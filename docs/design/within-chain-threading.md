# Within-chain threading for large-n single-chain sweeps: engineering design

Status: CLOSED - NO-GO on every tested hardware (x86 and Apple Silicon), 2026-07-21.

Summary: Single-chain large-n sweeps (n >= 1e5) are DRAM/latency-bound over a
handful of O(n) passes per tree. This note parallelizes those passes (suffstat
gather, fit scatter, residual roll, totalFits rebuild, SSR reduction) across a
persistent std::barrier worker pool, with a fixed-block reduction scheme that
makes every draw bitwise-identical across n.threads in {1, 2, 8} by
construction (section 3). An a-priori Amdahl model (section 6) projected a
modest ceiling of ~1.5-1.7x at n=1e5 and ~1.8-2.1x at n=1e6 for 4-8 threads.
The mechanism was built and measured on the real engine on two different
memory systems and came in below even that modest projection both times:
0.91x at the target size on an x86 Ryzen/DDR4 box (section 8), and 1.10x best
case with an 8-thread SLOWDOWN on Apple Silicon (section 10). A same-session
representative-kernel microbench had suggested ~3x scaling on M1 and proposed
reopening the question (section 10); that reversal was refuted once the
actual prototype was built and run on the real engine. The binding constraint
is not memory bandwidth but the PARALLEL FRACTION: only ~47% of the sweep is
the O(n) passes above, the rest (structure-move scans, per-leaf RNG draws,
bookkeeping) is inherently serial, and Amdahl's law over that fraction caps
any within-chain scheme at ~1.1-1.2x regardless of cores or bandwidth
(sections 6 and 10). The correctness half of the design works without
qualification: draws are byte-identical across n.threads in {1, 2, 8} on both
architectures. The reusable LESSON (section 11): a kernel microbench that
omits the serial fraction overestimates threading speedup by roughly
1/parallel-fraction; only an in-situ, real-engine measurement is trustworthy
for a threading claim. Multi-chain parallelism remains the effective way to
use additional cores. Companion document: docs/plans/blocked-jacobi-trees.md,
the competing noise-split mechanism for the same regime - independently
KILLED, and strictly dominated by this one wherever either shows a nonzero
gain.

Grounded in the current engine (bartcore HEAD, block-fusion excised: the
legacy setNodeAverages / computeLeafStats suffstat path is again the only
one). The estimates below are against the CURRENT u16 hot layer only. The
per-column-width (u8) layer is parked -- it is not in the tree (a standalone
width retrofit measured a regression, and per-column widths instead fold
structurally into the owned container, docs/design/data-ownership.md), so
there is no "on top of u8" variant to model. All numbers below assume u16
codes.

## 1. Which passes parallelize, which stay serial

The steady-state sweep (chain.hpp run(), the tree loop at :1457-1526 plus the
sigma/latent tail :1533-1545). Backfitting is sequential across trees; parallelism
lives INSIDE each pass, not across trees. Per tree t, in order:

- RESIDUAL ROLL (:1464-1470). A fused O(n) streaming pass: read totalFits and the
  tree's old/prev fit slabs, write the running residual treeY. PARALLELIZE
  (flat obs-order, split [0,n) into blocks). Streaming.
- SUFFSTAT GATHER (setNodeAverages -> computeLeafStats, :1469). Per leaf an
  O(n_leaf) GATHER-reduce of treeY over the leaf's shuffled member slice; O(n)
  over the tree. PARALLELIZE (a reduction, needs the fixed-block scheme,
  section 3). x86's #1 hotspot (~32%); latency-bound on the shuffled index
  buffer (block-fusion.md section 1, section 10).
- MOVE (metropolisJumpForTree, :1481). O(n_leaf) partition of one leaf's members
  plus 2 child suffstats, wrapped in an RNG proposal and accept/reject. SERIAL:
  small (one leaf), RNG-coupled, and the tree mutation is inherently sequential.
- LEAF-MEAN DRAW + SCATTER (sampleParametersAndSetFits, :1513). Draw one mu per
  leaf from a single RNG stream (SERIAL, O(#leaves)), then SCATTER-write each mu
  into the tree's fit slab (O(n), x86's #3 hotspot ~15%). The draw is serial;
  the SCATTER PARALLELIZES (writes only, no reduction ordering to preserve --
  partition leaves/blocks across threads).

Once per sweep, after the tree loop:

- TOTALFITS REBUILD (:1528-1530). O(n) streaming. PARALLELIZE.
- SIGMA DRAW (drawSigma, :1544-1545). An O(n) SSR reduction (fixed-block scheme) plus
  one gamma draw. PARALLELIZE the reduction; the draw is serial. Once/sweep.
- LATENT REFRESH (refreshLatents, :1535). O(n), embarrassingly parallel but
  RNG-CONSUMING per observation. DEFER (section 5). Once/sweep.

Per-pass work and traffic (u16 hot layer, m trees, doubles = 8 B):

  pass                per-tree traffic   x per sweep   regime by n
  residual roll       ~24 n B            m             streaming
  suffstat gather     ~16-24 n B         m             latency (shuffled gather)
  fit scatter         ~16 n B            m             latency (shuffled write)
  totalFits rebuild   ~24 n B            1             streaming
  SSR (sigma)         ~16 n B            1             streaming
  latent refresh      ~16-24 n B         1 (deferred)  streaming

The tree loop's ~3m O(n) passes dominate; the once/sweep passes are ~1/m of that
frequency. Per-tree working set (treeY + one index buffer + one fit slab) is
~24 n B, which sets the cache regime and thus what threading can buy:

  n       per-tree WS   resident in     dominant limit    matches
  1e4     ~0.24 MB      L2              cache latency     ~2.1 ms/sweep
  1e5     ~2.4 MB       L3/LLC          LLC latency       ~29 ms/sweep, ~8 GB/s
  1e6     ~24 MB        spills to DRAM  DRAM lat + bw      ~498 ms/sweep

The n=1e5 and n=1e6 anchors are block-fusion.md's measured b=1 sweep times on
the x86 bench box (section 10); n=1e4 is the prior per-sweep estimate. The full
n*m treeFits slab (60 MB at n=1e5, 600 MB at n=1e6) is only touched one tree
slab at a time in the loop, so the per-tree WS, not the slab, governs.

## 2. Substrate decision: a new std::barrier pool

Three candidates: revive the dormant hierarchical thread manager
(src/misc/hierarchicalThreadManager.c, zero live callers); extend the live
misc_mt pool (src/misc/thread.c, the substrate testFitPool_ already uses); or a
new C++20 std::thread + std::barrier pool. The workload is ~3 fan-out/join
regions per tree x m trees = ~3m per sweep (a few hundred at defaults), so the
decision hinges on per-sync overhead and on pool-persistence ergonomics.

Both misc_mt and the htm are CONDITION-VARIABLE (blocking) fork-join: each
dispatch locks a mutex, signals workers, and waits on a condvar for completion
(thread.c misc_mt_runTasks; the htm waits the same way). That is the wrong tool
for hundreds of fine-grained syncs. A standalone microbench measured the cost of
one fan-out+join round-trip (trivial per-worker work, isolating the sync), three
substrates, arm64 dev box (Apple Silicon, 10 hardware threads). AMBIENT LOAD was
present, so these are RELATIVE comparisons only, not absolute ceilings:

  threads   std::barrier   condvar(misc_mt model)   atomic spin
  2         ~0.55-0.6 us    ~4.7-5.8 us              ~0.22-0.28 us
  4         ~1.7 us         ~11-13 us                ~0.53-0.62 us
  8         ~4.8-7.3 us     ~27.5-27.9 us              ~1.45-2.0 us
  (us/round-trip; one round-trip = one parallel region)

The condvar fork-join is 5-10x more expensive per sync than std::barrier, which
is itself ~3-8x a pure spin. Projected onto ~225 regions/sweep (m=75) at n=1e5
(~29 ms sweep): condvar burns ~2.5 ms (4 threads) to ~6.3 ms (8 threads) =
9-22% of the sweep in pure sync; std::barrier burns ~0.4-1.6 ms = 1.3-5%; spin
~0.3-0.45 ms = ~1-1.5%. On the x86 bench box (real cores, no oversubscription at
8 threads) the std::barrier and spin numbers should fall further; the condvar
floor is a syscall/wakeup cost that stays. Condvar overhead alone would eat most
of a modest win, so reviving the htm or extending misc_mt is disqualified for the
hot loop.

RECOMMENDATION: a new persistent worker pool built on std::thread +
std::barrier. A pure spin barrier is cheaper still but burns T-1 cores through
the serial move/draw phase between the suffstat and scatter regions, and is
fragile under oversubscription (cross-chain workers times within-chain workers);
std::barrier (libc++'s spin-then-block atomic wait) keeps the low-contention
cost while parking cleanly when a phase runs long. Keep the misc_mt pool exactly
where it is -- testFitPool_ (chain.hpp:4063) is a COLD, coarse, once-per-run
test-fit fan-out, for which condvar dispatch is fine; this note does not touch
it. A spin-vs-block toggle on the new pool is a prototype tuning knob, not a
design fork.

## 3. Fixed-block reduction: bitwise thread-count invariance by construction

Floating-point addition is not associative, so a reduction whose grouping
depends on thread assignment gives different bytes at different n.threads. The
scheme that removes the dependence: partition each reduced range into
CONSTANT-SIZE blocks of B elements; reduce each block into a partial with a
scalar accumulator; combine the partials in BLOCK-INDEX ORDER. The block
boundaries and the combine order are fixed functions of the data layout, not of
how threads are handed blocks, so the sum is identical whether 1, 2, or 8
threads compute it. This applies to the reductions -- the suffstat gather (per
leaf, over its member slice) and the SSR. The roll, totalFits rebuild, and fit
scatter are not reductions (elementwise writes / independent stores) and only
need their index range partitioned; no ordering contract.

THE SERIAL PATH USES THE SAME SCHEME. n.threads = 1 walks the identical blocks in
the identical order, so the one-thread result equals the eight-thread result
bitwise. This is the extension of the existing cross-dispatch bitwise-invariance
component test to the threaded paths, and it is the hard gate (not a perf gate).

SCALAR BLOCK INTERIORS. The draw-path reductions are scalar today, which is what
makes results bitwise-identical across SIMD/ISA dispatch WITHIN a host. Each
block therefore reduces with a scalar accumulator; SIMD stays out of the draw
path. The reproducibility contract is unchanged and explicit: WITHIN-host bitwise
across thread count AND across SIMD dispatch; CROSS-host remains statistical
only (a data-generation and libm/compiler difference, not a reduction-order one),
exactly as today.

BLOCK SIZE. Pick B = 1024 doubles (8 KB, comfortably L1-resident). It is a
constant independent of n and of thread count -- the invariance contract forbids
tying B to either. At n=1e5 a flat range yields ~98 partials (a per-leaf slice
yields fewer); the extra combine adds << 0.5% arithmetic, blocks stay large
enough that per-block bookkeeping is amortized, and ~98 blocks over <=8 threads
load-balances well. B is a knob to tune when prototyping, but the CONTRACT
(constant B, block-index-order combine, scalar interior) is the load-bearing
part, not the value.

ONE-TIME RNG SHIFT. Versus today's straight serial accumulation, the block
grouping regroups the floating-point sum, so draws shift -- the "shifting" class.
This is a one-time change requiring snapshot regeneration (replay whole test
files) and an equivalence re-record, taken once when the threaded path lands.
If block fusion or any other shifting-class change lands in the same window they
should share ONE re-record rather than each paying its own.

## 4. Worker-pool lifecycle

The flagship consumer is an embedded single-chain Gibbs sampler calling
run(0, 1) once per outer sweep. Per-run thread startup (spawn + join ~tens of us
each) would swamp a per-sweep win, so the pool MUST persist across run() calls.

PRECEDENT. testFitPool_ (chain.hpp:4063, :895, :4058-4075) is exactly this
lifecycle: a pool held as a Chain member, lazily created on first use, resized
only when the budget changes, reused across calls, destroyed in the Chain
destructor. The new pool follows it structurally -- a persistent Chain member,
lazily created when within-chain threading engages, torn down in the destructor
-- but is a std::barrier pool (section 2), not a misc_mt one. Workers spin/park on
the barrier between sweeps; run(0, 1) reuses them with no spawn.

WORKERS NEVER TOUCH R. Within-chain workers only do arithmetic on the residual,
fit, and index buffers; they never call into R. This is the existing worker rule
(the ProgressSink pattern, sampler.hpp), satisfied trivially here.

INTERACTION WITH THE CROSS-CHAIN LAYER AND THE BUDGET SPLIT. Cross-chain
parallelism (sampler.hpp:374-375, :1264-1265) fans numChains across
numWorkers = min(numThreads, numChains) raw std::thread workers. Two cases:

- SINGLE CHAIN (numChains = 1): numWorkers = 1, the cross-chain layer is inert
  and chains run on the main thread, so the FULL numThreads budget is free for
  within-chain parallelism. This is the flagship case and the clean one.
- MULTIPLE CHAINS: each chain already occupies a core. The within-chain budget
  per chain is numThreads / numChainWorkers -- the same arithmetic testFitPool_
  uses (chain.hpp:4063, budget = numThreads / chains). Within-chain threading
  engages only when that per-chain budget is > 1 (i.e. numThreads > numChains);
  otherwise the chain runs its sweep serially. This keeps total live threads
  ~ numThreads and avoids oversubscription (cross-chain workers nesting
  within-chain workers). The single n.threads control governs both layers, so no
  new user-facing knob is added.

## 5. Latent refresh: deferred

refreshLatents (chain.hpp:1535) rewrites the working response (and, for logistic,
the Polya-Gamma weights) once per sweep. It is embarrassingly parallel over
observations but RNG-CONSUMING: each observation draws from the chain's single
RNG stream in obs order. Parallelizing it without breaking thread-count
invariance requires per-block RNG SUBSTREAMS -- each fixed block seeded
deterministically from the chain RNG plus its block index, so observation i's
draw depends only on its block and offset, never on which thread ran the block.

DEFER, for two reasons. First, cost/benefit: it is ONE O(n) pass per sweep
against ~3m O(n) passes in the tree loop, so it is well under 1% of the
parallelizable work at m=75 -- parallelizing it moves the ceiling by a rounding
error. Second, substreams are a separate RNG-architecture change: the engine's
ext_rng is a single sequential stream, not a counter-based / splittable
generator, so per-block substreams would themselves shift the stream and would
need their own reproducibility design. The tree-loop reductions (section 3) carry
none of that -- they regroup an ADDITION, they do not resplit the RNG. Keep
latent refresh serial; revisit only if profiling a bound problem shows it
material.

## 6. Expected-speedup ceiling

Combine two models. AMDAHL over the parallelizable O(n) fraction p, with the
parallel part scaling by s_par (capped by memory, not cores), minus the barrier
fraction:

  speedup(n, T) ~ 1 / [ (1 - p) + p / s_par(n, T) ]  -  barrier_frac(n, T)

The DRAM/latency accounting is taken from block-fusion.md section 10 (the
corrected model; the earlier ~6x block-fusion projection over-counted the
amortizable traffic by ~3x), NOT re-derived:

- Single-core already sustains ~8 GB/s effective at n=1e5, near the per-core
  streaming ceiling. Aggregate socket bandwidth on the x86 bench box gives the
  streaming passes headroom of roughly (socket bw / 8 GB/s) with enough threads.
- The two dominant passes -- suffstat gather (~32%) and fit scatter (~15%),
  together ~47% -- are LATENCY-bound on the shuffled index buffer and
  LLC-resident at n=1e5, not streaming. This is where threading beats block
  fusion's pessimism: bandwidth amortization could not touch latency, but
  multiple cores multiply outstanding misses (memory-level parallelism), so the
  latency-bound gather scales with threads until the LLC/DRAM ceiling, whereas
  block fusion left it untouched.

Parameterize: p ~ 0.6 (suffstat + scatter + roll + totalFits + SSR are the O(n)
passes; the serial remainder is the move logic, the per-leaf RNG draws, and tree
bookkeeping). barrier_frac from section 2 (std::barrier, x86-adjusted).

  n=1e5 (LLC latency-bound; s_par ~2.5 at 4T, ~3 at 8T; barrier ~1.3-2.3%):
    4 threads:  ~1.55x    8 threads:  ~1.6x
  n=1e6 (DRAM lat+bw; s_par ~3.5 at 4T, ~5 at 8T; barrier < 0.3%):
    4 threads:  ~1.75x    8 threads:  ~1.95x

So the ceiling is ~1.5-1.7x at n=1e5 and ~1.8-2.1x at n=1e6. This is honestly
MODEST, and it is a CEILING -- LLC contention at 8 threads, a higher serial
fraction, or false sharing on the partials would pull it down. The upside over
block fusion is not magnitude but risk: threading runs the same kernels with no
O(mn) bookkeeping, so its worst case is "neutral," and it addresses the
latency-bound passes block fusion could not.

SIZE CUTOFF. At n=1e4 the per-tree working set is L2-resident (little latency to
hide, s_par ~1.2-1.5) while barrier overhead is large relative to the ~2.1 ms
sweep (~225 regions x std::barrier: ~18% at 4T, ~50%+ at 8T) -- net neutral to a
LOSS. Threading must therefore be gated OFF below a size cutoff, with the
within-chain thread count forced to 1 there (the same pattern as the SIMD size
toggles). The cutoff sits somewhere between n=3e4 and n=1e5; its precise value
is a measurement output. A conservative default of engage-at-n >= 1e5 matches
the target regime and is safe to ship pending the measured crossover.

(As it turned out, section 8 measured the prototype below this cutoff too, and
found it byte-identical to a threads=1 run and cost-neutral -- see section 8's
n=1e4 row.)

## 7. Prototype gates and the head-to-head

PROTOTYPE. Thread the suffstat gather and the fit scatter only (the two
hotspots), behind the fixed-block scheme, and measure at n in {1e4, 1e5, 1e6} x
threads {2, 4, 8} on the x86 bench box, quiet, against the u16 hot layer. The
roll / totalFits / SSR passes are cheaper streaming variants of the same
machinery and can wait for the landing pass. The FIRST prototype action is to
re-run the section-2 barrier microbench on the bench box and on the real kernels
(false sharing on the partial array and NUMA effects only show on the true
workload), since the arm64 numbers here are relative-only.

GO / NO-GO.
- GO if, at n=1e5, net speedup is >= ~1.4x at 4 threads and the bitwise
  invariance across n.threads in {1, 2, 8} holds (the hard gate). The invariance
  gate is non-negotiable regardless of speed.
- NO-GO (record it, as block-fusion.md section 10 did for its bet) if the
  n=1e5 4-thread net is < ~1.3x, or the latency-bound gather fails to scale
  (s_par < 2 at 4 threads): that means the memory wall dominates and the
  complexity is not repaid. "Still not worth it" is a valid, recorded outcome.
- No regression at n=1e4 with threading compiled in but below its cutoff (the
  cutoff logic must keep it off cleanly).

HEAD-TO-HEAD vs BLOCKED-JACOBI. The competing mechanism (noise-splitting
augmentation, parallel-bart-frontier.md section 3.5) attacks the same
single-chain regime differently: augment per-tree pseudo-responses summing to
the batch residual, then update b trees INDEPENDENTLY -- m/b barriers per sweep
instead of ~3m, each worker a cache-coherent single-tree update rather than a
reduction slice. The two mechanisms are exclusive enough to bench against each
other on the SAME hardware (x86 bench box, quiet), same data (friedman, n in
{1e4, 1e5}, m in {75, 200}), at matched thread counts {2, 4, 8}. The common
metric is ESS-PER-SECOND, not wall-clock per sweep: this note preserves the exact
posterior (thread-count invariant, identical draws), so its ESS/sweep is
unchanged and ESS/sec = its speedup x baseline; blocked-jacobi changes the kernel
(structure moves see precision b/sigma^2, so mixing degrades with b) so its
ESS/sweep falls while its sweeps parallelize more cleanly. ESS/sec is the only
honest comparison. Whichever mechanism wins gets the implementation; the loser is
recorded as closed. Secondary axis: correctness risk -- this note is a mechanical
parallelization with a one-time RNG shift, whereas blocked-jacobi is a new exact
kernel needing exact-posterior gates at b in {2, 8}; a tie on ESS/sec breaks
toward the lower-risk mechanism.

(How this played out: sections 8-10 below record this mechanism's own measured
outcome on x86 and Apple Silicon. It no-goed on its own gates before the
head-to-head could run under this protocol. blocked-jacobi-trees.md records
that mechanism's independent Phase 0/Phase 1 evaluation and an eventual
head-to-head against this one's results -- see section 10 and section 11.)

## 8. Measured outcome: NO-GO (2026-07-13)

The section-7 prototype was built (persistent std::barrier pool in
wcpool.hpp, fixed-block gather, partitioned scatter; +306 lines) and
measured on the x86 bench box (AMD Ryzen 3700X, 8 cores across 2 CCXs with
a split L3; friedman, m = 75, single chain). The correctness half of the
design WORKED: draws byte-identical across n.threads in {1, 2, 8} at
n = 1e5 on both arm64 and x86, the below-cutoff path byte-identical to a
threads=1 run, component tests green, and the serial fixed-block regroup
cost nothing measurable at any n. The performance half did not:

  msec/iteration        tip 1T   proto 1T   2T     4T     8T
  n = 1e4 (gated off)   2.09     2.07       2.07   2.00   1.94
  n = 1e5               26.4     25.2       23.7   27.6   39.7
  n = 1e6               359      359        319    347    403

- Net at n = 1e5 x 4T: 0.91x against the >= 1.4x go / < 1.3x no-go gate.
  Best anywhere: 1.12x at n = 1e6 x 2T.
- Gather s_par at 4T: 1.67 against the < 2 no-go gate (2T: 1.36). At 8T the
  gather COLLAPSES to 0.86x once threads span both CCXs and leaf sharing
  turns into cross-CCX L3 traffic.
- The fit scatter (bandwidth-bound writes) was NET-NEGATIVE to parallelize.

Extrapolation does not rescue the landing pass: threading the remaining
streaming passes lifts the parallel fraction from ~0.47 to ~0.6, but at the
measured s_par -- and with the scatter already net-negative -- the ceiling
stays ~1.2-1.3x at 4 threads, under the bar the complexity was priced at.
The section-6 model was still optimistic: it assumed s_par ~2.5 for the
latency-bound gather where the machine delivered 1.67. Synchronization is
NOT the culprit: the re-run microbench put std::barrier at ~1.1/1.9/5.0 us
per round-trip at 2/4/8 threads on real cores (condvar 3-8x worse, spin
~0.2-1.1 us), ~1% of a n = 1e5 sweep. The wall is memory, as it was for
block fusion (block-fusion.md section 10).

The head-to-head never ran under the section-7 protocol: this mechanism
no-goed on its own gates first. At the time, blocked-jacobi remained
unevaluated, open on its own merits, with the section-7 ESS/sec protocol
still the right yardstick if it were ever tried -- it subsequently was; see
blocked-jacobi-trees.md and section 11 below.

The prototype is archived, buildable, on archive/within-chain-threading.
Revival preconditions: hardware whose memory system actually scales with
cores for this footprint (large unified LLC, materially higher bandwidth
per core), or routine workloads at n >> 1e6, or a sweep profile reshaped by
other layers; re-run the barrier microbench and the full section-7 protocol
on the actual target before believing anything.

## 9. fp32 re-check (2026-07-20): still NO-GO - fp32 does NOT revive it

The reduced-precision storage arc (uint32 index + opt-in fp32 residual, both
landed) is exactly the "sweep profile reshaped by other layers" precondition
section 8 flagged. Re-checked with a multi-threaded ISOLATED gather microbench
(fp32 vs fp64 storage, s_par at T=1,2,4,8, n in {1e5,1e6,1e7}, M1 + x86 box).
Two decisive findings:
- The isolated gather ALREADY parallelizes fine for BOTH precisions: s_par(4T)
  ~3.8x at n=1e5, ~2.3-3.8x at n=1e6, ~2x at n=1e7. So the gather was NEVER the
  parallelism bottleneck - the section-8 in-situ 1.67x came from BARRIER overhead
  + the net-negative fit-scatter + the serial fraction, none of which fp32 or the
  index narrowing touch.
- fp32 LOWERS s_par, it does not raise it (e.g. M1 n=1e7: fp64 s_par 2.54x vs
  fp32 1.92x). A faster, more-cache-resident single core has LESS memory latency
  to hide across threads, so reduced-precision storage makes within-chain
  parallelism LESS attractive, not more - the opposite of the revival hypothesis.
  (fp32's absolute throughput is still best at every thread count; only the
  parallel headroom shrinks.)

Verdict: fp32 does not flip the NO-GO. The revival preconditions that remain are
the ORIGINAL ones (a memory system that scales with cores for this footprint, or
n >> 1e6 routine single-chain), NOT the storage reductions. blocked-jacobi
(the exact noise-split kernel) remained, at this point, separately unevaluated
on its ESS/sec merits (section 7) -- see section 11 for how that concluded.

## 10. Apple Silicon: a microbench suggested revival; the real engine refuted it (2026-07-21)

Section 8's NO-GO was measured ONLY on the x86 bench box (Ryzen 3700X,
dual-channel DDR4, split L3 across CCXs) -- a bandwidth-bound,
cross-CCX-penalized machine, and section 8 named the revival precondition
precisely: "hardware whose memory system actually scales with cores -- large
unified LLC, materially higher bandwidth per core." That describes Apple
Silicon, which had never been tested.

A representative-kernel microbench (spin barrier, m=200 friedman, gather+scatter
field kernels; bj-wallclock-probe.cpp, job b073bb28) run ON an M1 showed
straight within-chain threading (data-parallel gather with a per-worker bucket
reduction, plus data-parallel scatter; EXACT/bitwise) SCALING:

    n=1e5:  T=2 1.48x  T=4 2.25x  T=8 3.04x
    n=1e6:  T=2 1.19x  T=4 2.02x  T=8 3.08x

i.e. ~3x ESS/sec at 8 threads on M1 (exact, so wall-clock IS ESS/sec, no tax),
versus the 0.91x LOSS measured on the Ryzen. On the same M1 this microbench also
beat blocked-jacobi head-to-head (blocked-jacobi-trees.md: 1.56-1.62x ESS/sec at
8T), because the straight mechanism is exact and carries none of noise-splitting's
extra RNG/scratch traffic. This looked like a genuine reopening: the revival
precondition section 8 had named appeared to be met, so the obvious next step
was to build the archived real-engine prototype (this branch: chain.hpp +
wcpool.hpp, +306 lines) on M1, run bench-sampler at n.threads in {1, 4, 8} on a
quiet Mac (n in {1e5, 1e6}, single chain), confirm the ~2-3x, and re-check the
correctness gate (byte-identical across thread count). A spin barrier
(wcpool-spin.hpp), which beat std::barrier even on M1 in the probe, was flagged
as part of any revival substrate.

That test was run, and it refuted the microbench. The archived prototype
(commit 54a60aa, the real engine, byte-identical draws reconfirmed) was built
to a temporary lib and benched on a QUIET M1, single chain, timing ONLY the
sampling loop (sampler$run, ingestion excluded):

    n=1e5:  1T 64.2 ms/iter   4T 1.10x   8T 0.53x
    n=1e6:  1T 691  ms/iter   4T 0.98x   8T 0.67x

BEST anywhere was 1.10x (4T, n=1e5); 8T is SLOWER than 1T at both sizes (M1's 4
performance + 4 efficiency cores mean the barrier ends up waiting on an
efficiency core once threads spill past the performance cluster). This is
essentially the same ceiling as the x86 result (best 1.12x, section 8). The
section-10 microbench's ~3x reading OVERESTIMATED BY ~3x because it modeled the
sweep as two clean, fully parallel gather+scatter passes (parallel fraction
~1.0); the REAL sweep is only ~47% parallel (section 6) -- the structure-move
scans, the residual/totalFits/SSR bookkeeping, and the sigma and leaf draws are
all SERIAL and were not threaded in the microbench's model. Amdahl on p=0.47
caps within-chain threading at 1/(1-0.47) = 1.9x even with INFINITE s_par, and
at the actually-measured s_par (~1.67, memory-bound, section 8) the ceiling
comes out to ~1.2x -- exactly what both machines deliver in the real-engine
test.

VERDICT: high-bandwidth memory does NOT revive within-chain threading. The
binding wall is the SERIAL FRACTION combined with s_par, not bandwidth alone.
Section 8's NO-GO stands on Apple Silicon as well as x86 -- there is no
revival. LESSON: a representative-kernel microbench that omits the serial
fraction overestimates threading speedup by the reciprocal of the parallel
fraction; trust the IN-SITU measurement (section 8 was right all along). The
one theoretical path left is threading the remaining serial passes too, which
section 8 already modeled (lifting p toward ~0.6): even then the ceiling is
only ~1.2-1.3x, under the bar the complexity was priced at. Multi-chain
parallelism remains the answer for using cores.

## 11. Lesson and final state

A representative-kernel microbench that omits a workload's serial fraction
overestimates threading speedup by roughly the reciprocal of the parallel
fraction -- here, omitting the ~53% serial remainder inflated the apparent M1
win by about 3x. This is not specific to this mechanism: the same flaw
independently inflated the blocked-jacobi-trees.md wall-clock probe numbers
(see that document's final verdict). The methodological lesson to carry
forward: gate any threading claim on an in-situ, real-engine measurement (the
actual sampler loop, actual sweep, actual bookkeeping), never on a microbench
of the hot kernels alone, however representative the kernels look in
isolation.

Final state:
- VERDICT: NO-GO on every tested hardware (x86 Ryzen/DDR4, Apple Silicon M1).
  Both revival preconditions section 8 named (a memory system that scales with
  cores, or n >> 1e6 routine single-chain) were retested where feasible
  (Apple Silicon) and did not change the outcome.
- BANKED: the correctness mechanism (fixed-block reduction, section 3) works
  exactly as designed -- draws are byte-identical across n.threads in {1, 2,
  8} on both architectures, and the below-cutoff serial path is byte-identical
  to a threads=1 run. This substrate (wcpool.hpp / wcpool-spin.hpp) is
  archived and buildable on archive/within-chain-threading, and was reused
  as-is for the blocked-jacobi-trees.md prototype.
- The head-to-head against blocked-jacobi under the section-7 protocol never
  ran on this mechanism's own numbers, because this mechanism no-goed on its
  own gates first (section 8). A later, informal comparison using the
  section-10 microbench's M1 numbers (this mechanism at ~3x vs blocked-jacobi
  at 1.56-1.62x) is superseded by the same refutation and carries no standing;
  the comparison that matters is both mechanisms' real-engine results, which
  are both NO-GO.
- blocked-jacobi-trees.md is independently KILLED (memory-bandwidth-bound on
  typical x86, and in any case dominated by whatever this mechanism delivers,
  which is itself NO-GO); see that document for its own banked result (the
  noise-split augmentation is proven exact with ~no per-sweep ESS tax --
  reusable knowledge, not a shippable win).
- Multi-chain parallelism remains the effective way to use additional cores
  for this workload.
