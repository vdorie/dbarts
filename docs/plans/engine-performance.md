# engine-performance

agent: opus (S1 the build flag and its report; S2 the vector suffstat
  kernel; S3 the wait loop; S4 within-chain threading; S5's engine
  plumbing and settings); sonnet (S1's CI axis, S5's scripts, manual
  table and NEWS). S1 before S2; S3 is file-disjoint; S4 after S2, since
  it re-blocks the same accumulation; S5 after front-door's
  family-object slice (below).
rng: S1 neutral on both builds - nothing moves yet, and both builds
  reproducing the baselines bitwise is what proves the flag inert. S2
  shifting in the SHIPPED build, neutral in the REFERENCE build: the
  vector kernel re-associates the leaf sufficient statistic, the
  posterior does not move. dec-B90 licenses an exception to the shifting
  class's usual gates - NO snapshot replay, NO baseline re-record -
  because the gates split by build: the three equivalence baselines and
  the four snapshot files stay recorded on the reference build, whose
  reproducing them bitwise IS the neutrality proof, while the shipped
  build compares against those same unchanged baselines in statistical
  (z) mode with the exact-posterior gates as its oracle under the P17
  rule (a draw-changing re-record names the oracle establishing the new
  values; see MANIFEST). The baselines have always been recorded on
  arm64 macOS and the reference build must stay bitwise there; the
  shipped build's cross-ISA variation is accepted by dec-B90, narrowable
  by fork 2. S3 neutral. S4 neutral at the default (opt-in off keeps
  today's serial association); the opt-in carries its own fixed-block
  association, bit-identical across thread counts, with its own
  invariance test (fork 5). S5 neutral at every default; a moved default
  is posterior-changing.
window: engine-only, so it runs beside the R-only front-door arc
  (docs/plans/front-door.md), whose family-object slice moves the
  person-period row cap onto a `hazard()` constructor - S5 reads the cap
  from that new home, so it lands after. S1 edits configure.ac, as does
  docs/plans/interfaces-and-dependencies.md's configure-stub slice; both
  regenerate configure on the "no AC_ARG today" premise, so S1 lands
  first and the stub slice rebases. S4's default waits on
  docs/plans/memory-footprint-audit.md.
budget: S1 ~40 configure and config headers, ~60 bridge, ~40 R, ~40
  tinytest, ~90 CI; S2 ~350 misc kernels and selection, ~40 stats.h, ~80
  tests/cpp, ~40 docs; S3 ~60 sampler.hpp, ~40 tests/cpp; S4 ~120 note,
  ~330 engine, ~60 R, ~90 tests, ~60 Rd; S5 ~200 scripts, ~150 R, ~120
  engine, ~150 tests, ~250 note and Rd.

Decisions: dec-B88 (the run loop waits on completion), dec-B90 (vector
draw path, scalar reference build), dec-B89 (within-chain threading as
an opt-in), dec-B91, B93 and B110 (the constants audit, the calibrated
predict cutoff, the GP fallback warning), in
[docs/decisions.md](../decisions.md); they supersede dec-A37, A38, A42,
A43, A45, A46 and the never-built dec-B73.

## Goal

The shipped build vectorizes the leaf sufficient-statistic gather and
reduction; a configure flag builds the scalar fixed-order kernel
instead, and R reports which build it is running. The bitwise gates run
on the reference build, the suite and statistical gates on the shipped
one. A multi-chain run returns when the last chain finishes, not at the
next 100 ms tick. Within-chain threading is back as an explicit opt-in,
re-measured at the tip first. Every fixed engine constant is documented,
measured where possible, and the ones that bind become settings.

## Context

- The run loop: [`Sampler::run`](../../src/bartcore/sampler.hpp) spawns
  `min(numThreads, numChains)` raw threads that stride chains and
  decrement `numChainsRunning` while the caller spins on a 100 ms sleep,
  polling the interrupt and flushing `QueuedProgressSink` per wake; misc
  has [`misc_mt_runTasksWithInfo`](../../src/include/misc/thread.h).
- The draw-path reduction is NOT one kernel. The default gaussian
  constant-leaf unweighted path takes the in-header fused roll plus
  node-average pass, [`rollAndSetNodeAveragesFused`](../../src/bartcore/chain.hpp),
  whose association is already fixed at four positional banks
  ([`fusedSuffstatBanks`](../../src/bartcore/chain.hpp),
  [`combineFusedSuffstatBanks`](../../src/bartcore/chain.hpp)) that the
  engine documents as part of the draw law - independent of ISA, lane
  width and worker count, and explicitly refused a knob. What that pass
  declines (weighted families, fp32 residuals, a stale leaf map, linear
  and GP leaves) and every per-move child statistic
  ([`computeLeafStats`](../../src/bartcore/tree.hpp)) goes instead to the
  eight scalar single-accumulator entry points of
  [`misc_computeSufficientStatisticsFast`](../../src/include/misc/stats.h)'s
  family - plain, indexed, weighted and indexed-weighted, each with a
  Float twin - which sit outside the
  [`misc_simd_init`](../../src/misc/simd.c) table and are S2's target.
  The 32 percent share and net 3 to 4 percent estimate in
  [RANKED IMPLEMENTATION PLAN (x86-measured)](x86-simd-plan.md#ranked-implementation-plan-x86-measured)
  predate the fused pass, so S2 re-profiles before it writes anything.
- Vectorizing a reduction ends host independence: every dispatched double
  kernel today is elementwise or a permutation
  ([THE KEY ARCHITECTURAL INVARIANT (drives every verdict below)](simd-survey.md#the-key-architectural-invariant-drives-every-verdict-below)),
  and that is enforced at run time:
  ["C_dbarts_setSIMDInstructionSet"](../../inst/tinytest/test-simd.R) walks
  every dispatch level and pins yhat.train, sigma and varcount bitwise
  against the level-0 scalar fit, which check-standard.yaml's
  windows-arm64-neon job - the only NEON gate - asserts too. A suffstat
  kernel inside the level switch breaks both on every host, so S2 selects
  it by build mode alone.
- Build configuration: configure.ac carries no `AC_ARG_ENABLE` and four
  `AC_CONFIG_HEADERS`, so one `AC_DEFINE` reaches src/config.hpp (the
  bridge) and src/misc/config.h (the kernels) together. Windows has no
  configure, and [check-win-drift.R](../../tools/check-win-drift.R) fails
  on a `.in` macro absent from its `.win` counterpart unless its table
  records it. R learns nothing about the build today.
- Within-chain threading: the prototype is one commit on
  origin/archive/within-chain-threading, +306 lines over chain.hpp and a
  new wcpool.hpp - persistent barrier pool, fixed-block gather,
  partitioned fit scatter, draws byte-identical across thread counts by
  construction. The design doc is CLOSED NO-GO on speed
  ([8. Measured outcome: NO-GO (2026-07-13)](../design/within-chain-threading.md#8-measured-outcome-no-go-2026-07-13),
  [11. Lesson and final state](../design/within-chain-threading.md#11-lesson-and-final-state));
  this arc does not reopen that, it ships the opt-in because losing it
  regresses 0.9-34, and re-measures under section 11's rule.
- On the R side `n.threads` is a [`dbartsControl`](../../R/dbarts.R)
  slot defaulting to `dbarts::guessNumCores()` uncapped, capped to
  `n.chains` only inside [`bart`](../../R/bart.R) (front-door.md S1
  renamed `bart2` to `bart`; the capping logic moved with the body), so
  above four cores
  the control's own default exceeds its default chain count. The
  calibrated-cutoff precedent is threaded predict
  ([8. Measurement, honestly stated](../design/threaded-predict.md#8-measurement-honestly-stated)):
  [`predictParallelCutoff`](../../src/bartcore/sampler.hpp) has a test-only
  `cutoffOverride`, the seam a setting reuses.

## Decision

All six forks were put to VD on 2026-09-08 and are recorded with the
choice; none remains open.

1. The flag's name and the Windows route (VD 2026-09-08, "Use your
   recommendation"): `--enable-reference-build` defining
   `DBARTS_REFERENCE_BUILD`; on Windows an environment variable read by
   src/Makevars.win sets the same macro, the `.win` headers stay
   untouched and the macro joins the check-win-drift table; the
   windows-arm64-neon job stays on the shipped build with its bitwise
   assertion. CRAN ships the vector build everywhere.
2. The vector kernel's accumulator layout: decided by measurement (VD
   2026-09-08, "Shouldn't we do research and evaluation to decide?").
   S2 builds both layouts of each vectorized kernel, a fixed four-bank
   split matching [`fusedSuffstatBanks`](../../src/bartcore/chain.hpp)
   and the natural per-ISA width, and times them in the real sampler
   loop on the arm64 Mac (and on the x86 bench box if VD grants it) at
   n in {1e4, 1e5, 1e6} and 75 and 200 trees. Rule: if the natural
   width buys under 2 percent of a whole fit at every cell, the fixed
   split ships and the shipped build stays bitwise across machines,
   thread counts and instruction sets; if it buys more anywhere, the
   fork returns to VD with the table and the reproducibility cost
   stated (cross-host CI falls to its statistical tier). The design
   note records the table either way.
   Continued (VD 2026-09-10): "Why don't we see how it works on weights
   first before deciding. It doesn't seem that complicated (we had it
   before), and if it makes a difference on weights then we'll already
   have an instance where bit-wise equivalence won't apply." Step 5
   therefore narrows to the two WEIGHTED double entry points, which the
   re-profile puts at 31 to 54 percent of a weighted fit against 5 to
   12 percent for the default path's single indexed kernel; the
   unweighted pair and the four Float twins stay scalar, and the layout
   rule above is settled on the weighted cells.
   MEASURED (2026-09-10, arm64 M1 Max, macOS 26.6.2, Apple clang 21.0.0,
   R 4.6.1, one chain, one thread, otherwise idle; `sampler$run` only,
   7 repeats interleaved round by round across the three builds, median
   seconds, p = 20 throughout):

   | cell | reference | fixed split | natural width | fixed vs reference | natural vs fixed |
   |---|---|---|---|---|---|
   | weighted, n = 1e4, 200 trees, 300 + 300 | 2.937 | 2.914 | 2.978 | +0.78% | -2.18% |
   | weighted, n = 1e5, 200 trees, 100 + 100 | 11.494 | 11.472 | 11.583 | +0.19% | -0.97% |
   | weighted, n = 1e5, 75 trees, 100 + 100 | 4.294 | 4.290 | 4.302 | +0.09% | -0.28% |
   | default (no weights), n = 1e4, 200 trees, 300 + 300 | 2.790 | 2.786 | 2.783 | +0.14% | +0.11% |

   Spread was tight - the widest min-to-max across the 7 repeats of any
   cell is 2 percent of that cell's median, and every percentage above
   is a fraction of a whole fit, not of the kernel.

   OUTCOME: the natural width buys under 2 percent at every cell (it
   LOSES at every weighted cell), so by the rule above the FIXED SPLIT
   ships and the natural-width code is deleted. The shipped build's
   weighted draws therefore do not depend on ISA, lane width or
   dispatch level, only on the build mode.
3. Settled with fork 1.
4. Within-chain threading's opt-in spelling (VD 2026-09-08: "No,
   n.threads shouldn't mean n.chains"): `n.threads` keeps 0.9-34's
   meaning, a total thread budget ("used for various internal
   calculations, as well as the number of chains" in 0.9-34's manual);
   chains run in parallel across it and, when the budget exceeds
   `n.chains`, the surplus is divided among the chains and used inside
   each. Supplying more threads than chains is the explicit opt-in; no
   new argument and no warning; `dbartsControl`'s default budget stays
   the core count as in 0.9-34. Step 11 changes accordingly.
5. The summation law under within-chain threading (VD 2026-09-08, "I
   don't care about recreating the exact result, but I do care about
   preserving its speed"): two implementations. At or below the chain
   count the fused four-bank pass runs unchanged, keeping its measured
   speed; above it the fixed-block sum runs, identical at any surplus
   size. The manual states that shipped-build draws depend on the
   thread budget as they depend on the machine's vector width, and
   that the reference build is where exact reproduction lives; the
   chain.hpp comment on `fusedSuffstatBanks` is amended to say the
   budget selects the implementation. Gates for the threaded path:
   identical draws across surplus sizes and statistical equivalence
   against the serial path.
6. The GP fallback warning threshold (VD 2026-09-08, "Use your
   recommendation"): 25 percent of leaf evaluations that fell back to a
   constant leaf because the node exceeded
   [`maxLeafSize_`](../../src/bartcore/model.hpp); the count itself is
   exposed on the fit so any threshold can be checked by hand.

## Constraints

- Gates by class per [RNG classes and their gates](README.md#rng-classes-and-their-gates);
  an equivalence run is counted per scenario, not off the summary line
  ([Gate hygiene](README.md#gate-hygiene)). Bench compares are
  maintainer-run on the quiet machine, never concurrent.
- S2 selects the suffstat kernel by build mode only, never through
  [`misc_stat_setSIMDInstructionSet`](../../src/misc/simd.c), so
  test-simd.R's level-invariance pin and the windows-arm64-neon job keep
  passing on both builds. S2 must not touch the fused pass. S4 does:
  fork 5(a) overrides that pass's standing refusal of a knob under
  dec-B89 and amends its comment in the same slice.
- The reference build must reproduce equivalence-deb144d2.rds,
  bcf-equivalence-fbff1989.rds and multinomial-equivalence-fbff1989.rds
  bitwise on arm64 macOS at every slice; if it cannot, the slice stops
  and the re-record owes a P17 oracle
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)). The four RNG-locked
  snapshot files [regenerate-snapshots.R](../../tools/regenerate-snapshots.R)
  names skip off reference mode; their values do not move here.
- Out of scope: the relayout of per-node fits into index-buffer order
  (the real gather lever, a redesign); FMA (measured zero); the
  Mean/Variance family (no callers); widening the cut code past 16 bits,
  which [`misc_xint_t`](../../src/bartcore/tree.hpp) static-asserts and
  five per-ISA partition units carry, so it is its own item; and
  vectorizing [`misc_computeSumOfSquaredResiduals`](../../src/include/misc/stats.h),
  a different reduction feeding the sigma draw, which would move more
  than this arc's rng line claims.

## Steps

S1, the reference build and its report:

1. `AC_ARG_ENABLE` in configure.ac for the chosen flag, default off,
   with `AC_DEFINE(DBARTS_REFERENCE_BUILD, 1, ...)` so it lands in
   src/config.hpp and src/misc/config.h together; `#undef` in both `.in`
   templates, regenerated with `autoreconf -i`. Windows per fork 3:
   src/Makevars.win reads an environment variable of the same name and
   appends `-DDBARTS_REFERENCE_BUILD` to `PKG_CPPFLAGS` and the sub-make
   line, the macro going on
   [check-win-drift.R](../../tools/check-win-drift.R)'s expected-absent
   table so the `.win` headers stay untouched.
2. A bridge entry registered as `dbarts_buildInfo`, reached as
   `dbarts:::buildInfo()`, returning the mode (reference or shipped),
   the compiled ISA set and the dispatch level chosen at load; one R
   wrapper, one tinytest assertion on names and mode. The four snapshot
   files gain a mode-keyed `exit_file()` guard, and one test asserts it
   by count: each file yields zero assertions on the shipped build and
   its full recorded count on the reference build.
3. CI build axis: cpp-tests.yaml today installs `--preclean` on arm64
   and runs the C++ component tests only, so ADD a reference arm that
   also runs the four snapshot files and the three equivalence compares,
   drop `benchmarks/**` from its paths-ignore since those compares read
   benchmarks/, and raise timeout-minutes from 15 to cover a second
   `--preclean` install plus the compares. exact-gates.yaml runs its two
   cross-host compares on the reference build (a second install in that
   job), because their tier-1 verdict is a tight deviation bound the
   shipped build cannot meet and this workflow is the per-push half of
   that gate; the rest of exact-gates, equivalence.yaml,
   check-standard.yaml and sanitizers.yaml stay shipped. Prove the flag
   inert here: both builds bitwise.

S2, the vector suffstat kernel:

4. Re-profile all eight suffstat entry points at the tip with the fused
   pass live, on the quiet machine and the x86 box, recording the share
   carried on a default gaussian fit and on a weighted family. If the
   default-path share is under one percent the slice says so and stops.
5. Write per-ISA kernels under the fork-2 layout, no FMA, no hardware
   gather (measured slower than scalar), for the four double entry
   points and the four Float twins - the weighted and Float variants
   being exactly what the fused pass declines. Install them by build
   mode, where `DBARTS_REFERENCE_BUILD` is visible and NOT in the
   [`misc_simd_init`](../../src/misc/simd.c) level switch; under the flag
   the scalar bodies stay, so the reference build is today's code.
6. tests/cpp: a deterministic test that vector and scalar node sums
   agree within a mixed absolute/relative bound on a fixed fixture,
   modelled on [`testFusedSuffstatMatchesStock`](../../tests/cpp/test_sampler.cpp),
   plus prologue-residue coverage per
   [`testGatherTailShapes`](../../tests/cpp/test_sampler.cpp). The
   MANIFEST header records which build each gate uses.

Step 5 (weighted) outcome, measured at 80ff32b4: the vectorization is
worth 0.78, 0.19 and 0.09 percent of a whole weighted fit at the three
cells above - UNDER the two percent bar everywhere, and under one
percent everywhere. The default (unweighted) cell does not move, as
that path calls neither kernel. The mechanism is that
`misc_computeIndexedWeightedSufficientStatisticsFast`, which the
re-profile puts at 31 to 54 percent of a weighted fit, is entirely
gather-bound: a kernel-level A/B of the two bodies on this host times
them at 1.00x for every node of 128 observations or more (1.27x at 8,
1.07x at 32), while the contiguous twin - about 1 percent of a fit -
gains 1.5x at mid sizes and 1.01x at 1e5. So the cost of the change,
which is that every weighted family's draws move on the shipped build,
buys well under the "about 3 to 4 percent of runtime" dec-B90 was
weighed against. Shipping it at all is VD's call and is NOT settled
here; what is settled is the layout, if it ships. Two consequences to
weigh: bcf-equivalence.R and multinomial-equivalence.R gate their
point-in-time snapshot channels (mu, tau, glue; forestFits) bitwise
with no statistical fallback, so on the shipped build both harnesses
FAIL on all 12 and all 11 scenarios even though every draws-axis
channel reports max |z| = 0.00 - only their reference-build and
`--cross-host` runs stay meaningful, which is where every workflow
already runs them; and benchmarks/baselines/MANIFEST's build-mode
header and the equivalence-deb144d2 row still say the shipped build
reproduces the baselines bitwise, which stops being true the moment
this lands.

AVX2 addendum (VD 2026-09-10, "Try the AVX2"). A four-wide body of the
SAME split - one 256-bit register holding all four banks, lane b being
bank b, the prologue into bank 0, the combine ((b0 + b1) + b2) + b3,
products named before accumulation, no FMA, no hardware gather - was
written into its own `-mavx2` translation unit and installed by CPU
detection, the only draw-path reduction the package has ever put behind
a dispatch pointer. Byte-identity is what makes that safe, and it was
asserted directly: the tests/cpp arm compares the entry points at
dispatch level 0 against the host's own maximum, bitwise, over 27
shapes. On real AVX2 hardware it passed at level 8 against level 0, and
a mutation - the AVX2 combine regrouped pairwise - failed it, so the
assertion discriminates.

Timed on a 4-core x86-64 Linux box with AVX2 and FMA, gcc 13.3, R 4.6.1,
otherwise idle; same four cells, same protocol as the table above
(7 repeats interleaved, median seconds, one chain, one thread):

| cell | reference | split (SSE2) | AVX2 | AVX2 vs split |
|---|---|---|---|---|
| weighted, n = 1e4, 200 trees, 300 + 300 | 1.957 | 1.962 | 1.960 | +0.10% |
| weighted, n = 1e5, 200 trees, 100 + 100 | 9.817 | 9.755 | 9.911 | -1.59% |
| weighted, n = 1e5, 75 trees, 100 + 100 | 3.736 | 3.768 | 3.730 | +1.02% |
| default (no weights), n = 1e4, 200 trees, 300 + 300 | 1.548 | 1.535 | 1.545 | -0.65% |

Nothing there is a signal: the sign flips between cells, and the default
cell - which calls neither kernel and cannot move - drifts by the same
magnitude as the weighted ones, which fixes the noise floor at about
1 percent on that box. The kernel-level A/B says why, in nanoseconds per
call against the scalar body:

| node rows | indexed: split / AVX2 | contiguous: split / AVX2 |
|---|---|---|
| 8 | 0.96x / 0.85x | 1.09x / 1.18x |
| 32 | 1.16x / 0.94x | 1.52x / 1.86x |
| 128 | 0.93x / 0.92x | 1.24x / 2.15x |
| 1024 | 0.95x / 1.01x | 1.40x / 1.99x |

The indexed weighted kernel - 35 to 54 percent of a weighted fit on that
box - is gather-bound: neither width beats the scalar body at any node
size, and AVX2 is the slower of the two more often than not, because
packing four scattered doubles costs more than the three adds it saves.
Where wide vectors do win, 2.15x, is the contiguous twin, which the
re-profile puts at 1.5 percent of a fit.

VERDICT: not worth keeping. Under the same rule the layout fork uses -
under 2 percent of a whole fit at every cell - the shipped result stands
and the AVX2 body is deleted; it was a new translation unit and a build
rule on every platform, plus a standing obligation to keep two bodies
byte-identical forever, for a measured zero. It is recoverable from
813630e6 with its test if the premise ever changes. What survives is the
level-invariance assertion, which now pins the constraint positively:
these entry points are selected by build mode and by nothing else.

S3, the run loop (dec-B88):

7. Replace the sleep loop in [`Sampler::run`](../../src/bartcore/sampler.hpp)
   with a mutex and condition variable: a worker takes the mutex, does
   the final `numChainsRunning` decrement under it, then notifies, so
   the caller cannot miss the wake and the latency assertion is not
   racy; the caller waits with a 100 ms timeout, and a timeout wake does
   what the tick does today. The inline single-chain branch is untouched.
8. tests/cpp: extend [`testRunCancellation`](../../tests/cpp/test_sampler.cpp)
   with a multi-worker arm asserting a short run returns well under 100
   ms and a cancel still lands, plus a verbose-sink arm.

S4, within-chain threading (dec-B89):

9. Measure first, before revival code lands: build the archived
   prototype's two files onto a scratch branch off the tip into a
   private library; time the sampler loop only (`sampler$run`, ingestion
   excluded), single chain, friedman at n in {1e5, 1e6}, m = 75, at 1, 2,
   4 and 8 threads - 8 is where both recorded real-engine runs lose, so
   not optional - on the quiet machine, alternating rounds against a tip
   build, per-round minima. Record msec/iteration and the ratio per cell,
   the machine and its load, and whether draws stayed byte-identical
   across thread counts; no microbench number is admissible. The result
   becomes a new section of docs/design/within-chain-threading.md, whose
   Status line then says the verdict stands and the opt-in ships under
   the regression rule.
10. Revive: rebase the block scheme and wcpool.hpp onto the tip,
    re-deriving the gather blocking over the fused pass, keeping the pool
    parked between sweeps and the size cutoff unchanged. The scatter half
    stays out unless step 9 measures it positive.
11. R surface per fork 4: no new argument. The bridge passes the
    budget through; the sampler runs `min(n.threads, n.chains)` chains
    at once and gives each chain `n.threads %/% n.chains` within-chain
    workers (1 means serial, the law unchanged), so the surplus is
    spent as 0.9-34 spent it. `dbartsControl`'s default stays the core
    count. The memory audit's number, once that plan reports, decides
    only the manual's recommendation, not a default. One test asserts
    that a budget at or below the chain count takes the serial path.
12. Tests: byte-identical draws at 2, 4 and 8 within-chain workers at
    a fixed seed; a budget at or below the chain count reproduces the
    serial draws; a below-cutoff run takes the serial path.
13. Manual: the control page states that `n.threads` is a budget,
    that a budget above `n.chains` turns within-chain threading on, what
    it measured (0.9-34 about 10 percent at four threads on
    one chain; the barrier prototype 12 percent at best, and SLOWER than
    serial at eight threads on both hosts measured), that multi-chain
    parallelism is the effective use of cores, and (fork 5) that
    enabling it changes draws.

S5, the constants audit (dec-B91, B93, B110):

14. Document each at its definition, with origin and limit. dec-B91
    names the first seven rows; max.leaf.size joins under dec-B110, and
    the row cap under front-door's family-object slice, which moves it
    onto `hazard()` in R/model.R - S5 cites that home.

    | constant | file | what it limits | measurable |
    |---|---|---|---|
    | [`maxNumCutsRepresentable`](../../src/bartcore/data.hpp), [`maxCategories`](../../src/bartcore/data.hpp), [`maxLevelsForKind`](../../src/bartcore/data.hpp) | data.hpp, over [`xint_t`](../../src/bartcore/data.hpp) and its `XINT_TYPE` mirror | 65533 cuts, 65535 categorical levels, 65534 for an ordered factor, which spends one code per cut | what breaks at and above each cap and whether a real design reaches one; widening the code is out of scope |
    | [`categoricalExhaustiveCap`](../../src/bartcore/scan.hpp) | scan.hpp | above ten present levels the exact partition enumeration degrades to prefix splits | time and acceptance of exact enumeration versus the scan path at 8, 10, 12, 14 present levels |
    | [`LinearGaussianLeaf::maxNumCovariates`](../../src/bartcore/model.hpp) | model.hpp, refused in facade.hpp | a ninth leaf-regression column | fit time and leaf conditioning at 4, 8, 12, 16 designated columns |
    | [`perturbWidth`](../../src/bartcore/moves.hpp) | moves.hpp | the perturb proposal moves one grid position | acceptance and ESS at widths 1, 2, 4 on the standard designs |
    | [`testFitParallelCutoff`](../../src/bartcore/chain.hpp) | chain.hpp | test fits below 65536 rows stay serial | threaded versus serial wall time across n.test |
    | [`predictParallelCutoff`](../../src/bartcore/sampler.hpp) | sampler.hpp | predict below 1e7 cells stays serial | threaded versus serial wall time across rows x trees x draws; calibrated per dec-B93 |
    | [`sparseDensityThreshold`](../../src/bartcore/data.hpp) | data.hpp | a column below 0.2 density is stored sparse | memory and gather time either side of 0.2 |
    | [`maxLeafSize_`](../../src/bartcore/model.hpp), [`gp`](../../R/model.R)'s max.leaf.size (dec-B110) | model.hpp, R/model.R | a GP leaf above 256 observations falls back to a constant leaf | fit time and fallback share at 128, 256, 512, 1024 |
    | the person-period row cap, today [`expandDiscreteTimeHazard`](../../R/dbarts.R)'s max.rows | R/dbarts.R, moving to `hazard()` in R/model.R | the expansion refuses above 1e7 rows | expansion time and peak memory at 1e6, 1e7, 3e7 rows |

15. One measurement script per measurable row under benchmarks/R, each
    printing a table and runnable in a quick mode, run on the quiet
    machine. Record every result in a new docs/design/engine-constants.md
    - the note dec-B91 owes for a moved default - stating the value, the
    measurement and whether it binds.
16. Expose as `dbartsControl` arguments, with validity checks and bridge
    plumbing, neutral at today's values: the categorical enumeration cap,
    the test-fit cutoff, the predict cutoff (reusing the
    `cutoffOverride` seam), plus any row step 15 shows binding. Each gets
    a test that the setting reaches the engine and the default reproduces
    today's draws. The predict cutoff is calibrated, not left at 1e7.
17. GP fallback counter and warning (dec-B110): count leaf evaluations
    that took the constant-leaf fallback at
    [`maxLeafSize_`](../../src/bartcore/model.hpp)'s four call sites,
    report the share through the fit object, and warn above the fork-6
    threshold, share included. One test fires it, one does not.
18. Manual: one control-page table listing every constant, its default,
    recommended range and whether it is settable; the GP page keeps its
    10 to 25 tree guidance and gains the warning's meaning; NEWS carries
    the settings and the warning.

## Verification

```
# per slice, from the slice's own worktree and private library
R CMD INSTALL --preclean -l <lib> .                      # shipped build
R CMD INSTALL --preclean -l <ref> . --configure-args=--enable-reference-build
cd tests/cpp && make && ./test_bartcore
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
R_LIBS=<ref> Rscript -e 'tinytest::test_package("dbarts")'
  # per snapshot file, by count: length(run_test_file(f)) is 0 under <lib>
  # (exit_file leaves no skip attribute), its recorded count under <ref>
R_LIBS=<ref> Rscript benchmarks/R/equivalence.R compare \
  benchmarks/baselines/equivalence-deb144d2.rds --strict-coverage
  # 52 "identical draws (same RNG stream)", no "max |z|"; bcf 12, multi 11
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare \
  benchmarks/baselines/equivalence-deb144d2.rds
  # shipped, S2 onward: statistical mode everywhere, none at |z| > 4
for g in bd-balance change-balance perturb-balance backfit-exact linear-exact \
  categorical-exact heteroscedastic-exact multinomial-exact hazard-exact \
  bcf-exact bcf-exact-weak bcf-exact-restricted
do R_LIBS=<lib> Rscript benchmarks/R/$g.R quick; done   # exact gates, quick
R_LIBS=<lib> Rscript -e 'dbarts:::buildInfo()'   # shipped; <ref> reads reference
Rscript tools/check-win-drift.R && Rscript tools/check-doc-freshness.R
cd tests/cpp && make clean && make OPT="-O2 -g -fsanitize=address,undefined" && \
  ASAN_OPTIONS=detect_container_overflow=0 ./test_bartcore
```

Maintainer-run, quiet machine, never concurrent: `Rscript
benchmarks/R/bench-sampler.R compare
benchmarks/baselines/bench-sampler-127f04ee.csv`. Expected: S1 both
builds bitwise on all three baselines, test-simd.R green on both; S2 the
shipped build no slower on any bench arm, the reference build unchanged,
test-simd.R still bitwise across every dispatch level; S3 a four-chain
single-sweep run under 100 ms, every gate identical; S4 byte-identical
draws across within-chain thread counts, the default bitwise against the
tip, the measured ratios in the design note whatever they say; S5 every
default unchanged, each setting reaching the engine, the GP warning
firing only on a degenerate fit.

## Landing note, S1 (2026-09-09)

LANDED at 167e2c62dff4576b35134518d7bdc026c2baf8db, six commits:

- 7641bdd5d3754d3ecba4b633229797ee5ed2c58b Add the --enable-reference-build configure flag and its Windows route
- eb7b137a08dff3a600c8c4486152fb3aee314131 Report the build mode, compiled ISA set and dispatch level to R
- d748f6e1cdbc34245a441b025c276f5be38f7a82 Key the seeded-drift snapshots on the reference build
- 6b4474b288cea158bc3481fb2dcfac130cb44867 Split the CI gates by build mode
- 7ec2aa58b96f2b99952dbe1407442f32c09da67e Assert the mode on each CI reference install
- 167e2c62dff4576b35134518d7bdc026c2baf8db Refuse snapshot regeneration off the reference build, and match the CI notes

`--enable-reference-build` in configure.ac defines
`DBARTS_REFERENCE_BUILD` in src/config.hpp and src/misc/config.h alike,
regenerated with autoconf 2.73, byte-for-byte reproducible; Windows
routes through an environment variable of the same name read by
src/Makevars.win, joining
[check-win-drift.R](../../tools/check-win-drift.R)'s expected-absent
table. The bridge entry `dbarts_buildInfo`, reached as
`dbarts:::buildInfo()`, reports mode, compiled ISA set and dispatch
level from the source `simd.c` uses. The four snapshot files exit on
the shipped build and run in full on the reference build (14, 3, 7, 3
assertions), asserted by count in
[test-build-info.R](../../inst/tinytest/test-build-info.R);
regenerate-snapshots.R refuses a shipped build. cpp-tests.yaml gains a
reference arm running the snapshot files and the three equivalence
compares (benchmarks/** no longer ignored, timeout 60, no skip of the
reference install on a C++ failure); exact-gates.yaml runs its two
cross-host compares on a reference install, a mode assertion after
each. NEWS's 1.0-0 entry and README's [CI](README.md#ci) section gain
the reference arm. Real diff: 18 files, +396/-15 (362 insertions
excluding the regenerated configure), inside budget.

Gates, independently on both builds: equivalence 50/12/11 identical, 0
skipped, no "max |z|" line either build - the inertness proof; tinytest
7859/0 shipped (snapshots exit), 7886/0 reference; tests/cpp all
passed; check-win-drift.R, check-doc-freshness.R and check-rc-codoc.R
exit 0; `R CMD check --as-cran` OK from a clean tarball; NEWS parses;
`air format --check` clean; regenerate-snapshots.R on the reference
build rewrites the four files byte-identically. Mutation probe: a
reference build forced to report "shipped" fails the CI mode assertion.

Review findings fixed before landing: the regeneration tool did not
stop on a shipped build; README's CI list was stale for cpp-tests; the
timeout raised from 45 to 60; the reference install no longer skipped
after a C++ failure. Remaining: S2, the vector suffstat kernel; S3, S4, S5.

x86 leg, run 2026-09-09 at tip 777d0f1c on an x86-64 Linux host
(AVX2/FMA, dispatch level 8, R 4.6.1) after this slice and S3 landed:
both builds install; tests/cpp 283/283 plain, under ASan/UBSan and
under ThreadSanitizer with zero diagnostics; tinytest 0 failures on
both builds (8182 shipped, 8209 reference); the four RNG-locked
snapshot files RUN AND PASS in full on the reference build there, so
the reference build reproduces the arm64-recorded draws bitwise across
ISA; and the equivalence trio against the arm64 baselines reports
max |z| = 0.00 on all 50 gaussian scenarios and tier-1 PASS on all 12
BCF and 11 multinomial scenarios under the cross-host compare, on the
shipped build as well as the reference build. Until S2 lands the two
builds are the same kernel, so this is the cross-ISA baseline S2's
layout rule (Decision 2) is measured against.

## Landing note, S3 (2026-09-09)

LANDED at 7030557578a9a6dd225ed8abc23d1dc654a531ff, three commits:

- 6c6582bf889438d65b67c9c3dc3206f0afe96cee Wait on chain completion in the multi-worker run loop
- 1e895684169badbd87f1cb049f632b0ad4d1dd26 Cover the run loop's latency, cancel and flush paths
- 7030557578a9a6dd225ed8abc23d1dc654a531ff Name the predicate re-check as the wait's guarantee, and size the verbose arm to the host

[`Sampler::run`](../../src/bartcore/sampler.hpp)'s multi-worker branch waits
on a condition variable instead of sleeping in 100ms ticks: a worker takes
the mutex, does the final `numChainsRunning` decrement under it and
notifies; the caller waits with the predicate overload of `wait_for` and a
100ms timeout, so a completion landing while the lock is released for the
flush and the interrupt poll is seen rather than slept through. A timeout
wake still does what the tick did - flush, then poll - except the flush now
precedes the poll, a benign reorder the plan did not call for, stated here
because it can surface a queued line up to 100ms sooner. The mutex and
condition variable are locals of `run()`, joined with every worker before
going out of scope; the inline single-chain branch is untouched and the
engine stays R-agnostic. [`testRunCancellation`](../../tests/cpp/test_sampler.cpp)
gains a latency arm (a short four-chain run returns well under 50ms), a
cancel arm and a verbose-sink arm that grows its run until it clears 100ms
on the host; [Threading model](../architecture.md#threading-model) describes
the wait.

Measured: 200 calls of a 4-chain, 4-thread `run(0, 1)`: median 105ms
before, 0.12ms after; the reviewer measured 0.09-0.20ms per call, a 250x
margin under the bound, 2.4ms under ThreadSanitizer. Real diff:
sampler.hpp +44/-12, test_sampler.cpp +104/-10, architecture.md +8/-4; the
test line runs about 2.5x its ~40-line budget, inside the slice's 1.5x stop.

Gates, independently: tests/cpp all passed three times in a row,
ThreadSanitizer 0 diagnostics, ASan/UBSan 0 diagnostics, 25 repeats of the
sampler suite 0 failures; tinytest 8069/0 at its base, 8105/0 on the merged
tree with front-door S3; equivalence 50/12/11 identical, 0 skipped, no
"max |z|" line (gaussian against equivalence-c42b72af.rds on the merged
tree); `R CMD check --as-cran` OK from a clean tarball; doc-freshness,
rc-codoc and check-win-drift exit 0; `air format --check` clean. Mutation:
removing the notify fails the latency and cancel arms at about 101ms,
removing the in-loop flush fails the verbose arm.

Review findings fixed before landing: the comment above the mutex claimed
the count could not reach zero between the test and the wait, false as
written, and was reworded to name the predicate re-check as the guarantee;
the verbose arm's fixed sample count had only 3.6x headroom on this host
and now grows until the run clears a timeout - later replaced outright by
a deterministic split rather than a clock-sized run, cb999b60,
[Follow-up, S3 (2026-09-09)](#follow-up-s3-2026-09-09). Remaining: S2, the
vector suffstat kernel, still open; S4, S5.

## Follow-up, S3 (2026-09-09)

The verbose-sink arm above was flaky, and failed once on the cpp-tests runner
on a tip that changed no engine file. Its premise measured the whole `run`
call - setup, the post-join flush and the terminal summary included - against
100ms, so a run whose chains all finished inside the first wait could still
clear the bound with nothing flushed before the join, and the assertion then
read the flag at false. Reproduced at 7 failures in 25 runs by pinning the
sample count so the call lands just past 100ms.

The arm no longer reads the clock. It splits in two. A short run that
completes checks that every chain's lines reach the console. A run asking for
far more samples than any host draws in the time the arm takes checks the
timeout flush: the poll records the capture's length at the first poll, which
follows the loop's first flush, and any later poll seeing a new "iteration:"
line has seen one flushed on a later pass - reachable only through a
`wait_for` timeout, since the predicate re-tests the count under the lock and
the chains are nowhere near done. That poll cancels the run, so the arm costs
one tick, and load can only delay the pass, never skip it. The latency and
cancel arms' bound moves from 50ms to 75ms: a reintroduced tick wait costs its
full 100ms whatever the host does, so any bound under 100ms discriminates the
same, while an honest return - worst seen 18ms, under ThreadSanitizer with the
box oversubscribed - gains margin against CI noise.

Gates: tests/cpp 283 ok, 0 failures; 50 repeats of the sampler suite under 4
spinners and 50 more under 16, 0 failures; ThreadSanitizer 0 diagnostics on
the full suite and over 50 sampler-suite repeats, 0 failures; ASan/UBSan 0
diagnostics on the full suite. Mutation: widening
the wait's timeout to 100 seconds fails both of the new checks. tests/cpp
only, no engine file touched, so the equivalence and speed compares are not
owed and `air format --check` is n/a.
