# blocked-jacobi-trees: implementation plan (falsify-first)

Design lives in docs/design/parallel-bart-frontier.md sec 3.5 (the noise-split
movement-budget kernel) + sec 1-2 (why exact tree-parallelism is hard, the
conservation law) and docs/design/within-chain-threading.md sec 7 (the ESS/sec
protocol) + sec 8-9 (the straight-mechanism NO-GO + the fp32 re-check). This doc
is the ORDERED experiment plan. Status: UNEVALUATED, greenlit to try 2026-07-20.

## What it is (one paragraph)
Single-chain parallelism. Conditional on structures, the forest is a tiny joint
Gaussian; backfitting is Gauss-Seidel with one color per tree, and at default
depth the cross-tree leaf graph admits no coloring with < m colors, so exact
tree-parallelism is unavailable. Noise-splitting augmentation MANUFACTURES it:
draw per-tree pseudo-responses y_k ~ N(g_k, sigma^2/b) summing to the batch
residual, so b trees become conditionally independent and update in PARALLEL
(m/b barriers/sweep instead of ~3m; each worker a cache-coherent single-tree
update). The precision budget is CONSERVED (a law, sec 2): the b released trees
share sigma^2, so their structure moves see precision b/sigma^2 -> mixing
DEGRADES with b (the tax). Sec 3.5's refinement: the variance ALLOCATION is free,
so RELEASE a rotating handful at near-full variance and PIN the rest - exact for
any STATE-INDEPENDENT rotation schedule; the uniform tax becomes an aimable
movement-budget. Arbiter is ESS-PER-SECOND (ESS/sweep falls from the tax, sweep
wall-clock falls from parallelism).

## Phase 0 - SERIAL kernel, statistics only (cheapest, decisive; DO FIRST)
NO threading, NO wcpool. Implement the noise-split movement-budget kernel
single-threaded: each sweep, pick the released batch by a state-independent
rotation (e.g. b=2 of m round-robin), draw the sum-constrained pseudo-responses,
update the released trees against THEIR pseudo-residual and the pinned trees
classically (or hold pinned), advance the rotation. Two gates:
1. EXACTNESS (the make-or-break): does the kernel target the SAME posterior as
   serial BART? The design claims exact for a state-independent schedule - VERIFY,
   do not assume. Use SBC / an exact-posterior check (compare posterior summaries
   - f, sigma, coverage - vs serial BART over many seeds on a small problem; and
   an unbiased-coupling or rank-uniformity check if feasible). A mis-derived
   augmentation silently shifts the posterior -> this is where it most likely dies.
2. ESS/SWEEP TAX: how much does mixing degrade vs serial backfitting at b=2 (then
   b=8)? Measure integrated autocorrelation / ESS per sweep on f and sigma,
   signal-concentrated synthetic (friedman) at n in {1e4,1e5}, m in {75,200}.
KILL if exactness fails, or if the ESS/sweep tax is so large no plausible
parallel speedup (<= b-fold) recovers it (grouped-mixing precedent: a corner
that looked huge was estimator-unstable and not worth a heavy kernel). Phase 0
needs only R-level or a minimal engine hack; it decides go/no-go before any
systems work.

## Phase 1 - PARALLEL, ESS/sec + head-to-head (only if Phase 0 passes)
Wire the kernel to the worker pool. SUBSTRATE TO REUSE:
`git show origin/archive/within-chain-threading:src/bartcore/wcpool.hpp` -
`WithinChainPool`, a general persistent std::barrier pool with `forRange(total,
fn)` running fn(begin,end) per worker (owner participates, workers park between
regions, never call into R). The straight prototype (same branch, commit
54a60aa, +206 lines in chain.hpp) shows the integration and ALREADY MEASURED the
per-barrier overhead - real data going in, since blocked-jacobi's whole bet is
m/b barriers beating the straight mechanism's ~3m. Dispatch INDEPENDENT tree
updates via forRange(b, updateReleasedTree) instead of the straight prototype's
reduction slices. Measure ESS/SEC vs serial, and the sec-7 head-to-head vs the
straight mechanism (x86 box, friedman n in {1e4,1e5}, m in {75,200}, threads
{2,4,8}). Whichever wins ships; the loser is recorded closed.

## Gates & process
- Phase 0 is a throwaway measurement (like the fp32 falsification) - no blind
  critique needed to RUN it, but its EXACTNESS finding must be rock-solid before
  believing any speed number.
- If Phase 0 passes -> design-first + INDEPENDENT BLIND critique of the exactness
  argument and the Phase-1 kernel BEFORE landing (it is a NEW EXACT KERNEL, not
  the shifting/equivalence class: needs exact-posterior gates at b in {2,8}, not
  just the bitwise trio). One re-record is NOT enough - the kernel changes the
  transition operator, so the gate is statistical (SBC/coupling), like fp32 but
  stricter (the posterior claim is the whole point).
- Model tiers: Opus for the kernel/numerics/exactness; the ESS machinery is
  numerics too. bench/ESS on the quiet x86 box.

## Prior art (checked 2026-07-20): none direct
Existing parallel BART is DATA-parallel (Pratola 2014 / OpenBT: partition obs,
reduce suffstats - the orthogonal axis dbarts already threads; bartMachine:
parallel chains; XBART: different grow-from-root algorithm). The noise-split
movement-budget tree kernel appears novel to this program; lineage is general
auxiliary-variable / augmented parallel-Gibbs. A literature/code sweep before
committing is optional insurance, not done.
