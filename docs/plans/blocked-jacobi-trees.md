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

## Phase 0 results (2026-07-21): EXACT, no per-sweep tax at realistic m -> GO Phase 1
Method: an ORACLE HARNESS - m single-tree dbarts bartcore handles (the real tested
tree kernel: grow/prune/change/swap + leaf Gibbs), backfit from R; the only novel
code is the R-level noise split + schedule + sigma draw. Serial-oracle validated
against bart2 (posterior-mean fit corr 0.996, recovers f as well, sensible sigma).
Throwaway scripts: $CLAUDE_JOB_DIR/tmp/bj-*.R (job b073bb28). Not landed; a
measurement of ESS-PER-SWEEP, not wall-clock.

EXACTNESS: PASS, by four INDEPENDENT lines (not one derivation):
1. Closed-form MVN match, constant stumps (0a) and fixed multi-leaf structure
   (0a-plus), b=1,2,4: posterior mean/var/cov within ~1.5 MCSE of the analytic
   posterior; b=1 noise-split reduces to serial exactly (max|z-R|~1e-15).
2. Closed-form match b=1,2,4,8 fixed-structure (isolates BIAS from mixing - no
   structure moves, converges fast even at b=8): every b within ~2.4 MCSE, var
   ratios ~1 -> NO b-growing bias, even at b=8.
3. An INDEPENDENT blind derivation (Opus agent, given only the model + mechanism,
   not my algebra): verdict EXACT via the conservation law sum_{k in B} v_k =
   sigma^2 (v_k = sigma^2/|B| meets it); structure moves provably safe (the
   marginalization over z_B uses g_k(theta_k) only as a value, agnostic to
   constant-leaf vs full tree); the standard sigma full-conditional is correct
   (z's are auxiliary, discarded, not in the Markov state); and it enumerated the
   off-spec bias hazards (b not dividing m with a rigid sigma^2/b; state-dependent
   schedule; pinned trees written inside a barrier; leaf prior recalibrated to the
   reduced noise) - ALL of which the harness avoids (v = sigma^2/|B| by actual
   batch size, fixed round-robin, pinned frozen, leaf prior fixed and empirically
   sigma-independent).
4. Real-tree b=2 vs serial within the serial-vs-serial MC-noise control (drawn
   sigma, m=30).

TAX: the per-sweep MIXING tax at realistic m is MILD (fit retains ~85-100% of
serial's per-sweep mixing). Warm-start gold-standard (all methods from an identical
converged state via storeState, so neither non-convergence nor the stuck-chain-
inflates-ESS artifact can bias it), m=75 n=1e4: fmin (worst pointwise-f coord)
ESS/sweep ratio 0.94-1.13 across b=2,4,8; global-fit ratio 1.06-1.47 (blocked mixes
BETTER - noise-injection DECOUPLES the trees, offsetting the precision penalty).
Consistent across cold seed30, warm seed30, and an independent replication (data +
tree seed) seed77 (fmin 0.93-0.99). CONFIRMED at m=200 n=4000 (the realistic BART
size): fmin ratio 0.86-0.95, gmean 0.80-0.89, break-even (ESSratio*b) 1.8-6.9; and
b=8 CONVERGES FROM COLD at m=200 (sigma windows flat ~1.006 vs serial 0.998) - the
dilution is stronger at larger m, so the tax and the burn-in both ease as m grows.
Across m=75/200 x b=2/4/8, break-even (ESSratio * b) = 1.8-8.5, always well above 1.
sigma ESS ratio is noisier (0.49-2.22 across runs) but non-binding (sigma mixes
fast, ESS ~1e3-3e3; the fit governs BART). The conservation-law precision tax is
REAL but manifests as a BURN-IN/transient cost (theorist: learning rate down 1/b,
structure SNR down sqrt(b)), NOT a stationary per-sweep cost: a COLD b=8 chain at
SMALL m (30) climbs out of underfit only over ~40k+ sweeps (sigma 1.22 -> 0.94),
while WARM b=8, or cold b=8 at m=200, mixes fine. Practical: b=2-4 fine from cold at
any m; b>=8 wants a warm start only at small m.

VERDICT: GO to Phase 1. Neither kill condition tripped (exactness holds; the
ESS/sweep tax does not block - it is ~zero at realistic m). Phase 1 decides the
REMAINING and decisive question, which Phase 0 cannot: does the b-fold parallel
WALL-CLOCK speedup survive the overheads - m/b barriers per sweep (the wcpool
substrate; amortizes only when per-tree work is large, i.e. large n) AND the noise
split's own O(n*m)/sweep cost (roughly doubles per-sweep field work, but
parallelizes)? The straight-threading NO-GO (net 1.67x, docs/design/within-chain-
threading.md) is the bar to beat. Caveats to carry into Phase 1: ESS/sweep is not
ESS/sec; only Friedman/n=1e4 tested (ratio is n-robust statistically; add a second
DGP); sigma mixing mildly degraded (non-binding). Recommended de-risk before the
full engine kernel: a cheap wall-clock probe of the barrier + noise-split overhead
vs the b-fold speedup, since that overhead is exactly what killed straight
threading.
