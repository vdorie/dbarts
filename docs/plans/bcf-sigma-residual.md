# BCF sigma burn transient: anatomy, candidate remedies, recommendation

Read-only design memo for TODO item `bcf-sigma-residual`. Repo dbarts @ bartcore
6944811 (a-move / interweaveGlueRidge LANDED). All numbers measured against the
installed shared library via benchmarks/R/sbc.R. Scratch, scripts, and rds
output were run out-of-repo and are not preserved.

## 0. TL;DR

The sigma SBC flag is a slow FOREST-MIXING transient, not a glue-scale
init problem. The memo's routed remedy (amplitude-aware a/mu init) is
INEFFECTIVE (measured); so is a warm start. The bias decays only with more
BURN. Winner: candidate (c), longer burn (+ a doc note). The complete fix for
the frozen extreme-|a0| tail is faster tree-structure mixing at high SNR, an
engine change (separate item), NOT another glue move.

## 1. Anatomy of the transient -- where the time goes

Init path (cited): BCFState starts a=1, b0=0, b1=1, aVariance=1
(chain.hpp:329-330); both forests start all-zero fits (buildBCFForest,
chain.hpp:2095-2096); sigma_ = initialSigma() = sigest/range (chain.hpp:468,
model.hpp:1908). Across SBC reps the sampler is reused: sampleTreesFromPrior /
sampleNodeParametersFromPrior redraw only the forests (chain.hpp:899-957) and
setResponse(updateScale=FALSE) reuses the fixed build scale (model.hpp:1936-1946)
-- so a, aVariance, b0, b1 and sigma CARRY OVER between reps; only the mu/tau
forests are re-drawn shallow-from-prior each rep.

Traced a, ||mu|| (sd of mu total), a*mu amplitude (sd), sigma over 20k-40k
sweeps on strong-|a0| reps (characterize.R, characterize2.R):

- a*mu AMPLITUDE is correct by sweep ~10 for every a0 (a0=6: amuSd 1.13@sw10 vs
  target 1.20; a0=14: 1.93 vs 1.98; a0=100: 12.3 vs 13.2). The scalar a jumps to
  ~a0 in one conjugate drawGlue (chain.hpp:2147-2157); the a-ridge move then just
  trades a<->||mu|| at fixed product. So glue and leaf MAGNITUDE are FAST.
- SIGMA is the slow variable: it starts high (unexplained signal) and decays over
  tens of thousands of sweeps. sigFirst ~= sigLast within an 18k window
  (instrument.R) -- within one window it looks stationary, i.e. the decay is
  slower than the sampling window.
- The residual that sigma absorbs is the mu-forest SHAPE misfit: a fresh prior
  forest has the right amplitude but the wrong partition, so y - a*mu - b_z*tau
  keeps a chunk of signal proportional to |a0|. Time goes into the mu forest's
  Metropolis tree moves slowly discovering mu0's split STRUCTURE.
- It is rate-limited by SNR: injecting a large a (a-init) raises the mu-forest
  weight a^2 -> higher effective precision -> LOWER structural-move acceptance ->
  structure freezes FASTER and bias is slightly WORSE mid-burn. Direct evidence
  the bottleneck is tree structure, not scale.

Ground truth, real reused-sampler harness, R=80, burn=150, thin=120
(instrument.R): sigma rank mean 65.1/75 (reproduces the flag). Bias =
mean(recorded sigma)/sig0 by |a0|:

    |a0| in (0,2]    n=45  bias 0.999  rankMean 77.0
    |a0| in (2,5]    n=16  bias 1.027  rankMean 67.6
    |a0| in (5,10]   n=10  bias 1.159  rankMean 36.9
    |a0| in (10,20]  n= 6  bias 1.188  rankMean 41.2
    |a0| in (20,Inf] n= 3  bias 1.217  rankMean 15.7

The flag is entirely the |a0|>5 minority (~24% of reps, Cauchy(0,2) tail); the
bulk are perfectly calibrated. Bias grows with |a0| and, when sig0 is small,
pushes sigma above all L draws -> rank ~0.

## 2. Candidate remedies (measured: bias = mean window sigma / sig0 vs burn)

burncurve.R, 10 representative strong reps (|a0| in 5-25), cumulative burn
sweeps (burn units at thin=120 in parens):

    burn sweeps:  2000    9000   18000   36000   72000
                  (~17u)  (75u)  (150u)  (300u)  (600u)
    COLD (a=1):   2.75    2.20   1.52    1.26    1.07
    (a) A-INIT:   2.40    1.90   1.66    1.21    1.13
    (b) WARM:     1.73    1.76   1.44    1.15    1.21

- (a) data-informed a-init (inject a = sd(y0)/sd(yBuild) via storeState/setState;
  an always-legal, data-dependent init): NO asymptotic benefit, slightly WORSE at
  the standard burn (1.66 vs 1.52) -- raising a raises mu SNR and freezes
  structure. REJECT.
- (b) short warm start (grow mu W sweeps before freeing a): only front-loads a
  little (1.73 vs 2.75 at 2k) and has NO asymptotic gain (1.21 vs 1.07 at 72k).
  REJECT as a fix.
- (c) longer burn: the ONLY lever that works. Time-to-settle (mean bias -> ~1.05):
  ~72k sweeps = burn ~600 units at thin=120. burn=150 leaves mean bias 1.52 on
  the strong minority; burn=300 (the memo's burn-2x) still 1.26; burn=600 ~1.07.

Per-rep, settle time tracks the rep's SNR (|a0| x amplitude / sig0), not |a0|
alone: reps pairing strong signal with small sig0 are slowest, and 3/10 strong
reps stay >15% biased even at 72k sweeps. Extreme tail (a0=40, 100,
characterize2.R): sigma plateaus ~5x high with NO decay through 40k sweeps --
frozen structure. Cauchy(0,2) tail mass: P(|a0|>5)=0.24, P(|a0|>20)=0.063,
so ~6% of SBC reps are effectively frozen and contribute rank ~0 each; that
alone puts ~0.05-0.06 excess in the low-rank ecdf -- the R=400 band is ~0.066,
so even burn=600 sits AT the edge by construction.

USER-RELEVANCE REFRAME: |a| >> 1 means prognostic signal >> the response scale
the sampler was BUILT with. With the data scale set from the actual y (normal
use), a is O(1) and the transient is mild; the strong-|a0| regime is real only
when the build scale is stale relative to a swapped-in response
(setResponse(updateScale=FALSE) inside a larger Gibbs sampler -- the
dbartsSampler use case) or in SBC's own prior tail. The doc note should say
exactly that.

Per-candidate draw-change class and smallest diff (for the record; (a),(b)
rejected on measurement):
- (a) a-init: opt-in possible (default a=1 preserved). Smallest diff: BCFSpec
  gains aInit (chain.hpp:314), BCFState init reads it (chain.hpp:482-487),
  bridge bcfParams grows to 9 (R_interface_bartcore.cpp:1202,1228-1241),
  bartcoreBCFSampler gains a.init (R/bartcore.R:419). Data-dependent init is
  MCMC-legal; reproducibility rides the existing seed (estimator sd(y)/sd
  (yBuild) is deterministic). All existing draws unchanged if default a=1.
  Prototyped read-only via storeState/setState (bcf slot, probe2.R).
- (b) warm start: needs NO engine change at all -- reachable today R-side:
  build update.a=FALSE, run W sweeps, storeState, setState into an
  update.a=TRUE sampler (burncurve.R does exactly this). As a shipped feature
  it would be a run-schedule wrapper; opt-in, changes no default draws.
- (c) longer burn: changes NO draw law -- more of the same chain. Smallest
  diff: benchmarks/R/sbc.R -- runSbcBCF's burn default (line 938, currently
  50L) and/or the isBCF CLI branch (line 1382, which passes no burn); pin
  ABSOLUTE sweeps, e.g. burn = ceiling(72000 / thin), since burn units are
  thinned (chain.hpp:642 totalIterations = (burn + samples) * thin). Plus a
  burn-guidance note in docs/design/bcf.md. ZERO engine change.

## 3. Acceptance (winner = longer burn), n=200 glue-on, thin=120, L=150

Real runSbcBCF, burn=600 thinned units (72000 sweeps), chunked 50/seed and
pooled (accept.R + poolreport.R; chunk rds accept-b600s*.rds).

Interim R=200 (seeds 101,202,303,404; band 0.0924):
    sigma  rankMean 75.2  chisqP 0.772  ecdf 0.0463  PASS (inside the R=400
    band ~0.066 already; burn=300 at R=200 was 72.05 / 0.0888)
    abs.a 77.2 PASS, abs.diff 76.8 PASS, prog1-5 71.6-75.6 all PASS,
    eff1,3,4,5 PASS; eff2 0.0951 vs 0.0924 FLAG (chisqP 0.063, marginal --
    R=400 adjudicates); a / b1.minus.b0 FLAG by design (sign-ill-posed).

FULL R=400 (8 x 50, seeds 101..808; band 0.0656; ~5.6 min/chunk, 2 at a time):

    functional   rankMean   chisqP  ecdfDiff   band    verdict
    sigma            74.5    0.913    0.0325   0.0656   PASS   <- the target
    abs.a            77.1    0.145    0.0590   0.0656   PASS
    abs.diff         76.4    0.061    0.0355   0.0656   PASS
    prog1..5      72.7-75.6  >=0.04  <=0.0544  0.0656   PASS (all)
    eff1..5       72.6-75.5  >=0.03  <=0.0468  0.0656   PASS (all)
    a                80.8    0.000    0.2613   0.0656   FLAG (sign-ill-posed by design)
    b1.minus.b0      72.2    0.000    0.2068   0.0656   FLAG (sign-ill-posed by design)

Reference: burn=150 at R=400 was sigma 67.6 / ecdf 0.1003 FLAG (both builds,
a-move landing); burn=300 at R=200 was 72.05 / 0.0888 (would have exceeded the
R=400 band). burn=600 clears with 2x margin. The single residue of the frozen
extreme tail is the lowest bin (28 vs 20 expected, ~= the P(|a0|>40)=3% mass)
-- inside the band: SBC integrates over the prior, and the truly frozen tail
is small enough that burn=600 absorbs it at R=400 resolution.

## 4. Recommendation

1. bcf-sigma-residual: ROUTE to burn/doc, NOT to init and NOT to a ridge move.
   Implementer diff: benchmarks/R/sbc.R only -- runSbcBCF burn default (line
   938) and/or the isBCF CLI branch (line 1382): pin absolute burn sweeps,
   burn = ceiling(72000/thin) (=600 units at thin=120). Plus a burn note in
   docs/design/bcf.md: burn scales with the fitted-|a| regime (rule of thumb
   measured here: |a|~O(1) needs ~2k sweeps; |a|~5-25 needs ~72k; SNR =
   signal/sig0 is the real driver), and the strong-|a| regime arises when the
   build scale is stale vs the swapped-in response (setResponse
   updateScale=FALSE inside a larger Gibbs loop) -- normal single-fit use with
   the scale set from y keeps |a|~O(1). Correct the prior memo's routing: the
   amplitude-aware a-init it proposed is measured INEFFECTIVE (section 2).
   No engine change; no draw-law change; all existing gates untouched.
2. OPTIONAL follow-up (new item, not required for acceptance): the extreme
   tail (|a0| >= ~40, ~3% of the prior) is frozen-structure mixing that no
   burn fixes -- a principled engine remedy would target tree-structure mixing
   at high SNR (e.g. tempered early sweeps or restart-style grow proposals).
   Decoupled from the shelved b-ridge move, which remains a mixing nicety.
