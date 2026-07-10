# architecture-numerical-review

agent: two opus readers (architecture critic; numerical auditor),
       read-only; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: two reports recorded here; findings become fix items,
        design-note updates, or explicit accepted-risk notes

Review 6 of the retrospective program (retrospective-reviews.md),
run before review 5 (the roadmap survey needs sustained web access;
this one is purely local).

## Goal

Fresh-eyes stress of the engine now that it is functionally frozen:
(a) does the facade/concept decomposition actually carry the next
five features without deforming, and (b) where does the arithmetic
break first under extreme but reachable inputs - scale, sample
count, weight ratios, and long single chains.

## Method

Two independent read-only Opus readers, one dimension each.

A. ARCHITECTURE CRITIC: stress the decomposition (facade.hpp
   concept selection; Chain/Forest split; ResponseModel decoration;
   LeafModel concepts constant/vector/function; the mutation
   surface; the C bridge) against the five likely-next features,
   each already specified in TODO/design notes:
   1. grow-from-root-warm-start (a second tree-construction kernel
      beside metropolisJumpForTree, sharing a cut-scan primitive);
   2. multi-forest-models (multinomial/heteroscedastic/hurdle on
      the BCF two-forest template);
   3. negative-binomial + ordinal-outcomes (new response families:
      PG-with-r-update, cumulative probit with a cutpoint block);
   4. data-ownership (owned predictor storage replacing borrowed
      REAL(x), per-column widths);
   5. the parallel seam (blocked-jacobi trees / gpu-bart backend).
   For each: what bends, what breaks, what is already provisioned.
   Verdict per feature: FITS (seams named) / STRAINS (what deforms
   and the smallest re-cut) / BLOCKED (what must change first).
   Known couplings to weigh: BCF is constant-leaf-only by
   static_assert; grouped and BCF are mutually exclusive by
   construction; the empty-leaf veto is the single surrogate target
   all moves share; the proposal-equals-prior identity underpins
   birth/death cancellation.

B. NUMERICAL AUDITOR: where does the arithmetic degrade first.
   Seed list (verify, quantify the onset threshold, and judge
   severity; add what you find beyond it):
   - the constant-leaf centered sum-of-squares catastrophic
     cancellation under a large shared residual offset (block-5
     fragility note; score side only);
   - totalFits telescoped rounding drift in a pure run (~1e-12
     relative at 1e5 sweeps measured; single-chain n >= 1e5 is a
     standing workload - where does it matter first);
   - the GP path: 1e-6 conditioning jitter vs the nugget; the
     predict-vs-stored-fits discrepancy (~2e-3 at training rows -
     re-kriging through the cached factorization is not the
     recorded fit; review-4 finding) - defect, doc note, or fix;
   - the grouped tau slice sampler: hangs under extreme ranef
     scale (review-3 poison side effect; review 4 found the
     half-Cauchy tail makes prior draws stall it) - bound the
     stall region, propose the cheapest guard (step-out cap,
     bracket limit);
   - the 1e-9 BCF multiplier floor (block-6 note: magic constant
     unrelated to data scale);
   - log-space robustness of the acceptance ratios (change move's
     counted-ratio logs, chi-k gamma draws at extreme shapes -
     review 3 saw k run away to 1e153 under a poison; is the
     UNPOISONED k feedback loop bounded for extreme nu/M);
   - range scaling at extreme response scales (y ~ 1e15 or 1e-15:
     where do fitScale/fitShift/sigma round-trips lose precision);
   - setState occupancy validation accepting an n-mismatched
     gaussian state (block-7 note) - quantify what a mismatch
     corrupts.
   For each: REAL (onset threshold + severity + smallest fix),
   ACCEPTED RISK (why it cannot bite at reachable inputs), or
   DOC NOTE. Small numerical experiments in R against the
   installed package are encouraged; no engine edits.

## Constraints

- Findings only; the tree does not change under this review.
- Readers are read-only, FOREGROUND only, and spawn nothing.
- docs/design/*.md may be read for INTENT here (unlike the
  correctness audit, this review judges design against goals, so
  the design notes are evidence, not contamination).

## Verification

- Both reports recorded here; each REAL finding becomes a fix item
  or an accepted-risk note with VD visibility.

## Status

- 2026-07-09: plan authored; both readers dispatched.
- 2026-07-09: both reports in; adjudicated. Review COMPLETE.

READER A (architecture) verdicts:
- grow-from-root-warm-start FITS: rides installForests/
  sampleTreesFromPrior, never enters the MH kernel (so the
  proposal-equals-prior identity is untouched by construction).
  Net-new: the per-cut scan primitive (benchmark-only today), which
  must be OCCUPANCY-AWARE - the empty-leaf veto lives in
  logLikelihoodForBranch's finite penalty and does not transfer to
  a scan-based builder (noted on the TODO entry).
- multi-forest-models STRAINS (deepest): the Forest split moved
  state, not combining - BCF is a hardcoded bcf_ special case
  (sweep branches, drawGlue, combinedFits, formForestResponse all
  two-forest-literal; no ForestCombiner abstraction), sigma_ is a
  chain-level scalar, and 36 forests_[0] reporting hardcodings are
  per-forest quantities for K-way models. Heteroscedastic
  additionally forces a FOURTH leaf kind (multiplicative-positive,
  non-integrable) bundled with a non-conjugate MoveStrategy that
  does not exist. Filed: forest-combiner (TODO).
- negative-binomial + ordinal FIT engine-side (ResponseModel is the
  right seam; GroupedResponse's weighted update composes with PG/
  probit latents). The strain is the FROZEN dbarts_results struct:
  it cannot gain fields (consumers stack-allocate against the
  compiled layout), so NB's r trace / ordinal cutpoints need a
  decided growth path BEFORE 1.0-0 (VD decision window).
- data-ownership FITS (borrow isolated in ColumnStore; C ABI
  unaffected) but there is NO column-subset view - a hard
  prerequisite for BCF moderators and per-forest predictors (noted
  on the TODO entry).
- blocked-jacobi FITS behind the move surface; gpu-bart BLOCKED
  (representational; correctly sequenced behind the scan
  primitive).
Cross-cutting: GroupedResponse decoration is the right pattern on
the wrong side of the Chain/combiner boundary - grouped x BCF (and
grouped x any multi-forest) is impossible until combining becomes a
ResponseModel-side object. Ranked debts: (1) ForestCombiner,
(2) dbarts_results growth path [pre-release decision],
(3) additive/optional-block state format [pre-release; flat-format-
v2 is the existing item], (4) column-subset views, (5) per-
observation residual variance convention, (6) occupancy-aware scan.

READER B (numerical) verdicts, measured (scripts in session tmp):
- REAL, high: grouped tau slice sampler step-out is UNCAPPED in
  BOTH directions (model.hpp:2299-2301; only shrinkage has the
  1000 cap; verified by orchestrator). Step-outs ~ tau/width:
  >1e5 iterations at tau > ~2.5e5, indefinite hang at tau > ~1e8.
  Reachable legally: empty groups draw effects ~N(0, tau^2) with
  no mean-reversion (mostly-empty groupings let tau random-walk
  up); heavy-tail excursions at small J. Fix: Neal's m cap both
  directions (~1e4) - never engages in normal runs, bit-identical.
  Filed: tau-slice-stepout-cap.
- REAL, high (headline): chi-k Gibbs runaway - chi()'s scale=Inf
  DEFAULT (R/model.R:451, verified) is an improper prior; when
  leaves are prior-dominated the fixed-point factor
  (1 + df/numLeaves) > 1 always. Legal chi(100, Inf) runs away
  deterministically (k -> 1e25; leaf sd -> 0; forest fits NOTHING;
  sigma -> marginal sd; NO error raised); few-tree df=1.5 runs away
  stochastically (k to 1e21 observed). chi(100, scale=2) stabilizes
  healthily at k~21. The default 200-tree chi(1.5) survives only
  because data usually overcome a 1.0037 per-sweep factor. Remedy
  is VD's call (sentinel cap / finite default scale / warn) -
  filed chi-k-runaway; interacts with chi-default-research.
- ACCEPTED RISK: constant-leaf centered-SS cancellation (garbage
  at shared offset ~1e6-1e8, but the response is pre-standardized
  to [-0.5, 0.5], 0 decision flips through offset 1e6; the BCF
  resid/m scaling cancels in the suffstat products); totalFits
  drift (sqrt-growth, ~1e-13 at 1e5 sweeps, invisible vs sigma;
  optional re-anchor not needed); BCF 1e-9 floor (m cancels to
  full precision in w*z products; sub-floor leaves draw from the
  prior, immaterial); extreme uniform response scales (1e15/1e-15
  round-trip clean); setState gaussian n-mismatch (accepted but
  benign: repartition + re-accumulation yields healthy fits, only
  donor sigma/scale linger and burn-in absorbs them; the guard IS
  load-bearing for latent/function-leaf/grouped states, which are
  n-checked).
- DOC NOTES: GP predict at training rows re-krigs WITHOUT the
  self-term nugget, landing nugget*|K^-1 f| ~ 2-3e-3 below the
  recorded (jittered) fit - predict interpolates the jitter-free
  mean; affects predict/test only, never MCMC; huge-offset/
  tiny-range responses quantize in DOUBLE PRECISION before the
  engine sees them (y in [1e15, 1e15+1e-3] collapses to 1 distinct
  value; engine copes, models a degenerate response silently) -
  R-boundary warning filed as input-precision-validation.
- FORWARD GUARDRAIL (ties A to B): seeds 1 and 5 are safe ONLY
  because responses reach the constant leaf pre-standardized to
  O(1); the data-ownership container or any family bypassing
  rescale() re-exposes the cancellation and the floor's scale
  assumption (noted on the data-ownership TODO entry).
Overall: numerically trustworthy at the documented envelope; the
two uncapped/improper objects (slice step-out, improper chi) are
the trust boundary, both cheap to guard.
