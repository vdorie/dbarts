# robust-errors

agent: Opus (C1-C2), Sonnet candidate for C3's docs half
rng: neutral for every EXISTING anchor (the gaussian stream gains no
  draws; the t path is new); the new equivalence scenario is recorded
  with a neutrality trail, grouped_aft precedent
budget: C1 ~350 lines; C2 ~200; C3 ~250

## Goal

Continuous responses gain Student-t residuals via the scale-mixture
augmentation on the workingWeights hook. resid.dist = student()
estimates nu on a capped grid; student(df = x) fixes it; default
gaussian(). Design: docs/design/robust-errors.md (all sections
resolved; section 8 records the decision and the precedent survey).

## Context

The design note carries the math and seams; implement against it, not
this summary. Key anchors: TResponse as an AFTResponse-pattern
decorator ([[model.hpp:2405@8e4e7140]]) over a contained GaussianResponse via
setWeights(w_i * lambda_i) ([[model.hpp:1968@8e4e7140]]); the lambda conditional
Gamma((nu+1)/2, (nu + w_i r_i^2/sigma^2)/2) drawn with
ext_rng_simulateGamma; the nu grid draw from (sum log lambda,
sum lambda) with precomputed per-point constants (DartPrior pattern,
[[model.hpp:1655-1707@8e4e7140]]); the weighted sigma path already exact
([[model.hpp:1758-1764@8e4e7140]]).

## Commits

C1 (Opus, engine): TResponse with both modes. Fixed mode: nu constant.
  Grid mode: nu on {3, 4, 5, 6, 8, 10, 12, 15, 20}, prior gamma(2,
  0.1) evaluated at the points and normalized (the survey's
  convention; the cap per the note's section 4). refreshLatents draws
  lambda then nu (grid mode), composes c_i, delegates.
  computeLogLikelihood reports dt. latents() exposes lambda (and the
  nu grid index in grid mode) for the state block; cold-init
  lambda = 1, nu = the grid median (probit-latent convention).
  tests/cpp: lambda conditional moments on fixed residuals (the
  PG-moments precedent); grid full-conditional against a
  hand-computed log posterior on a 3-point grid.
C2 (Opus, facade + bridge + state): response selection wiring, lambda
  and nu in the by-name state block (new names; old states load
  unchanged), bridge argument plumbing. NO dbarts.h change: the C API
  does not expose t in v1 (stan4bart has no need); record the door.
C3 (R surface + gates): resid.dist argument on dbarts()/bart()/
  bart2() taking gaussian() or student(df = NULL); parsePriors-style
  resolution; student() refused informatively for non-gaussian
  families. Rd + NEWS + the two documented caveats (flexible mean
  confounds tail inference - posterior-check estimated nu; sigma is
  the conditional scale). Exact gate: fixed-nu single-tree quadrature
  reference per the note's section 6 (marginal t, never augments).
  New equivalence scenario (contaminated normal, student()) ->
  re-record the equivalence baseline WITH a neutrality trail (the
  existing 22 scenarios must reproduce ac6ec2c bitwise inside the
  compare; MANIFEST entry states it). air format + lintr; pkgdown
  reference entry if a new Rd topic lands.

## Constraints

- Existing anchors bitwise until C3's re-record, whose trail proves
  the 22 untouched; bcf + multinomial fixtures untouched throughout.
- Serialize after the data-review-remediation arc; one implementer at
  a time.
- The gaussian path gains no draw and no per-iteration branch (the
  decorator composes outside; workingWeightsVaryPerSweep gates the
  per-sweep work exactly as logistic does).
- setResponse/setData cold-reinit lambda to 1 (design note section 2).
- Out of scope: robit (recorded follow-up, needs its own design
  pass); asymmetric/quantile errors; a bart2 resid.df scalar alias
  (add later only if asked).
- Stop conditions per docs/plans/README.md; any existing-anchor
  divergence = stop, report, never re-record without the trail.

## Verification

Per commit: install --preclean, tests/cpp from clean, full tinytest,
three equivalence compares (bitwise; after C3 against the NEW
baseline, trail verified). C3 additionally: the exact gate passes at
tolerance like logistic-reference; air format --check . exits 0;
pkgdown check on any new Rd topic.

## Landings

C1 045a058 (2026-07-18). TResponse + ResidualDfPrior in model.hpp;
one-line Rf_dt decl in external/stats.h. Both modes; grid {3..20}
under normalized gamma(2,0.1); the grid test verifies against a
reference keeping the term the code drops. Deviations accepted: nu
via dedicated accessors (latents() stays length n); zero-weight
lambdas drawn for stream regularity, excluded from nu statistics.
All six gates green, anchors bitwise (nothing instantiated it).

C2 9b58662 (2026-07-18). SamplerOptions.residualDf (NaN/Inf =
gaussian, > 0 fixed, <= 0 estimated), gaussian construction selects
TResponse when finite; lambda rides the EXISTING latents block, nu is
the one new by-name scalar ("resid.df"), written only for t; old
states load into gaussian samplers unchanged, cross-family restores
refused; bridge reads the model's resid.df attribute and refuses
non-gaussian families and negative values host-side; GroupedResponse
forwards the nu channel. State round-trips fixed + estimated, RNG-
insulated; continuation statistical per the file's own convention
(lambda/nu themselves bitwise). No dbarts.h diff. Gates green.

C3 31dc05a (2026-07-18). resid.dist = gaussian()/student(df) on
dbarts/bart/bart2, resolved in the parsePriors environment (no
stats::gaussian shadowing), df validated at construction, family
refusal R-side mirroring the bridge backstop; Rd + NEWS with both
caveats (confounded tail inference; sigma conditional, marginal
sigma^2 nu/(nu-2)). benchmarks/R/t-exact.R: fixed-nu single-tree
quadrature reference (never augments) - max arm gap 0.0010 vs
tolerance 0.012. New student equivalence scenario (contaminated
normal, estimated nu); baseline re-recorded as equivalence-31dc05a
(23 scenarios) with the neutrality trail: all 22 pre-existing
scenarios reproduce ac6ec2c bitwise. equivalence.yaml re-pinned.
Suite 3088/0. One implementer stall (backgrounded gates), resumed by
message; reviewer re-ran install, trail, t-exact, and suite
independently. ARC CLOSED; implementer file overage (483 vs ~375)
was the two mandated new files.
