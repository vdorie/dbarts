# survival-models

agent: opus
rng: posterior-changing (new families)
budget: design note + ~500 lines for the first model

## Goal

Time-to-event responses fit natively, in two steps the design already
provisioned: AFT log-normal (censoring via truncated-normal latents)
and discrete-time hazard (person-period expansion + binary family).
The design note picks which lands first; AFT is the plan's default
recommendation.

## Context

- Provision: "Survival (AFT, discrete-time hazard), ordinal, quantile:
  ResponseModel latents; person-period expansion at ingest"
  (docs/design/core-generalization.md, extensions table).
- AFT log-normal: log T_i = f(x_i) + sigma * eps_i; right-censored
  observations draw latent log-times from the normal truncated below
  at the censoring time - the probit truncation machinery generalized
  (src/bartcore/model.hpp ProbitResponse::refreshLatents is the
  template). Uncensored data reduce exactly to the gaussian path.
- Named consumer: riAFTBART wraps AFT-with-random-intercepts around
  dbarts today; the grouped decorator composes, replacing its outer
  loop like rbart_vi's.
- Discrete-time hazard: ingestion sugar (person-period expansion with
  interval index as a predictor) + probit/logistic; mostly a data and
  surface item once AFT's response plumbing exists.
- Surface sketch: response as Surv-like (time, status) - accept
  survival::Surv without importing survival (inherits() check),
  family = "aft" or auto-dispatch on the response class; the note
  decides.

## Constraints

- Exact-posterior gate: single-tree AFT with a censored observation
  admits enumeration + quadrature (truncated-normal marginals);
  uncensored-data fits must match gaussian fits exactly (strong free
  gate).
- Prediction semantics defined in the note: survival curves need
  sigma draws + f draws; decide what predict() returns (median time,
  linear predictor, curve evaluations).
- setResponse under censoring (latent re-draw) defined, per the
  probit pattern.
- Out of scope: nonparametric baseline hazards (BART-package-style
  heteroscedastic survival), competing risks, left/interval censoring
  in v1 (right-censoring only; the latent code should not preclude
  them).

## Steps

1. Design note (docs/design/survival.md): which model first, surface,
   prediction semantics, riAFTBART composition check.
2. AFTResponse (censored-latent refresh; sigma update includes latent
   times); ingestion of (time, status); R surface + transforms + Rd.
3. Gates: uncensored == gaussian equality, exact-posterior censored
   gate, recovery under known censoring rates, mutation smoke,
   grouped-decorator composition test.
4. Discrete-time hazard as a follow-up section landed against the
   same note (person-period ingestion + existing binary families).

## Verification

- The equality and exact-posterior gates; component tests; full
  tinytest; bench on existing families unchanged.

## Status

- 2026-07-10: AFT log-normal LANDED as f0efc03 (squash of wt/survival-aft;
  design note docs/design/survival.md rides the commit). Engine: AFTResponse
  composes a contained GaussianResponse over a mutable log-time buffer;
  right-censored latents redraw each sweep from the lower-truncated normal;
  zero facade changes (status arrives via SamplerOptions.survivalStatus off
  the bartcore.survival control attribute, the groupIndices pattern); state
  rides the existing latents + fit.scale blocks, no format bump. Weights and
  setData refused; setResponse keeps the censoring structure and redraws.
- Surface (revised to a three-agent design panel's synthesis, VD-approved):
  Surv or two-column (time, status) response through the matrix interface;
  Surv auto-dispatches to aft only from family "auto" - an explicit
  conflicting family errors, never a silent override; the formula interface
  refuses aft/Surv with a matrix-interface hint (both guards) instead of
  failing hostilely downstream. survivalProbabilities is an S3 GENERIC
  returning DRAWS (draws x times x observations; chain margin under
  combineChains = FALSE) per the extract/fitted/ci.level tiers; the rbart
  method errors (grouped AFT unreachable until rbart_vi grows a family
  argument - filed on TODO as survival-grouped-surface). predict/extract
  stay on the LOG-TIME scale (riAFTBART-source-verified); the Rd documents
  the scale and the inert response/link aliases.
- Gates (worktree, private library; the landed tree is byte-identical to
  the gated tip): preclean install; tests/cpp all pass (bitwise
  uncensored == gaussian reduction, censored-latent truncated-normal
  moments, censored state round trip); tinytest 2661/0 (2610 + 51);
  equivalence 21/21 identical vs equivalence-de67cbb.rds; aft-exact.R
  max gap 0.0006 with a poisoned-truncation-bound proof (0.1153 fail);
  air + lintr clean; pkgdown::check_pkgdown clean. Reviewer re-ran
  install + tests/cpp + tinytest on the landed tree from the shared lib.
- Follow-ups: discrete-time hazard (person-period ingestion over the
  binary families, same note) and survival-grouped-surface (rbart_vi
  family = "aft" + the rbart survivalProbabilities method) remain open.

## Discrete-time hazard landing

Design f494156 (2026-07-19; the "Discrete-time hazard" section of
docs/design/survival.md). Independently critiqued: REJECT round caught
the false no-new-R-code thesis (the hazard token cannot survive the
family-keyed switches - node.scale, control@binary, fixedUnitScale,
weight policy - so the design became remap-before-switches with a
packaged-fit marker; $family keeps the binary token so every link
generic is correct unchanged), the Surv-guard collision, and the
ragged-training-design prediction path. VD then requested the
precedent survey (surv.bart verified source-level: internal expansion,
unique-times default grid, K quantile coarsening, probit default;
discSurv/pammtools conventions; link traditions), which moved the
time-grid fork to unique-observed-times default with breaks
coarsening. VD resolved all three forks 2026-07-19: that grid, the
two tokens (hazard = probit, hazard.logistic; hazard.probit shipped
as a self-documenting alias), probit default.

Implementation e58ab73 (2026-07-19). R-only: person-period expander
(resolveHazardGrid, expandDiscreteTimeHazard, appendHazardPeriodColumn;
breaks = NULL distinct times / integer quantile count / boundary
vector; max.rows refusal guard), shared survival parsing refactored so
the hazard sibling keeps raw time, the Surv guard admits the hazard
tokens, remap precedes every family-keyed switch, $periods is the one
new packaged field and carries dispatch, survivalProbabilities gains
the re-expanding keepTrees-requiring branch, offsets/weights replicate
per subject. Formula interface and in-fit test sets refused
informatively (matrix interface + predict-time expansion). Gates:
benchmarks/R/hazard-reduction.R - the fit reduces BITWISE to the
hand-expanded binary fit on both links (yhat/varcount/trees), the
family's primary correctness gate; suite 3230 -> 3267; new hazard
equivalence scenario, baseline re-recorded as equivalence-f494156.rds
(engine binary unchanged since 1d9388d) with the trail: all 25
pre-existing scenarios bitwise; equivalence.yaml re-pinned. Doors
recorded in the design section 6: cloglog, competing risks, left
truncation, grouped frailty, time-varying covariates in the
auto-expander, continuous-time models. ITEM CLOSED except
survival-grouped-surface (its own plan).
