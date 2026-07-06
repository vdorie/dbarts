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
