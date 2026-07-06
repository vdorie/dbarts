# interaction-constraints

agent: opus
rng: posterior-changing (prior change)
budget: memo on demand; consumer-gated

## Goal

Structured split-variable control on the SplitSelector seam when a
consumer needs it: grouped sparsity (DART over variable groups) and
interaction limits (variables that may not co-occur on a path, or
max-interaction-depth constraints a la BART-based causal-inference and
tabular-ML practice).

## Context

- The seam was designed for exactly this: "Grouped/structured DART,
  interaction constraints: SplitSelector"
  (docs/design/core-generalization.md, extensions table); SplitSelector
  log-probabilities already enter acceptance ratios and have a
  per-iteration hyper-update hook.
- Path-dependent constraints (no-co-occurrence) touch move validity,
  not just selection weights - the change/swap good-rule flows would
  need the constraint in their satisfiability walks; the memo must
  scope that honestly.

## Constraints

- Consumer-gated: no code without a named workload (candidates:
  causal setups wanting treatment-effect additivity, grouped
  genomics-style predictors).
- Out of scope: monotone constraints (monotone-bart), soft splits.

## Steps

1. On demand: memo covering the constraint vocabulary, which seam each
   piece lands on (selection prior vs move validity), and prior
   math for the acceptance ratios; exact-posterior gate design for a
   constrained toy problem.
2. Fresh implementation plan on approval.

## Verification

- Memo reviewed; TODO updated with the outcome.
