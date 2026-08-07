# grouped-equivalence

Status: CLOSED, LANDED (2026-08-07 correction) - the grouped scenario
shipped in 0e9ccca (2026-07-09, gate-hardening battery: group.by,
gamma tau prior, non-unit weights; equivalence.R result$grouped) and
grouped_aft followed in ac6ec2c (result$grouped_aft); both sit in the
current baseline (MANIFEST, equivalence-7903855.rds, 27 scenarios).
This doc predates the landing and was never updated; kept as history.

agent: sonnet
rng: neutral (harness only)
budget: ~80 lines in equivalence.R

## Goal

Grouped rbart_vi fits are covered by the equivalence gate if and when a
grouped refactor is planned; until then this item records why they are
excluded.

## Context

- Exclusion reason (TODO, 2026-07-06): grouped fits' multi-chain result
  shape does not fit the harness's single-fit summary.
- The in-core grouped path landed with its own statistical-equivalence
  gates at the time (docs/design/grouped-random-effects.md), so the gap
  is regression coverage, not landing validation.

## Constraints

- Conditional: only worth doing ahead of a change that touches the
  grouped Gibbs blocks (in-core decorator, tau updates). Otherwise skip.
- Out of scope: reshaping the harness's summary format for all
  scenarios.

## Steps

1. Add a grouped scenario that summarizes per-chain then pools
   (mean/sd of ranef, tau, sigma, train fit), fitting the existing
   summary-column comparison.
2. Record a new baseline including it; update the MANIFEST
   (see equivalence-ci).

## Verification

- equivalence.R record + compare round-trips the new scenario;
  z-mode passes across two independent record runs (seed check).
