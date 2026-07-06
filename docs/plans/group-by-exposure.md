# group-by-exposure

agent: sonnet
rng: neutral (surface only)
budget: decision memo first

## Goal

A recorded decision on exposing grouped random effects beyond rbart_vi
(a group.by argument on dbarts()/bart2), currently internal by design.

## Context

- docs/design/grouped-random-effects.md: grouping is an attribute on
  the control object, deliberately internal; "no public group.by
  exposure beyond rbart_vi yet"; naming reserved.
- The in-core path already handles the built-in tau priors
  multi-chain; the R loop remains for custom priors and callbacks.

## Constraints

- Blocked on demand: rbart_vi covers the known use. Exposure is an API
  commitment (argument name, formula semantics for bart2) that cannot
  be walked back after release - but adding it later is additive, so
  there is no pre-release window pressure.
- Out of scope: any engine change.

## Steps

1. When demand appears: memo covering argument shape (bare column vs
   formula term), interaction with family (probit-only today via
   rbart_vi's construction - the in-core decorator composes wider,
   which changes the user-facing story), and what rbart_vi becomes.
2. Fresh implementation plan on approval.

## Verification

- Memo reviewed by VD; TODO updated with the outcome.
