# weighted-binary

agent: opus
rng: posterior-changing (new likelihood support)
budget: decision memo exists below; decision-gated (2026-08-06: the
  post-1.0 framing is retired - no item is deferred to post-release,
  VD 2026-07-17; "Neither lands in 1.0-0" below is history, a
  statement of 1.0-0's content, not a standing gate)

## Goal

The two weighted-binary follow-ups stay decision-gated with their
analysis preserved: integer-weight probit, and arbitrary real logistic
weights. Neither lands in 1.0-0.

## Context (analysis carried from the pre-rewrite TODO, 2026-07-06)

- Integer-weight probit: Phi(f)^w is exactly tractable as w replicated
  latent draws - the tree side stays O(n) because the latent mean is a
  sufficient statistic (zbar_i ~ N(f_i, 1/w_i), the existing
  weighted-gaussian machinery); only the latent-update phase scales as
  Theta(sum w_i), linear in the UNAGGREGATED count, which silently
  undoes the aggregation that motivates count weights. The shipped
  weighted logistic has the same structure (fbd2168 draws PG(w, psi)
  as a sum of w PG(1, psi) variates), so replication would not be a
  regression - both families exact but linear in total counts. If
  large counts matter, the standard fix is a threshold hybrid: exact
  summation for small w, a large-w approximation above it (CLT normal
  approximation for the truncated-normal sum; Windle et al. large-b
  Polya-Gamma samplers) - O(1) per observation but perturbs the
  stationary distribution, so it must be a deliberate, documented
  switch.
- Arbitrary real logistic weights: needs a general real-shape PG
  sampler; the survey/IPW use case wants it but is a pseudo-likelihood
  trap (docs/design/weighted-logistic.md - IPW point estimates with
  badly overconfident uncertainty; AIPW/TMLE keep weights out of the
  BART fit deliberately).

## Constraints

- Post 1.0; each needs its own design note and exact-posterior gate on
  landing (the weighted-logistic.md pattern).
- Out of scope now: everything but keeping this analysis findable.

## Steps

1. When picked up: design note per docs/design/weighted-logistic.md's
   structure; implementation plan follows approval.

## Verification

- On landing: exact-posterior gate (single-tree enumeration +
  quadrature) for the new likelihood, per the established pattern.
