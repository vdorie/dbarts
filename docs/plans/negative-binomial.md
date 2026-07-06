# negative-binomial

agent: opus
rng: posterior-changing (new family)
budget: design note + ~500 lines

## Goal

Count responses fit natively: family = "nbinom" running the
Polya-Gamma negative-binomial augmentation on the existing weighted
conjugate path, with offsets acting as log-exposure.

## Context

- The machinery mostly ships: NB augmentation (Polson-Scott-Windle /
  Zhou et al.) is omega_i ~ PG(y_i + r, psi_i) with working response
  and precision entering exactly like the logistic port did
  (docs/design/core-generalization.md phase 3; workingWeights hook,
  src/bartcore/model.hpp LogisticResponse is the template).
- The r (dispersion) update is the design decision: PG shape y + r is
  non-integer whenever r is continuous, and the shipped PG sampler is
  exact for integer shape only (it sums PG(1, psi) draws - the same
  real-shape gap recorded in weighted-binary). Options for the note:
  (a) fixed user-supplied r; (b) integer-sampled r (CRT/gamma per Zhou
  & Carin, rounded support); (c) implement a real-shape PG sampler
  (unblocks weighted-binary's real-weights half too - if the note
  picks this, split it out as its own shared item).
- Exposure/offset semantics: psi_i = f(x_i) + offset_i with offset as
  log exposure; document against the gaussian offset behavior.

## Constraints

- Exact-posterior gate mandatory (single tree, one predictor, small
  counts: enumeration + quadrature over the leaf, fixed r), per the
  weighted-logistic.md pattern; plus PG moment component tests at the
  shapes used.
- Prior calibration on the latent scale needs its own paragraph in
  the note (node.scale under NB's latent variance, following the
  logistic pi*sqrt(3) precedent).
- The embedded-Gibbs mutation surface must define setResponse for
  counts (latent re-draw semantics) - the probit port's pattern.
- Out of scope: zero-inflation (multi-forest-models' hurdle), Poisson
  (the r -> infinity limit; revisit after NB lands).

## Steps

1. Design note (docs/design/negative-binomial.md): r update choice,
   latent scale calibration, offset semantics, surface (family value,
   r argument or prior).
2. NBResponse riding workingWeights; family plumbing through control/
   model per the family-on-model outcome; R surface + Rd.
3. Gates: exact-posterior, PG shape moment tests, recovery test
   (dispersion and log-mean recovered on synthetic), mutation smoke.

## Verification

- Exact-posterior gate to MC error; component tests; full tinytest.
- bench: gaussian/binary paths unchanged; NB cost recorded (PG draw
  per observation per iteration dominates, as logistic).
