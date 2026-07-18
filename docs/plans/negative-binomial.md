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

## Landings

Design 13bdfa4 (2026-07-18). docs/design/negative-binomial.md,
Opus-written then independently critiqued: REJECT with two blocking
findings, both fixed pre-implementation - (1) the drafted sweep order
(omega -> r -> rebuild working) was proven non-invariant (the CRT/grid
r-draw is collapsed over omega; PCG requires regenerating the
marginalized variable before next conditioning on it), corrected to
r-first; (2) the real-shape PG exactness story was internally
inconsistent (the truncated gamma-sum is approximate, the ecosystem's
own tools approximate non-integer shape, and no posterior-level gate
can resolve the ~1e-3 bias), restructured into the exactness fork.
VD resolved fork (A): integer dispersion only, fixed or estimated on
a capped grid by a closed-form discrete conditional (the
ResidualDfPrior pattern exactly); real dispersion behind a recorded
door carrying fork (B)'s spec; forward-compat accommodations binding
(neutral real-valued "dispersion" state slot, numeric surface
argument with informative integer-only refusal, shape-parameterized
PG helper and r-update seams). The prototype
(benchmarks/R/negbin-r-update-mixing.R) validated both r-update
kernels against a grid posterior and measured CRT-Gamma at 2.5-4x
MH's ESS; fork (A)'s grid conditional supersedes both.

C1 a3ce9fd (2026-07-18). NBResponse + NBDispersionPrior +
simulatePolyaGammaShape in model.hpp; logit-p parameterization,
kappa = (y - r)/2, working = kappa/omega - offset; r-first sweep
order locked by a bitstream-replay component test; grid conditional
from a per-response count-histogram lgamma kernel (one O(n) reduction
per sweep); carriesDispersion virtual trio with the
restoreDispersion-before-restoreLatents contract; PG helper bit-equal
to the shipped stream at shape 1 (tested). One review round:
setResponse initially redrew r through refreshLatents, contradicting
the note's kept-r embedded-Gibbs clause - fixed by extracting
drawOmega, with the corrected contract locked by a stream-replay
assertion (an r draw would consume a uniform and break it). Five
RNG-insulated component tests. All anchors bitwise.

C2 1d9388d (2026-07-18). ResponseFamily::nbinom; construction in
chain.hpp off SamplerOptions.dispersion (positive fixed, non-positive
estimated, NaN non-count); the by-name "dispersion" scalar state slot
(conditional write, absence-tolerant decode, stateIsValid refuses
non-finite/non-positive; restore ordered before latents per the
contract); bridge count shape via the bartcore.dispersion control
attribute (presence marks the shape, value is the spec); host
refusals: weights (exposure belongs in the offset), grouped NB
(recorded door), off-shape family names. No dbarts.h diff. State
round-trip component test, fixed and grid modes plus cross-family
refusals. All anchors bitwise.

C3 67d93d2 (2026-07-18). family = "nbinom" on dbarts/bart2,
explicit-only (counts carry no unambiguous class; auto on counts
stays gaussian); non-negative-integer validation in the numeric
branch; dispersion = NA_real_ estimates, positive integer fixes,
real fixed values refused with the forward-compat wording;
bartNegbin fit class (mean-count draws, latent psi, per-draw
dispersion; fitted/extract/predict/print/residuals; predict replays
bitwise); offset documented as log-exposure; rbart_vi/xbart refuse;
Rd + NEWS. Gates: benchmarks/R/negbin-exact.R omega-free quadrature
over leaf log-odds and the r grid, both arms (estimated: mean-count
gap 0.0012 vs 0.070, grid-posterior gap 0.0100 vs 0.025; fixed r = 5:
0.0042 vs 0.070); suite 3168 -> 3230; new nbinom equivalence scenario,
baseline re-recorded as equivalence-1d9388d.rds (named for the engine
HEAD; C3 is R-only) with the trail: all 24 pre-existing scenarios
bitwise vs 227f46a; equivalence.yaml re-pinned. Scenario records
psi/mean/r/varcount - omega locked transitively through the trees,
the ordinal-z precedent. ARC CLOSED; doors recorded: real dispersion
(TODO negbin-real-dispersion, the fork-B spec; gates weighted-binary
fractional weights), Poisson family, grouped NB.
