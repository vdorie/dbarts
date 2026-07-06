# ordinal-outcomes

agent: opus
rng: posterior-changing (new family)
budget: design note + ~400 lines

## Goal

Ordered categorical responses fit via cumulative probit: the existing
truncated-normal latent machinery plus a cutpoint block, surfaced as
family = "ordinal" (or ordered-factor response auto-dispatch - the
note decides).

## Context

- Probit latents already exist (ProbitResponse,
  src/bartcore/model.hpp); ordinal generalizes the truncation from
  {(-inf,0), (0,inf)} to per-category intervals (gamma_{k-1},
  gamma_k].
- The cutpoint update is the known hard part: Albert-Chib's Gibbs
  cutpoints mix badly; the note should evaluate the standard fixes
  (Cowles' MH block update, or fixing K-2 cutpoints via a latent-scale
  identification and sampling a scale). Identification also interacts
  with the response scaling and node.scale calibration.
- Ingestion: ordered factors currently fit as ordinal predictors;
  as a response they error - dbartsData response typing gains a path.

## Constraints

- Exact-posterior gate (single tree, few categories: enumeration +
  quadrature over latents and cutpoints on a grid) before shipping.
- rbart_vi/grouped decorator composition is desirable (ordinal
  mixed models) but not v1; refuse cleanly.
- Out of scope: unordered multinomial (multi-forest-models);
  ordinal-scale loss functions in xbart (follow-up).

## Steps

1. Design note (docs/design/ordinal.md): identification scheme,
   cutpoint sampler choice with a small mixing study, surface and
   prediction semantics (category probabilities vs latent scale -
   probabilityFromLatents generalizes).
2. OrdinalResponse + cutpoint block; ingestion typing; R surface,
   predict/fitted transforms, Rd.
3. Gates: exact-posterior, component tests (truncation intervals,
   cutpoint update), recovery, mutation smoke (setResponse re-draw
   semantics defined).

## Verification

- Exact-posterior gate to MC error; component tests; full tinytest;
  bench on existing families unchanged.
