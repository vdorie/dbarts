# forest-split-bcf

agent: opus
rng: neutral for the split (bitwise-gated refactor);
     posterior-changing for BCF (new model)
budget: split ~600 lines; BCF ~800 lines + design note; separable PRs

## Goal

Forest becomes the composable unit the design promised (trees + leaf
model + split selector + working response + backfitting state), with
Sampler orchestrating one or more of them - and BCF (a prognostic
forest mu(x) plus a treatment forest tau(x), y = mu + z * tau + eps)
lands as the first two-forest sampler. bartCause consumes it.

## Context

- The deferral is recorded: "Forest is not yet split out of Sampler
  (do this when the facade lands in phase 2)" -
  docs/design/core-generalization.md:265-266; it never happened.
- Multi-forest was the designed provision for BCF/multinomial/
  heteroscedastic/hurdle (core-generalization.md:138-144).
- Chain currently owns trees/fits/response state directly
  (src/bartcore/chain.hpp); the split moves the per-forest members
  into a Forest struct Chain holds one-or-more of.
- BCF reference: Hahn, Murray, Carvalho (2020). Design decisions the
  note must fix: propensity score as a prognostic-forest covariate
  (bartCause already fits propensities); separate tree-count/prior
  defaults per forest (BCF convention: smaller treatment forest,
  tighter prior); moderator subset for tau(x); how z enters (binary
  treatment first, continuous later).

## Constraints

- The split alone must be draw-neutral: single-forest samplers produce
  bitwise-identical results before and after (equivalence exact mode
  is the gate). Land it as its own PR.
- BCF surface: internal first (bartcore helpers, like the data handle),
  bartCause drives from R; public dbarts-level exposure is a follow-up
  decision. The embedded-Gibbs mutation surface must work per forest
  (bartCause swaps response-side quantities).
- Exact-posterior gate for BCF: two single-tree forests on a
  one-predictor problem admit the same enumeration + quadrature
  treatment as the existing gates; build it.
- Out of scope: multinomial/heteroscedastic/hurdle
  (multi-forest-models); continuous treatment.

## Steps

1. Design note (docs/design/bcf.md): the decisions above + the Forest
   member split, reviewed by VD before code.
2. Refactor: Forest struct; Chain holds std::vector<Forest> (size 1
   everywhere today); mutation fan-out and state serialization become
   per-forest loops. Bitwise gate.
3. BCF ResponseModel combining forests; per-forest ModelParameters;
   creation path taking two model specs + treatment vector.
4. State/flat formats: per-forest tree channels (coordinate with
   flat-format-v2 and state-format-policy version bumps).
5. bcf exact-posterior gate + component tests + a bartCause smoke
   driver in inst/tinytest.

## Verification

- Step 2: equivalence compare reports exact; full tinytest unchanged.
- Steps 3-5: the new exact-posterior gate to MC error; component
  tests; bench-sampler no regression on single-forest paths.
