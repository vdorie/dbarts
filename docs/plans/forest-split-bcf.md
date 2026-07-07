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

## Status (2026-07-07)

Step 1 landed: docs/design/bcf.md reviewed by VD; resolutions recorded
there (range-anchoring kept with the approximate sd(y) map onto
per-forest k hyperpriors; bcf's 200/50 tree convention for the BCF
entry point). Open questions 3 (muscale expansion vs a = 1) and 4
(tau seeing pihat) block only steps 3-5.

Step 2 landed (e6e5e7a): Forest<L> struct in chain.hpp holding the
per-forest members per the design note's member split; Chain holds a
size-1 std::vector<Forest>; sweep, mutation fan-out, and prediction
loop over it. ChainStateData and every serialized format unchanged
(getState/setState address forests_[0]; the forest-dimension bump is
step 4's). One file, +628/-511, all relocation and loop-wrapping.
Gates: component tests pass, tinytest 2470/0, equivalence exact
18/18 identical draws vs 235bebc, bench-sampler vs 235bebc clean
(no metric regressed; ratios 0.91-0.97). Deviations: Forest lives in
chain.hpp rather than a new header; the param-carrying mutation
methods (revalidateTrees and kin) address forests_[0] directly until
TreeParameters gains a forest dimension; per-forest option fields are
duplicated on Forest while SamplerOptions stays the unsplit bridge
struct - Sampler-level options_ was already stale after setModel
pre-split (kIsSampled), unchanged and noted for step 3.

Step 3 landed (d1ffb92): the BCF two-forest sampler behind an
internal-only surface. Chain gains a BCF constructor (constant leaf,
Gaussian response): per-forest backfitting against the residual net
of the other forest's scaled contribution (response r/m with weight
w m^2, the exact weighted reduction), the combined a mu + b_z tau
fit feeding latents/sigma, and the glue draws - both expansions per
the resolutions (b0/b1 Gaussian conditionals over the treatment
partition; the prognostic a with its half-Cauchy via the t_1
inverse-gamma scale mixture). setTreatment, per-forest fits and glue
accessors fan through Sampler/facade/bridge to internal .Call entry
points and R wrappers; no dbarts.h change, nothing exported.
Multi-forest state serialization refuses cleanly (step 4's).
Interims recorded: both forests read the full store (the moderator
subset waits on data-ownership's views); per-forest scales are
placeholder constants pending step 5's exact map; single-forest
queries address forest 0. Review fixed one defect: the scale-mixture
rate carried 1/scale^2 instead of scale^2 (invisible at the default
scale 1; wrong anywhere else). Gates: component tests incl. the new
BCF test, tinytest 2484/0 (12 new), equivalence exact 18/18 - the
single-forest sweep is pointer-identical when the BCF state is null.

Steps 4-5 remain: serialization formats, then the exact-posterior
gate and the bcf-constant calibration map (which also replaces the
placeholder scales).

Step 4 landed (4e6b206): ChainStateData splits into chain-level state
plus a per-forest vector (ForestStateData: tree/saved/param/mask
channels and k); the serialized object gains a per-chain "forests"
list (length 1 off BCF) and an optional "bcf" glue vector
(a, aVariance, b0, b1); stateFormatVersion 2 -> 3, version-2 states
refused by name; the step-3 multi-forest refusals replaced by working
round trips. Forest-count mismatches rejected in stateIsValid and
setState. Gates (implementer's, re-run by reviewer): component tests
incl. a 5-seed fuzzer pass over the new layout, tinytest 2492/0
(8 new), equivalence exact 18/18, air/lintr clean. Deviation: the
BCF fit round-trip test asserts to 1e-5 (setState is semantic
continuation); glue equality exact. Step 5 remains: the two-forest
exact-posterior gate + bcf calibration map. Reviewer note: run rchk
locally over the new serialization code before release (zero-findings
status last confirmed at c2d591a).

## Step 5 spec (2026-07-07)

Budget ~550 lines. Three deliverables plus tests:

1. Calibration map, replacing step 3's placeholder scales. At BCF
   creation compute s = sd of the range-scaled response (y mapped to
   [-0.5, 0.5]). Defaults per the bcf-parity resolution: mu forest
   nodeScale = s with aPriorScale = 2 (half-Cauchy median of |a| = 2
   puts the prognostic total at bcf's 2 sd(y); the map overrides the
   host model's k/node prior for mu - note that on the wrapper); tau
   forest nodeScale = s / 0.674 with bPriorVariance = 0.5 (b1 - b0 ~
   N(0, 1), half-normal median 0.674, so the effect scale's median is
   one sd(y) - bcf's use_tauscale correction). k stays 1. The R
   wrapper's scale arguments become sd(y)-unit knobs (sd.control = 2,
   sd.moderate = 1) converted internally. Record the derivation as a
   short Calibration subsection in docs/design/bcf.md.
2. Glue-fixing switches: BCFSpec gains updateA/updateB (default
   true), threaded through the bridge params and R wrapper; false
   skips the matching drawGlue block (and the aVariance refresh when
   a is fixed).
3. Exact-posterior gate, benchmarks/R/bcf-exact.R, following
   categorical-exact.R's conventions (quick mode, tolerance, exit
   status). One ordinal predictor over a handful of distinct values,
   fixed 0/1 z balanced within each x cell, Gaussian response, two
   single-tree forests (host n.trees = 1, n.trees.treatment = 1).
   Enumerate each forest's tree space independently (the ordinal
   analog, logistic-reference.R:63); the joint space is the product.
   Conditional on (a, b0, b1, sigma) the leaf parameters integrate in
   closed form (Gaussian block marginal); 1-D quadrature over sigma
   against the sampler's actual prior on the scaled response -
   reproduce dbarts's range scaling and chisq(df, quantile) sigma
   calibration exactly (the one tricky reproduction), and reproduce
   the deliverable-1 map for the leaf scales so the gate validates
   the map's implementation end to end. Fits from bartcoreForestFits
   are internal-scale; match on that scale at the distinct x cells,
   to MC error over a few seeds. Modes: (1) glue fixed a = 1, b0 = 0,
   b1 = 1, match E[mu] and E[tau]; (2a) a free (quadrature against
   Cauchy(0, aPriorScale)), b fixed, match E[a mu] and E[tau];
   (2b) b0/b1 free (N(0, bPriorVariance) grid), a fixed, match E[mu]
   and E[(b1 - b0) tau]. A fully-free mode is not required.
4. Tests: a component test for the fixed-glue path; tinytest
   additions including a bartCause-style driver (setTreatment swap
   plus a pihat column update through setPredictor between runs).

Single-forest paths must stay untouched: equivalence exact 18/18 is
the gate, plus component tests, full tinytest, the new gate in both
modes, and air format + lintr on touched R files.
