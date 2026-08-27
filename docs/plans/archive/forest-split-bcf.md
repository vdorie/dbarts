# forest-split-bcf

agent: opus
rng: neutral for the split (bitwise-gated refactor);
     posterior-changing for BCF (new model)
budget: split ~600 lines; BCF ~800 lines + design note; separable PRs

## Status

LANDED, in two phases. Phase 1 - the Forest split plus the two-forest BCF
sampler (steps 1-5) - closed at c9fd2fe on 2026-07-07. Phase 2 - a
user-facing `moderators` restriction on the tau forest, sequenced after the
later data-ownership-4 refactor (C1-C4) - closed at f6804f1.

## Summary

This arc gave dbarts a composable Forest unit (trees + leaf model + split
selector + working response + backfitting state) with Sampler orchestrating
one or more of them, and landed BCF - a prognostic forest mu(x) plus a
treatment forest tau(x), y = mu + z * tau + eps - as the first two-forest
sampler, consumed by bartCause. Phase 1 did the split, the BCF response
model, the state/flat-format work, and the exact-posterior gate and
calibration map. Phase 2, picked up after data-ownership-4 built the
underlying per-forest column mask, added the `moderators` argument that lets
callers restrict tau to a covariate subset, plus its own exact-posterior
validation for the restricted case. Both phases gate on bitwise equivalence
for every single-forest and default-path run, and on exact-posterior
quadrature wherever a run actually samples a different posterior (BCF itself
in Phase 1, a restricted tau in Phase 2).

## Goal

Forest becomes the composable unit the design promised (trees + leaf
model + split selector + working response + backfitting state), with
Sampler orchestrating one or more of them - and BCF (a prognostic
forest mu(x) plus a treatment forest tau(x), y = mu + z * tau + eps)
lands as the first two-forest sampler. bartCause consumes it.

## Context

- The deferral is recorded: "Forest is not yet split out of Sampler
  (do this when the facade lands in phase 2)" -
  docs/design/core-generalization.md; it never happened.
- Multi-forest was the designed provision for BCF/multinomial/
  heteroscedastic/hurdle (core-generalization.md, the multi-forest provision).
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

## Phase 1 - Forest split and the BCF sampler

### Steps

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

### Verification

- Step 2: equivalence compare reports exact; full tinytest unchanged.
- Steps 3-5: the new exact-posterior gate to MC error; component
  tests; bench-sampler no regression on single-forest paths.

### Landing

(Recorded 2026-07-07. The step-3 full-store interim noted below is
discharged once Phase 2 landed through f6804f1 - see Phase 2's Landing.)

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
exact-posterior gate + bcf calibration map. Reviewer note: rchk was
re-run locally 2026-07-07 over the step-4 serialization code (tarball
at db96254): zero dbarts findings, same tool-noise-only output as the
c2d591a baseline run.

Step 5 spec (2026-07-07). Budget ~550 lines. Three deliverables plus
tests:

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

Step 5 landed (c9fd2fe), closing the item. The calibration map
replaces the placeholder scales: at creation s = sd of the
range-scaled response anchors mu's leaf scale at s (aPriorScale =
sd.control = 2 puts the prognostic total's half-Cauchy median at
2 sd(y)) and tau's at sd.moderate s / 0.674 (so the effect
(b1 - b0) tau, with b0/b1 ~ N(0, 1/2), sits at a median of one
sd(y)); k is fixed at 1 and the derivation is bcf.md's Calibration
subsection. BCFSpec gains updateA/updateB glue switches
(creation-time config, re-supplied on restore - no state-format
change needed). benchmarks/R/bcf-exact.R is the two-forest
exact-posterior gate: ordinal enumeration per forest over a 3-cell
design, closed-form Gaussian block marginals over the leaf
parameters via a whitened eigendecomposition vectorized across the
sigma grid, the chisq sigma calibration reproduced exactly, and
three modes (glue fixed; a free against its Cauchy prior; b0/b1
free on a Gaussian grid). One review fix round removed the dead
nodeScale/k knobs the map had orphaned (BCFForestSpec, the bridge
params, the R wrapper's treatment.k). Deviation: mode 2b's E[mu] is
a slow-mixing target (mu couples to the weakly identified
(b0 + b1) direction of the bscale ridge), so it carries a looser
tolerance with more seeds and heavier thinning; the identified
E[(b1 - b0) tau] stays tight, and the multi-seed mean matches the
closed form, confirming mixing rather than bias. Gates (implementer
ran all, reviewer re-ran all but the full-mode gate): component
tests incl. the new fixed-glue test, tinytest 2497/0 (5 new),
equivalence exact 18/18 identical draws vs 235bebc, bcf-exact.R
quick (max gap 0.0009) and full (max gap 0.017 vs 0.035 on the
2b E[mu], all others <= 0.0003), air/lintr clean.

## Phase 2 - moderator-restricted tau forests (post data-ownership-4)

Phase 1 shipped BCF with both forests reading the full covariate store,
noting the moderator subset as an interim waiting on the data-ownership
work's views. This second phase, sequenced once data-ownership-4 had landed
that mechanism, closes the interim: a `moderators` argument restricts the
tau forest to a column subset, with its own exact-posterior validation for
the restricted posterior (which a run using it samples, and which therefore
cannot ride the bitwise equivalence gate).

One implementer session owns the whole phase (the R surface, the exact
gate, and the landing note stay coherent with one owner; the engine
mechanism already exists from data-ownership-4). Budget ~350-550 lines: C1
~60 lines of R; C2 ~200-300 (the restricted exact-posterior gate); C3 ~80
tinytest lines (+ ~120 C/C++/R if the forest-addressed query in Q-C is
taken); C4 small. The steps are separable PRs; C1+C3 can land before C2.

RNG posture: the default path stays NEUTRAL - moderators absent reproduces
the shipped BCF draw, byte-for-byte. A run that PASSES moderators samples a
different, restricted posterior BY DESIGN and cannot be equivalence-gated;
it is gated by the exact-posterior protocol instead (step C2). Gate for the
phase as a whole: equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds;
full tinytest 2771 + the C3 additions, no regen of existing snapshots.

Bench posture: no bench-sampler run is needed for the R-argument step. The
restriction reuses the same nullptr-guarded per-forest mask
data-ownership-4 already bench-confirmed neutral (bench-sampler compare vs
bench-sampler-32fc7c8.csv, 0.96-1.0x; baseline re-recorded as
bench-sampler-4008675.csv); C1-C3 add no hot-path code, only R index
resolution and a benchmark/test consumer. Skip the compare.

### Context (what data-ownership-4 discharged, what remains)

- The engine mechanism is BUILT and proven. data-ownership-4 step 2
  (4e1fb5b) gave BCFForestSpec columns/numColumns (the borrowed
  1-based-resolved column list), the BCF ctor installs it on the tau
  forest only (mu unrestricted), and bartcore_createBCF grew a trailing
  moderatorsExpr argument parsed at src/R_interface_bartcore.cpp:1339-
  1351 (range-checked, consumed at construction, -> spec.tau.columns).
  A tests/cpp test proves tau containment with mu reading the full
  store. The per-forest availability mask (step 1, f7763ba) is
  nullptr-guarded, so absent = the default draw byte-for-byte.
- What does NOT exist yet, and is assigned here by data-ownership-4 Q2:
  (a) any user-facing `moderators` surface - R/bartcore.R:485 passes a
  literal NULL placeholder to the moderatorsExpr slot; (b) the
  two-forest exact-posterior validation at a RESTRICTED moderator set.
  bcf-exact.R validates only the full-store two-forest posterior.
- The BCF user surface is INTERNAL. bartcoreBCFSampler
  (R/bartcore.R:453) is the sole entry; it is unexported (not in
  NAMESPACE) and undocumented (no man/ topic). Consumers reach it via
  dbarts::: - inst/tinytest/test-bcf.R and bartCause. Public
  dbarts-level BCF exposure stays a deferred decision
  (public-surface.md section 5), so this phase lands the
  argument on the internal function, not on an exported wrapper.
- Forest-addressed queries are partial. bartcore_getForestFits
  (R_interface_bartcore.cpp:1653, R helper bartcoreForestFits) reads
  fits for any forest index via forestTotalFits; but getTrees/varcount
  address forest 0 only (getTrees at :3907 loops sampler.numTrees()
  with no forest axis). So the tau forest's SPLIT VARIABLES are not
  yet R-observable - relevant to how C3 checks containment (Q-C).
- Post-step-5 work not in this plan's earlier text but now in the tree
  (context for the restricted runs, no interaction with the
  restriction): the glue draw carries a ridge-interweaving rescale
  move (docs/plans/archive/bcf-ridge-interweaving.md) under update.a = TRUE,
  and bcf-exact-weak.R adds a prior-dominated a-glue gate. The
  restricted gate inherits both unchanged.

### Steps

C1. The `moderators` argument on bartcoreBCFSampler (R only). Add
    `moderators = NULL` to bartcoreBCFSampler (R/bartcore.R:453).
    Default NULL = all columns = today's byte-for-byte draw (the
    literal NULL placeholder at :485 is replaced by the resolved
    value). Resolve against the sampler's built columns, which are the
    already-expanded design matrix sampler$data@x (dbarts()/bartCause
    ran makeModelMatrixFromDataFrame upstream; factor expansion is
    NOT re-done here - see Q-A): accept a character vector matched to
    colnames(sampler$data@x) and/or an integer/numeric vector of
    1-based column indices, coerce to a sorted unique integer vector,
    and pass it (as.integer) into the moderatorsExpr slot the bridge
    already validates. Validation is R-side (safe over fast): error on
    an unknown name, an index < 1 or > ncol, an empty vector, or names
    given when colnames(sampler$data@x) is NULL. Update the
    bartcoreBCFSampler roxygen/comment header to document the
    argument. rng: NEUTRAL (default path unchanged; the bridge and
    engine are untouched). Gate: equivalence 22/22; full tinytest 2771
    + C3; air/lintr on R/bartcore.R.

C2. Restricted-moderator exact-posterior gate (benchmarks/R). Extend
    the exact machinery so a tau forest restricted to a strict column
    subset is validated against the closed-form posterior. Scenario:
    TWO ordinal predictors, x1 over K1 = 3 cells and x2 over K2 = 2
    cells (6 crossed cells), a fixed 0/1 z balanced within each crossed
    cell, gaussian response; moderators = column 1 (tau may split on x1
    only), mu unrestricted (spans x1, x2). DGP gives mu a real
    dependence on BOTH predictors and tau a dependence on x1 only, so
    the restriction is correctly specified and the restricted posterior
    is well-defined and distinguishable from the full one. Exact side:
    reuse bcf-exact.R verbatim for the tau forest (its 1-D
    contiguous-cell enumerate() over x1), and extend it with a bounded
    2-D axis-aligned single-tree enumerator for the UNRESTRICTED mu
    forest over the 6-cell grid; the joint tree space is the product.
    Downstream is unchanged: sufficient statistics key on the
    (mu-cell in the 6 crossed cells, tau-cell in the 3 x1 cells, z)
    stratum, and the closed-form Gaussian block marginal, the whitened
    eigendecomposition vectorized over the sigma grid, the chisq sigma
    calibration, and the deliverable-1 leaf-scale map carry over
    directly. Statistic and mode: mode 1 (glue fixed a = 1, b0 = 0,
    b1 = 1) matching E[mu] over the 6 cells and E[tau] over the 3 x1
    cells - the free-glue modes (2a, 2b) are already validated
    full-store by bcf-exact.R and need not repeat under restriction.
    Tolerance/seed/thinning follow bcf-exact.R's quick/full split
    (quick tol 0.05, full 0.015; nSeeds 1/3), same max-gap report and
    exit-status style as categorical-exact.R. Deliver as a sibling
    benchmarks/R/bcf-exact-restricted.R (sourcing shared helpers where
    clean) rather than growing bcf-exact.R, keeping the full-store gate
    unmoved. This is the gate for the restriction: a run that uses
    moderators is not equivalence-gated, so its correctness is proven
    here. rng: posterior-changing (uses the restriction) - not
    equivalence-gated; the gate IS the exact match. Abort: a restricted
    E[tau]/E[mu] gap over tolerance = the mask, the backfit under
    restriction, or the R resolution is wrong.

C3. Tests (inst/tinytest/test-bcf.R additions). (a) Validation errors
    from the C1 argument: unknown moderator name, out-of-range index,
    empty vector, names-without-colnames - each a targeted
    expect_error. (b) A restricted run is accepted and stays sane:
    build with moderators = a subset, run, assert finite train/sigma
    and moving per-forest fits (bartcoreForestFits), matching the
    existing sanity style; this does NOT re-prove the posterior (C2
    owns that) - it proves the R plumbing carries the restriction
    without error. (c) Neutrality of the default from R:
    moderators = NULL reproduces the no-argument run's fits bitwise
    (a fixed-seed pair), a cheap R-side echo of the equivalence gate.
    (d) Containment observable from R: see Q-C - if the forest-
    addressed getTrees/varcount is taken, assert the tau forest's
    reported split variables are a subset of the moderator set
    (a hard containment check); otherwise this is left to the existing
    tests/cpp containment test and (b)'s sanity, noted in C4. rng:
    neutral for (a)/(c); (b)/(d) exercise a restricted posterior but
    assert only finiteness/containment, not draw values, so they add
    no snapshot. Gate: tinytest 2771 + these, no regen.

C4. Docs + landing note. No NEW man/ topic and no _pkgdown.yml change:
    bartcoreBCFSampler is unexported and undocumented, so the argument
    is documented only in its in-file roxygen/comment header; the same
    holds for bcf-exact-restricted.R (a benchmark, not a shipped
    function). If and when a public dbarts-level BCF wrapper is added
    (deferred, public-surface.md section 5), THAT gains an Rd topic and
    the moderators argument, which then triggers the _pkgdown.yml entry
    + pkgdown::check_pkgdown requirement - flag it there, not here.
    One caveat: if Q-C is taken AND the forest selector surfaces on the
    public getTrees/varcount methods (not an internal helper), their
    EXISTING Rd topic gains the argument line - still no new topic, so
    still no _pkgdown.yml work, but the Rd edit and R CMD check man do
    apply. Append a landing note recording the resolved open questions,
    the restricted-gate scenario, and (if Q-C deferred) that R-side tau
    containment leans on tests/cpp. Gate: R CMD check man (unaffected
    unless the Q-C Rd caveat applies); full tinytest 2771 + C3.

### Verification

- equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds at every
  commit - the default BCF and every single-forest path are untouched;
  any deviation on the default path is a defect, stop.
- Full tinytest 2771 + the C3 additions from a preclean install if any
  header changed (Q-C's forest query would; C1 alone would not);
  existing snapshots UNREGENERATED (moderators defaults to NULL, so no
  default draw moves).
- bcf-exact.R quick + full UNMOVED (the full-store gate is untouched);
  the new bcf-exact-restricted.R quick + full pass to MC error.
- tests/cpp: the existing tau-containment test still passes; if Q-C is
  taken, delete stale binaries after the header edit.
- No bench-sampler run (see the bench posture note above); dbarts.h
  unchanged, no stan4bart lockstep (the growth is R + internal .Call
  only).
- air format + lintr on touched R files; rchk on the next scheduled
  run only if Q-C touches the bridge.

### Open questions

- Q-A (moderators as names, indices, or both). RECOMMEND BOTH: accept
  a character vector resolved against colnames(sampler$data@x) OR a
  1-based integer vector, resolving to indices R-side before the
  bridge. Names are the ergonomic surface bartCause wants and indices
  are the unambiguous fallback when a matrix is unnamed; the cost is a
  few lines of match()/range validation, all safe-over-fast R with no
  engine impact. NOTE on factor expansion: the data reaching
  bartcoreBCFSampler is already the expanded numeric design (the caller
  ran makeModelMatrixFromDataFrame), so a factor moderator is named or
  indexed by its EXPANDED dummy columns, not the original variable -
  original-variable-name expansion belongs to a future public wrapper
  (or bartCause) that still owns the formula/term map, not to this
  internal function.
- Q-B (expose a symmetric mu restriction now). RECOMMEND keep mu
  full-store for now; wire only tau's moderators. The engine is already
  symmetric (BCFForestSpec.columns exists for both forests, only the
  bridge wiring is tau-only), so adding a prognostic-columns knob later
  is cheap and neutral, and no consumer (bartCause, the exact gate)
  needs a restricted mu today - the C2 gate keeps mu unrestricted,
  matching the shipped default via a bounded 2-D enumerator. What would
  change it: a consumer that wants mu confined (e.g. a partially-linear
  BCF variant), at which point expose it in lockstep.
- Q-C (creation-only, or add a post-creation query for tau's splits).
  RECOMMEND add a `forest` selector to the getTrees/varcount query
  (thread a forest index through bartcore_getTrees into the bridge
  getTrees loop, defaulting to forest 0 so single-forest callers are
  unchanged), so the tau forest's split variables become R-observable
  and C3(d) is a HARD containment check rather than a statistical
  proxy; the moderators argument itself stays creation-time only (it is
  structural). What would change it: if the tests/cpp containment test
  is deemed sufficient and VD wants zero new bridge surface, defer the
  selector and let C3 lean on tests/cpp + the C2 gate for containment.

### Drift between the plan text and the code

- The step-3 interim "both forests read the full store (the moderator
  subset waits on data-ownership's views)" (recorded in Phase 1's
  Landing above, chain.hpp:471) is DISCHARGED: data-ownership-4 step 2
  (4e1fb5b) landed the mask and the moderatorsExpr slot. This phase
  consumes it; the earlier deferral no longer holds.
- Gate baselines moved since step 5: the exactness/neutrality gates
  now read equivalence-ac6ec2c.rds (22/22) not equivalence-235bebc
  (18/18), tinytest 2771 not 2497, and bench-sampler-4008675.csv not a
  235bebc baseline. This phase uses the current baselines throughout.
- bcf-ridge-interweaving.md and bcf-exact-weak.R landed after step 5
  and are not named in this plan's original body; they are recorded in
  the Context above as inherited, not reopened.

### Landing

C1 landed (8552b68): `moderators` on the internal bartcoreBCFSampler.
Q-A resolved BOTH names and 1-based indices, resolved R-side (match()
against colnames(data@x), range-checked, sorted-unique integer to the
bridge); factor moderators address EXPANDED dummy columns (original-
variable expansion stays with a future public wrapper). The comment
header carries the Hahn/Murray/Carvalho (2020) canonical usage note
(propensity in the design, out of moderators). C3(a)-(c) tinytest
additions landed with it (validation errors, restricted-run sanity via
per-forest fits, fixed-seed default-neutrality echo); tinytest
2771 -> 2782. Equivalence 22/22 identical.

Q-C landed (f39f335): forest selector on the tree query. Internal
bartcoreGetTrees gained forest = 0L (0-based, matching
bartcoreForestFits); savedTree/savedTreeSlopes/savedTreeMasks/
flattenTree threaded with a defaulted trailing forestIndex through
chain/sampler/facade plus a new read-only numTreesInForest virtual;
DEF_FUNC arity 7 -> 8. The PUBLIC R5 getTrees signature, its Rd topic,
and dbarts.h are all UNCHANGED (both pass forest 0). C3(d) landed with
it: hard tau containment from R (tau splits within the moderator set,
mu proven to split outside it); tinytest 2782 -> 2786. Equivalence
22/22 identical; component tests green.

Q-B resolved: mu stays full-store; the engine is symmetric
(BCFForestSpec.columns exists on both specs) so a prognostic
restriction is cheap later if a consumer needs it; only tau is
bridge-wired.

C2 landed (f6804f1): benchmarks/R/bcf-exact-restricted.R - two
ordinal predictors (3 x 2 crossed cells), tau restricted to x1, mu
unrestricted via a new bounded 2-D axis-aligned single-tree
enumerator (rectangle recursion weighted exactly as CGMTreePrior:
uniform over available variables then uniform over that variable's
cuts; enumeration weights asserted to sum to 1), tau reusing the 1-D
machinery; mode-1 fixed glue; sufficient stats keyed on the
(mu-cell, tau-cell, z) stratum. Results: quick (tol 0.05) E[mu] gap
0.0003 / E[tau] gap 0.0002; full (tol 0.015, 3 seeds) gaps
0.0002/0.0002 - the restricted posterior matches the closed form to
MC error, proving the mask, the backfit under restriction, and the R
resolution end to end. A per-seed containment stopifnot (tau fit
constant across x2 within each x1 cell) guards the mask directly.
Harness note recorded: x2 has 2 levels so the script sets per-column
n.cuts = c(K1-1, K2-1) - the scalar default would place degenerate
cuts on x2 (data.hpp numCuts behavior for ordinal columns without
quantiles) and corrupt the enumerator's prior weights. bcf-exact.R
itself untouched (full-store gate unmoved).

Close: no bench-sampler run owed (per the bench posture note above);
no Rd/_pkgdown.yml change anywhere (all surface internal; public
getTrees unchanged); dbarts.h frozen throughout; Phase 2 is COMPLETE.
