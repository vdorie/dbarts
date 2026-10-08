# small-rulings-1008-engine: six engine rulings of 2026-10-08

Status: PLANNED 2026-10-08 (dec-B329, dec-B356, dec-B365, dec-B366, dec-B369, dec-B386, dec-B391); one
call pending the maintainer (Open calls).

agent: opus implementer, one; blind critique taken 2026-10-08; one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a hazard fit of several forests whose vars omit period and for a constant
gaussian or aft response; at the call only (shifting) after setData or a warm start that takes a
leaf covariate between one value and spread, and after setLeafPrior moves k.scale into a drawn prior;
NEUTRAL, bit for bit, for linear and gp leaves at 8 or fewer covariates and for n.perturb.cuts = 1.
window: before 1.0-0; engine slices stay serial.
budget: ~950 lines (code ~380, tests ~420, help and docs ~150).

## Goal

The six rulings are built: a slope is converted by the formula where a leaf covariate holds one value on
one side of a conversion, every forest of a hazard fit may split on period, leaf covariates have no cap,
the perturb window is the control setting n.perturb.cuts, setLeafPrior keeps the spread in force, and a
constant response's window is centred on its value. The tier is "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Items

### 1. leaf-slope-formula-single-value (dec-B329)

Ruling: at setData and a warm start a live linear-leaf slope is converted by the formula where the
covariate holds one value on one side; where it holds one value on both sides the leaf stays as built,
c1 to c2 included.

Change: [`convertLinearCoefficients`](../../src/bartcore/model.hpp) loses `guardNoSpread`. It keeps the
skip where neither side has spread and runs the formula otherwise. A single-valued side divides by its
placeholder 1, as the saved-draw path already does. [`convertFlatLinearTree`](../../src/bartcore/chain.hpp)
loses the argument. Its warm-start caller in
[`convertDonorStandardization`](../../src/bartcore/chain.hpp) drops it, and so does the direct live call to
`convertLinearCoefficients` in [`applyNewData`](../../src/bartcore/chain.hpp). The comments above
[`LeafStandardization`](../../src/bartcore/model.hpp) and on the conversion are rewritten.

Draws: a sampler with a linear leaf whose covariate goes between one value and spread at setData or
installTrees, from that call on. No equivalence scenario reaches it.

Tests:
- [`testLinearCoefficientConversion`](../../tests/cpp/test_model.cpp) and
  [`testLinearLeafSetDataConversion`](../../tests/cpp/test_moves.cpp) assert the formula in the two
  one-sided cases and keep the both-sides case.
- [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R) changes, asserting slopes and
  intercepts by the formula in place of zero slopes and unmoved fits:
  - section 1's block after ["setData onto a covariate with spread"](../../inst/tinytest/test-leaf-conversions.R);
  - section 6's blocks ["values made constant after creation"](../../inst/tinytest/test-leaf-conversions.R)
    (held) and ["rows appended that give a constant covariate spread"](../../inst/tinytest/test-leaf-conversions.R)
    (grown);
  - section 6's header comment.
- Kept, unchanged: section 6's ["neither side has spread"](../../inst/tinytest/test-leaf-conversions.R)
  block and every kept-draw pin.

Help: the newData and installTrees items of [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd):
with one value on one side the slope is converted by the formula and the fit may move until the sampler
runs; on both sides the coefficients are left as they are (kept). Design note in
[leaf-conversions.md](../design/leaf-conversions.md).

Size: ~140 lines.

### 2. hazard-period-every-forest (dec-B365)

Ruling: every forest of a hazard fit may split on its period column whatever its vars; a forest constant
over time has no spelling in 1.0-0.

Change: in [spec.R](../../R/spec.R) the union with the period column (`ncol(data@x)`) moves from the
`singleForest` branch to every forest's resolved columns. `firstColumns` loses the `singleForest`
condition, and `forestColumns` unions period into each restricted later forest. Both happen before the
blocks and interaction constraints are resolved against them. An unrestricted forest stays NULL. The record
in the control's `bartcore.forests` carries the union, so a copy and a reload agree.

Draws: a hazard fit of several forests whose vars omit period. No harness scenario has one.

Tests: in [test-single-forest-vars.R](../../inst/tinytest/test-single-forest-vars.R) the block
["with several forests each forest's 'vars' is taken as written"](../../inst/tinytest/test-single-forest-vars.R)
now asserts period splits in the first forest whether or not vars names it. Restricted caller columns
still report none. A second fit restricts the basis forest too, `forest(vars = "a", basis = ...)`, and
asserts period splits there.

Help: the vars item of [forest.Rd](../../man/forest.Rd) drops the several-forests exception; design note in
[single-forest-vars.md](../design/single-forest-vars.md).

Size: ~50 lines.

### 3. leaf-covariate-cap-removed (dec-B366)

Ruling: linear and gp leaves take any number of covariates, the scratch sized at creation; the R refusal
and the limits-table row go; no draw changes at 8 or fewer.

Change: retired: [`LinearGaussianLeaf::maxNumCovariates`](../../src/bartcore/model.hpp),
retired: [`GPGaussianLeaf::maxNumCovariates`](../../src/bartcore/model.hpp) and
retired: [`maxFunctionLeafCovariates`](../../src/bartcore/tree.hpp) go.
- LinearGaussianLeaf holds `mutable` vectors sized `(q + 1)^2`, `q + 1` and `q + 1` at the point it
  learns q. They replace the stack arrays in `logIntegratedLikelihoodForNode`, `drawFromPosteriorForNode`
  and [`accumulateNodeStatistics`](../../src/bartcore/model.hpp). They are used only on the chain's own
  thread, as the statistics cache already is; the test-row pool calls `fitForTestObservation`, which
  needs none.
- [`CachedNodeStatistics`](../../src/bartcore/model.hpp) holds `p * p` doubles in place of a fixed 81, and
  `statisticsCacheUsedBytes_` counts them beside the member lists, so the budget still bounds the cache
  at large q. Counting them caches fewer leaves where the budget binds, at any q: speed only, since a
  served value is bitwise a fresh scan's.
- [`addFlatFunctionPredictionsBelow`](../../src/bartcore/tree.hpp) takes its `uStar` scratch from the
  caller. The outer overload sizes one vector of q per call. `Chain::addFlatPredictions` calls the
  templated form directly, with a `uStar` vector added to [`PredictScratch`](../../src/bartcore/chain.hpp).
- [`leafCovariateDesignationIsValid`](../../src/bartcore/facade.hpp) drops the count test; the other
  designation rules stay.
- [`resolveLeafCovariates`](../../R/model.R) drops "at most 8 leaf covariates are supported".

Draws: none. The arithmetic order is unchanged, so q of 8 or fewer is bitwise, and q above 8 is newly reachable.

Tests: tests/cpp: the linear marginal and posterior draw at q = 9 and 16 against an independent dense
computation (log determinant and solve by plain elimination) within 1e-10 relative, and the cache's byte
count at q = 16. tinytest ([test-linear-leaves.R](../../inst/tinytest/test-linear-leaves.R),
[test-gp-leaves.R](../../inst/tinytest/test-gp-leaves.R)): 9 and 12 columns fit with finite draws, and with
keepTrees the replay at the training rows equals the train fits to 1e-12.

Help and docs:
- The leaf-regression row leaves the limits table of [dbartsControl.Rd](../../man/dbartsControl.Rd); the
  table's prose is item 4's.
- [engine-constants.md](../design/engine-constants.md) marks its cite of the cap `retired:` and rewords
  its REFUSED rows; its verbatim cite of
  [constant-linear-leaf-covariates.R](../../benchmarks/R/constant-linear-leaf-covariates.R)'s header
  follows that header's new text.
- [engine-performance.md](engine-performance.md) marks its two cites `retired:`.
- [linear-leaves.md](../design/linear-leaves.md) ("at most 8", twice) and
  [gp-leaves.md](../design/gp-leaves.md) (`maxFunctionLeafCovariates = 8`) are reworded.
- [memory-footprint.md](../design/memory-footprint.md) ("an inline 81-double crossproduct") is reworded.

Size: ~250 lines.

### 4. perturb-width-control (dec-B366, dec-B391)

Ruling: the perturb window is a dbartsControl setting beside proposal.probs, named `n.perturb.cuts`, the
most cut positions a perturb proposal moves a split either way; default 1, which changes no draw.
setControl accepts a changed value between runs, installed with the mixture.

Change: the value is today's half-width, [`perturbWidth`](../../src/bartcore/moves.hpp): a displacement of
up to w positions either side, clipped to the node's interval.
- R: a slot and argument `n.perturb.cuts` (slot name equal to the argument, so `controlArgumentFromSlot`
  reproduces it), after proposal.probs. It is a double, so that it can hold `Inf`. The class validity in
  [A_class.R](../../R/A_class.R) and [`dbartsControl`](../../R/dbarts.R) take a positive whole number or
  `Inf`, and refuse 0, negatives, `NA`, non-whole values and a length other than 1 by name. Values above
  `.Machine$integer.max` are taken.
- R, setControl: `mixtureMoved` in [`setControl`](../../R/dbarts.R) compares `n.perturb.cuts` beside
  `proposal.probs`, so a changed width goes through the same model install and rollback.
- Bridge: [`parseProposalProbs`](../../src/R_interface_bartcore.cpp) reads it and clamps anything at or
  above the cut cap (65533) to the cap. The clamp is exact, the window being clipped to the node's
  interval, and it is the bound that keeps `current + width` inside `int32_t`. The bridge backstops the
  R checks with its own refusal. Four fill sites carry it: the sampler options, `forest.forest`,
  `spec.forest`, and the `ModelParameters` of setModel.
- Engine: an `int32_t perturbWidth` follows `perturbProbability` everywhere it goes:
  - the structs: [`SamplerOptions`](../../src/bartcore/chain.hpp),
    [`ModelParameters`](../../src/bartcore/chain.hpp), the forest,
    [`ForestStructureSpec`](../../src/bartcore/combiner.hpp) and
    [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp);
  - the variance forest's `vf` copy in chain.hpp;
  - into [`MoveContext`](../../src/bartcore/moves.hpp), where [`perturbMove`](../../src/bartcore/moves.hpp)
    reads `ctx.perturbWidth`.

  moves.hpp keeps `perturbWidth` as the named default, so the design cites stay live.
- The verbose mixture line is unchanged.

Draws: none at 1; a positive perturb probability at another width proposes over that window.

Tests:
- tests/cpp: [`testPerturbMove`](../../tests/cpp/test_moves.cpp) at widths 1 and 3, with no accepted
  move farther than the width.
- tinytest beside [test-proposal-probs.R](../../inst/tinytest/test-proposal-probs.R):
  - the default and an explicit 1 draw bitwise as the control without the argument;
  - 0, -1, 1.5, `NA` and a length-2 value are refused by name;
  - `Inf`, `.Machine$integer.max + 1` and 1e12 are accepted and draw bitwise as 65533;
  - width 3 runs on a single forest, a bcf, a multinomial and a variance-forest sampler;
  - width 3 is kept by copy and reload, and taken by setControl between runs, with a refused install
    rolling the control back.
- Exact gate: [perturb-balance.R](../../benchmarks/R/perturb-balance.R) takes the width from an argument
  (default 1), and exact-gates.yaml adds a width-3 arm.

Help and docs:
- [dbartsControl.Rd](../../man/dbartsControl.Rd):
  - usage, and the new argument: "the most cut positions a perturb proposal moves a split either way";
  - the whole-number list;
  - the perturb sentence of proposal.probs;
  - the engine-limits section after items 3 and 4: its count of fixed values, "the four that a workload
    can profitably move", the perturb row (settable), and the closing paragraph on which settings change
    what is sampled (n.perturb.cuts changes the proposal, not the posterior).
- NEWS's tree-moves entry gains the argument.
- [engine-constants.md](../design/engine-constants.md) and [perturb-move.md](../design/perturb-move.md)
  (window width).
- [benchmarks/README.md](../../benchmarks/README.md)'s "compile-time with no knob" paragraph.
- [mutation-battery.R](../../benchmarks/R/mutation-battery.R)'s perturb mutants are re-matched to the new
  source text.
- [constant-perturb-width.R](../../benchmarks/R/constant-perturb-width.R) sets the width by control in
  place of a private build.

Size: ~300 lines.

### 5. leaf-prior-spelling-keeps-spread (dec-B356, dec-B369)

Ruling: a setLeafPrior that changes a drawn leaf spread's hyperprior, its spelling or its invchi() scale
keeps each chain's spread in force at the call, k becoming k_old x k.scale_new / k.scale_old; a fixed sd
stated sets the spread. Planned here, pending the maintainer: the same on every switch into a drawn
prior, a fixed k or sd before it included.

Change:
- Engine: a new primitive, `Chain::scaleDrawnK(f, factor)`, beside
  [`setForestFixedK`](../../src/bartcore/chain.hpp). On a forest that draws k after the write, it
  multiplies k by factor; a factor of exactly 1 writes nothing. It writes `Forest::k`, the field the k-move
  prototype's scale move also writes. It runs per chain through the facade and
  [`Sampler`](../../src/bartcore/sampler.hpp), and the bridge gets one new entry.
- R: the setLeafPrior method, and not [`writeLeafPrior`](../../R/dbarts.R), reads k.scale before the
  write, writes, reads it after, and calls the primitive with after / before. Both reads are under one
  transform.
- Untouched: [`reissueNamedLeafSd`](../../R/dbarts.R) after a re-anchor, and $setModel. Under the k
  spelling a changed chi() leaves k.scale alone, the factor is 1, and nothing moves.

If the maintainer rules to keep k from a fixed prior instead, the R side calls the primitive only where
the forest drew k before the call, and the fixed-to-drawn test asserts k kept.

Draws: a single-forest sampler drawing k after a setLeafPrior that moves k.scale. Multinomial (fixed k)
and amplitude forests (k pinned at 1) are untouched. No scenario calls setLeafPrior.

Tests: [test-calibration-midchain.R](../../inst/tinytest/test-calibration-midchain.R), after 200 sweeps on
two chains:
- chi to invchi and back, invchi scale 2 to 4, and a fixed sd to invchi each keep `k.scale / getK()` per
  chain to 1e-14;
- chi(1.5, 2) to chi(1.5, 4) leaves getK identical;
- a drawn k to a fixed sd sets the stated spread (kept);
- setResponse(updateScale = TRUE) under a drawn invchi() leaves getK identical (new);
- the "touches nothing else" pins stay as they are.

Help: the setLeafPrior docstring in [dbarts.R](../../R/dbarts.R) (rc-codoc) and its paragraph in
[dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) replace "keeps its current k ... jump" with the
rule. They also say that the write keeps the chains' spread, so on an unrun sampler it is not the sampler
creation under the new prior would make. Design note in [prior-defaults.md](../design/prior-defaults.md),
beside the named-sd reissue.

Size: ~140 lines.

### 6. constant-response-window (dec-B386)

Ruling: a constant response's transform is c - 0.5 to c + 0.5 at creation, at a re-derivation, and
when a stored (c, c) is read. (c, c) stays the record.

Change: [`GaussianResponse`](../../src/bartcore/model.hpp) keeps `min_ = max_ = c`, so `getScale` reports
(c, c) to the state, `installsScale` and the sampler's anchor as today. The window is applied inside the
range-0 handling: a separate low end of the window, the range of 1 below it set rather than computed,
equal to `min_` where the range is positive and to c - 0.5 where it is 0. The low end is set in
[`GaussianResponse::readRange`](../../src/bartcore/model.hpp) and
[`GaussianResponse::restoreScale`](../../src/bartcore/model.hpp). It is read wherever the transform reads
`min_`:
- the working-response builds in `rescale`, `restoreScale` and `setOffset`'s pinned branch;
- `fitShift`;
- `computeLogLikelihood`.

`rescale`'s empty-response branch follows the same rule: (0, 0) recorded, low end -0.5.
[`Chain::unitsOf`](../../src/bartcore/chain.hpp) reads (c, c) as multiplier 1 and shift c; the comments of
[`carriesUnits`](../../src/bartcore/chain.hpp), [`convertStateUnits`](../../src/bartcore/chain.hpp) and
`getScale` follow. combiner.hpp only carries `fitMin` and `fitMax` in the state and needs no change. The
count family's (c, c + 1) is a separate encoding and is untouched.

Draws: fits on a constant gaussian or aft response, which already warn. A state stored before this lands
on a constant response has its function shifted by half a unit at install, which the ruling's "when a
stored (c, c) is read" decides.

Tests:
- tinytest: on rep(5, n) in [test-boundary-inputs.R](../../inst/tinytest/test-boundary-inputs.R),
  getLeafPrior()'s prior.mean and response.shift are 5.
- tinytest: [test-state-not-model.R](../../inst/tinytest/test-state-not-model.R)'s
  ["a constant response's transform is units too"](../../inst/tinytest/test-state-not-model.R) keeps units
  c(2, 2); its comment changes.
- tests/cpp: a constant GaussianResponse reports (c, c), shift c and a working response of 0, before and
  after `restoreScale(c, c)`.
- tests/cpp: the constant-pair case in [`testStateRoundTripScaledOffset`](../../tests/cpp/test_state.cpp)
  pairs (0, 0) with a response spanning (-0.5, 0.5) in place of (0, 1).

Help: [state-not-model.md](../design/state-not-model.md)'s "spans 1 upward" sentence; one NEWS line under
changes from 0.9-34 (a fit on a constant response is centred on it, not half a unit above).

Size: ~120 lines.

## Order of work

One commit per item, each with its tests/cpp and tinytest files green: 3, 4 (neutral; after them the
three equivalence harnesses must still be bitwise), then 1, 5, 2, 6. `--preclean` on every install;
facade.hpp gains a virtual in item 5.

## Gates

On the slice tip, against its own library, run independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing):
- tests/cpp, plain and under `-fsanitize=address,undefined`; the R-loaded ASAN path over the touched test
  files (q above 8 is new numerics).
- The full tinytest suite (`at_home = TRUE`).
- Reference build (`--enable-reference-build`, `--preclean`): equivalence.R, bcf-equivalence.R and
  multinomial-equivalence.R `compare --bitwise` (equivalence.R also `--strict-coverage`) against the
  MANIFEST's current files; the four arm64 reproducibility snapshot files.
- Every gate of exact-gates.yaml in quick mode on the shipped build, plus the new perturb-balance width arm.
- `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift and doc-freshness.
- Hot path (linear-leaf scratch, cache accounting, the perturb read), quiet machine, maintainer-run:
  bench-sampler.R compare, and a same-machine A/B of linear fits against the base build at q = 1, 2, 4
  and 8, at one n under the 256 MB cache budget and one where it binds.
- stan4bart's and bairrtt's suites against the build.

Reviewer's mutants, each of which must fail a test:
- the guard restored for one-sided cases in item 1, or the skip dropped where neither side has spread;
- period unioned into the first forest only;
- the scratch sized q rather than q + 1 (fails under the ASAN runs), and the cache bytes left out of the
  budget;
- the width ignored in `perturbMove`, the clamp dropped, and the width left out of `mixtureMoved`;
- the k re-expression dropped, or applied in `reissueNamedLeafSd`;
- the low end left at c, and `getScale` reporting the window.

## Expected verdicts

Against the MANIFEST: equivalence-1b7d730c 55 of 55, bcf-equivalence-1b7d730c 15 of 15 and
multinomial-equivalence-80b1c8d4 11 of 11 "identical draws (same RNG stream)", with no |z| line, so no
re-record and no MANIFEST change; the four snapshot files pass unchanged; every exact gate passes;
benchmarks no arm past 1.05. The draws that move are those the items name, and only tinytest sees them.

## Stop conditions

Stop and report, without working around, when:
- the diff passes ~1400 lines, or one item passes 1.5x its size;
- any equivalence scenario or snapshot moves (no re-record in this slice);
- the small-q A/B shows an arm past 1.05 that a reserved scratch does not recover;
- a ruling's reading needs a state-format change, or the pending call is ruled in a way this plan does
  not already state.

## Interactions

- The uncommitted k-move prototype (probit-k-mixing, dec-B371) touches chain.hpp. This slice's chain.hpp
  hunks are:
  - the width field beside `perturbProbability` and in the `vf` copy;
  - `scaleDrawnK` beside `setForestFixedK`;
  - `unitsOf`, `PredictScratch`, and the two conversion call sites.

  None is in the k draw or the forest's k fields. Both write `Forest::k`, the prototype in its scale move
  and this slice once at a setLeafPrior call; that move changes how k is drawn, and there is no semantic
  overlap. Engine slices are serial: whichever lands second rebases, and if the k move lands first its
  re-recorded baselines replace the ones named above.
- state-install-keeps-spread (dec-B384, with the install surface's first slice) applies the same
  re-expression at setState, copy and reload. This slice leaves alone the install paths, the state
  format, the setState, copy and reload help, and getLeafPrior's "a state install leaves the prior"
  sentence. `scaleDrawnK` is the one primitive that slice calls with k.scale_recipient / k.scale_state,
  and the TODO entry gains that pointer.
- response-scale-rows slice C (dec-B364) keeps the scale in force when the rows in hold one value, so
  after it the window of item 6 is reached at creation, at setData and from a stored pair; no shared code
  beyond `readRange`.
- leaf-covariate-single-value-off-centre (with dec-B304) adds a floor to the standardization that item 1's
  conversion reads; it is later and touches test-leaf-conversions.R only.
- This sitting's R-surface slice may land while this one's CI runs; R/dbarts.R and dbartsSampler-class.Rd
  conflicts are textual.

## Open calls

Settled by the maintainer's rulings or the coordinator on 2026-10-08:
- Item 1 follows dec-B329 as recorded: the leaf stays where the covariate holds one value on both sides.
- Item 5 applies to setLeafPrior only.
- The name is n.perturb.cuts (dec-B391), and setControl takes it.
- (c, c) stays the record.

Pending: item 5, a fixed k or sd switched to a drawn prior. Keeping the spread, as planned, is
recommended, the leaves having been drawn under it; the alternative keeps k, with the one-line change
stated under item 5.
