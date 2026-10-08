# small-rulings-1008-engine: six engine rulings of 2026-10-08

Status: PLANNED 2026-10-08 (dec-B329, dec-B356, dec-B365, dec-B366, dec-B369, dec-B386).

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a hazard fit of several forests whose vars omit period and for a constant
gaussian or aft response; at the call only (shifting) after setData or a warm start onto or off a
single-valued leaf covariate, and after setLeafPrior moves k.scale under a drawn k; NEUTRAL, bit for bit,
for linear and gp leaves at 8 or fewer covariates and for the perturb window at its default of 1.
window: before 1.0-0; engine slices stay serial.
budget: ~850 lines (code ~330, tests ~380, help and docs ~140).

## Goal

The six rulings are built: slopes are converted by the formula on both sides of a single-valued
covariate, every forest of a hazard fit may split on period, leaf covariates have no cap, the perturb
window is a control setting, setLeafPrior keeps the spread in force, and a constant response's window is
centred on its value. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Items

### 1. leaf-slope-formula-single-value (dec-B329)

Ruling: at setData and a warm start a live linear-leaf slope is converted by the formula whether the
covariate holds one value before the call, after it, or neither; where nothing is to convert the leaf stays.

Change: [`convertLinearCoefficients`](../../src/bartcore/model.hpp) loses `guardNoSpread` and runs the formula
on every column whose centre, scale or mark differs (a single-valued side divides by its placeholder 1, as
the saved-draw path already does). [`convertFlatLinearTree`](../../src/bartcore/chain.hpp) loses the argument;
its callers in [`applyNewData`](../../src/bartcore/chain.hpp) and
[`convertDonorStandardization`](../../src/bartcore/chain.hpp) drop it. The doc comment above
[`LeafStandardization`](../../src/bartcore/model.hpp) and the conversion's comment are rewritten.

Draws: only a sampler with a linear leaf whose covariate holds one value on one side of setData or
installTrees, from that call on. No equivalence scenario reaches it.

Tests: [`testLinearCoefficientConversion`](../../tests/cpp/test_model.cpp) and
[`testLinearLeafSetDataConversion`](../../tests/cpp/test_moves.cpp) assert the formula in all four cases (spread or one value, before and
after);
in [test-leaf-conversions.R](../../inst/tinytest/test-leaf-conversions.R) the block after
["setData onto a covariate with spread"](../../inst/tinytest/test-leaf-conversions.R) asserts slopes equal to
the old slope times the new scale and the intercept shifted by the formula, live and after installTrees,
in place of zero slopes and unmoved fits; the kept draws stay as they are.

Help: the newData and installTrees items of [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) say
the formula applies in every case and the fit may move until the sampler runs; design note in
[leaf-conversions.md](../design/leaf-conversions.md).

Size: ~110 lines.

### 2. hazard-period-every-forest (dec-B365)

Ruling: every forest of a hazard fit may split on its period column whatever its vars; a forest constant
over time has no spelling in 1.0-0.

Change: in [spec.R](../../R/spec.R) the union with the period column (`ncol(data@x)`) moves from the
`singleForest` branch to every forest's resolved columns: `firstColumns` loses the `singleForest` condition
and `forestColumns` unions period into each restricted later forest, before the blocks and interaction
constraints are resolved against them. An unrestricted forest stays NULL. The record in the control's
`bartcore.forests` carries the union, so a copy and a reload agree.

Draws: a hazard fit of several forests whose vars omit period. No harness scenario has one.

Tests: in [test-single-forest-vars.R](../../inst/tinytest/test-single-forest-vars.R) the block
["with several forests each forest's 'vars' is taken as written"](../../inst/tinytest/test-single-forest-vars.R)
now asserts period splits in the first forest whether or not vars names it, and the restricted
caller columns still report none.

Help: the vars item of [forest.Rd](../../man/forest.Rd) drops the several-forests exception; design note in
[single-forest-vars.md](../design/single-forest-vars.md).

Size: ~40 lines.

### 3. leaf-covariate-cap-removed (dec-B366)

Ruling: linear and gp leaves take any number of covariates, the scratch sized at creation; the R refusal
and the limits-table row go; no draw changes at 8 or fewer.

Change: [`LinearGaussianLeaf::maxNumCovariates`](../../src/bartcore/model.hpp) and
[`maxFunctionLeafCovariates`](../../src/bartcore/tree.hpp) go.
- LinearGaussianLeaf holds `mutable` vectors sized `(q + 1)^2`, `q + 1` and `q + 1` at the point it learns q,
  replacing the stack arrays in `logIntegratedLikelihoodForNode`, `drawFromPosteriorForNode` and
  [`accumulateNodeStatistics`](../../src/bartcore/model.hpp). They are used only on the chain's own thread,
  as the statistics cache already is; the test-row pool calls `fitForTestObservation`, which needs none.
- [`CachedNodeStatistics`](../../src/bartcore/model.hpp) holds `p * p` doubles in place of a fixed 81, and
  `statisticsCacheUsedBytes_` counts them beside the member lists, so the budget still bounds the cache
  at large q.
- [`addFlatFunctionPredictionsBelow`](../../src/bartcore/tree.hpp)'s outer overload allocates one `uStar`
  of q doubles per call and passes it down.
- [`leafCovariateDesignationIsValid`](../../src/bartcore/facade.hpp) drops the count test; the other
  designation rules stay.
- [`resolveLeafCovariates`](../../R/model.R) drops "at most 8 leaf covariates are supported".

Draws: none. The arithmetic order is unchanged, so q of 8 or fewer is bitwise, and q above 8 is newly reachable.

Tests: tests/cpp: the linear marginal and posterior draw at q = 9 and 16 against an independent dense
computation (log determinant and solve by plain elimination) within 1e-10 relative, and the cache's byte
count at q = 16. tinytest ([test-linear-leaves.R](../../inst/tinytest/test-linear-leaves.R),
[test-gp-leaves.R](../../inst/tinytest/test-gp-leaves.R)): 9 and 12 columns fit with finite draws, and with
keepTrees the replay at the training rows equals the train fits to 1e-12.

Help: the leaf-regression row leaves the limits table of [dbartsControl.Rd](../../man/dbartsControl.Rd), whose
count of fixed values drops by one; [engine-constants.md](../design/engine-constants.md) and
[linear-leaves.md](../design/linear-leaves.md); the header comment of
[constant-linear-leaf-covariates.R](../../benchmarks/R/constant-linear-leaf-covariates.R).

Size: ~220 lines.

### 4. perturb-width-control (dec-B366)

Ruling: the perturb window is a dbartsControl setting beside proposal.probs, a whole number of grid
positions, default 1, which changes no draw.

Change: the value is today's half-width, [`perturbWidth`](../../src/bartcore/moves.hpp): a displacement of
up to w positions either side, clipped by the node's interval, with no upper bound.
- R: a control slot and argument after proposal.probs (name: Open calls), held to a whole number of at
  least 1 by the class validity in [A_class.R](../../R/A_class.R) and refused as fractional by
  [`dbartsControl`](../../R/dbarts.R) with its other whole-number arguments.
- Bridge: [`parseProposalProbs`](../../src/R_interface_bartcore.cpp) reads it into the parsed model beside
  the mixture.
- Engine: an `int32_t perturbWidth` follows `perturbProbability` through
  [`SamplerOptions`](../../src/bartcore/chain.hpp), [`ModelParameters`](../../src/bartcore/chain.hpp), the
  forest and [`ForestStructureSpec`](../../src/bartcore/combiner.hpp) into
  [`MoveContext`](../../src/bartcore/moves.hpp), and [`perturbMove`](../../src/bartcore/moves.hpp) reads
  `ctx.perturbWidth`. moves.hpp keeps the default as a named constant, so the design cites stay live.
- The verbose mixture line is unchanged.

Draws: none at 1; a positive perturb probability at another width proposes over that window.

Tests: [`testPerturbMove`](../../tests/cpp/test_moves.cpp) at widths 1 and 3, no accepted move farther than
the width; tinytest beside [test-proposal-probs.R](../../inst/tinytest/test-proposal-probs.R): the default
and an explicit 1 draw bitwise as the control without the argument; 0, -1, 1.5, NA and a length-2 value are
refused by name; width 3 runs and is kept by copy, reload and setControl as Open calls settle.
[perturb-balance.R](../../benchmarks/R/perturb-balance.R) takes the width from an argument (default 1) and
exact-gates.yaml adds a width-3 arm.

Help: [dbartsControl.Rd](../../man/dbartsControl.Rd) usage, a new argument, the whole-number list, the
perturb sentence of proposal.probs, and the limits row (settable); NEWS's tree-moves entry gains the
argument; [engine-constants.md](../design/engine-constants.md) and [perturb-move.md](../design/perturb-move.md)
(window width); [mutation-battery.R](../../benchmarks/R/mutation-battery.R)'s perturb mutants are re-matched
to the new source text, and [constant-perturb-width.R](../../benchmarks/R/constant-perturb-width.R) sets the
width by control in place of a private build.

Size: ~250 lines.

### 5. leaf-prior-spelling-keeps-spread (dec-B356, dec-B369)

Ruling: a setLeafPrior that changes a drawn leaf spread's hyperprior, its spelling or its invchi() scale
keeps each chain's spread in force at the call, k becoming k_old x k.scale_new / k.scale_old; a fixed sd
stated sets the spread.

Change: an engine primitive, `Chain::scaleDrawnK(f, factor)`, beside
[`setForestFixedK`](../../src/bartcore/chain.hpp): on a forest that draws k, multiply k by factor; a
factor of exactly 1 writes nothing. It runs per chain through the facade and
[`Sampler`](../../src/bartcore/sampler.hpp), and the bridge gets one new entry. In R the setLeafPrior method,
and not [`writeLeafPrior`](../../R/dbarts.R), reads k.scale before the write, writes, reads it after, and
calls the primitive with after / before. Both reads are under one transform. So
[`reissueNamedLeafSd`](../../R/dbarts.R) after a re-anchor, and $setModel, are untouched. Under the k
spelling a changed chi() leaves k.scale alone, the factor is 1, and nothing moves. The primitive does not
touch the k draw.

Draws: a single-forest sampler with a drawn k, from a setLeafPrior that moves k.scale. Multinomial (fixed
k) and amplitude forests (k pinned at 1) are untouched. No scenario calls setLeafPrior.

Tests: [test-calibration-midchain.R](../../inst/tinytest/test-calibration-midchain.R), after 200 sweeps on
two chains: chi to invchi and back, and invchi scale 2 to 4, keep `k.scale / getK()` per chain to 1e-14;
chi(1.5, 2) to chi(1.5, 4) leaves getK identical; a drawn k to a fixed sd sets the stated spread (kept);
setResponse(updateScale = TRUE) under a drawn invchi() leaves getK identical (new); the "touches nothing
else" pins stay as they are.

Help: the setLeafPrior docstring in [dbarts.R](../../R/dbarts.R) (rc-codoc) and its paragraph in
[dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) replace "keeps its current k ... jump" with the
rule; design note in [prior-defaults.md](../design/prior-defaults.md), beside the named-sd reissue.

Size: ~130 lines.

### 6. constant-response-window (dec-B386)

Ruling: a constant response's transform is c - 0.5 to c + 0.5 at creation, at a re-derivation, and
when a stored (c, c) is read.

Change: [`GaussianResponse::readRange`](../../src/bartcore/model.hpp) and
[`GaussianResponse::restoreScale`](../../src/bartcore/model.hpp) take a zero-width range as the low end
c - 0.5 with range exactly 1 (set, not computed), so the working response is 0 and the shift c.
[`Chain::unitsOf`](../../src/bartcore/chain.hpp) reads (c, c) as multiplier 1, shift c; the comments of
[`carriesUnits`](../../src/bartcore/chain.hpp) and [`convertStateUnits`](../../src/bartcore/chain.hpp) follow.
What is recorded stays (c, c) (Open calls). The count family's (c, c + 1) is a separate encoding and is
untouched.

Draws: fits on a constant gaussian or aft response, which already warn. A state stored before this lands
on a constant response has its function shifted by half a unit at install, which the ruling's "when a
stored (c, c) is read" decides.

Tests: tinytest ([test-boundary-inputs.R](../../inst/tinytest/test-boundary-inputs.R)): on rep(5, n),
getLeafPrior()'s prior.mean and response.shift are 5; [test-state-not-model.R](../../inst/tinytest/test-state-not-model.R)'s
["a constant response's transform is units too"](../../inst/tinytest/test-state-not-model.R) keeps units
c(2, 2) and its comment changes; tests/cpp's constant-pair case in
[`testStateRoundTripScaledOffset`](../../tests/cpp/test_state.cpp) pairs (0, 0) with a response spanning
(-0.5, 0.5) in place of (0, 1).

Help: [state-not-model.md](../design/state-not-model.md)'s "spans 1 upward" sentence; one NEWS line under
changes from 0.9-34 (a fit on a constant response is centred on it, not half a unit above).

Size: ~100 lines.

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
- Hot path (linear-leaf scratch, the perturb read): bench-sampler.R compare and a same-machine A/B of
  linear fits at q = 1, 2, 4, 8 against the base build, quiet machine, maintainer-run.
- stan4bart's and bairrtt's suites against the build.

Reviewer's mutants, each of which must fail a test: the guard restored in item 1; period unioned into
the first forest only; the scratch sized q rather than q + 1, and the cache bytes left out of the
budget; the width ignored in `perturbMove`; the k re-expression dropped, or applied in
`reissueNamedLeafSd`; the low end left at c.

## Expected verdicts

Against the MANIFEST: equivalence-1b7d730c 55 of 55, bcf-equivalence-1b7d730c 15 of 15 and
multinomial-equivalence-80b1c8d4 11 of 11 "identical draws (same RNG stream)", with no |z| line, so no
re-record and no MANIFEST change; the four snapshot files pass unchanged; every exact gate passes;
benchmarks no arm past 1.05. The draws that move are those the items name, and only tinytest sees them.

## Stop conditions

Stop and report, without working around, when:
- the diff passes ~1300 lines, or one item passes 1.5x its size;
- any equivalence scenario or snapshot moves (no re-record in this slice);
- the small-q A/B shows an arm past 1.05 that a reserved scratch does not recover;
- a ruling's reading needs a state-format change, or an Open call is answered otherwise than recommended
  in a way that changes the plan.

## Interactions

- The uncommitted k-move prototype (probit-k-mixing, dec-B371) touches chain.hpp. This slice's chain.hpp
  hunks are the perturb field beside `perturbProbability`, `scaleDrawnK` beside `setForestFixedK`,
  `unitsOf`, and the two conversion call sites. None is in the k draw or the forest's k fields. Engine
  slices are serial: whichever lands second rebases, and if the k move lands first its re-recorded
  baselines replace the ones named above. There is no semantic overlap: that move changes how k is drawn,
  and this slice writes k once at a call.
- state-install-keeps-spread (dec-B384, with the install surface's first slice) applies the same
  re-expression at setState, copy and reload. To keep the two apart, this slice leaves the install paths,
  the state format, the setState, copy and reload help, and getLeafPrior's "a state install leaves the
  prior" sentence alone. It leaves `scaleDrawnK` as the one primitive that slice calls with
  k.scale_recipient / k.scale_state, and the TODO entry gains that pointer.
- response-scale-rows slice C (dec-B364) keeps the scale in force when the rows in hold one value, so
  after it the window of item 6 is reached at creation, at setData and from a stored pair; no shared code
  beyond `readRange`.
- leaf-covariate-single-value-off-centre (with dec-B304) adds a floor to the standardization that item 1's
  formula reads; it is later and touches test-leaf-conversions.R only.
- This sitting's R-surface slice may land while this one's CI runs; R/dbarts.R and dbartsSampler-class.Rd
  conflicts are textual.

## Open calls

1. Item 1, a covariate single-valued on both sides at different values (setData from c1 to c2): the
   formula, shifting each intercept by slope x (c2 - c1), as dec-B328 converts a column setPredictor moved
   (recommended), or leave the leaf as today. The formula is the identity where the value is the same,
   which reads "nothing is to convert".
2. Item 5, a fixed k or sd to a drawn prior: keep the spread in force, k re-expressed (recommended, the
   leaves having been drawn under it), or keep k as today.
3. Item 5, scope: setLeafPrior only (recommended, as ruled), leaving $setModel's install and the named-sd
   reissue after a re-anchor as they are; the alternative extends it to $setModel.
4. Item 4, the name: `perturbWidth`, like the engine settings it joins in the limits table and the
   engine's own name (recommended), or `perturb.width`, like proposal.probs.
5. Item 4, setControl on a live sampler: takes a changed width, installed with the mixture as
   proposal.probs is (recommended, a proposal setting that moves no stationary distribution), or refuses
   it as the creation-fixed engine limits are refused.
6. Item 6, what is recorded: the state's fit.scale and the model's response.range stay (c, c)
   (recommended: the count family's design notes that a midrange pair does not round-trip exactly, and
   test-state-not-model.R pins c(2, 2)), or the window (c - 0.5, c + 0.5) is recorded.
