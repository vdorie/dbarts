# single-forest-vars: one declared forest honours its column restriction

Status: PLANNED.

agent: opus implementer, one; opus reviewer.
rng: POSTERIOR-CHANGING for a fit that states `vars` on a single declared forest: today the argument is
ignored and the fit draws the unrestricted model; afterwards it draws, seed for seed, what the same fit on
the named columns alone draws (given the same residual scale estimate). NEUTRAL for every other fit, whose
draws are bit for bit unchanged.
window: pre-release (dec-B241 makes one declared forest the ordinary way to state a single forest). Serial
with any other work in [`resolveSamplerSpec`](../../R/spec.R), the bridge file or chain.hpp.
budget: ~550 lines (C++ engine ~45, bridge ~25, R ~45, tests/cpp ~110, tinytest ~250, design note, manual
and TODO ~75), upper figure 900. Plans have run 1.5-2x low.

## Goal

`forest(vars = )` does on a single declared forest what it does on any forest of a multi-forest fit: the
forest splits on the named columns and no others, for the sampler's whole life. Such a fit is the fit on
those columns alone. A model whose restriction differs from the sampler's is refused by `setModel`.

## Context

- The defect. `dbarts(x, y, forests = list(forest(vars = c("a", "b"))))` on three columns splits on the
  third 10899 times in 300 sweeps (in every sweep) and is draw for draw the fit without `vars`; the same
  through `dbartsSpec`, by name or by index. The value is not even read: a name that is no column and an
  empty selection are both accepted. [`forest`](../../man/forest.Rd) says "Any forest may be restricted."
- Why. A multi-forest sampler carries each forest's columns in
  [`ForestStructureSpec`](../../src/bartcore/combiner.hpp), filled by
  [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) from the control's forests attribute, and
  [`Chain::buildSpecifiedForest`](../../src/bartcore/chain.hpp) turns them into the forest's `columnMask`,
  which every tree reads through [`Tree::collectAvailableVariables`](../../src/bartcore/tree.hpp). The
  variance forest has the same thing ([`applyVarianceAttributes`](../../src/R_interface_bartcore.cpp),
  [`Chain::buildVarianceForest`](../../src/bartcore/chain.hpp)). A single mean forest is built from
  [`SamplerOptions`](../../src/bartcore/chain.hpp), which has no such field, and
  [`resolveSamplerSpec`](../../R/spec.R) reads the first forest's `n.trees`, `interactions` and `blocks`
  but resolves `vars` only inside its multi-forest branch.
- The same argument works elsewhere, run: the first of two forests restricted to two of three columns
  never splits on the third (300 sweeps); `varianceForest(vars = )` and `variance = ~ a + b` confine the
  variance trees (0 splits on the third column in 300 sweeps against 5466 unrestricted). The variance
  forest has no defect and is not touched.
- Everything else one `forest()` can state for a single forest is honoured or refused by name already:
  `n.trees`, `base`, `power`, `interactions`, `blocks` reach the sampler; `sd`, `update.amplitude` and
  `amplitude.prior.variance` are refused ([`resolveForests`](../../R/model.R)). `vars` is the one dropped.
- `bart`, `xbart` and `rbart_vi` have no `forests` argument. A formula's `forest()` term always declares
  at least two forests. Neither is affected.

## Three routes

(i) The column mask, as on every other forest. (ii) Zeros in the tree prior's split probabilities for the
excluded columns, which needs no engine change. (iii) Refuse `vars` on a single forest by name.

From the user's side there is one answer: `forest(vars = )` on the only forest should mean what it means
on any forest, and what the manual says. That is (i) or (ii); (iii) keeps the manual false on the main road.

(ii) does not hold. Run:
- It confines a forest only while a positive-probability column is available at every node. With two
  allowed 0/1 columns at one split point each and a continuous excluded column, the zero-probability
  column is split on 31454 times in 2000 sweeps of 20 trees, in every sweep; the mask gives 0. When the
  allowed columns run out along a path, [`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp)
  returns the first available column whatever its probability. 0.9-34 does the same (32423).
- It makes structure ride on a parameter. dec-B254: a forest's columns are structure, split probabilities
  are a parameter `setModel` may change. `setModel` with a plain model on a forest confined by zeros is
  accepted and the forest splits on the excluded column again (10906 of 29680 splits). And a zero given to
  a column the trees use is accepted unchecked (dec-B254's finding, re-run: 5 such splits left after 850
  sweeps).
- DART has no split probabilities to zero, and a multi-forest fit refuses `split.probs` and DART outright,
  so the argument would have two mechanisms and a gap.
- The caller sees split probabilities they did not write on the stored tree prior.

(i) was prototyped outside the tree (about 100 added lines over the engine, the bridge and R) and run:
- 0 splits on the excluded columns in 1000 sweeps of 75 trees, where the unrestricted fit splits on them in
  every sweep; the same for probit, logistic, ordinal, negative binomial, aft, Student-t, linear and
  Gaussian-process leaves, a monotone column, every category forest of a multinomial fit; by name, by
  index, through `dbartsSpec` and the direct door; with a factor column; with allowed columns that run out
  (0 in 2000 sweeps).
- A hazard fit's design carries a `period` column the caller did not supply. In the prototype
  `vars = c("a", "b")` there excluded it, giving a hazard that is constant over periods; the slice keeps
  `period` allowed (Constraints).
- Draw for draw the fit on the sub-matrix at equal `sigest`: gaussian, probit, with `interactions`, with a
  caller's `split.probs` (the caller's ratios hold among the allowed columns).
- It composes with `blocks` (each block row is intersected with the mask,
  [`Chain::installBlockMasks`](../../src/bartcore/chain.hpp)), with `interactions`, with a variance forest
  on other columns.
- It survives `copy`, a save and reload, `setPredictor` (a column, the whole matrix), `setData`,
  `sampleTreesFromPrior` and `growFromRoot`. `setState` with a state that splits on an excluded column and
  `installTrees` from an unrestricted donor are refused with the messages a restricted variance forest
  gives ([`Chain::columnMaskStateFeasible`](../../src/bartcore/chain.hpp) reads the same mask).
- 20 kinds of fit that state no `vars` on a single forest are bitwise the tip's (families, DART,
  `split.probs`, `interactions`, `blocks`, linear and monotone leaves, a restricted variance forest, two
  forests with and without `vars`, multinomial, `bart`, `xbart`); the full tinytest suite (13897 results)
  and tests/cpp pass unchanged.
- Three things it showed that the steps settle. `setModel` accepts a model with another restriction, in
  both directions, and a copy then disagrees with the live sampler. `split.probs` zero on every allowed
  column gives a forest that splits on the first allowed column only. Under DART the Dirichlet is laid over
  every column, the excluded ones included: not the fit on the sub-matrix, and `varprobs` reports
  probabilities for columns that cannot be split on. A second prototype that lays it over the allowed
  columns is draw for draw DART on the sub-matrix, `varprobs` included, a reload agreeing with a copy, and
  leaves the 20 kinds, the suite and tests/cpp as they were.

Recommendation: (i). Cost against (iii): about 500 lines more, almost all tests.

## The rule

A single declared forest restricted to columns C is the fit on `x[, C]`: same trees, same draws under the
same seed and residual scale estimate, with zero counts reported for the other columns. It holds under
every tree prior, DART included: its Dirichlet is laid over C and the other columns report probability 0.

## Constraints

- A fit that states no `vars` on a single declared forest draws exactly what it draws now. `vars` naming
  every column is no restriction and stores none.
- The restriction is a model fact: it is stored on the model object beside the interaction and block
  constraints, so `copy`, a reload and the direct door rebuild it. When the control migration gives every
  model a forest record, it moves there with the rest.
- `vars` resolves as on any forest ([`resolveModerators`](../../R/model.R)): unknown name, empty
  selection and out-of-range index are refused with the same messages.
- `blocks` beside `vars` partitions the allowed columns, as the manual says, on a single forest and on the
  first forest of a multi-forest fit (run: the first forest of two is refused today for not naming the
  columns it may not split on).
- The residual scale estimate still comes from a linear fit on every column: it is the data object's, not
  a forest's.
- On a hazard fit the `period` column the expansion adds is always allowed, named in `vars` or not; `vars`
  restricts the caller's columns. The manual says so.
- No change to the multi-forest path's own masks, to the variance forest, to the flat C header, or to the
  state format (a mask is not state).
- No NEWS entry: `forest()` is new in 1.0-0.
- Out of scope, to TODO: a zero split probability does not exclude a column once the positive ones run out
  (inherited from 0.9-34; `split.probs` alone, no `vars` needed); all-zero `split.probs` fails with a raw
  R error.

## Steps

1. Engine: a column list on [`SamplerOptions`](../../src/bartcore/chain.hpp) for the single mean forest,
   built into the forest's `columnMask` in the single-forest constructor before DART is initialized and
   before [`Chain::installBlockMasks`](../../src/bartcore/chain.hpp), and set on every tree
   ([`Tree::setColumnMask`](../../src/bartcore/tree.hpp)); the same list on
   [`MultinomialForestSpec`](../../src/bartcore/combiner.hpp), installed by
   [`Chain::buildMultinomialForest`](../../src/bartcore/chain.hpp) on every category forest. tests/cpp: a
   restricted single-forest chain never splits outside its list over 500 sweeps, with and without blocks,
   on the constant, linear and monotone leaves; a state that splits outside is infeasible; a multinomial
   chain likewise; an empty list is bitwise the chain without one.
2. DART over the allowed columns: [`DartPrior::initialize`](../../src/bartcore/model.hpp) takes the mask;
   the dimension, the default `rho` and the concentration grid use the number of allowed columns,
   [`DartPrior::update`](../../src/bartcore/model.hpp) draws for allowed columns only, and an excluded
   column's probability is 0. With no mask every expression is what it is now. tests/cpp: an unrestricted
   DART chain is bitwise unchanged; a restricted one reports 0 for excluded columns.
3. Bridge: [`parseModel`](../../src/R_interface_bartcore.cpp) reads the model's resolved columns (integer,
   in range, else an error), [`optionsFromParsed`](../../src/R_interface_bartcore.cpp) hands them to the
   options, [`buildMultinomialSampler`](../../src/R_interface_bartcore.cpp) to the category spec.
4. R, in [`resolveSamplerSpec`](../../R/spec.R): the first forest's `vars` is resolved once; on a model
   with no bases it is attached to the model (omitted when it names every column) and it narrows what
   `blocks` must partition; on a multi-forest model it narrows the first forest's blocks the same way.
   `split.probs` that gives no positive probability to any allowed column is refused by name.
   [`dbartsSampler$setModel`](../../R/dbarts.R) refuses, before anything is stored, a model whose
   restriction differs from the sampler's, absence included (dec-B254).
5. tinytest, a new file. Fails today: split variables over 300 single-sweep runs stay within `vars`, read
   three ways (`getTrees`, `getForestVariableCounts`, `run()`'s `varcount`), beside the same fit without
   `vars`, which must split on the excluded column (the response depends on it). Then: the sub-matrix
   identity at equal `sigest` for gaussian, probit, `interactions`, a caller's `split.probs` and DART;
   allowed 0/1 columns at one cut each, 2000 sweeps; names, indices, a factor column, every column named,
   the three refusals; `dbartsSpec` and the direct door; the families and leaves listed above; a hazard
   fit, where `period` is split on though `vars` does not name it; `blocks` (accepted over the allowed columns, refused when it
   names another, each tree within its group), on a single forest and on the first of two; a variance
   forest with its own `vars`; `copy`, save and reload, `setState` (own state `TRUE`; a foreign state
   refused, sampler unchanged), `installTrees`, `setPredictor`, `setData`, `sampleTreesFromPrior`,
   `growFromRoot`; `setModel` with the sampler's own model edited (accepted, restriction kept in a copy)
   and with a model of another or no restriction, both directions (refused by name, live sampler and
   stored model unchanged).
6. Mutations: the mask not installed; the model attribute not written; DART initialized without the mask;
   `setModel`'s check removed. Each fails the new tests; report the counts.
7. A design note in docs/design/ with its index row (the defect, the rule, why not split probabilities,
   the measurements); [`forest`](../../man/forest.Rd)'s `vars` item (a single forest; DART and
   `split.probs` beside it; `blocks` over the allowed columns on the first forest too; a hazard fit's
   `period` column, always allowed); three TODO entries: the two defects left out above, and that the attribute moves
   into the first forest's record when the control migration gives every model one.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library (a `SamplerOptions` field is a layout change);
  full tinytest; `tests/cpp` builds and passes, clean under ASan and UBSan; the new tinytest file under
  ASan on the R-loaded path.
- On a reference build: the four seeded snapshot files pass unchanged, and the three equivalence compares
  are bitwise against the current baselines, every scenario identical (none states `vars` on a single
  forest), counted per scenario. Nothing is re-recorded.
- Every exact gate in `.github/workflows/exact-gates.yaml` in quick mode, unchanged. No new exact harness:
  the restricted fit is, bit for bit, a fit those gates already cover.
- One script on the base and slice builds digesting seeded draws of the 20 kinds above: equal.
- `lintr::lint_package()`, `air format --check .`, the three tools/ checks,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, and `R CMD check --as-cran` on a tarball from
  a clean copy.
- Not a hot-path change: unrestricted trees keep a null mask, and DART's update gains one branch per
  column per sweep.

## Calls made in planning

- DART beside `vars` on one forest is built (step 2), so the rule has no exception; the alternative was
  refusing the pair by name, as a multi-forest fit refuses DART.
- On a hazard fit the `period` column stays allowed whatever `vars` names: the caller did not supply it and
  cannot be expected to name it, and a discrete-time hazard that cannot vary over periods is not the model
  asked for. The first draft took `vars` literally there.
- The zero-probability defect is not fixed here; it goes to TODO with the `setModel` work of dec-B254, which
  already has to settle a zero split probability on a column in use.
