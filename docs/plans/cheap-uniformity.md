# cheap-uniformity

agent: S0 sonnet, S1 opus (data.hpp + bridge + R mutation are one contract),
  S2 sonnet, S3/S4 opus (Branch A, engine) or sonnet (Branch B, bridge).
  Serialized: each slice lands before the next starts.
rng: neutral, EVERY slice. No slice touches a proposal, a suffstat, a cut
  builder, or an RNG consumer. Piece 1/2 hand the engine a dense block
  byte-identical to today's R-densified one; piece 3 touches only the flat
  replay, which no draw reaches (see Context). Gate: the full trio bitwise,
  no snapshot regenerated. A draw move is falsifier F8 - stop and report.
window: builds on test-data-parity (7d781a8, the resident test container),
  typed-ingestion, and csc-code-validation. Precedes the dbarts.h reshape
  arc, which mirrors the surface in "Source-shaped mutation surface" below.
budget: S0 ~150 lines R+tests; S1 ~450 across data.hpp,
  R_interface_bartcore.cpp, R/bartcore.R, R/mixedMatrix.R, tests/cpp/
  test_data.cpp, a new tinytest file; S2 ~250; S3+S4 (Branch A) ~800;
  S3B (Branch B) ~250.

## Goal

Three dense-only asymmetries close: a sparse-valued argument mutates a
sparse-backed design without the caller densifying, a single column of a
container-backed `x.test` can be replaced, and a sparse test set predicts
without a full n x p materialization. Along the way three live defects die:
a transposed `dgCMatrix` silently reinterpreted as a mutation argument, a
missing-value predicate that is wrong in both directions on a container,
and a container whose declared sparse reference level means one thing to
`predict()` and another to `setTestPredictor()`. No `dbarts.h` change.

## Binding decisions inherited (do not reopen)

1. Piece 1 is scoped to CONTAINER / CSC-backed targets only. A plain-matrix
   `data@x` keeps today's R-side `as.double` densification and today's dense
   `.Call`: `rawColumnForRequantize (data.hpp)` makes the caller-supplied
   matrix the sole re-quantize source for a dense build, so `data@x` cannot
   stop being a dense matrix and a bridge scratch there is purely additive.
   Door 1B (per-column streaming inside `runPredictorTransaction`) stays
   shut.
2. Ingestion semantics gate on the STORE's `columnTypes[j]`, never on the
   source's declared reference. This is what the engine already computes
   (`resolveCscCategoricalReferences` skips every non-categorical store
   column; `mutateCscColumnFromDense` keys its implicit on `types[j]`).
3. The dimension narrowing is keyed on "`dim(x)` non-NULL and disagrees",
   NOT on a class allowlist. `dim(x)` NULL keeps today's length-only check
   verbatim, and every `Matrix` class that is not `dgCMatrix` keeps today's
   `as.double` densification verbatim.
4. Per-observation mutation of a CSC-backed column stays refused
   (`bartcoreSamplerSetPredictor`, `bartcore_updatePredictorPerObservation`).
   `forceUpdate = "partial"` is refused for any sparse-VALUED argument.
5. Piece 3 is an OPEN FORK (see "Piece 3 fork"). S0-S2 are branch-
   independent and land first.

## Context (seams, all read at 39451b1)

- `PredictorSource` (data.hpp) is the one borrowed view; `isDenseBlock()` is
  the predicate the mutation entrances refuse on. The mutation virtuals
  (`SamplerBase::setPredictor`/`updatePredictor`, facade.hpp) ALREADY take
  the view, so pieces 1-2 need no facade signature change.
- `bartcore_setPredictor` requires a 2-element `R_DimSymbol` and reads
  `REAL`; `bartcore_updatePredictor` checks length only. Neither sees an S4
  argument today because R densifies first.
- `bartcoreSamplerSetPredictor` gates its dimension checks on `is.matrix(x)`,
  which is FALSE for every `Matrix` class - the transposed-dgCMatrix hole.
  Its whole-matrix sparse-source branch already densifies to a dense `.Call`
  and maintains `data@x` separately through `installPredictorColumns`.
- `bartcoreSamplerSetTestPredictor` refuses ANY column update when
  `@x.test` inherits `dbartsMixedMatrix`. The dense flavor's per-column
  update is already copy-modify plus a whole-object install, so piece 2
  needs no engine or bridge writer. There is no per-column test writer and
  none is being added (test-data-parity Open Q2 resolved whole-object).
- Flat replay: `partitionFlatIndices` and the four `*Below` helpers (tree.hpp)
  are reached only by `Chain::addFlatPredictions`, `bartcore_getTrees`, and
  tests/cpp - never by the draw loop, which descends on quantized codes
  through `findBottomNodeForRow`. Row fits accumulate only from a row's own
  leaf, so row chunking is bitwise-exact.
- THREE raw reads, not one: `partitionFlatIndices` (~1617), and the
  leaf-covariate reads inside `addFlatLinearPredictionsBelow` (~1735) and
  `addFlatFunctionPredictionsBelow` (~1819). The latter two index by the
  STORE column index `columns[j]`, which under a MAPPED source is not the
  block position.
- SIX `addFlatPredictions` call sites in chain.hpp: the four `predictFrom*`
  entries plus `predictVarianceFromSavedSample` (~741) and `predictVariance`
  (~769). `tests/cpp/test_model.cpp` is the only exerciser of the variance
  replay.
- Design docs: docs/design/sparse-columns.md (the standing sparse contract),
  docs/design/data-store.md (read first - data-adjacent work).

## Defects fixed, and where

| defect | lands | failing-first test |
|---|---|---|
| transposed `dgCMatrix` mutation argument silently reinterpreted | S0 | new refusal case, `test-data-mixed-mutation.R` |
| `anyNA` predicate wrong in BOTH directions on a container | S0 | two cases (spurious refusal, missed NA) |
| container reference level means different things to `predict()` and `setTestPredictor()` | S2 | pre-fix disagreement captured, then both refuse |

The `anyNA` defect, precisely: `checkMissingPolicy(data, anyNA(x), ...)` is
called with the user's raw argument. Base `anyNA` on a list is TRUE iff some
element is a length-one atomic NA, so a container with exactly ONE sparse
ordinal column - `sparseReference = NA_integer_`, the ordinal sentinel
`wrapSparseTestMatrix` writes - reports TRUE and is refused under
`missing = "error"` for missingness it does not have. That is live at tip on
the landed `setTestPredictor` container path. In the other direction, a
container whose `sparse@x` holds a real NA reports FALSE. `recursive = TRUE`
fixes NEITHER: it descends correctly into the S4 `dgCMatrix` (verified,
R 4.6.1 / Matrix 1.7.5) but it also descends into `sparseReference`, so it
returns TRUE for EVERY container carrying an ordinal sparse column - a
strictly worse false positive. The metadata fields must be excluded by hand.

The reference-semantics defect, precisely: `as.matrix.dbartsMixedMatrix`
takes a sparse column's implicit rows to be `sparseReference[k]` whenever it
is non-NA, gated only on `is.na`. The engine takes them to be the reference
code only when the STORE types the column categorical, and 0 otherwise. The
two disagree for a container declaring a non-NA reference against a
store-ordinal column, which is reachable through documented calls: train
with an ORDERED factor column (typed ordinal, but carrying a `factor.levels`
entry), then supply that column as a `sparseFactor` in a test frame -
`mapFactorColumnsToTrainingLevels` remaps it rather than refusing, since its
numeric-training guard keys on a NULL level table. `predict()` then
densifies with the reference implicit while `setTestPredictor()` stays
resident with implicit 0.

## Constraints

- Refuse, never coerce or clamp. Reuse existing error texts where they fit;
  `parseCscMatrix`'s dimension text ("number of rows of 'x' must equal
  length of 'y'") is creation-flavored and must not surface from a mutation
  call - give the helper a caller-supplied message or pre-check `Dim`.
- Code-range validation bounds the MATERIALIZED block, not the stored
  entries: an implicit row's value is bounded like any other
  (`validateColumnValues`).
- A mutation argument's declared K is metadata only; the store's category
  count is creation-pinned and its `numCuts` is the bound.
- No factor-level alignment on the TRAIN mutation path. Alignment is the
  test side's job, at the `validateXTest` funnel only.
- Compute the new `data@x` / `@x.test` BEFORE the `.Call`; assign only on
  acceptance. A rejected transaction must leave the R-side source
  byte-untouched.
- Out of scope: per-observation mutation of a CSC-backed column; an engine
  `setTestColumn` writer; door 1B; any `inst/include/dbarts/dbarts.h` edit.

## S0. R validation hygiene (R only, no new capability)

1. `sourceAnyNA(x)` in R/utility.R, dispatching:
   `dbartsMixedMatrix` -> `any(vapply(x$dense, anyNA, NA))` over the dense
   list (NULL dense yields FALSE) `|| (!is.null(x$sparse) && anyNA(x$sparse@x))`;
   a `sparseMatrix` -> `anyNA(x@x)` where the class carries values, FALSE
   for a pattern matrix; anything else -> `anyNA(x)`. Never `anyNA` on the
   container list, and never `recursive = TRUE`. Replace the three
   `checkMissingPolicy` arguments (`setPredictor`, `setTestPredictor`,
   `setTestPredictorAndOffset`).
2. Dimension checks keyed on `dim()`, in both `bartcoreSamplerSetPredictor`
   branches and `bartcoreSamplerSetTestPredictor`'s column branch: when
   `dim(x)` is non-NULL it must agree with the target (existing texts
   "dimension of x must be equal to N" / "number of columns of new x does
   not match length of columns to replace"); when it is NULL, today's length
   check runs verbatim.
3. Nothing else changes: a plain numeric vector, a `sparseVector`, and every
   non-`dgCMatrix` `Matrix` class keep today's `as.double` path. This is
   load-bearing - `bench-sampler.R`'s `setPredictor-accept-n1000-t75` and
   `setPredictor-reject-n1000-t75` arms both pass a plain length-n vector.

Gates: new refusal/acceptance cases green and failing on the unpatched R;
full tinytest; `air format --check .`; trio bitwise (formality - no C moves).
Falsifier: if any dim check refuses a shape that today accepts AND today
produces correct results, the narrowing is not confined to the mis-shaped
case - re-scope (F2).

## S1. Piece 1: sparse-valued mutation onto a container/CSC-backed target

1. data.hpp gains a free `materializePredictorSource(const PredictorSource&,
   const ColumnType* storeTypes, size_t rowBegin, size_t rowEnd, double* out)`
   beside `PredictorSource`: fill each column with its implicit value, then
   scatter the stored entries over it. Implicit = `referenceCodeOf(j)` iff
   `storeTypes[j]` is categorical, else 0 (binding decision 2). Null
   `storeTypes` means all-ordinal. The row range exists so Branch B of piece
   3 reuses the same helper; S1 always passes the whole range.
2. Bridge: `bartcore_setPredictor` gains an `unwindProtect` wrapper - it has
   none today, and `validateColumnValues`'s `Rf_error` would longjmp past a
   `std::vector<double>` destructor. Both mutation entrances branch on
   `Rf_isReal` FIRST (today's path, byte-identical), then `dgCMatrix`, then
   `dbartsMixedMatrix`; the sparse branches build a view with the existing
   `parseCscMatrix` / `parseMixedContainerBlock` / `mapColumnSources` /
   `resolveCscCategoricalReferences` helpers (`parseTestContainer` is the
   template), refuse a container declaring a non-NA reference against a
   store-ordinal column, materialize, run the existing per-column
   `validateColumnValues` loop, and call the existing dense entry.
3. R5 `bartcoreSamplerSetPredictor`: when `predictorSourceIsSparse(data@x)`,
   pass the sparse argument straight to `.Call` instead of
   `matrix(as.double(x), ...)`. A plain-matrix target is untouched.
4. `installPredictorColumns` accepts sparse `values` in its CONTAINER branch
   only: move `values <- as.double(values)` to per-column so an untouched
   column is never materialized; a dense-backed target column materializes
   its one column (O(n)); a sparse-backed target column whose implicit
   matches the argument's splices `@i`/`@x` directly (O(nnz)) but must
   CANONICALIZE first - drop stored entries equal to the implicit, the
   `is.na(v) | v != implicit` rule `replaceSparseColumn` documents and
   `mutateCscColumnFromDense` mirrors - or `data@x`'s pattern diverges from
   the store's. The `is.matrix(x)` early-return branch is untouched.
5. Memory, stated honestly: a whole-matrix replacement onto a container
   target still peaks at O(n x p) - the scratch moves from an R matrix to a
   C++ vector, it does not shrink. The per-column entrance peaks at
   O(n x k). What drops is R-side container maintenance: O(n) per touched
   column instead of the whole dense argument. The slice buys uniformity and
   validation, not a lower whole-matrix peak.

Gates: tests/cpp `test_data.cpp` materializer cases {ordinal CSC,
categorical CSC with refCode != 0, mixed dense+CSC, all-implicit column,
stored NaN, mapped source whose dense column is not at its own index};
tests/cpp from clean plus the ASAN leg; new `test-mutate-sparse-valued.R`
(the oracle below); full tinytest; trio bitwise; `air format --check .`;
rchk on the next scheduled run. data.hpp is a header: `R CMD INSTALL
--preclean`, and delete the `benchmarks/kernels` binaries (no header
dependency tracking there).
Oracle: for {bare dgCMatrix design, mixed container design} x {whole matrix,
single column, multi-column} x {dgCMatrix, sparseVector, dbartsMixedMatrix}
argument, a sampler mutated with the sparse argument is `expect_identical`
in fits to a twin mutated with the dense equivalent, and `as.matrix(data@x)`
agrees. Plus: a rejected transactional update leaves `data@x` byte-
untouched; a fresh sampler built from the mutated `data@x` matches the
mutated sampler; a malformed reference-vs-store-type container is refused;
an NA in a sparse argument trips `missing = "error"`;
`setCutPoints`/`setState` after a sparse-valued mutation still work
(`rawPredictorMatrix` returns NULL and the store re-quantizes from its
retained slices).
Falsifier F1: if the sparse argument is not BITWISE equal to its dense
equivalent on any accepted shape, the materializer disagrees with the
engine's implicit rule. Stop; do not relax to a tolerance.

## S2. Piece 2: per-column x.test update, and the reference-semantics pin

1. In `bartcoreSamplerSetTestPredictor` the container refusal becomes a
   container branch: `installPredictorColumns(data@x.test, NULL, column,
   values)`, then the EXISTING whole-object `.Call`. The container's
   per-column storage decision is preserved; the existing optimistic-install-
   plus-`tryCatch` rollback already covers a bridge refusal. R only.
2. The A6 pin, one rule at two funnels. R: in `validateXTest`, after
   `alignContainerFactorLevels` and before the resident/densify decision,
   refuse a container whose sparse column declares a non-NA reference for a
   column `attr(x.train, "varTypes")` does not type categorical - the same
   attribute the funnel already trusts for `factor.levels`. Bridge:
   `refuseCscReferenceAgainstStore(types, parsed)` called immediately after
   every `parseTestContainer` (creation parse, `bartcore_setTestPredictor`,
   `bartcore_setTestPredictorAndOffset`), beside
   `validateTestContainerAgainstStore`; it takes the STORE's types at
   mutation and the parsed types at creation, so no attribute is trusted
   where the store is available. The reverse mismatch (an NA reference
   against a store-categorical column) already errors inside
   `resolveCscCategoricalReferences` and stays as is.
3. Cost, stated honestly: a per-column container update rebuilds the whole
   test store, O(nnz + n_test x p_dense) - identical in kind to today's
   dense per-column update, which also rebuilds wholesale. No speed win.

Gates: extend the test-side tinytest file - a per-column update of a
container-backed `x.test` equals a whole-object install of the same spliced
container AND the dense equivalent, for a sparse-backed and a dense-backed
column of the same container, by index and by name, with a bridge refusal
rolling `@x.test`/`@offset.test` back; a NEW test builds the ordered-factor
training design plus the `sparseFactor` test column, asserts the pre-fix
`predict()` vs `setTestPredictor()` disagreement on the unpatched build, and
asserts both refuse after. Full tinytest; trio bitwise; `air format --check .`;
rchk next run.
Falsifier F4: if a per-column container update's fits differ from a
whole-object install of the same spliced container, the R-side splice is not
building the container the engine would build.

## Piece 3 fork (OPEN - VD decides before S3 starts)

S3 alone delivers nothing user-visible, so there is no partial landing: the
choice is `{Branch A}` or `{Branch B}`, decided against the target shape
n_test = 1e5, p = 1e4, density 0.01. Branch B is bitwise-exact and cheap but
BOUNDS densification rather than removing it; Branch A removes it and is the
only branch that answers the directive.

Both branches share: `bartcore_predict` and `bartcore_getTrees` accept a
`dgCMatrix`/`dbartsMixedMatrix` `xTestExpr`/`newdataExpr` through
`parseTestContainer` + `refuseCscReferenceAgainstStore` +
`validateTestContainerAgainstStore`; a dense real matrix keeps
`validatePredictorMatrix` verbatim; the two R densification lines
(`dbartsSampler$predict` and `$getTrees(newdata=)`) are deleted; the
heteroscedastic route still returns `list(mean, variance)`; `n_test == 0`
reaches the same error a dense zero-row matrix does; a sparse-backed LEAF
COVARIATE column is refused with `installTestContainer`'s existing text.

### Branch B: row-chunked densification at the bridge (fallback)

`bartcore_predict` materializes the parsed view in row chunks (4096 rows)
through S1's `materializePredictorSource` row range, calls the existing
dense `predict` per chunk, and scatters each chunk's rows into the right row
range of every [row x sample x chain] slab of the result; `predictVariance`
takes the same treatment. `bartcore_getTrees`'s `newdata` replay
materializes WHOLE (one n_new x p scratch, no worse than today's R
`as.matrix`) - chunking it would need a count-accumulating variant of
`countFlatObservationsBelow` for no memory benefit on a diagnostic path.
Peak input memory becomes O(chunkRows x p); the O(n_test x p) time stays.
ZERO engine change, no new virtual, no `--preclean` beyond S1's.
Bitwise-exactness is the row-independence argument in Context, and the gate
is that every existing assertion in `test-predict-sparse.R` keeps passing
VERBATIM plus a chunk-boundary case (n_test not a multiple of the chunk).

### Branch A: engine-resident sparse predict (S3 then S4)

S3 - templating, DENSE ONLY, zero behavior change.
- `partitionFlatIndices` and the four `*Below` helpers become templates on a
  `Columns` type exposing `column(j)` -> a reader with `at(i)`. ALL THREE raw
  reads route through it: the hoist in `partitionFlatIndices` AND the
  leaf-covariate reads in `addFlatLinearPredictionsBelow` and
  `addFlatFunctionPredictionsBelow`, which index by the STORE column index
  and would otherwise read the wrong column under a mapped source - a silent
  wrong prediction, not a crash. Leaving a dangling `const double* x`
  parameter after the signature change is the same defect.
- The existing `(const double* x, size_t numRows, ...)` signatures stay as
  non-template wrappers over a `DenseColumns{x, numRows}`, so chain.hpp,
  `bartcore_getTrees`, `test_state.cpp`, and `test_tree.cpp` move by zero
  lines. The dense instantiation must compile to today's code.
- chain.hpp: ALL SIX `addFlatPredictions` call sites gain the templated
  form - the four `predictFrom*` entries AND
  `predictVarianceFromSavedSample` and `predictVariance`. The
  heteroscedastic sparse test in S4 is unbuildable without the last two.
- Gates: tests/cpp from clean INCLUDING `test_model.cpp` (the only
  exerciser of the variance replay) and a new `test_tree.cpp` case with a
  linear-leaf and a function-leaf model over a MAPPED source whose leaf
  covariate is dense-backed but NOT at index j; every existing
  `test-predict-sparse.R` assertion verbatim; the ASAN leg; the codegen A/B
  below. `--preclean` (tree.hpp and chain.hpp are headers) and delete the
  `benchmarks/kernels` binaries. Land it alone - this is the risk slice.

S4 - the sparse route.
- data.hpp gains `SparseRawColumn`: `SparseColumnData`'s shape (bits,
  wordRanks, O(1) `at(i)`) holding raw doubles and a `double implicitValue`
  instead of codes and a `zeroCode`. No cut grid, no quantization, so it
  touches no draw-law surface. BUILT PER CALL, NEVER CACHED - `predict`
  reads its argument, not the resident store, so there is no staleness risk
  and no implementer should introduce one. The rank bitmap is correct for
  duplicate values (rank indexes the pattern, not the values).
- A NEW internal virtual `predict(const PredictorSource&, size_t numRows,
  double* out)` and its `predictVariance` twin on `SamplerBase`; the frozen
  dense entries become non-virtual inlines building a dense view, exactly
  the `setTestPredictors` precedent. FACADE VIRTUAL CHANGE: `--preclean` is
  mandatory (a stale object bus-errors) and the `benchmarks/kernels`
  binaries must be deleted.
- Gates: the new tinytest bitwise oracle (resident-sparse `predict` and
  `getTrees(newdata=)` identical to the dense twin; pooled categorical
  K > 63, exercising the out-of-range clamp in `partitionFlatIndices`;
  heteroscedastic `list(mean, variance)`; the sparse leaf-covariate
  refusal); trio bitwise; tests/cpp plain and ASAN from clean; F7 below.
- Falsifier F5: if any existing `test-predict-sparse.R` assertion stops
  holding, the refactor leaked into the dense path - abort the slice.
- Falsifier F7 (PRE-REGISTERED): on n_test = 1e5, p = 1e4, density 0.01,
  resident-sparse predict must be at least 2x faster than the densify path
  AND at least 5x smaller in peak input memory. If not, the arc's most
  expensive piece has no payoff - land Branch B instead.

Codegen gate for S3 (replaces F6, which `bench-sampler.R` cannot fire - it
has no predict arm and adding one would force a baseline re-record). Use the
house same-machine A/B protocol (docs/design/memory-wall-frontier.md
sections 11-12): a throwaway harness timing dense `predict` on a fixed
sampler, base and rolled build BACK-TO-BACK, min of 9 reps, 3 self-gating
resumable rounds, 1-min loadavg stamped per round, valid at loadavg < 6
judged by CROSS-ROUND ratio spread rather than absolute quiet. Record the
per-round ratios in this plan's landing note. A dense regression beyond ~2%
that reproduces across rounds means the reader indirection did not compile
away - revert to Branch B.

## Verification (every slice)

- `R CMD INSTALL .`; `--preclean` for S1 (data.hpp) and for S3/S4 of Branch
  A (headers, and a facade virtual in S4). Delete `benchmarks/kernels`
  binaries after any header edit; tests/cpp tracks headers via `-MMD`.
- `cd tests/cpp && make clean && make && ./test_bartcore` - all pass. ASAN
  leg (`ASAN_OPTIONS=detect_container_overflow=0`) for every slice that
  makes new engine code reachable: S1 (the data.hpp materializer) and S3/S4
  under Branch A. S0, S2, and Branch B add no engine code.
- Full `tinytest::test_package("dbarts")` from a preclean install; new tests
  ADD, no snapshot regenerated.
- `benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-c8f661a.rds`
  -> 27/27 "identical draws (same RNG stream)".
  `bcf-equivalence-99205ee.rds` -> 5 scenarios x 6 channels bitwise.
  `multinomial-equivalence-ec2a3d0.rds` -> 3 scenarios x 5 channels bitwise.
  No max-|z| line anywhere. THE TRIO IS NECESSARY, NOT SUFFICIENT: its
  `sparse` scenario is a dgCMatrix TRAINING design with a DENSE `x.test` and
  no mutation, so it exercises no surface this arc changes. It proves only
  that the dense path did not move. The NEW tests are the real oracle -
  `test-mutate-sparse-valued.R` for S1, the extended test-side file plus the
  ordered-factor divergence test for S2, the extended
  `test-predict-sparse.R` for piece 3.
- `air format --check .` on any slice touching R/.
- rchk on the next scheduled run for S1 and S2 (both move the bridge).
- Speed: no slice is a hot-layer change, so `bench-sampler.R` compare
  against `bench-sampler-ab1dc52.csv` is a formality, not a gate. S3 of
  Branch A is gated by the A/B protocol above instead.

Stop conditions per docs/plans/README.md: a step fails twice, the diff
exceeds 1.5x budget, or a needed change is out of scope - report and stop.

## Edge cases the tests must name

Sparse argument onto a dense-backed column of a mixed container (legal, the
dense kernels run); a stored entry >= the store's K (refused with the
existing categorical text); an all-implicit sparse column (materializes
constant, and the ordinal cut refresh refuses degeneracy exactly as a dense
constant column does); explicit stored zeros in the argument (materialize to
the implicit; the R splice MUST canonicalize); NA (Matrix stores it as a
stored NaN, implicit rows are never NA); `sparseVector` with
`length(column) > 1` (refused - no unambiguous column layout);
`numColumns == 0` (already errors); a BCF or multinomial sampler
(`denseCreationPredictorSource` drops sparse storage at creation, so a
sparse-valued mutation materializes and runs the dense path); a data-handle
view (`refuseMutationOnView` fires first).

## Source-shaped mutation surface (for the dbarts.h reshape arc to mirror)

Nothing here changes dbarts.h. Recorded so the follow-on copies a signed-off
shape. The C mirror of `PredictorSource` is a POD with the same ten fields
(`numRows`, `numColumns`, `denseValues`, the CSC triple, `columnSources`,
`columnTypes`, `categoryCounts`, `referenceCodes`), consumed by
`setPredictor`, `updatePredictor`, `setTestPredictors`, and `predict`.
Semantics this arc pins:

1. `columnSources == NULL` means column j is dense column j of
   `denseValues`; a nonnegative entry names a dense column; a negative one
   names CSC column `~v`.
2. A CSC column's implicit rows read 0 UNLESS THE STORE types column j
   categorical, in which case they read that column's `referenceCodes`
   entry. A reference declared against an ordinal column is malformed and
   refused - it is NOT an alternative implicit value. (The engine ignores
   `referenceCodeOf(j)` outright for an ordinal column; the R-side
   `as.matrix` rule that gates only on `is.na` is the divergence S2 closes.)
3. The source is BORROWED for the call; the callee retains nothing from it.
4. A mis-shaped source is REFUSED, never reinterpreted.
5. Rule 2's malformed-reference refusal applies at EVERY funnel - creation,
   mutation, test ingestion, and predict - not only at mutation.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S0: a mis-shaped sparse predictor argument (for example a transposed
  `dgCMatrix`) is now refused rather than silently reinterpreted
  column-major.
- S0: `missing = "error"` no longer refuses a sparse container that holds no
  missing values, and no longer accepts one that does.
- S1: `setPredictor`/`updatePredictor` accept `dgCMatrix`, `sparseVector`,
  and `dbartsMixedMatrix` values against a sparse or mixed design without
  the caller densifying.
- S2: a single column of a sparse or mixed test matrix can be replaced.
- S2: a sparse test column declaring a reference level for a training column
  that is not categorical is refused; previously `predict` and
  `setTestPredictor` read such a column differently.
- Piece 3: prediction on a sparse test set no longer materializes the full
  dense matrix (Branch A) / bounds its peak memory (Branch B).

## Open items

- The piece-3 fork is unresolved and belongs to VD. Branch A is the only
  branch that answers "without densifying" and is the only one whose risk
  needs the A/B protocol; Branch B is bitwise-exact, ~250 lines, and bounds
  rather than removes the materialization. Deciding evidence: whether the
  target shape (n_test = 1e5, p = 1e4, density 0.01) is a shape a consumer
  actually runs. S0-S2 do not wait on it.
- Design artifacts (memo, blind critique, synthesis) are durable at
  `.claude/cheap-uniformity-design/`. The critique's verdict was STANDS
  WITH AMENDMENTS; A1-A10 are adopted here, with two corrections recorded
  in the synthesis (A5's mechanism and A6's reachability route).
