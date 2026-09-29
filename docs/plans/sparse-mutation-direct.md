# sparse-mutation-direct: a sparse replacement goes sparse to sparse

Status: PLANNED 2026-09-28

agent: opus (both slices; slice B's man, NEWS and design-doc edits may go to sonnet)
rng: neutral, bitwise (the codes, cut grids, stored pattern and missing flags a sparse replacement leaves are the
ones the dense path leaves today, so no draw moves)
window: before the RC (TODO sparse-mutation-direct)
budget: ~720 lines in two serial slices. A, engine: ~430 (data.hpp ~170, sampler.hpp ~50, facade.hpp ~6
comments, tests/cpp ~200). B, bridge and surface: ~290 (bridge ~60, C_interface.cpp ~6 comments, R ~10 comments,
man ~8, NEWS ~8, design docs ~30, tinytest ~170).

## Goal

A sparse replacement through `setPredictor` (whole matrix or named columns) reaches the engine as each sparse
column's stored rows and values, and the engine rebuilds that column's pattern and codes from them. No dense
rows-by-columns block is built for a sparse column, so the call's peak memory follows the replacement's
nonzeros. Every accepted replacement leaves the store bitwise as it does today. The five places that say a
replaced sparse column "densifies its storage permanently" say what is true.

## Context

- Ruling: dec-B135 (maintainer, 2026-09-28), replacing dec-A35. The engine stays R-agnostic and the bridge
  converts (dec-A50, dec-B85): the new engine entry takes plain rows, values and a count, and every R-shaped
  decision (reference metadata, level-code messages) stays in the bridge. The shipped header carries no
  predictor mutation entry (dec-B86 withholds the pair until a consumer asks), and this ruling does not reopen
  it. NEWS covers only changes against 0.9-x (dec-B128); sparse predictors are new in 1.0-0.
- Today's pipeline. The bridge parses a `dgCMatrix` or `dbartsMixedMatrix` argument into a borrowed view
  ([`parseMutationSource`](../../src/R_interface_bartcore.cpp)), then
  [`materializeMutationSource`](../../src/R_interface_bartcore.cpp) expands it into an n x columns block
  ([src/R_interface_bartcore.cpp:1026-1049](https://github.com/vdorie/dbarts/blob/91e3db8640a8c36cc9f2e208c9fd6845001ec838/src/R_interface_bartcore.cpp#L1026-L1049)),
  called from [`bartcore_setPredictor`](../../src/R_interface_bartcore.cpp)
  ([src/R_interface_bartcore.cpp:5563-5567](https://github.com/vdorie/dbarts/blob/91e3db8640a8c36cc9f2e208c9fd6845001ec838/src/R_interface_bartcore.cpp#L5563-L5567))
  and [`bartcore_updatePredictor`](../../src/R_interface_bartcore.cpp)
  ([src/R_interface_bartcore.cpp:5620-5627](https://github.com/vdorie/dbarts/blob/91e3db8640a8c36cc9f2e208c9fd6845001ec838/src/R_interface_bartcore.cpp#L5620-L5627)).
  The engine's [`ColumnStore::setPredictors`](../../src/bartcore/data.hpp) and
  [`ColumnStore::setColumns`](../../src/bartcore/data.hpp) send each CSC-backed column to
  [`ColumnStore::mutateCscColumnFromDense`](../../src/bartcore/data.hpp), which scans all n rows for the ones that
  differ from the column's implicit value and hands the result to the shared install-rebuild-quantize tail.
- The transaction: [`runPredictorTransaction`](../../src/bartcore/sampler.hpp) over
  [`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp) and [`SubsetUpdate`](../../src/bartcore/sampler.hpp), both
  holding a `const double*` block today; the view overloads
  [`SamplerBase::setPredictor`](../../src/bartcore/facade.hpp) and
  [`SamplerBase::updatePredictor`](../../src/bartcore/facade.hpp) accept only a dense block. Rollback mechanics:
  [Predictor mutation transaction](../design/data-store.md#predictor-mutation-transaction).
- Creation already builds a CSC column's cuts from its stored entries:
  [`ColumnStore::quantileGridForCscColumn`](../../src/bartcore/data.hpp) and
  [`ColumnStore::fillCutsUniformlyCsc`](../../src/bartcore/data.hpp), reached from
  [`ColumnStore::buildCutsForColumn`](../../src/bartcore/data.hpp).
- The R side is already sparse to sparse: [`installPredictorColumns`](../../R/mixedMatrix.R) splices each column
  through [`sparseEntriesForColumn`](../../R/mixedMatrix.R) and [`replaceSparseColumns`](../../R/mixedMatrix.R) in
  O(nnz), materializing a column only when its implicit value differs from the target's.
- The feature record: [Extension (i): sparse-column in-place mutation (landed)](../design/sparse-columns.md#extension-i-sparse-column-in-place-mutation-landed).

## Paths that can receive a sparse or mixed source

Established by reading the code at 91e3db86, by runs against a private library built from it, and by tests/cpp.

1. Whole matrix, R: [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) with `column` missing on a design whose
   `data@x` is a `dgCMatrix` or a container with a sparse column passes a sparse argument to the bridge unchanged
   and splices it into `data@x` on acceptance. The bridge materializes it (above). Densifies today.
2. Named columns, R: the same method with `column` given, whether the columns are dense-backed, sparse-backed or
   both, and whether the argument is a `dgCMatrix` or a container. Same bridge step. Densifies today (n x the
   number of named columns).
3. Per-observation (`forceUpdate = "partial"`): refused on a sparse-backed column by name, and
   [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) coerces the argument with `as.double`. A dense-backed column
   only, one column; no change.
4. A plain-matrix design given a sparse argument: `as.double` in R before the bridge. The store keeps no raw for
   such a design and re-quantizes from `data@x`, which must stay a dense matrix; out of scope.
5. `setData`: refuses a sparse or mixed design and argument
   ([`bartcore_setData`](../../src/R_interface_bartcore.cpp)); no change.
6. Test side: [`bartcore_setTestPredictor`](../../src/R_interface_bartcore.cpp) and
   `bartcore_setTestPredictorAndOffset` parse a container with
   [`parseTestContainer`](../../src/R_interface_bartcore.cpp) and rebuild the test store from the view through
   [`installTestContainer`](../../src/R_interface_bartcore.cpp); a bare `dgCMatrix` becomes a container in
   `validateXTest`, and the per-column R path splices into the container. None of them densifies. Measured at
   n = 1e5, p = 200, 1% density: a sparse `setTestPredictor` peaks at 259-278 MB against 245-274 MB without it.
   No change.
7. The C API: inst/include/dbarts/dbarts.h has no predictor mutation entry. The one entry
   that takes a `dbarts_predictor_source`, `dbarts_sampler_predict`, replays CSC columns in place. No change.
8. Transactional rollback: [`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp) snapshots every CSC column's
   source descriptor, rank storage and owned buffers when the store was built from CSC;
   [`SubsetUpdate`](../../src/bartcore/sampler.hpp) takes a
   [`ColumnStore::CscColumnRollback`](../../src/bartcore/data.hpp) per touched CSC column. Both restore the state
   the new entry writes, so they are reused unchanged.
9. Cut refresh: [`ColumnStore::refreshCutsForColumn`](../../src/bartcore/data.hpp) reads a dense column, and so do
   the precheck [`ColumnStore::cutsWouldRemainValid`](../../src/bartcore/data.hpp) and the degenerate test
   [`ColumnStore::valuesAreDegenerate`](../../src/bartcore/data.hpp). Each needs a CSC form (Design).
10. Categorical sparse columns (`sparseFactor`): the store's implicit value is the reference code fixed at
    creation. The source's implicit value is the container's declared reference, resolved by
    [`resolveCscCategoricalReferences`](../../src/R_interface_bartcore.cpp), or 0 for a bare `dgCMatrix`. The two
    may differ, and today both are accepted: the materialized block carries the source's value and the dense
    scan re-sparsifies against the store's. A declared reference on an ordinal store column is refused by
    [`refuseCscReferenceAgainstStore`](../../src/R_interface_bartcore.cpp) first.
11. Storage tiers: a CSC-backed column is rank-bitmap or densified-codes, fixed at build and never flipped by a
    mutation (the [`ColumnStore::mutateCscColumnFromDense`](../../src/bartcore/data.hpp) contract). Both tiers
    keep their retained slice in `ownedCscRows` and `ownedCscValues`, so both take the new entry.
12. Missing values: a stored `NA` is a stored NaN, which never equals the implicit value, so it stays stored and
    sets `hasMissing` in [`ColumnStore::quantizeCscColumnInto`](../../src/bartcore/data.hpp). The R missing-policy
    check reads the argument itself (`sourceAnyNA`) and is unaffected.
13. Storage kind differs between source and store: a mixed argument may hold a column dense that the store keeps
    CSC, or sparse that the store keeps dense (a dense-backed column of a mixed design). Both are accepted today,
    through the block (the twin tests in ["a sparse argument onto a DENSE-backed column"](../../inst/tinytest/test-mutate-sparse-valued.R)).

## Design

### The engine entry

`ColumnStore::mutateCscColumnFromCsc(j, rows, values, nnz, sourceImplicit, updateCuts)`, beside
[`ColumnStore::mutateCscColumnFromDense`](../../src/bartcore/data.hpp). It builds the column's new stored pattern
under the same rule, keep an entry iff it differs from the store's implicit value (so NaN is kept):

- `sourceImplicit` equals the store's implicit value (every ordinal and ordered-factor column, where both are
  0, and a categorical column whose source declares the store's reference): one pass over the nnz entries,
  dropping explicit implicit-valued entries (stored zeros, -0.0, reference-coded entries). O(nnz).
- They differ (a categorical column whose source reads its absent rows as another level): one merge pass over n
  rows, since every absent row of the source is now a stored entry of the store. O(n) for that column only.

Both branches then call a new private tail, `installCscColumn`, which is today's tail of
`mutateCscColumnFromDense` moved out: compare with the old pattern, move the vectors into `ownedCscRows[j]` and
`ownedCscValues[j]`, repoint the slice, refresh the cut grid if asked (below), rebuild the rank bitmap when the
pattern changed, and quantize. `mutateCscColumnFromDense` keeps its dense-column cut refresh before extraction
and calls the same tail.

### Routing by storage kind

The two strategies hold the caller's `PredictorSource` instead of a `const double*`. The store gains a
source-aware `setPredictors(const PredictorSource&, updateCuts)` and `setColumns(const PredictorSource&,
columns, count, updateCuts)`, keeping the dense spellings as wrappers (tests/cpp calls them), and one per-column
router, `mutateColumnFromSource(j, source, k, updateCuts, scratch)`:

| store column | source column | path |
|---|---|---|
| CSC-backed | CSC | `mutateCscColumnFromCsc`, no dense column at all |
| CSC-backed | dense | `mutateCscColumnFromDense` on the source's own column pointer (no copy) |
| dense | dense | today's dense path, unchanged (`setColumnJournaled` on the subset path) |
| dense | CSC | materialize that one column into an n-double scratch reused across columns, then the dense path |

The implicit value of a source CSC column is decided by the store's type, as
[`materializePredictorSource`](../../src/bartcore/data.hpp) decides it today: the source's reference code on a
categorical store column, 0 otherwise. The existing view overloads of `setPredictor` and `updatePredictor` on
[`SamplerBase`](../../src/bartcore/facade.hpp) and [`Sampler`](../../src/bartcore/sampler.hpp) widen to accept
mapped and CSC views: no new virtual, no vtable change. A plain dense block is the identity map, so every dense
caller takes the dense row of the table exactly as today.

### Cut refresh and feasibility

- Refresh runs after the install and reads the new slice with creation's builders,
  [`ColumnStore::quantileGridForCscColumn`](../../src/bartcore/data.hpp) and
  [`ColumnStore::fillCutsUniformlyCsc`](../../src/bartcore/data.hpp), through `refreshCutsForCscColumn(j)`, the
  CSC sibling of [`ColumnStore::refreshCutsForColumn`](../../src/bartcore/data.hpp). Its uniform-mode refusal
  needs `cscColumnIsDegenerate(j)`, the sibling of
  [`ColumnStore::valuesAreDegenerate`](../../src/bartcore/data.hpp): no two logical values differ, the implicit
  value counting when any row is absent. Factor columns keep no grid to refresh, as today.
- The precheck in [`runPredictorTransaction`](../../src/bartcore/sampler.hpp) reads a CSC source column through
  `cutsWouldRemainValidCsc(j, values, nnz, sourceImplicit)`: level validity of every stored value and of the
  implicit value when any row is absent, for a factor column; the quantile grid's induced count, for a numeric
  one. The grid collector of `quantileGridForCscColumn` is refactored to take (values, nnz, implicit present),
  so build, refresh and precheck share one collector over the entries.
- Why the grids are bitwise the dense path's: the logical values are the same set, one implicit value plus the
  stored entries, and the collectors already match at creation.

### Bridge

`materializeMutationSource` becomes `prepareMutationSource`: the reference refusal and resolution it does today,
and no block (`ParsedMutationSource::block` is deleted). Both entrances hand `parsed.view` to the view overloads,
or a `densePredictorSource` over `REAL(x)` for a plain matrix. The per-column level-code check becomes
`validateSourceColumnValues(store, j, view, k)`, calling
[`validateColumnValues`](../../src/R_interface_bartcore.cpp) on a dense column, or on the stored values plus the
implicit value when the column has any absent row, so the messages and the refusal order are unchanged. The R
methods need no change beyond comments.

### dbarts.h

No change: no entry, struct, enumerator or documented contract moves, and `DBARTS_C_API_HASH` is untouched.

### Peak memory, measured

Design: n = 1e5, p = 200, 1% density (200,000 stored entries), gaussian, 1 chain, 20 trees, one whole-matrix
`sampler$setPredictor(x2)` with default arguments. Whole-process peak footprint from `/usr/bin/time -l`, three
runs per arm, macOS arm64. The "after" arm is a scratch prototype of this design, not the landed code.

| arm | peak footprint | increment over no call |
|---|---|---|
| no mutation (both builds) | 249-277 MB | - |
| before: dense block | 414-421 MB | ~160 MB, the 1e5 x 200 double block |
| after: sparse to sparse | 259-266 MB | ~10-15 MB, within run-to-run noise |

Replacing 100 named columns peaks at 332 MB before and 277-287 MB after. The call itself drops from ~25 ms to
~11 ms. The prototype passed the full tinytest suite (9630 of 9630) and tests/cpp, and a scratch store-level
check matched `mutateCscColumnFromDense` bitwise (cuts, codes, slice, rank bitmap, missing flag) over both
tiers, both cut modes, with and without a refresh, stored zeros, -0.0, NaN, an all-absent column, and
categorical columns with a matching and a mismatched source reference.

## Doc corrections

All five claims are false today: the engine stores only the entries that differ from the implicit value and
keeps the column's layout, and the R splice canonicalizes the same way. No NEWS entry is added for the memory
change (dec-B128); the NEWS text describing the feature is corrected in place.

- man/dbartsSampler-class.Rd, the `column` argument ("Note that replacing a sparse-backed column densifies its
  storage permanently: ... it stores n entries from that point on."). Replace with: "A replaced sparse-backed
  column stays sparse: it stores the replacement's entries that differ from its implicit value (zero, or a
  sparse factor's reference level), in the storage layout chosen at creation, and a sparse replacement is read
  as those entries rather than expanded, so the call's memory follows the replacement's nonzeros."
- man/dbarts.Rd, the `formula` argument ("replacing a sparse-backed column densifies its storage
  permanently."). Replace with: "a replaced sparse-backed column stays sparse, storing only the replacement's
  nonzeros, and a sparse replacement is never expanded to a dense matrix."
- inst/NEWS.Rd, the sparse-predictor item of NEW FEATURES:
  - "(and the flat C API's \code{setPredictor} / \code{updatePredictor})": delete. The header has no such entry
    (dec-B86), and man/dbartsSampler-class.Rd's Mutation cost section already says so.
  - "a replaced column's storage densifies permanently, since every row then differs from the implicit zero."
    Replace with: "a sparse replacement is read as its stored entries, never expanded, and a replaced column
    stores only its nonzeros."
  - "a replaced sparse column's storage densifies permanently, as for an all-sparse design." Replace with: "a
    replaced sparse column stays sparse, as for an all-sparse design."
- R/bartcore.R, [`bartcoreSamplerSetPredictor`](../../R/bartcore.R): "Replacing a sparse column densifies its
  storage - every row now differs from the implicit value." Replace with: "A replaced sparse column stays
  sparse: the engine and the splice both keep only the entries that differ from its implicit value." The comment
  above the `predictorSourceIsSparse` test ("it materializes there") becomes "the bridge hands its sparse columns
  to the engine as stored entries, under the store's own implicit rule".

Other text this change makes stale, corrected in slice B:

- docs/design/sparse-columns.md, extension (i): its first bullet ("hands the engine a DENSE column even for
  sparse storage") and the 2b supersession bullet's last sentence ("densifies its storage permanently"). The
  section is a landing record, so append a "SUPERSEDED by sparse-mutation-direct" bullet saying both are no longer
  true rather than rewriting them.
- docs/design/data-store.md: "Only the MUTATION entrances still need the dense values as one block" (The
  borrowed view's two value channels), and a sentence under Predictor mutation transaction on the per-column
  routing.
- Comments: [`Sampler::setPredictor`](../../src/bartcore/sampler.hpp) and
  [`SamplerBase::setPredictor`](../../src/bartcore/facade.hpp) ("only a dense block is consumable"),
  [`ColumnStore::mutateCscColumnFromDense`](../../src/bartcore/data.hpp) ("the mutation surface hands a dense
  column even for sparse storage"), [`ColumnStore::setPredictors`](../../src/bartcore/data.hpp),
  [`ParsedMutationSource`](../../src/R_interface_bartcore.cpp), the two comments in
  [`TranslatedSource`](../../src/C_interface.cpp)'s neighborhood that name a mutation entrance the C API does not
  have, [`sparseEntriesForColumn`](../../R/mixedMatrix.R) and
  [`replaceSparseColumn`](../../R/mixedMatrix.R) (name the sibling beside `mutateCscColumnFromDense` as the rule
  they mirror), and the header of ["a sparse-valued mutation argument"](../../inst/tinytest/test-mutate-sparse-valued.R).

## Agent-made calls

The ruling settles sparse to sparse and a second engine entry. These are not settled by it; each is a
recommendation VD can overturn before implementation.

1. The entry takes a sixth argument, `sourceImplicit`, beside (j, rows, values, nnz, updateCuts). A bare
   `dgCMatrix` onto a categorical column whose reference is not the first level, or a container declaring
   another reference, reads its absent rows as a different level than the store, and that is accepted today.
   Alternatives: refuse a mismatched reference (a behavior change, and R's splice already handles it), or have
   the bridge materialize such a column (the engine keeps five arguments; the bridge grows an O(n) path).
2. No new virtual on [`SamplerBase`](../../src/bartcore/facade.hpp): the existing view overloads widen to mapped
   and CSC views. The alternative, a separate sparse overload, adds a vtable slot and a stale-object hazard for
   no caller the widened one cannot serve.
3. A dense store column fed a CSC source column is materialized one column at a time into a reused n-double
   scratch, then takes the dense path. The alternative, a CSC form of the journaled dense path, is more code for
   the rarest shape; peak stays O(n), not O(n x columns).
4. A mixed argument's dense part is still assembled by `parseMixedContainerBlock` as n x (its dense columns).
   Those columns are dense in the argument already; the ruling is about sparse columns.
5. A plain-matrix design given a sparse argument keeps R's `as.double` (path 4): the store re-quantizes from
   `data@x`, which must stay dense.
6. The CSC cut refresh runs after the install and reads the new slice with creation's builders; the precheck
   reads the not-yet-installed entries. Both share one refactored collector rather than two copies.
7. Canonical pattern: an entry equal to the store's implicit value is dropped, -0.0 included on an ordinal
   column, NaN kept, exactly the dense rule and [`sparseEntriesForColumn`](../../R/mixedMatrix.R).
8. The rank tier still never flips. A rank column replaced by a mostly nonzero column stays rank, with nnz up to
   n: existing behavior, out of scope.
9. No automated memory gate. tinytest cannot bound a C++ peak reliably; the landing note records the
   measurement above rerun on the landed code, and a discrimination mutation proves the sparse entry is the one
   reached (Verification).
10. Two serial slices: A (engine and tests/cpp) changes nothing a user sees, since the bridge still hands dense
    blocks; B switches the bridge and carries the docs. One slice of ~720 lines is the alternative.
11. docs/design/bart-as-a-component.md section 6 still names `dbarts_sampler_setPredictor` and
    `dbarts_sampler_updatePredictor` as the flat surface's predictor channel, stale since dec-B86 and the same
    false claim as the NEWS parenthetical. Fixed in slice B (two lines) unless VD wants it kept separate.

## Constraints

- RNG class neutral, bitwise: any drift in a code, cut, stored entry or missing flag is a defect of this plan,
  not a snapshot to regenerate.
- inst/include/dbarts/dbarts.h and [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) unchanged. No
  change to the facade's virtual set; `--preclean` on every engine commit regardless.
- docs/plans/INDEX.md gains this plan's row when the plan is committed (the freshness guard fails without it).
- Engine code takes no R type and no R message; reference resolution and level-code messages stay in the
  bridge.
- Out of scope: per-observation replacement of a sparse-backed column; `setData` on a sparse design; the rank
  tier flip; a plain-matrix design's R-side densification; any C mutation entry; the test-side paths, which
  already keep sparse storage.

## Steps

Slice A, engine:

1. Move the tail of [`ColumnStore::mutateCscColumnFromDense`](../../src/bartcore/data.hpp) into `installCscColumn`
   and add `mutateCscColumnFromCsc` (Design). tests/cpp passes unchanged.
2. Refactor the CSC grid collector; add `refreshCutsForCscColumn`, `cscColumnIsDegenerate`,
   `cutsWouldRemainValidCsc`. Creation's grids unchanged (tests/cpp's CSC build checks).
3. Add the source-aware `setPredictors`, `setColumns` and `mutateColumnFromSource`; the dense spellings become
   wrappers.
4. Switch [`WholeMatrixUpdate`](../../src/bartcore/sampler.hpp), [`SubsetUpdate`](../../src/bartcore/sampler.hpp)
   and the precheck in [`runPredictorTransaction`](../../src/bartcore/sampler.hpp) to the view; update the
   view-overload comments in sampler.hpp and facade.hpp.
5. tests/cpp, in test_model.cpp beside `testSparseMutation`:
   - store level: `mutateCscColumnFromCsc` against `mutateCscColumnFromDense` on the materialized column,
     comparing `cutPoints` bytes, `numCuts`, every `codeAt`, the retained slice's rows and value bytes, the rank
     bitmap, word ranks, `nzCodes` and `zeroCode`, and `hasMissing`; over both tiers, uniform and quantile cuts,
     `updateCuts` off and on, stored zeros, -0.0, NaN, an all-absent column, and categorical columns on both
     tiers with a matching and a mismatched `sourceImplicit`;
   - `cutsWouldRemainValidCsc` equals `cutsWouldRemainValid` on the materialized column, passing and failing
     cases for quantile counts and level codes;
   - sampler level: `setPredictor` and `updatePredictor` with an all-CSC view and a mixed view (a CSC source on
     a dense store column, a dense source on a CSC store column) against the dense block: codes, cuts, then
     sweeps and test predictions bitwise;
   - rollback: a three-column transactional update whose middle column is refused, once by the precheck
     (`invalidCutPoints`, store byte-identical, nothing written) and once by revalidation (`rolledBack`, codes,
     cuts, rank storage, owned slices and tree fits restored byte for byte), whole and subset.

Slice B, bridge and surface (stacked on A):

6. Replace `materializeMutationSource` with `prepareMutationSource`, add `validateSourceColumnValues`, pass the
   view in both entrances; delete `ParsedMutationSource::block`.
7. The doc and comment corrections above; the NEWS parse check.
8. tinytest, extending inst/tinytest/test-mutate-sparse-valued.R, every arm a sparse argument against its dense
   twin with `expect_identical`:
   - the existing twins, plus the stored cut grid (`attr(sampler$state, "cutPoints")` after `storeState()`),
     `predict` on a held-out set, and further sweeps;
   - `updateCutPoints = TRUE`, uniform and `useQuantiles`;
   - `forceUpdate = FALSE`, accepted and rolled back;
   - categorical: a container with the store's reference, and a bare `dgCMatrix` onto a `sparseFactor` column
     whose reference is "s2";
   - subset path: one column, several, by name, spanning dense- and sparse-backed store columns;
   - rollback with a refused middle column, whole and subset: an out-of-table level code (error; `data@x` and
     the next sweeps identical to an unmutated twin), a quantile precheck failure (error, same checks), and an
     all-absent middle column under `forceUpdate = FALSE` (`FALSE` returned, same checks);
   - the implicit value checked against the store's level count: a container whose declared reference is past
     the store's K is refused with the categorical message.

## Verification

- Each slice: `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e
  'tinytest::test_package("dbarts")'`: zero failures; `cd tests/cpp && make && ./test_bartcore`: "all tests
  passed".
- Neutrality on the reference build (`--configure-args=--enable-reference-build`): the four seeded-drift
  snapshot files and the three bitwise equivalence compares CI's cpp-tests job runs, each "identical draws" with
  no "max |z|".
- Sanitizers: tests/cpp under ASAN and UBSan, and slice B's tinytest file R-loaded under ASAN, per
  [Gate hygiene](README.md#gate-hygiene).
- Discrimination, each run and reverted with a `touch`: keeping an implicit-valued entry in the one-pass branch
  fails the stored-zeros twins and the store-level pattern checks; skipping the cut refresh in
  `installCscColumn` fails the `updateCutPoints = TRUE` arms; routing a CSC source back through the scratch
  column passes every twin, the two being bitwise by design, so a store-level test routes CSC onto CSC and
  requires the router's scratch vector to come back with capacity 0, which that mutation fails.
- Memory, recorded in the landing note: the Peak memory design rerun on the landed library, before and after.
- `air format --check .`; `lintr::lint()` on R/bartcore.R and R/mixedMatrix.R; `Rscript
  tools/check-rc-codoc.R .`, `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each
  on its own exit status; the NEWS parse check; `R CMD check --as-cran` on a tarball from a clean copy (slice B
  touches R/ and man/).
