# bartcore: the merge review

Current at bd91c5e2 (bartcore), 2026-09-06.

This is the case for merging the bartcore branch into main. Sections 1 to 6
are the decision; Appendix A is the tour, what to read and in what order,
for a reader who opens the code. Code is cited by symbol, not by line
number.

## 1. What bartcore replaces

bartcore replaces the classic engine - the one dbarts 0.9-x shipped -
outright. `dbarts/R_C_interface.hpp` and the C++ ABI behind it are deleted;
`inst/include/dbarts/dbarts.h` is the only shipped header, a flat C API.

The one structural idea. The leaf model - the prior over a terminal node's
parameter and its conjugate draw given the observations in that node - is a
compile-time template parameter `L`, because `accumulate` and
`logIntegratedLikelihoodForNode` must inline. A template has no
runtime-uniform handle, hence the `facade.hpp` type-erasure layer between
the C API and the engine. The response family - the likelihood of the
response given the sum of trees, with its own draws after each sweep, sigma
and any latents - is a runtime virtual, chosen once per chain.

Three counts are easy to conflate. The engine enumerates six response
families - gaussian, probit, logistic, aft, ordinal, nbinom - in
`src/bartcore/model.hpp`'s `ResponseFamily`. Everything else called a family
here (multinomial, bcf, heteroscedastic, student, hazard, hurdle)
composes or reduces to those six. And `docs/design/feature-matrix.md` scores
13 rows, counting each composition a user selects as its own model.

BCF's own R verb, `bcf()`/`bartBCF`, lives in bartCause on its `dbarts-1.0`
branch; dbarts carries only the multi-forest engine it is built on.

## 2. Breaking changes for R users

`inst/NEWS.Rd`'s 1.0-0 UPGRADING block is the authoritative list. The breaks
most likely to bite a real script:

- Sampling no longer advances R's random stream, so seeded draws differ from
  0.9-x.
- Saved sampler states and `dbartsData` objects need a version-matched
  rebuild, not a reload.
- `bart2` defaults to `combineChains = TRUE`.
- Unordered factors split on subsets of their levels, and an ordered factor
  becomes a single column split at the midpoints between its consecutive
  declared levels, where 0.9-x expanded both into indicator columns; either
  is posterior-changing for a fit carrying one
  (`docs/plans/column-kind-consolidation.md`, sections 1 and 6).
- A new `missing` argument keeps rows with missing predictors instead of
  dropping them.
- An argument name foreign to the method called is refused by name rather
  than silently discarded, across `predict`, `extract`, `fitted` and
  `residuals`.

## 3. Breaking changes for linked packages

The main branch ships a C++ ABI, `inst/include/dbarts/*.hpp`, with a
C-callable face over it, `R_C_interface.hpp`, whose sampler entries take a
`dbarts::BARTFit*` and whose setters return `void` apart from the predictor
setters' rollback flag, failing through `Rf_error`. Both are gone. `dbarts.h`
is the whole contract, and its head comment is the authority on the three
classes a non-void return can belong to: VALUE, TRANSACTION result,
CAPABILITY STATUS. Below is what a caller of `R_C_interface.hpp` meets.

- Six entries that were `void` on `R_C_interface.hpp` return `int` on
  `dbarts.h`, where every sampler entry is spelled `dbarts_sampler_<name>`:
  `setResponse`, `setOffset`, `setSigma`, `setTestPredictors` (singular
  `setTestPredictor` there), `setTestOffset`, `predict`. The `int` is a
  CAPABILITY STATUS: 1 means the call did its work, which on an ordinary
  gaussian sampler it always does; 0 means this sampler's model has nothing
  for the call to act on, and it was left untouched - `setSigma` on a probit
  sampler, whose residual scale is pinned at 1, or `setTestOffset` on a
  multi-forest sampler. A bad argument still raises. The answer is a fixed
  property of the sampler, so probe once at setup; a caller that ignores
  the return still compiles. `setWeights` is new to the C API - it had no C
  entry point, only `dbarts::BARTFit::setWeights` behind the C++ ABI.
- `dbarts_sampler_getTrees` and `dbarts_sampler_printTrees` take `forest` as
  argument 2; a single-forest caller passes 0. The ABI hash is the backstop
  against a stale call site that a C compiler only warns about.
- `dbarts_sampler_setWeights` answers CAPABILITY STATUS 0 for probit,
  ordinal, aft and nbinom, none of which carries a weight to change, and
  raises on an out-of-support logistic count or gaussian weight.
- `dbarts_predictor_source` is the predictor-input struct for every
  predictor-taking entry, `structSize`-versioned: the caller sets
  `structSize` and may leave `denseCodes` null. Build the dense case with
  `dbarts_dense_predictor_source()`.
- Unlike the deleted `R_C_interface.hpp`, `dbarts.h` neither includes
  `<Rversion.h>` nor defines `USE_FC_LEN_T`; a consumer that relied on that
  pull-in must include it itself.
- `DBARTS_C_API_MAJOR` is 1 and `DBARTS_C_API_MINOR` is 0.
  `DBARTS_C_API_HASH` is recomputed at every ABI change - a signature,
  struct field, enumerator or callback parameter, not a header edit alone -
  so read it from the header at the merge tip rather than from this
  document.

Migration runs in lockstep, once dbarts installs clean
(`docs/plans/capi-shape.md` section 11):

| consumer | mandatory source edits |
|---|---|
| stan4bart, branch `bartcore` | none: it already passes `forest`. It carries `DBARTS_REQUIRE_EXACT_ABI`, so a later hash change forces only a rebuild |
| treatSens, branch `dbarts-1.0` | none: it calls neither `getTrees` nor `printTrees`. Not R-API-only, though - its main branch links the deleted C++ ABI |
| bartCause, branch `dbarts-1.0` | none: R API only, no `src/`, no `dbarts_` symbols |
| bairrtt, no compat branch | none: R API only, with no dbarts linkage in its `src/` |

`TODO`'s `release` item re-verifies all four against the final header.

## 4. What is checked

A green gate proves what its row says and no more.

| gate | what it proves |
|---|---|
| `check-standard` | `R CMD check` clean of errors and warnings, plus NEON kernels checked against scalar on Windows arm64 |
| `cpp-tests` | the C++ component suite green; a seventh `ResponseFamily` enumerator is a compile error |
| `sanitizers` | ASAN and UBSAN over engine and bridge; any finding fails |
| `exact-gates` quick | 21 exact-posterior and move-balance scripts, against closed forms rather than snapshots |
| `exact-gates` cross-host | bcf and multinomial equivalence at tier 1 |
| `equivalence.R` gaussian | 50 scenarios reproduce bitwise on one host |
| `sbc.R` | simulation-based calibration (SBC) over five family arms (gaussian, ordinal, nbinom, Student-t, multinomial) and 30 functionals, Bonferroni-corrected, with nbinom's two dispersion functionals waived as an adjudicated mixing ridge |
| `rchk` | PROTECT balance |
| `valgrind` | leaks and out-of-bounds reads |
| `revdep-smoke` | reverse dependencies install and run |

A cross-host comparison has two tiers. Tier 1, a tight relative-deviation
bound on the draws themselves (`rtol = 1e-8`), is the gate; tier 2, a Welch
z over posterior summaries, is a weaker fallback that cannot gate on its own
(`docs/plans/bcf-cross-host.md`). Within one host, reproducibility is
bitwise across every SIMD dispatch path.

The rewritten engine matches the shipped one where the priors match: the
equivalence harness's statistical mode ran released 0.9-34 against this
branch over 16 scenarios, 4 at high precision, with zero unexplained
disagreements, every large z tracing to a documented change
(`docs/plans/review-2026-08-24/anchor-main.md`, sections 4 "Explained
differences" and 5 "Unexplained disagreements").

## 5. What is not checked

Five workflows - `equivalence`, `sbc`, `rchk`, `valgrind`, `revdep-smoke` -
are `schedule` plus `workflow_dispatch`, and GitHub binds both triggers to
the default branch, which does not carry them. On bartcore each fires only
from a push that touches its own file, and each has run once that way.
`equivalence` and `revdep-smoke` passed. `rchk` reported eight unprotected
uses of a data frame's names attribute in the model-matrix code and one
multi-argument slot read in the multinomial bridge, false positives on the
running program that are protected anyway; the same image now reports zero
findings. `valgrind` found a 48-byte leak on the C API's test-missingness
refusal, a C++ object destroyed by a longjmp, since fixed; the full suite is
clean under it. `sbc` flags the nbinom dispersion functionals on the
identifiability ridge `docs/plans/sbc-family-tiers.md` adjudicates as mixing
rather than miscalibration, so that arm waives those two by name and fails
on any other. Merging to main registers the schedules (`TODO`'s `release`
item).

Things that could be wrong and would not be caught:

- No equivalence scenario and no SBC coverage reaches a multi-forest
  amplitude sampler - one whose forests enter the fit through per-forest
  amplitude scalars, as BCF's `a*mu + b_z*tau` does - under a latent family,
  probit or logistic; the BCF calibration evidence is gaussian-only
  (`docs/plans/review-2026-08-24/calibration-sbc.md`).
- aft and heteroscedastic are uncovered at ensemble scale, and both carry
  sampling code that reduces to no covered family. hazard and hurdle are not
  scored directly either; their draws are checked to reproduce bitwise the
  draws a covered family makes on the corresponding data, so they inherit
  that family's calibration.
- Warm start and grow-from-root are unrefused and untested at two forests.
- The cross-host tier-2 bar is weak by construction: it tolerates a shift of
  about 1.4 posterior standard deviations; a 20 percent node-prior widening
  passes tier 2 and fails tier 1. Its fix, independent per-scenario seeds
  rather than one chain's autocorrelated draws, waits until after the
  release candidate (`TODO`'s `equivalence-harness-seeds-axis`).
- The C++ mutation record,
  `docs/plans/review-2026-08-24/mutation-B-findings.md`, which planted 80
  deliberate engine mutations and recorded which ones the C++ component
  suite missed, has not been re-run against this tip, so those gaps are not
  confirmed closed.
- Nothing tests that `setState` itself honours the containment verdict -
  that a restored state's splits stay inside the columns the model allows
  (`sampler.hpp`'s `allValid = columnMaskOk`).
- Three `benchmarks/R` harnesses run in no workflow, one of them calling
  itself a gate, so drift in what they measure goes undetected because
  nothing re-runs them.
- `setForestBasis(k, ~var)` evaluates the formula in its own environment
  with no data attached, so a column living only in a data frame is not
  found.
- A per-forest weight is not part of saved state, and an active-row mask is
  mirrored nowhere, so two states can compare `identical()` while their fits
  diverge (`docs/design/bcf.md`, `docs/design/bart-as-a-component.md`).

## 6. Decided, open, and more expensive after the merge

Four scope questions are decided: `updateScale` is refused on every
multi-forest sampler, whatever its family, by a guard that reads the forest
count; real-valued nbinom dispersion and weighted binary are scheduled
after 1.0-0; formal heredity is the first work after 1.0-0. `TODO` carries
the last three as `negbin-real-dispersion`, `weighted-binary` and
`interaction-constraints`.

One question is open: whether to declare the release candidate (`TODO`'s
`rc-gate` item).

No shipped surface still needs changing before the release. Four surfaces
would be expensive to change after it, and each is in its final form:
`gp()` is calibrated at 25 trees, inside the range its man page recommends
(a GP leaf earns its keep at tens of trees, not hundreds), in four
configurations including the one where a tree holds both GP and
constant-fallback leaves (`docs/plans/sbc-calibration.md`, Tier C); the
pointwise log-likelihood on a BCF fit is pinned against a hand computation
in all three families BCF supports, gaussian, probit and logistic; the
heteroscedastic swap under `updateScale = TRUE` is refused; and a sampled GP
lengthscale would be an additive state block, not a format break.

## Appendix A. The tour: what to read, in order

This is the reading order for a reviewer who opens the code and the
documents after sections 1 to 6. It is ordered by what a linked package can
be broken by, and each document is placed at the stop where its subject
comes up. Word counts are given where they are known, so you can budget by
them. In the four design documents, only the sections named below state
the current design - about 4,100 of their 16,000 words; the rest can be
skipped.

### 1. Orientation

Open: `docs/architecture.md` - the current state, not a history; prefer it
to any paraphrase where the two overlap. It is the one document to read
whole before any code.

Then: the code walk starts at the surface a linked package compiles
against.

### 2. The C API

Open: `inst/include/dbarts/dbarts.h`, `src/C_interface.cpp` - the head
comment's contract list, then the X-macro entry table - and
`docs/plans/capi-shape.md` sections 0 and 13.

Judge: whether every non-void entry says which of the three return classes
it is, and whether a discarded capability 0 is a failure mode you accept:
it leaves the sampler unchanged and the run conditioned on what it held
before, quieter than `R_C_interface.hpp`'s `Rf_error` longjmp.

Then: behind that surface is the engine those entries call into.

### 3. The engine

Open: `facade.hpp`, `sampler.hpp`, `chain.hpp` - `SamplerBase` and its pure
virtuals, `SamplerFacade`, the `create*Sampler` factories; `Sampler`, `run`,
`predictColumns` fanning out over `std::thread` workers via
`fanOutPredictSlabs`; `Chain`, `setActiveRows`, `columnMaskStateFeasible`.
Prefer `docs/architecture.md` on RNG and threading.

Judge: the exhaustive `ResponseFamily` switch, which carries no `default:`
arm anywhere, and that state restore is semantic, not bitwise.

### 4. Multiple forests

Open: the mutation-legality table first, then the code that enforces it,
then the one weight that code does not save.

- `docs/design/bart-as-a-component.md`, sections 2 "Which mutations are
  legal between sweeps" and 3 "What engine state does not carry, and who
  reinstalls it", about 850 words: which mutations a multi-forest sampler
  admits, and the two state gaps, the per-forest weight and the active-row
  mask.
- `docs/design/multiplier-combiner.md`, the preamble's first paragraph,
  then "The model", "The amplitude layout", "The reparameterization", "The
  amplitude conditional", "bcf as the K = 2 instance", "Surfaces" and "What
  this family does not do", about 1,490 words: what the basis-and-amplitude
  family is, and where BCF sits in it.
- `combiner.hpp`: `ForestCombiner`, `AmplitudeForestCombiner` (saved-state
  key `"glue"`, after the per-forest amplitude scalars that glue the
  forests into one fit), `MultinomialForestCombiner`. BCF's `a*mu + b_z*tau`
  is the two-forest instance of the amplitude-and-basis family
  `docs/design/multiplier-combiner.md` sets out.
- `docs/design/bcf.md`, the preamble's model equation and "The multiplier
  snap and the per-forest weight (2026-08-10)", about 355 words: why a row
  can carry an exact-zero weight in one forest, and why that weight is not
  saved state.

Judge: which mutations the combiner refuses and why.

### 5. The R bridge

Open: `src/R_interface_bartcore.cpp` - `bartcore_create`, `_run`, the
setters, `_storeState`, `_setState`, `_installForests`, `_predict`,
`_predictPerForest`, `_getTrees`, then the shared guards
`refusedAmplitudeFamilyReason`, `refuseMultiForestMutation`,
`refuseUndefinedTestFits`, `refusePinnedSigmaChange`, `refuseNonBinaryMask`.
`tests/cpp/test_facade.cpp` is the facade's conformance test, one row per
`SamplerBase` virtual driven through the base.

Judge: `refusePinnedSigmaChange`'s own comment, the source's clearest
statement of why a guard is keyed on family rather than an internal flag.

### 6. Tree moves and data

Open: `docs/design/empty-leaf-veto.md`, "Where the constant is read", then
"Is vetoed-vs-vetoed reachable? Yes; the veto is a RANK (2026-08-18)",
"What counts as empty: the weight law (2026-08-12)" and "Which weights the
predicate sees", about 1,410 words: the member-empty versus weight-empty
ranking. Then `moves.hpp`, `tree.hpp`, `scan.hpp`, `grow.hpp`, `data.hpp`:
`metropolisJumpForTree`; `Tree`, `columnMaskSubtreeIsValid`;
`scanOrdinalCuts`, `growTreeFromRoot`; `ColumnStore`, `ScopedCutGrid`,
`ColumnKind` and the derived `kindSplitsBySubset`.

Judge: change-move detailed balance; the ranked empty-leaf veto in
`Tree::leafVetoRank` and `resolveVetoRank`, where a member-empty leaf
vetoes absolutely and a weight-empty leaf is only penalized
(`docs/architecture.md`'s "Tree moves"); whether the semantic kind axis and
the mechanical `splitsBySubset` axis stay separate (only grid construction,
ingestion validation and reporting may read the kind); and the doubled
entry layout `scanOrdinalCuts` uses for a node holding missing members.

Then: with the mechanisms read, the grid that scores them can be judged.

### 7. The capability grid

Open: `docs/design/feature-matrix.md` - the one deep read: the per-model
capability grid, and a Gaps section collecting every missing cell as a
candidate work item.

Judge: the cell values. Its cites are machine-checked; its cell values are
judgments.

### 8. Build support

Open: the build support files - `configure`, `tools`, `src/misc`,
`src/external` - are skim-only; the one thing worth a look is `simd.c`'s
`cpuid`, which requests subleaf 0 so that AVX2 is not misdetected as AVX as
it is in 0.9-x.

### 9. Reference, not reading

Open as the questions arise, not in order:

- `docs/design/INDEX.md`, `docs/plans/INDEX.md` - complete manifests,
  refused and closed items included.
- `docs/plans/release-candidate-review.md` - the pre-release review's master
  log, newest first.
- root `TODO` - an alphabetical backlog, some items scheduled after 1.0-0.
  Its `release` item is the one ordered procedure.
