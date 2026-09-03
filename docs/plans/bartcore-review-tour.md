# bartcore: the merge review

Current at 849f08ea (bartcore), 2026-09-03.

This is the case for merging the bartcore branch into main. Sections 1 to 6
are the decision; the appendices are for a reader who opens the code. Code
is cited by symbol, not by line number.

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
here (multinomial, bcf, grouped, heteroscedastic, student, hazard, hurdle)
composes or reduces to those six. And `docs/design/feature-matrix.md` scores
13 rows, counting each composition a user selects as its own model.

BCF's own R verb, `bcf()`/`bartBCF`, lives in bartCause on its `dbarts-1.0`
branch. It never shipped in dbarts and is not being removed from it; dbarts
carries the multi-forest engine it is built on.

## 2. Breaking changes for R users

`inst/NEWS.Rd`'s 1.0-0 UPGRADING block is the authoritative list. The breaks
most likely to bite a real script:

- Sampling no longer advances R's random stream, so seeded draws differ from
  0.9-x.
- Saved sampler states and `dbartsData` objects need a version-matched
  rebuild, not a reload.
- `bart2` and `rbart_vi` default to `combineChains = TRUE`.
- Unordered factors split on subsets of their levels, and an ordered factor
  becomes a single ordinal column by default.
- A new `missing` argument keeps rows with missing predictors instead of
  dropping them.
- An argument name foreign to the method called is refused by name rather
  than silently discarded, across `predict`, `extract`, `fitted` and
  `residuals`.
- An ordinal fit's R-visible spelling is `thresholds`, matching the engine
  and C API's `ordinalThresholds`.

An ordered-factor predictor is its own column kind, not a bare ordinal
predictor: it splits on the midpoints between consecutive declared levels
rather than on `n.cuts` uniform cuts over the observed codes, which is
posterior-changing for any fit carrying one
(`docs/plans/column-kind-consolidation.md`, sections 1 and 6).

## 3. Breaking changes for linked packages

`dbarts.h` is the whole contract, and its head comment is the authority on
the three classes a non-void return can belong to: VALUE, TRANSACTION
result, CAPABILITY STATUS. Below is only what changed.

- Seven setters went `void` -> `int`: `setResponse`, `setOffset`,
  `setWeights`, `setSigma`, `setTestPredictors`, `setTestOffset`, `predict`.
  Source-compatible everywhere; the answer is a fixed property of the
  sampler, so probe once at setup.
- `forest` is argument 2 on `dbarts_sampler_getTrees` and
  `dbarts_sampler_printTrees`, the last two entries where it followed
  `useLiveTrees`. A pointer where an integer is now expected fails to
  compile in C++ but may only warn in C; the ABI hash is the backstop.
- `dbarts_sampler_setWeights` answers CAPABILITY STATUS 0 for probit,
  ordinal, aft and nbinom, none of which carries a weight to change, and
  raises on an out-of-support logistic count or gaussian weight.
- `dbarts_predictor_source` gains `denseCodes` and `numDenseCodeColumns`,
  appended after the fields a 1.0-0 consumer compiles against, so an
  existing caller's layout is unchanged and a null `denseCodes` stays
  valid.
- `dbarts.h` no longer defines `USE_FC_LEN_T` nor includes `<Rversion.h>`; a
  consumer that relied on that pull-in must include it itself.
- `DBARTS_C_API_MAJOR 1` and `DBARTS_C_API_MINOR 0` do not move.
  `DBARTS_C_API_HASH` is recomputed at every ABI change - a signature,
  struct field, enumerator or callback parameter, not a header edit alone -
  so read it from the header at the merge tip.

Migration runs in lockstep, once dbarts installs clean
(`docs/plans/capi-shape.md` s11):

| consumer | mandatory source edits |
|---|---|
| stan4bart, branch `bartcore` | two: the `getTrees` and `printTrees` call sites in `src/init.cpp`. It carries `DBARTS_REQUIRE_EXACT_ABI`, so a later hash change forces only a rebuild |
| treatSens, branch `dbarts-1.0` | none: it calls no reordered entry. Not R-API-only, though - its main branch links the deleted C++ ABI |
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
| `equivalence.R` gaussian | 51 scenarios reproduce bitwise on one host |
| `sbc.R` | simulation-based calibration (SBC) over five family arms and 30 functionals, Bonferroni-corrected |
| `rchk` | PROTECT balance |
| `valgrind` | leaks and out-of-bounds reads |
| `revdep-smoke` | reverse dependencies install and run |

A cross-host comparison has two tiers. Tier 1, a tight relative-deviation
bound on the draws themselves (`rtol = 1e-8`), is the gate; tier 2, a Welch
z over posterior summaries, is a weaker fallback that cannot gate on its own
(`docs/plans/bcf-cross-host.md`). Within one host, reproducibility is
bitwise across every SIMD dispatch path.

The rewritten engine has also been checked against the shipped one: the
equivalence harness's statistical mode ran released 0.9-34 against this
branch over 20 scenarios, 10 at high precision, with zero unexplained
disagreements, every large z tracing to a documented change
(`docs/plans/release-candidate-review.md`, the Calibration lane paragraph of
the second whole-branch review).

## 5. What is not checked

Five workflows are `schedule` plus `workflow_dispatch` only - `equivalence`,
`sbc`, `rchk`, `valgrind`, `revdep-smoke` - and have never run, by any
trigger, on any branch: GitHub binds both triggers to the default branch,
which does not carry these files. Every clean rchk, valgrind, equivalence
and SBC claim above rests on a manual local run, recorded in
`docs/plans/release-candidate-review.md`. Merging to main registers them
(`TODO`'s `release` item).

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
  about 1.4 posterior standard deviations, and a probe confirmed it passes a
  20 percent node-prior widening that tier 1 fails. Its fix, independent
  per-scenario seeds rather than one chain's autocorrelated draws, waits
  until after the release candidate (`TODO`'s
  `equivalence-harness-seeds-axis`).
- The C++ mutation record,
  `docs/plans/review-2026-08-24/mutation-B-findings.md`, has not been re-run
  against this tip, so the gaps it measured are not confirmed closed.
- Nothing tests that `setState` itself honours the containment verdict -
  that a restored state's splits stay inside the columns the model allows
  (`sampler.hpp`'s `allValid = columnMaskOk`).
- Five `benchmarks/R` harnesses run in no workflow; three of them call
  themselves gates and have one recorded run each, and the ratios
  `grouped-mixing.R` measures now already disagree with its own header's
  figures, undetected because nothing re-runs it.
- `setForestBasis(k, ~var)` evaluates the formula in `environment(basis)`
  with no data attached, so a column living only in a data frame is not
  found.
- A per-forest weight is not part of saved state, and an active-row mask is
  mirrored nowhere, so two states can compare `identical()` while their fits
  diverge (`docs/design/bcf.md`, `docs/design/bart-as-a-component.md`).

## 6. Decided, open, and more expensive after the merge

Four scope questions are settled and need no ruling: `updateScale` stays
refused on every multi-forest family, keyed on the sampler's forest count
rather than its family; real-valued nbinom dispersion and weighted binary are
scheduled after 1.0-0; and formal heredity is the first work after 1.0-0.
`TODO` carries the last three as `negbin-real-dispersion`, `weighted-binary`
and `interaction-constraints`.

One question is open: whether to declare the release candidate (`TODO`'s
`rc-gate` item).

Three pieces of shipped surface are cheap to change now and become a second
breaking change after release. Leaving each as it stands is a legitimate
ruling, made deliberately:

- `gp()` is documented in `man/bart2.Rd`, `man/dbarts.Rd`,
  `man/dbartsPriors.Rd` and `man/xbart.Rd`. Its calibration evidence is a
  25-tree run, not an ensemble-scale one.
- The heteroscedastic scale leaf is calibrated once at creation and never
  rescaled when `sigmaScale` moves, so a response or offset swap under
  `updateScale = TRUE` leaves the prior on the old scale (`TODO`'s
  `variance-forest-mutation-routing`).
- Sampling a GP lengthscale rather than fixing it needs per-sample channels
  in the state and flatten formats: a saved-state format change, breaking
  state compatibility a second time (`TODO`'s `gp-followups`).

## Appendix A. Reading the code

Two renames make a grep of the old tree mislead. `BCFForestCombiner` is now
`AmplitudeForestCombiner` in `combiner.hpp`, saved-state key `"glue"` for
`"bcf"`, after the per-forest amplitude scalars that glue the forests into
one fit; the old names find nothing, which does not mean BCF was removed.
And `bartcore_createMultinomial` and its `Counts` variant are retired,
multinomial creating through `bartcore_create` like every other family.

The walk is ordered by what a linked package can be broken by.

ABI - `inst/include/dbarts/dbarts.h`, `src/C_interface.cpp`. The head
comment's contract list, then the X-macro entry table. Judge whether every
non-void entry says which of the three return classes it is, and whether a
discarded capability 0 is a failure mode you accept: it leaves the sampler
unchanged and the run conditioned on what it held before, quieter than the
old longjmp. `docs/plans/capi-shape.md` s0, s13.

Engine - `facade.hpp`, `sampler.hpp`, `chain.hpp`. `SamplerBase` and its
pure virtuals, `SamplerFacade`, the `create*Sampler` factories; `Sampler`,
`run`, `predictColumns` fanning out over `std::thread` workers via
`fanOutPredictSlabs`; `Chain`, `setActiveRows`, `columnMaskStateFeasible`.
Judge the exhaustive `ResponseFamily` switch, which carries no `default:`
arm anywhere, and that state restore is semantic, not bitwise. Prefer
`docs/architecture.md` on RNG and threading.

Multi-forest - `combiner.hpp`: `ForestCombiner`, `AmplitudeForestCombiner`,
`MultinomialForestCombiner`. BCF's `a*mu + b_z*tau` is the two-forest
instance of the amplitude-and-basis family
`docs/design/multiplier-combiner.md` sets out. Judge which mutations the
combiner refuses and why; the mutation-legality table in
`docs/design/bart-as-a-component.md` comes first.

Bridge - `src/R_interface_bartcore.cpp`: `bartcore_create`, `_run`, the
setters, `_storeState`, `_setState`, `_installForests`, `_predict`,
`_predictPerForest`, `_getTrees`, then the shared guards
`refusedAmplitudeFamilyReason`, `refuseMultiForestMutation`,
`refuseUndefinedTestFits`, `refusePinnedSigmaChange`, `refuseNonBinaryMask`.
Judge `refusePinnedSigmaChange`'s own comment, the tree's clearest statement
of why a guard is keyed on family rather than an internal flag.
`tests/cpp/test_facade.cpp` is the facade's conformance test, one row per
`SamplerBase` virtual driven through the base.

Moves and data - `moves.hpp`, `tree.hpp`, `scan.hpp`, `grow.hpp`,
`data.hpp`: `metropolisJumpForTree`; `Tree`, `columnMaskSubtreeIsValid`;
`scanOrdinalCuts`, `growTreeFromRoot`; `ColumnStore`, `ScopedCutGrid`,
`ColumnKind` and the derived `kindSplitsBySubset`. Judge change-move
detailed balance; the ranked empty-leaf veto in `Tree::leafVetoRank` and
`resolveVetoRank`, where a member-empty leaf vetoes absolutely and a
weight-empty leaf is only penalized (`docs/architecture.md`'s "Tree moves");
whether the semantic kind axis and the mechanical `splitsBySubset` axis stay
separate (only grid construction, ingestion validation and reporting may
read the kind); and the doubled entry layout `scanOrdinalCuts` uses for a
node holding missing members.

The build support files - `configure`, `tools`, `src/misc`, `src/external` -
are skim-only; the one thing worth a look is `simd.c`'s fix for AVX2
misdetected as AVX.

## Appendix B. Other documents

- `docs/architecture.md` - the current state, not a history; prefer it to
  any paraphrase where the two overlap.
- `docs/design/feature-matrix.md` - the one deep read: the per-model
  capability grid, and a Gaps section collecting every missing cell as a
  candidate work item. Its cites are machine-checked; its cell values are
  judgments.
- `docs/design/INDEX.md`, `docs/plans/INDEX.md` - complete manifests,
  refused and closed items included.
- `docs/plans/release-candidate-review.md` - the pre-release review's master
  log, newest first.
- root `TODO` - an alphabetical backlog, some items scheduled after 1.0-0.
  Its `release` item is the one ordered procedure.

The four design documents are proposals with landing notes; the current
design is these sections, about 4,200 of their 15,100 words, and the rest
can be skipped.

- `docs/design/bart-as-a-component.md`, sections 2 "Which mutations are legal
  between sweeps" and 3 "What engine state does not carry, and who
  reinstalls it", about 890 words: which mutations a multi-forest sampler
  admits, and the two state gaps, the per-forest weight and the active-row
  mask.
- `docs/design/multiplier-combiner.md`, the preamble's first paragraph, then
  "The model", "The amplitude layout", "The reparameterization", "The
  amplitude conditional", "bcf as the K = 2 instance", "Surfaces" and "What
  this family does not do", about 1,600 words: what the basis-and-amplitude
  family is, and where BCF sits in it.
- `docs/design/empty-leaf-veto.md`, the preamble's "Correction (2026-07-15)",
  then "Is vetoed-vs-vetoed reachable? Yes; the veto is a RANK (2026-08-18)",
  "What counts as empty: the weight law (2026-08-12)" and "Which weights the
  predicate sees", about 1,370 words: the member-empty versus weight-empty
  ranking.
- `docs/design/bcf.md`, the preamble's model equation and "The multiplier
  snap and the per-forest weight (2026-08-10)", about 355 words: why a row
  can carry an exact-zero weight in one forest, and why that weight is not
  saved state.
