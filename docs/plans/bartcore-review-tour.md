# bartcore: the merge review

Current at d28b087b (bartcore); the INDEX manifests it points at are
machine-checked, this file is not.

Regenerate: re-derive every claim against the new tip, print counts as
commands rather than numbers, then `Rscript tools/check-doc-freshness.R .` -
which certifies `docs/design` anchors and both INDEX manifests. `docs/plans`
anchors, this file's included, are NOT machine-checked, which is why code is
cited by symbol here and docs by section.

Sections 0-5 are the decision; section 6 is optional. Nothing is inlined from
a document that is the authority on it.

---

## 0. The merge in one page

bartcore replaces the classic engine outright. `dbarts/R_C_interface.hpp` and
the C++ ABI behind it are deleted; `inst/include/dbarts/dbarts.h` is the only
shipped header, and it is a flat C API.

**The one structural idea**: the leaf model is a compile-time template
parameter `L` (`accumulate`/`logIntegratedLikelihoodForNode` must inline),
forcing the `facade.hpp` type-erasure layer; the response family is a runtime
virtual chosen once per chain. Almost every "why is this file shaped like
this" question resolves to that sentence.

**Three renames that make greps lie**: `BCFForestCombiner` is now
`AmplitudeForestCombiner` (`combiner.hpp`), one of 20 identifiers a "BCF
sheds its spelling" cleanup renamed (state key `"bcf"` -> `"glue"`) - the old
name finds nothing, which does not mean BCF was removed;
`bartcore_createMultinomial(Counts)` are retired, multinomial going through
the same `bartcore_create` as every other family; and BCF's own R verb,
`bcf()`/`bartBCF`, lives in **bartCause**'s `dbarts-1.0` branch while dbarts
keeps the general engine it was built on.

    git diff --shortstat main...bartcore
    git diff --shortstat main...bartcore -- docs

---

## 1. What a user's script sees

What breaks in a script that ran on main. Every row is an `\item` of
`inst/NEWS.Rd`'s 1.0-0 UPGRADING block, the user-facing contract and the
authority where this table is terser.

| change | evidence |
|---|---|
| Requires R 4.2.0 or newer; the sampler core is C++20 | UPGRADING item 1 |
| Seeded draws differ from 0.9-x; sampling no longer advances R's stream, so R-level draws after a fit change too | UPGRADING item 2 |
| Binary `k` prior becomes `chi(1.5, 2)` (`k = chi(1.5, Inf)` restores the old one); `xbart` now fixes `k = 2` on every family | UPGRADING item 3 |
| Sampler states saved by earlier versions are unrestorable - refit | UPGRADING item 4 |
| A `dbartsData` object saved by an earlier version keeps its old column types and fits differently from one freshly built over the same data frame; an ordered factor in particular has no upgrade path, so rebuild saved data objects | UPGRADING item 5; `docs/plans/column-kind-consolidation.md` section 1 |
| Packages linking against dbarts must port to the flat C API in `dbarts.h`; the old C++ headers are removed | UPGRADING item 6 |
| `bart2`/`rbart_vi` default to `combineChains = TRUE`, merging the chain and sample margins; pass `combineChains = FALSE` for the old shape | UPGRADING item 7 |
| Unordered factors split on level subsets and ordered factors become one ordinal column by default (`bart` unaffected); pass `factors = "indicators"` for the old expansion | UPGRADING item 8 |
| A new `missing` argument keeps and models rows with missing predictors instead of silently dropping them; `NA` in test predictors or `newdata` errors where the training column was complete | UPGRADING item 9; `docs/plans/composition-refusals.md` s7 |
| Test column names not covering the training design's are an error naming what's missing, not a positional-match warning | UPGRADING item 10 |
| An unordered factor response with three or more levels fits multinomial via `bart2` under `family = "auto"`, an ordered one fits ordinal via `dbarts`/`bart2`; the fitters that do not carry the family refuse it by name | UPGRADING item 11 |
| `dbartsControl` drops `rngKind`/`rngNormalKind` and renames `rngSeed` to `seed` | UPGRADING item 12 |
| `xbart` loses its `control` argument, renames `sigma` to `sigest`, and `n.burn` now takes two values | UPGRADING item 13 |
| A `probit` fit refuses `weights`; 0.9-x accepted them and fit a weighted probit | UPGRADING item 14 |
| A wrong-length `weights` vector is an error, not recycled; `resid.prior = fixed(value)` now holds the residual variance at `value`; `getSumsOfSquaredResiduals` returns the raw sum of squares | UPGRADING item 15 |
| The sampler's mutators refresh `$state` only when called with `updateState = TRUE`; `NA` no longer defers to the control's setting | UPGRADING item 16 |
| The sampler's `run` renames `numThreads` to `n.threads` | UPGRADING item 17 |
| `setResponse` gains `updateScale` (default `FALSE`) as its second positional argument, displacing `updateState` to third | UPGRADING item 18 |
| An `rbart_vi` fit with a built-in prior and no callback carries a length-one `$fit` list instead of one sampler per chain | UPGRADING item 19 |
| `predict.bart`/`predict.rbart` share their first five formals `(object, newdata, type, offset, weights)`; `group.by` moves after `...`, name-matched only; `offset.test` is refused by name | UPGRADING item 20; `docs/plans/predict-surface.md` s3 |
| `fitted`'s third positional argument is `ci.level` (was `sample`); `sample` moves to fourth | UPGRADING item 21; `docs/plans/surface-refusals.md` s8 |
| `predict.rbart`'s deprecated `value` alias and `type = "post-mean"` are removed; use `type = "ev"` | UPGRADING item 22 |
| Names foreign to the method called - on `predict.bart`, `extract`, `fitted`, `residuals` - are refused by name instead of silently discarded, and any OTHER unknown name warns rather than vanishing (class `dbartsUnusedArgsWarning`; a warning, not an error, so a subclass forwarding through `NextMethod()` is not refused) | UPGRADING item 23; `refuseUnusedGenericArgs`, `foreignArgsFor`, the five `*ForeignReasons` tables (R/generics.R); `warnUnusedDots` (R/utility.R); surface-refusals.md s4 |
| A fractional double is refused, naming the argument, on every count argument, index and column/group selector a caller can write; `$getSigmas`/`$getSumsOfSquaredResiduals` also refuse their vestigial `result` argument | UPGRADING item 24; `coerceOrError` (R/utility.R); surface-refusals.md s10 |

Decided since: an ordinal fit's R-visible spelling was renamed to
`thresholds`, aligning it with the engine and C API's `ordinalThresholds`.
Section 5.

Also since: an ordered-factor predictor is now its own column kind rather
than a bare ordinal column, and its split grid is the K - 1 midpoints
between consecutive DECLARED levels rather than `n.cuts` uniform cuts over
the observed codes - posterior-changing for any fit with such a predictor,
re-recorded as `equivalence-02d41365.rds` against
`benchmarks/R/categorical-exact.R`'s ordered-factor arm. `n.cuts` no
longer applies to a factor column of either kind, `setCutPoints` refuses
one, and a value that is not an existing level code is refused on both
kinds at every entrance.
`docs/plans/column-kind-consolidation.md`, sections 1 and 6 and its
Landing notes.

---

## 2. What a linked consumer sees

`dbarts.h` is the whole contract; its head comment is where the doctrine
lives. `wc -l inst/include/dbarts/dbarts.h src/C_interface.cpp`

- **A three-class return doctrine**, stated in the preamble and restated per
  entry: a VALUE (the number is the answer), a TRANSACTION result
  (`setPredictor`, `updatePredictor` - 0 is a completed rollback), or a
  CAPABILITY STATUS (0 means this sampler cannot, and nothing was touched).
  `Rf_error` means this call is wrong, and longjmps.
- **Seven setters went `void` -> `int`**: `setResponse`, `setOffset`,
  `setWeights`, `setSigma`, `setTestPredictors`, `setTestOffset`, `predict`.
  Source-compatible at every call site; the answer is a fixed property of the
  sampler, so probe once at setup.
- **`forest` is argument 2** on `dbarts_sampler_getTrees` and
  `dbarts_sampler_printTrees` - the last two that took it after
  `useLiveTrees`. Pointer-to-integer at position 2: fails to compile in C++,
  may only warn in C. The hash is the real backstop.
- `dbarts_sampler_setForestBasis` takes **`basisRowMajor`**, naming the one
  documented exception to the header's column-major rule.
- `dbarts_sampler_setVerbose` **refuses `printEvery == 0`** (the print
  condition is a modulo).
- `dbarts_sampler_setWeights` **validates at the flat entrance**: probit
  refused, logistic counts must be positive integers, gaussian case weights
  finite and non-negative (`enforceBinaryWeightPolicy`, shared with the R
  bridge) - previously R's job alone.
- `dbarts_drawLatents` takes **`ordinalThresholds`**, not `cutPoints`.
- `dbarts_predictor_source` gains **`denseCodes`/`numDenseCodeColumns`**, an
  int32 code channel appended below the 1.0-0 field boundary, so a caller
  whose factor columns are already level codes hands them over unwidened
  (`INT_MIN` marks a missing code). Both channels index at the same
  `columnSources[j]`, in the channel the sampler's kind for that column
  selects; leaving `denseCodes` NULL is what every earlier caller does and
  stays valid.
- `dbarts.h` no longer defines `USE_FC_LEN_T` nor includes `<Rversion.h>`; a
  consumer relying on that pull-in must include it itself.
- `DBARTS_C_API_MAJOR 1` / `DBARTS_C_API_MINOR 0` do not move. The exact-ABI
  token `DBARTS_C_API_HASH` **re-bakes at every ABI change** - a signature, a
  struct field, an enumerator or the callback's parameters, not a header edit
  by itself - and separately has moved several times: since `capi-shape.md`
  recorded its own value, again when the ordered-factor column type was
  appended, and again when the code channel was - so read it from the header
  at the merge tip, never from a doc.

Migration, from `docs/plans/capi-shape.md` s11 (lockstep, after dbarts
installs clean):

| consumer | mandatory edits |
|---|---|
| stan4bart (`bartcore`) | two: the `getTrees` and `printTrees` call sites in `src/init.cpp`, at 979a91a - move `forest` to position 2. Nothing further for the later re-bakes beyond a rebuild; it carries `DBARTS_REQUIRE_EXACT_ABI`, so each one forces one |
| treatSens (`dbarts-1.0`) | none - it calls no reordered entry. NOT R-API-only: its main worktree still links the deleted C++ ABI, which is why the compat branch exists |
| bartCause (`dbarts-1.0`) | none - R API only, no `src/`, zero `dbarts_` symbols. Its `alignForestBasisToSubset` failures were fixed at 765a596; suite green |

---

## 3. Gates, and what each one proves

Two adversarial whole-branch reviews and a calibration lane sit on top of the
mechanical apparatus, and their headline finding is that green checks prove
less than they look like they prove.

| gate | trigger | what it proves | enforced |
|---|---|---|---|
| `check-standard` | push/PR | `R CMD check` at the standing 0E/0W/1N | CI |
| `cpp-tests` | push/PR | the C++ component suite green, and a 7th `ResponseFamily` enumerator a hard compile error (`-Werror=switch`; every `default:` arm deleted) | CI |
| `sanitizers` | push/PR | ASAN+UBSAN over engine and bridge on r-hub instrumented containers; any finding fails | CI |
| `exact-gates` quick | push/PR | 21 exact-posterior and balance scripts - math against a closed form, not a recorded snapshot. Full mode has never run in CI | CI, quick only |
| `exact-gates` cross-host step | push/PR | bcf and multinomial equivalence under `--cross-host` tier 1 LOCKED: `rtol = 1e-8`, `atol = rtol * max|a|` continuous, `identical()` combinatorial | CI |
| `lint` | push/PR | `air format --check`, lintr, `tools/check-rc-codoc.R` | CI |
| `doc-freshness` | push/PR, own workflow, no `paths-ignore` | `docs/design` anchor identity and both INDEX manifests. Split out of `lint.yaml`, whose `paths-ignore` excluded this check's own inputs | CI |
| `equivalence.R` gaussian | schedule + dispatch | 51 scenarios bitwise same-host | dormant |
| `sbc.R` | schedule + dispatch | 83 functionals, Bonferroni band | dormant |
| `rchk` | schedule + dispatch | PROTECT balance | dormant |
| `valgrind` | schedule + dispatch | leaks, OOB reads | dormant |
| `revdep-smoke` | schedule + dispatch | reverse dependencies | dormant |
| `composition-matrix.R` | manual | every `S`/`?` cell of `feature-matrix.md` CONSTRUCTS or REFUSES | manual, one run |
| `bench-sampler.R` | manual, quiet host | speed vs `bench-sampler-ab1dc52.csv` | manual |

**Five are `schedule` + `workflow_dispatch` only - `equivalence`, `rchk`,
`sbc`, `valgrind`, `revdep-smoke` - and have never run, not once, by any
trigger, on any branch**: GitHub binds both triggers to the default branch,
which does not carry these files, so they cannot even be hand-dispatched.
Every clean rchk / valgrind / equivalence / SBC claim rests on a **manual
local run**: rchk clean; valgrind clean after fixing one real leak and one
pre-bartcore OOB read in R's own `makeModelMatrixFromDataFrame.c`; gaussian
SBC PASS. **Merging to main is what registers them** (`TODO`, `release`).

Cross-host is settled and honestly labelled (`docs/plans/bcf-cross-host.md`):
s2 shows tier 2 (a Welch z) cannot gate on its own, and s4 records the x86-box
validation before that box retired - bcf 12/12 tier 1 PASS (worst ratio
3.7e-06), multinomial 11/11 (2.1e-06), `equivalence.R` 43 compared / 0 skipped
all `|z| = 0.00`, tinytest 7478 ok, four discrimination probes firing as
designed. Within a host, reproducibility stays bitwise across SIMD dispatch.

Calibration lane: `sbc.R` at `n.trees = 200`, bracketing both shipped
defaults. **11 families carry an ensemble-scale verdict, 75/83 functionals
PASS, 3 flag, all 3 pre-adjudicated.**

---

## 4. Risks a merge carries

Things that could be wrong and would not be caught.

- **PROTECT balance and leak safety rest on manual local runs**; the two
  workflows that would prove them are among the five dormant ones (s3).
- **No equivalence scenario and no SBC arm reaches a latent K-forest** - BCF's
  sigma clears at 200, and that evidence is gaussian-only.
- **aft and heteroscedastic are uncovered at ensemble scale**; both have real
  sampling code that reduces to no covered family. hazard and hurdle are
  calibration-inherited by bitwise reduction, not directly SBC'd.
- **Warm start and grow-from-root are unrefused and untested at two forests.**
- **The cross-host tier-2 bar is weak by construction**: at the measured ESS
  the `|z| = 4` bar tolerates ~1.4 posterior sd (0.85 of the residual sd on
  `train`), and a follow-up probe confirmed it passes a 20 percent
  node-prior widening tier 1 fails. Tier 1 carries the gate; the honest
  fix, a real seeds axis, is an explicit post-RC door (bcf-cross-host.md
  s2.1, s2.3, landing note).
- **The C++ suite still has measured mutation gaps** in `facade.hpp`,
  `model.hpp` and `moves.hpp`, despite a whole-suite replant that closed most
  of them; counts in `review-2026-08-24/mutation-B-findings.md`.
- **A stale cross-reference in the column-mask coverage**: `test-blocks.R`'s
  warm-start block cites "the BCF columnMask refusal test in
  test-interactions.R", which that file no longer contains. The containment
  gate is still tested from `test-blocks.R`,
  `test-heteroscedastic-warm-start.R` and `tests/cpp/test_state.cpp`; a
  planted mutant making `setState` ignore the containment verdict was
  MISSED, then adjudicated an equivalent mutant and deferred.
- **Six of 32 `benchmarks/R` harnesses are wired to no workflow**, and four
  that call themselves "repeatable" or "permanent" gates have one recorded
  manual run each - `mutation-battery.R`, `grouped-mixing.R` (whose re-measured
  IACT numbers already diverge from its own header, undetected because nothing
  re-runs it), `geweke-mc.R`, `composition-matrix.R`.
- **`setForestBasis(k, ~var)` evaluates the formula in `environment(basis)`
  with `data = NULL`** (`$setForestBasis` -> `evaluateForestBasis`), so a
  data-frame-only column is not found.
- **A per-forest weight is not part of saved state**, and an active-row mask is
  mirrored nowhere: two states can compare `identical()` while fits diverge
  (`docs/design/bcf.md`; pinned only for the same-holder round trip).

---

## 5. Doors that need VD's ruling

Filtered by "does this need a decision, or is it merely unbuilt?" The rest is
in `TODO`.

| door | the undecided part | who decides |
|---|---|---|
| `updateScale` on a multi-forest sampler | **DECIDED**: refusal kept as-is under every family, keyed on bases not family - relaxing later breaks nobody (VD 2026-08-28). Listed so it is not re-opened | decided |
| Real-`r` nbinom | **DECIDED**: scheduled after 1.0-0, value acknowledged (VD 2026-08-28) (`TODO` negbin-real-dispersion). Listed so it is not re-opened | decided |
| Weighted binary | **DECIDED**: scheduled after 1.0-0, shares the decision above (VD 2026-08-28) (`TODO` weighted-binary). Listed so it is not re-opened | decided |
| The RC declaration | VD's own read of this document, then the declaration; both VD-held (`TODO` rc-gate) | VD |
| D7: same-host bcf baseline re-record at the RC tip | a refresh, not a prerequisite - the absent `summaries` field is a pure function of what is stored. Checklist at bcf-cross-host.md s7 | VD schedules |
| Formal heredity (ensemble lattice prior) | **DECIDED**: the first post-1.0 arc, not pre-RC (VD 2026-08-26). Listed so it is not re-opened | decided |

---

## 6. If you read code: the walk

Ordered by what a third party can be broken by. Symbols, not line numbers.

**ABI** - `inst/include/dbarts/dbarts.h`, `src/C_interface.cpp`. The only
thing a third party can be broken by, and the only file a reviewer can
finish. Head comment's contract list, then the X-macro entry table. *Judge*:
whether every non-void entry's doc says which of the three return classes it
is, and whether a DISCARDED capability 0 is a failure mode you accept - it
leaves the sampler unconditioned and the run continuing, quieter than the old
longjmp. `capi-shape.md` s0, s13.

**Engine** - `facade.hpp`, `sampler.hpp`, `chain.hpp`. `SamplerBase` and its
pure virtuals, `SamplerFacade`, the `create*Sampler` factories; `Sampler`,
`run`, `predictColumns` fanning out over `std::thread` workers through
`fanOutPredictSlabs`; `Chain`, `setActiveRows`, `columnMaskStateFeasible`.
*Judge*: the exhaustive `ResponseFamily` switch (no `default:` anywhere), and
that state restore is semantic, not bitwise. Prefer `docs/architecture.md`'s
own prose on RNG and threading.

**Multi-forest** - `combiner.hpp`: `ForestCombiner`,
`AmplitudeForestCombiner`, `MultinomialForestCombiner`.
BCF's `a*mu + b_z*tau` is the K=2 instance of
`docs/design/multiplier-combiner.md`. *Judge*: which mutations the combiner
refuses and why - the mutation-legality table in
`docs/design/bart-as-a-component.md` is required reading first.

**Bridge** - `src/R_interface_bartcore.cpp`: `bartcore_create`, `_run`, the
setters, `_storeState`/`_setState`/`_installForests`,
`_predict`/`_predictPerForest`/`_getTrees`, then the shared guards
`refusedAmplitudeFamilyReason`, `refuseMultiForestMutation`,
`refuseUndefinedTestFits`, `refusePinnedSigmaChange`, `refuseNonBinaryMask`.
*Judge*: `refusePinnedSigmaChange`'s own comment, the clearest statement in
the tree of why a guard is keyed on family rather than an internal flag. The
facade had no conformance test until the second review found 5 of 7 planted
forwarding defects passing the whole C++ suite; `tests/cpp/test_facade.cpp`
is the answer, one row per `SamplerBase` virtual driven through the base.

**Moves and data** - `moves.hpp`, `tree.hpp`, `scan.hpp`, `grow.hpp`,
`data.hpp`: `metropolisJumpForTree`; `Tree`, `columnMaskSubtreeIsValid`;
`scanOrdinalCuts`, `growTreeFromRoot`; `ColumnStore`, `ScopedCutGrid`,
`ColumnKind` and the derived `kindSplitsBySubset`.
*Judge*: change-move detailed balance, the empty-leaf veto at `-HUGE_VAL`,
whether the semantic kind axis and the mechanic `splitsBySubset` axis stay
separate at every site (only grid construction, ingestion validation and
reporting may read the kind), and the doubled entry layout `scanOrdinalCuts`
uses when a node holds missing
members. Build sediment (`configure`, `tools`, `src/misc`, `src/external`) is
skim-only; `simd.c`'s AVX2-misdetected-as-AVX fix is the one thing in it
worth a look.

---

## 7. Where the rest lives

- `docs/architecture.md` - current state, not history; prefer its prose to
  any paraphrase where they overlap.
- `docs/design/feature-matrix.md` - **the one deep read**: 13 family rows, 26
  capability columns across five tables, 49 footnotes, and a Gaps section
  collecting every `M` and `?` cell as a candidate work item. Its anchors are
  machine-checked; its cell VALUES are adjudicated separately.
- `docs/design/INDEX.md`, `docs/plans/INDEX.md` - complete manifests, NO-GO
  and closed items included. The plans index carries a new **Archived**
  section: landed records live under `docs/plans/archive/`, rows kept.
- Landing notes, each authoritative on its own slice: `capi-shape.md` (the
  return doctrine; s11 the migration list), `dbarts-h-freeze.md`,
  `surface-refusals.md` (s4 refusal lists, s17 residue), `rd-records.md`,
  `bcf-cross-host.md` (s2.1-2.2, the most transferable methodology here),
  `composition-refusals.md`, `predict-surface.md`, and
  `pre-review-cleanup.md` (the cleanup rulings and open doors, s5).
- `docs/plans/release-candidate-review.md` - the master log, newest-first;
  sections 1-2 are premises at writing, not current state.
- `docs/plans/review-2026-08-24/` - `consolidated-report.md` (the ledger),
  `gate-ledger.md`, `mutation-{A,B,C,D}-findings.md`, `calibration-sbc.md`,
  `decision-brief.md`; memos and blind critiques under `memos/`.
- root `TODO` - an alphabetical backlog of open items, some scheduled after
  1.0-0; `release` is the one ordered procedure, run when VD triggers it.

Gate commands:

    R CMD INSTALL .            # --preclean after headers/Makevars.in/configure.ac
                               # ALWAYS --preclean after facade.hpp virtuals
    cd tests/cpp && make && ./test_bartcore
    Rscript -e 'tinytest::test_package("dbarts")'

    # equivalence trio - the pre-landing check; only the two cross-host
    # compares run in CI
    Rscript benchmarks/R/equivalence.R compare \
      benchmarks/baselines/equivalence-d4bca4ce.rds
    Rscript benchmarks/R/bcf-equivalence.R compare \
      benchmarks/baselines/bcf-equivalence-00cfa108.rds --cross-host
    Rscript benchmarks/R/multinomial-equivalence.R compare \
      benchmarks/baselines/multinomial-equivalence-4d9a3337.rds --cross-host

    Rscript benchmarks/R/composition-matrix.R
    Rscript benchmarks/R/bcf-exact.R quick    # exact posterior; CI runs quick

    # speed - never concurrent with other load, never in CI
    Rscript benchmarks/R/bench-sampler.R compare \
      benchmarks/baselines/bench-sampler-ab1dc52.csv
