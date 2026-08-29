# bartcore: the merge review

Current at 8e8a63ad (bartcore); the INDEX manifests it points at are
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
| Seeded draws differ from 0.9-x; sampling no longer advances R's stream, so R-level draws after a fit change too | UPGRADING item 1 |
| Saved sampler states from earlier versions are unrestorable - refit | UPGRADING item 2 |
| Family vocabulary refused by name: `bart(family=)` redirects to `bart2`, DART shorthand collisions named | `bartRedirectedFamilies`, `refuseBartRedirectedFamily` (R/bart.R); `resolveDartShorthand` (R/model.R) |
| One `offset` spelling: `predict.bartNegbin`'s `offset.test` is now `offset`, and `offset.test` is refused by name on all six predict methods | UPGRADING, predict-order item |
| `predict`'s first four positions uniform - `(object, newdata, type, offset, ...)` on all six; `group.by` moves after `...`, name-matched only | same item; `docs/plans/predict-surface.md` s3 |
| `fitted`'s slot 3 is `ci.level` on every class; `sample` moves to slot 4 on `bart`/`rbart`, and `fitted.bartHurdle` loses it | `docs/plans/surface-refusals.md` s8 |
| A name foreign to the method called is refused by name instead of vanishing through `...`, positionally as well as named, across predict/extract/fitted/residuals/survivalProbabilities | `refuseUnusedGenericArgs`, `foreignArgsFor`, the five `*ForeignReasons` tables (R/generics.R); surface-refusals.md s4 |
| `residuals(ci.level=)` refused - it never carried a band and mislabelled two of three columns. `type = "class"` with `ci.level`, and `ci.level` with `type = "forest"`, likewise | surface-refusals.md s5-s7 |
| A fractional double is refused, naming the argument, wherever a count is coerced | `coerceOrError` (R/utility.R); surface-refusals.md s10 |
| `n.samples %/% n.thin == 0` refused by name on entry points that return draws; `printEvery` floored at 1 under thinning | `refuseZeroSamples` (R/utility.R); the `printEvery` floor in `bart`/`bart2` control assembly (R/bart.R) |
| `NA` in test predictors or `newdata` is an error where the training column was complete | `docs/plans/composition-refusals.md` s7 |
| A variance forest alongside grouped random effects is a validation error | composition-refusals.md s3 |
| `print.bart`/`print.rbart` gain `\alias`/`\usage`/`\value`; `summary` has methods for every fit class | `docs/plans/rd-records.md` s1; `NAMESPACE` `S3method(summary, *)` |
| Multi-forest prior defaults move; a length-one `forests` list declaring a `basis` is refused | UPGRADING, prior-defaults item |

Kept deliberately: an ordinal fit's R-visible spelling stays `cutpoints`,
though the engine and C API now say `ordinalThresholds`. Section 5.

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
- `dbarts.h` no longer defines `USE_FC_LEN_T` nor includes `<Rversion.h>`; a
  consumer relying on that pull-in must include it itself.
- `DBARTS_C_API_MAJOR 1` / `DBARTS_C_API_MINOR 0` do not move. The exact-ABI
  token is currently `DBARTS_C_API_HASH 0xb6c0e97dc0688991ULL` and **re-bakes
  at every header change** - several times since `capi-shape.md` recorded its
  own value, so read it from the header at the merge tip, never from a doc.

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
| `exact-gates` quick | push/PR | 20 exact-posterior and balance scripts - math against a closed form, not a recorded snapshot. Full mode has never run in CI | CI, quick only |
| `exact-gates` cross-host step | push/PR | bcf and multinomial equivalence under `--cross-host` tier 1 LOCKED: `rtol = 1e-8`, `atol = rtol * max|a|` continuous, `identical()` combinatorial | CI |
| `lint` | push/PR | `air format --check`, lintr, `tools/check-rc-codoc.R` | CI |
| `doc-freshness` | push/PR, own workflow, no `paths-ignore` | `docs/design` anchor identity and both INDEX manifests. Split out of `lint.yaml`, whose `paths-ignore` excluded this check's own inputs | CI |
| `equivalence.R` gaussian | schedule + dispatch | 43 scenarios bitwise same-host | dormant |
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
- **Seven of 32 `benchmarks/R` harnesses are wired to no workflow**, and five
  that call themselves "repeatable" or "permanent" gates have one recorded
  manual run each - `mutation-battery.R`, `grouped-mixing.R` (whose re-measured
  IACT numbers already diverge from its own header, undetected because nothing
  re-runs it), `backfit-exact.R`, `geweke-mc.R`, `composition-matrix.R`.
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
| Ordinal threshold spelling | engine and `dbarts_drawLatents` say `ordinalThresholds`; the R-visible fit and `summary` say `cutpoints`. Keep the split, or align one to the other | VD |
| `updateScale` on a multi-forest sampler | refused under *every* family, keyed on bases rather than family, "though its transform is the identity and the re-anchoring the refusal guards against cannot occur" - a stated divergence from the arc's own plan | VD |
| Real-`r` nbinom | needs a non-integer-shape Polya-Gamma draw, for which no exact sampler exists; shipping it means dbarts' first documented-approximate family behind an explicit opt-in (`TODO` negbin-real-dispersion) | VD, decision-gated |
| Weighted binary | integer-weight probit via replicated latents is built; real weights need a real-shape PG sampler and share the decision above (`TODO` weighted-binary) | VD, decision-gated |
| The RC declaration | the pre-RC slate is complete except the bcf baseline refresh; the declaration is VD-held (`TODO` rc-gate) | VD |
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
`AmplitudeForestCombiner`, `shippedShape()`, `MultinomialForestCombiner`.
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
`scanOrdinalCuts`, `growTreeFromRoot`; `ColumnStore`, `ScopedCutGrid`.
*Judge*: change-move detailed balance, the empty-leaf veto at `-HUGE_VAL`,
and the doubled entry layout `scanOrdinalCuts` uses when a node holds missing
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
- root `TODO` - open work only. `second-review-followups` and `rc-gate` are
  the "what is left" pointers; `release` is the one ordered procedure.

Gate commands:

    R CMD INSTALL .            # --preclean after headers/Makevars.in/configure.ac
                               # ALWAYS --preclean after facade.hpp virtuals
    cd tests/cpp && make && ./test_bartcore
    Rscript -e 'tinytest::test_package("dbarts")'

    # equivalence trio - the pre-landing check; only the two cross-host
    # compares run in CI
    Rscript benchmarks/R/equivalence.R compare \
      benchmarks/baselines/equivalence-736bfb05.rds
    Rscript benchmarks/R/bcf-equivalence.R compare \
      benchmarks/baselines/bcf-equivalence-6e3b9fb8.rds --cross-host
    Rscript benchmarks/R/multinomial-equivalence.R compare \
      benchmarks/baselines/multinomial-equivalence-4d9a3337.rds --cross-host

    Rscript benchmarks/R/composition-matrix.R
    Rscript benchmarks/R/bcf-exact.R quick    # exact posterior; CI runs quick

    # speed - never concurrent with other load, never in CI
    Rscript benchmarks/R/bench-sampler.R compare \
      benchmarks/baselines/bench-sampler-ab1dc52.csv
