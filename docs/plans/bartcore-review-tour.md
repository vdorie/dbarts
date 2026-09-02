# bartcore: the merge review

Current at 18f47358 (bartcore), 2026-09-02.

This is the case for merging bartcore into main. Sections 0-5 are what
you need to rule on it; the code is current at the tip named above.
Section 6 is optional. Code is cited by symbol; `tools/check-doc-freshness.R`
checks every cite.

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

**Three things that make a grep of the old tree mislead**: `BCFForestCombiner`
is now `AmplitudeForestCombiner` (`combiner.hpp`, state key `"bcf"` -> `"glue"`)
- the old name finds nothing, which does not mean BCF was removed;
`bartcore_createMultinomial(Counts)` are retired, multinomial going through
the same `bartcore_create` as every other family; and BCF's own R verb,
`bcf()`/`bartBCF`, lives in **bartCause**'s `dbarts-1.0` branch while dbarts
keeps the general engine it was built on.

Size of the change:

    git diff --shortstat main...bartcore
    git diff --shortstat main...bartcore -- docs

---

## 1. What a user's script sees

`inst/NEWS.Rd`'s 1.0-0 UPGRADING block is the full, authoritative list of
what breaks in a script that ran on main; this is not a copy of it. The
breaks most likely to bite a real script: sampling no longer advances R's
stream, so seeded draws differ from 0.9-x; saved sampler states and
`dbartsData` objects need a version-matched rebuild, not a reload; packages
linking against dbarts must port to the flat C API in `dbarts.h`;
`bart2`/`rbart_vi` default to `combineChains = TRUE`; unordered factors
split on level subsets and ordered factors become one ordinal column by
default; a new `missing` argument keeps rows with missing predictors
instead of dropping them; and a name foreign to the method called is
refused by name instead of silently discarded, across `predict`, `extract`,
`fitted` and `residuals`.

An ordinal fit's R-visible spelling is `thresholds`, aligning it with the
engine and C API's `ordinalThresholds`. An ordered-factor predictor is its
own column kind rather than a bare ordinal column; its split grid is the
K - 1 midpoints between consecutive declared levels rather than `n.cuts`
uniform cuts over the observed codes, which is posterior-changing for any
fit with such a predictor - the current baseline `equivalence-d4bca4ce.rds`
carries it, and `benchmarks/R/categorical-exact.R`'s ordered-factor case is
a closed-form check. `n.cuts` no longer applies to a factor column of either
kind, `setCutPoints` refuses one, and a value that is not an existing level
code is refused on both kinds at every entrance
(`docs/plans/column-kind-consolidation.md`, sections 1 and 6).

---

## 2. What a linked consumer sees

`dbarts.h` is the whole contract; its head comment is where the doctrine
lives.

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
  `dbarts_sampler_printTrees` - the last two entries whose `forest` argument
  still came after `useLiveTrees`. Pointer-to-integer at position 2: fails to
  compile in C++, may only warn in C. The hash is the real backstop.
- `dbarts_sampler_setForestBasis` takes **`basisRowMajor`**, naming the one
  documented exception to the header's column-major rule.
- `dbarts_sampler_setVerbose` **refuses `printEvery == 0`** (the print
  condition is a modulo).
- `dbarts_sampler_setWeights` **is a capability status for four families and
  a value check for two**: returns 0 for probit, ordinal, aft and nbinom
  (none of them carry a weight to change); raises on a logistic count that
  is not a positive integer or a gaussian weight that is negative or
  non-finite (`enforceBinaryWeightPolicy`, shared with the R bridge, which
  raises a named error on all four instead).
- `dbarts_drawLatents` takes **`ordinalThresholds`**, not `cutPoints`.
- `dbarts_predictor_source` gains **`denseCodes`/`numDenseCodeColumns`**, an
  int32 code channel appended after the fields a 1.0-0 consumer compiles
  against, so an older caller's struct layout is unchanged, and a caller
  whose factor columns are already level codes hands them over unwidened
  (`INT_MIN` marks a missing code). Both channels index at the same
  `columnSources[j]`, in the channel the sampler's kind for that column
  selects; leaving `denseCodes` NULL is what every earlier caller does and
  stays valid.
- `dbarts.h` no longer defines `USE_FC_LEN_T` nor includes `<Rversion.h>`; a
  consumer relying on that pull-in must include it itself.
- `DBARTS_C_API_MAJOR 1` / `DBARTS_C_API_MINOR 0` do not move. The exact-ABI
  token `DBARTS_C_API_HASH` **is recomputed at every ABI change** - a
  signature, a struct field, an enumerator or the callback's parameters, not
  a header edit by itself - and has moved several times since it was first
  recorded, so read it from the header at the merge tip, never from a doc.

Migration, from `docs/plans/capi-shape.md` s11 (lockstep, after dbarts
installs clean):

| consumer | mandatory edits |
|---|---|
| stan4bart (`bartcore`) | two: the `getTrees` and `printTrees` call sites in `src/init.cpp` - move `forest` to position 2 (capi-shape.md s11). Nothing further for later hash changes beyond a rebuild; it carries `DBARTS_REQUIRE_EXACT_ABI`, so each one forces one |
| treatSens (`dbarts-1.0`) | none - it calls no reordered entry. NOT R-API-only: its main branch still links the deleted C++ ABI, which is why the compat branch exists |
| bartCause (`dbarts-1.0`) | none - R API only, no `src/`, zero `dbarts_` symbols. Suite green at the current tip |

---

## 3. Gates, and what each one proves

A green gate proves exactly what its row below says, no more - several
enforce less than their name suggests. (Cross-host verdicts have two
tiers: tier 1 a tight relative-deviation bound, the actual gate; tier 2 a
weaker statistical fallback, described further down.)

| gate | trigger | what it proves | enforced |
|---|---|---|---|
| `check-standard` | push/PR | `R CMD check` clean of errors and warnings (the standing 1 NOTE is observed, not enforced); a second job builds and runs the NEON kernels natively on Windows arm64 and asserts scalar-vs-NEON agreement | CI |
| `cpp-tests` | push/PR | the C++ component suite green, and a 7th `ResponseFamily` enumerator a hard compile error (`-Werror=switch`; every `default:` arm deleted) | CI |
| `sanitizers` | push/PR | ASAN+UBSAN over engine and bridge on r-hub instrumented containers; any finding fails | CI |
| `exact-gates` quick | push/PR | 21 exact-posterior and balance scripts - math against a closed form, not a recorded snapshot. Full mode has never run in CI | CI, quick only |
| `exact-gates` cross-host step | push/PR | bcf and multinomial equivalence under `--cross-host` tier 1 LOCKED: `rtol = 1e-8`, `atol = rtol * max|a|` continuous, `identical()` combinatorial | CI |
| `lint` | push/PR | `air format --check`, lintr, `tools/check-rc-codoc.R`, `tools/check-win-drift.R`, `benchmarks/R/mutation-battery.R verify-anchors` | CI |
| `doc-freshness` | push/PR, own workflow, no `paths-ignore` | `docs/design` anchor identity and both INDEX manifests | CI |
| `equivalence.R` gaussian | schedule + dispatch | 51 scenarios bitwise same-host | dormant |
| `sbc.R` | schedule + dispatch | 5 arms, 30 functionals, alpha Bonferroni'd over all 30 | dormant |
| `rchk` | schedule + dispatch | PROTECT balance | dormant |
| `valgrind` | schedule + dispatch | leaks, OOB reads | dormant |
| `revdep-smoke` | schedule + dispatch | reverse dependencies | dormant |
| `composition-matrix.R` | manual | every `S` cell of `feature-matrix.md` actually CONSTRUCTS | manual, one run |
| `bench-sampler.R` | manual, quiet host | speed vs `bench-sampler-ab1dc52.csv` | manual |

**Five are `schedule` + `workflow_dispatch` only - `equivalence`, `rchk`,
`sbc`, `valgrind`, `revdep-smoke` - and have never run, not once, by any
trigger, on any branch**: GitHub binds both triggers to the default branch,
which does not carry these files, so they cannot even be hand-dispatched
(their own header comments claim manual dispatch works on bartcore, and it
does not, for the same reason). Every clean rchk / valgrind /
equivalence / SBC claim rests on a **manual local run**, recorded in
`docs/plans/release-candidate-review.md`: rchk clean
([[docs/plans/release-candidate-review.md#rchk gate: triage and defensive PROTECTs]]);
valgrind clean after fixing one real leak and one pre-bartcore OOB read in
R's own `makeModelMatrixFromDataFrame.c`
([[docs/plans/release-candidate-review.md#valgrind: BCF spec leak, model-matrix OOB read, and the clean full-suite pass]]);
the three equivalence harnesses reproduced bitwise on every scenario
([[docs/plans/release-candidate-review.md#orphaned scripts and the first equivalence-leg run]]);
gaussian SBC PASS
([[docs/plans/release-candidate-review.md#ensemble-scale gaussian SBC, the premise oracle]]).
**Merging to main is what registers these workflows** (`TODO`'s `release`
item).

Cross-host is settled, with tier 2's weakness stated rather than hidden
(`docs/plans/bcf-cross-host.md`): s2 shows tier 2 (a Welch z) cannot gate on
its own, and s4 records the x86-box validation before that box retired -
bcf 12/12 tier 1 PASS (worst ratio 3.7e-06), multinomial 11/11 (2.1e-06),
`equivalence.R` 43 compared / 0 skipped, all `|z| = 0.00` (the corpus's size
at that validation run; it defines 51 scenarios now), tinytest 7478 ok, four
discrimination probes firing as designed. Within a host, reproducibility
stays bitwise across SIMD dispatch.

---

## 4. Risks a merge carries

Things that could be wrong and would not be caught.

- **PROTECT balance and leak safety rest on manual local runs**; the two
  workflows that would prove them are among the five dormant ones (s3).
- **No equivalence scenario and no SBC coverage reaches a K-forest amplitude
  sampler under a latent family (probit, logistic)** - BCF's `sigma` SBC
  functional, which flagged as a calibration defect at the smaller m = 50
  ensemble, passes cleanly at the shipped default of 200 trees
  (docs/plans/review-2026-08-24/calibration-sbc.md), and that evidence is
  gaussian-only.
- **aft and heteroscedastic are uncovered at ensemble scale**; both have real
  sampling code that reduces to no covered family. hazard and hurdle are
  calibration-inherited by bitwise reduction, not directly SBC'd.
- **Warm start and grow-from-root are unrefused and untested at two forests.**
- **The cross-host tier-2 bar is weak by construction**: at the measured ESS
  the `|z| = 4` bar tolerates ~1.4 posterior sd (0.85 of the residual sd on
  `train`), and a follow-up probe confirmed it passes a 20 percent
  node-prior widening tier 1 fails. Tier 1 carries the gate; the honest
  fix, a real seeds axis, is deferred to after the RC (bcf-cross-host.md
  s2.1, s2.3, landing note).
- **`docs/plans/review-2026-08-24/mutation-B-findings.md` (dated 9cebb352)
  measured several C++ mutation gaps in `facade.hpp`, `model.hpp` and
  `moves.hpp`** - among them a forwarding defect in `SamplerFacade`'s
  `setForestWeights`/`savedTree`/`savedSlotForDraw`, `ordinalRuleIsValid`'s
  descendant-interval check, and the ordinal per-observation
  log-likelihood's lower tail. Most have since gained direct unit coverage
  (`tests/cpp/test_facade.cpp`, `test_moves.cpp`, `test_model.cpp`), but the
  mutation record itself has not been re-run against this tip to confirm
  none remain.
- **A stale cross-reference in the column-mask coverage**: `test-blocks.R`'s
  warm-start block cites "the BCF columnMask refusal test in
  test-interactions.R", which that file no longer contains. The containment
  gate is still tested from `test-blocks.R`,
  `test-heteroscedastic-warm-start.R` and `tests/cpp/test_state.cpp`;
  nothing tests that `setState` itself honours the containment verdict
  (`sampler.hpp`'s `allValid = columnMaskOk`); `Chain::stateIsValid`'s own
  per-forest and per-variance-tree checks already reject the same invalid
  states independently, which is why the review-2026-08-24 mutation pass
  ruled this one not blocking (mutation-B-findings.md).
- **Six of 32 `benchmarks/R` harnesses are wired to no workflow at all**:
  `bench-sampler.R`, `composition-matrix.R`, `forest-ranef-collapse-proto.R`,
  `geweke-mc.R`, `grouped-mixing.R` (five harnesses), plus
  `equivalence-common.R`, a shared helper the bcf and multinomial
  equivalence harnesses source rather than a harness of its own. Three of
  the five harnesses call themselves a gate and have one recorded manual
  run each - `bench-sampler.R` ("the zero-regression gate"),
  `composition-matrix.R` ("a local correctness gate"), `grouped-mixing.R`
  ("a fixed gate" to score future mixing work against, whose re-measured
  IACT numbers already diverge from its own header, undetected because
  nothing re-runs it). `forest-ranef-collapse-proto.R` is explicitly a
  prototype and `geweke-mc.R` explicitly disclaims being a gate in its
  poison mode. `mutation-battery.R`'s full run is unwired the same way,
  though it is not one of the six above: only its anchor check is wired,
  into `lint.yaml`.
- **`setForestBasis(k, ~var)` evaluates the formula in `environment(basis)`
  with `data = NULL`** (`$setForestBasis` -> `evaluateForestBasis`), so a
  data-frame-only column is not found.
- **A per-forest weight is not part of saved state**, and an active-row mask
  is mirrored nowhere: two states can compare `identical()` while fits
  diverge (`docs/design/bcf.md`; measured only for a state round trip back
  into the same sampler, not a transplant to another one).

---

## 5. Scope questions: decided and open

Four items are DECIDED and need no ruling: `updateScale` stays refused on
every multi-forest family, keyed on bases rather than family (VD
2026-08-28); real-valued nbinom dispersion and weighted binary
(integer-weight probit, arbitrary real-weight logistic) are both scheduled
after 1.0-0 (VD 2026-08-28; `TODO`'s `negbin-real-dispersion` and
`weighted-binary` items); and formal heredity (an ensemble lattice prior)
is the first piece of work after 1.0-0, not pre-RC (VD 2026-08-26).

One remains open: the RC declaration itself - VD's own read of this
document, then the declaration (`TODO`'s `rc-gate` item).

---

## 6. If you read code: the walk

Ordered by what a third party can be broken by.

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
the tree of why a guard is keyed on family rather than an internal flag.
`tests/cpp/test_facade.cpp` is the facade's conformance test, one row per
`SamplerBase` virtual driven through the base - without it, a forwarding bug
in `SamplerFacade` could pass the whole C++ suite undetected.

**Moves and data** - `moves.hpp`, `tree.hpp`, `scan.hpp`, `grow.hpp`,
`data.hpp`: `metropolisJumpForTree`; `Tree`, `columnMaskSubtreeIsValid`;
`scanOrdinalCuts`, `growTreeFromRoot`; `ColumnStore`, `ScopedCutGrid`,
`ColumnKind` and the derived `kindSplitsBySubset`.
*Judge*: change-move detailed balance, the ranked empty-leaf veto
(`Tree::leafVetoRank`, `moves.hpp#resolveVetoRank`: a member-empty leaf is
an absolute veto, a weight-empty leaf is penalized but survives a chain
whose current state already carries one - docs/architecture.md's "Tree
moves"), whether the semantic kind axis and the mechanic `splitsBySubset`
axis stay separate at every site (only grid construction, ingestion
validation and reporting may read the kind), and the doubled entry layout
`scanOrdinalCuts` uses when a node holds missing members. The build support
files (`configure`, `tools`, `src/misc`, `src/external`) are skim-only; `simd.c`'s
AVX2-misdetected-as-AVX fix is the one thing in it worth a look.

---

## 7. Where the rest lives

- `docs/architecture.md` - current state, not history; prefer its prose to
  any paraphrase where they overlap.
- `docs/design/feature-matrix.md` - **the one deep read**: 13 family rows, 26
  capability columns across five tables, 49 footnotes, and a Gaps section
  collecting every `M` cell as a candidate work item: its cites are
  machine-checked; the cell values are judgments.
- Four design docs sections 4 and 6 call required reading:
  `docs/design/bart-as-a-component.md` (the mutation-legality table),
  `docs/design/multiplier-combiner.md` (the K-forest basis/amplitude
  family), `docs/design/empty-leaf-veto.md` (the ranked veto), and
  `docs/design/bcf.md` (the per-forest weight and active-row-mask gaps).
- `docs/design/INDEX.md`, `docs/plans/INDEX.md` - complete manifests, NO-GO
  and closed items included. The plans index carries an **Archived**
  section: landed records live under `docs/plans/archive/`, rows kept.
- `docs/plans/release-candidate-review.md` - the master log, newest-first;
  sections 1-2 are premises at writing, not current state; its landing
  notes record where each review finding was resolved.
- `docs/plans/review-2026-08-24/calibration-sbc.md` - the calibration
  record behind section 3's SBC numbers.
- root `TODO` - an alphabetical backlog of open items, some scheduled after
  1.0-0; its `release` item is the one ordered procedure, run when VD
  triggers it.

Gate commands:

    Rscript benchmarks/R/equivalence.R compare \
      benchmarks/baselines/equivalence-d4bca4ce.rds
    Rscript benchmarks/R/bcf-equivalence.R compare \
      benchmarks/baselines/bcf-equivalence-3c81d6df.rds --cross-host
    Rscript benchmarks/R/multinomial-equivalence.R compare \
      benchmarks/baselines/multinomial-equivalence-4d9a3337.rds --cross-host
