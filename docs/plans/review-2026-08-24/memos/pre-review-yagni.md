One of four adversarial pre-review reports read ahead of the manual review; VD's rulings on its findings are recorded in docs/plans/pre-review-cleanup.md.

# Adversarial YAGNI review - dbarts @ bartcore 3080a9c5

Read-only. Each candidate is stated with the strongest value claim I could
construct for it, then a verdict. Settled calls (R5 swap surface, the flat C
API's existence and X-macro mechanism, multi-forest/BCF, the nine shipping
families, the multinomial mutation arc, per-observation predictor mutation, the
equivalence/SBC/exact harnesses as gates, docs/design + docs/plans,
run(n.threads), the two declared-not-read predictor-source fields, the
three-class C API return doctrine) are excluded and not relitigated.

Ordered most valuable first. Line counts approximate.

---

## 1. Three thread managers ship; two have zero callers - CUT (~2470 lines)

`src/misc/hierarchicalThreadManager.c` (828), `src/misc/blockingThreadManager.c`
(395), their declarations in `src/include/misc/thread.h` (10 `misc_btm_*` +
~28 `misc_htm_*`), and ~1200 lines of `misc_mt_*` / `misc_htm_*` moment wrappers
and task bodies in `src/misc/moments.c` (74 functions spanning [[src/misc/moments.c:96-2770@74e2e050]] of 2775).

Generalizes over: a two-level (chains x rows) work-stealing thread hierarchy, a
blocking pool, and a threaded variant of every moment function.

Evidence. `grep -rn misc_htm_ src/ tests/ benchmarks/ R/ inst/ | grep -v
'^src/misc/'` is **empty**. Same for `misc_btm_`, outside its own file and its
10 header lines. All 12 `misc_mt_compute*` and all 11 `misc_htm_compute*` public
moment entries in `src/include/misc/stats.h` have zero callers - not in
`src/bartcore`, not in `src/*.cpp`, not in `tests/cpp`, not in
`benchmarks/kernels`. All of it compiles into `misc.a` on every build:
`src/misc/Makefile:53-54` lists both managers and `moments.c` in
`LOCAL_SOURCES`/`LOCAL_OBJECTS`.

What is actually live: `misc_mt_create/runTasks/getNumThreadsForJob/destroy`
from `src/misc/thread.c`, driving the test-fit pool at
`[[src/bartcore/chain.hpp:4066-4083@3080a9c5]]`, plus the single-threaded moments.

Best value claim: "within-chain threading (`docs/design/within-chain-threading.md`)
will need a hierarchical manager." It does not hold. The htm API's shape -
`runTopLevelTasks` + `runSubTask` + `reserveThreadsForSubTask` + a `taskId`
threaded through every moment signature - is the *classic* engine's
chain-outer/row-inner model. bartcore threads chains with `std::thread`
(`[[src/bartcore/sampler.hpp:410@3080a9c5]], 636, 1271`) and rows with `misc_mt_`; a future
within-chain design would be written against those, not against a taskId-keyed C
API inherited from a deleted engine. Nothing tests any of it either, so it is not
validated machinery being held in reserve.

Blast radius: **none**. No test, no consumer, no R path, no benchmark, no doc.
Touches `src/misc/Makefile` and the `.win` object lists.

Adjacent, same header, same verdict (~200 more lines): `misc_computeMeanFast`,
`misc_computeIndexedMeanFast`, `misc_computeWeightedMeanFast`,
`misc_computeIndexedWeightedMeanFast`, the four `*VarianceForKnownMeanFast`
twins, `misc_computeIndexedVariance`, `misc_computeIndexedVarianceForKnownMean`,
`misc_computeIndexedWeightedMean`, `misc_computeWeightedVarianceForKnownMean` -
12 zero-caller public entries. Also `misc_vectorIsConstant`
(`src/include/misc/linearAlgebra.h`), zero callers.

## 2. Dead harness files - CUT (~1000 lines)

- `benchmarks/R/linear-leaf-xcheck.R` (30) and
  `benchmarks/kernels/linear_leaf.cpp` (319): mutually referential, cited
  otherwise only by two `docs/plans/review-2026-08-24/gate-ledger*.md` rows. The
  xcheck prints three numbers nobody compares to anything; the "timing
  denominators" claim fails because no denominator is recorded.
- `benchmarks/R/ordinal-cutpoint-mixing.R` (207) and
  `benchmarks/R/negbin-r-update-mixing.R` (278): prototypes for arcs that have
  landed, whose findings are already written into `docs/design/ordinal.md` and
  `docs/design/negative-binomial.md`. The reproducible-evidence claim is served
  by the design docs; nothing re-runs them.
- `benchmarks/kernels/partition_u8.c` + `.h` (109): zero references outside
  `bench.c`'s include and the Makefile.
- `tools/m4/ax_log1p_in_namespace_std.m4`: `AX_LOG1P_IN_NAMESPACE_STD` is never
  expanded in `configure.ac`, absent from the generated `configure`, and no
  `HAVE_LOG1P*` symbol exists in any `config.*.in`. Provably dead.
- `tools/build-aux/install-sh`: a 0-byte placeholder; `AC_PROG_INSTALL` is not in
  `configure.ac`.
- `[[inst/common/bartcoreHandle.R:425@3080a9c5]]` `bartcorePredictPerForest`: zero callers in
  `inst/tinytest/`, `benchmarks/`, or `R/`.

## 3. `engine=new` and its entire shim apparatus - CUT (~280 lines)

`benchmarks/R/bartcore-shim.R` (129), `tests/cpp/rshim.cpp` (~150 - notably the
one `tests/cpp` source NOT in the Makefile's `SOURCES`), the `engine=new` branch
in `benchmarks/R/equivalence.R`, and the same flag in
`benchmarks/R/bench-sampler.R` where it is parsed, stripped, and discarded.

Generalizes over: running the equivalence gate against a second engine build.

Best value claim: "a second engine (GPU, or a SIMD variant) would compare through
it." Does not hold. The classic engine it was written against is deleted; no
workflow, no doc, no baseline uses the flag (`benchmarks/README.md` prose only);
and a future second engine would compare through the same stored `.rds`
baselines the default path already writes, which is what storing them is for.

## 4. GP leaves - KEEP (value named), MAINTAINER-CALL on cost

`[[src/bartcore/model.hpp:1313-2109@3080a9c5]]` (`GPGaussianLeaf`, 773 lines, the sole
`FunctionLeafModel` instantiation), `[[R/model.R:1452-1490@3080a9c5]]` (`gp()`) plus
column/lengthscale resolution at `[[R/model.R:198-290@3080a9c5]]`, 6 refs in `R/xbart.R`,
`inst/tinytest/test-gp-leaves.R`, `docs/design/gp-leaves.md`,
`docs/plans/gp-followups.md`, `docs/plans/gp-cache-test-flake.md`, a DESCRIPTION
claim.

Generalizes over: leaf functions that are neither constant nor linear.

Best value claim: smooth-in-x leaves are a published BART variant, R-reachable
and tested, and the only reason the `FunctionLeafModel` concept exists. Under the
standing rule (gate only on failure to name value), that **holds** - this is not
YAGNI by the stated bar.

It earns a line anyway as the highest cost-per-user item in the branch: 773
engine lines carrying an O(n_leaf^3) Cholesky, a per-draw weight cache with its
own byte budget (`statisticsCacheBudget_`), an oversized-leaf constant-leaf
fallback with its own MH-coherence argument, a known test flake with its own plan
doc, and zero consumers on any branch. If one thing here warrants a deliberate
"is this in 1.0" conversation rather than a cut, it is this.

## 5. The flat C API's forest/amplitude block - KEEP the block; MAINTAINER-CALL on two members

`inst/include/dbarts/dbarts.h` X-macro entries [[CAPI:573@3080a9c5]] `numForests`, [[CAPI:575@3080a9c5]]
`setForestBasis`, [[CAPI:579@3080a9c5]] `getForestFits`, [[CAPI:582@3080a9c5]] `numForestAmplitudes`, [[CAPI:584@3080a9c5]]
`getForestAmplitudes`, [[CAPI:588@3080a9c5]] `setForestWeights`, [[CAPI:591@3080a9c5]] `getForestCalibration`,
[[CAPI:595@3080a9c5]] `setForestPriorScale`.

Zero calls across stan4bart (which uses 24 of 48 entries), treatSens dbarts-1.0
(10), bartCause, bairrtt - and zero mentions in stan4bart's own docs/TODO. The
only caller is `inst/tinytest/capi/consumer.c`, which calls 48/48 by construction
(`[[docs/plans/capi-shape.md:204@3080a9c5]]`).

Best value claim, and it is real: `[[docs/plans/dbarts-h-reshape.md:1264-1267@3080a9c5]]`
ships `setForestWeights` on *named enabled models* (multilevel BCF,
latent-treatment sensitivity, zero-inflated log-linear BART). That is the
enabling-value gate working as specified, and consumer absence is not the gating
fact. **The claim holds for the block.**

Weakest for two members:
- `dbarts_getForestCalibration` ([[docs/plans/dbarts-h-reshape.md:591@3080a9c5]]) - 39 lines of implementation, its own
  13-member ABI struct, its own `_INIT` macro, its own field-boundary contract,
  ~35 lines of Doxygen. Largest surface-per-caller ratio in the header.
- `dbarts_sampler_setForestPriorScale` ([[docs/plans/dbarts-h-reshape.md:595@3080a9c5]]) - its own documentation concedes a
  read-then-write of the calibration is inert.

If cut as a pair (~90 lines + struct): `test-capi.R` legs, `capi/consumer.c`,
`docs/plans/capi-shape.md` anchors, the 48-count cited in several docs.

## 6. `drawShippedGlue` - a duplicate amplitude conditional - FINISH-OR-CUT (~125 lines)

`[[src/bartcore/combiner.hpp:1141@3080a9c5]]` (impl), `[[src/bartcore/combiner.hpp:985@3080a9c5]]` (branch), `[[src/bartcore/combiner.hpp:1478@3080a9c5]]`
(`shippedShape()`), with `generalAmplitudeDraw` plumbing at `[[src/bartcore/combiner.hpp:348@3080a9c5]]`, `[[src/bartcore/combiner.hpp:759@3080a9c5]]`,
`[[src/bartcore/combiner.hpp:1518@3080a9c5]]`.

Generalizes over: nothing. It is the K=2 special case of `drawAmplitudes`, which
already handles it. Real setters of `generalAmplitudeDraw`: **0 hosts, 3 test
refs**. `shippedShape()` selects the duplicate whenever BCF's canonical basis
holds - i.e. always, on the shipped path.

The comment at `[[src/bartcore/combiner.hpp:960-972@3080a9c5]]` names its own deletion trigger verbatim ("at which
point `AmplitudeSpec::generalAmplitudeDraw` becomes the default and this branch
is deleted") and records the measurement behind it (21 accumulation variants
tried; the general path cannot be made bitwise on bcf).

Best value claim: it preserves the bcf equivalence baseline bitwise. Holds today
- which makes this a baseline re-record decision, not a design one. Blast radius:
`benchmarks/baselines/bcf-equivalence-6e3b9fb8.rds` + MANIFEST.

## 7. Gates written and never wired - FINISH

The inverse smell: a gate apparatus whose only path is the manual one.

- `.github/workflows/{equivalence,rchk,revdep-smoke,sbc,valgrind}.yaml` are
  `schedule` + `workflow_dispatch` only. GitHub binds `schedule` to the default
  branch, and **none of the five exists on `main`** - so on `bartcore` none has
  ever fired automatically. rchk, valgrind, SBC, revdep and the gaussian
  equivalence arm have zero automated coverage today. They arm on merge; until
  then, calling them gates is false.
- `benchmarks/R/backfit-exact.R` (345) is the only `*-exact.R` omitted from
  `exact-gates.yaml`'s 20-gate list, despite the same `quick` flag and in-script
  tolerance. Wire it.
- `benchmarks/R/composition-matrix.R` (766) is an executable oracle diffing
  family x capability against `docs/design/feature-matrix.md`, and
  `tools/check-doc-freshness.R` does not do this. Wire it, or it is 766 lines of
  oracle nobody consults.
- `tools/check-doc-freshness.R` Part 4 (baseline-hash resolvability) self-skips
  under CI's shallow checkout, so ~1/5 of a 1771-line linter has never executed
  in CI. Fetch the tag or drop Part 4.
- `benchmarks/R/mutation-battery.R` pins `equivalence-31e52644.rds`, a MANIFEST
  role=`historical` baseline - the only live reader of a superseded file. Repin.
- `benchmarks/R/bench-sampler.R` `biggrid` is checked *before* mode dispatch, so
  `biggrid record` silently writes a CSV no baseline slot exists for.
- `benchmarks/R/sbc.R` accepts 26 `which` values; `sbc.yaml` runs 6. Twenty
  (`probit`, `dart`, `dart-sparse`, `weighted`, `linear*`, `gp*`, `aft`, `bcf`,
  `bcf-weak`, `burn-<family>`) run nowhere automatically. MAINTAINER-CALL: a
  hand-runnable SBC menu is defensible; 20 modes unrun since they landed is a
  maintenance surface, not a gate.
- 6 of the 11 stored baselines are role `historical` and reproducible by the
  current engine (only `equivalence-5430fdb.rds`, the classic-engine cutover
  evidence, is genuinely unreproducible). ~1 MB of archive presented as gate
  material. MAINTAINER-CALL.

## 8. The data-handle column-subset view is dead end to end - CUT (~50 lines)

`[[R/bartcore.R:687-709@23b9cde7]]` (`bartcoreSamplerFromHandle(..., columns = NULL)`) and its
bridge at `[[src/R_interface_bartcore.cpp:3622@23b9cde7]], 3694-3720` - the 8th argument of
`C_dbarts_bartcore_createFromHandle`, its bounds checks, and the
leaf-covariate-to-view-column remap that exists *only* for the subset case
(`viewLeafCovariates`, `[[src/R_interface_bartcore.cpp:3707-3720@23b9cde7]]`).

Generalizes over: a fold-view restricted to a subset of the handle's predictor
columns.

Evidence: `bartcoreSamplerFromHandle` has exactly one production caller
(`[[R/xbart.R:702@23b9cde7]]`), which passes seven positional arguments and never reaches
`columns`; the six test call sites (`[[test-data-handle.R:32@23b9cde7]], 47, 114`,
`[[test-linear-leaves.R:190@23b9cde7]], 207`, `[[test-gp-leaves.R:237@23b9cde7]], 254`) never name it
either. So the whole path - R formal, C argument, remap - is unreachable.

Best value claim: xbart could one day cross-validate over predictor subsets as
well as over (k, power, base, n.trees) cells, and the handle exists precisely so
folds share one binning. That is a nameable enabled use, and it is a good one -
which makes this **FINISH-OR-CUT rather than a clean cut**: either wire a
`columns` axis into xbart's cell grid, or delete the argument and the remap and
let the view span the handle. What it must not stay is a general mechanism whose
only path is the full-span special case, with an untested remap in the middle.

## 9. `AmplitudeSpec`'s mu/tau spelling is a test-only fixture - FINISH-OR-CUT (~56 engine lines)

`[[src/bartcore/combiner.hpp:316-353@23b9cde7]]` (fields), `[[src/bartcore/combiner.hpp:362-378@23b9cde7]]` (`expandForestSpecs`),
`[[src/bartcore/combiner.hpp:1424@23b9cde7]]` (`synthesizeIndicatorBasis`), `[[src/bartcore/combiner.hpp:781@23b9cde7]]` (guard).

Its own comment at `[[src/bartcore/combiner.hpp:326@23b9cde7]]` says it: "a FIXTURE default that no consumer reaches -
`createAmplitudeSampler` has no flat-C entry point, and the bridge fills
`spec.forests` unconditionally." Confirmed at `[[src/R_interface_bartcore.cpp:2157@23b9cde7]]`.
Host writers of `mu`/`tau`/`z`/`aPriorScale`/`bPriorVariance`/`sdModerate`/
`ridgeA`/`ridgeB`/`updateB`: **zero**. Test writers: 25-34 refs each.

Best value claim: the readable spelling for the two-forest case in C++ fixtures.
That is a test-ergonomics claim, not a capability one - and it buys a second
entrance to what `spec.forests` already says, which is the shape that drifts. Cut
is right; the cost lands on ~90 C++ fixture sites.

## 10. `sliceSample`/`rejectionSample`: the `log = FALSE` arm is test-only - CUT (~70 lines)

`[[R/sliceSample.R:1-44@23b9cde7]]` (`rejectionSample`, two parallel implementations of the
same loop) and the `useLog` branches through `sliceSample` (313 lines total). The
single production call site is `[[R/rbart.R:676@23b9cde7]]`, at the default `log = TRUE`.
Every `log = FALSE` use is in `inst/tinytest/test-slice-sample.R` (5 sites).

Generalizes over: densities supplied on the natural rather than the log scale.

Best value claim: `rbart_vi(prior = )` takes a user closure
(`test-rbart-custom-prior.R`), so a user might hand back a natural-scale density.
Does not hold - `[[R/rbart.R:19@23b9cde7]]` documents the contract inline ("on log scale"),
the closure is composed with a log likelihood inside `posteriorClosure`, and a
natural-scale prior would be numerically wrong there regardless of what the
sampler can consume.

Blast radius: `test-slice-sample.R` only (internal, not exported).

## 11. `unpack` / `[<-.named_lval` - a destructuring facility with zero users - CUT

`[[R/multipleAssignment.R:118-139@74e2e050]]` (`unpack` and its method) has **zero** uses
outside `inst/tinytest/test-multipleAssignment.R`. Its sibling `massign` has
exactly four uses, all `[[R/partialDependence.R:44@23b9cde7]], 75, 215, 329`, all of the form
`massign[sampler, fit] <- ...` - two-element positional destructuring. The
named-argument branch of `[<-.lval` ([[R/partialDependence.R:55-97@23b9cde7]], with duplicate-name warnings,
multiple-match warnings, and pop-if-not-referenced-later bookkeeping) is likewise
reached only from the test file.

Generalizes over: arbitrary named/positional multiple assignment with R-level
name resolution and collision diagnostics.

Best value claim: a general internal convenience. Not exported, so no user
benefit exists; the one real use is served by two ordinary assignments; and
`[[R/partialDependence.R:214@23b9cde7]], 328` already carry `sampler <- fit <- NULL ## for R
CMD check` lines that exist only because `massign` hides the assignment from the
static checker - the abstraction taxes every call site it has.

CUT `unpack` outright (~25 lines); MAINTAINER-CALL on the named branch of
`massign` (~45 more). `test-multipleAssignment.R` (99 lines) largely disappears.

## 12. Formals that exist in order to be ignored - CUT (~30 lines)

- `[[R/bart.R:2937-2940@23b9cde7]]` **`makeind(x, all = TRUE)`**: the body is
  `ignored <- all ## for R check` and then `makeModelMatrixFromDataFrame(x,
  TRUE)`. `[[man/makeind.Rd:25@23b9cde7]]` documents `all` as "not currently implemented;
  retained for signature compatibility with `BayesTree::makeind`", and
  `[[inst/tinytest/test-makeModelMatrix.R:299@23b9cde7]]` asserts
  `makeind(df, all = FALSE)` equals `makeind(df)` - a test whose only purpose is
  to certify that the parameter does nothing. The value claim (signature compat)
  inverts: a BayesTree user passing `all = FALSE` gets silently *different*
  behavior rather than an error, which is worse than not accepting it. FINISH
  (implement) or CUT (drop the formal); the current state is the worst of the
  three.
- `[[R/dbarts.R:1675-1689@23b9cde7]]` **`getSigmas(result)`** and `[[R/dbarts.R:1695-1708@23b9cde7]]`
  **`getSumsOfSquaredResiduals(result)`**: both declare `result` solely to
  `stop()` when supplied, reasoning that the name "cannot be quietly
  repurposed." Dropping the formal makes R itself raise "unused argument
  (result)", naming the same thing - and the symmetry with `getLatents(result)`
  is what creates the confusion the refusal then corrects.
- `[[R/utility.R:899-905@23b9cde7]]` **`quoteInNamespace(name, character.only = FALSE)`**: 9
  call sites, none supplies `character.only`, so the `TRUE` branch is dead.

## 13. `dbarts_sampler_family` exists so its enumerator can exist - FINISH-OR-CUT (~25 lines)

`[[inst/include/dbarts/dbarts.h:614@23b9cde7]]`. `[[docs/plans/dbarts-h-freeze.md:42@23b9cde7]]` states it
outright: the MULTINOMIAL arm is "dead until multinomial creation opens, the
point being the ENUMERATOR existing pre-freeze."

Best value claim: the enumerator must be placed before the ABI freezes. Nothing
is frozen pre-release (standing rule), so that carries no weight here. And the
freeze doc records at `[[docs/plans/dbarts-h-freeze.md:198-199@23b9cde7]]` that the accessor fails at the diagnostic it was
meant to serve - a heteroscedastic sampler answers GAUSSIAN and `setSigma` still
refuses it. Strip the unreachable arm and the entry only reports back what the
caller passed to `create`, modulo AUTO. Either make it answer the questions a
host actually asks (heteroscedastic, multi-forest, which refusals apply) or drop
it.

## 14. `setCallback` + `getDispersion` - KEEP the hook; MAINTAINER-CALL on the reader

`[[inst/include/dbarts/dbarts.h:506@23b9cde7]]` and `[[inst/include/dbarts/dbarts.h:600@23b9cde7]]`.
`[[docs/plans/archive/capi-callbacks.md:9-13@23b9cde7]]` names the outer-Gibbs host (stan4bart), then
`[[docs/plans/archive/capi-callbacks.md:84-85@23b9cde7]]` puts migrating stan4bart out of scope; `[[docs/plans/archive/capi-callbacks.md:23-24@23b9cde7]]` records that callbacks
had previously been deferred until a consumer existed.

`setCallback`'s claim holds: a per-sweep hook is strictly cheaper than
`run(0,1)` + `storeState` for a host driving one sweep at a time, which is this
package's headline use case.

`dbarts_sampler_getDispersion` is weaker. Its own header text concedes it returns
"the same scalar `dbarts_results::dispersion` records"; its only non-redundant
use is a mid-sweep read without serializing state, reachable **only from inside a
callback**. Speculative conditional on a surface that is itself unexercised.
(~20 lines; the R5 `$getDispersion()` is separate and kept.)

## 15. `SamplerOptions::forestColumns` - a per-forest mask nothing writes - CUT (~18 lines)

`[[src/bartcore/chain.hpp:127-132@23b9cde7]]` (decl), `[[src/bartcore/chain.hpp:694-701@23b9cde7]]` (consume), `[[src/bartcore/chain.hpp:764@23b9cde7]]` and `[[src/bartcore/chain.hpp:878@23b9cde7]]`
(null-outs). R's `vars =` restriction routes exclusively through
`ForestStructureSpec::columns` (`[[R/spec.R:726@23b9cde7]]` ->
`[[src/R_interface_bartcore.cpp:2197@23b9cde7]]`). Writers of the `SamplerOptions` twin: 0 in
`src/`, 0 in `R/`, 3 in `tests/cpp`.

The two null-out lines exist only to defend against a field nothing sets - the
incomplete-generality smell in miniature: a second entrance that must be defended
against precisely because it is never used.

## 16. Knobs whose branches no test reaches - FINISH

Not YAGNI (each names a value), but each is a code path with no exerciser, which
is how a general mechanism quietly becomes a broken one.

- `[[R/model.R:1633@23b9cde7]]` **`dart(rho = NULL, b = 1)`**: never supplied anywhere in
  `R/`, `inst/tinytest`, vignettes, `man`, or the four consumers (`dart(a = )`
  *is* tested). `rho <= 0` falls back to `numPredictors`, so the non-default
  branch at `[[src/bartcore/model.hpp:2446@23b9cde7]], 2454` never executes in any test.
  Value: they are Linero (2018)'s own sparsity hyperparameters. Add coverage.
- `[[R/validateComposition.R:52@23b9cde7]], 142` **`seed`, `alpha`**: never supplied. `alpha`
  drives the nSim floor ([[R/validateComposition.R:81@23b9cde7]]), the band quantile ([[R/validateComposition.R:93@23b9cde7]]), and the Bonferroni
  correction ([[R/validateComposition.R:225@23b9cde7]]), and its own error message ([[R/validateComposition.R:172@23b9cde7]]) is untested; `seed` is the
  only thing making a verdict reproducible. This is a *calibration validator* -
  its own knobs going unexercised is the sharpest version of the problem.
- `[[R/dbarts.R:1986@23b9cde7]]` **`printTrees(chainNums)`** has no default and is never
  supplied; resolved through `match.call()` at [[R/dbarts.R:1988@23b9cde7]].
- `[[R/bartcore.R:740@23b9cde7]]` **`mu.blocks`**: reached only from
  bartCause's `tests/testthat/test-14-bcf.R` line 353; absent from `inst/tinytest`
  entirely.

## 17. The per-generic "foreign argument" reason tables - KEEP (value named)

`[[R/generics.R:2050@23b9cde7]]` (`foreignArgsFor`) plus five tables at `[[R/generics.R:2109-2164@23b9cde7]]`
(`predictForeignReasons`, `extractForeignReasons`, `fittedForeignReasons`,
`residualsForeignReasons`, `survivalProbabilitiesForeignReasons`), composed at 22
call sites, ~130 lines of prose strings.

Generalizes over: every (generic, argument) pair that is a formal on one method of
a generic and foreign on another.

Best value claim: error quality is a stated commitment
(`docs/design/error-style.md`, `docs/plans/robust-errors.md`,
`inst/tinytest/test-error-quality.R`), and `extract(object, newdata = )` being
silently swallowed by `...` is a genuinely bad failure. **Holds.**
MAINTAINER-CALL only on whether the table formulation earns itself over an inline
`stop()` per method; it is proportionate to 6 fit classes x 5 generics. Same
verdict for the 8 one-line refusal methods at `[[R/generics.R:3069-3078@23b9cde7]],
3096-3105`.

## 18. `PredictorUpdateResult::unsupportedSource` - produced, never distinguished - CUT (~4 lines)

`[[src/bartcore/sampler.hpp:86@23b9cde7]]`, produced at `[[src/bartcore/sampler.hpp:1383@23b9cde7]]` and `[[src/bartcore/sampler.hpp:1397@23b9cde7]]`. Readers: zero.
`[[src/C_interface.cpp:238@23b9cde7]]` collapses it into the `accepted ? 1 : 0` mapping;
`[[inst/tinytest/test-capi.R:1212@23b9cde7]]` documents it as counterfactual ("Handing the
engine a non-dense view instead *would* return..."). The refusal is
load-bearing; only the distinct enumerator is dead. Fold into `rolledBack` or an
assert.

## 19. `SamplerShape::varianceLeafPrior` - one reader, and it is tautological - CUT (~5 lines)

`[[src/bartcore/facade.hpp:86@23b9cde7]]`, filled at `[[src/bartcore/facade.hpp:443@23b9cde7]]`. Sole reader:
`[[tests/cpp/test_facade.cpp:521-522@23b9cde7]]`, which asserts the shape field equals the
`impl_` accessor it was copied from one line earlier. No R reader, no `dbarts.h`
reader, no behavioral test. (+ `[[chain.hpp:918@23b9cde7]]` / `[[sampler.hpp:818@23b9cde7]]` accessors if
nothing else reads them.)

## 20. `bart()`'s lowercase formals and the plotting knobs - KEEP (value named)

`[[R/bart.R:2736@23b9cde7]]` `sigquant`, `printevery`, `numcut`, `proposalprobs`; also
`[[R/dbarts.R:2115@23b9cde7]]` `treePlotPars`, `[[R/plot.R:476@23b9cde7]], 505` `cols`/`justmedian`,
`[[R/partialDependence.R:200@23b9cde7]], 314` `plquants`, `[[R/diagnostics.R:177@23b9cde7]]+` `vars` at the
`as_draws_*` entrances. None is ever supplied by name anywhere - not in `R/`, not
in `inst/tinytest`, not in vignettes, not in `man` examples, not in any of the
four consumers.

Best value claim: `[[man/bart.Rd:17@23b9cde7]]` states that `bart` **is** "the
BayesTree-compatible surface", and DESCRIPTION claims dbarts "serves as a
drop-in replacement for package 'BayesTree'". `sigquant`, `printevery` and
`numcut` are BayesTree's own spellings; the pdbart plotting knobs are BayesTree's
too. A drop-in claim that silently drops arguments is not a drop-in claim.
**Holds.**

Two caveats worth a decision rather than a cut: `proposalprobs` is a dbarts-0.9
name, not a BayesTree one, so the compat claim does not cover it; and `sigquant`
is the one on the list with genuinely *no* coverage anywhere, which for a
prior-calibration argument is a testing gap more than a design one.

## 21. `dbartsControl@storage = "single"` - KEEP (value named)

`[[R/A_class.R:236@23b9cde7]], 256, 332-338`; `[[R/dbarts.R:218@23b9cde7]], 231, 238`; consumed at
`[[src/R_interface_bartcore.cpp:170@23b9cde7]]` (an opt-in fp32 running residual); refused at
`[[R/spec.R:450@23b9cde7]], 678`.

Best value claim: memory-wall relief at large n, a standing fact here
(single-chain n >= 1e5 is common), and a genuine second `ResidT` instantiation in
the engine rather than a stub. Holds. Named only because it is the one control
slot whose two values change numerical output, so its two `spec.R` refusal sites
are where an incomplete-generality regression would first appear.

## 22. `family = "twopart"` - KEEP (value named)

`[[R/dbarts.R:388@23b9cde7]], 408-410`; `[[R/bart.R:705@23b9cde7]], 721-725, 2698-2704`; three Rd files;
three test files. Generalizes over two spellings of `hurdle.lognormal`.

Best value claim: "two-part model" is the econometrics/health-economics name for
exactly this model, so the alias meets users at their own vocabulary; it resolves
to the canonical token immediately, costing no downstream branching. Holds.
Flagged only because it is the sole alias in the family vocabulary and
pre-release is when alias policy is cheapest to revisit.

## 23. Small stale or incorrect surface

- **`[[R/sparseFactor.R:1-7@23b9cde7]]`**: the header comment says the class "is recognized by
  ingestion but refused at data construction until sparse-categorical (CSC over
  level codes) engine support lands; it exists so data sets can be assembled
  ahead of engine support." `man/sparseFactor.Rd` and `[[R/utility.R:437-450@23b9cde7]]` say
  the opposite - it rides to the engine unexpanded and bins bitwise-identically
  to a dense factor. **The comment is stale**, and it is precisely the comment
  that would have made `sparseFactor` the headline YAGNI finding here. Fix it;
  the feature's own value (high-cardinality categorical memory) is named and
  holds.
- **`DESCRIPTION`**: `person("Armon", "Dadgar", role = "ctb", comment = "adaptive
  radix tree")`. There is **no radix-tree code anywhere** in `src/` or `R/` -
  `grep -rln "radix\|art_tree" src/ R/` is empty; it went with the classic
  engine. A contributor credited for code that does not ship. CUT the credit (or
  restore the code).
- **`[[R/model.R:1556-1567@23b9cde7]]`**: `normal(k = "chi(1.5)")` string-spec parsing,
  introduced as "compatibility with string specifications." Reached from
  `[[R/rbart.R:110@23b9cde7]]` and two tests. Pre-release there is no compatibility to
  preserve, and the call-shaped NSE path `normal(chi(1.5))` already exists.
  MAINTAINER-CALL, ~12 lines.
- **`src/include/external/random.h`**: `ext_rng_createAndSeed`,
  `ext_rng_seedsAreEqual`, `ext_rng_getAlgorithmName`, `ext_rng_getState`,
  `ext_rng_setState`, `ext_rng_setSeedFromClock`, `ext_rng_usesNativeRNG`,
  `ext_rng_simulateNormal`, `ext_rng_simulateLowerTruncatedStandardNormal` -
  zero callers outside `src/external/`. ~120 lines of public API.
  MAINTAINER-CALL (a vendored library keeping its own API is normal; this one is
  no longer vendored from anywhere upstream).
- **`inst/common/`**: `pdData.R` (1 test), `almostLinearBinaryData.R` (2),
  `multithreadData.R` (3) - shared fixtures with fewer users than the
  indirection costs. Trivial, MAINTAINER-CALL.
- **Not YAGNI, but found on the way and worth its own look**: `[[R/dbarts.R:959@23b9cde7]]`
  `run(..., n.threads = control@n.threads)` is never read in the body -
  `bartcoreSamplerRun(.self, numBurnIn, numSamples)` (`[[R/bartcore.R:145@23b9cde7]]`) takes
  no thread count - while `[[man/dbartsSampler-class.Rd:56@23b9cde7]], [[man/dbartsSampler-class.Rd:154@23b9cde7]]` describe it on
  the same footing as `predict`'s, which *does* honor it (`[[R/dbarts.R:1145@23b9cde7]]`).
  The reserved status is settled; the doc wording appears not to match it.

## 24. Checked and clean (recording the negative results)

- **No dead R functions.** A parse-based sweep of all 381 top-level definitions
  in `R/` against `R/`, `inst/tinytest/`, `inst/common/`, `man/`, `vignettes/`
  and `NAMESPACE` found exactly two names occurring once: `.onLoad` and
  `.onUnload`, both hooks R itself calls. An independent 146k-call-site sweep
  including the four consumers agrees: zero dead internal functions, zero
  exported-and-uncalled functions.
- **No zero-user R5 methods.** All 31 `dbartsSampler` mutators/readers have
  tinytest coverage; `reapplyForestWeights` is internal (called from
  `getPointer`/`setState`). Five have live consumer callers.
- **No dead `run()` result channels.** All 18 - `sigma`, `train`, `test`,
  `varcount`, `k`, `varprobs`, `tau`, `ranef`, `cutpoints`, `dispersion`,
  `resid.df`, `variance`, `varianceTest`, `forestFits`, `glue`, ... - have R
  readers and test assertions.
- **No zero-registrant engine hooks.** `SweepCallback` has 2 registrants
  (`rbart_vi` at `[[src/R_interface_bartcore.cpp:4594@23b9cde7]]`, flat C at
  `[[src/C_interface.cpp:605@23b9cde7]]`); `pollInterrupt` 1 with 3 call sites;
  `drawForestGlue` and `combinedTestFits` 1 override each.
- **Every engine policy template has >= 2 instantiations**: `L`/`ResidT` (5),
  `Merge` (2), `Columns` (2), `Ensemble` (2), `Observer` (2), `GoesRight` (6),
  `FitsOf` (2). The single-instantiation *concepts* (`VectorLeafModel`,
  `FunctionLeafModel`, `ScaleLeafModel`, `TreeDrawLeafModel`,
  `ParamScoringLeafModel`) each have a live R path, and the last two exist so the
  constant-leaf path stays bitwise - a present-tense reason.
- **No `#if 0`, no experimental `#ifdef`s, no dead `constexpr bool` flags** other
  than `generalAmplitudeDraw` (finding 6). The unreachable enumerator arms in
  `AmplitudeChain`'s family switch (`[[src/bartcore/chain.hpp:775-789@23b9cde7]]`) are
  deliberate: adding a family fails the build rather than silently acquiring the
  gaussian law. Correct as written.
- **`dbartsSpec(parentEnv = )`** looks unused in-repo but **is** supplied by
  stan4bart (`R/stan4bart_fit.R` line 556). Keep.
- `misc/partition.h`, `misc/simd.h`, `misc/io.h`, and `misc/linearAlgebra.h` (bar
  `misc_vectorIsConstant`) are fully used.

---

## Overall judgement

Relative to a typical 1.0, this branch carries **less** speculative generality
where one expects it and **more** in one place nobody looked. The engine is
unusually tight: every policy template has multiple live instantiations, every
single-instantiation concept has an R path, all 18 `Results` channels and 30 of
31 `SamplerShape` fields have real readers, there are no zero-registrant hooks,
no dead feature flags, no `#if 0`. What generality the engine does carry is
mostly *legacy duplication kept for bitwise equivalence* (findings 6, 9, 15) - a
cost the gate discipline imposed rather than an imagined future - and it is
self-documenting, with the deletion trigger written into the comment. The flat C
API's speculative residue (5, 13, 14) is bounded and deliberate; the plans show
the enabling-value gate biting in both directions, reserving flat multinomial
create and forest-indexed predict as doors rather than building them. The R
surface has essentially no dead code at all: 381 function definitions, zero
unreachable. The genuine outlier is `src/misc`, roughly 2500 lines - two entire
thread managers and 74 threaded moment functions - compiled into `misc.a` on
every build with zero callers anywhere, untested, undocumented, and structurally
matched to the deleted classic engine rather than to bartcore; that one item is
larger than every other cut on this list combined, and it survived only because
it sits below the layer everyone was reviewing. The second-order pattern is not
excess machinery but *unwired* machinery (findings 7, 8, 16): gates written and
never scheduled, an oracle nobody consults, a linter part that self-skips, a
column-subset view whose only path is the full span, a calibration validator
whose own `alpha` is untested. That is the inverse smell - a general mechanism
whose only path is the specific case - and there the right move is almost always
FINISH, not CUT.
