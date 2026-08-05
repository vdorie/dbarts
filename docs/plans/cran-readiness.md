# cran-readiness

agent: opus
rng: neutral - every fix must be behavior-neutral (equivalence
  identical); doc/example/build-flag changes only. Anything that
  cannot be fixed neutrally is REPORTED, not fixed.
budget: findings-driven; expect Rd/example/DESCRIPTION touch-ups plus
  compiled-code warning hygiene.

## Goal

Find the CRAN rejections early: a clean R CMD check --as-cran run on
the built tarball, a full ASAN+UBSAN pass over the compiled core, and
a triaged findings list - mechanical fixes landed, judgment calls
reported. Release itself waits (VD); this buys the lead time.

## Context

- R CMD check --as-cran has never run cleanly in the session harness
  (it breaks there); run it in a plain shell against a built tarball
  (R CMD build first). The season added many Rd examples
  (survivalProbabilities, samplePriorPredictive, growFromRoot,
  n.grow.sweeps) that have never seen the example-timing limits.
- CRAN's extra checks run gcc-ASAN/clang-UBSAN on Linux; this machine
  is arm64 macOS, where tests/cpp sanitizes directly and the R-loaded
  .so sanitizes with more ceremony. A sanitizer caught a real defect
  this season (the callback landing), so the pass has teeth.
- Sanitizer CI gating is a separate follow-up decision (TODO soft
  call), not this item. Manual pass only.

## Steps

1. R CMD build, then R CMD check --as-cran on the tarball in a clean
   shell. Triage every NOTE/WARNING/ERROR: fix the mechanical ones
   (Rd nits, cross-refs, example runtimes via smaller inputs or
   \donttest, DESCRIPTION fields, non-ASCII, stray files in the
   tarball via .Rbuildignore), report the judgment ones.
2. Compiled-code warning hygiene: build src/ with CRAN-like strict
   flags (-Wall -Wextra -pedantic) and clear what is clearable
   without behavior change.
3. Sanitizers: tests/cpp under -fsanitize=address,undefined
   (mandatory, full run); the R-loaded package under UBSAN and, if
   macOS tooling permits, ASAN over the full tinytest suite
   (best-effort; document exactly what ran and what remains for
   CRAN's Linux images).
4. Findings memo appended to this plan: what was fixed, what was
   found clean, what needs VD, what could not be exercised on macOS.

## Verification

- check --as-cran status recorded verbatim (target: no ERROR/WARNING,
  NOTEs enumerated and justified or fixed).
- Sanitized tests/cpp: zero reports. Sanitized R suite: zero reports
  on whatever subset ran.
- Behavior neutrality: full tinytest 2727+ / 0; equivalence vs
  equivalence-ac6ec2c.rds 22/22 IDENTICAL; tests/cpp clean rebuild.

## Findings memo (2026-07-11)

Env: macOS Tahoe 26.5, aarch64-apple-darwin23, R 4.6.1, Apple clang
21.0.0. Built with R CMD build (vignettes OK), checked with
R CMD check --as-cran --no-multiarch on the tarball. Base commit
7c0c7e1; two fix commits landed (10ee74c, 16d00c7) plus this memo.

### Final R CMD check --as-cran status (verbatim)

Re-run after the fix; only the compiled-code stderr NOTE remains:

```
* checking R code for possible problems ... OK
* checking compiled code ... NOTE
File 'dbarts/libs/dbarts.so':
  Found '___stderrp', possibly from 'stderr' (C)
    Object: 'misc.a'

Compiled code should not call entry points which might terminate R nor
write to stdout/stderr instead of to the console, nor use Fortran I/O
nor system RNGs nor [v]sprintf.

See 'Writing portable packages' in the 'Writing R Extensions' manual.
* checking examples ... OK
* checking examples with --run-donttest ... OK
* checking tests ... OK
  Running 'tinytest.R'
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
Status: 1 NOTE
```

Everything else (incoming feasibility, DESCRIPTION, top-level files,
non-ASCII, Rd cross-refs / metadata / usage / contents, example
timings incl. --run-donttest, tests, vignette rebuild, PDF+HTML
manual) checks OK. No ERROR, no WARNING. Example timings never tripped
the 5s limit, so no example was shrunk or moved to \donttest.

### Findings and dispositions

1. R NOTE "no visible global function definition for 'dbinom'"
   (pointwiseLogLikelihood, R/generics.R:54). dbinom sits beside dnorm
   (already imported); it was the one missing stats import and would
   fail in a clean session where stats is not attached.
   -> FIXED in 10ee74c: added dbinom to the stats importFrom in
   NAMESPACE. Behavior-neutral (equivalence 22/22 identical).

2. Compiled-code NOTE "Found '___stderrp', possibly from 'stderr'"
   in misc.a. Source: src/misc/io.c default hooks printToStderr /
   flushStderr, the documented R-independent fallback (io.h). In an R
   session R_init_dbarts (src/R_interface.cpp:362-363) UNCONDITIONALLY
   repoints misc_printf -> Rprintf and misc_flushOutput ->
   R_FlushConsole at DLL load, before any .Call, so no stderr write is
   reachable from R; the symbol persists only because the fallback is
   compiled in.
   -> JUSTIFIED NOTE / REPORTED to VD. Not fixed because every neutral
   removal is out of scope or non-neutral: (a) misc.a is shared -
   tests/cpp links the very same ../../src/misc.a - so gutting the
   stderr default changes the standalone/tests-cpp contract (not
   behavior-neutral for non-R consumers); (b) routing the default
   through Rprintf needs an R-only build guard + a -D macro, i.e. a
   build-system change this task is told not to make. Recommended
   durable fix for a later, deliberate pass: guard printToStderr/
   flushStderr behind a "built for R" macro so the R build defaults the
   hooks to NULL (they are set in R_init before use) while tests/cpp
   keeps the stderr fallback. This is a common, non-blocking NOTE for
   packages whose C core is intentionally host-agnostic.

3. Top-level file TODO (33 KB dev notes) ships in the tarball.
   "checking top-level files ... OK" - R does not flag it, so it is not
   a check finding. -> REPORTED (optional): if VD would rather not ship
   internal notes, add ^TODO$ to .Rbuildignore (trivial, neutral). Left
   as-is since the check is clean and it is the author's call.

### Strict-flag warning inventory (-Wall -Wextra -pedantic)

Full src/ rebuild with the sanitizer-free strict flags injected via a
scratch CPPFLAGS (not committed to the build system). Our-code
warnings, all cleared in 16d00c7, behavior-neutral:

- R_interface_bartcore.cpp:1237,1294,1470 missing field
  'ownedTreatment' initializer [-Wmissing-field-initializers]. The
  BartcoreHolder aggregate has 8 members; the three sites listed 7,
  leaving ownedTreatment implicitly value-initialized. FIX: added the
  explicit trailing {} (identical value-initialization).
- misc/moments.c:1330 unused parameter 'i' [-Wunused-parameter].
  The whole body of misc_stat_setSIMDInstructionSet is under
  #ifdef COMPILER_SUPPORTS_SSE2; on arm64 (no NEON moments kernel) i is
  unused. FIX: (void) i; (emits no code).

Remaining warnings are NOT our code and are left untouched:

- R_ext/Boolean.h:66 "enumeration types with a fixed underlying type
  are a C23 extension" [-Wc23-extensions] x6 - inside R's own installed
  header, triggered by -pedantic on C TUs that include R. Unfixable by
  us; disappears without -pedantic.
- tests/cpp/test_scan.cpp:34 'lwz2'/'rwz2' set but not used
  [-Wunused-but-set-variable] x2 - test-harness only (never compiled
  into the package), pre-existing, out of CRAN scope.

### Sanitizer coverage

- MANDATORY tests/cpp under ASAN+UBSAN: DONE, zero reports. The three
  archives (misc.a, external.a, rc.a) were rebuilt with the sanitizer
  baked into the compiler drivers (CC/CXX/CXX20 =
  clang[++] -fsanitize=address,undefined -fno-sanitize-recover=undefined
  -g) via a scratch R_MAKEVARS_USER, and test_bartcore was compiled and
  linked against those instrumented .a with the same flags. Run to
  completion: 109 test cases, "all tests passed", zero ASan/UBSan
  reports (ASAN_OPTIONS=abort_on_error=1, UBSAN_OPTIONS=halt_on_error).
  Note: macOS LeakSanitizer is unsupported, so at-exit leak detection
  did not run; heap-overflow/use-after-free/UB detection did.
- BEST-EFFORT R-loaded package under UBSAN: DONE. .so built with
  -fsanitize=undefined (dynamic UBSan runtime, so it dlopens into a
  plain framework R with no ceremony); test-load OK; full
  tinytest::test_package("dbarts") ran on the instrumented .so ->
  2727/2727 TRUE, zero "runtime error:" lines (UBSAN_OPTIONS=
  print_stacktrace=1:report_error_type=1).
- GAP - ASAN under R on macOS: NOT run. ASan must own malloc before
  main, which on an unsanitized framework R means DYLD_INSERT_LIBRARIES
  of libclang_rt.asan_osx_dynamic.dylib plus fighting R's allocator -
  genuine heroics with high false-positive risk from the unsanitized
  runtime. This path is covered by CRAN's Linux gcc-ASAN/clang-ASAN
  images (a sanitized R) and, at the core level here, by the mandatory
  tests/cpp ASAN pass above.

### Gate outputs

1. R CMD INSTALL --preclean -l <lib>: clean (exit 0).
2. tests/cpp normal rebuild + run: 109 tests, all passed.
3. tinytest::test_package: 2727 results, all TRUE, 0 fails.
4. equivalence vs equivalence-ac6ec2c.rds: 22 compared / 0 skipped,
   all "identical draws (same RNG stream)" -> behavior neutrality
   proven.
5. air/lintr/pkgdown: no .R or .Rd files were touched (only NAMESPACE
   and two compiled files), so these are out of scope / N-A.

## Revdep sweep (2026-07-22)

CLEAN: a fresh CRAN reverse-dependency query returned 24 packages, zero
dbarts-caused breaks. stan4bart in particular was verified GENUINELY
green via a real compile this time, not the earlier structSize-masked
pass (a stale build had been reusing an old dbarts_results layout and
so never exercised the current ABI).

METHOD (no revdep script is checked into the repo): a fresh
tools::package_dependencies("dbarts", reverse = TRUE, which =
c("Depends", "Imports", "LinkingTo", "Suggests")) call enumerates the
set; a fault-tolerant per-package loop then does
download.packages(..., type = "source") -> remotes::install_deps into
a throwaway temp library -> R CMD check, catching and recording
failures per package rather than aborting the sweep; a guard asserts
dbarts itself stays at 1.0-0 for the whole run (a stray reinstall
mid-sweep would silently change what is being tested). The three
MAINTAINED reverse deps - stan4bart, bartCause, treatSens - are pulled
from their COMPAT BRANCH, not the stale CRAN tarball, since that branch
is what will actually ship lockstep.

CAVEAT: remotes::install_deps(dependencies = NA) installs hard
dependencies only, NOT Suggests. A CRAN-grade sweep needs Suggests
installed too - 6 packages in this run ERRORed solely on a missing
Suggests package, none on a dbarts issue; treat those 6 as unresolved
by this method, not as dbarts breakage.

EQUIVALENCE-TRIO RESIDUAL: the 2026-07-22 engine commits (warm-starts
0a27207, sparse 343dd4c) were gated on the tinytest reproducibility
snapshots (passed) but not an explicit benchmarks/R/equivalence.R
baseline compare - that gate runs on the weekly cron, not push-gated.
Both changes are draw-neutral by construction, so the weekly run is
expected to confirm; recorded here so the gap is visible until it
does.

## Pre-submission gate run (2026-07-25)

Env: macOS Tahoe 26.5.2, aarch64-apple-darwin23, R 4.6.1, Apple clang.
Base commit 51edbf6 plus the fixes below. The 2026-07-11 memo above
covers the previous pass; this is the re-run on the finished branch.

### Fixes landed in this pass

1. `inst/CITATION` built `textVersion` from base R's `version` object
   instead of `meta$Version`, so `citation("dbarts")` printed FOURTEEN
   entries - one per element of `R.version` ("... Sampler. 6.1.",
   "... Sampler. Happy Hop."). Now a single correct entry.
2. `DESCRIPTION` `Depends: R (>= 3.5.0)` -> `R (>= 4.0.0)`. `CXX_STD =
   CXX20` is only recognised from R 4.0.0 (R 4.0.0 NEWS, PACKAGE
   INSTALLATION: "the beginnings of support for the recently approved
   C++20 standard"); the old floor promised installs on releases that
   cannot parse the standard the engine requires.
3. `DESCRIPTION` `Date` refreshed to the run date (CITATION derives its
   year from it).
4. `TODO` added to `.Rbuildignore` - the dev notes no longer ship in the
   tarball. `configure~` (a tracked autoconf backup) deleted, and `*~`
   added to `.gitignore` so it cannot come back. R CMD build already
   excluded both via its own `~$` and `config.(log|status)` rules, so
   neither had ever shipped.
5. Four unreachable functions deleted from `R/bartcore.R`:
   `bartcoreSetSigma`, `bartcoreSampleTreesFromPrior`,
   `bartcoreSampleNodeParametersFromPrior`, `bartcorePrintTrees`. The
   dbartsSampler R5 methods `.Call` the C entry points directly
   (R/dbarts.R), nothing in R/, src/, NAMESPACE, man/, tests, vignettes
   or benchmarks names them, and the package uses no dynamic dispatch
   (`do.call(paste0(...))`, `match.fun`, `get(paste0(...))`) that could
   reach them. Surfaced by the covr audit below.
6. Tests added for `makeind` and `makeTestModelMatrix` - both exported
   and documented on `makeind.Rd`, both with zero test coverage before
   this pass (appended to inst/tinytest/test-makeModelMatrix.R, which
   goes 65 -> 75 assertions).

### R CMD check --as-cran

Status: OK, ZERO NOTEs - run twice, once on a `git archive` of 51edbf6
and once on the tree with the fixes. Both clean, including the checks
the 2026-07-11 pass could not clear:

```
* checking compiled code ... OK          <- the ___stderrp NOTE is gone
* checking examples ... OK
* checking examples with --run-donttest ... OK
* checking tests ... OK
  Running 'tinytest.R' [36s/47s]
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... OK
Status: OK
```

Tarball 1.2 MB, no stray files. tinytest under check: 3402 results, all
pass. The 50+ R warnings the test suite emits are the intentional
"'test' is unnamed but 'x' had named predictors" warning (51 of them);
verified it does NOT fire on the plain `bart(x, y, x.test)` unnamed-matrix
idiom, only when the training data is genuinely named - correct, not noise.

### Equivalence

`equivalence.R compare baselines/equivalence-7903855.rds
--strict-coverage`: 27 compared / 0 skipped, EVERY scenario "identical
draws (same RNG stream)". This closes the EQUIVALENCE-TRIO RESIDUAL
recorded above: the 07-21..07-24 work (interaction-constraints P4,
warm-starts, sparse, the audit fixes, and the phase-B copy-paste dedup)
is draw-neutral, not merely statistically indistinguishable.

### SBC

`sbc.R gaussian 100 200 30`: prior moment check PASS; all seven
functionals (avg.f, f.star1-5, sigma) PASS on both chi-square and KS,
max ecdfDiff 0.079 against a 0.132 band. 351s.

### Coverage audit (the one-shot covr item)

R-level only (`options(covr.flags = list())` - gcov attribution over the
header-only C++ engine is the noise the item says to ignore), driven by
the tinytest suite. TOTAL 87.1%.

Lowest files: hooks.R 0%, guessNumCores.R 14.7%, sliceSample.R 57.0%,
plot.R 75.8%, partialDependence.R 77.0%,
updatePredictorPerObservationJointly.R 78.1%, model.R 80.9%,
rbart.R 85.5%. Everything else >= 87%.

11 of 283 functions are fully uncovered. Triaged:

- GENUINE GAPS, now fixed: `makeind`, `makeTestModelMatrix` (exported,
  documented, no test) - see fix 6.
- DEAD CODE, now deleted: the four R/bartcore.R wrappers - see fix 5.
- ATTRIBUTION NOISE, no action: `hooks.R::.onLoad`/`.onUnload` (run
  before instrumentation, so covr can never see them); `model.R::dart`
  (tests DO exercise it, via `tree.prior = dart()` in
  test-dart-mixed-columns.R and `dbartsPriors$dart(...)` in
  test-model-priors.R - the prior list holds an uninstrumented copy);
  `guessNumCores.R`'s 14.7% is the per-platform branch set, of which one
  host can only ever run one.
- LEFT AS REAL BUT LOW-VALUE GAPS: `dbarts.R::show` (S4 show method) and
  `mixedMatrix.R::dimnames.dbartsMixedMatrix`. Both are display/accessor
  surface, untested; flagged rather than fixed.

sliceSample.R at 57% is the largest remaining genuine R-level gap (it is
exercised only through the hyperprior paths that use it); worth a look if
a later pass wants the number up, not a release blocker.

### Not run, and why

rchk and valgrind still have never run: macOS arm64 supports neither, and
the workflows carrying them (plus equivalence, sbc, revdep-smoke) 404 on
GitHub because both `schedule` and `workflow_dispatch` bind to the default
branch and these files exist only on bartcore. VD held the merge on
2026-07-25 pending their own review of the branch; revisit the checks
then. win-builder was deliberately skipped this pass.

### Revdep sweep WITH Suggests (2026-07-25)

Closes the caveat on the 07-22 pass (which installed hard dependencies
only, leaving 6 packages ERRORing on a missing Suggests). Method as
before, plus each package's own Suggests installed into the sweep
library first; dbarts pinned at 1.0-0 and asserted between packages;
stan4bart and bartCause built from their local compat branches. 24
reverse dependencies (up from 23 - `lorax` is new, `treatSens` no longer
declares dbarts).

21 OK. Notably stan4bart OK on a real compile, tidytreatment OK, and
bartXViz OK - the last was "not checkable locally" in every previous
pass. Three results need a decision:

1. `lorax` - 2 ERRORs, GENUINE 1.0-0 BREAK. It calls
   `dbarts::bart()` with a 3-level factor response. 0.9-x coerced any
   factor response with `as.double(as.integer(data) - 1L)` (old
   R/data.R:255), i.e. silently fit levels as the numbers 0/1/2; 1.0-0's
   `resolveClassificationFamily` refuses and points at
   `bart2(family = "multinomial")`. The refusal is right - the old
   behavior was a silent miscoding - but lorax's examples and tests fail
   against this release.
2. `insight` - 1 ERROR, GENUINE 1.0-0 BREAK. Its tests call
   `dbartsData(bill_len ~ ., penguins)`, and palmerpenguins has NA
   responses. 0.9-x set `na.action = stats::na.omit` on the model frame
   (old R/data.R:202) and silently dropped those rows; 1.0-0 passes NAs
   through (R/data.R:766, deliberate, comment records the change) and
   then errors "response contains missing values".
3. `bartCause` - 1 NOTE, not a dbarts issue: the tarball built from the
   local working copy carried `.claude/`, which its `.Rbuildignore` did
   not exclude (dbarts excludes it via `^\.claude$`). RESOLVED the same
   day by VD (bartCause 695c603); re-checked against that commit,
   Status OK. With bartCause resolved the sweep stands at 22 of 24 OK;
   lorax and insight are the only non-OK results.

Both breaks are the same class - 1.0-0 refuses input that 0.9-x silently
mishandled - and both are consistent with the NEWS UPGRADING section. VD
decides whether to let them break or soften; softening (2) is the more
arguable of the two, since dropping NA-response rows was at least a
defined behavior. Neither is a defect in dbarts.

### Defects found by the coverage work

Writing tests for the uncovered branches turned up two real bugs, both
fixed in this pass:

1. `R/sliceSample.R` - the slice loop's exhaustion check tested the loop
   counter (`if (j == maxIter) stop(...)`) rather than whether a draw was
   accepted, so a draw accepted on the maxIter'th shrink pass raised
   "slice sampler failed: maxIter reached" over a perfectly good sample.
   Now tracks acceptance. RNG-neutral: the draws inside the loop are
   unchanged, only the post-loop test differs. Reachable from
   rbart_vi's custom-prior tau updates, the sampler's only caller.
2. `R/plot.R` and `R/generics.R` (three sites) - the guards that name a
   missing `keeptrees` / `keeptrainfits` compared
   `as.character(call[[1L]])` against "bart2". For a namespace-qualified
   fit (`dbarts::bart(...)`, the idiomatic form inside other packages)
   that expression is `c("::", "dbarts", "bart")`, so the comparison is a
   length-3 condition - an error since R 4.2. Users who forgot
   `keeptrees` got "the condition has length > 1" from `predict`,
   `extract`, or `plot` instead of the message telling them what to do.
   Now routed through a `callName()` helper in R/utility.R that strips
   the qualifier. The existing test missed it because it asserted
   `expect_error()` with no pattern; that assertion is now pinned to the
   message.

### Coverage added

- `inst/tinytest/test-slice-sample.R` 11 -> 22: the maxIter-acceptance
  regression, two finite boundaries, the low-density rejection restart
  and the natural-scale branch that cannot do it, the non-finite-start
  mode-finding failure, boundary rejection in `rejectionSample`, and
  envelope exhaustion on both scales.
- `inst/tinytest/test-plot-generics.R` 19 -> 24: multi-chain
  `plot.bart` and `plot.rbart` (the matrix-shaped sigma trace, 28 lines
  that had never executed), the binary grouped fit's
  probability-interval branch, and the two namespace-qualified error
  paths.
- `inst/tinytest/test-generics-errors.R`: `predict` and `extract` guards
  pinned to their messages for both qualified `bart` and `bart2` fits.
- `inst/tinytest/test-makeModelMatrix.R` 65 -> 75: `makeind` and
  `makeTestModelMatrix`.

Measured after the fact, same covr configuration: TOTAL 87.09% ->
88.23%, sliceSample.R 57.0 -> 80.5, plot.R 75.8 -> 83.8, generics.R 88.0
-> 88.9, bartcore.R 89.0 -> 91.2 (the dead wrappers leaving). Fully
uncovered functions 11 -> 5, the survivors being the two covr artifacts
(`.onLoad`, `.onUnload`), the `dart` misattribution, and the two real
but low-value ones (`dbarts.R::show`,
`mixedMatrix.R::dimnames.dbartsMixedMatrix`). Suite runtime 3402 -> 3438
assertions, no meaningful change in wall clock.

Not covered, deliberately: the slice sampler's underflow guard
(`getInterval` evaluates the target first and dies with R's own
"missing value where TRUE/FALSE needed", so the guard is unreachable by
a NaN target - a defensive-path robustness gap, not worth contorting a
test around at release time), and the hand-rolled gradient-ascent
fallback in `findMode`, which is exercised only up to its own failure.
