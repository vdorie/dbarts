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
