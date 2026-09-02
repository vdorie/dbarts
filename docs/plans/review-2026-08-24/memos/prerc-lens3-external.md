# Pre-RC audit, lens 3: external eyes (CRAN reviewer / first-time user / LinkingTo consumer)

Repo `/Users/vdorie/Repositories/dbarts`, branch `bartcore`, HEAD **7a8c7286** (verified, clean). Evidence from a tarball
(`R CMD build`, 54 s) installed into a private lib; R 4.6.1, aarch64-apple-darwin23. Read-only: no edits, no write-side git,
private worktrees untouched.

**Verdict: no CRAN-blocking defect.** Every `R CMD check` leg run here is `Status: OK`. But two *silent wrong-answer* traps sit
on the main user path, and the shipped consumer header omits the one fact every `LinkingTo` consumer needs first.

## A. User-facing correctness - highest severity found

**A1. A one-sided formula fits silently on an all-zero response - pre-RC, ~6 lines + test.** `bart2(~ x1 + x2, df)` returns a
`bart` object with no error, warning or message: `all(fit$y == 0)` is TRUE, `mean(fit$sigma)` is 1.9e-09. Cause
`[[R/data.R:1164-1166@7a8c7286]]`, `if (is.null(y)) { y <- rep(0, NROW(modelFrame)) }` - right for `dbartsData()` inside a composed sampler,
wrong for `bart`/`bart2`, which fit and return immediately. A typo yields plausible-looking output. Fix: refuse a response-free
formula on those entries only.

**A2. `family = "ordinal"` has no response-type guard - pre-RC, ~5 lines + test.** `bart2(ycont ~ x1, df, family = "ordinal")`
on 60 continuous values silently builds a **60-category** ordinal fit (levels `"-3.73309396896429"`, ...). `[[R/data.R:648-650@7a8c7286]]`
only refuses `K < 2`; every sibling family guards its response (`[[R/spec.R:188@7a8c7286]]` nbinom, `[[R/spec.R:206@7a8c7286]]` probit).

## B. LinkingTo consumer

**B1. `LinkingTo: dbarts` alone does not work, and the header never says so - pre-RC, ~6 lines.** I built a real minimal
consumer (C stub TU + C++ prototype TU). It **compiles clean**, but the runtime handshake worked only on the third try:
- `LinkingTo: dbarts` only -> `Error: function 'dbarts_apiHash' not provided by package 'dbarts'`
- plus `Imports: dbarts` **in DESCRIPTION** -> *same error*; `isNamespaceLoaded("dbarts")` is `FALSE`
- `importFrom(dbarts, dbarts)` in **NAMESPACE** -> works (`api major/minor: 1 0`, `structSize: 96`)

`R_GetCCallable` needs dbarts's DLL loaded, and R loads a namespace from NAMESPACE import directives, not the DESCRIPTION
`Imports` field. `inst/include/dbarts/dbarts.h` documents stubs, versioning and the ABI hash exhaustively, but
`Imports`/`Depends`/`loadNamespace`/`NAMESPACE` appear nowhere in it (only "Imports Rinternals.h", line 73). The error names
neither cause nor cure. Fix: Doxygen near lines 8-16. Pre-RC because the header is the only shipped consumer doc and the
surface freezes at 1.0-0.

**B2. `Depends: R (>= 4.0.0)` cannot be honored by a C++20-concepts package - pre-RC, 1 line.** `[[src/bartcore/model.hpp:40-145@7a8c7286]]`
defines 8 `concept`s; concepts need gcc >= 10 / clang >= 10, but the reference Windows toolchain for R 4.0/4.1 is Rtools40,
gcc **8.3.0**. Rtools42 (R 4.2, gcc 10.3) is the first that can build this. CRAN will not catch it (r-devel/release/oldrel
only), so it is a truthfulness issue. Raise to `R (>= 4.2.0)`, or `4.3.0` for margin.

**B3. Compat branches - all three bartCause failures are the consumer's.**
- **bartCause** (`dbarts-1.0`, `7ae6e83`, clean/pushed): the `alignForestBasisToSubset` failures at
  `[[tests/testthat/test-04-bartc.R:301@7a8c7286]]`, `[[test-14-bcf.R:187@7a8c7286]]`, `[[test-15-subset-pscore.R:38@7a8c7286]]` are **consumer's to fix**. dbarts
  deliberately inverted the basis/subset contract *after* bartCause's tip - `4aa52288` (08-20), `47cdb96a` (08-24), both in
  NEWS.Rd - and `[[R/model.R:955-978@47cdb96a]]` now refuses a basis whose row count matches the subset rather than the full data;
  `bartCause [[R/bcf.R:144@47cdb96a]]` still pre-subsets on the formula path. **Two edits the recorded count misses**:
  `[[test-14-bcf.R:220-224@47cdb96a]]` asserts the *old* direction, and the oracle at `[[test-14-bcf.R:16-17@47cdb96a]]` pre-subsets its basis.
  **`[[TODO:306-309@47cdb96a]]` still calls bartCause "ready and PUSHED ... Status OK"** - stale, contradicted by
  `[[docs/plans/bartcore-review-tour.md:463-467@47cdb96a]]`. Correcting it is pre-RC.
- **treatSens** (`1d7c697`): no recorded failure; needs only a `--preclean` rebuild (compile-time hash predates the re-bake).
- **stan4bart** (`bartcore`, `2185866`): **no signature mismatch** - all 28 `dbarts_*` symbols match, every `dbarts_results`
  field it writes is above the 1.0-0 boundary. The live risk closed itself: dbarts `c44fcbc5` added `numThreads` to
  `dbarts_sampler_predict` and re-baked the hash (08-25 03:21); stan4bart `2185866` updated its call site (04:19).
- All three still owe the `[[TODO:302-305@c44fcbc5]]` re-verify against the final re-baked header.

## C. CRAN reviewer

**C1. `Date: 2026-07-25` starts emitting a NOTE on 2026-08-26 - at-RC, 1 line.** `tools:::.check_package_CRAN_incoming` has
`if (dd < Sys.Date() - 31) "The Date field is over a month old."` Today `Sys.Date() - 31` is exactly 2026-07-25, so it is false
by one day; tomorrow it fires, adding a second NOTE beside the known "Days since last update". Set `Date` at submission.
(Dropping the field instead needs a companion edit - `inst/CITATION` derives `year` from `meta$Date`.)

**C2. Installed size 6.1Mb is an INFO, not a NOTE, on current R - post-1.0.** Check reports `installed size is 6.1Mb` with
`libs 1.9Mb` and `tinytest 1.7Mb` as the sub-directories over 1Mb. `tools:::.check_packages` branches
`if (R_check_use_log_info) infoLog else noteLog`; `_R_CHECK_LOG_USE_INFO_` defaults TRUE and `--as-cran` forces it TRUE, so
status is unaffected on R >= 4.5 (it would be a NOTE on older R). CRAN policy still asks for < 5 Mb, so a human may raise it.
Breakdown: libs 1.9 (`dbarts.so` 1.97 MB, `-g -O2`), tinytest 1.7 (**169 test scripts**, no large fixtures), doc 0.9 (3 PDFs),
R 0.8, help 0.5.

**C3. `configure.ac` selects the default-standard compiler, not the C++20 one - post-1.0, ~8 lines + `autoreconf -i`.** WRE
calls this "essential": use `R CMD config CXX20`/`CXX20STD`/`CXX20FLAGS`; `[[configure.ac:29-34@c44fcbc5]]` uses plain `CXX`/`CXXFLAGS`.
**Impact today is nil** - on R 4.6.1 `R CMD config CXX` already returns `clang++ -arch arm64 -std=gnu++20` (C++20 became R's
default in 4.6.0), and the only C++ probes are `AC_C_RESTRICT` (43-45) and `AC_FUNC_STRERROR_R` (66-68), neither
standard-sensitive. Latent on R 4.3-4.5, and a trap the first time a C++20-sensitive probe is added.

## D. First-time-user rough edges

- **`citation("dbarts")` prints two near-duplicate entries - pre-RC, 2-3 lines.** `inst/CITATION:1` is a bare
  `citation(auto = meta)` and 8-21 add a hand-written `bibentry`; both are collected, both with empty BibTeX keys and
  *conflicting* authors (auto: all three + github.io; manual: Dorie only + CRAN). Header reads "To cite the 'dbarts' in
  publications use:" (stray "the") and prints *after* the first entry, which has none. Neither cites Chipman/George/McCulloch
  (2010), whose doi is already in DESCRIPTION.
- **NEWS.Rd 1.0-0 is 1992 of 2551 lines** (78%; all 20 prior releases = 555). Renders clean (`checkRd`, `Rd2txt`, news DB
  silent). Flat: 4 `\itemize`, 200 items, no nesting - UPGRADING 11 (5-77), NEW FEATURES 68 (78-955), C API 9 (956-1023), BUG
  FIXES 112 (1024-1995). Worst passages: **131-195**, one 65-line bullet carrying ~15 changes, with the `sd` renormalization
  that *moves Gaussian draws* buried at 156-158; **1018-1023**, the only mandatory consumer action, last in its section, in
  internal vocabulary ("re-baked", "the lockstep handshake as designed"), absent from UPGRADING; **1024-1995**, 112
  undifferentiated items of which **>=13 are not bug fixes** (1518 is a 28-line feature; 1545 `xbart` loses `control` is a
  breaking removal). UPGRADING misses >=6 breaking changes documented in the body (1545, 1542-1544, 1794, 1902,
  1615/1625/1660/1710, 1018-1023); the `R_C_interface.hpp` removal is called out (17-20, 958-960) but with no call mapping, and
  `LinkingTo` never appears there. `news(package=)` drops the lead-in at 6-7, so CRAN's news view shows UPGRADING with no
  statement of what it is. **Minimum fix: repair UPGRADING, ~40 lines, pre-RC.** Full restructure to ~360 lines with overflow
  moved to a migration vignette is ~400 lines + ~250 of prose, roughly a day - the largest discretionary item in this report.
- **README never shows how to fit a model**, and documents priors that are not reachable - pre-RC (surface). `[[README.md:12-27@c44fcbc5]]`
  is a bare feature list; the only R is `install.packages`/`install_github` (38-46), and `bart()`/`bart2()`/`xbart()`/
  `rbart_vi()` are never named (line 27 points at `inst/NEWS.Rd`, a source path). Worse, `resid.dist = student(3)`,
  `tree.prior = dart()`, `node.prior = linear(...)`/`gp(...)` (20-21) work only typed literally inline (resolved by
  `substitute`, `[[R/model.R:1608-1633@c44fcbc5]]`): `p <- dart()` gives `could not find function "dart"`, and `dbartsPriors` carries
  `cgm dart normal linear gp chisq fixed chi` - **not `student`** - so no non-inline route to Student-t errors exists.
- **5 of the 10 advertised families have no runnable example anywhere**: multinomial, ordinal, hazard/hazard.probit/
  hazard.logistic, hurdle.lognormal/twopart. probit and nbinom have one only via the low-level `dbarts()`, never
  `bart()`/`bart2()`. Live: gaussian (many), logistic (`man/bart2.Rd`), aft (`man/survivalProbabilities.Rd`), variance forest,
  BCF (`man/forest.Rd`), grouped (`man/rbart.Rd`). Relatedly, **no vignette covers the new families** - the three vignettes
  call `dbarts()` 17x and `bart2()` 2x, all gaussian, so the headline 1.0 feature has no narrative doc (the natural home for
  the NEWS overflow). `man/dbarts-package.Rd` has no `\examples`.
- **Objects that dump raw lists**: `pdbart()` prints **1607 lines** of raw `$fd` matrices - no `print.pdbart`; same for
  `forest()` (32) and `interactions()`/`blocks()`, while sibling `varianceForest()` *does* have `print.dbartsVarianceForest`
  (`NAMESPACE:48-49`). `dbartsData()` prints **135 lines including the entire y vector and x matrix** via default S4 `show`.
  `dbartsSampler`'s `show` prints only `dbarts sampler` + the call - no n.obs, n.predictors, n.trees or family.
- **A bad family token gives a bare `match.arg`**: `'arg' should be one of "auto", "gaussian", ...` - never names the argument,
  function or offending token; `family = binomial` (the glm habit, unquoted) gives `'arg' must be NULL or a character vector`.
  Worse at narrow surfaces: `xbart(family = "nbinom")` and `rbart_vi(family = "logistic")` tell a user a *real* family is not a
  family, with no hint that `bart2` fits it. `[[R/bart.R:2621-2640@c44fcbc5]]` (`bartRedirectedFamilies`) already does this right, but only
  for `bart()`. By contrast the response/family *mismatch* messages are best-in-class and actionable (`[[R/data.R:571-614@c44fcbc5]]`).
  Separately, the formula-refusal message (`[[R/bart.R:2266@c44fcbc5]]`; also `[[R/dbarts.R:495@c44fcbc5]], [[R/dbarts.R:553@c44fcbc5]], [[R/dbarts.R:587@c44fcbc5]]`) says "use the matrix interface -
  `bart2(x.train, y.train, ...)`", but `bart2`'s formals are `formula, data`, so the working call deparses back as
  `bart2(formula = x, data = yh, ...)` in the fit's own `Call:` block.
- **`_pkgdown.yml` coverage is complete** (all 29 topics grouped, every export documented) **but the grouping misleads**:
  `xbart`/`pdbart` sit under *Model fitting* (15-20); `interactions`, `blocks`, `forest`, `varianceForest`, `dbartsPriors`,
  `sparseFactor` are buried under *The mutable sampler* (21-40) though they are `bart2()` arguments; there is **no
  prediction/extraction group** at all (`predict`/`extract`/`fitted`/`residuals` exist only as aliases in `bart.Rd`);
  `dbarts-package.Rd` and `dbarts.Rd` carry *identical* titles; *Survival* (42) omits discrete-time hazard.
- **Minor**: `xbart` supports 4 of the 13 family tokens, so the new families cannot be cross-validated (capability gap, not a
  message problem). `Biarch: yes` is obsolete (32-bit Windows went in R 4.2). `[[gibbs_sampler_mixture_model.Rmd:398@c44fcbc5]]` links
  `blob/master/...` while the default branch is `main` (200 only via GitHub's rename redirect).

## Checked and clean - no action

- **All `R CMD check` legs run here are `Status: OK`**: structural, examples **including `--run-donttest`** (each under 5 s;
  only `man/pdbart.Rd` and `man/summary.bart.Rd` use `\donttest`), vignette rebuild **11 s**, size/FF/registration. Build 54 s,
  install 31 s. **All 29 Rd files have `\value`** - the most common CRAN-reviewer rejection, avoided.
- **The C++20 declaration is exactly right.** WRE (R 4.6.1) says a package "should" put the standard in `SystemRequirements`
  *and* "If it has a Makevars file ... this should include the line `CXX_STD = CXX20`". dbarts has both, in
  `[[src/Makevars.in:33@c44fcbc5]]` **and** `[[src/Makevars.win:23@c44fcbc5]]`; they agree, and check reports `specified C++20` with no mismatch. The
  premise that current advice is "SystemRequirements only, no CXX_STD" is outdated - that was C++11/14-era guidance.
- **Portability + Windows + registration**: zero `#pragma` in `src/`; 5 `__builtin_*`, all behind configure/compiler guards
  (`[[misc/partition.c:28@c44fcbc5]]`, `[[include/misc/intrinsic.h:36@c44fcbc5]], [[include/misc/intrinsic.h:37@c44fcbc5]], [[include/misc/intrinsic.h:60@c44fcbc5]]`, `[[include/misc/alloca.h:19@c44fcbc5]]`); no `-O3`, `-march`, `-mtune`,
  `-flto`, `-ffast-math`, `-funroll` anywhere; no OpenMP; Windows pthreads via `$(SHLIB_PTHREAD_FLAGS)`.
  `tools/check-win-drift.R` OK (6 version literals, 4 config headers, 38 expected-absent entries, 0 UNKNOWN) - it checks
  `PACKAGE_*` literals in the `*.win` config headers against DESCRIPTION and each `*.win` macro nameset against its `*.in`
  template, against a table of deliberate omissions each with a recorded reason; `tools/check-rc-codoc.R` OK (42 methods).
  `R_registerRoutines` + `R_useDynamicSymbols(FALSE)` (`[[src/R_interface.cpp:285-287@c44fcbc5]]`), the flat C API registered by expanding
  `DBARTS_C_API_LIST` so the table cannot drift from the header, and `tools::checkFF(registration = TRUE)` returns empty.
- **Suggests all conditional** - every `posterior`/`Matrix`/`survival` use is behind `requireNamespace`, and `R/hooks.R`
  registers the `posterior` S3 methods correctly for **both** load orders (`setHook(packageEvent("posterior", "onLoad"), ...)`
  plus an `isNamespaceLoaded` catch-up). No `BART` dependency, declared or used. No `.onAttach`, so no startup message.
- **Metadata and hygiene**: all 11 distinct URLs in DESCRIPTION/Rd/vignettes/CITATION/README return **200** (`urlchecker` not
  installed; checked with `curl -L`). Non-ASCII is one char, `[[inst/NEWS.Rd:2080@c44fcbc5]]`, under a declared `\encoding{UTF-8}` (line 3);
  `R/`, `src/`, `man/`, `vignettes/` are pure ASCII. No `.DS_Store` or stray dotfiles in the tarball (`.Rbuildignore` catches
  the one in `inst/include/`); installed `include/` holds only `dbarts.h`. No `data/`, so `LazyData` correctly absent;
  `Authors@R` parses with ORCIDs; `License: GPL (>= 2)` needs no LICENSE file. **Every DESCRIPTION feature claim is backed by
  code** (discrete-time hazard, GP leaves, monotonicity, interaction constraints, DART, grouped random intercepts,
  missing-value incorporation, sparse predictors, Student-t, variance forests, two-part).
- **Consumer compile**: the minimal `LinkingTo` package built clean with `-Wall` in **both** C (stub path) and C++ (prototype
  path) - zero warnings. `dbarts.h` is C-clean (`extern "C"`, `stdint`, `R_NO_REMAP`-safe with a documented include-order
  caveat), so a pure-C consumer needs no C++20 declaration of its own. Its Doxygen is unusually thorough - every entry
  documented, plus per-setter pointer-retention contracts, thread rules, and the `structSize` forward-compatibility scheme.

Known-green items were not re-run: `--as-cran` (1 NOTE, "Days since last update" - a submission-cadence note, not a package
defect), tinytest, codoc, pkgdown, lintr, air, the six push workflows.
