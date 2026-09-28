# rbart-vi-port: rbart_vi back for one release as 0.9-34's R loop

agent: sonnet (R port, tests, manual, NEWS; no engine, bridge or header change)
rng: neutral (no existing draw moves; rbart_vi's own draws are new against both 0.9-34 and the tombstone)
budget: ~3000 added lines, nearly all ported: R ~1850 (the 0.9-34 loop, its slice sampler and six methods after
air formatting), tests ~900, man ~200, NEWS and docs ~150; ~120 removed (stubs, stub tests, stub Rd). Two slices.

## Goal

`rbart_vi` and its `predict`, `extract`, `fitted`, `residuals`, `plot` and `print` methods run again in 1.0-0,
as 0.9-34's pure-R Gibbs loop over the plain sampler, warn once per session that they are deprecated in favour
of `stan4bart::stan4bart`, and are removed in 1.1-0 with the rest of the tombstone registry (dec-B130).

## Context

- Ruling: dec-B130 (maintainer, 2026-09-28) replaces the tombstone half of dec-A01 and the R-only rejection in
  dec-B105 for 1.0-0. No engine work, no change to the shipped C header.
- Removal record: [The decision](../design/retire-grouped-random-effects.md#the-decision). Its "Keep the R loop,
  drop the engine path" alternative is what this plan builds.
- Source, main at 0.9-34: `rbart_vi`, `rbart_vi_run`, `rbart_vi_fit`, `packageRbartResults` and `rbart.priors`
  in [R/rbart.R:1-537](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/rbart.R#L1-L537)
  (lines 538-659 are two `if (FALSE)` stubs, not ported); `sliceSample` and `rejectionSample` in
  [R/sliceSample.R:1-179](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/sliceSample.R#L1-L179);
  the methods in [R/generics.R:143-406](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/generics.R#L143-L406),
  [R/generics.R:469-473](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/generics.R#L469-L473)
  and [R/plot.R:47-93](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/plot.R#L47-L93);
  the manual [man/rbart.Rd:1-175](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/man/rbart.Rd#L1-L175);
  nine `test-rbart-*.R` files, 769 lines.
- The bartcore-era engine version (removed by 1e5f80b2) is not ported. Two things are taken from it: the fix to
  the loop's first draw (adaptation 7) and the n.burn = 0 handling (adaptation 6).
- Today: seven error stubs in [tombstones.R](../../R/tombstones.R) with their rows in
  [`dbartsTombstones`](../../R/tombstones.R), the stub checks in
  [test-tombstones.R](../../inst/tinytest/test-tombstones.R), the aliases in
  [dbarts-deprecated.Rd](../../man/dbarts-deprecated.Rd), and two 1.0-0 NEWS items that call rbart_vi removed.

## Evidence: spike

A private build of this tree with the port applied (0.9-34's files, the adaptations below, stubs deleted), next to
a private build of main at 0.9-34, both installed outside the user library. Results:

- The ported loop runs end to end, gaussian and probit, one and several chains, serial and PSOCK, with a callback,
  the gamma prior, a dbartsData object as formula, and n.burn thinned to zero. Every method runs.
- Six seeds per family, n = 2000, 20 groups, tau = 1, two chains at the default 1500/1500/thin 5. Posterior mean
  tau, 0.9-34 against the port: gaussian 1.34/1.23, 1.06/1.08, 0.86/0.85, 0.80/0.76, 1.00/1.03, 1.38/1.08;
  probit 1.24/1.15, 1.05/1.08, 0.82/0.86, 0.78/0.76, 0.99/1.04, 1.07/1.07. Centered ranef RMSE agrees
  within 0.02 on every cell; sigma within 0.02. The port is about 30 percent faster (2.8 s against 4.1 s per fit).
- main's nine rbart test files against the port: 64 of 65 expectations pass unchanged. The one failure is
  main's hard-coded ranef snapshot (expected: draws differ). A new-level correlation check in
  test-rbart-groupby.R failed in one earlier spike run. It is seed-fragile on 0.9-34 too: below its 0.9 bar on
  3 of 12 seeds there, and on 0 of 12 on the port.
- Nothing needed a bridge or engine change. Every adaptation is R-side.

## Sampler calls the loop makes

Each is supported by the 1.0 sampler as-is unless the row says otherwise.

| 0.9-34 call or access | 1.0 | adaptation |
|---|---|---|
| `sampler$setOffset(offset, isWarmup)`, second argument `updateScale` | [`dbartsSampler$setOffset`](../../R/dbarts.R) keeps `(offset, updateScale)`; [`bartcoreSamplerSetOffset`](../../R/bartcore.R) refuses `updateScale = TRUE` only on amplitude-glued samplers, which rbart never builds | none |
| `sampler$run(0L, 1L)`; reads `$train`, `$test`, `$sigma`, `$k`, `$varcount` | same fields; `$train` still carries the offset; saved trees accumulate across calls | none |
| `sampler$getLatents(state$y.st)`, write in place | [`bartcore_getLatents`](../../src/R_interface_bartcore.cpp) fills a supplied numeric of length n | none |
| `sampler$predict(sampler$data@x)` right after `sampleTreesFromPrior` | refuses with keepTrees on and no saved draw (n.burn thins to 0) | 6, 7 |
| `sampler$sampleTreesFromPrior()`, `$setControl` | unchanged | none |
| `sampler$startThreads()`, `$stopThreads()` | no-op tombstones that warn ([`noOpThreadMethod`](../../R/tombstones.R)) | 5 |
| `model@node.hyperprior`, class `dbartsChiHyperprior` | slot keeps node (dec-A95) | none |
| `data@sigma` | slot unchanged, filled at creation | none |
| `control@keepTrees`, `n.thin`, `n.chains`, `n.threads`, `updateState` | unchanged | none |
| `.Call(C_dbarts_assignInPlace, ...)` | [`assignInPlace`](../../src/R_interface.cpp) still registered | none |
| `.Call(C_rbart_fitted, ...)` in `fitted.rbart` | removed from the bridge | 8 |
| `redirectCall` into `dbartsControl` | `dbartsControl` now has `seed` and `...` | 1 |
| `dbarts(node.prior =, resid.prior =, sigma =)` | all three are warning tombstones on `dbarts()` | 2, 3 |
| `makeCluster`, `clusterExport(rbart_vi_fit, rbart_vi_run)`, `clusterMap` | unchanged; workers load the installed namespace, so the slice sampler resolves there | none |
| `combineChains`, `convertSamplesFromDbartsToBart`, `sampleFromPPD` | combined order is chain-major now; [`sampleFromPPD`](../../R/generics.R) takes `n.chains` | 9 |

## Adaptations

Numbered for the table and the calls below; each one is small.

1. Drop `seed` and `...` from the forwarded control call. `seed` is a
   [`dbartsControl`](../../R/dbarts.R) formal now, and passing it would give every serial chain the same engine
   seed; the loop keeps 0.9-34's own `set.seed` handling. `...` would otherwise be copied in as an empty argument,
   which the spike hit as an error on every call.
2. Build the sampler as [`bartBT`](../../R/bart.R) does: `leaf.prior = normal(k)` only when `k` is supplied, the
   residual prior `chisq(sigdf, sigquant)` on an `"auto"` family through [`withResidPrior`](../../R/family.R),
   and `sigest =`. No tombstone warning fires.
3. Build the data with `factors = "indicators"` and `na.action = stats::na.omit`, 0.9-34's dummy expansion and
   row rule, as `bartBT` does.
4. Deprecation: the first statement of `rbart_vi` calls [`warnOnce`](../../R/utility.R) under key
   `tombstone.rbart_vi`: "'rbart_vi' is deprecated and is removed in dbarts 1.1-0; grouped random effects live in
   stan4bart (stan4bart::stan4bart), whose group-spread prior differs, so results move." The methods do not warn.
5. Delete the `$startThreads()`/`$stopThreads()` pair.
6. With n.burn thinned to zero and keepTrees on, take the initial fit with keepTrees off, then turn it back on.
7. Take the initial training fit from one sweep, `sampler$run(0L, 1L)$train`, not from `predict` on
   prior-drawn trees. 0.9-34 bug: a prior forest predicts the response midpoint, and on a skewed response the gap
   to the mean goes into the intercepts and stays. Spike, normal predictors, n = 1000, 10 groups, true tau 1:
   0.9-34 gives mean tau 105, 128 and 107 over three seeds, the port with the fix gives 1.56, 1.10 and 0.78. This
   is the bartcore-era fix, recorded in the archived plan rbart-custom-prior-divergence.md.
8. `fitted.rbart` averages in R rather than calling the removed C helper: `colMeans` of yhat plus the matched
   ranef columns, `pnorm` first for a binary fit. An unmatched group gives NA, not the out-of-bounds read the C
   helper did.
9. The methods keep 0.9-34's shapes: a private copy of 0.9-34's `combineOrUncombineChains` (a one-chain fit gets
   no chain margin) instead of [`combineOrUncombineChains`](../../R/generics.R) (dec-A79 adds one).
   `sampleFromPPD` gets `n.chains`. `browser()` becomes `stop()`.
10. Small 0.9-34 bugs the spike hit, fixed in R:
    - `group.by` (and `group.by.test`) named as a symbol that is not a column of `data` silently took the first
      column
      ([R/rbart.R:84](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/rbart.R#L84):
      `which.max` of all-FALSE is 1). On 0.9-34 the spike fit 200 groups where 5 were meant. Fix: look the name up in
      `data` only when it is there.
    - `predict` with a new level, several chains and `combineChains = TRUE` stops with "subscript out of bounds"
      on 0.9-34 as well: it tests the stored ranef's dimension rather than the combined one
      ([R/generics.R:227](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/generics.R#L227)).
      Fix that test, and draw the new-level ranef with tau in the same layout as the ranef it is bound to.
    - A saved fit read back into a new session: 0.9-34's `predict` silently gives different values; the 1.0
      sampler refuses because the fit stored no state. Fix: `$storeState()` on each kept sampler at the end of
      `rbart_vi_fit` when `updateState` is on. With it, the reloaded fit predicts identically in the spike.
11. Refuse any family other than gaussian and probit after the sampler is built, naming rbart_vi. `"auto"` would
    otherwise resolve a Surv response to aft. A 3-level factor is already refused by `dbarts()`.
12. `importFrom(stats, dcauchy, dgamma)`, since the tau priors use them; the graphics imports are already there.

Kept from 0.9-34 as they are: the warmup scale re-anchoring (`setOffset(updateScale = TRUE)` during burn-in), the
unweighted group mean in the ranef draw, the formals and their defaults (`seed = NA_integer_`,
`combineChains = FALSE`, `k = 2.0` shown but applied only when supplied), the "verbose output disabled" warning,
the result fields, and a `k` that is evaluated rather than parsed (`chi()` does not resolve in either build).

## Consumers

**bartCause, released (main = CRAN 1.0-10).** Calls `rbart_vi` from `getBartTreatmentFit` and
`getBartResponseFit` when `group.by` is set and `use.ranef` is TRUE (the default). Arguments: the formula (or a
`dbartsData` object as `formula` with `data` dropped), `data`, `subset`, `weights`, `group.by` (a symbol, or a
literal factor), `group.by.test`, `verbose = FALSE`, `n.chains` (10 when unset), `seed`, and any user argument that
passes a filter on the formals of rbart_vi and dbartsControl. Methods and fields it reads: `extract(fit, sample =,
combineChains =)`; `predict(fit, newdata, group.by, combineChains = FALSE)` with `group.by` third by position;
`fit$yhat.train` dims (for n.chains), `$sigma`, `$first.sigma`, `$fit` (non-NULL) and `$y`. Against the port, its
`rbart_vi fit matches manual call` treatment test passes (2 of 2, n.burn = 3 thinning to 0), with the
deprecation warning. But bartCause's own suite on 1.0 errors in 27 tests, rbart and non-rbart alike, before any fit:
`responseData@x.test[, treatmentName] <-` fails because `dbartsData@x` is a `dbartsMixedMatrix` in 1.0. So
released `bartc()` with `group.by` does not run on 1.0 with or without this port. Only a direct treatment-model
call does. The dbarts-1.0 branch has no rbart call left and needs nothing. See call 12.

**stan4bart, released (main).** Its at_home() test `nonlinearities are estimated well` calls
`rbart_vi(y ~ . - g.2, df.train, test = df.test, group.by = g.2, group.by.test = df.test$g.2, verbose = FALSE,
n.samples = 1000, n.burn = 1000)` and `fitted(fit, sample = "test")`. Replayed on four seeds, the test-set deviance
is 1.44, 1.39, 1.40, 1.41 on 0.9-34 and 1.42, 1.44, 1.39, 1.43 on the port. The bartcore branch already dropped that
comparator and needs nothing.

## Constraints

- No change under src/ or inst/include/. No change to any existing fit's draws. No rbart equivalence scenario.
- Surface: 0.9-34's argument list and six methods plus print, and nothing more. The bartcore-era extras
  (`plotTree`, `summary`, `survivalProbabilities`, `as_draws_*`, the family and factors arguments) stay out.
- The port stays one file, so the 1.1-0 removal deletes it with R/tombstones.R.
- Out of scope: bartCause and stan4bart edits, the vignette, and any mixing or prior work on rbart.

## Steps

Slice 1, code and tests (~2800 lines):

1. Add R/rbart.R: 0.9-34's rbart.R (without the `if (FALSE)` stubs), sliceSample.R, the six methods and print,
   with adaptations 1-12. Format it with air and make it lintr-clean. The spike left 7 lints: 3
   object_usage, 2 quotes, 2 return.
2. R/tombstones.R: delete the stub bodies. Keep the seven registry rows. Amend the header to say that rbart_vi runs
   0.9-34's implementation after warning once, and that R/rbart.R goes with the file at 1.1-0.
3. NAMESPACE: adaptation 12. The export and the seven S3method lines already stand.
4. Tests. Port main's nine files with these changes:
   - add patterns to the bare `expect_error` calls;
   - seed each fit in test-rbart-groupby.R so its correlation checks are deterministic;
   - drop the hard-coded ranef snapshot in test-rbart-reproducibility.R and keep its same-seed identity checks
     (call 10);
   - keep test-rbart-performance.R's lme4 comparison guarded as it is.

   Add test-rbart-port.R:
   - the warning fires once per session, names stan4bart and 1.1-0, and does not fire from the methods. Count
     warnings with `withCallingHandlers` after resetting the key in `dbarts:::onceWarnState`;
   - one regression check for each fix in adaptations 7 and 10, each small and seeded;
   - gaussian and probit recovery on a simulated design: posterior mean tau inside a bound around the true value,
     and `cor(ranef.mean, b)` above a floor;
   - bartCause's call shape: a dbartsData formula, a literal group.by, and `extract(sample = "test",
     combineChains = FALSE)` against `predict(type = "bart")`;
   - a non-gaussian, non-probit response refused (adaptation 11).

   In test-tombstones.R, replace the stub block with a call to the live function.

Slice 2, manual, NEWS and docs (~250 lines):

5. man/rbart.Rd back from main. Add a Deprecated paragraph first in the description, `print` usage, and a pointer
   to `bartBT` for `sigdf`, `sigquant`, `k`, `power` and `base`, which `bart` no longer takes. Wrap the example in
   `suppressWarnings`. In man/dbarts-deprecated.Rd, drop the rbart aliases, usage and paragraph, and keep one
   sentence naming rbart_vi as deprecated with a link to its page. Add `rbart` to `_pkgdown.yml` beside
   `dbarts-deprecated`.
6. inst/NEWS.Rd 1.0-0:
   - replace "rbart_vi and its methods are removed" with a deprecation item naming stan4bart and 1.1-0;
   - reword the rbart clause of the tombstone-list item to "deprecated, runs as in 0.9-x" (the names stay, so the
     registry check in test-tombstones.R still finds them);
   - add BUG FIXES items for the four 0.9-34 defects in adaptations 7 and 10;
   - add one line on what the port inherits from 1.0 (call 3).
7. docs:
   - the retire-grouped-random-effects.md Status line, plus a short "Reversed in part" paragraph under "What would
     reverse this";
   - its row in docs/design/INDEX.md and in docs/plans/INDEX.md, and a row for this plan;
   - the three rbart_vi statements in bartcore-review-tour.md (VD-read, so plain language);
   - the rbart_vi sentence in docs/design/bart-as-a-component.md;
   - the TODO item rbart-vi-port, removed at landing.

## Verification

- `R CMD INSTALL -l <lib> .`, then with `R_LIBS=<lib>`: `tinytest::test_package("dbarts")` green;
  `run_test_file` on each test-rbart-*.R file and on test-tombstones.R.
- The lint workflow's gates, each on its own exit status: `lintr::lint_package()`,
  `air format --check .`, `Rscript tools/check-rc-codoc.R .`, `tools/check-win-drift.R`,
  `tools/check-doc-freshness.R`. Also `pkgdown::check_pkgdown(".")`, the NEWS parse gate in docs/plans/README.md,
  and `R CMD check --as-cran` on a tarball from a clean copy.
- No existing snapshot moves. tests/cpp and the equivalence compares do not apply (no src/ change).
- Sister replay, recorded in the landing note but not committed:
  - bartCause main's test-02-treatmentFit.R against the slice library: the rbart test passes;
  - stan4bart main's at_home rbart_vi call on four seeds: test deviance within 0.05 of 0.9-34's;
  - the six-seed gaussian and probit comparison above against a main build: posterior mean tau and ranef RMSE
    agree within seed noise.

## Agent-made calls

1. Two defaults come from 1.0's `dbarts()`, not from 0.9-34:
   - binary k is chi(1.5, 2), where 0.9-34 used chi(1.25, Inf) (dec-B106);
   - the proposal mix is 0.6/0/0.4, where 0.9-34 used 0.5/0.1/0.4.

   Alternative: pin both to 0.9-34. That would reproduce 0.9-34's model more closely, but it keeps a k prior the
   package's study rejected. The spike shows the probit posterior agreeing either way.
2. Four 0.9-34 bugs are fixed rather than carried: the skewed-response first draw, the `group.by` lookup, the
   combined-chain new-level `predict`, and saved-fit reload. Alternative: port them verbatim. That is faithful,
   but a skewed response then gives tau near 100 against a truth of 1. It needs VD's eye because "0.9-34's
   implementation carried over" could be read as verbatim.
3. 0.9-34's formals, output shapes and chain handling are kept, not dec-A79's one-chain margin. The combined
   scalar order becomes chain-major through the 1.0 helpers. Alternative: the 1.0 shapes, which would break
   bartCause's `dim(fit$yhat.train)` reads and 0.9-x scripts.
4. The deprecation warning comes from `rbart_vi` only, once per session; the methods are silent. Alternative:
   warn from the methods too, which only repeats itself.
5. `warnOnce` with a plain warning, not `.Deprecated`, which warns on every call and has no once-per-session form.
6. The seven registry rows stay, under the same 1.1-0 expiry, with the header amended. Alternative: a new
   "deprecated" kind. That means editing the registry test and the NEWS list check for a single entry.
7. One file, R/rbart.R, with the slice sampler folded in. Alternative: restore R/sliceSample.R as its own file, one
   more file to delete at 1.1-0.
8. The man page comes back whole, and dbarts-deprecated.Rd keeps a pointer. Alternative: document rbart_vi on
   dbarts-deprecated.Rd, which would put 30 arguments on the deprecation page.
9. Data built with indicators and `na.omit`, as in 0.9-34 and `bartBT`. Alternative: 1.0's categorical columns
   and missing-value handling, which would change which rows and columns a 0.9-x script fits.
10. The hard-coded ranef snapshot is dropped. Alternative: a reference-build test-reproducibility-rbart.R. That
    adds a fifth file to tools/regenerate-snapshots.R and to the cpp-tests workflow list (a .github edit) for a
    function that is gone next release. The engine draws underneath are already pinned by the other four files.
11. No vignette change and no committed comparison harness. The vignette's own random-intercept example stays.
    The comparison numbers go in the landing note.
12. For VD: dec-B130 says the port makes "the released bartCause's grouped fits run again". On 1.0, released
    `bartc()` fails before reaching rbart_vi, on every call and not only grouped ones, because
    `dbartsData@x` is a `dbartsMixedMatrix`. The TODO already expects cross-checks against the released sister
    packages to fail. So the port serves 0.9-x scripts and direct `rbart_vi` callers, not released `bartc()`.
    The register text may want correcting. Nothing in this plan depends on it.

## Landing
