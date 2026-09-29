# rbart-vi-port: rbart_vi back for one release as 0.9-34's R loop

Status: LANDED 2026-09-28 (S1 ac5d8ea2, S2 e975b367, review fixes 94294de5)

agent: sonnet (R port, tests, manual, NEWS; no engine, bridge or header change)
rng: neutral (no existing draw moves; rbart_vi's own draws are new against both 0.9-34 and the tombstone)
budget: ~3150 added lines, nearly all ported: R ~1900 (the 0.9-34 loop, its slice sampler and six methods after
air formatting), tests ~1000, man ~200, NEWS and docs ~150; ~120 removed (stubs, stub tests, stub Rd). Two slices.

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
- The spike's PSOCK runs never called `predict`. An independent review's probe did, and found adaptation 10. The
  same review found D5-D7 and M9, each rechecked here on the spike build and on 0.9-34.

## Sampler calls the loop makes

Each is supported by the 1.0 sampler as-is unless the row says otherwise.

| 0.9-34 call or access | 1.0 | adaptation |
|---|---|---|
| `sampler$setOffset(offset, isWarmup)`, second argument `updateScale` | [`dbartsSampler$setOffset`](../../R/dbarts.R) keeps `(offset, updateScale)`; [`bartcoreSamplerSetOffset`](../../R/bartcore.R) refuses `updateScale = TRUE` only on amplitude-glued samplers, which rbart never builds | none |
| `sampler$run(0L, 1L)`; reads `$train`, `$test`, `$sigma`, `$k`, `$varcount` | same fields; `$train` still carries the offset; saved trees accumulate across calls | none |
| `sampler$getLatents(state$y.st)`, write in place | [`bartcore_getLatents`](../../src/R_interface_bartcore.cpp) fills a supplied numeric of length n | none |
| `sampler$predict(sampler$data@x)` right after `sampleTreesFromPrior` | refuses with keepTrees on and no saved draw (n.burn thins to 0) | 6, call D1 |
| `sampler$sampleTreesFromPrior()`, `$setControl` | unchanged | none |
| `sampler$startThreads()`, `$stopThreads()` | no-op tombstones that warn ([`noOpThreadMethod`](../../R/tombstones.R)) | 5 |
| `model@leaf.hyperprior`, class `dbartsChiHyperprior` | the slot is `leaf.hyperprior` | 11 |
| `data@sigma` | slot unchanged, filled at creation | none |
| `control@keepTrees`, `n.thin`, `n.chains`, `n.threads`, `updateState` | unchanged | none |
| `.Call(C_dbarts_assignInPlace, ...)` | [`assignInPlace`](../../src/R_interface.cpp) still registered | none |
| `.Call(C_rbart_fitted, ...)` in `fitted.rbart` | removed from the bridge | 7 |
| `redirectCall` into `dbartsControl` | `dbartsControl` now has `seed` and `...` | 1 |
| `dbarts(node.prior =, resid.prior =, sigma =)` | all three are warning tombstones on `dbarts()` | 2, 3 |
| `makeCluster`, `clusterExport(rbart_vi_fit, rbart_vi_run)`, `clusterMap` | unchanged; workers load the installed namespace, so the slice sampler resolves there | none |
| `combineChains`, `convertSamplesFromDbartsToBart`, `sampleFromPPD` | combined order is chain-major now; [`sampleFromPPD`](../../R/generics.R) takes `n.chains` | 8 |
| a sampler returned from a PSOCK worker | carries no stored state, so `predict` and `extract(type = "trees")` refuse | 10 |

## Adaptations

What the 1.0 sampler forces, numbered for the table. Changes to what 0.9-34 computed are not here; they are the
calls at the end, for the maintainer.

1. Drop `seed` and `...` from the forwarded control call. `seed` is a
   [`dbartsControl`](../../R/dbarts.R) formal now, and passing it would give every serial chain the same engine
   seed; the loop keeps 0.9-34's own `set.seed` handling. `...` would otherwise be copied in as an empty argument,
   which the spike hit as an error on every call. Force `keepFits = TRUE` alongside the `keepTrainingFits = TRUE`
   0.9-34 already forces: the loop reads `$train` every sweep, and `keepFits` is a 1.0 control setting a caller
   can pass through.
2. Build the sampler as [`bartBT`](../../R/bart.R) does: `leaf.prior = normal(k)` only when `k` is supplied, the
   residual prior `chisq(sigdf, sigquant)` on an `"auto"` family through [`withResidPrior`](../../R/family.R),
   and `sigest =`. No tombstone warning fires. Which k prior and move mix apply when `k` is not supplied is
   call M1.
3. Build the data with `factors = "indicators"` and `na.action = stats::na.omit`, 0.9-34's dummy expansion and
   row rule, as `bartBT` does.
4. Deprecation: the first statement of `rbart_vi` calls [`warnOnce`](../../R/utility.R) under key
   `tombstone.rbart_vi`: "'rbart_vi' is deprecated and is removed in dbarts 1.1-0; grouped random effects live in
   stan4bart (stan4bart::stan4bart), whose group-spread prior differs, so results move." The methods do not warn.
5. Delete the `$startThreads()`/`$stopThreads()` pair.
6. With n.burn thinned to zero and keepTrees on, take the initial fit with keepTrees off, then turn it back on.
7. `fitted.rbart` averages in R rather than calling the removed C helper: `colMeans` of yhat plus the matched
   ranef columns, `pnorm` first for a binary fit. An unmatched group gives NA, not the out-of-bounds read the C
   helper did.
8. The methods keep 0.9-34's shapes: a private copy of 0.9-34's `combineOrUncombineChains` (a one-chain fit gets
   no chain margin) instead of [`combineOrUncombineChains`](../../R/generics.R) (dec-A79 adds one).
   `sampleFromPPD` gets `n.chains`. `browser()` becomes `stop()`.
9. Refusals happen in `rbart_vi`, from the built `dbartsData`, before any chain starts. Two checks:
   - A response other than continuous or binary (a Surv response would resolve to aft under `"auto"`) is refused
     with an rbart_vi message. A 3-level factor is already refused by `dbarts()`.
   - A binary response with weights other than 0 and 1 is refused, naming rbart_vi. See call M9.

   Raised inside a PSOCK worker instead, the error does reach the caller, but only after 0.9-34's handler warns
   "error running multithreaded, defaulting to single" and retries serially. Checked on the spike with weighted
   probit.
10. `$storeState()` on every kept sampler at the end of `rbart_vi_fit`, always, before the sampler leaves a
    worker. A sampler returned from a PSOCK worker otherwise carries no stored state. So a default multi-chain
    fit (n.threads above one) cannot `predict` or `extract(type = "trees")`: it fails with "samplers cannot be
    re-created without a stored state". The reviewer's probe hit this, and 0.9-34 is fine there. The same call
    fixes serial reload (call D4).
11. Detect a modeled k from the sampler's leaf-prior reader, `getLeafPrior()[, "k.has.hyperprior"]`, not from
    `model@leaf.hyperprior`. The run's `$k` is NULL for a fixed k too, but
    the sample buffers are sized before the first run.
12. `importFrom(stats, dcauchy, dgamma)`, since the tau priors use them; the graphics imports are already there.

Kept from 0.9-34 as they are: the warmup scale re-anchoring (`setOffset(updateScale = TRUE)` during burn-in), the
formals and their defaults (`seed = NA_integer_`, `combineChains = FALSE`, `k = 2.0` shown but applied only when
supplied), the "verbose output disabled" warning, the result fields, and a `k` that is evaluated rather than parsed
(`chi()` does not resolve in either build).

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
call does. The dbarts-1.0 branch has no rbart call left and needs nothing. A released bartCause call with
binary treatment and non-0/1 weights now meets call M9's refusal.

**stan4bart, released (main).** Its skip_on_cran() test `nonlinearities are estimated well` calls
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
   with adaptations 1-12 and the calls as ruled. Format it with air and make it lintr-clean. The spike left 7 lints: 3
   object_usage, 2 quotes, 2 return.
2. R/tombstones.R: delete the stub bodies. Keep the seven registry rows. Amend the header to say that rbart_vi runs
   0.9-34's implementation after warning once, and that R/rbart.R goes with the file at 1.1-0.
3. NAMESPACE: adaptation 12. The export and the seven S3method lines already stand.
4. Tests. Port main's nine files with these changes:
   - add patterns to the bare `expect_error` calls;
   - seed each fit in test-rbart-groupby.R so its correlation checks are deterministic;
   - drop the hard-coded ranef snapshot in test-rbart-reproducibility.R and keep its same-seed identity checks
     (agent-made call 8);
   - keep test-rbart-performance.R's lme4 comparison guarded as it is.

   Add test-rbart-port.R:
   - the warning fires once per session, names stan4bart and 1.1-0, and does not fire from the methods. Count
     warnings with `withCallingHandlers` after resetting the key in `dbarts:::onceWarnState`;
   - one regression check for each defect fix the maintainer takes (calls D1-D7), each small and seeded;
   - `predict` and `extract(type = "trees")` on a two-chain PSOCK fit (adaptation 10), and a serial fit that
     predicts identically after `saveRDS` and `readRDS`;
   - refusals before any chain (adaptation 9), including under n.threads = 2, with no "defaulting to single"
     warning;
   - gaussian and probit recovery on a simulated design: posterior mean tau inside a bound around the true value,
     and `cor(ranef.mean, b)` above a floor;
   - bartCause's call shape: a dbartsData formula, a literal group.by, and `extract(sample = "test",
     combineChains = FALSE)` against `predict(type = "bart")`;

   In test-tombstones.R, replace the stub block with a call to the live function.

Slice 2, manual, NEWS and docs (~250 lines):

5. man/rbart.Rd back from main. Add a Deprecated paragraph first in the description, `print` usage, and a pointer
   to `bartBT` for `sigdf`, `sigquant`, `k`, `power` and `base`, which `bart` no longer takes. Wrap the example in
   `suppressWarnings`. The page says data handling follows 0.9-x (indicator columns, `na.omit`), not `bart`'s,
   and that a binary response takes only 0/1 weights (call M9). In man/dbarts-deprecated.Rd, drop the rbart aliases, usage and paragraph, and keep one
   sentence naming rbart_vi as deprecated with a link to its page. Add `rbart` to `_pkgdown.yml` beside
   `dbarts-deprecated`.
6. inst/NEWS.Rd 1.0-0:
   - replace "rbart_vi and its methods are removed" with a deprecation item naming stan4bart and 1.1-0;
   - reword the rbart clause of the tombstone-list item to "deprecated, runs as in 0.9-x" (the names stay, so the
     registry check in test-tombstones.R still finds them);
   - add BUG FIXES items for the defect fixes the maintainer takes (calls D1-D7);
   - add the weighted-binary refusal (call M9), and one line on the defaults as ruled (call M1).
7. docs:
   - the retire-grouped-random-effects.md Status line, plus a short "Reversed in part" paragraph under "What would
     reverse this";
   - its row in docs/design/INDEX.md and in docs/plans/INDEX.md, and a row for this plan;
   - the three rbart_vi statements in bartcore-review-tour.md (VD-read, so plain language);
   - the rbart_vi sentences in docs/design/bart-as-a-component.md, docs/design/heteroscedastic.md and
     docs/design/correlated-outcomes.md;
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
  - stan4bart main's skip_on_cran rbart_vi call on four seeds: test deviance within 0.05 of 0.9-34's;
  - the six-seed gaussian and probit comparison above against a main build: posterior mean tau and ranef RMSE
    agree within seed noise.

## Calls

### For the maintainer, before implementation

M1. Defaults when `k` is not supplied. Ruled (a), maintainer 2026-09-28: "Use 1.0 defaults." Two options:
    (a) 1.0's `dbarts()` defaults: binary k ~ chi(1.5, 2) (dec-B106) and moves birth/death 0.6, swap 0,
        change 0.4. The fit then matches every other 1.0 binary fit and uses the k prior the package's study
        chose; it differs from 0.9-34's model.
    (b) 0.9-34's: binary k ~ chi(1.25, Inf), and moves 0.5/0.1/0.4. A 0.9-x script then fits 0.9-34's model; the
        port carries a k prior the study rejected, set by rbart_vi alone.

    Gaussian fits are unaffected: k is fixed at 2 in both. The spike ran (a); its probit tau and intercepts agree
    with 0.9-34 within seed noise, so the choice moves the prior and mixing, not a visible bias on that design.

Each defect fix below states what 0.9-34 does, what the fix does, and the verbatim alternative.

D1. First draw on a skewed response. Ruled fix, maintainer 2026-09-28: "Fix it." 0.9-34 starts the intercept draw from `predict` on prior trees, which returns
    the response midpoint; on a skewed response the midpoint-mean gap goes into the intercepts and stays. Spike,
    n = 1000, 10 groups, true tau 1: mean tau 105, 128, 107 on three seeds. Fix: start from one sweep's
    `$train`; the port then gives 1.56, 1.10, 0.78 (the bartcore-era fix, archived plan
    rbart-custom-prior-divergence.md). Verbatim: tau near 100.
D2. `group.by` lookup. Ruled fix, maintainer 2026-09-28 ("Fix it, unless it somehow appeared to be intentional in docs or examples"): it does not; 0.9-34's NEWS promises the fall-through, and main's test-rbart-groupby.R second case exercises it while the bug silently groups by x_1. main's test-rbart-error.R not_a_symbol case passes only through the bug; ported, it expects "'group.by' not found". 0.9-34 takes a symbol that is not a column of `data` as `data`'s first column
    ([R/rbart.R:84](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/rbart.R#L84),
    `which.max` of all-FALSE is 1), silently: the spike fit 200 groups where 5 were meant. Same for
    `group.by.test`. Fix: take the column only when it is there, then fall through to the caller's scope as
    0.9-34 does. Verbatim: the silent wrong grouping.
D3. New level in `predict` with several chains and `combineChains = TRUE`. Ruled fix, maintainer 2026-09-28: "Fix it." 0.9-34 stops with "subscript out of
    bounds": it tests the stored ranef's dimension rather than the combined one
    ([R/generics.R:227](https://github.com/vdorie/dbarts/blob/cb29055019449614b4085eb7d30431aa39f790f5/R/generics.R#L227)).
    Fix: test the combined one, and draw the new-level intercepts with tau in the same layout. Verbatim: the
    error.
D4. Saved fit read into a new session. Not a fork: the store is forced on 1.0; noted to the maintainer 2026-09-28. 0.9-34's `predict` silently returns different values. Adaptation 10's
    `$storeState()` makes the reloaded fit predict identically. There is no verbatim option on 1.0: without the
    store, reload and every PSOCK fit refuse to predict. NEWS files it as a fix to reloaded fits, not as the
    whole reason for the store.
D5. Weights in the intercept draw. Ruled fix, maintainer 2026-09-28: "fix it to do proper precision weights" - the lme4 and dbarts convention (row variance sigma^2 / w_i), the same update the engine gives a leaf; the maintainer asked what a varying-intercept hierarchical linear model does before ruling. 0.9-34 draws each group intercept with precision n_j / sigma^2 around the
    unweighted residual mean, while BART's sigma is on the weighted scale. The intercepts' posterior then
    depends on the weights' overall scale. Reviewer's probe, same data, constant weights 1 against 1/n: posterior
    sd of the intercepts 0.366 against 0.175 on 0.9-34, 0.262 against 0.125 on the port. Fix: the weighted
    conditional, precision sum of w_i over group j divided by sigma^2, around the weighted residual mean.
    Verbatim: intercept uncertainty that moves with a rescaling of the weights.
D6. Test offset double count. Ruled fix, maintainer 2026-09-28: "Fix it." When a caller gives an offset and a test set with the same row count, 0.9-34's
    `setOffset` also copies each sweep's training offset, intercepts included, onto the test rows. The test fits
    then carry the training rows' intercepts, and extract adds the test groups' intercepts again. Probe, test =
    train rows: up to 2.78 between test and train fits on 0.9-34, 1.60 on the port, where they should agree.
    Fix: set `testUsesRegularOffset` to FALSE on the sampler's data for the loop, then restore it. Verbatim: the
    double count.
D7. `seed =` and the caller's random stream. Ruled (b), maintainer 2026-09-28: "Fix it." A stream absent before the call is removed after it. 0.9-34 means to restore the stream after a seeded fit, but
    `.Random.seed <- oldSeed` assigns a local, so the stream is left where the fit ended; checked on the port. The
    options:
    (a) delete the dead restore, and document that `seed =` resets the global stream;
    (b) assign in the global environment, as `packageRbartResults` already does for `result$seed`, so a seeded
        fit leaves the caller's stream untouched.

    Verbatim is (a) without the documentation.

M9. Weighted binary fits. Ruled refuse up front, maintainer 2026-09-28: "Refuse up front." (consistent with dec-B13). 0.9-34 fit probit with any weights; 1.0's `dbarts()` refuses probit weights other than
    0 and 1. The port refuses before any chain (adaptation 9), naming rbart_vi and stan4bart, and says so in
    NEWS and the manual. The alternative is working around the refusal inside rbart_vi. That would need a
    weighted probit the 1.0 sampler deliberately does not provide, which is engine work outside dec-B130.

### Agent-made, noted for the register

1. 0.9-34's formals, output shapes and chain handling are kept, not dec-A79's one-chain margin. The combined
   scalar order becomes chain-major through the 1.0 helpers. Alternative: the 1.0 shapes, which would break
   bartCause's `dim(fit$yhat.train)` reads and 0.9-x scripts.
2. The deprecation warning comes from `rbart_vi` only, once per session; the methods are silent. Alternative:
   warn from the methods too, which only repeats itself.
3. `warnOnce` with a plain warning, not `.Deprecated`, which warns on every call and has no once-per-session form.
4. The seven registry rows stay, under the same 1.1-0 expiry, with the header amended. Alternative: a new
   "deprecated" kind. That means editing the registry test and the NEWS list check for a single entry.
5. One file, R/rbart.R, with the slice sampler folded in. Alternative: restore R/sliceSample.R as its own file, one
   more file to delete at 1.1-0.
6. The man page comes back whole, and dbarts-deprecated.Rd keeps a pointer. Alternative: document rbart_vi on
   dbarts-deprecated.Rd, which would put 30 arguments on the deprecation page.
7. Data built with indicators and `na.omit`, as in 0.9-34 and `bartBT`. Alternative: 1.0's categorical columns
   and missing-value handling, which would change which rows and columns a 0.9-x script fits.
8. The hard-coded ranef snapshot is dropped. Alternative: a reference-build test-reproducibility-rbart.R. That
   adds a fifth file to tools/regenerate-snapshots.R and to the cpp-tests workflow list (a .github edit) for a
   function that is gone next release. The engine draws underneath are already pinned by the other four files.
9. No vignette change and no committed comparison harness. The vignette's own random-intercept example stays.
   The comparison numbers go in the landing note.

10. Refusal test for "continuous or binary" reads the built `dbartsData`: a Surv response (its `survivalStatus`
    attribute) or a multinomial count matrix is refused. Alternative: refuse from the formula, which would miss a
    `dbartsData` passed as `formula`.
11. The caller's random stream is restored by an `on.exit` in `rbart_vi`, so an error mid-fit restores it too, and by
    NULL-safe helpers that also fix `packageRbartResults`, whose restore assigned NULL when no stream existed (an
    error on the next draw). Alternative: restore on success only.
12. D5 covers zero-weight rows: a group whose weights total 0 draws its intercept from the prior, where 0.9-34's
    `mean` of an empty group gave NaN. A binary fit's 0/1 weights live in the sampler's active rows, not in
    `data@weights`, so the loop reads them from there (review fix; the first landing missed it and let excluded rows
    inform their group).
13. D3's new-level draws are not the same numbers in the combined and split layouts for one seed (each layout fills
    its array in its own order); only the measured levels agree, and the scales follow each draw's own tau.
14. The eight ported test files mark the deprecation key as warned at their top, so the suite's warning list is not
    the deprecation repeated; test-rbart-port.R and test-tombstones.R reset it and count.
15. The weighted-binary and Surv refusals are tested through direct calls; D2 is tested through direct calls too,
    since `do.call` evaluates a symbol before `rbart_vi` sees it.
16. TODO item rbart-vi-port is left in place until the landing record, since the plan says it is removed at landing.
17. `predict` on an rbart fit with an unnamed matrix `newdata` passes through the 1.0 sampler's own
    "'test' is unnamed but 'x' had named predictors" warning; not changed here.

## Landing

### Landing note (2026-09-28, wt/rbart-port off 25911928, not pushed)

What landed. Slice 1 (code and tests): `R/rbart.R` (1894 lines: `rbart_vi`, its loop and slice sampler, six methods
and print, adaptations 1-12, fixes D1-D7, refusals M9 and adaptation 9, M1 as ruled (a)), the seven stubs deleted
from `R/tombstones.R` with the registry rows kept and the header amended, `importFrom(stats, dcauchy, dgamma)`, the
nine 0.9-34 test files ported (patterns added, groupby fits seeded, ranef snapshot dropped), `test-rbart-port.R`
(421 lines, 51 expectations), and the live-function block in `test-tombstones.R`. 3350 added, 75 removed. Slice 2:
`man/rbart.Rd` back with the Deprecated paragraph, print usage, bartBT pointer, data and weights notes and a
`suppressWarnings` example; `dbarts-deprecated.Rd` reduced to a pointer; `_pkgdown.yml`; `inst/NEWS.Rd` (a
deprecation item, the tombstone-list reword, seven BUG FIXES items, M1 and M9 in the deprecation item); the
design and plan records, the review tour, and the design-doc sentences (bart-as-a-component, correlated-outcomes, heteroscedastic; the last was missed at
first). 262 added, 48 removed.

Gates, macOS arm64, R from this host, private library, tree at the slice 2 commit. Full tinytest suite: 9499
expectations, 0 failures. `lintr::lint_package()`: 0 lints. `air format --check .`: clean. `check-rc-codoc.R`,
`check-win-drift.R`, `check-doc-freshness.R`: exit 0. `pkgdown::check_pkgdown(".")`: no problems. NEWS parses,
332 entries (325 before). `R CMD build` and `R CMD check --as-cran --no-manual` on the tarball of a clean copy:
Status 1 NOTE, the Date field over a month old, present before this change; no WARNING. The 25 gates of
exact-gates.yaml in quick mode: 25 of 25 pass (the loop only; the two cross-host compares need the reference
build and were not run, and nothing under src/ changed). Mutation check: test-rbart-port.R's sections for D1, D2,
D3, D5, D6, D7 and M9 fail against a build of 0.9-34 and pass on the port; the D4 section fails there in one
expectation; the recovery and call-shape sections pass on both.

Sister replay. stan4bart main's `nonlinearities are estimated well` rbart_vi call (probit, 100 rows, default
1000/1000), read from main with git show and run from a scratch script against this library, seeds 1 to 4:
test-set deviance 1.45, 1.41, 1.41, 1.42 on the port against 1.38, 1.40, 1.38, 1.40 on 0.9-34. Seed 1 differs by
0.07, above the plan's 0.05, and all four sit far inside the test's own 1.35x slack. The plan's own spike numbers
(1.44, 1.39, 1.40, 1.41 and 1.42, 1.44, 1.39, 1.43) used other seeds. stan4bart's bartcore branch has no rbart_vi
call, so its tests were not run. bartCause and the six-seed tau comparison were not replayed.

What the plan got wrong or left open. Slice 1 ran 3350 lines against about 2800, from air's formatting of the
ported tests. Nothing else in the plan was wrong. Main's
`test-rbart-weights.R` bare `expect_warning` would have passed on the deprecation warning alone; it now
carries the pattern of the warning it means.

### Review fixes (third commit, same day)

Seven review findings applied. Binary zero-weight rows now read the active rows (a group of all-zero weights draws
from the prior, sd near tau; test added, 0.41 of tau before, 1.15 after). `rbart_vi` builds one sampler with the
chain's arguments before any chain and discards it, so a 3-level factor response and `k = -1` refuse once under
`n.threads = 2` (both added to the M9 test); the build leaves the caller's random stream as it found it (the seeded-chain identity test needs that). `keepFits` is forced
on, with a test. NEWS drops the claim that PSOCK fits could not predict on 0.9-34 and keeps the reload sentence.
The manual's shared-argument pointer names `bart`, with `sigdf`, `sigquant`, `power` and `base` as in `bartBT`.
heteroscedastic.md no longer says grouped intercepts were removed entirely. `predict` and `extract(type = "trees")`
on a fit saved by 0.9-x stop with a message to refit; the test is a sampler environment lacking the `activeRows`
field, which 0.9-x samplers lack (a real 0.9-x fit was checked by hand). New agent-made call: the 0.9-x detector is
that field's absence; the state-based `refuseLegacyState` cannot be used, since reading `$state` on such a sampler
already fails.

Gates after the review fixes, same environment: full tinytest suite 9508 expectations, 0 failures (run without
suppressWarnings); `lintr::lint_package()` 0 lints; `air format --check .` clean; codoc, win-drift and
doc-freshness exit 0; `check_pkgdown` clean; NEWS 332 entries; `R CMD check --as-cran --no-manual` on a clean
tarball, one NOTE (the Date field), no WARNING. The exact gates were not rerun: nothing they read changed since the
25 of 25 pass. Line counts of this commit are in its diffstat, about 190 added and 40 removed.

