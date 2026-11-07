# Review 4, slice S-D: the R surface

Status: LANDED on main 2026-11-06 (branch review4-r-surface).

Each item is a finding of the fourth whole-repo review. RNG class: neutral. No fix changes a sampler
draw. families-01 changes what an `nbinom()` fit reports, and so its predictive draws and
log-likelihood, which had read a drawn shape as held; the rest refuse a call that failed or gave NaN,
name the family the caller wrote, or move a warning. The one engine edit is a debug assert, compiled
out of the shipped build.

## Fixes as built

**families-01: `nbinom()` written as a call draws the shape as the token does.** `fitDescriptors`
reads the model's resolved shape (`isTRUE(model@shape > 0)`) where it read the family object's
setting, so `family = nbinom()`, `nbinom(shape = NULL)` and `dbartsFamilies$nbinom()` carry the shape
draws, keep `shape.raw` under keepTrees, print no "Fixed, not sampled" line, and draw the ppd and
loglik at each draw's shape. `nbinom(4)` reports its shape in `fit$fixed` alone. Tests:
test-nbinom.R, each call form against `"nbinom"` under one seed (yhat.train, shape, shape.raw, fixed,
`extract(type = "ppd")`, `predict(type = "ppd")`, loglik, the summary text with the call held aside),
and `nbinom(4)`.

**families-02: an `offset.test` forwarded through a wrapper's dots.** `dbartsData` replaces a `..N`
`offset.test` by the expression it was written as (`recoverForwardedArgument`) and reads its names
where it was written, as the training offset already does. Tests: test-data-testOffset.R, a name,
`offset + 1`, a column of `test`, one and two wrapper layers and `dbarts()`, each identical to the
direct call, and a wrapper whose caller holds the test offset as a local of its own, on the matrix
and formula interfaces.

**families-05.1: a refused factor response names "student".** `resolveClassificationFamily` takes
`shownFamily`, the caller's token, for its two "cannot fit a" messages; `resolveSamplerSpec` passes
`requestedFamily`. Tests: test-family.R, the token and `dbartsFamilies$student(df = 3)` on a 4-level
and a 2-level factor.

**families-05.4, families-V-02 and dec-B602: no rescale advice on a count or a category.**
`withBinaryResponsePrecision` holds the precision warning for the ordinal, negative binomial,
multinomial and hazard families, for a binary family fit to a single class, and for any categorical
response encoding; a binary family keeps it beside its refusal of a near-constant response that is
not one class. `bart2Hurdle` holds it during its data build (`holdingResponsePrecision`) and raises it
only once `splitHurdleResponse` returns, so a refusal drops it. An ordinal or multinomial fit whose
response shows one category is fitted with one warning, `warnSingleCategory`: "family \"ordinal\" is
fitted to a response with a single observed category; the prior alone drives the fitted
probabilities". Tests: test-data-precision-warning.R, warnings counted: ordinal and auto on an
ordered factor with one level at each of a, b and c, through `bart()` and `dbarts()`; multinomial
likewise and on a count matrix with one counted category; two observed categories of three silent;
nbinom constant 3 silent; hurdle constant 2.5 refused with no warning; probit constant 2 refused with
the warning; gaussian constant 5 warns. test-family.R keeps the near-constant probit case.

**models-01 and dec-B597: a variance forest of two or more trees refuses df 2 or below.**
`refuseVarianceCalibration(residPrior, numTrees)` stops on a `chisq()` residual prior with df 2 or
below beside two or more variance trees: "a variance forest of 40 trees does not support a residual
prior df of 1.5: its leaves are calibrated to the prior mean of sigma^2, which is infinite at df 2 or
below; give chisq() a df above 2, or use varianceForest(n.trees = 1)". `resolveSamplerSpec` calls it
once the variance tree count is resolved, `$setModel` with the sampler's variance record. The bridge
repeats it in the same words, from the model's variance record, in `parseSamplerSpecification` and in
`bartcore_setModel` before anything is written. `Chain::calibrateVarianceLeaf` asserts a finite leaf
df in a debug build. Help: dbarts.Rd's `variance`, varianceForest.Rd's `n.trees`,
dbartsSampler-class.Rd's `setModel`. Tests: test-heteroscedastic.R, df 1.5 and 2 refused through
`dbarts(gaussian(sigma = chisq()))`, `dbartsSpec()` (no sampler built, the R check alone),
`bart(variance = varianceForest(n.trees = 2))` and `bart(sigdf = )`; one variance tree at df 1.5
fitted; `$setModel` to df 1.5 refused, the sampler then identical to an untouched twin, df 2.5 taken;
a model handed straight to `new("dbartsSampler", ...)` refused by the bridge, and one handed
straight to the bridge's `setModel` entry, the sampler then running as its twin.
test-aft-heteroscedastic.R, `aft(sigma = chisq())` refused and one tree accepted.

**models-V-01, models-V-02, group B and dec-B598: the prior constructors check their scalars.**
`priorScalarRefusal(value, name, inRange, kind)` refuses a value that is not numeric, is an object,
is not of length 1, is NA or is out of range, naming the argument as written; `refusePriorScalar`
stops with it. The constructors call it before `newValidated`, the class validities share it, length
first, for a direct `new()`:

| constructor | argument | refused unless |
|---|---|---|
| `chisq` | `df` | positive finite; "to hold sigma at a value, use fixed()" |
| `chisq` | `quant` | in (0, 1) |
| `fixed` | `value` | positive finite |
| `chi` | `df`, `scale` | positive (an infinite scale is the improper limit) |
| `cgm` | `power`; `base` | positive; in (0, 1) |
| `dart` | `a`, `b`, `alpha`, a stated `rho` | positive finite |
| `dart` | a stated `update.delay` | non-negative |
| `dart` | `update.alpha` | TRUE or FALSE |

`dart(rho = NA)` and `dart(update.delay = NA)` keep their own message. A matrix or array of split
probabilities with at most one dimension longer than 1 is the vector it holds, named by that
dimension's names, as a named vector is; `cgm` refuses any other in one line. A forest's stated
`power` and `base`, in `forest()` at creation, through `$setForest`, in `varianceForest()` at creation
and through `$setVarianceForest`, are held to `cgm`'s checks (`validateForestTreePrior`), so a
string, a logical or a factor is refused there as `bart()` refuses it, where they had been coerced.
`normal` already refused each such value by name and is unchanged. The bridge's residual prior parse refuses a non-finite df in dec-B598's words. Help:
dbartsPriors.Rd's `chisq` and `fixed` items. NEWS: `bart(sigdf = Inf)` refused, where 0.9-34
returned NaN. Tests: test-model-priors.R, each constructor at NA, NA_real_, NaN, c(1, 2), "a" and its
range ends, the accepted ends (`chisq(1e6)`, `chi(1.5, Inf)`, `cgm(power = 3L)`), the validities
through `new()`, split probabilities as a one-column, one-row, named and three-dimensional array
taken and a 2 x 2 matrix refused through `cgm` and `bart`, the forest doors' power and base refusals
at each of the four, an infinite df through the
`dbarts()` door, the aft constructor, `bart(sigdf = Inf)` and `bart(sigdf = NA)`, and the bridge
backstop through `new("dbartsSampler", ...)` and `$setModel`. test-model-errors.R,
test-setmodel-parameters.R and test-forest-arguments.R take the new words.

**NEWS.** BUG FIXES: an `offset.test` passed through a wrapper's `...` is read as written (0.9-34
stopped with "cannot find test offset '..1'"). USER-VISIBLE CHANGES: a residual prior of infinite df
is refused, naming `fixed()`.

## rfit-10: k = Inf stays accepted

`normal(k = Inf)` was already accepted, documented and tested as a value: the class help's
multinomial `setLeafPrior` takes a fixed k "Inf included", the multiforest-leaf-prior-writer plan
states it, the forest-door plan extends the reader's NULL exception to it, and
test-multiforest-leaf-prior-writer.R, test-setmodel-forests.R, test-forest-lone-no-basis.R,
test-forest-door.R and test-binaryResponse-hyperprior.R's k.infinite kind fit with it. It stays
accepted, and the planned refusal is not built. bart.Rd's `k`, dbartsPriors.Rd's `normal` and
dbartsSampler-class.Rd's `setLeafPrior` say what it fits: every tree held at zero, so the fit is the
offset, plus on a continuous response the center of the range `y - offset` spans, about which the
response is scaled.

## Calls made

1. The precision warning is held by family and by encoding: ordinal, nbinom, multinomial and hazard
   families, a single-class binary fit, and any categorical encoding. Alternative: the plan's
   ordinal and nbinom alone, which passes it on for "auto" over a categorical response, the family
   being unresolved where the warning is held. A first build also held it for every binary family;
   test-family.R pins the warning beside a probit refusal of a near-constant response that is not one
   class, so a binary family holds it for a single class only.
2. The hurdle re-raise stays after `splitHurdleResponse`, as the plan places it, though every
   constant response the warning flags is refused there today. Alternative: drop the re-raise.
3. The bridge's `$setModel` backstop reads the model record's variance tree count, since the
   sampler's shape knows only whether a variance forest exists; `keepModelRecords` keeps that record
   equal to the sampler's. Alternative: a new accessor on the sampler shape (a facade change).
4. The R `$setModel` check and the bridge backstop are redundant layers with one message: removing
   the R check alone leaves the tests passing, removing both fails them. Alternative: drop the R
   check (the plan names both).
5. The message wording for dart's `a`, `b`, `alpha`, `rho`, `update.delay` and `update.alpha`, and
   chi's "positive number" for an accepted infinite scale, follow the plan's shape for the arguments
   it does not quote.
6. The split-probability matrix is refused in `cgm` before the slot, in its own line, instead of
   passing the validity's message through `rethrowValidityError`.
7. A numeric ordinal response with one value stays refused (it states no categories); a factor with
   one observed level is fitted per dec-B602.
8. The debug assert is on the computed leaf df (finite), which df 2 itself passes; the refusal at df 2
   is the R surface's and the bridge's, per dec-B597.
9. A split-probability matrix or array with at most one dimension longer than 1 is taken as its
   vector, named by that dimension's names; any other is refused. Alternative: refuse every matrix,
   as the first build did, which refused an equal-probability one-column matrix that main fit.
10. The forest doors' power and base refusals keep their "forest 'power'" and "forest 'base'" words,
    the variance forest's at creation included, which had reached the bridge's "variance power must
    be positive". Alternative: the bridge's words at creation and the forest words live.

## Gates

Run on arm64 macOS, at most two cores.

- Full tinytest suite, one process, shipped build: 41395 tests, 0 failed, after the review's fix
  round (41366 before it).
- Equivalence trio, reference build, `--bitwise --strict-coverage`, against the current baselines:
  gaussian 55 of 55 scenarios identical draws, BCF 19 of 19, multinomial 11 of 11; none skipped.
- The four seeded-drift snapshot files, reference build: 27 tests, 0 failed, unchanged.
- Exact gates in `quick`, the main job's list: 44 of 44 passed (33 gates and 11 arms).
- The monotone job of the same workflow in `quick`: the successive-conditional check and the
  enumeration gate both pass.
- `R CMD check --as-cran` on a tarball built from a clean copy: 1 NOTE, the incoming-feasibility one
  (the Date field, and an expired certificate on an external URL in bartBT.Rd and pdbart.Rd, pages
  this slice does not touch); the tests pass.
- tests/cpp, which compiles the debug assert in: 431 checks, all passed.
- The exact gates, `R CMD check` and tests/cpp ran before the fix round, which touches R, tests and
  help only; the suite, the trio, the snapshot files and the checks below ran after it.
- `lintr::lint_package()` against the slice's build: no lints. `air format --check .`: clean.
  check-rc-codoc: 53 methods. check-win-drift: OK. check-doc-freshness: OK.
  `mutation-battery.R verify-anchors`: every anchor resolves.

## Mutants

Each mutant was installed on its own and run against the test file named.

| mutant | fix | outcome |
|---|---|---|
| `fitDescriptors` reads the family object's shape setting again | families-01 | killed: test-nbinom.R, 21 failures in the call-form block |
| the forwarded `offset.test` is not recovered | families-02 | killed: test-data-testOffset.R stops, "cannot find test offset '..5'" |
| the remapped family is named in the refusal | families-05.1 | killed: test-family.R, 4 failures |
| ordinal and nbinom pass the precision warning on | families-05.4 | killed: test-data-precision-warning.R, the constant count |
| the hurdle raises the warning before its refusal | families-V-02 | killed: the hurdle case |
| no single-category warning, multinomial; then ordinal | dec-B602 | killed: test-data-precision-warning.R, each |
| no R check at creation | models-01 | killed: the `dbartsSpec()` case |
| the R check at df below 2 only | dec-B597 | killed: `dbartsSpec()` at df 2 |
| neither `$setModel` layer refuses | models-01 | killed: the refusal and the twin comparison |
| no bridge check at creation | models-01 | killed: the `new("dbartsSampler", ...)` case |
| no R check in `$setModel`, the bridge's kept | models-01 | survives: the bridge refuses in the same words, before anything is written (Call 4) |
| `chisq()` checks nothing of df | models-V-01, V-02 | killed: test-model-priors.R, the NA rows and `bart(sigdf = NA)` |
| the chisq validity takes an infinite df | models-V-01 | killed: the direct `new()` case |
| the bridge takes an infinite df | models-V-01 | killed: the bridge cases at creation and `$setModel` |
| `cgm()` checks nothing of power; `fixed()` of value; `chi()` of scale; `dart()` of a, b, alpha, rho | group B | killed: test-model-priors.R, each |
| `cgm()` passes a split-probability matrix on | group B | killed: the `cgm` and `bart` cases |
| the forwarded `offset.test`'s names read where the wrapper is | families-02 | killed: the caller's local, "cannot find test offset 'localOffset'" |
| no bridge check in `$setModel`, the R check kept | models-01 | killed: the direct bridge `setModel` case and its twin |
| every split-probability matrix refused | group B | killed: the one-column matrix |
| every split-probability matrix taken | group B | killed: the 2 x 2 through `cgm` and `bart` |
| a forest's power and base coerced again | group B | killed: 16 failures over the four doors |
| the variance forest's stated power unchecked at creation | group B | killed: the `varianceForest()` creation case |

## Landing note

Landed by fast-forward at 862424a55 after two reviews (the second over the fix round). On that commit: the nine
per-push workflows green by dispatch (R CMD check, cpp-tests, sanitizers, equivalence, exact-gates, lint,
doc-freshness, pkgdown, revdep-smoke); the x86 leg (Linux, x86-64): tests/cpp 431 checks plain and under
ASan/UBSan with the leak sanitizer, no report; the full tinytest suite at home, 41523 tests, 0 failed; the
equivalence trio in statistical mode, 55 scenarios compared and indistinguishable, BCF and multinomial within the
cross-host bound. forest() and varianceForest() now refuse a string, logical or factor power or base that the
development build took; both are new in 1.0-0, so no NEWS line.
