# small-r-batch: an unstated thread count falls back to one, a forced setPredictor returns NULL, a NULL power or base is the default, two forest-selection texts are corrected

Status: PLANNED (dec-B306, dec-B310; the last two items are root TODO backlog).

agent: sonnet implementer, one (R, tinytest, manual); opus reviewer.
rng: NEUTRAL, bit for bit, for every call accepted before and still accepted. A call that stopped with an
error before (an unstated thread count with uncountable cores, `power = NULL`, `base = NULL`) now runs and
draws as the same call with the stated default does. No draw of an accepted call moves.
window: pre-release, before the merge to main.
budget: ~300 lines beyond this file.

## Goal

A fit that says nothing of threads runs on one thread where the cores cannot be counted. A forced
`setPredictor` returns `NULL` invisibly, as 0.9-34's did. A `NULL` `power` or `base` is the default. Two
wordings that no longer fit the argument the caller wrote are corrected.

## Context

- [`dbartsControl`](../../R/dbarts.R), [`bart`](../../R/bart.R), [`rbart_vi`](../../R/rbart.R) and
  [`xbart`](../../R/xbart.R) default `n.threads` from `guessNumCores()`, which is `NA` where the cores cannot
  be counted; `bartBT` defaults `nthread` to 1 and is not affected. These four are every door that does so.
- [`naThreadsMessage`](../../R/A_class.R) is the refusal of an `NA` count, raised by the control's validity
  and by `rbart_vi` and `xbart` before it.
- [`bartcoreSamplerSetPredictor`](../../R/bartcore.R) ends `if (!forceUpdate) updateSuccessful else
  invisible(TRUE)`; the sampler method wraps it and keeps its visibility. The only caller that reads the
  value is the method; no other R code calls it.
- [`buildSamplerPriors`](../../R/bart.R) puts `power` and `base` into the quoted `cgm()` call by element assignment; a
  `NULL` removes the element and the next assignment is out of bounds. `bart` takes both through `...` as
  retired spellings (`consolidatedScalar`), `bartBT` as formals; every family door reaches the same builder.

## The rules

1. dec-B306. "Where n.threads is not stated and the cores cannot be counted, the fit runs on one thread with
   no message, at every function whose default is the guess (dbartsControl, bart, bartBT, rbart_vi, xbart); a
   stated n.threads that is NA stays refused, with dec-A163's message reworded to what a caller who wrote NA
   needs; guessNumCores() still returns NA; and the help for n.threads and for guessNumCores says the default
   is then one."
2. dec-B310. "The forced forms of setState, setPredictor and installTrees return NULL invisibly; the
   unforced forms return TRUE where the update was made and FALSE where it was not and the sampler is as it
   was. setPredictor's forceUpdate = "partial" keeps its vector by observation." Only setPredictor is built
   here; setState and installTrees have no forced form yet.
3. Root TODO, control-carried-count-edges: "bart(power = NULL) or bart(base = NULL) stops with 'subscript out
   of bounds' on any model, where a NULL is the default."
4. Root TODO, forest-selection-wordings: a list given forest by forest that carries a name it may not is told
   to select by position as `forest = 3`; the remedy must fit the argument the caller wrote. In the
   BayesTree-style help's Extracting Trees the sentence on the trailing forest margin must not cover `trees`,
   which return a data frame.

## Constraints

- R, tinytest and manual pages only; no C/C++ change; no draw of an accepted call moves.
- `guessNumCores()` still returns `NA`. A stated `NA` (and a stated 0 or negative) is refused as before.
- No subprocess in a test.
- `inst/NEWS.Rd` changes only where behaviour differs from the released 0.9-34; the root TODO is updated.

## Steps

1. The four formals keep their exported-name defaults. One internal helper, `fallBackToOneThread`, called from
   the bodies of `dbartsControl`, `rbart_vi` and `xbart` and, for `bart`, where it merges the defaults of the
   arguments it shares with the control, makes an unstated count whose default is `NA` one. A stated `NA` is
   refused with a message that tells its writer to leave the argument out. The help of `n.threads` at the
   four doors and of `guessNumCores` says the default is one where the cores cannot be counted.
2. The forced branch of `bartcoreSamplerSetPredictor` returns `invisible(NULL)`; its comment, the method's
   docstring and the Value of `setPredictor` in the sampler help say so.
3. `buildSamplerPriors` and `rbart_vi`, which builds its `cgm()` call by the same lines, take a `NULL` `power` or `base` as 2 and 0.95, the values `bart` and `bartBT`
   give them; checked at `bart`, `bartBT` and the family doors. A `NULL` still counts as named against a supplied
   `tree.prior` (`refuseColliding`'s rule, unchanged).
4. The selection refusal for a list given forest by forest names the argument written (`forests` of
   `setLeafPrior`, `bases` of `predict`); the Extracting Trees sentence is limited to the array it describes.

## Tests

- [test-excess-threads.R](../../inst/tinytest/test-excess-threads.R): with `guessNumCores` replaced in the
  package namespace for the length of the test, an in-process mock restored on exit, the four doors fit
  with `n.threads` unstated and no message, the control holds 1, and a stated `NA` is refused.
- [test-sampler-predictors.R](../../inst/tinytest/test-sampler-predictors.R): a forced `setPredictor`,
  whole matrix and column, returns `NULL` with `withVisible` showing it invisible; the unforced value stays
  `TRUE` or `FALSE`. Existing tests that read a forced `TRUE` are changed.
- [test-tombstones.R](../../inst/tinytest/test-tombstones.R), which covers the retired spellings: `power = NULL`, `base = NULL` equal the default fit, at `bart` and `bartBT`.
- [test-forest-selection.R](../../inst/tinytest/test-forest-selection.R): the new remedy text, at
  `setLeafPrior` and `predict`.

## Verification

Install into a private library; each new test is shown to fail with its fix reverted. Full tinytest suite
with `at_home = TRUE`, serial; the lint set of `docs/plans/README.md`; `R CMD check --as-cran` on a tarball
built from a clean `git archive`.

## Calls made in planning

- The thread default is one internal helper, not four inline `if`s, so the doors cannot drift apart.
- The fallback applies to the default only: a stated `n.threads` (including one computed by the caller from
  `guessNumCores()`) that is `NA` is refused.
- `NULL` `power` or `base` is read as the default in `buildSamplerPriors`, the one builder every door reaches,
  rather than at each door.
- No NEWS entry for item 3 or 4 where the released package behaved no differently; item 1 and 2 are decided
  in the report.
- The formals are not changed to a helper call, so the usage and `args()` show exported names only. The
  fallback sits in the bodies and reads `missing(n.threads)`; a wrapper that hands on its dots leaves the
  count unstated and gets one, while one that states a count of its own, `NA` included, is refused, as R's
  `missing()` leaves a wrapper's defaulted formal not missing.
- The orchestrator's calls: `cgm()` and `dart()` read a `NULL` `power` or `base` as their default, so
  the remedy the retired spelling's warning names works; on a model whose every forest has a basis, `bart(power =
  NULL)` and `bart(base = NULL)` stay refused as naming the plain forest's prior (dec-A179), pinned by a
  test; `forceUpdate = NA` at `setPredictor` and a control whose `n.threads` slot is `NA` given to `dbarts()`
  are refused by name.
