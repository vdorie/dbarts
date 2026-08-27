# Scoping: pre-RC documentation/record corrections (bartcore tip 9d0ee10f)

Read-only scope. Nine verified items plus a spot-check sweep. For each: file,
exact current text, proposed replacement (verbatim), risk. CI note per item.
Status: LANDED 2026-08-26 at 52c10e02 (design record 37ab6ea9; see the landing note); scoped at 9d0ee10f. Revised per an independent blind critique; all 15 adjudications
below applied - see the per-item notes for what changed and why.

---

## Cross-slice coordination (read first)

This scoping work touches `man/bart.Rd`, `man/rbart.Rd`, `man/xbart.Rd`, and
`man/dbartsSampler-class.Rd` - all four are ALSO edited by the concurrent
surface-refusals slice. This slice's edits should land AFTER surface-refusals
and be rebased onto its landed tip, not merged independently; a naive
parallel land risks silently clobbering whichever side lands second.

Specifically for item 2: surface-refusals homes the general whole-number
/fractional-count argument clause in `dbartsControl.Rd` plus the fitting
pages (`bart`/`bart2`/`rbart`/`xbart`) - not `dbarts.Rd`, the file item 2
edits. Item 2's proposed text below is therefore kept THIN and deferring
(a cross-reference to `dbartsControl.Rd`'s `n.samples` item for the full
contract, not a restated rationale) specifically so the two texts cannot
drift apart once surface-refusals lands its own version of that clause.

---

## 1. `print.bart`/`print.rbart`: no `\alias`, no `\usage`, no `\value`

**Cause.** `R/generics.R:3014-3018` and `:3020-3024` define
`print.bart`/`print.rbart` (both `{ printCall(x); fitSynopsis(x);
invisible(x) }`), registered in NAMESPACE:46/50. Neither has an `\alias`
anywhere in `man/`, so `?print.bart` fails, neither has a `\usage` entry, and
neither Rd's `\value` section mentions them. `man/bart2.Rd` shows the
sibling pattern done right: `print.bartMultinomial` has both
`\alias{print.bartMultinomial}` (bart2.Rd:6) AND a `\usage` entry
(`\method{print}{bartMultinomial}(x, \dots)`, bart2.Rd:88), sitting between
its `predict` and `residuals` usage blocks - every other documented method
on `bart.Rd`/`bart2.Rd` gets a `\usage` block, so alias-only would leave
`?print.bart` on a page whose `\usage` never shows it. Fixed here: alias
PLUS `\usage`, matching the sibling pattern.

**File: `man/bart.Rd`** (alias block, lines 1-8)

Current:
```
\name{bart}
\alias{bart}
\alias{plot.bart}
\alias{extract}
\alias{predict.bart}
\alias{extract.bart}
\alias{fitted.bart}
\alias{residuals.bart}
```

Proposed:
```
\name{bart}
\alias{bart}
\alias{plot.bart}
\alias{print.bart}
\alias{extract}
\alias{predict.bart}
\alias{extract.bart}
\alias{fitted.bart}
\alias{residuals.bart}
```

**File: `man/bart.Rd`**, `\usage` block: insert a `print.bart` entry between
the `predict.bart` block (ends line 54) and `extract(object, \dots)`
(line 56), matching bart2.Rd's predict-then-print adjacency.

Current (lines 53-56):
```
    n.threads,
    \dots)

extract(object, \dots)
```

Proposed:
```
    n.threads,
    \dots)

\method{print}{bart}(x, \dots)

extract(object, \dots)
```

**File: `man/bart.Rd`**, `\value` section, end (after line 340, before the
closing `}` at 341)

Current (line 340, last paragraph of `\value`):
```
  For continuous response fits, the plot method sets mfrow to c(1, 2) and makes two plots. ... [existing text, unchanged]
}
```

Proposed (new paragraph appended before the closing brace):
```
  For \code{print.bart}, the fit itself (\code{x}), returned invisibly; the call and a short synopsis (family, chain/tree/burn counts, kept-draw count) print to the console as the side effect.
}
```

Checked (critique note 3): this new paragraph lands immediately after the
line-340 `mfrow`/plotting paragraph, which was flagged as possibly
contradicting item 9's par()-restore fix. Read closely, line 340 only
describes what panels `plot` draws and that it sets `mfrow` to lay them out
- it makes no claim about restoring (or not restoring) the caller's
graphical parameters afterward, so it was never inaccurate and needs no
change; the new `print.bart` paragraph sits next to it without conflict.

**File: `man/rbart.Rd`** (alias block, lines 1-7): same treatment.

Current:
```
\name{rbart}
\alias{rbart_vi}
\alias{plot.rbart}
\alias{fitted.rbart}
\alias{extract.rbart}
\alias{predict.rbart}
\alias{residuals.rbart}
```

Proposed:
```
\name{rbart}
\alias{rbart_vi}
\alias{plot.rbart}
\alias{print.rbart}
\alias{fitted.rbart}
\alias{extract.rbart}
\alias{predict.rbart}
\alias{residuals.rbart}
```

**File: `man/rbart.Rd`**, `\usage` block: `rbart.Rd` uses `\S3method{...}`
(not `\method{...}`) throughout its own usage entries - matched here.
Insert between the `predict.rbart` block (ends line 63) and
`\S3method{residuals}{rbart}` (line 65).

Current (lines 63-65):
```
    \dots, group.by)

\S3method{residuals}{rbart}(object, type = "ev", \dots)
```

Proposed:
```
    \dots, group.by)

\S3method{print}{rbart}(x, \dots)

\S3method{residuals}{rbart}(object, type = "ev", \dots)
```

**File: `man/rbart.Rd`**, `\value` section, end (after line 155, before the
closing `}` at 156):

Proposed (new paragraph appended):
```
  For \code{print.rbart}, the fit itself (\code{x}), returned invisibly, printed the same way as \code{print.bart} (see \code{\link{bart}}).
}
```

**Risk.** man/ only - fires `R CMD check`/`check-standard` (missing-alias,
missing-\usage-entry, and missing/incomplete-\value are all check-standard
NOTEs), not `check-doc-freshness.R` (that script never walks `man/`). No
`_pkgdown.yml` change needed: an alias creates no pkgdown topic - pkgdown
indexes Rd FILES, and `print.bart`/`print.rbart` land on the existing
`bart`/`rbart` topics automatically. codoc exposure from the new `\usage`
entries is nil: both signatures (`x, \dots`) match the actual formals
exactly, so `R CMD check`'s codoc pass has nothing to flag. Also see the
cross-slice coordination note above - both files are shared with
surface-refusals.

---

## 2. `man/dbarts.Rd:54` - `n.samples` "positive integer" claim

**File: `man/dbarts.Rd`**, lines 53-55

Current (line 54):
```
  \item{n.samples}{
    A positive integer setting the default number of posterior samples to be returned for each run of the sampler. Can be overriden at run-time. See \code{\link{dbartsControl}}.
  }
```

**Cause.** `R/A_class.R:375-376`'s validity check refuses only
`object@n.samples < 0L` (not `<= 0`), so `0` is a valid value - and it is
given operational meaning at `man/dbartsControl.Rd:33`: "A non-negative
integer... `0` is accepted here (and by `dbarts`) - a sampler meant to be
driven by a host loop's own `run()` calls." (The length check at
`R/A_class.R:311-312` is a separate, unrelated constraint - length 1, not a
magnitude bound; that is the only other `n.samples`-adjacent check in
`A_class.R`, and it is not implicated here.) `dbarts.Rd`'s copy of the same
argument was never updated to match `dbartsControl.Rd`'s.

Proposed - kept thin per the cross-slice coordination note above (a
cross-reference, not a restated rationale, so it cannot drift from
`dbartsControl.Rd`'s own text once surface-refusals lands its version):
```
  \item{n.samples}{
    A non-negative integer (\code{0} accepted) setting the default number of posterior samples to be returned for each run of the sampler. Can be overriden at run-time. See \code{\link{dbartsControl}}'s \code{n.samples} item for the full contract.
  }
```

**Risk.** man/ only - fires check-standard. No `check-doc-freshness.R`
exposure. Sequence AFTER surface-refusals per the coordination note; if that
slice's `dbartsControl.Rd` wording changes materially, re-check this
cross-reference still resolves to the right item.

---

## 3. `xbart`'s `\value` understates the returned array shape

**File: `man/xbart.Rd`**, lines 129-135

Current:
```
\value{
  An array of dimensions \code{n.reps} \eqn{\times}{*} \code{length(n.trees)} \eqn{\times}{*} \code{length(k)} \eqn{\times}{*} \code{length(power)} \eqn{\times}{*} \code{length(base)}. If \code{drop} is \code{TRUE}, dimensions of length 1 are omitted. If all hyperparameters are of length 1, then the result will be a vector of length \code{n.reps}. When the result is an array, the \code{dimnames} of the result shall be set to the corresponding hyperparameters.

  For method \code{"k-fold"}, each element is an average across the \eqn{K} fits. For \code{"random subsample"}, each element represents a single fit.

  The result is a bare array with no class, so the fit generics - \code{predict}, \code{extract}, \code{fitted}, \code{residuals} - do not apply to it; it is a table of losses, not a fit.
}
```

**Cause** (`R/xbart.R:519-577`): the internal array is always allocated with
SIX dimensions - `rep, n.trees, k, power, base, loss` - not five, and which
ones survive to the returned object is not simply "drop iff length 1":
- `n.trees`, `power`, `base`: dropped iff `drop` and length 1 (matches the
  current text). `n.trees` is coerced to integer at `xbart.R:197`.
- `k`: dropped whenever `k` was NOT given as a numeric grid (i.e. resolved to
  a hyperprior; `kIsGrid <- is.numeric(k)`) - **regardless of `drop`**
  (`xbart.R:549`, `if (!kIsGrid) FALSE`; `kLength` is 1 off-grid,
  `xbart.R:437`). This differs from the other three: a hyperprior `k` never
  gets even a length-1 cell, unlike what `drop = FALSE` would suggest.
- `loss`: a trailing dimension, appended **only when the `loss` function's
  own return value has length > 1** (`numResults <- ncol(lossValues)`,
  `xbart.R:516`, `:558`) - independent of `drop` entirely (present with
  `drop = FALSE` too when the loss is scalar; the current text implies
  `drop = FALSE` keeps every dimension, which is false for this one). All
  three built-in losses (`xbart.R:606-633`) return scalars, so this axis
  never appears with a built-in loss.
- Collapse to a vector happens when zero of the above survive
  (`length(newDims) == 1`, `xbart.R:561`), not simply "all hyperparameters
  length 1" (loss counts too).
- `dimnames` (`xbart.R:566, 572-576`): only the double-valued axes (`k`,
  `power`, `base`) are rounded (`signif(x, 2L)`); `n.trees`, an integer, is
  `as.character`'d directly and so labeled exactly. `loss` and `rep` carry
  no per-slot names at all.

Proposed:
```
\value{
  An array with up to six dimensions, in order \code{rep} (length \code{n.reps}), \code{n.trees}, \code{k}, \code{power}, \code{base}, and \code{loss}. \code{rep} is always present. \code{n.trees}, \code{power}, and \code{base} are each omitted when \code{drop} is \code{TRUE} and the corresponding argument has length 1; with \code{drop = FALSE} they are always present. \code{k} follows a different rule: it is omitted whenever \code{k} was NOT given as a numeric grid (i.e. it resolved to a hyperprior; see the \code{k} argument above) REGARDLESS of \code{drop}, and otherwise follows the same drop-if-length-1 rule as the other three. The trailing \code{loss} dimension, sized to however many values a single call to \code{loss} returns, is present only when that count is greater than one - independent of \code{drop} entirely; the default losses and an ordinary scalar-returning custom \code{loss} never contribute it. When none of the above survive, the result collapses to a plain vector of length \code{n.reps}. When the result remains an array, its \code{dimnames} name the swept values on each surviving hyperparameter axis - exact integer labels for \code{n.trees}, values rounded to 2 significant digits for the double-valued \code{k}/\code{power}/\code{base} axes; the \code{loss} axis, when present, carries no per-slot names, and neither does \code{rep}.

  For method \code{"k-fold"}, each element is an average across the \eqn{K} fits. For \code{"random subsample"}, each element represents a single fit.

  The result is a bare array with no class, so the fit generics - \code{predict}, \code{extract}, \code{fitted}, \code{residuals} - do not apply to it; it is a table of losses, not a fit.
}
```

**Companion fix, now required** (the `\value` text above cross-references
"the `k` argument above" for why a hyperprior drops the axis; that
cross-reference must itself be correct, so this is no longer optional):

**File: `man/xbart.Rd`**, line 67 (the `k` `\item`)

Current:
```
     A vector of positive real numbers, setting the BART hyperparameter for the node-mean prior standard deviation. If \code{NULL}, the grid default of 2 is used for every response family. Binary responses do not inherit \code{bart2}'s Chi hyperprior default: a hyperprior is not a grid, so taking it here would leave the \code{k} axis a single cell. Hyperprior crossvalidation not possible at this time. A hyperprior \code{k} is held, not swept, and is DRAWN every sweep in every cell, so the reported loss is computed under a shrinkage that moves within each fit rather than under the named value. A numeric grid is always swept largest to smallest and the reported \code{k} axis un-permuted back to the order given, so results do not depend on the order \code{k} is listed in; cells still warm-start off the previous one, so this is order-invariance, not an unbiased estimate for each cell taken alone.
```

Proposed (fix only the "single cell" clause, rest unchanged):
```
     A vector of positive real numbers, setting the BART hyperparameter for the node-mean prior standard deviation. If \code{NULL}, the grid default of 2 is used for every response family. Binary responses do not inherit \code{bart2}'s Chi hyperprior default: a hyperprior is not a grid, so taking it here would drop the \code{k} axis from the result array entirely, rather than leave it a swept dimension (see \sQuote{Value}). Hyperprior crossvalidation not possible at this time. A hyperprior \code{k} is held, not swept, and is DRAWN every sweep in every cell, so the reported loss is computed under a shrinkage that moves within each fit rather than under the named value. A numeric grid is always swept largest to smallest and the reported \code{k} axis un-permuted back to the order given, so results do not depend on the order \code{k} is listed in; cells still warm-start off the previous one, so this is order-invariance, not an unbiased estimate for each cell taken alone.
```

**Risk.** man/ only, check-standard. No freshness-tool exposure. Two edits
in the same file now (the `\value` rewrite and the `k`-item wording fix)
must land together - landing only the `\value` rewrite would leave it
cross-referencing a sentence that still says the opposite thing.

---

## 4. `dbartsSampler$predict` return shape undocumented

**File: `man/dbartsSampler-class.Rd`**, `\value` section, the `predict`
paragraph (line 414)

Current:
```
  \code{predict} keeps the current test matrix in place and uses the current set of tree splits. This function has two use cases. The first is when \code{keepTrees} of \code{\link{dbartsControl}} is \code{TRUE}, in which case the sampler should be run to completion and the function can be used to interrogate the existing fit. When \code{keepTrees} is \code{FALSE}, the function can be used to obtain the likelihood as part of a proposed new set of covariates in a Metropolis-Hastings step in a full-Bayes sampler. This would typically be followed by a call to \code{setPredictor} if the step is accepted.
```

**Cause** (`src/R_interface_bartcore.cpp`, `predictFromSource` at 5669,
`bartcore_predict` at 5789): unlike `predictForests` (documented at Rd
line 430 with a full shape), the `predict` paragraph never states a return
shape at all - only behavior. Traced shape rules, all confirmed directly
against the source (`numSamples` at 5681-5682; multi-location branch
5715-5730; saved-tree branch 5731-5737; live-tree branch 5738-5744;
variance-list branch 5759-5772; `refuseEmptyTreeStore`,
`src/bartcore/data.hpp` area 2877-2883, makes the "at least one recorded
draw" gloss below safe - `keepTrees` alone with nothing run yet is refused,
not answered with an empty result):
- An amplitude-coupled (multi-forest, including a Bayesian causal forest)
  sampler refuses `predict` outright (`refuseUndefinedTestFits`,
  `src/R_interface_bartcore.cpp:5804`) - this is not a new fact: the SAME
  `\value` section already states it, one paragraph later, under
  `getFitsWithoutOffset` (line 426): "\code{predict} is refused" on such a
  sampler. The shape rules below apply only to a single-forest or
  multinomial sampler.
- Single-forest sampler, saved trees to replay (`keepTrees` true, capacity
  > 0): n.test x n.samples x n.chains, chain axis dropped at `n.chains == 1`
  (line 5731-5737) - mirrors `run()`'s own `train`/`test` channel shape,
  which line 408 of this same `\value` section already documents directly
  (n.obs x n.samples x n.chains).
- Single-forest, no saved trees (live/current trees replayed once per
  chain): NO samples axis at all - a plain vector of length n.test at one
  chain, an n.test x n.chains matrix otherwise (line 5738-5744).
- Multi-location (`family = "multinomial"`, `numLocations > 1`,
  `numReportedLocations() > 1` exists only on `MultinomialForestCombiner`,
  `src/bartcore/combiner.hpp:1598`): a K axis inserted between rows and
  samples - n.test x K x n.samples (x n.chains) - and here the samples axis
  is **never dropped**: it is sized 1 (not omitted) when there are no saved
  trees to replay (line 5715-5730 computes `numSamples` the same way
  regardless of `capacity`, unlike the single-forest branch). NOTE: this
  layout is stated outright below from the source directly, not by
  cross-reference to `run()`'s own entry - line 408 documents the K-widening
  explicitly only for `varcount`, and only says the SAME widening pattern
  "is what a multinomial sampler's per-category counts ride" as an aside;
  it does not itself spell out the `train`/`test` K-layout, so leaning on it
  for that claim would be circular.
- Heteroscedastic sampler (variance forest) with saved trees to replay:
  returns not a bare array but `list(mean = ..., variance = ...)`, both
  members the shape above (line 5759-5772); without saved trees to replay,
  there is nothing for the variance forest to replay, so the bare mean array
  returns alone, "backward-compatible" per the code's own comment
  (line 5757-5758).

Proposed (append to the existing paragraph, same location):
```
  \code{predict} keeps the current test matrix in place and uses the current set of tree splits. This function has two use cases. The first is when \code{keepTrees} of \code{\link{dbartsControl}} is \code{TRUE}, in which case the sampler should be run to completion and the function can be used to interrogate the existing fit. When \code{keepTrees} is \code{FALSE}, the function can be used to obtain the likelihood as part of a proposed new set of covariates in a Metropolis-Hastings step in a full-Bayes sampler. This would typically be followed by a call to \code{setPredictor} if the step is accepted.

  Its return shape mirrors \code{run()}'s own \code{train}/\code{test} channel for a single-forest sampler, with two exceptions. An amplitude-coupled (multi-forest, including a Bayesian causal forest) sampler refuses \code{predict} outright rather than answering here - see \code{getFitsWithoutOffset} below, which already records the same refusal - so everything that follows describes a single-forest or multinomial sampler only. For a single-forest sampler: with saved trees to replay (\code{keepTrees} \code{TRUE} and at least one recorded draw), the result is an n.test x n.samples x n.chains array, the trailing chain dimension dropped to a plain n.test x n.samples matrix when \code{n.chains} is \code{1}; without saved trees to replay, the CURRENT (live) trees are replayed once per chain instead, with no samples axis at all - a plain vector of length n.test at one chain, an n.test x n.chains matrix otherwise. A \code{family = "multinomial"} sampler's surface reports the K category probabilities instead: an n.test x K x n.samples (x n.chains) array, the K axis inserted between the rows and the samples. Unlike the single-forest case, this samples axis is never dropped: it is sized \code{1} rather than omitted when there are no saved trees to replay. A heteroscedastic (variance forest) single-forest sampler with saved trees to replay returns not a bare array but a named list \code{list(mean = ..., variance = ...)}, both members the single-forest shape above; without saved trees there is nothing for the variance forest to replay, so the bare mean array returns alone, exactly as a non-heteroscedastic sampler's would.
```

**Risk.** man/ only, check-standard. No freshness-tool exposure (this Rd is
not in `docs/design`). Shared with surface-refusals - see the cross-slice
coordination note above.

---

## 5. `rbart`'s `$fit`/`$n.chains` vs the "same elements as `bart`" claim

**Ruling: fix in code.** `$n.chains` should always be present, matching
`bart()`'s own unconditional behavior. Below is the code fix, the two test
sites it requires, and the `\value` text describing the FIXED behavior (no
absence clause - documenting the gap was rejected once the gap is closed
the same commit).

**Cause.** Traced `R/rbart.R`'s `packageRbartResults` (defined 1045,
called 418 and 539) against `bart.R:389-392`:
- `bart()`'s own `$fit`/`$n.chains`: `result$fit <- fit` (a SINGLE sampler
  object) only under `keepSampler`; `result$n.chains <- n.chains` is
  UNCONDITIONAL (bart.R:389-392) - `$n.chains` is always present.
- `rbart_vi()`'s `$fit`, when present, is ALWAYS A LIST, never a bare
  sampler - a different shape from `bart`'s (this part stays; it is
  intentional, not the bug). Two sub-cases: (a) the default,
  built-in-tau-prior/no-callback fast path (line 381, "in-core fast path")
  builds ONE multi-chain sampler; its own wrapper (line 428-432) sets
  `result$fit <- list(fitResult$sampler)` - a length-1 list wrapping one
  sampler that ran every chain; (b) the custom-prior/callback R-loop path
  (line 1293-1297) sets `result$fit <- lapply(chainResults, ...)` - a
  length-`n.chains` list, one sampler PER chain.
- `$n.chains` is present ONLY in the `else` branch of
  `packageRbartResults`'s own `if (keepSampler) {result$fit <- ...} else
  {result$n.chains <- n.chains}` (line 1293-1297) - i.e. only when `$fit`
  is absent. The fast path's wrapper (418-427) forces `keepSampler = FALSE`
  into the INTERNAL `packageRbartResults` call specifically to get
  `$n.chains` set, then bolts `$fit` on afterward (430-432) - so on that
  path BOTH are present when `keepSampler` is requested. On the
  custom-prior/callback path there is no such compensation: with
  `keepSampler = TRUE` there, `$n.chains` is simply ABSENT - an
  unintentional asymmetry between the two internal paths, not a documented
  design choice.
- Every reader already tolerates the absence
  (`fitNChains` `R/diagnostics.R:16-21`, `fitSynopsis` `R/generics.R:
  2958-2968`, `predict.rbart` `:2226-2230`, `extract.rbart` `:2435-2438,
  :2449`, `plotTree.rbart` `:2706-2709, :2722`) via `is.null(n.chains) ?
  length(fit) : n.chains`, and the two values already agree on the general
  path (`length(fit) == n.chains` there), so making `$n.chains`
  unconditional changes no reader's answer - it only removes the silent
  gap. `bartCause` reads only its own `$n.chains` (its own object, not
  rbart's), so it is unaffected either way.

**Proposed fix (code), three sites, land together:**

`R/rbart.R:1293-1297`

Current:
```r
  if (keepSampler) {
    result$fit <- lapply(chainResults, function(x) x$sampler)
  } else {
    result$n.chains <- n.chains
  }
```

Proposed:
```r
  if (keepSampler) {
    result$fit <- lapply(chainResults, function(x) x$sampler)
  }
  result$n.chains <- n.chains
```

`R/diagnostics.R:10-14` - retire the now-stale compensating parenthetical
(it singles out the in-core path as the one that "always keeps n.chains,"
which was true only because the general path did not; once both paths do,
the qualifier is misleading rather than merely redundant):

Current:
```r
# n.chains survives on the object whether or not the sampler was kept (see
# packageBartResults/packageRbartResults); fit is a single dbartsSampler for
# bart/bart2 and a list of them (length n.chains, or 1 for the in-core
# multi-chain rbart path, which always keeps n.chains on the object too) for
# rbart
```

Proposed:
```r
# n.chains survives on the object whether or not the sampler was kept (see
# packageBartResults/packageRbartResults); fit is a single dbartsSampler for
# bart/bart2 and a list of them for rbart (length n.chains on the general
# path, length 1 - one multi-chain sampler standing in for every chain - on
# the in-core path)
```

`inst/tinytest/test-rbart-bartcore.R:92` and `:125` - these currently PIN
the bug as expected behavior on the exact path the fix changes (a custom
`prior` function, default `keepTrees`/`keepSampler`, i.e. the
general/R-loop path) and must be updated in the same commit or the fix
breaks tinytest:

Current (line 92, `fit.custom` built with `n.chains = 2L`):
```r
expect_true(is.null(fit.custom$n.chains))
```
Proposed:
```r
expect_equal(fit.custom$n.chains, 2L)
```

Current (line 125, `fit.named` built with `n.chains = 1L`):
```r
expect_true(is.null(fit.named$n.chains))
```
Proposed:
```r
expect_equal(fit.named$n.chains, 1L)
```

**File: `man/rbart.Rd`**, `\value` section opening (line 139) - describes
the FIXED behavior only; no absence clause.

Current:
```
  An object of class \code{rbart}. Contains all of the same elements of an object of class \code{\link{bart}}, as well as the elements:
```

Proposed:
```
  An object of class \code{rbart}. Contains all of the same elements of an object of class \code{\link{bart}}, with one exception: \code{fit}, when present (\code{keepTrees}/\code{keepSampler} \code{TRUE}), is always a LIST rather than \code{bart}'s bare sampler object - a length-\code{n.chains} list of per-chain samplers on the general (custom \code{prior}/\code{callback}) fitting path, or a length-one list wrapping a single sampler that ran every chain together on the built-in-prior, no-\code{callback} path (the common case). \code{n.chains} is always present alongside it, whether or not \code{fit} is kept. \code{rbart_vi} also has the elements:
```

Note the spelling: `keepTrees`/`keepSampler` (camelCase), matching
`rbart_vi`'s own formals (`R/rbart.R:42,45`; `man/rbart.Rd:29,31`) - `bart`'s
lowercase `keeptrees`/`keepsampler` is a spelling split unique to `bart`'s
own page (`R/generics.R:270`), not shared by `rbart`.

**Risk.** `R/rbart.R` and `R/diagnostics.R` changes fire `check-standard`
and the C/C++-adjacent recompile discipline does not apply (pure R), but
`inst/tinytest` runs against the INSTALLED package - reinstall before
running `test-rbart-bartcore.R`. `man/rbart.Rd` is shared with
surface-refusals (cross-slice coordination note above). This does not
change any RNG draw (no sampling path touched), so no hardcoded-value
snapshot elsewhere should shift - confirmed no other test file asserts
`is.null(...$n.chains)` beyond the two sites above (`grep` swept every
`inst/tinytest/*.R` referencing `$n.chains`).

---

## 6. `man/xbart.Rd:99-101` missingness deferral - fold-view scope exclusion

**File: `man/xbart.Rd`**, lines 99-101

Current:
```
   \item{missing}{
     How missing values in the predictors enter the model; as in \code{\link{dbarts}}.
   }
```

**Cause.** `?dbarts`'s `missing` item (the contract this text defers to)
promises: "a column complete in training has no [missing] route, and an
\code{NA} there is refused, naming the column" - enforced by
`refuseTestMissingness` (`src/C_interface.cpp:254`), whose ONLY call site is
`validateTestSource` (`src/C_interface.cpp:277`, the "flat entrances" used
by `setTestPredictor`/`predict`; the R-side twin is `R/data.R:354`).
`xbart`'s internal per-fold sampler creation
(`R/xbart.R:687` -> `R/bartcore.R:687` -> `:699` ->
`C_dbarts_bartcore_createFromHandle` -> `bartcore_createFromHandle`,
`src/R_interface_bartcore.cpp:3575`) calls NEITHER function - the refusal
does not exist on this path at all, not a narrower version of it. Compounding
this: `ColumnStore::buildFromParent` (`src/bartcore/data.hpp:1392`) resets
`hasMissing` (`.assign`, line 1087) then recomputes it FRESH FOR EACH VIEW
from that view's own `rows` (the fold's TRAIN rows only, data.hpp:1459-1466)
- not inherited from the parent's dataset-wide flag. So a column that
carries `NA` somewhere in the full data, but happens to be complete within
one particular fold's training subset, gets NO learned missing-direction
rule for that fold; a held-out row of that same fold that IS `NA` on that
column is never refused - no refusal call exists on this path - and instead
silently takes the rule's unset default direction (left; `tree.hpp:138`'s
`setMissingGoesRight`, unset = left).

Proposed (one added sentence; reworded from the prior draft, which wrongly
implied a narrower REFUSAL still fires - the verified mechanism is that NO
refusal exists on this path at all, so a fold-complete column's held-out
`NA` row is never caught, the opposite of what "scoped to the fold" would
suggest):
```
   \item{missing}{
     How missing values in the predictors enter the model; as in \code{\link{dbarts}}, with one exception: \code{xbart}'s internal per-fold sampler construction does not go through the entry points that check it, so the test-\code{NA}-on-a-training-complete-column refusal described there does not apply here - a held-out row that is \code{NA} on a column complete within its own fold's training rows is never refused; it silently takes that rule's default (left) branch instead.
   }
```

**Risk.** man/ only, check-standard. No freshness-tool exposure. Shared with
surface-refusals (`xbart.Rd`) - see the cross-slice coordination note above.

---

## 7. `benchmarks/R/composition-matrix.R` - the two recorded disagreements

Both are pre-existing (per the composition-refusals landing note,
`docs/plans/composition-refusals.md:642-644`), unrelated to that slice.
Investigated cause and a proposed code fix for each (harness bugs, not real
model-composition regressions) - both fixes, and the landing-note mark
below, must land in the SAME commit, or the committed record keeps
asserting a gate finding that no longer reproduces.

### 7a. `multinom dbarts5` - "no base fixture recipe"

**Cause.** `docs/design/feature-matrix.md:85` claims `(multinom, dbarts5) =
S` (`dbarts.R:381` - the `family = "multinomial"` token on `dbarts()`'s own
`family=` formal, a real, shipped, direct R5 construction route -
`dbarts(x, factor(y), family = "multinomial")` is a real route,
`R/dbarts.R:584-592`). But `table1Probes$dbarts5`
(composition-matrix.R:545) is a BLANKET dispatcher:

```r
table1Probes <- list(
  bart = probeBart1,
  bart2 = probeBart2,
  dbarts5 = function(family, seed) buildBase(family, seed),
  rbart_vi = probeRbart,
  xbart = probeXbart,
  flatc = probeFlatC
)
```

`buildBase`'s own switch (composition-matrix.R:210-296) DELIBERATELY has no
`multinom`/`hurdle` arm - its header comment says so explicitly ("Everything
but multinom/hurdle: those two host-shell families ... are special-cased
directly in their probes instead of forcing a 13th shape into this
recipe") - and falls through to:
```r
    stop("composition-matrix: no base fixture recipe for '", family, "'")
```
Every OTHER `table1Probes` entry (`bart2`, `xbart`) already special-cases
`multinom` via `bart2Args` (composition-matrix.R:385: `multinom = list(d$x,
factor(d$label), family = "multinomial")`); only `dbarts5`'s one-line
wrapper never got the same routing-around that `runProbe` already gives
`multinomActiveRows`/`hurdleDart` for the OTHER two special-cased cells
(composition-matrix.R:644-649). This is a harness gap, not a genuine model
refusal - `attempt()` classifies the generic `stop()` as REFUSES (nothing in
`genericFailure` matches "no base fixture recipe"), so it reads as a
plausible-looking refusal in the report rather than what it is.

**Proposed fix** (composition-matrix.R, two spots):

Add a helper near `multinomActiveRows`/`hurdleDart` (after line 561), with
the same `invisible(sampler$run(0L, 1L))` sweep `buildBase` performs on
every other family (composition-matrix.R:301), for parity:
```r
multinomBase <- function(seed) {
  d <- mkXY(seed)
  sampler <- dbarts(d$x, factor(d$label), family = "multinomial", control = ctl(seed))
  invisible(sampler$run(0L, 1L))
  sampler
}
```

Add a special case in `runProbe`, alongside the other two (before the
`capability %in% names(table1Probes)` dispatch, composition-matrix.R:643-650):
```r
runProbe <- function(family, capability, seed) {
  if (family == "multinom" && capability == "activeRowsMask") {
    return(attempt(multinomActiveRows(seed)))
  }
  if (family == "hurdle" && capability == "dart") {
    return(attempt(hurdleDart(seed)))
  }
  if (family == "multinom" && capability == "dbarts5") {
    return(attempt(multinomBase(seed)))
  }
  if (capability %in% names(table1Probes)) {
    return(attempt(table1Probes[[capability]](family, seed)))
  }
  ...
```

### 7b. `logistic setWeights` - integer-weights

**Cause.** `docs/design/feature-matrix.md:137` (`[f10]` at `:315`, which
already states the count semantics) claims `(logistic, setWeights) = S`
(`MOD:3600`). The shared `mutate$setWeights` probe
(composition-matrix.R:620):
```r
  setWeights = function(s, d) s$setWeights(runif(length(s$data@y), 0.5, 1.5)),
```
draws continuous, virtually-never-integer weights - fine for gaussian and
every other weighted family, but `logistic`'s `setWeights` specifically
requires positive integers by MODEL design (weights are observation-count
replication under the PG(w, psi) augmentation), enforced at
`enforceBinaryWeightPolicy` (`src/R_interface_bartcore.cpp:2756-2761`, from
the call site at `:4902`):
```cpp
  if (family == bartcore::ResponseFamily::logistic)
    for (size_t i = 0; i < numObservations; ++i)
      if (!(weights[i] > 0.0) || weights[i] != std::floor(weights[i]))
        Rf_error("logistic weights are observation counts and must be "
                 "positive integers; drop zero-count rows, and use a gaussian "
                 "model for continuous weights");
```
`runProbe`'s generic dispatch (composition-matrix.R:686-690) already passes
`d$family <- family` into the closure before calling the probe, so the fix
is entirely local - the data needed to special-case logistic is already
there, unused.

**Proposed fix** (composition-matrix.R:620, in place):

Current:
```r
  setWeights = function(s, d) s$setWeights(runif(length(s$data@y), 0.5, 1.5)),
```

Proposed:
```r
  setWeights = function(s, d) {
    w <- if (identical(d$family, "logistic")) {
      sample(1:3, length(s$data@y), replace = TRUE)
    } else {
      runif(length(s$data@y), 0.5, 1.5)
    }
    s$setWeights(w)
  },
```

**Companion mark, required, same commit:**

**File: `docs/plans/composition-refusals.md`**, line 644 only (642-643
untouched verbatim) - this paragraph is the committed record that names
both disagreements by exact wording ("multinom dbarts5 'no base fixture
recipe'; logistic setWeights integer-weights"); landing the harness fix
without marking it leaves a committed record asserting a gate finding that
no longer reproduces, for the human RC review to re-derive from scratch.
Checked for line-number citers first (per house practice for this kind of
edit): `grep -rn "composition-refusals.md" docs/ TODO` finds only
`docs/plans/prerc-surface-freeze.md:7` and `docs/plans/INDEX.md:205`, and
NEITHER cites it by line number (both are bare filename mentions) - so
strictly nothing requires line-count invariance here, but the edit is kept
line-count-invariant anyway (single line, in place) as the lower-risk
default.

Current (line 644):
```
dbarts5 "no base fixture recipe"; logistic setWeights integer-weights); multinomial-mutation-arc.md sections 5-8 and
```

Proposed:
```
dbarts5 "no base fixture recipe"; logistic setWeights integer-weights - fixed at HEAD, both harness bugs, not model gaps); multinomial-mutation-arc.md sections 5-8 and
```

**Risk.** `benchmarks/` is not `R CMD check`-covered and not walked by
`check-doc-freshness.R`; it is its own gate, run manually
(`Rscript benchmarks/R/composition-matrix.R`). The `docs/plans/` mark fires
nothing either (not `docs/design`, so Part 3 never walks it; no hash token
added, so Part 4 is not implicated). Together the two edits remove two false
"DISAGREEMENTS" from that gate's report and close the paper trail on them in
the same commit, so neither is re-derived from scratch by the RC human
review.

---

## 8. `TODO`'s `rc-gate` entry claims completeness while D7 is open

**File: `TODO`**, lines 188-189

Current:
```
  rc-gate: the pre-RC slate is complete; the RC declaration itself is
    VD-held.
```

**Cause.** `docs/plans/prerc-surface-freeze.md:3-7` (Status: DECIDED
2026-08-25): "D6 and D8 LANDED 936825d7 (docs/plans/composition-refusals.md);
D7 remains, at the RC tip." `D7` is `bcf-baseline-cross-host`, still an open
item in TODO itself (lines 33-36: "exempt the BCF snapshot channels under a
cross-host flag and re-record bcf-equivalence at the RC tip..."). The
`rc-gate` entry's "the pre-RC slate is complete" is accordingly false at the
current tip (9d0ee10f).

Proposed (same 2 physical lines, so nothing after line 189 in TODO shifts):
```
  rc-gate: the pre-RC slate is complete except bcf-baseline-cross-host,
    still open at the RC tip; the RC declaration itself is VD-held.
```

**Risk.** TODO edits fire nothing in CI directly, but `TODO` IS walked by
`tools/check-doc-freshness.R`: Part 4 (hash resolvability) scans it for
commit-hash-looking tokens (the proposed text adds none, so no exposure
there). More importantly, several `docs/design/*.md` files cite `TODO` by
exact line number AFTER this entry - re-derived directly by grep (correcting
the file attributions from an earlier pass of this scoping, which had two
wrong):
- `docs/design/multinomial-mutation-arc.md:812` -> `TODO:309-322`
- `docs/design/tree-mixing-proposals.md:975` -> `TODO:302`
- `docs/design/tree-mixing-proposals.md:2067` -> `TODO:312`
- `docs/plans/archive/multiforest-extension-surface.md:71` -> `TODO:190-214` (one
  line past this entry; itself hash-qualified "at 934a02d5", so already
  self-marked as a snapshot reference rather than a live one)

(`TODO:54-57`, `:150-168`, `:162-168` - all cited from
`multinomial-mutation-arc.md` - sit BEFORE line 188 and are unaffected by an
in-place edit at 188-189 regardless of line count.) `tools/
check-doc-freshness.R` Part 3 walks every `docs/design/*.md` file's
`file:line` anchors, including bare `TODO:N` ones, and would FAIL if any of
the four post-188 citations above moved. The proposed edit is deliberately
line-count-invariant (2 lines in, 2 lines out) for exactly this reason -
confirmed by re-deriving the citer list and the line-count math directly,
not by executing the script (read-only scope). Also note (Part 1): TODO
itself is not index-checked by name, only `docs/design`/`docs/plans`
`INDEX.md` entries are.

---

## 9. `docs/plans/release-candidate-review.md` - stale plot par() claims

**Background.** Commit `fcbbc478` ("Refuse silently-wrong formula and
argument inputs at fit entry points", 2026-08-25, IS an ancestor of the RC
tip 9d0ee10f) fixed exactly this: "plot.bart and plot.rbart changed
par(mfrow) to lay out the sigma-trace and interval panels without saving or
restoring the caller's graphical parameters; both now save and restore
around the [plot]." Two passages in `release-candidate-review.md` still
assert stale claims about this - and, per re-examination below, the two
passages are not actually parallel: 9a's claim was accurate when written and
later became stale; 9b's claim was already overbroad when written.

**Citation check** (per house practice, before proposing an edit that could
change line counts): `grep -rn "release-candidate-review.md" docs/ TODO`
turns up citations by exact line number from
`docs/plans/review-2026-08-24/memos/prerc-lens2-backlog.md:29`
("release-candidate-review.md:846-850" - close to, though already offset
from, the second passage below; likely already-stale from unrelated earlier
edits, not something this edit should try to fix) and
`threaded-predict-critique.md:119,248` ("release-candidate-review.md:
2529-2531", ":2529", well after both target passages). No `docs/design/*.md`
file cites this doc at all (confirmed empty grep), so
`tools/check-doc-freshness.R` never walks or validates any citation into it
(Part 3 only walks `docs/design/*.md`) - this file's line-number citations
are a human-consistency concern only, not an automated-gate one. Given
existing citations after both target passages, **both proposed edits below
are single-line, in-place edits that add zero net lines**, to avoid adding a
new failure mode to those other docs regardless.

**Note on style.** `docs/design/model-space-survey.md` section 6 marks a
fixed-at-HEAD passage by PREPENDING a new, clearly-labeled sentence ("Fixed
at HEAD (not a live defect): ...") before the original defect description,
leaving the original wording 100% untouched below it - marking, not
rewriting. That pattern adds lines. Given the citation risk above, the
edits below adapt the same spirit (append a marker, minimize touching
original wording) but keep it to appended/adjusted clauses on the SAME line
rather than a new prepended sentence, to stay line-count-invariant. If VD
prefers the purer model-space-survey form (a new prepended line) instead,
that is a reasonable alternative - see the alternative below - but it
requires re-verifying/re-anchoring `prerc-lens2-backlog.md:29`'s `:846-850`
citation (and confirming no other citation beyond line 610 elsewhere) as a
coordinated follow-up, the same kind of pass commit `936825d7` did for
`docs/design` after the composition-refusals landing.

### 9a. Lines 609-610

Current:
```
Doors left (unchanged unless noted): fit-time test-basis; logistic
weights never ride saved state; DESCRIPTION Date bumps at the RC; the 70
freshness advisories; plot.bart/plot.rbart still leak par() through
plotSigmaTrace (the six own-class sites restore).
```

This claim was ACCURATE when written (own-class sites already restored;
plot.bart/plot.rbart did not, until fcbbc478) and only later became stale -
so it is marked, not reworded.

Proposed (edit only line 610; lines 607-609 untouched verbatim):
```
Doors left (unchanged unless noted): fit-time test-basis; logistic
weights never ride saved state; DESCRIPTION Date bumps at the RC; the 70
freshness advisories; plot.bart/plot.rbart still leak par() through
plotSigmaTrace (the six own-class sites restore; plot.bart/plot.rbart fixed fcbbc478).
```

### 9b. Lines 839-840 (task's "line 842 vicinity" - same paragraph, `residue`
list at 837-848)

Current (lines 839-840; 841 unchanged, shown for context):
```
one run(n.burn, n.samples), so plot draws no burn-in segment); par(mfrow) is not
restored by any plot method (plot.bart does not either - a shared change);
as_draws_array.bartMultinomial(vars =) ignores a non-meanProb value (documented);
```

Re-examined per critique: this claim's LEADING clause - "par(mfrow) is not
restored by ANY plot method" - was already partly false when written, not
merely later-outdated: line 610 (9a, above, in this same file) already
records that the six own-class sites restored par() at that point in the
document's own timeline. The prior draft of this fix only appended a marker
inside the trailing parenthetical ("plot.bart does not either"), leaving
the false leading claim itself unmarked. Corrected to mark the leading
claim, with precise attribution: overbroad when written (six sites already
restored, recorded at line 610), and the true remaining half (plot.bart/
plot.rbart) closed by fcbbc478.

Proposed (edit only line 840; line 839 and 841 untouched verbatim):
```
one run(n.burn, n.samples), so plot draws no burn-in segment); par(mfrow) is not
restored by any plot method (overbroad when written - the six own-class sites already restored, line 610; plot.bart/plot.rbart's half fixed fcbbc478);
as_draws_array.bartMultinomial(vars =) ignores a non-meanProb value (documented);
```

**Alternative (not proposed, for completeness): true model-space-survey
form.** Insert, immediately before the "Doors left (unchanged unless
noted):" paragraph (i.e. before line 607) and again before the "Doors and
residue recorded here, not in code:" paragraph (before line 837), a
standalone marker sentence such as: "Fixed at HEAD (fcbbc478, an ancestor of
the RC tip): plot.bart and plot.rbart now save and restore \code{par()}
around their plotting, closing the gap the next passage records as still
open (see also line 610's own record that the other six sites already
did)." This is closer in spirit to the design-doc convention (100%
additive, original wording untouched) but adds 2 lines total before line
610 and before line 840 respectively, shifting every line number after each
insertion point - which would require re-verifying (and likely re-stamping)
every downstream numeric citation into this file listed above. Recommend the
in-place, single-line version (9a/9b) unless VD specifically wants the
line-shift and the coordinated re-anchor that follows it.

**Risk.** docs/ only - fires nothing in CI (not `check-doc-freshness.R`, as
established above; not `check-standard`). Both proposed edits are
single-line, additive-only, zero net line change.

---

## Sweep: other `\value` vs. code, on predict-surface / composition-refusals
surfaces (spot-check, not exhaustive)

Checked, found already consistent (no action needed):
- `man/dbarts.Rd`'s `\item{variance}` (line 81) already states the
  grouped-random-effects-with-variance-forest refusal fe0b3292 added at
  `R/spec.R:532` ("resid.dist = student() residuals and a grouped
  (rbart_vi) fit are also refused together with variance") - matches the
  current code exactly; no man/ update needed for that landing, despite
  `git show fe0b3292 --stat` touching no `man/` file (the wording predates
  or was folded into an earlier pass; confirmed current, not stale).
- `man/rbart.Rd`: no `value=` argument was ever documented for
  `predict.rbart`, so `befc8f45`'s removal of the `value=`/`"post-mean"`
  deprecation shims (a pure `R/generics.R` change, no `man/` diff) left
  nothing stale to fix.
- `man/bart.Rd`/`man/bart2.Rd`/`man/rbart.Rd`'s `offset`/`offset.test`
  wording (predict-surface's "one offset spelling" unification, `7b3ac6bf`):
  internally consistent - `offset.test` (fit-time formal) vs. `predict`'s
  own `offset` (new-data-time) are each named and distinguished explicitly
  in `bart.Rd:189`; no stray `offset.test` references found where `predict`'s
  `offset` should be.
- `man/bart.Rd:340`'s `plot` `mfrow` paragraph (checked again for item 1's
  benefit, see above): describes panel layout only, makes no
  restore/non-restore claim, so it is not stale and needs no par()-related
  edit.

Not independently re-verified beyond the above (explicitly out of scope per
"spot-check, do not boil the ocean"): the full `man/dbartsSampler-class.Rd`
and `man/bart2.Rd` diffs under `71cc7133`/`7b3ac6bf` themselves (63 and 12
line changes respectively) were read for shape but not re-derived line by
line against current `src/`; items 1-6 above already found the concrete
bugs on this surface (dbartsSampler's bare `predict` shape, xbart's `\value`
and its `k`-item cross-reference, rbart's `$fit`/`$n.chains`) that a fuller
sweep would likely have surfaced anyway.

---

## Summary table

| # | File(s) | Kind | CI on edit | Line-count sensitive |
|---|---|---|---|---|
| 1 | man/bart.Rd, man/rbart.Rd | alias + usage + \value | check-standard | no |
| 2 | man/dbarts.Rd | \value/\arguments wording (thin, defers to dbartsControl.Rd) | check-standard | no |
| 3 | man/xbart.Rd | \value rewrite + k-item wording fix (2 sites, same file) | check-standard | no |
| 4 | man/dbartsSampler-class.Rd | \value addition | check-standard | no |
| 5 | R/rbart.R, R/diagnostics.R, inst/tinytest/test-rbart-bartcore.R, man/rbart.Rd | code fix (3 sites) + \value wording | check-standard + tinytest | no |
| 6 | man/xbart.Rd | \arguments, one sentence | check-standard | no |
| 7 | benchmarks/R/composition-matrix.R, docs/plans/composition-refusals.md | code fix (2 spots) + landing-note mark, same commit | none (manual gate) | no |
| 8 | TODO | wording, line-count-invariant | none directly; TODO is freshness-walked via docs/design citations | YES - kept invariant |
| 9 | docs/plans/release-candidate-review.md | marking, line-count-invariant | none | YES - kept invariant |

Cross-slice: items 1, 2, 4, 5, 6 touch `man/` files also edited by the
concurrent surface-refusals slice (`bart.Rd`, `rbart.Rd`, `xbart.Rd`,
`dbartsSampler-class.Rd`) - land this work after surface-refusals, rebased
on its tip; see the coordination section at the top.

All proposed text is ASCII, carries no attribution, and (items 1-8) no
docs/ path or slice codename inside a shipped file - item 7's fix lives in
`benchmarks/` and `docs/plans/`, neither a shipped directory, so existing
docs/plans citations in comments are fine as-is and untouched by the
proposed fix.

## Landing note (2026-08-26)

LANDED at 52c10e02 (design record 37ab6ea9), implemented on the
surface-refusals slice d48aef8a as its cross-slice section orders and
cherry-picked from its gated worktree build; both gate batteries
green on that base: tinytest 7458/0, equivalence trio bitwise
43/12/11, composition-matrix disagreements 2 -> 0 with confirmations
181 -> 183 and the multinom dbarts5 cell now probing construction
plus one sweep, R CMD check --as-cran 1 NOTE with codoc/usage OK over
the new print aliases, rc-codoc green, freshness byte-identical to
its base (17 FAIL / 77 WARN, all pre-existing stacked-slice drift),
line counts invariant on TODO (347), composition-refusals.md (649),
release-candidate-review.md (3894), multinomial-mutation-arc.md
(1258). One design gap found and closed at implementation: item 1's
man/bart.Rd insertions shift the Saving subsection by +3 lines and
docs/design/multinomial-mutation-arc.md:233 cites it by exact line -
the design's risk note covered man/ only as a freshness source, not
as an anchor target - re-anchored :251 -> :254 in the same commit,
single token, line-count neutral, outside that doc's frozen sections.
The rbart n.chains fix is the unconditional assignment in
packageRbartResults with the R/diagnostics.R comment retired and the
two pinned tests flipped to expect the count.
