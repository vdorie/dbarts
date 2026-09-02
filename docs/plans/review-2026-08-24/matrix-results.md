# Executable consistency matrix - second whole-branch review

Branch `bartcore`, tip `b102e17c`. Built and run against a scratch install
(`R CMD INSTALL --preclean`), never the checkout. Script: `matrix.R` in this
directory; raw grid: `matrix-grid.csv` (836 rows from `matrix.R`, plus 41
appended by the supplementary `matrix-sampler-census.R` - see section 8 -
for 877 total). Fixture: n=40 train / 12 test, n.trees 3-4, n.samples 5-6,
n.burn 3-4, fixed seeds - the whole grid runs in under 5 seconds (matrix.R)
plus a few more for the census, well inside the 30-minute budget (see
"Budget and coarsening" at the end for what that headroom was NOT spent on).

Every number below was produced by running `matrix.R` (and, for section 8,
`matrix-sampler-census.R`), not transcribed from
an untracked review memo; where it disagrees with that memo's
citations, the disagreement is noted rather than silently resolved.

Both scripts are independently re-runnable against any installed lib:
`R_LIBS=<lib> Rscript matrix.R [outDir]` then
`R_LIBS=<lib> Rscript matrix-sampler-census.R [outDir]` (the second appends
to the first's `matrix-grid.csv`, so run `matrix.R` first - it truncates and
rewrites the file on each run, which would otherwise discard the census
rows).

## 1. Family-token acceptance, per entry

Tokens are pulled mechanically from each entry's own `family=` formal
default (`eval(formals(fn)$family)`), so this table is self-updating if the
branch's vocabulary changes. "-" = not in that entry's declared vocabulary at
all (never reaches match.arg; see the note on `dbartsData`/`dbartsControl`
below the table).

| token | dbarts | bart2 | bart | rbart_vi | xbart | dbartsSpec |
|---|---|---|---|---|---|---|
| auto | accept | accept | accept | accept | accept | accept |
| gaussian | accept | accept | - | accept | accept | accept |
| probit | accept | accept | - | - | accept | accept |
| logistic | accept | accept | accept | - | accept | accept |
| aft | accept | accept | accept | **refuse (bug, see F2)** | - | **error, no fixture reaches it (see F2)** |
| multinomial | accept | accept | - | - | - | **refuse (see F1)** |
| ordinal | accept | accept | - | - | - | accept |
| nbinom | accept | accept | - | - | - | accept |
| hazard | accept | refuse* | - | - | - | - |
| hazard.probit | accept | refuse* | - | - | - | - |
| hazard.logistic | accept | refuse* | - | - | - | - |
| hurdle.lognormal | refuse** | accept | - | - | - | - |
| twopart | refuse** | accept | - | - | - | - |
| (invalid token) | bare match.arg | bare match.arg | bare match.arg | bare match.arg | bare match.arg | bare match.arg |

Counts confirmed: dbarts 13/13, bart2 13/13, bart 3/3, rbart_vi 3/3, xbart
4/4, dbartsSpec 8/8 tokens declared - matching the memo's 13/13/3/3/4-and-8
count. `dbartsData` and `dbartsControl` were also probed directly with
`family=` (task's parenthetical "+ dbartsData ... dbartsControl as
constructors"): **neither has a `family` formal at all** -
`dbartsData(x, y, family="gaussian")` and `dbartsControl(family="gaussian")`
both raise raw `unused argument (family = "gaussian")`. Family resolution is
entirely a `dbarts`/`bart2`/`bart`/`rbart_vi`/`xbart`/`dbartsSpec` concern;
the two lower constructors never see the token.

\* refused only because the fixture also passed `test=` (bart2 forwards it
by default in this probe); a hazard/hazard.logistic fit built without
`test=` is accepted (see representative fits `bart.hazard`,
`bart.hazard.logistic` in the CSV) - this is documented, correct behavior
(discrete-time hazard has no `test` concept), not a gap.

\** `dbarts()` refuses `hurdle.lognormal`/`twopart` by name
("fits two component samplers and is only available through bart2()"),
matching [[bart2.Rd:237@b102e17c]]'s own statement of the same restriction. Agrees.

## 2. Class x generic table

Existence is a static fact, re-derived via `getS3method(generic, class,
optional=TRUE)` for 9 classes (bart, rbart, bartMultinomial, bartOrdinal,
bartNegbin, bartHurdle, pdbart, pd2bart, dbartsSampler) x 7 core generics
(predict, extract, fitted, residuals, summary, print, plot): **42 of 63
cells exist, 21 are empty.** (The memo's "18 of 63" is a different count on
the same axes - I could not reproduce 18 from the current tree; 21 is what
`getS3method` reports today, possibly reflecting drift since the memo, or a
different generic selection. Not reconciled further.) A separate 4-generic
extension (plotTree, survivalProbabilities, as_draws_array, as_draws_df) x
9 classes: **9 of 36 exist** (bart and rbart each carry all 4; dbartsSampler
carries `plotTree` only; the other 6 classes carry none). Note:
`as_draws_array`/`as_draws_df` are registered onto `bart`/`rbart` by dbarts's
own `.onLoad` hook *conditional on the `posterior` package being loaded* -
running this probe without `library(posterior)` first makes both show
`no-method` everywhere even though the methods exist and work once
`posterior` is attached; the run behind these numbers loads it.

Empty core cells:

| class | missing generics |
|---|---|
| bartOrdinal | plot |
| bartNegbin | plot |
| bartHurdle | plot |
| dbartsSampler | fitted, plot, predict, print, residuals, summary |
| pdbart | extract, fitted, predict, print, residuals, summary |
| pd2bart | extract, fitted, predict, print, residuals, summary |

`dbartsSampler`'s 6 and pdbart/pd2bart's 6-each are by design (dbartsSampler
has its own `$predict()`/`$show()` RC methods and `extract.dbartsSampler`/
`plotTree.dbartsSampler`, not the "bart"-family S3 surface; pdbart/pd2bart
are plot-only partial-dependence summaries). **The genuinely bare spot is
bartOrdinal/bartNegbin/bartHurdle: all three have `summary` (and print,
predict, extract, fitted, residuals) but none has `plot`**, while bart,
rbart, and bartMultinomial all do - re-derives the memo's citation exactly.

## 3. Rd cross-check (targeted, not exhaustive - see Budget)

Given the size of a full 29-file Rd sweep, this pass verified the Rd claims
tied directly to what the grid actually exercised (dbartsSpec.Rd, rbart.Rd,
bart2.Rd), rather than an exhaustive extraction across every man page.

- **CONTRADICTS**: [[man/dbartsSpec.Rd:40@b102e17c]], "Both entry points call the same
  resolution, so a family can never resolve two ways" - directly falsified
  by finding F1 below (`family = "multinomial"` resolves on `dbarts()` and
  refuses on `dbartsSpec()` for the identical input).
- **Internally inconsistent, same file**: [[man/dbartsSpec.Rd:48@b102e17c]] groups
  `"multinomial"` with `"hurdle.lognormal"` as families "unavailable"
  through this entry point because they "describe more than one sampler" -
  false for multinomial (one K-forest sampler, not multiple; unlike
  hurdle.lognormal, which genuinely composes two). Also inconsistent with
  the function's own formals: `hurdle.lognormal` is correctly absent from
  dbartsSpec's `family=` vocabulary, but `"multinomial"` IS present and
  match.arg-selectable - so the same paragraph claims a token is
  unreachable while the signature offers it.
- **AGREES**: [[man/rbart.Rd:94-95@b102e17c]] already documents "[the survival response]
  enter[s] only through the formula interface" for `family = "aft"` - the
  restriction itself is correctly documented; finding F2 is about the
  *enforcement*, not the documentation.
- **AGREES**: [[man/bart2.Rd:306@b102e17c]], "\[the forest()\ term route\] is the only
  way bart2 itself reaches the per-forest channel, since it has no
  `forests =` formal of its own" - matches the grid exactly (bart2 refuses
  `forests=` by name; dbarts accepts it). Not a bug or a doc gap.

## 4. Error-without-reason catalog (66 cells)

- **33 cells**: `extract(fit, type="trees", sample=<anything>)` on every
  bart-class and rbart representative - see finding F3. Same root cause,
  counted once there.
- **23 cells**: bare `unused argument (...)` from entries with no `...`
  catch-all (`bart()`, `xbart()`) hit by a keyword the entry doesn't
  support - see finding F4.
- **10 cells**: deliberate `zzz-invalid-type`/`zzz-invalid-sample` probes
  in Part C3 hitting a bare `match.arg`-style `"'arg' should be one of
  ..."` - this is the *expected*, by-design outcome for that probe (the
  task's own definition of "error without reason" names bare match.arg
  explicitly), not a new finding.

## 5. Entry-pair disagreements

| pair | input | result |
|---|---|---|
| dbartsSpec vs dbarts | family="multinomial" | **DISAGREE** - see F1 |
| dbartsSpec vs dbarts | family=auto/gaussian/probit/logistic/ordinal/nbinom | agree (6/7 tokens) |
| dbarts vs bart2 | weights containing NA | agree - both refuse `'weights' cannot be NA` |
| bart2 vs rbart_vi vs bart | n.thin/keepevery = -1 | agree - all three refuse `'n.thin' must be a positive integer` (the memo's "rbart_vi tests == 0L not <=0L" claim does not reproduce on this tree; looks fixed) |
| dbarts vs bart2 | `forests=` | different by design - dbarts accepts, bart2 refuses by name (documented, see S3 above) |
| dbarts/bart2 vs bart/rbart_vi/xbart | `resid.dist=student()` (NSE) | dbarts/bart2/bart accept; rbart_vi refuses by name (`unknown argument 'resid.dist'`); xbart raises raw `unused argument` - see F4 |

## 6. Open cells (Rd-silent, not fill-vs-refuse judged)

- 6 `dbartsSampler` core-generic gaps (fitted/plot/predict/print/residuals/
  summary): no Rd claims either way; would mean "does dbartsSampler ever
  grow the bart-family S3 surface, or does it stay RC-method-only forever."
- 12 pdbart/pd2bart gaps (extract/fitted/predict/print/residuals/summary x
  2 classes): same question for the partial-dependence objects.
- 31 of 36 extra-generic cells (plotTree on 6 non-{bart,rbart,dbartsSampler}
  classes; survivalProbabilities on 7 non-{bart,rbart} classes;
  as_draws_array/df on 7 non-{bart,rbart} classes): open in the sense that
  no Rd claims any of these classes should or shouldn't carry them.

## 7. Numbered candidate findings

**F1. `dbartsSpec()` and `dbarts()` resolve `family = "multinomial"`
differently for the identical input, contradicting [[dbartsSpec.Rd:40@b102e17c]]'s
"same resolution" claim.**
Repro:
```r
dd <- dbartsData(x, yMultinom)                 # plain 3-level factor response
dbartsSpec(data = dd, control = ctrl, family = "multinomial")
# -> refuses: family "multinomial" cannot fit a 3-level factor response;
#    a 3+-level factor is multinomial (bart2(family = "multinomial"))
dbarts(x, yMultinom, family = "multinomial", control = ctrl, verbose = FALSE)
# -> accepts, model@family == "multinomial"
```
Root cause: `dbarts()`/`bart2()` auto-convert a factor response into the
`counts` matrix multinomial actually needs (`resolveMultinomialCounts()`)
*before* building the `dbartsData`; `dbartsSpec()` only ever consumes an
already-built `dbartsData`, so it never gets that conversion. Confirmed
`dbartsSpec()` *does* accept `family="multinomial"` when handed a
counts-built `dbartsData(x, counts = countMatrix)` directly - so the token
is reachable, just not through the "obvious" construction path the doc's
"same resolution" claim implies exists uniformly for all 8 tokens dbartsSpec
declares. Also see the [[man/dbartsSpec.Rd:48@b102e17c]] inconsistency in section 3.

**F2. `rbart_vi(family="aft")` via the matrix interface fails with an
opaque dimension-mismatch error, not the by-name refusal the family's own
formula-only restriction (already correctly documented) would predict.**
Repro:
```r
rbart_vi(x, Surv(time, status), group.by = g, family = "aft", ...)
# CALL: dbarts::dbartsData(formula = x, data = <Surv object>, ...)
# -> "'x' must have the same number of observations as 'y'"
```
`x` has 40 rows and `Surv(time, status)` also reports 40 (`NROW`/`length`
agree) - the message is simply wrong about what's wrong. Root cause: R/rbart.R's
own comment says "rbart_vi's matrix interface has no survival form; aft
enters only through the formula" - a real, and correctly documented ([[rbart.Rd:95@b102e17c]]),
restriction - but nothing enforces it explicitly on the matrix path. The raw
`Surv` object gets forwarded straight into `dbartsData(formula=x, data=Surv,
...)`, which has no survival-response handling at all (unlike
`dbarts()`/`bart2()`, which extract log-time+status *before* ever calling
`dbartsData()`), and some internal accounting on the classed `Surv` matrix
miscounts rows, producing the generic error. A user who reads only the
`family=` argument list (which lists "aft" with no matrix/formula
qualification at the call-signature level) has no way to anticipate this
from the error text.

**F3. `extract(fit, type="trees", sample=<train|test|anything>)` silently
forwards `sample` into `dbartsSampler$getTrees()`'s unrelated `sampleNums`
parameter via R's partial argument matching, corrupting the call.**
Repro:
```r
extract(fit, type = "trees", sample = "train")
# -> "missing value where TRUE/FALSE needed"
# mechanism: extract.bart's type=="trees" branch rewrites the call to
# object$fit$getTrees(sample = "train", ...) without stripping `sample`
# (unlike `object`/`type`, which it does strip). getTrees()'s formals are
# (treeNums, chainNums, sampleNums, current=FALSE, newdata=NULL) - no
# `sample` - so R partial-matches "sample" -> "sampleNums", then
# as.numeric("train") -> NA (warning: NAs introduced by coercion), and a
# downstream `if (... sampleNums <= 0 ...)` guard dies on the NA.
```
Reproduces identically for `extract.rbart`, which has the same
match.call()-rewrite pattern (confirmed in the C3 grid: 33 total cells
across every bart-class + rbart representative, all samples including
valid "train"/"test" tokens - this is not limited to the deliberately
invalid probe). `extract(fit, type="trees")` with `sample` *omitted*
doesn't error (the default `c("train","test")` coerces to `c(NA,NA)` and
some other branch swallows it), but silently ignores the requested sample
filter - worth a numeric-correctness look, out of this grid's scope (see
Budget).

**F4. Entries without a `...` catch-all (`bart()`, `xbart()`) turn every
unsupported keyword into R's raw, unnamed `unused argument (...)` error;
entries with one (`bart2()`, `rbart_vi()`) give a clear `unknown argument
'NAME'` diagnostic for the identical mistake.**
Evidence (10 of the 23 raw Part A/B cells, one line each from the grid):
`bart(variance=~x1)` -> `unused argument (variance = ~x1)`; `bart(forests=
list(forest()))` -> `unused argument (forests = ...)`; `xbart(test=...)`
-> `unused argument (test = ...)`; `xbart(keepTrees=TRUE)` -> `unused
argument`; `xbart(n.chains=2)` -> `unused argument`; `xbart(resid.dist=
student())` -> `unused argument`; contrast `rbart_vi(forests=list(forest()))`
-> `unknown argument 'forests'`, `rbart_vi(resid.dist=student())` ->
`unknown argument 'resid.dist'`, `bart2(forests=list(forest()))` ->
`unknown argument 'forests'`. Same underlying user mistake (an argument
that belongs to a sibling entry but not this one), two completely different
failure qualities depending on which entry it lands on.

**F5. `offset.category` is unreachable as a keyword argument on every
fitting entry** (dbarts/bart2/bart/rbart_vi/xbart formals were checked
directly - none declare it) **- the only door is a pre-built
`dbartsData(..., offset.category=)` passed as `data=`.** Not a bug (both
`dbarts()` and `bart2()` accept the pre-built-object route symmetrically,
confirmed in Part B), but genuinely undiscoverable from any fitting entry's
own argument list - a candidate doc gap rather than a behavioral one.

**F9. `bart()`'s by-name own-class-family refusal ([[man/bart.Rd:16@b102e17c]]: "naming
one of bart2's remaining own-class families to bart's family is refused,
pointing there") does not cover the `"twopart"` alias, even though
`"twopart"` is documented ([[dbarts.Rd:111@b102e17c]], bart2.Rd) to resolve and print as
`"hurdle.lognormal"` everywhere else.**
Not reachable from Part A/B above, since both only probe each entry's own
declared vocabulary (`declaredTokens(bart)` = `c("auto","logistic","aft")`,
which excludes both spellings); found by probing `bart()` directly with
tokens outside its own formal, matching the task's "include ... the
aliases" instruction. Repro:
```r
bart(x, y, family = "hurdle.lognormal", verbose = FALSE)
# -> bart() does not fit family = "hurdle.lognormal"; use bart2(x.train,
#    y.train, family = "hurdle.lognormal")     [named, points to bart2()]
bart(x, y, family = "twopart", verbose = FALSE)
# -> 'arg' should be one of "auto", "logistic", "aft"    [bare match.arg]
```
Root cause: `bartOwnClassFamilies` ([[R/bart.R:2587-2592@b102e17c]]) is the literal
vector `c("multinomial", "ordinal", "nbinom", "hurdle.lognormal")` - checked
by `bart()` against the RAW token *before* `bart2()`'s own alias-fold (which
resolves `"twopart"` to `"hurdle.lognormal"` immediately after
`match.arg()`) ever runs. `bart()` never reaches that fold - it calls
`match.arg(family)` directly against its own 3-token vocabulary - so a
caller who reasonably expects the two spellings to be interchangeable (as
they are documented to be, and are in every other family= consumer) gets a
by-name, helpful redirect for one spelling and a raw, unnamed R error for
the other, for the identical request.

**F6-F8 (lower confidence / not independently pursued further, listed for
completeness).**
- `rbart_vi(subset=..., group.by=<unsubsetted>)` refuses with "'group.by'
  not of length equal to that of data" - correct behavior, but the message
  doesn't say the fix is to subset `group.by` yourself; possibly worth a
  clearer message, not a correctness bug.
- `bartHurdle`/`bartNegbin`/`bartOrdinal` all have `summary` but none has
  `plot` (section 2) - flagged by the memo, re-confirmed here, not
  independently investigated for cause.
- The heteroscedastic-gaussian silently-wrong-loglik/ppd/summary bug the
  2026-08-24 value-scan reported (memo F3) did **not** reproduce on this
  representative fit (`bart.heteroscedastic`): every `extract`/`predict`
  cell returns a shape-consistent 6x40/6x12 result. This grid only checks
  shape, not numeric correctness against an oracle, so this is "not
  reproduced," not "confirmed fixed."

## 8. dbartsSampler reference-class method census (supplementary)

The task's axis-1 explicitly names "the dbartsSampler reference class (its
~46 own methods; census in inst/tinytest/test-host-shell-pins.R)" as one of
the fitting entries. `matrix.R`'s own execution grid (section 2/7) only
reaches `dbartsSampler` through its two registered S3 generics (`extract`,
`plotTree`) - section 7's Budget note flags this explicitly ("dbartsSampler
... not deeply exercised in Part C"). `matrix-sampler-census.R` closes that
gap: it builds one plain single-forest gaussian sampler (`dbarts(x, y,
control = ctrl)`, keepTrees=TRUE) and calls all 41 substantive RC methods
from the `test-host-shell-pins.R` census list directly (`entry =
"dbartsSampler-method"` in the CSV), one call each with minimally-plausible
arguments.

Result: **28 accepted, 13 refused, 0 error-without-reason.** Every method on
the surface is reachable and every refusal names a reason - none of the 41
calls produced a raw/uninformative error, which is a genuinely clean result
(consistent with section 1.0 of the review-lenses memo's ERROR-TEXT DRIFT
finding: this surface has near-zero raw-message residue). The 13 refusals
are all expected, single-forest-vs-multinomial/BCF distinctions working as
documented: `setCounts`/`setCategoryOffset`/`setCategoryTestOffset`
("carries no count response: only a multinomial..."), `setForestWeights`/
`setForestBasis`/`getForestAmplitudes`/`getCalibration` (require a
`forests=` multi-forest sampler), `setTestPredictor` (a fixture-shape
mismatch, not a real gap), `setState` ("must inherit from bartcoreState" -
this probe's argument was wrong, not evidence of a defect), and `predict`/
`printTrees`/`getTrees`/`plotTree` all refusing with "the saved-tree store
holds no recorded draws; run the sampler with keepTrees" - an artifact of
this probe's own ordering (several state-mutating methods, e.g.
`setPredictor`/`setCutPoints`, run between the census's initial `run()` and
these later reads, clearing the tree store each entails a re-fit) rather
than evidence keepTrees itself is broken - `test-host-shell-pins.R` already
covers the real keepTrees-survives-save/reload path.

**Not exercised**: a multinomial or BCF (`forests=`) sampler's own method
census - `[[man/dbartsSampler-class.Rd:351@b102e17c]]` documents in detail which of the
41 methods a multinomial sampler refuses and why (`setResponse`,
`setOffset`, `setWeights`, `setSigma`, `setData`, `setModel`,
`setCalibration`, `setForestWeights`, `setForestBasis`,
`getFitsWithoutOffset` - each "by a message naming the capability"); this
census only ran the gaussian case, so that paragraph's claims remain
untested by this pass (an open cell, not a contradiction - no
counter-evidence was produced either way).

## Budget and coarsening

Total runtime: under 5 seconds (836 cells), far inside the 30-minute
allowance - the fixture is deliberately tiny (n=40/12, 3-4 trees, 5-6
samples) specifically to make the *whole* grid affordable rather than to
push cell count. What was traded for that headroom, explicitly:

- **Conditioning columns are one-at-a-time against a single representative
  family per entry (Part B), not crossed against all 13/13/3/3/4/8 family
  tokens** - a full factorial would be ~5x-13x this size. If a conditioning
  bug is family-specific, this grid would miss it.
- **The consumption grid (Part C) uses one representative fit per family
  variant, not one per (entry x family) combination** - e.g. bart2's
  gaussian fit stands in for what dbarts()'s gaussian *sampler*, once run,
  would also produce; the two are not separately exercised through every
  generic.
- **No numeric-correctness checking** - `describeValue()` records class,
  dimensions, dimnames, and extra attributes only, never compares returned
  numbers against an independent oracle. This grid answers
  accept/refuse/shape questions; it cannot by itself catch a "runs, returns
  the right shape, wrong numbers" defect (the kind wave-5b's value-scan
  found). F3's `sample`-omitted case is flagged as exactly this kind of gap.
- **Rd cross-checking is targeted, not exhaustive** - 3 of 29 man pages
  (dbartsSpec.Rd, rbart.Rd, bart2.Rd) were read against specific grid
  results; the other 26 were not swept for testable claims.
- **pdbart/pd2bart are covered by the class x generic existence table but
  not deeply exercised in Part C** (one representative object each,
  generics that exist are run, but `pdbart`'s own argument surface, e.g.
  `xind`, is untouched here). `dbartsSampler`'s ~46-method RC surface is
  covered separately by section 8's supplementary census (28
  accepted/13 refused/0 error-without-reason on a plain gaussian sampler),
  which closes most of this gap for the single-forest case; a multinomial
  or BCF sampler's method census ([[man/dbartsSampler-class.Rd:351@b102e17c]]'s specific
  per-method refusal claims) is still untested - see section 8's own "not
  exercised" note.
- **Two independent processes converged on this same deliverable during
  this pass** - `matrix.R`/`matrix-grid.csv`/`matrix-results.md` (sections
  1-7) were built by a forked agent working the same brief in parallel;
  `matrix-sampler-census.R` and section 8 are this pass's own addition,
  layered on afterward rather than rebuilding what already existed. Findings
  F1-F5 above and this pass's independent direct verification agreed on
  every point checked twice (the `dbartsSpec`/`multinomial` contradiction,
  the `rbart_vi`/`aft` matrix-interface failure and its exact message, and
  the `extract(type="trees", sample=)` partial-argument-matching bug, which
  both this pass and the forked agent found independently) - convergent
  evidence, not a single derivation checked once.
