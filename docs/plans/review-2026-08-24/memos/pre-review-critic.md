One of four adversarial pre-review reports read ahead of the manual review; VD's rulings on its findings are recorded in docs/plans/pre-review-cleanup.md.

# Pre-review completeness critique - dbarts bartcore @ 3080a9c5

Read-only. Verified against a fresh install of this tip in a private library;
every "verified" claim below was run, not inferred.
NEWS db count confirmed at **346** rows (218 in 1.0-0: 20 UPGRADING / 68 NEW
FEATURES / 11 C API / 119 BUG FIXES) - matches expectation.

The claim "ready for a human to read end to end and find nothing that should
have been fixed first" is **refuted**. Twelve PRE-REVIEW items follow, the first
five of them user-visible defects rather than polish, all cheap.

---

## PRE-REVIEW

### 1. `n.thin` greater than `printEvery` hard-fails every front end, naming an argument the user never set. ~3 lines.

`R/bart.R:811`, `R/bart.R:2849`, `R/rbart.R:100` each do

```r
control@printEvery <- control@printEvery %/% control@n.thin
```

With the default `printEvery = 100L`, any `n.thin >= 101` drives it to 0, and
the bridge refuses 0 at `src/R_interface_bartcore.cpp:432`
(`rc_getInt(..., RC_VALUE | RC_GEQ, 1, ...)`). Verified live, all three:

```
bart(x, y, ndpost=400, nskip=200, keepevery=200L)   -> "print every must be greater than or equal to 1"
bart2(x, y, n.samples=400, n.burn=200, n.thin=200L) -> same
rbart_vi(..., n.samples=400, n.burn=200, n.thin=200L) -> same
```

Nothing fractional, nothing undocumented, nothing the caller typed - heavy
thinning is an ordinary request and the error names `print every`. The R5 path
does *not* divide (`dbarts()` + `$run()` with `n.thin = 200L` runs fine,
verified), so `printEvery` also means different things on the two paths.
Cure: `max(1L, ... %/% ...)` at the three sites.

Same knot, second contradiction: `man/dbartsControl.Rd:56` says printEvery
"Must be a positive integer"; `R/A_class.R:356-357` accepts 0 ("must be a
non-negative integer"); the bridge refuses 0. Verified: `dbartsControl(printEvery
= 0L)` passes `validObject`, and every `dbarts()` built from it is refused.
Three contracts, one argument. Cure: `< 1L` + "positive" in the validity message.

### 2. The new dbartsSampler `extract` catch-all is not a catch-all, and can print an empty name. 2 lines.

`R/generics.R:2864-2872`:

```r
refuseSamplerExtractArgs <- function(dots) {
  if (length(names(dots)) > 0L) {
    stop("'", names(dots)[1L], "' is not used by extract on a dbartsSampler: ...")
```

Verified: `extract(s, "predictors", 1)` is **silently accepted** (all-unnamed
dots give `names(dots) == NULL`), and `extract(s, "predictors", 1, foo = 2)`
raises

```
'' is not used by extract on a dbartsSampler: this method returns the sampler's coded predictor matrix
```

Both contradict the comment three lines above it ("a caller-supplied name of any
kind is refused") and `inst/NEWS.Rd:139` ("gain their own catch-all refusals").
Cure: `if (length(dots) > 0L)`, and pick the first *named* element for the
quoted name.

Same root, second instance, on the shipped fit surface:
`refuseUnusedGenericArgs` (`R/generics.R:2028-2043`) matches on
`names(dots)` only, so every *positional* extra is still discarded. Verified:
`residuals(ordinalFit, "ev")` returns silently, while
`residuals(ordinalFit, type = "ev")` raises the designed refusal. `type` is the
sibling spelling `residuals.bart` takes positionally, so this is the exact
migration a user makes. This is the residue at `docs/plans/surface-refusals.md:768-770`
("additive, so post-1.0") - but the named half already shipped today, and the
unnamed half is ~6 lines in one helper, not a new surface.

### 3. NEWS claims a `plotTree` catch-all that does not exist. 1 line.

`inst/NEWS.Rd:139-143` - "`extract` on a `dbartsSampler` and `plotTree` on a
`dbartsSampler` gain their own catch-all refusals". `refusePlotTreeArgs`
(`R/generics.R:2986-2999`) refuses exactly two names, `sample` and `chain`.
Verified: `plotTree(s, foo = 1)` falls through to R's own
`argument "treeNum" is missing, with no default`. The surface-refusals landing
note (`docs/plans/surface-refusals.md:816-818`) already narrowed this claim once
and still overshot. Cure: drop "catch-all" for the plotTree half.

### 4. `dbarts_sampler_setVerbose(s, 1, 0)` divides by zero in the engine. 2 lines.

`src/C_interface.cpp:969-972` forwards `printEvery` unvalidated;
`src/bartcore/chain.hpp:1393` evaluates `(iteration / numThin + 1) % options_.printEvery == 0`
under `options_.verbose`. `printEvery == 0` is integer modulo by zero: SIGFPE on
x86-64, a silently-true condition on arm64 (prints every iteration). "0 = never
print" is the obvious consumer reading, and `dbarts.h:1082-1083` is the one
entry on that page carrying **no** doc comment at all while its neighbours carry
paragraphs. Recorded at `docs/plans/capi-shape.md:613` as "still does not guard
0", without noting that it is a crash. (`setNumThin(0)` is benign - checked:
`totalIterations` goes to 0 and the loop never runs.) Cure: clamp in
`C_interface.cpp`, or `options_.printEvery != 0 &&` at the engine condition, plus
one header sentence.

### 5. The `fitted()` slot-3 reorder gives migrating users a message that names neither argument. ~4 lines.

`fitted(object, type, sample, ...)` is the **released** signature
(`main:R/generics.R:117-120`), so `fitted(fit, "ev", "test")` is live 0.9-x code
in the wild. On this tip, verified:

```
'ci.level' must be a single number in (0, 1)
```

`sample` is never mentioned, and `ci.level` is an argument the caller has never
heard of. This is the single most-used accessor in the package and the reorder
is the flagship UPGRADING entry (`inst/NEWS.Rd:152-158`). Cure: in the ci.level
validator, special-case a character value in `c("train", "test")` and point at
`sample`'s new fourth position.

### 6. `predict()`'s "uniform" order still has one silent positional collision. 2 lines (doc) or ~12 (formal).

Live formals, slot 4: `offset` on bart / rbart / bartMultinomial / bartNegbin,
`combineChains` on bartOrdinal / bartHurdle. Verified:
`predict(ordinalFit, x, "ev", FALSE)` binds `combineChains = FALSE` and returns a
20x60x3 array; `predict(ordinalFit, x, "ev", offset = z)` correctly raises
`'offset' is not used by predict on a bartOrdinal fit: ...`. This is the only
remaining positional collision on the six predicts, and it is the same defect
class `fitted()` was reordered to remove today.

Recorded as residue at `docs/plans/predict-surface.md:186`. What makes it a
pre-review item is `inst/NEWS.Rd:76-79`, which opens "Every `predict` method's
argument order is now uniform" - a reader takes that as licence to pass
positionally. Cheapest honest cure is two clauses (NEWS + the bart2 Value
paragraph) saying the order is uniform in *relative* order only and that absolute
slots shift where a class lacks a channel; the thorough cure is an `offset = NULL`
formal in slot 4 on both, refusing non-NULL.

### 7. A misleading doc anchor introduced today. 1 token.

`docs/design/multinomial-mutation-arc.md:754-755` - "`R/generics.R:1688` and
`:1869` now say so directly" (about `$bc`'s removal). `R/generics.R:1686-1690` is
`print.bartOrdinal`'s `cat("n.trees: ", ...)` block. The real comment is
`R/generics.R:1590-1592`. `bfc4571d` re-targeted this today, 1501 -> 1688, and
missed; the paired `:1869` cite *is* correct, which is what hides it. This is
one of the four new freshness advisories (69 -> 73) and the only one a reviewer
would trip on - the other three are benign nearest-identifier drift on
`feature-matrix.md:358` (RIB:2791, RIB:2928, still landing on the right
predicates) and `multinomial-mutation-arc.md:654`
(`multinomial-equivalence.R:189-196`, real anchor now ~472-475, also today's).

### 8. A landing record written today whose own anchors are all wrong. 4 tokens.

`docs/plans/bcf-cross-host.md:38`, `:126`, `:536-537` cite the
*pre-implementation* file and were never re-targeted by `de67cf07`. Concretely:
`bcf-equivalence.R:128-140` is now the `statChannels`/`snapshotChannels` vectors
(`recordChannels` moved to `:412`); the "Welch z at `bcf-equivalence.R:511-518`
(and `multinomial-equivalence.R:603-612`)" now lands in scenario bodies - the z
is at `:826-831` / `:911-916`. Below the freshness checker's radar (it walks
`docs/design`, and it exits 0 on WARN regardless: `.github/workflows/lint.yaml:89`).

### 9. Today's new per-push gate can pass while comparing nothing. ~4 lines.

`benchmarks/R/bcf-equivalence.R:764-767` (and the same shape in
`multinomial-equivalence.R`):

```r
if (is.null(b)) {
  cat(sprintf("%-14s skipped (not produced this run)\n", name))
  next
}
```

`anyFailure` is untouched. Drop or rename every scenario and both scripts print
the full `OK: every gated channel within the cross-host tier 1 bound` and exit 0.
`equivalence.R` has `--strict-coverage` for exactly this; neither sibling does -
and `exact-gates.yaml` now runs both on every push, so the coverage hole is
newly load-bearing. Cure: set `anyFailure` in the skip arm, or assert a produced
count.

### 10. The only test pinning today's `plot.bartMultinomial` reorder is vacuous. 1 line.

`inst/tinytest/test-plot-generics.R:372`:
`expect_silent(plot(fit.multinomial, c(0.1, 0.9)))`. Verified: `plot(fm, cols =
c(0.1, 0.9))` is **also** silent (numeric colours coerce), so the assertion
passes under both the old and the new argument order and cannot detect the change
it exists to pin. Cure: `expect_error(plot(fm, c("blue", "black")))` - verified to
raise `'probs' outside [0,1]` under the new order and to be silent under the old.

### 11. Error-style non-conformance inside today's own hunks. ~6 lines.

Governed by `docs/design/error-style.md`. Not the whole corpus - only messages a
reviewer meets inside a diff hunk landed today:

- `R/generics.R:2883` `stop("type must be 'predictors'")` sits two lines below the
  new conforming refusal in the same function: bare `type`, quoting inverted onto
  the value. Should be `"'type' must be one of 'predictors'"`.
- `R/utility.R:179-180` (new today) guards `mc[[2L]]` with `is.symbol()`, while
  `:164` and `:184` in the *same 20-line function* interpolate it unguarded. The
  asymmetry is visible in one screen. (`docs/plans/surface-refusals.md:771-774`
  defers these to error-style slice L - defensible for the corpus, not for two
  lines bracketing today's insertion.)
- `R/utility.R:180` `" must be a whole number, not ", deparse(x)[1L]` uses R10's
  `, not <type>` shape to carry a *value*; R12's shape is `; got '<value>'`.
- `R/generics.R:2092-2095` `bartUnusedArgs$contribution` is a verbatim copy of the
  `forest` reason and never mentions contribution; the same name now carries three
  different explanations (`:518`, `:2092`, `:2114`).
- `R/generics.R:2112` leaves `bases` unquoted while its list-mates quote
  `'sample'`, `'group.by'`.
- `R/generics.R:952` uses `;` for an explanation where the new siblings at
  `:1234/:1562/:1843` use `:`, so one method now emits both separators for
  adjacent argument names.

### 12. `summary`'s empty-table message names the wrong parameter set. ~3 lines.

`R/diagnostics.R:415`: `"No scalar parameters (sigma, k, tau) to summarize."` -
printed for `bartOrdinal` (default `vars` is `c("cutpoints", "sigma", "k",
"tau")`) and `bartNegbin` (`c("dispersion", ...)`) too. Low reachability, but it
is in the printed output a reviewer eyeballs first. Cure: interpolate the
requested `vars`.

---

## RESIDUE-OK (correctly deferred - priced and dismissed)

- `setResponse`/`setOffset`/`setWeights`' flat-unreachable 0 arms
  (`capi-shape.md:587-591`): stated contract, pinned by the BCF legs. Correct.
- `setWeights` finiteness at the flat entrance (`:609-610`): C API is
  fast-over-safe by house rule and R validates. Correct.
- `setActiveRows` scanning its mask twice (`:614`): perf only.
- The R5 twins of `getTrees`/`printTrees` carrying no forest selector (`:611-612`):
  real surface gap, but well over 30 lines and additive. Correct to defer.
- `extract(type = "trees")` raising R's own "unused argument"
  (`surface-refusals.md:775-778`): disagrees in voice, not outcome.
- `fitted(type = "ppd")` consuming RNG (`:779-781`): documented, not a defect.
- `as_draws_*.bartMultinomial`'s inert `vars` (`:782-784`): VD's call, and it *is*
  documented - `man/summary.bart.Rd:97-105` says its `vars` "only ever means
  `meanProb`". Verified inert live. Correct.
- `plot.pdbart`'s fourth shape (`:791-793`), `makeind(all=)`, `run(n.threads)`:
  correct.
- `plot`'s `cols` default splitting `NULL` (multinomial/ordinal) from
  `c("blue","black")` (the other four): **documented** at `man/bart2.Rd:353-355`.
  Not a finding.
- `plotTree.bart`'s `chainNum` having no default while `plotTree.rbart` defaults to
  `1L`: **documented** at `man/plotTree.Rd:53-57`. Not a finding.
- `fitted`'s type vocabulary lacking `"ppd"` on multinomial/ordinal where predict
  has it: a mean over category codes is meaningless; implicitly documented.
- D8 residues (`composition-refusals.md:424-435`), xbart fold-view scope: the
  xbart carve-out is now stated in `man/xbart.Rd`. Correct.

## NOTE

- `docs-site/reference/bart.md:84-89` still shows the pre-today `fitted` order.
  Untracked and `.Rbuildignore`d (line 23), so it ships nowhere - but the whole
  `docs-site/` tree predates today and must be rebuilt before any publish.
- 8 of 11 hash-pinned files in `benchmarks/baselines/` name commits reachable from
  no ref - including all three current ones (`6e3b9fb8`, `4d9a3337`, `ab1dc52`;
  verified `git branch -a --contains` empty for each). They are loose objects in
  this clone only, so "compare against baseline-`<hash>`" is unreproducible from a
  fresh clone. D7's pending RC-tip re-record fixes one of the three.
- `benchmarks/R/composition-matrix.R` runs in no workflow (zero hits in
  `.github/`), while `docs/plans/INDEX.md:217` records a 0-disagreement state CI
  will never re-check. Its oracle is genuinely external (tokens parsed out of
  `feature-matrix.md`), so today's two closures are real; but 116 of 299 cells are
  never probed and `CONSTRUCTS` means only "threw nothing".
- A tier-1 FAIL adjudicated by tier 2 exits 0
  (`bcf-equivalence.R:896-899`), and nothing in `exact-gates.yaml` greps for
  `decoupled:` although `docs/plans/bcf-cross-host.md:252-253` calls such a line
  "a finding to report". Enforced only by a human opening a collapsed group.
- `inst/tinytest/test-capi.R:73-89`'s six retired-token `expect_false`s are
  dominated by the live-token `expect_identical` at `:92` - diagnostic naming, not
  coverage. Not stale (all six literals were genuinely live), just not load-bearing.
- The R5 `getTrees` range messages (`R/dbarts.R:2059/:2065/:2068`) name
  `'chainNums'`/`'sampleNums'` when reached through `plotTree`, whose own
  arguments are `chainNum`/`sampleNum`. Verified: `plotTree(s, 1L, 99)` ->
  `'chainNums' must be in [1, 1]`.

## Everything else - one paragraph

Verified clean, and worth saying so because these were the likeliest places to
find rot. Every `\value`/argument claim edited today that I could execute is
**true**: `print.bart`'s synopsis and invisible return; `xbart`'s six-axis rule
(vector collapse at all-length-1, `rep x n.trees x k` at 2x2x2, `drop = FALSE`
giving 2x1x1x1x1, integer dimnames on `n.trees` and 2-sig-fig on `k`, `rep`
unnamed); every branch of `dbartsSampler$predict`'s new Value paragraph
(n.test x n.samples with the chain axis dropped at one chain, the live-tree
vector/matrix, the heteroscedastic `list(mean, variance)` collapsing to a bare
array without saved trees, multinomial's never-dropped K axis at 5x3x20);
`rbart`'s `$fit` list-of-one on the in-core path with `$n.chains` present either
way; the `getSigmas`/`getSumsOfSquaredResiduals` `result` refusals; the
`class`+`ci.level`, `forest`-outside-forest-arm, `residuals(ci.level=)`,
`summary.bartMultinomial(vars=)` and `group.by`-by-name refusals, each with the
designed wording. `inst/include/dbarts/dbarts.h` is internally consistent with
its new signatures (capability paragraphs on all seven converted setters, the
forest-position rule matching `getTrees`/`printTrees`, `basisRowMajor` named as
the column-major exception, no `USE_FC_LEN_T`); `inst/tinytest/capi/consumer.c`
is fully converted; `vignettes/` and `README.md` carry no stale surface;
`NAMESPACE`/`DESCRIPTION` are correct, with `as_draws_*` registered dynamically
in `R/hooks.R:18-24` because posterior is Suggests-only. The NEWS 1.0-0 section
reads coherently top to bottom - I found no entry reversed by a later one.
Findings 3 and 6 are the only two places where the section overclaims. One
near-miss for the record: `inst/NEWS.Rd:2067-2070` ("a `keepSampler` fit ... now
always carries `$n.chains`") was written at `08db8da1`, and the kept-sampler path
that makes it true only landed today at `52c10e02` - the entry described the code
before the code existed. It is accurate as of this tip; noted only because it is
the one place NEWS ran ahead of `R/`.
