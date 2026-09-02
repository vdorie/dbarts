# Error-message style: design

Status: ADOPTED for new messages 2026-08-17 (release-candidate-review
wave 0a; VD resolved Fork 5 as "follow best practices from highly
regarded R packages" and refined it the same day: tidyverse practices
carry no authority weight and are adopted only on their own scrutinized
merits. This revision applies both - the draft codified the in-repo
majority, the measured practice of base, stats, Matrix, survival, lme4,
and mgcv takes precedence wherever it clearly speaks, and each
guide-originated adoption states its own merit case). The wave-2 L
sweep of the existing corpus runs after VD sees the delta summary.
Governs every `stop()` in `R/` and every
`Rf_error`/`ext_throwError` in `src/`. Every NEW message from this point on
follows this rule; slice L (repo-wide, ~120 strings) rewrites the existing
corpus against it in one sweep, deliberately last among the message-touching
slices so nothing is reworded twice. Text only - no function is renamed, no
predicate changes (predicate reconciliation, e.g. `any(w < 0)` vs
`!(w >= 0.0)`, is slice K6, scheduled after this lands).

Evidence: 610 `stop()` in `R/`, 365 `Rf_error()` in `src/R_interface_bartcore.cpp`
+ `src/C_interface.cpp` + `src/R_interface.cpp`, 3 `ext_throwError` in
`src/bartcore/chain.hpp:1886,2050` and `src/bartcore/sampler.hpp:682`. Each
rule states the majority it codifies, or
says plainly that it invents (no majority existed). Frequencies in the
appendix. Every rule below was additionally checked against: a direct
source-level survey of `stop()` in base, stats, Matrix, survival, lme4,
and mgcv (the observable practice of highly-regarded base-style packages -
the tradition dbarts actually belongs to, since it ships no
rlang/cli/tidyverse dependency and never will); Writing R Extensions
(WRE); and, scrutinized rather than deferred to, the tidyverse style
guide's error chapter and rlang's `abort()` documentation. Full method and
citations are in the **External evidence appendix** at the end of this
document. Evidence precedence, per the project lead's direction
(2026-08-17): the measured practice of the highly-regarded base-style
packages wins over the in-repo majority wherever it clearly speaks; the
tidyverse/rlang guides carry NO authority weight of their own - a practice
that originates there is adopted only when it stands on its own
demonstrated merits under scrutiny (each such adoption below states that
merit case explicitly); where all of that is silent or conflicting, the
in-repo majority stands. Each rule below is tagged **CONFIRMED**,
**OVERRIDDEN**, or **SILENT-KEPT**.

## R1. Argument/formal-name quoting: single quotes — CONFIRMED

Wrap an argument, formal, or slot name in `'...'`. Never backtick, never bare,
never `sQuote`/`dQuote`/`gettextf` (0 uses of any of these on either side -
do not introduce them now).

- Conformance: `"'weights' must have the same length as 'y'"` (`R/data.R:772`).
- Violation, since reworded: `"chainNum must be a single chain index in
  [1, ...]"` was bare; `R/generics.R:3030` now reads `stop("'chainNum'
  must be a single chain index in [1, ", n.chains, "]")`.

A value drawn from a closed set of choices (family name, class name used as an
echoed value) is quoted the same way; a descriptive category noun used as the
sentence's own subject stays bare - `"probit models do not support weights"`
(`src/R_interface_bartcore.cpp:2784`), not `"'probit' models"`.

**External evidence.** The two external sources disagree, and the disagreement
is the crux of this rule. Published best practice (tidyverse style guide):
"Surround the names of arguments in backticks, e.g. `` `x` ``. Use 'column'
to disambiguate columns and arguments: `` Column `x` ``." That is an explicit,
unambiguous rule for backticks. Observable practice of highly-regarded
base-style packages says the opposite, and says it almost unanimously: across
~3,400 sampled `stop()` calls in base, stats, Matrix, survival, lme4, and mgcv,
backtick-quoted names appear in exactly 26 messages (0 in base, 0 in stats, 0
in Matrix, 1 in survival, 0 in lme4, 25 in mgcv - and even within mgcv,
single-quoted names still outnumber backtick-quoted ones 67 to 25). Base R
itself - the most authoritative "highly regarded package" of all, and the one
whose error-message idiom dbarts already imitates throughout this document -
uses single quotes in 703 of its own quote-pairs and backticks in 0.
Representative: `stop("'input' must be a character vector or 'NULL'")`
(`base::system`); `stop("'shape' must be one of \"g\", \"t\", \"s\"")`
(`Matrix::.sparseDiagonal`).

Resolution: the backtick convention in the tidyverse guide is not a
free-standing typographic preference - it exists because rlang/cli render
backtick-quoted spans specially (monospacing, sometimes color) in a
multi-line bulleted display. Lifted out of that rendering context into a
plain `stop()` string, a backtick is just a character with a second, older,
competing meaning in R itself: it denotes a non-syntactic name *in code*
(`` `my var` <- 1 ``), which is exactly the kind of confusion a reader
skimming a one-line error should not have to resolve. Per the task framing
("adopt textual/structural conventions only insofar as they translate to
plain stop()/Rf_error strings"), the backtick rule does not survive
translation; single-quote is the shape that does, and it is also what every
sampled highly-regarded base-style package - most authoritatively base R
itself - actually writes. **R1 is CONFIRMED unchanged.**

One more corroborating data point: base R's own quoting helpers, `sQuote()`/
`dQuote()`, are used 119 times across the six sampled packages (70 in base
alone) but were already excluded by the in-repo majority (0 uses in dbarts).
That exclusion is now doubly justified: `sQuote`/`dQuote` render *fancy*
(locale-dependent, non-ASCII) quote glyphs by default
(`getOption("useFancyQuotes")`), which would violate R4. Literal ASCII
`'...'` is the correct base-R-flavored quoting *without* the fancy-quote
liability - keep it exactly as R1 already specifies.

## R2. Sentence case: lowercase-initial — CONFIRMED

Begin the message lowercase, unless the first token is an established acronym
or proper noun always capitalized on its own (BCF, DART, R, PG) - keep its
natural casing rather than force `"bcf does not support"`.

- Conformance: 546/546 `stop()` bodies are lowercase-initial
  (`grep 'stop("[A-Z]' R/*.R` -> 0).
- Violation: none left on the C side. The prior outlier ("Student-t
  residuals ...") went with its message; `:7753` (DART) is the acronym
  exception and stays as-is; the two BCF-initial messages that sat beside it
  were reworded when the amplitude family took amplitude-rooted names,
  leaving DART the only one.

**External evidence.** Published best practice (tidyverse): "Errors should be
written in sentence case ... make sure to capitalise the first word (unless
it's an argument or column name)." That directly contradicts R2. Observable
practice again disagrees with the published rule and agrees with the draft:
lowercase-initial is the overwhelming majority in every sampled package -
base 559/584 (95.7%), stats 880/916 (96.0%), Matrix 177/193 (91.7%), lme4
186/211 (88.2%), mgcv 391/509 (76.8%), survival 595/921 (64.6%, the weakest
but still a clear majority). Base R and stats, the two most authoritative
samples, are both above 95%. **R2 is CONFIRMED unchanged**, on the same
tradition-membership reasoning given in full under R3 below (R2 and R3 are
one design decision, not two).

## R3. Terminal period: never — CONFIRMED

No message ends in `.`: 348/348 on the C side, 549/549 on the R side. The
former exception - one refusal restated three times, each copy ending in a
period - is fixed: `R/spec.R:54-59`, `:73-77`, `:82-86` (probit/ordinal/nbinom
weight refusals).

- Conformance: any `Rf_error` call (all 348 are period-free).

**External evidence.** This is the rule where the two traditions most
plainly conflict, and where the task asks for a one-paragraph resolution.
Published best practice (tidyverse): errors "should end in a full stop,"
matching sentence case (R2) and a bulleted, hint-terminated, multi-line
structure - full-sentence punctuation is part of one coherent rendering
package (rlang's structured condition + cli's colored bullets). Observable
practice of highly-regarded base-style packages is the mirror image, and
overwhelmingly so: no terminal period in 581/585 base messages (99.3%),
903/915 stats (98.7%), 191/193 Matrix (99.0%), 906/913 survival (99.2%),
206/211 lme4 (97.6%), 465/508 mgcv (91.5%) - the *weakest* of the six
packages still clears 91%, and the two most authoritative samples (base,
stats) both clear 98.7%. **Resolution (R3 CONFIRMED unchanged):** dbarts
belongs to the base-R tradition of terse, lowercase, unpunctuated error
fragments that lean on R's own `Error in call: ...` framing to supply the
sentence's subject and full stop, not the tidyverse/rlang/cli tradition of
self-contained, fully punctuated, multi-bullet messages - and a package with
no rlang/cli dependency (and no plan to add one) cannot half-adopt the
tidyverse convention, because full-sentence punctuation is only legible as
its *own* tradition when paired with the bulleted, hint-terminated structure
it was designed for; borrowing the period without the bullets produces a
message that reads as neither house style, and every highly-regarded
base-style package sampled - independently, consistently, above 91% in the
weakest case - has made the same choice dbarts's own corpus already made.
R2 and R3 are CONFIRMED together on this reasoning.

## R4. ASCII only — CONFIRMED, reinforced

Restates the house rule (`CLAUDE.local.md`): `-` not an en/em dash, `->` not
an arrow glyph, no smart quotes. Measured 0 non-ASCII bytes in any message
body today - already fully conformant; stated here only because the sweep is
about to retype ~120 of these strings by hand, the easiest way to reintroduce
one.

**External evidence.** Neither the tidyverse guide, rlang docs, nor WRE state
an ASCII rule (this is a house convention, not an ecosystem one - published
best practice is silent). Observable practice, however, corroborates it
independently: a byte-level scan of all ~3,400 sampled `stop()` calls in
base, stats, Matrix, survival, lme4, and mgcv found 0 lines with any
non-ASCII byte in any of the six packages. **R4 is CONFIRMED unchanged**,
now on two legs instead of one: it was already dbarts's own house rule, and
it also turns out to be exactly what R-core's own error-message corpus does.
(See also the sQuote/dQuote note under R1: the one place base R's own
machinery *could* emit non-ASCII output - fancy quotes - is a runtime,
locale-dependent rendering choice, not literal source text, and dbarts's
avoidance of `sQuote`/`dQuote` sidesteps it entirely.)

## R5. Placeholder / value formatting — SILENT-KEPT

R: interpolate via `stop()`'s own `...` concatenation, with the placeholder's
quote marks (R1) typed directly into the adjacent string literal - the
majority form (43 of ~120 sampled interpolating calls) over `sprintf()` (15).
Keep `sprintf()` as the accepted alternate when one clause interpolates two or
more values (`R/augmentation.R:18`); do not introduce `paste0()`/`gettextf()`
(0 uses of either today). Conformance: `stop("'chainNum' must be a single
chain index in [1, ", n.chains, "]")` (`R/generics.R:3030`). The example
this rule first cited, `"invalid monotone direction '", value, "'; use
-1, 0, or +1"`, did not survive: R12's rewording left `R/model.R:557-560`
a fixed two-literal message that interpolates nothing.

C: `Rf_error`'s only mechanism is its own printf placeholders (`%s`, `%d`,
`%zu`) - no alternative exists. Quote `%s` in `'...'` when it echoes a name or
a user-supplied choice; leave it bare for a descriptive phrase (`"a treatment
forest does not support %s"`, `src/R_interface_bartcore.cpp:2340`, where
`refused` is a phrase like `"a DART tree prior"`, not a name).

**External evidence.** Neither tidyverse nor rlang docs address the choice
between comma-concatenation and `sprintf()` for a base `stop()` call - that
choice only exists outside the rlang/glue-string world they're written for
(rlang messages are glue-interpolated, a mechanism dbarts does not use and
this document does not propose). WRE is silent on interpolation mechanism
entirely. Observable practice is not silent, but it answers a different
question than the one this rule asks: `gettextf()` is the dominant
interpolation mechanism in base (161 hits), stats (110), Matrix (82), and
lme4 (33) - overwhelmingly paired with `domain = NA`, which marks the string
for translation via R's own message-catalog (`.po`) machinery. That is an
internationalization infrastructure decision (a translatable-message
catalog, `po/` directory, `tools::update_pkg_po()`), not a text-style
decision, and it is out of scope for this document by the same
translatability filter used throughout: dbarts ships no message catalog and
this document does not propose adding one. survival, notably, uses
`gettextf()` almost never (1 hit) and leans on `paste()` (35 hits) instead -
so even setting the i18n question aside, there is no single external
convention for *which* interpolation mechanism to use, only for *whether* to
route it through a translation catalog. **R5 is SILENT-KEPT**: the in-repo
majority (comma-concatenation, `sprintf()` for the multi-value case) stands.

## R6. Multi-clause shape: main clause + one explanatory/remedy clause, max — CONFIRMED

One main clause stating the refusal, optionally followed by ONE more clause
(`:` for an explanation, `;` for a remedy) - not both.

- Conformance: `"%s: a multi-forest sampler fixes its data at creation; make
  a new sampler instead"` (`src/R_interface_bartcore.cpp:2642`) - main +
  one remedy clause.
- Both prior violations were fixed while they stood: `refuseHostMutation`'s
  message was refusal plus one remedy clause (the function and every call site
  are since deleted, multinomial-mutation-arc.md S4), and `R/spec.R:54-59`
  (probit weights) collapsed to the same shape, matching the C twin's
  already-shorter form (`src/R_interface_bartcore.cpp:2784-2786`).

**External evidence.** Published best practice's general philosophy agrees
with R6's spirit: "An error message should start with a general statement of
the problem then give a concise description of what went wrong" - a lead
clause plus one (or a small number of) elaborations, not a sprawl. But its
*mechanism* for the elaboration is a multi-line bulleted list (`x`/`i`
markers, an optional hint bullet) - a structure that has no plain-`stop()`
equivalent and is explicitly out of scope by the translatability filter (a
bulleted list concatenated into one string reads as a run-on, not as
bullets). With the structural half of the tidyverse rule excluded, what's
left to check against is observable practice, which independently supports
brevity: average length of the first string literal across the ~3,400
sampled calls was 36-40 characters per package (base 36.5, stats 36.6, Matrix
38.3, survival 36.4, lme4 36.6, mgcv 39.5) - short, and a spot check of the
samples (`sample_*.txt`) shows the overwhelming majority are single-clause.
**R6 is CONFIRMED unchanged**: the general concision principle is reinforced
by both sources; the specific "max one extra clause" threshold remains
dbarts's own reasonable translation of that principle into a single string,
since no external source dictates an exact clause count for a non-bulleted
message.

## R7. C-side caller-name prefixing: three call-site kinds, one convention each — SILENT-KEPT

Measured: of 351 `Rf_error` calls, 259 (74%) carry no caller-name prefix, 69
(20%) carry a dynamic `"%s..."` prefix fed a `caller` parameter, 23 (7%) carry
a hardcoded literal entry-point name. These track three distinct call-site
kinds, confirmed by reading the call graph, not three competing dialects at
the same sites:

- **Default - no prefix.** Raised directly in a `.Call` bridge entry point's
  own body (or a helper it alone reaches); R's `Error in .Call(...)` frame
  already attributes it. E.g. `"forest weight length must match the number of
  observations"` (`src/R_interface_bartcore.cpp:4072`).
- **Flat C API entry points** (`src/C_interface.cpp`, `src/R_interface.cpp`,
  the `dbarts.h` ABI) **- hardcoded literal self-name.** Reachable by a
  `LinkingTo: dbarts` consumer with no R call frame to consult, so the message
  must self-identify. All 23 hardcoded instances are a `dbarts_sampler_*`
  function naming itself (`src/C_interface.cpp:534`: `"dbarts_sampler_run:
  results.structSize is 0..."`). Always follow the name with `:`.
- **Shared bridge helpers reached from multiple `.Call` entry points -
  dynamic `"%s: ..."`.** `caller` can't be hardcoded since it varies per call
  site (`refuseMultiForestMutation`, called from `bartcore_setData` and
  `bartcore_setModel`, `src/R_interface_bartcore.cpp:4908`/`:5214`).
  Standardize on `%s: ` (colon) - the majority sub-style at the survey
  date (38 of 69 vs 31 without), e.g. `:2816`. Reword the no-colon
  instances to add the colon. The one cited here has been: `:133` now
  reads `"%s: requires a numeric matrix with matching columns"`, and the
  `Rf_error("%s", ...)` calls that remain in that file substitute a whole
  pre-built message rather than a caller prefix.

A helper reached from both the .Call bridge and the flat C API takes the
dynamic prefix and is handed the flat entry's own name, since that route
has no R call frame at all; three refusals are the outstanding exceptions,
prefixed by none of them - the factor level-code texts and the level-count
ceilings (src/R_interface_bartcore.cpp:1793-1800, :1826-1829), raised from
the shared validateCategoricalPredictors and
validateTestContainerAgainstStore, and the sparse leaf-covariate test
refusal raised in dbarts_sampler_setTestPredictors' own body
(src/C_interface.cpp:891).

**External evidence.** WRE's C-API chapter has a section titled "Error
signaling" in its table of contents, but the fetched manual text did not
include that section's body, and the surrounding chapters that *were*
fetched contain no guidance on caller-name prefixing or on identifying the
raising function in a C-level message (this is flagged, not glossed over -
see the appendix for exactly what was and wasn't reachable). None of the six
sampled R packages ship a public, `LinkingTo`-callable flat C ABI header
comparable to `inst/include/dbarts/dbarts.h` reachable with no R call frame
on the stack, so there is no comparable observable practice to check R7
against - this is architecture dbarts alone has among the sampled set, not a
convention any of them chose for or against. **R7 is SILENT-KEPT**: the
in-repo, call-graph-derived three-kind convention stands as written.

## R8. Cross-language policy: R canonical, C backstops deliberately distinct — SILENT-KEPT

When the same rule is guarded on both sides (R validates, the bridge
re-validates as a backstop for direct `.Call`/C-API callers bypassing R), the
R message is canonical - the one users see, the one docs/tests quote; the
C-side message is independently worded against this rule, not copied.
Verbatim-identical strings across the boundary already failed once: `"forest
weights must be finite and non-negative"` is byte-identical at
`R/dbarts.R:1435` / `src/R_interface_bartcore.cpp:4079`, but the very next
guard in the same pair, the length check, drifted apart with no R counterpart
to stay in sync with (`"forest weight length must match the number of
observations"`, C-only). A string shared across two languages by hand is
identical by discipline, not by construction, and the discipline already
lapsed once. Every C-side backstop should carry a one-line comment naming it
as direct-API defense in depth, on the model of `src/C_interface.cpp:788-794`
(`"defense in depth, since validateTestSource has already raised it"`). This
is the message policy K6 reconciles predicates toward; K6 remains free to
decide, guard by guard, whether a given C-side backstop is worth keeping.

**K6 addendum (predicate parity).** Ownership of the *message* (above) is a
separate question from parity of the *predicate*: whichever side is not
canonical must still accept and reject the identical input set, or a value
can cross the boundary having passed one guard's wording only to be caught
(or worse, missed) by the other's. This is not hypothetical: `dbartsSampler$
setForestWeights` checked `anyNA(w) || any(w < 0)` R-side while the bridge
checked `!R_FINITE(w) || w < 0.0` - an `Inf` forest weight passed the R guard
and was refused three frames deep in a `tryCatch` rethrow instead of cleanly
at the surface. Two rules follow. First, the C-side guard stays even when the
R-side front door already makes it unreachable on the normal path: dbarts
exposes low-level `bartcore*` `.Call` entries (see `R/bartcore.R`'s header)
that construct and mutate a sampler directly, skipping every R-level check,
so the C guard is load-bearing for that route regardless of what the front
door does. Second, an R-side predicate touching a value that can be `NA`/
`NaN` must be written to refuse it deterministically in one guard -
`weights < 0.0` alone is not equivalent to the bridge's `!(w >= 0.0)`,
because `NA_real_ >= 0.0` is itself `NA` in R and `if (NA)` raises an
unnamed error rather than the intended message; write `is.na(w) | w < 0.0`
(or `!is.finite(w) | w < 0.0` when the bridge also requires finiteness), not
a separate `anyNA` pre-check relying on evaluation order to paper over it.

**External evidence.** None of the four external sources (tidyverse, rlang,
WRE, the six-package survey) address a package with dual-language validation
of the same argument, an R layer plus an independently linkable C ABI - this
is, again, dbarts's own architecture, not a question the R error-message
ecosystem has occasion to answer. **R8 is SILENT-KEPT**: the draft's
recommendation stands as the final rule (this also closes open point 2
below - see there for the explicit sign-off).

## R9. Missing required argument — CONFIRMED (both sub-shapes), one open point now closed

Two sub-shapes, tracking two different R predicates. **Value present but
NULL** (`is.null(x)` / C `ptr == NULL`): `"'<name>' cannot be NULL"`.
Conformance: `"x.test cannot be NULL"` (`R/dbarts.R:1090`). **Argument omitted
from the call** (R `missing(x)`, no C analogue): `"'<name>' must be
specified"`. Conformance: `"'group.by' must be specified to use rbart_vi"`
(`R/rbart.R:127`). The former violation, `"'group.by' must be supplied when
'newdata' is given"`, is retired: `R/bart.R:2605` now refuses a different
predicate, `"'group.by' must be given by name when 'newdata' is given"`.

**External evidence, omitted-argument sub-case.** Neither tidyverse nor rlang
prescribe a specific verb here. Observable practice across the six packages
tallies `"must be specified"` 12 times (base 7, stats 1, lme4 1, mgcv 3),
`"must be given"` 6 times (stats 4, survival 2), `"must be supplied"` 3 times
(stats 1, mgcv 2) - `specified` is the external plurality too, and it's
attested repeatedly in base itself (`base::seq.Date` x3, `seq.POSIXt`,
`globalCallingHandlers`, `tryCatch`). **CONFIRMED, reinforced** - `specified`
was already the in-repo choice; external practice independently agrees.

**External evidence, NULL-present sub-case (was open point 1).** This is
where the external record has nothing to say. The literal phrase `"cannot be
NULL"` occurs 0 times across all six sampled packages (~3,400 calls). Related
phrasings exist but don't converge on one: `lme4::getStart`/`updateStart` use
`"'start' is not NULL, a numeric vector, or a list"` (a different predicate -
enumerating an allowed set that happens to include NULL, not rejecting a bare
NULL); `mgcv::uniquecombs`/`uniquecombs0` use the bare, unquoted, informal
`"x is null"` (2 hits, one package, lowercase `null`). Neither published
guidance nor observable practice speaks clearly enough to move the needle.
**Open point 1 is CLOSED, not merely recommended: external evidence is
silent, so the draft's own tie-breaking reasoning governs** - `"'<name>'
cannot be NULL"`, matching the `must be`/`cannot be` spine every other rule
in this document uses. This sub-shape of R9 is final.

## R10. Type/class mismatch — OVERRIDDEN (enrichment added, base shape unchanged)

`"'<name>' must be a/an <type description>[, not <actual>]"` - the bare shape
without the bracketed clause is dominant already (~25 sampled type-mismatch
messages, only 2 use a different verb) and stays the required minimum.
Conformance: `"'weights' must be a numeric vector"` (`R/data.R:765`,
`R/data.R:1290`) - both guards in the same file now read exactly alike.

The bracketed `, not <actual>` clause is new: append it when the actual
type/class is already in hand as a short noun (`class(x)[1]`, `typeof(x)`) at
essentially no extra cost - encouraged, not mandated (see R11/R12 for why the
same qualifier applies to all three "state the actual value" enrichments
uniformly).

**External evidence and the merit case.** The enrichment is adopted ON ITS
OWN MERITS, not on any guide's authority: a type-mismatch report that names
what was actually received converts "what did I pass?" from a debugging
round trip into a one-read diagnosis, at essentially zero cost when the
class is already in hand at the call site - that diagnostic value is the
whole argument, and it stands without citing anyone. Genuine base-style
precedent exists and shows the clause reads naturally in this tradition:
`stop("'data' must be a data.frame, not a matrix or an array")`, used
twice, verbatim, in stats (`get_all_vars`/`model.frame.default`) - real,
but rare (2 hits in ~3,400 calls), which is why the clause is encouraged
rather than required. (The tidyverse guide recommends the same "must be X,
not Y" shape; under the project lead's direction that recommendation
carries no weight of its own and is recorded here only as corroboration.)
**Verdict: OVERRIDDEN, on merits.** The base shape (bare `"must be a/an
<type>"`) is unchanged, because that's what both the majority of dbarts's
own corpus and the majority of the external corpus already write; the
enrichment is "encouraged when in hand" - not mandatory (observable
practice doesn't demand it and slice L's budget doesn't afford re-deriving
`class(x)` at every site that lacks it already), but no longer merely a
maybe.

- Old: `"'<name>' must be a/an <type description>"` only.
- New: same, `+ ", not <actual>"` when the actual type is already available.

## R11. Length mismatch — OVERRIDDEN (mandatory `got` dropped to encouraged)

Against a fixed/derived count: `"'<name>' must have length <N>"`. Against
another argument's length: `"'<name>' must have the same length as
'<reference>' (<N>)"`. Closest precedent (names the reference but not the
numbers): `sprintf("'%s' must have length %d, that of '%s'", name, n,
reference)` (`R/augmentation.R:18`) - keep its shape, the sweep's one named
exception. The prior furthest outlier, `R/bartcore.R`'s `"length of new x
does not match old"` (naming neither argument nor either length), is fixed:
`R/bartcore.R:284` now reads `"'x' must have length <N>"`, matching this
rule.

Appending the actual (got) length - `"; got <got>"` - is now encouraged, not
mandated, when the value is already in hand; see the shared reasoning under
R12, which found the identical pattern.

**External evidence.** The draft's own R11 already flagged this as invented
("no existing [in-repo] template names both the expected length and the
actual (got) length"), so there was no in-repo majority to protect here in
the first place - this rule was always going to be set by outside judgment,
and under the evidence precedence that judgment is the measured base-style
practice. Observable practice answers the length question directly and
says no: of ~90 length-related
messages sampled across all six packages (`"must have [the same] length"`,
`"length mismatch"`, `"lengths ... match"`, etc.), not one names the actual
received length - representative: `stop("'x' and 'y' must have the same
length")` (`base::formatDL`), `stop("'model$order' must be of length 3")`
(`stats::arima.sim`), `stop(sprintf("length mismatch in %s (%d != %d)", nm,
length(x), expected_len))` (`lme4::setParams` - the one place a package
*does* echo two numbers, and even there it echoes both current *and*
expected, not a bare "got"). A direct grep for any "length ... got" pattern
across the full corpus returned zero hits. **Verdict: OVERRIDDEN** - the
invented mandatory `"; got <got>"` suffix is dropped as a requirement; the
expected-length-only shape (already what most of dbarts's own corpus does
where it names a number at all) becomes the required minimum, with the
`got` value permitted and encouraged when cheap, exactly mirroring R10 and
R12. This also folds in and generalizes what was open point 3 (out-of-range
enrichment scope) - see that entry below.

- Old (invented): `"'<name>' must have length <N>; got <got>"` (mandatory).
- New: `"'<name>' must have length <N>"` required; `"; got <got>"` encouraged
  when the value is already in hand, never required.

## R12. Enum/choice rejection — OVERRIDDEN (shape confirmed, mandatory `got` dropped)

`"'<name>' must be one of <c1>, <c2>, ..."` - retires all four competing
adjectives (`invalid`/`unknown`/`unsupported`/`unrecognized`, no majority:
4/2/5/3 sampled hits, and `unsupported` is already claimed by R13). `"must be
one of"` matches base R's own `match.arg()` wording, which already governs 30
enum checks here with no hand-written message at all - the closed-choice
family's real majority convention is "let `match.arg` say it"; this template
is for checks that can't use it (non-character enums, C-side class dispatch).
Appending `"; got '<value>'"` is now encouraged, not mandated (see below).

Conformance (shape): `"'forest' must name one of '", paste0(..., collapse =
"', '"), "'"` (`R/generics.R:645-649`). The monotone-direction violation is
retired: `R/model.R:555` now reads `"'direction' must be one of -1, 0, 1"`,
with no `got` value - this rule's shape exactly. Still open:
`"unrecognized response family for a binary response"`
(`src/R_interface_bartcore.cpp:1744`) - names no choices at all.

**External evidence.** The `"must be one of"` shape itself is directly
attested in base R: `stop("'origin' must be one of 'start', 'current' or
'end'")` (`base::seek.connection`); also `stop("'shape' must be one of \"g\",
\"t\", \"s\"")` (`Matrix::.sparseDiagonal`), `stop(gettextf("view variables
must be one of %s", ...))` (`mgcv::vis.gam`) - **CONFIRMED** for the base
shape, now with a base-R precedent the original draft didn't have (it
inferred the shape only from `match.arg()`'s wording, not from a hand-written
`stop()` using it). The four competing rejected adjectives remain genuinely
mixed externally too (`"unknown"` and `"unrecognized"` both appear repeatedly
across the six packages, e.g. `stats::density.default`: `"unknown bandwidth
rule"`; `survival::cox.zph`: `"Unrecognized transform"`) - no external
convention picks a winner among them either, so retiring all four in favor of
`"must be one of"` remains the right call. But none of the three external
`"must be one of"` examples just quoted appends the value that was actually
received - not `seek.connection`, not `.sparseDiagonal`, not `vis.gam`. A
targeted search across the full six-package corpus for any `"must be one
of"`/`"should be one of"` call that also echoes the offending value found
zero. **Verdict: OVERRIDDEN on the `got` clause only** - same reasoning and
same resolution as R11: mandatory got-value is dropped, `"must be one of
<choices>"` is the required minimum, `"; got '<value>'"` is encouraged, not
mandated, when cheap.

- Old: `"'<name>' must be one of <c1>, <c2>, ...; got '<value>'"` (mandatory).
- New: `"'<name>' must be one of <c1>, <c2>, ..."` required; `"; got
  '<value>'"` encouraged when the value is already in hand, never required.

## R13. Unsupported-under-family/composition refusal — SILENT-KEPT

`"<subject> does not support <feature>[: <reason>]"` - clear majority (26 of
37 R + 16 C "not support" hits use this exact verb). Conformance: `"a
treatment forest does not support %s"`
(`src/R_interface_bartcore.cpp:2340`); `"probit models do
not support weights; fit integer count weights with family = \"logistic\",
or model continuous weights' latents directly"` (`R/spec.R:54-59`).
Violation, since reworded to exactly this rule's proposal: `"sample =
\"test\" is not available for type = \"forest\": an ..."` is now
`"type = \"forest\" does not support sample = \"test\": no test-sample
per-forest channel is stored, ..."` (`R/generics.R:633-638`).

Sub-shape, foreign-argument refusal: an argument a method does not take,
refused because a sibling method of the same surface does take it, reads
"'<name>' is not used by <generic> on a <class> fit: <reason>" - the
subject is the argument, not the method, so the caller reads the
offending name first (R/generics.R:2062-2075). Where the refusal has no
single offending name, the plain R13 shape stands: "predict on a bart fit
does not support unnamed arguments: ..." (R/generics.R:2086-2095).

**External evidence.** Checked hard, because this looked like a plausible
override going in. Neither tidyverse nor rlang address "unsupported feature"
phrasing at all - published best practice is silent. Observable practice is
not silent, but it's genuinely split and, if anything, leans the other way on
raw frequency: across the six packages, the passive constructions `"not
available"` (46 hits: base 2, stats 8, survival 9, lme4 2, mgcv 25) and `"not
supported"` (31 hits: base 4, stats 2, survival 17, mgcv 8) both outnumber
the active `"does not support"` (2 hits total: `stats::logLik.lm`,
`survival::survreg`). That's a real signal, but not a decisive one: no
package has an internal majority for any single verb (survival alone has both
17 `"not supported"` hits and its own `"does not support"` instance), no
published rule addresses the choice, and in-repo already has an unusually
strong, self-consistent 42-hit majority for the active form across both R and
C. Active-vs-passive voice for this one phrase is exactly the kind of
question the task's "split -> in-repo majority stands" clause is for.
**R13 is SILENT-KEPT**: `"does not support"` remains the rule, on the
strength of the in-repo majority, with the caveat now on record that
`"not available"`/`"not supported"` are the more common constructions in the
wild if a future revision wants to reopen this.

## R14. Mutation-guard refusal — CONFIRMED, reinforced

R-canonical shape (per R8): `"<operation> is not available on <context>[:
<reason>]"` - was the shape `refuseHostMutation` gave 22 call sites, before
the function and every call site were deleted (multinomial-mutation-arc.md
S4). The shape survives in its replacements, `refuseCountsMutation`
(`R/bartcore.R:62`) and `refuseAmplitudeMutation` (`R/bartcore.R:36`), and in
`$setControl`'s own per-slot refusal (`R/dbarts.R:1207-1213`). C-side
backstops keep their own idiom under R7/R8 (`"%s: a multi-forest sampler
fixes its data at creation; make a new sampler instead"`,
`src/R_interface_bartcore.cpp:2642`). Conformance: `setModel`'s DART-tree-
prior refusal (`R/dbarts.R:1249-1256`) reads `"changing a DART tree prior is
not available on an existing sampler: recreate it instead"`, matching this
rule.

**External evidence.** R14 happens to already use the phrase the R13 survey
just found to be the *more* common external construction: `"not available"`
was attested 46 times across the six packages, more than `"does not
support"`'s 2. That's a coincidental but welcome alignment - R14's existing
wording needed no change to already match the more common external idiom for
this family of refusal. **R14 is CONFIRMED unchanged**, reinforced by the
same data set that left R13 SILENT-KEPT.

---

## Measured current state (appendix)

Corpus: 610 `stop()` (`R/`), 365 `Rf_error()` (`src/R_interface_bartcore.cpp`
+ `src/C_interface.cpp` + `src/R_interface.cpp`), 3 `ext_throwError()`
(`src/bartcore/chain.hpp:1886,2050` and `src/bartcore/sampler.hpp:682`, the
only users of that fourth mechanism, not part of this drift). 0 uses
anywhere of `gettextf`,
`sQuote`, `dQuote`, backtick-quoted argument names, or non-ASCII bytes.

| dimension | finding |
|---|---|
| argument quoting | single-quote opens the message (177 R, 14 C); backtick/sQuote/dQuote: 0/0 |
| terminal period | period-free: 348/348 C, 549/549 R - the sweep's three restated-refusal exceptions are gone |
| sentence case | lowercase-initial: 534/534 R bodies that carry a literal (15 more rethrow a caught error object and carry none), 347/348 C; C's one remaining outlier is the DART acronym exception - the prior genuine miss is fixed, and the two BCF-initial messages were reworded with the amplitude rename |
| hyphenation | `non-negative` 24 (R 16 + C 8) vs `nonnegative` 0 - `non-negative` is now unanimous |
| out-of-range shape | `"... out of range"` 42 combined (R 8 + C 34) vs `"must be between"` 1 - bare `out of range` remains majority for index/discrete bounds |
| composition-refusal verb | `"does not support"` 33 (R 28 + C 5) vs `"incompatible with"` 3 (C) / `"is not available for"` 0 |
| missing-arg (NULL sub-case) | `"cannot be NULL"` 13 (R 8 + C 5) vs `"is NULL"` 4 (R 3 + C 1) - `cannot be NULL` is now the clear majority (settled by R9) |
| missing-arg (omitted sub-case) | `"must be specified"` 5, `"must be supplied"` 0, `"must be given"` 0 - `specified` is now unanimous |
| enum-rejection verb | messages opening on `invalid` 2 (C), `unknown` 3 (R), `unsupported` 2 (R 0 + C 2, claimed by R13), `unrecognized` 5 (R 2 + C 3); `"must be one of"` 2 - still no majority among the four rejected adjectives, migration to the settled shape is gradual |
| C caller-name prefix | no-prefix 259 (74%), dynamic `%s`+caller 69 (20%, colon 47 / no-colon 22), hardcoded literal 20 (6%) - all three map to distinct call-site kinds (R7) |
| length-mismatch templates | unified: 1 template (R11's `'<name>' must have the same length as '<reference>'` / `'<name>' must have length <N>`), applied at 14 call sites across `R/data.R`, `R/A_class.R`, `R/dbarts.R`, and `R/partialDependence.R`; `R/augmentation.R:18`'s `sprintf()` form is R11's one named, kept exception. Was 6 competing templates at this table's prior revision, drifted to 7 before this sweep |
| R interpolation | comma-concatenated `stop()` args majority (43 sampled) vs `sprintf()` (15); `paste0()`/`gettextf()` as the sole message wrapper: 0 (inline `paste0(..., collapse = ...)` for an enum/choice list, e.g. R12's own template, is not this and is unaffected) |

## External evidence appendix

Method: (1) fetched the tidyverse style guide's error-messages chapter and
rlang's `abort()` reference documentation directly; (2) fetched Writing R
Extensions (CRAN r-release build) and searched it for message-style guidance,
including a direct attempt at its "Error signaling" section; (3) extracted
every `stop()` call's source text from six installed, highly-regarded
base-style packages by walking each package namespace and deparsing every
function body (not `grep` over installed sources, which are lazy-loaded
databases with no plain-text `.R` files on disk for an installed package) -
script and raw per-package call dumps are in the scratchpad
(`extract_stops.R`, `analyze_stops.R`, `stops_<pkg>.txt`,
`sample_<pkg>.txt`). Versions: R 4.6.1; Matrix 1.7.5; survival 3.8.6; lme4
2.0.1 (Rcpp/RcppEigen-backed, no rlang/cli message dependency); mgcv 1.9.4;
rlang 1.2.0 (docs only, not adopted). All package/version pairs are from the
CRAN builds installed in this environment's R library.

1. **Tidyverse style guide, error-messages chapter**
   (https://style.tidyverse.org/errors.html, fetched 2026-08-17). Governs
   rlang/cli-based packages; source of the structure this document
   deliberately does not adopt wholesale.
   > "An error message should start with a general statement of the problem
   > then give a concise description of what went wrong... If the cause of
   > the problem is clear (e.g. an incorrect type or size), use 'must': `n`
   > must be a numeric vector, not a character vector... Errors should be
   > written in sentence case, and should end in a full stop... Surround the
   > names of arguments in backticks, e.g. `x`."

2. **rlang `abort()` reference** (https://rlang.r-lib.org/reference/abort.html,
   rlang 1.2.0, fetched 2026-08-17). Confirms what does *not* translate to
   plain `stop()`: the bullet-marker system and the message/body/footer
   split have no base-R equivalent and are excluded by this document's
   translatability filter.
   > "Elements named `\"*\"`, `\"i\"`, `\"v\"`, `\"x\"`, and `\"!\"` are
   > formatted as regular, info, success, failure, and error bullets
   > respectively... The first element is displayed as an alert bullet
   > prefixed with `!` by default."
   The docs give no prose guidance on phrasing, capitalization, or
   punctuation beyond the bullet mechanics themselves - that guidance lives
   in the tidyverse style guide (source 1), not here.

3. **Writing R Extensions** (CRAN r-release,
   https://cran.r-project.org/doc/manuals/r-release/R-exts.html, fetched
   2026-08-17). Searched for style guidance on `stop()`/`error()`/`warning()`
   message composition, including a direct fetch attempt at the "Error
   signaling" anchor (present in the manual's table of contents but not
   retrievable as fetched text in this environment - flagged as **not
   reachable**, not glossed over). Conclusion, stated plainly per the
   instruction not to invent guidance that isn't there: **WRE contains no
   style guidance on error-message wording, capitalization, punctuation, or
   quoting anywhere in the sections that were reachable.** Its message-related
   content is confined to translatability/internationalization mechanics
   (`gettext`, message catalogs), which is why R5's `gettextf()` question is
   resolved as an infrastructure decision, not a style one.

4. **base + stats** (R 4.6.1, installed base distribution). 602 + 927 =
   1,529 `stop()` calls extracted; 585 + 917 = 1,502 carried a literal
   string. The single most authoritative sample - base R's own error-message
   corpus governs every other rule in this document by extension, since
   dbarts's whole house style already imitates it.
   > `stop("'x' and 'y' must have the same length")` (`base::formatDL`)

5. **Matrix 1.7.5**. 199 `stop()` calls, 193 with a literal string. Chosen
   as a highly-regarded numerical/linear-algebra base-style package, the
   closest analog to dbarts's own domain among the four candidates.
   > `stop("'shape' must be one of \"g\", \"t\", \"s\"")`
   (`Matrix::.sparseDiagonal`).

6. **survival 3.8.6**. 926 `stop()` calls, 922 with a literal string. The
   package with the weakest (but still clear) majority on both R2 and R3,
   and the most frequent user of Title Case among the six.
   > `stop("Argument lengths do not match")` (`survival::coxph.wtest`).

7. **lme4 2.0.1**. 224 `stop()` calls, 211 with a literal string.
   > `stop(sprintf("length mismatch in %s (%d != %d)", nm, length(x),
   > expected_len))` (`lme4::setParams`) - the one package-internal instance
   found anywhere in the six-package survey of a length-mismatch message
   naming both the expected and the actual length.

8. **mgcv 1.9.4**. 511 `stop()` calls, 509 with a literal string. The
   heaviest external user of backtick-quoted names (25 hits) and of
   terminal periods (43 hits) among the six - still a minority within its
   own corpus on both counts, but the closest any sampled package comes to
   the tidyverse convention.
   > `` stop(gettextf("view variables must be one of %s", paste(v.names, `` \
   > `` collapse = ", "))) `` (`mgcv::vis.gam`).

## Open points for sign-off

1. **CLOSED. Missing-argument NULL-case verb: `"cannot be NULL"` vs `"is
   NULL"`.** External evidence is silent (see R9): 0/6 sampled packages use
   `"cannot be NULL"` or any single converged alternative for this exact
   predicate. With no external signal to weigh, the draft's own tie-breaking
   reasoning stands as final: `"cannot be NULL"` - matches the `must be`/
   `cannot be` spine every other rule in this document uses, rather than `is
   NULL`'s descriptive voice, which reads oddly once R9-R14 are all phrased
   as requirements.

2. **CLOSED. Cross-language canonical side (R8).** External evidence is
   silent - no sampled package has a comparable dual-language (R +
   independently linkable C ABI) validation architecture to check this
   against (see R8). The draft's recommendation is adopted as final: R
   canonical, C-side deliberately distinct and comment-justified, over
   verbatim-identical strings kept in sync by hand. This gates how K6 is
   written; it now has an explicit yes.

3. **CLOSED, and generalized. Out-of-range enrichment scope (R11).** The
   draft's own instinct here - "bare is the required minimum... enrichment
   is encouraged, not mandated" - turned out, once R10 and R12 were checked
   against external evidence, to be the *general* answer for every "should
   this message echo the actual/received value" question in this document,
   not just this one. Echoing the received value has real diagnostic merit
   (it is the enrichment's own case, R10); observable practice across six
   highly-regarded base-style packages (~3,400 `stop()` calls) essentially
   never does it (0 instances for length or enum rejection, 2 for type).
   Rather than resolve that tension differently in four different places -
   merit says sometimes, tradition says rarely - R10, R11, and R12 all
   now use the identical resolution this open point already proposed for
   out-of-range: state the expected value/type/choices as the required
   minimum; append the actual/received value when it's already in hand,
   never as a mandatory re-derivation. `"'forest' index must be between 1
   and <N>"` (`R/generics.R:613`) remains the model example; it was right
   the first time.
