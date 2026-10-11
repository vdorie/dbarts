# stan4bart-creation-mapping: stan4bart fixes its response range at creation, from a fit holding the forest's columns

Status: PLANNED (dec-B437, dec-B441, dec-B444 to dec-B447; dec-B442 and dec-B443 are revised by dec-B446). No open
call. Revised after the blind critique. The code is stan4bart's; this file is its record.

agent: one opus implementer (R; two flags in C++); one opus reviewer who runs the mutants.
rng: POSTERIOR-CHANGING for every continuous stan4bart fit (the leaf prior's center and spread change, and one
range replaces one per chain). SHIFTING for a binary fit (its starting values). NEUTRAL for dbarts: no dbarts
code changes, so no dbarts baseline, snapshot or exact gate moves.
window: before the merge to main (dec-B437). Either side of [response-scale-rows.md](response-scale-rows.md);
see [Against response-scale-rows](#against-response-scale-rows).
budget: ~1050 lines, almost all in stan4bart (R ~330, C++ ~6, tinytest ~470, help and NEWS ~80, benchmarks
~110; records here ~55). Stop at 2 times; steps 1 and 2 can land as their own slice.

## Goal

A continuous stan4bart fit has one response range for its forest, computed once before sampling, the same for
every chain and never re-derived: by default from a linear mixed model that holds the forest's columns as linear
terms, fitted by lme4 where the model has random terms, and never wider than the response's own range; else as
`bart_range` says. The fit records the pair and how it was obtained. A build made beside lme4 runs without it.
The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Why, with the measurements: [stan4bart-response-range.md](../design/stan4bart-response-range.md). Ruled, each
covering what it names: dec-B437 (one range, fixed at creation, shared, never re-derived); dec-B441 (the forest's
columns enter linearly); dec-B444 (where the fit cannot be computed the range is the response's own, less any
offset); dec-B445 (`bart_range`: a pair, "mixed", "lm" or "response"); dec-B446 (lme4 stays a suggested package
and does the mixed-model fits; without it, a warning and a fallback, for continuous and binary alike; nothing of
ours fits a mixed model); dec-B447 (the default is never wider than the response's own range, and the help lists
no cases). dec-B419's restore of one sampler per chain stays, with dec-A197's calls.

"Run" is this plan's or its critique's: stan4bart's bartcore branch at 3d2a295 against a build of dbarts bartcore
1c64c477, arm64 macOS, R 4.6.1, lme4 2.1.0, on a shared machine.

Today, read in stan4bart's source at that commit and run:
- `stan4bart()` evaluates the call again as `lme4::lmer` (binary: `glmer`) on the formula less its `bart()`
  term; without lme4, or where that stops, `lm` (`glm`) with the grouping terms as fixed effects; then the fixed
  terms alone. The fitted values are the forest's starting offset, on the probability scale for a binary fit.
  With the formula held in a variable both families stop ("invalid type (language)") after a full BART fit run
  inside the initial fit.
- The range is derived when each chain's sampler is created, through
  [`dbarts_sampler_setOffset`](../../inst/include/dbarts/dbarts.h) with its re-derive flag true, and through
  warm-up, the same entry with the flag true on a thinning schedule. Nothing else derives it.
- The restore puts each chain's sampler on its recorded pair by
  [`dbartsSampler$setResponse`](../../man/dbartsSampler-class.Rd) with `updateScale = TRUE` on a vector
  alternating the pair, then the real response with the range held.
- Weights and the offset are not cut to the kept rows when a row is dropped for a missing value: with 4 of 60
  rows dropped the response and the designs have 56 rows, the weights and offset 60. Weights then stop the fit
  ("weights must be of length equal to 56"); an offset fits misaligned.
- Two operators replace stan4bart's copies of lme4's and reformulas's helpers with those packages' own when the
  package is installed. A build made beside lme4 holds 21 and 24 of their functions and, with lme4 hidden, stops
  at the first call (`could not find function "chk.cconv"`).

Run, and used below:
- The install. A sampler built as stan4bart builds it takes a pair exactly from a made-up response (the lower
  end on the first row of positive weight, the upper on every other row): width identical in 1000 of 1000 random
  pairs, with zero weights on the first rows and with one row of positive weight; draws bit-identical to a
  sampler whose creation derived the pair under 5 of 5 leaf priors; held through `setOffset` with the range
  held. In a patched stan4bart every chain's state held the pair in 17 of 17 configurations and a reload
  replayed the stored fits to 1.7e-12 in 16 (gp leaves differ, as on today's build: not this slice's).
- lme4 on the parsed pieces. `lme4::mkLmerDevfun` and `lme4::optimizeLmer`, given the fixed columns with the
  forest's added, the random-effect structure cut to the rows of positive weight, and a frame of the response
  less the offset with its weights, equal `lmer` on the written-out formula (fixed effects within 2e-10, the
  own part on every row within 2e-9; levels left with no row get an effect of zero). `lme4::mkMerMod` refuses
  such a frame ("can't find formula"), so the fit is read from the criterion's environment. Rank-deficient
  columns stop it ("Downdated VtV is not positive definite"). One random effect per row, which `lmer` refuses,
  runs to a flat criterion. The probit twin (`mkGlmerDevfun`, `optimizeGlmer` twice around
  `updateGlmerDevfun`) equals `glmer` within 1.4e-5.
- With the two operators removed, a build made beside lme4 holds none of its functions and fits a correlated
  slope with a crossed term, and a binary model with nested terms, with lme4 hidden.
- Near zero residual degrees of freedom the fit's residual sd runs from 0.008 to 6.5 for a generating 1, also
  where dec-B447's bound does not take over (18 residual degrees of freedom and fewer, 20 seeds a cell).

## Design

A. The argument. `bart_range = c("mixed", "lm", "response")`, last in `stan4bart()`'s formals, checked before
the `iter = 0` return.
- A character value is matched as `offset_type` and `store` are (`match.arg`: abbreviations pass, anything else
  is its error); NULL is "mixed".
- A numeric value is a pair: two finite numbers, the first below the second, their difference finite; else
  `'bart_range' must be c(lower, upper), lower below upper, with a finite width`. It is never reordered.
- A pair that does not overlap the range of `y - offset - own` on the rows in (D) is used, with a warning:
  `'bart_range' (<lower>, <upper>) lies wholly outside the values BART is fitted to (<a> to <b>)`.
- For a binary response a value other than the default warns `'bart_range' is not used for a binary response`.

B. Rows and columns, from what `stan4bart()` has parsed.
1. Kept rows. Weights and the offset are cut to the rows the model keeps, where `stan4bart()` reads them: a
   defect fix, in its own commit.
2. Rows in: the kept rows of positive weight; n is their count. The fit and the range read no other row.
3. Added columns, from the forest's training matrix and its column kinds. A numeric or ordered column enters as
   the value the forest splits on. An unordered factor enters as indicators of its levels present among the rows
   in, less the first; unless every such level lies within one level of some grouping factor of the model, when
   it is left out and named (rules 1 and 3).
4. Aliasing. The fixed block is [1, the model's own columns, the added columns] on the rows in. A pivoted QR at
   `lm`'s tolerance (1e-7) keeps the first of any dependent set, so the model's own columns are kept before the
   added ones; the dropped are named (rule 2). p counts the kept columns. lme4 never chooses: it would stop.
5. Where no added column is left (each aliased or left out), the fit is the model's own; a verbose line says so:
   `no BART predictor enters the initial fit, so the range is that of the model's own terms; see 'bart_range'`.

C. The fits. z is the response less the offset on the rows in, w their weights.
1. Mixed (the model has random terms, "mixed", lme4 available by H.2): `mkLmerDevfun` on a frame of z and w,
   the kept columns, and the model's own random-effect structure with its design cut to the rows in, under
   restricted likelihood and `lmerControl`'s defaults with the convergence checks off; `optimizeLmer` without
   derivatives; the coefficients, predicted effects and residual sums read from the criterion's environment
   (`pp$beta`, `pp$b`, `pp$sqrL`, `resp$wrss`). Warnings are muffled; a singular fit and an optimizer's
   complaint change nothing (rule 10). Any error, or a non-finite value read, is a failure (E).
2. Before C.1, lme4's own refusal, which its builder does not make: a random term with at least n effects makes
   the fit undefined (E).
3. Linear ("lm", or "mixed" on a model with no random term): `lm.wfit` of z on the kept columns.
4. The own part, on EVERY kept row: the model's own columns centered at the means the parametric sampler centers
   them at, times their coefficients, plus the predicted effects through the full random-effect design (zero for
   a level with no row in; absent under C.3). The intercept and the added columns' part are left out.

D. Routes and starting values.
- "mixed" and "lm": the range is the smallest and largest of `y - offset - own` over the rows in. Where its
  width exceeds that of `y - offset` over the rows in, the response's own range replaces both ends (dec-B447;
  E). A pair is never bounded.
- "response": own is zero and the range is that of `y - offset` over the rows in. Nothing is fitted.
- A pair: the range is the pair; own is "lm"'s with its fallbacks, for the start only. No mixed model is fitted.
- Starting values: the forest's first offset is own, the user's offset being added once, where the sweep loop
  adds it (rule 7); the first residual sd is the sd of `y - offset - own` over the rows in, the spread of what
  the forest is first given, or 1, today's default, where that is zero or not finite. Any fallback of E takes
  the "response" route whole: both ends, own zero, its sd.
- A fitted range of zero width is undefined (E); if the response's own is zero too, nothing is installed and
  dbarts's own warning for a constant response stands.
- `offset_type` other than "default" (a debugging aid): own keeps the parts the forest's offset holds in the
  sweeps (random alone under "fixef", fixed alone under "ranef", both under "bart", where the offset is not
  subtracted for the range); under "parametric" own is zero and no fit is run. The C++ creation is not changed:
  under these types the forest's first offset is the user's alone, as today.

E. When the range is the response's own though another was asked, and how the user is told.

| `range.fallback` | when | told by |
|---|---|---|
| "lme4" | random terms, "mixed", lme4 not available | a warning |
| "undefined" | n - p < 1; a random term with n effects or more; a fitted width of zero | a verbose line |
| "failed" | C.1 or C.3 raised an error or gave a non-finite value | a warning |
| "wider" | the fitted range is wider than the response's own, under "mixed" or "lm" | a verbose line |

1. Record, on the fit: `range.bart`, the pair (NULL for a binary fit); `range.route`, what was asked, one of
   "mixed", "lm", "response", "supplied"; `range.fallback`, present exactly in the table's cases;
   `range.info`, a list of n, p, the fit that ran ("lme4", "linear" or "none"), the added, dropped and left-out
   column names, the variance parameters and the optimizer's message.
2. Warnings: `lme4 is not installed, so the BART response range is the response's own, not the mixed model's;
   install lme4, or set 'bart_range'`; `the initial fit for the BART response range failed (<message>); the
   response's own range is used; see 'bart_range'`.
3. Verbose lines (`verbose` above 0), from the calling process and never from a chain's worker: before a fit,
   `estimating the BART response range by a linear <mixed >model (lme4) on <n> rows and <p> fixed columns; see
   'bart_range'`; on "undefined" or "wider", `the BART response range is the response's own, <lower> to <upper>:
   <the initial fit is not defined on <n> rows | the fitted range was wider>; see 'bart_range'`.
4. `summary()` prints `BART response range: <lower> to <upper> (<route>[; the response's own: <fallback>])`, and
   no line for a fit that has no `range.bart`.

F. Installing the range.
1. One helper takes a sampler holding no offset and a pair. It makes the vector of Context's first run from the
   sampler's own weights, calls `setResponse` with `updateScale = TRUE`, then `setResponse` with the sampler's
   own response and the range held, both with `updateState = FALSE`. It then reads the sampler's
   [`getLeafPrior`](../../man/dbartsSampler-class.Rd): `response.scale` must be identical to the pair's width
   and `response.shift` within 1e-12 of its middle, else it stops with an internal error. The vector is its own
   function.
2. Creation: each chain's worker installs the pair after it creates its sampler and before the C++ creation.
   The C++ creation passes the re-derive flag false; the sweep loop passes it false and loses its schedule. Those
   are the two places the range is derived today, and the only C++ edits.
3. Restore: the fit's `range.bart` for every chain. A fit saved before this change has none; its chains take the
   pair each state records, as today, through the same helper. Data carrying an offset is refused, as today.
4. No dbarts change. The engine holds an entry that sets the pair
   ([`Sampler::setAnchor`](../../src/bartcore/sampler.hpp), reached by a re-creation through
   [`applyAnchor`](../../R/dbarts.R)); it is not a documented method. See Calls made.

G. Binary fits. No range. The start is built from the parsed pieces, never by evaluating the call again.
- Random terms and lme4 available: the probit mixed model of Context's second run on [1, the model's own
  columns], the kept rows the 0/1 weights select, and the offset.
- Otherwise `glm.fit` with a probit link on the same columns, with the warning `lme4 is not installed, so the
  binary fit starts from a probit regression on the fixed effects alone; install lme4 to start from the mixed
  model` where the model has random terms. Where it does not converge or reports fitted probabilities of 0 or 1
  (separation), and where the mixed model fails (a warning), the start is zero.
- The starting offset is the own part of C.4, on the linear-predictor scale, where today it is the fitted
  probability.

H. lme4 as a suggested package.
1. The two operators go; every helper is the package's own copy. reformulas leaves Suggests, nothing else naming
   it; lme4 stays, for the fits, the tests and the help's links. No call names lme4 outside C.1, G and H.2.
2. One internal function answers whether lme4 is used: the option `stan4bart.lme4` is not FALSE and
   `requireNamespace("lme4", quietly = TRUE)`. The option is how a test runs the paths without lme4 on a machine
   that has it.

## Change

stan4bart, by file: a new R file for B to E and G; `stan4bart()` (the argument, B.1, the initial-fit block
replaced); the fit function and the chain worker (the pair handed down and installed); the restore and
`summary`; the two C++ flags; the lme4-derived source file (H); DESCRIPTION and NAMESPACE (the stats functions
newly called); tests; help; NEWS; its TODO (per-chain-leaf-prior-anchor closes). Here: this plan's Status and
Landing note, the design note's Status, TODO, the ledger entry, and the amendment of response-scale-rows named
below.

What a user sees:
- Every continuous fit's draws change at a given seed, and at 200 rows the posterior moves by the amounts in the
  design note's section 3.7. Chains share one prior. The estimated range is never wider than the response's own.
- `bart_range`; `range.bart`, `range.route`, `range.fallback`, `range.info`; the summary's line; the verbose
  lines; the warnings of A, E and G.
- A model with random terms fitted without lme4 warns, and its draws differ from the same seed's with lme4; they
  can also differ between versions of lme4. A model with no random term never uses lme4.
- `bart = list(...)` no longer reaches `bart_args` by partial matching; R reports it as ambiguous.
- Weights beside a row dropped for a missing value no longer stop the fit, and an offset beside one is aligned
  with its rows (such a fit's results change).
- A formula held in a variable fits, and the initial fit no longer runs a hidden BART fit.
- `extract(type = "trees")` reports leaf values in one internal unit for every chain: at a row, a tree's value
  times the range's width, summed over a draw's trees, plus the range's middle, is that draw's fit.
  `sampler.bart` is still one sampler per chain; `predict`, `fitted` and `extract` keep their shapes.
- A fit with rows of weight zero or an offset starts from corrected values; a binary fit starts on the
  linear-predictor scale. A fit saved before the change reloads as today and has no `range.bart`.

Constraints: no dbarts code; no new dependency; no size or time constant chooses a route; `mvbart` untouched.
Not in scope: one sampler for every chain, splines, a learned range, a dbarts method setting the range; and two
defects that predate the slice, for TODO: a reload of a fit with gp leaves replays 0.013 off its stored fits,
and `extract(type = "ranef")` stops ("subscript out of bounds") on a `||` term.

## Tests

Mechanisms: "without lme4" is `options(stan4bart.lme4 = FALSE)`, restored on exit; a failing fit is a function
that stops, passed as the route function's fitter argument. A test that calls lme4 itself sits inside
`requireNamespace("lme4", quietly = TRUE)`. Fits are 2 chains of 20 to 40 sweeps on 1 core unless said; the route
function is called directly on a `stan4bart(iter = 0)` object where no sampling is needed.

- Kept rows (B.1): a fixture with a missing forest value, a missing fixed value, a missing response, a zero
  weight and an offset, with and without `subset`: weights, offset, response and designs have one length; the
  offset and weights equal the kept rows' by position; the fit runs.
- Definition, with lme4: two random terms, one a correlated slope, weights with zeros, an offset, a factor and an
  exact copy of a model column inside `bart()`: `range.bart`, the starting offset and the coefficients equal,
  within 1e-6 of the width, those computed in the test from `lme4::lmer` on the written-out formula over the
  rows of positive weight. A grouping level with no row of positive weight starts at zero.
- Columns (rules 1 to 3), on the column builder: a logical, an ordered factor, a log term; a constant, a sum of
  two columns and a copy of the model's own column dropped, the model's own kept; an unordered factor as levels
  - 1 indicators; the grouping factor inside `bart()` and a factor finer than it left out and named; a coarser
  one kept. B.5: a forest holding only the model's own column prints its line, and equals the model's own fit.
- Shared and held: with a supplied pair far from the data's own, and under "mixed" with zero weights and an
  offset, every chain's state records a pair identical to `range.bart` after warm-up and every sampler's
  `getLeafPrior` reports it; once more at home on 2 workers.
- The helper: its vector holds the lower end exactly once, on the first row of positive weight, for weights zero
  on the even rows, on the odd rows, on all but two rows, and absent; installed, it lands exactly in each.
- Routes: "lm" equals `lm.wfit`'s range computed in the test within 1e-10 of the width; "response" and each
  fallback are `identical` to `range((y - offset)[weights > 0])`; a pair is `identical` to itself. Traced: lme4's
  builder is called once under "mixed" and never under the other three, nor for a model with no random term. The
  starting sd equals the sd of `y - offset - own` over the rows in, under each route.
- The table of E, each row asserting the pair, the two record fields, the warnings counted with
  `withCallingHandlers` and the lines captured under `verbose` 1 and 0: without lme4 (one warning naming lme4;
  none under "lm", "response", a pair, or a model with no random term); 30 rows and 40 forest columns; one
  random effect per row; a failing fitter (one warning); a forest column 0.99995 correlated with the model's own
  ("wider", no warning, the response's own pair, start zero); a singular fit keeps "mixed" with no record.
- The argument: "resp" is "response"; "splines", a function, length 3, NA, the upper end first and
  `c(-1e308, 1e308)` stop by their words, also at `iter = 0`; a pair of 1000 to 1001 on data near 50 warns once
  and is used; a binary response with "lm" warns once and has a NULL `range.bart`; `bart = list()` stops.
- `offset_type`: under each of the five the route function's pair equals the D formula from the test's own lme4
  parts; "parametric" calls no fitter.
- Binary: with lme4 the starting offset equals `glmer`'s linear predictor less offset and intercept terms within
  1e-4; without it, `glm.fit`'s, with one warning; a separated fixture starts at zero; a formula held in a
  variable fits, for both families, and no BART output is printed by the initial fit.
- Restore: a fit with zero weights on every second row, saved and read back, predicts its stored fits within
  1e-10; so does a fit rebuilt without `range.bart` from two chains whose states record different pairs;
  `summary` of it prints no range line.
- lme4 as suggested: no function of the namespace has lme4's or reformulas's namespace as its environment.

Edits: the serialization tests' precondition that no two chains share a range becomes that all do, and its
leaf-sum identity uses the one pair; the weights-and-offset tests gain a zero weight; the at-home accuracy
thresholds are rerun over five seeds and reported, and a threshold that fails is a stop, not a moved seed.

Reviewer's mutants, each of which must fail a test: weights or offset left uncut; the added columns left out; a
factor's guard removed, or turned on a coarser factor; indicators for an ordered factor; the model's own column
dropped for its copy; lme4 handed the unpruned columns; the offset subtracted twice, or not at all; zero-weight
rows left in the fit, or in the range; the own part holding the intercept, or the added columns' part, or
uncentered; the starting sd taken from the fit; either C++ flag back to true; the helper's vector alternating,
or its lower end on row 1 whatever the weight; the install skipped; "lm" calling lme4; "response" reading every
row; the bound removed, applied to a pair, or replacing one end; a fallback keeping the fit's start; no warning
without lme4, or one for a model with no random term; a warning on "undefined"; none on "failed"; a fallback on
a singular fit; the n-effects check removed; a pair sorted; an infinite width accepted; the argument unchecked at
`iter = 0`; the binary start on the probability scale, or from the call; an operator left on one helper.

## Verification

In stan4bart, on the slice's own library `<lib>` holding the dbarts under test and lme4, at most 2 cores:

    R CMD INSTALL --preclean -l <lib> .
    NOT_CRAN=true R_LIBS=<lib> Rscript -e 'tinytest::test_package("stan4bart")'
    R_LIBS=<lib-without-lme4> Rscript -e 'tinytest::test_package("stan4bart")'
    R CMD build <clean copy>; R_LIBS=<lib> R CMD check --as-cran stan4bart_*.tar.gz
    _R_CHECK_DEPENDS_ONLY_=true R_LIBS=<lib> R CMD check --as-cran stan4bart_*.tar.gz

Each ends with no failure, error or warning of the check's own; the last two differ only in the tests that sit
behind lme4, and the help's example, whose model has random terms, runs in the last with the lme4 warning. Then
the two gate scripts as stan4bart's benchmarks README gives them, and in dbarts
`Rscript tools/check-doc-freshness.R .`.

## Gates and baselines

- stan4bart's CI: check-standard, sanitizers, and gates. Its exactness gate must pass unchanged. Its posterior
  gate compares five recorded tiers: the binary tier must pass against the recorded baseline (its start changed,
  its posterior did not); the four continuous tiers are reported against it, not gated; then all five are
  re-recorded and the MANIFEST row says why.
- Consumers, the three CRAN packages that suggest stan4bart. bartCause: its test suite on its dbarts-1.0 branch
  with the built stan4bart and lme4 installed (its `group.by` fits build `y ~ (1 | g) + bart(...)`); a seeded
  pin that moves is reported with old and new value. The bartCause leg of dbarts's
  [revdep-smoke.yaml](../../.github/workflows/revdep-smoke.yaml) runs without stan4bart and is not this gate.
  tidytreatment and WeightIt: their CRAN sources searched for calls of `stan4bart(` and of `bart =`, then
  `R CMD check` of each tarball with the built stan4bart installed; a failure is reported, not repaired here.

## After the build

A rerun of the design's deciding cells on the built package (lme4 installed) against the branch before the
slice, same data and seeds, 4 chains of 1000 + 1000, 2 cores: many small groups at 200 rows (20 seeds) and 1000
(10), the control at 200 (20), the group-level factor at 1000 (10), weights and an offset at 1000 (10). About
140 fits. Rule, fixed here; a miss on any is a stop for the orchestrator, with the table:
1. In every fit the built range equals, within 1e-6 of a width, the range computed in the script from
   `lme4::lmer` on the written-out formula, and dec-B447's bound takes over in none.
2. At 1000 rows no accuracy quantity is worse under the design's rule (interval excluding zero on the bad side
   and a mean beyond 3 percent).
3. Many small groups at 200 rows: expected-value error between +0.4 and +5.0 percent and split error between
   -14.3 and -0.3 percent against today's (the design's intervals widened by their own half-width).
4. The control at 200 rows: neither quantity worse.

Reported beside it and gated on nothing, the cost of running without lme4: a random intercept on groups of 4 and
of 25 at 500 rows, binary and continuous, 10 seeds, with lme4, without it, and with it at another seed (120
fits). With 10 seeds a difference of about 1 percent in the fitted mean's error shows (the half-widths were 0.6
and 1.6 percent at 5). The help's sentence on fitting without lme4 quotes what it finds.

## Help and docs

- The `bart_range` item: what the range is (the prior on the forest's function has its mean at the middle and,
  at `k = 2`, a standard deviation of a quarter of the width); the default and its definition; that an estimated
  range is never wider than the response's own, and that a user who wants another gives a pair, for example
  `range(y)`; "lm" for large crossed grouping factors, where the mixed model can outlast the sampler, at a cost
  in small samples; that too narrow a range leaves group-level signal in the random effects and too wide a one
  costs accuracy and mixing in small samples; that a strongly curved group-level effect is a case for a pair.
  No list of cases in which the estimate runs wide (dec-B447), and no timing.
- lme4, in Details: it fits the initial mixed model when installed; without it a model with random terms takes
  the response's own range (a binary fit a simpler start) with a warning, and results differ from those with
  it. The `verbose` item names the two lines. Value: the four elements of E.1. The generics' help: the unit of
  extracted leaf values.
- NEWS, under upgrading: the range and its argument, with "results differ from earlier versions, including
  under a fixed seed"; `bart =` no longer abbreviates `bart_args`; the kept-rows fix; the binary start's scale;
  that a package built where lme4 is installed now runs where it is not.

## Steps

Each ends with the suite green on the slice's own library. Steps 1 and 2 change no draw where lme4 is installed
and may land first, together; steps 3 to 5 land together (a range computed and not installed, or installed and
not documented, is no state to leave the branch in).
1. Kept rows (B.1) with its fixture. ~45 lines.
2. The operators out, reformulas out of Suggests, H.2, the namespace test. ~95 lines, 60 of them mechanical.
3. The range: A to F and the continuous start, with their tests and the edits. ~560 lines.
4. The binary start (G) with its tests. ~120 lines.
5. Help, NEWS, stan4bart's TODO. ~80 lines.
6. After review: the confirming measurement and the table beside it, the posterior baseline's re-record, the
   consumer runs. ~110 lines.
7. Records here. ~55 lines.

Stop when: the diff passes ~2100 lines; the helper's read-back fails on any fixture; lme4's functions refuse the
pieces on the lme4 under test; the exactness gate or the binary posterior tier moves; a test fails that Tests
does not name; the change needs dbarts code; a rule of After the build is missed.

## Against response-scale-rows

Neither slice needs the other. This one stops calling the re-deriving entries during a run and overwrites
creation's range, so afterwards dbarts's choice of rows reaches a continuous stan4bart fit only inside the
helper, whose vector is exact whether dbarts reads every row or the rows of positive weight (run for the first;
by construction for the second, both ends lying on rows of positive weight whenever two exist; with one such row
that plan's rule 4 reads every row).
- This first (recommended: it is ruled and smaller): response-scale-rows rebases the stan4bart paragraph of its
  [What a consumer sees](response-scale-rows.md#what-a-consumer-sees). stan4bart no longer calls
  `setOffset` with the flag true, a gaussian stan4bart fit with a zero weight changes with this slice and not
  with that one, and its verification pair (one seeded fit with zero weights differing, one without identical)
  becomes both identical. Its rule 3 stays for other callers.
- That first: nothing here changes but Context's sentence on the alternating vector, which then fails in fact.
- Its slice C (dec-B364: a re-derivation whose rows in hold fewer than two distinct values keeps the range and
  warns) would leave the install undone on a fit with one row of positive weight; the read-back turns that into
  a stop, and slice C's plan must name the helper.

## Calls made

The orchestrator's, written in by its direction: the fitted width is accepted as it comes up to dec-B447's
bound, with no factor toward today's method; the bound holds for every estimated route, "mixed" and "lm"; rows of
weight zero leave the initial fit and the range now, not with response-scale-rows; the restore of one sampler
per chain stays; the design's rules 1 to 11 stand, but for the starting sd below; a fallback inherent in the data
is a verbose line, one a user can repair (lme4 missing) or a failure is a warning; names are matched with
`match.arg`; a pair that misses the data wholly warns and any other pair is the user's number; `bart =` stops
abbreviating `bart_args`.

The planner's:
- The starting sd is the spread of what the forest is first given, not the fit's residual sd (rule 4's second
  half): Context's last run, and it is what today's fit supplies in effect, its residual holding the forest's
  signal. Rejected: the fit's sd bounded by a constant; the fit's sd kept where the bound does not take over.
- A fallback takes the "response" route whole. Rejected: the response's ends with the fit's start, whose own part
  is what made the range wide.
- The install is the made-up response, by documented methods, with a read-back. Deferred: a documented dbarts
  method over the engine's existing entry, about 60 lines there. Rejected: calling the undocumented bridge entry.
- lme4's fit is read from the criterion's environment, as its own constructor reads it; that constructor refuses
  our frame. An error there, on any version of lme4, is the "failed" fallback with its warning, and the
  Definition test shows it on CRAN's checks.
- The check for a random term with n effects is ours, in lme4's words' meaning; lme4's builder skips it.
- The option of H.2 is internal and not in the help. Offered to the orchestrator: documenting it would let a user
  reproduce on one machine the fit of a machine without lme4.
- The binary guard is `glm.fit`'s own report of non-convergence or of fitted probabilities of 0 or 1.
- A pair fits no mixed model; its start is the linear fit's. `range.route` is what was asked; a model with no
  random term records "mixed" with `range.info`'s fit "linear", so `range.fallback` means one thing.
- Under an `offset_type` other than the default the range follows what the forest's offset holds and the C++
  creation is left alone. Rejected: a C++ edit to start those types from the fit, with no reader to test it by.
- reformulas leaves Suggests with its operator: the same defect, and identical parses on six formulas (run).
- The critique's items on the withdrawn routine (its optimizer, its tolerance, its boundary stop, one random
  effect per row as its problem) are moot; the last became C.2.

## Not verified

Nothing was built in stan4bart. lme4's functions were run on 2.1.0 only, and not on a model with more than two
random terms beyond the timing cases. Not run: the route function as a whole, the bound inside a sampler, the
starting-sd rule in a fit, the `offset_type` rule, the zero-width rule, the old-fit restore, the binary start on
the linear-predictor scale, `glm.fit`'s guard on a separated fixture, the helper on a build with
response-scale-rows, the cost of a many-level forest factor's indicator columns, a check with
`_R_CHECK_DEPENDS_ONLY_`, bartCause's suite, tidytreatment's and WeightIt's use of stan4bart (not read),
stan4bart's posterior gate, Windows, forked workers.
