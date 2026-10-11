# stan4bart-creation-mapping: stan4bart fixes its response range at creation, from a mixed model it fits itself

Status: PLANNED (dec-B437, dec-B441 to dec-B445). Two open calls. The code is stan4bart's; this file is its record.

agent: a blind critique of this plan; one opus implementer (new numerics in R, ten lines of C++); one opus
reviewer who runs the mutants.
rng: POSTERIOR-CHANGING for every continuous stan4bart fit (the leaf prior's center and spread change, and one
range replaces one per chain). SHIFTING for a binary fit (its starting values). NEUTRAL for dbarts: no dbarts
code changes, so no dbarts baseline, snapshot or exact gate moves.
window: before the merge to main (dec-B437). Either side of [response-scale-rows.md](response-scale-rows.md);
see [Against response-scale-rows](#against-response-scale-rows).
budget: ~1150 lines, almost all in stan4bart (R ~390, C++ ~10, tinytest ~500, help and NEWS ~70, benchmarks
~125; records here ~55). Stop at 1.5 times.

## Goal

A continuous stan4bart fit has one response range for its forest, computed once before sampling, the same for
every chain and never re-derived: by default from a linear mixed model that stan4bart fits itself and that holds
the forest's columns as linear terms; else as `bart_range` says. The fit records the pair and how it was
obtained. lme4 has no use at run time. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Why, with the measurements: [stan4bart-response-range.md](../design/stan4bart-response-range.md). Ruled, each
covering what it names: dec-B437 (one range, fixed at creation, shared, never re-derived), dec-B441 (the forest's
columns enter linearly), dec-B442 (the mixed model is fitted inside stan4bart over Matrix, and matches lme4 on
large crossed models in the tests before it replaces lme4), dec-B443 (no run-time use of lme4; a binary fit
starts from glm with the grouping terms as fixed effects, once a check shows binary results after warm-up do not
move), dec-B444 (where the fit cannot be computed the range is the response's own, less any offset), dec-B445
(`bart_range`: a pair of numbers, "mixed", "lm" or "response"). dec-B419's restore of one sampler per chain
stays, with dec-A197's calls.

"Run" is this plan's: stan4bart's bartcore branch at 3d2a295 against a build of dbarts bartcore 1c64c477, arm64
macOS, R 4.6.1, R's own BLAS, lme4 2.1.0, Matrix 1.7.5, on a shared machine.

Today, read in stan4bart's source at that commit:
- `stan4bart()` evaluates `lme4::lmer` (binary: `glmer`) on the formula without its `bart()` term; where lme4 is
  not installed or the fit stops, `lm` (`glm`) with the grouping terms as fixed effects; then the fixed terms
  alone. The fitted values are the forest's starting offset, the fit's sigma the starting residual sd.
- The range is derived twice: when each chain's sampler is created, through
  [`dbarts_sampler_setOffset`](../../inst/include/dbarts/dbarts.h) with its re-derive flag true (the response
  less the starting offset); and through warm-up, the same entry with the flag true on a thinning schedule (every
  sweep of the first eighth, every second of the next, and so on). Nothing else derives it.
- The restore builds one sampler per chain and puts each on its chain's recorded pair by
  [`dbartsSampler$setResponse`](../../man/dbartsSampler-class.Rd) with `updateScale = TRUE` on a vector
  alternating the pair, then the real response with the range held.
- At the top of its lme4-derived source file two operators replace stan4bart's copies of lme4's and reformulas's
  helpers with those packages' own.

Run, and used below:
- The install. On the tip, a sampler built as stan4bart builds it takes a pair exactly from a made-up response
  (the lower end on the first row of positive weight, the upper on every other row): the model's recorded pair is
  identical to the pair with and without zero weights on every second row, and stays so through `setOffset` with
  the range held. A sampler created on other data and then installed gives bit-identical draws over 30 sweeps to
  one whose creation derived the same pair, under the default leaf prior, a drawn k and a named sd. A state
  forced into an installed sampler predicts identically. Today's alternating vector lands too, because dbarts
  reads every row; on the rows of positive weight alone it would hold one value where zero weights fall on one
  parity.
- The routine. The design's prototype against lme4 on fifteen models, three of the design's four timing cases
  among them: criterion equal within 5e-10, range ends within 2e-5 of the width, time 0.28 s for 0.36 (25,000 nested
  groups at 1e5 rows), 39 s for 46 (5000 by 2000 crossed, 1e5 rows), 76 s for 96 (3000 by 3000, 3e4 rows). At
  three crossed factors with slopes and 1e5 rows one evaluation takes 6.0 s for lme4's 9.5 and neither ends in
  ten minutes. On a second draw of seven models one singular fit stopped on the boundary 0.017 above lme4's
  criterion (ends off by 2.7e-4 of the width). The numbers are in the design note's section 6.
- Three facts against the record: [Found against the record](#found-against-the-record).

## Design

A. The argument. `bart_range = "mixed"`, last in `stan4bart()`'s formals.
- A string is exactly "mixed", "lm" or "response", in lower case; NULL is "mixed". Anything else that is not a
  pair stops: `'bart_range' must be c(lower, upper), "mixed", "lm" or "response"`.
- A pair is two finite numbers, the second above the first; else
  `'bart_range' must be two finite numbers, the lower first`. It is never reordered.
- For a binary response a value other than "mixed" warns `'bart_range' is not used for a binary response`.

B. Rows and columns, from what `stan4bart()` has parsed (the response, the model's fixed columns, the
random-effect structure, the forest's data object, weights, offset).
1. Rows in: the rows the model keeps, of positive weight. n is their count. The fit and the range read no other
   row (the orchestrator's call; the design's rule 6).
2. Added columns, from the forest's training matrix and its column kinds. A numeric or ordered column enters as
   the value the forest splits on (a logical 0/1, an ordered factor its level score, a transformed term as
   transformed). An unordered factor enters as indicators of its levels present among the rows in, less the
   first; unless every one of those levels lies within one level of some grouping factor of the model, when the
   factor is left out and its name recorded (rules 1 and 3).
3. Aliasing. The fixed block is [1, the model's own columns, the added columns] on the rows in. A pivoted QR at
   `lm`'s tolerance (1e-7) keeps the first of any dependent set, so the model's own columns are kept before the
   added ones; the dropped names are recorded (rule 2). p counts the kept columns.

C. The routine, one function over Matrix and stats.
1. Inputs: z, the response less the offset on the rows in; X, the kept fixed columns; the random-effect
   structure's transposed design `Zt` cut to the rows in, its factor template `Lambdat` with the index `Lind`
   that places the variance parameters theta in it, their starting values and lower bounds; the weights w.
2. With h = sqrt(w), Zw = Zt diag(h), Xw = h X and zw = h z: the factor L of `Lambda' Zw Zw' Lambda + I` is a
   sparse Cholesky whose fill-reducing ordering is found once and whose values are updated at each theta
   (`Matrix::Cholesky`, `Matrix::update`). Then, in Bates, Maechler, Bolker and Walker's notation (2015, section
   3): cu = L^-1 P Lambda' Zw zw, RZX = L^-1 P Lambda' Zw Xw, RX the Cholesky factor of Xw'Xw - RZX'RZX, beta
   from RX'RX beta = Xw'zw - RZX'cu, u from L'P'u = cu - RZX beta, b = Lambda u, and the penalized residual sum
   r2 = sum(w (z - X beta - Zt'b)^2) + u'u. The criterion is
   `2 log|L| + 2 log|RX| + (n - p) (1 + log(2 pi r2 / (n - p)))`. The crossproducts of X and z are formed once.
3. Optimizer: `stats::nlminb` from the structure's starting values with its lower bounds (0 on a factor's
   diagonal), relative tolerance 1e-8, iteration and evaluation limits 10,000. It is in base R, takes bounds and
   needs no gradient; in the runs it took 0.4 to 3.3 times lme4's evaluations. A trial theta at which RX cannot
   be formed returns Inf, a rejected step.
4. A singular fit is an optimum on a bound: that term's column of Lambda is zero, L stays positive definite, its
   predicted effects are zero, and the result is used as any other. So is a fit that ends at a limit or with a
   nonzero convergence code (rule 10); the code and the evaluation count are recorded.
5. Returns the kept columns' coefficients, b, the residual sd `sqrt(r2 / (n - p))`, theta, and the own part on
   EVERY row: the model's own columns centered at the means the parametric sampler centers them at, times their
   coefficients, plus Zt'b (zero for a level with no row in). The intercept and the added columns' part are left
   out.
6. A failure is an error at the optimum (RX not positive definite) or a non-finite criterion, coefficient or
   own part there.
7. Each evaluation returns to R, so an interrupt stops the fit between evaluations.

D. Routes. In each the range is the smallest and largest of `y - offset - own` over the rows in.
- "mixed": own and the starting residual sd from C. A model with no random term has no C; it takes "lm" and
  records "lm".
- "lm": weighted least squares of z on the same kept columns, no random term; own is the model's centered
  columns times their coefficients; the sd is that fit's on n - p degrees of freedom.
- "response": own is zero, the range is that of `y - offset` over the rows in, the starting sd their sd. Nothing
  is fitted.
- A pair: the range is the pair. The starting offset and sd are "lm"'s (or "response"'s where "lm" is not
  defined); no mixed model is fitted. Recorded as "supplied".
- The forest's starting offset is own; the user's offset is added once, where the sweep loop adds it (rule 7).
  Under an `offset_type` that swaps part of the model for the offset, own keeps only the parts the forest's
  offset will hold in the sweeps (random part alone under "fixef", fixed alone under "ranef", none under
  "parametric", both under "bart", where the user's offset is not subtracted for the range either); the C++
  creation, which adds the starting offset under the default type only, follows.
- A range of zero width (the rows in hold one value) is not installed; creation's own stands, as today.

E. Where the fit cannot be computed (dec-B444).
1. Undefined: n - p < 1 under "mixed" or "lm". Failed: C.6. Either way the route is "response".
2. Record, on the fit: `range.bart`, the pair (NULL for a binary fit); `range.route`, one of "mixed", "lm",
   "response", "supplied"; `range.fallback`, "undefined" or "failed", present only where the route is not the
   one asked for; `range.info`, a list of n, p, the added, dropped and left-out column names, theta, the
   convergence code and the evaluation count.
3. Verbose (`verbose` above 0), from the calling process and never from a chain's worker: before a fit,
   `estimating the BART response range by a linear mixed model on <n> rows, <p> fixed columns and <q> random
   effects; see 'bart_range'` ("by a linear model" under "lm"); on a fallback, `the BART response range is the
   response's own, <lower> to <upper>: <the fit has no residual degrees of freedom in <n> rows | the fit
   failed>; see 'bart_range'`.
4. A warning only for "failed": `the initial fit for the BART response range failed (<message>); the response's
   own range is used; see 'bart_range'`. None for "undefined", a singular fit or a convergence code.
5. `summary()` prints one line: `BART response range: <lower> to <upper> (<route><, fallback>)`.

F. Installing the range.
1. One helper takes a sampler holding no offset and a pair. It makes the vector of Context's first run from the
   sampler's own weights, calls `setResponse` with `updateScale = TRUE`, then `setResponse` with the sampler's
   own response and the range held, both with `updateState = FALSE`. It then reads the sampler's
   [`getLeafPrior`](../../man/dbartsSampler-class.Rd): `response.scale` must be identical to the pair's width
   and `response.shift` within 1e-12 of its middle, else it stops with an internal error. The vector is its own
   function, so its property is tested without a sampler.
2. Creation: each chain's worker installs the pair after it creates its sampler (and after a binary mask, which
   takes none) and before the C++ creation. The C++ creation passes the re-derive flag false; the sweep loop
   passes it false and loses its schedule. Those are the two places the range is derived today, and the only
   C++ edits besides D's `offset_type` line.
3. Restore: the restore takes the fit's `range.bart` for every chain. A fit saved before this change has none;
   its chains take the pair each state records, as today, through the same helper. Data carrying an offset is
   refused, as today.
4. No dbarts change. A dbarts entry that sets the pair directly exists inside the engine
   ([`Sampler::setAnchor`](../../src/bartcore/sampler.hpp), reached by a re-creation through
   [`applyAnchor`](../../R/dbarts.R)) and is not a documented method: see Calls made.

G. Binary fits (dec-B443). No range. The start is `glm` with a probit link on the formula with its grouping
terms as fixed effects, as today without lme4; on an error, on the fixed terms alone; on another, zero. The
starting offset is the linear predictor less the user's offset, centered over the rows in, where today it is the
fitted probability. See Open call 2 before this half is built.

H. lme4 at run time. The calls to `lme4::lmer`, `glmer` and `lmerControl` go. The two operators go and every
helper is the package's own copy. DESCRIPTION: lme4 stays under Suggests, for the tests and the help's links;
reformulas leaves it, nothing else naming it.

## Change

stan4bart, by file: a new R file for B to E; `stan4bart()` (the argument, the initial-fit block replaced, G);
the fit function and the chain worker (the pair handed down and installed); the restore and `summary`; the C++
creation and sweep loop (F.2, D); the lme4-derived source file (H); DESCRIPTION; tests; help; NEWS; its TODO
(per-chain-leaf-prior-anchor closes). Here: this plan's Status and Landing note, the design note's Status, TODO,
the ledger entry, and the amendment of response-scale-rows named below.

What a user sees:
- Every continuous fit's draws change at a given seed, and at 200 rows the posterior moves by the amounts in the
  design note's section 3.7. Chains share one prior.
- `bart_range`; `range.bart`, `range.route`, `range.fallback`, `range.info`; the summary's line; the verbose
  lines; one new warning.
- `predict`, `fitted` and `extract` return what they did in shape. `extract(type = "trees")` reports leaf values
  in one internal unit for every chain: at a row, a tree's value times the range's width, summed over a draw's
  trees, plus the range's middle, is that draw's fit. `sampler.bart` is still one sampler per chain.
- A fit with rows of weight zero, or with an offset, starts from corrected values and takes its range from the
  rows of positive weight.
- A binary fit's draws change at a given seed (its start). A seed gives the same draws with lme4 installed or
  not, and a build made beside lme4 runs without it.
- A fit saved before the change reloads as it does today, each chain at its own range; it has no `range.bart`.

Constraints: no dbarts code; no new dependency; no size or time constant chooses a route; `mvbart` untouched.
Out of scope: one sampler for every chain, splines, a learned range, a drawn k, a dbarts method setting the range.

## Tests

New file for the range. Fixtures of 150 to 400 rows; fits of 2 chains and 20 to 40 sweeps on 1 core unless said.
The route function is called directly on a `stan4bart(iter = 0)` object where no sampling is needed.
- Definition: two grouping terms, a fixed column, weights with zeros, an offset and a factor inside `bart()`.
  `range.bart` equals, within 1e-4 of its width, the range over the rows of positive weight of the response less
  the offset less the own part of `lme4::lmer` fitted to those rows with the forest's columns added (skipped
  without lme4); the starting offset and sd agree likewise.
- Columns (rules 1 to 3), on the column builder: a logical, an ordered factor, a log term; a constant, a sum of
  two columns and a copy of the model's own fixed column dropped, the model's own kept; an unordered factor as
  levels - 1 indicators; the grouping factor inside `bart()` and a factor finer than it left out and named; a
  factor coarser than it kept.
- Shared and held: after a fit with warm-up and `keepTrees`, every chain's state records a pair identical to
  `range.bart`, and every sampler's `getLeafPrior` reports it.
- The helper: its vector holds the lower end exactly once, on the first row of positive weight, for weights zero
  on the even rows, on the odd rows, on all but two rows, and absent; installed on a sampler it lands exactly in
  each case; a pair is installed on a sampler created at another range and the same state predicts identically.
- Routes: "lm" equals `lm.wfit`'s range computed in the test within 1e-10 of the width; "response" is
  `identical` to `range((y - offset)[weights > 0])`; a pair is `identical` to itself, route "supplied". Traced:
  the routine is called once under "mixed" and never under the other three. A model with no random term
  records "lm".
- dec-B444: 30 rows and 40 forest columns give route "response", fallback "undefined", zero warnings (counted
  with `withCallingHandlers`), the fallback line once under `verbose = 1` and nothing under 0. With the routine
  stubbed to stop: fallback "failed" and exactly one warning. A singular fit (no group effect) keeps route
  "mixed" with no warning.
- The argument: each refusal by its words (length 3, NA, the upper end first, "MIXED", "splines", a function);
  NULL is "mixed"; a binary response with "lm" warns once and has a NULL `range.bart`.
- Restore: a fit with zero weights on every second row, saved and read back, predicts its stored fits within
  1e-10; so does a fit object rebuilt without `range.bart` from two chains whose states record different pairs.
- lme4: no function of the namespace has lme4's or reformulas's namespace as its environment; a seeded fit is
  `identical` before and after `loadNamespace("lme4")`.

New file for the routine, skipped without lme4; lme4's criterion is `mkLmerDevfun`'s on the same pieces.
- Criterion: at lme4's optimum and two other settings the routine's criterion equals lme4's plus the sum of the
  log weights within 1e-8 relative, and its coefficients, b and sd at lme4's theta agree within 1e-8.
- Optimum: range ends within 1e-3 of the width, sd within 1e-3 relative, criterion no more than 0.05 above
  lme4's (the boundary stop of Context was 0.017 and 2.7e-4 of the width).
- Models, on CRAN: a correlated intercept and slope; three correlated effects on 15 levels and on 6; a slope
  with a crossed intercept; `||`; no group effect; a slope proportional to its intercept; nested terms with
  weights and an offset; crossed 500 by 200 at 1e4 rows; 25,000 nested groups at 1e5 rows. At home: crossed 1000
  by 1000 at 1e4 rows, three crossed factors with slopes at 1e4 rows, crossed 5000 by 2000 at 1e5 rows.

Edits: the serialization tests' precondition that no two chains share a range becomes that all do, and its
leaf-sum identity uses the one pair; the weights-and-offset tests gain a zero weight; the at-home accuracy
thresholds of the continuous and binary tests are rerun over five seeds (a seed moves if a threshold holds on
four of five; else stop).

Reviewer's mutants, each of which must fail a test: the added columns left out; a factor's guard removed, or
turned on a coarser factor; indicators for an ordered factor; the model's own column dropped for its copy; the
offset subtracted twice, or not at all; zero-weight rows left in the fit, or in the range; the own part holding
the intercept, or the added columns' part, or uncentered; either C++ flag back to true; the helper's vector
alternating, or its lower end on row 1 whatever the weight; the read-back removed with the install skipped; "lm"
running the mixed model; "response" reading every row; a pair sorted; an unknown string read as "mixed"; a
warning on "undefined", none on "failed", a fallback on a singular fit, the fallback taken from the model's own
fit; an old fit restored at creation's range; an operator left on one helper; in the routine, the weights' root
dropped from Zw, log|RX| dropped, n for n - p, the penalty u'u dropped, the permutation dropped, Lambda not
refreshed.

## Gates and baselines

- Before anything else replaces lme4: the routine's test file green, its at-home models included (dec-B442).
- stan4bart's tinytest at home (`NOT_CRAN`) and `R CMD check --as-cran` from a tarball built from a clean copy,
  each on a library holding lme4 and on one that does not, the package installed there; at most 2 cores.
- stan4bart's CI: check-standard, sanitizers, and gates. Its exactness gate must pass unchanged. Its posterior
  gate compares five recorded tiers; the four continuous tiers are reported against the recorded baseline, not
  gated (the posterior changed), then all five are re-recorded and the MANIFEST row says why. The binary tier
  must pass against the old baseline before the re-record: that is part of dec-B443's check.
- dec-B443's check, rule fixed here: the documented binary model and a random intercept on groups of 4 and of
  25, at 500 and 2000 rows, 5 seeds, 4 chains of 1000 + 1000; the branch before the slice with lme4 installed,
  the built package, and the former at another seed. The new start passes if, in every cell, its paired
  difference from today's in the error of the fitted probability and in split error has an interval containing
  zero or a size no larger than the rerun's.
- Consumers: bartCause's test suite on its dbarts-1.0 branch against the result (its `group.by` fits build `y ~
  (1 | g) + bart(...)`, and it reads only whether `sampler.bart` is NULL); a seeded pin that moves is reported
  with old and new value. dbarts's [revdep-smoke.yaml](../../.github/workflows/revdep-smoke.yaml) dispatched
  once on the result. treatSens and bairrtt do not use stan4bart.
- dbarts: `Rscript tools/check-doc-freshness.R .` for the records; nothing else here moves.

## After the build

A rerun of the design's deciding cells on the built package against the branch before the slice, same data and
seeds, 4 chains of 1000 + 1000, 2 cores: many small groups at 200 rows (20 seeds) and 1000 (10), the control at
200 (20), the group-level factor at 1000 (10), weights and an offset at 1000 (10). About 140 fits. Rule, fixed
here:
1. In every fit the built range equals, within 1e-3 of a width, the range computed in the script from
   `lme4::lmer` with the forest's columns added (the Definition test's oracle).
2. At 1000 rows no accuracy quantity is worse under the design's rule (interval excluding zero on the bad side
   and a mean beyond 3 percent).
3. Many small groups at 200 rows: expected-value error between +0.4 and +5.0 percent and split error between
   -14.3 and -0.3 percent against today's (the design's intervals widened by their own half-width).
4. The control at 200 rows: neither quantity worse.
A miss on any is a stop for the orchestrator, with the table.

## Help and docs

- The `bart_range` item of stan4bart's help: what the range is (the prior on the forest's function has its mean
  at the middle and, at `k = 2`, a standard deviation of a quarter of the width); the default and its definition;
  the three other values and when each serves ("lm" for large crossed grouping factors, where the mixed model
  can outlast the sampler, at a cost in small samples; a pair, for example `range(y)`, where the range is known
  or a strongly curved group-level effect is expected); that too narrow a range leaves group-level signal in
  the random effects and too wide a one costs accuracy and mixing in small samples; the two cases of the design
  note's [7. What does not work](../design/stan4bart-response-range.md#7-what-does-not-work) (as many added
  columns as rows; as many group-level covariates as groups) as cases for a pair; that the fit can be
  interrupted between evaluations and a verbose fit says what it will compute. No timing is quoted.
- Value: the four elements of E.2; `sampler.bart` unchanged. The generics' help: the unit of extracted leaf
  values, as in Change. Details: the initial fit is the package's own; lme4 is not used.
- NEWS, under upgrading: the range and its argument, with "results differ from earlier versions, including
  under a fixed seed"; that zero-weight rows and an offset no longer enter the range twice or at all; that lme4
  is no longer used when fitting and a binary fit's start changed.

## Steps

Each ends with the suite green on the slice's own library.
1. The routine (C) and its test file, lme4 untouched. ~230 lines. STOP for the orchestrator if any model misses
   a tolerance (dec-B442's condition); its time beside lme4's on the at-home models is reported.
2. Rows, columns and routes (B, D, E), the argument (A), the record and the lines, with their tests; the fit
   still samples as today. ~330 lines.
3. The install (F): the helper, the worker, the two C++ flags, the restore, the summary; the shared-and-held,
   helper and restore tests and the edits. ~210 lines.
4. Binary start (G) and lme4 out (H), with dec-B443's check. ~150 lines. STOP if the check fails, or before
   building G at all if Open call 2 is unanswered.
5. Help, NEWS, stan4bart's TODO. ~70 lines.
6. After review: the confirming measurement, the posterior baseline's re-record, the consumer runs. ~125 lines.
7. Records here. ~55 lines.

Stop when: the diff passes ~1700 lines; a tolerance of the routine's file is missed; the helper's read-back
fails on any fixture; the exactness gate or the binary posterior tier moves; a test fails that Tests does not
name; the change needs dbarts code; a rule of After the build is missed.

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

The orchestrator's, written in by its direction: the enriched fit's width is accepted as it comes, with no
factor toward today's method; rows of weight zero leave the initial fit and the range now, not with
response-scale-rows; the restore of one sampler per chain stays, and only how a sampler is put on the range
changes; the design's rules 1 to 11 stand; the binary starting-offset fault is fixed here; a fallback is
announced under verbose and recorded, with a warning only for a numerical failure.

The planner's:
- The install is the made-up response, by documented methods, with a read-back. Deferred: a documented dbarts
  method setting the pair through the engine's existing entry, about 60 lines of R, help and tests there and a
  new public method before the release candidate, with no change to the flat C header unless a C caller wants it
  (then a minor version and a new hash). Rejected: calling the undocumented bridge entry from stan4bart.
- nlminb at 1e-8. Rejected: `optim`'s L-BFGS-B (run on seven models: no better criterion, 6 to 416 evaluations
  for nlminb's 6 to 285, and it refuses a non-finite value); nlminb's default tolerance (criterion within 1e-7
  of lme4's for 4.5e-5, at 6 to 340 evaluations; at 1e-8 the range is already within 2e-5 of a width of lme4's).
- A boundary stop at a singular fit is accepted (1 of 22 comparisons, 2.7e-4 of a width): rule 10.
- A pair fits no mixed model; its start is the linear fit's. Rejected: the mixed fit for the start, which the
  user supplying a pair to avoid it would still pay for; a zero start.
- A model with no random term records "lm", what ran; "undefined" is no residual degree of freedom, "lm" too;
  a zero-width range is left to dbarts's handling of a constant response, as today.
- The record is four elements beside `sampler.bart`, not attributes of the pair. `range.bart` is the name the
  help carried until dec-A197 removed it.
- A binary response with a named range warns and fits (rejected: stopping, for a wrapper that passes one value
  to both of its fits; silence).
- reformulas leaves Suggests with its operator: the same defect, and identical parses on six formulas (run).
- Under an `offset_type` other than the default the start and the range follow what the forest's offset holds.

## Open calls

### 1. Two corners where the default range comes out too wide

Background. The initial fit gives the forest's covariates straight-line terms so that they can compete with the
random effects for signal. Where there are almost as many such terms as there is data to estimate them from,
they take too much, and the range comes out wider than it should. Two cases. (i) Few groups and many covariates
that are constant within a group: with 8 groups and 10 such covariates the covariates reproduce the 8 group
means exactly, the fit hands the forest all the variation between groups (range 1.32 to 1.35 times the right
width), and the sampler's random-intercept sd came out 1.62 where the data were generated with 4, its interval
covering 4 in 7 of 19 fits; with 6, 4 and 2 such covariates the range is 1.22, 1.17 and 1.07 times the right
width. (ii) Nearly as many forest covariates as rows: with 159 at 200 rows the range is 1.42 times the right
width and expected-value error 6.8 percent [5.5, 8.1] above today's method; with 59, 1.16 times. The design's
idea for both, falling back to the response's own range whenever the fit puts a random effect's variance at zero
that the fit without the forest's columns does not, was screened for this plan on 20 data sets a cell: it would
have fired in 0 of 20 in each corner (1 of 20 with 6 group-level covariates), and the range it falls to is wider
still (1.63 and 1.78).
- (a) Leave the default as designed; the help names both cases and says to give a pair there. Nothing to build.
  Cost: the numbers above for a user who does not read it. How realistic: 8 groups with several group-level
  covariates is common (regions, sites), 10 for 8 extreme; 150 covariates on 200 rows with random effects is
  rare.
- (b) Where the fixed columns reproduce a grouping factor's indicators exactly (case (i) at its extreme, found by
  a rank comparison with no tuning number), use the range of the model's own fit, without the forest's columns:
  about 25 lines and a second fit of the same cost. In the 8-group cell that range gave a random-intercept sd of
  4.06, covering in 20 of 20. Cost: it does nothing until the covariates reach the number of groups, so the
  range jumps from about 1.3 to 0.62 of the right width between 6 and 7 covariates; it does nothing for case
  (ii); and that range is the one that loses 10 to 20 percent of split accuracy with many small groups.
- (c) The design's trigger with the response's own range. Cost: a second fit for a rule that did not fire.

Recommended: (a). A user in either corner has a way out in one argument, and neither rule on offer helps a user
who is near a corner without being in it. A measurement could decide between (a) and (b) only: 8 groups with 2,
6, 7 and 10 group-level covariates at 200 and 1000 rows, 20 seeds, chains of 5000 + 5000 (nothing with 8 groups
converged at the default length), (a) against (b) against today's method. Rule: adopt (b) if, in both cells where
it acts, the random-intercept sd's interval covers the generating value in more fits than (a)'s by more than the
difference between two runs of (a), and its expected-value error is not worse under the 3 percent rule.

### 2. What a binary fit starts from when a grouping factor has many levels

dec-B443 names the route; this is a cost its record does not carry, so it comes back as a question.

Background. A binary fit takes no range from the initial fit, only the values its first sweep starts from. With
lme4 no longer used, that start is a probit regression with each grouping factor entered as ordinary indicator
columns. That regression is dense in the number of groups. Run for this plan: at 8,000 rows in 2,000 groups of 4
it took 353 seconds where lme4's probit mixed model took 0.7; at 2,000 rows in 500 groups, 6 seconds for 0.15;
it holds rows times groups numbers, 20 GB at 100,000 rows in 25,000 groups. And once the start is put on the
linear-predictor scale, as it should be, groups whose outcomes are all 0 or all 1 are fitted exactly: with
groups of 4 the start is beyond 4 in size (a probability under 0.0001 or over 0.9999) on 41 to 44 percent of
rows, where the mixed model's runs from -0.7 to 1.3. Also run: the start hardly matters. Over 30 fits at 500
rows, swapping lme4's start for the indicator regression's (both on the probability scale, as today) moved each
row's fitted probability by a median 0.026 to 0.046 posterior sd, where two runs with the same start differ by
0.025 to 0.040.
- (a) As ruled. Cost: the minutes and the memory above for a user with thousands of groups, who today (with lme4
  installed) waits under a second; the extreme start, unmeasured.
- (b) Start from a probit regression on the fixed columns only, the random effects not in the start: about 5
  lines fewer than (a) (it is today's last resort), 0.01 seconds at any number of groups. Cost: a start that
  ignores the groups, which the run above suggests costs nothing after warm-up; not measured for (b) itself.
- (c) Fit a probit mixed model inside stan4bart by repeating the new routine on a working response: about 40
  lines and a second numerical routine to hold to lme4 in the tests. Cost: that upkeep, for one sweep's start.

Recommended: (b), if dec-B443's check passes for it. The maintainer wanted no run-time lme4 and a start that does
not depend on what is installed; (b) gives both at no cost in time. The check (Gates and baselines) is run with
(b)'s start in place of (a)'s, under the same rule; if (b) fails it and (a) passes, (a) is built with the help
saying how long it can take.

## Found against the record

- dec-B443 says stan4bart prefers lme4's helpers "at load". The choice is made when the package is installed.
  Run: a build made beside lme4 holds lme4's and reformulas's own functions; loaded where neither is installed it
  attaches, then stops at the first call (`could not find function "chk.cconv"`). A build made without them
  holds stan4bart's copies. H removes both operators, which also repairs this.
- The design found the fitted-probability start on the lme4 route. The route without lme4 has it too: `fitted`
  on a `glm` ignores `type = "link"` (run: values 1.6e-9 to 1 where the linear predictor runs -5.9 to 5.9).
- The design's fallback trigger for its two corners does not fire in them (Open call 1).

## Not verified

Nothing was built in stan4bart. The routine's numbers are the prototype's, to 5000 by 2000 crossed levels. Read
and not run: nlminb's handling of an Inf trial value, L-BFGS-B's refusal of one. Not run: the `offset_type` rule,
the zero-width rule, the old-fit restore, Open call 2's option (b), the helper on a build with
response-scale-rows, bartCause's suite and stan4bart's posterior gate, an older Matrix, Windows, a threaded BLAS.
