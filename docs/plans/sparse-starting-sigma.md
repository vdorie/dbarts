# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403); rewritten 2026-10-10 for dec-B433 and dec-B434. Not built. Three calls
are open ([Open calls](#open-calls)); the steps build the side of each that adds nothing.

agent: one sonnet implementer (R only, no engine code); one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a continuous response under a chisq residual prior with no sigest, in two cases: a
caller's sparse x whose design leaves residual degrees of freedom, at any size (dec-B403, dec-B433), and any
fit that keeps a factor categorical, the default (dec-B422). SHIFTING for an indicator expansion R built
sparse (the sparse QR's value to `lm.fit`'s). NEUTRAL for a dense design with no factor, a design with no
residual degrees of freedom, a fixed-unit family, a fit given a number, draws under a fixed residual prior.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~880 lines (R/utility.R ~200 net of the QR's removal, R/spec.R and R/tombstones.R ~50, R/xbart.R ~50,
packers, summary and R/rbart.R ~40, tinytest ~350, man ~60, docs ~60, NEWS ~6, MANIFEST one row).

## Goal

`sigest` is a number, one of "auto", "dense" and "sparse", or in xbart a function. Each string gives the
residual sd of the linear regression of the response on the predictors, factors as indicators, at any size,
and the response's sd where that leaves no residual degrees of freedom, announced under verbose and in
summary. "auto", the default, is `lm.fit` for a matrix or a data frame and the exact sparse routine for a
sparse x. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Ruled, each covering what it names: dec-B370 (an R-built indicator expansion "takes the same linear-model
estimate whichever the storage"); dec-B382, group 10 (xbart chooses the "route once from all rows"); dec-B403
(a sparse x takes the sparse regression, no warning; how it is computed is left to measurement); dec-B421 (its
rank tolerance may differ from lm's on nearly dependent columns); dec-B422 (factors enter as indicators);
dec-B427 (a zero-filled column beside the indicator of its unstored rows is re-centered); dec-B429 (a verbose
line and a summary line, no warning); dec-B430 (LAPACK's pivoted Cholesky); dec-B433 ("do the regression when
it is defined with no size limit, sd(y) otherwise"; xbart may take a function for each fold; the help names
alternatives "in case a user has to cancel at the initial estimate phase"); dec-B434 ("Rather than try to
guess, it would be better to give the user the option and let them decide", and "Sure, do the "auto" fix":
`sigest` takes "auto", "dense" or "sparse", "auto" follows what the caller passed, and there is no size
constant). Neither dec-B413 nor dec-B425 says whether an estimate is made under a fixed residual prior.

Today, on the worktree's build (bartcore da9bff44, shipped build, arm64 macOS, R 4.6.1, reference BLAS); "run"
is this revision's (scratch/sspfold4/), "carried" an earlier author's or the critique's, not rerun:
- [`estimateSigmaFromLinearModel`](../../R/utility.R) returns sd(y - offset) with a warning for a source with
  a CSC column, sends an expansion R built sparse to Matrix's sparse QR
  ([`sparseResidualStandardError`](../../R/utility.R)) and anything else to `lm.fit`, a categorical factor as
  its codes; [`floorSigmaEstimate`](../../R/utility.R) turns a non-finite result into sd(y), warning.
- Run (g-doors.out): at two chains `bart` estimates once and `rbart_vi` three times; `rbart_vi`'s defaults
  turn the control's verbose off; a fixed residual prior still estimates (slot 0.473 where it fixes 0.7).
- Run (i-auto.out): a string for `sigest` stops in `as.double` ("must be coercible to numeric type");
  `density(bw = "foo")` stops with "unknown bandwidth rule".

## Design

`sigest` at every door (`bart`, `bartBT`, `dbarts`, `dbartsSpec`, `rbart_vi`, `xbart`): NULL, and NA where a
door takes it, reads as "auto"; a string must be exactly one of the three names, any other stopping with
`unknown 'sigest' rule "<string>"; use "auto", "dense" or "sparse"`, as `density` refuses an unknown `bw`. A
resolver beside [`resolveSigestArg`](../../R/tombstones.R) returns the route, handed on to
[`resolveSamplerSpec`](../../R/spec.R). A string is a rule for an estimate that is made: `dbartsSpec` keeps a
data object's own sigma under any string. The retired spelling `sigma` takes a number only. The slot and
`fit$sigest` hold the number used; the call keeps the string. `estimateSigmaFromLinearModel(data, route)`
returns the sigma, with the attribute `fallback = TRUE` where no regression was defined (E).

A. The design, built by one new function, `startingSigmaDesign(x, route)`, for the all-rows estimate and
xbart's chunk runner; `data@x` never sees it. For the dense route it is a dense matrix, indicators built
dense, so Matrix plays no part; for the sparse route columns keep the caller's order and record their factor
term and whether the caller gave them dense (for C.6 and D.1).
1. A plain matrix with no factor marked: [`sigmaDesignMatrix`](../../R/utility.R) as today. One whose
   `varTypes` mark a factor is rebuilt as a frame (codes to factors over its level table) and takes case 2.
2. A dense container with a factor: [`makeIndicatorModelMatrix`](../../R/utility.R) on its frame.
3. A caller's sparse source (a bare dgCMatrix wrapped by [`wrapSparseTestMatrix`](../../R/mixedMatrix.R)): a
   factor column, dense-backed or a sparseFactor, emits case 2's columns, the reference level included; other
   columns as in [`sparseDesignMatrix`](../../R/utility.R). An indicators-route container as it is.

A missing value takes its column's observed mean (an indicator's, the level's frequency), as today; on the
dense route after the design is made dense, so there a sparse x with NA entries shares its dense twin's bits.

B. Route. n is the rows with a response, weight and offset and a positive weight; p the built columns.
1. An infinite entry is refused first, for every design and before any routine, naming the column.
2. "dense": `lm.fit` ([`residualStandardError`](../../R/utility.R)) on the dense indicator form, a sparse
   source made dense first. "sparse": C and D on the assembled sparse form, which needs Matrix.
3. "auto" is "sparse" where the caller supplied a sparse-stored column, that is a dgCMatrix, or a data frame
   or list holding a Matrix sparseVector, a dgCMatrix or a sparseFactor column (today's test:
   [`predictorSourceIsSparse`](../../R/utility.R) without the attribute `sparse.from.indicators`, which
   stays), and "dense" otherwise: a matrix or a data frame, factors included, whatever storage R gave their
   indicators and whether or not Matrix is installed. dec-B370 holds by construction, and without Matrix no
   caller can pass a sparse column ([`assembleMixedMatrix`](../../R/mixedMatrix.R) requires it).
4. Under every string a rank that leaves no residual degrees of freedom is the fallback (E).

Run (b-route2.out, i-auto.out): the data frame with n 5000, a 4,600-level factor and a numeric column gets
0.480200 from "auto" in 51.6 s (0.9-34: 0.4802 after 68 s; today's sparse QR 0.0 s) and from "sparse" in 0.17
s; n 3000 with 2,200 levels 0.519528 in 7.7 s and 0.15 s. On 40 small caller-sparse designs D is within
1.1e-15 of `lm.fit` and bit-identical on 7: a sparse fit shares its dense twin's draws only under "dense".

C. Front end of the sparse routine, in this order.
1. Rows: drop rows with a missing response, weight or offset and rows of weight 0; z = y - offset.
2. Stored means nonzero: explicit zeros are dropped from the design, not from `data@x` (ten defeat step 4).
3. Constants: a column whose entries over the kept rows are all equal (zeros included) is dropped.
4. Re-centering: a column with unstored kept rows and unequal stored entries has their mean taken off them
   when another column is constant on exactly its unstored rows and unstored elsewhere, or constant on exactly
   its stored rows (found by stored count and sum of row indices, confirmed row for row).
5. Same pattern: a column whose stored rows are exactly an earlier column's (found and confirmed likewise) has
   its projection on that column taken off its stored entries, b - a (a'b / a'a); the span is unchanged.
6. Centering: each numeric column the caller gave dense and each column stored in every kept row is centered
   at its weighted mean. Then every column is divided by its largest absolute entry.

D. The sparse routine, on B = diag(sqrt(w)) [1 X] with each column divided by its norm.
1. Block elimination. Among the factor terms, widest first, take the first whose columns share no kept row
   other than rows imputed for that factor. Its block of the crossproduct is diagonal (plus c v v' where rows
   were imputed, v their common row and c their weight; inverted by the Sherman-Morrison formula), so with Z
   the block and C the other columns, S = C'C - C'Z (Z'Z)^-1 Z'C. With no such term C is B.
2. Narrow (ncol(C) < n): `chol(S, pivot = TRUE, tol = 1e-10)` inside `suppressWarnings` (its one warning is
   the rank deficiency); r is the block's width plus S's rank. C's pivoted columns take their coefficients b
   from the factor, the block (Z'Z)^-1 (Z'z - Z'C b), and the residual e is formed directly.
3. Wide (otherwise, nothing eliminated): the pivoted Cholesky of `as.matrix(Matrix::tcrossprod(B))` at the
   same tolerance; Q = `qr.Q(qr(L, LAPACK = TRUE))`, L the factor's first r columns, unpivoted, and e = z
   sqrt(w) - Q Q' z sqrt(w). Either way sigma = sqrt(sum(e^2) / (n - r)), and none when r >= n.

Cost. "sparse", with m the smaller of n and one more than the columns left after D.1: about 1 s at m 2,000, 2
minutes and 1.8 GB at 10,000, at least 7 hours at 50,000, an allocation error at 100,000 (dec-B433's figures;
0.8 s at 2,000 and 10.2 s at 4,602 run); a design whose width is one factor does not reach them (run,
h-absorb*.out): 0.4 s at n 2e5 with 50,000 levels, 13.7 s with two crossed 5,000-level factors. "dense" costs
what `lm.fit` costs and holds n by p doubles, 80 GB for that 50,000-level frame. Near n 3e6 D's rank can
exceed `lm.fit`'s by an exact dependency (2 of 4 designs; run), moving sigma by 1.7e-7.

Where D differs from `lm.fit` ("dense" gives lm's handling), at any size: in dec-B421's band, a near-copy of a
column that carries signal and rows weighted 1e12 times the rest (carried); by the plan's call, a value
recorded on the rows of two or more levels of a factor with a large offset and a small spread (1.325 for
0.501; run). D can be lower: `lm.fit` drops a time within 60 s beside its indicator, D keeps it (run).

E. No regression defined (dec-B429). Under a fixed residual prior none of E.2 to E.5 applies.
1. Value: sd(y - offset) through [`floorMarginalSigma`](../../R/utility.R), as today. No warning.
2. Record: the sigma value carries the attribute `fallback = TRUE`, which rides the `sigma` slot through
   `dbarts`, [`setData`](../../R/bartcore.R) and a saved sampler, changes no draw (run) and goes when the slot
   is rewritten. The fits store `sigest` without it and a logical `sigest.fallback`.
3. Verbose: a dbartsStartingSigmaMessage in [`announceAutoFamily`](../../R/utility.R)'s form, `starting sigma
   is the sd of the response (<value>): the linear model on <p> columns (factors as indicators) leaves no
   residual degrees of freedom in <n> rows; supply 'sigest' to set it`, sent once from `resolveSamplerSpec`
   under the door's verbose flag; by `rbart_vi` from its own `verbose` for its validation sampler, so its
   default call prints it ([`rbart_vi_fit`](../../R/rbart.R) muffles the class); by `xbart` for all rows only.
4. Before every estimate, by the same flags and class: `estimating the starting sigma by a <dense|sparse>
   linear regression on <n> rows and <p> columns; see 'sigest'`.
5. Summary: [`printSummaryBartBody`](../../R/diagnostics.R) prints `(Starting sigma: the sd of the response;
   the linear model had no residual degrees of freedom)` on `sigest.fallback` (hurdle: positive part).

F. xbart.
1. Route, once from all rows (dec-B382): the all-rows estimate runs once, raising any refusal and sending E.3
   and E.4, and names the routine (none, `lm.fit` or D, by B on what the caller passed). Each fold runs that
   routine on its training rows, the design built once per chunk in [`xbartRunUnits`](../../R/xbart.R), or
   takes its rows' sd where the route is none or its own rows leave no residual degrees of freedom.
2. `sigest` may also be a function, `function(x, y, weights, offset)`, called by name once per unit on the
   worker that runs it, never from an engine thread. `x` is the training rows of A's design for "auto" (a
   numeric matrix, or a dgCMatrix where the caller passed sparse; no NA), `y` the response, `weights` and
   `offset` the training rows' or NULL where the call has none, all after C.1's row rule. It returns one
   positive finite number; anything else stops the run with `'sigest' function must return one positive finite
   number`. At the door a function with neither `...` nor all four names is refused, as is one beside a fixed
   residual prior; under a binary family it is not called. No all-rows estimate then.
3. Workers. A socket or cluster worker gets the function serialized, as `loss` is
   ([`xbartLossFunction`](../../R/xbart.R)), so it names packages with `::`. R's generator is seeded per unit
   ([`withFixedSeed`](../../R/validateComposition.R)) from one more seed drawn after xbart's existing ones, so
   a function that samples rows reproduces at any `n.threads`.

## Change

By file: R/tombstones.R and each door take the route; R/utility.R takes A to D and E's helpers,
`sparseResidualStandardError` keeping its name and signature as C and D and returning NA where no regression
is defined; `floorSigmaEstimate` loses its warning. R/spec.R: [`nonFinitePredictorNames`](../../R/spec.R)
reads sparse sources and runs before the estimate. R/bart.R, R/rbart.R, R/diagnostics.R: E. R/xbart.R: F, the
returned value checked as [`validateSigest`](../../R/dbarts.R) checks a number, plus finiteness. Both warning
classes go: never released, and no branch of stan4bart, bartCause, treatSens or bairrtt names them (carried).

Constraints: no engine, bridge, C API, state or class change; no new dependency; no size constant. A dense
design with no factor under a chisq prior is bit for bit unchanged. Out of scope: a string "sd", a summary
method for `rbart_vi`, one estimate shared by its samplers.

## Tests

Edits, six files: test-starting-sigma.R (three pins on a factor's codes take the indicator design);
test-data-mixed.R and [test-indicator-storage.R](../../inst/tinytest/test-indicator-storage.R) (the
sparse-fallback blocks: no warning, the dense `sigest` within 1e-10; the latter's direct calls to
`sparseResidualStandardError` stay at 1e-8, except the loop under
["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R), whose near-copies D drops: it
goes through the estimate under "auto", `expect_identical` to the dense fit); test-boundary-inputs.R and
test-xbart-fold-oracle.R (no warning, `as.numeric(data@sigma)` as pinned, the attribute, the message under
verbose); test-sampler-splitProbabilities.R (a seeded comparison of split counts on a factor design fails:
report the counts over five seeds; move the seed if it holds, else stop; the same for any seeded expectation
on a caller-sparse fit, whose sigest moves in its last digits). test-data-code-channel.R pins A.1 unchanged.

New, in test-sparse-starting-sigma.R, against `lm.fit` on the dense design within 1e-10 unless said:
- dec-B422: a frame with a 5-level factor and a numeric column, default against `factors = "indicators"` and
  its matrix of codes, `sigest` identical and equal to `summary(lm(y ~ f + x))$sigma`; the factor ordered,
  with two levels, with 40, with a missing value; each also "sparse".
- dec-B434: "auto" runs `lm.fit` for a matrix and for a data frame (traced; `identical` under
  `options(dbarts.sparseIndicators = "dense")` and "auto") and D for a bare dgCMatrix, a sparseVector column
  and a sparseFactor column, through `dbarts`, `bart` and `xbart`, with no warning (counted, Gate hygiene);
  "dense" on those three is `identical` to the dense twin, NA entries included; "sparse" on the frame is
  `identical` under both storages. By each door: NULL is "auto"; "Dense" and "sd" stop naming the string;
  `sigma = "dense"` stops; `fit$sigest` is a number, the call keeps the string; `dbartsSpec` keeps a carried
  sigma.
- "sparse": one-hot 3 x 40 plus three duplicated columns with zero weights and an offset; a column scaled by
  1e160; p > n with rank below n - 1 and of full rank (the fallback); a fully stored timestamp-like column
  with signal, weighted and not, within 1e-8 (`lm.fit` errs by 3.4e-9; carried).
- dec-B421's band edges, "sparse", n 400 (run, f-tol.out): a near-copy differing by 1e-4 of its norm with the
  signal on the difference equals `lm.fit` (0.5085; at tolerance 1e-6, 0.8155); one differing by 1e-6 equals
  `lm.fit` without the copy (0.7933; at tolerance 1e-16, 0.4818); "dense" equals `lm.fit` on both.
- Front end, "sparse", within 1e-8: the timestamp pair in both orders, with ten explicit zeros stored in the
  time column and ten in the indicator; a value stored on a first level's rows; start and end times sharing a
  missing pattern. A partner matching in count and index sum but not in rows leaves the design untouched.
- D.1 and A.3, "sparse": one to three factors with numeric and sparse columns, a duplicated factor, weights,
  an offset and missing factor values; a sparseFactor with a first and with a middle reference level; a level
  per row is the fallback; `chol` sees only the columns left.
- dec-B429, dense and sparse, by `dbarts`, `bart`, `bartBT`, `rbart_vi` (two chains on one thread; four on two
  at home only) and `xbart`: zero warnings; one message with the flag TRUE and none with it FALSE; `sigest` is
  `sd(y)` with no attribute; `sigest.fallback` and the summary's line, absent for a linear estimate, a
  supplied sigest and a fit with the field removed; after `setData` and after `data@sigma <- 1` the record
  follows the value; none of it under a fixed prior. E.4's message once before each estimate.
- An Inf entry in a dgCMatrix and in a sparseVector column: the dense path's error, no routine called.
- xbart: each fold's sigma equals `lm.fit`'s on its rows' indicator design; every fold runs the all-rows
  routine (traced, under each string); an all-rows fallback fits no fold. `sigest` as a function: each unit's
  sigma is its value on that unit's training rows (arguments recorded and compared), `weights` and `offset`
  NULL where none; one drawing `sample` gives identical losses in two seeded runs and at `n.threads` 1 and 2
  (at home only) and leaves `.Random.seed` alone; a return of NA, 0, Inf, a character or length 2 stops;
  `function(y)` and a function beside a fixed prior are refused.

Reviewer's mutants, each of which must fail a test: the expansion removed for a container, a matrix of codes,
a fold; A.3 without its first level's indicator; "auto" sending an R-built sparse expansion to D, or a
dgCMatrix to `lm.fit`; "dense" ignored; an unknown string read as "auto"; explicit zeros kept; C.4 without its
partner test or its row-for-row confirmation; C.5 removed; tolerance 1e-10 to 1e-16 and to 1e-6; D.1's block
left out of the rank, its rank-one term dropped, and applied to a term whose columns share a row; the Inf
check after the routine; the message sent with the flag FALSE, or once per chain in `rbart_vi`; the attribute
left on `fit$sigest`; the function called on all rows, its return unchecked, its seed not set.

## Baselines and gates

- What moves (carried, scratch/sspnew/eq.out, and re-derived for dec-B434, not rerun in the harness): twelve
  of [equivalence.R](../../benchmarks/R/equivalence.R)'s 55 scenarios. POSTERIOR-CHANGING: sparse,
  mixedmatrix, sparsefactor, testswap, leaffactormixed, factorpartial, xbartmixed (a caller's sparse source,
  sd(y) to D's value, 0.16 to 0.62 of sd(y) and `lm.fit`'s to rounding) and categorical, leaffactor, nafactor,
  ordfactor (a dense frame's factor, codes to `lm.fit` on indicators, 0.9-34's value). SHIFTING:
  wideFactorIndicators (the sparse QR to `lm.fit`, 6.6e-16). NEUTRAL: the other 43, bcf's 15, multinomial's
  11, the snapshot files, every exact gate. Another mover is a stop.
- Re-record the twelve on the reference build into a copy of equivalence-e4faed5c
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)). Partition: 43 of 55 identical, the movers the twelve, no
  |z| above 4 on wideFactorIndicators; 55 of 55 bitwise from a second `--preclean` install.
- Oracle (MANIFEST rule P17): the identity is agreement with `lm.fit` on the indicator design: the tinytest
  pins, and a sweep of "sparse" over 400 random designs (n 30 to 2000, p 5 to 1530, dependent columns, zero
  weights, factors) and 30 tall ones (n 1e5 to 3e6, one to three factors of 6 to 5,000 levels).
- Independently of the implementer: the posterior-changing battery of
  [RNG classes and their gates](README.md#rng-classes-and-their-gates), `R CMD check --as-cran`, the lint set.
- Speed and memory, within 1.5x of: D at n 2600, p 1999 (0.8 s), at n 2e5 with a 50,000-level factor (0.4 s,
  300 MB) and with two 5,000-level factors (13.7 s, 770 MB).

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd, man/rbart.Rd and man/xbart.Rd: a number, or how to
  estimate one, the residual standard deviation of a linear regression of the response on the predictors,
  factors as indicator columns whatever `factors` says. `"dense"` is `lm`'s own fit on the dense design.
  `"sparse"` is an exact sparse regression that never forms the dense design; it agrees with `lm` to rounding
  except on nearly dependent columns (a near-copy of another column, a tiny spread beside a large mean), where
  it may drop a column `lm` keeps, and its last digit follows the BLAS and thread count (dec-B430). `"auto"`,
  the default, is `"dense"` for a matrix or a data frame and `"sparse"` for a sparse x. Pass `"sparse"` for a
  data frame with a factor of thousands of levels (0.17 s for 51.6 s at 4,600 levels; run); pass `"dense"` for
  `lm`'s handling of nearly dependent columns or for a sparse fit's draws to match its dense twin's. Where the
  regression leaves no residual degrees of freedom the response's standard deviation is used, which a verbose
  fit prints and `summary` shows. No size turns the estimate off: `"sparse"` takes about a second at 2,000
  rows and columns, minutes at 10,000, hours at 50,000, and can fail to allocate; `"dense"` holds rows times
  columns doubles; R cannot interrupt either (Open call 1). For an expensive estimate pass a number: `sd(y)`,
  as several BART packages do, too high where there is linear signal; where rows far outnumber columns a
  regression on a random subset of rows, `i <- sample(nrow(x), 10 * ncol(x)); summary(lm(y[i] ~ as.matrix(x[i,
  ])))$sigma`; for a wide or sparse x a cross-validated lasso, `sqrt(min(glmnet::cv.glmnet(x, y)$cvm))` (each
  run as written, d-examples.out: 2.949, 0.5033, 0.5012 for the default's 0.4997).
- man/xbart.Rd: `sigest` may be a function (F.2, F.3), with the alternatives as functions. man/bart.Rd: both
  warning classes go, the value lists `sigest.fallback`. man/dbartsData.Rd: the `sigma` slot's attribute.
- inst/NEWS.Rd: the classes go; one item for the fallback's line, one for xbart's function. docs/design: a
  dated paragraph in sparse-columns.md's [R surface](../design/sparse-columns.md#r-surface) and a dated line
  in starting-sigma-sensitivity.md's [Recommendation](../design/starting-sigma-sensitivity.md#recommendation);
  error-style.md and memory-footprint.md follow. At landing: TODO, the ledger entry, Status, the MANIFEST row.

## Steps

1. Design A to F with the tinytest edits and new tests; the suite green against `R CMD INSTALL -l <lib> .`.
2. The sweep; the speed points. 3. Help, NEWS, docs. 4. After review: the re-record, in its own commit.

Stop and report when: the diff passes ~1150 lines; a scenario other than the twelve, a snapshot or an exact
gate moves, or wideFactorIndicators shows a |z| above 4; the sweep shows a rank below `lm.fit`'s, a rank above
it at n under 1e6 or by more than the design's exact dependencies, or a sigma off by more than 1e-10 plus
(rank excess) / (n - r) above 30 residual degrees of freedom; a speed or memory point is past 1.5x; a tinytest
fails that Tests does not name; the change needs engine, bridge or C API code.

## Calls made

- D.1, block elimination (about 25 lines; dec-B403 leaves how the fit is computed to measurement), serves the
  case the help sends to "sparse", a data frame with one wide factor, and a sparseFactor column: 0.17 s for
  10.4 s at 4,600 levels (run), and by D's cost for 7 hours at 50,000. Rejected: today's sparse QR as the
  routine for factor designs (0.1 s there, unfinished after 170 s on crossed factors D does in 14 s; run).
- NULL stays the formals' default and is "auto". Rejected: `sigest = "auto"` in six signatures, moving codoc
  and every test of an absent sigest. A string is a rule: a fixed residual prior does not refuse it.
- C.5 is not in dec-B427's words: 8 lines and 0.02 s at n 2e5 take start and end times that share a missing
  pattern from 1.760 to `lm.fit`'s 0.4995 (run). Rejected: accepting 3.5 times, for spreads up to 6 hours.
- A.3 emits every present level, as A.2 does, so a value recorded on a first level's rows finds its partner
  (0.5002, `lm.fit`'s, against 1.005 with that level dropped; run).
- The record is an attribute of the sigma value: about 6 lines, no class change, never stale when stan4bart
  writes the slot. Rejected: a logical slot with a setter and a guarded reader, stale after a direct write.
- D.3 has no row equilibration: none of 150 wide designs, weights ranging over 1e14, needs it (run).
- xbart's function is `function(x, y, weights, offset)`, called by name: `lm.wfit`'s and `glm.fit`'s argument
  names, as `glm` calls a function given as its `method`. Rejected: `function(data)` handed the fold's data
  object. It is seeded per unit because two of the help's three alternatives draw random numbers.
- xbart announces for all rows only. E.4 prints before every estimate, no size constant being left to gate it
  (`family = "auto"` already messages such fits). Under a fixed residual prior E is silent.

## Open calls

### 1. Is an estimate that cannot be interrupted acceptable?

Background. dec-B433 has the help name alternatives "in case a user has to cancel at the initial estimate
phase". Run (e-interrupt.out): each long step was sent Ctrl-C's signal two seconds in, and in a second run
given a 2 second `setTimeLimit`. None stopped: `lm.fit`'s QR ran its 14 s (n 4000, p 2500), the sparse
crossproduct its 15 s (4e5 rows, 3,000 columns), LAPACK's pivoted Cholesky its 22 s (6,000 columns); a control
loop stopped in 0.1 s. R acts on an interrupt only between compiled calls, as with `lm`, so cancelling means
killing the R session. The steps that grow: the Cholesky ("sparse"; 2 minutes at m 10,000, at least 7 hours at
50,000, for a sparse matrix with that many rows and columns) and `lm.fit` ("dense"; 51.6 s on B's data frame).

Options. (a) Leave it: nothing built; the help says so and E.4's line shows what is running. (b) For "sparse",
the crossproduct in 256-column blocks and the pivoted Cholesky written in R in 96-column blocks (LAPACK's
algorithm), which R interrupts between blocks: about 40 lines and 25 of tests. Run (e-blocked.out,
e-check.out, e-mem.out): the signal stopped both 0.0 to 1.1 s later; a time limit stopped the Cholesky in 0.2
s and once ran through the crossproduct, which needs an explicit poll; rank equal and sigma within 7e-16 of
LAPACK's on 33 designs; the Cholesky 0.9 s against LAPACK's 0.8 at 2,000 columns, 3.6 against 10.2 at 4,602,
7.2 against 22.8 at 6,000; the crossproduct 18.8 s against 16.1; memory as prototyped 1,330 MB against 285 MB
over a 275 MB matrix. It replaces the routine dec-B430 kept on an estimate of "tens of seconds at 2000 columns
against about one" that this run does not bear out. Not run: another BLAS, thread counts, 10,000 columns.
"dense", as `lm`, and the wide branch's QR stay uninterruptible. (c) The estimate in a forked child the
session can kill: about 12 lines, no numerical change, both routes, 22.5 s against 22.8, stopped in 0.0 to 1.2
s; not on Windows, nor under RStudio, Positron or R.app, where `xbart` already declines to fork. Recommended:
(b). The ruling's remedy assumes a user can cancel, and (b) costs no time on this BLAS.

### 2. Does `sigest` take a function at the other doors?

Background. dec-B433 names xbart, where a function is called once per fold. At `bart`, `bartBT`, `dbarts`,
`dbartsSpec` and `rbart_vi` it would be called once, on all rows, and `bart(x, y, sigest = f(x, y))` does the
same in one line; a caller lacks only the design with factors expanded, which the exported
`makeModelMatrixFromDataFrame` builds. Options. (a) xbart only: nothing more. (b) Every door: about 25 lines,
help in four more files, about 30 lines of tests; `sigest` then means one thing everywhere, and `rbart_vi`
would call it once per sampler (five times at its defaults). Recommended: (a).

### 3. Is the estimate made under a fixed residual prior?

Background. Under `gaussian(sigma = fixed(value))` the estimate calibrates nothing and is made today: the fit
reports it as `sigest` (0.473 where the fixed sigma is 0.7, run) and an Inf in x is refused by it. With no
size limit it can cost minutes to hours there, and a caller cannot skip it: a number beside a fixed prior is
accepted only where it equals the fixed sigma, and refused from 1.1-0 (dec-B413). Options. (a) As today, which
the steps build: the wasted estimate, by the route `sigest` names. (b) No estimate: about 10 lines; `sigest`
and the slot report the fixed sigma (0.7 for 0.473 above); an Inf in x is accepted, as it is today when sigest
is given; draws are unchanged (carried: identical with 0.7 forced into the slot). Recommended: (b).
