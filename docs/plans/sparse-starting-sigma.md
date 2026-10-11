# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: LANDED 2026-10-10 as 8c27f225..d820097c (dec-B403, dec-B421, dec-B422, dec-B429, dec-B430, dec-B433, dec-B434,
dec-B438, dec-B439, dec-A200). See the [Landing note](#landing-note).

agent: one sonnet implementer (R only, no engine code); one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a continuous response under a chisq residual prior with no sigest, in two cases: a
caller's sparse x whose design leaves residual degrees of freedom (dec-B403, dec-B433), and any fit that keeps a
factor categorical, the default (dec-B422). SHIFTING for an indicator expansion R built sparse (the sparse QR's value
to `lm.fit`'s). NEUTRAL for a dense design with no factor, a design with no residual degrees of freedom, a fixed-unit
family, a fit given a number, draws under a fixed residual prior.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~1000 lines (R/utility.R ~270, doors ~70, R/xbart.R ~50, E ~30, tinytest ~420, man ~70, docs ~60, NEWS ~6).

## Goal

`sigest` is a number, one of "auto", "dense" and "sparse", or in xbart a function. Each string gives the residual sd
of the linear regression of the response on the predictors, factors as indicators, at any size, and the response's sd
where that leaves no residual degrees of freedom, announced under verbose and in summary. "auto", the default, is
`lm.fit` for a matrix or a data frame and the exact sparse routine for a sparse x. The tier is "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Context

Ruled, each covering what it names and cited where applied: dec-B370 (an R-built indicator expansion takes the
linear-model estimate whichever its storage), dec-B382 group 10, dec-B403 (a sparse x takes the sparse regression, no
warning; how is left to measurement), dec-B421 (the rank band), dec-B422 (factors as indicators), dec-B427 (the
re-centering), dec-B429, dec-B430, dec-B433 (the regression "with no size limit, sd(y) otherwise"; xbart may take "a
function to use for each fold"; the help names alternatives "in case a user has to cancel") and dec-B434 (the three
strings; "auto" follows what the caller passed). None says whether an estimate is made under a fixed residual prior.

"Run" is this revision's, in scratch/sspfold5/, on a build of bartcore 8288ac6a (the estimate's and the doors' sigest
code is unchanged at the tip), arm64 macOS, R 4.6.1, R's own BLAS unless said; "carried" is an earlier round's or the
fifth critique's (scratch/sspcrit5/), not rerun. Today's estimate is
[`estimateSigmaFromLinearModel`](../../R/utility.R), [`sparseResidualStandardError`](../../R/utility.R) and
[`floorSigmaEstimate`](../../R/utility.R). Run (h-gdoors.out, h-doors.out, h-doors2.out, h-cran.out): `bart` estimates
once at two chains, `rbart_vi` three times; a string or a function for `sigest` fails inside `as.double` or `is.na`, a
numeric string is coerced, as in 0.9-34; a number beside a binary response warns at `bart` by the argument's name (a
written-out NULL too) and is silently ignored elsewhere; `dbarts` and `bart` estimate over a data object's own sigma,
`dbartsSpec` keeps it unless given a number; a verbose fit sends one message, the `family = "auto"` line, at `bart`,
`dbarts`, `xbart` and one-thread `rbart_vi`, none at `bart` with the family named, at `bartBT` (the engine's lines, on
standard output) or at several-thread `rbart_vi` (verbose off).

## Design

`sigest` at every door (`bart`, `bartBT`, `dbarts`, `dbartsSpec`, `rbart_vi`, `xbart`) is resolved after
[`resolveSigestArg`](../../R/tombstones.R), by [`resolveSigestRule`](../../R/utility.R), into a number or a route
handed to [`resolveSamplerSpec`](../../R/spec.R).
The slot and `fit$sigest` hold the number used, with no attribute; the call keeps the string.
- NULL, and NA where a door takes it, is "auto"; the retired spelling `sigma` takes a number only.
- A string is exactly "auto", "dense" or "sparse"; one `as.double` reads as a number is that number, as today; any
  other stops with `unknown 'sigest' rule "<string>"; use "auto", "dense" or "sparse"`.
- A function, except at xbart (F), stops with `'sigest' argument to <door> must be a number, "auto", "dense" or
  "sparse"; only xbart takes a function`. "sparse" where [`matrixAvailable`](../../R/utility.R) is FALSE stops with
  `sigest = "sparse" requires the Matrix package`.
- Where `sigest` has no use (a fixed-unit family) a known string is treated as a number is today: `bart` warns by the
  argument's name, the other doors say nothing.
- A data object's own sigma: `dbartsSpec` keeps it under NULL and "auto" and estimates over it under "dense" or
  "sparse", as a number replaces it; the other doors estimate over it, as today.

A. The design: `startingSigmaDesign(x, route)`, new, for the all-rows estimate and xbart's chunks, never `data@x`.
1. Dense route: a dense matrix, with no use of Matrix. A plain matrix with no factor marked is
   [`sigmaDesignMatrix`](../../R/utility.R)'s. A source holding a factor (a container's factor column, a column whose
   `varTypes` mark one, a sparseFactor) is rebuilt as a data frame, codes to factors over the level table, and
   expanded by [`makeIndicatorModelMatrix`](../../R/utility.R) with `storage = "dense"`; an indicators-route container
   and a caller's sparse columns are made dense. A missing value then takes its column's observed mean.
2. Sparse route: a list of `X` (a dgCMatrix with no NA, columns in the caller's order), `term` (per column, its
   factor, else NA), `dense` (per column, whether a numeric column the caller gave dense) and `imputed` (per factor
   with a missing value, its missing rows and a frequency per column). By source column:
   - A factor: with counts the `tabulate` of its codes, the levels emitted are those of positive count, the reference
     too, or the higher alone where two are present ([`sparseFactorIndicatorSlices`](../../R/utility.R)'s rule); the
     columns are one `Matrix::sparseMatrix(i, j, x = 1)`, i the rows observed at an emitted level and j that level's
     position, never `makeIndicatorModelMatrix` (45 s for 0.01 s at 50,000 levels; carried). A row with a missing code
     stores nothing; `imputed` holds it, with count / sum(counts) per level.
   - An indicators-route container's indicator columns: its factors and their counts are attribute `drop`'s integer
     entries, a factor's columns numbering its positive counts (one where two); stored 1s are kept; rows holding NA
     are the factor's missing rows, not stored; frequencies as above.
   - Any other column: zeros unstored, a missing value at the mean of the column's observed entries, zeros counted, as
     retired: [`sparseDesignMatrix`](../../R/utility.R) did.

B. Route. n counts the rows with a response, weight and offset and a positive weight.
1. An infinite entry is refused first, for every design and before any routine, naming the column.
2. "dense" is `lm.fit` ([`residualStandardError`](../../R/utility.R)) on A.1; "sparse" is C and D on A.2; under
   either, a rank leaving no residual degrees of freedom is E.
3. "auto" is "sparse" where the caller supplied a sparse-stored column - a dgCMatrix, or a data frame or list holding
   a Matrix sparseVector, a dgCMatrix or a sparseFactor (today's test: [`predictorSourceIsSparse`](../../R/utility.R)
   without the attribute `sparse.from.indicators`, which stays) - and "dense" otherwise, whatever storage R gave a
   factor's indicators and whether or not Matrix is installed (without it no caller can pass a sparse column:
   [`assembleMixedMatrix`](../../R/mixedMatrix.R)).

C. Front end of the sparse routine, in this order. "Held" columns are those of a factor that still has imputed rows:
C.3 (c) to C.6 leave them alone, and they are never the other column in C.4 or C.5.
1. Rows: keep the rows B counts (in a fold, among its training rows, the design having been built on all rows); z is y
   less the offset. Each factor's imputed rows are cut to the kept rows, its frequencies unchanged.
2. Stored means nonzero: explicit zeros are dropped from the design, not from `data@x`.
3. Emptied levels, then constants. (a) A factor with imputed rows and a column left with no stored entry (a level
   whose rows were all dropped): the first such column becomes the indicator of the imputed rows and leaves the
   factor, the others are dropped, and the factor has no imputed rows from here. (b) A factor with imputed rows and
   one column has its frequency written into the column on those rows, and has none from here. (c) A column whose
   entries over the kept rows are all equal (zeros included) is dropped.
4. Re-centering: a column with unstored kept rows and unequal stored entries has their mean taken off them when
   another column is constant on exactly its unstored rows and unstored elsewhere, or constant on exactly its stored
   rows (found by stored count and sum of row indices, confirmed row for row).
5. Same pattern: a column with unstored kept rows whose stored rows are exactly an earlier column's (the earliest;
   found and confirmed likewise) has its projection on that column taken off its stored entries, `b - a (a'b / a'a)`.
   Where what is left is at most 1e-7 of b's norm, lm's tolerance, the column is dropped.
6. Centering: each numeric column the caller gave dense and each column stored in every kept row is centered at its
   weighted mean. Then every column is divided by its largest absolute entry (none is zero after C.3 and C.5).

D. The sparse routine, with `h = sqrt(w)` and `zs = h z`; it ends in e, and `sigma = sqrt(sum(e^2) / (n - rank))`.
1. Block. The factor with the most columns left, at least two, the first on a tie, is eliminated. Z0 is its columns
   times h and d their squared norms; its true columns are `Z = Z0 + u v'`, u being h on its imputed rows and 0
   elsewhere, v its frequencies and `c = u'u` (both 0 with none). Ct is h times [1, the other columns], plus one
   column per other factor with imputed rows, h on them. The true other columns are `C = Ct W`, W being the identity
   over one row per such factor that holds its frequencies in its columns; their squared norms are
   `N_j^2 = |Ct_j|^2 + W_tj^2 |Ct_t|^2`, and `F = W diag(1 / N)`. With `a = v / d`, `s = v'a`, `den = 1 + c s`,
   `M0 = Z0'Ct`, `t = Ct'u` and `g = M0'a`,
   `St = Ct'Ct - M0' diag(1 / d) M0 + (c / den) g g' - (g t' + t g') / den - (s / den) t t'` and `S = F' St F`. If no
   factor has two columns, or S would have n or more, none is eliminated: each imputed row left is written in at its
   frequency and B is `h [1 X]`. If B has under n columns step 2 takes `S = D B1'B1 D`, B1 being B without its
   intercept, the crossproduct formed first and D the inverse roots of its diagonal, and orders the intercept last;
   else step 3 takes B with unit-norm columns.
2. Narrow. The rank r of S is 0 when no diagonal entry exceeds 1e-10, and otherwise that of
   `chol(S, pivot = TRUE, tol = 1e-10)` inside `suppressWarnings` (its one warning is the rank deficiency; LAPACK
   never tests its first pivot). K is the first r pivots, R the factor's leading r by r block. The design's rank is
   the block's columns plus r; at n or more no regression is defined (E). With `Ai(x) = x / d - c a (a'x) / den`,
   `hz = Z0'zs + v (u'zs)` and `k = Ai(hz)`: `b = R^-1 R^-T (F'(Ct'zs - M0'k - t (v'k)))[K]`, `beta = F[, K] b`,
   `bZ = Ai(hz - M0 beta - v (t'beta))` and `e = zs - Z0 bZ - u (v'bZ) - Ct beta`. With nothing eliminated
   `[ez e1] = [zs h] - P [zs h]`, P projecting on the kept columns `B1 D[, K]` through R; the intercept joins, one
   more in the rank, where `|e1|^2 > 1e-10 |h|^2`, and then `e = ez - e1 (e1'ez / e1'e1)`, else `e = ez`. An
   intercept that is exactly the sum of a complete set of indicators is so found whatever n and the weights: the
   difference of sums `S_00 - |R^-T S_K0|^2` it replaces grew with n past the tolerance (fix round, B1).
3. Wide. The same rule and call on `as.matrix(Matrix::tcrossprod(B))`; with L the factor's first r columns, unpivoted,
   `Q = qr.Q(qr(L, LAPACK = TRUE))` and `e = zs - Q Q'zs`; no regression at r >= n.

scratch/sspfold5/proto5.R is A.2 for a data frame or a sparse matrix, C and D, in base R and Matrix
(`sparseRouteDesign`, `planFrontEnd`, `planRoutine`); where it and these words differ, it is what was run. Run against
`lm.fit` on the dense indicator form: rank equal and sigma within 1.5e-13 on the critique's 300 adversarial designs
(c-sweep.out), no difference on 63 named cases (c-named.out), rank equal on 300 designs of one factor alone
(c-firstpivot.out); dec-B433's two frames give `lm.fit`'s 0.480200 and 0.519528 in 0.2 s for its 56 s and 8 s
(e-cost.out); 40 small caller-sparse designs are within 1.1e-15 and bit-identical on 7, so a sparse fit shares its
dense twin's draws only under "dense" (f-premises5.out).

Cost (run: e-cost.out, g-blocked-blas.out). "sparse" costs the Cholesky of S, 117 s at 10,000 columns under R's BLAS
and 3.5 s under Accelerate, and twice 8 bytes times the columns squared; one wide factor never reaches it (0.8 s and
240 MB from the data frame at n 2e5 with 50,000 levels, missing values or not). "dense" costs what `lm.fit` costs and
holds n by p doubles. D differs from `lm.fit`, at any size: in dec-B421's band, on a near-copy of a column that
carries signal and on rows weighted 1e12 times the rest (carried); on a value recorded on the rows of two or more
levels of a factor with a large offset and a small spread (1.325 for 0.501; run); on two complete sets of a caller's
own one-hot columns under weights that are not whole numbers, by one rank and 8e-8 in sigma on 1 of 30 designs at n
6e6 and none of 36 at 1e6 and 3e6 (fix round; one set, or unit weights, is exact at every n run, to 3.2e7; a frame's
factor does not show this, but where a factor is eliminated and the weights take few values (0.1, 1, 7) what is left
of the intercept grows with n, 7.4e-11 at 1.2e7 rows, and would cross the 1e-10 tolerance near 1.6e7, the rank right
on 30 of 30 designs up to 1.2e7; before that round one set alone gave a rank too many from n 21,000, on 31 of 120
designs, and this sentence put it near n 3e6); and downward where `lm.fit` drops a time within 60 s beside its
indicator and D keeps it (run).

E. No regression defined (dec-B429). Under a fixed residual prior none of E.2 to E.6 applies.
1. Value: `sd(y - offset)` through [`floorMarginalSigma`](../../R/utility.R), as today. No warning.
2. Record: `estimateSigmaFromLinearModel(data, route)` returns a plain number and, where no regression was defined,
   first signals with `signalCondition` a condition of class `dbartsSigmaFallback` (not a message class, or a caller's
   handler muffling messages finds no restart: h-record.out). `bart`, `bartBT` and `rbart_vi` (around its validation
   sampler) set `sigest.fallback` on the fit from a calling handler around sampler creation; xbart reads the class
   where it reads the warning's today.
3. Verbose, the fallback: a message of class `dbartsStartingSigmaMessage` in
   [`announceAutoFamily`](../../R/utility.R)'s form, `starting sigma is the sd of the response (<value>): the linear
   model on <p> columns (factors as indicators) leaves no residual degrees of freedom in <n> rows; supply 'sigest' to
   set it`, sent by [`estimateStartingSigma`](../../R/spec.R) under the verbose flag its caller hands it: the
   control's in `resolveSamplerSpec`, xbart's own for its all-rows estimate.
4. Verbose, before every estimate, by the same class and flag: `estimating the starting sigma by a <dense|sparse>
   linear regression on <n> rows and <p> columns; see 'sigest'`.
5. Where the two lines print: at `bart`, `dbarts`, `dbartsSpec` and `xbart` (for all rows, never per fold) as
   messages, as the `family = "auto"` line does; at `bartBT` on standard output, its handler printing the text with
   `cat` and muffling the message; at `rbart_vi` once, from its validation sampler ([`rbart_vi_fit`](../../R/rbart.R)
   muffles the class for the chains), and so not at all on several threads.
6. Summary: [`printSummaryBartBody`](../../R/diagnostics.R) prints `(Starting sigma: the sd of the response; the
   linear model had no residual degrees of freedom)` on `sigest.fallback` (hurdle: positive part).

F. xbart.
1. Route (dec-B382): the all-rows estimate runs once, raises any refusal, sends E's lines and fixes the routine: none,
   `lm.fit` or D. Each fold runs it on its rows, the design built per chunk in [`xbartRunUnits`](../../R/xbart.R), or
   takes their sd under none or if they leave no residual degrees of freedom. Where all rows had a regression and
   some units' rows have none (dec-B429's "wherever"), each chunk returns their count and one line of E.3's class
   is sent for the run under verbose: `starting sigma is the sd of the response in <k> of <units> (replication,
   fold) units: the linear model leaves no residual degrees of freedom in their training rows; supply 'sigest' to
   set it`. A quiet run and the draws are as without it.
2. `sigest` may be a function, on `loss`'s conventions: `function(x, y, weights, offset)` of exactly four arguments
   (else `supplied sigest function must take exactly four arguments`), or `list(function, environment)`, called
   positionally once per unit on the worker that runs it. `x` is the unit's training rows of A's design for "auto" (a
   numeric matrix, or a dgCMatrix with imputed rows written in), built on first use (`delayedAssign`); `y`, `weights`
   and `offset` are those rows', NULL where the call has none. Any return but one positive finite number stops the run
   with `'sigest' function must return one positive finite number`. A binary family does not call it; a fixed residual
   prior refuses it as it does a differing number. There is no all-rows estimate and no E line.
3. It reaches a worker as `loss` does ([`xbartLossFunction`](../../R/xbart.R)): serialized with its environment, so it
   names packages with `::` and finds nothing of the global environment on a socket worker. It is not seeded: one that
   draws random numbers draws on the worker, as the help says of `loss`.

## Change

By file: R/tombstones.R and the doors; R/utility.R (A to D, `sparseResidualStandardError` keeping its name as C and D
and returning NA where no regression is defined, `floorSigmaEstimate` losing its warning); R/spec.R
([`nonFinitePredictorNames`](../../R/spec.R) reads sparse sources and runs before the estimate; E.2 to E.4); R/bart.R,
R/rbart.R, R/diagnostics.R (E); R/xbart.R (F). Both warning classes go: never released, named by no consumer's branch
(carried). Constraints: no engine, bridge, C API, state or class change; no new dependency; no size constant; a dense
design with no factor under a chisq prior bit for bit unchanged. Out of scope: a string "sd", a summary method for
`rbart_vi`, one estimate shared by its samplers, a function at another door.

## Tests

Edits: test-starting-sigma.R (three pins on a factor's codes take the indicator design); test-data-mixed.R and
[test-indicator-storage.R](../../inst/tinytest/test-indicator-storage.R) (the sparse-fallback blocks: no warning, the
dense `sigest` within 1e-10; the latter's direct calls to `sparseResidualStandardError` stay at 1e-8, except the loop
under ["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R), whose near-copies D drops: it runs
under "auto", `expect_identical` to the dense fit); test-boundary-inputs.R and test-xbart-fold-oracle.R (no warning,
`data@sigma` as pinned, the message under verbose); test-sampler-splitProbabilities.R (a seeded comparison of split
counts on a factor design fails: report the counts over five seeds; move the seed if it holds, else stop; likewise any
seeded expectation on a caller-sparse fit or a default fit on a factor frame).

New, in test-sparse-starting-sigma.R, against `lm.fit` on the dense indicator design within 1e-10 unless said;
"traced" counts the calls of `lm.fit` and of `sparseResidualStandardError`.
- dec-B422: `summary(lm(y ~ f + x))$sigma` for a 5-level factor and a numeric column: by default, under
  `factors = "indicators"` and as codes, each also "sparse"; also ordered, two-level, 40-level and incomplete factors.
- dec-B434, traced, by `dbarts`, `bart` and `xbart`: "auto" runs `lm.fit` alone for a matrix, a default frame and
  `factors = "indicators"` with a 6-level factor under `options(dbarts.sparseIndicators =)` "sparse" and "dense"
  (values `identical`), and D alone, with no warning (counted, Gate hygiene), for a bare dgCMatrix, a sparseVector
  column and a sparseFactor column; "dense" on those three runs `lm.fit` alone, `identical` to the dense twin, NA
  entries included; "sparse" on a matrix runs D alone.
- The doors: each rule of Design's list at each door it names, the refusals by their words; `matrixAvailable` stubbed
  FALSE, where "auto" and "dense" on a factor frame still run; "sd" stops under a binary response too; "AUTO",
  "Dense" and "SPARSE" are refused by name; `list(1.5)` is 1.5 at `dbarts`, `bart` and `dbartsSpec`; NULL at
  `bartBT` is "auto"; a failed estimate keeps the generic words and appends the cause.
- No size constant, at the door: a dgCMatrix of 2,100 rows and 2,004 columns (six stored) and its dense six
  columns equal `lm`'s value; 2,050 levels in 2,100 rows, as a factor under "sparse" and as a sparseFactor under
  "auto", equal the within-level regression's.
- D.2 with nothing eliminated: a dgCMatrix of 400 complete one-hot columns and a numeric one at n 40,000 has
  rank 401 and the within-level regression's sigma (rank 402 before the fix round); with a second complete set of
  7, rank 407.
- c-named.R's 63 cases, "sparse", its folds through `xbart`: rows dropped beside a missing factor value (one-row
  levels at weight 0 and at a missing response, five folds, two such factors, a two-level factor with a level held out
  whole, the factor missing on every kept row); C.5 on proportional and duplicated columns, alone, beside their rows'
  indicator in both orders and beside a factor, on nested factors with one-child levels in both orders and on a
  duplicated factor; C.4 and C.5 beside D.1 (within 1e-8); D.3 at rank below n and at full rank (the fallback).
- A.2 stores no missing value: with a 500-level factor and 40 missing at n 2000, its columns hold at most n entries. A
  sparseFactor with a first and a middle reference level; a level per row; `chol` sees only the columns left. D.2's
  rule (c-firstpivot.R): a weighted factor alone has `lm.fit`'s rank; a 4-level factor beside a value near 1.7e9 with
  a spread of 3600 on two levels' rows equals `lm.fit` on the factor alone (1.2684; 0.4755 without the rule).
- f-premises5.R's cases, "sparse": ten explicit zeros in a time column and in its indicator; start and end times
  sharing a missing pattern; a value on a first level's rows; dec-B421's band edges, a near-copy at 1e-4 of its norm
  equal to `lm.fit` (0.5085), one at 1e-6 equal to `lm.fit` without the copy (0.7933), "dense" equal to `lm.fit` on
  both. Also a column scaled by 1e160, and a partner matching in count and index sum but not in rows (no change).
- dec-B429, dense and sparse, by `dbarts`, `bart`, `bartBT`, `rbart_vi` (two chains on one thread; on two at home
  only) and `xbart`: zero warnings; with verbose TRUE E.3's line once and E.4's once per call, where and how E.5 says,
  with no message condition at `bartBT`; with verbose FALSE nothing, and a fit under a caller's
  `message = function(m) invokeRestart("muffleMessage")` runs. `sigest` is `sd(y)` with no attributes, also in a
  `pdbart` result; `sigest.fallback` is TRUE whatever verbose, on an `rbart_vi` fit too, and the summary's line
  prints; both absent for a linear estimate, a supplied sigest and a fit with the field removed; none of it under
  a fixed prior.
- An Inf entry in a dgCMatrix and in a sparseVector column: the dense path's error, no routine called.
- xbart: every fold runs the all-rows routine (traced, under each string) and gets `lm.fit`'s sigma on its rows; an
  all-rows fallback fits no fold. 30 rows, 24 columns and 5 folds, on one thread and on two: the all-rows line and
  F.1's line once each under verbose, nothing when quiet, the losses those of `function(x, y, weights, offset)
  sd(y)` either way. A function, and a list with an environment: each unit's sigma is its value on that
  unit's training rows (the four arguments recorded, `weights` and `offset` NULL where none); one ignoring `x` builds
  no design; one drawing `sample` gives identical losses in two runs after one `set.seed` at `n.threads` 1, and with
  `seed =` leaves `.Random.seed` as it was; a return of NA, 0, Inf, a character or length 2 stops; `function(y)` and a
  function beside a fixed prior are refused; a binary response does not call it.

Reviewer's mutants, each of which must fail a test (the routine's, run on the prototype: c-mutants.out,
c-mutants-named.out): the expansion removed for a container, a matrix of codes, a fold; A.2 without a first level's
indicator, or storing imputed rows; "auto" sending an R-built sparse expansion to D or a dgCMatrix to `lm.fit`;
"dense" or "sparse" ignored; an unknown string read as "auto"; a function accepted at `dbarts`; the Matrix check
removed; explicit zeros kept; C.3 (a) removed; held columns taking part in C.4 and C.5; C.4 without its partner test
or row confirmation; C.5 or its drop removed; tolerance 1e-16 and 1e-6; D.1's block left out of the rank, its u v'
term or the other factors' imputed rows dropped; D.2's rule removed; the intercept back among the unit-normed columns
of D.2's crossproduct; a size cutoff on n or on min(n, p + 1); a rule's name read whatever its case; the Inf check
run late; the condition signaled
only under verbose, or as a message; the message sent with the flag FALSE, once per chain in `rbart_vi`, or as a
message at `bartBT`; the function called on all rows, its return unchecked, `x` built unused, a worker seeded.

## Baselines and gates

- What moves (carried, scratch/sspnew/eq.out, and re-derived for dec-B434, not rerun in the harness): twelve of
  [equivalence.R](../../benchmarks/R/equivalence.R)'s 55 scenarios. POSTERIOR-CHANGING: sparse, mixedmatrix,
  sparsefactor, testswap, leaffactormixed, factorpartial, xbartmixed (a caller's sparse source: sd(y) to D's value)
  and categorical, leaffactor, nafactor, ordfactor (a dense frame's factor: codes to `lm.fit` on indicators, 0.9-34's
  value). SHIFTING: wideFactorIndicators (the sparse QR to `lm.fit`, 6.6e-16). NEUTRAL: the other 43, bcf's 15,
  multinomial's 11, the snapshot files, every exact gate. Another mover is a stop.
- Re-record the twelve on the reference build into a copy of equivalence-e4faed5c
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)). Partition: 43 of 55 identical, the movers the twelve, no |z|
  above 4 on wideFactorIndicators; 55 of 55 bitwise from a second `--preclean` install.
- Oracle (MANIFEST rule P17): agreement with `lm.fit` on the indicator design: c-sweep.R's 300 designs and c-named.R's
  cases on the package's routine; a sweep over 400 random designs (n 30 to 2000, p 5 to 1530, dependent columns, zero
  weights, factors with missing values) and 30 tall ones (n 1e5 to 3e6, one to three factors of 6 to 5,000 levels).
- Independently of the implementer: the posterior-changing battery of
  [RNG classes and their gates](README.md#rng-classes-and-their-gates), `R CMD check --as-cran`, the lint set.
- Speed and memory: `estimateSigmaFromLinearModel` under "sparse", timed end to end from the data object, one
  configuration per process, on e-cost.R's frames: n 2e5 with a 50,000-level factor, complete and with 2 percent
  missing, at most 2 s and 400 MB (run: 0.8 s, 240 MB); two crossed 5,000-level factors with 2 percent missing in
  each, at most 25 s and 1.5 GB (16 s, 980 MB); a dgCMatrix of n 2600, p 1999, at most 2 s under R's own BLAS.

## Help and docs

- `sigest` in the five man pages that take it: a number or one of the three strings, in Goal's terms. `"sparse"` needs
  Matrix; it agrees with `lm` to rounding except on nearly dependent columns (a near-copy of another column, a tiny
  spread beside a large mean), where it may drop a column `lm` keeps; its last digit follows the BLAS and thread count
  (dec-B430); it reads a sparseFactor, not a caller's own one-hot columns, as a factor. Pass `"sparse"` for a data
  frame with a factor of thousands of levels, `"dense"` for `lm`'s handling of nearly dependent columns or for a
  sparse fit's draws to match its dense twin's. No size turns the estimate off and no timing is quoted: `"sparse"`
  factors a square matrix as wide as the columns left beside the widest factor and can fail to allocate it, `"dense"`
  holds rows times columns doubles, and neither can be interrupted (Open call 1). For an expensive estimate pass a
  number, by the code of scratch/sspfold4/d-examples.R (carried): `sd(y)`, too high where there is linear signal; a
  regression on a random subset of rows; a cross-validated lasso (glmnet) for a wide or sparse x.
- man/xbart.Rd: `sigest` may be a function, pointing at `loss` and `seed` for its calling form, workers and random
  numbers (F.2, F.3), with the alternatives above as functions. man/bart.Rd: both warning classes go, the value lists
  `sigest.fallback`. man/rbart.Rd: the estimate is made once per chain and once more.
- inst/NEWS.Rd: one item each for the strings, the fallback's line and xbart's function. docs/design: a dated
  paragraph in sparse-columns.md's [R surface](../design/sparse-columns.md#r-surface) and a dated line in
  starting-sigma-sensitivity.md's [Recommendation](../design/starting-sigma-sensitivity.md#recommendation);
  error-style.md and memory-footprint.md follow. At landing: TODO, the ledger entry, Status, the MANIFEST row.

## Steps

1. Design A to F with the tinytest edits and new tests; the suite green against `R CMD INSTALL -l <lib> .`.
2. The sweeps; the speed and memory points.
3. Help, NEWS, docs.
4. After review: the re-record, in its own commit.

Stop when: the diff passes ~1300 lines; a scenario other than the twelve, a snapshot or an exact gate moves, or
wideFactorIndicators shows a |z| above 4; a sweep shows a rank below `lm.fit`'s, a rank above it on a frame's factors
or at n under 1e6, or a sigma off by more than 1e-10 plus `(rank excess) / (n - r)` above 30 residual degrees of
freedom; a speed or memory limit is passed; a tinytest fails that Tests does not name; the change needs engine code.

## Calls made

- D.1, C.3 (a), C.5 with its drop, and D.2's rule are in no ruling's words (dec-B403 leaves the how to measurement).
  D.1 serves one wide factor, the case the help sends to "sparse"; C.5 takes start and end times sharing a missing
  pattern from 1.760 to `lm.fit`'s 0.4995 (run); the rest repair a wrong rank (c-mutants.out, c-firstpivot.out).
- A factor's imputed rows are never stored. Rejected: storing them and gating the cost (19 s and 4.0 GB for 0.3 s and
  220 MB at n 1e5, 20,000 levels and 2 percent missing; carried, run); storing those of the factors not eliminated (31
  s for 15 s at 0.1 percent missing in two crossed 5,000-level factors; run, e-cost-stored.out).
- A.2 emits every present level, so a value on a first level's rows finds its partner (0.5002 for 1.005; run). D.3 has
  no row equilibration: none of 60 wide designs, weights over 1e14, needs it (run).
- NULL stays the formals' default (rejected: "auto" in six signatures, moving codoc and every test of an absent
  sigest). A numeric string stays a number, as in 0.9-34 (rejected: refusing "1.5" as `density` does an unknown `bw`).
- A function is taken by xbart only, the door dec-B433 names; elsewhere `bart(x, y, sigest = f(x, y))` does the same.
  `dbartsSpec` estimates over a carried sigma under "dense" or "sparse", as a number replaces it (rejected: keeping
  it, the argument ignored without a word; refusing).
- A data frame or list holding a sparse column is "the caller passed sparse" (dec-B434 names "a caller's sparse
  matrix" and "a data frame"; rejected: dense for every data frame). A string where `sigest` has no use follows the
  number's handling (rejected: silence for "auto" at `bart`, which a written-out NULL does not get).
- The record is a signaled condition and a field of the fit. Rejected: an attribute on the sigma value, which surfaces
  in `pdbart`'s and `rbart_vi`'s `sigest`, `rbart_vi`'s state, a one-draw prior-predictive sigma (h-record.out) and
  any consumer reading the slot; a logical slot, a class change and stale after a direct write.
- xbart's function follows `loss`: positional, a fixed argument count, the list form, no seeding. Rejected: named
  arguments and a seed per unit, against the help's "draws on the worker running it" and `xbartRunUnits`'s rule that
  no worker calls `set.seed`, and varying with the worker's `RNGkind` (carried, seed.out).
- E.4's line and where E's lines print are the orchestrator's call, not a ruling: where each door's verbose lines
  print today, nothing where a door is silent. `bart` with a named family gets its first message.

## Open calls

Both are ruled. Call 1: (a), the estimate is left uninterruptible and the help says so (dec-B438). Call 2: (b), no
estimate is made under a fixed residual standard deviation, whose value the fit reports as `sigest` (dec-B439).

### 1. Can a user stop the starting-sigma estimate once it has begun?

Background. Before a fit starts, dbarts estimates the residual standard deviation by a linear regression. Since
dec-B433 nothing limits its size, and the help is to name alternatives "in case a user has to cancel at the initial
estimate phase". A user cannot cancel it. Measured on an arm64 Mac, R 4.6.1, under R's own BLAS and under Apple's
Accelerate: Ctrl-C sent 2 seconds in, and separately a 2 second `setTimeLimit`, stopped none of the estimate's three
long steps. Under R's BLAS `lm.fit` (the default for a matrix or data frame; 4,000 rows by 2,500 columns) ran its full
14 seconds, the sparse crossproduct (400,000 rows by 3,000 columns) its 16 and LAPACK's Cholesky (6,000 columns) its
23; under Accelerate, 7, 16 and 7 seconds (the Cholesky at 13,000 columns). R acts on the interrupt only when the
compiled call returns; until then the one way out is to kill R, losing the session. `lm` itself behaves so, and so
does 0.9-34, whose estimate on the same matrix took Ctrl-C after 14 seconds. What 1.0 adds is the sparse estimate,
whose Cholesky takes 117 seconds at 10,000 columns under R's BLAS and 3.5 seconds under Accelerate.
- (a) Leave it. Nothing to build. The help says it cannot be interrupted, and a verbose fit prints what it will
  compute, on how many rows and columns, before starting. Cost: a user in a long estimate waits or loses the session.
- (b) Make the sparse estimate interruptible by running its crossproduct and Cholesky in blocks from R code: about 40
  lines, 25 of tests, and a poll for time limits not yet written (as prototyped a time limit stopped the Cholesky in
  0.2 seconds and ran through the crossproduct's 18 to 19). Ctrl-C then stops either within 0.1 to 1.5 seconds under
  both BLAS libraries. Cost: the R Cholesky takes the same time under either, 0.6 seconds at 2,000 columns, 7 at 6,000
  and 28 at 10,000, so 4 times faster than LAPACK under R's BLAS (117 seconds at 10,000) and 8 times slower under
  Accelerate (3.5 seconds); it needed 1.3 GB of working memory for LAPACK's 0.3 GB at 6,000 columns; it replaces the
  routine dec-B430 kept; and `lm.fit`, the default for a matrix or data frame, stays uninterruptible.
- (c) Run the estimate in a forked child process that Ctrl-C kills: about 12 lines, the same numbers, dense and sparse
  alike, in 23.3 seconds against 23.1 under R's BLAS and 8.1 against 7.5 under Accelerate, stopped within 0.1 to 1.2
  seconds. Cost: no fork on Windows, and none inside RStudio, Positron or R.app, where xbart already declines to fork;
  elsewhere every fit forks, and a fork after a multithreaded BLAS has started can hang, as xbart's help warns.

Recommended: (a). A user needs the session back, and no option gives that for the default estimate where most users
work: (b) covers only sparse input and slows it 8 times under a fast BLAS, (c) covers neither Windows nor the editors.
What every user can be given is the warning before the wait, which (a) has.

### 2. Under a fixed residual standard deviation, is the estimate made at all?

Background. With `family = gaussian(sigma = fixed(value))` sigma is not sampled, so the starting estimate calibrates
nothing. dbarts makes it anyway. Run on 80 rows with the standard deviation fixed at 0.7: the fit reports `sigest`
0.4734, lm's value, and a predictor holding Inf stops the fit with "unable to obtain a starting estimate of sigma".
0.9-34 did the same through `resid.prior = fixed()`: 0.4734 in the slot, the Inf refused (run). With no size limit
left, that unused estimate can take minutes, and a user cannot skip it with a number: dec-B413 accepts one beside a
fixed value only where the two agree, and refuses any from 1.1-0.
- (a) As today (the steps build it): the estimate is made by `sigest`'s rule and only reported, at full cost.
- (b) Make no estimate under a fixed value: about 10 lines. A user sees three changes, the first two also against
  0.9-34: `sigest` in the fit is the fixed standard deviation (0.7 for 0.4734 above); a predictor holding Inf no
  longer stops the fit there (run: accepted); "dense" or "sparse" beside a fixed value does nothing. Draws are
  identical (run, two chains of 50).

Recommended: (b). A user who fixed sigma gets nothing from the estimate and can wait minutes for it.

## Landing note

Landed 2026-10-10 as 8c27f225..d820097c: one opus implementer, one opus review that ran the plan's mutants and 49 of
its own, a fix round, the same reviewer's check of it with 27 more mutants, and two x86 runs of tests/cpp and tinytest.
The twelve scenarios the plan names moved and no other; they are re-recorded as equivalence-2b48939b, and the
MANIFEST row gives the old and new starting sigma of ten of them. What was built is the Design, with these departures
and additions (dec-A200):

- The step that eliminates no factor is not the plan's. The review found a caller's complete set of one-hot columns
  given one rank too many from about 400 levels (31 of 120 designs; sigma off by at most 2.5e-5): after unit norm the
  intercept tied the indicators at the first pivot. The crossproduct is now taken before unit norm, and the intercept
  stays out of the factorization and joins by what the other columns leave of it (0 of 120, and 0 of 216 adversarial
  designs in the review's check).
- Two limits are accepted. Two complete sets of a caller's one-hot columns under non-integer weights gave one rank too
  many in 1 of 30 designs at 6e6 rows in the fix round, and in 0 of 155 in the review's check. Where a factor is
  eliminated under weights of few distinct values the intercept's leftover grows with the number of rows and would
  cross the tolerance near 1.6e7 (rank right in 30 of 30 up to 1.2e7).
- Under dec-B429's "wherever": a verbose xbart run prints one line counting the (replication, fold) units that started
  from the response's sd where the all-rows design allowed a regression, and an rbart_vi fit records `sigest.fallback`.
- `bartBT(sigest = NULL)` is "auto"; `dbarts(sigest = list(1.5))` is 1.5; a rule's name is matched in lower case only;
  the estimate's generic failure carries its cause; the rule helpers sit beside the starting-sigma code.
- Two tests the plan did not name moved from the old coercion message to the refusal of an unknown rule.
- Size: about 2,200 lines added beside this plan, 1,140 of them tests, against a stop the orchestrator lifted to about
  2,000.

Left open, as before the slice, and in TODO: the fallback's sd counts rows of weight 0
(starting-sigma-fallback-zero-weights), and an integer column's NA in an indicator expansion
(integer-na-in-model-matrix). Not verified: an older Matrix, a threaded BLAS, Windows, more than 1.2e7 rows.
