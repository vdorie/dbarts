# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403); revised after its blind critique; revised 2026-10-10 for dec-B421,
dec-B422 and dec-B426 and for the sigest-beside-fixed landing (dec-B413). Not built.

agent: blind critique of this revision first (opus; it adds a routine and widens what moves); one sonnet
implementer (R only, no engine code); one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING, for a continuous response under a chisq residual prior with no sigest, in two
cases: a caller's sparse x (dec-B403), and any fit that keeps a factor as a categorical predictor, the
default, dense or sparse (dec-B422). Rounding only for an indicator expansion R built sparse (dec-B370).
NEUTRAL for a dense design with no categorical column, every fixed-unit family, every fit given sigest, and
the draws of a fit under a fixed residual prior.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~1000 lines (R/utility.R ~260 net of the QR's removal, R/spec.R ~35, R/xbart.R ~35, tinytest ~380,
benchmarks/R ~160, man ~40, docs ~90, NEWS ~2, MANIFEST one row).

## Goal

A sparse x gets the starting sigma its dense equivalent gets, with no warning for being sparse: from an exact
routine while the smaller of its row and column counts is at most 2,000, and from LSQR above that. Every
factor, in a dense design too, enters that regression as indicator columns. One routine serves every sparse
design and xbart's per-fold estimate; the sparse QR and the class dbartsSparseSigmaFallbackWarning go. An
infinite entry in a sparse x is refused as in a dense one, and a fixed residual prior makes no estimate. The
tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Rulings, each covering what it names:

- dec-B403: "Yes, just use the sparse regression, no warning." A sparse x takes its dense equivalent's
  linear-model starting sigma; a caller can pass sigest and the help says so; how it is computed is measured.
- dec-B421: "Band is fine." The exact routine's rank tolerance (below) may differ from lm's on near-copies.
- dec-B422: "What would keeping the codes mean? Treating them as integers? That doesn't make sense. I think
  indicators are necessary, right?" Every factor, unordered or ordered, enters the starting-sigma regression
  as indicator columns; the trees' representation is unchanged; the slice re-records what moves.
- dec-B426: "I'm ok with C, but do we want to evaluate it in more contexts?" then "Go ahead." The exact
  routine up to m 2,000 (m the smaller of n and p + 1); above it LSQR with structural degrees of freedom and
  no message; above the cutoff a design with under 10 percent residual degrees of freedom takes sd(y). Two
  conditions of the build ([Build conditions](#build-conditions)). The measurements are in
  [Part 2: LSQR](../design/starting-sigma-sensitivity.md#part-2-lsqr).

Today, run on this tree (2026-10-10, shipped build, arm64 macOS, R 4.6.1, Matrix 1.7-5):

- [`estimateSigmaFromLinearModel`](../../R/utility.R) returns sd(y - offset) with the class for any source
  with a CSC column: a bare dgCMatrix (1.053 where lm gives 0.500), a frame with a sparseVector column, a
  frame with a sparseFactor. An indicator expansion R built sparse (`sparse.from.indicators`, set by
  [`makeIndicatorModelMatrix`](../../R/utility.R)) goes to [`sparseResidualStandardError`](../../R/utility.R),
  Matrix's sparse QR; anything else to [`residualStandardError`](../../R/utility.R) (`lm.fit`).
  [`floorSigmaEstimate`](../../R/utility.R) turns a non-finite result into sd(y) with
  dbartsSigmaFallbackWarning: a dense design with p >= n warns; one with 2 residual degrees of freedom of
  40 rows gives an estimate and no warning.
- A factor kept categorical enters as its level codes, ordered or not
  ([`as.matrix.dbartsMixedMatrix`](../../R/mixedMatrix.R)): on a 40-level factor and one numeric column
  (n 300, noise sd 0.5) the default starts sigma at 1.088 (sd(y) 1.095) and `factors = "indicators"` at 0.499.
  An ordered factor: 0.534 on its codes, 0.489 as indicators, which is `lm`'s value under its polynomial
  contrasts.
- The frame rebuilt from a dense container's column list gives, through
  [`makeIndicatorModelMatrix`](../../R/utility.R), the same design values as the `factors = "indicators"`
  fit of the original frame. That builder stores an indicator column sparse where its density is at most
  [`sparseIndicatorDensity`](../../R/utility.R), so a 5-level factor at n 300 already takes the sparse
  routine under "indicators". A missing factor value is NA in every indicator of its block and is imputed
  with the level's frequency, so the block still sums to one.
- An infinite stored entry in a dgCMatrix is accepted (sd(y), the sparse class); a dense one is refused,
  naming the column ([`nonFinitePredictorNames`](../../R/spec.R)).
- A fixed residual prior with no sigest still runs the estimate ([`resolveSamplerSpec`](../../R/spec.R);
  [`xbart`](../../R/xbart.R) before it resolves the prior): the slot and `bart`'s `sigest` hold the linear
  estimate (0.293 under `fixed(4)`), and a sparse x warns about an estimate nothing uses. Draws are identical
  whether the slot holds that estimate or the fixed sigma. Since dec-B413
  ([`refuseSigestUnderFixedPrior`](../../R/family.R)) a sigest within 4 machine epsilons of the fixed sigma
  is accepted with a once-per-session message, then sits in the slot with no estimate made; a differing one
  is refused (dec-B425). The slot is read later: a `setModel` to a chisq prior, `rbart_vi`'s start and
  `bart`'s `sigest`, so it is never left NA.
- xbart on a sparse x takes the per-fold route "marginal" (each fold its own sd); a dense design "linear".
  A dense design whose folds have no residual degrees of freedom raises dbartsSigmaFallbackWarning.

The exact routine, measured 2026-10-09 (R's reference BLAS; the first prototype's seconds, the routine of
[Algorithm](#algorithm) D timing within 10 percent of it at four of these points):

| n | p | design | dense lm.fit | sparse QR (today) | exact routine |
|---|---|---|---|---|---|
| 1e3 | 1000 | 1% | 0.49 s | 0.28 | 0.14 |
| 1e4 | 100 | 1% / one-hot 10 x 10 | 0.09 / 0.09 | 0.03 / 0.86 | 0.04 / 0.04 |
| 1e4 | 1000 | 1% / one-hot 20 x 50 | 7.8 / 7.8 | 7.4 / 157 | 0.18 / 0.20 |
| 2e4 | 2000 | 1% | 64 | 65 | 1.2 |
| 1e5 | 1000 | 0.1-5% / one-hot 20 x 50 | 83 | not run | 0.16-1.6 / 0.45 |

- Its cost is the crossproduct plus a pivoted Cholesky of an m x m matrix: 1.0 s and about 100 MB over R's
  own 250 MB at m 2,000, the largest it now runs at (14 s at 5,000 and 133 s at 10,000 are why it stops
  there). Where both take under a second it can be the slower (0.03 s against 0.01 s at n 1e3, p 100).
- Agreement with `lm.fit` on the dense design: 400 random designs (n 30 to 2000, p 5 to 1530, numeric,
  one-hot and mixed, a third with dependent columns, a third weighted with zero weights, a third with an
  offset; 143 with p >= n). Rank equal to lm.fit's on 398; of the 292 with an estimate, sigma within 1e-11
  relative above 30 residual degrees of freedom and 3e-8 at 30 or fewer. On the other 2 (wide), lm.fit keeps
  columns with relative pivots near 1e-16, the routine's rank equals the SVD's, and sigma is 0.26 and 0.64
  percent apart.
- The band (dec-B421). A column is dropped when its residual after the kept columns is below 1e-5 of its norm
  (pivot 1e-10 on the unit-diagonal crossproduct; lm's is 1e-7). An exact dependency's pivot measured 0 to
  1.5e-13, so 1e-10 holds a margin of about 700 and lm's 1e-14 none. Two near-copies differing by 1e-7 to
  1e-5 of their norm both stay in lm's fit and one stays here (0.01 to 4.7 percent apart, equal to the dense
  fit without one copy); a dense column whose spread is below 1e-7 of its mean is dropped by lm and kept here.

LSQR, from the prototype behind dec-B426 (base R and Matrix; rerun here): two runs on one design return the
same bits; weighted fits (uniform, 20 percent zeros, fold-like 0/1, skewed, a few at 1e-12, with an offset)
agree with `lm.wfit` to 7e-11 in 11 to 50 iterations at n 6000, p 400, but 1 percent of rows at weight 1e6
takes 785 iterations (2.7e-6 apart), near the cap; a constant, fully stored column under weights survives
the prototype's sum-of-squares test and costs a degree of freedom; a full one-hot block left uncounted costs
one (9e-5), counted or built without one level it agrees to 3e-12; a zero response divides by zero. A row
subset of a 1.5e6-entry design costs 2.8 iterations' products. Costs:
[Cost and memory](../design/starting-sigma-sensitivity.md#cost-and-memory).

What moves (every call of [`estimateSigmaFromLinearModel`](../../R/utility.R) logged through a load hook,
2026-10-10, beside the value `lm.fit` gives on the indicator design):

- equivalence.R, quick mode, 132 calls in 44 of 55 scenarios. Twelve move:
  - seven on a caller's sparse source, sd(y) to the linear estimate:
    ["sparse <- list("](../../benchmarks/R/equivalence.R) 2.959 to 1.370,
    ["mixedmatrix <- list("](../../benchmarks/R/equivalence.R) 1.670 to 0.271,
    ["sparsefactor <- list("](../../benchmarks/R/equivalence.R) 3.609 to 2.122,
    ["testswap <- list("](../../benchmarks/R/equivalence.R) 3.664 to 2.265,
    ["leaffactormixed <- list("](../../benchmarks/R/equivalence.R) 3.941 to 1.833,
    ["factorpartial <- list("](../../benchmarks/R/equivalence.R) 3.483 to 2.125 and
    ["xbartmixed <- list("](../../benchmarks/R/equivalence.R) 4.124 to 2.386 (all rows);
  - four on a dense frame with a factor, codes to indicators:
    ["categorical <- list("](../../benchmarks/R/equivalence.R) 1.9855 to 1.9870,
    ["leaffactor <- list("](../../benchmarks/R/equivalence.R) 1.863 to 1.830,
    ["nafactor <- list("](../../benchmarks/R/equivalence.R) 2.469 to 2.389 and
    ["ordfactor <- list("](../../benchmarks/R/equivalence.R) 2.234 to 2.330;
  - ["wideFactorIndicators <- list("](../../benchmarks/R/equivalence.R), by rounding (the QR's
    2.6995748335729584 against lm.fit's ...566).
  The other 32 with a call are dense with no categorical column; 11 make none.
- bcf-equivalence.R: 12 calls, all dense, no factor. multinomial-equivalence.R: none.
- The four test-reproducibility files (build guard bypassed): 5 calls in three files, none sparse and none
  with a factor; the binary file stopped under the bypass with nothing logged (a binary family makes no
  estimate).
- exact-gates.yaml's list in quick mode, bcf-latent-exact.R aside: 35 calls, 30 dense with no factor and 5 on
  a factor design (mask-redraw-exact.R 4, change-balance.R 1), all five under a fixed residual prior (read
  at the call), which makes no estimate after Change 4. No exact gate moves.
- tinytest (19709 results, 0 failures): 4 files reach the caller-sparse branch, 3 the QR, 37 a dense
  container with a factor (409 calls). With the new values substituted at the hook, six expectations fail in
  four files: test-starting-sigma.R 3, test-indicator-storage.R 1, test-data-mixed.R 1,
  test-sampler-splitProbabilities.R 1 ([Tests](#tests)).

## Algorithm

`sparseResidualStandardError(y, x, weights, offset)` keeps its name; the cutoff (2000) and LSQR's cap (1000)
are formals with those defaults, so a test reaches either route on a small design.

A. The design, built by one new function that every caller uses (`startingSigmaDesign(x)`: the all-rows
estimate and xbart's chunk runner). The expansion of dec-B422 happens here and nowhere else; `data@x`, the
cut grid and the trees never see it.

1. A plain matrix: [`sigmaDesignMatrix`](../../R/utility.R) as today, for `lm.fit`.
2. A dense container (no CSC column) with a factor column: the frame is rebuilt from its column list and
   names and handed to [`makeIndicatorModelMatrix`](../../R/utility.R) with `drop = TRUE` and storage
   "auto", so it is the design `factors = "indicators"` builds and the two fits' sigest are identical. The
   result is a matrix (case 1) or a container with sparse-built indicators (case 4). Without a factor
   column: as today.
3. A caller's sparse source. A bare dgCMatrix is wrapped ([`wrapSparseTestMatrix`](../../R/mixedMatrix.R)).
   In a container, a dense-backed factor column becomes one indicator per present level but its first; a
   sparseFactor column (non-NA `sparseReference`) one indicator per present level but its reference, from
   its stored entries, a stored entry at the reference code dropped. A stored NA is an NA entry in each of
   the factor's indicators. Other columns as in [`sparseDesignMatrix`](../../R/utility.R) today (dense
   numerics to CSC), then its imputation (the mean of the observed entries, implicit zeros counted: for an
   indicator the level's frequency). No indicator of the reference level is built, so the columns span with
   the intercept what the dense equivalent's full block spans and no block is exactly dependent.
4. An indicators-route container (built by R, or beside a caller's sparse column): taken as it is.
   [`makeIndicatorModelMatrix`](../../R/utility.R) records per column the input term that emitted it
   (`indicator.term`, in place of `sparse.from.indicators`, whose one reader goes).

B. Front end, shared by both routines, in this order.

1. Inf: an infinite entry is an error that [`estimateStartingSigma`](../../R/spec.R) reports as for a dense
   design, naming the column: [`nonFinitePredictorNames`](../../R/spec.R) learns sparse sources.
2. Rows: drop rows with a missing response, weight or offset and rows of weight 0, as `lm.wfit` does;
   z = y - offset; n the rows kept.
3. Constants: a column whose entries over the kept rows are all equal (implicit zeros included; an exact
   comparison of values, not of a sum of squares) is dropped. p is the number of columns left and
   m = min(n, p + 1).
4. Centering: each dense-backed column and each column stored in every kept row is centered at its weighted
   mean over the kept rows. Then every column is divided by its largest absolute entry.
5. Blocks: the columns of one factor term (`indicator.term`) are a full block when at least two are left
   and they sum to one in every kept row (within 1e-8; checked on the design, so a caller's drop pattern
   or a fold's rows cannot miscount). b is the number of full blocks.

C. Route. m <= 2000: the exact routine (D). Otherwise the structural residual degrees of freedom are
df = n - 1 - p + b; when df < 0.1 n there is no estimate (F); else LSQR (E). m is the column count after the
expansion, so a dense frame with a 5,000-level factor can be an LSQR design.

D. The exact routine.

1. B = diag(sqrt(w)) [1 X]; each column divided by its norm. Sparse columns are not centered: a centered
   crossproduct cancels for a large-mean column, which is why the intercept stays a column.
2. Narrow (ncol(B) < n): `chol(as.matrix(Matrix::crossprod(B)), pivot = TRUE, tol = 1e-10)` with
   `suppressWarnings` around that one call (its one warning is the rank deficiency); r its rank, K its first
   r pivots, R its leading r x r block; b = R^-1 R^-T B_K' (z sqrt(w)), e = z sqrt(w) - B_K b, computed
   directly. No iterative refinement (one step moved nothing by more than 1e-12 on 30 near-tolerance designs).
3. Wide: K = `as.matrix(Matrix::tcrossprod(B))`, equilibrated by D = sqrt(diag(K)), the pivoted Cholesky of
   D^-1 K D^-1 at the same tolerance (unequilibrated, a row at relative weight 1e-14 reads as dependent).
   r >= n gives no estimate; else Q an orthonormal basis of D L (L the first r columns of the factor,
   unpivoted) and e = z sqrt(w) - Q Q' z sqrt(w).
4. sigma = sqrt(sum(e^2) / (n - r)); no estimate when n - r <= 0.

A full block is an exact dependency (pivot at most 1.5e-13), so the rank found is lm's on the indicator
design and the band of dec-B421 concerns near-copies only, as before the expansion.

E. LSQR (Paige and Saunders), in R in R/utility.R, on base R and Matrix's sparse products alone. The
design is never copied or made dense past the front end.

1. Operator. With h = sqrt(w), s = sum(w), mu_j the weighted mean of column j and d_j the reciprocal of its
   weighted centered norm, A = diag(h) [1 / sqrt(s), (X - 1 mu') diag(d)], applied as two sparse products a
   step (`X %*% v`, `Matrix::crossprod(X, u)`) with the centering as a rank-one correction. Every column has
   unit norm and the intercept is orthogonal to the rest. Unit column norms are the only preconditioning. A
   column whose centered sum of squares is not above 1e-24 of its uncentered one is dropped as constant
   (weights can make one so) and leaves p.
2. Right side b = h z, start at zero. If the norm of b is 0, sigma is 0; if A'b is 0, the residual is b.
3. Iteration: Golub-Kahan bidiagonalization with the standard updates, no reorthogonalization, no damping.
4. Stop at the first of: the normal-equations estimate at or below 1e-6 of the product of the running
   operator norm and the residual norm; the residual norm at or below 1e-10 of the norm of b (the fit
   reproduces the response, where the first test cannot fire and the recurrence drifts); a zero alpha or
   beta; 1,000 iterations. The tolerance is a constant, not an argument a caller sets.
5. At the cap the iterate is used like any other, with no message: its residual is at least the minimum,
   so the estimate is high, never low. [Build conditions](#build-conditions) (a) measures by how much.
6. sigma = sqrt(sum(r^2) / df), r = b - A x computed from the returned coefficients, not the recurrence.
7. Weights enter as the row scale h in the operator and the right side and in every mean and norm; n counts
   rows of positive weight. This is `lm.wfit`'s fit and `summary.lm`'s degrees of freedom.
8. Degrees of freedom are C's count, never the iteration's. A dependency the structure does not show (a
   near-copy, an interaction block containing its main effects) is counted as a fitted column and
   overestimates by sqrt((df + d) / df), 0.5 percent for 50 in 4550; dec-B426 names structural degrees of
   freedom.
9. Determinism. No random number, fixed constants, one operation order: sums are R's `sum`, the products
   Matrix's sparse kernels, and nothing on this route goes through LAPACK. The same inputs give the same
   bits on one platform and Matrix build (run on the prototype). The exact routine calls LAPACK and so
   holds its bits for one BLAS and thread count, as `lm.fit` does on the dense path.

F. No estimate (D with no residual rank, C under 10 percent) returns NA, which
[`floorSigmaEstimate`](../../R/utility.R) takes to sd(y - offset) with dbartsSigmaFallbackWarning, as the
dense path does. Whether the 10 percent case warns is [Open calls](#open-calls) 1; the steps build the
recommended option, and the other is about 8 lines. An allocation failure is an error, which
[`estimateStartingSigma`](../../R/spec.R) already reports as "unable to obtain a starting estimate of sigma".

## Change

1. R/utility.R: `startingSigmaDesign` (A); [`sparseDesignMatrix`](../../R/utility.R) gains the wrap and the
   two expansions; [`sparseResidualStandardError`](../../R/utility.R) rewritten as B to F with the exact
   routine and LSQR as two internal functions (the batching, the refactor loop and the `grepl` muffler of
   Matrix's warning go); [`estimateSigmaFromLinearModel`](../../R/utility.R) loses the fallback branch and
   its warning and routes by what the design builder returned;
   [`makeIndicatorModelMatrix`](../../R/utility.R) records `indicator.term` and no longer sets
   `sparse.from.indicators`.
2. R/spec.R: [`nonFinitePredictorNames`](../../R/spec.R) for sparse sources (the dense list and the CSC
   block's stored entries, named through the container's map, or positions for a bare dgCMatrix).
3. R/xbart.R: the design is built once per chunk by `startingSigmaDesign` in
   [`xbartRunUnits`](../../R/xbart.R); `foldData` fits its training rows by the same routes, each fold by
   its own m and its own degrees of freedom, a dense design keeping
   [`residualStandardError`](../../R/utility.R). A fold with no estimate falls back as a dense fold does
   today; where the all-rows estimate fell back every fold takes its marginal sd, as today. An indicator
   expansion stops being made dense there.
4. A fixed residual prior makes no estimate. In [`resolveSamplerSpec`](../../R/spec.R), where `residPrior`
   is a dbartsFixedPrior and the family is not on a fixed unit scale, the slot takes the square root of the
   fixed variance, whatever it held; [`xbart`](../../R/xbart.R) resolves its prior before its all-rows
   estimate and does the same, with no per-fold route. [`refuseSigestUnderFixedPrior`](../../R/family.R)
   still runs first and is untouched: a differing sigest is refused, an agreeing one accepted with its
   message, and the slot then holds the same value either way. Draws are unchanged; `bart`'s `sigest` and
   the slot report the fixed sigma, and a later `setModel` to a chisq prior, or `rbart_vi`'s start, reads it.
5. The class dbartsSparseSigmaFallbackWarning is retired: new in 1.0, nothing raises it, and no consumer
   branch names it (stan4bart, bartCause, treatSens, bairrtt).

## Constraints

- No engine, bridge, C API or state change. No new dependency: Matrix (Suggests) as today; without it no
  sparse source exists and every indicator builds dense. Base R calls within R 4.2, Matrix within 1.4-1.
- A dense design with no categorical column under a chisq prior is bit for bit unchanged.
- Out of scope: a cutoff or LSQR for the dense path; the 10 percent rule below the cutoff (Open calls 2); a
  per-fold crossproduct downdate or warm start in xbart (a warm start saved 2 to 8 percent).

## Build conditions

Both are dec-B426's. The rules are fixed here, before any number exists.

(a) LSQR on badly conditioned realistic designs. A tracked script, benchmarks/R/starting-sigma-lsqr.R, calls
the implemented routine; two jobs at a time.

- Designs, each at m 10,000, 20,000 and 50,000 with n = 1.5 m and n = 1.12 m: crossed factors (two factors
  with skewed level frequencies, main effects and their interaction as one-hot columns of a dgCMatrix, so no
  dependency is structural; 100, 140 and 220 levels each); word counts (term frequencies Zipf, document
  lengths lognormal, raw counts); correlated numeric columns (groups of 20 sharing one sparse pattern at 5
  percent, pairwise correlation 0.99, beside 200 fully stored columns correlated at 0.999). At m 10,000 each
  also with skewed weights (a cubed exponential) and with 1 percent of rows at weight 1e6, and with the
  signal on 1 percent of the columns and on all of them.
- Reference: the exact routine at m 10,000 (cutoff raised, one job at a time); above it LSQR at tolerance
  1e-10 and a cap of 20,000, admitted only if it stops on tolerance and if the same setting agrees with the
  exact routine within 1e-6 on that family at m 10,000. The crossed design's rank is its count of non-empty
  cells.
- Recorded per run: iterations, seconds, stop reason, the residual sum of squares against the reference's,
  and sigma against the reference's.
- Large means more than 10 percent above the reference. The note's evidence: sd(y), 15 to 75 percent above
  the linear estimate in its settings, was within seed noise on every measure at n 1000 and 5000
  ([Results](../design/starting-sigma-sensitivity.md#results)), and LSQR runs only where n is above about
  2,200 (m over 2,000 with a tenth of n left over); 10 percent is under the smallest error measured there
  and found invisible.
- PASS: every run stops on tolerance within the cap, or stops at the cap no more than 10 percent above its
  reference. A run more than 1 percent above, from the solver or the count, is named in the Landing note.
- BACK TO THE MAINTAINER, before landing: any run at the cap and more than 10 percent above, with the table
  and three priced options (a higher cap, sd(y) at the cap, a block preconditioner). A run more than 10
  percent above from the degrees-of-freedom count alone goes back the same way. A design with no admitted
  reference goes back as unmeasured.
- DEFECT, stop: any run more than 1e-6 below its reference (an early stop and a structural count can only
  be high).

(b) Weighted fits against the exact routine, in tinytest with LSQR forced (cutoff 0) on designs under the
real cutoff: sparse numeric, an indicators-route design with full blocks, a frame with a sparseFactor; under
uniform weights, 20 percent zero weights, 0/1 fold weights, skewed weights, a few rows at 1e-12, and weights
with an offset. PASS: sigma within 1e-6 relative and the same degrees of freedom on every one. FAIL, stop: a
weighted case past 1e-6 whose unweighted design is within it. The 1e6-weight case is (a)'s, where it is
expected near the cap.

## Tests

Edits the substitution run found:

- [test-starting-sigma.R](../../inst/tinytest/test-starting-sigma.R) pins `expect_identical` against `lm` on
  the extracted predictors, a factor's codes among them; its design becomes
  `makeModelMatrixFromDataFrame` of the same frame (3 levels at n 200 build dense, so the pins stay
  bitwise).
- [test-data-mixed.R](../../inst/tinytest/test-data-mixed.R), the block counting
  ["dbartsSparseSigmaFallbackWarning"](../../inst/tinytest/test-data-mixed.R): zero
  dbartsSigmaFallbackWarning for the sparse frame and its dense equivalent, the two `sigest` within 1e-10.
- [test-indicator-storage.R](../../inst/tinytest/test-indicator-storage.R):
  ["a sparse column the caller supplied still falls back"](../../inst/tinytest/test-indicator-storage.R)
  becomes no warning and the dense fit's value within 1e-10; the weights, offset and missing-value loop keeps
  its 1e-8; ["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R) pins every
  perturbation (0, 1e-9, 5e-8, 2e-7, 1e-6, 1e-5, 3e-5, 1e-4, 1e-3) within 1e-10 of the nearest of three dense
  fits (all columns, without `a`, without `c`), and 1e-6 asserts one rank fewer than lm.fit's.
- [test-sampler-splitProbabilities.R](../../inst/tinytest/test-sampler-splitProbabilities.R): a seeded
  comparison of split counts on a design with a factor fails at its seed. The implementer reports the
  counts before and after over five seeds; if the comparison holds in distribution the seed moves, and if
  not it is a stop.

New, in test-sparse-starting-sigma.R (under 5 s), against `lm.fit` on the dense indicator design within
1e-10 unless said:

- dec-B422, dense: a frame with a 5-level factor and a numeric column, default against
  `factors = "indicators"`: `sigest` identical, and equal to `summary(lm(y ~ f + x))$sigma`; the same with
  the factor ordered, with a missing value, with two levels and with 40 (built sparse); `data@x` and its
  `varTypes` identical to the fit given sigest (the trees' side is untouched).
- dec-B422, sparse: a frame with a sparseVector column and a dense factor; a sparseFactor with its first
  and a middle level as reference, and a middle reference with missing values.
- dec-B403: a bare dgCMatrix through `dbarts`, `bart` and `xbart` with no warning (counted, Gate hygiene);
  one-hot 3 x 40 plus three duplicated columns with zero weights and an offset; NA entries in a dgCMatrix; a
  dense timestamp-like column (mean 1.7e9, spread an hour) carrying signal, weighted and not, and the same
  stored in every row of a dgCMatrix; a column scaled by 1e160; p > n with rank below n - 1; p > n of full
  rank (NA, and through `dbarts` sd(y) with one dbartsSigmaFallbackWarning); that design with one row at
  weight 1e-14 (NA).
- Route, at n 2600 and 1 percent: p 1999 takes the exact routine and p 2000 LSQR (read from the route the
  routine reports), and the p 2000 design with the cutoff raised gives the exact value within 1e-6 of
  LSQR's; n 2200, p 2100: no estimate, sd(y), the warning of Open calls 1.
- LSQR forced: condition (b); a full block counted (degrees of freedom equal to the exact rank's), a block
  under a drop pattern that removes a present level not counted; a constant fully stored column under
  weights dropped; a zero response; the cap at 3 iterations returns a larger sigma than the converged one
  and raises nothing; two calls `identical`, for each route.
- An Inf entry in a dgCMatrix and in a sparseVector column: the dense path's error, naming the column.
- xbart on a sparse frame: each fold's sigma equals the dense indicator design's for the same rows.
- A fixed residual prior, dense and sparse: `estimateStartingSigma` is not called, the slot and `sigest`
  equal the fixed sigma with and without an agreeing sigest, draws identical to the build before; then
  `setModel` to a chisq prior draws finite sigmas. test-sigest-fixed-agree.R passes unchanged.

The 37 files that reach a factor design must otherwise pass unchanged; a failure there that is not a pinned
sigest, slot or seeded draw of a default fit with a factor is a stop.

Reviewer's mutants, each of which must fail a test: the expansion removed (codes); an ordered factor left
as codes; a reference-level indicator added to a sparseFactor; imputation before the expansion; cutoff 2000
to 20000, and `<=` to `<`; the 10 percent rule removed, and taken against p; b forced to 0, and counted
without the row-sum check; the constant test replaced by the sum-of-squares one; weights left out of the
operator, of the means, of the residual; the stopping tolerance at 1e-2; sigma from the recurrence's
residual norm; the intercept left out of the operator; zero-weight rows counted in n; exact tolerance 1e-10
to 1e-16 and to 1e-6; centering removed, and at the unweighted mean; the max-abs step removed; the
equilibration removed; the wide basis without D; `n - r` replaced by `n - p`; the Inf check removed; the
fixed-prior skip leaving the slot NA.

## Baselines

- Current: equivalence-e4faed5c, bcf-equivalence-1b7d730c, multinomial-equivalence-80b1c8d4
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)).
- Moves, by class: POSTERIOR-CHANGING the seven caller-sparse and the four dense-factor scenarios (Context);
  SHIFTING wideFactorIndicators (rounding; bitwise if the routine happens to round to the QR's value).
  NEUTRAL the other 43, bcf's 15 and multinomial's 11, the four snapshot files and every exact gate. Any
  other mover is a defect: stop.
- Re-record on the reference build (`--preclean --configure-args=--enable-reference-build`) with
  `EQUIVALENCE_SCENARIOS` set to the twelve and `EQUIVALENCE_CORES=2`, merged into a copy of e4faed5c in its
  scenario order, named after the slice's code commit; e4faed5c demoted to historical.
- Partition against e4faed5c in z mode: 43 of 55 identical and the movers exactly the twelve.
  wideFactorIndicators shows no |z| above 4. The four dense-factor scenarios' anchors move 0.08 to 4 percent
  and the seven sparse ones' fall by factors of 0.16 to 0.62 on 150 to 500 rows, so |z| above 4 is expected
  among the seven and possible among the four; the verdict there is which scenarios moved. The merged file
  reproduces 55 of 55 under `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17): the change is the value a prior is calibrated against, not the sampler; the
  identity is agreement with `lm.fit` on the indicator design (the tinytest pins, and the 400-design sweep
  rerun against the implemented function with factor columns added to a third of its designs).

## Gates

On the slice tip against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing):

- tests/cpp: unchanged, green (no C++ touched).
- Full tinytest suite: green with the edits above.
- The four seeded-drift snapshot files on the reference build: pass unchanged.
- The equivalence trio: as Baselines; bcf and multinomial bitwise against their current files.
- exact-gates.yaml's list in quick mode: all pass, output as before.
- `R CMD check --as-cran` from a clean tarball: no new NOTE; lintr, air, rc-codoc, win-drift,
  doc-freshness. Sanitizers are not owed: no compiled code changes.
- Speed, same machine, within 1.5x: the exact routine at n 2e4, p 2000 (1.2 s) and n 1e4, one-hot 20 x 50
  (0.20 s); LSQR at m 1e4, 1 percent (0.15 s) and m 5e4, 1 percent (5.4 s); xbart's 200 per-fold estimates
  at n 1e4, p 1000, 1 percent (28 s). No bench-sampler.R compare (no hot path).
- Build conditions (a) and (b) with their verdicts.

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd and man/xbart.Rd: the clause on sparse-backed
  columns becomes: a sparse x takes the same estimate as its dense equivalent; while the smaller of its row
  and column counts is at most 2,000 the fit is exact, above that iterative and never low, and above that
  with under 10 percent residual degrees of freedom the marginal standard deviation is used; supply
  `sigest` to skip the estimate. A factor enters the estimate as indicator columns whatever `factors` says.
  Under a fixed residual prior no estimate is made. man/xbart.Rd also gains the sentence the other pages
  carry on an agreeing `sigest` beside a fixed prior, which it lacks.
- man/bart.Rd's warning-class paragraph drops the class; `sigest` under a fixed prior is the fixed sigma.
- man/sparseFactor.Rd: the starting-sigma paragraph becomes one sentence: a sparse column leaves the default
  starting sigma as the same column stored dense gives it.
- inst/NEWS.Rd: the class leaves the warning-class list (never released, so no entry of its own; dec-B422
  restores 0.9-34's value, so none either).
- docs/design: sparse-columns.md's [R surface](../design/sparse-columns.md#r-surface) gets a dated paragraph
  with the rule and Algorithm in brief; starting-sigma-sensitivity.md a dated section with condition (a)'s
  table and verdict; error-style.md drops the class; memory-footprint.md's starting-sigma row gains the
  sparse routes (about 16 m^2 bytes, at most 64 MB, for the exact routine; vectors only for LSQR), per
  worker under xbart.
- At landing: TODO's item goes; a ledger entry for the calls below; this plan's Status and Landing note;
  the MANIFEST row.

## Steps

1. Change 1 to 5 with Algorithm A to F; the tinytest edits and new tests; the suite green against
   `R CMD INSTALL -l <lib> .`.
2. The 400-design sweep against the implemented function; build condition (a), its table and verdict
   written into the note. A verdict other than PASS stops the slice here.
3. The speed points; help, NEWS and docs.
4. After review: the re-record, MANIFEST row and partition, in their own commit.

## Stop conditions

Stop and report when: a build condition's verdict is not PASS; the diff passes ~1400 lines; an equivalence
scenario other than the twelve, a snapshot or an exact gate moves, or wideFactorIndicators shows a |z|
above 4; the sweep shows a rank different from the SVD's or a sigma off by more than 1e-10 above 30 residual
degrees of freedom outside the band; a speed point is past 1.5x; a tinytest outside the four files fails;
the change needs engine, bridge or C API code.

## Interactions

- Whichever lands second of this and any other slice re-recording equivalence.R re-records against the
  other's file and partitions against it.

## Calls made

- LSQR is written in R. Its time is Matrix's two sparse products a step, already compiled (3e-9 to 4e-9 s
  an entry a step); the R loop adds nothing measurable and can be interrupted. In C: about 250 lines, a
  bridge entry and its Windows twin, sanitizers owed, bits independent of Matrix's build, and no faster.
- At the cap the iterate is used silently, read from dec-B426's "with no message" and from condition (a)
  being the check on it. A warning at the cap would be about 6 lines and a test.
- The intercept is a column of the operator and every column is centered and scaled inside it, as measured;
  no other preconditioner. A block preconditioner is condition (a)'s remedy if it is needed.
- Two stops the prototype lacks: a reproduced response, and a zero norm at the start. A constant column is
  found by comparing values. All three from runs in Context.
- LSQR takes the exact routine's front end, a fold subsetting its rows (2.8 iterations' products) where the
  prototype passed held-out rows at weight 0. One front end; a few percent of a fold's time.
- m is taken after the expansion and after constant columns are dropped, on the rows kept; each xbart fold
  routes by its own m. A fold near the cutoff can take the other routine from the all-rows fit (they agree
  to 4e-5 or better where structure is known).
- A dense frame's factors are expanded by the indicators route's own builder, so default and "indicators"
  fits share one sigest bit for bit; the cost is that a default fit with a factor of five or more levels
  takes the sparse routine, 1e-11 from `lm.fit`, as an "indicators" fit does today.
- A caller's sparse source builds every factor without one level; full blocks come only from the indicators
  route and are recognized by `indicator.term` and a row-sum check. The alternative, reading widths off the
  `drop` attribute, miscounts: a one-level factor emits a column under forced sparse storage and none dense.
- "Large" is 10 percent, and a count-only overestimate past it also goes back (Build conditions).
- Condition (a)'s script is tracked in benchmarks/R so the reviewer reruns it (about 160 lines, unshipped).
- One routine for every sparse design (the QR goes); the smaller Gram side with the intercept as a column;
  centering of dense-backed and fully stored columns; max-abs then unit-norm scaling; the wide side
  equilibrated; tolerance 1e-10; residual computed directly; no refinement.
- An infinite entry in a sparse x is refused as the dense path refuses it.
- Under a fixed residual prior no estimate is made and the slot holds the fixed sigma, also beside an
  agreeing sigest (Change 4).
- dbartsSparseSigmaFallbackWarning retired.

## Open calls

1. Does the 10 percent case warn? Background: dec-B426 says LSQR runs "with no message" and that under 10
   percent residual degrees of freedom above the cutoff the design "takes sd(y)"; it does not say whether
   that sd(y) is announced. Every other sd(y) fallback raises dbartsSigmaFallbackWarning (run: a dense
   design with p >= n; a sparse one with no residual rank at or below the cutoff will too), and xbart reads
   that warning to choose its per-fold route. A wide sparse design, the common large one, lands here.
   Options: (a) warn with the existing class, the message naming the rule and `sigest`: one rule, sd(y) is
   always announced; the cost is a warning on every such fit until sigest is given, as wide dense fits
   have. (b) Silent: about 8 lines (the routine reports its route so xbart can still choose); the cost is
   that a wide sparse fit warns at m 2,000 and not at 2,001, and `sigest` reports sd(y) unmarked. (c) Warn
   only with no residual degrees of freedom, silent between 0 and 10 percent: about 10 lines; matches the
   dense path where it warns and hides the new rule where it is new. Recommended: (a).
2. The 10 percent rule stops at the cutoff. Background: by dec-B426 it applies above m 2,000 only. A design
   with 3 percent residual degrees of freedom gets the exact estimate at m 2,000 (relative sd about 9
   percent on 60 degrees of freedom) and sd(y) at m 2,001 (15 to 164 percent above the exact estimate in
   the note's row). Options: (a) as ruled, the step stays: nothing to build. (b) The rule for every sparse
   design: about 4 lines, and a sparse design then differs from its dense equivalent below the cutoff,
   against dec-B403. (c) The rule for dense designs too: about 10 lines, posterior-changing for every
   design with few residual degrees of freedom, a wider re-record, and its own evidence (the note measured
   sd(y) against the linear estimate at 19 degrees of freedom only at n 200). Recommended: (a) for this
   slice, and (c) as its own TODO item if wanted; it does not block the build.

## Estimate

Blind critique of this revision, half a day. Implementer about two days: the two routines and the design
builder one, tests and docs one. Condition (a) about three hours of machine time at two jobs (the exact
reference at m 10,000 is two minutes and 1.8 GB a design, one at a time; the tight LSQR reference at m
50,000 up to 15 minutes a design). Gates about two hours of machine time, the re-record minutes. Review
with mutants half a day, and one fix round.
