# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403); revised after its first blind critique; revised 2026-10-10 for
dec-B421, dec-B422 and dec-B426 and for the sigest-beside-fixed landing (dec-B413); revised twice more from
its second and third critiques. Not built. Five calls are open ([Open calls](#open-calls)); the steps marked
"OC1" build the recommendation of the first and wait on its ruling.

agent: one sonnet implementer (R only, no engine code); one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING, for a continuous response under a chisq residual prior with no sigest, in two
cases: a caller's sparse x (dec-B403), and any fit that keeps a factor as a categorical predictor, the
default, dense or sparse (dec-B422). An indicator expansion R built sparse (dec-B370) moves by rounding.
NEUTRAL for a dense design with no categorical column, every fixed-unit family, every fit given sigest, and
the draws of a fit under a fixed residual prior.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~1200 lines (R/utility.R ~360 net of the QR's removal, R/spec.R ~30, R/xbart.R ~35, tinytest ~460,
benchmarks/R ~200, man ~40, docs ~90, NEWS ~2, MANIFEST one row).

## Goal

A sparse x gets the starting sigma its dense equivalent gets, with no warning for being sparse: `lm.fit`'s
own while the design is small (Open calls 1), an exact sparse routine while the smaller of its row and
column counts is at most 2,000, and LSQR above that. Every factor, in a dense design too, enters that
regression as indicator columns. The sparse QR and the class dbartsSparseSigmaFallbackWarning go. An
infinite entry in a sparse x is refused as in a dense one, and a fixed residual prior makes no estimate. The
tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Rulings, each covering what it names:

- dec-B403: "Yes, just use the sparse regression, no warning." A sparse x takes its dense equivalent's
  linear-model starting sigma; a caller can pass sigest and the help says so; how it is computed is measured.
- dec-B421: "Band is fine." The exact routine's rank tolerance may differ from lm's on near-copies (shown
  0.01 to 4.7 percent, and 14 percent; what was measured since is Open calls 1).
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
  with a CSC column (a bare dgCMatrix: 1.053 where lm gives 0.500). An indicator expansion R built sparse
  (`sparse.from.indicators`) goes to [`sparseResidualStandardError`](../../R/utility.R), Matrix's sparse QR;
  anything else to [`residualStandardError`](../../R/utility.R) (`lm.fit`).
  [`floorSigmaEstimate`](../../R/utility.R) turns a non-finite result into sd(y) with
  dbartsSigmaFallbackWarning: a dense design with p >= n warns; one with 2 residual degrees of freedom of
  40 rows gives an estimate and no warning.
- A factor kept categorical enters as its level codes, ordered or not: a 40-level factor and a numeric
  column (n 300) start sigma at 1.088 by default (sd(y) 1.095) and 0.499 under `factors = "indicators"`.
  `dbarts` (bcf is its `forests` argument), `bart`, and `dbartsData` with `dbartsSpec` give a frame and
  the same data as a matrix of codes carrying `varTypes` one value; a multinomial fit makes no estimate;
  `bartBT` refuses a matrix of codes; `rbart_vi` builds through `dbarts` (not run).
- The frame rebuilt from a dense container's column list gives, through
  [`makeIndicatorModelMatrix`](../../R/utility.R), the design values of the `factors = "indicators"` fit,
  whose indicator columns are stored sparse at a density of at most
  [`sparseIndicatorDensity`](../../R/utility.R). A missing factor value is imputed with each level's
  frequency; the block still sums to one.
- An infinite stored entry in a dgCMatrix is accepted; a dense one is refused, naming the column. A fixed
  residual prior with no sigest still runs the estimate ([`resolveSamplerSpec`](../../R/spec.R);
  [`xbart`](../../R/xbart.R) before it resolves the prior), and draws are identical whatever the slot
  holds; `dbartsSpec` documents that a value the data object carries survives, and stan4bart writes one.
- xbart on a sparse x takes the per-fold route "marginal" (each fold its own sd); a dense design "linear",
  where `foldData` fits [`sigmaDesignMatrix`](../../R/utility.R) of `data@x`, a factor's codes included.

Measured 2026-10-10 on builds of this plan's text (R's reference BLAS unless said); the comparison of the
candidate rules on one set of designs is in Open calls 1.

- `lm.fit` on a dense design: 0.09 s at n 1e4, p 100; 0.74 s at p 300; 1.15 s at p 380; 0.92 s at n 2000,
  p 860; 0.59 s at n 1e5, p 80; 0.13 s at n 1e6, p 7, the last two peaking at 195 and 225 MB of heap with
  the design. With n(p + 1) <= 6e6 and n(p + 1) min(n, p + 1) <= 1.2e9 it stays near 1 s and under 200 MB.
  A narrow design past that bound has more than 1,000 rows.
- The exact routine against `lm.fit`, 400 random designs (n 30 to 2000, p 5 to 1530, numeric, one-hot and
  mixed, a third with dependent columns, a third weighted with zero weights, a third with an offset): rank
  equal to lm.fit's on 398, sigma within 9.1e-12 above 30 residual degrees of freedom and 2.6e-8 at 30 or
  fewer; the other 2 are wide, 0.26 and 0.64 percent apart. 373 of the 400 fall under the dense bound.
  Neither the shift nor the recheck changes any of the 400; the largest direct residual of a dropped
  column is 1.4e-15.
- Its cost at m 2,000: 1.0 to 1.1 s and 90 to 160 MB over R's heap narrow (n 2e4, 1 percent); 3.9 s and
  170 MB wide and rank-deficient (n 2000, p 3600, rank 901). The recheck adds 0.3 s for 20 dropped columns
  and 3.0 s for 1,000 exactly dependent ones; with 999 near-copies all re-admitted it takes 67 s and 1 GB
  against 1.1 s and 110 MB, and 4.9 s and 250 MB when only the 64 largest are.
- Without the recheck a column is dropped when under 1e-5 of its norm is left after the kept ones (pivot
  1e-10 on the unit-diagonal crossproduct; an exact dependency's pivot reaches 1.5e-13). Measured
  directly, an exact dependency leaves at most 1.4e-15 and a doubtful column 1e-9 or more.
- LSQR: two runs give the same bits. Weighted fits (uniform, 20 percent zeros, 0/1, skewed, a few at
  1e-12, with an offset) agree with `lm.wfit` to 7e-11 in 11 to 50 iterations at n 6000, p 400. A full
  one-hot block left uncounted costs a degree of freedom, and is left uncounted when the row sums are
  taken after a dense-backed level is centered. Run past convergence on an exact dependency the iterate
  drifts (31 times too high at 1,000 iterations); the standard condition-estimate stop at 1e10 ends it at
  iteration 27. A reproduced-response stop taken against the uncentered norm stops a response of mean
  1e10 at the first iteration.
- LSQR's tolerance: at 1e-6 a nearly dependent design stops on tolerance far from the fit (53 iterations,
  143 percent high on a near-copy at m 10,000); at 1e-10 that one converges in 234. Slow but well-posed
  designs (20-column groups correlated at 0.99: 704 iterations at 1e-6) reach the cap at 1e-10 with 1e-8
  left. Under weights (n 6000, p 400): 1 percent of rows at 1e6, 797 iterations and exact; at 1e12, 105
  percent high at 1e-6 and 5 percent at the cap at 1e-10; lognormal(0, 8) weights, the cap either way,
  1.6 percent high.
- Bits are in Open calls 5; LSQR's costs in
  [Cost and memory](../design/starting-sigma-sensitivity.md#cost-and-memory).

What moves (every call of [`estimateSigmaFromLinearModel`](../../R/utility.R) logged through a load hook,
beside the value `lm.fit` gives on the indicator design):

- equivalence.R, quick mode, 132 calls in 44 of 55 scenarios. Twelve move:
  - seven on a caller's sparse source, sd(y) to the linear estimate (factors of 0.16 to 0.62):
    ["sparse <- list("](../../benchmarks/R/equivalence.R), ["mixedmatrix <- list("](../../benchmarks/R/equivalence.R),
    ["sparsefactor <- list("](../../benchmarks/R/equivalence.R), ["testswap <- list("](../../benchmarks/R/equivalence.R),
    ["leaffactormixed <- list("](../../benchmarks/R/equivalence.R), ["factorpartial <- list("](../../benchmarks/R/equivalence.R),
    ["xbartmixed <- list("](../../benchmarks/R/equivalence.R);
  - four on a dense frame with a factor, codes to indicators (0.08 to 4 percent):
    ["categorical <- list("](../../benchmarks/R/equivalence.R), ["leaffactor <- list("](../../benchmarks/R/equivalence.R),
    ["nafactor <- list("](../../benchmarks/R/equivalence.R), ["ordfactor <- list("](../../benchmarks/R/equivalence.R);
  - ["wideFactorIndicators <- list("](../../benchmarks/R/equivalence.R), by rounding.
  The other 32 with a call are dense with no categorical column; 11 make none.
- bcf-equivalence.R: 12 calls, all dense, no factor. multinomial-equivalence.R: none.
- The four test-reproducibility files (build guard bypassed, shipped build): 5 calls in three files, none
  sparse and none with a factor; the binary file stopped under the bypass with nothing logged.
- exact-gates.yaml's list in quick mode, bcf-latent-exact.R aside: 35 calls, 5 of them on a factor design
  (mask-redraw-exact.R 4, change-balance.R 1), all five under a fixed residual prior. Rerun with the new
  values substituted, both scripts print what they printed, two elapsed-time lines aside.
- tinytest (19709 results, 0 failures): 37 files reach a dense container with a factor (409 calls). With
  the new values substituted for every factor design, six expectations fail in four files
  ([Tests](#tests)); with a matrix of codes left as codes, three more fail in test-data-code-channel.R.

## Algorithm

`sparseResidualStandardError(y, x, weights, offset)` keeps its name; the dense bound, the cutoff (2000),
LSQR's cap (1000) and its tolerance are formals with their shipped defaults, so a test and the build's
reference run reach any route at any size.

A. The design, built by one new function that every caller uses (`startingSigmaDesign(x)`: the all-rows
estimate and xbart's chunk runner). The expansion of dec-B422 happens here and nowhere else; `data@x`, the
cut grid and the trees never see it.

1. A plain matrix with no factor marked: [`sigmaDesignMatrix`](../../R/utility.R) as today, for `lm.fit`.
   A matrix whose `varTypes` mark a factor is rebuilt as a frame (codes to factors over its
   `factor.levels`, or over the codes present where no table rides) and takes case 2.
2. A dense container (no CSC column) with a factor column: the frame is rebuilt from its column list and
   names and handed to [`makeIndicatorModelMatrix`](../../R/utility.R) with `drop = TRUE` and storage
   "auto", so it is the design `factors = "indicators"` builds, and a frame, its matrix of codes and its
   "indicators" fit share one sigest bit for bit. The result is a matrix (case 1) or a container with
   sparse-built indicators (case 4). Without a factor column: as today.
3. A caller's sparse source. A bare dgCMatrix is wrapped ([`wrapSparseTestMatrix`](../../R/mixedMatrix.R)).
   In a container, a dense-backed factor column becomes one indicator per present level but its first; a
   sparseFactor column (non-NA `sparseReference`) one indicator per present level but its reference, from
   its stored entries, a stored entry at the reference code dropped. A stored NA is an NA entry in each of
   the factor's indicators. Other columns as in [`sparseDesignMatrix`](../../R/utility.R) today, then its
   imputation (for an indicator, the level's frequency). The columns span with the intercept what the dense
   equivalent's full block spans. Each indicator records its factor (`indicator.term`): a fold or zero
   weights can remove the level left out, and the rest are then a full block. New columns are assembled
   from slots as [`assembleMixedMatrix`](../../R/mixedMatrix.R) does, with no Matrix coercion.
4. An indicators-route container (built by R, or beside a caller's sparse column): taken as it is.
   [`makeIndicatorModelMatrix`](../../R/utility.R) records per column the input term that emitted it
   (`indicator.term`, in place of `sparse.from.indicators`, whose one reader goes).

B. Route by size (OC1). A design that is a dense matrix goes to `lm.fit` at any size, as today. A
sparse-stored design (cases 3 and 4) with n(p + 1) <= 6e6 and n(p + 1) min(n, p + 1) <= 1.2e9, n its rows
with a response, weight and offset and a positive weight and p its columns, is made dense and goes to
[`residualStandardError`](../../R/utility.R): the fit a dense caller gets, with `lm.fit`'s rank rule and
bits (to rounding where A.3 built a factor without one level). Anything larger takes C to G. An infinite
entry is refused before either, as for a dense design, naming the column
([`nonFinitePredictorNames`](../../R/spec.R) learns sparse sources).

C. Front end of the sparse routine, in this order.

1. Rows: drop rows with a missing response, weight or offset and rows of weight 0, as `lm.wfit` does;
   z = y - offset; n the rows kept.
2. Constants: a column whose entries over the kept rows are all equal (implicit zeros included; an exact
   comparison of values) is dropped. p is the number of columns left and m = min(n, p + 1).
3. Blocks, on the values as they stand: the columns of one factor term (`indicator.term`) are a full block
   when at least two are left and they sum to one in every kept row (within 1e-8). b is their number.
4. OC1, the shift. A column with unstored kept rows and stored entries that are not all equal, in a design
   that also holds a column constant on exactly its unstored rows and unstored elsewhere (or constant on
   exactly its stored rows), has the mean of its stored entries taken off them. The partner is found by
   stored count and the sum and sum of squares of row indices, then confirmed row for row. The span with
   the intercept is unchanged, and a zero-filled column beside its missing indicator stops being a
   near-copy of the two. With no such partner the column is left alone.
5. Centering: each dense-backed column and each column stored in every kept row is centered at its weighted
   mean over the kept rows. Then every column is divided by its largest absolute entry.

D. Route. m <= 2000: the exact routine (E). Otherwise the structural residual degrees of freedom are
df = n - 1 - p + b; when df < 0.1 n there is no estimate (G); else LSQR (F).

E. The exact routine.

1. B = diag(sqrt(w)) [1 X]; each column divided by its norm. Sparse columns are not centered: a centered
   crossproduct cancels for a large-mean column, which is why the intercept stays a column.
2. Narrow (ncol(B) < n): `chol(as.matrix(Matrix::crossprod(B)), pivot = TRUE, tol = 1e-10)` with
   `suppressWarnings` around that one call (its one warning is the rank deficiency); r its rank, K its first
   r pivots, R its leading r x r block; e = z sqrt(w) - B_K R^-1 R^-T B_K' (z sqrt(w)), computed directly.
3. Wide: K = `as.matrix(Matrix::tcrossprod(B))`, equilibrated by D = sqrt(diag(K)), the pivoted Cholesky of
   D^-1 K D^-1 at the same tolerance. r >= n gives no estimate; else Q = `qr.Q(qr(D L, LAPACK = TRUE))`, L
   the first r columns of the factor, unpivoted: LAPACK's pivoted QR applies no rank tolerance, so the
   basis keeps all r columns and the rank is not decided twice. e = z sqrt(w) - Q Q' z sqrt(w).
4. sigma = sqrt(sum(e^2) / (n - r)); no estimate when n - r <= 0.
5. OC1, the recheck. Each vector the Cholesky dropped (a column of B when narrow, a row of D^-1 B when
   wide) is measured again where it lives: its residual against the kept vectors, projected three times
   through R, in chunks of at most 8e6 entries. One with more than 1e-7 of its unit norm left is a
   candidate; only the 64 with most left are held. They are orthogonalized in order of size, twice, and
   one still above 1e-7 is re-admitted: it joins the fit (narrow: e loses its component along it after one
   more projection through R; wide: B v joins the basis for a re-admitted row direction v) and r grows by
   one. An exact dependency leaves a residual at rounding and is never a candidate. The aim past the dense
   bound is the least-squares fit, not lm's rank: a column lm would drop by its order can be re-admitted.

F. LSQR (Paige and Saunders), in R in R/utility.R, on base R and Matrix's sparse products alone.

1. Operator. With h = sqrt(w), s = sum(w), mu_j the weighted mean of column j and d_j the reciprocal of its
   weighted centered norm, A = diag(h) [1 / sqrt(s), (X - 1 mu') diag(d)], applied as two sparse products a
   step (`X %*% v`, `Matrix::crossprod(X, u)`) with the centering as a rank-one correction. Unit column
   norms are the only preconditioning. A column whose centered sum of squares comes out 0 takes d_j = 0:
   it adds nothing and stays counted, so p, m and the route are those of C and D.
2. Right side b = h z, start at zero. If the norm of b is 0, sigma is 0; if A'b is 0, the residual is b.
3. Iteration: Golub-Kahan bidiagonalization with the standard updates, no reorthogonalization, no damping.
4. Stop at the first of:
   - the normal-equations estimate at or below the tolerance times the running operator norm times the
     residual norm. OC1: the tolerance is 1e-10; the ruling was shown 1e-6;
   - the residual norm at or below 1e-10 of the norm of h (z - zbar), zbar the weighted mean of z: the fit
     reproduces the response, where the first test cannot fire;
   - the standard condition estimate (the running operator norm times the norm of the accumulated update
     directions) at or above 1e10: the iteration has begun on a dependency, where the iterate drifts;
   - a zero alpha or beta; 1,000 iterations.
5. At the cap the iterate is used like any other, with no message.
6. sigma = sqrt(sum(r^2) / df), r = b - A x computed from the returned coefficients, not the recurrence.
7. Weights enter as the row scale h in the operator and the right side and in every mean and norm; n counts
   rows of positive weight. This is `lm.wfit`'s fit and `summary.lm`'s degrees of freedom.
8. Degrees of freedom are D's count, never the iteration's. A dependency the structure does not show is
   counted as a fitted column and overestimates by sqrt((df + d) / df).
9. Determinism. No random number, fixed constants, one operation order; the bits held across two BLAS and
   two thread counts. The exact routine's last digit follows the BLAS and, under a threaded one, its
   thread count (Open calls 5); `lm.fit`'s follows the BLAS only.

G. No estimate (E with no residual rank, D under 10 percent) returns NA, which
[`floorSigmaEstimate`](../../R/utility.R) takes to sd(y - offset) with dbartsSigmaFallbackWarning, as the
dense path does. Whether the 10 percent case warns is Open calls 3; the steps build its recommended option,
and the other is about 8 lines. An allocation failure is an error, reported as today.

## Change

1. R/utility.R: `startingSigmaDesign` (A); [`sparseDesignMatrix`](../../R/utility.R) gains the wrap and the
   two expansions; [`sparseResidualStandardError`](../../R/utility.R) rewritten as B to G with the size
   route, the shift, the exact routine, its recheck and LSQR as internal functions (the batching, the
   refactor loop and the `grepl` muffler of Matrix's warning go); [`estimateSigmaFromLinearModel`](../../R/utility.R) loses the fallback
   branch and its warning and routes by what the design builder returned;
   [`makeIndicatorModelMatrix`](../../R/utility.R) records `indicator.term` and no longer sets
   `sparse.from.indicators`.
2. R/spec.R: [`nonFinitePredictorNames`](../../R/spec.R) for sparse sources (the dense list and the CSC
   block's stored entries, named through the container's map, or positions for a bare dgCMatrix).
3. R/xbart.R: the design is built once per chunk by `startingSigmaDesign` in
   [`xbartRunUnits`](../../R/xbart.R), for dense sources too, so a fold of a dense frame with a factor fits
   indicators; `foldData` fits its training rows by the same routes, each fold by its own size and its own
   degrees of freedom. A fold with no estimate falls back as a dense fold does today; where the all-rows
   estimate fell back every fold takes its marginal sd, as today.
4. A fixed residual prior makes no estimate. In [`resolveSamplerSpec`](../../R/spec.R), where `residPrior`
   is a dbartsFixedPrior and the family is not on a fixed unit scale, a slot that is NA takes the square
   root of the fixed variance; a slot that holds a value keeps it (an agreeing sigest, or what the data
   object carried). [`xbart`](../../R/xbart.R) resolves its prior before its all-rows estimate and does the
   same, with no per-fold route. [`refuseSigestUnderFixedPrior`](../../R/family.R) still runs first and is
   untouched. Draws are unchanged; with nothing given, `bart`'s `sigest` and the slot report the fixed
   sigma, and a later `setModel` to a chisq prior, or `rbart_vi`'s start, reads it.
5. The class dbartsSparseSigmaFallbackWarning is retired: new in 1.0, nothing raises it, and no consumer
   branch names it (stan4bart, bartCause, treatSens, bairrtt).

## Constraints

- No engine, bridge, C API or state change. No new dependency: Matrix (Suggests, no version floor declared)
  as today. The Matrix calls made (`crossprod`, `tcrossprod`, `rowSums`, `cbind2`, `Diagonal`, `%*%`,
  subsetting) are all exported by Matrix 1.4-1, R 4.2's (read from its NAMESPACE, not run).
- A dense design with no categorical column under a chisq prior is bit for bit unchanged.
- Out of scope: a cutoff or LSQR for the dense path; the 10 percent rule below the cutoff (Open calls 4); a
  per-fold crossproduct downdate or warm start in xbart.

## Build conditions

Both are dec-B426's. The rules are fixed here, before any number exists; the threshold T in them is Open
calls 2 (recommended: 10 percent). They are written for Open calls 1's recommendation; what changes under
its other options is said at the end.

(a) LSQR on badly conditioned realistic designs. A tracked script, benchmarks/R/starting-sigma-lsqr.R, calls
the implemented routine; two jobs at a time. Each family has its own reference, the least-squares fit of
every column obtained the way that family allows; a size with no admitted reference is reported, not gated.

- Crossed factors, by occupancy: two factors of L levels, K of the L^2 cells occupied, every occupied cell
  one row and the rest of the rows spread over them with Zipf frequencies; columns are the two main
  effects and the occupied cells as one-hot columns of a dgCMatrix. (L, K) = (150, 9,699), (200, 19,599),
  (300, 49,399) give m 10,000, 20,000 and 50,000 exactly, at n = 1.5 m and n = 1.12 m. Reference: the
  cell-means fit. The count alone is 3.0, 2.0 and 1.2 percent high at n = 1.5 m and 11.8, 8.0 and 4.9 at
  n = 1.12 m.
- Word counts (term frequencies Zipf, document lengths lognormal, raw counts) and correlated numeric
  columns (groups of 20 sharing one sparse pattern at 5 percent, pairwise correlation 0.99, beside 200
  fully stored columns correlated at 0.999), at the same six (m, n). Reference: the exact routine at m
  10,000 (cutoff and dense bound raised, one job at a time); above it LSQR at tolerance 1e-13 and a cap of
  30,000, admitted if it stops on tolerance and the same setting is within 1e-6 of the exact routine at m
  10,000.
- Ten zero-filled timestamp columns (spreads from a minute to a month), each beside its missing indicator
  and carrying signal, in a sparse numeric design, at the six (m, n) and in both column orders. Reference:
  LSQR at 1e-13 on the shifted design (134 iterations at m 10,000), checked on the m 2,000 twin against a
  dense QR at tolerance 1e-13 (both 0.4840, where `lm.fit` gives 0.693 or 0.842 by column order).
- Nearly dependent columns with no partner to shift against: 50 near-copies at 1e-8 to 1e-3, half
  carrying signal; ten pairs of zero-filled start and end times sharing one missing pattern. Reference:
  the dense QR at tolerance 1e-13 on an m 2,000 twin, the exact routine with its recheck at m 10,000, the
  tight LSQR run above that where it is admitted.
- At m 10,000, the first three families also under skewed weights (a cubed exponential) and with 1 percent
  of rows at weight 1e6, against the exact routine. One percent of rows at 1e12 is reported, not gated.
- Recorded per run: iterations, seconds, stop reason, the residual sum of squares against the reference's,
  and sigma against the reference's.
- PASS: every gated run's sigma is within T of its reference, in either direction, whatever its stop
  reason.
- BACK TO THE MAINTAINER, before landing: any gated run further than T from its reference, with the table
  and three priced options (a higher cap, sd(y) there, a block preconditioner). Known now: the crossed
  design at m 10,000, n 11,200 is 11.8 percent high by its count.
- DEFECT, stop: a run at the shipped tolerance whose residual sum of squares is more than 1e-6 below that
  of the tolerance-1e-13 run on the same design; or an m 2,000 twin further than 1e-6 from its dense QR.
- A run more than 1 percent from its reference, or at the cap, is named in the Landing note.
- Under Open calls 1's options 1 to 4 the fourth family is reported and not gated (it is the band as
  ruled); under options 1 to 3 the timestamp family has no reference the options can reach and is
  reported the same way.

(b) Weighted fits against the exact routine, in tinytest with LSQR forced (cutoff 0, dense bound 0) on small
designs: sparse numeric, an indicators-route design with a full block that has a level above 20 percent
density, a frame with a sparseFactor; under uniform weights, 20 percent zero weights, 0/1 fold weights,
skewed weights, a few rows at 1e-12, and weights with an offset. Each weighted case is read beside its
unweighted twin: both within 1e-6 relative with equal degrees of freedom, PASS; the weighted one past it and
the twin within, FAIL as a defect in how weights enter; the twin past it, the design is wrong for this test
and is replaced (these three are within 1e-10 unweighted).

## Tests

Edits the substitution run found:

- [test-starting-sigma.R](../../inst/tinytest/test-starting-sigma.R) pins `expect_identical` against `lm` on
  the extracted predictors, a factor's codes among them; its design becomes
  `makeModelMatrixFromDataFrame` of the same frame, and the pins stay bitwise.
- [test-data-mixed.R](../../inst/tinytest/test-data-mixed.R), the block counting
  ["dbartsSparseSigmaFallbackWarning"](../../inst/tinytest/test-data-mixed.R): zero
  dbartsSigmaFallbackWarning for the sparse frame and its dense equivalent, the two `sigest` within 1e-10.
- [test-indicator-storage.R](../../inst/tinytest/test-indicator-storage.R):
  ["a sparse column the caller supplied still falls back"](../../inst/tinytest/test-indicator-storage.R)
  becomes no warning and the dense fit's value; the weights, offset and missing-value loop keeps its 1e-8;
  ["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R) pins every perturbation (0,
  1e-9, 5e-8, 2e-7, 1e-6, 1e-5, 3e-5, 1e-4, 1e-3) to the dense fit with `expect_identical`: the design is
  under the dense bound (OC1), so no pin sits near a threshold.
- [test-sampler-splitProbabilities.R](../../inst/tinytest/test-sampler-splitProbabilities.R): a seeded
  comparison of split counts on a design with a factor fails at its seed. The implementer reports the
  counts over five seeds; if the comparison holds in distribution the seed moves, and if not it is a stop.
  test-data-code-channel.R passes unchanged and is the pin for A.1.

New, in test-sparse-starting-sigma.R (under 10 s), against `lm.fit` on the dense indicator design, missing
values imputed as the estimate imputes them, within 1e-10 unless said; "forced" means the dense bound at 0:

- dec-B422, dense: a frame with a 5-level factor and a numeric column, default against
  `factors = "indicators"` and against its matrix of codes: `sigest` identical; equal to
  `summary(lm(y ~ f + x))$sigma` where nothing is missing; the same with the factor ordered, with two
  levels, with 40, and with a missing value against the imputed design; `data@x` and its `varTypes`
  identical to the fit given sigest.
- dec-B422, sparse: a frame with a sparseVector column and a dense factor; a sparseFactor with its first
  and a middle level as reference, and a middle reference with missing values; each also forced.
- dec-B403: a bare dgCMatrix through `dbarts`, `bart` and `xbart` with no warning (counted, Gate hygiene)
  and `sigest` identical to the dense matrix's (OC1). Forced: one-hot 3 x 40 plus three duplicated columns
  with zero weights and an offset; NA entries; a fully stored timestamp-like column carrying signal,
  weighted and not; a column scaled by 1e160; p > n with rank below n - 1; p > n of full rank (NA, and
  through `dbarts` sd(y) with one dbartsSigmaFallbackWarning); that design with one row at weight 1e-14.
- OC1, under the bound: the zero-filled timestamp with its indicator (an hour and five minutes, both
  column orders), a near-copy with signal, a start and end time, a row at weight 1e12: each `identical`
  to `lm.fit` on the dense design.
- OC1, forced: ten zero-filled timestamps with indicators at 300 rows, both column orders, by the exact
  routine and by LSQR, within 1e-8 of a dense QR at tolerance 1e-13; a column with unstored rows and no
  partner is not shifted, nor the indicator itself, nor a column whose partner matches in count and sums
  but not in rows; a near-copy at 1e-6 and two zero-filled times sharing one pattern are re-admitted,
  within 1e-8; of 100 near-copies 64 are.
- Route: a sparse design just under and just over each limit of the dense bound takes `lm.fit` and the
  sparse routine (read from the route the routine reports). At n 2600 and 1 percent, p 1999 takes the
  exact routine and p 2000 LSQR, and the p 2000 design with the cutoff raised gives the exact value within
  1e-6 of LSQR's; n 2200, p 2100: no estimate, sd(y), the warning of Open calls 3. With the cutoff formal,
  forced: n 40, p 100 at cutoff 45 is exact (m is n); 20 live and 30 constant columns at n 200, cutoff 25,
  is exact (m counts live columns).
- LSQR forced (cutoff 0): condition (b); a full block with a level above 20 percent density counted; a
  block under a drop pattern that removes a present level not counted; a sparseFactor whose reference
  level has only zero-weight rows counted as a block; a constant fully stored column under weights
  dropped; a zero response; a response of mean 1e10 within 1e-6 of the exact routine; a response the
  design reproduces (sigma under 1e-8 of sd(y)); a one-hot design at tolerance 0 stops on the condition
  estimate within 1e-6; the cap at 3 iterations returns a larger sigma than the converged one and raises
  nothing; two calls `identical`, for each route.
- An Inf entry in a dgCMatrix and in a sparseVector column: the dense path's error, naming the column.
- xbart: on a sparse frame and on a dense frame with a factor, each fold's sigma equals `lm.fit`'s on the
  indicator design for the same rows; on n 3300, p 2100 at 1 percent each fold takes LSQR (traced).
- A fixed residual prior, dense and sparse: `estimateStartingSigma` is not called; with nothing given the
  slot and `sigest` equal the fixed sigma; an agreeing sigest and a sigma the data object carried stay;
  draws identical to the build before; then `setModel` to a chisq prior draws finite sigmas.
  test-sigest-fixed-agree.R passes unchanged.

The 37 files that reach a factor design must otherwise pass unchanged; a failure there that is not a pinned
sigest, slot or seeded draw of a default fit with a factor is a stop.

Reviewer's mutants, each of which must fail a test: the expansion removed (codes), for a container, for a
matrix of codes, and in `foldData` alone; an ordered factor left as codes; a reference-level indicator
added to a sparseFactor; `indicator.term` not recorded by the caller's-sparse builder; imputation before
the expansion; the size route removed, and each of its limits inverted; cutoff 2000 to 20000, and `<=` to
`<`; m taken as p + 1, and before constants are dropped; the 10 percent rule removed, and taken against p;
b forced to 0, counted without the row-sum check, and counted after centering; the constant test replaced
by a sum-of-squares one; the shift without its partner test, without the row-for-row confirmation, and
applied to a column whose stored entries are all equal; weights left out of the operator, of the means, of
the residual; the LSQR tolerance at 1e-2; the reproduced-response stop against the uncentered norm; the
condition stop removed; sigma from the recurrence's residual norm; the intercept left out of the operator;
zero-weight rows counted in n; exact tolerance 1e-10 to 1e-16 and to 1e-6; the recheck removed, its
threshold at 1e-5, its repeat projections removed, its bound of 64 removed; centering removed, and at the
unweighted mean; the max-abs step removed; the equilibration removed; the wide basis without D, and from
`qr` at its default tolerance; `n - r` replaced by `n - p`; the Inf check removed; the fixed-prior skip
leaving the slot NA, and overwriting a carried value.

## Baselines

- Current: equivalence-e4faed5c, bcf-equivalence-1b7d730c, multinomial-equivalence-80b1c8d4
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)).
- Moves, by class: POSTERIOR-CHANGING the seven caller-sparse and the four dense-factor scenarios (Context);
  SHIFTING wideFactorIndicators (rounding, the QR's value to `lm.fit`'s). All twelve are under the dense
  bound, so each takes `lm.fit`'s value on its indicator design. NEUTRAL the other 43, bcf's 15 and
  multinomial's 11, the four snapshot files and every exact gate. Any other mover is a defect: stop.
- Re-record the twelve on the reference build (`--preclean --configure-args=--enable-reference-build`,
  `EQUIVALENCE_CORES=2`), merged into a copy of e4faed5c in its scenario order and named after the slice's
  code commit; e4faed5c demoted to historical. Partition against it in z mode: 43 of 55 identical and the
  movers exactly the twelve, with no |z| above 4 on wideFactorIndicators; |z| above 4 is expected among
  the seven and possible among the four. The merged file reproduces 55 of 55 under
  `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17): the change is the value a prior is calibrated against, not the sampler; the
  identity is agreement with `lm.fit` on the indicator design (the tinytest pins, and the 400-design sweep
  rerun against the implemented function, forced, with factor columns added to a third of its designs).

## Gates

On the slice tip against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing):

- tests/cpp unchanged and green (no C++ touched); the full tinytest suite green with the edits above; the
  four seeded-drift snapshot files on the reference build pass unchanged.
- The equivalence trio as Baselines, bcf and multinomial bitwise against their current files;
  exact-gates.yaml's list in quick mode all pass, output as before.
- `R CMD check --as-cran` from a clean tarball: no new NOTE; lintr, air, rc-codoc, win-drift,
  doc-freshness. Sanitizers are not owed: no compiled code changes.
- Speed and memory, same machine, time within 1.5x and peak heap over the design within 1.5x: `lm.fit` at
  the dense bound (n 1e4, p 380: 1.2 s, 90 MB); the exact routine at n 2e4, p 1999, 1 percent (1.1 s, 160
  MB), with 999 near-copies of which 64 are re-admitted (4.9 s, 250 MB), and on the wide design n 2000,
  p 3600, rank 901 (6.2 s, 180 MB); LSQR at m 1e4 and m 5e4, 1 percent, timed at the shipped tolerance by
  the implementer and held by the reviewer; xbart's 200 per-fold estimates at n 1e4, p 1000, 1 percent
  (28 s). No bench-sampler.R compare (no hot path).
- Build conditions (a) and (b) with their verdicts.

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd and man/xbart.Rd: a sparse x takes the same
  estimate as its dense equivalent: `lm`'s own while the design is small, a sparse regression past that,
  iterative above 2,000 rows or columns, where under 10 percent residual degrees of freedom gives the
  marginal standard deviation; supply `sigest` to skip it. A factor enters as indicator columns whatever
  `factors` says. Under a fixed residual prior no estimate is made. man/xbart.Rd gains the other pages'
  sentence on an agreeing `sigest` beside a fixed prior; man/dbartsSpec.Rd's "an unset value is still
  estimated" gains "under a chisq prior"; man/bart.Rd's warning-class paragraph drops the class;
  man/sparseFactor.Rd's starting-sigma paragraph becomes one sentence (a sparse column leaves the default
  starting sigma as the same column stored dense gives it).
- inst/NEWS.Rd: the class leaves the warning-class list (never released; dec-B422 restores 0.9-34's value,
  so no entry for either). docs/design: sparse-columns.md's
  [R surface](../design/sparse-columns.md#r-surface) gets a dated paragraph with the rule;
  starting-sigma-sensitivity.md a dated section with condition (a)'s table and verdict; error-style.md
  drops the class; memory-footprint.md's starting-sigma row gains the sparse-source routes (up to 200 MB
  with the design under the dense bound, about 250 MB over the design for the exact routine, vectors only
  for LSQR), per worker under xbart.
- At landing: TODO's item goes; a ledger entry for the calls below; this plan's Status and Landing note;
  the MANIFEST row.

## Steps

0. Open calls 1 is ruled, or the implementer is told which option to build. The steps are written for its
   option 5. Under 4, E.5 goes and F.4's tolerance is 1e-6; under 1 to 3, route B and C.4 go as well (2
   keeps E.5), and the OC1 tests pin what the chosen option gives.
1. Change 1 to 5 with Algorithm A to G; the tinytest edits and new tests; the suite green against
   `R CMD INSTALL -l <lib> .`.
2. The 400-design sweep against the implemented function, forced; build condition (a), its table and
   verdict written into the note. A verdict other than PASS stops the slice here.
3. The speed and memory points; help, NEWS and docs.
4. After review: the re-record, MANIFEST row and partition, in their own commit.

## Stop conditions

Stop and report when: a build condition's verdict is not PASS; the diff passes ~1600 lines; an equivalence
scenario other than the twelve, a snapshot or an exact gate moves, or wideFactorIndicators shows a |z|
above 4; the sweep shows a rank different from lm.fit's on a design other than its two wide ones, or a
sigma off by more than 1e-10 above 30 residual degrees of freedom; a speed or memory point is past 1.5x; a
tinytest outside the four files fails; the change needs engine, bridge or C API code. Whichever lands
second of this and any other slice re-recording equivalence.R re-records against the other's file.

## Calls made

- LSQR is written in R: its time is Matrix's two compiled sparse products a step, and the loop can be
  interrupted. In C: about 250 lines, a bridge entry and its Windows twin, sanitizers owed, no faster.
- At the cap the iterate is used silently, read from dec-B426's "with no message" and from condition (a)
  being the check on it. A warning at the cap would be about 6 lines and a test.
- The intercept is a column of the operator and every column is centered and scaled inside it; no other
  preconditioner. Three stops the prototype lacks: a reproduced response (against the centered norm), a
  zero norm at the start, the condition estimate.
- A factor is expanded wherever the estimate meets one, a matrix of codes included, by the indicators
  route's own builder, so the doors share one sigest. A caller's sparse source builds every factor
  without one level; full blocks are found by a recorded term and a row-sum check on the raw values.
- LSQR takes the exact routine's front end, a fold subsetting its rows. m is taken after the expansion
  and constant drops, on the rows kept; each xbart fold routes by its own size.
- Under a fixed residual prior no estimate is made and only an empty slot is filled. The alternative,
  overwriting whatever the slot holds, gives one value on every path but breaks what `dbartsSpec`
  documents (a carried value survives) and what stan4bart writes; it is one line either way.
- Within Open calls 1's recommendation: the dense bound's two numbers are where `lm.fit` measured about
  1 s and 200 MB; the recheck holds the 64 columns with most left (unbounded it took 67 s and 1 GB on a
  thousand); the shift needs its partner column in the design and is confirmed row for row; past the
  bound the aim is the least-squares fit of every column, so no step imitates lm's order.
- Condition (a) has a reference per family and gates only where one exists; its verdicts are two-sided,
  and a miss from the count alone is treated as the solver's would be. Its script is tracked.
- Carried from the first revision: the smaller Gram side with the intercept as a column; centering,
  max-abs and unit-norm scaling; the wide side equilibrated; a sparse Inf refused; the class retired.

## Open calls

Each is the maintainer's; none is settled here.

### 1. What rule gives a sparse design its starting sigma when columns are nearly dependent?

What the number is for. The starting sigma is the residual sd of a linear regression of the response on
the predictors; the default prior on the residual variance is centered with it. A wrong one matters at
small n and fades: an estimate 2 to 4 times too high moved the posterior mean of sigma by 23 to 57 points
of the truth at n 200, by 2 to 6 points at n 1000, and not measurably at n 5000
([Results](../design/starting-sigma-sensitivity.md#results)).

What was ruled, and what was measured since. dec-B403 gives a sparse design its dense equivalent's sigma.
dec-B421 ("Band is fine.") accepted that the sparse routine drops a column lm keeps when between 1e-7 and
1e-5 of the column is independent of the others, shown differences of 0.01 to 4.7 percent and one of 14.
Three critiques then measured ordinary-looking designs where the difference is 130 to 270 percent:

- a timestamp holding 0 where the value is missing, beside the indicator of those rows, the usual way to
  carry a missing value in a sparse matrix: lm 0.536, the routine 1.796, sd(y) 1.899;
- the same two columns in a plain dense data frame that also holds a factor, which dec-B422 sends through
  the same routine: 1.714 against lm's 0.465, where 0.9-34 gives 0.531 and today's build 0.533 on like
  data, so as ruled this is a regression;
- a near-copy of a column that carries signal (140 percent high), and a few rows weighted 1e12 times the
  rest (175 percent high).

Above 2,000 columns the iterative solver, stopped at 1e-6, does the same (168 percent high at 10,000).

lm is not a fixed target on such designs. Its answer depends on the order of the columns: the same data
with a timestamp spread over five minutes gives 0.536 with the timestamp before its indicator and 1.796
with the indicator first, and a design with ten such timestamps gives 0.693 or 0.842. The least-squares
fit that keeps every column gives 0.536 and 0.484. The other way round, with a start and an end time a
few minutes apart lm drops the second and gives 1.059 where the full fit is 0.511. So there are two
targets: lm's answer, which is what a dense caller gets, and the full least-squares fit, which is the
better estimate.

Options. Percentages are above the full least-squares fit.

1. As ruled. A sparse design gets a sparse regression that calls a column redundant when under 1e-5 of it
   is independent of the rest; above 2,000 columns an iterative solver stops at 1e-6. A user gets lm's
   sigma on ordinary designs and 132 to 269 percent too much on the designs above, at every size, a
   default data frame with a factor included. Nothing more to build.
2. As ruled, plus a recheck. A column the routine drops is measured again directly and kept when more
   than 1e-7 of it is independent; the iterative solver stops at 1e-10. The one-hour timestamp, the
   near-copy and the weights come out right; a five-minute timestamp does not (235 percent), nor ten
   timestamps together (115 percent at 2,000 columns, 91 at 10,000, where the solver reaches its cap).
   About 45 lines; 67 s and 1 GB when a thousand columns are re-admitted, per fold in cross validation.
3. As ruled, plus `lm.fit` on doubtful designs. When a dropped column is doubtful the design is made dense
   and given to `lm.fit`, while it fits in 1 GB. A user gets lm's own answer there, its dependence on
   column order included (43 or 74 percent on the ten timestamps); it takes 8 s at 3,000 rows and 64 s at
   20,000 with 2,000 columns; above 2,000 columns nothing changes. About 25 lines.
4. Dense when small; larger, as ruled with zero-filled columns shifted. A design whose dense form
   `lm.fit` handles in about a second and 200 MB gets `lm.fit`'s sigma, whatever its storage; a larger one
   gets the sparse regression after a zero-filled column beside its missing indicator is re-centered. A
   user gets exactly what a dense caller and 0.9-34 get on small designs, lm's order dependence included,
   and the full fit on larger timestamp designs. Left, past the bound and so above 1,000 rows: a near-copy
   with signal (132 to 143 percent), two zero-filled times sharing one missing pattern with no indicator
   (120 percent, measured with the routine forced on a small design), rows weighted 1e12 (35 percent).
   About 40 lines; a fit under the bound takes up to 1.2 s where the sparse routine takes 0.01 s, so cross
   validation's 200 folds take up to 4 minutes against 2 s, as a dense design of that size does today.
5. Option 4, plus the recheck and the tighter solver past the bound. Small designs as in 4; larger ones
   get the full fit on every design measured but rows weighted 1e12 (16 percent at 2,000 columns, 5 at
   the solver's cap). About 90 lines; the recheck holds at most 64 columns, 4.9 s and 250 MB at worst
   measured; the solver runs to its cap on slow well-posed designs, with 1e-8 left.

| option | small designs | 2,000 columns | 10,000 columns | cost |
|---|---|---|---|---|
| 1 | +235 to +269 | timestamps +174, near-copy +132 | +168, +143 | none |
| 2 | 0; +235 at five minutes | +115, 0 | +91, 0 | 45 lines; 67 s, 1 GB worst |
| 3 | lm's | lm's (+43 or +74), 0; 8 s | +168, +143 | 25 lines; 8 to 64 s when triggered |
| 4 | lm's, to the bit | 0, +132 | 0, +143 | 40 lines; 1.2 s a fit under the bound |
| 5 | lm's, to the bit | 0, 0 | 0, 0 | 90 lines; 4.9 s, 250 MB worst |

"lm's" is the full fit in all but two of the small designs run: the five-minute timestamp with its
indicator first (235 percent) and the start and end times (107 percent), where a dense caller gets the
same today. Options 1, 2, 4 and 5 do not depend on column order past the bound; 3 does where it falls back.

How 4 and 5 sit with the rulings. dec-B403's rule, the dense equivalent's sigma, is met exactly on small
designs, which is where the sigma matters. Its words, "just use the sparse regression", are not followed
there: a small sparse design is made dense for this one estimate (at most 6e6 numbers, 48 MB). dec-B422
stands: a factor enters as indicators, and a default data frame with a factor gets `lm.fit` on the
indicator design, 0.9-34's value.

Recommended: 5. A user expects the dense equivalent's sigma, and 4 and 5 give exactly that wherever the
fit can show the difference. Past that size lm is neither affordable nor stable, and the correct target
is the least-squares fit; 5 reaches it on everything measured but weights spanning twelve orders of
magnitude, at a bounded cost. Option 4 is the smaller build and leaves the nearly dependent designs past
the bound as dec-B421 has them. The steps marked OC1 are written for 5; under 4, E.5 goes and F.4's
tolerance is 1e-6.

### 2. How far from the reference may LSQR be before a design comes back to the maintainer?

Background. dec-B426 lets LSQR ship on the condition that a badly conditioned design "that reaches the cap
with a large error comes back to the maintainer before landing". "Large" has no number yet, and the build's
verdict turns on it. At n 1000, sd(y) in place of the linear estimate (15 to 75 percent above it in the
sensitivity note's settings) moved RMSE and coverage by less than a reseed does, and the posterior mean of
sigma by 1.6 points of the truth where a reseed moves it 0.8; at n 5000 nothing exceeded a reseed. LSQR
runs only where n is above about 2,200. What is already known: the crossed-factor design at 10,000 columns
and 11,200 rows is 11.8 percent high from its degrees-of-freedom count alone (8.0 and 4.9 percent at the
two larger sizes). The cap is reached at 1e-10 by slow well-posed designs with under 1e-6 left, and by
weights spanning twelve orders of magnitude with 2 to 21 percent left. The timestamp family reaches it 86
to 165 percent high unless its columns are shifted (Open calls 1, options 4 and 5), and then converges in
about 100 iterations onto the full fit, inside any threshold; under options 1 to 3 it comes back under any.

Options. (a) 10 percent: under the smallest error the note measured (15 percent); the one crossed design
comes back whatever the solver does. (b) 1 percent, the note's "to a percent or better": every design with
dependencies the structure does not show comes back (all six crossed designs), so the review is of the
count more than of the solver. (c) A factor of 2, the largest error the note found close to a reseed at n
1000 (sigma 1.8 points): the crossed designs and the weights pass, and what comes back is only what Open
calls 1 leaves in the band (168 percent on the timestamps under its options 1 to 3).

Recommended: (a).

### 3. Does a design that takes sd(y) under the 10 percent rule get a warning?

Background. dec-B426 says LSQR runs "with no message" and that above the cutoff a design with under 10
percent residual degrees of freedom "takes sd(y)"; it does not say whether that is announced. Every other
fallback to sd(y) raises dbartsSigmaFallbackWarning (a dense design with p >= n; a sparse one with no
residual rank at or below the cutoff will too), and xbart reads that warning to choose its per-fold route.
A wide sparse design lands here. So does a dense frame with a high-cardinality factor, whatever this call
decides: at n 1000 with a factor of 1,000 levels each present once, `dbarts` and `bart` give 1.145 today
with no warning (the factor as codes) and `bartBT`, which expands factors, already gives sd(y) 1.533 with
the warning; after dec-B422 every door does, because no residual degrees of freedom are left. At n 5000
with 4,600 levels 399 structural degrees of freedom are left, under the 500 line: the 10 percent case.

Options. (a) Warn with the existing class, the message naming the rule and `sigest`: sd(y) is always
announced. Cost: a warning on every such fit until sigest is given, as wide dense fits have, now also on
default fits with an identifier-like factor that today fit silently. (b) Silent: about 8 lines (the routine
reports its route so xbart can still choose). Cost: a wide sparse fit warns at 2,000 columns and not at
2,001, and `sigest` reports sd(y) unmarked. (c) Warn only with no residual degrees of freedom, silent
between 0 and 10 percent: about 10 lines; matches the dense path where it warns and hides the new rule.

Recommended: (a).

### 4. Should the 10 percent rule stop at the cutoff?

Background. By dec-B426 the rule applies above m 2,000 only. A design with 3 percent residual degrees of
freedom gets the exact estimate at m 2,000 (relative sd about 9 percent on 60 degrees of freedom) and
sd(y) at m 2,001 (15 to 164 percent above the exact estimate in the note's row). On nearly dependent
designs the two sides of the cutoff can also disagree with each other: as ruled, a timestamp design at
the cutoff gets 1.33 from the exact routine, and LSQR, which a design one column wider would get, gives
0.93 to 1.32 by its tolerance, where the full fit is 0.48. Under Open calls 1's options 4 and 5 both
sides give 0.48.

Options. (a) As ruled; the step at the cutoff stays; nothing to build. (b) The rule for every sparse design
past the dense bound: about 4 lines, and such a design then differs from its dense equivalent below the
cutoff, against dec-B403. (c) The rule for dense designs too: about 10 lines, posterior-changing for every
design with few residual degrees of freedom, a wider re-record, and its own evidence (the note compared
sd(y) with the linear estimate at 19 degrees of freedom only at n 200).

Recommended: (a) for this slice, and (c) as its own TODO item if wanted; it does not block the build.

### 5. May the starting sigma's last digit follow the BLAS thread count?

Background. A one-digit change in the starting sigma changes a seeded fit's draws from the first. `lm.fit`
gives the same bits at any thread count and different ones under another BLAS (R's against Accelerate, run
here). The exact routine calls LAPACK's pivoted Cholesky and differs under another BLAS and, under
Accelerate, between 1 and 2 threads (5 of 14 designs in the third critique); two runs at one setting agree.
LSQR's bits follow neither. As ruled this reaches every sparse design up to 2,000 columns and every default
data frame with a factor level at or under 20 percent density, which `lm.fit` serves today. Under Open
calls 1's options 4 and 5 it reaches only designs past the dense bound and up to 2,000 columns.

Options. (a) Accept; the help says a seeded fit is reproducible at a fixed BLAS and thread count. (b) A
pivoted Cholesky written in the package: in R about 40 lines and, by estimate, tens of seconds at 2,000
columns; in C out of this slice. (c) LSQR for every design past the dense bound: about 5 lines; the exact
rank goes, so dependencies the structure does not show cost degrees of freedom (0.5 percent for 50 in
4550).

Recommended: (a), with Open calls 1's option 4 or 5 keeping small designs on `lm.fit`.

## Estimate

Implementer about three days (the size route, front end with the shift, exact routine with recheck, LSQR
and the design builder two; tests and docs one). Condition (a) about six hours of machine time at two
jobs: the exact reference at m 10,000 is two minutes and 1.8 GB a design; a tight LSQR reference is up to
30 minutes a design at m 50,000. Gates about two hours, the re-record minutes. Review with mutants half a
day, and two fix rounds.
