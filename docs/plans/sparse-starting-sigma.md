# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403); revised after its first blind critique; revised 2026-10-10 for
dec-B421, dec-B422 and dec-B426 and for the sigest-beside-fixed landing (dec-B413); revised again from its
second blind critique. Not built. Four calls are open ([Open calls](#open-calls)); the steps marked "OC1"
build the recommendation of the first and wait on its ruling.

agent: one sonnet implementer (R only, no engine code); one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING, for a continuous response under a chisq residual prior with no sigest, in two
cases: a caller's sparse x (dec-B403), and any fit that keeps a factor as a categorical predictor, the
default, dense or sparse (dec-B422). An indicator expansion R built sparse (dec-B370) moves by rounding,
and by more where a column sits in the rank band of Open calls 1. NEUTRAL for a dense design with no
categorical column, every fixed-unit family, every fit given sigest, and the draws of a fit under a fixed
residual prior.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~1100 lines (R/utility.R ~310 net of the QR's removal, R/spec.R ~30, R/xbart.R ~35, tinytest ~430,
benchmarks/R ~180, man ~40, docs ~90, NEWS ~2, MANIFEST one row).

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
- dec-B421: "Band is fine." The exact routine's rank tolerance may differ from lm's on near-copies. It was
  shown 0.01 to 4.7 percent, and 14 percent; what was measured since is Open calls 1.
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
  (`sparse.from.indicators`, set by [`makeIndicatorModelMatrix`](../../R/utility.R)) goes to
  [`sparseResidualStandardError`](../../R/utility.R), Matrix's sparse QR; anything else to
  [`residualStandardError`](../../R/utility.R) (`lm.fit`). [`floorSigmaEstimate`](../../R/utility.R) turns a
  non-finite result into sd(y) with dbartsSigmaFallbackWarning: a dense design with p >= n warns; one with
  2 residual degrees of freedom of 40 rows gives an estimate and no warning.
- A factor kept categorical enters as its level codes, ordered or not: a 40-level factor and a numeric
  column (n 300) start sigma at 1.088 by default (sd(y) 1.095) and 0.499 under `factors = "indicators"`.
  Every door reads one estimate: `dbarts` (bcf is its `forests` argument), `bart`, and `dbartsData` with
  `dbartsSpec` give a frame and the same data as a matrix of codes carrying `varTypes` one value (0.5538
  on a 4-level and an ordered 3-level factor, where lm on indicators gives 0.5255); a multinomial fit makes
  none; `bartBT` refuses a matrix of codes. `rbart_vi` builds through `dbarts` (not run).
- The frame rebuilt from a dense container's column list gives, through
  [`makeIndicatorModelMatrix`](../../R/utility.R), the design values of the `factors = "indicators"` fit.
  That builder stores an indicator column sparse at a density of at most
  [`sparseIndicatorDensity`](../../R/utility.R), so a 5-level factor at n 300 already takes the sparse
  routine there. A missing factor value is imputed with each level's frequency; the block still sums to one.
- An infinite stored entry in a dgCMatrix is accepted; a dense one is refused, naming the column.
- A fixed residual prior with no sigest still runs the estimate ([`resolveSamplerSpec`](../../R/spec.R);
  [`xbart`](../../R/xbart.R) before it resolves the prior), and a sparse x warns about an estimate nothing
  uses. Draws are identical whatever the slot holds. Since dec-B413 an agreeing sigest is accepted with a
  message and sits in the slot; a differing one is refused (dec-B425). The slot is read later (a `setModel`
  to a chisq prior, `rbart_vi`'s start, `bart`'s `sigest`); `dbartsSpec` documents that a value the data
  object carries survives, and stan4bart writes one under `fixed(1)`.
- xbart on a sparse x takes the per-fold route "marginal" (each fold its own sd); a dense design "linear",
  where `foldData` fits [`sigmaDesignMatrix`](../../R/utility.R) of `data@x`, a factor's codes included.

The exact routine ([Algorithm](#algorithm) D), built from this plan's text and run 2026-10-10 on R's
reference BLAS:

- Agreement with `lm.fit` on the dense design: 400 random designs (n 30 to 2000, p 5 to 1530, numeric,
  one-hot and mixed, a third with dependent columns, a third weighted with zero weights, a third with an
  offset; 143 with p >= n). Rank equal to lm.fit's on 398; 106 with no estimate on both; of the 292 with
  one, sigma within 9.1e-12 relative above 30 residual degrees of freedom and 2.6e-8 at 30 or fewer. The
  other 2 are wide (lm.fit keeps columns with relative pivots near 1e-16), 0.26 and 0.64 percent apart.
  With the recheck of D.5 the sweep is the same on every design (no column re-admitted; 14.2 s against
  12.7 s in all).
- Cost at m 2,000: 1.0 to 1.1 s and 90 to 160 MB over R's own heap on a narrow design (n 2e4, 1 percent);
  0.8 s on a wide one of full row rank; 3.9 s and 170 MB on a wide rank-deficient one (n 2000, p 3600,
  rank 901). The recheck adds nothing when no column is dropped, 0.3 s for 20 dropped columns, 3.0 s for
  1,000 (3.8 s against 0.8 s, 220 MB) and 2.3 s on the wide rank-deficient design. For scale: dense
  `lm.fit` at n 2e4, p 2000 took 64 s, the sparse QR 65 s there and 157 s on one-hot 20 x 50 at n 1e4.
- Bits: the sparse products, the crossproduct and LSQR give the same bits under R's BLAS and Accelerate and
  at 1 and 2 threads. The exact routine's sigma differs in its last digit between the two BLAS and, under
  Accelerate, between 1 and 2 threads; `lm.fit` differs between the BLAS and not with the thread count.
- Reach. The exact routine serves a caller's sparse x; an indicators-route design with a sparse-built
  column, which today takes the sparse QR at lm's tolerance; and, through dec-B422 and the shared builder,
  every default dense frame with a factor level at or under 20 percent density.
- The band, as planned without D.5: a column is dropped when its residual after the kept columns is below
  1e-5 of its norm (pivot 1e-10 on the unit-diagonal crossproduct; an exact dependency's pivot measured up
  to 1.5e-13, above the 1e-14 that lm's 1e-7 would need). What that costs where such a column carries the
  signal is in Open calls 1. Weights do the same through the crossproduct: one row at 1e12 times the
  others loses 4 of 21 ranks (175 percent high), three rows 388 percent, where `lm.wfit` keeps full rank
  through 1e15. In the other direction a centered dense column whose spread is under 1e-7 of its mean is
  dropped by lm and kept here (53 to 61 percent below lm where it carries the signal).

LSQR, built from this plan's text (base R and Matrix):

- Two runs give the same bits. Weighted fits (uniform, 20 percent zeros, 0/1, skewed, a few at 1e-12, with
  an offset) agree with `lm.wfit` to 7e-11 in 11 to 50 iterations at n 6000, p 400; 1 percent of rows at
  weight 1e6 takes 785 iterations. A full one-hot block left uncounted costs a degree of freedom (9e-5),
  and is left uncounted when the row sums are taken after a dense-backed level is centered (1.67e-4).
- At tolerance 1e-6 a stop on tolerance is no evidence of accuracy: the band's designs stop at iteration 4
  to 7, 136 to 258 percent high (Open calls 1, with what 1e-10 does).
- Slow designs: 20-column groups correlated at 0.99 (n 2240, p 2000) take 704 iterations at 1e-6 (3e-6 from
  exact) and reach the cap at 1e-10 (1e-8 from exact); a tolerance of 1e-10 run to the end takes 1,612.
- Run past convergence on a design with an exact dependency, the iterate drifts: 6 times too high at 100
  iterations and 31 times at 1,000. The standard condition-estimate stop at 1e10 ends it at iteration 27,
  2e-14 from lm.
- A stop on "the response is reproduced" taken against the norm of the response stops a response of mean
  1e10 at the first iteration, 104 percent high; against the centered norm it runs to tolerance.
- Costs: [Cost and memory](../design/starting-sigma-sensitivity.md#cost-and-memory).

What moves (every call of [`estimateSigmaFromLinearModel`](../../R/utility.R) logged through a load hook,
beside the value `lm.fit` gives on the indicator design):

- equivalence.R, quick mode, 132 calls in 44 of 55 scenarios. Twelve move:
  - seven on a caller's sparse source, sd(y) to the linear estimate (2.959 to 1.370, 1.670 to 0.271, 3.609
    to 2.122, 3.664 to 2.265, 3.941 to 1.833, 3.483 to 2.125, and 4.124 to 2.386 on all rows):
    ["sparse <- list("](../../benchmarks/R/equivalence.R), ["mixedmatrix <- list("](../../benchmarks/R/equivalence.R),
    ["sparsefactor <- list("](../../benchmarks/R/equivalence.R), ["testswap <- list("](../../benchmarks/R/equivalence.R),
    ["leaffactormixed <- list("](../../benchmarks/R/equivalence.R), ["factorpartial <- list("](../../benchmarks/R/equivalence.R),
    ["xbartmixed <- list("](../../benchmarks/R/equivalence.R);
  - four on a dense frame with a factor, codes to indicators (1.9855 to 1.9870, 1.863 to 1.830, 2.469 to
    2.389, 2.234 to 2.330): ["categorical <- list("](../../benchmarks/R/equivalence.R),
    ["leaffactor <- list("](../../benchmarks/R/equivalence.R), ["nafactor <- list("](../../benchmarks/R/equivalence.R),
    ["ordfactor <- list("](../../benchmarks/R/equivalence.R);
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

`sparseResidualStandardError(y, x, weights, offset)` keeps its name; the cutoff (2000), LSQR's cap (1000)
and LSQR's tolerance are formals with their shipped defaults, so a test and the build's reference run reach
either route at any setting.

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

B. Front end, shared by both routines, in this order.

1. Inf: an infinite entry is an error that [`estimateStartingSigma`](../../R/spec.R) reports as for a dense
   design, naming the column: [`nonFinitePredictorNames`](../../R/spec.R) learns sparse sources.
2. Rows: drop rows with a missing response, weight or offset and rows of weight 0, as `lm.wfit` does;
   z = y - offset; n the rows kept.
3. Constants: a column whose entries over the kept rows are all equal (implicit zeros included; an exact
   comparison of values) is dropped. p is the number of columns left and m = min(n, p + 1).
4. Blocks, on the values as they stand: the columns of one factor term (`indicator.term`) are a full block
   when at least two are left and they sum to one in every kept row (within 1e-8). b is their number.
5. Centering: each dense-backed column and each column stored in every kept row is centered at its weighted
   mean over the kept rows. OC1: one whose centered norm is under 1e-7 of its uncentered norm is dropped,
   lm's verdict against the intercept, and leaves p and m. Then every column is divided by its largest
   absolute entry.

C. Route. m <= 2000: the exact routine (D). Otherwise the structural residual degrees of freedom are
df = n - 1 - p + b; when df < 0.1 n there is no estimate (F); else LSQR (E). m is the column count after the
expansion, so a dense frame with a 5,000-level factor can be an LSQR design.

D. The exact routine.

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
   through R (the two repeats take the error from 1e-6 to rounding), in chunks of at most 8e6 entries. One
   whose residual is above 1e-7 of its unit norm, lm's tolerance, is a candidate. Candidates are
   orthogonalized in order of size, twice, and one still above 1e-7 is accepted: it joins the fit (narrow:
   e loses its component along it after one more projection through R; wide: B v joins the basis for an
   accepted row direction v) and r grows by one. At most 2.5e7 / n are accepted (200 MB); the rest stay
   dropped. An exact dependency leaves a residual at rounding and is never a candidate.

E. LSQR (Paige and Saunders), in R in R/utility.R, on base R and Matrix's sparse products alone.

1. Operator. With h = sqrt(w), s = sum(w), mu_j the weighted mean of column j and d_j the reciprocal of its
   weighted centered norm, A = diag(h) [1 / sqrt(s), (X - 1 mu') diag(d)], applied as two sparse products a
   step (`X %*% v`, `Matrix::crossprod(X, u)`) with the centering as a rank-one correction. Unit column
   norms are the only preconditioning. A column whose centered sum of squares comes out 0 takes d_j = 0:
   it adds nothing and stays counted, so p, m and the route are those of B and C.
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
8. Degrees of freedom are C's count, never the iteration's. A dependency the structure does not show is
   counted as a fitted column and overestimates by sqrt((df + d) / df). The estimate is never below the
   fit of every column, but against `lm.fit` or the exact routine, which drop columns by tolerance, it
   can be lower.
9. Determinism. No random number, fixed constants, one operation order; the bits held across two BLAS and
   two thread counts (Context). The exact routine's last digit follows the BLAS and, under a threaded
   one, its thread count; the reference build's baselines are recorded on R's own.

F. No estimate (D with no residual rank, C under 10 percent) returns NA, which
[`floorSigmaEstimate`](../../R/utility.R) takes to sd(y - offset) with dbartsSigmaFallbackWarning, as the
dense path does. Whether the 10 percent case warns is Open calls 3; the steps build its recommended option,
and the other is about 8 lines. An allocation failure is an error, reported as today.

## Change

1. R/utility.R: `startingSigmaDesign` (A); [`sparseDesignMatrix`](../../R/utility.R) gains the wrap and the
   two expansions; [`sparseResidualStandardError`](../../R/utility.R) rewritten as B to F with the exact
   routine, its recheck and LSQR as internal functions (the batching, the refactor loop and the `grepl`
   muffler of Matrix's warning go); [`estimateSigmaFromLinearModel`](../../R/utility.R) loses the fallback
   branch and its warning and routes by what the design builder returned;
   [`makeIndicatorModelMatrix`](../../R/utility.R) records `indicator.term` and no longer sets
   `sparse.from.indicators`.
2. R/spec.R: [`nonFinitePredictorNames`](../../R/spec.R) for sparse sources (the dense list and the CSC
   block's stored entries, named through the container's map, or positions for a bare dgCMatrix).
3. R/xbart.R: the design is built once per chunk by `startingSigmaDesign` in
   [`xbartRunUnits`](../../R/xbart.R), for dense sources too, so a fold of a dense frame with a factor fits
   indicators; `foldData` fits its training rows by the same routes, each fold by its own m and its own
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
calls 2 (recommended: 10 percent).

(a) LSQR on badly conditioned realistic designs. A tracked script, benchmarks/R/starting-sigma-lsqr.R, calls
the implemented routine; two jobs at a time.

- Crossed factors, defined by occupancy so the size is what is stated: two factors of L levels, K of the
  L^2 cells occupied, every occupied cell one row and the rest of the rows spread over them with Zipf
  frequencies; columns are the two main effects and the occupied cells as one-hot columns of a dgCMatrix,
  so no dependency is structural. (L, K) = (150, 9,699), (200, 19,599), (300, 49,399) give m 10,000,
  20,000 and 50,000 exactly (built and counted at n = 1.5 m and n = 1.12 m). The reference is the
  cell-means fit. The count alone is 3.0, 2.0 and 1.2 percent high at n = 1.5 m and 11.8, 8.0 and 4.9 at
  n = 1.12 m, so the first of those is past 10 percent before the solver runs.
- At the same six (m, n): word counts (term frequencies Zipf, document lengths lognormal, raw counts);
  correlated numeric columns (groups of 20 sharing one sparse pattern at 5 percent, pairwise correlation
  0.99, beside 200 fully stored columns correlated at 0.999); a sparse numeric design with ten zero-filled
  timestamp columns (spreads from a minute to a month) each beside its missing indicator and carrying
  signal; the same design with 50 near-copies at 1e-8 to 1e-3, half carrying signal.
- At m 10,000 each also with skewed weights (a cubed exponential) and with 1 percent of rows at weight 1e6.
- Reference: the cell means for the crossed family. Elsewhere the exact routine at m 10,000 (cutoff
  raised, one job at a time), checked against `lm.fit` on a twin of each family at m 2,000; above m
  10,000 LSQR at tolerance 1e-12 and a cap of 20,000, admitted only if it stops on tolerance and if the
  same setting agrees with the exact routine within 1e-6 on that family at m 10,000.
- Recorded per run: iterations, seconds, stop reason, the residual sum of squares against the reference's,
  and sigma against the reference's.
- PASS: every run's sigma is within T of its reference, in either direction, whatever its stop reason.
- BACK TO THE MAINTAINER, before landing: any run further than T from its reference, with the table and
  three priced options (a higher cap or tighter tolerance, sd(y) there, a block preconditioner); a design
  with no admitted reference, as unmeasured.
- DEFECT, stop: a run at the shipped tolerance whose residual sum of squares is more than 1e-6 below that
  of the tolerance-1e-12 run on the same design (the iteration cannot lose ground by running longer); or a
  family whose m 2,000 twin is further than T from `lm.fit` with OC1 built.
- A run more than 1 percent from its reference, from the solver or the count, is named in the Landing note.

(b) Weighted fits against the exact routine, in tinytest with LSQR forced (cutoff 0) on designs under the
real cutoff: sparse numeric, an indicators-route design with a full block that has a level above 20 percent
density, a frame with a sparseFactor; under uniform weights, 20 percent zero weights, 0/1 fold weights,
skewed weights, a few rows at 1e-12, and weights with an offset. Each weighted case is read beside its
unweighted twin: both within 1e-6 relative with equal degrees of freedom, PASS; the weighted one past it and
the twin within, FAIL as a defect in how weights enter; the twin past it, the design is wrong for this test
and is replaced (these three are within 1e-10 unweighted). The 1e6-weight case is (a)'s.

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
  perturbation (0, 1e-9, 5e-8, 2e-7, 1e-6, 1e-5, 3e-5, 1e-4, 1e-3) to the dense fit within 1e-8 with lm.fit's
  rank (OC1; under option (a) there, to the nearest of the dense fits with and without each copy).
- [test-sampler-splitProbabilities.R](../../inst/tinytest/test-sampler-splitProbabilities.R): a seeded
  comparison of split counts on a design with a factor fails at its seed. The implementer reports the
  counts over five seeds; if the comparison holds in distribution the seed moves, and if not it is a stop.
  test-data-code-channel.R passes unchanged and is the pin for A.1.

New, in test-sparse-starting-sigma.R (under 8 s), against `lm.fit` on the dense indicator design, missing
values imputed as the estimate imputes them, within 1e-10 unless said:

- dec-B422, dense: a frame with a 5-level factor and a numeric column, default against
  `factors = "indicators"` and against its matrix of codes: `sigest` identical; equal to
  `summary(lm(y ~ f + x))$sigma` where nothing is missing; the same with the factor ordered, with two
  levels, with 40 (built sparse), and with a missing value against the imputed design; `data@x` and its
  `varTypes` identical to the fit given sigest.
- dec-B422, sparse: a frame with a sparseVector column and a dense factor; a sparseFactor with its first
  and a middle level as reference, and a middle reference with missing values.
- dec-B403: a bare dgCMatrix through `dbarts`, `bart` and `xbart` with no warning (counted, Gate hygiene);
  one-hot 3 x 40 plus three duplicated columns with zero weights and an offset; NA entries in a dgCMatrix; a
  fully stored timestamp-like column carrying signal, weighted and not; a column scaled by 1e160; p > n
  with rank below n - 1; p > n of full rank (NA, and through `dbarts` sd(y) with one
  dbartsSigmaFallbackWarning); that design with one row at weight 1e-14 (NA).
- OC1: the zero-filled timestamp with its indicator, spread an hour, on a dgCMatrix and in a dense frame
  with a 5-level factor, within 1e-8; a near-copy carrying the signal at 1e-6; a centered column with
  spread 3e-8 of its mean (lm's rank); one row, and 1 percent of rows, at weight 1e12, within 1e-8; a wide
  design with a heavy row.
- Route: at n 2600 and 1 percent, p 1999 takes the exact routine and p 2000 LSQR (read from the route the
  routine reports), and the p 2000 design with the cutoff raised gives the exact value within 1e-6 of
  LSQR's; n 2200, p 2100: no estimate, sd(y), the warning of Open calls 3. With the cutoff formal: n 40,
  p 100 at cutoff 45 is exact (m is n); 20 live and 30 constant columns at n 200, cutoff 25, is exact (m
  counts live columns).
- LSQR forced: condition (b); a full block with a level above 20 percent density counted (degrees of
  freedom equal to the exact rank's); a block under a drop pattern that removes a present level not
  counted; a sparseFactor whose reference level has only zero-weight rows counted as a block; a constant
  fully stored column under weights dropped; a zero response; a response of mean 1e10 within 1e-6 of the
  exact routine; a response the design reproduces (sigma under 1e-8 of sd(y)); a one-hot design at
  tolerance 0 stops on the condition estimate within 1e-6; the cap at 3 iterations returns a larger sigma
  than the converged one and raises nothing; two calls `identical`, for each route.
- An Inf entry in a dgCMatrix and in a sparseVector column: the dense path's error, naming the column.
- xbart: on a sparse frame and on a dense frame with a factor, each fold's sigma equals `lm.fit`'s on the
  indicator design for the same rows; on n 3300, p 2100 at 1 percent each fold takes LSQR (the route
  traced).
- A fixed residual prior, dense and sparse: `estimateStartingSigma` is not called; with nothing given the
  slot and `sigest` equal the fixed sigma; an agreeing sigest and a sigma the data object carried stay;
  draws identical to the build before; then `setModel` to a chisq prior draws finite sigmas.
  test-sigest-fixed-agree.R passes unchanged.

The 37 files that reach a factor design must otherwise pass unchanged; a failure there that is not a pinned
sigest, slot or seeded draw of a default fit with a factor is a stop.

Reviewer's mutants, each of which must fail a test: the expansion removed (codes), for a container, for a
matrix of codes, and in `foldData` alone; an ordered factor left as codes; a reference-level indicator
added to a sparseFactor; `indicator.term` not recorded by the caller's-sparse builder; imputation before
the expansion; cutoff 2000 to 20000, and `<=` to `<`; m taken as p + 1, and before constants are dropped;
the 10 percent rule removed, and taken against p; b forced to 0, counted without the row-sum check, and
counted after centering; the constant test replaced by a sum-of-squares one; weights left out of the
operator, of the means, of the residual; the LSQR tolerance at 1e-2; the reproduced-response stop against
the uncentered norm; the condition stop removed; sigma from the recurrence's residual norm; the intercept
left out of the operator; zero-weight rows counted in n; exact tolerance 1e-10 to 1e-16 and to 1e-6; the
recheck removed, its threshold at 1e-5, its repeat projections removed; lm's intercept verdict removed;
centering removed, and at the unweighted mean; the max-abs step removed; the equilibration removed; the
wide basis without D, and from `qr` at its default tolerance; `n - r` replaced by `n - p`; the Inf check
removed; the fixed-prior skip leaving the slot NA, and overwriting a carried value.

## Baselines

- Current: equivalence-e4faed5c, bcf-equivalence-1b7d730c, multinomial-equivalence-80b1c8d4
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)).
- Moves, by class: POSTERIOR-CHANGING the seven caller-sparse and the four dense-factor scenarios (Context);
  SHIFTING wideFactorIndicators (rounding; bitwise if the routine happens to round to the QR's value).
  NEUTRAL the other 43, bcf's 15 and multinomial's 11, the four snapshot files and every exact gate. Any
  other mover is a defect: stop.
- Re-record the twelve on the reference build (`--preclean --configure-args=--enable-reference-build`,
  `EQUIVALENCE_CORES=2`), merged into a copy of e4faed5c in its scenario order and named after the slice's
  code commit; e4faed5c demoted to historical.
- Partition against e4faed5c in z mode: 43 of 55 identical and the movers exactly the twelve, with no |z|
  above 4 on wideFactorIndicators. The dense-factor anchors move 0.08 to 4 percent and the sparse ones fall
  by factors of 0.16 to 0.62, so |z| above 4 is expected among the seven and possible among the four. The
  merged file reproduces 55 of 55 under `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17): the change is the value a prior is calibrated against, not the sampler; the
  identity is agreement with `lm.fit` on the indicator design (the tinytest pins, and the 400-design sweep
  rerun against the implemented function with factor columns added to a third of its designs).

## Gates

On the slice tip against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing):

- tests/cpp unchanged and green (no C++ touched); the full tinytest suite green with the edits above; the
  four seeded-drift snapshot files on the reference build pass unchanged.
- The equivalence trio as Baselines, bcf and multinomial bitwise against their current files;
  exact-gates.yaml's list in quick mode all pass, output as before.
- `R CMD check --as-cran` from a clean tarball: no new NOTE; lintr, air, rc-codoc, win-drift,
  doc-freshness. Sanitizers are not owed: no compiled code changes.
- Speed, same machine, within 1.5x: the exact routine at n 2e4, p 1999, 1 percent (1.1 s), with 1,000
  dependent columns (3.8 s) and on the wide design n 2000, p 3600, rank 901 (6.2 s); LSQR at m 1e4, 1
  percent (0.15 s at tolerance 1e-6; re-timed at the shipped one) and m 5e4; xbart's 200 per-fold estimates
  at n 1e4, p 1000, 1 percent (28 s). No bench-sampler.R compare (no hot path).
- Build conditions (a) and (b) with their verdicts.

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd and man/xbart.Rd: a sparse x takes the same
  estimate as its dense equivalent, exact while the smaller of its row and column counts is at most 2,000
  and iterative above that, where under 10 percent residual degrees of freedom gives the marginal standard
  deviation; supply `sigest` to skip it. A factor enters as indicator columns whatever `factors` says.
  Under a fixed residual prior no estimate is made. man/xbart.Rd gains the sentence the other pages carry
  on an agreeing `sigest` beside a fixed prior; man/dbartsSpec.Rd's "an unset value is still estimated"
  gains "under a chisq prior".
- man/bart.Rd's warning-class paragraph drops the class; `sigest` under a fixed prior is the fixed sigma.
  man/sparseFactor.Rd's starting-sigma paragraph becomes one sentence: a sparse column leaves the default
  starting sigma as the same column stored dense gives it. inst/NEWS.Rd: the class leaves the
  warning-class list (never released; dec-B422 restores 0.9-34's value, so no entry for either).
- docs/design: sparse-columns.md's [R surface](../design/sparse-columns.md#r-surface) gets a dated paragraph
  with the rule and Algorithm in brief; starting-sigma-sensitivity.md a dated section with condition (a)'s
  table and verdict; error-style.md drops the class; memory-footprint.md's starting-sigma row gains the
  sparse routes (up to about 220 MB at m 2,000 for the exact routine; vectors only for LSQR), per worker
  under xbart.
- At landing: TODO's item goes; a ledger entry for the calls below; this plan's Status and Landing note;
  the MANIFEST row.

## Steps

0. Open calls 1 is ruled, or the implementer is told which option to build. Under its option (a), B.5's
   lm verdict, D.5 and the OC1 tests go, E.4's tolerance is 1e-6, and the band's cases are pinned as they
   fall.
1. Change 1 to 5 with Algorithm A to F; the tinytest edits and new tests; the suite green against
   `R CMD INSTALL -l <lib> .`.
2. The 400-design sweep against the implemented function; build condition (a), its table and verdict
   written into the note. A verdict other than PASS stops the slice here.
3. The speed points; help, NEWS and docs.
4. After review: the re-record, MANIFEST row and partition, in their own commit.

## Stop conditions

Stop and report when: a build condition's verdict is not PASS; the diff passes ~1500 lines; an equivalence
scenario other than the twelve, a snapshot or an exact gate moves, or wideFactorIndicators shows a |z|
above 4; the sweep shows a rank different from lm.fit's on a design other than its two wide ones, or a
sigma off by more than 1e-10 above 30 residual degrees of freedom; a speed point is past 1.5x; a tinytest
outside the four files fails; the change needs engine, bridge or C API code.

Whichever lands second of this and any other slice re-recording equivalence.R re-records against the
other's file and partitions against it.

## Calls made

- LSQR is written in R: its time is Matrix's two compiled sparse products a step, and the loop can be
  interrupted. In C: about 250 lines, a bridge entry and its Windows twin, sanitizers owed, no faster.
- At the cap the iterate is used silently, read from dec-B426's "with no message" and from condition (a)
  being the check on it. A warning at the cap would be about 6 lines and a test.
- The intercept is a column of the operator and every column is centered and scaled inside it; no other
  preconditioner. Four stops the prototype lacks, each from a run in Context: a reproduced response
  (against the centered norm), a zero norm at the start, the condition estimate, and a tolerance formal.
- A factor is expanded wherever the estimate meets one, a matrix of codes included, by the indicators
  route's own builder: one sigest across the doors, at the cost that a default fit with a factor level at
  or under 20 percent density takes the sparse routine, as an "indicators" fit does today.
- A caller's sparse source builds every factor without one level and records each indicator's factor; full
  blocks are found by that record and a row-sum check on the raw values, not by widths read off `drop`.
- LSQR takes the exact routine's front end, a fold subsetting its rows (2.8 iterations' products). m is
  taken after the expansion and constant drops, on the rows kept; each xbart fold routes by its own m.
- Under a fixed residual prior no estimate is made and only an empty slot is filled. The alternative,
  overwriting whatever the slot holds, gives one value on every path but breaks what `dbartsSpec`
  documents (a carried value survives) and what stan4bart writes; it is one line either way.
- The exact routine's last digit follows the BLAS and a threaded BLAS's thread count; `lm.fit`'s follows
  the BLAS only. Accepted: the alternative is a hand-written pivoted Cholesky, too slow in R and out of
  this slice in C.
- Extreme weights (a row at 1e12 times the rest) are the band's mechanism, and `lm.wfit` is sound there.
  With OC1 built they are pinned at 1e12; under its option (a) they are left unpinned.
- Condition (a)'s script is tracked in benchmarks/R; its verdicts are two-sided, and a miss from the
  count alone is treated as the solver's would be.
- Carried from the first revision: one routine for every sparse design; the smaller Gram side with the
  intercept as a column; centering, max-abs and unit-norm scaling; the wide side equilibrated; a sparse
  Inf refused; the sparse warning class retired.

## Open calls

Each is the maintainer's; none is settled here.

1. How closely must the sparse routine follow lm when a column is almost, but not exactly, a combination
   of others?
   Background. The starting sigma is the residual sd of a linear regression of the response on the
   predictors. The dense path (`lm.fit`) drops a column as redundant when less than 1e-7 of its length is
   left after the columns before it. The planned exact routine works from the design's crossproduct, where
   rounding makes an exactly redundant column look as if up to 4e-7 of its length were left, so it drops
   at 1e-5; LSQR, above the cutoff, stops at a tolerance with the same effect. A column between the two
   thresholds is kept by lm and dropped here. dec-B421 accepted that ("Band is fine.") when shown
   differences of 0.01 to 4.7 percent, and 14 percent for one dense column. The second critique found a
   design where it is 250 percent, reproduced here: a timestamp column holding 0 where the value is
   missing, beside the indicator of those rows, which is how a missing value is usually carried in a sparse
   matrix. With the times within an hour of each other the column is, to 6e-7 of its length, a combination
   of the intercept and the indicator (6e-6 at ten hours, still dropped on a sparse x). Where the response
   depends on the time, lm gives 0.536 and the routine 1.796 (sd(y) is 1.899): 235 and 258 percent high at
   n 300 and 3000. The same columns in an ordinary dense data frame that also holds a five-level factor
   take the same routine (dec-B422 and the shared builder) and come out 269 percent high within an hour's
   spread: 1.714 against lm's 0.465, where 0.9-34 gives 0.531 and today's build 0.533 on like data, so as
   planned this is a regression against both. A near-copy of a column that carries the signal is 136 to
   148 percent high between 2e-7 and 3e-6. LSQR at its tolerance of 1e-6 does the same above the cutoff.
   Options.
   (a) As ruled. Nothing to build. Cost: such designs start sigma 2.4 to 3.7 times too high, near sd(y),
   with no message, default dense frames with a factor among them. The sensitivity note found a 2 to 4
   times overestimate moves the posterior mean of sigma by 23 to 57 points of the truth at n 200, by 2 to
   6 at n 1000, and not measurably at n 5000.
   (b) Fall back to a QR whenever the Cholesky drops a column. It cannot be narrowed to the doubtful
   columns: an exactly redundant column and the timestamp look alike in the crossproduct (1.5e-13 against
   3.6e-13), so every design with a redundant column falls back, which is every indicators-route design.
   Cost: about 15 lines; when triggered, the dense fit (64 s at n 2e4, p 2,000, and the design held dense)
   or the sparse QR this plan removes (65 s there, 157 s on one-hot 20 x 50 at n 1e4). Where triggered it
   is lm.
   (c) Recheck what the Cholesky dropped. Each dropped column is measured again in the original space,
   where rounding does not square: its residual against the kept columns, compared with lm's own 1e-7; a
   survivor joins the fit. With it, a centered column gets lm's verdict against the intercept, and LSQR's
   tolerance goes from 1e-6 to 1e-10 with the standard stop on the condition estimate. Cost: about 45
   lines. At m 2,000: no time when nothing is dropped, 0.3 s for 20 dropped columns, 3.8 s against 0.8 s
   for 1,000, 6.2 s against 3.9 s on a wide rank-deficient design; up to 220 MB against 140. On the
   400-design sweep nothing changes (398 ranks equal to lm's, the same largest differences, no column
   re-admitted). The timestamp pair gives lm's 0.536, sparse and in the dense frame; near-copies have lm's
   rank from 1e-9 to 1e-3; a row at 1e12 to 1e15 times the weight of the rest has lm's rank. LSQR at
   1e-10 agrees with lm to 3e-9 on the timestamp design and to under 0.05 percent on the near-copies, in
   9 to 27 iterations where 1e-6 took 4 to 20, and a slow design goes from 704 iterations to the cap with
   1e-8 left. What remains: LSQR then keeps a column down to about 1e-8 of its length where lm stops at
   1e-7, and in that decade is below lm (72 percent on a timestamp spread over 36 seconds); the two wide
   designs of the sweep (0.26 and 0.64 percent); and nothing here was measured above m 2,240.
   Recommended: (c). A user expects the dense equivalent's sigma (dec-B403), and (c) makes the rule lm's
   own tolerance, measured where it can be measured, in one routine. The steps marked OC1 are written for
   it and wait on this ruling.
2. How far from the reference may LSQR be before a design comes back to the maintainer?
   Background. dec-B426 lets LSQR ship on the condition that a badly conditioned design "that reaches the
   cap with a large error comes back to the maintainer before landing". "Large" has no number yet, and
   the build's verdict turns on it. The sensitivity note measured how much a fit depends on the starting
   sigma: at n 1000, sd(y) in place of the linear estimate (15 to 75 percent above it in the note's
   settings) moved RMSE and coverage by less than a reseed does, and the posterior mean of sigma by 1.6
   points of the truth where a reseed moves it 0.8; at n 5000 nothing exceeded a reseed. LSQR runs only
   where n is above about 2,200. One number is known in advance: the crossed-factor design at m 10,000 and
   n 11,200 is 11.8 percent high from its degrees-of-freedom count alone (8.0 and 4.9 percent at the two
   larger sizes).
   Options. (a) 10 percent: under the smallest error the note measured (15 percent), which at n 1000 moved
   sigma by at most 1.6 points; that one crossed design comes back whatever the solver does. (b) 1
   percent, the note's "to a percent or better": every design with dependencies the structure does not
   show comes back (all six crossed designs), so the review is of the count more than of the solver.
   (c) A factor of 2, the largest error the note found close to a reseed at n 1000 (sigma 1.8 points):
   nothing measured so far comes back but a tolerance stop on the band's designs, which Open calls 1
   decides.
   Recommended: (a).
3. Does a design that takes sd(y) under the 10 percent rule get a warning?
   Background. dec-B426 says LSQR runs "with no message" and that above the cutoff a design with under 10
   percent residual degrees of freedom "takes sd(y)"; it does not say whether that is announced. Every
   other fallback to sd(y) raises dbartsSigmaFallbackWarning (a dense design with p >= n; a sparse one
   with no residual rank at or below the cutoff will too), and xbart reads that warning to choose its
   per-fold route. A wide sparse design lands here. So does a dense frame with a high-cardinality factor,
   whatever this call decides: at n 1000 with a factor of 1,000 levels each present once, today's default
   gives 1.145 with no warning (the factor as codes), and after dec-B422 it follows
   `factors = "indicators"`, sd(y) 1.533 with the warning, because no residual degrees of freedom are
   left; at n 5000 with 4,600 levels 399 structural degrees of freedom are left, under the 500 line, so
   that one is the 10 percent case.
   Options. (a) Warn with the existing class, the message naming the rule and `sigest`: one rule, sd(y)
   is always announced. Cost: a warning on every such fit until sigest is given, as wide dense fits have,
   now also on default fits with an identifier-like factor that today fit silently. (b) Silent: about 8
   lines (the routine reports its route so xbart can still choose). Cost: a wide sparse fit warns at m
   2,000 and not at 2,001, and `sigest` reports sd(y) unmarked. (c) Warn only with no residual degrees of
   freedom, silent between 0 and 10 percent: about 10 lines; matches the dense path where it warns and
   hides the new rule where it is new.
   Recommended: (a).
4. Should the 10 percent rule stop at the cutoff?
   Background. By dec-B426 the rule applies above m 2,000 only. A design with 3 percent residual degrees
   of freedom gets the exact estimate at m 2,000 (relative sd about 9 percent on 60 degrees of freedom)
   and sd(y) at m 2,001 (15 to 164 percent above the exact estimate in the note's row).
   Options. (a) As ruled; the step at the cutoff stays; nothing to build. (b) The rule for every sparse
   design: about 4 lines, and a sparse design then differs from its dense equivalent below the cutoff,
   against dec-B403. (c) The rule for dense designs too: about 10 lines, posterior-changing for every
   design with few residual degrees of freedom, a wider re-record, and its own evidence (the note compared
   sd(y) with the linear estimate at 19 degrees of freedom only at n 200).
   Recommended: (a) for this slice, and (c) as its own TODO item if wanted; it does not block the build.

## Estimate

Implementer about two and a half days (routines, recheck and design builder one and a half; tests and docs
one). Condition (a) about five hours of machine time at two jobs: the exact reference at m 10,000 is two
minutes and 1.8 GB a design; the tight LSQR reference took 1,612 iterations at small scale, 20 to 30
minutes a design at m 50,000 if the count holds. Gates about two hours, the re-record minutes. Review with
mutants half a day, and two fix rounds.
