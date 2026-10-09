# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403); revised the same day after its blind critique.

agent: sonnet implementer, one (R only, no engine code); blind critique of this plan first (opus); one opus
reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a fit with a caller's sparse x (a dgCMatrix, or a sparseFactor, sparseVector or
dgCMatrix column in a data frame), a continuous response and no sigest: the default residual prior is
calibrated against the linear-model estimate where it was sd(y). Rounding only, possibly none, for an
indicator expansion R built sparse (dec-B370), whose estimate moves from the sparse QR to the routine below.
NEUTRAL for every dense design, every binary family and every fit given sigest.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~450 lines (R/utility.R ~110 net of the QR's removal, R/spec.R ~25, R/xbart.R ~20, tinytest ~200,
man ~30, NEWS ~2, docs ~60, MANIFEST one row).

## Goal

A caller's sparse x gets the starting sigma its dense equivalent gets, with no warning, at a cost that grows
with the smaller of its row and column counts rather than with refactoring a QR per dependent column. One
routine serves every sparse design, an indicator expansion R built sparse included; the sparse QR and the
dbartsSparseSigmaFallbackWarning class go. xbart's per-fold estimate takes the same routine on each fold's
training rows. An infinite entry in a sparse x is refused as in a dense one, and a fit under a fixed
residual prior makes no estimate at all. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

- The ruling (dec-B403, TODO sparse-starting-sigma): same linear-model starting sigma as the dense
  equivalent, no warning; a caller can pass sigest and the help says so; how the fit is computed is an
  implementation call, measured here.
- Today: [`estimateSigmaFromLinearModel`](../../R/utility.R) returns sd(y - offset) with the warning for a
  caller's sparse source; an indicator expansion (the `sparse.from.indicators` attribute set by
  [`makeIndicatorModelMatrix`](../../R/utility.R)) goes to [`sparseResidualStandardError`](../../R/utility.R),
  Matrix's sparse QR on unit-norm columns, batched at n columns, one dependent column dropped per refactor;
  a dense design goes to [`residualStandardError`](../../R/utility.R) (`lm.fit` on `cbind(1, x)`).
  [`floorSigmaEstimate`](../../R/utility.R) turns a non-finite result into sd(y) with
  dbartsSigmaFallbackWarning. Callers: [`estimateStartingSigma`](../../R/spec.R) (every sampler and fit) and
  [`xbart`](../../R/xbart.R) (once on all rows, then per fold in `foldData`, whose "linear" route densifies
  the design with [`sigmaDesignMatrix`](../../R/utility.R) once per chunk).
- Measurements (2026-10-09, arm64 macOS M1 Max, R 4.6.1 with R's reference BLAS and LAPACK, two jobs at a
  time; scripts in scratch/ssp/ of the main checkout, untracked). Designs: random sparse numeric columns at
  0.1, 1 and 5 percent; full one-hot factor sets (every level a column, so dependent with the intercept by
  one column a factor), written as F x L; response a sparse linear signal plus noise of sd 0.5. Seconds:

  | n | p | design | dense lm.fit | sparse QR (today) | candidate (v3) |
  |---|---|---|---|---|---|
  | 1e3 | 100 | 0.1-5% | 0.01 | 0.01 | 0.03 |
  | 1e3 | 1000 | 1% | 0.49 | 0.28 | 0.14 |
  | 1e3 | 5000 | 0.1% / 1% | 3.5 (1%) | 0.73 / 0.27 | 1.45 / 0.15 |
  | 1e3 | 5000 | one-hot 50 x 100 | 3.6 | 3.6 | 0.15 |
  | 1e4 | 100 | 1% / one-hot 10 x 10 | 0.09 / 0.09 | 0.03 / 0.86 | 0.04 / 0.04 |
  | 1e4 | 1000 | 1% / one-hot 20 x 50 | 7.8 / 7.8 | 7.4 / 157 | 0.18 / 0.20 |
  | 2e4 | 2000 | 1% | 64 | 65 (dec-B403) | 1.2 |
  | 1e4 | 5000 | 0.1-5% / one-hot 50 x 100 | ~400 (2np^2 at the measured rate) | not run | 14.5-17.1 / 15.4 |
  | 1e5 | 100 | 1% / one-hot 10 x 10 | 1.0 / 1.0 | 1.0 / 11.7 | 0.04 / 0.15 |
  | 1e5 | 1000 | 0.1-5% / one-hot 20 x 50 | 83 (1%) / 83 | not run | 0.16-1.6 / 0.45 |
  | 1e5 | 5000 | 0.1% / 1% / 5% / one-hot 50 x 100 | ~4000 and a 4 GB matrix (extrapolated) | over 20 min (dec-B403) | 15 / 17 / 38 / 17 |

  The candidate column is the first prototype (v3); the revised routine (v4, Algorithm) times within 10
  percent of it at four of these points (17.5 s at n 1e5, p 5000, 1 percent; 0.22 s at n 1e4, one-hot
  20 x 50).

  The candidate's cost is the crossproduct (0.2, 2.6 and 19 s at n 1e5, p 5000 and 0.1, 1 and 5 percent)
  plus a pivoted Cholesky of an m x m matrix, m the smaller of n and p + 1. Measured, the whole routine:
  1.0 s at m 2000, 14 s at 5000, 133 s at 10000; peak resident memory 351 MB, 654 MB and 1.84 GB, of which
  about 250 MB is R with Matrix and dbarts loaded, so about 16 m^2 bytes (the dense crossproduct and the
  factor `chol` returns). The Cholesky alone: 13 s at 5000, 37 s at 7000, 118 s at 10000 and 592 s at
  14000 (3.2 GB), steeper than cubic past 10000 on this reference BLAS (a local exponent of 4.8).
  Extrapolated: m 2e4, 0.5 to 1 hour (exponent 3 to 4.8) and 6.4 GB; m 5e4, at least 7 hours and about
  40 GB, more than this 32 GB machine holds (an allocation error or swapping); m 1e5, about 160 GB, an
  allocation error. Neither LAPACK call can be interrupted. An optimized BLAS is many times faster; the
  memory does not change. The dense path at the same design costs at least 2 m^3 flops in `lm.fit`'s
  LINPACK QR, six times the Cholesky's m^3 / 3 (27 minutes at n = p = 1e4 at the measured rate), and holds
  n p doubles against 2 m^2. Where both take a second or two the routine can be the slower: 0.03 s against
  0.01 s at n 1e3, p 100; 2.0 s against 0.6 s on a tall mostly dense frame (n 2e5, 50 dense columns and
  one sparse, the dense block costing about 12 bytes an entry as CSC).
- xbart: the per-fold estimate runs n.reps x the fold count times, 200 at the defaults (today a sparse
  design's per-fold value is a marginal sd). At n 1e4, p 1000, 1 percent and one-hot 20 x 50: 28 and 32 s
  for the 200 (0.14 to 0.16 s each), where a dense design of the same size pays 5.6 to 5.9 s a fold, about
  19 minutes. Each worker holds its own crossproduct, so memory multiplies by the worker count.
- Agreement with the dense fit (`lm.fit` on `cbind(1, as.matrix(x))`, through
  [`residualStandardError`](../../R/utility.R)), the revised routine (v4: Algorithm below): 400 random
  designs (n 30 to 2000, p 5 to 1530, numeric, one-hot and mixed, a third with three dependent columns
  added, a fifth with a column scaled by 1e5, a third weighted with 5 percent zero weights, a third with an
  offset; 143 with p >= n). The rank equals lm.fit's on 398. Of those, 106 have no residual degrees of
  freedom and both give no estimate; of the 292 with an estimate, sigma agrees within 1e-11 relative at
  more than 30 residual degrees of freedom (245 designs) and within 3e-8 at 30 or fewer (47; 2.6e-8 at 10
  or fewer). On the other 2, wide designs, lm.fit's rank is above the SVD's numerical rank (975 against
  970, 935 against 933; it keeps columns with relative pivots near 1e-16), the routine's equals the SVD's,
  and sigma is 0.26 and 0.64 percent apart; the shipped sparse QR sides with the SVD there too. The
  critique reran the sweep independently on the v3 routine and reproduced these counts.
- Hand cases (v4, scratch/ssp/verify4.R and band4.R): duplicated, constant and all-zero columns, exact
  linear combinations, zero weights, offsets, NA entries in a dgCMatrix, sparseFactor columns with each
  reference level and with missing values, a mixed frame (dense factor, numeric and sparseVector columns),
  a column scaled by 1e160 or 1e-160, a dense timestamp column (mean 1.7e9, spread an hour) carrying
  signal, weighted and not, a fully stored sparse column of the same kind, a year column, dense columns of
  mean 2e3 and 1e7 beside an exact duplicate: all within 1.3e-11. A wide full-rank design with one row at
  relative weight 1e-14 gives no estimate, as lm.fit does (v3 gave 3.4e-8 at rank 49).
- The tolerance band. A column is dropped when its residual after the kept columns is below 1e-5 of its
  norm (pivot 1e-10 on the unit-diagonal crossproduct); lm.fit's is 1e-7. The pivot of an exactly
  dependent column measured 0 to 1.5e-13 (one-hot 50 x 100 at n 30000; 5e-14 at 20 x 100; 1.7e-14 for
  100 linear combinations of 1000 numeric columns), so 1e-10 holds a margin of about 700 and a tolerance
  matching lm's (1e-14) none. Two near-copies whose difference is between 1e-7 and 1e-5 of their norm
  keep both under lm and one here, and sigma moves by whatever the second copy fits: on
  test-indicator-storage.R's design (n 120) and two redraws, perturbations from 1e-7 to 1e-5 give one rank
  fewer and sigma 0.01 to 4.7 percent apart, each time equal within 1e-15 to the dense fit without one of
  the two copies, the routine choosing which. Since dense-backed and fully stored columns are centered
  (Algorithm), a column of large mean and small spread is no longer in that band against the intercept:
  a timestamp spanning an hour, 42 percent apart under the uncentered v3, agrees within 6e-13. The band
  turns over for such a column instead: where the spread is below 1e-7 of the mean (a mean of 1e9 with an
  sd of 6), lm drops the column as constant and the routine keeps it (14 percent apart on that design,
  where the column carried signal). See Open calls.
- What moves (an instrumented build logging each call of
  [`estimateSigmaFromLinearModel`](../../R/utility.R), 2026-10-09; xbart's per-fold fits were not logged,
  but a sparse xbart also makes the logged all-rows call):
  - equivalence.R in quick mode: the seven scenarios built on a caller's sparse source, each without sigest,
    ["sparse <- list("](../../benchmarks/R/equivalence.R),
    ["mixedmatrix <- list("](../../benchmarks/R/equivalence.R),
    ["sparsefactor <- list("](../../benchmarks/R/equivalence.R),
    ["testswap <- list("](../../benchmarks/R/equivalence.R),
    ["leaffactormixed <- list("](../../benchmarks/R/equivalence.R),
    ["factorpartial <- list("](../../benchmarks/R/equivalence.R) and
    ["xbartmixed <- list("](../../benchmarks/R/equivalence.R); their starting sigma goes from sd(y) to the
    linear estimate: 2.959 to 1.370, 1.670 to 0.327, 3.609 to 2.387, 3.664 to 2.508, 3.941 to 1.867,
    3.483 to 2.235 and (all rows) 4.124 to 2.626, each within 9e-16 of lm.fit's. wideFactorIndicators
    takes the QR route today and moves by rounding: v4 gives 2.6995748335729592 against the QR's
    2.6995748335729584 (lm.fit 2.6995748335729566), so its draws change and its posterior does not. The
    other 47 make dense estimates or none, and none uses a fixed residual prior.
  - bcf-equivalence.R: 12 estimates, all dense. multinomial-equivalence.R: none (no residual scale).
  - The four test-reproducibility files (run with the build guard bypassed): 5 estimates, all dense, none
    under a fixed residual prior.
  - exact-gates.yaml's gate list in quick mode, bcf-latent-exact.R aside: 35 estimates, all dense, every
    gate passing; bcf-latent-exact.R builds no Matrix, sparseFactor or indicator design (grep). No exact
    gate moves.
  - tinytest, the whole suite (19630 results, 0 failures): the caller-sparse branch is reached in
    test-data-mixed.R, test-indicator-storage.R, test-predict-na-action.R and test-row-names.R; the QR in
    test-dart-mixed-columns.R, test-indicator-storage.R and test-sampler-splitProbabilities.R. Every other
    sparse test passes sigest.
- A fixed residual prior (`gaussian(sigma = fixed(v))`) still runs the estimate today
  ([`dbartsSpec`](../../R/spec.R) calls [`estimateStartingSigma`](../../R/spec.R) unless the family has a
  fixed unit scale), and the creation discards it for sqrt(v): draws are bit for bit the same with the
  slot left NA. The slot is read later, though: a `setModel` to a chisq prior on a sampler whose slot is NA
  draws NA sigmas without an error (measured), `rbart_vi` starts its sigma from the slot, and `bart`
  records it as `sigest`. Hence the rule in Change 4.

## Algorithm

`sparseResidualStandardError(y, x, weights, offset)` keeps its name and arguments; `x` is a predictor source
(a dgCMatrix or a mixed container).

1. Design, in this order. A bare dgCMatrix is first wrapped with
   [`wrapSparseTestMatrix`](../../R/mixedMatrix.R). Then a sparse categorical column (non-NA
   `sparseReference`) has its reference code subtracted from its stored entries, so the column equals
   [`as.matrix.dbartsMixedMatrix`](../../R/mixedMatrix.R)'s (implicit rows at the reference code) up to a
   constant, which the intercept absorbs. Then [`sparseDesignMatrix`](../../R/utility.R) as today (dense
   block to CSC first, NAs mean-imputed with the implicit zeros counted). The shift comes before the
   imputation, since the imputation's mean counts the implicit rows as 0: imputing first is 0.4 to 2.2
   percent off on a sparseFactor with a middle reference and missing values. Without the shift at all, a
   sparseFactor whose reference is not its first level fits another model (0.2 percent off). The builder
   records which columns are dense-backed.
2. Inf. An infinite stored entry is an error, which [`estimateStartingSigma`](../../R/spec.R) reports as
   it does for a dense design, naming the column: [`nonFinitePredictorNames`](../../R/spec.R) learns sparse
   sources (the dense list and the CSC block's stored entries, mapped to names through the container's map
   and column names, or positions for a bare dgCMatrix). Today a sparse Inf is accepted with the sd(y)
   fallback; the uncentered routine would have dropped the column silently (13 percent off).
3. Rows. Drop rows with a missing response, weight or offset and rows of weight 0, as `lm.wfit` does;
   z = y - offset.
4. Centering. Each dense-backed column, and each column with a stored entry in every kept row, is
   centered at its weighted mean over the kept rows (they are dense already; centering the values leaves
   the span with the intercept unchanged). Sparse columns are not centered: a centered crossproduct
   (X'WX - s m m') cancels for a large-mean column (pivot noise 2.5e-11 at a mean of 2000 and sd 6, a margin
   of 4), which is why the intercept stays a column.
5. Scaling. B = diag(sqrt(w)) [1 X]; each column divided by its largest absolute entry, then by its norm
   (so a column scaled by 1e160 neither overflows nor drops); all-zero columns dropped.
6. Narrow. If ncol(B) < n, `chol(as.matrix(Matrix::crossprod(B)), pivot = TRUE, tol = 1e-10)` (LAPACK
   dpstrf) with `suppressWarnings` around that one call only (its one warning is the rank deficiency; no
   matching on message text), r its rank, K its first r pivots, R its leading r x r block;
   b = R^-1 R^-T B_K' (z sqrt(w)), e = z sqrt(w) - B_K b, computed directly (no normal-equation RSS). The
   error in b is orthogonal to the true residual, so its effect on the RSS is second order: one step of
   iterative refinement changed no result by more than 1e-12 on 30 designs with a kept column near the
   tolerance and two columns of mean 1e4. It is not taken.
7. Wide. Otherwise K = `as.matrix(Matrix::tcrossprod(B))` (n x n), equilibrated: D = sqrt(diag(K)) (the row
   norms of B, positive since the intercept column has an entry in every row), the pivoted Cholesky of
   D^-1 K D^-1 at the same tolerance. Its pivots are then relative to each row's norm, as B'B's are to
   each column's; unequilibrated, a row at relative weight 1e-14 is taken as dependent and a full-rank
   design gets sigma 3.4e-8 (B2). r >= n gives no estimate; else Q an orthonormal basis of D L (L the first
   r columns of the pivoted factor, unpivoted; `qr.Q(qr(.))`, n x r) and e = z sqrt(w) - Q Q' z sqrt(w).
   On n 1000, p 5000 at 0.1 percent (rank 993) this takes 1.5 s where the narrow side takes about 14.
8. sigma = sqrt(sum(e^2) / (n - r)), NA when n - r <= 0, for [`floorSigmaEstimate`](../../R/utility.R) to
   take to sd(y) with dbartsSigmaFallbackWarning, exactly as the dense path does with more columns than rows.

No size cap, no message (Open calls). An allocation failure is an error, which
[`estimateStartingSigma`](../../R/spec.R) already turns into "unable to obtain a starting estimate of sigma;
provide one instead".

## Change

1. R/utility.R: [`sparseResidualStandardError`](../../R/utility.R) and its comment rewritten as above (the
   batching, the refactor loop and the `grepl` muffler of Matrix's structural-rank warning go);
   [`sparseDesignMatrix`](../../R/utility.R) gains the dgCMatrix wrap, the reference shift ahead of its
   imputation and the dense-backed record; [`estimateSigmaFromLinearModel`](../../R/utility.R) loses the
   fallback branch and its warning and sends every sparse source to the sparse routine;
   [`makeIndicatorModelMatrix`](../../R/utility.R) no longer sets `sparse.from.indicators`, its one reader
   gone.
2. R/spec.R: [`nonFinitePredictorNames`](../../R/spec.R) for sparse sources (Algorithm 2).
3. R/xbart.R: [`xbart`](../../R/xbart.R)'s comment on the per-fold route drops "a sparse design"; in the
   chunk runner a sparse source's `sigmaDesign` is built once per chunk by the design builder, and
   `foldData` fits the sparse routine on its training rows (centering at each fold's own weighted means)
   where a dense source keeps `residualStandardError`. An indicator expansion stops being densified there.
4. A fixed residual prior runs no linear fit: in [`dbartsSpec`](../../R/spec.R), where `residPrior` is a
   dbartsFixedPrior, the data's sigma slot takes the sigma it fixes, sqrt of the fixed variance, in place
   of the estimate, for dense and sparse designs alike; xbart's all-rows estimate the same. Creation draws
   are unchanged (they read the fixed value). What changes: `bart`'s `sigest` and the slot under a fixed
   prior report the fixed sigma, not an estimate the fit never used, and a later `setModel` to a chisq
   prior, or `rbart_vi`'s start, reads that value. The slot is never left NA (Context). This is the
   critique's S3 in the form that keeps every reader of the slot finite.
5. The class dbartsSparseSigmaFallbackWarning is retired: it was new in 1.0 and nothing raises it. No
   consumer branch references it (the critique's grep of stan4bart, bartCause, treatSens and bairrtt).

## Constraints

- No engine, bridge, C API or state change. No new dependency: Matrix (Suggests) as today, base `chol`,
  `backsolve`, `qr` and `tapply`, all within R 4.2 and Matrix 1.4-1.
- Dense designs under a chisq prior untouched: [`residualStandardError`](../../R/utility.R) and every
  dense caller bit for bit.
- Out of scope: how a categorical column enters the linear fit (Open calls); a per-fold crossproduct
  downdate in xbart (the dense path refits per fold too).

## Tests

- [test-data-mixed.R](../../inst/tinytest/test-data-mixed.R), the block counting
  ["dbartsSparseSigmaFallbackWarning"](../../inst/tinytest/test-data-mixed.R): zero
  dbartsSigmaFallbackWarning for the sparse frame and its dense equivalent, and the two fits' `sigest`
  equal within 1e-10.
- [test-indicator-storage.R](../../inst/tinytest/test-indicator-storage.R):
  - ["a sparse column the caller supplied still falls back"](../../inst/tinytest/test-indicator-storage.R)
    becomes no warning and `data@sigma` equal within 1e-10 to `residualStandardError` on
    `sigmaDesignMatrix` of the same source.
  - The weights, offset and missing-value loop and the "other units" and more-columns-than-rows cases keep
    their 1e-8 against the dense fit.
  - ["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R): every perturbation in
    the series (0, 1e-9, 5e-8, 2e-7, 1e-5, 3e-5 and 1e-3 today, with 1e-6 and 1e-4 added) is pinned within
    1e-10 to the nearest of three dense fits: with all columns, without `a` and without `c`. The routine
    picks which near-copy it drops, so no pin names one; 5e-8, at 5.9e-9 from the full dense fit, is no
    longer pinned at 1e-8 against it. Below 1e-7 and from 3e-5 up the nearest is the full fit (the
    existing 3e-5 case stays kept, 2.4e-14 apart); 1e-6 additionally asserts one rank fewer than lm.fit's.
- New, in test-indicator-storage.R or a new test-sparse-starting-sigma.R (under 2 s), each against the
  dense fit within 1e-10 unless said:
  - a bare dgCMatrix through `dbarts`, `bart` and `xbart`, with no warning (counted, Gate hygiene);
  - one-hot 3 x 40 plus three duplicated columns, weights with zeros, an offset;
  - a dgCMatrix with NA entries (the mean imputation);
  - sparseFactor columns with reference the first level and a middle level, and a middle reference with
    missing values, in a frame with a numeric column (the shift and its order);
  - a frame with a dense timestamp-like column (mean 1.7e9, spread an hour) carrying signal, weighted and
    unweighted, and a dgCMatrix with such a column stored in every row (the centering);
  - a column scaled by 1e160 (the max-abs scaling);
  - p > n with rank below n - 1 (n 150, p 240, rank 41: the wide side); p > n of full rank: NA from the
    routine and, through `dbarts`, sd(y) with one dbartsSigmaFallbackWarning, as the dense design gives;
    the same full-rank wide design with one row at weight 1e-14 of the others: NA (the equilibration);
  - an Inf entry in a dgCMatrix and in a sparseVector column of a frame: the dense path's error, naming
    the column;
  - xbart on a sparse frame: each fold's sigma equals the dense design's per-fold value for the same rows
    (the fold oracle's route, test-xbart-fold-oracle.R's helper or a direct `foldData` probe);
  - a fixed residual prior on a dense and a sparse design: no linear fit (a traced or mocked
    `estimateStartingSigma` is not called), the slot and `sigest` equal sqrt of the fixed variance, draws
    identical to the build before; then `setModel` to a chisq prior: finite sigma draws.
- test-predict-na-action.R and test-row-names.R reach the routine and pin nothing about sigma; they must
  pass unchanged. The 21 test files that use a fixed residual prior are checked for a pinned `sigest` or
  sigma slot under it; such a pin moves to the fixed sigma, any other change is a stop.

Reviewer's mutants, each of which must fail a test: tolerance 1e-10 to 1e-16 (exact dependencies kept) and
to 1e-6 (the 1e-4 column dropped); the reference shift removed; imputation before the shift; the
intercept column dropped; centering removed (the timestamp cases, 42 percent off); centering at the
unweighted mean on the weighted case; the max-abs step removed (the 1e160 case); the equilibration removed
(the weight-1e-14 case); weights left out of B or of the residual; the wide side's basis taken from the
unpivoted factor or without D; zero-weight rows kept; `n - r` replaced by `n - p`; the Inf check removed;
the fixed-prior skip leaving the slot NA (the `setModel` case).

## Baselines

- Moves: equivalence.R's seven caller-sparse scenarios and wideFactorIndicators (Context; the last by
  rounding, unless the implemented routine happens to round to the QR's value). Current file
  equivalence-734441f1 ([MANIFEST](../../benchmarks/baselines/MANIFEST)).
- Does not move: the other 47 equivalence scenarios, bcf-equivalence-1b7d730c (15),
  multinomial-equivalence-80b1c8d4 (11), the four test-reproducibility files, every exact gate (Context;
  none fits under a fixed residual prior or records its slot). Any other mover is a defect: stop.
- Re-record on the reference build (`--preclean --configure-args=--enable-reference-build`):
  `EQUIVALENCE_SCENARIOS=sparse,mixedmatrix,sparsefactor,testswap,leaffactormixed,factorpartial,xbartmixed,wideFactorIndicators`
  (wideFactorIndicators dropped if it came out bitwise), `EQUIVALENCE_CORES=2`, merged into a copy of
  734441f1 in its scenario order, named after the slice's code commit; 734441f1 demoted to historical.
  Partition against 734441f1 in z mode: 47 (or 48) of 55 identical, the movers exactly those named,
  wideFactorIndicators with no |z| above 4. With the prior's anchor falling by factors of 0.2 to 0.7 on
  150 to 500 rows, sigma and fit summaries past |z| 4 are expected in the seven and are not a failure; the
  partition's verdict is which scenarios moved. The merged file reproduces 55 of 55 under
  `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17), as deb3fe50's row: the change is the starting value a prior is calibrated
  against, not the sampler; the identity is the routine's agreement with `lm.fit` (the tinytest pins and
  the 400-design sweep above, rerun by the implementer from scratch/ssp/sweep4.R's recipe against the
  implemented function).

## Gates

On the slice tip against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing), with the expected
verdict:
- tests/cpp: unchanged, green (no C++ touched).
- Full tinytest suite: green with the edits above.
- The four seeded-drift snapshot files on the reference build: pass unchanged (27 results).
- The equivalence trio: as Baselines.
- exact-gates.yaml's list in quick mode: all pass, output as before (none reaches the routine).
- `R CMD check --as-cran` from a clean tarball: no new NOTE; lintr, air, rc-codoc, win-drift,
  doc-freshness. Sanitizers: not owed (no compiled code changes); the routine's numerics are base R's and
  Matrix's.
- Speed, same machine, within 1.5x of Context: the routine at n 1e5, p 5000 (1 percent; 17.5 s) and n 1e4,
  p 1000 (one-hot 20 x 50; 0.22 s); xbart's 200 per-fold estimates at n 1e4, p 1000, 1 percent (28 s,
  timed around the per-fold calls or as the difference to a run given sigest). No bench-sampler.R compare
  (no hot path).
- The design note: [R surface](../design/sparse-columns.md#r-surface) in sparse-columns.md says "the sigma
  estimate falls back to sd(y)"; a dated paragraph there states the rule and Algorithm in brief, citing
  dec-B403 (the landing notes stay as written).

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd and man/xbart.Rd: the clause "a design with
  sparse-backed predictor columns skips the linear model altogether and falls back the same way (class
  dbartsSparseSigmaFallbackWarning, a dbartsSigmaFallbackWarning)" becomes: a sparse x takes the same
  estimate as its dense equivalent, fitted from its crossproduct without making it dense; its time grows
  at least with the cube, and its memory with twice the square, of the smaller of its row and column counts
  (seconds up to a few thousand; at ten thousand about two minutes and 2 GB, at twenty thousand most of an
  hour and over 6 GB with R's reference BLAS), and it cannot be interrupted, so for a large sparse design
  supply `sigest`. Under a fixed residual prior no estimate is made. man/xbart.Rd adds that the estimate is
  repeated for every fold of every repetition, each worker holding its own crossproduct.
- man/bart.Rd's warning-class paragraph drops dbartsSparseSigmaFallbackWarning; its `sigest` value under
  a fixed prior is the fixed sigma.
- man/sparseFactor.Rd: the paragraph on the default starting sigma is replaced by one sentence: a sparse
  column leaves the default starting sigma as the same column stored dense gives it.
- inst/NEWS.Rd: dbartsSparseSigmaFallbackWarning leaves the warning-class list (never released, so no
  entry of its own).
- docs/design/error-style.md: the class leaves the list and the table.
- docs/design/memory-footprint.md: the starting-sigma row, which prices only the dense `lm` path
  (2 n (p + 1)), gains the sparse route's transient, about 16 m^2 bytes, m the smaller of n and p + 1, per
  worker under xbart.
- At landing: TODO's item goes; ledger entry for the calls below; this plan's Status and Landing note;
  MANIFEST row.

## Steps

1. Algorithm and Change 1 to 5; the tinytest edits and new tests; the full suite green against
   `R CMD INSTALL -l <lib> .`.
2. The 400-design sweep against the implemented function; the speed points; help, NEWS and docs.
3. After review: the re-record, MANIFEST row and partition, in their own commit.

## Stop conditions

Stop and report when: the diff passes ~700 lines; an equivalence scenario other than the eight, a
snapshot or an exact gate moves, or wideFactorIndicators shows a |z| above 4; the sweep shows a rank
different from the SVD's or a sigma off by more than 1e-10 at more than 30 residual degrees of freedom
outside the band; a speed point is past 1.5x; a test under a fixed prior changes for any reason but the
slot's value; the change needs engine, bridge or C API code.

## Interactions

- Whichever lands second of this and any other slice re-recording equivalence.R re-records against the
  other's file, and partitions against it.
- The categorical question below: if it is ruled before this lands, the sparseFactor scenarios move once.

## Calls made in planning

- One routine for every sparse design: the sparse QR goes, so an indicator expansion takes this path too
  (157 s against 0.2 s on one-hot 20 x 50 at n 1e4).
- The smaller Gram side with the intercept as a column; dense-backed and fully stored columns centered at
  the weighted mean over the kept rows; max-abs then unit-norm column scaling; the wide side equilibrated;
  tolerance 1e-10 on unit-diagonal pivots; residual computed directly; no refinement.
- The design's order: wrap, shift, impute.
- An infinite entry in a sparse x is refused as the dense path refuses it, the column named (coordinator,
  for VD's rule that the estimate follow the dense equivalent).
- xbart's per-fold estimate uses the routine per fold, as the dense path refits per fold.
- Under a fixed residual prior no estimate is made and the slot holds the fixed sigma (Change 4).
- dbartsSparseSigmaFallbackWarning retired.

## Open calls

- The band. Two near-copies of a column whose difference is between 1e-7 and 1e-5 of their norm both enter
  `lm`'s fit and one enters this one; sigma differs by what the second copy fits, 0.01 to 4.7 percent on
  the n 120 test design and two redraws, and equals the dense fit without one of the two. A dense column
  whose spread is below 1e-7 of its mean, which lm drops as constant, is kept here (14 percent apart where
  it carries signal). Accept? Recommended: in a crossproduct an exact dependency's pivot reaches 1.5e-13 at
  p 5000, above lm's 1e-14, so no tolerance near lm's separates the two, and matching lm exactly needs a
  QR, the cost this replaces. Centering has already removed the case the critique found largest, a
  timestamp-like column against the intercept (42 percent apart before, 6e-13 now).
- Large designs. With m the smaller of n and p + 1, on R's reference BLAS: 14 s and 0.65 GB at m 5000,
  2.2 minutes and 1.8 GB at 1e4, 10 minutes and 3.2 GB for the Cholesky alone at 1.4e4; extrapolated, 0.5 to
  1 hour and 6.4 GB at 2e4, at least 7 hours and about 40 GB at 5e4, an allocation error at 1e5 (about
  160 GB). Not interruptible. Today such a fit starts at once from sd(y) with a warning; the dense path at
  the same design runs at least six times longer, with no message. xbart repeats it 200 times by default.
  Options: no cap and no message (the ruling's text; the help names `sigest`; recommended for consistency
  with the dense path, though at m above about 2e4 a fit can sit for an hour before its first sweep);
  a message once m passes 1e4 naming `sigest`, then the fit (about 6 lines and a test); or sd(y) with the
  old warning past it (a size at which the prior changes).
- A categorical column enters the linear fit as its codes on the dense path in 1.0 (a factor kept as a
  factor; [`as.matrix.dbartsMixedMatrix`](../../R/mixedMatrix.R)), and this plan gives a sparseFactor the
  same. 0.9-34 expanded factors into indicators: on a 5-level factor and one numeric column (n 300, noise sd
  0.5) 0.9-34's dbarts starts sigma at 0.537 (lm with the factor) and 1.0's at 1.713 (lm on the codes). A
  separate call; if indicators, a sparse categorical column is expanded to indicator columns in the design
  (about 15 lines here) and the dense path changes with it.

## Estimate

Implementer about a day; gates about two hours of machine time, the re-record minutes.
