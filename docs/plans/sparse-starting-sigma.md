# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403).

agent: sonnet implementer, one (R only, no engine code); blind critique of this plan first (opus); one opus
reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a fit with a caller's sparse x (a dgCMatrix, or a sparseFactor, sparseVector or
dgCMatrix column in a data frame), a continuous response and no sigest: the default residual prior is
calibrated against the linear-model estimate where it was sd(y). Rounding only, possibly none, for an
indicator expansion R built sparse (dec-B370), whose estimate moves from the sparse QR to the routine below.
NEUTRAL for every dense design, every binary family and every fit given sigest.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~350 lines (R/utility.R ~90 net of the QR's removal, R/xbart.R ~15, tinytest ~150, man ~25, NEWS ~2,
docs ~50, MANIFEST one row).

## Goal

A caller's sparse x gets the starting sigma its dense equivalent gets, with no warning, at a cost that grows
with the smaller of its row and column counts rather than with refactoring a QR per dependent column. One
routine serves every sparse design, an indicator expansion R built sparse included; the sparse QR and the
dbartsSparseSigmaFallbackWarning class go. xbart's per-fold estimate takes the same routine on each fold's
training rows. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

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

  | n | p | design | dense lm.fit | sparse QR (today) | candidate |
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

  The candidate's cost is the crossproduct (0.2, 2.6 and 19 s at n 1e5, p 5000 and 0.1, 1 and 5 percent)
  plus a pivoted Cholesky of an m x m matrix, m the smaller of n and p + 1: 0.8 s at m 2000, 13 s at 5000,
  37 s at 7000, 118 s at 10000, memory 8 m^2 bytes (0.8 GB at 10000). At the same design the dense path
  costs at least 2 m^3 flops in `lm.fit`'s LINPACK QR, six times the Cholesky's m^3 / 3, and holds n p
  doubles against m^2; the candidate is never the slower or larger of the two.
- Agreement with the dense fit (`lm.fit` on `cbind(1, as.matrix(x))`, through
  [`residualStandardError`](../../R/utility.R)): 400 random designs (n 30 to 2000, p 5 to 1530, numeric,
  one-hot and mixed, a third with three dependent columns added, a fifth with a column scaled by 1e5, a
  third weighted with 5 percent zero weights, a third with an offset; 143 with p >= n): the rank equals
  lm.fit's on 398. Of those, 106 have no residual degrees of freedom and both give no estimate; of the 292
  with an estimate, sigma agrees within 1e-11 relative at more than 30 residual degrees of freedom (245
  designs) and within 3e-8 at 30 or fewer (47). The other 2 are wide designs where lm.fit keeps a column
  whose relative pivot is 1.4e-16, one above the SVD's numerical rank, which the routine's rank equals
  (sigma 0.26 and 0.64 percent apart); the shipped sparse QR sides with the SVD there too. Hand cases:
  duplicated, constant and all-zero columns, exact linear combinations, zero weights, offsets, NA entries in
  a dgCMatrix, sparseFactor columns with each reference level and with missing values, a mixed frame (dense
  factor, numeric and sparseVector columns), columns with means of 2e3 to 1e6 beside exact duplicates: all
  within 1e-13.
- The tolerance band. A column is dropped when its residual after the kept columns is below 1e-5 of its
  norm (pivot 1e-10 on the unit-diagonal crossproduct); lm.fit's is 1e-7. The pivot of an exactly
  dependent column measured 0 to 1.5e-13 (one-hot 50 x 100 at n 30000; 5e-14 at 20 x 100; 1.7e-14 for 100
  linear combinations of 1000 numeric columns), so 1e-10 holds a margin of about 700 and a tolerance
  matching lm's (1e-14) none. A column between 1e-7 and 1e-5 is kept by lm and dropped here: on
  test-indicator-storage.R's collinearity series, perturbations 2e-7, 1e-6 and 3e-6 give one rank fewer and
  sigma 0.9, 0.4 and 0.08 percent apart; 5e-8 and below, and 1e-5 and above, agree within 6e-9. A column
  whose spread is below 1e-5 of its mean (a mean of 1e7 with an sd of 6) is the same case against the
  intercept (7e-5 apart). See Open calls.
- What moves (instrumented build logging every starting-sigma call, 2026-10-09):
  - equivalence.R in quick mode: the seven scenarios built on a caller's sparse source, each without sigest,
    ["sparse <- list("](../../benchmarks/R/equivalence.R),
    ["mixedmatrix <- list("](../../benchmarks/R/equivalence.R),
    ["sparsefactor <- list("](../../benchmarks/R/equivalence.R),
    ["testswap <- list("](../../benchmarks/R/equivalence.R),
    ["leaffactormixed <- list("](../../benchmarks/R/equivalence.R),
    ["factorpartial <- list("](../../benchmarks/R/equivalence.R) and
    ["xbartmixed <- list("](../../benchmarks/R/equivalence.R); their starting sigma goes from sd(y) to the
    linear estimate: 2.959 to 1.370, 1.670 to 0.327, 3.609 to 2.387, 3.664 to 2.508, 3.941 to 1.867,
    3.483 to 2.235 and (all rows) 4.124 to 2.626. wideFactorIndicators takes the QR route today and moves
    by rounding: the prototype gives 2.6995748335729579 against the QR's 2.6995748335729584 (lm.fit
    2.6995748335729566), so its draws change and its posterior does not. The other 47 make dense
    estimates or none.
  - bcf-equivalence.R: 12 estimates, all dense. multinomial-equivalence.R: none (no residual scale).
  - The four test-reproducibility files (run with the build guard bypassed): 5 estimates, all dense.
  - exact-gates.yaml's gate list in quick mode, bcf-latent-exact.R aside: 35 estimates, all dense, every
    gate passing; bcf-latent-exact.R builds no Matrix, sparseFactor or indicator design (grep). No exact
    gate moves.
  - tinytest, the whole suite (19630 results, 0 failures): the caller-sparse branch is reached in
    test-data-mixed.R, test-indicator-storage.R, test-predict-na-action.R and test-row-names.R; the QR in
    test-dart-mixed-columns.R, test-indicator-storage.R and test-sampler-splitProbabilities.R. Every other
    sparse test passes sigest.

## Algorithm

`sparseResidualStandardError(y, x, weights, offset)` keeps its name and arguments; `x` is a predictor source
(a dgCMatrix or a mixed container).

1. Design. [`sparseDesignMatrix`](../../R/utility.R) as today (dense block to CSC, NAs mean-imputed with the
   implicit zeros counted), plus: a bare dgCMatrix is first wrapped with
   [`wrapSparseTestMatrix`](../../R/mixedMatrix.R); a sparse categorical column (non-NA `sparseReference`)
   has its reference code subtracted from its stored entries, so the column equals
   [`as.matrix.dbartsMixedMatrix`](../../R/mixedMatrix.R)'s (implicit rows at the reference code) up to a
   constant, which the intercept absorbs. Without the shift a sparseFactor whose reference is not its first
   level fits another model (0.2 percent off on the test design).
2. Rows. Drop rows with a missing response, weight or offset and rows of weight 0, as `lm.wfit` does;
   z = y - offset. B = diag(sqrt(w)) [1 X], each nonzero column scaled to unit norm, all-zero columns
   dropped. No centering: a centered crossproduct cancels for a column with a large mean (pivot noise 2.5e-11
   at a mean of 2000 and sd 6, a margin of 4 to the tolerance).
3. Rank. If ncol(B) < n, `chol(as.matrix(Matrix::crossprod(B)), pivot = TRUE, tol = 1e-10)` (LAPACK dpstrf;
   its rank-deficiency warning muffled), r its rank, K its first r pivots, R its leading r x r block;
   b = R^-1 R^-T B_K' (z sqrt(w)), e = z sqrt(w) - B_K b, computed directly (no normal-equation RSS, so no
   cancellation). One step of iterative refinement changed no result by more than 1e-12 on 30 designs with
   a kept column near the tolerance and two columns of mean 1e4; it is not taken.
4. Wide. Otherwise the same on `as.matrix(Matrix::tcrossprod(B))` (n x n; its nonzero eigenvalues are
   B'B's, so the tolerance carries over); r >= n gives no estimate; else an orthonormal basis Q of the
   first r columns of the pivoted factor (`qr.Q(qr(.))`, n x r) and e = z sqrt(w) - Q Q' z sqrt(w).
   On n 1000, p 5000 at 0.1 percent (rank 993) this takes 1.5 s where the column side takes about 14.
5. sigma = sqrt(sum(e^2) / (n - r)), NA when n - r <= 0, for [`floorSigmaEstimate`](../../R/utility.R) to
   take to sd(y) with dbartsSigmaFallbackWarning, exactly as the dense path does with more columns than rows.

No size cap, no message (Open calls). An allocation failure is an error, which
[`estimateStartingSigma`](../../R/spec.R) already turns into "unable to obtain a starting estimate of sigma;
provide one instead". Neither LAPACK call checks for an interrupt, as `lm.fit` does not.

## Change

1. R/utility.R: [`sparseResidualStandardError`](../../R/utility.R) and its comment rewritten as above (the
   batching, the refactor loop and Matrix's structural-rank muffler go);
   [`sparseDesignMatrix`](../../R/utility.R) gains the dgCMatrix wrap and the reference shift;
   [`estimateSigmaFromLinearModel`](../../R/utility.R) loses the fallback branch and its warning and sends
   every sparse source to the sparse routine; [`makeIndicatorModelMatrix`](../../R/utility.R) no longer sets
   `sparse.from.indicators`, its one reader gone.
2. R/xbart.R: [`xbart`](../../R/xbart.R)'s comment on the per-fold route drops "a sparse design"; in the
   chunk runner a sparse source's `sigmaDesign` is `sparseDesignMatrix(data@x)` once per chunk, and
   `foldData` fits `sparseResidualStandardError` on its training rows (`sigmaDesign[trainRows, ]`) where a
   dense source keeps `residualStandardError`. An indicator expansion stops being densified there.
3. The class dbartsSparseSigmaFallbackWarning is retired: it was new in 1.0 and nothing raises it.

## Constraints

- No engine, bridge, C API or state change. No new dependency: Matrix (Suggests) as today, base `chol`,
  `backsolve` and `qr`, all within R 4.2 and Matrix 1.4-1.
- Dense designs untouched: [`residualStandardError`](../../R/utility.R) and every dense caller bit for bit.
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
  - ["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R): perturbations 0, 1e-9,
    5e-8, 1e-4 and 1e-3 agree with the dense fit within 1e-8; 1e-6 pins the band: rank one fewer than
    lm.fit's and sigma equal within 1e-8 to the dense fit with the perturbed column removed. 2e-7 and 1e-5
    go (too near a boundary to pin).
- New, in test-indicator-storage.R or a new test-sparse-starting-sigma.R (under 2 s), each against the
  dense fit within 1e-10 unless said:
  - a bare dgCMatrix through `dbarts`, `bart` and `xbart`, with no warning (counted, Gate hygiene);
  - one-hot 3 x 40 plus three duplicated columns, weights with zeros, an offset;
  - a dgCMatrix with NA entries (the mean imputation);
  - sparseFactor columns with reference the first level, a middle level, and with missing values, in a
    frame with a numeric column (the shift);
  - p > n with rank below n - 1 (n 150, p 240, rank 41: the wide side); p > n of full rank: NA from the
    routine and, through `dbarts`, sd(y) with one dbartsSigmaFallbackWarning, as the dense design gives;
  - columns of mean 1e5 beside an exact duplicate: rank equal to lm.fit's;
  - xbart on a sparse frame: each fold's sigma equals the dense design's per-fold value for the same rows
    (the fold oracle's route, test-xbart-fold-oracle.R's helper or a direct `foldData` probe).
- test-predict-na-action.R and test-row-names.R reach the routine and pin nothing about sigma; they must
  pass unchanged.

Reviewer's mutants, each of which must fail a test: tolerance 1e-10 to 1e-16 (exact dependencies kept) and
to 1e-6 (the 1e-4 column dropped); the reference shift removed; the intercept column dropped; centering
reinstated in place of the intercept column (the mean-1e5 case); weights left out of B or of the residual;
the wide side's basis taken from the unpivoted factor; zero-weight rows kept; `n - r` replaced by `n - p`.

## Baselines

- Moves: equivalence.R's seven caller-sparse scenarios and wideFactorIndicators (Context; the last by
  rounding, unless the implemented routine happens to round to the QR's value). Current file
  equivalence-734441f1 ([MANIFEST](../../benchmarks/baselines/MANIFEST)).
- Does not move: the other 47 equivalence scenarios, bcf-equivalence-1b7d730c (15),
  multinomial-equivalence-80b1c8d4 (11), the four test-reproducibility files, every exact gate (Context). Any other mover is a defect: stop.
- Re-record on the reference build (`--preclean --configure-args=--enable-reference-build`):
  `EQUIVALENCE_SCENARIOS=sparse,mixedmatrix,sparsefactor,testswap,leaffactormixed,factorpartial,xbartmixed,wideFactorIndicators`
  (wideFactorIndicators dropped if it came out bitwise), `EQUIVALENCE_CORES=2`, merged into a copy of
  734441f1 in its scenario order, named after the slice's code commit; 734441f1 demoted to historical.
  Partition against 734441f1 in z mode: 47 (or 48) of 55 identical, the movers exactly those named,
  wideFactorIndicators with no |z| above 4. With the prior's anchor
  falling by factors of 0.2 to 0.7 on 150 to 500 rows, sigma and fit summaries past |z| 4 are expected and
  are not a failure; the partition's verdict is which scenarios moved. The merged file reproduces 55 of 55
  under `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17), as deb3fe50's row: the change is the starting value a prior is calibrated
  against, not the sampler; the identity is the routine's agreement with `lm.fit` (the tinytest pins and
  the 400-design sweep above, rerun by the implementer from scratch/ssp/sweep3.R's recipe against the
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
- Speed: the routine at n 1e5, p 5000 (1 percent) and n 1e4, p 1000 (one-hot 20 x 50) within 1.5x the
  Context table on the same machine; no bench-sampler.R compare (no hot path).
- The design note: [R surface](../design/sparse-columns.md#r-surface) in sparse-columns.md says "the sigma
  estimate falls back to sd(y)"; a dated paragraph there states the rule and Algorithm in brief, citing
  dec-B403 (the landing notes stay as written).

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd and man/xbart.Rd: the clause "a design with
  sparse-backed predictor columns skips the linear model altogether and falls back the same way (class
  dbartsSparseSigmaFallbackWarning, a dbartsSigmaFallbackWarning)" becomes: a sparse x takes the same
  estimate as its dense equivalent, fitted from its crossproduct without making it dense; its time grows
  with the cube and its memory with the square of the smaller of its row and column counts (seconds up to a
  few thousand, a couple of minutes and most of a gigabyte at ten thousand), so for a large sparse design
  supply `sigest`.
- man/bart.Rd's warning-class paragraph drops dbartsSparseSigmaFallbackWarning.
- man/sparseFactor.Rd: the paragraph on the default starting sigma is replaced by one sentence: a sparse
  column leaves the default starting sigma as the same column stored dense gives it.
- inst/NEWS.Rd: dbartsSparseSigmaFallbackWarning leaves the warning-class list (never released, so no
  entry of its own).
- docs/design/error-style.md: the class leaves the list and the table.
- At landing: TODO's item goes; ledger entry for the calls below; this plan's Status and Landing note;
  MANIFEST row.

## Steps

1. Algorithm and Change 1 to 3; the tinytest edits and new tests; the full suite green against
   `R CMD INSTALL -l <lib> .`.
2. The 400-design sweep against the implemented function; the speed points; help, NEWS and docs.
3. After review: the re-record, MANIFEST row and partition, in their own commit.

## Stop conditions

Stop and report when: the diff passes ~600 lines; an equivalence scenario other than the eight, a
snapshot or an exact gate moves, or wideFactorIndicators shows a |z| above 4; the sweep shows a rank different from the SVD's
or a sigma off by more than 1e-10 at more than 30 residual degrees of freedom outside the band; either speed
point is past 1.5x; the change needs engine, bridge or C API code.

## Interactions

- Whichever lands second of this and any other slice re-recording equivalence.R re-records against the
  other's file, and partitions against it.
- The categorical question below: if it is ruled before this lands, the sparseFactor scenarios move once.

## Calls made in planning

- One routine for every sparse design: the sparse QR goes, so an indicator expansion takes this path too
  (157 s against 0.2 s on one-hot 20 x 50 at n 1e4).
- The smaller Gram side, uncentered with the intercept as a column, tolerance 1e-10 on the unit-diagonal
  pivots, residual computed directly, no refinement.
- xbart's per-fold estimate uses the routine per fold, as the dense path refits per fold.
- dbartsSparseSigmaFallbackWarning retired.

## Open calls

- The band: accept that a column whose residual after the others is between 1e-7 and 1e-5 of its norm is
  dropped here and kept by `lm` (one residual degree of freedom each, and what that column fit; 0.08 to 0.9
  percent on the test design)? Recommended: in a crossproduct an exact dependency's pivot reaches 1.5e-13 at
  p 5000, above lm's 1e-14, so no tolerance near lm's separates the two, and matching lm exactly needs a QR,
  the cost this replaces.
- Large designs: no cap and no message (recommended: the ruling's text, and the dense path at the same size
  runs at least six times longer with no message), against a message once the Gram dimension passes 10000
  (about two minutes and 0.8 GB here; about 6 lines and a test), or sd(y) with the old warning past it.
- A categorical column enters the linear fit as its codes on the dense path in 1.0 (a factor kept as a
  factor; [`as.matrix.dbartsMixedMatrix`](../../R/mixedMatrix.R)), and this plan gives a sparseFactor the
  same. 0.9-34 expanded factors into indicators: on a 5-level factor and one numeric column (n 300, noise sd
  0.5) 0.9-34's dbarts starts sigma at 0.537 (lm with the factor) and 1.0's at 1.713 (lm on the codes). A
  separate call; if indicators, a sparse categorical column is expanded to indicator columns in the design
  (about 15 lines here) and the dense path changes with it.

## Estimate

Implementer about half a day; gates about two hours of machine time, the re-record minutes.
