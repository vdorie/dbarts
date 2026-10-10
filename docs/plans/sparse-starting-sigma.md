# sparse-starting-sigma: a sparse x takes the linear-model starting sigma

Status: PLANNED 2026-10-09 (dec-B403); rewritten 2026-10-10 for dec-B427 to dec-B430. Not built. One call
is open ([Open calls](#open-calls)); the steps build its recommendation and it does not block the build.

agent: one sonnet implementer (R only, no engine code); one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for a continuous response under a chisq residual prior with no sigest, in two
cases: a caller's sparse x up to the cutoff (dec-B403), and any fit that keeps a factor as a categorical
predictor, the default, dense or sparse (dec-B422). SHIFTING for an indicator expansion R built sparse.
NEUTRAL for a dense design with no categorical column, a sparse design above the cutoff, a design with no
residual degrees of freedom, a fixed-unit family, a fit given sigest, and draws under a fixed residual prior.
window: before the merge (dec-B403). R only; may run beside an engine slice.
budget: ~800 lines (R/utility.R ~200 net of the QR's removal, R/spec.R ~35, R/xbart.R ~35, the data class,
packers and summary ~40, tinytest ~330, man ~40, docs ~60, NEWS ~4, MANIFEST one row).

## Goal

A sparse x gets the starting sigma its dense equivalent gets: `lm.fit`'s own while its dense form is small,
an exact sparse routine past that while the smaller of its row and column counts is at most 2,000, and the
response's sd above that. Every factor, in a dense design too, enters that regression as indicator columns.
Where the starting sigma is the response's sd because no linear estimate was made, the fit prints a line
under verbose and its summary says so; nothing warns, and both fallback warning classes and the sparse QR
go. An infinite entry in a sparse x is refused as in a dense one, and a fixed residual prior makes no
estimate. The tier is "Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

Ruled, each covering what it names:

- dec-B403: "Yes, just use the sparse regression, no warning." A sparse x takes its dense equivalent's
  linear-model starting sigma; a caller can pass sigest and the help says so.
- dec-B421: "Band is fine." The sparse routine's rank tolerance may differ from lm's on nearly dependent
  columns; after dec-B427 that reaches only designs past the dense bound.
- dec-B422: "What would keeping the codes mean? Treating them as integers? That doesn't make sense. I think
  indicators are necessary, right?" Every factor enters this regression as indicator columns.
- dec-B427: "Dense when small is fine." A design whose dense indicator form `lm.fit` handles in about a
  second and 200 MB takes `lm.fit`'s sigma as a dense caller does; a larger one takes the sparse routine, a
  zero-filled column beside the indicator of its unstored rows re-centered first. Nothing more is built.
- dec-B428: "Use sd(y) above the cutoff." A sparse design past the dense bound whose m is above 2,000
  starts from sd(y). No iterative solver, none of its build conditions, no 10 percent rule.
- dec-B429: "Then just use a verbose print, not a warning. And a line in the summary." Wherever the
  starting sigma is sd(y) because no linear estimate was made, dense and sparse alike.
- dec-B430: "Accept." The sparse routine keeps LAPACK's pivoted Cholesky; the help says a seeded fit
  reproduces at a fixed BLAS and thread count.

Today, run on bartcore 91173fe4 (2026-10-10, shipped build, arm64 macOS, R's reference BLAS):

- [`estimateSigmaFromLinearModel`](../../R/utility.R) returns sd(y - offset) with
  dbartsSparseSigmaFallbackWarning for any source with a CSC column; an indicator expansion R built sparse
  goes to Matrix's sparse QR; anything else to [`residualStandardError`](../../R/utility.R) (`lm.fit`), a
  factor kept categorical as its level codes. [`floorSigmaEstimate`](../../R/utility.R) turns a non-finite
  result into sd(y) with dbartsSigmaFallbackWarning (dense n 20, p 30: 0.866 and the warning).
  [`xbart`](../../R/xbart.R) catches that warning from its all-rows estimate to choose the per-fold route.
- The estimate runs once, in the calling process, for `dbarts`, `bart` (any chains and threads), `bartBT`,
  `xbart` and `dbartsSpec`, and three times for `rbart_vi` with two chains on one thread (a validation
  sampler and one per chain). The control there holds the door's verbose flag, except in `xbart`, which
  sets its control's to FALSE and keeps its own formal. `rbart_vi` has no summary method.
- No branch of stan4bart, bartCause, treatSens or bairrtt names either warning class or its text (`git
  grep` over every local branch); stan4bart writes the sigma slot itself and never reaches the estimate.

Measured the same day, the new rule emulated from `residualStandardError` and an R build of C and D:

- The dense bound, sparse designs just under it, densifying included: 0.98 s and 88 MB over R's heap at n
  1e4, p 345; 0.58 s at n 1100, p 1043; 0.42 s and 142 MB at n 1e5, p 59; 0.13 s and 183 MB at n 1e6, p 5.
  The sparse routine on the same four agrees within 4e-14, in 0.06 to 0.17 s.
- Each side of each boundary (today all are sd(y) with the warning). Dense bound: 0.508 from `lm.fit` just
  under it; its neighbor two columns wider 0.508 from the sparse routine, `lm.fit`'s value there. Cutoff,
  n 2600 at 1 percent: p 1999 (m 2,000) gets 0.496 in 0.9 s, `lm.fit`'s value; p 2000 gets sd(y) 1.121
  where `lm.fit` gives 0.508. No residual degrees of freedom: n 60, p 90 gets sd(y), its neighbor at p 40
  gets 0.386; past the bound n 1500, p 3000 of full rank gets sd(y), and at rank 1,001 0.480 in 2.9 s.
- The sparse routine against `lm.fit` on 400 random designs (n 30 to 2000, p 5 to 1530, thirds with
  dependent columns, zero weights, an offset): equal rank on 398, sigma within 9.1e-12 above 30 residual
  degrees of freedom. At m 2,000 it takes 1.0 to 1.1 s and 90 to 160 MB narrow, 3.9 s and 170 MB wide.
- Re-centering (C.4): a timestamp spread over an hour, 0 where missing, beside the indicator of those
  rows, at n 1e6 (past the bound): 1.761 unshifted, 0.4995 shifted, `lm.fit` 0.4995, in either column
  order. With no indicator, or one that differs in one row, nothing is shifted and the design is
  identical. A partner test keyed on squared row indices is inexact past n 3e5 and missed that design.

## What moves

Every call of [`estimateSigmaFromLinearModel`](../../R/utility.R) logged through a load hook beside
`lm.fit`'s value on the indicator design; the suite then rerun with that value substituted and both
warnings muffled. Every design reached is under the dense bound.

- [equivalence.R](../../benchmarks/R/equivalence.R), quick mode: 132 calls in 44 of 55 scenarios. Twelve
  move, the same twelve as before dec-B427, each to `lm.fit` on its dense indicator design:
  - seven on a caller's sparse source, sd(y) to the linear estimate (0.16 to 0.62 of today's value):
    sparse, mixedmatrix, sparsefactor, testswap, leaffactormixed, factorpartial, xbartmixed;
  - four on a dense frame with a factor, codes to indicators, 0.9-34's value (0.08 to 4.1 percent):
    categorical, leaffactor, nafactor, ordfactor;
  - wideFactorIndicators, the sparse QR's value to `lm.fit`'s (6.6e-16 apart).
- tinytest: 20055 results in 247 files, 4,643 calls. Six expectations fail in four files and a fifth file
  stops; a sixth file is found by reading ([Tests](#tests)).
- Not rerun since the rebase, which changed none of the files they read: bcf-equivalence.R (12 calls, all
  dense, no factor), multinomial-equivalence.R (none), the four test-reproducibility files (no sparse or
  factor design), exact-gates.yaml's list (5 calls on a factor design, all under a fixed residual prior).

## Design

`estimateSigmaFromLinearModel(data)` returns `list(sigma, route)`, route one of "dense", "sparse", "df" (no
residual degrees of freedom) and "size" (above the cutoff). Three named constants in R/utility.R, each
the default of a formal of the routing function so a test reaches any route at any size:
`sigmaDenseCells <- 6e6`, `sigmaDenseWork <- 1.2e9`, `sigmaExactCutoff <- 2000L`.

A. The design, built by one new function, `startingSigmaDesign(x)`, used by the all-rows estimate and by
xbart's chunk runner. dec-B422's expansion happens here only; `data@x` and the trees never see it.

1. A plain matrix with no factor marked: [`sigmaDesignMatrix`](../../R/utility.R) as today. A matrix whose
   `varTypes` mark a factor is rebuilt as a frame (codes to factors over its `factor.levels`, or over the
   codes present where no table rides) and takes case 2.
2. A dense container with a factor column: the frame is rebuilt from its column list and names and handed
   to [`makeIndicatorModelMatrix`](../../R/utility.R) with `drop = TRUE` and storage "auto", so a frame,
   its matrix of codes and its `factors = "indicators"` fit share one sigest bit for bit. The result is a
   matrix (case 1) or a container with sparse-built indicators (case 4).
3. A caller's sparse source. A bare dgCMatrix is wrapped ([`wrapSparseTestMatrix`](../../R/mixedMatrix.R)).
   A dense-backed factor column becomes one indicator per present level but its first; a sparseFactor
   column one per present level but its reference, from its stored entries, a stored entry at the
   reference code dropped; a stored NA is an NA in each of the factor's indicators. Other columns as in
   [`sparseDesignMatrix`](../../R/utility.R) today, then its imputation (for an indicator, the level's
   frequency). Columns are assembled from slots as [`assembleMixedMatrix`](../../R/mixedMatrix.R) does.
4. An indicators-route container: taken as it is. Every indicator column, from either builder, records the
   term that emitted it (`indicator.term`, replacing `sparse.from.indicators`, whose one reader goes).

B. Route. n is the rows with a response, weight and offset and a positive weight; p the built columns.

1. An infinite entry is refused first, as for a dense design, naming the column.
2. A design that is a dense matrix goes to [`residualStandardError`](../../R/utility.R) at any size.
3. A design with a sparse-stored column, R's own indicators included, with n (p + 1) <= `sigmaDenseCells`
   and n (p + 1) min(n, p + 1) <= `sigmaDenseWork` is made dense and goes there too: `lm.fit`'s rank rule
   and bits (to rounding where A.3 built a factor without one level). Route "dense".
4. Otherwise the front end (C), then m = min(n, p + 1) with p counted after C.2. m <= `sigmaExactCutoff`:
   the exact routine (D), route "sparse". Above it no fit is made: route "size".
5. A non-finite sigma from 2, 3 or D is route "df".

A design with m above 2,000 is always past the dense bound (2001^3 is above 1.2e9), so the regions nest.
Across the dense bound two neighbors differ by rounding, or by dec-B421's band where columns are nearly
dependent; across the cutoff the step is from the linear estimate to sd(y).

C. Front end of the sparse routine, in this order.

1. Rows: drop rows with a missing response, weight or offset and rows of weight 0; z = y - offset.
2. Constants: a column whose entries over the kept rows are all equal (implicit zeros included; an exact
   comparison) is dropped.
3. Blocks, counted before any centering: the columns of one `indicator.term` are a full block when at
   least two are left and they sum to one in every kept row (within 1e-8).
4. Re-centering. On the assembled design, where a dense-backed column's zeros are unstored: a column with
   unstored kept rows and stored entries that are not all equal has the mean of its stored entries taken
   off them when the design also holds a column constant on exactly its unstored rows and unstored
   elsewhere, or constant on exactly its stored rows. Partners are found by stored count and the sum of
   stored row indices (exact in doubles), then confirmed row for row. Otherwise the column is left alone.
5. Centering: each dense-backed column and each column stored in every kept row is centered at its weighted
   mean over the kept rows. Then every column is divided by its largest absolute entry.

D. The exact routine.

1. B = diag(sqrt(w)) [1 X], each column divided by its norm. Sparse columns are not centered: a centered
   crossproduct cancels for a large-mean column, which is why the intercept stays a column.
2. Narrow (ncol(B) < n): `chol(as.matrix(Matrix::crossprod(B)), pivot = TRUE, tol = 1e-10)` inside
   `suppressWarnings` (its one warning is the rank deficiency); r its rank, K its first r pivots, R its
   leading r x r block; e = z sqrt(w) - B_K R^-1 R^-T B_K' (z sqrt(w)), computed directly.
3. Wide: K = `as.matrix(Matrix::tcrossprod(B))`, equilibrated by D = sqrt(diag(K)), the pivoted Cholesky of
   D^-1 K D^-1 at the same tolerance. r >= n gives no estimate; else Q = `qr.Q(qr(D L, LAPACK = TRUE))`, L
   the first r columns of the factor, unpivoted (LAPACK's pivoted QR applies no rank tolerance, so the rank
   is not decided twice). e = z sqrt(w) - Q Q' z sqrt(w).
4. sigma = sqrt(sum(e^2) / (n - r)); no estimate when n - r <= 0.

Accepted past the dense bound (dec-B421, dec-B427): a near-copy of a column that carries signal (132 to 143
percent high), two zero-filled times sharing one missing pattern with no indicator (120 percent), rows
weighted 1e12 times the rest (35 percent).

E. No estimate (routes "df" and "size"; dec-B429).

1. Value: sd(y - offset) through [`floorMarginalSigma`](../../R/utility.R), as today. No warning.
2. Record. `dbartsData` gains a slot `sigma.fallback` beside `sigma`: NULL, "df" or "size". One setter
   writes both: the estimate sets the record, a supplied sigest clears it, and
   [`setData`](../../R/bartcore.R)'s copy of the sigma carries it. A data object saved before the slot
   existed lacks it, so every read goes through a guarded reader, as [`dataRowNames`](../../R/data.R)
   does. `bart`, `bartBT` and `rbart_vi` copy it to the fit as `sigest.fallback`, beside `sigest`.
3. Verbose. A helper in [`announceAutoFamily`](../../R/utility.R)'s form: a message of class
   dbartsStartingSigmaMessage (inheriting dbartsMessage), sent when the flag is TRUE:
   - "df": `starting sigma is the sd of the response (<value>): the linear model has no residual degrees
     of freedom; supply 'sigest' to set it`
   - "size": `starting sigma is the sd of the response (<value>): no linear estimate is made for a sparse
     design above 2000 rows and columns; supply 'sigest' to set it`

   | door | flag (default) | prints |
   |---|---|---|
   | `dbarts`, `dbartsSpec` | the control's `verbose` (FALSE) | once, from [`resolveSamplerSpec`](../../R/spec.R) |
   | `bart`, `bartBT` | `verbose` (TRUE) | once, the same place, whatever the chains and threads |
   | `rbart_vi` | `verbose` (TRUE) | once, for its validation sampler; [`rbart_vi_fit`](../../R/rbart.R) muffles the class; nothing where `rbart_vi` turns verbose off (several chains on several threads) |
   | `xbart` | its own `verbose` (FALSE) | once, in the calling process, for the all-rows estimate; folds print nothing |

4. Summary. [`summary.bart`](../../R/diagnostics.R) carries `sigest.fallback` and
   [`printSummaryBartBody`](../../R/diagnostics.R) prints, after the fixed line, `(Starting sigma: the sd of
   the response; the linear model had no residual degrees of freedom)` or `(Starting sigma: the sd of the
   response; no linear estimate for a sparse design above 2000 rows and columns)`. A hurdle fit prints it
   under its positive part. A fit saved without the record prints no line, as one saved without `fixed`
   names nothing fixed. `rbart_vi` and `xbart` have no summary to print in.

F. A fixed residual prior makes no estimate. In [`resolveSamplerSpec`](../../R/spec.R), where the residual
prior is a dbartsFixedPrior and the family is not on a fixed unit scale, an NA slot takes the square root
of the fixed variance and a slot that holds a value keeps it (an agreeing sigest, or what the data object
carried). [`xbart`](../../R/xbart.R) resolves its prior before its all-rows estimate and does the same,
with no per-fold estimate. [`refuseSigestUnderFixedPrior`](../../R/family.R) still runs first, untouched.

G. xbart. The all-rows estimate runs once, to raise a refusal and send the message; its warning handler
goes. `startingSigmaDesign` runs once per chunk in [`xbartRunUnits`](../../R/xbart.R), for dense sources
too, and `foldData` routes each fold's training rows through B by that fold's own n, p and m, taking the
fold's sd under "df" or "size". Where the all-rows route is "df" every fold takes its sd without a fit (a
subset of independent rows has no residual degrees of freedom either). A fold may take another route than
all rows do: fewer rows can put it under the dense bound or under the cutoff.

## Change

By file: R/utility.R takes A to D and the helper of E.3 ([`sparseResidualStandardError`](../../R/utility.R)
is rewritten, its batching, refactor loop and `grepl` muffler going; [`floorSigmaEstimate`](../../R/utility.R)
loses its warning). R/spec.R: [`nonFinitePredictorNames`](../../R/spec.R) reads sparse sources (the dense
list and the CSC block's stored entries, named through the container's map, or positions for a bare
dgCMatrix); [`estimateStartingSigma`](../../R/spec.R) passes the route on; `resolveSamplerSpec` records,
announces and applies F. R/A_class.R and R/data.R: the slot, its setter and reader, through which
R/bartcore.R and R/dbarts.R write. R/bart.R and R/rbart.R: `sigest.fallback` on the fit, the muffle.
R/diagnostics.R: E.4. R/xbart.R: F and G. Both warning classes are retired: never in a release, and no
consumer names them, so the consumers change nothing.

Constraints: no engine, bridge, C API or state change. No new dependency: Matrix (Suggests) as today; the
calls made (`crossprod`, `tcrossprod`, `rowSums`, `cbind2`, `Diagonal`, `%*%`, subsetting) are exported by
Matrix 1.4-1, R 4.2's (read from its NAMESPACE, not run). A dense design with no categorical column under
a chisq prior is bit for bit unchanged, at any size. Out of scope: a size rule for dense matrices, a
summary method for `rbart_vi`, any iterative solver.

## Tests

Edits:

- [test-starting-sigma.R](../../inst/tinytest/test-starting-sigma.R): three `expect_identical` pins against
  `lm` on a factor's codes; the design becomes `makeModelMatrixFromDataFrame` of the same frame, bitwise.
- [test-data-mixed.R](../../inst/tinytest/test-data-mixed.R), the block counting
  ["dbartsSparseSigmaFallbackWarning"](../../inst/tinytest/test-data-mixed.R): no warning for the sparse
  frame or its dense equivalent, the two `sigest` within 1e-10.
- [test-indicator-storage.R](../../inst/tinytest/test-indicator-storage.R):
  ["a sparse column the caller supplied still falls back"](../../inst/tinytest/test-indicator-storage.R)
  becomes no warning and the dense fit's value; ["collinear to 1e-9 and 5e-8"](../../inst/tinytest/test-indicator-storage.R)
  pins every perturbation to the dense fit with `expect_identical` (under the dense bound).
- [test-boundary-inputs.R](../../inst/tinytest/test-boundary-inputs.R), the n 2 and n 1 fits: no fallback
  warning (the n 1 fit keeps its other one); `data@sigma` as pinned; the reader returns "df".
- [test-xbart-fold-oracle.R](../../inst/tinytest/test-xbart-fold-oracle.R),
  ["falls back to the marginal response sd"](../../inst/tinytest/test-xbart-fold-oracle.R): silent by
  default, the message under `verbose = TRUE`.
- [test-sampler-splitProbabilities.R](../../inst/tinytest/test-sampler-splitProbabilities.R): a seeded
  comparison of split counts on a factor design fails at its seed. The implementer reports the counts
  over five seeds; if it holds in distribution the seed moves, and if not it is a stop.

New, in test-sparse-starting-sigma.R (under 10 s), against `lm.fit` on the dense indicator design, missing
values imputed as the estimate imputes them, within 1e-10 unless said; "forced" is both dense constants 0:

- dec-B422: a frame with a 5-level factor and a numeric column, default against `factors = "indicators"`
  and its matrix of codes: `sigest` identical and equal to `summary(lm(y ~ f + x))$sigma`; the same with
  the factor ordered, with two levels, with 40, with a missing value; `data@x` and its `varTypes` as in
  the fit given sigest. A frame with a sparseVector column and a dense factor; a sparseFactor with its
  first and a middle level as reference, and a middle reference with missing values; each also forced.
- Under the bound: a bare dgCMatrix through `dbarts`, `bart` and `xbart` with no warning (counted, Gate
  hygiene), no message under verbose, route "dense", `sigest` identical to the dense matrix's; the
  zero-filled timestamp with its indicator in both orders, a near-copy with signal and a row at weight
  1e12, each `identical` to `lm.fit`.
- Forced: one-hot 3 x 40 plus three duplicated columns with zero weights and an offset; NA entries; a fully
  stored timestamp-like column carrying signal, weighted and not; a column scaled by 1e160; p > n with rank
  below n - 1; p > n of full rank (route "df").
- Re-centering, forced: the timestamp pair in both orders and a value stored on one level's rows beside
  that level's indicator, within 1e-8; a column with no partner, the indicator itself, a column whose
  stored entries are all equal, and one whose partner matches in count and index sum but not in rows leave
  the design `identical` to the one with re-centering off. Unforced at n 1e6 with 7 columns: shifted.
- Route, read from the value returned: a sparse design just under and just over each dense constant is
  "dense" and "sparse", within 1e-10; a dense matrix past both is "dense"; n 2600 at 1 percent is "sparse"
  at p 1999 and "size" with `sd(y)` at p 2000. Forced, with the cutoff formal: n 40, p 100 at cutoff 45 is
  exact (m is n); 20 live and 30 constant columns at n 200, cutoff 25, is exact.
- dec-B429, for "df" (dense and sparse) and "size", by `dbarts`, `bart`, `bartBT`, `rbart_vi` (two chains,
  one thread) and `xbart`: zero warnings; no message with the flag FALSE and exactly one
  dbartsStartingSigmaMessage naming the reason with it TRUE; `sigest` is `sd(y)`; `sigest.fallback` and
  the summary's line. No line for a linear estimate, a supplied sigest, or a fit with the field removed;
  a hurdle fit prints it under its positive part; after `setData` the record follows the sigma.
- An Inf entry in a dgCMatrix and in a sparseVector column: the dense path's error, naming the column.
- xbart: on a sparse frame and on a dense frame with a factor each fold's sigma equals `lm.fit`'s on the
  indicator design for its rows; through the constants, a design "size" on all rows and "sparse" per fold,
  and one "sparse" on all rows and "dense" per fold; an all-rows "df" design fits no fold (traced).
- A fixed residual prior, dense and sparse: `estimateStartingSigma` is not called; the slot and `sigest`
  equal the fixed sigma; an agreeing sigest and a carried sigma stay; draws identical to the build before.

test-data-code-channel.R passes unchanged and is the pin for A.1. The 37 files that reach a factor design
must otherwise pass unchanged; a failure there that is not a pinned sigest, slot or seeded draw of a
default fit with a factor is a stop.

Reviewer's mutants, each of which must fail a test: the expansion removed (codes) for a container, for a
matrix of codes, and in `foldData` alone; an ordered factor left as codes; a reference-level indicator
added to a sparseFactor; `indicator.term` not recorded by the caller's-sparse builder; imputation before
the expansion; the size route removed; each dense constant's `<=` inverted; a dense matrix sent through
the size route; cutoff 2000 to 20000, and `<=` to `<`; m taken as p + 1, and before constants are dropped;
a fold routed by the all-rows size; the constant test replaced by a sum-of-squares one; blocks counted
after centering; the re-centering without its partner test, without the row-for-row confirmation, keyed
on squared indices, and applied to a column whose stored entries are all equal; weights left out of the
crossproduct, of the means, of the residual; zero-weight rows counted in n; tolerance 1e-10 to 1e-16 and
to 1e-6; centering removed, and at the unweighted mean; the max-abs step removed; the equilibration
removed; the wide basis without D, and from `qr` at its default tolerance; `n - r` replaced by `n - p`; the
Inf check removed; a warning left in `floorSigmaEstimate`; the message sent with the flag FALSE, and once
per chain in `rbart_vi`; the record not written, written for a linear estimate, kept when sigest is
supplied, and dropped by `setData`; the summary line printed without the record; the fixed-prior skip
leaving the slot NA, and overwriting a carried value.

## Baselines and gates

- Current: equivalence-e4faed5c, bcf-equivalence-1b7d730c, multinomial-equivalence-80b1c8d4
  ([MANIFEST](../../benchmarks/baselines/MANIFEST)). POSTERIOR-CHANGING the seven caller-sparse and the
  four dense-factor scenarios; SHIFTING wideFactorIndicators; NEUTRAL the other 43, bcf's 15,
  multinomial's 11, the four snapshot files and every exact gate. Any other mover is a defect: stop.
- Re-record the twelve on the reference build (`--preclean --configure-args=--enable-reference-build`,
  `EQUIVALENCE_CORES=2`), merged into a copy of e4faed5c in its scenario order and named after the slice's
  code commit; e4faed5c demoted to historical. Partition against it in z mode: 43 of 55 identical, the
  movers exactly the twelve, no |z| above 4 on wideFactorIndicators. The merged file reproduces 55 of 55
  under `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17): the change is the value a prior is calibrated against, not the sampler; the
  identity is agreement with `lm.fit` on the indicator design (the tinytest pins, and the 400-design sweep
  rerun against the implemented function, forced, with factor columns added to a third of its designs).
- On the slice tip against its own library, independently of the implementer
  ([RNG classes and their gates](README.md#rng-classes-and-their-gates), posterior-changing): tests/cpp;
  the full tinytest suite; the four snapshot files on the reference build unchanged; bcf and multinomial
  bitwise; exact-gates.yaml's list in quick mode, output as before (the fit carries a new field);
  `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness. No
  sanitizers and no bench-sampler.R compare (no compiled code, no hot path).
- Speed and memory, same machine, within 1.5x: `lm.fit` on a sparse design at the dense bound (n 1e4, p
  345: 1.0 s, 90 MB over the heap); the exact routine at n 2e4, p 1999, 1 percent (1.1 s, 160 MB) and wide
  at n 2000, p 3600, rank 901 (3.9 s, 170 MB).

## Help and docs

- `sigest` in man/bart.Rd, man/bartBT.Rd, man/dbarts.Rd, man/rbart.Rd and man/xbart.Rd: a sparse x takes
  its dense equivalent's estimate, `lm`'s own while the design is small and a sparse regression past that;
  a sparse design above 2,000 rows and columns, like any design with no residual degrees of freedom,
  starts from the response's standard deviation, which a verbose fit prints and `summary` shows; supply
  `sigest` to set it. A factor enters as indicator columns whatever `factors` says. Under a fixed residual
  prior no estimate is made. A seeded fit reproduces at a fixed BLAS and thread count (dec-B430).
- man/xbart.Rd: each fold is routed by its own size. man/bart.Rd: the warning-class paragraph drops both
  classes, the value section lists `sigest.fallback`, `summary`'s text names the line. man/dbartsData.Rd
  lists the slot. man/dbartsSpec.Rd's "an unset value is still estimated" gains "under a chisq prior".
  man/sparseFactor.Rd's starting-sigma paragraph becomes one sentence.
- inst/NEWS.Rd: both classes leave the warning-class list; one item: a design with no residual degrees of
  freedom starts from the response's standard deviation, printed under verbose and shown by `summary`.
- docs/design: sparse-columns.md's [R surface](../design/sparse-columns.md#r-surface) gets a dated
  paragraph with the rule; error-style.md drops both classes and lists the message class;
  memory-footprint.md's starting-sigma row gains the sparse-source routes; starting-sigma-sensitivity.md's
  [Recommendation](../design/starting-sigma-sensitivity.md#recommendation) a dated line for dec-B428.
- At landing: TODO's item goes; a ledger entry for the calls below; this plan's Status and Landing note;
  the MANIFEST row.

## Steps

1. The changes above with Design A to G; the tinytest edits and new tests; the suite green against
   `R CMD INSTALL -l <lib> .`.
2. The 400-design sweep against the implemented function, forced; the speed and memory points.
3. Help, NEWS and docs.
4. After review: the re-record, MANIFEST row and partition, in their own commit.

Stop and report when: the diff passes ~1100 lines; an equivalence scenario other than the twelve, a
snapshot or an exact gate moves, or wideFactorIndicators shows a |z| above 4; the sweep shows a rank
different from `lm.fit`'s on a design other than its two wide ones, or a sigma off by more than 1e-10 above
30 residual degrees of freedom; a speed or memory point is past 1.5x; a tinytest outside the six files
fails; the change needs engine, bridge or C API code. Whichever lands second of this and any other slice
re-recording equivalence.R re-records against the other's file.

## Calls made

- A dense matrix keeps `lm.fit` at any size: dec-B427's "a larger one takes the sparse routine" is read
  with dec-B428's "a sparse design past the dense bound", so the size rule applies to a design with a
  sparse-stored column. Sending large dense matrices to the sparse routine is 3 lines and moves the last
  digit of every dense fit past the bound.
- The dense bound is tested on the built design's columns and the cutoff after constant columns are
  dropped: the first prices the matrix `lm.fit` is handed, the second is the ruled m.
- The record is a slot on the data object, where the sigma it describes lives and every door's packer
  reads; it names the reason so the two lines can. A fit-only field is about 10 lines fewer and is lost by
  `dbarts` and by `setData`.
- The verbose line is a classed message in `announceAutoFamily`'s form, so `rbart_vi`, which builds a
  sampler per chain, mutes the repeats by class. A plain `cat` prints three times at two chains.
- `rbart_vi` gets the record and the verbose line but no summary line: it has no summary method, and one
  is its own slice (about 60 lines with its help).
- xbart announces for its all-rows estimate only. Reporting folds that fall back alone means returning
  each unit's route from the workers and printing a count, about 20 lines.
- Each xbart fold routes by its own size. Routing every fold as all rows are routed is 2 lines fewer and
  gives sd(y) to folds under the cutoff when all rows are above it.
- The re-centering takes either partner, the indicator of a column's unstored rows (dec-B427's words) or a
  column constant on exactly its stored rows; both leave the span unchanged, and a value recorded for one
  level beside that level's indicator needs the second (1.017 against `lm.fit`'s 0.505 without it).
- A factor is expanded wherever the estimate meets one, a matrix of codes included, by the indicators
  route's own builder, so the doors share one sigest. Under a fixed residual prior only an empty slot is
  filled: overwriting it breaks what `dbartsSpec` documents and what stan4bart writes.

## Open calls

### 1. Does a plain data frame whose factors expand past the cutoff start from sd(y)?

Background. dec-B422 expands every factor into indicators for this estimate, and R stores an indicator
column sparse when its level holds at most 20 percent of the rows. dec-B427 sends a design "sparse or not"
past the dense bound to the sparse routine; dec-B428 names "a sparse design" for sd(y) above the cutoff. A
plain data frame with more than 2,000 rows and factors of more than 2,000 levels in all (an identifier, a
postal code) is not sparse to its caller. Run at n 5000 with one factor of 4,600 levels and a numeric
column, noise sd 0.5: `lm.fit` on the dense indicator matrix, 0.9-34's route, gives 0.480 from 399 residual
degrees of freedom in 53 s on 176 MB; today's build, regressing on the level codes, gives 1.124; sd(y) is
1.513.

Options. (a) sd(y), as the steps are written: one rule by the size of the built design, with the line
under verbose and in summary; three times 0.9-34's value on the design above. (b) `lm.fit` on the dense
indicator matrix at any size for a frame the caller gave dense, as 0.9-34: about 3 lines; 53 s there, and
by the same rate about 10 minutes and 2.4 GB at n 1e5 with 3,000 levels (not run), again in every xbart
fold. (c) The exact routine at any m for such a frame: about 3 lines; about 2 minutes and 1.8 GB at m
10,000 (dec-B426's measurement) and hours past that.

Recommended: (a). The study behind dec-B428 found a starting sigma 2 to 4 times too high moving the fit at
n 200 and not measurably at n 5000 ([Results](../design/starting-sigma-sensitivity.md#results)), and
every such design has more than 2,000 rows.

## Estimate

About 800 lines. Implementer one to two days: the design builder with the caller's-sparse expansion is the
largest and riskiest part (half a day); the front end and exact routine half a day against the R build in
hand; the record, message and summary a few hours; tests and docs half a day. Gates about two hours, the
re-record minutes. Review with the mutants half a day, and one or two fix rounds.
