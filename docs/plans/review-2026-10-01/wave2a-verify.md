wave2a-verify - independent verification of ingest.md (9) and aux.md (6), review-3, pinned 01dee4b4

Method. Own probes scratchpad/r3-verify-wave2a-p1..p16.R against r3-lib (1.0-0 at 01dee4b4) and
r3-rfit-cranlib (CRAN 0.9-34; bart2 is its formula/x door). Code read in src/bartcore/data.hpp,
sampler.hpp, src/R_interface_bartcore.cpp, R/data.R, R/utility.R, R/xbart.R, R/partialDependence.R,
R/rbart.R; checked root TODO, docs/decisions.md (through dec-B164), docs/design, man/. None of the 15
is filed in TODO or ruled in the ledger; dec-B156 (setter refuses numbers on a categorical column) and
dec-B155 (indicator route maps labels) are the nearest rulings and cover neither ingest-04's entrance
nor any other finding here.
Result: 12 CONFIRMED, 3 QUALIFIED (ingest-08 broader, aux-04 broader and 0.9-34 detail wrong,
ingest-05 severity), 0 refuted.

Summary
  id         verdict    sev      vs 0.9-34                  surface decision?          moves draws?
  ingest-01  CONFIRMED  BLOCKER  regression                 no (option A)              no (option A)
  ingest-05  QUALIFIED  MINOR    65535 part new; NaN same   no                         no
  ingest-02  CONFIRMED  BLOCKER  new feature                no                         only the broken fits
  ingest-03  CONFIRMED  MAJOR    regression                 no (restores 0.9-34)       no
  ingest-07  CONFIRMED  MINOR    int refusal same as 0.9-34 yes, small (accept more)   no
  ingest-04  CONFIRMED  MAJOR    0.9-34 errored             yes (new refusal)          no
  ingest-06  CONFIRMED  MINOR    same                       yes, small                 only fits with Inf
  ingest-08  QUALIFIED  MINOR    partly same                no                         no
  ingest-09  CONFIRMED  MINOR    behaviour same; doc new    yes, small                 no
  aux-01     CONFIRMED  MAJOR    regression                 no (restores 0.9-34)       no
  aux-02     CONFIRMED  MAJOR    new (categorical new)      yes (levs content)         no
  aux-03     CONFIRMED  MAJOR    same silent; loud new      no                         no
  aux-04     QUALIFIED  MINOR    same                       no                         no
  aux-05     CONFIRMED  MINOR    same (0.9-34 worse)        no                         no
  aux-06     CONFIRMED  MINOR    same                       yes (xbart results move)   xbart snapshot

=============================================================================================
G1 - the cut-grid validator (ingest-01, ingest-05): one helper, three call sites

ingest-01  CONFIRMED  BLOCKER (loud, but a saved fit is unusable)  regression
  Code: data.hpp fillCutsOverRange writes xMin + (k+1) * (xMax - xMin)/(n+1), so a zero range gives
  n.cuts equal values; sampler.hpp Sampler::setState refuses `cutPoints[j][k] <= cutPoints[j][k-1]`
  (and installForests' cross-grid branch does the same). The engine already knows the uniform
  degenerate grid repeats a value (refreshCutsForColumn's comment); test-sampler-degenerate-cuts.R's
  comment ("the single degenerate cut a uniform grid would place") is false.
  Own probe p1/p2 (bart, 2 chains, keepTrees, storeState, saveRDS/readRDS, predict):
    1.0-0  constant 0 column        live 40x3, reloaded: "state is not consistent with this sampler"
           all-TRUE logical         same refusal
           constant on subset 1:50  same refusal
           x = 1e15 + 4 runif       same refusal (ulp-collapsed cuts)
           useQuantiles = TRUE      reload ok
           dbarts(cbind(a, b = 0)); run; storeState: $copy() and $setState($state) refused;
           with updateState = TRUE (the default) $copy() after any run is refused
    0.9-34 every case reloads and predicts 40x3.
  Options. A: the validator accepts a non-decreasing grid (equal neighbours allowed), which is
  exactly what the store builds; no draws move, no surface. B: build a degenerate uniform column
  as one cut, as the quantile path does; moves draws for every fit with a constant column (the
  cut-index draw changes), and does not cover the ulp case, which needs a dedupe as well.
  Recommendation A.
ingest-05  QUALIFIED  MINOR (filed MAJOR)
  Confirmed (p6, 600 rows, 60 NA, y = 5 where x < .1 or NA): 65534 cuts: hi/lo/mid/NA 0.01 5.02
  0 4.99; 65535 cuts: 2.51 4.85 0.17 2.59 (top bin pooled with NA). setCutPoints(c(.2, NaN, .6))
  and setCutPoints(NaN) accepted; 0.9-34 accepts both too but had no NA code. Silent, but only an
  explicit setCutPoints of exactly 65535 cuts or a NaN cut reaches it: MINOR by reachability.
Fix (one helper, e.g. bartcore::cutGridIsValid(const double*, size_t n, bool strict)): n in
  [1, maxNumCutsRepresentable]; no NaN (test `!(c[k] >= c[k-1])` / `!(c[k] > c[k-1])`, and isnan on
  c[0]). Call non-strict from Sampler::setState and installForests' cross-grid check; strict from
  bartcore_setCutPoints (replacing its 65535 literal). Tests: test-sampler-degenerate-cuts.R adds
  the uniform constant column (save/load predict, copy, setState) and the 1e15 + runif column, and
  fixes the comment; setCutPoints refuses 65534+1 cuts and NaN; tests/cpp setState on a repeated grid.

=============================================================================================
ingest-02  CONFIRMED  BLOCKER  new in 1.0-0 (sparse frame columns are new)
  Code: R/data.R dbartsData formula branch, pos <- match(rownames(modelFrame), rownames(data));
  model.frame renames repeated subset rows "1.1", which match nothing.
  Own probe p3 (n 200, bootstrap subset, 77 repeats): sparse column NA cells 77, dense 0; values
  agree elsewhere; unique, logical subsets 0; row-named data 77; bart() fitted values differ from
  the dense-column fit by up to 4.01 (truth: 5 where s > 0).
  Fix: carry a row index through model.frame instead of matching names: add an extra model.frame
  argument (e.g. `(dbartsRow)` = seq_len(NROW(data))), which subset and na.action then shape exactly
  as weights, and index subsetSparseColumn by it. Tests: bootstrap subset, sparseVector, dgCMatrix
  and sparseFactor columns, x equal to the dense frame's; with na.omit as well.

=============================================================================================
G2 - column kinds refused at ingestion (ingest-03, ingest-07): one coercion pass in the R doors

ingest-03  CONFIRMED  MAJOR  regression
  Code: R/utility.R makeCategoricalModelMatrix tests is.numeric(column), FALSE for Date and
  difftime (is.numeric.Date/.difftime) and POSIXct; a classed integer passes.
  Own probe p4: 1.0-0 Date, POSIXct, difftime: formula and x/y doors "column 'z' cannot be converted
  to a predictor"; factors = "indicators" and bartBT fit; 0.9-34 fits all three (cor .96-.98).
  No doc intends it; lm() takes a Date as its number.
  Fix: in makeCategoricalModelMatrix's numeric branch accept any non-factor atomic column of type
  double or integer, as.double(unclass(column)); POSIXlt stays refused (message could suggest
  as.POSIXct). The test/predict frame goes through the same builder, so a Date test column follows;
  test fit and predict with Date, POSIXct, difftime on both doors.
ingest-07  CONFIRMED  MINOR  integer refusal also in 0.9-34 ("x must be of type real")
  Own probe p8: 1.0-0 integer matrix/vector "'x' must be numeric" (bart, dbarts, bartBT); logical
  matrix, dgTMatrix, dgRMatrix "unrecognized 'formula' type"; dgCMatrix ok; predict() takes an
  integer matrix but refuses a logical one ("test matrix must be numeric").
  Decision (small widening): accept integer and logical matrices (storage.mode<- "double") and any
  Matrix class (as(x, "CsparseMatrix") then to double) in dbartsData's x branches and predict's
  validateXTest, as lm's model.matrix accepts them; the cost is a slightly wider R input surface.
  Recommendation: accept.

=============================================================================================
ingest-04  CONFIRMED  MAJOR  0.9-34 errored on the same call; silent miscoding is new
  Code: R/utility.R mapFactorColumnsToTrainingLevels skips any column that is not factor or
  character, so numbers reach the engine as 0-based codes.
  Own probe p5 (levels a-d, effects 0/10/20/30): factor a,b,c -> 0 10 20; integer 1:3 -> 10 20 30;
  double 0:2 -> 0 10 20; logical TRUE -> 10; 1:4 or 1.5 error; bart(test = 1:3) and an ordered factor
  given 1:3 shift the same way. predict.lm refuses this (.checkMFClasses: fitted with type "factor"
  but type "numeric" was supplied).
  Decision. A: refuse a numeric or logical data-frame column whose training column has a level
  table (categorical or ordered), naming the column - the frame-entrance twin of dec-B156, which
  covers only $setPredictor/$setTestPredictor. B: read numbers as 1-based codes (as.integer(factor)).
  C: keep 0-based, documented. A costs nothing a user should rely on; B and C keep a silent trap for
  the other convention. Recommendation A; a raw numeric matrix (no frame) is left as is.
  Fix: in mapFactorColumnsToTrainingLevels, after the sparseFactor branch, stop() when the column is
  neither factor nor character; tests for predict, bart(test =), dbartsData(test =), $setTestPredictor
  with a frame, categorical and ordered.

=============================================================================================
ingest-06  CONFIRMED  MINOR  same in 0.9-34
  Own probe p7 (z step at .5, z[1] = Inf, sigest given): uniform cor .14 / .22 (two cols), sigma 1.0;
  quantiles cor .975, sigma .23. Without sigest both versions stop with "unable to obtain a starting
  estimate of sigma" (lm on Inf), naming the wrong input.
  Decision (small). A: take the uniform range over finite values (fillCutsUniformly and the CSC twin
  skip !isfinite when seeding and folding), so +-Inf code past the ends, as the quantile path already
  behaves; moves draws only for fits containing Inf. B: refuse non-finite predictors by column, as
  lm.fit refuses. Recommendation A (trees are rank-based; the two grid modes then agree), plus naming
  the non-finite column in the sigest failure.

ingest-08  QUALIFIED  MINOR  partly pre-existing
  Broader than filed (p9): with a dbartsData object, factors, subset, weights and na.action are all
  ignored without the warning (weights = rep(2, n), subset = 1:10: 50 rows, no warning); 0.9-34 also
  ignored subset and weights silently. Fix: add every supplied argument (missing() tests for subset,
  weights, factors, na.action) to the existing warning in dbartsData's dbartsData short-circuit.
  Under dec-B161 that warning should become a plain warning() (it carries dbartsIgnoredArgWarning).

ingest-09  CONFIRMED  MINOR  doc sentence false; recycling and truncation as in 0.9-34
  p10: n.cuts = 5:7 on two columns -> 5 6, silently. Decision (small): A, correct bart.Rd to
  dbartsControl.Rd's wording; B, also refuse a vector longer than ncol(x) in dbarts()/bart()/xbart()
  (its extra entries can never apply). Recommendation A plus B's refusal of the too-long case only.

=============================================================================================
aux-01  CONFIRMED  MAJOR  regression
  Code: R/xbart.R wraps a bare function as list(loss, evalEnv); xbartLossFunction then runs
  environment(result) <- loss[[2L]], dropping the closure's own scope (introduced dd092cb0). 0.9-34
  built a call to the closure and evaluated it in the environment, leaving the closure intact.
  Own probe p11: makeLoss(10) factory: 1.0-0 "object 'mult' not found"; with mult <- 1 in the calling
  function 0.2178 (= plain rmse, factor 10 lost, silent); list(makeLoss(10), globalenv()) also fails.
  0.9-34: 2.2444 in all three.
  Fix: never re-parent. Bare function: use as is. List form: function(y.test, s, w)
  eval(as.call(list(fn, y.test, s, w)), env). Tests: factory closure and a same-name shadow inside a
  function, n.threads 1 and 2 (PSOCK serializes the closure's environment).

=============================================================================================
G3 - partial dependence (aux-02, aux-03): one pass over R/partialDependence.R

aux-02  CONFIRMED  MAJOR  new (categorical route is the 1.0 default)
  Own probes p12, p16: default levs are quantiles of 0-based codes (15 levels: 1 3 5 6 7.5 ...),
  so pdbart(fit, xind = "f") errors "categorical predictor values must be existing category codes";
  a 12-level ordered factor errors "ordered factor predictor values must be existing level codes";
  levs = list(c("A", "B")) errors "test matrix must be numeric". Ordered factors are affected too.
  Decision. A: for a factor column default to every level, report levs as level names, accept names
  in levs, and plot per level (points, no line); pd2bart's grid likewise. B: refuse factor columns,
  pointing to predict(). C: all codes, labelled by names only. A matches how partial-dependence tools
  treat factors; B loses a working feature for the default fit. Recommendation A.
  Fix: pdbart.defaultLevs reads the column type and level table from sampler$data; map names to
  codes before x.test[, j] <- code; carry names to buildResult and plot.pdbart/plot.pd2bart.
aux-03  CONFIRMED  MAJOR  silent collapse same as 0.9-34; the named-input error is new
  Code: pd2bart, both ncol == 2 branches pass the whole grid as x.test and then pdbart.drawMeans
  averages over it, giving one column.
  Own probe p13: unnamed 2 columns fd 20 x 1 (levs 11, 11); named x1, x2: 1.0-0 error "column names
  of 'test' do not match ... (whose columns are 'Var1, Var2')", 0.9-34 20 x 1 with a warning;
  3 columns 20 x 121.
  Fix: in those branches fd is the per-draw prediction itself, t(pred) with chains combined
  (draws x grid), not its mean; set colnames(x.test) to the training names. Test: a two-predictor
  model's fd equals a direct predict at each grid point, both keepTrees settings, and plot runs.

=============================================================================================
G4 - rbart_vi (deprecated, removed in 1.1-0)

aux-04  QUALIFIED  MINOR  pre-existing
  Code: R/rbart.R rbart_vi computes builtinTauPrior correctly (line ~213) and resolves built-ins,
  then a second block (~297) runs rbart.priors[[which(names == matchedCall$prior)]] for ANY symbol.
  Own probe p14, same results on 1.0-0 and 0.9-34: prior = myPrior fails (get1index); lst$p, gamma,
  "gamma" run; "nope" fails with "could not find function 'prior'". Broader than filed: a wrapper
  passing its own argument (prior = pr) fails even for "gamma". The filed claim that 0.9-34 failed
  on lst$p and gamma is wrong here.
  Fix: delete the second block (the first already resolved built-ins); refuse an unknown string by
  name listing names(rbart.priors). Test: symbol, wrapper argument, unknown string.
aux-05  CONFIRMED  MINOR  pre-existing (0.9-34 also fails the character ev case)
  p15: character group.by: type = "ranef" 2 x 0, unseen level "subscript out of bounds"; a factor
  works. Fix: group.by <- as.factor(group.by) after the length check in predict.rbart.

=============================================================================================
aux-06  CONFIRMED  MINOR  pre-existing (0.9-34 ran the same lm on all rows)
  Code: R/xbart.R `data@sigma <- estimateStartingSigma(data)` once on all rows before the folds;
  every fold calibrates its chisq prior from it, so a fold's fit reads its own held-out y. Mild:
  one scalar calibration, not the fit. dec-A05 settled the chain-carry leakage only.
  Decision. A: when sigest is NULL estimate it per fold from the fold's training rows (one lm per
  fold); moves every default xbart result and the test-reproducibility-xbart.R snapshot. B: keep,
  and say in xbart.Rd's sigest entry that the default is estimated once on all rows.
  Recommendation A (cheap, and it is what cross-validation promises); B if the snapshot move is
  unwelcome this close to the RC.

Not re-verified: ingest-01's quantile-midpoint overflow near 1.7e308 (p15 of the ingest lens) and
aux-06's own probe (the code argument is direct).
