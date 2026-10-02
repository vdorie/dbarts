# Review 3 - lens: aux (wave two)

Tree: review3 worktree at 01dee4b4, library r3-lib. Old behaviour from CRAN 0.9-34 (scratch lib
r3-rgen-cranlib). Probe scripts: scratchpad/r3-aux-p*.R, small n, 3-30 trees, n.threads <= 2.

Covered, with real fits:
- xbart: fold construction (k-fold, LOO, random subsample by fraction and by count, one held-out row),
  offset and factor/logical response reaching the loss, k-order invariance, full grid
  (n.trees x k x power x base) shape and labels vs 0.9-34, drop, multi-valued loss, sd axis vs k axis
  (probit and logistic anchors, leaf.prior normal(sd =)), control precedence for n.samples/n.thin, subset,
  NA predictors/response, n.threads 1 vs 2 under fork and socket, custom loss environments, leakage of the
  held-out response.
- pdbart/pd2bart: against a direct predict()-based computation per draw (1 and 2 chains, gaussian and
  probit, keepTrees TRUE and FALSE, with offsets), categorical and sparseFactor columns, two-predictor
  models, reloaded fits, a keepTrees = FALSE bart fit.
- rbart_vi: fitted/extract/predict agreement on train and test (1 and 2 chains, combineChains both ways,
  offset, probit), unseen groups (factor, character, numeric), type = "ranef", user tau priors, seeded
  stream restore (serial and PSOCK), keepTrees = FALSE, ppd sigma pairing across layouts, saveRDS/readRDS
  into a new session then predict, extract(type = "trees"), plot.
- summary(): sigma and varcount rhat/ess against posterior 1.7.0 summarise_draws (1 and 3 chains, odd and
  even n.samples, combined and split storage) - exact to 4e-15.
- extract(type = "loglik") against a hand computation from extract(type = "ev") and the stored sigma:
  weighted gaussian with offset, probit with offset, weighted logistic, student; combined and split.
- varcount against split counts read off extract(type = "trees").

Not covered: plot.bart and the other family plot methods (wave one rendered them), aft/hazard/hetero
loglik, statistical calibration of xbart losses, any C++ path.

## aux-01 - MAJOR (regression vs 0.9-34)
Location: R/xbart.R, xbart (`loss <- list(loss, evalEnv)`) and xbartLossFunction
(`environment(result) <- loss[[2L]]`).
Claim: a custom loss function has its enclosing environment replaced by the caller's frame, so any
variable it captured from a factory or local scope is lost: the call fails, or, when the caller's frame
has a variable of the same name, the loss silently uses that value instead.
Probe (r3-aux-p20.R):
  makeLoss <- function(mult) function(y.test, s, w) mult * sqrt(mean((y.test - rowMeans(s))^2))
  xbart(x, y, ..., loss = makeLoss(10), seed = 1)        # 1.0: ERROR: object 'mult' not found
  mult <- 1; xbart(x, y, ..., loss = makeLoss(10), seed = 1)   # 1.0: 0.2363561 (= plain rmse, factor 10 lost)
  0.9-34, same script: 2.033674 both times (10 x its rmse 0.2033674).
0.9-34 evaluated the call in the given environment and left the closure intact
(R_interface_crossvalidate.cpp, CustomLossFunctor).
Why gates missed it: every test loss is defined at the top level of its test file, which is also the
caller's frame.
Fix: call a bare function as given, keeping its own environment; for the list form, evaluate the call in
the supplied environment rather than re-parenting the function.

## aux-02 - MAJOR (new in 1.0-0: categorical columns are new)
Location: R/partialDependence.R, pdbart.defaultLevs, pdbart and pd2bart (x.test[, xind[j]] <- levs[[j]][i]);
plot.pdbart.
Claim: on a categorical column (bart()'s default factors = "categorical", and sparseFactor columns)
pdbart and pd2bart treat the 0-based category codes as numbers: default levs are quantiles of the codes,
which either fail inside the engine or silently skip levels, and results and plots are labelled by codes,
never by level names; levs cannot be given as level names.
Probe (r3-aux-p8.R, r3-aux-p11.R, r3-aux-p7.R):
  15-level factor, n = 300: default levs 0 1 2 3 5 7 9 10 12 13 14 - levels 4, 6, 8, 11 are absent, and
    plot.pdbart draws a line through the codes as if ordered. PD values at the codes it did pick match a
    direct predict() computation.
  12-level balanced factor, n = 60: quantile levs 0 1 2 3 4 5.5 7 ...;
    pdbart(fit, xind = "f") -> ERROR: categorical predictor values must be existing category codes
    (same with keepTrees = FALSE).
  3-level factor lo/mid/hi: pd$levs[[1]] is 0 1 2.
Why gates missed it: test-pdbart*.R use numeric predictors only.
Fix: for a categorical column default levs to every category code, carry the level names into the
result (levs labels, plot axis), and accept levs given as level names.

## aux-03 - MAJOR (pre-existing silent collapse; the loud failure for named inputs is new)
Location: R/partialDependence.R, pd2bart, both `ncol(sampler$data@x) == 2L` branches.
Claim: on a model with exactly two predictors pd2bart averages over the whole grid, so fd is
n.draws x 1 instead of n.draws x (length(levs[[1]]) * length(levs[[2]])), and plot.pd2bart then errors;
in 1.0 a named predictor matrix, a formula, a data frame or a bart() fit now fails earlier, because the
grid is passed with expand.grid's Var1/Var2 names.
Probe (r3-aux-p10.R, r3-aux-p9.R):
  pd2bart(unnamed X (2 cols), y, xind = c(1, 2), pl = FALSE, keeptrees = TRUE, ...)
    dim(fd) = 50 1, levs lengths 11 11; plot(): ERROR: dimensions of z are not length(x)(-1) times
    length(y)(-1). Same on 0.9-34. A third column gives 50 x 121 and the right surface.
  pd2bart(X with colnames x1, x2, y, ...), pd2bart(bart fit, ...):
    1.0: ERROR: column names of 'test' do not match those of 'x': 'x1, x2' present in 'x' but not in
    'test' (whose columns are 'Var1, Var2'); 0.9-34 returned the collapsed 20 x 1.
With pl = FALSE and unnamed predictors the wrong result is silent.
Why gates missed it: test-pdbart.R runs pd2bart only on data with more than two predictors.
Fix: in that branch each grid point is an observation, so fd is the draws x grid prediction itself
(chains combined), not its per-draw mean; give x.test the training column names.

## aux-04 - MINOR (pre-existing; 7ad0bbea had fixed it, the port restored 0.9-34's test)
Location: R/rbart.R, rbart_vi (`if (is.symbol(matchedCall$prior) || ...) prior <- rbart.priors[[which(...)]]`).
Claim: prior = myPrior, a user function referred to by name, fails because every symbol is looked up in
rbart.priors; rbart.Rd documents "A function or symbolic reference to built-in priors".
Probe (r3-aux-p5.R):
  prior = myPrior            -> ERROR: attempt to select less than one element in get1index
  prior = function(x, rel.scale) ..., prior = lst$p, prior = gamma, prior = "gamma" -> run
  prior = "nope"             -> ERROR: could not find function "prior"
  0.9-34 failed on the first four (only the string form ran).
Why gates missed it: test-rbart-options.R passes built-ins and a literal function.
Fix: take the built-in only when the symbol or string names one (7ad0bbea's builtinTauPrior test);
otherwise evaluate the argument; refuse an unknown string by name.

## aux-05 - MINOR (pre-existing)
Location: R/rbart.R, predict.rbart (`ranefNames.test <- levels(group.by)`).
Claim: with a character or numeric group.by at predict (both accepted by rbart_vi at fit), type = "ranef"
returns an empty matrix and an unseen group errors instead of drawing from tau as a factor does.
Probe (r3-aux-p3.R):
  predict(f, nd, group.by = c("a","b","c"), type = "ranef")  -> num[1:100, 0]
  predict(f, nd, group.by = c("a","b","zz"))                  -> ERROR: subscript out of bounds
  predict(f2 (numeric group.by fit), nd, group.by = c(1, 2, 9)) -> ERROR: subscript out of bounds
  predict(f, nd, group.by = factor(c("a","b","zz")))          -> warns, draws "zz" from tau
Why gates missed it: rbart predict tests pass factors.
Fix: coerce group.by with as.factor at entry.

## aux-06 - MINOR (pre-existing)
Location: R/xbart.R, xbart (`data@sigma <- estimateStartingSigma(data)` on the full data); every fold's
sampler reads data@sigma for its chisq residual prior (bartcore_createFromHandle, bartcore_setModel).
Claim: the default sigest is one least-squares fit over all rows, held-out responses included, so a
fold's posterior depends on its own held-out y through the residual prior's calibration - the leakage
xbart.Rd's n.burn entry says the design avoids.
Probe (r3-aux-p30.R), k-fold, n.reps = 1, same seed, only y[1] changed (+50):
  fold 3 holds out row 1; its test draws identical after changing y[1] only: FALSE
  same, with sigest = 0.2 fixed: TRUE
Why gates missed it: test-xbart-fold-oracle.R checks row membership, not what a fold's fit reads.
Fix: estimate sigest per fold from its training rows when the caller gave none (or document it).

## Checked and found correct
- xbart: fold sizes and membership; offset reaches both the fit and the test draws (gaussian and probit);
  factor response coded 0/1 at the loss; k-order invariance; grid placement and labels match 0.9-34;
  sd axis equals the k axis at anchor 3 (probit) and pi*sqrt(3) (logistic); subset equals pre-subsetting;
  NA predictors kept, NA response dropped; LOO and one-row subsample; identical results at 1 thread,
  2 fork and 2 socket workers; control n.samples precedence.
- pdbart/pd2bart (three or more predictors): per-draw values equal a direct predict() computation to
  1e-14 for 1 and 2 chains, gaussian and probit; keepTrees TRUE and FALSE agree with an offset; sparseFactor
  column PD at 0/1/2 matches direct predict.
- rbart_vi: predict(train rows) == extract to 6e-15 for 1 and 2 chains, both combineChains; test-set
  extract == predict with an explicit offset.test (gaussian, probit); saveRDS/readRDS then predict identical;
  PSOCK fit predicts; seeded stream restored serial and PSOCK; ppd noise sd tracks the paired sigma
  (<= 3 percent) for every fit/predict combineChains pair; extract(type = "trees") with chainNums/treeNums.
- summary(): sigma and varcount statistics equal posterior 1.7.0 summarise_draws to 4e-15.
- extract(type = "loglik"): exact against hand computation for weighted gaussian + offset, probit +
  offset, weighted logistic, student; combined and split.
- varcount equals split counts from the saved trees.

Not reported (pre-existing in 0.9-34, deprecated port kept as 0.9-34): plot.rbart fails ("invalid 'xlim'")
on a multi-chain fit with n.burn = 0; fitted.rbart fails on a keepTrainingFits = FALSE fit; pdbart on a
keepTrees = FALSE bart fit runs the user's sampler further and resets its test predictor without the warning
the dbartsSampler branch gives. A reloaded keepTrees fit needs fit$fit$storeState() before saving for
pdbart as for predict (docs-03).
