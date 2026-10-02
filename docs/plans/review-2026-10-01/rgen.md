# Review 3, lens "rgen": post-fit generics, data/predict path, sparseFactor, tombstones, sampler user methods

Tree: review3 worktree at 01dee4b4, library r3-lib. Old behaviour from CRAN 0.9-34 in a private scratch lib.
Probe scripts: scratchpad/r3-rgen-p*.R (all small n, 5-20 trees, n.threads = 1).

Covered (with real fits):
- predict / fitted / residuals / extract (every type) / summary / print / plot / plotTree on bart (gaussian,
  weighted, offset, student, heteroscedastic, probit, probit + offset, weighted logistic, aft, bcf,
  monotone, linear leaf, hazard), bartMultinomial, bartOrdinal, bartNegbin, bartHurdle; 1 and 2-3 chains;
  combineChains TRUE/FALSE at fit time crossed with TRUE/FALSE at call time; ci.level; row names.
- predict(train rows) == stored draws for every class above (max diff <= 2e-14).
- newdata: reordered columns, missing column, character for factor, reordered/subset levels, unseen level,
  ordered factors, logical, sparseFactor test columns (reference and level order differing), na.action
  (default, na.omit, na.exclude, na.pass) with missing predictors and offsets, fit-time na.exclude padding.
- save/load (saveRDS in one process, predict/run in another) for bart fits incl. sparse columns; sampler
  copy() deep/shallow; sampler storeState/save/load continuation.
- Sampler user methods: run, setResponse, setOffset, setPredictor, setTestPredictor,
  setTestPredictorAndOffset, predict.
- Every tombstone in R/tombstones.R's registry exercised once (message and successor).
- sparseFactor vs factor semantics (36 operations) and dense vs sparse fits bitwise.
- Formula data path: subset with vector and column offset/weights, x/y interface with subset/na.action.

Not covered: xbart, rbart_vi internals beyond the tombstone, pdbart/pd2bart, diagnostics.R internals beyond
summary(), as_draws, the dbartsMixedMatrix S3 methods beyond what predict exercises.

## Findings

### rgen-01 - BLOCKER (pre-existing in 0.9-34, not a regression)
Location: R/data.R, dbartsData (formula branch, "## offset, when in data frame": model.offset is read only
when offsetGivenAsScalar is FALSE, i.e. only when an 'offset =' argument was also given).
Claim: an offset() term in a formula (the standard R idiom) is silently dropped on bart, dbarts and
dbartsData - neither applied as an offset nor used as a predictor - while the refusal text at the ':'
check says "poly(), ns(), log(), and offset() are supported".
Probe (r3-rgen-p3.R, r3-rgen-p29.R; y = x1 + 3 * off + noise):
```
fa <- bart(y ~ x1 + x2, df, offset = off, ...)          # rmse of fitted 0.95, sigma 1.06
fb <- bart(y ~ x1 + x2 + offset(off), df, ...)          # rmse 1.39, sigma 1.50
fc <- bart(y ~ x1 + x2, df, ...)                        # rmse 1.39, sigma 1.50  (fb == fc)
dbartsData(y ~ x1 + x2 + offset(off), df)@offset        # EMPTY
dbarts(y ~ x1 + x2 + offset(off), df)$data@offset       # EMPTY
dbartsData(y ~ x1 + offset(off), df, offset = rep(0, n))@offset  # = off: the term counts only beside an offset= argument
```
0.9-34 gives the same fb == fc. Why gates missed: no test fits a formula offset() term without an offset
argument; docs/plans/review-2026-08-24/memos/prerc-lens2-backlog.md P1 records "offset() all work" without
a probe.
Fix: read model.offset(modelFrame) whenever the frame has an offset term (sum with the argument as
model.offset already does), or refuse an offset() term by name; drop "offset()" from the "supported" text
if refused.

### rgen-02 - BLOCKER (new in 1.0-0)
Location: R/generics.R, predict.bartNegbin (type = "ppd" arm: negbinPpd(means, object$dispersion));
R/bart.R negbinPpd (as.vector(array(r, dim(mu)))).
Claim: predict(type = "ppd") on a negative-binomial fit pairs each mean-count draw with another draw's
dispersion whenever the caller's combineChains differs from the fit's stored layout (fit stored split +
default predict, or default fit + predict(combineChains = FALSE)); the means are paired correctly, only
the noise draw is not, so the posterior predictive is silently wrong.
Probe (r3-rgen-p15.R; fit with combineChains = FALSE, 2 chains, 4 draws):
```
fF$dispersion            # chains x samples: [5 5 5 4; 6 5 5 5]
set.seed(3); p <- predict(fF, nd, type = "ppd")
predict ppd == correctly paired draw: FALSE
predict ppd == mispaired draw:        TRUE
r paired with row 2 (chain 1, draw 2): used 6 should be 5
combined-stored fit, predict cC=F ppd correctly paired: FALSE
```
extract.bartNegbin(type = "ppd") is correct (it normalizes through scalarDrawVec); the default-fit +
default-predict path is also correct. Why gates missed: test-nbinom.R checks predict ppd shape and
non-negativity at one chain only; no cross-layout pairing check.
Fix: build the ppd in predict.bartNegbin from the split (chains x samples x obs) means with
dispersion.raw (as the means loop already does) and reshape afterward, or pass
scalarDrawVec(object$dispersion, n.chains, length(means)) after forcing means to the split layout.

### rgen-03 - MAJOR (new in 1.0-0: categorical columns are new)
Location: R/bartcore.R, bartcoreSamplerSetTestPredictor (column branch: new.x.test[, column] <-
as.double(x.test)) and bartcoreSamplerSetPredictor (same coercion for column updates).
Claim: giving $setTestPredictor or $setPredictor a factor for a categorical column installs its 1-based
R codes as the engine's 0-based category codes, so every value shifts up one level silently (and the top
level is refused with "categorical predictor values must be existing category codes"); a character
vector becomes NA with only a coercion warning. The man page says x is "a numeric predictor vector" and
documents no encoding for a categorical column, so a user has no documented way to do this right.
Probe (r3-rgen-p17.R, r3-rgen-p30.R; y = x1 + 3*(g == "b") - 3*(g == "c")):
```
setTestPredictor(df g = a,b,c)            -> test: -0.05  3.14 -3.01
setTestPredictor(factor a,a,b, column 'g') -> test:  3.10  3.10 -2.91     # predicted as b,b,c
setTestPredictor(char a,a,b, column 'g')   -> test: -1.71 -1.71 -1.71     # NA + coercion warning
setPredictor(newg, "g", forceUpdate = TRUE): labels b a a b b b -> codes 2 1 1 2 2 2  (0-based: c b b c c c)
```
Why gates missed: test-sparse-factor.R's column updates use numeric columns; whole-frame
setTestPredictor goes through validateXTest and is correct.
Fix: in both column branches, map a factor or character value through the training level table
(mapFactorColumnsToTrainingLevels / the sparse remap) for categorical columns, refusing unknown labels;
document the accepted forms.

### rgen-04 - MAJOR (new in 1.0-0)
Location: R/generics.R, extract.bartOrdinal / ordinalLogLik; R/bart.R, bart2Ordinal packaging (no
'active' element recorded).
Claim: an ordinal fit with 0/1 case weights installs the active-row mask (fit$fit$activeRows shows it),
but its packaged fit records no 'active' and extract(type = "loglik") reports finite log-probabilities
for the rows outside the likelihood, where a probit fit reports NaN; the shipped comment in
pointwiseLogLikelihood states ordinal's masked rows are "not in the model and [have] no likelihood to
report". loo/WAIC over this matrix silently counts rows the model never saw.
Probe (r3-rgen-p26.R; weights = rep(c(1, 0), 20)):
```
probit names: ... active
ordinal names: ... (no active)
probit loglik NaN cols: 20  of 40
ordinal loglik NaN cols: 0  of 40
```
Why gates missed: test-binary-weight-mask.R checks the NaN columns on the probit fit only; ordinal is
tested at the sampler level.
Fix: package result$active on ordinal fits (as bart's single-forest path does) and set those columns NaN
in ordinalLogLik.

### rgen-05 - MINOR
Location: R/data.R, validateXTest (the "does not match training's indicator columns" stop).
Claim: on a bart(factors = "indicators") fit, a test factor with a subset of the training levels is
refused with "... use bart() or dbarts(), which track levels across predict by default", though the
caller did use bart(); the remedy is factors = "categorical" (or a test factor carrying all training
levels).
Probe (r3-rgen-p7b.R):
```
fit <- bart(y ~ x1 + o + b + f, df, factors = "indicators", ...)
predict(fit, nd7)   # nd7$o <- factor("mid", ordered = TRUE)
Error: 'test' factor 'o' does not match training's indicator columns ('test' levels: mid); use bart() or
dbarts(), which track levels across predict by default
```
The refusal itself matches 0.9-34 (which failed with a bare column-count error). Why gates missed: the
message is pinned only for the bartBT/x-y route. Fix: say "keep the training levels on the test factor
(factor(x, levels = ...)), or fit with factors = \"categorical\"".

## Checked and found correct
- predict on training rows reproduces yhat.train (ev and bart scale) for 16 fit kinds x 2 chains;
  fitted == column means of extract; residuals == y - fitted, with fit-time na.exclude/default padding.
- Shapes and dimnames of fitted/residuals/extract/predict (ev, ppd, bart, loglik, sigma, k, varcount, ci)
  across all five classes at 1 and 2 chains, including the dec-A79 one-chain chain margin.
- extract/fitted/residuals/ppd/loglik identical between fits stored combined and split, all families
  (gaussian, weighted, student, hetero, probit, weighted logistic, multinomial, ordinal, nbinom, aft, bcf,
  k = chi()); predict ppd identical across storage for every family but nbinom (rgen-02).
- newdata column reordering, character for factor, reordered or subset levels, unseen-level refusal,
  ordered factors, logical columns, sparseFactor test columns, missing-column message.
- predict na.action arms match ?na.keepPredictors "At Prediction"; offset/weights NA routing (dec-A89).
- Multinomial predict offset: by-name permutation, data frame, unnamed, required-when-trained refusal.
- Formula subset with vector/column offset and weights, x/y subset and na.omit: bitwise to hand-subset fits.
- save/load: bart fits (incl. sparseFactor and sparseVector columns) predict and run after reload once
  storeState() is called (documented requirement); sampler copy() and reload continue to ~1e-15.
- sampler $predict / setTestPredictor(data frame) / setTestPredictorAndOffset / setOffset agree with run.
- sparseFactor: 36 factor operations match base factor; dense vs sparse fits bitwise (sigest fixed).
- Tombstones: every registry entry warns once, names the right successor and forwards or refuses as
  NEWS says; bart(x, y, x.test) note, fourth-positional refusal, BT-name forwarding.
- survivalProbabilities on hazard fits: test channel == newdata replay, chain shapes.
- 0.9-34 parity: predict without offset on an offset fit gives the offset-free surface in both versions;
  bartBT binaryOffset predict/yhat.test behaviour identical to 0.9-34.
