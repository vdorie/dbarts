rfit lens - R fitting entry points and argument handling (review-3, pinned 01dee4b4)

Covered: bart() (formula and x/y doors, legacy forwarding to bartBT, bart2 alias,
retired flat names power/base/split.probs/sigdf/sigquant/resid.prior/proposal.probs,
control precedence), formula/data/subset/weights/offset/offset.test/na.action handling
against lm/glm, factor/ordered/character/logical/Date predictors, unseen and
declared-unused levels, NA in x and y at fit and predict, test-set column matching,
data-dependent formula terms, every family token (gaussian, student, probit, logistic,
binomial()/gaussian() objects, multinomial, ordinal, nbinom, aft, hazard, hurdle,
heteroscedastic, monotone, blocks, interactions, dart, linear leaves) and their
refusals, seed/set.seed reproducibility across n.threads 1 vs 2 for every family plus
xbart and rbart_vi, keepTrees/keepSampler/keepCall/keepTrainingFits/keepFits/
combineChains and the generics on the result, xbart argument sweep, rbart_vi basics,
output shapes against CRAN 0.9-34 (installed in a scratch lib). Every probe ran
against the r3-lib build of 01dee4b4; old-behaviour claims ran against CRAN 0.9-34.

Not covered: statistical correctness of any posterior (no SBC/exact work), the
sampler reference-class methods, dbartsSpec/multi-forest formula terms beyond a
smoke fit, sparse columns, warm.start internals, plot output, pdbart.

Probe scripts: scratchpad r3-rfit-p*.R (header r3-rfit-hdr.R defines A(), a
match.call wrapper that fills n.trees = 5, n.samples = 10, n.burn = 5, n.chains = 1,
verbose = FALSE when not given; tr() turns errors/warnings into strings).

-----------------------------------------------------------------------------

rfit-01  BLOCKER
Location: R/data.R dbartsData (formula branch, "offset, when in data frame" block,
  offsetGivenAsScalar); same path serves bart, bartBT, dbarts, xbart.
Claim: an offset() term in the formula is silently dropped - the model frame carries
  it, but model.offset() is read only when the 'offset' argument was also supplied, so
  y ~ x + offset(log(exposure)) fits as if no offset were written.
Probe (r3-rfit-p6.R; O = 100 * rnorm, y = a + O + noise):
  bart(y ~ a + offset(O), data = x, ...)   cor(fitted, y) 0.347, sigma 86.4
  bart(y ~ a, offset = O, data = x, ...)   cor(fitted, y) 0.9999, sigma 1.04
  xbart(y ~ a + b + offset(o), ...) returns the same loss as xbart(y ~ a + b, ...)
  (r3-rfit-p19.R: offsetTerm mean 1.128 == basic mean 1.128); an NA inside the
  offset() column is also accepted silently (r3-rfit-p5.R), where offset = o with an
  NA errors "'offset' contains missing values".
  0.9-34 drops it the same way (not a regression), but 1.0 now advertises it: the
  interaction refusal in dbartsData says "poly(), ns(), log(), and offset() are
  supported", and the nbinom docs describe the offset as a log-exposure, the case
  where offset(log(exposure)) is the idiomatic glm spelling.
Why gates missed it: inst/tinytest/test-bart-formula.R ("poly()/log()/offset() terms
  are unaffected - each still fits") and test-formula-terms.R
  (expect_silent(fit(y ~ a + offset(o)))) assert only that the call returns; no test
  compares an offset() fit with the offset= fit.
Fix: when 'offset' is missing, take model.offset(modelFrame) (and add it to an
  explicit offset when both are present, as lm does); or refuse offset() by name.

rfit-02  BLOCKER
Location: R/data.R validateXTest, inner replayTerms (rebuilds a formula from
  term.labels and calls model.frame on the new data).
Claim: data-dependent terms (poly(), ns(), scale(), bs()) are recomputed on the test
  set / predict newdata instead of reusing the training basis (terms' predvars), so the
  prediction at a given x depends on which other rows are in newdata - silently wrong.
Probe (r3-rfit-p17.R; y = 3 c^2 + N(0, 0.1^2), truth f(0.5) = 0.75; fitted near
  c = 0.5 is about 0.70-0.76):
  term        predict f(0.5): newdata nd1   newdata nd2   lm predict nd1/nd2
  poly(c, 2)                  1.478         -0.022        0.753 / 0.753
  scale(c)                    0.797          0.194        0.967 / 0.967
  ns(c, 3)                    0.762          0.013        0.735 / 0.735
  log(c)                      0.794          0.794        (stateless: correct)
  (nd1 = c(0.5, 0.9, 0.1, 0.3), nd2 = c(0.5, 0.51, 0.52, 0.53)); test = nd2 at fit
  time gives the same wrong value. A one-row newdata with poly(c, 2) errors
  "'degree' must be less than number of unique points". 0.9-34 behaves the same (not a
  regression), but 1.0 names poly() and ns() as supported in its own refusal message.
Why gates missed it: the only poly() test (test-bart-formula.R) checks the fit's class;
  the predict tests replay the training rows or the fit-time test, which reproduces
  yhat.test bit for bit while both are wrong.
Fix: keep attr(terms(modelFrame), "predvars") (or the terms object) on the data and
  build the test frame with model.frame(delete.response(terms), newdata), as predict.lm
  does.

rfit-03  MAJOR
Location: R/data.R dbartsData, explicit offset.test branch
  (offset.test <- rep_len(offset.test, nrow(x.test))) and getTestOffset symbol lookup.
Claim: an explicit offset.test of the wrong length is silently recycled or truncated
  to the test row count, and a bare column name resolves in the TRAINING data first, so
  offset.test = o picks the first m training offsets even when 'test' carries its own
  column o.
Probe (r3-rfit-p30n.R; 60 training rows, 5 test rows, test$o = 100):
  offset.test = c(7, 8, 9)   -> stored 7 8 9 7 8, no error
  offset.test = O (a 60-row training column)  -> stored O[1:5], no error
  offset.test = o (test has o = 100)          -> stored training o[1:5]; yhat.test
                                                 ~1, not ~100
  0.9-34 does the same. bart.Rd's offset.test item says a length that does not match
  "is refused by name rather than silently recycled" (for the default), and
  predict()'s offset, by contrast, is length-checked - the fit-time door is the
  lenient one.
Why gates missed it: no test passes an explicit offset.test of the wrong length or a
  column name present in both data and test.
Fix: refuse a non-scalar explicit offset.test whose length is not nrow(test); resolve a
  bare name in 'test' before 'data' (or refuse the ambiguity).

rfit-04  MAJOR
Location: R/data.R factor coding on the factors = "indicators" route
  (makeModelMatrix / indicator drop) and bartBT, which always takes that route.
Claim: under factors = "indicators" a level with no training rows (declared-unused, or
  emptied by subset) loses its indicator column, so predict at that level is silently
  predicted as another level - the 0.9-34 behaviour NEWS ("A factor level that training
  never saw is an error (0.9-34 silently predicted it as another level)") and dbartsData
  ("A level declared but unobserved in training therefore keeps its own bin") say is
  gone.
Probe (r3-rfit-p38.R):
  bart(y ~ a + b, subset = b != "c", factors = "indicators"): varcount columns a, b.b;
    predict at b = a, b, c (a = 0) -> 0.920 2.095 0.920  (c == a exactly)
  bartBT on the same rows: columns a, b.b -> 0.903 1.977 0.903
  declared levels a,b,c,z with z unused, indicators: columns a, b.a, b.b, b.c; z is
    accepted and routed as the all-zero pattern -> 2.938
  The categorical default keeps the bin (documented), and a character column refuses
  the unseen level; only the indicators route still mispredicts silently.
Why gates missed it: the unseen-level tests exercise the categorical route and
  character columns; none subsets away a level or declares an unused one under
  "indicators".
Fix: on the indicators route, keep the level table and refuse (by name) any predict
  level whose indicator column was dropped or never existed in training.

rfit-05  MINOR
Location: R/generics.R predict.bart, the keepTrees-free "current trees" path.
Claim: a multi-chain fit with keepSampler = TRUE and keepTrees = FALSE fails predict
  with "internal error: row names do not match the observations"; a one-chain fit
  silently returns a single current-state draw (a length-n vector) instead of draws.
Probe (r3-rfit-p23.R, r3-rfit-p24.R):
  n.chains = 1: predict -> Named num [1:3] (one draw); n.chains = 2: "internal error:
  row names do not match the observations" (also for type = "bart" and ci.level);
  r$fit$predict returns an n x n.chains matrix that nameObservationMargin rejects.
  0.9-34 returned uninitialized memory (1.55e-313) for the two-chain case, so 1.0 is
  better, but the message is an internal error.
Why gates missed it: no test predicts from a keepSampler-only multi-chain fit.
Fix: refuse as the ppd and amplitude arms already do ("predict requires the fit's
  saved trees") whenever keepTrees is FALSE, or at least when n.chains > 1.

rfit-06  MINOR
Location: R/rbart.R rbart_vi (PSOCK/FORK cluster branch).
Claim: rbart_vi results depend on n.threads and, with n.threads > 1 (the default on a
  multi-core machine, min(guessNumCores(), n.chains)), set.seed does not make the run
  reproducible; rbart.Rd says n.threads and seed are "Same as in bart", whose
  Reproducibility section promises both.
Probe (r3-rfit-p11.R): set.seed(11) twice, n.threads = 2 -> identical FALSE;
  seed = 7, n.threads 2 vs 1 -> identical FALSE (seed = 7 twice at 2 threads -> TRUE).
  0.9-34 port behaviour (dec-B130), deprecated function.
Why gates missed it: rbart tests run n.threads = 1.
Fix: state in rbart.Rd that its thread count changes draws and that only an explicit
  seed reproduces a multi-threaded run (or draw per-chain seeds from R's stream when
  seed is NULL).

rfit-07  MINOR
Location: R/data.R restateMatrixResponseError.
Claim: bart(cbind(time, status) ~ x, family = "aft" or "hazard") is refused with a
  message that contradicts itself, while the x/y door accepts the same two-column
  matrix under those families.
Probe (r3-rfit-p40.R): formula cbind(tm, st) ~ a + c, family = "aft" -> 'family = "aft"
  takes a single-column response but 'y' is an n x 2 matrix - ... a (time, status) pair
  needs a survival::Surv response ... or family = "aft" / "hazard"'; bart(X, cbind(time,
  status), family = "aft") fits.
Why gates missed it: refusal tests match the message prefix only.
Fix: accept the two-column LHS on the formula door under an explicit survival family,
  or drop 'or family = "aft" / "hazard"' from the reading when that family was given.

rfit-08  MINOR
Location: R/bart.R bart (combineChains is not validated before use).
Claim: combineChains = NA fails with the bare base error "missing value where
  TRUE/FALSE needed"; every sibling logical (verbose, keepTrees) is refused by name.
Probe (r3-rfit-p14.R): combineNA : ERROR: missing value where TRUE/FALSE needed;
  keepTreesNA : 'keepTrees' must be TRUE/FALSE.
Fix: validate combineChains with the same TRUE/FALSE check.

rfit-09  MINOR
Location: R/generics.R summary.bart call preamble; bart keepCall = FALSE stores
  call("NULL").
Claim: summary() of a keepCall = FALSE fit prints "Call: `NULL`()" (print omits the
  call correctly), and update() fails with 'could not find function "NULL"'.
Probe (r3-rfit-p21.R): summary : |Call:|`NULL`()|...; update : could not find
  function "NULL".
Fix: have summary use the same "no call kept" test print uses; refuse update() by
  name on a call-less fit.

-----------------------------------------------------------------------------

Checked and found correct
- Seeds: bart with seed = and with set.seed gives identical draws at n.threads 1 vs 2
  for gaussian, student, logistic, multinomial, ordinal, nbinom, aft, hazard, hurdle,
  heteroscedastic, monotone, 3 chains/2 threads, n.grow.sweeps, warm.start; xbart too;
  seed = leaves .Random.seed untouched; R's stream after a NULL-seed fit does not
  depend on n.samples; predict(n.threads) and ppd under set.seed thread-invariant.
- predict reproduces yhat.test (combined and uncombined) for gaussian, student,
  multinomial, ordinal, nbinom, aft, heteroscedastic, monotone, dart, linear leaves,
  offset fits.
- subset (index, negative, logical with NA, expression), weights/offset columns and
  vectors with subset, NA weights/offset refused or dropped per na.action, na.exclude
  pads fitted to n as lm does, na.omit/na.fail, predict na.action arms, test NA
  refusal naming the columns, training+test NA routing.
- Test/newdata column matching by name (matrix and data frame), reordered levels,
  level subsets, character test column for a factor, missing column named.
- Legacy door: BayesTree names forward to bartBT with a warning, including from a
  function's local frame; fourth positional argument refused; bart2 forwards; retired
  power/base/split.probs/sigdf/sigquant/resid.prior/proposal.probs give draws
  identical to their new spellings; flat-over-control precedence as documented.
- Family resolution: binomial()/binomial(probit)/gaussian()/student(df) objects,
  unsupported links and families refused by name; probit/logistic/nbinom response and
  weight validation; survival type refusals; monotone (numeric, ordered factor,
  transformed term, NA column, probit, dart, indicators) holds on a predict grid,
  refusals for unordered factors, linear leaves, variance forest; blocks and
  interactions (max.order, forbid, groups, with NA) give exactly additive predictions.
- Argument validation of counts, threads, k, sigest, seed, n.cuts, factors.
