r-verify - independent verification of rfit.md and rgen.md (review-3, pinned 01dee4b4)

Method. Own probes (scratchpad r3-verify-r-p1..p11*.R) against r3-lib (1.0-0 at 01dee4b4) and,
for old behaviour, r3-docs-mainlib (0.9-34 built from main; bart2 there is the formula door).
Code paths read in R/data.R, R/generics.R, R/bart.R, R/bartcore.R, R/rbart.R; checked root TODO,
docs/decisions.md (dec-A85, dec-B78/A06, dec-B130), inst/NEWS.Rd, man/bart.Rd, man/dbarts.Rd,
man/dbartsSampler-class.Rd. None of the 14 findings is filed in TODO or ruled in the ledger.
Result: 14 CONFIRMED (one with a severity qualification); 0 refuted. Grouped by shared fix below.

Summary table
  id                  verdict    sev     vs 0.9-34        surface change for VD?          moves draws?
  rfit-01 = rgen-01   CONFIRMED  BLOCKER same (not regr.) yes (offset() starts to count)   only fits that use offset()
  rfit-02             CONFIRMED  BLOCKER same             no (predict.lm behaviour)       no (test/predict only)
  rfit-03             CONFIRMED  MAJOR   same             small (refusals; name lookup)   no (test channel only)
  rfit-04             CONFIRMED  MAJOR   same; NEWS false yes (new refusals)              no
  rgen-05             CONFIRMED  MINOR   message new      yes if fixed via G2 (widening)  no
  rgen-02             QUALIFIED  MAJOR   new feature      no                              no, if fixed as below
  rgen-03             CONFIRMED  MAJOR   new feature      yes (accepted input forms)      no
  rgen-04             CONFIRMED  MAJOR   new feature      small (new field, NaN cols)     no
  rfit-05..09         CONFIRMED  MINOR   see entries      rfit-09 maybe                   no

=============================================================================================
Group G1 - formula terms are not carried to the offset, test and predict frames
(rfit-01/rgen-01, rfit-02; one mechanism: keep the training terms object and replay it)

rfit-01 / rgen-01  CONFIRMED  BLOCKER  not a regression
  Code: R/data.R dbartsData, formula branch. If 'offset' is missing, modelFrameArgs omits it,
  offsetGivenAsScalar stays NA, and "## offset, when in data frame" reads model.offset() only
  when offsetGivenAsScalar is FALSE. The term is in terms(modelFrame) attr "offset" but never read.
  Own probe p1 (y = sin(3 x1) + 5 log(e) + N(0, .1^2), lo = 5 log(e); rmse(fitted - y), sigma):
    1.0-0  term offset(lo)        2.057 2.129   == no offset at all (2.057 2.129)
           offset = lo argument   0.099 0.109
           term + offset = 0      2.057 2.129   NEW sub-case: a scalar offset also drops the term
           term + rep(0, n)       0.099 0.109   (term counts only beside a vector offset=)
           dbartsData(...)@offset length 0; dbarts(...)$data@offset length 0
    0.9-34 identical pattern (2.063 / 0.101 / 2.063 / 2.063 / 0.101).
    lm: coef(y ~ x1 + offset(lo)) == coef(y ~ x1, offset = lo).
  1.0's own refusal text ("poly(), ns(), log(), and offset() are supported") and the tinytest
  cases ([inst/tinytest/test-bart-formula.R:85](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/tinytest/test-bart-formula.R#L85), [inst/tinytest/test-formula-terms.R:77](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/tinytest/test-formula-terms.R#L77), [inst/tinytest/test-hazard-factors.R:212](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/tinytest/test-hazard-factors.R#L212) with
  offset(log(period)) under family = "hazard") only check that the call returns.
  Fix (surface change - VD decides; base R is unambiguous): in dbartsData's formula branch set
  offset <- as.vector(model.offset(modelFrame)) whenever attr(terms(modelFrame), "offset") is
  non-NULL, summing with a scalar argument too (lm/model.offset semantics: every offset() term plus
  the offset= argument). Then route through the same offsetGivenAsScalar = FALSE path, so every
  family that refuses a flat offset (multinomial: refuseFlatOffsetOnMultinomial) refuses the term
  by the same message, and check the hazard person-period expansion carries it
  ([inst/tinytest/test-hazard-factors.R:212](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/inst/tinytest/test-hazard-factors.R#L212) currently passes because the term is ignored). Test default: offset.test must be
  the term evaluated on 'test' (model.offset of the G1 test frame), not the training vector -
  the current "tracks offset" default would apply training rows to test rows. predict: copy
  predict.lm, which evaluates offset() terms on newdata (model.offset(model.frame(Terms, newdata)));
  that departs from dbarts' current "predict without offset = offset-free surface" for the
  offset= argument, so VD should rule whether the term and the argument behave alike. Alternative
  with no draw movement: refuse an offset() term by name and drop it from the "supported" text.
  Tests: offset(o) fit == offset = o fit bitwise (same seed); offset(o) + offset = c sums;
  term on multinomial refused; yhat.test and predict(newdata) carry the term from test data.
  Draws: move only for fits whose formula has an offset() term (none in the reproducibility
  snapshots; grep finds no offset( there).

rfit-02  CONFIRMED  BLOCKER  not a regression
  Code: R/data.R validateXTest -> replayTerms pastes term.labels into "~ a + b" and re-runs
  model.frame on the new data; the training terms' predvars (poly coefs, scale center/scale,
  ns knots) are never stored, so the basis is refit on the test rows.
  Own probe p2 (y = 4 (c - .4)^2, truth f(.7) = .360; prediction at c = .7):
    term         nd1 = (.7, 0, 1)  nd2 = (.7, .71, .69)  one-row newdata
    poly(c, 2)   0.470             0.062                 error "'degree' must be less than..."
    scale(c)     0.128             0.054                 refused as "missing values in
                                                          'scale(c).1', which carried none in
                                                          training" (misleading; sd of 1 row)
    log(c + 1)   0.372             0.372                 0.372 (stateless, correct)
    lm poly      0.368             0.368
    test = nd2 at fit time reproduces the wrong 0.054. 0.9-34: same (0.749 / 0.063).
  Fix (no surface change; copies predict.lm): store the training terms (terms(modelFrame), which
  model.frame already stamped with predvars via makepredictcall) on data@x as an attribute beside
  "term.labels"; in replayTerms build model.frame(delete.response(trainTerms), newdata,
  na.action = na.pass) - column names stay the term labels, so downstream matching is unchanged.
  Keep the label-based replay as fallback for objects without the attribute. Same frame gives the
  offset() term for G1's test/predict side. Tests: predict at a fixed x is invariant to the other
  rows of newdata for poly/scale/ns (splines is base), and equals the fit-time test channel; a
  one-row newdata works. Draws: none (training untouched).

=============================================================================================
Group G2 - the indicators route keeps no level table (rfit-04, rgen-05)

rfit-04  CONFIRMED  MAJOR (a documented claim is false)  not a regression
  Own probe p4 (g in a/b/c, declared levels a,b,c,z; effect a 0, b +2, c -2; x = .5):
    bart(subset = g != "c", factors = "indicators")  cols x g.b   predict a,b,c,z:
      0.508 2.554 0.508 0.508  (c and z silently predicted as a)
    bartBT on the same rows                          0.548 2.571 0.548 0.548
    indicators, all rows, z declared unused          cols x g.a g.b g.c; z -> 1.341 (all-zero
      pattern, a level that does not exist)
    categorical default                              c 1.798, z 0.850 (own bin; documented in
      dbartsData.Rd, so not part of this finding)
    lm(subset = g != "c")                            refuses "factor g has new levels c, z"
    0.9-34 bart2 subset: 0.539 2.542 0.539 0.539 - same.
  inst/NEWS.Rd: "A factor level that training never saw is an error (0.9-34 silently predicted
  it as another level)"; dec-A85 states the same rule. False on this route.
rgen-05  CONFIRMED  MINOR
  Own probe: bart(factors = "indicators") then newdata g = factor("b") (levels b only) ->
  "... use bart() or dbarts(), which track levels across predict by default" though the caller
  used bart(). factor("b", levels = training levels) predicts fine; lm accepts the subset-level
  factor (xlevels).
  Shared fix: record the training level table (levels observed in the kept training rows, as
  lm's .getXlevels with drop.unused.levels) on the indicators route too, recode a test factor
  over it by label (factor(as.character(x), levels = trainLevels)) before expansion, and refuse by
  name any level with no training rows. That makes subset-level test factors work (a widening,
  lm's behaviour - VD) and makes the NEWS sentence true. If VD prefers no widening: keep the
  refusal and reword it ("give the test factor the training levels, or fit with factors =
  \"categorical\""), and still add the unseen-level refusal for rfit-04. Covers bartBT, which
  always takes this route. Tests: subset-away level and declared-unused level refused by name
  under indicators and in bartBT; subset-level test factor predicts == full-level factor.
  Draws: none.

=============================================================================================
rfit-03  CONFIRMED  MAJOR  not a regression
  Code: R/data.R dbartsData explicit branch: offset.test <- rep_len(offset.test, nrow(x.test));
  getTestOffset resolves a bare symbol in names(data) (training) and never in 'test'.
  Own probe p3 (50 train rows, 4 test rows, test$o = 50, y = x + o):
    offset.test = c(7, 8, 9)  -> yhat.test 7.76 8.10 9.37 7.17 (recycled, no error)
    offset.test = o           -> identical to offset.test = df$o (training o[1:4]); the test
                                 frame's own o = 50 ignored; offset.test = te$o gives ~50
    predict(fit, te, offset = c(7, 8, 9)) refuses on length - the fit-time door is the lax one.
    0.9-34 identical.
  Fix: after getTestOffset, refuse a flat offset.test whose length is neither 1 nor nrow(test)
  (same wording predict uses); resolve a bare name in 'test' first (predict.lm evaluates offset
  in newdata), then the caller's frames; drop the training-'data' lookup, or keep it only when
  the length matches. Tests: both cases above. Draws: none (test channel feeds no likelihood).

rgen-02  QUALIFIED  MAJOR (rgen says BLOCKER)  new in 1.0-0 (nbinom family)
  Code: R/generics.R predict.bartNegbin: means go through convertSamplesForCaller (the caller's
  layout), then negbinPpd(means, object$dispersion) recycles the dispersion in the FIT's stored
  layout (flat chain-major if stored combined; chains x samples chain-fastest if split).
  Own probe p9 (2 chains, 25 draws, dispersions vary 3..8; reference built independently):
    fit cC  predict cC  ppd correctly paired
    TRUE    TRUE        TRUE
    TRUE    FALSE       FALSE (equals the mispaired reference)
    FALSE   TRUE        FALSE (equals the mispaired reference)
    FALSE   FALSE       TRUE
  Qualification: means are correct and each noise draw uses a genuine posterior dispersion, so
  the output is drawn from the product of the (mu, r) marginals rather than the joint; wrong but
  usually second-order, and only when the caller's combineChains differs from the fit's. Keep
  BLOCKER if the orchestrator reads "silently wrong" strictly.
  Fix that moves no draws: build the dispersion array exactly as the means are built - fill
  rArr[, s, chain] <- disp[s, chain] beside the means loop, pass it through the same
  convertSamplesForCaller(rArr, n.chains, combineChains), and call rnbinom over that; the RNG
  order and the two currently-correct cells stay bitwise. (extract.bartNegbin's approach -
  draw in split layout then reshape - is also correct but changes the default path's ppd draws.)
  Test: the 4-cell pairing check above (reference built by replaying set.seed).

rgen-03  CONFIRMED  MAJOR  new in 1.0-0 (categorical columns)
  Code: R/bartcore.R bartcoreSamplerSetTestPredictor column branch new.x.test[, column] <-
  as.double(x.test); bartcoreSamplerSetPredictor's column update is the same.
  Own probe p10 (g lo/mid/hi stored as 0/1/2, effects -3/0/+3; baseline test -2.61 0.42 3.43):
    factor(lo,lo,mid) on column g -> 0.40 0.40 3.41   (each shifted one level up)
    factor(hi,hi,hi)              -> "categorical predictor values must be existing category codes"
    codes 0,0,1                   -> -2.63 -2.63 0.37 (the only working form; undocumented)
    character lo,lo,mid           -> 0.86 x3 with only a coercion warning: NA routed although g
                                     had no NA in training, which the whole-frame path refuses
    setPredictor(factor, "g")     -> same code-range refusal
  man/dbartsSampler-class.Rd says x is "a numeric predictor vector"; no encoding is documented.
  Fix: in both column branches, when the target column's varType is categorical, map a factor or
  character value by label through the training level table (mapFactorColumnsToTrainingLevels on
  a one-column frame), refusing unknown labels and introduced NAs by name; decide (VD) whether a
  numeric value on a categorical column is accepted as a 0-based code or refused, and document it.
  Tests: factor/character column update == whole-frame setTestPredictor result; setPredictor too.
  Draws: none for numeric columns.

rgen-04  CONFIRMED  MAJOR  new in 1.0-0 (ordinal family)
  Own probe p11/p11b (60 rows, weights rep(c(1, 0))): mask is really installed (flipping y at
  the zero-weight rows leaves active-row fitted P(H) bitwise 0.4597), yet the ordinal fit has no
  'active' element and extract(type = "loglik") has 0 NaN columns; the probit twin has 'active'
  and 30 NaN columns. pointwiseLogLikelihood's comment and bart.Rd's weights item say ordinal
  follows probit.
  Fix: in packageOrdinalResults (R/bart.R) copy result$active <- sampler$activeRows as the
  single-forest packager does (R/bart.R ~1631), and NaN those columns in ordinalLogLik
  (R/generics.R) as pointwiseLogLikelihood does. Test: mirror test-binary-weight-mask.R's
  NaN-column check for ordinal. Draws: none.

=============================================================================================
Minors (each CONFIRMED with my own probe; no draw movement for any fix)

rfit-05  keepSampler = TRUE, keepTrees = FALSE predict. Probe p5: 1 chain -> length-3 vector
  (one current-state draw; generics.R calls this the "long-standing keepTrees-free reading", so
  intended); 2 chains -> "internal error: row names do not match the observations" for ev and
  ppd (the ppd refusal sits after the failing naming). 0.9-34 2-chain: 1.97e-313 garbage, so
  1.0 is better. Fix: refuse with the existing "predict requires the fit's saved trees" text when
  n.chains > 1 and keepTrees is FALSE (or always - VD; refusing at 1 chain is a surface change).
rfit-06  rbart_vi: set.seed twice at n.threads = 2 -> not identical; 1 thread -> identical; seed =
  7 at 2 vs 1 threads -> differ. Same code as 0.9-34 (dec-B130 port, deprecated). Default
  n.threads = min(cores, n.chains), so set.seed does not reproduce a default run. Fix: say so
  in rbart.Rd ("Same as in bart" is false for these two); changing seeding would move draws.
rfit-07  bart(cbind(tm, st) ~ a, family = "aft"/"hazard") refused with a message ending "or
  family = \"aft\" / \"hazard\"", which the caller gave; x/y door with the same matrix fits. Fix:
  under an explicit survival family, say "write the response as survival::Surv(time, status)";
  accepting cbind on the formula door is a surface widening (VD).
rfit-08  combineChains = NA -> "missing value where TRUE/FALSE needed"; "yes" -> "invalid
  argument type"; c(TRUE, FALSE) -> "the condition has length > 1"; keepTrees = NA is refused
  by name. 0.9-34 same. Fix: the same TRUE/FALSE validator keepTrees uses.
rfit-09  keepCall = FALSE: print omits the call, summary prints "Call: `NULL`()", update ->
  'could not find function "NULL"', formula -> "invalid formula". 0.9-34 stored the same
  call("NULL") (summary.bart is new). Fix: summary uses print's test; update/formula refuse by
  name. Base mimicry would store call = NULL (update.default then says "need an object with call
  component" as for lm) - a change to the fit object's field, VD.
