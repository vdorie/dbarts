Review 3, wave 2, lens: families (tree 01dee4b4, library r3-lib)

Covered: every response family end to end through bart()/dbarts() - gaussian, student (df estimated
and fixed), probit, logistic (integer weights), aft, hazard / hazard.logistic (default grid, quantile
and explicit breaks, offset, weights, test set), hurdle.lognormal, nbinom (estimated and fixed r,
offset), ordinal (offset, K = 2, numeric response, unobserved levels), multinomial (factor and count
matrix, zero-trial rows, K = 2 as binomial-with-trials). Per family: extract/predict/fitted types
(ev, response, link/bart, ppd, loglik, scalar channels) checked against hand computations, including
multi-chain pairing under combineChains TRUE/FALSE; survivalProbabilities for aft and hazard;
refusals and family()/dbartsFamilies constructors, base-R family objects; edge data (all censored,
all events, single period, all-zero counts, counts near the 1e6 cap, single observed category,
unobserved middle category, constant positive part, NA in response/status). Exact checks written
here: Student-t nu grid posterior plus E[mu] (fixed sigma, root-only tree, 2-D quadrature) and AFT
(mu, sigma) posterior with sigma DRAWN under 37% censoring (2-D quadrature) - both pass; these are
the two pieces t-exact.R and aft-exact.R leave ungated (both pin sigma, t-exact also pins nu).
Mixing of the per-family scalar parameters (nbinom r, student nu, ordinal thresholds) on realistic n.

Not covered: heteroscedastic aft beyond ppd pairing, monotone/linear/gp leaves per family, BCF and
multi-forest families, xbart families, dbartsSampler mutation methods per family ($setCounts,
$setCategoryOffset, ...), flat C API families, Windows.

Probes: scratchpad r3-families-p1.R .. p30.R (header r3-families-hdr.R: B() fills n.trees = 10,
n.samples = 50, n.burn = 20, n.chains = 1, n.threads = 1; tr() turns errors/warnings into strings).

Findings: 2 BLOCKER, 2 MAJOR, 6 MINOR.

-----------------------------------------------------------------------------

families-01  BLOCKER
Location: src/bartcore/model.hpp NBResponse::refreshLatents / NBDispersionPrior::drawIndex (the r
  step conditions on psi = f + o); R/bart.R bart2Negbin (reports r$dispersion as the posterior).
Claim: the estimated negative-binomial dispersion r essentially never moves after burn-in: given
  psi the full conditional of r is a point mass (mean = r exp(psi), so changing r at fixed psi
  rescales every mean), so each chain freezes near wherever burn-in left it, usually the cold
  start r = 8; the reported dispersion draws, ppd and loglik are not posterior draws.
Probe (p11, p13, p14, p22; n = 200, mu = 8 exp(x1), default bart() settings):
  true r 50: four chains, n.burn 1000, n.samples 2000 each:
    chain 1 : 10:394 12:1606 ; switches: 2      chain 2 : 12:1528 15:472 ; switches: 1
    chain 3 : 12:2000        ; switches: 0      chain 4 : 10:151 12:1849 ; switches: 1
  true r 1: every chain 2:2000, switches 0.
  n.chains = 2 defaults, r draws by truth:  r=1 -> all 2;  r=5 -> all 8;  r=50 -> 8/10/12
    90% ppd coverage of fresh y: 0.85, 0.90, 0.97 (nominal 0.90)
  grid full conditional of r at a stored draw's psi (my own computation of drawIndex's weights):
    r=1 5.4e-23, r=2 1.0, r=3 3.7e-23, r=4 3.5e-69, ...
  yet the r = 1 fit of the same data has essentially the same log-likelihood
    (sum over obs, posterior mean): estimated-r fit -708.9, nbinom(dispersion = 1) fit -709.4,
    so the marginal posterior has real mass at r = 1 that the chain can never reach.
  The sampler's own state agrees: s$getDispersion() after each of 200 single-sweep runs: 8 (200 of 200).
Why gates missed: negbin-exact.R is a single tree on a two-cell design with tiny n, where p(r | psi)
  is broad; it checks the stationary law (correct), not mixing. test-nbinom.R checks that draws lie
  on the grid. docs/design/negative-binomial.md measured r mixing with psi held fixed and flagged
  the f/log r confounding as making ESS "optimistic"; it did not foresee a frozen chain.
Fix: move r along the mean-preserving ridge - propose r' with a compensating level shift
  log(r / r') of the forest (the level-fibre move drawLevelShift is the existing template) and
  accept by MH on the collapsed NB likelihood; until then document that r does not mix, or default
  to a fixed r.

families-02  BLOCKER
Location: R/bart.R hazardSurvivalProbabilities (training branch: predict(object, bigX, type = "ev")
  with no offset; newdata branch: codedRowDraws(..., NULL offset)); survivalProbabilities.bart (aft
  newdata branch: predict(..., type = "bart") with no offset); man/survivalProbabilities.Rd details.
Claim: on a hazard fit with an offset, survivalProbabilities(fit) silently drops the training offset
  (the documented h(k | x) = g(f(x, k) + o) loses o), and neither family can carry an offset to
  newdata (survivalProbabilities refuses 'offset' by name), so curves for an offset fit are wrong
  in-sample (hazard) and out of sample (hazard and aft).
Probe (p4, p5, p15):
  hazard, offset = rep(c(-1, 1)), subject 1 (offset -1, t = 2):
                 [,1]      [,2]
    manual   0.9274745 0.7945967     # cumprod(1 - h) from extract(fit, "ev") on its rows
    survProb 0.7044648 0.4074754     # survivalProbabilities(fit)[, 1:2, 1], column means
    hazards: stored 0.0725 0.1479 | predict, no offset 0.2955 0.4439 | predict, offset 0.0725 0.1479
  hazard with test = x[1:4, ], offset.test = -1: S at period 4, same four subjects
    storedTest 0.0608 0.1652 0.1152 0.0770 | newdata = x[1:4, ] 0.0017 0.0113 0.0064 0.0025
  aft, offset = rep(c(-2, 2)): S(1 | x) mean by offset group
    train   -2: 2.4e-06   2: 1.0000      newdata = x (same rows)  -2: 0.934   2: 0.901
Why gates missed: no exact gate, tinytest or reduction gate puts an offset on aft or hazard
  (hazard-exact.R, hazard-reduction.R, aft-exact.R: zero mentions of offset); the aft training path
  is right only because it reads the stored channel, which already carries the offset.
Fix: give survivalProbabilities an offset argument (per subject, replicated over periods for
  hazard), and in the hazard training branch replay with the fit's own per-subject offset (the
  period-1 rows' offsets).

families-03  MAJOR
Location: R/generics.R predict.bartOrdinal (refusePredictOffsetChannel); R/bart.R bart2Ordinal /
  dbarts(family = "ordinal") offset acceptance; man/bart.Rd \item{offset} ("neither family has an
  out-of-sample offset channel").
Claim: an ordinal fit accepts offset and offset.test and uses them (the stored train and test
  probabilities include them), but predict refuses an offset and replays the offset-free surface,
  so predict at the training rows silently disagrees with fitted().
Probe (p17): offset = rep(c(-3, 3)), y from cut(off + 2 x1 + e):
  fit$latent.train means by offset: -3: -1.99   3: 4.12      (offset is in the fit)
  extract(test) vs extract(train) on the same 4 rows with offset.test: max diff 7.8e-16
  predict(f, x, type = "ev", offset = off)
    ERROR: 'offset' is not used by predict on a bartOrdinal fit: this fit has no out-of-sample
    offset channel; predict replays the offset-free surface
  max |predict(f, x, type = "ev") - extract(f, "ev")| = 0.866
Why gates missed: ordinal-exact.R and test-ordinal.R fit without an offset; the predict refusal is
  tested as a refusal (dec on the shared predict signature), never against an offset fit.
Fix: add the offset channel to predict.bartOrdinal (the latent is f + o, as probit's), or refuse
  offset/offset.test at fit time for ordinal so the documented "no offset channel" is true.

families-04  MAJOR
Location: R/generics.R multinomialPpdFromProbs (used by extract.bartMultinomial and
  predict.bartMultinomial type = "ppd").
Claim: on a count-matrix multinomial fit (rows with n_i trials) the posterior predictive draws one
  category per draw instead of a count vector of n_i trials, and a zero-trial row still gets a
  category; this contradicts the fit's own likelihood unit (loglik is the whole count row,
  residuals are proportions y / n_i).
Probe (p16): every row 10 trials, K = 3:
  head(f$y): 2 4 4 / 2 2 6 / ...      str(extract(f, type = "ppd")): int [1:50, 1:40] 2 3 2 2 3 ...
  predict(f, x[1:2, ], type = "ppd"): int [1:50, 1:2] category codes
Why gates missed: multinomial-exact.R gates probabilities only; tinytest checks ppd shapes; bart.Rd
  documents the single-category draw without distinguishing count rows.
Fix: for a count response draw rmultinom(1, n_i, p) per (draw, row) and return a K-widened count
  array (zero vector at a zero-trial row); keep category codes for the single-trial factor case.

families-05  MINOR
Location: R/generics.R sampleFromPPD ("posterior predictive sampling does not support student
  residuals").
Claim: extract/predict type = "ppd" on a student() fit errors, though the fit stores the per-draw
  scale and nu (resid.df) that the draw needs; no help page says ppd is unavailable for student.
Probe (p17): extract(B(x, y, family = "student", weights = w), "ppd")
  ERROR: posterior predictive sampling does not support student residuals
Why gates missed: test-pointwise-loglik.R asserts the refusal.
Fix: draw f + sigma / sqrt(w) * rt(1, nu) per draw with nu from resid.df (paired as loglik already
  pairs it), or document the refusal on bart.Rd/bartBT.Rd type = "ppd".

families-06  MINOR
Location: R/dbarts.R binary response validation; R/spec.R / bart variance check; hazard ingestion;
  R/dbarts.R parseSurvivalResponse.
Claim: several refusals name the wrong thing - the internal component family, or a base error.
Probe (p6, p26, p27):
  hazard, every subject censored:  ERROR: family "probit" requires a response coded 0/1
  probit, y all 1:                 ERROR: family "probit" requires a response coded 0/1
                                   (the response is coded 0/1; the problem is one class)
  hurdle.lognormal + variance = ~a: ERROR: a variance forest requires family = "gaussian" or
                                   "aft"; family "probit" routes precision through its own latent channel
  hazard or aft, NA in status:     ERROR: missing value where TRUE/FALSE needed
  (and an NA survival time is refused, "survival times must be finite and positive", where every
  other family drops an NA response under na.action - e.g. ordinal or nbinom with one NA y fit 39 rows)
Why gates missed: message text for these paths is untested.
Fix: say "every subject is censored" / "the response has a single class" / "hurdle.lognormal does
  not take variance"; anyNA-check status; route NA survival rows through na.action like NA y.

families-07  MINOR
Location: R/bart.R bartBT family.spec construction; R/family.R specifiedFamily.
Claim: family() of a binary bartBT fit is probit(sigma = chisq(3, 0.9)), a call probit() cannot
  build ("Printing one names the call that would build it", dbartsFamilies.Rd); passing it back to
  bart() is accepted silently.
Probe (p20, p27):
  family(bartBT(x, yb, ...))                              probit(sigma = chisq(3, 0.9))
  dbartsFamilies$probit(sigma = dbartsPriors$chisq(3, .9)) ERROR: unused argument (sigma = ...)
  bart(x, yb, family = family(fb), ...)                   accepted; family(): probit(sigma = chisq(3, 0.9))
Why gates missed: family() is tested on bart fits, whose binary spec carries no sigma.
Fix: drop the residual-prior setting from a binary family's spec in bartBT (and validate settings
  against the token in newValidated("dbartsFamily")).

families-08  MINOR
Location: R/bart.R bart2Hurdle (refuseHurdlePositiveMissingness runs on the raw x; the positive
  fit's test = full x).
Claim: hurdle.lognormal ignores na.action for missing predictors: with na.action = na.omit a column
  NA only on zero rows is still refused, and NA on both row sets fails with a message about 'test'
  predictors the user never supplied.
Probe (p3, p30):
  NA in 'a' on two zero rows, na.action = na.omit: ERROR: ... 'a' carry missing values only on the
    zero (y == 0) rows ...
  NA on one positive and one zero row, na.omit: ERROR: test predictors have missing values in 'a',
    which carried none in training ...
Why gates missed: test-hurdle*.R exercise the default na.keepPredictors only.
Fix: apply na.action to (x, y) before the split, or refuse na.action other than the default for
  hurdle by name.

families-09  MINOR
Location: R/generics.R sampleFromPPD (weighted binary branch).
Claim: predict(type = "ppd", weights = ) on a logistic fit accepts non-integer weights and returns
  NA draws with base R's "NAs produced" warning, where the fit itself refuses non-integer weights.
Probe (p17): predict(f, x[1:3, ], type = "ppd", weights = c(1.5, 2, 3))
  WARN: NAs produced;  [1,] NA 2 0 / [2,] NA 1 0
Why gates missed: weighted ppd is tested with integer weights only.
Fix: validate predict weights with the fit-time logistic policy (positive integers).

families-10  MINOR
Location: R/dbarts.R nbinom weight refusal.
Claim: nbinom refuses 0/1 weights ("exposure belongs in the offset as a log-exposure term") though
  NBResponse supports the active-row mask and probit/ordinal accept 0/1 weights as exactly that mask;
  the message answers a question the caller did not ask.
Probe (p26): B(x, rpois(n, 3), family = "nbinom", weights = rep(c(0, 1), n / 2))
  ERROR: nbinom (count) models do not support weights: exposure belongs in the offset ...
Why gates missed: test-nbinom.R asserts the refusal for real-valued weights only.
Fix: accept all-0/1 weights as the mask (as probit/ordinal), or say why the mask is not offered.

-----------------------------------------------------------------------------

Checked and found correct
- Student-t: nu grid posterior and E[mu] match a 2-D quadrature (fixed sigma, root tree, 80k draws,
  all nine grid z-scores within +-0.4); nu mixes (110-120 switches per 1000 draws at nu = 20);
  loglik (marginal t at sigma / sqrt(w), per-draw resid.df) exact under combineChains TRUE/FALSE.
- AFT: (mu, sigma) posterior with sigma drawn under censoring matches quadrature (E[mu] 1.1523 vs
  1.1524, E[sigma] 0.7962 vs 0.7959, P(sigma < 0.6) 0.0352 vs 0.0356); loglik uses the upper tail
  for censored rows; fixed(v) pins sigma = sqrt(v); all-censored and all-event fits run.
- ppd noise pairs with each draw's sigma (gaussian and aft, 3 chains, extract and predict,
  combineChains TRUE/FALSE: cor(var(noise), sigma^2) 0.71-0.80 vs ~0 if mispaired).
- Ordinal: yhat = Phi-differences of latent and per-draw thresholds, loglik, predict replay
  (2e-15) under 3 chains both layouts; thresholds mix (lag-50 autocorrelation ~0); K = 2, numeric
  and unobserved-level responses fit.
- nbinom: loglik equals dnbinom(y, r_draw, mu_draw) (combined and split); mean = r exp(psi);
  count cap 1e6 refused by name; fixed r honored (the r mixing defect aside).
- Hurdle: ev = pi exp(f + sigma^2 / 2), loglik with the -log y Jacobian, predict at training rows
  equals the in-sample channel (2 chains); weights/offset/subset/test refused by name.
- Probit/logistic: ev = link(latent) with offset, predict with offset reproduces ev, loglik
  (weighted logistic = w * Bernoulli); multinomial K = 2 count matrix fits binomial-with-trials;
  zero-trial rows give loglik 0 and fitted probabilities with the documented warning.
- Hazard: expansion, default/quantile/explicit grids, breaks validation, horizons before the
  first period give S = 1, data-frame and matrix newdata agree, stored-test path without offset.
- family(): every bart family returns its dbartsFamily; stats binomial(probit/logit) map;
  poisson, cloglog, quasibinomial refused by name; constructor argument validation.

Exact-gate coverage gaps (what the shipped gates cannot see)
- No exact gate puts an offset on aft, hazard, ordinal, nbinom or student; only logistic-reference
  and multinomial arm 6 do (families-02, -03 live there).
- t-exact pins both nu and sigma, aft-exact pins sigma: the residual-scale draw under the t mixture
  and under censoring is ungated (probed here and correct for nu and for aft sigma).
- negbin-exact is single-tree, tiny n: it validates the stationary law but cannot see r mixing.
- No gate covers R-side reporting: ppd per family, survivalProbabilities, count-row ppd.
