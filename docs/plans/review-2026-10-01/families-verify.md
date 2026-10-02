Review 3, families lens - independent verification (tree 01dee4b4, library r3-lib)

Method: every finding re-probed with my own scripts (scratchpad r3-verify-families-v1.R .. v7.R), code
read at the cited symbols, and checked against TODO, docs/decisions.md (incl. dec-B154..B165),
docs/design/negative-binomial.md, docs/plans/sbc-family-tiers.md and .github/workflows/sbc.yaml.
Regression vs 0.9-34: none of the ten - every family here except gaussian/probit is new in 1.0-0, and
the probit one (06) only changes a message.

Verdicts: 01 CONFIRMED (known mechanism, much worse than recorded) BLOCKER; 02 QUALIFIED - BLOCKER
in-sample (hazard), MAJOR out of sample; 03 CONFIRMED MAJOR; 04 QUALIFIED MAJOR (extract only);
05-09 CONFIRMED MINOR; 10 QUALIFIED MINOR. None refuted.

-----------------------------------------------------------------------------

families-01  CONFIRMED, BLOCKER.  Known mechanism, new evidence that it is a defect, not slow mixing.

Is it the adjudicated ridge? Same mechanism: the r-vs-psi level ridge of the logit-p parameterization
(mean = r exp(psi)), named in docs/design/negative-binomial.md section 1 ("a real (if mitigated)
mixing cost. If that bites, the log-mean surface is a documented follow-up") and waived in sbc.yaml
as "MIXING, an identifiability ridge, not a defect". But the record describes slow mixing (ACF > 0.1
past lag 200 at n = 150, crossing into the band at 5x thin). What the review found, and I reproduce,
is a chain that never leaves its cold start at ordinary n: not slow, effectively reducible. The
adjudication rested on avg.mu calibrating; it did not look at what users read off r, ppd and loglik.

Mechanism (code): NBResponse::refreshLatents draws r from NBDispersionPrior::drawIndex, already
COLLAPSED over the Polya-Gamma omegas (weights = L_k + r_k S + log prior, S = sum log(1 - p_i)). So
"marginalize the PG latents" is already done and is not the fix. The problem is conditioning on psi:
at fixed psi, r -> r' rescales every mean by r'/r, so with sum(y) in the thousands the conditional is
a point mass on the current r. The forest can only follow r by a coordinated level shift of every
tree, which tree-local moves never propose. r freezes wherever the forest settled during early
burn-in, which is the cold start r = 8 unless the truth is far below it.

Probes (v1, v2; mu = 8 exp(x1), 50 trees, n.burn 500, n.samples 1000, 2 chains):
  n 200,  r0 3  -> chain 1 8:1000 (0 switches), chain 2 8:1000 (0 switches)
  n 1000, r0 3  -> 8:1000 / 8:1000, 0 switches
  n 1000, r0 30 -> 8:1000 / 8:1000, 0 switches
  n 1000, r0 3: oracle NB log-lik at true mu, relative to r = 3:
    r=1 -168, r=2 -17.2, r=3 0, r=4 -22.3, ..., r=8 -191.3   (data decisively say r = 3)
  engine's own conditional at a stored psi: r=8 1.0, r=6 6.5e-45, r=10 1.3e-93, r<=4 ~0
  log-mean conditional p(r | mu, y) at the SAME stored means: r=3 0.999 (one Gibbs step would move)
  90% ppd coverage of fresh y: estimated r 0.78; r fixed at 3 0.915
  posterior-mean total loglik: estimated r -3557.5; r fixed at 3 -3463.1
So extract(type = "dispersion"), ppd, loglik (and anything built on them: loo/WAIC, plot) silently
report the cold start, not the posterior. The mean channel (ev) is identified and roughly fine.
Gates missed it: negbin-exact.R is single-tree, tiny n (r | psi broad there); the SBC arm uses
n = 150 with a tight psi prior and its flag was waived; test-nbinom.R checks grid membership only.

Fix options (all exact; none needs a new augmentation):
  A. Log-mean parameterization (the design doc's own recorded escape): forest fits eta = log mu - o,
     psi = eta + o - log r. r is then a grid Gibbs step on p(r | mu, y), still collapsed over omega,
     still using the precomputed L_k; per sweep O(13 n) logs, negligible beside the PG draw's
     sum(y + r) loop. The working response re-anchors by log r when r moves (setOffset's reshift).
     Costs: type = "bart"/link on nbinom becomes the log mean (cleaner; the brms/Stan neg_binomial_2
     convention; new in 1.0, no tombstone); leaf prior needs a log-mean calibration (a fitShift such
     as log(mean y) and a node.scale on that scale), a small design item; every nbinom draw moves;
     negbin-exact.R, the SBC nbinom arm and the nbinom docs re-record.
  B. Keep logit-p, add a joint ridge move: propose a grid neighbour r' with a forest level shift of
     total c = log(r / r') split across trees from drawLevelShift's per-tree Gaussian conditional
     (non-zero-sum version), accept by MH (mean is unchanged, so the ratio is the dispersion part of
     the NB likelihood x r prior x the level-sum prior ratio), run before the omega draw. Costs:
     chain.hpp work across the response/forest seam (refreshLatents has no forest access), constant
     leaves only (as drawLevelShift), ~200-300 lines plus component tests; only estimated-r draws
     move; surface unchanged.
  C. For 1.0-0 refuse estimated dispersion (nbinom() requires dispersion =) and document. Cheap,
     loses a feature, a surface change.
  Recommendation: A - it is a smaller, exact change than B, removes the confound instead of
  working around it, and gives users f on the scale they reason in. Maintainer decision: yes (the
  logit-p choice was an agent pick, not a VD ruling - dec-B15 ruled integer r only; A changes what
  the link means and moves draws; the sbc.yaml waiver should be withdrawn either way).
Tests: a mixing gate (n 500-1000, r0 in {2, 30}, default settings, 2 chains: every chain visits the
r0 cell, split-Rhat on r < 1.05, 90% ppd coverage of fresh y within +-0.04); drop SBC_EXPECTED_FLAGS
r,agg.psi; keep negbin-exact.R's stationary-law check (re-derived on the new scale under A).

-----------------------------------------------------------------------------

families-02  QUALIFIED.  Hazard in-sample: BLOCKER. Out of sample (hazard, aft): MAJOR.
Probe (v3; hazard, offset rep(c(-2, 2)), 20 trees, 200 draws):
  S(period 1) by offset group: survivalProbabilities(fit) 0.679 / 0.916;
  1 - stored period-1 hazard 0.870 / 0.705 (direction reversed: the forest learned around o)
  same data, no offset: 0.8009 vs 0.8009 (path correct when there is no offset)
  survivalProbabilities(f, newdata, offset =) -> refused, "takes 'times' and 'newdata' alone"
  aft S(1) by group: training 0.0004 / 1.000; newdata = x 0.900 / 0.913
Cause confirmed: hazardSurvivalProbabilities' training branch calls predict(object, bigX, "ev") with
no offset although the fit carries it (fit$data@offset, per expanded row). The aft training branch
reads the stored channel and is right. Out of sample, dropping an offset= argument is base R's
default (dec-B154 restates it), so that part is not wrong by itself; what is wrong is that there is
no way to supply one, and under dec-B154 a formula offset() term must be evaluated on newdata.
Draws: none move (R-side reporting). Surface: survivalProbabilities gains an offset argument -
additive, on a function new in 1.0; small maintainer item (see shared group G1).
Fix: training branch replays with the period-1 rows' offsets replicated over K; add offset (per
subject, replicated over periods for hazard) passed through to predict / codedRowDraws; evaluate a
formula offset() term per dec-B154. Tests: hazard and aft with offset, survivalProbabilities(fit)
equal to cumprod(1 - stored ev) at at-risk periods; newdata = training x with offset reproduces it.

families-03  CONFIRMED, MAJOR.
Probe (v5; ordinal, offset rep(c(-2, 2)), 2 chains): max |mean predict(f, x, "ev") - mean extract(f,
"ev")| = 0.646; predict(f, x, offset = off) refused ("no out-of-sample offset channel"). The fit
accepts offset and offset.test, so bart.Rd's "neither family has an out-of-sample offset channel" is
false for ordinal, and the pending dec-A114 entry ("those fits have no out-of-sample offset") rests
on the same false premise. dec-B154 (offset() term evaluated on predict's newdata for any family
that takes an offset) now requires the channel. Draws: none. Surface: offset on ordinal predict goes
from refused to accepted; the ledger's A114 text needs correcting - maintainer to confirm, but the
recommendation follows from B154. Fix: predictCodedTest(object$fit, rows$x, offset, ...) exactly as
probit (latent eta + o before the threshold differences). Test: predict at training rows with the
training offset equals extract(ev) to 1e-12.

families-04  QUALIFIED, MAJOR for extract on count-matrix fits; predict defensible.
Probe (v6; every row 10 trials, row 1 zero trials): extract(f, "ppd") int [1:50, 1:40] category
codes; the zero-trial row still draws categories; its loglik is 0. bart.Rd documents "one category
per posterior draw", so this is documented, but for count rows it is not a posterior predictive of
the response the fit modelled. For predict on newdata the trial count is unknown, so a one-trial
draw is a defensible definition if documented.
Surface options: (a) on a count-matrix fit, extract ppd returns rmultinom(n_i) counts, K-widened
(zero vector at n_i = 0); predict keeps one-trial category draws, documented. (b) Refuse ppd on
count-matrix fits with any n_i != 1. (c) Document only. Recommend (a): the shape then matches the
response, as every other family's ppd does. Draws: ppd draws on count fits change. Maintainer: yes
(shape of a documented output). Test: count-matrix ppd row sums equal n_i, zero rows all zero.

families-05  CONFIRMED, MINOR.  extract(..., "ppd") on a student fit: "posterior predictive sampling
does not support student residuals" (sampleFromPPD, deliberate per its comment, undocumented on
bart.Rd/bartBT.Rd). Fix: f + sigma / sqrt(w) * rt(nu) with per-draw nu paired as loglik pairs it;
new draws only where refused today. No surface decision.

families-06  CONFIRMED, MINOR.  All four messages reproduce (hazard all censored -> 'family "probit"
requires a response coded 0/1'; probit all-1 -> same, and it is preceded by a misleading warning
"response values are indistinguishable ... center and/or rescale the response", which the review did
not list; hurdle + variance names "probit"; aft NA status -> "missing value where TRUE/FALSE
needed"; NA time refused instead of routed through na.action). Message-only except the NA routing
(a behaviour change: NA survival rows drop under na.action like NA y). No draws move.

families-07  CONFIRMED, MINOR.  family(bartBT binary) prints probit(sigma = chisq(3, 0.9));
dbartsFamilies$probit(sigma = ...) errors "unused argument"; bart(x, yb, family = family(fb))
accepted. Fix as the review says (drop the residual prior from a binary spec; validate settings
against the token). No draws move.

families-08  CONFIRMED, MINOR.  Both hurdle na.omit refusals reproduce. Fix: apply na.action to
(x, y) before the zero/positive split. Not documented as a limitation anywhere I found.

families-09  CONFIRMED, MINOR.  predict(logistic, type = "ppd", weights = c(1.5, 2, 3)) -> base
"NAs produced" warning and NA draws, while the fit refuses 1.5 by name. Fix: run the fit-time
logistic weight check on predict's weights. No draws move.

families-10  QUALIFIED, MINOR.  Reproduces (nbinom refuses 0/1 weights; ordinal and probit accept
them as the active-row mask, and NBResponse::supportsActiveRows is true). Not filed or ruled. It is a
consistency choice rather than a defect: either accept 0/1 weights as the mask (reusing the probit
policy) or reword the refusal to point at subset/$setActiveRows. Recommend accepting; light
maintainer item since it widens the nbinom surface.

-----------------------------------------------------------------------------

Shared fixes
G1 offset replay (02, 03, plus dec-B154's implementation): one rule - every replay path that a fit
   with an offset can reach takes an offset argument with predict's semantics and evaluates a formula
   offset() term on newdata; in-sample replays use the stored training offset. Touches
   survivalProbabilities (hazard and aft), predict.bartOrdinal, bart.Rd offset item, dec-A114 text.
G2 per-family ppd (04, 05, 09): sampleFromPPD / multinomialPpdFromProbs - count-row multinomial,
   student-t noise, weight validation. One tinytest file covering ppd per family would have caught
   all three; the review's gate-gap note (no gate covers R-side reporting) stands.
G3 refusals and messages (06, 07, 08, 10): R-only, no draws.
01 stands alone and is the one release-relevant engine decision.
