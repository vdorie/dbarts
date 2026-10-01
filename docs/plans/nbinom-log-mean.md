# nbinom-log-mean: the forest models log(mean), and r is drawn given the mean

Status: PLANNED (dec-B170 in [decisions.md](../decisions.md), 2026-10-01). Blocks 1.0-0: the third whole-branch
review's finding families-01 is a blocker.

agent: opus for C1-C2 (engine numerics, the r step, the anchor, the exact-gate re-derivation); sonnet for C3-C5
(SBC arm, equivalence re-record, records). Serialized: one implementer, each commit gated before the next.
rng: POSTERIOR-CHANGING for family = "nbinom" only - every nbinom draw moves, fixed r included (see
[Fixed r](#fixed-r)). NEUTRAL for every other family: no other ResponseModel, kernel or RNG call order is touched,
so the equivalence compare must show the 54 non-nbinom scenarios identical and bcf/multinomial harnesses bitwise.
window: pre-release, before the 1.0-0 merge.
budget: ~920 lines (engine ~150, bridge ~15, R ~60 with Q1 (d), Rd and dbarts.h comments ~40, tests/cpp ~200,
tinytest ~90, negbin-exact.R ~70 changed, negbin-mixing.R ~170 new, sbc.R and sbc.yaml ~50, design doc ~120, records
~40). Plans have run 1.5-2x low: expect up to ~1800. Compute: the SBC nbinom arm (46 min at R=200 today) plus its
burn ladder, the exact gates in FULL mode, the equivalence trio; bench-sampler.R compare is maintainer-run.

## Goal

A negative-binomial fit's forest models eta(x) = log E[y | x] less the offset, as MASS::glm.nb and Stan's
neg_binomial_2 do. The dispersion r is a grid Gibbs step on p(r | mu, y), collapsed over the Polya-Gamma latents,
and it mixes: every chain of the mixing gate leaves the cold start, the chains agree and the pooled posterior covers
the true r, predictive intervals cover at nominal rate, and the SBC arm passes with its waiver withdrawn. type =
"link" is the log mean; the leaf prior is centered at the data's log rate and its spread recalibrated for the
log-mean scale.

## Context

### The defect

Under the shipped logit-p parameterization ([`NBResponse`](../../src/bartcore/model.hpp)) the forest fits
psi = f + o, p = plogis(psi), mu = r exp(psi), and [`NBResponse::refreshLatents`](../../src/bartcore/model.hpp)
draws r from p(r | psi, y). The review's probes (n 200 and 1000, mu = 8 exp(x1), defaults): every chain froze
within burn-in, mostly at the cold start r = 8 (true r 3, 5 and 30 all read 8; true r 1 read 2; true r 50 read 8
to 15 by chain, one or two switches in 2000 draws). At n = 1000, r0 = 3: the engine's own conditional at a stored
psi put mass 1.0 on r = 8 (r = 10 at 1.3e-93); 90% ppd coverage of fresh y was 0.78 (0.915 with r fixed at 3);
posterior-mean total log-likelihood -3557.5 against -3463.1 at r = 3. extract(type = "dispersion"), ppd, loglik,
and anything built on them report the cold start. The mean channel is roughly fine.

Mechanism. At fixed psi, r -> r' multiplies every mean by r'/r. The NB2 information on a common log-mean shift is
sum_i r mu_i / (mu_i + r) (about 5000 at the probe's n = 1000 and r = 8), so the grid neighbour 8 -> 10 (a shift
of log 1.25) costs of order 10^2 nats. The forest could follow r only by a coordinated level shift of every tree, which
tree-local moves never propose. negbin-exact.R (one tree, n = 50) checks the stationary law, not mixing; the SBC
arm's r flag was waived as a ridge (["SBC_EXPECTED_FLAGS: r,agg.psi"](../../.github/workflows/sbc.yaml)).

### The new conditionals

Write eta_i = f(x_i) + c + o_i = log mu_i, with c the prior center ([Leaf prior](#leaf-prior)), and
psi_i = eta_i - log r, so p_i = mu_i / (mu_i + r). The count law is unchanged:

    NB(y | r, p) = G(y + r) / (G(r) y!) p^y (1 - p)^r,    p^y (1 - p)^r = e^{y psi} / (1 + e^psi)^{y + r}.

PG step (Polson-Scott-Windle). e^{y psi} / (1 + e^psi)^b = 2^-b e^{kappa psi} E[e^{-omega psi^2 / 2}] with
omega ~ PG(b, 0), b = y + r, kappa = (y - r)/2. So omega_i | r, f ~ PG(y_i + r, psi_i), and given (omega, r) the
f-likelihood is Gaussian, exp(-omega_i/2 (psi_i - kappa_i/omega_i)^2). Since psi_i = f_i + a_i with the anchor
a_i = o_i + c - log r, the trees see working response z_i = kappa_i/omega_i - a_i at precision omega_i: the logistic
seam unchanged, with -log r and c entering exactly as an offset does.

r step. Integrating omega out of p(y, omega | f, r) returns the NB likelihood, so

    p(r | f, y) ~ pi(r) prod_i NB(y_i | r, mu_i)
    log w_k = sum_i [lgamma(y_i + r_k) - lgamma(r_k)] + sum_i [y_i log mu_i + r_k log r_k - (y_i + r_k) log(mu_i + r_k)] + log pi_k
            = K_k - sum_i (y_i + r_k) log(1 + e^{eta_i - log r_k}) + log pi_k + (r-free),
    K_k = L_k - Y log r_k,   L_k = sum_c n_c [lgamma(c + r_k) - lgamma(r_k)],   Y = sum_i y_i,

using (y + r) log(mu + r) = (y + r) log r + (y + r) log(1 + mu/r) and n r_k log r_k - sum_i (y_i + r_k) log r_k =
-Y log r_k. K_k depends on y alone and precomputes where L_k does
([`NBDispersionPrior::computeKernel`](../../src/bartcore/model.hpp)). The rest no longer separates: under logit-p
p_i was r-free and the eta-dependence collapsed to one statistic S; now p_i moves with r_k, so the per-sweep cost is
13 n evaluations of [`logOnePlusExp`](../../src/bartcore/model.hpp) (one per row and grid point), beside the PG
draw's sum_i (y_i + r) unit draws.

Exactness. The sweep is a two-block Gibbs sampler: block 1 draws the trees given (omega, r, y) (the Gaussian
working-response conditional); block 2 draws (r, omega) jointly given the trees, as p(r | f, y) p(omega | r, f, y):
r from its omega-marginal, then omega ~ PG(y + r, psi(r)) at the new r. A joint draw of a block is a valid Gibbs
step. The order is the one the code already uses (r first, never reading omega; omega regenerated before anything
conditions on it), so the ordering argument of
[5. State and mutation](../design/negative-binomial.md#5-state-and-mutation) carries over; what changes is what
is held fixed while r moves: mu, not psi. Probe (not checked in): a cell-means stand-in (fixed partition, constant leaves, these
two blocks, omega from [`dbartsDrawLatents`](../../R/augmentation.R)), n = 60 over 3 cells, 58000 kept draws: the
grid posterior of r matches the quadrature reference within 0.0015 at every grid point.

Why it mixes. In the NB2 (mu, r) parameterization the Fisher information is diagonal: d2l / dmu dr =
(y - mu) / (mu + r)^2, mean zero. Conditioning r on mu is conditioning on an orthogonal parameter, so a move of r
costs no likelihood through the mean. Under logit-p the cross term is the full log-mean information divided by r,
which is the ridge. The same stand-in at n = 1000, 10 cells, mu = 8 exp(-1..1), start r = 8, 1200 kept sweeps:

    true r   logit-p (shipped)                      log-mean (this plan)
    2        visits {4,5,6,8}, never 2, 3 switches  P(r = 2) 1.000 (exact 1.000)
    30       frozen at 8, 0 switches                574 switches, P(r = 30) 0.618 (exact 0.603), lag-1 ACF -0.01

### Fixed r

With r fixed, the PG draw and the tree block are the logit-p ones with anchor o + c - log r in place of o: the
fixed-r log-mean model IS the logit-p model with offset o + c - log r. The likelihood family is the same; the
reported link moves by log r (eta = psi + log r); the posterior moves only through the leaf prior. Logit-p centered
the prior on log mu at log r + o (mean count r whatever the data); this plan centers it at c + o and recalibrates its
spread (Q1). The two coincide only when c = log r and the leaf prior is kept, so under the defaults fixed-r draws move
as well, and fixed-r fits get a data-centered prior they lacked.

## Design

### Engine

[`NBResponse`](../../src/bartcore/model.hpp) gains two scalars, the shift c and a cached log r, and one anchor
a_i = o_i + c - log r used by every working-response build:

- c at construction and in [`NBResponse::setData`](../../src/bartcore/model.hpp):
  c = log(max(sum_i y_i, 1/2) / sum_i e^{o_i}), sum_i e^{o_i} = n with no offset. This is the intercept MLE of the
  Poisson model with offset (and of the NB model when o is constant). The 1/2 floor keeps an all-zero response
  finite. Computed over all rows: the response transform is the full-data one by design, the aft precedent.
- `fitShift()` returns c; fitScale and sigmaScale stay 1. The train and test channels, the fits without offset, the
  predict replay and `forestCalibration().priorMean` then carry c with no chain edit
  ([`Chain::storeSample`](../../src/bartcore/chain.hpp), [`Chain::fitsWithoutOffset`](../../src/bartcore/chain.hpp),
  [`Chain::forestCalibration`](../../src/bartcore/chain.hpp)). The internal fits keep their zero-centered prior,
  as gaussian's do.
- [`NBResponse::refreshLatents`](../../src/bartcore/model.hpp): (1) when estimating, eta_i = totalFits_i + c + o_i
  and r from [`NBDispersionPrior::drawIndex`](../../src/bartcore/model.hpp), now taking (y, eta, active) and
  forming w_k above; update log r; (2)-(3) [`NBResponse::drawOmega`](../../src/bartcore/model.hpp) with
  psi_i = totalFits_i + a_i and z_i = kappa_i/omega_i - a_i. `collapsedStatistic` is deleted.
  [`NBDispersionPrior::computeKernel`](../../src/bartcore/model.hpp) folds -Y log r_k (Y over the active rows) into
  the kernel; log r_k is a constant table. The weight loop is scalar and fixed-order (rows outer, grid inner), so the
  reference and shipped builds stay bitwise.
- [`NBResponse::coldStart`](../../src/bartcore/model.hpp), [`NBResponse::restoreLatents`](../../src/bartcore/model.hpp)
  and [`NBResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) use psi_i = totalFits_i + a_i.
- Re-anchoring. Whenever c, r or o changes without an omega draw, the working response is REBUILT from omega as
  kappa_i/omega_i - a_i, never shifted by a delta: a shift is not bitwise the expression creation evaluates, and the
  setOffset parity test needs bitwise. Sites: [`NBResponse::setOffset`](../../src/bartcore/model.hpp) (replacing
  its [`reshiftWorkingForOffset`](../../src/bartcore/model.hpp) call), restoreScale. Inside refreshLatents r moves
  only before the omega draw, which rebuilds anyway.
- updateScale. The ResponseModel contract ("updateScale re-anchors the transform to the new response, as setOffset
  does; false locks it") now has a transform to act on: setResponse and setOffset with updateScale = TRUE
  recompute c from all rows, then rebuild; FALSE keeps c, the embedded-Gibbs default. setData always recomputes, as
  gaussian's does. setResponse keeps r and rebuilds the kernel, as today.
- State. c rides the existing fit.scale block as (c, c + 1); restoreScale decodes c = min, an exact round trip (a
  midrange encoding such as (c - 1/2, c + 1/2) does not round-trip for about 1% of values). fitScale stays the
  constant 1 and never reads the width. getScale writes the pair; restoreScale sets c and rebuilds, which
  [`Chain::installForest`](../../src/bartcore/chain.hpp) needs (it restores the scale but not the latents).
  [`Chain::stateIsValid`](../../src/bartcore/chain.hpp) refuses an nbinom state with fitMax <= fitMin, which is
  every pre-change nbinom state ((0, 0)) and the case setState would otherwise skip silently. No new block and no
  version bump: no state format has shipped ([`stateFormatVersion`](../../src/R_interface_bartcore.cpp)). Copies go
  through the same state struct. The restore contract (dispersion before latents) is unchanged; restoreDispersion
  also updates log r.
- Threading. NBResponse is per chain and c is chain-invariant; nothing is shared across chains, and the r step runs
  in the chain's own thread like the PG loop.

### Bridge and C API

- [`drawAugmentationLaws`](../../src/R_interface_bartcore.cpp) nbinom arm: psi = fit + offset - log(dispersion).
  [`computeWorkingResponse`](../../src/R_interface_bartcore.cpp) nbinom arm: kappa/omega + log(dispersion) - offset.
  `fit` there is what getFitsWithoutOffset reports, which now includes c, so the helpers need no shift argument.
- dbarts.h: the "(gaussian only)" in the updateScale text of
  [`dbarts_sampler_setResponse`](../../inst/include/dbarts/dbarts.h) and
  [`dbarts_sampler_setOffset`](../../inst/include/dbarts/dbarts.h) becomes "(gaussian and nbinom)". Comments only:
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) hashes declarations, so it does not move. The train
  channel's meaning for nbinom changes (log mean); no consumer reads it (below).

### Leaf prior

Center: c. The data-centered families (gaussian, aft, the hurdle positive part) center at the response's midrange
on their own scale; the latent binary families center at 0 because their link has a natural zero (p = 1/2). A log
mean has no natural zero, so the gaussian/aft precedent applies on the link scale, using the log rate rather than a
midrange because counts include zeros. dbarts has no other log-link family (hazard is a binary link on person-periods;
aft and hurdle.lognormal are gaussian on a log response), so there is no in-package spread to copy.

Spread: the sigma-free families fix an anchor A, with prior sd of f(x) = A/k. probit and logistic default to
k ~ chi(1.5, 2) (the binary hyperprior, [`isBinaryFamily`](../../R/spec.R) covers only those two); ordinal and
nbinom default to a fixed k = 2, so nbinom today has sd pi sqrt(3)/2 = 2.72. The anchor in
[`defaultLeafScale`](../../R/model.R) (and its C backstop) and, under Q1's recommendation, the k default in
[`resolveLeafHyperprior`](../../R/model.R) move; the values are Q1.

Probe for Q1 (not checked in): a scratch build of this design (the engine change above, anchor settable), r
estimated, one chain, 500 + 500 sweeps, default trees, 500 test rows, 8 replicates at each of n = 300 and n = 2000,
six designs (log-mean signal sd 0.06 to 2, one a step of 3 on the log scale; true r 2 to 30). Fresh-y coverage is
by randomized PIT. Means over designs and replicates:

    option                      n     RMSE eta  cov 90% mu  worst design  cov 90% y  r0 in 95%  best RMSE
    (a) A = pi sqrt(3), k = 2   300   0.354     0.943       0.893         0.886      0.62       0 of 6
                                2000  0.191     0.904       0.807         0.895      0.92       0 of 6
    (b) A = 3, k = 2            300   0.289     0.954       0.887         0.896      0.90       1 of 6
                                2000  0.175     0.921       0.828         0.897      0.92       0 of 6
    (c) A = 2, k = 2            300   0.282     0.927       0.781         0.897      0.96       2 of 6
                                2000  0.163     0.928       0.813         0.896      0.94       3 of 6
    (d) A = 3, k ~ chi(1.5, 2)  300   0.251     0.957       0.898         0.892      0.92       3 of 6
                                2000  0.152     0.934       0.826         0.896      0.94       3 of 6

The worst design is always the strongest signal (log mu = 1 + 2 z). (d) has the lowest RMSE at both n, beating (b)
in 78% and (c) in 73% of paired fits, with coverage of mu at or near the best; (c) under-covers the strong signal
at n = 300 (0.78); (a) is worst throughout and, at n = 300, its wide prior lets the forest absorb over-dispersion
(r0 inside the 95% set in 62% of fits). A named sd (prior.scale) is stated on the log-mean scale with no
conversion (fitScale = 1).

### R surfaces afterwards

| surface | value |
|---|---|
| link / "bart" | eta + o = log mean, per draw (latent.train / latent.test) |
| response / "ev" | exp(link) per draw; r no longer enters ([`negbinMeanCounts`](../../R/bart.R) drops it) |
| ppd | rnbinom(size = r_s, mu = mu_s), unchanged in form, now with r mixing |
| loglik | dnbinom(y, size = r_s, mu = mu_s), unchanged ([`negbinLogLik`](../../R/generics.R)) |
| dispersion | the r draws, unchanged |
| predict | ev = exp(replayed link + offset); ppd pairs each draw with its r ([`predict.bartNegbin`](../../R/generics.R)) |
| fitted(type = "link") | posterior-mean log mean |
| print | "negative binomial (log link)" ([`print.bartNegbin`](../../R/generics.R)) |
| plot, summary | unchanged ([`plot.bartNegbin`](../../R/plot.R)) |
| sampler | run()$train, getFitsWithoutOffset, predict: log mean; getLeafPrior()$prior.mean and response.shift: c |
| augmentation helpers | psi = fit + offset - log(dispersion) ([`dbartsWorkingResponse`](../../R/augmentation.R)) |

[`bart2Negbin`](../../R/bart.R) keeps its per-sample loop (its stated reasons are gone, but collapsing it is a
separate draw-neutral change, not this plan's); the loop body computes exp(train). dispersion.raw stays as the
keepTrees marker predict checks.

r prior: unchanged - the 13-point grid, gamma(2, 0.1) weights renormalized, cold start at the median 8. The ruling
does not reopen it, and both gates below exercise it.

### Gates

- Mixing gate, new: benchmarks/R/negbin-mixing.R. mu = 8 exp(x1), 5 uniform predictors, default bart() settings
  except n.chains = 2, fixed seeds, on cells where r is identified: r0 = 5 at n = 2000 and r0 = 2 at n = 500
  (the scratch build: 5:488 6:12 and 5:469 6:31 per chain; 2:500 on both). Not r0 = 30: at these means r = 30 and
  50 are barely distinguishable, the forest absorbs the variance difference, and a correct sampler puts most mass
  on 50 depending on the seed (the scratch build: 98%; an independent emulation 74-100%). Pass when, per cell:
  (i) each chain left the cold start, under half its draws at r = 8; (ii) the chains agree: split-Rhat on r below
  1.05, where a split half with zero variance is handled by rule - every half constant at one shared value passes
  (Rhat taken as 1), any other zero-variance case fails; (iii) the pooled central 95% set of r contains r0;
  (iv) the 90% ppd coverage of 1000 fresh y by randomized PIT, u = F(y - 1) + V (F(y) - F(y - 1)) with F the
  draws' ecdf and V uniform, lies in 0.90 +- 0.04 (the probe's randomized coverage sat at 0.886 to 0.897). quick
  and full differ in seeds only. Listed in exact-gates.yaml. Discrimination: it must FAIL at the parent commit
  (r frozen at 8), the [Gate hygiene](README.md#gate-hygiene) rule. Estimated at about a minute from the probe's
  fit times.
- [negbin-exact.R](../../benchmarks/R/negbin-exact.R), re-derived: the leaf m is the cell log mean with prior
  N(c, tau^2), tau = A/(k sqrt(numTrees)); the cell log-likelihood in (m, r) is sum lgamma(y + r) - n lgamma(r) -
  sum lgamma(y + 1) + sum y (m + o) + n r log r - sum (y + r) log(r + e^{m + o}); the gated mean is exp(m) per cell
  (train channel less the offset); the quadrature grid centers on c. The estimated arm gains a two-level exposure
  offset (o in {0, log 2} within each cell), so the anchor's o, c and log r terms are all exercised; the fixed arm
  (r = 5) stays offset-free. The script reads c from getLeafPrior()$response.shift and fails unless it equals the
  formula to 1e-12. Tolerances unchanged.
- SBC arm ([`sbcFamilyConfig`](../../benchmarks/R/sbc.R)): draw eta0 from the prior via the generator, y0 ~
  NB(r0, mu = exp(eta0)); functionals r, avg.mu = mean(exp(eta0)), agg.eta = mean(eta0 at the test rows). c is
  data-derived, so the fit must share the generator's: build the fit sampler from the generator's placeholder
  response (same c) and install y0 with setResponse(updateScale = FALSE), rebuilding per replication as today so r
  restarts. Keep the prior sd of f at the current 0.68 (k = A/0.68) and a placeholder mean count near 5, which keeps
  the PG cost inside the budget k = 8 was sized for. Re-measure the burn ladder with
  [`sbcBurnLadder`](../../benchmarks/R/sbc.R) (24000 sweeps was set by the ridge;
  [`sbcBurnSweeps`](../../benchmarks/R/sbc.R)), then thin and timeout. Drop SBC_EXPECTED_FLAGS from the arm and
  its waiver comment; M stays 77 (still three functionals). Record the R=200 verdict in
  [sbc-family-tiers.md](sbc-family-tiers.md).
- Equivalence re-record: [`fitViaNbinom`](../../benchmarks/R/equivalence.R) needs only its comments (yhat.test is
  now the log mean). DRAW-CHANGING on nbinom only; partition at the tip on the reference build: 54 of 55
  identical, nbinom the mover. P17 oracle (MANIFEST rule): negbin-exact.R FULL (both arms) and negbin-mixing.R at
  the tip on the shipped build, plus the SBC arm's R=200 verdict. Neutrality: bcf-equivalence and
  multinomial-equivalence bitwise; the four seeded-drift snapshot files carry no nbinom and must pass unchanged.
- tests/cpp: rewrite [`testNBDispersionGridConditional`](../../tests/cpp/test_model.cpp) against w_k above, with K_k
  checked against a direct lgamma sum and the weights against a per-row dnbinom-form sum;
  [`testNBSweepOrderAndRestore`](../../tests/cpp/test_model.cpp) and
  [`testActiveRowsNBKernels`](../../tests/cpp/test_model.cpp) for the anchor and the active-row K_k. New: c's formula
  (offset, all-zero floor), the working response equal to kappa/omega - a_i after refresh, setOffset at both
  updateScale values, restoreScale and dispersion-then-latents restore; the fit.scale round trip (exact) and the
  fitMax <= fitMin refusal.
- tinytest: [test-nbinom.R](../../inst/tinytest/test-nbinom.R) (ev = exp(link) draw by draw; leaf.scale = A; the
  recovery smoke's r band tightened to what the mixing gate supports; pre-change state refusal),
  [test-dispersion-channel.R](../../inst/tinytest/test-dispersion-channel.R) (the same identity),
  [test-augmentation.R](../../inst/tinytest/test-augmentation.R) (helper parity with the engine, psi less log r),
  [test-family-mutation-parity.R](../../inst/tinytest/test-family-mutation-parity.R) (creation with offset equals
  setOffset(updateScale = TRUE), bitwise; updateScale = FALSE keeps c and differs),
  [test-active-rows-pins.R](../../inst/tinytest/test-active-rows-pins.R) (the nbinom substitution must preserve the
  inactive rows' total count, since c is full-data, as the aft arm keeps its extremes),
  [test-calibration-prior-draws.R](../../inst/tinytest/test-calibration-prior-draws.R) and
  [test-fits-without-offset.R](../../inst/tinytest/test-fits-without-offset.R) (center c).
- Per the [RNG classes and their gates](README.md#rng-classes-and-their-gates): all exact gates in the
  exact-gates.yaml list (only nbinom's may move), sanitizers on tests/cpp and the R-loaded path (new numerics), and
  bench-sampler.R compare (hot path) including an nbinom configuration (added if absent), maintainer-run.

### Docs

- [negative-binomial.md](../design/negative-binomial.md): rewrite section 1 to the log-mean parameterization with
  logit-p kept as the recorded alternative and its measured defect; 2A's conditional (w_k, O(13 n)); section 3's
  caveat; section 4's reporting; section 5's state (fit.scale, the anchor, rebuild-not-shift); section 6's gate. In
  section 7 the log-mean door closes and the real-r door changes: CRT-Gamma conjugacy needs an r-free p_i, so under
  log mean a real r would move by a slice or MH step on the same collapsed likelihood, still waiting on a
  real-shape PG draw. Status line amended at landing.
- Rd: [bart.Rd](../../man/bart.Rd) (bartNegbin components and generics, the predict offset text),
  [dbartsFamilies.Rd](../../man/dbartsFamilies.Rd) (nbinom item: log link),
  [dbartsPriors.Rd](../../man/dbartsPriors.Rd) (anchor table: an nbinom row, "log mean", A),
  [dbartsAugmentation.Rd](../../man/dbartsAugmentation.Rd) (psi less log r), and
  [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) only if its updateScale text names gaussian.
- NEWS (nbinom is new in 1.0-0, so release state only): the families item reads "\code{"nbinom"}
  (negative-binomial counts on a log link, with an integer dispersion fixed or sampled on a grid)".
- TODO: negbin-real-dispersion's "the CRT-Gamma machinery drops in as specified" becomes the slice/MH note above;
  nbinom-dispersion-name and second-review-followups unchanged (the burn-in dispersion channel is now worth more,
  since r moves; noted, not done).
- Consumers: neither stan4bart nor bartCause uses nbinom (stan4bart builds gaussian and probit samplers; its
  neg_binomial_2 is unreachable rstanarm vocabulary; bartCause has no count path). No downstream change.

### Not in this plan

Collapsing bart2Negbin's per-sample loop; the dispersion name (nbinom-dispersion-name); real or sub-1 r; the
large-shape PG cost ([nbinom-large-shape-pg.md](nbinom-large-shape-pg.md), unaffected: it draws PG(y + r, psi)
whatever psi's parameterization).

## Steps

Each commit passes its gates before the next starts; landing per [Landing](README.md#landing).

1. C1, gate script only: benchmarks/R/negbin-mixing.R, unlisted. Gate: it FAILS at this commit (shipped engine),
   recorded in the commit body with its numbers.
2. C2, the model change, one commit because exact-gates runs negbin-exact.R on every push and tinytest pins the
   R identities: engine, bridge, R, defaultLeafScale and its C backstop, tests/cpp, tinytest, negbin-exact.R,
   exact-gates.yaml gains negbin-mixing.R, Rd, dbarts.h comments, NEWS, TODO, design doc sections. Gates:
   `--preclean` install; tests/cpp; full tinytest; negbin-exact.R and negbin-mixing.R FULL; every other exact gate
   quick; equivalence compare (54/55 identical, only nbinom moves, |z| expected large on its latent channel);
   bcf/multinomial bitwise; ASan/UBSan both halves; lint, air, check-rc-codoc, check-win-drift,
   check-doc-freshness; R CMD check --as-cran; bench-sampler.R compare (maintainer).
3. C3, SBC: sbc.R nbinom arm, burn ladder, sbc.yaml (waiver dropped, burn/thin/timeout), R=200 verdict recorded in
   sbc-family-tiers.md. Gate: the arm PASSES locally at R=200 with every functional ranked (no flag).
4. C4, equivalence re-record: new current baseline, MANIFEST row with partition, P17 oracle (C2's exact gates, the
   mixing gate, C3's verdict), neutrality and build mode; cpp-tests.yaml/equivalence.yaml filenames.
5. C5, records: this plan's Landing note, negative-binomial.md Status, the dec-B170 record mark per
   [decisions.md](../decisions.md) practice.

## Verification

    R CMD INSTALL --preclean -l $LIB .
    cd tests/cpp && make && ./test_bartcore
    R_LIBS=$LIB Rscript -e 'tinytest::test_package("dbarts")'
    R_LIBS=$LIB Rscript benchmarks/R/negbin-exact.R        # PASS, both arms
    R_LIBS=$LIB Rscript benchmarks/R/negbin-mixing.R       # PASS here; FAIL at C2's parent
    R_LIBS=$LIB Rscript benchmarks/R/equivalence.R compare # 54 identical, nbinom moves
    R_LIBS=$LIB Rscript benchmarks/R/sbc.R nbinom 200 150 <thin>   # no FLAG

## Risks

- r near the grid's edges: a truth below 1 still piles at r = 1 (the integer envelope, dec-B15), now visibly rather
  than hidden by the freeze.
- Residual slowness between r and tree structure at small n (a small r absorbing heterogeneity the forest could
  fit). The level ridge is gone by orthogonality, but this one is not ruled out; the SBC arm at n = 150 is where it
  would show, and a flag there is a finding to bring back, not a waiver to restore.
- r is weakly identified at large r and moderate means (30 against 50 above): the posterior there leans on the
  grid's cap and prior. The mixing gate avoids that cell on purpose; the help page should say r above about 20 is
  hard to tell apart at ordinary counts.
- Coverage of discrete counts: equal-tailed integer quantiles over-cover at small means, which the randomized PIT
  removes; a seed failing anyway is a finding, never a reason to widen the band.
- Cost, unmeasured: 13 n logOnePlusExp per sweep against sum_i (y_i + r) Devroye draws. A count of
  transcendental calls puts it well under the PG loop at the probe's means (y + r near 20) and comparable or above
  it when y + r is near 2. bench-sampler.R decides; the cheap exact fallback tabulates exp(eta_i) once per sweep
  and takes one log per row and grid point, halving the calls.
- c is full-data and fixed under updateScale = FALSE: an embedded user who creates on a placeholder response gets
  that placeholder's center until they re-anchor, as with gaussian.
- Pre-change nbinom states and saved fits from development builds stop loading (refused by fitMax <= fitMin). No
  release carried them.

## Open questions for the maintainer

Q1. The leaf prior on the log-mean scale: anchor A and the k default (prior sd of f is A/k). Numbers in
[Leaf prior](#leaf-prior).
- (a) Keep A = pi sqrt(3), k = 2 (sd 2.72): the width logit-p had. Worst on RMSE at both n and in r placement at
  n = 300.
- (b) A = 3, k = 2 (sd 1.5, probit's anchor): beats (a) everywhere; second-best coverage.
- (c) A = 2, k = 2 (sd 1.0): good RMSE at n = 2000, but under-covers the strong signal at n = 300 (0.78).
- (d) A = 3, k ~ chi(1.5, 2), the binary families' default form: lowest RMSE at both n, coverage at or near the
  best. Cost: nbinom joins probit and logistic in drawing k, so the ordinal-negbin-k-gap TODO's nbinom half (a k
  channel packageNegbinResults drops) becomes part of C2, about 20 more lines.
Recommendation: (d). It wins on accuracy without the coverage loss of (c), and it matches what the other
fixed-scale binary-link families already do. What would change it: a realistic count design where the k draw
under-covers, which would argue for (b).
