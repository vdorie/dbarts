# Negative-binomial count outcomes: design

Status: LANDED 2026-07-18 (9c28b31); sections 1, 2A and 3-7 AMENDED by
nbinom-log-mean ([Landing](../plans/nbinom-log-mean.md#landing)), LANDED 2026-10-01 (fdfc1fe4): the forest models the log
mean and r is drawn given the means (dec-B170). Section 4 is also AMENDED by
[front-door](../plans/front-door.md#front-door) S2, LANDED 2026-09-09
(44b3fa6d): shape is not a `bart()`/`dbarts()` formal; it is the
`nbinom(shape = NA)` [`dbartsFamily`](../../R/family.R) constructor's
argument, and `family = "nbinom"` resolves to `nbinom()`'s default. The name
shape for r is nbinom-dispersion-name ([Landing](../plans/nbinom-dispersion-name.md#landing)), LANDED 2026-10-02 (32ea44b2;
dec-B189). Plan: docs/plans/archive/negative-binomial.md (this is
its step 1). Non-negative integer counts fit natively by the Polya-Gamma
negative-binomial augmentation (Polson-Scott-Windle 2013; Zhou-Li-Dunson-Carin
2012), riding the per-observation working weights the LogisticResponse port
already carries ([`LogisticResponse`](../../src/bartcore/model.hpp)). The forest fits the log mean
(section 1); a shape parameter r governs over-dispersion. Surfaced as
`family = "nbinom"`. The load-bearing resolution (section 2): exact PG draws
exist only for INTEGER shape, so v1 ships the exact envelope - r a positive
integer, fixed or estimated on a capped grid by a closed-form conditional (the
robust-errors nu-grid pattern) - with continuous r behind a recorded door.
Poisson (the r -> infinity limit) and zero-inflation are out of scope
(section 7). This note touches NO data layer
(docs/design/data-store.md, "Response family implementer": a family owns only
the working response/weights/latents channel); it is a pure ResponseModel
addition plus its family plumbing, the robust-errors and ordinal precedent.

## 1. The model and link (the parameterization fork)

A negative binomial with shape (size) r > 0 and a mean set by the forest.
RESOLVED twice: logit-p shipped first (2026-07-18) and was replaced by the
log-mean parameterization (dec-B170, 2026-10-01) after the third whole-branch
review measured that r never moved under logit-p. Both are written down here,
since they share every piece of machinery and differ in what f means and in
what the r step holds fixed.

**Log-mean (shipped): log mu_i = eta_i = f(x_i) + c + o_i.** The count law is
NB2, the MASS::glm.nb and Stan `neg_binomial_2` convention:

    y_i ~ NB(size = r, mu = mu_i),   variance mu_i + mu_i^2 / r,

with o_i = log(exposure_i) a log-exposure offset and c the response transform
below. In log-odds terms p_i = mu_i / (mu_i + r) and psi_i = logit(p_i) =
eta_i - log r, and as a function of psi_i the likelihood is the Polya-Gamma
form

    p_i^{y_i} (1 - p_i)^r = e^{y_i psi_i} / (1 + e^{psi_i})^{y_i + r},

so omega_i ~ PG(y_i + r, psi_i) and, with kappa_i = (y_i - r)/2,
psi_i | omega_i ~ N(kappa_i/omega_i, 1/omega_i). Since psi_i = f_i + a_i with
the anchor a_i = o_i + c - log r, the trees see the working response
z_i = kappa_i/omega_i - a_i under per-sweep precisions omega_i - the
LogisticResponse seam, with c and -log r entering exactly as an offset does.
sigma is fixed at 1. f is the log mean less c and the offset, so the reported
link (the train channel, type = "link") is eta = f + c + o, and the mean count is
exp(link) with no r in it ([`negbinMeanCounts`](../../R/bart.R)).

**Why log-mean: r mixes.** In the NB2 (mu, r) parameterization the Fisher
information is diagonal (d2l / dmu dr = (y - mu)/(mu + r)^2, mean zero), so
drawing r given the means costs nothing through them. Under logit-p, moving r
at fixed psi multiplies every mean by r'/r; the information on a common
log-mean shift is sum_i r mu_i / (mu_i + r), so at n = 1000 and r = 8 the grid
neighbour 8 -> 10 cost of order 10^2 nats, and the forest could follow r only
by a coordinated level shift of every tree, which tree-local moves never
propose. The review's probes (n 200 and 1000, mu = 8 exp(x1)) found every chain
frozen within burn-in, mostly at the cold start r = 8, whatever the true r
(3, 5, 30), with 90% predictive coverage of fresh counts 0.78. The same probe
under log-mean: the mixing gate (section 6) passes, r recovered at 5 and 2 and
coverage within 0.89-0.92.

**Logit-p (the recorded alternative): psi_i = f(x_i) + o_i, mu_i = r exp(psi_i).**
The form Zhou-Carin and Polson-Scott-Windle write down: the PG tilt is the fit
itself, p_i is free of r, and so the r conditional separates into a precomputed
kernel plus one O(n) statistic, and real r would have CRT-Gamma conjugacy. f is
then a log-odds, log mu = log r + f + o, and the level of f and log r are the
same direction in the likelihood - the ridge above, which made the cheap r
update useless in practice.

**Response transform and leaf prior.** c = log(max(sum_i y_i, 1/2) /
sum_i exp(o_i)) is the intercept of the Poisson model with offset, computed
over all rows (the aft precedent for a full-data transform), with the floor
keeping an all-zero response finite. It is [`NBResponse::fitShift`](../../src/bartcore/model.hpp);
fitScale and sigmaScale stay 1. The latent binary families center their leaf
prior at 0 because their link has a natural zero; a log mean has none, so the
gaussian and aft precedent applies on the link scale, with the log rate rather
than a midrange because counts include zeros. Spread: anchor A = 3 on the
log-mean scale ([`defaultLeafScale`](../../R/model.R)) with k drawn under
chi(1.5, 2) by default, as for probit and logistic (dec-B183; the probe and the
alternatives are in the plan's [Leaf prior](../plans/nbinom-log-mean.md#leaf-prior)). A named
sd is stated on the log-mean scale with no conversion.

## 2. The r update and the exactness fork (the load-bearing decision)

RESOLVED (VD 2026-07-18): fork (A) - integer shape, fully exact; r
fixed or estimated on the capped grid by the closed-form discrete
conditional; real shape stays behind the recorded door carrying
fork (B)'s spec. VD's rider: design decisions should accommodate a
later real-r expansion where possible. Binding accommodations:
- The by-name state slot stores r as a real-valued scalar under a
  parameterization-neutral name ("shape"); grid mode writes
  integer-valued doubles, so a real-r mode later loads and saves with
  no state-format change.
- The R surface takes shape as a positive number; v1 refuses a
  non-integer fixed value informatively ("real shape is not yet
  supported"), so admitting it later is a validation relaxation, not a
  signature change.
- In the engine, the omega loop calls a shape-parameterized PG helper
  (b = y_i + r; integer-sum implementation behind it today) and the r
  update sits behind its own small seam (grid conditional today, a
  CRT-or-real strategy later); neither the sweep body nor the working
  rebuild assumes integrality anywhere but inside those two seams.
- The exact gate quadratures r over the prior's support set, which is
  any finite grid unchanged; the r-first sweep order (section 5) is
  the valid order for BOTH modes.

Two SEPARATE questions hide here, and the plan/TODO conflate them.

**2A - how r itself moves** (fixed / grid / CRT-Gamma / Metropolis-on-log-r).
**2B - the RNG requirement for the MEAN update** (the real-shape PG gap).

They interlock through one fact: the mean update draws omega_i ~ PG(y_i + r,
psi_i) every sweep, and y_i + r is non-integer whenever r is - INDEPENDENT of
how (or whether) r is updated. **A fixed real r = 2.5 needs a non-integer-shape
PG draw every sweep exactly as an estimated real r does.** So the exactness
boundary is not fixed-vs-estimated; it is INTEGER-vs-REAL r, and that boundary
is what VD must pick (the three-way fork below).

### 2B: the real-shape Polya-Gamma gap, stated plainly

The shipped sampler is Devroye PG(1, psi) only (ext_rng_simulatePolyaGamma,
[`ext_rng_simulatePolyaGamma`](../../src/external/random.c); declared [`ext_rng_simulatePolyaGamma`](../../src/include/external/random.h)), EXACT. LogisticResponse
handles an integer trial count w by SUMMING w independent PG(1, psi) draws
([`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)) - exact because PG(n, z) = sum of n PG(1, z) for integer
n. Non-integer shape has no such reduction, and the facts are:

(i) With real r (fixed OR estimated), EVERY omega draw has non-integer shape
    y_i + r. There is no rare path: real r means approximate draws n times per
    sweep, every sweep.
(ii) The candidate fractional primitive, the Devroye/Zhou gamma-sum

         PG(a, z) = (1/(2 pi^2)) sum_{k>=1} g_k / ((k - 1/2)^2 + z^2/(4 pi^2)),
         g_k ~ Gamma(a, 1) iid,

     truncated at K terms, is APPROXIMATE: truncation drops a nonnegative tail,
     so the draw is systematically biased LOW. The bias is one-sided and
     bounded - the omitted tail's mean is (a/(2 pi^2)) sum_{k>K} 1/((k-1/2)^2 +
     z^2/(4 pi^2)) < a / (2 pi^2 (K - 1/2)), i.e. < 2.6e-4 absolute at K = 200
     for a < 1, relative bias < 2/(pi^2 (K-1/2)) ~ 1.0e-3 of E[PG(a, z)] ~ a/4
     at small z - and it shrinks only linearly in K (bias < eps needs K ~
     2/(pi^2 eps): 1e-6 costs K ~ 2e5 gamma draws per fractional part).
(iii) No established EXACT sampler exists for general real shape b. The
     ecosystem's own real-shape tools are all approximations: BayesLogit's
     rpg.gamma is this truncated sum, rpg.sp is a saddlepoint approximation,
     and its hybrid rpg routes large shapes to a normal approximation (section
     3). Windle-Polson-Scott (2014, arXiv 1405.0506) generalize Devroye's
     alternating series but BayesLogit itself does not use it as a general
     exact real-shape path; treating "the ecosystem has real-shape PG" as
     "exact real-shape PG exists" was this note's original error.
(iv) No gate would catch the bias. The exact-posterior gate (section 6) is
     omega-free - it quadratures the closed-form NB likelihood - and in an
     integer-r configuration it never exercises a fractional draw at all; in a
     real-r configuration its MC tolerance (~1e-2) dwarfs a ~1e-3 one-sided
     bias. The equivalence/bitwise gates check REPRODUCIBILITY, not
     correctness. The only honest witness is a PG-moment component test whose
     tolerance is set to the truncation bound (accepting the bias), not to
     exactness.

Prototype Part B (same script as below) CHARACTERIZES the truncated
composition's error rather than validating exactness: for non-integer b in
{0.5, 2.5, 3.3, 10.7}, z in {0.3..1.5}, K = 200, the sampled mean sits within
1.2% of the exact E[PG(b, z)] = (b/2z) tanh(z/2) at 4000 draws - consistent
with the ~0.1% one-sided truncation bias bound plus ~0.5-1% MC noise; the
numbers bound the approximation, they do not certify an exact sampler.

Cost model, honestly. Whatever the fork, the integer part alone is
floor(y_i + r) Devroye rejections per observation per sweep - **NB's PG cost
scales with the counts**, unlike logistic's one-draw-per-observation. The plan's
"PG draw per observation per iteration dominates, as logistic" is true in
STRUCTURE but understates the per-draw cost at large counts. The O(1)-per-draw
escape (saddlepoint) is itself an approximation, so within the exact envelope
there is NO large-count escape; that is a real, documented cost cliff.

### 2A: the three-way fork (VD's call)

**(A) v1 ships the EXACT envelope only: integer r, fixed or grid-estimated
(recommend).** Working the fixed-real-r fact through, the exact envelope is
precisely: r restricted to positive INTEGERS, whether fixed (user-supplied) or
estimated. Estimation cannot be CRT-Gamma (its Gamma full conditional yields
real draws, leaving the envelope), so the exact estimated mode is a **discrete
grid full conditional against the closed-form NB likelihood - the robust-errors
ResidualDfPrior pattern** ([`ResidualDfPrior`](../../src/bartcore/model.hpp): precomputed per-grid-point
kernel, one discrete draw). Under the log-mean model r is drawn given the means
eta_i = log mu_i, collapsed over omega; for grid values r_k the log full
conditional is

    log w_k = K_k - sum_i (y_i + r_k) log(1 + mu_i / r_k) + log prior_k,
    K_k = sum_c n_c [lgamma(c + r_k) - lgamma(r_k)] - Y log r_k,   Y = sum_i y_i,

with n_c the count histogram, so K_k precomputes once per response
([`NBShapePrior::computeKernel`](../../src/bartcore/model.hpp)). The rest
does not separate - p_i moves with r_k - so a sweep costs one exp per row and
one log1p per row and grid point, 13 n in all
([`NBShapePrior::drawIndex`](../../src/bartcore/model.hpp)), beside the PG
draw's sum_i (y_i + r) unit draws. (Under logit-p the conditional was
L_k + r_k S + log prior_k with one O(n) statistic S = sum_i log(1 - p_i), cheaper,
and useless for the reason section 1 gives.) No tuning. All PG shapes stay integer, the shipped
integer-sum Devroye path serves every draw bit-exactly, and NO new RNG
primitive is needed. Real r is deferred behind a recorded door (section 7)
pending either an exact real-shape primitive or an explicit project-level
decision to admit approximate MCMC. Cost, stated plainly: **r in (0, 1) - the
heavy-over-dispersion regime, variance > mu + mu^2 - is unrepresentable**, and
r between grid points is rounded to the grid; the ecosystem estimates
continuous shape everywhere (section 3), so integer-r is a genuine
modeling restriction, not just a discretization.

**(B) Real r estimated by CRT-Gamma, with the truncated fractional primitive,
DOCUMENTED as approximate MCMC.** The Zhou-Carin update: L_i ~ CRT(y_i, r) =
sum_{j=1}^{y_i} Bernoulli(r/(r + j - 1)) (0 when y_i = 0), then under an r ~
Gamma(a0, b0) prior (shape, rate) the full conditional is conjugate,

    r | {L_i}, {p_i} ~ Gamma(a0 + sum_i L_i, b0 - sum_i log(1 - p_i)),

one Gamma draw; CRT needs only integer table counts and holds for real r, and
the conjugacy requires p_i independent of r - logit-p only (section 1), so
under the shipped log-mean model this fork's r step would be a slice or
Metropolis move on the same collapsed likelihood instead (section 7). The
mean update then draws PG(y_i + r, psi) by integer-sum + truncated gamma-sum
fractional part, with K sized so the one-sided bias is provably below a stated
threshold (the (ii) bound: relative bias < 2/(pi^2 (K-1/2)); K = 200 -> ~1e-3,
documented in the family's docs and enforced by the component-test tolerance).
Costs: the CRT draw is O(sum_i y_i) Bernoullis per sweep with NO exact
large-count escape (the honest cliff: a normal/Poisson approximation to the
Bernoulli sum exists for large y_i but stacks a second approximation); the
fractional PG adds K gamma draws per observation; and - the deep cost - **this
would be the codebase's FIRST approximate-MCMC family**, breaking the exact-MCMC
uniformity the equivalence and exact-posterior gates are built around. The
honest argument for it: the entire PG ecosystem (BayesLogit, every published
NB-PG Gibbs including Zhou-Carin's and Neelon's own code) runs on exactly these
approximations, and the bias bound is smaller than any posterior feature a user
can resolve.

**(C) Exact-with-correction: investigated and EXCLUDED.** The candidate scheme:
augment only the integer part of the shape - omega_i ~ PG(y_i + floor(r),
psi_i), exact by integer-sum - leaving a leftover likelihood factor
(1 + e^{psi_i})^{-frac(r)} that the conjugate machinery does not see, and
Metropolis-correct every move against it. Written down, it fails structurally:
the leftover factor depends on psi_i and so multiplies into EVERY tree-stage
move - each leaf-mean draw and every grow/prune acceptance would need a
per-node MH correction over its member observations, converting the engine's
conjugate backfitting scan and closed-form marginal-likelihood ratios into
per-node Metropolis, an engine-wide change that cannot be expressed through the
ResponseModel seam (the data-store role contract confines a family to the
working response/weights/latents channel). A second candidate - independence-MH
on the fractional omega with the truncated gamma-sum as proposal - fails
because the proposal density (a K-fold convolution of scaled gammas) is
intractable, so the acceptance ratio cannot be computed. Recorded as excluded
with these reasons; if an exact real-shape PG primitive ever lands (a Devroye-
style alternating-series sampler with proven bounds for real b), the (A)->real
door opens without any of this.

**Evidence (prototype).** No-trees comparison of CRT-Gamma vs Metropolis-on-
log-r on synthetic NB data, psi held fixed (the logit-p isolation where both
target the same 1-D posterior), n in {200, 2000}, r in {0.5, 2, 10}, shared
Gamma(2, 0.1) prior, 6000 iterations / 1500 burn-in, seed 20260718, ESS by the
Geyer initial-positive-sequence estimator on 4500 kept draws. Script:
benchmarks/R/negbin-r-update-mixing.R.

    n     r    | ESS_CRT  ESS_MH  MHacc | post_CRT  post_MH  grid   grid_sd
    200   0.5  |   3339     754   0.43  |  0.607    0.607   0.607   0.079
    2000  0.5  |   3758    1023   0.47  |  0.474    0.475   0.474   0.022
    200   2.0  |   3298     947   0.40  |  2.013    2.011   2.012   0.152
    2000  2.0  |   2996     718   0.56  |  1.998    2.000   1.999   0.048
    200  10.0  |   2529    1161   0.46  | 10.524   10.534  10.529   0.359
    2000 10.0  |   2900    1039   0.46  |  9.995   10.004   9.994   0.109

- CORRECTNESS: CRT, MH, and the deterministic fine-grid posterior agree to ~3
  decimals in every cell, including r < 1: the r-UPDATE schemes themselves are
  exact. (The approximation in fork (B) lives in the PG mean update, which this
  isolation never draws - see the caveats.)
- MIXING: CRT-Gamma delivers 2.5-4x the ESS of Metropolis (2500-3760 vs
  720-1160) with no tuning knob. This ranking is why (B)'s r update is CRT, not
  MH; under fork (A) the grid conditional supersedes both (a direct draw from
  the discrete full conditional, no kernel to mix).

HONEST CAVEATS. (i) psi is held fixed, removing the f<->r level confounding a
real BART fit induces (mean = r exp(psi)); the ESS numbers are an OPTIMISTIC
bound on in-sampler r mixing, not a dbarts prediction - the transferable
findings are the RANKING (CRT > MH) and the correctness agreement. The
ordinal-mixing-study caveat verbatim. (ii) The isolation never draws omega at
all (the r comparison runs against the collapsed likelihood directly), so the
prototype can witness NEITHER the PG truncation bias NOR sweep-ordering bugs
between the r and omega draws (section 5's invariance argument) - it validates
the r-update kernels in isolation, nothing about their composition into the
sweep. (iii) Each cell is a single replicate; the ESS gap is order-of-magnitude
signal. (iv) grid_sd shrinks sharply because psi is known; real posteriors on r
are wider (the mean absorbs part of r).

**Decision.** Recommend **(A): v1 ships the exact envelope - r a positive
integer, fixed (user-supplied) or estimated on a capped integer grid by the
closed-form full conditional under a proper prior (the ResidualDfPrior pattern,
now an EXACT analogy: capped grid, precomputed kernel, per-sweep scalar
statistic), estimated by default; real r behind a recorded door.** It preserves
the project's exact-MCMC uniformity - the equivalence gates, the exact-posterior
gates, and the never-widen-tolerances discipline all presume the sampler
targets its posterior exactly, and admitting the first approximate family is a
project-level identity decision that deserves its own arc, not a rider on a
family note. It also ships with ZERO new RNG primitives and the cheapest r
update on the table. Strongest argument against: **integer r cannot represent
r < 1**, heavy over-dispersion (variance > mu + mu^2), a genuinely common
count-data regime - this note's own prototype headline includes r = 0.5 - and
no mainstream package restricts the shape's support, so (A) may be judged
too weak to ship as "negative binomial"; if VD weighs that regime above exact-
MCMC uniformity, (B) is the coherent alternative and its error budget is
specified above, ready to implement.

## 3. Survey (the r prior and the estimate-vs-fix default)

**PG real-shape samplers.** `pgdraw` (Makalic-Schmidt, CRAN) is Devroye and its
man page states `b` is an integer scalar/vector - it does NOT ship a real-shape
draw, so dbarts cannot lean on it. `BayesLogit` (Windle/Polson/Scott, CRAN) IS
the real-shape reference: `rpg.devroye(h, z)` (integer h only), `rpg.gamma(h, z,
trunc)` (the truncated gamma-sum series above, any real h, documented "slow"),
`rpg.sp(h, z)` (saddlepoint, any real h, O(1) per draw independent of h), and the
hybrid `rpg` that ROUTES on h: h in {1, 2} -> Devroye, ~13 < h <= 170 ->
saddlepoint, h > 170 -> normal approximation, else -> gamma-sum. The real-shape
methods are Windle-Polson-Scott 2014 ("Sampling Polya-Gamma random variates:
alternate and approximate techniques," arXiv 1405.0506): an alternating-series
generalization of Devroye plus the saddlepoint approximation. The design
implication (section 2B facts (iii)): the ecosystem's practical real-shape
paths are ALL approximations - the truncated gamma-sum, the saddlepoint, the
large-shape normal - and BayesLogit's own routing concedes it; there is no
established exact sampler for general real b to import. Cost side: naive
Devroye-summation is O(y_i + r) per observation, and the only O(1) escapes
(saddlepoint, normal) are approximate, so the exact envelope has no large-count
escape (the 2B cost cliff).

**Zhou-Carin lineage (the direct CRT precedent).** Zhou-Li-Dunson-Carin (2012,
ICML, arXiv 1206.6456) and Zhou-Carin (2015, IEEE TPAMI 37:307, arXiv 1209.3442)
introduce the CRT-Gamma r update for NB regression under logit-p: L_i ~ CRT(y_i,
r) (Stirling-number PMF, integer y_i, real r), and under r ~ Gamma(a0, rate h0)
the conjugate full conditional is Gamma(a0 + sum L_i, h0 - sum ln(1 - p_i)) =
Gamma(a0 + sum L_i, h0 + sum ln(1 + mu_i/r)) - verified against the paper's own
`Gamma(e0 + L, 1/(f0 - ln(1-p)))`. They ESTIMATE r under a DIFFUSE Gamma (a0, h0
~ 0.01) by default. Neelon (2019, Bayesian Analysis 14:849, "Bayesian
Zero-Inflated Negative Binomial Regression Based on Polya-Gamma Mixtures") is the
closest methodological template: the EXACT NB-PG + CRT machinery this note
builds, for GLM/spatiotemporal regression rather than BART.

**Standard Bayesian NB regression.** All use the log-mean NB2 convention (mu =
exp(eta), variance mu + mu^2/r) and ESTIMATE the dispersion: Stan
`neg_binomial_2(mu, phi)`; brms `negbinomial()` shape, log link, default prior
gamma(0.01, 0.01) (with a known move toward a tail-bounding PC-prior,
approx inv_gamma(0.4, 0.3), brms issue 1614, because the diffuse gamma is weakly
identifying near the Poisson limit); rstanarm `reciprocal_dispersion`; MASS
`glm.nb` theta (ML); PyMC `alpha`. WARNING for reporting: statsmodels
`NegativeBinomial` uses `alpha = 1/r` - the one inverse convention; dbarts's r is
the "size" (variance mu + mu^2/r), matching R's rnbinom `size`, Stan `phi`, brms
`shape`, MASS `theta`, PyMC `alpha` - NOT statsmodels alpha.

**NB-BART precedent.** Murray (2021, JASA 116:756, "Log-Linear BART," arXiv
1701.01503) fits Poisson / NB / multinomial / zero-inflated count BART - but
via a POISSON/GAMMA augmentation (the "gamma trick": phi_i ~ Gamma(n_i, sum_j
f^(j)), rendering each log-intensity a Gaussian-response BART problem), NOT
Polya-Gamma, with the NB dispersion kappa an explicit estimated parameter. So
Murray is the count-BART precedent but a DIFFERENT augmentation; its existence
means the log-mean-count-BART niche is filled by a non-PG method, and a PG+CRT NB
is the augmentation dbarts already speaks (LogisticResponse) rather than a new
gamma-trick engine. No shipped BART package (`BART` Sparapani-McCulloch has
wbart/pbart/lbart/mbart/surv but NO count; `bartMachine`, `stochtree`) fits an
NB count outcome at all. A dbarts NB-BART reusing its weighted-Gaussian tree
sampler through the PG working response with CRT for r appears to be a genuine
gap (Zhou-Carin, Neelon, and Murray each hold two of the three ingredients but
not this combination) - not a proven-novelty claim, a read of the searchable
literature; a very recent (2024-2026) paper cannot be ruled out.

**Prior and default, characterized like robust-errors' nu.** As with nu, no BART
sets the precedent and the general-Bayesian convention ESTIMATES the dispersion
under a proper prior. The historical default is the DIFFUSE Gamma(0.01, 0.01)
(Zhou, brms) - but the r likelihood FLATTENS toward the Poisson (large-r) limit,
so a diffuse prior leaves the upper tail data-undetermined, exactly robust-
errors' nu weak-identification pathology, and the ecosystem is moving to a
tail-bounding prior (brms's PC-prior direction). Precision matters here: with a
proper Gamma prior the posterior is ALWAYS proper (conjugacy or the bounded
likelihood guarantees it) - the risk is not impropriety but a posterior that
silently LEANS ON THE PRIOR where the data cannot distinguish r = 30 from
r = 300. Under fork (A) the mechanism is a CAPPED integer grid with a proper
renormalized prior over it - which makes the robust-errors nu analogy EXACT
(ResidualDfPrior: capped grid, gamma-kernel prior weights, [`ResidualDfPrior`](../../src/bartcore/model.hpp))
and turns the tail question into a grid-cap question. Proposed default,
PROVISIONAL pending recovery-gate calibration: grid {1, 2, 3, 4, 5, 6, 8, 10,
12, 15, 20, 30, 50} (dense where overdispersion matters, sparse toward the
Poisson-like cap; cap justified because at r >= 50 the NB is practically
Poisson for moderate mu) with prior weights the gamma(2, 0.1) kernel
renormalized on the grid - noting honestly that gamma(2, 0.1)'s mean of 20 is
a real upward pull if left unexamined, which the capped renormalization tames
but does not justify; the location is provisional, marked for the recovery
gate to calibrate. Estimate r on the grid by default; allow a user-fixed
integer r. (Under fork (B) the same gamma(2, 0.1) serves as the continuous
CRT-conjugate prior; the prototype used it and recovered r across {0.5, 2, 10}.)
BART-specific caveat: under the log-mean model r is orthogonal to the means
(section 1), so the level ridge logit-p had is gone, but r stays weakly
identified when counts are small (the variance mu + mu^2/r is then close to mu)
and at large r: at ordinary counts r = 30 and r = 50 are barely
distinguishable, the posterior there leans on the grid's cap and prior, and the
help says so. A flexible mean can still absorb part of the over-dispersion at
small n; the SBC arm at n = 150 is where that would show.

## 4. Surface

**Family value: `family = "nbinom"` (recommend).** The plan's choice, and it
matches R's own distribution vocabulary - `stats::dnbinom`/`rnbinom` - the same
token dbarts otherwise leans on. The ecosystem alternatives are `"negbinomial"`
(brms, VGAM) and `"negative.binomial"` (MASS); `"negbin"` (the task's working
name) is an abbreviation no major package uses. Recommend `"nbinom"` for the
R-core `dnbinom` alignment; `"negbinomial"` is the runner-up if brms-alignment is
valued over R-core. Added to the dbarts and bart2 family vectors ([`dbarts`](../../R/dbarts.R),
[`bart2`](../../R/bart.R)) and to the A_class whitelist ([`dbartsModel`](../../R/A_class.R)).

**Response validation.** y must be non-negative integers. The check belongs in
the numeric-response branch of the family resolution ([`resolveSamplerSpec`](../../R/spec.R), where
the binary families run their 0/1 test): a "nbinom" arm beside the binary check,
refusing a non-integer or negative response by name with a message like
`family "nbinom" requires a non-negative integer (count) response`. This is a
REQUIREMENT, not a nicety: the NB pmf puts zero mass on non-integers, the grid
kernel's count histogram (section 2A) presumes integer y, and the (B)-door CRT
draw is only defined for integer y_i. Unlike the binary check, counts are
unbounded, so validation is integrality + non-negativity, not a two-value test.

**A third response-shape channel (the ordinal precedent).** resolveFamily
([`resolveFamily`](../../src/R_interface_bartcore.cpp)) branches on control.responseIsBinary and
control.numOrdinalCategories; a count response is neither, so - exactly as ordinal
added a K-level channel (docs/design/ordinal.md section 4) - NB needs a `count`
response-shape flag plumbed through ParsedControl beside responseIsBinary, with
resolveFamily accepting `"nbinom"` only on that shape and refusing it by name
everywhere else. The engine enum gains ResponseFamily::nbinom (shipped under
that name, not the "negbin" this note originally proposed;
[`ResponseFamily`](../../src/bartcore/model.hpp)),
and the chain family switch ([`Chain::Chain`](../../src/bartcore/chain.hpp)) a case constructing
NBResponse(y, offset, numObservations, rSpec), with the r spec (a fixed integer,
or the grid-estimate flag; the residualDf convention of "positive fixes,
non-positive estimates" carries over) threaded through the options struct as
residualDf and numCategories are ([`optionsFromParsed`](../../src/R_interface_bartcore.cpp)).

**Offset.** o_i = log(exposure_i), entering the mean multiplicatively (section 1);
a fixed-unit-scale family keeps its zero offset meaningful ([`resolveSamplerSpec`](../../R/spec.R)),
as probit/logistic/ordinal do. The offset enters c (section 1), so an offset
swap with updateScale re-derives c, as a gaussian swap re-derives its range.

**Weights: refused, except 0/1 as the row mask (dec-B179).** Weights of 0
and 1 name the rows in the data set and install as the active-row mask, as on
probit and ordinal. Any other weight is refused. The original reasoning: a frequency/case weight replicating an observation
w_i times would draw PG(w_i (y_i + r), psi) - inside fork (A)'s integer
envelope that stays exact for integer w (w (y_i + r) is integer), but it also
multiplies the count histogram into the grid kernel and the exposure question
into the likelihood, and the usual "weight" a count modeler reaches for is
EXPOSURE, which belongs in the offset (log-exposure), not in replication. v1
refused weights by name at ingestion, beside the probit/logistic/ordinal weight
policy ([`enforceWeightPolicy`](../../R/spec.R)), keeping the surface honest rather than guessing
which weighting the user meant; dec-B179 kept that refusal for every weight
but 0/1. Door: integer frequency weights are EXACT under
fork (A) and cheap to add later (weight the grid statistics and the PG shape);
continuous weights inherit the real-shape question and wait on the section 7
weighted-binary fork.

**Prediction / reporting.** type = "bart"/"link" returns the log mean
eta = f + c + o per draw. type = "ev"/"response" returns the mean counts
mu = exp(link), which no longer read r. The r draws are still a first-class
posterior output, the `shape` field (the count analog of gaussian's sigma
and ordinal's thresholds; section 5), which ppd (rnbinom(size = r_s, mu = mu_s))
and loglik (dnbinom at the same pair) read. A drawn k rides `k` (or `sd`), as on
a bart fit. fitted()/predict() mean shapes match the gaussian single-column ev
shape. predict requires keepTrees (the predict.bart guard).

**xbart refusal.** xbart's mechanism is match.arg over its family
vector `c("auto", "gaussian", "probit", "logistic")` ([`xbart`](../../R/xbart.R), matched at
[`xbart`](../../R/xbart.R)): omitting "nbinom" from the vector makes match.arg itself the refusal,
BEFORE resolveClassificationFamily ([`xbart`](../../R/xbart.R)) ever sees the value - and its losses
are misclassification/continuous, so a count loss (NB deviance / log-loss) is a
separate xbart pass. The refusal is the vector-omission mechanism, the
ordinal precedent.

## 5. State and mutation

**r in a new by-name scalar state block - the resid.df pattern EXACTLY.** Add the
virtual trio carriesR() / r() / restoreR() to ResponseModel (default false / 0 /
no-op - retired: shipped as carriesShape() / shape() /
restoreShape() instead, [`NBResponse::carriesShape`](../../src/bartcore/model.hpp), [`NBResponse::shape`](../../src/bartcore/model.hpp), [`NBResponse::restoreShape`](../../src/bartcore/model.hpp)), mirroring carriesResidualDf() / residualDf() / restoreResidualDf()
([`TResponse::carriesResidualDf`](../../src/bartcore/model.hpp), [`TResponse::residualDf`](../../src/bartcore/model.hpp), [`TResponse::restoreResidualDf`](../../src/bartcore/model.hpp)). r is a scalar, so it needs no length (the residualDf
analog, not the thresholds vector analog); in grid mode the stored value is a
grid member, the TResponse estimatesResidualDf convention (retired: [`TResponse::estimatesResidualDfForTesting`](../../src/bartcore/model.hpp)).
ChainStateData gains a scalar field near its residualDf field, named
`shape` as shipped (retired: proposed as `r`; [`ChainStateData::shape`](../../src/bartcore/combiner.hpp), NaN
marking absent);
getState writes it when carriesR() ([`Chain::getState`](../../src/bartcore/chain.hpp), the residualDf line);
stateIsValid refuses an NB state with a non-finite/non-positive r
([`Chain::stateIsValid`](../../src/bartcore/chain.hpp)); setState restoreR()s it ([`Chain::setState`](../../src/bartcore/chain.hpp)). The bridge adds a
SLOT_SHAPE enum (retired: renamed from SLOT_R) + name to slotNames
([`storeState`](../../src/R_interface_bartcore.cpp)), a
conditional write when finite ([`storeState`](../../src/R_interface_bartcore.cpp), the resid.df line), and a by-name
read tolerating absence ([`setState`](../../src/R_interface_bartcore.cpp)). Old states omit the slot and load
unchanged - the whole point of the additive by-name block; no
state-format-version bump (additive, per the [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) rule).

**c rides the existing fit.scale block as (c, c + 1).** The response
transform is state: under updateScale = FALSE it is not recoverable from the
data. [`NBResponse::getScale`](../../src/bartcore/model.hpp) writes the pair and
[`NBResponse::restoreScale`](../../src/bartcore/model.hpp) decodes c = min, an
exact round trip (a midrange encoding such as (c - 1/2, c + 1/2) does not
round-trip for about 1% of values); fitScale stays 1 and never reads the width.
[`Chain::stateIsValid`](../../src/bartcore/chain.hpp) refuses an nbinom state
whose pair is not increasing, which is every state written under logit-p
((0, 0)) and the case the chain's restore would otherwise skip silently. No new
block and no version bump: no state format has shipped.

**Rebuild, never shift.** Whenever c, r or the offset change without an omega
draw - setOffset, restoreScale, restoreLatents - the working response is
REBUILT from omega as kappa_i/omega_i - a_i through the one expression the draw
uses, so creation with an offset and setOffset(updateScale = TRUE) to it are
bitwise the same state.

**omega rides the existing latents slot.** The per-observation PG draws omega_i
are the latents, serialized through the existing `latents` slot exactly as
LogisticResponse's omega does ([`NBResponse::latents`](../../src/bartcore/model.hpp)); latents() returns omega_.data().
**Restore-ordering REQUIREMENT: restoreR runs before restoreLatents.** The
working response is ((y_i - r)/2)/omega_i - a_i, so restoreLatents rebuilds
working from omega AND the current r and c; a restore that installs latents before r
rebuilds working against the stale r. setState must sequence the r block ahead
of the latents block (or restoreLatents must be the sole working-rebuild site
and restoreR must re-trigger it); this ordering is a stated contract of the
implementation, tested by a state round-trip. workingWeightsVaryPerSweep() is
true (per-sweep omega), dropping the sufficient-statistic caches each sweep
([`Chain::run`](../../src/bartcore/chain.hpp)), the logistic behavior.

**refreshLatents order - r FIRST, then omega (the invariance requirement).**
Per sweep: (1) update r from its full conditional given the log means
(section 2A), COLLAPSED over omega: it conditions on (y, f) only and never reads
the omega draws; (2) draw omega_i ~ PG(y_i + r_new, psi_i) at the NEW r, with
psi_i = f_i + a_i at the new log r; (3) rebuild the working response with r_new.
The sweep is then a two-block Gibbs sampler: block 1 draws the trees given
(omega, r, y), block 2 draws (r, omega) jointly given the trees as
p(r | f, y) p(omega | r, f, y), and a joint draw of a block is a valid Gibbs
step. What the r step holds fixed is the means, not psi. The reverse order (omega first, then r,
then rebuild) is NOT invariant: omega would carry shape y + r_old while the
tree stage consumes kappa built from r_new - the trees then condition on an
omega that has the wrong distribution given the state they see. The first draft
of this note had exactly that bug, rationalized as "conditioning the CRT on the
fresh fit" - vacuous, since the collapsed r update never reads omega. Note the
prototype could not have caught this: its isolation never draws omega at all
(section 2 caveat (ii)), so sweep-ordering errors are invisible to it; the
exact-posterior gate (section 6), which runs the full composed sweep against an
augmentation-free reference, is the gate that would. sigma is ignored (fixed
at 1) in all three steps.

**setResponse / setData cold-init.** setResponse (same n, new y - the embedded-
Gibbs count swap): KEEP the current r (a slow-moving global the outer sampler
wants persisted across a small y perturbation - the ordinal kept-cutpoints
clause), RECOMPUTE the grid kernel L_k (it derives from the count histogram,
which the new y changes - the ordinal computeScales-on-setResponse precedent,
[`OrdinalResponse::setResponse`](../../src/bartcore/model.hpp)), and re-draw omega under the new y, rebuild working. setData
(n changes, everything stale): cold-init r to the grid median (the
ResidualDfPrior medianIndex convention, [`ResidualDfPrior`](../../src/bartcore/model.hpp); or the user's fixed
value), rebuild the kernel, and cold-start omega at its PG(y+r, 0) mean
(y_i + r)/4 - the LogisticResponse coldStart generalization
([`LogisticResponse::coldStart`](../../src/bartcore/model.hpp), which uses w/4) - so the working response starts deterministic and
the first sweep's draw replaces it; setData recomputes c, setResponse and
setOffset recompute it only under updateScale = TRUE (FALSE, the
embedded-Gibbs default, keeps it, as gaussian keeps its range). setWeights is a
no-op (weights refused; 0/1 weights reach the mask instead); setSigmaPrior a
no-op (sigma fixed); setOffset keeps omega and kappa and rebuilds the working
response under the new anchor.

## 6. Gates

**Exact-posterior gate (single tree, small n, small counts).** In the single-tree
enumeration style of the logistic and ordinal gates. The NB category likelihood
is CLOSED FORM in (leaf log mean, r) - the augmentation omega integrates out -
so the reference is omega-FREE and the gate quadratures only over the leaf log
mean and sums over the r grid, never over omega, exactly as ordinal quadratures over
leaf means + gamma_2 and never over z. Concretely, under fork (A): use the
shipped grid and prior weights; enumerate the tree structures a single predictor
with a few cuts admits (root, or one split into two leaves); for each structure
the marginal is

    sum over grid r_k of  prior_k  x  integral over (m_leaf...) of
      [ prod_i NB(y_i; size = r_k, mu = exp(m_{node(i)} + o_i)) ]
      x prod_leaf N(m_leaf; c, (nodeScale/(k sqrt(numTrees)))^2),

a 1-2-D quadrature per grid point (the r dimension is a FINITE SUM - cleaner
than ordinal's continuous cutpoint integral). Renormalize each structure's
tree-prior x marginal over the enumeration. Match the sampler's posterior means
of the identified quantities - the cell mean counts exp(m) and the posterior
distribution over grid r - to the reference to Monte Carlo error; the estimated
arm carries a two-level exposure offset so the anchor's offset, c and log r
terms are all exercised, and the script checks the engine's c against its
formula first ([negbin-exact.R](../../benchmarks/R/negbin-exact.R));
tolerances bound MC plus quadrature error and are never widened to pass.
Agreement validates the PG mean augmentation, the grid r update, AND their
composition into the sweep (the section 5 ordering; an invalid scan shifts the
stationary law, which this gate CAN see) - the robust-errors / ordinal
reference-never-augments logic. STATED LIMITATION: because the reference is
omega-free and fork (A) draws only integer-shape (exact) PG variates, this gate
exercises NO approximate code path; if the (B) door ever opens, the gate gains
a real-r cell but its MC tolerance (~1e-2) CANNOT resolve the ~1e-3 truncation
bias - the gate does not certify the fractional primitive, only the PG-moment
component test below does, at a tolerance honestly sized to the truncation
bound (section 2B(iv)). The bitwise/equivalence gates likewise witness
reproducibility, never correctness.

**Component tests.**
- Integer-shape PG moment (the mainline): on fixed integer b = y + r at the
  shapes the gate uses, the integer-sum draw's mean and variance match the
  analytic PG moments (mean (b/2z) tanh(z/2); the PG-moments precedent) - plus
  the stream identity that NB at b = 1 consumes exactly the shipped PG(1)
  Devroye stream (the ordinal K = 2 identity analog: NB's PG path IS
  logistic's, generalized only in the summation count).
- Grid r conditional: on a tiny fixed (y, {p_i}) the sampled grid-index
  histogram matches the hand-computed discrete full conditional
  w_k (section 2A) against a per-row dnbinom-form sum - the ResidualDfPrior
  drawIndex test pattern, with the K_k kernel checked against a direct lgamma
  evaluation less Y log r_k.
- Behind the (B) door only: the fractional PG moment test (tolerance = the
  K-truncation bound of 2B(ii), NOT exactness - the test that makes the
  approximation's size a checked contract) and the CRT moment / conditional
  test (E[L_i] = sum_{j=1}^{y_i} r/(r + j - 1); the r histogram against the
  collapsed Gamma conditional, for which the prototype's grid agreement is the
  pilot).

**Mixing gate.** The exact gate checks the stationary law on one tree and
n = 50, which cannot see a chain that never leaves its cold start: logit-p
passed it while frozen. [negbin-mixing.R](../../benchmarks/R/negbin-mixing.R)
fits the default forest with two chains at mu = 8 exp(x1), r0 = 5 at n = 2000
and r0 = 2 at n = 500, and requires each chain to leave the cold start, split-Rhat
on r below 1.05, the pooled 95% set to cover r0, and 90% predictive coverage
of fresh counts by randomized PIT within 0.90 +- 0.04. It fails under logit-p
(every chain at r = 8, coverage 0.71-0.87) and passes under log-mean.

**Recovery.** Simulated NB counts over a nonlinear f at moderate n, checking
mean-count calibration and r recovery against truth across grid values r in
{2, 5, 10} plus an off-grid truth (r = 7, recovered to the bracketing grid
mass) and - the honest boundary probe - an r = 0.5 truth, DOCUMENTING what fork
(A) does when the data are more dispersed than the grid can say (mass piles at
r = 1 and the mean fit absorbs what it can): the failure mode is recorded, not
hidden. The family-level smoke beyond the exact gate.

**Equivalence fixture and neutrality.** A NEW scenario in
benchmarks/R/equivalence.R (count response, family = "nbinom") recording the NB
channels: mean counts, r draws, and the omega latents. Existing anchors are
untouched - NB is a new family behind a new enum value and a new response-shape
flag, adds NO draw to any existing family's stream, and does not touch the
gaussian/probit/logistic paths - so the frozen baselines and every RNG-locked
snapshot stay stable; the neutrality trail is verified by re-running
equivalence.R compare and expecting IDENTICAL draws for the existing families (no
re-record), the robust-errors and ordinal precedent. Component C++ tests
(tests/cpp) cover the integer-shape PG-moment and grid-conditional checks; the
cross-ISA PG stream gate extends to the integer-sum path (and to any fractional
primitive only if the (B) door opens).

## 7. Out of scope, and the doors

- **Zero-inflation / hurdle NB.** A two-component (structural-zero + count) model
  is a multi-forest hurdle (docs/design/multinomial.md machinery; the Neelon ZI/
  hurdle lineage), not a single-forest family. OUT OF SCOPE. The door: once NB
  ships as a count family, a hurdle wrapper composing a binary "at-risk" forest
  with an NB count forest is a coherent follow-up, recorded, not designed here.

- **Poisson.** The r -> infinity limit. Deferred until NB lands (the plan's
  sequencing); a Poisson family is a separate log-mean augmentation (no PG shape
  parameter), revisited after NB. Recorded.

- **Mixed-model NB.** Not dbarts's: multilevel structure is stan4bart's
  ([The decision](retire-grouped-random-effects.md#the-decision)).

- **Real (continuous) r.** THE door this note's fork creates. Fork (A) defers
  real r pending one of two unlocks: an exact real-shape PG primitive (a
  Devroye-style sampler with proven series bounds for real b - if one is
  published), or an explicit project-level decision to admit approximate MCMC,
  for which fork (B) specifies the PG side (error budget, K sizing, bias-aware
  component-test tolerances). Its r step changes under the log-mean model:
  CRT-Gamma conjugacy needs an r-free p_i, which logit-p had and log-mean does
  not, so a real r would move by a slice or Metropolis step on the same
  collapsed likelihood section 2A's grid draws from, still waiting on a
  real-shape PG draw for the mean update.

- **Log-mean surface.** CLOSED: it is the shipped model (section 1,
  dec-B170).

- **dbarts.h exposure.** NONE in v1, the robust-errors / ordinal precedent. The
  flat C API (inst/include/dbarts/dbarts.h) is unchanged; NB is reachable only
  through the R surface and the internal bartcore .Call path, so no LinkingTo
  consumer (stan4bart) sees an ABI change. Door: a future dbarts.h entry could
  expose the count family for embedded use, deferred until demand.

- **The weighted-binary implication (state it explicitly).** The real-shape
  PG(b, z) gap (section 2B) is the SAME gap weighted-binary's real-weights half
  faces: a weighted logistic with a non-integer case weight w needs PG(w, psi)
  for real w, which the shipped integer-sum cannot draw. What this note's
  resolution implies has INVERTED from its first draft: because fork (A)
  restricts to integer shapes, **NB lands NO real-shape primitive** - so the
  shared item between NB and weighted-binary is no longer a primitive one of
  them builds for the other, but the DECISION both are gated on: does the
  project admit approximate PG draws (no exact real-shape sampler exists,
  2B(iii)), or does it hold the integer-exact line? Weighted-binary's
  real-weights half inherits fork (A)'s answer verbatim: integer weights exact
  and shipped (LogisticResponse already does them, [`LogisticResponse::refreshLatents`](../../src/bartcore/model.hpp)),
  fractional weights deferred behind the SAME door as real r, and the two
  doors should open together (one primitive, one bias budget, one component-
  test contract serves both - fork (B)'s specification is written to be that
  shared design). If VD overrides to (B) here, NB pays for the primitive and
  weighted-binary becomes its thin consumer, the original framing.
