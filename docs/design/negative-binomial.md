# Negative-binomial count outcomes: design

Status: LANDED 2026-07-18 (9c28b31). Plan: docs/plans/negative-binomial.md (this is
its step 1). Non-negative integer counts fit natively by the Polya-Gamma
negative-binomial augmentation (Polson-Scott-Windle 2013; Zhou-Li-Dunson-Carin
2012), riding the per-observation working weights the LogisticResponse port
already carries (src/bartcore/model.hpp:3517). The forest fits a log-odds latent
psi; a dispersion parameter r governs over-dispersion. Surfaced as
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

A negative binomial with dispersion (shape) r > 0 and per-observation success
odds set by the forest. Write the fit eta_i = f(x_i) + o_i for the sum-of-trees
f and offset o_i. Two parameterizations put eta in different places; they share
all machinery below but differ in what f MEANS and in whether the r update is
conjugate (the coupling to section 2).

**Logit-p (recommend): psi_i = logit(p_i) = f(x_i) + o_i.** The count law is

    y_i ~ NB(r, p_i),   P(y_i = k) = C(k + r - 1, k) (1 - p_i)^r p_i^k,

with p_i = plogis(psi_i). As a function of psi_i the likelihood is the exact
Polya-Gamma form

    p_i^{y_i} (1 - p_i)^r = e^{y_i psi_i} / (1 + e^{psi_i})^{y_i + r},

so (PSW/Zhou) omega_i ~ PG(y_i + r, psi_i) and, with kappa_i = y_i - (y_i+r)/2 =
(y_i - r)/2, the conditional is psi_i | omega_i ~ N(kappa_i/omega_i, 1/omega_i).
The backfitting engine therefore sees working response z_i = kappa_i/omega_i - o_i
= ((y_i - r)/2)/omega_i - o_i under per-iteration precision weights omega_i - the
LogisticResponse seam exactly, with two changes: the PG shape is y_i + r (real,
not 1) and kappa is (y_i - r)/2 (not y_i - 1/2). sigma is fixed at 1.

What f means to the user. The MEAN is

    E[y_i] = mu_i = r p_i/(1 - p_i) = r exp(psi_i) = r exp(f(x_i) + o_i),

so **f is a log-odds latent, NOT the log-mean**: log mu_i = log r + f(x_i) + o_i.
The offset enters the mean multiplicatively (mu_i scales by exp(o_i)), so
o_i = log(exposure_i) is a clean log-exposure offset - the plan's requirement
and the standard count-model convention, matching the gaussian offset's additive
role on ITS scale (here the log-mean scale). The honest cost: the level of f and
log r both shift the log-mean, so f's grand intercept and r confound (the
mean-level ridge; section 2's mixing caveat, robust-errors' lambda-sigma and
ordinal's f-vs-cutpoint analogs). Reporting mean counts requires r: a fitted or
predicted mean is r_draw * exp(f + o) per posterior draw, so r must be a reported
output (section 4).

**Log-mean (the alternative): log mu_i = f(x_i) + o_i, i.e. psi_i = f + o - log r.**
Now f IS the log-mean (the brms/Stan `neg_binomial_2` convention, section 3) and
the offset is log-exposure with the identical meaning. The mean no longer rides
r's level, so the interpretation is cleaner. Cost: psi_i now depends on r (through
- log r), which (a) makes the PG working response's shift move every time r moves
(a per-sweep re-anchor, cheap) and (b) degrades every r update: the CRT-Gamma
conjugacy BREAKS ((1 - p_i)^r = (r/(mu_i + r))^r is no longer const^r, so r's
full conditional is not Gamma), and the integer-grid conditional (section 2A)
loses its precompute economics (p_i moves with each candidate r_k, so the
per-sweep cost becomes O(n x grid) likelihood evaluations instead of one O(n)
reduction). Under log-mean the r update is Metropolis on log r or a
full-likelihood grid, both strictly worse than their logit-p forms.

Recommendation: **logit-p.** It is the Zhou-Carin / PSW native form: the PG tilt
psi IS the fit f + o, so the working-response identity is the LogisticResponse
seam unchanged and nothing re-anchors when r moves. Decisively, it keeps p_i
independent of r, which every good r update needs (section 2A): the integer-grid
conditional's lgamma kernel precomputes once and its per-sweep cost collapses to
one O(n) reduction (the ResidualDfPrior economics) only because the r-dependent
likelihood term separates, and the CRT-Gamma conjugacy behind the real-r door
requires it outright. The mean-interpretation cost (f is log-mean minus log r)
is a reporting concern, discharged by exposing r and computing mu = r exp(f + o)
at report time - not a modeling defect. The strongest argument against: users
reason about counts on the log-mean scale, and logit-p makes the directly-fit
quantity (log-odds) one step removed from the reported mean, with the
f<->log r level confounding a real (if mitigated) mixing cost. If that bites,
the log-mean surface is a documented follow-up reachable with the Metropolis r
update and no new augmentation.

**Response scaling / node.scale.** NB is a fixed-unit-scale family like probit
and logistic: it does NOT center or [-0.5, 0.5]-rescale the response (y are
counts, entering kappa directly), so fitScale = 1, fitShift = 0, sigmaScale = 1,
initialSigma = 1, and sigma stays fixed at 1 (drawSigma returns sigma), exactly
LogisticResponse (model.hpp:3581,3640-3643). Leaf-prior calibration therefore rides
node.scale, not a response range. psi is a log-odds, so v1 reuses **logistic's
node.scale = pi*sqrt(3)** (R/model.R:406; the logistic-latent sd pi/sqrt(3)
times probit's 3.0), giving a total-fit prior sd node.scale/k = pi*sqrt(3)/2 ~
2.72 that admits a plausible several-unit swing in the log-odds. Honest caveat
(the ordinal scheme-C / robust-errors k analog): the induced prior on the MEAN
depends on r, since mu = r exp(f); a fixed node.scale on the psi scale does not
fix the prior on counts, and at very small r (heavy over-dispersion) or extreme
exposure the log-odds swing per unit count changes. v1 documents the induced
prior and reuses the logistic constant rather than deriving an r-aware
node.scale; an r-aware calibration is a follow-up.

## 2. The r update and the exactness fork (the load-bearing decision)

RESOLVED (VD 2026-07-18): fork (A) - integer dispersion, fully exact; r
fixed or estimated on the capped grid by the closed-form discrete
conditional; real dispersion stays behind the recorded door carrying
fork (B)'s spec. VD's rider: design decisions should accommodate a
later real-r expansion where possible. Binding accommodations:
- The by-name state slot stores r as a real-valued scalar under a
  parameterization-neutral name ("dispersion"); grid mode writes
  integer-valued doubles, so a real-r mode later loads and saves with
  no state-format change.
- The R surface takes dispersion as a positive number; v1 refuses a
  non-integer fixed value informatively ("real dispersion is not yet
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
src/external/random.c:607; declared random.h:135), EXACT. LogisticResponse
handles an integer trial count w by SUMMING w independent PG(1, psi) draws
(model.hpp:3555-3557) - exact because PG(n, z) = sum of n PG(1, z) for integer
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
ResidualDfPrior pattern EXACTLY** (model.hpp:3994: precomputed per-grid-point
kernel, per-sweep scalar statistics, one discrete draw). Under logit-p it is
just as economical: for grid values r_k the log full conditional is

    log w_k = L_k + r_k * S + log prior_k,   S = sum_i log(1 - p_i),
    L_k = sum_c n_c [lgamma(c + r_k) - lgamma(r_k)]

with n_c the count histogram (y is FIXED, so every L_k precomputes ONCE at
construction, the ResidualDfPrior kernel_ move; the y_i log p_i term is
r-free and cancels in the normalization). Per sweep: ONE O(n) reduction for S
plus O(gridSize) multiply-adds - CHEAPER than CRT's O(sum y_i) Bernoullis, with
no count-scaling and no tuning. All PG shapes stay integer, the shipped
integer-sum Devroye path serves every draw bit-exactly, and NO new RNG
primitive is needed. Real r is deferred behind a recorded door (section 7)
pending either an exact real-shape primitive or an explicit project-level
decision to admit approximate MCMC. Cost, stated plainly: **r in (0, 1) - the
heavy-over-dispersion regime, variance > mu + mu^2 - is unrepresentable**, and
r between grid points is rounded to the grid; the ecosystem estimates
continuous dispersion everywhere (section 3), so integer-r is a genuine
modeling restriction, not just a discretization.

**(B) Real r estimated by CRT-Gamma, with the truncated fractional primitive,
DOCUMENTED as approximate MCMC.** The Zhou-Carin update: L_i ~ CRT(y_i, r) =
sum_{j=1}^{y_i} Bernoulli(r/(r + j - 1)) (0 when y_i = 0), then under an r ~
Gamma(a0, b0) prior (shape, rate) the full conditional is conjugate,

    r | {L_i}, {p_i} ~ Gamma(a0 + sum_i L_i, b0 - sum_i log(1 - p_i)),

one Gamma draw; CRT needs only integer table counts and holds for real r, and
the conjugacy requires p_i independent of r - logit-p only (section 1). The
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
no mainstream package restricts the dispersion's support, so (A) may be judged
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
(ResidualDfPrior: capped grid, gamma-kernel prior weights, model.hpp:3994-4023)
and turns the tail question into a grid-cap question. Proposed default,
PROVISIONAL pending recovery-gate calibration: grid {1, 2, 3, 4, 5, 6, 8, 10,
12, 15, 20, 30, 50} (dense where dispersion matters, sparse toward the
Poisson-like cap; cap justified because at r >= 50 the NB is practically
Poisson for moderate mu) with prior weights the gamma(2, 0.1) kernel
renormalized on the grid - noting honestly that gamma(2, 0.1)'s mean of 20 is
a real upward pull if left unexamined, which the capped renormalization tames
but does not justify; the location is provisional, marked for the recovery
gate to calibrate. Estimate r on the grid by default; allow a user-fixed
integer r. (Under fork (B) the same gamma(2, 0.1) serves as the continuous
CRT-conjugate prior; the prototype used it and recovered r across {0.5, 2, 10}.)
BART-specific caveat, verbatim from robust-errors: a flexible mean confounds
the dispersion (mean = r exp(f) shares its level with r); posterior-check the
estimated r, and document that r is weakly identified when counts are small or
the fit very flexible (the prototype's grid_sd shrinks only because psi is held
fixed there).

## 4. Surface

**Family value: `family = "nbinom"` (recommend).** The plan's choice, and it
matches R's own distribution vocabulary - `stats::dnbinom`/`rnbinom` - the same
token dbarts otherwise leans on. The ecosystem alternatives are `"negbinomial"`
(brms, VGAM) and `"negative.binomial"` (MASS); `"negbin"` (the task's working
name) is an abbreviation no major package uses. Recommend `"nbinom"` for the
R-core `dnbinom` alignment; `"negbinomial"` is the runner-up if brms-alignment is
valued over R-core. Added to the dbarts and bart2 family vectors (R/dbarts.R:383,
R/bart.R:700) and to the A_class whitelist (R/A_class.R:455-465).

**Response validation.** y must be non-negative integers. The check belongs in
the numeric-response branch of the family resolution (R/spec.R:199-210, where
the binary families run their 0/1 test): a "nbinom" arm beside the binary check,
refusing a non-integer or negative response by name with a message like
`family "nbinom" requires a non-negative integer (count) response`. This is a
REQUIREMENT, not a nicety: the NB pmf puts zero mass on non-integers, the grid
kernel's count histogram (section 2A) presumes integer y, and the (B)-door CRT
draw is only defined for integer y_i. Unlike the binary check, counts are
unbounded, so validation is integrality + non-negativity, not a two-value test.

**A third response-shape channel (the ordinal precedent).** resolveFamily
(src/R_interface_bartcore.cpp:1581-1614) branches on control.responseIsBinary and
control.numOrdinalCategories; a count response is neither, so - exactly as ordinal
added a K-level channel (docs/design/ordinal.md section 4) - NB needs a `count`
response-shape flag plumbed through ParsedControl beside responseIsBinary, with
resolveFamily accepting `"nbinom"` only on that shape and refusing it by name
everywhere else. The engine enum gains ResponseFamily::negbin (model.hpp:2580),
and the chain family switch (chain.hpp:597-631) a case constructing
NBResponse(y, offset, numObservations, rSpec), with the r spec (a fixed integer,
or the grid-estimate flag; the residualDf convention of "positive fixes,
non-positive estimates" carries over) threaded through the options struct as
residualDf and numCategories are (R_interface_bartcore.cpp:1759, 1762).

**Offset.** o_i = log(exposure_i), entering the mean multiplicatively (section 1);
a fixed-unit-scale family keeps its zero offset meaningful (R/spec.R:218-237),
as probit/logistic/ordinal do. No response de-scaling of the offset (fitShift = 0).

**Weights: refused in v1.** A frequency/case weight replicating an observation
w_i times would draw PG(w_i (y_i + r), psi) - inside fork (A)'s integer
envelope that stays exact for integer w (w (y_i + r) is integer), but it also
multiplies the count histogram into the grid kernel and the exposure question
into the likelihood, and the usual "weight" a count modeler reaches for is
EXPOSURE, which belongs in the offset (log-exposure), not in replication. v1
refuses weights by name at ingestion, beside the probit/logistic/ordinal weight
policy (R/spec.R:45-88), keeping the surface honest rather than guessing
which weighting the user meant. Door: integer frequency weights are EXACT under
fork (A) and cheap to add later (weight the grid statistics and the PG shape);
continuous weights inherit the real-shape question and wait on the section 7
weighted-binary fork.

**Prediction / reporting.** type = "bart"/"link" returns the single latent
column eta = f (as probit/logistic return their latent). type = "ev"/"response"
returns the MEAN counts mu = r exp(f + o) - which REQUIRES r, so the r draws are a
first-class posterior output, an n.samples-length `r` field (the count analog of
gaussian's sigma and ordinal's cutpoints; section 5). fitted()/predict() mean
shapes match the gaussian single-column ev shape (n x n.chains x n.samples where
extract does), computed per draw as r_draw * exp(f + o). predict requires
keepTrees (the predict.bart guard). A "prob-like" latent p = plogis(eta) may also
be exposed for diagnostics but the reported deliverable is the mean count.

**xbart / rbart_vi refusals.** xbart's mechanism is match.arg over its family
vector `c("auto", "gaussian", "probit", "logistic")` (R/xbart.R:26, matched at
:104): omitting "nbinom" from the vector makes match.arg itself the refusal,
BEFORE resolveClassificationFamily (:130) ever sees the value - and its losses
are misclassification/continuous, so a count loss (NB deviance / log-loss) is a
separate xbart pass. rbart_vi likewise omits "nbinom" from its family vector
(R/rbart.R:48) so its match.arg refuses, ahead of the group attribute build -
grouped NB is out of scope (section 7). Both refusals are the vector-omission
mechanism, the ordinal precedent.

## 5. State and mutation

**r in a new by-name scalar state block - the resid.df pattern EXACTLY.** Add the
virtual trio carriesR() / r() / restoreR() to ResponseModel (default false / 0 /
no-op), mirroring carriesResidualDf() / residualDf() / restoreResidualDf()
(model.hpp:4166-4169). r is a scalar, so it needs no length (the residualDf
analog, not the cutpoints vector analog); in grid mode the stored value is a
grid member, the TResponse estimatesResidualDf convention (model.hpp:4168).
ChainStateData gains a scalar field near its residualDf field, named
`dispersion` as shipped (retired: proposed as `r`; combiner.hpp:81-85, NaN
marking absent);
getState writes it when carriesR() (chain.hpp:3096-3100, the residualDf line);
stateIsValid refuses an NB state with a non-finite/non-positive r
(chain.hpp:3272-3275); setState restoreR()s it (chain.hpp:3746-3747). The bridge adds a
SLOT_DISPERSION enum (retired: renamed from SLOT_R) + name to slotNames
(R_interface_bartcore.cpp:6513-6515), a
conditional write when finite (:6642-6644, the resid.df line), and a by-name
read tolerating absence (:7117-7123). Old states omit the slot and load
unchanged - the whole point of the additive by-name block; no
state-format-version bump (additive, per the :6390-6392 rule).

**omega rides the existing latents slot.** The per-observation PG draws omega_i
are the latents, serialized through the existing `latents` slot exactly as
LogisticResponse's omega does (model.hpp:4438); latents() returns omega_.data().
**Restore-ordering REQUIREMENT: restoreR runs before restoreLatents.** The
working response is ((y_i - r)/2)/omega_i - o_i, so restoreLatents rebuilds
working from omega AND the current r; a restore that installs latents before r
rebuilds working against the stale r. setState must sequence the r block ahead
of the latents block (or restoreLatents must be the sole working-rebuild site
and restoreR must re-trigger it); this ordering is a stated contract of the
implementation, tested by a state round-trip. workingWeightsVaryPerSweep() is
true (per-sweep omega), dropping the sufficient-statistic caches each sweep
(chain.hpp:1538-1542), the logistic behavior.

**refreshLatents order - r FIRST, then omega (the invariance requirement).**
Per sweep: (1) update r from its full conditional given the fit - the grid
draw against {p_i} = plogis(psi_i) (fork (A)), or CRT-Gamma (fork (B)) - both
of which are COLLAPSED over omega: they condition on (y, f) only and never read
the omega draws; (2) draw omega_i ~ PG(y_i + r_new, psi_i) at the NEW r;
(3) rebuild the working response with r_new (kappa = (y_i - r_new)/2 over the
fresh omega). This order is what makes the scan a valid partially-collapsed
Gibbs sampler (van Dyk-Park): a draw made with a variable marginalized out
requires that variable be REGENERATED from its conditional under the new value
before anything conditions on it again. The reverse order (omega first, then r,
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
model.hpp:3397), and re-draw omega under the new y, rebuild working. setData
(n changes, everything stale): cold-init r to the grid median (the
ResidualDfPrior medianIndex convention, model.hpp:3998; or the user's fixed
value), rebuild the kernel, and cold-start omega at its PG(y+r, 0) mean
(y_i + r)/4 - the LogisticResponse coldStart generalization (model.hpp:
3682-3692, which uses w/4) - so the working response starts deterministic and
the first sweep's draw replaces it. setWeights is a no-op (weights refused);
setSigmaPrior a no-op (sigma fixed); setOffset shifts the working response by
the offset delta, keeping omega and kappa (the logistic setOffset,
model.hpp:3608-3613).

## 6. Gates

**Exact-posterior gate (single tree, small n, small counts).** In the single-tree
enumeration style of the logistic and ordinal gates. The NB category likelihood
is CLOSED FORM in (leaf log-odds, r) - the augmentation omega integrates out - so
the reference is omega-FREE and the gate quadratures only over the leaf log-odds
and sums over the r grid, never over omega, exactly as ordinal quadratures over
leaf means + gamma_2 and never over z. Concretely, under fork (A): use the
shipped grid and prior weights; enumerate the tree structures a single predictor
with a few cuts admits (root, or one split into two leaves); for each structure
the marginal is

    sum over grid r_k of  prior_k  x  integral over (mu_leaf...) of
      [ prod_i NB(y_i; r_k, plogis(mu_{node(i)} + o_i)) ]
      x prod_leaf N(mu_leaf; 0, (nodeScale/(k sqrt(numTrees)))^2),

a 1-2-D quadrature per grid point (the r dimension is a FINITE SUM - cleaner
than ordinal's continuous cutpoint integral). Renormalize each structure's
tree-prior x marginal over the enumeration. Match the sampler's posterior means
of the identified quantities - the leaf mean counts mu = r exp(f + o) and the
posterior distribution over grid r - to the reference to Monte Carlo error;
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
  w_k proportional to exp(L_k + r_k S + log prior_k) (section 2A) - the
  ResidualDfPrior drawIndex test pattern, with the L_k kernel checked against
  a direct lgamma evaluation.
- Behind the (B) door only: the fractional PG moment test (tolerance = the
  K-truncation bound of 2B(ii), NOT exactness - the test that makes the
  approximation's size a checked contract) and the CRT moment / conditional
  test (E[L_i] = sum_{j=1}^{y_i} r/(r + j - 1); the r histogram against the
  collapsed Gamma conditional, for which the prototype's grid agreement is the
  pilot).

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

- **Grouped / mixed-model NB (rbart_vi + nbinom).** Refused cleanly at the R
  layer (R/rbart.R:49,56) before the group attribute builds; rbart_vi's family
  vector omits "nbinom" for v1. The door is real and FEASIBLE: GroupedResponse is
  a base-response decorator whose conjugate group update needs a Gaussian working
  response on the latent scale (model.hpp:4706+), and NB HAS one - the PG working
  response z_i = kappa_i/omega_i is Gaussian given omega, exactly as logistic's
  is - so grouped NB is a coherent future once the r block and the group block are
  shown to interleave (the group draw would condition on the current omega and r).
  v1 refuses; the composition is recorded as feasible, the ordinal precedent.

- **Real (continuous) r.** THE door this note's fork creates. Fork (A) defers
  real r pending one of two unlocks: an exact real-shape PG primitive (a
  Devroye-style sampler with proven series bounds for real b - if one is
  published, the CRT-Gamma machinery of fork (B) drops in exactly, prototype-
  validated), or an explicit project-level decision to admit approximate MCMC,
  for which fork (B) is the fully specified plan (error budget, K sizing,
  bias-aware component-test tolerances). Opening it later needs no design
  revisit; this note IS the (B) design.

- **Log-mean surface.** The alternative parameterization (section 1), f =
  log-mean directly, is a documented follow-up reachable with a Metropolis
  update on log r (the prototype's MH arm) or an O(n x grid) full-likelihood
  grid, and NO new augmentation - only a re-anchoring of the working response's
  shift as r moves. Recorded as the escape if users need f on the count scale.

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
  and shipped (LogisticResponse already does them, model.hpp:3555-3557),
  fractional weights deferred behind the SAME door as real r, and the two
  doors should open together (one primitive, one bias budget, one component-
  test contract serves both - fork (B)'s specification is written to be that
  shared design). If VD overrides to (B) here, NB pays for the primitive and
  weighted-binary becomes its thin consumer, the original framing.
