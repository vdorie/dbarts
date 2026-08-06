# Correlated errors and correlated outcomes around a BART mean

Status: RESOLVED 2026-07-22 (decision-gated door, TODO correlated-outcomes).
Fable scoping investigation, load-bearing API claims re-verified against
the code by the orchestrator; the covariance-structure hypothesis confirmed
and since VALIDATED empirically. (B) multivariate/SUR shipped as mvbart()
in stan4bart (commit e27a7c3, branch bartcore) - ZERO dbarts engine change.
(A) AR-1 stays deferred; see the landing note below.

## The capability and its causal motivation

dbarts fits a BART mean f(x) under scalar (per-observation, optionally
heteroscedastic) Gaussian error. Two causal settings want richer error/outcome
covariance:

- (A) AR-1 / serially-correlated errors: longitudinal, panel, interrupted-time-
  series, and difference-in-differences designs. Ignoring serial correlation
  mis-states treatment-effect standard errors (the Bertrand-Duflo-Mullainathan
  2004 serial-correlation-in-DiD critique).
- (B) Multivariate / multiple correlated outcomes: joint treatment effects with
  coherent joint credible regions and multiplicity-aware co-primary endpoints;
  causal mediation (the outcome+mediator bivariate system, with the error
  correlation rho as the sequential-ignorability sensitivity parameter);
  surrogate / principal outcomes (the joint law of surrogate and true outcome).

## The conduit dbarts already exposes (verified)

dbarts is designed to be a conditional BART mean inside an outer Gibbs/MH
sampler, and both the dbartsSampler reference class (R/dbarts.R) and the shipped
flat C API (inst/include/dbarts/dbarts.h) carry a complete DIAGONAL-precision
conditioning conduit:

- setResponse / setOffset / setWeights / setSigma (dbarts_sampler_setResponse /
  _setOffset / _setWeights / _setSigma): swap response, additive offset,
  per-observation PRECISION weights (y_i ~ N(f(x_i)+offset_i, sigma^2 / w_i),
  gaussian family only), and the held residual sd between sweeps. resid.prior =
  fixed() (dbartsFixedPrior) suppresses the internal sigma draw so the outer
  sampler owns sigma.
- run(0, 1) / dbarts_sampler_run: one sweep, returns the train fits the outer
  sampler differences into residuals.
- A per-sweep conditioning callback: dbarts_sampler_setCallback in dbarts.h, and
  the internal R-level C_dbarts_bartcore_runWithCallback (R_interface_bartcore.cpp,
  registered in R_interface.cpp) that rbart_vi already drives (R/rbart.R: a
  preStep that draws random intercepts and calls setOffset before each sweep) -
  one engine run with an R closure per sweep instead of a run(0,1) round trip.

What the conduit CANNOT express: a non-diagonal precision. Nothing accepts one,
and the constant-leaf sufficient statistics (sumWeights / residualVariance) are
per-leaf scalar. So the design rule is: anything whose CONDITIONAL precision is
diagonal composes today; off-diagonal correlation must be conditioned away in
the outer sampler.

## (B) Multivariate: pure composition, exact, today - NO new dbarts code

Seemingly-unrelated BART: one forest per outcome, e_i ~ N_q(0, Sigma). The
coupling is across EQUATIONS, not across observations of one f, so chained
Gaussian conditionals make each outcome's step ordinary: conditional on the
other outcomes' residuals,

    y_ik | . ~ N(f_k(x_i) + m_ik, v_k),  m_ik = Sigma_{k,-k} Sigma_{-k,-k}^-1 e_{i,-k},
    v_k = Sigma_kk - Sigma_{k,-k} Sigma_{-k,-k}^-1 Sigma_{-k,k}.

m_ik is a per-observation OFFSET, v_k a scalar sd - both already conduit inputs.
So the recipe is q dbartsSampler objects with resid.prior = fixed(), per-sweep
setOffset(m_k) / setSigma(sqrt(v_k)) / run(0,1), and a conjugate inverse-Wishart
Sigma draw in the outer sampler. Mediation is the triangular special case (Y's
equation conditions on observed m; rho on a sensitivity grid or with a prior).
Exact WITHOUT augmentation - observations stay independent, so no whitening
problem arises. dbarts exposes nothing new for correctness; native q-coupled
forests (riding the multi-forest / multinomial machinery) would buy only speed.

## (A) AR-1: cannot whiten a nonparametric mean

L^-1 y has mean L^-1 F, which mixes f across observations - BART's per-node
marginal likelihood and leaf draw assume the mean enters per observation, so you
cannot pre-whiten and hand the result to dbarts. Three routes:

- RECOMMENDED (outer-sampler feature, no new dbarts code): latent AR-1 state +
  white-noise nugget. y_t = f(x_t) + u_t + nu_t, u_t = rho u_{t-1} + eta_t. The
  outer sampler draws u by FFBS/HMC and (rho, variances) conjugately; the dbarts
  step is an ordinary iid-error BART fit conditioned on u via setOffset(u) and
  the nugget sd via setSigma/fixed(). This is literally the rbart_vi loop with
  the random intercept replaced by an AR process (docs/design/grouped-random-
  effects.md is the in-engine precedent). Caveat: this is AR-1-PLUS-NUGGET (the
  measurement-error model BSTS/CausalImpact use), not pure AR-1 errors, and f/u
  compete for the residual (the mixing hazard forest-ranef-interweaving.md
  studied for intercepts) - expect slow mixing as the nugget shrinks or rho -> 1.
- REJECTED: a pseudo-Gibbs conditional offset using the previous sweep's f - not
  a valid MCMC scheme (f in its own conditioning offset).
- EXACT pure-AR-1 (engine door, deferred): a GLS leaf update. For a fixed tree
  the conjugate leaf draw and move marginals become Z' Sigma^-1 Z (L x L,
  couples the leaves within a tree) and Z' Sigma^-1 R; AR-1 makes Sigma^-1
  tridiagonal so these stay O(n), but it needs a new banded/sparse
  setResidualPrecision channel plus cross-leaf sufficient-statistic kernels and
  a joint L-dimensional leaf draw in every move - breaking the per-leaf
  independence the constant-leaf kernels and their SIMD reductions assume.
  Generalizes to any banded/sparse precision (spatial CAR). Large; defer unless
  a consumer demands exact AR-1 errors.

## Priority and recommendation

- By raw demand, A leads (panel/DiD are ubiquitous) - but rbart_vi / stan4bart
  random intercepts already absorb the exchangeable share of within-unit
  correlation, blunting the worst of the BDM problem; AR-1 sharpens it mainly for
  longer panels. By value-per-new-code, B leads decisively (free, and unlocks
  mediation-with-sensitivity, a recurring bartCause/bairrtt ask).
- Sequencing: (B) ship as a validated composition recipe (a bartCause/bairrtt-
  side function or vignette, gated against a linear-SUR oracle) now; (A) add the
  AR-1 latent state to stan4bart/WALNUTS's parametric vocabulary (dbarts
  unchanged); leave exact pure-AR-1 as the recorded engine door above.
- Smallest dbarts-side enabler (only if profiling demands it): export the
  per-sweep R callback rbart_vi uses internally (C_dbarts_bartcore_runWithCallback)
  as a callback argument on dbartsSampler$run - it removes the per-sweep round-
  trip cost for every R-side outer-Gibbs composition (the multivariate recipe,
  custom rbart loops, AR prototypes) without touching the engine or the shipped C
  ABI, whose analogue dbarts_sampler_setCallback already exists.

Hypothesis CONFIRMED: the existing weighted-fit + between-sample-swap conduit
covers everything whose conditional precision is diagonal, and both models
condition into diagonal form (B exactly; A after latent augmentation). The
covariance belongs in the WALNUTS outer sampler.

## Open questions (verify at build time)

- setSigma composed with resid.prior = fixed(): the dbarts.h contract holds
  sigma "until the next call or gaussian draw"; confirm a per-sweep setSigma
  cleanly overrides a fixed() creation value.
- Scale anchoring: the internal [-0.5, 0.5] response rescale keys off
  creation-time min/max; a latent u or offset m_k growing beyond that range
  could compress the working scale. rbart_vi uses updateScale during warmup only
  - the recipes should adopt that convention; the failure mode is untested.
  COST (2026-08-06): a BCF host cannot adopt it. setOffset(updateScale = TRUE)
  is now refused on any multi-forest sampler, because both forests' leaf scales
  are calibrated against the creation-time transform and do not ride a rescale.
  The alternative - recomputing s and both leaf scales on every rescale - is
  rejected as a data-dependent prior refresh mid-run.
- BCF-mode mutability: whether a two-forest BCF sampler is driveable through the
  same per-sweep mutation conduit for multivariate-BCF (dbarts_results.logLikelihood
  is NaN for BCF, which matters for joint WAIC reporting).
- Literature (unverified, stated from memory): LongBet (panel causal BART),
  seemingly-unrelated / multivariate-Gaussian BART, AR-error tree models,
  tsBART (mean-side dependence). BDM 2004 is the standard DiD reference.
- The latent-AR mixing claim is by analogy to forest-ranef, not measured.

## Landing note (2026-07-22): (B) validated and productized

mvbart() in stan4bart (commit e27a7c3, branch bartcore) implements the
Priority section's (B) recipe exactly: q dbartsSampler objects with
resid.prior = fixed(), a pure-R outer Gibbs doing per-sweep setOffset
(conditional mean) + setSigma (conditional sd) + run(0,1), and a
conjugate inverse-Wishart draw of Sigma over the residual cross-products.
ZERO dbarts engine change - the conditional-mean conduit this doc
verified is complete and stable as specified.

VALIDATED against the conjugate MNIW / SUR oracle across rho in
{0, 0.5, 0.9}: recovers it, with a small explainable positive bias at
high rho. Beats independent per-outcome BART fits on held-out joint
log-score, by a margin that grows with rho - the coupling the recipe
exists to capture.

FOOTGUN recorded: dbartsFixedPrior's fixed(v) takes a VARIANCE, but
setSigma (R method and dbarts_sampler_setSigma alike) takes an SD.
Silently mixing the two fits the wrong residual scale with no error.
mvbart's outer loop passes sqrt(v_k) at the setSigma call site; any
other composition built on this recipe must do the same.

The per-sweep run callback (dbarts_sampler_setCallback /
C_dbarts_bartcore_runWithCallback, the "smallest dbarts-side enabler"
noted above) stays a PERF-ONLY optimization - it only removes the
run(0,1) round-trip cost - and is not freeze-gated: it may land any
time relative to 1.0-0 without an ABI change, since
dbarts_sampler_setCallback already exists in the shipped header.

(A) AR-1 stays deferred per the Priority section: its latent-state
recipe is destined for stan4bart/WALNUTS's outer loop, not dbarts.
