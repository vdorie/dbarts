# Robust (Student-t) errors: design

Status: LANDED 2026-07-17 (b4f818e). Plan: docs/plans/archive/robust-errors.md (this is its
step 1). Continuous responses gain a Student-t residual law by the classic
Gaussian scale-mixture augmentation, riding the per-observation precisions the
workingWeights hook already carries for logistic.

## 1. Augmentation and the lambda draw

Internal (rescaled) scale: working response z_i = yRescaled_[i]
([[src/bartcore/model.hpp#GaussianResponse::rescale]]), fit f_i = totalFits[i], residual
r_i = z_i - f_i, residual variance sigma^2, user precision w_i. The Student-t
error is the mean-1 Gamma scale mixture

  r_i | lambda_i ~ N(0, sigma^2 / (w_i lambda_i)),   lambda_i ~ Gamma(nu/2, nu/2),

(shape, rate; mean 1), marginal sqrt(w_i) r_i / sigma ~ t_nu. The lambda_i full
conditional is conjugate (likelihood proportional to
lambda_i^{1/2} exp(-w_i lambda_i r_i^2 / (2 sigma^2))):

  lambda_i | z, f, sigma ~ Gamma( (nu+1)/2,  (nu + w_i r_i^2 / sigma^2)/2 )

in (shape, rate), drawn via ext_rng_simulateGamma(rng, (nu+1)/2, scale) with
SCALE = 2 / (nu + w_i r_i^2/sigma^2) (the DartPrior.update idiom,
[[src/bartcore/model.hpp#TResponse::refreshLatents]]). Its mean (nu+1)/(nu + w_i r_i^2/sigma^2) is ~1 in-range, small
for a gross outlier -- that downweighting IS the robustness.

lambda enters the tree stage exactly as logistic's omega does: the engine weights
node sufficient statistics by ResponseModel::workingWeights() ([[src/bartcore/chain.hpp#Chain::run]];
into setNodeAverages [[src/bartcore/tree.hpp#Tree::setNodeAverages]] and MoveContext.weights [[src/bartcore/moves.hpp#MoveContext::weights]]), where
the t family returns the composite precision c_i = w_i lambda_i.
**No engine change, no GaussianResponse edit (item 2).** The constant leaf
recomputes (sumWeights, sumWeightedResponse, sumWeightedResponseSq) each sweep
from the current weights ([[src/bartcore/tree.hpp#Tree::setNodeAverages]]); the vector-leaf crossproduct cache
drops when workingWeightsVaryPerSweep() is true ([[src/bartcore/chain.hpp#Chain::run]]). Both were
built for logistic's per-sweep omega, so this is a pure ResponseModel addition.

## 2. Composition, families, scope

User weights compose multiplicatively into c_i = w_i lambda_i -- the rule
logistic uses spelled differently: it folds its integer trial-count weight INTO
omega by drawing PG(w_i, psi) as a sum of w_i PG(1, psi) variates
([[src/bartcore/model.hpp#LogisticResponse::refreshLatents]]), so workingWeights() returns count * latent
([[src/bartcore/model.hpp#LogisticResponse::workingWeights]]). For continuous Gaussian weights the product is a plain
multiply and the working RESPONSE stays z_i (unlike logistic's kappa/omega,
[[src/bartcore/model.hpp#LogisticResponse::refreshLatents]]). Offsets need no handling: rescale() subtracts the offset before
min/max ([[src/bartcore/model.hpp#GaussianResponse::rescale]]), so r_i is already offset-adjusted.

Recommended construction: a decorator TResponse holding a GaussianResponse plus
lambda_ and a composite buffer; refreshLatents draws lambda then calls
gaussian_->setWeights(composite) ([[src/bartcore/model.hpp#TResponse::refreshLatents]]) so BOTH workingWeights() and
drawSigma() see c_i unedited, with workingWeightsVaryPerSweep() true. This is the
AFTResponse pattern (a contained GaussianResponse over a mutable buffer, latents
redrawn each sweep, sigma delegated;
[[src/bartcore/model.hpp#TResponse::refreshLatents, TResponse::drawSigma]]) -- a closer precedent
than the plan's "grouped decorator." Only computeLogLikelihood needs a real
override (report dt, not the conditional dnorm); latents() exposes lambda_ for
the state block as probit exposes omega; setResponse/setData cold-init lambda_ to
1 and let the next sweep draw (the probit-latent convention).

Scope answer: **gaussian-only for v1;** probit/logistic have a fixed unit scale
and are untouched. A t-latent probit (robit; Liu 2004) is a legitimate
follow-up, NOT a trap: the same augmentation applies to the truncated latent and
ProbitResponse's null workingWeights ([[src/bartcore/model.hpp#ProbitResponse::workingWeights]]) would carry lambda. Caveat
-- a t_nu latent changes the LINK (t-cdf, not probit) and confounds latent scale
with nu under the sigma=1 identification, so robit needs its own default nu and
calibration, a separate design pass.

**Heteroscedasticity is refused.** A Student-t sampler cannot also carry a
heteroscedastic variance forest: the variance forest routes through the same
per-observation weight channel the scale mixture already owns, so the two
would multiply into one another. The refusal is one predicate,
`varianceForestIsRefused` (src/bartcore/facade.hpp), asked by both sampler
factories so the full-data and view paths cannot drift; a finite
`residualDf` is its own term there rather than a family test, since the t
augmentation reports the gaussian family. The host message names every term
at once: "invalid sampler specification: either the leaf covariate
designation is not one the leaf models can take, or a variance forest is
combined with a non-gaussian family, Student-t residuals, leaf covariates,
or a monotone constraint".

**The exported augmentation helpers dispatch on a law, not a family.**
`dbartsDrawLatents`/`dbartsWorkingResponse` and their flat twins
`dbarts_drawLatents`/`dbarts_workingResponse` draw the same six per-sweep
augmentation variables the engine draws, at a caller's fit. Keying that on
`ResponseFamily` would put the Student-t mixture on the `gaussian` arm - the
family a t sampler genuinely reports - where it is indistinguishable from a
plain gaussian response, which draws no augmentation variable at all.
`AugmentationLaw` (`src/R_interface_bartcore_common.hpp`) is a separate enum,
total over probit, logistic, ordinal, aft, nbinom and studentT, and is the
only thing either surface dispatches on. The family survives on that path at
exactly one point, `supportFamily`, which states the response support a law
constrains: studentT answers gaussian there because a t response is
unconstrained, not because the two share a law. Each surface maps into the
enum at its own boundary - a family name for the R helpers, whose argument is
a name, and `DBARTS_FAMILY_STUDENT` for the flat entries, the family list
being that surface's vocabulary for its create AND augmentation entries. A
list-derived law count backs a `static_assert` in each translation unit that
enumerates the laws, so adding one fails the build.

## 3. sigma^2 under the mixture

The sigma draw is the scaled-inverse-chi-square
[[src/bartcore/model.hpp#ChiSquaredScalePrior::drawSigmaSqFromPosterior]]: posterior
scale = df0 s0 + weighted SSR, df = df0 + numPositiveWeights, weighted
SSR = sum_i weights_i (z_i - f_i)^2 (misc_computeWeightedSumOfSquaredResiduals,
[[src/bartcore/model.hpp#ChiSquaredScalePrior::drawSigmaSqFromPosterior]]). Feed it c_i and it computes exactly the mixture draw:
SSR = sum c_i r_i^2, df = df0 + n (lambda_i > 0 always, so the positive count is
the positive-user-weight count, [[src/bartcore/model.hpp#GaussianResponse::countPositiveWeights]]). **The existing weighted sigma
path already does exactly this**, given c_i (item 2's setWeights delegation).
sigma is now the CONDITIONAL scale; the marginal residual variance is
sigma^2 nu/(nu-2) (item 6; a reporting trap in item 7).

## 4. nu: the fork (do not pick)

Fixed nu with a default. The augmentation is Geweke (1993, exact Gamma(nu/2,
nu/2) Gibbs); robust-t modeling is Lange-Little-Taylor (1989, nu typically 3-8,
often estimated) and West (1984). Defensible band nu in [3, 10]. nu = 4 (the
plan's claim) gives finite variance 2 sigma^2 and moderate tails; caveat:
kurtosis is infinite for nu <= 4, so a finite fourth moment wants nu > 4 (nu = 5
-> variance 5/3). Cost: one mis-specification knob (too small over-downweights
good points and inflates sigma; too large loses robustness), zero mixing cost,
no extra state.

Sampled nu on a grid. The DartPrior alpha-grid pattern
([[src/bartcore/model.hpp#DartPrior::initialize, DartPrior::update]]
precompute and discrete draw): a fixed nu grid, per-point lgamma constants
of Gamma(nu/2, nu/2) precomputed once, then per-iteration the log full conditional
  sum_i [(nu/2) log(nu/2) - lgamma(nu/2) + (nu/2) log lambda_i - (nu/2) lambda_i]
      + log p(nu)
from two scalar statistics (sum log lambda_i, sum lambda_i) at O(grid)
multiply-adds -- DART's alpha cost. Cost: one extra scalar reduction per sweep
plus a weakly-identified parameter (item 7). Out of scope for v1 per the plan;
the grid is its natural home.

## 5. R surface: candidates (recommendation per option; no global pick)

(a) A family value, family = "t". Precedent: family already resolves
gaussian/probit/logistic/aft ([[R/dbarts.R#dbarts]], [[R/bart.R#bart2]]). Cost: a scalar
match.arg cannot carry nu, so it needs a companion resid.df arg, splitting one
concept in two. Recommend only paired with (b)/(c) for the df.

(b) A residual-distribution constructor, resid.dist = student(df = 4), default
gaussian(). Precedent: the priors-as-objects vocabulary (tree.prior = cgm,
node.prior = normal, resid.prior = chisq; [[R/dbarts.R#dbarts]]) resolved by
[[R/model.R#parsePriors]] from [[R/model.R#dbartsPriors]]. student() composes df
now and a sampled-nu spec later with no new top-level argument. Do NOT overload
resid.prior -- that is the sigma^2 PRIOR (chisq/fixed, [[R/model.R#dbartsPriors]]),
orthogonal to the error law. Recommend: cleanest and most future-proof, a new
resid.dist slot.

(c) A bart2 scalar, resid.df = Inf (Inf = gaussian, finite = t). Precedent:
bart2's bare scalars sigdf/sigquant/k ([[R/bart.R#bart2]]). Cost: bart2-only
(dbarts() still needs the object form) and no room for sampled-nu without an
overload. Recommend as the convenience alias over (b), not the primitive.

## 6. Interactions

Range scaling's outlier sensitivity (the motivation). rescale() keys the
[-0.5, 0.5] map off the raw min/max of the offset-adjusted response
([[src/bartcore/model.hpp#GaussianResponse::rescale]]) and initialSigma = sigmaEstimate / range_
([[src/bartcore/model.hpp#GaussianResponse::GaussianResponse]]). One extreme point inflates range_, compressing the bulk into a
sliver and mis-anchoring both the leaf prior (scale/k on the compressed scale)
and the sigma prior. Robust errors let the model attribute that point to the tail
(lambda_i small) so it stops dragging f and inflating sigma. HONEST CAVEAT / plan
contradiction: they do NOT fix the compression -- range_ stays set by the raw
min/max at construction, which the t-model never revisits, so the plan's "direct
mitigation for range scaling's outlier sensitivity" (plan line 12) is only partly
true: it mitigates the outlier's LEVERAGE, not the SCALING. Follow-up if it bites:
a robust (quantile) range.

k / node.scale calibration. The k = 2 signal-to-noise argument assumes residual
variance sigma^2; under t_nu the marginal is sigma^2 nu/(nu-2), so a sample-sd
sigmaEstimate overstates the conditional sigma -- document sigma as the
conditional scale or derive sigmaEstimate robustly. Robustness shields leaves
from outliers, so default k is if anything more appropriate here.

Exact-posterior gate. Single tree, one predictor, few cuts, as the logistic and
categorical gates (benchmarks/R/logistic-reference.R, categorical-exact.R). Fix
nu, condition on sigma; enumerate tree structures and compute each leaf marginal
by 1-D quadrature over the leaf mean mu against the EXACT product of t_nu
densities times the N(0, (scale/k)^2) leaf prior -- the reference uses the
marginal t and never augments, so agreement validates the augmentation (free
sigma -> a 2-D quadrature, still tractable). Plus a component test: on fixed r_i
and sigma the lambda_i draw's moments match
Gamma((nu+1)/2, (nu + w_i r_i^2/sigma^2)/2), as PG moments were checked.

Equivalence fixture. A NEW scenario in benchmarks/R/equivalence.R (t- or
contaminated-normal data, resid.dist = student). Existing anchors untouched: the
frozen classic baselines cannot be re-recorded, and the t path adds no draws to
the gaussian stream (GaussianResponse::refreshLatents stays a no-op,
[[src/bartcore/model.hpp#GaussianResponse::refreshLatents]]), so the nine z-stats and every RNG-locked snapshot stay stable.

## 7. Costs and risks

- Lambda block mixing. lambda_i are conditionally independent (one vectorized,
  PG-cheap sweep), but lambda and sigma both scale the residual and mix slowly
  together when nu is small or outliers many -- the known scale-mixture Gibbs
  pathology. Remedy (out of scope): parameter expansion / interweaving (Liu-Wu,
  van Dyk-Meng).
- Prior sensitivity at small nu. As nu -> 1 (Cauchy) the variance vanishes and
  sigma is weakly identified, so sigquant dominates. Cap the grid (and any fixed
  value) at nu >= 2, ideally >= 3, to keep the marginal finite.
- Sampled-nu weak identification. The likelihood in nu flattens for large nu
  (30 vs 100 indistinguishable with few outliers); cap the grid, the prior
  matters.
- Reporting trap. getSigmas reports the CONDITIONAL scale; sigma^2 is NOT the
  marginal residual variance (nu/(nu-2) sigma^2). Consumers (stan4bart, rbart)
  may misread it -- document, or expose the marginal.
- Zero user weights compose cleanly: c_i = 0, lambda_i draws from its prior, and
  numPositiveWeights excludes the row ([[src/bartcore/model.hpp#GaussianResponse::countPositiveWeights]]), so the sigma df is right.

## 8. Resolution (VD, 2026-07-18)

A software/literature survey (brms, Bambi, R-INLA, PyMC, Stan wiki,
JAGS corpus, bayesreg, gamlss, heavy, MASS, statsmodels, the BART
package family) settled the forks:

- No BART implementation ships t errors; expectations import from
  general Bayesian regression, where the dominant packages (brms,
  Bambi, INLA, the BUGS tutorial corpus) ESTIMATE nu by default under
  a proper tail-bounding prior - gamma(2, 0.1) is the de facto
  convention (Stan prior-choice wiki, Juarez-Steel). Fixed-nu lineage
  exists (PyMC's canonical example nu = 3; heavy's default family is
  Student(df = 4); Lange-Little-Taylor advise fixed 4-6 when
  estimation is unstable) but no major Bayesian package defaults to
  it.
- The estimation pathologies are all about HOW: the nu likelihood is
  unbounded at the boundary (Fernandez-Steel 1999), flat priors give
  improper posteriors (Fonseca-Ferreira-Migon 2008), and the moment
  estimator is inconsistent (Wang-Ip). A capped grid under a proper
  prior - section 4's sampled arm - is the sanctioned form.
- Surface: constructor objects are the ecosystem idiom (brms
  student(), gamlss TF(), heavy Student(df = 4)); nobody splits a
  family string from a bare df argument.

DECIDED: both modes ship in one arc behind the section 5(b) surface.
resid.dist = student() estimates nu on a capped grid (proper prior;
the section 4 machinery); student(df = x) fixes it. Default remains
gaussian(). Docs carry the BART-specific caveat (a flexible mean
confounds tail inference; posterior-check estimated nu) and document
sigma as the conditional scale in both modes.
