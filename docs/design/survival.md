# Survival models: design

Status: proposed 2026-07-10. First survival family (AFT log-normal) for
1.0-x. Companion to the roadmap in docs/plans/survival-models.md and the
extensions table in docs/design/core-generalization.md ("Survival (AFT,
discrete-time hazard), ordinal, quantile: ResponseModel latents;
person-period expansion at ingest").

Discrete-time hazard is a follow-up section against this same note; it is
person-period ingestion sugar over the existing binary families and lands
after AFT. This note covers AFT only.

## Which model first

AFT log-normal, per the plan's default. It is the survival family the
ResponseModel latent seam already provisions for, it reduces exactly to
the Gaussian path (a free correctness gate), and it is riAFTBART's model,
the named consumer whose outer loop the grouped decorator replaces.

## Model

Accelerated failure time, log-normal error:

  log T_i = f(x_i) + offset_i + sigma * eps_i,   eps_i ~ N(0, 1).

The forest models the conditional mean of log T on the log scale; sigma
is the log-scale residual sd, drawn conjugately exactly as Gaussian's.
Right-censored observations (status_i = 0) observe only a lower bound
C_i on T_i: their latent log T_i is redrawn each sweep from
N(f(x_i) + offset_i, sigma^2) truncated BELOW at log C_i, the truncated
normal already in external/random.h
(ext_rng_simulateLowerTruncatedNormal). Uncensored observations
(status_i = 1) carry log T_i fixed at the observed log event time and
enter as ordinary Gaussian data. Right-censoring only in v1; the latent
code must not preclude left/interval censoring later.

## Engine: AFTResponse (src/bartcore/model.hpp)

AFTResponse contains a GaussianResponse over a mutable log-time buffer
logT_ (the Gaussian borrows logT_.data() as its response). Each sweep,
refreshLatents redraws logT_ for the censored observations from the
truncated normal, using the public fitScale()/fitShift()/offset()
accessors to map the internal-scale fit back to the log scale, then
rebuilds the Gaussian working response via setResponse(updateScale =
false). drawSigma, workingResponse, and the scale accessors all delegate
to the contained Gaussian, so the freshly drawn latent times are data
for the conjugate sigma draw (loop order: refreshLatents then drawSigma,
chain.hpp).

REDUCTION PROPERTY (the free gate): with zero censored observations
refreshLatents is a complete no-op (early return), construction feeds the
Gaussian the identical log-time buffer, and every other virtual
delegates. So an AFT fit on all-uncensored log T is BIT-IDENTICAL to a
Gaussian fit on the same log T with the same seed/settings. R passes the
engine log-times already transformed (log() at ingest), so the reduction
target and the AFT fit start from the identical buffer -- no in-C log to
diverge from R's.

State: AFT rides the EXISTING serialization blocks -- `latents` (logT_,
via latents()/restoreLatents, the probit/logistic channel) and
`fit.scale` (the Gaussian min/max, via getScale/restoreScale). No new
block, no stateFormatVersion bump; the registry rule
(docs/design/public-surface.md 2) is satisfied trivially (an additive
change needs no new name when existing blocks already carry the state).

Family plumbing: a new ResponseFamily::aft; sigmaIsFixed_ treats aft like
gaussian (drawn unless the user fixes it); status threads to the chain as
a borrowed SamplerOptions.survivalStatus vector (copied at construction
into a censored-index list, the groupIndices precedent), delivered by the
bridge from an internal control attribute (bartcore.survival), the
applyGroupAttribute pattern.

## R surface (stage 4, kept separable)

TASTE CALLS, each with a recommendation; the reviewer may revise with VD
before merge, so stage-4 surface code stays cleanly separable from the
engine.

1. Response ingestion. RECOMMEND accepting both: a survival::Surv object
   (detected by inherits(y, "Surv"), NOT importing survival -- Surv is a
   two-column matrix with a "type" attribute; we read time = column 1,
   status = column 2, and refuse types other than "right"), and a plain
   two-column (time, status) matrix / cbind for users without survival.
   Both transform to (log(time), status) at ingest.

2. Family naming. RECOMMEND auto-dispatch on a Surv response (a Surv y
   selects AFT with no family argument), plus explicit family = "aft" for
   the plain two-column path (a bare 2-column matrix is otherwise
   ambiguous). "aft" is the only survival family name in v1; the
   discrete-time follow-up will not need a new family (it reuses probit /
   logistic on the person-period expansion).
   [Review 2026-07-10: auto-dispatch fires only from family "auto" (or an
   explicit "aft"); a Surv response with a conflicting explicit family
   (e.g. "gaussian") errors instead of being silently overridden. A factor
   status (survival codes it as type "mright") is detected with a hint;
   status must be 0/1 or logical. A Surv-like object with no type
   attribute is documented as treated right-censored.]

3. Prediction semantics. The recorded fit is f(x) + offset = E[log T | x]
   on the LOG scale (storeSample de-scales through fitScale/fitShift like
   Gaussian). RECOMMEND the minimal v1 surface:
   - predict()/extract() return the linear predictor E[log T | x] (log
     scale) by default -- posterior samples of f, no probability
     transform (unlike the binary families). This is the raw quantity
     every downstream survival summary is built from.
   - median survival time is exp(E[log T | x]) (the log-normal median);
     document it as a one-line transform of the returned lp rather than a
     separate predict type in v1.
   - survival curve S(t | x) = 1 - Phi((log t - E[log T | x]) / sigma)
     needs the per-draw sigma, already exposed (fit$sigma / getSigmas).
     RECOMMEND shipping one small R helper that combines the lp draws
     with the sigma draws over a user time grid, rather than overloading
     predict with a type argument in v1. Keep it out of the engine.
   The `family` element on the packaged fit records "aft" so
   predict/extract skip the probability transform (probabilityFromLatents
   is probit/logistic only).
   [Review 2026-07-10: the helper shipped as the S3 generic
   survivalProbabilities (bart method; an rbart method errors - grouped
   AFT is unreachable until rbart_vi grows a family argument, and a curve
   dropping the intercepts would be wrong). It returns DRAWS (draws x
   times x observations; a chain margin under combineChains = FALSE) per
   the package's extract-draws / fitted-mean / ci.level-interval tiers,
   computed in the uncombined convention where sigma aligns with the fit
   draws. The log-scale predict default stands (riAFTBART consumes
   yhat.train.mean on the log scale); the Rd states the scale plainly and
   documents the response/link type aliases as inert for aft.]

4. setResponse under censoring. RECOMMEND: setResponse(newLogTimes)
   replaces the observed log-times, keeps the status fixed at creation
   (which observations are censored is structural, like the group
   assignment), and redraws the censored latents against the current fit
   (the probit pattern; delegated through the contained Gaussian's
   setResponse). Changing the censoring pattern needs a new sampler.

5. Weights. RECOMMEND unsupported in AFT v1 (reject at creation, as
   probit does): a weighted truncated-latent draw is not a coherent
   likelihood at this scale, and riAFTBART does not need it. The
   uncensored path could carry Gaussian weights, but mixing weighted and
   truncated draws is deferred rather than half-specified.

## riAFTBART composition check

riAFTBART's model is AFT with random intercepts per cluster; its outer
loop alternates a BART AFT fit with a Gibbs draw of the cluster
intercepts. That is exactly GroupedResponse decorating AFTResponse:
GroupedResponse::refreshLatents calls base_->refreshLatents against
f + b (so censored latents are drawn against the group-adjusted fit),
draws the intercepts against the fresh working response with null
weights (drawGroupEffects handles null = unit weights), then drawSigma
delegates to the base against f + b. The composition is generic over the
base family (it already carries probit and Gaussian), so AFT slots in
with no GroupedResponse change, replacing riAFTBART's R outer loop like
rbart_vi's. Confirmed by construction; gated by a grouped-composition
test in stage 3.

## Out of scope (v1)

Nonparametric / heteroscedastic baseline hazards (BART-package-style),
competing risks, and left/interval censoring. The discrete-time hazard
family is a separate follow-up section, not this landing.
