# Survival models: design

Status: AFT log-normal LANDED 2026-07-10 (f0efc03; grouped rbart_vi support
ac6ec2c). Discrete-time hazard PROPOSED 2026-07-18 (the section below).
First survival families for 1.0-x. Companion to the roadmap in
docs/plans/survival-models.md and the extensions table in
docs/design/core-generalization.md ("Survival (AFT, discrete-time hazard),
ordinal, quantile: ResponseModel latents; person-period expansion at
ingest").

The AFT model is documented first; the discrete-time hazard section follows
it below. Discrete-time hazard is person-period ingestion sugar over the
existing binary families - it adds no engine code - proposed here after AFT
landed.

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
family is the separate section below (proposed 2026-07-18), not part of
this AFT landing.

# Discrete-time hazard (proposed 2026-07-18)

RESOLVED (VD 2026-07-19): all three forks as re-grounded by the survey.
Time grid: default = the distinct observed times (the surv.bart and
pammtools convention), breaks coarsens (an explicit boundary vector or
an integer quantile count), the expansion-size guard is the refusal
point. Naming: the two tokens, family = "hazard" (probit) and
"hazard.logistic". Default link: probit - the house binary default and
surv.bart's own default (type = "pbart") align.

The plan's second survival family, and the lightest possible one: it adds
no ENGINE code - no bridge, ResponseModel, or C++ change of any kind. A
discrete-time hazard fit IS an ordinary probit or logistic fit on a
person-period-expanded design. The R surface, though, gains real code
beyond the expander, and this section says so plainly: a family REMAP (the
hazard token resolves to its binary link token before any family-keyed
switch runs, section 2), a hazard MARKER on the packaged fit (reporting
dispatch, section 4), and the survival-curve reporting branch. The data
layer is untouched (docs/design/data-store.md): the expanded rows are
ordinary rows and the period index is an ordinary ordinal predictor, so
ColumnStore, the cut grid, the bridge, and the binary ResponseModels all
run bit-for-bit as they do on any binary fit. Contrast aft, a REAL engine
family: it adds ResponseFamily::aft, threads a status vector to C++, and
takes its own arm in every family-keyed switch (node.scale, sigma
handling, weights, log-likelihood); hazard adds NONE of that, because
after the remap the engine-facing family IS probit or logistic and there
is nothing left to switch on. That is why this section is short - most
decisions below reuse what AFT and the binary families already ship, and
the genuine forks (the time grid, section 1; naming, section 2; the
default link, section 3) are surface taste, not sampler design.

## 1. Model

K ordered discrete periods 1..K. The discrete hazard is the conditional
probability of failing in period t given survival through t-1,

    h(t | x) = P(T = t | T >= t, x) = g(f(x, t) + o),

for a binary link g (probit or logistic, section 3), the sum-of-trees fit
f, and an offset o on the LINK scale (log-odds / probit, NOT the log-time
scale AFT's offset lives on). Survival and the discrete density follow by
the usual product,

    S(t | x) = prod_{k<=t} (1 - h(k | x)),   P(T = t | x) = h(t) S(t-1).

Person-period expansion is the standard device (Cox 1972; Allison 1982;
Singer-Willett 2003) that turns this into a binary regression. A subject i
observed to period t_i (event or censoring) contributes t_i rows, one per
at-risk period k = 1..t_i, each carrying the subject's covariates x_i, the
period index k, and a binary outcome y_ik = status_i * 1{k = t_i} (1 only
in the terminal period of an event; a censored subject's rows are all
zero). The per-subject likelihood

    L_i = h(t_i) prod_{k<t_i} (1 - h(k))    (event, status_i = 1)
        = prod_{k<=t_i} (1 - h(k))          (censored, status_i = 0)

is in BOTH cases exactly the Bernoulli product over the person-period rows,
prod_{k=1}^{t_i} h_ik^{y_ik} (1 - h_ik)^{1 - y_ik}, so posterior inference
on the expanded binary data IS inference on the discrete hazard. The
same-period convention, stated explicitly: a subject whose event falls in
period t and a subject censored in period t are BOTH at risk in period t
and both contribute a row for it - the event's terminal row contributes
h(t) (y = 1), the censored subject's terminal row contributes 1 - h(t)
(y = 0). That is, censoring within a period is read as surviving the
period (censoring at interval end, the Singer-Willett convention); the
expansion above encodes it by construction, since a censored subject's
rows are all zero. No censoring parameter and no latent survival time
survive the expansion: right-censoring is encoded entirely by a subject's
rows ending without a terminal 1. (Contrast AFT, where censoring drives a truncated-latent redraw
in the engine; here it is pure data shape, which is why discrete-time
hazard needs no engine support.)

**How time enters (settled: an ordinal period column).** The period index k
is appended to x as one ORDINAL predictor column (integer 1..K) the trees
split on with ordinary threshold cuts (period <= c). This is the whole
mechanism by which the baseline hazard is learned: a tree splitting on the
period column and nowhere else reproduces a step-function baseline; a tree
splitting on period AND a covariate represents a time-varying
(non-proportional) covariate effect for free. The ordinal type is
load-bearing - a threshold cut groups CONTIGUOUS periods, so adjacent
periods share leaves and the leaf prior borrows strength across time,
giving a smooth baseline rather than K independent hazards. Leaf-prior
calibration is inherited from the chosen binary family (node.scale 3.0
probit / pi*sqrt(3) logistic, R/dbarts.R:502/509); nothing new is
calibrated.

The classical model instead uses a saturated set of K per-period intercepts
(a free baseline parameter per period; the cloglog/logit discrete-hazard
GLM). BART does not fit per-period intercepts - there is no such parameter
block - and should not: exchangeable per-period effects throw away the
ordering the ordinal column exploits, and the tree-structured baseline is
the more parsimonious, strength-borrowing object BART exists to provide. A
dummy-indicator expansion of the period (K-1 binary columns) is a third
option, dominated: it lets a LINEAR model represent free per-period effects
but for BART just fragments one clean ordinal column into K sparse ones with
no benefit, since the trees already represent any function of period from
the single column. Settled: one ordinal period column.

**DECISION (open for VD) - ties and the time grid.** The model lives on a
discrete grid, so continuous event times must be placed on one before
expansion; ties WITHIN a period are automatic (same period index = same
time-column value = same risk set {i : t_i >= k}). The open question is
who defines the grid and how. The survey below is decisive here and MOVED
this fork's recommendation: an earlier draft refused continuous times
outright unless `breaks` was given, a stance no surveyed package shares.

- Recommend (revised per the survey): the default grid is the sorted
  DISTINCT OBSERVED TIMES - exactly surv.bart's default (events <-
  unique(sort(times)), all observed times, censoring times included) and
  pammtools' (cut = NULL -> all unique event times). A `breaks` argument
  coarsens on request: a numeric vector gives explicit interval
  boundaries (the discSurv contToDisc convention), a single integer K
  bins at the quantiles (1:K)/K of the observed times with times snapped
  to the grid (the surv.bart K convention, verbatim). Coarsening is thus
  never silent - the default grid is the data's own resolution - and the
  refusal point is the N' row guard (section 2): heavily continuous times
  on the default grid trip the cap with an actionable message naming
  breaks/K, so the pathology the earlier refusal targeted is still
  caught, while a BART user's surv.bart-formed expectation (the grid is
  my observed times) is honored.
- Alternative (the earlier draft): refuse continuous times without
  explicit breaks. Maximally explicit about the fact that the grid IS
  the model, but the survey shows no package behaves this way -
  surv.bart and pammtools default to the unique-times grid, and the one
  explicit-limits requirement (discSurv contToDisc) lives in a dedicated
  binning utility, not at the model entry. Refusing what the direct BART
  precedent accepts trades user expectation for a purity the N' guard
  already delivers.
- Silent auto-coarsening (deciles or a fixed K when times look
  continuous) stays REJECTED: no surveyed package coarsens without being
  asked either, and it hides a modeling choice - the number and placement
  of periods change the model, not just its resolution.

Strongest counter to the revised recommendation: the unique-times default
makes the model depend silently on the data's time RESOLUTION (the same
phenomenon measured finer gets a different grid), where an
explicit-breaks requirement forces the user to own that choice; and the
N' guard catches only the memory pathology, not a
statistically-too-fine grid that fits but mixes badly. VD to confirm.

## Survey: surv.bart and the expansion ecosystem (added 2026-07-18)

What a user arriving at a discrete-time-hazard BART would expect, per
VD's ask, surveyed before the forks were settled. Package behavior
verified against current documentation and source on 2026-07-18 except
where flagged as memory.

**BART::surv.bart (Sparapani-Logan-McCulloch; BART 2.9.10, 2026-01-24) -
the direct precedent.** The one shipped discrete-time survival BART, and
the strongest user-expectation signal for BART users specifically:

- Ingestion: (times, delta) vectors beside x.train; the person-period
  expansion happens INTERNALLY, in surv.pre.bart - exactly this
  section's ingestion-sugar shape. surv.pre.bart emits y.train (the
  binary indicators) and tx.train with the time PREPENDED as the first
  column of the covariates (X.train[k, ] <- c(events[j], x.train[i, ]),
  source-verified) - the same single-time-column device as section 1,
  modulo column position.
- Grid: default events <- unique(sort(times)) - ALL distinct observed
  times, censoring times included (source-verified). An optional K
  coarsens to the quantiles (1:K)/K (events <- unique(quantile(times,
  probs = (1:K)/K))) with times snapped to the grid; an optional
  `events` argument supplies an explicit grid. As-given by default,
  coarsened on request, never a refusal.
- Link: type = "pbart" (PROBIT) is the DEFAULT; "lbart" (logit) is the
  alternative. Link selection by ARGUMENT on one entry point, not by
  separate families.
- Returns: surv.test / surv.test.mean are SURVIVAL probabilities
  S(t | x) cumulated from the per-period draws, with yhat.* the
  f(t, x) latent draws alongside - the package-level deliverable is the
  survival curve plus the latent channel, section 4's shape
  independently arrived at.
- Siblings in the same package (man index verified): crisk.bart /
  crisk2.bart + crisk.pre.bart (competing risks, two causes) and
  recur.bart / recur.pre.bart (recurrent events) - so the section 6
  competing-risks door has a shipped design to imitate when it opens.

**Expansion utilities.** discSurv (dataLong + contToDisc): dataLong
expands ALREADY-DISCRETE times; binning continuous times is a separate,
explicit step - contToDisc(dataShort, timeColumn, intervalLimits), whose
intervalLimits is REQUIRED (explicit right borders, or with equi = TRUE
a count of equidistant intervals), intervals right-bounded
[a_{k-1}, a_k). pammtools (as_ped): the piecewise-exponential
(Poisson, log-duration offset) cousin of the cloglog tradition;
cut = NULL defaults to "all unique event times", with max_time for
administrative truncation, and ped objects carry newdata methods that
expand prediction data onto the training intervals. casebase fits
smooth-in-time hazards by case-base sampling (logistic regression over
sampled person-moments) - a different device with no person-period
grid, noted for completeness (memory, not re-verified). Net finding: NO
surveyed package refuses continuous times, and none coarsens silently;
the expansion ecosystems default to the data's own unique times and
coarsen only on request.

**Link expectations by tradition.** The applied discrete-time-hazard
literature is logit-dominant - Allison (1982) and Singer-Willett (2003),
the canonical textbook treatment, are logistic throughout
(standard-knowledge claim, not re-verified this pass). The
grouped-continuous-PH tradition expects CLOGLOG (Prentice-Gloeckler
1978; the mgcv/brms person-period recipes use binomial(cloglog), and
pammtools' Poisson PED is the same limit). Probit's precedent in this
exact space is surv.bart itself: the default link of the only shipped
discrete-time BART. Each candidate default has a real constituency -
logit (applied survival analysts), cloglog (grouped-PH), probit (BART
users via surv.bart, plus the dbarts house default).

## 2. Ingestion surface

**What the user passes (reuse AFT's response shapes).** The response is the
SAME (time, status) AFT takes - a survival::Surv object or a two-column
(time, status) matrix / data frame - parsed by the existing
extractSurvivalResponse (R/dbarts.R:17), or a thin sibling of it that skips
the log() transform (R/dbarts.R:56) and returns the raw times and status for
expansion. Everything extractSurvivalResponse already enforces carries over
(right-censoring only, status 0/1 or logical, the factor-status "mright"
hint, positive times). A Surv response cannot AUTO-dispatch to the hazard
model - family = "auto" with a Surv already selects aft (the auto-dispatch
at R/dbarts.R:252-266) - so the discrete-time model is always requested
explicitly by family, and the Surv conflict guard must admit the hazard
tokens (the remap block below).

**The remap (settled; required whatever naming wins).** A hazard token
CANNOT flow through the family-keyed resolution unchanged - verified
against the code: node.scale is a switch(family, ...) with NO default
(R/dbarts.R:498-510), so an unknown token yields NULL; control@binary
keys on family %in% c("probit", "logistic") (:359), so the binary
machinery would stay off; fixedUnitScale excludes it (:368-370), so sigma
would be ESTIMATED and the 0/1 response fit as gaussian; and the weight
policy keys on the literal tokens (:382-403). The design is therefore an
early REMAP: immediately after the ingestion resolves the hazard request
and the expander runs, the R surface rewrites family to the underlying
binary token ("probit" or "logistic"), BEFORE the node.scale /
control@binary / fixedUnitScale / weight-policy resolution - so every one
of those switches, the model object (model@family records the binary
token), the bridge, and the engine see an ordinary binary fit and are
untouched. The hazard provenance survives R-side only: the ingestion
parks the period grid where the wrapper packaging can reach it (the
bartcore.survival -> $status precedent, R/bart.R:247-251), and the
packaged fit carries it as the hazard marker (section 4). aft is the
contrast: a real engine family with its own arm in each of those switches
and a status vector threaded to C++; hazard takes an arm NOWHERE because
after the remap there is nothing left to switch on.

The Surv conflict guard extends to admit the hazard tokens: today
responseIsSurv errors for any explicit family outside "auto"/"aft"
(R/dbarts.R:229-235), so Surv + a hazard family would be refused as a
conflict. The guard's whitelist gains the hazard tokens - a Surv response
declares survival, and the family then picks WHICH survival model
(family = "auto" keeps selecting aft, the landed default; a hazard token
selects discrete-time). Every naming form below needs this same guard
edit, so it differentiates none of them.

**DECISION (open for VD; entangled with section 3) - naming.** probit and
logistic are already two family tokens for one Bernoulli model under two
links (R/dbarts.R:198-199); the hazard model, being sugar over exactly
those, can be named the same way:

- Recommend: two tokens, family = "hazard" (probit link, the house binary
  default) and family = "hazard.logistic" (logit link), with
  "hazard.probit" an accepted alias. No new argument is introduced (dbarts
  has no link argument; ordinal and nbinom added tokens, not arguments -
  docs/design/ordinal.md, negative-binomial.md), the form is symmetric
  with the existing probit/logistic split, and a future binary family
  (cloglog, section 3) yields "hazard.cloglog" for a remap-table entry.
- Alternative A: a single family = "hazard" plus a new
  `hazard.link = c("probit", "logistic")` argument. Tidier as one token,
  but a link argument is a new surface primitive for one family's
  benefit. The survey STRENGTHENS this alternative: surv.bart selects its
  link by exactly this shape (one entry point, type = "pbart"/"lbart"),
  so the direct BART precedent is argument-shaped, not token-shaped.
- Alternative B: an orthogonal `discrete.time = TRUE` flag with the link
  carried by the ordinary family = "probit"/"logistic". Cleanly separates
  expansion from link and composes with any future binary family
  automatically, at the cost of a new flag argument whose meaning is
  entangled with the response shape (discrete.time = TRUE on a
  non-survival response must be refused). An earlier draft rejected B for
  a claimed Surv-ambiguity against aft; that argument does not withstand
  scrutiny - ALL three forms need the identical guard edit at :229, so
  the collision differentiates nothing, and B stands or falls on the
  argument-count trade alone.

Recommend the two-token form, HELD after the survey on internal-
consistency grounds: dbarts's families are tokens (aft, ordinal, nbinom
all extended the family vector, never an argument or a dedicated entry
point), a dbarts user's expectation is set by dbarts's own surface, and
surv.bart's other shape signal - a DEDICATED FUNCTION (surv.bart beside
pbart/lbart) - is a non-option here for the same reason. Merits: no new
argument, family-vector symmetry with probit/logistic, one token to grep
for. Strongest counter, upgraded by the survey: the one shipped
discrete-time BART picks its link by argument (type = "pbart"/"lbart"),
so a surv.bart migrant would find Alternative A the familiar shape. The
tokens are added to the dbarts and bart2 family vectors
(R/dbarts.R:195-203, R/bart.R:378-386) and remapped as above; xbart and
rbart_vi omit them (their match.arg vectors, R/xbart.R:27, R/rbart.R:48,
ARE the refusal, the ordinal/nbinom precedent) - grouped hazard is
section 6.

**Where expansion happens (settled: an R helper at ingest).** The
person-period expander is a pure R transform invoked during ingestion,
exactly where extractSurvivalResponse is (R/dbarts.R:252-266), BEFORE
dbartsData builds the model matrix. It consumes (x, time, status, grid) and
emits an ordinary (X', y') binary problem: X' is x replicated down each
subject's at-risk periods with the ordinal period column appended, y' the
per-period indicators. With the remap done (family now reads "probit" or
"logistic"), that design flows through the dbarts -> dbartsData -> bridge
-> binary-family path with no change downstream of the remap point. The
engine, the bridge (resolveFamily, applyGroupAttribute,
applySurvivalAttribute - src/R_interface_bartcore.cpp:1110/1393/1449), and
the ResponseModels never learn the rows are person-periods; there is no
bartcore.hazard attribute and no status vector to C++ (unlike aft,
R/dbarts.R:530-537) - the censoring is already baked into y'.

**The memory story (be honest).** The expanded design is
N' = sum_i t_i = n * mean(t_i) rows by (p + 1) columns, and x is REPLICATED
mean(t_i)-fold down the rows. For n = 1000 over a mean 20 periods this is
20000 rows - trivial; for n = 1e5 over a mean 100 periods it is 1e7 rows
with a 100-fold predictor blowup - a real memory event. Practical guard:
refuse (or require an explicit opt-in) when N' exceeds a cap, the message
naming the two coarsening levers (`breaks` boundaries or integer K,
section 1). Under the revised grid fork this cap IS the refusal point for
over-fine grids: continuous times run on the unique-times default until
they trip it, and the error is actionable rather than a blanket refusal.
Genuinely discrete data with a modest K never sees the guard.

**Offsets and weights under expansion (settled: replicate per subject,
inherit the binary policy).** A per-subject offset o_i is a link-scale
hazard shift replicated to all of subject i's rows; a per-subject weight
likewise replicates and then follows the CHOSEN binary family's existing
policy unchanged - probit refuses non-unit weights (R/dbarts.R:382-393),
logistic requires positive integer counts (:394-403). No new offset or
weight semantics are introduced; they are whatever the underlying binary fit
already does, applied to the replicated rows. A time-VARYING offset (a
distinct o_ik per period row) is the natural person-period generalization,
tied to the time-varying-covariate door (section 6).

**Reused vs new.** Reused wholesale: the (time, status)/Surv parsing, the
probit/logistic families and their entire engine/bridge/prediction stack
(the prediction generics stay correct BECAUSE of the marker design,
section 4), the matrix interface, dbartsData, and the cut grid (the period
column is an ordinary ordinal column, and a test subject's period column
bins identically to training by the shared-grid construction,
data-store.md). New, all R-side: the expander, the family token(s) and
their remap, the Surv-guard extension, the `breaks`/grid argument, the N'
guard, the packaged $periods marker, and the survival-curve reporting
branch (section 4).

## 3. Link (probit vs logistic; cloglog)

**DECISION (open for VD) - default link.** Both binary links give a valid,
standard discrete-time hazard model; the choice is convention, not
correctness, and the reduction gate (section 5) holds for either.

- Recommend probit as the default (family = "hazard" = probit link),
  UPHELD and strengthened by the survey: surv.bart's default is
  type = "pbart" - the PROBIT link - so the only shipped discrete-time
  BART defaults to probit, and a BART user arriving from it gets exactly
  what they expect; probit is also dbarts's house binary default
  (auto -> probit, R/dbarts.R:350; node.scale 3.0, :502), aligning the
  family with the package's other fixed-unit-scale defaults.
- Counter: the applied discrete-time-hazard literature is logit-dominant
  (Allison 1982; Singer-Willett 2003, logistic throughout - the survey's
  link-traditions finding), and the grouped-PH tradition expects cloglog
  (unavailable in v1, section below) - so a survival analyst from OUTSIDE
  the BART world expects logit (or cloglog), not probit; logistic's
  node.scale is already provisioned (pi*sqrt(3), :509).

Recommend probit: the two constituencies split (BART users -> probit,
applied survival -> logit), the tie-breaks are the house default and the
direct BART precedent, and logit stays one token away. VD may still
prefer logit to court the applied-survival audience. Low stakes - the two
differ only in the link, and the section-5 gate certifies whichever wins.

**cloglog (settled: not needed for v1, recorded as a door).** The
complementary-log-log link is the CLASSICAL discrete-hazard link, for a
specific reason: if continuous T follows a proportional-hazards (Cox) model,
the grouped/interval discrete hazard satisfies cloglog(h_k) = alpha_k +
beta'x exactly (Prentice-Gloeckler 1978), so cloglog is the discrete link
that preserves continuous-time hazard ratios. That motivation LARGELY
DISSOLVES under BART: f(x, t) is a fully flexible, non-proportional surface
(the period column interacts with covariates freely through the trees), so
there is no proportional coefficient to preserve and no hazard-ratio
interpretation to protect. logit (or probit) is the standard, defensible
substitute in exactly this nonparametric setting. The engine has no cloglog
family, and adding one is a genuine ResponseModel item (a new binary link
with its own latent/augmentation), NOT ingestion sugar - so it is out of v1.
The door is clean and cheap once wanted: because discrete-time hazard is
sugar over whatever binary families exist, a future cloglog family becomes
family = "hazard.cloglog" for one remap-table entry plus its link's
probability transform (section 2's naming recommendation is what keeps it
this cheap).

## 4. Prediction and reporting

**The marker, and why the link transforms stay correct (settled).** The
packaged fit's $family element records the BINARY token ("probit" or
"logistic") - the engine truth, and what the packaging already does for
free, since the family element is read off fit$model@family
(R/bart.R:218/:230), which the remap set. The hazard provenance rides a
separate marker: $periods, the ordered period grid, present only on
hazard fits (the aft $status precedent, R/bart.R:247-251, and the
bartOrdinal $levels precedent). This split is load-bearing, not
bookkeeping: every generic that keys behavior on $family then does the
right thing with ZERO modification. probabilityFromLatents selects its
inverse link solely by identical(family, "logistic") and defaults to
pnorm (R/generics.R:13-19) - had $family recorded a hazard token instead,
a logit-link hazard fit would take the probit transform SILENTLY, a wrong
number with no error - and pointwiseLogLikelihood dispatches on $family
and errors on families it does not know (R/generics.R:53; its
probit/logistic Bernoulli branch is exactly the per-row discrete-hazard
log-likelihood, by the section-1 identity). The rejected alternative -
recording a hazard token in $family and teaching each consumer the new
tokens' links - is real new code in every $family consumer, and any
missed one fails silently in the transform rather than loudly. So the
rule is: LINK dispatch keys on $family (the binary token); HAZARD
dispatch (survivalProbabilities, any hazard-aware printing) keys on the
$periods marker, never on $family. No new S3 class: a plain bart object
with $periods.

**Hazards are the native output; the survival curve is the new code.**
Because the fit is a binary fit on the expanded rows, the ordinary
predict/extract channels return per-(subject, period) hazard draws
h_ik = g(f(x_i, k) + o) unchanged - type = "ev"/"response" the hazard
probability (through the correct link, per the marker design above),
type = "bart"/"link" the latent f (the binary convention, R/generics.R).
The new reporting step is the survival curve
S(t | x) = prod_{k<=t} (1 - h(k | x)), a cumulative product of
(1 - hazard) over the periods up to t.

**survivalProbabilities: same entry point and shape, DIFFERENT evaluation
path than aft's.** The generic (R/generics.R:7) gains a hazard branch in
survivalProbabilities.bart, keyed on the $periods marker beside the aft
$family gate (R/bart.R:1779), computing S by cumulative hazard products
instead of the log-normal tail. The SHAPE convention matches aft's
exactly - draws x times x observations, a chain margin under
combineChains = FALSE, the extract-draws / fitted-mean / ci.level tiers -
so the two survival families share ergonomics. The EVALUATION cannot
mirror aft's training path, though: aft's newdata = NULL path reads one
stored linear predictor per subject (extract type = "bart",
sample = "train", R/bart.R:1793-1794) because every subject has exactly
one row, but the hazard training design is RAGGED - subject i has only
t_i rows, so its stored training fits stop at its own event/censoring
period and cannot supply h(k | x_i) beyond it. A full-horizon S(t | x_i)
needs hazards at EVERY requested period, and only tree replay provides
them: the hazard method ALWAYS re-expands its subjects to the requested
grid and predicts - training data included - and therefore requires
keepTrees UNCONDITIONALLY (the predict.bart guard, R/generics.R:153),
where aft's training-data path never needs it. Requested `times` are grid
periods; S cumulates (1 - h) through each. Both quantities are thus
available: per-period hazards via the ordinary binary predict/extract,
survival curves via survivalProbabilities.

**Test data expands too (settled).** A test subject has no observed event
time, so it expands to the FULL grid - one row per requested period - the
period column binned against the shared training cut grid (identical bins
by construction, data-store.md). survivalProbabilities' `times` argument
sets the horizon; the default is the training grid's periods.

## 5. The reduction gate (primary correctness gate)

**Bitwise equality with the by-hand binary fit.** The single load-bearing
gate, and it is free: a hazard fit must equal, draw for draw, the ordinary
probit/logistic fit a user would get by expanding the data by hand and
calling dbarts(family = "probit"/"logistic") on it with the same seed and
control. Since the family adds no engine code and remaps to the binary
token before the sampler is built, this equality is exactly the claim
"hazard is sugar" made testable: with the same seed and the same expanded
design (same column order - period column appended last - same cut grid,
same offset/weights), the two fits consume the identical RNG stream. The
equality is over the DRAW components - expect_identical on the trees,
latents, yhat.train/yhat.test, and varcount channels - NOT on the packaged
objects, which differ by construction in the hazard-only fields ($periods,
and $family reads the same binary token on both, section 4); the test
compares the draw channels and asserts the marker fields are the ONLY
difference. This is the discrete-time analog of AFT's uncensored ==
gaussian reduction (the "Engine: AFTResponse" section above), and like it
certifies the feature by tying it to already-gated code.

**No engine gate is needed.** The engine, bridge, and ResponseModels are
untouched, so tests/cpp does not change and no exact-posterior gate is
added: the binary exact-posterior gate
(benchmarks/R/logistic-reference.R, which certifies both probit and logistic
to the exact posterior) ALREADY certifies the sampler's correctness on
binary data, and the reduction gate proves the expander feeds the sampler
exactly that. The new correctness surface is all R-side: the expander and
the remap (covered bitwise by the reduction gate) and the marker-dispatched
reporting (covered by the surface tests below).

**The usual surrounding gates.** (i) A new equivalence scenario in
benchmarks/R/equivalence.R (a small hazard fit) recording the per-period
hazard and survival draws; existing anchors stay IDENTICAL - hazard adds no
draw to any existing family's stream, it just drives the probit/logistic
path on expanded data - so the neutrality trail is a re-record that only
ADDS the scenario, the ordinal/nbinom precedent (docs/design/ordinal.md
section 7, negative-binomial.md section 6). (ii) tinytest surface tests:
family routing and refusals (including Surv + hazard admitted, Surv + a
non-survival family still refused), the breaks/grid argument in both
forms (boundary vector and integer K; the grid defaulting to the
distinct observed times), the N' guard refusal, offset/weight
replication, a
logit-hazard fit's probability channel equal to plogis of its latent
(locking the marker design's link dispatch, section 4), and the
survivalProbabilities shape (draws x times x obs) and monotonicity (S in
[0, 1], non-increasing in t). (iii) A recovery test: simulate discrete-time
survival data from a known baseline hazard and covariate effect, fit, and
check the recovered hazards / S(t | x) against truth - the family-level
smoke beyond the reduction gate.

## 6. Scope and doors

- **Time-varying covariates (OUT of the auto-expander; IN by construction).**
  The person-period form makes them natural: each (subject, period) row can
  carry period-specific covariates rather than the replicated baseline x_i.
  The v1 AUTO-expander assumes time-constant covariates (it replicates x_i),
  so time-varying covariates are out of the sugar path - but they are
  already fully supported by the MANUAL path, because a user who builds the
  person-period design themselves (time-varying columns plus a period column)
  and fits family = "probit"/"logistic" gets exactly the model, and that
  manual fit is precisely the reduction-gate target (section 5). The escape
  hatch is a first-class, gated use, not a missing feature. Door: a
  long-format ingestion mode accepting pre-expanded (subject, period, x, y)
  is the natural convenience follow-up.

- **Competing risks (OUT; door).** J > 1 event types make the per-period
  outcome MULTINOMIAL (survive, or fail from cause 1..J), i.e. discrete-time
  hazard over the multinomial forest (docs/design/multinomial.md), not a
  binary family. Out of v1; once multinomial ships, competing-risks discrete
  hazard is the same expansion with a multinomial outcome, and BART's
  crisk.bart / crisk2.bart (the survey) are the shipped design to compare
  against. Recorded.

- **Left truncation / delayed entry (OUT of v1; the cheapest door).** A
  subject entering the risk set at period e_i > 1 is handled by STARTING its
  expansion at e_i (omitting the pre-entry rows) - a one-argument extension
  of the expander (an entry-time input), strictly cheaper than the others.
  Out of the minimal v1 surface, but flagged as near-term and low-cost.

- **Grouped / frailty (OUT of v1 surface; composition already exists).** A
  random-intercept (shared-frailty) discrete-time hazard is GroupedResponse
  decorating the binary family on the expanded data - and that decorator is
  generic over its base family and already ships and is gated for probit
  (grouped-random-effects.md; the AFT note's riAFTBART composition check
  above, "generic over the base family ... already carries probit and
  Gaussian"), with grouped logistic feasible on the same seam (logistic has
  a Gaussian working response). So grouped hazard is reachable the moment the
  expander feeds rbart_vi - whose family vector (R/rbart.R:48) would grow the
  hazard token as it grew "aft". The only work is surface, exactly
  paralleling the AFT grouped follow-up
  (docs/plans/survival-grouped-surface.md). v1 ships the ungrouped hazard on
  dbarts/bart2; grouped hazard (rbart_vi + hazard expansion +
  survivalProbabilities.rbart cumulating hazards) is a recorded follow-up,
  feasible by construction.

- **cloglog link (OUT; door).** Section 3: not needed for a valid v1, and a
  future cloglog binary family becomes a hazard link for free.

- **Continuous-time models (OUT; that is AFT).** Genuinely continuous
  survival is AFT's territory (shipped) and the future nonparametric
  baseline-hazard / Cox-type families in the AFT out-of-scope list. Discrete-
  time hazard is deliberately the discrete-grid model; a user with continuous
  times either bins them (section 1) or uses family = "aft".
