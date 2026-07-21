# Hurdle / two-part model: design

Status: LANDED 2026-07-20 (901581e..7903855; records docs/plans/hurdle.md
Landing). Plan: docs/plans/hurdle.md. A hurdle model fits a zero-inflated / semicontinuous outcome by factoring
its likelihood into two conditionally-independent parts: an OCCUPANCY model of
1{y > 0} over all n observations, and a POSITIVE-PART model of y restricted to the
subset {i : y_i > 0}. The load-bearing finding (section 0): because the parts share
no parameters and the split is OBSERVED (not latent), the two forests are
conditionally independent given the data - a hurdle fit is two ordinary dbarts fits
glued at prediction time, not a coupled multi-forest engine model. This inverts the
framing forest-combiner.md and multi-forest-models.md carried (hurdle as the
model that breaks Chain's single response_) and lands hurdle R-side, the discrete-
time hazard precedent (docs/design/survival.md section 5), NOT the heteroscedastic
engine precedent. rng class: BYTE-neutral in the engine unconditionally (no C++
code; a hurdle fit drives the existing binary and gaussian families verbatim).

Landscape (general-knowledge survey, not re-verified this pass): the semicontinuous
two-part model is Duan-Manning-Morris-Newhouse (1983, JBES 1:115) - a logistic
any-expenditure part plus a lognormal (or smearing-retransformed) positive-cost
part - the canonical healthcare-cost / insurance model. The count hurdle is Mullahy
(1986, J.Econometrics 33:341) and Cragg (1971): binary at-zero plus a
zero-TRUNCATED count. Distinct from ZERO-INFLATION (Lambert 1992), which mixes a
structural-zero mass with a count that can also emit zeros - a LATENT split that
does couple (section 0). No maintained BART package ships a hurdle/two-part family:
BART, bartMachine, and stochtree have none; Murray (2021, Log-Linear BART) fits
zero-inflated counts by a coupled augmentation, not a hurdle. negative-binomial.md
recorded the count-hurdle door. So there is a real gap and no shipped BART
convention to imitate - the choice is grounded in the applied two-part / hurdle
literature, not in a BART precedent.

## 0. The load-bearing finding: the parts are CONDITIONALLY INDEPENDENT

A hurdle factors the likelihood of a non-negative outcome with exact zeros as

    p(y_i) = p(y_i > 0 | x_i) * p(y_i | y_i > 0, x_i),

    L = prod_{i: y_i = 0} (1 - pi_i)  x  prod_{i: y_i > 0} pi_i g(y_i | y_i > 0),

with pi_i = P(y_i > 0 | x_i) the occupancy probability and g the positive-part
density. The occupancy parameters (pi, a binary forest over all n) and the positive
parameters (g's forest + its sigma over the y > 0 subset) appear in DISJOINT factors
and share nothing. Under independent priors the joint posterior therefore FACTORS
into the occupancy posterior times the positive-part posterior. This is the central
architectural fact, and it separates hurdle from every landed multi-forest model:

- BCF couples through a mu + b tau: the forests blend into ONE per-observation
  location a shared response scores (combiner.hpp combinedFits/drawGlue).
- Multinomial couples through the softmax: K forests share one likelihood, an
  interleaved PG cycle, and a level-centering move (multinomial.md).
- Heteroscedastic couples through the precision channel: the variance forest sets
  the mean forest's per-observation weight w_i / s^2(x_i) EVERY sweep, and vice
  versa (heteroscedastic.md). That within-sweep coupling is exactly what
  forced the Chain-level integration, the multiplicative roll, and the m'=2 closing
  gate.

Hurdle has NO such coupling. The two forests never read each other within a sweep;
there is no combined location, no glue draw, no interleaving, no shared latent. The
engine's entire multi-forest apparatus - the ForestCombiner surface, the coupling
draws - has nothing to own.

**Hurdle is not zero-inflation.** In a hurdle the split is OBSERVED: every zero
comes from the occupancy part, and the positive part is zero-truncated (emits no
zeros), so there is no latent class-membership to draw. In zero-inflation the zeros
are a MIXTURE of structural and sampling zeros, so each zero carries a latent
indicator drawn each sweep against both components - a genuine coupling that a
factored R composition CANNOT express (it is the coupled cousin, a recorded door,
section 9). negative-binomial.md filed "zero-inflation / hurdle NB" as one item;
they are different models and only the hurdle is conditionally independent. This
note is the hurdle; zero-inflation stays out of scope precisely because it couples.

## 1. The model and its reduction

**Occupancy part (both variants).** A binary regression of z_i = 1{y_i > 0} over
ALL n observations, fit by a shipped binary family: probit (house default) or
logistic. This is an ordinary dbarts binary fit on the indicator - no new machinery.

**Positive part (the variant fork, section 3).** A model of y over the subset
S = {i : y_i > 0}:

- Semicontinuous (v1 recommend): lognormal, i.e. an ordinary GAUSSIAN fit on
  log(y_i) for i in S. E[y_i | y_i > 0] = exp(f(x_i) + o_i + sigma^2 / 2) (the
  lognormal mean; the retransformation, section 6). Nearly free - it is the AFT
  reduction verbatim (uncensored AFT on log-times == gaussian, survival.md):
  R logs the positive values at ingest and hands the engine an ordinary gaussian
  problem over |S| rows.
- Count (door): a zero-TRUNCATED Poisson or negative binomial over S. dbarts has an
  untruncated NB count family (negative-binomial.md) but NO zero-truncated
  likelihood; a truncated count is real new single-forest family work (section 9).

**Combined mean (prediction only).** E[y_i | x_i] = pi_i * E[y_i | y_i > 0, x_i],
evaluated per posterior draw at report time (section 6). The two fits are combined
NOWHERE else - not in any sweep.

## 2. Decision (fork 1, the gating question) - COMPOSE IN R, do not build in the engine

The decision this note turns on. Two routes:

- **Compose in R (recommend).** A wrapper fits two ordinary samplers - an occupancy
  binary fit over all n and a positive-part fit over the y > 0 subset - and combines
  their posterior draws at prediction time. Zero engine code. The discrete-time
  hazard shape exactly (survival.md): a family token, an R ingestion split,
  a packaged fit that reuses every existing generic.
- **Build in the engine.** A Chain-level two-response generalization: two independent
  (response, forest-group, sigma) tuples, the positive forest scoring only the y > 0
  subset (section 4, the rejected alternative).

Honest weighing of what the engine route BUYS, candidate by candidate:

- One sampler object / shared data handle + cut grid. Weak: the two parts operate on
  DIFFERENT observation sets, so their natural cut grids differ anyway (the positive
  part's cuts are quantiles of the subset's x). Forcing a shared grid is a modeling
  choice, expressible R-side by passing explicit cutpoints to both fits - not an
  engine necessity. The handle-view mechanism data-ownership-4-views.md
  names hurdle as a future consumer of is real, but it serves a subset VIEW the
  engine route would need; R composition sidesteps it entirely (two separate stores).
- One RNG stream. No statistical requirement - the parts are independent, so any
  pairing of their draws is a valid joint draw. Worse, an engine stream would force
  us to DEFINE an interleaving order and then build a gate to match it; R composition
  makes the reduction gate free by construction (section 5), because the hurdle IS
  the two fits.
- Combined prediction / reporting / state serialization. Pure packaging convenience,
  which R gives for free in the bartMultinomial / bartNegbin idiom (a few dozen lines,
  R/bart.R:1269,1746). The engine route would instead DOUBLE the state wire format
  (two response latent/sigma blocks, two forest lists) - a format bump, not a saving.
- LinkingTo / dbarts.h reach. None in v1: nbinom, ordinal, hazard, and heteroscedastic
  all ship with zero dbarts.h exposure (negative-binomial.md). No gain.

Against the engine route: it is the WIDEST Chain surgery in the queue - two response_,
two sigma_, and a subset view threaded through the positive forest's per-observation
kernels (section 4) - and it buys nothing statistically, because nothing couples. Every
landed multi-forest engine model earned its Chain work with a within-sweep coupling an
R outer loop cannot express exactly (BCF's glue, multinomial's interleaved PG,
heteroscedastic's precision feedback). Read carefully, the heteroscedastic precedent
argues AGAINST engine hurdle: the thing that made the engine pay there (coupling) is
absent here. Hurdle is to the binary+gaussian families what discrete-time hazard is to
the binary family - a factored construction over shipped families, right-sized R-side.

**Recommendation: compose in R.** It ships the model at a fraction of the cost, keeps
the engine byte-neutral unconditionally (section 7), and gives the correctness gate for
free (section 5). The strongest argument against - that a user expects "a hurdle
sampler," one object, one .Call, the shape BCF/multinomial took - is answered by the
packaged fit class (section 6): the user sees one `bart2(family = "hurdle")` call and
one fit object; only the internals are two samplers. The case where the engine route
becomes justified is real but is a DIFFERENT model: a sample-selection / correlated
two-part (Heckman) model, where the occupancy and amount errors are correlated and the
parts do NOT factor - that coupling is what would finally pay for the Chain
two-response generalization (section 9). It is not this model.

RESOLVED (pending VD): recommend R composition for v1; the engine architecture is
sketched (section 4) so VD can fork on it with the cost asymmetry explicit.

## 3. Decision (fork 2) - which hurdle variant for v1

- **(a) Semicontinuous two-part, lognormal positive part (recommend).** Occupancy
  binary + gaussian on log(y) over S. The positive part is NEARLY FREE - an existing
  family on a log-transformed subset, the AFT reduction (survival.md). It is the
  Duan et al. healthcare-cost model, the canonical two-part model with the largest
  applied audience (costs, insurance claims, semicontinuous rainfall/consumption).
- **(b) Count hurdle, zero-truncated Poisson / NB positive part.** The zero-inflated-
  count competitor. Needs a zero-TRUNCATED count likelihood dbarts does NOT have: the
  shipped NB (negative-binomial.md) fits untruncated NB(r, p) with mass at zero, and
  truncating at y >= 1 changes the normalization and the PG augmentation - real new
  single-forest-family work, not a composition. Poisson is not shipped at all.

**Recommendation: (a), the lognormal semicontinuous two-part.** The positive-part
family question decides it: lognormal-via-gaussian costs nothing, while a
zero-truncated count is a separate engine arc that must land as its own single-forest
family BEFORE any composition can use it. Under the R-composition architecture the
wrapper is family-agnostic on the positive part (it just calls a sampler over S), so
(b) becomes available the moment a truncated-count family ships - v1 need not choose
the harder positive part to keep the door open. Gamma is a third positive-part family
(new family, door). The occupancy link mirrors hazard: probit default (house default,
R/dbarts.R node.scale 3.0), logistic one token away (section 6).

## 4. The engine architecture (the rejected alternative, sketched for the fork)

Presented so VD can fork on the engine route with substance; NOT the recommendation.
Were hurdle built in the engine, the shape would be two independent tuples in one
Chain, and it must be distinguished sharply from the ForestCombiner:

- **NOT a ForestCombiner.** A combiner assumes shared-response coupling: formForestResponse
  forms each forest's residual against a SHARED response, and combinedFits blends the
  forests into ONE per-observation location that ONE response_ scores
  (combiner.hpp; forest-combiner.md). Hurdle has no shared response and no
  combined location - two forests feed two DIFFERENT responses over two DIFFERENT
  observation sets. The combiner's whole virtual surface is inert. Hurdle is a second,
  orthogonal Chain-level axis, not a combiner instance.
- **Reuses the heteroscedastic nullable-second-object precedent as far as it goes.** A
  distinctly-typed, nullable member beside forests_ (the varianceForest_ shape,
  chain.hpp:304,549-556), null off hurdle so the single-forest path is byte-neutral by
  construction (heteroscedastic.md section 8).
- **Where hurdle needs MORE than heteroscedastic.** Heteroscedastic's second forest (i)
  saw ALL n observations, (ii) SHARED the one response_ (it routed into that response's
  weight channel), and (iii) had NO sigma of its own (the variance forest IS the
  variance). Hurdle's positive-part forest (i) sees a DIFFERENT observation set - the
  y > 0 subset, requiring a subset mask or compacted view threaded through every
  per-observation kernel of that forest (a genuinely new plumbing item; the
  data-ownership-4-views.md handle-view is the mechanism); (ii) has its OWN response_
  of a DIFFERENT family (the single-response_ break forest-combiner.md flagged);
  and (iii) has its OWN sigma_ (the single-sigma_ break). So hurdle breaks BOTH Chain
  invariants heteroscedastic kept, plus adds the subset view. This is why it is the
  widest surgery: two response_, two sigma_, a subset-restricted forest, doubled state.
- **Reporting / predict / state.** Reporting channels: occupancy probability pi(x) over
  all n (train + test); positive-part mean E[y | y > 0, x] (the positive forest replays
  at ALL x, though trained on S); combined E[y] = pi * E[y | y > 0]. Predict replays
  both forests' saved trees at newdata and combines. State serializes two forest lists
  plus two responses' latent/sigma blocks - a wire-format widening (a non-BCF glue bump,
  forest-combiner.md, but here it is a second full response block, not glue).
- **dbarts.h:** no change in v1 either way (universal precedent).

The engine route thus costs the queue's largest Chain change for a model with zero
coupling to exploit - the case for R composition, restated structurally.

## 5. The correctness gate: bitwise reduction to two components

**Primary gate (the hazard-reduction precedent, benchmarks/R/hazard-reduction.R).** A
hurdle fit must reduce BITWISE to its two independently-fit components:

- Component 1, occupancy: bart2(family = "probit"/"logistic") on z = 1{y > 0} over
  all n, same seed/control.
- Component 2, positive part: bart2(family = "gaussian") on log(y[S]) over the y > 0
  subset S, same seed/control.
- Channels: component 1's trees / latents / varcount identical to the wrapper's
  occupancy fit; component 2's trees / sigma / varcount identical to the wrapper's
  positive fit; the packaged objects differ ONLY in the hurdle marker fields (the
  markerOnly assertion, hazard-reduction.R:85). The wrapper's two internal fits ARE
  the two standalone fits at their seeds, so the equality is bitwise by construction -
  the "hurdle is two glued fits" thesis made testable, exactly hazard's "hurdle is
  sugar" gate (survival.md).

**Why NO separate joint exact-posterior gate is needed.** Because the parts are
conditionally independent (section 0), the joint posterior IS the product of the two
marginal posteriors, and each marginal is a shipped family already certified to its
exact posterior (the binary logistic-reference.R / probit gate; the gaussian conjugate
path). Certifying each part certifies the whole. This is stronger than heteroscedastic,
where the m'=1 reduction was only a FLOOR and an independent m'=2 quadrature was the
CLOSING gate (heteroscedastic.md section 9) - because there the coupling could be
self-consistently wrong. Hurdle has no joint coupling to get wrong, so the reduction
gate is not a floor but the COMPLETE correctness argument. (If VD picks the engine
route, section 4, the subset-view plumbing DOES introduce a way to get the joint wrong -
the positive forest scoring the wrong rows - and then an exact gate confirming the
positive posterior matches a subset-only quadrature becomes load-bearing, as
heteroscedastic's divisor gate was. R composition needs none of it.)

**Surrounding gates.** (i) A recovery / simulation smoke: simulate semicontinuous data
(a probit occupancy surface and a lognormal amount surface), fit, and check the
recovered pi(x), E[y | y > 0, x], and combined E[y | x] against truth. (ii) A NEW
hurdle scenario in benchmarks/R/equivalence.R recording the two component channels plus
the combined predict; every existing anchor stays IDENTICAL with NO re-record, since
hurdle adds no draw to any engine family's stream (the ordinal/nbinom/hazard neutrality
trail, survival.md). (iii) tinytest surface tests: family routing, the y >= 0
validation and zero handling, the retransformation of E[y | y > 0], draw alignment
across the two fits, and predict on new data.

## 6. The R surface

- **How the user asks.** family = "hurdle" (v1: probit occupancy + lognormal positive
  part), added to the dbarts and bart2 family vectors (R/dbarts.R:349-360, R/bart.R
  around :413). Following the dbarts token convention (families are tokens, not
  arguments - aft, ordinal, nbinom, hazard all extended the vector, survival.md),
  future variants are further tokens: "hurdle.logistic" (logistic occupancy),
  "hurdle.nbinom" (count positive part once a truncated-count family ships). The
  surv.bart-flavored alternative - a single "hurdle" token plus a
  `hurdle.positive = c("lognormal", ...)` argument (the type = "pbart"/"lbart" shape) -
  is recorded but not recommended, on the same internal-consistency grounds hazard
  settled (survival.md section 2). xbart / rbart_vi omit the token (their match.arg
  vectors are the refusal, the nbinom precedent).
- **The response.** A non-negative outcome with exact zeros (semicontinuous: a spike at
  0 plus a continuous positive part). The wrapper splits it R-side at ingest: the
  occupancy response is z = 1{y > 0} over all n; the positive response is log(y[S]) over
  S. Validate y >= 0 (refuse negatives by name, the nbinom validation precedent,
  negative-binomial.md). No latent split is drawn - the zeros are observed.
- **The fit object + generics.** A dedicated class `bartHurdle` (the bartMultinomial /
  bartNegbin idiom, R/bart.R:1269,1746) holding the two component fits (or their packaged
  draws) plus a marker recording the variant. Generics combine the two components' draws
  by sample index (any pairing is a valid joint draw, section 0):
  - predict / fitted / extract, type = "ev"/"response": combined E[y | x] =
    pi(x) * E[y | y > 0, x]. For the lognormal part, E[y | y > 0, x] =
    exp(f(x) + o + sigma^2 / 2) per draw - which REQUIRES the positive part's sigma
    draws, so those are a first-class output (the AFT survival-curve-needs-sigma shape,
    survival.md). Flag the retransformation honestly: the parametric lognormal
    mean assumes normal log-residuals; Duan's smearing estimator is the classic
    distribution-free alternative, a documented option/door.
  - type = "prob": the occupancy probability pi(x) (through the correct link).
  - a positive-part channel: E[y | y > 0, x], and both components' trees / varcounts
    reachable per part.
  - print reports "family: hurdle (probit occupancy + lognormal)".
  Because each component fit carries the ordinary $family element ("probit"/"gaussian"),
  every existing $family-dispatched generic (probabilityFromLatents, pointwiseLogLikelihood,
  R/generics.R:42-58) stays correct on the components unchanged; the hurdle-level combine
  is new code keyed on the bartHurdle class / marker, the hazard $periods split
  (survival.md section 4). predict requires keepTrees (the predict.bart guard).

## 7. Base-speed neutrality

Trivial and total under R composition: the engine has NO new code, so every non-hurdle
fit is byte-identical by definition - stronger than the nullable-member discipline
heteroscedastic needed (it had a real engine member to gate behind if(varianceForest_)).
All existing equivalence scenarios stay bitwise identical with no re-record (section 5
(ii)); bench-sampler is untouched (no engine path moves). This neutrality is itself an
argument for fork 1's R route: it removes the whole class of "did the second-forest
member perturb the single-forest hot path" risk that heteroscedastic's equivalence
22/22 gate had to guard (heteroscedastic.md section 8).

## 8. Scope: minimal v1

- Semicontinuous two-part: probit occupancy on 1{y > 0} over all n + gaussian-on-log-y
  positive part over the y > 0 subset, combined at predict.
- family = "hurdle" on dbarts()/bart2(); the bartHurdle class + fitted/predict/extract/
  print; y >= 0 validation.
- The reduction gate (section 5), the recovery smoke, the equivalence scenario, tinytest.
- No engine code, no bridge, no state-format change, no dbarts.h change.

## 9. Doors

- **Count hurdle (zero-truncated Poisson / NB positive part).** Needs a zero-truncated
  count single-forest family first (its own arc, section 3); then the composition wrapper
  accepts it family-agnostically as "hurdle.nbinom". negative-binomial.md filed it.
- **Gamma positive part.** A new continuous positive family; token door.
- **Logistic occupancy link.** "hurdle.logistic", one token away (section 6).
- **Grouped hurdle (random effects on either part).** Each part is a shipped family that
  GroupedResponse already decorates (rbart_vi); a grouped hurdle composes two grouped
  fits R-side, exactly as v1 composes two ungrouped fits. Surface-only follow-up.
- **Smearing retransformation.** Duan's distribution-free E[y | y > 0] estimator as a
  predict option beside the parametric lognormal mean (section 6).
- **Zero-INFLATION and sample-selection (the coupled cousins - where the engine finally
  pays).** A zero-inflated model (latent structural-vs-sampling zeros) and a
  Heckman-type correlated two-part / sample-selection model (correlated occupancy and
  amount errors) do NOT factor - the parts couple, R composition breaks, and a coupled
  sampler is required. THAT is the model that would justify the Chain two-response
  generalization sketched in section 4, because then there is a coupling for it to own.
  Recorded as the natural occasion for the engine work, distinct from this hurdle.
- **dbarts.h exposure.** None in v1 (universal precedent).

## 10. Budget

R-only, the hazard shape (survival.md's discrete-time landing was R-only). Approx:
the wrapper (ingestion split, two-fit orchestration, draw-combining), the family token +
validation + Surv-guard-style routing, the bartHurdle class + generics, the reduction
gate script + equivalence scenario + tinytest. Estimate ~400-700 lines R + gate scripts,
across ~2-3 commits (ingestion + wrapper; generics + packaging; gates + docs). NO engine,
bridge, or state work.

Contrast the rejected engine route (section 4), budgeted like heteroscedastic's C2 -
the two-response Chain surgery with the subset view and doubled state - at ~1500-2500+
lines of engine C++, the widest review in the queue, plus a state-format bump and its own
exact gate. The ~4x cost asymmetry, for a model with no coupling to exploit, is the
concrete form of fork 1's recommendation.

## 11. Plan-vs-note reconciliation

forest-combiner.md and multi-forest-models.md frame hurdle as the model
whose "per-forest response families break Chain's single response_" - i.e. as an engine
model requiring the two-response Chain break, and (heteroscedastic.md) as "the other
two-leaf-type model" that would sharpen heteroscedastic's combiner-vs-Chain-integration
sub-fork. FINDING: the break is REAL but AVOIDABLE - hurdle's conditional independence
(section 0) means it should be composed R-side, not built in the engine, so the invariant
break stays unbuilt. Consequences to record at landing (the orchestrator owns the edits):

- forest-combiner.md's hurdle bullet ("Response families differing per forest breaks
  Chain's single response_ - Chain-level, not combiner-API") is CORRECT about what an
  engine hurdle would need, but this note supersedes its implicit assumption that hurdle
  WILL be an engine model: v1 hurdle adds no Chain code. Cite this note.
- The "second consumer that shapes the two-leaf-type split" heteroscedastic.md deferred
  to hurdle DOES NOT ARRIVE via hurdle. Heteroscedastic's dedicated Chain-level integration
  therefore stays a single-consumer design; the engine two-leaf-type generalization waits
  for a model with genuine coupling AND two leaf types (the sample-selection cousin,
  section 9), not hurdle.
- core-generalization.md lists hurdle under "multi-forest Sampler"; accurate as a
  statistical description (two forests) but the IMPLEMENTATION is R composition of two
  single-forest samplers, not a multi-forest Chain. Note the distinction at landing.

## 12. Risks and confidence

- **Confidence in the factorization: HIGH.** The hurdle likelihood factoring into
  conditionally-independent occupancy and positive parts is textbook (Cragg 1971, Mullahy
  1986); the split is observed, not latent, so there is no coupling term to have missed.
  The single sharpest way to be wrong here is to conflate hurdle with zero-inflation
  (which DOES couple) - section 0 pins the distinction, and it is the first thing a critic
  should check.
- **Confidence in the R-composition sufficiency: HIGH.** The discrete-time hazard arc
  proved a factored two-part survival construction ships cleanly R-side with a bitwise
  reduction gate (survival.md); hurdle is the same pattern with two families instead of
  one on expanded data.
- **Top open question a critic should hammer #1: is the shared-cut-grid / one-object
  convenience worth engine work?** The note argues no (section 2). A critic who values a
  single unified sampler object, or who wants the positive forest to split on the FULL
  data's cut grid (not the subset's), should press whether the R route's two-store split
  costs prediction consistency. Rebuttal: explicit cutpoints in R force a shared grid
  without engine work; but VD should confirm the default (subset-derived cuts) is
  acceptable.
- **Top open question #2: the retransformation.** E[y | y > 0] = exp(f + sigma^2/2)
  assumes normal log-residuals; on heavy-tailed cost data the parametric lognormal mean
  is biased and Duan's smearing is the applied-standard fix (section 6). v1 ships the
  parametric mean and documents the assumption; a critic may argue smearing should be v1,
  not a door.
- **Top open question #3: does "hurdle" without a count variant undersell the name?** The
  count hurdle (section 3b) is the zero-inflated-count competitor many users mean by
  "hurdle"; shipping only the semicontinuous lognormal two-part may read as incomplete.
  Rebuttal: the wrapper is family-agnostic on the positive part, so the count variant is
  a pure follow-up once a truncated-count family lands - but VD may want the count arc
  sequenced alongside, or the family named "twopart" to set expectations.

## 13. Resolutions (2026-07-20; survey- and critique-grounded, authoritative)

An independent blind critique returned SOUND WITH FIXES and a software/literature
survey (brms, GLMMadaptive, pscl, statsmodels, twopartm/twopm, Stata churdle, mhurdle;
Duan 1983, Manning 1998, Manning-Mullahy 2001) grounded the two open forks. These
resolutions supersede the tentative in-body defaults of sections 3 and 6.

- CORE (confirmed): compose in R (section 2), semicontinuous lognormal positive part
  for v1 (section 3). Both the critique and the survey confirm; no BART/tree package
  ships either variant, so there is no in-family token to match.

- NAMING (supersedes section 6's bare "hurdle"): the canonical token is
  `hurdle.lognormal`, matching brms (`hurdle_lognormal`) and GLMMadaptive
  (`hurdle.lognormal`) verbatim - the closest Bayesian/mixed analogs, whose audience
  overlaps dbarts'. The `.lognormal` qualifier disambiguates from the count "hurdle" of
  pscl/statsmodels (a bare token would set the wrong truncated-count expectation).
  `twopart` is an accepted alias (health-econ audience) but prints as `hurdle.lognormal`.
  Future siblings: `hurdle.gamma`, `hurdle.nbinom` (count, once a zero-truncated count
  family ships).

- REPORTING (supersedes section 6's naive parametric-mean default; dissolves critique
  finding 1): predict/fitted default to the NATURAL (response) scale, computed
  model-consistently and heteroscedasticity-aware via posterior-predictive Monte Carlo -
  the convention every shipped two-part tool follows (twopartm/twopm/churdle natural
  scale; brms posterior_epred returns exp(mu + sigma^2/2)*(1 - hu) with sigma allowed
  per-observation). E[y|y>0,x] = exp(f(x) + sigma^2/2) consuming whatever variance the
  positive fit carries (single sigma^2 homoscedastic; the per-observation s(x)^2 surface
  when heteroscedastic BART is enabled on the positive part); E[y|x] = Pr(y>0|x) *
  E[y|y>0,x]. The MCMC posterior-predictive integrates the retransformation and per-obs
  sigma captures the log-scale heteroscedasticity Manning (1998) warns of - the case a
  single-sigma retransformation cannot, which dbarts's just-landed variance forest can.
  The retransformation plumbing MUST consume per-observation sigma from v1 even though the
  v1 DEFAULT fit stays homoscedastic, so enabling `variance = ~x` on the positive part
  "just works." Opt-ins: `type = "link"`/`"log"` (the log-scale linear predictor); Duan
  (1983) smearing as a distribution-free robustness estimator (a door - v1's lognormal
  model does not admit the non-normal-log-residual misspecification smearing guards).

- CRITIQUE-DRIVEN HARDENING (folded into the plan): (a) the gate is NOT the tautological
  "each component reduces to its standalone fit at the same seed" (true by construction
  for an R wrapper) - it adds a first-class analytic oracle for the COMBINE +
  retransformation (a hand-checked known-pi/mu/sigma example) plus explicit
  predict-on-newdata and save/load tests, the only genuinely new code; (b) INDEPENDENT
  deterministic per-component RNG seeds (a shared seed correlates the two chains and
  biases the combined credible interval, not the mean); (c) the positive fit receives the
  full-n x as x.test so in-sample fitted()/extract() has E[y|y>0,x] at the zero rows it
  never trained on; (d) a proper bimodal `type="ppd"` (draw Bernoulli(pi), then the
  lognormal), which the plain gaussian ppd path cannot produce; (e) the R route forecloses
  a SHARED variable-selection prior across the parts (a coupling the default independent
  spec lacks) - documented, with zero-inflation and Heckman selection staying the
  engine-justifying coupled doors (section 9).
