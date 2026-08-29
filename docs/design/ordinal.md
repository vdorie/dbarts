# Ordinal (cumulative-probit) outcomes: design

Status: LANDED 2026-07-18 (e0c4982). Plan: docs/plans/archive/ordinal-outcomes.md (this is its
step 1). Ordered categorical responses fit by a cumulative probit: the
truncated-normal latent machinery ProbitResponse already carries
(src/bartcore/model.hpp:3057), generalized from a single threshold at 0 to K-1
ordered cutpoints, plus a cutpoint block sampled by a marginal Metropolis
update. Surfaced as family = "ordinal". The unordered case is out of scope and
stays with the multinomial forest (docs/design/multinomial.md); the dispatch key
is is.ordered().

## 1. The model and the probit reduction

K ordered categories, response y_i in {1..K}. A latent utility

    z_i = f(x_i) + o_i + eps_i,   eps_i ~ N(0, 1),

for the sum-of-trees fit f, offset o_i, and cutpoints

    -inf = gamma_0 < gamma_1 < ... < gamma_{K-1} < gamma_K = +inf,

with y_i = k iff gamma_{k-1} < z_i <= gamma_k. Writing eta_i = f(x_i) + o_i the
category likelihood marginalizes the latent in closed form,

    P(y_i = k | eta_i, gamma) = Phi(gamma_k - eta_i) - Phi(gamma_{k-1} - eta_i),

the cumulative probit. This is exactly ProbitResponse with two changes: the
truncation interval per observation is (gamma_{y_i-1}, gamma_{y_i}] instead of
the fixed {(-inf, 0], (0, inf)} split (model.hpp:3085-3097), and the cutpoint
vector gains sampled values - gamma_1 stays pinned at 0 (section 2), so of the
K-1 finite cutpoints exactly K-2 are free. Everything else in the probit
seam is unchanged: sigma is fixed at 1 (drawSigma returns sigma,
model.hpp:3106), the tree stage sees working response z_i - o_i under unit
weights (workingWeights() null, rebuildWorking model.hpp:3074,3166), and the
log-likelihood channel is the Phi-difference above.

At K = 2 there is one cutpoint, gamma_1, fixed at 0 (below): the free-cutpoint
set is empty, the cutpoint block is skipped entirely, and the latent draw is the
one-sided truncation at 0 - i.e. **ordinal K = 2 is probit, bit for bit**,
PROVIDED the boundary categories route through probit's one-sided rejection
primitives and its NaN fallback (the two invariants of section 4); with them the
K = 2 path consumes the identical rng stream and the identity is a load-bearing
gate anchor (section 7). Without them it degrades to distributional equivalence.

## 2. Decision - identification

RESOLVED (VD 2026-07-18): scheme A. VD's condition - that pinning gamma_1 = 0
not complicate the latent distribution - holds: the pin is a pure location
anchor (Albert-Chib 1993, the MCMCoprobit normalization); the latent
conditionals keep their truncated-normal form and probit at K = 2 is this
scheme exactly.

The likelihood depends on gamma and eta only through gamma_k - eta_i, so it is
invariant to (i) a common location shift of all gamma and eta and (ii) a common
positive rescaling of all gamma, eta, and the latent sd sigma. A cumulative
probit needs one anchor for each. The BART mean f is a flexible location that
floats, softly shrunk toward the offset by the leaf prior (total-fit prior sd
nodeScale/k; probit sets nodeScale = 3.0, R/model.R:405), so f carries the
location that an explicit intercept would in a linear model.

**Scheme A (recommend): sigma = 1 fixed, gamma_1 = 0 fixed; free interior
cutpoints gamma_2 < ... < gamma_{K-1} (K-2 of them).** The unit latent sd anchors
scale exactly as probit does; gamma_1 = 0 anchors location, leaving f to carry
the rest. Identified quantities: eta_i on the unit-variance latent scale and the
K-2 free cutpoints relative to gamma_1 = 0; the reported deliverable, the K
category probabilities, is identified at any leaf-prior scale. The prior on the
latent mean is probit's: f ~ N(0, (nodeScale/k)^2) at the total-fit level
(nodeScale = 3.0, k = 2 default -> sd 1.5), reused unchanged. This is the exact
generalization of the shipped probit (K = 2 identical), keeps sigma fixed so
sigmaIsFixed_ stays true (chain.hpp:666-668) with no sigma draw and no new sigma
path, and confines all new machinery to the free cutpoints and their sampler.
Cost: as K grows the free cutpoints must mix (section 3, the hard part); and
two DISTINCT posterior correlations need naming. First, the Albert-Chib gamma-z
coupling - cutpoints pinned between adjacent latent order statistics - which the
marginal cutpoint update of section 3 removes outright by integrating z out.
Second, the f-vs-cutpoint ridge: f's grand level and a uniform shift of the free
cutpoints move the likelihood only weakly. gamma_1 = 0 removes the exact flat
direction and the leaf prior softly anchors f's level, but the ridge lives in
the MEAN, so it survives any cutpoint update that conditions on f - marginalizing
z does not touch it. v1 accepts it as mitigated, not eliminated; the true remedy,
if it ever bites, is a joint (f-level, gamma) shift move v1 does not provide
(section 9).

**Scheme B: sigma free, fix gamma_1 = 0 and gamma_2 = 1; free sigma plus K-3
cutpoints.** Estimates the latent scale, which some econometric treatments want.
Cost: needs K >= 3 to fix two cutpoints, and at K = 3 there are zero free
cutpoints - all shape rides sigma - so there is no smooth K = 2 reduction and
ordinal no longer unifies with probit; sigma free re-introduces the
sigmaIsFixed_ = false path, a sigma draw, and a scale that confounds with the
leaf-prior scale (both multiply eta), forcing a nodeScale recalibration; and the
unit gap gamma_2 - gamma_1 = 1 is an arbitrary scale the data must absorb. The
mixing study (section 3) shows the marginal cutpoint update fixes mixing without
a free sigma, so B buys machinery for no measured benefit.

**Scheme C: anchor via the response scaling / node.scale calibration.** Not a
competing anchor - it is the leaf-prior-scale sub-decision inside A. With sigma =
1 and gamma_1 = 0 fixed, the only remaining freedom is how diffuse f is, set by
nodeScale. v1 reuses probit's nodeScale = 3.0. Honest caveat (the multinomial
K-dependence, docs/design/multinomial.md): the induced prior on any single
category probability depends jointly on the leaf-prior sd of f and the cutpoint
spacing, and tightens toward 1/K as K grows; the anchor is exact only in the
probit sense at K = 2. A K-aware nodeScale is a follow-up; v1 documents the
induced prior rather than reshaping it.

Precedent survey. Scheme A (unit scale, one cutpoint pinned) is Albert-Chib
(1993, JASA 88:669), MCMCpack::MCMCoprobit (which normalizes gamma_1 = 0 with
standard-normal errors), and - in a variant that drops the intercept and frees
all K-1 cutpoints instead of pinning one - brms cumulative() and the Stan
ordered vector. dbarts pins gamma_1 = 0 rather than dropping an intercept
because BART's mean IS the (flexible) intercept-carrying location; there is no
scalar intercept to drop, and pinning gamma_1 = 0 is what makes K = 2 collapse
to probit. Scheme B (free sigma, two cutpoints fixed) is Nandram-Chen (1996) and
bqror's ORII model (used only at exactly 3 outcomes, where sigma is of interest).
No BART implementation ships a cumulative ordinal model at all: BART::mbart /
mbart2 and MPBART (Kindo-Wang-Pena 2016) are unordered; the one ordered-probit
BART, OPBART (Lee-Hwang 2024, Stat 13:e643), uses scheme A with a plain
Albert-Chib Gibbs cutpoint draw - which section 3 improves on.

## 3. Decision - cutpoint sampler and draw order

RESOLVED (VD 2026-07-18): as recommended - Cowles-style marginal MH, one
cutpoint at a time, plain random-walk proposal with out-of-bounds rejection,
fixed count-derived proposal scale.

**Draw order (decided): cutpoints, then latents.** The response family's
refreshLatents fires once per sweep after every tree updates, with the current
totalFits (chain.hpp:1535). Ordinal splits it into (1) a marginal Metropolis
update of the free cutpoints against eta = totalFits + offset, with the latents
integrated out (the Phi-difference likelihood of section 1), then (2) a fresh
draw of z_i from the doubly-truncated N(eta_i, 1) on (gamma_{y_i-1}, gamma_{y_i}]
under the new cutpoints, then rebuildWorking. The order matters: the cutpoint
move must NOT condition on the current z (that is precisely the Albert-Chib
pathology below); drawing z afterward hands the tree stage a working response
consistent with the accepted cutpoints. At K = 2 step (1) is empty and step (2)
is probit's one-sided draw, so the convention degenerates to probit exactly.

The candidates, and a mixing study (no trees - see the honest caveat) comparing
cutpoint effective sample size (ESS) across samplers, at
benchmarks/... (script currently untracked, not yet
promoted to benchmarks/R/ordinal-cutpoint-mixing.R; promote it if
kept). Ordinal probit, sigma = 1, gamma_1 = 0, an intercept+slope linear mean
drawn by a conjugate Gibbs step so the location floats; 6000 iterations, 1500
burn-in, seed 20260718. ESS is the mean over free cutpoints (initial-positive-
sequence estimator):

    n     K   AC(min/mean)   Cowles   log-gap
    200   3      30 / 30        472       319
    2000  3      10 / 10        543       402
    200   5      -- (stalled)  7 / 269   8 / 291
    2000  5       8 / 14      267 / 312  111 / 224

- **Albert-Chib (1993) Gibbs.** gamma_k ~ Uniform(max latent in category k, min
  latent in category k+1), latents pinned by adjacent order statistics. ESS
  DEGRADES as n grows (K = 3: 30 -> 10) because the order statistics crowd
  together, and at K = 5, n = 200 it STALLS: a sparse boundary cell leaves the
  window one-sided (Uniform(lo, +inf)), which a flat-prior AC cannot draw from;
  the prototype skips such updates, freezing that cutpoint (ESS undefined). This
  is the documented pathology (correlation of cutpoints with the latent data,
  worsening in n). Rejected as the primary sampler; it is the reference the gate
  reproduces distributionally (its target is correct, only its mixing is bad).

- **Marginal Metropolis, Cowles-style (recommend).** A PLAIN (untruncated)
  normal random-walk proposal on each cutpoint; a proposal outside
  (gamma_{k-1}, gamma_{k+1}) is rejected outright, and an in-bounds proposal is
  accepted on the marginal category likelihood ratio times the prior ratio (only
  the two adjacent cells' Phi-differences change), the latents integrated out;
  then z is redrawn. The plain-RW kernel is SYMMETRIC, so that ratio is the whole
  acceptance - exactly what the prototype validates. Pinned deliberately: a
  literal truncated-normal proposal under the same acceptance would be an INVALID
  sampler (asymmetric; it needs the truncation-normalizer ratio in the
  acceptance), so the spec is plain-RW-with-rejection, not Cowles' original
  truncated proposal. In the study it delivers 15-54x AC's ESS and is stable or
  improving in n (K = 3: 543 at n = 2000 vs AC's 10; K = 5 survives n = 200 where
  AC stalls). Recommended as the primary update, one cutpoint at a time on the
  ordered scale (simplest, needs no proposal covariance, adapts per-cutpoint).

- **Log-gap reparameterization (documented alternative / proposal mechanism).**
  Map the ordered cutpoints to unconstrained reals delta_j = log(gamma_{j+1} -
  gamma_j), a joint random-walk MH with the Jacobian sum(delta) in the acceptance,
  z redrawn after. This is the Stan/brms ordered-vector internal transform and
  bqror's map; it enforces ordering automatically and blocks the cutpoints. In the
  study it tracks Cowles closely (K = 3: 319/402) and is the better structure when
  cutpoints are strongly correlated at large K, at the cost of a proposal
  covariance to tune. v1 ships one-at-a-time Cowles; the log-gap block is the
  recorded upgrade path if K is routinely large.

HONEST CAVEAT. The study has NO trees: a 2-parameter linear mean, not a BART
ensemble. A flexible tree mean competes with the cutpoints for the location far
more aggressively, so these ESS numbers are an OPTIMISTIC bound on the in-sampler
mixing, not a BART prediction. The transferable finding is the RANKING (AC
collapses in n and at sparse cells; marginal-MH and log-gap are an order of
magnitude better and stable), which matches the qualitative literature (Cowles
1996; Nandram-Chen 1996). The absolute ESS is not a claim about dbarts. Further
prototype limits: the MH acceptances use a FLAT cutpoint prior (no prior ratio),
so the study validates the kernel ranking, not the shipped log-gap prior; the
post-cutpoint z redraw in the script is overwritten unused at the next loop head
(each iteration redraws z before the mean update), so the study never exercises
the shipped cutpoints -> latents -> trees within-sweep coupling directly - the
coupling flows only across iterations; and each (n, K) cell is a single
replicate, so the Cowles-vs-log-gap differences are noise - only the
order-of-magnitude AC gap is signal.

**Sub-decision - fixed vs adaptive proposal scale.** The optimal random-walk
scale for gamma_k shrinks like 1/sqrt(n_k + n_{k+1}) as the adjacent-cell counts
grow (the marginal-likelihood curvature scales with the count near the cutpoint),
so a single hard-coded scale cannot be robust across n. Two ways to get an
n-robust scale: (a) a FIXED scale computed once at construction from the observed
cell counts, e.g. c / sqrt(n_k + n_{k+1}) with c ~ 2-3 - a fixed Markov kernel,
fully reproducible, nothing extra to serialize, gate-friendly; or (b) diminishing
adaptation during warmup toward a target acceptance, then frozen - better tuned
but non-Markov during warmup, and the tuned scale must serialize and the warmup
must be deterministic to stay bitwise reproducible. Given dbarts's reproducibility
priority (within-host bitwise across any sample count), **recommend (a), the
count-derived fixed scale**; (b) is the fallback if fixed scales mis-tune at
extreme cell imbalance.

**Cutpoint prior.** A weakly-informative proper prior on the K-2 free gaps keeps
the posterior proper even when a boundary or interior cell is empty (the case that
breaks flat-prior AC): concretely a normal prior on the log-gaps delta_j with
mean and sd set so the default spacing is O(1) on the unit latent scale.
Shipped constants (C1): log-gap prior N(0, 1.5^2) (default spacing exp(0) = 1
on the unit latent scale), proposal-scale constant c = 2.5, and the interior
categories' NaN fallback is the interval midpoint (the boundary categories keep
probit's sign * DBL_EPSILON). The exact gate (section 7) and any reference
computation must use these same constants. The
classic improper flat prior on the ordered region (Albert-Chib, Cowles) is the
recorded alternative - fine when every cell is occupied, and it is the limit the
proper prior approaches as its sd grows. A proper prior is also required for the
exact-posterior gate's quadrature to integrate (section 7).

## 4. Seams

**OrdinalResponse (new, src/bartcore/model.hpp).** Not a decorator - a sibling of
ProbitResponse holding the same per-observation latents_ and working_ buffers
plus a K-1 cutpoint vector and the fixed per-cutpoint proposal scales. It reuses
probit's rebuildWorking (working = latents - offset), workingWeights() = nullptr,
drawSigma() = sigma, initialSigma/fitScale/fitShift/sigmaScale = 1/1/0/1. The
overrides are refreshLatents (the two-step cutpoints-then-latents draw of section
3), computeLogLikelihood (the Phi-difference, generalizing probit's two-tail
form at model.hpp:3148-3163), and the new cutpoint state hooks (section 6). Ordinal
carries K on construction (levels count).

**A doubly-truncated-normal primitive is required.** random.h shipped only
one-sided truncations (ext_rng_simulateLowerTruncatedNormalScale1 /
Upper..., random.h:86-107); probit uses the lower form at 0 (model.hpp:3092).
That primitive has since been added, as specified here:
ext_rng_simulateTruncatedNormalScale1 is declared at random.h:108-115.
Interior categories need z ~ N(eta, 1) on a finite interval (gamma_{k-1},
gamma_k], so add ext_rng_simulateTruncatedNormalScale1(rng, mean, lower, upper):
inverse-CDF (mean + qnorm(Phi(lower-mean) + u (Phi(upper-mean) - Phi(lower-mean))))
in the bulk, Robert (1995) rejection when the interval sits in a tail where the
CDF gap underflows. Boundary categories (y = 1 and y = K) MUST keep the existing
one-sided rejection primitives, never the new inverse-CDF one: the upper draw
ext_rng_simulateUpperTruncatedNormalScale1(mean, bound) is
mean - LowerTruncatedStandardNormal(mean - bound) (random.c:486-493), the same
underlying primitive and argument folding as probit's sign-flipped draw
(model.hpp:3092), and the NaN fallback must match probit's sign * DBL_EPSILON
(model.hpp:3093). These two invariants are what make the K = 2 identity of
section 1 BITWISE rather than merely distributional - a different primitive
consumes a different rng stream - and the K = 2 component test (section 7) locks
both. The doubly-truncated primitive itself is a self-contained external/
addition with its own component test (section 7).

**Family selection (chain.hpp:597-637, facade.hpp).** Add ResponseFamily::ordinal
to the enum (model.hpp:2580) and a case to the chain's family switch constructing
OrdinalResponse(y, offset, K, ...). No leaf-model change - ordinal is a
ConstantGaussianLeaf single-forest model like probit, so it composes with the
existing SamplerFacade instantiations (facade.hpp:413+) unchanged; K threads in
through the same options struct the survival status and group indices use
(chain.hpp:774-791). The bridge is NOT a string addition: resolveFamily
(src/R_interface_bartcore.cpp:1581-1619) branches on the boolean
control.responseIsBinary and refuses every family but gaussian/aft for a
non-binary response, so ordinal needs a third response-shape channel - a K-level
categorical flag plus K itself - plumbed through ParsedControl beside
responseIsBinary, with resolveFamily accepting "ordinal" only on that shape and
refusing it by name everywhere else.

## 5. Decision - surface and prediction semantics

RESOLVED (VD 2026-07-18): auto-dispatch with announcement, family = "ordinal"
as the explicit primitive. VD's condition - detection by a concrete class, not
level-name inference - holds: orderedness is the class attribute
c("ordered", "factor") set by factor(ordered = TRUE)/ordered(), tested by
is.ordered(); level names are never consulted. Prediction shape (n x K labeled
probabilities) and the state design (section 6) were presented alongside and
drew no objection.

**Family surface (recommend family = "ordinal" as the primitive, with
auto-dispatch on ordered factors, announced).** family = "ordinal" is added to the
bart2 and dbarts family vectors (R/bart.R:699, R/dbarts.R:382) and routed through
resolveClassificationFamily (R/data.R:588). The response must be an ordered
factor (is.ordered); the level ORDER defines the category order. An unordered
factor or character under family = "ordinal" is accepted with an informational
message that order is taken from the level order (factor default: alphabetical -
usually not what is meant, hence the warning); an integer/numeric response uses
sort(unique(y)) as the ordered levels.

BEHAVIOR AT PROPOSAL TIME (the baseline every option changes; the option that
won has landed, so this paragraph is the pre-landing baseline, not the live
tree). An ordered K >= 3 factor then took two different paths.
bart2(family = "auto") routed it through detectAutoMultinomial, whose type
match included "ordered factor" at n.levels >= 3, and SILENTLY fit an unordered
multinomial, discarding the ordering. The single-forest entries (dbarts, xbart,
rbart_vi) errored in resolveClassificationFamily's K >= 3 refusal
(R/data.R:619-673). So EVERY option below - explicit-only included - had to
split ordered factors out of detectAutoMultinomial, and bart2's behavior on
ordered responses changed whichever option won; there was no
zero-behavior-change choice. The detection hook already existed:
classifyResponse tags "ordered factor" as its own response type
(R/data.R:498-510).

LIVE: the split is made on the disjoint is.ordered() key (stated at
R/bart.R:1344-1346). detectAutoMultinomial (R/bart.R:1391) matches unordered
factors and characters only; detectAutoOrdinal (R/bart.R:1405-1421) catches the
3+-level ordered factor and selects family = "ordinal". The single-forest
entries still refuse, and R/data.R:619-621 names ordinal as the reason.

The fork on auto-dispatch - a genuine decision point for VD, with the internal
and external precedents on OPPOSITE sides:

- **Auto-on-ordered with announcement (recommend).** family = "auto" routes
  is.ordered(y) -> ordinal and announces it ("ordered response -> family =
  'ordinal'; set family to override"); unordered/character K >= 3 stays
  multinomial; explicit family = "ordinal" is the primitive and forces the model
  on any response. The INTERNAL precedent is squarely here: dbarts's own binary
  handling auto-dispatches probit from the unambiguous 2-level signal via
  announceAutoFamily (R/utility.R:23-33), and bart2's auto already routes
  unordered K >= 3 factors to multinomial - under explicit-only, ordered factors
  would become the ONE factor response family = "auto" refuses, incoherent with
  the package's own auto concept. is.ordered(y) is an unambiguous type signal,
  exactly what auto exists to consume. Against: the ordinal ECOSYSTEM
  (ordinal::clm, brms cumulative(), rstanarm::stan_polr, MCMCoprobit) always
  requires an explicit family/function - none auto-dispatches on is.ordered -
  and ordered() is sometimes set for contrasts or printing, not modeling, so
  auto makes a type attribute load-bearing for MODEL choice; the announcement is
  the mitigation, not an elimination.

- **Explicit-only.** family = "auto" errors on an ordered K >= 3 factor with a
  hint to set family = "ordinal" (a louder, better error for the single-forest
  entries; for bart2 it converts today's silently-wrong multinomial fit into an
  error). Matches the ecosystem convention above and the spirit of multinomial's
  never-auto-detect rule (docs/design/multinomial.md) - though that rule is
  already breached by bart2's own auto peek, which weakens it as precedent.
  Against: it leaves the unambiguous is.ordered signal inert, makes the user
  restate in family what the type declares, and is inconsistent with the
  package's binary auto-dispatch.

Recommend auto-with-announcement: once bart2 must change behavior anyway, the
internal-consistency argument dominates, and the announcement keeps the dispatch
visible. Whichever wins, family = "ordinal" is the primitive and ordinal never
overlaps unordered multinomial: the disjoint key is is.ordered().

**Prediction semantics.** fitted()/predict(type = "ev"/"response") return the
n x K category-probability array (x n.chains x n.samples where extract does),
P(y = k | x) = Phi(gamma_k - eta) - Phi(gamma_{k-1} - eta), columns labeled by the
ordered levels - the multinomial K-column shape (R/generics.R:1584-1598), computed
through the cutpoints rather than a softmax. type = "bart"/"link" returns the
single-column latent eta = f (probit returns its latent likewise,
R/generics.R:400-404). A class-prediction convenience maps the argmax category
back to the ordered level (multinomial's max.col + levels idiom). The K-1 threshold
draws are posterior output too: expose an n.samples x (K-1) "thresholds" field, the
ordinal analog of gaussian's sigma, so users can reconstruct probabilities at
arbitrary eta and inspect the thresholds. predict requires keepTrees (the
predict.bart guard, R/generics.R:302-304). This K-column path deliberately
DIVERGES from the plan's suggestion that probabilityFromLatents generalizes
(docs/plans/archive/ordinal-outcomes.md step 1): probabilityFromLatents
(R/generics.R:12-19) is a scalar link-inverse of a single latent column and
stays binary-only; an ordinal probability is a Phi-DIFFERENCE against per-sample
cutpoints, a two-argument transform that does not fit that seam, so ordinal gets
its own transform beside it rather than a generalization of it.

**Levels round-trip.** Single-forest probit currently DISCARDS the factor levels
(codeResponse, R/data.R:552-566). Ordinal must persist the ordered levels on the
fit object to label the K probability columns and map argmax back - as
bartMultinomial persists $levels (R/bart.R:1710). Give ordinal its own S3
class bartOrdinal (mirroring bartMultinomial, R/bart.R:1729) carrying $levels and
$thresholds, so the bart generics do not misread the K-widened arrays as a
single-forest fit.

## 6. Decision - state and mutation

**Cutpoints in a new by-name state block - the resid.df / C2 pattern.** Add the
virtual trio carriesCutpoints() / cutpoints() / restoreCutpoints() to ResponseModel
(default false / nullptr / no-op), exactly as carriesResidualDf() /residualDf()/
restoreResidualDf() was added for Student-t (model.hpp:4169). serializeState
writes a "thresholds" slot only when carriesCutpoints()
(the residualDf pattern at chain.hpp:3084-3086); the R<->C++ state list gains a
named "thresholds" slot beside "resid.df" (slotNames, R_interface_bartcore.cpp:
6457-6466), written only when present (the conditional write at 6612-6619) and
read by name with absence tolerated (the getListElement decode at 7087-7095).
stateIsValid refuses an ordinal sampler whose state lacks a length-(K-1) cutpoint
block (the residualDf check, chain.hpp:3262); setState restores it
(chain.hpp:3753). **Old states (gaussian/probit/etc.) omit the slot and load
unchanged** - the whole point of the by-name additive block. The per-observation
latent z_i rides the EXISTING latents slot (probit's z already serializes there,
model.hpp:3136-3141); restoreLatents rebuilds working from z. Choosing the
count-derived FIXED proposal scale (section 3) means nothing else serializes; an
adaptive scale would add a second scalar-vector block.

**setResponse / setData cold-init.** setResponse (same n, new y - the embedded-
Gibbs swap): KEEP the current cutpoints (a slow-moving global parameter the outer
sampler wants persisted across a small y perturbation) and re-draw z from the
current fits under the new y's intervals - probit's minimal-disruption re-draw
(model.hpp:3110-3114), plus leaving gamma alone. setData (n changes, everything
stale): cold-init the cutpoints to the default prior spacing AND re-draw z from
scratch, matching probit/TResponse cold-init on a data swap (model.hpp:3121-3134,
4122-4132). State this rule in the class doc as the probit convention plus a
kept-cutpoint clause.

**setSigma / weights.** sigma is fixed at 1 (scheme A), so sigmaIsFixed_ is true
and the sigma prior/setSigmaPrior is a no-op, as for probit. Weights are
unsupported for the same reason probit refuses them (model.hpp:3053-3056): there
is no coherent weighted truncated-normal latent likelihood, so a weighted ordinal
probit is not a real model. Refuse weights at ingestion, by name, like probit.

## 7. Gates

**Exact-posterior gate (single tree, K = 3, small n).** In the single-tree
enumeration style of the logistic and multinomial gates
(benchmarks/R/multinomial-exact.R, categorical-exact.R). The key simplification:
the cumulative-probit likelihood MARGINALIZES the latents analytically (section 1
is a closed-form Phi-difference), so the exact reference is z-FREE and the gate
needs quadrature only over the leaf means and the one free cutpoint gamma_2
(gamma_1 = 0 fixed) - never over z. Agreement therefore validates BOTH the
truncated-normal augmentation AND the cutpoint sampler, exactly as robust-errors'
reference validates the scale mixture by never augmenting. Computation: enumerate
the tree structures a single predictor with a few cuts admits (root, or one split
into two leaves); for each, the marginal is

    integral over (mu_leaf..., gamma_2) of
      [ prod_i ( Phi(gamma_{y_i} - mu_{node(i)}) - Phi(gamma_{y_i - 1} - mu_{node(i)}) ) ]
      x prod_leaf N(mu_leaf; 0, (nodeScale/(k sqrt(numTrees)))^2)
      x p(gamma_2),

a 3-D quadrature at K = 3 with a two-leaf split (mu_left, mu_right, gamma_2),
2-D at the root. Three details keep the target exact: p(gamma_2) is the log-gap
prior pushed forward to the gamma scale INCLUDING the Jacobian - the prior is
normal on delta = log gamma_2 (scheme A's single free cutpoint at K = 3), so
p(gamma_2) carries the 1/gamma_2 factor and the gate integrates the sampler's
actual prior, not a flat stand-in; a leaf holding no observations contributes a
factor of 1 (its empty likelihood product integrates its mu prior to unity); and
each enumerated structure's posterior weight is its tree prior TIMES its computed
marginal likelihood, renormalized over the enumeration. Match the sampler's
posterior mean of the identified category probabilities and of gamma_2 to the
quadrature, to Monte Carlo error; tolerances bound MC plus quadrature error and
are never widened to pass.

**Component tests.**
- Truncation interval: on fixed gamma and eta the z_i draw's support is
  (gamma_{y_i-1}, gamma_{y_i}] and its mean/variance match the analytic
  doubly-truncated-normal moments (the new primitive tested in isolation,
  including a deep-tail interval where the inverse-CDF path must hand off to
  Robert rejection).
- Cutpoint conditional: on a tiny fixed (y, eta) the sampled gamma_2 histogram
  matches the exact 1-D posterior p(gamma_2 | y, eta) proportional to the
  Phi-difference product times the prior; and the marginal-MH acceptance ratio is
  checked against a direct evaluation.
- K = 2 identity: ordinal with K = 2 and gamma_1 = 0 reproduces the probit
  sampler BIT FOR BIT, locking the two section 4 invariants that make it true:
  the empty free-cutpoint set consumes no rng, boundary draws route through the
  one-sided rejection primitives (never the inverse-CDF primitive), and the NaN
  fallback is probit's sign * DBL_EPSILON. A strong equivalence anchor tying the
  new family to the shipped one.

**Recovery.** The parent plan's recovery gate (docs/plans/archive/ordinal-outcomes.md
step 3): simulated ordinal data over a nonlinear f at moderate n, checking
category-probability calibration and cutpoint recovery against the truth - the
family-level smoke beyond the exact gate, in the style of the existing recovery
tests.

**Equivalence fixture and neutrality.** A NEW scenario in
benchmarks/R/equivalence.R (an ordered-factor response, family = "ordinal", K in
{3, ...}) recording the ordinal channels: category probabilities, cutpoints, and
latents. Existing anchors are untouched - ordinal is a new family selected by a
new enum value, adds no draw to any existing family's stream, and reduces to
probit at K = 2 - so the frozen baselines and every RNG-locked snapshot stay
stable; the neutrality trail is verified by re-running equivalence.R compare and
expecting IDENTICAL draws for the existing families (no re-record), the
robust-errors precedent (docs/design/robust-errors.md section 6).

## 8. Out of scope, and the doors

- **Unordered multinomial.** Owned by docs/design/multinomial.md (the K-forest
  softmax). Ordinal is is.ordered() only; an unordered K >= 3 factor stays
  multinomial (or its refusal on the single-forest entries). The two never
  overlap - the dispatch key is is.ordered() - which is why ordinal is a
  single-forest cumulative probit and not a second multi-forest model.

- **Grouped / mixed-model ordinal (rbart_vi + ordinal).** Refused cleanly at the
  R layer, at the resolveClassificationFamily call in rbart_vi (R/rbart.R:338),
  before the group attribute is built; rbart_vi's family vector (R/rbart.R:48)
  omits "ordinal" for v1. The door is real: GroupedResponse is a base-response
  decorator whose conjugate group update needs a Gaussian working response with
  the group intercept on the latent scale (model.hpp:4668+, drawGroupEffects over
  base_->workingResponse()), and ordinal HAS exactly that - z is Gaussian given
  the cutpoints - so grouped ordinal is a coherent future once the cutpoint block
  and the group block are shown to interleave cleanly. v1 refuses; the composition
  is recorded as feasible.

- **xbart ordinal loss.** Out of scope, follow-up. The natural cross-validation
  metric for ordered categories is an order-aware loss (ranked probability score),
  not the misclassification/log-loss the current xbart carries; deferred to a
  separate pass.

- **dbarts.h exposure.** None in v1, the robust-errors precedent. The flat C API
  (inst/include/dbarts/dbarts.h) is unchanged; ordinal is reachable only through
  the R surface and the internal bartcore .Call path, so no LinkingTo consumer
  (stan4bart) sees an ABI change. The door: a future dbarts.h entry could expose
  the cumulative-probit family for embedded use, deferred until there is demand.

## 9. Costs and risks

- Cutpoint-mean confounding. Two correlations, only one of which v1 fixes
  (section 2). The Albert-Chib gamma-z coupling is removed by the marginal
  cutpoint update (z integrated out). The f-vs-uniform-cutpoint-shift ridge is
  NOT: it lives in the mean and survives any update that conditions on f;
  gamma_1 = 0 removes only the exact flat direction and the leaf prior softly
  anchors f's level, so the ridge is mitigated, not eliminated. At large K with
  a very flexible mean it could slow mixing; the true remedy is a joint
  (f-level, gamma) shift move - the centering-move analog - which v1 does not
  provide. The log-gap block (section 3) helps only the cutpoint-internal
  correlation, not this ridge.
- Sparse cells. An empty or near-empty category makes a cutpoint weakly
  identified; the flat-prior AC has no valid draw there at all (section 3 study,
  the one-sided window), which is why v1 uses a proper log-gap prior and the
  marginal-MH update. Document that a category with no observations cannot be
  placed from data alone.
- The K-dependent induced prior (section 2, scheme C caveat): the default
  nodeScale is calibrated at the probit K = 2 sense, so the prior on a single
  category probability tightens toward 1/K as K grows. Documented, not reshaped in
  v1.
- Reproducibility. The fixed count-derived proposal scale keeps the cutpoint
  kernel Markov and adds no serialized state, preserving within-host bitwise
  reproducibility across sample counts; an adaptive scale would trade that away.
- The new truncated-normal primitive is a numerical-robustness surface: the
  inverse-CDF underflows for tail intervals and must hand off to Robert rejection
  (section 4), covered by its own deep-tail component test.
