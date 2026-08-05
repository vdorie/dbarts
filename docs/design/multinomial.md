# Multinomial softmax: design

Status: LANDED, 2026-07-15 (commit bb8855e). The first multi-forest MODEL, and
the second consumer of the forest combiner (docs/design/forest-combiner.md).
K symmetric constant-leaf forests couple through a single-trial softmax
likelihood, augmented by an interleaved one-vs-rest Polya-Gamma cycle and pinned
by a likelihood-invariant level-centering move. The combiner
(MultinomialForestCombiner<L>, src/bartcore/combiner.hpp), the response family
(MultinomialResponse, src/bartcore/model.hpp), the K-forest chain constructor
and buildMultinomialForest (src/bartcore/chain.hpp), and the exact-posterior
gate (benchmarks/R/multinomial-exact.R) are all posterior-defining; creation is
internal, mirroring BCF's .Call surface. Scope: single-trial classification
(labels 0..K-1) and a fixed, sigma-free likelihood.

## The model

Category k gets its own forest f_ik = f_k(x_i); the K forests couple only
through the likelihood

    P(y_i = k) = softmax(f_i)_k = exp(f_ik) / sum_j exp(f_ij).

The forests are symmetric - one spec builds all K (MultinomialForestSpec,
src/bartcore/combiner.hpp) - rather than K-1 fit against a held-out reference
category. Symmetry buys exchangeability of the levels (no category is
privileged) and the ecosystem's shape: BART::mbart2, the closest precedent, is
fully symmetric with K-length variable counts and tree draws and a labeled
n x K probability output, and BART carries no coefficient table in which a
reference category's effect could hide. Every category then gets its own forest
and its own variable-count channel, so reporting is symmetric too. The price is
a per-observation additive non-identification of the raw f_ik (below), handled
by the centering move and sidestepped in all reporting by comparing the
identified probabilities. K-1-with-reference removes the invariance but makes
the reference category's reporting implicit and its variable counts asymmetric -
the wrong trade for a classification model whose deliverable is the K-vector of
probabilities.

## The interleaved Polya-Gamma augmentation

Category k's ONE-VS-REST conditional, given the other K-1 forests, is a binomial
logistic with linear predictor

    eta_ik = f_ik - C_ik,   C_ik = log sum_{j != k} exp(f_ij),

where C_ik is the log-sum-exp margin of the other categories. So
omega_ik ~ PG(1, eta_ik), and forest k sees working response
(y_ik - 1/2)/omega_ik + C_ik under precision omega_ik - the shipped logistic PG
machinery verbatim (ext_rng_simulatePolyaGamma; the construction in
LogisticResponse, src/bartcore/model.hpp), one binomial per category, with the
margin C_ik playing the role logistic's offset plays.

The augmentation is INTERLEAVED, not joint. The product of the K one-vs-rest
kernels integrates to a product of binomials, NOT the multinomial: each
omega_ik is a temporary latent valid only for category k's conditional, against
the margins it was drawn under. So omega_ik is drawn against the CURRENT margins
immediately before forest k's tree update, cycling the categories. A single
post-loop draw of all K omegas is not a valid Gibbs blocking: the one-vs-rest
kernels do not compose into a joint augmentation, and the batched form gives
forest k its working-response MEAN from fresh margins but its PRECISION from
stale ones - a Jacobi-style simultaneous update targeting the wrong invariant
distribution. Each omega_k is valid only against the margins it was drawn under.

The draw therefore lives in the per-forest pre-update hook drawForestGlue(f)
(src/bartcore/combiner.hpp), fired inside the sweep just before forest f's tree
update with the partially updated forests (0..f-1 new this sweep, f..K-1 old);
formForestResponse(f), called immediately after, reads the same margin and
omega through the combiner's own scratch. The base hook is a no-op consuming no
rng, so every additive combiner (BCF included) is bitwise unchanged - which is
how the hook landed without moving any existing draw.

Attribution: the interleaved one-vs-rest conditional cycle is Held and Holmes
(2006) and Polson, Scott and Windle (2013), section 4. Murray (2021),
"Log-Linear Bayesian Additive Regression Trees", is related work, not the
source of this construction - it augments with a per-observation normalizer
latent, a different device.

## Identification and the level-centering move

The softmax is invariant to a common per-observation shift of all f_ik (add the
same constant to every category's fit at observation i and the probabilities do
not move), so the raw f_ik are non-identified along a flat additive direction
the prior pins only weakly; left alone it mixes as a slow random walk. The
landed move pins ONLY the global, dataset-wide flat direction: a single scalar
shift c added to every f_ik at once, and ABSORBED UNIFORMLY - c/m_k onto every
occupied leaf of every one of forest k's m_k trees, plus c onto every totalFits
entry, which the next sweep's margins read. The residual roll's total =
sum-of-tree-fits invariant is therefore preserved to rounding
(m_k * fl(c/m_k) is not exactly c).

The conditional is exact in LEAF space, which is where the prior lives: each
leaf value is a priori N(0, s_k^2) with per-leaf sd s_k = nodeScale/(k sqrt(m_k))
(NOT N(0, tau^2) on the total fit - that would treat the n*K fits as independent
prior draws and both over-count the precision and, worse, hand the move a mean
that subtracts the whole level every sweep). With L_k and S_k the count and value
sum of forest k's occupied leaves over all its trees,

  prec = sum_k L_k / (m_k^2 s_k^2),   num = sum_k S_k / (m_k s_k^2),
  c    = -num/prec + N(0, 1)/sqrt(prec).

Empty leaves are skipped (they carry no fit). Uniform absorption is also the
better mixing device than loading c onto a single tree, whose own constant-leaf
conditional is not shift-equivariant and pulls back against the move; at the
intercept-only configuration the draw reduces to an exact independence sampler
from the level's marginal N(0, tau^2/K). The move lives in afterCombine
(MultinomialForestCombiner<L>, src/bartcore/combiner.hpp), the post-loop
combiner move, the BCF-ridge-interweave analog for the softmax's flat direction.

The shift is GLOBAL, not per-observation, by necessity. A per-observation shift
- the naive reading of the invariance - is not representable by shared-leaf
trees: a per-observation vector added to the fits cannot be carried exactly by a
piecewise-constant forest, so the next backfit PROJECTS the mismatch, leaking
spurious variance into the identified log-odds and biasing the reported
probabilities (a Jensen-type bias, measured at 2-4 percent against the exact
gate). A dataset-wide shift is the one flat direction a forest carries exactly -
add c to every observation uniformly - so it moves only the non-identified grand
level and leaves every identified log-odds f_ij - f_ik untouched. The move is a
mixing device, not a correctness requirement: the identified probabilities are
unbiased at any variance of the non-identified level, and the exact gate passes
with the move disabled. It only collapses the level direction's random walk.

## Leaf-scale calibration

The identified quantity is a DIFFERENCE of forests: a pairwise log-odds
f_ik - f_ij has sd sqrt(2)*s for per-forest total-fit sd s. Anchoring that at
the shipped logistic calibration (log-odds prior sd pi*sqrt(3)/2 at k = 2) fixes
s = pi*sqrt(3)/sqrt(2). Concretely (src/bartcore/combiner.hpp,
src/bartcore/chain.hpp): the per-forest node scale is
nodeScale = pi*sqrt(3)/sqrt(2) ~ 3.8476 (MultinomialSpec::nodeScale), the
per-leaf scale is nodeScale/sqrt(numTrees) (buildMultinomialForest), and the
per-forest total-fit prior sd is tau = nodeScale/k, with k = 2 the default, so
tau = pi*sqrt(3)/sqrt(2)/2 ~ 1.9238 - the sd the exact gate reads. The centering
conditional reads the per-LEAF sd s = nodeScale/(k sqrt(numTrees)) = tau/sqrt(m)
instead, since that is the scale its prior term is written on.

The margin C_ik makes the effective prior K-DEPENDENT: the softmax is a coupled
nonlinear map of the K forests, so no single per-forest scale matches the binary
calibration at every K. The K = 2 pairwise-log-odds anchor is exact only at
K = 2; at K >= 3 the same per-forest scale spreads its probability mass over
more categories and the implied prior on any one category's probability tightens
toward 1/K. The induced prior on a category probability p_1 - simulating f_k iid
N(0, tau^2), tau = pi*sqrt(3)/sqrt(2)/2, and softmaxing (4e6 draws, fixed seed):

    K   p1: q05    q25    q50    q75    q95    P(p1 > 0.9)
    2      0.011  0.138  0.501  0.863  0.989      0.210
    3      0.005  0.050  0.216  0.588  0.937      0.076
    5      0.002  0.019  0.084  0.294  0.785      0.020

The mean of p_1 is 1/K by symmetry at every K; the median falls below it (the
distribution is right-skewed) and the tail mass above 0.9 thins as K grows. The
K-dependence is inherent to the coupled prior, not a defect of the anchor.

## The exact-posterior gate

benchmarks/R/multinomial-exact.R, in the single-tree enumeration style of
benchmarks/R/bcf-exact.R and categorical-exact.R. Three arms, each matching the
sampler's posterior mean of the IDENTIFIED softmax probabilities to a closed-form
quadrature, to Monte Carlo error; the raw f_ik are compared never, and no arm
fixes a forest to zero (that would zero Var(f_K), break exchangeability, and
compare against the wrong posterior).

Arm 1, intercept-only K = 3. A constant predictor admits no valid cut points, so
every tree stays a root and each forest is a single leaf. The quadrature
integrates the sampler's actual symmetric N(0, tau^2)^K prior with the level
marginalized analytically: the induced prior on the differences
d_k = f_k - f_K (k = 1..K-1) is a CORRELATED Gaussian with covariance
tau^2 (I + 11') - each Var(d_k) = 2 tau^2, each Cov(d_k, d_j) = tau^2 - integrated
by 2-D quadrature. This gates the margin formation and the coupling in isolation
from tree growth.

Arm 2, K = 2 == logistic, and INTERCEPT-ONLY. The two-category multinomial
matches the shipped logistic sampler DISTRIBUTIONALLY (the two consume the rng
differently - K forests plus interleaved PG plus centering versus one forest - so
the match is distributional, not draw-for-draw). This is the only exact equality
with logistic, and only at intercept: the K = 2 multinomial log-odds is f1 - f0,
a DIFFERENCE of two m-tree ensembles, while logistic's is a SINGLE m-tree
ensemble; the two share the prior covariance at the sqrt(2) calibration but not
the function-space prior, so with a covariate they differ. Covariate K = 2 tree
growth is gated by arm 3, not by a logistic equivalence.

Arm 3, covariate-dependent K = 3. One binary categorical predictor (two cells),
one tree per forest; the joint tree space (each forest a root or a two-cell
split, 2^K combinations) is enumerated and integrated by nested per-cell
quadrature over the difference-space Gaussian - the only exact gate on tree
growth under softmax for K >= 3. Failure of any arm means the softmax coupling,
the interleaved PG draw, the centering move, or the margin formation is wrong;
the tolerances bound MC plus quadrature error and are never widened to pass.

## The surface

Creation runs through bartcoreMultinomialSampler (R/bartcore.R) ->
C_dbarts_bartcore_createMultinomial -> createMultinomialHolder
(src/R_interface_bartcore.cpp), which builds a MultinomialSpec and calls
createMultinomialSampler (a single ConstantGaussianLeaf instantiation, as BCF).
The public entry is bart2(family = "multinomial") (R/bart.R): a factor
response, K from levels(y), the level names threaded onto every K-shaped
output, probability-scale generics with an argmax class convenience, and a
fit class of its own (bartMultinomial) so no single-forest generic misreads
the K-widened arrays; its fit path reproduces the internal pattern bit for
bit (test-multinomial-surface.R pins it). Unsupported surface is refused by
name, never silently reshaped: weights, offset, and the latent type. Test
data at creation and out-of-sample prediction under keepTrees are both
supported (below).

- The response is an n x K nonnegative integer count matrix with row
  sums n_i >= 1 (trials): category k's one-vs-rest conditional is
  binomial(n_i, sigmoid(eta_ik)), augmented by omega_ik ~ PG(n_i, eta_ik)
  drawn as the sum of n_i PG(1, .) draws - exact, because the shape is
  observed integer data, never sampled, so the real-shape gap that
  constrains negative binomial's dispersion has no analog here (landed
  2bd34db, docs/plans/multinomial-counts.md). Single-trial labels
  (0..K-1 codes) enter as a one-hot counts matrix with unit trials, the
  byte-identical n_i = 1 reduction that anchors every recorded
  equivalence baseline; K defaults to one past the largest code on the
  label entry and to the column count on the count entry. Empty rows
  (n_i = 0) are refused at ingestion - PG(0, .) is a point mass at zero
  and would poison the working response. K = 2 counts reduce to
  binomial(n_i, p) distributionally (not bitwise - two forests, a
  different draw stream). The counts (or labels) are the response, so
  the host sampler's own response is ignored. Real-shape (non-integer)
  counts stay out of scope, and case weights are refused, both for the
  same real-shape PG gap.
- The run's train channel is the n x K x n.samples softmax probabilities, the
  identified deliverable (combinedFits writes the K location-major channels,
  log-sum-exp-safe). Per-category function values and split counts read through
  the per-forest queries bartcoreForestFits / bartcoreForestVariableCounts
  (0-based forest = category). The fit also carries a per-sample per-category
  variable-count channel: bart2(family = "multinomial")$varcount is
  (n.chains x) n.samples x p x K, levels on the K margin and predictor names
  on the p margin, mirroring every other K-shaped fit field.
- The formula interface (bart2(y ~ x1 + x2, data = df, family =
  "multinomial")) is accepted beside the matrix interface (landed
  2026-07-17, docs/plans/multinomial-formula.md); family = "multinomial"
  is never auto-detected from a factor response - a multi-level factor
  left-hand side under any other family setting is untouched by this and
  still whatever it did before (an error, from dbartsData's own
  response-to-numeric coercion). The response is pulled via
  model.frame/model.response with NO type coercion: dbartsData's own
  formula ingestion cannot be reused for this, since it coerces the
  response to numeric and would discard a factor's levels or a
  cbind(...) matrix's column names. A factor left-hand side routes to
  the label path (K and levels(y) as above); a cbind(c1, ..., cK) ~ x
  left-hand side (the glm binomial idiom) routes to the count path
  above, K and levels from the cbind column names. The right-hand side
  is handed to the host sampler build as the term-labeled predictor data
  frame, uncoded: dbartsData's own data-frame-as-x.train branch codes
  it, choosing the categorical or indicators builder via 'factors'
  exactly as any other family's formula fit does, which is what threads
  term.labels/factor.levels onto the host sampler's data@x and lets
  predict.bartMultinomial's validateXTest code a data.frame newdata - no
  separate terms/xlevels retention was needed. A formula fit reproduces
  the equivalent matrix fit bit for bit at the same seed, the
  reproduction gate extended to the formula surface
  (test-multinomial-surface.R).
- Test fits are DEFINED (testFitsAreDefined() true): test data supplied at
  creation (x.test, wired by createMultinomialHolder before the run) flows to
  all K forests, and combinedTestFits (MultinomialForestCombiner<L>,
  src/bartcore/combiner.hpp) blends their totalTestFits into the same K
  softmax combinedFits already applies to totalFits, reported in the run's
  test channel (yhat.test, same shape as train). The level-centering shift c
  is common to every forest and is NOT applied to totalTestFits (afterCombine
  leaves it alone); softmax invariance to a common per-observation shift (the
  identification section, above) makes the blend correct without touching it.
  A test offset is refused at creation (the softmax carries no offset), and a
  mixed/sparse test store cannot arise (multinomial requires dense
  predictors).
- Out-of-sample prediction (predict.bartMultinomial, R/generics.R) replays the
  K forests' saved trees at newdata into the same K-location slab and
  softmaxes it through the same map, coding newdata through the same
  training-level validation as creation; it requires the fit built with
  keepTrees (refused otherwise, mirroring predict.bart's own guard) - the
  multinomial chain constructor gained the saved-tree initialization it had
  silently lacked, so keepTrees now actually allocates storage on this path.
  Predicting at the fit-time test matrix reproduces the run's test channel
  bitwise, since both sum the same saved per-forest leaf fits through the
  same softmax. The bitwise fixture (benchmarks/R/multinomial-equivalence.R)
  records five channels: train, per-category forestFits, per-category
  varcount, the per-sample per-category run varcount, and test.
- The log-likelihood channel stays flagged undefined
  (logLikelihoodIsDefined() false): storeSample scores one forest's fits
  through response_->computeLogLikelihood, which cannot see the K-blend -
  BCF's exact NaN-flag choice, untouched by this arc.
- Whole-data mutation is refused, as for every multi-forest sampler:
  setData/setResponse/setWeights refuse on any sampler with >= 2 forests
  (refuseMultiForestMutation), and the transactional (non-force) predictor
  paths - including the per-observation sessions, which have no force variant
  - refuse likewise (refuseMultiForestTransactionalUpdate), because
  applyNewData and revalidateTrees rebuild only forest 0. A forced
  whole-matrix setPredictor, which refreshes every forest, stays available;
  the TEST predictor family (setTestPredictor and friends) is a separate
  gate, refuseBCFTestSurface, keyed on testFitsAreDefined() rather than forest
  count, so it now passes multinomial through (above) and refuses only BCF.
- State carries NO combiner wire blocks. The K forests serialize through the
  existing per-forest tree list (any length), and the multinomial combiner
  overrides neither serializeGlue nor restoreGlue - it serializes nothing. omega
  is a per-sweep latent redrawn against whatever margins the restored forests
  present (the interleaved draw seeds omega_k against restored forest k on the
  first restored sweep), so restore is STRUCTURAL, not bitwise - matching the
  interleaved draw's own semantics and the structural-restore contract every
  sampler follows. This is why a non-BCF combiner need not always carry glue
  state, correcting the blanket implication in
  docs/design/forest-combiner.md.

## No probit path

BART's multinomial precedent (BART::mbart) exposes both a probit (pbart) and a
logit (lbart) latent path as a user choice. This model is PG-softmax only - a
deliberate performance carve-out, a one-way door taken knowingly. The softmax is
one coupled likelihood served by a single augmentation mechanism (the
interleaved one-vs-rest PG), reusing the shipped PG(1, .) sampler and the
weighted-conjugate kernels verbatim. A multinomial probit would instead need
truncated multivariate-normal latent utilities and a covariance sampler -
per-category latent truncation machinery reusing none of the PG code. A probit
path is not precluded: it could live behind the same combiner seam, a probit
combiner drawing its own latents in drawForestGlue, if it is ever wanted. It is
simply not built.
