# Prior defaults

Every default on the current surface, its source, and what it interacts
with. Values are unchanged from BayesTree/BART/dbarts history; this is a
record, not a re-derivation.

## k (node-prior scale)

Continuous responses default to 2 (`dbartsPriors$normal()`, `bart2`);
`bart` matches it verbatim. CGM's argument: the node prior is
`N(0, sigma_mu^2)` with `sigma_mu = 0.5 / (k sqrt(m))` for `m` trees.
Summed over the forest, `Var(f(x)) = m sigma_mu^2 = 0.25 / k^2`
independent of `m`, so `k` prior standard deviations of `f(x)` always
span the coded response range `[-0.5, 0.5]`; `k = 2` puts that range at
roughly a 95% prior interval. Binary responses (probit/logistic) default
instead to the `chi(1.25, Inf)` hyperprior - an empirical choice
(improper, mildly penalizes small `k`), not independently derived. See
"Response scaling" below for the interaction with the response transform
and `node.scale`.

## power, base (tree prior, `cgm()`)

2 and 0.95. CGM's empirical recommendation for the split-probability
decay `base * (1 + depth)^-power`, tuned in their experiments to favor
shallow trees without forbidding deeper ones. Not derived; adopted as-is
across BayesTree/BART/bartMachine/dbarts.

## sigdf, sigquant (residual variance prior, `chisq()`)

3 and 0.9. CGM's calibration: the inverse-chi-squared prior's scale is
picked so a rough sigma estimate sits at the `sigquant` quantile with
`sigdf` degrees of freedom - an aggressive default (substantial prior
mass below the naive estimate), not a derived one. The engine's
`ChiSquaredScalePrior` (`src/bartcore/model.hpp`) inherits the mechanics
verbatim from the classic engine; the calibration itself is CGM's.

## node.scale (latent reference range)

3.0 for probit - a bare anchor (no formula ties it to anything else; it
is simply the assumed spread of the latent index), inherited unchanged
from the classic engine. Logistic's `pi * sqrt(3)` is mechanically
derived from it: multiply by the ratio of the logistic and normal latent
standard deviations (`(pi / sqrt(3)) / 1`), so the logistic leaf prior
spans the same number of latent standard deviations as probit's rather
than picking an independent constant (`R/dbarts.R`, `node.scale`
assignment).

## n.trees

75 (dbarts's own historical default, `dbartsControl`); BayesTree's and
`bart`'s default is 200. Both are unargued round numbers; `k`'s
derivation above holds for either since `sigma_mu` is defined in terms
of `m`.

## dart update.delay

Half of `n.burn`. The BART package's "startdart" convention: hold the
Dirichlet split-probability update until the forest has had time to
become informed by the data, because a cold, uniform-probability forest
under an immediately-sampled concentration is bistable.

## tau slice steps (grouped random effects)

`n.thin`. `rbart_vi`'s R loop reused the thinning interval as the slice
sampler's step count - a convention (matching tau's refresh cadence to
however many tree sweeps separate two kept draws elsewhere in the
chain), not a derived mixing requirement. The in-core engine reproduces
the count exactly (see `docs/design/grouped-random-effects.md`).

## Response scaling

Continuous responses are range-scaled: `y` (net of offset) is mapped to
`[-0.5, 0.5]` by its observed min/max, and every other constant above is
calibrated against that anchor. This is the convention of the entire
BART software lineage (BayesTree, BART, bartMachine), which is what lets
`k` (and `sigma_mu = 0.5 / (k sqrt(m))`) transfer across packages and
papers without re-derivation; changing it would invalidate every
cross-package comparison at once.

The known failure mode is outlier sensitivity: two extreme `y` values
stretch the range and compress the effective leaf prior for everything
else. bartMachine's JSS paper names the same issue and recommends the
fix dbarts offers no automation for: log-transform or winsorize extreme
values before fitting. The in-package alternative is the `chi(1.25,
Inf)` hyperprior on `k` (default for binary responses, available for
continuous ones too via `node.prior = normal(chi(1.25))`) - letting the
leaf scale adapt some of the outlier's effect away rather than letting a
fixed `k` absorb it.

Because the anchor is fixed from the data at creation, a Gibbs sampler
that swaps `y` or the offset between draws (the `dbartsSampler` use
case) would otherwise let the effective prior drift as the range
changes. `setResponse` and `setOffset` both take `updateScale` (default
`FALSE`): locked reuses the scale fixed at creation, so ordinary
sampling and between-draw substitution never drift; `updateScale =
TRUE` re-anchors and is documented as burn-in only, since re-anchoring
mid-run makes fits across iterations no longer comparable.
