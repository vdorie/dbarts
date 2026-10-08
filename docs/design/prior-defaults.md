# Prior defaults

Status: reference, current. Documents shipped defaults as-is; not a design
proposal and carries no landing date - update in place whenever a default
changes.

Every default on the current surface, its source, and what it interacts
with. Values are unchanged from BayesTree/BART/dbarts history; this is a
record, not a re-derivation.

## k (leaf-prior scale)

Continuous responses default to 2 (`dbartsPriors$normal()`, `bart2`);
`bart` matches it verbatim. CGM's argument: the leaf prior is
`N(0, sigma_mu^2)` with `sigma_mu = 0.5 / (k sqrt(m))` for `m` trees.
Summed over the forest, `Var(f(x)) = m sigma_mu^2 = 0.25 / k^2`
independent of `m`, so `k` prior standard deviations of `f(x)` always
span the coded response range `[-0.5, 0.5]`; `k = 2` puts that range at
roughly a 95% prior interval. Binary responses (probit/logistic) default
instead to the `chi(1.5, 2)` hyperprior - an empirical choice, proper,
not independently derived. See "Response scaling" below for the
interaction with the response transform and `leaf.scale`.

## power, base (tree prior, `cgm()`)

2 and 0.95. CGM's empirical recommendation for the split-probability
decay `base * (1 + depth)^-power`, tuned in their experiments to favor
shallow trees without forbidding deeper ones. Not derived; adopted as-is
across BayesTree/BART/bartMachine/dbarts.

## df, quant (residual variance prior, `chisq()`, `family = gaussian(sigma = )`)

3 and 0.9. CGM's calibration: the inverse-chi-squared prior's scale is
picked so a rough sigma estimate sits at the `quant` quantile with
`df` degrees of freedom - an aggressive default (substantial prior
mass below the naive estimate), not a derived one. The engine's
`ChiSquaredScalePrior` (`src/bartcore/model.hpp`) inherits the mechanics
verbatim from the classic engine; the calibration itself is CGM's.

## leaf.scale (latent reference range)

3.0 for probit - a bare anchor (no formula ties it to anything else; it
is simply the assumed spread of the latent index), inherited unchanged
from the classic engine. Logistic's `pi * sqrt(3)` is mechanically
derived from it: multiply by the ratio of the logistic and normal latent
standard deviations (`(pi / sqrt(3)) / 1`), so the logistic leaf prior
spans the same number of latent standard deviations as probit's rather
than picking an independent constant (`R/dbarts.R`, `leaf.scale`
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
values before fitting. The in-package alternative is the `chi(1.5, 2)`
hyperprior on `k` (default for binary responses, available for
continuous ones too via `leaf.prior = normal(chi(1.5, 2))`) - letting
the leaf scale adapt some of the outlier's effect away rather than
letting a fixed `k` absorb it.

Because the data's scale is fixed at creation under a k-named prior, a Gibbs sampler
that swaps `y` or the offset between draws (the `dbartsSampler` use
case) would otherwise let the effective prior drift as the range
changes. `setResponse` and `setOffset` both take `updateScale` (default
`FALSE`): locked reuses the scale fixed at creation, so ordinary
sampling and between-draw substitution never drift; `updateScale =
TRUE` re-anchors and is documented as burn-in only, since re-anchoring
mid-run makes fits across iterations no longer comparable.

## Naming the leaf prior: k or sd

Locking the scale is not the same as choosing it. A composed model -
one whose driving R program hands the sampler latents, residuals or
another block's offsets - still inherits whatever spread the
CONSTRUCTION vector implied, which is an accident of how the outer loop
was initialized rather than a modelling statement. So the leaf prior is
named one of two ways, never both: `k`, relative to the scale the data
fixes (a number, or a law on k, `chi(df, scale)`), or `sd`, the prior sd
of the forest total's leaf parameter on the scale the family's forest
fits (a number, or a law on the sd itself, `invchi(df, scale)`). For the
constant leaf `sd` is the prior sd of f(x); for linear leaves it is each
coefficient's, per standardized covariate; for GP leaves the amplitude.

Only the ratio of k.scale to k enters a draw law. Under k ~ s chi_df the
spread is (k.scale / s) / chi_df, a scaled inverse chi on the sd, so a
named k.scale and the k hyperprior's scale are not separately identified:
`chi(df, s)` is `invchi(df, k.scale / s)`, and `chi(df, Inf)`, the
improper sd^-(df + 1), is `invchi(df, 0)`. That family is the only law on
the sd offered because it is the one conjugate to the normal leaves.

The translation rides a reference k of 2: `sd = x` reaches the engine
as k.scale 2x with k fixed at 2, and `invchi(df, c)` as k.scale 2c with
k ~ chi(df, 2); `invchi(df, 0)` states no scale and runs as chi(df, Inf) against the data's. A drawn
k starts at 2, so the chain starts at the named spread, and the binary
default and the old k spellings at their defaults keep bitwise engine
inputs. A named sd's k.scale is the dbartsModel slot `prior.scale` (NA under a k-named prior), which
overrides the family-keyed `leaf.scale` above and is converted
engine-side against the transform in force, at creation and on every
model install.

A named sd is absolute. The sampler restates the named sd's `prior.scale` after
every channel that re-anchors the response transform
([`reissueNamedLeafSd`](../../R/dbarts.R)), using the latest
`$setLeafPrior` write, which the R5 model records; a k moves with the
data. A `$setLeafPrior` write into a drawn prior keeps each chain's spread
in force at the call, from a fixed or a drawn prior in either spelling
(dec-B356, dec-B369, dec-B392, dec-B393): it reads k.scale before and after
the write and multiplies each chain's drawn k by their ratio
([`Chain::scaleDrawnK`](../../src/bartcore/chain.hpp)), so the new prior acts
from the next draw of k; a fixed k or sd stated sets the spread as written.
A re-anchor's reissue moves no k. The reader reports in the terms the prior was named in: `prior.sd`,
the sd law in force while k is drawn, and k relative to the data's
scale, reported as `k.scale` ([`reportLeafPrior`](../../R/dbarts.R)); a fit named by an sd
hyperprior carries draws of the sd in place of k. The two-forest and
multinomial models have their own calibration maps and refuse a named
sd rather than drop it, as does a hurdle fit, whose two parts are on
different scales.
