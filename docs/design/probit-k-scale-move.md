# The probit rescaling step: a parameter expansion for k under probit

Status: BUILT 2026-10-08 (dec-B371, dec-B394, dec-B397), not yet landed; plan
[probit-k-scale-move.md](../plans/probit-k-scale-move.md). The finding it
answers: [probit-k-calibration.md](probit-k-calibration.md).

## What it does

Once a sweep, ahead of the trees, a single-forest probit fit with a drawn k
multiplies every active latent and every occupied leaf by one factor alpha and
divides k by it, alpha drawn from its exact conditional
([`Chain::drawForestRescale`](../../src/bartcore/chain.hpp)). The posterior is
unchanged; k and the fit's overall scale move far faster. It has no switch
(dec-B397): it runs wherever it applies.

## Why k mixes slowly

One probit forest: k ~ s chi(nu) (chi(1.5, 2) by default); each occupied leaf
mu ~ N(0, (c / k)^2); latent z_i ~ N(o_i + f_i, 1), y_i = 1{z_i > 0}. Two
couplings hold the fit's scale still. Given the leaves, k is pinned to a
relative spread near 1 / sqrt(2 M) over M leaves (the funnel). Given z, the
fit's scale is pinned by z, and given the fit z sits within about one unit of
it; when the fit is large, near separation, their common scale moves by a
vanishing amount a sweep - the slowness of probit data augmentation that
parameter expansion was invented for. Both couplings lie along one direction:
multiply every leaf, every latent and 1 / k by one factor.

## The step

The group alpha > 0 acts as (z, mu, k) -> (alpha z, alpha mu, k / alpha) on the
active latents, the occupied leaves and k. The draw is the generalized Gibbs
step of Liu and Sabatti (2000): alpha has density proportional to the target at
the moved state times the Jacobian, against the group's Haar measure
d alpha / alpha. With n active rows, R = sum (z_i - f_i)^2, Q = sum o_i (z_i - f_i)
and C = k^2 / (2 s^2) (zero under an infinite scale):

    Jacobian                          alpha^(n + M - 1)
    leaf prior at (alpha mu, k/alpha) alpha^-M times its value at (mu, k)
    k prior at k / alpha              alpha^-(nu - 1) exp(-C / alpha^2)
    likelihood at alpha z             exp(-alpha^2 R / 2 + alpha Q)

so in v = log alpha

    log p(v) = (n - nu) v - (R / 2) e^(2v) + Q e^v - C e^(-2v).

Carrying k keeps eta = k mu / c invariant, so the leaf prior drops out and M
with it. The truncation 1{sign z_i = y_i} is invariant for alpha > 0. Empty
leaves sit at zero and stay there; an inactive row is not in the model, so its
latent is neither counted nor scaled; a leaf whose members are all inactive is
drawn from its prior and scales like any other. Without an offset alpha^2 is
generalized inverse Gaussian; with one it is not, so v is moved by one slice
step from v = 0 ([`sliceFromZero`](../../src/bartcore/chain.hpp); Neal 2003,
stepping out at width 0.25 with a step limit of 1000 split at random, then
shrinkage). One slice step is a reversible kernel for the conditional, which is
all a Gibbs component needs.

Exactness: the step is a Gibbs draw of the orbit coordinate given the
orbit-invariant coordinates (k z, k mu, the trees). Every decline - another
family, a fixed or infinite k, more than one forest or a combiner, a non-constant
or monotone leaf, a stale tree map, R = 0, C = 0 with n <= nu where the
conditional is improper - depends only on orbit-invariant quantities, so mixing
the step with the identity keeps the posterior. A decline consumes no generator
draw, so every fit the step does not apply to draws exactly as before.

Placed beside the level step, ahead of the tree loop, every recorded channel -
training and test fits, saved trees, k - is written after it. totalFits is
multiplied by alpha in place while the forest's accumulated factor since its
last re-derivation stays within [1/2, 2], and re-summed from the scaled leaves
in tree order when the factor would leave that range, so the cache gap a factor
multiplies stays bounded ([Forests and
combiners](../architecture.md#forests-and-combiners)). A re-sum every sweep
would cost 9 to 14 percent of a sweep (below).

## Alternatives measured on the prototype

- (a) Interweaving on k (Yu and Meng 2011): a non-centred draw of k given
  eta = k mu / c and z. Its conditional's precision in the fit's scale is
  sum f^2, so it moves freely only where the fit is small, where today's
  sampler already mixes; median tau of k 84 sweeps but a 99th percentile of
  31,000, and nothing on the separated toy (P(k < 0.135) 0.184 against an exact
  0.103). Dropped.
- (b0) The expansion with k held: median tau 129 against (b)'s 66 on the same
  25 datasets, worst 1711 against 682. Carrying k is worth a further 1.5 times
  at the median and halves the tail.
- (c) The collapsed move, k given eta and y with the latents integrated out:
  about three times faster per sweep than (b) but 46 percent more time a sweep
  at n = 5000 and 200 trees (23 percent by two-evaluation Metropolis), gaining
  mainly where k is not slow today. A door: it generalizes to any family with a
  cheap likelihood, which (b) does not.

## Other families (doors)

- logistic and nbinom (Polya-Gamma; k drawn by default): (b) does not apply, the
  augmentation's variates being precisions; (c) does. Whether k is slow there is
  TODO k-mixing-pg-families.
- ordinal: (b) extends by scaling the free cutpoints with the latents, their
  count in the Jacobian and their prior in the density.
- linear and gp leaves under probit with a drawn k: every coefficient's sd is
  scale / k, so (b) extends by scaling the parameter blocks and the fit slab.
- gaussian-family fits with a drawn k: (a) is the classical interweave beside a
  drawn sigma.

## Evidence

Measured on the built step, arm64 laptop unless said. Every figure without the
step was measured while the build carried a switch (dec-B394), since removed
(dec-B397); the scripts now run the step alone.

**Mixing** ([`runCensus`](../../benchmarks/R/probit-k-mixing.R)): the SBC
probit-k arm's generator (n = 150, 50 trees), 100 prior-drawn datasets, the
arm's start, 30,000 sweeps of burn-in and 200,000 recorded. tau is k's
integrated autocorrelation time in sweeps.

    step       tau(k): median   90th   99th    max   tau(avg f): median   90th
    off (prototype)       524   6108  23278  52751                  162   5354
    on                     74    429    830   1129                   26    369

By the dataset's true k0, median tau(k) with the step: 609 (k0 <= 0.3, two
datasets), 550 (0.3 to 0.6, four), 217 (0.6 to 1, 14), 124 (1 to 1.5, 17), 49
(above 1.5, 63); without it the prototype read 12797, 4744, 2711, 1749 and 293.

**Invariance** ([`runTruth`](../../benchmarks/R/probit-k-mixing.R)): a state
drawn exactly from the posterior run on for 3000 sweeps, R = 2000: the drift of
log k over lags 300 to 3000 is -0.022 (se 0.011), z -2.0, within the |z| 3
criterion. The prototype's Jacobian error of alpha^2 drifted it at z -16.5.

**Exact gate** ([probit-k-scale-exact.R](../../benchmarks/R/probit-k-scale-exact.R)):
full mode, |z| <= 2.1 on all 20 decile statistics and every leaf mean, the
separated arm's batch-means error 0.0043; with the step off that arm's error was
0.0238 and the gate failed (quick: z -15.2, error 0.019), as it must. Quick
mode's mixing bound, 0.014, sits between the step's 0.0064 (arm64) and 0.0099
(x86-64) and the 0.019 without it, so a step that silently declines fails it.
The gate cannot see k left unscaled after a correct draw of the latents and
leaves: that mutant passes quick mode, worst |z| 4.3. Only the mapping check in
[`runForestRescaleTests`](../../tests/cpp/test_ensemble.cpp) catches it.

**The SBC arm**: at thin 1000 every one of the 100 census datasets keeps 50 or
more effective draws of the 99 kept (median 95, least 54; an AR(1) reading of
the lag-1000 autocorrelation), against a quarter under 50 without the step.
The arm runs at thin 700 after a 10,000-sweep burn (dec-A190 as marked): the
fewest effective draws are then 49 by the lag-1000 autocorrelation read at lag
700 and 55 by k's autocorrelation time, the burn about nine times k's slowest
autocorrelation time (1,129 sweeps), at about 60% of thin 1000's cost. Its first run at those
settings, in the full matrix on 2026-10-09, passed every functional: k at
chi-square p 0.383, KS p 0.907 and ECDF difference 0.022 against the 0.080
band, its end bins 33 and 30 against 30 expected, in a 44-minute job.

**Cost.** Per-sweep time with the step on over off, the same build, data and
chain state (both samplers restored from one stored state before every timed
segment), ten interleaved repetitions in alternating order, one chain, one
thread; median ratio and range, on an idle x86-64 host (4 cores, load 0.5 to
0.8 over the as-built run, 0.1 to 1.7 over the others):

    n       trees   as built (bounded in place)   re-sum every sweep, four trees a pass   one tree a pass
    500       75    1.027 (1.014-1.047)           1.095 (1.081-1.100)                     1.126
    500      200    1.017 (1.000-1.077)           1.088 (1.059-1.109)                     1.121
    2000      75    1.008 (1.008-1.023)           1.104 (1.092-1.111)                     1.154
    2000     200    1.013 (1.004-1.031)           1.119 (1.020-1.145)                     1.169
    5000      75    1.016 (1.008-1.029)           1.119 (1.115-1.131)                     1.179
    5000     200    1.000 (0.968-1.033)           1.114 (1.060-1.133)                     1.189
    20000     75    1.004 (0.996-1.013)           1.108 (1.103-1.135)                     1.179
    20000    200    1.009 (1.004-1.025)           1.139 (1.115-1.144)                     1.220

The step's own arithmetic is a few percent at most at the smallest n and within
noise from n = 2000. The prototype's figures (+19 percent at n = 500, +1
percent at n = 5000) timed fresh samplers at different chain states.
