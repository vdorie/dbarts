# benchmark-surfaces: test problems for measuring dbarts, and the battery to build

Status: COMPLETE (survey), 2026-09-06

Every measurement this package has taken of its own sampler has been taken
on Friedman's five-dimensional function or a near relative of it: the
move-set A/B and the response-swap recovery run both
([13. Move-set A/B (2026-09-06)](tree-mixing-proposals.md#13-move-set-ab-2026-09-06),
[14. Recovery after a response swap (2026-09-06)](tree-mixing-proposals.md#14-recovery-after-a-response-swap-2026-09-06)),
the grow-from-root default study, the composition probe's own generator
(`inst/common/friedmanData.R`), and most of `benchmarks/R`. Friedman's
function has two linear terms, one quadratic and one bounded smooth
interaction on five of ten uniform predictors; shallow trees fit it easily,
and a battery built on it cannot tell a kernel that mixes better from one
that does not. This document surveys the test problems the tree,
nonparametric-regression, causal-inference and MCMC literatures actually
use, verifies each against its primary source, and proposes a battery.

The framing constraint the battery is built to. dbarts is used as the
nonparametric residual component of semi-parametric models - stan4bart and
bartCause put a parametric block on the structured terms and hand BART what
is left - so the battery measures the AVERAGE case first: realistic
residual surfaces with moderate nonlinearity, a few interactions,
correlated mixed-type covariates, realistic n and p, moderate noise.
Pathologies are a secondary set, chosen to cover distinct failure modes
rather than to be maximally cruel. In a causal or embedded problem the
quantity judged is what the OUTER model reports - a treatment effect, a
random-effect scale, a coefficient - and BART's own fit is secondary.

The decision rule this battery exists to serve: **a kernel or prior change
is accepted only if it is neutral or better on the average-case core AND
better on at least one pathology. Never the reverse.** A change that wins
a pathology and loses the core is refused.

Scope: survey only. Nothing here is scheduled, no source is touched, and
no default moves. Section 6 is the deliverable; sections 1 to 5 are the
evidence, and section 8 lists what could not be fetched.

---

## 1. Tree literature

### 1.1 Chipman, George and McCulloch: what BART was measured on

The 2010 BART paper's synthetic case is Friedman's function,
`y = 10 sin(pi x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5 + eps`, with
`x1..xp` iid uniform, `eps ~ N(0,1)`, at **n = 100** and p = 10, then
p = 20, 100 and 1000 with n held at 100 [verified: arXiv 0806.3286 sec 5.2,
eqs 26-27]. That is the whole of the paper's controlled-truth evidence.
Its own verdict on mixing is a visual one: "We see that the BART MCMC
burns-in quickly and mixes well" (sec 5.2.1), and, in general,
"Compared to the single tree model MCMC approach of CGM98, our backfitting
MCMC algorithm mixes dramatically better... we have found that restarts of
the backfitting MCMC algorithm give remarkably similar results even in
difficult problems. Consequently, we run one long chain with BART rather
than multiple starts. Although mixing does not appear to be an issue..."
(sec 3.1). Every measurement in sections 1.6, 1.7 and 5 below contradicts
that last clause on larger data than the paper ran.

The **42-dataset bake-off** (sec 5.1) is a subset of the 52 sets of Kim et
al. (2007): "between 3 and 28 numeric predictors and 0 to 6 categorical
predictors", categorical predictors expanded to indicators, sample sizes
"from 96 to 6806 observations", 20 random 5/6-1/6 train/test splits each
(840 splits), scored by relative RMSE - "the RMSE divided by the minimum
RMSE obtained by any method" on that split, so 1.0 is the per-split winner.
Competitors were the Lasso, gradient boosting, random forests and a
one-hidden-layer neural net. Burn-in was "determined by inspection of a
single long run. Typically, 200 burn-in steps and 1000 iterations".

Two further facts from that paper bear directly on this battery. First,
variable inclusion is a small-m readout, not a default-m one: "this
strategy is less effective when m is large because the redundancy offered
by so many trees tends to mix many irrelevant predictors in with the
relevant ones" (sec 3.2) - which is the same self-averaging effect this
house measured at 75 trees. Second, the paper's one large binary example
(sec 5.3, a drug-discovery set with p = 266 molecular descriptors,
n = 29374 of which 542 active, probit, m = 50) is where the authors
actually reached for a convergence check: "we performed four independent
repetitions of 250,000 MCMC iterations and obtained essentially the same
results each time." A quarter of a million sweeps is the honest price of
that check.

CGM 1998, the Bayesian CART paper whose stochastic search the 2010 paper
compares against, could **not be fetched** (section 8).

### 1.2 Pratola 2016: the two problems built to be sticky

Pratola's two motivating examples are the only published pair designed
specifically so that tree mixing, not fit, is what fails.

- **The low-noise Friedman emulator** (sec 2.2): the deterministic
  Friedman function treated as a simulator, n = 5000, m = 200, 5000 burn-in
  plus 5000 kept. At `sigma^2 = 1` the birth/death-only sampler "was found
  to mix reasonably well, having an acceptance rate around 18% and the 90%
  credible interval having an empirical coverage of 81%." At
  `sigma^2 = 0.1` the acceptance rate is "just around 4%", the coverage
  "around 53%", and "the tree structure became stuck in a local mode with,
  for all practical considerations, zero chance of moving to a different
  area of tree-space" [verified: arXiv 1312.1895 sec 2.2]. Rotation at 20%
  of proposals took acceptance to about 25% and coverage to about 96%;
  adding perturb took acceptance to 65% at coverage 92% (sec 5.2, sec 6).
- **The Wu-Tjelmeland-West confounded step function** (sec 2.3): p = 3,
  n = 300, `y = 1 + eps` if `x1 <= .5, x2 <= .5`; `3 + eps` if
  `x1 <= .5, x2 > .5`; `5 + eps` if `x1 > .5`, with `eps ~ N(0, 0.25)`.
  The covariates are drawn in blocks so that **x1 and x3 are confounded**:
  `x1 ~ U(0.1,0.4)` for i <= 200 and `U(0.6,0.9)` after, while
  `x3 ~ U(0.6,0.9)` for i <= 200 and `U(0.1,0.4)` after. "We fit BART to
  this dataset using only m = 1 trees and found that the acceptance rate of
  tree moves (after the initial few steps of the sampler) was 0."
  [verified: arXiv 1312.1895 sec 2.3]

The second is the cheapest sticky problem in the whole survey: 300 rows,
three columns, an exactly known two-mode answer (x1 or x3, by symmetry of
the design), and a reported acceptance rate of zero. It is a
representation-multimodality probe with a checkable null, which is what
the XOR probe in `grow-from-root-default.md` was built by hand to be.

### 1.3 Linero and Yang 2018 (SoftBart): where smoothness beats trees

SoftBart's simulation section supplies three things this battery wants.
Its Friedman cell is n = 250, `sigma^2` in {1, 10}, "p from 5 to 1000 along
an evenly-spaced grid on the scale of log p" [verified: arXiv 1707.09461
sec 4.1]. Its variable-selection cell rescales the linear part,
`f = 10 sin(pi x1 x2) + 20 (x3 - .5)^2 + lambda (10 x4 + 5 x5)` for
`lambda` in [0.1, 1], and scores precision, recall and F1 on posterior
inclusion probability above 50% - a **signal-strength sweep**, which is a
better inclusion probe than a fixed design because it walks the detection
boundary rather than sitting on one side of it.

Its section 4.2 supplies two pathologies at the opposite extremes of
smoothness: a pure step, `f(x) = 2 - 4 I(x1 < 0.5)` at n = 250 and
`sigma = 0.1`, and "a highly localised Daubechies wavelet of smoothness
order 10", where "the fit of BART, by contrast, possesses many artifacts
outside the support of the wavelet, and possesses generally wider credible
bands" [verified: arXiv 1707.09461 sec 4.2]. The wavelet cell is the
inhomogeneous-smoothness case of section 2 in tree clothing.

Its benchmark table (sec 4.3) is ten datasets - ais, abalone, bbb, cpu,
diamonds, hatco, servo, tecator, triazines, wipp, mostly again from Kim et
al. (2007) - scored by 5-fold predictive error normalized to SBART-CV.
One row matters here: **tecator**, where BART-CV scores 1.87 and 1.63 for
DART against 0.98 for SBART, and the authors say "leveraging smoothness for
this dataset is essential to attaining good performance". That is the
cleanest published instance of a real dataset where a smooth method beats
trees by a wide margin.

Linero 2018 (DART itself) could **not be fetched** (section 8); its test
functions are known to this survey only through the SoftBart paper's
description of them.

### 1.4 He, Yalov and Hahn 2019, He and Hahn 2023: a DGP factory

XBART's suite is four fixed mean functions on `d = 30` iid standard normal
predictors, `sigma = kappa Var(f)` for `kappa` in {1, 10}, at
n in {10000, 50000, 250000} [verified: arXiv 1810.02215 sec 4.1, Table 1]:

    Linear         x^T gamma,  gamma_j = -2 + 4(j-1)/(d-1)
    Single index   10 sqrt(a) + sin(5a), a = sum_{j=1}^{10} (x_j - gamma_j)^2
    Trig + poly    5 sin(3 x1) + 2 x2^2 + 3 x3 x4
    Max            max(x1, x2, x3)

The 2023 journal version turns this into a full factorial, and it is the
single most reusable generator in this survey [verified: arXiv 2002.03375v4
sec 4.1]. The mean functions are chosen "to cover a range of important
special cases: linearity, additive models, models with interactions,
nonlinear smooth functions, and functions with discontinuities". The
predictor matrix is drawn either independent standard normal or
**correlated with factor structure**: `k = p/5` factors, loadings 0/1 with
exactly five ones per column and one per row so `BB^T` is block diagonal,
`X = (BF)^T + e` with `e ~ N(0, 0.01k)`, columns rescaled to unit standard
deviation. Errors are Gaussian or `t_3` rescaled to the same variance,
with `sigma^2 = kappa^2 Var(f)` and `kappa` in {1, 10}. The (n, p) grid is
p = 30 at n in {10000, 50000, 250000}; p = 100 at n = 1000; p = 500 at
n = 300; p = 1000 at n in {500, 1000}.

The decisive number for this package is in their section 5. At n = 10000
and `kappa = 1`, 95% pointwise coverage of the true mean function, averaged
over 100 replications, is **0.77 for BART on Linear, 0.78 on Max, 0.74 on
Trig+Poly**, against 0.91 to 0.96 for the same BART sampler warm-started
from 25 XBART forest draws, and BART is the slowest of the three
[verified: arXiv 2002.03375v4 Table 4]. Their reading: warm starting
"yields considerable improvement in the estimation, which may indicate
inadequate chain length of BART (that is, poor mixing)." That is a
30-percentage-point coverage deficit at MODERATE noise on a plain n = 10000
fit - not in a low-noise tail, not at one tree - and coverage of the true
mean is exactly the statistic section 6.3 of the mixing survey nominated.

### 1.5 Kapelner and Bleich 2016 (bartMachine)

bartMachine's synthetic case is Friedman again, at n = 500 with p = 5
signal columns and **p0 = 95 noise columns**, `sigma = 1`, used to show the
effect of an informed split prior [verified: JSS 70(4) sec 4.10]. Its
interaction-detection readout counts, per draw, pairs of variables
appearing on a root-to-leaf path together (sec 4.11) - a structural
functional, so it inherits every caveat this house has recorded about
structural readouts at large m. Its shipped convergence panel is worth
noting as prior art: sigma by iteration, **percent acceptance of MH
proposals per iteration**, average leaves per tree and average tree depth
(Figure 4). Its nine-dataset bake-off is boston, triazine, ozone,
baseball, wine.red, ankara, wine.white, pole and compactiv, 10-fold CV
averaged over 20 replicates; random forests beat both BART implementations
on triazine, ozone and pole (Table 3).

### 1.6 Ronen, Saarinen, Tan, Duncan and Yu 2022

This is the paper that ran dbarts itself. Two experiments, and the survey
record in `tree-mixing-proposals.md` carries only the second.

**Experiment 1 is the one this battery should take.** Four PMLB datasets -
Breast Tumor (n = 116640, p = 9), California Housing (20640, 8), Echo
Months (17496, 9), Satellite Image (6435, 36) - subsampled to n = 200,
n = 2000 and full; **full BART at ntree = 200**, 8 chains, nskip = 5000, 20
replicates; the statistic is the Gelman-Rubin diagnostic of the held-out
test RMSE across the 8 chains. The finding: "even when running the chain
with more burn-in samples than recommended, we observe that BART fails to
mix on all of the original datasets", and the trend is monotone in n
[verified: arXiv 2210.09352 sec 4.1]. One number is worth carrying whole:
"the between-chain variation for the Breast Tumor dataset is no more than
0.1% of the RMSE value, while the GR value is 1.74." A diagnostic that
fires at a between-chain spread of one part in a thousand is measuring
structure, not fit.

**Experiment 2** is Bayesian CART at one tree with the full move set: "the
root split changes in less than 0.2% of the samples on average across 160
chains", and "for full datasets, an overwhelming majority of the root
splits occur on the same feature, and furthermore, this feature is
different for different chains" [verified: arXiv 2210.09352 sec 4.2].

The authors' own non-transfer sentence - "We did not find strong evidence
that this bottleneck affects the BART algorithm to the same degree" -
attaches to the ROOT-SPLIT bottleneck of Experiment 2, not to Experiment
1's Gelman-Rubin failure, which they report for BART at m = 200 without
qualification. `tree-mixing-proposals.md` section 3.3 uses the sentence
correctly for the mechanism it is about; nothing in this package's record
yet carries Experiment 1.

### 1.7 Tan, Ronen, Saarinen and Yu 2024

The follow-up keeps the four PMLB sets and adds two synthetic DGPs
"borrowed from recent papers on heterogeneous treatment effect modeling"
[verified: arXiv 2406.19958 sec 9.1.3]:

- **Low-Dimensional Smooth** (Lei and Candes 2021): `x in R^10`,
  `x ~ N(0, Sigma)` with `Sigma_ii = 1` and `Sigma_ij = 0.01`;
  `y = g(x1) g(x2) + eps` with `g(x) = 2 / (1 + exp(-12(x - 0.5)))`; noise
  calibrated to a signal-to-noise ratio of 3.
- **Piecewise Linear** (Kunzel et al. 2019): `x in R^20` a Gaussian copula
  with the same weak equicorrelation; three regimes selected by `x20` at
  cuts -0.4 and 0.4; within each regime the mean is a linear combination of
  a DIFFERENT disjoint block of five covariates, coefficients drawn once
  from `U(-15, 15)`; noise variance 1.

Their protocol: m = 200, `pi_g = pi_p = 0.25, pi_c = 0.4, pi_s = 0.1`
(dbarts' own mixture), 8 chains, 1000 burn-in, 10000 kept, 25 replicates,
n in {200, 500, 1000, 10000}, R-hat on held-out RMSE plus 95% coverage.
Four of their findings matter here:

- R-hat increases with n on every dataset tested.
- Raising the temperature dampens the trend and improves coverage.
- "Increasing the number of trees consistently dampens the trend in R-hat."
- "Initializing the chains from a fitted XGBoost ensemble... or simply
  increasing the number of burn-in iterations can increase R-hat values and
  exaggerate their increasing trend. We believe this suggests that the
  chains become stuck around different local maxima for the target
  posterior."

The last is the sharpest warning in this survey for anyone reading a
between-chain statistic: on these problems a BETTER-fitting initializer
makes the diagnostic WORSE, because it drops each chain into its own mode.
It is also the mechanism behind this package's own KILLED grow-from-root
warm-start default, seen from the diagnostic side rather than the accuracy
side.

Their two theorems name the two structural pathologies precisely. Theorem
5.1: for an additive `f = f1(x1) + ... + fm'(xm')` with independent
components fitted at `m <= m'` trees, the hitting time for an optimal
representation is `Omega_P(n^{1/2})`, worsening to `Omega_P(n^{qmin/2-1})`
if change and swap are disallowed and `m < m'`. Theorem 5.2: if `f*`
contains a **pure interaction** - XOR over binary features is their
canonical case - and change moves are disallowed, the hitting time is
`Omega_P(n^{1/2})` with a suboptimality gap wide enough that "the MSE of
the BART sampler output does not even converge to zero with the training
sample size" [verified: arXiv 2406.19958 sec 5.1, 5.2]. Both theorems are
conditional in a way that matters to dbarts: the additive bound bites at
FEW trees, and the XOR bound assumes `pi_c = 0`, which the shipped default
is not.

And their Experiment 7 is the null this package already reproduced:
"restricting the move set does not substantially affect R-hat, coverage,
or RMSE" [verified: arXiv 2406.19958 appendix L.6].

### 1.8 Breiman's oblique problems, and the checkerboard

Breiman's three synthetic classification problems, at 300 training cases
each, are the canonical axis-alignment stressors [verified: Berkeley TR
460, sec 2.3]:

    twonorm    20-dim, 2 class, unit covariance, means +/-(a,...,a),
               a = 2/sqrt(20).            Optimal boundary: an oblique plane.
    threenorm  class 1 an equal mixture of N(+a...), N(-a...); class 2
               centered at (a,-a,a,-a,...).  Boundary: two joined oblique
               hyperplanes.
    ringnorm   class 1 N(0, 4I); class 2 unit covariance at (a,...,a),
               a = 1/sqrt(20).            Boundary: a sphere.

"These problems are difficult for CART. For instance, in twonorm the
optimal separating surface is an oblique plane. This is hard to
approximate by the multidimensional rectangles used in CART... Threenorm is
the most difficult". The caveat that must travel with them is in the same
paragraph: "Yet in all examples CART has low bias. The problem is its
variance." An ensemble is a variance reducer, so an oblique boundary is a
weaker mixing probe than it looks; it is a fit-and-coverage probe.

The **checkerboard** in current use is not Breiman's. It is scenario 3 of
Zhu, Zeng and Kosorok (2015), reused widely since:
`X ~ N(0, Sigma)` with `Sigma_jk = 0.9^|j-k|`,
`f = 2 x5 x10 + 2 x15 x20`, `eps ~ N(0,1)`, run at n in {800, 1600} and
p in {20, 40} [verified: arXiv 2012.10737 sec 4.1]. It is a pure
two-way-interaction surface on strongly autocorrelated predictors, which
puts Tan's Theorem 5.2 pathology and a correlated design in the same cell
and gives a clean inclusion answer (columns 5, 10, 15, 20 and nothing
else). The same source verifies three more standard surfaces used with it:
van der Laan et al. (2007), `f = xt1 xt2 + xt3^2 + xt8 xt10 - xt6^2` on
`xt = 2(x - 0.5)` with uniform x and `eps ~ N(0, 0.5)`; and Meier, van de
Geer and Buhlmann (2009) 1 and 2, both **purely additive** in four
components, Meier 2's fourth component strongly oscillatory.

---

## 2. Nonparametric regression test functions

### 2.1 The Donoho-Johnstone four

The standard inhomogeneous-smoothness battery, and the reason the wavelet
literature never uses a single global smoothness assumption. All four are
one-dimensional on `t` in [0,1], at n = 2048, with Gaussian noise "rescaled
to have signal-to-noise ratio, SD(f)/sigma = 7" [verified: Donoho and
Johnstone 1994, figure legends and the appendix "Formulas for Test
Functions"]:

    Blocks     f(t) = sum_j h_j K(t - t_j),  K(t) = (1 + sgn(t))/2
               t_j = (.1,.13,.15,.23,.25,.40,.44,.65,.76,.78,.81)
               h_j = (4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,5.1,-4.2)
    Bumps      f(t) = sum_j h_j K((t - t_j)/w_j),  K(t) = (1 + |t|^4)^-1
               t_j as Blocks; h_j = |h_j| of Blocks
               w_j = (.005,.005,.006,.01,.01,.03,.01,.01,.005,.008,.005)
    HeaviSine  f(t) = 4 sin(4 pi t) - sgn(t - .3) - sgn(.72 - t)
    Doppler    f(t) = (t(1-t))^(1/2) sin(2 pi (1 + eps)/(t + eps)), eps = .05

Their relevance here is not the wavelet comparison but the shape of the
difficulty, and it is not uniform across the four:

- **Blocks** is a piecewise constant with eleven knots at irregular
  spacings, four of them within 0.05 of each other. It is the one surface
  a tree ensemble should fit exactly, and the failure mode it probes is
  therefore not bias but the sampler: getting eleven cut points right
  requires eleven correct splits and there is only one correct answer per
  knot, so this is the cleanest **cut-point placement** probe available,
  with an exactly known truth for a `varcount`-style readout.
- **Bumps** is eleven spikes with widths spanning 0.005 to 0.03 - a
  6-to-1 scale range - and is the locally-adaptive-bandwidth case.
- **HeaviSine** is smooth plus two jumps: BART's easy part and hard part
  in one function, and the natural place to look for over-smoothing at the
  jump and artifacts away from it.
- **Doppler** is a chirp: the local frequency diverges as `t -> 0`, so
  there is a region where no finite partition suffices and a region where
  a constant is right. Of the four, it is the one where a tree must
  place its splits at wildly unequal densities.

A caveat worth stating plainly: these are one-dimensional and equispaced.
Ported to dbarts they become n rows and one column, which is a regime the
package's defaults were never tuned for (m = 75 trees on one predictor).
The honest port is to embed them: one signal column carrying the DJ
function plus a realistic set of nuisance columns, which is what section
6's pathology P3 does.

### 2.2 Rotated boundaries and smooth ridges

Section 1.8's twonorm/threenorm/ringnorm are the canonical axis-alignment
stressors, with Breiman's own warning that they are variance problems and
therefore partly cured by any ensemble. The regression analogue that is
NOT cured that way is a RADIAL surface: XBART's
**Single index** function, `10 sqrt(a) + sin(5a)` with
`a = sum_{j=1..10} (x_j - gamma_j)^2`, makes the response depend on a
squared distance in ten predictors at once, so its level sets are spheres
about `gamma`, no hyperplane arrangement of finite size represents them,
and every axis-aligned split captures a vanishing share of the variation.
It is not a ridge and not a rotated coordinate: `a` is invariant under any
rotation of those ten axes, and on the independent-normal design so is the
predictors' joint law, so the cell is unchanged in distribution by every
rotation and can neither be helped by one nor discriminate one. It is the
cheapest published non-axis-aligned regression surface, it comes with
(n, p) settings and a coverage number for BART, and it sits in the same
factory as the other
three He and Hahn functions.

The other useful ridge is Tan et al.'s **Low-Dimensional Smooth**
(section 1.7): a product of two logistic ridges in ten weakly correlated
normals at a signal-to-noise ratio of 3. It is not pathological - it is
close to the average case - but it is a genuine two-way interaction with
no axis-aligned representation of finite size, and it is already
instrumented in a published R-hat study.

### 2.3 Additive with many weak effects

Two sources put this in the battery for different reasons.

Empirically, Meier, van de Geer and Buhlmann's two surfaces are purely
additive in four components on uniform predictors with `eps ~ N(0, 0.5)`
[verified: arXiv 2012.10737 sec 4.1]:

    Meier 1   f = -sin(2 xt1) + xt2^2 + xt3 - exp(xt4)
    Meier 2   f = -xt1 + (2 xt2 - 1)^2
                  + sin(2 pi xt3)/(2 - sin(2 pi xt3))
                  + 2 cos(2 pi xt4) + 4 cos^2(2 pi xt4)

with `xt = 2(x - 0.5)` and the remaining p - 4 columns pure noise; run at
n in {800, 1600}, p in {20, 40}. Meier 2's oscillatory fourth component is
the part trees handle worst.

Theoretically, Tan et al.'s Theorem 5.1 says exactly when an additive
surface becomes a MIXING problem rather than a fit problem: with `m'`
independent additive components and `m <= m'` trees, hitting an optimal
representation takes `Omega_P(n^{1/2})` sweeps, and `Omega_P(n^{qmin/2-1})`
if change and swap are removed and `m < m'`. The lever is the tree count,
and it points the opposite way from the usual advice: **more trees mixes
better here**, which Tan et al.'s Experiment 2 also finds empirically
("increasing the number of trees consistently dampens the trend in
R-hat"). Any additive-many-components cell in this battery must therefore
be run at more than one m, or it measures the tree count rather than the
kernel.

### 2.4 Regime switches keyed on a categorical

The tgp package ships the one mixed-type surface in this survey that a
consumer could call in a single line. `fried.bool` augments ten uniform
real covariates with a four-level categorical indicator, binary-encoded
into three columns (13 columns in all), and switches the mean function
entirely by level [verified: JSS 33(6) eq 1]:

    I = 1   10 sin(pi x1 x2)
    I = 2   20 (x3 - 0.5)^2
    I = 3   10 x4 + 5 x5
    I = 4   10 x1 + 5 x2 + 20 (x3 - 0.5)^2 + 10 sin(pi x4 x5)

with unit-normal noise, n = 500 train and 1000 test in the vignette's own
call. "Irrespective of I, the response depends only on {x1, ..., x5}, thus
combining nonlinear, linear, and irrelevant effects." It is a good
average-case candidate precisely because it is not cruel: mixed types,
five irrelevant columns, a categorical that genuinely modulates
everything, and a truth that supports a clean inclusion readout.

Two more one-dimensional regime-switch surfaces from the same package are
worth naming for the pathology set. The 1-d nonstationary example is
`z(x) = sin(pi x/5) + 0.2 cos(4 pi x/5)` for `x <= 9.6` and `x/10 - 1`
above it, on [0, 20] at n = 100 with `sd = 0.1` [verified: JSS 19(9)
sec 4.2] - oscillation abutting a line, with the changepoint interior.
And the 2-d exponential, `z(x) = x1 exp(-x1^2 - x2^2)` on
`[-6,6] x [-6,6]` at n = 400 [verified: JSS 33(6) sec 5.3], is a small
active bump in a large flat plane; it is the example on which Gramacy
demonstrates the failure directly, the untempered chain "almost never
visits trees of height less than five after burn-in and instead makes
rather lengthy excursions into deeper trees, exploring a local mode in
the posterior", while the tempered chain "frequently prunes back to the
tree root". That is a **tree-height trace** as the discriminating
statistic, which is cheap, needs no truth, and reads the random-walk-in-
size mode directly.

Silverman's motorcycle data belongs here too and is discussed in section 5.

---

## 3. Causal inference

bartCause is a dbarts consumer, so this section's problems are ones where
the quantity judged is the outer model's, not BART's.

### 3.1 Hill 2011's IHDP construction

The construction, as it is actually implemented in `vdorie/npci`
(`examples/ihdp_sim/data.R`, the code the maintainer wrote and the one this
survey read rather than the paywalled paper - section 8):

- Covariates are the real IHDP measurements, `covariates = "select"` giving
  25 columns - bw, b.head, preterm, birth.o, nnhealth, momage, sex, twin,
  b.marr, mom.lths, mom.hs, mom.scoll, cig, first, booze, drugs, work.dur,
  prenatal and seven site indicators - standardized. Six are continuous,
  the rest binary, and they are real-data correlated. `"full"`, `"reduced"`
  (a random 5) and `"junk"` (real plus fabricated normals) are also
  supported.
- **The imbalance is induced by deletion**: `subset(ihdp, treat != 1 |
  momwhite != 0)` drops every treated child with a nonwhite mother, so the
  treated group is a non-random subset of covariate space. That single
  line is the whole overlap pathology, and it is why IHDP is a weak-overlap
  problem rather than an ordinary one.
- `sigma.y = 1`, `tau = 4`. Setting A: `mu0 = exp((x + w) beta)` with
  `w` typically 0.5 and `mu1 = x beta - omega` - a nonlinear control
  surface against a linear treated surface, so the treatment effect is
  strongly heterogeneous by construction. Settings B and C build a design
  matrix of main effects plus quadratics and pairwise interactions and draw
  sparse coefficients from {0, 1, 2} for main effects and {0, 0.5, 1} for
  the quadratic block, with the inclusion probability shrinking in p.
  Setting C additionally re-randomizes Z from a logistic propensity on a
  random subset of those terms, shifted so the median linear predictor is
  -1.35.
- `omega` is set from the treated units under `overlap = TRUE` and from
  the CONTROL units under `overlap = FALSE`, which is the switch that turns
  the estimand into one requiring extrapolation.

For this battery IHDP's value is that it is small (well under a thousand
rows; the exact count is whatever the shipped `ihdp.RData` yields after the
deletion), real-covariate, mixed-type, weak-overlap, and reports an outer
quantity
(SATT bias, RMSE, interval coverage) rather than a fit.

### 3.2 ACIC 2016, and what made the hard cells hard

The 2016 challenge is the best-documented causal benchmark in existence and
its covariate matrix is exactly the average case this battery wants: real
data from the Collaborative Perinatal Project, "4802 observations and 58
covariates remained. Of these covariates, 3 are categorical, 5 are binary,
27 are count data, and the remaining 23 are continuous" [verified: Dorie,
Hill, Shalit, Scott and Cervone, arXiv 1707.02641v5 sec 4.2].
Response surfaces and assignment mechanisms are randomly generated
"generalized additive functions" - transformed covariates added and
multiplied, optionally passed through a link - under six knobs: degree of
nonlinearity, percentage treated, overlap, alignment, treatment effect
heterogeneity and effect magnitude. 216 knob combinations were pruned to
77 scenarios, 100 replications each, 7700 realizations.

Three findings from it bear on how to build a battery at all.

- **The hard knobs.** "two of the three that created the most difficulty
  across the board for achieving low bias were nonlinear response surfaces
  and treatment effect heterogeneity", and separately "lack of alignment
  across the assignment mechanism and the response surface emerged as one
  of the most challenging features of the data" - alignment meaning how
  much the covariates driving assignment overlap those driving the
  response.
- **Coverage is the binding metric, not bias.** "good coverage was
  difficult for most methods to achieve even when bias was low." The
  numbers: "the original BART implementation had average coverage around
  82%", raised to nearly nominal by a TMLE adjustment at about a 50%
  increase in interval length, to a little over 90% by including the
  propensity score at no length cost, and similarly by symmetric rather
  than percentile intervals.
- **Difficulty is not predictable from the data.** "Across methods the R2
  from the predictive models rarely exceeds 0.10 when predictive models
  include only non-oracle measures", and even with oracle knowledge it
  rarely exceeds 0.10 for the strong methods. A battery cannot be
  hand-selected to be hard; it has to be sampled over knobs.

The generator ships as an R package, so a cell costs a call rather than a
reimplementation.

### 3.3 ACIC 2017: targeted selection, group-correlated errors

The 2017 challenge is smaller, sharper and closer to stan4bart's own
structure [verified: arXiv 1905.09515]. Covariates are eight IHDP columns
held fixed across all 8000 datasets - mother's age, cigarettes per day
(continuous), bilirubin (continuous), four binaries, and mother's birth
place as a **16-level categorical** - at n = 4302, with pairwise
correlations at most 0.20. The generator is fully specified:

    f(x)   = x1 + x43 + 0.3(x10 - 1)
    pi(x)  = (1 + exp(kappa1 f(x) + kappa2))^-1
    mu(x)  = -sin(Phi(pi(x))) + x43
    tau(x) = xi (x3 x24 + (x14 - 1) - (x15 - 1))
    sigma(x) = 0.4 + (x21 - 1)/15

with `xi` in {1/3, 2} (effect magnitude), `eta` in {1/4, 5/4} (noise, via
`sigma_y = eta sqrt(Var(mu + pi tau))`), and `(kappa1, kappa2)` in
{(0.5, 0), (3, -1)} (selection strength). Four error types cross those
eight cells to give 32 DGPs at 250 replicates each.

The reason to single this one out: **mu is a function of pi**, which is
targeted selection stated as an equation, and one of the four error types
is `sigma_y (0.9 eps_i + 0.1 eps_{x21})` - "10% of the error term is shared
among variables with common values of x21" - which is a 16-group random
intercept written into the DGP. That is stan4bart's own model shape with a
known truth, and the outer quantity is the CATE. Note the paper's own
erratum: the heteroskedastic-error results are declared incorrect and
withdrawn, so use the iid and group-correlated arms.

Hahn, Murray and Carvalho's own reading of the 2016-to-2017 change is worth
carrying: the 2016 sets had "unrealistically large average treatment
effects and similarly unrealistic degrees of heterogeneity", with SATT
interquartile range 0.57 to 0.79 standard deviations of Y, and "the 2017
competition explicitly incorporated targeted selection (unlike the 2016
datasets)" [verified: arXiv 1706.09523v4 sec 6.2].

### 3.4 ACIC 2019 and weak overlap

The 2019 challenge targets the population ATE across two tracks, low
dimensional at roughly 500 x 20 and high dimensional at roughly
1000 x 200 and 2000 x 200, 3200 datasets per track from 32 DGPs, covariates
"from publicly available data and also simulated", scored on the point
estimate and a 95% confidence interval. Its stated built-in difficulties
are "non-linearity of the response surface, treatment effect heterogeneity,
varying proportion of true confounders among the observed covariates, and
near violations of the positivity assumption" [verified: the challenge
site]. The last of those is the explicit weak-overlap arm this battery
wants, and the 1000 x 200 / 2000 x 200 track is the only causal cell in
this survey with p in the hundreds.

### 3.5 Hahn, Murray and Carvalho 2020: targeted selection as a design

Two constructions, both cheap.

**The two-dimensional shelf** (their Example 1): `d = 2`, `n = 250`,
`x1, x2 ~ U(0,1)`, homogeneous effect `tau = -1`, `eps ~ N(0,1)`,
propensity
`pi = 0.8 Phi(mu(x1,x2) / (0.1(2 - x1 - x2) + 0.25)) + 0.025(x1 + x2) + 0.05`.
The prognostic `mu` is described rather than printed: it has a "shelf" at
the line `x1 = x2` and ranges from -3 to 3. Over 200 datasets, standard
BART records **bias 0.27, 95% coverage 65%, RMSE 0.31**, against BCF's
0.14, 95% and 0.21 [verified: arXiv 1706.09523v4 Table 1]. The mechanism
is stated exactly: "it takes many axis-aligned splits to approximate the
'shelf' across the diagonal ... At the same time, due to the strong
confounding in this example a single split in Z can stand in for many
splits on x1 and x2 that would be required to approximate mu(x)." A
diagonal boundary and regularization-induced confounding in one 2-column,
250-row cell, with a published BART coverage number to reproduce. The
authors label it "somewhat stylized in that we designed it specifically to
be difficult to learn for tree-based models", which is exactly what a
pathology is for.

**The eight-way factorial** (their section 6.1): five covariates, "the
first three are continuous, drawn as standard normal random variables, the
fourth is a dichotomous variable and the fifth is unordered categorical,
taking three levels"; `tau(x)` either 3 or `1 + 2 x2 x5`; `mu(x)` either
`1 + g(x4) + x1 x3` or `-6 + g(x4) + 6|x3 - 1|` with
`g = (2, -1, -4)`; propensity
`pi(x) = 0.8 Phi(3 mu(x)/s - 0.5 x1) + 0.05 + u/10` with `s` the sample
standard deviation of `mu`; n in {250, 500}; 200 replications; scored by
RMSE, coverage and average interval length for both ATE and CATE
[verified: arXiv 1706.09523v4 sec 6.1]. Mixed covariate types, small n,
targeted selection, and a homogeneous-versus-heterogeneous toggle. It is
the closest thing in the literature to an average-case causal cell that
still has a known truth, and dbarts already ships the BCF sampler it was
written for.

---

## 4. The parametric MCMC canon, and what transfers

For each: does the pathology appear in a sum-of-trees posterior, and can a
BART analogue be built.

### 4.1 Neal's funnel - YES, and dbarts already has three of them

`v ~ N(0, 3^2)`, and `x1..x9 | v ~ N(0, e^v)` independent [verified: Neal,
"Slice sampling", arXiv physics/0009028 sec 8]. Neal's own framing is the
reason it is in this document: "Such a distribution is typical of priors
for components of Bayesian hierarchical models - x1 to x9 might, for
example, be random effects for nine subjects, with v being the log of the
variance of these random effects." Multivariate Metropolis on it gives
results that "are grossly incorrect... Moreover, there is little in the
plot to indicate that anything is wrong", and Neal notes the remedy that
would have caught it: "Running several chains from different starting
states might have revealed the problem".

**Appears in a sum-of-trees posterior, three ways.** (i) A grouped
random-intercept scale over the forest - which is stan4bart's tau, and
which is why that scale is the outer quantity to report. (ii) The
per-forest amplitude of the multiplier combiner: the response is a scale
times a sum of leaf values, and both are sampled, which is a level/scale
funnel by construction. This package has already measured its narrow end:
at `a0 = 40` and `100` the sampler sits with "sigma plateaus ~5x high with
NO decay through 40k sweeps - frozen structure", and injecting a large `a`
made the bias WORSE mid-burn than a cold start
([3.2 Structure freezes when the noise level is low (ESTABLISHED)](tree-mixing-proposals.md#32-structure-freezes-when-the-noise-level-is-low-established)).
(iii) The leaf-prior scale itself, which is the same object at m = 1.
**A BART analogue is cheap to build**: a grouped design with few
observations per group and a genuinely small group-level scale, fitted with
the group scale sampled, reporting the scale's own ESS and a
prior-predictive calibration check rather than any fit statistic.

### 4.2 Eight schools - YES, as the parameterization pair

`J = 8`, `y = (28, 8, -3, 7, -1, 1, 18, 12)`,
`sigma = (15, 10, 16, 11, 9, 11, 10, 18)`; the centered model puts
`theta ~ N(mu, tau)` with `tau ~ Cauchy(0, 5)` and the non-centered model
reparameterizes it [verified: posteriordb, `eight_schools_centered.stan`
and the shipped `eight_schools.json` data]. Its whole content as a
benchmark is that two algebraically identical posteriors mix completely
differently, so it is a **parameterization** test, not a surface.

**The transfer is a warning rather than a problem to add.** BART's own
centered/non-centered fork is the forest-plus-random-effect
interweaving question, and this package already investigated it and
recorded NO-GO (`docs/design/forest-ranef-interweaving.md`), and grouped
random intercepts are being retired from dbarts entirely
(`docs/design/retire-grouped-random-effects.md`). So the eight-schools
analogue belongs to stan4bart, not to dbarts, and the thing dbarts owes it
is that any embedded-use cell in this battery reports the OUTER scale
parameter's effective sample size.

### 4.3 Banana and Rosenbrock ridges - YES, but already measured here

A curved narrow valley in a continuous parameter space. The primary
sources for the twisted-Gaussian construction could **not be fetched**
(section 8), so nothing numeric is claimed for it here.

**The sum-of-trees analogue is not in tree space at all**: it is the ridge
between a parametric block and the forest when both can explain the same
signal, and this house has it at 6x
([4.1 Move signal out of the forest, and let BART fit the remainder](tree-mixing-proposals.md#41-move-signal-out-of-the-forest-and-let-bart-fit-the-remainder),
sourced from `forest-ranef-interweaving.md`). Adding a synthetic banana
would measure nothing that a composed-model cell does not measure better,
because the composed cell has the ridge AND the right estimand.

### 4.4 Multimodal mixtures and label switching - YES, exactly

The `k!` symmetry of a mixture posterior under relabelling, and the family
of remedies - identifiability constraints, relabelling algorithms, and
**label-invariant loss functions** (Jasra, Holmes and Stephens, Statistical
Science 20(1); record and abstract only, full text **not fetched**).

**This is the closest structural match in the canon.** The sum-of-trees
map from ensembles to functions is massively many-to-one, and the modes are
exchangeable in the same way: permute tree labels, or split one main effect
across two trees.
[3.1 Many tree arrangements, one fitted function (ESTABLISHED)](tree-mixing-proposals.md#31-many-tree-arrangements-one-fitted-function-established)
records the measurement (between-chain standard deviation of the
root-on-x1 fraction 0.3619 against a mixing null near 0.05).

The methodological transfer is the one worth taking: **the mixture
literature's answer is not to fix the sampler but to report only
label-invariant functionals**. Applied here that says the fitted function,
sigma and any outer-model estimand are legitimate one-chain readouts;
`varcount`, interaction counts and `plotTree` are not, and must either be
symmetrized or read between chains. Every structural statistic in section
6 is specified that way.

### 4.5 Correlated-design logistic regression - PARTLY

The canonical cell is Girolami and Calderhead's collection of five
logistic-regression datasets, "where n ranges between 250 and 1000 and d
ranges between 7 and 25", German credit at `d = 25, n = 1000` the largest,
scored by "the minimum ESS (across the d dimensions of pi) per second
CPU-time", with ESS from Geyer's initial monotone sequence estimator
[verified through Kleppe, arXiv 1501.07454 sec 5, which replicates their
experiment; the Girolami and Calderhead paper itself was **not fetched**].

**The pathology - a strongly correlated Gaussian-ish posterior in a
moderate number of coefficients - does not transfer**: BART has no
coefficient vector to be correlated. What transfers is the **metric**.
"Minimum ESS across dimensions per unit of run time" is the right shape for
this battery's ranking column, and it is stricter than a median: it is the
worst coordinate, which is where a sticky sampler shows.

The correlated DESIGN does transfer, separately, and is already covered:
the He and Hahn factor-structure predictor matrix (section 1.4), the
`0.9^|j-k|` checkerboard (section 1.8), and ACIC's real covariate matrices
(section 3).

### 4.6 Spike-and-slab and horseshoe model-space mixing - YES, measured

The horseshoe's difficulty is documented as exactly the funnel of section
4.1: "the problem arises due to posterior having an extreme funnel shape
which is challenging for Markov chain Monte Carlo (MCMC) methods. The
problem was revealed with the help of the divergence diagnostics of the
NUTS algorithm... when fitting the models in Stan." The fix - half-t with
small `nu` in place of half-Cauchy - works, but "the drawback is that the
prior becomes less sparsifying" [verified: Piironen and Vehtari, arXiv
1707.01694 sec 2.4].

**dbarts' analogue is DART, and the trade has already been measured in the
same direction.** Tan et al.'s Experiment 3 compares a Dirichlet split
prior at `alpha = 1` against the uniform one: "Use of a Dirichlet prior on
split feature probabilities either exacerbates the increasing trend for
R-hat, or at best has ambiguous effect. This is despite improving
prediction performance on datasets exhibiting sparsity such as Low
Dimensional Smooth and Piecewise Linear" [verified: arXiv 2406.19958
appendix L.2]. A sparsity prior buys accuracy on a sparse surface and
costs mixing, in trees as in coefficients. That is a live instance of this
battery's decision rule pointing at a shipped feature, and it is the
strongest argument in this survey for measuring `sparse = TRUE` cells
separately rather than folding them into a default cell.

### 4.7 posteriordb as a curated source

148 posteriors, of which 57 carry reference draws (a gold-standard sample
to compare against), each pairing a Stan model with a JSON data set
[verified: the `stan-dev/posteriordb` repository tree]. It contains
eight_schools in both parameterizations and two motorcycle-data posteriors
(`mcycle_gp-accel_gp`, `mcycle_splines-accel_splines`); it contains no
funnel, no German credit, no horseshoe and no Rosenbrock under those names,
and no tree model at all.

**Its value here is the convention, not the problems.** A reference-draw
set - a long, trusted run stored alongside the problem, against which a
candidate sampler is scored - is what would let this battery answer
"is dbarts' posterior right" rather than only "did it move". dbarts has one
already in a different currency: `benchmarks/baselines` holds bitwise
equivalence baselines, which answer "did the draws change" and cannot
answer "are they correct". The SBC harness (`benchmarks/R/sbc.R`) is the
correctness gate this package does own, and it is the natural place to
attach a funnel-shaped prior tail.

---

## 5. Real datasets with known hardness for trees

Real data has no truth, so it can only carry between-chain statistics,
held-out error and coverage of the posterior PREDICTIVE. That is a real
restriction, and it is why sections 1 to 4 carry the weight. Three
clusters are worth having anyway.

**The four PMLB sets, because two groups have already published dbarts'
failure on them.** Breast Tumor (n = 116640, p = 9), California Housing
(20640, 8), Echo Months (17496, 9), Satellite Image (6435, 36), each
subsampled to n = 200 and 2000 as well as full, with a fixed 10% held-out
set. Ronen et al. report BART at m = 200 failing Gelman-Rubin on all four
at full size; Tan et al. reuse the same four and report R-hat rising with
n on every one [verified: arXiv 2210.09352 sec 4.1, arXiv 2406.19958 sec
9.1.3]. They are the only real datasets in this survey with a published,
reproducible dbarts mixing failure attached.

One warning that comes with them. Tan et al.'s Experiment 4 explains why
R-hat rises with longer burn-in: "having a longer burn-in leaves the
between-chain differences largely unchanged because each individual chain
does not escape its local mode. On the other hand, a longer burn-in removes
early, highly variable samples, thereby reducing within-chain variance. As
a result, R-hat increases" [verified: arXiv 2406.19958 appendix L.3].
R-hat can therefore worsen for a reason that is not a regression, so a
battery that gates on R-hat must gate on R-hat AT A FIXED SWEEP BUDGET and
never compare arms that differ in burn-in.

**The bake-off sets, for the average case and for the cases trees lose.**
The 42-set CGM bake-off (section 1.1) and SoftBart's ten (section 1.3)
overlap heavily because both draw on Kim et al. (2007). Three specific
losses are worth keeping as named cells:

- **tecator**: BART-CV 1.87 and DART 1.63 against SBART 0.98, with the
  authors' own reading that smoothness is essential here [verified: arXiv
  1707.09461 Table 1]. The strongest published "smooth beats trees" case.
- **ozone**, **triazine** and **pole**: random forests beat both BayesTree
  and bartMachine [verified: JSS 70(4) Table 3, values 4.068 vs 4.105 and
  4.129; 0.119 vs 0.130 and 0.130; 10.699 vs 11.731 and 12.764]. These are
  small margins on small data and are worth having as a "BART is not
  uniformly best" reminder rather than as a diagnostic.
- **boston**: bartMachine 3.003 against random forests 4.581, a
  significant win, and the set Gramacy and Lee also used for treed GPs. It
  belongs in the average-case group, not the hard one.

**Gramacy's nonstationary examples**, for a failure mode nothing else in
this survey reads directly. The motorcycle data (Silverman 1985) is
one predictor, and "Many authors have commented on the existence of
two - perhaps three - regimes in the data over time where the
characteristics of the mean process and noise level change (i.e., a
nonstationarity and heteroskedasticity, respectively)" [verified: JSS 19(9)
sec 1.1]. dbarts has a heteroscedastic variance forest, so motorcycle is
a natural end-to-end check of it on real data, and posteriordb carries two
independent posteriors on the same data (`mcycle_gp-accel_gp`,
`mcycle_splines-accel_splines`) to compare a fit against.

The 2-d exponential (section 2.4) is the one with the diagnostic attached:
the untempered chain "almost never visits trees of height less than five
after burn-in and instead makes rather lengthy excursions into deeper
trees", while a tempered chain "frequently prunes back to the tree root".
Tree height as a trace is the cheapest readout of section 3.5's
random-walk-in-size mode in the mixing survey, needs no truth, and is
already recoverable from a dbarts fit.

The LGBB rocket-booster computer experiment (Gramacy, Lee and Macready) is
named in the tgp papers but its data were **not fetched** here.

---

## 6. The battery

### 6.1 The rule, stated operationally

A kernel or prior change is accepted only if it is **neutral or better on
every average-case core cell** and **better on at least one pathology**.
Never the reverse. Two consequences that are easy to lose:

- The core is a gate, not a score. A change that improves the core is not
  thereby accepted; a change that regresses one core cell is refused even
  if it wins every pathology.
- Winning a pathology means winning THAT pathology's pre-registered
  primary statistic, chosen before the run, not whichever of its metrics
  happened to move.

House law from `grow-from-root-default.md` applies unchanged: per-cell
checks and never a pooled aggregate; thresholds frozen against a pilot
before any confirmatory contrast; a mandatory fresh-seed re-run of any
single flagged cell; a null control whose failure voids the estimator
family. Every cell below is paired on matched seeds in the
`benchmarks/R/grouped-mixing.R` idiom - data seed indexed by cell and
replicate, sampler seed the replicate, both shared across arms.

Cost is quoted in **unit fits**, where one unit fit is the n = 2000,
p = 10, m = 200, 3000-sweep, single-threaded fit that section 13 of the
mixing survey measured at 4.2 to 6.7 seconds on arm64 macOS. Costs below
are for one confirmatory run of one contrast (two arms), not for the pilot.

### 6.2 The average-case core

| id | problem | source | truth | n, p, family | stresses | dbarts use | primary statistic | cost |
|---|---|---|---|---|---|---|---|---|
| C1 | He and Hahn factorial, correlated-factor arm, Trig+poly and Single index | arXiv 2002.03375v4 sec 4.1 | yes | n 10000, p 30, and n 1000, p 100; gaussian, kappa 1 | correlated continuous design, one true interaction, a rotated ridge, realistic n | one-shot | 95% pointwise coverage of true f on 1000 held-out rows | ~120 unit fits |
| C2 | `fried.bool` mixed-type regime switch | tgp, JSS 33(6) eq 1 | yes | n 500 train / 1000 test, p 13 (10 real + 3 dummy for a 4-level factor); gaussian | mixed types, categorical modulation, 5 irrelevant columns | one-shot | coverage of true f; summed inclusion share on x1..x5 read between chains | ~10 unit fits |
| C3 | ACIC 2016 moderate cells | arXiv 1707.02641; `aciccomp2016` | yes | n 4802, p 58 (3 categorical, 5 binary, 27 count, 23 continuous); gaussian | real correlated mixed-type covariates, generated additive-plus-interaction surfaces, confounding | causal (bartCause) | SATT bias and 95% interval coverage | ~500 unit fits |
| C4 | ACIC 2017 group-correlated arm, weak selection, low effect | arXiv 1905.09515 sec 4.1 | yes | n 4302, p 8 (3 continuous, 4 binary, one 16-level factor); gaussian, 16-group shared error | embedded use: an outer block moves the response between sweeps; targeted selection | embedded, moving response; causal | CATE RMSE and coverage; ESS of the outer group-scale parameter | ~250 unit fits |

**C1** is the only core cell with a published BART number to reproduce, and
it is the number that most changes the picture: 95% pointwise coverage of
0.74 to 0.78 at n = 10000 and moderate noise (section 1.4). Run the
correlated-factor predictor arm, not the independent one, because the
independent arm is not the average case for a residual surface. Two mean
functions are enough - Trig+poly for the interaction, Single index for the
rotated ridge; Linear and Max belong to the pathology set if anywhere.
Report coverage on a fixed 1000-row held-out set drawn once per cell, and
carry minimum ESS across 25 fixed points as the secondary.

**C2** is the cheap mixed-type cell and the only one in the battery where a
categorical genuinely changes the function rather than shifting it. It is
also the natural inclusion cell for a core gate, because the truth
(x1..x5 signal, x6..x10 irrelevant, the factor relevant) is exact. Read
inclusion BETWEEN chains, per section 4.4: at m = 75 the within-chain
average self-averages and detects nothing.

**C3** is the average case as a consumer meets it - real covariates, real
correlation, real mixed types, p = 58 - and it is judged on what bartCause
reports. Do not use the whole 77-scenario grid as a gate; select a
moderate stratum once (medium nonlinearity, overlap on, high alignment,
heterogeneous effects) and freeze the scenario list. The reason to gate on
coverage and not bias is section 3.2's own finding: BART's ACIC-2016
coverage was about 82%, so there is room to lose and room to gain, whereas
its bias was already near the floor.

**C4** is the cell that exists because dbarts' distinguishing use is a
sampler inside an outer loop, and no published measurement covers it. The
group-correlated error arm is a 16-group random intercept written into the
DGP, so the composed model is correctly specified and the outer scale has a
known truth. This is the cell that inherits the response-swap recovery
question from
[14. Recovery after a response swap (2026-09-06)](tree-mixing-proposals.md#14-recovery-after-a-response-swap-2026-09-06):
alongside the estimand, record recovery sweeps after the first few
`setResponse` calls, so the battery measures the moving-response regime and
not just the fixed one. Its honest caveat is that it needs stan4bart or an
equivalent outer loop, so it is the one core cell that is not
self-contained in this repository.

### 6.3 The pathologies

Each covers a mode the others do not.

| id | problem | source | truth | n, p, family | failure mode | dbarts use | primary statistic | cost |
|---|---|---|---|---|---|---|---|---|
| P1 | Low-noise Friedman | arXiv 1312.1895 sec 2.2 | yes | n 5000, p 10, m 200, sigma^2 0.1 | structure freezes as sigma falls; acceptance collapse | one-shot | 90% pointwise coverage of true f; per-move acceptance rate | ~80 unit fits |
| P2 | Confounded predictors (Wu, Tjelmeland and West step function) | arXiv 1312.1895 sec 2.3 | yes, plus an exact symmetry | n 300, p 3, m 1, 8 chains | representation multimodality with two exactly equiprobable modes | one-shot, structural | between-chain sd of the fraction of draws with the root on x1, against the 0.5 symmetry null | ~3 unit fits |
| P3 | PMLB R-hat ladder (California Housing, Echo Months) | arXiv 2210.09352 sec 4.1 | no | n 200 / 2000 / full, p 8 and 9; 8 chains, m 200 | mixing degrades with n on real data | one-shot | R-hat of held-out RMSE across 8 chains at a FIXED sweep budget, and its slope in n | ~4000 unit fits |
| P4 | Doppler embedded in nuisance columns | Donoho and Johnstone 1994 | yes | n 2048, p 10 (1 signal + 9 noise), SD(f)/sigma 7 | inhomogeneous smoothness: required cut density varies by orders of magnitude | one-shot | pointwise coverage split by region (high-frequency third vs flat third) | ~10 unit fits |
| P5 | Checkerboard on an autocorrelated design | Zhu, Zeng and Kosorok 2015 scenario 3, via arXiv 2012.10737 | yes | n 1600, p 40, Sigma_jk 0.9^abs(j-k) | pure two-way interaction plus a design where inclusion is ambiguous | one-shot, structural | between-chain sd of inclusion on {x5,x10,x15,x20} and their immediate neighbours | ~50 unit fits |
| P6 | The diagonal shelf with targeted selection | arXiv 1706.09523v4 sec 4, Table 1 | yes | n 250, d 2 plus Z; 200 replications | rotated boundary plus regularization-induced confounding | causal | ATE bias and 95% interval coverage (published BART: 0.27 and 65%) | ~40 unit fits |
| P7 | IHDP with overlap off | `vdorie/npci` `examples/ihdp_sim/data.R`, setting A, `overlap = FALSE` | yes | n under a thousand after the deletion rule, p 25 mixed; tau 4, sigma 1 | extrapolation into a region with no treated support | causal | SATT bias and 95% coverage; interval width inside the unsupported region | ~20 unit fits |
| P8 | Hierarchical variance funnel | Neal, arXiv physics/0009028 sec 8, realized on ACIC 2017's group-correlated arm at low noise | yes | n 4302, 16 groups; gaussian | Neal's funnel in the outer scale, at the narrow end | embedded | SBC rank uniformity and ESS of the group scale; NOT any fit statistic | ~500 unit fits |

**P1** is the failure this package has already reproduced twice, and it
must stay in the battery as the known-positive control. Run it on two
rungs: Pratola's own cell (n = 5000, m = 200, `sigma^2 = 0.1`, published
90% coverage 53% under birth/death only) and section 13's cheaper cell
(n = 2000, m = 200, `sigma = 0.25`, measured here at 90% coverage 0.71
under all three mixtures section 13 ran). If the second rung does not come back
near 0.71 in the control arm, the harness is broken and no verdict from any
other cell is valid. Its weakness as a discriminator is established: three
proposal mixtures were indistinguishable on it.

**P2** is the highest-value pathology per second in this survey. Three
columns, three hundred rows, one tree, a published acceptance rate of zero,
and - because x1 and x3 are confounded by construction - an exact 0.5
symmetry that serves as its own oracle. It reads representation
multimodality directly, which no core cell can. It needs the duplicate-
column null control alongside it (`x1 == x2`, where both the likelihood and
prior ratios are exactly 1, so both arms must return pooled p1 within Monte
Carlo error of 1/2 with non-zero switch counts in every chain).

**P3** is the only real-data cell with a published dbarts mixing failure,
which makes it the battery's external credibility, and the most expensive
thing in it. Restrict it to California Housing and Echo Months; Breast
Tumor at n = 116640 costs roughly six times as much for the same finding.
Gate at a FIXED sweep budget, never comparing arms with different burn-in,
for the reason section 5 gives.

**P4** ports the wavelet literature's inhomogeneity into dbarts' own shape:
Doppler on one column with nine nuisance columns, so the fit is a realistic
p = 10 problem and the failure is local. Split the coverage readout by
region - the chirp's high-frequency third against the flat third - because
a pooled coverage number averages the two failures away. The natural
secondary is the empirical distribution of cut points along the signal
column against the local frequency, which is a structural readout and must
be read between chains.

**P5** is the pure-interaction mode that Tan et al.'s Theorem 5.2 names,
placed on the correlated design that makes it a variable-selection problem
too. It is the one cell where inclusion has an exactly right answer
({x5, x10, x15, x20}) AND a set of near-decoys (x4, x6, x9, x11, ... at
correlation 0.9 and 0.81 with the true columns). Their theorem's condition
is `pi_c = 0`, which the default is not, so this cell tests whether the
change move actually rescues the pure interaction that the theory says
grow/prune cannot reach.

**P6** is the cheapest cell in the battery with a published BART failure in
an OUTER quantity: bias 0.27 and 65% coverage against BCF's 0.14 and 95%,
over 200 replications of a 250-row, 2-column problem. It is the only cell
that isolates a rotated boundary from everything else, and dbarts ships
both arms of the published comparison, so the reproduction is a script.

**P7** is extrapolation, which nothing else here covers: `overlap = FALSE`
changes `omega` so that the estimand is defined on treated units whose
covariate region has no controls, on top of the deletion rule that already
removed the treated nonwhite-mother stratum. The right readouts are
whether the interval widens where it should and whether coverage survives,
not RMSE.

**P8** is the funnel, realized on a DGP that already contains it rather
than as a synthetic. It is the only cell judged purely on the outer scale
parameter, and the only one where SBC rank uniformity - not coverage, not
error - is the primary. It is expensive and it depends on stan4bart, so it
is last to build; but it is the only cell that reads the pathology
stan4bart's users actually hit.

**Modes NOT separately covered, deliberately.** Additive-with-many-weak-
effects (section 2.3) is a tree-count question rather than a kernel
question, so it belongs in a `m` sweep attached to C1 rather than in its
own cell. Oblique classification boundaries (twonorm/threenorm/ringnorm)
are subsumed by P6, which carries the battery's one genuinely LINEAR
boundary - the shelf at the line `x1 = x2` - with Breiman's own caveat that
an ensemble already fixes most of them. The Single index arm of C1 does not
belong in that list: it is radial, hence rotation-invariant on the
independent design (section 2.2), so it neither tests nor rewards an
oblique rule. P7's setting A is the battery's other oblique structure, both
of its surfaces being functions of one linear index. The banana ridge is
subsumed by C4, which has the ridge and an estimand.

### 6.4 What "no regression on the core" means numerically

Twenty matched pairs per cell. All checks per cell, Holm-corrected across
cells within a metric, with a mandatory fresh-seed re-run of any single
flagged cell before a flag counts.

The margins are set from the resolution section 13 actually achieved at
five pairs, scaled to twenty: paired standard error about 0.012 on RMSE,
0.02 on coverage, 6 on median f ESS and 0.01 on summed inclusion share at
five pairs, so about 0.006, 0.010, 3 and 0.005 at twenty. A regression is
flagged when the paired mean difference is worse than the margin AND its
one-sided 95% bound excludes the null; both conditions, so that a noisy
cell cannot flag on its point estimate alone.

    metric                                  margin (worse than this fails)
    95% pointwise coverage of true f        -0.010 absolute
    held-out RMSE against true f            ratio > 1.02
    minimum ESS over 25 fixed points,       ratio < 0.90
      per second
    summed inclusion share on true columns  -0.010 absolute
    outer estimand RMSE (C3, C4)            ratio > 1.02
    outer estimand interval coverage        -0.010 absolute
    ESS of the outer scale (C4)             ratio < 0.90
    wall time per sweep                     ratio > 1.05

Two absolute gates on top of the paired ones. First, **the null control**:
any arm that is meant to be inert - a new kernel at weight zero, a new
prior at its identity setting - must be BITWISE identical to the control
under `benchmarks/R/equivalence.R`, not merely statistically equal. A
statistical pass with a bitwise failure means the change is not what it
claims to be. Second, **P1 must still fail**: if the
n = 2000, `sigma = 0.25` rung's 90% coverage does not sit near 0.71 in the
control arm, the harness is mismeasuring and no verdict is valid (measured
0.725 held-out at the shipped default,
[10.8 P1, the low-noise Friedman emulator (2026-09-07)](#108-p1-the-low-noise-friedman-emulator-2026-09-07)).

"Better on at least one pathology" is the mirror image, and is
deliberately harder than the no-regression bar: a paired improvement on
that pathology's pre-registered primary statistic exceeding **four times**
the measured per-replicate standard error, per cell, surviving a fresh-seed
re-run. Four times, not two, because
[6.4 Kill criteria, pre-registered](tree-mixing-proposals.md#64-kill-criteria-pre-registered)
already fixed that margin for this package and there is no reason to
loosen it.

### 6.5 Ranking, and the three to run first

Ranked by discriminating power per unit of run time - published failure to
reproduce, exactness of the oracle, and distinctness of the mode, divided
by cost.

    rank  cell  cost (unit fits)  why it ranks here
    1     P2       3    published acceptance rate of 0, an exact symmetry
                        oracle, three columns, and the only direct read on
                        representation multimodality
    2     P6      40    published BART bias 0.27 / coverage 65% against a
                        shipped alternative at 0.14 / 95%, in an outer
                        quantity, on 250 rows
    3     P5      50    exact inclusion truth with near-decoys at
                        correlation 0.9; tests the one theorem whose
                        premise dbarts' default violates
    4     P4      10    cheap, exact truth, and the only local-adaptivity
                        read; ranks below P5 only because the cut-density
                        readout needs new tooling
    5     P7      20    cheap and consumer-relevant, but its statistic
                        (interval width off-support) is soft
    6     P1      80    known-positive control, but three proposal mixtures
                        were indistinguishable on it
    7     C1     120    the largest published coverage deficit in this
                        survey (0.74-0.78 at 95% nominal), but it is a
                        core gate rather than a discriminator
    8     C2      10    cheap; low power, since it is not a hard problem
    9     C4     250    the only embedded cell, but it needs an outer loop
                        that lives in another repository
    10    C3     500    the truest average case and the slowest gate
    11    P8     500    the right pathology for stan4bart, but SBC replicate
                        counts dominate the cost
    12    P3    4000    the strongest external anchor and by far the most
                        expensive; R-hat also carries the burn-in artifact

**Run first: P2, P6, P5.** Together they cost under 100 unit fits - a
single afternoon - and they cover three disjoint modes (representation
multimodality, rotated boundary in an outer estimand, pure interaction with
ambiguous inclusion). Two of the three have published numbers that a
correct harness must reproduce before any arm contrast is believable, which
makes them the pilot as well as the first measurement.

**Build the core in parallel but gate on C1 first.** C1 is the only core
cell that carries a published BART number, and reproducing 0.74 to 0.78
coverage at n = 10000 would establish, in this house, that the coverage
deficit is not a low-noise curiosity. If dbarts does NOT reproduce it, that
is the more interesting result and it changes what the rest of the battery
is for.

### 6.6 What this battery does not measure

- **Non-gaussian families.** Every cell above is gaussian or a causal
  gaussian. dbarts ships nine response models; the coverage and inclusion
  statistics carry over, but no cell here exercises them, and section 13's
  one probit design is the only published hint (change helped the worst
  point's ESS there, swap did not).
- **p in the thousands.** The largest p here is 200 (ACIC 2019's high
  track, not adopted above) and 58 in C3. CGM's own p = 1000 at n = 100
  cell is not represented, because a residual surface with p >> n is not
  the average case this battery is built for.
- **Anything at m = 1 except P2.** That is deliberate - the package's own
  finding is that structural statistics need one tree to be visible and
  that one tree is not the shipped configuration, so P2 quarantines that
  regime rather than spreading it.
- **Correctness.** Coverage of the true mean is a frequentist calibration
  statistic, not a proof the sampler targets the right posterior; only
  `benchmarks/R/sbc.R` answers that, and only P8 uses it.

---

## 7. References

Every row was fetched and read in this arc unless the "read" column says
otherwise. "Full text" means the PDF or source was downloaded and the
cited passage read directly, not summarized.

| # | source | read | url |
|---|---|---|---|
| 1 | Chipman, George, McCulloch, "BART: Bayesian additive regression trees", AOAS 4(1):266-298, 2010 | full text: sec 3.1, 3.2, 5.1, 5.2, 5.3, 6 | arxiv.org/abs/0806.3286 |
| 2 | Chipman, George, McCulloch, "Bayesian CART model search", JASA 93(443):935-948, 1998 | NOT FETCHED - the reachable copy is a scanned JSTOR image with no text layer | jstor.org/stable/2669832 |
| 3 | Pratola, "Efficient Metropolis-Hastings proposal mechanisms for Bayesian regression tree models", Bayesian Analysis 11(3):885-911, 2016 | full text: sec 2.2, 2.3, 5.1, 5.2, 6 | arxiv.org/abs/1312.1895 |
| 4 | Linero, "Bayesian regression trees for high-dimensional prediction and variable selection", JASA 113(522):626-636, 2018 | NOT FETCHED - paywalled, no preprint located | doi.org/10.1080/01621459.2016.1264957 |
| 5 | Linero, Yang, "Bayesian regression tree ensembles that adapt to smoothness and sparsity", JRSS-B 80(5):1087-1110, 2018 | full text: sec 4.1, 4.2, 4.3 and Table 1 | arxiv.org/abs/1707.09461 |
| 6 | He, Yalov, Hahn, "XBART: accelerated Bayesian additive regression trees", AISTATS 2019 | full text: sec 4.1 (Table 1), 4.2, 4.4, 4.5 | arxiv.org/abs/1810.02215 |
| 7 | He, Hahn, "Stochastic tree ensembles for regularized nonlinear regression", JASA 2023 | full text: sec 4.1 (DGP factory, Tables 1-2), sec 5 (Table 4) | arxiv.org/abs/2002.03375 |
| 8 | Kapelner, Bleich, "bartMachine: machine learning with Bayesian additive regression trees", JSS 70(4), 2016 | full text: sec 4.3, 4.10, 4.11, Figure 4, appendix B (Tables 3-4) | doi.org/10.18637/jss.v070.i04 |
| 9 | Hill, Linero, Murray, "Bayesian additive regression trees: a review and look forward", Annu. Rev. Stat. Appl. 7:251-278, 2020 | NOT FETCHED - the SSRN copy returns 403, no open copy located | doi.org/10.1146/annurev-statistics-031219-041110 |
| 10 | Ronen, Saarinen, Tan, Duncan, Yu, "A mixing time lower bound for a simplified version of BART", 2022 | full text: sec 1.2, 4.1 (Table 2, Figures 2-3), 4.2, 5 | arxiv.org/abs/2210.09352 |
| 11 | Tan, Ronen, Saarinen, Yu, "On the computational efficiency of Bayesian additive regression trees: an asymptotic analysis", 2024 | full text: sec 5.1, 5.2, 9.1-9.2, appendices L.2, L.3, L.6 | arxiv.org/abs/2406.19958 |
| 12 | Breiman, "Bias, variance, and arcing classifiers", Berkeley Statistics TR 460, 1996 | full text: sec 2.3 and Table 1 | stat.berkeley.edu/~breiman/arcall96.pdf |
| 13 | Feng, Baumgartner, "(Decision and regression) tree ensemble based kernels", 2020 - used here as the verifying secondary for the Checkerboard (Zhu, Zeng, Kosorok 2015 scenario 3), van der Laan et al. 2007 and Meier et al. 2009 surfaces | full text: sec 4.1 | arxiv.org/abs/2012.10737 |
| 14 | Donoho, Johnstone, "Ideal spatial adaptation by wavelet shrinkage", Biometrika 81(3):425-455, 1994 | full text: sec 2, figure legends, appendix "Formulas for Test Functions", Tables 2-3 | imjohnstone.su.domains/WEBLIST/1994/isaws.pdf |
| 15 | Gramacy, "tgp: an R package for Bayesian nonstationary, semiparametric nonlinear regression and design by treed Gaussian process models", JSS 19(9), 2007 | full text: sec 1.1, 4.1, 4.2, 4.4 | doi.org/10.18637/jss.v019.i09 |
| 16 | Gramacy, Taddy, "Categorical inputs, sensitivity analysis, optimization and importance tempering with tgp version 2", JSS 33(6), 2010 | full text: eq 1, sec 5, sec 5.3 | doi.org/10.18637/jss.v033.i06 |
| 17 | Hill, "Bayesian nonparametric modeling for causal inference", JCGS 20(1):217-240, 2011 | NOT FETCHED - paywalled; the IHDP construction was instead read from the maintainer's own implementation, `vdorie/npci`, `examples/ihdp_sim/data.R` and its README, in full | doi.org/10.1198/jcgs.2010.08162 ; github.com/vdorie/npci |
| 18 | Dorie, Hill, Shalit, Scott, Cervone, "Automated versus do-it-yourself methods for causal inference", Statistical Science 34(1), 2019 | full text: sec 4.2, 4.3, 6, 7.2, 7.3, 8 | arxiv.org/abs/1707.02641 |
| 19 | Hahn, Dorie, Murray, "Atlantic Causal Inference Conference (ACIC) data analysis challenge 2017", 2019 | full text: sec 2, 3, 4.1, 4.2, 6.2 | arxiv.org/abs/1905.09515 |
| 20 | ACIC 2019 data challenge, official site | fetched: estimand, tracks, dataset counts, covariate sources, stated difficulties | sites.google.com/view/acic2019datachallenge |
| 21 | Hahn, Murray, Carvalho, "Bayesian regression tree models for causal inference", Bayesian Analysis 15(3):965-1056, 2020 | full text (arXiv v4): sec 4 and Table 1, sec 6.1 and Tables 2-3, sec 6.2 | arxiv.org/abs/1706.09523 |
| 22 | Neal, "Slice sampling", Annals of Statistics 31(3):705-767, 2003 | full text: sec 8 | arxiv.org/abs/physics/0009028 |
| 23 | posteriordb (Stan development team) | repository tree and files read directly: the posterior list, `eight_schools_centered.stan`, `eight_schools.json` | github.com/stan-dev/posteriordb |
| 24 | Girolami, Calderhead, "Riemann manifold Langevin and Hamiltonian Monte Carlo methods", JRSS-B 73(2):123-214, 2011 | NOT FETCHED - paywalled; the dataset collection and the min-ESS-per-second metric were verified through Kleppe's replication instead | doi.org/10.1111/j.1467-9868.2010.00765.x |
| 25 | Kleppe, "Adaptive step size selection for Hessian-based manifold Langevin samplers", 2015 - verifying secondary for row 24 | full text: sec 5 and Table 2 | arxiv.org/abs/1501.07454 |
| 26 | Piironen, Vehtari, "Sparsity information and regularization in the horseshoe and other shrinkage priors", EJS 11(2), 2017 | full text: sec 2.4 | arxiv.org/abs/1707.01694 |
| 27 | Jasra, Holmes, Stephens, "MCMC methods and the label switching problem in Bayesian mixture modeling", Statistical Science 20(1):50-67, 2005 | NOT FETCHED - record and abstract only | doi.org/10.1214/088342305000000016 |
| 28 | Haario, Saksman, Tamminen, "An adaptive Metropolis algorithm", Bernoulli 7(2):223-242, 2001 - the twisted-Gaussian "banana" | NOT FETCHED - paywalled | doi.org/10.2307/3318737 |
| 29 | Kim, Loh, Shih, Chaudhuri, the 52-dataset collection behind both the CGM bake-off and SoftBart's ten | NOT FETCHED - named only, through rows 1 and 5 | (see rows 1 and 5) |

---

## 8. What could not be fetched

Stated plainly, because a survey that quietly recalls a paper it could not
read is worse than one that says so.

- **Chipman, George and McCulloch 1998** (Bayesian CART). The reachable
  copy is a scanned JSTOR PDF whose only text layer is the cover page; no
  OCR was available. Nothing in this document rests on it. The 2010
  paper's characterization of its predecessor's sampler ("tends to quickly
  gravitate toward a single large tree and then gets stuck in a local
  neighborhood of that tree") is quoted from the 2010 paper, which WAS
  fetched, and the restart advice is quoted at second hand from Gramacy's
  JSS 19(9), which was also fetched.
- **Linero 2018**, the DART paper itself. No arXiv preprint was located
  and the JASA copy is paywalled. DART's test functions are described here
  only as Linero and Yang 2018 describes them.
- **Hill, Linero and Murray 2020**, the review. The SSRN delivery URL
  returns 403 and no open copy was found. No claim in this document is
  attributed to it, so the item the task asked for - "Hill, Linero, Murray
  2020's failure modes" - is NOT answered here.
- **Hill 2011**, the IHDP paper. Paywalled. The construction in section
  3.1 comes from the maintainer's own `vdorie/npci` implementation, read in
  full, which is a stronger source for what the simulation actually does
  but does not carry the paper's own labelling of response surfaces A, B
  and C; the letters used in section 3.1 are the code's.
- **Girolami and Calderhead 2011**. Paywalled; verified through Kleppe's
  replication, which recoded their method and reports their dataset
  collection and metric.
- **Jasra, Holmes and Stephens 2005**. Record and abstract only, through
  the Project Euclid landing page; the full text is behind access control.
  Section 4.4's methodological point (label-invariant functionals) is
  supported in-house rather than by that paper.
- **Haario, Saksman and Tamminen 2001**, the banana. Paywalled. No
  numeric claim about it appears above.
- **Kim, Loh, Shih and Chaudhuri 2007**, the 52-dataset collection.
  Named through rows 1 and 5 of the reference table and not fetched; the
  dataset-level details quoted are CGM's and SoftBart's descriptions of it.
- **The LGBB rocket-booster data**, named in both tgp papers, was not
  located.
- **ACIC 2018.** No 2018 data challenge documentation was found; the
  documented challenges are 2016, 2017 and 2019, and section 3 covers those
  three. If a 2018 challenge exists this survey did not find it.

---

## 9. Provenance

```
repo          /Users/vdorie/Repositories/dbarts, branch bartcore
written at    6e2fcc47
in-repo       docs/design/tree-mixing-proposals.md sections 1-3, 6.3, 13,
reading       14; docs/design/grow-from-root-default.md and
              docs/design/forest-ranef-interweaving.md through the mixing
              survey's own verified record rather than re-derived;
              benchmarks/README.md, benchmarks/R/ and inst/common/ for what
              already exists
scope         survey only - no source touched, no default moved, nothing
              scheduled, no cell built
method        each external claim fetched and read from the primary source
              in this arc, or explicitly marked NOT FETCHED in section 8.
              Where a summarizing fetch returned a paraphrase, the PDF was
              re-extracted locally and the passage read directly; that is
              how every "verified" tag above was obtained
cost model    "unit fit" = the n 2000, p 10, m 200, 3000-sweep single-
              threaded fit that tree-mixing-proposals section 13 measured
              at 4.2 to 6.7 s on arm64 macOS. Unit-fit counts in section 6
              are estimates scaled from that anchor by n and sweeps, not
              measurements
```

---

## 10. Pilot results (2026-09-06)

**Every record in this section was measured at the mixture in force on
2026-09-06, `birth_death 0.5, swap 0.1, change 0.4, birth 0.5`. The default
has since dropped swap to zero, `birth_death 0.6, swap 0, change 0.4, birth
0.5`
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)),
which is this section's own `noswap` arm.** Nothing below was re-measured:
"default", "shipped default" and "shipped mixture" in this section name the
swap-carrying mixture these numbers were taken at, not the mixture that
ships. The swap-carrying arm is still reachable and still runs, as the
`swap` arm of both P2 scripts.

Section 6.5 named P2, P6 and P5 as the three to run first and C1 as the core
cell to gate on. This is that run, against the then-shipped default sampler.
Sections 1 to 9 are the survey as written and are unchanged. The cells live
in `benchmarks/R/surfaces`, one script each plus a shared
`surfaces-common.R`; every generating process there is transcribed from the
primary source and carries its citation.

Two corrections to the survey's transcriptions came out of re-reading the
sources before the run, and both change what gets generated.

- **P2's covariates.** Section 1.2 gives x1 and x3 in two blocks each and is
  silent on x2. The paper draws x2 in THREE blocks:
  `x2 ~ unif(0.1,0.4)` for i <= 100, `unif(0.6,0.9)` for i = 101..200 and
  `unif(0.1,0.9)` for i = 201..300 [verified: arXiv 1312.1895 sec 2.3]. A
  uniform x2 is a different problem: the mean function's second level is
  defined by a cut on x2 at 0.5, and the block structure is what aligns it
  with the x1/x3 confounding. The cell uses the paper's version.
- **C1's Single index gamma.** Section 1.4 gives one gamma, Linear's. Table
  1 gives Single index its own, `gamma_j = -1.5 + (j-1)/3`, against Linear's
  `gamma_j = -2 + 4(j-1)/(d-1)` [verified: arXiv 2002.03375v4 Table 1]. The
  cell uses the Table 1 pair. Table 4's noise levels are also
  `kappa` in {1, 2}, not section 4.1's {1, 10}.

### 10.1 P2, the confounded step function

Five matched data seeds, 8 chains, m = 1, 1000 burn-in and 2000 kept at
`n.thin = 1`, two arms through `proposal.probs`. Mean over seeds, min-max in
parentheses. The acceptance proxy is the share of sweeps on which the tree
structure differs from the previous sweep's, which at one tree and no
thinning is exactly the share on which a structural move was accepted.

    design      arm          between-chain sd    pooled p(root x1)   share on x1 given {x1,x3}
    confounded  default      0.501(0.470-0.518)  0.431(0.361-0.625)  0.513(0.375-0.625)
    confounded  birth/death  0.455(0.354-0.535)  0.350(0.125-0.625)  0.504(0.167-0.833)
    duplicate   default      0.051(0.017-0.075)  0.499(0.481-0.512)  -
    duplicate   birth/death  0.492(0.354-0.535)  0.425(0.125-0.625)  -

    design      arm          acceptance proxy       root switches/chain  min switches
    confounded  default      0.0584(0.0520-0.0632)  0.7(0.0-1.2)         0
    confounded  birth/death  0.0298(0.0222-0.0355)  0.0(0.0-0.0)         0
    duplicate   default      0.1030(0.0907-0.1132)  78.5(70.2-85.0)      38
    duplicate   birth/death  0.0097(0.0061-0.0140)  0.0(0.0-0.0)         0

| statistic | published | measured, shipped default | verdict |
|---|---|---|---|
| acceptance rate of tree moves after burn-in | 0 | 0.058 (0.052-0.063) of sweeps change the tree | does not, literally |
| root variable moves per chain per 2000 sweeps | implied 0 | 0.7 (0.0-1.25); at least one chain never moves in all 5 seeds, no chain moves at all in 1 of them | reproduces |
| between-chain sd of p(root on x1), primary | mixing null near 0 | 0.501 (0.470-0.518) | reproduces |
| pooled share on x1 among {x1, x3} draws | 0.5 exactly | 0.513 (0.375-0.625) | reproduces, pooled only |

**Verdict: reproduces, with the published number corrected.** The acceptance
rate is not zero - the tree changes on about one sweep in seventeen - but
what it never changes is the root variable. A chain visits 1.3 root
variables on average over two thousand sweeps, and the between-chain sd of
0.501 is 94 percent of the 0.535 a 0/1 statistic can reach at eight chains.
The exact 0.5 symmetry is satisfied only after pooling the eight. That is
representation multimodality read directly, and it is the failure Pratola's
sentence describes even though the number in it does not survive.

**Arm contrast.** Birth/death only halves the structural acceptance rate
(0.030 against 0.058) and moves the root variable zero times in 40 of 40
chains, against 0.7 per chain for the shipped mixture. Neither arm mixes
between representations on the confounded design at this length, so the
pathology itself does not separate them. The null control does. With two
exactly duplicated columns - likelihood and prior ratios both exactly 1 -
the shipped mixture returns a pooled 0.499 with 78.5 switches per chain and
never fewer than 38, while birth/death only returns 0 switches in every
chain and a between-chain sd of 0.49. **Birth/death only fails the null
control this cell requires; the shipped mixture passes it.** The change move
is what supplies variable switching at one tree, and death-then-rebirth does
not substitute for it even when the two representations are exactly
exchangeable. Note that the null control's columns are drawn on a
four-value grid: with continuous duplicates a rule-changing proposal has to
hit the twin column at the same cut out of the whole grid, which makes the
switch rare for a reason unrelated to the kernel.

**No-swap arm - now the shipped mixture.** Arm C of
[14.2 Design](tree-mixing-proposals.md#142-design) (birth_death 0.6,
swap 0, change 0.4, birth 0.5) was added to this cell and run on the same
five seeds; it is the mixture that ships. Null control, pooled p(root x1),
switches per chain (min-max), minimum switches, between-chain sd, chains
parked off the pair (of 40):
default 0.499 (0.481-0.512), 78.5 (70.2-85.0), 38, 0.051 (0.017-0.075), 0;
birth/death only 0.425 (0.125-0.625), 0.0 (0.0-0.0), 0, 0.492 (0.354-0.535),
12; no-swap 0.457 (0.359-0.555), 70.8 (57.9-81.2), 0, 0.149 (0.036-0.251), 5.
No-swap's mean switches per chain falls inside the default's seed-to-seed
range, but 5 of its 40 chains never switch at all: every one sits at an x3
root with 3 to 4 interior nodes and a child already split on x1, a
representation no default chain ever visits. Change at an x3 root is vetoed
once a child splits on x1, and death of that child loses the signal rather
than returning it to the root ([`changeMove`](../../src/bartcore/moves.hpp)); swap
is the only move that rotates a child's rule up to the root
([`swapMove`](../../src/bartcore/moves.hpp)). On this null control change alone
does not carry representation switching; swap does the rule rotation. This is
the one measured cost of the default's zero, and it is confined to m = 1: the
paragraph below finds no stuck tree at 50 or 200, and a single-tree caller can
set swap positive.

**At production tree counts.** `P2-null-at-scale.R` re-runs the
duplicate-column null at m = 50 and m = 200 trees, swap-carrying and no-swap
arms only, 5 seeds x 8 chains x 2000 kept sweeps, reading a tree-level
root share instead of a chain-level one. Mean x3 share (min-max over
chains), between-chain sd of the x3 share, total stuck-on-x3 trees over
all chains: default 50 0.281(0.266-0.299), 0.008(0.005-0.010), 0; default
200 0.323(0.313-0.331), 0.003(0.002-0.005), 0; no-swap 50
0.283(0.268-0.305), 0.009(0.007-0.012), 0; no-swap 200 0.323(0.317-0.330),
0.003(0.002-0.004), 0. The duplicated pair's own mean share sits at
0.307-0.330 across the same four cells. The single-tree stuck-on-x3
representation does not survive at either tree count: zero trees, in
every arm and both tree counts, hold an x3 root for every kept draw of
their chain.

### 10.2 P6, the diagonal shelf with targeted selection

200 replications, n = 250, true effect -1, shipped defaults, ATE read off
the counterfactual contrast. The BCF arm is dbarts' own: an estimated
propensity added as a covariate and a `forest()` term the treatment
indicator modulates.

    reconstruction  arm   bias   coverage  rmse   interval length
    figure          bart  0.314  0.590     0.366  0.695
    figure          bcf   0.100  0.885     0.214  0.729
    shelf 0.15      bart  0.592  0.260     0.643  0.821
    shelf 0.40      bart  0.172  0.840     0.235  0.672

| statistic | published | measured, figure reconstruction | verdict |
|---|---|---|---|
| BART ATE bias, primary | 0.27 | 0.314 | partial |
| BART 95% coverage, primary | 0.65 | 0.590 | partial |
| BART RMSE | 0.31 | 0.366 | partial |
| BCF ATE bias | 0.14 | 0.100 | partial |
| BCF 95% coverage | 0.95 | 0.885 | partial |
| BCF RMSE | 0.21 | 0.214 | partial |

**Verdict: partial, and the reason is the source, not the sampler.** Every
measured value sits within 0.07 of its published counterpart, and the
BART-to-BCF gap - bias 0.31 to 0.10, coverage 0.59 to 0.89 - is the gap the
paper reports. But the paper never prints its prognostic function. Section
3.5 recorded that; what this run adds is what it costs. Across three
reconstructions that all satisfy the paper's stated constraints, BART's bias
runs from 0.17 to 0.59 and its coverage from 0.84 to 0.26. The
figure-calibrated one agrees with the published row because it was built to
agree with the figure the published row came from, so the agreement
calibrates the reconstruction and cannot on its own confirm the sampler.

Two source findings worth carrying. First, the printed propensity is
incomplete. The equation reads
`0.8 Phi(mu / (0.1(2 - x1 - x2) + 0.25)) + 0.025(x1 + x2) + 0.05`, but the
paper's own LaTeX source carries the generating expression one line above it
as a comment, `0.8*pnorm(m/3, 0, 0.1*(2-xtilde)+0.25) + 0.025*xtilde + 0.05`,
which divides mu by 3 first. Only the commented form reproduces the paper's
Figure 4: it puts the propensity near 0.21 at mu = -1 and near 0.10 at
mu = -1.9, where the printed form puts it near 0.08 and 0.05. Second,
Figure 4 pins mu harder than the captions do - a realized range of about
-1.9 to 3.0, not the symmetric -3 to 3 the Figure 3 caption states, and a
median well below zero. That is why no member of the symmetric near-step
family the captions describe can be made to fit, and why the `figure`
reconstruction is asymmetric.

**Move-set arms (2026-09-07).** The table above was measured at the
mixture in force on 2026-09-06 (`birth_death 0.5, swap 0.1, change 0.4`),
the former default section 10's preamble names, not the mixture that
ships: its three `bart` rows reproduce digit for digit under this build's
`swap` arm and not under `default` (`birth_death 0.6, swap 0, change
0.4`), so "measured, figure reconstruction" above names the swap-carrying
mixture at every `bart` row, not the shipped one. `P6-diagonal-shelf.R
default birthdeath swap` reruns all three, paired, on the same 200
replications per reconstruction:

    reconstruction  moveset     bias   coverage  rmse
    figure          default     0.314  0.600     0.366
    figure          birthdeath  0.315  0.580     0.366
    figure          swap        0.314  0.590     0.366
    shelf 0.15      default     0.594  0.225     0.646
    shelf 0.15      birthdeath  0.593  0.250     0.644
    shelf 0.15      swap        0.592  0.260     0.643
    shelf 0.40      default     0.173  0.850     0.236
    shelf 0.40      birthdeath  0.173  0.840     0.236
    shelf 0.40      swap        0.172  0.840     0.235

Paired against `default`
([`surfacesPairedDifference`](../../benchmarks/R/surfaces/surfaces-common.R),
[`surfacesMarginVerdict`](../../benchmarks/R/surfaces/surfaces-common.R)),
every `|bias|` and RMSE ratio sits at 0.993 to 1.003 against the 1.02
margin - within margin in all six cells. Coverage's point estimate crosses
the -0.010 margin in three cells (figure birthdeath -0.020 t -1.42, figure
swap -0.010 t -0.71, shelf 0.40 swap -0.010 t -0.82) but none separates
from zero at the required one-sided 95% bound, so each reads "past margin,
not separated" rather than FLAG; the largest |t| anywhere is 2.36 (shelf
0.15 swap coverage, +0.035), short of the four-times-SE improvement bar
([`surfacesImprovementVerdict`](../../benchmarks/R/surfaces/surfaces-common.R))
every cell reads below. Nothing is flagged and nothing clears the bar: the
move set does not move this cell.

`bcf` cannot take an arm - a treatment forest refuses a non-default
`proposal.probs` outright - so it always runs at whatever
`proposal.probs` resolves to be the default, and its row moves with the
default's own change: the pilot's 0.100/0.885/0.214/0.729
(bias/coverage/rmse/interval length, also the swap-carrying mixture) reads
0.101/0.890/0.213/0.726 on this build. The BART-to-BCF gap the verdict
above reads from is the same gap either way.

### 10.3 P5, the checkerboard on an autocorrelated design

Twenty matched seeds, n = 1600, p = 40, 8 chains, 1000 burn-in and 2000
kept, shipped defaults. There is no published dbarts number; the oracle is
the inclusion truth and its near-decoys.

    statistic                                  measured
    between-chain sd of inclusion, true cols   0.0109(0.0080-0.0154)
    its mixing null, sd/sqrt(ESS)              0.0064(0.0055-0.0074)
    ratio to the null, true cols               1.72(1.19-2.52)
    ratio to the null, immediate decoys        1.61(1.29-1.93)
    summed inclusion share, 4 true cols        0.510(0.490-0.534)
    summed inclusion share, 8 decoys           0.176(0.152-0.199)
    largest non-true column's inclusion        0.030(0.024-0.040)
    non-true columns above the weakest true    0.0(0.0-0.0)
    95% pointwise coverage of true f           0.984(0.975-0.991)
    held-out RMSE                              0.704(0.609-0.909)

| statistic | reference | measured | verdict |
|---|---|---|---|
| between-chain sd of inclusion, primary | 1.0x its own mixing null if chains agree | 1.72x (1.19-2.52) | partial |
| inclusion oracle, {x5, x10, x15, x20} | exact | all four top the ranking in 20 of 20 seeds; no other column outranks the weakest true one | reproduces the truth |

**Verdict: the pathology does not appear at the shipped default.** Tan et
al.'s Theorem 5.2 bounds the hitting time for a pure interaction when the
change move is disallowed; the shipped mixture is not that, and this cell is
where the gap gets tested. It comes out clean. Each true column carries
about 0.128 of the splits against 0.022 for its correlated-0.9 neighbours
and 0.010 for the far columns; the four true columns take 0.51 of all splits
between them; and in all twenty replicates none of the other thirty-six
columns outranks the weakest true one. Coverage is 0.98 against a nominal
0.95,
so the intervals are conservative rather than short. The one thing not
clean is the between-chain spread: chains disagree about inclusion by 1.7
times their own Monte Carlo resolution, and by about the same factor on the
decoys as on the true columns, so the disagreement is a shared-splits effect
across the correlated block rather than a wrong answer.

**Move-set arms (2026-09-07).** The table above was measured at the same
mixture 10.2 was, the former default section 10's preamble names
(`birth_death 0.5, swap 0.1, change 0.4`), not the mixture that ships: all
ten of its statistics reproduce digit for digit under this build's `swap`
arm and not under `default` (`birth_death 0.6, swap 0, change 0.4`), so
every row above is the swap-carrying mixture, not the shipped one.
`P5-checkerboard.R default birthdeath swap` reruns all three, paired, on
the same twenty seeds:

    arm         between-chain sd  ratio to null  inclusion share  coverage  RMSE
    default     0.0117            1.84           0.501            0.984     0.710
    birthdeath  0.0123            1.82           0.478            0.984     0.737
    swap        0.0109            1.72           0.510            0.984     0.704

(between-chain sd and its ratio to the mixing null are of inclusion on the
four true columns; inclusion share is their summed share of all splits.)

Paired against `default`
([`surfacesPairedDifference`](../../benchmarks/R/surfaces/surfaces-common.R),
[`surfacesMarginVerdict`](../../benchmarks/R/surfaces/surfaces-common.R)),
`birthdeath` FLAGS on two of the three gated metrics: inclusion share
-0.023 +/- 0.010 (0/20 seeds positive) t -10.14, past the -0.010 margin
and separated; and held-out RMSE +0.028 t 4.44 (ratio 1.039), past the
1.02 margin and separated. Its coverage is unchanged (+0.000). `swap`
stays within margin on all three: inclusion +0.009 (15/20) t 3.16,
coverage +0.000, RMSE -0.005 (ratio 0.993). On the primary - the
between-chain sd itself - neither arm clears the four-times-SE improvement
bar
([`surfacesImprovementVerdict`](../../benchmarks/R/surfaces/surfaces-common.R)):
`birthdeath` moves it +0.0006 (12/20) t 1.40, `swap` moves it -0.0008
(7/20) t -1.21, both below 4x SE. The inclusion oracle holds in every arm -
each true column carries 0.117 to 0.128 of the splits, and no decoy
outranks the weakest true column in 20 of 20 seeds under any of the three.

`birthdeath` is `pi_c = 0`, the premise
[6.3 The pathologies](#63-the-pathologies) names for Tan et al.'s Theorem
5.2, which bounds a pure interaction's hitting time when the change move
is disallowed. This cell is where that premise's cost turns into a number:
cutting the change move to zero does not break the inclusion oracle, but
it costs 0.023 of the true columns' summed share and 3.9 percent of
held-out RMSE, both separated from noise. On this correlated design the
change move is what buys back the inclusion share and accuracy the
theorem's premise gives up.

Per [6.1 The rule, stated operationally](#61-the-rule-stated-operationally),
a flag counts only after a mandatory fresh-seed re-run of the flagged
cell. That re-run has been done: `P5-checkerboard.R` gained a second seed
block (`block2`, seeds 21 to 40, the default block's outputs unchanged) and
`default` against `birthdeath` was run again on it, 20 matched pairs, 8
chains.

    statistic                    default                birthdeath
    between-chain sd, true cols  0.0115(0.0088-0.0139)  0.0124(0.0081-0.0164)
    ratio to the mixing null     1.91(1.30-2.40)        1.81(1.17-2.75)
    between-chain sd, decoys     0.0060(0.0044-0.0072)  0.0076(0.0064-0.0099)
    inclusion share, true cols   0.503(0.481-0.549)     0.476(0.453-0.504)
    95% coverage of true f       0.984(0.971-0.992)     0.985(0.974-0.994)
    held-out RMSE                0.705(0.622-0.853)     0.733(0.660-0.888)

    paired, birthdeath - default   mean     sd      positive  t       one-sided 95% bound
    primary, between-chain sd     +0.0009  0.0022  15/20     +1.87   no margin; not an improvement
    decoy between-chain sd        +0.0016  0.0012  17/20     +6.01   no margin
    inclusion share, true cols    -0.0267  0.0117   0/20    -10.21   0.022 worse vs margin 0.010: FLAG
    95% coverage                  +0.0003  0.0039  10/20     +0.29   within margin
    held-out RMSE                 +0.0276  0.0228  19/20     +5.43   ratio 1.039 vs 1.02: FLAG

Both flags REPRODUCE on fresh seeds under 6.4's two conditions (the paired
mean worse than the margin and the one-sided bound excluding the null), so
they count: on this cell birth/death only loses 0.027 of the true columns'
inclusion share and 3.9 percent of held-out RMSE against the shipped
mixture. Coverage stays within margin; the primary again fails to separate
as an improvement and leans the wrong way (+0.0009, t +1.87, after +0.0006
in the first block), so `birthdeath` has twice failed the four-times-SE bar
it would need to win this pathology; the decoy spread is worse (t +6.01).
The inclusion oracle holds in 20 of 20 fresh seeds in both arms.

### 10.4 C1, the He and Hahn factorial

Twenty matched seeds per cell, n = 10000, p = 30, kappa = 1, one chain, 1000
burn-in and 2500 kept - the paper's own chain length. Coverage, length and
RMSE are of the true f at the 10000 training rows, which is the readout the
published table uses; the held-out column is the survey's own pre-registered
thousand-row version, and it agrees with the in-sample one to within 0.01
throughout. Section 5 of the paper, which carries Table 4, fixes the sample
size, the noise level and the chain length but says neither which of section
4.1's two predictor arms it used nor how many trees, so both are varied as
diagnostic arms beside the pre-registered one.

    mean fn      arm             95% coverage        interval length     RMSE                min ESS
    trigpoly     correlated 75   0.850(0.763-0.912)  3.76(3.45-4.04)     1.26(1.17-1.39)     2(1-5)
    trigpoly     correlated 200  0.957(0.916-0.979)  5.19(4.96-5.46)     1.26(1.13-1.36)     5(1-13)
    trigpoly     independent 75  0.822(0.786-0.862)  3.32(3.13-3.56)     1.25(1.20-1.30)     2(1-4)
    trigpoly     independent 200 0.922(0.902-0.943)  4.13(3.99-4.29)     1.20(1.15-1.28)     2(1-4)
    singleindex  correlated 75   0.893(0.852-0.926)  8.96(8.60-9.43)     2.60(2.46-2.79)     2(1-4)
    singleindex  correlated 200  0.965(0.945-0.982)  11.25(10.69-11.77)  2.54(2.36-2.80)     3(2-4)
    singleindex  independent 75  0.822(0.771-0.858)  5.78(5.47-6.22)     2.10(2.05-2.21)     2(1-3)
    singleindex  independent 200 0.924(0.894-0.935)  6.96(6.74-7.20)     1.93(1.84-2.09)     4(2-6)

| statistic | published | measured, pre-registered arm (correlated, 75) | measured, closest arm (independent, 75) | verdict |
|---|---|---|---|---|
| Trig+poly 95% coverage, primary | 0.74 | 0.850 | 0.822 | partial |
| Trig+poly interval length | 2.89 | 3.76 | 3.32 | partial |
| Trig+poly RMSE | 1.27 | 1.26 | 1.25 | reproduces |
| Single index 95% coverage, primary | 0.73 | 0.893 | 0.822 | partial |
| Single index interval length | 4.62 | 8.96 | 5.78 | does not / partial |
| Single index RMSE | 2.08 | 2.60 | 2.10 | does not / reproduces |

**Verdict: partial, and the pre-registered predictor arm is the wrong one.**
The evidence is Single index's RMSE. On the correlated-factor design dbarts
returns 2.60 against a published 2.08 and an interval nearly twice the
published length; on the independent design it returns 2.10 and 5.78.
Trig+poly's RMSE matches on both designs (1.25 to 1.26 against 1.27), so it
does not discriminate, but Single index does, and it says Table 4 was run on
the independent standard-normal predictor arm. Section 6.2's instruction to run
the correlated arm should be read as a choice this battery makes for
realism, not as a reproduction of the paper's setting.

At the setting closest to the paper's - independent design, shipped 75 trees
- the coverage deficit reproduces in direction and about half in size: 0.82
against a nominal 0.95 on both mean functions, where the paper reports 0.73
and 0.74. dbarts is better calibrated than the published BART at the same
point accuracy: its intervals are 15 percent longer on Trig+poly and 25
percent longer on Single index at the same RMSE.
Raising the tree count to 200 buys 10 points of coverage on both functions
and costs nothing in RMSE.

The secondary is worth stating plainly because it is the same number in
every shipped-default arm: **the minimum effective sample size over 25 fixed
points is 2, out of 2500 kept draws.** Not 2 percent - two draws. Coverage
near nominal at 200 trees and an ESS of 2 to 5 are both true of the same
chain.

**Grow-from-root arm.** `correlated75grow` is the pre-registered
`correlated75` arm's settings plus `n.grow.sweeps = 5` (k = 5, stochtree's
own default; neither `man/bart2.Rd` nor `grow-from-root.md` names a study
value of its own - `grow-from-root-default.md`), the same twenty seeds.
Trig+poly: 95% coverage 0.856(0.792-0.917), interval length
3.80(3.59-4.12), RMSE 1.26(1.17-1.38), min ESS 2(1-5). Single index: 95%
coverage 0.881(0.839-0.915), interval length 8.80(8.19-9.26), RMSE
2.63(2.45-2.82), min ESS 2(2-5). Paired `correlated75`: Trig+poly
0.850(0.763-0.912), 3.76(3.45-4.04), 1.26(1.17-1.39), 2(1-5); Single index
0.893(0.852-0.926), 8.96(8.60-9.43), 2.60(2.46-2.79), 2(1-4). Mean wall
time per fit is 9.4s against 8.3s for the no-grow arm. He and Hahn report
0.91-0.96 coverage for their own warm start (Table 4, WS-BART, k = 2). The
grow-from-root arm moves coverage, RMSE and minimum ESS by less than the
no-grow arm's own seed-to-seed range on both mean functions.

**Independent-design grow arm.** `independent75grow` is the `independent75`
arm's settings plus `n.grow.sweeps = 5`, the same k = 5 rationale as
`correlated75grow`, the same twenty seeds. Trig+poly: 95% coverage
0.797(0.747-0.838), interval length 3.20(3.00-3.37), RMSE 1.28(1.21-1.37),
min ESS 2(1-5). Single index: 95% coverage 0.803(0.758-0.830), interval
length 5.65(5.44-5.95), RMSE 2.13(2.05-2.28), min ESS 2(1-3). Paired
`independent75`: Trig+poly 0.822(0.786-0.862), 3.32(3.13-3.56),
1.25(1.20-1.30), 2(1-4); Single index 0.822(0.771-0.858), 5.78(5.47-6.22),
2.10(2.05-2.21), 2(1-3). Mean wall time per fit is 9.3s against 8.0s for
the no-grow arm. The grow-from-root arm moves coverage, RMSE and minimum
ESS by less than the no-grow arm's own seed-to-seed range on both mean
functions.

**Chain configuration (2026-09-07).** Every arm above ran ONE chain of 2500
kept draws, the paper's configuration; `bart2`'s shipped default is four
chains of 500 burn-in and 500 kept, pooled. Three arms were added to the
independent design at 75 trees on the same twenty seeds, at the kernel that
ships (swap at zero; the recorded arms above were measured with swap at
0.1, and the 2500-draw prefix of the long chain re-reads 0.818 and 0.817
against their 0.822, inside the seed range): `independent75pool4`, four
chains of 500 + 500 pooled, the shipped default exactly;
`independent75pool4long`, four chains of 1000 + 2500 pooled; and
`independent75long`, one chain of 1000 + 25000. Coverage, interval length
and RMSE are of the pooled draws; minimum ESS is per chain and summed;
"between" is the median over the 25 ESS points of the between-chain
standard deviation of each chain's posterior mean divided by the pooled
posterior standard deviation, near 0 when the chains agree and near 1 when
each sits in its own place.

    mean fn      arm                     95% coverage        length  RMSE   min ESS (sum)  per chain  between
    trigpoly     independent75pool4      0.961(0.945-0.977)  4.61    1.12   15(8-31)       2(1-2)     0.78
    trigpoly     independent75pool4long  0.959(0.937-0.975)  4.19    1.04   18(9-46)       2(1-2)     0.65
    trigpoly     independent75long       0.902(0.858-0.928)  3.72    1.16   2(1-3)         2(1-3)     -
    singleindex  independent75pool4      0.895(0.878-0.915)  6.45    1.92   21(9-32)       2(1-2)     0.68
    singleindex  independent75pool4long  0.905(0.886-0.930)  6.62    1.94   18(10-27)      2(1-2)     0.49
    singleindex  independent75long       0.888(0.860-0.914)  6.45    2.00   3(2-10)        3(2-10)    -

Within the long chain, read on the same fit: Trig+poly 0.818 (0.770-0.866)
at 2500 kept draws and 0.902 (0.858-0.928) at 25000; Single index 0.817
(0.759-0.871) and 0.888 (0.860-0.914).

Two readings, both against the single-chain deficit of 0.82. First, the
deficit is chain exploration, and the shipped default already carries most
of the remedy: four pooled chains at the default length lift coverage to
0.96 on Trig+poly and 0.90 on Single index at equal or better RMSE, the
between-chain ratio of 0.5 to 0.8 says each chain sits in its own place so
that pooling is what widens the interval, and the 500 + 500 configuration
matches the 1000 + 2500 one at a fifth of the compute, so chain count and
not chain length is the active ingredient. Second, the posterior is right
and the sampler is slow: within one fit coverage climbs from 0.82 to 0.90
as the chain runs ten times longer while its minimum ESS stays at 2 to 3,
and that tenfold chain buys less than four short chains buy at a fifth of
its cost. The tree-count contrast above is the same mechanism seen from
the prior's side and is not a default question; the measure a kernel
change has to move on this cell is the per-chain ESS. The host carried a
1-minute load of 8 to 15 throughout, so the wall times in the script's
output carry no timing claim.

**Move-set arms in the four-chain configuration (2026-09-07).** Four
further arms hold the shipped four-chain configuration fixed and vary only
`proposal.probs`, re-reading
[13. Move-set A/B (2026-09-06)](tree-mixing-proposals.md#13-move-set-ab-2026-09-06)'s
contrast at the chain count that ships rather than at that grid's one
chain: `independent75pool4bd` (birth_death 1, birth/death only),
`independent75pool4swap` (birth_death 0.5, swap 0.1, change 0.4, the
former default), `independent75pool4perturbMixed` (birth_death 0.5, change
0.34, perturb 0.16, the 0.16 taken from change and birth/death both) and
`independent75pool4perturbB`
([5.2 Arms, dosage, and cells](perturb-move.md#52-arms-dosage-and-cells)'s
own arm B, birth_death 0.6, change 0.24, perturb 0.16 at `w = 1, d =
0.16`, the 0.16 taken from change alone). `independent75pool4` is re-run
in the same session as their paired control and reproduces the Chain
configuration paragraph's own row digit for digit. Twenty matched seeds,
both mean functions; columns as above.

    mean fn      arm                             95% coverage        length  RMSE   min ESS (sum)  per chain  between
    trigpoly     independent75pool4              0.961(0.945-0.977)  4.61    1.12   15(8-31)        2(1-2)     0.78
    trigpoly     independent75pool4bd            0.963(0.950-0.978)  4.98    1.20   15(7-30)        1(1-2)     0.81
    trigpoly     independent75pool4swap          0.958(0.935-0.974)  4.70    1.18   14(7-23)        1(1-2)     0.81
    trigpoly     independent75pool4perturbMixed  0.959(0.936-0.974)  4.63    1.14   18(6-32)        2(1-2)     0.67
    trigpoly     independent75pool4perturbB      0.962(0.948-0.977)  4.66    1.14   15(9-24)        2(1-2)     0.67
    singleindex  independent75pool4              0.895(0.878-0.915)  6.45    1.92   21(9-32)        2(1-2)     0.68
    singleindex  independent75pool4bd            0.905(0.892-0.924)  6.56    1.88   19(10-28)       2(1-2)     0.74
    singleindex  independent75pool4swap          0.894(0.864-0.913)  6.45    1.92   16(9-35)        2(1-2)     0.67
    singleindex  independent75pool4perturbMixed  0.888(0.859-0.913)  6.42    1.98   20(10-32)       2(2-3)     0.49
    singleindex  independent75pool4perturbB      0.889(0.856-0.912)  6.44    1.98   23(8-36)        2(1-3)     0.50

Paired differences against `independent75pool4`, mean +/- sd (seeds
positive of 20); a t is given only where it separates.

    mean fn      arm                             d min ESS (sum)              d 95% coverage                 d RMSE
    trigpoly     independent75pool4bd            -0.0 +/- 8.4 (10/20)         +0.002 +/- 0.008 (12/20)       +0.083 +/- 0.099 (16/20)
    trigpoly     independent75pool4swap          -0.6 +/- 9.2 (11/20)         -0.003 +/- 0.009 (6/20)        +0.062 +/- 0.133 (15/20)
    trigpoly     independent75pool4perturbMixed  +3.1 +/- 8.7 (12/20) t 1.57  -0.002 +/- 0.008 (9/20)        +0.023 +/- 0.077 (12/20)
    trigpoly     independent75pool4perturbB      +0.1 +/- 8.1 (12/20) t 0.06  +0.001 +/- 0.007 (10/20) t 0.43 +0.023 +/- 0.069 (14/20) t 1.45, ratio 1.020
    singleindex  independent75pool4bd            -2.1 +/- 8.6 (10/20)         +0.010 +/- 0.009 (18/20) t 5.2 -0.041 +/- 0.036 (3/20) t -5.0, ratio 0.979
    singleindex  independent75pool4swap          -4.3 +/- 8.4 (7/20)          -0.000 +/- 0.010 (8/20)        +0.004 +/- 0.025 (13/20)
    singleindex  independent75pool4perturbMixed  -0.5 +/- 9.6 (12/20)         -0.007 +/- 0.010 (4/20) t -2.9 +0.057 +/- 0.035 (20/20) t 7.2, ratio 1.030
    singleindex  independent75pool4perturbB      +1.9 +/- 11.3 (10/20) t 0.74 -0.006 +/- 0.010 (7/20) t -2.71 +0.058 +/- 0.034 (19/20) t 7.63, ratio 1.030

[13.5 What separated and what did not](tree-mixing-proposals.md#135-what-separated-and-what-did-not)'s
null - no move set separates from another - reproduces inside this
battery on its own primary cell and statistic: neither birth/death-only
nor the former swap default moves summed minimum ESS on Trig+poly (t
about 0 and -0.3). The two secondaries that do clear
[6.4 What "no regression on the core" means numerically](#64-what-no-regression-on-the-core-means-numerically)'s
margin both sit on the ungated Single index: birth/death-only is better, coverage
+0.010 (t 5.2) and RMSE ratio 0.979 (t -5.0); both perturb arms are
worse, RMSE ratio 1.030 on Single index (t 7.2 Mixed, t 7.63 B), with
Trig+poly close behind at a ratio of 1.020 for B (t 1.45). The two perturb
arms land at the same Single index ratio whichever move the 0.16 share is
drawn from - change alone for B, change and birth/death both for Mixed -
so the penalty reads as the size of the discarded share, not its source.
Both perturb arms also pull the between-chain ratio down without buying
ESS: 0.67 and 0.50 against the control's 0.78 and 0.68, while their own
summed minimum ESS differences from control (+3.1 and +0.1 on Trig+poly,
-0.5 and +1.9 on Single index) stay inside noise.

The paired standard error of summed minimum ESS runs 1.8 to 2.5 across
the eight cells above, at or below
[5.1 The chain configuration, and what it makes the primary statistic](perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic)'s
derived 2.0 bound on seven of the eight; the eighth, Single index against
`independent75pool4perturbB`, reads 2.5, above the bound, off two large
offsetting seeds rather than a shift in location. This is the calibration
the +8 bar was set against.

**Sham and fresh-seed re-run (2026-09-07).** perturb-move.md's slice 3 run adds a sham arm, `independent75pool4sham` - the
control against itself at sampler seeds offset 1000, same twenty data seeds - and a fresh-seed block, seeds 21 to 40, both
run through [`C1-he-hahn.R`](../../benchmarks/R/surfaces/C1-he-hahn.R)'s new `samplerOffset` and `seedBlock` options. Paired
differences, mean +/- sd (seeds positive of 20), against each row's own control on its own seeds:

    mean fn      arm                             seeds    d min ESS (sum)             d per-chain min ESS   d 95% coverage                   d RMSE, ratio
    trigpoly     independent75pool4sham          1-20     -2.3 +/- 9.7 (8/20) t -1.07  -0.04 +/- 0.26        +0.001 +/- 0.009 (13/20)         -0.010 +/- 0.060, ratio 0.991
    trigpoly     independent75pool4perturbB      21-40    -1.9 +/- 6.5 (9/20) t -1.27  +0.06 +/- 0.28 t 0.92 -0.004 +/- 0.009 (5/20) t -2.00  +0.046 +/- 0.068 t 3.02, ratio 1.042
    singleindex  independent75pool4perturbB      21-40    +9.5 +/- 8.0 (19/20) t 5.35  +0.37 +/- 0.26 t 6.47 -0.005 +/- 0.011 (6/20)          +0.051 +/- 0.034 t 6.72, ratio 1.027

Held-out RMSE ratio is 0.991 on the sham, 1.032 on the fresh-block Trig+poly contrast (t 2.36) and 1.028 on Single index
(t 6.27); wall ratio is 0.997, 1.041 and 0.996. The sham's summed minimum ESS, 14.82 against 12.51, sits inside the +8 bar
and is the harness's own calibration of it. **The fresh block reproduces the pilot's null on Trig+poly**: the pilot's own
+0.1 +/- 8.1 (t 0.06) against the same control, above, sharpens to -1.9 +/- 6.5 (t -1.27) at fresh seeds. It does not
reproduce on Single index, where the pilot read +1.9 +/- 11.3 (t 0.74) and the fresh block reads +9.5 +/- 8.0 (t 5.35),
separated from noise. Full verdict at
[5.3 What arm B must produce, and the kill](perturb-move.md#53-what-arm-b-must-produce-and-the-kill).

**rule_gibbs arms (2026-09-07).** nog-gibbs.md's slice 3 run adds two
arms on the same shipped four-chain configuration and the same twenty
seeds, varying only `proposal.probs`: `independent75pool4ruleGibbsB`
(birth_death 0.6, swap 0, change 0.24, rule_gibbs 0.16, the 0.16 taken
from change alone, which is
[6. Benefit, pre-registered](nog-gibbs.md#6-benefit-pre-registered)'s own
arm B at `d` = 0.16) and `independent75pool4ruleGibbs32` (change 0.08,
rule_gibbs 0.32, the second dosage of that grid). `independent75pool4` is
re-run in the same session as their paired control and reproduces the
Chain configuration paragraph's own rows digit for digit on both mean
functions. Both cells flagged, so arm B took the mandatory fresh-seed
re-run on seeds 21 to 40, control and all. Columns as above.

    mean fn      arm                             seeds  95% coverage        length  RMSE  min ESS (sum)  per chain  between
    trigpoly     independent75pool4              1-20   0.961(0.945-0.977)  4.61    1.12  15(8-31)       2(1-2)     0.78
    trigpoly     independent75pool4ruleGibbsB    1-20   0.939(0.911-0.958)  3.98    1.09  36(19-53)      3(2-4)     0.58
    trigpoly     independent75pool4ruleGibbs32   1-20   0.939(0.917-0.957)  3.96    1.09  52(23-97)      3(2-6)     0.55
    trigpoly     independent75pool4              21-40  0.964(0.950-0.977)  4.58    1.10  16(9-28)       2(1-2)     0.80
    trigpoly     independent75pool4ruleGibbsB    21-40  0.938(0.917-0.951)  3.96    1.08  39(15-55)      2(2-3)     0.59
    singleindex  independent75pool4              1-20   0.895(0.878-0.915)  6.45    1.92  21(9-32)       2(1-2)     0.68
    singleindex  independent75pool4ruleGibbsB    1-20   0.900(0.872-0.921)  6.57    1.97  35(16-57)      2(2-3)     0.40
    singleindex  independent75pool4ruleGibbs32   1-20   0.904(0.872-0.923)  6.61    1.96  36(18-74)      3(2-5)     0.37
    singleindex  independent75pool4              21-40  0.891(0.868-0.911)  6.43    1.93  15(9-27)       2(1-2)     0.67
    singleindex  independent75pool4ruleGibbsB    21-40  0.898(0.879-0.914)  6.55    1.97  39(20-60)      3(2-6)     0.39

Paired differences against each row's own control on its own seeds, mean
+/- sd (seeds positive of 20):

    mean fn      arm                             seeds  d min ESS (sum)               d per-chain min ESS    d 95% coverage                     d RMSE, ratio
    trigpoly     independent75pool4ruleGibbsB    1-20   +21.5 +/- 12.8 (19/20) t 7.51  +1.14 +/- 0.59 t 8.66  -0.022 +/- 0.008 (0/20) t -12.46   -0.028 +/- 0.071, ratio 0.975
    trigpoly     independent75pool4ruleGibbs32   1-20   +37.0 +/- 19.6 (20/20) t 8.45  +1.37 +/- 0.96 t 6.37  -0.023 +/- 0.008 (0/20) t -13.20   -0.026 +/- 0.070, ratio 0.977
    trigpoly     independent75pool4ruleGibbsB    21-40  +22.1 +/- 12.2 (19/20) t 8.10  +0.86 +/- 0.53 t 7.21  -0.026 +/- 0.012 (0/20) t -9.93    -0.016 +/- 0.045, ratio 0.986
    singleindex  independent75pool4ruleGibbsB    1-20   +14.4 +/- 13.2 (17/20) t 4.87  +0.73 +/- 0.47 t 6.94  +0.005 +/- 0.010 (13/20) t 2.38    +0.046 +/- 0.035, ratio 1.024
    singleindex  independent75pool4ruleGibbs32   1-20   +15.1 +/- 12.8 (18/20) t 5.29  +1.09 +/- 0.81 t 5.99  +0.009 +/- 0.010 (17/20) t 3.93    +0.042 +/- 0.036, ratio 1.022
    singleindex  independent75pool4ruleGibbsB    21-40  +23.9 +/- 11.0 (20/20) t 9.74  +1.14 +/- 0.86 t 5.91  +0.006 +/- 0.008 (18/20) t 3.51    +0.043 +/- 0.027, ratio 1.022

Held-out RMSE ratio, in the order of the paired table: 0.979, 0.980,
0.987, 1.028, 1.027, 1.023. Wall ratio: 3.12, 6.20, 3.05, 3.28, 6.86,
4.15. Paired standard error of the summed minimum ESS runs 2.5 to 4.4
across the six cells.

**This is the first arm on this cell to move the summed minimum ESS on
Trig+poly**, the gated mean function, and it moves it by two to five
times the +8 bar the perturb arms above were measured against, on both
mean functions and at both seed blocks; the per-chain minimum, which
every other four-chain arm here leaves at 2 on Trig+poly, moves with it. It costs 95 percent coverage on Trig+poly: -0.022 and
-0.026, past
[6.4 What "no regression on the core" means numerically](#64-what-no-regression-on-the-core-means-numerically)'s
-0.010 with the one-sided bound excluding the null in both blocks, and
the fresh block confirms the flag rather than dissolving it. The
interval is 14 percent shorter and the between-chain ratio falls from
0.78 to 0.58, so what narrows is the width pooling was supplying; the
control over-covers at 0.961 against a nominal 0.95 and arm B
under-covers at 0.939, the same miss with the sign reversed. Single
index, ungated, moves the other way on coverage (+0.005, +0.006) and
regresses held-out RMSE to 1.028 and 1.023, past the 1.02 margin at both
blocks. Full verdict at
[6. Benefit, pre-registered](nog-gibbs.md#6-benefit-pre-registered).

Chain overlap, `independent75pool4`: the median over the 25 ESS points of
the mean pairwise shared-over-spanned interval length is 0.48
(0.42-0.57) on Trig+poly and 0.58 (0.53-0.63) on Single index; every pair
of chains overlaps at every point of every seed. Host load ran 9 to 67
throughout, so these wall times carry no timing claim either.

**Level-fibre arms (2026-09-07).** level-fibre.md's slice 3 run adds two
more arms on the same shipped four-chain configuration and the same
twenty seeds: `independent75pool4level`, which varies no proposal
probability at all and only switches on the level-fibre Gibbs step
([1. The draw](level-fibre.md#1-the-draw)), and
`independent75pool4ruleGibbsBlevel`, that step on top of the rule_gibbs
arm above. `independent75pool4` is re-run in the same session as their
paired control and reproduces the Chain configuration paragraph's own
rows digit for digit on both mean functions. The level arm's Single
index wall ratio flagged, so it took a fresh-seed re-run on seeds 21 to
40, control and all, whose control rows reproduce the rule_gibbs
paragraph's own fresh-block rows digit for digit. Columns as above.

    mean fn      arm                                seeds  95% coverage        length  RMSE  min ESS (sum)  per chain  between
    trigpoly     independent75pool4                 1-20   0.961(0.945-0.977)  4.61    1.12  15(8-31)       2(1-2)     0.78
    trigpoly     independent75pool4level            1-20   0.963(0.930-0.976)  4.65    1.13  14(6-24)       2(1-2)     0.79
    trigpoly     independent75pool4ruleGibbsBlevel  1-20   0.942(0.914-0.960)  4.01    1.08  35(16-64)      2(2-3)     0.57
    trigpoly     independent75pool4                 21-40  0.964(0.950-0.977)  4.58    1.10  16(9-28)       2(1-2)     0.80
    trigpoly     independent75pool4level            21-40  0.964(0.951-0.981)  4.68    1.12  14(6-28)       1(1-2)     0.79
    singleindex  independent75pool4                 1-20   0.895(0.878-0.915)  6.45    1.92  21(9-32)       2(1-2)     0.68
    singleindex  independent75pool4level            1-20   0.896(0.880-0.920)  6.48    1.91  15(8-26)       2(1-2)     0.65
    singleindex  independent75pool4ruleGibbsBlevel  1-20   0.898(0.867-0.923)  6.59    1.98  36(17-63)      3(2-5)     0.39
    singleindex  independent75pool4                 21-40  0.891(0.868-0.911)  6.43    1.93  15(9-27)       2(1-2)     0.67
    singleindex  independent75pool4level            21-40  0.891(0.874-0.909)  6.46    1.93  17(9-25)       2(1-2)     0.70

Paired differences against each row's own control on its own seeds, mean
+/- sd (seeds positive of 20):

    mean fn      arm                                seeds  d min ESS (sum)                d per-chain min ESS     d 95% coverage                    d RMSE, ratio
    trigpoly     independent75pool4level            1-20   -0.9 +/- 6.4 (9/20) t -0.61    +0.01 +/- 0.19 t 0.22   +0.002 +/- 0.008 (13/20) t 1.14   +0.012 +/- 0.114, ratio 1.011
    trigpoly     independent75pool4level            21-40  -2.5 +/- 8.0 (9/20) t -1.41    -0.06 +/- 0.21 t -1.21  -0.000 +/- 0.007 (11/20) t -0.27  +0.026 +/- 0.069, ratio 1.024
    trigpoly     independent75pool4ruleGibbsBlevel  1-20   +19.9 +/- 15.4 (19/20) t 5.76  +0.63 +/- 0.43 t 6.57   -0.019 +/- 0.009 (0/20) t -9.47   -0.038 +/- 0.072, ratio 0.966
    singleindex  independent75pool4level            1-20   -5.3 +/- 8.2 (6/20) t -2.90    -0.01 +/- 0.28 t -0.23  +0.001 +/- 0.011 (10/20) t 0.53   -0.008 +/- 0.031, ratio 0.996
    singleindex  independent75pool4level            21-40  +1.7 +/- 6.6 (12/20) t 1.16    +0.07 +/- 0.25 t 1.23   -0.000 +/- 0.008 (9/20) t -0.20   +0.004 +/- 0.031, ratio 1.002
    singleindex  independent75pool4ruleGibbsBlevel  1-20   +15.0 +/- 13.4 (18/20) t 5.00  +0.86 +/- 0.66 t 5.86   +0.004 +/- 0.011 (12/20) t 1.45   +0.056 +/- 0.038, ratio 1.029

Held-out RMSE ratio, in the order of the paired table: 1.012, 1.014,
0.970, 0.994, 1.002, 1.033. Wall ratio: 0.966, 0.933, 3.629, 1.086,
0.962, 4.432. Paired standard error of the summed minimum ESS runs 1.43
to 3.45 across the six cells.

**The level arm moves nothing on this cell**, in any column. Its summed
minimum ESS reads -0.9 and -2.5 on Trig+poly against the +8 bar and
inside the sham's own -2.3 +/- 9.7, its per-chain minimum +0.01 and
-0.06, its coverage flat to the third digit, its held-out RMSE 1.012 and
1.014 against a 1.02 margin, and its between-chain ratio 0.79 against
0.78 and 0.79 against 0.80 - so this is not the perturb arms' pattern of
a between-chain ratio falling without ESS rising, it is a null
everywhere. Wall per sweep is the formality the design predicted: 0.966
and 0.933 in the same-session pairs, the one reading past 1.05 being
Single index at 1.086, and 0.986 by the means (0.915-1.024 per seed) on
a re-measure of the two arms alone at five seeds on Trig+poly with the
host at a load of 4.2 to 6.3. The stacked arm is the rule_gibbs arm to
within the pairing noise - +19.9 against that arm's +21.5, 35(16-64)
against 36(19-53), coverage -0.019 against -0.022, between 0.57 against
0.58, held-out RMSE 0.970 against 0.979 - so the two kernels do not
stack, and the coverage flag the structural kernel carries is unmoved by
the level step. Host load ran 4.2 to 9.7 through the arm runs, so those
wall times carry no timing claim beyond the quiet re-measure. Full
verdict at
[6. Benefit, pre-registered](level-fibre.md#6-benefit-pre-registered).

**Frozen-structure ESS (2026-09-07).**
[`C1-frozen-ess.R`](../../benchmarks/R/surfaces/C1-frozen-ess.R) separates
the leaf half of the kernel from the structural half on the recorded C1 arm
(Trig+poly, independent design, n = 10000, 75 trees, one chain, 1000
burn-in and 2500 kept): at the first five seeds, the structural chain's own
25-point minimum and median ESS and median lag-1 autocorrelation, against a
chain frozen (all four structural probabilities zero) at the last kept
draw and, separately, at the 1250th, then run for 2500 further draws with
no burn-in of its own. Each cell is minimum ESS / median ESS / median
lag-1 autocorrelation.

    seed    structural                 frozen at 2500              frozen at 1250
    1       1.2 / 10.8 / 0.71          19.6 / 841.9 / 0.32         4.1 / 794.9 / 0.39
    2       1.8 / 14.9 / 0.70          5.5 / 694.1 / 0.38          3.1 / 580.8 / 0.30
    3       1.9 / 16.2 / 0.68          25.1 / 623.3 / 0.33         4.2 / 532.6 / 0.40
    4       1.6 / 6.1 / 0.72           21.2 / 710.6 / 0.34         160.8 / 723.4 / 0.41
    5       1.4 / 28.2 / 0.63          162.7 / 603.3 / 0.37        81.3 / 646.8 / 0.45
    median  1.6 / 14.9 / 0.70          21.2 / 694.1 / 0.34         4.2 / 646.8 / 0.40

At a typical point the deficit is overwhelmingly structural: median-point
ESS rises from 14.9 unfrozen to about 670 of 2500 frozen, and the lag-1
autocorrelation drops from 0.70 to 0.34-0.40. At the worst point it is not:
minimum ESS rises only from 1.6 to 4-21, still three orders below the 2500
kept, and the worst frozen point carries two to three times the median
posterior spread, so the leaf Gibbs itself is slow at that coordinate and
[16.3 Ranking](tree-mixing-proposals.md#163-ranking)'s row 1 (the exact draw
on the level fibre) bears on the ranked (minimum-ESS) statistic rather than
the typical one. The two frozen start points agree in order of magnitude
throughout.

**Frozen-structure ESS with the level step (2026-09-07).** The same script
at its `levelGibbs` flag branches a paired second frozen chain at each
freeze point, identical to its partner but for the level-fibre Gibbs step
([1. The draw](level-fibre.md#1-the-draw)). The step is fixed when a
sampler is created, so it is enabled only AFTER the freeze: the arm is
a fresh sampler carrying the flag with the recorded chain's stored
state transplanted in, refused unless its trees, leaf values and sigma
are identical to the recorded sampler's, so both arms stand at one
forest when their chains begin. Its generators ride the state, so it
starts on the stream the recorded continuation starts on and parts from
it only as the step consumes draws. The off columns are the run above, reproduced digit
for digit.

    minimum ESS of 2500 kept
                frozen at 2500                  frozen at 1250
    seed    off      on       on - off      off      on       on - off
    1       19.6     293.0    +273.4        4.1      244.9    +240.8
    2       5.5      194.9    +189.4        3.1      249.9    +246.8
    3       25.1     142.4    +117.3        4.2      235.8    +231.6
    4       21.2     241.2    +220.0        160.8    273.4    +112.6
    5       162.7    183.7    +21.0         81.3     247.0    +165.7
    median  21.2     194.9    +189.4        4.2      247.0    +231.6

Each column's median is taken down that column, so the last entry of a
difference column is the median of the paired differences and not the
difference of the medians. The paired median rise is +189.4 at the 2500
freeze (bootstrap SE of the median 58.3; paired mean +164.2, SE 43.8) and
+231.6 at the 1250 one (SE 39.1; mean +199.5, SE 26.1), with all ten pairs
positive. The median point is unharmed: its frozen ESS medians 715.1 and
743.0 with the step on against 694.1 and 646.8 without, the worst of the
ten arms 632.7. What changes at the worst point is its character as well
as its ESS - the sd at the minimum over the median falls from 2.49 and
2.45 to 0.96 and 1.27, so the ranked coordinate is no longer the
high-spread one. Read against the paragraph above, the leaf Gibbs's own
slowness at the worst coordinate is a level-fibre deficit and the exact
draw on that fibre removes most of it; this is the frozen chain's residual
channel alone, and the structural channel it cannot see is
[6. Benefit, pre-registered](level-fibre.md#6-benefit-pre-registered)'s
slice 3. Host load ran 5 to 6 throughout, so these wall times carry no
timing claim: 12.4 to 13.3 s per seed for the recorded fit and its two
frozen chains, 2.9 s for each chain of the paired arm.

### 10.5 What the four cells say about the kernel that was measured

Facts only, against the rule in section 6.1.

- **The pilot establishes levels, not verdicts.** The rule accepts a change
  only on a paired contrast, and only P2 was run paired. Nothing here
  accepts or refuses anything.
- **P2 is the first design in this house where two shipped move-set arms
  separate at all.** Section 13 ran three mixtures over five Friedman-family
  designs and found no separation on any of seven metrics; section 13.7
  named the missing statistic as a structural one read between chains, and
  that is exactly the statistic that separates them here. What separates
  them is not the pathology - neither arm escapes the confounded design's
  two modes - but the null control, where the shipped mixture switches
  representation 78.5 times per chain and birth/death only switches zero.
- **The move set is not the missing term on the pathology itself.** Both
  arms lock. That extends section 13.6's reading of Pratola's low-noise
  collapse from coverage to representation: restoring or removing change and
  swap does not repair either.
- **P5 is currently a null.** It has an exact oracle and the shipped kernel
  answers it in 20 of 20 seeds. Theorem 5.2's premise is `pi_c = 0` and the
  default is not that, so this is the outcome the theory allows; but a cell
  the control arm passes cleanly cannot discriminate until some arm fails
  it. Its live statistic is the between-chain spread at 1.7x the mixing
  null.
- **The largest reproduced deficit is C1's coverage, and it is smaller than
  published.** 0.82 against a nominal 0.95 at n = 10000 and moderate noise,
  where He and Hahn report 0.73 to 0.74. Section 6.5 said that reproducing
  0.74 would establish in this house that the deficit is not a low-noise
  curiosity, and that failing to reproduce it would be the more interesting
  result. It is the second: the deficit is real and it is not low-noise, but
  at 13 points rather than 21, and 10 of those 13 are bought back by raising
  the tree count to 200.
- **The one large within-package gap is P6's, and it is not a kernel gap.**
  On the same 200 data sets the shipped default returns bias 0.31 and
  coverage 0.59 where dbarts' own propensity-augmented `forest()` surface
  returns 0.10 and 0.89.
- **Two internal validity checks fired correctly**: P2's duplicate-column
  null under the shipped mixture, and P5's inclusion oracle. Both are
  cheap, both have exact answers, and both should stay attached to their
  cells.

### 10.6 What the pilot could not do

- **No paired contrast on P6, P5 or C1, at the time.** One arm each;
  section 6.4's twenty matched pairs per cell is a paired-contrast count
  and is not met by twenty matched seeds of a single arm. All three now
  carry one, `default` rerun beside the added arms so the contrast is
  paired within one session: P6 and P5 each took two more arms on their
  own matched seeds
  ([10.2 P6, the diagonal shelf with targeted selection](#102-p6-the-diagonal-shelf-with-targeted-selection),
  [10.3 P5, the checkerboard on an autocorrelated design](#103-p5-the-checkerboard-on-an-autocorrelated-design)),
  and C1 took four in its own four-chain configuration
  ([10.4 C1, the He and Hahn factorial](#104-c1-the-he-and-hahn-factorial)).
  P5's two flags are measured findings pending the fresh-seed re-run 6.1
  requires before either counts, which has not been run.
- **P1 was not run** at the time, so section 6.4's absolute gate - the
  n = 2000, sigma = 0.25 rung's 90 percent coverage sitting near 0.71 in the
  control arm - was not in force behind any verdict above; each cell's own
  published-number check and the two internal oracles stood in for it. It
  has since run (10.8) and reads 0.725, so the gate is in force from
  2026-09-07 on.
- **P6's prognostic function is a reconstruction** and cannot be made
  otherwise from the published record.
- **C1's predictor arm and tree count are inferred**, from RMSE agreement,
  not stated by the source.
- **No between-chain statistic on C1**, which ran one chain per fit, so the
  coverage deficit is not attributed to mixing here - only measured.
- **No fresh-seed re-run** of any cell, which section 6.1 requires before a
  flag counts. Nothing was flagged, because nothing was contrasted.
- **Eight cells remain unbuilt**: P1, P3, P4, P7, P8, C2, C3 and C4.
- **Wall times are indicative only.** The host carried a load average
  between 6 and 136 across the run, none of it this measurement's.

### 10.7 Provenance

```
repo          /Users/vdorie/Repositories/dbarts, worktree on bartcore
measured at   fdba3809
build         private library installed from the worktree; dbarts 1.0.0,
              R 4.6.1, posterior 1.7.0, arm64 macOS, single-threaded
grid          P2   2 designs x 2 arms x 5 seeds x 8 chains          7 s
              P5   20 seeds x 8 chains                          5 m 19 s
              P6   3 reconstructions x 200 reps, plus a BCF arm  13 m 22 s
              C1   2 mean functions x 4 arms x 20 seeds             43 m
scope         measurement only - no source change, no default change,
              nothing scheduled. Sections 1 to 9 unchanged.
sources       every generating process was re-read from the primary source
              in this arc before the run, not taken from sections 1 to 5.
              The P6 propensity correction comes from the arXiv LaTeX
              source of 1706.09523; the P6 mu constraints from that
              paper's figures, read directly
scripts       benchmarks/R/surfaces, one per cell plus surfaces-common.R
              and a README; results written outside the working tree
```

### 10.8 P1, the low-noise Friedman emulator (2026-09-07)

Built and run after the pilot, as the battery's known-positive control
(6.3) and absolute gate (6.4). `P1-friedman.R` runs two rungs on twenty
matched seeds, one chain per fit, three arms through `proposal.probs`:
`default` (the shipped mixture, birth_death 0.6, swap 0, change 0.4),
`birthdeath` (Pratola's own arm) and `swap` (the former default, swap 0.1).
The house rung is section 13's cell, n = 2000, 200 trees, sigma = 0.25,
1000 burn-in and 2000 kept; Pratola's rung is n = 5000, 200 trees, 5000
burn-in and 5000 kept at sigma^2 in {1, 0.1}. Pratola's text prints
`10 sin(2 pi x1 x2)` where Friedman's function has `pi`, so his rung runs
his frequency, with one diagnostic arm (`variance01friedman`) at Friedman's;
his y is the function plus noise, not the deterministic simulator, and his
coverages are read as in-sample. In the published Bayesian Analysis text
the emulator is section 2.3 and the confounded step function 2.2, the
reverse of the arXiv numbering section 1.2 cites. Mean over seeds (min-max);
90 percent pointwise coverage of the true f.

    house rung (n 2000, sigma 0.25)
    arm         train                held-out             RMSE   min ESS
    default     0.782(0.743-0.812)   0.725(0.682-0.760)   0.260   1
    birthdeath  0.780(0.755-0.813)   0.709(0.657-0.737)   0.268   2
    swap        0.787(0.752-0.818)   0.728(0.684-0.775)   0.259   1

    Pratola's rung (n 5000)
    design              arm         train                held-out             RMSE   min ESS
    variance1           default     0.802(0.764-0.816)   0.791(0.753-0.830)   0.487   2
    variance1           birthdeath  0.799(0.761-0.824)   0.784(0.739-0.822)   0.488   2
    variance1           swap        0.805(0.779-0.842)   0.797(0.754-0.833)   0.481   2
    variance01          default     0.675(0.649-0.713)   0.622(0.569-0.681)   0.291   2
    variance01          birthdeath  0.674(0.637-0.708)   0.622(0.573-0.661)   0.292   2
    variance01          swap        0.667(0.633-0.689)   0.612(0.569-0.650)   0.295   2
    variance01friedman  birthdeath  0.710(0.671-0.734)   0.680(0.615-0.732)   0.207   2

**The absolute gate is in force.** Section 6.4 requires the house rung's
90 percent coverage to sit near 0.71 in the control arm; section 13 read
0.714, 0.725 and 0.714 held-out for its three mixtures, and the shipped
default now reads 0.725 (0.682-0.760) held-out on the same design. The
three mixtures remain indistinguishable, as section 13 found.

**Per-move acceptance**, from `P1-friedman-census.R` on a
`-DBARTCORE_MOVE_CENSUS` build, one seed, 200 sampled sweeps after each
rung's burn-in, as a share of proposals made: sigma^2 = 1, birth/death
arm, birth 14.3 / death 17.1 / all 15.6 percent (published about 18);
sigma^2 = 0.1, birth/death arm, birth 5.1 / death 5.6 / all 5.3 (published
about 4); sigma^2 = 0.1, default arm, birth 5.7 / death 6.0 / change 2.1 /
all 4.3; the swap arm's swap move accepts 1.0 percent at 74 percent no-op;
house rung, default arm, birth 10.2 / death 11.8 / change 4.7 / all 8.4.

| statistic (birth/death arm) | published | measured | verdict |
|---|---|---|---|
| acceptance, sigma^2 = 1 | ~18% | 15.6% | reproduces |
| acceptance, sigma^2 = 0.1 | ~4% | 5.3% | reproduces |
| 90% coverage, sigma^2 = 1 | 0.81 | 0.799 train, 0.784 held-out | reproduces |
| 90% coverage, sigma^2 = 0.1 | 0.538 | 0.674 train, 0.622 held-out | direction reproduces; the collapse is 0.13 shallower here |

The frequency transcription does not explain the shallower collapse:
Friedman's own frequency reads 0.710 train at RMSE 0.207, no closer to
0.538 than the doubled one. Wall times are not recorded as measurements;
the host carried a 1-minute load of 12 to 18 throughout.
