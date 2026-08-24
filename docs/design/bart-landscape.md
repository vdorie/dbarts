# The BART implementation landscape

Snapshot: 2026-08-12. In-repo anchors were read in `wt/zero-weight` between
`b2290614` and `797d1527` (the tip advanced mid-write when `feature-matrix.md`
landed; no anchor cited here is touched by the difference). DESCRIPTIVE only -
no ranking, no scheduling, no recommendation, and no maintenance promise: this
is a dated photograph.

Companion to two documents, duplicating neither.
`docs/plans/roadmap-survey.md` (Review 5, 2026-07-10) RANKED dbarts's backlog
by external demand; this document ranks nothing and is raw material a later
ranking pass can cite, one implementation at a time (two of that survey's
figures are corrected below from measurement taken today: embarcadero is no
longer a CRAN reverse dependency, and SoftBart's binary `k = 1` is now
sourced). `docs/design/feature-matrix.md` is the INTERNAL-completeness
companion, cell by cell; this is the EXTERNAL one. The rubric comes from
`docs/design/r-c-division.md` (ACCEPTED 2026-08-11).

## Verification protocol

Every claim about an external artifact was checked first-hand today: 27 CRAN
source tarballs downloaded and unpacked, six GitHub repositories cloned, two
papers fetched as PDF and converted to text (SoftBart's and stochtree's) plus
two abstracts read, PyPI JSON read directly, download counts pulled from
cranlogs, GitHub metadata through the API. Clones and tarballs
live in the session scratchpad, never in this repo. Each entry ends with a
`Verified:` line naming what was read; nothing is carried from memory. Claims
that could not be established first-hand are marked **UNVERIFIED** inline -
there are three, one in the artifacts section and two in the appendix. A
subsequent independent agent re-fetched every load-bearing source and audited
the document adversarially; its findings are folded in, and the entries it
changed say so on their `Verified:` line (pymc-bart, nftbart - whose grade
moved - nbbart, VCBART, and synthesis point 5).

Coverage was derived, not recalled: the CRAN `src/contrib` index (24,724
tarballs) grepped for `bart` and for each engine name, plus an r-universe
full-text search for "Bayesian additive regression trees" to catch packages
whose name lacks "bart" (that is how `ShrinkageTrees` and `dlmtree` surfaced),
plus PyPI and the GitHub search API for the Python and research tail.

CRAN traction, one call, total downloads 2026-07-13..2026-08-11: dbarts 6043,
bartMachine 3575, BART 2617, stan4bart 1806, BayesTree 906, bartCause 417,
stochtree 387, dlmtree 366, bcf 356, SoftBart 284, bartcs 281, flexBART 199,
ShrinkageTrees 178, riAFTBART 173, ridgeBART 161, nftbart 158, EBcoBART 157,
VCBART 135, uddbart 110, sparseVCBART 100.

## The rubric

For each artifact: what would its author have had to do to build it on a
public engine instead of writing tree code?

- **(A)** Composable on dbarts as it stands, from R, by driving a sampler
  between sweeps - `setResponse` / `setOffset` / `setWeights` / `setPredictor`
  / `setSigma` / `setModel`, zero-weight row subsetting, `run(0L, 1L)`,
  `getLatents`, warm starts. Name the recipe.
- **(B)** Composable once the committed adoption slate lands (r-c-division,
  "Adoption slate"): nameable calibration, the latent-family active-rows mask,
  the composition validator, augmentation helpers, the M4 basis family. Name
  the unlocking item.
- **(C)** Needs a CONTRIBUTED engine family or leaf model - an integrand
  change (likelihood family, link, leaf model, tree prior, hyperprior law,
  split-rule representation), which this architecture takes as a C++
  contribution selected from R (`facade.hpp` concepts).
- **(D)** Structurally engine-only in ANY implementation: a likelihood that
  does not factorise over observations given the forest, a move that writes
  engine state, a coupling not additive on a carriable quantity. A fork or a
  new engine was the honest choice. PLATFORM ports (another language or
  execution model) are graded (D) too, with "platform, not model" stated - a
  new engine was genuinely right, for a reason unrelated to the integrand.

"As it stands" means the surface in THIS worktree (the pre-1.0 bartcore line:
gaussian, Student-t, probit, logistic, ordinal, nbinom, multinomial, AFT,
hazard, hurdle, BCF, grouped, heteroscedastic - `feature-matrix.md`), not CRAN
0.9-33, which is gaussian plus binary probit (verified: its `R/A_class.R`
carries a single `binary` logical). Most authors below had only the CRAN line,
so each (A) says whether the recipe would have worked on what they had. Where
another public engine's surface would have served an author better than
dbarts's, the entry says so.

---

## R packages: engines

**BayesTree 0.3-1.5** (Chipman, McCulloch; CRAN, published 2024-01-30;
906/mo). The original reference implementation of Chipman, George and
McCulloch (2010). *Unique:* it is the origin - there was nothing to reuse.
*Built:* hand-written C++ with its own class hierarchy (`Node.cpp`,
`Prior.cpp`, `Likelihood.cpp`, `EndNodeModel.h`, `MuS.cpp`, `Sdev.h`), no
Rcpp, `Imports: nnet`. *Counterfactual:* none - it is what the others reuse or
replace; dbarts began as a rewrite of exactly this code. *Verified:* CRAN
tarball DESCRIPTION and `src/` listing.

**BART 2.9.10** (Sparapani, McCulloch, Spanbauer; CRAN, 2026-01-24; 2617/mo).
The model-suite package: `wbart`, `pbart`, `lbart` logistic, `mbart`/`mbart2`
multinomial, `abart` AFT, `surv.bart` / `recur.bart` / `crisk.bart` /
`crisk2.bart` survival, recurrent-event and competing-risk suites, `heterbart`
heteroscedastic. *Unique:* breadth of outcome types under one prior and one
code base, plus the person-period survival machinery (`surv.pre.bart`) that
made discrete-time hazard BART routine. *Built:* own C++ (`tree.cpp`,
`bd.cpp`, `bartfuns.cpp`, headers marked "Copyright (C) 2017 Robert McCulloch
and Rodney Sparapani"), a rewrite in the BayesTree line rather than a reuse of
it; Rcpp plus OpenMP for the `mc.*` entry points. *Counterfactual:* **(A) for
most of the suite.** The survival, recurrent and competing-risk models are a
data reshape plus probit BART - the reshape is R, the fit is a stock probit
sampler, and the tip ships that reshape itself as the `hazard` row; `lbart`'s
logistic is latent augmentation, the shape `nbbart` demonstrates on dbarts
below. Two parts are not: `mbart`'s softmax margin needs a log-sum-exp over
the other K-1 forests, r-c-division's named non-carriable coupling - **(D)** -
and `heterbart`'s variance forest is **(C)**, an integrand change dbarts has
since taken as one. *Verified:* NAMESPACE export lines 3-16, `src/` listing
and copyright headers, `R/lbart.R`, DESCRIPTION.

**dbarts (this package)** - self-placement. Among engines it is the one whose
distinguishing feature is not a model but a SURFACE: predictors, response,
offset, weights, sigma, treatment and (on the tip) model swap between sweeps,
so the forest sits inside somebody else's Gibbs sampler. Ten CRAN packages
import it and one (`stan4bart`) links to it - the largest reverse-dependency
footprint of any BART engine, at ~1.7x bartMachine's downloads. The 1.0 line
has absorbed much of the model space too (13 rows in `feature-matrix.md`,
including BCF, ordinal, negative binomial, AFT, grouped intercepts, variance
forests, and four leaf models - constant, monotone, linear, GP - a combination
no other engine carries; stochtree pairs constant with regression leaves, but
nobody else ships the GP one). Gaps against this landscape: no soft trees, no
graph-structured categorical splits, no horseshoe leaf prior, no JSON
serialization, no Python. *Verified:* CRAN reverse-dependency lists; this
worktree's `man/` and `feature-matrix.md`.

**bartMachine 1.4.2** (Kapelner, Bleich; CRAN, 2026-05-03; 3575/mo). BART on
the JVM. *Unique:* the substrate - Java, with `serialize`, multithreaded
chains, tree illustrations - plus a built-in inference layer (permutation
variable selection, `cov_prior_vec`, `interaction_constraints`) and
missing-data handling (`use_missing_data`, MIA) years before other R engines;
`SystemRequirements` now names Java >= 21.0 and TornadoVM for optional GPU.
*Built:* from scratch in Java (`inst/java/bart_java.jar`, via
`bartMachineJARs` + `rJava`); no C++ at all. *Counterfactual:* **(A) on the
model, (D) on the platform.** Nothing in the model leaves a stock BART
integrand - MIA is the dbarts DEFAULT today (`missing = "incorporate"`), and
permutation variable selection is an R layer over any engine's `varcount`.
What the authors wanted was a JVM-native artifact with its own memory and
threading story, which no R/C++ engine supplies. *Verified:*
`R/bart_package_builders.R` formals, DESCRIPTION, `inst/java` listing; dbarts
`man/dbarts.Rd:123-124`.

**SoftBart 1.0.3** (Linero; CRAN, 2025-11-23; 284/mo). Soft BART: leaves are
reached probabilistically through logistic gates rather than by hard
partition. *Unique:* two things - the soft trees, and `MakeForest()`, an
`Rcpp_Forest` module with `do_gibbs(X, Y, X_test, i)`,
`do_gibbs_weighted(X, Y, weights, X_test, i)`, `set_sigma`, `get_sigma_mu`,
`set_s`/`get_s`, `do_predict`, which is the only R surface besides dbarts's
built to be driven from somebody else's Gibbs loop. The package paper states
the ambition in as many words: "as to the best of my knowledge SoftBart is the
only BART package that allows users to embed BART within a larger model
without having to modify the underlying C++ code" (arXiv 2210.16375, "Our
Motivation") - false as to dbarts, which that same paper cites. *Built:* own
C++ (`soft_bart.cpp/.h`, RcppArmadillo) including its own slice sampler.
*Counterfactual:* **(C).** Soft routing changes the integrand: every
observation contributes to every leaf with a weight, so leaf sufficient
statistics stop being within-leaf sums. Note also that
`Hypers()` NAMES the leaf-prior scale (`sigma_mu <- 0.5 / (k * sqrt(num_tree))`,
`R/SoftBart.R:42`) and the probit entry overrides it
(`hypers$sigma_mu = 3 / k / sqrt(num_tree)`, `k = 1` default,
`R/softbart_probit.R:79,111`) - a peer engine already ships a nameable
calibration, on the internal scale, with the user responsible for normalizing
X and Y. *Verified:* `R/MakeForest.R`, `src/soft_bart.h` class `Forest`,
`R/SoftBart.R`, `R/softbart_probit.R`; arXiv 2210.16375 PDF converted to text.

**bcf 2.0.2** (Hahn, Murray, Carvalho, Starling; CRAN, 2024-02-27; 356/mo).
Bayesian causal forests: separate prognostic and treatment forests, the
estimated propensity as a covariate (targeted selection), and an amplitude
parameterization. *Unique:* the asymmetric prior pair - `ntree_control = 200`,
`base 0.95 / power 2`, `sd_control = 2*sd(y)` against `ntree_moderate = 50`,
`base 0.25 / power 3`, `sd_moderate = sd(y)` - plus `use_muscale` /
`use_tauscale`, the `bcf_clean_overpar.cpp` half-Cauchy scale
parameterization. *Built:* own C++ in the Hahn/Murray line (`tree.cpp`,
`bd.cpp`, `funs.cpp`, `TreeSamples.cpp`), RcppParallel, own logging and RNG.
*Counterfactual:* **(D), but only for the amplitude move.** The additive
two-forest structure is (A) - two dbarts samplers each carrying the other's
fit as an offset target the same posterior, measured twice in the r-c-division
arc. What does not port is the `b0`/`b1` rescale, which reinterprets
already-drawn leaf values and so writes engine state. dbarts has since taken
BCF as an engine model for exactly that reason. *Verified:* `R/bcf.R:224-242`
formals, `src/` listing, `tree.h` header.

**stochtree 0.4.5** (Herren, Hahn, Murray, Carvalho; CRAN, 2026-06-30;
387/mo; also PyPI). A C++ library with R and Python bindings for BART, XBART,
BCF and XBCF. *Unique:* three things no other package pairs. (i)
grow-from-root as the headline sampler, not an option (`num_gfr = 5` default
in `bart()`); (ii) a genuinely low-level R6 surface - `ForestDataset`,
`Outcome` with `add_vector` / `subtract_vector` / `update_data`,
`ForestModel$sample_one_iteration(..., gfr = TRUE)`,
`sampleGlobalErrorVarianceOneIteration`, `resetActiveForest` - so users write
their own outer sampler against the same objects `bart()` uses; (iii) JSON
serialization with fits passing between R and Python. Leaf models are an enum,
not a hard-coded model: `kConstantLeafGaussian`,
`kUnivariateRegressionLeafGaussian`, `kMultivariateRegressionLeafGaussian`,
`kLogLinearVariance`, `kCloglogOrdinal`. *Built:* new C++ core
(`src/include/stochtree/*`, vendored Eigen), cpp11 for R. *Counterfactual:*
**(D), by design** - the library IS the deliverable; grow-from-root writes
engine state, and the two-binding goal cannot be met by a single-language
engine. Its own paper's Table 1 gives the "C/C++ API" row exactly one check
(its own) and the "Random effects" row exactly one, and its prose calls dbarts
"a limited interface for interchanging a forest MCMC step with other
samplers"; dbarts ships `inst/include/dbarts/dbarts.h` (39 entries in
`DBARTS_C_API_LIST`) and `rbart_vi`. *Verified:* `R/data.R:229-318`,
`R/model.R:93-133`, `R/bart.R:187-207`,
`src/include/stochtree/leaf_model.h:352-358`, NAMESPACE; arXiv 2512.12051v1
Section 1.2 and Table 1 fetched and converted to text.

**flexBART 2.0.3** (Deshpande, Perrett; CRAN, 2026-02-12; 199/mo). BART with
first-class unordered categorical predictors. *Unique:* a categorical split
sends an arbitrary SUBSET of levels left, and when levels carry a network the
subset is drawn by partitioning a random spanning tree of that network
(`src/graph_funs.h`: `wilson`, `find_components`, `fiedler_split`,
`floydwarshall`) - the split rule itself is a random graph partition. Version
2 added a formula interface with `family`/`link` parsing (gaussian, binomial
probit and logit, Poisson), heteroscedastic mean-and-variance ensembles, and
varying-coefficient models. *Built:* own C++ with RcppArmadillo, one fit file
per model shape (`single_ensm_probit.cpp`, `multi_ensm_poisson.cpp`, ...).
*Counterfactual:* **(C)** - a split-rule representation change is clause 2 by
name. The plain categorical half is no longer a gap: this tip splits unordered
factors by level subset natively (`factors = "categorical"` is the DEFAULT,
`man/dbarts.Rd:105`); the graph-structured half remains a contribution nobody
else has. *Verified:* `R/flexBART.R:1-30`, `src/graph_funs.h`, `src/structs.h`,
DESCRIPTION.

**nftbart 2.3** (Sparapani, McCulloch, Pratola; CRAN, 2025-12-03; 158/mo).
Nonparametric failure time BART: `Y = mu + f(x) + sd(x) E`, `f` a BART forest,
`sd` an HBART variance forest, `E` a Dirichlet-process mixture (Low
Information Omnibus prior). *Unique:* a nonparametric error distribution
joined to a heteroscedastic forest. *Built:* own C++ in the
Pratola/McCulloch/Chipman OpenBT line (`brt.h`, `ambrt.h`, `psbrt.h`,
`DPMLIOneal8.h` implementing Neal's algorithm 8), header-only, Rcpp.
*Counterfactual:* **(C).** The pieces exist - a variance forest and a weights
channel coexist on the tip - but the composition does not close. nftbart's own
R reconstructs the observation-level moments as
`mu. = dpmu * s.train + f.train` and `sd. = dpsd * s.train`
(`R/nft2.R:425-426`), so given a DPM component the mean forest's conditional
carries offset `dpmu_c * s(x_i)` and precision `1 / dpsd_c^2` ON TOP of the
`s(x_i)^2` the variance forest already supplies - both scaled by a quantity
the engine REDRAWS inside `run()`. A host driving the loop from R can only
install a one-sweep-stale `s(x)`, so the recipe is approximate, not an exact
Gibbs kernel; making it exact means the scaled DPM location entering within
the sweep, which is a contributed response family. (The DPM block itself is
C++ in nftbart - `DPMLIOneal8.h`, `drawDPM` passed in at `R/nft2.R:292`; it is
the counterfactual port's R program that would draw it.) *Verified:*
DESCRIPTION Description field, `src/` listing and copyright headers,
`R/nft2.R:292,425-426`; the offset-scaling defect was found by an independent
adversarial re-fetch of the package, which is why this grade moved from (A).

**bartcs 1.3.0** (Yoo, Kim; CRAN, 2025-04-08; 281/mo). Confounder selection:
an exposure model and an outcome model fit jointly, sharing one Dirichlet
prior over splitting probabilities so a variable used by one is more likely to
be offered to the other. *Unique:* the coupling runs through the SPLIT PRIOR,
not the mean - `single_bart.cpp:63,98` draws `var_prob` from `rdirichlet` and
updates it by Metropolis against the union of both models' variable counts.
*Built:* own C++ (`BART.cpp`, `Node.cpp`, `grow.cpp`, `prune.cpp`,
`change.cpp`), Rcpp, calling back into MCMCpack from C++. *Counterfactual:*
**(B).** `$setModel` already accepts a fresh `cgm(power, base, split.probs)`
mid-run, so a two-sampler R loop that reads `varcount`, draws the Dirichlet
and re-installs the probabilities is expressible today - but `setModel` also
re-pins sigma (engine-side, not in the R layer: `Chain::setModel` sets
`sigmaIsFixed_` and calls `setSigma` or `setSigmaPrior` for gaussian and aft
samplers with no variance forest, `src/bartcore/chain.hpp:1304-1314`) and is
refused outright on a DART prior, so the composition carries side effects its
author never named. Unlocking items: a
split-probability channel narrower than `setModel` (dbarts issue #67,
"Feature request: more convenient updates of splitting probabilities", open
since 2024-03-14) and the validator. *Verified:* `R/separate_bart.R:4-21`,
`src/single_bart.cpp`, `src/separate_bart.cpp`; this worktree's
`R/dbarts.R:1229-1238` and `man/dbartsSampler-class.Rd:103`.

**VCBART 1.2.5** (Deshpande, Bai, Balocchi, Starling, Weiss; CRAN,
2026-04-21; 135/mo). Varying-coefficient models `y = sum_j beta_j(x) z_j`,
one forest per coefficient function. *Unique:* per-forest BASES - each
ensemble's output multiplies its own modifier rather than entering additively.
*Built:* own C++ (`update_trees.cpp`, `vcbart_ind_fit.cpp`,
`vcbart_cs_fit.cpp` for a compound-symmetry error), RcppArmadillo.
*Counterfactual:* **(A) for the independent-error arm**, by a
weighted-least-squares change of variable that is an IDENTITY, not an
approximation. Block `j`'s conditional is
`(y - sum_{l != j} beta_l z_l) = beta_j(x) z_j + e`; divide by `z_j` and it is
an ordinary Gaussian regression on `beta_j(x)` with case weights `z_j^2`, so
the leaf sufficient statistics are `sum z_j^2` and `sum z_j r` - exactly the
two sums VCBART's own C++ accumulates (`src/funs.cpp:206-207`, inside
`compute_p_theta_ind`), and the leftover `sum r^2` is partition-independent
and cancels in the Metropolis ratio, so the tree moves agree too. Two clauses
the recipe needs. (i) `z_j == 0` rows: `resid / z_j` is NaN or Inf there, and
NaN is refused by `setResponse` while Inf silently poisons the leaf's weighted
mean - so write any FINITE placeholder at those rows, where the zero weight
makes the value irrelevant, and the shipped "excluded from the likelihood,
fitted value retained" semantics apply (`man/dbartsSampler-class.Rd:137`).
(ii) The compound-symmetry arm (`vcbart_cs_fit.cpp`, `update_rho`,
`compute_p_theta_cs`) is OUT of scope: its within-subject covariance is
non-diagonal and case weights are diagonal. The standing caveat is the
calibration wall - `beta_j`'s leaf prior is inherited from the range of
whatever `resid / z_j` the sampler was CONSTRUCTED on, so the composition is
correct only by care. Nameable calibration makes it **(B)** and correct by
construction; the M4 basis family would make it one sampler instead of K. The
recipe works on CRAN 0.9-33 too. *Verified:* DESCRIPTION, `src/` and `R/`
listings, `R/VCBART_ind.R`; `src/funs.cpp:206-207` and the zero-weight
failure modes come from an independent adversarial re-check.

**sparseVCBART 1.0.0** (Ghosh, Deshpande; CRAN, 2026-05-19; 100/mo). VCBART
with global-local (horseshoe-type) shrinkage across coefficient functions.
*Unique:* the hyperprior law on the leaf scales, not the mean structure.
*Built:* a sibling code base to VCBART, own C++. *Counterfactual:* **(C)** - a
hyperprior law over leaf scales is clause 2, and the leaf draw changes.
*Verified:* DESCRIPTION and `src/` listing, side by side with VCBART's.

**ridgeBART 1.0.2** (Yee, Deshpande; CRAN, 2026-05-27; 161/mo). BART whose
leaves output ridge functions of random features rather than constants.
*Unique:* the leaf model - `rand_basis_funs.cpp`, `get_bases.cpp`,
`eval_bases.cpp` draw random feature maps per leaf, so the fit is smooth
within a leaf. *Built:* own C++ (flexBART's file layout, same maintainer),
RcppArmadillo. *Counterfactual:* **(C).** dbarts's `linear()` and `gp()`
leaves are the same genus, and a FIXED feature map supplied as designated leaf
covariates would come close - but ridgeBART redraws the bases inside the
sampler, an integrand change. *Verified:* DESCRIPTION, `src/` listing.

**ShrinkageTrees 2.0.2** (Jacobs et al.; CRAN, 2026-04; 178/mo). Tree
ensembles with horseshoe-type shrinkage on the leaf parameters, aimed at
survival and causal inference (`CausalHorseForest`). *Unique:* the leaf
shrinkage prior in a high-dimensional survival setting. *Built:* own C++
(`HorseTrees.cpp`, `Forest.cpp`, `ForestEngine.h`, `probit-HorseTrees.cpp`,
`OuterGibbsFunctions.cpp`), Rcpp. *Counterfactual:* **(C)**, sparseVCBART's
class. *Verified:* CRAN tarball DESCRIPTION, `src/` and `R/` listings.

**monbart** (Murray; GitHub `jaredsmurray/monbart`, DESCRIPTION dated
2022-04-15, packaged 2018; not on CRAN). BART for a binary outcome constrained
monotone in one binary covariate. *Unique:* the constraint couples two
surfaces - the fitted probability under treatment must dominate the one under
control - a truncation of the joint leaf draw, not a per-forest constraint.
Motivating application is a bankruptcy-forecast sensitivity analysis (README
cites Papakostas, Hahn, Murray, Zhou and Gerakos, arXiv 2106.04503). *Built:*
own C++ sharing the bcf lineage's file names (`bd.cpp`, `funs.cpp`,
`tree.cpp`, `TreeSamples.cpp`, `slice.cpp`) plus `bartRcppMono.cpp`.
*Counterfactual:* **(C).** dbarts ships monotone leaves, but for monotonicity
in a NUMERIC predictor within one forest; a cross-forest dominance constraint
on a probit latent is a different constrained draw. *Verified:* cloned the
repository - DESCRIPTION, README.md, `src/` listing.

## R packages: consumers and layers

These write no tree code; they are evidence about the surface, not the
integrand, and are graded only where they drive a sampler iteratively.

**bartCause 1.0-10** (Dorie, Hill; 417/mo) - causal estimands over dbarts
fits; in-house, and the home r-c-division assigns `bcf()`. **stan4bart 0.0-12**
(Dorie; 1806/mo) - the only CRAN package that links to dbarts: a BART forest
and a Stan-sampled parametric block alternating through the offset channel,
`mvbart()` routing around the response channel deliberately; in-house, and the
reference implementation of the (A) pattern. **riAFTBART 0.3.3** (Hu, Ji,
Zhang; 173/mo) - random-intercept AFT with multiple treatments, and the
measured defect: `R/sampling_funs.R:13-17` calls `dbarts::bart(...)` FRESH
every Gibbs iteration and takes `mod$yhat.train.mean` and `mean(mod$sigma)` as
the draw. Posterior means used as draws, plus a cold-start refit per sweep.
**(A)** - the correct version is one persistent sampler,
`setResponse(logT - b[cluster])`, `run(0L, 1L)`, `$train` and `getSigmas()` -
and it is exactly the defect class the composition validator exists for.
**EBcoBART 1.1.2** (Goedhart et al.; 157/mo) - empirical-Bayes co-data
learning: refits `dbarts::bart` in an outer loop and reads
`dbarts::extract(fit, "trees")` to re-estimate split probabilities and the
tree prior. **(A)**, and legitimately refit-based (hyperparameter estimation,
not a Gibbs sweep). **bmabart 2.0** (Yu, Li) - mediation analysis on BART.
**cjbart 0.3.2** (Robinson, Duch) - heterogeneous effects in conjoint
experiments, on BART. **uddbart 0.2.0** (Pan) - interval-censored dynamic risk
prediction, `NeedsCompilation: no`, BART only in `Suggests`. **bartMan 0.2.1**
(Inglis et al.) and **bartXViz 1.0.11** (Lee, Lee) - visualization and SHAP
layers reading fits from BART, dbarts and bartMachine alike; the only
artifacts here that treat three engines as interchangeable. *Verified:* each
package's CRAN DESCRIPTION plus the cited source lines from the unpacked
tarballs.

## Python

**pymc-bart 0.12.0** (PyMC developers; PyPI 2026-05-13; 154 GitHub stars, the
most-starred actively maintained BART artifact - `JakeColtman/bartpy` has 232
but was last pushed 2023-12-22). *Unique:* BART as a random VARIABLE inside a
probabilistic programming language - `pmb.BART('bart', X, y)` inside a
`pm.Model()` block, composable with any other PyMC distribution, trees sampled
by a dedicated particle-Gibbs step. *Built:* the RELEASED 0.12.0 is pure
Python on numba - the sdist ships `pymc_bart/pgbart.py`, `tree.py` and
`split_rules.py`, its `__init__.py` reads `from pymc_bart.pgbart import
PGBART`, and `requires_dist` names `numba` and not `bartrs`. A Rust rewrite
has landed on unreleased `main` (commit `f743d904`, 2026-06-11, "Rust sampler
(#276)"), where the package becomes a thin layer over **bartrs** imported for
its side effect (`import bartrs  # registers PGBART with PyMC.`); `main` still
self-reports `__version__ = "0.12.0"`, so reading HEAD makes the Rust layer
look shipped when `pip install pymc-bart` today gets the numba sampler.
*Counterfactual:* **(D), platform** - no R engine can be a node in a PyTensor
graph. The direction of travel is worth recording, softly: the Python flagship
is moving its sampler out of Python into a compiled core, arriving where the R
packages already were - though `GStechschulte/bart-rs` (2024-07-14) carries
the same description, so `pymc-devs/bartrs` (2026-05-08) dates the repository,
not the first Rust sampler. *Verified:* PyPI JSON and the 0.12.0 sdist
contents, `main`'s `__init__.py` and `requirements.txt`, `pymc-devs/bartrs`
metadata, README - re-checked against an independent adversarial re-fetch that
REFUTED an earlier draft attributing the Rust sampler to the release.

**bartz 0.12.0** (Petrillo; PyPI 2026-07-24; 96 stars). BART in JAX.
*Unique:* the execution model - the sampler is vectorized and jit-able, so it
runs on GPU; the README claims up to 200x on an A100 at 2M+ observations, and
CPU parity with dbarts ("the fastest implementation I know of") at lower
memory. *Built:* from scratch in JAX. Its module layout is itself a landscape
observation: `bartz.mcmcstep` exposes `init` / `step` over explicit `State`
and `Forest` objects (a functional step API a custom sampler can drive),
`bartz.BART` provides `gbart` / `mc_gbart` "that mimic the R BART3 package",
and `bartz.stochtree` is a "`stochtree`-compatible interface to bartz" - one
artifact carrying two rivals' API shapes. *Counterfactual:* **(D), platform** -
a data-parallel accelerator rewrite cannot borrow a pointer-chasing C++ tree
sampler. *Verified:* PyPI JSON, README, three submodule docstrings from raw
githubusercontent.

**stochtree (Python) 0.4.5** - the same C++ core as the R package, same
version, same low-level objects. *Counterfactual:* **(D), platform**, and the
reason the core-plus-two-bindings design exists. *Verified:* PyPI JSON
(documentation `stochtree.ai/python_docs`, repository
`StochasticTree/stochtree`).

**xbart 0.1.8** (He, Yalov, Hahn; PyPI, last released 2020-03-02). The
original accelerated-BART (grow-from-root) implementation, R and Python faces
from one repository; effectively superseded by stochtree, whose authors
overlap - its README still says "The CRAN version will be submitted soon".
*Counterfactual:* **(D)** - a tree-construction algorithm that writes engine
state; dbarts has since taken it as `growFromRoot`. *Verified:* PyPI JSON,
GitHub README and repository listing.

**BART-Survival** (Tiegs et al., CDC; GitHub `CDCgov/BART-Survival`, Python).
Discrete-time survival with BART. *Unique:* nothing in the integrand -
person-period expansion plus BART, packaged for public-health users. *Built:*
**reuse** - `pyproject.toml` reads "Survival analyses with Bayesian Additive
Regression Trees using PyMC-BART as BART backend", with pymc-bart a hard
dependency. *Counterfactual:* **(A)**, and the clearest instance in the
landscape of an author treating a BART engine as a library. *Verified:*
`pyproject.toml` through the GitHub API.

## Research artifacts and forks

**princeBART 0.2.0** (Godoy Garraza; GitHub `AlkemaLab/prince_BART`).
Principal stratification with BART for a binary instrument and binary outcome
(LATE / complier effects). *Unique:* the row set each component is fit to is
LATENT and moves every outer iteration - compliers, always-takers and
never-takers are imputed per sweep. *Built:* reuse of dbarts, awkwardly. Six
probit samplers in `initialize_samplers` (`R/fit_psbart_binary.R:257-302`),
five carrying a `subset =` stratum; `update_samplers` calls
`$setData(dbartsData(X, Y, subset = <new stratum>))` on five of them EVERY
outer iteration (`:339-386`), for `n_warmup + n_samples` = 2000 iterations by
default (`:35-46`). Nine dbarts internals are reached by `getFromNamespace`
(`R/utils.R:145-156`, plus `normal` again at `R/fit_psbart_ordinal.R:385`), of
which four are used - to hand-build a probit sampler with a chosen node prior,
because the convenience layer is not reusable. The row subsetting itself uses
the PUBLIC `dbartsData(subset =)`. *Counterfactual:* **(B).** Gaussian row
subsetting already ships through zero weights, but `setWeights` is refused for
probit on identification grounds, so the probit path has no per-sweep row
channel; the unlocking item is the latent-family active-rows mask
(`setActiveRows`, v1 scope gaussian, Student-t, probit, ordinal). With it,
princeBART writes one sampler per outcome model and one
`$setActiveRows(co == 1 & Z == 1)` per iteration, keeping its trees, its
latents and its cut grid instead of rebuilding all three. Its ordinal path
already needs nothing - it drives GAUSSIAN samplers on host-drawn latents,
subset by an OBSERVED instrument that never moves. *Verified:* cloned the
repository and read the cited files; agrees with the first-hand read in
`.claude/latent-subset-mask-design/memo.md` Part 2.

**clbart 0.1.0** (Englert; GitHub `jacobenglert/clbart`). Heterogeneous
exposure effects in the case-crossover design: a conditional logistic
regression whose exposure coefficient is a BART function of time-invariant
covariates. The method paper is Englert, Ebelt and Chang, arXiv 2311.12016,
whose abstract says "CL-BART utilizes reversible jump Markov chain Monte Carlo
to bypass the conditional conjugacy requirement". *Unique:* the likelihood is
stratum-conditional, so it does not factorise over observations - it
factorises over STRATA, whose rows can land in different leaves - and the leaf
parameter is not conjugate. *Built:* the tree moves are PURE R
(`R/tree_proposals.R` implements grow/prune/change as data-frame surgery,
`R/update_forest.R` runs backfitting, `R/update_mu.R` draws leaves); C++
appears only for the conditional-logistic log-likelihood (`src/clr.cpp`), a
log-gamma variate and a variance update. Node parameters come from a Laplace
approximation (`get_node_map`) inside a reversible-jump acceptance.
*Counterfactual:* **(D)** - the non-factorising case in the principle's own
list; no engine's leaf kernel can serve it. It is also the counterexample that
keeps clause 2 honest: the tree MOVES here are R, slow but correct, so "C++
owns the moves" is a home-when-shipped claim, not a capability claim.
*Verified:* cloned the repository and read `R/update_forest.R`,
`R/tree_proposals.R`, `src/` listing, DESCRIPTION, README; arXiv 2311.12016
abstract fetched (the repository carries no citation).

**nbbart 0.1.0 / spanbbart / soft_spanbbart** (Englert; GitHub
`jacobenglert/nbbart`). Negative-binomial BART, and a spatial version with CAR
random effects. *Unique:* nothing in the integrand - Polya-Gamma augmentation
plus a Gaussian forest. *Built:* **reuse**, by the same author who abandoned
every engine for clbart, which makes the pair a natural experiment.
`R/nbbart.R:32-38` builds one `dbarts::dbarts(x, z, offset = offset,
weights = omega, resid.prior = fixed(1), sigma = 1)`; the loop draws
`omega ~ PG(y + xi, eta)`, forms `z = (y - xi) / (2 omega)`, then
`setResponse(z)` / `setWeights(omega)` / `run(0L, 1L)` at `:81-83`.
`spanbbart` builds its sampler with NO offset argument (`:69`), residualizing
spatial and fixed effects out of the response instead, and mutates at
`:150-152`. *Counterfactual:* **(A)**, demonstrated - and it would have run on
CRAN 0.9-33 unchanged. Two free observations. (i) `R/soft_spanbbart.R` is the
SAME outer loop written against SoftBart's `MakeForest` (`:78`), with
`G <- forest$do_gibbs_weighted(x2, r, omega, x2, 1)[1,]` at `:153` standing
line for line where `spanbbart.R:150-152` puts `setResponse` / `setWeights` /
`run` - one author ported one model across both embeddable engines, differing
only in the forest call, so the surfaces are substitutable for this model
class. (ii) `nbbart` appears to carry a live instance of r-c-division
defect 5: `offset` is passed at `:32` and `:87` computes `eta <- offset + G`
from `G <- samples$train` (`:84`), while `run()$train` already carries the
offset - the CRAN 0.9-33 gaussian branch adds it at
`src/dbarts/bartFit.cpp:2449` (the binary branch does the same at `:2433`).
Inert at the default `offset = rep(0, n)`; double-counts the moment a user
supplies the log-exposure a count model exists for. *Verified:* cloned the
repository and read all three files, plus both CRAN dbarts source lines.

**TBAFTcure 0.1.0** (Sun and Song; GitHub `roxiesun/TBAFTcure`). Tree-based
Bayesian AFT cure model for heterogeneous treatment effects: an incidence
component fit on the whole sample, a latency component on the uncured
subgroup, membership latent and redrawn every sweep. *Unique:* per-component
row sets that differ and move. *Built:* a standalone C++ implementation whose
README states the tree code is "constructed following the framework of R
package bcf" - a reimplementation, not a fork (own `tree.cpp`, `bd.cpp`,
`changerule.cpp`, `swaprule.cpp`, `TBAFTcure.cpp`). *Counterfactual:* **(B).**
Incidence is a probit over all rows and latency an AFT over the uncured, both
shipped families here; what is missing is a per-sweep row channel for a
non-gaussian family - `setWeights` is refused for AFT, and the per-forest
weight channel (`dbarts:::bartcoreSetForestWeights`) is internal, absent from
`dbarts.h` and from the R5 class. Unlocking items: the active-rows mask (AFT
is its S2 scope) or a public per-forest weight. *Verified:* cloned the
repository - README.md, DESCRIPTION, `src/` listing; the model shape agrees
with `docs/design/model-space-survey.md:164-183`, which verified the paper and
the per-sweep `glabel` draw first-hand.

**BPCF** (Kim et al.; GitHub `lit777/BPCF`, C++). Principal stratification
with Bayesian causal forests. *Unique:* imputed potential intermediates enter
the outcome ensemble as SPLITTING COVARIATES, not as a row filter. *Built:*
standalone C++. *Counterfactual:* **(A)** - a per-sweep predictor column is
`setPredictor(x, column = j)`, which ships, with transactional roll-back and a
per-observation install mask. *Verified:* the repository exists and is C++
(GitHub API); the "every component sees all n rows every sweep" reading is
carried from `model-space-survey.md:190-196`, which verified paper and code -
**UNVERIFIED here first-hand**.

**embarcadero** (Carlson) - a dbarts-based species-distribution package the
2026-07 roadmap survey listed as a CRAN reverse dependency. It is NOT on CRAN
today: the `src/contrib` index has no `embarcadero` tarball and the dbarts
CRAN page's reverse-dependency lists do not name it. Treat that reference as
stale. *Verified:* CRAN index grep, dbarts CRAN page.

## Synthesis: forced versus chosen

Thirty-five entries: 15 R engines, 9 R consumers and layers, 5 Python, 6
research artifacts and corrections. Twenty-five carry a primary grade; ten do
not (BayesTree as origin, dbarts itself, the seven consumers that never drive
a sampler between sweeps, the embarcadero correction).

- **(A) 8**: BART's survival/recurrent/competing-risk suites, bartMachine's
  model, VCBART (independent-error arm), riAFTBART, EBcoBART, BART-Survival,
  nbbart, BPCF.
- **(B) 3**: bartcs, princeBART, TBAFTcure.
- **(C) 7**: SoftBart, flexBART, nftbart, sparseVCBART, ridgeBART,
  ShrinkageTrees, monbart.
- **(D) 7**: bcf (the amplitude move only - its additive structure is (A)),
  stochtree, clbart, pymc-bart, bartz, stochtree-Python, xbart. Named
  secondary (D)s without their own entry: BART's multinomial softmax,
  bartMachine's platform.

Five things the entry-by-entry pass shows that a headline count does not.

1. **Reimplementation was usually FORCED, but not always by the model.** The
   (C) and (D) entries split on two unrelated causes. One is the integrand:
   soft routing, graph-partition splits, horseshoe leaves, random-feature
   leaves, non-conjugate stratum likelihoods - genuine clause-2 and clause-D
   territory where no surface would have helped anyone. The other is the
   PLATFORM: JVM, JAX, Rust, Python, a C++ library with two bindings. Four of
   the seven (D)s are platform, not model, and for those authors an engine's
   composition surface was simply irrelevant.
2. **Where reuse was CHOSEN it worked, and one author chose both ways.**
   Englert reused dbarts for negative-binomial and spatial negative-binomial
   BART - and ported the identical loop to SoftBart, showing the barrier was
   the model, not the engine - then wrote pure-R tree moves the moment the
   leaf stopped being conjugate. That is the cleanest evidence available that
   the boundary people actually navigate is conjugacy of the leaf, not
   familiarity with an API.
3. **The three (B) entries all want one missing thing: a row channel for a
   non-gaussian family.** princeBART (probit strata), TBAFTcure (an AFT cure
   subgroup) and, obliquely, bartcs (a coupling re-installed per sweep) are
   blocked not by the integrand but by dbarts's per-observation channels
   stopping at the gaussian family. Gaussian row subsetting already ships
   through zero weights and nobody needed telling; the gap is exactly one
   family boundary wide.
4. **The composition failure mode is real, silent and unattended.** Of the
   artifacts here that drive a dbarts sampler iteratively, riAFTBART uses
   posterior MEANS as draws and nbbart appears to double-count its offset.
   Neither is a capability gap: one is an SBC-detectable error in the outer
   loop, the other a documentation-and-accessor problem in the fit channel.
   The landscape's composition failures are failures of the CONTRACT, not of
   the surface.
5. **Peer PROSE lags dbarts's capability; the peer TABLE is fairer than the
   prose.** SoftBart's paper claims to be the only embeddable BART, and
   stochtree's calls dbarts "a limited interface for interchanging a forest
   MCMC step" while giving the "C/C++ API" and "Random effects" rows one check
   each - but the same Table 1 DOES check dbarts for "Custom sampler
   interface", one of only four packages to get it. The two unchecked rows are
   refutable against the contemporaneous CRAN releases the papers cite
   (SoftBart 0.9-22, stochtree 0.9-28): that header layout predates both -
   0.9-33 still installs `inst/include/dbarts/`, `R_C_interface.hpp` included,
   which is exactly how stan4bart links to it - and `rbart_vi` ships. So the
   honest reading is narrow: the prose is grudging, the table is mostly fair,
   and two rows are wrong. Meanwhile SoftBart already exposes the named
   leaf-prior scale this program's calibration wall is about, and stochtree
   already exposes the residual as a mutable object with `add_vector` /
   `subtract_vector` - on the two surfaces treated here as distinctive, the
   competition is closer than either the perception or the tables suggest.

What the landscape does NOT show, and must not be read to show: that any
particular dbarts work item is worth more than another. Seven of the ten
most-downloaded artifacts in the traction list are model-suite packages rather
than composition surfaces (two of the other three are in-house dbarts
consumers), and traction does not track model completeness - dbarts leads at
6043/mo while shipping neither soft trees nor serialization nor a Python face.
What either observation implies for scheduling belongs to the ranking pass.

## Appendix: the long tail

One line each, from the CRAN index, package DESCRIPTION or the GitHub API
today.

- **dlmtree 1.2.0** (Mork, Im; 366/mo) - Bayesian treed distributed lag
  models; own C++ (Rcpp/RcppArmadillo/RcppEigen); tree ensembles over exposure
  lag structure, adjacent to BART rather than an instance of it.
- **bark 1.0.5** (Ouyang, Clyde) - Bayesian Additive Regression KERNELS; a
  different basis, listed only because full-text search returns it.
- **bartMachineJARs 1.2.2** - bartMachine's shipped JVM dependency.
- **BayesTreePrior 1.0.1** - simulation from the CGM tree prior; no sampler.
- **tidytreatment, voi, countSTAR, MatchIt, WeightIt, marginaleffects,
  insight, bundle, butcher, tmle, mcmcsae, adrftools, lorax, tidyAML,
  twoStageDesignTMLE, glossa, funcml, nlfh** - CRAN packages importing or
  suggesting dbarts without touching its sampler surface (from the dbarts CRAN
  page's reverse-dependency lists, 2026-08-12).
- **BCFM 1.0.0** (Bayesian Clustering Factor Models) and **ssBartik 0.1.1**
  (shift-share Bartik instruments) - name collisions, not BART; recorded so a
  later sweep need not re-investigate them.
- **CDCgov/BART-Survival** - entered above; the JOSS paper referenced by the
  2026-07 roadmap survey was **not re-verified** here - **UNVERIFIED**.
- **theodds/BART4RS, nillen0/SurvBART, twj8CDC/BART_SURVIVAL** - further
  survival-BART repositories returned by GitHub search; contents **not read** -
  **UNVERIFIED** beyond existence and description.
- GitHub-only artifacts entered above: **jacobenglert/nbbart**,
  **jacobenglert/clbart**, **jaredsmurray/monbart**, **roxiesun/TBAFTcure**,
  **AlkemaLab/prince_BART**, **lit777/BPCF**, **JingyuHe/XBART**.
