# How sampler comparisons measure the time to reach the stationary distribution

Status: SURVEY, 2026-11-07. A research survey for pricing the logistic collapsed scale move (dec-B622); it
recommends nothing on its own. The measure pricing adopts is in the plan's
[The convergence measure](../plans/k-move-pricing.md#the-convergence-measure).

## The question

A Markov chain started anywhere but its stationary distribution draws, at sweep t, from a distribution that is
not yet the posterior: the marginal of its state after t sweeps from the start, here written pi_t. Two different
questions are asked of a sampler. The first is how many sweeps it takes before pi_t is close enough to the
posterior, from a given start: the convergence or burn-in question. The second is how many sweeps, once there,
buy one independent draw: efficiency, read as effective draws or effective draws per second. A sampler can be
good at one and poor at the other, and two samplers compared on efficiency alone can rank differently on
convergence. This survey covers how published comparisons read the first, and what they need to read it.

The answer needs three things: a start (one point, or a distribution of starting points), a distance between
pi_t and the posterior, read on some set of quantities, and a threshold on that distance. Each method below
supplies them differently, and each fails in its own way. The reviews agree that no method can certify that a
finite run is representative: Cowles and Carlin (1996) compared 13 diagnostics on two simple models and found
that every one could miss the failure it was designed to catch, and recommend several in combination, run on a
few parallel chains, with autocorrelations watched and parameterizations or algorithms changed where they mix
poorly. Roy (2020) reaches the same conclusion a generation later, adding the rigorous methods developed since.

## The methods

### Bounds from theory

**What it measures.** A bound on the total variation distance between pi_t and the posterior, as a function of t
and the start, or a rate at which that distance shrinks. Drift and minorization conditions give explicit bounds
(Rosenthal's technique, reviewed by Roy 2020); Jones and Hobert (2001), as Roy reports them, define an honest
burn-in as the smallest t at which the total variation distance falls below 0.01. Rates come from the spectral
gap, the conductance, or, for data augmentation, the fraction of missing information and the largest
autocorrelation over linear functions (Meng and van Dyk 1999; van Dyk and Meng 2001). For samplers over
combinatorial structures, canonical paths bound the spectral gap (Kim and Rockova 2025).

**What it needs.** Analysis of the specific sampler on the specific model; no runs.

**Cost.** Mathematical work, not compute. The bounds hold for every dataset the hypotheses cover.

**Where it fails.** The constants are usually unknown or far too loose to say how many sweeps a given dataset
needs, so the results order samplers and describe how the time grows with the sample size or the dimension,
rather than giving a burn-in. Most results cover simple models, and none covers a sum of trees under a logistic
likelihood with a drawn leaf scale.

**Who used it on which samplers.**
- Data augmentation, with working parameters and parameter expansion: Meng and van Dyk (1999) is a theoretical
  study of conditional and marginal augmentation, with the univariate t model as its running example, judged by
  the rate of convergence. Liu and Wu (1999) and Liu and Sabatti (2000) were not reached; van Dyk and Meng (2001)
  describe their results as establishing the validity of parameter-expanded and group-orbit moves and the choice of
  working prior that gives the fastest rate.
- Polya-Gamma augmentation for logistic regression: Choi and Hobert (2013), not reached, are cited by later work
  as proving the sampler of Polson, Scott and Windle (2013) uniformly ergodic under a normal prior on the
  coefficients, through a minorization condition. That is a qualitative guarantee, not a burn-in.
- Imbalanced categorical data: Johndrow, Smith, Pillai and Dunson (2019) show that in the limit of a fixed number of
  successes and a growing sample, the spectral gap of data augmentation, Polya-Gamma for logit and truncated normal
  for probit, shrinks polynomially in the sample size, while adaptive Metropolis's shrinks only logarithmically;
  rare outcomes at large n are where augmentation samplers slow down.
- Tree samplers: Ronen, Saarinen, Tan, Duncan and Yu (2022) prove a lower bound on the mixing time of a simplified
  BART (one tree, grow and prune) that grows with the number of rows; Kim and Rockova (2025) give upper bounds for
  Bayesian CART under connectivity conditions on the signal, show grow and prune cannot reach deep isolated signals
  in polynomial time, and show that proposals adding and removing whole twigs mix polynomially.

### Coupled chains

**What it measures.** An upper bound, estimated by simulation, on the total variation or Wasserstein distance
between pi_t and the posterior, for every t, from a stated distribution of starts. Two chains with the same
kernel and the same start distribution are run with one ahead of the other by a lag L, their randomness coupled
so that they meet exactly and stay together; the distribution of the meeting time bounds the distance (Biswas,
Jacob and Vanetti 2019, generalizing the lag of one in Jacob, O'Leary and Atchade 2020). Jacob, O'Leary and
Atchade use the meeting times to choose the burn-in of their unbiased estimators, as a large quantile of the
meeting time. Johnson (1996), as Biswas and coauthors describe it, bounded the distance from many chains coupled
without a lag, given a constant that is rarely known. Mixing times are read as the first t at which the
estimated bound falls below a threshold: 0.25 in Biswas and coauthors' comparison of exact and approximate
Langevin samplers, the convention of the Markov chain literature.

**What it needs.** A coupled version of every update in the sampler, built so the two chains meet exactly in
finite time: maximal couplings of the conditional draws, common random numbers in Metropolis steps, and for large
or structured spaces a coupling designed to contract.

**Cost.** One pair of chains per estimate, run for the lag plus the meeting time, repeated independently some
hundreds to thousands of times; the repeats run in parallel. Biswas and coauthors used up to 10,000 repeats in
one example.

**Where it fails.** The bound is an upper bound and can be loose, more so at a lag of one; a pair can take
impractically long to meet when the coupling does not contract, as in high dimensions or on discrete
structures; and with a multimodal target a few hundred repeats can all meet in the mode nearest the start and
mislead (Biswas and coauthors show this with 1,000 against 10,000 repeats). Each sampler needs its own coupling,
and a new coupling for a complicated sampler is research: Biswas, Bhattacharya, Jacob and Johndrow (2022) needed a
two-scale coupling for shrinkage-prior regression, and Kelly, Ryder and Clarte (2023) a purpose-built contractive
coupling for phylogenetic trees.

**Who used it on which samplers.**
- Biswas, Jacob and Vanetti (2019): the Polya-Gamma Gibbs sampler against Hamiltonian Monte Carlo at four tuning
  settings on a Bayesian logistic regression (German credit, 1,000 rows, 49 coefficients), all started from the
  normal prior on the coefficients, the bound plotted against iteration for each sampler on one figure; also
  single-site Gibbs against parallel tempering on an Ising model, and exact against approximate Langevin samplers
  as the dimension grows.
- Jacob, O'Leary and Atchade (2020): a Gibbs sampler for nuclear pump failure counts, an Ising model near its
  critical temperature, high-dimensional variable selection and the cut distribution.
- Biswas, Bhattacharya, Jacob and Johndrow (2022): blocked Gibbs samplers for regression under the horseshoe and
  half-t priors, finding fewer than 1,000 iterations enough at 100,000 covariates and that the prior's degrees of
  freedom change how fast the chains meet.
- Kelly, Ryder and Clarte (2023): Metropolis-Hastings samplers over phylogenetic trees with branch lengths,
  scalar parameters and latent variables, on trees with up to 200 leaves.

No coupling of a BART sampler has been published that this survey found.

### Ensembles read against a reference

**What it measures.** The distance between pi_t and the posterior on chosen quantities, read directly: many
independent chains are started from the start distribution, their states at sweep t form a sample from pi_t, and
that sample is compared with a reference sample from the posterior, sweep by sweep. Readings in use are
histograms of the t-th marginal set against the target (Biswas and coauthors, as the check on their bounds), the
squared error of an estimate against the target as the run lengthens (Margossian, Hoffman, Sountsov,
Riou-Durand, Vehtari and Gelman 2024), the 1-Wasserstein distance between a method's draws and reference draws
from long thinned runs (Zhang, Carpenter, Gelman and Vehtari 2022, who include short chains of Stan's dynamic
Hamiltonian sampler), and,
in phylogenetics as Kelly and coauthors review it, the largest discrepancy in a clade's posterior support against
extremely long chains, convergence declared below a fixed threshold (Hohna, Defoin-Platel and Drummond 2008), the
root-mean-square deviation of split frequencies against chains of a billion iterations (Whidden and Matsen 2015),
and a Euclidean distance between split distributions below 0.02 against a reference built from multiple long
runs (Meyer 2021).

**What it needs.** A reference that can be trusted, on every dataset read: a very long run, ideally of a sampler
other than the one under test, or several samplers agreeing; and enough chains that the sample from pi_t is
precise.

**Cost.** The number of chains times the length of the stream, on every dataset; the reference once per dataset.

**Where it fails.** Where no reference can be trusted, which is where the sampler is slowest. A reference made by
the sampler under test reads whether it agrees with itself. The distance reads only the quantities chosen, so a
slow direction outside them goes unseen, and an empirical distance from a few chains is noisy.

**Who used it on which samplers.** Margossian and coauthors on Hamiltonian samplers run as many short chains;
Zhang and coauthors on Hamiltonian Monte Carlo and variational approximations over the posteriordb models;
the phylogenetic studies on tree samplers; van Dyk and Meng (2001), who set the histograms of 5,000 draws of two
probit augmentation samplers against the target and find the standard sampler too narrow after that many draws.

### Agreement between chains: R-hat and its successors

**What it measures.** Whether chains started from different points have forgotten their starts, by comparing the
variance between chains with the variance within them. Gelman and Rubin (1992) run several chains of 2n
iterations from a starting distribution overdispersed relative to the target, discard the first n, and read the
potential scale reduction factor, R-hat. Brooks and Gelman (1998) read it as a trajectory: R-hat, the pooled
variance and the within-chain variance computed in batches over the second half of ever longer runs and plotted
against the iteration, convergence read when all three settle. Vehtari, Gelman, Simpson, Carpenter and Burkner
(2021) rank-normalize and fold the draws and split each chain, so the statistic also sees differences in scale and
in the tails, recommend at least four chains, diffuse starting values, R-hat below 1.01 with bulk and tail
effective draws of at least 400, and add local measures of efficiency for quantiles and small intervals. Vats and
Knudson (2021) tie R-hat one to one to the effective number of draws, so a threshold can be set by the precision
wanted. Margossian and coauthors (2024) group many short chains into superchains that share a start and read a
nested R-hat between them, for runs of many chains with few draws each.

**What it needs.** Several chains from starts dispersed relative to the posterior; no reference.

**Cost.** Nothing beyond the chains; it is read on the run a user makes.

**Where it fails.** It reads agreement, not stationarity. Chains that share a start, or share a slow direction,
agree long before they reach the posterior: Brooks and Gelman's Weibull example shows R-hat falling under 1.1
between 800 and 1,300 iterations while every chain had explored only one part of the space, the empirical 90
percent interval after 1,000 iterations covering 79.8 percent of the target, and Cowles and Carlin's bimodal
example shows chains started around one mode passing. Vehtari and coauthors report R-hat dipping under 1.1 well
before convergence in their simulations, and much more rarely under 1.01. And, as Vats and Knudson show, a
threshold on R-hat is also a demand for a number of effective draws, not only for stationarity.

**Who used it on which samplers.** van Dyk and Meng (2001) plot R-hat against the iteration from three starting
values for each of two augmentation samplers on a probit regression of 55 patients with strong covariates and no
separation, beside autocorrelations, time series and lag-one scatterplots; the improved sampler's R-hat settles
near one much sooner. Ronen and coauthors (2022) read Gelman and Rubin's statistic on the held-out root-mean-square
error of BART fits made with dbarts, eight chains of 200 trees after 5,000 burn-in sweeps, as the number of rows
grows, and find it often above 1.1. Phylogenetic software reads it on scalar parameters and on split frequencies.

### Single-chain diagnostics and trace plots

Geweke's comparison of the early and late parts of one chain, Heidelberger and Welch's stationarity test, Raftery
and Lewis's run-length calculation and visual inspection of trace plots read one chain's output (Cowles and Carlin
1996; Roy 2020). They are cheap and cannot see what the chain never visited. Pratola (2016) judged his tree
proposals for Bayesian CART and BART this way: trace plots of the error standard deviation (convergence read by
eye), effective draws and effective draws per second, Geweke's statistic and acceptance rates, by the share of
rotation proposals, and the coverage of 90 percent intervals, for Bayesian regression trees and BART.

### Hitting times and run-length studies of tree samplers

Kim and Rockova (2025) read, as a stand-in for the mixing time, the number of iterations a chain started at the
root takes to first contain the true tree, over 50 chains each of Bayesian CART, its twiggy variant and a
spike-and-slab sampler; on a signal four levels deep none of 50 Bayesian CART chains arrived within a million
iterations. Ronen and coauthors (2022) count how often the root split of Bayesian CART changes along a chain,
including its burn-in: under 0.2 percent of iterations on the full datasets, across 160 chains. He and Hahn (2023) compare BART from its default start (every tree a single leaf,
1,000 burn-in sweeps, 2,500 draws) with 25 chains started at draws of their fast approximation (XBART) and run
for 100 sweeps with no burn-in, by the coverage and length of 95 percent intervals for the mean function, the
root-mean-square error and the running time; the started chains have the better coverage and error, and take
less time even counted conservatively. Both read the
effect of the start on what users estimate, at fixed run lengths, rather than a distance.

## The nearby samplers, by what was read

| study | sampler and model | convergence from a start | efficiency once there |
|---|---|---|---|
| Meng and van Dyk 1999 | data augmentation with working parameters; t model | theoretical rate | |
| van Dyk and Meng 2001 | data augmentation; t, probit, mixed models | R-hat against iteration from three starts; histograms after 5,000 draws against the target | lag-one and lag-two autocorrelation; lag at which it falls under 0.05 |
| Liu and Wu 1999; Liu and Sabatti 2000 | parameter expansion; moves along group orbits | not reached | not reached |
| Polson, Scott and Windle 2013 | Polya-Gamma Gibbs; logit and negative binomial benchmarks | none in the main text | effective draws and effective draws per second after burn-in, median over coefficients |
| Choi and Hobert 2013 | Polya-Gamma Gibbs; logistic regression | uniform ergodicity (not reached) | |
| Johndrow, Smith, Pillai and Dunson 2019 | Polya-Gamma and truncated normal augmentation; imbalanced data | spectral gap as n grows with successes fixed | autocorrelations to lag 100 after 30,000 burn-in; effective draws per second |
| Biswas, Jacob and Vanetti 2019 | Polya-Gamma Gibbs against Hamiltonian Monte Carlo; logistic regression | coupling bound on total variation against iteration, start from the prior | |
| Ronen, Saarinen, Tan, Duncan and Yu 2022 | BART and one-tree BART (dbarts); Bayesian CART | mixing-time lower bound; Gelman and Rubin's statistic on held-out error across 8 chains as n grows; root-split changes | |
| Pratola 2016 | Bayesian regression trees and BART with rotate and perturb proposals | trace plots of the error standard deviation; Geweke | effective draws and effective draws per second; acceptance rates |
| Kim and Rockova 2025 | Bayesian CART, its twiggy variant, spike-and-slab | mixing-time bounds; hitting time of the true tree from the root over 50 chains | |
| He and Hahn 2023 | BART from its default start against chains started at XBART draws | coverage, interval length and error at fixed run lengths | running time |

Under data augmentation and Polya-Gamma samplers, convergence from a start has mostly been argued from theory and
from R-hat and autocorrelation; the one published start-based distance for the Polya-Gamma logistic sampler is
Biswas and coauthors' coupling bound, on a logistic regression with a linear predictor (49 coefficients, the German
credit data), not on trees. For BART, the published readings are R-hat across a few
chains, trace plots, hitting times and coverage at fixed run lengths from different starts.

## What carries over to a BART sampler under logistic

- **Couplings** would give the reading closest to the question, a bound on the distance from the default start at
  every sweep, but none exists for BART. Two chains of 75 trees would have to meet exactly in every tree, leaf, k
  and latent, and the published tree-space coupling needed a construction of its own. A coupling of this sampler,
  with its tree moves, Polya-Gamma draws and k's draws, would be new work with its own proof of validity.
- **A reference** read sweep by sweep is the standard where one can be trusted. In the regimes where the shipped
  sampler is slowest, long runs of the shipped sampler and of the move mostly disagree, so the reference is missing
  where it matters most (the k study's long runs that failed acceptance were concentrated in logistic fits with a
  small true k).
- **R-hat** costs nothing and is what a user reads, but chains that start together and share a slow direction, as
  k's slow drift from its starting value is shared, agree before they arrive.
- **A simulated truth.** When the truth is drawn from the model's prior and the data from the truth, the truth given
  the data is itself one exact draw from the posterior, independent of every chain run on those data. Calibration
  checks of samplers rest on this at stationarity (Cook, Gelman and Rubin 2006; Talts, Betancourt, Simpson,
  Vehtari and Gelman 2018). It supplies a reference draw in every replication, including those where no long run
  can be trusted. Reading it sweep by sweep, as the rank of the truth among independent chains' states at each
  sweep, is not a method this survey found in a published comparison; a design using it is the study's own.
- **Starts.** Published comparisons start from what the user gets or from a distribution overdispersed relative to
  the posterior: Gelman and Rubin and their successors ask for diffuse starts, Biswas and coauthors start the
  logistic comparison from the prior, Kim and Rockova start at the root, and He and Hahn compare the default start
  with an informed one. A draw from the model's prior is the overdispersed start these support.

## Models for the report

pricing's report follows published reports of its kind (dec-B626). Among the sources here, these show what such
a comparison contains:
- one figure per comparison with the distance against iteration for each sampler (Biswas, Jacob and Vanetti 2019,
  their figure 4; Biswas, Bhattacharya, Jacob and Johndrow 2022);
- R-hat, the pooled variance and the within-chain variance as trajectories (Brooks and Gelman 1998, figure 1), and
  R-hat against iteration beside autocorrelation, time series and histograms set against the target, for two
  samplers from the same starts (van Dyk and Meng 2001, figures 3 to 6);
- squared error of an estimate against the target as the run lengthens, with the precision of 100 independent
  draws marked, beside R-hat over the same lengths (Margossian and coauthors 2024, figure 1);
- tables of effective draws and effective draws per second by dataset and sampler (Polson, Scott and Windle
  2013; Pratola 2016, table 2);
- rank-normalized R-hat with bulk and tail effective draws and local efficiency plots (Vehtari and coauthors 2021);
- coverage, interval length, error and running time of each start or sampler at fixed run lengths (He and Hahn
  2023, table 4).

## Sources

Read for this survey (copies kept with the project's papers): Biswas, Bhattacharya, Jacob and Johndrow (2022),
Journal of the Royal Statistical Society B; Biswas, Jacob and Vanetti (2019), Advances in Neural Information
Processing Systems 32; Brooks and Gelman (1998), Journal of Computational and Graphical Statistics 7; Cowles and
Carlin (1996), Journal of the American Statistical Association 91; Gelman and Rubin (1992), Statistical Science 7
(a scan, read in part); He and Hahn (2023), Journal of the American Statistical Association; Jacob, O'Leary and
Atchade (2020), Journal of the Royal Statistical Society B 82; Johndrow, Smith, Pillai and Dunson (2019), Journal
of the American Statistical Association 114 (the 2017 preprint); Kelly, Ryder and Clarte (2023), Annals of Applied
Statistics 17; Kim and Rockova (2025), Electronic Journal of Statistics 19; Margossian, Hoffman, Sountsov,
Riou-Durand, Vehtari and Gelman (2024), Bayesian Analysis; Meng and van Dyk (1999), Biometrika 86 (a scan, read in
part); Polson, Scott and Windle (2013), Journal of the American Statistical Association 108; Pratola (2016),
Bayesian Analysis 11; Ronen, Saarinen, Tan, Duncan and Yu (2022), preprint; Roy (2020), Annual Review of
Statistics and Its Application 7; van Dyk and Meng (2001), Journal of Computational and Graphical Statistics 10;
Vats and Knudson (2021), Statistical Science 36; Vehtari, Gelman, Simpson, Carpenter and Burkner (2021), Bayesian
Analysis 16; Zhang, Carpenter, Gelman and Vehtari (2022), Journal of Machine Learning Research 23.

Not reached: Liu and Sabatti (2000), Biometrika 87, and Liu and Wu (1999), Journal of the American Statistical
Association 94, both behind a paywall; Choi and Hobert (2013), Electronic Journal of Statistics 7, open access but
refusing automated download. What is said of them here is how the papers above cite them.

Cited only as other papers describe them, not read: Johnson (1996) and Jones and Hobert (2001), as Biswas and
coauthors and Roy describe them; Hohna, Defoin-Platel and Drummond (2008), Whidden and Matsen (2015) and Meyer
(2021), as Kelly and coauthors review them; Cook, Gelman and Rubin (2006) and Talts and coauthors (2018), cited
for the property of simulated truths that pricing's step 1 already uses.
