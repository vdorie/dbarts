# dbarts 1.0-0: what it is, what it breaks, what is left

This memo answers four questions the maintainer confirmed on 2026-09-13:

1. What is dbarts 1.0-0?
2. What breaks for someone using 0.9-34 as they use it today, in dbarts
   and in the sister packages we own?
3. What is missing or wrong that should be fixed before a release
   candidate and before 1.0-0, and does each item go before or after the
   merge to main?
4. Where 1.0-0 and 0.9-34 fit the same model, do they give the same
   posterior?

## What dbarts 1.0-0 is

dbarts 1.0-0 replaces the sampling engine with a rewrite that keeps the
sampler object's interface and extends what it can fit. It adds eight
response families that 0.9-34 did not have: Student-t residuals, logistic,
ordinal, multinomial, negative binomial, log-normal accelerated failure
time, discrete-time hazard, and two-part log-normal. It adds a variance
forest for heteroscedastic residuals under the gaussian and survival
families, a multi-forest family in the shape of Bayesian causal forests,
monotone and interaction constraints, linear and Gaussian-process leaves,
the DART sparsity prior, sparse predictor matrices, in-model handling of
missing predictors, per-draw callbacks from R or C, an active-row mask,
and a flat C header that other packages compile against in place of the
old C++ headers. The cross-validation function is now an R driver over the
sampler, so it reaches the sampler's own priors and settings and
reproduces at any thread count. It repairs defects of 0.9-34: the tree
change move lacked half of its acceptance ratio and biased splits toward
low-cardinality predictors; the residual variance counted zero-weight rows
in its degrees of freedom; the sampled k used the wrong chi shape; a fixed
residual variance was silently ignored; rows with a missing predictor were
silently dropped; a wrong-length weight vector was silently recycled;
seeded multi-chain fits depended on the thread count; and cross-validation
carried one chain across folds, so each fold started from a forest that
had seen its own held-out rows. On a four-core machine with 200 trees and
1000 draws after 500, 1.0-0 fits a continuous response 1.7 times faster
than 0.9-34 at 1000 observations and 2.6 times faster at 10,000; four
chains on four threads and a probit fit are each about 1.7 and 1.6 times
faster. A four-chain fit of 100,000 rows also holds 1.6 GB less in memory.

## What breaks for a 0.9-34 user

Every fit gives different draws. That is silent and unavoidable: the
engine draws a different stream, the change move is repaired, the initial
forest is drawn without empty leaves, and the modern doors' default
proposal mixture sets the swap move to zero, where bartBT keeps
BayesTree's. A script that compares against draws saved under
0.9-34 will not match. Beyond that, in roughly the order of how common the
usage is:

- Calling bart with the old argument names. Loud, once: a call that names
  a BayesTree-style argument is forwarded to bartBT, which carries
  0.9-34's 31 arguments and defaults unchanged, after a once-per-session
  warning. Calling bart positionally with a matrix and a vector is silent
  and changes the answer: bart is now the modern door, formerly bart2, and
  fits four chains of 500 draws with 75 trees where 0.9-34 fit one chain
  of 1000 draws with 200 trees.
- Factor predictors. Silent. Unordered factors split on subsets of levels
  as one column and ordered factors become one ordinal column, where
  0.9-34 expanded both into indicators. The design matrix, the variable
  counts and their names, and the draws all change. bartBT expands as
  before, and the modern doors accept factors = "indicators".
- Missing predictors. Silent for bart, dbarts and the sampler, which keep
  the row and model the missingness, so the number of rows and the length
  of the fitted vector grow. bartBT drops the row, as 0.9-34 did. Test
  data with a missing value in a column that was complete in training is
  now an error where 0.9-34 answered.
- Multi-chain results. Silent shape change. bart now merges the chain and
  sample margins by default, and the merged order is chain-major where it
  was sample-major, so code that reshapes combined draws back into chains
  pairs them with the wrong chain.
- Binary responses. Silent. The leaf scale k is sampled under a proper chi
  prior with scale 2 rather than the improper one; passing k = chi(1.5,
  Inf) restores the old prior. Weights on a probit fit are loud: only 0
  and 1 are accepted, marking rows in and out of the likelihood, and the
  message points weighted counts to the logistic family. 0.9-34 fit
  weighted probit incorrectly.
- Weights with zeros on a gaussian fit. Silent and large: the residual
  scale's posterior no longer counts the zero-weight rows, and on the
  comparison design 0.9-34 reported a residual scale of 0.29 where 1.0-0
  and an unweighted fit on the remaining rows both report about 0.7.
- Residual prior shorthand. A warning, once per session, and the value is
  honoured: bart's sigdf, sigquant and resid.prior, and dbarts's and
  xbart's resid.prior, all still work, spelled now on the family object,
  family = gaussian(sigma = chisq(...)) or fixed(...), removed in dbarts
  1.1-0. Loud where the two spellings collide: giving the prior both ways
  is refused if they disagree, and sigest beside a fixed residual scale is
  now an error, where 0.9-34 accepted it and silently ignored it.
- Cross-validation with xbart. Loud for a three-element burn-in, which is
  refused with an explanation. Silent otherwise: reported losses rise
  because no fold is warm-started on its own rows, and sigma is renamed
  sigest.
- Sampler methods. Loud: a thread count passed to run is ignored with a warning naming setControl. Silent:
  setResponse gains updateScale as its second positional argument, so
  setResponse(y, TRUE) now rescales instead of storing state, with a once-
  per-session warning; the mutators refresh the stored state only when
  asked, so a fit saved after a mutation may not reflect it.
- rngSeed. A warning, and the value is honoured; but a package that
  filters control arguments against the control constructor's formals
  drops it silently and runs unseeded.
- Saved objects. Loud for a saved sampler state, which cannot be restored.
  Silent for a saved data object, which keeps its old design and fits
  differently from a fresh one; an ordered factor in it has no upgrade
  path.
- Removed functions. Loud: rbart_vi points to stan4bart and its plot and
  print methods are gone; a three-or-more-level factor response through
  bartBT is refused with the multinomial route named; the control
  constructor loses its RNG-kind arguments; startThreads and stopThreads
  warn and do nothing.
- Smaller silent changes: NULL components are dropped from the fit object;
  the third positional argument of fitted is the interval level.
- Installation. Loud: R 4.2.0 and a C++20 compiler are required, and any
  package that compiled against the old C++ headers fails to compile.

The sister packages we own each need their ported branch:

- stan4bart. The CRAN release fails to load against 1.0-0, since it
  resolves symbols the flat header no longer exports. The ported branch
  builds its sampler in R, drives it through the header, seeds correctly,
  and gains a slice move that fixes the group-spread mixing the removal of
  rbart_vi was conditioned on.
- bartCause. The CRAN release installs, then errors on every response
  route, because it assigns into the predictor matrix of a data object in
  a way 1.0-0 refuses; its grouped route fails earlier on rbart_vi. The
  ported branch fits grouped models through stan4bart, adds a causal
  forest fit, and matches the chain-major layout.
- treatSens. Not on CRAN. Its current source compiles against deleted
  headers; the ported branch uses the flat header and bartBT for its null
  surface. It still reaches two dbarts internals by namespace.
- bairrtt installs and runs unchanged and silently returns a different
  posterior.

## What is missing or wrong before a release candidate and before 1.0-0

Before the merge to main, because each is shipped surface or a decision
the merge is made on:

- At the freeze: re-verify the four sister packages against the final
  header. They passed on 2026-09-23, so this repeats only if the header or
  the R functions they call change before then. Bump the DESCRIPTION date.
  Then the maintainer's own items: contacting lorax's maintainer about its
  three-level factor response, contacting WeightIt's and MatchIt's
  maintainer about bart2's removal in 1.1-0 and the rngSeed argument their
  filters drop, and closing the setResponse issue.

After the merge, because they cannot run before it or the maintainer
placed them there:

- The scheduled workflows, equivalence, rchk, the reverse-dependency smoke
  test, SBC and valgrind, bind to the default branch and cannot run until
  the merge. rchk and valgrind have each run once by hand, rchk clean, and
  they are the checks CRAN itself runs.
- Placed there by ruling, before the 1.0-0 submission: the engine stops
  calling R for densities, printing and error reporting, and the C
  interface gains an entry that creates a sampler from a plain-C
  specification, so a host without R, such as Python, can build one. The
  sister packages are re-verified after it.
- Deferred by ruling: the fused residual pass loses up to 8 percent below
  its crossover size and ships without a size gate; the rule-Gibbs move
  ships at zero weight and is adopted later; the binary prior is revisited
  after the mixing research, since the sampled k never converges at any
  affordable length even though the fitted probabilities do; and C entries
  to shift a constant between the forest and a host's intercepts, which
  stan4bart's remaining mixing failure needs, land as a minor header
  addition.
- Gaps with no decision yet: xbart reaches only the gaussian and binary
  families; the negative binomial has no continuous dispersion; the
  alternate families have no warm start; the multinomial family has no
  per-observation log-likelihood channel; the negative binomial SBC arm
  passes only with two functionals waived as an identifiability ridge.

## Whether 1.0-0 and 0.9-34 give the same posterior

Where both releases fit the same model they were compared directly. Both
were installed side by side and one script in 0.9-34's vocabulary was run
under each on 26 scenarios covering continuous and probit responses,
offsets, test sets, weights, factor predictors in both encodings, ordered
factors, cut-point settings, tree counts from one to 200, 5000
observations, four chains, an embedded Gibbs loop, predictor swaps and
cross-validation. Each scenario ran 20 seeds a side with every moved
default pinned on both sides, and posterior summaries were compared by a
two-sample z statistic. Twenty-two scenarios agree at the rate the null
predicts. The other four differ in three ways, and each difference was
traced to a defect in 0.9-34 with a control that isolates it: zero
weights, where dropping the rows makes both releases agree with 1.0-0's
weighted answer; the change move under unequal cut counts, where removing
the move makes them agree and 0.9-34's answer drifts with the move's share
while 1.0-0's does not; and cross-validation, where a fold-by-fold loop
through the plain fitting function agrees with 1.0-0 and not with 0.9-34.
Two expected shifts, the empty-leaf initial forest and the chi degrees-of-
freedom relabel, are below what this comparison can resolve, which is
about a third of a posterior standard deviation on a fitted value. The
comparison does not reach the sampler's accessors, state saving or
prediction from a saved sampler.

For everything 0.9-34 cannot fit, the evidence is internal. Every family
with sampling code of its own has an exact check against a closed-form or
brute-force posterior on a small design, and the Student-t, ordinal,
negative binomial, multinomial, survival and both heteroscedastic
compositions have simulation-based calibration arms that rank prior-drawn
truths against posterior draws at 200 replications. The logistic family
and monotone constraints have reference checks but no calibration arm. The
multi-forest family has exact checks and its own bitwise baseline but no
calibration arm, because its amplitude chains decorrelate too slowly for
one. Hazard and two-part fits are compositions with no sampling code and
are checked by reduction to the probit and logistic fits they expand into.
