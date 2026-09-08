# tree-mixing-proposals: how the BART posterior is sticky, and what could move it

Status: COMPLETE (survey with an adjudicated evidence base), 2026-08-09;
**section 12 is an addendum, 2026-08-10**; **section 13 is a measured
move-set A/B, 2026-09-06**; **section 14 measures recovery after a response
swap, 2026-09-06**; **section 15 is a four-lens proposal brainstorm with a
refutation pass, 2026-09-07**; **section 16 is a novelty-gated
first-principles brainstorm with a refutation pass, 2026-09-07**. Section
12 refutes one of this document's
recorded inferences (erratum in sec 5.4), amends sec 3.1, and adjudicates
fourteen new candidates. Section 13 reproduces Tan et al.'s null on this
package's own kernel and, with a one-paragraph addendum in sec 6.1,
supplies the per-move acceptance rates the move census asked for.
Section 14 finds that section 13's null does NOT extend to a moving
response - dropping the change move costs 20 to 25 percent more sweeps to
re-adapt after a large swap, while dropping swap costs nothing. That last
finding was acted on: the shipped default now sets swap to zero, the move
itself staying in the kernel for the single-tree case
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)),
so every "shipped default" in this document describes the mixture in force
when it was written, not the one that ships. Section 15 ranks ten proposal
mechanisms across four lenses, marks each with a refutation verdict, and
recommends none of them; two of its findings landed as sentence corrections
in `benchmark-surfaces.md`. Section 16 runs the first-principles sequel
section 15 asks for, ranks eight more mechanisms under a novelty gate, and
recommends none of them either. Nothing else here is proposed
or scheduled. TODO
`tree-mixing-proposals` (VD 2026-08-09: "I'm interested in the ways in
which the posterior is sticky and if we can come up with some other
proposal mechanisms to help BART explore and mix better"). Likely
post-release; timing free.

This document is the durable record. It was produced by a survey pass, an
adversarial critique that refuted the survey's top recommendation, and an
adjudication pass that re-opened every in-repo claim at `d3cb94b` and
spot-checked every external proof against the primary source - including,
for the tempering literature, the authors' released experiment code. The
working papers are untracked (gitignored):
`memo.md` (the survey), `critique.md` (the review), `synthesis.md` (the
per-finding adjudication, ADOPT/OVERTURN with evidence). Section 4.1 was
added at VD's direction during the adjudication pass and is not in either
earlier paper.

Summary: BART's sampler is sticky in five mechanically distinct ways, and
this document separates the two we have measured in this house (the
ensemble sits at many equally-good tree arrangements and never enumerates
them; tree structure freezes when the noise level is low) from the three
that are, today, arguments from the source code with no identifying
measurement behind them. The survey's leading recommendation - warming the
sampler up during burn-in and cooling it back down - **does not survive**:
not one of the five sources cited for it evaluates the construction
proposed, the released code behind the one empirical result does something
different, and the version proposed is pointed the wrong way for the one
failure this package has actually reproduced. What survives first is much
duller and much cheaper: **dbarts has no move that nudges an existing split
point without also redrawing which variable it splits on**, that move's
Metropolis-Hastings correction is provably 1 under machinery already
shipped, and it is the only candidate whose step size can be tuned.

A second candidate, on a different axis, ranks alongside it: instead of
making a sticky chain move better, **give the forest less to do** - fit the
smooth, separable share of the signal in a parametric block (sampled
jointly by an HMC-type sampler, or held in the leaves) and let BART target
the residual. Its falsifier is the cheapest in the document because every
arm already ships. Its hazard - the two components competing for the same
signal, and a two-block sampler crawling along the resulting ridge - is
established in print three times over and measured in this house at 6x, and
is discussed head-on in section 4.1. What has never been measured anywhere
is the benefit: no published work reports a mixing diagnostic for a
composed model, and stan4bart's own paper poses the question and says
"More research will need to be performed to confirm this."

The recommended next step is two pure measurements, neither of which
changes a draw or needs an engine change: the per-move acceptance rate,
which nobody has ever taken on this sampler, and the composition probe.

---

## 1. The question, and why the answers are hard to get

**Mixing** is how quickly a Markov chain Monte Carlo (MCMC) sampler moves
around the set of models the data support. A **sticky** posterior is one
where the sampler settles into one region and stays: the draws it reports
are correlated with each other, so a thousand draws carry much less
information than a thousand independent ones. BART's sampler visits tree
structures by proposing a small local edit and accepting or rejecting it
with a probability that keeps the long-run distribution correct
(**Metropolis-Hastings**, MH). The **acceptance rate** is the fraction of
proposals accepted; a rate near zero means the chain is frozen.

Three facts make this question unusually hard to answer empirically in
this package, all established by the grow-from-root default study
(`grow-from-root-default.md`, KILLED 2026-08-08):

- **At the shipped default of 75 trees, no structural statistic can detect
  structural mode collapse.** The ensemble self-averages any label like
  "the root splits on x1" at a rate of `1/sqrt(#splits)` with no mixing
  required to produce that average. Structure probes therefore have to run
  at one tree, on a purpose-built scenario.
- **R-hat is not gateable on this package's slow functional.** R-hat is the
  standard cross-chain convergence diagnostic (values near 1 mean the
  chains agree). On held-out prediction error, both a cold and a warm arm
  sat at absolute R-hat 1.28 (n = 2000) and 1.53 (n = 20000) after 1000
  draws, with per-replicate standard deviation of the *difference* 0.198
  and 0.482. A 0.05 margin at 4x the standard error would need roughly 250
  and 1500 replicates. (Those are Stage-0 probe cells; the study's own
  m = 75 cells show cold R-hat of 1.003-1.026, so this is a statement about
  one functional at those probe sizes, not a claim that dbarts is broadly
  non-converged at defaults.)
- **The house design law that follows.** Per-cell kill checks, never a
  pooled aggregate - the study's aggregates all passed while the per-cell
  check killed it. Thresholds frozen against a pilot before any
  confirmatory contrast. A mandatory fresh-seed re-run of any single
  flagged cell. A null control whose failure voids the estimator family.
  Every falsifier below is bound by these.

---

## 2. What the sampler can and cannot do today

`metropolisJumpForTree` ([`metropolisJumpForTree`](../../src/bartcore/moves.hpp)) draws one uniform per tree
per sweep and dispatches to exactly one of three kernels:

```
u < birthOrDeathProbability          -> birthOrDeathMove
u < birthOrDeath + swapProbability   -> swapMove
else                                 -> changeMove
```

Shipped mixture `birth_death = 0.6, swap = 0, change = 0.4, birth = 0.5`
(`defaultProposalProbs` [`defaultProposalProbs`](../../R/model.R), `dbarts()`'s
`proposal.probs` [`dbarts`](../../R/dbarts.R), engine defaults
[`SamplerOptions`](../../src/bartcore/chain.hpp)).
`StepType` is `{birth, death, swap, change}` ([`StepType`](../../src/bartcore/moves.hpp)): the
"four-move set" is three kernels with four labels. The sweep is
Gauss-Seidel over trees; the variance forest runs the identical kernel
([`sweepVarianceForest`](../../src/bartcore/chain.hpp)).

Everything below that discusses a swap move at a positive weight was written
while the shipped default carried `birth_death = 0.5, swap = 0.1, change =
0.4, birth = 0.5`. The move is still in the kernel and still reachable
through `proposal.probs`; only the default moved.

- **Birth** picks a leaf uniformly among those with a usable variable and
  draws the new splitting rule from the prior, so the rule's prior density
  cancels in the acceptance ratio. One split, at the fringe.
- **Death** removes one split, at the fringe, among nodes with no
  grandchildren.
- **Swap** exchanges the rules of a parent and one child. It preserves the
  tree's shape, is symmetric, and carries no proposal correction
  ([`swapMove`](../../src/bartcore/moves.hpp)). It cannot lift a rule more than one level.
  It ships at probability zero, being nearly all no-op at production forest
  sizes, and is the only move that rotates a rule up the tree, which a
  single-tree fit needs.
- **Change** picks uniformly among all non-leaf nodes *including the root*,
  redraws the split variable from the prior, then the cut point uniformly
  over the descendant-valid set. **The entire skeleton below the node is
  held fixed** ([`changeMove`](../../src/bartcore/moves.hpp)) and every observation is rerouted
  through it, with a hard veto if any descendant leaf empties
  ([`resolveVetoRank`](../../src/bartcore/moves.hpp)).

Three consequences.

1. **Only birth and death change the number of splits, and only at the
   fringe.** Pratola states the same for the classic move set: "of all the
   proposal mechanisms that have been developed in the literature, only the
   birth/death move changes dimensionality of the model. Because these
   moves can only alter the very bottom of the tree, it is very unlikely in
   a practical amount of time for a regression tree MCMC algorithm to fully
   explore the space of nearly equivalent trees that have high posterior
   probability." [verified: arXiv 1312.1895 sec 4]
2. **No kernel re-derives the rules below a changed ancestor.** Change
   installs a new rule and reroutes members through an unchanged skeleton,
   so at a deep node it is near-certain to ruin a descendant leaf.
3. **dbarts has no same-variable-only cut move.** `changeMove` always
   redraws the variable, so a pure cut-point move happens only when the
   redraw lands back on the incumbent variable: probability `~1/p_avail`
   per change proposal (about 2% at p = 50 under the default
   `split.probs = NULL`, uniform over the available variables,
   [`CGMTreePrior::drawSplitVariable`](../../src/bartcore/model.hpp)). The machinery for such a
   move already exists and is already gated (section 5.1).

Adjacent machinery already landed, relevant to cost:
`scanOrdinalCuts` [`scanOrdinalCuts`](../../src/bartcore/scan.hpp) (leaf-templated full-cut scan);
`growTreeFromRoot` [`growTreeFromRoot`](../../src/bartcore/grow.hpp); `growForestFromRoot`
[`growForestFromRoot`](../../src/bartcore/chain.hpp) (opt-in
`n.grow.sweeps`, init only, with a reset/regrow/rebuild/redraw loop in the
same function); `SubtreeSnapshot` [`SubtreeSnapshot`](../../src/bartcore/tree.hpp) (restores node
*contents* for a fixed set of node ids - it cannot undo a shape change);
`collapseEmptyNodes` / `collapseSubtreeToLeaf`
[`collapseEmptyNodes`](../../src/bartcore/tree.hpp), [`collapseSubtreeToLeaf`](../../src/bartcore/tree.hpp).

---

## 3. Five ways the posterior is sticky

For each: the mechanism, what it corrupts, and - kept strictly separate -
what is **measured** versus what is **inferred from the code**.

### 3.1 Many tree arrangements, one fitted function (ESTABLISHED)

**Mechanism.** The map from tree ensembles to fitted functions is massively
many-to-one. The same partition of the data arises from several split
orders; the same ensemble fit arises from permuting tree labels, or from
splitting one main effect across two trees instead of one. Moving between
two such arrangements under fringe-only moves means pruning one back to a
stump and regrowing the other - crossing a valley whose two ends are
equally probable and whose middle is not.

**Amendment (2026-08-10, forest-specialization addendum).** The second
clause above - splitting one main effect across two trees instead of one -
also has a *fixed-structure* form: how the fitted signal is apportioned
among the trees' leaf values, a slow direction with a derived timescale of
`n_leaf nodeScale^2 / (m k^2 s^2)` and no label-invariant functional
reading it. That is not a sixth mode; it belongs here, and section 12.2
(B5) carries the derivation, the corrected `k^2` scaling, and the reason
its marginal autocorrelation time is the wrong thing to measure.

**Corrupts.** Everything built on tree structure - variable inclusion
proportions (`varcount`), interaction reporting, `plotTree`, DART's
split-count feedback loop ([`Chain::run`](../../src/bartcore/chain.hpp) recomputes split counts
from the current forest each sweep, so a locked structure feeds a locked
prior), and structural readouts in bartCause and treeSens. It does not
obviously corrupt the fitted function itself.

**Measured, in this house.** At one tree on an XOR scenario
(`y = 4*XOR(x1>.5, x2>.5) + eps`, n = 5000, 24 replicates x 8 chains, 2000
draws), the between-chain standard deviation of the fraction of draws with
the root on x1 was **0.3619** against a mixing null near 0.05, with pooled
p1 = 0.3707 where the correct answer is 0.5 by symmetry
(`grow-from-root-default.md` sec 4.8). The chains move their roots freely
(cold mean 123.5 switches per chain) and *still* disagree that much. That
is representation multimodality, cleanly.

**Independent.** Deshpande makes both halves explicit: "despite BART's
failure to mix over tree space, we often obtain accurate point estimates
and reasonably well-calibrated uncertainty intervals"; but "the fact that E
is not identified complicates the use of these heuristics" for variable
importance and interaction detection. "It is, in our opinion, essentially
hopeless to expect our local grow/prune transition kernel to navigate
efficiently between such representations."
[verified: arXiv 2211.04459 appendix B3]

### 3.2 Structure freezes when the noise level is low (ESTABLISHED)

**Mechanism.** Every structural acceptance carries `exp(dLogL)` where the
integrated log-likelihood difference scales like
`explainedSumOfSquares / residualVariance` - a fact documented in the engine
itself ([`logIntegratedLikelihood`](../../src/bartcore/model.hpp)). As the residual standard deviation sigma falls,
or as a multiplicative forest weight rises, that exponent's magnitude grows
and every proposal that is not an improvement is rejected outright. Leaf
values converge in a handful of sweeps; what remains is a partition-shape
misfit that sigma has to absorb, and absorbing it keeps sigma high, which
is the only thing keeping acceptance non-zero.

**Corrupts.** Sigma's own posterior, interval coverage for the fitted
function, and every structural readout. Point estimates survive.

**Measured, in this house** (`docs/plans/archive/bcf-sigma-residual.md`, measured
at `bartcore 6944811`). In the causal-forest sampler's prior tail, where a
scale parameter `a` multiplies one forest's contribution and the engine
hands that forest weight `w_i * a^2` ([`formForestResponse`](../../src/bartcore/combiner.hpp), [`forestMultiplier`](../../src/bartcore/combiner.hpp)):

- At `a0 = 40` and `100`: "sigma plateaus ~5x high with NO decay through
  40k sweeps - frozen structure."
- The decisive causal test: injecting a large `a` (raising the forest's
  effective signal-to-noise ratio) made things **worse** mid-burn (bias
  1.66 at 18k sweeps versus 1.52 cold), recorded as "direct evidence the
  bottleneck is tree structure, not scale."
- Burn curve on 10 strong replicates, bias = mean window sigma / true
  sigma: 2k sweeps 2.75, 9k 2.20, 18k 1.52, 36k 1.26, 72k 1.07. Longer burn
  was the only lever that worked; a warm start was not (1.21 at 72k against
  cold 1.07).

Scope, stated honestly: the item is RESOLVED for acceptance by pinning
absolute burn sweeps, and the doc's own reframe says the strong-`|a0|`
regime "is real only when the build scale is stale relative to a swapped-in
response (`setResponse(updateScale=FALSE)` inside a larger Gibbs sampler)
or in SBC's own prior tail". That is not obscure for this package - a
`dbartsSampler` inside a larger Gibbs loop is dbarts' distinguishing use
case - but it is not ordinary single-fit `bart()` either.

**Independent.** Pratola's Friedman example: n = 5000, m = 200,
sigma^2 = 0.1, birth/death acceptance ~4%, empirical coverage of the 90%
interval 53%: "the tree structure became stuck in a local mode with, for
all practical considerations, zero chance of moving to a different area of
tree-space that could give an equally good fit." At sigma^2 = 1 the same
setup mixed "reasonably well", acceptance ~18%, coverage 81%.
[verified: arXiv 1312.1895 sec 2.2] **Two scope facts that must travel with
this citation**: that baseline sampler is verbatim "the BART MCMC algorithm
(with birth/death proposals only)" - no change, no swap, against dbarts'
0.4 and 0.1 - and the paper states no tree-prior hyperparameters anywhere,
so the match to dbarts is on `(n, m, sigma)` and on nothing else.

### 3.3 A high split cannot be changed once the tree is deep (CODE-DERIVED HYPOTHESIS)

**Mechanism.** Changing the rule at a node of depth `d` requires the
skeleton below it to remain sensible under a completely different routing
of that node's members. `changeMove` keeps that skeleton fixed by
construction and vetoes any emptied descendant leaf unconditionally, so
acceptance should fall combinatorially with the depth of the subtree below.
Swap can only lift a rule one level. Death cannot reach a node with
grandchildren. So a tree's root variable should be pinned until the tree
collapses back to a stump.

**This is an argument from the source, and this pass demotes it from
"measured" to "hypothesized".** The datum previously cited for it - 147 of
192 warm chains recording zero root switches on duplicate columns
`x1 == x2`, against 10 of 192 cold with a cold mean of 3.8 switches
(`grow-from-root-default.md` sec 4.8) - **cannot identify the mechanism,
because the study recorded switches and never recorded proposals.** The
proposal rate for that event is computable and small: change kernel (0.4)
times the root drawn uniformly from the non-leaf nodes (`1/|notBottom|`)
times the partner column drawn by `drawSplitVariable` (`~1/p`). At p = 10
that reproduces the cold arm's 3.8 switches in 2000 draws at
`|notBottom| ~ 20`, an unremarkable interior-node count for one tree on
n = 5000 - i.e. **the cold datum is equally consistent with acceptance near
1 and with acceptance near 0.** The warm-cold contrast does not rescue it:
deeper warm trees lower the proposal rate and the acceptance rate at once.

The design record itself (`grow-from-root-default.md` sec 4.8) states the
acceptance-collapse reading as established. That reading is plausible and
may well be right; it is not measured, and this document supersedes it on
that point.

**Independent, with the authors' own limit.** Ronen, Saarinen, Tan, Duncan
and Yu ran Bayesian CART **with the full move set** precisely to test
whether change and swap rescue the root: "the root split changes in less
than 0.2% of the samples on average across 160 chains", and "for full
datasets, an overwhelming majority of the root splits occur on the same
feature, and furthermore, this feature is different for different chains".
Four PMLB datasets, California Housing (n = 20640) among them.
[verified: arXiv 2210.09352 sec 4.2] Three scope facts:

- **They ran it with dbarts itself**, named: "We use the dbarts R package
  (Dorie, 2022)", version 0.9-22, `nskip=5000`, `nchain=8`. That raises the
  relevance of the result to this package considerably.
- Their Table 1 fixes Bayesian CART at **one tree**, and section 1.3
  declines the ensemble transfer in their own words: "We did not find
  strong evidence that this bottleneck affects the BART algorithm to the
  same degree."
- 0.9-22 is the classic engine and **predates `change-move-balance.md`**
  (LANDED 2026-07-08), which fixed a since-origin detailed-balance defect
  biasing the change move toward low-cardinality split variables. Their
  change move carried that defect. It biases which variable is chosen
  rather than whether the root moves at all, so the finding stands - but
  the number has never been re-measured on a correct change move, and
  dbarts is the only package that could do it.

**Detection trap, house-established.** Duplicate columns prove nothing
about the sampler: with `x1 == x2` the change move's likelihood ratio and
prior ratio are both exactly 1, so the two "modes" are trivially connected.
Their correct use is as a **null control** - both arms must return pooled
p1 within Monte Carlo error of 1/2 with non-zero switch counts in every
chain (`grow-from-root-default.md` sec 3).

### 3.4 Signals hidden behind a neutral first step (THEORY ONLY, DOES NOT BIND dbarts)

**Mechanism.** Birth proposes one split and scores it against the current
leaf. If the true structure needs two splits before any improvement
appears - XOR is the canonical case, where splitting on x1 alone leaves
both children with the same mean - the first split is a coin flip against
the prior and the second is never reached in a directed way.

**Theory.** Kim and Rockova prove, for one-dimensional dyadic Bayesian CART
with only grow and prune, that it "cannot reach deep isolated signals in
faster than superpolynomial mixing time" (published EJS wording;
the preprint says exponential; both conditional on `L = L_max ~ log(n/2)`).
[verified: EJS 19(2):3041-3067, DOI 10.1214/25-EJS2397, Theorem 5.1;
arXiv 2306.00126] Ronen et al. prove the multi-dimensional single-tree
analogue, also grow/prune only. **Neither binds dbarts' kernel**: single
tree, grow/prune only, and Ronen et al. concede "If either Change or Swap
moves are allowed, the conductance computations would become more
complicated and we may not be able to use the same bottleneck set". Both
papers name a data-fitted initialization as their own remedy - which this
package measured, and killed as a default.

**And the negative datum that has to travel with this.** Tan et al.'s
Experiment 7 found that "restricting the move set [to grow and prune] does
not substantially affect R-hat, coverage, or RMSE"
[verified: arXiv 2406.19958 appendix L.6] - on their own Python
implementation, comparing `{grow .5, prune .5}` against
`{grow .25, prune .25, change .4, swap .1}`. That is a null result about
dbarts' *existing* extra moves. It is the single most directly relevant
disconfirming datum for any new move that is a variation on change, and it
belongs beside those candidates, not quarantined.

### 3.5 Tree size moves by a random walk (CODE-DERIVED HYPOTHESIS)

**Mechanism.** Only birth and death change the number of splits, one node
at a time, at the fringe, against a depth-penalizing prior. Tree size
therefore executes a random walk with O(1) steps; traversing k levels of
size costs O(k^2) sweeps at best.

**No identifying in-repo evidence exists.** The reading previously offered
for this - that the study's autocorrelation-time table shows structural
functionals slow and sigma much faster - does not survive checking. Cold
integrated autocorrelation time on held-out prediction error across the ten
cells is 6.5, 14.6, 15.4, 23.2, 40.2, 45.4, 95.9, 108.0, 177.6, 949.4; the
per-cell ratio of sigma's autocorrelation time to it is 0.89, 0.54, 0.33,
0.76, 0.74, 0.62, 0.40, 0.51, 0.94, 0.51 - median about 0.6, and 888.9 in
absolute terms in the one-tree cell. Sigma is a downstream functional of
the structure and inherits its autocorrelation; a ratio near 0.6 separates
nothing (`grow-from-root-default.md` sec 4.4).

What the table does establish is that **the slow coordinate is slow**:
autocorrelation times in the hundreds on the primary functional at ship
defaults, and a two-long-chain invariance check whose disagreement shrank
only from 0.5693 to 0.3904 posterior standard deviations at 4x length
(ideal would be 0.5x) while the sigma difference **flipped sign** - a
systematic bias cannot flip sign with chain length; a slowly-mixing chain
can (sec 4.9).

**Why it matters here specifically.** Time to re-equilibrate tree size
after the data change is exactly the cost dbarts pays inside a larger Gibbs
sampler (`setData`, `setResponse` between sweeps) - the package's
distinguishing use case.

---

## 4. Candidates that survive, ranked

Ranked by (verified evidence + mechanism) divided by cost, with a bias
toward candidates whose failure is cheap to establish. Cost scale:
XS ~ tens of lines; S ~ a hundred lines plus a surface knob; M ~ a few
hundred lines plus a gate arm; L ~ a new subsystem.

**Two different axes.** Sections 4.2-4.7 are all *tree-space proposals*:
they try to make a sticky chain move better. Section 4.1 is not a proposal
at all - it tries to reduce how much work the sticky part has to do, by
moving part of the signal out of the forest entirely. It ranks first
because its falsifier needs **no engine change and every arm already
ships**, and because a positive result would change what the other
candidates are worth. It is not a substitute for them: it cannot touch the
interaction-discovery mode at all, and its own hazard is measured in this
house at 6x.

### 4.1 Move signal out of the forest, and let BART fit the remainder

**Erratum (2026-08-10, composition mixing probe).** This section's
first-overall rank is WITHDRAWN by measurement. The pre-registered probe
(docs/plans/archive/composition-mixing-probe.md, run to verdicts the day after this
survey landed) returned YELLOW with its registered HARM clause fired and
fresh-seed confirmed, which by the registration kills the blanket
composition recommendation whatever the mixing gates say. What the probe
measured: the representation-transfer leg below is REAL (absorbing the
smooth share robustly shrinks the forest's own job - this section's
mechanism holds); the mixing payoff does NOT reliably follow (no arm
reached the inclusion-efficiency margin at the frozen replication, and the
one fresh-seed pass disagrees with the main block at z 2.5); and the
accuracy guardrail bites where this section did not predict - linear
leaves cost 18% held-out RMSE when there is NOTHING to absorb, while
outer composition buys 15% when there is. What survives is a CONDITIONAL
tool, not a recommendation: outer composition where absorbable structure
is known present. The cross-block ridge is material in every outer arm
(IACT 850-880 of a 2000-sweep window; cor(a_t, b_t) -0.98 to -0.995), the
first such measurement on a composed parametric-plus-BART sampler. The
falsifier this section proposed has therefore RUN; section 7's
recommendation to run it first is discharged.

**What it is.** Put a parametric component - a linear predictor `Z beta`,
random effects, splines, a latent vector - *outside* the forest, sample it
with a gradient-based sampler (Hamiltonian Monte Carlo, NUTS, WALNUTS) that
updates all of its coordinates jointly, and give BART only the residual.
VD's framing: "use a parametric latent vector which can traverse modes more
freely using WALNUTS or some other HMC-type sampler and have BART target
the residual."

**The mechanism, stated so it can be false.** This is **representation
transfer, not mode-hopping.** Hamiltonian Monte Carlo does not traverse
separated modes either - it moves efficiently through a *correlated but
connected* region, which is a different and easier problem. The claim that
has to be true for this to pay is:

> Signal that is multimodal in tree space is unimodal in coefficient space.

A smooth additive surface has exactly one representation as a coefficient
vector and combinatorially many as a deep tree ensemble. That is section
3.1's mode, precisely. Move that share of the signal into a block where it
is one point, and the tree-space multimodality it was generating goes away
rather than being navigated.

**Which stickiness it targets, and which it cannot.**

- **3.1 (many arrangements, one function): directly, and this is the whole
  case.** The equivalent renderings that the chain cannot enumerate are
  renderings *of the smooth share*. Remove the share, remove the
  renderings.
- **3.3 (deep-node lock) and 3.5 (tree size random walk): indirectly, via
  depth.** Both pathologies are functions of how deep the trees are - the
  lock is worse the larger the subtree below a node, and the size walk is
  longer the larger the equilibrium size. A forest with less to explain
  needs less depth for the same fit, so both shrink. This is the most
  defensible of the indirect claims.
- **3.2 (high-signal-to-noise freeze): ambiguous, and possibly backwards.**
  Acceptance turns on the size of a proposal's fit change relative to
  sigma. Absorbing signal shrinks the fit changes the forest needs to make
  (helps) *and* shrinks sigma, because the model as a whole fits better
  (hurts). Worse, the second effect arrives *first*: a parametric block
  converges in a handful of draws while a forest takes thousands, so
  composition rapidly removes the signal that was holding sigma high and
  then asks the forest to discover the remaining interaction structure in a
  low-noise regime - which is precisely the frozen regime. This is not
  speculation about a different model; it is the same mechanism this house
  measured when injecting a large scale into the causal forest "raises mu
  SNR and freezes structure" and made mid-burn bias *worse* (section 3.2).
  **So a real possible outcome of composition is that it improves the fit
  and degrades the forest's structural mixing at the same time**, and the
  falsifier must be able to see that.
- **3.4 (myopia / XOR): not at all.** Interactions are exactly what a
  main-effects parametric block cannot absorb. They stay in tree space and
  the two-splits-before-any-gain problem is untouched. Absorbing them would
  mean specifying them, at which point BART is not doing the work.

**Two variants, and they are genuinely different candidates.**

- **Outer composition**: a separate parametric block alternating with the
  forest. This is stan4bart's shape, of which a random intercept is the
  degenerate case.
- **Inner composition**: the parametric part lives in the *leaves* - dbarts'
  landed linear leaves (`linear-leaves.md`) and GP leaves (`gp-leaves.md`);
  MOTR-BART is the external name for the linear case. **The inner variant
  has no cross-block ridge in the structural move at all**, verified: the
  linear leaf's `logIntegratedLikelihoodForNode` ([`LinearGaussianLeaf::logIntegratedLikelihoodForNode`](../../src/bartcore/model.hpp))
  integrates the leaf coefficients out in closed form (a ridge-regression
  marginal reducing exactly to the constant leaf at q = 0), so the
  acceptance decision never conditions on a realized coefficient. That is a
  strictly better mixing story than the outer variant, and it is already
  shipped. Its honest limit: what it absorbs is *locally* linear (or
  locally smooth) structure per leaf, not a global term - there is no
  coefficient to report, shrink, or put a substantive prior on, and it
  cannot hold random effects. So the two variants are not
  interchangeable, and the depth-reduction claim is the only thing they
  share.

**The hazard, head-on: additive competition makes a ridge, and alternating
updates crawl along ridges.** If both the parametric block and the forest
can absorb the same smooth share, the posterior has a ridge in the
combination, and a two-block Gibbs sampler that updates one conditional on
the other moves along it in small steps. Naive subdivision can therefore
*trade* tree-space multimodality for cross-block correlation that is just
as sticky.

**This house has measured that hazard, in the one composition it already
ships.** `forest-ranef-interweaving.md` (NO-GO, recorded door, 2026-07-20;
its section 9 is a critique-hardened authoritative correction layer):

- The dominant weak-signal bottleneck in the grouped model is *forest-ranef
  confounding* - "the mean forest f and the group intercepts b competing to
  explain group-level structure" - established by measurement, and larger
  than the variance-parameter kernel or the funnel.
- The cleanest number, from an isolated surrogate with a closed-form
  collapse (`benchmarks/R/forest-ranef-collapse-proto.R`), at 3 groups
  where confounding is strong: **alternating conditional blocks give an
  integrated autocorrelation time of 56.1; a joint (collapsed) move gives
  9.3.** A ~6x penalty purely for alternating. At 10 groups, where the
  confounding is mild, the two are level (13.1 versus 11.7).
- On the real engine, removing the forest drops the variance parameter's
  autocorrelation time ~3x at 3 groups, with single-seed excursions to 310
  that the no-forest control never produces (max 42) - corrected by the
  authoritative section 9 to "several-fold, ~3-9x, estimator-unstable,
  worst case unmeasurable by this contrast".

**And the standard remedy does not work here, for a structural reason that
generalizes.** Interweaving / ancillarity-sufficiency (ASIS) and parameter
expansion travel a ridge by alternating two parameterizations. Yu and
Meng's Theorem 1 bounds the interwoven chain's convergence rate by
`R_1,2 sqrt(r_1 r_2)`, where `R_1,2` is the *maximal posterior correlation
between the two augmentation schemes* - so the method's power comes from
having two schemes that are nearly posterior-independent, which is why the
sufficient/ancillary pair is the canonical choice. Two things follow. The
theorem is about two data augmentations *for the same parameter*, linked by
a map, and it does not license treating "parametric block" and "forest
block" as such a pair - that construction would have to be invented. And
operationally, the method needs a *shared scalar coordinate* that both
blocks own. A forest has none: its contribution is emergent from leaf
values organized by leaf, not by group or by covariate, so "shift a
constant out of f and into b" is not expressible without redrawing leaves -
which is the tree sweep. Measured consequence: ASIS on the grouped model
made things **worse** at 3 groups (56.1 -> 114.6), confirmed by the
authoritative critique section as a real mixing effect and not a bug, with
in-repo precedent. Contrast BCF, where the ridge is between two explicit
scalars (`a` and the leaf values) and the rescale move is one line with a
closed-form conditional.

The only remedy that worked in the prototype is a **collapse** - draw the
forest with the parametric block marginalized out. For a per-group
intercept that is priced in-house at ~800-1300 lines on the hot path with
the ungrouped-bitwise guarantee to protect, and it is a declined door.

**Does an HMC block fix the cross-block ridge? No - and this is the
precision point the framing invites.** HMC updates all coordinates of the
parametric block jointly, so it removes ridges *inside* that block (slopes
versus intercepts versus variance components). It does nothing to the ridge
*between* the parametric block and the forest, because that is still an
alternating two-block Gibbs. Verified in stan4bart's own loop
(its `src/init.cpp`, an out-of-repo file this doc cannot bracket-cite):
`dbarts_sampler_run(sampler.bartSampler, 0, 1, ...)`
runs exactly **one** BART sweep, the fit is copied a few lines later into `stanOffset`,
the parametric sampler is given it (`setOffset`), and WALNUTS
draws conditional on it. One-to-one alternation - exactly the configuration
the surrogate measured at 6x. (`forest-ranef-interweaving.md`'s landscape
paragraph describes stan4bart as marginalizing "the ridge in one joint
proposal"; that is accurate for the ridge *within* the parametric block and
overstated for the forest-versus-parametric ridge, which stan4bart
alternates like any two-block Gibbs. Recorded as a correction.)

**Evidence quality: the hazard is well established in print, the benefit is
measured nowhere, and the field openly disagrees about the sign.**

*The hazard is established, independently, three times over.* Hahn,
Carvalho, Puelz and He named it "regularization-induced confounding": a
regularizing prior on one component makes the model prefer to shift signal
into the other, "over-stating the magnitude of the treatment effect
parameter ... while simultaneously attenuating the control variable
coefficients", with a closed-form finite-sample bias. Their remedy is a
reparameterization, and - the one mixing admission in this whole
literature - their appendix adds a step "to improve mixing over the
parameter of interest", noting it "is not possible in the naive
parametrization". BCF diagnoses the same thing for BART specifically, with
a tree-flavoured mechanism: "due to the strong confounding in this example
a single split in Z can stand in for many splits on x1 and x2 that would be
required to approximate mu(x). These simpler structures are favored by the
BART prior, leading to RIC" - and its own fix is again a
reparameterization plus a propensity-score covariate, with mu and tau
stated to "alias one another". For the semiparametric case, CSP-BART is
explicit that a parametric and a BART block sharing covariates hits
"non-identifiability issues", and Bhandari et al. measure the resulting
attenuation directly: "the flexible BART component can absorb variability
that might otherwise be attributed to the linear predictor".

*The remedies in print are reparameterizations and constraints, never
better samplers.* Keep the covariate sets disjoint (Zeldow et al., who
"found that modeling a covariate in both h* and omega sometimes led to bias
and undercoverage"); reparameterize (Hahn et al., BCF); or constrain the
proposal (CSP-BART, below). Orthogonality constraints are proposed and not
built.

*The benefit is measured nowhere.* A dedicated pass searched the full text
of every relevant paper for effective sample size, autocorrelation, R-hat
or Gelman-Rubin: CSP-BART has zero hits in 2592 lines; Zeldow et al. show
trace plots of single coefficients with no accompanying text; MOTR-BART
claims "faster convergence, when we look at the overall log-likelihood"
with no number; the SoftBart vignette's partial linear model asserts "the
chain mixes well" from a trace plot; BCF never uses the words. **No paper
reports a mixing diagnostic for a composed parametric-plus-BART model, let
alone against BART alone.** This is a documented gap, not a search failure.

*And the field disagrees about the sign, in this package's own ecosystem
paper.* stan4bart's section 4.5 takes the opposite view from CSP-BART -
that putting a covariate in both components "may have computational
benefits because it may simplify the nonparametric model", is "an example
of parameter expansion, a technique often employed in Gibbs samplers to
reduce dependence between parameters and increase the efficiency of the
sampler", and that although "neither the parametric nor the nonparametric
components would be directly identifiable ... crucially their sum would
still be". It then flags the claim as unverified in as many words: **"More
research will need to be performed to confirm this."** That was 2022. It is
still unconfirmed, and this probe is what would confirm or refute it.

*One thing that IS quantified, and it is the depth claim.* MOTR-BART -
the inner variant, a linear model in every leaf - reports both halves:
"fewer trees are required to achieve equal or better performance than
BART", and "the trees from MOTR-BART tend to be shallower than those from
BART (10 trees)". Their headline runs use 10 trees against BART's default
200, and on Friedman at n = 1000, p = 50 they estimate 391,193 parameters
against BART-200's 2,371,140 with lower error. So the *transfer* half of
the mechanism - absorbing signal shrinks the forest - has published
support. The *mixing* half does not.

**A remedy in print that is itself a proposal mechanism, and therefore
belongs in this document.** CSP-BART's fix for the shared-covariate
non-identifiability is not a reparameterization but a modified tree kernel:
paired "double-grow" and "double-prune" moves accepted or rejected as a
single Metropolis step, a near-zero-variance leaf prior on the opposing
branch, and rejection of trees whose branch splits only on a shared
variable - so that "the double-grow move ensures that the linear component
estimates only main effects and forces the BART component to work
specifically on interactions and non-linearities". They also flag an
intercept trap worth recording: a leading column of ones in the parametric
part "would conflate the linear component's constant with the constant
node-level mu parameters in the BART component". If dbarts ever grows a
first-class semiparametric surface, this is the design to start from - and
note that a paired-move-accepted-as-one-step kernel is structurally the
same idea as the twig move in section 4.6.

**dbarts fit: this is the cheapest candidate to falsify by a wide margin,
because every arm already exists.**

| arm | what runs it | status |
|---|---|---|
| BART alone | `bart()` | shipped |
| outer composition, HMC parametric block | `stan4bart` (WALNUTS + dbarts exchanging offsets) | exists and runs; 0.0.14 installed here |
| outer composition, arbitrary block, user-driven | `dbartsSampler$setOffset` ([`dbartsSampler$setOffset`](../../R/dbarts.R)), and `dbarts_sampler_setOffset` in the shipped C API ([`dbarts_sampler_setOffset`](../../inst/include/dbarts/dbarts.h)) | shipped, supported |
| inner composition | `node.prior = linear(columns)` / `gp(columns)` ([`dbartsModel`](../../R/model.R)) | shipped |

The probes exist too: the XOR scenario from the grow-from-root study, the
pooled
inclusion-dispersion statistic, and plateau prediction error. So does the
substrate: `inst/common/friedmanData.R` ships a DGP that is exactly one
interaction term (`10 sin(pi x1 x2)`) plus three separable terms
(`20 (x3-0.5)^2 + 10 x4 + 5 x5`), so the absorbable share can be dialed by
changing coefficients while the un-absorbable interaction is held fixed.

**Cost.** The falsifier is XS-to-S: a harness, no engine change, no new
gate, no draw-law change. What a *positive* result costs is the honest
question. If composition helps, the immediate deliverable is guidance plus
documentation of a surface that already exists. If it helps *and* the
cross-block ridge turns out to be the binding constraint, the deliverable
is a collapse - the widest change in the queue, already priced and already
declined once.

**Falsifier sketch, with kill criterion.** Matched-seed paired arms on a
Friedman family whose absorbable share is dialed in three settings
(separable share 0%, 50%, 100% of the signal variance), n = 5000, at ship
defaults. Arms: (A) BART alone; (B) inner composition, linear leaves on the
separable columns; (C) outer composition through `setOffset` with a
conjugate linear block; (D) outer composition with an HMC block
(`stan4bart`). Matched exposure is the design's spine, as everywhere else
here: arms must be compared at matched *BART sweeps*, and arm D's
per-iteration cost must be reported, because the parametric block is not
free. Same threshold discipline as everything else in this document: a
Stage-0 pilot measures the per-replicate standard error of every readout
and the thresholds are frozen against it before any confirmatory contrast
is looked at, with a mandatory fresh-seed re-run of any single flagged
cell.

Primary readout - and it must be a **tree-space** readout, or the study
measures fit quality and calls it mixing: realized tree depth and leaf
count; the pooled between-chain standard deviation of time-averaged
variable inclusion; structural acceptance rate by move type; and root-
switch counts on the one-tree XOR scenario. A **required** companion
readout, because of the freeze mechanism above: recovery of the
un-absorbable interaction term (`10 sin(pi x1 x2)`) and the realized sigma
trajectory, so that "the fit got better and the forest got more frozen" is
visible rather than hidden behind the total error. Secondary:
autocorrelation time on held-out error and on sigma, 90% interval coverage,
plateau error.

**KILL** if, at the 50% absorbable-share setting, no composition arm
reduces mean tree depth *and* improves at least one tree-space mixing
readout beyond 4x the per-replicate standard error over at least 20 matched
pairs - i.e. the transfer either does not shrink the forest's job or does
not translate into better tree-space behaviour when it does. **Also kill
the outer variant specifically** if arm C or D shows a *worse*
autocorrelation time than arm B on the same functional at matched sweeps:
that is the ridge eating the transfer, and it would say the inner variant
is the only form worth pursuing. The 0% and 100% settings are the controls
- at 0% every arm must agree within Monte Carlo error (nothing to absorb),
and at 100% the composition arms must collapse the forest to near-stumps
(the transfer is working) or the harness is wrong.

### 4.2 A same-variable cut move ("perturb") - the first tree-space candidate

**What it is.** A fourth kernel that picks a non-leaf node, **keeps its
splitting variable**, and moves only the cut point. Pratola calls this
"perturb"; his valid-cut interval (his equation 2) accounts for
"constraints from the ancestral ... and descendant ... parts of the tree
about node 5 to avoid any such spurious rejections" - exactly what dbarts'
`findGoodOrdinalRules` already computes.

**Which stickiness it targets.** Freezing at low noise (3.2), and
indirectly the tree-size walk (3.5), because a well-placed cut makes a
subsequent birth acceptable. It does not touch high-node lock (3.3) or
myopia (3.4).

**The mechanism argument, which is the strongest thing here and is new.**
Every move dbarts has makes a *large* change to the partition: birth adds a
split drawn from the prior, change redraws both the variable and the cut,
swap exchanges whole rules. At low noise the acceptance exponent scales
like `dSS / sigma^2`, so a large partition change is a large exponent and
acceptance collapses. **A cut move is the only candidate whose step size is
a free parameter.** Shrink the step and `dSS` shrinks with it, so
acceptance can be held at a workable rate no matter how small sigma is -
the ordinary random-walk Metropolis tuning argument, which none of dbarts'
current structural moves can make. This is a first-principles argument, not
a citation, and it is the reason this candidate is first.

It comes with a fork the argument itself forces: **the property belongs to
a local-window form, not to a uniform redraw over the whole valid
interval.** A uniform redraw makes a large change and its acceptance
collapses like everything else.

**The evidence, read honestly - it is weaker than it first appears.**

- Pratola's ensemble arm reaching 65% acceptance and 92% coverage is
  verbatim "using all MH proposal mechanisms: birth/death proposals, tree
  rotation proposals **and perturb within change-of-variable** proposals",
  and coverage went **96% -> 92%**, i.e. DOWN, when that arm was added to
  the rotation-only arm. Perturb is never isolated there, and it is bundled
  with the correlated-variable proposal that is a separate candidate below.
  [verified: arXiv 1312.1895 sec 5.2, fig 12]
- The one table that does isolate it (Mohammadi, Pratola and Kaptein,
  single tree, n = 300, 100 replications) shows x1/x3 effective sample size
  going 1419/1482 to **13134/13144** when perturb is added on top of
  rotation. But: **x2's effective sample size falls 2899 -> 1306 (-55%)**
  and its per-second rate collapses 3221 -> 759 (-76%); sigma^2's
  per-second rate halves at flat sigma^2 sample size, implying perturb
  costs about 1.9x per iteration in their implementation; and **x1 and x3
  are a deterministic mirror pair** (x1 uniform on (0.1,0.4) then (0.6,0.9)
  by index, x3 exactly reversed) whose splits at 0.5 induce the *identical*
  partition. So the 10x is measured on a purpose-built
  representation-equivalence functional, not on general mixing. In the same
  table the best rotation-containing arm reaches 40041/37925 on that pair,
  so "this puts perturb ahead of rotation" is false.
  [verified: JMLR 21(201) Table 1]
- **Shipped precedent, at a very different dose.** OpenBT runs one
  structural move per tree per sweep and then calls `pertcv`, which loops
  over **every interior node** with no random selection and no early exit
  (OpenBT's `brtmoves.cpp`, an out-of-repo file this doc cannot bracket-cite), sending 90% of those attempts to perturb
  (`pchgv = 0.1` shipped in both the R and Python wrappers). At a typical
  75-tree interior-node count that is roughly 3-6 perturb attempts per tree
  per sweep - 20 to 40 times the throughput the survey proposed testing.
- **Window width is untuned and the two published settings are 8.5x
  apart.** OpenBT's window is 10% of the admissible range centred on the
  current cut (its `brtmoves.cpp`, out-of-repo), with the boundary asymmetry
  corrected by a simple count ratio a few dozen lines later. Pratola's own paper fixes
  it at **85%** ("we have fixed alpha to cover 85% of the range defined in
  the interval (2)"), and his headline numbers were produced at 85%. No
  comparison between the two exists anywhere.
- The nearest disconfirming datum is Tan et al.'s Experiment 7 (section 3.4
  above): restricting away change and swap changed nothing measurable on
  their battery. Perturb is a restriction of the change move.

**Why it is first anyway: the fit to dbarts is exceptional.**

- **The Metropolis-Hastings correction is provably 1** for the uniform
  form, verified by construction in [`changeMove`](../../src/bartcore/moves.hpp): with the new
  variable equal to the old one on an ordinal column, the forward and
  reverse interval and valid-set counts are produced by the *same two
  functions on the same unmodified tree*, and both ignore the node's own
  rule, so `log(rI) - log(fI) + log(fV) - log(rV) = 0` identically - not
  approximately. Both-categorical falls through all three branches of the same
  function and leaves the correction at its `0.0` initializer, which is also
  correct. `change-move-balance.md` records the same: "a same-variable or
  equal-cut ordinal change gives correction 1".
- For the local-window form the correction is `|W(c)| / |W(c')|`, the ratio
  of the two clipped window sizes - still exact, still computed by the same
  two counters, just not identically 1.
- The subtree-prior term below the node must still be computed in general
  (a same-variable cut can flip a descendant's variable availability, which
  changes `growthProbability` and `splitVariableLogProbability` at
  [`CGMTreePrior::growthProbability`](../../src/bartcore/model.hpp), [`CGMTreePrior::splitVariableLogProbability`](../../src/bartcore/model.hpp)), and `changeMove` already computes it.
- Default weight 0 makes the landing bitwise-identical, so the equivalence
  gate stays green by construction and the falsifier can run before
  anything changes for users.

**Cost.** S for the kernel: about 120-160 lines across `moves.hpp`,
`chain.hpp`, `R/model.R`, `R/A_class.R`, `R/dbarts.R`,
`R_interface_bartcore.cpp`, plus ~15 more for the window count ratio.
**The kernel is not the cost.** The correctness gate is (section 5.2), and
so is the study.

### 4.3 Tree rotation

**What it is.** A rotation at an interior node swaps that node's rule with
its parent's, duplicates the parent's other subtree, cuts each copy along
the rotated rule to remove branches that can no longer be reached, and
merges where possible. The partition of the predictor space is preserved:
"the actual decomposition of the X-space has not otherwise changed between
the original and rotated trees". It changes the number of splits at *any*
interior node, which no dbarts move does.
[verified: arXiv 1312.1895 sec 4]

**Which stickiness it targets.** High-node lock (3.3) directly - it is the
only candidate that can change a high node's rule while keeping the
descendants sensible - and representation multimodality (3.1) directly.

**Evidence.** The best-matched published result: Pratola's Friedman example
with rotation at 20% of proposals and birth/death at 80%, acceptance
4% -> 25% and 90%-interval coverage 53% -> 96%. Shipped in OpenBT
(its `brtmoves.cpp`, out-of-repo). Two honest limits: **the number of trees is never
restated** in that section, its figure captions, or the discussion, so
"m = 200 governs it" is inference (it is the only reasonable reading, and
the section says "the same dataset", but it is not text); and the baseline
it improves on is a birth/death-only sampler, which dbarts is not. In the
one controlled single-tree comparison, adding rotation to plain
reversible-jump moves lifted the mirror-pair sample size only 1037 -> 1419,
though it lifted unique trees visited 1.83 -> 3.07.

**dbarts fit - the hard part, and the reason it is not first.**

- *Rollback.* `SubtreeSnapshot` restores node *contents* for a fixed set of
  node ids; it cannot undo a shape change. Rotation needs a full save of
  the subtree at the parent: the node id set, the mask-pool words, and the
  index segment. (OpenBT copies the whole tree per proposal.)
- *Categorical predictors.* Pratola's cut-and-merge assumes a totally
  ordered rule. dbarts' categorical rules are canonical-gauge direction
  masks whose validity invariant is exactly what a rotation would break.
  Ordinal-only in a first version, the same scoping `growTreeFromRoot`
  ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)) used when this was written. That precedent has since
  been spent: the builder scans categoricals too
  (`scanCategoricalPartitions` [`scanCategoricalPartitions`](../../src/bartcore/scan.hpp), the
  winning rule built by `growCategoricalRule` [`growCategoricalRule`](../../src/bartcore/grow.hpp)), so an ordinal-only rotation
  would now be scoping narrower than the builder rather than matching it.
- *Interaction constraints.* A rotation lifts one variable above another -
  the same class of break that `swapMove` already guards with
  `tree.interactionSubtreeIsValid` ([`swapMove`](../../src/bartcore/moves.hpp)). That guard is
  written; note that swap's symmetry additionally relies on a parent's rule
  never equalling a non-selected child's, which holds because
  `splitInterval` gives a child a strictly interior interval. Rotation
  would have to re-establish the analogous property from scratch.
- *The acceptance ratio.* Counting admissible merge arrangements forward
  and reverse, plus the reverse factor for rotations invertible from either
  of two nodes. Pratola enumerates seven merge types and counts recursively.

**Cost.** L. Realistically 400-600 lines plus a merge enumerator plus a new
exact-posterior gate arm; weeks. Highest defect risk of anything here, in a
package whose last change-move defect survived its entire history until an
exact-posterior gate caught it. **Its saving grace is a genuinely cheap
first stage** (section 4.7).

### 4.4 Heated companion chains with private ladders

**What it is.** Metropolis-coupled MCMC (MC^3, also called parallel
tempering). Alongside each chain the user asked for, run a small ladder of
"hotter" companion chains that see a flattened version of the posterior and
therefore move more freely, and periodically propose swapping states
between adjacent rungs. Only the cold chain's draws are kept. The swap is
itself a Metropolis move on the joint space, so the cold chain's kept draws
still target the exact posterior.

**Which stickiness it targets.** Freezing at low noise (3.2) and
representation multimodality (3.1); in principle high-node lock too, since
a hot rung can dismantle and rebuild a tree.

**Evidence.** MC^3 is the reference-quality baseline for tree posteriors in
phylogenetics. Angelopoulos and Cussens applied it to Bayesian
classification and regression trees: four chains, ladder
`beta_i = 1/(1 + 0.2(i-1))`, a swap proposed after every iteration, "Only
trees visited by the cold chain are collected to form the MCMC sample";
"our results show that a clear improvement is achieved using tempering".
The robust part of their result is that the across-seed standard deviation
of accuracy was smaller with tempering in 15 of 16 cases - **over three
seeds** - while mean accuracy was mixed (PIMA fell 76.5 -> 73.4 and
76.9 -> 73.6 in two of three settings). [verified: ICML 2005]

**Why the previously recorded architectural objection does not hold.** The
survey declined this on the grounds that swapping breaks dbarts'
per-chain RNG reproducibility, forces a synchronization barrier across the
thread-parallel chain layout ([`bart2`](../../R/bart.R)), and destroys the
diagnostic value of multiple chains. All three assume swaps happen *among
the user's chains*. Under the private-ladder construction - each cold chain
owns its own rungs and cold chains never exchange with each other - each
cold chain plus its ladder is one deterministic unit driven by that chain's
own Mersenne Twister, the swaps happen inside that unit rather than across
the parallel layout, and the cold chains remain independent starts so
R-hat means exactly what it means today.

**The real cost is compute and tuning.** A ladder of L rungs multiplies
per-chain work by L (forest state also multiplies by L, which is small).
Ladder spacing is the classic MC^3 weakness: too wide and swaps never
accept, too narrow and the hot rung is not hot. And an *unbounded* flatten
of the likelihood alone runs into the same problem as section 5.1 below -
if the depth prior is not flattened too, the hot rung relaxes toward
shallow trees.

**Cost.** M-L for the engine (a swap kernel, per-rung forest state, a
per-unit scheduler), plus a real user-facing decision about what to do with
`n.threads`. Position: this is the *expensive but well-founded* option -
strictly better founded than annealed burn on both validity and evidence.

### 4.5 Informed birth/death over the shared cut scan

**What it is.** One scan per variable over a leaf's members
([`scanOrdinalCuts`](../../src/bartcore/scan.hpp)) yields the collapsed marginal likelihood for *every*
candidate cut at once. Instead of drawing the cut from the prior, propose
over the whole birth/death neighbourhood with weights proportional to
`sqrt(posterior ratio)` (Zanella's locally-balanced construction), and the
acceptance collapses algebraically to `min(1, Z(T)/Z(T'))`, a ratio of two
scan sums. Sketched in-repo at `parallel-bart-frontier.md` sec 3.1.

**Which stickiness it targets.** Freezing at low noise (3.2) - an informed
proposal picks a cut that *is* an improvement, so acceptance survives small
sigma - and the tree-size walk (3.5). It explicitly does not touch
high-node lock or myopia.

**Evidence and its qualifiers.** Zanella characterizes locally-balanced
proposals as optimal *within the class of pointwise-informed proposals*,
asymptotically, under a bounded-degree conditional-independence condition,
and says explicitly that no single balancing function dominates within the
class. On tree posteriors specifically, Zhang, Huelsenbeck and Ronquist
replaced random topology proposals with parsimony-guided ones and found
"single chains using parsimony-guided moves usually converge an order of
magnitude faster"; their own caveat is that "relative performance ...
depends strongly on the data set". Kim and Rockova cap what it can deliver:
informedness alone "does not solve the myopic problem of Bayesian CART" -
footnote-conditioned on the proposal neighbourhood being unchanged, which
is the whole argument. Deshpande warns that a naive informed proposal
*deflates* acceptance; the locally-balanced construction is the standard
repair, but that gloss is this program's reading, not Deshpande's.

**Cost, corrected upward.** The in-repo measurement (10.4x at p = 10,
53-56x at p = 50) is for *one node expansion* against *one classic move on
the same members*. A locally-balanced proposal over the whole birth/death
neighbourhood must score every leaf times every variable, i.e. `O(p n)` per
proposal since the leaves partition the data - roughly `p x (#leaves)`
classic moves. (The reverse normalizing constant is cheaper than it looks:
only the split leaf changes, and every leaf is scored against the same
residual within one move, so the unchanged leaves' sums carry forward.) The
break-even bar is therefore far above "p times the classic kernel", which
is where the survey set its kill line.

**Cost.** M-L, ~300 lines plus a gate arm; the scan itself is already paid
for, and the categorical analogue - scheduled when this was written - has
LANDED (`docs/plans/archive/grow-from-root-categorical-scan.md`; the TODO door is
recorded closed and removed from `TODO`). **The half worth as much as the
mixing gain** is the free by-product: because the scan scores the whole
neighbourhood, posterior functionals (variable inclusion, DART split
counts) can be averaged over the neighbourhood instead of the single
realized move - unbiased, variance never larger. That attacks the variance
of exactly the readouts that 3.1 and 3.3 corrupt, without needing the chain
to move at all.

### 4.6 Multi-split grow and prune ("twigs")

**What it is.** Grow attaches a whole chain of splits to a leaf rather than
a single split; prune removes an entire such chain. Kim and Rockova
penalize depth in the proposal with a geometric layer weight so it does not
always reach for the deepest layer.

**Which stickiness it targets.** Myopia (3.4) - it is the only candidate
that changes which states are reachable in one step, which is what the
theory says is required.

**Evidence.** Kim-Rockova's Theorem 5.3: under the assumption where plain
Bayesian CART is superpolynomial, twig-augmented Bayesian CART mixes in a
bound polynomial in n; Theorem 5.4: the informed twig variant is at most
linear in n. Crucially their Remark 9 says informedness alone does not
achieve this - the twig is doing the work. **But** the model is
one-dimensional dyadic single-tree Bayesian CART, and every method in their
experiments is single-tree: **there is no evaluation of twig moves in a
BART ensemble in print.**

**dbarts fit.** A depth-2 twig birth is two chained `tree.birth` calls with
both rules drawn from the prior, so both rule densities cancel exactly as a
single birth's does. The transition ratio needs the twig-length
distribution, the forward leaf-selection probability, and a reverse count
of nodes whose entire branch below is a twig - a new traversal, but a
simple one (`fillNoGrand` is the length-1 case). The practical worry is the
empty-leaf veto plus the depth prior (`base` 0.95, `power` 2): a depth-2
twig has three leaves that must all be occupied and is charged twice for
depth, so acceptance may be dominated by the prior penalty except in
exactly the XOR-like cases the move exists to catch. That is cheap to test.

**Cost.** M, ~200 lines plus a gate arm. Its falsifier is the cheapest and
most purpose-built of any candidate (hitting time from a stump to the
correct split pair on XOR at one tree), so it can be pulled forward if the
census says myopia dominates.

### 4.7 The cheap first stage that belongs to rotation

Worth recording separately because it is the right way to buy down the
largest implementation risk in the set: implement rotation as a **proposal
generator only**, with no acceptance ratio, and instrument on a stock
low-noise fit (i) the fraction of proposals that are structurally
admissible under dbarts' ordinal-only, interaction-guarded,
empty-leaf-vetoed constraints and (ii) the distribution of the
log-likelihood change for those. If the admissible fraction is very low or
the likelihood changes are systematically large and negative, dbarts'
rotation is not the near-neutral move Pratola's is and the premise fails
before the merge enumerator is written. Days, not weeks.

---

## 5. Set aside, and honestly why

### 5.1 Warming the sampler during burn-in and cooling it back down

**What it was.** Divide the likelihood part of the structural acceptance by
a temperature `T` that decays to 1 by the end of burn-in, so that early
sweeps accept more freely. Kept draws all run at `T = 1`, so the reported
posterior is exactly the posterior and the hot phase is just an
initializer - formally the identical validity argument
`grow-from-root.md` sec 3(a) makes for the shipped warm start. It was the
survey's top recommendation, at XS cost and bitwise-neutral by default.

**Why it leaves the recommended set. Four reasons, in order of weight.**

1. **It is pointed the wrong way for the failure it targets.** Acceptance
   is `priorRatio * transitionRatio * exp(dLogL / T)`. For a proposal that
   *is* an improvement, dividing by `T` **shrinks** the improvement while
   the depth-penalizing tree prior stays at full strength. So
   likelihood-only tempering makes good structural moves *less* likely to
   be accepted, and at high `T` its target is the prior over structures -
   which concentrates on depth 1-2, which is where a cold start already
   sits. The one house failure this was aimed at (section 3.2) is a chain
   that starts shallow and never *grows* the right structure; that is a
   directed-growth failure, and this construction trades directed growth
   for aimless diffusion. The most likely realized outcome is **inertness**:
   a cold start with fewer effective sweeps.
2. **Not one of the five sources cited for it evaluates the construction.**
   Tan et al.'s Theorem 7.3 analyses a *fixed* tempered chain, and the
   paper says so in the same paragraph: "the fact that we use a fixed
   temperature ... means that the stationary distribution for the sampler
   analyzed in Theorem 7.3 **is not the posterior**". Their Proposition 7.4
   is a concentration statement about the *tempered* posterior. Their
   Experiment 1 - the only empirical support - does **not** confine its
   schedule to burn-in: the authors' released code computes the schedule
   over `ndpost + nskip` and evaluates it on every iteration, gating only
   *recording* on the burn-in count. At their published settings, the
   temperature at the first kept draw is **2.82**, the mean over the 10,000
   kept draws is **1.91**, and only **5.5%** of kept draws sit at
   `T <= 1.1`. Their reported coverage gain is measured on tempered draws,
   where a flattened likelihood widens intervals that were under-covering.
   Angelopoulos and Cussens is verbatim "tempering (aka Metropolis-coupled
   MCMC)" - it is section 4.4 above, not this. `tgp` ships *importance*
   tempering, which keeps the tempered draws with weights, and is off by
   default.
3. **The temperature cannot be chosen from the literature.** The BCF
   combiner hands a forest weight `w_i * a^2` ([`formForestResponse`](../../src/bartcore/combiner.hpp), [`forestMultiplier`](../../src/bartcore/combiner.hpp)), so
   its acceptance exponent carries `a^2 / sigma^2`. At `a0 = 100`, `a^2` is
   1e4 and the recorded 5x-high sigma claws back only 25x, leaving the
   exponent roughly 400x too large - about 2.6 orders of magnitude. Divide
   by `T = 3` and nothing happens. Even the primary strata (`|a0|` 5-25)
   are 10-280x out. Tan's condition scales `T` with *n*; the dbarts failure
   is a signal-to-noise pathology at fixed *n*. **The right temperature is
   a function of the measured log-likelihood-difference distribution**,
   which nobody has measured - so this candidate cannot be *designed*, let
   alone built, before the census in section 6.
4. **Its falsifier was unrunnable and could not have answered the
   question.** The named harness (`burncurve.R`, `characterize.R`,
   `instrument.R`, `accept.R`, `poolreport.R`) is not in the repo and never
   was - `bcf-sigma-residual.md` says "Scratch, scripts, and rds output
   were run out-of-repo and are not preserved". `sbc.R`'s burn ladder is
   scoped to four response-family tiers and BCF is not one of them; there
   is no bias-versus-burn metric for BCF and no stratification by `|a0|`.
   No per-replicate standard error was ever recorded, so the proposed kill
   line ("beat it by more than the measured SE") had no SE. The control
   constant 1.52 was measured at `bartcore 6944811`, before the change-move
   balance fix and the variance-forest arc both changed the draw law. And
   holding total sweeps fixed makes "tempering does not help" and
   "tempering wasted 30% of the burn" read identically.

**What would have to be true to bring it back.** (a) The census shows that
high-signal-to-noise structural rejections are *close calls* - proposals
with small negative log-likelihood differences that a modest flattening
would let through - rather than proposals that are simply wrong. (b) If
built, temper the **whole** structural log-ratio (tree prior and
likelihood, transition ratio untempered), because that is what the one
supporting experiment actually ran and it is the form that escapes reason
(1). (c) The falsifier gains a matched-`T = 1`-sweeps arm so a null can be
attributed. (d) The harm check is treated as mandatory, not confirmatory -
this package has a confirmed +11.10% plateau posterior-mean error
regression in one cell from the last burn-phase change that carried the
same formal validity argument.

**One thing that is genuinely in its favour, and should be said.** The
error direction is benign. The shipped warm start errs *greedily against
the data* - deep, data-committed trees the fringe-only move set cannot
unwind. Tempering errs toward the prior - shallow trees, which births
repair cheaply. So the honest characterization is "no stationarity
obligation, a benign error direction, and therefore a likely inert rather
than harmful failure". That is a reason not to fear it, not a reason to
build it first.

### 5.2 Proposing a wholly different rule set for the same leaf partition ("restructure")

Wu, Tjelmeland and West propose replacing a tree's internal rules with a
different set inducing the *same* partition of observations into leaves,
by enumerating admissible (variable, threshold-interval) pairs at the root
and recursing. On a two-mode synthetic problem the non-restructure chain
produced "3,997 out of the 4,000 samples [with] exactly the same structure
and splitting variables as the starting tree"; with restructure the chain
converged in under 500 iterations versus more than 4,000 - and the
iteration counts are cost-normalized by construction, though the authors
note that figure is the *most* favourable of the several scalar functions
they tried.

**This is the weakest of the three declines, and the previously recorded
reason was wrong.** It was declined on the grounds that its cost grows with
the number of predictors and its mixing degrades with sample size -
"exactly where dbarts is hard, on both axes". The source contradicts the
first: section 4.2 supplies a subset-of-predictors variant (map only a
random subset `C`, with in-tree variables held at probability 1, plus the
matching acceptance correction) precisely to answer it, and the discussion
says "scale-up in terms of numbers of potential predictor variables is
immediate". The second is an untested conjecture in their closing remarks
("it is reasonable to believe"), on a paper whose largest dataset is
n = 683, p = 9, in a paragraph that also frames it as open and offers a
parallel multi-try mitigation.

**The honest reasons to set it aside.** It targets representation
multimodality, which the 75-tree ensemble self-averages away - so its
payoff is structural interpretation at small ensembles, a narrower prize.
And Pratola positions rotation as its cheaper local form in as many words:
the rotation's cost, "while greater than a simple birth/death proposal, is
much reduced compared to a more drastic restructure move such as the
proposal of Wu et al. (2007)". Recorded as a genuinely open door for
small-m interpretability work.

### 5.3 Replacing the per-tree kernel with a particle sampler

Particle Gibbs / conditional sequential Monte Carlo (Lakshminarayanan, Roy
and Teh) proposes a *complete* tree to fit the residual instead of editing
the current one.

The previously recorded reason - that their own tables show the local
sampler winning in the shallow-tree regime dbarts' prior produces - rests
on a table run at **one tree** with ten particles, so it cannot speak to an
ensemble. Their *ensemble* evidence (three real datasets, n = 2000, at the
paper's default of 200 trees) has the particle sampler ahead **3 of 3** on
raw effective sample size and 2 of 3 on effective samples per second,
losing the lowest-dimensional dataset to grow/prune by about 2.9x. There is
**no accuracy table anywhere in the paper**; the only predictive statement
is that all three samplers are "very similar".

**The honest reasons to set it aside** are architecture and cost, not
evidence: it replaces the per-tree kernel wholesale rather than adding a
proposal, it needs particle state of order (trees x particles x
observations), and its ensemble advantage is on raw sample size rather than
per second. If ever probed, the pre-registered read is effective samples
*per second* at 75 trees.

### 5.4 Grow-from-root as a periodic move rather than an initializer

Reset one tree per sweep and rebuild it with `growTreeFromRoot`, accepting
with a proper Metropolis-Hastings ratio. The builder already computes its
per-node candidate weights ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)), so the forward density is
nearly free; the reverse density requires replaying the builder's candidate
assembly along the *current* tree's construction path, the same cost again.
The reset/regrow/rebuild/redraw loop already exists ([`growForestFromRoot`](../../src/bartcore/chain.hpp))
and lacks only the acceptance filter. Reachability limited it to
ordinal-only forests until the scheduled categorical scan landed; that scan
has landed, and `growTreeFromRoot` ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)) now emits categorical
candidates of its own via `scanCategoricalPartitions` [`scanCategoricalPartitions`](../../src/bartcore/scan.hpp), so the limit is gone.

**Low priority, with a nearly free pre-check.** Instrument the realized
acceptance rate on a stock 75-tree fit. The residual-conditional posterior
of one tree inside a 75-tree ensemble is close to its prior, so an
independence proposal is likely to land rarely; below a couple of percent
it is pure overhead. That check costs a proposal generator and an hour.

**Erratum (2026-08-10, forest-specialization synthesis).** The middle
clause above is unsound and is withdrawn. "The target is diffuse" does
not imply "an independence proposal lands rarely": that inference is
correct for a proposal drawn from the *prior*, and `growTreeFromRoot` is
not one - its candidate weights are the prior factors *times* the
integrated likelihood ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)). The governing object
is the ratio, and it factorizes:
`pi(T)/q(T) = Z_root * prod_{w != root} [(1 - g_w) + g_w B_w]`, where
`g_w` is the CGM growth probability at `w` and `B_w` is the prior-averaged
split Bayes factor there. "Conditional close to prior" is exactly `B_w`
near 1, whence the ratio is near 1 and acceptance is near 1 - the opposite
of the recorded reading. Two new receipts, both taken after this section
was written: an exact enumeration of the single-predictor tree space under
the shipped arithmetic gives realized independence-MH acceptance of
**0.53-0.76** (the reviewing pass) and **0.84-0.97** (an independent
re-enumeration by the adjudicating pass on a different data-generating
process), in both cases *not* decaying as the tree space grows from 5 to
2950 trees; and Lakshminarayanan, Roy and Teh adjudicate the distinction
in one sentence - "proposing complete trees from the tree prior, however
these moves would be rejected, leading to slow mixing... The PG-BART
sampler succeeds not only because non-local moves are considered, but
because those non-local moves have high posterior probability."
[verified: AISTATS 2015 primary PDF, sec 1] Scope: both enumerations are
`m = 1`, one predictor, `n ~ 400`; they establish the mechanism and
predict nothing at 75 trees and `p = 10`.

**The conclusion nevertheless stands.** This section's *action* - low
priority, run the nearly free pre-check first - is unchanged, because the
acceptance rate at ship scale is still unmeasured and the pre-check is
still an hour. What changes is the recorded reason, which is now "we do
not know the acceptance rate" rather than "it will be low". Section 12.2
(B2) carries the derivation and section 12.6 schedules the pre-check as
Stage R0.

### 5.5 Continuous-time birth-death

Replace accept/reject with a continuous-time jump process where every jump
is accepted and implausible states simply have short waiting times.
Declined - on corrected grounds.

The previously recorded reason (the authors call it "too expensive") quotes
a naive intermediate construction they discard in the very next sentence.
The better reasons are: (i) the implemented substitute is a mixture of a
birth/death part and a rotate part at a **fixed** constant `alpha`, while
the value that makes the mixture exact is state-dependent and computing it
would restore the cost the split removed - and the paper never discloses
the `alpha` its benchmarks used; (ii) every published evaluation is single
tree at n = 300, with the authors scoping their methods to "reasonably
sized problems (e.g., thousands of observations, tens of variables)"; and
(iii) the change to time semantics (waiting-time-weighted, Rao-Blackwellized
estimators) would ripple through every consumer of this package. Recorded
fairly: in the one head-to-head table, the continuous-time arms beat every
reversible-jump arm on the mirror-pair functional.

### 5.6 Choosing the new split variable by correlation

Pratola's preconditioner proposes the replacement split variable with
probability proportional to its absolute correlation with the incumbent
(with correlations at or below 0.30 zeroed), and he reports it "leads to
much higher acceptance rates". It needs one p x p correlation matrix and a
proposal-density correction, and the per-side correction composition in
`changeMove` already has the right shape to carry it.

**Resolved in one direction, because it was previously both dismissed and
cited as half of another candidate's support.** It is a *variable*
proposal; perturb is a *cut* proposal. They are separable, they should be
built and tested separately, and Pratola's 65%/92% arm bundles them - so
that number cannot be credited to either alone. Low priority on its own
merits: dbarts already ships DART, which adapts split probabilities from
realized usage and occupies the adjacent design space from the prior side,
and the interaction between the two has never been measured.

### 5.7 The two levers already shipped, which are the baseline everything must beat

Neither is engine work, and both should be in the documentation regardless.

- **More trees.** Tan et al.: "Increasing the number of trees consistently
  dampens the trend in R-hat. Its effect on coverage and RMSE is
  ambiguous." dbarts defaults to `n.trees = 75L` ([`bart2`](../../R/bart.R)), below
  BART's classic 200. **Carry the caveat with it**: more trees dampens
  R-hat partly *because* the ensemble self-averages structural labels
  harder, so an improved R-hat at larger m is not by itself evidence that
  tree-space mixing improved.
- **More chains.** Ronen et al.'s own recommendation is to "increase the
  number of chains with the number of data points"; dbarts defaults to
  `n.chains = 4L` ([`bart2`](../../R/bart.R)).

No engine candidate should be measured against a single-chain, 75-tree
straw man.

---

## 6. The first tree-space falsifier, pre-registerable as written

This is the falsifier for the first *engine* candidate. The composition
probe (section 4.1) is a separate, cheaper study whose sketch lives in its
own section; the two are independent and section 7 recommends both.

One study, three stages, testing the same-variable cut move (section 4.2).
Its Stage 0 is the measurement four separate findings in this survey
resolve to, and it also freezes every threshold the later stages use - so
the house's pilot-then-freeze discipline is satisfied by construction
rather than bolted on.

**What has to be built.** Nothing reusable exists for the benefit contrast.
The nearest precedent is `benchmarks/R/change-fix-instrumentation.R` (285
lines), which did exactly this shape of engine instrumentation before:
environment-variable-gated CSV logging from data the move already computes,
the RNG stream untouched, logging switched on only after a silent burn-in,
and the engine patch reverted before commit. The correctness gate has two
in-repo templates in `swap-balance.R` (407 lines) and `bd-balance.R` (237
lines). The grow-from-root battery, which the *default* decision would
eventually need, is not in `benchmarks/` and per its own section 8 must be
reconstructed from the pre-registration rather than recovered.

### 6.1 Stage 0 - the move census (pilot; no kill criterion)

Instrument `moves.hpp` behind an environment gate and log, per structural
proposal: move type, target node depth, tree depth, interior-node count,
the integrated log-likelihood difference, the log prior ratio, the proposal
correction, and accepted/rejected. Run on a grid: dbarts defaults
(n = 5000, p = 10, m = 75, sigma = 1); the low-noise cell (sigma^2 = 0.1);
a wide cell (p = 50, 45 noise columns); and the causal-forest strong-scale
cell.

Four things nobody has, that everything downstream needs:

- **Per-move acceptance rates on this sampler.** Pratola's 4% / 18% / 25% /
  65% are the only numbers of this kind in print and they come from a
  different implementation with a birth/death-only baseline and an unstated
  tree prior.
- **The distribution of the log-likelihood difference among *rejected*
  structural proposals.** This is the fork that decides the whole program.
  If rejections cluster at small negative differences, the bottleneck is
  scale and the temperature family (sections 4.4, 5.1) is live. If they are
  overwhelmingly large and negative, the bottleneck is proposal accuracy
  and only better-aimed proposals help - which is the cut move's premise.
- **Change-move proposals *and* acceptances as a function of node depth.**
  This is the measurement section 3.3 has never had, and it decides whether
  rotation's motivation is real.
- **The cut displacement that puts a same-variable move at a target
  acceptance rate**, from which the window-width grid is set. Sample a few
  displacements per interior node rather than profiling the whole interval.

Stage 0 output freezes, before any Stage 2 contrast is looked at: the
window-width grid, the dosage grid, the per-replicate standard errors, and
every threshold.

Both addenda below were measured at the former default, where swap carried
0.1; their swap rows are the evidence the default dropped it to zero
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)).

**Addendum (2026-09-06): the first bullet is now measured.** A scaffold
build that instrumented the four [`resolveVetoRank`](../../src/bartcore/moves.hpp)
call sites for a different purpose - the occupancy veto's rejection budget -
also classified every structural proposal by move type and outcome, and
those counts carry the per-move acceptance rates this bullet asks for. At
the default mixture, the default prior, n = 2000, p = 10, 200 trees, 200
burn plus 500 sampled sweeps, one chain
([Measured occupancy rejection rate (2026-09-06)](empty-leaf-veto.md#measured-occupancy-rejection-rate-2026-09-06)):
birth accepts 22.56 percent of its proposals, death 25.95, change 14.79 and
swap 6.36, pooling to 18.61 percent over 140000 proposals. Two denominator
facts travel with those numbers. Change and swap make proposals that never
reach a score - 7.4 percent of change proposals and 73.0 percent of swap
proposals are no-ops - so per *scored* proposal the rates are birth 22.56,
death 25.95, change 15.98 and swap 23.55: swap's headline weakness is
almost entirely its no-op rate, not its acceptance ratio. The run was also
200 trees, not this section's m = 75 grid cell. The comparison to Pratola's
4 / 18 / 25 / 65 percent still cannot be pinned, because that section does
not record the noise level of its configurations and Pratola's own
birth/death rate moves 4.5x between his sigma^2 = 1 and sigma^2 = 0.1 cells.
The other three bullets - the log-likelihood difference among rejected
proposals, change acceptance against node depth, and the
displacement-versus-acceptance curve - remain unmeasured; none of them
falls out of a counter at the veto call site.

**Addendum (2026-09-07): the census ran, and the other three bullets are
now measured.** Instrumentation that survives in the tree -
[`cutProbe`](../../src/bartcore/moves.hpp) and the hooks around it, compiled
only under a macro - records one line per structural proposal, and the runner
[`runCell`](../../benchmarks/R/move-census.R) fits this section's four cells:
Friedman n = 5000, p = 10, m = 75 at sigma = 1 (`default`) and at
sigma^2 = 0.1 (`lownoise`); Friedman p = 50 with 45 noise columns (`wide`);
and the causal-forest cell (`bcf`), two forests with the prognostic surface
standardized to 8 sigma, which is the strong-|a| regime
[Burn-in under strong prognostic signal (2026-07-10)](bcf.md#burn-in-under-strong-prognostic-signal-2026-07-10)
records. Every cell ran 200 burn plus 500 sampled sweeps, one chain, one
thread, fixed seeds: 52500 proposals per single-forest cell, 87500 for BCF,
whose treatment forest carries 50 trees to the mean forest's 75. Sampling
takes 2.2 to 3.1 seconds per cell. Every table below is the 500 SAMPLED
sweeps; the burn column is acceptance per proposal over the 200 burn-in
sweeps, which is higher than the sampled rate in every cell pooled and in
every move but `default`'s swap.

Per move: proposals made, the share that never reached a score, and
acceptance on both denominators (percent).

    cell      move   proposals  no-op  accept  scored  burn
    default   birth       9613   0.00   10.05   10.05  14.95
              change     14736   2.72    4.07    4.19   5.98
              death       9403   0.00   10.24   10.24  12.59
              swap        3748  71.48    4.30   15.06   3.18
              all        37500   8.21    7.17    7.82   9.62
    lownoise  birth       9483   0.00    5.25    5.25  13.85
              change     14819   0.39    1.61    1.62   4.64
              death       9518   0.00    5.04    5.04   8.78
              swap        3680  70.14    1.71    5.73   2.71
              all        37500   7.04    3.41    3.67   7.77
    wide      birth       9886   0.00   11.64   11.64  21.88
              change     15079   4.99    5.37    5.65  12.66
              death       8826   0.00   12.70   12.70  25.26
              swap        3709  74.55    2.10    8.26   5.09
              all        37500   9.38    8.43    9.30  17.16
    bcf       birth      20985   0.00    5.91    5.91   9.58
              change     24948  34.32    2.80    4.27   4.63
              death      10361   0.00   11.84   11.84  15.29
              swap        6206  76.97    1.77    7.70   2.23
              all        62500  21.34    5.24    6.67   7.82

Swap's no-op share, 70 to 77 percent, reproduces the 73.0 percent the
2026-09-06 counts recorded. Change's no-ops are stumps rather than
unsatisfiable draws - 97 percent of them at `default` are trees with no
interior node to change - and BCF's 34.32 percent belongs to the treatment
forest alone: 8454 of that cell's 8562 change no-ops are forest 1, every one
a stump, which its shallow prior leaves common. The pooled rate, 7.17
percent at `default` against the 18.61 percent recorded at n = 2000 with
200 trees, moves with n and with the tree count in the direction that
record's own rows already show (18.61 at n = 2000, 10.92 at n = 10000,
9.63 at 50 trees).

The log-likelihood difference among REJECTED proposals, in log units.
Quantiles and shares are over the finite differences; the veto's -Inf, which
pooled over the moves is 0.07 to 0.25 percent of a cell's rejections, is
excluded from both.

    cell      move    rejected      q05      q25     q50    q75
    default   birth       8647    -3.33    -2.85   -2.33  -1.66
              change     13735  -402.60  -136.88  -56.61 -22.35
              death       8440  -342.32  -104.92  -43.00 -17.36
              swap         908  -335.04   -53.61  -22.53  -5.92
              all        31730  -307.19   -85.09  -25.78  -2.98
    lownoise  birth       8985    -3.98    -3.43   -2.81  -2.10
              change     14522 -1073.72  -330.77 -135.19 -52.47
              death       9038  -526.13  -182.60  -76.02 -27.71
              swap        1036  -727.59  -200.17  -79.03 -30.10
              all        33581  -677.65  -187.96  -54.87  -3.78
    wide      birth       8735    -3.38    -2.92   -2.43  -1.78
              change     13516  -351.32  -152.60  -80.14 -28.61
              death       7705  -282.57  -126.46  -65.04 -23.07
              swap         866  -246.36  -116.90  -34.75 -14.01
              all        30822  -277.43  -109.35  -32.70  -2.96
    bcf       birth      19744    -3.43    -2.98   -2.48  -1.74
              change     15687  -364.85  -160.31  -71.97 -27.76
              death       9134  -272.86  -111.78  -47.55 -19.87
              swap        1319  -306.78  -106.15  -41.19 -13.86
              all        45884  -262.38   -74.99   -9.84  -2.59

Share of those rejections within 1, 2 and 5 log units of zero (percent).

    cell        birth              change      death       swap
                1     2     5      1    2   5   1    2   5   1     2     5
    default   9.75 36.43 99.98  0.78 2.00 6.27 1.20 2.51 7.97 4.19 10.58 22.05
    lownoise  5.50 22.32 99.94  0.27 0.70 2.27 0.03 0.30 3.78 1.35  3.29  8.02
    wide      8.40 31.80 99.98  1.05 2.92 7.18 1.53 3.27 7.24 1.73  5.20 11.66
    bcf      10.14 32.06 99.83  0.85 1.99 5.53 0.34 0.92 5.35 1.67  4.10 12.06

Pooled over the four move types the shares are 3.43 / 11.76 / 32.69 percent
at `default`, 1.64 / 6.44 / 28.93 at `lownoise`, 3.27 / 11.25 / 33.60 at
`wide` and 4.77 / 14.76 / 46.23 at `bcf`.

Change proposals and acceptances by the depth of the node whose rule is
redrawn, scored proposals only. Depth is the target node's own depth; these
trees reach depth 5.

    cell      depth  proposals  accepted  accept
    default       0       9869       287    2.91
                  1       3483       238    6.83
                  2        795        60    7.55
                  3        152        14    9.21
                  4         35         0    0.00
                  5          1         1  100.00
    lownoise      0       7237        37    0.51
                  1       4906        97    1.98
                  2       1873        73    3.90
                  3        653        26    3.98
                  4         84         5    5.95
                  5          8         1   12.50
    wide          0      11152       489    4.38
                  1       2482       264   10.64
                  2        542        53    9.78
                  3        149         3    2.01
                  4          1         1  100.00
    bcf           0      10093       291    2.88
                  1       4348       255    5.86
                  2       1474       107    7.26
                  3        401        44   10.97
                  4         70         2    2.86

The same-variable cut move, priced but never run: at each interior node a
change proposal visited, the MH log ratio the move would have had with the
variable held and the cut displaced, evaluated on a snapshot that is
restored exactly and with no draw of its own. Acceptance is
mean min(1, exp(log ratio)) over the probes at that displacement; the two
signs are pooled by magnitude and are within a few hundred probes of each
other everywhere.

    cell      |1|   |2|   |4|   |8|   probes at |1|  median log ratio at |1|
    default  38.34 23.76 12.46  7.69          28578                    -1.83
    lownoise 26.43 14.75  6.88  3.83          29201                    -4.44
    wide     47.50 30.74 17.40 10.20          28514                    -1.04
    bcf      34.08 21.02 10.91  7.17          32439                    -2.40

Magnitudes 3, 5, 6 and 7 also appear, at 139 to 1140 probes against ~28000
for each power of two, because a requested displacement is clipped to the
node's descendant-valid interval. They come only from nodes whose interval
is narrower than the step asked for, so they are a biased subsample and are
not comparable to the unclipped magnitudes.

**Reading, against section 7's fork.** The rejected differences do not
answer as one distribution: birth's are small - median -2.33 log units at
`default`, and 99.98 percent of them within 5 - while change, death and swap
are large, change's median -56.61 at `default` and -135.19 at `lownoise`
with 0.78 and 0.27 percent within one log unit. Lowering the noise moves
every stream further out and birth's least. Change acceptance does not fall
with the depth of the node whose rule is redrawn: it rises from depth 0 to
depth 1 in every cell - 2.91 to 6.83 percent at `default`, 0.51 to 1.98 at
`lownoise`, 4.38 to 10.64 at `wide`, 2.88 to 5.86 at `bcf` - and keeps
rising to depth 3 in three of them, falling back at `wide`'s depth 3 on 149
proposals; depth 4 and below carries 84 proposals or fewer in every cell and
cannot be read. What concentrates at depth 0 is the proposal mass, 69
percent of scored change proposals at `default`, since most interior nodes
are the root. Within the depths these trees reach, that is the opposite sign
to the one section 3.3 predicts, and it leaves 3.3 untested below depth 5. A
cut displacement of one position sits in the 25 to 50 percent band in every
cell (38.34, 26.43, 47.50, 34.08); two positions holds only at `wide`
(23.76, 14.75, 30.74, 21.02).

The census build is not the shipped one: the instrumentation compiles only
under `-DBARTCORE_MOVE_CENSUS`, which reaches the compile through a user
Makevars appending to `CPPFLAGS` under `R_MAKEVARS_USER`, installed to a
private library so the ordinary one is untouched. The runner's header
carries the exact commands and the record format.

**Addendum (2026-09-07): the census re-ran at the two-move kernel.** The
swap move is removed
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)),
so `metropolisJumpForTree` now dispatches to birth/death or change only, at
`birth_death = 0.6, change = 0.4, birth = 0.5`. Same instrumentation, same
runner, same four cells, same seeds, same 200 burn plus 500 sampled sweeps:
52500 proposals per single-forest cell, 87500 for BCF. Sampling takes 2.2 to
3.1 seconds per cell, same as the three-move run. Every table below drops
the swap row and column outright; there is no swap move to report.

Per move: proposals made, the share that never reached a score, and
acceptance on both denominators (percent).

    cell      move   proposals  no-op  accept  scored    burn
    default   birth      11573   0.00   10.39   10.39   15.52
              change     14919   2.46    3.77    3.87    5.85
              death      11008   0.00   10.88   10.88   13.17
              all        37500   0.98    7.90    7.98   10.95
    lownoise  birth      11291   0.00    5.46    5.46   12.26
              change     14969   0.73    1.66    1.67    4.54
              death      11240   0.00    5.26    5.26    8.72
              all        37500   0.29    3.88    3.89    8.15
    wide      birth      12136   0.00   12.62   12.62   21.29
              change     14874   7.56    6.06    6.55   13.17
              death      10490   0.00   14.43   14.43   23.93
              all        37500   3.00   10.53   10.85   18.70
    bcf       birth      25525   0.00    5.63    5.63    8.75
              change     25086  36.96    2.88    4.57    4.53
              death      11889   0.00   12.13   12.13   15.04
              all        62500  14.84    5.76    6.77    8.27

Pooled (scored) against the three-move figure: default 7.98 vs 7.82,
lownoise 3.89 vs 3.67, wide 10.85 vs 9.30, bcf 6.77 vs 6.67. Burn is higher
than sampled in every cell and every move now, with no exception - the
three-move run's one exception was `default`'s swap, which is gone.

The log-likelihood difference among REJECTED proposals, in log units
(cell/move as above; the veto's -Inf is excluded, 0.13 to 0.27 percent of a
cell's rejections pooled, against 0.07 to 0.25 before).

    cell      move    rejected         q05      q25       q50     q75
    default   birth      10371     -3.35    -2.84    -2.31   -1.64
              change     13989   -299.21  -122.47   -62.34  -23.84
              death       9810   -247.06  -103.28   -48.93  -18.33
              all        34170   -229.08   -84.68   -24.68   -2.82
    lownoise  birth      10675     -4.05    -3.51    -2.87   -2.15
              change     14611   -809.48  -336.50  -143.45  -52.71
              death      10649   -520.59  -201.74   -77.13  -27.04
              all        35935   -589.19  -188.75   -46.95   -3.64
    wide      birth      10604     -3.47    -2.97    -2.48   -1.83
              change     12848   -375.24  -149.74   -63.16  -23.12
              death       8976   -279.13  -123.94   -49.63  -21.21
              all        32428   -275.35   -91.68   -23.32   -2.78
    bcf       birth      24088     -3.62    -3.19    -2.63   -1.91
              change     15091   -443.94  -164.05   -75.90  -27.69
              death      10447   -293.45  -108.22   -52.28  -19.51
              all        49626   -272.02   -67.93    -4.04   -2.60

Pooled median against the three-move figure: default -24.68 vs -25.78,
lownoise -46.95 vs -54.87, wide -23.32 vs -32.70, bcf -4.04 vs -9.84 -
every cell's rejected mass sits closer to zero with swap's long left tail
gone. The fork still answers the same way: birth's rejections stay small
(median -2.31 to -2.87 across cells) while change's and death's stay large
(change's median -62.34 to -143.45).

Share of those rejections within 1, 2 and 5 log units of zero (percent).

    cell        birth              change            death
                1     2     5      1    2    5      1    2    5
    default   10.44 37.45  99.96  0.87 2.12 6.52  1.19 2.55  8.44
    lownoise   5.67 21.28 100.00  0.25 0.72 2.20  0.07 0.46  3.66
    wide       7.60 30.69  99.92  1.34 3.72 9.67  2.31 4.66 10.15
    bcf        8.36 27.33  99.43  0.93 2.14 5.66  0.38 1.18  6.20

Pooled over the three move types, against the four-move figure: 3.86 /
12.94 / 35.37 percent at `default` (was 3.43 / 11.76 / 32.69), 1.80 / 6.73 /
31.61 at `lownoise` (was 1.64 / 6.44 / 28.93), 3.65 / 12.79 / 39.27 at
`wide` (was 3.27 / 11.25 / 33.60), 4.42 / 14.16 / 51.26 at `bcf` (was 4.77 /
14.76 / 46.23).

Change proposals and acceptances by the depth of the node whose rule is
redrawn, scored proposals only.

    cell      depth  proposals  accepted  accept
    default       0      10125       249    2.46
                  1       3404       233    6.84
                  2        711        60    8.44
                  3        268        19    7.09
                  4         43         1    2.33
                  5          1         1  100.00
    lownoise      0       8187        37    0.45
                  1       4197       123    2.93
                  2       1668        49    2.94
                  3        547        24    4.39
                  4        184        11    5.98
                  5         72         2    2.78
                  6          2         0    0.00
                  7          2         2  100.00
    wide          0      10282       618    6.01
                  1       2866       221    7.71
                  2        511        51    9.98
                  3         76         8   10.53
                  4         11         2   18.18
                  5          3         1   33.33
    bcf           0       9947       347    3.49
                  1       3863       240    6.21
                  2       1475        98    6.64
                  3        464        28    6.03
                  4         65        10   15.38

Depth 0 to depth 1 rise, against the three-move figure: default 2.46 -> 6.84
(was 2.91 -> 6.83), lownoise 0.45 -> 2.93 (was 0.51 -> 1.98), wide 6.01 ->
7.71 (was 4.38 -> 10.64), bcf 3.49 -> 6.21 (was 2.88 -> 5.86). The rise from
depth 0 to depth 1 still holds in every cell.

The same-variable cut move, priced but never run, unchanged in method from
the three-move addendum.

    cell         |1|    |2|    |4|    |8|  probes at |1|  median log ratio at |1|
    default    40.65  24.06  11.94   6.81          28682                    -1.53
    lownoise   27.15  13.73   5.63   3.07          29221                    -4.83
    wide       51.38  34.82  20.77  12.91          27281                    -0.76
    bcf        33.34  21.09  12.31   8.09          31266                    -2.46

Displacement-1 acceptance and median log ratio, against the three-move
figure: default 40.65 / -1.53 (was 38.34 / -1.83), lownoise 27.15 / -4.83
(was 26.43 / -4.44), wide 51.38 / -0.76 (was 47.50 / -1.04), bcf 33.34 /
-2.46 (was 34.08 / -2.40).

Nothing moved beyond what the mixture change alone predicts: birth and
death's combined proposal share rises by 9.4 to 10.4 points in every cell,
matching the 9.8 to 10.0 points swap held (change's own share is flat, 39.3
to 40.2 percent under either kernel), and every move's conditional numbers
above - accept/scored, the rejected-difference quantiles, the depth curve,
the |1| cut acceptance - sit within a couple of points of the three-move
figures in every cell, `bcf`'s birth taking more of the gain than death only
because its two forests already carried different tree counts (75 against
50), not because the kernel changed.

**Addendum (2026-09-07): the generator-only probes.** Four measurements the
two proposal-brainstorm rounds named as their one-day falsifiers -
[15.3 Cross-lens ranking](#153-cross-lens-ranking) rows 1 to 3 and
[16.3 Ranking](#163-ranking) row 2 - are now taken, all generator-only:
nothing draws, nothing changes the RNG stream, and every hook computes, logs
and restores. [`census::nogProbe`](../../src/bartcore/moves.hpp) rides
`changeMove`'s target node and prices the closed rule neighbourhood a
collapsed Gibbs draw would score there - every (available ordinal variable,
admissible cut) pair weighted by the cut scan's collapsed marginal, the
node's own rule prior `1/|SI|` and the two `log(1 - growth(child))` terms
that are exactly `changeMove`'s below-node prior, the split-variable factor
itself cancelling since it is uniform over the available set - and logs the
weight entropy, the incumbent's share and its rank, both jointly over the
available variables and restricted to the incumbent's own.
[`census::deathProbe`](../../src/bartcore/moves.hpp) rides
`birthOrDeathMove`'s death branch and weighs every nog node by the
merged-leaf marginal ratio against the uniform pick the kernel actually
made. [`census::perturbProbe`](../../src/bartcore/moves.hpp) logs the
perturb move's node and its signed displacement (target minus current), and
[`census::treeShape`](../../src/bartcore/moves.hpp), hooked at both forest
sweep loops in [chain.hpp](../../src/bartcore/chain.hpp), logs every tree's
settled leaf count once its move has settled. All four compile only under
`-DBARTCORE_MOVE_CENSUS`; the default build carries none of it, and the
existing `p` and `d` records are byte-identical before and after the
instrumentation, 208126 lines.

The runner grows a fifth cell, `c1`: the He and Hahn independent design at
n = 10000, p = 30, Trig+poly, kappa = 1, 75 trees
([10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)),
the surface battery's primary-benefit cell, at the same 200 burn plus 500
sampled sweeps, one chain. It ran twice: once at the shipped mixture
(birth_death 0.6, change 0.4), which is what the nog, death and shape
tables below read, and once under
[`censusProposalProbs`](../../benchmarks/R/move-census.R)'s
perturb-carrying mixture (birth_death 0.5, change 0.34, perturb 0.16) - the
only way a signed displacement is ever recorded, since perturb never fires
at its shipped zero - which the perturb table reads; `bcf`'s treatment
forest refuses a non-default `proposal.probs` outright, so it has no
perturb row. The `default` cell reproduces the second 2026-09-07 addendum's
per-move table exactly.

Closed rule neighbourhood at a nog node, sampled sweeps only (`nog%` the
share of interior nodes that are nog, `target-nog%` the share of change
proposals landing on one, `H` the median weight entropy in nats, `P(inc)`
the median incumbent probability under the normalized weights, `top%` the
share of nog proposals where the incumbent already holds the top weight;
"joint" ranges over the available variables and their cuts, "cut" over the
incumbent variable's cuts alone):

    cell        nog%  target-nog%   joint H  P(inc)  top%   cut H  P(inc)  top%
    default     57.6      71.9        0.911  0.361   50.4   0.894  0.379   52.0
    lownoise    48.5      62.7        0.334  0.737   60.1   0.328  0.746   60.6
    wide        64.1      76.5        1.066  0.079   27.3   1.059  0.110   29.9
    bcf f0      52.2      65.0        0.689  0.522   57.1   0.682  0.539   57.9
    bcf f1      98.3      99.1        1.659  0.130   37.6   1.579  0.164   39.5
    c1          65.2      77.3        6.425  0.0017  27.4   3.426  0.034   31.9

`c1`'s incumbent has a median rank of 26.5 of up to 3000 candidates jointly
and 4 of 100 on the cut axis alone - the Gibbs step has the most to buy
there, in entropy and in how far the current rule sits from the mode.
`lownoise` is the opposite pole: the smallest joint entropy and the highest
`P(inc)` and top-share of any single-forest cell, so the incumbent already
sits close to what a Gibbs draw would pick - the least there is to buy. That
is where 15.3's falsifier ("`P(incumbent)` near 1 at `lownoise` kills it")
was aimed; measured, it is 0.737, not near 1, so the falsifier's kill
condition does not fire anywhere in this grid, though the ordering it
predicted holds. A temporary assertion cross-checking the scan's incumbent
entry against the kernel's own cached branch score agreed to 3.7e-13 over
the run, so the probe is pricing the same rule the kernel actually holds.

Informed death (15.3 row 3): among death proposals seeing two or more nog
nodes - 5.7 (`c1`) to 32.3 (`lownoise`) percent of death proposals by cell,
the rest (68 to 94 percent) seeing exactly one, where the uniform pick has
nothing to be wrong about - the normalized merged-marginal weight vector is
a point mass: median entropy 3e-5 nats at `c1` down to 5e-20 at `lownoise`,
median max weight 1.000000, the uniform pick's median rank 1 in every cell.
It already sits at the weighted mode 83.1 to 97.1 percent of proposals
overall, falling to 47.7 to 49.9 percent - indistinguishable from choosing
at random - once there are two or more candidates. There is nothing here
for a weighting to inform: one nog node's merged statistic dominates the
others so completely that the uniform draw already lands on it almost every
time it could matter.

Perturb signed runs, from the perturb-carrying mixture, among consecutive
accepted displacements at one node (`pairs` the count of such consecutive
pairs; `reversible%` the null a reversible walk implies):

    cell      pairs  same-direction%  reversible%
    default    1266             36.6           50
    lownoise    528             38.1           50
    wide       1265             38.7           50
    c1         1611             43.1           50

Every cell sits below the reversible null, not above it - accepted
displacements tend to REVERSE, not continue - and the extension hazard (the
chance a streak already k long extends by one more) is flat in streak
length rather than rising, so there is no persistence to exploit.
16.3 row 2 pitched a lifted cut displacement as a gain if same-direction
runs exist; measured, the opposite holds, and a lift - which spends its
whole saving forcing continuation in one direction - would spend it
fighting a chain that already prefers to turn around. Perturb's own scored
acceptance in this run is 29.7 (`default`), 13.2 (`lownoise`), 34.2
(`wide`) and 54.2 (`c1`) percent.

Leaves per tree, sampled sweeps, shipped mixture (mean of the per-sweep
means; quantiles over trees and sweeps together):

    cell      mean   q05  q50  q95   max  stump%
    default   2.83     2    3    5     8     2.4
    lownoise  3.79     2    3    8    13     0.5
    wide      2.53     1    2    5     7     7.5
    bcf f0    3.36     2    3    6     9     1.6
    bcf f1    1.11     1    1    2     4    89.5
    c1        2.52     1    2    4     8     6.1

Every mean sits just above 16.2's Jensen bound for its cell (2.44 / 2.82 /
2.34 / 2.59) - shipped trees really do carry two to three leaves, `bcf`'s
treatment forest closer to one.

`c1` as a census cell, shipped mixture, sampled sweeps: scored acceptance is
25.9 (birth), 20.7 (change), 29.2 (death), pooled 24.9 percent, against the
`default` cell's 7.98. `c1` is not a sticky regime by this measure - it is
the loosest of the five cells on every move - which means
[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s
minimum effective sample size of 2 out of 2500 kept draws cannot be a
structural-acceptance deficit: the trees move freely and the deficit
survives it untouched. That is what elevates 16.3 row 1's frozen-structure
ESS test over the rest of this round's ranking - if leaf values alone,
structure held fixed, already carry an ESS this low, the bottleneck is the
fibre row 1 targets, not which moves the kernel proposes.

Not run: 16.3 row 4's pair-transfer probe. Pricing a transfer needs the
residual net of the other `m - 2` trees and the partner tree's own leaf
statistics, neither of which `changeMove` or `birthOrDeathMove` sees; it
needs a [chain.hpp](../../src/bartcore/chain.hpp) block scoring the joint
`L_j + L_k` system row 4's refutation already describes, not a hook inside
one of the moves this round's other three probes reused.

### 6.2 Stage 1 - correctness (`perturb-balance.R`, new)

A per-kernel exact-posterior gate on the **within-variable cut
distribution**. `change-balance.R` cannot serve: its pass/fail statistics
are z-tests on root-split-*variable* marginals, and its cut distribution is
computed but only printed (`cutReport`, lines 574-593) with no threshold,
no z, and no verdict. A cut move never changes a variable, so a defect in
its interval, its cut law, or its window symmetry lives precisely in the
quantity that gate does not gate - in a package whose last change-move
defect survived its entire history because the existing gates could not see
it.

Design: a single-tree enumerable problem, exact posterior over (variable,
cut) from the region dynamic program `change-balance.R` already implements,
run with the mixture set perturb-dominant so the realized cut distribution
*is* the kernel's stationary law. Gate on per-cut z against the exact
probabilities, Holm-corrected across cuts. **The test problem must place
cuts near the interval edges**, because the local-window form's correction
`|W(c)| / |W(c')|` is exercised only when the window is clipped. Failure
means fix the kernel, not drop the candidate.

### 6.3 Stage 2 - benefit, with matched exposure

Four arms, matched seeds, paired (data seed `BASE_SEED + s`, sampler seed
`s`, shared `s` across a pair).

| arm | what it runs | purpose |
|---|---|---|
| A | today's mixture, N sweeps | control |
| B | perturb takes weight `w` from change, N sweeps | same proposal count as A |
| C | today's mixture plus an OpenBT-style perturb pass, N sweeps | extra throughput |
| D | today's mixture, N x (measured cost ratio of C to A) sweeps | C's extra compute spent on plain sweeps |

This is the design the previous falsifier lacked. It separates the three
outcomes that otherwise read identically:

- **B beats A**: the cut move is worth more than the change move's marginal
  proposal. Help, at no extra cost.
- **B is flat and C beats D**: it helps, but only as extra throughput.
  Help, at a cost.
- **C is flat against D**: **inert** - any apparent gain was just more
  compute.
- **B is worse than A**: it displaces something better. Harm.

Metrics, all frozen at Stage 0. Primary: 90% pointwise interval coverage of
the true mean on held-out points, in the low-noise cell - the regime where
the failure is established and coverage is the published failure. Coverage
is preferred to R-hat, which this house has established is not gateable.
Secondary: per-move acceptance rates (Stage 0's instrumentation kept on),
integrated autocorrelation time on held-out error and on sigma, and the
pooled between-chain standard deviation of time-averaged variable inclusion
- the statistic that replaced R-hat in the grow-from-root study, whose
probe standard deviation per replicate was 0.0011-0.0015 against a cold
level of 0.0075, i.e. it has power at 12 replicates.

Cells: (i) low-noise Friedman, n = 5000, p = 10, m = 75, sigma^2 = 0.1;
(ii) the ship default at sigma = 1, which must not regress; (iii) Friedman
p = 50 with 45 noise columns, the regime that killed the warm-start
default; (iv) the causal-forest strong-scale cell.

Grids: window width at 10% (OpenBT's shipped value), 85% (Pratola's own),
and 100% (the uniform form whose correction is identically 1). Dosage from
the survey's 0.04-0.16 attempts per tree per sweep up to OpenBT's
throughput (~0.9 attempts per interior node per sweep, i.e. roughly 3-6 per
tree per sweep at 75 trees). Both grids frozen against Stage 0's measured
acceptance-versus-displacement curve; the confirmatory contrast runs at ONE
(width, dosage) point chosen by a rule written down before Stage 1.

Null control: an all-categorical cell, where a first version's ordinal-only
cut move provably cannot act. Every metric must sit within Monte Carlo
error between arms. If it does not, the estimator family is void.

### 6.4 Kill criteria, pre-registered

**KILL** if, at the Stage-0-selected (width, dosage) in the low-noise cell,
arm B does not improve coverage over arm A by more than 4x the measured
per-replicate standard error, AND arm C does not beat arm D on the same
metric by the same margin - over at least 20 matched pairs, with a
mandatory fresh-seed re-run of any single flagged cell before a flag counts.

**KILL the default question independently** on any per-cell harm rejection
(Holm-corrected, point estimate beyond its frozen margin, in at least two
cells) on plateau prediction error in the noise-heavy or the large-n
stratum - mandatory, not confirmatory. The last change of this shape that
carried a clean formal argument cost this package a confirmed +11.10%
plateau posterior-mean error in one cell.

**Asymmetry to state plainly.** Passing Stage 2 justifies shipping the move
**opt-in**, at default weight 0, bitwise-neutral. Flipping any default is a
separate decision that needs the grow-from-root harm battery, which does
not exist in `benchmarks/` and must be reconstructed. Pricing the diff is
not pricing the decision.

### 6.5 Cost, honestly

| piece | size | note |
|---|---|---|
| Stage 0 instrumentation | ~100 lines, reverted before commit | pattern exists in `change-fix-instrumentation.R` |
| Stage 0 driver | ~250 lines | new |
| the kernel | ~120-160 lines, +15 for the window ratio | 6 files; default weight 0 |
| `perturb-balance.R` | ~400-600 lines | the single largest artifact; not skippable |
| Stage 2 harness | ~400 lines | reuses the matched-pair seed idiom |
| compute | Stage 0 hours; Stage 2 on the order of a day | 4 arms x 4 cells x 20 pairs plus re-runs, against the grow-from-root study's 1.93 h floor |

Total: **a week of implementation, not a day**. The kernel is the cheap
part and the survey priced only the kernel.

---

## 7. Recommended next step

**Run two measurements first, both of which change no draw and need no
engine change, and decide afterwards.** They are independent and can run in
either order or together.

**Measurement 1: the move census (section 6.1).** It has no kill criterion,
it costs hours, and four separate open questions in this document resolve
to it:

- whether high-signal-to-noise rejections are close calls (temperature is
  live) or wrong proposals (only better-aimed proposals help);
- whether the change move's acceptance really collapses with node depth,
  which is rotation's entire motivation and has never been measured;
- what the per-move acceptance rates actually are on this sampler, a number
  that does not exist anywhere;
- what window width puts a cut move at a workable acceptance rate.

The strong recommendation of this document is that **no candidate should be
designed, not merely built, before that measurement exists**. The survey's
top pick failed partly because its temperature could not be chosen without
it.

**Measurement 2: the composition probe (section 4.1's falsifier sketch).**
Also pure measurement, also no engine change, and every arm already ships -
`bart()`, linear leaves, `setOffset` with a conjugate block, and stan4bart.
It answers a question that changes what the engine candidates are worth:
does moving the smooth share of the signal out of the forest actually make
the forest's own sampling behave better, or does the ridge between the two
blocks eat the transfer? It is worth running first for two reasons beyond
cost. A positive result on the *inner* (leaf-level) variant would be an
improvement users can have today, with no new kernel and no new gate. And
it would settle a question this package's own ecosystem paper left open:
stan4bart's section 4.5 argues that letting the two components overlap is
beneficial parameter expansion and then says, in as many words, "More
research will need to be performed to confirm this." Four years on, no one
has - and the semiparametric-BART literature has meanwhile taken the
opposite position and built constraints against overlap. The probe
adjudicates that disagreement, and dbarts plus stan4bart is the only place
it can be run cheaply.

**Then decide.** If the census says proposals are wrong rather than close:
build the same-variable cut move and run sections 6.2 and 6.3. If it says
rejections are close calls: reopen the temperature family, and reopen it at
**heated companion chains with private ladders** (section 4.4) rather than
at annealed burn - the kept draws stay exact, the published support is
real, and the architectural objection does not hold. If it says change-move
acceptance collapses sharply with depth: run rotation's cheap first stage
(section 4.7) next instead. If the composition probe says depth transfer
works and survives the ridge, then every tree-space candidate is worth less
than it looks, because the cheapest way to fix a sticky forest is to give
it less to do.

---

## 8. Things this program should not re-derive

- **The change move's root sensitivity is not a bug.** It is the
  fixed-skeleton design working as specified ([`changeMove`](../../src/bartcore/moves.hpp)), and the
  same mechanism would pin a cold-started tree's root too; a cold start
  only escapes it by passing through shallow states.
- **Duplicate columns cannot separate a locked chain from a mixing one**
  (their likelihood ratio is exactly 1). Use the XOR construction; keep
  duplicates as a null control whose failure is itself informative.
- **At 75 trees no structural statistic can detect mode collapse** at
  feasible replicate counts. Structure probes run at one tree on
  purpose-built scenarios.
- **R-hat on held-out error is not gateable** at those probe sizes.
- **Neither formal mixing lower bound binds dbarts' kernel**, and both name
  a data-fitted initialization as the remedy - which this package measured
  and killed as a default. The empirical half of Ronen et al. *does* apply
  (full move set, real data, this package by name). And Tan et al. found
  change and swap add nothing measurable over grow/prune on their battery -
  so the case for a *new* move cannot rest on "dbarts has four moves".
- **A same-variable ordinal cut redraw carries acceptance correction
  exactly 1** under shipped machinery ([`changeMove`](../../src/bartcore/moves.hpp)). This is the
  single most load-bearing code fact in the document and it has now been
  verified twice, independently.

## 9. What this survey could not settle

- **Whether 200 trees governs Pratola's 25% / 96% rotation result.** The
  section says only "the same dataset"; a full-text search of that section,
  its figure captions, and the discussion finds no restatement of the tree
  count. It is the only reasonable reading, but it is inference, and
  rotation's rank leans on it.
- **What "ESS x1" measures in the one table that isolates the cut move.**
  Its 10x sits on a mirror pair whose two splits induce the identical
  partition, so it could be dominated by cut-value randomization that
  changes no fit.
- **The per-replicate standard error of the causal-forest burn curve** and
  **which 10 replicates it was measured on.** Both unpreserved; any arm
  there must re-measure both.
- **Whether dbarts' rotation would be near-likelihood-neutral** the way
  Pratola's is, given the ordinal-only scoping, the interaction guard and
  the empty-leaf veto. Section 4.7 is the cheap way to find out.
- **Ronen et al.'s root-bottleneck number on a corrected change move.**
  Theirs ran on dbarts 0.9-22, which predates the change-move
  detailed-balance fix. dbarts is the only package that could re-measure it.
- **Whether composing with a parametric block improves tree-space mixing at
  all.** Nobody has measured it - not in print (a full-text audit of ten
  papers found no ESS, autocorrelation or R-hat for any composed
  parametric-plus-BART sampler) and not here. The house has measured the
  *hazard* (the 6x alternation penalty in the grouped surrogate) but never
  the *benefit*, and stan4bart, which is the pattern, reports no mixing
  diagnostic at all: its vendored sampler surfaces no divergence, treedepth
  or energy hook, and the Stan-layout rows it still writes are constant
  zero placeholders (`stan4bart/docs/design/walnuts.md`, "Diagnostics:
  dropped, not conditionally"). The literature does not even agree on the
  sign - stan4bart's section 4.5 argues overlap is beneficial parameter
  expansion and says "More research will need to be performed to confirm
  this"; CSP-BART and Bhandari et al. argue overlap is harmful
  non-identifiability and build constraints against it. This is the single
  largest open question the survey found, the cheapest to close, and the
  one where closing it would settle a question this package's own
  ecosystem paper left open in 2022.
- **Two blocked sources that could contain a hit.** The published BCF text
  with its discussion and rejoinder (Project Euclid 403), where a
  discussant might raise mu/tau mixing; and Zhang et al., Statistics in
  Medicine 45 (2026) e70593, whose abstract describes near-collinearity and
  instability from cluster-level covariates in a BART-plus-mixed-effects
  model. Neither was read. The item-7 NOT-FOUND is stated against
  everything that could be fetched, not against everything that exists.
- **Whether the forest-versus-parametric ridge is as bad as the
  forest-versus-group-intercept ridge.** The 6x is measured on a per-group
  intercept, where the forest can alias each group's mean with a free
  constant. A smooth global term is a more constrained competitor, so the
  ridge could be milder - or, because a forest can render a smooth surface
  many ways, worse. The in-house record's own authoritative critique
  declines to bound this in either direction for the ranef case, and the
  same uncertainty applies here.

---

## 10. Citation ledger

Every source below was fetched and read during this arc. "Full text" means
the PDF or source was downloaded and the cited passage read directly. Every
load-bearing claim was verified twice: once by the survey, once
independently by the adjudication pass, which also fetched the released
experiment code where a paper's prose was ambiguous.

| # | Source | What was verified | Where |
|---|---|---|---|
| 1 | Pratola, Bayesian Analysis 11(3):885-911, 2016 | Full text (arXiv v1/v4): sec 2.2 stuck Friedman ("with birth/death proposals only"), no tree-prior hyperparameters anywhere, sec 3 perturb + the 85% window, sec 3.1 the 0.30 correlation cutoff, sec 4 rotation, sec 5.1/5.2 results, sec 6 discussion | arXiv 1312.1895 |
| 2 | Wu, Tjelmeland, West, JCGS 16(1):44-66, 2007 | Full text: sec 4.1 restructure, sec 4.2 subset-of-predictors variant + its MH correction, sec 5.1/5.2 results and their cost normalization, Discussion scaling conjecture | www2.stat.duke.edu/~mwest/MWextrapubs/Wu2007.pdf |
| 3 | Mohammadi, Pratola, Kaptein, JMLR 21(201), 2020 | Full text: Table 1 transcribed in full, the mirror-pair construction (eqs 13-16), Appendix B and the fixed-alpha substitute, scope sentence, M=1 n=300 | jmlr.org/papers/v21/19-307.html |
| 4 | Lakshminarayanan, Roy, Teh, AISTATS 2015 | Full text: Tables 1/2 ("We fix m = 1", C = 10), Tables 4/5 (m = 200 default, n = 2000), absence of any accuracy table | proceedings.mlr.press/v38/lakshminarayanan15.pdf |
| 5 | Kim, Rockova, EJS 19(2):3041-3067, 2025 | Full text of BOTH the published EJS PDF and arXiv 2306.00126: Theorem 5.1 (superpolynomial), 5.3/5.4 (twigs), sec 3.1 geometric layer weight, Remark 9 + its footnote, the sec 5.1 initialization preamble, the EJS-only DP-initialization paragraph, and the renumbering | doi.org/10.1214/25-EJS2397 ; arXiv 2306.00126 |
| 6 | Ronen, Saarinen, Tan, Duncan, Yu, arXiv 2210.09352 | Full text: sec 4.2 (<0.2%, full move set, 1 tree, 4 PMLB sets), the dbarts 0.9-22 implementation statement, sec 1.3's non-transfer sentence, Appendix A.3 (figure only), sec 5 caveats and recommendations | arXiv 2210.09352 |
| 7 | Tan, Ronen, Saarinen, Yu, arXiv 2406.19958v2 | Full text: eq 10, eq 13, Thm 7.3 + its "is not the posterior" paragraph, Prop 7.4, the T pincer, sec 9 findings, Appendix L.6. **Plus the released code** (both repos out-of-repo here, so cited by file and function rather than by line): `bart-comp-efficiency` @ 0589240 (`runner.py`, `utils.py`) and `bart-playground` @ e686d34 (`samplers.py`) | arXiv 2406.19958 ; github.com/yanshuotan/bart-comp-efficiency ; github.com/yanshuotan/bart-playground |
| 8 | Zanella, JASA 115(530):852-865, 2020 | Abstract + the balancing-function conditions and the no-dominance statement | arXiv 1711.07424 |
| 9 | Angelopoulos, Cussens, ICML 2005 | Full text: sec 5 "tempering (aka Metropolis-coupled MCMC)", the 4-chain ladder, swap-every-iteration, cold-chain-only collection, Table 6 (3 seeds, 15/16), PIMA accuracy, the CGM-1998 block quote | icml.cc/Conferences/2005/proceedings/papers/003_Tempering_AngelopoulosCussens.pdf |
| 10 | Deshpande, flexBART, arXiv 2211.04459 | Full text: sec 2.2 informed-proposal objection, appendix B2/B3 | arXiv 2211.04459 |
| 11 | Zhang, Huelsenbeck, Ronquist, Syst Biol 69(5):1016-1032, 2020 | Record + verbatim abstract (order-of-magnitude faster convergence; dataset-dependence caveat) | academic.oup.com/sysbio/article/69/5/1016/5716338 |
| 12 | Gramacy, Taddy, JSS 33(6), 2010 | Record + abstract; importance tempering, `itemps = NULL` by default | jstatsoft.org/article/view/v033i06 |
| 13 | OpenBT | Source @ main (out-of-repo, so cited by file and function rather than by line): `brt.cpp` (`drawvec`), `brtmoves.cpp` (`pertcv`, every interior node; the 10% window; the window correction; the live `rot`); a later stretch of `brt.cpp` is one unclosed block comment, so its own `rot` copy is dead code; `pchgv = 0.1` in `misc/openbt.R` and `misc/openbt.py` | github.com/jcyannotty/OpenBT |
| 14 | He, Hahn, arXiv 2002.03375 | Relied on through `grow-from-root-default.md`'s verified record (sec 5, Table 4) | arXiv 2002.03375 |
| 15 | Chipman, George, McCulloch, JASA 93(443), 1998 | Read only as a verbatim block quotation inside #9 (JASA text paywalled). Cite accordingly. | (quoted in #9) |

**Composition block (section 4.1), sources and status.** Section 4.1 was
added late and its evidence base is deliberately weighted toward what could
be verified directly rather than toward citations.

| source | status |
|---|---|
| `docs/design/forest-ranef-interweaving.md` sec 0, 2, 5, 6, 9 | Read in full at `d3cb94b`. The 56.1 / 9.3 / 114.6 prototype table, the with-f/no-f attribution, the "no cheap ASIS/PX" structural argument, and section 9's authoritative corrections are all quoted from it directly. This is the load-bearing evidence for the hazard. |
| `LinearGaussianLeaf::logIntegratedLikelihoodForNode` [`LinearGaussianLeaf::logIntegratedLikelihoodForNode`](../../src/bartcore/model.hpp) (linear leaf integrated likelihood) | Read. Confirms the inner variant marginalizes leaf coefficients out of the structural score. |
| `dbartsSampler$setOffset` [`dbartsSampler$setOffset`](../../R/dbarts.R), `dbarts_sampler_setOffset` [`dbarts_sampler_setOffset`](../../inst/include/dbarts/dbarts.h), `dbartsModel`'s linear/gp node-prior checks [`dbartsModel`](../../R/model.R) | Read. Confirms the composition surface is public on both the R and C sides. |
| `inst/common/friedmanData.R` | Read. Confirms the probe DGP decomposes into one interaction plus three separable terms. |
| stan4bart's `src/init.cpp` (out-of-repo), `docs/design/walnuts.md` | Read in the live tree (0.0.14 installed). Confirms the 1:1 two-block Gibbs alternation and that no mixing diagnostic is reported. |
| Hahn, Carvalho, Puelz, He, Bayesian Analysis 13(1):163-182, 2018 | Full text (arXiv 1602.02176v3): the RIC definition, the competing-criteria mechanism, the closed-form bias (2.3), the reparameterization (2.5)-(2.6), and the appendix's extra alpha step added "to improve mixing". |
| Hahn, Murray, Carvalho, Bayesian Analysis 15(3):965-1056, 2020 | Full text (arXiv 1706.09523v4): RIC named for BART, the "single split in Z can stand in for many splits" mechanism, the propensity-covariate motivation, the mu/tau reparameterization, and "mu and tau alias one another". The PUBLISHED version with its discussion and rejoinder was 403-blocked and NOT checked - if a discussant raises mu/tau mixing, this survey did not see it. |
| Prado, Parnell, Murphy, McJames, O'Shea, Moral, CSP-BART, AOAS (arXiv 2108.07636v7) | Full text: the shared-covariate non-identifiability, the double-grow / double-prune paired kernel and its supporting constraints, the intercept-conflation trap, and appendix B's "identifiable but ... the individual components are not". Grepped for mixing diagnostics: **zero hits in 2592 lines**. |
| Bhandari, Bhatti, Chiu, Ji, arXiv 2605.20143v2 | Full text: "the flexible BART component can absorb variability that might otherwise be attributed to the linear predictor"; measured coefficient attenuation (their Table 6) that more data reduces but does not remove; orthogonality constraints proposed, not built. |
| Zeldow, Lo Re, Roy, AOAS 13(3):1989-2010, 2019 | Full text (arXiv 1806.04200): the "modeling a covariate in both ... sometimes led to bias and undercoverage" aside and the disjoint-covariate-sets remedy. Trace plots only, no ESS/ACF, no comparison arm. |
| Tan, Roy, Statistics in Medicine 38(25):5048-5069, 2019 | Full text (arXiv 1901.07504): section 4.1 is descriptive; **NOT-FOUND** on competition, identification, or block mixing. |
| Prado, Moral, Parnell, MOTR-BART, Statistics and Computing 31:20, 2021 | Full text (arXiv 2006.07493v5): "fewer trees are required", "the trees from MOTR-BART tend to be shallower than those from BART (10 trees)", and the parameter counts. Mixing only as an unquantified log-likelihood-convergence claim. |
| Linero, Yang, SoftBART, JRSS-B 80(5):1087-1110, 2018, and the SoftBart package vignette (arXiv 2210.16375) | Full text of both: **NOT-FOUND** in the paper (no semiparametric variant, no competition discussion). The vignette's section 4.3 partial linear model is a literal residual-swapping two-block Gibbs and asserts "the chain mixes well" from a trace plot, with no number and no comparison. |
| Yu, Meng, JCGS 20(3):531-570, 2011 (ASIS) | Full text (author PDF via Internet Archive): Theorem 1's bound `r_1&2 <= R_1,2 sqrt(r_1 r_2)`, the beauty-and-beast premise, and the scope - two augmentations of the same parameter linked by a map, which is not what a parametric block and a forest are. |
| Dorie, Perrett, Hill, Goodrich, Entropy 24(12):1782, 2022 (stan4bart) | Full text via EuropePMC: the Gibbs composition and offset exchange, and section 4.5's parameter-expansion framing with its own "More research will need to be performed to confirm this". Note the venue: Entropy, not Observational Studies. **No mixing diagnostic is reported for the composed sampler.** |
| The item-7 negative audit | Every paper above was searched full-text for effective sample size, autocorrelation, R-hat and Gelman-Rubin, case-corrected for the "regression" false positive. **No published measurement of a composed parametric-plus-BART sampler's mixing, against BART alone or otherwise, was found.** Two near-misses were blocked and NOT checked: the published BCF discussion/rejoinder, and Zhang et al., Statistics in Medicine 45 (2026) e70593, whose abstract describes cluster-level covariates and random intercepts inducing "near-collinearity and instability in selection" in a BART-plus-mixed-effects framework. That second one is the closest thing to a hit and should be retried from a network with journal access. |

Not independently re-fetched by this arc, and relied on only for claims
`docs/design/grow-from-root-default.md` already carries with its own
verification tags: Hill, Linero and Murray (Annu. Rev. Stat. Appl. 7);
Gelman and Rubin (Statist. Sci. 7); Carnegie (Stat. Sci. 34). The Bayesian
Analysis 11(3) discussion of Pratola is paywalled, was not read, and
nothing from it appears here.

## 11. Provenance

```
repo          /Users/vdorie/Repositories/dbarts, branch bartcore
code anchors  d3cb94b (all re-read at this tip by the adjudication pass;
              moves.hpp is unchanged since c637506, chain.hpp and tree.hpp
              are not, so earlier working-paper line numbers for those two
              files may be off by a few). Re-checked at 81df361, the tip
              when this document closed: d3cb94b..81df361 touches only
              TODO and docs/plans/, no source, so every anchor is live.
scope         research only - no source change, no commit, nothing scheduled
seeded by     TODO: tree-mixing-proposals (VD 2026-08-09)
working papers untracked {memo,critique,synthesis}.md
              (gitignored; synthesis.md carries the per-finding
              ADOPT/OVERTURN record and the evidence for each)
in-repo data  docs/design/grow-from-root-default.md sec 3, 4.4, 4.8, 4.9, 8
              docs/plans/archive/grow-from-root-default-study.md
              docs/plans/archive/bcf-sigma-residual.md sec 1-4
              docs/design/change-move-balance.md
              docs/design/grow-from-root.md sec 2, 3, 5
              docs/design/parallel-bart-frontier.md sec 3.1, 3.3
              docs/design/forest-ranef-interweaving.md sec 0, 2, 5, 6, 9
              docs/design/linear-leaves.md, docs/design/gp-leaves.md
out of repo   ~/Repositories/stan4bart: src/init.cpp (the two-block Gibbs
              loop, offset exchange), docs/design/walnuts.md (the
              parametric target, and the dropped diagnostics) - read-only,
              nothing written there
```

**Process note.** The survey pass that opened this arc ranked warming the
sampler during burn-in first, on the strength of five external sources.
An adversarial review found that none of the five evaluates that
construction, and an adjudication pass confirmed it by reading the
authors' released experiment code, which settles what the paper's prose
leaves ambiguous. Recorded because the failure mode is general and cheap
to repeat: **for an empirical claim about a schedule, a parameter, or a
dosage, the paper's prose is not the primary source; the code that
produced the number is.** Three of this document's most consequential
corrections - the temperature schedule spanning the kept draws, the
perturb throughput being 20-40x higher than proposed, and the perturb
window differing 8.5x between the paper and the shipped implementation -
came from reading code, not papers.

---

## 12. Addendum: forest specialization, whole-tree regrow, and eight orchestrator candidates (2026-08-10)

Status: ADDENDUM, adjudicated. No source touched, nothing scheduled. A
second research cycle ran on two VD directions commissioned 2026-08-09:
(a) whether individual trees can be made to target different parts of the
posterior, kept apart by penalties; (b) whether XBART's approximate
builder plus a "parametric, orthogonal trick" can generate proposals far
from the current trees. Same pipeline as this document's own: a research
memo, a blind refuting critique, then this adjudication, which re-opened
every derivation against the live tree at `ef7335d`, re-ran the numerics
independently, and re-fetched every citation it carries. The working
papers are untracked {memo,critique}.md files
(gitignored), so every load-bearing fact is carried here rather than
referenced.

The critique's verdict on the memo was STANDS WITH AMENDMENTS, five
blocking findings. This addendum adopts all five, on independent receipts
in four of them.

### 12.1 The four constructions the memo separated

VD's sketch ("multiple functions, each one of which describes a different
distinct region of the posterior, kept apart from each other [by]
penalties") is not one construction. Sum-of-functions is a decomposition
of the fitted FUNCTION, which is what BART already is; coverage of a
posterior DISTRIBUTION by an ensemble of whole explanations is a mixture
over explanations, a different object.

| # | construction | components indexed by | composition | exactness |
|---|---|---|---|---|
| 1 | input-space specialization (mixture of experts) | region of covariate space | sum or mixture over a gate | exact MCMC on a changed model |
| 2 | posterior-space specialization (repulsive particles) | posterior mode / explanation | mixture over particles | (2a) repulsion-as-prior exact; (2b) particle flow approximate at finite particle count |
| 3 | within-forest diversity | tree, within one forest | sum (unchanged) | exact MCMC on a changed prior |
| 4 | whole-tree regrow as a far-jump proposal | nothing - it is a kernel | n/a | exact, model and target unchanged |

The memo also proposed a sixth stickiness mode, **"mode F", apportionment
stickiness**: how the fitted signal is divided *among* trees at fixed
structure, as opposed to the five structural modes of section 3. B5 below
rules that this is a quantification of section 3.1's second clause, not a
sixth mode.

### 12.2 The five blocking findings, adjudicated

#### B1. The memo's leaf-prior no-go is REFUTED, and BART's own prior rewards specialization

**The memo's claim** (advertised as "the sharpest negative result in this
memo"): leaf values are exchangeable iid Gaussians across trees, so
conditional on the fitted function the prior is a minimum-norm penalty on
the apportionment, maximized at the equal split; hence "no exchangeable
Gaussian prior over tree contributions can reward specialization". That
was the memo's stated reason for ranking construction 3 near-last.

**ADOPTED: refuted.** Three defects, checked here.

1. *The exchangeability step is a non sequitur.* Exchangeability gives the
   conditional **mean** of `f_j` given the sum; the conclusion is about
   where the conditional puts its **mass**. `(1/3, 2/3)` and `(2/3, 1/3)`
   with equal weight is exchangeable, has mean `(1/2, 1/2)`, and puts zero
   mass at the equal split.
2. *The pointwise minimum-norm claim fails as soon as the structures
   differ.* The prior is iid Gaussian over **leaf values**, not over
   per-observation contributions. Take `m = 2` and two observations; tree 1
   a stump (value `a`, contributing `(a, a)`), tree 2 two leaves (`b`, `c`).
   Condition on the total `S = (s1, s2)`: the conditional prior is Gaussian
   on `{a + b = s1, a + c = s2}` with mode at the minimum of
   `a^2 + b^2 + c^2`, i.e. `a = (s1 + s2)/3`. At `S = (1, 0)` the mode
   apportions `1/3 : 2/3` at observation 1, not `1/2 : 1/2` - and the equal
   split is not merely improbable but **infeasible**, since a constant tree
   cannot equal `S/2` at both points. The tree whose partition resolves the
   signal is *given* the larger share, by the shipped leaf prior. The
   general reason: the memo's "pointwise norm" implicitly weights each leaf
   by its occupancy (`sum_leaves n_leaf mu^2`) whereas the prior is
   `sum_leaves mu^2`.
3. *On the memo's own block-additive example the shipped prior prefers the
   specialized arrangement by hundreds of nats.* "Specialized" (trees carry
   `f_A` or `f_B`) and "spread" (every tree carries `(f_A + f_B)/m`) cannot
   have the same structures: a tree carrying the sum of two functions on
   disjoint variable sets must partition the **product** grid. The factor
   that prices structures is the CGM tree prior. Recomputed here from
   [`CGMTreePrior`](../../src/bartcore/model.hpp) (`base = 0.95`, `power = 2`, root at depth 0),
   branching factors only, balanced trees:

   ```
   leaves  2   log p(T) =  -0.5936
   leaves  4   log p(T) =  -3.3727
   leaves  8   log p(T) = -12.4102
   leaves 16   log p(T) = -35.1314
   ```

   A 2-leaf specialized tree against a 2x2 product grid costs **2.78 nats
   per tree, 208 at `m = 75`**; against a 2x4 grid, **11.82 per tree, 886
   at `m = 75`**; a 4-leaf specialized tree against a 4x4 grid, **31.76 per
   tree, 2382 at `m = 75`**. (The critique labelled the 11.82 figure "4
   cells per block"; 4 cells per block is the 16-leaf grid, i.e. the 31.76
   row. Direction and magnitude are unaffected. These count branching
   factors only - the deeper tree also pays more split-rule factors, each
   below 1 - so every figure is a lower bound.) The leaf prior's
   factor-of-2 exponent penalty on apportionment at *fixed identical*
   structures is not in the same universe.

**What survives, and is worth recording.** Exactly one statement: *at
fixed, identical tree structures*, reapportioning a fitted function
unevenly across trees costs prior mass, with equality iff the per-component
tree budget is proportional to that component's magnitude
(`f_A/k = f_B/(m-k)`). That is a real design statement about
`blocks(trees.per.group = )`. It is not a no-go, because BART samples
structures, and specialization in BART arrives *as* a structural change -
smaller, differently-shaped trees - which is the direction the prior
rewards.

**Consequence for this record.** Nothing in sections 1-11 rested on the
no-go. What changes is downstream: a within-forest diversity prior is now
an open question rather than a pre-argued decline, and `blocks()` is
**prior-favoured** for a block-additive truth rather than merely
prior-neutral, which makes a null result on the `blocks()` arm *more*
informative, not less.

#### B2. `q = pi` at depth 1 is FALSE, and the correct object is a product over nodes

**The memo's claim** (advertised as "the most useful new argument in the
memo"): for a tree capped at depth 1 the builder's root candidate weights
are exactly `pi`'s factors, so `q = pi` identically and the regrow's
acceptance ratio is 1 - whence BART's shallow trees are the favourable
regime for an independence regrow.

**ADOPTED: false, on an independent derivation and an independent exact
enumeration.**

*The cap does not exist.* `CGMTreePrior::growthProbability`
([`CGMTreePrior::growthProbability`](../../src/bartcore/model.hpp)) returns 0 **only** when
`!tree.hasAnyAvailableVariable(...)`; otherwise `base/(1+depth)^power`,
strictly positive at every finite depth. The memo states this itself
elsewhere and then assumes its negation.

*The correct general formula*, re-derived here from [`growTreeFromRoot`](../../src/bartcore/grow.hpp)
directly. At a node `v` the builder's candidate set is `{no-split}` with
weight `(1-g_v) L(v)` and one entry per legal cut with weight
`g_v P(j) P(c) L(l) L(r)`, so `Z_v = L(v) [(1-g_v) + g_v B_v]` with
`B_v = sum_{j,c} P(j) P(c) L(l) L(r) / L(v)` the prior-averaged split
Bayes factor at `v`. Collecting `q(T)` over visited nodes and cancelling
against `pi(T) propto p(T) prod_leaves L(leaf)`:

```
pi(T) / q(T)  =  Z_root * prod_{w != root} [ (1 - g_w) + g_w B_w ]
```

the product running over every **non-root node** of `T`, internal and leaf.
`Z_root` does not depend on `T`. A depth-1 tree has two non-root nodes, so
`q = pi` there would require `g_child = 0` - the cap that does not exist.
The memo's diagnosis of the mechanism ("the divergence is the builder's
as-if-terminal scoring") is right; its claim about where the divergence
starts is wrong. It starts at the first split, and it scales with **node
count**, not with depth.

*Receipt, independent of the critique's.* Exact enumeration of the
single-predictor tree space under the shipped arithmetic - constant
Gaussian leaf exactly as [`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp), CGM prior exactly as
[`CGMTreePrior`](../../src/bartcore/model.hpp), builder weights exactly as [`growTreeFromRoot`](../../src/bartcore/grow.hpp),
`nodeScale = 0.5`, `k = 2`, `m = 1`, `n ~ 420` (an even split over the
codes), `q` asserted to normalize to 1 and the closed form above asserted
to reproduce `log pi - log q` for every enumerated tree to 1e-9:

```
cuts  signal  sigma  trees  accept   log-wt spread   spread over depth<=1
  4    0.0     1.0      51   0.964        0.26              0.23
  4    0.0     0.3      51   0.913        0.50              0.42
  4    0.2     0.3      51   0.919       35.47             14.36
  4    0.5     1.0      51   0.842       15.77              6.32
  4    1.0     0.3      51   0.941     1023.72            457.15
  2    0.0     1.0       5   0.966        0.89              0.89
  6    0.0     1.0     731   0.965        0.30              0.20
  7    0.0     1.0    2950   0.966        0.37              0.18
  7    0.3     1.0    2950   0.875        2.23              1.25
```

The depth<=1 spread is never zero, so `q != pi` there in every cell. Two
things the table also settles, both used below: acceptance does **not**
decay as the tree space grows 600-fold at fixed signal (0.966, 0.964,
0.965, 0.966 across 5, 51, 731 and 2950 trees), and the weight spread
grows by four orders of magnitude with signal while acceptance does not
fall with it. The critique's own enumeration, on a different
data-generating process, gave acceptance 0.384-1.000 and a depth<=1 spread
of 0.99-1323 nats. Two independent enumerations agree on the conclusions,
not on the digits.

**What survives, and it is worth keeping.** The log importance weight is a
*sum over non-root nodes*, so a prior that keeps trees small keeps the
number of terms small. That is a derived version of "BART's shallow trees
are the favourable regime", and it gives Stage R0 a **per-node** quantity
to log rather than a whole-tree acceptance rate.

**The SNR direction needs re-derivation, not deletion.** The memo argued
high SNR is the *best* regime for the regrow. The weight spread says the
opposite - it explodes with signal, 0.26 nats at zero signal against 1024
at signal 1 in the table above - yet realized acceptance stays high,
because `pi` and `q` concentrate on the same tree. Both facts are real.
Acceptance governs the move's throughput; the spread governs the
**variance** of anything built on it. Log both.

#### B3. The reallocation census cannot do the two jobs it was assigned

The memo's #1 recommendation was a "reallocation census" measuring the
per-tree apportionment autocorrelation time. **ADOPTED with a replacement
statistic**; the census is still worth running, but not as designed.

- *(a) The "discriminating" sigma cell does not discriminate.* The memo's
  own timescale for mode F is `n_leaf tau^2 / s^2`, which **grows** as
  sigma falls; mode B (sec 3.2) also gets worse as sigma falls. Same
  direction, so the sigma cell is the least discriminating in the design,
  not the most - as the memo says itself two pages earlier ("If both are
  real, low noise is doubly bad").
- *(b) The `blocks()` arm cannot move the statistic the kill criterion
  reads.* Verified by reading the code: `forest.leaf.scale =
  resolvedNodeScale(options.nodeScale, options.priorScale) / sqrt(forest.numTrees)` ([`Chain::Chain`](../../src/bartcore/chain.hpp)) uses the
  **total** tree count, and `installBlockMasks` ([`installBlockMasks`](../../src/bartcore/chain.hpp),
  called at [`Chain::Chain`](../../src/bartcore/chain.hpp)) installs per-tree column masks and nothing else - no
  per-group rescaling anywhere in the function. So the per-leaf prior sd
  `tau` is identical in both arms and the within-block apportionment
  timescale is the same number. What `blocks()` removes is not the
  timescale but the **dimension** of the slow subspace, `m - 1` free
  apportionment directions down to `m - G`: at `m = 75` with 4 groups, 3
  directions out of 74. Pre-registering a kill on a test the theory
  predicts will fail is not a falsifier.
- *(c) The arm contrast is confounded even if (b) is wrong.* In the
  unrestricted arm a tree's fitted vector can wander across the whole
  function space; under a mask it cannot. The restricted arm's series has
  strictly smaller support by construction, so any reduction it shows is
  partly definitional. A fair contrast needs a statistic invariant to the
  restriction - e.g. apportionment *within* a fixed variable block,
  measured in both arms.
- *(d) The kill ratio's denominator is pathological in the cells of
  interest.* "Within 2x the total-fit IACT" can fire in the low-noise cell
  because the total-fit IACT is itself huge there (this record measures 6.5
  to 949.4 across ten cells, sec 3.5); and in the memo's own stumps
  caricature the apportionment subspace has posterior **equal to** the
  prior, so a large ratio is analytically guaranteed and measuring it
  confirms nothing.

**Replacement statistic: measure the coupling, not the marginal.** Tree
`j`'s structural question is scored against
`treeY = y - sum_{k != j} f_k`, a running residual rolled tree by tree
inside the sweep ([`rollTreeResidual`](../../src/bartcore/chain.hpp)), so it depends on the other trees
**only through their total**. Apportionment can therefore reach structure
only through `f_j`'s own drift. The pre-registerable claim is: *tree `j`'s
structural acceptance pattern is autocorrelated at the `f_j` timescale, and
that autocorrelation is what apportionment stickiness costs.* Same
instrumentation, falsifiable in both directions.

**Three further corrections to the census design**, all adopted.

- **No matching or alignment step is needed.** The sweep updates tree `t`
  in place ([`Chain::run`](../../src/bartcore/chain.hpp)) and nothing in `chain.hpp` shuffles, permutes
  or relabels trees (verified by search). The *posterior* is
  label-exchangeable; the *chain* never exercises the symmetry.
- **`getTrees` is necessary but not sufficient.** It returns "a data.frame
  containing the internal state of the trees" ([`dbartsSampler$getTrees`](../../R/dbarts.R)) - flat
  node structure with leaf values, decoded categorical directions and
  missing routes - not per-tree fitted vectors. The census must walk trees
  in R itself (the package walks trees in R only in
  `plotTree`/`getTreeDepthAndSize`, neither of which evaluates one) or take
  a small C-side per-tree readout. That is real code, not zero.
- **Vary `k`, and drop the `m = 1` control.** At `m = 1` the per-tree
  fitted vector *is* the total fit, so the two autocorrelation times
  coincide by identity - that tests the plumbing. And raising `m` at fixed
  `n` raises per-leaf occupancy, so the two `m` effects partly cancel in
  the prediction; `k` enters **squared** (B5), is a pure prior knob, and
  does not move the tree geometry.

#### B4. The support enumeration is incomplete; the refusal predicate must not be a hand list

Construction 4's exactness rests on refusing the move when the incumbent is
outside the builder's support. The abstract argument is sound (12.3 below),
but the memo's four-item support list is incomplete.

**A fifth gap: monotone / `ParamScoring` leaves.** Verified here:
`growForestFromRoot` refuses only `hasVectorParams || hasFunctionParams`
([`growForestFromRoot`](../../src/bartcore/chain.hpp)), so the **monotone constant leaf is in the
builder's scope**. `MonotoneConstantGaussianLeaf::logIntegratedLikelihood`
delegates to the *unconstrained* `ConstantGaussianLeaf` and is documented
"never on the constrained hot path" ([`MonotoneConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp)). Meanwhile the
structural target under that leaf is not a leaf-marginalized posterior at
all: `logLikelihoodForBranch` dispatches `ParamScoringLeafModel` to
`leaf.logLikelihoodForBranchWithParams(...)`, reading frozen neighbour
parameters ([`logLikelihoodForBranch`](../../src/bartcore/moves.hpp)). So for a monotone forest `q` is computable
but `pi` is **not** `p(T) prod_leaves L(leaf)`, and "every term in the
acceptance ratio has a shipped implementation" is false there. Monotone
forests must be scoped out of any regrow v1 explicitly.

**The predicate itself.** A hand-written list (categorical split | zero-
growth node | missing-only side | non-constant leaf | ...) will rot. The
equivalent, always-correct predicate is a property of the replay the move
needs anyway: **refuse iff the reverse replay returns `-inf`**. It is
decidable from the incumbent alone because the builder's support is
residual-independent - occupancy depends on member counts and availability
on the cut grid and masks, never on `y` ([`scanOrdinalCuts`](../../src/bartcore/scan.hpp),
[`variableAvailable`](../../src/bartcore/tree.hpp)) - it is exactly the support indicator, and it costs
nothing extra because the replay runs regardless.

#### B5. Mode F is not a sixth mode, and its timescale is off by `k^2`

**ADOPTED on both halves. This addendum amends section 3.1.**

*Novelty.* Section 3.1's own first paragraph already names the dynamic:
"the same ensemble fit arises from permuting tree labels, **or from
splitting one main effect across two trees instead of one**". And section
4.1 makes the *cross-block* instance of the same ridge its central hazard,
measured in this house at 6x. Mode F is a **quantification** of section
3.1's second clause, specialised to leaf-value space at fixed structure.
**Section 3.1 is amended to say so; no sixth letter is added to the
taxonomy**, and the memo's claim that the dynamic is absent from this
record is withdrawn. Section 3.2's opposing datum - "leaf values converge
in a handful of sweeps; what remains is a partition-shape misfit" - stands
and must be engaged by anything that measures it.

*Timescale.* The shipped leaf prior sd is `scale/k`, not `scale`:
`priorPrecision = (k/scale)^2` ([`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp)) and
`drawFromPrior = (scale/k) * z` ([`ConstantGaussianLeaf::drawFromPrior`](../../src/bartcore/model.hpp)). So

```
tau = nodeScale / (k sqrt(m)),   timescale ~ n_leaf nodeScale^2 / (m k^2 s^2)
```

a factor `k^2 = 4` smaller at the default `k = 2`. With `nodeScale = 0.5`
([`ParsedModel`](../../src/R_interface_bartcore.cpp)), `tau = 0.0289` at `m = 75`, not
`0.0577`. `k` is the shipped knob that enters the prediction **squared**
while `n`, `m` and `sigma` enter linearly, and it is itself sampled when
`updateK` is on ([`Chain::run`](../../src/bartcore/chain.hpp)), which the caricature assumes
fixed.

### 12.3 Substantive advisories, adjudicated

- **The variable-ban regrow is reversible - HOLDS, with two conditions the
  memo omits.** The construction: grow the proposal against the unmodified
  residual but with the incumbent's realized split variables banned, so
  `used(T)` and `used(T')` are disjoint and `T` is in the reverse support.
  Checked mechanically and it does not fail: every internal node of `T`
  splits on a variable in `used(T)`, which the reverse ban leaves alone and
  which is available at that node by construction, so `growthProbability`
  is positive at every internal node of `T` under the reverse ban
  ([`CGMTreePrior::growthProbability`](../../src/bartcore/model.hpp)); and banning changes `numAvailable` and hence
  `P(var)`, but both directions compute their own normalizer, which is all
  MH needs ([`growTreeFromRoot`](../../src/bartcore/grow.hpp) mirrors [`CGMTreePrior`](../../src/bartcore/model.hpp)). **Condition
  one**: the ban does not relieve B4's gaps, it stacks on them, so the
  refusal predicate must be evaluated **under the reverse ban**, not on the
  incumbent in isolation. **Condition two**: `q` depends on `used(T)`, so
  the kernel is *not* an independence sampler and must not inherit
  independence-MH analysis. Scope the memo already states and this pass
  confirms: the move can only reach disjoint-variable representations, so
  at `p = 10` with a tree using 4 variables it explores a much-reduced
  space, and at high `p` with sparse trees it is nearly the plain regrow.
- **Mixture-of-kernels invariance needs no citation.** If `K_i pi = pi` for
  each `i` and `K = sum_i a_i K_i` with `a_i >= 0` summing to 1, then
  `K pi = pi` by linearity. Reversibility of the refusing regrow is equally
  short: on the support set the MH ratio gives detailed balance pairwise, a
  proposal never leaves the support, off it the kernel is the identity, and
  there is zero flow across the boundary in both directions. The memo's
  unverified Tierney / Roberts-Rosenthal citation should be dropped rather
  than chased; the genuinely unproven part of that argument was the support
  characterization, which is B4.
- **Three corrections to the builder's `q`, all verified here.**
  1. The dropped `sum w z^2` term cancels only when the candidate classes
     share a member set. On a **missing-capable column they do not**: the
     no-split candidate is scored from the node's cached statistic over
     *all* members ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)) while `scanOrdinalCuts` skips
     `naCode` members outright ([`scanOrdinalCuts`](../../src/bartcore/scan.hpp)). This does not break
     exactness - `q` is whatever the builder's realized normalized weights
     are, and that is what an accumulator records - but it does break the
     memo's stated reason, and it is a structural `q`/`pi` mismatch rather
     than a data-driven one.
     SUPERSEDED PREMISE: the scan no longer skips them - [`scanOrdinalCuts`](../../src/bartcore/scan.hpp)
     accumulates the `naCode` rows into a missing bin that every candidate
     adds to one of its two children, and [`scanOrdinalCuts`](../../src/bartcore/scan.hpp) states that
     the scan's scores therefore agree with the leaf statistics
     `tree.birth` caches, missing rows included (TODO item
     `ordinal-scan-missing-rows`, DISCHARGED), so the structural mismatch
     concluded from the old premise no longer holds.
  2. **Missing-capable columns are halved twice.** `logCut` already
     subtracts `log 2` for `data.hasMissing[j]` ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)), matching
     `ruleForVariableLogProbability`'s `+log 2` for *one* rule
     ([`ruleForVariableLogProbability`](../../src/bartcore/model.hpp)); the builder then draws the direction coin
     separately ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)). Since the candidate stands for the
     *pair* of direction-rules, its weight carries one rule's prior mass,
     so `q` under-weights splits on missing-capable columns by 2x per
     split. **This is already known and scheduled**: it is exactly the
     pre-registered two-arm falsifier in
     `docs/plans/archive/grow-from-root-categorical-scan.md` sec S0, with an OPEN
     VD FORK on whether to change the shipped weight. What a regrow adds is
     stakes: for a warm start it is a start-quality bias, for a regrow it
     becomes a systematic term in the importance weight.
     SUPERSEDED PREMISE: that `log 2` is gated on `routesMissing`
     ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)) - whether the scan emitted both directions for THIS
     node's members ([`growTreeFromRoot`](../../src/bartcore/grow.hpp)), each direction then its own candidate -
     and not on `data.hasMissing[j]`; the direction coin ([`growTreeFromRoot`](../../src/bartcore/grow.hpp))
     fires only where a candidate did NOT already name its direction. The
     double-halving concluded from the old premise no longer holds.
  3. The reverse replay must scan **every node of the incumbent, leaves
     included** - a leaf contributes `(1-g) L(u) / Z_u` and `Z_u` needs the
     full candidate assembly there. Cost is `2 x (nodes)` candidate
     assemblies, not `2 x (levels)`.
- **The categorical support gap is confirmed, and the `P(var)` accounting
  is right.** [`growTreeFromRoot`](../../src/bartcore/grow.hpp) skips a subset-splitting column outright
  while `numAvailable` counts categoricals even though they generate no
  candidates - which is *correct* for matching [`CGMTreePrior`](../../src/bartcore/model.hpp)'s
  `-log(numAvailable)`, since both sides count the same set and the missing
  mass is absorbed by `Z_node`. Worth recording because it is the one place
  the builder's `P(var)` and the prior's agree exactly and it is not
  obvious from either file alone.
  SUPERSEDED PREMISE: [`growTreeFromRoot`](../../src/bartcore/grow.hpp) now ENTERS the categorical branch and
  emits one candidate per admissible partition
  (`scanCategoricalPartitions`, [`growTreeFromRoot`](../../src/bartcore/grow.hpp), [`scanCategoricalPartitions`](../../src/bartcore/scan.hpp)), so
  categoricals do generate candidates and the support gap concluded from
  the old premise is closed.
- **Construction 2(a)'s exactness needs a DART scope condition.** A
  cross-tree repulsion prior has a normalizing constant over the joint tree
  configuration space; it cancels in the MH ratio only while everything it
  depends on is fixed. dbarts ships DART, which resamples the split
  probabilities every sweep from realized split counts
  ([`Chain::run`](../../src/bartcore/chain.hpp)), and that constant depends on those
  probabilities. So fixed-mask `blocks()`-style constraints stay exact, but
  **DART plus any cross-tree split-usage repulsion is doubly intractable in
  the `s` draw** - which lands on the memo's own "cheap and exact" verdict
  for the cheapest repulsion variant, since that variant's cheapness came
  from reusing DART's counts. The `k` hyperprior
  ([`Chain::run`](../../src/bartcore/chain.hpp)) is unaffected; only variable-selection parameters
  enter the constant.
- **Census cost is understated.** To estimate an autocorrelation time of
  order `10^2-10^3` needs chains far longer than `10^3` draws, with
  `keeptrees = TRUE` storing per-sample forests for every cell; "hours of
  compute" holds only for the small cells.

### 12.4 Citation corrections carried into this record

Every row was fetched and read by this pass. Where a locator or referent in
an earlier paper was wrong, the correction is here and the earlier form
must not be repeated.

| # | source | correction | receipt |
|---|---|---|---|
| 1 | Lakshminarayanan, Roy, Teh, AISTATS 2015 | Both quotes sit in **sec 3.3 "PG sampler for BART"**, not sec 4; sec 4 is the experimental evaluation. The introduction (sec 1) additionally carries the sentence that adjudicates sec 5.4's erratum. | primary PDF, `pdftotext`: sec 3.3 "The conditional-SMC algorithm is an MH kernel with pi(T_j) as its stationary distribution"; sec 1 "proposing complete trees from the tree prior, however these moves would be rejected, leading to slow mixing... because those non-local moves have high posterior probability" |
| 2 | He, Hahn, arXiv 2002.03375 | The "not a proper full conditional" sentence is in **sec 3.3 "Prediction"** and "this estimator" is their point-wise posterior-mean *predictor*, not the sampler - so it must be quoted with a scoping clause. Two adjacent facts make the point better: sec 6.2 opens "This section proves that a **slightly modified** version of GrowFromRoot generates draws from a Markov chain with a stationary distribution. The slight modification is that all leaf parameters are drawn jointly...", and Theorem 2 states only "a finite-state Markov chain with stationary distribution" - the abstract's "**unique**" does not appear in the theorem. | primary PDF, `pdftotext` |
| 3 | Chipman, George, McCulloch, BART, sec 3.1 | The famous fragment is half a sentence. In full: "Although mixing does not appear to be an issue, **the recently proposed modifications of Blanchard (2004) and Wu, Tjelmeland and West (2007) might well provide additional benefits.**" BART's own authors point at the non-local-move literature in the same breath, which is evidence *for* sections 5.2 and 5.4, not against them. | arXiv 0806.3286v2 PDF, `pdftotext` |
| 4 | Du, Linero, DP-Forests, AISTATS 2019, PMLR 89:108-117 | **A published within-BART-ensemble diversity prior exists**, so any claim of the form "no published within-ensemble diversity prior in BART" is false. What survives is the narrower "no published *repulsion between trees* in a BART ensemble". Its cost is ~2x and it reports no mixing diagnostic at all. | published PMLR PDF, `pdftotext`: sec 3 "we specify a prior which clusters the trees into non-overlapping groups such that each cluster constructs splits using different subsets of the predictors"; sec 4 "SBART and DP-Forests took 118 seconds and 241 seconds respectively to obtain 40,000 samples from the posterior. By comparison, iRF took 279 second, HL-LS took 91 seconds, and additive groves took 4966 seconds"; **zero** occurrences of mixing, effective sample, autocorrelation, R-hat, Gelman, convergence, Markov chain, burn-in or MCMC in the whole PDF |
| 5 | Thakkar et al., Quantum Machine Intelligence 6 (2024), arXiv 2306.12965 | A determinantal point process HAS been applied to a tree ensemble, so "no DPP work on tree ensembles" is false as written. It is a DPP over **data**, not a repulsive prior over trees, so the argument it was cited against is untouched. | arXiv PDF, `pdftotext`: "we propose an extension of the Random Forest, called the DPP-Random Forest (DPP-RF), which utilizes Determinantal Point Processes (DPPs) instead of uniform sampling to subsample rows and columns for individual decision trees" |
| 6 | Neal, "Sampling from multimodal distributions using tempered transitions", Statistics and Computing 6:353-366, 1996 | Bibliographic record verified from the reference list of a citing paper; **the primary was not fetched**. The citing paper's abstract is the load-bearing caveat: "Unfortunately the improved movement between modes comes at a high computational cost with a low acceptance rate of expensive proposals", and its body records that "Neal (1996) demonstrates that the algorithm satisfies detailed balance with respect to the target". | Behrens, Friel, Hurn, "Tuning Tempered Transitions", arXiv 1010.0842 PDF, `pdftotext`, abstract and reference [14] |
| 7 | Liu, Liang, Wong, "The multiple-try method and local optimization in Metropolis sampling", JASA 95(449):121-134, 2000 | Bibliographic record only; **body not read**. Named for R3(b), which is not scheduled. | publisher record and author-hosted PDF listing |

Incidental, and it belongs beside section 5.3: the PG-BART paper **does**
report effective sample size and ESS per second, and PG *loses* on ESS/s in
the shallow-tree regime BART's prior produces (their Table 2, Hypercube-D:
at `D = 2`, CGM 157.67 against PG 7.69; at `D = 3`, 93.01 against 11.03),
winning by an order of magnitude only at `D >= 4`. That is the same
shallow-tree caveat section 5.3 already records, now with numbers.

**Not verified by this pass, and therefore not carried:** the memo's whole
PASS-VERIFIED tier for constructions 1-3 (the treed-model and treed-GP
ancestors, SoftBART, the SVGD and repulsive-mixture line, the label-
switching literature, Wood et al., Breiman). Those claims stay in the
gitignored memo. Mengersen and Tweedie remains unretrieved by anyone in
this program; after the reversibility ruling nothing here depends on it.

### 12.5 The eight orchestrator candidates, item by item

None carried authority. Each is adjudicated on the merits.

**R1. Ridge-traversal move for composition - SET ASIDE (subsumed), with a
corollary worth keeping.** The proposal: shift smooth mass between a
parametric block and the trees while holding the fitted function fixed
(`theta -> theta + delta`, leaves absorb `-delta`), so the likelihood
cancels and acceptance is prior-ratio only. The open question it names -
"work out what IS exactly absorbable" - has a short linear-algebra answer.
Holding structures fixed, the forest can absorb `-Z delta` exactly iff
`Z delta` lies in the span of **all** leaf indicators over all `m` trees,
i.e. in the space of functions representable as a sum of the current
trees' piecewise constants. The constant function is always in that span
(a tree's leaf indicators sum to 1), and a continuous column never is, so
for a design whose non-intercept columns are continuous the **only**
exactly absorbable direction is the intercept - shift `c`, subtract `c/m`
from every leaf of every tree. Once absorption is inexact the likelihood
stops cancelling and the prior-ratio-only tractability - the whole point of
the move - is gone. Note also that the span is state-dependent: which
directions are absorbable changes as the trees move.

That is decisive for the case that motivated it. The 6x hazard section 4.1
records is a forest-versus-**group-intercept** ridge, and a per-group shift
is a step function on group membership, which lies in the leaf-indicator
span only if some tree splits on the group indicator - and in dbarts'
grouped model the group is not a predictor column at all. So the move does
not exist for the one composition this house has measured. This is
`forest-ranef-interweaving.md`'s recorded "the forest exposes no per-group
scalar" argument, now as an exact statement about spans rather than an
intuition, and it is why R1 does not reopen that NO-GO. One genuine edge,
recorded rather than built: if a user puts group dummies in `x` *and* a
tree splits on them, the shift is absorbable - i.e. the move exists exactly
when the forest is already rendering the parametric direction.

**R2. Identified specialization (multi-resolution forest; block-structured
DART) - ADMITTED, gated behind the `blocks()` arm.** The premise is sound
and is the one B1 leaves standing: components made non-exchangeable *by
construction* carry no label-switching problem, which is the same
"labelled" column `blocks()` already occupies. Two corrections to the
stated mechanism, and one cost fact.

- *A depth prior does not confine structure scale.* A depth-2 tree is
  already a two-way interaction, and a shallow tree's cells can be
  arbitrarily small if its cuts sit at the edge of the grid. Depth caps the
  number of cells, not their size, so "shallow block = coarse/smooth niche"
  is not what the prior buys. If the intent is a smooth-versus-local split,
  the identifying constraint is on *variables* (which `blocks()` and
  `interactions()` already express) or on interaction *order* (which
  `interaction-constraints.md` already ships), not on depth.
- *Per-block depth priors are not a knob.* `CGMTreePrior` is a per-**forest**
  member ([`Forest`](../../src/bartcore/combiner.hpp)), so per-block `base`/`power` needs a per-tree
  prior indirection through every scoring path.
- *(b) is DP-Forests with fixed labels.* Per-block variable-inclusion
  priors are exactly Du and Linero's construction with user-supplied rather
  than sampled clusters (12.4 row 4), i.e. published, buildable, ~2x, and
  reporting no mixing benefit. It is also the soft form of what `blocks()`
  does hard.

Both variants are downstream of the question the `blocks()` arm asks, so
they are gated on it rather than scheduled beside it. Documented absence:
no published BART variant with per-tree-block depth priors was found by
this pass (one web search over BART plus varying/multi-resolution depth
priors; nearest hits are SoftBART's smoothness adaptation and DP-Forests).

**R3(a). Partial regrow (rebuild only below a chosen internal node) -
ADMITTED, and promoted above the full regrow.** Better shaped on three
counts, all following from B2's closed form. The importance weight is a
product over the **rebuilt subtree's** non-root nodes only, so its spread
is smaller by exactly the terms outside. The support requirement is
confined to the subtree, so a tree with a categorical split *above* the
rebuild point is still eligible - which makes it strictly less blocked by
the `grow-from-root-categorical-scan` dependency than the full regrow is.
And the rebuild depth is a genuine jump-size dial interpolating the change
move and the full regrow, making it the only far-jump candidate in this
document with a tunable step size other than perturb (sec 4.2) - the
property that made perturb first. It does not dodge the largest missing
piece: the correction needs the node-selection probability in both
directions (the incumbent's and the proposal's internal-node counts
differ), and rollback still needs the subtree save/restore primitive
rotation needs (sec 4.3), since `SubtreeSnapshot` cannot undo a shape
change.

**R3(b). Multiple-try - SET ASIDE pending Stage R0, not scheduled.** Its
value is entirely conditional on acceptance being the binding constraint.
Both enumerations put single-try acceptance far from the floor (0.84-0.97
here, 0.53-0.76 in the critique's cells), in `m = 1` caricatures. If R0
reproduces anything of that order at ship scale, a K-fold proposal cost
buys little. Gate it on R0's number.

**A1. Marginal-sigma structural scoring - SET ASIDE.** This was the
orchestrator's top pick, and it fails on two independent grounds, either
sufficient.

1. *Its premise does not hold for dbarts' prior.* The Student-t marginal it
   names requires the leaf prior to scale with sigma - the classical
   conjugate normal-inverse-gamma CART setup. dbarts' does not:
   `priorPrecision = (k/scale)^2` with `scale = nodeScale/sqrt(m)`
   ([`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp), [`Chain::Chain`](../../src/bartcore/chain.hpp)), fixed and sigma-free. That is
   BART's design, not an oversight. Integrating sigma against a
   fixed-variance leaf prior gives no closed form in the same sufficient
   statistics.
2. *Even granting the algebra, the marginal cannot move the number it
   targets.* Under a scaled-inverse-chi-squared sigma prior with
   `nu = sigmaDf` (3 by default, [`ModelParameters`](../../src/bartcore/chain.hpp)), marginalizing replaces
   the exponent `dSS / (2 s^2)` with
   `((nu + n)/2) log(1 + dSS / (nu lambda + SS))`, and the two agree to
   first order whenever `dSS` is a small fraction of the total residual sum
   of squares. Arithmetic in the cell that matters: at `n = 5000`,
   `sigma^2 = 0.1`, a proposal at `dLogL = -50` has `dSS = 10` against a
   total near 500, so the marginal softens it to -49.5 - **1% relief**
   against a freeze this document prices at 10-280x out, and roughly 400x
   in the causal-forest tail (sec 5.1 point 3). The softening becomes real
   only when one tree's structural move changes the whole residual sum of
   squares by an `O(1)` fraction, which does not happen in a 75-tree
   ensemble.

Recorded rather than dropped, because the mechanism is real in single-tree
Bayesian CART and inert in BART for a reason worth knowing: BART's leaf
prior is deliberately sigma-free.

**A2. Tempered transitions (Neal) - ADMITTED into the section 4.4
temperature family, ranked below it, gated on the census fork.** Exact by
construction, and the within-chain member of the family whose
between-chain member is section 4.4. Nothing here reverses section 5.1,
which declined a *third* thing - a burn-only anneal with no stationarity
claim. Two additions from this pass. The cost is the published weakness
rather than a guess (12.4 row 6: "high computational cost with a low
acceptance rate of expensive proposals"), which is why it ranks below 4.4,
whose kept draws are exact by a simpler argument and which has verified
empirical support on tree posteriors while tempered transitions has none
this pass could find. And the *likelihood* half has an unexpectedly cheap
route in this engine: raising a Gaussian likelihood to a power `beta` is
exactly scaling every observation weight by `beta`, and the leaf marginal
reads weights only through `(sum w, sum w z)` ([`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp)), so the
intermediate distributions' leaf draws and structural scores need no new
leaf math. The *prior* half does - and section 5.1 point 1 establishes that
tempering the likelihood alone points the wrong way - so a valid
construction must temper the CGM factors too, and that is where its cost
lives.

**A3. Subtree crossover / graft within one forest - SET ASIDE, on a
mechanism error and on target.** "Computable acceptance (additive forest)"
does not hold in the form this engine scores moves in. A joint proposal on
trees `i` and `j` must be scored against `y` minus every *other* tree,
under the sum `f_i + f_j`, and the leaf-marginalized likelihood of a sum of
two trees is not the product of their per-tree marginals: both trees' leaf
values are integrated jointly against overlapping design columns - the
common refinement of the two partitions - so the per-branch factorized
score at [`logLikelihoodForBranch`](../../src/bartcore/moves.hpp) cannot serve and the cost is a joint solve
rather than a sum of per-leaf scalars. A leaf-*conditional* variant is
computable directly from fits, but it gives up the marginalization that is
why BART's structural moves accept at all. Independently: its target is
mode F, which per B3 and B5 is diffusion along a direction the posterior
does not constrain and that no label-invariant functional reads.

**A4. Importance-weighted warm starts - SET ASIDE, and the item's own
doubt is wrong.** The doubt ("this is a Markov-chain init bias, not a
simple importance bias") does not hold. If the initial state carries weight
`w = pi/q` and every subsequent kernel is `pi`-invariant, the weighted
estimator is unbiased at every later sweep:
`E[w f(X_t)] = int q (pi/q) K^t f = int pi K^t f = E_pi[f]`. The
construction is valid in principle. It dies on weight variance, and B2's
closed form prices it: a forest's weight is a **product over 75 trees** of
per-tree weights whose log spread is a sum over that tree's non-root nodes,
and the per-tree spread is 0.18-457 nats even in the `m = 1` caricature
above. Any per-tree spread of order 1 nat gives a forest log-weight spread
of order 75 nats, so the effective sample size of a weighted ensemble
collapses to one member at any feasible chain count. There is also an
output-contract precedent: importance-tempered draws are what `tgp` ships,
and it ships them off by default (sec 5.1 point 2). Finally, the study this
reopens was killed on a per-cell plateau posterior-mean RMSE cost measured
at 1000+ draws, and its reopen clause is not "weight the init".

**A5. Mode atlas - SET ASIDE, subsumed.** Its proposal densities are
computable only for atlas entries inside the builder's support, so it
inherits B4's refusal discipline wholesale. Its real cost is validity, not
bookkeeping: an atlas refreshed from the chain's own history makes the
kernel **adaptive**, and an adaptive kernel is not `pi`-invariant under the
standard MH argument - it needs diminishing adaptation and containment, or
a frozen atlas. A frozen atlas built before the chain runs is stale by
construction, because every entry was fit against a residual the other 74
trees have since moved. What remains once both are handled is R3(b) over a
cached candidate set.

### 12.6 Ranked disposition

Three of the eight orchestrator candidates are admitted and five are set
aside; all three admissions are gated on a measurement, none is scheduled.
Constructions 1 (soft or treed gate) and 2 stay research programs,
unchanged: section 4.4 continues to dominate construction 2 on validity and
evidence, and construction 1's value remains bounded above by the
`blocks()` result. Nothing in sections 4-7 is re-ranked - section 5.4 keeps
its rank and its action, and loses only its recorded reason - and this
addendum schedules no code.

| rank | item | what it is | cost | gate |
|---|---|---|---|---|
| 1 | **Stage R0** - generator-only regrow census | build a proposal, compute `log alpha`, discard it; no state change, no draw-law change | days; rides sec 6.1's Stage-0 instrumentation as extra logged columns | none (pilot) |
| 2 | **reallocation census, re-specified** | the coupling statistic of B3, not the apportionment IACT | ~300-400 lines of R plus a tree walk or a small C-side readout | none (pilot) |
| 3 | **the `blocks()` arm** | hard, free, exact variable specialization as the premise test for constructions 1 and 3 | zero engine work | restated kill (B3(c)) |
| 4 | R3(a) partial-regrow dial | rebuild below a chosen node; the jump-size dial | M plus a gate arm plus the shared subtree save/restore | conditional on R0 |
| 5 | A2 tempered transitions | within-chain member of the sec 4.4 family | M-L | live only if the census says rejections are close calls |
| 6 | R2 identified specialization | fixed-label DP-Forests / per-block priors | M-L, plus a per-tree prior indirection for the depth variant | conditional on (3) |
| - | set aside | R1, R3(b), A1, A3, A4, A5 | - | receipts in 12.5 |

**Stage R0, specified.** At each tree each sweep, build a candidate `T'`
from a separate RNG stream, compute
`log alpha = log pi(T') + log q(T) - log pi(T) - log q(T')`, log it, and
discard `T'`. Log, per proposal: the **per-node** terms
`log[(1-g_w) + g_w B_w]` (B2's decomposition - the quantity that governs
both the acceptance and the variance, and the reason to log per node rather
than per tree); the realized log acceptance; whether the incumbent was in
support and which reverse-replay node returned `-inf` (B4); the structural
distance between `T` and `T'` (split count, variable-set Jaccard, partition
disagreement); and build-plus-replay wall cost. Run the R3(a) rebuild-depth
dial and the variable-ban variant as logged variants of the same generator.
Scope monotone forests out (B4). Cells as in sec 6.1, plus one
all-categorical cell where the move must refuse 100% of the time. No kill
criterion: this is a pilot, and it freezes the thresholds anything
downstream uses.

**The re-specified reallocation census.** Statistic: the coupling (B3), not
the marginal apportionment autocorrelation time. Varied knob: `k`, which
enters squared. Dropped: the sigma cell as "discriminating", the `m = 1`
control, and the premise kill, which is a foregone answer in the informative
direction and can fire spuriously in the low-noise cell. Added: an R-side
tree walk or a small C-side per-tree readout, because `getTrees` does not
return fitted vectors. Its cheapest form is as a rider on the pre-registered
composition probe (`docs/plans/archive/composition-mixing-probe.md`, RE-REGISTERED
v2, not run), which already runs the arms and has no timing metric - and
which would then also pick up the prediction that absorbing the smooth
share shrinks apportionment stickiness as well as section 3.1's mode.

**Section 7's recommendation stands.** The move census and the composition
probe still go first. This addendum adds one measurement (Stage R0) that
rides the first of them, and re-specifies a third.

### 12.7 What this addendum could not settle

- **End-to-end acceptance at ship scale.** Both enumerations are `m = 1`,
  one predictor, `n ~ 400`, at most 2950 trees. They establish the
  mechanism and refute the depth-1 identity; they predict nothing at 75
  trees and `p = 10`. Stage R0 is the right measurement.
- **Whether the `q`/`pi` weight spread degrades with `p`.** Both
  enumerations are single-predictor. The formula says `p` enters only
  through `B_w`'s dispersion. Worth one more R0 cell.
- **Whether a within-forest diversity prior is worth building.** B1 removes
  the argument that said no. It does not supply one that says yes, and the
  only published within-ensemble construction (DP-Forests) reports no
  mixing diagnostic at all.
- **Whether the missing-column double-halving should be repaired.** Owned
  by `grow-from-root-categorical-scan` S0 with an open VD fork. This
  addendum only records that a regrow raises the stakes from start-quality
  to importance-weight bias.
- **Everything in the memo's PASS-VERIFIED external tier.** Not re-fetched
  by this pass and therefore not carried into this record (12.4).

### 12.8 Provenance

```
repo          /Users/vdorie/Repositories/dbarts, branch bartcore
tip           ef7335d; working tree clean at the start of this pass
scope         research only - no source change, no commit, nothing scheduled
seeded by     TODO: tree-mixing-proposals, two further VD directions
              commissioned 2026-08-09
working papers untracked {memo,critique,
              orchestrator-refinements}.md (gitignored; the memo is the
              survey, the critique is the blind review, this section is the
              adjudication)
code re-read  src/bartcore/grow.hpp (growTreeFromRoot 63-167, candidate
                assembly 73-130, missing coin 158-159)
              src/bartcore/scan.hpp (scanOrdinalCuts 75-120, naCode skip 94,
                occupancy sentinel 105-110)
              src/bartcore/chain.hpp (leaf scale 548-549; prior defaults 51;
                blocks install 621 / installBlockMasks 3882-3905; sweep 985-
                1010; k hyperprior 1106-1109; DART 1110-1116; sigmaDf 247;
                growForestFromRoot 1332-1345; regrow loop 1365-1395;
                per-forest CGMTreePrior 330)
              src/bartcore/model.hpp (constant leaf 155-215; monotone leaf
                520-540; CGM prior 2050-2135)
              src/bartcore/moves.hpp (logLikelihoodForBranch 50-80)
              src/bartcore/tree.hpp (depthOf 346-353, availability 545-600)
              src/R_interface_bartcore.cpp (nodeScale 255)
              R/dbarts.R (getTrees 1329, setOffset 1004), R/model.R
                (blocks 1061), inst/include/dbarts/dbarts.h (setOffset 366)
in-repo docs  docs/plans/archive/grow-from-root-categorical-scan.md (S0),
              docs/plans/archive/composition-mixing-probe.md,
              docs/design/forest-ranef-interweaving.md, TODO entries
              tree-mixing-proposals and grow-from-root-categorical-scan
numerics      independent exact enumeration of the single-predictor tree
              space under the shipped CGM prior, the shipped constant-
              Gaussian leaf marginal and the builder's exact candidate
              weights; `q` asserted to normalize to 1 and the closed form
              `pi/q = Z_root prod_{w != root} [(1-g_w) + g_w B_w]` asserted
              against `log pi - log q` for every enumerated tree to 1e-9;
              plus the CGM balanced-tree log-prior table. Scripts were run
              out of repo and are not preserved; every input is named above
              and the enumeration is reproducible from them.
citations     every source in 12.4 fetched and read in this session; two
              carried as bibliographic record only, marked as such
```

---

## 13. Move-set A/B (2026-09-06)

### 13.1 The question

Should dbarts follow the BART package, which proposes birth and death only,
or bartMachine, which proposes grow, prune and change but no swap? Section 2
records that dbarts dispatched, when this was measured, to three kernels
under four labels at `birth_death 0.5, swap 0.1, change 0.4, birth 0.5`, and
section 5.7 requires
that every candidate be measured against what is already shipped. This is
the cheapest possible instance of that requirement: the two rival move sets
are already reachable at runtime through `proposal.probs`
([`dbarts`](../../R/dbarts.R), [`defaultProposalProbs`](../../R/model.R)) with no engine
change at all, so the contrast costs a grid of fits and nothing else.

Three cited sources disagree about what the answer should be. Pratola's
birth/death-only sampler collapses at low noise - acceptance ~4 percent,
90 percent coverage 53 percent - and mixes at sigma^2 = 1 [verified: arXiv
1312.1895 sec 2.2, as section 3.2 carries it]. Tan et al. found that
"restricting the move set [to grow and prune] does not substantially affect
R-hat, coverage, or RMSE" on their own Python implementation [verified:
arXiv 2406.19958 appendix L.6, as section 3.4 carries it]. bartMachine's
retention of change is attributed to variable-inclusion proportions rather
than to posterior means; **that attribution was not verified in this arc and
bartMachine is not in section 10's ledger**, so it is carried here as the
question's motivation and not as a citation. The three quantities they
disagree about - predictive accuracy, coverage, and variable inclusion - are
what this measurement takes.

### 13.2 Arms

Runtime only, everything else identical.

    A   birth_death 0.5, swap 0.1, change 0.4, birth 0.5   the then-shipped default
    B   birth_death 1.0, swap 0.0, change 0.0, birth 0.5   birth/death only
    C   birth_death 0.6, swap 0.0, change 0.4, birth 0.5   change, no swap

200 trees, 1000 burn-in, 2000 kept, one chain, one thread, `n.thin = 1`.
Five matched seeds per arm and design, the `grouped-mixing.R` idiom section
6.3 specifies: data seed indexed by design and replicate, sampler seed the
replicate, both shared across arms so every contrast is paired. Note that
200 trees and 4 chains are section 5.7's own baseline on the tree count but
not on the chain count; one chain is what makes a 75-fit grid affordable,
and it is why nothing between-chain appears below.

### 13.3 Designs

    (1) Friedman, n = 2000, p = 10 (5 noise), sigma = 1
    (2) as (1) at sigma = 0.25            the low-noise regime of section 3.2
    (3) 10 * 1[x1>.5] * 1[x2>.5] * 1[x3>.5] + 5 x4, n = 2000, p = 10,
        sigma = 1                         a depth-3 interaction, where a rule
                                          at a high node has to change
    (4) Friedman on the first 5 of p = 100, n = 1000, sigma = 1
                                          the variable-inclusion design
    (5) probit, n = 2000, p = 10, latent = standardized Friedman

The test set is 1000 rows, drawn once per design and held fixed across every
arm and seed; 25 evenly spaced rows of it carry the ESS. RMSE and 90 percent
pointwise coverage are of the true f on all 1000 rows - the latent index
under (5). ESS is `posterior::ess_basic`, median over the 25 points for f,
and of the sigma draws where sigma exists. Inclusion is `varcount` converted
to per-draw proportions and averaged over draws; "incl true" below is the
summed proportion on the signal columns.

### 13.4 Tables

Mean over the five seeds, min-max in parentheses.

    (1) Friedman, sigma = 1
    arm  RMSE                cover             ESS f         ESS sigma    incl true            wall s
    A    0.517(0.479-0.552)  0.904(0.88-0.93)   20( 11- 28)   75(  8-244)  0.693(0.671-0.707)  6.17
    B    0.515(0.493-0.534)  0.905(0.89-0.93)   21( 18- 23)   96(  6-246)  0.701(0.682-0.713)  5.16
    C    0.511(0.485-0.533)  0.915(0.90-0.94)   23( 20- 25)   16(  7- 24)  0.697(0.682-0.718)  6.25

    (2) Friedman, sigma = 0.25
    arm  RMSE                cover             ESS f         ESS sigma    incl true            wall s
    A    0.265(0.257-0.277)  0.714(0.69-0.75)    8(  3- 18)   19(  2- 84)  0.912(0.893-0.921)  5.95
    B    0.263(0.251-0.275)  0.725(0.70-0.75)    9(  6- 15)    5(  2- 13)  0.918(0.914-0.924)  5.31
    C    0.267(0.248-0.278)  0.714(0.68-0.74)   10(  5- 13)   25(  2- 81)  0.908(0.896-0.918)  6.04

    (3) deep interaction
    arm  RMSE                cover             ESS f         ESS sigma    incl true            wall s
    A    0.753(0.635-0.859)  0.949(0.94-0.97)   48( 30- 68)   36(  2-112)  0.482(0.458-0.512)  6.12
    B    0.762(0.660-0.825)  0.955(0.95-0.96)   49( 42- 57)  118(  2-333)  0.494(0.481-0.508)  5.34
    C    0.722(0.653-0.785)  0.959(0.95-0.97)   43( 35- 48)   15(  2- 32)  0.473(0.456-0.500)  6.13

    (4) sparse p = 100
    arm  RMSE                cover             ESS f         ESS sigma    incl true            wall s
    A    0.925(0.826-0.996)  0.900(0.88-0.92)   15(  9- 25)    2(  1-  3)  0.277(0.270-0.288)  4.60
    B    0.945(0.875-0.994)  0.885(0.85-0.92)   14( 10- 20)    6(  1- 18)  0.279(0.260-0.290)  4.24
    C    0.960(0.881-1.035)  0.881(0.85-0.92)   12(  9- 15)    9(  1- 30)  0.278(0.268-0.287)  4.63

    (5) probit
    arm  RMSE                cover             ESS f         ESS sigma    incl true            wall s
    A    0.316(0.288-0.353)  0.950(0.91-0.98)  234(117-294)      -         0.540(0.535-0.547)  6.63
    B    0.316(0.288-0.351)  0.948(0.90-0.98)  186(143-222)      -         0.541(0.533-0.548)  6.01
    C    0.318(0.289-0.352)  0.946(0.92-0.97)  241(206-265)      -         0.540(0.532-0.545)  6.65

Design 4 is the inclusion design, so it gets its own readout. Per-column
inclusion proportion, averaged over the five seeds:

    arm    x1      x2      x3      x4      x5     95 noise cols     noise cols
                                                   mean     max   above min true
    A    0.0726  0.0664  0.0497  0.0548  0.0336  0.0076  0.0112         0
    B    0.0684  0.0725  0.0503  0.0525  0.0350  0.0076  0.0117         0
    C    0.0739  0.0674  0.0493  0.0526  0.0349  0.0076  0.0106         0

Zero noise columns clear the weakest signal column's inclusion proportion,
in every arm and every one of the fifteen fits.

### 13.5 What separated and what did not

Sixty-eight paired contrasts - five designs, two contrasts against A, seven
metrics, less the two sigma cells probit does not have. The largest is
|t| = 2.53 on 4 degrees of freedom; the smallest Holm-adjusted p over the
whole family is 1.000. **No arm separates from any other on RMSE, on 90
percent coverage, on median f ESS, on sigma ESS, or on variable inclusion,
in any of the five designs.**

The resolution actually achieved, read off the paired standard errors: about
0.012 on RMSE, 0.02 on coverage, 6 on median f ESS, 0.01 on the summed
true-column inclusion share. Section 6.4's own pre-registered margin is 4x
the per-replicate standard error, which is 0.05 to 0.08 on coverage here. No
contrast reaches a third of it. Five pairs per cell rules out effects of the
size the three sources disagree about; it does not rule out small ones.

Two raw signals are worth recording even though neither survives
multiplicity, because they point the same way and both sit in the probit
design: B loses median f ESS against A (-48 +/- 22, four of five seeds
negative) and C raises the *worst*-point f ESS against A (+48 +/- 19, five
of five seeds positive). That is the change move helping, and the swap move
not, on the probit latent. It is a mixing signal, not the inclusion signal
bartMachine's retention of change is attributed to, and this grid cannot
establish it.

Sigma ESS does not discriminate and should not be carried forward at this
replicate count: it swings from 1 to 333 across seeds inside a single arm
and design. That is the slowness of the sigma chain at 200 trees, and it is
consistent with section 3.2's mechanism, not an arm effect.

Wall time is the one place the arms separate cleanly, consistently and in
the same direction in all five designs. Pooled, B costs 0.884x arm A and C
costs 1.008x; per design B/A is 0.836, 0.891, 0.873, 0.923, 0.907 and C/A is
1.013, 1.015, 1.002, 1.006, 1.003. Dropping change and swap buys 8 to 16
percent of the sweep. Dropping swap alone buys nothing measurable, which
section 6.1's addendum explains: 73.0 percent of swap proposals are no-ops,
so the move is already nearly free.

### 13.6 Reading against the three sources

- **Tan et al. is reproduced, on this package's own kernel.** Their
  Experiment 7 compared `{grow .5, prune .5}` against
  `{grow .25, prune .25, change .4, swap .1}` on a Python implementation and
  found no substantial effect on R-hat, coverage or RMSE. Arm B against arm
  A is that same contrast on dbarts, and it finds no effect on coverage or
  RMSE either, in five designs including two their battery did not have (a
  depth-3 interaction and a probit response). Section 8's standing
  conclusion - "the case for a *new* move cannot rest on 'dbarts has four
  moves'" - now rests on a dbarts measurement rather than on transfer from
  another implementation.
- **Pratola's low-noise collapse reproduces, and is not caused by the move
  set.** Design 2 drops 90 percent coverage from 0.90 to 0.71, which is
  section 3.2's established failure appearing on schedule. But it drops to
  0.714, 0.725 and 0.714 in arms A, B and C alike. Pratola's own sampler was
  verbatim "with birth/death proposals only", and section 3.2 already warns
  that the match to dbarts is on `(n, m, sigma)` and nothing else; this
  measurement adds that the move set is not the missing term. Restoring
  change and swap to a birth/death-only sampler does not repair low-noise
  coverage, so whatever repairs it is not in the shipped mixture.
- **The bartMachine premise is not reproduced on the design built to test
  it.** Design 4 puts Friedman signal on 5 of 100 columns, which is where
  inclusion proportions are supposed to degrade if change is dropped. The
  summed signal share is 0.277, 0.279 and 0.278 across arms; the 95 noise
  columns average 0.0076 in every arm; and no noise column in any arm or
  seed outranks the weakest signal column. If dropping change hurts
  variable-inclusion proportions, it does not do so here at n = 1000,
  p = 100, 200 trees and 2000 draws.

### 13.7 What this does not settle

- One chain per fit. Every between-chain statistic section 6.3 prefers - the
  pooled between-chain standard deviation of time-averaged inclusion, above
  all - is absent, and section 8 already records that at 75 trees no
  structural statistic detects mode collapse at feasible replicate counts.
  A null here is a null on within-chain readouts.
- Five pairs. The margins are stated above; small effects are not excluded.
- The three arms are mixtures over the *shipped* kernels. Nothing here
  speaks to a new kernel, and in particular nothing here weakens or
  strengthens the section 4.2 cut move, whose premise is that the shipped
  change move is badly aimed rather than that there are too few moves.
- Acceptance rates were not instrumented for this grid; section 6.1's
  addendum supplies them from a separate run at a separate configuration.

The practical consequence, stated plainly: there is no measured reason to
change the shipped default, and no measured reason to fear either
alternative. A user who sets `proposal.probs` to birth/death only gets the
same answers 8 to 16 percent faster on these five designs. That is a fact
worth documenting; it is not an argument for moving the default, which would
need the harm battery section 6.4 requires and does not have.

### 13.8 Provenance

```
repo          /Users/vdorie/Repositories/dbarts, branch bartcore
measured at   916271f3, re-checked at c585ba2c (docs only in between, no
              source, so the measurement is live at that tip)
build         private library installed from a clean `git archive HEAD`
              export; dbarts 1.0.0, R 4.6.1, posterior 1.7.0, arm64 macOS
grid          3 arms x 5 designs x 5 seeds = 75 fits, run in the foreground,
              single-threaded; 4.2 to 6.7 s per fit
scope         measurement only - no source change, no default change, nothing
              scheduled
caveat        the host carried a load average of 9 to 16 throughout. Absolute
              seconds are indicative; the B/A ratio holds in all five designs
              and all 25 pairs, which is what carries the wall-time claim.
scripts       run out of repo and not preserved; every input is named above
              and the grid is reproducible from `proposal.probs` alone
```

---

## 14. Recovery after a response swap (2026-09-06)

### 14.1 The question

Section 13 measured three proposal mixtures on five one-shot designs and
found no difference. Every one of those designs fits a fixed response once.
dbarts' distinguishing use is not that: it is a `dbartsSampler`
([`dbartsSampler`](../../R/dbarts.R)) inside a larger Gibbs loop whose response
moves between sweeps, which is what stan4bart and bartCause do and what
section 3.2's scope note already names as the regime where a stale build
scale bites ("`setResponse(updateScale=FALSE)` inside a larger Gibbs
sampler"). How fast the trees re-adapt after the response moves is measured
by none of the three cited sources and by none of section 13's designs. It
is the criterion on which change and swap should be kept or dropped for
THIS package, and it is reachable the same cheap way section 13 was:
through `proposal.probs` ([`dbarts`](../../R/dbarts.R),
[`defaultProposalProbs`](../../R/model.R)), with no engine change.

### 14.2 Design

Section 13's three arms, unchanged:

    A   birth_death 0.5, swap 0.1, change 0.4, birth 0.5   the then-shipped default
    B   birth_death 1.0, swap 0.0, change 0.0, birth 0.5   birth/death only
    C   birth_death 0.6, swap 0.0, change 0.4, birth 0.5   change, no swap

n = 2000, p = 10 uniform predictors, 200 trees, one chain, one thread,
`keepTrees` FALSE. Burn 1000 sweeps on y1 = Friedman(x) + N(0, sigma^2),
then `setResponse(y2, updateScale = FALSE)`
([`setResponse`](../../R/dbarts.R)) and 400 sweeps ONE AT A TIME, reading the
current fit after each. Four shifts, y2 carrying the SAME noise draw as y1
so the response moves because the outer block moved and not because the
data were redrawn:

    control  f2 = f1                                 no swap; the level to recover to
    small    f2 = f1 + 2 x6                          a new linear term on a noise variable
    large    f2 = Friedman with x6..x10 in x1..x5's roles
    offset   f2 = f1 + b_g, 20 groups, b_g ~ N(0, 1) a random-intercept block's shape

Two noise levels, sigma = 1 and sigma = 0.25, and section 6.3's matched-seed
idiom: data seed indexed by replicate, sampler seed the replicate, both
shared across arms, so every contrast is paired. The test set is 1000 rows
drawn once and held fixed across every arm, seed and shift; 25 evenly spaced
rows of it carry the ESS.

Recovery sweeps is the first sweep t >= 10 at which the 10-sweep median
window of the RMSE trajectory comes within 10 percent of a target. Two
targets are reported: the arm's OWN mean over sweeps 300-400, and a COMMON
one pooled over the three arms within (shift, sigma, seed). The common
target is the honest one, since an arm that settles at a worse level would
otherwise look fast for recovering to its own worse level; the two agree on
every ranking below. ESS is `posterior::ess_basic` over post-swap sweeps
200-400, median over the 25 points. Wall time is 400 sweeps issued in one
`run` call after the measured loop, so the per-call R overhead is out of the
ratio.

### 14.3 The confirmation set

The first grid - four shifts, two sigmas, three arms, five matched seeds,
120 fits - found exactly one cell where the arms can differ: the large
shift, the only shift whose post-swap trajectory has a transient at all.
There arm B was slower and worse than A in 5 of 5 seeds on every level
metric at both noise levels, largest |t| 6.37 on 4 df, and nothing survived
Holm over the 144-contrast family (smallest adjusted p 0.448) or over the
36-contrast large-shift family (0.112).

That cell was then re-run on 15 FRESH seeds, independent of the five that
selected it, with the primary contrasts fixed in writing before the run:
B - A on recovery-to-common-target and on RMSE at post-swap sweep 50, at
each sigma, four tests, Holm over the four, predicted direction positive.
C - A on the same four, the ESS and the remaining level metrics were
secondary and not gating.

### 14.4 Tables

Mean over seeds, min-max in parentheses. Control, small and offset are the
five seeds; large is all twenty.

    recovery sweeps to a common target, train RMSE
    sigma shift     n     A                B                C
    1     control   5   10.0(10-10)      10.0(10-10)      10.0(10-10)
    1     small     5   10.0(10-10)      10.0(10-10)      10.0(10-10)
    1     offset    5   10.0(10-10)      10.0(10-10)      10.0(10-10)
    1     large    20   72.0(52-89)      89.1(58-124)     80.8(55-121)
    0.25  control   5   10.0(10-10)      10.0(10-10)      10.0(10-10)
    0.25  small     5   10.4(10-12)      10.6(10-13)      10.4(10-11)
    0.25  offset    5   10.0(10-10)      10.0(10-10)      10.0(10-10)
    0.25  large    20  213.5(163-264)   238.7(177-299)   214.4(166-289)

    post-swap ESS, sweeps 200-400, median over the 25 fixed test points
    sigma shift     n     A                B                C
    1     control   5   11.0(7.9-12.6)   11.1(8.3-13.4)   11.5(9.3-16.0)
    1     small     5    9.2(5.0-11.7)    8.7(4.5-12.6)   10.9(7.8-15.1)
    1     offset    5   10.0(5.9-14.4)   13.0(6.8-17.0)    9.8(5.2-14.8)
    1     large    20   10.3(5.4-16.6)   11.3(6.9-16.8)    8.8(3.8-11.5)
    0.25  control   5   12.4(9.3-14.7)   13.1(8.2-21.9)    9.9(7.8-11.9)
    0.25  small     5    9.2(6.1-12.5)    9.9(5.5-13.1)   11.2(4.6-14.4)
    0.25  offset    5    9.1(4.6-13.2)    9.2(4.7-12.0)   13.2(8.1-17.8)
    0.25  large    20    6.8(2.9-14.3)    8.3(4.1-15.0)    7.6(4.2-12.3)

Test RMSE against the true f2 along the transient, and the steady state the
matching no-swap control holds:

    sigma shift  arm    @10    @50   @100   @400   control steady state
    1     large  A     1.479  0.897  0.790  0.747  0.734
    1     large  B     1.539  0.925  0.800  0.745  0.728
    1     large  C     1.493  0.896  0.806  0.744  0.747
    0.25  large  A     1.447  0.686  0.505  0.354  0.321
    0.25  large  B     1.605  0.802  0.558  0.378  0.322
    0.25  large  C     1.452  0.706  0.521  0.365  0.315

The four primary contrasts on the confirmation seeds, 14 df, Holm over the
four:

    contrast                       sigma    diff       se     t       Holm p
    B - A, recovery(common target)  1      +16.07    4.10   3.92     0.005  *
    B - A, RMSE at sweep 50         0.25    +0.108   0.025  4.28     0.003  *
    B - A, recovery(common target)  0.25   +24.67   10.53   2.34     0.069
    B - A, RMSE at sweep 50         1       +0.017   0.021  0.79     0.443

C - A on the same four: +7.67 (p 0.17), -3.53 (0.81), +0.0018 (0.93),
+0.019 (0.32).

### 14.5 What separated and what did not

**Three of the four shifts have no transient at all.** The small shift, the
offset shift and the no-swap control fire the recovery criterion in the
first available 10-sweep window, in every arm, at both noise levels, in
every seed. A new linear term on a variable the forest was not splitting on
is absorbed by leaf redraws; the group offsets are not a function of any
predictor, so the forest cannot represent them and absorbs them into sigma
(which rises to 1.25 at sigma = 1 and 0.89 at sigma = 0.25) rather than into
the trees. There is no move-set question in these cells, because there is no
structural work for a move to do.

**The large shift separates, and B is the arm that loses.** Retargeting the
signal onto five variables the forest was not splitting on is the only shift
that makes structural work, and there birth/death only takes 16 more sweeps
to recover at sigma = 1 (Holm p 0.005) and sits 0.108 higher in RMSE at
post-swap sweep 50 at sigma = 0.25 (Holm p 0.003). The direction is the same
in every level metric at both noise levels and in the recovery metric under
either target.

**The gain belongs to change, not to swap.** Arm C carries change at the
same 0.4 weight with swap set to zero and is indistinguishable from A on all
four primary contrasts and on every secondary one, in either direction; B,
which drops both moves, is the only arm that separates. Swap contributes
nothing measurable to recovery, which is what section 6.1's addendum already
predicts from the other side: 73.0 percent of swap proposals are no-ops
([`swapMove`](../../src/bartcore/moves.hpp)).

**Recovered mixing does not separate.** Post-swap ESS is flat across arms in
all eight shift x sigma cells; the largest |t| on ESS anywhere in the
confirmation set is 1.38 on 14 df. The move set changes how fast the trees
get back, not how well they mix once back.

**Wall time reproduces section 13 exactly.** Over the 120-fit grid B costs
0.887x arm A and C costs 0.993x, and the per-pair B/A ratio is below 1 in
all 40 pairs. So B trades an 11 percent sweep saving for a 20 to 25 percent
longer recovery on the one shift that needs one.

### 14.6 Reading

Section 13's null holds for a fixed response and does not extend to a moving
one. On the criterion that matters for this package's distinguishing use,
the change move earns its 0.4 and the swap move does not earn its 0.1: the
one place the arms separate is a large response swap, and the arm that
carries change without swap matches the shipped default there. That is
consistent with section 2's mechanics - change is the only shipped kernel
that can install a new split VARIABLE at an existing interior node
([`changeMove`](../../src/bartcore/moves.hpp)), which is precisely the edit a
retargeted signal demands, while swap only exchanges a parent's rule with a
child's and cannot introduce a variable the tree does not already carry
([`metropolisJumpForTree`](../../src/bartcore/moves.hpp)).

**This is what dropped swap out of the default.** Arm C is the mixture that
ships: change keeps its 0.4 and swap's 0.1 goes to birth/death. The move
itself stays in the kernel - at one tree it is the only proposal that rotates
a rule up, which no default at production forest sizes can speak to
([Removing the swap tree-proposal](swap-removal.md#removing-the-swap-tree-proposal)).
Change stays too, which this section is the evidence for. What the ledger now
says about dropping change: a user who sets `proposal.probs` to birth/death
only gets section 13's same answers 8 to 16 percent faster on a fixed
response, and pays 20 to 25 percent more sweeps to re-adapt when the response
moves under a sampler embedded in an outer loop.

### 14.7 What this does not settle

- One chain per fit, so nothing between-chain, exactly as in section 13.
- 400 post-swap sweeps is not enough at sigma = 0.25: no arm returns to its
  control's steady-state RMSE by sweep 400 (0.354, 0.378, 0.365 against
  0.321, 0.322, 0.315), so that cell ranks an unfinished transient. At
  sigma = 1 every arm has returned.
- The large shift is a TOTAL retarget of the signal. A real Gibbs loop's
  sweep-to-sweep response move is small, and the two shifts here that are
  closer to that size show no arm effect because they show no transient.
  What is established is that the move set matters when the response moves
  far enough to require new split variables, not that stan4bart's or
  bartCause's own sweep-to-sweep moves are that large. Measuring the size of
  move those samplers actually produce is the obvious next question and is
  not asked here.
- Five seeds outside the large shift. Those cells rule out effects the size
  of the large-shift effect, not small ones.
- The arms are mixtures over the SHIPPED kernels. Nothing here speaks to a
  new kernel, and in particular nothing here bears on the section 4.2 cut
  move, whose premise is that change is badly aimed rather than absent.
- Recovery is measured on RMSE against the true f2. No coverage, no
  variable inclusion, no sigma-chain readout is taken along the transient.

### 14.8 Provenance

```
branch        bartcore
measured at   edf5a85b (section 13's own landing tip; no source touched
              since, so the measurement is live at that tip)
build         private library installed from a clean `git archive HEAD`
              export; dbarts 1.0.0, R 4.6.1, posterior 1.7.0, arm64 macOS
grid          120 fits (3 arms x 4 shifts x 2 sigmas x 5 seeds) plus a
              90-fit confirmation set (3 arms x 1 shift x 2 sigmas x
              15 fresh seeds), run in the foreground, single-threaded;
              1801 sweeps and about 3.5 s per fit
scope         measurement only - no source change, no default change,
              nothing scheduled
caveat        the host carried a load average of 8 to 28 throughout.
              Absolute seconds are indicative; the wall-time claim rests
              on the per-pair B/A ratio, below 1 in all 40 pairs
scripts       run out of repo and not preserved; every input is named
              above and the grid is reproducible from `proposal.probs`,
              `setResponse` and the four shift definitions alone
```

---

## 15. Proposal brainstorm: geometry, latents, informed moves, rotated inputs (2026-09-07)

### 15.1 The question and the four lenses

VD, 2026-09-07: "I still wonder if there aren't proposals that use
geometry, latents, and least-favorable-directions or something like that."
Four independent lens surveys were run against sections 2 to 7 and 12, the
move census in sec 6.1, and `benchmark-surfaces.md` sec 10: tree and
partition geometry; latent and auxiliary constructions; informed and
gradient-like proposals; rotations of the input space. Each had to state a
proposal precisely enough to write its Metropolis-Hastings correction, say
why it is valid, price it against one cut scan, name a measured deficit,
and supply a falsifier runnable in a day.

**This was a literature-anchored pass**, and the novelty column of sec 15.3
reports what that produced honestly: most rows are known elsewhere, several
are standard MCMC constructions with no tree-ensemble instance the four
searches could find, and none is new outright. A first-principles round -
mechanisms derived from dbarts' own structure with the literature held out -
is the natural sequel, and it is
[16. Proposal brainstorm, second round: first principles under a novelty gate (2026-09-07)](#16-proposal-brainstorm-second-round-first-principles-under-a-novelty-gate-2026-09-07). A refutation pass then checked
each lens's load-bearing claim against source or arithmetic, and its
verdicts ride the last column. **No recommendation to build follows**; the
ranking is by (evidence + mechanism) / cost and nothing here is scheduled.

### 15.2 The two facts the lenses agreed on

**(A) The posterior lives on the partition of the rows, so a
partition-preserving move has likelihood ratio exactly 1.**
[`ConstantGaussianLeaf`](../../src/bartcore/model.hpp) scores a leaf through
`(sumWeights, sumWeightedResponse)` and nothing else, and the raw sum of
squares it drops is additive over the members, so it cancels under any
repartition. The refutation pass extended this rather than assuming it, and
it holds for every other integrable leaf on its own statistic:
[`ConstantVarianceLeaf`](../../src/bartcore/model.hpp) on `(n, ssr)`,
[`LinearGaussianLeaf`](../../src/bartcore/model.hpp) on `(U'WU, U'Wz, z'Wz)`,
[`GPGaussianLeaf`](../../src/bartcore/model.hpp) on the leaf's own rows - all
functions of the member set alone, the response families changing only the
working response and weights those statistics read. **The exception is
the monotone leaf**: [`MonotoneConstantGaussianLeaf`](../../src/bartcore/model.hpp)
adds a truncation term whose bounds come from
[`monotoneNeighborBounds`](../../src/bartcore/model.hpp), which reads leaf BOXES
through [`monotoneLeafBox`](../../src/bartcore/model.hpp) and the frozen values of
the other leaves - both of which a rule set can move while holding the row
partition fixed. A fibre move must scope monotone forests out, the exclusion
[12.6 Ranked disposition](#126-ranked-disposition) already writes for the
regrow census.

So representation multimodality
([3.1 Many tree arrangements, one fitted function (ESTABLISHED)](#31-many-tree-arrangements-one-fitted-function-established))
and the rooting lock of
[10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)
are connectivity failures on the fibre, not likelihood barriers, and a fibre
move's acceptance is prior ratio times transition ratio, both already in
[`CGMTreePrior::treeLogProbability`](../../src/bartcore/model.hpp). dbarts owns one
such move and ships it at zero: in [`swapMove`](../../src/bartcore/moves.hpp)'s
both-children-share-a-rule branch the four grandchild cells are preserved as
a set with the two cross cells exchanging slots, so where those slots are
leaves the partition and the likelihood are untouched. It fires only when
[`Tree::rulesAreEqual`](../../src/bartcore/tree.hpp) holds for two siblings, which
is rare by construction and counted nowhere.

**(B) dbarts gets exact neighbourhood scores where the discrete-MCMC
literature pays a Taylor surrogate for them.** Gibbs-with-gradients and the
discrete Langevin family spend a first-order expansion to avoid `O(d)` exact
evaluations; [`scanOrdinalCuts`](../../src/bartcore/scan.hpp) returns the collapsed
marginal for EVERY cut of one variable in one pass over a node's members,
scoring the missing-direction bit rather than drawing it, so importing the
approximation would be strictly worse. Three boundaries the refutation pass
fixed. It is ordinal-only:
[`scanCategoricalPartitions`](../../src/bartcore/scan.hpp) is marked INIT-ONLY and
"not a valid Metropolis-Hastings neighborhood and must not be reused as
one", so on a mixed design the conditional over rules is not enumerable in
`p` scans at all. It is scalar-leaf-only, being templated on
`ScalarLeafModel`. And it is exact only at a node whose two children are
leaves (a "nog" node), where the score is a two-way partition of the node's
members; below that frontier a moved cut reroutes members through a fixed
skeleton and no prefix scan exists. Neither this document nor
`perturb-move.md` had drawn that line.

### 15.3 Cross-lens ranking

Cost is per proposal in cut-scan units (one
[`scanOrdinalCuts`](../../src/bartcore/scan.hpp) pass over a node's members for one
variable); a sweep makes `m` proposals, 0.4 of them changes. Novelty is
"known elsewhere" (named), "new to BART" (a standard construction with no
tree-ensemble instance found) or "new outright"; nothing here is new
outright.

| # | mechanism | what it is | validity | cost | deficit | novelty | one-day falsifier | refutation |
|---|---|---|---|---|---|---|---|---|
| 1 | Rule Gibbs at a nog node | scan all `p` variables over a nog node's members, draw the rule from prior x marginal | Gibbs, acceptance 1 - the neighbourhood is closed | `p` scans; sec 4.5's measured 10.4x at `p` = 10, 53-56x at `p` = 50 | change's aim (3.77 percent accepted, median rejected -62.34 default, -143.45 lownoise) and 3.2's low-noise freeze | new to BART: collapsed Gibbs on a closed discrete neighbourhood | extend [`cutProbe`](../../src/bartcore/moves.hpp) to log the weight entropy, P(incumbent) and the nog share; P(incumbent) near 1 at `lownoise` kills it | CONFIRMED for ordinal availability sets; REFUTED where any available variable is categorical, and for every leaf model the scan cannot serve. MEASURED (sec 6.1's third 2026-09-07 addendum): joint weight entropy and P(incumbent) at a nog node range from 0.334 nats / 0.737 at `lownoise`, the least there is to buy, to 6.425 nats / 0.0017 (median rank 26.5 of up to 3000) at `c1`, the most; the falsifier's near-1 kill condition does not occur in any cell measured. BUILT and BENEFIT-TESTED as `rule_gibbs` (nog-gibbs.md): on the C1 four-chain cell it raises the summed minimum ESS by +21.5 and +22.1 over twenty matched pairs at two seed blocks against a +8 bar, at 2.21 cut-scan sweep-equivalents, so its pre-registered kill does not fire; it regresses that cell's 95 percent coverage by -0.022 and -0.026 against a -0.010 margin, which blocks a nonzero default share ([6. Benefit, pre-registered](nog-gibbs.md#6-benefit-pre-registered)) |
| 2 | Cut Gibbs at a nog node | item 1 restricted to the incumbent variable | Gibbs, acceptance 1 | 1 scan - the pass change already makes | same | new to BART, as item 1 | rides item 1's probe | as item 1; it dominates perturb at `w` = 1 only if the conditional is not a point mass. MEASURED: the cut-restricted entropy is never near zero (0.328 to 3.426 nats across cells) and P(incumbent) never near 1 (0.034 to 0.746), so the conditional is not a point mass anywhere in this grid, and item 2 dominates a `w` = 1 perturb everywhere measured, most narrowly at `lownoise` |
| 3 | Informed death, uniform birth | weight each nog node by `sqrt` of its pruned posterior ratio; leave birth alone | ordinary MH: the forward gains `w(v)/W(T)`, the reverse is unchanged | 0 scans; `O(#nog)` arithmetic plus three availability walks per candidate | death's median rejected -48.93 default, -77.13 lownoise, at 10.88 percent; and 3.5's size walk | new to BART: one-sided locally-balanced weighting | log the full nog weight vector and the realized uniform pick at every death; dead if the uniform pick already sits near the weighted mode | CONFIRMED for the constant leaf; the cited "[`Tree::computeLeafStats`](../../src/bartcore/tree.hpp) already forms a parent as left + right" is REFUTED - it accumulates over the index span. MEASURED, and DEAD by its own kill criterion: the normalized weight vector over nog nodes is a near point mass everywhere - median entropy 3e-5 nats at `c1` to 5e-20 at `lownoise` - and the uniform pick already sits at that mode 83.1 to 97.1 percent of proposals overall, falling to 47.7 to 49.9 percent (indistinguishable from chance) once there are two or more candidates; the weighting has nothing left to inform |
| 4 | DART-on census re-run | re-run the census with `dart = TRUE`, which already routes `splitProbabilities` into change's variable redraw | shipped and valid | 0 | the variable half of change's aim | known elsewhere: Linero 2018, shipped here | one census re-run, no engine work | code claim CONFIRMED; the inference is REFUTED - DART changes the target, so this is not a decomposition |
| 5 | Same-variable rotation | rotate a parent-child interior pair splitting on one ordinal variable | fibre move: likelihood exactly 1, the correction is the rotatable-node count | 0 scans; one index permutation plus `O(#subtree)` prior factors | the rooting lock | known elsewhere: tgp's `rotate`, which Pratola dismisses without a number | count parent-child interior pairs sharing a split variable per sweep in the four cells; kill under ~1 percent of proposals | CONFIRMED as arithmetic - three cuts on one axis re-bracket the same intervals, so no merge enumeration exists; the rate is unmeasured |
| 6 | Partition-preserving restructure | hold the leaf blocks, re-derive the whole rule set top-down | fibre move; `q` factorizes over nodes and is re-run on `T'` | `O(p)` interval scans per interior node, no likelihood evaluation | the rooting lock at `m` = 1, and 3.1 | known elsewhere: Wu, Tjelmeland and West, already set aside in sec 5.2 | on P2's five seeds compute both rootings' leaf partitions and their variation of information; VI > 0 kills the fibre program for that cell | CONFIRMED as a construction; whether the fibre is larger than a point at production scale is unmeasured |
| 7 | Rows-crossed step dial | size a displacement window so a target number of ROWS crosses, not a target number of grid positions | valid if the width is a deterministic function of the current state, recomputed at `T'`; a chain-history average is adaptation | free rider on the scan | change's aim, and sec 6.1's window grid, which is in grid positions | new to BART: state-dependent step size | add a rows-reassigned column to [`cutProbe`](../../src/bartcore/moves.hpp) and re-run the four cells | the metric is CONFIRMED to leading order in `1/sigma^2`; its pooling prediction across birth, death and change is REFUTED - a leaf-count change leaves an Occam term that does not cancel. The leaf-count spread behind that term is now MEASURED directly (sec 6.1's third 2026-09-07 addendum) rather than asserted - mean leaves per tree runs 2.5 to 3.8 across the four cells, 1.1 on the BCF treatment forest - which is the size the pooling prediction would have to absorb |
| 8 | Conditional SMC as a mixture component | with probability `p_pg` one tree's step is a conditional-SMC sweep with the incumbent clamped | a mixture of pi-invariant kernels; no correction at the mixture level | `C` x (#interior) scans on the sweeps it fires, about 40 at `C` = 10 | change's aim, and the He-Hahn ESS | known elsewhere: PG-BART, PyMC-BART; the mixture framing is new to BART | generator-only: run the build, log the returned particle against the clamped one, discard; above ~95 percent retention it is overhead | CONFIRMED, with two conditions the lens omits: `p_pg` must not depend on the state, and the build must respect the veto's support |
| 9 | Fixed oblique augmentation | fit on `[X, XW]` with `W` fixed before the run | exact - the kernel is untouched; the MODEL changes | 0 in the proposal; per-node availability goes `O(p)` to `O(p+K)` | P6's outer failure, measured here at 0.314 bias and 0.590 coverage | known elsewhere: rotation forests, feature augmentation | add oracle and random-20 arms to the P6 script on matched seeds, with a mandatory harm clause on P2's duplicate-column null | premise CONFIRMED: P6's mu has a shelf at the line `x1 = x2`; but P6 is not the only oblique cell - P7's setting A is a function of a linear index |
| 10 | Outer projection step | propose `W'` under a Stiefel prior, install `XW'` through `$setPredictor`, accept on the fit | valid: `W` to partition is deterministic, `q` symmetric, rollback exact | `O(m n)` per outer step; budget one extra sweep, not a proposal | none measured; a model extension at the package's design centre | known elsewhere as Bayesian single-index / projection pursuit; new to BART as a forest-conditional | grid small `W` perturbations on an oblique surface, recording what fraction `$setPredictor(forceUpdate = FALSE)` accepts | "prior ratio 1 iff the grid is pinned" is REFINED - the CGM prior reads no predictor value at all, so pinning gives 1 for every column - and the "only if" is REFUTED: what breaks under a re-derived grid is a non-invertible remap |

Runners-up, below the ten on cost or on target and not re-argued here:
rotated-PAIR split rules under an angle prior (GP-BART's construction, the
only rotation item with a mixing story); a Jain-Neal restricted-Gibbs launch
for categorical rules, where
[`drawCategoricalRuleFromPrior`](../../src/bartcore/moves.hpp) is the worst-aimed
draw in the kernel; a nog-shaped surrogate scan with exact MH below the nog
frontier; multiple-try Metropolis where no scan exists; locally-balanced
birth over the whole tree, which is
[4.5 Informed birth/death over the shared cut scan](#45-informed-birthdeath-over-the-shared-cut-scan)
unchanged and which the census de-prioritizes because birth is already a
near miss; one private rotation per tree; full oblique rules in the engine;
a soft-gate bridge; and Linero's retrospective auxiliary leaf values.

### 15.4 The readings taken of "least-favorable directions"

- **Steepest-change / locally-informed.** Three lenses took it; it yields
  items 1, 2, 3 and 7. Its exact content is that a deterministic argmax edit
  has no reverse density, so the reversible softening is a balancing
  function - and on a CLOSED neighbourhood that function is the identity,
  which makes the move a Gibbs step rather than an informed proposal at all.
- **Flat directions, the fibre.** The geometric lens inverted the phrase and
  asked which directions the posterior does not see. This produced the most:
  items 5 and 6, fact (A), and the observation that a shipped move already
  sits on the fibre. It is the opposite of a least-favorable direction in
  the information-geometry sense.
- **The semiparametric least-favorable submodel.** No proposal: it describes
  an estimator's fluctuation, not a transition kernel, and MH restores the
  target whatever direction one aims at. What it supplies is a READOUT, and
  [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s
  minimum ESS over 25 fixed points is that idea already, by brute force.
- **Adversarial mode targeting.** Requires knowing the modes, i.e. sec
  12.5's A5 mode atlas, adaptive and stale by construction.
- **Directions in the INPUT space**, the rotation lens's reading: a model
  question rather than a kernel one, producing items 9 and 10.

### 15.5 Refuted along the way

- **BHV tree space, geodesics and Riemann-manifold methods.** BHV
  coordinates are interior edge lengths over a fixed leaf-label set, which
  regression trees do not have, and the collapsed marginal is piecewise
  constant between adjacent order statistics, so the Fisher information is
  degenerate. Tree kernels and leaf-assignment embeddings fail next door: a
  point in embedding space cannot be decoded to an admissible tree without
  an enumerator, and once that exists it is item 6.
- **Mask complement as a move.** Complementing a categorical direction mask
  and swapping the two child subtrees preserves partition, depth and prior
  exactly - a free `Z2` relabelling the canonical gauge does not quotient,
  changing no fit. A caution that an informed categorical scan must not
  double-count, not a move.
- **C1's Single index as the rotation falsifier**, and **a learned rotation
  as a cure for the He-Hahn ESS.** The first is radial and exactly
  rotation-invariant on the independent design, which is what the two
  sentence corrections in `benchmark-surfaces.md` sec 2.2 and 6.3 record;
  the second misreads a chain-exploration deficit, since a
  reparameterization leaving an equally multimodal posterior does not move a
  chain that sits in one place.
- **Per-observation latent leaf membership, and a shadow forest sampled from
  a flattened target.** Not every membership vector is realizable by an
  axis-aligned tree and the realization map's density is intractable; the
  exchange algorithm needs an exact draw from the auxiliary model, and no
  exact tree-posterior sampler exists.
- **A latent per-tree temperature, and pseudo-marginal for the Gaussian
  leaf.** Trees share sigma and are exchangeable, so a per-tree exponent has
  no Gibbs conditional and breaks the target; the Gaussian marginal is exact
  and closed form, so an unbiased estimator only adds variance.
- **Gibbs-with-gradients' central trick** (unnecessary here, for fact (B));
  **delayed rejection for change** (no cost control where the first stage
  rejects 96 percent of what it scores - delayed ACCEPTANCE is the live
  direction, for the expensive leaf models); **a latent continuous cut
  location** (the one-position probe already accepts 40.65 / 27.15 / 51.38 /
  33.34 percent, so the grid is not the binding constraint); **treed models
  with linear splits as prior art** (there are none); **informed hyperplane
  proposals**, whose one published attempt reports the opposite of the
  intuition; and **tempering**, again.

### 15.6 What this could not settle

- **The nog share of interior nodes at production settings**, on which items
  1, 2 and 3 all live: for a tree with `L` leaves it lies between 1 and
  `L/2`, and nobody has counted. **Whether the nog conditional is a point
  mass at the incumbent** decides items 1 and 2, and the perturb probe's
  median log ratio is a displacement statistic, not the conditional.
- **Whether P2's two rootings induce the SAME leaf partition**, on which
  item 6 rests and where the record says only "exactly equiprobable"; and
  the rate of same-variable parent-child interior pairs, which prices item 5
  and which Pratola calls "usually not satisfied" with no number and no
  stated `p`.
- **Whether a valid categorical neighbourhood can be enumerated cheaply at
  all.** Fact (B) says the shipped scan cannot supply one and the exact
  space is `2^R - 2`.
- **Whether the leading-order partition metric survives the `lownoise`
  cell**, where the exponent is largest and the sampler worst; until item
  7's falsifier runs, items 5 to 7 share an unverified first-order argument.
- **Whether PG's ensemble advantage survives at 75 trees.** Its evidence is
  `n` = 2000 at 200 trees; sec 5.3's pre-registered read, effective samples
  per second, is still right and still unmeasured.
- **Whether the empty-leaf veto locks an outer `W` step**, the veto being a
  hard support constraint on `W` given the forest; whether GP-BART's rotated
  splits were ever ablated from its GP leaves; and whether any published
  work puts a prior on a full rotation MATRIX sampled jointly with a forest.
  Both rotation search legs found none: apparently nonexistent, not verified
  absent.

### 15.7 Provenance

```
inputs        four independent lens surveys, 2026-09-07 (tree and partition
              geometry; latents and auxiliary constructions; informed and
              gradient-like proposals; rotations of the input space), then
              a refutation and synthesis pass over the same records
base          cd32e6d0
scope         ranking and evidence only. No code, no default, no schedule.
              Two sentence corrections land with this section, in
              benchmark-surfaces.md sec 2.2 and 6.3, both consequences of
              the Single index finding
```

Literature verified during this arc; a source that could not be fetched is
named at the end and nothing above rests on it alone.

| # | Source | Verified | Where |
|---|---|---|---|
| 1 | He and Hahn, XBART factorial | Full text, sec 4.1 Table 1: "10 sqrt(a) + sin (5a); a = sum (x_j - gamma_j)^2", and "each element of X is drawn independently from a standard Gaussian distribution". The basis of the radial finding | arXiv 2002.03375v4 |
| 2 | Wu, Tjelmeland and West, Bayesian CART | "the likelihood function depends on T only through the induced partition of observations to leaves"; the move leaves "unchanged the partition of observations into terminal nodes"; "The computational cost of this is proportional to the number of candidate predictor variables p" | www2.stat.duke.edu/~scs/Projects/Trees/BayesianCART/WuWestPaper.pdf |
| 3 | Gramacy, Bayesian treed Gaussian process models | "Since the partitions at the leaves remain unchanged, the likelihood ratio of a proposed rotate is always 1. The only 'active' part of the MH acceptance ratio is the prior on T" | bobby.gramacy.com/prepo/gra2005-02.pdf |
| 4 | Pratola 2016 | rotate "efficiently traverses disparate regions of the model space along contours of equal probability"; on tgp's version, "it requires that all 3 internal nodes involved in a rotation split on the same variable. This constraint in general will usually not be satisfied" | arXiv 1312.1895 |
| 5 | Zanella, informed proposals | "Qg,s(x,dy) = g(pi(y)/pi(x))Ks(x,dy) / Zg(x)"; the condition "g(t) = t g(1/t) for all t > 0"; the "'naively informed' choice ... when g(t) = t"; "also the computational cost of sampling from the pointwise informed proposals increases" | arXiv 1711.07424 |
| 6 | Jain and Neal, split-merge | "Split-merge moves are produced by exploiting properties of a restricted Gibbs sampling scan"; the launch rule, quoted from a citing source, "It is only the last restricted Gibbs sweep which is used to compute the transition probability" | projecteuclid.org 10.1214/07-BA219 ; ar5iv 1406.0071 |
| 7 | Kim and Rockova, mixing rates for Bayesian CART | "locally informed proposal schemes leverage posterior information in the vicinity of the current state to propose the next state"; the twig move "attaches/detaches entire twigs (not just single nodes)" | ar5iv 2306.00126 |
| 8 | The discrete-MCMC surrogate family | Grathwohl et al., "We use a Taylor series computed on the underlying continuous function to estimate likelihood ratios of making discrete moves"; Rhodes and Gutmann, "the cost is O(d) evaluations of f for a d-dimensional problem"; Zhang, Liu and Liu, "DLP is able to update all coordinates in parallel in a single step"; Sun et al., an informed proposal "requires evaluating all energy changes in the neighborhood" | arXiv 2102.04509 ; 2208.00040 ; 2206.09914 ; iclr.cc/virtual/2022/spotlight/7061 |
| 9 | Lakshminarayanan, Roy and Teh, PG-BART | "it can mix faster since it can propose a completely different tree that explains the data"; "a change in an internal node that leaves any of the nodes in the subtree below empty will be rejected"; "We set the number of particles C = 10" | arXiv 1502.04622 |
| 10 | Linero, Generalized BART | "these auxiliaries end up canceling in the Metropolis-Hastings acceptance ratio"; "data augmentation can slow down mixing substantially, especially in cases where the outcome distribution is highly imbalanced"; and against itself, the RJMCMC algorithm "should be inferior in terms of mixing to the algorithm of Chipman et al. (2010)" | arXiv 2202.09924 |
| 11 | Linero and Yang, SoftBART | "SBART requires computing a likelihood contribution for each leaf-observation pair, whereas BART only requires a single likelihood contribution for each tree"; the reversible-jump latent over the number of trees "resulted in poor mixing" | arXiv 1707.09461 |
| 12 | Nguyen, Yee and Deshpande, oblique BART | "we choose to propose decision rules in grow moves from the prior"; "drawing overly-informed proposals can result in even slower MCMC exploration than drawing proposals from the prior"; "obliqueBART was about twice as slow as BART"; "the smallest overall average SMSE (0.296)"; "the smallest classification accuracy (0.846)". No ESS or mixing diagnostic anywhere | arXiv 2411.08849 |
| 13 | Maia, Murphy and Parnell, GP-BART | "an angle theta is sampled with equal probability from a predefined grid of 20 equally spaced values within the interval [0,pi]", through grow-rotate and change-rotate. The only published prior on a rotation sampled jointly with a forest | arXiv 2204.02112 |
| 14 | Blaser and Fryzlewicz, random rotation ensembles | "The best overall average rank of 6.10 (of 15)"; "in 67.8% of cases, random forests without rotation in 64.3% of cases"; "For categorical variables, rotation is unnecessary and ill-defined" | jmlr.org/papers/volume17/blaser16a/blaser16a.pdf |
| 15 | Jauch, Hoff and Dunson | "we parametrize the Stiefel and Grassmann manifolds ... using the Cayley transform. We derive the necessary Jacobian terms for change of variables formulas" | arXiv 1810.02881 |

NOT VERIFIED, and load-bearing for nothing above: Liu, Liang and Wong -
both reachable scans lack a text layer, so the multiple-try weight and
acceptance form come from a secondary source; Christen and Fox (2005)
delayed acceptance - no primary source fetched; Hohna and Drummond, Syst
Biol 61 - abstract only, the guiding weight paywalled; Rodriguez, Kuncheva
and Alonso, Rotation Forest - abstract only.

## 16. Proposal brainstorm, second round: first principles under a novelty gate (2026-09-07)

### 16.1 The rule of the round

VD's verdict on section 15: "These don't really feel novel - they all feel
like retreads of tree proposal mechanisms used elsewhere." Section 15's own
novelty column agrees: nothing there was new outright. A second round
therefore ran under the inverted rule - **derive from the posterior's
structure and from what this engine alone can do, with the literature held
shut until the end**, then run a non-existence check on each survivor and
report what was searched. Every mechanism section 15 names was excluded by
fiat, its runners-up included. Section 15's five elements still bind - a
writable correction, a validity argument, a cost in cut scans, a measured
deficit, a one-day falsifier - and so does its last clause: **no
recommendation to build**.

### 16.2 The two structural facts this round rests on

**(S1) The tree prior never reads a predictor value. CONFIRMED, with two
additions.** [`CGMTreePrior::treeLogProbability`](../../src/bartcore/model.hpp) descends through
[`CGMTreePrior::growthProbability`](../../src/bartcore/model.hpp), [`CGMTreePrior::splitVariableLogProbability`](../../src/bartcore/model.hpp) and
[`CGMTreePrior::ruleForVariableLogProbability`](../../src/bartcore/model.hpp), and every leaf of that descent is a
walk over ANCESTOR RULE INDICES intersected with the grid's SHAPE:
[`Tree::splitInterval`](../../src/bartcore/tree.hpp) reads `numCuts`, [`Tree::reachableCategories`](../../src/bartcore/tree.hpp) reads
`categoryCounts` and `hasMissing`, and the branch between them is
[`ColumnStore::splitsBySubset`](../../src/bartcore/data.hpp). No `x` value is touched anywhere, so the
prior is a function of rule indices and grid metadata alone. **Addition
one**: [`Tree::variableAvailable`](../../src/bartcore/tree.hpp) also reads the forest's column mask and the
interaction constraint set, so a prior cancels between two states only when
those agree too - the predicate [`Sampler::installForests`](../../src/bartcore/sampler.hpp) and
[`Sampler::setState`](../../src/bartcore/sampler.hpp) already enforce. **Addition two**: under DART
`splitProbabilities` points at a per-chain Dirichlet draw, so the prior does
not cancel ACROSS chains there, and row 3 below is not valid under `dart`.

**(S2) A tree induces the identical row partition in every chain of one
sampler. CONFIRMED.** `Sampler` owns one `ColumnStore data_` and builds
every chain against it from one forest spec, so codes, cut grid, column
masks and case weights are shared and a chain differs only in trees, leaf
values, sigma and family latents; `installForests` records the consequence
in words, "a same-grid donor installs verbatim". A donor tree therefore
installs without repartition, and the empty-leaf veto, which counts members,
cannot fire in one chain and not another. **What (S2) does not buy is a free
move.** [`Sampler::run`](../../src/bartcore/sampler.hpp) runs each chain's whole burn-and-sample loop to
completion - serially, or one chain per worker thread - so any cross-chain
step needs a per-sweep barrier, a designated RNG so the answer does not
depend on `n.threads`, and the loss of per-chain stream independence.

**A third fact, derived here and used three times below.** [`changeMove`](../../src/bartcore/moves.hpp)
picks its node uniformly among interior nodes ([`Tree::fillNotBottom`](../../src/bartcore/tree.hpp)), so
the depth-0 share of the census's change proposals is `E[1/I]` for `I` the
interior-node count, and by Jensen the mean leaves per NON-STUMP tree is at
least `1 + 1/E[1/I]`: **2.44 at `default`, 2.82 at `lownoise`, 2.34 at
`wide`, 2.59 at `bcf`**, the unconditional mean being lower still. Shipped
trees carry two to three leaves, not sixteen to thirty-two, and no item
below may price itself off a depth-5 tree.

### 16.3 Ranking

Cost is per proposal in cut-scan units, as in
[15.3 Cross-lens ranking](#153-cross-lens-ranking); one full pass over `n`
is about `L` units, and `L` is 2.3 to 2.8 by the fact above.

| # | mechanism | what it is | validity | cost | deficit | novelty | one-day falsifier | refutation |
|---|---|---|---|---|---|---|---|---|
| 1 | Exact draw on the level fibre | add `c_t` to every leaf of tree `t` with `sum_t c_t = 0`; `f` is unchanged exactly, so the conditional on that subspace is the leaf prior alone. Draw `u_t ~ N(-S_t/L_t, tau^2/L_t)` with `S_t` the tree's leaf sum, then `c = u - v (1'u)/(1'v)`, `v_t = tau^2/L_t` | exact Gibbs, acceptance identically 1: no likelihood is evaluated, nothing is tuned | `sum_t L_t` additions, no data pass - under 1/100 of a cut scan | [3.1 Many tree arrangements, one fitted function (ESTABLISHED)](#31-many-tree-arrangements-one-fitted-function-established)'s second clause, quantified as mode F in sec 12.2 and never given a fix | new to BART: the construction is Hastie and Tibshirani's backfitting centering, whose own version CHANGES the target | RUN 2026-09-07: [`C1-frozen-ess.R`](../../benchmarks/R/surfaces/C1-frozen-ess.R) branches the recorded C1 arm (Trig+poly, independent design, 75 trees) with [`storeState`](../../R/dbarts.R) and [`copy`](../../R/dbarts.R), zeroes all four structural probabilities at the last kept draw and separately at the 1250th, five seeds, runs leaf and sigma draws alone for 2500 more sweeps, and reads ESS of `f` at [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s 25 points | derivation CONFIRMED, dimension `m-1` and generically ALL of `ker(Z)` (col(Z) is a sum of `m` subspaces each holding `1_n`, so `rank Z <= sum L_t - (m-1)`). The second-order claim is CONFIRMED and exact: [`ConstantGaussianLeaf::logIntegratedLikelihood`](../../src/bartcore/model.hpp) reduces to `0.5 log(P/(P+Q)) + 0.5 b^2/(s^4 (P+Q))`, so shifting a residual by `c_t` moves every `b` and no birth, death or change ratio is invariant. The "1.7 sweeps at default" timescale is REFUTED as a scale mismatch - `tau = 0.0289` is on the internal response scale ([`GaussianResponse::fitScale`](../../src/bartcore/model.hpp) returns the range) while the quoted `s = 1` is the original one; corrected, the number can only be larger, and it is unmeasured. MEASURED (2026-09-07), by point - see [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial): at the median point the deficit is overwhelmingly structural, ESS 14.9 unfrozen to about 670 of 2500 frozen and lag-1 autocorrelation 0.70 to 0.34-0.40, so leaf space is not the bottleneck there. At the worst point it is: minimum ESS rises only from 1.6 to 4-21 frozen, three orders below 2500, so the leaf Gibbs itself is slow at that coordinate and this row's construction bears on the ranked (minimum-ESS) statistic specifically, not the typical point |
| 2 | Lifted cut displacement | carry a bit `d` per interior node; propose into the ONE-SIDED window `{j in [lo,hi] : 0 < d(j-c) <= 1}` from [`findGoodOrdinalRules`](../../src/bartcore/moves.hpp); keep `d` on acceptance, flip it on rejection and at an interval end | skew-detailed balance against `(T,d) -> (T,-d)`; the mixture stays invariant because [`changeMove`](../../src/bartcore/moves.hpp) redraws the node's `d` from its uniform marginal and [`birthOrDeathMove`](../../src/bartcore/moves.hpp) draws a new node's bit uniformly, whose `1/2` cancels the extended target's `2^-I(T)` exactly | zero scans, one bit per interior node - CHEAPER than the reversible perturb | change's aim (sec 3.3) and sec 3.5's cut axis | new to BART: lifting is canonical general MCMC and no tree sampler carries it | already priced; the open question is whether same-direction RUNS exist - add a signed column to [`cutProbe`](../../src/bartcore/moves.hpp) and count consecutive same-sign improvements | the WINDOW correction is CONFIRMED identically 1 at `w = 1`: the one-sided window holds exactly one index each way and `[lo,hi]` is identical on `T` and `T'`. The claim that the PRIOR ratio is 1 is REFUTED - [2.2 Acceptance, the veto, and the grid](perturb-move.md#22-acceptance-the-veto-and-the-grid) already records that the subtree strictly below the node contributes a non-zero prior difference. The `1/(1-a)` gain reproduces exactly (1.68 / 1.37 / 2.06 / 1.50) but it is this house's own arithmetic for a homogeneous walk at state-independent acceptance, not a literature bound, and the probe measured `a` in a chain where perturb never fires. MEASURED AGAINST: same-direction runs do not exist - consecutive accepted displacements at one node continue in the same direction only 36.6 to 43.1 percent of the time, below the reversible null of 50, with a streak-extension hazard flat in streak length. REFUTED as a source of gain: a lift spends its saving forcing continuation in a direction the chain already prefers to reverse |
| 3 | Same-temperature one-tree exchange between chains | draw chains `A != B` and indices `j, k`, swap the two structures with leaves integrated out and redrawn after; `alpha = L_A(T_{B,k}) L_B(T_{A,j}) / (L_A(T_{A,j}) L_B(T_{B,k}))`, the four collapsed marginals of each tree against each chain's own residual and sigma | ordinary symmetric-proposal MH on the product target `prod_c pi(theta_c)`; the swap is a deterministic involution and the two tree priors cancel by (S1) | 2 to 4 passes over `n`, 3 to 13 percent of a sweep; the real price is (S2)'s barrier | 10.4's per-chain ESS of 2 with a between-chain ratio of 0.5 to 0.8 - the deficit nothing in section 15 targets | FOUND for a single tree, unfound for an ensemble | RUN 2026-09-07: 8 chains in one sampler, 500 burn sweeps then 100 states 5 apart, uniform and matched-index exchanges scored offline against each chain's own collapsed marginal (acceptance table in the paragraph below); the MANDATORY harm clause never gets invoked - acceptance is already under half the kernel's own scored rate at `c1` and four to eight orders of magnitude short of it at `lownoise` - KILLED | algebra CONFIRMED. Prior art is CLOSER than either lens allowed: Rigat's cross-chain CART sampler carries a component swap AND "a whole tree swap between chains", so only the ENSEMBLE instance is unfound. Harm CONFIRMED verbatim - 10.4 says "pooling is what widens the interval", and an exchange is a coupling that drives "between" toward 0. The sizing "one tree's worth = `2^depth` = 16-32x a leaf's worth" is REFUTED by 16.2's third fact: it is 2.3 to 2.8x. MEASURED (2026-09-07): acceptance is 10.1 to 11.7 percent at `c1`'s two seeds, 4.4e-8 to 2e-4 at `lownoise`; it lives entirely on the one- and two-split trees birth/death already reaches and vanishes on the multi-split trees the move is sized for, and it is highest exactly on the couplings that would collapse the between-chain spread pooling converts into 10.4's coverage. KILLED |
| 4 | Pairwise-collapsed split transfer | pick trees `(j,k)`, delete the children of a nog node `v` of `T_j` and install `v`'s rule at an admitting leaf of `T_k`; score by the JOINT collapsed marginal of the pair, which needs only the `L_j x L_k` weight contingency table plus the cached leaf statistics | MH on the block `(T_j, T_k)` with both leaf vectors integrated out; the empty-leaf veto applies at both ends | one pass over `n` for the table (`~L` units) plus one leaf scan; the `(L_j+L_k)^3` Cholesky is noise at `L ~ 2.5` | 3.1's FIRST clause - a split leaves one tree and enters another without either passing through a stump | new outright as far as two search legs reach | generator-only, the shape of [`cutProbe`](../../src/bartcore/moves.hpp): each sweep score one candidate transfer, log the pair-collapsed ratio, RESTORE. A median where change's sits kills it | computability CONFIRMED - the table is one pass over the [`rebuildLeafOf`](../../src/bartcore/chain.hpp) maps, and `r'Wr` cancels because both states share the residual net of the other `m-2` trees. The stated `q` ratio `(N_j M_k)/(N'_k M'_j)` is REFUTED as incomplete: the pair is drawn from the set of trees SHARING a split variable, which the move itself changes, so that selection density does not cancel and both normalizers must be computed. Repairable, and cheap over variable-usage bitmasks. NOT RUN this round: pricing a transfer needs the residual net of the other `m-2` trees and the partner tree's own leaf statistics, which neither `changeMove` nor `birthOrDeathMove` sees, so it needs a dedicated `chain.hpp` block scoring the joint `L_j + L_k` system rather than a hook inside one of the moves this round's other three probes reused |
| 5 | Lifted birth/death | carry a bit `v` per tree; at `+1` the move is a birth, at `-1` a death; flip on rejection. Against [`birthOrDeathMove`](../../src/bartcore/moves.hpp) this deletes both [`probabilityOfBirthStep`](../../src/bartcore/moves.hpp) factors from the transition ratio and replaces the Bernoulli with a read of `v` | row 2's argument with the birth/death involution; the boundary cases (stump, saturated tree) become reflecting and leave the ratio | negative - one bit per tree, one fewer Bernoulli | 3.5, tree size as a random walk | new to BART | recompute `1/(1-a)` per cell from the shipped census; below ~1.3 everywhere, do not write it | code claims CONFIRMED. The price is REFUTED as quoted: the reports read the SUPERSEDED three-move census. Recomputed at the two-move kernel the gain is 1.12 / 1.06 / 1.14 / 1.06 on birth and 1.12 / 1.06 / 1.17 / 1.14 on death - 1.06 to 1.17, a rider on an aim improvement and nothing alone. The lift is on the SIZE axis only; at fixed size the walk over WHICH nodes is untouched. UNCHANGED by this round's probes: they price the cut-displacement axis (row 2) and the nog neighbourhood (15.3 rows 1-3), not the birth/death step itself |
| 6 | Variance-forest-weighted birth leaf | where a variance forest is in the model, choose the birth leaf with weight proportional to the leaf mean of the fitted `s^2(x)`, corrected by the ratio of normalizers | an informed weight that is a deterministic function of the current state; the reverse is the uniform nog draw, computable at `T'`. Fitting a variance forest AS a proposal device off chain history would be adaptation | zero - one accumulator in the leaf-stat pass the sweep already makes | the steepest-change reading of least-favorable directions | new outright | fit a He-Hahn cell with a variance forest and correlate `s^2(x_i)` against `abs(f - fhat)`; near-zero correlation and the weight carries no signal | validity CONFIRMED and the HBART contrast is exact - there `s^2(x)` enters the mean forest only as a precision weight. Two objections. The direction is UNDERDETERMINED: that same precision weight already FLATTENS the mean likelihood exactly where this would aim more proposals. And 10.4 is homoscedastic and fitted without a variance forest, so at the deficit it names this is a MODEL change, not a kernel change |
| 7 | Outer Metropolis on a per-column warp | put a prior on a per-column monotone warp of `x`, install through `$setPredictor`, accept on the fit, roll back on rejection | valid as a model extension; by (S1) the tree prior ratio is exactly 1 whenever the grid SHAPE is held, for every column and even though every cut VALUE moved | one full forest refit per proposal - a whole sweep, so a low-rate outer move | the rooting lock of [10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function) and P6's shelf | new to BART; Snoek's input warping is the GP construction | R-only: apply a fixed monotone warp to P2's confounded columns on the five stuck seeds and see whether the two rootings stop being equiprobable | prior-ratio-1 CONFIRMED and it RECONCILES with 15.3 item 10's verdict (h) rather than extending it: both say pinning the shape gives 1. The claim that a re-derived grid costs only "one `treeLogProbability` call per tree" is REFUTED - [`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp) is a many-to-one nearest-cut remap that collapses starved subtrees, so the proposal is not a bijection and no reverse density exists at all. The bounding negative CONFIRMED: under a re-derived quantile grid codes are ranks, the partition is exactly invariant, and `alpha` collapses to the warp prior |
| 8 | Per-leaf per-variable histograms | cache the binned `(count, sum w, sum wz)` per leaf per variable; a node's histogram is the sum of its children's, a sibling's the parent's minus the other child's | not a proposal - it changes how a number is obtained, no draw | claims a (#leaves)-fold cut in the scan surface's cost | nothing directly; it decides whether the informed constructions of 15.3 are affordable | new to Bayesian tree MCMC (the subtraction trick is LightGBM's docs, not Ke et al.) | microbenchmark against `benchmarks/kernels`, half a day, no engine change | the SIZE of the saving is REFUTED. A scan is over a NODE's members ([`scanOrdinalCuts`](../../src/bartcore/scan.hpp)), so scanning every leaf for every variable already costs `O(n p)` in total, not `O(n p)` per leaf; the histogram wins only on the scan-EVERY-NODE workload and the factor there is depth-fold, about 2 at 16.2's measured tree size, not 16-32. The associativity objection CONFIRMED and unavoidable: a summed or subtracted histogram rounds differently from a member pass, so `equivalence.R` breaks either way |

Novelty column, search record: lifting and non-reversible MCMC against BART,
Bayesian CART, decision trees, treed GP and phylogenetic topology;
population-MCMC crossover and cross-chain tree exchange against BART and
sum-of-trees; joint or blocked updates of a PAIR of trees; centering and
interweaving against BART leaf values; warped, latent and
errors-in-variables covariates under trees; variance forests as proposal
devices; and histogram sufficient statistics against XBART, bartz,
flexBART, stochtree, PyMC-BART, bartMachine and PG-BART.

**Addendum (2026-09-07): row 3 priced at `m = 75`.**
[`C1-cross-chain-probe.R`](../../benchmarks/R/surfaces/C1-cross-chain-probe.R)
runs eight chains in one sampler on `c1` (He-Hahn independent design,
Trig+poly, n = 10000, p = 30, 75 trees, two data seeds) and on the
move census's `lownoise` cell (Friedman, n = 5000, p = 10, 75 trees,
sigma^2 = 0.1): 500 burn sweeps, then 100 states five sweeps apart, 50
uniform (chain pair, tree pair) exchanges plus 25 matched-index
(`k = j`) exchanges scored per state from the four collapsed marginals
in closed form
([`probeLogMarginal`, `probeLogMarginalNumeric`, `probeAssignLeaves`, `probeStructureKeys`](../../benchmarks/R/surfaces/C1-cross-chain-probe.R)
do the scoring, the numerical cross-check and the structure
bookkeeping) - generator-only throughout: nothing draws, nothing is
proposed, the sampler runs its shipped moves untouched.

Acceptance of min(1, alpha) (`excl.` drops exchanges between
identical structures, `ident` is that share, all proportions):

    cell      design    n     accept   excl.    medianLogAlpha  ident
    c1 seed1  uniform  5000   0.1074   0.1034      -53.64       0.0044
    c1 seed1  matched  2500   0.1136   0.1111      -47.86       0.0028
    c1 seed2  uniform  5000   0.1037   0.1012      -49.87       0.0028
    c1 seed2  matched  2500   0.1169   0.1134      -49.33       0.0040
    lownoise  uniform  5000   2.0e-4   4.4e-8     -1374.00       0.0002
    lownoise  matched  2500   4.1e-8   4.1e-8     -1394.00       0.0000

By leaf-count pair (`c1` seed 1, uniform; flat across the run's four
sweep blocks, and only smaller with more splits on either side):

    pair    accept
    1-2     0.331
    1-3     0.170
    2-2     0.164
    2-3     0.073
    3-3     0.031
    4-4     0.028
    5+-5+   4e-22

Shared structure, mean over states and chains:

    cell      shared  non-stump  stump   meanLeaves
    c1 seed1   0.290    0.248    0.057     2.59
    c1 seed2   0.283    0.241    0.055     2.55
    lownoise   0.056    0.054    0.002     4.06

consistent with 16.2's third fact and with 6.1's fourth addendum,
where `lownoise` also carries the deepest trees of any single-forest
cell. Checks: the closed form against
numerical integration on all 15 leaves scored agrees to max |diff|
2.2e-11; the leaf assignment reproduces `getTrees`' own `n` column for
all 600 trees of a state; summed per-tree internal leaf values
reproduce the fitted values to 2.2e-15 to 6.6e-15, pinning the
internal scale, `tau = 0.028868`; and one whole-tree brute-force log
alpha (90.21401) reproduces the closed-form value digit for digit.

The positive tail is real - P(log alpha > 20) is 0.0052 (`c1` seed 1,
uniform) against a stationary bound of 2.06e-9 - and the brute-force
check's own donor and acceptor chains differ in internal sigma by
about 5 percent (0.0846 against 0.079-0.081 for the other six), the
same likelihood-side signal that the chains sit apart.

Reading: not worth building. At `c1` the exchange accepts at under
half the kernel's own scored acceptance (24.9 percent, 6.1's fourth
addendum) in the same cell; the acceptance lives on the one- and
two-split trees birth/death already reaches and vanishes on the
multi-split trees that are the move's own justification; at
`lownoise`, the regime where structure freezing is established, it is
unavailable; and the accepted exchanges are exactly the couplings
that would collapse the between-chain spread that pooling converts
into
[10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)'s
coverage.

### 16.4 Merged view across the two rounds

Generator-only, on the census build and nothing else: 15.3 items 1 to 3's
nog probe (weight entropy, `P(incumbent)`, the nog share); row 4's
pair-collapsed ratio, the same shape; row 2's signed run-length column; and
the leaves-per-tree count itself, which is logged nowhere, which 16.2 could
only bound by Jensen, and which prices rows 3, 4 and 8 at once. R-only, no
engine change: row 1's frozen-structure ESS, which reorders every item in
all six reports and which nothing on disk answers; row 3's per-tree fit
pre-check across chains; row 7's fixed warp on P2. Needing code before any
evidence exists: rows 2, 3 and 5 as kernels (row 1's kernel: slice 1 landed
cbe80534 and the benefit run KILLED it 2026-09-07, the step staying at
default off - see
[level: an exact Gibbs draw on the level fibre](level-fibre.md#level-an-exact-gibbs-draw-on-the-level-fibre)),
row 8 as a cost model, and row 4 beyond its probe. Nothing here is
scheduled. Of the four, three have
now run - the nog probe, the signed run-length column and the
leaves-per-tree count (sec 6.1's third 2026-09-07 addendum) - leaving row
4's pair-collapsed probe as the one generator-only item still not built.

### 16.5 Discarded across both reports

- **Aiming at the maximum-posterior-variance directions of `f`**: they lie
  in `col(Z)`, which the leaf Gibbs already draws exactly, which is why the
  round read "least-favorable" as its complement.
- **Full joint leaf draw over all `m` trees**: `m(m-1)/2` Gram increments
  per row plus a `(sum L_t)^3` Cholesky, 30-100x a sweep. Rows 1 and 4 are
  the affordable members of that family.
- **Centering each tree's fit on the fly**, Hastie and Tibshirani's own fix:
  BART's leaf prior is proper, so a constraint moves the target; row 1 is
  the target-preserving form.
- **Prior-only tempering as the fibre temperature**: it flattens exactly
  `ker(Z)` and nothing else - the clean answer to "flatten only the fibre" -
  but it is tempering, dormant for cause, and row 1 makes it unnecessary
  wherever the fibre-restricted conditional is closed form.
- **Pseudo-prior dormant subtrees** (Carlin-Chib in the product space):
  valid, and one accept for a `k`-level size change, but dead on arrival at
  the PRIOR pseudo-prior - a prior-drawn rule is what `changeMove` draws and
  the census prices it at -62 to -143 log units, so a `k`-rule subtree is
  about `k` times worse. Its repair needs a fitted pseudo-prior.
- **Informed choice of WHICH tree gets the move**: the correction is right,
  the cost REFUTED - every tree's residual moves when any tree's fit does,
  so the normalizer costs a pass per tree, and reading the sweep's own
  running caches instead is history-dependent, hence adaptation.
- **Lifted variable choice**: `p` split variables carry no order, and
  ordering them by DART's `s` is the adaptation trap.
- **Whole-state exchange between same-temperature chains**: a relabelling of
  exchangeable chains, acceptance 1, effect nil. (What makes Rigat's
  whole-tree swap real is that his chains differ in PROPOSAL.)
- **Block exchange of `k > 1` trees, and subtree crossover between chains**:
  the first reshuffles the representation without moving the fit; the second
  loses the exact prior cancellation and pays all of rotation's revalidation.
- **Rank-one subtree transfer WITHIN one tree**: reversibility needs a merge
  enumeration, which row 4 avoids by checking the destination's constraints.
- **A second copy of the family latents to aim while the first scores**: `z`
  is already conditionally sufficient, so the copy carries strictly less.
- **Per-observation latent coordinates in the design matrix**: `n` free
  parameters against a leaf that charges nothing for a singleton - an
  identifiability failure, not a proposal.
- **Discrete zig-zag on tree size** (row 5 at more machinery); **a global
  shift on (forest, probit latents)** as an interweaving step (null for
  gaussian, not shift-invariant for probit); **deterministic
  prune-and-remember** (history-dependent; its valid form is the discarded
  pseudo-prior above).

### 16.6 What this could not settle

- ~~**Whether the leaf-value half or the structural half carries the
  He-Hahn ESS.** Row 1's frozen-structure experiment answers it in a day
  and would reorder both rounds; section 3.2's standing datum ("leaf
  values converge in a handful of sweeps") points the other way from mode
  F's timescale.~~ SETTLED (2026-09-07): both, at different points - the
  median point's deficit is overwhelmingly structural and the worst
  point's is the leaf half; see
  [10.4 C1, the He and Hahn factorial](benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial).
- **The corrected magnitude of that timescale.** The formula is 12.2 B5's
  and stands; its instantiation was on the wrong scale, and the internal
  `sigma` of the census cells is recorded nowhere.
- **Whether `ker(Z)` exceeds the `m-1` level directions.** `m-1` is now
  known to be the generic floor, reached exactly unless two trees'
  partitions coincide or nest, so the excess is a count of coincidences
  rather than the full Gram the geometry lens assumed.
- ~~**Whether exchange acceptance survives at `m = 75`**, the premise
  being measured and the acceptance not.~~ SETTLED (2026-09-07): it does
  not, at either cell tried - see the addendum after
  [16.3 Ranking](#163-ranking)'s table.
- **Whether same-direction cut RUNS exist**, without which row 2 collapses
  onto the reversible perturb at no loss and no gain; and whether tree size
  is a slow coordinate at all, on which row 5's whole value rests.
- **The leaves-per-tree distribution**, bounded here by Jensen off a table
  built for another purpose. Three rows are priced against that bound.

### 16.7 Provenance

```
inputs        two independent lens surveys under a novelty gate, 2026-09-07
              (function-space geometry; constructions from general MCMC and
              from dbarts' own affordances), then a refutation and synthesis
              pass over both records
base          15c908cb
scope         ranking and evidence only. No code, no default, no schedule.
```

Verified in THIS pass; a claim carried from the round-two records without a
re-fetch is named at the end.

| # | Source | Verified | Where |
|---|---|---|---|
| 1 | Rigat, Parallel hierarchical sampling | "the equilibrium distributions of all chains is the same but the proposal distribution used to update each chain is different"; "since all temperatures have value 1, the Metropolis swap acceptance ratio (5) is equal to one"; "the cross-chains version of the insert, graft and change transitions, swapping the elements of the tree structure"; and "The second class of cross-chains transitions includes a whole tree swap between chains". A SINGLE-tree CART sampler: "The leaves are the final nodes of a single-rooted binary partition of the covariates space". At that URL the author is Rigat alone; the Rigat and Mira (2012) journal version was not fetched | arXiv 0812.1484 |
| 2 | Liang and Wong, Evolutionary Monte Carlo | "The population is updated by mutation (Metropolis update), crossover (partial state swapping) and exchange operators (full state swapping)" | statistica.sinica.edu.tw A10n21 |
| 3 | Gagnon and Doucet, nonreversible jump | "By lifting this model indicator variable, we obtain non-reversible jump algorithms"; restricted to "nested models, a class of models for which the model indicator variable is an ordinal random variable". **No acceptance-rate bound of any kind appears**; the paper's quantitative claim is a scaling limit, time accelerated "by a factor of only sqrt(n) ... comparatively to n", and its Remark 1 says the Diaconis et al. orders are "K^2_max and K^2_max log Kmax steps ... for the non-reversible and reversible ... samplers" - a log factor, on that target | arXiv 1911.01340 |
| 4 | Gagnon and Maire | "the asymptotic variances cannot increase by a factor of more than 2, regardless of the target distribution" - a guarantee against DEGRADATION, not of improvement, and hedged in the original by "essentially" | arXiv 2405.15952 |
| 5 | Diaconis, Holmes and Neal | "nonreversibility can indeed lead to improvements over the diffusive behavior of simple Markov chain sampling schemes"; "order n steps are necessary and sufficient for convergence in total variation distance of the non-reversible walk" against `n^2` reversible - and the paper flags this as "one of the few natural instances where total variation and chi^2 relaxation times differ" | projecteuclid.org 10.1214/aoap/1019487508 |
| 6 | Hastie and Tibshirani, Bayesian Backfitting | "in the Bayesian backfitting algorithm, we have to center the fits after smoothing and generation", because the constraints "are necessary to ensure that the posterior distribution of alpha and the fj is not singular" - identifiability, not an improper prior | projecteuclid.org 10.1214/ss/1009212815 |
| 7 | LightGBM feature docs | "It then can get histograms of its neighbor by histogram subtraction with small cost (O(#bins))". The NeurIPS paper of Ke et al. does not contain it | lightgbm.readthedocs.io Features |
| 8 | Pratola, Chipman, George and McCulloch, HBART | the mean forest's likelihood differs from BART only by "replacing a scalar variance s^2 with a vector variance s^2(x_i)", and "The draws of T_j and T_0j are done using Metropolis-Hastings steps as in Chipman et al. (2010) and Pratola (2016)" - the variance forest is a precision weight, never a proposal | arXiv 1709.07542 |

CARRIED, NOT RE-FETCHED in this pass and load-bearing for nothing above:
Park and van Dyk on partially collapsed Gibbs; Jasra, Stephens and Holmes,
and Drugan and Thierens, on temperature-free crossover; Turitsyn, Chertkov
and Vucelja; Koskela's zig-zag; Mohammadi, Pratola and Kaptein's
continuous-time birth-death; Snoek et al. on input warping; CGM 2010 on
monotone invariance; He, Yalov and Hahn on XBART's cumulative sums; and
Carlin and Chib, on which 16.5's pseudo-prior verdict is provisional.
