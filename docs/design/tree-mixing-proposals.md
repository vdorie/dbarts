# tree-mixing-proposals: how the BART posterior is sticky, and what could move it

Status: COMPLETE (survey with an adjudicated evidence base), 2026-08-09;
**section 12 is an addendum from a second research cycle, 2026-08-10**,
which refutes one of this document's recorded inferences (erratum in
sec 5.4), amends sec 3.1, and adjudicates fourteen new candidates.
No work proposed, nothing scheduled, no source touched. TODO
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

`metropolisJumpForTree` (`moves.hpp:844-865`) draws one uniform per tree
per sweep and dispatches to exactly one of three kernels:

```
u < birthOrDeathProbability          -> birthOrDeathMove   (moves.hpp:206)
u < birthOrDeath + swapProbability   -> swapMove           (moves.hpp:734)
else                                 -> changeMove         (moves.hpp:475)
```

Shipped mixture `birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5`
(`R/model.R:5-9`, `R/dbarts.R:365`, engine defaults `chain.hpp:59-62`).
`StepType` is `{birth, death, swap, change}` (`moves.hpp:838`): the
"four-move set" is three kernels with four labels. The sweep is
Gauss-Seidel over trees; the variance forest runs the identical kernel
(`chain.hpp:4232-4248`).

- **Birth** picks a leaf uniformly among those with a usable variable and
  draws the new splitting rule from the prior, so the rule's prior density
  cancels in the acceptance ratio. One split, at the fringe.
- **Death** removes one split, at the fringe, among nodes with no
  grandchildren.
- **Swap** exchanges the rules of a parent and one child. It preserves the
  tree's shape, is symmetric, and carries no proposal correction
  (`moves.hpp:825`). It cannot lift a rule more than one level.
- **Change** picks uniformly among all non-leaf nodes *including the root*,
  redraws the split variable from the prior, then the cut point uniformly
  over the descendant-valid set. **The entire skeleton below the node is
  held fixed** (`moves.hpp:595-596`) and every observation is rerouted
  through it, with a hard veto if any descendant leaf empties
  (`moves.hpp:79`).

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
   `split.probs = 1/num.vars`, `R/model.R:293-297`). The machinery for such a
   move already exists and is already gated (section 5.1).

Adjacent machinery already landed, relevant to cost: `scan.hpp:105`
`scanOrdinalCuts` (leaf-templated full-cut scan); `grow.hpp:179`
`growTreeFromRoot`; `chain.hpp:1967` `growForestFromRoot` (opt-in
`n.grow.sweeps`, init only, with a reset/regrow/rebuild/redraw loop at
`chain.hpp:2009-2022`); `tree.hpp:999-1021` `SubtreeSnapshot` (restores node
*contents* for a fixed set of node ids - it cannot undo a shape change);
`tree.hpp:1065`/`:1257` `collapseEmptyNodes` / `collapseSubtreeToLeaf`.

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
split-count feedback loop (`chain.hpp:1569-1574` recomputes split counts
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
itself (`model.hpp:178-181`). As the residual standard deviation sigma falls,
or as a multiplicative forest weight rises, that exponent's magnitude grows
and every proposal that is not an improvement is rejected outright. Leaf
values converge in a handful of sweeps; what remains is a partition-shape
misfit that sigma has to absorb, and absorbing it keeps sigma high, which
is the only thing keeping acceptance non-zero.

**Corrupts.** Sigma's own posterior, interval coverage for the fitted
function, and every structural readout. Point estimates survive.

**Measured, in this house** (`docs/plans/bcf-sigma-residual.md`, measured
at `bartcore 6944811`). In the causal-forest sampler's prior tail, where a
scale parameter `a` multiplies one forest's contribution and the engine
hands that forest weight `w_i * a^2` (`combiner.hpp:877-894`):

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
(docs/plans/composition-mixing-probe.md, run to verdicts the day after this
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
  forest. This is stan4bart's shape and rbart_vi's random intercept is its
  degenerate case.
- **Inner composition**: the parametric part lives in the *leaves* - dbarts'
  landed linear leaves (`linear-leaves.md`) and GP leaves (`gp-leaves.md`);
  MOTR-BART is the external name for the linear case. **The inner variant
  has no cross-block ridge in the structural move at all**, verified: the
  linear leaf's `logIntegratedLikelihoodForNode` (`model.hpp:1090-1116`)
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
(`src/init.cpp:642`): `dbarts_sampler_run(sampler.bartSampler, 0, 1, ...)`
runs exactly **one** BART sweep, the fit is copied into `stanOffset`
(`:654`), the parametric sampler is given it (`setOffset`), and WALNUTS
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
| outer composition, conjugate block | `rbart_vi` / in-engine `GroupedResponse` | shipped |
| outer composition, arbitrary block, user-driven | `dbartsSampler$setOffset` (`R/dbarts.R:1316`), and `dbarts_sampler_setOffset` in the shipped C API (`dbarts.h:755`) | shipped, supported |
| inner composition | `node.prior = linear(columns)` / `gp(columns)` (`R/model.R:37-40`, `:51`) | shipped |

The probes exist too: `benchmarks/R/grouped-mixing.R` (the autocorrelation
harness), `benchmarks/R/forest-ranef-collapse-proto.R` (the isolated ridge
surrogate), the XOR scenario from the grow-from-root study, the pooled
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
  (`brtmoves.cpp:62-73`), sending 90% of those attempts to perturb
  (`pchgv = 0.1` shipped in both the R and Python wrappers). At a typical
  75-tree interior-node count that is roughly 3-6 perturb attempts per tree
  per sweep - 20 to 40 times the throughput the survey proposed testing.
- **Window width is untuned and the two published settings are 8.5x
  apart.** OpenBT's window is 10% of the admissible range centred on the
  current cut (`brtmoves.cpp:218-225`), with the boundary asymmetry
  corrected by a simple count ratio at line 266. Pratola's own paper fixes
  it at **85%** ("we have fixed alpha to cover 85% of the range defined in
  the interval (2)"), and his headline numbers were produced at 85%. No
  comparison between the two exists anywhere.
- The nearest disconfirming datum is Tan et al.'s Experiment 7 (section 3.4
  above): restricting away change and swap changed nothing measurable on
  their battery. Perturb is a restriction of the change move.

**Why it is first anyway: the fit to dbarts is exceptional.**

- **The Metropolis-Hastings correction is provably 1** for the uniform
  form, verified by construction at `moves.hpp:566-580`: with the new
  variable equal to the old one on an ordinal column, the forward and
  reverse interval and valid-set counts are produced by the *same two
  functions on the same unmodified tree*, and both ignore the node's own
  rule, so `log(rI) - log(fI) + log(fV) - log(rV) = 0` identically - not
  approximately. Both-categorical falls through all three branches at lines
  567-581 and leaves the correction at its `0.0` initializer, which is also
  correct. `change-move-balance.md` records the same: "a same-variable or
  equal-cut ordinal change gives correction 1".
- For the local-window form the correction is `|W(c)| / |W(c')|`, the ratio
  of the two clipped window sizes - still exact, still computed by the same
  two counters, just not identically 1.
- The subtree-prior term below the node must still be computed in general
  (a same-variable cut can flip a descendant's variable availability, which
  changes `growthProbability` and `splitVariableLogProbability` at
  `model.hpp:2129-2145`), and `changeMove` already computes it.
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
(`brtmoves.cpp:305`). Two honest limits: **the number of trees is never
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
  (`grow.hpp:179`) used when this was written. That precedent has since
  been spent: the builder scans categoricals too
  (`scanCategoricalPartitions`, `grow.hpp:224-263`, `scan.hpp:337`, the
  winning rule built at `grow.hpp:332-335`), so an ordinal-only rotation
  would now be scoping narrower than the builder rather than matching it.
- *Interaction constraints.* A rotation lifts one variable above another -
  the same class of break that `swapMove` already guards with
  `tree.interactionSubtreeIsValid` (`moves.hpp:804`). That guard is
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
thread-parallel chain layout (`R/bart.R:669-670`), and destroys the
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
(`scan.hpp:105`) yields the collapsed marginal likelihood for *every*
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
LANDED (`docs/plans/grow-from-root-categorical-scan.md`; the TODO door is
recorded closed at `TODO:302`). **The half worth as much as the
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
   combiner hands a forest weight `w_i * a^2` (`combiner.hpp:877-894`), so
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
per-node candidate weights (`grow.hpp:195-247`), so the forward density is
nearly free; the reverse density requires replaying the builder's candidate
assembly along the *current* tree's construction path, the same cost again.
The reset/regrow/rebuild/redraw loop already exists (`chain.hpp:2009-2022`)
and lacks only the acceptance filter. Reachability limited it to
ordinal-only forests until the scheduled categorical scan landed; that scan
has landed, and `growTreeFromRoot` (`grow.hpp:179`) now emits categorical
candidates of its own (`grow.hpp:224-263`), so the limit is gone.

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
integrated likelihood (`grow.hpp:187-190`, `:195-247`). The governing object
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
  ambiguous." dbarts defaults to `n.trees = 75L` (`R/bart.R:666`), below
  BART's classic 200. **Carry the caveat with it**: more trees dampens
  R-hat partly *because* the ensemble self-averages structural labels
  harder, so an improved R-hat at larger m is not by itself evidence that
  tree-space mixing improved.
- **More chains.** Ronen et al.'s own recommendation is to "increase the
  number of chains with the number of data points"; dbarts defaults to
  `n.chains = 4L` (`R/bart.R:669`).

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

Four arms, matched seeds, paired (the `grouped-mixing.R` idiom: data seed
`BASE_SEED + s`, sampler seed `s`, shared `s` across a pair).

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
| Stage 2 harness | ~400 lines | reuses the `grouped-mixing.R` seed idiom |
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
  fixed-skeleton design working as specified (`moves.hpp:595-596`), and the
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
  exactly 1** under shipped machinery (`moves.hpp:566-580`). This is the
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
| 7 | Tan, Ronen, Saarinen, Yu, arXiv 2406.19958v2 | Full text: eq 10, eq 13, Thm 7.3 + its "is not the posterior" paragraph, Prop 7.4, the T pincer, sec 9 findings, Appendix L.6. **Plus the released code**: `bart-comp-efficiency` @ 0589240 (`runner.py:40`, `utils.py:274-276`) and `bart-playground` @ e686d34 (`samplers.py:113-118`, `:268-273`) | arXiv 2406.19958 ; github.com/yanshuotan/bart-comp-efficiency ; github.com/yanshuotan/bart-playground |
| 8 | Zanella, JASA 115(530):852-865, 2020 | Abstract + the balancing-function conditions and the no-dominance statement | arXiv 1711.07424 |
| 9 | Angelopoulos, Cussens, ICML 2005 | Full text: sec 5 "tempering (aka Metropolis-coupled MCMC)", the 4-chain ladder, swap-every-iteration, cold-chain-only collection, Table 6 (3 seeds, 15/16), PIMA accuracy, the CGM-1998 block quote | icml.cc/Conferences/2005/proceedings/papers/003_Tempering_AngelopoulosCussens.pdf |
| 10 | Deshpande, flexBART, arXiv 2211.04459 | Full text: sec 2.2 informed-proposal objection, appendix B2/B3 | arXiv 2211.04459 |
| 11 | Zhang, Huelsenbeck, Ronquist, Syst Biol 69(5):1016-1032, 2020 | Record + verbatim abstract (order-of-magnitude faster convergence; dataset-dependence caveat) | academic.oup.com/sysbio/article/69/5/1016/5716338 |
| 12 | Gramacy, Taddy, JSS 33(6), 2010 | Record + abstract; importance tempering, `itemps = NULL` by default | jstatsoft.org/article/view/v033i06 |
| 13 | OpenBT | Source @ main: `brt.cpp:1278-1299` (`drawvec`), `brtmoves.cpp:62-73` (`pertcv`, every interior node), `:218-225` (10% window), `:266` (window correction), `:305` (live `rot`); `brt.cpp:1952-2490` is one unclosed block comment, so the `rot` at 2198 is dead code; `pchgv = 0.1` in `misc/openbt.R:53-54` and `misc/openbt.py:29-30` | github.com/jcyannotty/OpenBT |
| 14 | He, Hahn, arXiv 2002.03375 | Relied on through `grow-from-root-default.md`'s verified record (sec 5, Table 4) | arXiv 2002.03375 |
| 15 | Chipman, George, McCulloch, JASA 93(443), 1998 | Read only as a verbatim block quotation inside #9 (JASA text paywalled). Cite accordingly. | (quoted in #9) |

**Composition block (section 4.1), sources and status.** Section 4.1 was
added late and its evidence base is deliberately weighted toward what could
be verified directly rather than toward citations.

| source | status |
|---|---|
| `docs/design/forest-ranef-interweaving.md` sec 0, 2, 5, 6, 9 | Read in full at `d3cb94b`. The 56.1 / 9.3 / 114.6 prototype table, the with-f/no-f attribution, the "no cheap ASIS/PX" structural argument, and section 9's authoritative corrections are all quoted from it directly. This is the load-bearing evidence for the hazard. |
| `model.hpp:1090-1116` (linear leaf integrated likelihood) | Read. Confirms the inner variant marginalizes leaf coefficients out of the structural score. |
| `R/dbarts.R:1316`, `dbarts.h:421`, `R/model.R:37-40,51`, `R/rbart.R:947` | Read. Confirms the composition surface is public on both the R and C sides, and that rbart_vi's loop is 1:1 alternation. |
| `inst/common/friedmanData.R` | Read. Confirms the probe DGP decomposes into one interaction plus three separable terms. |
| stan4bart `src/init.cpp:642,654`, `docs/design/walnuts.md` | Read in the live tree (0.0.14 installed). Confirms the 1:1 two-block Gibbs alternation and that no mixing diagnostic is reported. |
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
              docs/plans/grow-from-root-default-study.md
              docs/plans/bcf-sigma-residual.md sec 1-4
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
   `model.hpp:2117-2131` (`base = 0.95`, `power = 2`, root at depth 0),
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
(`model.hpp:2129-2132`) returns 0 **only** when
`!tree.hasAnyAvailableVariable(...)`; otherwise `base/(1+depth)^power`,
strictly positive at every finite depth. The memo states this itself
elsewhere and then assumes its negation.

*The correct general formula*, re-derived here from `grow.hpp:179-247`
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
Gaussian leaf exactly as `model.hpp:165-183`, CGM prior exactly as
`model.hpp:2117-2145`, builder weights exactly as `grow.hpp:179-247`,
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
  resolvedNodeScale(options.nodeScale, options.priorScale) / sqrt(forest.numTrees)` (`chain.hpp:648-650`) uses the
  **total** tree count, and `installBlockMasks` (`chain.hpp:5116-5138`,
  called at `:722`) installs per-tree column masks and nothing else - no
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
inside the sweep (`chain.hpp:4500-4518`), so it depends on the other trees
**only through their total**. Apportionment can therefore reach structure
only through `f_j`'s own drift. The pre-registerable claim is: *tree `j`'s
structural acceptance pattern is autocorrelated at the `f_j` timescale, and
that autocorrelation is what apportionment stickiness costs.* Same
instrumentation, falsifiable in both directions.

**Three further corrections to the census design**, all adopted.

- **No matching or alignment step is needed.** The sweep updates tree `t`
  in place (`chain.hpp:1457`) and nothing in `chain.hpp` shuffles, permutes
  or relabels trees (verified by search). The *posterior* is
  label-exchangeable; the *chain* never exercises the symmetry.
- **`getTrees` is necessary but not sufficient.** It returns "a data.frame
  containing the internal state of the trees" (`R/dbarts.R:1939`) - flat
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
(`chain.hpp:1968-1970`), so the **monotone constant leaf is in the
builder's scope**. `MonotoneConstantGaussianLeaf::logIntegratedLikelihood`
delegates to the *unconstrained* `ConstantGaussianLeaf` and is documented
"never on the constrained hot path" (`model.hpp:536-542`). Meanwhile the
structural target under that leaf is not a leaf-marginalized posterior at
all: `logLikelihoodForBranch` dispatches `ParamScoringLeafModel` to
`leaf.logLikelihoodForBranchWithParams(...)`, reading frozen neighbour
parameters (`moves.hpp:91-94`). So for a monotone forest `q` is computable
but `pi` is **not** `p(T) prod_leaves L(leaf)`, and "every term in the
acceptance ratio has a shipped implementation" is false there. Monotone
forests must be scoped out of any regrow v1 explicitly.

**The predicate itself.** A hand-written list (categorical split | zero-
growth node | missing-only side | non-constant leaf | ...) will rot. The
equivalent, always-correct predicate is a property of the replay the move
needs anyway: **refuse iff the reverse replay returns `-inf`**. It is
decidable from the incumbent alone because the builder's support is
residual-independent - occupancy depends on member counts and availability
on the cut grid and masks, never on `y` (`scan.hpp:105-140`,
`tree.hpp:572-610`) - it is exactly the support indicator, and it costs
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
`priorPrecision = (k/scale)^2` (`model.hpp:174`) and
`drawFromPrior = (scale/k) * z` (`model.hpp:198`). So

```
tau = nodeScale / (k sqrt(m)),   timescale ~ n_leaf nodeScale^2 / (m k^2 s^2)
```

a factor `k^2 = 4` smaller at the default `k = 2`. With `nodeScale = 0.5`
(`src/R_interface_bartcore.cpp:280`), `tau = 0.0289` at `m = 75`, not
`0.0577`. `k` is the shipped knob that enters the prediction **squared**
while `n`, `m` and `sigma` enter linearly, and it is itself sampled when
`updateK` is on (`chain.hpp:1565-1567`), which the caricature assumes
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
  (`model.hpp:2129-2132`); and banning changes `numAvailable` and hence
  `P(var)`, but both directions compute their own normalizer, which is all
  MH needs (`grow.hpp:207-222` mirrors `model.hpp:2135-2145`). **Condition
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
     *all* members (`grow.hpp:187-190`) while `scanOrdinalCuts` skips
     `naCode` members outright (`scan.hpp:128`). This does not break
     exactness - `q` is whatever the builder's realized normalized weights
     are, and that is what an accumulator records - but it does break the
     memo's stated reason, and it is a structural `q`/`pi` mismatch rather
     than a data-driven one.
     SUPERSEDED PREMISE: the scan no longer skips them - `scan.hpp:128`
     accumulates the `naCode` rows into a missing bin that every candidate
     adds to one of its two children, and `scan.hpp:100-103` states that
     the scan's scores therefore agree with the leaf statistics
     `tree.birth` caches, missing rows included (`TODO:312`,
     `ordinal-scan-missing-rows`, DISCHARGED), so the structural mismatch
     concluded from the old premise no longer holds.
  2. **Missing-capable columns are halved twice.** `logCut` already
     subtracts `log 2` for `data.hasMissing[j]` (`grow.hpp:290-291`), matching
     `ruleForVariableLogProbability`'s `+log 2` for *one* rule
     (`model.hpp:2179-2186`); the builder then draws the direction coin
     separately (`grow.hpp:340-342`). Since the candidate stands for the
     *pair* of direction-rules, its weight carries one rule's prior mass,
     so `q` under-weights splits on missing-capable columns by 2x per
     split. **This is already known and scheduled**: it is exactly the
     pre-registered two-arm falsifier in
     `docs/plans/grow-from-root-categorical-scan.md` sec S0, with an OPEN
     VD FORK on whether to change the shipped weight. What a regrow adds is
     stakes: for a warm start it is a start-quality bias, for a regrow it
     becomes a systematic term in the importance weight.
     SUPERSEDED PREMISE: that `log 2` is gated on `routesMissing`
     (`grow.hpp:291`) - whether the scan emitted both directions for THIS
     node's members (`grow.hpp:288`), each direction then its own candidate -
     and not on `data.hasMissing[j]`; the direction coin (`grow.hpp:341-342`)
     fires only where a candidate did NOT already name its direction. The
     double-halving concluded from the old premise no longer holds.
  3. The reverse replay must scan **every node of the incumbent, leaves
     included** - a leaf contributes `(1-g) L(u) / Z_u` and `Z_u` needs the
     full candidate assembly there. Cost is `2 x (nodes)` candidate
     assemblies, not `2 x (levels)`.
- **The categorical support gap is confirmed, and the `P(var)` accounting
  is right.** `grow.hpp:224` skips `ColumnType::categorical` outright
  while `numAvailable` counts categoricals even though they generate no
  candidates - which is *correct* for matching `model.hpp:2135-2145`'s
  `-log(numAvailable)`, since both sides count the same set and the missing
  mass is absorbed by `Z_node`. Worth recording because it is the one place
  the builder's `P(var)` and the prior's agree exactly and it is not
  obvious from either file alone.
  SUPERSEDED PREMISE: `grow.hpp:224` now ENTERS the categorical branch and
  emits one candidate per admissible partition
  (`scanCategoricalPartitions`, `grow.hpp:224-263`, `scan.hpp:337`), so
  categoricals do generate candidates and the support gap concluded from
  the old premise is closed.
- **Construction 2(a)'s exactness needs a DART scope condition.** A
  cross-tree repulsion prior has a normalizing constant over the joint tree
  configuration space; it cancels in the MH ratio only while everything it
  depends on is fixed. dbarts ships DART, which resamples the split
  probabilities every sweep from realized split counts
  (`chain.hpp:1569-1574`), and that constant depends on those
  probabilities. So fixed-mask `blocks()`-style constraints stay exact, but
  **DART plus any cross-tree split-usage repulsion is doubly intractable in
  the `s` draw** - which lands on the memo's own "cheap and exact" verdict
  for the cheapest repulsion variant, since that variant's cheapness came
  from reusing DART's counts. The `k` hyperprior
  (`chain.hpp:1565-1567`) is unaffected; only variable-selection parameters
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
  member (`combiner.hpp:153`), so per-block `base`/`power` needs a per-tree
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
   (`model.hpp:174`, `chain.hpp:648-650`), fixed and sigma-free. That is
   BART's design, not an oversight. Integrating sigma against a
   fixed-variance leaf prior gives no closed form in the same sufficient
   statistics.
2. *Even granting the algebra, the marginal cannot move the number it
   targets.* Under a scaled-inverse-chi-squared sigma prior with
   `nu = sigmaDf` (3 by default, `chain.hpp:256`), marginalizing replaces
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
reads weights only through `(sum w, sum w z)` (`model.hpp:165-183`), so the
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
score at `moves.hpp:66-95` cannot serve and the cost is a joint solve
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
composition probe (`docs/plans/composition-mixing-probe.md`, RE-REGISTERED
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
in-repo docs  docs/plans/grow-from-root-categorical-scan.md (S0),
              docs/plans/composition-mixing-probe.md,
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
