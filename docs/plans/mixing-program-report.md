# The mixing program: what was measured, what was found, what to do

Status: REPORT, 2026-09-08. A record; no default is changed here.

## 1. Summary

BART's tree sampler moves by small local edits, and on hard problems it
settles into one arrangement of trees and stays there. The visible symptom
is under-coverage of the fitted function; the underlying quantity is
effective sample size, which on the one core problem with a published BART
coverage number sits at two draws out of two thousand five hundred. The
package now has an evaluation battery with a pre-registered accept rule,
and five candidate remedies were measured against it. Four failed: the swap
proposal earns nothing at production forest sizes, a same-variable cut move
earns nothing anywhere, an exchange of trees between chains is priced out,
and an exact draw on the forest's level fibre helps only when the trees are
held fixed. One succeeded. Replacing the Metropolis change proposal at a
node whose two children are both leaves with an exact draw from that rule's
full conditional raises the effective sample size on the core cell by two
to five times the pre-registered bar, at both seed blocks and on both mean
functions, and a restricted version that redraws only the cut keeps about
half that gain for a thirtieth of the arithmetic. The recommendation is to
consider that restricted draw for a nonzero default share, at a dose the
dose-response run has still to name, and to keep every other kernel at
zero.

## 2. The question

BART's intervals do not always cover, and colleagues report it in settings
where the fit looks fine. The literature has one clean number: He and
Hahn's factorial study reports 95 percent pointwise coverage of the true
mean function of 0.73 to 0.74 at ten thousand rows, thirty predictors and
moderate noise. That setting is this battery's core cell, called C1 below.
dbarts reproduces the deficit in direction and about half in size - 0.82
against a nominal 0.95 at the setting closest to the paper's, with
intervals 15 to 25 percent longer at the same point accuracy, so this
package is better calibrated than the published BART and still short of
nominal
([10.4 C1, the He and Hahn factorial](../design/benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)).

Two levers could close that gap: more trees, or better mixing. The
tree-count lever works and is cheap - raising the forest from 75 trees to
200 buys about ten points of coverage on both of C1's mean functions and
costs nothing in error. It was nevertheless parked, because the same run
showed the coverage arithmetic is not a prior question: four short chains
pooled read 0.96 where one long chain of the same total cost reads 0.90,
and the chains disagree enough that pooling is what supplies the width. The
posterior is right and the sampler is slow, so the statistic a change has
to move is per-chain effective sample size, not coverage
([C1 chain configuration: the coverage deficit is exploration (fef6dca6, 2026-09-07)](release-candidate-review.md#c1-chain-configuration-the-coverage-deficit-is-exploration-fef6dca6-2026-09-07)).
That ruling is what made this a mixing program.

## 3. How things were measured

The battery is a fixed set of test problems, each with a stated failure
mode, a truth to score against and one statistic chosen before the run. It
divides into an average-case core, which gates, and pathologies, each
isolating a mode the others do not
([6. The battery](../design/benchmark-surfaces.md#6-the-battery)). Five
cells are built:

- **C1**, the He and Hahn factorial: ten thousand rows, thirty continuous
  predictors, moderate noise, two mean functions - a trigonometric
  polynomial carrying one true interaction, and a single-index rotated
  ridge. The only core cell with a published BART coverage number, and the
  only cell where a kernel change has moved anything.
- **P1**, a low-noise Friedman emulator: Pratola's setting, where structure
  freezes as the residual variance falls. The known-positive control.
- **P2**, the confounded step function: three columns, three hundred rows,
  one tree, two exactly equiprobable representations of the same fit, so
  the answer is a 0.5 symmetry the sampler finds or does not.
- **P5**, a checkerboard interaction on columns correlated at 0.9, whose
  inclusion truth is exact and whose near-decoys are strong.
- **P6**, a diagonal shelf with targeted selection: 250 rows, an outer
  causal estimand, a published BART bias of 0.27 at 65 percent coverage.

The accept rule is asymmetric. A change is accepted only if it is neutral
or better on every core cell and better on at least one pathology: the core
is a gate, not a score, so improving it earns nothing and regressing one
cell refuses the change however many pathologies it wins
([6.1 The rule, stated operationally](../design/benchmark-surfaces.md#61-the-rule-stated-operationally)).

The no-regression margins are per cell, at twenty matched pairs, a flag
counting only when the paired mean difference is worse than the margin and
its one-sided 95 percent bound excludes the null:

    95% pointwise coverage of true f          -0.010 absolute
    held-out RMSE against true f              ratio above 1.02
    minimum ESS over 25 fixed points, /second ratio below 0.90
    summed inclusion share on true columns    -0.010 absolute
    outer estimand interval coverage          -0.010 absolute
    wall time per sweep                       ratio above 1.05

Winning a pathology is harder than not regressing the core: the improvement
must exceed four times the measured per-replicate standard error. Three
controls sit on top of the paired ones. A kernel added at weight zero must
be bitwise identical to the control, not merely statistically equal; a
statistical pass with a bitwise failure means the change is not what it
claims to be. Any flagged cell takes a mandatory re-run on fresh seeds
before the flag counts. And a sham arm - the control against itself at a
different sampler seed - calibrates the bar: on C1 it reads -2.3 in summed
minimum effective sample size with a standard deviation of 9.7, which is
what the +8 improvement bar was set against. The absolute gate is P1: if
its 90 percent coverage does not sit near 0.71 in the control arm, the
harness is mismeasuring and no verdict is valid. It reads 0.725 held-out at
the shipped default
([10.8 P1, the low-noise Friedman emulator (2026-09-07)](../design/benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)).

One measurement changed how that coverage margin is read. C1's shipped
configuration - four chains of 500 burn-in and 500 kept, pooled -
over-covers, and the reason is now measured rather than inferred. A
well-mixed reference arm, four long chains under the best-mixing kernel
available, reads 0.941 against a nominal 0.95 where the shipped
configuration reads 0.961, and its coverage has converged in length: half
its kept draws read 0.939. Doubling the chain count at one kernel and one
length moves coverage from 0.939 to 0.956, where quintupling the length
moves it from 0.939 to 0.941. Chain count, not chain length, inflates the
pooled interval, so the extra coverage is unfinished mixing, and coverage
on this cell is read against the reference rather than the
shipped-configuration control, at the same -0.010 margin
([6.4 What "no regression on the core" means numerically](../design/benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)).
The same substitution applies to held-out root mean squared error on the
single-index mean function, where the reference reads 1.024 against the
control.

## 4. The shipped kernel today

The engine holds five structural tree moves: birth and death, change,
swap, perturb and the nog-node rule draw named `rule_gibbs`. One uniform
per tree per sweep selects among them, at the default mixture

    birth_death 0.6   swap 0   change 0.4   perturb 0   rule_gibbs 0

with birth taking 0.5 of the birth/death share. So the kernel a user runs
is fringe birth and death plus change, which redraws an interior node's
split variable and cut with the skeleton below held fixed. The other three
are reachable through `proposal.probs` and ship at zero. Setting all five
to zero freezes the forest: no move is proposed, none of its randomness
drawn, while leaf values, sigma and the family's latents keep sampling
([Tree moves](../architecture.md#tree-moves)).

One leaf-value step sits beside the tree moves. `levelGibbs` adds a
constant to every occupied leaf of a tree, the constants summing to zero
across the forest, so the fitted function is unchanged exactly. Its default
`NA` means automatic: the step is taken for a forest exactly where that
forest's structural mixture is frozen; `TRUE` and `FALSE` override
([One sweep](../architecture.md#one-sweep)).

## 5. What was tried

### The swap proposal: removed, then restored at zero

Swap exchanges the rules of a parent and one child. The move census priced
it: 70 to 77 percent of its proposals never reach a score at all, so it
accepts 1.7 to 4.3 percent of what it proposes though 5.7 to 15.1 percent
of what it scores. The weakness is the no-op rate, not the acceptance
ratio. On the one criterion where the shipped mixtures separate - how fast
the trees re-adapt after the response is swapped under them, this package's
distinguishing use - the arm carrying change with swap at zero matches the
shipped default on every contrast, while the arm dropping both moves is the
only one that loses
([1. The decision, and its evidence](../design/swap-removal.md#1-the-decision-and-its-evidence)).
Swap's 0.1 was accordingly moved to birth/death and the move deleted.

It came back at a default of zero, on an exact-posterior gate rather than a
mixing argument. One gate fits a single tree on two live columns and
compares the sampler's subject-level hazard against a brute-force
enumeration of all 62 reachable trees. Without swap it fails, maximum
hazard gap 0.0120 against a tolerance of 0.004; with it, 0.0008, and twenty
of twenty-one gates pass against twenty-one of twenty-one. The mechanism is
specific to one tree: once the splits beneath the root depend on the root's
own variable, no new variable has a valid cut there and change no-ops, and
swap alone rotates a child's rule up. At fifty and two hundred trees the
ensemble averages the trap away, zero stuck trees in either arm at either
count, so no default moved
([9. Reversal: the move returns at default zero](../design/swap-removal.md#9-reversal-the-move-returns-at-default-zero)).
What remains in the code is the move, its name in `proposal.probs`, its
detailed-balance gate, and a positive swap share in every one-tree exact
gate that accepts a caller mixture.

### Perturb, a same-variable cut move

Change always redraws the split variable, so a pure cut displacement
happens only when the redraw lands back on the incumbent. The census made
the case for building one: change accepts 1.7 to 6.1 percent of its
proposals by cell, and its rejections are not close calls, the median
log-likelihood difference among them running -62 to -143. A one-position
displacement of the same variable's cut accepts 27 to 51 percent - the
random-walk step-size argument made concrete
([1. What change does, and why a displacement is a different move](../design/perturb-move.md#1-what-change-does-and-why-a-displacement-is-a-different-move)).

The move landed at weight zero, bitwise neutral, with a detailed-balance
gate carrying two poisoned variants that both fail as designed. Its benefit
stage ran on C1 at a share of 0.16 taken from change and did not move the
primary: +0.1 with a standard deviation of 8.1 in the pilot and -1.9 with
6.5 on the fresh seed block, against a +8 bar and inside the sham arm's own
reading. Held-out error regressed to 1.032 against a 1.02 margin. KILLED at
that setting; the kernel stays at zero, and whether it is removed before
release is an open call
([5.3 What arm B must produce, and the kill](../design/perturb-move.md#53-what-arm-b-must-produce-and-the-kill)).

### The move census and its probes

Before any kernel was built, the engine gained a compile-time census
recording one line per structural proposal and restoring every byte it
touches. No shipped build carries it and it draws nothing, so a census run
at an arm's seeds reproduces that arm's chain exactly. It priced four
mechanisms without writing any of them
([6.1 Stage 0 - the move census (pilot; no kill criterion)](../design/tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion)).

- **An exact rule draw at a nog node**, an interior node whose two children
  are both leaves. Nog nodes are 48.5 to 98.3 percent of interior nodes by
  cell and take 62.7 to 99.1 percent of change's proposals, so they are not
  a corner case. Asked how far the incumbent rule sits from the node's full
  conditional, the probe reads 0.737 for the incumbent at the low-noise
  cell - the least there is to buy - against 0.0017 at C1, a median rank of
  26.5 among up to three thousand candidates. This justified the kernel.
- **Informed death.** Weighting which nog node to prune by its merged-leaf
  marginal has nothing to inform: the weight vector is a point mass and the
  uniform pick already has median rank 1 in every cell.
- **Perturb displacement.** Consecutive accepted displacements at one node
  continue in the same direction 36.6 to 43.1 percent of the time, below
  the 50 percent a reversible walk implies.
- **Lifting**, which carries a direction bit and flips it on rejection, was
  ranked on the premise that same-direction runs exist. They do not, so a
  lift would spend its saving forcing continuation in a direction the chain
  prefers to reverse. Refuted before any code was written.

The same run recorded what shipped trees look like - two to three leaves on
average, not sixteen to thirty-two - and that C1 is not a low-acceptance
regime: pooled scored acceptance there is 24.9 percent against the default
Friedman cell's 8.0. The trees move freely and the deficit survives, which
redirected attention from acceptance to what accepted moves are worth.

### Cross-chain exchange

Exchanging one tree between two chains at the same temperature is ordinary
Metropolis on the product target, the two tree priors cancelling because
the prior never reads a predictor value. It was priced without drawing
anything, from four collapsed marginals in closed form. Acceptance is 10.4
to 11.7 percent on C1 and 4e-8 to 2e-4 on the low-noise cell, and it lives
entirely on the one- and two-split trees birth and death already reach: by
leaf-count pair it runs 0.331 at one leaf against two, down to 4e-22 at
five against five. The accepted exchanges are also the
couplings that would collapse the between-chain spread pooling converts
into this cell's coverage. KILLED
([16.3 Ranking](../design/tree-mixing-proposals.md#163-ranking)).

### The exact rule draw at a nog node, full version

At a nog node the whole conditional over split rules is available for the
price of a scan, and that scan's marginal is exact, the two children being
a two-way partition of the node's members with no skeleton below to reroute
through. The neighbourhood does not depend on the incumbent rule, so the
candidate set is the same from every state in it, the normalizer is
unchanged by the draw, and no proposal count survives into the acceptance:
a Gibbs step at acceptance one, not an informed proposal
([2. The move](../design/nog-gibbs.md#2-the-move)).

Correctness is gated twice. A component test assembles the neighbourhood
once by the kernel and once by an independent reference that installs each
candidate rule and re-scores it from the prior and the branch likelihood;
they agree candidate by candidate. A detailed-balance script then runs the
kernel at a dominant share against three closed-form target statistics,
with two poisons sized in advance - dropping the below-node growth factors,
and dropping the rule prior whose omission was once a real defect in this
package's change move - and both poisons fail while the kernel passes
([5. Correctness: rule-gibbs-balance.R](../design/nog-gibbs.md#5-correctness-rule-gibbs-balancer)).

On C1's four-chain configuration, at a share of 0.16 taken from change,
over twenty matched pairs and then twenty fresh ones:

    statistic (Trig+poly)      control   rule_gibbs 0.16   fresh block
    summed minimum ESS             15         36              +22.1
    per-chain minimum ESS           2          3              +0.86
    95% coverage                0.961      0.939             -0.026
    held-out RMSE ratio          1.00      0.979              0.987

The paired improvements of +21.5 and +22.1 are about eight paired standard
errors, and the per-chain minimum moves with them, which no other arm on
this cell has done. Doubling the share to 0.32 buys more effective sample
size (+37.0) and nothing else. The price was predicted to within one
percent: 2.21 sweeps of arithmetic where a sweep costs one
([6. Benefit, pre-registered](../design/nog-gibbs.md#6-benefit-pre-registered)).

The coverage secondary flagged, at -0.022 and -0.026, and the fresh-seed
re-run confirmed it. That flag has since dissolved: read against the
well-mixed reference rather than the over-covering control, the same arm is
-0.002, inside the margin, because what it narrows is width that unfinished
mixing was supplying. The single-index error regression goes the same way,
1.028 against the control being 1.004 against the reference.

### The same draw restricted to the cut axis

The variable axis is thirty times the price of the cut axis, so a
restricted variant holds the node's incumbent variable and enumerates that
variable's cuts alone - one scan a proposal instead of thirty. It remains
an exact Gibbs step for the same reason the full draw is: the restricted
set is a deterministic function of state the move cannot change. It sits
behind a compile-time switch, off in every shipped build, and the
detailed-balance script passes on it with both poisons still failing.

    statistic (C1, 0.16 share)      full draw        cut-only
    cut scans per sweep               337.32           10.97
    sweep-equivalents of cost           2.21            1.04
    Trig+poly summed min ESS      +21.5, +22.1     +11.3, +9.1
    Single index summed min ESS   +14.4            +16.2
    Trig+poly coverage            -0.022, -0.026   -0.007, -0.012

So the variable axis costs thirty times the scans, buys about half of one
mean function's gain and none of the other's, and buys most of the coverage
movement that flagged the full draw. The restricted kernel clears the +8
bar at both seed blocks on the gated mean function, at a cost
indistinguishable from a plain sweep
([2.4 Cost, and the two restricted variants](../design/nog-gibbs.md#24-cost-and-the-two-restricted-variants)).

[DOSE RESPONSE: pending, filled from nog-gibbs.md section 6 when the run lands]

### The level-fibre step

Adding a constant to every occupied leaf of one tree, the constants summing
to zero across the forest, leaves the fitted function unchanged exactly, so
the conditional on that subspace is the leaf prior alone and the draw is
closed form. It costs about a five-thousandth of a sweep and changes no
rule
([1. The draw](../design/level-fibre.md#1-the-draw)).

Its pilot was spectacular. Freezing C1's structures at a recorded draw and
running leaf and sigma draws alone, the minimum effective sample size over
the cell's twenty-five points rises by a paired median of +189.4 at one
freeze point and +231.6 at the other, all ten pairs positive.

Under the live kernel it is nil: -0.9 and -2.5 at the two seed blocks
against the +8 bar, inside the sham arm's own reading, with coverage flat
to the third digit, held-out error at 1.012 and 1.014, and the
between-chain ratio unmoved. Stacked on the nog-node draw it reads as that
kernel alone in every column. KILLED as a general default. The explanation
is in the pilot's own scope: it measures one channel with the structure
held fixed, and the structural half of the sweep moves the level fibre
faster than the exact draw pays for itself
([6. Benefit, pre-registered](../design/level-fibre.md#6-benefit-pre-registered)).

What was kept is the regime where the pilot's gain is the whole story. The
control slot became tri-state, its default meaning automatic: the step runs
for a forest exactly where that forest's mixture is frozen. The decision is
per sweep and per forest, because the mixture is mutable between samples
while the slot is fixed at sampler creation, so a driver loop that freezes
structure between response swaps would otherwise never reach it. A forest
that is not frozen consumes no generator draw, so the shipped engine is
bitwise unchanged
([4. The surface](../design/level-fibre.md#4-the-surface)).

### The tree count, and grow-from-root

The tree-count contrast is the same mechanism seen from the prior's side
and is not a kernel question. It is parked: four pooled chains at the
shipped length buy more coverage than 200 trees do.

Grow-from-root as a default was killed earlier, on its own harm battery.
The early-iteration benefit is real and every aggregate test passed; the
per-cell check caught a plateau posterior-mean error cost of +11.10 percent
in one noise-heavy cell against a 6 percent margin, confirmed on fresh
seeds, that pooling averaged away. It stays opt-in, and its design law -
per-cell checks, frozen thresholds, mandatory fresh-seed re-runs, a null
control that voids the family - is what the battery inherits
([5. Verdict and consequences](../design/grow-from-root-default.md#5-verdict-and-consequences)).

### The two brainstorm rounds

Two rounds of proposal generation ran: the first anchored in the literature
across four lenses (tree and partition geometry, latent constructions,
informed proposals, rotations of the input space), the second deriving from
this engine's own structure with the literature held shut. Every candidate
had to state a proposal precisely enough to write its Metropolis-Hastings
correction, price it against one cut scan, name a measured deficit and
supply a one-day falsifier
([16.4 Merged view across the two rounds](../design/tree-mixing-proposals.md#164-merged-view-across-the-two-rounds)).

What survived is what the sections above report. What was refuted is worth
as much: tree-space geodesics fail because regression trees carry no fixed
leaf-label set and the collapsed marginal is piecewise constant; a learned
rotation misreads a chain-exploration deficit, a reparameterization leaving
an equally multimodal posterior not moving a chain that sits in one place;
a per-tree temperature has no Gibbs conditional; and a per-leaf
per-variable histogram cache wins a depth-fold, about two at the measured
tree size, rather than the leaf-fold claimed
([15.5 Refuted along the way](../design/tree-mixing-proposals.md#155-refuted-along-the-way)).

Two constructions remain unbuilt with their arguments intact. A
pairwise-collapsed split transfer, moving a split from one tree to another
without either passing through a stump, is the only candidate that
addresses representation multimodality directly; its probe was not built
because pricing a transfer needs the residual net of the other trees and
the partner's leaf statistics, which no existing move sees. A lifted
birth/death is cheap and valid but recomputes at a gain of 1.06 to 1.17, a
rider on an aim improvement rather than anything alone
([16.5 Discarded across both reports](../design/tree-mixing-proposals.md#165-discarded-across-both-reports)).

## 6. What it adds up to

**The deficit has two halves and they sit at different points.** Freezing
the structures on C1 and running the leaf draws alone, the median point's
effective sample size rises from 14.9 to about 670 of 2500 kept and its
lag-one autocorrelation falls from 0.70 to about 0.35: at a typical
coordinate the deficit is overwhelmingly structural. At the worst
coordinate it is not. The minimum rises only from 1.6 to between 4 and 21,
and that coordinate carries two to three times the median posterior spread.
The leaf draw is itself slow there, which is why the level-fibre step helps
in a frozen chain and why a structural kernel alone cannot be the whole
answer.

**Chain count is doing work that looks like coverage.** Four short chains
each sit in their own place - the between-chain ratio runs 0.5 to 0.8 - and
pooling them widens the interval to 0.961 on a nominal 0.95, where the
well-mixed reference reads 0.941. Any arm that makes the chains agree
narrows that interval and looks like a coverage regression until it is read
against the reference. That is what happened to the nog-node draw, and why
the core's coverage margin now names the reference.

**The only kernel that has moved either half is the exact rule draw**,
which raises the summed minimum from 15 to 36 at the low dose and 52 at the
high one and moves the per-chain minimum, which every other arm left at 2.

**The candidate for a default is the cut-only rule draw.** It clears the
improvement bar on the gated mean function at both seed blocks, gains the
whole of the full draw's benefit on the ungated one, moves coverage only to
the edge of the margin, and costs 1.04 sweep-equivalents rather than 2.21,
close enough to free that the cost-adjusted question stops mattering. The
dose is not yet named.

[DOSE RESPONSE: pending, filled from nog-gibbs.md section 6 when the run lands]

**What a default change costs.** Any nonzero share is a stream shift: every
draw moves from the first sweep of the first chain, so the three
equivalence baselines, the documentation counts computed from them, the
gate ledger's baseline footnote and four seeded-drift regression files
re-record in one bundle, and the consumer packages on their lockstep
branches re-record with them. That is a day of mechanical work, not a risk;
the bitwise oracle is the usual one.

The decision is whether to spend that re-record and about one percent of
sweep cost to roughly double the worst-coordinate effective sample size on
the average-case cell, on one cell at two seed blocks, with the pathologies
untouched and the equal-cost arm and plateau-error gate owed.

## 7. What was not measured

- **Response-swap recovery inside the battery.** How fast the trees
  re-adapt after the response moves under them was measured once, on its
  own grid, and decided swap's default. It has not been re-run against any
  new kernel, and the cell that would carry it is unbuilt
  ([14. Recovery after a response swap (2026-09-06)](../design/tree-mixing-proposals.md#14-recovery-after-a-response-swap-2026-09-06)).
- **The equal-cost arm**, in which the control is given the extra sweeps
  the new kernel's cost ratio buys it, on a quiet machine. It answers the
  per-second question and is owed. Wall time was not readable on the host
  these runs used, so no cost claim here is anything but a cut-scan
  count.
- **The plateau-error gate.** The house rule requires a per-cell harm check
  on plateau prediction error in a noise-heavy or large-n stratum before
  any default flips. No built cell is either, and none measures it
  ([6.4 Kill criteria, pre-registered](../design/tree-mixing-proposals.md#64-kill-criteria-pre-registered)).
- **Whether composing BART with a parametric block helps tree-space mixing
  at all.** Nobody has measured it, in print or here; the hazard of that
  composition is measured at a factor of six and the benefit never
  ([9. What this survey could not settle](../design/tree-mixing-proposals.md#9-what-this-survey-could-not-settle)).
- **The embedded cells.** The core cell that runs a sampler inside an outer
  loop and the funnel pathology both need an outer sampler living in
  another repository.
- **Seven cells remain unbuilt**: the real-data convergence ladder, the
  inhomogeneous-smoothness cell, the no-overlap extrapolation cell, the
  variance funnel, the mixed-type regime switch, the real-covariate causal
  cell and the embedded moving-response cell
  ([10.6 What the pilot could not do](../design/benchmark-surfaces.md#106-what-the-pilot-could-not-do)).
- **Non-gaussian families.** Every built cell is gaussian or causal
  gaussian, and coverage of a true mean is a calibration statistic, not a
  proof the sampler targets the right posterior.

## 8. Where the evidence lives

- `docs/design/benchmark-surfaces.md` - the battery: its cells, accept rule
  and margins, and every run of the five built cells.
- `docs/design/tree-mixing-proposals.md` - how the posterior is sticky, the
  ranked candidates, the move census and its probes, the brainstorm rounds.
- `docs/design/nog-gibbs.md` - the exact rule draw at a nog node: cost,
  correctness gate, benefit verdict, cut-only pilot.
- `docs/design/level-fibre.md` - the level-fibre step: derivation, pilot,
  live-kernel kill, automatic mode.
- `docs/design/perturb-move.md` - the same-variable cut move and its kill.
- `docs/design/swap-removal.md` - swap's removal and the gate that brought
  it back at zero.
- `docs/design/grow-from-root-default.md` - the earlier warm-start study,
  whose per-cell design law the battery inherits.
- `docs/plans/release-candidate-review.md` - landing notes, in date order.
- `docs/architecture.md` - what one sweep does and what each move is.
