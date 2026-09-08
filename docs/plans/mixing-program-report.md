# The mixing program: what was measured, what was found, what to do

Status: REPORT, 2026-09-08. A record; no default is changed here.

## 1. Summary

BART's tree sampler moves by small local edits, and on hard problems it
settles into one arrangement of trees and stays there. The visible symptom
is under-coverage of the fitted function; the underlying quantity is
effective sample size, which on the one core problem carrying a published
BART coverage number sits at two draws out of 2500.

The package now has an evaluation battery with a pre-registered accept
rule. Five candidate remedies were assessed: three through its cells, one
on evidence predating it, one priced in closed form. Four do not earn a
default. Swap earns nothing at production forest sizes and ships at zero. A
same-variable cut move was killed on its registered statistic, leaving an
unclaimed gain on the ungated mean function. An exchange of trees between
chains was killed on acceptance. An exact draw on the forest's level fibre
is nil under the live kernel and is kept only where structure is frozen,
now its default.

One succeeded on its pre-registered statistic: replacing the Metropolis
change proposal at a node whose two children are both leaves with an exact
draw from that rule's full conditional. It raises effective sample size on
the core cell by 1.8 to 3.0 times the pre-registered bar, at both seed
blocks and on both mean functions, and a restricted version redrawing only
the cut keeps about half that gain on the gated function and all of it on
the other, for a thirtieth of the arithmetic. The proposal is to
consider that restricted draw for a nonzero default share at `d` = 0.16,
the dose the dose-response run names, and to keep every other kernel at
zero. Section 6 states what the proposal lacks: the accept rule asks for a
pathology win and this kernel has none, three of the battery's four core
cells do not exist, and the cost answer rests on a cut-scan count.

## 2. The question

BART's intervals do not always cover, and colleagues report it. The
literature's cleanest number for the average case is He and Hahn's
factorial study: 95 percent pointwise coverage of the true mean function of
0.73 to 0.74 at ten thousand rows, thirty predictors and moderate noise.
That setting is this battery's core cell, C1, the only one carrying a
published BART coverage number. dbarts reproduces the deficit in direction
and about half in size - 0.82 against a nominal 0.95
at the setting closest to the paper's, at intervals 15 to 25 percent longer
for the same point accuracy, so the package is better calibrated than the
published BART and still short of nominal
([10.4 C1, the He and Hahn factorial](../design/benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial)).

Two levers could close that gap: more trees, or better mixing. The
tree-count lever works - 200 trees against 75 buy about ten points of
coverage on both mean functions at no cost in error - and was parked
anyway, on a later measurement of this cell's coverage arithmetic. Three
chain arms on the same twenty seeds read

    four chains, 500 burn-in and 500 kept, pooled   0.961
    four chains, 1000 burn-in and 2500 kept         0.959
    one chain, 1000 burn-in and 25000 kept          0.902

the long chain costing about six times the first arm's sweeps, and the four
disagreeing enough that pooling supplies the width. The posterior is right
and the sampler is slow. The maintainer's ruling followed: the coverage
deficit is a mixing symptom and mixing is the lever, so the statistic a
change has to move here is per-chain effective sample size, not coverage
([5.1 The chain configuration, and what it makes the primary statistic](../design/perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic)).
Section 6 returns to what that commits the battery to.

## 3. How things were measured

The battery is a fixed set of test problems, each with a failure mode, a
truth and one statistic chosen before the run: a core of four cells, which
gates, and eight pathologies, each isolating its own mode
([6. The battery](../design/benchmark-surfaces.md#6-the-battery)). Five are
built.

- **C1**, the He and Hahn factorial, on two mean functions: a trigonometric
  polynomial with one true interaction, written Trig+poly and the only one
  the accept rule gates, and a single-index rotated ridge, Single index,
  reported beside it.
- **P1**, a low-noise Friedman emulator, where structure freezes as the
  residual variance falls: the known-positive control. The rung in force is
  the house rung at n = 2000 and sigma 0.25, not Pratola's n = 5000 cell
  and not its published 53 percent coverage.
- **P2**, the confounded step function: one tree, two exactly equiprobable
  representations of one fit, so a 0.5 symmetry the sampler finds or does
  not.
- **P5**, a checkerboard interaction on correlated columns, with an exact
  inclusion truth and near-decoys.
- **P6**, a diagonal shelf with targeted selection: an outer causal
  estimand at a published BART bias of 0.27 and 65 percent coverage.

The accept rule is asymmetric. A change is accepted only if it is neutral
or better on every core cell and better on at least one pathology: the core
is a gate, not a score, so improving it earns nothing and regressing one
cell refuses the change however many pathologies it wins
([6.1 The rule, stated operationally](../design/benchmark-surfaces.md#61-the-rule-stated-operationally)).

The no-regression margins are per cell, at twenty matched pairs,
Holm-corrected across cells within a metric, a flag counting only when the
paired mean difference is worse than the margin and its one-sided 95
percent bound excludes it:

    95% pointwise coverage of true f          -0.010 absolute
    held-out RMSE against true f              ratio above 1.02
    minimum ESS over 25 fixed points, /second ratio below 0.90
    summed inclusion share on true columns    -0.010 absolute
    outer estimand RMSE                       ratio above 1.02
    outer estimand interval coverage          -0.010 absolute
    ESS of the outer scale                    ratio below 0.90
    wall time per sweep                       ratio above 1.05

Winning a pathology is harder: the improvement must exceed four times the
measured per-replicate standard error on that cell's own statistic, and any
flagged cell takes a mandatory fresh-seed re-run first. Two absolute gates
sit on top. A kernel added at weight zero must be bitwise identical to the
control; and P1 checks the harness, since unless its 90 percent coverage
sits near 0.71 in the control arm no verdict is valid. It reads 0.725
held-out at the shipped default
([10.8 P1, the low-noise Friedman emulator (2026-09-07)](../design/benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07)).

The improvement bar on C1 is +8 in summed minimum effective sample size:
four times a paired standard error of 2.0, from that cell's own recorded
spread and matched by the 1.8 to 2.5 its eight move-set cells measure. A
sham arm - the control against itself at a different sampler seed - ran
once as a check on that bar, not its source, and reads -2.3 +/- 9.7 at a
paired standard error of 2.17, four times which is 8.7, so +8 is marginally
optimistic
([5. Benefit, pre-registered](../design/perturb-move.md#5-benefit-pre-registered)).
The sham has not been re-run since; every arm below rests on that reading.

**The baseline C1's coverage is read against.** C1's shipped four chains of
500 burn-in and 500 kept, pooled, over-cover, so coverage here is read
against a well-mixed reference arm rather than that control, at the -0.010
margin
([6.4 What "no regression on the core" means numerically](../design/benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically)).
The reference has to be named, being not independent of what it judges: it
is section 5's candidate kernel in full, at 0.32 - twice the share of the
arm whose coverage flag it dissolves - at four chains of 1000 burn-in and
2500 kept. The eight-chain arm separating count from length is the same
kernel again.

    arm                                    95% coverage, Trig+poly
    four chains 500 + 500, shipped                   0.961
    reference, four chains 1000 + 2500               0.941
    that reference fit at half its length            0.939
    the same kernel at eight chains of 500 + 500     0.956

Chain count, not chain length, is what inflates the pooled interval.

The reference is credible on three grounds: its between-chain ratio of 0.48
is the lowest recorded here against the control's 0.78; its lengths agree;
and the kernel producing it is an exact Gibbs step on the posterior the
shipped kernel targets, so a bias not running through mixing would have to
be a defect both its gates missed - a candidate-by-candidate identity test
against an independent reference assembly, and a prior-only
detailed-balance gate with two poisons that both fail. What it is not is
independent: no kernel but the candidate mixes well enough here to serve as
a reference, which is the finding itself. Residual disagreement widens a
pooled interval, so 0.941 is an upper bound. The same substitution applies
to Single index held-out error, at 1.024.

## 4. The shipped kernel today

The engine holds five structural tree moves: birth and death, change, swap,
perturb and the nog-node rule draw `rule_gibbs`, a nog node being an
interior node whose two children are both leaves. One uniform per tree per
sweep selects among them, at the default mixture

    birth_death 0.6   swap 0   change 0.4   perturb 0   rule_gibbs 0

with birth taking 0.5 of the birth/death share
([9. Reversal: the move returns at default zero](../design/swap-removal.md#9-reversal-the-move-returns-at-default-zero)).
So the kernel a user runs is fringe birth and death plus change, which
redraws an interior node's split variable and cut with the skeleton below
held fixed. The other three ship at zero, reachable through
`proposal.probs`, the R argument carrying the mixture; setting all five to
zero freezes the forest, leaving leaf values, sigma and the latents
sampling ([Tree moves](../architecture.md#tree-moves)).

One leaf-value step sits beside them: `levelGibbs` adds a constant to every
occupied leaf of a tree, the constants summing to zero across the forest,
leaving the fitted function unchanged exactly. Its default `NA` means
automatic - taken for a forest exactly where that forest's mixture is
frozen - and `TRUE` and `FALSE` override
([One sweep](../architecture.md#one-sweep)).

## 5. What was tried

### The swap proposal: removed, then restored at zero

Swap exchanges the rules of a parent and one child, and its evidence is not
a battery cell. The census priced it: 70 to 77 percent of its proposals
never reach a score, so it accepts 1.7 to 4.3 percent of what it proposes
though 5.7 to 15.1 percent of what it scores. On the one criterion where
the shipped mixtures separate - re-adaptation after the response is swapped
under the trees - change with swap at zero matches the default on every
contrast, and only the arm that drops both loses
([1. The decision, and its evidence](../design/swap-removal.md#1-the-decision-and-its-evidence)).
Swap's 0.1 moved to birth/death and the move was deleted, then returned at
a default of zero on an exact-posterior gate: one tree, two live columns,
the subject-level hazard against a brute-force enumeration of all 62
reachable trees, where the maximum gap is 0.0120 against a 0.004 tolerance
without swap and 0.0008 with it. Swap alone rotates a child's rule up; at
fifty and two hundred trees the ensemble averages the trap away, so no
default moved
([9. Reversal: the move returns at default zero](../design/swap-removal.md#9-reversal-the-move-returns-at-default-zero)).
Every one-tree exact gate that takes a caller mixture now sets a positive
swap share; five cannot - four BCF gates whose two-forest path refuses a
non-default mixture, and one monotone gate that rewrites it - and those
five are the standing exposure.

### Perturb, a same-variable cut move

Change always redraws the split variable, so a pure cut displacement
happens only when the redraw lands back on the incumbent. The census made
the case for building one: change accepts 1.7 to 6.1 percent of its
proposals by cell against a one-position displacement's 27 to 51. It landed
at weight zero, bitwise neutral, behind a detailed-balance gate whose two
poisoned variants both fail as designed, and its benefit stage ran on C1 at
a share of 0.16 taken from change without moving the gated statistic:

    Trig+poly summed minimum ESS   +0.1 +/- 8.1 (pilot)   -1.9 +/- 6.5 (fresh)

against the +8 bar and inside the sham arm's own reading, with held-out
error at 1.032 against the 1.02 margin. KILLED at that setting. On the
ungated function, recorded as residue rather than survival, it gains +9.5
summed on the fresh block, 19 of 20 seeds at t 5.35, with its own 1.027
held-out error regression
([5. Benefit, pre-registered](../design/perturb-move.md#5-benefit-pre-registered)).
The kernel stays at zero, its removal a decision the maintainer takes after
this report.

### The move census and its probes

A compile-time census logs one line per structural proposal and restores
every byte it touches; it draws nothing, so a census run at an arm's seeds
reproduces that chain exactly
([6.1 Stage 0 - the move census (pilot; no kill criterion)](../design/tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion)).
Nog nodes are 48.5 to 98.3 percent of interior nodes and take 62.7 to 99.1
percent of change's proposals, and the incumbent rule sits far from the
node's full conditional where it matters:

    P(incumbent), the low-noise cell   0.737, the least there is to buy
    P(incumbent), C1                   0.0017, median rank 26.5 of up to 3000

That justified the kernel, and killed informed death on a point-mass weight
vector and a lifted cut displacement on the absence of same-direction runs.
It also put C1's pooled scored acceptance at 24.9 percent against the
default Friedman cell's 8.0, so its minimum effective sample size of 2 is
no acceptance deficit - the reading that elevated the frozen-structure
test.

### Cross-chain exchange

An exchange of one tree between two chains at one temperature is ordinary
Metropolis on the product target, the two tree priors cancelling, and was
priced in closed form without drawing anything.

    acceptance at C1's two seeds                     10.1 to 11.7 percent
    acceptance at the low-noise cell                 4e-8 to 2e-4
    by leaf-count pair, one against two              0.331
    by leaf-count pair, five-plus against five-plus  4e-22

It lives on the one- and two-split trees birth and death already reach and
vanishes on the multi-split trees it was sized for; the accepted exchanges
are the couplings that would collapse the between-chain spread pooling
turns into coverage. KILLED, on acceptance
([16.3 Ranking](../design/tree-mixing-proposals.md#163-ranking)).

### The exact rule draw at a nog node, full version

At a nog node the whole conditional over split rules costs one scan, and
that scan's marginal is exact, the two children being a two-way partition
of the node's members with no skeleton below. The neighbourhood does not
depend on the incumbent rule, so the candidate set is the same from every
state in it and no proposal count survives into the acceptance. The law has
two cases and they differ: where the incumbent sits in the candidates' own
empty-leaf-veto rank stratum - the only case with no weight mask and no
missing value at the node - the step is an exact Gibbs draw at acceptance
one; where the stratum entered is strictly better it is
Metropolis-within-Gibbs, also at acceptance one, valid because that stratum
is absorbing, with no stationarity claimed for the one it leaves
([2. The move](../design/nog-gibbs.md#2-the-move)). Its two correctness
gates are the ones section 3 names
([5. Correctness: rule-gibbs-balance.R](../design/nog-gibbs.md#5-correctness-rule-gibbs-balancer)).

The benefit stage ran on C1's shipped four-chain configuration at a share
of 0.16 from change, twenty matched pairs then twenty fresh. Absolute
readings on the gated function, first block:

    Trig+poly, seeds 1-20      control    rule_gibbs 0.16
    summed minimum ESS            15            36
    per-chain minimum ESS          2             3
    95% coverage               0.961         0.939
    interval length             4.61          3.98
    between-chain ratio         0.78          0.58

Paired differences on the same seeds, first block then the fresh one:

    summed minimum ESS      +21.5 +/- 12.8, t 7.5   +22.1 +/- 12.2, t 8.1
    per-chain minimum ESS   +1.14, t 8.7            +0.86, t 7.2
    95% coverage            -0.022, t -12.5         -0.026, t -9.9
    held-out RMSE ratio       0.979                   0.987

Those improvements are 2.7 and 2.8 times the +8 bar, and the per-chain
minimum moves with them, which no other arm here has done. Doubling the
share to 0.32 buys +37.0 and nothing else. Cost is 2.21 sweeps where a
sweep costs one. On the ungated function the arm gains +14.4 and +23.9
summed minimum ESS and pays 1.028 and 1.023 in held-out error, past the
1.02 margin at both blocks
([6. Benefit, pre-registered](../design/nog-gibbs.md#6-benefit-pre-registered)).

The coverage secondary flagged and the fresh-seed re-run confirmed it. That
flag has since dissolved: read against section 3's reference rather than
the over-covering control, the arm is -0.002, inside the margin, because
what it narrows is width unfinished mixing was supplying. The Single index
error regression goes the same way, 1.028 against the control being 1.004
against the reference.

### The same draw restricted to the cut axis

The variable axis is thirty times the price of the cut axis, so a
restricted variant holds the incumbent variable and enumerates its cuts
alone - one scan a proposal instead of thirty. It is still an exact Gibbs
step, its candidate set a deterministic function of state the move cannot
change, and sits behind a compile-time switch, off in every shipped build;
the balance script passes with both poisons failing.

    C1, 0.16 share                  full draw        cut-only
    cut scans per sweep               337.32           10.97
    sweep-equivalents of cost           2.21            1.04
    Trig+poly summed min ESS      +21.5, +22.1     +11.3, +9.1
    Single index summed min ESS   +14.4            +16.2
    Trig+poly coverage            -0.022, -0.026   -0.007, -0.012
    Single index held-out RMSE       1.028            1.033

So the variable axis costs thirty times the scans, buys about half the
gated function's gain and none of the other's, and accounts for most of the
coverage movement that flagged the full draw
([2.4 Cost, and the two restricted variants](../design/nog-gibbs.md#24-cost-and-the-two-restricted-variants)).

What it has not been tested against: its runs carried the bitwise null and
the balance script and nothing else - no P1 control reading, no P2 arm, no
sham, no fresh-seed block on the ungated function, all of which the full
draw's run carried. The cut-scan figures come from a census cell at its own
data and seed on one chain of 200 + 500 sweeps, not a replay of the arms'
twenty, and that cell takes its 0.16 share from change and birth/death
both, neither departure scaled for.

**Dose response.** A third dose below the pilot's and both of its own, over
the same twenty matched pairs, each taking a fresh block:

    paired against the control        d 0.08   d 0.16   d 0.32
    Trig+poly summed min ESS, 1-20      +8.1    +11.3    +13.4
    Trig+poly summed min ESS, 21-40     +1.8     +9.1     +9.3
    Single index summed min ESS, 1-20   +8.3    +16.2    +27.5
    Trig+poly held-out RMSE, 1-20      0.999    0.992    1.057
    Trig+poly held-out RMSE, 21-40     0.994    1.012    1.048
    sweep-equivalents of cost           1.02     1.04     1.08

Only the 0.16 scan count is measured; the other two scale it linearly in
the share, holding the tree population fixed. `d` = 0.08 does not survive
its fresh block, t 1.29 against the +8 bar; 0.16 and 0.32 confirm. Coverage
flags nowhere - against section 3's reference every Trig+poly dose sits
0.011 to 0.015 above 0.941 - but 0.32 regresses held-out error past 1.02 at
both blocks, and the restriction is the reason. The cut-only share comes out of
change, the only move that changes a node's variable, so at 0.32 change is
left at 0.08, the pooled interval goes back up to 4.58 and 4.62 from 4.33
and 4.36, and the chains agree better on a worse fit. The full draw at that
dose does the opposite, +37.0 at a held-out 0.980, carrying the variable
axis itself; Single index, where that axis has little to buy, stays
monotone at 1.03
([6. Benefit, pre-registered](../design/nog-gibbs.md#6-benefit-pre-registered)).

### The level-fibre step

Because the constants of section 4's step leave the fitted function
unchanged, the conditional on that subspace is the leaf prior alone and the
draw is closed form, at about a five-thousandth of a sweep
([2. Placement and cost](../design/level-fibre.md#2-placement-and-cost)).
Its pilot was spectacular: with C1's structures frozen and only leaf and
sigma draws running, the minimum effective sample size over the cell's
twenty-five points rises by a paired median of +189.4 at one freeze point
and +231.6 at the other, all ten pairs positive. Under the live kernel it
is nil.

    paired difference             seeds 1-20        fresh seeds
    Trig+poly summed min ESS      -0.9              -2.5
    Single index summed min ESS   -5.3, t -2.90     +1.7
    Trig+poly held-out RMSE        1.012             1.014
    Trig+poly between-chain        0.79 vs 0.78      0.79 vs 0.80

against the +8 bar and inside the sham arm's own reading, with coverage
flat to the third digit. Two readings ran against the step and neither
survived: the Single index loss at the first block is separated from noise,
and its wall ratio of 1.086 triggered the fresh-seed re-run, which
reproduced neither. Stacked on the nog-node draw the step reads as that
kernel alone in every column. KILLED as a general default, the explanation
being the pilot's scope: it measures one channel with the structure fixed,
where the sweep's structural half moves the level fibre faster than the
exact draw pays for itself
([6. Benefit, pre-registered](../design/level-fibre.md#6-benefit-pre-registered)).

What was kept is the regime where the pilot's gain is the whole story: the
control slot became tri-state, its default automatic, so the step runs for
a forest exactly where that forest's mixture is frozen. That is decided per
sweep and per forest, the mixture being mutable between samples while the
slot is fixed at creation, so a driver loop freezing structure between
response swaps would never reach it. A non-frozen forest draws nothing, so
the shipped engine is bitwise unchanged
([4. The surface](../design/level-fibre.md#4-the-surface)).

### The tree count, and grow-from-root

The tree-count contrast is the same mechanism from the prior's side, not a
kernel question, and it is parked. Against chain count it is not one-sided:
on Trig+poly four pooled chains at 75 trees read 0.961 where one chain at
200 trees reads 0.922, while on Single index the ordering reverses, 0.895
against 0.924.

Grow-from-root as a default was killed earlier in both strata: every
aggregate test passed, but the per-cell check caught plateau posterior-mean
error costs of +11.10 percent in a noise-heavy small-n cell and of +4.47
and +10.66 in a large-n one, each past its frozen margin, confirmed on
fresh seeds and averaged away by pooling. Its design law - per-cell checks,
frozen thresholds, mandatory fresh-seed re-runs, a null control voiding the
family - is what the battery inherits
([5. Verdict and consequences](../design/grow-from-root-default.md#5-verdict-and-consequences)).

### The two brainstorm rounds

Two rounds of proposal generation ran, each candidate obliged to write its
own Metropolis-Hastings correction, price itself against one cut scan, name
a measured deficit and supply a one-day falsifier. Tree-space geodesics, a
learned rotation and a per-tree temperature each fail on their own terms
([15.5 Refuted along the way](../design/tree-mixing-proposals.md#155-refuted-along-the-way)),
while a per-leaf per-variable histogram cache was refuted on size in the
second round's ranking instead - a depth-fold, about two at the measured
tree size, not the leaf-fold claimed. Two constructions of that ranking
remain unbuilt with their arguments intact
([16.3 Ranking](../design/tree-mixing-proposals.md#163-ranking)). A
pairwise-collapsed split transfer, moving a split from one tree to another
without either passing through a stump, is the only candidate addressing
representation multimodality directly; its probe was not built because
pricing a transfer needs the residual net of the other trees and the
partner's leaf statistics, which none of the moves this round's probes rode
inside can see. A lifted birth/death is cheap and valid but recomputes at a
gain of 1.06 to 1.17, a rider rather than anything alone.

## 6. What it adds up to

**The deficit has two halves and they sit at different points.** With C1's
structures frozen and the leaf draws running alone, the median point's
effective sample size rises from 14.9 to about 670 of 2500 kept and its
lag-one autocorrelation falls from 0.70 to 0.34-0.40: at a typical
coordinate the deficit is structural. At the worst it is not - the minimum
rises only from 1.6 to between 4 and 21, at a coordinate carrying two to
three times the median posterior spread. The leaf draw is slow there, so
the level-fibre step helps a frozen chain and a structural kernel cannot be
the whole answer.

**Chain count is doing work that looks like coverage.** Four short chains
each sit in their own place, and pooling them widens the interval to 0.961
on a nominal 0.95 where the reference reads 0.941. Any arm that makes the
chains agree narrows that interval and looks like a coverage regression
until read against the reference - which is what happened to the nog-node
draw, and why the core's margin now names it.

**The exact rule draw is the only kernel measured here that has moved
either half.** It raises the summed minimum from 15 to 36 at the low dose
and 52 at the high one and moves the per-chain minimum, which every other
arm here left at 2. Elsewhere a kernel change has moved things too, though
not this one: on P5, dropping change costs 0.027 of the true columns'
inclusion share and 3.9 percent of held-out error on fresh seeds; on P2 at
one tree it stops the root variable moving.

**Under the accept rule as written, this kernel has earned nothing.** The
rule wants a pathology win and treats the core as a gate; every gain
`rule_gibbs` has recorded is on C1, a core cell. P2 ran as a
must-not-degrade control, and by the design's own account the move cannot
reach the recorded parked chains: at one tree the root is a nog node only
while the tree has one split, and every parked chain sits at three or four
interior nodes with a child split on the confounded column. P1's rung was
re-run as the absolute gate at the shipped mixtures and carries no
`rule_gibbs` arm; nor do P5 and P6. And the core gate is checkable on one
cell of four.

**The case for the kernel rests on a departure the design pre-registered.**
Its kill criterion names five departures from the program's own kill
criteria, two of which matter here: the statistic is minimum effective
sample size rather than coverage, which has no headroom at 0.961, and the
cell is C1 rather than the low-noise cell, on which three shipped mixtures
are indistinguishable. That is section 2's ruling in arithmetic, accepted
when the program was framed around mixing rather than tree count.
Defaulting this kernel accepts that departure: it treats the core cell's
mixing statistic as the target rather than as a gate, and leaves the
pathologies untouched.

**The equal-cost arm is now run, and the kernel wins on it.** The cut-only
draw's 1.04 sweep-equivalents is four percent more arithmetic, not one, and
every cost figure here is a cut-scan count, not a wall time. The shipped
kernel at four chains of 1105 burn-in and 1105 kept is 2.21 times the
control's sweeps, the full draw's own budget:

    at 2210 cut-scan units      summed minimum ESS
    the full draw                      36.3
    2.21 times the sweeps              22.9

so length buys 63 percent of what the kernel buys. Its coverage falls from
0.961 to 0.952 on length alone, a third of the way to the reference, and
its between-chain ratio only to 0.73 from 0.78, the kernel's being 0.58.
Per thousand cut-scan units the cut-only draw at 0.16 reads 25.1 against
the control's 14.8, the full draw's 16.4 and this arm's 10.4. Per second it
reads the other way, on a host at load 100 to 122, so carries nothing
([6. Benefit, pre-registered](../design/nog-gibbs.md#6-benefit-pre-registered)).

**The records now recommend the cut-only rule draw at `d` = 0.16, as a
recommendation and not a decision.** It is the only cut-only dose clearing
the +8 bar at both seed blocks, with every secondary clean at the
reference-read margins - coverage 0.955 and 0.952 against the reference's
0.941, held-out error 0.992 and 1.012 - and it takes the whole of the full
draw's gain on the ungated function. **Three things stand between that
candidate and a default**, none a further measurement of benefit: the
decision, which is the maintainer's and carries the surface question with
it, the restricted kernel being a private compile-time build rather than a
mode anything can select; the plateau-error gate, which does not exist; and
the re-record. The item is scheduled post-release
([8. Slices](../design/nog-gibbs.md#8-slices)).

**What a default change costs.** Any nonzero share is a stream shift
([7. RNG and baselines](../design/nog-gibbs.md#7-rng-and-baselines)), so
the three equivalence baselines, the counts computed from them, the gate
ledger's footnote and four seeded-drift regression files re-record in one
bundle, and the consumer packages on their lockstep branches with them:
mechanical work against a bitwise oracle, not a risk. Defaulting the
restricted draw means either deleting the variable axis from `rule_gibbs`
or adding a sixth name to `proposal.probs` with its own fill rule, formals,
slot, validity term, Rd pages and tests - twenty-four files for the last
name added.

The decision is whether to spend that re-record and four percent of sweep
cost to raise the worst-coordinate effective sample size on the
average-case cell from 15 to about 26, a factor of 1.75, on one cell at two
seed blocks, with the pathologies untouched and the plateau-error gate
owed.

**A standing rule.** No exploration kernel is removed before this report is
read; a kill parks it at zero.

## 7. What was not measured

- **Wall time, and anything per second.** No host was quiet - load ran 4 to
  269 across the arms and 100 to 122 on the equal-cost machine - so every
  cost figure here is a cut-scan count
  ([6. Benefit, pre-registered](../design/nog-gibbs.md#6-benefit-pre-registered)).
- **The plateau-error gate.** The house rule wants a per-cell harm check on
  plateau prediction error in a noise-heavy or large-n stratum before a
  default flips; no built cell is either
  ([6.4 Kill criteria, pre-registered](../design/tree-mixing-proposals.md#64-kill-criteria-pre-registered)).
- **The sampler's own speed since the default moved.** A bench-sampler
  comparison on a quiet machine has been owed since swap's 0.1 went to
  birth/death
  ([8. Landing](../design/swap-removal.md#8-landing)).
- **Response-swap recovery inside the battery.** Measured once on its own
  grid, where it decided swap's default; not re-run against any new kernel,
  its cell unbuilt
  ([14. Recovery after a response swap (2026-09-06)](../design/tree-mixing-proposals.md#14-recovery-after-a-response-swap-2026-09-06)).
- **Whether composing BART with a parametric block helps tree-space
  mixing.** Nobody has measured it, in print or here; the hazard is a
  factor of six, the benefit unmeasured
  ([9. What this survey could not settle](../design/tree-mixing-proposals.md#9-what-this-survey-could-not-settle)).
- **Seven cells remain unbuilt**: the real-data convergence ladder, the
  inhomogeneous-smoothness cell, the no-overlap extrapolation cell, the
  variance funnel, the mixed-type regime switch, the real-covariate causal
  cell and the embedded moving-response cell. The record's list of eight
  includes P1, since built
  ([10.6 What the pilot could not do](../design/benchmark-surfaces.md#106-what-the-pilot-could-not-do)).
  The funnel and the embedded cell need an outer sampler from another
  repository
  ([6.2 The average-case core](../design/benchmark-surfaces.md#62-the-average-case-core)).
- **Non-gaussian families.** Every built cell is gaussian or causal
  gaussian, and coverage of a true mean calibrates rather than proves
  ([6.6 What this battery does not measure](../design/benchmark-surfaces.md#66-what-this-battery-does-not-measure)).

## 8. Where the evidence lives

- `docs/design/benchmark-surfaces.md` - the battery: cells, rule, margins,
  every run of the built cells.
- `docs/design/tree-mixing-proposals.md` - the stickiness survey, ranked
  candidates, move census, brainstorms.
- `docs/design/nog-gibbs.md` - the exact rule draw and its cut-only pilot.
- `docs/design/level-fibre.md` - the level-fibre step, automatic mode.
- `docs/design/perturb-move.md` - the same-variable cut move and its kill.
- `docs/design/swap-removal.md` - swap's removal and restoration.
- `docs/design/grow-from-root-default.md` - the warm-start study whose
  design law the battery inherits.
- `docs/plans/release-candidate-review.md` - landing notes by date.
- `docs/architecture.md` - one sweep, and each move.
