# The tree-sampler mixing program

Status: RECORD, 2026-09-08. The cut-only exact rule draw is adopted, to land
after the first release; no pre-release default changes.

## 1. Summary

BART's tree sampler moves by small local edits, and on hard problems it settles
into one arrangement of trees and stays there. The visible symptom is
under-coverage of the fitted function; the underlying quantity is effective
sample size. On the one average-case problem with a published BART coverage
number, the worst of that cell's 25 evaluation points has an effective sample
size of about two draws out of 2500 kept, in every chain.

The package now has an evaluation battery of test problems with an accept rule
fixed before any run, and five candidate kernels were measured against it. Four
earn no share of the per-tree draw over move types: swap, a same-variable cut
move killed on its registered statistic, an exchange of trees between chains
killed on its acceptance rate, and an exact draw on the forest's level fibre,
the set of leaf-value shifts that leave the fitted function exactly unchanged,
which does nothing measurable while the tree moves run and was kept only for a
frozen forest. One succeeded: replacing the Metropolis change proposal at a
node whose two children are both leaves with an exact draw from that rule's
full conditional. It raises effective sample size on the average-case cell by
1.8 to 3.0 times the threshold registered there, at both blocks of twenty seeds
and on both of the cell's mean functions, and a restricted version that redraws
only the cut keeps about half that gain for a thirtieth of the added
arithmetic.

**The decision.** On 2026-09-08 the maintainer adopted that restricted draw at
a mixture weight of `d` = 0.16, sixteen percent of the per-tree draw over move
types, taken out of the change move's share, which falls from 0.40 to 0.24. It
lands after the first release; until then the shipped kernel is unchanged and
every candidate stays at weight zero. Adoption still owes what section 6 sets
out: the kernel wins no pathology, three of the four average-case cells do not
exist, one secondary regression stands, and the harm check required before a
default flips has nowhere to run.

## 2. The coverage deficit

The literature's cleanest average-case number is He and Hahn's factorial study:
95 percent pointwise coverage of the true mean function of 0.73 to 0.74 at ten
thousand rows, thirty predictors and moderate noise. dbarts reproduces that
deficit in direction and about half in size, 0.82 against a nominal 0.95 at the
setting closest to the paper's, at intervals 15 to 25 percent longer than the
published BART's for the same point accuracy.

Two levers could close the gap: more trees, or better mixing. The tree-count
lever works, 200 trees against 75 buying about ten points of coverage on both
mean functions at no cost in error. It was parked anyway, meaning kept as a
comparison and not adopted as a default, because it is not one-sided: four
pooled chains at 75 trees read 0.961 on the first mean function where one chain
at 200 trees reads 0.922, and 0.895 on the second where that chain reads 0.924.

The measurement that parked it is three arms on twenty matched seeds. An arm is
one sampler configuration; every arm in a contrast runs on the same seeds as
the others, so differences are paired.

| arm | sweeps, total over four chains | 95% coverage |
|---|---|---|
| four chains, 500 burn-in and 500 kept | 4000 | 0.961 |
| four chains, 1000 burn-in and 2500 kept | 14000 | 0.959 |
| one chain, 1000 burn-in and 25000 kept | 26000 | 0.902 |

The four short chains disagree enough that pooling them is what supplies the
interval's width, while one chain run ten times longer still under-covers and
its per-chain minimum effective sample size stays at 2 to 3. The posterior is
right and the sampler is slow. The maintainer ruled on that reading: the
coverage deficit is a mixing symptom, mixing is the lever, and the statistic a
change has to move here is per-chain effective sample size, not coverage, which
has no headroom at 0.961.

The deficit has two channels. Freeze the cell's tree structures, so only leaf
values and the residual scale are still drawn, and at the median evaluation
point effective sample size rises from 14.9 to about 670 of 2500 kept while its
lag-one autocorrelation falls from 0.70 to between 0.34 and 0.40: at a typical
coordinate the deficit is in the structural channel, which tree structures the
chain visits. At the worst coordinate it is not, the minimum rising only from
1.6 to between 4 and 21, and there the deficit is in the leaf channel, which
leaf values the chain draws given a structure. Nor is that coordinate an
acceptance problem: among the proposals this cell actually evaluates, 24.9
percent are accepted, against 8.0 percent in a low-noise cell that mixes worse.
So a structural kernel cannot be the whole answer.

## 3. The shipped tree kernel

The engine holds five structural tree moves: birth and death, change, swap,
perturb, and the nog-node rule draw `rule_gibbs`, a nog node being an interior
node whose two children are both leaves. One uniform draw per tree per sweep
selects among them, at the default mixture `birth_death 0.6, swap 0, change
0.4, perturb 0, rule_gibbs 0`. So the kernel a user runs is birth and death at
the tree's fringe, the leaf pairs at the bottom, plus change, which redraws an
interior node's split variable and cut while leaving the subtree beneath that
node in place. The other three ship at weight zero and are reachable through
`proposal.probs`, the R argument carrying the mixture; all five at zero freezes
the forest while leaf values, the residual scale and any latent variables the
response family carries keep sampling.

One leaf-value step sits beside them. `levelGibbs` adds a constant to every
occupied leaf of a tree, the constants summing to zero across the forest, so
the fitted function is exactly unchanged; those shifts are the level fibre. Its
slot takes three values: `TRUE` and `FALSE` force it, and the default `NA`
resolves per forest and per sweep, running the step exactly where that forest's
mixture is frozen.

## 4. The battery and the accept rule

The battery is a fixed set of test problems, each with a known truth, a named
failure mode and one statistic chosen before the run: four average-case core
cells, coded C1 to C4, and eight pathologies, P1 to P8. Five are built.

| code | problem | what it stresses | statistic |
|---|---|---|---|
| C1 | He and Hahn factorial: 10000 rows, 30 predictors, moderate noise, two mean functions (Trig+poly, a trigonometric polynomial with one true interaction; Single index, a rotated ridge) | a correlated design at a realistic size | 95% pointwise coverage of the true mean function, on 1000 held-out rows |
| P1 | low-noise Friedman emulator: 2000 rows, sigma 0.25 | structure freezes as the noise falls | 90% coverage; per-move acceptance rate |
| P2 | confounded step function: 300 rows, 3 columns, one tree | two exactly equiprobable representations of one fit | between-chain sd of the fraction of draws with the root split on x1, against a 0.5 null |
| P5 | checkerboard on an autocorrelated design: 1600 rows, 40 columns | a two-way interaction with ambiguous inclusion | between-chain sd of inclusion on the four true columns and their near neighbours |
| P6 | diagonal shelf with targeted selection: 250 rows, two covariates and a treatment | a rotated boundary plus confounding | treatment-effect bias and 95% coverage; published BART reads 0.27 and 65% |

The accept rule covers only C1's first mean function; Single index is reported
beside it and is called the ungated function below.

The rule is asymmetric. The core cells are a gate: a change must not regress
any of them, and improving one earns it nothing. The pathologies are where a
change has to win, by beating at least one on that cell's own registered
statistic by more than four times the measured per-replicate standard error.
Two absolute gates sit on top: a kernel added at weight zero must be bitwise
identical to the control, the arm at the shipped mixture; and P1's 90 percent
coverage must come back near 0.71 in the control arm, or the harness is broken
and no verdict from any cell is valid. It reads 0.725.

### 4.1 Margins, and what stands in for wall time

A secondary metric is flagged, counted as a regression, when the paired mean
difference over twenty seeds is worse than that metric's fixed margin and its
one-sided 95 percent bound also excludes the margin: both conditions, so a
noisy cell cannot flag on a point estimate alone. Margins are per cell,
Holm-corrected within a metric, and a flagged cell is re-run on a fresh block
of twenty seeds before the flag counts.

| metric | margin: worse than this fails |
|---|---|
| 95% pointwise coverage of the true mean function | -0.010 absolute |
| held-out RMSE against the true mean function | ratio above 1.02 |
| minimum ESS over 25 fixed points, per second | ratio below 0.90 |
| summed inclusion share on the true columns | -0.010 absolute |
| wall time per sweep | ratio above 1.05 |
| outer estimand RMSE | ratio above 1.02 |
| outer estimand interval coverage | -0.010 absolute |
| ESS of the outer group scale | ratio below 0.90 |

The last three govern causal and embedded cells, none of them built; an outer
estimand is the quantity an enclosing model reports, a treatment effect or a
group-level scale, as opposed to BART's own fit.

Two of these margins, the per-second ESS ratio and wall time per sweep, are on
quantities the program never measured. No host was quiet, so no timing taken
here would mean anything, and a kernel priced at more than twice a sweep's
arithmetic would fail both by construction whatever it bought. The rule draw's
design therefore replaced them, before its run, with a count of cut scans as
the cost conjunct of its kill criterion, and with one equal-cost arm, which is
where the per-second question is answered. A cut scan is one pass over a node's
members for a single variable, so a full pass over the data is about L such
units at L leaves per tree; costs are quoted as sweep-equivalents, multiples of
a sweep's own three passes. Every cost figure below is a scan count.

### 4.2 The improvement bar, the sham arm and the reference arm

The primary statistic on C1 is the summed minimum effective sample size: within
a chain, the smallest effective sample size over the cell's 25 fixed evaluation
points, added across the four chains. A change has to raise it by +8, four
times that cell's paired standard error of 2.0. That threshold is the +8 bar
below.

A sham arm checked the bar rather than setting it: the control against itself
at fresh sampler seeds, so its paired difference should be zero. It reads -2.3
+/- 9.7. Every "mean +/- sd" here is the mean of the twenty paired differences
and their standard deviation across seeds, so the sham's paired standard error
is 2.17 and four times it is 8.7; the +8 bar is therefore marginally
optimistic. The sham has not been re-run, and every arm below rests on that one
reading.

Coverage on C1 is not read against the control. The shipped four short chains
over-cover, so the -0.010 coverage margin is read against a well-mixed
reference arm instead, on the maintainer's ruling of 2026-09-08 after a
candidate's coverage secondary flagged against the control and its fresh-seed
re-run confirmed the flag. Two long-run readings appear in this report and they
are different arms. Section 2's long single chain is the shipped kernel, one
chain of 1000 burn-in and 25000 kept, reading 0.902: one chain run ten times
longer does not fix the mixing. The reference arm is a different kernel at a
different chain count, the exact rule draw of section 5.1 at `d` = 0.32, four
chains of 1000 burn-in and 2500 kept, reading 0.941: that is what a well-mixed
pooled four-chain interval covers at. Pooling four chains that disagree is what
adds the extra width, so the two readings are consistent.

| arm, C1 Trig+poly | 95% coverage |
|---|---|
| four chains, 500 burn-in and 500 kept, the shipped control | 0.961 |
| the reference arm, four chains, 1000 burn-in and 2500 kept | 0.941 |
| that same reference fit read at half its kept length | 0.939 |
| the same kernel at eight chains | 0.956 |

Chain count, not chain length, is what inflates the pooled interval.

The reference is credible on three grounds. Its between-chain ratio is 0.48,
the lowest recorded here, against the control's 0.78; that ratio is, at each of
the 25 points, the standard deviation across chains of a chain's posterior mean
divided by the pooled posterior standard deviation, then the median over
points, near zero when the chains agree and near one when each sits in its own
place. Its two lengths agree. And on a cell with no observation weights masked
out and no missing values, the kernel producing it is an exact Gibbs step on
the posterior the shipped kernel already targets. What it is not is independent
of what it judges: no kernel but the candidate mixes well enough here to serve,
which is the finding itself. Residual disagreement still widens a pooled
interval, so 0.941 is an upper bound. The same substitution gives the ungated
function's held-out error a reference ratio of 1.024.

## 5. The kernels tried

### 5.1 The full exact rule draw

This is the only kernel that moved the structural channel under the shipped
sampler, and it moved it by 2.7 and 2.8 times the +8 bar.

Two counts justified building it, both from a census build: an instrumented
compile that logs one line per structural proposal and consumes no random draw,
so it reproduces the uninstrumented chain exactly. Nog nodes are 48.5 to 98.3
percent of interior nodes and take 62.7 to 99.1 percent of change's proposals,
so they are no corner case. And the probability the node's own conditional puts
on the rule already in place is 0.0017 on C1, at a median rank of 26.5 among
3000 candidate rules, so there is much to gain where it matters; on the
low-noise cell that probability is 0.737 and there is almost nothing.

At a nog node the whole conditional over split rules costs one scan and that
scan's marginal is exact, the two children being a two-way partition of the
node's members with nothing below them. The candidate set does not depend on
the rule currently in place, so no proposal count survives into the acceptance.
In the ordinary case, no weights masked out and no missing value routed at the
node, the step is an exact Gibbs draw at acceptance one; where a mask or a
missing value puts the incumbent rule in a lower stratum of the veto that
forbids empty leaves, it is Metropolis-within-Gibbs, also at acceptance one and
valid because the stratum it enters is absorbing. Correctness rests on the two
absolute gates of section 4, bitwise identity with the control at weight zero
and P1's control reading, and on a prior-only detailed-balance script carrying
two poisons, deliberately broken variants that the script must reject, which is
what shows it can detect an error at all. Both fail as designed.

The benefit run took C1's four-chain configuration at `d` = 0.16, twenty
matched seeds and then a fresh block of twenty. Absolute readings on Trig+poly,
seeds 1 to 20, four chains of 500 burn-in and 500 kept; ESS figures are draws.

| statistic | control | rule draw at `d` = 0.16 |
|---|---|---|
| summed minimum ESS | 15 | 36 |
| per-chain minimum ESS | 2 | 3 |
| 95% coverage | 0.961 | 0.939 |
| interval length | 4.61 | 3.98 |
| between-chain ratio | 0.78 | 0.58 |

Paired differences against the control on the same seeds.

| statistic | seeds 1 to 20 | seeds 21 to 40 |
|---|---|---|
| summed minimum ESS | +21.5 +/- 12.8, t 7.5 | +22.1 +/- 12.2, t 8.1 |
| per-chain minimum ESS | +1.14, t 8.7 | +0.86, t 7.2 |
| 95% coverage | -0.022, t -12.5 | -0.026, t -9.9 |
| held-out RMSE, ratio | 0.979 | 0.987 |

The per-chain minimum moves with the sum, which no four-chain arm on this mean
function had achieved before. Doubling `d` to 0.32 buys +37.0 and nothing else.
The arm costs 2.21 sweep-equivalents. On the ungated function it gains +14.4
and +23.9 and pays 1.028 and 1.023 in held-out error, past the 1.02 margin at
both blocks. Its coverage secondary flagged and the fresh-seed re-run confirmed
it; against the reference arm it reads -0.002, inside the margin, because what
the kernel narrows is width unfinished mixing was supplying. The ungated error
regression goes the same way, 1.028 against the control being 1.004 against
that arm.

### 5.2 The cut-only restriction

The variable axis costs thirty times what the cut axis costs, so a restricted
variant holds the incumbent variable and enumerates its cuts alone: one scan
per proposal instead of thirty. It is still an exact Gibbs step, it sits behind
a compile-time switch off in every shipped build, and the balance script passes
with both poisons failing. C1 at `d` = 0.16; the paired rows give the first
seed block then the fresh one.

| | full draw | cut-only |
|---|---|---|
| cut scans per sweep | 337.32 | 10.97 |
| cost, sweep-equivalents | 2.21 | 1.04 |
| Trig+poly summed minimum ESS, paired | +21.5, +22.1 | +11.3, +9.1 |
| Trig+poly coverage, paired | -0.022, -0.026 | -0.007, -0.012 |

The variable axis costs thirty times the scans for about half of the gated
function's gain, and accounts for most of the coverage movement that flagged
the full draw.

**Dose response.** A dose is a value of `d`. Three were run, each over twenty
matched pairs with a fresh block: a half dose at 0.08, and the full draw's own
0.16 and 0.32.

| paired against the control | `d` = 0.08 | `d` = 0.16 | `d` = 0.32 |
|---|---|---|---|
| Trig+poly summed minimum ESS, seeds 1 to 20 | +8.1 | +11.3 | +13.4 |
| Trig+poly summed minimum ESS, seeds 21 to 40 | +1.8 | +9.1 | +9.3 |
| Single index summed minimum ESS, seeds 1 to 20 | +8.3 | +16.2 | +27.5 |
| Trig+poly held-out RMSE ratio, seeds 1 to 20 | 0.999 | 0.992 | 1.057 |
| Trig+poly held-out RMSE ratio, seeds 21 to 40 | 0.994 | 1.012 | 1.048 |
| cost, sweep-equivalents | 1.02 | 1.04 | 1.08 |

Only the 0.16 scan count is measured; the other two scale it linearly in the
share at a fixed tree population, an assumption and not a measurement. `d` =
0.08 does not survive its fresh block, at t 1.29 against the +8 bar. Coverage
flags nowhere: against the reference arm every Trig+poly dose sits 0.011 to
0.015 above 0.941. But 0.32 regresses held-out error at both blocks, and the
restriction is why. The kernel's share comes out of change, the only move that
changes a node's split variable, so at 0.32 change is left with 0.08; the
interval widens, the chains agree no better and the fit is worse, where the
full draw at that dose does the opposite because it carries the variable axis
itself. Change is load-bearing elsewhere too: dropping it costs P5 0.027 of the
true columns' inclusion share and 3.9 percent of held-out error on fresh seeds,
and stops P2's root variable moving at all. On the ungated function the
cut-only kernel is monotone at held-out ratios of 1.033, 1.033 and 1.032, past
the 1.02 margin at every dose, and 1.009 against that function's reference
ratio of 1.024.

### 5.3 The level-fibre step

Frozen, this step is the largest effect in the program; alongside the tree
moves it is nothing. Because it leaves the fitted function exactly unchanged,
its conditional is the leaf prior alone and the draw is closed form, at about a
five-thousandth of a sweep. With C1's structures frozen and only leaf and
residual-scale draws running, the minimum effective sample size over the cell's
25 points rises by a paired median of +189.4 at one freeze point and +231.6 at
the other, all ten pairs positive. Live, the paired differences are these.

| paired difference | seeds 1 to 20 | fresh seeds |
|---|---|---|
| Trig+poly summed minimum ESS | -0.9 | -2.5 |
| Single index summed minimum ESS | -5.3, t -2.90 | +1.7 |
| Trig+poly held-out RMSE, ratio | 1.012 | 1.014 |
| Trig+poly between-chain ratio | 0.79 against 0.78 | 0.79 against 0.80 |

That is against the +8 bar and inside the sham arm's own reading, with coverage
flat to the third digit. Two adverse readings were taken and neither stood: the
Single index loss separates from noise on the first block but the fresh block
reads +1.7, and a wall-time ratio of 1.086 was taken with three arms
interleaved on a loaded host, so the two arms were re-measured alone on a quiet
one and read 0.986. Run together with the nog-node rule draw, the step adds
nothing: every column matches that kernel on its own. KILLED as a general
default, the explanation being the frozen run's scope: it measures one channel
with the structure held fixed, and when the tree moves run they move the level
fibre faster than the exact draw pays for itself. What was kept is that one
regime, through the automatic `levelGibbs` default of section 3; a non-frozen
forest draws nothing, so the shipped engine is bitwise unchanged.

### 5.4 Perturb, the same-variable cut move

Change always redraws the split variable, so a pure cut displacement happens
only when the redraw lands back on the variable already there. The census made
the case for building it: change accepts 1.7 to 6.1 percent of its proposals by
cell, against 27 to 51 percent for a one-position cut displacement. The move
landed at weight zero, bitwise neutral, behind a detailed-balance script whose
two poisons both fail as designed. Its benefit run on C1 at `d` = 0.16 did not
move the gated statistic: +0.1 +/- 8.1 on the first block and -1.9 +/- 6.5 on
the fresh one, against the +8 bar and inside the sham arm's own reading, with
held-out error at 1.032 against the 1.02 margin. KILLED at that setting. On the
ungated function it gains +9.5 summed on the fresh block, 19 of 20 seeds at t
5.35, carrying its own 1.027 error regression; that gain is unclaimed, the
accept rule not covering that mean function. The kernel stays at weight zero.

### 5.5 Swap

Swap exchanges the rules of a parent and one child, and its evidence is not a
battery cell. A census priced it as mostly wasted work, and on the one
criterion where the shipped mixtures separate, re-adaptation after the response
is swapped under the trees, change with swap at zero matches the default on
every contrast. Swap's share of 0.1 moved to birth and death and the move was
deleted, then restored at a default of zero on an exact-posterior gate: one
tree, two live columns, the fitted quantity against a brute-force enumeration
of all 62 reachable trees, where the largest absolute gap in the tree
probabilities is 0.0120 without swap and 0.0008 with it, against a tolerance of
0.004. Swap alone rotates a child's rule up the tree, which is what the
one-tree gate sees; at fifty and two hundred trees the ensemble averages the
effect away, so no default moved.

### 5.6 Cross-chain exchange

KILLED on acceptance. An exchange of one tree between two chains at one
temperature is ordinary Metropolis on the product target, the two tree priors
cancelling, so its acceptance rate is closed form; it was evaluated on states
taken from a running sampler, the move itself never proposed.

| acceptance rate | |
|---|---|
| at C1's two seeds | 10.1 to 11.7 percent |
| at the low-noise cell | 4.4e-8 to 2e-4 |
| by leaf-count pair, one leaf against two | 0.331 |
| by leaf-count pair, five-plus against five-plus | 4e-22 |

The move lives on the one- and two-split trees birth and death already reach,
and vanishes on the multi-split trees it was sized for.

## 6. The decision and what it rests on

**What was decided.** On 2026-09-08 the maintainer adopted the cut-only exact
rule draw at `d` = 0.16, to land after the first release. It is the only
cut-only dose that clears the +8 bar at both seed blocks, +11.3 and +9.1, with
every gated secondary clean at the reference-read margins. Its one adverse
secondary is the ungated function's held-out error at 1.033, past the 1.02
margin against the control and 1.009 against that function's reference ratio.

**Cost.** An equal-cost arm gives the shipped kernel the extra sweeps a
competitor's cost ratio buys, so the two run at one budget. The arm that was
run tested the full draw, not the adopted kernel: four chains of 1105 burn-in
and 1105 kept against the shipped 500 and 500, at the full draw's 2.21.

| arm, at a 2210 cut-scan-unit budget | summed minimum ESS |
|---|---|
| the full rule draw, 1000 sweeps at 2.21 | 36.3 |
| the shipped kernel, 2210 sweeps at 1.00 | 22.9 |

Length buys 63 percent of what the kernel buys. It also moves coverage from
0.961 to 0.952 with no kernel change at all, a third of the way to the
reference, while barely moving what the kernel moves: between-chain 0.73 from
0.78, against the kernel's 0.58. So on coverage the kernel and simply running
longer are hard to tell apart, and the case for the kernel is a cost case, not
a coverage one.

A budget below is a run's sweep count times its cost in sweep-equivalents, so
the control's 1000 sweeps at 1.00 are 1000 units and its rate per thousand
units is numerically its own summed minimum ESS.

| arm | budget, cut-scan units | summed minimum ESS | ESS per thousand units |
|---|---|---|---|
| the control | 1000 | 14.8 | 14.8 |
| cut-only at `d` = 0.16 | 1040 | 26.1 | 25.1 |
| the full draw | 2210 | 36.3 | 16.4 |
| the equal-cost arm | 2210 | 22.9 | 10.4 |

The adopted kernel leads that ranking, and the ranking is weaker than it looks:
it never ran at equal cost, and this statistic grows sub-linearly in chain
length, so a per-unit ratio across four different budgets favours the smallest
budget by construction. The design record bridges the two kernels by
extrapolation instead: on the equal-cost arm's measured scaling the cut-only
kernel at the full draw's budget would read about 40 against 36.3, extrapolated
off a scaling measured on the shipped kernel and not this one. Cut-only at 0.32
does not match the full draw at 0.16, so the choice is not closed by dominance.

**What the kernel has not earned.** Under the accept rule as written, nothing.
The rule wants a pathology win and treats the core as a gate, and every gain
recorded for the rule draw is on C1, a core cell. P2 ran as a must-not-degrade
control, and the move cannot reach the chains that get stuck there. P1 was
re-run as the absolute gate and carries no rule-draw arm; nor do P5 and P6. And
the core gate is checkable on one cell out of four.

**The departure the case rests on.** The kernel's kill criterion was registered
with departures from the program's own, two of which matter: the statistic is
minimum effective sample size rather than coverage, and the cell is C1 rather
than the low-noise cell, on which three shipped mixtures are indistinguishable.
Both follow the ruling of section 2. Adoption accepts it: the core cell's
mixing statistic becomes the target rather than a gate, and the pathologies
stay untouched.

**What adoption owes, all of it after the release.** First, the design
amendment that chooses the surface, the restricted kernel being a private
compile-time build today and not a mode anything can select; the choice is
between deleting the variable axis from the rule draw and adding a sixth name
to `proposal.probs`, twenty-four files for the last name added. Second, the
harm controls listed in the kernel's own design and never run for the
restricted draw. Third, the plateau-error gate, a gap and not a plan: the house
rule wants a per-cell check that posterior-mean prediction error has not
worsened once the sampler has reached its plateau, measured in a noise-heavy or
a large-n cell, and no built cell is either, so it has nowhere to run. Fourth,
the default flip and its baselines: any nonzero share moves every draw from the
first sweep of the first chain, so the three equivalence baselines, the counts
computed from them, the gate ledger's baseline footnote and the four
seeded-drift snapshots are regenerated, with the consumer packages on their
lockstep branches re-recording alongside. That work is mechanical, every
regenerated baseline being checked against a build that reproduces the draws
exactly.

### 6.1 A standing rule

No exploration kernel is removed from the engine before this report has been
read. A kill leaves the kernel in place at weight zero rather than deleting it.

## 7. What was not measured

- **Wall time, and anything per second.** No host was quiet; the one-minute
  load ran from 4 to 269 across the arms and 100 to 122 on the equal-cost
  machine, so every cost figure here is a scan count.
- **The plateau-error gate.** No built cell is noise-heavy or large-n, so the
  house rule's harm check has nowhere to run.
- **The sampler's own speed since the default last moved.** A bench-sampler
  comparison on a quiet machine, owed since swap's 0.1 went to birth and death.
- **The adopted kernel's own controls.** No P1 control reading, no P2 arm, no
  sham arm, no fresh seed block on the ungated function and no equal-cost arm,
  all of which the full draw's run carried. Its cut-scan figures come from a
  census cell at its own data and seed on one chain of 200 burn-in and 500
  kept, not from a replay of the arms' twenty seeds, and that census cell takes
  its 0.16 share out of change and out of birth and death together, where the
  benefit arms take it out of change alone. Neither difference was corrected
  for in the scan count.
- **Any rule-draw arm on a pathology.** P1, P5 and P6 carry none.
- **A swap share in five one-tree exact gates.** Every one-tree exact gate that
  accepts a caller mixture now sets a positive swap share; four two-forest
  causal gates and one monotone gate cannot, so those five run without one.
  That is the one place swap's argument is not covered.
- **Response-swap recovery inside the battery.** Measured once on its own grid,
  where it decided swap's default; never re-run against any new kernel, its
  cell unbuilt.
- **Whether composing BART with a parametric block helps tree-space mixing.**
  Unmeasured in the literature and here; the cost of finding out is a factor of
  six.
- **Seven cells.** The real-data convergence ladder (P3), the
  inhomogeneous-smoothness cell (P4), the no-overlap extrapolation cell (P7),
  the hierarchical variance funnel (P8), the mixed-type regime switch (C2), the
  real-covariate causal cell (C3) and the embedded moving-response cell (C4).
  The funnel and the embedded cell each need an outer sampler from another
  repository.
- **Non-gaussian families.** Every built cell is gaussian or causal gaussian,
  and coverage of a true mean function calibrates a sampler, it does not prove
  it correct.

## Appendix A. Refuted and unbuilt proposals

Two rounds of proposal generation ran under a common bar: a written
Metropolis-Hastings correction, a price in cut scans, a named deficit and a
falsifier runnable in a day. Tree-space geodesics, a learned rotation and a
per-tree temperature each fail on their own terms. A per-leaf per-variable
histogram cache was refuted on size, the saving being a depth-fold, about two
at the measured tree size, not the leaf-fold claimed.

The move census refuted two more mechanisms. An informed death proposal,
choosing which leaf pair to prune by weight rather than uniformly, failed its
own kill criterion: the weights are effectively a point mass, which makes the
proposal the uniform one. A lifted cut displacement, which would have given a
cut move a persistent direction, was refuted as a source of gain, accepted
displacements reversing rather than continuing.

Two constructions remain unbuilt with their arguments intact. A
pairwise-collapsed split transfer is the only candidate addressing
representation multimodality directly; its probe was not built because pricing
a transfer needs the residual net of the other trees and the partner's leaf
statistics, which no existing move sees. A lifted birth and death is cheap and
valid but recomputes at a gain of only 1.06 to 1.17.

The battery's design law comes from an earlier study. Grow-from-root as a
default was killed in both strata: every aggregate test passed, but per-cell
plateau posterior-mean error costs of +11.10 percent in a noise-heavy small-n
cell and of +4.47 and +10.66 percent in a large-n one, each past its frozen
margin and each confirmed on fresh seeds, were averaged away by pooling. What
the battery inherits is that law: per-cell checks, thresholds frozen before the
run, mandatory fresh-seed re-runs and a null control that voids the family.

## Appendix B. Source records

Every claim above is recorded in one of these sections.

| what it holds | record |
|---|---|
| The C1 cell, and every arm run on it | [10.4 C1, the He and Hahn factorial](../design/benchmark-surfaces.md#104-c1-the-he-and-hahn-factorial) |
| What the battery's pilot could not build | [10.6 What the pilot could not do](../design/benchmark-surfaces.md#106-what-the-pilot-could-not-do) |
| P1's control reading of 0.725 | [10.8 P1, the low-noise Friedman emulator (2026-09-07)](../design/benchmark-surfaces.md#108-p1-the-low-noise-friedman-emulator-2026-09-07) |
| The battery's cells and their statistics | [6. The battery](../design/benchmark-surfaces.md#6-the-battery) |
| The four average-case core cells | [6.2 The average-case core](../design/benchmark-surfaces.md#62-the-average-case-core) |
| The accept rule as a procedure | [6.1 The rule, stated operationally](../design/benchmark-surfaces.md#61-the-rule-stated-operationally) |
| The margins, and the coverage baseline ruling | [6.4 What "no regression on the core" means numerically](../design/benchmark-surfaces.md#64-what-no-regression-on-the-core-means-numerically) |
| What the battery does not cover | [6.6 What this battery does not measure](../design/benchmark-surfaces.md#66-what-this-battery-does-not-measure) |
| The rule draw's construction and correctness law | [2. The move](../design/nog-gibbs.md#2-the-move) |
| Its cost, and the two restricted variants | [2.4 Cost, and the two restricted variants](../design/nog-gibbs.md#24-cost-and-the-two-restricted-variants) |
| Its benefit run, doses, equal-cost arm and cost instrument | [6. Benefit, pre-registered](../design/nog-gibbs.md#6-benefit-pre-registered) |
| Its slices, and the adoption ruling | [8. Slices](../design/nog-gibbs.md#8-slices) |
| The level-fibre step's placement and cost | [2. Placement and cost](../design/level-fibre.md#2-placement-and-cost) |
| Its three-valued surface | [4. The surface](../design/level-fibre.md#4-the-surface) |
| Its frozen and live benefit runs | [6. Benefit, pre-registered](../design/level-fibre.md#6-benefit-pre-registered) |
| The re-record bundle a default share pays | [7. RNG and baselines](../design/level-fibre.md#7-rng-and-baselines) |
| Perturb's benefit run, the +8 bar and the sham arm | [5. Benefit, pre-registered](../design/perturb-move.md#5-benefit-pre-registered) |
| The chain configuration, and the primary statistic | [5.1 The chain configuration, and what it makes the primary statistic](../design/perturb-move.md#51-the-chain-configuration-and-what-it-makes-the-primary-statistic) |
| Swap's removal and the evidence for it | [1. The decision, and its evidence](../design/swap-removal.md#1-the-decision-and-its-evidence) |
| The sampler-speed comparison still owed | [8. Landing](../design/swap-removal.md#8-landing) |
| Swap's return at weight zero, and the default mixture | [9. Reversal: the move returns at default zero](../design/swap-removal.md#9-reversal-the-move-returns-at-default-zero) |
| The move census and its probes | [6.1 Stage 0 - the move census (pilot; no kill criterion)](../design/tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion) |
| The program's kill criteria and the plateau clause | [6.4 Kill criteria, pre-registered](../design/tree-mixing-proposals.md#64-kill-criteria-pre-registered) |
| What the stickiness survey could not settle | [9. What this survey could not settle](../design/tree-mixing-proposals.md#9-what-this-survey-could-not-settle) |
| Recovery after a response swap | [14. Recovery after a response swap (2026-09-06)](../design/tree-mixing-proposals.md#14-recovery-after-a-response-swap-2026-09-06) |
| Proposals refuted in the first round | [15.5 Refuted along the way](../design/tree-mixing-proposals.md#155-refuted-along-the-way) |
| The ranked candidates, built and unbuilt | [16.3 Ranking](../design/tree-mixing-proposals.md#163-ranking) |
| The warm-start study whose design law the battery inherits | [5. Verdict and consequences](../design/grow-from-root-default.md#5-verdict-and-consequences) |
| The five structural tree moves | [Tree moves](../architecture.md#tree-moves) |
| What one sweep does | [One sweep](../architecture.md#one-sweep) |
