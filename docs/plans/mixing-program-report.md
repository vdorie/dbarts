# The tree-sampler mixing program

Status: RECORD, 2026-09-08. The cut-only exact rule draw is adopted, to land
after the first release; no pre-release default changes.

## 1. Summary

BART's tree sampler moves by small local edits, and on hard problems it settles
into one arrangement of trees and stays there. The symptom is under-coverage of
the fitted function; the quantity behind it is effective sample size. On the
one average-case problem with a published BART coverage number, the worst of
that problem's 25 evaluation points has an effective sample size of about two
draws out of 2500 kept, in every chain.

The package now has an evaluation battery with an accept rule fixed before any
run, and five candidate kernels were measured against it, one in two variants.
Four earn no default weight in the mixture of proposal probabilities over tree
moves, drawn once per tree per sweep: swap; a same-variable cut move, killed on
its registered statistic; an exchange of trees between chains, killed on
acceptance; and an exact draw on the forest's level fibre, the leaf-value
shifts that leave the fitted function exactly unchanged, which does nothing
measurable while the tree moves run and was kept only for a frozen forest.

One succeeded: replacing the Metropolis change proposal at a node whose two
children are both leaves with an exact draw from that rule's full conditional.
It raises effective sample size on the average-case problem by 1.8 to 3.0 times
the threshold registered there, spanning both mean functions and both blocks of
twenty seeds. A restricted version redrawing only the cut keeps about half that
gain on the mean function the accept rule covers and all of the full draw's
first-block gain on the other, for a thirtieth of the added arithmetic.

**The decision.** On 2026-09-08 the maintainer adopted that restricted draw at
a mixture weight of `d` = 0.16, taken out of the change move's share, which
falls from 0.40 to 0.24. It lands after the first release; until then the
shipped kernel is unchanged and every candidate stays at weight zero. Adoption
still owes what section 6 sets out: the kernel wins no pathology, three of the
four average-case cells are not built, one secondary regression stands, and the
harm check required before a default flips has nowhere to run.

## 2. The coverage deficit

The literature's cleanest average-case number is He and Hahn's factorial study:
95 percent pointwise coverage of the true mean function of 0.73 to 0.74 at ten
thousand rows, thirty predictors and moderate noise. dbarts reproduces that
deficit in direction and about half in size, 0.82 against a nominal 0.95, at
intervals 15 to 25 percent longer for the same point accuracy.

Two levers could close the gap: more trees, or better mixing. The tree-count
lever works, 200 trees against 75 buying about ten points of coverage at no
cost in error, and it was parked anyway, meaning kept as a comparison and not
made a default, because it is not one-sided: four pooled chains at 75 trees
read 0.961 on Trig+poly, the trigonometric polynomial this problem fits, where
one chain at 200 trees reads 0.922, and 0.895 on Single index, its
rotated-ridge companion, where that chain reads 0.924.

The measurement that parked it is three arms on twenty matched seeds. An arm is
one sampler configuration, and every arm in a contrast runs on the same seeds,
so differences are paired.

| arm | sweeps, total over four chains | 95% coverage |
|---|---|---|
| four chains, 500 burn-in and 500 kept | 4000 | 0.961 |
| four chains, 1000 burn-in and 2500 kept | 14000 | 0.959 |
| one chain, 1000 burn-in and 25000 kept | 26000 | 0.902 |

The four short chains disagree enough that pooling them is what supplies the
interval's width, while one chain run ten times longer still under-covers and
its effective sample size at the worst of the problem's 25 evaluation points
stays at 2 to 3 draws. The posterior is right and the sampler is slow. The
maintainer ruled on that reading: the coverage deficit is a mixing symptom,
mixing is the lever, and the statistic a change has to move here is per-chain
effective sample size, not coverage, which has no headroom at 0.961.

The sampler is not refusing to move, its proposals being accepted at 24.9
percent here against 8.0 percent on a low-noise problem that mixes worse, which
is why the next test froze the tree structures and read the deficit's two
channels apart. Frozen, at the median evaluation point effective sample size
rises from 14.9 to about 670 of 2500 kept, so at a typical coordinate the
deficit is in the structural channel, which tree structures the chain visits.
At the worst coordinate the minimum rises only from 1.6 to between 4 and 21, so
there it is in the leaf channel, which leaf values the chain draws given a
structure. A structural kernel cannot be the whole answer.

## 3. The shipped tree kernel

The engine holds five structural tree moves: birth and death, change, swap,
perturb, and the nog-node rule draw `rule_gibbs`, a nog node being an interior
node whose two children are both leaves. One uniform draw per tree per sweep
selects among them, at the default mixture `birth_death 0.6, swap 0, change
0.4, perturb 0, rule_gibbs 0`. So the kernel a user runs is birth and death at
the tree's fringe, the leaf pairs at the bottom, plus change, which redraws an
interior node's split variable and cut while leaving the subtree beneath that
node in place. The other three ship at weight zero, reachable through
`proposal.probs`, the R argument carrying the mixture; all five at zero freezes
the forest while leaf values, the residual scale and any latents the response
family carries keep sampling.

One leaf-value step sits beside them, not in the mixture. `levelGibbs` adds a
constant to every occupied leaf of a tree, the constants summing to zero across
the forest, so the fitted function is exactly unchanged; those shifts are the
level fibre. Its slot takes three values: `TRUE` and `FALSE` force it, and the
default `NA` resolves per forest and per sweep, running the step exactly where
that forest's mixture is frozen.

## 4. The battery and the accept rule

The battery is a fixed set of test problems, each with a known truth, a named
failure mode and one statistic chosen before the run: four average-case core
cells, C1 to C4, and eight pathologies, P1 to P8. Five are built.

| code | problem | what it stresses | registered statistic |
|---|---|---|---|
| C1 | the He and Hahn factorial of section 2, on two mean functions: Trig+poly, with one true interaction, and Single index, a rotated ridge | a correlated design at a realistic size | 95% coverage of the true mean function, on held-out rows |
| P1 | a low-noise Friedman emulator | structure freezes as the noise falls | 90% coverage; acceptance rate |
| P2 | a confounded step function, fitted with one tree | two exactly equiprobable representations of one fit | between-chain sd of the root-on-x1 fraction, against a 0.5 null |
| P5 | a checkerboard on an autocorrelated design | a two-way interaction with ambiguous inclusion | between-chain sd of inclusion on the four true columns |
| P6 | a diagonal shelf with targeted selection | a rotated boundary plus confounding | treatment-effect bias and coverage; bias 0.27 and coverage 65 percent in the published BART |

C1's registered statistic is coverage. Section 4.2 explains why the rule-draw
arms were judged on effective sample size instead, and section 6 treats that
substitution as the departure it is.

The rule covers only Trig+poly; Single index is reported beside it and is
called the ungated function below. It is asymmetric: the core cells are a gate,
so a change must not regress one and gains nothing by improving one, while the
pathologies are where it has to win, by beating one on that cell's statistic by
four times the measured per-replicate standard error. Two absolute gates sit on
top: a kernel added at weight zero must be bitwise identical to the control,
the arm at the shipped mixture, and the low-noise cell's 90 percent coverage
must come back near 0.71 in the control arm or no verdict is valid; it reads
0.725.

### 4.1 Margins, and what stands in for wall time

Each cell has one primary statistic, the one a change has to move, and
secondary metrics that must not regress. A secondary is flagged, counted as a
regression, only when the paired mean difference over twenty seeds is worse
than that metric's margin and its one-sided 95 percent bound also excludes it;
a flagged cell is then re-run on a fresh block of twenty seeds. Three margins
carry the verdicts below: coverage of the true mean function at -0.010
absolute, held-out RMSE against it at a ratio of 1.02, and inclusion share on
the truly relevant columns at -0.010 absolute. The rest are in the record.

Two further registered margins are per-second quantities, minimum effective
sample size per second and wall time per sweep, and no host was ever quiet, so
neither was read. In their place the design put a count of cut scans, one scan
being a single pass over a node's members for one variable. Costs are quoted as
sweep-equivalents, multiples of a sweep's own three passes over the data, and
one equal-cost arm was run to answer the per-second question directly.

### 4.2 The improvement bar, the sham arm and the reference arm

The primary statistic on C1 is the summed minimum effective sample size: within
a chain, the smallest effective sample size over the cell's 25 fixed evaluation
points, added across the four chains. A change has to raise it by +8, four
times that cell's paired standard error; that is the +8 bar below. A sham arm,
the control against itself at fresh sampler seeds, reads -2.3 +/- 9.7 where it
should read zero; its own paired standard error is 2.17, four times which is
8.7, so the +8 bar is marginally optimistic. That arm was run once, and every
verdict below rests on the one reading. A "mean +/- sd" in this report is the
mean of the twenty paired differences and their standard deviation across
seeds.

Coverage on C1 is read not against the control but against a well-mixed
reference arm, and the substitution was made after the fact: the rule draw's
coverage secondary flagged against the shipped control, the mandatory
fresh-seed re-run confirmed the flag, and only then was the margin's baseline
replaced, on the maintainer's ruling of 2026-09-08, the shipped four short
chains having been measured to over-cover. That reference arm is the exact rule
draw of section 5.1 at `d` = 0.32, four chains of 1000 burn-in and 2500 kept:
the control reads 0.961 and the reference 0.941, and at one kernel and one
length doubling the chain count moves coverage from 0.939 to 0.956 where
quintupling the chain length moves it only from 0.939 to 0.941, so it is chain
count and not chain length that inflates the pooled interval.

Section 2's long single chain is a third arm, and the two are easy to confuse.
That one is the shipped kernel, one chain of 1000 burn-in and 25000 kept,
reading 0.902, which says that running one chain ten times longer does not fix
the mixing; the reference's 0.941 is what a well-mixed pooled four-chain
interval covers at.

The reference is credible on three grounds. Its between-chain ratio is 0.48
against the control's 0.78, the lowest recorded here, that ratio being the
median over the 25 points of the across-chain standard deviation of a chain's
posterior mean over the pooled posterior standard deviation. Its reading has
converged in length, the same fit read at half its kept draws giving 0.939
against 0.941. And the kernel producing it is, on this cell, an exact Gibbs
step on the posterior the shipped kernel targets. What it is not is independent
of what it judges, and residual disagreement still widens a pooled interval, so
0.941 is an upper bound; the same substitution gives the ungated function's
held-out error a reference ratio of 1.024.

## 5. The kernels tried

### 5.1 The full exact rule draw

This is the only kernel that moved the structural channel under the shipped
sampler, and it moved it by 2.7 and 2.8 times the +8 bar.

A census build made the case for it. That is an instrumented compile which logs
every structural proposal without consuming a random draw, so a run reproduces
its arm's chain exactly. It shows that nog nodes are most of a tree's interior
nodes and take most of change's proposals, and that the probability the node's
own conditional puts on the rule already in place is 0.0017 on C1 against 0.737
on the low-noise cell (P1): much to gain where it matters.

At a nog node the whole conditional over split rules costs one scan, that
scan's marginal is exact, and the candidate set does not depend on the rule
currently in place, so the step is an exact Gibbs draw at acceptance one. Where
a masked observation weight or a routed missing value strands part of the
candidate set and the step enters a strictly better stratum of the veto that
forbids empty leaves, it is Metropolis-within-Gibbs, also at acceptance one,
valid because the better stratum is absorbing; no stationarity is claimed for
the stratum it leaves. Correctness rests on two gates of the kernel's own: an
assertion that its candidate set and scores agree candidate by candidate with
an independently written reference assembly, and a prior-only detailed-balance
script carrying two poisons, deliberately broken variants the script must
reject, both of which fail as designed. Separately, at weight zero the kernel
is bitwise identical to the control, which is the battery's own first absolute
gate.

The benefit run took C1's four-chain configuration at `d` = 0.16, twenty
matched seeds then a fresh block of twenty; Trig+poly, four chains of 500
burn-in and 500 kept. Control and rule-draw columns are absolute readings on
seeds 1 to 20; ESS is in draws.

| statistic | control | rule draw | paired, seeds 1 to 20 | paired, seeds 21 to 40 |
|---|---|---|---|---|
| summed minimum ESS | 15 | 36 | +21.5 +/- 12.8, t = 7.5 | +22.1 +/- 12.2, t = 8.1 |
| per-chain minimum ESS | 2 | 3 | +1.14, t = 8.7 | +0.86, t = 7.2 |
| 95% coverage | 0.961 | 0.939 | -0.022, t = -12.5 | -0.026, t = -9.9 |
| held-out RMSE, ratio | - | - | 0.979 | 0.987 |

The per-chain minimum moves with the sum, which no four-chain arm here had done
before. Doubling `d` to 0.32 buys +37.0 and moves no other gated reading. At
`d` = 0.16 the full draw costs 2.21 sweep-equivalents. On the ungated function
it gains +14.4 and +23.9 but pays 1.028 and 1.023 in held-out error, past the
1.02 margin at both blocks. Read against the reference arm rather than the
over-covering control, both regressions dissolve, coverage to -0.002 and the
ungated error to 1.004, because what the kernel narrows is width unfinished
mixing was supplying.

### 5.2 The cut-only restriction

The variable axis costs thirty times what the cut axis costs, so a restricted
variant holds the incumbent variable and enumerates its cuts alone. It is still
an exact Gibbs step, it sits behind a compile-time switch off in every shipped
build, and the balance script passes with both poisons failing. On C1 at `d` =
0.16 it costs 1.04 sweep-equivalents against the full draw's 2.21, and it keeps
about half the gated function's gain. So the variable axis buys the other half
of that gain for thirty times the scans, and it accounts for most of the
coverage movement that flagged the full draw.

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
0.08 does not survive its fresh block, at t = 1.29 against the +8 bar, and
coverage flags nowhere. But 0.32 regresses held-out error at both blocks, and
the restriction is why: the kernel's share comes out of change, the only move
that changes a node's split variable, so at 0.32 change is left with 0.08 and
the fit gets worse, where the full draw at that dose improves it. Change is
load-bearing elsewhere too: dropping it costs the checkerboard cell (P5) 0.027
of the true columns' inclusion share against a margin of 0.010, and stops the
confounded step function's root variable moving at all. On the ungated function
the cut-only kernel is past the 1.02 margin at every dose, 1.033 at worst, and
1.009 against that function's reference ratio.

**The other four kernels.** The level-fibre step, perturb, swap and the
cross-chain exchange earned no weight. Each verdict and the evidence behind
it are in Appendix A.

## 6. The decision and what it rests on

**What was decided.** On 2026-09-08 the maintainer adopted the cut-only exact
rule draw at `d` = 0.16, to land after the first release. It is the only
cut-only dose that clears the +8 bar at both seed blocks with every gated
secondary clean: it reads +11.3 and +9.1 at coverage of 0.955 and 0.952 against
a reference of 0.941 that this same kernel produced at twice the dose, while
0.08 fails the bar on its fresh block and 0.32 regresses held-out error at
both. Its one adverse secondary is the ungated function's held-out error at
1.033, past the 1.02 margin against the control and 1.009 against that
function's reference ratio.

**Cost.** An equal-cost arm gives the shipped kernel the extra sweeps a
competitor's cost ratio buys, so the two run at one budget; the arm that was
run tested the full draw, not the adopted kernel. At one budget the full draw
reads 36.3 summed minimum ESS against the shipped kernel's 22.9, so length buys
63 percent of what the kernel buys. Length also moves coverage from 0.961 to
0.952 with no kernel change, while barely moving the between-chain ratio, 0.73
from 0.78 against the kernel's 0.58. On coverage the kernel and running longer
are hard to tell apart, and the case for the kernel is a cost case.

A budget below is a run's sweep count times its cost in sweep-equivalents, so
the control, at 1000 sweeps a chain and a cost of 1.00, is 1000 units, and its
rate per thousand units is numerically its own summed minimum ESS. The 14.8 and
36.3 here are section 5.1's 15 and 36 unrounded.

| arm | budget, cut-scan units | summed minimum ESS | ESS per thousand units |
|---|---|---|---|
| the control | 1000 | 14.8 | 14.8 |
| cut-only at `d` = 0.16 | 1040 | 26.1 | 25.1 |
| the full draw | 2210 | 36.3 | 16.4 |
| the equal-cost arm | 2210 | 22.9 | 10.4 |

The adopted kernel leads that ranking, and the ranking is weaker than it looks:
it never ran at equal cost, and the statistic grows sub-linearly in chain
length, so a per-unit ratio across four budgets favours the smallest by
construction. The record bridges the two kernels by extrapolation instead, off
a scaling measured on the shipped kernel and not this one, and its own reading
is that cut-only at 0.32 does not match the full draw at 0.16; the choice is
not closed by dominance.

**What the kernel has not earned.** Under the accept rule as written, nothing.
The rule wants a pathology win and treats the core as a gate, and every gain
recorded for the rule draw is on C1, a core cell. The confounded step function
(P2) ran as a must-not-degrade control, and the move cannot reach the chains
that get stuck there: with one tree the root is a nog node only while the tree
has a single split, and every recorded stuck chain sits deeper. The low-noise
cell carries no rule-draw arm; nor do the checkerboard or the diagonal shelf.
And the core gate is checkable on one cell out of four.

**The departure the case rests on.** The kernel's kill criterion was registered
with departures from the program's own, two of which matter: the statistic is
minimum effective sample size rather than coverage, and the cell is C1 rather
than the low-noise cell, on which three shipped mixtures are indistinguishable.
Both follow section 2's ruling, and adoption accepts it: a core cell's mixing
statistic becomes the target rather than a gate, and the pathologies stay
untouched.

**What adoption owes, all of it after the release.** First, the design
amendment that settles what a user sees, the restricted kernel being a private
compile-time build today and not a mode anything can select; the choice is
between deleting the variable axis from the rule draw and adding a sixth name
to `proposal.probs`, which cost twenty-four files the last time a name was
added. Second, the harm controls listed in the kernel's own design and never
run for the restricted draw. Third, the plateau-error gate, a gap and not a
plan: the house rule wants a per-cell check that posterior-mean prediction
error has not worsened once the sampler has reached its plateau, in a
noise-heavy or a large-n cell, and no built cell is either. Fourth, the default
flip itself, since any nonzero share moves every draw from the first sweep of
the first chain: the three equivalence baselines and everything derived from
them are regenerated, with the consumer packages on their lockstep branches
re-recording alongside. That is mechanical work against a bitwise oracle, not a
risk.

### 6.1 A standing rule

No exploration kernel is removed from the engine before this report has been
read. A kill leaves the kernel in place at weight zero rather than deleting it.

## 7. What was not measured

- **Wall time, and anything per second.** No host was quiet; the one-minute load
  ran from 4 to 269 across the arms and 100 to 122 on the machine that ran the
  equal-cost arm, so every cost figure here is a scan count.
- **The plateau-error gate.** No built cell is noise-heavy or large-n, so the
  house rule's harm check has nowhere to run.
- **The sampler's own speed since the default last moved.** A bench-sampler
  comparison on a quiet machine, owed since swap's share moved.
- **The adopted kernel's own controls.** No low-noise control reading, no
  confounded-step arm, no sham arm, no fresh seed block on the ungated function
  and no equal-cost arm, all of which the full draw's run carried. Its cut-scan
  figures come from a census cell at its own data and seed, not a replay of the
  arms' seeds, and that cell takes its 0.16 out of change and birth/death
  together where the benefit arms take it out of change alone; neither
  difference was corrected for.
- **Any rule-draw arm on a pathology.** P1, P5 and P6 carry none.
- **A swap share in five exact single-tree tests.** These check the sampler's
  draws against a brute-force enumeration of every reachable one-tree state.
  Each one that lets its caller set the mixture now sets a positive swap share;
  four two-forest causal tests and one monotone test cannot, so those five never
  exercise swap, which is the one gap in the evidence that swap belongs in the
  engine.
- **Response-swap recovery inside the battery.** Measured once on its own grid,
  where it decided swap's default; never re-run against a new kernel.
- **Whether composing BART with a parametric block helps tree-space mixing.**
  Unmeasured in the literature and here.
- **Seven cells.** The real-data convergence ladder (P3), the
  inhomogeneous-smoothness cell (P4), the no-overlap extrapolation cell (P7),
  the hierarchical variance funnel (P8), the mixed-type regime switch (C2), the
  real-covariate causal cell (C3) and the embedded moving-response cell (C4).
  The funnel and the embedded cell each need an outer sampler from another
  repository.
- **Non-gaussian families.** Every built cell is gaussian or causal gaussian,
  and coverage of a true mean function calibrates a sampler, it does not prove
  it correct.

## Appendix A. Kernels that earned no weight

Section 1 gives each verdict in a sentence; this appendix gives the evidence.

### A.1 The level-fibre step

With the forest's structures held fixed this step is the largest effect in the
program; alongside the tree moves it is nothing. Because it leaves the fitted
function exactly unchanged, its conditional is the leaf prior alone and the
draw is closed form, at about a five-thousandth of a sweep. With C1's
structures frozen, the minimum effective sample size over that cell's 25 points
rises by a paired median of +189.4 at one of the two sweeps where the freeze
was taken and +231.6 at the other, over ten matched pairs, all ten positive.
Live, the same statistic reads -0.9 on the first seed block and -2.5 on the
fresh one, inside the sham arm's own reading, and stacked on the nog-node rule
draw it adds nothing.

Two adverse readings were taken and neither stood: a Single index loss of -5.3
on the first block reads +1.7 on the fresh one, and a wall-time ratio of 1.086,
taken on a loaded host, fell to 0.986 when the two arms were re-measured alone
on a quiet one. Killed as a general default, the explanation being the frozen
run's scope: it measures one channel with the structure held fixed, where the
tree moves move the level fibre faster than the exact draw pays for itself.
What was kept is that one regime, through the automatic `levelGibbs` default of
section 3; a non-frozen forest draws nothing, so the shipped engine is bitwise
unchanged.

### A.2 Perturb, the same-variable cut move

Change always redraws the split variable, so a pure cut displacement happens
only when the redraw lands back on the variable already there, and the census
made the case for building one: change's acceptance rate is far below a
one-position displacement's. The move landed at weight zero, bitwise neutral,
behind a detailed-balance script whose two poisons both fail as designed. Its
benefit run on C1 at `d` = 0.16 did not move the gated statistic, +0.1 +/- 8.1
on the first block and -1.9 +/- 6.5 on the fresh one, with held-out error at
1.032 against the 1.02 margin. Killed at that setting. On the ungated function
it gains +9.5 summed on the fresh block, at t = 5.35, carrying its own 1.027
error regression; that gain is unclaimed, the accept rule not covering that
mean function.

### A.3 Swap

Swap exchanges the rules of a parent and one child, and its evidence is not a
battery cell. A census priced it as mostly wasted work, and on the one
criterion where shipped mixtures do separate, how fast the sampler re-adapts
after the response is swapped under the trees, the mixture with swap at zero
matched the default on every contrast, and only the arm that dropped change as
well lost ground. Swap's 0.1 moved to birth and death and the move was deleted,
then restored at a default of zero on an exact-posterior test: one tree, two
live columns, against a brute-force enumeration of all 62 reachable trees,
where the largest absolute gap in the tree probabilities is 0.0120 without swap
and 0.0008 with it, against a tolerance of 0.004. Swap alone rotates a child's
rule up the tree, which is what that test sees; at fifty and two hundred trees
the ensemble averages the effect away, so no default moved.

### A.4 Cross-chain exchange

Killed on acceptance. An exchange of one tree between two chains at one
temperature is ordinary Metropolis on the product target, the two tree priors
cancelling, so its acceptance rate is closed form; evaluated on states from a
running sampler, without the move ever being proposed, it is 10.1 to 11.7
percent on C1 and effectively zero on the low-noise cell, and alive only on the
one- and two-split trees birth and death already reach. The exchanges that do
accept are, moreover, exactly the ones that would make the chains agree, and it
is the chains' disagreement that pooling turns into interval width.

## Appendix B. Refuted and unbuilt proposals

Two rounds of proposal generation ran under a common bar: a written
Metropolis-Hastings correction, a price in cut scans, a named deficit and a
falsifier runnable in a day. Tree-space geodesics, a learned rotation and a
per-tree temperature each fail on their own terms, and a per-leaf per-variable
histogram cache was refuted on size, the saving being a depth-fold rather than
the leaf-fold claimed. The census refuted two more: an informed death proposal,
choosing which leaf pair to prune by weight rather than uniformly, whose
weights turn out to be effectively a point mass, so the proposal is the uniform
one; and a lifted cut displacement, which would have given a cut move a
persistent direction, where accepted displacements reverse rather than
continue.

Two constructions remain unbuilt with their arguments intact. A
pairwise-collapsed split transfer is the only candidate addressing
representation multimodality directly, and its probe was not built because
pricing a transfer needs the residual net of the other trees and the partner's
leaf statistics, which no existing move sees. A lifted birth and death is cheap
and valid but recomputes at a gain of only 1.06 to 1.17.

The battery's design law comes from an earlier study. Grow-from-root as a
default was killed: every aggregate test passed, but per-cell plateau
posterior-mean error costs of +11.10 percent in a noise-heavy small-n cell and
of +4.47 and +10.66 percent in a large-n one, each past its frozen margin and
each confirmed on fresh seeds, were averaged away by pooling. What the battery
inherits is that law: per-cell checks, thresholds frozen before the run,
mandatory fresh-seed re-runs and a null control that voids the family.

## Appendix C. Source records

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
