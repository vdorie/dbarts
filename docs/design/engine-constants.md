# Engine constants

Status: MEASURED (2026-09-10) on arm64 macOS (10 cores, 32 GB), one
measurement at a time on an otherwise quiet machine, against the shipped
build at 650510f6. Every row below comes from a script under benchmarks/R
named after its constant; each is re-runnable, prints a table, and has no
baseline and no pass/fail exit.

Nine fixed constants shape what the engine will represent and where it
switches strategies. None of them was measured when it was chosen. This note
records, per constant, its value, where the value came from, what a
measurement on this machine says, and one verdict: whether the value BINDS -
whether a workload anyone would actually run reaches it and would be better
off with it moved. It stops there. Which of them become settings is the next
step's question, not this note's.

Two caveats apply throughout. The measurements are single-machine and
single-ISA; a cutoff calibrated here is calibrated for this class of machine,
which is why the ones that bind want a setting rather than a new hard-coded
number. And two of the nine are compile-time constants with no knob at all,
so their arms came from scratch builds of an archived copy of the tree, one
per value, never from the worktree.

## Categorical exhaustive cap

Value 10, at [`categoricalExhaustiveCap`](../../src/bartcore/scan.hpp).
Origin, from the comment at its definition: 2^(P-1) - 1 = 511 candidates
"already cover every factor an R user constructs", and raising the cap
"relocates the boundary rather than removing it". Above the cap the scan
emits the P - 1 sorted prefixes instead of the full balanced-partition set,
in [`categoricalNumEmitted`](../../src/bartcore/scan.hpp) and the two
branches of the scan below it.

["Present-category count above which the exact enumeration gives way"](../../src/bartcore/scan.hpp)
is the definition; the sweep is
["Measures categoricalExhaustiveCap"](../../benchmarks/R/constant-categorical-cap.R),
n = 2000, 200 trees, 300 iterations, one factor and two continuous columns.

| present levels | candidates | msec/iter | structure acceptance | share of splits on the factor |
|---|---|---|---|---|
| 8 (exact) | 127 | 0.960 | 0.2173 | 0.4143 |
| 10 (exact) | 511 | 0.980 | 0.2980 | 0.3271 |
| 12 (prefix) | 11 | 0.980 | 0.2679 | 0.3838 |
| 14 (prefix) | 13 | 0.987 | 0.2716 | 0.3904 |

Enumerating 511 candidates costs nothing measurable over emitting 11
prefixes: the per-proposal cost is dominated by the O(node members) histogram
that precedes the enumeration, not by the candidate loop. A factor-only
control design (one column, 200 trees, same grid) is flat the same way and if
anything rises with P: 0.642, 0.648, 0.666, 0.680 msec/iter at 8, 10, 12, 14.
Acceptance is highest at the cap and drops when the prefix family takes over,
but not by much, and the level count is confounded with the path.

Verdict: does NOT bind on time - the enumeration is free at the cap - and the
acceptance evidence for raising it is weak; a cap arm at P = 12 or 14 exact
needs a scratch build and was not run.

## Linear leaf covariate cap

Value 8, at
[`LinearGaussianLeaf::maxNumCovariates`](../../src/bartcore/model.hpp).
Origin, from the comment at its definition: the per-node sufficient-statistic
scratch is a fixed-size stack array sized for it, and the factory rejects any
designation above it. The refusal is user-visible in
[`leafCovariateDesignationIsValid`](../../src/bartcore/facade.hpp), which is why the sweep can
only measure up to the cap.

The sweep is
["Measures LinearGaussianLeaf::maxNumCovariates"](../../benchmarks/R/constant-linear-leaf-covariates.R),
n = 4000, 25 trees, 300 iterations, 16 available columns.

| designated columns | msec/iter | leaves with n <= q+1 | median leaf n / (q+1) | outcome |
|---|---|---|---|---|
| 4 | 3.687 | 0.0000 | 800.00 | fit |
| 8 | 3.890 | 0.0017 | 444.44 | fit |
| 12 | - | - | - | REFUSED: at most 8 leaf covariates are supported |
| 16 | - | - | - | REFUSED: at most 8 leaf covariates are supported |

Doubling the designation from four to eight costs 5 percent of a sweep, so
the O(q^3) leaf draw is not what the cap is protecting. Conditioning is not
what it is protecting either: at eight columns the median leaf still holds
444 times q + 1 members and only 0.17 percent of leaves are rank-deficient
without the prior's ridge, so a wider leaf regression would still be
identified on this design.

Verdict: BINDS as a refusal - nine columns is a stop, not a slowdown - but
nothing measured says eight is the right place to stop; neither time nor
conditioning is near a limit there.

## Perturb width

Value 1, at [`perturbWidth`](../../src/bartcore/moves.hpp). Origin, from the
comment at its definition: a compile-time constant rather than a knob,
because "acceptance falls off steeply with the displacement and no caller can
set it from evidence", and a width arm therefore needs a private build. That
is what was done: an archived copy of the tree per width, the constant edited
there, installed into a private library, the worktree untouched.

The sweep is
["Measures perturbWidth"](../../benchmarks/R/constant-perturb-width.R), a
perturb-dominant mixture (birth_death 0.10, change 0.10, perturb 0.80), 50
trees, 2000 kept draws after 500 burn, Friedman design.

| width | n | msec/iter | cut acceptance | ESS sigma | ESS per point (median of 20) |
|---|---|---|---|---|---|
| 1 | 500 | 0.099 | 0.5679 | 12.5 | 48.7 |
| 2 | 500 | 0.102 | 0.5414 | 55.6 | 44.1 |
| 4 | 500 | 0.103 | 0.4941 | 25.3 | 40.0 |
| 1 | 2000 | 0.304 | 0.4721 | 29.8 | 23.1 |
| 2 | 2000 | 0.310 | 0.4298 | 31.2 | 19.0 |
| 4 | 2000 | 0.308 | 0.3876 | 13.8 | 30.5 |

Acceptance falls monotonically with the window at both sizes, exactly as the
comment predicts. Nothing gains: the per-point ESS falls with width at
n = 500 and is unordered at n = 2000, and the sigma ESS is unordered at both.
Two chains would be needed to separate the sigma column from noise.

Verdict: does NOT bind. The measurement reproduces the comment's premise and
finds no mixing the wider window buys back.

## Test-fit parallel cutoff

Value 65536, at [`testFitParallelCutoff`](../../src/bartcore/chain.hpp).
Origin: none recorded beyond the comment naming it a cutoff below which the
pool is never created. It is the count of test rows below which a chain
routes its test matrix on its own thread rather than borrowing its share of
the thread budget.

The sweep is
["Measures testFitParallelCutoff"](../../benchmarks/R/constant-testfit-parallel-cutoff.R),
n.train = 2000, 75 trees, one chain, 100 iterations per round, 5 rounds,
minimum taken, at one thread and at four.

| n.test | serial msec/iter | threaded msec/iter | speedup | path taken |
|---|---|---|---|---|
| 8192 | 4.630 | 4.600 | 1.007 | serial (both) |
| 16384 | 8.680 | 8.650 | 1.003 | serial (both) |
| 32768 | 16.500 | 16.500 | 1.000 | serial (both) |
| 65536 | 31.000 | 20.080 | 1.544 | threaded |
| 131072 | 64.260 | 23.630 | 2.719 | threaded |
| 262144 | 132.070 | 44.430 | 2.973 | threaded |

Below the cutoff the two columns are the same code and agree to a fraction of
a percent, which is the check that the sweep measures the routing and nothing
else. At the cutoff itself threading already wins 1.54x, and the serial cost
is linear in n.test: at 32768 a caller pays 16.5 msec per iteration on one
thread where the threaded path, extrapolating its own 131072-to-262144 slope
back, would cost roughly 8. The cutoff is above the point where fanning out
pays.

Verdict: BINDS. The cutoff sits past the crossover; the serial path is
charged to test sets that would already profit from the pool.

## Predict parallel cutoff

Value 1e7 traversals, at
[`predictParallelCutoff`](../../src/bartcore/sampler.hpp). Origin, from the
comment at its definition: an arithmetic estimate, "the measured cost of one
is ~2.6 ns, so 1e7 of them is ~26 ms of work", chosen to sit comfortably
above spawn and join cost. dec-B93 asks for it calibrated rather than left
fixed.

The sweep is
["Calibrates predictParallelCutoff"](../../benchmarks/R/constant-predict-parallel-cutoff.R).
Both arms are the same build: the test seam
[`bartcore_setPredictParallelCutoff`](../../src/R_interface_bartcore.cpp)
replaces the constant, so cutoff 1 forces the fan-out and a cutoff above the
whole sweep forces the inline path. Four threads, 7 rounds, minimum taken.

| trees | draws | rows | traversals | serial ms | threaded ms | speedup |
|---|---|---|---|---|---|---|
| 10 | 10 | 10 | 1.0e3 | 0.059 | 0.130 | 0.455 |
| 10 | 10 | 30 | 3.0e3 | 0.063 | 0.127 | 0.495 |
| 10 | 10 | 100 | 1.0e4 | 0.075 | 0.129 | 0.582 |
| 10 | 10 | 300 | 3.0e4 | 0.119 | 0.154 | 0.772 |
| 10 | 10 | 1000 | 1.0e5 | 0.470 | 0.311 | 1.512 |
| 10 | 10 | 3000 | 3.0e5 | 2.090 | 0.856 | 2.442 |
| 10 | 10 | 10000 | 1.0e6 | 7.352 | 2.773 | 2.651 |
| 75 | 200 | 10 | 1.5e5 | 0.537 | 0.275 | 1.951 |
| 75 | 200 | 100 | 1.5e6 | 2.608 | 0.864 | 3.018 |
| 75 | 200 | 1000 | 1.5e7 | 52.749 | 15.460 | 3.412 |
| 75 | 200 | 10000 | 1.5e8 | 638.001 | 179.229 | 3.560 |
| 75 | 200 | 100000 | 1.5e9 | 6133.886 | 1727.984 | 3.550 |

The crossover is between 3e4 and 1e5 traversals; the fan-out is already at
1.5x by 1e5 and saturates near 3.5x on four threads. The spawn and join cost
the sweep isolates is about 60 to 70 microseconds, two orders of magnitude
below the 26 msec the constant budgets for. The two arms agree where they
overlap (1.5e5 against 1e5 and 3e5), so the crossover is a property of the
traversal count and not of the fit shape.

Verdict: BINDS, and by about 200x. A calibrated value on this machine is
roughly 5e4 traversals; at the shipped 1e7 every replay up to 1e7
traversals - 46 msec of avoidable serial work at the top of that range - runs
inline on one thread with the fan-out available.

## Sparse density threshold

Value 0.2, at [`sparseDensityThreshold`](../../src/bartcore/data.hpp).
Origin: none recorded; the comment states the rule ("CSC-built columns at or
below this nonzero fraction take rank-bitmap storage") without attribution,
and the sparse-columns design note states the 20 percent number the same way.
The two layouts are [`SparseColumnData`](../../src/bartcore/data.hpp) - a
bitmap, a word-rank index and the packed nonzero codes - and a dense
`xint_t` per row.

The sweep is
["Measures sparseDensityThreshold"](../../benchmarks/R/constant-sparse-density-threshold.R),
n = 1e5, 100 columns, 50 trees, each density in a fresh process, exactly
round(n * density) distinct rows per column so the realized density is the
requested one.

| density | storage | rss growth MB | dense model MB | sparse model MB | msec/iter |
|---|---|---|---|---|---|
| 0.050 | sparse | 48.1 | 20.0 | 2.9 | 9.277 |
| 0.100 | sparse | 48.7 | 20.0 | 3.9 | 10.991 |
| 0.150 | sparse | 49.7 | 20.0 | 4.9 | 12.580 |
| 0.190 | sparse | 46.7 | 20.0 | 5.7 | 13.336 |
| 0.210 | dense | 63.7 | 20.0 | 6.1 | 9.640 |
| 0.250 | dense | 44.6 | 20.0 | 6.9 | 9.654 |
| 0.350 | dense | 45.3 | 20.0 | 8.9 | 9.739 |
| 0.500 | dense | 25.6 | 20.0 | 11.9 | 9.963 |

The time column is the measurement. It steps cleanly across the threshold -
13.34 msec per sweep at 0.19 against 9.64 at 0.21 - and the dense side is
flat at about 9.7 regardless of density, so the whole 38 percent is the rank
decode. Sparse gather cost rises with density inside the sparse regime and
meets the dense cost only well below 0.2. The resident-set column does NOT
isolate the store: everything else a sampler allocates costs the same in
either layout and swamps a 14 MB difference, so the memory side of the
tradeoff rests on the model columns, which say sparse still saves 3.5x at
0.19.

Verdict: BINDS as a real tradeoff at its stated value - 38 percent of sweep
time bought for 3.5x of predictor memory - and the threshold is a memory
choice, which is what the sparse-columns note claims. Whether 0.2 is the
right price is a per-workload question a fixed number cannot answer.

## GP leaf max leaf size

Value 256, at [`maxLeafSize_`](../../src/bartcore/model.hpp), set from
[`gp`](../../R/model.R)'s max.leaf.size argument, so unlike the other engine
constants this one is already settable. Origin: the four call sites in
[`GPGaussianLeaf`](../../src/bartcore/model.hpp) abandon the Gaussian process
and score a leaf as a constant leaf above it, the GP draw being cubic in the
leaf's member count. dec-B110 asks how much of a fit that silently disables.

The sweep is
["Measures maxLeafSize_"](../../benchmarks/R/constant-gp-max-leaf-size.R),
n = 1200, 20 trees, 60 iterations, two designated columns. The grid is small
because the 1024 arm is cubic in a leaf that can hold most of the design.

| max.leaf.size | msec/iter | share of leaves fallen back | share of ROWS in a fallen-back leaf | rmse |
|---|---|---|---|---|
| 128 | 2.702 | 0.7076 | 0.9522 | 0.2082 |
| 256 | 21.984 | 0.4447 | 0.8277 | 0.1961 |
| 512 | 276.354 | 0.2558 | 0.5629 | 0.1876 |
| 1024 | 2232.528 | 0.1058 | 0.2840 | 0.1900 |

At the shipped default, 83 percent of training rows sit in a leaf that took
the constant-leaf path: a fit asked for GP leaves is mostly not getting them,
and nothing today says so. Raising the ceiling fixes that and costs 8 to 12x
per doubling, the cubic law plus the shrinking fallback share. The rmse gain
is real but small and flattens by 512.

Verdict: BINDS, in both directions. The default silently disables the model
on most of the data, and moving it is unaffordable by a factor of ten per
step - which is the case for dec-B110's counter and warning: a user needs to
be told which regime the fit landed in.

## Person-period row cap

Value 1e7, the max.rows argument of [`hazard`](../../R/family.R), enforced by
[`expandDiscreteTimeHazard`](../../R/dbarts.R) before it allocates. Origin:
the comment at the expander calls it a guard on an over-fine time grid, and
the refusal names the two levers - coarsen with hazard(breaks =), or raise
the cap with hazard(max.rows =).

The sweep is
["Measures the person-period row cap"](../../benchmarks/R/constant-person-period-rows.R),
10 predictor columns, 20 periods, each row count in a fresh process.

| target rows | subjects | expanded rows | seconds | peak resident MB |
|---|---|---|---|---|
| 1e6 | 95239 | 998506 | 0.074 | 191.6 |
| 1e7 | 952381 | 10004962 | 0.571 | 1896.7 |
| 3e7 | 2857143 | 29989717 | 1.624 | 5692.5 |

The expansion is linear and cheap in time: even three times over the cap it
is under two seconds. Memory is the whole cost, about 190 MB per million
rows at ten columns, and it scales with the column count. At the cap a fit
holds 1.9 GB of expanded design before the sampler has allocated anything of
its own; the 3e7 arm, which only runs with the cap raised, holds 5.7 GB.

Verdict: BINDS at a defensible place. 1e7 rows is where the expansion stops
being something a 16 GB machine absorbs, and the cap is already settable per
fit, which is the right shape for a number that depends on the column count
and the host.

## The xint caps

Values 65533 cuts at
[`maxNumCutsRepresentable`](../../src/bartcore/data.hpp), 65535 categorical
levels at [`maxCategories`](../../src/bartcore/data.hpp), and 65534 for an
ordered factor via [`maxLevelsForKind`](../../src/bartcore/data.hpp). Origin,
from the comment over [`xint_t`](../../src/bartcore/data.hpp): all three fall
out of the 16-bit predictor code, which reserves
[`naCode`](../../src/bartcore/data.hpp) for missing, and an ordered factor
spends one more code because its K - 1 midpoint grid must itself fit.
Widening the code is out of scope for this arc.

This one is a probe, not a timing sweep:
["Probes the xint_t caps"](../../benchmarks/R/constant-xint-caps.R).

| channel | request | outcome |
|---|---|---|
| control n.cuts | 65532, 65533, 65534, 100000 | accepted, every one |
| sampler setCutPoints length | 65533, 65534 | accepted, both |
| categorical predictor levels | 65534, 65535 | accepted |
| categorical predictor levels | 65536 | REFUSED: factor 'g' has more than 65535 levels |
| ordered factor levels | 65533, 65534 | accepted |
| ordered factor levels | 65535 | REFUSED: ordered factor 'g' has more than 65534 levels |

The two caps behave differently at the boundary. A level count past its cap
is a named refusal at the bridge, which is what a caller wants. A cut request
past its cap is accepted and silently clamped where the store sizes its grid,
so a caller who asks for 100000 cuts gets 65533 and is never told.

Whether a real design reaches one: the per-column cut count is min(n.cuts,
distinct values - 1) and n.cuts defaults to 100, so reaching 65533 takes an
explicit request on a design with more rows than that. Asking is not
expensive - at n = 70000 with five trees, a 65533-cut grid runs at 0.684
msec per iteration against 0.722 for the default 100 - so the cap is reached
only by a caller who went looking for it. A 65535-level factor is a stranger
object still: it needs at least that many rows to be non-degenerate.

Verdict: does NOT bind. No design that arrives on its own gets near either
cap. The one defect the probe finds is the asymmetry: the cut cap clamps
silently where the level caps refuse by name.
