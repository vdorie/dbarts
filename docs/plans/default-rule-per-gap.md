# default-rule-per-gap: the default rule gives a few-valued column one cut point per gap

Status: PLANNED 2026-10-09 (dec-B313). Every claim below is marked (ran), run on the reference build of
bartcore 8ad1fa58 in a private library with the probes under scratch/drplan, or (read), read in the tree or
the source named.

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for every fit under the default rule (`useQuantiles = FALSE`) with a numeric column
whose distinct values fall within the boundary of Open call 1, at creation, `setData`, a refresh and in
xbart's folds: the tree prior over that column's thresholds changes (Context). NEUTRAL, bit for bit, for
every fit under the quantile rule (but the columns of Open call 3), every default fit whose numeric columns
all hold more distinct values than the boundary, constant and all-missing columns (one point already),
factor columns of either kind, grids set by `setCutPoints`, and states, copies and reloads, which install
the grid they carry.
window: before the merge to main (dec-B313, TODO default-rule-per-gap). An engine slice of its own after
repeated-cut-restore, which has landed; engine slices stay serial.
budget: ~550 lines (data.hpp ~60, tests/cpp ~170, tinytest ~130, bd-balance arm and workflow ~50, help ~40,
NEWS ~10, design docs ~70, MANIFEST ~15), plus machine-written snapshot values if Open call 1 goes as
recommended.

## Goal

Under the default rule a numeric column whose distinct values leave no more gaps than it can be given points
gets one cut point in each gap between neighbouring values, the grid the quantile rule gives it; any other
column keeps its `n.cuts` evenly spaced points. A 0/1 column holds one point at 0.5 where it holds 100 today.
The change is in the one routine that derives a numeric grid, so creation, `setData`, a refresh and xbart's
folds follow it together. No state format, C API, facade or R surface change. The tier is "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Context

The ruling, dec-B313 (VD 2026-10-07). Shown the plan's reading of dec-B297, which kept the evenly spaced
points, against one point per gap, VD asked: "What do other packages do? It seems inefficient from a mixing
perspective to have so many possible proposals that don't change anything." and ruled: "Use option B, One
point per gap when a column has fewer distinct values than n.cuts; evenly spaced otherwise." The entry's
cost, the orchestrator's text: every default fit with a binary or few-valued numeric column draws
differently from 0.9-34, and a gap's chance of being split no longer follows its width. It revises dec-A178's
second call. Related rulings, none changed here: dec-B297 (a grid holds each value once), dec-B298 (a refresh
derives "Up to n.cuts, from the new values"), dec-B299 (a refresh onto fewer distinct values shrinks), dec-B300
(no grid repeats a point; weights later as probabilities), dec-B231 (a grid is data, kept until a call says
otherwise), dec-B235 (the quantile grid spread, VD: "Fix it everywhere, bartBT included."), dec-B311 and
dec-B312 (where splits go when a grid is replaced), dec-B172 (the uniform range is over finite values).

Where the grid comes from (read). [`ColumnStore::deriveNumericCuts`](../../src/bartcore/data.hpp) is the one
derivation; creation ([`ColumnStore::buildCutsForColumn`](../../src/bartcore/data.hpp)), the whole-data
replacement ([`ColumnStore::setData`](../../src/bartcore/data.hpp), which rebuilds every numeric grid), a
refresh ([`ColumnStore::refreshCutsForColumn`](../../src/bartcore/data.hpp),
[`ColumnStore::refreshCutsForCscColumn`](../../src/bartcore/data.hpp), reached from `setPredictor` and the
per-observation forms with `updateCutPoints = "position"` or `"value"`) and the xbart data handle
(`bartcoreDataHandle`, which builds a store; its row-subset views copy the parent's grid) all call it. Its
uniform branch is [`ColumnStore::fillCutsUniformly`](../../src/bartcore/data.hpp) and
[`ColumnStore::fillCutsUniformlyCsc`](../../src/bartcore/data.hpp) into
[`ColumnStore::fillCutsOverRange`](../../src/bartcore/data.hpp); its quantile branch collects the sorted
distinct finite values ([`ColumnStore::quantileGridForColumn`](../../src/bartcore/data.hpp),
[`ColumnStore::quantileGridForEntries`](../../src/bartcore/data.hpp),
[`ColumnStore::finishQuantileGrid`](../../src/bartcore/data.hpp)) and fills
([`ColumnStore::fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp)), taking every midpoint when the
midpoints number no more than `n.cuts`.

Paths that do not derive (read). `setCutPoints` installs the caller's grid
([`ColumnStore::setCutPointsForColumn`](../../src/bartcore/data.hpp)). `setState`, a copy and a reload
install the grid the state carries ([`setCutPointsForColumn`](../../src/bartcore/sampler.hpp), from the state's
`cutPoints`), so a state or saved fit from an earlier build continues on its 100-point grid until a call
re-derives. `predict` codes new values against the held grid
([`ColumnStore::codeFor`](../../src/bartcore/data.hpp)). A warm start maps the donor's splits onto the
recipient's own grid by value ([`Tree::mapOldCutPointsOntoNew`](../../src/bartcore/tree.hpp)). The flat C
API has no entry that builds or takes a grid ([`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h));
stan4bart builds its samplers through `dbartsData` and `new("dbartsSampler")` and calls no grid setter
(read, stan4bart branch bartcore 963956b), so its fits follow creation.

What counts as a distinct value (read): the quantile collectors count finite values only (NA, NaN and
`Inf` drop out), add `+0.0` so -0 and 0 are one value, count a CSC column's implicit zeros as one 0, and read
every training row of the store, so weights, an active-row mask, offsets and test rows do not enter. The
default rule uses the same set; the grid is derived data (dec-B231), not a function of the likelihood.

Why it changes the posterior (read, ran). A node may split on a column while its interval holds a point
([`Tree::variableAvailable`](../../src/bartcore/tree.hpp)). After a split on a 0/1 column at one of its 100
points both children still hold points, so they count the column as available, it takes a share of their
variable draws, and a proposal on it empties a side. With one point the children hold none. In one tree on
a 0/1 column, noise response, 40,000 draws (scratch/drplan/prior-mass.R), the share of draws holding a split
is 0.177 at one point, 0.151 at two, 0.138 at five and 0.131 at 100 (ran). On a column of uneven values the
prior over gaps also changes from proportional to width to uniform.

How far a fit moves (ran, scratch/drplan/move.R): n = 300, 200 trees, four 0/1 columns, two of five values,
four uniform, the per-gap grid installed by `setCutPoints` on a fresh sampler, ten seeds: test RMSE 0.440
today and 0.441 per gap (paired difference 0.001, se 0.002), sigma 0.872 and 0.870, effective sample size of
sigma and of a test fit not distinguishable (differences 28 +- 35 and -29 +- 27), splits on the 0/1 columns
112.5 and 116.5 per draw (+4.0, se 0.3), sweep time equal. A continuous column on fewer rows than `n.cuts`
is also within the rule; there it costs a little test error (Open call 2).

0.9-34 (ran, scratch/libs/cran): the default rule gives 100 points to a 0/1 column, a 13-value column and a
continuous one; the quantile rule 1, 12 and 100. The branch today gives the same (ran). Other packages: the
BART package's `bartModelMatrix` takes the midpoints when a column has fewer distinct values than `numcut`,
under either rule unless the column is flagged continuous, and `numcut` points otherwise (ran on a 0/1, a
13-value, a six-value and a continuous column: 1, 12, 5 and 100; source read). BayesTree spaces its points
evenly whatever the values, bartMachine draws a threshold from the distinct values in the node, and bcf takes
midpoints under 8 distinct values (as recorded in dec-B313; not re-read).

## Constraints

- No state format, C API ([dbarts.h](../../inst/include/dbarts/dbarts.h)), facade virtual, bridge or R surface
  change. `DBARTS_C_API_HASH` does not move.
- States, copies and reloads keep the grid they carry; nothing converts an old grid.
- One derivation for every path: the change sits in `deriveNumericCuts` and what it calls, never in a caller.
- No sort of a column with more distinct values than the boundary: today's uniform path reads a column
  once, and the new count must stop once it passes the boundary (Speed).
- Out of scope: factor columns (their grids follow the level table), `setCutPoints`, which midpoints the
  quantile rule's thinning picks (only the point it places in a gap gains the guard), cut-point weights
  (dec-B300).

## The rule

Under the default rule, for numeric column j with U distinct finite values and `requestedNumCuts[j]` = n:

- U of 2 or more within the boundary: the U - 1 gap points, point k in the gap between sorted values u_k and
  u_(k+1). The boundary is Open call 1: U - 1 <= n (recommended, the quantile rule's) or U < n (the words).
- Otherwise: n evenly spaced points over the finite range, as today. U of 0 or 1 is today's one point.

A gap's point is its midpoint, 0.5 (u_k + u_(k+1)), the quantile rule's expression, so a few-valued column's
grid is the quantile rule's to the bit. Where the midpoint is not strictly below u_(k+1), which happens only
when the two are adjacent doubles or their sum overflows, the point is u_k, which a row at u_k does not
exceed and one at u_(k+1) does (Open call 3). With that, every gap point separates its pair and the points
strictly increase. Elsewhere the guard changes no double: the midpoint of two doubles that are not adjacent
rounds strictly between them (read; checked on adjacent doubles and on an overflowing pair, ran).

## Change

All in [`ColumnStore`](../../src/bartcore/data.hpp).

1. A bounded distinct collector, dense and CSC, returning the sorted distinct finite values (with `+0.0`
   and the CSC implicit zero, as the quantile collectors) when there are at most n + 1 of them, and stopping
   with "more" as soon as there are n + 2 (n - 1 and n under the words). Expected linear time with no sort before the stop; scratch kept
   in the store, sized to n + 2, never to the row count. The implementer may share it with
   `quantileGridForColumn` only if the quantile path's doubles do not move.
2. `deriveNumericCuts`'s uniform branch: run the collector; within the boundary, set `numCuts[j]` to U - 1
   and fill through the quantile fill; otherwise the uniform fill as today. A `QuantileGrid` built from the
   collected values with `inducedNumCuts = U - 1` reaches `fillCutsFromQuantileGrid` unchanged.
3. The gap point: one helper used by `fillCutsFromQuantileGrid` for every point it places, thinning
   included: m = 0.5 (a + b); m if a <= m and m < b, otherwise a. `dropRepeatedCuts` stays as a backstop.
4. Comments: `deriveNumericCuts` (the two rules and the boundary), `fillCutsOverRange` ("the uniform rule's
   count and fill" now applies past the boundary only), `fillCutsFromQuantileGrid` (the guard; "a cut is a
   midpoint" gains its exception), the class's header paragraph on per-column counts.

Nothing else changes: the refresh, `setData` and the handle reach the rule through `deriveNumericCuts`, and
the splits on a refreshed column move by position or value as dec-B311 rules
([`Tree::rescaledSplitIndex`](../../src/bartcore/tree.hpp)).

A column whose values grow later (read): with `updateCutPoints = "none"` the grid is kept, as today, so a
column created 0/1 and later given many values keeps its one point; with `"position"` or `"value"` the grid
is derived again for the new values, n evenly spaced points once they pass the boundary, and back to gap
points if they fall within it. Today the same calls keep or re-derive the same way; what differs is the grid
a few-valued column starts with. setPredictor's help says so (Docs).

## Tests

tests/cpp ([test_data.cpp](../../tests/cpp/test_data.cpp)):
- [`testDistinctCutGrids`](../../tests/cpp/test_data.cpp): the uniform counts become {4, 1, 1, 1, 1, asked}
  and the quantile counts {4, 1, 1, 1, 1, asked} (the five adjacent doubles hold four points under the guard,
  where today they hold five under the uniform rule, the last separating nothing (ran), and three under the
  quantile rule (read, the test's own counts)); the 0/1 grid is {0.5} under both. The guarded points on four
  adjacent doubles are 1, 1 + eps and 1 + 2 eps, bd-balance's narrow cuts exactly (ran).
- New `testDefaultRulePerGap`, at n = 10 on unevenly spaced values: U from n - 1 to n + 2 straddles the
  boundary, gap points exactly where it says and n evenly spaced points otherwise; within it the grid equals
  the quantile rule's bit for bit; dense and CSC builds equal (implicit zeros once, -0 and 0 one value); NaN
  and `Inf` do not count (two finite values beside `Inf` give one point); per-column counts follow each
  column's n; a 0/1 column gives {0.5} at n = 1, 2 and 100.
- The guard: four adjacent doubles under both rules hold three points and code to four distinct codes; a
  pair near the largest double whose sum overflows gets the lower value.
- [`testRefreshDerivesAsCreation`](../../tests/cpp/test_data.cpp) gains a 0/1 column and a five-value column:
  creation, `setData` and a refresh give one grid; a refresh from 0/1 values onto 200 distinct values gives n
  evenly spaced points, and back gives {0.5}.

tinytest:
- [test-cut-grid-distinct.R](../../inst/tinytest/test-cut-grid-distinct.R): the uniform counts
  ["c(100L, 5L, 1L, 1L, 100L, 100L)"](../../inst/tinytest/test-cut-grid-distinct.R) become those above, the
  quantile ones likewise, and the narrow refresh's `c(5L, 100L)` its new count.
- A new section in [test-quantile-grid.R](../../inst/tinytest/test-quantile-grid.R) beside
  ["spreadMidpoints"](../../inst/tinytest/test-quantile-grid.R): through `dbarts`, `bart`, `bartBT`, a
  BayesTree-spelled `bart` call and `rbart_vi`, a 0/1, a six-value and a continuous column give the same
  grid under either rule wherever the column is within the boundary; xbart builds through the same store
  and returns no sampler to read one from.
- A 0/1 column draws bitwise the same at `n.cuts` 1, 2 and 100 (one point at 0.5 each); under the literal
  boundary of Open call 1, `n.cuts = 2` gives two points and drops out of the test.
- `"none"` keeps {0.5} after a 0/1 column is given continuous values; `"position"` and `"value"` re-derive.
- A state carrying a 100-point grid on a 0/1 column (set with `setCutPoints`, then stored) installs that
  grid on `setState`, a copy and a reload, and the restored chain continues it.
- Every other expectation the suite changes is a pinned count or draw on a few-valued column; the
  implementer lists each with its reason.

## Baselines

Counted with a tracer on the entry points that hand a store to the engine (scratch/drplan/hook.R: creation,
re-creation, the xbart handle, `setData` and a refresh), comparing each numeric column's grid today with the
per-gap grid under both boundaries (ran).

- Equivalence, current equivalence-51201107 ([MANIFEST](../../benchmarks/baselines/MANIFEST)), 55 scenarios:
  under either boundary 8 move: factorpartial, hazard (its period column, 6 values), hurdle (95 distinct values
  in every column), leaffactormixed, mixedmatrix, sparse, wideFactorIndicators (0/1 indicator columns) and xbartmixed.
  Under the recommended boundary bart2twoforest moves as well (columns of 100 distinct values), 9 in all. The
  other 46 or 47, quants among them (quantile rule), do not.
- bcf-equivalence-1b7d730c (15 scenarios) and multinomial-equivalence-80b1c8d4 (11): no affected column
  (ran); any mover is a defect.
- The four test-reproducibility-*.R files: under the words none moves (binaryResponse and multithreaded have
  no few-valued column; singleThreaded and xbart fit columns of 100 distinct values). Under the recommended
  boundary singleThreaded and xbart move, regenerated by tools/regenerate-snapshots.R replaying each whole
  file on the reference build; binaryResponse and multithreaded stay.
- Exact gates (exact-gates.yaml's list, quick, ran under the tracer): draws move in backfit-exact,
  heteroscedastic-exact and multinomial-exact (0/1 columns at 100 points) and hazard-reduction (both sides
  derive one grid, so it stays bitwise); bd-balance's narrow arm keeps its grid under the guard and its
  draws. Under the recommended boundary the gates that set n.cuts to one fewer than equally spaced cell
  values (bd-balance and its zeroweight arm, bcf-exact, hazard-exact, linear-exact, mask-redraw-exact,
  monotone-reference, monotone-exact-enumeration, monotone-successive-conditional at 100 values) get the same
  partitions at midpoints in place of evenly spaced values: draws unchanged, split values changed, and
  bd-balance and monotone-exact-enumeration, which match split values to their own evenly spaced cuts,
  restate their cuts as midpoints. bcf-latent-exact was cut off by the probe's limit and not counted; the
  implementer counts it.
- SBC (sbc.R arms run for a few replications under the tracer): the arms on fewer than 100 rows move,
  bcf-probit-weak and bcf-logistic-weak at 40 rows, gp, gp-weighted and gp-mixed at 80, monotone-1-leaf,
  monotone-1-joint and monotone-bd at 20 (ran). None is in sbc.yaml's matrix: of its twelve,
  discrete-selfcheck builds no sampler and the other eleven do not move (ran).
- heteroscedastic-exact and multinomial-exact weigh a two-cell tree on their 0/1 column by `base` alone
  (read), which is the prior only where the children hold no point; under the change it is exact, where today
  it misses a factor of (1 - base / 4)^2 for the children (inferred from prior-mass.R, not run against these
  gates). The reviewer runs both in full mode on the base and the slice and reports the two gaps.
- Re-record: the movers fresh on the reference build (`EQUIVALENCE_SCENARIOS=<the 8 or 9>`,
  `EQUIVALENCE_CORES=2`), merged into a copy of 51201107 in its scenario order, named after the slice's code
  commit, with a MANIFEST row; 51201107 demoted to historical. Partition against 51201107 in z mode: the
  others identical, the movers' |z| reported as sizing, not as the oracle. The new file reproduces 55 of 55
  under `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17): a new bd-balance arm, `fewvalued`, its four cells at uneven values (0, 1, 3, 10)
  with the default `n.cuts = 100`, asserting the grid {0.5, 2, 6.5} first and then the enumeration the other
  arms use, whose three cuts are now one per gap; on the base build the grid holds about 10, 20 and 70 points
  in the three gaps and the children of a split stay available, so the arm fails there, by its grid
  assertion and, with that assertion removed for the run, by its statistic. It joins the
  [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) loop
  ["for arm in narrow constmissing"](../../.github/workflows/exact-gates.yaml). Beside it, one SBC arm that
  moves (monotone-bd, one tree on 20 rows) passing its band.
- `Rscript benchmarks/R/mutation-battery.R verify-anchors` passes, and each mutant whose killer is a moved
  scenario still kills against the new file.

## Speed

The change is at grid derivation, not in a sweep: no bench-sampler.R arm reaches it (its columns are uniform
draws on 1000 rows or more and it refreshes no grid; read). In its place, same machine, shipped builds of
the base and the slice, median of five interleaved creations: `dbarts()` creation at n = 1e6 with 50 uniform
columns (the collector stops early; accept 1.05) and with 50 0/1 columns (the collector reads every row;
report it beside the quantile rule's creation on the same data, which sorts). A refresh loop on one uniform
column at n = 1e5, 200 calls, accept 1.05.

## Docs

- [cut-grid.md](../design/cut-grid.md), the design note the tier requires: [The rule](../design/cut-grid.md#the-rule)
  states the default rule's two cases and the gap point; [What moves](../design/cut-grid.md#what-moves)'s
  0/1 and six-value rows and its closing forward pointer to dec-B313 are restated; a dated section records
  the mixing argument, prior-mass.R, move.R and smalln.R, and the boundary as ruled.
  [Checks](../design/cut-grid.md#checks) names the `fewvalued` arm.
- [quantile-grid.md](../design/quantile-grid.md): [What does not move](../design/quantile-grid.md#what-does-not-move)'s
  "The uniform grid, which is the default." points here, and "A cut is a midpoint, never an observed value"
  gains the adjacent-doubles exception.
- [data-store.md](../design/data-store.md)'s `cutPoints` bullet (it still says a grid can repeat a value) and
  [sparse-columns.md](../design/sparse-columns.md)'s sentence on uniform grids over CSC columns.
- [classic-compare.md](classic-compare.md): one entry under
  [What differs, and why](classic-compare.md#what-differs-and-why); its numbers wait for the release-candidate
  run (manual-gate-runs).
- Help: man/dbartsControl.Rd (`useQuantiles`, `n.cuts`), man/bart.Rd (`n.cuts`, `useQuantiles`),
  man/bartBT.Rd (`usequants`, `numcut`, Decision Rules), man/xbart.Rd (`n.cuts`, `useQuantiles`) and
  man/dbartsSampler-class.Rd (`setCutPoints`'s "with useQuantiles, a column with no more distinct values"
  sentence now holds under either rule; `updateCutPoints = "none"` keeps the grid a column got at creation,
  one point for a 0/1 column). The text: with either setting, a numeric predictor whose distinct values
  leave no more than `n.cuts` gaps (or: has fewer distinct values than `n.cuts`, per Open call 1) gets one
  cut point halfway between each pair of neighbouring values.
- NEWS, under USER-VISIBLE CHANGES beside the item on distinct cut points: under the default
  `useQuantiles = FALSE`, a predictor with few distinct values, a 0/1 column, a small count or any column on
  fewer rows than `n.cuts`, gets one cut point halfway between each pair of neighbouring values, as
  `useQuantiles = TRUE` gives it, where 0.9-34 placed `n.cuts` equally spaced points; a 0/1 column has one
  cut point where it had 100, and fits with such a column change. The existing quantile item's "and the
  default useQuantiles = FALSE, are unaffected" stays true of that item.
- TODO's item goes and the ledger records the calls at landing.

## Steps

1. Change 1 to 4 with the tests/cpp tests; `make && ./test_bartcore` green, and under
   `-fsanitize=address,undefined`.
2. `R CMD INSTALL --preclean -l <lib> .`; the tinytest changes; the full suite green.
3. The `fewvalued` arm and its workflow line; run it on the base build (fails) and the slice (passes).
4. Docs, help and NEWS.
5. After review, in their own commit: the re-record and MANIFEST row, the regenerated snapshot files if
   the boundary moves them, the gates restating their cuts, the oracle runs.

## Gates

On the slice tip against its own library, independently of the implementer, posterior-changing
([RNG classes and their gates](README.md#rng-classes-and-their-gates)):
- tests/cpp plain and sanitized; the R-loaded ASAN path over test-cut-grid-distinct.R and
  test-quantile-grid.R (the new collector reads raw columns).
- Full tinytest suite; `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift,
  doc-freshness, the NEWS parse.
- Reference build: the three equivalence harnesses as in Baselines; the four snapshot files.
- exact-gates.yaml's list in quick mode with the new arm, and full mode for backfit-exact,
  heteroscedastic-exact, multinomial-exact and the `fewvalued` arm.
- Consumers against the slice's library: stan4bart's suite and its posterior baselines, bartCause,
  treatSens and bairrtt. bartCause's response fit and treatSens's carry the 0/1 treatment as a predictor
  column, and bairrtt its 0/1 `z` (read: bartCause responseFit.R, treatSens branch dbarts-1.0
  treatSensBART.R, bairrtt irt_causal_bart.R), so their draws move; bairrtt's latent columns take its own
  grid through `setCutPoints` and do not. A failure that is not a pinned draw is a defect.
- Speed as above.

Reviewer's mutants, each of which must fail a test:
- the boundary off by one either way (`U <= n + 1` read as `U <= n`, or `U < n` as `U <= n`);
- the gap guard dropped (four adjacent doubles hold two points, and bd-balance's narrow arm fails its grid
  assertion);
- the per-gap branch on the dense path only (the CSC twin differs);
- the collector counting non-finite values, or -0 and 0 apart;
- the collector stopping one distinct value early (at n + 1 in place of n + 2 under the recommended
  boundary);
- the rule placed in `buildCutsForColumn` and not in `deriveNumericCuts` (a refresh differs from creation).

## Stop conditions

Stop and report when: the diff passes ~900 lines before machine-written values; a scenario, snapshot or
harness outside the lists above moves; an exact gate fails on two seeds; a creation or refresh arm stays
above 1.05; the change needs a state, facade, bridge or C API change.

## Interactions

- gp-copy-continuation (planned, engine): disjoint code; both re-record the equivalence baseline, on
  disjoint scenarios. Whichever lands second records against the other's file.
- repeated-cut-restore (landed): this slice keeps its every invariant (distinct points, refresh as creation,
  splits by position or value) and changes only which points the default rule places.

## Open calls

1. The boundary. VD's words, "fewer distinct values than n.cuts", are the BART package's test (U < n);
   the quantile rule takes every midpoint while they number no more than n (U - 1 <= n). They differ at U =
   n and U = n + 1. Recommended: the quantile rule's. At U = n the words put n points in n - 1 gaps, so some
   gap holds two points that make one split, the proposals VD named; at U = n + 1 both give n points, but
   evenly spaced ones can leave one gap empty and double another, where one per gap separates every pair;
   and with it a few-valued column has one grid under either rule, and a 0/1 column one point at every
   `n.cuts` (under the words, two at `n.cuts = 2`). Cost: bart2twoforest, two snapshot files and the
   midpoint restatement in about nine exact gates also move, all by machine or mechanically. Evidence that
   would change it: none on the model; VD preferring the BART package's test.
2. Continuous columns on fewer rows than `n.cuts`. The words reach them: a column of 50 uniform draws has
   50 distinct values. Measured (ran, scratch/drplan/smalln.R, Friedman's function, 200 trees, 8 seeds):
   test RMSE at 30 rows 2.944 today and 3.033 per gap (+0.089, se 0.042), at 60 rows 1.956 and 2.004
   (+0.048, se 0.017), at 90 rows 1.805 and 1.815 (+0.011, se 0.024); evenly spaced points spread a gap's
   threshold over its width, which interpolates between training values. The alternative is to apply the
   rule only to a column with a repeated value. Recommended: as ruled. It is the BART package's rule, the
   quantile rule already does this to every column within its count, the loss is a few percent on fits of
   under ~90 rows and nil after, and a rule keyed on a repeated value would make a column's grid depend on
   whether one value happens to repeat. Evidence that would change it: a larger loss on real small data.
3. The gap point's guard, and the quantile rule. Two adjacent doubles have no double strictly between them,
   so their midpoint rounds onto one of them; today the quantile rule gives four adjacent doubles two points
   and leaves a pair inseparable (ran). Recommended: the guard of The rule, in the one fill both rules use, so
   the quantile rule also changes, on such columns only (no recorded baseline holds one). The alternatives:
   the guard under the default rule only, leaving the quantile rule as it is; or no guard, under which the
   default rule also gives four adjacent doubles two points and bd-balance's narrow arm fails.
4. bartBT. It shares the store, so `bartBT` and `bart`'s forwarded BayesTree spelling follow, departing
   from BayesTree, whose evenly spaced points this rule replaces; with `factors = "indicators"`, its default,
   every factor's indicator columns move. Recommended: follow, as dec-B235 did for the quantile grid ("Fix it
   everywhere, bartBT included."): one spelling, one grid.
5. bench-sampler.R compare, owed for a hot-path change, waived for this slice in favour of the creation and
   refresh timings in Speed, no bench arm reaching it. Recommended: waive.

## Estimate

Implementer about half a day; the re-record, exact gates in full mode for the four arms and the consumer
suites about four hours of machine time, the timings about half an hour on a quiet machine.
