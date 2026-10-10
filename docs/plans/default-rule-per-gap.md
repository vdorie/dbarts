# default-rule-per-gap: one cut point per gap, each gap weighted by its width

Status: LANDED 2026-10-09 (f946d728..650cdecb; Landing note below). Planned under dec-B406, which revises
dec-B313; replanned the same day after the blind critique. Claims below marked (ran) were run on a
prototype outside the tree.

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants below.
rng: POSTERIOR-CHANGING for every fit under the default rule (`useQuantiles = FALSE`) with a numeric column
of at least two and fewer than `n.cuts` distinct finite values, at creation, `setData`, a refresh and in
xbart's folds: the tree prior over that column's thresholds changes (Context). NEUTRAL, bit for bit, for
every other fit: the quantile rule (but the adjacent-doubles columns of Open call 1), columns with at least
`n.cuts` distinct values, constant and all-missing columns, factor columns, grids set by `setCutPoints`, and
states, copies and reloads, which install the grid and weights they carry. The weighted pick is new code on
the hot path; an unweighted column runs today's code and draws.
window: before the merge to main (dec-B313, dec-B406), if VD accepts the size; an engine slice of its own,
serial with the queue (Sizing).
budget: ~1000 lines planned (Sizing), forecast 1500 to 2000.

## Goal

Under the default rule a numeric column with fewer distinct values than `n.cuts` holds one cut point in each
gap between neighbouring values, and the tree prior and every proposal choose among a node's gaps with
probability in proportion to each gap's width. A column with `n.cuts` or more distinct values keeps its
`n.cuts` evenly spaced, equally likely points. The engine gains weighted cut choice, stored per column beside
the grid and carried by a state, in a form TODO cut-point-weights can reuse. The tier is "Changes draws"
([Process by risk](README.md#process-by-risk)).

## Context

The rulings. dec-B313 (VD 2026-10-07): "Use option B, One point per gap when a column has fewer distinct
values than n.cuts; evenly spaced otherwise", after "It seems inefficient from a mixing perspective to have so
many possible proposals that don't change anything." dec-B406 (VD 2026-10-09), shown the critique's
measurements, the rule as ruled, a cutoff at about 8 distinct values and width-weighted gaps: "We can do
width-weight gaps, but yes, I'd want it planend and sized first." The entry's rule, its text and not VD's:
one point per gap where a column has fewer distinct values than `n.cuts`, "each gap chosen with probability in
proportion to its width, in the prior and the proposal; planned and sized before it is built, and brought
back if too large for 1.0." Unchanged and kept: dec-B297 (a grid holds each value once), dec-B298 and dec-B299
(a refresh derives "Up to n.cuts, from the new values" and shrinks), dec-B300 (no repeats; "Weights later
come as probabilities beside the grid, never as repeats"), dec-B231 (a grid is data), dec-B311 and dec-B312
(where splits go when a grid is replaced), dec-B172 (finite values only). TODO cut-point-weights: weights as
`setCutPoints(cuts, column, probs = )`, after 1.0-0 (read).

Why weights (ran). The critique emulated equal per-gap grids: on three zero-inflated columns (12% nonzero,
n = 200) test RMSE rose 5.5% (+0.0145, se 0.0039), on continuous columns of 30 and 60 rows by 3% and 2.5%, and
the equivalence movers shifted by up to |z| 32 (mixedmatrix). An equal per-gap prior moves the zero/nonzero
gap's share from about a third to 1/U. The prototype weights gaps by width; with the same seeds and designs,
default fit on base against prototype (scratch/drplan/acc.R), paired difference in test RMSE (se):

| design | base | weighted per gap | difference | equal per gap (earlier) |
|---|---|---|---|---|
| zero-inflated, n 200, 10 seeds | 0.2646 | 0.2651 | +0.0006 (0.0018) | +0.0145 (0.0039) |
| Friedman, 30 rows, 8 seeds | 2.944 | 2.954 | +0.010 (0.017) | +0.089 (0.042) |
| Friedman, 60 rows | 1.956 | 1.962 | +0.005 (0.014) | +0.048 (0.017) |
| Friedman, 90 rows | 1.805 | 1.797 | -0.008 (0.019) | +0.011 (0.024) |
| 0/1 and five-value columns, n 300, 10 seeds | 0.4397 | 0.4408 | +0.0011 (0.0022) | +0.001 (0.002) |

The losses go. Against equivalence-51201107 (z mode, the eight movers, ran) the prototype's largest shifts are
hurdle 8.95, sparse 5.40, mixedmatrix 3.35 (32.4 under equal gaps), and wideFactorIndicators 2.52,
factorpartial 2.29, leaffactormixed 2.21, hazard 2.91, xbartmixed 0.63.

What the prior keeps and what it changes (ran). At a node the chance of a gap is its width over the summed
widths of the node's gaps, which is the share of evenly spaced points falling in it to within one point in
`n.cuts`. What changes is availability: a node may split on a column while its interval holds a point
([`Tree::variableAvailable`](../../src/bartcore/tree.hpp)). Today a node holding one value of a 0/1 column
still holds points, counts the column as available, gives it a share of its variable draws, and every such
proposal empties a side; per gap it holds none. In one tree on a 0/1 column, noise response, the share of
draws holding a split is 0.177 at one point and 0.131 at 100 (scratch/drplan/prior-mass.R). That change is
the rule's purpose and stays; the residual hurdle and sparse shifts above are of this kind.

Where the grid comes from (read). [`ColumnStore::deriveNumericCuts`](../../src/bartcore/data.hpp) is the one
derivation, called from [`ColumnStore::buildCutsForColumn`](../../src/bartcore/data.hpp) (creation,
[`ColumnStore::setData`](../../src/bartcore/data.hpp), the xbart data handle) and
[`ColumnStore::refreshCutsForColumn`](../../src/bartcore/data.hpp) /
[`ColumnStore::refreshCutsForCscColumn`](../../src/bartcore/data.hpp) (a whole-column `setPredictor` with
`updateCutPoints = "position"` or `"value"`; the per-observation forms refuse a cut-point update, "partial
updates cannot also update cut points" in [`bartcoreSamplerSetPredictor`](../../R/bartcore.R)). The other
writers of a grid (critique finding 8, read): [`ColumnStore::setCutPointsForColumn`](../../src/bartcore/data.hpp)
(`setCutPoints`, and [`setState`](../../src/bartcore/sampler.hpp)'s install from the state), a row-subset view
copying its parent's grid (xbart folds), [`ScopedCutGrid`](../../src/bartcore/data.hpp) (a warm start's
donor grid, temporary), and the restores on a refusal
([`revalidateAllChains`](../../src/bartcore/sampler.hpp), `setState`'s rollback). The flat C API has no entry
that sets a grid ([`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h)); stan4bart (bartcore 963956b)
calls no grid setter.

Where the engine counts or picks a cut (read). The rule's prior is uniform over the positions of the node's
interval ([`Tree::splitInterval`](../../src/bartcore/tree.hpp)) in
[`CGMTreePrior::ruleForVariableLogProbability`](../../src/bartcore/model.hpp), which
[`CGMTreePrior::treeLogProbability`](../../src/bartcore/model.hpp) sums, so every move scoring a subtree
(swap, change, perturb) reads it; the draw is
[`CGMTreePrior::drawRuleForVariable`](../../src/bartcore/model.hpp) (birth, and `sampleTreesFromPrior` through
`drawRuleAndVariable`). Four places count positions themselves:
[`changeMove`](../../src/bartcore/moves.hpp) (its forward and reverse interval and valid counts and its draw
over the valid range from [`findGoodOrdinalRules`](../../src/bartcore/moves.hpp)),
[`perturbMove`](../../src/bartcore/moves.hpp) (symmetric over window positions; the node's own rule factor
cancels only while positions are equally likely),
[`enumerateNogRuleNeighbourhood`](../../src/bartcore/moves.hpp) (rule_gibbs, `logRulePrior`) and
[`growTreeFromRoot`](../../src/bartcore/grow.hpp) (`logCut`). Availability, the monotone leaf's moves and the
interaction constraint read only whether an interval is empty or which variable a node holds, and do not
change. The scan scores positions and does not change.

What counts as a distinct value (read): as the quantile collectors count, finite values only, `+0.0` so -0 and
0 are one, a CSC column's implicit zeros once, over every training row (weights, an active-row mask, offsets
and test rows do not enter; the grid is data, dec-B231).

0.9-34 and others (ran unless marked). 0.9-34 gives the default rule 100 points on 0/1, 13-value, six-value and
continuous columns, the quantile rule 1, 12, 5 and 100, as the branch does today. The BART package's
`bartModelMatrix` takes midpoints, unweighted, when a column has fewer distinct values than `numcut` (source
read; its boundary is the words'). BayesTree spaces points evenly, bartMachine draws a threshold from the
node's distinct values, bcf takes midpoints under 8 distinct values (as recorded in dec-B313).

## The rule

For numeric column j under the default rule, with U distinct finite values u_1 < ... < u_U and n =
`requestedNumCuts[j]`:

- 2 <= U < n: U - 1 points, point k in the gap (u_k, u_(k+1)), with weight w_k = u_(k+1) - u_k. Where every
  w_k is equal (a 0/1 column, equally spaced integers) the column carries no weights: equal weights are the
  uniform choice, and the unweighted path draws as today's does.
- Otherwise: today's grid, unweighted (n evenly spaced points; one point over U of 0 or 1).

The point is the midpoint 0.5 (u_k + u_(k+1)), the quantile rule's expression, so the grid is the quantile
rule's to the bit; where that midpoint is not strictly below u_(k+1) (adjacent doubles, or a sum that
overflows) the point is u_k (Open call 1). A difference never rounds to zero for u_k < u_(k+1); where one
overflows, every weight of the column is taken as 0.5 u_(k+1) - 0.5 u_k instead, the same proportions.

At a node whose interval holds positions a to b of column j, rule k has prior probability
w_k / (w_a + ... + w_b), and an unweighted column 1 / (b - a + 1) as today. Every proposal that draws a
position draws it from the same weights over the range it draws from, and every count in an acceptance ratio
becomes the matching sum of weights.

What the ruling covers and what it does not:
- Covered: the default rule's derivation at creation, `setData`, a refresh and the xbart handle.
- Not covered, kept unweighted: the quantile rule (its midpoints already stand for equal shares of the
  distinct values; Open call 3), a grid from `setCutPoints` (until cut-point-weights adds `probs =`;
  Open call 4), columns at or past the boundary.
- Carried, not derived: a state's weights (State), a view's (copied with its grid). A warm start maps the
  donor's splits onto the recipient's grid by value and uses the recipient's weights.

## Change

Engine ([`ColumnStore`](../../src/bartcore/data.hpp) unless named):

1. Storage: per column, weights and their prefix sums (`numCuts[j] + 1` entries from 0), empty for an
   unweighted column; accessors for the mass of positions [a, b] (the count when unweighted), one position's
   weight, and the position at a given mass. `setCutPointsForColumn` takes an optional weights pointer (null
   unweighted), the entry cut-point-weights will call; the view constructor copies weights with the grid;
   `ScopedCutGrid` installs a donor grid unweighted and restores the column's weights with its grid.
2. Derivation: a bounded distinct collector, dense and CSC (stops past n - 1 distinct values; expected linear
   time, scratch sized to n, never to the row count); within the boundary the quantile fill
   ([`ColumnStore::fillCutsFromQuantileGrid`](../../src/bartcore/data.hpp) with `inducedNumCuts = U - 1`)
   and the widths; otherwise today's uniform fill.
3. The gap point's guard in the one fill both rules use (Open call 1).
4. Prior and draws: `CGMTreePrior::ruleForVariableLogProbability` adds log w_k - log(sum) on a weighted
   column; a shared position draw replaces `ext_rng_simulateIntegerUniformInRange` in
   `CGMTreePrior::drawRuleForVariable` and `changeMove`, unweighted columns keeping that call and its
   generator use.
5. Moves: `changeMove`'s four counts become sums; `perturbMove` adds log w_target - log w_current;
   `enumerateNogRuleNeighbourhood` and `growTreeFromRoot` add log w_c to each candidate and take the
   interval's sum as normalizer.
6. Rollbacks: every snapshot that saves `cutPoints` and `numCuts` saves the weights
   (`setState`'s, `revalidateAllChains`'s, the refresh undo in [`Chain::undoSplitMoves`](../../src/bartcore/chain.hpp)'s
   caller); `setState`'s skip of an equal grid compares weights too.

Bridge ([`storeState`](../../src/R_interface_bartcore.cpp), [`setState`](../../src/R_interface_bartcore.cpp),
[`SamplerStateData`](../../src/bartcore/sampler.hpp)): State, below. No facade virtual, no R surface, no
`dbarts.h` change.

The prototype carries changes 1, 2 (unbounded collector) and 4 to 5 in about 100 changed lines (ran:
`diff` against the base engine), not 3, 6, the state or the CSC-specific paths beyond the shared collector.

## State

A state's grid comes from its sampler and is installed by `setState`, a copy and a reload (read), so weights
must travel with it: a column created few-valued and later given other values with `updateCutPoints = "none"`
keeps the old grid and weights, and a reload re-derives from the new values, so neither the grid nor the
weights can be re-derived. `storeState` writes a top-level attribute `cutWeights`, a list with one entry per
column, NULL for an unweighted one; `setState` installs it with the grid and refuses, naming the column, an
entry of the wrong length, or one not positive and finite. Under the registry rule at
[`stateFormatVersion`](../../src/R_interface_bartcore.cpp) this is an additive attribute: no version bump; a
state without it installs unweighted grids, as every state written before this slice holds. No released
format exists.

## Tests

tests/cpp:
- [test_data.cpp](../../tests/cpp/test_data.cpp): [`testDistinctCutGrids`](../../tests/cpp/test_data.cpp)'s
  counts become {4, 1, 1, 1, 1, asked} under both rules (the five adjacent doubles: four points under the
  guard, today five under the uniform rule, the last separating nothing (ran), and three under the quantile
  rule (read)), the 0/1 column {0.5} unweighted. New `testDefaultRuleWeightedGaps`: U from n - 2 to n + 1 at
  n = 10 straddles the boundary; weights equal the widths and an equally spaced column carries none; dense and
  CSC equal; NaN and `Inf` do not count; -0 and 0 one value; per-column n; overflow widths; the guard
  (four adjacent doubles: three points, four codes). [`testRefreshDerivesAsCreation`](../../tests/cpp/test_data.cpp)
  gains a weighted column across creation, `setData`, a refresh and back.
- [test_tree.cpp](../../tests/cpp/test_tree.cpp): `ruleForVariableLogProbability` and `treeLogProbability`
  on a weighted column against hand sums; the position draw's frequencies against its weights (chi-square).
- [test_moves.cpp](../../tests/cpp/test_moves.cpp): `changeMove`'s and `perturbMove`'s log corrections and the
  rule_gibbs candidate weights on a weighted column against enumeration; [test_grow.cpp](../../tests/cpp/test_grow.cpp)'s
  `logCut` check gains a weighted case.
- [test_state.cpp](../../tests/cpp/test_state.cpp) or test_sampler.cpp: weights survive store, `setState`,
  a refused update's rollback and a warm start's scoped grid.

tinytest:
- [test-cut-grid-distinct.R](../../inst/tinytest/test-cut-grid-distinct.R): the counts
  retired: ["c(100L, 5L, 1L, 1L, 100L, 100L)"](../../inst/tinytest/test-cut-grid-distinct.R) become
  `c(100L, 4L, 1L, 1L, 1L, 5L)` and the quantile counts the same; the narrow pins at the creation, refresh and
  forced-refresh sections go from 5 to 4 (critique finding 16).
- [test-cut-points-undo.R](../../inst/tinytest/test-cut-points-undo.R): its narrow count 5 becomes 4.
- [test-xbart-fold-oracle.R](../../inst/tinytest/test-xbart-fold-oracle.R): its premise, that a `dbarts()`
  fit on a fold's training rows rebuilds the fold's full-data grid, holds only while the grid depends on the
  range alone; restate it with a fixture whose columns hold at least `n.cuts` distinct values on the fold's
  rows, and add a check that a fold's grid and weights are the full data's. A fold's grid then holds gaps
  between held-out values, points that move no training row, as today's evenly spaced grid does.
- test-sampler-prior.R (a seeded tolerance on hillData, whose z is 0/1) and test-joint-update-factor.R (a
  `declined > 0` on a 40-row continuous fixture) move with the draws (critique finding 4, ran on the equal-gap
  emulation); the implementer re-derives each pin and lists every changed expectation with its reason.
- New: a 0/1 column draws bitwise the same at `n.cuts` 3 and 100 (one point, unweighted); a weighted column's
  state, copy and reload continue the chain; a state without `cutWeights` installs unweighted; `setCutPoints`
  drops weights except for the held grid (Open call 4); a refresh re-derives them; `"none"` keeps them.
- [test-quantile-grid.R](../../inst/tinytest/test-quantile-grid.R): through `dbarts`, `bart`, `bartBT`, a
  BayesTree-spelled `bart` call and `rbart_vi`, a few-valued column's grid equals the quantile rule's.

## Baselines

Counted with a tracer on every entry that hands a store to the engine (scratch/drplan/hook.R) under the words'
boundary (ran):
- equivalence-51201107 ([MANIFEST](../../benchmarks/baselines/MANIFEST)), 55 scenarios: 8 move
  (factorpartial, hazard, hurdle, leaffactormixed, mixedmatrix, sparse, wideFactorIndicators, xbartmixed;
  sizes above); 47 do not, bart2twoforest and quants among them. bcf-equivalence-1b7d730c (15) and
  multinomial-equivalence-80b1c8d4 (11): no affected column.
- The four test-reproducibility-*.R files: none moves (singleThreaded and xbart fit columns of exactly 100
  distinct values, binaryResponse and multithreaded no few-valued column).
- Exact gates: draws move in backfit-exact, heteroscedastic-exact and multinomial-exact (0/1 columns, one
  point), and in hazard-reduction on both of its sides alike, so it stays bitwise; bd-balance's narrow arm
  keeps its grid (three points on adjacent doubles under the guard, equal widths, unweighted; ran) and its
  draws. Not moving: bcf-latent-exact (K = 2 and 3 cells with `n.cuts` K - 1; critique finding 7) and
  monotone-successive-conditional (100 distinct values at `n.cuts = 100`, which the critique's finding 3 had
  moving under the boundary this plan no longer takes). The balance gates (change, swap, perturb,
  rule-gibbs) use the quantile rule and do not move. No existing gate holds a weighted column.
- SBC: bcf-probit-weak and bcf-logistic-weak (40 rows), gp, gp-weighted and gp-mixed (80), monotone-1-leaf,
  monotone-1-joint and monotone-bd (20) move (ran); none is in sbc.yaml's matrix, whose twelve do not move
  (discrete-selfcheck builds no sampler).
- Re-record: the 8 movers fresh on the reference build, merged into a copy of 51201107, named after the code
  commit, MANIFEST row, 51201107 demoted; partition in z mode, the |z| reported as sizing; 55 of 55 under
  `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (P17), all new or re-run, each on the base and the slice:
  - bd-balance gains a `weighted` arm: cells at 0, 1, 3 and 10, default `n.cuts`, the grid {0.5, 2, 6.5}
    asserted, and the enumeration with log w_j - log(sum) for its rule factor, run under five mixtures
    (birth/death, change, swap, perturb, rule_gibbs) so each kernel's weighted path faces the exact posterior.
    The prototype passes all five (largest |z| 1.7, ran, scratch/drplan/bd-weighted.R); against the
    equal-gap oracle it fails at |z| up to 162 (ran), so the arm sees the weights. On the base build it fails
    at its grid assertion.
  - monotone-exact-enumeration gains an uneven-values arm, its rule factor weighted likewise.
  - multinomial-exact arm 8 and heteroscedastic-exact part (b) in full mode, against the gate's oracle and
    the one with the children's factor (critique finding 5): on the slice the gate's own oracle becomes the
    exact one.
  - growTreeFromRoot's weighted candidates: the C++ test above (no R gate reaches grow-from-root).
  - The monotone-bd SBC arm runs as a no-breakage check, not as the oracle (calibration holds on any grid).
- `Rscript benchmarks/R/mutation-battery.R verify-anchors` passes; the mutants whose killer is a moved
  scenario still kill.

## Speed

The weighted pick is in the hot path. Measured on the prototype (ran, scratch/drplan/speed.R, shipped
builds, 75 trees, 5000 sweeps after 200, five interleaved rounds on a machine at load 15 to 19, so a 5% grain):
uniform columns (unweighted path) 0.99 base by minimum and 1.04 by median; few-valued uneven columns
(weighted path) 1.005 and 0.97; 0/1 columns 1.05 and 1.05, where the grid and posterior differ and the code
does not: the prototype's 0/1 and equally spaced columns draw bit for bit what the base draws once given the
same per-gap grid by `setCutPoints` (ran, scratch/drplan/ident.R). Owed on a quiet machine: bench-sampler.R
compare (its arms are uniform columns, the unweighted path; accept 1.03), a weighted arm (accept 1.05 against
the base's evenly spaced fit of the same data), creation at n = 1e6 with 50 uniform columns (the collector
stops early; 1.05) and with 50 few-valued ones (reported), and a refresh loop of 200 calls at n = 1e5 on a
uniform and on a few-valued column (critique finding 15; 1.05 and 1.5).

## Docs

- [cut-grid.md](../design/cut-grid.md), the design note: [The rule](../design/cut-grid.md#the-rule) states
  the per-gap grid, the weights, the guard and the state attribute; [What moves](../design/cut-grid.md#what-moves)'s
  0/1 and six-value rows and its forward pointer to dec-B313 are restated; a dated section records the
  measurements above and the availability effect; [Checks](../design/cut-grid.md#checks) names the new arms.
- [quantile-grid.md](../design/quantile-grid.md): [What does not move](../design/quantile-grid.md#what-does-not-move)'s
  "The uniform grid, which is the default." points here; "A cut is a midpoint, never an observed value" gains
  the guard's exception.
- [data-store.md](../design/data-store.md)'s `cutPoints` bullet (stale already: it says a grid can repeat a
  value) gains the weights; [sparse-columns.md](../design/sparse-columns.md)'s uniform-grid sentence.
- [classic-compare.md](classic-compare.md): an entry under
  [What differs, and why](classic-compare.md#what-differs-and-why); its numbers at the release-candidate run.
- Help: man/dbartsControl.Rd (`useQuantiles`, `n.cuts`), man/bart.Rd, man/bartBT.Rd (Decision Rules),
  man/xbart.Rd, man/dbartsSampler-class.Rd (`setCutPoints`: a grid given there is unweighted; `setPredictor`'s
  `"none"` keeps the grid and weights a column has; `storeState`'s attribute). The text: under the default
  rule a numeric predictor with fewer distinct values than `n.cuts` gets one cut point halfway between each
  pair of neighbouring values, chosen in proportion to the gap's width.
- NEWS, beside the item on distinct cut points: under the default `useQuantiles = FALSE`, a predictor with
  fewer distinct values than `n.cuts` (a 0/1 column, a small count, a sparse column, any column on fewer rows
  than `n.cuts`) gets one cut point halfway between each pair of neighbouring values, each chosen in
  proportion to the gap's width, where 0.9-34 placed `n.cuts` equally spaced points; a 0/1 column has one
  cut point where it had 100, no proposal moves among points that split the same rows, and fits with such a
  column change.

## Records

At landing: the ledger entry for the calls made; Status and Landing note; TODO's default-rule-per-gap item
out; TODO cut-point-weights restated to what remains (the `probs =` argument, its validation and help); and
two new TODO items for the gate oracles the critique ran (finding 5), whether or not this slice lands, since
each oracle holds only while its column keeps no point in a split's children:
- heteroscedastic-exact-split-prior: part (b) weighs a two-cell variance or mean tree by `base` alone; on a
  0/1 column of 100 points the children stay available and the prior carries (1 - base / 4)^2 more (gaps
  0.0027 and 0.0012 at 100 points against 0.0013 and 0.0001 at one, critique run). Fix: assert the column's
  grid is one point, or put the factor in the oracle.
- multinomial-exact-arm8-split-prior: arm 8 (`armConstrained`, two 0/1 numeric columns, `max.order = 1`,
  base 0.95) weighs a split by `base / 2` alone; at `n.cuts = 100` the sampler misses that oracle by 0.0017
  and matches the one with the children's factor to 0.0003. Same fix.

## Steps

1. Changes 1 to 6 with the tests/cpp tests; `make && ./test_bartcore` green, and under
   `-fsanitize=address,undefined`.
2. The state attribute and its tests; `R CMD INSTALL --preclean -l <lib> .`; the tinytest changes; the full
   suite green.
3. The `weighted` bd-balance arm and the monotone arm; each run on the base (fails) and the slice (passes);
   the workflow loop gains them.
4. Docs, help, NEWS.
5. After review, in their own commit: the re-record and MANIFEST row, the oracle runs, the speed runs.

## Gates

Posterior-changing ([RNG classes and their gates](README.md#rng-classes-and-their-gates)), on the slice tip in
its own library, independently of the implementer: tests/cpp plain and sanitized; the R-loaded ASAN path over
the grid, state and xbart test files; the full tinytest suite; `R CMD check --as-cran`; lintr, air, rc-codoc,
win-drift, doc-freshness, the NEWS parse; the three equivalence harnesses and the four snapshot files on the
reference build; exact-gates.yaml's list in quick mode with the new arms, full mode for the new arms,
backfit-exact, heteroscedastic-exact and multinomial-exact; Speed; consumers against the slice's library
(stan4bart's suite and posterior baselines, bartCause, treatSens, bairrtt: bartCause's response fit and
treatSens carry the 0/1 treatment as a predictor and bairrtt its 0/1 `z`, read, so their draws move; bairrtt's
latent columns take its own grid through `setCutPoints` and do not).

Reviewer's mutants, each of which must fail a test:
- the weight dropped from each site in turn: the prior's rule factor, the birth draw, change's counts, change's
  draw, perturb's ratio, rule_gibbs's candidates, grow-from-root's candidates;
- weights as equal (uniform choice on a weighted column): the `weighted` arm fails;
- the boundary off by one (`U <= n` for `U < n`);
- the guard dropped (four adjacent doubles hold two points; bd-balance's narrow arm fails its grid assertion);
- weights not copied to a view, not saved by a rollback, not written or not read by the state;
- the per-gap branch on the dense path only (the CSC twin differs).

## Sizing

Planned lines, from the prototype (about 100 engine lines for the core, ran) and the list above:

| part | lines |
|---|---|
| data.hpp: storage, collector, derivation, guard, view, scoped grid, setter | 130 |
| model.hpp, moves.hpp, grow.hpp: prior, draw, change, perturb, rule_gibbs, grow | 60 |
| sampler.hpp, chain.hpp: rollbacks, equal-grid skip | 40 |
| bridge: state attribute, refusal | 70 |
| tests/cpp | 300 |
| tinytest | 170 |
| gate arms (bd-balance five mixtures, monotone) and workflow | 120 |
| docs, help, NEWS, MANIFEST | 130 |
| total | ~1020 |

Forecast at the usual 1.5 to 2 times: 1500 to 2000 lines; implementer two to three days, review and fix
round one to two, gates about eight hours of machine time, two of them quiet. For scale: repeated-cut-restore
was planned at 1400 and built near its 2600 upper figure; leaf-conversions planned at 1270, built at 2900.

Does it fit before 1.0-0: yes, as one engine slice in the serial queue, about the size of one install-surface
slice and smaller than the forest arc. Order: after gp-copy-continuation, before state-frame-prior's grid
half and state-grid-replaces-sampler, which change how a state carries a grid and would otherwise be
re-planned around `cutWeights`. monotone-exact-birth-death (design first) would touch the same moves;
whichever comes second rebases.

The fallback, bcf's cutoff (one point per gap, unweighted, only where a column has at most about 8 distinct
values): the derivation change alone, no weighted machinery, state or new kernel paths; about 500 lines
planned, forecast 750 to 1000, one engine day. It moves 0/1 and small-count columns and leaves alone the
zero-inflated, the sparse and the small-n continuous columns where the equal-gap prior lost accuracy; within
2 to 8 values an uneven count column's gaps become equally likely (0, 1, 2, 5, 10, 20: the 10-to-20 gap
from about half the prior to a fifth), not measured. It builds nothing cut-point-weights can reuse.

## Stop conditions

Stop and report when: the diff passes 2000 lines before machine-written values; a kernel's weighted arm fails
on two seeds after review; a scenario, snapshot or harness outside Baselines moves; bench-sampler's compare
on the unweighted path stays above 1.03 or the weighted arm above 1.05 after one round of levers (a branch
hoisted per column, the prefix lookup inlined); the change needs a facade, C API or state format version
change.

## Interactions

- gp-copy-continuation (in flight): disjoint code; both re-record the equivalence baseline on disjoint
  scenarios; the second records against the other's file.
- The install-surface slices and state-frame-prior's grid half: `setState` is touched by both; land this
  first or rebase it over them with the attribute carried.
- cut-point-weights (after 1.0-0): reuses the storage, the setter's weights pointer, the state attribute, the
  weighted kernels and the `weighted` arm; adds only `setCutPoints(probs = )`, its checks and help.

## Open calls

Settled by dec-B406 and the prototype, recorded for the reader of the earlier plan:
- The boundary: the words, U < n. With widths as weights a per-gap grid and today's evenly spaced one carry
  nearly the same gap prior, so the earlier case for U - 1 <= n (two values of U) is small, and the words
  leave the 100-row fixtures (friedmanData.R in 54 test files, two snapshot files, monotone-successive-
  conditional) and `useQuantiles` on 101-row data alone (critique finding 1).
- Continuous and zero-inflated columns under the boundary: covered, the accuracy loss gone (table above).

Open:
1. The guard, and whether the quantile rule takes it. Two adjacent doubles have no double between them, so
   their midpoint rounds onto one; today the quantile rule gives four adjacent doubles two points and leaves a
   pair inseparable (ran). The guard places the point on the lower value, so a cut equals an observed value,
   against quantile-grid.md's "never an observed value". Recommended: the guard in the one fill, so the
   quantile rule changes on such columns only (no recorded baseline holds one; tinytest pins do). The
   default rule needs it either way, or bd-balance's narrow arm and the narrow pins lose a point. Alternative:
   the default rule only.
2. bartBT. It shares the store, so `bartBT` and a forwarded BayesTree-spelled `bart` follow, departing from
   BayesTree's evenly spaced points; with its default `factors = "indicators"` every factor's indicator
   columns move. dec-B235's "Fix it everywhere, bartBT included." fixed a defect and does not decide this
   (critique finding 14). Recommended: follow. Under width weights the gap prior is BayesTree's to within a
   point; what differs is proposals that move no row, which no caller relies on; a BayesTree grid in bartBT
   would need a rule flag through the data object (about 60 lines) to keep a behaviour the ruling calls
   wasteful. Alternative: bartBT keeps evenly spaced points.
3. The quantile rule unweighted. The ruling names the default rule. Recommended: unweighted; its midpoints
   are chosen as equal shares of the distinct values, which is what `useQuantiles = TRUE` asks for.
   Alternative: weight its every-midpoint case too, making the two rules agree on few-valued columns.
4. `setCutPoints` and weights. Recommended: a grid given to `setCutPoints` is unweighted until `probs =`
   exists, except the grid the column holds bit for bit, which keeps its weights, so that call stays a
   no-op. Cost: undoing a `setCutPoints` by handing back the derived grid loses its weights; `setState`
   restores both. Alternative: always unweighted.
5. Go, or the fallback. Recommended: build the weighted rule now (about 1000 lines, forecast 1500 to 2000):
   it is what was ruled, it measured no loss on any design tried, and it builds cut-point-weights' machinery
   once. The cutoff is half the size and leaves the 5.5% and 3% losses' columns unchanged, which avoids them
   too, but changes uneven small counts unmeasured and builds nothing reusable.

## Landing note

Landed 2026-10-09 as f946d728..650cdecb after the plan's blind critique, one opus review with 31 mutants
and one fix round the same reviewer checked. Under the default rule a numeric column with at least two and
fewer than n.cuts distinct values holds one cut point per gap between its values, each gap chosen in
proportion to its width in the tree prior and every proposal (birth, change, perturb, rule_gibbs,
grow-from-root); states carry the weights as an additive cutMass attribute, older states reload unchanged.
Equivalence re-recorded as equivalence-e4faed5c: the eight planned movers (hurdle max |z| 8.95, sparse
5.40, mixedmatrix 3.35, hazard 2.91, wideFactorIndicators 2.52, factorpartial 2.29, leaffactormixed 2.21,
xbartmixed 0.63), the other 47 bitwise; oracle in its MANIFEST row (bd-balance's weighted arm under five
move mixtures, monotone design u1, multinomial-exact arm 8 and heteroscedastic-exact part (b), tests/cpp
testWeightedCutKernels). bcf-equivalence-1b7d730c 15/15 and multinomial-equivalence-80b1c8d4 11/11
bitwise; the four snapshot files unchanged (27 results). After the final rebase: equivalence 55/55 bitwise
on the reference build, the suite in one process 19693 results and 0 failures, lint chain clean,
stan4bart's suite 582 results and 0 failures. tests/cpp plain and under address/undefined sanitizers,
R-loaded sanitizers on the touched files, every exact gate quick and the touched ones full, R CMD check
--as-cran one NOTE (Date). Mutants: 30 of 31 killed by tests/cpp or tinytest, change's forward side as
counts only by bd-balance until the fix round reshaped testWeightedCutKernels so tests/cpp kills it too.
Speed on the Mac under load: bench-sampler 0.976 to 1.026 (bound 1.03), weighted columns 1.003, creation
with uniform columns 1.00, refresh 1.015, few-valued creation 1.21 and refresh 1.13 (reported, bound 1.5);
the quiet x86 A/B after landing (one thread, A and B interleaved, seven rounds): bench-sampler geometric mean 0.999 both directions (worst metric 1.041), weighted columns 1.002, creation with uniform columns 1.014, refresh 1.025, few-valued creation 1.156 (reported) and refresh 1.196 (bound 1.5); every bounded arm passes. The bart-as-a-component vignette's recipe 4 was corrected (a
revert can be declined) and its handling queued (TODO embedding-recipe-declined-revert). Calls: dec-A194.
