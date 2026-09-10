# Memory footprint

Status: VALIDATED against measured peak resident set size by
[benchmarks/R/memory-footprint.R](../../benchmarks/R/memory-footprint.R) on
arm64 macOS, first at 80ff32b4 and re-validated after the transient copies
came out of packaging. Every one of the 30 grid cells falls within
max(10 pct, 20 MB) of the model, and the median absolute relative residual
over the whole grid is 3.0 to 4.7 pct across runs against a 5 pct limit. The
top of that spread, and the one linear-cell miss seen with it, came from a
run that overlapped another job on the host; on a quiet machine the median is
3.0 to 4.2 pct with no cell outside. The worst cell on a quiet run is
n = 1e5, p = 20, T = 200, C = 1, S = 200, keepTrees, 39.6 to 39.7 MB high
against a 55.8 MB tolerance, which is the collector churn the model no longer
carries a row for. The margin on
the median is thin, and it is thin for a stated reason: sixteen of the thirty cells sit at n = 1e4, where a
prediction of 20 to 60 MB is scored against page-level allocator behaviour a
byte model cannot reach. "What the measurement moved" records the rows the
first run changed and "What the removal moved" the rows the second did.
Two rows have moved since that validation without a re-recording - the
ingestion allowance and the starting-sigma row, both cut by the follow-ons
in "The largest avoidable allocations, ranked" - so a re-record of
benchmarks/baselines/memory-footprint-b184b6b2.csv against the model as it
now stands is OWED. Both moves are downward and neither touches the engine
rows, so the grid's residuals can only have grown more negative.

The closed-form model of what a fit allocates, per component, in the units the
allocations have. Every byte count is read off the element type and the
allocation site, not measured. Symbols: n rows, p predictors, T trees per
forest, C chains, S kept draws, nTest test rows, K categories, q designated
leaf covariates, nnz sparse nonzeros, m the mean live node count of a tree (the
one non-derived input, below). Recurring widths, each confirmed by `sizeof` on
arm64: [`xint_t`](../../src/bartcore/data.hpp) 2,
[`index_t`](../../src/bartcore/data.hpp) 4, [`Node`](../../src/bartcore/tree.hpp)
56, [`FlatNode`](../../src/bartcore/tree.hpp) 24. No row differs by build mode:
the reference build re-associates draw-path kernels and allocates nothing
differently.

## The table

| component | symbol | scope | bytes/unit | units | exists when |
| --- | --- | --- | --- | --- | --- |
| train cut codes | [`CodeBlock::codes`](../../src/bartcore/data.hpp) under [`ColumnStore::train`](../../src/bartcore/data.hpp) | sampler | 2 | n*p | dense-stored columns (always, off CSC rank storage) |
| test cut codes | [`ColumnStore::test`](../../src/bartcore/data.hpp) | test row | 2 | nTest*p | a test set |
| cut grid | [`ColumnStore::cutPoints`](../../src/bartcore/data.hpp) | sampler | 8 | up to n.cuts (default 100) per numeric column | always |
| column metadata | [`ColumnStore::numCuts`](../../src/bartcore/data.hpp), [`ColumnStore::categoryCounts`](../../src/bartcore/data.hpp), [`ColumnSource`](../../src/bartcore/data.hpp) | sampler | ~70 per column per row set | p | always |
| sparse column | [`SparseColumnData`](../../src/bartcore/data.hpp) | sampler | 0.1875*n + 2*nnz | per rank-stored column, replacing its 2*n codes | a CSC column at or below 20 pct density |
| owned dense raw | [`ColumnStore::ownedDenseValues`](../../src/bartcore/data.hpp), [`ColumnStore::ownedTestValues`](../../src/bartcore/data.hpp) | sampler | 8 | n (nTest) per real dense-backed column | CSC/mixed build; test side on any test build |
| gathered leaf raw | [`ColumnStore::gatheredRawValues`](../../src/bartcore/data.hpp) and the leaf's standardized copy of it | sampler | 16 | n*q | a designated-covariate leaf (linear, gp) |
| owned conditioning vectors | [`BartcoreHolder::ownedResponse`](../../src/R_interface_bartcore_common.hpp), [`BartcoreHolder::ownedWeights`](../../src/R_interface_bartcore_common.hpp), [`BartcoreHolder::ownedOffset`](../../src/R_interface_bartcore_common.hpp) | sampler | 8 | 3*n | always, sized in [`createHolder`](../../src/R_interface_bartcore.cpp) whether filled or not |
| owned test offset | [`BartcoreHolder::ownedTestOffset`](../../src/R_interface_bartcore_common.hpp) | test row | 8 | nTest | always, same site |
| predictor-update cache | [`UpdateSessionImpl`](../../src/bartcore/sampler.hpp) | sampler | 4 | n per tree splitting on the column | an open update session, over every chain |
| index buffer | [`Forest::indexBuffer`](../../src/bartcore/combiner.hpp), sliced into [`Tree::indices`](../../src/bartcore/tree.hpp) | chain | 4 | n*T per forest | always |
| leaf map | [`Forest::leafOf`](../../src/bartcore/combiner.hpp), from [`Chain::initForestFitStorage`](../../src/bartcore/chain.hpp) | chain | 4 | n*T per forest | constant leaf |
| dense tree fits | [`Forest::treeFits`](../../src/bartcore/combiner.hpp) | chain | 8 | n*T per forest | non-constant leaf, replacing the leaf map |
| total fits | [`Forest::totalFits`](../../src/bartcore/combiner.hpp) | chain | 8 | n per forest | always |
| residual store | [`Forest::treeY`](../../src/bartcore/combiner.hpp) | chain | 8, or 4 under storage="single" | n per forest | always |
| move scratch | [`Tree::SubtreeSnapshot`](../../src/bartcore/tree.hpp)'s [`indexSegment`](../../src/bartcore/tree.hpp), in [`Forest::scratch`](../../src/bartcore/combiner.hpp) | chain | 4 | high-water n per forest, reached at a root-level change proposal | always |
| live trees | [`Tree::nodes`](../../src/bartcore/tree.hpp) | tree | 56 | T*m per forest, plus 296 per tree for the Tree object | always |
| live leaf values | [`Forest::muByTree`](../../src/bartcore/combiner.hpp) | tree | 8 | T*m per forest, plus 24 per tree | constant leaf |
| rescaled response | [`GaussianResponse::yRescaled_`](../../src/bartcore/model.hpp) | chain | 8 | n | gaussian, aft |
| latent channels | [`ProbitResponse::latents_`](../../src/bartcore/model.hpp) and its working twin | chain | 8 | 2*n | probit, ordinal, logistic, nbinom |
| row mask | [`GaussianResponse::activeRows_`](../../src/bartcore/model.hpp), [`GaussianResponse::composite_`](../../src/bartcore/model.hpp) | chain | 8 | 2*n | an active-row mask installed |
| variance forest | [`VarianceForest::indexBuffer`](../../src/bartcore/chain.hpp) 4*n*Tv, [`VarianceForest::factorByTree`](../../src/bartcore/chain.hpp) 8*n*Tv | chain | 12 | n*Tv | heteroscedastic |
| variance scratch | [`VarianceForest::combinedVariance`](../../src/bartcore/chain.hpp) and three siblings, plus the chain's mean weights | chain | 8 | 5*n | heteroscedastic |
| multinomial glue | [`MultinomialForestCombiner::omega_`](../../src/bartcore/combiner.hpp), [`MultinomialForestCombiner::suffix_`](../../src/bartcore/combiner.hpp), [`MultinomialForestCombiner::combined_`](../../src/bartcore/combiner.hpp) | chain | 8 | 3*n*K, plus n*K more under a category offset | multinomial; and every per-forest row above runs K times |
| ordinal cutpoints | [`OrdinalResponse::gamma_`](../../src/bartcore/model.hpp) | chain | 8 | 2*(K-1) | ordinal |
| test fits | [`Forest::totalTestFits`](../../src/bartcore/combiner.hpp), [`Forest::currTestFits`](../../src/bartcore/combiner.hpp), from [`Chain::resizeTestStorage`](../../src/bartcore/chain.hpp) | test row | 16 | nTest per forest per chain | a test set |
| saved trees | [`Forest::savedTrees`](../../src/bartcore/combiner.hpp) | saved sample | 24 | S*T*m per forest per chain, plus about 40 per saved tree of vector object and allocator header | keepTrees |
| stored state trees | [`storeFlatTrees`](../../src/R_interface_bartcore.cpp) | saved sample | 13 (4 variable + 8 value + 1 flags) | T*m per forest per chain, plus 4 per tree; times S again under keepTrees | storeState called |
| result train channel | [`allocChannel`](../../src/R_interface_bartcore.cpp), [`installChannel`](../../src/R_interface_bartcore.cpp) | saved sample | 8 | n*L*S*C, L the reported locations (K for multinomial, else 1) | keepTrainingFits |
| result test channel | [`testExpr`](../../src/R_interface_bartcore.cpp), filled from [`Results::testFits`](../../src/bartcore/chain.hpp) | saved sample | 8 | nTest*L*S*C | a test set |
| result varcount | [`varcountExpr`](../../src/R_interface_bartcore.cpp), filled from [`Results::variableCounts`](../../src/bartcore/chain.hpp) | saved sample | 4 | p*F*S*C, F the forest axis | always |
| result sigma, k, thresholds, dispersion | [`sigmaExpr`](../../src/R_interface_bartcore.cpp), [`ordinalThresholdsExpr`](../../src/R_interface_bartcore.cpp), filled from [`Results::sigma`](../../src/bartcore/chain.hpp), [`Results::ordinalThresholds`](../../src/bartcore/chain.hpp) | saved sample | 8 | S*C each; (K-1)*S*C for ordinal thresholds | per family |
| result variance channel | [`varianceTrainExpr`](../../src/R_interface_bartcore.cpp), [`varianceTestExpr`](../../src/R_interface_bartcore.cpp) | saved sample | 8 | n*S*C, plus nTest*S*C | heteroscedastic |
| yhat.train transient copy | [`convertSamplesFromDbartsToBart`](../../R/bart.R), [`packageBartResults`](../../R/bart.R) | saved sample | 8 | 1 extra n*L*S*C live at peak, on every family and either setting of combineChains | keepTrainingFits |
| yhat.test transient copy | [`convertSamplesFromDbartsToBart`](../../R/bart.R), [`packageBartResults`](../../R/bart.R) | saved sample | 8 | 1 extra nTest*L*S*C live at peak, on the same arithmetic | a test set |
| ordinal probability array | [`probsTrain`](../../R/bart.R) | saved sample | 8 | n*K*S*C, beside the n*S*C latent channel | ordinal, built R-side |
| raw predictors, caller's | the matrix the caller passes | sampler | 8 | n*p | always, beside the store's 2*n*p codes |
| raw predictors, ingestion high-water | [`dbartsData`](../../R/data.R)'s subset copy | sampler | 8 | up to n*p, plus nTest*p | always: the subset copy is persistent and fires on a matrix with no subset. The transient complete-cases copy beside it fired unconditionally too until this arc, and now only when a row is actually dropped. Measured per host and shape, not derived - see below |
| starting-sigma linear model | [`estimateSigmaFromLinearModel`](../../R/utility.R), [`residualStandardError`](../../R/utility.R) | sampler | 8 | about 2*n*(p+1) transiently - the design matrix and the QR's own copy of it | no `sigest` given and the family estimates a residual sd (not binary) |
| leaf statistics cache | [`LinearGaussianLeaf`](../../src/bartcore/model.hpp)'s per-tree crossproduct cache | chain | 4 | n per tree per chain - one live partition, the draw's prune having released everything else - plus up to as much again in retained capacity; 1.1 levels measured | the linear leaf; the gp leaf keeps its own kernel cache, arena-indexed and unpruned for dead slots |
| fit-path warm-up | the R session itself | sampler | not a byte count | one-off: byte-compiling the fit closures and populating the S4 dispatch tables | the first fit of a session |
| ingestion transients | [`makeModelMatrixFromDataFrame`](../../R/data.R) | sampler | 8 | up to 2*n*p live at once - the model frame's columns and the numeric matrix built from them, before the store quantizes | the formula and data-frame doors |
| quantile collector | [`QuantileGrid`](../../src/bartcore/data.hpp)'s [`sortedUnique`](../../src/bartcore/data.hpp) | sampler | 8 | n reserved per column, one column at a time | usequants |
| retained sampler | [`keepSampler`](../../R/bart.R) | sampler | the whole engine total above | 1 | keepTrees or keepSampler |

The engine writes recorded draws straight into the R vectors the bridge
allocates, so no engine-side copy of any result channel exists. Packaging
holds exactly two full-size arrays of a prediction channel at once, on every
family and either setting of combineChains: the engine's own, and the
permuted one the fit returns. It reached three until the two transient copies
came out (below): `matrix()` then `t()` built the combined layout in two
steps, and `apply` permuted the whole channel again to take the posterior
mean. Getting below two would need the bridge to allocate the channel in the
returned layout and the engine to write into it draw-major, which is a
strided write on the draw path and no part of this note.

Which INSTANT a total describes matters, and the measurement settled which
one wins. Two instants compete: ingestion, where the predictor copies and the
starting-sigma linear model are live, and packaging, where the result channels
and their copies are. On the matrix door a short run (n.burn = 0, small
n.samples, as the audit's cells are) peaks at PACKAGING, not at ingestion: the
per-chain n*T pair is already resident by then and dominates both. Ingestion
wins only when no `sigest` is given, and then by the starting-sigma least
squares, not by the store: on this host that estimate raises the peak by
18.2 MB at p = 10, 16.2 MB at p = 20 and 101.3 MB at p = 50 (n = 1e5,
T = 200, S = 10, each the paired difference between the same call with and
without `sigest`). It was 48.0, 70.8 and 178.0 MB while the estimate went
through `lm` (below). The formula and data-frame doors add
their own model-frame transients on top; the audit measures the matrix door,
where they do not exist.

Two rows above are not byte counts and are measured per host rather than
derived. The script reports each in its own column beside the closed form, so
the gate never hides behind them.

- The fit-path warm-up, 5.8 MB on this host: the R session's one-off growth
  from byte-compiling the fit closures and populating the S4 dispatch tables.
- The ingestion high-water of the predictor copies.
  [`dbartsData`](../../R/data.R) made two beyond the caller's matrix - a
  persistent subset copy and a transient complete-cases copy, both taken even
  when nothing was subset and nothing was missing - and the derived 16*n*p
  was an UPPER bound on what reached the peak, since the transient's pages
  are the collector's to reclaim before packaging. Measured then, over and
  above the caller's own matrix and response: 0.6 MB at n = 1e4 p = 10,
  8.2 MB at n = 1e4 p = 50, 22.4 MB at n = 1e5 p = 10, 38.3 MB at n = 1e5
  p = 20, 83.5 MB at n = 1e5 p = 50 and 376.5 MB at n = 1e6 p = 20 - between
  0.8 and 2.8 copies, never the derived 2 flat. The transient copy is now
  taken only when a row is actually dropped, so those figures are an upper
  bound on the current allowance, and the derived bound is 8*n*p; they have
  not been re-measured over the grid (below).

The measured p slope of the fit's own peak carries the same spread: 14.9
bytes per n*p between p = 10 and p = 20 at n = 1e5, 30.0 between p = 20 and
p = 50 there, 18.8 to 28.5 between p = 10 and p = 50 at n = 1e4. Two copies
plus the codes would be 18 flat and three plus codes 26; neither holds
everywhere, because whether the transient copy is still resident at packaging
is a collector outcome and not a byte count. Each allowance is the high-water
of its OWN instant, and the model sums them, which over-counts wherever a
later instant recycles an earlier one's pages: that is what the -5 to -7 MB
residuals at n = 1e4, p = 50 and the -24 MB at the two-chain cell are.

## The non-derived inputs

m, the mean live node count of a tree, is the first quantity not fixed by the
source (the leaf statistics cache's level count, below, is the second). It
enters the live-tree rows (negligible: 0.09 MB at T = 200, m = 8) and the
saved-tree and stored-state rows (where it is the whole term). The
current estimate is 3.8 at n = 2e3, growing slowly with n; the reference cases
below assume 8 at n >= 1e5. The audit measures it by summing the `tree.sizes`
blocks [`storeFlatTrees`](../../src/R_interface_bartcore.cpp) writes and
dividing by C*T: 5.73 at n = 2e4, T = 75, C = 2. That is a SHORT-RUN count -
the audit's cells run n.burn = 0, so the trees never reach their stationary
size - and so it neither confirms nor moves the reference cases' 8; it is a
lower bound on it. Two node counts exist and differ: the flattened
live count m that storeState and saved trees write, and the arena length
[`Tree::nodes`](../../src/bartcore/tree.hpp) holds. The arena only grows -
a death recycles its pair through [`freePairs`](../../src/bartcore/tree.hpp)
rather than shrinking the vector - so it is a high-water mark at or above the
live count, and 56*T*m is a LOWER bound on the live-tree row, not an upper one.
Step 2 can measure only the flattened count; the arena high-water is not
observable from R, and no channel reports it.

The leaf statistics cache's level count was the second, and is no longer one:
the cache is pruned, and what stays resident is one partition per tree. The
mechanism, not a node count: a cached entry is one leaf's ordered member list
plus an inline 81-double crossproduct, and the cache is indexed by ARENA slot
([`TreeStatisticsCache`](../../src/bartcore/model.hpp)). Leaf memberships
partition the observations, so the lists live in one tree sum to at most 4*n -
and that, [`statisticsEntryBytes`](../../src/bartcore/model.hpp) over
`members.size()`, is the only thing the 256 MiB budget counts.

What the audit found resident was three times that, for two reasons the budget
did not see. `assign` leaves each slot's member vector at the capacity of the
largest membership that slot ever held, and nothing released a slot when its
node stopped being a leaf, so a run that started from stumps parked 4*n in the
root's slot and 4*n more across each depth level below it; separately, each
populated slot carries 648 bytes of inline crossproduct plus its vector header
whatever q is. The first dominates and the second is noise: instrumented at
n = 1e5, T = 200, C = 1 on the linear leaf, the member lists held 248.3 MB
resident against 178.1 MB counted, over 1097 populated slots whose inline
crossproducts came to 0.8 MB - 99.7 pct of the cache was member lists, live
and stale together, and 0.3 pct was the inline entry. An honest budget alone
would not have helped: it would have started refusing at a ceiling the cache
was already near, trading the megabytes for rescans.

Both are now bounded at the source.
[`drawFromPosteriorForNode`](../../src/bartcore/model.hpp) runs on the settled
tree - a rejected move has already rolled back - and releases the slots that
are not live leaves, which is what the interior nodes an accepted grow leaves
behind and the pairs an accepted death frees have in common;
[`storeCrossproduct`](../../src/bartcore/model.hpp) reallocates rather than
reuse a capacity that has run past twice its membership. Neither can move a
value: every lookup re-validates its member list, so a released entry costs a
rescan and a rescan is bitwise what the entry would have served. What remains
resident is the live partition, 4*n per tree per chain, plus at most as much
again in retained capacity and `sizeof(CachedNodeStatistics)` per arena slot
ever touched. The same instrumented cell reads 86.1 MB of member lists against
80.0 MB counted (exactly 4*n*T) over 447 populated slots, so the model's
multiplier is 1.1 levels rather than 3.2, and process peak RSS for that fit
fell from 689.0 MB to 518.2 MB, medians of nine and six runs. The six-cell
excursion behind the old 3.2 (n in {1e4, 1e5, 2e5}, T in {75, 200}, C in
{1, 2}) has not been re-run, and neither has the grid; the multiplier
[`memory-footprint.R`](../../benchmarks/R/memory-footprint.R) carries is the
post-prune measurement at the one cell. The budget still never bound at any
cell measured, and now bounds residency within the factor above when it does.

## Reference cases

MB is 1e6 bytes. Both are `bart` calls on the defaults except where named:
S = 500, C = 4 (case 1) are defaults, T = 200 is set. Gaussian, constant leaf,
no test set, keepTrees FALSE, dense matrix, m = 8.

Case 1: n = 1e5, p = 20, T = 200, C = 4, S = 500.

| component | scope | MB |
| --- | --- | --- |
| cut codes 2*n*p, cut grid 8*100*p | sampler | 4.02 |
| owned conditioning vectors 24*n | sampler | 2.40 |
| index buffer 4*n*T | chain | 80.00 |
| leaf map 4*n*T | chain | 80.00 |
| total fits, residual, rescaled response 24*n | chain | 2.40 |
| move scratch 4*n | chain | 0.40 |
| live trees, leaf values, Tree objects | chain | 0.17 |
| per chain | | 162.97 |
| engine total, 4.02 + 2.40 + 4*162.97 | | 658.3 |
| x, two live copies 16*n*p, y 8*n, sigma, varcount | R | 33.0 |
| yhat.train, two copies live at peak, 16*n*S*C | R | 3200.0 |
| peak, of which the engine is 17 pct | | 3891.3 |

Case 2: n = 1e6, p = 50, C = 1, everything else as above.

| component | scope | MB |
| --- | --- | --- |
| cut codes and cut grid | sampler | 100.04 |
| owned conditioning vectors 24*n | sampler | 24.00 |
| index buffer 4*n*T | chain | 800.00 |
| leaf map 4*n*T | chain | 800.00 |
| total fits, residual, rescaled response 24*n | chain | 24.00 |
| move scratch 4*n | chain | 4.00 |
| live trees, leaf values, Tree objects | chain | 0.17 |
| per chain | | 1628.17 |
| engine total | | 1752.2 |
| x, two live copies 16*n*p, y, sigma, varcount | R | 808.1 |
| yhat.train, two copies live at peak, 16*n*S | R | 8000.0 |
| peak, of which the engine is again 17 pct | | 10560.3 |

Both cases were measured on this host, each as a single `bart` call under
`/usr/bin/time -l` with n.burn = 0 and `sigest` supplied, once against the
library built before the packaging copies came out and once after. Case 1:
5847.0 MB before, 4234.8 MB after. Case 2: 15172.7 MB before, 11161.5 MB
after. Both readings predate the ingestion guard above, which takes a
further 8*n*p off each - 16.0 MB in case 1 and 400.0 MB in case 2 - and
neither has been re-measured since. Each drop is one whole prediction array
to within a megabyte - 1612.2 MB against a derived 1600.0, and 4011.2 MB
against 4000.0 - which is
what a copy count, and not a coefficient, predicts. Both peaks read 4 to 6 pct
above the derived totals above, the same direction and size as the grid's own
residuals.

What the consuming arcs read off this: the three owned conditioning vectors
plus the test offset (dec-B87, [docs/decisions.md](../decisions.md)) are
24*n + 8*nTest, 2.4 MB in case 1 (0.4 pct of the engine) and 24.0 MB in case 2
(1.4 pct); one more chain costs the per-chain line, 163.0 MB and 1628.2 MB,
keepTrees costs 24*S*T*m + 40*S*T = 23.2 MB per chain in both. The cost of one
more thread inside a chain is a FORWARD claim this note does not make:
within-chain threading is not in the tree, and dec-B89's arc measures it.

## Where this disagrees with the plan's Context

- The four owned conditioning vectors are NOT confined to the data-handle path
  today: [`createHolder`](../../src/R_interface_bartcore.cpp) sizes all four
  unconditionally on the ordinary creation route, so they are already in the
  baseline rather than added by dec-B87, and the total is 3*8*n + 8*nTest,
  not 4*8*n. The reference-case engine totals rise by that amount.
- The per-chain block is 163.0 MB, not 162.5: the move scratch's high-water
  4*n and the per-tree Tree/muByTree objects are real and were unlisted.
- dec-B104's "about 2.8 times the training predictions" is FIXTURE-SPECIFIC,
  and the disagreement is with that decision's wording alone (the plan's step 4
  already carries the caveat). What storeState writes under keepTrees is
  13*T*m per chain per saved draw against an 8*n*S*C prediction array, so the
  ratio is 13*T*m/(8*n) = 2600/n at T = 200 and m = 8: it crosses 2.8 near
  n = 930 and falls like 1/n, reaching 0.03 in case 1. Its provenance is the
  auto-store that fired only for keepTrees fits, gated on keepTrees, so the
  ratio was measured on a review fixture's small n; the concatenated per-forest
  vectors storeState writes were already the layout then, so small n is the
  whole explanation and R list overhead is not part of it. That auto-store is
  gone, reverted by dec-B33 and confirmed by dec-B104. retired:
  [R/bart.R:144-150](https://github.com/vdorie/dbarts/blob/6c740a41e5bb7da90d2e58494e8f622e8131ce87/R/bart.R#L144-L150)
  The audit anchors the arithmetic: at n = 2e4, T = 75, m = 5.73, C = 2 and
  S = 100 one draw's stored state is 0.00035 of the run's whole prediction
  array (13*T*m/(8*n*S)) and 0.035 of a single draw's (13*T*m/(8*n)). Neither
  reaches 2.8 anywhere but at n in the hundreds, so the manual carries the
  per-draw form with its condition and not the decision's bare number.

## What the measurement moved

The run at 80ff32b4 changed five rows and added none that a byte count
could have been read off without it.

- The single raw-predictor row split in two. The caller's own matrix, 8*n*p,
  is derived; the copies [`dbartsData`](../../R/data.R) makes on top of it -
  the persistent subset copy and, when a row is dropped, a transient
  complete-cases copy (taken unconditionally until this arc) - are a
  measured allowance, because the derived bound is an upper bound on how
  much of them survives to the peak. The measured p slope of the fit's peak
  runs 14.9 to 30.0 bytes per n*p depending on the pair, straddling the 18 a
  two-copy model gives and the 26 a three-copy one does; the single-copy
  model the note started with gave 10 and was wrong everywhere.
- The starting-sigma linear model is a new row and, on a short run with no
  `sigest`, the largest single R-side term. It is a least-squares fit over
  the whole design, not a sampler allocation, which is why the audit's cells
  supply `sigest` and price it on a paired excursion instead. It was `lm`
  when the row was written, at about 4*n*(p+1) for the model frame, the na
  filter, the model matrix and the QR; it is now the QR alone, over a design
  matrix built directly (below).
- The leaf statistics cache is a new row and was the largest single term of
  a designated-covariate fit: 276 MB at n = 1e5, T = 200, C = 1, where the
  whole rest of the fit is 350 MB. Without it the linear cell missed by
  274 MB, the only cell that missed at all. The first formula tried, 4*n per
  cached NODE, was wrong and only looked right because it saturated the
  256 MiB budget; leaf memberships partition the observations, so the live
  lists are at most 4*n per tree, and what made the cache large was retained
  vector capacity across an arena-indexed store that was never pruned. The
  budget counts the live lists only and bound no cell measured. The prune
  above took the row to one partition per tree and the cell's peak RSS from
  689.0 MB to 518.2 MB.
- Gathered leaf raw doubled, 8*n*q to 16*n*q: the leaf keeps its own
  standardized copy beside the store's gather.
- Two rows are not byte counts at all and are measured per host: the
  fit-path warm-up and the ingestion high-water. Naming them, and reporting
  each in its own column beside the closed form, is the honest form; folding
  them into a per-unit coefficient would have hidden R-collector effects
  inside an allocation model. A third, the training-fit mean's collector
  churn, was carried while the mean went through `apply` and is gone with it
  (below).

Nothing in the engine's own rows moved. The n*T pair measured 8.02 to 8.06
bytes per n*T under a constant leaf across every cell, against the derived
8, and the per-chain and per-sampler rows carried the chain excursion to
within 2 pct.

## What the removal moved

The second ranked row landed inside this arc, and it moved three rows of the
model above.

- The transient-copy rows lost their extra copy and their conditional. The
  combined-chain reshape is one `aperm` whose flattening is a `dim<-`, and
  the posterior mean is reduced over the returned layout, where an
  observation's draws are one contiguous slab. Measured on the duplicates
  probe at a 32 MB channel, the reshape's heap high-water went from 2.00
  times the array to 1.01 and the mean's from 2.48 to 1.47 - each exactly one
  array apart.
- The order the mean sums in is what fixes where the reduction can live.
  `mean()` adds in the order it is handed and applies a second correction
  pass over the same order, so a reduction over the engine's
  observation-major array would change the last bits of every reported mean.
  Reducing after the permutation keeps them identical, which is what lets the
  removal be gated by comparing whole packaged fits.
- The collector-churn allowance is gone from the model, not because the R
  heap stopped churning - the duplicates probe still shows 1.47 times the
  array under the reduction against `apply`'s 2.48 - but because none of it
  reaches peak RSS in a fit any more. Scored with the allowance kept, the
  grid's median absolute relative residual was 8.6 pct against a 5 pct limit
  and every residual was negative; scored without it the median is 4.2 pct
  and the residuals fall either side of zero. Its former size, 16.4 MB at
  n = 1e5 and 121.2 MB at n = 1e6 (10 draws), was `apply`'s.

## The largest avoidable allocations, ranked

By bytes saved per line of change, at the two reference cases. The key is
UNCONDITIONAL value first: the two rows that pay at the reference cases as
they stand lead, and the rows below them are ordered by bytes per line within
their condition, which the recommendation column names. Everything below the
second row is its own TODO entry rather than a change in this arc.

| item | saved, case 1 | saved, case 2 | change | recommendation |
| --- | --- | --- | --- | --- |
| name `keepTrainingFits = FALSE` (legacy `keeptrainfits`) in the manual as the large-n lever | 3200 MB | 8000 MB | one sentence | taken, in the manual's Memory section; the figures are the two live copies, down from three |
| take the column means over the returned layout and build that layout in one permutation, so neither extra copy exists | 1600 MB | 4000 MB | the R reshape and mean | TAKEN, measured 1612.2 MB and 4011.2 MB; the last copy would need the bridge to allocate the channel draw-major and the engine to write into it strided, which is its own item |
| drop the transient complete-cases copy of the predictor matrix when nothing is missing | 16 MB | 400 MB | one branch in [`dbartsData`](../../R/data.R) | TAKEN, one guard in [`dbartsData`](../../R/data.R): the row selection runs only when a row is actually dropped. Measured at n = 1e5, p = 50, T = 200, S = 10 with `sigest` supplied, 410.0 MB to 367.6 MB - 42.4 MB against the derived 8*n*p of 40.0 |
| prune the leaf statistics cache, or bound it by resident bytes rather than by tracked member bytes | 0 (constant leaf) | 0 (constant leaf) | the store and the draw | TAKEN, as the prune: [`drawFromPosteriorForNode`](../../src/bartcore/model.hpp) releases the slots that are not live leaves and [`storeCrossproduct`](../../src/bartcore/model.hpp) caps retained capacity at twice the membership. CONDITIONAL on a designated-covariate leaf, where the cache was 4*n per arena level per tree per chain and is now one partition per tree. Measured at n = 1e5, T = 200, C = 1 on the linear leaf: peak RSS 689.0 MB to 518.2 MB, member lists 248.3 MB to 86.1 MB, the fit 1.8 pct slower on the same cell. The honest-budget alternative was declined: it would have refused entries at a ceiling the cache was already near instead of giving the bytes back |
| a cheaper starting sigma than an `lm` over the whole design | 0 today (packaging peaks higher) | 0 today | a few lines in [`estimateSigmaFromLinearModel`](../../R/utility.R) | TAKEN: [`residualStandardError`](../../R/utility.R) calls the QR routine `lm` calls, on the design matrix `model.matrix` would have built, and takes `summary.lm`'s own expression for sigma over it, so the estimate is bitwise unchanged with no model frame and no second design. Measured at n = 1e5, T = 200, S = 10 with no `sigest`: 348.9 to 318.5 MB at p = 10, 411.8 to 357.3 MB at p = 20, 588.5 to 511.2 MB at p = 50. Zero at either reference case, which supplies `sigest` |
| a flat arena for saved trees instead of a vector per tree (keepTrees only) | 4 MB/chain at keepTrees TRUE | 4 MB/chain at keepTrees TRUE | one engine struct | own TODO entry, post-release |
| drop the raw x when no mutation surface is in use | 16 MB | 400 MB | ingestion and predict both touched | not recommended; re-quantization needs it |
| leafOf as uint16 | 40 MB/chain | 400 MB/chain | declined by measurement | do not reopen |

What the consuming arcs take from this. Within-chain threading (dec-B89,
[docs/decisions.md](../decisions.md)) weighs the per-chain block, 163.0 MB in
case 1 and 1628.2 MB in case 2, against a per-thread block this note contains
no row for at all: no allocation here is indexed by thread, so at large n
threads inside a chain are the cheap parallelism and chains are the expensive
one. The header arc (dec-B87) takes its owned conditioning vectors as a share
of the per-sampler total: 24*n + 8*nTest is 2.4 MB of case 1's 658.3 MB engine
(0.4 pct) and 24.0 MB of case 2's 1752.2 MB (1.4 pct), and the measurement
leaves both unchanged - the rows the audit moved are all R-side or leaf-model
rows, none of them on the conditioning path.
