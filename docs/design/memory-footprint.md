# Memory footprint

Status: OPEN until step 2 of [memory-footprint-audit.md](../plans/memory-footprint-audit.md)
validates it against measured peak RSS.

The closed-form model of what a fit allocates, per component, in the units the
allocations actually have. Every byte count below is read off the element type
and the allocation site, not measured. Symbols: n rows, p predictors, T trees
per forest, C chains, S kept draws, nTest test rows, K categories, q designated
leaf covariates, nnz sparse nonzeros, and m the mean live node count of a tree
(the one non-derived input, below). Widths that recur:
[`xint_t`](../../src/bartcore/data.hpp) 2, [`index_t`](../../src/bartcore/data.hpp)
4, [`Node`](../../src/bartcore/tree.hpp) 56, [`FlatNode`](../../src/bartcore/tree.hpp)
24 - each confirmed by `sizeof` on arm64, and none of them changes with the
reference build, which re-associates draw-path kernels and allocates nothing
differently. No row below differs by build mode.

## The table

| component | symbol | scope | bytes/unit | units | exists when |
| --- | --- | --- | --- | --- | --- |
| train cut codes | [`CodeBlock::codes`](../../src/bartcore/data.hpp) under [`ColumnStore::train`](../../src/bartcore/data.hpp) | sampler | 2 | n*p | dense-stored columns (always, off CSC rank storage) |
| test cut codes | [`ColumnStore::test`](../../src/bartcore/data.hpp) | test row | 2 | nTest*p | a test set |
| cut grid | [`ColumnStore::cutPoints`](../../src/bartcore/data.hpp) | sampler | 8 | up to n.cuts (default 100) per numeric column | always |
| column metadata | [`ColumnStore::numCuts`](../../src/bartcore/data.hpp), [`ColumnStore::categoryCounts`](../../src/bartcore/data.hpp), [`ColumnSource`](../../src/bartcore/data.hpp) | sampler | ~70 per column per row set | p | always |
| sparse column | [`SparseColumnData`](../../src/bartcore/data.hpp) | sampler | 0.1875*n + 2*nnz | per rank-stored column, replacing its 2*n codes | a CSC column at or below 20 pct density |
| owned dense raw | [`ColumnStore::ownedDenseValues`](../../src/bartcore/data.hpp), [`ColumnStore::ownedTestValues`](../../src/bartcore/data.hpp) | sampler | 8 | n (nTest) per real dense-backed column | CSC/mixed build; test side on any test build |
| gathered leaf raw | [`ColumnStore::gatheredRawValues`](../../src/bartcore/data.hpp) | sampler | 8 | n*q | a designated-covariate leaf (linear, gp) |
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
| result test channel | same | saved sample | 8 | nTest*L*S*C | a test set |
| result varcount | same | saved sample | 4 | p*F*S*C, F the forest axis | always |
| result sigma, k, thresholds, dispersion | same | saved sample | 8 | S*C each; (K-1)*S*C for ordinal thresholds | per family |
| result variance channel | same | saved sample | 8 | n*S*C, plus nTest*S*C | heteroscedastic |
| yhat.train transient copies | [`convertSamplesFromDbartsToBart`](../../R/bart.R), [`packageBartResults`](../../R/bart.R) | saved sample | 8 | 2 extra n*L*S*C live at peak | keepTrainingFits, non-binary |
| ordinal probability array | [`probsTrain`](../../R/bart.R) | saved sample | 8 | n*K*S*C, beside the n*S*C latent channel | ordinal, built R-side |
| raw predictors | [`dbartsData`](../../R/A_class.R) | sampler | 8 | n*p, plus nTest*p | always, beside the store's 2*n*p codes |
| retained sampler | [`keepSampler`](../../R/bart.R) | sampler | the whole engine total above | 1 | keepTrees or keepSampler |

The engine writes recorded draws straight into the R vectors the bridge
allocates, so there is no engine-side copy of any result channel.

## The one non-derived input

m, the mean live node count of a tree, is the only quantity not fixed by the
source. It enters the live-tree rows (negligible: 0.09 MB at T = 200, m = 8)
and the saved-tree and stored-state rows (where it is the whole term). The
current estimate is 3.8 at n = 2e3, growing slowly with n; the reference cases
below assume 8 at n >= 1e5. Step 2 measures it per grid cell by summing the
`tree.sizes` blocks [`storeFlatTrees`](../../src/R_interface_bartcore.cpp)
writes and dividing by C*T. Two node counts exist and differ: the flattened
live count m that storeState and saved trees write, and the arena size
[`Tree::nodes`](../../src/bartcore/tree.hpp) holds, a high-water mark that
retains dead pairs. The 56-byte row uses the arena, so it is an upper bound of
the same order.

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
| x 8*n*p, y 8*n, sigma, varcount | R | 17.0 |
| yhat.train, three copies live at peak, 24*n*S*C | R | 4800.0 |
| peak, of which the engine is 12 pct | | 5475.3 |

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
| x, y, sigma, varcount | R | 408.1 |
| yhat.train, three copies live at peak, 24*n*S | R | 12000.0 |
| peak, of which the engine is again 12 pct | | 14160.3 |

What the consuming arcs read off this: the three owned conditioning vectors
plus the test offset (dec-B87, [docs/decisions.md](../decisions.md)) are
24*n + 8*nTest, 2.4 MB in case 1 (0.4 pct of the engine) and 24.0 MB in case 2
(1.4 pct); one more chain costs the per-chain line, 163.0 MB and 1628.2 MB,
while one more thread inside a chain costs under 1 MB (dec-B89); keepTrees
costs 24*S*T*m + 40*S*T = 23.2 MB per chain in both.

## Where this disagrees with the plan's Context

- The four owned conditioning vectors are NOT confined to the data-handle path
  today: [`createHolder`](../../src/R_interface_bartcore.cpp) sizes all four
  unconditionally on the ordinary creation route, so they are already in the
  baseline rather than added by dec-B87, and the total is 3*8*n + 8*nTest,
  not 4*8*n. The reference-case engine totals rise by that amount.
- The per-chain block is 163.0 MB, not 162.5: the move scratch's high-water
  4*n and the per-tree Tree/muByTree objects are real and were unlisted.
- dec-B104's "about 2.8 times the training predictions" holds for no
  configuration of what storeState writes. The stored live state is 13*T*m per
  chain against an 8*n*S*C prediction array: 5e-5 of it in case 1. Saved trees
  over all S draws reach 2.8 only near n = 1700 and fall like 1/n. Step 2
  measures both ratios; the manual carries whichever survives, with its
  condition.
