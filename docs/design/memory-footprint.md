# Memory footprint

Status: VALIDATED against measured peak resident set size by
[benchmarks/R/memory-footprint.R](../../benchmarks/R/memory-footprint.R) at
80ff32b4 on arm64 macOS. Every one of the 29 grid cells falls within
max(10 pct, 20 MB) of the model and the median absolute relative residual is
2.0 pct over the twelve cells predicting more than 100 MB; the worst cell
(n = 1e5, p = 10, T = 200) is 13.5 MB high against a 24.0 MB tolerance. The
relative half of the tolerance is scored only above 100 MB: below it the
residual is set by page-level allocator behaviour rather than by any row
here, so a relative criterion there measures the host. "What the measurement
moved" below records the rows the run changed.

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
| yhat.train transient copies | [`convertSamplesFromDbartsToBart`](../../R/bart.R), [`packageBartResults`](../../R/bart.R) | saved sample | 8 | 2 extra n*L*S*C live at peak; 1 under combineChains = FALSE on a binary fit | keepTrainingFits |
| yhat.test transient copies | [`convertSamplesFromDbartsToBart`](../../R/bart.R), [`packageBartResults`](../../R/bart.R) | saved sample | 8 | 2 extra nTest*L*S*C live at peak, on the same arithmetic | a test set |
| ordinal probability array | [`probsTrain`](../../R/bart.R) | saved sample | 8 | n*K*S*C, beside the n*S*C latent channel | ordinal, built R-side |
| raw predictors, three live copies | [`dbartsData`](../../R/A_class.R) | sampler | 24 | n*p, plus nTest*p | always, beside the store's 2*n*p codes: the caller's matrix, the one the data object keeps, and the ingestion copy, all resident at the peak |
| starting-sigma linear model | [`estimateSigmaFromLinearModel`](../../R/utility.R) | sampler | 8 | about 4*n*(p+1) transiently, at the `summary.lm` instant - base R's model frame, na filter, model matrix and QR | no `sigest` given and the family estimates a residual sd (not binary) |
| leaf statistics cache | [`LinearGaussianLeaf`](../../src/bartcore/model.hpp)'s per-tree crossproduct cache | chain | 4 | n per cached node, over every chain, capped at the leaf's 256 MiB total budget | a designated-covariate leaf (linear, gp) |
| training-fit mean churn | [`packageBartResults`](../../R/bart.R)'s `apply` | saved sample | not a byte count | two small R objects per observation, resident until the collector's next cycle | keepTrainingFits on a non-binary fit |
| fit-path warm-up | the R session itself | sampler | not a byte count | one-off: byte-compiling the fit closures and populating the S4 dispatch tables | the first fit of a session |
| ingestion transients | [`makeModelMatrixFromDataFrame`](../../R/data.R) | sampler | 8 | up to 2*n*p live at once - the model frame's columns and the numeric matrix built from them, before the store quantizes | the formula and data-frame doors |
| quantile collector | [`QuantileGrid`](../../src/bartcore/data.hpp)'s [`sortedUnique`](../../src/bartcore/data.hpp) | sampler | 8 | n reserved per column, one column at a time | usequants |
| retained sampler | [`keepSampler`](../../R/bart.R) | sampler | the whole engine total above | 1 | keepTrees or keepSampler |

The engine writes recorded draws straight into the R vectors the bridge
allocates, so no engine-side copy of any result channel exists. The transient
copies are not a non-binary phenomenon: `bart`'s default combineChains = TRUE
reaches three live full-size arrays on every fit, since `matrix()` and `t()`
each allocate while the engine's array is still bound; the non-binary path
reaches the same three through `apply`'s `aperm`, which copies even when the
permutation is the identity. Only a binary fit under combineChains = FALSE
peaks at two.

Which INSTANT a total describes matters, and the measurement settled which
one wins. Two instants compete: ingestion, where the predictor copies and the
starting-sigma linear model are live, and packaging, where the result channels
and their copies are. On the matrix door a short run (n.burn = 0, small
n.samples, as the audit's cells are) peaks at PACKAGING, not at ingestion: the
per-chain n*T pair is already resident by then and dominates both. Ingestion
wins only when no `sigest` is given, and then by base R's `lm`, not by the
store: on this host that estimate raised the peak by 21.0 MB at p = 10,
58.3 MB at p = 20 and 147.7 MB at p = 50 (n = 1e5, T = 200, S = 10), about
30 further bytes per n*p as p grows. The formula and data-frame doors add
their own model-frame transients on top; the audit measures the matrix door,
where they do not exist.

Two rows above are not byte counts and are measured per host rather than
derived: the fit-path warm-up (6.5 MB on this host) and the training-fit
mean's collector churn. `apply` over the observation margin allocates two
small R objects per observation and R's collector leaves them resident until
its next cycle, so the peak carries 2.8 MB at n = 1e4, 16.2 MB at n = 1e5 and
120.4 MB at n = 1e6 (10 draws) that no byte count predicts. Both vanish from
the model the moment the mean is taken with a colMeans-style reduction
instead.

## The one non-derived input

m, the mean live node count of a tree, is the only quantity not fixed by the
source. It enters the live-tree rows (negligible: 0.09 MB at T = 200, m = 8)
and the saved-tree and stored-state rows (where it is the whole term). The
current estimate is 3.8 at n = 2e3, growing slowly with n; the reference cases
below assume 8 at n >= 1e5. The audit measures it by summing the `tree.sizes`
blocks [`storeFlatTrees`](../../src/R_interface_bartcore.cpp) writes and
dividing by C*T: 5.73 at n = 2e4, T = 75, C = 2. That is a SHORT-RUN count -
the audit's cells run n.burn = 0, so the trees never reach their stationary
size - and so it neither confirms nor moves the reference cases' 8; it is a
lower bound on it, and it is the value the leaf statistics cache row was
checked against. Two node counts exist and differ: the flattened
live count m that storeState and saved trees write, and the arena length
[`Tree::nodes`](../../src/bartcore/tree.hpp) holds. The arena only grows -
a death recycles its pair through [`freePairs`](../../src/bartcore/tree.hpp)
rather than shrinking the vector - so it is a high-water mark at or above the
live count, and 56*T*m is a LOWER bound on the live-tree row, not an upper one.
Step 2 can measure only the flattened count; the arena high-water is not
observable from R, and no channel reports it.

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
| x, three live copies 24*n*p, y 8*n, sigma, varcount | R | 49.0 |
| yhat.train, three copies live at peak, 24*n*S*C | R | 4800.0 |
| peak, of which the engine is 12 pct | | 5507.3 |

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
| x, three live copies 24*n*p, y, sigma, varcount | R | 1208.1 |
| yhat.train, three copies live at peak, 24*n*S | R | 12000.0 |
| peak, of which the engine is again 12 pct | | 14960.3 |

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

- Raw predictors went from 8*n*p to 24*n*p. Three copies of the predictor
  matrix are resident at the peak, not one: the caller's, the one the data
  object keeps beside the codes, and a third the ingestion path makes. The
  p slope of the measured peak is 25 to 28 bytes per n*p against the 10 the
  single-copy model gave.
- The starting-sigma linear model is a new row and, on a short run with no
  `sigest`, the largest single R-side term. It is base R's `lm`, not a
  sampler allocation, which is why the audit's cells supply `sigest` and
  price it on a paired excursion instead.
- The leaf statistics cache is a new row and the largest single term of a
  designated-covariate fit: 268 MB at n = 1e5, T = 200, C = 1, where the
  whole rest of the fit is 350 MB. Without it the linear cell missed by
  274 MB, the only cell that missed at all; with it the same cell lands
  1.7 MB high. It is capped, so it does not scale past 256 MiB, but it
  reaches the cap at every cell above n*T*m*C = 6.7e7.
- Gathered leaf raw doubled, 8*n*q to 16*n*q: the leaf keeps its own
  standardized copy beside the store's gather.
- Two rows are not byte counts at all and are measured per host: the
  fit-path warm-up and the training-fit mean's collector churn. Naming them
  is the honest form; folding them into a per-unit coefficient would have
  hidden an R-collector effect inside an allocation model.

Nothing in the engine's own rows moved. The n*T pair measured 8.02 to 8.06
bytes per n*T under a constant leaf across every cell, against the derived
8, and the per-chain and per-sampler rows carried the chain excursion to
within 2 pct.
