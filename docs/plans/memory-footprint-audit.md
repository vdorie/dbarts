# memory-footprint-audit

agent: opus (steps 1 and 5, the derivation and the ranked list); sonnet
  (steps 2 to 4 and 7, the script, the fit run, the manual and the
  two-copy removal). Step 1 lands the note first: the script's predicted
  column is that model; step 7 lands last, after the script has recorded
  the before row.
rng: neutral. No draw path is touched; the script is a benchmark, not a
  gate; step 7 changes how the result array is laid out and reduced, not
  any draw, so the three bitwise equivalence baselines (gaussian, BCF,
  multinomial) expect IDENTICAL and every packaged result is compared
  element-wise against the pre-change fit.
window: independent of the other pre-release arcs and wanted early: it sets
  the within-chain threading default (dec-B89) and prices the owned
  conditioning vectors the header arc adds ([pure-c-header.md](pure-c-header.md), dec-B87).
budget: ~150 lines benchmarks/R/memory-footprint.R, ~100 lines
  docs/design/memory-footprint.md, ~40 lines Rd; step 7 ~40 lines R (and
  bridge dimension changes if the implementer takes that route) plus
  ~60 lines of tests.

Decisions in [docs/decisions.md](../decisions.md): dec-B87 (this audit prices
the owned copies), dec-B89 (it informs the threading default), dec-B104 (the
stored state is about 2.8 times the training predictions).

## Goal

A closed-form footprint model, per component, in the five units the
allocations have: once per sampler, per chain, per tree, per saved sample,
per test row. It is validated against measured peak resident set size over a
grid in n, p, trees and chains, with and without keepTrees, with and without
a test set, for the constant leaf and one non-constant leaf, on gaussian and
probit; the R-layer duplicates are priced separately. The manual states it in
words with one worked table, and the arc closes with the largest avoidable
allocations ranked by bytes saved per line of change, each becoming its own
TODO entry rather than a change here.

## Reference cases, from the code

Computed from the sizes derived in Context. MB is 1e6 bytes. Both cases are
`bart` calls (front-door.md S1 renamed `bart2` to `bart`) on its defaults
except where named: S = 500 kept draws per
chain, and 4 chains in case 1, are defaults; 200 trees is set explicitly (75
is the default). Both assume a mean of 8 nodes per live tree (measured 3.8 at
n = 2e3; it grows slowly with n and step 1 measures it) and the dense-matrix
path, whose store owns codes but no raw doubles.

Case 1: n = 1e5, p = 20, 200 trees (T), 4 chains (C), 500 draws (S),
gaussian, constant leaf, no test set, keepTrees FALSE.

| component | scope | bytes/unit | units | MB |
| --- | --- | --- | --- | --- |
| train cut codes and cut grid | sampler | 2, 8 | n*p = 2.0e6, 100*p | 4.02 |
| index buffer | chain | 4 | n*T = 2.0e7 | 80.0 |
| leaf map | chain | 4 | n*T = 2.0e7 | 80.0 |
| total fits, residual store, rescaled response | chain | 24 | n | 2.4 |
| live trees | chain | 56 | T*nodes = 1600 | 0.09 |
| per chain | | | | 162.5 |
| engine total | | | 4.0 + 4*162.5 | 654.0 |
| x, y, sigma, varcount | R | | | 17.0 |
| yhat.train, three copies live at peak | R | 24 | n*S*C = 2.0e8 | 4800.0 |
| R total, at peak | | | | 4817.0 |
| peak, of which the engine is 12 pct | | | | 5471.0 |

Case 2: n = 1e6, p = 50, 200 trees, 1 chain, everything else as above.

| component | scope | bytes/unit | units | MB |
| --- | --- | --- | --- | --- |
| train cut codes and cut grid | sampler | 2, 8 | n*p = 5.0e7, 100*p | 100.04 |
| index buffer | chain | 4 | n*T = 2.0e8 | 800.0 |
| leaf map | chain | 4 | n*T = 2.0e8 | 800.0 |
| total fits, residual store, rescaled response | chain | 24 | n | 24.0 |
| live trees | chain | 56 | T*nodes = 1600 | 0.09 |
| engine total | | | 100.0 + 1624.1 | 1724.1 |
| x, y, sigma, varcount | R | | | 408.1 |
| yhat.train, three copies live at peak | R | 24 | n*S = 5.0e8 | 12000.0 |
| R total, at peak | | | | 12408.1 |
| peak, of which the engine is again 12 pct | | | | 14132.2 |

| what the consuming arcs read off this | case 1 | case 2 |
| --- | --- | --- |
| the header arc's four owned vectors, 24*n + 8*nTest | 2.4 MB, 0.4 pct of engine | 24.0 MB, 1.4 pct |
| one more chain | 162.5 MB | 1624.1 MB |
| one more thread inside a chain | under 1 MB | under 1 MB |
| keepTrees, per chain, 24*S*T*nodes plus per-tree overhead | 23.2 MB | 23.2 MB |

## Context

- Widths: [`xint_t`](../../src/bartcore/data.hpp) codes two bytes,
  [`index_t`](../../src/bartcore/data.hpp) indices four,
  [`Node`](../../src/bartcore/tree.hpp) 56 and
  [`FlatNode`](../../src/bartcore/tree.hpp) 24 by layout. Once per sampler: packed codes ([`CodeBlock::codes`](../../src/bartcore/data.hpp)
  under [`ColumnStore::train`](../../src/bartcore/data.hpp)) at 2*n*p, the
  same again at 2*nTest*p for a test set, the cut grid
  ([`ColumnStore::cutPoints`](../../src/bartcore/data.hpp)) at 8 bytes per
  cut, O(p) metadata, and a sparse column
  ([`SparseColumnData`](../../src/bartcore/data.hpp)) at 0.19*n + 2*nnz
  instead. [`ColumnStore::ownedDenseValues`](../../src/bartcore/data.hpp)
  is sized only on the container/CSC path and
  [`ColumnStore::gatheredRawValues`](../../src/bartcore/data.hpp) only under a
  designated-covariate leaf.
- Per chain, the dominant pair:
  [`Forest::indexBuffer`](../../src/bartcore/combiner.hpp) at 4*n*T, sliced
  per tree into [`Tree::indices`](../../src/bartcore/tree.hpp) and sized by
  the forest initializer in [`Chain`](../../src/bartcore/chain.hpp)'s
  constructor, and, under a constant leaf,
  [`Forest::leafOf`](../../src/bartcore/combiner.hpp) at 4*n*T from
  [`Chain::initForestFitStorage`](../../src/bartcore/chain.hpp); a
  non-constant leaf drops leafOf for
  [`Forest::treeFits`](../../src/bartcore/combiner.hpp) at 8*n*T, half again
  as much. Also [`Forest::totalFits`](../../src/bartcore/combiner.hpp) at
  8*n, [`Forest::treeY`](../../src/bartcore/combiner.hpp) at n times the
  residual width (4 under storage="single"), and the live trees.
- Per chain, per family: [`GaussianResponse::yRescaled_`](../../src/bartcore/model.hpp)
  (8*n, always) and [`ProbitResponse::latents_`](../../src/bartcore/model.hpp)
  (two 8*n); a row mask adds two more 8*n, heteroscedasticity a whole
  [`VarianceForest::indexBuffer`](../../src/bartcore/chain.hpp) family. Test
  slabs are nTest-length, not nTest*T -
  [`Forest::totalTestFits`](../../src/bartcore/combiner.hpp) and its current
  twin, 16*nTest per forest per chain, from
  [`Chain::resizeTestStorage`](../../src/bartcore/chain.hpp).
- Per saved sample: [`Forest::savedTrees`](../../src/bartcore/combiner.hpp),
  C*S*T vectors of FlatNode - 24*S*T*nodes per chain plus about 40 bytes per
  saved tree of header and allocator overhead, 4 MB at S*T = 1e5.
- The bridge allocates each recorded channel directly in R memory
  ([`allocChannel`](../../src/R_interface_bartcore.cpp),
  [`installChannel`](../../src/R_interface_bartcore.cpp)) and the engine
  writes into it, so no engine-side yhat.train exists. The R duplicates:
  [`dbartsData`](../../R/A_class.R) keeps raw x at 8*n*p beside the store's
  2*n*p codes;
  [`convertSamplesFromDbartsToBart`](../../R/bart.R) permutes n x S x C into
  C x S x n, allocating the destination before releasing the source, and
  [`packageBartResults`](../../R/bart.R) then takes the column means with
  `apply`, which allocates a third full-size array (it transposes its
  argument even when the permutation is the identity) while the engine's
  array is still bound - peak is three times yhat.train. It also retains the
  whole sampler under [`keepSampler`](../../R/bart.R), which keepTrees
  implies. [`BartcoreHolder::ownedResponse`](../../src/R_interface_bartcore_common.hpp)
  and its three siblings exist only on the data-handle path today; dec-B87
  makes them unconditional at 4*8*n per sampler.
- Standing records: reduced-precision-storage
  [6. What landed](../design/reduced-precision-storage.md#6-what-landed) is
  the only existing formula set (index narrowing, about 400 MB saved at
  n = 5e5; leafOf uint16 declined by measurement);
  [The adopted design](../design/data-ownership.md#the-adopted-design) and
  [data-store.md](../design/data-store.md) describe the store. Harness: [`genFriedman`](../../benchmarks/R/bench-sampler.R) and
  the record/compare CLI beside it, row schema scenario, metric, value, rev,
  date, quick; [benchmarks/baselines/MANIFEST](../../benchmarks/baselines/MANIFEST)
  (the provenance file giving each baseline's commit, machine and coverage)
  carries a new baseline's rules.

## Decision

Forks 1 and 2 were put to VD on 2026-09-08 and are recorded with the
choice; fork 3 is the plan's own call.

1. Whether the audit also runs on the x86 bench box (VD 2026-09-08,
   "Use your recommendation"): arm64 only. Every quantity above is a
   byte count fixed by the source; the host-dependent terms are
   allocator overhead and page granularity, which the baseline
   subtraction mostly removes. A residual above the step 2 tolerance
   with no code-side explanation returns to VD as a grant request for
   the second host.
2. Whether any reduction lands pre-release (VD 2026-09-08, "Use your
   recommendation" on adding the two-copy removal): the audit lands its
   model and manual section, and the removal of the two transient copies
   of the training-prediction array joins this arc as step 7; every
   other ranked item stands on its own afterwards. The regression rule
   binds nothing here, verified rather than assumed: 0.9-34's classic
   engine allocated per chain a size_t index array and a double
   tree-fits array at n*trees each, 16 bytes per observation-tree pair
   against today's 8 under a constant leaf. retired:
   [src/dbarts/state.cpp:51-60](https://github.com/vdorie/dbarts/blob/edcdf735855ad331f785c3577f0305d7b2dc224f/src/dbarts/state.cpp#L51-L60),
   in the since-deleted classic engine.
3. Grid size against run time. Recommended: a spine, not the full crossing
   of the eight axes (384 cells, several hours). The model is additive and
   separable, so a base cell plus one-axis-at-a-time excursions (about 40
   cells) identifies every coefficient, with a small full crossing at the
   smallest n as the interaction check. Each cell runs n.burn = 0 and
   n.samples = 10, since only saved trees and the R channels track S; one
   cell at S = 200 checks those. Under ten minutes.

## Constraints

- Gates: the script runs to completion and the residual meets its stated
  tolerance. No engine or C++ change, so no tinytest, equivalence or
  sanitizer run is owed; but the new script is R source, so `air format
  --check .` and `lintr::lint()` cover it, and the step 4 Rd edits need
  `R CMD check --as-cran` from a tarball plus `tools/check-rc-codoc.R` for
  the reference-class methods `codoc` cannot see. Maintainer-run on a quiet
  machine: a peak RSS beside other load is not evidence.
- Out of scope: any change to what the engine allocates (the leafOf
  narrowing included, declined by measurement); the R result layout; the
  threading default, which is dec-B89's arc reading these numbers.

## Steps

1. Derive the model into a new design note, `docs/design/memory-footprint.md`,
   as one table over component, symbol, scope (sampler, chain, tree, saved
   sample, test row), bytes per unit, unit count and the condition under
   which the component exists; every row cites its symbol. Rows are the
   Context list plus the variance forest, the per-category multinomial and
   ordinal channels, and the R objects. It carries the two reference cases and
   its one non-derived input, the mean live node count. Status OPEN until step
   2 validates it.
2. Write `benchmarks/R/memory-footprint.R`: one subprocess per cell under
   `/usr/bin/time -l` on macOS and `-v` on Linux, parsing maximum resident
   set size (bytes on macOS, kilobytes on Linux) and subtracting a baseline
   measured the same way from a subprocess that loads the package and fits
   nothing. Grid: base cell n = 1e5, p = 20, 200 trees, 1 chain, gaussian,
   constant leaf, no test set, keepTrees FALSE; excursions in n over 1e4,
   1e5, 1e6; p over 10, 20, 50; trees over 75, 200; chains over 1, 2, 4;
   keepTrees TRUE; a test set at nTest = n/5; [`linear`](../../R/model.R) as
   `node.prior` over three columns; family probit; the S = 200 cell. Output
   is a CSV in the bench-sampler schema, one row per cell per metric
   (peak_rss_mb, predicted_mb, residual_mb), plus a MANIFEST row naming the
   host and the baseline subtracted. It then fits measured against predicted,
   to a tolerance of every cell within 10 percent or 20 MB, whichever is
   larger, and a median absolute relative residual under 5 percent; a miss is
   a missing term, a defect in the note. Flip its Status to LANDED.
3. Price the R-layer duplicates in the same script: `gc(reset = TRUE)` then
   `gc()`'s max-used column around each reshape and the `apply` mean, with
   `object.size` on x, y, on yhat.train before and after
   [`convertSamplesFromDbartsToBart`](../../R/bart.R), and on the sampler
   retained under `keepSampler`, and a check of whether the combined-chain
   reshape adds a fourth full-size copy. Their own CSV rows.
4. Manual: `man/dbarts-package.Rd` gains `\section{Memory}` beside its
   existing sections - the model in words (two arrays of 4 bytes per
   observation per tree per chain dominate, so it is linear in n, trees and
   chains; the predictor store is 2 bytes per cell once, whatever the chain
   count; a test set costs its own store plus 16 bytes per test row per
   chain), a table of worked examples at the grid corners, and that a
   large-n fit is usually limited by the returned yhat.train array, not by
   the sampler, naming the lever `keepTrainingFits` (`bart` and the control
   object; `keeptrainfits` on the legacy door, `bartBT`).
   [`dbartsControl`](../../man/dbartsControl.Rd)'s `keepTrees` entry gains a
   sentence pricing saved trees;
   [`dbartsSampler$storeState`](../../man/dbartsSampler-class.Rd)'s entry one
   with the stored-state ratio. Caution: dec-B104's 2.8x is not scale-free,
   and which ratio it means must be settled first. Saved trees over all S
   draws against the S-draw prediction array is 24*S*T*nodes / (8*n*S) =
   3*T*nodes/n, which is 2.8 at n near 1700 and falls like 1/n (0.05 at
   n = 1e5); what `storeState` writes is one draw's state, 24*T*nodes per
   chain, so against the same array it is 3*T*nodes/(n*S) and reaches 2.8 at
   no plausible cell. Step 2 measures both; the Rd carries whichever the
   measurement supports, with its condition.
5. The ranked list, in the design note, by bytes saved per line of change.
   Expected order, from the model:

   | item | saved, case 1 | saved, case 2 | change | recommendation |
   | --- | --- | --- | --- | --- |
   | name `keepTrainingFits = FALSE` (legacy `keeptrainfits`) in the manual as the large-n lever | 4800 MB | 12000 MB | one sentence | take it, in step 4 |
   | have the bridge allocate the result in the layout the R side returns, and take the column means over it, so neither extra copy exists | 3200 MB | 8000 MB | bridge dims plus the R reshape and mean | step 7 of this arc (VD 2026-09-08) |
   | a flat arena for saved trees instead of a vector per tree (keepTrees only; zero in both cases, which run keepTrees FALSE) | 4 MB/chain at keepTrees TRUE | 4 MB/chain at keepTrees TRUE | one engine struct | own TODO entry, post-release |
   | drop the raw x when no mutation surface is in use | 16 MB | 400 MB | ingestion and predict both touched | not recommended; re-quantization needs it |
   | leafOf as uint16 | 40 MB/chain | 400 MB/chain | declined by measurement | do not reopen |

   Then what the consuming arcs take: within-chain threading (dec-B89) gets
   the per-chain block against the near-zero per-thread block, making threads
   the cheap parallelism at large n; the header arc (dec-B87) its owned
   copies as a share of the per-sampler total.
6. Records: close the TODO entry, naming the note and the script; add the
   CSV and its MANIFEST row; one NEWS line under 1.0-0. The feature matrix is
   untouched.
7. Remove the two transient copies of the training-prediction array in
   [`packageBartResults`](../../R/bart.R) and
   [`convertSamplesFromDbartsToBart`](../../R/bart.R): the fitted means are
   taken over the engine's own layout before any permutation (a
   `colMeans`-style reduction over the observation-major array, not
   `apply`), and the array reaches its returned layout either by the bridge
   allocating it in that layout or by one permutation into a
   preallocated array, so no moment holds more than one full-size copy
   beyond the one returned. Same treatment for the test predictions.
   Every returned element is identical to the pre-change fit (a test
   compares whole result objects at one seed for gaussian, probit and one
   multi-chain uncombined case); the script's before and after rows show
   the peak drop, expected from about 5.5 GB to about 2.2 GB in case 1 and
   14 GB to 6 GB in case 2.

## Verification

```
Rscript benchmarks/R/memory-footprint.R record benchmarks/baselines/memory-footprint-<rev>.csv
Rscript benchmarks/R/memory-footprint.R fit benchmarks/baselines/memory-footprint-<rev>.csv
air format --check . && Rscript -e 'lintr::lint("benchmarks/R/memory-footprint.R")'
Rscript tools/check-doc-freshness.R && Rscript tools/check-rc-codoc.R
R CMD build . && R CMD check --as-cran dbarts_*.tar.gz
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'                      # step 7
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-c42b72af.rds   # step 7: 50 identical
```

Expected: record writes one row per cell per metric; fit prints the residual
summary and exits non-zero on a missed step 2 tolerance; format, lint,
freshness, codoc and check clean.
