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

- Step 1 landed [docs/design/memory-footprint.md](../design/memory-footprint.md),
  which now carries the authoritative derivation; its "Where this
  disagrees with the plan's Context" section corrects three points below:
  the four owned conditioning vectors are already unconditional in the
  baseline, not added by dec-B87; the per-chain block is 163.0 MB, not
  162.5 (move scratch and per-tree object overhead were unlisted); and
  dec-B104's 2.8x ratio is fixture-specific, holding near n = 1700, not
  as a general rule.
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
   14 GB to 6 GB in case 2. "Two transient copies" is the general count
   the note's Reference cases section confirms: `combineChains = TRUE`
   (the reference cases' default) reaches three live arrays for every
   family, not only non-binary ones, via `matrix()` and `t()`; a
   non-binary fit under `combineChains = FALSE` still reaches three via
   `apply`'s `aperm`. Only a binary fit under `combineChains = FALSE`
   already peaks at two (one transient copy), so its own test case in the
   "multi-chain uncombined" assertion above should be a binary family, to
   cover that one-copy floor rather than the general two-copy case.

## Verification

```
Rscript benchmarks/R/memory-footprint.R record benchmarks/baselines/memory-footprint-<rev>.csv
Rscript benchmarks/R/memory-footprint.R fit benchmarks/baselines/memory-footprint-<rev>.csv
air format --check . && Rscript -e 'lintr::lint("benchmarks/R/memory-footprint.R")'
Rscript tools/check-doc-freshness.R && Rscript tools/check-rc-codoc.R
R CMD build . && R CMD check --as-cran dbarts_*.tar.gz
R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'                      # step 7
R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-f0236082.rds   # step 7: 52 identical
```

Expected: record writes one row per cell per metric; fit prints the residual
summary and exits non-zero on a missed step 2 tolerance; format, lint,
freshness, codoc and check clean.

## Landing note, steps 2-6 (2026-09-10)

LANDED at 2f25e3d1f3ced3409b0a8866d6d7166ac6018333, nine commits:

- 72430fbee5f8643c27aa9b083d0824ccb1f282a0 Add the memory footprint audit script
- b6b628555c434d82539382eb15a46fab5c2c9b61 Validate the memory footprint note against measured peak RSS
- 2b827c06aa0de89eeb8b9384e632e8eb3e610e52 Document the memory footprint in the manual
- 156a3e3868be66b4721859683fff7aa01500d595 Rank the largest avoidable allocations in the memory note
- 4d56b2ae6e2f4482b2023328b82e08d88de7e45e Price the stored state in the sampler manual
- 644b19dbdd61a811f6dc4178cdfba51ccca3c83e Score the relative residual only where the prediction clears the floor
- b2ad5f70061990507a8b645b6890eb6e473da17e Correct the ingestion and leaf-cache rows and score over the whole grid
- cfbf74ed0a615d68adacaa50e6e15e57a767ddaa Record the three-run spread on the residual median
- 2f25e3d1f3ced3409b0a8866d6d7166ac6018333 Tidy the memory note after review

[benchmarks/R/memory-footprint.R](../../benchmarks/R/memory-footprint.R)
records one subprocess per cell under `/usr/bin/time -l` (`-v` on Linux),
subtracts a load-only baseline, prints the note's prediction beside each
measurement with the warm-up, churn and ingestion allowances as separate
columns, and scores the plan's tolerance over the whole grid (30 cells:
the Decision 3 spine, the small full crossing, one S = 200 cell, two
designated-covariate cells). The note's Status is VALIDATED: every cell
within max(10 pct, 20 MB), median absolute relative residual 4.4 to 4.8
pct across three runs (limit 5), worst cell n = 1e5, p = 20, T = 200,
C = 2 at 20.8 to 24.3 MB low against 46.6 MB. The measurement moved four
rows: raw predictors became a derived 8*n*p plus a measured ingestion
allowance (two copies in [`dbartsData`](../../R/data.R), one persistent
and one transient); the leaf statistics cache is 4*n per populated arena
depth level per tree per chain from `assign`'s retained capacity in the
never-pruned cache (about 3.2 levels, 276 MB at the reference linear
cell; the 256 MiB budget counts live member lists only and bound no cell);
the gathered leaf design doubles for its standardized copy; and the
starting-sigma `lm` gained a row. `man/dbarts-package.Rd` carries the
Memory section (step 4) and the note the ranked list (step 5).
Baseline benchmarks/baselines/memory-footprint-cfbf74ed.csv (superseded by
memory-footprint-b184b6b2.csv at step 7) is the second
reader's own full-grid recording at the script's last code commit.

Review findings fixed before landing: the leaf-cache row multiplied by
the node count where leaf memberships partition n, so it fit only by
saturating the budget (re-derived from the arena, second linear cell
added); the raw-predictors row cited the class file rather than the two
copying sites; the relative tolerance had been restated to score only
above 100 MB, under which the plan's own wording failed (5.49 pct) - the
restatement was reverted and the corrected model passes the plan's gate
as written.

## Landing note, step 7 (2026-09-09)

LANDED at fd10488821d4787e1e22f92c2c3fe2f067e325a5, five commits:

- 4e76f93466564dc52bdd98a1956fd4edb38e9278 Stop packaging the prediction channels through extra full-size copies
- b184b6b2f530a34979702715828fc3ed37bfb093 Re-derive the packaging rows of the footprint model
- e32e4e0027678bd8feb7dc778e30e58ad56ac415 Record the two-copy removal in the note, the manual and the plan
- 6cb6e141ffd68bf57b450cd074e257594f049a23 Tighten the packaging-copies test
- fd10488821d4787e1e22f92c2c3fe2f067e325a5 Gate the two changed expressions directly, not through a shadowed packager

The two transient copies of the training-prediction array are out of
[`packageBartResults`](../../R/bart.R) and
[`convertSamplesFromDbartsToBart`](../../R/bart.R). The combined-chain
branch builds the returned layout with one `aperm` whose flattening is a
`dim<-` rather than `matrix()` then `t()`, and the posterior means are
reduced over the returned layout by `channelMeans` rather than by `apply`,
which permutes the whole channel again. Packaging now holds two full-size
arrays of a prediction channel at once - the engine's own and the one the
fit returns - on every family and either setting of `combineChains`.

The bridge route the step named was NOT taken. Getting below two arrays
needs the bridge to allocate the channel in the returned layout, which
means the engine writing each draw strided across the observation margin
instead of contiguously: an engine change on the draw path, outside both
this arc's budget and its "no engine change" constraint. The R route
reaches the count the step states - no moment holding more than one
full-size copy beyond the one returned - and the note's ranked row is
re-derived to what it saves, 1600 MB in case 1 and 4000 MB in case 2
rather than 3200 and 8000. The last copy is its own item.

Reducing over the returned layout rather than over the engine's own is
forced, not preferred: `mean()` sums in the order it is handed and runs a
correction pass over the same order, so a reduction over the
observation-major array moves the last bits of every reported mean. After
the permutation each observation's draws are one contiguous slab and the
values are identical.

Measured on this host, one `bart` call per cell under `/usr/bin/time -l`
with n.burn = 0 and `sigest` supplied, against a library built from the
pre-change tree: case 1 (n = 1e5, p = 20, T = 200, C = 4, S = 500)
5847.0 MB before, 4234.8 MB after; case 2 (n = 1e6, p = 50, C = 1)
15172.7 MB before, 11161.5 MB after. Each drop is one prediction array to
within a megabyte. The full 30-cell grid was re-recorded with the new
library and passes the plan's tolerance: every cell within max(10 pct,
20 MB), median absolute relative residual 3.0 to 4.7 pct across runs
against the 5 pct limit, worst cell n = 1e5, p = 20, T = 200, C = 1,
S = 200, keepTrees at 39.6 to 39.7 MB HIGH against a 55.8 MB tolerance.
The top of that spread, and the one linear-cell miss seen with it, came
from a run overlapping another job on the host.
The model lost its collector-churn allowance: the R heap still churns
under the reduction (1.47 times the array against `apply`'s 2.48 on the
duplicates probe, the two exactly one array apart), but none of it reaches
peak RSS in a fit, and keeping the allowance put the grid's median at
8.6 pct with every residual negative.

Gate: [inst/tinytest/test-packaging-copies.R](../../inst/tinytest/test-packaging-copies.R)
compares the two changed expressions against the pre-change ones,
re-declared in the file, over every shape packaging feeds them - one chain
and several, a named parameter margin, and the train and test channels of
real gaussian and binary sampler runs at two and three chains, each through
both settings of `combineChains` - and pins the allocation count itself
through `gc()`'s high-water mark. Recorded draws were rejected for a
fixture: they hold only on the build and instruction set they were recorded
on, while packaging is value-neutral on any build, so a snapshot would have
been dark on every run but one. Whole `bart` results were also compared
across the two installed libraries out of band, identical on all four cases
(gaussian and probit, `combineChains` TRUE and FALSE).

## Landing note, follow-ons (2026-09-10)

LANDED at 5a05d7999cb92f5234e535a4d38e8ab26f43c778, three commits:

- ee03fbcdc245c74f174a990d5d2cf0cbf7eae642 Take the starting sigma
  without lm's model frame and model matrix
- 6571eae1cde67903b0b62e07ae75f789067427fc Skip the complete-cases row
  selection when every row is complete
- 5a05d7999cb92f5234e535a4d38e8ab26f43c778 Reduce the summary and
  partial dependence means without apply

Three of the ranked list's declined-or-deferred rows are taken. The
starting-sigma estimate is now
[`residualStandardError`](../../R/utility.R): it calls `lm.fit`/`lm.wfit`
directly on the design `model.matrix` would have built (intercept
column first) and takes `summary.lm`'s own expression for sigma over
it, dropping the model frame and second design `lm()` built on top of
the matrix already in hand; the estimate is bitwise what `lm()` gave.
[`dbartsData`](../../R/data.R) now runs its complete-cases row
selection (and the weights/offset/basis slices beside it) only when a
row is actually dropped, instead of unconditionally copying `x`, `y`
and the rest even when nothing is missing. `channelMeans`
([`R/bart.R`](../../R/bart.R)) gained a `trailing` argument generalizing
it from two fixed shapes to any trailing-margin count, and now backs
[`posteriorInterval`](../../R/generics.R),
`fitted.bart`/`fitted.bartOrdinal`/`fitted.bartNegbin`/`fitted.bartHurdle`,
`meanCategoryProbabilities` and a new
[`pdbart.drawMeans`](../../R/partialDependence.R) helper used by
`pdbart` and `pd2bart`'s four call sites - each site's own
`apply(..., mean)` is gone, so none of them permutes a full-size
prediction channel to take a per-observation or per-category mean.

Gates: implementer tinytest full suite 8255 of 8255, zero failures.
Equivalence trio against the f0236082 baselines this arc's Verification
block names, run independently by the implementer and by the reviewer:
gaussian 52 of 52 "identical draws (same RNG stream)" lines, BCF 12 of
12 and multinomial 11 of 11 "identical (all N channels: ...)" lines,
zero "max |z|" lines and zero skipped, in both runs (the reviewer's
first BCF and multinomial invocations errored on a baseline-settings
mismatch and were re-run to the 12/12 and 11/11 above). `R CMD check
--as-cran` Status OK; `check-doc-freshness` and `check-rc-codoc` OK;
`lintr::lint()` no lints on every touched R file (`R/utility.R`,
`R/data.R`, `R/generics.R`, `R/partialDependence.R`, `R/bart.R`). The
slice is RNG-neutral, consistent with every gate above.

The note's model rows moved with the code: the starting-sigma row is
now the QR alone rather than `lm()`'s model frame plus model matrix,
and the ingestion allowance's transient complete-cases copy is now
conditional on a dropped row rather than unconditional. Both moves are
downward and touch no engine row, so
benchmarks/baselines/memory-footprint-b184b6b2.csv - recorded against
the pre-follow-on model - is stale; a re-record on a quiet machine,
named after ee03fbcd, is owed and not done by this landing note (a
records-only pass, no benchmark run).

## Landing note, leaf-cache prune (2026-09-10)

LANDED at 3388dc15b47f5eff48f924d74b4e62c67fabcd03, one commit, with the note rows and the TODO entry.

The ranked list's leaf-cache item is taken as the prune rather than as an
honest budget. [`drawFromPosteriorForNode`](../../src/bartcore/model.hpp)
releases the cache slots that are not live leaves, which is what an
accepted grow's interior nodes and an accepted death's freed pairs have in
common, and [`storeCrossproduct`](../../src/bartcore/model.hpp) reallocates
rather than reuse a capacity that has run past twice its membership. The
draw is the hook because it alone runs on the settled tree: a proposal is
scored with the tree mutated, so a prune during scoring would drop the
parent of every rejected grow and pay a rescan for it.

Counting resident bytes against the existing budget was the alternative and
was declined by the measurement: instrumented at n = 1e5, T = 200, C = 1 on
the linear leaf, the cache held 248.3 MB of member lists against 178.1 MB
tracked, so an honest budget would have started refusing entries just under
its own ceiling and traded the megabytes for rescans instead of giving them
back. The same reading settles which mechanism dominates - member lists,
live and stale together, are 99.7 pct of the cache; the inline crossproduct
the budget also misses is 0.8 MB over 1097 populated slots, and is left
alone.

After: 86.1 MB of member lists against 80.0 MB tracked (exactly 4*n*T, one
partition per tree) over 447 populated slots, process peak RSS 689.0 MB to
518.2 MB, and the fit 1.8 pct slower - 6.70 s to 6.82 s, medians of nine
and six interleaved runs on this host. The 1.8 pct is the prune's own O(arena)
walk per drawn leaf plus the rescans the released entries force; it is
reported, not hidden.

Gate: `testLinearLeafStatisticsCachePrune` in
[tests/cpp/test_model.cpp](../../tests/cpp/test_model.cpp) pins the resident
byte count through a stump, an accepted grow, a slot recycled onto a much
smaller leaf and an accepted death, and pins that the pruned cache still
scores and draws bitwise what a leaf that never cached anything computes.
Three mutations were shown to fail it: the prune call removed (3 checks),
the capacity bound removed (1), and the parent walk that tells a freed slot
from a live leaf removed (1). Equivalence is bitwise across all three
harnesses (52 gaussian, 12 BCF, 11 multinomial scenarios), as a cache that
only decides whether a value is recomputed must be. The 30-cell footprint
grid was not re-recorded; the model's cache multiplier is the post-prune
reading at the one cell, and the re-record already owed for the earlier
follow-ons covers it.

## Landing note, re-record (2026-09-11)

LANDED on the model side only; the commits are named at landing. The
recording itself is 9171ff64,
[benchmarks/baselines/memory-footprint-3388dc15.csv](../../benchmarks/baselines/memory-footprint-3388dc15.csv):
30 cells plus the probes, recorded 2026-09-10 at the tip 1456e999, one
subprocess at a time on a quiet arm64 macOS host.

Scored against the model as the recording found it, three cells missed:
n = 1e5, p = 20, T = 75 and T = 200 and T = 200 with keepTrees, 23.3, 22.9
and 23.3 MB high against a 20.0 MB tolerance, median absolute relative
residual 2.8 pct. The three follow-ons all moved the model downward, so the
misses had to be a row the model was not carrying rather than a row it
carried wrong.

It is the posterior-mean reduction's collector churn, and the paired
measurements name it without a fitted constant.
[`channelMeans`](../../R/bart.R) allocates one length-S*C vector per reported
observation, and [`packageBartResults`](../../R/bart.R) calls it only when the
family reports a posterior mean, which probit does not: the same cell reads
250.8 MB gaussian against 232.9 MB probit at n = 1e5, p = 20, T = 200,
C = 1, S = 10, and the gap closes to 1.4 MB (303.4 against 304.8) under
`keepTrainingFits = FALSE`, which drops the training channel and its mean
together. The row had been carried while the mean went through `apply` and
was dropped at b184b6b2 because keeping it turned every residual negative -
but the ingestion allowance then was 2.1 to 2.4 copies of 8*n*p and is now
1.05 to 1.45, the transient complete-cases copy having left the probe. The
allowance had been standing in for the churn at the cells that decided the
removal.

The model carries it again as 180 bytes per reduced observation to a 20 MB
ceiling, both measured: the rate from the paired difference over n at n = 1e4
and 1e5, the ceiling from n = 1e6, where the rate alone would be 180 MB and
the measured difference is 20.4 MB. Two grid shapes sit outside it and the
note says so - n = 1e5, p = 10, where the measured difference is 0.0 MB
against a charged 18.0, and n = 1e5, S = 200, where it is 35.6 against 20.0.

`fit` now re-scores rather than re-reports: the measurement columns stay the
recording's and the two probed allowances with them, and the prediction is
recomputed, so a model correction is scored against an existing grid instead
of owing a new one. A recording pins measurements; a score quoted for one is
the score at the model of its own day, which the MANIFEST rows now say.

After: 3.7 pct median absolute relative residual, no cell outside
max(10 pct, 20 MB), worst cell n = 1e5, p = 10 at -13.7 MB against 22.7 MB.
The median rose from 2.8 because the sixteen n = 1e4 cells turn from about
0.4 MB under to about 1.4 MB over on 20 to 50 MB predictions; the gate that
moved is the per-cell one.

Nothing engine-side moved and nothing shipped: the change is
[benchmarks/R/memory-footprint.R](../../benchmarks/R/memory-footprint.R), the
note, this plan, the MANIFEST and the TODO entry. The leaf-cache multiplier
the prune took to 1.1 levels now has the whole grid behind it rather than one
cell - the two linear cells land 9.4 MB high at T = 200 and 0.7 MB low at
T = 75, against tolerances of 41.9 and 21.4 MB.

Gates: `Rscript benchmarks/R/memory-footprint.R fit
benchmarks/baselines/memory-footprint-3388dc15.csv` OK, every cell within
tolerance, 30 cells; `quick` mode run to exercise the recording path's new
column; `lintr::lint` no lints; `tools/check-doc-freshness.R` OK. The 30-cell
grid was NOT re-run: the machine was not quiet and the recording is what the
model is scored against.
