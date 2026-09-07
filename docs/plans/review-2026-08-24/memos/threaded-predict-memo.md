# Threaded predict: design memo

dbarts, branch bartcore, tip 0045507c (read-only). Ruling taken as given: predict's `n.threads` gets
wired, the shipped flat C API changes, consumer lockstep is a cost to enumerate.

## 0. Summary

It existed. 93f354a8 "Add native predict multithreading" (2025-03-03) shipped it on the classic engine
in 0.9-31; it was never ported to bartcore (f3d01938, 2026-07-02, deliberately serial), and the classic
branch carrying it went in edcdf735 (2026-07-03). Dead formal since.

The recommended decomposition is NOT the classic one. Classic split across CHAINS; bartcore's loop's
natural outer axis is the (chain, draw) SLAB, of which there are numChains * numSavedDraws - orders of
magnitude more parallelism, perfectly balanced, no cross-thread reduction, bitwise-identical by
construction. Measured on this build predict is 99.4-99.5% engine time: Amdahl ceiling 3.93x at 4
threads, 7.66x at 8 - the opposite of docs/design/within-chain-threading.md's NO-GO, whose binding
constraint was the SWEEP's 47% parallel fraction. That no-go does not transfer.

## 1. Archaeology

**93f354a8b763000df66d166a7fd50fb74b06949e**, Vincent Dorie, 2025-03-03, "Add native predict
multithreading", 18 files +249/-27. Its NEWS entry is still in the tree at [inst/NEWS.Rd:1914-1919](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L1914-L1919):
"`predict` now accepts a `n.threads` argument and will use native threads to parallelize across chains."

What it added (all reachable on branch `main`):
- `BARTFit::predict(x_test, n, testOffset, numThreads, result)`, main:[src/dbarts/bartFit.cpp:514-604](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L514-L604).
  Decomposition: ONE top-level task per CHAIN via `misc_htm_runTopLevelTasks(threadManager,
  &predictThreadFunction, threadDataPtr, control.numChains)` ([src/dbarts/bartFit.cpp:582-587](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L582-L587)), plus nested
  `misc_htm_reserveThreadsForSubTask` ([src/dbarts/bartFit.cpp:409](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L409), [src/dbarts/bartFit.cpp:415](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L415)) handing leftover threads down inside a chain when
  numChains == 1. Threads were started and stopped around the call ([src/dbarts/bartFit.cpp:515-521](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L515-L521), [src/dbarts/bartFit.cpp:599-603](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L599-L603)).
- `predictThreadFunction` ([src/dbarts/bartFit.cpp:385-508](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/dbarts/bartFit.cpp#L385-L508)): per chain, loop samples, loop trees, `SavedTree::getPredictions`,
  `misc_addVectorsInPlace`, rescale, add offset.
- Flat-ABI twin `dbarts_predictMultiThreaded` on the OLD C++ ABI header
  (inst/include/dbarts/R_C_interface.hpp), defined [src/R_C_interface.cpp:221](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/R_C_interface.cpp#L221). Note the commit's own
  typo: header `...MultiThreaded`, definition `...Multithreaded`. Both gone; zero hits anywhere on
  bartcore.
- Bridge arity 3 -> 4 with rc_getInt validation ([src/R_interface_sampler.cpp:343-352](https://github.com/vdorie/dbarts/blob/93f354a8b763000df66d166a7fd50fb74b06949e/src/R_interface_sampler.cpp#L343-L352) then), plus the two
  R formals still standing today: R5 `predict(x.test, offset.test, n.threads = control@n.threads)` and
  predict.bart's `n.threads = object$fit$control@n.threads`.
- Test: tests/testthat/test-09-generics.R "predict gives same result when single or multi-threaded",
  asserting `predict(fit, x, n.threads = 2L)` equals `fit$yhat.train` for a 1-chain and a 4-chain fit.

**Where it went.** That test became inst/tinytest/test-generics-multithreaded.R verbatim in b9b9b77f
"Replace testthat with tinytest"; c4391bf2/bcf33186 only reformatted it with air. It still passes -
vacuously, since n.threads is discarded before the .Call.

The loss was a fork that never got it, not a removal. f3d01938 "Add tree storage, prediction, and state
serialization to bartcore" introduced `C_dbarts_bartcore_predict` with the comment "the engine runs
prediction serially" and a 3-argument .Call, keeping the R5 `n.threads` formal beside it while classic
callers still took the 4-argument path. edcdf735 "Remove the classic engine from the R surface" deleted
that branch, taking `.Call(C_dbarts_predict, ptr, x.test, offset.test, n.threads)` with it (pre-rebase
twin 7ff59ea4, under tag bartcore-pre-cran-rebase). Neither message mentions the thread count: the loss
is collateral, not a decision.

The disclosure came later, in c5f85804/8744e778 "Fix the interface review's bug list and close the
shipped-docs gaps": [man/dbartsSampler-class.Rd:152-157](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/dbartsSampler-class.Rd#L152-L157) now says n.threads "Currently has no effect:
`run` and `predict` both execute serially regardless of the value passed here ... reserved for a future
per-call override."

**Two doc defects this exposes.** [man/bart.Rd:149-151](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/bart.Rd#L149-L151)'s shared `\item{nthread, n.threads}` still reads
"Integer specifying how many threads to use" with no caveat, and predict.bart's usage block
([man/bart.Rd:49](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/bart.Rd#L49)) lists `n.threads` under it. And [inst/NEWS.Rd:1914-1919](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L1914-L1919) stands unretracted while 1.0-0's
UPGRADING section (:5-20) never mentions the regression.

**Not related:** origin/archive/within-chain-threading (54a60aa8) and
docs/design/within-chain-threading.md are the SWEEP, "CLOSED - NO-GO on every tested hardware,
2026-07-21". `git log -S/-G` over all three refs finds no other threaded-predict work.

## 2. The current path

1. `predict.bart` [R/generics.R:207](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L207); formal [R/generics.R:214](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L214); `n.threads <- as.integer(n.threads)[1L]` [R/generics.R:294](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L294);
   `object$fit$predict(newdata, offset, n.threads)` [R/generics.R:297](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L297).
2. R5 `predict` [R/dbarts.R:1079](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L1079). n.threads is NEVER READ in the body; `.Call(C_dbarts_bartcore_predict,
   ptr, x.test, offset.test)` at [R/dbarts.R:1140](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L1140).
3. `DEF_FUNC("dbarts_bartcore_predict", bartcore_predict, 3)` [src/R_interface.cpp:224](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface.cpp#L224) ->
   `bartcore_predict` [src/R_interface_bartcore.cpp:5760](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L5760) -> shared tail `predictFromSource` [src/R_interface_bartcore.cpp:5641-5745](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L5641-L5745),
   which allocates the result, calls `sampler.predict(source, n, categoryOffset, REAL(result))` at
   [src/R_interface_bartcore.cpp:5718](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L5718), adds the flat offset at [src/R_interface_bartcore.cpp:5723](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L5723), and clones the shape for `predictVariance` at [src/R_interface_bartcore.cpp:5733](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L5733).
4. Facade virtual `SamplerBase::predict` [src/bartcore/facade.hpp:258-260](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L258-L260), forwarded [src/bartcore/facade.hpp:545](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L545).
5. `Sampler::predictColumns` [src/bartcore/sampler.hpp:564-593](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L564-L593) - THE LOOP:

```
for c in chains:                         # numChains
  if savedTreeCapacity > 0:
    for i in 0 .. recordedDraws_:        # numSavedDraws
      dst = out + (c * numDraws + i) * slab
      chains_[c]->predictFromSavedSample{,Multi}(savedSlotForDraw(i), columns, n, dst)
  else:
    chains_[c]->predictFromCurrentTrees{,Multi}(columns, n, out + c * slab)
```

`slab = numTestObservations * numReportedLocations`. Inside `Chain::predictFromSavedSample`
([src/bartcore/chain.hpp:2800-2820](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L2800-L2820)): zero `out`, loop the forest's numTrees saved flattened trees,
`addFlatPredictions` each into `out`, then `out[i] = scale*out[i] + shift`. Full iteration: chains x
draws x trees x observations, observations innermost and index-partitioned inside
`addFlatPredictionsBelow` ([src/bartcore/tree.hpp:1831](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/tree.hpp#L1831)).

**Read-only-ness.** `predictFromSavedSample` (2800), `...Multi` (2865) and
`predictPerForestFromSavedSample` (2938) are `const`: they read
`forest.savedTrees/savedTreeParams/savedTreeMasks`, `data_` and `response_->fitScale()`, writing only
the caller's `out` plus function-local scratch. The Columns adapters are const-read -
`DenseColumns::column` ([src/bartcore/tree.hpp:1716](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/tree.hpp#L1716)) and `PredictorSourceColumns::column` ([src/bartcore/data.hpp:434](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/data.hpp#L434)) return by
value and build everything in the constructor - so ONE `columns` object is safely shared.
`predictVarianceFromSavedSample` ([src/bartcore/chain.hpp:904](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L904)) is read-only in fact but not marked const: a one-word
fix buying compiler-enforced proof. `predictFromCurrentTrees` (2825), `...Multi` (2898) and
`predictPerForestFromCurrentTrees` (2962) call the non-const `flattenTree` (2683), but that is the
capacity == 0 arm - one slab per chain, so no two threads land on one Chain.

**No RNG, no R API.** Grep over [src/bartcore/chain.hpp:2757-3010](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L2757-L3010) finds no `rng`, `unif_rand` or `norm_rand` - and no
`Rf_error`, `ext_printf` or `R_alloc` either. Predict is a pure function of (saved trees, columns):
determinism is a property to preserve, not one to engineer. [inst/include/dbarts/dbarts.h:45-47](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/include/dbarts/dbarts.h#L45-L47) warns
`dbarts_sampler_predict` is main-R-thread-only because it is "R_alloc-backed internally", but that
allocation is `translateSource` ([src/C_interface.cpp:185-241](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L185-L241)), which runs BEFORE `engine.predict(...)`
at [src/C_interface.cpp:792](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L792); the R bridge likewise parses the source before predictFromSource reaches the engine. Threading
strictly inside `predictColumns` touches no R API.

**Every variant funnels through that loop.** Multi-chain is the `c` loop. Multinomial:
`numReportedLocations > 1` selects `predictFromSavedSampleMulti` - K forests summed per slab, then
`softmaxLocationMajor` per row; same slab structure, wider slab. Amplitude/BCF: plain predict is REFUSED
(`refuseUndefinedTestFits` [src/R_interface_bartcore.cpp:2858](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L2858), called [src/R_interface_bartcore.cpp:5768](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L5768) and [src/C_interface.cpp:780](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L780)), so
those take `predictPerForestColumns` [src/bartcore/sampler.hpp:619-637](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L619-L637) - the same slab loop with a forest margin.
Heteroscedastic adds `predictVarianceColumns` [src/bartcore/sampler.hpp:674-680](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L674-L680). Ordinal ([R/generics.R:1221](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1221)) and negbin ([R/generics.R:1354](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1354))
reach the same .Call via `bartcorePredict` ([R/bartcore.R:1457-1473](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bartcore.R#L1457-L1473)); hurdle via two nested `predict()`
at [R/generics.R:1531](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1531); rbart at [R/generics.R:1706](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1706), [R/generics.R:1712](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1712), [R/generics.R:1722](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1722); partial dependence at [R/partialDependence.R:242-365](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/partialDependence.R#L242-L365). And
[R/dbarts.R:849](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L849) `draw$predict(xt, offset.test, n.threads)` is a SECOND site already forwarding a thread
count into the void.

**How fitting already threads, for reuse.** Across chains: `std::thread` workers, `numWorkers =
min(numThreads, numChains)`, spawned per run with SIGINT blocked so only the main thread runs R's
interrupt handler, then joined ([src/bartcore/sampler.hpp:349-420](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L349-L420); grow-from-root twin [src/bartcore/sampler.hpp:1072-1099](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L1072-L1099)). Within a chain
over TEST ROWS: `Chain::routeTestRows` ([src/bartcore/chain.hpp:4033-4066](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4033-L4066)) over a persistent `misc_mt_manager_t
testFitPool_` ([src/bartcore/chain.hpp:5466](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L5466)), budget `numThreads / numChains`, serial below `testFitParallelCutoff = 65536`
([src/bartcore/chain.hpp:5467](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L5467)). Its comment is this memo's argument, already in-tree: "Routing draws no rng and each row writes
its own output slot, so splitting the range across this chain's share of the thread budget yields
byte-identical results at any thread count."

## 3. Decomposition

Reduction order resolves trivially, and the reason should be stated rather than asserted. The only
summation is per (slab, row) over trees: `addFlatPredictions` accumulates into `out[row]` in tree order
t = 0 .. numTrees-1 ([src/bartcore/chain.hpp:2806-2812](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L2806-L2812)), so every (slab, row) pair owns its accumulator. There is NO
cross-slab and NO cross-row reduction anywhere; the final rescale and the multinomial softmax are both
per-row maps. Any partition keeping each (slab, row) pair whole in one thread therefore reproduces the
serial floating-point sequence exactly. Both candidate axes do.

**(a) Over slabs, i.e. (chain, draw) pairs. RECOMMENDED.** Task count numChains * numSavedDraws: 1000 in
the benchmark below, 4000 for a `bart2(n.chains = 4, n.samples = 1000)` fit. Perfectly balanced - every
slab is the same numTrees traversals over the same rows. Zero extra memory (each task writes a disjoint,
already-allocated `out` range) and zero engine surgery below `predictColumns` (per-task scratch is
already function-local). Contiguous slot blocks per worker keep each on a contiguous region of the
slot-major store (`forest.savedTrees[slot * numTrees + t]`). Degenerate case capacity == 0 gives
numChains slabs - the classic decomposition, and cheap (one forest pass, not ndpost).

**(b) Over observation chunks** works when slabs are few, but the row range must thread through
`predictFromSavedSample` into `addFlatPredictions`'s `indices` seeding while `out` stays indexed
absolutely: real surgery, for an advantage that is the capacity == 0 arm, i.e. one forest pass. **(c)
Over chains only**, the classic scheme, caps at numChains, and a `bart()` fit has nchain = 1. Reject
both.

Recommend (a) alone, with a work-based serial cutoff: `routeTestRows` uses a row count (65536), but the
right predicate is total traversals, numChains * numDraws * numTrees * numTest, and at ~2.5 ns per
traversal a 1e7 cutoff is ~25 ms of serial work, comfortably above any spawn cost. State the constant
with its derivation in the comment. Use `std::thread`, spawned once per predict call and joined,
mirroring `Sampler::run` ([src/bartcore/sampler.hpp:385-420](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L385-L420)) including the `pthread_sigmask` SIGINT block under the
same `#ifndef _WIN32` - a worker must never run R's interrupt handler, since a longjmp across threads is
fatal. Predict is one coarse call of milliseconds to seconds, so a one-time spawn is noise, while a
persistent pool adds a lifetime to manage. Do NOT reuse `testFitPool_`: its budget is `numThreads /
numChains` and it is per-Chain, while a slab partition is cross-chain.

## 4. Surface

### 4.1 Header

No per-call override convention exists in inst/include/dbarts/dbarts.h; the only thread token is the
persistent setter `dbarts_sampler_setNumThreads` (X-macro [src/bartcore/sampler.hpp:473-474](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L473-L474), prototype [src/bartcore/sampler.hpp:939](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L939)). Telling C consumers
to bracket predict with it is wrong three ways: it mutates state that also drives `run`'s worker count
and `routeTestRows`'s budget; the restore leg is not exception-safe, since every refusal on the predict
path raises `Rf_error` and longjmps past it; and it is a lie about a read-only replay.

```c
void dbarts_sampler_predict(dbarts_sampler* sampler,
                            const dbarts_predictor_source* xTest,
                            const double* offsetTest, size_t numThreads,
                            double* out);
```

numThreads before `out`, matching the header's rule that outputs come last; the X-macro entry ([src/bartcore/sampler.hpp:449-452](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L449-L452))
and doc block ([src/bartcore/sampler.hpp:841-855](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L841-L855)) change in lockstep. Semantics: `numThreads == 0` means "use the sampler's own
count", `>= 1` is a per-call override that does not persist. Zero-as-defer makes the consumer migration
a literal `0` and reads as deferral; requiring >= 1 forces every consumer to invent a number, or adds a
`dbarts_sampler_numThreads()` query the header lacks. `dbarts_sampler_forestFits` ([src/bartcore/sampler.hpp:993](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L993)) is IN-sample
and unaffected, and there is no flat entry for predictPerForest or predictVariance: one prototype, one
X-macro entry, one doc block.

`DBARTS_C_API_HASH` ([inst/include/dbarts/dbarts.h:142](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/include/dbarts/dbarts.h#L142), 0x6c9776ae1197e8f5ULL) is an FNV-1a over the entry table,
`static_assert`ed against `dbarts_apiToken()` at [src/C_interface.cpp:465](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L465). It MUST be re-baked; the build
fails loudly until it is. Bump DBARTS_C_API_MINOR alongside.

### 4.2 Consumers

**stan4bart** (branch bartcore, tip 54e157b), `LinkingTo: dbarts (>= 1.0-0)` DESCRIPTION:41, built
`-DDBARTS_USE_STUBS`. ONE C call site, its `src/init.cpp` line 342: `dbarts_sampler_predict(fit, &x_test,
testOffset, REAL(result));` inside `predictBART` (`src/init.cpp` line 294). Minimal migration is one token, `0`; it builds
its samplers with `n.threads = 1L` (`R/stan4bart_fit.R` line 515, `R/mvbart.R` line 119) anyway. It hard-equality
checks `dbarts_apiHash()` at `src/init.cpp` lines 962-977, so every installed binary errors at load until
rebuilt: intended lockstep pre-1.0-0. A full version - a 4th SEXP on `predictBART` - additionally
touches `src/init.cpp` line 1056 and [R/generics.R:205](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L205), [R/generics.R:275](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L275), [R/generics.R:965](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L965): 1 line minimum, ~8 at most.

**bartCause** (branch dbarts-1.0, tip 7ae6e83), `Imports: dbarts (>= 1.0-0)`, `NeedsCompilation: no`, no
LinkingTo, no src/. ZERO lines: it reaches predict only by S3 dispatch, forwarding `...` at
[R/generics.R:143-176](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L143-L176) and [R/generics.R:193-230](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L193-L230), so a user-supplied `n.threads` already reaches predict.bart through
it and will simply start working. tests/testthat/test-08-predict.R would catch a numeric change.

### 4.3 Bridge and R

`DEF_FUNC(..., bartcore_predict, 3)` -> 4 ([src/R_interface.cpp:224](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface.cpp#L224)); `bartcore_predict` ([src/R_interface_bartcore.cpp:5760](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5760)) and
`predictFromSource` ([src/R_interface_bartcore.cpp:5641](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5641)) gain the argument, validated with the rc idiom of section 1;
`bartcore_predictPerForest` ([src/R_interface_bartcore.cpp:5851](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5851), DEF_FUNC [src/R_interface.cpp:225](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface.cpp#L225)) and `predictPerForestFromSource` ([src/R_interface_bartcore.cpp:5813](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5813)) should take
it too, being the amplitude family's only replay path.

R5 `predict` ([R/dbarts.R:1079](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L1079)): the formal exists; pass it at [R/dbarts.R:1140](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L1140). Default `control@n.threads`: KEEP,
being what the formal and [man/dbartsSampler-class.Rd:63](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/dbartsSampler-class.Rd#L63) already promise, and a default that probes cores
is exactly what `_R_CHECK_LIMIT_CORES_` exists to catch. Consequence to accept: `bart()` defaults
`nthread = 1L` ([R/bart.R:2634](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L2634)) so a BayesTree-style fit predicts serially unless asked, while `bart2`
defaults `min(guessNumCores(), n.chains)` ([R/bart.R:657](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/bart.R#L657)) so a 4-chain fit gets up to 4. Note the
mismatch for section 7: that budget was sized for CHAIN parallelism, while predict's is over draws.

`predict.bart` ([R/generics.R:294](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/generics.R#L294)) currently does `as.integer(n.threads)[1L]`, silently accepting NA, 0 and
negatives. Replace with the house pattern ([R/rbart.R:82-85](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/rbart.R#L82-L85), [R/xbart.R:398-401](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/xbart.R#L398-L401)): `n.threads <-
coerceOrError(n.threads, "integer")[1L]; if (is.na(n.threads) || n.threads < 1L) stop("'n.threads' must
be a positive integer")`.

Family generics `predict.bartMultinomial` ([R/generics.R:1013](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1013)), `predict.bartOrdinal` ([R/generics.R:1197](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1197)), `predict.bartNegbin`
([R/generics.R:1330](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1330)), `predict.bartHurdle` ([R/generics.R:1614](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1614)) and `predict.rbart` ([R/generics.R:1647](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L1647)) all take `...` and none names
`n.threads`. Add the formal to each and forward it, at the sites listed in section 2, plus
`predictForest` ([R/xbart.R:603](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/xbart.R#L603)), `samplePriorPredictive` ([R/dbarts.R:849](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L849), already carrying it) and the five
partialDependence sites - partial dependence is predict-in-a-loop and the natural beneficiary.
`bartcorePredict` needs the argument too.

**`run`: NOT NOW.** Its R5 formal ([R/dbarts.R:958](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L958)) is equally discarded - `bartcoreSamplerRun(.self,
numBurnIn, numSamples)` takes no thread count - and the Rd reserves the override for both. But run's
budget is not a pure per-call knob: it also sizes `routeTestRows` ([src/bartcore/chain.hpp:4040](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4040)), and a run mutates
state, so an override moves progress reporting and interrupt behaviour mid-object even though the draws
are per-chain-generator and thread-count-independent ([inst/NEWS.Rd:11-12](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/NEWS.Rd#L11-L12)). Separate, smaller slice.

Rd: [man/dbartsSampler-class.Rd:152-157](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/dbartsSampler-class.Rd#L152-L157) loses the no-effect paragraph for the per-call semantics and default;
[man/bart.Rd:149-151](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/bart.Rd#L149-L151) gains a sentence that predict's `n.threads` parallelizes the replay across posterior
draws; NEWS gets a 1.0-0 entry.

## 5. Tests and gates

**test-generics-multithreaded.R, rewritten.** Today it asserts `predict(fit, x, n.threads = 2L)` equals
`fit$yhat.train` - true whatever the argument does. Replace with `identical()` (bitwise, not
`expect_equal`'s tolerance) between `n.threads = 1L` and 2L/3L/4L on a fixed fit, across: gaussian at 1
and 4 chains (draw-axis vs chain-axis partition); binary probit with `type = "bart"`; multinomial
(K-forest softmax replay); heteroscedastic (both legs of `list(mean, variance)`); BCF glue via
`predictForests` / `type = "forest"`; per-forest replay; a dgCMatrix and a dbartsMixedMatrix test set
(the `PredictorSourceColumns` adapter rather than `DenseColumns`); and one test set under the serial
cutoff and one over it. Use thread counts that do not divide the slab count evenly (3L against 4 chains
x 20 draws) - where an off-by-one partition shows.

**A probe that threads actually run,** since the identity test passes if the wiring is a no-op - the
failure mode being fixed. Do both, they are cheap. (a) tests/cpp: call the engine predict at partition
counts 1, 2, 3, 7 on a fixture of known slab count and `memcmp` the buffers, plus assert the slab ->
worker map is a partition (covers, disjoint); -MMD header tracking makes this a normal add. (b) R:
`system.time(...)[["user.self"]] / [["elapsed"]] > 1.2` on a large predict at n.threads = 4 (user CPU
above wall time requires real concurrency), gated on `NOT_CRAN` and `guessNumCores() >= 4`.

**Equivalence trio stays bitwise.** `benchmarks/R/equivalence.R compare
benchmarks/baselines/equivalence-<hash>.rds` must be UNCHANGED - predict consumes no RNG and this slice
touches no sampling code - as must bcf-equivalence and multinomial-equivalence. If any moves, the slice
touched something it should not have. exact-gates is untouched too, but several of its scripts call
predict: free coverage.

**CRAN.** `R CMD check --as-cran` with `_R_CHECK_LIMIT_CORES_=TRUE`. Every existing inst/tinytest fit
passes 1L or 2L, so nothing exceeds 2 today; the new file must respect it too - cap thread values at
`min(4L, dbarts::guessNumCores())` and skip the >2 arms when it is set. The repo has no
`_R_CHECK_LIMIT_CORES_` anywhere (.github/ and tools/ grep empty): new discipline, not a regression.

**Sanitizers.** sanitizers.yaml is ASAN+UBSAN only; no TSAN leg exists. A slab partition writes disjoint
output and reads immutable state, so a race would be a bug in the partition arithmetic, which the C++
partition test catches deterministically and more cheaply than TSAN. Recommend no TSAN job, and that the
ASAN leg run the new multithreaded file instead.

**Speed, same machine.** This build, Apple M1 Max 8P+2E; fit n = 2000, p = 10, ntree = 200, ndpost =
250, nchain = 4, keepTrees, so 1000 slabs; min of 5 reps, gc between:

| nTest | engine .Call | R reshape | predict.bart | engine share |
|-------|--------------|-----------|--------------|--------------|
|  2000 | 0.962 s      | 0.003 s   | 0.967 s      | 99.5%        |
| 10000 | 5.476 s      | 0.017 s   | 5.511 s      | 99.4%        |

Throughput 2.40-2.74 ns per (row x tree x slab) traversal, single-threaded. Amdahl ceiling over that
parallel fraction: 1.99x @2, 3.93x @4, 7.66x @8.

THAT IS A CEILING, NOT A SPEEDUP - no threaded predict exists to measure. The honest statement:
predict's serial fraction does not cap it, unlike the sweep's; whether memory bandwidth does is unknown,
and the 10000-row case writes 76 MB of output while streaming the whole saved-tree store, which is the
shape of a bandwidth-bound problem. The gate before landing is a real scaling curve at n.threads in {1,
2, 4, 8} on a quiet box, as a benchmarks/R/bench-predict.R counterpart to bench-sampler.R. Enable the R
default only if it beats ~2x at 4 threads; below that the honest outcome is the within-chain note's, and
the surface change is still worth landing for the C consumer.

## 6. Cost

Dense lines, excluding comments (house style roughly doubles the engine):

| Layer | Lines | What |
|---|---|---|
| engine | ~60 | slab partition + std::thread fan-out + sigmask in `predictColumns` ([src/bartcore/sampler.hpp:565](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L565)), shared with `predictPerForestColumns` ([src/bartcore/sampler.hpp:621](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L621)) and `predictVarianceColumns` ([src/bartcore/sampler.hpp:674](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L674)) via one helper; `numThreads` on 3 facade virtuals ([src/bartcore/facade.hpp:258](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L258), [src/bartcore/facade.hpp:267](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L267), [src/bartcore/facade.hpp:272](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L272)) and 3 forwarders ([src/bartcore/facade.hpp:545](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L545), [src/bartcore/facade.hpp:549](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L549), [src/bartcore/facade.hpp:554](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L554)); `const` on [src/bartcore/chain.hpp:904](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L904) |
| flat API | ~10 | [src/C_interface.cpp:773](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/C_interface.cpp#L773) signature + pass-through; [inst/include/dbarts/dbarts.h:449](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/include/dbarts/dbarts.h#L449), [inst/include/dbarts/dbarts.h:853](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/include/dbarts/dbarts.h#L853) + doc; re-bake DBARTS_C_API_HASH ([inst/include/dbarts/dbarts.h:142](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/include/dbarts/dbarts.h#L142)), bump minor |
| bridge | ~20 | [src/R_interface.cpp:224-225](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface.cpp#L224-L225) arity; [src/R_interface_bartcore.cpp:5641](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5641), [src/R_interface_bartcore.cpp:5760](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5760), [src/R_interface_bartcore.cpp:5813](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5813), [src/R_interface_bartcore.cpp:5851](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5851) + one rc_getInt block |
| R | ~45 | [R/dbarts.R:1140](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L1140) (+ [R/dbarts.R:1153](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/R/dbarts.R#L1153)); 6 generics gain a formal and forward it; bartcorePredict ([R/bartcore.R:1457](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/bartcore.R#L1457)); predictForest ([R/generics.R:603](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/R/generics.R#L603)); 5 partialDependence sites; one validation block |
| Rd | ~20 | [man/dbartsSampler-class.Rd:152-157](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/dbartsSampler-class.Rd#L152-L157), [man/bart.Rd:149-151](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/man/bart.Rd#L149-L151), 6 usage blocks, NEWS |
| tests | ~180 | rewritten test-generics-multithreaded.R (~120) + tests/cpp partition test (~60) |
| bench | ~90 | benchmarks/R/bench-predict.R + baseline |
| consumers | 1-8 | stan4bart `src/init.cpp` line 342 (+7 optional); bartCause 0 |

Risks. (1) The facade vtable change forces a full rebuild of every consumer of the header-only engine;
CLAUDE.local.md's `--preclean` rule applies and benchmarks/kernels binaries must be deleted by hand (no
header dep tracking). (2) The hash re-bake breaks every installed stan4bart binary at load until
rebuilt: intended pre-release, and the error message says so. (3) Bandwidth, not correctness, is the
live risk to the payoff. (4) Six generics gaining a formal is six usage blocks that can drift from the
Rd. (5) `predictFromCurrentTrees`'s non-const `flattenTree` is safe under a slab partition ONLY because
capacity == 0 gives one slab per chain: an invariant, not an accident, wanting a comment and an assert.

**Slice plan: TWO.** S1, surface + engine: header change, hash re-bake, bridge arity, engine partition,
R wiring, Rd, rewritten test-generics-multithreaded.R, tests/cpp partition test, stan4bart one-liner -
this must land before 1.0-0 locks the public face. S2, measurement + tuning: bench-predict.R, its
baseline, the tuned cutoff, the default decision; it can follow, since S1 with a conservative cutoff is
correct at any constant. Do NOT split the header change out of S1 - landing the parameter without the
parallelism ships a lie in the doc block and burns the consumer rebuild for nothing.

## 7. Open decisions

1. **`numThreads == 0` = defer to the sampler, or require >= 1?** RECOMMEND 0-as-defer: stan4bart's
   migration becomes a literal `0`, it reads as deferral, and the header has no thread-count query to
   defer to otherwise.
2. **predict's R default: `control@n.threads`, or something predict-shaped?** RECOMMEND
   `control@n.threads` for 1.0-0 - what the formal and Rd already promise, never a CRAN surprise -
   accepting that `bart()` fits predict serially by default. Revisit in S2; if the win is large the
   honest alternative is a `dbartsControl` slot or an option, not a silent core probe.
3. **Does `run` get the same override now?** RECOMMEND no, separate slice (4.3). Decide whether the Rd
   item splits in two now or keeps a joint wording with a "predict only" clause.
4. **Does `predictPerForest` / R5 `predictForests` take it in S1?** RECOMMEND yes: it is the amplitude
   family's ONLY out-of-sample path (plain predict is refused there), so leaving it serial leaves BCF
   users with no threaded predict at all. ~6 lines.
5. **[inst/NEWS.Rd:1914-1919](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/inst/NEWS.Rd#L1914-L1919), the unretracted 0.9-31 entry.** RECOMMEND leave it and say nothing: once
   this slice lands the statement is true again, and retracting a regression that never shipped in a
   release is noise.
6. **Cutoff predicate.** RECOMMEND total traversals (numChains * numDraws * numTrees * numTest >= 1e7)
   over `routeTestRows`'s row count - a 200-row / 4000-slab predict is worth threading and a 200000-row
   / 1-slab one may not be. Tune in S2.
7. **Interrupt polling during a long predict.** Absent today (a 5-second predict is already
   un-Ctrl-C-able); this slice does not change it. RECOMMEND deferring, but recording it as a decision.
