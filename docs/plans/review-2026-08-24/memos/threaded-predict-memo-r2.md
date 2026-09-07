# Threaded predict: design memo, REVISION 2

Anchor RE-CUT: origin/bartcore **33837351** (was 0045507c; 42b12ac7 and 33837351 landed during
review). Every line reference below is at 33837351 and was re-verified there. Revision 1 is
superseded; the critique's 35 findings are dispositioned in section 1.

## 1. Disposition of the 35 findings

| # | Sev | Disposition | Effect on the design |
|---|-----|-------------|----------------------|
| 1 | ok | ACK | Archaeology stands as written; no change. |
| 2 | MAJ | **ACCEPTED** | Verified: `main:DESCRIPTION` = 0.9-34 (2026-08-20); 93f354a8 bumped to 0.9-31; `main:src/dbarts/bartFit.cpp` still carries the numThreads overload (7 hits). 0.9-31 SHIPPED and 1.0-0 is a live regression against four releases. Precedent f04e8686 (2026-08-20) edits NEWS in place. Decision 5 REVERSED - the entry is corrected in place (4.6). |
| 3 | min | **ACCEPTED** | docs/plans/archive/interface-review.md F10 ([docs/plans/interface-review.md:198-205](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/interface-review.md#L198-L205)) records both formals "fully inert, not merely serial"; [docs/plans/interface-review.md:543](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/plans/interface-review.md#L543) lists "threaded prediction" on the 2.0-WISHLIST. Both cited; the wishlist line is struck by this slice. |
| 4 | ok | ACK | Independent `mutable` / SIMD-dispatch audit strengthens the no-shared-state claim; folded into 2.1. |
| 5 | ok | ACK | Bitwise-identity-by-construction survives an adversarial read. Unchanged. |
| 6 | MAJ | **ACCEPTED** | Scratch is per-SLAB, not per-task ([src/bartcore/chain.hpp:2805-2825](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/chain.hpp#L2805-L2825) constructs `indices`/`blockOffsets`/`raw` on every call). Design now hoists scratch per WORKER before the spawn, and wraps each worker body in try/catch with a rethrow after the join (2.2). |
| 7 | MAJ | **ACCEPTED** | `Sampler::setNumThreads` ([src/bartcore/sampler.hpp:1011-1014](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/sampler.hpp#L1011-L1014)) and `dbarts_sampler_setNumThreads` ([src/C_interface.cpp:874-877](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/C_interface.cpp#L874-L877)) store unvalidated; [R/A_class.R:349-350](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/A_class.R#L349-L350) guards only the R path. Resolved count is now floored at 1 (2.3). |
| 8 | min | **ACCEPTED** | The one-slab-per-chain invariant is asserted in the capacity == 0 arm and stated in the comment (2.2). |
| 9 | min | **ACCEPTED** | Cutoff is now per-entry-point traversal count (2.4). |
| 10 | ok | ACK | No RNG/R API in the threaded region; [inst/include/dbarts/dbarts.h:45-47](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/include/dbarts/dbarts.h#L45-L47)'s main-thread contract is NOT relaxed. |
| 11 | min | **ACCEPTED** | `bartcore_checkInterrupt` ([src/R_interface_bartcore.cpp:4238-4244](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/R_interface_bartcore.cpp#L4238-L4244)) is R_ToplevelExec-wrapped and reached only by run, so a Ctrl-C during join cannot longjmp. Deferring the poll is SAFE, not merely unchanged (door D3). |
| 12 | ok | ACK | Mirror run's fan-out; do not reuse `testFitPool_`. Unchanged. |
| 13 | min | **ACCEPTED** | `numWorkers <= 1` runs inline on the caller's thread, and numWorkers is capped at the slab count (2.2). |
| 14 | MAJ | **ACCEPTED** | Verified in base R: `parallel:::.check_ncores` is the only reader, and it refuses on "%d simultaneous processes spawned" - processes, never `std::thread`. The safety claim is STRUCK; the Rd states the fact affirmatively (4.5). |
| 15 | MAJ | **ACCEPTED** | Re-derived independently with a paren matcher: 307 `bart2(`/`rbart_vi(` calls in inst/tinytest, **117** with no thread argument, **50** of those with n.chains > 1 or defaulted; worst files test-multinomial-surface.R (17), test-hurdle-surface.R (7), test-argument-surface.R (3) - the critique's figures exactly. "Nothing exceeds 2 today" is struck. |
| 16 | min | **ACCEPTED** | TWO literals: [src/C_interface.cpp:461](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/C_interface.cpp#L461) `dbarts_apiSignatureToken == 0x85bd1ef04beb3848ULL` and [inst/include/dbarts/dbarts.h:142](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/include/dbarts/dbarts.h#L142) `DBARTS_C_API_HASH 0x6c9776ae1197e8f5ULL`. Both in the budget (4.1). |
| 17 | min | **ACCEPTED** | No default argument on any facade virtual (static-type binding). The dense convenience spellings ([src/bartcore/facade.hpp:277](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L277), [src/bartcore/facade.hpp:283](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L283)) take the parameter too, or the dbarts.h-shaped calls stay serial. |
| 18 | MAJ | **ACCEPTED** | tests/cpp/test_facade.cpp is an edit site at three places (FacadeVirtual list [src/bartcore/facade.hpp:40-41](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L40-L41), SPY_VOID [src/bartcore/facade.hpp:140-146](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L140-L146), conformance [src/bartcore/facade.hpp:752-785](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L752-L785)). The spy now RECORDS the forwarded numThreads and asserts pass-through. |
| 19 | MAJ | **ACCEPTED** | Four direct 3-argument `.Call` sites verified: [inst/tinytest/test-predict-sparse.R:164](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-predict-sparse.R#L164), [inst/tinytest/test-multinomial-test-offset.R:406](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-multinomial-test-offset.R#L406) and [inst/tinytest/test-multinomial-test-offset.R:415](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-multinomial-test-offset.R#L415), [inst/tinytest/test-multinomial-category-offset.R:431](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-multinomial-category-offset.R#L431). Plus [R/bartcore.R:1472](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/bartcore.R#L1472), [R/bartcore.R:1490](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/bartcore.R#L1490) and [R/dbarts.R:1140](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/dbarts.R#L1140), [R/dbarts.R:1153](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/dbarts.R#L1153). All budgeted. |
| 20 | MAJ | **ACCEPTED** | bartCause's `R/generics.R` line 162 calls `predict(object$fit.trt, x.new, group.by, combineChains = FALSE, ...)` POSITIONALLY. The formal is APPENDED LAST on all six generics, and each gets its own default expression (4.5). |
| 21 | MAJ | **ACCEPTED** | `return(predictForest(...))` at [R/generics.R:265](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/generics.R#L265) and `return(predictBlend(...))` at [R/generics.R:272](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/generics.R#L272) both precede the coercion at [R/generics.R:294](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/generics.R#L294). Validation moves ABOVE both; predictBlend ([R/generics.R:727-738](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/R/generics.R#L727-L738)) joins the site list. |
| 22 | MAJ | **ACCEPTED** | Measurement re-taken paired, alternating, per-round, loadavg stamped (section 3). No share above 100% is reported. |
| 23 | MAJ | **ACCEPTED** | The offset add ([src/R_interface_bartcore.cpp:5736-5739](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/R_interface_bartcore.cpp#L5736-L5739)) and `Rf_duplicate` ([src/R_interface_bartcore.cpp:5745](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/R_interface_bartcore.cpp#L5745)) are serial work inside the `.Call` the partition never touches. Section 3 now measures the offset and heteroscedastic configurations and gives the structural O(1/numTrees) bound. |
| 24 | min | **ACCEPTED** | The bandwidth presumption is dropped: per slab the engine re-touches a 32 KB output numTrees times; the shape is latency/branch bound. The scaling curve remains the gate. |
| 25 | - | AGREE | Decision 1 kept, with finding 7's clamp. |
| 26 | - | **SPLIT** | The two supporting facts ARE refuted (14, 15 accepted, both claims struck). The DECISION is settled the other way on different grounds by the maintainer: the fit's thread count is the constant a user expects. Recorded, not relitigated; the unguarded exposure is stated plainly (4.5) and the new tests pass explicit values <= 2. |
| 27 | - | ACCEPTED | The Rd item cannot keep joint wording; splitting it is forced. |
| 28 | - | ACCEPTED | predictBlend routes a BCF fit's PLAIN predict through predictForest, so predictPerForest in S1 is required, not optional. |
| 29 | - | ACCEPTED | See 2. |
| 30 | - | ACCEPTED | Named constexpr carrying the ns/traversal derivation, per-entry counts. |
| 31 | - | ACCEPTED | Recorded as door D3 with finding 11's stronger reason. |
| 32 | min | **ACCEPTED** | `git show --stat 63df524e`: +651/-18 across 15 files (202 test lines + 112 tests/cpp), with no ABI change. Budget re-cast as diff lines (4.7). |
| 33 | MAJ | **ACCEPTED** | The timing ratio is gone. The observable is deterministic: resolved worker count and slab->worker coverage, surfaced through a test-only channel. |
| 34 | min | **ACCEPTED** | [.github/workflows/sanitizers.yaml:72-82](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/.github/workflows/sanitizers.yaml#L72-L82) runs `tinytest::test_package("dbarts")` - the whole suite - so there is nothing to "run instead", and ASAN/UBSAN detect no races. A `-fsanitize=thread` arm over tests/cpp replaces "no TSAN". |
| 35 | min | **DISPUTED (as stated), remedy ACCEPTED** | At the memo's stated anchor 0045507c the cited numbers are EXACT: `git show 0045507c:src/R_interface_bartcore.cpp` gives predictFromSource [src/R_interface_bartcore.cpp:5641](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5641), refuseUndefinedTestFits [src/R_interface_bartcore.cpp:2858](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L2858), bartcore_predict [src/R_interface_bartcore.cpp:5760](https://github.com/vdorie/dbarts/blob/0045507c17f93c9fcc738645ce6586484874c4e0/src/R_interface_bartcore.cpp#L5760) - the memo's values. The +13 drift is 42b12ac7 and 33837351 landing during review (+35/-10 in that file). The critique measured against a moved tree. Remedy accepted regardless: R2 is anchored at 33837351 ([src/R_interface_bartcore.cpp:5654](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5654), [src/R_interface_bartcore.cpp:2871](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L2871), [src/R_interface_bartcore.cpp:5773](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5773)). |

**Score: 33 accepted, 1 split (26), 1 disputed as stated with the remedy accepted (35).**

## 2. The amended design

### 2.1 What is unchanged

The partition axis. `Sampler::predictColumns` ([src/bartcore/sampler.hpp:565-593](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L565-L593)) iterates chains x
recordedDraws_, each (chain, draw) SLAB writing a disjoint `out + (c * numDraws + i) * slab`
range. The only accumulation in the whole replay is `fits[indices[k]] += leafValue` inside
`addFlatPredictionsBelow` ([src/bartcore/tree.hpp:1831](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/tree.hpp#L1831)), once per row per tree, with the tree loop
t = 0..numTrees-1 in every entry point: each (slab, row) pair owns its accumulator and sees the
same addend order. There is no cross-slab and no cross-row reduction, so any partition keeping a
(slab, row) pair whole in one thread is BITWISE identical at every thread count. No RNG, no
`Rf_error`, no `R_alloc` inside the threaded region; `translateSource`'s R_alloc
([src/C_interface.cpp:185-241](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/C_interface.cpp#L185-L241)) runs strictly before the engine call, so [inst/include/dbarts/dbarts.h:45-47](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/inst/include/dbarts/dbarts.h#L45-L47)'s
main-R-thread-only contract is unchanged and must not be relaxed.

### 2.2 Worker bodies: hoisted scratch, caught exceptions, inline single worker

```
resolve numWorkers                       # 2.3
if numWorkers <= 1: run inline on the caller's thread, return   # finding 13
scratch[w] for w in 0..numWorkers        # allocated HERE, main thread
failed[numWorkers] = {false}; firstError[numWorkers]
block-partition slabs contiguously over workers
spawn with SIGINT blocked (#ifndef _WIN32), restore the mask, join
each worker body: try { for slab in myBlock: replay(slab, scratch[w]) }
                  catch (const std::exception& e) { failed[w]=true; firstError[w]=e.what(); }
                  catch (...)                     { failed[w]=true; firstError[w]="unknown"; }
after the join: if any failed[w] -> Rf_error with the first message
```

Two amendments, both from finding 6. **Scratch is hoisted per worker.** Today
`predictFromSavedSample` ([src/bartcore/chain.hpp:2805](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2805)) constructs `std::vector<size_t> indices(numTest)` plus
`blockOffsets`, and the multinomial twin ([src/bartcore/chain.hpp:2870](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2870)) additionally `raw(numTest * K)`, on EVERY slab -
numChains * numSavedDraws times, which the partition would make concurrently contended. Passing a
per-worker scratch struct in is a serial-path win as well. **Exceptions are caught and rethrown on
the calling thread.** An escaping `std::bad_alloc` in a `std::thread` body is `std::terminate` -
it aborts the R session, and a large predict is exactly where bad_alloc is plausible. Nothing in
the worker may call `Rf_error`; the raise happens after the join, on the main thread.

The capacity == 0 arm gets an assert: one slab per chain means no two workers touch one Chain,
which is what makes the non-const `flattenTree` ([src/bartcore/chain.hpp:2683](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2683), reached from [src/bartcore/chain.hpp:2830](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2830), [src/bartcore/chain.hpp:2903](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2903),
[src/bartcore/chain.hpp:2967](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2967)) safe. The comment states that a row-axis partition would have to privatize the flatten
buffers (finding 8).

`numWorkers = min(resolvedThreads, numSlabs)`, so a 4-thread predict of a 1-chain no-keepTrees fit
does not spawn three idle workers.

### 2.3 Thread-count resolution: floored at 1

`numThreads == 0` means "the sampler's own count". That count is unvalidated on the C path
(`Sampler::setNumThreads` [src/bartcore/sampler.hpp:1011-1014](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L1011-L1014), `dbarts_sampler_setNumThreads`
[src/C_interface.cpp:874-877](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/C_interface.cpp#L874-L877) both store as given; [R/A_class.R:349-350](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/A_class.R#L349-L350) guards only the R path), so the
resolution is:

```
size_t resolved = numThreads != 0 ? numThreads : options_.numThreads;
if (resolved < 1) resolved = 1;                 // never zero workers
```

Without the floor, `setNumThreads(0)` followed by `predict(..., 0, out)` yields
`min(0, numSlabs) == 0` workers, nothing writes `out`, and predict hands R the contents of an
uninitialized `Rf_allocVector(REALSXP, ...)`. The doc block says what 0 means when the sampler's
own count is 0, and a test pins `setNumThreads(0)` + predict.

### 2.4 Cutoff: per entry point

A named `constexpr` in the style of `testFitParallelCutoff` ([src/bartcore/chain.hpp:5481](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L5481)), carrying its
derivation: at the measured ~2.6 ns per (row x tree x slab) traversal, 1e7 traversals is ~26 ms of
serial work, comfortably above a few tens of microseconds of spawn cost. The traversal count is
THIS entry's, not forest 0's (finding 9):

- `predictColumns`, single location: numChains * numDraws * numTrees * numTest
- `predictColumns`, multinomial: numChains * numDraws * (sum_f numTrees_f) * numTest
- `predictPerForestColumns`: numChains * numDraws * (sum_f numTrees_f) * numTest
- `predictVarianceColumns`: numChains * numDraws * numVarianceTrees * numTest

## 3. Measurement, re-taken

Protocol per the ruling: ONE R process; 9 rounds; every round times all six arms; arm ORDER
ALTERNATES between rounds (forward on odd, reverse on even) so slow drift cancels; `gc(FALSE)`
before each timing; ratios computed PER ROUND and then summarized. Host Apple M1 Max (8P+2E),
dbarts 1.0.0 built from 33837351. **loadavg start 8.40 / 13.33 / 30.42, end 6.72 / 12.22 / 29.03 -
this box was NOT quiet.** Fit: n = 1000, p = 10, ntree = 200, ndpost = 200, nchain = 2 (400 slabs),
keepTrees; nTest = 4000.

| config | median T (predict.bart) | median E (.Call) | per-round E/T median [min, max] | additive T - E median [min, max] |
|---|---|---|---|---|
| best (no offset, dense, gaussian) | 0.829 s | 0.839 s | 1.000 [0.958, 1.100] | 0.000 s [-0.089, +0.038] |
| offset supplied | 0.848 s | 0.830 s | 0.995 [0.962, 1.021] | 0.004 s [-0.017, +0.033] |
| heteroscedastic (25 variance trees) | 0.935 s | 0.930 s | 0.973 [0.877, 1.015] | 0.025 s [-0.014, +0.131] |

**The ratio statistic is not reportable and no share above 100% is quoted.** T and E are separate
calls, so their ratio can exceed 1 on noise, and it does (max 1.100 in the best case). What the
paired design DOES establish is the additive complement, whose best-case median is 0.000 s against
a 0.829 s call - i.e. **the R-side serial work is at or below this box's noise floor.**

The one figure that cannot be negative, measured directly:
`convertSamplesFromDbartsToBart(raw, 2L, TRUE)` = **0.0030 s of a 0.829 s predict.bart, 0.36%**.
Throughput 2.62 ns per (row x tree x slab), consistent with revision 1's 2.40-2.74.

Finding 23's correction stands and is answered structurally rather than by a noisy delta. Serial
work INSIDE the `.Call` that the partition never touches is one pass over the whole output per
non-parallel step: the flat-offset add ([src/R_interface_bartcore.cpp:5736-5739](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5736-L5739)) and, on the
heteroscedastic path, `Rf_duplicate(resultExpr)` ([src/R_interface_bartcore.cpp:5745](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5745)). The replay is numTrees passes over that
same output. So the bridge's serial share is **O(1 / numTrees)**: ~0.5% at ntree = 200, ~1.3% at
`bart2`'s default n.trees = 75L ([R/bart.R:653](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/bart.R#L653)). It shrinks as the forest grows, and it never
scales with the draw count.

**Honest bound.** Parallel fraction f ~ 0.98-0.995 depending on configuration and tree count, so
the ceiling is roughly 3.8-3.9x at 4 threads and 7.0-7.7x at 8. **These are ceilings computed from
a serial fraction, on a loaded box; they are not a measured speedup, and no threaded predict exists
to measure.** The gate before the R default is enabled is the real scaling curve from
`benchmarks/R/bench-predict.R` at n.threads in {1, 2, 4, 8} on a quiet machine, as
CLAUDE.local.md already requires for bench-sampler. Per finding 24 the bandwidth presumption is
dropped: each slab re-touches a 32 KB output numTrees times and the streamed working set is the
saved-tree store at single-digit MB/s, so the shape is latency/branch bound - which threads well -
but that is a prediction, not a result.

## 4. Implementation slice spec (S1)

Anchor 33837351. Executable without reading the critique.

### 4.1 Header and ABI

Exact new prototype ([inst/include/dbarts/dbarts.h:853](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/inst/include/dbarts/dbarts.h#L853), and the X-macro at [inst/include/dbarts/dbarts.h:449-452](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/inst/include/dbarts/dbarts.h#L449-L452) in lockstep):

```c
/// ... numThreads is a PER-CALL override that does not persist: 0 means the
/// sampler's own count (dbarts_sampler_setNumThreads), and a resolved count
/// below 1 - including a sampler whose own count was set to 0 - is treated as
/// 1. The replay is bitwise identical at every value.
void dbarts_sampler_predict(dbarts_sampler* sampler,
                            const dbarts_predictor_source* xTest,
                            const double* offsetTest, size_t numThreads,
                            double* out);
```

Hash re-bake, in order: (1) make the signature change; (2) build - `static_assert` at
[src/C_interface.cpp:461](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/C_interface.cpp#L461) fails and PRINTS the new `dbarts_apiSignatureToken`; paste it over
`0x85bd1ef04beb3848ULL`; (3) rebuild - `static_assert` at [src/C_interface.cpp:465](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/C_interface.cpp#L465) fails and prints
the new `dbarts_apiToken()`; paste it over `DBARTS_C_API_HASH 0x6c9776ae1197e8f5ULL` at
[inst/include/dbarts/dbarts.h:142](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/inst/include/dbarts/dbarts.h#L142); (4) bump `DBARTS_C_API_MINOR` ([inst/include/dbarts/dbarts.h:104](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/inst/include/dbarts/dbarts.h#L104)) from 0 to 1. Both asserts fail loudly
with instructions, so the sequence is self-correcting.

### 4.2 Engine

- [src/bartcore/facade.hpp:258](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L258), [src/bartcore/facade.hpp:267](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L267), [src/bartcore/facade.hpp:272](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L272) - add `std::size_t numThreads` to the three virtuals.
  **No default argument** (a default on a virtual binds from the static type; every real call goes
  through `SamplerBase&`). Overrides at [src/bartcore/facade.hpp:545](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L545), [src/bartcore/facade.hpp:549](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L549), [src/bartcore/facade.hpp:554](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L554). Dense convenience spellings at [src/bartcore/facade.hpp:277](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L277) and
  [src/bartcore/facade.hpp:283](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/facade.hpp#L283) take the parameter and forward it, or the dbarts.h-shaped calls silently stay serial.
- [src/bartcore/sampler.hpp:565](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L565) `predictColumns`, [src/bartcore/sampler.hpp:621](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L621) `predictPerForestColumns`, [src/bartcore/sampler.hpp:674](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L674)
  `predictVarianceColumns` - one shared fan-out helper (2.2), plus the entry overloads at [src/bartcore/sampler.hpp:539](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L539),
  [src/bartcore/sampler.hpp:554](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L554), [src/bartcore/sampler.hpp:602](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L602), [src/bartcore/sampler.hpp:614](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L614), [src/bartcore/sampler.hpp:655](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L655), [src/bartcore/sampler.hpp:662](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L662). Mirror `Sampler::run`'s spawn/mask/join ([src/bartcore/sampler.hpp:349-420](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/sampler.hpp#L349-L420)).
- [src/bartcore/chain.hpp:2805](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2805), [src/bartcore/chain.hpp:2830](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2830), [src/bartcore/chain.hpp:2870](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2870), [src/bartcore/chain.hpp:2903](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2903), [src/bartcore/chain.hpp:2943](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2943), [src/bartcore/chain.hpp:2967](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L2967) - take a per-worker scratch
  struct instead of constructing vectors; mark `predictVarianceFromSavedSample` ([src/bartcore/chain.hpp:909](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L909)) `const`.
- Do NOT reuse `testFitPool_` ([src/bartcore/chain.hpp:5481](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/bartcore/chain.hpp#L5481)): per-Chain, budget numThreads/numChains, row cutoff.

### 4.3 Bridge

- [src/R_interface.cpp:224](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface.cpp#L224), [src/R_interface.cpp:225](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface.cpp#L225) - `DEF_FUNC(..., 3)` -> 4 for both entries.
- [src/R_interface_bartcore.cpp:5654](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5654) `predictFromSource`, [src/R_interface_bartcore.cpp:5773](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5773) `bartcore_predict`, [src/R_interface_bartcore.cpp:5826](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5826)
  `predictPerForestFromSource`, [src/R_interface_bartcore.cpp:5864](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/R_interface_bartcore.cpp#L5864) `bartcore_predictPerForest`, and the two decls in
  src/R_interface_bartcore.hpp - add the argument, validated with
  `rc_getInt(expr, "number of threads", RC_LENGTH|RC_EQ 1, RC_VALUE|RC_GEQ 1, RC_NA|RC_NO, RC_END)`.
- [src/C_interface.cpp:773](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/src/C_interface.cpp#L773) - signature and pass-through.

### 4.4 R

- [R/dbarts.R:1140](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/dbarts.R#L1140) (pass the existing formal), [R/dbarts.R:1153](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/dbarts.R#L1153) (`predictForests` gains the formal and passes).
- [R/bartcore.R:1457](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/bartcore.R#L1457) `bartcorePredict` and [R/bartcore.R:1486](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/bartcore.R#L1486) `bartcorePredictPerForest` gain the argument.
- [R/partialDependence.R:242](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/partialDependence.R#L242), [R/partialDependence.R:244](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/partialDependence.R#L244), [R/partialDependence.R:346](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/partialDependence.R#L346), [R/partialDependence.R:363](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/partialDependence.R#L363), [R/partialDependence.R:365](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/partialDependence.R#L365) - forward the caller's value.
- [R/dbarts.R:849](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/dbarts.R#L849) already forwards; no change.

### 4.5 The six generics

The formal is **APPENDED LAST**, after every existing positional formal and before `...`, on all
six. bartCause's `R/generics.R` line 162 calls `predict(object$fit.trt, x.new, group.by, combineChains =
FALSE, ...)` positionally, so inserting `n.threads` before `group.by` in predict.rbart would pass a
factor as a thread count. Each gets its own default expression:

| generic | R/generics.R | default expression |
|---|---|---|
| predict.bart | [R/generics.R:207](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L207) | `object$fit$control@n.threads` (already present at [R/generics.R:214](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L214) - move it to last) |
| predict.bartMultinomial | [R/generics.R:1013](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1013) | `object$fit$control@n.threads` |
| predict.bartOrdinal | [R/generics.R:1197](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1197) | `object$fit$control@n.threads` |
| predict.bartNegbin | [R/generics.R:1330](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1330) | `object$fit$control@n.threads` |
| predict.bartHurdle | [R/generics.R:1614](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1614) | `object$occupancy$fit$control@n.threads` (there is no `object$fit`) |
| predict.rbart | [R/generics.R:1647](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1647) | `object$fit[[1L]]$control@n.threads` (`object$fit` is a LIST of samplers) |

Validation moves **ABOVE** both early returns - `return(predictForest(...))` at [R/generics.R:265](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L265) and
`return(predictBlend(...))` at [R/generics.R:272](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L272) - so the amplitude family's value is validated too:

```r
n.threads <- coerceOrError(n.threads, "integer")[1L]
if (is.na(n.threads) || n.threads < 1L)
  stop("'n.threads' must be a positive integer")
```

Forwarding sites: [R/generics.R:265](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L265) predictForest, [R/generics.R:272](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L272)/[R/generics.R:727-738](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L727-L738) **predictBlend** (a BCF fit's PLAIN predict
reaches the replay here), [R/generics.R:297](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L297), [R/generics.R:603](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L603), [R/generics.R:1056](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1056), [R/generics.R:1221](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1221), [R/generics.R:1354](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1354), [R/generics.R:1531](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1531) (hurdle's two inner `predict()`),
[R/generics.R:1706](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1706), [R/generics.R:1712](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1712), [R/generics.R:1722](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/R/generics.R#L1722).

Rd. [man/dbartsSampler-class.Rd:152-157](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/man/dbartsSampler-class.Rd#L152-L157) SPLITS - the joint "run and predict both execute serially"
becomes half false the moment this lands. predict's item states the per-call override, the default,
and **affirmatively that `_R_CHECK_LIMIT_CORES_` does not govern native threads** - it is read only
by `parallel:::.check_ncores`, whose own message is "%d simultaneous processes spawned", so it sees
`makeCluster`/`mclapply` and never a `std::thread`; run's item
keeps the inert-formal wording. Usage lines [man/dbartsSampler-class.Rd:56](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/man/dbartsSampler-class.Rd#L56), [man/dbartsSampler-class.Rd:63](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/man/dbartsSampler-class.Rd#L63). [man/bart.Rd:49](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/man/bart.Rd#L49) and the shared
`\item{nthread, n.threads}` at [man/bart.Rd:149-151](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/man/bart.Rd#L149-L151) gain a sentence that predict's `n.threads` parallelizes the
replay across posterior draws, plus the six generics' usage blocks.

Exposure stated plainly, not hedged: 117 of 307 `bart2`/`rbart_vi` calls in inst/tinytest declare
no thread argument, and 50 of those carry n.chains > 1, so their `control@n.threads` is
`min(guessNumCores(), n.chains)` - up to 4 on a CI box - and under the chosen default every
`predict()` on them becomes 4-worker. That posture already exists in `run`; no check tooling will
flag it either way. The NEW tests pass explicit values <= 2.

### 4.6 NEWS

- [inst/NEWS.Rd:1914-1919](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/inst/NEWS.Rd#L1914-L1919), the 0.9-31 entry, is corrected IN PLACE (precedent f04e8686). It
  describes a feature that shipped in 0.9-31 and still stands in released 0.9-34.
- The 1.0-0 entry states that predict threading is BACK and that it partitions **per (chain,
  draw)**, not "across chains" - the 0.9-31 axis is not this one, so the old wording would stay
  false after the slice.
- [docs/plans/interface-review.md:543](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/docs/plans/interface-review.md#L543) - strike "threaded prediction" from the 2.0-WISHLIST; F10
  ([docs/plans/interface-review.md:198-205](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/docs/plans/interface-review.md#L198-L205)) gets a landed note.

### 4.7 Tests, with the mutation each proves

| test | file | mutation it fails on |
|---|---|---|
| forwarded-numThreads spy | tests/cpp/test_facade.cpp (FacadeVirtual list [tests/cpp/test_facade.cpp:40-41](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/tests/cpp/test_facade.cpp#L40-L41), `SPY_VOID` [tests/cpp/test_facade.cpp:140-146](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/tests/cpp/test_facade.cpp#L140-L146), conformance [tests/cpp/test_facade.cpp:752-785](https://github.com/vdorie/dbarts/blob/42b12ac7f23255a2575c466e7c5a92b405ea54c0/tests/cpp/test_facade.cpp#L752-L785)) | a forwarder at [src/bartcore/facade.hpp:545](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L545)/[src/bartcore/facade.hpp:549](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L549)/[src/bartcore/facade.hpp:554](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/src/bartcore/facade.hpp#L554) that drops the argument - "the wiring is a no-op" |
| worker-count / coverage | tests/cpp, new case | resolution not `min(requested, slabs)`; a partition that skips or double-covers a slab, at counts 1, 2, 3, 7 |
| `setNumThreads(0)` then predict | tests/cpp | zero workers returning an uninitialized buffer (finding 7) |
| bitwise identity | inst/tinytest/test-generics-multithreaded.R (rewritten) | any partition-dependent result. `identical()`, not `expect_equal`: 1L vs 2L/3L/4L across gaussian at 1 and 4 chains, probit `type = "bart"`, multinomial, heteroscedastic (both list legs), BCF via `predictForests` AND plain predict through predictBlend, per-forest replay, dgCMatrix and dbartsMixedMatrix sources, and one test set each side of the cutoff. Thread counts that do not divide the slab count (3L against 4 x 20). |
| arity repair | [inst/tinytest/test-predict-sparse.R:164](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-predict-sparse.R#L164), [inst/tinytest/test-multinomial-test-offset.R:406](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-multinomial-test-offset.R#L406) and [inst/tinytest/test-multinomial-test-offset.R:415](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-multinomial-test-offset.R#L415), [inst/tinytest/test-multinomial-category-offset.R:431](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/inst/tinytest/test-multinomial-category-offset.R#L431) | the four direct 3-argument `.Call`s - a hard error in the first, a wrong-message failure in the other three |

No timing assertion anywhere. The worker count and the slab->worker map are surfaced through a
test-only channel so both "the argument reached the engine" and "the engine used it" are proved
deterministically (findings 18 and 33).

### 4.8 Gates

- **Equivalence trio bitwise, unchanged**: `equivalence` 43 scenarios, `bcf-equivalence` 12,
  `multinomial-equivalence` 11 (counts read from the current baselines). Predict changes no draws;
  any movement means the slice touched sampling code.
- **tests/cpp under ASAN and TSAN**: it is a standalone binary with its own Makefile and no R, so a
  `-fsanitize=thread` arm over the predict partition costs a make flag, not an r-hub container.
  ASAN/UBSAN see no races, and [.github/workflows/sanitizers.yaml:72-82](https://github.com/vdorie/dbarts/blob/f04e86864eaa3fdd18bf4842e5f1be9da6ec55b5/.github/workflows/sanitizers.yaml#L72-L82) already runs the whole tinytest suite, so
  there is nothing to add there.
- **`R CMD check --as-cran`**, with a run under `_R_CHECK_LIMIT_CORES_=TRUE` for the record. It
  constrains nothing native; the new tests bound themselves at <= 2 by construction.
- exact-gates: untouched, but several of its scripts call predict - free coverage.

### 4.9 Budget

Diff lines (finding 32: the comparable per-forest predict slice 63df524e was +651/-18 across 15
files with NO ABI change; this one adds the ABI, the generic fan-out and a bench, so ~700-900):

| layer | diff lines | notes |
|---|---|---|
| engine (facade, sampler, chain) | ~180 | fan-out helper, scratch struct, three virtuals + overrides + dense spellings, per-entry cutoff, invariant assert |
| flat API | ~25 | signature, pass-through, TWO hash literals, minor bump, doc block |
| bridge | ~45 | two DEF_FUNC arities, four entries, two decls, rc_getInt |
| R | ~90 | R5 x2, six generics + predictForest + predictBlend + bartcorePredict x2, five partialDependence sites, one validation block |
| Rd + NEWS + docket | ~70 | Rd item split, six usage blocks, 0.9-31 correction, 1.0-0 entry, wishlist strike |
| tests | ~280 | rewritten tinytest file ~170, tests/cpp ~90, four arity repairs ~20 |
| bench | ~100 | benchmarks/R/bench-predict.R + baseline |
| stan4bart | 1 | stan4bart's `src/init.cpp` line 342: `dbarts_sampler_predict(fit, &x_test, testOffset, 0, REAL(result));` plus a rebuild (it pins the hash by equality at stan4bart's `src/init.cpp` line 972) |
| bartCause | 0 | S3 `...` forwarding only - but ONLY under the append-last rule (4.5) |

**S2** (after S1): bench-predict.R's scaling curve on a quiet box, the tuned cutoff constant, and
the decision on whether the R default stays as ruled. S1 is correct at any cutoff constant.

## 5. Doors

- **D1** `run`'s per-call override. Its formal ([R/dbarts.R:958](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/R/dbarts.R#L958)) is equally inert. Deferred: run's
  count also sizes `routeTestRows` ([src/bartcore/chain.hpp:4048](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/bartcore/chain.hpp#L4048)) and a run mutates state, so the override moves
  progress reporting and interrupt behaviour, not just worker count.
- **D2** A predict-shaped thread default (a `dbartsControl` slot or an option) if S2's curve
  justifies it.
- **D3** Interrupt polling during a long predict. Deferred, and SAFE to defer: R's unix SIGINT
  handler only sets `R_interrupts_pending`, and the jump happens in `R_CheckUserInterrupt`, which
  predict never calls - the bridge's poll (`bartcore_checkInterrupt`, [src/R_interface_bartcore.cpp:4238](https://github.com/vdorie/dbarts/blob/63df524e8b54948c6bdca1820a949f2c5eb17385/src/R_interface_bartcore.cpp#L4238)
  -4244) is `R_ToplevelExec`-wrapped and reached only by run. A Ctrl-C during `worker.join()`
  cannot longjmp out. The SIGINT mask mirror is defence in depth, not necessity.
