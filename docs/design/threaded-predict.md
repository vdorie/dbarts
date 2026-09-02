# Threaded predict

Status: design settled 2026-08-25, after an independent blind critique
of 35 findings (33 accepted, 1 settled on different grounds, 1
disputed as stated with its remedy accepted; section 9). LANDED
2026-08-25; section 11 records what shipped and what differs from the
proposal. Anchor: bartcore 70b30c18 - every citation below was
re-verified live against that tree.

## 1. The ruling

`predict`'s `n.threads` formal existed on every generic and on
`dbartsSampler$predict`/`$predictForests` but was inert:
man/dbartsSampler-class.Rd then said plainly that "`run` and `predict`
both execute serially regardless of the value passed here" (that item is
now split - section 6, and Rd:154-170 live).
Vincent ruled that this formal gets wired to real threading. Because
dbarts is pre-1.0-0, the shipped `dbarts.h` C API is free to change to
carry it - no consumer's ABI is frozen yet, and consumer compatibility
is a migration cost to enumerate, not a design input.

## 2. Archaeology

Native predict multithreading shipped once already, in 0.9-31 (commit
`93f354a8`, +249/-27 over 18 files), parallelizing **across chains**;
`inst/NEWS.Rd`'s 0.9-31 entry still says so, and `DESCRIPTION` on
`main` is at Version 0.9-34, so the feature is live across four
released versions, not dead work. The bartcore C++20 rewrite forked
past it during the engine cutover: `main`'s classic engine still
carries the `numThreads` overload of `BARTFit::predict`, but
bartcore's replay entry points (`Sampler::predictColumns` at
sampler.hpp:705, `predictPerForestColumns` at :783,
`predictVarianceColumns` at :850) never gained one, while
`Sampler::run` (:349-420) kept its fan-out.
`docs/plans/archive/interface-review.md`'s F10 item (:199) already records both
formals as "fully inert, not merely serial," and lists "threaded
prediction" on the 2.0-WISHLIST (:543); this design discharges both,
striking the wishlist line. The NEWS entry is corrected in place
rather than left standing (`f04e8686` is precedent for editing NEWS
pre-release): its "across chains" wording would stay false even after
landing, since the axis below is (chain, draw) slabs, not chains.

## 3. The partition axis

`Sampler::predictColumns` (sampler.hpp:705-729) iterates chains x
`recordedDraws_`; each (chain, draw) pair is a SLAB writing a disjoint
`out + (c * numDraws + i) * slab` range. The only accumulation anywhere
in the replay is `fits[indices[k]] += leafValue` inside
`addFlatPredictionsBelow` (tree.hpp:1862), once per row per tree, tree
loop `t = 0..numTrees-1` identical at every entry point: each (slab,
row) pair owns its accumulator and sees the same addend order, so a
partition keeping a (slab, row) pair whole in one thread is bitwise
identical at every thread count. No RNG, no `Rf_error`, no `R_alloc`
runs inside the threaded region - `translateSource`'s `R_alloc`
(C_interface.cpp:187,196,221) runs strictly before, so dbarts.h's
main-R-thread-only contract (:45-47) stays unrelaxed. The same shape
covers `predictPerForestColumns`, `predictVarianceColumns`, and
`predictBlend` (R/generics.R:867) - the route by which a BCF fit's
*plain* `predict` reaches the replay, via `predictForest` (:718).

## 4. Worker bodies and thread-count resolution

`predictFromSavedSample` (chain.hpp:2801) constructs an `indices`
vector plus `blockOffsets` on every slab call - numChains x
numSavedDraws times, contended allocator traffic under a naive
partition - so the design hoists one scratch struct per worker,
allocated on the main thread before the spawn, a serial-path win too.
Each worker body runs inside `try`/`catch`, recording rather than
raising any exception; `Rf_error` fires only
after the join, on the main thread, since an escaping `std::bad_alloc`
inside a bare `std::thread` body is `std::terminate` and would abort
the whole R session. `numWorkers <= 1` runs inline with no spawn,
mirroring `run`'s own `sampler.hpp:377`/`:1266` arms, and is capped at the slab
count (the one-slab-per-chain invariant, asserted where there are no
saved trees, is what makes the non-const `flattenTree` safe without
synchronization).

`numThreads == 0` means "the sampler's own count," stored unvalidated
on the C path - `Sampler::setNumThreads` (sampler.hpp:1197-1200) and
`dbarts_sampler_setNumThreads` (C_interface.cpp:960-962) both store
the value as given, and `R/A_class.R:349-350` guards only the R path -
so the resolved count floors at 1, never 0 workers: without the floor,
`setNumThreads(0)` then `predict(..., 0, out)` would resolve to zero
workers and hand R an uninitialized `Rf_allocVector(REALSXP, ...)`; a
test pins this case. A named `constexpr`, sized from a measured ~2.6
ns per (row x tree x slab) traversal, gates each entry point on its
own traversal count, not forest 0's, so a small predict runs without
spawning threads.

## 5. Header and ABI

The prototype, as shipped (inst/include/dbarts/dbarts.h:899, X-macro in
lockstep at :449-452):

```c
/// numThreads is a PER-CALL override that does not persist: 0 means the
/// sampler's own count; a resolved count below 1 is treated as 1. The
/// replay is bitwise identical at every value.
void dbarts_sampler_predict(dbarts_sampler* sampler,
                            const dbarts_predictor_source* xTest,
                            const double* offsetTest, size_t numThreads,
                            double* out);
```

The change moved two hash literals, both failing loudly and
self-correcting: `dbarts_apiSignatureToken` (C_interface.cpp:525, then
`0x85bd1ef04beb3848ULL`) fires first and prints the new signature token
to paste over itself; a rebuild then fails `dbarts_apiToken() ==
DBARTS_C_API_HASH` (:467) and prints the new layout hash to paste over
`DBARTS_C_API_HASH` (dbarts.h:189). Both were re-signed - live,
C_interface.cpp:525 asserts `0x0b33edcf638a3cd3ULL` and dbarts.h:189
bakes `0xb6c0e97dc0688991ULL`. This section also called for
`DBARTS_C_API_MINOR` (dbarts.h:138) to bump from 0 to 1; it did not, on
the header's own pre-release rule - see section 11.
`facade.hpp` took the parameter on all three predict virtuals -
`predict` (:263), `predictPerForest` (:273), `predictVariance` (:278) -
with no default argument (a default on a virtual binds from the static
type, and every real call goes through `SamplerBase&`); the overrides
(:553, :559, :564) and the dense convenience spellings (:285, :291) take
and forward it too, or those calls would silently have stayed serial.
`tests/cpp/test_facade.cpp`'s virtual enumerator, spy, and conformance
checks (:40-41, :140-146, :771-822) record the `numThreads` a forwarder
was handed and assert it passed through, disproving "the wiring is a
no-op".

## 6. Bridge and R surface

The bridge took the argument on both `.Call` entries -
`dbarts_bartcore_predict`, `dbarts_bartcore_predictPerForest`
(`DEF_FUNC` arity 3 -> 4, live at 4, R_interface.cpp:225-226) - and on
`predictFromSource` (R_interface_bartcore.cpp:5727), `bartcore_predict`
(:5791), `predictPerForestFromSource` (:5850), and
`bartcore_predictPerForest` (:5890), validated with
`rc_getInt(..., RC_VALUE|RC_GEQ 1, ...)`. Four inst/tinytest sites
(test-predict-sparse.R:164, test-multinomial-test-offset.R:411 and
:421, test-multinomial-category-offset.R:436) call these `.Call`s
directly and needed the fourth argument; all four pass it.
`dbartsSampler$predict`/`$predictForests` (R/dbarts.R:1084-1169) pass
their formal through; `bartcorePredict` (R/bartcore.R:1085) and the
test harness's `bartcorePredictPerForest` in `inst/common/bartcoreHandle.R`
(no per-forest wrapper ships in `R/`; `$predictForests` calls the entry
point directly) gain the argument; `R/partialDependence.R`'s five
`sampler$predict(x.test)` sites (:242, :244, :346, :363, :365) take it
from the caller.

The formal is appended **last**, after every existing positional
formal and before `...`, on all six generics: bartCause's
`predict.bartcFit` calls `predict.rbart` positionally, ending in
`group.by, combineChains = FALSE, ...`, so inserting `n.threads`
before `group.by` would pass a factor as a thread count.

| generic | R/generics.R | default expression |
|---|---|---|
| predict.bart | :289 | `object$fit$control@n.threads` (already present at :331, moved last) |
| predict.bartMultinomial | :1283 | `object$fit$control@n.threads` |
| predict.bartOrdinal | :1593 | `object$fit$control@n.threads` |
| predict.bartNegbin | :1870 | `object$fit$control@n.threads` |
| predict.bartHurdle | :2412 | `object$occupancy$fit$control@n.threads` (no `object$fit`) |
| predict.rbart | :2463 | `object$fit[[1L]]$control@n.threads` (a list of samplers) |

`predict.bart`'s validation moved above the two early returns that
preceded it - `return(predictForest(...))` and `return(predictBlend(...))`
- so the amplitude family's value is validated too; live,
`validatePredictThreads(n.threads)` sits at R/generics.R:350, above both,
and `predictBlend` (:867-879) is a forwarding site.
man/dbartsSampler-class.Rd:154-170's joint "run and predict both execute
serially" item is split: predict's half states the per-call override and
the default, and affirms that `_R_CHECK_LIMIT_CORES_` (read only by
`parallel:::.check_ncores`) does not govern native threads; run's half
keeps the inert wording (door D1).

The chosen R default is the fit's own thread count, and that exposure
is stated plainly: 117 of 307 `bart2`/`rbart_vi` calls across
inst/tinytest declare no thread argument, 50 with `n.chains > 1` or
defaulted, so their `control@n.threads` is `min(guessNumCores(),
n.chains)`, up to 4 on a CI box (R/bart.R:670), and under this default
every `predict()` becomes multi-worker on those fits, a posture
already true of `run`. The new tests pin explicit counts of 2 or
fewer.

## 7. Verification and consumer migration

No test asserts a timing ratio. `predictColumns` and its siblings
surface the resolved worker count and slab-to-worker map through a
test-only channel; `tests/cpp` asserts the count is `min(requested,
slabs)` and the map covers every slab exactly once, at counts 1, 2, 3,
and 7. Combined with the facade spy (section 5), both "the argument
reached the engine" and "the engine used it" are proved without a
clock. `inst/tinytest/test-generics-multithreaded.R` compares
`identical()` (not a tolerance) across thread counts 1 through 4,
every response family and predictor source, BCF through both
`predictForests` and plain `predict` via `predictBlend`, and test sets
on both sides of the cutoff.

stan4bart's single call site is one line - the extra `numThreads`
argument on its `dbarts_sampler_predict` call - plus a rebuild forced
by its hash pin (`dbarts_apiHash() == DBARTS_C_API_HASH`, checked at
load). bartCause needs zero lines, since `predict.bartcFit` forwards
`...` into dbarts' `predict.bart`/`predict.rbart` - but only under the
append-last rule above, since its positional call would otherwise
misread an inserted argument as `group.by`.

## 8. Measurement, honestly stated

Predict changes no draws, so the equivalence trio (43/12/11 scenarios)
stays bitwise unchanged. A naive ratio of `predict.bart` against its
underlying `.Call`, two independently-timed series, is not reportable
as a percentage on a loaded box; the one directly measured, always
non-negative figure, `convertSamplesFromDbartsToBart`, is 0.30% of a
representative `predict.bart` call.

The serial work predict's threaded region does **not** cover - the
flat-offset add over the whole output (R_interface_bartcore.cpp:5804-
5808) and, on the heteroscedastic path, a full `Rf_duplicate` before
the second (variance) fan-out (:5807) - is one pass over the output
per call against `numTrees` passes inside the replay, so the bridge's
serial share is structurally `O(1 / numTrees)`: about 0.5% at
ntree=200, about 1.3% at `bart2`'s default `n.trees = 75L`
(R/bart.R:666). A parallel fraction near 0.98-0.995 gives an Amdahl
ceiling of roughly 3.8-3.9x at 4 threads and 7.0-7.7x at 8 - **a
ceiling from a serial fraction, not a measured speedup**. No scaling
curve has been measured, and the script that would produce one at
`n.threads` in `{1, 2, 4, 8}` on a quiet machine, under the standing
`bench-sampler` rule, has not been written; the R-level default ships
enabled, so that curve gates the cutoff constant, not the default.

## 9. Critique disposition

An independent blind critique, run against a read-only checkout with
the ruling taken as given, produced 35 findings. 33 were accepted
outright and are folded into sections 3-8 above. One finding split:
its two supporting facts for keeping the fit's thread count as
predict's R default were both refuted (`_R_CHECK_LIMIT_CORES_` does
not cover native threads, and "nothing exceeds 2 threads today" is
false), but the decision was reaffirmed on different grounds - the
fit's own count is the constant a user already expects from `run` -
so it stands while the unguarded CRAN exposure (section 6) is stated
plainly rather than argued away. The last finding was disputed as
stated (a citation-drift claim measured against a moved tree), but its
remedy - re-anchoring every citation at one fixed commit - was
accepted regardless and is the anchor used throughout this document.

## 10. Doors

- **D1**: `run`'s own per-call thread override (R/dbarts.R:960-964),
  equally inert today. Deferred: `run`'s count also sizes
  `routeTestRows`, and `run` mutates state, so wiring it moves
  progress reporting and interrupts too - a larger change.
- **D2**: A predict-shaped thread default distinct from the fit's own
  count (a `dbartsControl` slot, or an option), if section 8's
  scaling curve justifies trading away the current CRAN exposure.
- **D3**: Interrupt polling during a long predict. Deferred, and safe
  to defer: a Ctrl-C during `worker.join()` cannot longjmp out of a
  `std::thread` body - R's Unix `SIGINT` handler only sets
  `R_interrupts_pending`, and the jump happens in
  `R_CheckUserInterrupt`, which predict never calls (the bridge's
  poll is `R_ToplevelExec`-wrapped and reached only by `run`) - so
  the mirrored SIGINT mask is defense in depth, not necessity.

## 11. Landing note

Landed on bartcore. What shipped follows sections 3-7, with three
recorded departures.

**Version constants did not move.** Section 5 called for bumping
`DBARTS_C_API_MINOR` from 0 to 1 alongside the hash re-bake. The
header's own rule forbids it: the minor-version block states that the
constants "become a compatibility contract at the first release ...
and they do not move before it", and the struct field boundaries say a
pre-1.0-0 change "moves no version constant". A signature that GAINS a
parameter is not additive anyway - it is source- and ABI-breaking, so
a minor bump would misdescribe it. Both baked literals were re-signed;
they live in `DBARTS_C_API_HASH` (inst/include/dbarts/dbarts.h) and
`src/C_interface.cpp`, re-baked at every header change, and section 5
has them live. inst/tinytest/test-capi.R pins the current hash, keeps
the superseded ones as negative assertions, and asserts `c(1L, 0L)`.

**partialDependence.R was not touched.** Section 6 listed its five
`sampler$predict(x.test)` sites as forwarding sites. They need no
edit: the R5 method's own default IS `control@n.threads`, the
sampler's own count, so an explicit forward would pass the identical
value. They thread at the fit's count today, through the default.

**The observable carries a cutoff override.** Section 7's deterministic
channel reports the resolved count, the worker count and the
slab-to-worker map. It also carries a writable cutoff override, because
without one every fixture small enough for a unit test falls below the
traversal cutoff, collapses to one worker, and makes an
identity-across-thread-counts assertion vacuous. Both halves are
test-only: `bartcore::predictPartition` in src/bartcore/sampler.hpp,
reached from R through `dbarts_bartcore_lastPredictPartition` and
`dbarts_bartcore_setPredictParallelCutoff`, registered under
src/R_interface.cpp's "below: testing" section beside the SIMD
dispatch knob. Nothing in the R surface reads either.

Mutation proofs run against the pins: a worker that skips its last slab
fails 11 tinytest assertions and both engine identity checks; a
partition-dependent tree-sum order fails 6 and the same two; forcing
the resolved count to 1 inside the engine fails 10 tinytest assertions,
4 engine partition checks and all three facade pass-through checks.

Consumer migration: stan4bart src/init.cpp:342 inserts a `0` before
`REAL(result)` and rebuilds (it pins the hash by equality at
src/init.cpp:972). bartCause needs no edit - its
`predict(object$fit.trt, x.new, group.by, combineChains = FALSE, ...)`
at R/generics.R:162 stays correct precisely because the formal is
appended last.

Open, unchanged: section 8's scaling curve, whose benchmark script has not
been written; the tuned cutoff constant (`Sampler::predictParallelCutoff`,
1e7 traversals); and whether the R-level default stays the fit's own count.
The implementation is correct at any cutoff value.
