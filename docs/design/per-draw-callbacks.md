# Per-draw callbacks

Status: LANDED 2026-09-11 (S1 711942ad, S2 e4063140, S3 db3bd084, S4
4cf56b69, S5 1e24bd55), docs/plans/per-draw-callbacks.md. Proposed 2026-09-10
and revised the same day after an independent verification against the code
found sixteen defects in the first draft; the core proposal stands and the
specifics below replace it. All nine forks of section 9 are settled - see
Decision below and dec-B114 in [docs/decisions.md](../decisions.md). Section
7's per-draw cost is UNMEASURED: the bench-sampler.R scenario landed with S4
but its numbers await a maintainer run on a quiet machine ("measured at
landing:" placeholder, section 7). Anchor: bartcore 4cf56b69.

## Decision

Settled by VD, 2026-09-10, verbatim; recorded as dec-B114 and carried into
the plan's Decision section.

- The arc's motivation: "That strikes me as an argument for more and better
  callback support." / "we should bump up the priority. I also think it would
  be good to have an Rcpp example showing how to have a C callback write to a
  preallocated array, if that makes sense. That would be the default way to do
  it in R, since that wouldn't require blocking. The callback would of course
  use raw pointers or lists or void\*, just not SEXPs."
- Fork 4 (storage opt-out): a new logical `keepFits` over every
  per-observation channel (training, test, variance, per-forest), default
  TRUE, set FALSE automatically when a callback is supplied unless the user
  overrides; `keepTrainingFits` stays as the narrower existing switch;
  variable counts and scalar channels are always kept; `keepTrees` stays the
  recompute path. This mirrors stan4bart's `keep_fits` ("Intended to be used
  with callback"). VD: "Sounds good." The spelling supersedes section 4's
  `keepFits`.
- Fork 2 (registration): VD: "Ship the header entry and its two types now."
  The flat C header gains the setter entry and the two types, one ABI hash
  re-bake, dec-B86 argued as section 3 argues it.
- Fork 5 (stop): VD: "Sure, option 1." The callback returns `int`, nonzero
  aborts the run: a shared cancel flag any worker can set, the worker loop
  reading it, the bridge distinguishing a callback stop from an interrupt, the
  sampler's inconsistent-after-abort state documented.
- Forks 1 and 3 are settled by those words: the hook runs on worker threads as
  a const-pointer observer forbidden to call R, reversing dec-B62's premise as
  section 2 argues; no blocking R-closure variant in this arc. If one is ever
  revived, VD's shape for it is that the machinery calls the closure once to
  learn its return length, preallocates, and copies each draw's result in
  itself (stan4bart's `callbackResultLength` pattern) - the design for a
  post-release item, not this arc.
- Forks 6 through 9 are decided by the orchestrator on the recommendations
  above, and are agent-made: the example is a recipe in
  [dbarts-as-a-component.Rmd](../../vignettes/dbarts-as-a-component.Rmd) (Rcpp
  form, unevaluated) plus a plain-C copy under `inst/tinytest` the suite
  compiles; `xbart` does not get the argument; the hazard family's
  expanded-row meaning is documented only; and undefined channels are handed
  over as null pointers.

## 1. Why now

The memory audit priced a large fit and found the engine is not the problem.
At its two reference cases ([Reference cases](memory-footprint.md#reference-cases))
`yhat.train` and its packaging copy are 3200 MB of a 3891 MB peak (n = 1e5,
p = 20, C = 4, S = 500) and 8000 MB of a 10560 MB peak (n = 1e6, p = 50,
C = 1) - five sixths of the peak in the first case and three quarters in the
second, against an engine share of 17 pct in both. Three copies were live
before the audit's step 7; it removed ONE, measured at 1612.2 MB and
4011.2 MB ("each drop is one whole prediction array to within a megabyte"),
leaving the two priced above. The last packaging copy is its own TODO item;
the returned array is inherent to returning an n x S x C object at all.

VD, 2026-09-10, reading that result: "That strikes me as an argument for more
and better callback support." And:

> we should bump up the priority. I also think it would be good to have an
> Rcpp example showing how to have a C callback write to a preallocated
> array, if that makes sense. That would be the default way to do it in R,
> since that wouldn't require blocking. The callback would of course use raw
> pointers or lists or void\*, just not SEXPs.

The argument is exact. A user who wants posterior means, quantiles over a row
subset, a running loss, or their own strided layout pays 8*n*L*S*C bytes to
materialize every draw only so R can reduce over it afterwards. A hook that
hands each draw to caller code as it is produced lets that user allocate what
the ANSWER needs instead, and the S axis never exists.

## 2. What exists today

**An internal R-level entry.**
[`bartcore_runWithCallback`](../../src/R_interface_bartcore.cpp) drives a run
with an R closure evaluated once per SWEEP. It refuses more than one chain
outright, so the closure runs inline on the main thread and blocks the run.
Its RNG contract is load-bearing: no `GetRNGstate`/`PutRNGstate` bracket,
because the chain's generator never touches R's stream while the closure may
draw from it, so R owns `.Random.seed` throughout. Its error contract is
`R_tryEval` - an error cannot longjmp across `Chain::run`'s C++ frames, so it
becomes a cooperative stop the bridge re-raises. It is a CONDITIONING hook,
firing at the top of each iteration, unthinned and including burn-in. No
caller in the tree since rbart_vi was retired, but not unwanted: it is the
named "smallest dbarts-side enabler" of
[Priority and recommendation](correlated-outcomes.md#priority-and-recommendation).

**A retired flat-C entry.** `dbarts_sampler_setCallback`
(retired: [`dbarts_sampler_setCallback`](../../inst/include/dbarts/dbarts.h))
took `(userData, sampler, chainIndex, sweepIndex, isBurnIn)`, returned 0 to
stop, and was refused whenever chains would run on worker threads
([5. The sweep-boundary hooks that exist](bart-as-a-component.md#5-the-sweep-boundary-hooks-that-exist),
[6. C API and callbacks](public-surface.md#6-c-api-and-callbacks)). dec-B86
([docs/decisions.md](../decisions.md)) trimmed it out with the other entries
no consumer calls, narrowing the enabling-value rule for the C header alone
to "ships when a consumer or a named host design calls for it".

**dec-B62, which this note argues against.** dec-B62 shipped that callback
"refused where chains would run on worker threads"; this proposal reverses
that premise, the note's central claim and not a detail. dec-B62 governs a
CONDITIONING hook whose purpose is to re-enter the host between sweeps -
mutating sampler state, in the R case evaluating a closure - so it must run
where the host can be called. A per-draw OBSERVER over const pointers,
forbidden to call R at all, carries no such requirement, so the refusal's
reason does not transfer; what does transfer is that a worker-side hook
cannot be an R closure. The two coexist: the conditioning hook keeps the
inline-only refusal, the observer is worker-safe by construction. Fork 1.

**The engine machinery.** [`SweepCallback`](../../src/bartcore/chain.hpp) is a
`std::function` the chain calls before every sweep;
[`Sampler::run`](../../src/bartcore/sampler.hpp) forwards it only inline,
because with more than one worker chains run on spawned `std::thread`s that
must never call into R - progress queues through
[`QueuedProgressSink`](../../src/bartcore/sampler.hpp), cancellation is a
relaxed atomic, and since dec-B88 the caller blocks on a condition variable.
[`storeSample`](../../src/bartcore/chain.hpp) writes each kept draw into the
caller-owned `Results` slabs
[`allocChannel`](../../src/R_interface_bartcore.cpp) sized, which is the
firing point a per-draw hook wants. Sampler internals stay R-agnostic and the
bridge converts (dec-A50, superseded by dec-B85), so the hook is a plain C
function pointer over plain arrays.

## 3. The C contract

A per-draw callback, distinct from the pre-sweep `SweepCallback` and not a
replacement for it.

```c
typedef struct dbarts_draw_t {
  size_t structSize;   /* library sets; read through DBARTS_HAS_FIELD */
  size_t chainIndex, drawIndex;   /* 0-based; drawIndex over THIS run */
  size_t numObservations, numTestObservations, numPredictors;
  size_t numReportedLocations;    /* L: 1, or K for multinomial */
  size_t numVariableCountForests, numForests, numAmplitudes;
  size_t numOrdinalThresholds;
  const double *train, *test;                    /* n x L; nTest x L */
  const double *varianceFits, *varianceTestFits; /* n; nTest */
  const double* forestFits;         /* n x numForests, forest-major */
  const double* glue;               /* numAmplitudes, ragged, forest-major */
  const double *splitProbabilities, *logLikelihood;  /* numPredictors; n */
  const double* ordinalThresholds;  /* numOrdinalThresholds */
  const uint32_t* varcount;  /* numPredictors x numVariableCountForests */
  double sigma, k, dispersion, residualDf;  /* NaN where inapplicable */
} dbarts_draw;

typedef int (*dbarts_draw_callback)(void* context, const dbarts_draw* draw);
```

**Every channel storeSample settles is there** - the whole of `storeSample`,
not a selection. A pointer is null wherever the fit does not carry that
channel (the variance pair off heteroscedastic, the forest pair off a
multi-forest coupling, thresholds off ordinal, split probabilities off DART)
and the scalars are NaN rather than absent, so a callback tests the channel,
never the family. Two layout facts: `train` carries L channels with any
offset folded in at L = 1, and `varcount` is
`numPredictors * numVariableCountForests`, forest-major within a draw - one
slab for a single-forest model, K for multinomial and for a multi-forest
amplitude model, `Sampler::run` clamping the count to what the combiner can
report before striding by it.

**The two NaN channels.** `storeSample` NaN-fills `testFits` where the
combiner's `testFitsAreDefined()` is false (a BCF test blend a*mu + b_z*tau
is ill-defined with no test treatment vector) and `logLikelihood` where
`logLikelihoodIsDefined()` is false. For the callback both are NULLED
instead: it has no other channel in which to tell "absent" from "present but
NaN". The R channels keep the fill, their shape being part of the returned
object. Fork 9.

**When it fires.** Once per SAVED draw per chain, on the thread owning that
chain, immediately after `storeSample` settles that draw; `drawIndex` counts
this call's saved draws from 0. Sweeps discarded as burn-in inside one
`run(numBurnIn, numSamples)` call never reach `storeSample` and never fire it
- but see section 4, because `bart()` does not use that form.

**Pointer validity: the call, and no longer.** A kept channel's pointer is
into the R result array the bridge allocated, offset by `drawIndex`, and
outlives the call - but the callback cannot know that, because an opted-out
channel's is a per-chain scratch buffer that chain's NEXT draw overwrites.
Copy or reduce inside the call. The struct is size-first like
`dbarts_results`, so an older consumer reads only the fields it knows.

**No R API inside the callback, ever.** Three independent reasons. *Worker
threads*: R's evaluator, allocator and protection stack are single-threaded,
so `Rf_allocVector`, `Rf_eval`, `R_alloc`, even `PROTECT`, corrupt state the
main thread owns. *GC*: any R allocation may collect, and a worker allocating
while the main thread sits inside `Sampler::run` can free objects nothing has
protected on its behalf. *longjmp*: `Rf_error` longjmps to a context the MAIN
thread established - from a worker it unwinds the wrong stack, and even on
the main thread it skips every C++ destructor between raise and catch, the
defect ["capi-unwind-protect"](../../TODO) records for the existing entry. On
POSIX, SIGINT is blocked in workers, so R's interrupt handler cannot run
there either. For C++ the rule is sharper than "do not allocate": do not
CONSTRUCT OR DESTROY an Rcpp proxy type inside the callback -
`Rcpp::NumericVector` and its siblings touch `Rcpp_PreserveObject` and the
protection stack on construction and destruction, the prohibited operation
even where no allocation is visible. Take the raw `double*` out beforehand.

**Interrupts are dead for the duration of a call.** The main thread's poll
sets a cancel flag workers read at the TOP of the next sweep, so a hung
callback never reaches that read while the main thread stays in
`chainsDone.wait_for`: the session hangs with no Ctrl-C. Inline, the callback
runs ON the main thread and no poll happens at all.

**Return value, errors, and what a stop costs.** `int`, nonzero to stop -
with the plumbing specified, because today none exists: the worker lambda in
`Sampler::run` DISCARDS `Chain::run`'s return value, `cancelFlag` is written
only by the main thread's `pollInterrupt`, and a worker holding several
chains (`c += numWorkers`) moves on regardless. Making the flag real takes
one cancel flag shared by both arms, the worker lambda storing into it when
`Chain::run` returns true, and the bridge telling a callback stop from an
interrupt as `bartcore_runWithCallback` tells `closureStopped` from
`cancelled` - about 30 lines. Two properties are stated rather than fixed:
another chain sees the flag only at its next sweep boundary, and a stop is an
ABORT, since `Sampler::run` returns before advancing `currentSampleNum_` and
`recordedDraws_` while `storeSavedTreeRecord` has already written saved trees
into the slots those cursors count - so results and saved trees are
discarded, the property the interrupt path already has. A C callback cannot
signal an error at all: it records a status in its own context, read from R
after the run, and returns nonzero only to abort. Fork 5.

**Concurrency.** Calls for different chains may run CONCURRENTLY; calls
within a chain are ordered by `drawIndex`. The engine takes no lock, so a
callback touching shared state owns its synchronization; the discipline that
needs none, and the one the example uses, is a write addressed by
`(chainIndex, drawIndex)`, disjoint by construction. An engine mutex would
make every chain wait on the slowest callback: rejected, said in the header.

**Registration: both sites, together** (fork 2). The R bridge's `.Call` run
entry takes a function pointer and a context pointer read out of external
pointers; `dbarts.h` gains `dbarts_sampler_setDrawCallback(sampler, fn,
context)`, null `fn` clearing. The types can neither be held back nor ship
alone: the ABI token folds the `DBARTS_C_API_LIST` signatures, the three
enums and the layouts of `dbarts_results` and `dbarts_predictor_source`, so a
folded-in `dbarts_draw` re-bakes the hash while one left out is an
unversioned ABI layout - and a type with no entry is exactly the irreversible
inclusion with no consumer dec-B86 refuses. dec-B86 is therefore argued
head-on: the named consumers are the documented example, which cannot compile
against a declaration that does not ship, and stan4bart.

## 4. The R surface

`bart()`, `dbarts()` and `dbartsSampler$run` gain one argument, `callback`,
taking either of two things. **(a) A C callback, non-blocking** - VD's
default: a list of two external pointers, one over the function address and
one over the caller's context. The run never re-enters R, so it composes with
`n.threads > 1` and every chain at once, and the R side only checks that both
are external pointers and passes the addresses down. **(b) An R closure,
blocking**, evaluated on the main thread between sweeps and forcing the run
inline; single-chain is `bartcore_runWithCallback` today, and multi-chain
would need workers to enqueue each finished draw and the main thread to drain
the queue - so either buffering draws (reintroducing the array this note is
about) or blocking each worker until its draw is consumed, serializing the
chains behind R's evaluator. Fork 3: do NOT ship (b) here.

**Burn-in is an R-layer problem, not an engine one.**
[`runWithBurnIn`](../../R/bart.R) passes no burn-in count to the engine: it
calls `sampler$run(0L, control@n.burn)` for the burn phase and
`sampler$run(0L, control@n.samples)` for the kept phase, so EVERY burn-in
sweep is a recorded draw at `sampleNum` 0..n.burn-1 and `drawIndex` restarts
at 0 for the second call. A hook fired after `storeSample` therefore sees
every burn-in draw at the `bart()` defaults, and a naive running mean
averages burn-in in. The fix is one line in the R layer: install the callback
on the kept-sample run only, leaving the engine rule ("fires on saved draws")
as stated. The alternative, stopping `bart()` from splitting the run, gives
up the burn phase's draw-neutral speedup (it forces `keepTrainingFits` FALSE
and drops the test set). Driving
[`dbartsSampler$run`](../../man/dbartsSampler-class.Rd) with a real
`numBurnIn` never sees it; calling `$run` twice restarts `drawIndex` each
time, which the example's context must survive.

**The storage opt-out needs its own control.** `keepTrainingFits` cannot
carry it: the test channel is gated only on `numTestObservations > 0` (see
[`installChannel`](../../src/R_interface_bartcore.cpp)), so there is no
`keepTestFits` to set, and the heteroscedastic variance channels
([`varianceTrainExpr`](../../src/R_interface_bartcore.cpp)) and multi-forest
[`forestFitsExpr`](../../src/R_interface_bartcore.cpp) are gated only on the
model carrying them. Those are n*S*C and n*F*S*C - as large as `yhat.train`
or F times larger - so a train-only opt-out fails where the arrays are
biggest. Proposal: one new logical `keepFits`, FALSE when `callback` is
supplied and TRUE otherwise, over every PER-OBSERVATION channel (train, test,
variance train and test, forest fits); scalars and `varcount` are kilobytes
and stay, and `keepTrainingFits` remains the legacy spelling for the train
half. Fork 4.

**What the opt-out costs the returned object**, stated plainly because the
first draft understated it. With no train channel,
[`plot.bart`](../../R/plot.R) and `extract(sample = "train")`
([`extract`](../../R/generics.R)) both ERROR by name, `yhat.train.mean` is
absent ([`packageBartResults`](../../R/bart.R)), and `fitted`/`residuals`
lose their source; with no test channel `yhat.test` and `yhat.test.mean` go
too, which `keepTrainingFits = FALSE` does NOT do today; and `keepTrees`
defaults FALSE, so `predict` recovers none of them. The default callback fit
therefore returns no fitted value at all - the accumulator IS the answer -
which is why fork 4's simpler alternative is not absurd.

**The mechanism, and where the stride lives.** The first draft had this
wrong. `Sampler::run` computes each chain's base as
`results.trainingFits + c * numSamples * n * numLocations` BEFORE any draw is
stored, and `Results` is one struct the sampler slabs, so the bridge cannot
point a single chain's buffer elsewhere and a one-draw buffer would have
chains 1..C-1 writing past its end. Both halves change: the bridge allocates
C one-draw buffers CONTIGUOUSLY per opted-out channel, and `Sampler::run`
takes a per-draw stride (default `n * numLocations`) that both its per-chain
base arithmetic and `storeSample`'s offset read, so a stride of zero lands
every draw of a chain in that chain's one buffer. `bartcore::Results` is
engine-internal, so the field is additive, not an ABI event.

## 5. The Rcpp example

A running posterior mean over the training rows - the reduction that most
often motivates the full array.

```r
Rcpp::sourceCpp("drawMean.cpp")     # makeMeanContext, meanCallbackPtr
acc <- numeric(nrow(x) * nChains)   # the ONLY per-observation allocation
ctx <- makeMeanContext(acc, nrow(x), nChains)
fit <- bart(x, y, n.chains = nChains,
            callback = list(fn = meanCallbackPtr(), context = ctx))
means <- rowMeans(matrix(acc, nrow(x), nChains))
stopifnot(contextStatus(ctx) == 0L)
```

```cpp
#include <Rcpp.h>
#include <dbarts/dbarts.h>   // dbarts_draw; PKG_CPPFLAGS -I the include dir

struct MeanContext {
  double* out;           // REAL() of the R vector, taken BEFORE the run
  std::size_t n, numChains;
  std::size_t* counts;   // per chain: no chain reads another's counter
  int status;            // set here, read from R after the run
};

extern "C" int meanCallback(void* context, const dbarts_draw* draw) {
  MeanContext* ctx = static_cast<MeanContext*>(context);
  if (draw->train == nullptr || draw->numObservations != ctx->n ||
      draw->chainIndex >= ctx->numChains) {
    ctx->status = 1;     // record and continue; nonzero would ABORT the run
    return 0;
  }
  double* out = ctx->out + draw->chainIndex * ctx->n;    // disjoint slice
  const double m = static_cast<double>(++ctx->counts[draw->chainIndex]);
  for (std::size_t i = 0; i < ctx->n; ++i)
    out[i] += (draw->train[i] - out[i]) / m;             // running mean
  return 0;
}

// freeMeanContext (not shown) deletes counts and the context

// [[Rcpp::export]]
SEXP meanCallbackPtr() {   // ...PtrFn: casting a function pointer to void* is
                           // what -Wpedantic flags under R CMD check
  return R_MakeExternalPtrFn(reinterpret_cast<DL_FUNC>(&meanCallback),
                             R_NilValue, R_NilValue);
}

// [[Rcpp::export]]
SEXP makeMeanContext(Rcpp::NumericVector acc, std::size_t n, std::size_t C) {
  MeanContext* ctx = new MeanContext{REAL(acc), n, C, new std::size_t[C](), 0};
  // acc goes in the PROTECTED slot (third argument): that is what keeps the
  // array alive as long as the context is
  SEXP ptr = PROTECT(R_MakeExternalPtr(ctx, R_NilValue, acc));
  R_RegisterCFinalizerEx(ptr, &freeMeanContext, TRUE);
  UNPROTECT(1);
  return ptr;
}
```

A strided layout instead of a reduction is the same shape with a `numDraws`
field added and the body replaced by a `memcpy` into
`ctx->out + (draw->chainIndex * ctx->numDraws + draw->drawIndex) * ctx->n` -
exactly the layout the memory note's second array existed to produce, built
once, in place. It is where the per-run caveat bites: a context reused across
two `$run` calls sees `drawIndex` restart at 0 and overwrites the first run's
slots, so a context is per run unless the caller carries a draw offset.

**Lifetime rules, all the caller's.** `acc` must stay reachable and unmoved
for the run; the external pointer's protected slot ties it to the context.
`REAL(acc)` is read ONCE, before the run: a reallocation on assignment makes
a DIFFERENT vector and the callback keeps writing the old one. No R entry
point, no constructing or destroying an Rcpp proxy type, no throwing. Chains
need disjoint slices, as here, or atomics - which on doubles give up bitwise
reproducibility. Nonzero ABORTS the run, so a callback that merely disagrees
with a draw records a status and returns 0. Getting any of it wrong CRASHES
the session: no condition to catch, a traceback into the engine.

stan4bart's embedding pattern - per-iteration caller-owned buffers, one draw
at a time, chains in separate single-chain samplers - is this shape reached
from C, and is the second consumer. Rcpp itself is not a dbarts dependency:
the repo's precedent for a compiled example under test is
["consumer source not installed"](../../inst/common/capiConsumer.R), plain C
through `R CMD SHLIB`, self-gating on the toolchain. Fork 6.

## 6. Memory consequence

Case 1 (n = 1e5, p = 20, T = 200, C = 4, S = 500), gaussian single forest.
Every `today` row below is a row of
[Reference cases](memory-footprint.md#reference-cases) as it now stands, after
the ingestion guard that note records took 8*n*p off the R-side predictor
block - so its R row reads 33.0 MB here, not the 49.0 MB the pre-guard
revision carried:

| | today | callback, `keepFits = FALSE` |
| --- | --- | --- |
| engine | 658.3 MB | 658.3 MB |
| R, predictors and scalars | 33.0 MB | 33.0 MB |
| yhat.train, two copies | 3200.0 MB | 0 |
| per-draw scratch, 8*n*L*C | 0 | 3.2 MB |
| caller's accumulator, 8*n*C | 0 | 3.2 MB |
| peak | 3891.3 MB | 697.7 MB |

Case 2 (n = 1e6, p = 50, C = 1) carries a much larger R-side predictor block
(808.1 MB): 8000.0 MB of `yhat.train` becomes 8.0 MB of scratch plus 8.0 MB
of accumulator and the peak falls from 10560.3 MB to 2576.3 MB. The engine's
share rises from 17 to 94 pct in case 1 and 17 to 68 pct in case 2, which is
where the audit's remaining ranked items point.

Two families are worth more than this: a heteroscedastic fit adds an n*S*C
variance channel and a BCF fit n*F*S*C forest fits, neither reachable by
`keepTrainingFits` - the argument for `keepFits` - and both now carry their
formula-derived sizes under [Reference
cases](memory-footprint.md#reference-cases), 1600.0 and 3200.0 MB at case 1's
shape. The honest comparison is the manual's cheaper lever:
`keepTrainingFits = FALSE` alone saves the same 3200 and 8000 MB for free, but
DISCARDS the draws.

## 7. Threading interaction

The per-draw cost is an indirect call plus whatever the callback does, once
per saved draw. A sweep at n = 1000, p = 10, T = 75 measures 0.174 ms
(benchmarks/baselines/bench-sampler-127f04ee.csv) and the example's callback
touches 8 KB, so at that size it is a fraction of a percent. How it SCALES is
not asserted: the sweep is cache-resident at small n while the callback is an
O(n) streaming pass over 8 MB per draw at n = 1e6, so the ratio could move
either way. `benchmarks/R/bench-sampler.R`'s `callback` scenario times a
no-op C callback and the vignette's running-mean recipe against no callback
at all, at both reference shapes (`callback-none-*`, `callback-noop-*`,
`callback-mean-*`); maintainer-run on a quiet machine, per Verification.
measured at landing: [pending - the orchestrator records the compared
numbers here once bench-sampler.R has run on a quiet machine].

The load-bearing point is the one VD named: a C callback is the only kind
that can run under `n.threads > 1` without serializing the chains, an R
closure making every worker wait on the main thread's evaluator.

## 8. Scope and sequencing

Minimum shippable surface, pre-release (VD raised the priority, 2026-09-10):

1. Engine: a `DrawCallback` fired from `Chain::run` after `storeSample` and
   forwarded by `Sampler::run` on BOTH paths; the per-draw stride read by the
   per-chain base arithmetic and by `storeSample`; the shared cancel flag the
   worker lambda writes when `Chain::run` returns true.
2. `dbarts_draw`, `dbarts_draw_callback` and
   `dbarts_sampler_setDrawCallback` in `dbarts.h`, one hash re-bake, and the
   bridge's per-run function-pointer argument.
3. R: the `callback` argument on `bart()`, `dbarts()` and
   [`dbartsSampler$run`](../../man/dbartsSampler-class.Rd); `keepFits` and
   its C-draw buffers; `runWithBurnIn` installing the callback on the kept
   run only; the crash, interrupt and lost-fitted-value warnings in the Rd.
4. The example, a seventh recipe in `vignettes/dbarts-as-a-component.Rmd`.
5. Tests: a plain-C callback compiled with `R CMD SHLIB`, its accumulated
   mean checked against the same fit's `yhat.train` mean and self-gating like
   ["consumer source not installed"](../../inst/common/capiConsumer.R); a
   multi-chain run asserting per-chain call counts and ordering; a run
   asserting exactly `n.samples` calls at the `bart()` defaults (the burn-in
   regression); a stop-flag run; heteroscedastic and BCF fits asserting the
   extra channels arrive; and a tests/cpp test that the hook fires once per
   recorded draw.

Waits: the blocking R closure (fork 3), and any per-draw TREE handoff, which
would hand a callback a variable-length structure the engine does not hold
outside `keepTrees`. Budget: roughly 550 lines over the engine hook and
stride, the cancel plumbing, the header entry, the bridge argument, the R
surface and its Rd, the vignette recipe and the tests.

RNG class NEUTRAL: the hook consumes no generator draw and mutates no sampled
state, so a run with no callback registered stays bitwise identical, the
property `SweepCallback` and `shouldCancel` already carry. Gates owed:
tests/cpp; the full tinytest suite; all three equivalence harnesses IDENTICAL
on every scenario; and `bench-sampler.R compare` on a quiet machine, the hook
sitting on the draw path. The new entry is ADDITIVE and re-bakes
`DBARTS_C_API_HASH`; under dec-B111 the version pair is the guard and CI
asserts a changed hash comes with a minor bump, so pre-release the constants
stay at 1/0 with one re-bake. treatSens calls neither the new type nor entry.

stan4bart is dec-B86's second named consumer, and its source confirms the
claim at this landing. It calls neither `dbarts_draw_callback` nor
`dbarts_sampler_setDrawCallback`, so there is nothing to port. Its load-time
handshake bakes in the major/minor pair its own `dbarts.h` copy defines
(1/0) and refuses to load on a mismatch, so S2's re-baked
`DBARTS_C_API_HASH` with the pair held at 1/0 under dec-B111 admits
stan4bart's existing compiled binary unchanged. Its embedding already has
the shape this entry serves: `stan4bart_fit_worker` builds one
single-chain, single-sample `dbartsSampler` per chain (`n.chains = 1L,
n.samples = 1L`) and drives it a draw at a time from R into its own
caller-owned buffers via its per-sweep R closure - so it may adopt
`dbarts_sampler_setDrawCallback` in that closure's place post-release,
without being obliged to. Its `keep_fits` argument - "Logical that, when
false, prevents the sampler from storing each draw. Intended to be used
with `callback`." - is the precedent `keepFits` mirrors (Decision, above).

## 9. Open forks for VD

All nine are settled; the recommendations below are kept as the argument
behind each ruling, and Decision above records the rulings themselves.

1. **Reverse dec-B62's worker-thread refusal for an observer hook?**
   Recommend yes, scoped in words to the observer: dec-B62 refused a hook
   that must re-enter the host, and a const-pointer observer forbidden to
   call R does not. Alternative, hold the refusal for both - which makes the
   callback unusable at the default `n.chains = 4`, the very fit carrying the
   3200 MB array.
2. **Ship the header entry now, or the R route alone?** Recommend the entry
   and the two types together, arguing dec-B86 head-on with the documented
   example and stan4bart as named consumers; one irreversible entry and one
   hash re-bake, free pre-release. Alternatives: types with no entry, which
   dec-B86 refuses and the hash either folds or silently omits; or no header
   change, leaving every example hand-declaring an unversioned ABI layout.
3. **The blocking R-closure variant: ship or not?** Recommend NOT here: no
   multi-chain without a draw buffer (the array we are removing) or
   serialized workers, plus `R_UnwindProtect` and a second error contract, to
   serve a user `$run(0, 1)` in a loop already serves. Reopen it if the
   crash-on-mistake cost of the C route is what users hit.
4. **`keepFits`, or reuse `keepTrainingFits`?** Recommend the new logical
   over every per-observation channel. Cost: one more control frozen at 1.0,
   and a default that silently drops `yhat.test`, which nothing drops today.
   Alternative, reuse `keepTrainingFits`: simpler, surprising nobody, and
   leaving heteroscedastic (n*S*C) and BCF (n*F*S*C) fits paying the whole
   array - the two cases where the saving is largest.
5. **The stop flag: `int` with new plumbing, or `void`?** Recommend `int`,
   with the shared cancel flag, the worker lambda storing into it and the
   bridge telling a callback stop from an interrupt, about 30 lines. The stop
   is an ABORT and other chains see it only at their next sweep boundary.
   Alternative `void`: no plumbing, and a callback that has detected a
   problem can do nothing - less bad than it sounds, since a hung callback
   kills Ctrl-C either way.
6. **The example's home.** Recommend a seventh recipe in
   `vignettes/dbarts-as-a-component.Rmd` (Rcpp form, unevaluated) plus a
   plain-C copy under `inst/tinytest/` the suite compiles. Alternatives:
   `inst/examples/` alone (nothing tests it), or a companion package.
7. **Does `xbart()` get the argument?** [`xbart`](../../R/xbart.R) builds its
   own control with `keepTrainingFits = FALSE` and never surfaces the fit
   arguments, so a callback would thread through the CV driver and fire once
   per fold per rep. Recommend NOT here: xbart holds loss values, not fits.
   Alternative, thread it through for a per-fold custom loss, at the price of
   fold and rep indices in the struct.
8. **What does `numObservations` mean under `family = "hazard"`?** The
   person-period expander turns N subjects into N' rows before the sampler
   sees anything, so a callback's `train` is over expanded rows and its
   indices are not subjects. Recommend documenting it and nothing more - the
   mapping is R-side. Alternative: pass the expansion index, a
   family-specific field on a family-neutral struct.
9. **NaN channels: null the pointer, or hand over the NaN buffer?** Recommend
   nulling `test` and `logLikelihood` where the combiner declares them
   undefined, keeping the R channels' NaN fill. Alternative: expose the NaN
   buffer, a line cheaper and a trap in every consumer that skips the test.
