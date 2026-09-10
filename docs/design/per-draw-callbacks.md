# Per-draw callbacks

Status: DESIGN PROPOSAL (2026-09-10). Nothing is built. Anchor: bartcore
09a9c2a5; every citation below was read live against that tree.

## 1. Why now

The memory audit priced a large fit and found the engine is not the problem.
At its two reference cases
([Reference cases](memory-footprint.md#reference-cases)) the R-side prediction
array and its packaging copy are 3200 MB of a 3907 MB peak (n = 1e5, p = 20, C
= 4, S = 500) and 8000 MB of a 10960 MB peak (n = 1e6, p = 50, C = 1) - five
sixths of the peak in the first case and three quarters in the second, against
an engine share of 17 and 16 pct. One of the two copies came out in the
audit's step 7; the other is inherent to returning an n x S x C array at all.

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
the ANSWER needs instead: the reduction happens in place and the S axis never
exists.

## 2. What exists today

Three things, none of them the hook above.

**An internal R-level entry.**
[`bartcore_runWithCallback`](../../src/R_interface_bartcore.cpp) drives a run
with an R closure evaluated once per SWEEP. It refuses more than one chain
outright, so the closure runs inline on the main thread and blocks the run.
Its RNG contract is load-bearing: no `GetRNGstate`/`PutRNGstate` bracket,
because the chain's generator never touches R's stream while the closure may
draw from it, so R owns `.Random.seed` throughout. Its error contract is
`R_tryEval` - an error cannot longjmp across `Chain::run`'s C++ frames, so it
becomes a cooperative stop the bridge re-raises afterwards. It is a
CONDITIONING hook, firing at the top of each iteration, unthinned and
including burn-in, and it has had no caller since rbart_vi was retired.

**A retired flat-C entry.** `dbarts_sampler_setCallback`
(retired: [`dbarts_sampler_setCallback`](../../inst/include/dbarts/dbarts.h))
took `(userData, sampler, chainIndex, sweepIndex, isBurnIn)`, returned 0 to
stop, and was refused whenever chains would run on worker threads
([5. The sweep-boundary hooks that exist](bart-as-a-component.md#5-the-sweep-boundary-hooks-that-exist),
[6. C API and callbacks](public-surface.md#6-c-api-and-callbacks)). dec-B86
([docs/decisions.md](../decisions.md)) trimmed it out with the other entries
no consumer calls, narrowing the enabling-value rule for the C header alone to
"ships when a consumer or a named host design calls for it" - inclusion there
being irreversible. Re-opening it needs a named consumer.

**The engine hook both are built on.**
[`SweepCallback`](../../src/bartcore/chain.hpp) is a `std::function` the chain
calls before every sweep; [`Sampler::run`](../../src/bartcore/sampler.hpp)
forwards it only on the inline path, because with more than one worker chains
run on spawned `std::thread`s that must never call into R - progress lines
queue through [`QueuedProgressSink`](../../src/bartcore/sampler.hpp),
cancellation is a relaxed atomic, and since dec-B88 the caller blocks on a
condition variable. Writing results is separate:
[`storeSample`](../../src/bartcore/chain.hpp) writes each kept draw
contiguously into the caller-owned `Results` slabs that
[`allocChannel`](../../src/R_interface_bartcore.cpp) sized as n x L x S x C,
which is the natural firing point for a per-draw hook.

One rule governs everything below: sampler internals stay R-agnostic and the
bridge converts (dec-A50, superseded by dec-B85,
[docs/decisions.md](../decisions.md); it holds for conversion, not linkage -
the model header still reaches Rmath). The hook is therefore a plain C
function pointer over plain arrays, and no SEXP crosses into `chain.hpp`.

## 3. The C contract

A per-draw callback, distinct from the pre-sweep `SweepCallback` and not a
replacement for it.

```c
typedef struct dbarts_draw_t {
  size_t structSize;   /* library sets; read through DBARTS_HAS_FIELD */
  size_t chainIndex, drawIndex;  /* 0-based; drawIndex over KEPT draws */
  size_t numObservations, numTestObservations, numPredictors;
  size_t numReportedLocations;   /* L: 1, or K for multinomial */
  const double* train;      /* numObservations x L, or null */
  const double* test;       /* numTestObservations x L, or null */
  const uint32_t* varcount; /* numPredictors, or null */
  double sigma, k, dispersion, residualDf;  /* NaN where inapplicable */
} dbarts_draw;

typedef int (*dbarts_draw_callback)(void* context, const dbarts_draw* draw);
```

**When it fires, and for how long the pointers are good.** Once per SAVED
draw per chain, from the thread that owns that chain, immediately after
`storeSample` has settled every channel of that draw; thinned-away and burn-in
sweeps do not fire it (fork 5). Every pointer is valid for the duration of the
call and not one instruction longer (they address engine buffers the next
sweep overwrites), so a callback that wants a draw afterwards copies it. The
struct is size-first for the same reason `dbarts_results` is.

**No R API inside the callback, ever.** Three independent reasons, each
sufficient:

- *Worker threads.* R's evaluator, allocator and protection stack are
  single-threaded; `Rf_allocVector`, `Rf_eval`, `R_alloc` and even `PROTECT`
  from a spawned thread corrupt state the main thread owns.
- *GC.* Any R allocation may collect, which walks the protection stack. A
  worker allocating while the main thread sits inside `Sampler::run` can free
  objects nothing has protected on the worker's behalf.
- *longjmp.* `Rf_error` longjmps to a context the MAIN thread established.
  From a worker it unwinds the wrong stack; even on the main thread it skips
  every C++ destructor between raise and catch, the defect
  ["capi-unwind-protect"](../../TODO) already records for the existing entry.
  On POSIX, SIGINT is blocked in workers so R's interrupt handler cannot run
  there either.

The alternative, firing only on the inline path as `SweepCallback` does, is
rejected: it is unusable at the default `n.chains = 4`, exactly the
configuration whose result array the memory argument is about.

**Concurrency.** Calls for different chains may run CONCURRENTLY; calls
within one chain are ordered by `drawIndex`. The engine takes no lock, so a
callback touching shared state owns its own synchronization; the discipline
that needs none, and the one the example uses, is a write into a slot
addressed by `(chainIndex, drawIndex)`, disjoint by construction. Serializing
the calls behind an engine mutex would make every chain wait on the slowest
callback and turn a lock-free run into a contended one; rejected, and stated
in the header instead.

**What it may do.** Write into caller-owned memory allocated before the run;
accumulate (a running mean, a sum of squares, a held-out loss, a quantile
sketch); call plain C and C++ that allocates nothing R knows about. Not: call
R, throw (an exception escaping a worker body is `std::terminate`), block for
long, or retain a pointer from the struct.

**Return value: `int`, nonzero to stop, and errors.** Recommended over
`void`, because the callback cannot raise and a stop flag is its only escape;
the engine already carries the machinery (the per-chain `shouldCancel`
inline, the shared atomic on the worker path), so it costs one predictable
branch per saved draw. The caveat is that a cancelled run today discards its
results wholesale (`Sampler::run` returns true, the bridge raises "sampler
run interrupted"), so the flag means ABORT and not "I have enough", and the
bridge must distinguish "stopped by the callback" from "interrupted" exactly
as `bartcore_runWithCallback` distinguishes `closureStopped` from
`cancelled`. `void` is the alternative: simpler, and an infinite run then
cannot be ended from inside. A callback that fails records a status in its
own context, which the caller reads in R after the run; the engine reports
only that the callback stopped it, never why.

**Registration.** Two candidate sites, not exclusive:

- A per-run argument on the R bridge: the `.Call` run entry gains a function
  pointer and a context pointer, both read out of external pointers. No
  header change, no hash motion, no new frozen entry - and it is all the R
  surface needs.
- A sampler-level setter on the flat C header,
  `dbarts_sampler_setDrawCallback(sampler, fn, context)`, null `fn` clearing,
  re-opening the retired entry under dec-B86. Cost: one irreversible entry,
  one hash re-bake, and sampler state that persists across runs. Benefit: a
  `LinkingTo` consumer (stan4bart) registers once and drives many runs.

Recommendation (fork 2): ship the first now, hold the second for a named
consumer, and declare the TYPES in `dbarts.h` either way - an example that
hand-declares a matching struct is an unversioned copy of an ABI layout.

## 4. The R surface

`bart()`, `dbarts()` and `dbartsSampler$run` gain one argument, `callback`,
taking either of two things.

**(a) A C callback, non-blocking.** A list of two external pointers: one over
the function address, one over the caller's context. This is VD's default. The
run never re-enters R, so it composes with `n.threads > 1` and with every
chain running at once, and the R side's whole job is to check that both are
external pointers and pass the addresses down.

**(b) An R closure, blocking.** Evaluated on the main thread between sweeps,
which forces the run inline; single-chain is `bartcore_runWithCallback` today.
Multi-chain would need workers to enqueue each finished draw and the main
thread to drain the queue - so either buffering draws (reintroducing the array
this note is about) or blocking each worker until its draw is consumed,
serializing the chains behind R's evaluator. Recommendation (fork 1): do NOT
ship (b) in this arc. It costs a queue, a back-pressure policy, an
`R_UnwindProtect` around the run (["capi-unwind-protect"](../../TODO)) and a
second error contract, to serve a user who already gets the same effect from
`$run(0, 1)` in a loop - the pattern the component vignette teaches - at a
measured 0.256 ms per Gibbs step at n = 1000, T = 75
(benchmarks/baselines/bench-sampler-127f04ee.csv). The one thing that loop
cannot do is run chains concurrently, and neither can (b).

**How the storage opt-out composes.** Supplying a callback sets the
per-observation channels to not-kept unless the caller asks otherwise:
`keepTrainingFits` ([`keepTrainingFits`](../../R/bart.R)) defaults to FALSE
when `callback` is supplied, TRUE otherwise, and the test channel follows it.
That is what turns the hook into a memory lever; the alternative - never opt
out, let the user pass `keepTrainingFits = FALSE` themselves - is one line
simpler and makes the default callback fit pay for both. The mechanism is
small and belongs in the bridge: for a not-kept channel it allocates ONE
draw's worth per chain and points that chain's `Results` slab at it, so
`storeSample` writes somewhere the callback can read and the next draw
overwrites it. The write offset today is `sampleNum * n * L`, so `Results`
needs a per-draw stride the bridge can set to zero - an additive field on an
engine-internal struct, not an ABI event.

**What the returned object contains.** Everything a fit returns today except
the opted-out channels - sigma, varcount, k, the family scalars, the call, the
sampler under `keepSampler` - with `yhat.train` and `yhat.test` NULL, exactly
the shape `keepTrainingFits = FALSE` already produces and
[`packageBartResults`](../../R/bart.R) already handles. The families that
refuse `keepTrainingFits = FALSE` because their packaging reads the training
fits (ordinal probabilities, nbinom, hurdle) keep refusing it, and so refuse
the opt-out, not the callback.

## 5. The Rcpp example

The complete sketch the implementation slice will validate: a running
posterior mean over the training rows, the reduction that most often motivates
the full array. R side:

```r
Rcpp::sourceCpp("drawMean.cpp")   # makeMeanContext, meanCallbackPtr
acc <- numeric(nrow(x) * nChains)   # the ONLY per-observation allocation
ctx <- makeMeanContext(acc, nrow(x), nChains)
fit <- bart(x, y, n.chains = nChains,
            callback = list(fn = meanCallbackPtr(), context = ctx))
means <- rowMeans(matrix(acc, nrow(x), nChains))
```

C++ side, `drawMean.cpp`:

```cpp
#include <Rcpp.h>
#include <dbarts/dbarts.h>   // dbarts_draw; PKG_CPPFLAGS -I the include dir

struct MeanContext {
  double* out;          // REAL() of the R vector, taken BEFORE the run
  std::size_t n, numChains;
  std::size_t* counts;  // per chain, so no chain reads another's counter
};

extern "C" int meanCallback(void* context, const dbarts_draw* draw) {
  MeanContext* ctx = static_cast<MeanContext*>(context);
  if (draw->train == nullptr || draw->numObservations != ctx->n) return 1;
  double* out = ctx->out + draw->chainIndex * ctx->n;   // disjoint slice
  const double m = static_cast<double>(++ctx->counts[draw->chainIndex]);
  for (std::size_t i = 0; i < ctx->n; ++i)
    out[i] += (draw->train[i] - out[i]) / m;            // running mean
  return 0;
}

static void freeMeanContext(SEXP ptr) {
  MeanContext* ctx = static_cast<MeanContext*>(R_ExternalPtrAddr(ptr));
  if (ctx == nullptr) return;
  delete [] ctx->counts;
  delete ctx;
}

// [[Rcpp::export]]
SEXP meanCallbackPtr() {
  return R_MakeExternalPtr(reinterpret_cast<void*>(&meanCallback),
                           R_NilValue, R_NilValue);
}

// [[Rcpp::export]]
SEXP makeMeanContext(Rcpp::NumericVector acc, std::size_t n, std::size_t C) {
  MeanContext* ctx = new MeanContext{REAL(acc), n, C, new std::size_t[C]()};
  // acc rides the tag, so the array cannot be collected before the context
  SEXP ptr = PROTECT(R_MakeExternalPtr(ctx, R_NilValue, acc));
  R_RegisterCFinalizerEx(ptr, &freeMeanContext, TRUE);
  UNPROTECT(1);
  return ptr;
}
```

Writing a strided layout instead of a reduction is the same shape with one
extra context field and the body replaced by
`memcpy(ctx->out + (draw->chainIndex * ctx->numDraws + draw->drawIndex) * ctx->n, draw->train, ctx->n * sizeof(double))`.
That is exactly the layout the memory note's second array existed to produce,
built once, in place, by the consumer that wanted it.

**Lifetime rules, all of them the caller's.**

- `acc` must stay reachable and unmoved for the run. The external pointer's
  tag (above) ties it to the context's lifetime. `REAL(acc)` is read ONCE,
  before the run: R never moves a vector's data in place, but a reallocation
  on assignment makes a DIFFERENT vector and the callback keeps writing into
  the old one.
- The callback must not allocate through R, call any R entry point, or throw.
- It must be thread-safe across chains. The example is, writing a disjoint
  slice per chain with a per-chain counter. A SHARED accumulator needs
  atomics or per-chain partials combined afterwards - and atomics on doubles
  give up bitwise reproducibility.
- A callback that gets any of this wrong CRASHES the R session: no condition
  to catch, and a traceback pointing into the engine. That is the price of
  the non-blocking route and it belongs in the argument's documentation in
  those words.

**A second consumer.** stan4bart's embedding pattern - per-iteration
caller-owned buffers, one draw at a time, chains in separate single-chain
samplers - is this same shape reached from C rather than R, and it is what
would name the flat-header setter in fork 2.

**Rcpp is not a dbarts dependency.** The example lives in documentation and
its `Rcpp` use is the user's. The repo's precedent for a compiled example
under test is
["consumer source not installed"](../../inst/tinytest/test-capi.R), plain C
through `R CMD SHLIB`, self-gating on the toolchain and adding no Suggests. VD
asked for Rcpp specifically, so: show the Rcpp form in the vignette, test the
plain-C form (fork 4).

## 6. Memory consequence

Case 1 (n = 1e5, p = 20, T = 200, C = 4, S = 500), against
[Reference cases](memory-footprint.md#reference-cases):

| | today | with a callback, channels not kept |
| --- | --- | --- |
| engine | 658.3 MB | 658.3 MB |
| R, predictors and scalars | 49.0 MB | 49.0 MB |
| yhat.train | 3200.0 MB | 0 |
| per-draw scratch, 8*n*L*C | 0 | 3.2 MB |
| caller's accumulator, 8*n*C | 0 | 3.2 MB |
| peak | 3907.3 MB | 713.7 MB |

Case 2 (n = 1e6, p = 50, C = 1) has a much larger R-side predictor block
(1208.1 MB): 8000.0 MB of `yhat.train` becomes 8.0 MB of scratch plus 8.0 MB
of accumulator, and the peak falls from 10960.3 MB to 2976.3 MB. The engine's
share rises from 17 to 92 pct in case 1 and from 16 to 59 pct in case 2, which
is where the audit's remaining ranked items point. The honest comparison is
against the cheaper lever already in the manual: `keepTrainingFits = FALSE`
alone saves the same 3200 and 8000 MB for free, but DISCARDS the draws. The
callback is what buys the reduction back.

## 7. Threading interaction

The per-draw cost is an indirect call plus whatever the callback does, once
per SAVED draw - never during burn-in, and at `n.thin > 1` less often than per
sweep. A sweep at n = 1000, p = 10, T = 75 measures 0.174 ms
(benchmarks/baselines/bench-sampler-127f04ee.csv) and the example's callback
touches 8 KB, well under a microsecond: a few tenths of one percent at the
SMALLEST benchmarked size, and less at every larger one, where the sweep grows
with n*T and the callback with n. An estimate to be measured, not a claim -
the slice owes a bench-sampler scenario with a no-op callback.

The load-bearing point is the one VD named: a C callback is the only kind that
can run under `n.threads > 1` without serializing the chains. An R closure
runs on the main thread, so every worker that produced a draw waits on R's
evaluator, collapsing a four-chain run to one chain's throughput plus
queueing. The non-blocking route costs a branch; the blocking route costs the
parallelism.

## 8. Scope and sequencing

Minimum shippable surface, pre-release (VD raised the priority, 2026-09-10):

1. The engine hook: a `DrawCallback` fired from `Chain::run` after
   `storeSample` and forwarded by `Sampler::run` on BOTH paths, plus the
   per-draw stride on `bartcore::Results` that lets the bridge point an
   opted-out channel at a one-draw buffer.
2. The `dbarts_draw` and `dbarts_draw_callback` types in `dbarts.h`, and the
   bridge's per-run function-pointer argument.
3. The R `callback` argument on `bart()`, `dbarts()` and
   [`dbartsSampler$run`](../../man/dbartsSampler-class.Rd), the
   `keepTrainingFits` opt-out of section 4, and documentation carrying the
   crash warning verbatim.
4. The worked example, a seventh recipe in
   `vignettes/dbarts-as-a-component.Rmd` (fork 4).
5. Tests: a tinytest compiling a plain-C callback with `R CMD SHLIB` and
   checking its accumulated mean against the same fit's `yhat.train` mean,
   self-gating like ["consumer source not installed"](../../inst/tinytest/test-capi.R);
   a multi-chain run asserting per-chain call counts and ordering; a
   stop-flag run; and a tests/cpp test that the hook fires once per kept draw
   and never during burn-in.

Waits: the flat-header setter (fork 2, on a named consumer); the blocking R
closure (fork 1); any per-draw TREE handoff, which would hand a callback a
variable-length structure the engine does not hold outside `keepTrees` -
stated as not offered rather than half-offered.

Budget: roughly 450 lines across the engine hook, the bridge argument, the R
surface and its Rd, the vignette recipe and the tests.

RNG class NEUTRAL: the hook consumes no generator draw and mutates no sampled
state, so a run with no callback registered stays bitwise identical, the
property `SweepCallback` and `shouldCancel` already carry. Gates owed:
tests/cpp; the full tinytest suite; all three equivalence harnesses IDENTICAL
on every scenario; and, the hook sitting on the draw path,
`bench-sampler.R compare` on a quiet machine with a no-op-callback scenario.

ABI and sister packages. Declaring the two types moves nothing a consumer
links against; adding the setter is ADDITIVE (a name-looked-up entry) and
re-bakes `DBARTS_C_API_HASH`. Under dec-B111 the version pair is the guard and
CI asserts a changed hash comes with a minor bump, so pre-release the
constants stay at 1/0 with one re-bake and a post-release addition is a bump
and a re-bake together. Neither breaks a consumer: stan4bart and treatSens
call neither type nor entry.

## 9. Open forks for VD

1. **The blocking R-closure variant: ship or not?** Recommend NOT in this
   arc: no multi-chain without either a draw buffer (the array we are
   removing) or serialized workers, plus `R_UnwindProtect` and a second error
   contract, to serve a user `$run(0, 1)` in a loop already serves. Reopen it
   if the crash-on-mistake cost of the C route is what users actually hit.
2. **Registration site: sampler setter, run argument, or both?** Recommend
   the run argument plus the TYPES in `dbarts.h` now, holding
   `dbarts_sampler_setDrawCallback` until stan4bart names it - dec-B86 makes
   header inclusion irreversible, and the R layer does not need the entry.
3. **The stop-flag return.** Recommend `int`, nonzero to stop, the bridge
   distinguishing "stopped by callback" from "interrupted". Alternative
   `void`: one less branch and one less contract, and a callback that cannot
   end a run it knows has gone wrong.
4. **The example's home.** Recommend a seventh recipe in
   `vignettes/dbarts-as-a-component.Rmd` (Rcpp form, unevaluated) plus a
   plain-C copy under `inst/tinytest/` the suite compiles. Alternatives:
   `inst/examples/` alone (nothing tests it), or a companion package (a
   second release to keep in step, for one file).
5. **Does it also fire during burn-in?** Recommend NO: a burn-in sweep has no
   slot and no draw index, the pre-sweep hook covers anyone who wants every
   sweep, and firing makes `drawIndex` mean two things. Alternative: an
   `isBurnIn` flag and a second counter - one more field, and a running mean
   silently wrong for anyone who ignores it.
6. **The family scalars by value or by pointer?** Recommend by value with NaN
   for inapplicable: one cache line, cannot dangle, no null check per field.
   Alternative, null pointers, matches `dbarts_results` and costs a check per
   read.
7. **Does the opt-out extend to varcount and the scalars?** Recommend NO -
   p*F*S*C and S*C are kilobytes at both reference cases, and dropping them
   breaks every diagnostic for no gain. Alternative: a `keep =` character
   vector naming the channels retained, a cleaner surface and a wider one to
   freeze.
