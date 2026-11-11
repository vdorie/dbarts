# threaded-chains: stan4bart's chains on dbarts' threads, one sampler for every chain

Status: PLANNED 2026-11-11, revised after its blind critique (round 1). Route (b), the cores default and the
user's R callback ruled 2026-11-11, and per-chain values on a copy ruled (a) the same day; the compiled
callback planned at the maintainer's request (S4). One conditional call remains ([Open calls](#open-calls)).

agent: one opus implementer per repository line (dbarts D1 and D2, stan4bart S1), effort high; one opus
reviewer who runs mutants on D1, D2 and S1; a sonnet implementer and reviewer for D3, S2, S3, S4 and the
bartCause slice B1, with an opus read of S4's header and its worker-side C++ before commit.
rng: NEUTRAL for dbarts: with no hook registered every run takes the path it takes today and no baseline,
snapshot or exact gate moves; D3 changes thread counts only under R CMD check, and draws do not depend on the
thread count. SHIFTING for chains 2 and on of a stan4bart fit (each chain's seeds derive from the seed and
the chain number, the same at every thread count; chain 1 is the one-chain fit); a one-chain fit is expected to keep its draws bit for bit
([Draws at a fixed seed](#draws-at-a-fixed-seed)). Not posterior-changing anywhere: the Gibbs order inside a
chain is today's.
window: before the 1.0-0 release (scheduled by the maintainer 2026-11-11). The dbarts branch stays open until
S1's gates pass against it, then merges once.
budget: dbarts ~1610 lines planned (forecast 2900 to 3800); stan4bart ~1590 (forecast 2500 to 3300); bartCause
~50. See [Sizes](#sizes).

## Goal

A stan4bart fit holds one dbarts sampler with every chain in it. Each chain carries its own offset and sigma,
and stan4bart makes one dbarts run per chunk of sweeps. Before each group of sweeps, the dbarts worker thread
that owns chain c calls a compiled stan4bart hook. The hook draws that chain's parametric part and sets that
chain's offset and sigma through a view of that chain alone, which never calls R and never raises. When the
group's draw settles, the same worker calls stan4bart's draw callback, which records the draw and hands it to
the parametric part. dbarts keeps all thread machinery. Chains run on min(cores, chains) threads under the
ruled default. Seeded draws are bitwise identical at every thread count. The PSOCK and FORK code goes. A user's
R callback runs on the main thread between chunks; a compiled one runs on the workers at every draw.

## Context

The ruling (2026-11-11, recorded in the backlog rulings under "Threaded chains route"): "Ok, (b) it is, where
it is possible to do so. If we need to do (a) in some cases, we can." Route (b) is expert A's route 3 in the
three memos and synthesis (scratch/pcthreads/, not tracked). Route (a), one sampler per chain with each whole
chain on a thread of the consumer's, needs off-main-thread entries in dbarts. It is built only where (b)
cannot serve, for example treatSens with `trt.model = bart()`, where two samplers run in lockstep. It is
a later, separate item, not built here ([Later items](#later-items)).

The cores default, ruled 2026-11-11 ("Do (b)"): stan4bart's `cores` defaults to `getOption("mc.cores")` when
that option is set, else to `min(dbarts::guessNumCores(), chains)`, and runs on threads. bartCause's own rule
for when to give stan4bart cores goes in a bartCause slice (B1). Tests and examples cap at 2. NEWS says the
default changed from 0.0-13, which was `getOption("mc.cores", 1L)` on processes. S1 confirms that by running
0.0-13.1; this plan read the tag.

The user's R callback, ruled 2026-11-11: "(a), although if there's a way to support C++ callbacks through
something like Rcpp, ideally compiled in an R script, that would be ideal." So the R callback is evaluated on
the main thread between chunks over buffered draws (S3). The compiled callback is S4. It was answered feasible
and is planned here as its own slice after the threading slices; the maintainer may still object to it.

Earlier rulings this rests on: dec-B419 (the maintainer's question whether one sampler could hold every
chain); dec-B437 (stan4bart's response range is fixed at creation and shared, so chains have no transform to
agree on and need not meet inside a run); dec-A197 (today's restore, one sampler per chain, which S2
replaces).

Code read for this plan: dbarts main 6a34f43b1; stan4bart main 89af617; bartCause main ec1ee3e; treatSens main
268b2db. The blind critique checked the hook position, the absence of races, the one-chain seed order and the
posterior's identity against the code (scratch/pcthreads/crit/critique.md, not tracked).

What exists:
- [`Sampler::run`](../../src/bartcore/sampler.hpp) starts min(numThreads, numChains) workers at each call
  through [`spawnWorkers`](../../src/bartcore/sampler.hpp), with SIGINT masked on POSIX, a cancel flag the main
  thread sets after its 100 ms interrupt poll, progress queued through
  [`QueuedProgressSink`](../../src/bartcore/sampler.hpp), and the first exception a chain throws rethrown after
  the join. It joins every worker before it returns, so chains meet only at the end of a run.
- The per-sweep hook, [`SweepCallback`](../../src/bartcore/chain.hpp), fires in
  [`Chain::runSweeps`](../../src/bartcore/chain.hpp) before every sweep, after the cancel check and before the
  augmentation draw, so a change made there reaches the sweep's latents. `Sampler::run` passes it only on the
  inline path. Its one installer, [`bartcore_runWithCallback`](../../src/R_interface_bartcore.cpp), evaluates
  an R closure and refuses more than one chain. The flat API once had a hook of this kind, which it refused with
  more than one thread (retired: [`dbarts_sampler_setCallback`](../../inst/include/dbarts/dbarts.h); see
  [5. The sweep-boundary hooks that exist](../design/bart-as-a-component.md#5-the-sweep-boundary-hooks-that-exist)).
- The per-draw observer [`dbarts_draw_callback`](../../inst/include/dbarts/dbarts.h), registered by
  [`dbarts_sampler_setDrawCallback`](../../inst/include/dbarts/dbarts.h), already runs on workers with no lock,
  under a no-R contract. It fires once per saved draw on the chain's own thread, right after the draw settles,
  with every channel the fit carries (train, test, sigma, k, varcount and the rest) whether or not the run
  stores them, and carries no latents. A run with null results still fires it.
- [`Chain::setOffset`](../../src/bartcore/chain.hpp) and [`Chain::setSigma`](../../src/bartcore/chain.hpp)
  touch only that chain when the transform is held: the response model, its working response and sigma are
  per chain; the saved-draw restatement returns at once when no slot is named or the scale did not move, and a
  variance forest's re-anchor returns when its scale did not move. The sharing is in
  [`Sampler::setOffset`](../../src/bartcore/sampler.hpp), which fans one vector out to every chain, and in the
  holder's one `ownedOffset` ([`dbarts_sampler_t`](../../src/R_interface_bartcore_common.hpp)).
- A loop of one-sweep runs is bitwise identical to one batched run
  ([1. The per-sweep driver loop is bitwise identical to a batched run](../design/bart-as-a-component.md#1-the-per-sweep-driver-loop-is-bitwise-identical-to-a-batched-run)).
  That identity is what lets chunked runs keep stan4bart's draws.
- [`dbarts_sampler_run`](../../src/C_interface.cpp) raises R conditions and polls R through
  [`userInterrupted`](../../src/R_interface_bartcore.cpp); the setters raise on a non-finite value, and
  `dbarts_sampler_setOffset` raises, at any `updateScale`, on a multi-forest sampler whose coupling admits no
  response mutation ([`refuseMultiForestResponseMutation`](../../src/R_interface_bartcore.cpp)). These stay
  main-thread entries; nothing here changes them.
- The token ([`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)) folds the list's stringized signatures,
  the ABI enums and each crossing struct's layout by field name and offset. A struct member's function-pointer
  type is outside that fold, and so is the `dbarts_draw_callback` typedef's parameter list.
- The thread defaults of [`bart()`](../../R/bart.R), [`rbart_vi()`](../../R/rbart.R),
  [`dbartsControl()`](../../R/dbarts.R) and [`xbart()`](../../R/xbart.R) read
  [`guessNumCores`](../../R/guessNumCores.R), which ignores `_R_CHECK_LIMIT_CORES_`.

stan4bart today, read at 89af617: `stan4bart_fit_worker` builds one `dbartsSampler` per chain
(`n.chains = 1`, `n.threads = 1`, `n.thin = skip.bart`) and calls `C_stan4bart_run` once per phase with result
type `"both"`; the `"stan"` and `"bart"` types are never passed. Per iteration that call:
1. runs the parametric step: the ridge move in sampling, then `skip` WALNUTS transitions;
2. reads the parametric mean into the BART offset, sets sigma and the offset;
3. runs one dbarts run of one kept draw;
4. stores the training fits less that draw's offset (`bart_train`), which become the parametric part's
   offset, and hands a binary fit's latents to the parametric part;
5. evaluates the user's R callback.

Under `keep_fits = FALSE` it keeps one draw. The parametric code draws from its own `std::mt19937_64` and
reports failure by C++ exceptions; it calls R only through stan4bart's data readers at creation and a debugging
print in `double_writer.hpp`. `cores > 1` starts PSOCK (else FORK) workers, seeds each chain from `sample.int`
after `set.seed(seed)`, and falls back to serial with a warning when the cluster fails to start or a chain
fails. `cores = 1` draws every chain's seeds from one running stream, so seeded draws depend on `cores` (run
2026-11-11, synthesis fact 6). `.Random.seed` is restored after a fit only when `seed` was given.
`C_stan4bart_predictBART` passes 0 threads, meaning the sampler's own count.

## Design

### The hook and the chain view (dbarts.h)

```c
#define DBARTS_CHAIN_VIEW_OPS(X) \
  X(int, setOffset, (dbarts_chain_view* view, const double* offset)) \
  X(int, setSigma, (dbarts_chain_view* view, double sigma))
struct dbarts_chain_view_t {
  size_t structSize;       /* the library's sizeof, as dbarts_draw */
  size_t chainIndex;
  size_t numObservations;
  void* internal;          /* the library's; never read by the caller */
  /* one function-pointer member per DBARTS_CHAIN_VIEW_OPS entry */
};
typedef int (*dbarts_sweep_hook)(void* context, dbarts_chain_view* chain,
                                 size_t sweepIndex, int isBurnIn);
int dbarts_sampler_setSweepHook(dbarts_sampler* sampler,
                                dbarts_sweep_hook fn, void* context);
```

The sketch fixes the shape; the implementer settles names in review. One entry joins `DBARTS_C_API_LIST`
(registration). The view's operations are function pointers that the library fills, valid only during the
call they are passed to (call A4). The members are generated from `DBARTS_CHAIN_VIEW_OPS`, whose full text,
return type and parameter list included, is stringized into the token. In the same change the token gains the
text of both hook typedefs, `dbarts_draw_callback` and `dbarts_sweep_hook`, so a changed parameter type moves
it (call A13).

- Registration copies (fn, context) into the holder on the copy-on-set rule and allocates one n-double offset
  buffer per chain, not seeded: a chain draws against the sampler-wide offset until its hook first writes its
  own through the view, which copies all n values. It is a capability
  status: 0, touching nothing, on a sampler where a chain's offset can never be written: a fixed response
  conduit, or a multi-forest coupling that admits no response mutation, the cases where
  `dbarts_sampler_setOffset` returns 0 or raises whatever its arguments. 1 otherwise. A null fn clears the hook
  and returns every chain to the shared offset. Registration runs on the main thread and may allocate; no
  chain operation does. The hook fires on runs through `dbarts_sampler_run` only, as the draw callback does,
  and is not part of a saved state or a copy (call A9).
- The hook fires before every sweep of every chain, on the thread that owns the chain: a worker, or the calling
  thread when the run is inline. It gets the 0-based sweep index within this run and the burn-in flag.
  0 continues. Nonzero stops the run, a positive value as a deliberate stop and a negative one as the hook's
  failure, the sign rule both callbacks follow (dec-B581): every chain stops at its next sweep boundary, the run's results are
  discarded and its saved trees dropped, exactly as a draw callback's stop. `dbarts_sampler_run` then returns
  normally, and the caller reads why from its own context (call A5).
- The contract the header states: no R call of any kind (allocation, printing, errors, RNG, interrupt checks)
  on any path, inline included, so one rule holds at every thread count. The hook may touch only its own
  chain, through the view, and must not call any `dbarts_sampler_*` entry on this sampler; a debug build
  asserts that (call 12). No C++ exception may escape either hook (call A6).
- `setOffset(view, offset)`: n finite values copied into that chain's buffer, with the transform held, never
  re-derived (no updateScale argument). Returns `DBARTS_OK` (0) on a write and `DBARTS_INVALID_ARGUMENT` (-2),
  touching nothing, on a non-finite value. Registration has already ruled out a sampler where the write is
  never legal. `setSigma(view, sigma)`: `DBARTS_OK` on a write; `DBARTS_UNSUPPORTED` (-1) where the model fixes
  sigma, as `dbarts_sampler_setSigma` decides; `DBARTS_INVALID_ARGUMENT` on a value that is not finite and
  positive. The names, values and words are the host-neutral plan's statuses, from its three-argument
  `DBARTS_STATUS_LIST` (call A7).
- `dbarts_sampler_setDrawCallback` takes the channels the callback reads, `unsigned int channels`: one
  `dbarts_draw_channel` flag per pointer channel of `dbarts_draw` (the variance pair one flag), or
  `DBARTS_DRAW_ALL`. A run computes the channels it stores and the named ones; one neither stored nor named
  is NULL and not computed, the log-likelihood pass included. The flags fold into the token, and the R route
  passes `DBARTS_DRAW_ALL` (ruled 2026-11-10; call A16).
- `dbarts_draw` gains `latents`: the chain's n augmentation values this draw's trees were drawn from, null
  where the family augments none (the values `dbarts_sampler_getLatents` reports after a run). It is a
  pre-release field, so it goes above the 1.0-0 boundary marker, not below it, and moves no version constant.
  It sits among the pointers, after `varcount` and before `sigma`, and joins
  [`DBARTS_DRAW_POINTER_FIELDS`](../../src/C_interface.cpp). The token folds pointer fields by offset in
  pointer units, which is one number everywhere only ahead of the four doubles. The token re-bakes.
  "The draw struct carries none" changes in the header and in architecture's
  [Callbacks and errors](../architecture.md#callbacks-and-errors) (call A8).
- `dbarts_sampler_run`'s comment "A caller makes its changes between runs: this API has no hook inside a
  sweep" becomes: between runs, or inside a sweep hook through the chain view, which reaches that chain alone.
  The threading notes at the top of the header gain the one exception: the view's operations are called from
  whatever thread runs the hook.

The view, not the handle, is what the hook receives. A hook cannot reach another chain or a sampler-wide
setter. Its operations are never resolved through `DBARTS_USE_STUBS`, whose first-call
`R_GetCCallable`, and the handshake `Rf_error` beside it, would otherwise run on a worker.

### Engine (D1)

- `SweepCallback` becomes one plain hook type, a function pointer plus context, called with a `ChainView&`: a
  small non-template interface that `Chain` implements over its own state. It is legal on both paths.
  [`Sampler::run`](../../src/bartcore/sampler.hpp) passes it to the workers, as it passes the draw relay
  today. The inline-only restriction leaves the engine. The bridge entry that installs an R closure keeps its
  own refusal of more than one chain and adapts the closure into the hook type through its context. One hook
  type serves both, so there is no second call site (call A1).
- `ChainView` offers two operations. One sets the offset with the transform held: a `Chain` method that calls
  `response_->setOffset(ptr, false, ...)` and names no saved slot, so it never touches the shared tree store.
  The other sets sigma.
- `DrawInfo` gains the chain's latents pointer, which the bridge's draw fill passes on.
- The hook's stop takes the path the draw relay's stop takes: a set-once relaxed flag read after the join, so
  `Sampler::run` can tell a hook's stop from an interrupt. A hook's exception is caught in the worker body that
  already catches chain exceptions, and at the call on the inline path, as the draw callback's is.
- `running_` becomes atomic so the debug assertion (call 12) can read it from a worker. No draw-path code
  changes, but the hook's type and its passing change inside `Chain::runSweeps` and `Sampler::run`. Those are
  on the hot path, so D1 takes a bench-sampler compare (maintainer-run, quiet machine) against main.

### Bridge (D2)

- The holder gains `sweepHook`, beside `drawHook`: the shipped hook (fn, context) and its
  `std::vector<std::vector<double>> chainOffsets`, nested, as `ownedForestWeights` is, so a pointer a chain
  borrows survives later growth. The buffers are sized at registration, empty otherwise.
- The view's operations live in the shared bridge header, so tests/cpp reaches them, with no R API inside
  (no `R_finite`; `std::isfinite`).
  The bridge fills one view per chain per run on the worker's stack.
- A sampler-wide `dbarts_sampler_setOffset` (or the R `$setOffset`) after chains hold their own offsets puts
  every chain back on the shared buffer, as it does today. A sampler-wide re-derivation of the transform
  (`updateScale` true, or [`Sampler::rederiveScale`](../../src/bartcore/sampler.hpp), which reads chain 0's
  offset) while chains hold different offsets is refused with an error naming the cause. The rule it would
  need, one transform read from which chain's offset, is not one a user means (call A10).
- `bartcore_runWithCallback` moves onto the single hook type unchanged in behavior.
- The token re-bakes. The version constants
  ([`DBARTS_C_API_MINOR`](../../inst/include/dbarts/dbarts.h)) do not move before the first release, by the
  header's own rule. If this lands after 1.0-0 ships, the minor moves to 1 and stan4bart's floor with it (call
  A11).

### Default thread count under R CMD check (D3)

dbarts exports one small function, `defaultNumThreads(n.tasks)`. It returns `min(guessNumCores(), n.tasks)`,
capped at 2 when `_R_CHECK_LIMIT_CORES_` is set and not `"false"` in any case, and NA when the cores cannot be counted, so
dec-B306's one-thread fallback stays where it is. The thread defaults of `bart()`, `bart2()`, `rbart_vi()` and
`dbartsControl()` (`n.tasks = n.chains`) and `xbart()` (`n.tasks = Inf`) call it. stan4bart's default calls it
through dbarts' exported surface, not an internal: `getOption("mc.cores")` when set, else
`defaultNumThreads(chains)`, with NA read as 1. `guessNumCores` keeps its meaning, the cores on the machine. An
explicit count, as an argument or as `mc.cores`, is never capped. The cap applies to the default alone, since
parallel's own check is what refuses an explicit count above 2 under the variable (call 14). R-surface tier,
one sonnet review.

### stan4bart (S1 to S4)

- S1, the loop. `stan4bart_fit` builds one `dbartsSampler` with `n.chains = chains` and `n.threads =
  min(cores, chains)` and installs the shared range once. `C_stan4bart_create` reads the handle and holds C
  parametric samplers, each with its own buffers (BART offset, parametric offset, latents, parameter writer)
  and seed. It runs the one creation run (every chain once, no hook) and does each chain's first hand-off on
  the main thread, as today. Then it registers two compiled functions over one context that holds the C
  per-chain records, the raw pointers of the result arrays (taken on the main thread, since `REAL()` is an R
  call) and an error slot per chain:
  - the sweep hook, `noexcept`, acting only when `sweepIndex % skip.bart == 0`. It runs the parametric step
    (the ridge move in sampling, then `skip` transitions, one in warmup), writes the chain's parameter row, and
    sets the chain's offset (the parametric mean, plus a user offset by `offset_type`) and, for a continuous
    fit, sigma. These are steps 1 and 2 above.
  - the draw callback, `noexcept`, at each kept draw, which ends each group. It writes train less this chain's
    offset into `bart_train` at (chain, chunk base plus draw index), writes test, sigma, k and varcount, and
    hands the draw to the parametric part: the fits as its offset, and a binary fit's latents from the draw.
    These are step 4 above. At `keep_fits = FALSE` it writes the one slot, as today. Under S4 it also calls a
    compiled user callback.

  Every step of today's iteration runs in today's order inside one group, on the chain's thread, so a run
  boundary is not a special case: the warmup-to-sampling boundary, where `C_stan4bart_disengageAdaptation`
  runs, finds every chain as today's loop leaves it. dbarts runs with null result channels, so it allocates
  nothing per chunk beyond its own scratch for the draw callback (call A2). A WALNUTS exception, or a -1 from
  a view operation (a diverged parametric step's NaN), is caught; its message goes into the chain's slot and
  the function returns nonzero. Everything the two functions touch is owned by the external-pointer object,
  never by a `.Call` frame, so an interrupt's jump strands nothing. Runs keep `numBurnIn = 0`, and the Gibbs
  coupling depends on it: the draw callback never sees a sweep dbarts discards as burn-in, so under dbarts
  burn-in the BART draws would silently stop reaching the parametric part (and the warmup trace would stop
  recording). stan4bart's warmup is its own phase, recorded, never dbarts' burn-in. `C_stan4bart_run` loses its dead `"stan"` and `"bart"`
  result types. The adaptation, eval-count and initial-summary entries loop over the chains.
- Chunks: a phase runs in one dbarts run unless something stops it sooner. With verbose above 1, chunks end
  every `refresh` iterations, where `iter i / N` prints; with an R callback, every 100 iterations or at
  `refresh` (S3). Chunk length cannot move draws; an internal option (`stan4bart.chunk`, unexported,
  test-only) forces a length, so the gate can compare 1 against a whole phase. A compiled callback needs no
  chunks.
- Progress: printed by the main thread between chunks. Every chain ends a chunk at the same iteration, so one
  `iter i / N` line describes them all; at one chain the lines are today's (call 5). Verbose 1's "fitting
  chain" lines become one line naming the chain count and thread count. The "starting multithreaded fit,
  further output silenced" line goes.
- Errors: a chain's failure stops every chain at its next sweep boundary and raises one R error after the join,
  naming the chain and its message. It replaces today's warning and silent serial rerun (call 6). An interrupt
  is dbarts': polled every 100 ms, it stops the workers at their next sweep boundary and raises R's interrupt
  condition after the join. Latency is at most one sweep plus one parametric step. dbarts' slow-count warning
  comes once per fit rather than once per chain's sampler: a NEWS line.
- Seeds (call 2): each chain's seeds, parametric and BART, derive from the user's seed and the chain number
  alone, so chain 1 of a multi-chain fit equals the one-chain fit at the same seed, as in dbarts' `bart()`
  and today's stan4bart. After `set.seed(seed)` (or R's stream as it stands when `seed` is NA), stan4bart
  draws a fixed number of values whatever `chains` is, chain 1's parametric seed first, as today's one-chain
  order does, and derives chain c's parametric seed from them and c; the sampler's chain seeds come from its
  own generator given one of them, which by dbarts' rule gives each chain a seed set by that value and the
  chain's index. S1 settles the derivation within that rule. `.Random.seed` is put back afterwards only when
  `seed` was given, as now; an unseeded fit advances R's stream, so two in a row differ. Every generator is a
  chain's own and the hook touches only its chain, so draws do not depend on `cores`.
- The cluster code goes: `makeCluster`, `clusterMap`, `clusterExport`, the FORK fallback, both fallback
  warnings, the per-worker seeding and `parallel` from Imports if nothing else uses it.
- The cores default is as ruled, through `defaultNumThreads` (D3). `cores` keeps its name and is validated as
  one positive whole number; `min(cores, chains)` threads run.
- S1 runs stan4bart 0.0-13.1 against dbarts 0.9-34, in private libraries, to confirm what NEWS will say: the
  old cores default, and that `bart_args` seeding gave every chain one BART stream.
- S1's help and NEWS: `cores` describes threads, the new default, the 2-thread cap under R CMD check and the
  draw change. `bart_args = list(seed = )` giving each chain its own BART stream is listed as a correctness
  fix ([Draws at a fixed seed](#draws-at-a-fixed-seed)). NEWS diffs against 0.0-13.
- S2, restore, predict, extract and save. `sampler.bart` becomes one `dbartsSampler`, and `state.bart$state`
  the one sampler's state. `restoreBartSampler` builds one sampler of C chains, installs the one range once
  and sets the state. `bart_pointer_is_live`, `getBartSampler` and `predict_bart_chains` (one
  `C_stan4bart_predictBART` call returning draws x chains) follow. `extract(type = "trees")` reads one
  sampler's chains by `chainNums`. A saved fit holds one sampler's data copy instead of C. Predictions and
  recomputed blocks thread across chains. They take their thread count from the session's default rule when
  they run, not from the fit's saved control (call A14). Help and NEWS say so, and the tests cap it at 2.
  `restoreBartSampler` keeps a branch that installs each chain's own `fit.scale` for fits saved before the one
  range. S2 first establishes which shipped version, if any, wrote such fits. A 0.0-13.1 fit, written through
  `C_stan4bart_exportBARTState` against dbarts 0.9-34, is loaded on main and on the branch (call A12).
- S3, the R callback, as ruled. While `callback` is given, chunks end every 100 iterations, or at `refresh`
  when verbose. After each run the main thread evaluates the callback chain by chain, iteration by iteration,
  over the chunk's buffered draws: train and test fits and the parameter row, from the fit's arrays when
  `keep_fits` is TRUE, else from a chunk buffer of 100 x chains draws. The results keep today's shape and
  per-chain order. Only the interleaving across chains differs: serial fits used to run chain 1 to the end
  first.
- S4, the compiled callback (the maintainer's request), after S3:
  - stan4bart ships `inst/include/stan4bart/callback.h`, plain C, versioned as dbarts.h versions its structs.
    It holds a `STAN4BART_CALLBACK_MAJOR` and `_MINOR`, and a `stan4bart_draw` struct whose leading
    `structSize` stan4bart fills, with fields appended below a marked boundary: chain, iteration, the warmup
    flag, the three lengths, `yhatTrain`, `yhatTest`, the parameter row and the parameter names. The
    signature is `typedef int (*stan4bart_callback)(void* context, const stan4bart_draw* draw, double* result,
    size_t resultLength)`. Under R headers the header also gives a macro that builds the external pointer with
    `R_MakeExternalPtrFn` and tags it with the major version, since a function pointer is not a `void*`.
    The pointer is valid only while the user's `sourceCpp` library stays loaded, so the help says to create
    it right before the fit and not to keep it across sessions.
  - The user writes the function in an R script with `Rcpp::sourceCpp(code = )` and
    `// [[Rcpp::depends(stan4bart)]]`, exports a function returning that pointer, and passes `callback =
    compiled_callback(ptr, length = 1)`. stan4bart refuses a pointer whose tag's major differs.
  - It runs inside stan4bart's draw callback, on the chain's worker, at every draw, after the record and the
    hand-off. It sees what the R callback sees, plus chain and iteration, and writes `length` values into a
    slot stan4bart gives it. The results land in `fit$callback` in the R callback's shape.
  - Contract: no R and no Rcpp R types inside (RcppParallel's rule). Exceptions follow call A6. stan4bart
    catches anything the call throws, and a nonzero return stops the fit with one R error naming the chain and
    iteration.
  - A test-only callback compiled into stan4bart (an unexported entry returning its pointer) drives the tests,
    so CRAN tests compile nothing. A worked example in the help is under `\dontrun{}`; it needs a compiler and
    Rcpp, which becomes Suggests.
  - A helper taking only the body and writing the boilerplate (signature, `noexcept` wrapper, the export) is
    an agent call (call A15).
  - No dbarts change: everything the callback needs reaches stan4bart's draw callback already.

### bartCause (B1)

`R/argParse.R` gives stan4bart `cores = min(n.threads, chains)` only when the fit is long enough
(`chains * iter >= 4000`), because a cluster costs a fixed startup. When it does, it also injects a seed so a
clustered fit stays reproducible under `set.seed`. Under threads neither reason holds. B1 passes `cores =
n.threads` when the user gave `n.threads`, and otherwise leaves `cores` unset so stan4bart's default applies,
capped under R CMD check by D3. The injected seed goes, since stan4bart now seeds from R's stream at every core
count. bartCause tests that fit stan4bart with a seed are re-recorded, and its tests and examples pass at
most 2.

## Draws at a fixed seed

- dbarts: none move. No hook, no new path.
- stan4bart, one chain: today's serial order is "parametric seed by `sample.int(max, 1)`, then the sampler's
  one chain seed", and the new order at C = 1 is the same. Today's iteration order is kept step for step
  inside a group. So a one-chain fit is expected to reproduce stan4bart main 89af617 bit for bit, and S1 gates
  on it. If something else moves it, the implementer names the cause and the reviewer judges. The equality
  would not come back by accident.
- stan4bart, more than one chain: chains 2 to C move once, at `cores = 1` and at `cores > 1` alike; chain 1
  equals the one-chain fit at the same seed. From then on, `cores` = 1, 2 and 4 give bitwise-identical fits.
- A correctness fix, by the code: today `bart_args = list(seed = s)` gives every chain's own sampler the same
  control seed, so every chain draws the same BART stream. The chains were not independent. In one sampler the
  seed generator gives each chain its own seed, which is what `seed` means in dbarts (a one-chain run with
  seed s reproduces chain 1 of any run with seed s). S1 confirms the old behavior by a run of 0.0-13.1.
- Snapshots: stan4bart has no reference-build snapshot files. Its seeded expectations that hold values
  (test-05-rng holds identities, not values) are re-recorded by replaying whole files. bartCause's
  stan4bart-backed seeded tests are re-recorded in B1.

## Gates

Tiers by [Process by risk](README.md#process-by-risk). D1 and D2 change the engine: plan with blind critique,
an opus review that runs mutants, the RNG class's battery for NEUTRAL, the full tinytest suite, sanitizers on
tests/cpp and on the R-loaded path (the chain view is new code reachable from a `.Call`), `R CMD check`, the
arm64 reproducibility snapshots, which must not move, and bench-sampler.R compare for the hot-path change. S1
changes stan4bart's draws: the same tier on stan4bart's side. D3, S2, S3, S4 and B1 are R-surface tier.

| gate | slice | expected |
|---|---|---|
| tests/cpp: a hook that sets each chain's offset and sigma from its own last draw (a toy Gibbs step) through the view, against m one-sweep runs with the same changes made by `Chain::setOffset` and `Chain::setSigma` between them | D1 | bitwise equal, every channel, at 1, 2 and 4 threads |
| tests/cpp: the same sampler at 1, 2, 3 and 4 threads (3 leaves one worker with two chains) | D1 | bitwise equal |
| tests/cpp: a hook stop on chain 1 of 4 at sweep k; a hook that throws, inline and on a worker | D1 | the run reports the stop (not an interrupt); the store's cursor and counts as after a draw callback's stop; the throw comes back after the join |
| tests/cpp: chain c's offset write leaves every other chain's working response and fits bit for bit | D1 | equal |
| tests/cpp: `DrawInfo` latents for probit equal the sampler's latents after a one-sweep run | D1 | equal |
| bench-sampler.R compare, D1 against main, no hook (maintainer-run) | D1 | within the compare's noise |
| the draw callback's own cost: a one-chain probit sampler at n = 5000, sweeps timed with no callback and with a no-op one naming train alone, same machine, alternating | D2 | recorded (call A16) |
| consumer: the recording callback registered for each channel flag alone, results NULL; every array stored and no channel named; every channel named | D2 | the named channel where the fit carries it and no other, its draw the every-channel one; the stored channels; today's channels |
| test-capi.R plus capi/consumer.c: a compiled hook and draw callback through the flat API doing the toy Gibbs step; chunked against one-sweep flat runs with `dbarts_sampler_setOffset` or setSigma between, at one chain (where the two routes are the same model) | D2 | bitwise equal |
| the same, 4 chains at 1 and 2 threads | D2 | bitwise equal |
| consumer: NaN offset, sigma 0, a pinned-sigma sampler (probit); registration on a fixed-conduit sampler | D2 | -2, -2, -1 with the chain unchanged; registration 0 and nothing registered; nothing raised |
| a C++ consumer compiled in test-capi.R whose hook and draw callback throw, inline and at 2 threads, on all three CI systems | D2 | the safety net: an R error carrying the what(), both hooks; the sampler runs on (call A6) |
| consumer: a hook returning 1 and -1, a draw callback returning -1; an interrupt injected through dbarts' count hooks at 2 chains and 2 threads, in one-sweep runs so a poll lands inside every run | D2 | normal return with the context's reason; the interrupt condition, then a run equal to an uninterrupted one |
| sampler-wide `setOffset` after chain offsets; `updateScale = TRUE` while chains differ | D2 | every chain back on the shared offset; a refusal naming the cause |
| `capi_hash` expectations, including a probe that a changed view-op or hook-typedef parameter moves the token; the header compiles as C99 `-pedantic` | D2 | the token moved; clean |
| ThreadSanitizer: tests/cpp built with `-fsanitize=thread`, running the threaded-hook tests | D1, CI | no report (call 12) |
| mutants: the hook fired after the augmentation draw; the hook skipped on workers; the view writing chain 0; latents of the wrong chain in the draw | D1, D2 | each caught by a gate above |
| `defaultNumThreads` with `_R_CHECK_LIMIT_CORES_` unset, `"false"`, `"FALSE"`, `"TRUE"`, `"warn"`, `"WARN"`; an explicit `n.threads` and `mc.cores` under it | D3 | the machine's count, the same twice, then 2 three times; the explicit values untouched |
| stan4bart: chains = 1, seeded, continuous and binary, against a 89af617 install | S1 | bitwise equal (see above) |
| stan4bart: chains = 3, `cores` = 1 and 2 (4 locally) | S1 | bitwise equal `bart_train`, `stan`, `sigma` |
| stan4bart: chain 1 of chains = 4 against chains = 1, a fixed seed, at cores = 1 and 2 | S1 | bitwise equal |
| stan4bart: chains = 3 seeded; every pair of chains | S1 | differ in `stan`, `sigma` and `bart_train` |
| stan4bart: two unseeded fits in a row | S1 | differ |
| stan4bart: `stan4bart.chunk` = 1 against a whole phase, with and without `verbose = 2`, continuous and binary | S1 | bitwise equal |
| stan4bart: `keep_fits = FALSE` | S1 | the one stored draw equal to the last draw of a `keep_fits = TRUE` fit |
| stan4bart: a forced parametric failure on chain 2 | S1 | one R error naming chain 2; no warning; no rerun |
| stan4bart: test-24's interrupt through dbarts' count hooks, extended to chains = 2, cores = 2, with `stan4bart.chunk = 1`: on the worker path the hooks count main-thread polls (one at each run's entry, then one per 100 ms), so one-sweep runs make the injected poll land inside the fit however fast it is; the checks do not depend on the iteration it lands in | S1 | R's interrupt condition; a fresh fit equal to one taken before |
| stan4bart: 20 seeds, a continuous and a binary fit, against 89af617 | S1 | posterior means of fixed effects, sigma, random-effect sd and the BART fit mean agree under the threshold the equivalence harnesses' statistical mode uses; a fit time at n = 500 and 5000 recorded, no bar |
| stan4bart: a one-chain, `cores = 1` binary fit at n = 5000 against 89af617, alternating runs on one machine | S1 | no slowdown above noise; the draw callback computes every carried channel at each kept draw, the log-likelihood with a `pnorm` per row included, where today's loop asked for neither |
| stan4bart mutants: every chain given the first parametric seed; the draw callback's offset subtraction using the next group's offset; the hand-off moved after the next parametric step; the warmup run made with `numBurnIn = warmup` | S1 | caught by the pairwise, the one-chain and the chunk gates |
| stan4bart: R-loaded ThreadSanitizer (or ASan) run of test-01 and test-02 at cores = 2, macOS or the bench box | S1, local | no report |
| stan4bart: audit of the parametric code and its headers for mutable statics, `Rprintf` and R API on the hook's and the draw callback's paths | S1, review | none, or each removed |
| stan4bart: predict, extract (all types), fitted, saveRDS round trip, `getBartSampler` rebuild | S2 | equal to the stored blocks to 1e-10, as test-19 does today |
| stan4bart: an R callback at cores = 1 and 2; with and without `keep_fits` | S3 | results equal to each other and in today's shape |
| stan4bart: the test-only compiled callback at cores = 1 and 2; a wrong-major pointer; a throwing and a nonzero-returning test callback | S4 | equal results in the R callback's shape; a refusal at registration; one R error each |
| CRAN: every test and example asks at most 2 threads; test-01's `cores = 4` becomes 2 | S1 to S4, B1 | `R CMD check --as-cran` clean |
| timing thresholds: dbarts, stan4bart and bartCause checked with `_R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5`, `_R_CHECK_TEST_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5` and `_R_CHECK_LIMIT_CORES_=TRUE` | D3, S1, B1 | no timing NOTE |
| bartCause: its stan4bart tests | B1 | pass; re-recorded values only where a seed moved |

Gate hygiene ([Gate hygiene](README.md#gate-hygiene)): `--preclean` on every engine commit (the facade's
virtuals change), mutation reverts touched, warnings counted with `withCallingHandlers`.

CI. A ThreadSanitizer leg on tests/cpp is one more job in the cpp-tests workflow on Linux, where gcc and clang
both carry it, costing minutes. An R-loaded TSan leg would need R's own frames uninstrumented and is noisy,
so it stays a local gate at S1's landing (call 12).

## Sizes

Planned lines, then the forecast by the house's record. Recent dbarts slices ran 1.4 to 2.3 times their plans,
with half to four fifths of the diff in tests (the records in forest-sd-unit and forest-kind-by-class). As
there, tests are assumed at 2.5 times and the rest at 1.3.

| slice | package code | tests | docs, help, records |
|---|---|---|---|
| D1 engine | ~140 (hook type, ChainView, worker pass, latents in DrawInfo, atomic flag) | tests/cpp ~380 | architecture's Threading model and Callbacks and errors ~40 |
| D2 bridge and header | bridge ~200 (view, registration, draw fill, token fold), header ~170 (mostly the contract's text) | consumer.c, a C++ consumer and test-capi.R ~450 | bart-as-a-component section 5, NEWS ~80 |
| D3 thread default | R ~50 | ~40 | Rd ~60 |
| S1 stan4bart loop | C++ ~420, R ~-170 net | ~300 | help, NEWS ~60 |
| S2 restore and readers | R ~120 | ~100 | ~20 |
| S3 R callback | ~120 | ~50 | ~20 |
| S4 compiled callback | header ~80, C++ ~120, R ~60, test-only callback ~30 | ~150 | help and example ~60 |
| B1 bartCause | ~15 | ~30 | NEWS ~5 |

dbarts: package ~560, tests ~870, the rest ~180; ~1610 planned, forecast ~3140 (2900 to 3800). stan4bart:
package ~830 net, tests ~600, the rest ~160; ~1590 planned, forecast ~2790 (2500 to 3300). S4 alone is ~290
plus the header against the ~250 plus header the maintainer was told. Stop and report when a repository's diff
passes the top of its forecast, or its package code passes 1.5 times its figure (dbarts ~840, stan4bart
~1250). An overrun in tests is expected; one in package code means the design moved.

## Order

1. D1 on branch `threaded-chains` in dbarts; review; gates.
2. D2 on the same branch; review; gates; TSan job added.
3. D3 on the same branch; one sonnet review.
4. S1 on a stan4bart branch building against the dbarts branch; review with mutants; gates. A defect S1 finds
   in the hook or the view is fixed on the dbarts branch and re-gated there.
5. The dbarts branch merges to main once, green, with its records, after S1's gates pass against it. A release
   cut never freezes a hook no consumer has run.
6. S2, S3, then S4; one R-surface review each; stan4bart merges once.
7. B1 in bartCause after stan4bart merges.

## Per-chain values on a copy (ruled)

Ruled 2026-11-11 (the maintainer: "Ok, fair enough. Do (a) but write (b) down as a post release
enhancement."). The hook's offsets live in the holder's per-chain buffers. Per-chain sigma already rides the
state, since each chain's state records its sigma, so only the offset was in question.
- Built here, (a): a copy, a stored state and `sampler$data` carry the sampler-wide offset. The chains' own
  offsets are the larger sampler's state, set again before every sweep. D2 states this in the header beside
  the hook and in the sampler's help. stan4bart needs no more: after a fit it only predicts, and predict takes
  its offset as an argument.
- After the release, (b), an enhancement on the same buffers: per-chain offsets recorded in copies and states,
  `$setOffset(chain = )` on the R object, and every reader of the offset choosing a chain (~420 lines by the
  per-chain responses design note, [per-chain-responses.md](../design/per-chain-responses.md), section 2). It
  is additive, and nothing here forecloses it. Its TODO entry comes with the records batch.

## Open calls

For the maintainer; the design fits either answer.

None. The exception rule for compiled hooks, once conditional here, is settled by dec-B581 (call A6).

## Agent calls

Made here, for the close-out (dec-B629). Each names its alternative.

2. **Seed rule.** One rule at every `cores`: each chain's parametric and BART seeds derive from the user's seed
   and the chain number alone (pass-2 finding C135), so chain 1 of C chains is the one-chain fit. Rejected
   draft: C parametric seeds, then the sampler's C chain seeds, from one running stream, under which chain 1
   depends on C. Alternative: keep today's `cores > 1` draws, which would need a dbarts
   creation entry taking per-chain seeds (an R surface change) and would still move `cores = 1`. One-chain
   fits keep their draws; chains 2 and on of a multi-chain fit move once. `bart_args$seed` stops giving every chain one BART
   stream.
5. **Progress.** One main-thread line per `refresh` between chunks. Alternative: bart()'s per-chain lines
   through dbarts' queued sink, whose counts are per run and so read wrong under chunks, and which needs a
   report operation on the view. Trivial either way.
6. **Chain failure.** One R error after every chain stops. Alternative: today's warning and serial rerun,
   which hides a diverged chain and doubles the time.
7. **Entry shape.** Chain operations in a library-filled view, separate from the sampler-wide entries.
   Alternatives: a mode flag on the existing entries (hidden state), or the existing entries promised R-free
   (they raise on NaN).
8. **Public.** In dbarts.h. Alternative: private callables in lock-step with stan4bart. A private accessor
   needs a known saving across packages we maintain in lock-step, and treatSens is a second consumer.
   Cost: after the release the view is part of the contract.
12. **Guards.** A tests/cpp ThreadSanitizer job in CI, and a debug-build assertion that no `dbarts_sampler_*`
    entry runs on a sampler while its run is in flight. Alternative: local-only TSan, which a later engine
    change would not meet.
14. **The check cap** (the coordinator's call). The default honors `_R_CHECK_LIMIT_CORES_` through an exported
    `defaultNumThreads(n.tasks)` that every default calls, stan4bart's included. Alternatives: an argument on
    `guessNumCores` (blurs its meaning), or stan4bart copying the rule (two homes for one rule); or no cap,
    consistent with bart() today, leaving every test and example to pass a count. Scope beyond the three
    fitters named: `dbartsControl()` and `xbart()`, whose defaults read the same function.

A1. **One engine hook type** for the R closure and the compiled hook. Alternative: keep `SweepCallback`
    beside a new worker hook, which gives two call sites in `Chain::runSweeps`.
A2. **Record and hand off in the draw callback** (from the critique). The draw callback fires on the chain's
    thread right after each group's draw, with the channels stan4bart stores. Recording there subtracts the
    offset that draw was made under, keeps memory to what `keep_fits` asks, and keeps today's order with no
    run-boundary special case. Alternative, the first draft: copy the run's results after each run and do
    the last hand-off on the main thread. That cannot subtract per-draw offsets, which have been overwritten
    by then, without a store of n x iterations x chains; it allocates iterations x chains x (n + n.test) even
    under `keep_fits = FALSE`; and it would need a `getFits` and a test-fits operation on the view.
A4. **Function pointers in the view** rather than two named entries in the list. Named entries resolve
    lazily through the stubs, which call `R_GetCCallable` and may `Rf_error` on first use, here on a worker.
    Cost: the view's member types enter the token through a list macro (A13).
A5. **Stop protocol.** Nonzero return and the caller's own context carry the reason, as the draw callback does.
    Alternative: a message buffer in the view that dbarts raises. That would put stan4bart's error text into
    dbarts' raise and still need a status.
A6. **No exception escapes either hook**, by dec-B581's rule: a consumer catches its own in its wrapper, keeps
    it, returns a negative value and re-raises after the run returns. dbarts still catches one at the call on
    every thread as a safety net, tested across libraries in D2, and the header says the net is no portable
    contract across separately built C++ runtimes (Windows DLLs). stan4bart's slices follow the rule.
    Alternative: the net as the contract, which fails where unwinding does not cross libraries.
A7. **The host-neutral statuses** (dec-B119, dec-B571, dec-B572): 0 `DBARTS_OK`, -1 `DBARTS_UNSUPPORTED`, -2
    `DBARTS_INVALID_ARGUMENT`, defined now from `DBARTS_STATUS_LIST` in that plan's three-argument form (name,
    value, words, all folded into the token); the host-neutral work adds its other codes to the same list.
    Alternative: a separate three-valued chain status, which the host-neutral header would carry beside its
    own; or the capability 0/1 with a raise on bad values, which is not allowed on a worker.
A8. **Latents in `dbarts_draw`.** Alternative: keep a `getLatents` on the view and hand a binary fit's latents
    over at the next group's hook. That splits the hand-off across two calls and, at a run's last draw,
    pushes it past the adaptation freeze at the warmup boundary, out of today's order.
A9. **Flat runs only, not saved.** Alternative: the hook fires on R's `$run()` too. No R caller can register
    one, and a hook firing on a run its registrant did not make would surprise it.
A10. **Refuse a sampler-wide re-derivation while chains hold different offsets.** Alternative: read chain 0's,
    as `Sampler::rederiveScale` does today, which silently gives chains 2 to C a transform read from another
    chain's data.
A11. **Version.** The token re-bakes; the minor stays 0 while no release has shipped, by the header's rule.
    Moves to 1 if 1.0-0 ships first.
A12. **Pre-range restores.** Kept only if a shipped stan4bart wrote such fits. Fits from development builds
    are not owed compatibility (VD 2026-10-17). If none shipped, the per-chain `fit.scale` branch and its test
    go. Alternative: keep it regardless, which keeps a one-sampler-per-chain path alive beside the one
    sampler.
A13. **Fold member and typedef signatures into the token** (from the critique). Alternative: accept the
    same-width blind spot the header already admits for struct fields, which here would cover every
    parameter of every hook.
A14. **Predict threads from the session's rule**, so a fit saved on a large machine does not ask for its
    count on a small one. Alternative: the saved control's count, which is what passing 0 gives today.
A15. **A body-only helper for the compiled callback** (`compile_callback(body, length)`, calling
    `Rcpp::sourceCpp` with the signature, a `noexcept` wrapper that catches inside the user's own library, and
    the export written for the user). Cost ~50 lines and a Suggests on Rcpp, which S4 takes anyway for its
    example. Alternative: the header and a documented template only. Recommended: build it. A user's
    exception then never crosses a library boundary, which makes call A6's doubt moot for this path. Marked
    for close-out.
A16. **The callback names the channels it reads** (ruled 2026-11-10, "Use your recommendation for this item"):
    `dbarts_sampler_setDrawCallback` takes a channel mask, so stan4bart, which reads train and the latents,
    pays for no log-likelihood or test pass it never reads. The scalars are always filled, being free and
    having no NULL. Alternative, the plan's first course: time the log-likelihood scratch and add a second
    registration entry only for a cost above noise.
A17. **One log-gamma, the C library's reentrant form** (`misc_logGamma`), after the TSan job caught chains
    racing on `signgam` through `lgamma`. Bitwise `lgamma`'s value wherever R builds: `lgamma_r` on glibc,
    musl, macOS and the BSDs, declared in the C file since each library hides it behind its own feature
    macro, and libmingwex's `__lgamma_r` on Windows, the function its `lgamma` wraps, exported since 2007
    (its `lgamma` passing `&signgam` to it since 2012) but declared in no header. No configure check: every
    Unix C library R builds on has `lgamma_r`, and one without fails the install's load test by name.
    Alternatives: a lock around `lgamma` (bitwise, but a lock on the per-draw negative binomial path and the
    variance forest's moves), or an own implementation (lock-free, but last-bit changes, so every draw moves
    on every build).
A18. **The Student-t log-likelihood is the engine's own density** (`StudentTLogDensity`), since R's `dt`
    reaches the C library's `lgamma` through its `stirlerr` at a non-integer df up to about 10.5, so two
    chains at a fixed df of 3.5 raced on `signgam` inside libR, where the symbol read cannot look. The plain
    form, a constant per df through `misc_logGamma` and one `log1p` per row, with the Stirling series'
    difference for the constant from df 15, where the two log-gammas cancel, agrees with `dt` to a relative
    1e-14 over df in [1e-6, 1e10]; the channel's last bits move, no draw does. Alternatives: a port of R's
    `dt`, `stirlerr` and `bd0` with `stirlerr`'s `lgamma` swapped (about 200 lines and `lgamma1p`'s
    dependencies, to keep `dt`'s last bits), or a lock around `dt` (on the per-draw channel). Of R's routines a chain then calls only `pnorm`, `dnorm` and
    `qnorm`, pure arithmetic; `qt`, `qchisq` and `lchoose`, which can reach `lgamma` (the last never at
    whole counts), run on R's thread, at construction and in setters.

## Risks and the objections answered

- **"Foreign code mutating a chain inside a run reverses the header's between-runs rule."** It narrows rather
  than reverses it. The hook gets a view of one chain with two writes, both with the transform held, at a
  sweep boundary the engine already treats as between sweeps. The augmentation draw follows the hook, so the
  change reaches the sweep as it would after a run. The per-sweep identity of section 1 of
  bart-as-a-component is the gate, chunked against one-sweep.
- **"A no-throw, no-R contract on code dbarts calls."** That contract already exists for the draw callback on
  workers. The view's statuses mean the hook never needs R to refuse. D2 tests the exception path across
  libraries, and TSan and the debug assertion catch the rest in testing.
- **"Per-chain offsets a copy would silently drop."** Not silently: under the ruled (a) the header and the
  help say what a copy carries, and (b) is the post-release enhancement.
- **"Route 3 costs route 1's work plus a hook."** It does not. No R surface, no per-chain response, no
  per-chain transform: offset buffers, one hook, two view operations and one appended draw field (synthesis
  fact 3).
- **Recording the wrong quantity** (the critique's blocking item). The draw callback subtracts the offset
  that draw was made under; the one-chain bitwise gate and the offset mutant pin it.
- **Identical chains passing every equality gate** (the critique's mutant). The pairwise-differ gate catches
  it.
- **A static in the parametric code shared across chains** (WALNUTS, Eigen, stan4bart's `activeSamplers`).
  The audit item and the R-loaded sanitizer run catch it. `activeSamplers` is touched only at creation and
  finalization, on the main thread.
- **Oversubscription under `mclapply`.** A fit inside forked children threads by default, as bart() does. That
  follows from the ruling and is not changed here.
- **Gain unmeasured.** About 3.5x on 4 cores is the estimate from measured parts (synthesis fact 5). The fit
  time is recorded at S1 with no bar.

## Later items

- **Route (a) off-main-thread entries**, built only when a consumer needs them: a run taking a caller's stop
  poll and returning a status, and a main-thread call for the deferred
  slow-count warning. The first candidate is treatSens with `trt.model = bart()`, whose two samplers advance in
  lockstep and so cannot share one sampler's chains. A TODO entry, not built here.
- **treatSens' default treatment models** (`probitEM()`, `probit()`) as a later consumer of this hook. By the
  ruling's account, each sensitivity grid cell becomes a chain carrying its own offset, which is the outcome
  model's term in the unobserved confounder. The treatment-model step runs in the hook and the draw callback.
  It would need what this plan builds: a per-chain offset, sigma for a continuous outcome, the draw's fits,
  and latents for a binary outcome. It would also need seeds per cell and a grid that fits in one sampler's
  chain count. Not verified against treatSens' loop here; its plan reads that loop first.
- **Per-chain responses** (TODO per-chain-responses) stay as they are: additive later, for R-level loops.

## Verification

- dbarts:
  `R CMD INSTALL --preclean -l <lib> .`, then
  `cd tests/cpp && make && ./test_bartcore`, then the TSan build of tests/cpp
  (`OPT="-O1 -g -fsanitize=thread"`) and its run. After that,
  `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts", lib.loc = "<lib>")'` in one process,
  `lintr::lint_package()`, `air format --check .`, and `tools/check-rc-codoc.R`, `tools/check-win-drift.R`,
  `tools/check-doc-freshness.R` and `benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit
  status. The reference-build snapshots must be unchanged, the ASan R-loaded run clean, and bench-sampler.R
  compare (maintainer) within noise. `R CMD check --as-cran` runs with the three timing and core variables
  named in the gates.
- stan4bart: its tinytest suite against the dbarts branch library, then the one-chain bitwise script against a
  89af617 install in a separate library, then the thread-count, pairwise, unseeded and chunk gates, then the
  0.0-13.1 run, then `R CMD check --as-cran` with the same variables.
- bartCause: its suite against both branch libraries, then its check with the same variables.
