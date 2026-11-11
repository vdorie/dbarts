# threaded-chains: stan4bart's chains on dbarts' threads, one sampler for every chain

Status: PLANNED 2026-11-11, before its blind critique. Route (b) ruled 2026-11-11; the cores default ruled
2026-11-11; the user's R callback and per-chain values on a copy still before the maintainer
([Open calls](#open-calls)).

agent: one opus implementer per repository line (dbarts D1 and D2, stan4bart S1 to S3), effort high; one opus
reviewer who runs mutants on D1 and D2 and on S1; a sonnet reviewer for S2, S3 and the bartCause slice B1.
rng: NEUTRAL for dbarts: with no hook registered every run takes the path it takes today and no baseline,
snapshot or exact gate moves. SHIFTING for stan4bart fits of more than one chain (chain seeds are drawn in a new
order, the same at every thread count); a one-chain fit is expected to keep its draws bit for bit
([Draws at a fixed seed](#draws-at-a-fixed-seed)). Not posterior-changing anywhere: the Gibbs order inside a
chain is today's.
window: before the 1.0-0 release (scheduled by the maintainer 2026-11-11). dbarts lands first; stan4bart builds
against the dbarts branch until it merges.
budget: dbarts ~1690 lines planned (forecast 2900 to 3800); stan4bart ~1210 (forecast 1900 to 2500); bartCause
~50. See [Sizes](#sizes).

## Goal

A stan4bart fit holds one dbarts sampler with every chain in it. Each chain carries its own offset and sigma,
and stan4bart makes one dbarts run per chunk of sweeps. Before each sweep, the dbarts worker thread that owns
chain c calls a compiled stan4bart hook. The hook draws that chain's parametric part and sets that chain's
offset and sigma through entries that reach that chain alone, never call R and never raise. dbarts keeps all
thread machinery. Chains run on min(cores, chains) threads; the cores default is the ruled one. Seeded draws are
bitwise identical at every thread count. The PSOCK and FORK code goes.

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
default changed from 0.0-13, which was `getOption("mc.cores", 1L)` on processes (read in stan4bart's 0.0-13.1
tag, not run).

Earlier rulings this rests on: dec-B419 (the maintainer's question whether one sampler could hold every
chain); dec-B437 (stan4bart's response range is fixed at creation and shared, so chains have no transform to
agree on and need not meet inside a run); dec-A197 (today's restore, one sampler per chain, which S2
replaces).

Code read for this plan: dbarts main 6a34f43b1; stan4bart main 89af617; bartCause main ec1ee3e; treatSens main
268b2db.

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
- The per-draw observer [`dbarts_draw_callback`](../../inst/include/dbarts/dbarts.h) already runs on workers
  with no lock, under a no-R contract. It is the precedent for compiled consumer code running on dbarts'
  workers.
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
  [`userInterrupted`](../../src/R_interface_bartcore.cpp); the setters raise on a non-finite value. Both stay
  main-thread entries; nothing here changes them.

stan4bart today, read at 89af617: `stan4bart_fit_worker` builds one `dbartsSampler` per chain
(`n.chains = 1`, `n.threads = 1`, `n.thin = skip.bart`) and calls `C_stan4bart_run` once per phase. Per
iteration that call runs the parametric step (the ridge move in sampling, then `skip` WALNUTS transitions),
reads the parametric mean into the BART offset, sets sigma and the offset, runs one dbarts run of one kept draw,
subtracts the offset from the training fits to give the parametric part its offset, hands a binary fit's
latents to the parametric part, and evaluates the user's R callback. The parametric code draws from its own
`std::mt19937_64` and reports failure by C++ exceptions; it calls R only through stan4bart's data readers at
creation and a debugging print in `double_writer.hpp`. `cores > 1` starts PSOCK (else FORK) workers, seeds each
chain from `sample.int` after `set.seed(seed)`, and falls back to serial with a warning when the cluster fails
to start or a chain fails. `cores = 1` draws every chain's seeds from one running stream, so seeded draws
depend on `cores` (run 2026-11-11, synthesis fact 6).

## Design

### The hook and the chain view (dbarts.h)

```c
typedef struct dbarts_chain_view_t dbarts_chain_view;
struct dbarts_chain_view_t {
  size_t structSize;       /* the library's sizeof, as dbarts_draw */
  size_t chainIndex;
  size_t numObservations;
  void* internal;          /* the library's; never read by the caller */
  int (*setOffset)(dbarts_chain_view* view, const double* offset);
  int (*setSigma)(dbarts_chain_view* view, double sigma);
  int (*getLatents)(const dbarts_chain_view* view, double* out);
  int (*getFits)(const dbarts_chain_view* view, double* out);
};
typedef int (*dbarts_sweep_hook)(void* context, dbarts_chain_view* chain,
                                 size_t sweepIndex, int isBurnIn);
void dbarts_sampler_setSweepHook(dbarts_sampler* sampler,
                                 dbarts_sweep_hook fn, void* context);
```

The sketch fixes the shape; the implementer settles names in review. One entry joins `DBARTS_C_API_LIST`
(registration). The four chain operations are function pointers in a view that the library fills, valid only
during the call it is passed to (call A4).

- Registration copies (fn, context) into the holder on the copy-on-set rule. It allocates one n-double
  offset buffer per chain, each seeded from the sampler-wide offset in force (zeros if none). A null fn clears
  the hook and returns every chain to the shared offset. Registration runs on the main thread and may allocate;
  no chain operation does. The hook fires on runs through `dbarts_sampler_run` only, as the draw callback does,
  and is not part of a saved state or a copy (call A9).
- The hook fires before every sweep of every chain, on the thread that owns the chain: a worker, or the calling
  thread when the run is inline. It gets the 0-based sweep index within this run and the burn-in flag.
  0 continues. Nonzero stops the run: every chain stops at its next sweep boundary, the run's results are
  discarded and its saved trees dropped, exactly as a draw callback's stop. `dbarts_sampler_run` then returns
  normally, and the caller reads why from its own context (call A5).
- The contract the header states: no R call of any kind (allocation, printing, errors, RNG, interrupt checks)
  on any path, inline included, so one rule holds at every thread count. No C++ exception may leave the hook;
  dbarts catches one defensively at the call and treats it as a stop that raises an R error after the join,
  but unwinding across two shared libraries is not promised on every toolchain (call A6). The hook may
  touch only its own chain, through the view. It must not call any `dbarts_sampler_*` entry on this sampler;
  a debug build asserts that (call 12).
- `setOffset(view, offset)`: n finite values copied into that chain's buffer; the transform is held, never
  re-derived (no updateScale argument). Returns 1 on a write; 0, touching nothing, where the sampler has no
  offset conduit (the capability `dbarts_sampler_setOffset` reports as 0); -1, touching nothing, on a
  non-finite value. `setSigma(view, sigma)`: 1 on a write; 0 where sigma is pinned, as
  `dbarts_sampler_setSigma` decides; -1 on a value that is not finite and positive. One status enum,
  `dbarts_chain_status`, generated from a list macro like the other four (call A7).
- `getLatents(view, out)`: that chain's n augmentation values, the ones the previous sweep's trees were drawn
  from, with `dbarts_sampler_getLatents`'s return rule. `getFits(view, out)`: that chain's current training fits
  with exactly the values and layout the `train` channel would record for the last sweep (n x F, the chain's
  offset included), so it adds no new quantity (call A8).
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
- `ChainView` offers four operations: set the offset with the transform held (a `Chain` method that calls
  `response_->setOffset(ptr, false, ...)` and names no saved slot, so it never touches the shared tree store);
  set sigma; read the latents; copy the training fits as `storeSample` would write them (the same code, not a
  copy of it).
- The hook's stop takes the path the draw relay's stop takes: a set-once relaxed flag read after the join, so
  `Sampler::run` can tell a hook's stop from an interrupt. A hook's exception is caught in the worker body that
  already catches chain exceptions, so nothing new is needed there.
- `running_` becomes atomic so the debug assertion (call 12) can read it from a worker. No draw-path code
  changes; the no-hook path stays one null test per sweep, as today.

### Bridge (D2)

- The holder gains `std::vector<std::vector<double>> ownedChainOffsets`: nested, as `ownedForestWeights` is,
  so a pointer a chain borrows survives later growth. It is sized at registration, empty otherwise. It also
  gains the shipped hook (fn, context), beside `drawHook`.
- The view's operations live in the flat entry file, with no R API inside (no `R_finite`; `std::isfinite`).
  The bridge fills one view per chain per run on the worker's stack.
- A sampler-wide `dbarts_sampler_setOffset` (or the R `$setOffset`) after chains hold their own offsets puts
  every chain back on the shared buffer, as it does today. A sampler-wide re-derivation of the transform
  (`updateScale` true, or [`Sampler::rederiveScale`](../../src/bartcore/sampler.hpp), which reads chain 0's
  offset) while chains hold different offsets is refused with an error naming the cause. The rule it would
  need, one transform read from which chain's offset, is not one a user means (call A10).
- `bartcore_runWithCallback` moves onto the single hook type unchanged in behavior.
- The token re-bakes ([`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)). It folds the new list entry,
  the hook typedef's parameters, the view's layout and the status enum. The version constants
  ([`DBARTS_C_API_MINOR`](../../inst/include/dbarts/dbarts.h)) do not move before the first release, by the
  header's own rule. If this lands after 1.0-0 ships, the minor moves to 1 and stan4bart's floor with it (call
  A11).

### stan4bart (S1 to S3)

- S1, the loop. `stan4bart_fit` builds one `dbartsSampler` with `n.chains = chains` and `n.threads =
  min(cores, chains)` and installs the shared range once. `C_stan4bart_create` reads the handle and holds C
  parametric samplers, each with its own buffers (BART offset, parametric offset, latents, parameter writer)
  and seed. It runs the one creation run (every chain once, no hook), does each chain's first hand-off, then
  registers the hook. The hook is a `noexcept` C++ function over a context that holds the C per-chain records
  and an error slot per chain. It acts only when `sweepIndex % skip.bart == 0`:
  1. if this is not the run's first group, it hands the previous group's draw to the parametric part (fits
     less this chain's offset through `getFits`; a binary fit's latents through `getLatents`), as today's loop
     does after its BART run;
  2. it runs the parametric step (ridge move in sampling, then `skip` transitions, one in warmup) and writes the
     chain's parameter row;
  3. it sets the chain's offset (the parametric mean, plus a user offset by `offset_type`) and, for a
     continuous fit, sigma.

  A WALNUTS exception, or a -1 from a view operation (a diverged parametric step's NaN), is caught, its
  message kept in the chain's slot, and the hook returns nonzero. After each run the main thread does the last
  group's hand-off for every chain, so each run boundary leaves every chain exactly where today's loop leaves
  it after an iteration. That includes the warmup-to-sampling boundary, where
  `C_stan4bart_disengageAdaptation` runs. It then copies the chunk's BART channels into the per-chain result
  arrays (the flat results have no chain stride beyond the chunk). The adaptation, eval-count and
  initial-summary entries loop over the chains. Everything the hook touches is owned by the external-pointer
  object, never by a `.Call` frame, so an interrupt's jump strands nothing.
- Chunks: a phase runs in one dbarts run unless a reason to stop sooner applies. Verbose above 1 prints
  `iter i / N` every `refresh` iterations, and chunks end there. A user callback ends them too (open call 4).
  Chunk length cannot move draws; an internal option (`stan4bart.chunk`, unexported, test-only) forces a length
  so the gate can compare 1 against a whole phase.
- Progress: printed by the main thread between chunks. Every chain ends a chunk at the same iteration, so one
  `iter i / N` line describes them all; at one chain the lines are today's (call 5). Verbose 1's "fitting
  chain" lines become one line naming the chain count and thread count. The "starting multithreaded fit,
  further output silenced" line goes.
- Errors: a chain's failure stops every chain at its next sweep boundary and raises one R error after the join,
  naming the chain and its message. It replaces today's warning and silent serial rerun (call 6). An interrupt
  is dbarts': polled every 100 ms, it stops the workers at their next sweep boundary and raises R's interrupt
  condition after the join. Latency is at most one sweep plus one parametric step.
- Seeds (call 2): after `set.seed(seed)` (or R's stream as it stands when `seed` is NA), stan4bart draws the C
  parametric seeds with one `sample.int(.Machine$integer.max, chains)` and then builds the sampler, which
  draws the C BART chain seeds from R's stream in chain order. `.Random.seed` is put back afterwards, as now.
  Every generator is a chain's own and the hook touches only its chain, so draws do not depend on `cores`.
- The cluster code goes: `makeCluster`, `clusterMap`, `clusterExport`, the FORK fallback, both fallback
  warnings, the per-worker seeding and `parallel` from Imports if nothing else uses it.
- The cores default is as ruled. `cores` keeps its name and is validated as one positive whole number;
  `min(cores, chains)` threads run.
- S2, restore, predict, extract and save. `sampler.bart` becomes one `dbartsSampler`, and `state.bart$state`
  the one sampler's state. `restoreBartSampler` builds one sampler of C chains, installs the one range once
  and sets the state. `bart_pointer_is_live`, `getBartSampler` and `predict_bart_chains` (one
  `C_stan4bart_predictBART` call returning draws x chains) follow. `extract(type = "trees")` reads one
  sampler's chains by `chainNums`. A saved fit holds one sampler's data copy instead of C.
  `restoreBartSampler` keeps a branch that installs each chain's own `fit.scale` for fits saved before the
  one range. S2 first establishes which shipped version, if any, wrote such fits. A 0.0-13.1 fit, written
  through `C_stan4bart_exportBARTState` against dbarts 0.9-34, is loaded on main and on the branch (call A12).
- S3, the user callback, by the maintainer's answer to open call 4.
- S1's help and NEWS: `cores` describes threads, the new default, the 2-thread cap CRAN asks of examples, and
  the draw change. NEWS diffs against 0.0-13.

### bartCause (B1)

`R/argParse.R` gives stan4bart `cores = min(n.threads, chains)` only when the fit is long enough
(`chains * iter >= 4000`), because a cluster costs a fixed startup. When it does, it also injects a seed so a
clustered fit stays reproducible under `set.seed`. Under threads neither reason holds. B1 passes `cores =
n.threads` when the user gave `n.threads`, and otherwise leaves `cores` unset so stan4bart's default applies. The
injected seed goes, since stan4bart now seeds from R's stream at every core count. bartCause tests that fit
stan4bart with a seed are re-recorded, and its tests and examples pass at most 2.

## Draws at a fixed seed

- dbarts: none move. No hook, no new path.
- stan4bart, one chain: today's serial order is "parametric seed by `sample.int(max, 1)`, then the sampler's
  one chain seed", and the new order at C = 1 is the same. The hook keeps the per-chain Gibbs order, and the
  run-boundary hand-off keeps every run boundary where today's leaves it. So a one-chain fit is expected
  to reproduce stan4bart main 89af617 bit for bit, and S1 gates on it. If something else moves it, the
  implementer names the cause and the reviewer judges. The equality would not come back by accident.
- stan4bart, more than one chain: every seeded draw moves once, at `cores = 1` and at `cores > 1` alike. From
  then on, `cores` = 1, 2 and 4 give bitwise-identical fits.
- A side effect, by the code: today `bart_args = list(seed = s)` gives every chain's own sampler the same
  control seed, so every chain draws the same BART stream. In one sampler the seed generator gives each chain
  its own seed, which is what `seed` means in dbarts (a one-chain run with seed s reproduces chain 1 of any
  run with seed s).
- Snapshots: stan4bart has no reference-build snapshot files. Its seeded expectations that hold values
  (test-05-rng holds identities, not values) are re-recorded by replaying whole files. bartCause's
  stan4bart-backed seeded tests are re-recorded in B1.

## Gates

Tiers by [Process by risk](README.md#process-by-risk). D1 and D2 change the engine: plan with blind critique,
an opus review that runs mutants, the RNG class's battery for NEUTRAL, the full tinytest suite, sanitizers on
tests/cpp and on the R-loaded path (the chain view is new code reachable from a `.Call`), `R CMD check`, and
the arm64 reproducibility snapshots, which must not move. S1 changes stan4bart's draws: the same tier on
stan4bart's side. S2, S3 and B1 are R-surface tier.

| gate | slice | expected |
|---|---|---|
| tests/cpp: a hook that sets each chain's offset and sigma from its own fits (a toy Gibbs step); one run of m sweeps against m one-sweep runs with the same changes made by `Chain::setOffset` and `Chain::setSigma` between them | D1 | bitwise equal, every channel, at 1, 2 and 4 threads |
| tests/cpp: the same sampler at 1, 2, 3 and 4 threads (3 leaves one worker with two chains) | D1 | bitwise equal |
| tests/cpp: a hook stop on chain 1 of 4 at sweep k; a hook that throws | D1 | the run reports the stop (not an interrupt); the store's cursor and counts as after a draw callback's stop; the throw comes back after the join |
| tests/cpp: chain c's offset write leaves every other chain's working response and fits bit for bit | D1 | equal |
| test-capi.R plus capi/consumer.c: a compiled hook through the flat API doing the toy Gibbs step; chunked against one-sweep flat runs with `dbarts_sampler_setOffset` or setSigma between, at one chain (where the two routes are the same model) | D2 | bitwise equal |
| the same, 4 chains at 1 and 2 threads | D2 | bitwise equal |
| consumer: NaN offset, sigma 0, a pinned-sigma sampler (probit), a no-offset-conduit sampler | D2 | -1, -1, 0, 0 and the chain unchanged; nothing raised |
| consumer: a hook returning nonzero; an interrupt, in tests/cpp through the poll function (a CRAN test cannot send one) | D2 | normal return with the context's reason; the sampler runs on |
| sampler-wide `setOffset` after chain offsets; `updateScale = TRUE` while chains differ | D2 | every chain back on the shared offset; a refusal naming the cause |
| `capi_hash` expectations; the header compiles as C99 `-pedantic` | D2 | the token moved; clean |
| ThreadSanitizer: tests/cpp built with `-fsanitize=thread`, running the threaded-hook tests | D1, CI | no report (call 12) |
| mutants: the hook fired after the augmentation draw; the hook skipped on workers; chain view writing chain 0; the run-boundary flag not reset | D1, D2 | each caught by a gate above |
| stan4bart: chains = 1, seeded, continuous and binary, against a 89af617 install | S1 | bitwise equal (see above) |
| stan4bart: chains = 3, `cores` = 1 and 2 (4 locally) | S1 | bitwise equal `bart_train`, `stan`, `sigma` |
| stan4bart: `stan4bart.chunk` = 1 against a whole phase, with and without `verbose = 2` | S1 | bitwise equal |
| stan4bart: a forced parametric failure on chain 2 | S1 | one R error naming chain 2; no warning; no rerun |
| stan4bart: 20 seeds, a continuous and a binary fit, against 89af617 | S1 | posterior means of fixed effects, sigma, random-effect sd and the BART fit mean agree under the threshold the equivalence harnesses' statistical mode uses; a fit time at n = 500 and 5000 recorded, no bar |
| stan4bart: R-loaded ThreadSanitizer (or ASan) run of test-01 and test-02 at cores = 2, macOS or the bench box | S1, local | no report |
| stan4bart: audit of the parametric code and its headers for mutable statics, `Rprintf` and R API on the hook's path | S1, review | none, or each removed |
| stan4bart: predict, extract (all types), fitted, saveRDS round trip, `getBartSampler` rebuild | S2 | equal to the stored blocks to 1e-10, as test-19 does today |
| CRAN: every test and example asks at most 2 threads; test-01's `cores = 4` becomes 2 | S1 | `R CMD check --as-cran` clean |
| bartCause: its stan4bart tests | B1 | pass; re-recorded values only where a seed moved |

Gate hygiene ([Gate hygiene](README.md#gate-hygiene)): `--preclean` on every engine commit (the facade's
virtuals change), mutation reverts touched, warnings counted with `withCallingHandlers`.

CI. A ThreadSanitizer leg on tests/cpp is one more job in the cpp-tests workflow on Linux, where gcc and clang
both carry it, costing minutes. An R-loaded TSan leg would need R's own frames uninstrumented and is noisy,
so it stays a local gate at S1's landing (call 12).

## Sizes

Planned lines, then the forecast by the house's record. Recent dbarts slices ran 1.4 to 2.3 times their plans,
with half to four fifths of the diff in tests (the records in forest-sd-unit and forest-kind-by-class). As there, tests
are assumed at 2.5 times and the rest at 1.3.

| slice | package code | tests | docs, help, records |
|---|---|---|---|
| D1 engine | ~170 (hook type, ChainView, worker pass, fits copy, atomic flag) | tests/cpp ~380 | architecture's Threading model and Callbacks and errors ~40 |
| D2 bridge and header | bridge ~260, header ~170 (mostly the contract's text) | consumer.c and test-capi.R ~420 | bart-as-a-component section 5, NEWS ~80 |
| S1 stan4bart loop | C++ ~450, R ~-150 net | ~300 | help, NEWS ~60 |
| S2 restore and readers | R ~120 | ~100 | ~20 |
| S3 callback | ~120 under option (a), ~15 under (b) | ~50 | ~20 |
| B1 bartCause | ~15 | ~30 | NEWS ~5 |

dbarts: package ~600, tests ~800, the rest ~290; ~1690 planned, forecast ~3160 (2900 to 3800). stan4bart:
package ~690, tests ~450, the rest ~100; ~1210 planned, forecast ~2160 (1900 to 2500). Stop and report when a
repository's diff passes the top of its forecast, or its package code passes 1.5 times its figure (dbarts
~900, stan4bart ~1040). An overrun in tests is expected; one in package code means the design moved.

## Order

1. D1 on branch `threaded-chains` in dbarts; review; gates.
2. D2 on the same branch; review; gates; TSan job added; records. The branch merges to main once, green (dbarts
   additive: nothing that builds against main breaks).
3. S1 on a stan4bart branch building against the dbarts branch, then main; review with mutants; gates.
4. S2, then S3 once open call 4 is answered; one R-surface review each; stan4bart merges once.
5. B1 in bartCause after stan4bart merges.

Open call 9 does not gate any slice. Under option (b) an R surface slice follows in dbarts later.

## Open calls

For the maintainer; each design above fits either answer.

**4. The user's R callback.** `callback` runs R once per iteration. It cannot run on a worker.
- (a) The main thread evaluates it between chunks over the chunk's buffered draws, chain by chain and
  iteration by iteration. Chunks end every L iterations (L = 100, or `refresh` when verbose) while a callback
  is given. Within a chain the order is today's; across chains, calls interleave by chunk where serial runs
  today ran chain 1 to the end first. The cost is ~120 lines in stan4bart and a buffer of L x chains draws
  (n + n.test + parameters each) even when `keep_fits = FALSE`. Callback users get the speed.
- (b) A fit given a callback runs on one thread, with chunks of one iteration and the callback evaluated after
  each. A message says so when `cores` was set above 1 explicitly; under the new default it is silent, or it
  would speak on every callback fit. The cost is ~15 lines, and callback users get no speed.
- Recommendation: (a). The callback's results have the same shape and per-chain order. Only the interleaving
  across chains, which no reader depends on, differs.
- Either way the loop has chunks and a main-thread step after each run; (b) only sets threads to 1 and the
  chunk to 1.

**9. Per-chain values on a copy or reload.** The hook's offsets live in the holder's per-chain buffers.
Per-chain sigma already rides the state (each chain's state records its sigma), so only the offset is open.
- (a) Not recorded: a copy, a stored state and `sampler$data` carry the sampler-wide offset. The chains'
  own offsets are the larger sampler's state, set again before every sweep. Documented in the header and the
  sampler's help. This costs nothing, and stan4bart never needs more: after a fit it only predicts, and
  predict takes its offset as an argument.
- (b) Recorded: an R `chain` argument on `$setOffset`, per-chain vectors in copies and states, and every reader
  of the offset choosing a chain (~420 lines by the per-chain responses design note,
  [per-chain-responses.md](../design/per-chain-responses.md), section 2). It is additive later on the same
  buffers.
- Recommendation: (a) now, (b) when an R-level consumer asks.

## Agent calls

Made here, for the close-out (dec-B629). Each names its alternative.

2. **Seed rule.** One rule at every `cores`: parametric seeds, then the sampler's chain seeds, from R's
   stream after `set.seed(seed)`. Alternative: keep today's `cores > 1` draws, which would need a dbarts
   creation entry taking per-chain seeds (an R surface change) and would still move `cores = 1`. One-chain
   fits keep their draws; multi-chain fits move once. `bart_args$seed` stops giving every chain one BART
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

A1. **One engine hook type** for the R closure and the compiled hook. Alternative: keep `SweepCallback`
    beside a new worker hook, which gives two call sites in `Chain::runSweeps`.
A4. **Function pointers in the view** rather than four named entries in the list. Named entries resolve
    lazily through the stubs, which call `R_GetCCallable` and may `Rf_error` on first use, on a worker
    here. They would also need the caller to resolve them on the main thread first. Cost: a struct in the
    token's fold.
A5. **Stop protocol.** Nonzero return and the caller's own context carry the reason, as the draw callback does.
    Alternative: a message buffer in the view that dbarts raises. That would put stan4bart's error text into
    dbarts' raise and still need a status.
A6. **No exception across the ABI.** The hook is `noexcept` by contract, and dbarts catches one anyway.
    Alternative: promise that a thrown exception is carried back, as for the draw callback. Unwinding from
    one shared library through another is not known to hold under every toolchain CRAN builds with (Windows
    in particular), and nothing here has run it.
A7. **A three-valued status** (1, 0, -1). Alternative: reuse the capability 0/1 and raise on bad values, which
    is not allowed on a worker.
A8. **`getFits` returns the `train` channel's values**, offset included. Alternative: the forest's part
    alone, which saves the consumer one subtraction but adds a quantity the header does not define today.
    Alternative two: the hook reads the previous draw from the run's results buffer, which breaks at chunk
    boundaries and under null results.
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

## Risks and the objections answered

- **"Foreign code mutating a chain inside a run reverses the header's between-runs rule."** It narrows rather
  than reverses it. The hook gets a view of one chain with two writes, both with the transform held, at a
  sweep boundary the engine already treats as between sweeps. The augmentation draw follows the hook, so the
  change reaches the sweep as it would after a run. The per-sweep identity of section 1 of
  bart-as-a-component is the gate, chunked against one-sweep.
- **"A no-throw, no-R contract on code dbarts calls."** That contract already exists for the draw callback on
  workers. The view's statuses mean the hook never needs R to refuse. The defensive catch turns a broken
  contract into an R error, not a crash, wherever unwinding works, and TSan and the debug assertion catch the
  rest in testing.
- **"Per-chain offsets a copy would silently drop."** Open call 9 says it plainly either way, and under (a)
  the header and help say what a copy carries.
- **"Route 3 costs route 1's work plus a hook."** It does not. No R surface, no per-chain response, no
  per-chain transform: offset buffers, one hook, four view operations (synthesis fact 3).
- **Chunk equivalence fails at the run boundary.** The run-boundary hand-off is the one place the order could
  differ from today's. The one-chain bitwise gate against 89af617 and the chunk = 1 against whole-phase gate
  both cover it.
- **A static in the parametric code shared across chains** (WALNUTS, Eigen, stan4bart's `activeSamplers`).
  The audit item and the R-loaded sanitizer run catch it. `activeSamplers` is touched only at creation and
  finalization, on the main thread.
- **Gain unmeasured.** About 3.5x on 4 cores is the estimate from measured parts (synthesis fact 5). The fit
  time is recorded at S1 with no bar.

## Later items

- **Route (a) off-main-thread entries**, built only when a consumer needs them: a run taking a caller's stop
  poll and returning a status, status forms of setOffset and setSigma, and a main-thread call for the deferred
  slow-count warning. The first candidate is treatSens with `trt.model = bart()`, whose two samplers advance in
  lockstep and so cannot share one sampler's chains. A TODO entry, not built here.
- **treatSens' default treatment models** (`probitEM()`, `probit()`) as a later consumer of this hook. By the
  ruling's account, each sensitivity grid cell becomes a chain carrying its own offset, which is the outcome
  model's term in the unobserved confounder. The treatment-model step runs in the hook. It would need what this
  plan builds: a per-chain offset, sigma for a continuous outcome, the fits, and latents for a binary outcome.
  It would also need seeds per cell and a grid that fits in one sampler's chain count. Not verified against
  treatSens' loop here; its plan reads that loop first.
- **Per-chain responses** (TODO per-chain-responses) stay as they are: additive later, for R-level loops.

## Verification

- dbarts:
  `R CMD INSTALL --preclean -l <lib> .`, then
  `cd tests/cpp && make && ./test_bartcore`, then the TSan build of tests/cpp
  (`OPT="-O1 -g -fsanitize=thread"`) and its run. After that,
  `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts", lib.loc = "<lib>")'` in one process,
  `lintr::lint_package()`, `air format --check .`, and `tools/check-rc-codoc.R`, `tools/check-win-drift.R`,
  `tools/check-doc-freshness.R` and `benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit
  status. The reference-build snapshots must be unchanged, and the ASan R-loaded run must be clean.
- stan4bart: its tinytest suite against the dbarts branch library, then the one-chain bitwise script against a
  89af617 install in a separate library, then the thread-count and chunk gates, then
  `R CMD check --as-cran`.
- bartCause: its suite against both branch libraries.
