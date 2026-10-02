# Review 3 - lens: bridge

Pinned tree: .claude/worktrees/review3 at 01dee4b4. Library: scratchpad r3-lib. Probes:
scratchpad r3-bridge-p*.R (R level), r3-bridge-c*.R plus r3-bridge-shim/shim.cpp (a LinkingTo-style
consumer built with DBARTS_USE_STUBS and DBARTS_REQUIRE_EXACT_ABI against inst/include). The
0.9-34 comparisons ran against the CRAN 0.9-34 build in the scratchpad.

Covered: src/R_interface.cpp, src/R_interface_bartcore{,_common}.hpp/.cpp (the entry points, the
exception and unwind wrappers, the holder and finalizer, input validation), src/C_interface.cpp
and inst/include/dbarts/dbarts.h (every entry, the stubs, the token and layout asserts),
R/bartcore.R and the dbartsSampler methods that call them; stan4bart bartcore's use of the flat
API (src/init.cpp, R/generics.R). Probes: malformed arguments through every setter, gctorture over
about 40 sampler methods (gaussian, BCF, multinomial, ordinal, heteroscedastic), flat callbacks
that raise or throw (inline, worker, and the R route), 60,000 raising runs as a protect-stack leak
test, destroy then R methods, save/load/copy round trips, setTimeLimit interrupts at 1 and 2
threads, out-of-bounds canaries on flat run and predict per sampler kind, and malformed CSC sources.

Not covered: a reading PROTECT audit of all 8.5k bridge lines (gctorture and the existing rchk
pass stand in); fuzzing of malformed bartcoreState objects through setState; the xbart handle path
(createFromHandle) and updatePredictorPerObservationJointly; Windows. One unprobed portability
note, not a finding: ShippedDrawHook::callProtected throws UnwindJump out of R_UnwindProtect's
cleanup, so the C++ unwind crosses R_UnwindProtect's own C frame in libR; it works on macOS arm64
here (probed), but Rcpp and cpp11 deliberately longjmp to their own setjmp instead, so a platform
whose libR lacks unwind tables would terminate.

---

## bridge-01 - BLOCKER - a negative burn-in returns uninitialized memory as draws

Location: src/R_interface_bartcore.cpp bartcore_run (also bartcore_runWithCallback),
`static_cast<size_t>(Rf_asInteger(numBurnInExpr))`; R/bartcore.R bartcoreSamplerRun and
R/dbarts.R resolveRunCount, which check neither sign.

Claim: `sampler$run(numBurnIn = -3L, 5L)` succeeds and returns uninitialized heap memory as
sigma and fitted draws, and the saved-tree store then reports draws that were never recorded,
so `predict` serves the zero-leaf slots that refuseEmptyTreeStore exists to refuse.

Probe (r3-bridge-p2.R and inline):

```r
s <- dbarts(x, y, control = dbartsControl(n.chains = 1L, n.threads = 1L,
  keepTrees = TRUE, n.samples = 5L, verbose = FALSE))
r <- s$run(-2L, 5L); r$sigma
# [1] 4.880590e-313 6.153788e-313 6.153788e-313 6.153788e-313 8.063584e-313
s$predict(x)[1:2, ]
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.4255184 0.4255184 0.4255184 0.4255184 0.4255184   (5 "draws", all the empty store)
```

Mechanism: -2 becomes 2^64 - 2; Chain::runSweeps computes (numBurnIn + numSamples) * numThin,
which wraps to 3 sweeps, all inside the burn-in, so nothing records; Sampler::run then advances
recordedDraws_ by numSamples. A negative numSamples instead fails with R's "negative length vectors
are not allowed" from Rf_allocVector, and run(0L, 0L) silently returns NULL.

Regression: 0.9-34 refuses ("number of burn-in steps must be greater than or equal to 0",
probe r3-bridge-p16.R), and its run also refused 0 + 0.

Why gates missed: dbartsControl validates n.burn and n.samples, so every test that runs through
the defaults passes; no test passes an explicit negative count to $run.

Fix: refuse negative and NA counts at the boundary (in bartcore_run, as 0.9-34's rc_getInt did)
and in resolveRunCount, with the 0.9-x message.

## bridge-02 - BLOCKER - flat run and predict overrun the caller's buffers on a multinomial handle

Location: src/C_interface.cpp dbarts_sampler_run (`engineResults.numReportedLocations =
shape.numReportedLocations`) and dbarts_sampler_predict; inst/include/dbarts/dbarts.h
dbarts_results and dbarts_sampler_predict layout docs.

Claim: on a handle to an R-built `family = "multinomial"` sampler, which the header says a handle
may name (dbarts_sampler_family is "total over every sampler any construction path can build" and
reports DBARTS_FAMILY_MULTINOMIAL), dbarts_sampler_run writes numObservations x K x numSamples x
numChains into `train` and `test` and dbarts_sampler_predict writes numRows x K x ... into `out`,
K times the sizes the header documents, a heap overflow of the consumer's buffers.

Probe (r3-bridge-c4.R, r3-bridge-c5.R; buffers sized per the header with a NaN canary tail of
the same size):

```
gaussian       family=1 canary changed: sigma=0 train=0 test=0 varcount=0
multinomial    family=8 canary changed: sigma=0 train=180 test=30 varcount=0   (whole canary)
bcf            family=1 canary changed: sigma=0 train=0 test=0 varcount=0
ordinal        family=5 canary changed: sigma=0 train=0 test=0 varcount=0
hetero         family=1 canary changed: sigma=0 train=0 test=0 varcount=0
multinomial keepTrees = FALSE : status 1  documented size 8  canary entries overwritten 8
multinomial keepTrees = TRUE : status 1  documented size 24  canary entries overwritten 24
```

The run comment says the stride is "1 for every dbarts.h-created sampler, since the flat C API
builds no multi-location model", but the flat API builds no sampler at all: every handle is an R
object's pointer, multinomial included. predict's only refusal is testFitsAreUndefined, and the
multinomial blend is defined.

Why gates missed: inst/tinytest/test-capi.R and capi/consumer.c drive only gaussian and binary
handles; no gate runs a multinomial sampler through a flat entry.

Fix: refuse multi-location samplers in dbarts_sampler_run (a raise, as it has no status) and
return 0 from dbarts_sampler_predict when shape.numReportedLocations > 1, documenting both;
or state the n x K x S layout in the header and size by it.

## bridge-03 - MAJOR - a refused installTrees has already replaced earlier chains

Location: src/bartcore/sampler.hpp Sampler::installForests (the per-chain
`chains_[c]->installForest(...)` loop) and Chain::installForest (restoreScale before the
forest rebuild); src/R_interface_bartcore.cpp bartcore_bridge::installForests, which maps the
result to the shapeMismatch message.

Claim: when a cross-grid donor's tree fails to rebuild on chain c > 0, installTrees raises a
refusal but chains 0..c-1 (and chain c's restored response scale) already carry the donor's
trees, contrary to the engine's "On any mismatch nothing is touched" (sampler.hpp, above
installForests); the refusal names "number of trees, forests, or predictors differ", which is
not the cause.

Probe (r3-bridge-p13.R): a donor whose x was changed by setPredictor plus setData, installed
into a fresh 2-chain sampler with `samples = c(1L, 3L)`:

```
warm-start donor is not shape-compatible with this sampler (
live predictions unchanged after refusal: FALSE
live trees unchanged: FALSE
           [,1]        [,2]       [,3]        [,4]      <- before chain 1, chain 2 | after
[1,] -1.0642017 -0.19920680 -0.5436500 -0.19920680        chain 1 replaced, chain 2 kept
```

Why gates missed: the refusal tests check the message, not the sampler's state afterwards, and
use donors that fail the up-front checks (tree counts, DART, grid), which do run before any
mutation; only the rebuild failure arrives mid-loop.

Fix: rebuild every chain into scratch (or snapshot and restore) before committing any, and give
the rebuild failure its own WarmStartResult and message.

## bridge-04 - MAJOR - C++ exceptions still escape several entries and abort the session

Location: src/R_interface_bartcore.cpp bartcore_create / createHolder (bartcore::createSampler),
bartcore_setControl (setTreeStorage), bartcore_growFromRoot (no captureExceptions; its threaded
arm in Sampler::growFromRoot has no worker catch), and src/C_interface.cpp
dbarts_sampler_setTreeStorage; Sampler::fanOutPredictSlabs and Sampler::run spawn threads with
`workers.emplace_back` outside any guard.

Claim: an allocation failure (or a thread-creation failure) in these paths leaves a C++ exception
crossing extern "C" into R and the process terminates, though R_interface_bartcore_common.hpp
says of captureExceptions that "every entrance the engine can throw through therefore ends in
this".

Probe (r3-bridge-p15.R):

```
dbarts(x, y, control = dbartsControl(n.trees = 200000L, n.samples = 2000000000L, keepTrees = TRUE, ...))
libc++abi: terminating due to uncaught exception of type std::bad_alloc: std::bad_alloc
exit=134
s$setControl(ctl)   # keepTrees = TRUE, n.samples = 2e9 on an existing sampler: same abort, exit=134
```

For the spawn paths no probe was run (the brief caps n.threads at 2): a std::thread constructor
that cannot start a thread throws std::system_error, and unwinding past `workers` destroys
joinable threads, which calls std::terminate. predict accepts any n.threads >= 1 and starts
min(n.threads, recorded draws x chains) workers, so `predict(x, n.threads = 10000L)` on a fit
with 10,000 saved draws exceeds macOS's 8,192 threads per task.

Not a regression for allocation (0.9-34 aborts on the same creation call, r3-bridge-p16.R); on
Linux, without macOS's overcommit, an oversize fit fails this way at sizes that are merely large.

Why gates missed: no test drives an allocation failure outside the run path; the monotone
failNextCount hook covers only order counts.

Fix: wrap createSampler, setControl's setters, growFromRoot and the flat setTreeStorage in
captureExceptions; catch spawn failures inside fan-out and run (join what started, then rethrow).

## bridge-05 - MAJOR - setResponse and setOffset accept infinities that poison the sampler for good

Location: R/bartcore.R bartcoreSamplerSetResponse and bartcoreSamplerSetOffset (NA checks only);
src/R_interface_bartcore.cpp bartcore_setResponse and bartcore_setOffset (no finiteness check
for gaussian); src/C_interface.cpp validateResponseSupport, which constrains nothing for
gaussian.

Claim: creation refuses a non-finite response ("response contains non-finite values"), but the
conduits accept one, and a single Inf leaves sigma and every fit NaN even after the response or
offset is put back.

Probe (r3-bridge-p19.R, r3-bridge-p20.R):

```
create with Inf: response contains non-finite values
after setResponse(Inf): sigma NaN NaN  train[1:3] NaN NaN NaN
after restoring y: sigma NaN NaN               (100 more sweeps)
after restoring offset: sigma NaN NaN
```

0.9-34 behaves the same on the conduits, so this is not a regression, but the R layer already
refuses NA there, and an embedded Gibbs loop that feeds one overflowed value loses the chain
silently.

Why gates missed: conduit tests check NA and length only.

Fix: refuse non-finite values in both conduits (R and bridge) with the creation message; a
non-finite offset likewise.

## bridge-06 - MINOR - a refused setOffset leaves data@offset changed

Location: R/bartcore.R bartcoreSamplerSetOffset (`sampler$data@offset <- offset` before the
.Call, with no restore on error; setWeights and setTestOffset restore).

Claim: on a multi-forest sampler, setOffset with an updateScale that is neither TRUE nor FALSE
(NA, 1) skips the R pre-check (isTRUE), is refused by refuseMultiForestResponseMutation, and
leaves the mirror holding the new offset, so a save/load re-creates a sampler conditioned on an
offset the live one never had. On a single-forest sampler, updateScale = 1 rescales in the
engine (`updateScale == TRUE`) while the R side skips reissueNamedLeafSd.

Probe (r3-bridge-p22.R):

```
setOffset(updateScale = NA): $setOffset: a multi-forest sampler supports an offset swap only with updateScale = FALSE, ...
offset after refusal: 5 5 5
live vs reloaded train mean: 0.5406754 0.5999214
```

Why gates missed: the BCF refusal tests use updateScale = TRUE, which R refuses first.

Fix: validate updateScale as a single TRUE/FALSE in the R methods, and assign the mirror after
the .Call succeeds.

## bridge-07 - MINOR - flat predict does not check CSC structure

Location: src/C_interface.cpp translateSource and validateTestSource.

Claim: a CSC source with a row index outside [0, numRows) is replayed silently (wrong values,
out-of-bounds reads), and a column pointer past the stored entries crashes, though the header
says the entry refuses a source "whose declared shape disagrees".

Probe (r3-bridge-c3.R, a valid source then three broken ones):

```
mode 0 returned; matches R predict: TRUE
mode 1 (row index n + 100000) returned; matches R predict: FALSE
mode 2 (last column pointer past nnz): *** caught bus error, exit=138
mode 3 (row index -5) returned; matches R predict: FALSE
```

The flat API validates partially by design and no current consumer passes CSC, hence MINOR.

Why gates missed: capi tests build CSC only from well-formed dgCMatrix data.

Fix: check monotone pointers, last pointer equal to the entry count, and row indices in range
in translateSource (O(nnz)), or say in the header that CSC structure is unchecked.

## bridge-08 - MINOR - flat-set tree storage is lost on re-creation and blocks the restore

Location: src/C_interface.cpp dbarts_sampler_setTreeStorage (and setNumThreads, setVerbose);
R/dbarts.R getPointer, which re-creates from control, model, data and state.

Claim: storage turned on through the flat API is not mirrored into the R object's control, so a
state stored afterwards cannot be reinstalled after save/load.

Probe (r3-bridge-c6.R): keepTrees = FALSE sampler, flat setTreeStorage(1, 5), R run, storeState,
saveRDS/readRDS:

```
live predict dims: 3 5
reloaded predict: state is not consistent with this sampler
```

stan4bart restores through its own state.bart, so it is not affected; a generic consumer is,
and the header's handle section does not say which flat settings survive a re-creation.

Why gates missed: no test mixes a flat setter with an R-side save and load.

Fix: document it in dbarts.h, or have setState adopt the state's store capacity.

## bridge-09 - MINOR - misleading messages at the R boundary

Location: R/dbarts.R getTrees and setPredictor (column and chain index checks).

Claim: `s$getTrees(chainNums = NA_integer_)` and `s$setPredictor(x[, 1], NA_integer_)` fail with
R's "missing value where TRUE/FALSE needed" rather than naming the argument; run(5L, -1L) fails
with "negative length vectors are not allowed" (see bridge-01).

Probe (r3-bridge-p1.R):

```
setPredictor col NA                      ERR: missing value where TRUE/FALSE needed
getTrees chain NA                        ERR: missing value where TRUE/FALSE needed
```

Why gates missed: the range tests use out-of-range integers, not NA.

Fix: refuse NA indices by name before the range comparison.

---

## Checked and found correct

- Draw callbacks, flat and R route: a raise and a C++ throw from an inline callback become R
  errors with the sampler usable afterwards; a worker-thread throw is carried to the main thread;
  a nonzero return stops the run with "stopped by the callback". 60,000 raising runs on each route
  left the protect stack balanced, and gctorture over raise and throw was clean.
- gctorture(TRUE) over creation, run, predict, getTrees, store/setState, copy, every setter,
  setModel/setControl/setData, printTrees, getLatents; the BCF readers, forest weights and
  predictForests; multinomial setCounts/setCategoryOffset; ordinal latents; the heteroscedastic
  predict and getVariance. No faults.
- dbarts_sampler_destroy is idempotent; R methods afterwards re-create from the stored state, or
  refuse with the storeState message when there is none.
- save/load: predict is identical after a reload; active rows, offset, response and sigma all
  carry over.
- R-boundary validation of about 45 malformed inputs (lengths, NA, types, column and tree
  ranges, thread counts, unsorted cut points, zero-row predict) refuses with a message.
- ABI: the layout and token static_asserts match the header; DBARTS_RESULTS_INIT and
  DBARTS_PREDICTOR_SOURCE_INIT list every field; the stubs build clean under DBARTS_REQUIRE_EXACT_ABI.
  Every entry stan4bart bartcore calls exists with the signature it uses.
- setTimeLimit at 1 and 2 threads stops a run with "sampler run interrupted" (as documented). The
  sampler runs and predicts afterwards, and the cancelled run does not advance the recorded-draw
  count.
- Flat run on gaussian, BCF, ordinal and heteroscedastic handles writes nothing past the
  documented sizes (bridge-02's canary).
- setControl with an unchanged keepTrees does not discard the saved trees.
