# per-draw-callbacks

agent: opus (S1 the engine hook, the per-draw stride, the cancel plumbing
  and the bridge's draw fill; S2 the flat C header, its entry and the ABI
  hash re-bake); sonnet (S3 the R surface and its Rd, S4 the example, its
  test and the manual, S5 the sister-package record), with an opus read
  of S4's C source before it is committed - a defect there crashes a
  user's session with no condition to catch.
rng: NEUTRAL in every slice. The hook consumes no generator draw, writes
  no sampled state and moves no default, so a run with no callback
  registered is bitwise what it is today - the property
  [`SweepCallback`](../../src/bartcore/chain.hpp) and
  [`shouldCancel`](../../src/bartcore/sampler.hpp) already carry and S1's
  own tests/cpp assertion pins directly. Gate in every slice: the
  equivalence trio IDENTICAL on every scenario - 52 gaussian, 12 BCF, 11
  multinomial, against the f0236082 baselines, on the shipped build - no
  scenario added, no baseline re-recorded, exact-posterior gates untouched.
window: serial. S1 before S2 (the header's struct mirrors the fill S1
  writes), S2 before S3 (the R argument passes what the header declares,
  and S3's tests reuse S2's compiled counting consumer), S3 before S4 (the
  recipe calls the shipped argument); S5 last and docs-only. S1 re-blocks
  [`Chain::run`](../../src/bartcore/chain.hpp) and
  [`Sampler::run`](../../src/bartcore/sampler.hpp), which
  docs/plans/engine-performance.md's S4 and S5 also rewrite, so those do not
  run beside it.
budget: S1 ~90 engine + ~70 bridge + ~50 tests/cpp; S2 ~50 header + ~15
  src/C_interface.cpp + ~40 tests; S3 ~110 R + ~50 Rd + ~80 tinytest + ~10
  NEWS; S4 ~70 vignette + ~60 tinytest C + ~20 manual and memory note; S5
  ~25 docs. About 740 against the note's aggregate estimate of 550
  ([8. Scope and sequencing](../design/per-draw-callbacks.md#8-scope-and-sequencing));
  the difference is `keepFits` over four channels with its family refusals,
  and the plain-C test copy fork 6 asks for, neither itemized there.

Decisions: dec-B114 (the whole arc: motivation, storage opt-out, header
entry, stop flag, worker-thread observer) in
[docs/decisions.md](../decisions.md). It reverses dec-B62 for an OBSERVER
hook only, re-opens the surface dec-B86 trimmed with two named consumers,
and is governed by dec-B111's version-pair rule at the hash re-bake.

## Goal

A per-draw C callback fires once per saved draw per chain, on the thread
owning that chain, over const pointers into engine buffers, so a caller
reduces draws as they are produced instead of materializing the n x draws
x chains array the memory audit found to be five sixths of a large fit's
peak. `bart()`, `dbarts()` and the sampler's run method take it as a pair
of external pointers; a new `keepFits` control stops keeping the
per-observation channels when one is supplied. The flat C header ships the
callback signature, the draw struct and the setter, so a compiled consumer
registers one without hand-declaring an unversioned layout.

## Context

The argument, the priced memory consequence and the worked example are in
[Per-draw callbacks](../design/per-draw-callbacks.md#per-draw-callbacks),
whose figures come from
[Reference cases](../design/memory-footprint.md#reference-cases). What the
implementer needs from the code:

- The firing point. [`storeSample`](../../src/bartcore/chain.hpp) writes
  each kept draw into the caller-owned
  [`Results`](../../src/bartcore/chain.hpp) slabs the bridge sizes
  ([`allocChannel`](../../src/R_interface_bartcore.cpp),
  [`installChannel`](../../src/R_interface_bartcore.cpp)), so a hook fired
  right after it sees every channel the fit carries; a sweep discarded as
  burn-in inside one `run(numBurnIn, numSamples)` call never reaches it.
- The stride is in the sampler, not the bridge.
  [`Sampler::run`](../../src/bartcore/sampler.hpp) computes each chain's
  base as `results.trainingFits + c * numSamples * n * numLocations` before
  any draw is stored, and `Results` is one struct it slabs, so the bridge
  cannot point one chain's buffer elsewhere and a one-draw buffer would
  have chains 1..C-1 writing past its end. `Results` is engine-internal: a
  new field is additive, not an ABI event.
- `Sampler::run` clamps
  [`numVariableCountForests`](../../src/bartcore/sampler.hpp) to what the
  combiner can report before striding by it, and the bridge accumulates
  `uint32_t` counts before widening them to R integers.
  [`testFitsAreDefined`](../../src/bartcore/combiner.hpp) and
  [`logLikelihoodIsDefined`](../../src/bartcore/combiner.hpp) gate two
  channels `storeSample` NaN-fills rather than skips; the R channels keep
  that fill, their shape being part of the returned object.
- No stop plumbing exists today: the worker lambda in `Sampler::run`
  discards `Chain::run`'s return value,
  [`cancelFlag`](../../src/bartcore/sampler.hpp) is written only by the
  main thread's [`pollInterrupt`](../../src/bartcore/sampler.hpp), and a
  worker holding several chains moves on regardless. The precedent for
  telling a self-caught stop from an interrupt is
  [`closureStopped`](../../src/R_interface_bartcore.cpp) beside `cancelled`
  in [`bartcore_runWithCallback`](../../src/R_interface_bartcore.cpp).
- Burn-in is an R-layer problem: [`runWithBurnIn`](../../R/bart.R) passes
  no burn-in count to the engine, running `sampler$run(0L, control@n.burn)`
  then `sampler$run(0L, control@n.samples)`, so every burn-in sweep is a
  recorded draw there.
- [`keepTrainingFits`](../../R/A_class.R) is the only per-observation
  switch: the test channel is gated on `numTestObservations > 0` and the
  heteroscedastic
  ([`varianceTrainExpr`](../../src/R_interface_bartcore.cpp)) and
  multi-forest ([`forestFitsExpr`](../../src/R_interface_bartcore.cpp))
  channels only on the model carrying them, so no argument reaches the two
  largest cases. [`checkFamilyUnsupportedArgs`](../../R/bart.R) already
  refuses `keepTrainingFits = FALSE` by name for four families.
- The ABI token folds the `DBARTS_C_API_LIST` signatures, the three enums
  and the layouts of the structs that cross the ABI;
  [`dbarts_apiToken`](../../src/C_interface.cpp) static_asserts against the
  baked [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h), so a
  type folded in re-bakes the token and one left out is an unversioned
  layout. [`tools/check-api-hash.sh`](../../tools/check-api-hash.sh)
  enforces dec-B111's pair rule and prints "no release tag, skipped".
- The compiled-consumer precedent under test is
  ["consumer source not installed"](../../inst/common/capiConsumer.R): plain
  C through `R CMD SHLIB`, self-gating on the toolchain. Rcpp is not a
  dbarts dependency.

## Decision

Every fork the design note left open is settled. VD's rulings, 2026-09-10,
verbatim:

- The arc's motivation: "That strikes me as an argument for more and
  better callback support." / "we should bump up the priority. I also
  think it would be good to have an Rcpp example showing how to have a C
  callback write to a preallocated array, if that makes sense. That would
  be the default way to do it in R, since that wouldn't require blocking.
  The callback would of course use raw pointers or lists or void\*, just
  not SEXPs."
- Fork 4 (storage opt-out): a new logical `keepFits` over every
  per-observation channel (training, test, variance, per-forest), default
  TRUE, set FALSE automatically when a callback is supplied unless the
  user overrides; `keepTrainingFits` stays as the narrower existing
  switch; variable counts and scalar channels are always kept; `keepTrees`
  stays the recompute path. This mirrors stan4bart's `keep_fits`
  ("Intended to be used with callback"). VD: "Sounds good."
- Fork 2 (registration): VD: "Ship the header entry and its two types
  now." The flat C header gains the setter entry and the two types
  (callback signature, draw struct), one ABI hash re-bake, dec-B86 argued
  as the note does - dbarts's own R layer, the documented example and
  stan4bart are the consumers.
- Fork 5 (stop): VD: "Sure, option 1." The callback returns `int`, nonzero
  aborts the run: a shared cancel flag any worker can set, the worker loop
  reading it, the bridge distinguishing a callback stop from an interrupt,
  the sampler's inconsistent-after-abort state documented.
- Forks 1 and 3 are settled by those words: the hook runs on worker
  threads as a const-pointer observer forbidden to call R (reversing
  dec-B62's premise, argued in the note), and no blocking R-closure
  variant ships in this arc. If one is ever revived, VD's shape for it is
  that the machinery calls the closure once to learn its return length,
  preallocates, and copies each draw's result in itself (stan4bart's
  `callbackResultLength` pattern) - the design for a post-release item,
  not this arc.

Forks 6 through 9 are decided by the orchestrator on the note's
recommendations, and are agent-made: the example is a recipe in
[dbarts-as-a-component.Rmd](../../vignettes/dbarts-as-a-component.Rmd) (Rcpp
form, unevaluated) plus a plain-C copy under `inst/tinytest` the suite
compiles; `xbart` does not get the argument ([`xbart`](../../R/xbart.R)
holds loss values, not fits); the hazard family's expanded-row meaning for
`numObservations` is documented only; and undefined channels are handed
over as null pointers, not NaN buffers.

## Constraints

- Gates per slice are in Verification below; the class-wide ones are
  `tests/cpp`, the full tinytest suite and the equivalence trio at
  52/12/11, plus `pkgdown::check_pkgdown(".")` if S3 adds an Rd topic and
  ASan over both `tests/cpp` and the R-loaded path for S1 and S2, the hook
  making new pointer arithmetic reachable from a `.Call`.
- Contract freezes: `dbarts_draw` is size-first like `dbarts_results` and
  its fields append monotonically; pre-release the version constants stay
  at 1/0 with one hash re-bake (dec-B111). The engine never sees the
  shipped header - `Results` and the engine-side draw struct stay
  engine-internal and the entry files convert.
- No R API inside the callback, ever - not `Rf_allocVector`, not
  `PROTECT`, not `Rf_error`, and in C++ not the construction or
  destruction of an Rcpp proxy type. Stated in the header, the Rd and the
  vignette, in those terms.
- Out of scope: the blocking R-closure variant (fork 3); any per-draw tree
  handoff; `xbart` (fork 7); a family-specific expansion index on the
  struct (fork 8); an engine-held lock around the callback; and any change
  to `bartcore_runWithCallback`, which keeps its inline-only refusal as a
  CONDITIONING hook.

## Steps

S1, engine and bridge (src/bartcore/chain.hpp, src/bartcore/sampler.hpp,
src/R_interface_bartcore.cpp, tests/cpp/test_sampler.cpp; landed also
touching src/bartcore/facade.hpp, src/R_interface.cpp, R/bartcore.R and
tests/cpp/test_facade.cpp, the virtual `run` signature and its two R
call sites):

1. Give `Chain::run` a draw callback beside its sweep callback and cancel
   predicate: a function pointer and a context pointer, both null by
   default, over an engine-internal POD draw struct. Fire it once per SAVED
   draw, immediately after `storeSample` settles it, with a draw index
   counting this call's saved draws from 0 and the chain index the chain
   already holds; a sweep discarded as burn-in fires nothing.
2. Fill the struct per channel from `Results` and the chain's own shape:
   every channel `storeSample` settles, not a selection. A null pointer
   wherever the fit does not carry the channel and wherever the combiner
   declares it undefined (`testFitsAreDefined`, `logLikelihoodIsDefined`),
   NaN for an inapplicable scalar, and a varcount pointer covering
   `numPredictors * numVariableCountForests` forest-major within a draw, at
   the count `Sampler::run` clamped.
3. Add a per-draw stride to `Results`, defaulting to `numObservations *
   numReportedLocations` and read BOTH by `Sampler::run`'s per-chain base
   arithmetic and by `storeSample`'s per-draw offset, with the parallel
   field for the test, variance and forest channels; a stride of zero lands
   every draw of a chain in that chain's one buffer. The comment states why
   it is a field: the two readers must agree, which a computation repeated
   at each site does not guarantee.
4. Forward the callback from `Sampler::run` on BOTH paths, inline and
   spawned-worker, unchanged in either: no lock, no queue, no main-thread
   hop. Calls for different chains may run concurrently, calls within a
   chain are ordered by draw index, and the comment says why no mutex - it
   would make every chain wait on the slowest callback.
5. Make the stop real: one cancel flag shared by the interrupt arm and the
   callback arm, the worker lambda storing into it when `Chain::run`
   returns true, each chain reading it at its next sweep boundary. Document
   that a stop is an ABORT - `Sampler::run` returns before advancing the
   sample cursors while saved tree records are already in the slots those
   cursors count, so results and saved trees are discarded, as on the
   interrupt path.
6. Bridge: a per-run function-pointer and context argument on
   `bartcore_run` read out of two external pointers, and an adapter
   converting the engine's draw struct into the shipped `dbarts_draw` (S2)
   before the call, so the engine stays header-agnostic at a cost of one
   small stack copy per draw. Tell a callback stop from an interrupt on
   return, as `bartcore_runWithCallback` tells `closureStopped` from
   `cancelled`, and raise a different condition for each.
7. Bridge: honour a `keepFits` flag on the control the run entry already
   reads (its R slot arrives in S3; default TRUE here). When FALSE,
   allocate C one-draw buffers CONTIGUOUSLY per opted-out per-observation
   channel as scratch, set that channel's stride to zero and return no R
   array for it; counts and scalars are allocated as today.
8. tests/cpp: the hook fires exactly once per recorded draw per chain and
   in draw order; a multi-chain run under more than one worker gets its
   per-chain counts and a chain-disjoint write pattern; a stride of zero
   leaves each chain's scratch buffer holding that chain's last draw and no
   chain writing outside its own; a nonzero return stops the run; and a
   seeded no-op callback produces draws bitwise identical to no callback at
   all - the neutrality pin.

S2, the flat C header (inst/include/dbarts/dbarts.h, src/C_interface.cpp,
tests/cpp, inst/tinytest/test-capi.R; landed also touching
src/R_interface_bartcore_common.hpp, src/R_interface_bartcore.cpp,
tests/cpp/Makefile, tests/cpp/main.cpp and tests/cpp/common.hpp):

9. Declare `dbarts_draw` - size-first `structSize`, the indices, the shape
   counts, the channel pointers and the scalars - and
   `dbarts_draw_callback`, `int (*)(void* context, const dbarts_draw*
   draw)`. The header states: nonzero aborts the run; a null pointer means
   the channel is absent and a NaN scalar means inapplicable; pointer
   validity is the call and no longer, an opted-out channel's pointer being
   scratch the chain's next draw overwrites; calls run concurrently across
   chains with no engine lock, and why; and the no-R-API rule in the three
   terms of the Constraints above, the Rcpp proxy-type ban included.
10. Add `dbarts_sampler_setDrawCallback(sampler, fn, context)` to
    `DBARTS_C_API_LIST`, null `fn` clearing. It is a setter that copies
    into the sampler like every other, so the list's existing pure-C rules
    carry over unchanged.
11. Re-bake the token: update the literal in `dbarts_apiToken` and
    `DBARTS_C_API_HASH` together, the existing static_assert being the gate
    that they agree. Under dec-B111 the version pair stays at 1/0
    pre-release and `tools/check-api-hash.sh` still prints "no release tag,
    skipped", so nothing else moves.
12. Tests: a tests/cpp contract test that the entry registers, clears on
    null and survives a second registration; and a tinytest C API arm, a
    plain-C consumer compiled and self-gated the way the existing capi
    consumer is, registering a COUNTING callback and asserting the call
    count and per-chain indices. S3 reuses it, so give it a per-chain
    counter and a status field readable from R.

S3, the R surface (R/A_class.R, R/dbarts.R, R/bart.R, R/generics.R,
R/plot.R, man/, inst/NEWS.Rd, inst/tinytest; landed also touching
R/bartcore.R and inst/tinytest/test-argument-surface.R):

13. `callback` on `bart()`, `dbarts()` and the sampler's `run` method: a
    list of two external pointers, `fn` and `context`. Check that both are
    external pointers and pass the addresses down; nothing else is checked,
    and the Rd says why (the address is dereferenced as handed). `bartBT`
    does not get it - its formals are 0.9-34's.
14. `keepFits` on [`dbartsControl`](../../R/A_class.R): a length-1 logical
    slot, prototype TRUE, its validity check beside `keepTrainingFits`, and
    the matching argument on [`dbartsControl()`](../../R/dbarts.R) and
    `bart()`. FALSE drops the training, test, variance and per-forest
    channels; variable counts, scalars and `keepTrees` are unaffected.
15. The automatic default: when `callback` is supplied and the caller did
    not name `keepFits`, `keepFits` becomes FALSE; an explicit value always
    wins. Drive it off `missing()`/the matched call, not a sentinel, so
    `keepFits = TRUE` with a callback is honoured.
16. `checkFamilyUnsupportedArgs` rejects `keepTrainingFits = FALSE` for
    multinomial, ordinal, nbinom and hurdle.lognormal, so an AUTOMATIC
    `keepFits = FALSE` must raise there by name rather than silently
    breaking the fit - the same message shape, naming `keepFits` and the
    callback that set it.
17. `runWithBurnIn` installs the callback on the KEPT-sample run only,
    leaving the engine rule ("fires on saved draws") as stated. The comment
    states the constraint: every burn-in sweep is a recorded draw at this
    layer, so a hook installed for both calls would average burn-in in.
18. The returned object under `keepFits = FALSE`:
    [`packageBartResults`](../../R/bart.R) omits the absent channels and
    the means taken from them, and [`plot.bart`](../../R/plot.R),
    [`extract`](../../R/generics.R), `fitted`, `residuals` and `predict`
    each name the absent channel and the argument that dropped it rather
    than failing on a NULL. The default callback fit returns no fitted
    value at all - the accumulator is the answer - and the Rd says so.
19. Rd and NEWS: `callback` and `keepFits` on bart.Rd, dbartsControl.Rd and
    the sampler class Rd, carrying three warnings - a mistake in the
    callback crashes the session, an interrupt cannot land while a call is
    running, and the default drops the fitted values - plus the hazard note
    (fork 8: `numObservations` counts person-period rows, not subjects; the
    mapping is R-side). One NEWS entry.
20. tinytest: argument validation and its messages; the control slot, its
    validity and the automatic default including the explicit override; the
    four family refusals; the returned object under `keepFits = FALSE` with
    no callback; and, reusing S2's counting consumer, exactly `n.samples`
    calls at the `bart()` defaults - the burn-in regression - plus
    per-chain counts under `n.threads > 1` and a stop-flag run.

S4, the example, its test and the manual
(vignettes/dbarts-as-a-component.Rmd, inst/tinytest, man/dbarts-package.Rd,
docs/design/memory-footprint.md, benchmarks/; landed also touching
inst/common/capiConsumer.R, docs/design/per-draw-callbacks.md and
docs/plans/pure-c-header.md):

21. A seventh recipe in the component vignette: the running posterior mean
    of
    [5. The Rcpp example](../design/per-draw-callbacks.md#5-the-rcpp-example),
    Rcpp form, unevaluated, with `R_MakeExternalPtrFn` for the function
    pointer (casting one to `void*` is what `-Wpedantic` flags under `R CMD
    check`) and the accumulator in the external pointer's protected slot.
    The prose carries the lifetime rules: read the accumulator once as a
    raw pointer before the run, chains write disjoint slices, a context is
    per run because the draw index restarts, and nonzero ABORTS - so a
    callback that merely disagrees with a draw records a status, returns 0.
22. A plain-C copy of that callback under `inst/tinytest`, compiled and
    self-gated the way the existing capi consumer is, asserting the
    accumulated mean equals the same seeded fit's `yhat.train` mean to
    floating-point tolerance - the test that the example is correct, not
    merely that it compiles.
23. Manual: one sentence in the Memory section of
    [`dbarts-package`](../../man/dbarts-package.Rd) naming the callback
    route beside `keepTrainingFits`, with the reference-case numbers (3891
    MB to 698 MB at n = 1e5, p = 20, C = 4, S = 500; 10560 MB to 2576 MB at
    n = 1e6, p = 50, C = 1, both re-derived from the post-ingestion-guard
    rows); and the two rows
    [Reference cases](../design/memory-footprint.md#reference-cases) owes -
    the heteroscedastic variance channel and the BCF per-forest fits,
    neither reachable by `keepTrainingFits` today.
24. The per-draw cost the note owes: a bench-sampler scenario with a no-op
    callback and one with the example's, at both reference shapes, so
    [7. Threading interaction](../design/per-draw-callbacks.md#7-threading-interaction)
    - a fraction of a percent at small n, unasserted at large n - is
    measured rather than argued. Maintainer-run on a quiet machine; the
    scenario and its numbers land with this slice, the compare does not
    gate it.

S5, the sister-package record (docs/design/per-draw-callbacks.md):

25. Record stan4bart as the second consumer: nothing to port - it calls
    neither the new type nor the new entry, and the unchanged version pair
    admits its existing binary - while its per-iteration caller-owned
    buffer pattern, one draw at a time with chains in separate single-chain
    samplers, is the shape this entry serves, so it may adopt it
    post-release. Name its `keep_fits` as the precedent `keepFits` mirrors.

## Verification

Run from the slice's own worktree and private library, each gate on its
own exit status.

- `cd tests/cpp && make && ./test_bartcore` - all pass, S1's five new
  assertions included. Re-run for S1 and S2 under `OPT="-O2 -g
  -fsanitize=address,undefined"` with
  `ASAN_OPTIONS=detect_container_overflow=0`; a symbolization prompt means
  a diagnostic fired, so read the count.
- `R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'` - full suite
  passes, the compiled arms reporting run rather than skipped on this
  machine (a skip is not evidence; check the file's own gate line).
- The equivalence trio against the f0236082 baselines on the shipped build:
  `equivalence.R` 52 of 52, `bcf-equivalence.R` 12 of 12,
  `multinomial-equivalence.R` 11 of 11, counted as per-scenario "identical
  draws (same RNG stream)" lines with zero "max |z|" lines and zero
  skipped. Any deviation is a defect in the change.
- `Rscript tools/check-doc-freshness.R` and `Rscript
  tools/check-rc-codoc.R` exit 0 (S3 adds a reference-class formal);
  `air format --check .` tree-wide and `lintr::lint_package()` after S3;
  `R CMD check --as-cran` from a tarball staged outside the tree after S3
  and after S4.
- S2 only: the `dbarts_apiToken` static_assert compiles, which is the
  proof the re-baked hash matches the folded layout, and
  `tools/check-api-hash.sh` still reports "no release tag, skipped".
- S1's stop path is shown to discriminate: mutate the worker lambda to
  discard `Chain::run`'s return value again, confirm the tests/cpp stop
  assertion fails, then revert and `touch` the file before rebuilding.
- Maintainer-run, not a merge gate: `benchmarks/R/bench-sampler.R compare`
  against benchmarks/baselines/bench-sampler-127f04ee.csv on a quiet
  machine, plus step 24's two callback scenarios.

## Landing note, S1 (2026-09-10)

LANDED at 711942adbccbab4b7084d741bdf9483bcc72d089, three commits:

- b91cabb951feb1327a682620ff866423329d1078 Fire a per-draw observer on
  the chain's own thread
- 0bb5b52a16e94229aef6922151eefe35d6e0d874 Pin the per-draw observer's
  contract in tests/cpp
- 711942adbccbab4b7084d741bdf9483bcc72d089 Cover the storage opt-out
  and refuse an unknown keepFits

[`Chain::run`](../../src/bartcore/chain.hpp) takes a plain-C
[`DrawHook`](../../src/bartcore/chain.hpp) beside its sweep callback and
cancel predicate, firing it once per SAVED draw, immediately after
`storeSample` settles it, on the thread that already owns the chain - a
sweep discarded as burn-in never reaches it. The struct handed over,
[`DrawInfo`](../../src/bartcore/chain.hpp), carries every channel
`storeSample` settles, null wherever the fit does not carry the channel
or the coupling declares it undefined; `DrawCallback` returns `int` and a
nonzero return aborts the run. `Results` gained one per-draw stride
field per per-observation channel, read by both
[`Sampler::run`](../../src/bartcore/sampler.hpp)'s per-chain slab base
and `storeSample`'s per-draw offset so the two cannot disagree; a stride
of zero lands every draw of a chain in that chain's own one-draw buffer.
`Sampler::run` forwards the hook unchanged on both the inline and the
spawned-worker path via `DrawRelay`, which also records whether an
abort came from the callback rather than from an interrupt and reports
that through an out-parameter. Changing the virtual `run` signature
also touched [src/bartcore/facade.hpp](../../src/bartcore/facade.hpp) and its
[tests/cpp/test_facade.cpp](../../tests/cpp/test_facade.cpp) spy - both
beyond the plan's stated S1 file list, along with
[`R/bartcore.R`](../../R/bartcore.R) and
[src/R_interface.cpp](../../src/R_interface.cpp) below.

The bridge's `bartcore_run` now takes six `.Call` arguments: the
pointer, burn-in and sample counts, two external pointers for the
callback function and context, and `keepFits`. `keepFits` is a per-RUN
`.Call` argument, not a `dbartsControl` slot - `R/bartcore.R`'s two
existing call sites pass `NULL, NULL, TRUE`, so today's behaviour is
bitwise unchanged; S3 wires `bart()`/`dbarts()`/the sampler's `run`
method and the control slot to these three arguments.
[tests/cpp/test_sampler.cpp](../../tests/cpp/test_sampler.cpp) (465
lines, new) pins the firing count and order at one, two and four
chains, the stride arm's per-chain one-draw blocks, the stop arm at one
and four chains distinguishing a callback stop from an interrupt, and a
neutrality pin across no-hook, hooked and hooked-with-scratch runs.
[inst/tinytest/test-bartcore-keepfits.R](../../inst/tinytest/test-bartcore-keepfits.R)
(115 lines, new) drives `bartcore_run` directly at two chains, once
heteroscedastic with a test set and once two-forest, asserting the
opted-out channels come back null with the list's shape unchanged and
the sigma draws and variable counts bit-for-bit what a keeping run
drew. Review fix before landing: `Rf_asLogical(x) != FALSE` mapped an
NA `keepFits` to TRUE; 711942ad refuses NA explicitly with `Rf_error`,
since it decides how much the run allocates and which slots come back
null. The slice was reviewed by a second reader and independently
verified before landing.

Gates: tinytest full suite 8243 of 8243, zero failures. Equivalence
trio on the shipped build: gaussian 52 of 52 "identical draws (same RNG
stream)" lines, BCF 12 of 12 and multinomial 11 of 11 "identical (all N
channels: ...)" lines, zero "max |z|" lines and zero skipped across all
three. `tests/cpp` passed whole under both ASan and TSan, and the
compiled `test-bartcore-keepfits.R` arm passed 20 of 20 under ASan; the
reviewer's mutation run shows the ASan gate discriminates - reverting
the storage-opt-out stride guard produced a heap-buffer overflow
AddressSanitizer caught at [`Chain::storeSample`](../../src/bartcore/chain.hpp).
`check-doc-freshness`, `check-rc-codoc` and `check-win-drift` all report
OK. CI ran on 711942ad but those runs were cancelled when 5a05d799
(memory-footprint-audit's follow-ons, a sibling slice) pushed on top;
the CI evidence for this slice is therefore the 5a05d799 run.

Remaining: S2 through S5 (the flat C header, the R surface wiring
`callback` and the control slot, the vignette example and manual, the
stan4bart record) are open, per the plan's Steps.

## Landing note, S2 (2026-09-10)

LANDED at e40631405ea04a8deb5cb1b50418ce81bb5633b0, six commits:

- 886c25374ad34b6bee797a0eda8c1e6c5d9294a5 Stop the draw-callback stop test from chain 0 only so its count is pinned
- 658f9f25756d8246158e75b504e66ed418234a85 Ship the per-draw callback types and the setter in dbarts.h
- 34d2650ad93353709292b26b0841e5f21ba20212 Pin the shipped draw struct and the setter in tests/cpp
- 752dd2d1641626d2618583120c594edade535599 Drive the shipped draw callback from the plain-C consumer
- c5a77936d8465cce93625a866e2d84c3e0904846 Say what the ABI token folds for the per-draw struct
- e40631405ea04a8deb5cb1b50418ce81bb5633b0 Drop the field-order mirror test and the setter's restated abort contract

886c2537 precedes S2: CI failed the draw-stop test on 5a05d799 under
four workers (any chain could trip the abort, leaving chain 0's count
unpinned); chain 0 alone stops now, exactly.

[`dbarts_draw`](../../inst/include/dbarts/dbarts.h) - size-first
`structSize` filled by the library, the indices, shape counts, ten
channel pointers, four by-value scalars, `DrawInfo`'s order, an
optional field read through
[`DBARTS_DRAW_HAS`](../../inst/include/dbarts/dbarts.h) - and
`dbarts_draw_callback` and `dbarts_sampler_setDrawCallback` join
`DBARTS_C_API_LIST`: a copying setter, a null `fn` clears and drops
the context. Token re-baked: `DBARTS_C_API_HASH`
`0x6380bf095d5cae3f`, signature literal `0xb6f41cfcbd996897`, version
pair unchanged at 1/0. `dbarts_draw` is the first ABI struct with
by-value doubles: its scalars fold by name, position and width rather
than pointer-unit offset, so the token is one number on ILP32 too -
the header now says so.

One adapter, [`fillShippedDraw`](../../src/R_interface_bartcore_common.hpp)
and [`ShippedDrawHook`](../../src/R_interface_bartcore_common.hpp),
used by both routes: R's
[`bartcore_run`](../../src/R_interface_bartcore.cpp) (now casts
through `dbarts_draw_callback`) and the flat
[`dbarts_sampler_setDrawCallback`](../../src/C_interface.cpp).
[tests/cpp/test_capi.cpp](../../tests/cpp/test_capi.cpp) (new):
registration, clear, re-registration, the adapter copying every
field. The plain-C
[inst/tinytest/capi/consumer.c](../../inst/tinytest/capi/consumer.c)
gains `capi_draw_function`, `capi_draw_context`,
`capi_draw_reset(stopAfter)` and `capi_draw_report` (per-chain
counts, indices, sigma, status), driven by
[inst/tinytest/test-capi.R](../../inst/tinytest/test-capi.R).

Files beyond the plan's S2 list: src/R_interface_bartcore_common.hpp,
src/R_interface_bartcore.cpp, tests/cpp/Makefile, tests/cpp/main.cpp,
tests/cpp/common.hpp. About 660 added lines against the ~145 S2
budget (header prose, the shared adapter, tests - none itemized
there), reported by the implementer; the reviewer judged it not
padded but cut about 40 lines at landing (e4063140: a field-order
mirror test pinning a convention the code does not depend on, a
restated abort paragraph on the setter's doc).

Gates, the second reader's run on the rebased slice, all foreground:
`R CMD INSTALL --preclean` exit 0; tests/cpp 299 ok lines, all
passed; tinytest 8297 TRUE, 0 FALSE, 165 files, 0 skips, `test-capi.R`
143 tests ran; equivalence against the f0236082 baselines 52 of 52,
BCF 12 of 12, multinomial 11 of 11 "identical draws (same RNG
stream)", zero "max |z|", zero skipped; `check-doc-freshness`,
`check-rc-codoc`, `check-win-drift` exit 0; `check-api-hash.sh` "no
release tag, skipped"; `R CMD check --as-cran` Status OK; mutation
probe: dropping `numObservations` from the adapter fails "capi draw:
the adapter copies every channel and scalar". Review findings fixed
before landing: the header's token paragraph described the fold
wrongly for `dbarts_draw` (fixed in c5a77936); one over-long comment
line.

Two notes to S3: a hook via `dbarts_sampler_setDrawCallback` does not
fire on the R route - `bartcore_run` builds its own per-run hook from
the R callback argument, so the R argument is the R route's only
channel, and S3's Rd says so; the R route's shipped-struct path is
proven structurally (one adapter chain) and is first executed by
S3's tests.

Remaining: S3 through S5 (the R surface wiring `callback` and the
control slot, the vignette example and manual, the stan4bart record)
are open, per the plan's Steps.

## Landing note, S3 (2026-09-11)

LANDED at db3bd084db42141e986c961fc30884400c3983e1, three commits:

- 93b8bf36bd60496c80b764309615ce1b04ccdef1 Wire callback and keepFits
  through the R surface (bart, dbarts, sampler run)
- e98cd22472c911724b39b78a63f2439bbfe0ca7d Name keepFits when a dropped
  per-forest channel is asked for
- db3bd084db42141e986c961fc30884400c3983e1 Say what numObservations
  counts on a hazard sampler's callback

`callback` - two external pointers, `fn`/`context`, validated by
[`validateCallback`](../../R/bartcore.R) and passed to `bartcore_run`
per run - lands on `bart()`, `dbarts()` and the sampler's `run` method,
not `bartBT`. `dbarts()` validates but never persists its `callback`
(it never runs the sampler; a later bare `$run()` must pass one again,
`dbarts.Rd`) - the reviewer confirmed step 13 implies exactly this.
[`dbartsControl`](../../R/A_class.R) gains `keepFits` (validity check
beside `keepTrainingFits`), automatic default `keepFits =
is.null(callback)` resolved through `bart()`'s shared-formal path so an
explicit `TRUE` wins. `checkFamilyUnsupportedArgs` refuses an automatic
or explicit `keepFits = FALSE` for multinomial, ordinal, nbinom and
hurdle.lognormal, naming `keepFits` and `callback`; `runWithBurnIn`
installs the callback on the kept-sample run only.

`packageBartResults` omits the channels `keepFits = FALSE` drops and
their means; `plot`, `extract`, `fitted`, `residuals`, `predict` name
the absent channel and dropping argument rather than fail on a bare
NULL. A new `hasVariance` marker survives `keepFits = FALSE` where
`s.train` does not - a latent defect the implementer found: without it
a heteroscedastic fit run with `keepFits = FALSE` and no saved trees
silently skipped `s(x)` in `predict(type = "ppd")` instead of refusing.
e98cd224: an amplitude-coupled fit that way kept `n.forests` but not
`forestFits`, so `extract(type = "forest")`/`predict` fell to the
engine's basis error, naming neither; `refuseDroppedForestChannel` now
names both. db3bd084 adds the hazard row-counting note
(`numObservations` counts person-period rows, not subjects) to the
sampler class Rd's callback item, matching `bart.Rd`/`dbartsControl.Rd`.

Rd on `bart.Rd`, `dbartsControl.Rd`, `dbartsSampler-class.Rd`: the
three warnings (crash on a bad callback, no mid-call interrupt, the
automatic default drops the fitted values) and the hazard note; one
NEWS entry for the arc. Tests: `test-callback.R` (new), `test-capi.R`
(extended: exactly `n.samples` calls at `bart()` defaults, per-chain
counts under `n.threads = 2`, a stop run), `test-argument-surface.R`
(the two new formals).

Files beyond the plan's S3 list: R/bartcore.R (`validateCallback`),
inst/tinytest/test-argument-surface.R (formal-parity contract). 605
added / 35 removed against a budget of about 250; the reviewer judged
it load-bearing (validation, the latent-defect fix, the per-family
refusal tests, the arc's only NEWS entry) and cut four redundant
assertions.

Gates, the reviewer's run: preclean install exit 0; tinytest 8349
TRUE, 0 FALSE, 166 files, 0 skips, `test-capi.R` 154 and
`test-callback.R` 40 ran; equivalence 52/52, 12/12, 11/11 identical
draws, zero max |z|, zero skipped; no file under src/, inst/include or
tests/cpp touched; doc-freshness, rc-codoc, air, lintr clean; `R CMD
check --as-cran` OK. Mutation: installing the callback on the burn-in
run too fails six test-capi assertions (call counts 500 to 1000).

Remaining: S4 (the vignette example, its test and the manual) is open,
per the plan's Steps.

## Landing note, S4 (2026-09-11)

LANDED at 4cf56b69bca5f1855c856cdb6fbe524b876e4561, twelve commits, the
slice then its review fixes:

- 2ff936265fa6d6b762a7dfb8fdcbc14ed81af985 Drive the vignette's
  running-mean recipe from a plain-C copy under test
- 6a66a2d14750fb0a6063a9d5250fefd725704dd4 Add the per-draw callback
  recipe to the component vignette
- 1256078a0e7bd88e1f7692cd29cefc352baa0c7c Name the per-draw callback
  route beside keepTrainingFits, and the two rows the memory note owed
- 2f43401484f3f34b1f0c52f79b40f55c515a3105 Add no-op and running-mean
  callback scenarios to bench-sampler.R
- 1bfd93888d1fa947831dfcad0e1d47d0aecfde6d Fix the memory note's
  cross-reference link text to match its target heading
- 4e9f969a013c0e30221c5b6fa4cc470c715ec32e Reconcile the callback
  memory figures with the post-guard reference cases
- d22a0c01d9e8cd597674e9172683ceb0b55b738e Share the capi consumer
  compile between its two test files
- 85a36fefe7c0bc8016557f494304e954d9027acd Guard the running-mean
  recipe on the training channel's column count
- 55b4ac882084241b63a735f96ffd6f19c355f1dc State the no-R-API rule in
  the callback recipe's prose
- 683ed71e2a88ab35d2d1b5af56fd2142f9730dd2 Fold the callback bench
  scenarios into one loop and document the mode
- 8455f73ca26703d2b89a50f8c83336c44096d1e3 Retarget the
  compiled-consumer cites at the shared helper and the source
- 4cf56b69bca5f1855c856cdb6fbe524b876e4561 State the callback's
  pointer shape rather than quoting the request for it

Recipe 7 in
[dbarts-as-a-component.Rmd](../../vignettes/dbarts-as-a-component.Rmd) is
the running posterior mean of
[5. The Rcpp example](../design/per-draw-callbacks.md#5-the-rcpp-example),
Rcpp form, unevaluated - `Rcpp` is not a dbarts dependency, so the chunk
does not run as part of building the vignette. `R_MakeExternalPtrFn`
casts the callback's function pointer without the plain cast
`-Wpedantic` flags under `R CMD check`, and the caller-owned accumulator
goes in the external pointer's PROTECTED slot, the only thing keeping it
alive across the run. The prose states the lifetime rules - `REAL(acc)`
read exactly once before the fit starts, chains writing DISJOINT slices
so no lock is needed, a context built PER RUN because `drawIndex`
restarts at 0 on every call, and the return value an ABORT switch, not an
error channel - and, after review, the no-R-API rule in the Constraints'
own terms: no `Rf_allocVector`, no `PROTECT`, no `Rf_error`, and in C++
no construction or destruction of an `Rcpp` proxy type. The reviewer
compiled and ran the chunk out of band with `-Wall -Wpedantic`,
reproducing `yhat.train.mean` to 1.1e-16.

A plain-C mirror of the same callback,
[inst/tinytest/capi/consumer.c](../../inst/tinytest/capi/consumer.c)'s
`capi_meanDraw` (guarded, after review, on `numReportedLocations == 1`
beside the checks the recipe itself states), is driven by
[inst/tinytest/test-callback-example.R](../../inst/tinytest/test-callback-example.R)
at tolerance 1e-12 - the measured gap is 1.1e-16 absolute, so the
tolerance covers summation order (an incremental per-chain mean then
`rowMeans`, against `yhat.train.mean`'s own whole-array reduction) and
nothing else. [inst/common/capiConsumer.R](../../inst/common/capiConsumer.R)
is the compile/self-gate helper `test-capi.R` and this new file now
share (`compileCapiConsumer`, lifting what had been duplicated in both);
`exit_file` itself stays in each test file rather than moving into the
helper, since tinytest masks it only in the test file's own environment
and a call from the sourced helper would silently reach the namespace
version instead.

The manual gained one sentence:
[dbarts-package.Rd](../../man/dbarts-package.Rd)'s Memory section names
the callback route beside `keepTrainingFits`.
[memory-footprint.md](../../docs/design/memory-footprint.md)'s Reference
cases gained the two rows it owed - a heteroscedastic fit's variance
channel and a BCF fit's per-forest fits, neither reachable by
`keepTrainingFits` - both formula-derived, not measured on an actual fit.
Reconciling those against the post-ingestion-guard reference-case rows
mid-slice moved the design note's own figures (section 1 and
[6. Memory consequence](../design/per-draw-callbacks.md#6-memory-consequence))
and the manual sentence together: 3891.3 to 697.7 MB at n = 1e5, p = 20,
C = 4, S = 500 (engine 658.3 + R 33.0 + two `yhat.train` copies 3200.0
MB, against engine + R + 3.2 MB per-draw scratch + 3.2 MB accumulator)
and 10560.3 to 2576.3 MB at n = 1e6, p = 50, C = 1.

`benchmarks/R/bench-sampler.R` gained a `callback` mode (`Rscript
bench-sampler.R callback [record|compare ...]`, opt-in like the big
grid, own `sampler-callback.csv`), timing three variants - none, a
no-op C callback, and the running-mean recipe's own compiled copy, the
SAME source `test-callback-example.R` checks - at both reference
shapes. The consumer library is compiled once per invocation and the
mean context built once per shape, reused across that shape's timing
repetitions (only the elapsed time is read off; the correctness test
needs a fresh context per run, this scenario does not). NOT yet run:
the scenario lands with this slice but its numbers are maintainer-run
on a quiet machine, per Verification, and
[7. Threading interaction](../design/per-draw-callbacks.md#7-threading-interaction)'s
"measured at landing:" placeholder is still open.

Review fixes before landing: a docs/ path cite dropped from two shipped
comments (`inst/tinytest/capi/consumer.c` and
`inst/tinytest/test-callback-example.R` each named
`docs/design/per-draw-callbacks.md`, against the house rule that shipped
comments cite no docs/ path); the running-mean guard checks
`numReportedLocations` too, in both the vignette and its plain-C mirror,
since the recipe's reduction is written for a single reported location;
the no-R-API rule stated in the vignette's own prose, in the
Constraints' terms; a quoted VD request replaced by a statement of the
callback's actual pointer shape; and five doc cites retargeted from
`inst/tinytest/test-capi.R` to the new shared helper after its
extraction, two of them in `docs/plans/pure-c-header.md`.

Files beyond the plan's S4 list: `inst/common/capiConsumer.R` (new,
the shared helper, not itemized there), `docs/design/per-draw-callbacks.md`
(the memory reconciliation and the doc-cite retargeting) and
`docs/plans/pure-c-header.md` (two of the retargeted cites); the S4
header above is updated to list them. About 395 added lines before
review against a budget of about 150 (~70 vignette + ~60 tinytest C +
~20 manual and memory note); the reviewer judged the overage the
unbudgeted bench mode and the test file's boilerplate, the latter
removed at landing by the shared helper.

Gates (reviewer's run): `R CMD INSTALL --preclean` exit 0; tinytest
8351 TRUE, 0 FALSE, both `test-callback-example.R` and `test-capi.R`
ran (neither skipped); the equivalence trio against the f0236082
baselines 52 of 52, BCF 12 of 12, multinomial 11 of 11 "identical draws
(same RNG stream)", zero "max |z|", zero skipped; `R CMD build` and `R
CMD check --as-cran` from a tarball staged outside the tree OK, 0 NOTE;
`check-doc-freshness`, `check-rc-codoc`, `air format --check .` and
`lintr::lint_package()` all clean. Mutation: collapsing the per-chain
slice offset in the running-mean recipe fails
`test-callback-example.R`'s comparison against `yhat.train.mean`,
relative difference 0.5.

Remaining: the bench-sampler callback numbers owed to the design note's
section 7 ("measured at landing:" placeholder) are recorded separately,
by the maintainer, on a quiet machine. S5 (the sister-package record)
is open, per the plan's Steps.

## Landing note, S5 (2026-09-11)

LANDED, recorded in the commit that carries this note. One commit.

Records stan4bart as dec-B86's second named consumer in
[8. Scope and sequencing](../design/per-draw-callbacks.md#8-scope-and-sequencing):
verified against its current source, it calls neither
`dbarts_draw_callback` nor `dbarts_sampler_setDrawCallback`, so nothing
ports; the unchanged 1/0 version pair (dec-B111) admits its existing
compiled binary against S2's re-baked hash; its per-iteration,
one-draw-at-a-time embedding, each chain a separate single-chain sampler,
is the shape the entry serves, so it may adopt it post-release; and its
`keep_fits` argument, quoted in full, is the documented precedent
`keepFits` mirrors. Docs-only: no code, no engine gates.

Gate: `Rscript tools/check-doc-freshness.R` exit 0.

Remaining: S3 and S4 (the R surface wiring `callback` and the control
slot, the vignette example and manual) are open, per the plan's Steps.
