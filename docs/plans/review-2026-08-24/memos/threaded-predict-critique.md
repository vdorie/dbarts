# Threaded predict: adversarial critique

Target: threaded-predict-memo.md. Repo read-only at bartcore 0045507c; the ruling is taken as given, so everything
below is about HOW. Measurement build: `git archive HEAD` staged under scratchpad/tpcrit, `R CMD INSTALL --preclean -l
.../libs/tpcrit .`, R run under `R_LIBS=` that lib only.

## A. Archaeology

**1. CONFIRMED-OK - the three cited commits say what the memo says.** `git log -1` + `--stat` on each.
93f354a8b763000d, Vincent Dorie, 2025-03-03, "Add native predict multithreading", 18 files +249/-27. f3d019386a2c2619
2026-07-02, edcdf735855ad331 2026-07-03 (pre-rebase twin 7ff59ea4 under tag bartcore-pre-cran-rebase, identical
message). Neither of the last two mentions a thread count: "collateral, not a decision" holds. `git log -S` over
bartcore/main/the tag finds no other threaded-predict work.

**2. MAJOR - decision 5's premise is factually wrong: 0.9-31 SHIPPED.** Memo: "retracting a regression that never
shipped in a release is noise." `git show main:DESCRIPTION` -> Version 0.9-34, Date 2026-08-20; 93f354a8 bumps
DESCRIPTION under `Version: 0.9-31`; `git show main:src/dbarts/bartFit.cpp` still carries `BARTFit::predict(...,
size_t numThreads, ...)`. The feature went out in released 0.9-31 and still stands in released 0.9-34, so 1.0-0 as it
stands is a real regression against four releases. Two consequences the memo misses: (a) even AFTER this slice the
0.9-31 wording "parallelize across chains" is FALSE under the recommended (chain, draw) decomposition, so "once this
slice lands the statement is true again" is wrong; (b) the repo DOES fix NEWS in place pre-release - f04e8686 "Split
the released 0.9-34 notes back out of 1.0-0's NEWS section".

**3. MINOR - two records bearing on the serial choice, uncited.** No design doc records a reason (memo right). But
[[docs/plans/interface-review.md:200-205@f04e8686]] is the F10 item that renamed run's formal, with "(wiring threading itself is
wishlist)" and "Verified alongside D1 that both run() and predict()'s thread-count formals are fully inert, not merely
'serial'"; and [[docs/plans/interface-review.md:543@f04e8686]] puts "threaded prediction" on the explicit 2.0-WISHLIST. Cite both - the item is already docketed,
and the wishlist line must be struck when this lands.

## B. The partition

**4. CONFIRMED-OK - no shared mutable state a worker would touch. I tried hard to find one.** Read
predictColumns/predictPerForestColumns/predictVarianceColumns ([[sampler.hpp:565@f04e8686]], [[sampler.hpp:621@f04e8686]], [[sampler.hpp:674@f04e8686]]),
Chain::predictFromSavedSample / ...Multi / predictPerForest... and their live twins ([[chain.hpp:2797-2985@f04e8686]]),
addFlatPredictions ([[chain.hpp:2758@f04e8686]]), and the tree.hpp kernels addFlatPredictionsBelow / addFlatLinearPredictionsBelow /
partitionFlatIndices ([[chain.hpp:1723-1900@f04e8686]]).
- Adapters: DenseColumns/DenseColumnReader ([[tree.hpp:1700-1719@f04e8686]]), PredictorSourceColumns / SparseRawColumn /
  PredictorSourceColumnReader ([[data.hpp:346-442@f04e8686]]). All built in the constructor, `column(j) const` returns by value,
  `SparseRawColumn::at` is a pure popcount lookup with no search cursor and no memo. One shared `columns` object is
  safe.
- Accessors the replay reaches - [[model.hpp:984-988@f04e8686]], [[model.hpp:1347-1351@f04e8686]] (covariateColumns/Means/Sds/lengthscales/
  numParams/numCovariates) and fitScale/fitShift/sigmaScale ([[model.hpp:2718-2719@f04e8686]] + overrides) - are plain returns of stored
  members. No lazy init.
- `mutable` audit over src/bartcore: [[model.hpp:521@f04e8686]], [[model.hpp:1302-1303@f04e8686]], [[model.hpp:2093-2107@f04e8686]], [[model.hpp:2121-2127@f04e8686]] and [[tree.hpp:1604-1626@f04e8686]] are all
  sampling-path scratch (monotone neighbours, suffstat/GP kernel caches, pooled-mask and availability scratch); none
  is reachable from addFlatPredictions*, checked by reading the call graph, not by grep alone.
- SIMD dispatch ([[src/misc/simd.c:193-350@f04e8686]]): `misc_setVectorToConstant` / `misc_addVectorsInPlace` are file-scope function
  POINTERS assigned once by `misc_simd_init`, read-only after. No race.

**5. CONFIRMED-OK - bitwise identity at any thread count is structural.** The only accumulation is `fits[indices[k]]
+= leafValue` in addFlatPredictionsBelow, once per row per tree, with the tree loop t = 0..numTrees-1 in every entry
point: each (slab, row) pair owns its accumulator and sees the same addend order. The rescale, the multinomial softmax
and predictVariance's `out[i] *= treeFit[i]` are per-row maps or per-row products over the same j order. I looked
specifically for a running mean, an amplitude glue combine and a softmax across forests: the softmax is per (slab,
row) over K, inside one task. No cross-slab and no cross-row reduction exists.

**6. MAJOR - an escaping exception in a worker is std::terminate, and the design puts an allocation in every worker
iteration.** Memo: "zero engine surgery below predictColumns (per-task scratch is already function-local)". The
scratch is per-SLAB, not per-task: `std::vector<size_t> indices(numTest)` (plus blockOffsets, plus `raw(numTest*K)`
for multinomial) is constructed on EVERY call to predictFromSavedSample - numChains*numSavedDraws times, now
concurrently. Two costs: allocator traffic the partition makes contended, and a `std::bad_alloc` thrown inside a
std::thread body with no catch, which calls std::terminate and aborts the R session - and a large predict is where
bad_alloc is plausible. (Sampler::run has the same latent exposure: precedent, not novelty, but predict need not
inherit it.) Amendment: allocate scratch ONCE PER WORKER on the main thread before the spawn and pass it in; wrap each
worker body in `try { ... } catch (...) { failed[w] = true; }` and Rf_error after the join. Also a serial-path
speedup.

**7. MAJOR - "0 = defer" can resolve to zero workers and return uninitialized memory.** `Sampler::setNumThreads`
([[sampler.hpp:1011-1014@f04e8686]]) and `dbarts_sampler_setNumThreads` ([[C_interface.cpp:875-878@f04e8686]]) store the value with NO
validation; dbartsControl's check ([[R/A_class.R:349-350@f04e8686]]) guards only the R path. A C consumer that sets 0 and then
predicts with numThreads == 0 defers to 0; written as `numWorkers = min(numThreads, numSlabs)` that is zero workers,
nothing writes `out`, and predict hands R the contents of an uninitialized `Rf_allocVector(REALSXP, ...)`. Amendment:
clamp the resolved count to >= 1, pin it with a test, and say in the doc block what 0 means when the sampler's own
count is 0.

**8. MINOR - the flattenTree invariant needs the stronger statement.** The memo's reasoning is right
([[sampler.hpp:583-591@f04e8686]]), but the load-bearing half is that the partition is over SLABS: one slab per chain plus a slab
partition means no two workers touch one Chain. Under candidate (b) (row chunks) the invariant dies. Amendment: assert
it in the capacity == 0 arm and say in the comment that a row-axis partition must privatize the flatten buffers.

**9. MINOR - the cutoff predicate is wrong for three entry points.** `numChains * numDraws * numTrees * numTest` uses
forests_[0].numTrees. A multinomial replay traverses sum_f numTrees_f ([[chain.hpp:2871-2884@f04e8686]]), the per-forest replay
likewise, and the heteroscedastic second fan-out traverses numVarianceTrees. State the predicate as "the traversals
THIS entry performs" and compute it per entry.

**10. CONFIRMED-OK - no RNG, no R API, no Rf_error inside the threaded region.** Read `dbarts_sampler_predict`
([[C_interface.cpp:773-802@f04e8686]]) and `predictFromSource` ([[R_interface_bartcore.cpp:5654-5758@f04e8686]]). Every refusal and every
R_alloc is strictly BEFORE the engine call: refuseUndefinedTestFits ([[R_interface_bartcore.cpp:5781@f04e8686]] / [[C_interface.cpp:780@f04e8686]]),
refuseEmptyTreeStore, translateSource ([[C_interface.cpp:185-241@f04e8686]] - the R_alloc the header's main-thread warning is
about), validateTestSource, the offset-shape refusals ([[src/R_interface_bartcore.cpp:5631-5645@f04e8686]]), and the result allocation. No ext_rng/unif_rand
and no Rf_error/ext_printf/R_alloc under chain.hpp's replay block. [[dbarts.h:45-47@f04e8686]]'s main-R-thread-only contract is
unaffected and must NOT be relaxed - translateSource still runs R_alloc.

**11. MINOR - the interrupt argument is stronger than the memo makes it.** A Ctrl-C while the main thread is blocked
in `worker.join()` cannot longjmp out: R's unix SIGINT handler only sets R_interrupts_pending, and the jump happens in
R_CheckUserInterrupt, which predict never calls (the bridge's poll, `bartcore_checkInterrupt`,
[[R_interface_bartcore.cpp:4224-4231@f04e8686]], is R_ToplevelExec-wrapped and used only by run). So deferring interrupt polling
(decision 7) is SAFE, not merely unchanged, and the SIGINT mask mirror is defence-in-depth rather than necessity. Say
so.

## C. Fan-out and thread count

**12. CONFIRMED-OK - mirror run's fan-out; do not reuse testFitPool_.** [[sampler.hpp:349-420@f04e8686]] and the grow-from-root
twin [[sampler.hpp:1072-1099@f04e8686]]: std::thread, numWorkers = min(numThreads, numChains), pthread_sigmask(SIG_BLOCK, SIGINT) around the
spawn under `#ifndef _WIN32`, mask restored immediately, join at the end. Portability is configure-provided
([[configure.ac:53@f04e8686]] AX_PTHREAD; [[src/Makevars.in:18-34@f04e8686]] threads the flags through), and Windows needs no mask because R's
console Ctrl-C already lands on the main thread. testFitPool_ ([[chain.hpp:5353-5354@f04e8686]]) is per-Chain, budget
numThreads/numChains, 65536-ROW cutoff; a slab partition is cross-chain, so the rejection is right. Spawn cost (tens
of microseconds x a few threads) against a 25 ms cutoff is under 1%.

**13. MINOR - two arms the memo does not spell out.** run's `numWorkers <= 1` branch runs inline on the caller's
thread; predict must too, rather than spawn one worker and join. And numWorkers must be capped at the SLAB count, or a
4-thread predict of a 1-chain no-keepTrees fit spawns three workers with nothing to do.

**14. MAJOR - `_R_CHECK_LIMIT_CORES_` does not constrain native threads, and the memo leans on it.** Memo: "a default
that probes cores is exactly what `_R_CHECK_LIMIT_CORES_` exists to catch", and "The repo has no
`_R_CHECK_LIMIT_CORES_` anywhere ... new discipline, not a regression." The variable is read only by `parallel`'s
`.check_ncores` - makeCluster/mclapply, not std::thread. The repo's own record proves it:
[[docs/plans/release-candidate-review.md:2529-2531@f04e8686]] records the single trip as "xbart's auto n.threads trips CRAN's core
limit" (xbart uses `parallel::makeCluster`, [[man/xbart.Rd:61@f04e8686]]) "fixed with an explicit n.threads = 1L and the
as-cran-for-thread-spawning-tests lesson recorded in the runbook", and [[docs/plans/release-candidate-review.md:2210-2211@f04e8686]] records the P8 battery running tinytest
"with and without _R_CHECK_LIMIT_CORES_". The discipline exists; the mechanism does not protect dbarts's own threads.
Any "skip the >2 arms when it is set" clause must read the env var by hand - nothing errors for you.

**15. MAJOR - "Every existing inst/tinytest fit passes 1L or 2L" is false.** Parsed every `bart2(`/`rbart_vi(` call in
inst/tinytest/*.R with a paren matcher: 117 declare no n.threads, and 50 of those also have n.chains > 1 (or take the
default 4), so their control@n.threads is `min(guessNumCores(), n.chains)` - up to 4 on any CI box. Worst files:
test-multinomial-surface.R (17), test-hurdle-surface.R (7), test-argument-surface.R (3). Those fits already spawn up
to 4 workers in `run`, so the CRAN posture is pre-existing - but under decision 2's default every `predict()` on them
becomes 4-threaded too, and test-multinomial-surface.R is predict-heavy. The premise for "never a CRAN surprise"
fails.

## D. Header, bridge, R surface

**16. MINOR - the hash re-bake is TWO literals.** [[src/C_interface.cpp:461@f04e8686]] `static_assert(dbarts_apiSignatureToken ==
0x85bd1ef04beb3848ULL, ...)` fires on any SIGNATURE change; [[dbarts.h:142@f04e8686]] DBARTS_C_API_HASH on signature-or-layout.
Adding a parameter moves both. Both fail loudly with instructions, so it is self-correcting, but the flat-API cost
line names only one. stan4bart pins the hash by hard equality (its `src/init.cpp` line 972), so the lockstep break is real and
intended.

**17. MINOR - no default arguments on the facade virtuals.** [[facade.hpp:258@f04e8686]]/267/272 with overrides at [[facade.hpp:545@f04e8686]]/[[facade.hpp:549@f04e8686]]/[[facade.hpp:554@f04e8686]]
(verified). A default value on a virtual binds from the STATIC type: the base's for `SamplerBase&` calls, the
derived's for derived calls - silent divergence, and every real call goes through `SamplerBase&`. Amendment: no
default arg; explicit at every call site. Also give the non-virtual dense convenience spellings ([[facade.hpp:276-286@f04e8686]],
re-exposed by `using SamplerBase::predict;` at [[facade.hpp:542-543@f04e8686]]) the parameter, or the dbarts.h-shaped overloads quietly stay
serial while their view-taking twins thread.

**18. MAJOR - tests/cpp/test_facade.cpp is an unlisted edit site, and the cheapest answer to the memo's own named
failure mode.** It enumerates every virtual (FacadeVirtual list [[facade.hpp:40-41@f04e8686]]), spies them (`SPY_VOID(predict, ...)`,
`SPY_VOID(predictPerForest, ...)`, `SPY_VOID(predictVariance, ...)` at [[facade.hpp:140-146@f04e8686]]) and conformance-checks all three
predict virtuals at [[facade.hpp:752-785@f04e8686]]. A signature change forces edits in all three places; the cost table has no line for it.
Amendment: extend the spy to RECORD the numThreads it was handed and assert the forwarder passes it through. That
deterministically disproves "the wiring is a no-op" for ~10 lines instead of a timing assertion.

**19. MAJOR - the DEF_FUNC arity bump breaks four tinytest files the memo does not list.** Grepped every `.Call` of
the symbol. Besides [[R/dbarts.R:1140@f04e8686]] and [[R/bartcore.R:1354@f04e8686]], four test files call the entry point DIRECTLY with three
arguments: [[test-predict-sparse.R:163@f04e8686]] (inside `expect_identical` - a hard error), [[test-multinomial-test-offset.R:405@f04e8686]]
and [[test-multinomial-test-offset.R:414@f04e8686]], and [[test-multinomial-category-offset.R:430@f04e8686]] (inside `expect_error` with message patterns - they fail on the
arity message instead of the target refusal). Add [[R/bartcore.R:1490@0045507c]] and [[R/dbarts.R:1153@f04e8686]] for the per-forest arity.
Loud, but unbudgeted.

**20. MAJOR - bartCause is not "ZERO lines" without a placement rule.** bartCause dbarts-1.0 @7ae6e83
(its `R/generics.R` lines 141-230): predict.bartcFit forwards `...` into FOUR different predicts - dbarts predict.bart,
predict.rbart, stan4bart's, and stats predict.glm/lm/merMod (`predict(object$fit.trt, x.new.g, type = "response",
...)`) - and it calls predict.rbart POSITIONALLY: `predict(object$fit.trt, x.new, group.by, combineChains = FALSE,
...)`. So the new formal MUST be appended AFTER every existing positional formal on each of the six generics;
inserting `n.threads` before `group.by` in predict.rbart ([[R/generics.R:1647-1657@f04e8686]]) would pass a factor as a thread
count. The migration is still zero LINES, but only under that constraint, and the memo says only "add the formal to
each". Each generic also needs its own default expression: predict.bartHurdle ([[R/generics.R:1614@f04e8686]]) has no `object$fit` at all (it
is `object$occupancy$fit`). stan4bart is exactly as described: one C site, stan4bart's `src/init.cpp` line 342 in predictBART (line 293),
hash equality line 972, `-DDBARTS_USE_STUBS` (stan4bart's `src/Makevars.in` line 1).

**21. MAJOR - predict.bart's validation sits BELOW two early returns, and predictBlend is missing from the site
list.** [[R/generics.R:207-300@f04e8686]]: `type == "forest"` returns through `predictForest` at [[R/generics.R:660@f04e8686]] and an amplitude-coupled fit
through `predictBlend` at [[R/generics.R:672@f04e8686]], both BEFORE `n.threads <- as.integer(n.threads)[1L]` at [[R/generics.R:294@f04e8686]]. Putting the
house-pattern validation at [[R/generics.R:294@f04e8686]] leaves the amplitude family's n.threads unvalidated. And predictBlend ([[R/generics.R:727-738@f04e8686]]) is
how a BCF fit's PLAIN predict reaches the replay - it calls `predictForest(object, newdata, NULL, TRUE, NULL)` at
a line whose location could not be placed at this sha - `predictBlend` does not exist in
R/generics.R at f04e8686 (unresolved: [[R/generics.R:738@f04e8686]]). The memo lists predictForest but not predictBlend. Amendment: move the validation above both early returns; give
predictBlend the formal and forward it.

## E. Measurement

**22. MAJOR - could not reproduce, and the protocol cannot produce a valid share.** Same hardware class (Apple M1 Max,
8P+2E). Fit `bart2(n = 2000, p = 10, n.trees = 200, n.samples = 250, n.chains = 4, keepTrees, n.threads = 1)` = 1000
slabs; timed `sampler$predict(xt)` and `predict(fit, xt)` at nTest = 2000, 7 paired reps, gc between. Engine reps 5.41
7.28 6.47 6.36 6.48 6.97 2.33 s; full reps 6.70 6.96 5.50 6.43 6.36 6.44 1.73 s -> min(engine) = 2.330, min(full) =
1.734, "share" 134%. `full` STRICTLY CONTAINS `engine`, so a share above 100% is arithmetically impossible: it is what
dividing the minima of two INDEPENDENT rep series taken in different load regimes produces. That is the memo's
protocol ("min of 5 reps, gc between"), and it is the flaw, not just my noise - this host was at loadavg 125 and the
same call varied 3.1x across reps. CLAUDE.local.md already requires a quiet machine for bench-sampler; the memo's
table states no such condition. Best-case ns/traversal here was ~4.3 against the memo's 2.40-2.74, so the figure is
unconfirmed either way. Amendment: compute the share from PAIRED reps (per-rep ratio, then median) or time the R
reshape directly and subtract; state the box was quiet; re-take the table before quoting it.

**23. MAJOR - the Amdahl ceiling is over the wrong denominator.** The measured share is (whole
`.Call`)/(predict.bart), but the partition parallelizes only `predictColumns`. Serial work INSIDE the `.Call` that it
never touches: the flat-offset add over EVERY slab ([[R_interface_bartcore.cpp:5736-5739@f04e8686]] - one pass over the entire
output, 76 MB at the memo's nTest = 10000); `Rf_duplicate(resultExpr)` on the heteroscedastic path ([[R_interface_bartcore.cpp:5745@f04e8686]], a full copy
of that array) before a SECOND fan-out for predictVariance; the result allocation; translateSource; R-side
validateXTest/as.matrix. The measured configuration had no offset, no variance forest and a dense test set - the most
favourable case - and its ceiling is reported as general. Report the worst configuration too, or say plainly that
99.5% is the best case.

**24. MINOR - the bandwidth worry is overstated, in the memo's favour.** Per slab the engine touches only slab-sized
output (16-80 KB, L1/L2 resident) and re-touches it numTrees times; the streamed state is the saved-tree store, a few
tens of MB per FULL pass - single-digit MB/s at the measured runtime, nowhere near saturation. The shape is
latency/branch bound (per-tree `std::vector<FlatNode>` pointer chasing, a data-dependent partition in
partitionFlatIndices), which threads well. "76 MB of output ... the shape of a bandwidth-bound problem" conflates a
one-time output write with a streamed working set. Keep the scaling curve as the gate; drop the presumption of a wall.

## F. The seven open decisions

**25. Decision 1 (0-as-defer). AGREE**, with finding 7's clamp and a doc sentence for the sampler-count-is-0 case.

**26. Decision 2 (control@n.threads as predict's R default). REFUTED.** Both supporting facts are false (14, 15), so
the argument has nothing under it. Substantive case for 1L: predict is routinely called on small newdata, and
partialDependence calls `sampler$predict(x.test)` per grid level with NO n.threads ([[R/partialDependence.R:340@f04e8686]], [[R/partialDependence.R:357@f04e8686]],
[[R/partialDependence.R:359@f04e8686]]), taking the R5 default silently once per level; the fit-time budget was sized for CHAIN parallelism
(the memo says so in 4.3); the house rule is "safe over fast in R". Recommend `n.threads = 1L` for 1.0-0, revisited in
S2 against a real curve. If control@n.threads is kept, the memo must state that examples and tests are the exposure
and that no check tooling will catch it.

**27. Decision 3 (run deferred). AGREE.** Add: the Rd item cannot keep joint wording - after this slice "run and
predict both execute serially" ([[man/dbartsSampler-class.Rd:152-157@f04e8686]]) is half false, so splitting it is forced, not
optional.

**28. Decision 4 (predictPerForest in S1). AGREE, and the case is stronger than argued.** Because predictBlend routes
an amplitude fit's PLAIN predict through predictForest (finding 21), leaving predictPerForest serial leaves every BCF
fit's ORDINARY predict serial, not just `type = "forest"`.

**29. Decision 5 (leave the 0.9-31 NEWS entry). REFUTED** - finding 2. The entry describes a released feature, its
stated axis stays wrong after this slice, and f04e8686 shows the repo edits NEWS in place pre-release.

**30. Decision 6 (total-traversal cutoff). AGREE** on the predicate class; amend per finding 9 and make it a named
constexpr carrying the ns/traversal derivation, in the style of `testFitParallelCutoff` ([[chain.hpp:5354@f04e8686]]).

**31. Decision 7 (interrupt polling deferred). AGREE**, with finding 11's stronger reason.

## G. Budget, slicing, tests

**32. MINOR - the line budget is low by roughly 2x.** `git show --stat 63df524e` (the comparable per-forest predict
slice): +651/-18 across 15 files, with NO header/ABI change, no R generic fan-out and no bench - and it still needed
202 test lines plus 112 in tests/cpp. This slice adds all of those. Expect 700-900 diff lines, not ~425 dense. Not a
blocker; the S1/S2 split and "do not split the header out of S1" are both right.

**33. MAJOR - the R timing assertion is not a viable observable.** `system.time(...)[["user.self"]]/[["elapsed"]] >
1.2` is precisely the environment-fragile assertion the repo's own review waves have been removing (the P8
anti-vacuity pass; the xbart pin at [[release-candidate-review.md:2529@63df524e]]). On a loaded runner elapsed inflates and the
ratio collapses; on a single-core runner it can never pass. Amendment: make the observable DETERMINISTIC - have
predictColumns return (or record) the resolved worker count and the slab->worker map, and assert in tests/cpp that the
count is min(requested, slabs) and the map covers every slab exactly once at partition counts 1, 2, 3,
7. With finding 18's spy assertion, "the argument reached the engine" and "the engine used it" are both
proved without a clock. Keep the R side to `identical()` across thread counts; if an R-side liveness probe is still
wanted, gate it on NOT_CRAN and make it a message, not an assertion.

**34. MINOR - the sanitizer recommendation rests on a false premise.** [[sanitizers.yaml:77@63df524e]] runs
`tinytest::test_package("dbarts")` - the WHOLE suite - so the new file runs under ASAN automatically; there is nothing
to "run instead". ASAN+UBSAN detect no data races at all, so "the ASAN leg covers it" is not available as a claim.
Cheap real coverage: tests/cpp is a standalone binary with its own Makefile and no R, so a `-fsanitize=thread` arm
over the predict partition test costs a make flag, not an r-hub container. Recommend that in place of "no TSAN". (The
memo is right that the equivalence trio must be unchanged - predict touches no sampling code.)

**35. MINOR - citation drift in one region.** Against 0045507c the memo's R_interface_bartcore.cpp numbers run ~13 low
around the 2800/5600-5800 band: predictFromSource is [[src/R_interface_bartcore.cpp:5654@0045507c]] (memo [[src/R_interface_bartcore.cpp:5641@0045507c]]), refuseUndefinedTestFits defined [[src/R_interface_bartcore.cpp:2871@0045507c]] and
called [[src/R_interface_bartcore.cpp:5781@0045507c]] (memo [[src/R_interface_bartcore.cpp:2858@0045507c]], [[src/R_interface_bartcore.cpp:5768@0045507c]]). Everything else spot-checked is exact - facade.hpp, sampler.hpp, chain.hpp,
dbarts.h, man/bart.Rd, man/dbartsSampler-class.Rd, R/dbarts.R and R/generics.R all match.

## Verdict

**EXECUTE-WITH-AMENDMENTS.** The core technical claim - that a (chain, draw) slab partition is bitwise-identical by
construction and touches no shared mutable state - survived a deliberate attempt to break it, and choosing the slab
axis over the classic chain axis is right. What does not survive is the surrounding evidence: two CRAN-safety facts
are false, the measurement protocol cannot produce the share it reports, the NEWS decision rests on a wrong premise
about what shipped, and five edit sites are unbudgeted.

Amendments to fold in before implementation:

1. Clamp the resolved thread count to >= 1 after the 0-defer resolution; test setNumThreads(0).
2. Hoist per-slab scratch to per-worker, allocated before the spawn; wrap worker bodies in try/catch and re-raise
   after the join.
3. Move predict.bart's n.threads validation above the `type == "forest"` and forestFits early returns; give
   predictBlend the formal and forward it.
4. Append the new formal AFTER every existing positional formal on all six generics (critically after group.by on
   predict.rbart); give each its own default expression.
5. Replace the timing assertion with a deterministic worker-count / partition-coverage observable in tests/cpp, and
   extend test_facade.cpp's spy to record and assert the forwarded numThreads.
6. Add test_facade.cpp, the four direct `.Call` tinytest sites, [[R/bartcore.R:1490@0045507c]] and [[src/C_interface.cpp:461@0045507c]]'s second
   hash literal to the edit-site list and the budget.
7. Strike the `_R_CHECK_LIMIT_CORES_` safety claim and the "nothing exceeds 2 today" claim; re-argue decision 2 on
   merits (recommend 1L) or state the unguarded CRAN exposure explicitly.
8. Re-take the timing table on a quiet box with PAIRED reps, and report the Amdahl ceiling for the offset +
   heteroscedastic configuration alongside the best case.
9. Correct decision 5's premise (0.9-31 shipped; main is 0.9-34) and make the 1.0-0 NEWS entry state the (chain, draw)
   decomposition, since "across chains" stays false otherwise.
10. Per-entry-point cutoff traversal counts; no default arguments on the facade virtuals; give the dense convenience
   spellings the parameter; assert the one-slab-per-chain invariant in the capacity == 0 arm.
