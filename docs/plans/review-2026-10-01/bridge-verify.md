# Review 3 - verification of bridge.md

Tree .claude/worktrees/review3 at 01dee4b4, library r3-lib; 0.9-34 comparisons against the CRAN build in
r3-rfit-cranlib. Own probes: scratchpad r3-verify-bridge-*.R and a C shim r3-verify-bridge-shim/w.c
(DBARTS_USE_STUBS against the installed inst/include, identical to the tree's header). Every entry below was
re-derived, not re-run from bridge.md's scripts.

Summary: 9 findings, 7 CONFIRMED, 2 QUALIFIED (07 is documentation, 05 not a regression). None of the fixes
moves draws for a valid caller.

| id | verdict | severity | regression vs 0.9-34 |
|---|---|---|---|
| 01 | CONFIRMED | BLOCKER by the rubric (silent garbage), practically a MAJOR validation regression | yes |
| 02 | CONFIRMED | BLOCKER (heap overflow under the header's own contract); no shipped consumer reaches it | no (multinomial is new) |
| 03 | CONFIRMED | MAJOR | no (installTrees is new) |
| 04 | CONFIRMED (allocation arm probed; thread-spawn arm by code only) | MAJOR | partly (create aborted in 0.9-34 too) |
| 05 | QUALIFIED | MAJOR/MINOR boundary; pre-existing | no (0.9-34 identical) |
| 06 | CONFIRMED | MINOR | no |
| 07 | QUALIFIED | MINOR, a header-comment fix | no |
| 08 | CONFIRMED, slightly worse than stated | MINOR | no |
| 09 | CONFIRMED | MINOR | partly (run's negative-length message) |

---

## bridge-02 - CONFIRMED, BLOCKER, not a regression

Reachability. The header puts multinomial handles inside the contract in four places. THE HANDLE says
every handle is an R dbartsSampler's pointer, and dbarts() builds multinomial samplers. dbarts_sampler_family
is "total over every sampler any construction path can build" and returns DBARTS_FAMILY_MULTINOMIAL.
dbarts_sampler_getLatents documents its multinomial answer (0). The dbarts_results comment says logLikelihood
is NaN-filled for "the multinomial softmax", which describes a multinomial run through this entry, and
dbarts_draw documents L = K. The stride comment in dbarts_sampler_run ("1 for every dbarts.h-created sampler,
since the flat C API builds no multi-location model") was true when 45440d1f added it (2026-07-14). It went
stale on 2026-09-08 when d0403979 (pure-c-header) removed flat creation and made every handle an R object.
docs/design/feature-matrix.md [f4] and the TODO item multinomial-doors (Door 2) record flat multinomial
support as an unbuilt, consumer-gated door. So the intent is "not supported", but nothing refuses it.

Own probe (r3-verify-bridge-02.R plus shim w.c). n = 30, K = 3, 4 samples, 2 chains, nTest = 10. Buffers
are 6x the documented size and NaN-filled; the table shows the last written index against the documented
size:

```
multinomial kt=FALSE fam=8 | sigma 8/8 train 720/240 test 240/80 varcount 16/16 | predict st=1 60/20
multinomial kt=TRUE  fam=8 | sigma 8/8 train 720/240 test 240/80 varcount 16/16 | predict st=1 240/80
gaussian    kt=TRUE  fam=1 | sigma 8/8 train 240/240 test 80/80 varcount 16/16 | predict st=1 80/80
```

The train layout is n x K x S x C on the probability scale: row sums over K are exactly 1, and the R run
returns 30 x 3 x 1 x 2. A first version of the shim, with buffers 2x the documented size, overran them, and
the process segfaulted later in unrelated R code (exit 139). That is real heap corruption, not just a
canary hit. varcount is safe because the engine clamps it to the single documented slab.

A second, smaller defect: with a non-null offsetTest, flat predict adds the offset slab by slab at stride
numTestObservations (C_interface.cpp, dbarts_sampler_predict). On a multinomial sampler that adds an additive
offset to probability channels, at the wrong positions.

Why setting the stride to 1 is NOT a fix. Chain::storeSample writes n x combiner->numReportedLocations() per
draw whatever the caller declares. A stride of n would make the draws overlap, and the last draw would still
overflow. The current line is what keeps the R route correct.

Fix (recommended, no header change in any ABI sense):
- dbarts_sampler_run: when shape.numReportedLocations > 1 and results carries a non-null train or test
  (present by size), Rf_error before any sweep, naming the R $run route and the callback. Allow the run when
  both are null. A multinomial host can then still drive the sampler flat and read the K channels through
  dbarts_draw, which already documents L. sigma, varcount and logLikelihood are safe as sized.
- dbarts_sampler_predict: return 0 (capability status) when numReportedLocations > 1. The header's general
  rule ("0 means the SAMPLER cannot do this at all ... a capability the model does not carry") already
  covers this; only the entry's own sentence ("refused on any sampler whose blend is undefined") needs a
  second clause.
- Fix the stale stride comment, and amend the dbarts_results logLikelihood sentence and the
  train/test/predict layout lines to say a multi-location sampler refuses those buffers.
- These are comment-only edits to dbarts.h. DBARTS_C_API_HASH folds signatures, enumerators and layouts, not
  comments, so no hash re-bake, no version move and no consumer rebuild. It is still a shipped-contract
  wording change, so it goes to VD.

Alternative (an ABI event, maintainer decision): document the K-wide layout. That needs a flat accessor for
K, since no current entry reports it before a run (dbarts_sampler_numTrees refuses past the last forest),
and a new entry re-bakes the hash and the signature token. That is Door 2 of TODO multinomial-doors, which
is post-RC and gated on a consumer asking for it. Not recommended now.

Tests: add a multinomial arm to inst/tinytest/test-capi.R and capi/consumer.c:
- flat run with a train buffer raises;
- flat run with null train and test plus a callback sees numReportedLocations == K;
- flat predict returns 0;
- the gaussian canary stays clean.

Gates missed it because test-capi drives only gaussian and binary handles.

## bridge-01 - CONFIRMED, BLOCKER by the rubric, regression

Code: bartcore_run and bartcore_runWithCallback compute static_cast<size_t>(Rf_asInteger(...)) with no sign
check. bartcoreSamplerRun (R/bartcore.R) checks only for NA. 0.9-34 used rc_getInt with RC_GEQ 0 plus a
refusal of 0 + 0 (main:src/R_interface_sampler.cpp).

Own probe (r3-verify-bridge-01.R, keepTrees, one chain):
- run(-1, 4): sigma 0 0 0 0.
- run(-2, 5): sigma 4.881e-313 6.154e-313 ... For both, predict returns identical columns (the never-written
  store).
- run(0, 0): returns NULL.
- run(3, -1): "negative length vectors are not allowed".
- run(-5, 0): about 2^64 sweeps. It hit my 10 s setTimeLimit ("sampler run interrupted"); without the limit
  it does not finish.

0.9-34 on the same calls:
- run(-2, 5): "number of burn-in steps must be greater than or equal to 0".
- run(0, 0): "either number of burn-in or samples must be positive".

Fix: refuse NA-free negative counts in bartcore_run and runWithCallback, and refuse 0 + 0, using the 0.9-34
messages. Mirror the check in bartcoreSamplerRun. Draws do not move. Add a tinytest in test-sampler-run (or
similar) calling $run(-1L, 1L), $run(1L, -1L) and $run(0L, 0L) with expect_error.

## bridge-03 - CONFIRMED, MAJOR, not a regression

Code: Sampler::installForests (src/bartcore/sampler.hpp) runs every check before the commit loop except
Chain::installForest's own rebuild. installForest calls restoreScale and then rebuilds forest by forest, so
a rebuildLiveForest or rebuildLiveForestRemapped false on chain c returns shapeMismatch after chains
0..c-1 are installed. The variance loop runs after every mean install has committed, so a varianceMismatch
there leaves the same partial state. This contradicts the doc comment "On any mismatch nothing is touched".

Own probe (r3-verify-bridge-03.R), a 2-chain recipient with every pair of samples tried: all 12 refused
with "not shape-compatible (number of ...", and every one showed "chain1 changed TRUE, chain2 changed
FALSE".

A control in r3-verify-bridge-03b.R tried an ordinary cross-grid donor at 1, 2 and 3 chains over 6 seeds.
Every install succeeded, so the remap path itself works, and the defect is only the non-atomic refusal.

Side note, not chased: the donor that fails, from setPredictor plus setData after recording, fails on its
live trees as well. Why its trees no longer rebuild may deserve its own look.

Fix: build each chain's forests into scratch, or snapshot (fitMin/fitMax, the trees, sigma, k, glue) and
restore on failure, then commit. Give the rebuild failure its own WarmStartResult and message. Draws do not
move: installs consume no RNG and a successful install is unchanged. Add a tinytest that asserts the
recipient's predictions are identical after a refused install, not only the message.

## bridge-04 - CONFIRMED (allocation), MAJOR, partly pre-existing

Code: bartcore_growFromRoot, bartcore_setControl (sampler.setTreeStorage), createHolder
(createSampler/createAmplitudeSampler) and the flat dbarts_sampler_setTreeStorage have no captureExceptions.
That falsifies the comment on captureExceptions in R_interface_bartcore_common.hpp.

Own probe (r3-verify-bridge-04.R): flat setTreeStorage(1, 2^60) gives "libc++abi: terminating due to
uncaught exception of type std::length_error: vector", exit 134. The setControl probe instead hit the OS
memory killer (exit 137) because macOS overcommits, so I relied on the code reading there.

Thread spawn: workers.emplace_back at [src/bartcore/sampler.hpp:558](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/src/bartcore/sampler.hpp#L558), [src/bartcore/sampler.hpp:837](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/src/bartcore/sampler.hpp#L837) and [src/bartcore/sampler.hpp:1506](https://github.com/vdorie/dbarts/blob/01dee4b4f1a21565c91eceea391e666851662a0d/src/bartcore/sampler.hpp#L1506) is outside any guard. That is accepted
on code grounds only; not probed, given the thread cap.

Fix as proposed: wrap those entries in callConvertingExceptions or captureExceptions, and have each fan-out
join the threads it started before rethrowing. Draws do not move. The flat setTreeStorage size probe can
live in test-capi.

## bridge-05 - QUALIFIED, MAJOR/MINOR, not a regression

Own probe (r3-verify-bridge-0569.R):
- setResponse with an Inf entry is accepted; sigma is NaN from then on, including after the response is
  restored and 50 more sweeps.
- setOffset with an Inf entry gives the same result under both updateScale = TRUE and FALSE, and sigma stays
  NaN after setOffset(NULL).
- 0.9-34 also accepts it and stays NaN after the restore, so this is not a regression.

The flat side contradicts the header: the header requires setResponse values to "lie in the family's
support", but validateResponseSupport checks nothing for gaussian.

Fix: refuse non-finite values in both conduits, R and bridge, with the creation message (the bridge check
covers the flat route). This is a new refusal, not a draw move. Add a tinytest with expect_error on Inf for
both conduits.

## bridge-06 - CONFIRMED, MINOR

BCF setOffset(updateScale = NA): refused, and data@offset afterwards is 5 5 5.
bartcoreSamplerSetOffset assigns the mirror before the .Call. setResponse was already reordered to
validate first, with a comment explaining why, so setOffset is the inconsistent one.

Multinomial setOffset is refused by R before the mirror is touched, so it is safe.

Gaussian setOffset(updateScale = 1) is accepted.

Fix as proposed: check updateScale as a single TRUE or FALSE and assign the mirror after the .Call succeeds.

## bridge-07 - QUALIFIED, MINOR, documentation

translateSource checks the CSC triple's presence and the columnSources range, but not pointer monotonicity
or row-index range. That is consistent with the header's "validation is deliberately partial", and TODO
python-bindings records that CSC structure validation lives in the R bridge. "A source whose declared shape
disagrees" means numRows and numColumns, which are checked. So the claim that the header promises this check
overreads it.

Fix: one header-comment sentence saying the CSC triple's structure is unchecked. That is comment-only, so no
hash change. Alternatively add an O(nnz) check in translateSource.

## bridge-08 - CONFIRMED, MINOR, slightly worse than stated

Own probe (r3-verify-bridge-08.R): after a flat setTreeStorage(1, 4), an R run and storeState, the live
predict is 20 x 4. After saveRDS and readRDS, both predict AND run fail with "state is not consistent with
this sampler", so the reloaded object is unusable, not just its predict. The failure is loud, so it stays
MINOR.

Fix: have setState adopt the state's store capacity (bridge-side, no header change), or document in THE
HANDLE that flat settings do not survive a re-creation.

## bridge-09 - CONFIRMED, MINOR

Reproduced all three messages exactly:
- setPredictor(x, NA_integer_) and getTrees(chainNums = NA_integer_) give "missing value where TRUE/FALSE
  needed".
- run(5L, -1L) gives "negative length vectors are not allowed".

Fix as proposed. 01's fix covers the run message.
