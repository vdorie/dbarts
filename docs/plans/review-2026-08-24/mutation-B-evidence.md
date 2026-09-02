# Mutation testing leg B: long-form evidence (companion to mutation-B-findings.md)

## Ladder audit

Probe: a 7th ResponseFamily enumerator in a staged copy; tests/cpp rebuilt and
R_interface_bartcore.cpp compiled -fsyntax-only, -Wall -Wextra.  Exactly THREE
sites warned (-Wswitch): [[chain.hpp:582@658869ac]] (single-forest response factory),
R_interface:2268 (refusedAmplitudeFamilyReason), unresolved: [[chain.hpp:6181@658869ac]] (drawAugmentationLaws).
Four `default:` ladders stayed silent, plus one open-coded chain no -Wswitch
can reach.  (Only two of the four are in src/bartcore; the memo's list spans
both files.)
  [[chain.hpp:756@658869ac]]  AmplitudeSpec ctor        aft/ord/nbinom/new  -> GaussianResponse
  [[chain.hpp:5026@658869ac]] latentScaleAnchor         gauss/aft/ord/nb/new -> scaledResponseSd()
  R_interface:2291 defaultNodeScale        gaussian/aft/new    -> 0.5 (gaussian's)
  R_interface:2812 validateResponseSupport gaussian/aft/new    -> NO validation
  R_interface:6254 computeWorkingResponse  if/else, no switch  -> latent[i]
By call graph, [[chain.hpp:756@658869ac]] and [[chain.hpp:5026@658869ac]] take spec.family from the K-forest route, which
R_interface:2268 refuses for aft/ordinal/nbinom BEFORE either runs, so those
arms are gaussian-only today.  [[chain.hpp:2812@658869ac]] is the worst door - its own comment
enumerates the harm it prevents (a negative nbinom count underflowing into a
~1.8e19 allocation, "an uncatchable crash, not an error").  The de-facto 7th
family already exists: [[chain.hpp:862@658869ac]] sets `family_ = logistic` for multinomial
commenting "family() is not read on this path", yet SamplerShape::family
([[facade.hpp:432@658869ac]]) reports logistic and the bridge branches on it at unresolved: [[facade.hpp:2729@658869ac]],
:4611, unresolved: [[facade.hpp:4662@658869ac]], unresolved: [[facade.hpp:4702@658869ac]], unresolved: [[facade.hpp:4877@658869ac]], unresolved: [[facade.hpp:4949@658869ac]] - each unreachable only because a SEPARATE
multi-forest refusal stands in front, a coupling nothing asserts.
VERDICT the doors are shut by refusals in another file, not by the type system,
and one family already travels under another's name.

## Assertion audit

R compiles packages with -DNDEBUG (Makeconf:174), so every runtime assert() in
src/bartcore is DEAD in the shipped library; tests/cpp never defines it, so the
suite runs with them LIVE.  22 static_asserts (model 9, chain 5, combiner 3,
data 2, tree 2, facade 1) are compile-time - no gap.  The six runtime asserts:
  [[chain.hpp:4425@658869ac]] recoveredFactors || bottomNodesAreOccupied - counterpart YES,
    recoverVarianceLeafValuesBelow's numObservations() > 0 test returns 1.0.
  [[chain.hpp:4494@658869ac]]/[[chain.hpp:4522@658869ac]] leaf[j] < muByTree[t].size() - NO.  Its own comment
    concedes it: a stale index reads inside mu's CAPACITY, so it returns a
    stale value rather than faulting, invisible to ASAN too.
  [[chain.hpp:4628@658869ac]] leaf[j] < numNodes && isBottom - PARTIAL: the fused pass
    declines a map FLAGGED stale ([[chain.hpp:4617@658869ac]]), not one that is fresh and WRONG.
  [[grow.hpp:249@658869ac]], [[model.hpp:2224@658869ac]] numPresent in [2, numReachable] - NO, but the
    failure is a finite wrong prior mass, not a fault.
EMPIRICAL.  N1 makes leafOf fresh-but-wrong: plain build SIGABRTs in the GROW
suite, long before any leafOf test runs; rebuilt -DNDEBUG (R's own flag) it is
rc 1 with 53 failures, first "leafOf matches the derived map after every sweep"
([[test_moves.cpp:1456@658869ac]]).  M10, the only other SIGABRT, also reddens under -DNDEBUG
(70 failures).  VERDICT no catch here depends on an assertion the release build
lacks.  Residual exposure: the SHIPPED library has no runtime guard on the
(mu, leafOf) pairing at all - only a test does.  The sentinel shape, milder.

## Findings

1 BLOCKER [[chain.hpp:1518@658869ac]] (also [[chain.hpp:1578@658869ac]],[[chain.hpp:1616@658869ac]],[[chain.hpp:1665@658869ac]],[[chain.hpp:1678@658869ac]]) [C1].  Deleting the
  invalidateStatistics the latent refresh owes the vector leaf's U'WU cache
  passes the whole suite.  A missing CROSS: every LinearGaussianLeaf/
  GPGaussianLeaf fixture in tests/cpp is gaussian (test_model:388,
  test_moves:1341, test_sampler:385/1757/1761, test_fuzz:1054, test_shape:170),
  whose workingWeightsVaryPerSweep() is false, so the branch all five call
  sites live in is NEVER ENTERED; under logistic/nbinom/ordinal/Student-t the
  cache would serve every sweep against the PREVIOUS sweep's precisions.  The
  primitive is unit-tested (test_model:1393); the wiring is not.  REACHABLE:
  [[facade.hpp:816@658869ac]] applies no family gate, R/spec.R refuses a linear/GP node
  prior only on the multinomial ([[facade.hpp:445@658869ac]]) and K-forest ([[facade.hpp:643@658869ac]]) routes, and a
  +-25-line window around each of tinytest's 59 `node.prior = linear|gp` sites
  holds ZERO non-gaussian `family =`.  No backstop in either suite.
  ASSERTION one end-to-end fixture per vector leaf model under a
  workingWeightsVaryPerSweep family, draws compared with a cache-disabled run.
  FIX agent-fix (test); refusing the cross at the bridge instead is VD.

2 BLOCKER facade.hpp forwarding [F1,F3,F4,F6,F7].  5 of 7 facade mutations
  survive: numSavedDraws reporting capacity ([[facade.hpp:429@658869ac]]), setResponse forwarding
  !updateScale ([[facade.hpp:454@658869ac]]), savedSlotForDraw returning the identity ([[facade.hpp:514@658869ac]]),
  setForestWeights installing on forest 0 regardless ([[facade.hpp:619@658869ac]]), savedTree reading
  forest 0 regardless ([[facade.hpp:519@658869ac]]).  Only the two shape() fields were caught, by
  test_shape.cpp.  Every other test drives Sampler<L> DIRECTLY, so the
  36-virtual facade - the one dispatch layer between the shipped flat C API and
  the engine - is exercised only through shape(); the ENGINE-side twins of
  three of these WERE caught (A2, A7), isolating the gap to the forwarders.
  ASSERTION a facade conformance test: for each virtual with a selecting
  argument (forestIndex, slot, chainNum, updateScale) one call through
  SamplerBase whose answer differs per argument value.  FIX agent-fix (test).

3 MAJOR [[chain.hpp:4207@658869ac]] [C3].  Binding y + meanFits for y - meanFits - the
  variance forest fits the wrong quantity - passes.  The only VALUE-level
  variance assertion (testVarianceForestRecovery, test_model:689-725) uses
  y = s(x)*N(0,1) with NO MEAN FUNCTION, so y-f and y+f differ by 2f ~ 0 and
  the mutation is invisible by construction; its bound is loose too
  (highMean > 2.5*lowMean against a truth ratio of 44).  MAJOR not BLOCKER
  because tinytest backstops it ([[test-heteroscedastic.R:11-42@658869ac]]: fTrue = 2x, s(x)
  bounded in (0.15,0.6) and (0.9,2.2)).  [[docs/design/heteroscedastic.md:465-468@658869ac]] states gate
  (d) as "recovers f(x) AND s(x)"; the C++ fixture dropped f(x).  ASSERTION the
  same test on a strong non-constant mean, asserting s^2(x) tracks the true
  noise and NOT |f(x)|.  FIX agent-fix (test).

4 MAJOR [[sampler.hpp:575@658869ac]] (identically [[sampler.hpp:605@658869ac]], [[sampler.hpp:679@658869ac]]) [A8].  Giving predict's
  destination the CAPACITY stride while the loop runs over filledSavedDraws()
  passes.  Regressed that is an out-of-bounds heap WRITE: out is sized
  slab * filledSavedDraws() * numChains, so a partially filled store and any
  chain c >= 1 writes past the end.  Every saved-tree test is single-chain -
  testSavedDrawOrder (test_state:255) leaves numChains at 1, so
  c*numDraws == c*capacity == 0 identically, and the sanitizer cannot help
  because no fixture drives the shape.  ASSERTION repeat the partial-fill
  section with numChains >= 2, checking chain 1's slab starts at numDraws*slab,
  for all three readers.  FIX agent-fix (test).

5 MAJOR [[moves.hpp:626-706@658869ac]], the move-validity predicates [M8, M9].
  ordinalRuleIsValid's descendant interval off by one on both sides ([[moves.hpp:635@658869ac]]) and
  categoricalSubtreeIsValid's gauge dropping `directions == reachable` ([[moves.hpp:659@658869ac]])
  both survive.  ruleIsValid, ordinalRuleIsValid, categoricalSubtreeIsValid,
  categoricalSubtreeIsValidWide and findGoodOrdinalRules have ZERO direct
  callers in tests/cpp; they are reached only as a side condition of
  swap/change, where a wrong verdict moves a proposal between "no-op" and
  "scored" - visible in the acceptance ledger, in no structural check.  M6, the
  same family, was caught only by two unrelated tests reading shifted draws - a
  tripwire, not a measurement.  ASSERTION direct unit tests over hand-built
  trees: a descendant at exactly the ancestor bound accepted, one past it
  refused; a mask equal to the reachable set refused, a strict nonempty subset
  accepted; the pooled sibling on the same table.  FIX agent-fix (test).

6 MAJOR [[grow.hpp:222@658869ac]] [G5].  -log(numPredictors) for -log(numAvailable) -
  identical at a root where every variable is available, wrong at every node
  below - survives.  All three grow-from-root LAW tests measure ROOT rules only
  (testOrdinalMissingRowsAreRouted, testCategoricalExactDrawLaw,
  testCategoricalPrefixDrawLaw); below the root the suite asserts legality and
  gauge, never a probability.  ASSERTION one conditional law test at depth
  >= 1 - fix the root rule, chi-square a child's realized rules against the
  exact law - so the availability normalizer, the depth factor and the ancestor
  interval are measured off a root.  FIX agent-fix (test).

7 MAJOR [[model.hpp:3314@658869ac]] [D8].  Dropping the lower tail from
  OrdinalResponse::computeLogLikelihood - log(Phi(g_k - eta)) for the Phi
  difference - survives.  This is the EXPORTED per-observation log-likelihood
  channel of every ordinal fit; the 08-24 value scan already found four
  silently-wrong exported channels of this shape.  ASSERTION pin it against an
  independently coded Phi difference on a small ordinal fixture, as the file
  already does for the cutpoint acceptance.  FIX agent-fix (test).

8 MAJOR [[model.hpp:2185@658869ac]] [D4].  Dropping the `+ log 2` that
  ruleForVariableLogProbability owes a missing-bearing ordinal column survives.
  It enters treeLogProbability and so the change and swap ratios, but cancels
  on a same-variable redraw - which is why the grow-side twin (G1,
  [[grow.hpp:292@658869ac]]) is caught and this is not: nothing measures an MH ratio moving
  BETWEEN a missing-bearing and a plain column.  ASSERTION pin
  ruleForVariableLogProbability against a hand-computed value on both column
  kinds (it has no direct test today).  FIX agent-fix (test).

9 MAJOR [[model.hpp:3603@658869ac]] [D7].  Deleting LogisticResponse::setWeights' cold start
  of INACTIVE rows survives - the line that IS the logistic weight channel's
  landing claim ("a row that reactivates cannot carry an omega shaped by counts
  the sampler no longer holds").  No landed test constructs a count swap while
  a mask is installed.  ASSERTION install a mask, swap the counts, clear the
  mask, check the reactivated rows' latents equal a cold start against the NEW
  counts.  FIX agent-fix (test).

10 MINOR [[combiner.hpp:911@658869ac]] [B5].  formForestVetoWeights losing the near-zero
  multiplier snap survives, though testBCFZeroMultiplierSnap exists: it covers
  the response half ([[combiner.hpp:894-897@658869ac]], all three caught), not the veto half - the one
  the empty-leaf veto reads.  ASSERTION extend it to formForestVetoWeights,
  asserting a snapped row's veto precision is exactly 0.0.  FIX agent-fix.

11 MINOR diagnosis, not detection.  16 of 63 catches were crashes naming
  nothing (S2, S3, G7, M7, T1, T2, T5, N1 and the four reruns, on
  SIGSEGV/SIGBUS/SIGTRAP/SIGABRT, stdout lost to buffering), so CI reddens with
  no FAIL line.  `setvbuf(stdout, NULL, _IOLBF, 0)` in tests/cpp/main.cpp would
  make a crashing run name the check it died in.  FIX agent-fix (one line).

## Real defects found incidentally

None user-visible in pristine code; two code-health items, with reproductions.
R1 [[moves.hpp:280-281@658869ac]] - the birth-rejection restore of the parent's
  sumWeights/sumWeightedResponse is DEAD.  Tree::birth ([[tree.hpp:953@658869ac]]) writes
  the parent's rule and leftChild and the two CHILDREN's stats, never the
  parent's sums.  Reproduction: M3 deletes one and the suite passes (equivalent
  mutant); M3b replaces both with 1.0e6 and the suite reddens, so the fields
  are read later but the restore changes nothing.  The comment above reads as
  though it were load-bearing.  FIX agent-fix (delete two lines, fix comment).
R2 [[sampler.hpp:765@658869ac]] - `bool allValid = columnMaskOk;` is redundant:
  Chain::stateIsValid runs the same scratch build and the same
  columnMaskSubtreeIsValid per forest ([[chain.hpp:3229@658869ac]]) and per variance tree
  ([[chain.hpp:3318@658869ac]]), and *columnMaskRefused reads columnMaskOk directly, so `true` is
  behaviourally equivalent (A6, missed).  [[chain.hpp:3390@658869ac]] already calls
  stateIsValid "the invariant's backstop".  Noted so the MISS is not read as a
  hole.  FIX defer.

## Sanitizer leg

Pristine staged tree, `make clean && make OPT="-O2 -g -fsanitize=address,
undefined"`, run once under ASAN_OPTIONS=detect_container_overflow=0: build rc
0, run rc 0, "all tests passed", ZERO "runtime error" or "AddressSanitizer"
lines.  Clean.
