# Mutation testing leg B: long-form evidence (companion to mutation-B-findings.md)

## Ladder audit

Probe: a 7th ResponseFamily enumerator in a staged copy; tests/cpp rebuilt and
R_interface_bartcore.cpp compiled -fsyntax-only, -Wall -Wextra.  Exactly THREE
sites warned (-Wswitch): [src/bartcore/chain.hpp:582](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L582) (single-forest response factory),
R_interface:2268 (refusedAmplitudeFamilyReason), and a third,
[src/R_interface_bartcore.cpp:6244](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L6244) (`drawAugmentationLaws`).
Four `default:` ladders stayed silent, plus one open-coded chain no -Wswitch
can reach.  (Only two of the four are in src/bartcore; the memo's list spans
both files.)
  [src/bartcore/chain.hpp:756](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L756)  AmplitudeSpec ctor        aft/ord/nbinom/new  -> GaussianResponse
  [src/bartcore/chain.hpp:5026](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L5026) latentScaleAnchor         gauss/aft/ord/nb/new -> scaledResponseSd()
  R_interface:2291 defaultNodeScale        gaussian/aft/new    -> 0.5 (gaussian's)
  R_interface:2812 validateResponseSupport gaussian/aft/new    -> NO validation
  R_interface:6254 computeWorkingResponse  if/else, no switch  -> latent[i]
By call graph, [src/bartcore/chain.hpp:756](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L756) and [src/bartcore/chain.hpp:5026](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L5026) take spec.family from the K-forest route, which
R_interface:2268 refuses for aft/ordinal/nbinom BEFORE either runs, so those
arms are gaussian-only today.  [src/R_interface_bartcore.cpp:2812](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L2812) is the worst door - its own comment
enumerates the harm it prevents (a negative nbinom count underflowing into a
~1.8e19 allocation, "an uncatchable crash, not an error").  The de-facto 7th
family already exists: [src/bartcore/chain.hpp:882](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L882) sets `family_ = logistic` for multinomial
commenting "family() is not read on this path", yet SamplerShape::family
([src/bartcore/facade.hpp:438](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L438)) reports logistic and the bridge branches on it at six sites
(unresolved: [src/R_interface_bartcore.cpp:2729](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L2729), unresolved: [src/R_interface_bartcore.cpp:4611](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L4611), unresolved: [src/R_interface_bartcore.cpp:4662](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L4662), unresolved: [src/R_interface_bartcore.cpp:4702](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L4702), unresolved: [src/R_interface_bartcore.cpp:4877](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L4877), unresolved: [src/R_interface_bartcore.cpp:4949](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/R_interface_bartcore.cpp#L4949) - the exact line within the bridge that each reference names could not be re-derived) - each unreachable only because a SEPARATE
multi-forest refusal stands in front, a coupling nothing asserts.
VERDICT the doors are shut by refusals in another file, not by the type system,
and one family already travels under another's name.

## Assertion audit

R compiles packages with -DNDEBUG (Makeconf:174), so every runtime assert() in
src/bartcore is DEAD in the shipped library; tests/cpp never defines it, so the
suite runs with them LIVE.  22 static_asserts (model 9, chain 5, combiner 3,
data 2, tree 2, facade 1) are compile-time - no gap.  The six runtime asserts:
  [src/bartcore/chain.hpp:4425](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4425) recoveredFactors || bottomNodesAreOccupied - counterpart YES,
    recoverVarianceLeafValuesBelow's numObservations() > 0 test returns 1.0.
  [src/bartcore/chain.hpp:4494](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4494)/[src/bartcore/chain.hpp:4522](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4522) leaf[j] < muByTree[t].size() - NO.  Its own comment
    concedes it: a stale index reads inside mu's CAPACITY, so it returns a
    stale value rather than faulting, invisible to ASAN too.
  [src/bartcore/chain.hpp:4628](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4628) leaf[j] < numNodes && isBottom - PARTIAL: the fused pass
    declines a map FLAGGED stale ([src/bartcore/chain.hpp:4617](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4617)), not one that is fresh and WRONG.
  [src/bartcore/grow.hpp:249](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/grow.hpp#L249), [src/bartcore/model.hpp:2224](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/model.hpp#L2224) numPresent in [2, numReachable] - NO, but the
    failure is a finite wrong prior mass, not a fault.
EMPIRICAL.  N1 makes leafOf fresh-but-wrong: plain build SIGABRTs in the GROW
suite, long before any leafOf test runs; rebuilt -DNDEBUG (R's own flag) it is
rc 1 with 53 failures, first "leafOf matches the derived map after every sweep"
([tests/cpp/test_moves.cpp:1456](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/tests/cpp/test_moves.cpp#L1456)).  M10, the only other SIGABRT, also reddens under -DNDEBUG
(70 failures).  VERDICT no catch here depends on an assertion the release build
lacks.  Residual exposure: the SHIPPED library has no runtime guard on the
(mu, leafOf) pairing at all - only a test does.  The sentinel shape, milder.

## Findings

1 BLOCKER [src/bartcore/chain.hpp:1518](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L1518) (also [src/bartcore/chain.hpp:1578](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L1578),[src/bartcore/chain.hpp:1616](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L1616),[src/bartcore/chain.hpp:1665](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L1665),[src/bartcore/chain.hpp:1678](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L1678)) [C1].  Deleting the
  invalidateStatistics the latent refresh owes the vector leaf's U'WU cache
  passes the whole suite.  A missing CROSS: every LinearGaussianLeaf/
  GPGaussianLeaf fixture in tests/cpp is gaussian (test_model:388,
  test_moves:1341, test_sampler:385/1757/1761, test_fuzz:1054, test_shape:170),
  whose workingWeightsVaryPerSweep() is false, so the branch all five call
  sites live in is NEVER ENTERED; under logistic/nbinom/ordinal/Student-t the
  cache would serve every sweep against the PREVIOUS sweep's precisions.  The
  primitive is unit-tested (test_model:1393); the wiring is not.  REACHABLE:
  [src/bartcore/facade.hpp:816](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L816) applies no family gate, R/spec.R refuses a linear/GP node
  prior only on the multinomial ([src/bartcore/facade.hpp:445](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L445)) and K-forest ([src/bartcore/facade.hpp:643](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L643)) routes, and a
  +-25-line window around each of tinytest's 59 `node.prior = linear|gp` sites
  holds ZERO non-gaussian `family =`.  No backstop in either suite.
  ASSERTION one end-to-end fixture per vector leaf model under a
  workingWeightsVaryPerSweep family, draws compared with a cache-disabled run.
  FIX agent-fix (test); refusing the cross at the bridge instead is VD.

2 BLOCKER facade.hpp forwarding [F1,F3,F4,F6,F7].  5 of 7 facade mutations
  survive: numSavedDraws reporting capacity ([src/bartcore/facade.hpp:429](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L429)), setResponse forwarding
  !updateScale ([src/bartcore/facade.hpp:454](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L454)), savedSlotForDraw returning the identity ([src/bartcore/facade.hpp:514](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L514)),
  setForestWeights installing on forest 0 regardless ([src/bartcore/facade.hpp:619](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L619)), savedTree reading
  forest 0 regardless ([src/bartcore/facade.hpp:519](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/facade.hpp#L519)).  Only the two shape() fields were caught, by
  test_shape.cpp.  Every other test drives Sampler<L> DIRECTLY, so the
  36-virtual facade - the one dispatch layer between the shipped flat C API and
  the engine - is exercised only through shape(); the ENGINE-side twins of
  three of these WERE caught (A2, A7), isolating the gap to the forwarders.
  ASSERTION a facade conformance test: for each virtual with a selecting
  argument (forestIndex, slot, chainNum, updateScale) one call through
  SamplerBase whose answer differs per argument value.  FIX agent-fix (test).

3 MAJOR [src/bartcore/chain.hpp:4207](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L4207) [C3].  Binding y + meanFits for y - meanFits - the
  variance forest fits the wrong quantity - passes.  The only VALUE-level
  variance assertion (testVarianceForestRecovery, test_model:689-725) uses
  y = s(x)*N(0,1) with NO MEAN FUNCTION, so y-f and y+f differ by 2f ~ 0 and
  the mutation is invisible by construction; its bound is loose too
  (highMean > 2.5*lowMean against a truth ratio of 44).  MAJOR not BLOCKER
  because tinytest backstops it ([inst/tinytest/test-heteroscedastic.R:11-42](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/inst/tinytest/test-heteroscedastic.R#L11-L42): fTrue = 2x, s(x)
  bounded in (0.15,0.6) and (0.9,2.2)).  [docs/design/heteroscedastic.md:465-468](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/docs/design/heteroscedastic.md#L465-L468) states gate
  (d) as "recovers f(x) AND s(x)"; the C++ fixture dropped f(x).  ASSERTION the
  same test on a strong non-constant mean, asserting s^2(x) tracks the true
  noise and NOT |f(x)|.  FIX agent-fix (test).

4 MAJOR [src/bartcore/sampler.hpp:575](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L575) (identically [src/bartcore/sampler.hpp:605](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L605), [src/bartcore/sampler.hpp:679](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L679)) [A8].  Giving predict's
  destination the CAPACITY stride while the loop runs over filledSavedDraws()
  passes.  Regressed that is an out-of-bounds heap WRITE: out is sized
  slab * filledSavedDraws() * numChains, so a partially filled store and any
  chain c >= 1 writes past the end.  Every saved-tree test is single-chain -
  testSavedDrawOrder (test_state:255) leaves numChains at 1, so
  c*numDraws == c*capacity == 0 identically, and the sanitizer cannot help
  because no fixture drives the shape.  ASSERTION repeat the partial-fill
  section with numChains >= 2, checking chain 1's slab starts at numDraws*slab,
  for all three readers.  FIX agent-fix (test).

5 MAJOR [src/bartcore/moves.hpp:626-706](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/moves.hpp#L626-L706), the move-validity predicates [M8, M9].
  ordinalRuleIsValid's descendant interval off by one on both sides ([src/bartcore/moves.hpp:635](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/moves.hpp#L635)) and
  categoricalSubtreeIsValid's gauge dropping `directions == reachable` ([src/bartcore/moves.hpp:659](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/moves.hpp#L659))
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

6 MAJOR [src/bartcore/grow.hpp:222](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/grow.hpp#L222) [G5].  -log(numPredictors) for -log(numAvailable) -
  identical at a root where every variable is available, wrong at every node
  below - survives.  All three grow-from-root LAW tests measure ROOT rules only
  (testOrdinalMissingRowsAreRouted, testCategoricalExactDrawLaw,
  testCategoricalPrefixDrawLaw); below the root the suite asserts legality and
  gauge, never a probability.  ASSERTION one conditional law test at depth
  >= 1 - fix the root rule, chi-square a child's realized rules against the
  exact law - so the availability normalizer, the depth factor and the ancestor
  interval are measured off a root.  FIX agent-fix (test).

7 MAJOR [src/bartcore/model.hpp:3314](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/model.hpp#L3314) [D8].  Dropping the lower tail from
  OrdinalResponse::computeLogLikelihood - log(Phi(g_k - eta)) for the Phi
  difference - survives.  This is the EXPORTED per-observation log-likelihood
  channel of every ordinal fit; the 08-24 value scan already found four
  silently-wrong exported channels of this shape.  ASSERTION pin it against an
  independently coded Phi difference on a small ordinal fixture, as the file
  already does for the cutpoint acceptance.  FIX agent-fix (test).

8 MAJOR [src/bartcore/model.hpp:2185](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/model.hpp#L2185) [D4].  Dropping the `+ log 2` that
  ruleForVariableLogProbability owes a missing-bearing ordinal column survives.
  It enters treeLogProbability and so the change and swap ratios, but cancels
  on a same-variable redraw - which is why the grow-side twin (G1,
  [src/bartcore/grow.hpp:292](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/grow.hpp#L292)) is caught and this is not: nothing measures an MH ratio moving
  BETWEEN a missing-bearing and a plain column.  ASSERTION pin
  ruleForVariableLogProbability against a hand-computed value on both column
  kinds (it has no direct test today).  FIX agent-fix (test).

9 MAJOR [src/bartcore/model.hpp:3603](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/model.hpp#L3603) [D7].  Deleting LogisticResponse::setWeights' cold start
  of INACTIVE rows survives - the line that IS the logistic weight channel's
  landing claim ("a row that reactivates cannot carry an omega shaped by counts
  the sampler no longer holds").  No landed test constructs a count swap while
  a mask is installed.  ASSERTION install a mask, swap the counts, clear the
  mask, check the reactivated rows' latents equal a cold start against the NEW
  counts.  FIX agent-fix (test).

10 MINOR [src/bartcore/combiner.hpp:911](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/combiner.hpp#L911) [B5].  formForestVetoWeights losing the near-zero
  multiplier snap survives, though testBCFZeroMultiplierSnap exists: it covers
  the response half ([src/bartcore/combiner.hpp:894-897](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/combiner.hpp#L894-L897), all three caught), not the veto half - the one
  the empty-leaf veto reads.  ASSERTION extend it to formForestVetoWeights,
  asserting a snapped row's veto precision is exactly 0.0.  FIX agent-fix.

11 MINOR diagnosis, not detection.  16 of 63 catches were crashes naming
  nothing (S2, S3, G7, M7, T1, T2, T5, N1 and the four reruns, on
  SIGSEGV/SIGBUS/SIGTRAP/SIGABRT, stdout lost to buffering), so CI reddens with
  no FAIL line.  `setvbuf(stdout, NULL, _IOLBF, 0)` in tests/cpp/main.cpp would
  make a crashing run name the check it died in.  FIX agent-fix (one line).

## Real defects found incidentally

None user-visible in pristine code; two code-health items, with reproductions.
R1 [src/bartcore/moves.hpp:280-281](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/moves.hpp#L280-L281) - the birth-rejection restore of the parent's
  sumWeights/sumWeightedResponse is DEAD.  Tree::birth ([src/bartcore/tree.hpp:953](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/tree.hpp#L953)) writes
  the parent's rule and leftChild and the two CHILDREN's stats, never the
  parent's sums.  Reproduction: M3 deletes one and the suite passes (equivalent
  mutant); M3b replaces both with 1.0e6 and the suite reddens, so the fields
  are read later but the restore changes nothing.  The comment above reads as
  though it were load-bearing.  FIX agent-fix (delete two lines, fix comment).
R2 [src/bartcore/sampler.hpp:765](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/sampler.hpp#L765) - `bool allValid = columnMaskOk;` is redundant:
  Chain::stateIsValid runs the same scratch build and the same
  columnMaskSubtreeIsValid per forest ([src/bartcore/chain.hpp:3229](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L3229)) and per variance tree
  ([src/bartcore/chain.hpp:3318](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L3318)), and *columnMaskRefused reads columnMaskOk directly, so `true` is
  behaviourally equivalent (A6, missed).  [src/bartcore/chain.hpp:3390](https://github.com/vdorie/dbarts/blob/658869ac6a90fc23dc9d0860e0aa0890e743c6d7/src/bartcore/chain.hpp#L3390) already calls
  stateIsValid "the invariant's backstop".  Noted so the MISS is not read as a
  hole.  FIX defer.

## Sanitizer leg

Pristine staged tree, `make clean && make OPT="-O2 -g -fsanitize=address,
undefined"`, run once under ASAN_OPTIONS=detect_container_overflow=0: build rc
0, run rc 0, "all tests passed", ZERO "runtime error" or "AddressSanitizer"
lines.  Clean.
