# grow-from-root-categorical-scan

agent: S0/S1/S2 opus (each moves a draw law on the opt-in path); S3 sonnet
  (R surface + tests); S4 opus (it records a measurement); S5 sonnet (docs).
  Serialized: one implementer, each slice lands before the next starts.
rng: NEUTRAL on every default path, EVERY slice. The categorical branch sits
  INSIDE the existing `j` loop of `growTreeFromRoot`, so the ordinal emission
  order and count are unchanged and an all-ordinal design draws an identical
  stream; all new state lives in `GrowScratch`. Gate: the trio bitwise at every
  tip, with the verified fact that NO trio scenario reaches `growTreeFromRoot`
  at all (zero `n.grow.sweeps` / `growFromRoot` hits anywhere under
  `benchmarks/`). A trio deviation is a bug in the slice, NEVER a re-record.
  DRAW-LAW-CHANGING on the OPT-IN grow path only, and deliberately so: a
  categorical win spends `1 + A` new Bernoullis, and under `dart = TRUE` the
  grow sweep's own DART update now sees nonzero categorical split counts - see
  "The DART consumer", which is a stream-POSITION change, not only a value
  change.
window: after the multinomial counts/offset channel
  (docs/plans/multinomial-counts-mutation.md), whose S4 re-records the
  multinomial baseline - name the hash in `benchmarks/baselines/MANIFEST` at
  this arc's start, not `ec2a3d0`. Pre-release by default. Does NOT reopen the
  `n.grow.sweeps` default (KILLED in both strata) or the arithmetic-only
  renormalization (CLOSED by measurement as cell S7).
budget: S0 ~110 tests/cpp + 1 line grow.hpp + ~10 docs; S1 ~150 scan.hpp +
  ~30 model.hpp + ~8 tree.hpp + ~90 grow.hpp + ~110 tests/cpp; S2 ~380
  tests/cpp; S3 ~180 R test; S4 ~350 R harness (gitignored) + a landing note;
  S5 ~80 docs. ~5.5 sessions plus ~90 min compute. Budgets are sized to the
  MANDATED ORACLES, not the engine delta.

## Goal

`growTreeFromRoot` places REAL categorical split rules, weighted commensurably
with the ordinal branch, so the opt-in `n.grow.sweeps` init stops paying
categorical columns into `availableSplitProbability` while drawing nothing out.
Below a cap the node draw is EXACT over the full present-partition space; above
it a Fisher-prefix family carries the variable's full prior mass, so a
variable's total realized split mass is continuous in its level count. The v1
"categorical predictors are never split here" contract inverts.

## Binding decisions inherited (do not reopen)

1. **The arithmetic-only renormalization stays CLOSED** (D1 of
   docs/plans/grow-from-root-default-study.md; matched its pre-registered null
   as cell S7; docs/design/grow-from-root-default.md section 5). Real support
   supersedes it without relitigating it: the mass categorical columns already
   contribute stops being reserved-and-discarded and starts being spent.
2. **The `n.grow.sweeps` default question stays KILLED**, both strata, same
   document. This arc reports no verdict on it.
3. **`exhaustiveCap = 10`**, 12 rejected. Raising the cap relocates the mass
   boundary rather than fixing it, and 511 candidates already cover every
   factor an R user constructs. One compile-time constant.
4. **Deterministic mass-spread above the cap**, not Horvitz-Thompson. HT is
   unbiased in expectation and pathological per draw: at P = 20 the complement
   holds 524268 partitions, so with M = 16 a sampled complement partition
   carries 32767x a prefix's weight and dominates the node draw whenever it is
   sampled. It also adds `M*P` Bernoullis to the RNG contract.
5. **Draw the A absent-position coins; do not pin them left.** Pinning makes
   every grown mask a subset of the present set, biasing where categories
   unseen at that node route at PREDICT time, in a direction the MH sweeps do
   not correct.
6. **No new `CGMTreePrior` or `ColumnStore` members.** All new state in
   `GrowScratch`.
7. Design artifacts (memo, adversarial critique, synthesis) are durable at
   `.claude/grow-from-root-categorical-scan/`. **Read the synthesis before
   starting** - it carries the corrected enumeration recipe, the mass-repair
   derivation, and the re-aimed measured arm, none of which is in the memo.

## Context (seams, all read at 54c0b9f / engine 06c0254)

- `growTreeFromRoot` (grow.hpp) assembles one discrete draw over `{no-split} U
  {(variable, cut)}`; the whole defect is the `data.types[j] ==
  ColumnType::categorical` `continue` sitting downstream of the
  `availableSplitProbability` sum. The categorical branch goes INSIDE this loop.
- `scanOrdinalCuts` (scan.hpp) is NOT touched - that is what keeps ordinal-only
  designs bitwise. Its `cutScanEmptySentinel`, its `binScratch.assign()` clear,
  and its `naCode` skip are all ordinal-only behavior.
- The scan marginal omits `sum wz^2` deliberately (scan.hpp header): additive
  over any partition of a fixed member set, so it cancels cut-vs-cut AND
  split-vs-no-split. That cancellation is what puts categorical candidates,
  ordinal candidates and the no-split term on one scale with no extra
  bookkeeping.
- `Rule` (tree.hpp) is `{int32 variableIndex, uint64 bits}`; for a categorical
  column `bits` IS the direction mask, inline for K <= 63 and a `maskPool` word
  offset for pooled columns (`maskWordsForCount`, data.hpp). There is NO
  separate missing coin for a categorical: a missing value takes the reserved
  pseudo-category (`naCategory` = 63 inline, K pooled - `missingCategoryCode`,
  data.hpp) and its direction is another bit of the same mask.
- `CGMTreePrior::ruleForVariableLogProbability` (model.hpp) is CGM's uniform
  over the `2^R - 2` non-degenerate masks, with an `R > 54` branch because
  `pow(2,R) - 2` stops being exact. `drawRuleForVariable` /
  `drawCategoryPattern` draw the pattern BIT BY BIT deliberately (generator
  granularity pins low pattern bits for wide masks - docs/design/pooled-masks.md)
  and scatter it onto reachable positions ASCENDING; unreachable positions stay
  zero, the canonical gauge.
- `buildFromFlat` rejects a mask that is zero, not a subset of reachable, or
  equal to reachable (tree.hpp, inline and pooled). Every rule this arc builds
  must satisfy all three.
- `Tree::reachableCategories` / `reachableCategoriesWide`;
  `collectAvailableVariables` narrows every INLINE categorical's reachable mask
  on its single O(p + depth) walk (`availMaskScratch_`, private) and falls back
  to `variableAvailable` for pooled columns. `variableAvailable` needs
  `popcount(reachable) >= 2`.
- **Invariant this design leans on:** every category PRESENT among a node's
  members is REACHABLE there (a member arrives only by satisfying every
  ancestor rule, so its bit survives every AND). Hence `A = R - P >= 0`. Debug
  assert.
- `maxCategories = 0xFFFF` and `numCuts[j]` holds K for a categorical column,
  so **K in 1..65535**. But the size governing the scan is not K: it is `R`
  (reachable, <= K+1) and `P` (distinct categories PRESENT at the node,
  <= min(R, numMembers)). P collapses fast with depth. Every cost below is in P.
- Design docs: docs/design/grow-from-root.md (the standing contract, and where
  the v1 categorical paragraph lives), docs/design/pooled-masks.md (the two
  mask tiers), docs/design/grow-from-root-default.md (the closed study),
  docs/design/interaction-constraints.md (grow is the FIFTH split generator).

## The weight, and the mass repair

Fix a node and an available categorical variable j. `R` = popcount(reachable,
including the missing pseudo-category iff `hasMissing[j]`), `P` = distinct
present categories, `A = R - P >= 0`. The prior is uniform over the `2^R - 2`
non-degenerate masks. Grouping by the induced partition of the present set:
`2^A` masks per oriented side, `2^(A+1)` per unordered partition,
`2^(P-1) - 1` partitions, and the residual `2^(A+1) - 2` masks are exactly the
empty-child rules. A variable's whole enumerable family therefore carries

    M(R, P) = (2^(P-1) - 1) * 2^(A+1)/(2^R - 2) = (1 - 2^(1-P)) / (1 - 2^(1-R))

and the plan's single weight rule is **spread that family mass evenly over
whatever the branch emits**:

    logRuleCat = log1p(-exp2(1-P)) - log1p(-exp2(1-R)) - log(numEmitted)

- Exact branch, `numEmitted = 2^(P-1) - 1`: this reduces IDENTICALLY (<= 1 ulp,
  verified) to `(1-P) log 2 - log1p(-exp2(1-R))`, the per-partition group mass.
  Implement the exact branch in that closed form - it is better conditioned
  than the shipped `pow(2,R) - 2` and needs no `R > 54` branch (`exp2(1-R)`
  underflows to zero harmlessly at `R -> 65536`, where the shipped form returns
  `inf`).
- Fisher branch, `numEmitted = P - 1`: the variable's TOTAL mass is `M(R, P)`
  for every P, so it is continuous across the cap. Without this the total drops
  from 1.000 at P = 10 to 0.009775 at P = 11 - a 102x cliff at a compile-time
  constant the user cannot see.

Bias, recorded honestly: within the Fisher branch the retained greedy prefixes
inherit the family's whole mass, so relative to the ideal `sum over all
partitions of pi * L` the variable is over-weighted by at most
`(2^(P-1)-1)/(P-1)` - 51x at P = 11. That bound is attained only in the
strong-signal regime, where the likelihood ratio is `exp(O(n_node))` and a 51x
prior factor does not decide. The exact conditional errs by the SAME factor in
the opposite direction in the weak/no-signal regime, where the prior IS what
decides. That asymmetry is the whole reason for the choice.

The full candidate weight, exactly parallel to the ordinal branch:

    splitBase_j = logGrowth + logSplitVariable_j + logRuleCat_j
    logW(s)     = splitBase_j + logL(W_L(s), S_L(s)) + logL(W_R(s), S_R(s))

with `logSplitVariable_j` taken unchanged (categoricals already pay into
`availableSplitProbability`, so nothing there moves).

## The enumeration (verified exhaustively; do not re-derive)

Sort the P present categories by `(S_c / (a s2 + W_c), code)` ascending, where
`a = (k/scale)^2` and `s2 = residualVariance`. That key is exactly the
singleton-category leaf posterior mean `ConstantGaussianLeaf::drawFromPosterior`
would compute - the analog of LightGBM's `cat_smooth` and XGBoost's output leaf
value - and `a s2 > 0` always, so the `W_c == 0` case disappears. The code
tie-break makes the order a deterministic function of the inputs, which the
post-draw reconstruction needs. Both branches use this one order.

- `P < 2`: emit nothing, return. The variable still counts in `numAvailable`
  and `availableSplitProbability`, exactly as an ordinal column all of whose
  cuts are occupancy-empty already does. No new asymmetry.
- `P <= exhaustiveCap`: **the plain binary counter `s = 0 .. 2^(P-1) - 2` over
  the P-1 non-anchor categories, anchor (sorted position 0) always LEFT; left
  side = `{anchor} U s`, bit b of `s` naming sorted position b+1.** Recompute
  both sides in O(P) from the compact present array per candidate, so the
  suffstat is a pure function of `s` - no path-dependent add/subtract sequence,
  no last-ulp drift into a negative `sumWeights`, no NaN in
  `posteriorPrecision`, and a trivial decode. Do NOT use a Gray code: the memo's
  `G(g) = g ^ (g>>1)` over `g = 0..2^m-2` emits the FULL non-anchor mask (empty
  right child) for every `P >= 3` and drops one legitimate partition, and the
  measured O(P) recompute costs ~+1 us per (node, variable) at P = 10 against a
  5.62 us enumeration.
- `P > exhaustiveCap`: the `P - 1` prefixes `t = 1..P-1` of the sorted array.

**Keep a count-based sentinel anyway**, one compare per candidate:
`left.count == 0 || left.count == total.count -> cutScanEmptySentinel`. The
enumeration domain makes it unreachable, which is the point - it converts a
recipe bug into an undrawable candidate. It cannot be replaced by reasoning
about scores: `logIntegratedLikelihood` returns 0.0 at `sumWeights == 0`
(model.hpp), so a legal all-zero-weight side and an illegal empty side score
identically.

## Post-draw rule assembly (fixed order, part of the RNG contract)

1. Recompute the compact present-category array and its sorted order for the
   winning variable - one O(n + P log P) pass, deterministic, no draws, so no
   state need be carried through the discrete draw.
2. Draw 1 Bernoulli: orientation. `side = orientation ? complement(S) : S`
   within the present set.
3. Walk the reachable positions ASCENDING; draw 1 Bernoulli for each ABSENT
   one, matching `drawCategoryPattern`'s bit-by-bit convention and its stated
   generator-granularity reason.
4. Set the mask: present-side bits plus the drawn absent bits; unreachable
   positions stay 0. Pooled columns `tree.allocateMask(maskWordsForCount(K))`
   first and write through `mutableMaskWordsFor`. Grow never rejects, so no
   pool mark/truncate; `growForestFromRoot` already calls
   `compactMaskPoolIfNeeded` per grown tree under `hasPooledCategorical`.

Exactly `1 + A` Bernoullis, with NO rejection loop - unlike
`drawCategoryPatternWide`, because the present side already guarantees both
children non-empty. Assert the count. The resulting mask is nonzero, a strict
subset of reachable, and not equal to reachable, so it passes the gauge at
`buildFromFlat` and both children are occupied.

## The DART consumer

`Chain::growForestFromRoot` closes each sweep with, per forest,
`if (forest.useDart) { memset splitCounts; trees[t].countVariableUses(...);
forest.dart.update(rng_, ...) }`, and `countVariableUses` is type-agnostic.
Today categorical columns receive STRUCTURALLY ZERO split counts during every
warm-start sweep. After this arc they receive real ones, so
`bart2(..., dart = TRUE, n.grow.sweeps = k)` - a shipped, non-hypothetical
combination that no test anywhere covers - changes.

Stronger than a value change: `DartPrior::update` calls
`ext_rng_simulateGamma(rng, alpha/p + splitCounts[j], 1.0)` per predictor, and
that sampler is a REJECTION sampler whose uniform consumption depends on the
shape (src/external/random.c). A count moving off zero can shift the stream
POSITION for the remainder of the fit.

**Posture: do not gate it.** A DART posterior fed real categorical counts is
more correct than one fed structurally-zero counts, and gating would make the
init less complete under DART for no correctness reason. Test it (S3), document
it (NEWS plus the design note), name it here.

## Constraints

- No `facade.hpp` virtual moves, so the stale-object bus-error hazard does not
  apply - but `grow.hpp`, `scan.hpp`, `model.hpp` and `tree.hpp` are ALL
  headers pulled into the bridge TU by `src/bartcore/bartcore.hpp`, so **every
  `R CMD INSTALL` in this arc needs `--preclean`**.
- `scanOrdinalCuts` is not touched, and neither is `binScratch`. The
  categorical kernel gets its OWN `categoryBins` in `GrowScratch`.
- **Do not `assign()` the category bins.** Bin count is
  `hasMissing[j] ? missingCategoryCode(K) + 1 : K`, up to 65536 x 24 bytes; an
  O(K) clear at every node of every tree would dominate. Keep a `touched` list
  of first-seen codes plus a `seen` flag array cleared only over `touched`.
  Keep it unconditionally - one code path is worth more than a branch on K.
- **Do not materialize a mask per candidate**; a pooled column would need 1024
  words each. `candidateCut` carries `s` (exact) or `t` (Fisher), both well
  inside `int32`, and which one is implied by the winning variable's own
  recomputed P.
- **Grow owns its pooled reachable-mask scratch in `GrowScratch`.** It must not
  borrow `Tree::reachableScratch_`, which `variableAvailable` and
  `buildFromFlatBelow` also use. `collectAvailableVariables` narrows only
  INLINE masks, so grow needs its own `reachableCategoriesWide` + `maskPopcount`
  per available pooled column to get R.
- **Label the new kernel INIT-ONLY in its header comment.** A candidate set
  that is Fisher-truncated and mass-reweighted is not a valid MH neighborhood,
  and the present-partition grouping is not either (the absent-position
  assignment affects the reverse move's reachable sets). Without the label a
  future docs/design/parallel-bart-frontier.md section 3.1 implementer reuses it
  and breaks detailed balance silently.
- Out of scope: any `inst/include/dbarts/dbarts.h` edit; the renormalization;
  the `n.grow.sweeps` default; measuring whether the shrunken sort key beats the
  raw mean (Q5 - after the mass repair the family total does not depend on which
  prefixes are retained); `man/bart.Rd`'s stale "replaced with dummies" line for
  `bart2` (real, confirmed, its own ticket).

## S0. The ordinal `log 2` convention: falsifier, then pin

Decides a question, then acts on the answer. Lands first so S1 writes the
commensurability rule against a settled convention. **VD signs before the fix
half lands** (see Open items).

1. Two-arm falsifier in `tests/cpp/test_grow.cpp`, on a small missing-bearing
   ordinal fixture. Arm A = the shipped code. Arm B = an explicit
   `{no-split} U {(c, missing-left), (c, missing-right)}` enumeration where
   each of the `2*numRules` candidates carries its OWN exact likelihood with
   the missing rows PLACED on the stated side (no scan approximation). Arm B is
   the ideal law on the larger enumerated set by construction. Chi-square both
   arms' root-rule frequencies against arm B's own weights, ~2e5 seeds.
2. Pre-registered decision rule, written into the test file before it is run:
   arm A rejecting while arm B does not CONFIRMS that the shipped weight
   implements the smaller-set convention and deflates the group by exactly 2;
   arm A not rejecting REFUTES the memo's reading and the shipped weight stands.
3. If confirmed and signed: delete the `logCut -= std::log(2.0)` line. Keep the
   falsifier, re-pointed so arm A is the new code and now agrees with arm B.
4. Either way, PIN the convention: one sentence in `growTreeFromRoot`'s
   draw-discipline comment naming which rule set an enumerated candidate stands
   for, and a matching note in docs/design/grow-from-root.md. The convention is
   what S1's categorical branch is written against.

Gate: `cd tests/cpp && make && ./test_bartcore`; the trio bitwise;
full tinytest with NO snapshot regenerated. The regeneration prediction is
ZERO and is a HARD gate, not an expectation: the three tinytest files reaching
`growFromRoot` all source `inst/common/friedmanData.R`
(`matrix(runif(n*10), n, 10)`, no factors, no NAs) and
test-grow-from-root.R's heteroscedastic block is `cbind(runif, runif)`, so no
recorded value can move. A shift means the change leaked into the non-missing
ordinal path - STOP, never regenerate.

## S1. The kernel, the weight, and the install path

1. `scan.hpp`: `scanCategoricalPartitions`, the categorical sibling of
   `scanOrdinalCuts`, INIT-ONLY in its header comment. Pass 1 is one
   O(numMembers) sweep branching once on dense/sparse into `categoryBins` with
   the `touched`/`seen` bookkeeping, compacted into a dense P-entry array of
   `(code, count, sumWeights, sumWeightedResponse)`. No `naCode` skip - for a
   categorical column a missing value is a real code. Pass 2 is the sort, then
   the enumeration of "The enumeration" above with the count-based sentinel.
2. `model.hpp`: a `static` helper beside `ruleForVariableLogProbability`, e.g.
   `CGMTreePrior::categoricalGroupLogProbability(R, P, numEmitted)`, so the two
   derivations sit together. Precondition assert `2 <= P && P <= R`; the helper
   is public API for the component test and is meaningless outside that range
   (it returns a "probability" above 1 at `P < 2`).
3. `tree.hpp`: `const uint64_t* inlineReachableMasks() const`, valid
   immediately after `collectAvailableVariables`, so grow does not re-walk
   O(p * depth) for inline columns. Pooled columns still need grow's own
   per-node wide call.
4. `grow.hpp`: the categorical branch INSIDE the existing `j` loop, in the same
   `j` order, emitting after that variable's ordinal candidates would have
   been; the new `GrowScratch` members; post-draw assembly per "Post-draw rule
   assembly"; the debug assert `A >= 0`; the rewritten draw-discipline comment
   (`1 discrete` + ordinal-with-missing `1 Bernoulli` + categorical
   `1 + A Bernoullis`, no rejection loop), with the "categorical predictors are
   ordinal-only in v1" paragraph deleted.
5. Invert `testCategoricalNeverSplit` (test_grow.cpp) rather than deleting it,
   and refresh the file header comment and the two draw-count helper comments,
   which go stale with it.
6. `tests/cpp/test_scan.cpp`: the exhaustive enumeration obligation. For
   `P = 2..12`, assert the emitted set has exactly `2^(P-1) - 1` entries, no
   duplicate induced partition, both sides non-empty at every entry, and a
   bijection onto an independently brute-forced enumeration of the unordered
   nonempty-both-sides partitions built in the test. This is the check that
   would have caught the Gray-code recipe.

Gate: tests/cpp from clean, plain AND the ASAN+UBSAN leg (new reachable
pointer arithmetic: `mutableMaskWordsFor` writes, a bin array indexed to 65536,
the `touched`/`seen` bookkeeping); `R CMD INSTALL --preclean .`; full tinytest,
no snapshot regenerated; the trio bitwise against the MANIFEST hashes current
at the arc's start; `benchmarks/R/categorical-exact.R` unchanged-pass (it is
the MH-path analog of S2's O1 and must not move). Delete the
`benchmarks/kernels` binaries after the header edits - no dependency tracking
there. `benchmarks/kernels/grow_from_root.c` is standalone C including only
`misc/*`; if a categorical kernel variant is ever added it must keep that
property. `bench-sampler` is not a gate (the opt-in path is in no timed arm).

## S2. The oracles

Its own slice per the standing lesson that oracles are the budget.

1. **O1 (primary, correctness).** In test_grow.cpp: one 4-level categorical
   column (P = R = 4, so the cap is not engaged and the draw is EXACT), small
   n, one tree, depth-1 only. Compute the 7 present-partition weights plus the
   no-split weight independently INSIDE the test from the closed form and
   `logIntegratedLikelihood`; run `growTreeFromRoot` over ~2e5 seeds;
   chi-square the empirical root-rule frequencies. Falsifies the weight, the
   orientation coin, and the candidate ordering in one shot.
2. **O1F (Fisher branch).** The same shape at P = 12 with `exhaustiveCap`
   forced to 10: 11 candidates, closed-form conditional, chi-square. Plus the
   closed-form MASS INVARIANT, checked directly with no sampling: for a grid of
   `(R, P)` spanning both branches and the cap boundary, the sum of the
   emitted candidates' group masses equals `(1 - 2^(1-P))/(1 - 2^(1-R))`. Plus
   the `A = 0` cross-check `helper(R,R,2^(R-1)-1) - ruleForVariableLogProbability
   == log 2`, at a RELATIVE tolerance - it holds analytically in both branches
   of the shipped form but not bitwise, since one path goes through
   `pow(2,R) - 2` and the other through `log1p(-exp2(1-R))`.
3. **O2 (gauge and legality).** Over many seeds on a K = 6 inline and a K = 70
   pooled fixture (mirroring inst/tinytest/test-data-categorical-wide.R), each
   with an NA-bearing variant: for every internal node the mask is nonzero, a
   strict subset of `reachableCategories`, both children occupied; the A absent
   positions really are drawn (not all-zero across seeds) and the coin count is
   exactly `1 + A`; a `buildFromFlat` round trip of the grown tree is accepted.
   Plus one assertion on categorical x interaction-constraint x grow, newly
   reachable and today uncovered (grow is the FIFTH split generator the
   interaction design named). Do NOT copy `testGrownTreeWellFormed`'s gauge
   loop, which calls `splitInterval`/`splitIndex` unguarded by column type -
   safe on its one-ordinal fixture, wrong on a categorical rule where `bits` is
   the mask or a pool offset.
4. **Coverage gaps closed.** An `OP_GROW` op in `test_fuzz.cpp`'s `FuzzOp`,
   which already carries the `cat4` and `cat3Miss` configs - the cheapest
   possible oracle for the flat round trip and for BL-1's failure mode. A
   POOLED column added to `testMappedSourceReplay` (test_tree.cpp), which
   exercises ordinal and `categoricalInline` only and becomes load-bearing now
   that grow emits pooled categorical rules.

Gate: as S1. Record a K = 70 timing note from O2.

## S3. R surface, and the DART assertion

1. New `inst/tinytest/test-grow-from-root-categorical.R`: a categorical-heavy
   design through `dbarts(...)$growFromRoot(2L)` where the init forest MUST now
   contain categorical rules on the signal factors, read via `getTrees` /
   `extract(fit, "trees")`; and a `bart2(..., n.grow.sweeps = 2L)` early-RMSE
   assertion shaped like test-bart2-grow-from-root.R's but on a factor-signal
   design, where the current code cannot beat cold at all.
2. **The DART assertion**: `bart2(..., dart = TRUE, n.grow.sweeps = 2L)` on a
   factor-signal design gives the signal factors nonzero posterior split
   probability after the grow sweeps. This is the only place the interaction is
   observable, and nothing covers it today.
3. Record the `exhaustiveCap` sizing observation (the enumeration is not the
   cost driver - the O(numMembers) histogram is, and it is the same cost the
   ordinal path already pays) in the landing note. No slice is spent measuring
   the touched-list bookkeeping at K = 8; it is a few lines and immaterial
   either way.

Gate: as S1, plus `air format --check .`.

## S4. The measured arm (three arms, all constants pinned)

Harness under `.claude/grow-from-root-categorical-scan/` (gitignored). The
recorded S7 numbers (`B1@t1 = -0.7287` etc.) are the motivation for this arm
and appear in NO threshold: they were measured at `k* = 16` (this arc's arms
are k = 2/8) from a harness that exists at no commit, against a
pre-registration that pins only the scenario's shape.

**Arms, per replicate, all on the same paired seed:** cold
(`n.grow.sweeps = 0`); skip@k=2; cat@k=2; skip@k=8; cat@k=8. The SKIP arms are
the same binary with a compile-time `kGrowCategoricalScan = false`, so the two
warm arms differ in exactly one thing and the floor is measured IN SITU.

**Pinned DGP** (the pre-registration's largest hole; all of this is new):
n = 5000 train plus 200 held out from the same DGP; 10 predictors, x1..x5
unordered factors of 8 levels sampled uniform over levels (all inline - the
pooled path is O2's job, not this arm's), x6..x10 ordinal `runif(n)`;
`mu_i = sum_{j=1..5} E[j, level_ij]` where **E is a FROZEN LITERAL 5x8 matrix**
pasted into the harness, generated once by
`set.seed(20260809); round(matrix(rnorm(40), 5, 8), 4)`, never re-drawn per
replicate; `sigma = 1.0`; `factors = "categorical"`. Record the realized
`sd(mu)`. m = 75, 8 chains, 2000 draws, `n.burn = 0`, `combineChains = FALSE`,
`n.threads` pinned and stated. Seeds follow `grouped-mixing.R`: data
`set.seed(20260809 + s)`, sampler `seed = s`, shared `s` across all five arms
of replicate s. R = 12, `s = 1..12`.

**Metrics** exactly as defined in the study: B1 at `t in {1, 25, 100}`
(`(RMSE_cold - RMSE_warm)/RMSE_cold` on per-iterate test RMSE vs true mu), B2
(`log2(B_warm/B_cold)`), C1 and C2 (relative gaps over the last 500 iterates).

**Gates.** k = 2 is PRIMARY, k = 8 is DESCRIPTIVE with no gate, so there is no
combination rule and no alpha inflation.

    G1 PRIMARY   paired D_s = B1@t1(cat, s) - B1@t1(skip, s) > 0,
                 one-sided alpha 0.05, R = 12
    G2 REPORTED  B1@t1(cat) > 0 (beats cold) - secondary, NOT gating
    G3 GUARDRAIL C1, C2 for the cat arm inside the study's frozen SMALL
                 margins (0.03 and 0.06 relative)
    G4 REPORTED  paired B2(cat) vs B2(skip)

Power, at the re-aimed target: conservative planning SD for the paired
difference is `sqrt(2) * 0.1352 = 0.191` - an upper bound, since both warm arms
share the same cold arm (which enters B1 only through the common denominator)
and the same data realization. That gives min-detectable `0.091` at R = 12 and
80% power at `0.137`, i.e. a bar of a 12.5% recovery of the 0.73 gap. If
categorical support puts S7 where the study's other cells sit (B1@t1 in
0.21..0.87) then `D = 0.94..1.60`, z = 10 to 18 conservatively. **Report the
realized paired SD.** If it exceeds the planning bound, the result is reported
UNDERPOWERED; the threshold does not move.

**Mandatory fresh-seed re-run** (the study's step 5, carried explicitly rather
than inherited by reference - a premature-GREEN was its recorded failure mode):
any gate failing or any guardrail breached is re-run at `s = 101..112` before
the flag is treated as confirmed.

Wall clock enters no criterion, so this arm does not need a quiet machine.
~90 min compute.

## S5. Docs and landing

docs/design/grow-from-root.md: the v1 categorical contract replaced, the
enumeration and weight rule stated once, the INIT-ONLY label explained, the
pinned `log 2` convention from S0, and the DART consumer. A landing note on
docs/design/grow-from-root-default.md section 5's categorical bullet (real
support landed; the renormalization stays closed; no verdict on the default).
The `Status:` line bumped to `LANDED <date> (<commit>)`. TODO, `inst/NEWS.Rd`,
and the `## Landing` note in this file.

## Verification (every slice)

- `R CMD INSTALL --preclean .` - `grow.hpp`, `scan.hpp`, `model.hpp` and
  `tree.hpp` are all headers in the bridge TU. Delete the `benchmarks/kernels`
  binaries after any header edit; `tests/cpp` tracks headers via `-MMD -MP`, so
  plain `make` is correct there.
- `cd tests/cpp && make clean && make && ./test_bartcore` - all pass. ASAN+UBSAN
  leg for S0-S2 (each makes new engine or new reachable code); S3-S5 add none.
- Full `tinytest::test_package("dbarts")` from a preclean install. New tests
  ADD; NO snapshot is regenerated at any slice. A forced snapshot is a signal
  the slice changed more than intended - stop and report.
- The trio, EVERY slice, expecting no deviation:
  `benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-c8f661a.rds`
  -> 27/27 "identical draws (same RNG stream)";
  `bcf-equivalence-99205ee.rds` -> 5 scenarios x 6 channels bitwise; the
  multinomial baseline named in `benchmarks/baselines/MANIFEST` at the arc's
  start -> bitwise on every channel. No max-|z| line anywhere. **THE TRIO IS
  NECESSARY, NOT SUFFICIENT, and here it is unusually weak**: no trio scenario
  reaches `growTreeFromRoot` at all (verified - zero `n.grow.sweeps` /
  `growFromRoot` hits under `benchmarks/`), including `equivalence.R`'s
  `categorical` scenario, which is pure MH. The trio proves only that nothing
  leaked into the default path. O1/O1F/O2 and S4 are the real oracles.
- `benchmarks/R/categorical-exact.R` unchanged-pass at every slice.
- `air format --check .` on any slice touching R/. No bridge move in this arc,
  so no rchk obligation.

Stop conditions per docs/plans/README.md, plus:

1. Any trio scenario not reporting exact match -> STOP, the change leaked into
   the default path. Never re-record to clear it.
2. Any tinytest failure outside the new file, or any of the three grow tinytest
   files SHIFTING -> STOP. The zero-regeneration prediction is a gate.
3. O1 or O1F's chi-square rejecting, or the closed-form mass invariant failing
   -> STOP; the weight or the enumeration order is wrong. Do not tune the
   fixture.
4. G1 failing after the fresh-seed re-run -> STOP and report. The only ladder
   step is raising `exhaustiveCap`; HT augmentation is NOT a fallback
   (binding decision 4).
5. More than 2 sessions inside any single slice -> escalate rather than
   continue.

## Falsifiers

- **F1 (S0).** The two-arm falsifier must be written so arm B can disagree with
  arm A. If both arms pass the same chi-square, the fixture has no missing
  predictor on a splittable column and proves nothing - `testDeterminism-
  AndDrawCount` asserts `hasMissing[0] == 0` precisely because the existing
  fixtures cannot see this.
- **F2 (S1).** The exhaustive enumeration check must be shown RED against the
  memo's Gray recipe `G(g) = g ^ (g>>1)` over `g = 0..2^m-2`. A check that has
  never been red is not a check.
- **F3 (S1).** The count-based sentinel must be shown to fire when the
  enumeration range is deliberately widened to include the full non-anchor mask
  - that is the exact defect it exists to convert into an undrawable candidate.
- **F4 (S2).** O1 must be shown red under a deliberate `log 2` perturbation of
  `logRuleCat` and under a swapped orientation coin.
- **F5 (S2).** The mass invariant must be shown red against the memo's
  unrepaired Fisher weight (`(1-P) log 2` above the cap), which is where the
  102x cliff lives.
- **F6 (S2).** `testMappedSourceReplay`'s pooled arm must be shown red when the
  pooled mask channel is bypassed.
- **F7 (S4).** The skip arm must be verified to place ZERO categorical rules
  and the cat arm nonzero, on the same seed, before any contrast is computed. A
  three-arm harness whose two warm arms are secretly the same binary measures
  nothing.

## Edge cases the tests must name

A node where `R >= 2` but `P < 2` (available, emits nothing, still counts in
`availableSplitProbability` - the same behavior an all-empty-cut ordinal column
already has). `P = 2` (one candidate; the counter's `m = 1` case, the one the
Gray recipe was accidentally right about). Empty-but-reachable categories at a
node, `A > 0`, including deep nodes where `A` grows because absent categories
are routed down by ancestor masks. The cap boundary, P = 10 and P = 11, both
branches and the mass continuity across them. Pooled, K = 70, including a
pooled column whose node-level P is below the cap. The missing pseudo-category:
inline at position 63 and pooled at position K, present at the node (a real
histogram bin, routed by the partition) and absent at the node (routed by an
absent-position coin). A present category with count >= 1 but all-zero weights
(legal, scores 0.0, must NOT be sentineled - occupancy is on counts). A
categorical variable vetoed by an interaction constraint. `dart = TRUE` with
`n.grow.sweeps > 0`. A design where every categorical column has `P < 2` at the
root (grow must still terminate and produce a legal forest).

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S0: (only if the fix lands) the grow-from-root initializer no longer
  under-weights cut candidates on ordinal predictors that carry missing values.
- S1: the `n.grow.sweeps` initializer now places categorical split rules
  instead of skipping categorical predictors; below ten categories present at a
  node the node's choice is exact over all category subsets, and above it a
  sorted-prefix family carries the same total prior mass.
- S3: with `dart = TRUE`, categorical predictors now receive split counts
  during the grow sweeps, so their sampled split probabilities entering the
  first regular sweep differ from previous versions.

## Open items

- **The `log 2` convention needs VD's signature before S0's fix half lands.**
  The falsifier settles the fact; VD signs the consequence. Fix it and the two
  branches obey one rule - a candidate carries its rule group's total prior
  mass, post-draw coins pick uniformly within the group - at the cost of a
  deliberate draw-law change on the opt-in path for missing-bearing ordinal
  columns (zero snapshot burden, verified; magnitude `log 2 = 0.69` against
  likelihood ratios of `exp(O(n_node))`, so it bites only at the shallow-signal
  small-node margin). Or pin the shipped weight as the smaller-enumerated-set
  convention and document it, at the cost of the two branches obeying different
  rules. Recommendation: fix.
- **The mass-spread bias direction is RECORDED, not scheduled.** Above the cap
  the retained greedy prefixes inherit the family's whole mass, which is the
  ESL 9.2.4 many-level-predictor overfitting direction. S7's 8 levels sit below
  the cap, so this arc's measured arm does not exercise it. Reopening trigger:
  an observed C1/C2 guardrail breach traced to a factor with more than
  `exhaustiveCap` levels, or a user report of a high-cardinality factor being
  over-selected at init. The remedy if triggered is a shrunk spread (a factor
  between 1 and `(2^(P-1)-1)/(P-1)`), not HT.
- `man/bart.Rd`'s "is replaced with dummies" is stale for `bart2` (that is
  `bart()`'s forced `factors = "indicators"`; `bart2`, `dbarts`, `dbartsData`,
  `rbart_vi` and `xbart` all default to `factors = "categorical"`). Confirmed,
  out of this arc, worth its own ticket - an outside reader otherwise concludes
  dbarts dummy-expands.
- Pre-existing hazard found in the census, unrelated and unreached here:
  `storeVarianceSavedTrees` and `rebuildVarianceForest` (chain.hpp) call
  `flatten`/`buildFromFlat` with `masks == nullptr`, which would null-deref if a
  variance-forest tree ever held a pooled categorical rule. Grow does not sweep
  the variance forest, so this arc does not reach it. Its own ticket.
