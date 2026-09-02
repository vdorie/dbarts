# Pooled category masks

Lifts the 53-level categorical cap (public-surface.md section 7; the
reporting decision it waited on landed 2026-07-04 as the getTrees
directions column). After this, categorical predictors support up to
65535 levels - the code type's limit - with no change of any kind to
columns of 53 or fewer levels: those paths stay textually identical and
the equivalence gate must keep reporting identical draws.

LANDED 2026-07-04, as proposed. Deltas and notes:

- Building the wide component tests surfaced a pre-existing bitwise-restore
  fragility unrelated to masks: restoreScale re-anchors the variance
  prior's internal scale through a multiply-divide round trip
  (scale * range^2 / range^2), which is not a floating-point identity, so
  a restored gaussian chain's sigma draws could drift one ulp from the
  original's for data-dependent values. The state now carries the prior's
  internal scale exactly (ChainStateData.sigmaPriorScale, the third
  element of the R state's fit.scale slot; two-element slots from older
  states still restore, with the old last-ulp caveat).
- The R state slots are tree.masks/saved.masks as proposed; restore
  splits the concatenated words by walking each tree's vars against the
  store's category counts (maskWordsForFlatTree).
- The changeMove proposal reuses one scratch pool allocation across its
  64 rejection-sampling attempts and takes its pool mark before the draw
  so every exit path (abort, MH reject) truncates uniformly; the ordinal
  branch's truncate is a no-op.
- Component tests: testPooledMaskMechanics (partition across words,
  reachability, draw-margin uniformity, closed-form rule probability,
  compaction, flat round trip and its malformed rejections, NA-code
  composition) and testPooledMaskSampler (recovery over K = 70 pooled
  plus K = 60 inline-band columns, saved replay equal to recorded test
  fits, bitwise state round trip, live prediction). tinytest
  test-data-categorical-wide.R covers the R surface end to end.

## Two tiers

The 53 cap conflated two limits: the rule word holds 64 direction bits
(one reserved for the missing direction), while the flattened format's
double encoding holds 53. They are lifted separately.

- **Inline, K <= 63**: a rule's mask lives in its 64-bit word exactly as
  today, category c at bit c, the missing direction at the fixed
  naCategory position 63. The only engine change this band needs is
  drawCategoryPattern handling numReachable = 64 (63 categories plus the
  missing pseudo-category; `1ull << 64` is undefined), a guarded
  all-ones constant on the existing path.
- **Pooled, 64 <= K <= 65535**: the rule word holds an offset into a
  per-tree pool of mask words instead of the mask itself. A pooled
  column's mask is ceil((K + 1) / 64) words: category c at bit c and the
  missing direction at bit K, so the partition kernel tests bit `code`
  uniformly and the reachable-mask machinery routes missing values with
  no special case - the reserved missing code for a pooled column is K
  (naCategory's 63 would collide with a real category), assigned
  per-column by codeFor and buildFromParent.

A column is pooled iff maskWordsForColumn(categoryCounts[j]) > 1, i.e.
ceil((K + 1) / 64) > 1. Everything keys off that predicate; ordinal
columns and narrow categorical columns never touch the pool.

## Pool ownership and lifetime

Each Tree owns `std::vector<std::uint64_t> maskPool_`. Rules store word
offsets (stable across vector growth), pool entries are immutable once a
rule in the tree references them, and Rule copies alias freely - which
makes subtree snapshots, undoBirth, and node-struct restores correct
with no bookkeeping. Allocation discipline:

- A move records the pool size before drawing; every rejection path
  (draw abort, MH reject) truncates back to the mark. The change move's
  rejection-sampling loop overwrites one scratch allocation in place -
  nothing references it until the rule is installed.
- An accepted change or death strands the replaced rule's words. The
  chain compacts a tree's pool at a safe point (after its metropolis
  step returns, when the pool exceeds a high-water mark of
  max(256, 4x the live size found by the last compaction)): walk the
  arena nodes, copy live masks, rewrite offsets. No snapshot is live
  then, trees are small, and the walk is deterministic - no RNG
  involvement, so draws are unaffected.
- Tree::initialize and buildFromFlat clear and rebuild the pool; Tree
  copies carry the pool by value, keeping offsets valid.

Multi-word helpers (popcount, test/set bit, and/andnot into a buffer,
equality, is-zero) are free functions next to Rule. The recursive
categorical-validity walk gets a wide variant threading a depth-indexed
mask arena through MoveScratch; the right-child branch passes the rule's
own pooled words (immutable), only the left side needs a fresh slot per
level. reachableCategories gets a wide variant writing into a
caller-supplied buffer; Tree carries a mutable scratch for the
compute-check-discard call sites (variableAvailable, buildFromFlat's
gauge check). Rule equality for the swap move compares pooled words
through a Tree helper (offsets alone would distinguish equal masks).

The rule-count log-probability -log(2^R - 2) switches to the closed form
R log 2 + log1p(-2^(1 - R)) when R > 54, where pow(2, R) - 2 is no
longer exact; existing data only produces R <= 54, so existing draws are
bit-identical, and the formula is a deterministic function of R, so MH
ratios stay consistent. drawRuleForVariable takes Tree& non-const (a
pooled draw allocates); the wide pattern draw extends the bit-by-bit
Bernoulli scheme unchanged - one draw per reachable category ascending,
all-same patterns rejected - so the narrow stream shape is the wide
stream shape.

## Flat format

The flattened encoding is keyed at 53 - the double-exactness boundary -
independent of the engine tier:

- Rules on columns with K <= 53 keep today's encoding exactly (mask in
  value, missing direction in flags): byte-identical formats, states
  from f87a65c restore unchanged.
- Rules on wider columns (54 <= K, whether the engine mask is inline or
  pooled) store in value the rule's word offset into a per-tree side
  channel `std::vector<std::uint64_t>` - offsets are small integers,
  double-exact, and self-describing, so replay indexes the channel
  directly with no cursor threading through the pre-order recursion (the
  slopes side array needs (numOnLeft + 1) / 2 bookkeeping; masks do
  not). Each such rule contributes maskWordsForColumn(K) words holding
  category bits only, the missing-direction bit cleared: the missing
  direction stays in flags for every rule kind, single source of truth,
  and the getTrees missing column is untouched. flatten emits offsets
  sequentially in pre-order; buildFromFlat validates them as such,
  bounds-checks against the channel, and re-checks the canonical gauge
  against the reachable mask like the narrow path.
- Tree::flatten/buildFromFlat and the flat replay family
  (partitionFlatIndices, countFlatObservationsBelow,
  addFlat(Linear)PredictionsBelow, flatSubtreeIsWellFormed,
  printFlatSubtree) gain defaulted numCategories/maskWords parameters;
  scalar-narrow call sites compile the identical arithmetic. Chain
  threads the channel everywhere it holds flat trees: savedTreeMasks_
  parallel to savedTrees_, ChainStateData.treeMasks/savedTreeMasks
  parallel to treeParams/savedTreeParams, allocated only when the store
  has any wide categorical column (a bool computed at build; category
  counts are fixed for the life of the store).
- SamplerBase/SamplerFacade: flattenTree gains a defaulted masks out
  parameter; savedTreeMasks(chain, sample, tree) accessor mirrors
  savedTreeSlopes.

## State objects

Per-chain slots tree.masks/saved.masks (RAWSXP, each word written
explicitly little-endian - raw bytes serialize portably where NaN
payloads in a REALSXP might not), concatenated across trees in tree
order and split on restore by walking the vars arrays against the
store's category counts, exactly as tree.params splits by leaf counts.
Absent slots mean no wide rules; states holding wide rules require them
(stateIsValid). State objects remain opaque and engine-specific.

## Reporting

- getTrees: a narrow categorical rule's value stays the raw mask
  (documented, unchanged). A wide rule's value is NA and the bridge
  emits a directions string for it C-side - one L/R per *observed*
  category, from the store's count - in a directions column added
  whenever the store has wide categorical columns. The R decode then
  only fills narrow rules from value (today's math, byte-identical) and
  right-pads every string with L to the declared level count
  (unobserved trailing levels are canonically left under the gauge -
  exactly what the narrow decode already yields for them).
- plotTree reads the directions column and needs no change; printTrees
  and the live print dump read pooled words where they read the rule
  word.
- The C API's dbarts_sampler_getTrees returns the same data.frame and
  inherits the behavior; no signature changes.

## R surface

- makeModelMatrixFromDataFrame/dbartsData: the 53-level factor error
  becomes a 65535-level error (codes must fit the engine's uint16).
- Bridge validation ([0, 53) code range) widens to match; the xbart
  data-handle path needs nothing (views copy category counts, the pool
  is per-tree and data-independent).
- man/dbarts.Rd's "(at most 53)" and dbartsSampler-class.Rd's value
  paragraph gain the new cap and the wide-value-is-NA note; NEWS entry.

## Non-goals

- Per-column code widths (uint8 hot layer) and sparse columns stay
  separate (core-generalization data model).
- No new split vocabulary: pooled rules are the same subset splits with
  wider masks; DART, MIA, and linear leaves compose untouched.

## Gates

Component tests (wide fixture, K = 70 pooled and K = 60 inline-band):
partition vs brute-force routing, gauge invariants under the moves, pool
compaction preserving masks, pattern-draw uniformity smoke, flatten/
build and state round trips bitwise, saved replay equal to recorded test
fits, MIA composition (NA routes by bit K). tinytest: 60- and 70-level
factors end to end via factors = "categorical" (fit, getTrees directions
width and value NA, state round trip, plotTree, predict on newdata).
Full suite; equivalence identical draws on all nine scenarios; speed
compare; R CMD check --as-cran.
