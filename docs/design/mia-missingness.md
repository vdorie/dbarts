# MIA missingness

LANDED 2026-07-04, as proposed. Deltas discovered while landing: the
Rule union became a single uint64 with accessors outright (the moves and
tests read through splitIndex()/categoryDirections()); validateXTest's
model.frame needed na.pass like ingestion's (test data.frames were
silently na.omit-dropping rows); the starting sigma estimate mean-imputes
NAs for its linear model only (complete cases are scarce under scattered
missingness); all-NA predictor columns are rejected at ingestion (no
observed values to split on); and a nonexistent row-name subset now
surfaces as "response contains missing values" (na.pass keeps the all-NA
rows the old na.omit swallowed). Gates: component tests (ingestion,
mechanics, end-to-end signal recovery + bitwise state round trip),
tinytest test-data-missing.R, equivalence identical draws on all nine
complete-data scenarios, speed within noise.

Original proposal follows. It builds on public-surface.md section 2,
whose DECIDED paragraph reserved the representation room this design
spends: per-column cut counts cap one below the code type's maximum so a
reserved NA code always fits, and mask bits 53-63 are unused (the
categorical cap is 53 and ordinal rules' high union word is invariantly
zero).

## Model

Missing Incorporated in Attributes (Twala et al. 2008): every split rule
gains a missing direction, so an observation whose split variable is NA
follows that rule's chosen side. The direction is part of the rule,
proposed from a symmetric prior and accepted through the ordinary MH
ratio, so the chain learns where missingness should flow per split.
Missingness is a modeling feature only for predictors: NAs in the
response, weights, or offset stay rejected.

## Representation

- Ordinal codes: NA quantizes to the reserved code 0xFFFF (xint_t max).
  Enforce numCuts <= 0xFFFE - 1 at build so real codes (which reach
  numCuts) never collide. Cut construction (uniform ranges, quantile
  grids) skips NaN; a column must retain at least one observed value.
- Categorical codes: NA takes the fixed category position 63, i.e. code
  63 and mask bit 63. It deliberately does NOT take position K: the
  reachable-mask machinery (Tree::reachableCategories,
  categoricalSubtreeIsValid, drawCategoryPattern) then handles NA routing
  with no special cases - seed the root's reachable mask with bit 63 when
  the column has NAs and NA becomes one more "category" that ancestor
  rules filter down the tree, in canonical gauge, with the pattern draw
  covering it uniformly. Popcount-driven code paths count it correctly.
- Ordinal rules: the missing direction lives in bit 63 of the Rule union
  (the invariantly-zero high word). Do not rely on union aliasing to
  preserve it across splitIndex writes - that is endianness-dependent.
  Instead make the uint64 the sole storage with splitIndex()/
  setSplitIndex(index, missingGoesRight) accessors; rules change only by
  whole-Rule assignment already, and equals() compares the wide member,
  so rule identity picks up the missing direction for free.
- Per-column NA presence: ColumnStore grows hasMissing (one flag per
  column), set during build/quantize and refreshed by the mutation
  surface. Columns without NAs never draw or store a missing direction
  (the bit stays canonical zero), which is what keeps NA-free data
  bitwise identical to today - no extra RNG draws, no kernel changes on
  the fast path.

## Kernels

- Rule::sendsRight: ordinal becomes
  code == NA_CODE ? missingGoesRight : code > splitIndex; categorical is
  unchanged (NA is just code 63 against the mask).
- partitionChildren: the categorical mask partition is unchanged. The
  ordinal SIMD kernels (misc_partitionRange/misc_partitionIndices) order
  by code <= splitIndex, which happens to send NA_CODE right; when the
  column has NAs and the rule sends missing left, use a scalar MIA-aware
  partition instead. Branch per node on data.hasMissing[var] - NA-free
  columns keep the SIMD path unconditionally.
- Flat replay (countFlatObservationsBelow, predict on raw doubles): NaN
  routes by the flat node's missing direction; all comparisons against
  NaN are false today, which silently sends NaN left on the x <= value
  test - that path must consult the flag before comparing.

## Moves

- Birth (drawRuleAndVariable): ordinal rules on columns with NAs draw the
  missing direction as an extra Bernoulli(1/2); the proposal is from the
  prior so its density cancels, exactly like the cut draw. Categorical
  needs no extra draw (the pattern covers bit 63).
- Change: redraw the missing direction with the new cut (same
  cancellation). Swap moves rules wholesale, so directions travel with
  them; ordinal missing bits cannot invalidate a swap, and categorical
  NA bits flow through the existing reachable-mask validity checks.
- Identifiability: a subtree no NA can reach has meaningless ordinal
  missing bits; the prior is symmetric and the likelihood flat over
  them, so the chain wanders harmlessly (the categorical side is gauged
  by the reachable mask and does not have this slack). Accept that
  rather than gauge ordinal bits - change/swap validity would otherwise
  need NA-reachability tracking for no inferential gain.
- Variable availability and tree priors are untouched: availability
  still counts satisfiable cuts/categories, and the growth probability
  does not see missingness.

## Flat format and serialization

FlatNode grows a flags byte (bit 0 = missing goes right); value keeps
the cut double or the observed-category mask confined to bits 0-52. The
NA direction cannot ride in the categorical mask double: bit 63 (or any
bit above 52) pushes the mask past 2^53 and doubles stop round-tripping
integers exactly there. State objects are opaque and engine-specific, so
tree.vars/tree.values grow a parallel tree.flags (raw or integer)
element, absent meaning all-zero for restores of older states within the
same major version. buildFromFlat validates flags (only bit 0, zero when
the column lacks NAs).

## Bridge and R surface

- Ingestion: dbartsData accepts NAs in predictors (train and test) when
  the engine does. Add missing = c("incorporate", "error") to dbartsData
  and thread it through dbarts/bart2/xbart/rbart_vi like factors. The
  value is named for what happens - missingness is incorporated into the
  split rules - rather than for the MIA acronym (DECIDED 2026-07-04,
  Vincent: the name and the default are both "incorporate"; "error" is
  the strictness escape). bart keeps rejecting NAs (classic
  makeModelMatrix surface).
- Mutation surface: setPredictor/updatePredictor/setData/setCell accept
  NaN, refresh hasMissing, and re-quantize; a column that gains NAs
  mid-run routes them by the existing rules' canonical-zero bit (left)
  until moves revisit those rules - deterministic and documented.
  cutsWouldRemainValid and categoricalValueIsValid treat NaN as valid
  when missing = "incorporate".
- Views (buildFromParent) copy hasMissing with the rest of the grid; the
  reserved codes gather like any other.
- Reporting: getTrees gains a missing column ("L"/"R", NA where the rule
  has no missing direction), mirroring the directions decode; the
  categorical directions string does NOT grow a 54th character (the NA
  route reports in the missing column uniformly for both rule kinds).
  printTrees prints the direction with the rule. plotTree appends the
  route to the label ("NA->L").

## Gates

NA-free data must be bitwise identical to 1.0-0: component tests
comparing a fit before/after the feature on complete data, the full
equivalence compare (identical draws), and the speed benchmark (the
hasMissing branch should be unmeasurable; verify). New coverage: MIA
fits recover a signal carried by missingness (y depends on is.na(x1)),
NA routing consistency between getTrees counts, predict, and an R-side
replay, state round-trip with flags, views over NA data, and test-data
NAs. xbart/rbart_vi/bart2 smoke tests with NAs once the argument
threads through.

## Out of scope

Imputation-based alternatives (bartMachine-style hot-decking), NA in the
response, and sparse-aware storage (public-surface section 7 keeps
sparse columns separate; MIA's reserved codes are dense and orthogonal).
