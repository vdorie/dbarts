# Consolidating the predictor store's semantic-type axis

Status: DESIGNED AND RULED 2026-09-01, implementation pending. Every
alternative below has been weighed and every open question ruled (see
"Decisions"); nothing has landed yet. Slices land in the ruled order S0,
S1, S2, S3, S4a, S4c, S4b (see "Slice decomposition and sizing"); each
appends its own note under "Landing" at EOF as it lands. TODO:
column-kind-consolidation.

rng: draw-preserving (NEUTRAL) on S0, S1, S3, S4a, S4c and S4b, subject to
the four hazards enumerated in section 6's hot-path subsection.
POSTERIOR-CHANGING on S2 alone - the ordered-factor midpoint cut grid
changes the split-candidate set for every such column, hence the
stationary distribution over trees for any fit that has one - which
needs its own baseline re-record plus the exact-posterior gates
(.github/workflows/exact-gates.yaml) beside the usual equivalence trio.

Scheduled pre-release-candidate. Everything below is measured against the
bartcore branch at 2da4f101 and verified current through 48219484 (the
commits between are docs-only; three anchors into docs/design/data-store.md
moved and are corrected in place below).

Two things are settled going in and are treated as scope rather than as
questions: the **ordered-factor cut-grid fix lands pre-RC** (section 1, with its
bill in section 6), and the ABI cost of any enumerator or struct append - a
`DBARTS_C_API_HASH` re-bake and a lockstep stan4bart rebuild - **is accepted**
and is stated where it applies rather than argued.

## What is actually wrong

Four separate defects share one root, and they are worth separating because
they have different fixes and very different risk.

**1. There is no ordered-factor marker anywhere below R's model-matrix
builder.** `makeCategoricalModelMatrix` collapses `is.ordered()` into
`ORDINAL_VARIABLE` (`R/utility.R:510-514`), which is `0L`
(`R/data.R:1`), the same value a numeric column gets (`R/utility.R:530-532`).
The bridge states the consequence in as many words: "an ordered factor carries
a level table too, but enters as an ordinal column, whose grid is cut points
rather than categories" (`src/R_interface_bartcore.cpp:883-885`). The engine's
`ColumnType` has two values (`src/bartcore/data.hpp:93`), so nothing below the
bridge can tell an ordered factor from a numeric column.

This is exactly the gap the response side closed with `response.type`
(`R/A_class.R:500-505`), whose comment gives the reason verbatim: "0/1/2 codes
are indistinguishable from integer data".

The splitting *mechanic* is deliberate and correct -
`docs/design/public-surface.md:123-126` rules that ordered factors ingest as
ordinal on their integer codes "because their order is meaningful; threshold
splits are the right vocabulary" - and nothing here changes it. What the
missing marker costs is four concrete things:

  - *A silent data defect.* An ordered factor's grid is built by
    `fillCutsUniformly` -> `fillCutsOverRange` (`data.hpp:797-817`) with
    `numCuts[j] = maxNumCuts[j]` (`data.hpp:881`), over the observed range
    `[0, K-1]`. That cap is `n.cuts`, default `100L`
    (`R/A_class.R:258`, `R/data.R:15`), recycled per column with no
    cardinality clamp (`R/data.R:38`, `R/dbarts.R:684`), passed through
    `options.maxNumCutsPerVariable` (`R_interface_bartcore.cpp:1775`) into
    `maxNumCuts` (`data.hpp:1149-1150`) and reduced only by the
    representability cap (`:1155-1157`). Codes come from `lower_bound` over
    those cuts (`data.hpp:713-716`). With K distinct level codes and only
    `numCuts + 1` reachable code values, K > 101 guarantees by pigeonhole that
    two adjacent levels share a code and can never be separated by any rule.
    Nothing warns. An *un*ordered factor with too many levels errors and names
    the cap (`R/utility.R:500-507`); the ordered one does not, because the cap
    check is gated on `!is.ordered(column)`.

  - *A split-prior distortion, not merely wasted work.* Below that threshold
    the column carries `n.cuts` cut points where K-1 are distinct.
    `scanOrdinalCuts` scores every candidate (`src/bartcore/scan.hpp:89-98`),
    and grow-from-root pushes every non-sentinel entry with weight
    `splitBase + logLikelihood`, where `splitBase` carries
    `logCut = -log(right - left + 1)` normalized over the *inflated* index
    range (`src/bartcore/grow.hpp:270`, `:283-296`). Each distinct partition
    therefore appears with multiplicity equal to the number of index slots that
    map to it, and the induced prior over distinct thresholds is proportional
    to that multiplicity. The multiplicities are equal only when `K-1` divides
    `n.cuts`. The distortion is mild for small K (K=4 at `n.cuts=100` gives
    slot counts 33/34/33, about 3 percent) and becomes severe as K approaches
    `n.cuts`: at `K-1 = 51` and `n.cuts = 100` adjacent thresholds carry
    multiplicities 2 and 1, a factor-of-two prior tilt between neighbouring
    level boundaries. On top of the prior effect, the scan is O(members) +
    O(numCuts) per (node, variable), so the wasted scoring is invisible at the
    root and dominant in deep, small nodes - which is where grow-from-root
    spends most of its node visits.

  - *No test-side level check.* `validateCategoricalPredictors` gates on
    `columnTypes[j] != categorical` (`R_interface_bartcore.cpp:1704`,
    `:1728`), and `validateTestContainerAgainstStore` on
    `store.types[j] != categorical` (`:3028`), so an ordered-factor test column
    is never checked against the training level table. A test value of `7.5`,
    or code `40` on a 5-level factor, is quantized against the grid and routed
    without complaint. **Disposition: fixed, by decision 5 - strict refusal
    matching categorical**, at the sites and with the wording the categorical
    check already uses. The permissive alternative (treating a between-levels
    value as interpolation on an ordered scale) was considered and rejected;
    the reasoning is recorded with the decision.

  - *Degraded reporting.* `decodeCategoricalSplits` maps split values back to
    level names only for categorical columns (`R/plotTree.R:7-25`), which it
    selects via `which(varTypes == CATEGORICAL_VARIABLE)` (`:10`, `:25`). A
    split on an ordered factor reports a numeric threshold like `2.0198` where
    it could report the level boundary. **Disposition: deferred, not in any
    slice.** It is pure R, touches no engine surface, and is cosmetic; under
    the midpoint grid the thresholds become exact half-integers (`k + 0.5`),
    which makes the decode a one-line mapping whenever someone wants it. Worth
    a TODO line rather than scope here.

**2. `numCuts` is overloaded.** The declaration says so: "For categorical
columns numCuts holds the (fixed) category count, cutPoints stays empty"
(`data.hpp:508-510`, field at `:538`, written at `:872`). There are exactly 62
`numCuts[...]` references in `src/bartcore/*.hpp` and
`src/R_interface_bartcore.cpp`; TWENTY-SEVEN of them read it as a category
count: `data.hpp:555`, `:710`, `:723`, `:1466`; `tree.hpp:320`, `:330`,
`:445`, `:473`, `:503`, `:650`, `:1178`, `:1181`, `:1475`, `:1536`, `:1561`,
`:1611`, `:2069`, `:2083`, `:2148`; `scan.hpp:303`, `:304`;
`R_interface_bartcore.cpp:1521`, `:3029`, `:5483`, `:5663`, `:7620`, `:7624`,
plus the WRITE at `data.hpp:872` and one test-side read,
`tests/cpp/test_grow.cpp:860`, which sizes a category histogram from it and
so must move with them. (An earlier draft of this list named `moves.hpp:404`,
`:708`, `model.hpp:2158`, `:2349`, `grow.hpp:84` and `:230`; those lines hold
no `numCuts[` at all - they are `types[j] == categorical` comparisons, which
the `splitsBySubset` sweep carries instead. Corrected against the tree at
S1's implementation.) Every read is already inside a `types[j] ==
categorical` guard, which is what makes the split mechanical. Three deserve
attention: `scan.hpp:303-304` sizes the categorical scan's histogram and is on
the **hot path**, `data.hpp:1466` is inside `buildFromParent`, which section 1
returns to, and `R_interface_bartcore.cpp:1521` is the verbose printer's
untyped line, which the recommendation below makes kind-aware.

**3. A categorical column's raw double IS its code.** `codeFor` does an
unchecked `static_cast<xint_t>(value)` for a categorical (`data.hpp:708-711`),
and the flat replay reads the raw double and casts it to a category index
(`tree.hpp:1768-1769` pooled, `:1783-1784` inline). The bridge widens R's
`INTEGER` factor codes to doubles a cell at a time
(`R_interface_bartcore.cpp:543-555`) so that the engine can narrow them back.

**4. `FlatKind` duplicates type-plus-pooled** (`tree.hpp:171-176`), and both
`buildFromFlatBelow` (`tree.hpp:1530-1574`) and `flatSubtreeIsWellFormed`
(`tree.hpp:2065-2088`) re-derive the same thing from the store and refuse a
disagreement. See section 4 - this one is less redundant than it looks.

## Two premises worth settling before the sections

**The flat C API is not an unvalidated back door.** `dbarts_sampler_create`
takes SEXPs (`inst/include/dbarts/dbarts.h:484-486`; the header states
"control/model/data/state cross as SEXP" at `:199-200`) and forwards to the
same `bartcore_bridge::createHolder` the R entry uses
(`src/C_interface.cpp:542-548`), so it inherits `validateCategoricalPredictors`
in full. Every other flat entrance that can deliver predictor values validates
before `codeFor` casts: `dbarts_sampler_setPredictor` and
`dbarts_sampler_updatePredictor` sweep every cell through
`validateColumnValues` (`C_interface.cpp:757`, `:790`, the function defined at
`R_interface_bartcore.cpp:2965-2972`); `setTestPredictors` and `predict` go
through `validateTestSource` (`C_interface.cpp:275-279`); `setState`
re-quantizes only creation-validated values; the `getTrees` replay clamps
defensively (`tree.hpp:1760-1769`); and `translateSource` bounds any declared
`categoryCounts` and `referenceCodes` (`C_interface.cpp:146-160`). The residual
unvalidated surface is the header-only engine used directly - `tests/cpp`, and
any future non-R host. That is a real but much narrower exposure than "the flat
C API", and section 3 is scoped to it.

**The column kind does not need to be in the serialized state.**
`SamplerStateData` carries chains, cut points, a sample cursor and a
recorded-draw count, and nothing else (`sampler.hpp:89-97`); neither `types` nor
`numCuts` is ever written. `docs/design/data-store.md:348-350` states why:
"Predictors are not serialized in state; `getPointer()` re-creation rebuilds the
store from the stored data object". Every restore entrance re-creates from a
data object first - R's `getPointer`/`setState` (`R/dbarts.R:1927-1948`,
`:1959-1973`), the flat C `dbarts_sampler_setState` (`C_interface.cpp:953-959`),
and stan4bart's own re-creation (`~/Repositories/stan4bart/src/init.cpp:373-376`)
- so the kind rides `dbartsData@varTypes` (`R/A_class.R:487`), an S4 slot R
serializes with the object, and reaches the rebuilt store through the ordinary
creation parse. Section 1 develops this, including the two bridge sites that
currently prevent the kind from surviving that parse.

## 1. The kind axis

### Alternatives

**A. A flag beside the existing pair** - `std::vector<std::uint8_t>
isOrderedFactor` alongside `types`. Memory +1 byte per column. Zero existing
branches change meaning: every `== ColumnType::categorical` test keeps its
sense, and new reads appear only where the kind matters. It is the cheapest
option, and it is the wrong shape: it creates a second per-column typing array
that can disagree with the first, which is precisely the failure mode
`cutPoints[j].empty()`-as-a-type-channel already demonstrates.

**B. A three-valued enum** - `enum class ColumnKind : std::uint8_t { numeric,
orderedFactor, categorical }` replacing `ColumnType`, with `numeric = 0` and
`categorical = 1` keeping their current values and `orderedFactor = 2`
appended. Memory zero: one array, same element width. Sites: exactly 45
comparisons (`== ColumnType::` / `!= ColumnType::`) across
`src/bartcore/*.hpp`, `R_interface_bartcore.cpp` and `C_interface.cpp`, plus 27
`columnIsPooled(...)` call sites that read the type indirectly. Almost all of
the comparisons are `== categorical` or `!= categorical` and translate
mechanically. The ones that need reading are the `else` arms that today mean
"ordinal" and would come to mean "numeric or orderedFactor" - which is the
correct meaning in every one inspected, because all of them are about threshold
splitting (`tree.hpp:648-658`, `:1573-1583`, `:2086-2088`;
`grow.hpp:267-296`; `sampler.hpp:903-909`, `:924-933`).

**C. A richer per-column descriptor** - one `ColumnSpec` struct absorbing
kind, category count, cut count, cut cap and missingness. ~16 bytes per column,
replacing four vectors; it touches all 45 comparison sites, all 62 `numCuts`
sites, and the `CodeBlock`/`ColumnSource` split that the data-store
consolidation already settled - `docs/design/data-store.md:362-365` records
"why there is no nested `CutGrid` type". This relitigates a landed decision for
no gain the other two do not give.

### The distinction that makes B cheap

The 45 comparison sites and the 27 `columnIsPooled` sites do not care about
semantics. They care about *mechanic*: does this column split by threshold or by
subset mask. Only a handful of sites - grid construction, ingestion validation,
reporting - care about *semantics*. So B should ship as two names, one derived
from the other:

  - `ColumnKind { numeric, orderedFactor, categorical }`, the semantic axis,
    one `std::uint8_t` per column, stored.
  - `splitsBySubset(j)` (equivalently, its complement `splitsByThreshold(j)`),
    the mechanic axis, derived, replacing every `types[j] == ColumnType::categorical`.

That reading also disposes of section 4 cleanly, and it means the 45 mechanical
sites become predicate calls rather than enum comparisons - which is what they
should have been.

**One admissibility decision this makes silently, and which must be stated
rather than absorbed.** Three refusal predicates key on `categorical`:
monotone directions (`facade.hpp:770-773`), leaf-covariate designation from a
predictor source (`facade.hpp:827-831`), and the same from a pre-built store
(`facade.hpp:867-868`). Under `splitsBySubset(j)` all three keep an
**orderedFactor column admissible**, exactly as today. For monotonicity that is
right and already ruled: `docs/design/monotone.md:117-121` says in as many
words that "ORDERED factors and numeric columns ... DO accept the constraint -
their codes are ordered, so `splitInterval` and the neighbor test are
meaningful." For **leaf covariates** it is a modelling choice rather than a
settled ruling: an ordered factor designated as a linear- or GP-leaf covariate
enters as its integer level codes standardized by
`standardizationMomentsForColumn` (`data.hpp:63-84`), i.e. treated as
equally-spaced interval data. That is defensible and is what happens today, but
today it happens because nothing can tell the column apart; after this change it
happens because we decided so. **Decision 2 keeps all three admissible,
deliberately** - the same reading that makes threshold splits right for the kind
makes its codes usable as an interval covariate.

### Killing the numCuts overload

**Accessor only**: `categoryCount(j)` returns `numCuts[j]`. Memory zero,
twenty-two call sites renamed, and the overload survives under a better name. It
does not meet the mandate.

**Real field**: `std::vector<std::uint32_t> categoryCounts`, with `numCuts[j]`
set to 0 for a categorical column (which is *true* - a categorical column has no
cuts). +4 bytes per column; on a 1000-predictor problem, 4 KB. The twenty-two
K-reads listed above change spelling; all are already type-guarded, so none
needs new logic. Three sites need more than a rename:

  - `data.hpp:555` (`columnIsPooled`), which must read the category count to
    decide the mask tier.
  - `R_interface_bartcore.cpp:1520-1521`, the verbose printer's "Number of
    cutoffs" line, which prints `store.numCuts[j]` for **every** column with no
    type guard and would start reporting 0 for categoricals. (The cut-point
    listing further down is already safe: `:1528` skips categoricals before
    `:1532` tests `numCuts[j] == 0`.) The printer needs a kind-aware line -
    reporting the category count for a categorical column is better output than
    it produces today, so this is an improvement to make deliberately rather
    than a regression to avoid.
  - `ColumnStore::buildFromParent` (`data.hpp:1413-1427`), which copies
    `types`, `cutPoints`, `numCuts` and `maxNumCuts` explicitly in both its
    column-subset and whole-store arms. `categoryCounts` **must** join that copy
    list. Two things break immediately otherwise: `refreshCategoricalTiers()` at
    `:1434` reads the count to set the pooled tier, and the gather loop at
    `:1464-1466` computes `missingCategoryCode(numCuts[j])` per column. A view
    with an empty or zeroed count vector gets the wrong mask tier and an
    out-of-range read, and every xbart fold view is built this way.

**Recommendation: the field.** The 4 bytes buy the property that `numCuts[j]`
means one thing everywhere, which is the whole point of the exercise, and it is
what a later per-column width experiment wants to read.

### The cut grid for an ordered factor (settled: this lands)

An `orderedFactor` column takes **K-1 cut points at consecutive level
midpoints**, replacing the `n.cuts` uniform cuts over `[0, K-1]` it gets today.

The midpoint placement itself is what quantile mode already produces: K unique
values under a cap of **at least K-1** induce exactly K-1 cuts
(`data.hpp:749-750`) at consecutive midpoints (`:789-794`). So
`buildCutsForColumn` (`data.hpp:867-885`) gains a third arm that runs the
quantile path for an `orderedFactor` column regardless of the store's
`useQuantiles` flag.

**That is not sufficient on its own, and the cap is the whole difficulty.**
`finishQuantileGrid` has a third arm (`data.hpp:752-756`) that fires when
`numUnique > maxNumCuts[j] + 1` and **thins** the grid to `maxNumCuts[j]`
cuts. `maxNumCuts[j]` is `n.cuts`, default 100. So routing a 150-level ordered
factor through the quantile path unchanged gives 100 cuts at `k + 0.5` for
k = 0..99, which sends every level from 100 to 149 to code 100: the levels
still merge. The collision would change shape - from spread across the range to
concentrated in the tail - rather than disappear. Any spec that does not say
what happens to the cap does not fix the defect.

**The cap treatment: `maxNumCuts[j]` is raised to K-1 for an `orderedFactor`
column, and the column is refused above the code type's ceiling.** Concretely,
in the `orderedFactor` arm of `buildCutsForColumn`, before the grid is built:
`if (maxNumCuts[j] < K - 1) maxNumCuts[j] = K - 1;`. The in-tree precedent for
mutating the cap to fit an externally determined grid is
`setCutPointsForColumn`, which does exactly `if (maxNumCuts[j] < numCutPoints)
maxNumCuts[j] = numCutPoints;` (`data.hpp:933`). Raising the cap rather than
bypassing the thinning arm is the better of the two for one reason: it keeps
`numCuts[j] <= maxNumCuts[j]`. A bypass would leave a column carrying K-1 cuts
under a recorded cap of `n.cuts` - a cut count exceeding the column's own stated
ceiling, which is the invariant `build`'s clamp (`data.hpp:1155-1157`)
establishes and which `setCutPointsForColumn` preserves whenever an externally
chosen grid is wider than the cap.

**The cap raise is necessary but not sufficient: the mutation path must
dispatch on kind too.** Neither `refreshCutsForColumn` (`data.hpp:907-918`) nor
`cutsWouldRemainValid` (`:923-931`) reads `maxNumCuts[j]` at all. Both branch on
the **store-wide** `useQuantiles` flag after a single `categorical` early
return, so an `orderedFactor` column falls straight through to the ordinal arms.
On a default `useQuantiles = FALSE` store an `updateCuts` mutation therefore
reaches `fillCutsUniformly` (`:915`) and replaces the midpoint grid with K-1
uniform cuts over the **observed** range. Three things follow, and the third is
the serious one:

  - The cut *values* change - midpoints at `k + 0.5` become
    `(k+1)(max-min)/K` - and the flat tree format matches a saved cut value
    **exactly** (`cuts[k] != flat.value` refuses, `tree.hpp:1578-1580`), so a
    state saved before such a mutation no longer restores onto the column.
  - On a `useQuantiles = TRUE` store the other arm rebuilds from observed
    uniques (`:927`), silently degrading a declared-level grid to an
    observed-level one - the sub-decision above, undone by a mutation.
  - When the mutation's observed range does not span every declared level, the
    uniform grid's codes stop being level indices and levels outside that range
    share a code again. **The defect S2 exists to remove would be reintroduced
    by the first `updateCuts` mutation** - worse than leaving it in place,
    because it would then depend on call history.

The fix is the same dispatch in the same place: an `orderedFactor` column joins
`categorical` in the early arm of both functions. For both kinds the grid is a
property of the level table rather than of the values, so a refresh is a no-op
that keeps the existing grid, and the feasibility check is level membership -
the categorical arm's `categoricalValueIsValid` sweep (`:941-945`) - rather than
an induced-cut-count comparison. So S2 is three small edits, not one, and
omitting either of the latter two would make the grid fix conditional on the
sampler never being mutated.

The ceiling is the code type's, not `n.cuts`: `build` clamps `maxNumCuts[j]` to
`maxNumCutsRepresentable = 0xFFFD` = 65533 (`data.hpp:35`, applied at
`:1155-1157`), so a grid of K-1 cuts is representable exactly while
K <= 65534. **An ordered factor with more than 65534 levels is refused at
ingestion, naming the cap** - the ordered-side counterpart of the unordered
factor's refusal above 65535 levels (`R/utility.R:500-507`), which is where the
R-side check belongs so the message names the column. The two ceilings are one
apart (65534 against 65535) because the ordered side spends a code on the
grid's upper bin and the unordered side does not; the two ceilings are
stated as the numbers they are rather than claimed to cap identically.

**The cost this accepts.** A K-level ordered factor now carries K-1 cut points:
`8 * (K - 1)` bytes of `cutPoints[j]`, `2 * (K - 1) * 8` bytes of grow-from-root
scratch (`grow.hpp:276`), and an O(K) scoring tail per (node, variable) in
`scanOrdinalCuts`. At the ceiling that is roughly 524 KB of cuts and 1 MB of
scratch for one column, and a scan tail that will dominate deep-node cost. This
is expensive, and it is the honest expense: a K-level ordered factor genuinely
admits K-1 threshold splits, so a grid smaller than that is a silently
restricted model rather than a cheaper one. A user who wants 100 cuts over a
high-cardinality ordinal encoding should pass the column as numeric, where
`n.cuts` means what it says. The comparison worth noting is that a
high-cardinality *unordered* factor is already expensive in the same way -
pooled masks cost `maskWordsForCount(K)` words per rule (`data.hpp:45-47`) and
the category histogram is O(K) per node (`scan.hpp:301-304`) - so this does not
introduce a cliff the store did not already have for wide factors.

The rejected alternative is to keep thinning and warn. It preserves today's
performance for wide ordered factors at the cost of keeping a silent
model restriction behind a warning, which is the defect with a message attached
rather than a fix.

**One sub-decision inside that: observed levels or declared levels?**
`finishQuantileGrid` builds its grid from `sortedUnique`, the **observed**
values only (`data.hpp:735-740`). So a level that appears in the training data
zero times contributes no cut, its neighbours merge, and the column's codes
become ranks among observed levels rather than level indices. Two consequences:
a level absent from training but present in test data cannot be separated from
its neighbour, and the code-to-level correspondence is not stable under
subsetting. The alternative is to place K-1 cuts from the **declared** level
count, which the factor's own level table already carries - `factor.levels`
(`R/utility.R:515`, filled for ordered factors too). That table does not reach
the store today for exactly the reason this proposal exists:
`readDeclaredCategoryCounts` skips any column the bridge does not consider
categorical (`R_interface_bartcore.cpp:896`), so an ordered factor's declared
count is dropped. The kind axis is what makes it available, and using it is a
small extension of that function rather than new machinery.

That function has a **second** early-continue worth naming: `:897` skips any
column whose source is CSC-backed (`sourceOf(j) < 0`), because a container
declares its own K for those. An ordered factor cannot reach that path from R -
`makeCategoricalModelMatrix` only ever produces a `sparseFactor` as a
*categorical* column (`R/utility.R:491-493`) - but the flat C API and
`tests/cpp` can construct one, and under a declared-levels grid such a column
would fall through to an inferred count. Whichever way the sub-decision goes,
the CSC-backed ordered factor needs a stated answer rather than an inherited
`continue`.

**Decided: declared levels** (decision 4), because it makes the grid a property
of the factor rather than of the sample, which is what the level table is for.
`readDeclaredCategoryCounts` is extended past **both** gates, so a CSC-backed
ordered factor takes its container's declared K rather than an inferred count.
This also fixes what `K` means in the cap rule above - declared count, not
observed-unique count - so the two are one decision.

**Two residues of that ruling, stated rather than inherited.** Extending the
type gate is enough for every column R can build; the `sourceOf(j) < 0` gate
could not be dropped literally, because `resolveCscCategoricalReferences`
resolves a CSC-backed CATEGORICAL column's K and its reference code together
from the container's own metadata, and letting a `factor.levels` entry
override that would contradict "a container declares its own K". S1 therefore
substitutes an equivalent rule - a count the container already declared stays
authoritative - which selects the same set for categoricals (that resolution
zero-fills, then writes a strictly positive count for exactly the CSC-backed
categorical columns) and lets a CSC-backed ordered factor take its
`factor.levels` entry. What is left open, and what S2 must rule when the
declared count becomes the grid:

  - **Which channel wins for a hand-built CSC-backed ordered factor** that
    declares both a `factor.levels` entry and a `sparseCategoryCount` for the
    same column. R cannot construct one - `makeCategoricalModelMatrix` types
    every `sparseFactor` categorical whether or not it is ordered
    (`R/utility.R:491-493`) - so for anything R produces the two are the same
    number; a hand-built container or a header-only host can disagree.
  - **Such a column's implicit rows read the quantized zero, not its declared
    reference level.** `resolveCscCategoricalReferences` skips a
    non-categorical column, so its `refCode` is never resolved and
    `quantizeCscColumnInto` takes `codeFor(j, 0.0)` for the implicit rows.
    That is the pre-kind ordinal behavior and is unreachable from R, but it is
    a silent wrong answer for a host that declares one, and S2's
    declared-levels grid is what makes the code positions matter.

With the cap treatment specified, three consequences follow:

  - **The level-collision defect disappears**, for K up to the ceiling. With
    `maxNumCuts[j]` raised to K-1 the thinning arm never fires, every level
    boundary gets its own cut, and the pigeonhole failure has nowhere to occur.
    Above the ceiling the column is refused rather than silently merged, which
    is the point: no configuration ends in a quiet collapse.
  - **The prior-multiplicity tilt disappears.** Each distinct threshold gets
    exactly one index slot, so `logCut` (`grow.hpp:270`) normalizes over the
    real partition set and the induced prior over level boundaries is uniform.
  - **`n.cuts` no longer applies to an ordered factor** - it is overridden, not
    merely inapplicable, since `maxNumCuts[j]` is actively raised past it. That
    is the consistent completion rather than a new exception: `n.cuts` already
    does not apply to unordered categoricals, where `numCuts[j]` holds K and no
    cut grid exists at all (`data.hpp:868-873`). A factor's grid is determined
    by its levels; a numeric column's by `n.cuts`.

**What this does to the "never cap" test invariant.**
`inst/tinytest/test-data-categorical.R:40-41` carries the comment "ordered
factors are ordinal and never cap" over a 54-level ordered-factor case at
`:47-51`. The **assertion survives**, but not for the arithmetic reason: it
survives because `dbartsData` is pure R construction and never reaches
`buildCutsForColumn` at all, so no grid is built and no cap is consulted at the
point the test measures. What the new refusal changes is the *R-side* level
check that would sit beside `R/utility.R:500-507`, and 54 levels is far below
the 65534 ceiling, so that case passes too. The comment should be rewritten to
say what is now true - an ordered factor caps at 65534 levels rather than never
- and the file wants two additions: a case above `n.cuts` (say 150 levels)
asserting that adjacent levels remain separable in a **fitted** sampler, which
is the regression the old grid could not catch and which `dbartsData` alone
cannot express, and a case above the ceiling asserting the new refusal.

### Where the kind lives in the serialized state

**It does not, and it should not.**

The state's one bit of per-column typing is carried by absence: the writer skips
a categorical column's cut vector, leaving `R_NilValue`
(`R_interface_bartcore.cpp:6691-6692`), and both readers test emptiness against
the live store's type (`sampler.hpp:899-910` in `setState`, `:993-999` in
`installForests`). Under `ColumnKind` that bit still discriminates exactly what
it must - subset-mask rules versus threshold rules - because an `orderedFactor`
column keeps a non-empty cut grid, exactly as a numeric one does. Nothing about
the encoding changes, and the grid fix does not change it either: a state
carries its own cut vector and `setState` installs it (`sampler.hpp:929-932`),
so a state written under either grid restores onto its own cuts.

Therefore neither `stateFormatVersion` (`R_interface_bartcore.cpp:6437`) nor
`minReadableStateFormatVersion` (`:6446`) moves.

The silent-misread hazard is real but applies to a different move. The registry
rule at `R_interface_bartcore.cpp:6411-6436` (restated at
`docs/design/public-surface.md:163-185`) says a *new named* block or top-level
attribute is additive - an old reader ignores it, a new reader defaults it - and
bumps nothing. What is non-additive, and what bumps both constants, is
re-encoding an *existing* slot. Writing a non-null `cutPoints` entry for a
categorical column would be exactly that: an old reader would see a present cut
vector for a categorical column and refuse at `sampler.hpp:902`, and a new
reader looking at an old state could not distinguish "categorical" from "some
new kind with no cut grid". That is the same silent-rename hazard the
`bcf` -> `glue` bump was for (`docs/design/multiplier-combiner.md:44-48`). So:
do not touch the `cutPoints` encoding.

If a future capability genuinely needs the kind in state, the additive route is
a new top-level `columnKinds` attribute, INTSXP of length `numPredictors`, read
by name and defaulted when absent, with the current emptiness test kept as the
fallback. That would bump nothing either. `docs/design/data-ownership.md:226-230`
already anticipated this and recorded that the cutPoints-only encoding was kept
"on simplicity". Keep it.

### The changes that carry the kind from R to the store

The kind must survive four hops, and two of them currently drop it. Both are in
the bridge, not in R, and both are load-bearing enough that omitting either
makes the feature silently wrong rather than merely absent.

  - `R/data.R:1-2` gains a third constant.
  - `R/A_class.R:579-583` currently refuses any `varTypes` value that is not 0
    or 1; it must admit the third.
  - **`R_interface_bartcore.cpp:1060-1065` parses `varTypes` as "zero is
    ordinal, anything else is categorical"** (`if (variableTypes[j] == 0)
    continue; data.columnTypes[j] = categorical;`). Admitting `2L` on the R side
    without rewriting this loop makes every ordered factor enter the store as
    **categorical** - subset splits over identity codes - which is the exact
    opposite of the settled ruling at `public-surface.md:123-126`, and it is a
    silent draw change rather than an error. This is the single highest-risk
    edit in the proposal.
  - **`R_interface_bartcore.cpp:1066-1067` publishes the typing channel only
    when something is categorical**: `data.predictors.columnTypes =
    data.anyCategorical ? data.columnTypes.data() : NULL`. A data frame with an
    ordered factor and no unordered factor sets `anyCategorical` false, so the
    channel is null, `build` falls back to all-ordinal (`data.hpp:1142-1146`),
    and the kind dies at the store boundary. That is the *headline* case for
    this feature - an ordered factor on its own - so the gate must become "any
    column is not numeric", not "any column is categorical".
  - `inst/tinytest/test-data-categorical.R:25` asserts
    `varTypes == c(1L, 0L, 0L, 0L, 1L)` where the second column is the ordered
    factor `o` (`:8-12`); it becomes `c(1L, 2L, 0L, 0L, 1L)`.

**Provenance divergence, and why it cannot be repaired.** A `dbartsData` saved
by an older build carries `varTypes` 0 for an ordered factor, and 0 is ambiguous
- it also means numeric. There is no upgrade path: nothing in the stored object
distinguishes the two. (`factor.levels` is *nearly* a marker - `R/utility.R:515`
fills it for ordered factors too - but it is an attribute of `@x`, not a slot,
and it is equally present for unordered factors, so it identifies "was a factor"
rather than "was an ordered factor".) With the grid fix landing, the consequence
is concrete: **the same data frame yields a different grid, and therefore a
different fit, depending on whether its `dbartsData` was built before or after
this change.** Old objects keep working and keep behaving exactly as they do
today; they simply do not get the new grid. That is acceptable pre-release, but
it is observable by `load()`ing an old object beside a freshly built one, so it
gets one UPGRADING bullet in the scheduled NEWS rewrite (decision 7).

### Recommendation

Alternative B: `ColumnKind {numeric, orderedFactor, categorical}` with
`numeric = 0`, `categorical = 1`, `orderedFactor = 2`; a derived
`splitsBySubset(j)` predicate carrying the 45 mechanical sites; a real
`categoryCounts` field with `numCuts[j] = 0` for categoricals, copied in
`buildFromParent`; the kind-aware midpoint grid; no state channel and no version
movement; the two bridge sites above rewritten in the same commit as the R-side
widening, never after it.

**On the ABI, and on what an enumerator costs.** `DBARTS_COLUMN_TYPE_LIST`
(`dbarts.h:286-288`) is a list macro built for exactly this - its comment at
`:281-285` explains that enumerators are added to the list rather than the enum
body so the compile-time token cannot miss them. Appending
`DBARTS_COLUMN_ORDERED_FACTOR = 2` keeps both existing values, and forces
`C_interface.cpp:146-151`, which rejects any `columnTypes[j]` that is neither
`DBARTS_COLUMN_ORDINAL` nor `DBARTS_COLUMN_CATEGORICAL` by name in its error
text, to be extended. The enumerator list is folded into `DBARTS_C_API_HASH`
(`C_interface.cpp:521`, alongside the struct layouts at `:486-509`), so the
append re-bakes the token, and stan4bart compiles with
`-DDBARTS_REQUIRE_EXACT_ABI` (`~/Repositories/stan4bart/src/Makevars:3`), which
turns token movement into a hard failure at first stub resolution. That is a
lockstep stan4bart rebuild, and it is an accepted cost.

It costs two re-bakes across the plan, not one: the enumerator lands in S1 and
section 2's struct append in S4a. **Two is accepted.** The alternative - S1
carrying the struct fields too, appended and documented as reserved with
`translateSource` refusing a non-null value until S4a reads them - spends the
re-bake once but publishes a field shape before its consumer exists, which is
the larger commitment. Inside the pre-release window, where the
`-DDBARTS_REQUIRE_EXACT_ABI` guard is documented as temporary and the
coordinated rebuild is already scheduled, a second rebuild is cheaper than a
guessed shape. Revisit only if the lockstep rebuild proves expensive in practice
rather than in principle.

## 2. Native integer ingestion, and the typed-storage end state

This section is the design for the full rework, not a comparison against doing
less. The narrow alternative is stated once, at the end, for the record.

### The end state, if the store were designed today

The organizing principle is one sentence: **a column's stored representation is
the one its values actually have.** A factor column's values are level indices,
so its storage is integer codes; a numeric column's values are real numbers, so
its storage is doubles. Nothing is stored twice, and no channel carries a type
it is pretending to be.

Concretely, for a factor column - `categorical` or `orderedFactor` - the end
state is:

  - **No double representation anywhere in the store.** The authoritative
    per-cell datum is the quantized code already in `train.codes` /
    `test.codes` (`xint_t`, 2 bytes, `data.hpp:21`). For a `categorical` column
    that code *is* the level index by construction, unconditionally: `codeFor`
    is the identity cast (`data.hpp:709-711`). For an `orderedFactor` column
    under the midpoint grid the code equals the level index only when the grid
    is built from **declared** levels, or from observed ones with every level
    present - K-1 cuts at `k + 0.5` send level `v` to code `v` under
    `lower_bound` (`data.hpp:713-716`) - and it degrades to a rank among
    observed levels otherwise, per the sub-decision in section 1.

    **The rework does not depend on that identity.** What it depends on is the
    weaker and already-established fact that *no consumer needs a double back
    from a factor column*: routing reads codes, and the leaf-covariate loads -
    the only readers of `rawColumn`/`rawTestColumn` - never see one (see two
    bullets down). The identity matters only to a hypothetical code-to-double
    inverse, which the end state deliberately does not provide. Worth one
    `tests/cpp` assertion anyway, because it is the property a reader will
    assume when they encounter the storage.
  - **`ColumnSource` carries a typed source, not a `double*`.** Today
    `denseRaw` is a bare `double*` (`data.hpp:148-150`). In the end state the
    descriptor names which pool the column's raw lives in - the double pool, the
    code pool, or neither - so that "this column has no double" is representable
    rather than being a null a reader must know to expect.
  - **`ownedDenseValues` shrinks to the real-valued columns.** Today it is a
    single contiguous double block sized by the largest dense *source* index and
    copied wholesale from the caller (`data.hpp:1170-1181`), covering factor and
    numeric columns alike. In the end state it holds only columns whose values
    are genuinely real-valued, addressed through a per-column offset table
    rather than by source index. The same applies to `ownedTestValues`
    (`data.hpp:1314-1315`), which today owns every dense test column as doubles
    including the factors.
  - **`rawColumn` / `rawTestColumn` return null for factor columns, by
    design rather than by accident.** This is already safe: their only consumers
    are the linear- and GP-leaf covariate gathers (`model.hpp:1017`, `:1039`,
    `:1059`, `:1384`, `:1409`, `:1423`), the view gather
    (`data.hpp:1443`), and the designation validity check (`facade.hpp:868`) -
    and a categorical column is refused as a leaf covariate at the factory
    (`facade.hpp:827-831`, `:867-868`). The one live question the grid fix
    raises is settled by decision 2: an `orderedFactor` **stays admissible** as
    a leaf covariate, making it the single column kind that needs a double for
    a reason other than splitting. It is served the way every other leaf
    covariate is - gathered into `gatheredRawValues` (`data.hpp:596`) by
    designation - and never from a retained pool, so the no-factor-pool shape
    below holds without exception.
  - **`rawColumnForRequantize` returns null for factor columns.** Also already
    correct: a factor column never re-quantizes. `codeFor` for a categorical is
    grid-independent, `setCutPointsForColumn` skips categoricals
    (`sampler.hpp:924-925`), `buildCutsForColumn` clears their cuts
    (`data.hpp:873`), and `refreshCutsForColumn` returns immediately
    (`data.hpp:908`). Under the midpoint grid an `orderedFactor`'s cuts are a
    function of its level count alone, so it does not re-quantize either.
  - **The double-typed ingestion paths gain code-typed arms**, and the arms are
    selected by kind, not by a nullness accident:
      * `build`'s per-column loop selects `const double* column` and hands it to
        both `buildCutsForColumn` and `quantizeColumn` (`data.hpp:1229-1236`).
        This becomes a kind-dispatched pair.
      * `inferredCategoryCount` scans doubles for a max (`data.hpp:842-847`);
        its code-typed sibling scans `int32_t` skipping `NA_INTEGER`.
      * `quantizeDenseObserved` takes `const double* raw` and detects
        missingness with `isNA` on a NaN (`data.hpp:953-965`); its code-typed
        sibling narrows `int32_t` to `xint_t` and detects `NA_INTEGER`.
      * `buildCutsForColumn` (`data.hpp:867-885`) needs no double at all for a
        factor column: the categorical arm needs only the observed maximum and
        the declared count, and the `orderedFactor` arm needs only the level
        count.
      * `writeOwnedDenseColumn` / `writeOwnedDenseCell` (`data.hpp:1506-1521`),
        the mutation write-through, become kind-aware: a factor column's
        mutation writes codes, or writes nothing at all if the codes are the
        only representation.
  - **CSC columns are the one place doubles survive for a factor.** A
    `dgCMatrix`'s `@x` slot is double by Matrix's own definition, and the train
    side *borrows* it (`data.hpp:1211-1212`) rather than copying - the borrow is
    the point of the CSC path. So a CSC-backed factor column keeps a double
    slice on the train side, and the code channel does not reach it. The
    **owned** copies can be typed: `ownedTestCscValues`
    (`data.hpp:1331`, `:1350-1351`) and the post-mutation
    `ownedCscValues` (`data.hpp:620`) are the store's own buffers, and a factor
    column's could hold `int32_t`. Whether that is worth a second typed CSC pool
    is a real question and the honest answer is probably not - see the byte
    table.

### The predict and replay story

This is the half of the design where "the factor column has no double" has to
survive contact with a consumer that only knows how to read doubles.

**What actually reads doubles per row, and what does not.** The training-side
hot loops read **codes only** and are untouched by anything in this section:
`partitionChildren` (`tree.hpp:866-926`) reads `data.column(variable)`, a
`const xint_t*`, on every arm (`:888`, `:912`) or the `SparseColumnData` rank
storage (`:879`, `:900`); `scanOrdinalCuts` and the categorical histogram scan
read the same codes (`scan.hpp:310`). Recorded test fits during `run` route by
the `test` code block (`data.hpp:575`), also codes. So the doubles this section
removes are not on any sampling hot path - they are on the **replay** path only,
which serves `dbarts_sampler_predict` and saved-tree replay against a
caller-supplied matrix.

**The reader seam is already the right shape.** The flat replay takes a
`Columns` template parameter whose reader exposes `double at(size_t row)`
(`tree.hpp:1750-1754`), and there are two implementations:
`DenseColumns`/`DenseColumnReader` for a plain block (`tree.hpp:1724-1741`) and
`PredictorSourceColumns`/`PredictorSourceColumnReader` for a borrowed
`PredictorSource`, which already dispatches per column between a dense double
pointer and a sparse rank lookup (`data.hpp:382-389`, built at `:400-438`). Its
comment states the intent - "a source that maps store columns onto other storage
substitutes here without a second code path." An int-backed source adds a third
arm to that reader. `partitionFlatIndices` does not change, and there is no
second categorical loop.

**Two distinct sources, two distinct NA conventions.** A factor column's double
is not simply a view over its codes: the two sources carry different NA
conventions, and the distinction has to be explicit:

  1. **A caller's integer source** (what an int-backed `PredictorSource`
     carries) uses R's `NA_INTEGER`. Its reader maps `NA_INTEGER -> NA_REAL` and
     widens everything else straight across. The replay's own arms then behave
     exactly as they do today: `isNA(value)` routes the missing direction
     (`tree.hpp:1766`, `:1781`, `:1796`) and the categorical arms cast to an
     index (`:1768-1769`, `:1783-1784`). No knowledge of the column is needed,
     and the mapping is total.
  2. **The store's own codes** use the reserved missing code, which is *not*
     `NA_INTEGER` and *not* a level index: `codeFor` sends NA to
     `missingCategoryCode(numCuts[j])` (`data.hpp:710`), which is `K` for a
     pooled column (K >= 64) and the fixed position `63` otherwise
     (`data.hpp:50-54`). The inverse is therefore
     `code == missingCategoryCode(K) -> NA_REAL`, everything else
     `-> (double)code`, and **it requires the column's K** to know which code is
     the sentinel. It is unambiguous - an inline column has codes 0..K-1 <= 62
     with 63 free, a pooled one has 0..K-1 with K free (`data.hpp:37-41`) - but
     it is not a cast.

     The design consequence: **the store's codes are never served to a
     double-typed reader.** Nothing in the end state requires it - the replay
     always reads a caller's source, and the leaf-covariate loads never see a
     factor column - so the inverse in (2) exists as a documented property of
     the encoding rather than as a code path. If a future feature does need it
     (predicting on the training rows without re-supplying the matrix, say), it
     is a `rawFromCodes(j, i)` accessor on the store, which has K, and not
     something a bare pointer can do. That asymmetry is a reason to keep the two
     apart, not a detail to paper over.

**The test-side store.** `buildTest` owns every dense test column as doubles
today (`data.hpp:1314-1315`) and points each column's `denseRaw` into that block
(`:1339-1341`), then quantizes into `test.codes`. In the end state a factor test
column keeps only its codes; the doubles it owned served re-quantization (which
a factor never needs) and `rawTestColumn` for leaf covariates (which a
categorical cannot be). That is where the test-side saving in the table below
comes from.

### Byte accounting

Per column, `n` train rows and `m` test rows; `nnz` the nonzeros of a CSC
column. "Retained" means held for the sampler's lifetime; "transient" means
allocated for the duration of one creation or mutation call.

**Train side.**

| build path | column kind | today (retained) | end state (retained) | delta |
| --- | --- | --- | --- | --- |
| dense (unmapped) | factor | 2n codes | 2n codes | 0 |
| dense (unmapped) | numeric | 2n codes | 2n codes | 0 |
| mapped (mixed / all-dense container) | factor | 2n + 8n owned double | 2n | **-8n** |
| mapped | numeric | 2n + 8n owned double | 2n + 8n | 0 |
| CSC-backed | factor | 2n codes (+ borrowed slice) | 2n (+ borrowed slice) | 0 |
| CSC-backed, after first mutation | factor | 12*nnz owned (8 value + 4 row) | 8*nnz (4 code + 4 row) | -4*nnz |

The dense (unmapped) row is 0 because the engine already retains no raw there
(`build`'s unmapped arm packs code offsets and moves on, `data.hpp:1193-1198`);
a later re-quantize reads whatever matrix the caller supplies
(`rawColumnForRequantize`, `data.hpp:655-660`). **This is the common
factor-data-frame path** (`R_interface_bartcore.cpp:986-1008`), so the retained
saving does not apply to it. Its saving is transient instead.

**Test side.**

| build path | column kind | today (retained) | end state (retained) | delta |
| --- | --- | --- | --- | --- |
| dense or mapped `buildTest` | factor | 2m codes + 8m owned double | 2m | **-8m** |
| dense or mapped `buildTest` | numeric | 2m + 8m | 2m + 8m | 0 |
| CSC-backed test | factor | 2m + 12*nnz owned | 2m + 8*nnz | -4*nnz |

**Transient, at creation.** The bridge assembles a `p x n` double block via
`codeDenseColumn` (`R_interface_bartcore.cpp:543-555`, `:995-1008`). Per factor
column that is 8n today and 4n with a code channel: **-4n transient per factor
column**, plus one avoided widen-then-narrow per cell.

**Views** (`buildFromParent`) already hold no factor raw: they gather codes for
routing (`data.hpp:1431-1434`, `:1464-1471`) and raw only for designated leaf
covariates (`:1441-1449`), which cannot be categorical. No change.

**Worked total.** A mixed container with `n = 1e6` train rows, `m = 1e5` test
rows and 20 dense-backed factor columns retains
`20 * 8 * 1e6 = 160 MB` of train doubles and `20 * 8 * 1e5 = 16 MB` of test
doubles today, all of it duplicating information already in the codes. The end
state retains none of it: **176 MB**, against a store whose codes for those
columns total `20 * 2 * 1.1e6 = 44 MB`. The direction is a four-to-one shrink on
the affected columns, and nothing anywhere widens.

The honest qualification: **the common all-dense factor-data-frame path retains
nothing today**, so it gets the transient saving and no retained saving. The
retained numbers above are the mixed/CSC and test-side paths.

### Hot-path impact, and where bitwise preservation could silently fail

**Loops that do not change type at all** - byte-identical machine code, no
verification burden beyond recompiling:

  - `partitionChildren`, every arm (`tree.hpp:866-926`): codes and rank
    storage only.
  - `scanOrdinalCuts` and the categorical histogram (`scan.hpp:301-310`): codes
    only. (`scan.hpp:301-304` reads the category count, so it changes *spelling*
    under section 1's `categoryCounts` field, not type.)
  - The leaf-covariate loads (`model.hpp:1017`, `:1039`, `:1059`, `:1384`,
    `:1409`, `:1423`): doubles, and never a factor column.
  - Recorded test fits during `run`: `test.codes`, unchanged.

**The one loop that changes type**: `PredictorSourceColumnReader::at`
(`data.hpp:382-389`), on the replay path only. It already branches per row
(`dense != nullptr ? dense[row] : sparse->at(row)`); an int arm makes that a
two-test chain. The branch is perfectly predicted within a column - the reader
is fetched once per node outside the row loop (`tree.hpp:1754`) even though
`at()` is called inside it - so the cost is a predicted test per row on a path
that already does one. It is not free and it is not measurable against the
`isNA` test and the mask lookup in the same loop body. This is the only place a
`bench-sampler` predict arm could move, and it is worth a measurement rather
than an assumption.

**Where a pure representation change could silently stop being bitwise.** All
four of these produce correct-looking output and wrong draws:

  1. **`hasMissing[j]` flipping.** `quantizeDenseObserved` sets `anyMissing`
     from `isNA(raw[i])` (`data.hpp:960`). A code-typed arm must set it from
     `raw[i] == NA_INTEGER`. If it does not, a column with missing values goes
     from flagged to unflagged, and the invariant at `data.hpp:542-544` is
     explicit that an unflagged column "consumes no missing-direction draw from
     the rng" - so **every subsequent draw in the chain shifts**. This is the
     single most dangerous failure in the rework, and it is silent unless the
     gate corpus contains an NA-bearing factor column.
  2. **The NA code itself.** `NA_INTEGER` is `INT_MIN`; narrowed without a
     guard it becomes some `xint_t` in range, which is a legal-looking category.
     The forward map must send it to `missingCategoryCode(K)`, matching what
     `NA_REAL` produces today (`data.hpp:710`).
  3. **`inferredCategoryCount` on the code path.** The double version skips NaN
     and takes a max (`data.hpp:842-847`). A signed int version skipping
     `NA_INTEGER` is correct; an unsigned one would let `INT_MIN` become a huge
     maximum and inflate K, which changes the mask tier
     (`maskWordsForCount`, `data.hpp:45-47`) and therefore the rule encoding.
  4. **The CSC implicit-row rule.** A CSC categorical column's absent rows read
     the reference code (`data.hpp:985`, mirrored in
     `materializePredictorSource` at `:319-323` and in the replay reader at
     `data.hpp:425-428`). Any typed variant must apply the same rule in all
     three places or the implicit rows diverge between the store and the replay.

The mitigation for all four is the same and is cheap: the gate corpus needs an
NA-bearing factor column and a CSC-backed factor column, both of which section 6
already requires for other reasons.

### Alternatives within the end state

**int32 vs. int16 codes at the ABI and in `PredictorSource`.** `xint_t` is
`uint16_t` and `maxCategories` is `0xFFFF` (`data.hpp:21`, `:95`), so the store
narrows to 16 bits regardless. Taking `int16_t` at the boundary would halve the
transient block again. Against it: R integer vectors are int32 and `NA_INTEGER`
is `INT_MIN`, so an int16 channel has no representation for R's NA and would
force the bridge to do a narrowing pass R would not otherwise do - reintroducing
exactly the per-cell conversion the change removes. It would also couple the ABI
to the store's internal code width, which the backlogged per-column width
experiment (`docs/plans/archive/hot-layer-u8.md`) may make per-column. **int32
at the boundary, narrowed once at ingestion**, with `int32_t` already the
precedent at `dbarts.h:340`, `:346`.

**One typed pool per kind vs. per-column buffers vs. no factor pool.** Three
shapes for the retained raw:
  * *Two pools* - a double pool for real-valued columns and an int32 pool for
    factor columns, each with a per-column offset table. Keeps the current
    contiguous-block character; costs an offset table and a second pool that
    duplicates the codes.
  * *Per-column buffers* - `std::vector<std::vector<...>>`. Simplest to reason
    about; loses contiguity and adds an allocation per column.
  * *No factor pool at all* - factor columns retain nothing beyond their codes.
**Decided: no factor pool** (decision 3). It is what makes the byte table's
`-8n` real: a second int32 pool would only turn `-8n` into `-4n` while keeping a
duplicate representation, which is the thing the section exists to remove. The
double pool stays contiguous, addressed through a per-column offset table
instead of the source index, holding real-valued columns only. The accepted cost
is recorded with the decision: a factor column's raw width now follows the code
width, narrowing what a later per-column width experiment can vary
independently.

### Slice decomposition and sizing

Seven slices. S0-S3 are section 1, 3 and the gate work; S4a-c are this section.

  - **S0 - gate enabling.** Add an ordered-factor predictor case, an NA-bearing
    factor column, and a CSC-backed factor column to
    `benchmarks/R/equivalence.R`; re-record. No product code. Small. Must go
    first: it is what makes every later slice observable.
  - **S1 - kind axis.** `ColumnKind`, `splitsBySubset`, `categoryCounts`
    (including the `buildFromParent` copy), the two bridge sites, the printer,
    the ABI enumerator, the `FlatKind` predicate collapse, the
    `denseCallSupplied`/`denseResident` rename. Draw-preserving. Independently
    landable and gateable. Medium: 45 comparison sites, 22 count reads, 21
    rename references, plus R and bridge.
  - **S2 - the midpoint grid.** Draw-CHANGING; own re-record. Modest in code -
    a third arm in `buildCutsForColumn` reusing the quantile path, the
    `maxNumCuts[j]` raise, matching `orderedFactor` arms in
    `refreshCutsForColumn` and `cutsWouldRemainValid` so a mutation cannot undo
    the grid, and the R-side level refusal - and large in gate. Independently
    landable, and must follow S1 because it needs the kind. The mutation arms
    are not optional: without them the fix survives only until the first
    `updateCuts` call, so a mutation test over an ordered factor belongs in this
    slice.
  - **S3 - engine-side ingestion validation** (section 3), including the
    engine-side backstop for S2's level ceiling, since that refusal is
    otherwise R-side only. Draw-preserving.
    Small. Independent of S2; pairs naturally with S4a.
  - **S4a - the code channel, end to end through the bridge's validation.**
    `denseCodes` on `PredictorSource` and `dbarts_predictor_source`; the bridge
    fills it from an `INTSXP` instead of widening; code-typed
    `inferredCategoryCount` and `quantizeDenseInto`; `build`'s unmapped arm
    consumes it. **And** code-typed arms for the bridge's own pre-store
    validation sweep, which is double-typed throughout and reads the very block
    the slice is trying not to build: `rawViewColumn`
    (`R_interface_bartcore.cpp:1625-1629`) returns a `const double*` off
    `view.denseValues`, and `trainingCategoryBound` (`:1680-1690`),
    `validateCategoricalPredictors` (`:1702-1725`) and
    `validateTestContainerAgainstStore` (`:3024-3041`) all scan through it -
    as does `materializePredictorSource` (`data.hpp:304-332`) for the non-dense
    views. Draw-preserving subject to the four hazards above. **Medium-large**,
    once those arms are counted: the channel itself is small and the validation
    sweep is most of the work.

    The scoping matters because the two halves are not separable. The transient
    block *is* what the sweep reads, so a version of S4a that adds the channel
    and leaves the sweep double-typed must keep materializing the doubles and
    **delivers no transient saving at all** - it would be pure plumbing with no
    observable effect, which is neither gateable nor justifiable as a standalone
    slice. Either S4a carries the sweep, or S4a does not exist as a slice and
    the channel folds into S4b.
  - **S4b - typed sources.** `ColumnSource` names its pool; `ownedDenseValues`
    and `ownedTestValues` shrink to real-valued columns behind a per-column
    offset table; `rawColumn`/`rawTestColumn`/`rawColumnForRequantize` return
    null for factor columns by design; `writeOwnedDenseColumn`/`Cell` become
    kind-aware; `build`'s mapped arm and `buildTest` consume the code channel.
    Draw-preserving. **Large** - this is the data-ownership rework, and it is
    the slice that delivers the retained saving. Depends on S4a.
  - **S4c - the int-backed replay reader.** A third arm in
    `PredictorSourceColumnReader` and its construction in
    `PredictorSourceColumns`, so `predict` and `setTestPredictors` accept an
    integer source. Draw-preserving. Small-to-medium. Independent of S4b;
    depends on S4a only for the channel's existence.

**Ordering (decision 1): S0, S1, S2, S3, then S4a -> S4c -> S4b.** S4c precedes
S4b deliberately: it is small, independent of S4b, and exercises the code
channel end to end through `predict`, so the channel earns a second consumer and
a full gate run before the largest slice commits to it. S4b stays last - it is
the one slice that cannot be made small, it should carry its own review, and it
is the one to slip if anything must.

### Risks, and the honest case against doing this now

  - **The retained saving is narrower than the headline.** The common
    factor-data-frame path is the unmapped dense build, which retains no raw
    today (`data.hpp:1193-1198`). S4b's `-8n` applies to mixed and CSC
    containers and to the test side. A user fitting a plain data frame with
    factors and no test set gets the transient saving from S4a and nothing from
    S4b. If the argument for S4b were memory, that would be close to
    disqualifying; the argument is that the double-as-code representation is
    wrong, and the memory follows.
  - **S4b is a data-ownership rework, not a field addition.** It touches
    invariants `docs/design/data-ownership.md` records, both `build` arms,
    `buildTest`, the mutation entrances and the view gather. A partial landing
    mostly leaves `denseRaw` null where a reader expects a pointer - loud (a
    segfault, not a wrong number) - and the test matrix has to cover mixed and
    CSC builds with factor columns on both train and test sides.

    **Two exceptions where the failure is silent, and they are the two sites
    S4b rewrites.** `writeOwnedDenseColumn` and `writeOwnedDenseCell`
    (`data.hpp:1506-1521`) both open with a guard that returns quietly when the
    column is not `denseBorrowed` or its `denseRaw` is null. Under S4b a factor
    column legitimately has neither, so a mutation that should have written
    codes instead writes nothing and returns success, leaving the store's codes
    stale against the values the caller believes it installed. These guards must
    become assertions - or, better, the two functions must dispatch on kind so
    that "no double pool" routes to the code write rather than to the no-op
    return. This is the one place in S4b where "it will segfault if I get it
    wrong" is false, so it is the one place that needs a deliberate test rather
    than trust in the crash.
  - **The four silent-bitwise hazards above**, of which `hasMissing` is the
    serious one, because its failure mode is a shifted rng stream rather than a
    wrong value.
  - **Coupling to the width experiment.** "The codes are the raw" means a factor
    column's raw width becomes whatever `xint_t` is. If
    `docs/plans/archive/hot-layer-u8.md` later makes the code width per-column,
    the raw width becomes per-column too - a coupling that does not exist while
    a separate double pool holds the raw. This is a real narrowing of future
    freedom, and it is the strongest structural argument for the two-pool shape
    over the no-factor-pool one.
  - **Where it meets the backlogged dispatch reshape.**
    `docs/design/data-layout.md:319-354` names `partitionChildren`
    (`tree.hpp:866`) as "the one shared touchpoint" of the layout and dispatch
    axes and warns that landing them independently means reworking it twice.
    The good news is concrete: **this section does not touch
    `partitionChildren` at all**, because that function reads only codes and
    rank storage (`tree.hpp:888`, `:900`, `:912`) and never a double. The
    contact point is elsewhere - `PredictorSourceColumnReader`
    (`data.hpp:382-389`), where S4c adds a third arm and where a dispatch
    reshape would plausibly want to specialize the replay loop per column kind
    instead of branching per row. If that reshape lands first, S4c should
    express its int arm inside whatever tag mechanism the reshape introduces; if
    S4c lands first, the reshape inherits a three-arm reader rather than a
    two-arm one, which is more work but not a rewrite. Worth one sentence of
    coordination, not a blocking dependency.
  - **Opportunity cost.** S4b is the largest single item in this proposal and
    lands before an RC. If something has to slip, S4b is the slice to slip,
    because S1, S2, S3, S4a and S4c all stand without it.

### The narrow alternative, for the record

Landing **S4a only** - the code channel plus the bridge's code-typed validation
sweep, reaching the unmapped dense build - delivers the transient saving and the
avoided conversion, and leaves `codeFor(size_t, double)` and the mapped path's
`double* denseRaw` exactly as they are. It is a peak-memory optimization, not
the representational fix, and it permanently adds a second typed arm to both the
ingestion path and the validation sweep that S4b would then have to unify. It is
a reasonable place to *stop* if S4b slips; it is not a reasonable place to
*aim*. What it is emphatically not is *small*: once the sweep is counted, S4a is
medium-large, so "do the cheap half" is not what stopping here buys.

## 3. Invariant enforcement

### Where the check is today

`docs/design/data-store.md:330-331` assigns it to the host: "You own ingestion
validation and container assembly; the engine trusts its input beyond what you
check." `ColumnStore::build`'s own comment restates it (`data.hpp:1109-1111`:
"The host validates structure ... since the quantize trusts both"), and
`codeFor` acts on that trust with a bare `static_cast<xint_t>(value)`
(`data.hpp:708-711`). An out-of-range double-to-integer conversion is undefined
behavior in C++, and a code past K shifts past an inline mask
(`tree.hpp:153-155`) or over-reads a pooled one - the replay path clamps
defensively (`tree.hpp:1760-1769`) but `ColumnStore::codeAt` does not
(`data.hpp:489-492`).

As established above, every shipping host reaches `codeFor` through a validated
entrance. The unvalidated surface is the header-only engine used directly:
`tests/cpp`, and any future non-R host - which is precisely the audience the
R-agnosticism constraint is about.

### Alternatives

**A. Keep bridge-only trust.** Zero cost, zero change; leaves the engine
unsound standing alone, which matters more once the kind axis makes the store's
typing model something a non-R host would want to drive.

**B. Per-cell check inside `codeFor`.** `codeFor` has exactly six callers
(`data.hpp:959`, `:985`, `:997`, `:1005`, `:1769`; `sampler.hpp:1870`), all on
O(n) ingestion or mutation paths and none in the MCMC hot loop, so the marginal
cost is one compare-and-branch per cell on a cold path. The problem is not cost:
it turns a total function into a partial one, so it needs a sentinel return or a
throw, and `Rf_error` longjmps past destructors (`data.hpp:240-243` documents
that constraint) - so the failure has to travel back to a bridge that can raise.
That is a return-value protocol threaded through six call sites for a check
better placed one layer up.

**C. Check once per column at the store's own ingestion entrances**, fused into
a pass that already exists. `buildCutsForColumn` already calls
`inferredCategoryCount`, which touches every cell of a categorical column
(`data.hpp:842-847`); making it refuse a negative, non-integral or
unrepresentable value costs nothing extra. On the mutation side
`cutsWouldRemainValid` (`data.hpp:923-928`) already runs the correct per-cell
test via `categoricalValueIsValid` (`data.hpp:721-725`) and is already called by
the mutation transaction (`sampler.hpp:1821`); the gap is that `build`,
`buildTest` and the CSC quantize paths have no equivalent.

**D. Make the invariant unrepresentable on the code path.** With section 2's
S4a, ingestion performs exactly one checked narrow int32 -> `xint_t` per cell,
so the path that carries every factor after that lands is checked by
construction and the double path is the only one still needing C.

### Recommendation

**C, with D following from S4a.** The engine checks at its own ingestion
entrances and refuses by return value; the bridge keeps its checks unchanged.
This is defense in depth rather than a relocation, and the ordering matters: the
bridge's checks run *before* the store is touched, which is what preserves the
nothing-mutated-on-refusal property that `setState`'s rollback
(`sampler.hpp:945-955`) and the mutation transaction depend on. Message quality
for R users does not change, because the bridge still fires first with its own
wording (`R_interface_bartcore.cpp:1635-1638`).

S3 also carries the engine-side backstop for S2's level ceiling. The
K > 65534 refusal specified in section 1 sits beside `R/utility.R:500-507`,
which is R, so a header-only host reaches `build` with a wider ordered factor
and gets a grid silently clamped by `maxNumCutsRepresentable`
(`data.hpp:1155-1157`) - a merge, which is the defect S2 removes. `build` must
refuse the column itself. This is precisely the residue class S3 exists for:
the R message stays the one a user sees, and the engine stops being unsound
standing alone.

Cost: one fused per-cell test in `inferredCategoryCount` (free - the pass
already runs), one new refusal path out of `build`/`buildTest`, and a per-cell
range test on the code-channel narrow. No hot-loop cost. Draw-preserving for any
data that passes today.

## 4. FlatKind redundancy

### Correcting the premise

`FlatKind` is not a second copy of type-plus-pooled that exists by accident.
Its declaration says why it exists (`tree.hpp:169-170`): "kept in flags so the
replay family routes raw predictors *without the store that typed the columns*",
and the signature bears it out - `partitionFlatIndices` takes a `Columns` and a
`FlatNode` and no `ColumnStore` (`tree.hpp:1750-1753`). Storeless replay is a
shipped capability: `dbarts_sampler_getTrees` hands flat trees to consumers and
stan4bart calls it (`~/Repositories/stan4bart/src/init.cpp:502`).

What *is* redundant is the **cross-check**: `buildFromFlatBelow`
(`tree.hpp:1530-1574`) and `flatSubtreeIsWellFormed` (`tree.hpp:2065-2088`) each
re-derive the expected kind from `data.types[variable]` and
`data.columnIsPooled(variable)` and refuse a mismatch, in two separately
written three-arm ladders.

### Alternatives

**A. Keep both as they are.** Two parallel ladders to maintain; the tag must
stay in sync with the store at every write site.

**B. Drop `FlatKind`, pass the store to the replay.** Kills storeless replay,
which is a shipped feature with a live consumer. Refuse.

**C. Keep `FlatKind`, drop the cross-check, trust the tag.** `buildFromFlat`
accepts hand-built states - `sampler.hpp:1003-1013` shows the codebase already
reasoning about hand-buildable donors - so an untrusted tag against a store of
a different type would mis-route silently. Refuse.

**D. Keep both, unify the check, and keep `ColumnKind` out of the flat
format.** This is the interaction with section 1. An `orderedFactor` column
splits by threshold, so its flat node's tag is `FlatKind::ordinal`, unchanged;
`FlatKind` stays a four-valued *mechanic* tag while `ColumnKind` is a
three-valued *semantic* tag, and they are not two copies of one thing. The
cross-check then collapses to a single shared predicate -
`flatKindOf(node) == expectedFlatKind(data, variable)` - used by both ladders,
each of which keeps only its own payload validation.

### Recommendation

**D.** Memory zero: the tag already rides bits 1-2 of an existing flags byte
(`tree.hpp:198`, `:201-203`), three of four values are used, and bits 3-7 are
free. Sites: two three-arm ladders become one predicate plus their payload arms.

Headroom worth noting and not spending: if a future kind ever needed its own
flat tag, the flags byte has one spare `FlatKind` value and five free bits, and
the state's `tree.flags` block is a `RAWSXP` per node
(`R_interface_bartcore.cpp:5404-5434`), so widening the tag inside the byte
would be additive and would not move the state version -
`docs/design/public-surface.md:149-152` makes exactly that argument for the
MIA flags field.

## 5. The denseBorrowed / denseOwned naming inversion

### The fact

`ColumnSourceKind::denseBorrowed` (`data.hpp:133-134`) names the case where
`denseRaw` points into the store's **own** `ownedDenseValues`: the block is
assigned by the store at `data.hpp:1170-1181` and the descriptor is pointed at
it at `:1201-1203` (and at `:1339-1341` for the test side). The field comment
already concedes it - "denseBorrowed: the store's own dense column ... writable
so a mutation keeps the raw current" (`data.hpp:148-150`). `denseOwned`
(`:130-132`) names the opposite case, where the store owns nothing and
re-quantizes from whatever matrix the caller supplies
(`rawColumnForRequantize`, `data.hpp:655-660`). Twenty-one references across
`src/`.

### Alternatives

**A. Leave it.** Every reader of `data.hpp` must un-invert the name from the
comment, and `docs/design/data-ownership.md` - the one place a reader goes for
exactly this question - is contradicted by the code it describes.

**B. Swap the two spellings.** Refuse. Every existing reference would silently
mean the opposite of what it meant; a stale mental model is worse than an odd
name.

**C. Rename to the axis the code actually branches on.** The discriminator is
not ownership, it is *where the re-quantize source lives*:
`denseCallSupplied` (was `denseOwned` - the raw arrives with the call) and
`denseResident` (was `denseBorrowed` - the raw lives in the store), with
`denseRaw` becoming `residentRaw`. This is exactly the question
`rawColumnForRequantize` asks (`data.hpp:655-660`) and exactly what
`hasRequantizeSource`/`acceptsNewRawPredictors` (`:641-649`) are about. Neither
new spelling collides with an existing identifier.

**D. Rename ownership-honestly**: `denseUnowned` / `denseOwned`, so
`denseOwned` means store-owned. Correct, but it reuses an existing name for its
opposite - the same failure B has, in a smaller dose.

### Recommendation

**C**, in S1. Pure rename: four enumerator and field spellings, 21 references,
zero memory, zero behavior, no draw movement. Doing it in S1 rather than
alongside S4b matters: S4b rewrites what these names describe, and it is much
easier to review a rework of `denseResident` than a rework of something called
`denseBorrowed` that is not borrowed.

## 6. Bitwise and RNG impact

The equivalence gate expects identical draws while sampling is untouched;
baselines are re-recordable pre-release but at real cost. With the grid fix ruled
in, exactly one slice changes draws, and the bill for it is itemized below.

### Draw-preserving by construction

  - **S1**, the kind axis, *provided* the bridge's `varTypes` parse
    (`R_interface_bartcore.cpp:1060-1065`) is rewritten in the same commit as
    the R-side widening. If it is not, every ordered factor silently becomes
    categorical and every draw moves - the failure is loud in results and silent
    in review.
  - **S1's `categoryCounts` field**, *including* the `buildFromParent` copy.
    Every consumer reads the same value from a different array; omitting the
    copy is not a draw change but an out-of-range read.
  - **S1's `FlatKind` predicate collapse and the rename.**
  - **S3**, on data that passes today. A refusal does not move a draw.
  - **S4a/S4b/S4c**, subject to the four hazards enumerated in section 2's
    hot-path subsection - `hasMissing` flipping, the `NA_INTEGER` code, the
    signed max in `inferredCategoryCount`, and the CSC implicit-row rule. For an
    integral k in [0, 65535), `static_cast<xint_t>(double(k))` and
    `static_cast<xint_t>(std::int32_t(k))` are the same integer, so the codes
    themselves are identical; everything that could go wrong is at the
    boundaries.

### Draw-changing: S2, the midpoint grid

`numCuts[j]` changes for every ordered-factor column, so every code in the
column changes, `splitInterval`'s range changes (`tree.hpp:416-437`),
`logCut = -log(right - left + 1)` changes (`grow.hpp:270`), and
`scanOrdinalCuts` emits a different number of candidates into a
different-sized discrete draw. Every draw after the first rule on such a column
diverges. This is intended - it is the fix - and it is the only slice with this
property.

### The bill, and a gate that currently cannot see it

  - **`benchmarks/R/equivalence.R` will not move, and that is the problem.** Its
    factor cases are all built with plain `factor(...)` (`:185-186`, `:195-204`,
    `:871`, `:899`); the only `ordered` occurrences are `ordered_result = TRUE`
    at `:351` and `:523`, which is the *response*. There is no ordered-factor
    **predictor** anywhere in the corpus, so a grid change would leave the gate
    green while changing user-visible fits, and the gate cannot detect an
    ordered-factor regression at all. **S0 must add one, and must land and
    re-record before S2**, so S2 has a baseline to move. Section 2's hazards
    argue for two more cases in the same pass: an NA-bearing factor column
    (which is what makes a `hasMissing` regression visible) and a CSC-backed
    factor column.
  - **tinytest replay-regeneration - one file, not two.** The relevant question
    is not which files contain an ordered factor but which ones *assert* a
    quantity S2 moves, and R's rng stream is not one of those quantities.
    Sampling never advances R's stream: each chain gets its own Mersenne
    twister, seeded once at creation from `unif_rand` inside the only
    `GetRNGstate`/`PutRNGstate` bracket on the path
    (`R_interface_bartcore.cpp:1881-1889`), and the comment there states it
    outright - "sampling itself never advances R's stream" (`:1861-1862`), with
    the run entry restating it at `:4535`. The number of creation draws is one
    per chain regardless of any grid. So a file's seeded `rnorm()` sequence is
    **invariant** to S2.
      * `inst/tinytest/test-data-categorical.R` **does** need replaying. Not for
        stream reasons: the `o` column is an ordered factor in the `df` that
        feeds the samplers at `:86` and `:114`, so the fits asserted after
        `sampler$run(150L, 300L)` (`:93`) and `sampler.keep$run(50L, 5L)`
        (`:120`) change because the model changes. S1 separately changes the
        `varTypes` assertion at `:25` to `c(1L, 2L, 0L, 0L, 1L)`, and S2's cap
        rule changes the comment at `:40-41` and wants two new cases beside the
        54-level one at `:47-51`.
      * `inst/tinytest/test-sparse-factor.R` does **not**. Its ordered-factor
        block (`:1043-1075`) runs `sampler.ord$run(50L, 1L)` but asserts only
        error patterns from `predict`/`setTestPredictor`, never a fit. The
        file's own comment at `:1101-1103` - "a new `rnorm()` call here would
        shift every seeded draw that follows it" - is about adding an `rnorm()`
        call, which does advance the stream; `run()` does not, so the comment
        does not extend to S2. Re-verify the file after S2, but do not budget a
        regeneration for it.
      (`test-argument-surface.R:31` and `test-calibration-prior-draws.R:217` use
      ordered *responses*, not predictors, and are unaffected.)
  - **`tests/cpp`.** The C++ component tests build stores directly and assert
    against grids and codes; S1's `numCuts`-to-`categoryCounts` move and S2's
    grid both land there. Section 2's design also asks for one new assertion
    here: that an `orderedFactor` column's code equals its level index under the
    midpoint grid, with a companion case covering an unobserved interior level -
    which holds under a declared-levels grid and does not under an
    observed-levels one, so the assertion is also how that sub-decision gets
    pinned.
  - **`benchmarks/R/bench-sampler.R`.** S2 changes the candidate count in
    grow-from-root for ordered-factor columns, and S4c adds a branch to the
    replay reader. Neither baseline is invalidated by S1, but both want a
    measurement after S2 and after S4c.

### Recommended sequencing

1. **S0**: extend the equivalence corpus (ordered-factor predictor, NA-bearing
   factor, CSC-backed factor); re-record. No product change.
2. **S1**: the kind axis. Verify **identical** draws against the S0 baseline.
3. **S2**: the midpoint grid and the cap treatment. Re-record;
   replay-regenerate `test-data-categorical.R`; update `tests/cpp`.
4. **S3**, then **S4a**: engine-side validation (with the level-ceiling
   backstop), then the code channel. Identical draws against the S2 baseline.
5. **S4c**, then **S4b**: the int-backed reader, then typed sources. Identical draws
   again, with the four hazards as the specific things the corpus is there to
   catch.

## Interactions with the out-of-scope items

**Per-column code width** (`docs/plans/archive/hot-layer-u8.md`). That plan
records itself as rng-neutral ("codes are comparison operands; identical
comparisons, identical draws"), and nothing here disturbs that: the store's
cells stay `xint_t` and the code channel is a *host* width, narrowed once at
ingestion. `categoryCounts` helps the experiment - a width decision wants to
read K, not "the cut count" - and `scan.hpp:301-304` sizes its histogram from
that same value on the hot path, so the two changes meet there. The one genuine
tension is noted in section 2's risks: S4b's "the codes are the raw" makes a
factor column's raw width follow the code width, which a per-column width
experiment would then propagate further than it does today.

**Hot-path dispatch reshape of `partitionChildren`/`scan`.**
`docs/design/data-layout.md:319-354` names `partitionChildren` (`tree.hpp:866`)
as the shared touchpoint of the layout and dispatch axes. **Nothing in this
proposal touches it** - it reads codes and rank storage only (`tree.hpp:888`,
`:900`, `:912`). The contact point is the replay reader
(`PredictorSourceColumnReader`, `data.hpp:382-389`), where S4c adds a third arm;
see section 2's risks for how the two should be ordered. `splitsBySubset(j)`
becomes the natural dispatch key if that reshape lands a per-kind
specialization, and `ColumnKind` the natural tag parameter.

## Decisions

Every question this proposal raised is ruled below. The alternatives and their
costs are kept in the sections above rather than erased - the rulings conclude
that reasoning, they do not replace it.

1. **Section 2's scope: the FULL REWORK, in the order S4a -> S4c -> S4b.**
   DECIDED (2026-09-01, at orchestrator discretion under VD grant). S4a alone
   buys peak memory while leaving the double-as-code representation in place,
   which is the thing the item exists to remove, so stopping there was never the
   aim. The ordering puts S4c before S4b because S4c is small, is independent of
   S4b, and exercises the code channel end to end through `predict` - so the
   channel earns a second consumer and a gate run before the largest slice
   commits to it. S4b remains the one slice that cannot be made small and the
   one to slip if anything must.

2. **Ordered factors REMAIN ADMISSIBLE as leaf covariates and as monotone
   directions, deliberately.** DECIDED (2026-09-01, at orchestrator discretion
   under VD grant). Under the `splitsBySubset` predicate this is what
   `facade.hpp:770-773`, `:827-831` and `:867-868` already do; the ruling makes
   it intentional rather than inherited. An ordered factor designated as a
   linear- or GP-leaf covariate enters as standardized integer level codes,
   i.e. treated as equally-spaced interval data, which is the same reading that
   makes threshold splits right for the kind. Monotonicity was already settled
   in favour (`docs/design/monotone.md:117-121`). Consequence for section 2:
   an `orderedFactor` is the one factor kind that can need a double for a
   non-splitting reason, and it is served the way every other leaf covariate is
   - gathered into `gatheredRawValues` (`data.hpp:596`) by designation - not
   from a retained pool.

3. **Retained-raw shape: NO FACTOR POOL.** DECIDED (2026-09-01, at orchestrator
   discretion under VD grant). Factor columns retain nothing beyond their
   codes, which is what makes the `-8n` per column real; a second int32 pool
   would halve the saving while keeping a duplicate representation, which is
   the thing the section exists to remove. **Accepted cost:** a factor column's
   raw width now follows the code width, so if
   `docs/plans/archive/hot-layer-u8.md` later makes that width per-column, the
   per-column choice propagates further than it does today. That narrowing of
   the width door is accepted knowingly rather than discovered later.

4. **Grid levels: DECLARED.** DECIDED (2026-09-01, at orchestrator discretion
   under VD grant). The grid is a property of the factor rather than of the
   sample, which is what the level table is for, and it keeps an unobserved
   interior level splittable for a test row that carries it.
   `readDeclaredCategoryCounts` (`R_interface_bartcore.cpp:889-906`) is extended
   past **both** of its gates - the type gate at `:896` and the
   `sourceOf(j) < 0` gate at `:897`, so a CSC-backed ordered factor takes its
   container's declared K rather than an inferred count. `K` in the cap rule of
   section 1 means the **declared** K throughout.

5. **Test-side level check: STRICT REFUSAL, matching categorical.** DECIDED
   (2026-09-01, at orchestrator discretion under VD grant). An ordered-factor
   test value that is not an existing level code refuses, with the wording and
   at the sites the categorical check already uses
   (`R_interface_bartcore.cpp:1637-1638`, gates at `:1704`, `:1728`, `:3028`).
   The permissive reading - a between-levels value as interpolation on an
   ordered scale - is rejected: a column the user declared a factor is a factor
   on both sides of the fit, and a user who wants interpolation has the numeric
   column, where it is unambiguous. This closes the last of the four defects
   enumerated in section 1.

6. **`numCuts` de-overload: YES.** DECIDED (2026-09-01, at orchestrator
   discretion under VD grant). `numCuts[j] = 0` for a categorical column, a real
   `categoryCounts` field at 4 bytes per column carrying K, copied in
   `buildFromParent` (`data.hpp:1413-1427`), and the verbose printer's "Number
   of cutoffs" line (`R_interface_bartcore.cpp:1520-1521`) made kind-aware -
   reporting the category count for a factor column rather than 0.

7. **Provenance divergence: one UPGRADING bullet**, folded into the scheduled
   NEWS rewrite rather than written as a standalone note. DECIDED (2026-09-01,
   at orchestrator discretion under VD grant). The bullet says what the section
   1 paragraph establishes: a `dbartsData` saved before this change keeps
   `varTypes` 0 for its ordered factors, cannot be upgraded because 0 is
   ambiguous, and therefore fits differently from a freshly built object over
   the same data frame.

## Landing

### S0 - gate enabling (8cd10833, 95d1375e)

Adds three predictor-column shapes to benchmarks/R/equivalence.R -
ordfactor (an ordered-factor predictor at K = 150 levels against the
default n.cuts = 100), nafactor (two NA-bearing factor columns), and
sparsefactor (a CSC-backed sparseFactor column) - and records the new
baseline, equivalence-ee5ffe74.rds, 46 scenarios. Draw-preserving: all 43
predecessor scenarios reproduce equivalence-736bfb05.rds bitwise;
bcf-equivalence-00cfa108 stays 12/12 and multinomial-equivalence-4d9a3337
stays 11/11. No product code. S1 next.
