# The multiplier combiner: design

Status: LANDED, Gaussian-complete, 2026-08-13 to 2026-08-14, across M4.0
(562ee684), M4.1 (1458328c + e48fc5de), M4.2 (1a2aaedc) and M4.3 (9c63e9d8),
with the landing records at e7708b7c and 54e114ff
(docs/plans/multiforest-extension-surface.md, M4). The general K-forest
basis/amplitude family: each forest carries its own basis, its row contracts
with that forest's amplitude vector into a per-observation scalar, and the
forest enters the combination scaled by it. bcf's `a mu + b_z tau` is the
K = 2 instance. Posterior-defining: `BCFForestCombiner<L>`
(src/bartcore/combiner.hpp:685), the K-forest chain constructor
(chain.hpp:692-731), `expandForestSpecs` (combiner.hpp:325), the R resolution
(R/model.R `resolveForests`, R/spec.R), and the gates
(benchmarks/R/bcf-equivalence.R, bcf-exact{,-restricted,-weak}.R). Scope:
GAUSSIAN responses only (R/spec.R:423-431 refuses a basis off gaussian,
chain.hpp:702-705 builds `GaussianResponse` unconditionally; M4.4 is the
non-Gaussian slice), and CONSTANT leaves only (combiner.hpp:686-687).

This file documents ONE INSTANCE of the combiner hierarchy.
docs/design/forest-combiner.md owns the abstraction - why `ForestCombiner<L>`
exists beside `Forest<L>`, the null-short-circuit fast path, the
per-sweep-never-per-observation dispatch rule, the virtual surface as an
interface, and what stays Chain-level. docs/design/bcf.md owns the CAUSAL
instance - propensity, the moderator subset, tree-count and prior defaults, the
sd(y)-unit calibration, the exact-posterior gate, the burn-in study, and the
public creation history. Neither restates the math below.

**The naming debt, recorded and not acted on.** The family is general; the
SPELLING is bcf's. `BCFSpec`, `BCFForestSpec`, `BCFState`, `BCFForestCombiner`,
`createBCFSampler`, `ChainStateData::hasBCF` and the state format's `"bcf"`
block all read "BCF" where they mean "carries amplitudes", and the code's own
Doxygen has already re-described them as general (facade.hpp:774-776,
chain.hpp:687-689). A rename is ENGINE work, not a docs edit: it costs a
`--preclean`, a `structSize` move on the state fixtures (`common.cpp` already
tripped once at 272 -> 344, plan :3589-3591), and a state-block key change
needing the same in-place re-encode M4.3 used. Tracked as the root TODO's
`bcf-naming-generalization`. The mitigation that makes the debt survivable is
that no consumer-facing PROBE is keyed on the name or on a forest count:
capability is `totalAmplitudes() != 0` throughout (below, "Surfaces").

## The model

Forest f carries an n x q_f basis `B_f` and a length-q_f amplitude vector
`a_f`. The mean is

    E[y_i] = sum_f m_{f,i} f_f(x_i),   m_{f,i} = dot(a_f, B_f(i, .))

The per-row multiplier stays a SCALAR, which is why this is a generalization of
bcf's `a mu + b_z tau` and not a second mechanism: the reparameterization, the
combination and the reporting channels all read one number per forest per row.

`q_f >= 1` always, and there is NO implicit all-ones column. A forest whose
multiplier is a plain amplitude carries the ones column densely and reaches it
by the same contraction every other forest uses, which is what leaves exactly
one multiplier path (combiner.hpp:386-390).

Storage is ROW-major, row i at `i * numColumns`, because the contraction is the
only read the engine makes of a basis: a row is contiguous and the multiplier
costs one stream per forest (combiner.hpp:380-394). The plan's "The channel,
specified" section said COLUMN-major and was wrong about the shipped code;
M4.3 item 5 resolved the transpose IN THE DOC DIRECTION, the header agrees
(dbarts.h:707-712, src/C_interface.cpp:775-781), and M4.5 corrected the plan
line itself (docs/plans/multiforest-extension-surface.md:557-561).

## The amplitude layout

The amplitudes are ONE flat vector, forest f's block at `amplitudeOffset[f]`,
the offsets a pure prefix sum of the widths (combiner.hpp:412-418, :1388-1408).

The vector is RAGGED, and a TOTAL IS NOT A LAYOUT: `q = (1, 3)` and
`q = (2, 2)` both carry four amplitudes, so the per-forest widths travel with
them or a restore silently permutes the blocks. That one sentence is what the
whole persistence contract turns on (combiner.hpp:92-98, :666-672).

bcf's named accessors `a()`, `b0()`, `b1()` index THROUGH the offsets rather
than at 0/1/2 (combiner.hpp:454-459): forest 1's block starts wherever forest
0's ends, so a prognostic basis wider than one column moves them.

## The reparameterization

`formForestResponse` (combiner.hpp:815-835) hands forest f the pair its own
constant-leaf node sums need: response `r_i / m_{f,i}` and weight
`w_i m_{f,i}^2`, where `r_i` is the residual net of every OTHER forest's scaled
contribution.

A multiplier indistinguishable from zero at `0x1p-26` (combiner.hpp:787,
applied :823-827) leaves the row with exactly zero weight AND exactly zero
response. The zero response is required rather than cosmetic: the chain reads
this buffer arithmetically when it rolls the running residual and finalizes
total fits, and the node sufficient-statistic kernels accumulate `w * y`, where
a zero weight against an amplified response would be `0 * inf`. The tolerance
doubles as a condition-number cap - the division amplifies by at most 2^26, so
the cancellation downstream stays inside half the mantissa
(combiner.hpp:794-802).

The snap belongs to the REPARAMETERIZATION, not to the model: `combinedFits`
and the amplitude draw keep the exact multiplier, so a snapped row still
receives `m_{f,i} f_f(x_i)` in the combination and still informs the amplitude
conditional (combiner.hpp:804-806, :851-853).

The caller-settable per-forest, per-observation weight `s_{f,i}` composes as
one further multiplicative factor after this call returns, before the tree
loop. Its two edges - at `s = 0` with `m != 0` only the weight is zeroed, and
the weight lives on the chain rather than the serialized state - are
docs/design/bcf.md:159-169's.

## The amplitude conditional

Forest by forest in INDEX order, each block drawn jointly from its Gaussian
full conditional given the current value of every other block, so the pass is a
Gibbs scan and a block sees the blocks before it already updated
(combiner.hpp:1101-1127).

Forest f's design row is its basis row scaled by its own fit,
`x_i = B_f(i,.) f_f(x_i)`, against the residual net of every other forest, so
the conditional's precision is `P = I/priorVar + sum_i w_i x_i x_i' / sigma^2`
and its first moment `sum_i w_i x_i r_i / sigma^2`. The prior term is what
keeps P positive definite whatever the basis does, which is why the
factorization needs no failure path (combiner.hpp:1129-1136).

A scale-mixture prior's variance is refreshed straight after its own block,
from `IG((1 + q)/2, (scale^2 + ||a_f||^2)/2)` - at q = 1 exactly the shipped
path's shape 1.0 (combiner.hpp:1104-1107, :1116-1125).

**The LDL' fact, and why it is not the shipped Cholesky.** The factorization is
the square-root-free unit-lower `L D L'` because only its solve reduces to ONE
division per coordinate. Over an ORTHOGONAL basis - bcf's indicator pair, any
factor basis - the unit triangles are exactly identity, so the q-variate draw
is q scalar draws BITWISE, in coordinate order, one standard normal each
(combiner.hpp:1148-1152). The two-sqrt Cholesky solve gives
`x/sqrt(d)/sqrt(d) != x/d` and breaks the q = 1 reduction (M4.2 landing note,
plan :3488-3490). `testUnitLowerFactorization` (tests/cpp/test_model.cpp) is
its teeth, with a p = 1 arm asserting the Cholesky route DIFFERS and a p = 2
orthogonal arm (plan :3527-3532).

## Why bcf keeps a specialized draw

`drawGlue` selects `drawShippedGlue` on
`forests.size() == 2 && !generalAmplitudeDraw_ && shippedShape()`
(combiner.hpp:903-909). The two paths are the SAME conditional in exact
arithmetic - all four accumulators bitwise equal under `-ffp-contract=off` -
and differ only in where the compiler forms fused multiply-adds: the a block
accumulates in one statement and fuses, while the b block forms per-row
products before a branch and accumulates inside it, which fuses unevenly.

The MEASURED split, carried here verbatim from the code comment because it is
the trigger for deleting the branch: all four PRECISIONS reproduce bitwise,
weighted and unweighted; the divergence is in the MOMENTS - unweighted, n1
reproduces and n0 differs; weighted, both differ. No single accumulation shape
reproduces both blocks, over 21 variants tried (combiner.hpp:876-894). So the
general path CANNOT be bitwise on bcf, and the specialized one is kept until a
`bcf-equivalence` re-record is authorized. `BCFSpec::generalAmplitudeDraw` is
the one-line switch that re-record flips (combiner.hpp:307-312).

## The ASIS ridge

Per forest, at most one GIG draw each, in index order. The blocks are DISJOINT,
so the moves commute and each is an exact Gibbs update given the rest, which
makes the order a stream convention rather than a modelling choice
(combiner.hpp:911-923). "At most" is exact: `afterCombine` skips a forest
entirely on `!prior.update || !prior.ridge`, before any draw
(combiner.hpp:934), and `rescaleAmplitudeRidge` returns 1.0 consuming NO rng
below two occupied leaves or at a zero leaf sum (combiner.hpp:1242). Two
further 1.0 returns guard a non-finite or non-positive draw (:1258, :1260), but
those are reached only AFTER the GIG draw is taken, so :1242 is the sole
rng-free skip of the three. M4.0's pins hold all three inert on the stream.

Forest f's `L + q` scale coordinates travel the likelihood-invariant orbit
`(a_f, leaves) -> (a_f/c, c leaves)` with `c = sqrt(v)` and
`v ~ GIG((L - q)/2, M/leafVar, ||a_f||^2/priorVar)`, L and M the count and
squared sum of that forest's OCCUPIED leaves (combiner.hpp:1200-1217,
:1242-1268). The exponent is docs/plans/bcf-b-ridge.md:192-194's general rule
`p = (k - d)/2` for rescaling k leaf parameters against d glue scalars; the
naive move-map Jacobian's `(L - q + 1)/2` is off by one and that memo's
prototype rejects it at KS 1.6e-21 (:166-171, :329-357).

**One mechanism, not two.** Instantiated at q = 1 it IS bcf's shipped a-move
bitwise; at q = 2 with a fixed prior variance it IS the b-move
docs/plans/bcf-b-ridge.md derives (combiner.hpp:925-928).

B reads the LIVE prior variance, which for a scale mixture is the auxiliary
this move conditions on. Refreshing it here would re-randomize the coordinate
just conditioned on and throttle the mixing gain - measured, IACT 69 -> 196 on
`|a|` (docs/plans/bcf-ridge-interweaving.md:488-492); the one-sweep lag is
benign, the next `drawGlue` refreshing it given the new amplitude
(combiner.hpp:1207-1211).

The rescale-consistency set the move must carry, or a stored
`amplitude * leaf` stops being the identified product: `muByTree`, `totalFits`,
`totalTestFits`/`currTestFits` under `record`, and the keepTrees flattened slot
(combiner.hpp:1265-1285).

## ridgeB is code that is OFF

The b-move ships, but `BCFSpec::ridgeB = false` (combiner.hpp:306). Enabling it
costs a GIG draw per sweep, which re-records `bcf-equivalence`, and its own
acceptance gate (docs/plans/bcf-b-ridge.md:495-506 - IACT payoff, bcf-exact
mode-2b, keepTrees round trip) was not run. It is a DOOR, not a fork: it flips
only on a named measured mixing case, plus that gate, plus a re-record with the
`equivalence.yaml` bump in the same commit.

**No silent enablement is possible** (resolved at M4.3 review, plan
:3603-3610). On every creation route the scale mixture holds if and only if a
forest is basis-FREE. R writes the forest's amplitude prior SCALE as 0 whenever
that forest declares a basis (`R/model.R:874`,
`if (withBasis) 0 else declared(spec$sd, 2)`); the bridge reads it into
`ForestSpec::amplitudePriorScale` (src/R_interface_bartcore.cpp:2205) and
derives `forest.ridge = forest.amplitudePriorScale > 0.0` (:2207). So every
basis-carrying forest gets `ridge = false`, and the combiner's own
`halfCauchyScale` - the field `amplitudePriorScale` becomes - is zero there. The
one reachable q > 1 scale-mixture state, a post-creation widening of a
basis-free forest, consumes the SAME single GIG draw already taken at q = 1, at
the M4.2-validated exponent, so no stream moves.

## Draw-path selection is a VALUE predicate

Per forest, IS-CANONICAL: `basis[f]` is canonical when it is a dense all-ones
column, or a complementary two-column 0/1 pair (each entry in {0, 1}, each row
summing to 1) (combiner.hpp:1341-1363).

It is NOT a width test, and the reason is a silent-wrong-answer one:
`drawShippedGlue` never reads `basis[1]` as a DESIGN MATRIX. It borrows that
forest's column 1 and tests it for nonzero as a GROUP KEY, never multiplying by
the stored values, and forms two disjoint group-precision accumulators from the
partition (combiner.hpp:1045-1099). So a legal continuous 0.25/0.75 pair would
be drawn as a different MODEL - both entries are nonzero, so every row lands in
the treated group, and the values the general conditional would contract with
never enter. A non-canonical basis at ANY forest forces the general path for
the whole draw.

The flag is recomputed at install and at restore and never serialized, so
nothing can carry a stale answer across a round trip (combiner.hpp:421-429,
:757, :1016). Recorded as M4.3's doc-level residue (plan :3626-3629): the
RESTORE-side recompute is VACUOUS, because `restoreGlue` never touches basis
values - the load-bearing recompute is the install-time one. Verified at
combiner.hpp:1000-1017, which calls `refreshCanonical()` after writing
amplitudes only.

## One mutation route

Synthesis is CONSTRUCTION-ONLY: `synthesizeIndicatorBasis`
(combiner.hpp:1329-1338) is called once, from the constructor (:726).
`installForestBasis` is the SOLE mutator and therefore owns the guards nothing
else is left to apply - the index, `numColumns >= 1`, a non-null pointer, and
finiteness (combiner.hpp:747-759). It wins by being the only operation there
is, which is why no ordering between two mutators has to be specified
(combiner.hpp:694-700).

Ordering is LAST INSTALL WINS, per forest, and both orderings of a widen and a
swap collapse to it because `rebuildAmplitudeLayout` derives the offsets as a
pure prefix sum of the width vector and carries every block by position
(combiner.hpp:731-741, :1388-1408). Amplitudes PRESERVE and remap; surviving
coordinates keep their values at their new offsets and new ones enter at the
neutral 1.0. A width-preserving install is the BITWISE IDENTITY on every
amplitude (the :1396 early return) - that is bcf's mid-life z swap, and it is
baseline-gating.

`glue_.z` is DELETED. `drawShippedGlue` partitions from `basis[1]`'s treated
column, not from a borrowed indicator, because a width-preserving swap to a
different complementary pair would otherwise install the new pair, leave the
layout unmoved, and then draw `b0`/`b1` under the OLD partition while
`forestMultiplier` contracted the NEW basis (M4.3 item 1, plan :1903-1923;
landed per :3559-3560).

## Persistence

The AMPLITUDES are state; the BASES are not.

`serializeGlue`/`restoreGlue` carry the per-forest widths, the flat amplitude
vector, and each forest's prior variance (combiner.hpp:982-1017).
`glueIsValid` is the LAYOUT check `stateIsValid` routes through, because
`restoreGlue` writes THROUGH the live offsets and a state with the same total
over different widths would otherwise be admitted and silently permute the
blocks (combiner.hpp:666-672, :1018-1028).

The four named scalars `a`/`aVariance`/`b0`/`b1` survive in `ChainStateData` as
a hand-written K = 2 READING only, non-authoritative, read exactly when
`amplitudeWidths` is empty - which is exactly a state written by hand rather
than by a combiner (combiner.hpp:103-109, :1000-1011).

The bases ride CREATION, on `data@bases` (a LIST, R/A_class.R:514-524,
validated by `validateForestBases`, R/data.R:635-655), the way the design
matrix does. That is how RESTORE-THEN-WIDEN is met with no fourth reapply hook:
a widening applied after a restore preserves and remaps the RESTORED amplitudes
rather than the constructed ones (combiner.hpp:978-981; plan :1953-1973, landed
:3564-3568).

## Bitwise contracts

THREE accumulation directions, all observable once q > 1 or K > 2, each a
contract:

1. `combinedFits` accumulates WITHIN the row and from the LAST forest DOWN
   (combiner.hpp:854-869, seed at :863, descending loop :864-865). This is the
   load-bearing one and it is MEASURED. Under fused multiply-add contraction
   exactly one product in a sum escapes its own rounding - the one the closing
   add absorbs - and the two-forest `a mu + b_z tau` this replaces absorbed
   forest 0's. Seeding with the LAST forest leaves forest 0's product as the
   only bare multiply in the closing add. Accumulating FORWARD absorbs the last
   forest's product instead, moves ~30% of rows by one ulp and contaminates the
   trajectory within ~40 sweeps: all 12 `bcf-equivalence` scenarios red on mu,
   tau, glue, sigma and train (combiner.hpp:840-849; M4.1 landing note, plan
   :3434-3444). `testCombinedFitsAssociation`
   (tests/cpp/test_sampler.cpp:3210) is the ONLY in-process guard - the M4.0
   seam pin structurally CANNOT see association, its reference expression
   inheriting the test compiler's own contraction.
2. `formForestResponse`'s residual accumulates FORWARD, subtracting the other
   forests in increasing index order, and deliberately NOT combinedFits'
   reverse: this sum has no two-term fused expression to reproduce, and the
   amplitude conditional forms the same residual the same way, so the two agree
   by construction rather than by coincidence (combiner.hpp:808-814, :828-830).
3. `forestMultiplier` contracts FORWARD over the columns; at q > 2 a
   reassociation moves the multiplier by an ulp and every reader of it with it
   (combiner.hpp:1298-1303).

The standing lesson these pins carry, twice learned: a pin fixture must give
every factor in the pinned expression a DISCRIMINATING value - unit values
silently vacate pins (M4.0's required fix, plan :3381-3389; recurred at M4.3's
Arm 5(iii) dead pin, :3593-3598).

## The calibration map, general in K

Stated on `ForestSpec` rather than in the chain (combiner.hpp:258-287, applied
chain.hpp:709-727). The node scale is
`nodeScaleFactor * s / nodeScaleDivisor`, s the sample sd of the range-scaled
response.

Two branches, which are the two ways a magnitude can be carried. A
SCALE-MIXTURE amplitude carries it in `amplitudePriorScale` and leaves the node
scale at s - bcf's prognostic forest, factor and divisor both 1, giving
`mu ~ N(0, s^2)` with the half-Cauchy a at `aPriorScale sd(y)`. A
FIXED-VARIANCE amplitude carries it in the node scale, divided through the
half-normal median 0.674 - bcf's treatment forest, giving
`tau ~ N(0, (sdModerate s / 0.674)^2)`, so `(b1 - b0) tau` sits at
`sdModerate sd(y)`.

The divisor is CARRIED rather than derived from the prior, so a forest
declaring a zero scale keeps the node scale its host asked for. Both
expressions are written exactly as bcf's were, which is what keeps the K = 2
instance bitwise.

## bcf as the K = 2 instance

`expandForestSpecs` (combiner.hpp:319-340) is the thin adapter between bcf's
two-forest spelling and the K-length vector every other layer works in, and it
is LOAD-BEARING rather than courtesy: 25 `tests/cpp` fixtures and
benchmarks/R/bcf-equivalence.R drive through the `BCFSpec` spelling (plan
:1780-1786, :2177-2179). Forest 0 takes the half-Cauchy amplitude over the
implicit intercept and leaves its node scale at s; forest 1 takes the
fixed-variance pair over the treatment indicator basis and carries `sdModerate`
in its node scale.

`bcfGlue(a, b0, b1)` survives as a named READING that returns false on any
other layout, which is how a caller learns it is not looking at bcf
(combiner.hpp:761-770, chain.hpp:1040-1048).

## Surfaces

R creation: `forests = list(forest(basis = ...))` plus the per-forest knob map
(`resolveForests`, R/model.R), and `dbartsData(bases = )` for a numeric basis
(R/spec.R:406-431). R5: `$setForestBasis(forest, basis)` (R/dbarts.R:1223) and
`$getForestAmplitudes(forest)` (R/dbarts.R:1450), both 1-based via
`resolveForestIndex` (R/bartcore.R:1083). Flat C:
`dbarts_sampler_setForestBasis`, `dbarts_sampler_numForestAmplitudes` and
`dbarts_sampler_forestAmplitudes`, ragged and ROW-major
(inst/include/dbarts/dbarts.h:707-751).

Capability probes are `totalAmplitudes() != 0`, NEVER a forest count: a
K-forest multinomial defeats a `numForests` test, and the shipped code states
that rule in two independent places (src/C_interface.cpp:764-770,
src/R_interface_bartcore.cpp:3736-3741).

## What this family does not do

- Non-Gaussian responses. M4.4.
- The test surface. `testFitsAreDefined` and `logLikelihoodIsDefined` are both
  false (combiner.hpp:941-946), so `setTestPredictors`, `setTestOffset` and
  `predict` refuse - through `refuseBCFTestSurface`, gated on
  `testFitsAreDefined` rather than on the forest count
  (src/R_interface_bartcore_common.hpp:185-189).
- A per-draw amplitude channel in flat C. `dbarts_results` carries none
  (dbarts.h:139-152); DECLINED at plan :1521-1531, a `DBARTS_C_API_MINOR` bump
  binding decision 8 forbids.
- A variance forest. `createBCFSampler` refuses `numVarianceTrees > 0`
  (facade.hpp:782-785).
- Nameable leaf-prior calibration. The map owns it, so the write is refused on
  ANY combining sampler: `Chain::setForestPriorScale` returns `false` on
  `f >= forests_.size() || combiner_ != nullptr` (chain.hpp:1024), which the
  flat entry surfaces to a caller as a 0 return.

The classes the family ENABLES, and their evidence status, are
docs/design/model-space-survey.md's D4: continuous/dose-response exposure BCF
(unpublished preprint), VCBART's varying coefficients (which carries no
sampled amplitude), heterogeneous mediation, and the multiplier half of
principal stratification with BCF.

## Landing notes

**M4.0, 562ee684 (pins).** Component pins over BOTH `afterCombine` overrides -
BCF's GIG rescale with its three reachable 1.0 skips inert on the rng stream,
and the multinomial additive shift with its returns-1.0-having-moved
convention, whose pin is the SOLE guard on that convention. The review's
required fix became the slice's lesson: the fixture originally set
`numTrees = 1` on every forest, collapsing each per-tree factor to 1 and
leaving four sites undiscriminated, with the reviewer MEASURING a wrong
implementation staying green. Fixed with unequal per-forest tree counts. ~348
dense-equivalent of the TESTS band - the slice added no engine code at all, so
this figure is not comparable with the engine nets quoted for M4.1-M4.3 below.

**M4.1, 1458328c + follow-up e48fc5de (the multiplier generalization).**
`forestMultiplier` became `dot(a_f, B_f(i,.))` and `combinedFits` a K loop with
no K = 2 special case. The defining event was a 12/12 `bcf-equivalence`
divergence root-caused to FMA CONTRACTION ASSOCIATION, fixed by accumulating
from the last forest down (see "Bitwise contracts"). Engine +48 net
dense-equivalent, tests +98. The reviewer's trap note, standing: a perturbed
install WITHOUT `--preclean` reported bitwise identical, so every re-check must
`--preclean`. BENCH, on a granted quiet window: B/A per-sweep 1.0098 at
n = 20000 and 1.0105 at n = 2000, flat across 100x in n - per-element compute
in the combiner, not bandwidth - and the ~1% BCF-only cost was ACCEPTED as the
price of the general multiplier (plan :3470-3482).

**M4.2, 1a2aaedc (the q-variate conditional).** The per-forest q-variate
conditional through the new unit-lower LDL' helper, and ONE general per-forest
ASIS rescale at `p = (L - q)/2` whose q = 1 instance is bitwise-identical to
bcf's a-move. The in-slice rule's outcome: bitwise on the ASIS half,
SPECIALIZED-PATH-KEPT on the conditional half, with impossibility CONFIRMED by
the reviewer over 21 accumulation variants against the real engine (a first
standalone probe false-alarmed - vectorization split the multiplies in the
replica, a recorded methodology lesson). Engine ~180 dense-equivalent.

**M4.3, 9c63e9d8 (SUBSUME).** `setTreatment` retired as a mutator at all four
layers in favour of `setForestBasis(f, values, q)`; `installTreatment` became
the constructor-only `synthesizeIndicatorBasis`; `installForestBasis` became
the sole mutator; `glue_.z` DELETED; draw-path selection moved onto the
per-forest canonical VALUE predicate; `bcfGlue` re-signed to
`totalAmplitudes`/`numForestAmplitudes`/`amplitudes` with a K = 2 thin adapter;
`BCFSpec` gained its K-length `forests` vector with `expandForestSpecs`, all 25
BCFSpec fixtures untouched; `data@treatment` retired for a `data@bases` LIST
riding creation. Refusal ledger as the review re-derived it: 5 relax, 3
restate, 4 rewritten by the slot retirement, 33 stand, and ONE dropped outright
with licensed successors (the treatment-coding refusal, replaced by
length/finiteness refusals at creation). `consumer.c` `LEG_COUNT` 18 -> 19.
Engine ~165 dense-equivalent.

**Budget units, as a standing convention** (plan :3620-3624). Slice bands on
this arc are DENSE-EQUIVALENT lines - lines counted without blank-line and
formatting inflation, which roughly doubles raw counts. M4.3's implementer
reported raw nets against a dense band and appeared over budget when it was
85-135 UNDER. Never mix the two currencies in one report.

## Status

LANDED and Gaussian-complete. M4.4 (non-Gaussian responses) is the immediate
follow-on and owns dbarts.h:513's "Gaussian responses only." The ridgeB door is
shut. The naming debt is recorded above and in the root TODO.
