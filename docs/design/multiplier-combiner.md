# The multiplier combiner: design

Status: LANDED, 2026-08-13 to 2026-08-14, across M4.0 (eb340f4f), M4.1
(e3170da4 + 9459bef0), M4.2 (a35ff7df), M4.3 (983d7f0a) and M4.4 (e5e93f11),
with the landing records at 53d64798 and c0c9e12f
(docs/plans/archive/multiforest-extension-surface.md, M4). The general K-forest
basis/amplitude family: each forest carries its own basis, its row contracts
with that forest's amplitude vector into a per-observation scalar, and the
forest enters the combination scaled by it. bcf's `a mu + b_z tau` is the
K = 2 instance. Posterior-defining: `AmplitudeForestCombiner<L>`
([[src/bartcore/combiner.hpp#AmplitudeForestCombiner]]), the K-forest chain
constructor ([[src/bartcore/chain.hpp#Chain]]), `expandForestSpecs`
([[src/bartcore/combiner.hpp#expandForestSpecs]]), the R resolution
(R/model.R `resolveForests`, R/spec.R), and the gates
(benchmarks/R/bcf-equivalence.R, bcf-exact{,-restricted,-weak}.R). Scope:
GAUSSIAN, PROBIT and LOGISTIC responses ([[R/spec.R#resolveSamplerSpec]] admits
those three and refuses each other family by name,
[[src/bartcore/chain.hpp#Chain]] builds the matching response model,
[[src/bartcore/facade.hpp#createAmplitudeSampler]] the engine-side door), and
CONSTANT leaves only ([[src/bartcore/combiner.hpp#AmplitudeForestCombiner]]).
Under a latent family the combination is the
INDEX, on the link's own fixed scale: sigma is pinned, the response transform
is the identity, and every forest's prior scale is stated in latent sd units.

This file documents ONE INSTANCE of the combiner hierarchy.
docs/design/forest-combiner.md owns the abstraction - why `ForestCombiner<L>`
exists beside `Forest<L>`, the null-short-circuit fast path, the
per-sweep-never-per-observation dispatch rule, the virtual surface as an
interface, and what stays Chain-level. docs/design/bcf.md owns the CAUSAL
instance - propensity, the moderator subset, tree-count and prior defaults, the
calibration map's gaussian defaults, the exact-posterior gate, the burn-in
study, and the public creation history. Neither restates the math below.

**The naming debt, DISCHARGED.** The family is general and so, now, is the
spelling: `AmplitudeSpec`, `ForestStructureSpec` (the tree/structure block,
which carries no amplitude), `AmplitudeState`, `AmplitudeForestCombiner`,
`createAmplitudeSampler`, `ChainStateData::hasAmplitudes`, and the state
format's `"glue"` block - each of which used to read "BCF" where it meant
"carries amplitudes". R follows with `samplerCarriesAmplitudes` and
`refuseAmplitudeMutation`, and the control attribute is `bartcore.forests`.
Two costs this paragraph priced were wrong. `structSize` did NOT move: the
assert is `sizeof(ChainStateData) == 416` (`tests/cpp/common.cpp`), the M4.3
trip from 272 -> 344 was field ADDITION, and `hasAmplitudes` is the same `bool`
`hasBCF` was, so the only fixture edit was the field name in `statesAgree`. The
state-block change was a KEY rename, not M4.3's in-place re-encode, and a
rename is non-additive in the opposite direction: a reader looking up an absent
name defaults an OPTIONAL block in silence, leaving the amplitudes at their
construction values. So `stateFormatVersion` and
`minReadableStateFormatVersion` both went 1 -> 2, and a state carrying the old
key is refused by version before any block is read. No baseline was
re-recorded and all three equivalence suites stayed bitwise. The mitigation
that made the debt survivable in the meantime, and still holds, is that no
consumer-facing PROBE is keyed on the name or on a forest count: capability is
`totalAmplitudes() != 0` throughout (below, "Surfaces").

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
one multiplier path ([[src/bartcore/combiner.hpp#forestMultiplier]]).

Storage is ROW-major, row i at `i * numColumns`, because the contraction is the
only read the engine makes of a basis: a row is contiguous and the multiplier
costs one stream per forest ([[src/bartcore/combiner.hpp#forestMultiplier]]).
The plan's "The channel, specified" section said COLUMN-major and was wrong
about the shipped code; M4.3 item 5 resolved the transpose IN THE DOC
DIRECTION, the header agrees
([[inst/include/dbarts/dbarts.h#dbarts_sampler_setTreeStorage]],
[[src/C_interface.cpp#dbarts_sampler_printTrees]]), and M4.5 corrected the
plan line itself
([[docs/plans/multiforest-extension-surface.md:557-561@4c018187]]).

## The amplitude layout

The amplitudes are ONE flat vector, forest f's block at `amplitudeOffset[f]`,
the offsets a pure prefix sum of the widths
([[src/bartcore/combiner.hpp#AmplitudeState]],
[[src/bartcore/combiner.hpp#rebuildAmplitudeLayout]]).

The vector is RAGGED, and a TOTAL IS NOT A LAYOUT: `q = (1, 3)` and
`q = (2, 2)` both carry four amplitudes, so the per-forest widths travel with
them or a restore silently permutes the blocks. That one sentence is what the
whole persistence contract turns on ([[src/bartcore/combiner.hpp#ChainStateData]],
[[src/bartcore/combiner.hpp#glueIsValid]]).

bcf's named accessors `a()`, `b0()`, `b1()` index THROUGH the offsets rather
than at 0/1/2 ([[src/bartcore/combiner.hpp#a, b0, b1]]): forest 1's block
starts wherever forest 0's ends, so a prognostic basis wider than one column
moves them.

## The reparameterization

`formForestResponse` ([[src/bartcore/combiner.hpp#AmplitudeForestCombiner::formForestResponse]]) hands
forest f the pair its own constant-leaf node sums need: response
`r_i / m_{f,i}` and weight `w_i m_{f,i}^2`, where `r_i` is the residual net of
every OTHER forest's scaled contribution.

A multiplier indistinguishable from zero at `0x1p-26`
([[src/bartcore/combiner.hpp#zeroMultiplierTolerance]], applied inside
[[src/bartcore/combiner.hpp#AmplitudeForestCombiner::formForestResponse]]) leaves the row with exactly
zero weight AND exactly zero response. The zero response is required rather
than cosmetic: the chain reads this buffer arithmetically when it rolls the
running residual and finalizes total fits, and the node sufficient-statistic
kernels accumulate `w * y`, where a zero weight against an amplified response
would be `0 * inf`. The tolerance doubles as a condition-number cap - the
division amplifies by at most 2^26, whatever the family, the weights and K,
and the arithmetic downstream cancels that amplification EXACTLY rather than
merely bounding it: the node kernels accumulate `sumWeights = sum_i w_i
m_i^2` and `sumWeightedResponse = sum_i (w_i m_i^2)(r_i / m_i)`, whose exact
value is `sum_i w_i m_i r_i` ([[src/bartcore/combiner.hpp#AmplitudeForestCombiner::formForestResponse]]).

The snap belongs to the REPARAMETERIZATION, not to the model: `combinedFits`
and the amplitude draw keep the exact multiplier, so a snapped row still
receives `m_{f,i} f_f(x_i)` in the combination and still informs the amplitude
conditional ([[src/bartcore/combiner.hpp#combinedFits, AmplitudeForestCombiner::formForestResponse]]).

The caller-settable per-forest, per-observation weight `s_{f,i}` composes as
one further multiplicative factor after this call returns, before the tree
loop. Its two edges - at `s = 0` with `m != 0` only the weight is zeroed, and
the weight lives on the chain rather than the serialized state - are
[[docs/design/bcf.md#The multiplier snap and the per-forest weight]]'s.

## The amplitude conditional

Forest by forest in INDEX order, each block drawn jointly from its Gaussian
full conditional given the current value of every other block, so the pass is a
Gibbs scan and a block sees the blocks before it already updated
([[src/bartcore/combiner.hpp#drawAmplitudes]]).

Forest f's design row is its basis row scaled by its own fit,
`x_i = B_f(i,.) f_f(x_i)`, against the residual net of every other forest, so
the conditional's precision is `P = I/priorVar + sum_i w_i x_i x_i' / sigma^2`
and its first moment `sum_i w_i x_i r_i / sigma^2`. The prior term is what
keeps P positive definite whatever the basis does, which is why the
factorization needs no failure path
([[src/bartcore/combiner.hpp#drawForestAmplitude]]).

A scale-mixture prior's variance is refreshed straight after its own block,
from `IG((1 + q)/2, (scale^2 + ||a_f||^2)/2)` - at q = 1 exactly the 1.0
bcf's own spelling carries ([[src/bartcore/combiner.hpp#drawAmplitudes]]).

**The LDL' fact, and why it is not the shipped Cholesky.** The factorization is
the square-root-free unit-lower `L D L'` because only its solve reduces to ONE
division per coordinate. Over an ORTHOGONAL basis - bcf's indicator pair, any
factor basis - the unit triangles are exactly identity, so the q-variate draw
is q scalar draws BITWISE, in coordinate order, one standard normal each
([[src/bartcore/combiner.hpp#drawForestAmplitude]]). The two-sqrt Cholesky
solve gives `x/sqrt(d)/sqrt(d) != x/d` and breaks the q = 1 reduction (M4.2
landing note, [[docs/plans/multiforest-extension-surface.md:4855-4857@4c018187]]).
`testUnitLowerFactorization` (tests/cpp/test_model.cpp) is its teeth, with a
p = 1 arm asserting the Cholesky route DIFFERS and a p = 2 orthogonal arm
([[docs/plans/multiforest-extension-surface.md:4891-4896@4c018187]]).

## Why bcf no longer keeps a specialized draw

`drawGlue` ([[src/bartcore/combiner.hpp#drawGlue]]) is the general sweep and nothing else. bcf's
K = 2 shape landed with a two-scalar specialization beside it, selected on the
forest count, a basis-shape predicate and an `AmplitudeSpec` flag. The two were
the SAME conditional in exact arithmetic - all four accumulators bitwise equal
under `-ffp-contract=off` - and differed only in where the compiler forms fused
multiply-adds: the a block accumulated in one statement and fused, while the b
block formed per-row products before a branch and accumulated inside it, which
fused unevenly.

The MEASURED split, kept because it is what the deletion cost: all four
PRECISIONS reproduced bitwise, weighted and unweighted; the divergence was in
the MOMENTS - unweighted, n1 reproduced and n0 differed; weighted, both
differed. No single accumulation shape reproduced both blocks, over 21 variants
tried. So the general path could not be bitwise on bcf, and the specialization
was held until a `bcf-equivalence` re-record was authorized. This IS that
re-record: the branch, its predicate and its spec flag are gone, and every
shape draws through the one conditional.

## The ASIS ridge

Per forest, at most one GIG draw each, in index order. The blocks are DISJOINT,
so the moves commute and each is an exact Gibbs update given the rest, which
makes the order a stream convention rather than a modelling choice
([[src/bartcore/combiner.hpp#afterCombine]]). "At most" is exact: `afterCombine`
skips a forest entirely on `!prior.update || !prior.ridge`, before any draw
([[src/bartcore/combiner.hpp#afterCombine]]), and `rescaleAmplitudeRidge`
returns 1.0 consuming NO rng below two occupied leaves or at a zero leaf sum
([[src/bartcore/combiner.hpp#rescaleAmplitudeRidge]]). Two further 1.0 returns
guard a non-finite or non-positive draw, but those are reached only AFTER the
GIG draw is taken, so the leaf-count/leaf-sum guard is the sole rng-free skip
of the three ([[src/bartcore/combiner.hpp#rescaleAmplitudeRidge]]). M4.0's pins
hold all three inert on the stream.

Forest f's `L + q` scale coordinates travel the likelihood-invariant orbit
`(a_f, leaves) -> (a_f/c, c leaves)` with `c = sqrt(v)` and
`v ~ GIG((L - q)/2, M/leafVar, ||a_f||^2/priorVar)`, L and M the count and
squared sum of that forest's OCCUPIED leaves
([[src/bartcore/combiner.hpp#rescaleAmplitudeRidge]]). The exponent follows the
general rule `p = (k - d)/2` for rescaling k leaf parameters against d glue
scalars; the naive move-map Jacobian's `(L - q + 1)/2` is off by one and the
b-move's prototype (q = 2) rejects it at KS 1.6e-21 - derived and evidenced
below, "The exponent rule".

**One mechanism, not two.** Instantiated at q = 1 it IS bcf's shipped a-move
bitwise; at q = 2 with a fixed prior variance it IS the b-move
docs/plans/archive/bcf-b-ridge.md derives
([[src/bartcore/combiner.hpp#rescaleAmplitudeRidge]]).

B reads the LIVE prior variance, which for a scale mixture is the auxiliary
this move conditions on. Refreshing it here would re-randomize the coordinate
just conditioned on and throttle the mixing gain - measured, IACT 69 -> 196 on
`|a|` ([[docs/plans/bcf-ridge-interweaving.md:488-492@4c018187]]); the
one-sweep lag is benign, the next `drawGlue` refreshing it given the new
amplitude ([[src/bartcore/combiner.hpp#drawForestAmplitude]]).

The rescale-consistency set the move must carry, or a stored
`amplitude * leaf` stops being the identified product: `muByTree`, `totalFits`,
`totalTestFits`/`currTestFits` under `record`, and the keepTrees flattened slot
([[src/bartcore/combiner.hpp#rescaleAmplitudeRidge]]).

## The exponent rule

The rule is worked here at bcf's treatment forest (q = 2, `b0`, `b1` the two
glue scalars, `Lt`/`Mt` that forest's occupied-leaf count/squared sum,
`leafVar_tau` the leaf prior variance, `bPriorVariance` the fixed b prior
variance) - moved from docs/plans/archive/bcf-b-ridge.md secs 2.2-2.3 and 5a; its
Gibbs-block invariance argument (why rescaling `(b0, b1, tau)` on this orbit
is a legitimate sub-block update at all) stays there, sec 2.1, as BCF-specific
bookkeeping rather than exponent evidence. That argument ends at `b0`'s full
conditional given the orbit's invariants `(r = b1/b0, psi_l = b0 tau_l)`:

    q(b0) prop exp(-(b0^2+b1^2)/(2 sB2)) * exp(-S_psi/(2 leafVar_tau b0^2)) * |b0|^{1-Lt}.  (*)

**In terms of c (operational form).** Set `b0 = b0_0/c`, `c > 0`, `b0_0` =
current b0. Then `S_psi = b0_0^2 Mt`, `(1+r^2) b0_0^2 = b0_0^2 + b1_0^2`,
`db0 = -(b0_0/c^2) dc`. Substituting into (*) and folding constant `|b0_0|`
powers into the normalizer:

    exp(-(1+r^2)b0^2/(2 sB2))    -> exp(-(b0_0^2+b1_0^2)/(2 sB2 c^2))
    exp(-S_psi/(2 leafVar_tau b0^2)) -> exp(-Mt c^2/(2 leafVar_tau))
    |b0|^{1-Lt}                  -> |b0_0|^{1-Lt} c^{Lt-1}
    |db0/dc|                     -> |b0_0| c^{-2}

    q_c(c) prop c^{Lt-3} * exp( - Mt c^2/(2 leafVar_tau)
                                - (b0_0^2+b1_0^2)/(2 sB2 c^2) ).   (**)

The `c^{Lt-3}` (vs the a-move's `c^{L-2}`) is the operational fingerprint of
the second glue scalar: ONE fewer power of c. [The naive "move-map Jacobian"
shortcut -- `q(c) prop pi(T_c(x)) |det T_c'|` with `|det T_c'| = c^{Lt-2}` --
gives the WRONG `c^{Lt-2}`; it is off by one, exactly as it would be off by
one for the a-move (`c^{L-1}` vs the correct `c^{L-2}`). The prototype below
adjudicates decisively in favour of `c^{Lt-3}`.]

**Result: c^2 is Generalized Inverse Gaussian.** Substitute `v = c^2`
(`c^{Lt-3} dc = (1/2) v^{(Lt-4)/2} dv`):

    q_v(v) prop v^{(Lt-4)/2} exp(-alpha v - beta/v),
      alpha = Mt/(2 leafVar_tau),   beta = (b0_0^2+b1_0^2)/(2 sB2).

Matching the GIG density `prop v^{p-1} exp(-(A v + B/v)/2)`:

    -------------------------------------------------------------------------
    v = c^2 ~ GIG( p = (Lt-2)/2,  A = Mt/leafVar_tau,  B = (b0^2+b1^2)/bPriorVariance )
    -------------------------------------------------------------------------
      p = (Lt - 2) / 2
      A = Mt / leafVar_tau = Mt * (k/leaf.scale_tau)^2   (= Mt * P_tau)
      B = (b0^2 + b1^2) / bPriorVariance
    then  c = sqrt(v),  b0 <- b0/c,  b1 <- b1/c,  tau_l <- c * tau_l.

Contrast the a-move `GIG((L-1)/2, M/leafVar, a^2/aVariance)`. Two differences,
both structural and both prototype-confirmed:
  (i) p = (Lt-2)/2, NOT (Lt-1)/2 -- two glue scalars remove one more half-df
      than the one-scalar a-move. General rule: rescaling k leaf params by c
      against d glue scalars by 1/c gives p = (k-d)/2  (a: k=L,d=1; b: k=Lt,d=2).
  (ii) B = (b0^2+b1^2)/bPriorVariance with a FIXED prior variance -- no
      auxiliary, no conditioning, no lag. Both b coordinates enter B.

GIG generator: `ext_rng_simulateGeneralizedInverseGaussian(rng, p, A, B)`
ALREADY SHIPS (density `x^(p-1) exp(-(A x + B/x)/2)`; Dagpunar noshift
ratio-of-uniforms; the a-move added it, external/random.{h,c}). The b-move
reuses it verbatim -- no new RNG. `B=0 -> Gamma(p, rate A/2)`,
`A=0 -> inverse-gamma`. A single GIG draw covers every regime below.

**Pure-R prototype (adversarial check on the algebra) -- PASSED.**
An untracked prototype script, run and passed. Same logic as the a-memo:
the likelihood is constant along the orbit, so the move preserves the posterior
IFF it preserves the PRIOR's along-orbit conditional. Draw `(b0,b1,tau)` from
the prior (`b0,b1 ~ N(0,sB2)`, `tau_l ~ N(0,leafVar_tau)` iid), apply ONE move
(draw c from its conditional, rescale), test the pushed sample still has the
prior law (KS on marginals). Two independent parameterizations -- a vectorized
log-grid inverse-CDF on GIG(p=(Lt-2)/2, A=Mt/vT, B=(b0^2+b1^2)/sB2) via v=c^2,
and one on the operational form (**) in c directly. Results (N=20000):

    Lt=3  (p=0.5): GIG KS  b0=.113 b1=.061 tau1=.710 tauLast=.523 Mt=.340
                   oper KS b0=.051 b1=.208 tau1=.623   sign preserved TRUE
    Lt=8  (p=3.0): GIG KS  b0=.804 b1=.914 tau1=.900 tauLast=.497 Mt=.428
                   oper KS b0=.644 b1=.899 tau1=.995   sign preserved TRUE
    Lt=20 (p=9.0): GIG KS  b0=.785 b1=.887 tau1=.814 tauLast=.967 Mt=.752
                   oper KS b0=.773 b1=.978 tau1=.970   sign preserved TRUE

    DISCRIMINATION (wrong exponent -> prior NOT preserved, KS -> 0):
      p=(Lt-1)/2 [+0.5], Lt=3: tau1 KS = 1.6e-21 ; Lt=8: 2.7e-05
      p=Lt/2     [+1.0], Lt=8: tau1 KS = 7e-15
    combined-fit invariance: both arms 4.4e-16 ; all-treated 1.1e-16

All KS p-values non-significant across THREE Lt and TWO parameterizations
(cannot reject prior preservation); sign preserved; combined fit invariant to
machine precision including the all-treated edge. The discrimination arm shows
the check is SHARP: the off-by-one exponent p=(Lt-1)/2 (what the naive move-map
Jacobian gives) is rejected at KS=1.6e-21. This CONFIRMS p=(Lt-2)/2. The GIG
parameterization the implementer will code is the one validated. PASSES.

## ridgeB is code that is OFF

The b-move ships, but `AmplitudeSpec::ridgeB = false`
([[src/bartcore/combiner.hpp#AmplitudeSpec]]). Enabling it costs a GIG draw
per sweep, which re-records `bcf-equivalence`, and its own acceptance gate
([[docs/plans/archive/bcf-b-ridge.md:438-449@9cebb352]] - IACT payoff,
bcf-exact mode-2b, keepTrees round trip) was not run. It is a DOOR, not a fork:
it flips only on a named measured mixing case, plus that gate, plus a re-record
with the `equivalence.yaml` bump in the same commit.

**No silent enablement is possible** (resolved at M4.3 review,
[[docs/plans/multiforest-extension-surface.md:4967-4975@4c018187]]).
On every creation route the scale mixture holds if and only if a
forest is basis-FREE. R writes the forest's amplitude prior SCALE as 0 whenever
that forest declares a basis ([[R/model.R#forestParams]],
`if (withBasis) 0 else declared(spec$sd, 2)`); the bridge
reads it into `ForestSpec::amplitudePriorScale`
([[src/R_interface_bartcore.cpp#applyAmplitudeSpec]]) and derives
`forest.ridge = forest.amplitudePriorScale > 0.0`
([[src/R_interface_bartcore.cpp#applyAmplitudeSpec]]). So every
basis-carrying forest gets `ridge = false`, and the combiner's own
`halfCauchyScale` - the field `amplitudePriorScale` becomes - is zero there. The
one reachable q > 1 scale-mixture state, a post-creation widening of a
basis-free forest, consumes the SAME single GIG draw already taken at q = 1, at
the M4.2-validated exponent, so no stream moves.

## The basis is read, not classified

Nothing routes on basis SHAPE any more. `drawForestAmplitude`
([[src/bartcore/combiner.hpp#drawForestAmplitude]]) contracts forest f's design row - its basis row scaled by
its own fit - whatever the basis holds, so a continuous two-column pair is the
same conditional as an indicator pair rather than a different model.

That was not free before. The deleted two-scalar path never read `basis[1]` as
a DESIGN MATRIX: it borrowed that forest's column 1 and tested it for nonzero
as a GROUP KEY, never multiplying by the stored values, and formed two disjoint
group-precision accumulators from the partition. A legal continuous 0.25/0.75
pair would have been drawn as a different MODEL - both entries are nonzero, so
every row lands in the treated group - which is why a per-forest IS-CANONICAL
VALUE predicate, and not a width test, had to route it away.

The predicate is gone with the path it guarded: `AmplitudeState` carries no
canonical flag, `installForestBasis` recomputes nothing, and there is no
restore-side recompute left to explain as vacuous. What survives is the
discriminator the tests/cpp arm keeps: two width-2 bases sharing a treated
column and differing only in their control column must draw differently,
because the conditional contracts the whole row - where a three-amplitude
specialization keyed on the treated column alone would draw them
identically.

## One mutation route

Synthesis is CONSTRUCTION-ONLY: `synthesizeIndicatorBasis`
([[src/bartcore/combiner.hpp#synthesizeIndicatorBasis]]) is called once, from
the constructor ([[src/bartcore/combiner.hpp#AmplitudeForestCombiner]]).
`installForestBasis` is the SOLE mutator and therefore owns the guards nothing
else is left to apply - the index, `numColumns >= 1`, a non-null pointer, and
finiteness ([[src/bartcore/combiner.hpp#installForestBasis]]). It wins by
being the only operation there is, which is why no ordering between two
mutators has to be specified ([[src/bartcore/combiner.hpp#setForestBasis]]).

Ordering is LAST INSTALL WINS, per forest, and both orderings of a widen and a
swap collapse to it because `rebuildAmplitudeLayout` derives the offsets as a
pure prefix sum of the width vector and carries every block by position
([[src/bartcore/combiner.hpp#rebuildAmplitudeLayout]]). Amplitudes PRESERVE
and remap; surviving coordinates keep their values at their new offsets and
new ones enter at the neutral 1.0. A width-preserving install is the BITWISE
IDENTITY on every amplitude (the early return in
[[src/bartcore/combiner.hpp#rebuildAmplitudeLayout]]) - that is bcf's mid-life
z swap, and it is baseline-gating.

`glue_.z` is DELETED. The amplitude conditional contracts `basis[1]` itself,
not a borrowed indicator, because a width-preserving swap to a different
complementary pair would otherwise install the new pair, leave the layout
unmoved, and then draw `b0`/`b1` under the OLD partition while
`forestMultiplier` contracted the NEW basis (M4.3 item 1,
[[docs/plans/multiforest-extension-surface.md:1908-1928@4c018187]];
landed per
[[docs/plans/multiforest-extension-surface.md:4926-4927@4c018187]]).

## Persistence

The AMPLITUDES are state; the BASES are not.

`serializeGlue`/`restoreGlue` carry the per-forest widths, the flat amplitude
vector, and each forest's prior variance
([[src/bartcore/combiner.hpp#serializeGlue]]).
`glueIsValid` is the LAYOUT check `stateIsValid` routes through, because
`restoreGlue` writes THROUGH the live offsets and a state with the same total
over different widths would otherwise be admitted and silently permute the
blocks ([[src/bartcore/combiner.hpp#glueIsValid]], routed at
[[src/bartcore/chain.hpp#stateIsValid]]).

The four named scalars `a`/`aVariance`/`b0`/`b1` survive in `ChainStateData` as
a hand-written K = 2 READING only, non-authoritative, read exactly when
`amplitudeWidths` is empty - which is exactly a state written by hand rather
than by a combiner ([[src/bartcore/combiner.hpp#ChainStateData]],
[[src/bartcore/combiner.hpp#restoreGlue]]).

The bases ride CREATION, on `data@bases` (a LIST, [[R/A_class.R#dbartsData]],
validated by `validateForestBases`, [[R/data.R#validateForestBases]]), the way
the design matrix does. That is how RESTORE-THEN-WIDEN is met with no fourth
reapply hook: a widening applied after a restore preserves and remaps the
RESTORED amplitudes rather than the constructed ones
([[src/bartcore/combiner.hpp#serializeGlue]];
[[docs/plans/multiforest-extension-surface.md:1958-1978@4c018187]],
landed
[[docs/plans/multiforest-extension-surface.md:4931-4933@4c018187]]).

## Bitwise contracts

THREE accumulation directions, all observable once q > 1 or K > 2, each a
contract:

1. `combinedFits` accumulates WITHIN the row and from the LAST forest DOWN
   ([[src/bartcore/combiner.hpp#combinedFits]]). This is the
   load-bearing one and it is MEASURED. Under fused multiply-add contraction
   exactly one product in a sum escapes its own rounding - the one the closing
   add absorbs - and the two-forest `a mu + b_z tau` this replaces absorbed
   forest 0's. Seeding with the LAST forest leaves forest 0's product as the
   only bare multiply in the closing add. Accumulating FORWARD absorbs the last
   forest's product instead, moves ~30% of rows by one ulp and contaminates the
   trajectory within ~40 sweeps: all 12 `bcf-equivalence` scenarios red on mu,
   tau, glue, sigma and train ([[src/bartcore/combiner.hpp#combinedFits]];
   M4.1 landing note,
   [[docs/plans/multiforest-extension-surface.md:4798-4808@4c018187]]).
   `testCombinedFitsAssociation`
   ([[tests/cpp/test_sampler.cpp#testCombinedFitsAssociation]]) is the ONLY
   in-process guard - the M4.0 seam pin structurally CANNOT see association,
   its reference expression inheriting the test compiler's own contraction.
2. `formForestResponse`'s residual accumulates FORWARD, subtracting the other
   forests in increasing index order, and deliberately NOT combinedFits'
   reverse: this sum has no two-term fused expression to reproduce, and the
   amplitude conditional forms the same residual the same way, so the two agree
   by construction rather than by coincidence
   ([[src/bartcore/combiner.hpp#formForestResponse, drawForestAmplitude]]).
3. `forestMultiplier` contracts FORWARD over the columns; at q > 2 a
   reassociation moves the multiplier by an ulp and every reader of it with it
   ([[src/bartcore/combiner.hpp#forestMultiplier]]).

The standing lesson these pins carry, twice learned: a pin fixture must give
every factor in the pinned expression a DISCRIMINATING value - unit values
silently vacate pins (M4.0's required fix,
[[docs/plans/multiforest-extension-surface.md:4745-4753@4c018187]];
recurred at M4.3's Arm 5(iii) dead pin,
[[docs/plans/multiforest-extension-surface.md:4957-4960@4c018187]]).

## The calibration map, general in K

Stated on `ForestSpec` rather than in the chain
([[src/bartcore/combiner.hpp#ForestSpec]], applied
[[src/bartcore/chain.hpp#Chain]]). The node scale is
`nodeScaleFactor * s / (nodeScaleDivisor * c)`, s the family's own LATENT
SCALE and c the median nonzero row norm of the forest's basis - exactly `1.0`
on every shipped route, so c is inert there:

| family | s |
|---|---|
| gaussian | sample sd of the range-scaled response |
| probit | 1 |
| logistic | pi/sqrt(3) |

Two branches, which are the two ways a magnitude can be carried. A
SCALE-MIXTURE amplitude carries it in `amplitudePriorScale` and leaves the node
scale at s - bcf's prognostic forest, factor and divisor both 1, giving
`mu ~ N(0, s^2)` with the half-Cauchy a at `aPriorScale` units of s. A
FIXED-VARIANCE amplitude carries it in the node scale, divided through the
half-normal median 0.674 - bcf's treatment forest, giving
`tau ~ N(0, (sdModerate s / 0.674)^2)`, so `(b1 - b0) tau` sits at
`sdModerate` units of s. Under gaussian s IS `sd(y)`; under probit and
logistic s is the link's own fixed latent scale, sigma is PINNED, and there is
no response sd for either quantity to be stated against.

The divisor is CARRIED rather than derived from the prior, so a forest
declaring a zero scale keeps the node scale its host asked for. Both
expressions are written exactly as bcf's were, which is what keeps the K = 2
instance bitwise.

**Two further facts about the shipped prior, both load-bearing.** Adaptivity
is capped at one forest, for any K: `resolveForests`
([[R/model.R#resolveForests]], refusal also in
[[R/model.R#"needs a 'basis': the amplitudes multiplying it"]]) requires every forest past the first to carry a
basis, and `forestParams` writes the LITERAL `0` for `amplitudePriorScale`
whenever a basis is present ([[R/model.R#forestParams]]), from which the
bridge derives `forest.ridge = false`
([[src/R_interface_bartcore.cpp#applyAmplitudeSpec]]) - forests 2..K are
ALWAYS fixed-variance. There is also NO per-K renormalization anywhere in the
MAP, and `binary-kforest-prior-default` S2 added none: the map still disperses
as `sqrt(K)` by construction (exponent exactly 1/2), `1.04912 sqrt(K) s` at
unit row norms. What that slice moved is the DEFAULT `nodeScaleFactor`
(`forestParams`, R/model.R), from the literal 1 to `sqrt(2/K)`, whose product
with the map's own `sqrt(K)` is `sqrt(2)` at every count. So the shipped
all-basis index prior is `1.4837 s` at every K rather than `1.4837 s` at K = 2
growing to `2.9674 s` at K = 8, and the argument is compositionality rather
than magnitude: without it the prior on the combined location depends on how
the caller DECOMPOSED the mean rather than on what they said about it, the same
model written as two forests or as four differing by `sqrt(2)`.

**On the SHIPPED shape the law caps rather than pins, and no sentence anywhere
may state the stronger form.** With one basis-free forest and K - 1 basis
forests the fixed-variance channel under the default is
`1.049120 sqrt(K-1) sqrt(2/K) s` = `sqrt((K-1)/K) s / 0.674`: 1.0491199,
1.2114193, 1.2849042, 1.3878551 at K = 2, 3, 4, 8, which is 0.699, 0.808,
0.857, 0.925 of the `k = 2` binary index budget `1.5 s`, rising monotonically
toward the limit 0.989 and never reaching it. The reason is structural - a
Cauchy channel has no finite variance to enter a root-sum-square budget with,
so a shape carrying one cannot have its total pinned by any per-forest law -
and it is the same fact that killed index normalization. Against today's
unbounded `1.049120 sqrt(K-1) s`, which is 1.850x the budget at K = 8 and
passes 2x at K = 10.

**The contract, in the form the row-norm fix makes true, identical under
gaussian and under a latent family.** (1) `sd` and `amplitude.prior.variance`
are stated PER UNIT OF BASIS ROW NORM: a forest whose basis rows have median
nonzero norm c contributes the scale named, and the map divides c out.
(2) The induced prior sd on the index at row i is
`sqrt( sum_f prior.scale_f^2 v_f ||B_f(i,.)||^2 )` over the fixed-variance
forests, `prior.scale_f` read from `$getCalibration(f)[, "prior.scale"]` and
`v_f` the forest's `amplitude.prior.variance` (default 0.5); a basis-free
forest's own contribution is Cauchy with no sd. (3) Under probit and logistic
that index is in LATENT sd units and sigma is PINNED, so nothing absorbs a
mis-scaled basis; under gaussian it is in `sd(y)` units and a drawn sigma
partly does. (4) The MAP disperses the index as `sqrt(K)`, per the two facts
above, and the shipped DEFAULT `sd` cancels that: the all-basis index prior is
`1.4837 s` at every K, `0.989` of the `k = 2` binary budget, and on the
one-basis-free shape the fixed-variance channel is BOUNDED by it, rising from
`0.699` of that budget rather than pinned. A declared `sd` opts out and keeps
its per-forest reading at every K. (5) At the shipped defaults a K = 2 binary
K-forest prior puts `P(p < 0.01 or p > 0.99)` at 0.2468 under probit, against
0.2387 for the shipped single-forest binary default (`chi(1.5, 2)`), where
before `binary-kforest-prior-default` S2 it was 0.3764 - M4.4's documentation
mandate, discharged here in its successor form.

## bcf as the K = 2 instance

`expandForestSpecs` ([[src/bartcore/combiner.hpp#expandForestSpecs]]) is the
thin adapter between bcf's two-forest spelling and the K-length vector every
other layer works in, and it is LOAD-BEARING rather than courtesy: 25
`tests/cpp` fixtures and benchmarks/R/bcf-equivalence.R drive through the
`AmplitudeSpec` spelling
([[docs/plans/multiforest-extension-surface.md:1785-1791@4c018187]],
[[docs/plans/multiforest-extension-surface.md:2182-2184@4c018187]]).
Forest 0 takes the half-Cauchy amplitude over the implicit intercept and
leaves its node scale at s; forest 1 takes the fixed-variance pair over the
treatment indicator basis and carries `sdModerate` in its node scale.

`bcfGlue(a, b0, b1)` survives as a named READING that returns false on any
other layout, which is how a caller learns it is not looking at bcf
([[src/bartcore/combiner.hpp#bcfGlue]], [[src/bartcore/chain.hpp#forestTotalFits]]).

## Surfaces

The map is READABLE, at all three layers, and is the only route to any of its
five quantities (`binary-kforest-prior-default` S1):
`$getCalibration(f)` reports `amplitude.prior.variance` and
`amplitude.prior.scale` - exclusive per forest, the fixed-variance and
scale-mixture spellings of `ForestAmplitudePrior` - beside `node.scale.factor`,
`node.scale.divisor` and `basis.row.norm`, NaN on any forest with no map entry;
`bartcore_getCalibration` carries the five columns and
`dbarts_forest_calibration` the five appended fields. The anchor s has no
column and is recovered as `prior.scale * divisor * rowNorm / factor` whenever
`node.scale.factor` is not NaN, which is exactly when the calibration in force
is the map's: a `setState` or `installTrees` that brings a foreign leaf scale
clears both `node.scale` columns (the amplitude prior follows the state, the
row norm is unaffected - bases are not state) and `setForestBasis` re-imposes
the map and restores them.

R creation: `forests = list(forest(basis = ...))` plus the per-forest knob map
(`resolveForests`, R/model.R), and `dbartsData(bases = )` for a numeric basis
([[R/spec.R#resolveSamplerSpec]], the `validateForestBases` install; the
family gate it feeds is also in [[R/spec.R#resolveSamplerSpec]]). R5:
`$setForestBasis(forest, basis)` ([[R/dbarts.R#dbartsSampler$setForestBasis]])
and `$getForestAmplitudes(forest)`
([[R/dbarts.R#dbartsSampler$getForestAmplitudes]]), both 1-based via
`resolveForestIndex` ([[R/bartcore.R#resolveForestIndex]]). Flat C:
`dbarts_sampler_setForestBasis`, `dbarts_sampler_numForestAmplitudes` and
`dbarts_sampler_getForestAmplitudes`, ragged and ROW-major
([[inst/include/dbarts/dbarts.h#dbarts_sampler_setForestBasis, dbarts_sampler_numForestAmplitudes, dbarts_sampler_getForestAmplitudes]]).

Capability probes are `totalAmplitudes() != 0`, NEVER a forest count: a
K-forest multinomial defeats a `numForests` test, and the shipped code states
that rule in two independent places
([[src/C_interface.cpp#dbarts_sampler_numForests]],
[[src/R_interface_bartcore.cpp#bartcore_setForestBasis]]).

Per-forest SPLIT COUNTS are reported too (bcf-bartcause-relocation D3): the
combiner overrides `numVariableCountForests` to its own forest count and
`variableCountForest(j)` to `j`, so `run()$varcount` is
`n.predictors x n.forests x n.samples x n.chains`, forest-major and prognostic
first - the multinomial coupling's channel, at a coupling whose
`numReportedLocations` stays 1, which is exactly why the two axes are keyed
separately. The width a run actually writes is the CALLER's declared
`Results::numVariableCountForests`, clamped once to this count in
`Sampler::run`; a caller declaring one (the flat C API, rbart_vi's callback
loop) gets slot 0, the reported forest, byte for byte as before.

## What this family does not do

- `aft`, `ordinal` and `nbinom` responses - refused at creation, naming what
  each is missing (gaussian, probit and logistic landed at M4.4).
- The test surface, for EVERY family including the two M4.4 added.
  `testFitsAreDefined` and `logLikelihoodIsDefined` are both false
  ([[src/bartcore/combiner.hpp#testFitsAreDefined, logLikelihoodIsDefined]]),
  so `setTestPredictors`, `setTestOffset` and `predict` refuse - through
  `refuseUndefinedTestFits`, gated on `testFitsAreDefined` rather than on the
  forest count
  ([[src/R_interface_bartcore_common.hpp#refuseUndefinedTestFits, testFitsAreDefined]])
  - and no log-likelihood is reported. Unchanged by M4.4.
- A per-draw amplitude channel in flat C. `dbarts_results` carries none
  ([[inst/include/dbarts/dbarts.h#dbarts_results]]); DECLINED at
  [[docs/plans/multiforest-extension-surface.md:1521-1531@4c018187]],
  a `DBARTS_C_API_MINOR` bump binding decision 8 forbids.
- A variance forest. `createAmplitudeSampler` refuses `numVarianceTrees > 0`
  ([[src/bartcore/facade.hpp#createAmplitudeSampler]]).
- Nameable leaf-prior calibration. The map owns it, so the write is refused on
  ANY combining sampler: `Chain::setForestPriorScale` returns `false` on
  `f >= forests_.size() || combiner_ != nullptr`
  ([[src/bartcore/chain.hpp#setForestPriorScale]]), which the
  flat entry surfaces to a caller as a 0 return.

The classes the family ENABLES, and their evidence status, are
[[docs/design/model-space-survey.md#D4. The general per-forest multiplier]]:
continuous/dose-response exposure BCF (unpublished preprint), VCBART's
varying coefficients (which carries no sampled amplitude), heterogeneous
mediation, and the multiplier half of principal stratification with BCF.

## Landing notes

**M4.0, eb340f4f (pins).** Component pins over BOTH `afterCombine` overrides -
BCF's GIG rescale with its three reachable 1.0 skips inert on the rng stream,
and the multinomial additive shift with its returns-1.0-having-moved
convention, whose pin is the SOLE guard on that convention. The review's
required fix became the slice's lesson: the fixture originally set
`numTrees = 1` on every forest, collapsing each per-tree factor to 1 and
leaving four sites undiscriminated, with the reviewer MEASURING a wrong
implementation staying green. Fixed with unequal per-forest tree counts. ~348
dense-equivalent of the TESTS band - the slice added no engine code at all, so
this figure is not comparable with the engine nets quoted for M4.1-M4.3 below.

**M4.1, e3170da4 + follow-up 9459bef0 (the multiplier generalization).**
`forestMultiplier` became `dot(a_f, B_f(i,.))` and `combinedFits` a K loop with
no K = 2 special case. The defining event was a 12/12 `bcf-equivalence`
divergence root-caused to FMA CONTRACTION ASSOCIATION, fixed by accumulating
from the last forest down (see "Bitwise contracts"). Engine +48 net
dense-equivalent, tests +98. The reviewer's trap note, standing: a perturbed
install WITHOUT `--preclean` reported bitwise identical, so every re-check must
`--preclean`. BENCH, on a granted quiet window: B/A per-sweep 1.0098 at
n = 20000 and 1.0105 at n = 2000, flat across 100x in n - per-element compute
in the combiner, not bandwidth - and the ~1% BCF-only cost was ACCEPTED as the
price of the general multiplier
([[docs/plans/multiforest-extension-surface.md:4834-4846@4c018187]]).

**M4.2, a35ff7df (the q-variate conditional).** The per-forest q-variate
conditional through the new unit-lower LDL' helper, and ONE general per-forest
ASIS rescale at `p = (L - q)/2` whose q = 1 instance is bitwise-identical to
bcf's a-move. The in-slice rule's outcome: bitwise on the ASIS half,
SPECIALIZED-PATH-KEPT on the conditional half, with impossibility CONFIRMED by
the reviewer over 21 accumulation variants against the real engine (a first
standalone probe false-alarmed - vectorization split the multiplies in the
replica, a recorded methodology lesson). Engine ~180 dense-equivalent.

**M4.3, 983d7f0a (SUBSUME).** `setTreatment` retired as a mutator at all four
layers in favour of `setForestBasis(f, values, q)`; `installTreatment` became
the constructor-only `synthesizeIndicatorBasis`; `installForestBasis` became
the sole mutator; `glue_.z` DELETED; draw-path selection moved onto the
per-forest canonical VALUE predicate; `bcfGlue` re-signed to
`totalAmplitudes`/`numForestAmplitudes`/`amplitudes` with a K = 2 thin adapter;
`AmplitudeSpec` gained its K-length `forests` vector with `expandForestSpecs`,
all 25 AmplitudeSpec fixtures untouched; `data@treatment` retired for a
`data@bases` LIST riding creation. Refusal ledger as the review re-derived it:
5 relax, 3 restate, 4 rewritten by the slot retirement, 33 stand, and ONE
dropped outright with licensed successors (the treatment-coding refusal,
replaced by length/finiteness refusals at creation). `consumer.c` `LEG_COUNT`
18 -> 19. Engine ~165 dense-equivalent.

**M4.4, e5e93f11 (the latent family).** `probit` and `logistic` wired
through the K-forest constructor's response switch
([[src/bartcore/chain.hpp:747-762@e5e93f11]]), the
calibration map re-based on `latentScaleAnchor` at the settled Option L
anchor - probit s = 1, logistic s = pi/sqrt(3) - and closed by option E, the
basis row-norm divisor (median of nonzero row norms, re-derived on
`$setForestBasis` rather than left stale so a mid-sample swap cannot go
stale). `$getSumsOfSquaredResiduals` fixed to read the COMBINED location
rather than forest 0's bare total, per-family rather than by blind
substitution (a multinomial's slab is K x n, not length n); draw-neutral,
argued from full buffer overwrite, no in-sweep caller and no rng in the path,
not from the getter's name. M4.4's own gate: FA5 arm E, the K-forest probit
sampler, AGREES with leg P's INDEPENDENT reference - approximate by
construction, carrying its own convergence gate at max R-hat 1.0003 against
1.01, and never called exact - at max |z| = 1.46, and with decisive arm B at
max |z| = 1.73 (threshold 3.0 on both, 12 functionals,
the FA5 harness). The pinned sigma
costs no measurable mixing: arm E's IACT is at or below arm B's on 8 of 12
functionals. The anchor's own gate (checklist item 25, a prior-predictive
exactness check) was SUBSTANTIATED by mutation builds rather than asserted -
GREEN as shipped at 76/76, and RED under the naive cold-start anchor (14
assertions move), under the rejected 3x-wider anchor in both families (12),
and under a dropped row-norm divisor (4), each failing on the assertion built
to catch it. That last arm's coverage rests on ONE planted fixture cell whose
basis row norm is 4; remove it and the mutation goes green, so do not tidy it
away. Equivalence trio bitwise on all three suites (bcf-equivalence,
equivalence, multinomial-equivalence), no baseline re-records. Engine 67,
bridge 33, R 27 dense-equivalent (the bridge's 42 is its RAW net, and
reporting raw as dense would have put it at its own stop).

**binary-kforest-prior-default S1 (the map becomes readable).** Not an M4
slice, but it lands on the map this file owns, so its note is here.
`ForestCalibration` gained five fields - the two exclusive amplitude-prior
spellings, the node scale's factor and divisor, and the basis row norm the
constructor already computed and discarded - surfaced as five
`$getCalibration` columns and five appended `dbarts_forest_calibration`
fields. DRAW-NEUTRAL by construction and confirmed: the engine echoes its own
retained SPEC rather than reading the combiner, so no virtual joined
`ForestCombiner<L, ResidT>` and the vtable did not move, and the leaf-scale
expression keeps its written association with the row norm merely NAMED where
it was previously an argument. Equivalence trio bitwise on all three
baselines. The one non-echo is Fork 3b's truthfulness rule at the two state
install sites: a donor leaf scale differing BITWISE from the one in force
clears a per-forest `nodeScaleIsMapDerived_` flag, so the two `node.scale`
columns read NaN rather than describing a decomposition of a number the chain
no longer runs under, while the amplitude prior FOLLOWS the state under
`restoreGlue`'s own guard (the existing `glueIsValid` virtual, no new one) and
`setForestBasis` re-imposes the map. `sizeof(dbarts_forest_calibration)` moved,
so a `LinkingTo` consumer recompiles; `dbarts_apiHash()` did not, the append
being below the header's documented 1.0-0 boundary.

**Budget units, as a standing convention**
([[docs/plans/multiforest-extension-surface.md:4985-4989@4c018187]]).
Slice bands on
this arc are DENSE-EQUIVALENT lines - lines counted without blank-line and
formatting inflation, which roughly doubles raw counts. M4.3's implementer
reported raw nets against a dense band and appeared over budget when it was
85-135 UNDER. Never mix the two currencies in one report.

## Status

LANDED, gaussian, probit and logistic. M4.4 discharged the header's former
"Gaussian responses only" - [[inst/include/dbarts/dbarts.h:717-718@e5e93f11]]
now names gaussian, probit and logistic, with aft, ordinal and nbinom refused
by name at creation. The ridgeB door is shut. The naming debt is discharged;
see above.
