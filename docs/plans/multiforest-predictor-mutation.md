# multiforest-predictor-mutation

agent: SL (stop-loss) sonnet, standalone, lands out of band. S0 sonnet (pins,
  harness, baselines). S1, S2, S3 opus (each moves a live path and each can
  strand tree state). S4 sonnet (docs, records, residual coverage). Serialized:
  one implementer, each slice lands before the next starts.
rng: SL, S0, S4 NEUTRAL. S1, S2, S3 NEUTRAL on every path that exists today -
  the widened paths are refused at the bridge, so they are NEW streams with no
  baseline. A trio deviation is a LEAK and aborts the slice; it is never a
  re-record. S3 is additionally NEUTRAL on the heteroscedastic FORCED and
  setData paths, which the refactor must recompose byte-for-byte.
window: bcf-public-surface goes FIRST. VD 2026-08-09, recorded in TODO:72-73:
  "goes FIRST, ahead of multiforest-predictor-mutation." Not reopened here; the
  only interleave is SL, a ~155-line guard commit that must precede
  bcf-public-surface S1 (see "The live defect, and the stop-loss"). S0-S4 run
  after bcf-public-surface lands.
budget: SL ~35 bridge/C + ~120 test. S0 ~120 harness R + ~190 tinytest + three
  recorded baselines. S1 ~95 engine + ~25 bridge + ~265 test. S2 ~135 engine +
  ~20 bridge + ~285 test. S3 ~160 engine + ~15 bridge + ~250 test. S4 ~95 docs
  + ~60 test. Total ~1770.
tip: every anchor below re-read against the live tree at fed324e (branch
  bartcore). The zero-weight S2 weight channel landed at 153d1dd and moved
  chain.hpp and facade.hpp; the design memo's anchors (f05b604) and the blind
  critique's corrections to them are BOTH stale in places and are superseded by
  "Verified seams".

## Goal

The transactional predictor surface stops being single-forest-only. A BCF, a
multinomial and a heteroscedastic sampler accept `setPredictor` /
`updatePredictor` with rollback and the per-observation update session, on the
same acceptance criterion the single-forest sampler already ships - no leaf
empties in any tree of any ensemble of any chain - with the per-forest column
mask as the exact opt-out. `refuseMultiForestTransactionalUpdate` and
`refuseVarianceForestPredictorMutation` retire from the four transactional
entries. No `dbarts.h` symbol is added or changed and no API hash moves.

## Binding decisions inherited (do not reopen)

1. **Sequencing.** bcf-public-surface first (VD 2026-08-09, TODO:72-73).
2. **Authorization.** `docs/plans/multiforest-veto-rate-falsifier.md` ran
   2026-08-09, verdict YELLOW on both column types, final. The widening costs at
   most +0.0369 pp of install rate against a 2 pp GREEN line; the frozen-
   subpopulation share is 0.000 pp in all 18 gated cells; M1 (mask exactness)
   was exactly 0 over 14,400,000 install decisions. YELLOW obligates three
   things, discharged here: mask and j-split pruning ("Pruning"), a priced
   comparison against the alternatives ("Alternatives priced"), and the P4/T2
   attribution table with its caveat ("The attribution caveat").
3. **The veto is per SAMPLER, not per forest** (the falsifier's open fork 2,
   closed here; see "Fork 2"). The per-forest column mask is the opt-out.
4. **E1, the widened veto, is the semantics.** E2 (collapse-on-object) is
   rejected on record. E4 (rebuild the sampler) is priced and corrected below;
   it is NOT strictly dominant for the variance forest.
5. **Occupancy is count-based** (zero-weight-exactness binding decision 2). A
   per-forest weight, including an exact zero, changes no veto and no leaf's
   occupancy. Nothing in this arc makes per-forest zero weight an opt-out.
5b. **The large-K price is UNPRICED and is recorded, not assumed away**
   (multiforest-extension-surface, 2026-08-10). The veto rate scales with the
   j-splitting trees summed over ENSEMBLES, and the per-forest column mask
   does not save a VCBART-shaped model, where all p+1 ensembles are functions
   of the SAME modifiers, so the mask is vacuous there rather than an opt-out.
   No measurement exists above K = 4. If a general K-forest family lands, that
   arc owns a fresh veto-rate measurement at its own K.
6. **Consumer compatibility is a cost to enumerate, never a design constraint**
   (VD scope ruling 2026-08-10): everything is pre-release and the sister
   packages migrate in lockstep.

## Adjudication (design memo vs. blind critique)

The design cycle ran memo -> blind critique -> this synthesis; both artifacts
live under `.claude/`, which is gitignored, so every load-bearing fact is
carried here. The critique's verdict was STANDS WITH AMENDMENTS: the
architecture holds, the gate battery did not. Six blocking findings, adjudicated:

1. **The pruning-exactness falsifier was self-contradictory. UPHELD, repaired
   differently.** The memo forbade pruning forest 0's rebuild because
   `(T - f) + f != T` in IEEE, then demanded bitwise-identical draws between a
   pruned and an unpruned build while pruning `f >= 1` of the same loop.
   Measured: 3.1% of 2,000,000 subtract-then-add round trips at the loop's real
   magnitudes are inexact, max residual 4.44e-16; and in-engine,
   `forceRefreshTrees` versus `revalidateAllChains` on a homoscedastic sampler
   differ by 4.49e-14 on `train` from reassembly order alone. The critique's
   repair (declare the pruned build the record, compare partitions instead of
   draws) is REJECTED as unnecessarily weak. The repair taken here: **there is no
   second build.** For `f >= 1` the pruned path is the arithmetically exact one -
   an untouched tree's fits and partition are not recomputed at all, so its
   contribution to `forest.totalFits` is preserved bitwise, where the unpruned
   loop would perturb it by a ULP. The gates become F6 (untouched-tree bitwise
   invariance, asserted on the shipped build) and F5 (mask equality against an
   in-test oracle). Neither needs a build flag and both can go red.
2. **The mask oracle was not constructible. UPHELD.** The install mask is
   scan-order dependent (`commitObservation` mutates `leafCounts_`, which
   `observationWouldRemainValid` reads) and the order is drawn inside the engine
   by `ext_rng_drawPermutation(chains_[0]->rng(), ...)` (`sampler.hpp:1066`); no
   R oracle can reproduce it. The falsifier's own mask gate was distributional
   for exactly this reason. Replaced by an in-`tests/cpp` oracle driven through
   the public session API at a CALLER-CHOSEN scan order
   (`Sampler::beginPredictorUpdate`, `sampler.hpp:1077`, which the fuzzer already
   drives row by row at `test_fuzz.cpp:355-362`). This is exact, deterministic,
   and drops the arc's dependency on the gitignored falsifier harness entirely.
3. **The rollback and session-abandon gates were weak. UPHELD, closed
   structurally.** `FuzzSnapshot` (`test_fuzz.cpp:33-38`) carries codes, cuts,
   sigma and `Chain::treeFits()`, which is forest 0 only (`chain.hpp:2795-2797`).
   A rollback that restored mu and left tau routed by the proposal passes it
   green. Rather than adding tau and variance fields by hand - the
   statesAgree-omits-a-field pattern that has bitten this program before - S1
   rebuilds the snapshot so that omission is not expressible: see "The snapshot,
   closed structurally".
4. **E4's strict dominance is false for the variance forest. UPHELD, and it cuts
   the other way.** `Chain::stateIsValid`'s variance branch
   (`chain.hpp:2369-2380` at fed324e) checks tree count, flat well-formedness and
   strict leaf positivity and NOT occupancy; `setState` installs through
   `rebuildVarianceForest` (`:3269-3296`), which repartitions and scatters with
   no occupancy test. The occupancy check at `:2704` lives in
   `installVarianceForest`, which `setState` does not use. So for a
   heteroscedastic sampler the storeState / new sampler / setState route is not
   the same criterion at higher cost: it is a WEAKER criterion that admits
   variance states the veto refuses, and `refreshVarianceForest`'s occupancy
   assertion (`:3352`) is compiled out in a release build. That makes the refusal
   message's current advice ("make a new sampler instead") a pointer at an
   unchecked route, which strengthens the case for lifting the variance forest
   rather than weakening it. S3 closes the gap (item 5) so the dominance claim
   becomes true before the message retires.
5. **The window contradicted a settled ruling. UPHELD.** Not reopened; see
   "window" above. The critique's second half is also upheld: bcf-public-surface
   S0 (`docs/plans/bcf-public-surface.md:286-289`) pins the BCF refusals this
   arc's S1/S2 retire. Ownership is settled in S0 item 1.
6. **Fork 2's closure rested on a mis-stated premise. UPHELD.** Restated in
   "Fork 2". The conclusion (per-sampler veto) survives; the argument changes.

Advisories accepted: the multinomial margin is one and a half to two orders
under the GREEN lines, not three (binary per-observation `E_all` 0.087% against a
matched 0.0279% baseline is D = 0.059 pp against a 2 pp line, a factor of 34; the
largest margin anywhere is P1's 15 pp KILL line, a factor of 254); the memo's
attribution of Q1's "frozen forest" arm to a fresh sampler is wrong (Q1 froze the
trees, a fresh sampler regrows them - different kernels, so Q1's 0.488 sd prices
neither); "the session alone is not constructible" overstates (it requires the
whole-matrix slice's engine change - say that); the pruned session cache needs an
explicit survivor table, not offset arithmetic; the pruned multinomial tree count
is ~105, not ~115; `test_model.cpp:718-733` already asserts the E = 3 factory
refusals, so that is a standing assertion and not new work; the falsifier harness
is under `.claude/` and is not a repository artifact.

Advisory OVERTURNED: the claimed `dropStaleMissingDirections` ordering hazard is
not a hazard. `Tree::dropStaleMissingDirectionsBelow` (`tree.hpp:1213-1223`)
clears the bit only when `!data.hasMissing[j] && !data.columnIsPooled(j)`, i.e.
only when no observation is missing on that column, so the bit routes nothing and
clearing it moves no observation. The comment at `tree.hpp:1062-1067` states this
and the code confirms it. The drop may therefore sit on either side of the
validate phase's repartition; it matters only to `buildFromFlat`'s gauge check.
F10's negative half does not need to cover it.

## Verified seams (read at fed324e)

**The refusals.** `refuseVarianceForestPredictorMutation`
(`src/R_interface_bartcore.cpp:2118-2124`) keys on `shape().hasVarianceForest`;
`refuseMultiForestTransactionalUpdate` (`:2136-2144`) calls it when
`!forcedUpdate`, then errors on `numForests >= 2` when `!forcedUpdate`. Both live
in that file's ANONYMOUS namespace (`:35-2222`), unlike `refuseMultiForestMutation`
(`:2245`) and `refusePinnedSigmaChange`, which sit in `bartcore_bridge` and are
declared in `R_interface_bartcore_common.hpp:115,122`. Four call sites:
`bartcore_setPredictor` (`:3945`, forwards force), `bartcore_updatePredictor`
(`:3991`, forwards force), `bartcore_updatePredictorPerObservation` (`:4099`, NO
force argument), `bartcore_updatePredictorPerObservationJointly` (`:4164`, per
sampler). `Chain::numForests` is `forests_.size()`; a heteroscedastic sampler
reports `numForests == 1` and holds its variance trees in `varianceForest_`. The
unit of study is `E = numForests + (hasVarianceForest ? 1 : 0)`.

**The two-phase transaction.** `Sampler::runPredictorTransaction`
(`sampler.hpp:1334-1368`): cut-validity precheck; then either `applyForced` plus
`chain->forceRefreshTrees()`, or snapshot-apply plus `revalidateAllChains()`,
restoring and calling `chain->repartitionTrees()` on failure.
`Sampler::revalidateAllChains` (`:1195-1207`) validates every chain, then rebuilds
every chain, so a late chain's rejection never leaves an early chain's fits
overwritten. `Chain::revalidateTrees` (`chain.hpp:1571-1581`) opens
`Forest& forest = forests_[0]`, recovers leaf parameters, `repartitionSubtree`,
AND-folds `bottomNodesAreOccupied()` with early exit.
`Chain::rebuildFitsFromParameters` (`:1589-1612`) is forest 0 only:
`dropStaleMissingDirections()` (which DOES loop all forests, `:1733-1740`), then
per tree `subtractTreeFitsFromTotal` / `setTreeFits` /
`installLeafOfAndAddToTotal`. `Chain::repartitionTrees` (`:1616-1624`), the
rollback re-route, ALREADY loops every forest; it does not touch
`varianceForest_`. `Chain::forceRefreshTrees` (`:1743-1789`) loops every forest,
repartitions AND collapses, then calls `refreshVarianceForest(nullptr)` at
`:1789`.

**`totalFits` is per forest**, not per chain: `addTreeFitsToTotal` and
`subtractTreeFitsFromTotal` (`chain.hpp:3666-3691`) read and write
`forest.totalFits`, and for the constant leaf they gather `mu[leafOf[i]]`. The
BCF combiner consumes each forest's total separately. Widening the rebuild loop
is therefore a clean per-forest operation with no cross-forest arithmetic.

**`refreshVarianceForest`** (`chain.hpp:3337-3383`): per tree, recover
node-indexed factors through the LIVE partition (or take them from a caller's
earlier recovery), drop stale missing directions, optionally remap onto a new cut
grid, optionally re-point at a resized index buffer, repartition,
`collapseEmptyNodes<GeometricMerge>`, scatter surviving factors through the new
partition; then recompute `combinedVariance` as the per-row product. It has no
validate/rebuild split and it collapses rather than reporting. Its comment
records that recover-before-repartition is load-bearing.

**The session.** `Sampler::UpdateSessionImpl` (`sampler.hpp:1382-1459`):
`totalNumTrees = options_.numTrees * chains_.size()`; `observationLeaf_` is
`totalNumTrees * n` int32; `leafCounts_` is per tree, arena-id indexed;
`treeAt(t)` is `chains_[t / numTrees]->tree(t % numTrees)` (`:1446-1449`) and
`Chain::tree` is `forests_[0].trees[t]` (`chain.hpp:2792`). For a BCF,
`options_.numTrees` is mu's count alone and tau's differs.
`observationWouldRemainValid` iterates `leafCounts_.size()` with early exit and
descends with the candidate code overriding `column_`; `commitObservation` writes
exactly one cell (`data_.setCell`) plus count bookkeeping. There is no journal and
no abandon-rollback. `Sampler::updatePredictorPerObservation` (`:1060-1075`) draws
the scan permutation from `chains_[0]->rng()`, and `finalize()` IS
`revalidateAllChains()` (`:1443`). `updatePredictorPerObservationJointly`
(`facade.hpp:482`) runs one session per sampler on one shared permutation from
`samplers[0]->rng()` and installs in every sampler or none.

**The pruning primitives.** `Forest::columnMask` is per forest and pushed onto
every tree; `VarianceForest::columnMask` likewise.
`Tree::countVariableUses(uint32_t*)` (`tree.hpp:1072-1074`, recursive helper at
`:1213`) is an O(nodes) per-tree variable-use census. That is the j-split test and
it already exists.

**The criterion is not new.** `Chain::stateIsValid` loops every MEAN forest,
builds a scratch tree per stored tree, repartitions, and returns false unless
`bottomNodesAreOccupied()` (`chain.hpp:2324`). So "no leaf empties in any tree of
any mean forest of this chain, against the current store" is ALREADY the shipped
acceptance criterion of `setState` and `installForests`. This arc makes the
transactional path agree with the restore path; it does not invent a criterion.
The VARIANCE branch (`:2369-2380`) is the exception adjudication 4 turns on.

**Accessors that already exist** and this arc reuses: `Chain::numTreesInForest`
(`chain.hpp:719`), `Chain::numVarianceTrees` (`:733`), `Chain::forestTreeFits(f,
out)` (`:885`), `Chain::treeInForestForTesting(f, t)` (`:2830`),
`Chain::varianceTreeForTesting(j)` (`:2838`), `Chain::varianceFactorsForTesting`
(`:2841`), `Sampler::numTreesInForest(f)` (`sampler.hpp:1163`),
`Sampler::numForests()` (`:1118`), `Sampler::getState(SamplerStateData&)`
(`:588`). Only `Chain::totalFits()` (`chain.hpp:2800`) is forest-0 pinned and
needs a per-forest sibling.

**Factories.** `createBCFSampler` returns nullptr on `numVarianceTrees > 0` at
`facade.hpp:667`; `createMultinomialSampler` at `:680`. (The memo cited 631/644,
the critique corrected to 651/664; both are stale at fed324e - the weight channel
moved them again.) E = 3 stays unconstructible.

**Gate honesty, measured.** None of the three equivalence harnesses drives any
predictor mutation: `equivalence.R`, `multinomial-equivalence.R` and
`bcf-equivalence.R` all return zero hits for
`setPredictor|updatePredictor|PerObservation`. `equivalence.R`'s only mutation
surfaces are `setData`, `setWeights` and `setTestOffset` (`:646-665`). The trio as
it stands is necessary and NOT sufficient for this arc, on all three baselines.
S0 fixes that before any engine change; see "The harness this arc needs".
`bench-sampler.R` has `setPredictor` accept and reject arms and NO creation arm,
so no measured sampler-creation cost exists in this repository.

## The invariant asymmetry

The forced path (`forceRefreshTrees`) re-routes every forest and the variance
forest and COLLAPSES any leaf that empties, so it maintains I - every ensemble's
live partition agrees with the store's current codes, and no bottom node is
unoccupied - by changing tree structure. The transactional path maintains I for
forest 0 only, and by VETO. On ACCEPT, forests 1.. and the variance forest are
never repartitioned, so their `Node::begin/end` ranges still describe the
PRE-change partition while the store holds POST-change codes. That is the "routed
against stale codes" the guard comment names, and it is what the probe below
measured.

Two consequences the design must respect:

1. Rollback is already three-quarters correct: `repartitionTrees` loops all mean
   forests. It does not reach `varianceForest_`, which is safe only because
   nothing repartitions it during a transaction. The moment S3 repartitions
   variance trees in the validate phase, `repartitionTrees` MUST gain the variance
   arm or a rejected transaction leaves `s^2(x)` routed by the proposal.
2. **The cell guard and the revalidation must be widened TOGETHER.** Widening
   `revalidateTrees` alone lets an ensemble outside the session's cache empty a
   leaf; `finalize()` then returns false and
   `bartcore_updatePredictorPerObservation` responds with `Rf_error("...produced a
   tree with an empty leaf")` AFTER `commitObservation` has written cells with no
   journal. That converts a recoverable wrong answer into an unrecoverable state.
   F7 exists to falsify it, with a mandatory negative half.

## The live defect, and the stop-loss

`dbarts_sampler_setPredictor` (`src/C_interface.cpp:239-254`) and
`dbarts_sampler_updatePredictor` (`:256-275`) carry NO multi-forest and no
variance-forest guard. Their siblings `setResponse` (`:199`), `setOffset` (`:208`),
`setWeights` (`:215`) and `setTestOffset` (`:295`) all carry
`refuseMultiForestMutation`. This is reachable TODAY: `dbarts_sampler_create`
routes to `createHolder`, which calls `applyVarianceAttributes`, so a
heteroscedastic sampler is flat-creatable from a `dbartsSpec(variance = ...)`
triple.

CONFIRMED by probe, built through the supported `LinkingTo` path against an
installed archive of the tree: both flat entries return 1 (accepted) at
`forceUpdate = 0` on a heteroscedastic sampler while all three R-bridge
transactional entries refuse. Differentially, against the same sampler taking the
supported forced path on a row-reversal replacement (which preserves every
column's values and cut grid, so no leaf can empty and the forced collapse is a
certified no-op):

    HETEROSCEDASTIC   sigma    max|A-B| = 0
                      train    max|A-B| = 2.589   (sd 0.553)
                      variance max|A-B| = 6.222   (sd 0.986)
    HOMOSCEDASTIC     sigma 8.88e-16, train 4.49e-14 (reassembly order only)

The per-row combined variance is wrong by 6.3 times its own standard deviation.
No collapse occurred (variance flat-tree size 162 before and after on both arms),
the variance leaf VALUES are bitwise identical across the swap, and everything
stays finite: the entire difference is routing. It is a silent wrong answer, not
a crash - each variance tree owns its index buffer, so the stale partition is
self-consistent and merely wrong.

Scope: heteroscedastic ONLY at fed324e. `resolveFamily` has no multinomial branch
and `createHolder` never reaches `createBCFHolder`, so neither multi-forest shape
is flat-creatable yet. **bcf-public-surface S1 changes that**
(`bcf-public-surface.md:296-318`: `createHolder` reads the treatment slot and
routes to `createBCFSampler`), at which point the same hole strands tau on every
flat BCF.

R5 wrinkle worth recording, because a reader checking this by hand gets a false
negative: `sampler$setPredictor(x)` with `forceUpdate` MISSING defaults to TRUE
for a whole matrix (`R/bartcore.R:118-122`, `forceUpdate <- is.null(column)`), so
the transactional path must be requested explicitly.

**Disposition: a standalone guard commit (SL), landing out of band, BEFORE
bcf-public-surface S1.** It does not wait for this arc and it does not belong to
bcf-public-surface, whose S3 is where the flat surface is otherwise touched. It
refuses, at `dbarts_sampler_setPredictor` and `dbarts_sampler_updatePredictor`,
exactly what the R bridge refuses at `bartcore_setPredictor` and
`bartcore_updatePredictor`: with `forceUpdate == 0`, a sampler with a variance
forest, or a sampler with `numForests >= 2`. Forced flavors stay open. No
`dbarts.h` change, no hash re-bake - the entries already exist and only their
bodies gain a guard.

## Fork 2 - per-forest or per-sampler veto

The measurement REPORTED this and did not decide it. Decided here: **per-sampler
conjunction.** A row installs iff it empties no leaf in any tree of any ensemble
of any chain. The per-forest column mask is the opt-out and it is EXACT: a tree
with no split on the mutated column cannot veto, structurally - in the session the
override never fires so `newLeaf == oldLeaf`; in the revalidation the repartition
reproduces the same partition - and M1 measured it at exactly 0 over 14.4M
decisions.

The rejected alternative is a caller-declared per-forest veto opt-out. Stated
correctly (the memo's version described NOT REPARTITIONING the excluded forest,
which no coherent opt-out requires): exclude forest f from the VETO while still
repartitioning it. The shared-store objection then evaporates, and what makes the
state inadmissible is the empty-leaf invariant - an unoccupied bottom carries no
drawn parameter, and the engine treats that as invalid in `forceRefreshTrees`'
collapse, in `stateIsValid`'s mean branch, and in `installVarianceForest`. The
only repair is to collapse, and repartition-with-collapse for one forest IS E2
restricted to one forest, rejected below. The falsifier's shared-code-store
argument is retained where it was registered: against forest-LOCAL ROLLBACK, for
which it is genuinely load-bearing.

## Alternatives priced (YELLOW obligation 2)

**E1, the widened veto (taken).** Roll the whole transaction back, or decline the
row. Extends the shipped semantics unchanged. Price, DERIVED from the in-tree
bench arms and the tree counts - arithmetic, NOT measurement; no timing metric
exists in the falsifier and no benchmark was run for this plan:

| shape | trees touched on accept | derived cost, n = 1000 | in sweeps |
|---|---|---|---|
| single forest 75 (today) | 75 | 0.247 ms (MEASURED) | 1.4 |
| BCF 75/50, unpruned | 125 | ~0.41 ms | 2.4 |
| BCF 75/50, pruned (T_j(tau) ~ 1.2) | ~76 | ~0.25 ms | 1.4 |
| multinomial 4 x 75, unpruned | 300 | ~1.0 ms | 5.7 |
| multinomial 4 x 75, pruned | ~105 | ~0.35 ms | 2.0 |

Read as orders of magnitude, not measurements. The scaling assumes the accept path
is linear in trees x n, which is what the revalidate/rebuild pair is, and it
ignores the O(p n log n) cut work a transaction does not redo for an unchanged
grid. Measured inputs, all reproducing from
`benchmarks/baselines/bench-sampler-ab1dc52.csv` (n = 1000, p = 10, 75 trees, 1
chain): `setPredictor-accept` 0.247 ms/update, `setPredictor-reject` 0.136
ms/update, `run` 0.174 ms/iteration, `embedded-offset-run1` 0.248 ms per Gibbs
step.

**E2, collapse-on-object. REJECTED.** Forced semantics for the OBJECTING ensemble
only: when tau or the variance forest would empty a leaf, collapse that leaf and
install the row anyway. It is a THIRD posterior, neither the veto's nor the forced
path's. Wall-clock price comparable to E1 (a collapse is O(nodes)). Q1 (an arm-B
surrogate, REPORTED and never gated) puts per-row latent posterior means at 0.041
veto / 0.045 collapse / 0.029 frozen, sds 0.277 / 0.268 / 0.488, pooled KS 0.025
veto-vs-collapse against 0.039 veto-vs-frozen - i.e. on the surrogate the two
kernels this arc treats as categorically different are the two CLOSEST of the
three. Reject on the BENEFIT and the EVIDENCE, not on the description: the
benefit is at most 0.037 pp of install rate, the only evidence is one surrogate
that gates nothing, and there is no reversibility argument for the resulting
kernel while the one recorded consumer treats the install mask as an MH accept
mask (`bairrtt/R/irt_causal_bart.R:614`, rejected rows reverted as rejected MH
moves). A collapse is not an involution: it accepts the move AND changes tree
structure irreversibly, in scan order, and neither mandated BCF oracle
(`bcf-exact.R`, `bcf-exact-restricted.R`) can see it. The strongest form - E2
restricted to the VARIANCE forest alone, which would buy back most of the measured
cost and leave every mean-forest decision identical to E1 - is answered in open
decision 1, not here.

**E3, decline the widening for the objecting ensemble class.** Keep
`refuseVarianceForestPredictorMutation` on the transactional entries. Price: the
heteroscedastic latent-covariate class stays unbuildable; the tau half is
unaffected (tau is the sole objector in 0-3% of widened per-observation
rejections, with "The attribution caveat" attached). This is today's state, it is
a coherent shipping state, and it is S3's pre-registered fallback.

**E4, rebuild the sampler instead.** What the refusal message currently
recommends. For the MEAN forests it is the same acceptance criterion at higher
cost with the rollback removed: `stateIsValid` imposes occupancy on every tree of
every mean forest as an all-or-none hard refusal (no partial install, no rollback
to the previously accepted state, an R error rather than a FALSE), and the
rebuild's work is a strict superset of E1's - a full store rebuild including a
per-column cut grid, the flat-tree rebuild, plus the same repartitions. No
creation-cost benchmark exists here and none was run; the dominance argument for
the mean forests needs no constant. For the VARIANCE forest the dominance claim is
FALSE (adjudication 4): E4 is a weaker, unchecked route until S3 item 5 lands, and
the refusal message may not retire before that.

## Pruning (YELLOW obligation 1)

Two distinct benefits, and they are different.

**(a) The widening is exactly free for a mask-restricted consumer.** A tau forest
restricted to moderators, or a `variance = ~subset` forest, contributes zero
j-splitting trees for a latent outside its mask, hence zero vetoes. No new code:
the trees already carry `columnMask` and the argument is structural. M1 = 0 over
14.4M decisions is its receipt.

**(b) It removes most of the widening's time and memory cost.** Measured `T_j`
(j-splitting trees, summed over ensembles and chains) is ~10 of 75 for a single
forest at p = 10, ~1.1-1.4 of 50 for tau under `treatment.base = 0.25 /
treatment.power = 3`, and 5-6 of 40 for the variance forest. So a pruned session
caches ~11 trees where an unpruned widened session caches 125. **This is a DERIVED
arithmetic claim. The falsifier has no timing metric anywhere and licenses nothing
about speed.** It may not appear in NEWS or in a docstring without a
`bench-sampler` arm behind it.

Where pruning applies:

- **The session cache: prune, mask-exact.** Dropping non-j-splitting trees from
  `observationLeaf_` and `leafCounts_` cannot change the returned mask: for those
  trees the column override never fires during the descent, so `newLeaf ==
  oldLeaf` always, `valid` is untouched, and no move is staged. The permutation is
  drawn before the loop and the early-exit conjunction reaches the same boolean.
  Bitwise neutral on the single-forest path, INCLUDING forest 0's own
  non-j-splitting trees; F5 and the S0 `predpartial` scenario gate it.
- **The revalidate/rebuild loop for forest 0: DO NOT prune.**
  `rebuildFitsFromParameters` subtracts, rewrites and adds back per tree, and
  `(T - f) + f != T` on 3.1% of rows at the loop's magnitudes. Pruning forest 0
  would force an equivalence re-record. Forest 0's loop stays byte-for-byte
  today's; that is what keeps the trio bitwise.
- **The revalidate/rebuild loop for f >= 1 and the variance forest: prune, and the
  pruned form is the EXACT one.** Those paths do not exist today, so there is no
  baseline to preserve. Skipping a tree with no split on a touched column is
  exactly correct - its partition, its recovered parameters and its fits are
  unchanged - and its contribution to `forest.totalFits` is preserved bitwise
  rather than round tripped through a subtract and an add. Pruning here is not
  merely an optimization; it is the more accurate implementation. Forest 0 keeps
  the ULP noise because that noise is the recorded baseline.

The predicate is `Tree::countVariableUses` over the TOUCHED column set, computed
once per transaction. For a whole-matrix `setPredictor` the touched set is every
column and pruning is vacuous; for `updatePredictor(columns)` and the
single-column session it is real; under `updateCutPoints` it is the same set (a
rebuilt grid moves only the mutated columns' codes). Say that rather than
overclaim.

Two structural requirements, both load-bearing:

1. **One predicate, one survivor list, shared by the validate and the rebuild
   phase.** `revalidateAllChains` hands `params` from phase one to phase two; the
   survivor list travels with it. An empty parameter vector is NOT a legal skip
   marker - function leaves legitimately produce one.
2. **The session cache and the revalidation use the SAME predicate over the SAME
   touched set.** That is what makes F7 (the pairing invariant) provable: the
   guarded set is a subset of the revalidated set, and every tree in the
   difference is one whose partition the revalidation reproduces unchanged from a
   valid pre-state, so it cannot fail. The subset is proper precisely on forest
   0's non-j-splitting trees, which the session prunes and the rebuild does not.

## The attribution caveat (YELLOW obligation 3, standing)

The P4/T2 table reports that the variance forest is the sole objector in 21-36% of
widened per-observation rejections and tau in 0-3%. **No sentence in this plan, in
any commit message, or in any code comment cites that fact without this caveat
named beside it.**

1. Oracle-attributed, not engine-attributed, for the shapes that matter. Every
   shipped surface returns a conjunction and never an attribution:
   `setPredictor(partial)` returns one mask over forest 0 across all chains, the
   joint call ANDs the per-sampler validities before committing, and `setPredictor`
   returns one boolean. Only arm B's per-SAMPLER marginals are engine-measured, by
   seed-pinned subset replay; per-FOREST attribution inside one sampler is
   oracle-only.
2. Confounded with the PRIOR CONFIGURATION, not a property of "variance forests".
   tau is stumpy because `treatment.base = 0.25 / treatment.power = 3` gives it
   mean leaf depth ~0.4 and ~1.2 j-splitting trees of 50; the variance forest
   carries 5-6 of 40 at depth ~1.3 because it runs under the MEAN model's tree
   prior at the default `n.trees.variance = 40`. Change `base.variance`,
   `power.variance` or `n.trees.variance` and the number moves. The correct
   statement is "an ensemble's share of the veto is its share of `T_j`", the count
   law the same measurement confirmed to within 8% under a 4x widening.
3. Conditional on one response surface (deviation D4), one cut grid, and the n/p/
   move grid the pre-registration pinned.
4. T1's ratio clause, the only non-GREEN clause, is Poisson-noise dominated at the
   registered replication (0 to 11 whole-transaction vetoes per cell); an ungated
   40-seed supplement found every interval covering 1. Cite T1 as "YELLOW as
   registered, and the supplement says the ratio failures are noise", never as
   evidence of a cost.

## The snapshot, closed structurally

S1 REPLACES `FuzzSnapshot` rather than extending it, so that "the snapshot omits a
field" stops being expressible:

    struct FuzzSnapshot {
      SamplerStateData state;   // Sampler::getState - the persisted state, entire
      LiveGeometry geom;        // the live structures getState does not carry
    };

- `state` is compared with `statesAgree` (`tests/cpp/common.cpp:25`), which
  already compares every persisted field per chain and per forest and which the
  state round-trip fuzz keeps honest.
- `geom` is built by INDEX LOOPS ONLY, with no forest literal anywhere: for each
  chain, for `f` in `numForests()`, for `t` in `numTreesInForest(f)`, capture the
  tree's split rules, its `begin`/`end` per node and its `indices`, plus
  `forestTreeFits(f)` and that forest's `totalFits`; then, under
  `hasVarianceForest()`, the same per variance tree plus
  `varianceFactorsForTesting()` and `varianceFits()`. A new forest class is covered
  on arrival; nothing needs remembering.
- Equality is `statesAgree(a.state, b.state) && a.geom == b.geom`.
- Two tripwires: a `static_assert(sizeof(ChainStateData) == N)` beside
  `statesAgree`, with a comment saying a new persisted field must gain a comparison
  here; and F4, a table-driven test that perturbs ONE captured family at a time and
  requires each perturbation to make the comparison false. The `sizeof` guard is
  honest but not airtight - a small field can hide in padding - which is why F4
  exists beside it.
- One engine addition: a per-forest `Chain::totalFitsInForest(f)` beside the
  forest-0 `totalFits()` (`chain.hpp:2800`).

Related, and a real gap neither design document names:
`fuzzInvariantViolation`'s "totalFits != tree-order sum" check
(`test_fuzz.cpp:118-123`) sums forest 0's fits against forest 0's total. It is
self-consistent, so it will not misfire on a BCF - it simply does not cover forest
1. S1 generalizes it to loop forests. The rest of `fuzzInvariantViolation` is
already forest-generic (routing agreement over `ch.numForests()` at `:126-130`, the
variance forest and the `combinedVariance == prod(factors)` identity at
`:131-151`), which is why the BCF fuzz config is scheduled in S1, beside the code
it gates, and not in S4.

## The harness this arc needs

The trio is blind to predictor mutation (measured; see "Gate honesty"). S0 fixes
that BEFORE any engine change, so the recorded baseline is a true pre-change
baseline and every later slice compares against a harness that can see the paths
it moves.

`equivalence.R`: `fitViaSamplerApi`'s `scenario$mutate` block (`:654-662`) gains
predictor keys beside `weights` and `offset.test`. New scenarios, each burning in,
mutating, then running a short second leg so carried-over tree state matters:

- `predswap` - transactional whole-matrix `setPredictor(x2, forceUpdate = FALSE)`.
- `predcol` - transactional `setPredictor(v, column = j)` on a column subset.
- `predpartial` - the per-observation session, `setPredictor(v, column = j,
  "partial")`. The only scenario that consumes the scan permutation, so it is the
  one that pins S2's session refactor.
- `predreject` - an `updatePredictor` proposal built to be rolled back (a
  two-level replacement column at `updateCutPoints = FALSE`). It gates that a
  rejected transaction leaves the run bitwise; it cannot by itself distinguish
  accept from reject, and the comment says so.
- `predforce` - the forced whole-matrix swap, which routes through
  `forceRefreshTrees`.
- `hetforce` - a heteroscedastic sampler (`samplerArgs = list(variance = TRUE,
  n.trees.variance = 40L)`, spliced through the existing `samplerArgs` mechanism)
  taking a forced whole-matrix swap. There is NO heteroscedastic scenario in the
  harness today, so `refreshVarianceForest` - which S3 splits - has zero
  equivalence coverage. This scenario is S3's abort gate.

`bcf-equivalence.R` and `multinomial-equivalence.R` each gain one forced
whole-matrix swap scenario at S0 (the only predictor mutation those shapes accept
today), recorded per forest, so S1's loop change is gated on the forced path for
both shapes.

Baseline procedure, which the harness already supports: a scenario absent from a
baseline is reported SKIPPED, not failed, and "the anchor re-records at landing"
is the established pattern (`equivalence.R:427-428, 457-458, 482-483`). So S0
records all three baselines at its own tip - `equivalence-<S0>.rds`,
`multinomial-equivalence-<S0>.rds`, `bcf-equivalence-<S0>.rds` - and the S0 commit
must first show every PRE-EXISTING scenario bitwise identical to
`equivalence-c8f661a`, `multinomial-equivalence-ec2a3d0` and
`bcf-equivalence-c820227` (or to whatever bcf-public-surface leaves current) before
the new ones are added. S1, S2 and S3 then compare against the S0 baselines. Each
of S1/S2/S3 additionally ADDS the scenarios its own new paths make legal
(BCF/multinomial transactional at S1, the BCF/multinomial per-observation session
at S2, heteroscedastic transactional and per-observation at S3) and re-records at
that tip; those are new streams with no prior baseline and become the regression
floor for later work.

## Slices

### SL. Stop-loss. Standalone, out of band, before bcf-public-surface S1.

1. Move `refuseVarianceForestPredictorMutation` and
   `refuseMultiForestTransactionalUpdate` out of the anonymous namespace in
   `R_interface_bartcore.cpp` into `bartcore_bridge`, and declare them in
   `R_interface_bartcore_common.hpp` beside `refuseMultiForestMutation` - the
   one-shared-helper mechanism bcf-public-surface S3 prescribes, so the two
   surfaces cannot diverge.
2. Call `refuseMultiForestTransactionalUpdate(engine, "dbarts_sampler_setPredictor",
   forceUpdate != 0)` at `C_interface.cpp:239`, and the same at
   `dbarts_sampler_updatePredictor` (`:256`), before the column validation.
3. `inst/tinytest/capi/consumer.c` drives a heteroscedastic sampler built through
   `dbarts_sampler_create` and asserts the refusal at `forceUpdate = 0` and the
   success at `forceUpdate = 1`, for both entries.
4. A comment at both sites recording what is live and what is defensive: the
   variance clause is reachable today; the `numForests >= 2` clause is unreachable
   until bcf-public-surface S1 makes BCF flat-creatable and is load-bearing from
   that tip on.

rng: NEUTRAL. Gates: `R CMD INSTALL` into the private library; `test-capi.R`; full
tinytest; trio bitwise (a formality that catches an accidental edit);
`air format --check .`. NO `dbarts.h` change and no hash re-bake.
ABORT: any trio divergence.

### S0. Pin the surface; give the trio eyes. No engine change.

1. Pins for the refusals this arc will retire, so a later slice cannot open one
   silently. **Ownership, settled here:** bcf-public-surface S0 already pins the
   BCF transactional `setPredictor` and the per-observation session
   (`bcf-public-surface.md:286-289`) and lands FIRST, so it OWNS the BCF pins.
   This slice writes only the pins that arc does not: on a MULTINOMIAL and on a
   HETEROSCEDASTIC sampler, transactional `setPredictor` / `updatePredictor` and
   both per-observation entries REFUSE with their current messages, while forced
   `setPredictor` / `updatePredictor` and `setCutPoints` SUCCEED. Beside the
   existing cases in `inst/tinytest/test-multi-forest-seam.R` and
   `test-heteroscedastic-mutation.R`. The BCF pins are INVERTED, not duplicated
   and not deleted wholesale, by S1 and S2 - each retirement commit edits the
   assertion bcf-public-surface wrote, in place, from REFUSE to ACCEPT, so no arc
   deletes another's tests behind its back and no tinytest goes red across an arc
   boundary.
2. The harness extension and the three recorded baselines, per "The harness this
   arc needs".

rng: NEUTRAL. Gates: full tinytest; every pre-existing equivalence scenario bitwise
identical to the current baselines BEFORE the new scenarios are added;
`air format --check .`; `lintr::lint` on touched R files.
ABORT: any pre-existing scenario moving.

### S1. Widen the two-phase revalidation across `forests_`.

Lifts the whole-matrix and subset transaction for BCF and multinomial.

1. `Chain::revalidateTrees` and `Chain::rebuildFitsFromParameters` loop every
   forest. Introduce `using ForestParameters = std::vector<TreeParameters>` with a
   per-forest survivor list and give the transactional pair the new type. LEAVE
   `recoverTreeParameters` and `applyNewData` on their forest-0 `TreeParameters`
   signatures: `setData` is still refused at `numForests >= 2` by
   `refuseMultiForestMutation` and is not in scope.
2. Forest 0's inner loop is byte-for-byte today's. Forests f >= 1 are j-split
   pruned per "Pruning", with the survivor list shared between the two phases.
3. `Chain::repartitionTrees` needs no change (it already loops `forests_`). Assert
   that by test rather than by reading.
4. `FuzzSnapshot` rebuilt per "The snapshot, closed structurally", plus the
   per-forest `totalFitsInForest` accessor and the generalized totalFits invariant
   check.
5. A BCF `ConfigSpec` and a `fuzzRunBCF` sibling of `fuzzRunConstant`
   (`test_fuzz.cpp:506-546`) constructing `Sampler<ConstantGaussianLeaf>` with a
   `BCFSpec` and a treatment vector; op mask covering `OP_SET_PREDICTOR`,
   `OP_UPDATE_COLUMNS`, `OP_SET_CUTS`, `OP_RUN`, `OP_STATE` and EXCLUDING
   `OP_SET_DATA`, `OP_SET_RESPONSE`, `OP_SET_WEIGHTS`, `OP_SET_OFFSET` (refused
   for a multi-forest sampler). The per-observation ops join at S2.
6. Bridge: `refuseMultiForestTransactionalUpdate` drops its `numForests >= 2`
   clause at `bartcore_setPredictor` and `bartcore_updatePredictor` only. The
   per-observation entries keep it until S2; the variance clause stays until S3.
   The message text of the entries still refusing stops naming "make a new sampler
   instead" as the remedy and names the forced swap instead: E4 is priced out for
   the mean forests, and for the variance forest it is UNCHECKED until S3 item 5.
7. The BCF and multinomial transactional equivalence scenarios, recorded at this
   tip.

rng: NEUTRAL on every existing path. Gates: `R CMD INSTALL --preclean` into a
PRIVATE library at `.claude/multiforest-predictor-mutation-design/privlib`
(chain.hpp is a header); delete the `benchmarks/kernels` binaries; `tests/cpp` from
clean, plain AND ASAN/UBSAN; full tinytest with NO snapshot regenerated; the trio
bitwise on the S0 baselines; the mandated BCF oracle set (`bcf-exact.R` quick,
`bcf-exact-restricted.R`) UNMOVED; `air format --check .`.
ABORT: any trio divergence.

### S2. Widen `UpdateSessionImpl`.

Lifts the per-observation session for BCF and multinomial.

1. The session's flat index becomes an explicit survivor table of
   `(chain, forest, tree)` triples built in the constructor - NOT offset
   arithmetic, because the cached set is pruned and therefore sparse. `treeAt(t)`
   resolves through a new `Chain::treeInForest(f, t)` const accessor, promoted
   from `treeInForestForTesting`. `Chain::tree` keeps its forest-0 meaning for its
   other callers.
2. j-split pruning of the cache, on the SAME predicate the revalidation uses.
3. Bridge: drop the `numForests >= 2` clause at
   `bartcore_updatePredictorPerObservation` and at the joint entry.
4. The joint entry's cross-sampler contract is unchanged: one permutation from
   `samplers[0]->rng()`, install in every sampler or none, AND the `finalize()`
   returns.
5. `OP_PER_OBS` and `OP_SESSION_ABANDON` join the BCF fuzz config.
6. The BCF and multinomial per-observation equivalence scenarios, recorded at this
   tip.

rng: NEUTRAL on every existing path, INCLUDING the single-forest per-observation
session, whose returned mask the cache pruning cannot change (see "Pruning"); the
`predpartial` scenario is the gate on that claim. Gates: as S1, plus the widened
fuzz config and F5's oracle run.
ABORT: any trio divergence.

### S3. The variance forest.

1. Split `refreshVarianceForest` into `bool
   revalidateVarianceTrees(TreeParameters&)` - recover node-indexed factors through
   the LIVE partition FIRST (the ordering `refreshVarianceForest` documents),
   repartition, report `bottomNodesAreOccupied()`, collapsing nothing, dropping no
   directions and scattering nothing - and `void
   rebuildVarianceFactors(const TreeParameters&)` - drop stale directions, scatter
   through the new partition, recompute `combinedVariance`. The existing
   collapse-and-remap form STAYS for `forceRefreshTrees` and `applyNewData`,
   unchanged and byte-for-byte; `hetforce` is the gate on that.
2. `Chain::revalidateTrees` and `rebuildFitsFromParameters` call the two halves
   under `varianceForest_ != nullptr`, APPENDED after the `forests_` body - the
   variance arc's append-only discipline, so a homoscedastic chain's code path is
   unchanged. The neutrality claim is that the runtime branch is not taken, NOT
   that it is compiled out: `varianceForest_` is a `unique_ptr` and every guard on
   it is a runtime test.
3. `Chain::repartitionTrees` gains the variance arm. Load-bearing from this slice
   on. `combinedVariance` and `factorByTree` are untouched by the validate phase,
   so restoring the partition restores the state exactly.
4. The session cache admits the variance trees, through a promoted
   `Chain::varianceTree(j)`; the survivor table gains a variance arm (there is no
   `numTreesInForest` entry for it - use `numVarianceTrees()`).
5. **Close `stateIsValid`'s variance occupancy gap** (adjudication 4): the variance
   branch (`chain.hpp:2369-2380`) gains the scratch build, the repartition and the
   `bottomNodesAreOccupied()` test its mean sibling has at `:2324`. Without this
   the arc retires a refusal message that points users at a route WEAKER than the
   veto it is giving them. This CHANGES `setState` behavior - a variance state with
   an unoccupied bottom stops installing - so it carries its own tinytest and the
   NEWS bullet says so. If VD declines the variance forest (open decision 1), this
   item is re-homed to `docs/plans/variance-forest-mutation-routing.md` as
   unblocked work rather than dropped.
6. Bridge: `refuseVarianceForestPredictorMutation` retires from the four
   transactional entries. It keeps no other call site; DELETE it deliberately
   rather than leave it dead by accident.
7. The heteroscedastic config in `test_fuzz.cpp` re-admits `OP_SET_PREDICTOR |
   OP_UPDATE_COLUMNS | OP_PER_OBS | OP_SESSION_ABANDON`; the "for the whole arc"
   comment at `:643-648` is EDITED, not silently flipped.
8. The heteroscedastic transactional and per-observation equivalence scenarios,
   recorded at this tip.

Interaction with zero-weight-exactness, confirmed: occupancy is COUNT-based (that
arc's binding decision 2), so a per-forest weight - including an exact zero -
changes no veto and no leaf's occupancy. The two mechanisms are orthogonal by
construction.

rng: NEUTRAL for every homoscedastic path and for the heteroscedastic forced and
setData paths. Gates: as S1, with the ASAN/UBSAN leg MANDATORY (new engine code on
a newly reachable path), plus `test-heteroscedastic.R` and
`test-heteroscedastic-mutation.R` passing verbatim, plus `hetforce` bitwise.
FALLBACK, pre-registered: if the validate/rebuild split cannot preserve the
geometric-merge contract or the rollback identity, take E3 - keep the variance
refusal, land SL and S0-S2, and record the door as declined with the reason.
ABORT: any trio divergence, `hetforce` included.

### S4. Docs, records, residual coverage.

1. R docs for `setPredictor`, `updatePredictorPerObservation` and
   `updatePredictorPerObservationJointly`; `inst/NEWS.Rd`.
2. `docs/design/model-space-survey.md` secs 4-5; `docs/design/heteroscedastic.md`;
   `docs/design/empty-leaf-veto.md` (record that the criterion the transaction
   enforces is the one `stateIsValid` already enforced, and that S3 closed the
   variance asymmetry); bump each file's `Status:` line per the plans README
   landing rule.
3. `docs/plans/c-api-growth.md`: record that this arc added NO X-list entry and
   moved NO hash, so the queued reshape is unaffected.
4. TODO edits (below) and the Landing note in this file.

rng: NEUTRAL. Gates: `tests/cpp` from clean plus ASAN; full tinytest; trio bitwise;
`air format --check .`; `lintr::lint` on touched R files.

## Falsifiers (pre-registered)

- **F1 (every slice).** The trio is bitwise on all three S0 baselines. Any
  divergence on a pre-existing scenario is a LEAK and aborts the slice. From S0 on
  the trio DOES exercise predictor mutation, so it is a real gate here and not
  only a necessary one.
- **F2 (S1), the load-bearing one, BOTH halves.** On a BCF and on a multinomial
  sampler, after an ACCEPTED transactional `setPredictor`, every observation of
  every tree of every forest of every chain sits in the leaf
  `findBottomNodeForObservation` routes it to. Reuse `fuzzInvariantViolation`,
  which already loops `ch.numForests()`; do not write a second one. NEGATIVE HALF,
  mandatory: with the `f >= 1` arm of the widened revalidation removed, it must
  FAIL. A green suite that cannot go red is not evidence.
- **F3 (S1).** A REJECTED transactional update on a BCF leaves the sampler
  unchanged under the WIDENED snapshot comparison, not the old one.
- **F4 (S1), the structural gate on the snapshot.** Table-driven: perturb one
  captured family at a time - tau's fits, a tau partition range, a mean split rule,
  a variance factor, `combinedVariance`, sigma, a code, a cut - and require each to
  make `fuzzSnapshotsEqual` false. One row per family; a family with no row is a
  family the snapshot does not cover.
- **F5 (S2), mask exactness. Replaces the memo's un-constructible R oracle.**
  Driven from `tests/cpp` through `beginPredictorUpdate` at a CALLER-CHOSEN scan
  order, so nothing depends on the engine's drawn permutation: the engine's install
  decisions must equal, BITWISE, an in-test oracle that maintains its own per-tree
  leaf counts over every tree of every forest (and, from S3, every variance tree),
  over >= 1e5 decisions on BCF, multinomial and heteroscedastic configurations.
  NEGATIVE HALF: corrupt the pruning predicate to drop a j-splitting tree and it
  must fail.
- **F6 (S2), untouched-tree exactness. Replaces the memo's pruning falsifier.** On
  the shipped (pruned) build, after any accepted transaction, every tree with no
  split on a touched column has bitwise-unchanged rules, partition ranges,
  `indices`, recovered parameters and fit slab, and its forest's `totalFits`
  retains its contribution bitwise. Asserted for f >= 1 and the variance forest.
  NOT asserted for forest 0, which round trips through subtract-and-add by rule and
  differs at the ULP on ~3% of rows - state that in the test comment so no later
  reader reads the asymmetry as a bug. NEGATIVE HALF: skip a tree that DOES split
  on a touched column and F2 must go red.
- **F7 (S2), the pairing invariant.** `finalize()` returns true by construction:
  the cell guard and the revalidation quantify over sets related by the subset
  argument in "Pruning", so `bartcore_updatePredictorPerObservation`'s "produced a
  tree with an empty leaf" error is unreachable. Assert as a `tests/cpp` invariant
  across the fuzz surface. NEGATIVE HALF, mandatory: widen the revalidation WITHOUT
  widening the cell guard, show the hard error firing after cells have been written
  - the unrecoverable state - and revert.
- **F8 (S2).** `OP_SESSION_ABANDON` on the BCF config: a session dropped without
  commit or finalize leaves the sampler unchanged under the WIDENED snapshot.
- **F9 (S3), rollback.** Under a variance forest, a rolled-back transaction
  restores the variance trees' partitions and `combinedVariance` bitwise. NEGATIVE
  HALF: omit `repartitionTrees`' variance arm and it must fail.
- **F10 (S3), ordering.** The factors recovered in the validate phase equal the
  pre-transaction `varianceFactorsForTesting()` slab exactly. NEGATIVE HALF: move
  the recovery after the repartition and it must fail (it then reads the new
  partition's members out of the old leaves' slots). It need NOT cover the
  missing-direction drop, which is routing-neutral by construction
  (`tree.hpp:1062-1067`, `:1213-1223`).
- **F11 (S3), forced-path invariance.** `hetforce` and the heteroscedastic
  `setData` path are bitwise across the `refreshVarianceForest` split. Divergence
  aborts S3; it is never a re-record.
- **F12 (standing, not new work).** `tests/cpp/test_model.cpp:718-733` already
  asserts that `createBCFSampler` and `createMultinomialSampler` return nullptr on
  `numVarianceTrees > 0`. E = 3 stays unconstructible; no slice may make it
  reachable as a side effect.
- **F13 (SL), the live hole.** A heteroscedastic sampler created through
  `dbarts_sampler_create` refuses `dbarts_sampler_setPredictor(..., 0, 0)` and
  accepts `(..., 1, 0)`. Pre-fix, the refusal must be shown ABSENT in the same
  file, so the commit shows the inversion.

## Open decisions (VD)

### 1. The variance forest: in this arc, or permanently refused?

The class is the heteroscedastic sampler - one mean forest plus a variance forest
of scale trees, `numForests == 1` but `hasVarianceForest` true - taking
transactional and per-observation predictor mutation. It is S3.

The evidence changed since the memo. `Chain::stateIsValid` does NOT impose
occupancy on variance trees, and `setState` installs them through
`rebuildVarianceForest` with no occupancy test, so the "just rebuild the sampler"
alternative is not the same criterion at higher cost for this class: it is a WEAKER
route that admits variance states the veto refuses, with only a release-disabled
assertion behind it. This class also carries essentially the whole measured
whole-transaction cost of the widening (sole objector in 21-36% of widened
per-observation rejections, with "The attribution caveat" attached; tau 0-3%), and
it is the one whose own arc parked this door here explicitly.

- **(a) In, as S3, and close the `stateIsValid` gap in the same slice
  (RECOMMENDED).** Cost ~425 lines, the hardest mechanism in the arc, and one new
  refusal on `setState`. Buys the heteroscedastic latent-covariate class the
  variance arc named - IRT response-time or careless-responding variance keyed on
  latent ability, latent-confounder sensitivity where residual scale moves with the
  confounder - and makes the arc's retirement of the refusal message honest.
- **(b) In, gap left open.** Saves ~30 lines and the new refusal. Leaves a route by
  which a consumer refused by the veto still installs an unchecked variance state.
  Not recommended: the arc would add a criterion on one path while knowingly
  leaving the sibling path weaker.
- **(c) Decline permanently (E3).** Keep `refuseVarianceForestPredictorMutation` on
  the transactional entries; land SL and S0-S2. Coherent, and it is today's state.
  Costs the enabling class. Under VD's enabling-value gate, absence of a consumer
  today is not the gating fact, and this class is named. The `stateIsValid` gap
  should then still be closed, in the variance arc.

A fourth option exists and is priced in "Alternatives priced": E2 restricted to the
variance forest alone (collapse the objecting variance leaf, install the row,
leave every mean-forest decision identical to E1). It buys back most of the
measured cost and leaves the semantics bairrtt's MH filter depends on untouched.
It is rejected because it is a third posterior with no reversibility argument, no
oracle that can see it, and a benefit measured at 0.037 pp; it is recorded so it is
not re-proposed as a fresh idea.

Recommendation: (a). S3 already carries the pre-registered fallback to (c) if the
validate/rebuild split cannot hold the geometric-merge contract.

RESOLVED: (a), 2026-08-10 - VD granted "proceed at your discretion" over the
open recommendations; the variance forest is IN as S3 with the `stateIsValid`
occupancy gap closed in the same slice, the fallback to (c) staying
pre-registered.

### 2. Sequencing: what actually remains open

Nothing about the ORDER. VD ruled 2026-08-09 that bcf-public-surface goes first,
and this plan sequences inside that ruling. Pin ownership is settled in S0 item 1:
bcf-public-surface S0 owns the BCF pins, and this arc's S1/S2/S3 invert exactly
those assertions in the same commit that retires each refusal.

The one thing that needs a yes: **may SL land out of band, ahead of
bcf-public-surface S1?** It is ~155 lines, touches two flat entries and one header,
changes no `dbarts.h` symbol and moves no hash.

- **(a) Land SL now, as its own commit (RECOMMENDED).** Closes a live silent 6-sd
  corruption reachable today by any `LinkingTo` consumer that builds a
  heteroscedastic sampler, and forecloses the strictly larger BCF exposure that
  bcf-public-surface S1 opens.
- **(b) Fold the guard into bcf-public-surface S1.** One fewer commit; leaves the
  heteroscedastic hole open for the length of that arc and puts an unrelated guard
  inside its architectural slice.
- **(c) Wait for this arc's S0.** Leaves both holes open across the whole of
  bcf-public-surface. Not recommended.

### 3. Multinomial: named and gated, named and ungated, or excluded?

The class is a K-forest multinomial (softmax) sampler taking the same widened
mutation. It rides the same `forests_` loop as BCF; the existing predicate is
`numForests >= 2`, which multinomial trips, so EXCLUDING it costs code rather than
saving it.

Corrected facts. The flat-C hole is heteroscedastic-only today and
bcf-public-surface does not make multinomial flat-creatable, so no C consumer
reaches this shape. `FuzzSnapshot` is blind to every forest but forest 0 until S1
rebuilds it, so today's fuzz coverage of a multinomial rollback would be vacuous -
which is an argument for the S1 snapshot work, not against the shape. The measured
margin is a factor of 34 against the nearest GREEN line, not three orders (memo
error, corrected). Measured UNGATED at K = 4: per-observation reject 0.010%
continuous / 0.087% binary, against 0.002% / 0.028% at `E_0`. Highest price of any
shape in scope, and the shape where `T_j` is largest (~40).

- **(a) Named and ungated (RECOMMENDED).** Include it, and say in this plan, in
  NEWS and in the docs that multinomial's veto rate was measured as ungated stress
  rather than covered by a verdict. Note what IS gated: correctness - F2, F5, F6
  and F8 all run on a multinomial configuration. What is ungated is only the
  modelling cost of the rate.
- **(b) Named and gated.** Re-run the falsifier's arm at K = 4 with pre-registered
  thresholds before S1. The harness lives under `.claude/` and is not a repository
  artifact, so this costs a rebuild plus a measurement cycle, to gate a shape with
  no consumer today.
- **(c) Excluded.** Requires a NEW predicate distinguishing multinomial from BCF,
  and leaves a refusal whose stated reason has been removed.

Recommendation: (a).

RESOLVED: (a), 2026-08-10, same discretion grant - multinomial is named and
ungated, with the ungated status stated in NEWS and the docs as specified.
Question 2's SL interleave was likewise taken as recommended and SL has
LANDED (see Landing notes).

## Migration costs per consumer (enumerated, never constraining)

- **bairrtt** (`/Users/vdorie/Repositories/bairrtt`, live R5 driver: two
  `dbarts::dbarts` samplers at `irt_causal_bart.R:446, 458`,
  `updatePredictorPerObservationJointly` at `:609`,
  `setPredictor(forceUpdate = TRUE)` during burn at `:619, 624`, rejected rows
  reverted as rejected MH moves at `:614`). Existing code: ZERO change - it is
  homoscedastic single-forest and every path it uses is untouched and gated bitwise
  from S0 on. What it GAINS: the ability to make its outcome model a two-forest
  causal forest with the latent ability as a treatment moderator, the motivating
  class. What it still needs BEYOND this arc: the BCF test surface - its MH filter
  calls `$predict` (`:568, 571, 636, 637, 709, 713`), which `refuseBCFTestSurface`
  refuses (`R_interface_bartcore.cpp:4634`). This arc removes ONE of its two
  blockers; stated plainly so no reader over-reads "the named consumer is
  unblocked". Its measured price for the widening at its own shape: 0.043%
  (continuous) / 0.226% (binary) of moves rejected at n = 300 with 150 trees,
  against its own prose "under 1%".
- **stan4bart** (23 flat C entries, mutation every sweep). ZERO, plus SL's guard,
  which fires only on a sampler it chose to create with a variance forest. No
  `setPredictor`/`updatePredictor` call found in `src/bart_util.cpp`,
  `parametric_sampler.hpp` or `R/stan4bart_fit.R`; the scan was partial and is
  labelled partial.
- **treatSens** (10 flat C entries, mutation every sweep). Same posture; its
  recorded mutation is on the RESPONSE side.
- **bartCause.** ZERO - pure R, no `LinkingTo`, and no sampler-mutation call
  anywhere.

## Doors held open (recorded, not scheduled)

- **E = 3, heteroscedastic BCF.** Not lifted; no verdict authorizes it (the
  falsifier's E = 3 cells are an arm-B surrogate, reported and ungated). F12 keeps
  it unconstructible.
- **`setData` on a multi-forest sampler.** `refuseMultiForestMutation` untouched.
  It is a memory-safety guard as well as a staleness guard (the combiner indexes a
  borrowed z over the live observation count); its survey is recorded in
  `docs/plans/runsbcbcf-repair.md`.
- **A caller-declared per-forest veto opt-out.** Closed: it is E2 restricted to one
  forest. Reopen only on a reversibility argument for the resulting kernel.
- **A per-forest saved-tree replay / test-treatment surface** - bairrtt's second
  blocker.
- **Pooled categorical latents (> 64 levels).** Out of scope. Related and
  pre-existing, logged not fixed: a heteroscedastic sampler whose store holds a
  pooled categorical column null-dereferences in `Tree::flattenBelow` on
  `storeState`. Because the widened `FuzzSnapshot` calls `getState` on every
  captured op, the fuzz configs must not pair a pooled categorical with a variance
  forest until that is fixed; today's heteroscedastic config uses three ordinal
  columns and does not.
- **Occupancy-aware proposals** (`docs/design/empty-leaf-veto.md:115`, a 250-400
  line posterior-changing rewrite). Not needed - that was the KILL branch and no
  KILL clause fired.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- SL: the flat C `dbarts_sampler_setPredictor` and `dbarts_sampler_updatePredictor`
  now refuse a transactional update on a heteroscedastic or multi-forest sampler,
  matching the R interface; the forced flavors are unchanged.
- S1: `setPredictor` and `updatePredictor` accept transactional (rollback) updates
  on BCF and multinomial samplers; a row is installed only if no leaf of any tree
  of any forest of any chain would empty.
- S2: the per-observation update session, including the joint variant, accepts BCF
  and multinomial samplers.
- S3: heteroscedastic samplers accept both, and `setState` now validates variance
  tree occupancy (a state with an unoccupied variance leaf is refused where it
  previously installed).
- S4, beside the S1 bullet: multinomial's veto rate was measured as ungated stress
  rather than covered by a verdict (per open decision 3, if (a) is taken).

## TODO edits at landing

- `multiforest-predictor-mutation`: replace the entry with the landing record.
  Correct its two stale anchors - `R_interface_bartcore.cpp:1909-1923` ->
  `refuseMultiForestTransactionalUpdate` (`:2136-2144` at fed324e); `chain.hpp:1484`
  -> `Chain::revalidateTrees` (`:1571` at fed324e) - and strike "The BCF half is
  only consumable behind bcf-public-surface", which bcf-public-surface's landing
  retires. Record fork 2 as settled: the veto is per sampler.
- `variance-forest-mutation-routing`: close its transactional door, or record it as
  declined with the reason, per open decision 1; either way record that
  `stateIsValid`'s variance occupancy gap is closed (S3 item 5) or re-homed there.
- `multiforest-mutation-gaps`: record that the flat predictor entries are now
  guarded and share one helper with the R bridge.
- `bcf-public-surface`: note that its S0 BCF pins are inverted in place by this
  arc's S1 and S2, and where.

## Departures from the memo and the critique (record)

1. **The pruning falsifier is not split, it is replaced.** The critique's amendment
   kept a two-build comparison and softened its claim to partitions. This plan
   removes the second build: F6 asserts untouched-tree bitwise invariance on the
   shipped build, and F5 checks the mask against an in-test oracle. Reason: for
   f >= 1 the pruned build is the exact one, so there is no "definition of record"
   to negotiate.
2. **The mask oracle drops the R harness entirely** rather than downgrading to the
   falsifier's distributional pair. Reason: the C++ session API admits a
   caller-chosen scan order, which makes an EXACT gate available, and the R harness
   is not a repository artifact.
3. **`FuzzSnapshot` is rebuilt, not extended**, with F4 and a `sizeof` tripwire
   beside it. Reason: the omitted-field pattern had to be closed structurally, not
   per test.
4. **The `dropStaleMissingDirections` ordering advisory is OVERTURNED**, with the
   receipt at `tree.hpp:1213-1223`.
5. **Two findings neither document has:** `fuzzInvariantViolation`'s totalFits check
   covers forest 0 only and must be generalized for the BCF config; and
   `Chain::totalFits()` is forest-0 pinned and needs a per-forest sibling for the
   widened snapshot.
6. **The stop-loss is a standalone commit (SL), not an item of S0.** Reason: it
   must precede bcf-public-surface S1, and the rest of S0 must not.
7. **E4's correction is used as an argument FOR lifting the variance forest**,
   where the critique presented it as a reason to price E4 honestly and let VD
   weigh it. Both are in the plan; the recommendation follows from the fact that
   the "alternative" is unchecked, not merely more expensive.
8. **The equivalence harness gains a heteroscedastic scenario.** Neither document
   noticed that `refreshVarianceForest` - which S3 splits - has zero equivalence
   coverage today. `hetforce` is S3's abort gate.
9. **The E2 rejection is re-argued on benefit and evidence**, not on the semantic
   description alone, and the variance-only steelman is recorded in open decision 1
   so it is not re-proposed as new.

## Landing notes

SL LANDED 7299b8b (2026-08-10, out of band as sequenced). Both refusal
helpers promoted into bartcore_bridge and declared in the common header;
both flat entries guarded at forceUpdate == 0; reachability comments in
place. The capi consumer test extended in its existing shape (deviation:
two new forced-variant wrappers rather than parameterizing the existing
ones, which other tests rely on unchanged). Gates double-run
(implementer + orchestrator): install, capi 70/70, tinytest 3791/0, trio
bitwise (27/27 + 3x5 + 5x7), air 0, lintr 0. Refusal matrix as
pre-registered: both entries refuse a flat-created heteroscedastic
sampler unforced and accept it forced. No dbarts.h change, no hash
re-bake. The live silent-corruption window is CLOSED; the numForests
clause sits dormant until BCF becomes flat-creatable.

S0 LANDED 1357e7d (2026-08-11). The trio gained its predictor-mutation
eyes: the six scenarios of "The harness this arc needs" (predswap,
predcol, predpartial, predreject, predforce, hetforce) plus one forced
whole-matrix swap each in the bcf and multinomial harnesses, and the
three baselines re-recorded at the engine head 33f6fdc (R-only
landing, engine binary identical). Neutrality partition: every
pre-existing scenario reproduced its superseded baseline bitwise
first (27/27 identical draws, bcf 5/5 and multinomial 3/3 on every
channel, additions reported skipped), and each new baseline
reproduces itself (33/33 under --strict-coverage, 6x7, 4x5). Pins
landed beside the existing cases in test-multi-forest-seam.R and
test-heteroscedastic-mutation.R: multinomial and heteroscedastic
samplers refuse every transactional and per-observation entry with
their current messages while the forced entries and setCutPoints
succeed; the BCF pins stay where bcf-public-surface wrote them.
equivalence.yaml re-pinned to the new baseline (the established
re-recording pattern). Deviations, all deliberate: no S0 NEWS bullet
(the plan's NEWS section lists none and the slice is user-invisible;
precedent 994c161); an s2.test channel added so hetforce gates
refreshVarianceForest directly (S3's abort gate, F11's subject);
hetforce runs binary = TRUE because sigma is structurally pinned
under a variance forest and both sigma summaries would be degenerate
(ordinal/nbinom precedent). Findings for later slices: R `$`
partial-matches list keys (the new mutate keys are read with [[ ]]
and named prefix-free); the bcf/multinomial compare loops print no
line for a run-only scenario (uncovered additions were verified by
hand at recording); the refuse* helpers now live at
R_interface_bartcore.cpp:2605/2625 after SL's promotion. Gates
double-run (implementer + independent verifier): install, tinytest
3982/0, the neutrality and self-reproduction partitions above, air 0,
lintr no new findings.

S1 LANDED 938eb81 (engine) + 7a9c6f3 (baselines), 2026-08-11/12. The
two-phase revalidation transaction runs across forests_: forest 0's
inner loop byte-for-byte (33/33 + 6 + 4 pre-existing scenarios
bitwise vs the S0 baselines), forests f >= 1 j-split pruned with the
survivor list shared between phases; ForestParameters =
std::vector<TreeParameters> on the transactional pair only
(recoverTreeParameters/applyNewData keep forest-0 signatures; setData
stays refused). Bridge: the numForests clause dropped from the
transactional entries on BOTH surfaces - deviation from the plan's
R-only wording, deliberate: SL unified the helper so the surfaces
cannot diverge, and BCF has been flat-creatable since
bcf-public-surface S1, so the flat clause was live, not dormant; the
per-observation clause moved to a new
refuseMultiForestPerObservationUpdate (S2 retires it); variance
clause untouched (S3). Still-refusing messages now name the forced
swap, not "make a new sampler"; an R-level refuseBCFMutation on the
transactional path (outside the plan's seam list) removed. Pins
inverted IN PLACE in three files - test-bcf-mutation-pins.R (a BCF
pin file the plan did not name), test-bcf-r5-surface.R,
test-multi-forest-seam.R (BCF + multinomial transactional) - with
rollback/finite-run assertions beside each; per-observation and
heteroscedastic pins still refuse. Fuzz: BCF + multinomial
ConfigSpecs (four configs), snapshot rebuilt with THREE members
(statesAgree compares nothing above the per-chain level, so store
codes + SamplerStateData cutPoints/currentSampleNum are compared
explicitly - carry to S3); snapshot captures LEAF ASSIGNMENT, not raw
index order, because a rejected transaction permutes members within a
leaf (measured, pre-existing on single-forest configs too, rollback
re-routes from the root and the partition kernel is not order-stable;
draw-invisible, pinned by predreject at the draw level, recorded in
the capture comment). Falsifiers: F2 both halves green AND verified
able to go red (forest-0-only validation reddens exactly the
multi-forest configs); F3 green; F4 shipped as
testSnapshotCoversEveryFamily (26 families, drop-test verified);
repartitionTrees no-change asserted by test (plan item 3). Harness:
three transactional scenarios per shape (accept whole-matrix, accept
single-column, rollback) with the engine's accepted verdict recorded
as a channel; multinomial accept jitter 0.002 (0.005 rolls back at
K=3x40 - the count law visible in the fixture). Baselines
equivalence/bcf-equivalence/multinomial-equivalence-938eb81.rds:
neutrality partition first, self-reproduction 33/33 strict, 9/9, 7/7;
equivalence.yaml re-pinned. Gates double-run (implementer +
independent verifier): install --preclean, tests/cpp plain + ASAN
clean, tinytest 3985/0 no snapshot regenerated, trio + oracles
(bcf-exact quick 1e-4, restricted 2e-4), air 0, lintr no new
findings.

S2 LANDED 3cab553 (engine) + 5acd116 (baselines), 2026-08-12.
UpdateSessionImpl builds an explicit (chain, forest, tree) survivor
table in its constructor, pruned to the trees splitting on the
session's column and resolved through the promoted const
Chain::treeInForest (Chain::tree keeps its forest-0 meaning); one
predicate, Chain::collectSplittingTrees, feeds both consumers -
collectSurvivors (forest-0 exempt, the rebuild's
arithmetic-preservation rule) and the session's
treesSplittingOnColumn (not exempt) - so the guarded set is a subset
of the revalidated set in every forest, forest 0 included, as the
Pruning section mandates. Bridge:
refuseMultiForestPerObservationUpdate deleted; the retirement is
R-bridge ONLY - deviation from the plan's both-surfaces wording,
verified sound: the helper had no C_interface.cpp caller and
dbarts.h exposes no per-observation entry, so there was no flat
clause to drop; both per-observation entries now call
refuseVarianceForestPredictorMutation directly (S3 retires it); an
R-level refuseBCFMutation on the partial path removed. Joint entry
untouched (one permutation from samplers[0]->rng(), all-or-none
install, ANDed finalize returns). Fuzz: OP_PER_OBS +
OP_SESSION_ABANDON in fuzzMultiForestMask across all four
BCF/multinomial configs; snapshot discipline held (codes +
cutPoints/currentSampleNum compared beside statesAgree;
leaf-assignment capture at all three rejected/abandoned sites).
Falsifiers: F5 green at 110000 caller-ordered decisions (109456
installed / 544 declined) against an in-test per-tree leaf-count
oracle over every tree of every forest, negative half red on both
the mask and the finalizes; F6 green (82 pruned trees bitwise
including raw indices and fit slab; forest 0's unpruned round trip
moves totalFits on ~394 rows, reported not asserted, per the
invariant asymmetry), negative half reddened F2 on exactly the four
multi-forest configs; F7 green, negative half run in full - the
cache pinned to forest 0 with the revalidation widened raises the
unrecoverable "produced a tree with an empty leaf" hard error on a
multinomial sampler after cells were written, where the shipped
build installs 190-196/200; F8 green. Pins inverted in place; the
R5 pin asserts mask shape + finite run (that sampler is never run,
its trees are stumps, no row can empty a leaf - the decline half is
pinned on the burned-in handles in test-bcf-mutation-pins.R and
test-multi-forest-seam.R). Harness: per_observation 200/200 +
per_observation_partial 197/200 (BCF), k3perobs 198/200 +
k3perobspartial 194/200 (multinomial); the verdict channel is the
whole install mask, the session answering per row. Baselines
equivalence/bcf-equivalence/multinomial-equivalence-3cab553.rds:
neutrality partition first (33/33 identical draws, 9/9, 7/7 vs
938eb81, no max |z| anywhere), then self-reproduction 33/33 strict,
11/11, 9/9; equivalence.yaml re-pinned. Budget: ~820 changed lines
vs the ~440 estimate - engine 126 and bridge 51 within intent, the
overage entirely the mandated F5 oracle, F6 capture, fixtures, pins
and the harness lines item 7 mandates but the estimate omitted.
Gates double-run (implementer + independent verifier, fresh
privlibs): install --preclean, tests/cpp plain + ASAN/UBSAN clean,
tinytest 3993/0 no snapshot regenerated, trio + oracles (bcf-exact
quick to 4 dp, restricted 2e-4), air 0, lintr 0 new. Carried to S3:
Chain::treesSplittingOnColumn takes a forest index, so the variance
arm needs its own accessor beside varianceTree(j) (no
numTreesInForest entry); the heteroscedastic config's "for the
whole arc" comment in test_fuzz.cpp is still the pre-S2 text and is
S3's to edit; MultiForestFixture is the natural host for the
heteroscedastic F5/F6 configuration; the statesAgree
above-chain-level gap stands. Carried to S4 (residual coverage):
the bcf/multinomial equivalence legs are pinned by no CI workflow -
local-only gates, a pre-existing shape, not introduced here.

S3 LANDED a825263 (engine) + b174737 (baselines), 2026-08-12. The
E3 fallback was not needed: refreshVarianceForest split into
revalidateVarianceTrees + rebuildVarianceFactors for the
transactional paths - validate recovers node-indexed factors through
the LIVE partition first, repartitions, reports occupancy, collapses
nothing, drops no directions, scatters nothing; rebuild drops stale
directions, scatters survivors, recomputes combinedVariance - while
the collapse-and-remap form stays byte-for-byte for
forceRefreshTrees and applyNewData (hetforce bitwise gates it, F11
green). Deviation from the plan's signatures, verified sound: the
two halves take a ForestRevalidation handoff (one per chain, phase
one over all chains before phase two) plus the touched-column list
rather than a bare TreeParameters&, because the Pruning section's
requirement 1 mandates ONE survivor list travelling between phases
and the variance arm is pruned; the struct gained exactly two
members. Both calls appended after the forests_ body under the
runtime varianceForest_ guard; repartitionTrees gains the variance
arm; the session cache admits variance trees through promoted
Chain::varianceTree(j) + varianceTreesSplittingOnColumn (the
forest-index carry from S2). stateIsValid's variance branch gains
the scratch build, repartition and occupancy test its mean sibling
has - a setState behavior change with its own tinytest and NEWS
text; chain-level install stays ungated and the flatten identity
stays covered (test_model.cpp verdicts flipped in place). Bridge:
refuseVarianceForestPredictorMutation deleted AND
refuseMultiForestTransactionalUpdate deleted with it - the variance
clause was its whole body - so the flat C entries are unguarded
again BY DESIGN: the path SL's stop-loss protected is now legal, the
SL/F13 pins inverted in place at test-capi.R:167-215, and the flat
decline arm is reachable and pinned through a new test-support
capi_update_predictor_fixed_cuts wrapper in the tinytest consumer.c
(dbarts.h zero diff). Fuzz: the heteroscedastic config re-admits
OP_SET_PREDICTOR | OP_UPDATE_COLUMNS | OP_PER_OBS |
OP_SESSION_ABANDON, comment edited not flipped; snapshot discipline
held. Falsifiers: F9 green (5 rollbacks, 2 with the variance forest
sole objector on a mean-restricted fixture), negative half red on
exactly those 2; F10 green (12 accepted transactions, 22
observations changed variance leaves), negative half red on all 12;
F5 extension green at 215000 caller-ordered decisions, 105000 under
a variance forest, 0 mismatches, negative half red on mask + 413
finalizes; F6 extension green (21 pruned variance trees bitwise
incl. the factor slab, s2(x) bitwise where no variance tree is
reached), negative half reddened the routing invariant on the het
config; F12 untouched and intact. Harness: hetswap (verdict channel
15 accepts / 5 rollbacks over 20 seeds) + hetpartial (install mask
399-400/400, declines on 4 seeds), both carrying s2.test;
recordVerdict is per scenario so S0's summaries are unchanged. NEWS:
the SL "still refuse" sentences rewritten into the S3 bullet so the
release notes are not self-contradictory. Baselines
equivalence/bcf-equivalence/multinomial-equivalence-a825263.rds:
neutrality 33/33 (hetforce included) / 11/11 / 9/9 vs 3cab553, no
max |z| anywhere; self-reproduction 35/35 strict, 11/11, 9/9; the
bcf/multinomial re-records are pure re-anchors, stated in MANIFEST
with the reason (E = 3 unconstructible); equivalence.yaml re-pinned.
Budget: ~997 changed lines vs the ~640 checkpoint - engine 183,
bridge 92 (85 percent deletions), the remainder the plan-mandated
F9/F10 tests, het fixture, pin inversions, the item-5 setState
tinytest, harness scenarios and NEWS. Gates double-run (implementer
+ independent verifier, fresh privlibs): install --preclean, 0
compiler warnings, tests/cpp plain + ASAN/UBSAN clean, tinytest
4003/0 no snapshot regenerated, trio + oracles exact, air 0, lintr
0 new. Carried to S4: the dimnames drop on sampler$data@x after a
whole-matrix transactional setPredictor (pre-existing R-layer
behavior since S1) breaks by-name column resolution in
updatePredictorPerObservationJointly - the het pin orders the joint
call first as a workaround; recoverVarianceLeafValues' 1.0
fallback docstring still names setState as a reachable case and no
longer is, from any public route; the bcf/multinomial equivalence
legs remain un-pinned in CI; seven new comments cite slice/falsifier
identifiers (arc precedent from S1/S2) - whether S4 sweeps these
into self-contained rationale is an open hygiene choice; the
statesAgree above-chain-level gap remains worked around in the fuzz
snapshot rather than fixed in the helper.

S4 LANDED 031184c, 2026-08-12. Rd: dbartsSampler-class.Rd gains the
multi-forest/heteroscedastic mutation subsection (the per-sampler
conjunction, the install mask, the structural column-subset
exemption) and updatePredictorPerObservationJointly.Rd the joint
semantics; the sampler$data@x dimnames drop after a whole-matrix
setPredictor is DOCUMENTED with its exact blast radius (every later
by-name resolution fails - setPredictor's character column= and the
joint entry even under an integer column - while integer column= to
setPredictor itself survives) and workarounds; the fix is queued as
the new setpredictor-dimnames TODO entry, not shipped here. NEWS
needed nothing: the S1 commit already carried the plan's S4
sentence inside the S1 bullet - discovered, not decided. Design
docs: model-space-survey secs 4-5 + Status, heteroscedastic new sec
14 + Status, empty-leaf-veto records that the transaction enforces
the criterion stateIsValid already enforced and that S3 closed the
variance asymmetry (that file carries no Status line by its own
convention), c-api-growth records this arc added NO X-list entry
and moved NO hash. The VD-decided comment-anchor sweep landed:
about 20 comment/message-string edits across 13 files (chain.hpp,
sampler.hpp, common.cpp, test_fuzz.cpp incl. the check-message
falsifier tags, test_model.cpp, test_state.cpp, consumer.c, five
tinytest files) rewritten into self-contained rationale,
deliberately scoped to THIS arc's citations - sibling plans reusing
the S/F vocabulary keep theirs for their own arcs; the independent
verifier read every hunk: comment/string-only, no statement
touched, no technical claim lost, and the bitwise trio corroborates.
recoverVarianceLeafValues' docstring corrected (its 1.0 fallback is
unreachable from every public route since S3). TODO block landed
with the drift corrections (multiforest-mutation-gaps records the
flat guards' SL-to-S2 lifetime and S3 retirement; fork 2 recorded
settled per-sampler; bcf-public-surface pin citation verified to
:326-327) plus the statesAgree-gap record. Residual coverage: the
CI equivalence net widened - bcf-equivalence and
multinomial-equivalence jobs added parallel to the gaussian job in
equivalence.yaml, pinned to the a825263 baselines, compare commands
validated locally verbatim by implementer and verifier both.
Runner execution is DEFERRED, not confirmed: GitHub registers
workflows from the default branch, equivalence.yaml has never been
on main, so the workflow is neither schedule-fired nor dispatchable
from bartcore (checked at this landing - gh workflow list does not
know it); its first real execution, gaussian job included, comes
when bartcore lands on main. This is the workflow's own documented
shape, not a defect of this slice. Verification ran two rounds: the verifier's first verdict was
DO-NOT-LAND on four documentation defects (a phantom
updatePredictor R method named in the new Rd text and echoed in
three design docs, a false claim that setPredictor's by-name
column= survives the dimnames drop, a wrong TODO pin citation,
Status lines bumped by addendum rather than literally), all fixed
in the amended commit and re-checked (parse_Rd clean, clean-copy
tarball R CMD check --as-cran Status OK both rounds). Gates:
install --preclean, tests/cpp plain + ASAN/UBSAN clean, tinytest
4003/0, trio bitwise 35/35 strict / 11/11 / 9/9 vs a825263 with no
re-record (rng-neutral, and the trio doubles as proof the sweep
moved no draw), air 0, lintr 0.

THE ARC IS COMPLETE. Landed: SL 7299b8b; S0 1357e7d (+ records
430d4ee); S1 938eb81 + 7a9c6f3 (+ records 2e5114b); S2 3cab553 +
5acd116 (+ records e8d4812); S3 a825263 + b174737 (+ records
f5284cf); S4 031184c. Every predictor-mutation refusal is retired:
BCF, multinomial and heteroscedastic samplers accept transactional
and per-observation predictor mutation, joint variant included,
under the per-sampler empty-leaf conjunction with per-forest
j-split pruning, and setState validates variance occupancy.
