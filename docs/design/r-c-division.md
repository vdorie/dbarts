# The R / C++ division

Status: ACCEPTED (VD 2026-08-11), text amended (the text in "The
principle" below). Arc record, 2026-08-11. Artifacts (gitignored):
`.claude/r-c-division-design/{memo.md,critique.md,guide-memo.md}` - a
channel-by-family census memo, its adversarial critique (fifteen probe
scripts, run against a private build of b70b373), and a guide-first
memo written after VD twice corrected the framing ("a principle ...
should cover what we aspire to and not what we have currently
implemented"; "grounding the principle in code ... misses the entire
point of a principle which is as a guide"). The census and its
empirical claims were adversarially verified; the guide memo's two new
measurements are independently corroborated by the critique's probes;
the guide memo's own Part B demand census has since had its own
refuting critique (2026-08-11) - verdict NOT SETTLED AS WRITTEN,
amended below rather than treated as settled fact.

## VD's conjecture, adjudicated

The question: "any model that is built on mutating a sampler can be
done in R, but I think it's the case that any residualized/
mean-structure model assumes a Gaussian."

**Half 1 - nearly true, adoptable as aspiration.** Eleven of twelve
family configurations accept response/offset/sigma/predictor mutation
(census, ~60 cells adversarially re-run, none mismarked). The one real
hole is multinomial, whose response side is sealed (the response is
the combiner-owned counts matrix); the fix is the committed
multinomial-counts-mutation plan. A bridge nit: grouped samplers
refuse setResponse at the bridge while GroupedResponse::setResponse
delegates correctly.

**Half 2 - true of exactly one of three residualization routes.**
- Swapping the response to a residual (y minus another block) is
  Gaussian-only by arithmetic, as conjectured - and today it is a TRAP:
  latent-family samplers ACCEPT out-of-support responses (see
  "Defects" below).
- The OFFSET channel is family-generic residualization, measured, with
  production precedent (rbart_vi's own R loop for probit and aft;
  stan4bart's bernoulli path). Route-(b) compositions reach oracle
  accuracy for probit, logistic and negative binomial (critique probes
  B/D, E3, F; the nbinom route uses the shipped per-sweep state-read
  of the dispersion, the same idiom bart2Negbin itself uses).
- Splitting the mean across K samplers is NOT Gaussian-only and NOT
  broken: a host that draws the latents against the COMBINED fit (or
  lets each block draw its own latents against cross-offsets) targets
  the same posterior. Proven twice independently: pure-R Albert-Chib
  gold-vs-split arms agree to <2% of a posterior sd (Welch p
  0.145-0.963, ten chains per arm), and a two-forest additive probit
  model built from two dbarts samplers reaches engine accuracy (RMSE
  0.3116 vs 0.2918) once the standard sqrt(K) leaf-prior rescale is
  applied. The census memo's contrary claim, and the two-sampler
  latent-correlation probe supporting it, are REFUTED (the correlation
  is below the same sampler's own consecutive-sweep autocorrelation -
  augmentation noise, not model error).

**What is genuinely engine-only, in any implementation:** couplings
that are NOT additive on a quantity a channel can carry (the
multinomial softmax margin needs a log-sum-exp over the other K-1
forests); moves that write engine state (the per-forest ASIS rescale);
likelihoods that do not factorise over observations given the forest
(AR-1/spatial/copula errors - leaf statistics stop being within-leaf
sums); and the CALIBRATION (below). Mixing loss from separate blocking
is a price, not a wall. CORRECTED (critique, 2026-08-11): the arc's
first pass reported "up to ~6x effective samples in the bad case,
level in the common one", which conflated RMSE against the engine
target with actual mixing loss; the corrected measurement is that
separate blocking moves the truth by about ~1.5x in the bad case,
level in the common one.

**The wall that remains is calibration, and it is measured.** A
structurally correct pure-R probit composition inherits its leaf prior
from the range of whatever vector the sampler was CONSTRUCTED on: a
16x sweep of that accident moves the implied leaf sd 16x and the
posterior sd of f by 4.6x with no error or warning; correctness is
recovered only near one lucky range. stan4bart's mvbart() routes around
the response channel deliberately - the response is fixed and only the
offset moves - and warns when the per-sweep offset drift threatens the
anchor (`stan4bart/R/mvbart.R:40-44`, `:262-266`). The composed model's
BLOCKING ports to R; its PRIOR does not, unless it can be named.

## The principle (ACCEPTED AS AMENDED, VD 2026-08-11)

> R addresses the conditionals; C++ addresses the integrand.
>
> 1. R's side. Anything expressible as "one of the engine's own
>    models, fit to quantities the driving R program can compute
>    between sweeps" - residuals, latent draws, offsets from other
>    model blocks - should be reachable from R. The R program owns the
>    outer loop, the residual algebra, the latent draws and every
>    non-forest block, and changes response, offset, weights,
>    predictors and sigma between sweeps, on each channel where the
>    family's own identification permits it (a refusal on model
>    grounds is part of the model). Selecting an engine-provided
>    family or prior by argument is also R's side.
> 2. C++'s side. Anything that changes what the engine integrates - a
>    family, link, leaf model, tree prior, hyperprior law or
>    split-rule representation - belongs to C++ when it ships:
>    contributed to dbarts, then selected from R by argument. So does
>    any move that writes engine state, and any likelihood that does
>    not factorise over observations given the forest. The principle
>    adds reach to R; it never argues for moving what already runs
>    well in C++ out of it.
> 3. The promise. On R's side dbarts aspires to POSSIBLE and CORRECT,
>    and is honest about PRACTICAL. A composition targets the same
>    posterior as the equivalent engine model, or differs in a way
>    that is named, measured and testable - never silently.
> 4. What that asks. R should be able to name the calibration it
>    composes against, rather than inherit it from a construction
>    accident. Every engine capability should have an R route that
>    reaches the same posterior or a stated, tested difference, even
>    if slower. Every R affordance aims to ship with a tested recipe
>    and a validator its author can run. Read/write symmetry is the
>    default: whatever the engine conditions on, the R program may
>    write; whatever the engine draws, the R program may read.

Amendments made in accepting the text (VD 2026-08-11, relative to the
proposed draft): "owns" removed from both halves of the opening split;
clause 2 is a HOME-WHEN-SHIPPED claim, not a capability claim (the
clbart counterexample - its tree MOVES are pure R - would otherwise
read as a violation), plus the explicit sentence that the principle
adds reach to R and never argues for moving what already runs well in
C++ out of it; "host-computed data" is unpacked into what the R
program can actually compute between sweeps (residuals, latent draws,
offsets from other model blocks); the voice is aspirational throughout
("should"/"aims"/"asks", never "must"/"obliges"); clauses 3-4 read
"same posterior OR a named, measured, testable difference" rather than
a bare "must target". Clause 1's channel-permits-it parenthetical
carries forward the orchestrator's original absorption of the
critique's finding that sigma and per-observation weights are refused
for latent families BY IDENTIFICATION - the guide memo's wording
survives it; the census memo's "every channel for every family" did
not.

Arbitrations behind it (guide memo A12): velocity vs calibration ->
SPLIT THE OBJECT (R owns structure, the engine owns calibration, so
the calibration must be nameable from R); contract surface vs demand
-> general channels only, model-specific entries refused, impossible
specs refuse at creation; performance vs prototyping -> speed is
almost never the reason to cross (measured R-loop overhead 1.11-1.43x;
the honest reasons are joint moves, calibration, concurrency,
non-factorising likelihoods); expertise vs failure modes -> ship a
VALIDATOR, not a warning; introspection vs contract -> read/write
symmetry as a rule (readers carry no legality matrix); release latency
vs correctness -> the R shadow is mandatory and allowed to be worse,
provided the difference is stated and testable.

## Demand headlines (guide memo Part B; tool-verified, adversarially
critiqued 2026-08-11 - verdict NOT SETTLED AS WRITTEN, amendments
folded in below; no fabricated citation found anywhere in Part B,
including an independent re-check of Sun and Song)

- 24 of 24 BART-family artifacts that ship a tree sampler wrote or
  copied their own; zero reused an engine as a library; in every case
  the new thing was an integrand change or a non-factorising
  likelihood. One author under experimental control: nbbart reuses
  dbarts while the change is data; clbart abandons every engine the
  moment the leaf stops being conjugate. CRITIQUE: NEAR-TAUTOLOGICAL
  as stated - the honest ratio is 24 of 27, not 24 of 24; the guide
  memo's own list carries four integrand-change counterexamples,
  dbarts among them; and clbart's tree MOVES are themselves pure R,
  which cuts against the headline rather than for it.
- Of four CRAN packages driving a dbarts sampler iteratively, two are
  statistically invalid (posterior means used as draws) - a 50% defect
  rate on the response channel among peer-reviewed methodology. The
  composition failure mode is real and nobody catches it today.
  CRITIQUE: the honest count is 1 of 3 true sampler-drivers, not 50%
  of four; both defects are the SAME mistake (posterior means used as
  draws), which is exactly what the composition VALIDATOR (Adoption
  slate, below) catches - not something a calibration fix would have
  prevented.
- The top UNMET affordance is a row subset that changes between draws
  (principal stratification, mediation, IV): princeBART builds six
  probit samplers and re-runs `setData` on all six every outer
  iteration - five of them with a changed latent-stratum `subset` -
  for 2000 default iterations at n.thin 20. (Its nine getFromNamespace
  fetches, of which four are used, hand-build a probit sampler with a
  chosen node prior; they are NOT the row subsetting, which uses the
  public `dbartsData(subset =)`. The earlier "eight unexported
  internals for it" and the separate "AOAS 2024 / 6000 iterations"
  attribution are corrected and withdrawn: re-verified first-hand at
  princeBART 0.2.0, no such distinct artifact was found.) CRITIQUE:
  REFUTED as an unmet affordance - Gaussian row subsetting already
  ships from R via zero weights (`dbartsSampler-class.Rd:178`); the
  real gaps are the empty-leaf veto counting zero-weight rows as
  occupied, the latent-family refusal on sigma/weight mutation
  (clause 1's identification parenthetical), and a missing
  method-list doc - all priced in the Adoption slate rather than a new
  mutation shape. (This still reprices the multiforest-mutation-gaps
  row-subsetting door, which was value-lifted but unscheduled on a
  single consumer - the demand survey found the class, not just a
  consumer; what it opens is the latent-family subset mask below, not
  a general row-subset arc.)
- Peer packages misdescribe dbarts (SoftBart claims to be the only
  embeddable BART; stochtree's feature table denies dbarts a C API);
  the capability gap is smaller than the perception gap. Recipes and
  documentation, not architecture, are the competitive answer.

## Adoption slate (VD 2026-08-11: ALL items pre-release, framed by
UTILITY - price sizes a budget, it never ranks or gates; every item
below lands before the 1.0-0 freeze; committed plans carry their own
budgets)

- Nameable calibration (DISCHARGED; the named TODO entry is gone, the
  outcome is `docs/design/nameable-calibration.md`, PARTIAL - creation
  half LANDED c2a7e89b, mid-chain half LANDED, flat-C half an item
  inside the dbarts.h reshape's S1 - with the plan at
  `docs/plans/nameable-calibration.md`): design
  before the dbarts.h reshape's S1 (it may leave a header
  footprint); the single highest-value item; sized at ~1593 lines
  across four slices (budget expected above ~300 at design time);
  the excess is the flat-C half, the leaf-model contract, the
  three creation refusals and the oracle.
- A latent-family subset mask (DISCHARGED; the named TODO entry is gone,
  the outcome is `docs/plans/latent-subset-mask.md`): SHIPPED, S0-S4
  (dc11a805, 6db22aee, 87d370ea, 8b047f8b, 93afd635). This was the
  row-subset door the corrected demand survey opened (see "Demand
  headlines") - Gaussian row subsetting already shipped via zero weights;
  `$setActiveRows` extends the same "row i is not in the data set this
  sweep" membership channel to gaussian, Student-t, probit, ordinal,
  logistic, negative-binomial, aft and multinomial (global only), the
  latent families zero weights cannot reach. princeBART guarded each
  `setData` block with `if (sum(<stratum>) > 0)`, so a stratum that
  emptied was SKIPPED and that sampler silently kept the previous
  iteration's data; and `dbarts_binary` widens dbarts's binary detection
  to admit an all-0/all-1 stratum (R/utils.R:191-200, comment "can happen
  in principal stratification") because a constant all-1 response is NOT
  detected as binary today. Both were defects the mask removes by
  construction - an all-zeros mask is accepted and runs rather than
  needing to be guarded around - and both were stronger demand evidence
  than the count that was wrong.
- The empty-leaf veto fix (DISCHARGED; the named TODO entry is gone, the
  outcome is `docs/design/empty-leaf-veto.md` - the veto is now a
  lexicographic rank over leaves that separates "no member" from "no
  positive-weight member" - with the plan at
  `docs/plans/empty-leaf-veto.md`): its own
  measured slice, a draw-law change - the veto counted LEAF MEMBERS
  where it should count POSITIVE-WEIGHT members
  (`src/bartcore/moves.hpp`, the `numObservations() == 0` veto site);
  expect to re-record the zero-weight baseline.
- A composition validator (DISCHARGED; the named TODO entry is gone, the
  outcome is `docs/plans/adoption-slate.md` S5 and its landing note; SBC
  over a user-supplied one-sweep closure): the item that catches the measured
  posterior-mean-as-draws defect class (see "Demand headlines"), not a
  calibration fix.
- Exported augmentation helpers (TODO `augmentation-helpers`, LANDED
  890efd3d as adoption-slate S4) and named recipes cross-referenced from
  `?bart` (surface shakedown; DISCHARGED, the named TODO entry is gone
  and the outcome is `docs/plans/adoption-slate.md` S6 and its landing
  note).
- An offset-free fit accessor and a seeding contract for composed
  samplers: already committed inside the getLatents docs slice
  (defect 5, below); confirmed here, not double-booked. The claim that
  motivated it corrects from "four independent codebases hand-write
  train - offset today" to three packages under one author - the
  demand is real but narrower than first stated.
- Six small census items, itemized: G1 (setResponse support
  validation), G5 (the weight-refusal message) and G11 (the
  variance-forest updateScale guard) are DONE at 33f6fdc (see
  "Defects"). G2 (getLatents location-vs-precision docs, the wrong
  "(gaussian)" parenthetical in dbarts.h) IS the committed getLatents
  docs slice (defect 5) - no new entry. G3: export the nbinom
  dispersion r per draw as a first-class surface (a run-result slot or
  getter, dbarts.h plus R5) in place of the bart2Negbin per-sweep
  state-read idiom; small pre-release slice, scheduled after the
  getLatents docs slice. G10: relax the bridge-only refusal of
  setResponse on a GroupedResponse sampler - the model already
  delegates correctly, a bridge nit this arc recorded; small
  pre-release slice, same era. No VD fork: every remaining item
  carries a surface and a nameable value, so all schedule pre-release
  under the recorded frame.

## Fork implications (docs/plans/multiforest-extension-surface.md,
Open decisions) - ALL FOUR RESOLVED (VD 2026-08-11)

All four recommendations ADOPTED; fork 1's scheduling changed from the
recommendation, and one other justification changes:

1. Basis family: BUILD, PRE-RELEASE - not post-freeze as recommended.
   Scheduled after the dbarts.h reshape (M3 carries its header entries
   into reshape S1). The "non-Gaussian multiforest is otherwise
   unreachable from R" premise is falsified by measurement, so the
   ~300-line nameable-calibration affordance is built first (Adoption
   slate, above) and FA1/FA2 (the amplitude/ASIS mixing gates) become
   design-informing probes rather than go/no-go gates.
2. Non-Gaussian v1: probit+logistic, ADOPTED - the reason is
   correct-by-construction calibration vs correct-by-care, not
   impossibility. AMENDMENT REQUIRED at M4 scheduling: falsifier FA5
   as committed tests a strawman (K independent probit samplers); the
   decisive arm is K GAUSSIAN samplers with host-drawn latents against
   the combined fit, which measurement predicts will AGREE.
3. Flat rename: ADOPTED, re-sign at the dbarts.h reshape re-bake to
   setForestBasis/numForestAmplitudes/forestAmplitudes, + KEEP
   setForestWeights; stakes lowered (the flat C surface has exactly
   one consumer and it is in-house - 1 reverse LinkingTo vs 23
   R-level consumers).
4. bcf() home: bartCause (its dbarts-1.0 branch), ADOPTED;
   discoverability is answered by documentation, not vocabulary at the
   engine's front door.

Layering, sharpened: stan4bart is where the OTHER block gets a better
sampler; dbarts is where the FOREST does; everything else is a recipe,
and recipes are dbarts' documentation, not anyone's C++.

## Defects found by this arc (live at da37242; all adversarially
confirmed)

1. MEMORY-UNSAFE: setResponse on a negative-binomial sampler with a
   single negative element SEGFAULTS the R process (uncatchable;
   static_cast of a negative lround underflows to ~1.8e19 and sizes a
   histogram); y = 1e9 hangs on unbounded allocation. Same hole on the
   flat C path. Creation validates what mutation accepts. FIXED at
   33f6fdc: `bartcore_bridge::validateResponseSupport`, shared by both
   surfaces at creation, setResponse and setData, refuses the negative
   element. DELIBERATE DEVIATION: no nbinom magnitude cap - creation
   imposes none either, so y = 1e9 still hangs, on both creation and
   mutation; a magnitude-cap or sparse-tally decision was left a new open
   item, VD-facing, unscheduled. DISCHARGED: the named TODO entry is gone
   and the outcome is `docs/plans/adoption-slate.md` S7 and its landing
   note - a derived nbinom magnitude cap of 1e6 inside
   `validateResponseSupport`'s nbinom arm, so one edit covers creation and
   every y-swapping conduit on both surfaces
   (src/R_interface_bartcore_common.hpp:190).
2. SILENTLY WRONG: probit/ordinal setResponse accept out-of-support
   responses (non-0/1 y gives latents in the hundreds; constant 0.5
   collapses every latent to zero). setTreatment already shows the
   correct pattern (validates support post-creation). FIXED at
   33f6fdc: the same `validateResponseSupport` guard.
3. Variance forest + setResponse/setOffset(updateScale = TRUE): the
   sigma-pin is bypassed and the fit RUNS AWAY (mean variance 156 ->
   0.70 against truth ~171 over 750 sweeps) while getSigmas() reads
   unchanged - a fifth sigma door beyond variance-forest-mutation-
   routing's four, invisible from R. A refusal keyed on
   hasVarianceForest (mirroring the BCF one) needs no algebra. FIXED
   at 33f6fdc: `refuseVarianceForestScaleUpdate`, on setResponse and
   setOffset, both surfaces.
4. refuseBinaryWeightChange tells aft/ordinal/nbinom users about "a
   binary response" (~6 lines). FIXED at 33f6fdc: reworded per family.
   FOLLOW-ON CLOSED: the flat C `dbarts_sampler_setWeights` took the
   guard with its promotion into `bartcore_bridge`, and the value
   policy with the logistic channel below. The harm the follow-on
   stated was wrong either way - a non-positive logistic count never
   divided by zero, since the count is a loop bound AFTER one
   unconditional PG draw, so such a row carried a full PG(1, psi)
   precision instead: a phantom observation, not a crash.
5. getLatents is family-polymorphic with no documentation (probit/
   ordinal/aft return a LOCATION; logistic/nbinom/Student-t a
   PRECISION; dbarts.h's "(gaussian)" parenthetical is wrong), and
   run()$train CARRIES THE OFFSET - the natural reading of
   "getLatents() - train" as a conditional draw biases slopes by 7-10
   oracle SEs; correct usage is incremental. Doc + Rd + dbarts.h fix,
   plus the offset-free fit accessor. STILL OPEN: the one committed
   docs slice this arc's stop-loss landing did not close.
6. The multinomial bart2 host's $fit is a placeholder sampler that
   accepts every mutation as a silent no-op (the real K-forest sampler
   rides $bc); needs a guard or prominent documentation until the
   multinomial mutation arc replaces the area. FIXED at 33f6fdc for
   multinomial: a `hostFor` field plus `refuseHostMutation` on the R
   sampler class. NEW OPEN ITEM: bart2's ordinal and nbinom hosts
   retain the same disconnected-$fit shape; whether to extend the
   `hostFor` guard to them is a VD question BY NAME, unscheduled.
   CLOSED: the multinomial mutation arc replaced the area, as the
   guard's own comment anticipated. S1a pointer-adopts ordinal/nbinom's
   engine; S4 constructs multinomial's sampler directly. $fit is the
   sampler that ran for all three, $bc is gone, and the `hostFor`
   field/guard is deleted - nothing is left to guard.

Recommended sequencing, as landed: items 1-4 and 6 were one stop-loss
slice at 33f6fdc (the SL/VF S1 pattern: refusals first, no draw-law
change on valid paths); 5 remains a docs slice, the one item that
commit does not close.

## Errata this arc generates elsewhere

- docs/plans/multiforest-extension-surface.md's anchor-refresh pass is
  DONE, this commit (31 corrections applied; three stale line numbers
  its own header claimed were verified - :146, :360, :970 - plus the
  chain.hpp/combiner.hpp anchors drifted by the S4 landing, incl.
  fact 12's latent draw at chain.hpp:1125-1126, not :1109-1112).
- The census memo's tier (iii).1, its fork-1/fork-2 "strengthened"
  justifications, its two-sampler latent probe, and its "nbinom
  dispersion unreachable" gap are superseded per the critique
  (recorded here; the artifact is gitignored and left as written).

## The latent-family weight channel (G4), adjudicated 2026-08-15

The census's G4 - "decide the latent-family weight channel: either
build it or record a considered decline; leaving it as an accident is
what the principle forbids" - was never itemized into the adoption
slate and never adjudicated. It is adjudicated here as a PARTIAL
BUILD, ticketed outside that arc (TODO `latent-family-weight-channel`).
Logistic is built: its weights are observation counts, the
Polya-Gamma augmentation's shape is the sum of that many PG(1, psi)
draws, so a weight swap is a model change with a defined meaning and
varying exposure is a capability the channel would add. Probit,
ordinal, aft and nbinom are DECLINED BY IDENTIFICATION - a weighted
probit has no tractable latent form, ordinal inherits that, aft fixes
its censoring structure at creation, and nbinom's shape parameter is
y_i + r with no weight slot. Clause 1's identification parenthetical
covers all four, and this is the record that makes the coverage
deliberate rather than inherited.

LANDED d0701a6a (2026-08-24), logistic build. A swap replaces omega (the
count is the PG SHAPE, so a stale one is a WRONG draw, not an old one),
the working response, the masked composite, and the borrowed weight
pointer - the last making the family override mandatory, since the bridge
frees the previous vector as it retains the new one. `GroupedResponse`
delegates the swap against f + b without drawing its own group effects or
tau block, those staying a sweep concern. The refresh is IMMEDIATE: a
sweep draws its TREES first and refreshes after, so an outer sampler at
run(0, 1) per iteration would otherwise move every tree under the
previous counts. It draws from the CHAIN generators in a serial per-chain
fan-out, so consumption is independent of n.threads and R's stream never
moves. Under a mask only ACTIVE rows are drawn for - variates for an
inactive row would desynchronize the stream against a sampler built on
the retained rows - and the rest return to the deterministic cold start
against the NEW counts, so none reactivates on a shape the sampler no
longer holds. A zero count stays refused for a better reason than first
recorded: `leafVetoRank` counts positive WORKING weights and a zero-count
row still carries omega > 0, so it does not drop the row a zero gaussian
weight drops; the mask is the only spelling that does.
`enforceBinaryWeightPolicy` states the admissible values once, at
creation and on both mutation conduits, so the R bridge and the flat C
entry cannot disagree, and `setData` hands the replacement counts through
the SAME conduit, so its latents are drawn rather than cold-started;
replacement data given without weights is single-trial, as at creation.
CLOSED since: weights still do not ride the saved state, but a digest of
them does, so `setState` sees a count vector other than the one the stored
omega was shaped by and re-derives the latents through this same conduit.
The state seam now agrees with the live seam on the operation both serve;
a matched round trip re-derives nothing and stays the identity.

## Honest limits

The decisive composition experiments are single-DGP, single-machine
measurements - decisive about mechanisms (the augmentation composes;
the calibration does not; the runaway diverges), not about magnitudes.
The demand census counts what was findable on CRAN/GitHub/arXiv. The
IACT figures carried from earlier arcs keep their recorded caveats
(none is a forest-vs-forest measurement). The guide memo's Part B has
now been adversarially critiqued (2026-08-11): verdict NOT SETTLED AS
WRITTEN, amended as recorded in "Demand headlines" and the "Adoption
slate" above; its Part A load-bearing claims were independently
corroborated by the critique's probes before this pass.
