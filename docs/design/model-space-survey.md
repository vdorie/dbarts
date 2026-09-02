# Model-space survey: which model classes justify which update shapes

Code citations are at 60f899b4; they are not live.

Status: COMPLETE as a survey, 2026-08-08; D1 (section 4) subsequently LANDED
2026-08-10 through 2026-08-12 (docs/plans/archive/multiforest-predictor-mutation.md,
SL a13152c5 through S4) - see the landing note appended to section 4; D4
(section 4, the general per-forest multiplier channel) ADDED and CLOSED
2026-08-14 (docs/plans/archive/multiforest-extension-surface.md M4.0-M4.3), carrying
the four model classes the survey lacked as multiplier classes. Research
only - no code changed. Run
as the research gate VD placed in front of the two multi-forest update-shape
doors (TODO `multiforest-mutation-gaps`; docs/plans/archive/runsbcbcf-repair.md "setData
door survey"), on the standing rule that **the shape follows from the models**.
Carries the `backlog-value-scan` TODO entry as section 7. Method was design memo
-> blind refuting critique -> this synthesis; the critique's verdict was STANDS
WITH AMENDMENTS, and where memo and critique conflict this record follows the
critique except at two line-number claims re-checked here (section 8). Working
papers are untracked: `memo.md`, `critique.md`, `synthesis.md`.

Summary: door 1 (multi-forest whole-data `setData`, n free) and door 3
(`setData` on CSC/mixed stores) both stay UNDESIGNED - no verified model class
asks for either, and the growing-n absence claim, once corrected, still supports
the KEEP. Door 2 (per-forest row subsetting) has its VALUE gate lifted by a
single verified class (Sun & Song's mixture-cure BCF, whose authors forked the
`bcf` package to get it) but its SHAPE is UNRESOLVED: the memo's "mutable row
mask" recommendation rested on three classes that do not survive checking, and
the one class that does implemented physical compaction. The survey's most
valuable finding is not any of the three doors - it is a NEW door, D1
(multi-forest transactional and per-observation predictor mutation at fixed n),
which carries an existing consumer and four model classes and is strictly
cheaper than door 1. Two further new doors (D2, D3) and one memory-safety defect
came out of the same pass.

---

## 1. The question, and the method

The doors on record were all stated in terms of *data replacement*. The question
put to the survey was not "can we build them" but "which model classes are worth
enabling", so that the shape would be derived rather than guessed.

Method: (1) read the three door records and the surrounding engine code to
establish what the mutation surface already serves, because a door's value is
only its delta over that; (2) enumerate candidate classes from the causal/BART
literature, from the ecosystem dbarts already serves (stan4bart, bairrtt,
bartCause, treatSens, rbart_vi), and from the repo's own design docs, several of
which had already argued the model-space question without being read as survey
input; (3) per class, ask what update shape it needs, whether dbarts expresses it
today and how awkwardly, and its value under an R-user-need lens; (4) map the
surviving needs back onto the doors and record the shapes no model asked for.

Two methodological findings, both load-bearing:

- **Every model that wants these doors wants mutation of *membership* or of
  *values* at fixed n, not replacement of the data.** The doors were stated in
  the wrong currency. That single observation drives door 1's and door 3's
  verdicts and produces D1.
- **A model's structure is not its MCMC.** The memo read four classes' shapes
  off their abstracts; the critique read two of them at code level and found the
  implied update shape different in three cases out of four. Abstracts identify
  the *model*; only the methods section or the source identifies the *update*.
  Any future model-space pass must reach the sampler, not the summary.

---

## 2. What the engine serves today

The baseline every door is a delta over. Verified by reading at 60f899b4.

**Single-forest samplers.** `setData` already replaces predictors, response,
weights, offset and test data **with n free to change**
([[src/bartcore/sampler.hpp:1046-1056@60f899b4]], "with a possibly different
number of observations"; the predictor count is fixed). It rebuilds the cut
grid, remaps existing splits onto value-nearest new cuts via
`Tree::mapOldCutPointsOntoNew`, and collapses whatever is left invalid
([[src/bartcore/chain.hpp:1553-1611@60f899b4]]). Dense stores only
([[src/R_interface_bartcore.cpp:3346-3351@60f899b4]]); grouped and AFT
samplers refuse ([[src/R_interface_bartcore.cpp:3314-3319@60f899b4]]).

**Multi-forest samplers** are exactly BCF (2 forests) and multinomial (K
forests): `numForests` is `forests_.size()`
([[src/bartcore/chain.hpp:711@60f899b4]]). The heteroscedastic variance forest
is a *separate* nullable member `varianceForest_`, so a heteroscedastic
sampler reports `numForests == 1` and no multi-forest guard sees it - see
section 6.

- `setData` refuses for `numForests >= 2`
  ([[src/R_interface_bartcore.cpp:2649-2654@60f899b4]]).
- Response-side mutation rides `supportsResponseMutation`: **true for BCF**
  ([[src/bartcore/combiner.hpp:621@60f899b4]]), **false for multinomial**
  ([[src/bartcore/combiner.hpp:411@60f899b4]] default, not overridden).
- `setTreatment` swaps z mid-run at fixed n
  ([[src/bartcore/combiner.hpp:439@60f899b4]]).
- `setCutPoints` installs an arbitrary grid with no multi-forest guard.
- Predictor mutation: the **forced** whole-matrix `setPredictor` is supported
  multi-forest and is the documented swap; the **transactional** and
  **per-observation** paths refuse on `numForests >= 2`
  (`refuseMultiForestTransactionalUpdate`,
  [[src/R_interface_bartcore.cpp:1929-1943@60f899b4]]),
  because `Chain::revalidateTrees` opens `Forest& forest = forests_[0]`
  ([[src/bartcore/chain.hpp:1484@60f899b4]]) and revalidates the primary
  forest only.

**Row subsets already exist, at creation.** `bartcore_createFromHandle`
([[src/R_interface_bartcore.cpp:2609-2620@60f899b4]]) builds a sampler over a
row-subset view of a data handle, copying the handle's cut grid so folds bin
identically to the full data. That is xbart's fold mechanism and the mechanism
docs/plans/archive/data-ownership-4-views.md names for "hurdle/IRT-style embeddings". The
subset is **fixed at creation**; such a sampler refuses `setPredictor`,
`setData`, `setCutPoints` and `setState`.

**Zero weights are already exact.** A zero-weight row leaves the weighted SSR and
every leaf conditional, and since docs/plans/archive/sigma-df-zero-weights.md the sigma
posterior df counts only positive-weight rows. "Exact exclusion" as a
*likelihood* concept is served. What is not served is exclusion of a row from
*one forest of a coupled model*, and the O(|S|) cost that would come with it.

**The BCF near-exclusion was a floor, not a zero - CORRECTED, it is now an
exact zero.** As surveyed, `formForestResponse` divided the residual by the
forest multiplier and floored |m| at 1e-9, so a b0 = 0 control row entered
tau's suffstats with weight ~1e-18 and response ~1e9 x resid: effective, not
exact. The floor was REPLACED by an exact-zero snap at `0x1p-26` on 2026-08-10
([[src/bartcore/combiner.hpp#zeroMultiplierTolerance, formForestResponse]];
docs/plans/archive/zero-weight-exactness.md), the shipment door 2's own update
below records at "Update, 2026-08-10". A snapped row now leaves with exactly
zero weight and exactly zero response. What the snap does NOT change is the
cost: tau still pays O(n) per sweep, which is the half of this paragraph that
door 2 rests on.

---

## 3. Door verdicts, as amended

### Door 1 - multi-forest whole-data `setData`, n free, z required

**KEEP UNDESIGNED.** No model surfaced for the n-free multi-forest joint swap.

- Every multi-forest class the survey verified holds **n fixed** and mutates
  *values* or *membership*: latent-confounder/latent-treatment sensitivity
  [verified: https://pubmed.ncbi.nlm.nih.gov/27139250/], principal
  stratification with BCF [verified: https://arxiv.org/abs/2403.13256], the
  mixture-cure AFT BCF
  [verified: https://projecteuclid.org/journals/bayesian-analysis/advance-publication/A-Tree-based-Bayesian-Accelerated-Failure-Time-Cure-Model-for/10.1214/23-BA1402.full],
  and panel/longitudinal BCF (LongBet [verified: https://arxiv.org/abs/2406.02530];
  longitudinal BCF [verified: https://arxiv.org/abs/2407.11927]). None grows or
  shrinks the row set mid-chain.
- The one verified growing-n class - BART as a surrogate for sequential design of
  computer experiments [verified: https://arxiv.org/abs/1203.1078] - is
  single-forest and is **already served** by the shipped n-free single-forest
  `setData`.
- **Corrected absence statement.** The memo claimed no BART or BCF instance of
  adaptive-trial / response-adaptive-randomization design exists. That is FALSE
  as written and is narrowed here to its growing-n form. Three instances exist:
  "a covariate-adjusted response adaptive randomization, updating treatment
  allocation probabilities grounded on causal effect estimates using a random
  intercept accelerated failure time BART model"
  [verified: https://pubmed.ncbi.nlm.nih.gov/38688389/] - i.e. built on
  riAFT-BART, the same CRAN package this survey cites elsewhere; BFTS, "the first
  contextual bandit algorithm to integrate Bayesian Additive Regression Trees
  (BART) ... directly into the exploration loop", evaluated on a
  micro-randomized trial [verified: https://arxiv.org/abs/2602.07767]; and
  DOD-BART, "a novel seamless phase I/II design ... automatically updated with
  emerging data to allocate patients to the most promising dose levels"
  [verified: https://pubmed.ncbi.nlm.nih.gov/39611660/], with a 2025 successor
  DOD-PRO-BART. **The verdict survives and is strengthened**: none of them
  maintains a single fit while n grows - they refit at each interim or refresh,
  which needs nothing from dbarts. BFTS is affirmative evidence *for* the KEEP:
  "At each refresh, we use a cold-start MCMC and re-initialize the sampler from
  the prior to prevent the chains from becoming trapped in local modes"
  [verified: https://arxiv.org/html/2602.07767v1]. A leading BART adaptive-design
  method throws the fit away *on purpose*. The narrow claim that survives:
  **no verified BART or BCF model maintains one fit across an accruing sample.**
- The z-required constraint remains correct if the door is ever taken
  ([[src/bartcore/combiner.hpp:467@60f899b4]],
  [[src/bartcore/combiner.hpp:486@60f899b4]],
  [[src/bartcore/combiner.hpp:510@60f899b4]] index `glue_.z` over live
  `numObservations` - `glue_.z` was since deleted at M4.3), but no
  model was found to pay for it.

If any part of door 1 is ever taken, take D1 and D3 (section 4) first: strictly
cheaper, and they carry the models.

### Door 2 - per-forest row subsetting

**VALUE GATE LIFTED. SHAPE UNRESOLVED.**

*Value.* One verified class wants it: Sun and Song's tree-based Bayesian AFT
**mixture cure** model, which joins a BCF to an AFT cure model. The paper's body
confirms the shape - the cure-rate component is fit "based on the entire sample",
the latency components "based on the uncured subgroup (across the iterations)",
with the group label latent and unobservable for censored subjects
[verified: https://projecteuclid.org/journals/bayesian-analysis/advance-publication/A-Tree-based-Bayesian-Accelerated-Failure-Time-Cure-Model-for/10.1214/23-BA1402.pdf].
The per-sweep draw itself lives in the paper's Web Appendix, but is explicit in
the implementation: `glabel[k] = gen.binom(1, p_post)` per censored row, every
sweep [verified: https://github.com/roxiesun/TBAFTcure src/TBAFTcure.cpp].

To get it, the authors wrote a **standalone C++ implementation** whose tree code,
by their own README, is "constructed following the framework of R package bcf" -
they reimplemented rather than reused. (The critique called this a fork of `bcf`;
it is not, and the distinction matters only in that reimplementation is the
*stronger* evidence of an unmet need.) Under the standing enabling-value rule one
such consumer is sufficient to lift the value gate on its own.

*What did not survive.* The memo claimed four classes, "one shape". Three fall
away under checking, and this is the survey's largest correction:

- **Principal stratification with BCF** imputes the potential intermediates per
  sweep - "the value of the principal stratum membership of the r-th MCMC
  iteration" - and enters them into the outcome ensemble **as splitting
  covariates**: "The two potential intermediates ... along with X_i, serve as
  independent variables for splitting the tree during the growth of the tau_Y
  function." Every component sees all n rows every sweep
  [verified: https://arxiv.org/html/2403.13256v1; code
  https://github.com/lit777/BPCF]. That is the mutable-predictor-column channel -
  **door D1, not door 2**.
- **Zero-inflation / log-linear BART** keeps every row in every component of the
  augmented likelihood. Its likelihood is a product over all n, and the node
  sufficient statistics `r_ht = sum u_i`, `s_ht = sum f(x_i) v_i` run over all
  rows filtered only by node membership; in the ZINB augmentation a row assigned
  to the point-mass component has Z_i = 0, hence u_i = v_i = 0, hence contributes
  exactly nothing - **a zero weight, not a removed row**
  [verified: https://arxiv.org/abs/1701.01503, eqs. 4, 5, 12, 14; this is a
  derivation from those equations, not a verbatim statement of the paper]. The
  memo's self-declared strongest citation argues the existing weight channel
  suffices.
- **Latent-treatment sensitivity** was double-counted: the memo itself routes it
  to D1 ("needs no new data shape - only multi-forest revalidation") and then
  lists it among door 2's four. The pinned record is right and the memo was
  wrong here: that class **runs today** at fixed n via the fixed glue
  (a = 1, b0 = 0, b1 = 1) plus per-iteration `setTreatment`
  (docs/plans/archive/runsbcbcf-repair.md). Its residual is **efficiency** - O(treated)
  rather than O(n) tau cost - not mutability. (As written this also claimed
  "exact exclusion instead of the 1e-9 floor"; the floor became an exact zero on
  2026-08-10, so only the cost half survives. See section 2.)
- The IV/CACE "neighbour" cited as support has **no compliance classes at all**
  (no principal strata, no latent classes; it is a two-equation continuous-exposure
  IV with correlated errors, both ensembles on all n rows). Citation **dropped**
  from this record.

So door 2's support is one class, and its content for the *other* interested
class is efficiency, not mutability.

*Shape: three candidates, none settled.*

1. **Mutable per-forest row MASK.** A per-forest membership set, settable per
   sweep in the `setTreatment` idiom. Excluded rows keep their identity: full-n
   reporting and prediction, the combiner keeps indexing over live n, no
   renumbering anywhere. Argument for: the subset forest must still predict at
   all n (a cure model needs the latency surface at every row), so a renumbered
   view pays a scatter-back on every read. Argument against: the mask must be
   threaded through every per-observation kernel of that forest, and the O(|S|)
   win only materializes where kernels iterate an index list rather than test a
   flag. Cannot give the subset forest its own cut grid.
2. **Per-forest ZERO WEIGHT.** Excluded rows get weight 0 in that forest's
   `ForestResponse`. By far the cheapest, and SHIPPED since this was written
   (2026-08-10; see the update below): the combiner already formed a per-forest
   weight vector every sweep (now [[src/bartcore/combiner.hpp#formForestResponse]])
   and needed only an
   exact-zero path in place of the then-live 1e-9 multiplier floor, plus a
   caller-settable per-forest weight. Both landed. Exact exclusion from the leaf
   conditionals and the sigma df is an already-solved semantic
   (docs/plans/archive/sigma-df-zero-weights.md).
   This is the shape the surviving zero-inflation construction actually uses.
   Does **not** buy the O(|S|) cost reduction - the forest still pays O(n).
3. **Physical COMPACTION.** The subset forest holds a gathered |S|-row view,
   rebuilt when membership changes; predictions scatter back to n. Buys the
   O(|S|) cost honestly and gives the subset forest its own row space (and, if
   cuts are rebuilt, its own grid). Costs a per-sweep gather of every
   per-observation field, a scatter-back on every read, and renumbering through
   reporting, varcounts and the state format.

*The only real-world datapoint favours compaction, with a caveat.* TBAFTcure
draws the cure label per censored row each sweep, then physically allocates fresh
`nuncured`-sized arrays with a gather loop and sets the data-info row count to
`nuncured`, restoring n only for prediction - and the **zero-weight alternative
is visible, commented out** in the same source
(`weight[k] = w[k]*mscale*mscale/(sigma*sigma);//0.;`, in both the treated and
control branches). They contemplated zero weights and replaced them with
compaction. The caveat that keeps this from settling the fork: TBAFTcure follows
`bcf`'s framework, and `bcf` has neither mask machinery nor a caller-settable
per-forest weight surface, so compaction may be the cheapest thing to build in
*that* setting rather than the right thing for *this* one. dbarts already ships
the per-forest weight channel and a creation-time row-subset view; the local
costs are not the same.

*A finding in its own right: the two interested classes may want different
shapes.* The cure model wants (3) or (1) - it needs the cost reduction and
possibly a subset-specific grid. The latent-treatment sensitivity class wanted
(2) - exactness against the then-live 1e-9 floor - and would take the cost
reduction as a bonus. That last question is ANSWERED: (2) DID ship alone, on
2026-08-10, as the cheap exactness fix, leaving the cost reduction and the
one-mechanism-or-two question open for (1) and (3).

*What would settle the shape,* in the order that decides the most per unit cost:

- **Does the subset forest want its own cut grid?** hurdle.md raised this and
  declined to force a shared grid; TBAFTcure's compaction implicitly gives the
  latency component a subset grid. If the answer is yes, the mask is eliminated
  outright and the fork collapses to compaction. Evidence: the cure paper's
  methods section, and whether TBAFTcure rebuilds cuts on the compacted arrays.
- **Cost at realistic membership fractions.** Kernel-level A/B of mask vs
  compaction vs zero-weight at 10-50% membership, the existing benchmark
  methodology. Nobody has measured this.
- **Empty-leaf accounting under changing membership.** A leaf non-empty under S
  can be empty under S'. Either the change is transactional (the `setPredictor`
  posture: roll back a change that empties a leaf) or empty-under-membership
  leaves must be defined (prior draw). A modeling decision, and the sharpest way
  to get the door wrong.
- **Which residual sigma sees.** A row excluded from forest f is *not*
  zero-weight in the likelihood - it still contributes through the other forests -
  so the exclusion must not reach the positive-weight df count, but a cure model's
  latency component conditions the residual definition on membership. Undecided.
- **What the exclusion does to the split prior**, and to per-forest varcount
  reporting, which under changing membership describes a sample-dependent
  subpopulation and so changes meaning per sweep.
- **Whether one class is enough to build for.** The value gate lifts on one; the
  shape decision does not have a second consumer to triangulate against.

**Update, 2026-08-10 (docs/plans/archive/zero-weight-exactness.md).** Shape (2), per-
forest zero weight, SHIPPED: the combiner's exact-zero multiplier snap (S1) and
a caller-settable per-forest, per-observation weight on
`AmplitudeForestCombiner` (S2, `dbarts:::`-only). Shapes (1) (mutable row mask)
and (3) (physical compaction) stay open and unscheduled - the shipment resolves
only shape (2)'s branch of the fork above, not the fork itself. One fact
recorded, not decided: because occupancy is count-based, per-forest zero weight
and physical compaction are NOT the same draw law even in the limit - under
compaction the excluded rows do not exist, so their cut positions cannot be
split on and their leaves cannot be created; under zero weight they can, at the
prior. The three shapes are therefore not interchangeable implementations of
one semantic.

### Door 3 - `setData` data-type shapes (whole-data replacement on CSC/mixed)

**KEEP UNDESIGNED.** No model surfaced, and the corrections strengthen it.

- The classes needing sparse designs (high-dimensional text, genomic,
  administrative) do not swap whole data mid-chain.
- Every class that *does* mutate mid-chain mutates a **dense** latent column or a
  dense label: latent U, sequential-BART imputation
  [verified: https://academic.oup.com/biostatistics/article/17/3/589/1744410],
  latent ability, imputed principal strata, the cure label. Mixed containers
  already carry dense columns alongside sparse ones. That is independently the
  reason the per-observation-CSC refusal stays legitimate.
- The remaining sparse-uniformity work has a home in the in-flight
  `cheap-uniformity` arc.
- The one live sub-case - whole-data replacement with n free on a sparse store -
  inherits door 1's verdict.

---

## 4. New doors the survey surfaced

### D1. Multi-forest transactional and per-observation predictor mutation, fixed n

`refuseMultiForestTransactionalUpdate`
([[src/R_interface_bartcore.cpp:1929-1943@60f899b4]])
refuses the transactional `setPredictor`/`updatePredictor` and every
per-observation session on `numForests >= 2`, because `Chain::revalidateTrees`
([[src/bartcore/chain.hpp:1484@60f899b4]]) revalidates `forests_[0]` only. The
forced whole-matrix swap is the only supported multi-forest predictor
mutation.

**Models blocked - four, all fixed n, all needing no new data object:** IRT +
causal forest (bairrtt, an **existing consumer**, whose `irt_causal_bart` runs
two samplers sharing a latent ability column updated by
`updatePredictorPerObservationJointly`, and whose natural next model is theta as
a treatment *moderator*); treatSens-style latent-confounder sensitivity **on a
BCF** [verified: https://pubmed.ncbi.nlm.nih.gov/27139250/] - the very class the
setData door record names as its motivation; principal-stratification BCF
[verified: https://arxiv.org/abs/2403.13256], which the critique moved here from
door 2; sequential-covariate imputation inside a BCF
[verified: https://academic.oup.com/biostatistics/article/17/3/589/1744410].

**Why it is better-shaped than door 1:** fixed n, no new data object, no
calibration change, and the work is extending revalidation across forests rather
than inventing a shape.

**`forceUpdate = TRUE` is not a workaround.** Forced collapses emptied leaves
into their parents; transactional rolls the whole change back; the per-observation
session takes no force argument at all
([[src/R_interface_bartcore.cpp:3823-3829@60f899b4]],
[[src/R_interface_bartcore.cpp:3860-3868@60f899b4]]). These are **different
posteriors, and neither is the other** - the
partial-rollback session is itself a constraint-vetoed proposal, so "force is the
deviant one" overstates the case in the memo's direction. What is unambiguous is
that force destroys the per-observation partial-rollback semantics bairrtt
depends on.

**UNPRICED - name it before designing.** Lifting D1 widens the acceptance
criterion from "no leaf empties in any tree of any chain" to "no leaf empties in
any tree of any **forest** of any chain". A K-forest multinomial or a 2-forest
BCF will therefore reject strictly more transactions than the single-forest
surface does, and the per-observation session's partial-install rate - the exact
property bairrtt depends on - falls, by an amount nobody has measured. That is a
modeling consequence, not plumbing. A design pass must measure the per-forest
veto rate at realistic K on the models that want it, and decide whether the veto
is per-forest or per-sampler, **before** committing to the semantics. This is the
one open question that could change D1 from cheap to expensive.

**LANDED, 2026-08-10 through 2026-08-12** (docs/plans/multiforest-predictor-
mutation.md, SL 7299b8b through S4). The named question resolved to
**per-sampler**: a row installs only if it empties no leaf in any tree of any
forest (and, from that plan's S3, the variance forest) of any chain, and the
per-forest column mask is the exact opt-out - a forest whose trees never split
on the touched column cannot veto, structurally. The veto-rate falsifier this
section called for ran (`docs/plans/multiforest-veto-rate-falsifier.md`,
verdict YELLOW on both column types) and its obligations were discharged in the
landing plan rather than blocking it. Multinomial is included, named and
UNGATED: its veto rate is recorded as measured stress, not covered by a
verdict. The heteroscedastic (variance-forest) class turned out to carry the
largest share of the measured cost, not a side case, and landed alongside the
mean-forest classes; doing so required also closing a real gap in
`Chain::stateIsValid`'s variance branch, which had not imposed occupancy the
way its mean branch always had, so `setState` had been admitting variance
states the widened veto would refuse. `refuseMultiForestTransactionalUpdate`
and `refuseVarianceForestPredictorMutation` are retired from all four
transactional entries.

### D2. Multinomial response-side mutation - a counts/offset channel

`MultinomialForestCombiner` does not override `supportsResponseMutation`
([[src/bartcore/combiner.hpp:411@60f899b4]] default false), so `setResponse`,
`setOffset` and `setWeights` are all refused on a shipped public sampler. The
gap is real: multinomial BART cannot participate in an outer Gibbs sampler on
its **response** side.

**Re-scoped, because the memo's premise was wrong.** Flipping the flag would
enable **nothing**. `MultinomialForestCombiner::formForestResponse`
([[src/bartcore/combiner.hpp:785-796@60f899b4]]) takes `const double* /*y*/,
const double* /*w*/` - both ignored, as its own doc comment says ("the passed
chain y is ignored"). The response *is* the combiner-owned `counts_` matrix,
fixed at construction; `setResponse` on a multinomial would write a `y`
nothing reads. D2's real content is **a new counts (and n x K offset) mutation
channel on the combiner**, plus surface work at the creation end
([[R/bart.R:829-856@60f899b4]] refuses `weights` and `offset` for
`family = "multinomial"`). The price is higher than "audit the combiner".

**Two further corrections to the memo's framing.** (a) "A multinomial BART cannot
be a conditional inside a larger Gibbs sampler at all" is FALSE: the forced
whole-matrix `setPredictor` is open for multinomial today
(`bartcore_setPredictor` gates only `forceUpdate != TRUE`,
[[src/R_interface_bartcore.cpp:3703-3706@60f899b4]]), so the *predictor*
channel - the latent-covariate shape of D1 - already works there. It is the
*response* side that is closed. (b) The shipped entry is
`bart2(family = "multinomial")`; `bartMultinomial` is the fit **class** name
([[R/generics.R:772@60f899b4]]), not an exported constructor.

**Models blocked:** competing-risks discrete-time hazard, which
docs/design/survival.md section 6 already records as "the same expansion with a
multinomial outcome"; nominal-response IRT and discrete-choice models with
alternative-specific offsets (the recorded `n x K multinomial offset` door in its
motivated form); SUR-style multi-outcome composition, which
docs/design/correlated-outcomes.md shipped for gaussian outcomes in stan4bart and
which stops at the multinomial boundary.

### D3. Public BCF creation surface - CLOSED (2026-08-10 to 2026-08-11)

**CLOSED by bcf-public-surface.** BCF is now reachable from
`dbarts(treatment = z, ...)`/`dbartsSpec(...)` (S1, 60d7eb7c): the one flat
creation entry, `dbarts_sampler_create`
([[src/C_interface.cpp:107@60d7eb7c]]), routes through
`bartcore_bridge::createHolder`, which now dispatches to BCF creation itself
when the data carries a treatment vector
([[src/R_interface_bartcore.cpp:2600-2675@60d7eb7c]]) - `createBCFHolder` is
no longer the only path in, and `dbarts:::bartcoreBCFSampler`
([[R/bartcore.R:644@60d7eb7c]], corrected from this section's stale line 629) is
no longer the only R entry. `bcf()` stays a comment
([[R/model.R:1102@60d7eb7c]], [[R/model.R:1120@60d7eb7c]]), expected to ship
in bartCause instead (docs/plans/archive/multiforest-extension-surface.md
fork 4). The R5 surface (S2, 18657f8f), the flat C surface (S3, 0ea8af56) and
per-draw reporting (S4, f2d82bfe) followed; `treatment =`/`moderators =`/
`treatmentForest =` are PROVISIONAL, scheduled for replacement by
`forests = list(forest(basis = ...))` (same plan, M2).

**This conditioned every BCF verdict above; it no longer does** - door 1,
door 2, and D1 all now have a public creation route to be consumed through.

### D4. The general per-forest multiplier (basis/amplitude) channel - CLOSED (2026-08-13 to 2026-08-14)

The channel: each forest carries its own n x q_f row-major basis `B_f` and its
own length-q_f amplitude vector `a_f`, and enters the mean scaled by the
contraction of the two,

    E[y_i] = sum_f m_{f,i} f_f(x_i),   m_{f,i} = dot(a_f, B_f(i, .))

so a MODEL is a choice of basis columns per forest. bcf is the K = 2 instance
(an intercept basis with free amplitude `a`; the two-column indicator basis
`(1 - z, z)` with amplitudes `(b0, b1)`). The mechanism, its conditionals, its
mutation and persistence contracts and its bitwise contracts are
docs/design/multiplier-combiner.md; the arc is
docs/plans/archive/multiforest-extension-surface.md M4. **Gaussian, probit and logistic
responses** (M4.4) - the creation surface admits those three and refuses `aft`,
`ordinal` and `nbinom` by name ([[R/spec.R:504-526@e5e93f11]]), and the
K-forest chain builds the matching response model
([[src/bartcore/chain.hpp:710-725@e5e93f11]]) - so every class below is
expressible for a continuous outcome and for a binary one, where the forests
combine into the INDEX rather than the mean and sigma is pinned.

The survey already carried BCF itself (section 2, ledger below). The four
classes it did NOT carry, each a different basis, are:

- **Continuous / dose-response exposure BCF** (Woody, Carvalho, Hahn, Murray).
  bcf's shape with a CONTINUOUS exposure: forest 1 takes the single column `z`
  instead of the two-level indicator, so the treatment forest's contribution
  scales linearly in dose. **Evidence status, stated on its own because it does
  not meet this survey's standard: the only source is UNPUBLISHED** -
  arXiv:2007.09845, v1 only, no journal reference. Not merely venue-unverified;
  there is no venue. Recorded here because the phrase this section executes
  ("the four VERIFIED classes it lacks",
  docs/plans/archive/multiforest-extension-surface.md M4.5) called it verified, and the
  plan's own later corrected-citations pass
  ([[docs/plans/multiforest-extension-surface.md:488-493@4c018187]])
  downgraded it. Neither claim is left standing silently: the class is
  written, and it is written as unverified. dbarts expresses it TODAY - a
  one-column continuous basis is a legal `forest(basis = )`, and there is no
  shape predicate left to route it anywhere: `drawForestAmplitude`
  ([[src/bartcore/combiner.hpp#drawForestAmplitude]]) contracts whatever the
  basis holds, so a continuous column is drawn as the design matrix it is.
- **VCBART** (Deshpande, Bai, Balocchi, Starling, Weiss), Bayesian Analysis
  21(1):281-308, 2026. Varying coefficients: forest j carries the single column
  `X_j`, so `E[y_i] = sum_j beta_j(x_i) X_ij` with each `beta_j` its own
  ensemble. **The structure matches; the MCMC does not.** VCBART has NO
  amplitude - its leaf jumps are `N(0, tau_j^2)` with `tau_j` a FIXED
  hyperparameter, never sampled
  ([[docs/plans/multiforest-extension-surface.md:488-493@4c018187]]) -
  so dbarts expresses the model
  shape with an amplitude the original does not draw, which under this survey's
  own "a model's structure is not its MCMC" rule (section 1) is a different
  sampler, not a port. Holding a forest's amplitude fixed is one knob
  (`updateAmplitude = false`), so the exact reduction is reachable.
- **Heterogeneous mediation** (Ting, Linero), JASA 120(551):1400-1413, 2025,
  doi:10.1080/01621459.2025.2491155. Forests carrying columns `1`, the
  treatment `a`, and the OBSERVED mediator `m`, which decomposes the effect into
  direct and mediated parts. Paper verified at its published venue. Expressible
  today as a three-forest `forests = list(...)` spec; the mediator column is
  data, so nothing about it needs a mutation channel.
- **Principal stratification with BCF** (Kim, Zigler), Biometrics 81(1) article
  ujaf024, 2025. Forests carrying the column `a_i`. **PRESENT in this survey,
  but only as its other half:** section 3 door 2 records that the imputed
  strata enter the outcome ensemble as SPLITTING COVARIATES and routes that
  half to D1 (this document, section 3 door 2's "What did not survive" list
  and section 4's Models-blocked list), which was and remains the right
  routing for it. What was unrecorded is the MULTIPLIER half - the arm
  indicator entering as
  a basis column rather than as a split variable. Both halves are the same
  paper; neither subsumes the other, and the D1 entry is not duplicated here.

**CLOSED.** The channel shipped across M4.0 (eb340f4f, pins), M4.1 (e3170da4 +
9459bef0, the multiplier generalization), M4.2 (a35ff7df, the q-variate
amplitude conditional and the one general ASIS rescale) and M4.3 (983d7f0a, the
K-length spec surface, the sole mutation route, and `data@bases`), with the
landing records at 53d64798 and c0c9e12f. It is reachable from the PUBLIC R
surface at any K through `dbarts(..., forests = list(forest(basis = ...)))`,
from the R5 sampler through `$setForestBasis`/`$getForestAmplitudes`, and from
the flat C API through `dbarts_sampler_setForestBasis` and the ragged amplitude
pair. Since M4.4 it also reaches `probit` and `logistic`. What it does NOT
reach is `aft`, `ordinal` or `nbinom` - each refused at creation naming what
it is missing - nor the test surface, nor nameable per-forest leaf
calibration, which the calibration map owns on every combining sampler.

---

## 5. The most valuable model class the survey found

**Causal forests whose covariate or moderator is a latent quantity resampled
every sweep, reached through D1.** Not any of doors 1-3.

The case, in the order the evidence lands:

- **An existing consumer, an existing shape, and a named blocked extension.**
  bairrtt already runs `irt_causal_bart` on dbarts through the R API: two BART
  samplers sharing a latent ability column theta, committed by
  `updatePredictorPerObservationJointly` (one sequential sweep committing
  observation j in every sampler or none). That surface was built for this
  consumer. It stops working the moment the outcome model becomes a causal
  forest - and theta as a *moderator of a treatment effect* is the natural next
  model.
- **Three more verified classes want the same channel** (D1's list above),
  including the one the setData door record itself names as motivating, and
  including one the memo had misfiled under door 2.
- **Nothing about the data shape has to change.** Fixed n, dense latent column,
  no new object, no calibration question. What is missing is per-forest
  revalidation.
- **`forceUpdate = TRUE` is a different posterior, not a workaround** - it
  collapses emptied leaves rather than rolling back, and the per-observation
  session has no force variant, so the partial-rollback semantics bairrtt depends
  on cannot be recovered by any current call.

The one thing that would change this assessment is D1's unpriced veto-scope
consequence (section 4). It does not change the *value*; it could change the
*price*.

**D1 landed (section 4); this class is now buildable engine-side.** A BCF
sampler accepts `$setPredictor` (whole matrix, column subset, or
per-observation) and the per-observation session across samplers
(`updatePredictorPerObservationJointly`) on the same terms as a single-forest
sampler, so theta as a treatment moderator no longer needs new engine work.
bairrtt's OTHER named blocker is untouched by
this landing: its MH filter still calls `$predict` on the BCF test surface,
which `refuseUndefinedTestFits` still refuses
(`docs/plans/multiforest-predictor- mutation.md` "Migration costs per
consumer").

---

## 6. Defect found while establishing the baseline - pointer only
Fixed at HEAD (not a live defect): `applyNewData` resizes and re-anchors the variance-forest storage via `resizeVarianceStorage`, and both `forceRefreshTrees` and the donor-state restore path route the new data through `refreshVarianceForest` as well.
The survey found `setData` accepted on a heteroscedastic sampler:
`refuseMultiForestMutation` keys on `numForests >= 2`
([[src/R_interface_bartcore.cpp:2649-2654@60f899b4]]), but the variance forest
is the separate `varianceForest_` member, so such a sampler reports
`numForests == 1` and passes.

The critique escalated it on two axes, and the escalation is confirmed here:

- **Wider than `setData`.** No data-mutation path in the engine re-routes the
  variance forest. `revalidateTrees` ([[src/bartcore/chain.hpp:1484@60f899b4]])
  and `applyNewData` ([[src/bartcore/chain.hpp:1557@60f899b4]]) open on
  `forests_[0]`; `forceRefreshTrees` ([[src/bartcore/chain.hpp:1628@60f899b4]])
  and `dropStaleMissingDirections`
  ([[src/bartcore/chain.hpp:1615@60f899b4]]) iterate `forests_` and never
  touch `varianceForest_`. The only site that repartitions variance trees is
  `rebuildVarianceForest`, reached from state restore, never from mutation. So
  the corruption reaches the **supported, documented** forced `setPredictor`
  idiom and `setCutPoints`, not just `setData`.
- **Memory safety, not silent wrongness.** `meanWeights_` is `assign(n, 0.0)`
  once at `buildVarianceForest` ([[src/bartcore/chain.hpp:2898@60f899b4]]) and
  never resized, while `formMeanWeights` and `sweepVarianceForest` write over
  the *new* n. The n-growing case is an out-of-bounds **write** and segfaults;
  it is reachable from the public R surface (`dbarts(..., variance = ~1)` then
  `sampler$setData(...)`). Not reachable from the flat C API - there is no
  `dbarts_sampler_setData` in inst/include/dbarts/dbarts.h.

**This is recorded, not specified, here.** The fix - extend re-routing to
`varianceForest_` across all four mutation paths, versus refuse the whole
predictor-mutation surface on `hasVarianceForest` (a user-visible removal of
currently-accepted behaviour) - is its own arc, and it is a release blocker that
must not queue behind any door in this document. The narrow `setData`-only guard
the memo recommended does **not** close it. TODO entry:
`variance-forest-mutation-routing`.

---

## 7. Backlog value scan

Every TODO entry carrying `decision-gated`, `consumer-gated`, `no demand`,
`demand-driven`/`low-ranked`, or a DOCUMENTED REFUSAL, re-read under the standing
enabling-value rule. The memo scanned 21 rows over 15 entries and claimed 12
lifts; the critique cut that to **~8 effective**, and this record follows the
critique. The proposed entry edits are in
an untracked todo-draft memo.

**Lifted (8 rows, 6 entries), each naming the model that lifts it:**
`gp-followups` (GP-leaf BART as a sequential-design surrogate, where a fixed
lengthscale asks the user to know what they are estimating); `group-by-exposure`
(clustered/multi-site causal survival - riAFTBART does this on CRAN today
[verified: https://cran.r-project.org/web/packages/riAFTBART/index.html]);
`interaction-constraints`' formal-heredity half (functional-ANOVA BART
[verified: https://arxiv.org/abs/2509.03317]); `multiforest-mutation-gaps` door 2
(mixture-cure BCF - value only, shape unresolved); `multiforest-mutation-gaps`'
n x K multinomial offset (subsumed into D2); `sparse-extensions`' rbart_vi /
linear-leaf-on-sparse half (grouped random effects on high-dimensional sparse
designs); `survival-followups` competing risks (which moves with D2) and left
truncation (registry/EHR staggered entry, which survival.md already calls the
cheapest door).

**Gates that stand (with their reason, not a shrug):** `bcf-sigma-tail-mixing`
(a considered narrowness finding plus no principled remedy);
`forest-ranef-interweaving` (measured NO-GO with a benchmark gate and a reopen
clause); `interaction-constraints`' soft path-dependent penalties (a knob on an
existing constraint, not a model); `multiforest-mutation-gaps` door 1 and
`typed-ingestion` door 3 (both now gated on a *considered failure to find an
enabled model*, which is the licensed form, rather than on absent demand);
`python-bindings` (real value, but no model an R user cannot already express);
`sparse-extensions`' perf-only remainder; `survival-followups` cloglog and
continuous-time; `typed-ingestion`'s per-observation CSC refusal - now with a
model-space reason, since every class needing per-observation mutation mutates a
dense column.

**Four memo "lifts" rejected here**, because they do not move their entry's
operative gate: `negbin-real-dispersion` and `weighted-binary` (the memo itself
concedes the remaining approximate-PG posture decision stands - the value was
never the gate); `integer-predictor-storage` (self-described as a use, not a
model - a note, not a lift); `correlated-outcomes`' AR-1 engine door (already
labelled causal-motivated, and the entry states panel/ITS/DiD compose today via a
latent AR-state plus nugget in the outer sampler; the panel-BCF citations do not
show that composition insufficient, which is what lifting the *engine* door would
require).

---

## 8. What the survey could not settle

1. **Door 2's shape.** The largest open item; the fork and its deciding evidence
   are in section 3.
2. **D1's veto-scope price** (section 4). Measurable, unmeasured.
3. **Whether one consumer justifies building door 2.** The value gate lifts on
   the standing rule; the shape decision has no second consumer to triangulate
   against.
4. **Transfer / multi-task / co-data BART.** One candidate surfaced (co-data
   learning through prior weights) and was never fetched, so nothing was assessed
   there. It is prior-side, not data-swap, so no door rides on it.
5. **Online / streaming BART: not found, and the near-neighbour named.** An
   independent broad search found no BART instance maintaining a posterior under
   data arrival without refitting. The thing a reader will ask about is
   **dynaTree**: "a dynamic tree model whose state changes in time with the
   accumulation of new data", with "particle learning algorithms that allow for
   the efficient on-line posterior filtering of tree-states" and explicit
   sequential-design hooks [verified: https://arxiv.org/abs/0912.1586;
   https://cran.r-project.org/package=dynaTree]. It is a *single* tree carried by
   a particle cloud, not an additive ensemble, hence not a counterexample - but
   it is what a reader will ask about. Two naming traps: "Sequential BART" is
   sequential *conditionals* for imputation, and `stochtree`'s
   `previous_model_json` is same-data warm start, not accrual.

**Two process lessons, recorded so the next pass does not repeat them.** (a) The
memo declared that unverified claims would never support a recommendation, then
issued its only design verdict on exactly such a claim - the door 2 mutability
premise it had itself flagged as unresolved. Verification against source resolved
it against the memo for three of four classes. A shape recommendation needs the
sampler, not the abstract. (b) Build provenance: findings "confirmed by running
the installed package" must record the build hash. The memo's probes ran against
a library eight engine commits stale, including changes to the very file the
defect lives in; the critique caught this and showed the probes transfer by
diffing the two revisions, which the memo had not done.

**Three corrections against the critique**, each re-checked here. (a) The
grouped/AFT `setData` refusals are at
[[src/R_interface_bartcore.cpp:3334-3339@60f899b4]] - the memo was right and
the critique's 3316-3321 is off by two. (b) The dense-only refusal is at
3326-3331 - memo right, critique's 3325 off by one. (The critique's third
drift claim *is* correct: `applyNewData` runs to
[[src/bartcore/chain.hpp:1611@60f899b4]], not 1608.)
(c) TBAFTcure is **not** a fork of the `bcf` package, as the critique stated; its
README says the tree code is "constructed following the framework of R package
bcf" [verified: https://github.com/roxiesun/TBAFTcure]. The correction runs in
the critique's own favour - an independent reimplementation is stronger evidence
of an unmet need than a fork - but the record should say what is true.

---

## 9. Citation ledger

Kept citations, each with its verification URL. One citation the memo used was
**dropped**: the IV/BART paper offered as the principal-stratification
"neighbour", because it has no compliance classes, no principal strata and no
latent classes, and so supported nothing that survives.

- Hahn, Murray, Carvalho (2020), Bayesian Analysis 15(3) - BCF
  [verified: https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.full]
- Dorie, Harada, Carnegie, Hill (2016), Statistics in Medicine - treatSens
  [verified: https://pubmed.ncbi.nlm.nih.gov/27139250/]
- Kim, Zigler - Bayesian nonparametric trees for principal causal effects;
  Biometrics 81(1), article ujaf024, 2025, doi:10.1093/biomtc/ujaf024 (venue
  added 2026-08-14 from the extension-surface plan's corrected-citations pass,
  [[docs/plans/multiforest-extension-surface.md:483-484@4c018187]];
  this entry was previously venue-less)
  [verified: https://arxiv.org/abs/2403.13256; code https://github.com/lit777/BPCF]
- Sun, Song, Bayesian Analysis 20(2) 345-373 - tree-based Bayesian AFT cure model
  [verified: https://projecteuclid.org/journals/bayesian-analysis/advance-publication/A-Tree-based-Bayesian-Accelerated-Failure-Time-Cure-Model-for/10.1214/23-BA1402.pdf
  - the `.full` HTML page is bot-blocked; the gold-OA PDF fetches. Body confirms
  the latent label and the uncured-subgroup latency fit "across the iterations";
  the per-iteration draw is in the Web Appendix. Implementation, which carries
  the draw explicitly: https://github.com/roxiesun/TBAFTcure]
- Murray (2021), JASA 116(534) 756-769 - log-linear BART, zero-inflated counts
  [verified: https://arxiv.org/abs/1701.01503 - the all-n zero-weight reading is
  derived from eqs. 4, 5, 12, 14, not stated verbatim by the paper]
- LongBet - panel-data BCF [verified: https://arxiv.org/abs/2406.02530]
- Bayesian causal forests for longitudinal data
  [verified: https://arxiv.org/abs/2407.11927]
- Murray, Yuan, Thall (2018), JASA 113(523) - Bayesian ML for dynamic treatment
  regimes [verified: https://pmc.ncbi.nlm.nih.gov/articles/PMC6366650/]
- Chipman, Ranjan, Wang (2012), CJS - sequential design with BART
  [verified: https://arxiv.org/abs/1203.1078]
- Xu, Daniels, Winterstein (2016), Biostatistics 17(3) - sequential BART
  imputation
  [verified: https://academic.oup.com/biostatistics/article/17/3/589/1744410]
- Xiong, Roy, Liu, Hu (2024), Contemporary Clinical Trials 142:107547 -
  covariate-adjusted Bayesian adaptive randomization on a random-intercept AFT
  BART [verified: https://pubmed.ncbi.nlm.nih.gov/38688389/]
- Deng, Chakraborty, Chen, Tan - BFTS, Thompson sampling with BART; the
  cold-start line is in the full text, not the abstract
  [verified: https://arxiv.org/abs/2602.07767; https://arxiv.org/html/2602.07767v1]
- Zhao, Liu, Lin, Chi, Davies (2024), J. Biopharmaceutical Statistics - DOD-BART,
  seamless phase I/II dose optimization
  [verified: https://pubmed.ncbi.nlm.nih.gov/39611660/]; successor DOD-PRO-BART,
  Statistics in Biopharmaceutical Research 17(3) 347-356 (2025)
  [verified: https://www.tandfonline.com/doi/full/10.1080/19466315.2025.2491596]
- riAFTBART, CRAN
  [verified: https://cran.r-project.org/web/packages/riAFTBART/index.html]
- Taddy, Gramacy, Polson (2011), JASA 106(493) 109-123 - dynamic trees
  [verified: https://arxiv.org/abs/0912.1586; https://cran.r-project.org/package=dynaTree]
- Lakshminarayanan, Roy, Teh (2015), AISTATS - Particle Gibbs BART
  [verified: https://arxiv.org/abs/1502.04622]
- He, Yalov, Hahn (2019), AISTATS - XBART
  [verified: https://arxiv.org/abs/1810.02215]
- Taddy, Chen, Yu, Wyle (2015), ICML - Bayesian and empirical Bayesian forests
  [verified: https://arxiv.org/abs/1502.02312]
- Newton, Polson, Xu - weighted Bayesian bootstrap
  [verified: https://arxiv.org/abs/1803.04559]
- ANOVA-BART - functional ANOVA decomposition for BART
  [verified: https://arxiv.org/abs/2509.03317]

**Added 2026-08-14 with section 4's D4.** These three were verified by the
extension-surface plan's corrected-citations pass
([[docs/plans/multiforest-extension-surface.md:483-493@4c018187]]),
NOT by a fresh fetch in this file. Distinguishing the two is the point of this
ledger, so it is stated rather than blurred.

- Woody, Carvalho, Hahn, Murray (arXiv:2007.09845) - continuous /
  dose-response exposure BCF. **UNPUBLISHED: v1 only, no journal reference.**
  It does NOT meet this file's verification standard, and it is the SOLE source
  for D4's continuous-multiplier class - which the extension-surface plan
  called "verified" before its own later pass downgraded it
  ([[docs/plans/multiforest-extension-surface.md:488-493@4c018187]]).
  Kept as a named class with an unverified source, rather than dropped or
  relabelled.
- Deshpande, Bai, Balocchi, Starling, Weiss - VCBART, varying-coefficient BART;
  Bayesian Analysis 21(1):281-308, 2026 (the 2024-10-04 date is the advance
  publication). Structure verified; the MCMC differs - its leaf jumps are
  N(0, tau_j^2) with tau_j a FIXED hyperparameter, never sampled, so the model
  carries no amplitude at all.
- Ting, Linero - heterogeneous mediation with BART; JASA 120(551):1400-1413,
  2025, doi:10.1080/01621459.2025.2491155.

**Verified as a negative (searched, not found):** any BART or BCF model that
maintains a single fit across an accruing sample; any online/streaming *additive
ensemble* BART.

**Not verified, and load-bearing for nothing above:** Lambert (1992) ZIP; the
Dirichlet observation-weight mechanism in Taddy et al. and in the weighted
Bayesian bootstrap; random effects in the longitudinal BCF; co-data learning for
BART.
