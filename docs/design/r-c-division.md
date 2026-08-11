# The R / C++ division

Status: PRINCIPLE PROPOSED, awaiting VD acceptance or edit (the text in
"The principle" below). Arc record, 2026-08-11. Artifacts (gitignored):
`.claude/r-c-division-design/{memo.md,critique.md,guide-memo.md}` - a
channel-by-family census memo, its adversarial critique (fifteen probe
scripts, run against a private build of b70b373), and a guide-first
memo written after VD twice corrected the framing ("a principle ...
should cover what we aspire to and not what we have currently
implemented"; "grounding the principle in code ... misses the entire
point of a principle which is as a guide"). The census and its
empirical claims were adversarially verified; the guide memo's two new
measurements are independently corroborated by the critique's probes,
but the guide memo as a document has NOT had its own refuting critique
- run one before treating its Part B demand census as settled fact.

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
is a price, not a wall (measured up to ~6x effective samples in the
bad case, level in the common one).

**The wall that remains is calibration, and it is measured.** A
structurally correct pure-R probit composition inherits its leaf prior
from the range of whatever vector the sampler was CONSTRUCTED on: a
16x sweep of that accident moves the implied leaf sd 16x and the
posterior sd of f by 4.6x with no error or warning; correctness is
recovered only near one lucky range. stan4bart's mvbart() refuses the
response channel outright to avoid exactly this. The composed model's
BLOCKING ports to R; its PRIOR does not, unless it can be named.

## The principle (proposed text, VD to accept or edit)

> R addresses the conditionals; C++ owns the integrand.
>
> 1. R's side. Anything expressible as "one of the engine's own models
>    fit to data the host can compute" must be reachable from R. The
>    host owns the outer loop, the residual algebra, the latent draws
>    and every non-forest block, and changes response, offset,
>    weights, predictors and sigma between sweeps (each channel where
>    the family's own identification permits it - a refusal on model
>    grounds is part of the model). Selecting an engine-provided
>    family or prior by argument is also R's side.
> 2. C++'s side. Anything that changes what the engine INTEGRATES - a
>    family, link, leaf model, tree prior, hyperprior law or
>    split-rule representation - is C++, contributed to dbarts, then
>    selected from R by argument. So is any move that must write
>    engine state, and any likelihood that does not factorise over
>    observations given the forest.
> 3. The promise. On R's side dbarts promises POSSIBLE and CORRECT,
>    and prices PRACTICAL. A composition must target the same
>    posterior as the equivalent engine model; where a difference
>    remains it is named, measured and testable - never silent.
> 4. What that obliges. R must be able to NAME the calibration it
>    composes against, rather than inherit it from a construction
>    accident. Every engine capability has an R route that reaches the
>    same posterior, even if slower. Every R affordance ships with a
>    tested recipe and a validator the author can run. Read/write
>    symmetry is a rule: whatever the engine conditions on, the host
>    may write; whatever the engine draws, the host may read.

The parenthetical in clause 1 is an orchestrator amendment absorbing
the critique's finding that sigma and per-observation weights are
refused for latent families BY IDENTIFICATION - the guide memo's
wording survives it; the census memo's "every channel for every
family" did not.

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

## Demand headlines (guide memo Part B; tool-verified, not yet
adversarially critiqued)

- 24 of 24 BART-family artifacts that ship a tree sampler wrote or
  copied their own; zero reused an engine as a library; in every case
  the new thing was an integrand change or a non-factorising
  likelihood. One author under experimental control: nbbart reuses
  dbarts while the change is data; clbart abandons every engine the
  moment the leaf stops being conjugate.
- Of four CRAN packages driving a dbarts sampler iteratively, two are
  statistically invalid (posterior means used as draws) - a 50% defect
  rate on the response channel among peer-reviewed methodology. The
  composition failure mode is real and nobody catches it today.
- The top UNMET affordance is a row subset that changes between draws
  (principal stratification, mediation, IV): princeBART monkey-patches
  eight unexported dbarts internals for it; an AOAS 2024 sampler
  builds six fresh samplers inside each of 6000 iterations. (This
  reprices the multiforest-mutation-gaps row-subsetting door, which
  was value-lifted but unscheduled on a single consumer - the demand
  survey found the class, not just a consumer.)
- Peer packages misdescribe dbarts (SoftBart claims to be the only
  embeddable BART; stochtree's feature table denies dbarts a C API);
  the capability gap is smaller than the perception gap. Recipes and
  documentation, not architecture, are the competitive answer.

## Adoption cost (new work the principle creates; committed plans
carry their own budgets)

~1950-2200 lines: mutable row subset ~450-700 (own design arc);
nameable calibration ~300 (the single highest-value item); exported
augmentation helpers ~240; a composition validator (SBC over a
user-supplied one-sweep closure) ~350; named recipes cross-referenced
from ?bart ~150; an offset-free fit accessor ~110 (four independent
codebases hand-write train - offset today); six small census items
~290; a seeding contract for composed samplers ~60 doc.

## Fork implications (docs/plans/multiforest-extension-surface.md,
Open decisions)

All four recommendations UNCHANGED; two justifications change:

1. Basis family: still (a) build post-freeze entry-gated, but
   RE-ORDERED - the "non-Gaussian multiforest is otherwise
   unreachable from R" premise is falsified by measurement, so build
   the ~300-line nameable-calibration affordance first and let FA1/FA2
   (the amplitude/ASIS mixing gates) become the whole argument.
2. Non-Gaussian v1: still (a) probit+logistic - the reason is now
   correct-by-construction calibration vs correct-by-care, not
   impossibility. AMENDMENT REQUIRED at M4 scheduling: falsifier FA5
   as committed tests a strawman (K independent probit samplers); the
   decisive arm is K GAUSSIAN samplers with host-drawn latents against
   the combined fit, which measurement predicts will AGREE.
3. Flat rename: still (a) + KEEP setForestWeights; stakes lowered (the
   flat C surface has exactly one consumer and it is in-house - 1
   reverse LinkingTo vs 23 R-level consumers).
4. bcf() home: still (a) bartCause; discoverability is answered by
   documentation, not vocabulary at the engine's front door.

Layering, sharpened: stan4bart is where the OTHER block gets a better
sampler; dbarts is where the FOREST does; everything else is a recipe,
and recipes are dbarts' documentation, not anyone's C++.

## Defects found by this arc (live at da37242; all adversarially
confirmed)

1. MEMORY-UNSAFE: setResponse on a negative-binomial sampler with a
   single negative element SEGFAULTS the R process (uncatchable;
   static_cast of a negative lround underflows to ~1.8e19 and sizes a
   histogram); y = 1e9 hangs on unbounded allocation. Same hole on the
   flat C path. Creation validates what mutation accepts.
2. SILENTLY WRONG: probit/ordinal setResponse accept out-of-support
   responses (non-0/1 y gives latents in the hundreds; constant 0.5
   collapses every latent to zero). setTreatment already shows the
   correct pattern (validates support post-creation).
3. Variance forest + setResponse/setOffset(updateScale = TRUE): the
   sigma-pin is bypassed and the fit RUNS AWAY (mean variance 156 ->
   0.70 against truth ~171 over 750 sweeps) while getSigmas() reads
   unchanged - a fifth sigma door beyond variance-forest-mutation-
   routing's four, invisible from R. A refusal keyed on
   hasVarianceForest (mirroring the BCF one) needs no algebra.
4. refuseBinaryWeightChange tells aft/ordinal/nbinom users about "a
   binary response" (~6 lines).
5. getLatents is family-polymorphic with no documentation (probit/
   ordinal/aft return a LOCATION; logistic/nbinom/Student-t a
   PRECISION; dbarts.h's "(gaussian)" parenthetical is wrong), and
   run()$train CARRIES THE OFFSET - the natural reading of
   "getLatents() - train" as a conditional draw biases slopes by 7-10
   oracle SEs; correct usage is incremental. Doc + Rd + dbarts.h fix,
   plus the offset-free fit accessor.
6. The multinomial bart2 host's $fit is a placeholder sampler that
   accepts every mutation as a silent no-op (the real K-forest sampler
   rides $bc); needs a guard or prominent documentation until the
   multinomial mutation arc replaces the area.

Recommended sequencing: items 1-4 are one stop-loss slice (the SL/VF
S1 pattern: refusals first, no draw-law change on valid paths) and
should take the next code slot ahead of multiforest S0; 5 is a docs
slice; 6 folds into whichever lands first of the stop-loss slice or
the multinomial arc.

## Errata this arc generates elsewhere

- docs/plans/multiforest-extension-surface.md needs an anchor-refresh
  pass: three stale line numbers its own header claims were verified
  (:146, :360, :970) plus at least six chain.hpp/combiner.hpp anchors
  drifted by the S4 landing (fact 12's latent draw is
  chain.hpp:1125-1126, not :1109-1112, at the current tip).
- The census memo's tier (iii).1, its fork-1/fork-2 "strengthened"
  justifications, its two-sampler latent probe, and its "nbinom
  dispersion unreachable" gap are superseded per the critique
  (recorded here; the artifact is gitignored and left as written).

## Honest limits

The decisive composition experiments are single-DGP, single-machine
measurements - decisive about mechanisms (the augmentation composes;
the calibration does not; the runaway diverges), not about magnitudes.
The demand census counts what was findable on CRAN/GitHub/arXiv. The
IACT figures carried from earlier arcs keep their recorded caveats
(none is a forest-vs-forest measurement). The guide memo's Part B has
not been adversarially re-verified; its Part A load-bearing claims
have, via the critique's independent probes.
