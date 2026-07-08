# interface-review

agent: sonnet readers, one dimension each; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: findings recorded here; deliverable is a fix-before-release
        list, a 2.0-wishlist, and taste calls batched to VD

Review 2 of the retrospective program (retrospective-reviews.md).
Runs now because 1.0-0 freezes the R surface at submission: renames
and default changes cost nothing today and are breaking changes the
day after.

## Goal

Every exported entry point reads as one library rather than fifteen
years of accretion: shared concepts share names and defaults,
generics behave uniformly across fit types, failures speak the
user's language, and a newcomer can get from install to prediction
on the shipped docs alone.

## Scope

The exported surface per NAMESPACE: dbarts, bart, bart2,
dbartsControl, dbartsData, dbartsPriors, the dbartsSampler
reference class, rbart_vi, xbart, pdbart/pd2bart, predict/fitted/
extract/residuals/plot/print/summary methods, plotTree,
updatePredictorPerObservationJointly, makeind/
makeModelMatrixFromDataFrame/makeTestModelMatrix, guessNumCores.
dbarts.h is out of scope (its compat contract is already frozen);
internal BCF plumbing is out of scope (unexported, bartCause-gated).

## Method

Four read-only readers, one dimension each, findings only:

A. Argument and default coherence: the full argument matrix across
   entry points - names, casing, defaults, ordering; the same
   concept under different names (or worse, different concepts
   under one name); defaults that disagree between bart, bart2, and
   dbartsControl for the same knob; documented vs actual defaults.
B. Generics and methods: predict/fitted/extract/residuals/plot/
   print/summary across bart, rbart, and the sampler - argument
   conventions (type/value vocabularies), return shapes and
   dimension orders, combineChains coherence, behavior when the fit
   cannot support the request (keeptrees off, no test data), and
   whether the sampler's mutable surface matches its documentation.
C. Error-message and validation quality: the user-reachable stop()/
   warning()/Rf_error sites - do they fire early, name the
   argument, say what was received and what was expected; silent
   coercions and silently-ignored arguments; NA/type/length edge
   validation at the R boundary.
D. New-user walkthrough on shipped docs only (installed man pages,
   README, examples - no source reading, no docs/): regression,
   binary outcome, weights and offset, prediction on new data,
   the sampler as a conditional model, rbart_vi grouping, xbart
   cross-validation. Every stumble recorded verbatim at the moment
   it happens.

Fable adjudicates into: fix-before-release (cheap, surface-stable,
uncontroversial), 2.0-wishlist (breaking or opinionated), and taste
calls for VD (batched in chat with recommendations). Fixes become
their own items under the standard gates.

## Constraints

- Findings only; the tree does not change under this review.
- Readers are read-only, FOREGROUND only, and spawn nothing.
- Reader D judges only what ships in the installed package.

## Verification

- The three lists recorded in this file's Status; VD's decisions on
  the taste calls recorded beside them.

## Status

- 2026-07-08: started; readers A-D dispatched.
- 2026-07-08: all four readers reported; orchestrator verified every
  load-bearing claim in source or by live probe before recording.
  One reader finding REFUTED: the "continuous offset is a silent
  no-op" walkthrough claim came from a constant-offset probe, which
  the scale re-anchor absorbs exactly (correct behavior); a
  heterogeneous-offset probe shows offset fully applied (sigma
  0.49 -> 5.33). The docs are wrong in the opposite direction:
  dbarts.Rd calls offset "only useful for binary responses".

FIX-BEFORE-RELEASE, code (bugs; uncontroversial):
F1. generics.R:198 object$weigths typo: train-sample PPD ignores
    case weights silently (verified live: weighted halves both draw
    at unweighted sigma). Fix + regression test.
F2. extract(rbartFit, type = "trees") crashes for single-chain fits
    (generics.R:571-577 selects an absent chain column; the
    package's own rbart example configuration). Fix + test.
F3. bart()/dbarts() raw-matrix path silently drops incomplete rows
    BEFORE the missing = "error" check reads xHasNA (data.R
    ~556-580 vs ~703), so the safety net bart() itself requests is
    dead code. Evaluate NA state pre-filter; error under
    missing = "error", message row count otherwise.
F4. keepTrainingFits = FALSE fits: plot() dies deep in apply();
    fitted()/residuals() silently return NA vectors. Add early
    stops naming the control flag (mirror the keeptrees pattern).
F5. Numeric-at-train / factor-or-character-at-test columns are
    silently recoded against the TEST set's own alphabet
    (utility.R mapFactorColumnsToTrainingLevels skips non-factor
    training columns). Validate column types against training.
F6. validateXTest's commented-out warning branch (data.R ~102)
    falls back to positional matching with zero diagnostic;
    reinstate the warning.
F7. xbart n.test error prints the literal placeholder "[2, N]"
    (xbart.R:271).
F8. bart()'s legacy coercions bypass coerceOrError: bad ndpost/
    nskip/ntree first emit a generic NAs-introduced warning, then
    an error naming the INTERNAL slot (n.burn/n.trees) and blaming
    dbarts. Route through coerceOrError with the user-facing names.
F9. Formula-interface weights/offset length errors leak R's
    "(weights)" phrasing from model.frame; pre-validate lengths.
F10. dbartsSampler$run's thread formal is numThreads while the Rd
    and predict() say n.threads; rename to n.threads (wiring
    threading itself is wishlist; document both as currently
    serial).
F11. print.bart/print.rbart print only the deparsed call - with
    keepCall = FALSE, literally "NULL()". Print a synopsis
    (family, chains, trees, burn-in, draws) from object fields.

FIX-BEFORE-RELEASE, docs:
D1. dbartsSampler-class.Rd: no Examples and generic-style Usage
    (run(sampler) fails; the marquee embedded-sampler feature is
    undiscoverable). Add worked create/run/setResponse/run example,
    $-dispatch usage, run()'s return dimensions (currently
    undocumented and transposed vs bart's yhat convention), and
    plotTree's chainNum/sampleNum.
D2. dbarts.Rd offset: delete the "only useful for binary" claim,
    document gaussian semantics (y = f(x) + offset + eps); fix
    bart.Rd "offest" typo.
D3. rbart.Rd type doc: "value" -> "ev"; xbart.Rd example loss
    prototype arg naming (y.train vs y.test); document xbart's
    vector n.burn as distinct from the scalar n.burn elsewhere.
D4. bart.Rd "Extracting Trees" documents the rbart column order
    for the bart generic (see T2).
D5. residuals methods: document (binary fits return response-scale
    y - P(Y=1|x)); document binary fits' absent .mean fields.
D6. sigma (dbarts/xbart) vs sigest (bart family): add "same as"
    cross-references; acknowledge the intentional default
    differences (n.chains, combineChains, keepTrees, k) in Rd the
    way n.trees already is.

TASTE CALLS (VD; recommendations inline):
T1. combineChains stored-object default: bart TRUE vs bart2/
    rbart_vi FALSE, while every accessor defaults TRUE. Rec:
    document for 1.0-0, unify (FALSE) in 2.0.
T2. extract(type = "trees") column order: bart emits chain,sample;
    rbart emits sample,chain. Rec: chain-first everywhere (matches
    getTrees), fixed together with F2.
T3. rbart PPD generics never thread weights into sampleFromPPD.
    Rec: thread them (matches the documented estimand).
T4. updateState is dead in seven of ten mutators (documented as
    live). Rec: drop the parameter pre-freeze (honest signature;
    Gibbs loops explain why wiring the NA-default would cost);
    alternative: honor explicit TRUE only.
T5. rbart_vi fixes k = 2 where bart2 defaults k = NULL (binary
    gets the chi hyperprior in bart2 but never in rbart_vi, though
    rbart.Rd bundles k under "same as bart2"). Rec: decide with
    chi-default-research; doc-only for now.
T6. dbartsData's actual factors default is "categorical" while
    bart/bart2 hardcode "indicators". Rec: document now, unify 2.0.
T7. plot.pdbart cols default is the reverse of plot.bart's. Rec:
    flip pdbart to c("blue", "black").
T8. pdbart/pd2bart silently clobber sampleronly/samplerOnly. Rec:
    document (it is structurally required).
T9. Mismatched-name/same-count test columns predict by position
    with only a warning. Rec: upgrade to error pre-freeze.
T10. Vignette "Working with dbarts Saved Trees" is referenced from
    bart.Rd and DESCRIPTION but does not ship. Rec: ship or strip.
T11. dbartsControl's n.cuts rides as an attribute, not a slot
    (invisible to slotNames/validObject). Rec: document the
    exception now; slot in 2.0.

2.0-WISHLIST: retire the bart()/pdbart() legacy argument
vocabulary; align run()'s return orientation with bart's yhat
convention; casing normalization (sigdf/sigquant/sigest vs
resid.prior, rngSeed vs seed); fold rbart.priors (cauchy/gamma)
into the dbartsPriors vocabulary; echo offending values in
validation messages package-wide; standardize message register;
centralize deprecation shims; threaded prediction; announce
character-to-factor promotion.

TASTE CALLS DECIDED (VD, 2026-07-08): "Now is the time to do these
updates. It took 15 years to get to 1.0, there's no good reason to
defer." Nothing moves to 2.0.
T1. Unify combineChains NOW (stored-object default FALSE family-
    wide, accessors keep TRUE).
T2. Chain-first trees column order everywhere (extract.rbart
    changes to match bart/getTrees).
T3. Thread weights into rbart's PPD draws.
T4. Investigate the delayedAssign interaction first (an unforced
    state promise already materializes CURRENT state at first
    access, so mutate-then-first-save is covered; staleness arises
    only after the promise was forced/stored, then mutated). Then
    wire updateState in the seven mutators as an OPT-IN store
    (explicit TRUE stores; the default must not tax Gibbs loops).
    If the NA-default semantics get murky against run()'s
    NA -> control@updateState convention, ask VD before landing.
T5. Per recommendation: decide the rbart_vi k default with
    chi-default-research; doc-only for now (rbart.Rd stops claiming
    k matches bart2).
T6. bart2 defaults factors = "categorical", matching dbartsData.
    NOTE: posterior-changing for bart2 fits with factor predictors
    (model-matrix representation changes); needs the statistical
    equivalence verdict and snapshot care, not just identical-draw
    gates.
T7. Flip plot.pdbart cols to c("blue", "black").
T8. pdbart/pd2bart: stop silently clobbering samplerOnly - error
    (or warn) when the user explicitly supplies it, since the
    override is structurally required.
T9. Upgrade mismatched-name/same-count test-column positional
    matching from warning to error.
T10. The saved-trees vignette used to ship - investigate why it no
    longer does (build config or lost directory on the bartcore
    branch) and restore it.
T11. Make n.cuts a real dbartsControl slot now (pre-release, class
    change is free).
Implementation PAUSED at VD's request pending his usage
assessment; the in-flight uncontroversial-fixes implementer
completes but does not land until VD resumes.

Reader logs: walkthrough scripts under the session tmp; full
reports in the orchestration transcript. Positive findings worth
keeping: setPredictor rollback semantics verified exactly right
under adversarial probing; per-chain RNG documentation precise;
casual bart() path well documented.
