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

- 2026-07-08: interface-fixes-1.0 landed. All eleven code items
  (F1-F11) and all six doc items (D1-D6) done; nothing deferred.
  Taste calls (T1-T11) and the 2.0-wishlist untouched, as scoped.

  Notable calls made while implementing:
  - F2: the trees data.frame's column selection previously used a
    fixed varOrder match, which for fits with directions/missing/
    beta.* columns silently dropped them even in the working
    multi-chain case; fixed to reorder only the known columns
    present and keep any others trailing, incidentally fixing that
    too (not just the single-chain crash).
  - F3: scoped to the raw-matrix/data.frame (non-formula) branch
    per the finding; response-NA and predictor-NA under
    missing = "error" get separate, existing-wording stops; the
    default-path warning names a row count via warning() (the
    codebase's established idiom; no message() calls exist
    elsewhere in R/).
  - F4: fitted()/residuals() needed the guard in two places -
    extract.bart/extract.rbart (the shared path) and
    fitted.rbart's type = "ev" fast path, which calls the C
    routine directly on yhat.train and bypasses extract() entirely.
  - F8: scoped to the integer coercions (ntree/ndpost/nskip/
    nchain/nthread/keepevery/printevery/printcutoffs/numcut).
    Logicals (keeptrainfits, usequants, ...) share the
    wrong-internal-name failure mode but not the coercion
    mechanism - as.logical() never warns on bad input, so
    coerceOrError can't catch it there; left alone as out of scope
    for this fix.
  - F9: pre-validation only fires when 'data' is a literal
    data.frame, the case where the eventual row count is known
    without duplicating model.frame's own NSE resolution; weights
    is re-resolved early via the same findTermInFormulaData helper
    offset already uses, wrapped in tryCatch so any resolution
    mismatch falls back to the prior (leakier) behavior rather than
    introducing a new failure mode. list/environment 'data' still
    falls through to model.frame's message, unchanged.
  - F10/D1: verified both run() and predict()'s thread-count
    formals are fully inert, not merely "serial" - passing any
    value has zero effect; the sampler's real thread count is
    fixed at creation from control@n.threads. Documented the
    stronger claim.
  - F11: n.trees/n.burn only print when the sampler was kept
    (keepSampler/keepTrees = TRUE); dbarts does not otherwise
    retain them on the returned object. Kept-draws-per-chain falls
    back to varcount's dimensions when the sampler isn't kept.
  - D1: run()'s $train dimension order (n.obs x n.samples x
    n.chains, collapsing to a plain matrix at n.chains = 1)
    verified empirically on small fits before writing it up.

  Tests updated (deliberate, not silent workarounds):
  - test-bartcore.R: setTestPredictor's replacement-matrix regression
    (~line 345) built an unnamed matrix against a named x.engine,
    which used to pass expect_silent() only because F6's warning
    was dead code. Named the replacement matrix like x.engine so
    the test still isolates its actual subject (offset-sync on
    predictor replacement) rather than turning it into a naming
    test.
  - No existing test relied on F3's silent drop-then-missing="error"
    -never-fires bug; none needed updating for F3 itself.

  New regression tests added: F1 (weighted-PPD variance contrast,
  test-generics-posteriorPredictiveDistribution.R), F2 (single-chain
  rbart trees extraction, test-sampler-trees.R), F4 (keepTrainingFits
  = FALSE early stops for plot/fitted/residuals on both bart and
  rbart, plus a keepCall = FALSE print smoke test, test-plot-generics.R).

  Gates (from the worktree, R changes only, no C++ touched):
  - R CMD INSTALL .: clean.
  - tinytest::test_package("dbarts"): All ok, 2477 results (up from
    2465 at the pre-fix baseline; +12 from the new regressions).
  - benchmarks/R/equivalence.R compare
    benchmarks/baselines/equivalence-cf99a00.rds: 18/18 identical
    draws (same RNG stream); R-surface-only changes did not touch
    sampling.
  - air format --check .: clean. lintr on every touched R file:
    no lints found.
