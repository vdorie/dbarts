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

- 2026-07-08: interface-taste-1.0 landed. All eleven taste calls
  (T1-T11) done; two (T6, T10) turned out to already be resolved by
  prior work and needed verification plus doc polish rather than a
  functional fix - reported below rather than silently no-op'd.

  Per-item status:
  - T1 done. bart2/rbart_vi's combineChains default flipped
    FALSE -> TRUE (R/bart.R, R/rbart.R); man/bart.Rd and man/rbart.Rd
    updated. RNG-neutral as scoped (shapes only).
  - T2 done. extract.rbart's varOrder (R/generics.R) reordered to
    chain-first; getTrees' own C-level emission was already
    chain-first (verified in src/R_interface_bartcore.cpp's
    emitTreeDataFrame), so this was the one place still emitting
    sample-first. Single-chain crash fix (F2) untouched - "chain" is
    still dropped from knownOrder when absent from colnames.
  - T3 done. extract.rbart/predict.rbart/fitted.rbart's PPD paths now
    thread weights. rbart_vi's fit object did not store weights/
    weights.test at all (checked, per the task's instruction to
    verify) - added the same "needed to extract ppd" block bart2 has
    to packageRbartResults (R/rbart.R), fed by the `data` argument
    already passed in. extract.rbart/fitted.rbart (which shares a
    sample = "train"/"test" argument) now compute
    `weights <- if (sample == "train") object$weights else object$weights.test`
    exactly like extract.bart's F1 fix. predict.rbart has no
    train/test sample of its own (it always predicts fresh newdata),
    so mirroring extract.bart's sample-keyed weights would not have
    applied; instead gave it a `weights` argument mirroring
    predict.bart's own (user-supplies weights for the newdata,
    missing -> NULL) - the literal task wording ("weights.test for
    test") assumed a sample argument predict.rbart does not have;
    resolved as the closer parity target (predict.bart) rather than
    invented from scratch. man/rbart.Rd's usage block gained the new
    formal; no new \arguments entry needed since "weights" is already
    in the shared bart2-bundle item.
  - T4 done. Lazy-state analysis (delayedAssign at R/dbarts.R
    ~370-396, confirmed by reading and by the live check below):
    the "state" field is bound with delayedAssign; its expression
    reads control@updateState and, if TRUE, calls storeState() on
    first *access* of $state, capturing whatever the sampler's state
    is AT THAT MOMENT. So mutate-then-first-access is always covered:
    an unforced promise fires after the mutation and captures current
    state, promise or no opt-in wiring needed. Staleness is only
    possible once the promise has already been forced (by an earlier
    $state read, an explicit storeState(), or a run()/
    sampleTreesFromPrior()/sampleNodeParametersFromPrior() call that
    itself stores) and a later mutator call is not followed by
    another store - state then keeps referring to the pre-mutation
    snapshot. Wired all seven mutators (setData, setResponse,
    setOffset, setWeights, setSigma, setPredictor, setCutPoints in
    R/dbarts.R) as OPT-IN: `if (identical(updateState, TRUE))
    storeState()` after the mutation succeeds; NA (default) and
    FALSE store nothing, deliberately NOT following run()'s
    NA -> control@updateState convention (mutators run per-sweep in
    Gibbs loops; the default must stay free). No ambiguity surfaced
    that needed asking VD - the "opt-in only" reading was unambiguous
    once the delayedAssign behavior was confirmed. setPredictor's
    wrapper had to preserve its return value's visibility (invisible
    NULL vs a visible logical) explicitly, since capturing the result
    in a local and returning the local strips invisibility. Documented
    in man/dbartsSampler-class.Rd's shared updateState \item, split by
    method family. New test file test-sampler-updateState.R: since
    most of the seven mutators change `data`, not tree/RNG structure,
    $state's *content* is unaffected by them regardless of storing
    (verified empirically for setWeights); setCutPoints prunes empty
    leaves and so is the one used to demonstrate an actual content
    change, plus a smoke test that all seven accept updateState = TRUE
    without error.
  - T5 done, doc-only. man/rbart.Rd: removed k from the bart2 "same
    as" bundle; added its own \item{k} noting rbart_vi's fixed 2.0
    default vs. bart2's NULL (chi hyperprior on binary responses only
    in bart2). No default-value research beyond what's already in
    man/bart.Rd's existing k documentation was needed to write this.
  - T6 ALREADY RESOLVED before this task started; no functional
    change made. Traced the redirectCall/match.call chain (R/bart.R,
    R/dbarts.R, R/data.R): bart2's own `factors` formal default
    (`c("categorical", "indicators")`) is dead unless the caller
    passes it explicitly - redirectCall only forwards args present in
    the caller's *actual* call, so an omitted `factors` falls through
    to dbarts()'s default, which falls through to dbartsData()'s,
    both already `"categorical"` first. Verified live: `bart2(y ~ x
    + f, ...)` with a 3-level factor `f` and no `factors` argument
    produces a single categorical column (verified against explicit
    `factors = "indicators"` giving the expected 3 dummy columns).
    Git blame traces this to b33d40b3 ("Flip the default engine to
    bartcore", 2026-07-03), which predates the interface-review
    finding this task is based on - the review's finding text was
    stale by the time this task was written. bart()'s own hardcoded
    `factors = "indicators"` (R/bart.R, BayesTree-compat shim,
    non-overridable) is unaffected and correct as scoped ("xbart
    untouched" also holds - not touched). Only change: sharpened
    man/bart.Rd's `factors, family, missing` entry for bart2 to spell
    out the model-matrix-representation difference from bart(), since
    the doc gap (not the behavior) was real.
    Gate fallout (as instructed, reported even though empty): no
    tinytest failures traceable to bart2's factors default (there
    was nothing to change); equivalence gate's `categorical` scenario
    exercises factor predictors through `dbarts()` directly (
    samplerApi = TRUE, not bart2), so it does not directly probe this
    default, but it is included in the 18/18 identical-draws result
    below regardless, and no code touching bart2's factors resolution
    changed, so nothing could have diverged. No baseline touched.
  - T7 done. plot.pdbart's cols default flipped to c("blue", "black")
    (R/plot.R, man/pdbart.Rd). No test asserted the old order.
  - T8 done. pdbart.getAndInitializeSampler (R/partialDependence.R,
    shared by pdbart and pd2bart) now stops with "'<name>' is set
    internally by pdbart/pd2bart and cannot be overridden" if
    sampleronly/samplerOnly is already present in the redirected call
    before it is set - i.e. only when the immediate caller (or, on the
    bart-fit-object refit path, the original fit-time call) passed it
    explicitly. Untouched when not passed. Regression tests added to
    test-pdbart.R for both pdbart and pd2bart.
  - T9 done. validateXTest's (R/data.R) both-named/mismatched-names
    branch now stop()s, naming the 'x' columns missing from 'test'
    and 'test's full column set. The not-both-named branch's warning
    (reinstated by F6) is untouched. Updated test-data-errors.R's
    expect_warning to expect_error with the new message; no other
    test hit this path.
  - T10 ALREADY RESOLVED before this task started; no restoration
    needed. `git log --all --oneline -- vignettes` shows the
    directory is not lost: vignettes/working_with_saved_trees.Rmd
    (and gibbs_sampler_mixture_model.Rmd) are checked in, current,
    and NOT in .Rbuildignore; DESCRIPTION's VignetteBuilder/Suggests
    are correctly wired. Prior commits 4e6106a ("Re-seed the
    saved-trees vignette for the new engine", 2026-07-05), d2a658b
    ("Bring the vignettes up to the 1.0-0 surface", 2026-07-07), and
    b17a78f (2026-07-07) - all ancestors of this worktree's base
    (8744e77) and all predating the interface-review's recorded
    finding (2026-07-08) - already updated the vignette's content and
    fixed whatever build-config gap once caused it to drop (same
    stale-finding pattern as T6: the review read an earlier tree
    state). Verified directly: `R CMD build .` succeeds and the
    tarball contains vignettes/working_with_saved_trees.Rmd,
    inst/doc/working_with_saved_trees.{Rmd,R,pdf} - the code chunks
    already execute against the installed package (knitr runs during
    build) and already reflect current API, including this task's own
    T2 fix (the vignette's "Flattened Trees" section documents
    chain-first columns and its examples use chain-first subsetting).
    No stop-and-report-only situation arose since there was nothing
    stale to find; reporting the investigation as instructed.
  - T11 done. n.cuts is now a real integer slot on dbartsControl
    (R/A_class.R), prototype 100L, with a validity check
    (`anyNA(...) || any(... <= 0L)`) deliberately not constrained to
    length 1 (n.cuts can be a per-predictor vector, recycled later).
    dbartsControl()'s constructor (R/dbarts.R) now passes
    `n.cuts = coerceOrError(n.cuts, "integer")` straight into
    newValidated() like every other slot, instead of a manual
    stop() plus attr() stash; rethrowValidityError already strips the
    S4 "invalid class ... object:" prefix, so both prior error
    messages ("must be coercible to type: integer" from
    coerceOrError, "must contain only positive integers" from the new
    validity check) are byte-identical to before - confirmed against
    test-control-errors.R without editing it. Removed the attr(control,
    "n.cuts") <- NULL strip in both read sites (R/dbarts.R's dbarts(),
    R/xbart.R) since there is no longer a transient attribute to
    clear; both now read control@n.cuts directly. Grepped all of R/
    for `"n.cuts"` in quotes to confirm no other attr(., "n.cuts")
    reader was missed. Constructor signature unchanged. setControl's
    "settings fixed at creation" guard list (n.trees/n.chains/
    useQuantiles/rngSeed) deliberately NOT extended to n.cuts: it was
    never consulted after data-object construction even as an
    attribute (data@n.cuts, a real dbartsData slot, is what the
    engine actually reads), so a setControl() call with a different
    n.cuts was already a silent no-op before this change and remains
    one now - adding a guard would be a behavior change outside this
    item's scope. man/dbartsControl.Rd never mentioned the attribute
    exception, so nothing to strike there.

  Ambiguities resolved (see per-item notes above for the reasoning):
  T3's predict.rbart weights source (predict.bart parity, not a
  sample argument that does not exist); T6 and T10 being no-ops
  (verified live rather than assumed from the task text).

  Tests updated (deliberate, not silent workarounds), all because
  bart2/rbart_vi's combineChains default flip (T1) changed the shape
  of fields these tests index directly - each pinned
  combineChains = FALSE explicitly to keep asserting the *uncombined*
  shape it was actually testing, rather than rewriting the assertions
  to the new combined shape:
  - test-convergence-diagnostics.R: the bart2 fit feeding
    bartDrawsArray's reconstruction check, and the rbart_vi fit
    feeding the same for tau.
  - test-rbart-bartcore.R: the built-in-tau-prior fit whose
    dim(fit$tau)/dimnames(fit$ranef) are asserted directly.
  - test-rbart-generics.R: both rbart_vi fits that index $ranef/
    $yhat.train/$yhat.test as raw 3-d arrays (bart-type and
    fitted-type comparisons).
  - test-rbart-groupby.R: the missing-levels fit whose $ranef is
    indexed 3-d.
  - test-rbart-reproducibility.R: both fit1/fit2, whose
    per-chain-independence check indexes $yhat.train[1L,,] and
    [2L,,] directly.
  - test-model-priors.R: the DART fit whose $varprobs shape/sum
    are asserted directly.
  - test-rbart-example.R: the canonical two-chain rbart_vi example,
    whose entire point is asserting the raw per-field shapes.
  Also found and fixed a genuine latent bug this default flip
  unmasked (not a test-only issue): predict.rbart's unmeasured-level
  branch (R/generics.R, in the code building ranef for test levels
  absent from training) branched on `length(dim(object$ranef)) == 2L`
  - the fit's *stored* shape - to decide whether to cbind or
  array-rebuild the newly-drawn unmeasured-level ranef into the
  already-locally-reshaped `ranef` variable, instead of branching on
  `length(dim(ranef))` (the local variable, already reshaped earlier
  in the same call to match *this call's* combineChains argument).
  The two disagreed whenever a fit's own stored combineChains default
  differed from the value requested at predict() time - previously
  impossible to hit by accident because both defaulted the same way
  (FALSE); T1 unifying the family default to TRUE exposed the gap the
  first time a test called predict(..., combineChains = FALSE) against
  a now-combined-by-default fit with an unmeasured test level. Fixed
  by branching on the local `ranef` instead.

  New regression tests added: T3 (rbart_vi weighted-PPD variance
  contrast for extract/fitted/predict.rbart's `weights` argument,
  test-generics-posteriorPredictiveDistribution.R), T4
  (test-sampler-updateState.R, new file), T8 (sampleronly/samplerOnly
  override errors for both pdbart and pd2bart, test-pdbart.R), plus
  explicit chain-first column-order assertions added to
  test-sampler-trees.R for T2.

  Gates (from the worktree; C++ untouched - only R/, inst/tinytest/,
  man/, docs/):
  - R CMD INSTALL . (--preclean once, out of caution for the
    dbartsControl slot addition; ordinary R CMD INSTALL . confirmed
    clean afterward too): clean.
  - tinytest::test_package("dbarts"): All ok, 2496 results (up from
    2477 at the interface-fixes-1.0 baseline; +19 from the new/
    strengthened regressions listed above). Every deliberate test
    edit is listed above; no other file needed changes.
  - benchmarks/R/equivalence.R compare
    benchmarks/baselines/equivalence-cf99a00.rds: 18/18 identical
    draws (same RNG stream) - no scenario diverged; T6 required no
    RNG-affecting change (see above), so there was no
    bart2-with-factors statistical-verdict case to report. Baseline
    left untouched.
  - R CMD build .: succeeds; tarball contains vignettes/*.Rmd,
    vignettes/references.bib, and inst/doc/*.{Rmd,R,pdf} for both
    vignettes (T10 - already shipping, reconfirmed).
  - air format --check on every touched file: clean. lintr on every
    touched R file (via lint() with the project .lintr config): 0
    lints. tools::checkDocFiles(dir = "."): 0 issues (Rd usage vs.
    formals, including predict.rbart's new `weights` argument).
    Touched Rd files (man/bart.Rd, man/rbart.Rd,
    man/dbartsSampler-class.Rd, man/pdbart.Rd) contain no non-ASCII
    bytes and parse cleanly under tools::parse_Rd.
