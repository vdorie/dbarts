# interface-review

agent: sonnet readers, one dimension each; fable adjudicates
rng: neutral (findings only; any fix becomes its own gated item)
budget: findings recorded here; deliverable is a fix-before-release
        list, a 2.0-wishlist, and taste calls batched to VD

Review 2 of the retrospective program (retrospective-reviews.md).
Runs now because 1.0-0 freezes the R surface at submission: renames
and default changes cost nothing today and are breaking changes the
day after.

## Status

DONE. Four read-only readers (A: argument/default coherence, B:
generics/methods, C: error-message/validation quality, D: new-user
walkthrough on shipped docs) surveyed the exported R surface ahead
of the 1.0-0 freeze and produced eleven code findings (F1-F11), six
doc findings (D1-D6), eleven taste calls (T1-T11), and a
2.0-wishlist. Every load-bearing claim was verified in source or by
live probe before being recorded; one candidate finding was refuted
this way (see D2). All seventeen fix-before-release items (F1-11,
D1-6) landed under interface-fixes-1.0, nothing deferred. VD then
decided (2026-07-08) to land all eleven taste calls now rather than
defer any to 2.0 ("It took 15 years to get to 1.0, there's no good
reason to defer") - these landed under interface-taste-1.0, after a
brief pause for VD's usage assessment. Two taste items (T6, T10)
turned out to already be resolved by prior work and needed only
verification and doc polish. The 2.0-wishlist remains deferred,
untouched. Gates stayed clean throughout: tinytest grew from a
2465-result baseline to 2477 (interface-fixes-1.0) to 2496
(interface-taste-1.0); equivalence stayed 18/18 identical draws both
times (R-surface-only changes, no sampling touched); air/lintr/
checkDocFiles clean.

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

- The three lists (fix-before-release, 2.0-wishlist, taste calls)
  recorded in Findings and Taste calls below; VD's decisions on the
  taste calls recorded beside them.

## Findings

Every finding below was checked against source or by a live probe
before being recorded; reader D's initial claim that a continuous
offset is silently a no-op did not survive this check (see D2) - it
came from a constant-offset probe, which the response-scale
re-anchor absorbs exactly (correct behavior), while a heterogeneous-
offset probe shows the offset fully applied (sigma 0.49 -> 5.33).
The documentation turned out to be wrong in the opposite direction,
which is what D2 fixes.

All eleven code findings (F1-F11) and all six doc findings (D1-D6)
below landed under interface-fixes-1.0 (2026-07-08); nothing was
deferred. Where implementing a fix required a scoping or naming
decision beyond the finding text, that decision is folded into the
item below.

### Code (F1-F11)

F1. generics.R:198 object$weigths typo: train-sample PPD ignores
    case weights silently (verified live: weighted halves both draw
    at unweighted sigma). Landed: fixed, with a new weighted-PPD
    variance-contrast regression test
    (test-generics-posteriorPredictiveDistribution.R).

F2. extract(rbartFit, type = "trees") crashes for single-chain fits
    (generics.R:571-577 selects an absent chain column; the
    package's own rbart example configuration). Landed: the trees
    data.frame's column selection previously used a fixed varOrder
    match, which for fits with directions/missing/beta.* columns
    silently dropped them even in the working multi-chain case;
    fixed to reorder only the known columns present and keep any
    others trailing, incidentally fixing that too (not just the
    single-chain crash). New regression test: single-chain rbart
    trees extraction (test-sampler-trees.R).

F3. bart()/dbarts() raw-matrix path silently drops incomplete rows
    BEFORE the missing = "error" check reads xHasNA (data.R
    ~556-580 vs ~703), so the safety net bart() itself requests is
    dead code. Landed: scoped to the raw-matrix/data.frame
    (non-formula) branch per the finding; response-NA and
    predictor-NA under missing = "error" get separate, existing-
    wording stops; the default-path warning names a row count via
    warning() (the codebase's established idiom; no message() calls
    exist elsewhere in R/). No existing test relied on the silent-
    drop-then-missing="error"-never-fires bug, so none needed
    updating for this fix itself.

F4. keepTrainingFits = FALSE fits: plot() dies deep in apply();
    fitted()/residuals() silently return NA vectors. Landed: early
    stops added naming the control flag (mirroring the keeptrees
    pattern), needed in two places - extract.bart/extract.rbart
    (the shared path) and fitted.rbart's type = "ev" fast path,
    which calls the C routine directly on yhat.train and bypasses
    extract() entirely. New regression tests: keepTrainingFits =
    FALSE early stops for plot/fitted/residuals on both bart and
    rbart, plus a keepCall = FALSE print smoke test
    (test-plot-generics.R).

F5. Numeric-at-train / factor-or-character-at-test columns are
    silently recoded against the TEST set's own alphabet
    (utility.R mapFactorColumnsToTrainingLevels skips non-factor
    training columns). Landed: column types validated against
    training, as specified.

F6. validateXTest's commented-out warning branch (data.R ~102)
    falls back to positional matching with zero diagnostic. Landed:
    warning reinstated, as specified. Test fallout: test-bartcore.R's
    setTestPredictor replacement-matrix regression (~line 345) had
    built an unnamed matrix against a named x.engine, which used to
    pass expect_silent() only because this warning was dead code;
    renamed the replacement matrix to match x.engine so the test
    still isolates its actual subject (offset-sync on predictor
    replacement) rather than turning into a naming test.

F7. xbart n.test error prints the literal placeholder "[2, N]"
    (xbart.R:271). Landed: fixed, as specified.

F8. bart()'s legacy coercions bypass coerceOrError: bad ndpost/
    nskip/ntree first emit a generic NAs-introduced warning, then
    an error naming the INTERNAL slot (n.burn/n.trees) and blaming
    dbarts. Landed: routed through coerceOrError with user-facing
    names, scoped to the integer coercions (ntree/ndpost/nskip/
    nchain/nthread/keepevery/printevery/printcutoffs/numcut).
    Logicals (keeptrainfits, usequants, ...) share the wrong-
    internal-name failure mode but not the coercion mechanism -
    as.logical() never warns on bad input, so coerceOrError can't
    catch it there; left alone as out of scope for this fix.

F9. Formula-interface weights/offset length errors leak R's
    "(weights)" phrasing from model.frame. Landed: pre-validation
    added, but it only fires when 'data' is a literal data.frame -
    the case where the eventual row count is known without
    duplicating model.frame's own NSE resolution; weights is
    re-resolved early via the same findTermInFormulaData helper
    offset already uses, wrapped in tryCatch so any resolution
    mismatch falls back to the prior (leakier) behavior rather than
    introducing a new failure mode. list/environment 'data' still
    falls through to model.frame's message, unchanged.

F10. dbartsSampler$run's thread formal is numThreads while the Rd
     and predict() say n.threads. Landed: renamed to n.threads
     (wiring threading itself is wishlist). Verified alongside D1
     that both run() and predict()'s thread-count formals are fully
     inert, not merely "serial" - passing any value has zero effect;
     the sampler's real thread count is fixed at creation from
     control@n.threads. Documented the stronger claim (see D1).
     Superseded for predict: predict()'s n.threads is now wired on
     all six generics and on the sampler's predict/predictForests,
     partitioning the saved-tree replay per (chain, draw); the flat
     C API's dbarts_sampler_predict takes the count too. run()'s
     formal is still inert (its count also sizes routeTestRows and a
     run mutates state, so the override moves more than worker
     count) and the Rd now says so of run alone.

F11. print.bart/print.rbart print only the deparsed call - with
     keepCall = FALSE, literally "NULL()". Landed: prints a
     synopsis (family, chains, trees, burn-in, draws) from object
     fields. n.trees/n.burn only print when the sampler was kept
     (keepSampler/keepTrees = TRUE), since dbarts does not otherwise
     retain them on the returned object; kept-draws-per-chain falls
     back to varcount's dimensions when the sampler isn't kept.

### Docs (D1-D6)

D1. dbartsSampler-class.Rd: no Examples and generic-style Usage
    (run(sampler) fails; the marquee embedded-sampler feature is
    undiscoverable). Landed: worked create/run/setResponse/run
    example added, plus $-dispatch usage, run()'s return dimensions,
    and plotTree's chainNum/sampleNum. run()'s $train dimension order
    (n.obs x n.samples x n.chains, collapsing to a plain matrix at
    n.chains = 1) was verified empirically on small fits before
    writing it up - it was previously undocumented and transposed
    vs bart's yhat convention. See F10 for the joint thread-formal
    verification documented in the same pass.

D2. dbarts.Rd offset: delete the "only useful for binary" claim,
    document gaussian semantics (y = f(x) + offset + eps); fix
    bart.Rd "offest" typo. Landed: as specified. This is the doc
    gap the refuted reader-D finding (above) exposed - offset is
    fully applied for gaussian responses, and the shipped docs said
    the opposite.

D3. rbart.Rd type doc: "value" -> "ev"; xbart.Rd example loss
    prototype arg naming (y.train vs y.test); document xbart's
    vector n.burn as distinct from the scalar n.burn elsewhere.
    Landed: as specified.

D4. bart.Rd "Extracting Trees" documents the rbart column order for
    the bart generic (see T2). Landed: as specified.

D5. residuals methods: document (binary fits return response-scale
    y - P(Y=1|x)); document binary fits' absent .mean fields.
    Landed: as specified.

D6. sigma (dbarts/xbart) vs sigest (bart family): add "same as"
    cross-references; acknowledge the intentional default
    differences (n.chains, combineChains, keepTrees, k) in Rd the
    way n.trees already is. Landed: as specified.

### Gates (interface-fixes-1.0)

- R CMD INSTALL .: clean.
- tinytest::test_package("dbarts"): All ok, 2477 results (up from
  2465 at the pre-fix baseline; +12 from the new regressions listed
  above).
- benchmarks/R/equivalence.R compare
  benchmarks/baselines/equivalence-cf99a00.rds: 18/18 identical
  draws (same RNG stream); R-surface-only changes did not touch
  sampling.
- air format --check .: clean. lintr on every touched R file: no
  lints found.

### Also verified as correct

Positive findings worth keeping (reader logs and full reports live
in the orchestration transcript): setPredictor rollback semantics
verified exactly right under adversarial probing; per-chain RNG
documentation precise; the casual bart() path is well documented.

## Taste calls and the 2.0-wishlist

Recorded first as recommendations batched to VD, alongside a
separate 2.0-wishlist of larger/breaking changes not tied to a
specific taste call. Implementation of the taste calls was then
paused briefly at VD's request pending his usage assessment (the
in-flight fixes implementer completed its work but did not land
until VD resumed).

VD's decision (2026-07-08): "Now is the time to do these updates.
It took 15 years to get to 1.0, there's no good reason to defer."
Nothing moved to 2.0; all eleven taste calls landed under
interface-taste-1.0 the same day. Per-item: the original finding
and recommendation, VD's decision, and what landed (including any
deviation or ambiguity resolved along the way).

T1. combineChains stored-object default: bart TRUE vs bart2/
    rbart_vi FALSE, while every accessor already defaults TRUE.
    Recommended: document for 1.0-0, unify (to FALSE) in 2.0.
    Decided: unify now rather than defer (recorded then as "default
    FALSE family-wide, accessors keep TRUE"). Landed: bart2/
    rbart_vi's stored default flipped FALSE -> TRUE, unifying the
    family on TRUE (matching bart's original default and every
    accessor); R/bart.R, R/rbart.R, man/bart.Rd, man/rbart.Rd
    updated. RNG-neutral (shapes only). Consequence: seven existing
    tests indexed uncombined shapes directly and needed
    combineChains = FALSE pinned explicitly to keep testing what
    they were actually testing, not a silent workaround -
    test-convergence-diagnostics.R, test-rbart-bartcore.R,
    test-rbart-generics.R, test-rbart-groupby.R,
    test-rbart-reproducibility.R, test-model-priors.R, and
    test-rbart-example.R. The default flip also unmasked a genuine
    latent bug (not test-only): predict.rbart's unmeasured-level
    branch (R/generics.R) chose cbind vs. array-rebuild based on the
    fit's *stored* combineChains shape instead of the already-
    locally-reshaped ranef variable's shape, which only disagreed
    once a fit's stored default differed from the value requested
    at predict() time - impossible to hit before T1 (both defaulted
    FALSE). Fixed by branching on the local ranef instead.

T2. extract(type = "trees") column order: bart emits chain,sample;
    rbart emits sample,chain. Recommended: chain-first everywhere
    (matches getTrees); fix together with F2. Decided: adopt as
    recommended. Landed: extract.rbart's varOrder (R/generics.R)
    reordered to chain-first; getTrees' own C-level emission was
    already chain-first (verified in
    src/R_interface_bartcore.cpp's emitTreeDataFrame) - rbart's
    extract was the one place still emitting sample-first. F2's
    single-chain crash fix is untouched by this change ("chain" is
    still dropped from knownOrder when absent from colnames).
    Explicit chain-first column-order assertions added to
    test-sampler-trees.R.

T3. rbart PPD generics never thread weights into sampleFromPPD.
    Recommended: thread them (matches the documented estimand).
    Decided: implement as recommended. Landed: extract.rbart/
    predict.rbart/fitted.rbart's PPD paths now thread weights.
    rbart_vi's fit object did not store weights/weights.test at all
    (checked per the task's instruction to verify) - added the same
    "needed to extract ppd" block bart2 has to
    packageRbartResults (R/rbart.R), fed by the `data` argument
    already passed in. extract.rbart/fitted.rbart (sharing a
    sample = "train"/"test" argument) now compute
    `weights <- if (sample == "train") object$weights else object$weights.test`,
    exactly like extract.bart's F1 fix. predict.rbart has no
    train/test sample of its own (it always predicts fresh newdata),
    so mirroring extract.bart's sample-keyed weights would not
    apply; ambiguity resolved by giving it a `weights` argument
    mirroring predict.bart's own instead (user-supplies weights for
    the newdata, missing -> NULL) - the closer parity target than
    the literal task wording, which assumed a sample argument
    predict.rbart does not have. man/rbart.Rd's usage block gained
    the new formal; no new \arguments entry needed (weights already
    documented in the shared bart2-bundle item). New regression
    test: rbart_vi weighted-PPD variance contrast for extract/
    fitted/predict.rbart's weights argument
    (test-generics-posteriorPredictiveDistribution.R).

T4. updateState is dead in seven of ten mutators (documented as
    live). Recommended: drop the parameter pre-freeze (honest
    signature - Gibbs loops explain why wiring the NA-default would
    cost); alternative: honor explicit TRUE only. Decided:
    investigate the delayedAssign interaction first (an unforced
    state promise already materializes CURRENT state at first
    access, so mutate-then-first-save is already covered; staleness
    only arises once the promise has been forced/stored and is then
    mutated again without a following store), then wire updateState
    in the seven mutators as an OPT-IN store (explicit TRUE stores;
    the default must not tax Gibbs loops); if the NA-default
    semantics got murky against run()'s NA -> control@updateState
    convention, ask VD before landing. Landed: the lazy-state
    analysis (delayedAssign at R/dbarts.R ~370-396) confirmed the
    read above by inspection and live check - the "state" field is
    bound with delayedAssign whose expression reads
    control@updateState and, if TRUE, calls storeState() on first
    *access* of $state, capturing whatever the sampler's state is
    at that moment; so mutate-then-first-access is always covered,
    and staleness is only possible once the promise has already
    been forced (by an earlier $state read, an explicit
    storeState(), or a run()/sampleTreesFromPrior()/
    sampleNodeParametersFromPrior() call that itself stores) and a
    later mutator call is not followed by another store. All seven
    mutators (setData, setResponse, setOffset, setWeights, setSigma,
    setPredictor, setCutPoints in R/dbarts.R) wired as OPT-IN:
    `if (identical(updateState, TRUE)) storeState()` after the
    mutation succeeds; NA (default) and FALSE store nothing,
    deliberately not following run()'s NA -> control@updateState
    convention (mutators run per-sweep in Gibbs loops; the default
    must stay free). No ambiguity surfaced that needed asking VD -
    the "opt-in only" reading was unambiguous once the delayedAssign
    behavior was confirmed. setPredictor's wrapper had to preserve
    its return value's visibility (invisible NULL vs a visible
    logical) explicitly, since capturing the result in a local and
    returning the local strips invisibility. Documented in
    man/dbartsSampler-class.Rd's shared updateState \item, split by
    method family. New test file test-sampler-updateState.R: since
    most of the seven mutators change `data`, not tree/RNG
    structure, $state's *content* is unaffected by them regardless
    of storing (verified empirically for setWeights); setCutPoints
    prunes empty leaves and so is the one used to demonstrate an
    actual content change, plus a smoke test that all seven accept
    updateState = TRUE without error.

T5. rbart_vi fixes k = 2 where bart2 defaults k = NULL (binary
    responses get the chi hyperprior in bart2 but never in
    rbart_vi, though rbart.Rd bundles k under "same as bart2").
    Recommended: decide the actual default with chi-default-
    research; doc-only for now. Decided: doc-only for now, per the
    recommendation - rbart.Rd stops claiming k matches bart2.
    Landed: man/rbart.Rd's k removed from the bart2 "same as"
    bundle; given its own \item{k} noting rbart_vi's fixed 2.0
    default vs. bart2's NULL (chi hyperprior on binary responses
    only in bart2). No default-value research beyond what already
    exists in man/bart.Rd's k documentation was needed.

T6. dbartsData's actual factors default is "categorical" while
    bart/bart2 hardcode "indicators". Recommended: document now,
    unify in 2.0. Decided: make bart2 default factors =
    "categorical" to match dbartsData now - flagged as posterior-
    changing for bart2 fits with factor predictors (model-matrix
    representation changes), needing the statistical equivalence
    verdict and snapshot care rather than just identical-draw
    gates. Landed: ALREADY RESOLVED before this task started; no
    functional change made. Tracing the redirectCall/match.call
    chain (R/bart.R, R/dbarts.R, R/data.R) showed bart2's own
    `factors` formal default (`c("categorical", "indicators")`) is
    dead unless the caller passes it explicitly - redirectCall only
    forwards args present in the caller's *actual* call, so an
    omitted `factors` falls through to dbarts()'s default, which
    falls through to dbartsData()'s, both already "categorical"
    first. Verified live: `bart2(y ~ x + f, ...)` with a 3-level
    factor `f` and no `factors` argument produces a single
    categorical column (verified against explicit
    `factors = "indicators"` giving the expected 3 dummy columns).
    Git blame traces this to b33d40b3 ("Flip the default engine to
    bartcore", 2026-07-03), which predates the interface-review
    finding this task is based on - the finding text was stale by
    the time this task was written. bart()'s own hardcoded
    `factors = "indicators"` (BayesTree-compat shim, non-
    overridable) is unaffected and correct as scoped; xbart also
    untouched. Only actual change: sharpened man/bart.Rd's
    `factors, family, missing` entry for bart2 to spell out the
    model-matrix-representation difference from bart(), since the
    doc gap (not the behavior) was real. Gate fallout (reported
    even though empty, as instructed): no tinytest failures
    traceable to bart2's factors default; equivalence's
    `categorical` scenario exercises factor predictors through
    `dbarts()` directly (samplerApi = TRUE, not bart2), so it does
    not directly probe this default, but it is included in the
    18/18 identical-draws result regardless, and no code touching
    bart2's factors resolution changed, so nothing could have
    diverged. No baseline touched.

T7. plot.pdbart cols default is the reverse of plot.bart's.
    Recommended and decided: flip pdbart to c("blue", "black").
    Landed: done - R/plot.R, man/pdbart.Rd updated. No test
    asserted the old order.

T8. pdbart/pd2bart silently clobber sampleronly/samplerOnly.
    Recommended: document (it is structurally required). Decided:
    go further than documenting - stop with an error (or warn) when
    the user explicitly supplies it, since the override is
    structurally required. Landed: pdbart.getAndInitializeSampler
    (R/partialDependence.R, shared by pdbart and pd2bart) now stops
    with "'<name>' is set internally by pdbart/pd2bart and cannot
    be overridden" if sampleronly/samplerOnly is already present in
    the redirected call before it is set - i.e. only when the
    immediate caller (or, on the bart-fit-object refit path, the
    original fit-time call) passed it explicitly; untouched when
    not passed. Regression tests added to test-pdbart.R for both
    pdbart and pd2bart.

T9. Mismatched-name/same-count test columns predict by position
    with only a warning. Recommended and decided: upgrade to error
    pre-freeze. Landed: validateXTest's (R/data.R) both-named/
    mismatched-names branch now stop()s, naming the 'x' columns
    missing from 'test' and 'test's full column set; the not-both-
    named branch's warning (reinstated by F6) is untouched.
    test-data-errors.R's expect_warning updated to expect_error
    with the new message; no other test hit this path.

T10. Vignette "Working with dbarts Saved Trees" is referenced from
     bart.Rd and DESCRIPTION but does not ship. Recommended: ship
     or strip. Decided: investigate why it no longer ships (build
     config or lost directory on the bartcore branch) and restore
     it. Landed: ALREADY RESOLVED before this task started; no
     restoration needed. `git log --all --oneline -- vignettes`
     shows the directory is not lost: vignettes/
     working_with_saved_trees.Rmd (and
     gibbs_sampler_mixture_model.Rmd) are checked in, current, and
     not in .Rbuildignore; DESCRIPTION's VignetteBuilder/Suggests
     are correctly wired. Prior commits 4e6106a ("Re-seed the
     saved-trees vignette for the new engine", 2026-07-05), d2a658b
     ("Bring the vignettes up to the 1.0-0 surface", 2026-07-07),
     and b17a78f (2026-07-07) - all ancestors of this worktree's
     base (8744e77) and all predating the interface-review's
     recorded finding (2026-07-08) - already updated the vignette's
     content and fixed whatever build-config gap once caused it to
     drop (same stale-finding pattern as T6: the review read an
     earlier tree state). Verified directly: `R CMD build .`
     succeeds and the tarball contains
     vignettes/working_with_saved_trees.Rmd,
     inst/doc/working_with_saved_trees.{Rmd,R,pdf} - the code
     chunks already execute against the installed package (knitr
     runs during build) and already reflect the current API,
     including this task's own T2 fix (the vignette's "Flattened
     Trees" section documents chain-first columns and its examples
     use chain-first subsetting). No stop-and-report-only situation
     arose since there was nothing stale to find; reported the
     investigation as instructed.

T11. dbartsControl's n.cuts rides as an attribute, not a slot
     (invisible to slotNames/validObject). Recommended: document
     the exception now; slot in 2.0. Decided: make it a real slot
     now instead - pre-release, the class change is free. Landed:
     n.cuts is now a real integer slot on dbartsControl
     (R/A_class.R), prototype 100L, with a validity check
     (`anyNA(...) || any(... <= 0L)`) deliberately not constrained
     to length 1 (n.cuts can be a per-predictor vector, recycled
     later). dbartsControl()'s constructor (R/dbarts.R) now passes
     `n.cuts = coerceOrError(n.cuts, "integer")` straight into
     newValidated() like every other slot, instead of a manual
     stop() plus attr() stash; rethrowValidityError already strips
     the S4 "invalid class ... object:" prefix, so both prior
     error messages ("must be coercible to type: integer" from
     coerceOrError, "must contain only positive integers" from the
     new validity check) are byte-identical to before - confirmed
     against test-control-errors.R without editing it. Removed the
     attr(control, "n.cuts") <- NULL strip in both read sites
     (R/dbarts.R's dbarts(), R/xbart.R) since there is no longer a
     transient attribute to clear; both now read control@n.cuts
     directly. Grepped all of R/ for `"n.cuts"` in quotes to
     confirm no other attr(., "n.cuts") reader was missed.
     Constructor signature unchanged. setControl's "settings fixed
     at creation" guard list (n.trees/n.chains/useQuantiles/
     rngSeed) deliberately NOT extended to n.cuts: it was never
     consulted after data-object construction even as an attribute
     (data@n.cuts, a real dbartsData slot, is what the engine
     actually reads), so a setControl() call with a different
     n.cuts was already a silent no-op before this change and
     remains one now - adding a guard would be a behavior change
     outside this item's scope. man/dbartsControl.Rd never
     mentioned the attribute exception, so nothing to strike there.

2.0-WISHLIST (deferred, untouched by the taste-call landing above):
retire the bart()/pdbart() legacy argument vocabulary; align run()'s
return orientation with bart's yhat convention; casing
normalization (sigdf/sigquant/sigest vs resid.prior, rngSeed vs
seed); fold rbart.priors (cauchy/gamma) into the dbartsPriors
vocabulary; echo offending values in validation messages package-
wide; standardize message register; centralize deprecation shims;
announce character-to-factor promotion.

### Gates (interface-taste-1.0)

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
  RNG-affecting change (see above), so there was no bart2-with-
  factors statistical-verdict case to report. Baseline left
  untouched.
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
