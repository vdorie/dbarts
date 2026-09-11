# Evaluations of the not-mine decision rows

This file consolidates four neutral evaluations of every row in
docs/decisions.md marked `not mine`, plus an independent critic pass over
those evaluations. Method: an evaluator judged each row from the code and
the cited record on its own merits, independent of the record's framing;
a critic then checked the evaluator's factual claims a second time.
`evaluator` and `critic` are those two passes; `final` is the verdict that
stands - the evaluator's, unless the critic marked it OVERTURNED, in which
case the critic's; a critic finding that falls short of overturning
(WEAKENED, STANDS) leaves the evaluator's verdict as `final` and its
substance in `reason`. CLOSED means a later row in this register, or the
register's own adjudication process, already answers the question; that
row or process is named in `reason`.

Rows cover docs/decisions.md's `not mine` set as it stood when the evaluations
were made. Every row the register has gained since - dec-B74 onward, the
maintainer's own adjudication of this file's questions and of the six arcs that
followed - is marked `mine`, so the table gains no members there. One
earlier row has since turned `not mine` without gaining one: dec-B50, which
the maintainer disowned as a misattributed standing fact rather than a
decision, so it has nothing to evaluate. What the later rows changed is
recorded under the two headings below.

## All rows

| id | decision | evaluator | critic | final | reason |
|---|---|---|---|---|---|
| dec-A03 | remove startThreads/stopThreads | KEEP, add an UPGRADING line | OVERTURNED to CHANGE, citing dec-B76 | CLOSED | dec-B76 is a VD ruling dated 2026-09-08 that the thread methods stay as no-ops; the removal is reversed by that row, not by a merits argument |
| dec-A04 | drop rngKind/rngNormalKind, rename rngSeed to seed | DISCUSS a deprecated rngSeed alias | OVERTURNED to CHANGE, citing dec-B76 | CLOSED | dec-B76 rules rngSeed accepted with a warning for one release, closing the alias question the evaluator raised |
| dec-A09 | ordered factors get K-1 midpoint cuts above the cut cap | KEEP | not reviewed | KEEP | K-1 is the full candidate-split set an ordered factor has at any K; the large-K memory bound is worth revisiting only if bart later defaults to ordinal coding |
| dec-A10 | prior forest is rejection-sampled to carry no empty leaf | KEEP | not reviewed | KEEP | rejection sampling draws from the law the move kernels score; main's collapse targeted a different distribution |
| dec-A11 | residual-variance df counts positive-weight rows only | KEEP | not reviewed | KEEP | the conjugate update's df equals nu0 plus the positive-weight count; counting every row deflated the sigma posterior for fits carrying zero weights |
| dec-A13 | sampled k capped at 1e6, no warning path | KEEP | not reviewed | KEEP | the cap addresses a chain that is transient under an improper k prior; the shipped default chi(1.5, 2) keeps the cap off the ordinary path |
| dec-A14 | mutators store state only on explicit updateState = TRUE | CHANGE, re-arm the state promise on mutation | OVERTURNED to KEEP | KEEP | the two-convention hazard is documented in man/dbartsSampler-class.Rd; re-arming would discard a state an explicit storeState() had just captured under updateState = FALSE |
| dec-A15 | fitted.bart inserts ci.level as its third positional argument | KEEP | STANDS; register cost line disputed | KEEP | the 0.9-x positional call errors rather than changing meaning, matched against train/test and refused by name; register cell corrected below |
| dec-A16 | family = "auto" detects a factor response and announces | KEEP | not reviewed | KEEP | 0.9-x already dispatched on the data with no announcement and fit a 3-level factor response as gaussian on the level codes; auto replaces that with an announced route |
| dec-A17 | fit object drops NULL-valued elements | KEEP | WEAKENED | KEEP | matches base-R practice of varying the return component set; the proposed NEWS addition does not apply, since indexing a missing list element by name returns NULL rather than erroring |
| dec-A18 | NA indicator cells in the model matrix; guessNumCores returns NA | KEEP | STANDS; register cost line disputed | KEEP | main's NA path was an out-of-bounds write past the indicator array, not merely a wrong value; register cell corrected below to state the model-matrix NA case is tested |
| dec-A31 | errors are Rf_error longjmps only, no return codes | KEEP, name the leak sites | WEAKENED | KEEP | longjmp-only fits both known consumers' .Call frames; the leak claim was not reproduced in the C interface file but reproduces at the predict fan-out, which raises after building several owning containers |
| dec-A32 | a non-void return is one of three documented kinds | KEEP | not reviewed | KEEP | no entry returns 0 with two meanings; the cost is reading which of three kinds an entry is, not resolving an ambiguity |
| dec-A34 | validation is deliberately partial | KEEP | not reviewed | KEEP | matches Matrix, xts and data.table, none of which validate a consumer's arguments either; main's activeFits registry never guarded the C path |
| dec-A36 | drop the cross-repository sanitizer contract CI job | KEEP, name revdep-smoke as a mandatory pre-submission step | WEAKENED | KEEP | the named condition is already recorded as a release-procedure entry in TODO |
| dec-A39 | two threading mechanisms coexist | CHANGE, fold routeTestRows onto fanOutPredictSlabs | WEAKENED | CHANGE | routeTestRows runs once per sweep per chain, which the cached pool amortizes, unlike the once-per-call fanOutPredictSlabs; an A/B at the parallel cutoff should run before the fold, not after |
| dec-A40 | engine header-only, compiled into two translation units | DISCUSS a five-TU factory split | WEAKENED | CLOSED | settled as dec-B92: the instantiation shape is kept for the release and an independent deep dive is held for later research; the proposed parallel build did not land |
| dec-A44 | BCF treatment-forest ridge ships off | DISCUSS - unrun gate or baseline hygiene | OVERTURNED to KEEP | KEEP | the bridge states one rule - a forest travels the ridge exactly when its prior is a scale mixture - and names its own acceptance gate as unrun; the fixed-variance ridge flags are unreachable on any shipped path |
| dec-A51 | about 30 test-only accessors compile into the shipped engine | KEEP | not reviewed | KEEP | the accessors cost two vtable slots; a compile-time guard would mean the tested translation unit is not the shipped one |
| dec-A52 | a second R handle layer with its own validation ships in the package | CHANGE | self-critiqued | CHANGE | one BCF creation path: retire the test-only BCF sampler constructor and its C entry in favour of the public spec and forest() route, move the two multinomial shims to inst/common, keep the handle environment; at minimum build the BCF parameters through forestParams so the eight-slot layout exists once |
| dec-A53 | samplePriorPredictive re-derives the sigma calibration in R | KEEP | not reviewed | KEEP | the R-side derivation matches the bridge's own formula and reads the same resolved sigest; an engine accessor would need its own un-scaling |
| dec-A54 | export the composition validator and two augmentation primitives | DISCUSS - scalar sigma, augFamilies gaps | WEAKENED | CLOSED | settled in dec-A54's own VD cell (VD 2026-09-08): the helpers stay exported and the one-sweep contract is stated in the manual as the part that will not change |
| dec-A55 | eight refusal-only S3 methods plus a one-token type | CHANGE to two .default methods | OVERTURNED to KEEP | KEEP | a .default method answers every unrelated class with the same message; the per-class methods give a remedy specific to the class, which a shared default cannot |
| dec-A56 | blocks() ships alongside interactions(groups=) | KEEP | STANDS; register cost line disputed | KEEP | blocks() fixes each group's tree capacity where interactions(groups=) only forbids cross-group splits and lets allocation float; register cell corrected below |
| dec-A57 | ship Windows ARM64 NEON support before a native probe existed | KEEP, fix the stale Makevars.win comment | STANDS | KEEP | a missed ARCH hedge produces a link error, not a silent scalar fallback; the Makevars.win comment predates the native runner and should be fixed |
| dec-A60 | RNG-locked values in four files, labelled a tripwire, script-regenerated | KEEP | not reviewed | KEEP | the four files are the only per-push draw-drift detector; extending the MANIFEST's oracle-naming rule to regenerations is one process line |
| dec-A61 | equivalence only against bartcore-recorded baselines | CHANGE, adopt anchor-main.md sec 7 | OVERTURNED to KEEP | KEEP | benchmarks/R/change-balance.R and bd-balance.R are already permanent baseline-independent gates for the change move; the claim that no bartcore-descended gate can notice it reversed does not hold |
| dec-A62 | gate policy is prose: self-declared RNG class, unenforced oracle rule | CHANGE, script the oracle-naming half | WEAKENED | CHANGE | the MANIFEST header allows the oracle to be named in the commit body, so a MANIFEST-diff check would flag legitimate commits and would test for the token rather than for an oracle |
| dec-A63 | tinytest count floor of 5200 with six per-file floors | CHANGE to per-file non-emptiness plus a total floor near the measured count | STANDS on diagnosis, WEAKENED on remedy | CHANGE | 5200 sits below the 5478 static expect_ call sites; nine files can exit_file() with no skip marker, so a non-emptiness allowlist needs to name them or it unguards test-capi.R |
| dec-A64 | keep superseded baselines and the near-duplicate test file names | CHANGE the file names | OVERTURNED to KEEP | KEEP | each paired file's header already states what separates it and names its sibling; renaming would break those cross-references and the per-file floors in sanitizers.yaml |
| dec-A65 | a standing discretion grant settled user-facing forks | DISCUSS - ratify the adopted defaults, narrow the grant | WEAKENED | CLOSED | the register is itself the ratification mechanism the discussion asks for; rows dec-B74 through dec-B77 are VD rulings it has already produced |
| dec-A67 | four bare noun exports; priors bundled into one list | DISCUSS - one rule for both vocabularies | OVERTURNED to KEEP | KEEP | the prior constructors are captured unevaluated and resolved against a layered environment, which is why bundling costs them nothing; interactions and blocks are evaluated in the caller's frame before the function is entered, so export is the only way to write them inline |
| dec-A68 | documented-but-inert arguments ship | CHANGE, drop tau and align the multinomial vars formals | WEAKENED | CHANGE | man/summary.bart.Rd already documents tau as a token silently dropped when absent, so it is not undocumented dead vocabulary; the multinomial vars asymmetry is real but the cheaper fix gives summary.bartMultinomial the vars formal its siblings carry |
| dec-A69 | R5 forest index 1-based, C API 0-based | KEEP | WEAKENED | KEEP | 1-based R over 0-based C is the conversion at both named call sites; a third origin exists beyond those - sparseFactor takes 1-based positions and stores them 0-based in the exported S4 i slot |
| dec-A70 | NA response errors rather than dropping rows | DISCUSS a missing = "omit" token | STANDS, redirected | CLOSED | settled as dec-B108, which put the remedy on na.action as the evaluation said it belonged; landed, default na.keepPredictors (R/data.R) |
| dec-B11 | xbart's n.burn refuses a length past two | KEEP, name the removed phase in the message | not reviewed | KEEP | refusing beats honoring two of three elements silently; the message should name the removed per-replication burn-in phase |
| dec-B15 | Student-t nu grid floored at 3; NB dispersion integer-only | KEEP | not reviewed | KEEP | both grid full conditionals match the exact posterior under their parameterization; the floor and cap are documented in man/dbarts.Rd |
| dec-B16 | monotone leaves take shape B', the exact constrained marginal | KEEP | not reviewed | KEEP | the truncated marginal matches the constrained target over the cone posterior; the quadrature cost falls only on a leaf model users opt into |
| dec-B17 | heteroscedastic variance forest routed through the weight channel | KEEP | not reviewed | KEEP | the route is an exact chi-inverse-square leaf marginal with a mean-matched calibration; the refusals cover every family that owns the weight channel, not a silent misfit |
| dec-B18 | BCF prognostic forest takes a half-Cauchy amplitude | KEEP | not reviewed | KEEP | the half-Cauchy fixture is bcf's own prior and what its ridge needs; the disagreement with the R-side family-aware default is unreachable, since no flat-C entry creates the amplitude sampler |
| dec-B20 | ordinal: scheme A identification, Cowles-style updates, auto-dispatch | KEEP | STANDS; register cost line disputed | KEEP | main fit an ordered-factor response as as.integer(y) - 1, a continuous response on the level codes; register cell corrected below |
| dec-B25 | no 8-bit hot layer or standalone per-column widths | KEEP | not reviewed | KEEP | phase-1 measurement showed no partition win for 8-bit codes on arm64; per-column widths remain structural to the container work still ahead |
| dec-B26 | archive the support library's two thread managers, cut them | KEEP | not reviewed | KEEP | the removed managers had no callers and are preserved on a branch ref; the manager that remains (dec-A39) is the open question, not this row |
| dec-B28 | de-scale the sum of squared residuals by range squared | KEEP | not reviewed | KEEP | the fix multiplies by the response-scale factor squared, matching the residuals' own scale; main's single-range multiplier was a units error, fixed in both engines the same day |
| dec-B29 | getTrees reports NA plus directions for every categorical rule | KEEP | not reviewed | KEEP | the packed 32-bit mask 0.9-x used cannot represent 65535 levels, so one vocabulary (NA plus per-level directions) is the only choice stable across the new range |
| dec-B31 | saved states carry a format version, refuse on mismatch | DISCUSS a 0.9-x recognition message | WEAKENED | CLOSED | settled as dec-B109: a 0.9-x fit is recognized by the absent format field and refused by name at restore and predict |
| dec-B32 | state restore is semantic, not bitwise | KEEP | not reviewed | KEEP | classic and BayesTree never promised bitwise restore either; semantic restore stores less per save and matches what stan4bart's splice-and-predict path reads back |
| dec-B34 | new-data NA on a column with none in training is refused | KEEP, add a remedy clause to the message | not reviewed | KEEP | refusing beats inventing an answer with no basis in the fitted model, which is the fork the record rules on; a per-row NA option was not on that fork |
| dec-B38 | sum-to-one tolerance snaps at sqrt(machine epsilon) | KEEP | STANDS; register cost line disputed | KEEP | the tolerance is looser than main's fixed value, so nothing that validated under 0.9-x is newly refused; register cell corrected below |
| dec-B39 | accept and document the mutation path's per-update cost | KEEP | not reviewed | KEEP | the measured per-update cost is small next to a sweep; an opt-out flag would leave data@x stale for every other consumer reading it off the hot path |
| dec-B41 | document 10-25 trees for GP leaves rather than change the default | DISCUSS a message above n.trees ~50 | STANDS, strengthened | CLOSED | settled as dec-B110, a warning keyed on the realized constant-leaf fallback share, which is what the evaluation asked for; landed |
| dec-B42 | Student-t log-likelihood channel reports the t marginal | KEEP | STANDS; register cost line disputed | KEEP | the t marginal is the observation-level density WAIC and PSIS-LOO are defined on; register cell corrected below |
| dec-B43 | getSigmas/getSumsOfSquaredResiduals refuse result by name | KEEP for 1.0, schedule the formal's removal | WEAKENED | KEEP | the code comment's stated reason - the formal reserves the token against a silent future repurpose - holds on its own terms; dec-B76's one-release tombstone policy now governs the scheduled removal |
| dec-B56 | two-component version, X-macro entry list, in-header constexpr hash | KEEP, move the re-bake recipe beside the assert | not reviewed | KEEP | the X-macro half matches the xts idiom with the registration-key drift closed; the hash's value is narrow to git installs within one version, which is the case pre-1.0 |
| dec-B57 | reject the get-api dispatch table, keep per-symbol lookup | KEEP | not reviewed | KEEP | a dispatch-table struct would itself need the same layout-stability treatment dbarts_results gets, to save 48 one-time symbol lookups the per-symbol stubs already cache |
| dec-B58 | the exact-ABI hash check is opt-in, lockstep is not the default | DISCUSS - keep the flag past 1.0? | WEAKENED | CLOSED | settled as dec-B111: the flag comes off both consumers at the merge and stays off, with a CI assertion tying a hash change to a minor bump |
| dec-B59 | ABI structs lead with a size field, grow by appending pre-1.0 | KEEP | not reviewed | KEEP | size-first struct growth is the documented C-ABI pattern, applied in both read and write directions; the structSize == 0 refusal catches the mistake both known consumers are one slip from making |
| dec-B61 | item 1 covers seven entries; the basis setter is renamed, not transposed | KEEP item 1; CHANGE item 4 | STANDS | CLOSED | moot at the tip: dec-B86 took the seven entries and setForestBasis out of dbarts.h, so neither item has a header to shape |
| dec-B62 | a per-sweep callback ships inline-multi-chain only | KEEP | not reviewed | KEEP | the conditioning setters a callback exists to call are sampler-wide, so per-sweep conditioning under threaded multi-chain runs has no meaning to preserve |
| dec-B66 | multi-forest veto-rate thresholds ratified as judged, not measured | KEEP | not reviewed | KEEP | no benchmarks script, tinytest file or workflow re-runs the bands, so nothing can fail on a misjudged one; the provenance table already labels each as judged rather than measured |
| dec-B73 | bit-identity of the draw path is a toggle (never built) | CHANGE to superseded | WEAKENED | CLOSED | the shipped rule and its opposite are already cross-referenced in dec-A42's and dec-B73's own cost columns; both rows are now superseded by dec-B90, which made the reference build a configure flag rather than a runtime toggle |
| dec-C01 | register the C callables under the full symbol names | KEEP | not reviewed | KEEP | matches xts, Matrix and data.table, which all register callables under the symbol's own name; the registration key cannot drift from the symbol by construction |
| dec-C02 | bind the provider side at compile time | KEEP | not reviewed | KEEP | catches in-package ABI drift at dbarts's own build, before any consumer exists, at no registration or runtime cost |
| dec-C03 | rename every .Call registration to a bartcore-specific prefix | KEEP | not reviewed | KEEP | R_useDynamicSymbols(FALSE) means .Call by string never resolved externally, so the prefix change has no R-visible effect |
| dec-C06 | two research-only compile switches ship in the source | KEEP | not reviewed | KEEP | both switches are inert with the macro unset and are the instruments the mixing program used; re-running that program is a candidate post-1.0 action |
| dec-C07 | the stale-install detector and snapshot regenerator ship local-only | KEEP | not reviewed | KEEP | CI always builds clean, so the stale-install detector has nothing to detect there, and snapshot regeneration is meant to be a deliberate act a machine must not perform silently |

## Change without discussion

- dec-A39: fold routeTestRows onto fanOutPredictSlabs's pattern and delete the C thread manager, after an A/B measurement at the test-fit parallel cutoff.
- dec-A62: require that a commit touching a baseline file name an oracle in the MANIFEST row or the commit body; leave the RNG-class prose as prose.
- dec-A63: require every inst/tinytest/test-*.R to contribute at least one result, with a named allowlist for scripts that legitimately call exit_file(), behind a total floor near the measured count.
- dec-A68: drop "tau" from the vars defaults that carry it and give summary.bartMultinomial a vars formal matching its draws siblings (as_draws until dec-B99 replaced it). Not applied at the tip: R/diagnostics.R still carries "tau" in five defaults and summary.bartMultinomial takes no vars.
- dec-B61: moot, not owed. dbarts_sampler_setForestBasis and the seven entries left dbarts.h under dec-B86; if the multi-forest block returns for a consumer, the header goes column-major with no exception then.

- dec-A52: collapse the test-only BCF creation path onto the public spec and forest() route and move the multinomial shims out of the namespace; at minimum build the BCF parameters through forestParams once.

## Discuss with the maintainer

- dec-A40: settled 2026-09-08 under dec-B92, keep the instantiation shape. The parallel build this evaluation proposed as the latency mitigation is not in src/Makevars.in, whose sublibs loop is serial.
- dec-A54: settled 2026-09-08, keep exported; the one-sweep contract is documented as stable.
- dec-A70: settled 2026-09-08 as dec-B108, an na.action argument with a response-only default.
- dec-B31: settled 2026-09-08 as dec-B109, recognize a 0.9-x fit by the missing format field and refuse by name.
- dec-B41: settled 2026-09-08 as dec-B110, a warning keyed on the realized constant-leaf fallback share; landed with the engine-constants slice, and the warning fires through the cross-validation, multinomial, ordinal and negative-binomial doors too.
- dec-B58: settled 2026-09-08 as dec-B111, flag off after the merge with a CI assertion tying hash changes to minor bumps; the assertion landed, and both consumers have already dropped the flag.
- the four engine limits that became settings under dec-B91 (categoricalExhaustiveCap, testFitParallelCutoff, predictParallelCutoff, sparseDensityThreshold), no register row when raised: settled 2026-09-10 as dec-B116, control-only; bart() and xbart() gain control = dbarts::dbartsControl() reaching all four through the control rather than through new formals, matching the landing note's recommendation over putting categoricalExhaustiveCap on cgm().

## Register cost lines corrected

- dec-A15: the 0.9-x positional call errors rather than silently changing meaning - it is matched against train/test and refused by name.
- dec-A18: the model-matrix NA case is tested; the guessNumCores fallback is not.
- dec-A56: blocks() and interactions(groups=) express different priors over tree allocation, not one approximating the other.
- dec-B20: the 0.9-x behavior being changed is fitting an ordered-factor response as continuous on its integer level codes, not a silent change.
- dec-B38: the new tolerance is looser than main's fixed value, so nothing that validated under 0.9-x is newly refused.
- dec-B42: the t marginal is what stays comparable across residual laws; the conditional-given-augmentation alternative is the one that would not compare.
