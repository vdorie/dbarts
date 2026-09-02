# Plans doc index

Manifest of every `docs/plans/*.md` implementation plan (30 files; `README.md`
is the process/contract doc, indexed separately at the bottom, not listed as
a plan). Grouped by cluster/theme. STATUS reflects each doc's live
`Status:`/`## Status` section (or its equivalent closing Landing note) as of
8e8a63ad, 2026-08-28 -
see `docs/README.md` for how this index relates to the other navigation
surfaces. See `docs/design/INDEX.md`
for the paired design docs. 132 further plans that are LANDED/CLOSED/NO-GO
with no present-facing reader live under `docs/plans/archive/`; their rows
are kept, verbatim, in the "Archived" section at the bottom.

Columns: `file | STATUS | one-liner`.

## BCF follow-on cluster (pairs with docs/design/bcf.md)

| file | STATUS | purpose |
|---|---|---|
| bcf-cross-host.md | LANDED 3f532af2, 2026-08-26 | prerc-surface-freeze D7: a `--cross-host` compare-only flag on bcf-equivalence.R and multinomial-equivalence.R that exempts the six snapshot channels (mu, tau, glue, varcount.tau, accepted, installed - the decided list's `forestFits` is multinomial's name and `varcount` is draws-axis, kept gated) and gates the draws-axis channels under a two-tier verdict: tier 1 LOCKED (tight relative-deviation bound on continuous channels, combinatorial identical) is the gate, tier 2 `decoupled: statistical` a labeled weak fallback with its tolerated shift quantified; non-finite guard, partition assert in the compare loop, zero-RNG storeSample discrimination probe; exact-gates.yaml the live CI host; the current .rds already suffices so the RC-tip re-record is a refresh with a checklist covering the MANIFEST anchor shift. |

## Forest / multi-forest infrastructure

| file | STATUS | purpose |
|---|---|---|
| multiforest-veto-rate-falsifier.md | RUN AND REPORTED (YELLOW both column types), 2026-08-09 | Prices the acceptance-rate cost of widening the transactional/per-observation empty-leaf veto from `forests_[0]` alone to every ensemble before multiforest-predictor-mutation opens; no KILL clause fired, decides that arc's straightforward-extension fork. |

## Tau sampler cluster

| file | STATUS | purpose |
|---|---|---|
| tau-slice-review.md | REFERENCE (review) | Reviews the grouped tau slice sampler; verdict KEEP as default, exact-IG replacement recommended for the cauchy prior only. |

## SIMD / x86 cluster

| file | STATUS | purpose |
|---|---|---|
| simd-survey.md | REFERENCE (READ-ONLY survey) | arm64 SIMD candidate survey; ranks the fused residual-roll kernel and x86 dispatch gaps as safe wins, recommends leaving suffstat scalar for bitwise reproducibility. |
| x86-simd-plan.md | PARTLY LANDED (stats.h and AVX2 fixes shipped; R1b/R2/R3 open) | x86 SIMD investigation; its two fixes landed in src/misc/simd.c and stats.h (the AVX2 leaf-7 subleaf misdetection, the missing stats.h externs), remaining items open. |

## Mutation-surface cluster

| file | STATUS | purpose |
|---|---|---|
| latent-subset-mask.md | COMMITTED | Per-observation 0/1 active-row channel (`$setActiveRows`) extending between-draw row subsetting to the latent families zero weights cannot reach; v1 gaussian/Student-t/probit/ordinal, then logistic/nbinom/aft, then multinomial (global-only). |

## Data-store / predictor-storage cluster

| file | STATUS | purpose |
|---|---|---|
| column-kind-consolidation.md | PARTLY LANDED (S0, S1, S2, S3, S4a; S4c and S4b open) | Closes four defects sharing one root in the predictor store's semantic-type axis (no ordered-factor marker below the model matrix, an overloaded `numCuts`, an unchecked double-to-code cast, `FlatKind`'s duplicated cross-check) via a three-valued `ColumnKind` plus a derived `splitsBySubset` predicate, an ordered-factor midpoint cut grid (posterior-changing), engine-side ingestion validation, and a full native-integer-ingestion/typed-storage rework so a factor column's raw is its codes rather than a codes-plus-double pair; seven slices S0-S4c landing in the ruled order S0, S1, S2, S3, S4a, S4c, S4b. S0, S1, S2, S3 and S4a have landed: the three-valued `ColumnKind` with a real `categoryCounts` field, the ordered-factor midpoint grid (draw-changing, baseline `equivalence-02d41365.rds`), engine-side ingestion validation with a matching bridge-side level-count ceiling, and the int32 code channel through the bridge's validation sweep (`DBARTS_C_API_HASH` re-baked twice, now `0x37288e7c56449b34`). |

## SBC / calibration cluster

| file | STATUS | purpose |
|---|---|---|
| sbc-calibration.md | MIXED (Tier A complete; B/C outstanding) | Simulation-based-calibration harness; found the BCF glue-on sigma mixing issue and the cauchy-tau SBC-intractable tooling gap. 646-line running log, heaviest rewrite-signal doc in the corpus. |
| sbc-family-tiers.md | BUILT d094675 (2026-08-04) | Extends SBC to ordinal, nbinom (tightened k=8), robust-t, and multinomial: t + softmax ALL PASS, ordinal 9/10 (ridge mixing), nbinom identified-mu PASS with r/psi H-MIX; the raw-f_ik finding spun out to multinomial-level-centering.md. aft/hazard/hurdle excluded as ill-posed; heteroscedastic/monotone liftable. |

## Response-family / model-surface singletons

| file | STATUS | purpose |
|---|---|---|
| weighted-binary.md | ACTIVE (parked memo, post-1.0 only) | Preserves analysis for integer-weight probit and arbitrary-real-weight logistic; deliberately not implemented in 1.0-0. |

## Grouped-effects singletons

| file | STATUS | purpose |
|---|---|---|
| group-by-exposure.md | RESEARCH-OPEN | Decision memo (not yet written) on exposing grouped random effects beyond rbart_vi via a group.by argument; blocked on demand. |

## API-surface cluster

| file | STATUS | purpose |
|---|---|---|
| adoption-slate.md | SPECCED (S1-S8 pending), 2026-08-15 | Ships the r-c-division adoption slate: per-family `getLatents` semantics on both surfaces plus `$getFitsWithoutOffset` (an engine accessor through the facade, with `storeSample` refactored onto it); the nbinom dispersion as a per-draw channel; the grouped `setResponse` relaxation, gated two-way on gaussian/aft with the correctness analysis and a pre-registered SBC falsifier; exported augmentation helpers on R's own RNG stream; `dbartsValidateComposition`, SBC over a host's one-sweep step; six tested recipes plus an embedding page; a residue slice; and the arc's ONE `dbarts.h` re-bake. Three adversarial passes discharged; fork F1 settled by VD, the rest under a delegated grant, 2026-08-15; budgets in RAW ADDITIONS (a recorded arc-level deviation from multiplier-combiner.md's dense convention); zero baseline re-records expected on any slice. |
| capi-shape.md | LANDED 9df0cb50, 2026-08-26 | The pre-RC dbarts.h shape-freeze continuation of dbarts-h-freeze.md: seven refusing entries answer a fixed-property refusal in an `int` return instead of raising (a three-class return doctrine - VALUE, TRANSACTION, CAPABILITY STATUS - stated in the header), `forest` becomes the argument after the sampler on `getTrees`/`printTrees`, `setForestBasis` renames its parameter `basisRowMajor` (the header's one stated exception to column-major), `setActiveRows` raises on a non-binary mask leaving its 0 to mean the family implements no mask, the two declared-not-read predictor-source fields kept with the rationale in the header, and `USE_FC_LEN_T` no longer defined; one hash re-bake, two stan4bart migration lines, equivalence trio bitwise. |
| composition-refusals.md | LANDED 936825d7, 2026-08-25 | The pre-RC refusal slice for prerc-surface-freeze D6 + D8: grouped random effects with a variance forest - one model, two spellings, which constructed and ran with the group block drawing at a residual variance of 1 - becomes a validation error at resolveSamplerSpec with a createHolder backstop closing the four entrances that never reach it, formals kept and a door memo naming the engine/prior/oracle an adjudication would need; and an NA in test predictors on a column complete in training is refused by name at validateXTest, with a flat-entrance backstop in validateTestSource, since a split rule learns a missing direction only where the training column had one. |
| dbarts-h-freeze.md | LANDED 6446ddce, 2026-08-25 | The pre-RC dbarts.h freeze slice for prerc-surface-freeze D4 + D3: a `dbarts_family` enum replacing the stringly-typed family at create/drawLatents/workingResponse plus a `dbarts_sampler_family` accessor, `get` ADDED to the four C readers whose R twins carry it, `useLiveTrees` on `printTrees`, `const int32_t*`/`size_t printEvery` retypes, and the stub path checking major-equality + minor-floor with hash equality behind `DBARTS_REQUIRE_EXACT_ABI`; one hash re-bake, no version constant moves, four edits each in stan4bart and treatSens, both opting into `DBARTS_REQUIRE_EXACT_ABI` until the 1.0 merge. |
| nameable-calibration.md | LANDED ab3aa2fa, 2026-08-13 | Lets an R composition name the leaf-prior calibration (`prior.scale`, response units) at creation and mid-chain instead of inheriting it from the construction range; the per-chain getter is the authoritative reader of what is in force. |
| predict-surface.md | LANDED 78f334c1, 2026-08-25 | The pre-RC R-surface slice for prerc-surface-freeze D1 + D9 + D5 + D2: one `(object, newdata, type, ...)` order across the six predict methods with `offset` the single out-of-sample offset spelling and `group.by` named-only (after dots) on rbart's predict and survivalProbabilities, `forest = NULL` meaning every forest on the four R5 readers behind one new `bartcore_numForests` bridge entry and R-side stacking, fitted/predict type vocabularies equalized (`"class"` on the two categorical predicts, `"ppd"` second on negbin/hurdle fitted), the 15 bart2.Rd alias-only S3 methods given usage, predict.rbart's `value=`/`"post-mean"` shims deleted with `value` refused by name, and one saved-tree refusal wording naming `keepTrees = TRUE` across predict/extract(trees)/plotTree; four codoc-clean commits, one bartCause line to migrate, all baselines bitwise. |
| rd-records.md | LANDED 52c10e02, 2026-08-26 | Pre-RC documentation and record corrections ahead of the human review: `print.bart`/`print.rbart` aliases + usage + value, corrected `xbart` (6-D array with the loss axis) and sampler-`$predict` value claims with the amplitude-coupling refusal carve-out, the rbart `n.chains` kept-sampler code fix, the xbart fold-view missingness scope sentence, both composition-matrix harness bugs closed (multinom dbarts5 base recipe, logistic integer-weight probe), and line-count-invariant fixed-at marks on stale records (TODO rc-gate, release-candidate-review's par claim). |
| surface-refusals.md | LANDED d48aef8a, 2026-08-26 | The pre-RC R-surface refusal-coverage slice: a name that is a formal on a sibling method but foreign to the one called is refused by name across predict/extract/fitted/residuals/survivalProbabilities via formals-derived reason tables (extract's refusal sits below the trees branch so `extract(type = "trees", newdata =)` survives); `type = "class"` with `ci.level`, `ci.level` on `type = "forest"`, forest selection outside the forest arm, and `residuals(ci.level =)` all refused; `fitted`'s positional slot 3 unified to `ci.level` and `plot.bartMultinomial` aligned to `(x, plquants, cols, ...)`; fractional counts refused inside `coerceOrError`'s integer branch (32 call sites); `$getSigmas`/`$getSumsOfSquaredResiduals` refuse their inert `result` argument. |

## Build / infra singletons

| file | STATUS | purpose |
|---|---|---|
| repo-modernization.md | MIXED (recurring/standing item) | CI/tooling hygiene bundle; sub-items 1-2 landed (widened 2026-07-22: concurrency + paths-ignore, the sanitizer PR gate DROPPED not deferred), codecov closed as a documented no. |

## Review / retrospective programs

| file | STATUS | purpose |
|---|---|---|
| architecture-numerical-review.md | REFERENCE | Two-reader fresh-eyes stress review of the frozen engine; complete, no blocking issues. |
| correctness-audit.md | REFERENCE | Re-derives every acceptance ratio / conjugate update term-by-term across 7 blocks; all CONFIRMED. |
| release-candidate-review.md | SPECCED (2026-08-17, in execution) | Pre-RC review program: census-derived slate over the two families the slice gates cannot see (baseline-rightness and accumulation); six waves plus a parallel oracle lane under a freeze protocol; six VD forks. |
| prerc-surface-freeze.md | DECIDED, 2026-08-25 (nine rulings, work items in TODO) | The pre-RC public-surface freeze: predict() signature order, keepTrees refusal text, stub version check, dbarts.h type/naming fixes, deprecation shims, composition refusals, BCF baseline format, NA-at-predict refusal, small surface items; evidence in review-2026-08-24/memos/prerc-lens*.md. |
| bartcore-review-tour.md | Current at 8e8a63ad | The merge review: what a user's script and a linked consumer see, the gates and what each proves, the risks, the open doors, and an ordered code walk. |
| pre-review-cleanup.md | LANDED 7cd71f2d, 2026-08-26 | The four adversarial pre-review reviews (staleness, completeness, YAGNI, agent accumulation) commissioned ahead of VD's own manual read, tracked verbatim under review-2026-08-24/memos/pre-review-*.md; the rulings table, the ten-commit landing, coverage deltas (a weakened tests/cpp counter, a lost dead-entrance test path), and the open doors (the seeds-axis re-record, drawShippedGlue, the still-unwired gates; the cutpoints/ordinalThresholds R-vs-C-API spelling was later decided and renamed to `thresholds`). |

## Research doors / decision-gated (no or minimal code; open per TODO)

| file | STATUS | purpose |
|---|---|---|
| gp-followups.md | RESEARCH-OPEN | Placeholder for two deferred GP-leaf extensions (sampled lengthscales, low-rank kernels); blocked on a named consumer. |
| gpu-bart.md | DONE (memo complete) | Survey memo written (docs/design/gpu-bart.md): GO on cut-scan-as-warm-start (delivered), broader GPU experiment left open behind the informed-kernel and CG-leaf prototypes (parallel-bart-frontier next-actions 4-5). |
| python-bindings.md | RESEARCH-OPEN (no spike run recorded) | Feasibility spike memo for a Python binding over bartcore; decision-gated, likely out-of-repo. |
| sparse-extensions.md | MIXED (ext (i) LANDED 2026-07-22; rest consumer-gated) | In-place nonzero mutation on sparse columns landed 343dd4c; sparse x.test, streaming range kernel, dense-backed mixed-column mutation and the smaller extensions stay deferred. |

## Process doc (not a plan)

| file | purpose |
|---|---|
| README.md | Process/contract doc: roles (Fable/Opus/Sonnet), the plan-file template, RNG gate classes and their required gates, brevity rubric, review checklist. Does not enumerate the directory's contents - that is this INDEX's job. |

## Archived (landed; records only)

Docs below are LANDED/CLOSED/NO-GO with no present-facing reader; moved to
`docs/plans/archive/` verbatim, rows kept here for the historical manifest.

| file | STATUS | purpose |
|---|---|---|
| archive/bart2-argument-consolidation.md | LANDED (S1-S14 + closure), 2026-08-16 | Consolidates bart2's argument surface: all eight forks VD-decided plus the formula-level forest term; landed semantics that recur include bart()'s subset forwarding rule, per-forest packaging gated on coupling, xbart's deliberate grid-axes-override divergence, the S13 formula-path bases rule, and the S14 split man/bart2.Rd; the arc record - do not re-derive. |
| archive/bcf-b-ridge.md | NO-GO | Treatment-scale (b) ridge interweaving move, the named suspect for the sigma residual; both controls clear it - shelved as an unimplemented future mixing win. |
| archive/bcf-bartcause-relocation.md | LANDED (D0/D2/D3/B0/B0b/B1, ARC CLOSED), 2026-08-16 | Relocates bcf-public-surface's S5 - `bcf()` and the `bartBCF` fit class - to bartCause's dbarts-1.0 branch per multiforest-extension-surface fork 4 (VD 2026-08-11), and ships the dbarts side it needs: three R-surface guards on the multi-forest construction seam (D2), the per-draw per-forest varcount channel that makes S5's contract literal (D3, the arc's one engine slice, widening `varcount`'s forest axis with the caller's declared count authoritative so every caller-owned buffer keeps today's bytes), and in bartCause a snapshot refresh, a per-chain sigma fix and the fit function itself with its five S3 methods, the `method.rsp = "bcf"` arm and the propensity-score moderator exclusion; three blind critiques discharged, eight VD decisions recorded, one shape-only bcf-equivalence re-record expected and both other baselines bitwise. |
| archive/bcf-public-surface.md | LANDED (S0-S6) | BCF reachable through public `dbarts(treatment=, moderators=, treatmentForest=)` (an ordinary `dbartsSampler`, not `bartcoreBCFSampler`), a C consumer via `dbarts_sampler_create`, per-draw mu/tau/glue reporting, and `$setTreatment` mirroring `data@treatment`; landed aa6978b (2026-08-11), argument names PROVISIONAL pending multiforest-extension-surface.md M2. |
| archive/bcf-ridge-interweaving.md | LANDED | Prognostic (a) ridge interweaving move, landed 9617c94; confirmed mixing win, sigma SBC flag persists and routes to bcf-sigma-residual. |
| archive/bcf-sigma-residual.md | RESOLVED (burn routing adopted) | Diagnoses the BCF sigma burn-in transient as slow forest-structure mixing, not a glue-scale defect; its recommended `burn = ceiling(72000/thin)` is live in benchmarks/R/sbc.R, and the extreme-tail engine remedy is the TODO door bcf-sigma-tail-mixing. |
| archive/bcf-testfits-guard.md | LANDED | Guards BCF test-fit/predict entry points that silently reported the bare prognostic forest instead of the combined a*mu+b*tau surface. |
| archive/zero-weight-exactness.md | LANDED (S0-S3, ARC COMPLETE), 2026-08-10 | A per-forest multiplier at or below R's almost-equal tolerance produces an exact zero response/weight instead of the +/-1e-9 floor-and-divide, so a row carrying no information about a forest contributes nothing to that forest's sufficient statistics; adds a caller-settable per-forest, per-observation weight for callers the glue cannot reach. |
| archive/binary-kforest-prior-default.md | LANDED (ARC COMPLETE, S0-S2), 2026-08-15 | Family-aware aPriorScale default (2 gaussian / 1 probit+logistic), the sqrt(2/K) nodeScaleFactor dispersion default (all families, K = 2 the fixed point), five calibration observability columns at three layers, the nodeScaleFactor-times-anchor product pin (first tests/cpp calibration coverage), M4.4's diffuseness-doc debt in cap-not-pin form, and the K = 1 floor fix; argued on prior coverage after arm E refuted the mixing hypothesis; zero baseline re-records expected. |
| archive/facade-shape.md | LANDED 40082c7, 2026-08-05 | Collapses SamplerBase's 21 nullary count/capability virtuals into one SamplerShape POD filled on demand by a single `shape()` virtual, so a new capability costs one field instead of declaration+forward+override; bitwise-neutral, no dbarts.h change. |
| archive/forest-combiner.md | LANDED, 2026-07-14 | Extracts BCF's hardcoded glue into a polymorphic `ForestCombiner<L>` so multinomial (and future models) can plug in without re-forking Chain. |
| archive/forest-split-bcf.md | LANDED (two phases) | Splits `Forest<L>` out of Chain and lands BCF as the first two-forest sampler (steps 1-5); a later "Phase 2 (post data-ownership-4)" wires BCF's moderators restriction - both phases complete, one file. |
| archive/multi-forest-models.md | LANDED (tracker; historical value only) | Queue/tracker for the multi-forest family (multinomial, heteroscedastic, hurdle); all three now correctly marked landed in-file. |
| archive/multiforest-extension-surface.md | LANDED (ARC COMPLETE, M0-M4.5), 2026-08-13/14 | Replaces BCF's hardcoded a*mu+b*tau glue with the general K-forest basis/amplitude family: `forests = list(forest(basis = ...))` creation (M2), `$setForestBasis` mirroring (M1), the flat mean channel folded into dbarts-h-reshape S1 (M3), and the general engine (M4, sub-sliced M4.0-M4.5) wired for gaussian/probit/logistic with aft/ordinal/nbinom refused by name; wrote docs/design/multiplier-combiner.md. |
| archive/multiforest-mutation-gaps.md | LANDED (4 commits), 2026-08-06 | Closes multinomial setOffset/setTreatment silent no-ops, the vacuously guarded BCF setWeights, and the setState/installForests leaf-scale sibling door surfaced by runsbcbcf-repair's survey and the design + critique. |
| archive/multiforest-predictor-mutation.md | LANDED (SL, S0-S4, ARC COMPLETE), 2026-08-10/11 | Retires every multi-forest transactional predictor-mutation refusal: BCF, multinomial and heteroscedastic samplers accept setPredictor/updatePredictor with rollback and the per-observation session under "no leaf empties in any tree of any ensemble", with the per-forest column mask as the opt-out; no dbarts.h change. |
| archive/multinomial.md | LANDED (ARC CLOSED, C1-C7) | Base K-forest interleaved-PG-softmax multinomial model (RNG-neutral seams + the model + bart2 surface); hub of the cluster, wrote docs/design/multinomial.md as its C6 commit. |
| archive/multinomial-counts.md | LANDED (C1-C3) | Generalizes to an n x K count matrix via PG(n_i,.) as a sum of PG(1,.) draws; single-trial reduction bitwise. |
| archive/multinomial-counts-mutation.md | LANDED (S1-S5, ARC COMPLETE), 2026-08-12 | Multinomial stops fixing its response at creation: a counts mutation channel at fixed n/K plus an n x K category offset (train, test, predict) let a softmax chain be a conditional inside a larger Gibbs/MH sampler; every channel that cannot carry the offset refuses rather than silently omitting it. |
| archive/multinomial-formula.md | LANDED (C1-C2) | Formula interface (factor LHS, cbind(...) count matrix) for bart2(family="multinomial"); RNG-neutral. |
| archive/multinomial-margins.md | LANDED | Rewrites per-sweep margin computation O(nK^2) -> O(nK); rounding-level divergence confined to multinomial channels. |
| archive/multinomial-test-surface.md | LANDED (C1-C2; C3 docs commit unconfirmed in-file) | Test-at-creation (x.test softmax probabilities) + predict-on-newdata (K-forest saved-tree replay). |
| archive/multinomial-varcounts.md | LANDED (C1; C2 docs commit unconfirmed in-file) | Per-sample per-category variable-count channel (mbart2-style). |
| archive/block-fusion-stage-a.md | SUPERSEDED (historical) | Stage A (b=1 atom-vocabulary refactor) landed and shipped as default, then the whole block-fusion project was later reverted/excised. |
| archive/block-fusion-stage-b.md | NO-GO (historical) | Stage B (b>1 fusion, the intended win); CLOSED WONT-DO, never landed, machinery later excised entirely. |
| archive/blocked-jacobi-trees.md | NO-GO (KILLED as a build target, 2026-07-21) | Noise-split single-chain tree parallelism (b trees update independently per barrier). Phase 0 exact/GO, Phase 1 promising, then a head-to-head reversed the call and a real-engine test refuted an optimistic revival attempt - strictly dominated by straight within-chain threading, also NO-GO. |
| archive/within-chain-threading.md | NO-GO (CLOSED 2026-07-13) | Within-chain parallelism for large-n single-chain workloads; step 2 prototype passed bitwise-invariance but missed both speed thresholds. Full plan file with Landing notes; the deep analysis lives in docs/design/within-chain-threading.md section 8. |
| archive/chi-default-research.md | LANDED | 48-cell simulation study picking the binary-outcome default chi(1.5, Inf) -> chi(1.5, 2); equivalence baselines re-pinned. |
| archive/chi-hyperprior-df.md | LANDED | Fixes chi(df) to sample an actual chi distribution with the stated df instead of a silently-doubled shape; draw-neutral at the default. |
| archive/chi-k-empty-leaf-count.md | LANDED | Excludes leaves stranded empty by a data mutation from the chi-k hyperprior's count/sum-of-squares accounting. |
| archive/chi-k-runaway.md | LANDED | Caps the sampled end-node precision k so an improper/weak prior scale can't let the Gibbs draw run away; bit-identical below the cap. |
| archive/tau-cauchy-exact-ig.md | LANDED | Replaces the cauchy-prior grouped tau slice sampler with an exact Makalic-Schmidt IG two-block Gibbs draw. |
| archive/tau-slice-stepout-cap.md | LANDED | Bounds sliceSampleOnce's step-out loops to prevent hangs; reproduced a 2.5-minute pre-fix hang. |
| archive/rbart-loop-profile.md | LANDED (measurement) | Measures rbart_vi custom-prior R loop overhead (12-77% by cell); surfaces the divergence bug as an out-of-scope flag. |
| archive/rbart-loop-fix.md | LANDED | Replaces per-sweep run(0,1) round trips with one per-chain callback hook; draws unchanged. |
| archive/rbart-custom-prior-divergence.md | LANDED | Fixes 15-91x tau divergence between rbart_vi's custom-prior R loop and the in-core path. |
| archive/simd-flag-multiversioning.md | NO-GO | Evaluates AVX2-flag-built variants of unrolled Mean/Variance; no material end-to-end win. |
| archive/x86-simd.md | CLOSED/SUPERSEDED | Action plan to restore/bench x86 dispatch gaps; all four goal items settled by later work (NT-store deleted, SSR SIMD leaf NO-GO, Windows-arm64 NEON landed, memory-wall-frontier.md is the live perf frontier past SIMD). |
| archive/survival-models.md | LANDED (ITEM CLOSED except grouped-surface) | AFT log-normal + discrete-time hazard native survival families. |
| archive/survival-grouped-surface.md | LANDED | Grouped (random-intercept) AFT survival reachable from R; unblocks riAFTBART migration. |
| archive/mutation-fuzzing.md | LANDED | Randomized adversarial property-test fuzzer over setPredictor/setCell/sessions/state round-trips; found two latent state-serialization edge cases. |
| archive/mutation-journal.md | LANDED | Replaces full-array copy-then-restore with build-new-and-swap; journals only changed cells past a break-even threshold. |
| archive/gate-blindspot-audit.md | LANDED | Review-3 poison sweep (16 deliberate kernel breakages) + feature x gate coverage matrix; found several blind spots. |
| archive/gate-hardening-1.0.md | LANDED | Implements six new gates targeting the audit's blind spots (bcf-exact-weak, wtgp/grouped/chik2 equivalence, linear-exact, bd-balance). |
| archive/sbc-ci-gate.md | LANDED | Weekly/manual non-blocking CI run of the gaussian SBC baseline. |
| archive/multinomial-level-centering.md | LANDED ec2a3d0 (2026-08-05) | The SBC raw-f_ik arm implicated afterCombine's level-centering precision; fixed with uniform absorption plus an exact leaf-space conditional draw, per the Memo's five blocking corrections. |
| archive/runsbcbcf-repair.md | LANDED 62caed0, 2026-08-05 | Diagnoses and repairs `runSbcBCF` (benchmarks/R/sbc.R) via FIX-B (a setModel guard precondition); acceptance PASS across thin=30/90/120. |
| archive/data-ownership.md | COMPLETE | Hub/design-tracking doc for the 5-plan program; all five landed. |
| archive/data-ownership-1-container.md | LANDED | Plan 1 of 5: engine-owned quantized (u8/u16) predictor container replacing the borrow-and-alias apparatus. |
| archive/data-ownership-2-ingestion.md | LANDED | Plan 2 of 5: data.frame-direct ingestion, the sparseFactor class, narrowing data@x's role. |
| archive/data-ownership-3-mutation.md | LANDED | Plan 3 of 5: mutation-surface rewire to call-time raw supply. |
| archive/data-ownership-4-views.md | LANDED (CLOSED) | Plan 4 of 5: per-forest column views/restriction plus a standalone handle column axis; unblocks forest-split-bcf's moderators argument. |
| archive/data-ownership-5-sparse.md | LANDED (CLOSED) | Plan 5 of 5: CSC sparse-categorical storage kernel; declares the whole 5-plan program COMPLETE. |
| archive/data-review-remediation.md | LANDED (ARC CLOSED) | Fixes 3 verified bugs + 1 first-tier perf finding (validation gaps, a quantize loop). |
| archive/data-store-consolidation.md | LANDED (ARC COMPLETE) | Consolidates the C++ data layer (train/test twinning collapsed, explicit per-column storage kind); bench-neutral. |
| archive/data-store-residuals.md | LANDED | Discharges the last 4 deferred findings (shared quantize core, isView provenance-only, conditional raw gather, linear-leaf row-major layout). |
| archive/heteroscedastic.md | LANDED (C1-C4) | Second (variance) forest for y=f(x)+s(x)eps via conjugate scaled-inv-chi-sq leaves; homoscedastic fits byte-identical throughout. |
| archive/hurdle.md | LANDED (3 commits) | family="hurdle.lognormal" composed R-side from a probit occupancy fit + a gaussian fit on log(y>0); engine byte-neutral. |
| archive/monotone-bart.md | LANDED (3 commits) | Per-variable monotone constraints via a constrained constant-leaf model; a bias in an intermediate design was caught by the exact-posterior gate before shipping. |
| archive/negative-binomial.md | LANDED (C1-C3) | family="nbinom" via PG augmentation; sweep-order and real-shape-PG issues caught pre-implementation, corrected to integer-dispersion-only. |
| archive/ordinal-outcomes.md | LANDED (C1-C3 + follow-up) | family="ordinal" via cumulative probit with Cowles-style cutpoints; auto-dispatches on ordered factors (fixed a pre-existing silent bart2 bug). |
| archive/robust-errors.md | LANDED (ARC CLOSED) | Student-t residuals via scale-mixture augmentation (resid.dist); both estimated-nu and fixed-nu modes shipped. |
| archive/variance-forest-mutation-routing.md | LANDED (S1-S5), 2026-08-08/09 | Fixes three memory-safety holes (an out-of-bounds setData heap write, an out-of-bounds cut-grid read after setCutPoints, and after setData), stops the global sigma being unpinnable behind the variance forest's back, and routes the heteroscedastic predictor-mutation surface correctly instead of accepting it and computing silently wrong; no dbarts.h change. |
| archive/weighted-binary-ppd.md | LANDED | Fixes weighted-binary posterior-predictive draws (was degenerate two-point, now coherent binomial); closes the last live path of issue #79. |
| archive/grouped-equivalence.md | CLOSED, LANDED | The grouped scenario shipped in 0e9ccca and grouped_aft in ac6ec2c; both sit in the current equivalence baseline (MANIFEST). This doc predates the landing and is kept as history. |
| archive/collapse-merge.md | LANDED | Unifies the two subtree-collapse sites to merge leaf parameters by the same effective-observation-weighted rule. |
| archive/constant-leaf-fits.md | LANDED (compare discharged 2026-08-04) | Replaces the per-tree fits slab with node-indexed mu tables + uint32 leafOf maps for constant leaves; two commits landed bitwise. The x86 bench discharge measured the intended 14-18% n=1e4 win plus a +22-28% setPredictor-accept regression -> setpredictor-leafof-rebuild.md. |
| archive/constant-leaf-suffstat.md | LANDED | Switches ConstantGaussianLeaf to a one-pass, order-insensitive sufficient statistic; a bench amendment fixed a 12-18% regression via kernel unrolling. |
| archive/hot-layer-u8.md | NO-GO (phase 2) | Measured u8 column-width kernels vs u16; no partition-throughput win on arm64, folded per-column widths into data-ownership instead. |
| archive/kernel-cleanups.md | LANDED | Exposes fast serial moment accumulators without the null-thread-manager indirection; adds an injectable misc.a output hook. |
| archive/linear-leaf-reuse.md | LANDED | Per-(tree,node) crossproduct (U'WU) cache for linear leaves; 18-35% sweep speedup measured. |
| archive/parallel-falsifiers.md | LANDED | Measured go/no-go on three parallel-BART-frontier falsifiers (E1/E2/E3); all three survived. |
| archive/test-fit-parallel.md | LANDED | Parallelizes test-fit routing over rows above a cutoff; byte-identical at any thread count, 1.48-3.38x measured. |
| archive/cheap-uniformity.md | LANDED (S0-S4, ARC COMPLETE), 2026-08-08/09 | Closes three dense-only asymmetries: a sparse-valued setPredictor/updatePredictor mutates a sparse-backed design without densifying, a single x.test column of a container-backed test set can be replaced, and a sparse test set predicts without full n x p materialization (34.5x faster, 9.96x smaller peak RSS at n=1e5/p=1e4/density 0.01); fixes three live defects along the way (a transposed `dgCMatrix` misread as a mutation argument, a two-way-wrong missingness predicate, a reference-level container/predict() mismatch). |
| archive/csc-code-validation.md | LANDED df79f17, 2026-08-05 | Bounds categorical codes against the declared level count at every ingestion entrance - training and test, dense and CSC slice, reference code. |
| archive/cutpoints-shrink-orphan.md | LANDED | Fixes setCutPoints grid-shrink leaving an ordinal split indexing past the new grid. |
| archive/family-on-model.md | LANDED | Moves the response-family slot from dbartsControl to dbartsModel. |
| archive/flat-format-v2.md | LANDED | Replaces FlatNode's bit-cast mask encoding with a tagged value/mask/pool union; state format version bumped to 2. |
| archive/fuzz-state-roundtrip.md | LANDED | Fixes two getState/setState round-trip edges the mutation fuzzer found. |
| archive/range-scaling.md | LANDED | Decision memo (keep range scaling over sd-standardization) plus an updateScale flag. |
| archive/state-continuation.md | LANDED | Drops bitwise-continuation-only state fields in favor of semantic restore. |
| archive/state-format-policy.md | LANDED | Saved states carry a format version + package provenance; incompatible loads refuse cleanly. |
| archive/test-data-parity.md | LANDED (CLOSED) | Test-side data store gains resident sparse storage (no densification); 1.83x-6.98x memory/code-shrink measured. |
| archive/typed-ingestion.md | LANDED (Slices 1, 2a, 2b-pre, 2b), 2026-08-06 | One borrowed typed PredictorSource view (dense-double \| dense-integer \| CSC; declared K) replaces SamplerOptions' 8 ingestion fields, the 7-arg setTestData, and the dense-only mutation entries; slice 0 spun out to csc-code-validation.md. No dbarts.h change. |
| archive/c-api-growth.md | LANDED | Grows the C API with a size-first results struct, additive by-name state blocks, and a per-family loglik channel. |
| archive/capi-callbacks.md | LANDED | Adds a per-sweep conditioning callback and observer hook to the C API. |
| archive/capi-dispatch-table.md | LANDED (arc complete) | ABI-compatibility mechanism for LinkingTo consumers (X-macro single-source stubs + version/hash handshake); both dbarts and stan4bart landed. |
| archive/consumer-spec-surface.md | LANDED, 2026-07-25 | Exports `dbartsSpec()` over a shared internal resolution, so an embedding package builds a sampler specification without reaching into `dbarts:::parsePriors`. |
| archive/dbarts-h-reshape.md | LANDED | Re-signs the flat C API onto a self-describing predictor-source POD (dense or CSC, self-widthed) and a forest-indexed tree-query family; one hash re-bake, no version-constant move. |
| archive/interface-review.md | LANDED | Review-2 retrospective auditing the exported R surface; 11 code fixes + 6 doc fixes + 11 taste calls, all landed same day. |
| archive/pre-release-surface-fixes.md | LANDED | Fixes aft-loglik defect + freeze-regret paper-cuts from a pre-release surface audit. |
| archive/r-ingestion-cleanups.md | LANDED | Unifies duplicated classification-family routing, sparse test-matrix handling, missing-policy guards. |
| archive/r5-cleanup.md | LANDED | Removes startThreads/stopThreads no-ops, documents offset sync, hides internal class names from S4 validity errors. |
| archive/cran-readiness.md | CLEAN (2026-07-25) | R CMD check --as-cran plus ASAN/UBSAN sweep; the 2026-07-25 pre-submission gate run closes it at Status OK with zero NOTEs (the earlier 1-NOTE ___stderrp pass is superseded). |
| archive/equivalence-ci.md | LANDED | Scheduled CI run of the equivalence gate in statistical mode with baseline provenance (a MANIFEST). |
| archive/gp-cache-test-flake.md | LANDED | Fixed a build-layout-dependent flaky GP-leaf kernel-cache test. |
| archive/post-mutation-assertions.md | LANDED | Replaces near-tautological mutation-test asserts with statistical from-scratch-agreement checks. |
| archive/rchk-ci.md | LANDED | Scheduled rchk (PROTECT-balance) CI job on the bridge; first real run's ~50 findings fixed and re-verified zero. |
| archive/revdep-smoke-ci.md | LANDED | Monthly CI smoke test of stan4bart + bartCause against dev dbarts. |
| archive/small-validation-fixes.md | LANDED | dbartsData precision-degenerate-response warning + raw-matrix incorporate NA-row fix. |
| archive/snapshot-tests.md | LANDED | Per-family reproducibility snapshot files replace whole-file-history-dependent literals. |
| archive/test-suite-trim.md | LANDED | Trims tinytest suite wall time 33.8s -> 20.3s without losing coverage. |
| archive/tests-cpp-split.md | LANDED | Splits the monolithic tests/cpp TU into per-area TUs with dependency tracking. |
| archive/autoconf-dead-code.md | LANDED | Strips dead configure knobs (--with-xint-size, Solaris probes, orphaned m4); ~3700 net lines deleted. |
| archive/bridge-error-path-leaks.md | LANDED | Wraps bridge entry points so Rf_error longjmps don't leak C++ heap; confirmed zero-leak by valgrind CI. |
| archive/architecture-doc.md | LANDED | Commissions docs/architecture.md as the single orientation doc for the engine. |
| archive/package-review-remediation.md | LANDED | Remediation of the 2026-07-17 seven-reviewer package review; zero high-severity engine findings. |
| archive/readability-review.md | LANDED | Retroactive maintainer-readability pass over the whole bartcore branch diff. |
| archive/retrospective-reviews.md | LANDED (program complete) | Umbrella for a six-review retrospective program; totals: 4 engine defects, 1 calibration defect, roadmap re-ranked. |
| archive/review-perf-followups.md | LANDED (ARC CLOSED) | Lands Tier 5 perf findings + Tier 4 engine notes from the package review. |
| archive/roadmap-survey.md | LANDED (survey delivered) | 15-year BART ecosystem survey vs. dbarts's feature set/backlog; survival promoted to first tier. |
| archive/change-move-fix.md | LANDED | Fixes a detailed-balance violation in the change move for mixed ordinal/categorical splits. |
| archive/composition-mixing-probe.md | KILLED (harm clause fired), 2026-08-10 | Measures whether absorbing a smooth signal share into a parametric block improves the forest's own mixing (tree-mixing-proposals.md sec 4.1's top-ranked candidate); the registered harm clause fired, withdrawing that rank. v2 replaces v1 after eleven blind-critique findings were adopted. |
| archive/convergence-diagnostics.md | LANDED | Adds posterior-package-shaped draws plus an R-hat/ESS summary for multi-chain fits. |
| archive/empty-leaf-veto.md | LANDED (NO-GO on full removal) | Investigates replacing the -1e7 empty-leaf veto with occupancy-aware proposals; keep-and-document decision. |
| archive/grow-from-root.md | LANDED (memo phase) | XBART-style root-down tree sampling; GO on cut-scan-as-shared-primitive and warm-start use, NO-GO as a standalone posterior sampler. |
| archive/grow-from-root-categorical-scan.md | LANDED (S0a-S5, ARC CLOSED), 2026-08-12 | `growTreeFromRoot` places real categorical split rules weighted commensurably with the ordinal branch (exact below a cap, Fisher-prefix above it), inverting the v1 "categoricals never split here" contract; the pre-registered falsifier confirms the shipped/exact draw-law gap, and also finds scan.hpp's naCode split-vs-no-split asymmetry. |
| archive/grow-from-root-default-study.md | KILLED (measured, both strata), 2026-08-08 | Pre-registered study deciding whether bart2's grow-from-root init defaults on; KILLED in both size strata on a per-cell plateau RMSE cost concentrated in noise-heavy/large-n regimes a pooled aggregate hides, each confirmed by a mandated fresh-seed re-run; n.grow.sweeps stays opt-in. |
| archive/grow-from-root-warm-start.md | LANDED | Promotes the cut-scan kernel into an XBART-style root-down warm-start producer (n.grow.sweeps); byte-identical default path. |
| archive/interaction-constraints.md | LANDED f455d7c (2026-07-21) | Per-forest interaction constraints (max-order cap + hard co-occurrence deny/allow groups) via interactions(); post-landing BCF use-after-free fix 04ca425. |
| archive/interaction-constraints-p4.md | LANDED aadbbc8/103dbe2 (2026-07-21) | P4 follow-on: blocks() fixed-capacity per-tree block-additive prior plus the columnMask warm-start/setState feasibility gate (F1; follow-up 073d3db); soft path penalties and formal heredity stay deferred. |
| archive/monotone-prior-draw.md | LANDED 173a710 (2026-08-04) | samplePriorPredictive on a monotone fit drew unconstrained leaf values and left the chain monotone-infeasible; fixed by exact per-tree rejection from the constrained prior. Found by the SBC-extension blind critique. |
| archive/sampletreesfromprior-midchain.md | LANDED 1947b10, 2026-08-08 | Fixes mid-chain sampleTreesFromPrior leaving totalFits/leafOf out of sync with the reset muByTree (a permanent residual displacement and an ASAN-invisible out-of-capacity read); the reset is forest-only and lands the forest in the zero-fit state a freshly built chain carries. |
| archive/setpredictor-leafof-rebuild.md | CLOSED 2026-08-05 | All three leafOf-scoped mechanisms declined by measurement (ceiling +9.7%); the mu[leafOf]-gather SIMD door was taken instead, in a roll-only form, LANDED 9141274. |
| archive/moves-degenerate-root-guard.md | LANDED | Fixes a segfault when a root-only tree has no available split variable. |
| archive/pointwise-loglik.md | LANDED | R-side per-observation/per-draw log-likelihood extract for loo/waic. |
| archive/ppd-sigma-pairing.md | LANDED | Fixes multi-chain sampleFromPPD sigma pairing (chain-interleaved vs chain-blocked). |
| archive/prior-constants.md | LANDED | Provenance/derivation for every adopted prior default, recorded in docs/design/prior-defaults.md + Rd. |
| archive/prior-predictive.md | LANDED | R-level samplePriorPredictive verb (issue #31). |
| archive/sigma-df-zero-weights.md | LANDED | Gaussian sigma posterior df counts only positive-weight rows. |
| archive/vignette-refresh.md | LANDED | Vignettes updated for the 1.0-0 surface (family selection, DART, sparse, missingness, non-constant leaves). |
| archive/warm-starts.md | LANDED | R-level verb installs a saved fit's forest as another sampler's starting state; post-landing valgrind leak fix; cross-grid remap landed 2026-07-22, lifting the different-grid refusal. |
