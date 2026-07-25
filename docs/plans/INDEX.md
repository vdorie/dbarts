# Plans doc index

Manifest of every `docs/plans/*.md` implementation plan (121 files; `README.md`
is the process/contract doc, indexed separately at the bottom, not listed as
a plan). Grouped by cluster/theme. STATUS reflects each doc's live
`Status:`/`## Status` section (or its equivalent closing Landing note) as of
2026-07-21, corrected against a handful of docs whose status header was
missing or stale when an earlier inventory pass read them - see
`docs/README.md` for the reconciliation notes. See `docs/design/INDEX.md`
for the paired design docs.

Columns: `file | STATUS | one-liner`.

## BCF follow-on cluster (pairs with docs/design/bcf.md)

| file | STATUS | purpose |
|---|---|---|
| bcf-b-ridge.md | NO-GO | Treatment-scale (b) ridge interweaving move, the named suspect for the sigma residual; both controls clear it - shelved as an unimplemented future mixing win. |
| bcf-ridge-interweaving.md | LANDED | Prognostic (a) ridge interweaving move, landed 9617c94; confirmed mixing win, sigma SBC flag persists and routes to bcf-sigma-residual. |
| bcf-sigma-residual.md | MIXED (no Status header) | Diagnoses the BCF sigma burn-in transient as slow forest-structure mixing, not a glue-scale defect; recommends `burn = ceiling(72000/thin)`, now live in benchmarks/R/sbc.R, but the doc itself never closes with a Status/Landing note. |
| bcf-testfits-guard.md | LANDED | Guards BCF test-fit/predict entry points that silently reported the bare prognostic forest instead of the combined a*mu+b*tau surface. |

## Forest / multi-forest infrastructure

| file | STATUS | purpose |
|---|---|---|
| forest-combiner.md | LANDED, 2026-07-14 | Extracts BCF's hardcoded glue into a polymorphic `ForestCombiner<L>` so multinomial (and future models) can plug in without re-forking Chain. |
| forest-split-bcf.md | LANDED (two phases) | Splits `Forest<L>` out of Chain and lands BCF as the first two-forest sampler (steps 1-5); a later "Continuation (post data-ownership-4)" phase wires BCF's moderators restriction - both phases complete, one file. |
| multi-forest-models.md | LANDED (tracker; historical value only) | Queue/tracker for the multi-forest family (multinomial, heteroscedastic, hurdle); all three now correctly marked landed in-file. |

## Multinomial cluster (pairs with docs/design/multinomial.md)

| file | STATUS | purpose |
|---|---|---|
| multinomial.md | LANDED (ARC CLOSED, C1-C7) | Base K-forest interleaved-PG-softmax multinomial model (RNG-neutral seams + the model + bart2 surface); hub of the cluster, wrote docs/design/multinomial.md as its C6 commit. |
| multinomial-counts.md | LANDED (C1-C3) | Generalizes to an n x K count matrix via PG(n_i,.) as a sum of PG(1,.) draws; single-trial reduction bitwise. |
| multinomial-formula.md | LANDED (C1-C2) | Formula interface (factor LHS, cbind(...) count matrix) for bart2(family="multinomial"); RNG-neutral. |
| multinomial-margins.md | LANDED | Rewrites per-sweep margin computation O(nK^2) -> O(nK); rounding-level divergence confined to multinomial channels. |
| multinomial-test-surface.md | LANDED (C1-C2; C3 docs commit unconfirmed in-file) | Test-at-creation (x.test softmax probabilities) + predict-on-newdata (K-forest saved-tree replay). |
| multinomial-varcounts.md | LANDED (C1; C2 docs commit unconfirmed in-file) | Per-sample per-category variable-count channel (mbart2-style). |

## Single-chain parallelism & block-fusion

| file | STATUS | purpose |
|---|---|---|
| block-fusion-stage-a.md | SUPERSEDED (historical) | Stage A (b=1 atom-vocabulary refactor) landed and shipped as default, then the whole block-fusion project was later reverted/excised. |
| block-fusion-stage-b.md | NO-GO (historical) | Stage B (b>1 fusion, the intended win); CLOSED WONT-DO, never landed, machinery later excised entirely. |
| blocked-jacobi-trees.md | NO-GO (KILLED as a build target, 2026-07-21) | Noise-split single-chain tree parallelism (b trees update independently per barrier). Phase 0 exact/GO, Phase 1 promising, then a head-to-head reversed the call and a real-engine test refuted an optimistic revival attempt - strictly dominated by straight within-chain threading, also NO-GO. |
| within-chain-threading.md | NO-GO (CLOSED 2026-07-13) | Thin pointer: "see Landing notes and docs/design/within-chain-threading.md section 8" for the full record. Within-chain parallelism for large-n single-chain workloads; step 2 prototype passed bitwise-invariance but missed both speed thresholds. |

## Chi hyperprior cluster

| file | STATUS | purpose |
|---|---|---|
| chi-default-research.md | LANDED | 48-cell simulation study picking the binary-outcome default chi(1.5, Inf) -> chi(1.5, 2); equivalence baselines re-pinned. |
| chi-hyperprior-df.md | LANDED | Fixes chi(df) to sample an actual chi distribution with the stated df instead of a silently-doubled shape; draw-neutral at the default. |
| chi-k-empty-leaf-count.md | LANDED | Excludes leaves stranded empty by a data mutation from the chi-k hyperprior's count/sum-of-squares accounting. |
| chi-k-runaway.md | LANDED | Caps the sampled end-node precision k so an improper/weak prior scale can't let the Gibbs draw run away; bit-identical below the cap. |

## Tau sampler cluster

| file | STATUS | purpose |
|---|---|---|
| tau-slice-review.md | REFERENCE (review) | Reviews the grouped tau slice sampler; verdict KEEP as default, exact-IG replacement recommended for the cauchy prior only. |
| tau-cauchy-exact-ig.md | LANDED | Replaces the cauchy-prior grouped tau slice sampler with an exact Makalic-Schmidt IG two-block Gibbs draw. |
| tau-slice-stepout-cap.md | LANDED | Bounds sliceSampleOnce's step-out loops to prevent hangs; reproduced a 2.5-minute pre-fix hang. |

## Rbart custom-prior cluster

| file | STATUS | purpose |
|---|---|---|
| rbart-loop-profile.md | LANDED (measurement) | Measures rbart_vi custom-prior R loop overhead (12-77% by cell); surfaces the divergence bug as an out-of-scope flag. |
| rbart-loop-fix.md | LANDED | Replaces per-sweep run(0,1) round trips with one per-chain callback hook; draws unchanged. |
| rbart-custom-prior-divergence.md | LANDED | Fixes 15-91x tau divergence between rbart_vi's custom-prior R loop and the in-core path. |

## SIMD / x86 cluster

| file | STATUS | purpose |
|---|---|---|
| simd-survey.md | REFERENCE (READ-ONLY survey) | arm64 SIMD candidate survey; ranks the fused residual-roll kernel and x86 dispatch gaps as safe wins, recommends leaving suffstat scalar for bitwise reproducibility. |
| x86-simd-plan.md | REFERENCE (READ-ONLY investigation) | The x86-measured follow-up to simd-survey.md; finds AVX2 misdetected as AVX (real bug), ranks a suffstat SIMD toggle as the only worthwhile kernel win. |
| simd-flag-multiversioning.md | NO-GO | Evaluates AVX2-flag-built variants of unrolled Mean/Variance; no material end-to-end win. |
| x86-simd.md | CLOSED/SUPERSEDED | Action plan to restore/bench x86 dispatch gaps; all four goal items settled by later work (NT-store deleted, SSR SIMD leaf NO-GO, Windows-arm64 NEON landed, memory-wall-frontier.md is the live perf frontier past SIMD). |

## Survival cluster

| file | STATUS | purpose |
|---|---|---|
| survival-models.md | LANDED (ITEM CLOSED except grouped-surface) | AFT log-normal + discrete-time hazard native survival families. |
| survival-grouped-surface.md | LANDED | Grouped (random-intercept) AFT survival reachable from R; unblocks riAFTBART migration. |

## Mutation-surface cluster

| file | STATUS | purpose |
|---|---|---|
| mutation-fuzzing.md | LANDED | Randomized adversarial property-test fuzzer over setPredictor/setCell/sessions/state round-trips; found two latent state-serialization edge cases. |
| mutation-journal.md | LANDED | Replaces full-array copy-then-restore with build-new-and-swap; journals only changed cells past a break-even threshold. |

## Gate hardening cluster

| file | STATUS | purpose |
|---|---|---|
| gate-blindspot-audit.md | LANDED | Review-3 poison sweep (16 deliberate kernel breakages) + feature x gate coverage matrix; found several blind spots. |
| gate-hardening-1.0.md | LANDED | Implements six new gates targeting the audit's blind spots (bcf-exact-weak, wtgp/grouped/chik2 equivalence, linear-exact, bd-balance). |

## SBC / calibration cluster

| file | STATUS | purpose |
|---|---|---|
| sbc-calibration.md | MIXED (Tier A complete; B/C outstanding) | Simulation-based-calibration harness; found the BCF glue-on sigma mixing issue and the cauchy-tau SBC-intractable tooling gap. 646-line running log, heaviest rewrite-signal doc in the corpus. |
| sbc-ci-gate.md | LANDED | Weekly/manual non-blocking CI run of the gaussian SBC baseline. |

## Data-ownership program (pairs with docs/design/data-ownership.md)

| file | STATUS | purpose |
|---|---|---|
| data-ownership.md | COMPLETE | Hub/design-tracking doc for the 5-plan program; all five landed. |
| data-ownership-1-container.md | LANDED | Plan 1 of 5: engine-owned quantized (u8/u16) predictor container replacing the borrow-and-alias apparatus. |
| data-ownership-2-ingestion.md | LANDED | Plan 2 of 5: data.frame-direct ingestion, the sparseFactor class, narrowing data@x's role. |
| data-ownership-3-mutation.md | LANDED | Plan 3 of 5: mutation-surface rewire to call-time raw supply. |
| data-ownership-4-views.md | LANDED (CLOSED) | Plan 4 of 5: per-forest column views/restriction plus a standalone handle column axis; unblocks forest-split-bcf's moderators argument. |
| data-ownership-5-sparse.md | LANDED (CLOSED) | Plan 5 of 5: CSC sparse-categorical storage kernel; declares the whole 5-plan program COMPLETE. |

## Data-store review cluster (2026-07-17 four-reviewer review; pairs with docs/design/data-store.md)

| file | STATUS | purpose |
|---|---|---|
| data-review-remediation.md | LANDED (ARC CLOSED) | Fixes 3 verified bugs + 1 first-tier perf finding (validation gaps, a quantize loop). |
| data-store-consolidation.md | LANDED (ARC COMPLETE) | Consolidates the C++ data layer (train/test twinning collapsed, explicit per-column storage kind); bench-neutral. |
| data-store-residuals.md | LANDED | Discharges the last 4 deferred findings (shared quantize core, isView provenance-only, conditional raw gather, linear-leaf row-major layout). |

## Response-family / model-surface singletons

| file | STATUS | purpose |
|---|---|---|
| heteroscedastic.md | LANDED (C1-C4) | Second (variance) forest for y=f(x)+s(x)eps via conjugate scaled-inv-chi-sq leaves; homoscedastic fits byte-identical throughout. |
| hurdle.md | LANDED (3 commits) | family="hurdle.lognormal" composed R-side from a probit occupancy fit + a gaussian fit on log(y>0); engine byte-neutral. |
| monotone-bart.md | LANDED (3 commits) | Per-variable monotone constraints via a constrained constant-leaf model; a bias in an intermediate design was caught by the exact-posterior gate before shipping. |
| negative-binomial.md | LANDED (C1-C3) | family="nbinom" via PG augmentation; sweep-order and real-shape-PG issues caught pre-implementation, corrected to integer-dispersion-only. |
| ordinal-outcomes.md | LANDED (C1-C3 + follow-up) | family="ordinal" via cumulative probit with Cowles-style cutpoints; auto-dispatches on ordered factors (fixed a pre-existing silent bart2 bug). |
| robust-errors.md | LANDED (ARC CLOSED) | Student-t residuals via scale-mixture augmentation (resid.dist); both estimated-nu and fixed-nu modes shipped. |
| weighted-binary.md | ACTIVE (parked memo, post-1.0 only) | Preserves analysis for integer-weight probit and arbitrary-real-weight logistic; deliberately not implemented in 1.0-0. |
| weighted-binary-ppd.md | LANDED | Fixes weighted-binary posterior-predictive draws (was degenerate two-point, now coherent binomial); closes the last live path of issue #79. |

## Grouped-effects singletons

| file | STATUS | purpose |
|---|---|---|
| group-by-exposure.md | RESEARCH-OPEN | Decision memo (not yet written) on exposing grouped random effects beyond rbart_vi via a group.by argument; blocked on demand. |
| grouped-equivalence.md | RESEARCH-OPEN | Records why grouped rbart_vi fits are excluded from the equivalence gate; conditional on a future grouped-Gibbs refactor. |

## Performance & kernel singletons

| file | STATUS | purpose |
|---|---|---|
| collapse-merge.md | LANDED | Unifies the two subtree-collapse sites to merge leaf parameters by the same effective-observation-weighted rule. |
| constant-leaf-fits.md | MIXED (bench-sampler compare PENDING) | Replaces the per-tree fits slab with node-indexed mu tables + uint32 leafOf maps for constant leaves; two commits landed bitwise. |
| constant-leaf-suffstat.md | LANDED | Switches ConstantGaussianLeaf to a one-pass, order-insensitive sufficient statistic; a bench amendment fixed a 12-18% regression via kernel unrolling. |
| hot-layer-u8.md | NO-GO (phase 2) | Measured u8 column-width kernels vs u16; no partition-throughput win on arm64, folded per-column widths into data-ownership instead. |
| kernel-cleanups.md | LANDED | Exposes fast serial moment accumulators without the null-thread-manager indirection; adds an injectable misc.a output hook. |
| linear-leaf-reuse.md | LANDED | Per-(tree,node) crossproduct (U'WU) cache for linear leaves; 18-35% sweep speedup measured. |
| parallel-falsifiers.md | LANDED | Measured go/no-go on three parallel-BART-frontier falsifiers (E1/E2/E3); all three survived. |
| test-fit-parallel.md | LANDED | Parallelizes test-fit routing over rows above a cutoff; byte-identical at any thread count, 1.48-3.38x measured. |

## Data-layer / state-format singletons

| file | STATUS | purpose |
|---|---|---|
| cutpoints-shrink-orphan.md | LANDED | Fixes setCutPoints grid-shrink leaving an ordinal split indexing past the new grid. |
| family-on-model.md | LANDED | Moves the response-family slot from dbartsControl to dbartsModel. |
| flat-format-v2.md | LANDED | Replaces FlatNode's bit-cast mask encoding with a tagged value/mask/pool union; state format version bumped to 2. |
| fuzz-state-roundtrip.md | LANDED | Fixes two getState/setState round-trip edges the mutation fuzzer found. |
| range-scaling.md | LANDED | Decision memo (keep range scaling over sd-standardization) plus an updateScale flag. |
| state-continuation.md | LANDED | Drops bitwise-continuation-only state fields in favor of semantic restore. |
| state-format-policy.md | LANDED | Saved states carry a format version + package provenance; incompatible loads refuse cleanly. |
| test-data-parity.md | LANDED (CLOSED) | Test-side data store gains resident sparse storage (no densification); 1.83x-6.98x memory/code-shrink measured. |

## API-surface cluster

| file | STATUS | purpose |
|---|---|---|
| c-api-growth.md | LANDED | Grows the C API with a size-first results struct, additive by-name state blocks, and a per-family loglik channel. |
| capi-callbacks.md | LANDED | Adds a per-sweep conditioning callback and observer hook to the C API. |
| capi-dispatch-table.md | LANDED (arc complete) | ABI-compatibility mechanism for LinkingTo consumers (X-macro single-source stubs + version/hash handshake); both dbarts and stan4bart landed. |
| consumer-spec-surface.md | LANDED, 2026-07-25 | Exports `dbartsSpec()` over a shared internal resolution, so an embedding package builds a sampler specification without reaching into `dbarts:::parsePriors`. |
| interface-review.md | LANDED | Review-2 retrospective auditing the exported R surface; 11 code fixes + 6 doc fixes + 11 taste calls, all landed same day. |
| pre-release-surface-fixes.md | LANDED | Fixes aft-loglik defect + freeze-regret paper-cuts from a pre-release surface audit. |
| r-ingestion-cleanups.md | LANDED | Unifies duplicated classification-family routing, sparse test-matrix handling, missing-policy guards. |
| r5-cleanup.md | LANDED | Removes startThreads/stopThreads no-ops, documents offset sync, hides internal class names from S4 validity errors. |

## Testing / gates / CI singletons

| file | STATUS | purpose |
|---|---|---|
| cran-readiness.md | MIXED (1 NOTE) | R CMD check --as-cran plus ASAN/UBSAN sweep; clean except one persisting compiled-code NOTE. |
| equivalence-ci.md | LANDED | Scheduled CI run of the equivalence gate in statistical mode with baseline provenance (a MANIFEST). |
| gp-cache-test-flake.md | LANDED | Fixed a build-layout-dependent flaky GP-leaf kernel-cache test. |
| post-mutation-assertions.md | LANDED | Replaces near-tautological mutation-test asserts with statistical from-scratch-agreement checks. |
| rchk-ci.md | LANDED | Scheduled rchk (PROTECT-balance) CI job on the bridge; first real run's ~50 findings fixed and re-verified zero. |
| revdep-smoke-ci.md | LANDED | Monthly CI smoke test of stan4bart + bartCause against dev dbarts. |
| small-validation-fixes.md | LANDED | dbartsData precision-degenerate-response warning + raw-matrix incorporate NA-row fix. |
| snapshot-tests.md | LANDED | Per-family reproducibility snapshot files replace whole-file-history-dependent literals. |
| test-suite-trim.md | LANDED | Trims tinytest suite wall time 33.8s -> 20.3s without losing coverage. |
| tests-cpp-split.md | LANDED | Splits the monolithic tests/cpp TU into per-area TUs with dependency tracking. |

## Build / infra singletons

| file | STATUS | purpose |
|---|---|---|
| autoconf-dead-code.md | LANDED | Strips dead configure knobs (--with-xint-size, Solaris probes, orphaned m4); ~3700 net lines deleted. |
| bridge-error-path-leaks.md | LANDED | Wraps bridge entry points so Rf_error longjmps don't leak C++ heap; confirmed zero-leak by valgrind CI. |
| repo-modernization.md | MIXED (recurring/standing item) | CI/tooling hygiene bundle (action pins, sanitizers, optional codecov); sub-items 1-2 landed. |

## Review / retrospective programs

| file | STATUS | purpose |
|---|---|---|
| architecture-doc.md | LANDED | Commissions docs/architecture.md as the single orientation doc for the engine. |
| architecture-numerical-review.md | REFERENCE | Two-reader fresh-eyes stress review of the frozen engine; complete, no blocking issues. |
| correctness-audit.md | REFERENCE | Re-derives every acceptance ratio / conjugate update term-by-term across 7 blocks; all CONFIRMED. |
| package-review-remediation.md | LANDED | Remediation of the 2026-07-17 seven-reviewer package review; zero high-severity engine findings. |
| readability-review.md | LANDED | Retroactive maintainer-readability pass over the whole bartcore branch diff. |
| retrospective-reviews.md | LANDED (program complete) | Umbrella for a six-review retrospective program; totals: 4 engine defects, 1 calibration defect, roadmap re-ranked. |
| review-perf-followups.md | LANDED (ARC CLOSED) | Lands Tier 5 perf findings + Tier 4 engine notes from the package review. |
| roadmap-survey.md | LANDED (survey delivered) | 15-year BART ecosystem survey vs. dbarts's feature set/backlog; survival promoted to first tier. |

## Research doors / decision-gated (no or minimal code; open per TODO)

| file | STATUS | purpose |
|---|---|---|
| gp-followups.md | RESEARCH-OPEN | Placeholder for two deferred GP-leaf extensions (sampled lengthscales, low-rank kernels); blocked on a named consumer. |
| gpu-bart.md | RESEARCH-OPEN (memo DONE) | Survey memo written (docs/design/gpu-bart.md): GO on cut-scan-as-warm-start (delivered), broader GPU experiment left open behind the informed-kernel and CG-leaf prototypes (parallel-bart-frontier next-actions 4-5). |
| interaction-constraints.md | RESEARCH-OPEN | Consumer-gated memo placeholder for structured split-variable control (grouped DART, interaction limits). |
| python-bindings.md | ACTIVE (no Status/Landing section) | Feasibility spike memo for a Python binding over bartcore; unclear if the spike ever ran. |
| sparse-extensions.md | ACTIVE (open, consumer-gated) | Backlog of deferred sparse-column extensions (in-place mutation, streaming range kernel, mixed-column mutation). |

## Misc feature / fix singletons

| file | STATUS | purpose |
|---|---|---|
| change-move-fix.md | LANDED | Fixes a detailed-balance violation in the change move for mixed ordinal/categorical splits. |
| convergence-diagnostics.md | LANDED | Adds posterior-package-shaped draws plus an R-hat/ESS summary for multi-chain fits. |
| empty-leaf-veto.md | LANDED (NO-GO on full removal) | Investigates replacing the -1e7 empty-leaf veto with occupancy-aware proposals; keep-and-document decision. |
| grow-from-root.md | LANDED (memo phase) | XBART-style root-down tree sampling; GO on cut-scan-as-shared-primitive and warm-start use, NO-GO as a standalone posterior sampler. |
| grow-from-root-warm-start.md | LANDED | Promotes the cut-scan kernel into an XBART-style root-down warm-start producer (n.grow.sweeps); byte-identical default path. |
| moves-degenerate-root-guard.md | LANDED | Fixes a segfault when a root-only tree has no available split variable. |
| pointwise-loglik.md | LANDED | R-side per-observation/per-draw log-likelihood extract for loo/waic. |
| ppd-sigma-pairing.md | LANDED | Fixes multi-chain sampleFromPPD sigma pairing (chain-interleaved vs chain-blocked). |
| prior-constants.md | LANDED | Provenance/derivation for every adopted prior default, recorded in docs/design/prior-defaults.md + Rd. |
| prior-predictive.md | LANDED | R-level samplePriorPredictive verb (issue #31). |
| sigma-df-zero-weights.md | LANDED | Gaussian sigma posterior df counts only positive-weight rows. |
| vignette-refresh.md | LANDED | Vignettes updated for the 1.0-0 surface (family selection, DART, sparse, missingness, non-constant leaves). |
| warm-starts.md | LANDED | R-level verb installs a saved fit's forest as another sampler's starting state; post-landing valgrind leak fix. |

## Process doc (not a plan)

| file | purpose |
|---|---|
| README.md | Process/contract doc: roles (Fable/Opus/Sonnet), the plan-file template, RNG gate classes and their required gates, brevity rubric, review checklist. Does not enumerate the directory's contents - that is this INDEX's job. |
