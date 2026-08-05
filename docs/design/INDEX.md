# Design doc index

Manifest of every `docs/design/*.md` (37 files), grouped by theme. STATUS is
read from each doc's live `Status:` line or `## Status` section as of
2026-08-04, when a consistency pass re-verified every row against its doc;
a few standing-reference docs carry no status line by design
(noted as REFERENCE). See `docs/README.md` for how this index relates to the
root TODO and `docs/plans/README.md`.

Columns: `file | STATUS | one-liner`.

## Multi-forest / BCF

| file | STATUS | purpose |
|---|---|---|
| bcf.md | LANDED | Two-forest (prognostic mu + treatment tau) Bayesian Causal Forests sampler (Hahn/Murray/Carvalho 2020); landed 2026-07-07, mixing refinements through 07-10. Plans: bcf-b-ridge, bcf-ridge-interweaving, bcf-sigma-residual, bcf-testfits-guard, forest-split-bcf. |
| forest-combiner.md | LANDED, 2026-07-14 | Promotes BCF's hardcoded Chain-side glue into a polymorphic `ForestCombiner<L>` so multinomial/future multi-forest models reuse it; null-short-circuit costs single-forest chains nothing. Plan: forest-combiner.md. |

## Response families

| file | STATUS | purpose |
|---|---|---|
| heteroscedastic.md | LANDED, 2026-07-20 | HBART: second (variance) forest, s^2(x) via conjugate scaled-inv-chi-sq leaves (multiplicative-positive trees); reuses existing conjugate move machinery. Plan: heteroscedastic.md (+ multi-forest-models.md tracker). |
| monotone.md | LANDED, 2026-07-19 | Per-variable monotone constraints (mBART, Chipman-George-McCulloch-Shively 2022) via a constrained constant leaf; first maintained implementation (no CRAN mBART). Plan: monotone-bart.md. |
| ordinal.md | LANDED, 2026-07-18 | Ordered-categorical responses via cumulative probit, K-1 sampled cutpoints (Cowles marginal MH); K=2 collapses to bit-for-bit probit. Plan: ordinal-outcomes.md. |
| multinomial.md | LANDED, 2026-07-15/17 | First multi-forest model: K symmetric constant-leaf forests coupled by softmax via interleaved one-vs-rest Polya-Gamma augmentation + level-centering. Plan: multinomial.md + the multinomial-* cluster. |
| negative-binomial.md | LANDED, 2026-07-18 | family="nbinom" via Polya-Gamma NB augmentation; r restricted to positive-integer shape (exact PG only exists there), continuous r deferred (TODO: negbin-real-dispersion). Plan: negative-binomial.md. |
| survival.md | LANDED | AFT log-normal (07-10, reduces exactly to Gaussian path) + discrete-time hazard (07-18, person-period sugar over existing binary families). Plan: survival-models.md (+ survival-grouped-surface.md). |
| hurdle.md | LANDED, 2026-07-20 | Two-part/hurdle model; the zero/positive split is OBSERVED not latent, so the two forests are conditionally independent - composed R-side as two ordinary fits, zero engine/bridge/state-format code. Plan: hurdle.md. |
| weighted-logistic.md | LANDED, 2026-07-05 | Extends observation weights to logistic responses (previously refused) via PG(w,psi) "multiple copies"; positive-integer weights only. Plans: weighted-binary.md (declined follow-ups), weighted-binary-ppd.md. |
| grouped-random-effects.md | LANDED, 2026-07-04 | Moves rbart_vi's random-intercept Gibbs loop into the engine as a `GroupedResponse` decorator composable with any response family; R loop kept as custom-prior fallback. |
| forest-ranef-interweaving.md | NO-GO (recorded door, reconfirmed 2026-07-20) | Investigates a joint/interweave move for forest-ranef confounding; no cheap ASIS move exists (forest exposes no per-group scalar), only fix is a full marginal collapse - not worth it for a narrow weak-signal payoff. |
| correlated-outcomes.md | RESOLVED, 2026-07-22 (decision-gated door) | Richer error/outcome covariance around a BART mean, for causal DiD/panel and multi-outcome designs: (B) multivariate/SUR shipped as mvbart() in stan4bart with ZERO dbarts engine change (commit e27a7c3); (A) AR-1 serially-correlated errors deferred. TODO: correlated-outcomes. |

## Performance & parallelism frontier

| file | STATUS | purpose |
|---|---|---|
| memory-wall-frontier.md | CLOSED as idea map; recommended lever LANDED; re-profiled 2026-08-04 | Three-lens (algorithms/hardware/systems) ranked ideation map of the per-sweep memory-wall bottleneck; lands on fp32 residual storage as the answer, GPU and block-fusion re-tread both ruled dead. Sec 10 re-profile: the histogram-fused-suffstat flagship's premise survives (VD fork on a falsifier), Z-order/locality reordering measured dead. |
| parallel-bart-frontier.md | MIXED (research survey; 3 falsifiers measured 07-08) | Ranks every axis of BART parallelism beyond the settled observation-axis; structural facts, conservation laws, surviving constructions (block-fused atoms, delayed acceptance, coupled chains), CRAN GPU axis. Plan: parallel-falsifiers.md. |
| within-chain-threading.md | CLOSED, NO-GO on x86 and Apple Silicon, 2026-07-21 | std::barrier worker pool parallelizing single-chain large-n sweeps; ~1.5-2.1x ceiling, missed both pre-registered gates; reconfirmed under fp32 storage and after a real-engine test refuted an optimistic microbench revival. Plan: within-chain-threading.md. |
| reduced-precision-storage.md | LANDED / COMPLETE, 2026-07-20/21 | Optional compile-time narrowed hot-path storage: Track 1 (uint32 indices, bitwise-preserving) LANDED; Track 2 (opt-in fp32 residual, gaussian constant leaf) v1 LANDED (1.10-1.30x), v2 (fp32 scratch/fits bundle) MEASURED NO-GO. |
| block-fusion.md | CLOSED, WONT-DO | Block-fused sub-sweeps over a persistent atom map meant to amortize O(n) DRAM passes across b trees (~6x target). Stage A (b=1 refactor) landed; Stage B (b>1) benchmarked 4-9x SLOWER (cost model over-counted amortizable traffic ~3x); machinery excised to archive/block-fusion. Plans: block-fusion-stage-a.md, -b.md. |
| gpu-bart.md | NO-GO (survey; no direction earns a prototype yet) | Surveys 7 GPU-acceleration directions rated by ceiling/engine-fit/cost; recommends no GPU prototype now, ranks grow-from-root's cut-scan kernel as best future candidate. Plan: gpu-bart.md. |

## Leaf models

| file | STATUS | purpose |
|---|---|---|
| linear-leaves.md | LANDED, 2026-07-04 | Per-leaf small Bayesian linear regression (treed regression / BCF varying-slope) as the second leaf model beyond constant; all four stages + xbart/view deferral landed same day. |
| gp-leaves.md | LANDED (Part 1, stages 1-4); Part 2 unscheduled | GP leaves as a third leaf-parameter shape on the conjugate engine (Part 1: engine, formats, R surface, chi-k+kernel caching); Part 2 (non-conjugate MoveStrategy for non-Gaussian likelihoods) deliberately unscheduled pending a real consumer. Plan: gp-followups.md (open extensions). |

## Data layer

| file | STATUS | purpose |
|---|---|---|
| data-ownership.md | COMPLETE | Owned/quantized predictor container design (never-retain-raw-by-default, call-time raw supply, shared-handle views); all 5 sub-plans landed. Plans: data-ownership.md (hub) + data-ownership-1..5-*.md. |
| data-layout.md | CLOSED - SHELVED (re-evaluated 2026-08-04) | Contiguous per-node response-axis (y/fits/weights/residual) reorder; the ~10% standalone win was re-evaluated on the 06f73b0 re-profile after block-fusion's kill orphaned its substrate framing - the fit-scatter consumer is gone and the reorder pays a same-magnitude gather, so the realistic win falls below the marginal ~10%; dated close-out in the doc's post-mortem. |
| data-store.md | REFERENCE | Standing technical reference for ColumnStore/CodeBlock/ColumnSource, the predictor mutation transaction, and ownership/borrow rules; required reading before data-adjacent engine work. Plans: data-review-remediation.md, data-store-consolidation.md, data-store-residuals.md. |
| sparse-columns.md | LANDED, 2026-07-04 | Rank-bitmap CSC sparse predictor-column representation for wide mostly-zero designs; pure-sparse and mixed dense/sparse both landed, densifies above 20% density; in-place nonzero mutation (extension (i)) landed 2026-07-22. Plan: sparse-extensions.md (remaining deferred extensions). |
| pooled-masks.md | LANDED, 2026-07-04 | Lifts the categorical-predictor cap from 53 to 65535 levels via pooled masks; no change to columns of <=53 levels. |
| mia-missingness.md | LANDED, 2026-07-04 | Missing Incorporated in Attributes (Twala 2008): every split gains a learned missing-direction so predictor NAs route through splits. |

## Core/infra & surface

| file | STATUS | purpose |
|---|---|---|
| core-generalization.md | LANDED (mostly; phase 6 open) | Founding design for the bartcore C++20 rewrite: layering, concept decomposition (ResponseModel/LeafModel/MoveStrategy/SplitSelector/Forest), phase-by-phase log. Cutover complete 2026-07-03; phase 6 (GP leaves/non-conjugate moves) tracked via gp-leaves.md. Self-identifies as historical reference; docs/architecture.md is authoritative for current state. |
| consumer-spec-surface.md | LANDED, 2026-07-25 | Exports `dbartsSpec()`, which resolves the (control, model, data) triple and family token without constructing a sampler, so a LinkingTo consumer reaches the whole single-forest feature set through supported surface instead of an unexported internal; prior evaluation environment narrowed 2026-07-27 (b96d3bb, section 7). Plan: consumer-spec-surface.md. |
| public-surface.md | MIXED (reviewed 07-03, updated through 07-06; decisions recorded inline) | The major-version public R surface: engine cutover sequencing, factor/categorical ingestion, response-family exposure, DART, standalone data handle, C API/callbacks; section 7 catalogs features deferred here that have since landed elsewhere. |
| kernel-vocabulary.md | REFERENCE | Normative contract between the generic BART core and the compiled `misc.a` kernel library: dispatch mechanism, current vocabulary, planned additions. Plan: kernel-cleanups.md. |
| robust-errors.md | LANDED, 2026-07-17 | Student-t residuals for continuous responses via scale-mixture augmentation (resid.dist), riding the existing workingWeights hook. Plan: robust-errors.md. |
| prior-defaults.md | REFERENCE | Plain record of every current prior default (k, power/base, sigdf/sigquant, node.scale, n.trees, dart delay, tau slice steps) with source and rationale; a record, not a re-derivation. |
| change-move-balance.md | LANDED, 2026-07-08 | Fixes a since-origin detailed-balance defect in the tree change move (missing proposal-density correction, biased toward low-cardinality variables); hybrid ordinal/categorical fix. Plan: change-move-fix.md. |
| empty-leaf-veto.md | LANDED (keep-and-document) | Investigates occupancy-aware tree proposals vs. the existing finite-penalty empty-leaf veto; keeps the veto, corrected 2026-07-15 to -HUGE_VAL (was a beatable finite -1e7). Plan: empty-leaf-veto.md. |
| grow-from-root.md | MIXED (GO on cut-scan/warm-start; NO-GO on standalone sampler) | XBART-style root-down stochastic tree construction; ships only as a validity-free warm-start producer, rejected as a standalone approximate posterior sampler (not MH-exact). Plans: grow-from-root.md, grow-from-root-warm-start.md. |
| interaction-constraints.md | LANDED, 2026-07-21 (f455d7c; P4 aadbbc8/103dbe2/073d3db) | Per-forest limits on which predictors may jointly shape the fit: cap interaction ORDER (1 = additive/GAM, 2 = pairwise) and/or DENY named co-occurrences or ALLOW declared groups; motivates calibrated-additivity causal use (additive tau, free mu). Plans: interaction-constraints.md, interaction-constraints-p4.md. |
