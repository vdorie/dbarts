# Design doc index

Manifest of every `docs/design/*.md`, grouped by theme. STATUS is read from
each doc's live `Status:` line or `## Status` section. A few standing-
reference docs carry no status line by design (noted as REFERENCE); a
handful of landed-feature docs likewise carry none - their STATUS here is
drawn from the doc's own opening text instead. See `docs/README.md` for how
this index relates to the root TODO and `docs/plans/README.md`.

Columns: `file | STATUS | one-liner`.

## Multi-forest / BCF

| file | STATUS | purpose |
|---|---|---|
| bcf.md | LANDED | Adds a two-forest (prognostic + treatment) Bayesian Causal Forests sampler (Hahn/Murray/Carvalho 2020), reachable through `forests =`. |
| forest-combiner.md | LANDED, 2026-07-14 | Generalizes BCF's forest-combining logic so multinomial and other multi-forest models can reuse it, at no cost to single-forest fits. |
| multiplier-combiner.md | LANDED, gaussian/probit/logistic, 2026-08-13/14 | Generalizes BCF's two-forest coupling to K forests with a per-forest basis and amplitude, for gaussian, probit and logistic responses (aft/ordinal/nbinom refused by name); reachable via `forests =` or a `bart2()` formula `forest()` term. |
| model-space-survey.md | COMPLETE (survey), 2026-08-08; D1 subsequently LANDED | A research survey of which multi-forest mutation shapes (whole-data `setData`, per-forest row subsetting) are worth building; found and fixed a variance-forest mutation-routing defect along the way, left the rest open. |

## Response families

| file | STATUS | purpose |
|---|---|---|
| heteroscedastic.md | LANDED, 2026-07-20 | Adds a heteroscedastic variance forest (`variance =`), modeling s^2(x) as a second tree ensemble. |
| monotone.md | LANDED, 2026-07-19 | Adds per-variable monotonicity constraints (`monotone =`) via a constrained constant leaf. |
| ordinal.md | LANDED, 2026-07-18 | Adds ordered-categorical responses (`family = "ordinal"`) via cumulative probit with sampled cutpoints; K = 2 reduces bitwise to probit. |
| multinomial.md | LANDED, 2026-07-15/17 | Adds multinomial responses (`family = "multinomial"`): K constant-leaf forests coupled by a softmax link. |
| negative-binomial.md | LANDED, 2026-07-18 | Adds negative-binomial counts (`family = "nbinom"`) via Polya-Gamma augmentation; the dispersion `r` is a positive integer only, real-valued `r` deferred. |
| survival.md | LANDED | Adds two survival families: AFT log-normal (`family = "aft"`) and discrete-time hazard (`family = "hazard"`, a person-period expansion over the binary families). |
| hurdle.md | LANDED, 2026-07-20 | Adds semicontinuous two-part/hurdle responses (`family = "hurdle.lognormal"`), composed in R from two ordinary fits with no engine changes. |
| weighted-logistic.md | LANDED, 2026-07-05 | Lets logistic responses take observation weights (previously refused), as positive-integer replicate counts. |
| grouped-random-effects.md | LANDED, 2026-07-04 | Moves `rbart_vi()`'s random-intercept sampling into the engine, composable with any response family; the R-level loop stays as a fallback for a custom prior. |
| forest-ranef-interweaving.md | NO-GO (recorded door, reconfirmed 2026-07-20) | Investigated a mixing fix for forest/random-effect confounding; no cheap fix exists, and a full fix was not judged worth its cost. Nothing changed for users. |
| correlated-outcomes.md | RESOLVED, 2026-07-22 (decision-gated door) | Investigated richer error covariance around a BART mean; the multivariate/SUR case shipped as `mvbart()` in stan4bart with no dbarts engine change, AR-1 serial correlation stays deferred. |

## Performance & parallelism frontier

| file | STATUS | purpose |
|---|---|---|
| memory-wall-frontier.md | CLOSED as idea map; recommended lever LANDED; re-profiled 2026-08-04 | Surveyed the per-sweep memory-wall bottleneck; its recommended fix (fp32 residual storage) landed, GPU and a block-fusion retread were ruled out. |
| parallel-bart-frontier.md | MIXED (research survey; 3 falsifiers measured 07-08) | Surveyed BART parallelism beyond per-chain parallelism; ranks the surviving directions (block-fused atoms, delayed acceptance, coupled chains) without recommending any for a prototype yet. |
| within-chain-threading.md | CLOSED, NO-GO on x86 and Apple Silicon, 2026-07-21 | Tried a worker-pool for within-chain parallelism on large-n single-chain sweeps; measured too small a speedup on both x86 and Apple Silicon to ship, confirmed by a second attempt. |
| reduced-precision-storage.md | LANDED / COMPLETE, 2026-07-20/21 | Adds optional narrowed hot-path storage: bitwise-preserving uint32 indices shipped, and an opt-in fp32 residual (`storage = "single"`) shipped for the gaussian constant leaf; a further fp32 scratch bundle was tried and measured not worth it. |
| block-fusion.md | CLOSED, WONT-DO | Tried block-fused sub-sweeps to cut memory traffic; the single-tree refactor shipped as the default, but fusing multiple trees measured 4-9x slower than hoped and was removed. |
| gpu-bart.md | NO-GO (survey; no direction earns a prototype yet) | Surveyed seven GPU-acceleration directions; recommends none for a prototype yet, ranking grow-from-root's cut-scan kernel the best future candidate. |

## Leaf models

| file | STATUS | purpose |
|---|---|---|
| linear-leaves.md | LANDED, 2026-07-04 | Adds a per-leaf linear-regression leaf model (`node.prior = linear(...)`), the second leaf model beyond constant. |
| gp-leaves.md | LANDED (Part 1, stages 1-4); Part 2 unscheduled | Adds a per-leaf Gaussian-process leaf model (`node.prior = gp(...)`); a non-conjugate extension for non-Gaussian likelihoods is designed but not built, pending a real consumer. |

## Data layer

| file | STATUS | purpose |
|---|---|---|
| data-ownership.md | COMPLETE | Redesigned the predictor container to own and quantize its data rather than borrow and re-alias it. |
| data-layout.md | CLOSED - SHELVED (re-evaluated 2026-08-04) | Tried reordering per-node storage for a memory-bandwidth win; re-measured at about 10%, below the bar, and shelved. |
| data-store.md | REFERENCE | Standing technical reference for the predictor store's layout, mutation transaction, and ownership rules; required reading before data-adjacent engine work. |
| sparse-columns.md | LANDED, 2026-07-04 | Adds a sparse-column representation for wide, mostly-zero designs (dense/sparse mixed input accepted), densifying automatically above 20% nonzero density. |
| pooled-masks.md | LANDED, 2026-07-04 | Lifts the categorical-predictor level cap from 53 to 65535 via pooled masks; the shipped inline/pooled boundary is 63 categories. |
| mia-missingness.md | LANDED, 2026-07-04 | Adds Missing Incorporated in Attributes: every split learns a missing-value direction so predictor NAs route through splits instead of being dropped. |

## Core/infra & surface

| file | STATUS | purpose |
|---|---|---|
| multinomial-mutation-arc.md | LANDED, 2026-08-24 | Gives multinomial responses a real sampler surface, constructed directly rather than through a placeholder host object. |
| core-generalization.md | LANDED (mostly; phase 6 open) | The founding design for the bartcore C++20 engine rewrite; docs/architecture.md is now authoritative for current state. |
| bart-landscape.md | SNAPSHOT 2026-08-12 | A survey of 35 BART implementations (R, Python, research) and how each compares to dbarts's own feature set. |
| feature-matrix.md | LIVING REFERENCE, cites are symbols checked live by `tools/check-doc-freshness.R` | A living reference of what each of the nine shipped response models, plus four couplings and decorations, can and cannot do, kept current at every landing. |
| consumer-spec-surface.md | LANDED, 2026-07-25 | Exports `dbartsSpec()`, letting a `LinkingTo` consumer resolve a sampler specification without building an unexported internal one. |
| r-c-division.md | ACCEPTED (VD 2026-08-11) | States the rule for what belongs in R versus C++ ("R addresses the conditionals; C++ addresses the integrand") and prices the adoption slate that rule implies (docs/plans/adoption-slate.md). |
| bart-as-a-component.md | LANDED, 2026-08-19 | Documents the contract a driver loop can rely on when embedding dbarts inside a larger sampler: what a mutation does and does not carry, and which multi-forest mutations are legal. |
| public-surface.md | MIXED (reviewed 2026-07-03, updated through 2026-08-13; decisions recorded inline) | The major-version public R surface: engine cutover, factor/categorical ingestion, response-family exposure, DART, the standalone data handle, and the C API/callbacks. |
| error-style.md | ADOPTED for new messages, 2026-08-17 | Sets the error-message style rule (quoting, case, templates per refusal kind) for new messages and the warning-class taxonomy under `dbartsWarning`, following published and measured practice from base R and CRAN packages. |
| kernel-vocabulary.md | REFERENCE | Standing reference for the contract between the generic BART core and the compiled kernel library (`misc.a`). |
| robust-errors.md | LANDED, 2026-07-17 | Adds outlier-robust Student-t residuals (`resid.dist = student(...)`) via scale-mixture augmentation. |
| prior-defaults.md | REFERENCE | A plain record of every current prior default and its source. |
| active-rows-mask.md | REFERENCE | Standing reference for the per-observation 0/1 active-row mask (`$setActiveRows`): its contract, how each response family composes it, and what it leaves untouched. |
| nameable-calibration.md | ARC COMPLETE | Lets a fit name its per-forest leaf-prior scale directly in response units (`prior.scale`), with a matching `$getCalibration`/`$setCalibration` pair, also reachable through the flat C API. |
| change-move-balance.md | LANDED, 2026-07-08 | Fixes a detailed-balance defect in the tree change move that biased splits toward low-cardinality variables. |
| empty-leaf-veto.md | DECIDED (keep-and-document), 2026-07-07 | Keeps the existing empty-leaf veto rather than replacing it with occupancy-aware proposals; the veto was later made ranked rather than flat (docs/architecture.md, "Tree moves"). |
| grow-from-root.md | MIXED (GO on cut-scan/warm-start; NO-GO on standalone sampler) | XBART-style root-down tree construction; ships only as a warm-start producer, not as a standalone sampler (it is not MH-exact). |
| grow-from-root-default.md | KILLED (measured), 2026-08-08 | Measured whether the XBART grow-from-root warm start should default on; it costs accuracy in noisy/large-n settings, so `n.grow.sweeps` stays opt-in. |
| tree-mixing-proposals.md | COMPLETE (survey), 2026-08-09; ADDENDUM sec 12, 2026-08-10 | Surveys tree-space proposal mechanisms beyond the shipped four moves; ranks a same-variable cut move ("perturb") as the best candidate, without shipping one. |
| interaction-constraints.md | LANDED, 2026-07-21 | Adds per-forest interaction constraints (`interactions =`): a max order cap and/or named co-occurrence deny/allow rules. |
| threaded-predict.md | LANDED, 2026-08-25 | Wires `predict()`'s `n.threads` argument to a real thread fan-out over (chain, draw), bitwise identical at every thread count. |
