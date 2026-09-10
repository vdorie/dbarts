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
| model-space-survey.md | COMPLETE (survey), 2026-08-08 | A research survey of which multi-forest mutation shapes (whole-data `setData`, per-forest row subsetting) are worth building; the multi-forest predictor-mutation door it surfaced shipped, and both surveyed data-swap shapes stay undesigned. |

## Response families

| file | STATUS | purpose |
|---|---|---|
| heteroscedastic.md | LANDED, 2026-07-20 | Adds a heteroscedastic variance forest (`variance =`), modeling s^2(x) as a second tree ensemble. |
| monotone.md | LANDED, 2026-07-19 | Adds per-variable monotonicity constraints (`monotone =`) via a constrained constant leaf. |
| ordinal.md | LANDED, 2026-07-18 | Adds ordered-categorical responses (`family = "ordinal"`) via cumulative probit with sampled cutpoints; K = 2 reduces bitwise to probit. |
| multinomial.md | LANDED, 2026-07-15/17 | Adds multinomial responses (`family = "multinomial"`): K constant-leaf forests coupled by a softmax link. |
| negative-binomial.md | LANDED, 2026-07-18 | Adds negative-binomial counts (`family = "nbinom"`) via Polya-Gamma augmentation; the dispersion `r` is a positive integer only, real-valued `r` deferred. |
| survival.md | LANDED | Adds two survival families: AFT log-normal (`family = "aft"`) and discrete-time hazard (`family = "hazard"`, a person-period expansion over the binary families). |
| aft-variance-forest.md | LANDED, 2026-09-06 | Lifts the refusal of a variance forest under `family = "aft"`, so log-normal AFT carries a covariate-dependent dispersion s(x); the latent-channel reason the refusal gave does not hold for aft. |
| aft-status-setter.md | LANDED slices 1 and 2, 2026-09-07; slices 3-4 PROPOSED | Lets a live `family = "aft"` sampler take a new censoring status through `$setResponse(y, status =)`, so the structure stops being fixed at creation; enables the aft and heteroscedastic SBC arms. |
| hurdle.md | LANDED, 2026-07-20 | Adds semicontinuous two-part/hurdle responses (`family = "hurdle.lognormal"`), composed in R from two ordinary fits with no engine changes. |
| weighted-logistic.md | LANDED, 2026-07-05 | Lets logistic responses take observation weights as positive-integer replicate counts. |
| grouped-random-effects.md | RETIRED, 2026-09-06 | Historical record of the in-engine random-intercept sampler; grouped random effects are removed from dbarts (retire-grouped-random-effects.md). |
| forest-ranef-interweaving.md | SUPERSEDED, 2026-09-06 | Investigated a mixing fix for forest/random-effect confounding; NO-GO when recorded, and the door closed with the retirement of grouped random effects (retire-grouped-random-effects.md). |
| correlated-outcomes.md | RESOLVED, 2026-07-22 (decision-gated door) | Investigated richer error covariance around a BART mean; the multivariate/SUR case shipped as `mvbart()` in stan4bart with no dbarts engine change, AR-1 serial correlation stays deferred. |
| retire-grouped-random-effects.md | LANDED 2026-09-06 (1e5f80b2; stan4bart tau-mixing bar and bartCause group.by route remain release prerequisites) | Retires grouped random intercepts from dbarts entirely - `rbart_vi()`, the `GroupedResponse` decorator, its bridge and its two `dbarts_results` fields - and makes stan4bart the home for multilevel structure; the speed comparison fired its tau-mixing gate on two of three seeds and the decision stands, against two sister-repo release prerequisites. |

## Performance & parallelism frontier

| file | STATUS | purpose |
|---|---|---|
| memory-wall-frontier.md | CLOSED as idea map; recommended lever LANDED; re-profiled 2026-08-04 | Surveyed the per-sweep memory-wall bottleneck; its recommended fix (fp32 residual storage) landed, GPU and block fusion were ruled out. |
| parallel-bart-frontier.md | MIXED (research survey; 3 falsifiers measured 07-08) | Surveyed BART parallelism beyond per-chain parallelism; ranks the surviving directions (block-fused atoms, delayed acceptance, coupled chains) without recommending any for a prototype yet. |
| memory-footprint.md | VALIDATED, 2026-09-09 | The closed-form memory model: every engine and R-layer allocation by component, symbol, scope and bytes per unit, with two worked reference cases. |
| within-chain-threading.md | CLOSED, NO-GO on x86 and Apple Silicon, 2026-07-21 | Tried a worker-pool for within-chain parallelism on large-n single-chain sweeps; measured too small a speedup on both x86 and Apple Silicon to ship. |
| reduced-precision-storage.md | LANDED / COMPLETE, 2026-07-20/21 | Adds optional narrowed hot-path storage: bitwise-preserving uint32 indices shipped, and an opt-in fp32 residual (`storage = "single"`) shipped for the gaussian constant leaf; a further fp32 scratch bundle was tried and measured not worth it. |
| block-fusion.md | CLOSED, WONT-DO | Tried block-fused sub-sweeps to cut memory traffic; the single-tree refactor shipped as the default, and fusing multiple trees measured 4-9x slower than hoped, so nothing fused ships. |
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
| pooled-masks.md | LANDED, 2026-07-04 | Raises the categorical-predictor level cap to 65535 via pooled masks; the shipped inline/pooled boundary is 63 categories. |
| mia-missingness.md | LANDED, 2026-07-04; AMENDED 2026-09-09 by [front-door](../plans/front-door.md#front-door) S2 (`missing` argument retired, incorporation unconditional) | Adds Missing Incorporated in Attributes: every split learns a missing-value direction so predictor NAs route through splits instead of being dropped. |

## Core/infra & surface

| file | STATUS | purpose |
| engine-generics-review.md | MEMO for VD, 2026-09-08 | Independent review of what the engine's generic axes (leaf kind, family, forest coupling, host hooks) should be, derived from the model space; critiques the current shape. No code until VD reads it. |
|---|---|---|
| multinomial-mutation-arc.md | LANDED, 2026-08-24 | Gives multinomial responses a sampler surface constructed directly, like every other family. |
| core-generalization.md | LANDED (mostly; phase 6 open) | The founding design for the bartcore C++20 engine rewrite; docs/architecture.md is authoritative for current state. |
| bart-landscape.md | SNAPSHOT 2026-08-12 | A survey of 35 BART implementations (R, Python, research) and how each compares to dbarts's own feature set. |
| feature-matrix.md | LIVING REFERENCE, cites are symbols checked live by `tools/check-doc-freshness.R` | A living reference of what each of the nine shipped response models, plus four couplings and decorations, can and cannot do, kept current at every landing. |
| consumer-spec-surface.md | LANDED, 2026-07-25 | Exports `dbartsSpec()`, letting a `LinkingTo` consumer resolve a sampler specification without building an unexported internal one. |
| r-c-division.md | ACCEPTED (VD 2026-08-11) | States the rule for what belongs in R versus C++ ("R addresses the conditionals; C++ addresses the integrand") and prices the adoption slate that rule implies (docs/plans/adoption-slate.md). |
| bart-as-a-component.md | LANDED, 2026-08-19 | Documents the contract a driver loop can rely on when embedding dbarts inside a larger sampler: what a mutation does and does not carry, and which multi-forest mutations are legal. |
| per-draw-callbacks.md | PLANNED, docs/plans/per-draw-callbacks.md, 2026-09-10 | Proposes a per-draw C callback fired from each chain's worker thread over every channel `storeSample` settles, an R `callback` argument taking a function and context external pointer, a `keepFits` storage opt-out, and a worked Rcpp example, so a large fit reduces draws in place instead of materializing the n x draws x chains array the memory audit found to be five sixths of its peak; reverses dec-B62's worker-thread refusal for an observer hook, and all nine forks are settled (dec-B114). |
| public-surface.md | MIXED (reviewed 2026-07-03, updated through 2026-08-13; decisions recorded inline) | The major-version public R surface: engine cutover, factor/categorical ingestion, response-family exposure, DART, the standalone data handle, and the C API/callbacks. |
| error-style.md | ADOPTED for new messages, 2026-08-17 | Sets the error-message style rule (quoting, case, templates per refusal kind) for new messages and the warning-class taxonomy under `dbartsWarning`, following published and measured practice from base R and CRAN packages. |
| kernel-vocabulary.md | REFERENCE | Standing reference for the contract between the generic BART core and the compiled kernel library (`misc.a`). |
| robust-errors.md | LANDED, 2026-07-17; AMENDED 2026-09-09 by [front-door](../plans/front-door.md#front-door) S2 (`resid.dist` retired to a tombstone) | Adds outlier-robust Student-t residuals (`family = student(...)`) via scale-mixture augmentation. |
| prior-defaults.md | REFERENCE | A plain record of every current prior default and its source. |
| active-rows-mask.md | REFERENCE | Standing reference for the per-observation 0/1 active-row mask (`$setActiveRows`): its contract, how each response family composes it, and what it leaves untouched. |
| nameable-calibration.md | ARC COMPLETE | Lets a fit name its per-forest leaf-prior scale directly in response units (`prior.scale`), with a matching `$getCalibration`/`$setCalibration` pair, also reachable through the flat C API. |
| change-move-balance.md | LANDED, 2026-07-08 | Fixes a detailed-balance defect in the tree change move that biased splits toward low-cardinality variables. |
| empty-leaf-veto.md | DECIDED (keep-and-document), 2026-07-07 | Keeps the empty-leaf veto rather than replacing it with occupancy-aware proposals; the veto is ranked (docs/architecture.md, "Tree moves"). |
| grow-from-root.md | MIXED (GO on cut-scan/warm-start; NO-GO on standalone sampler) | XBART-style root-down tree construction; ships only as a warm-start producer, not as a standalone sampler (it is not MH-exact). |
| grow-from-root-default.md | KILLED (measured), 2026-08-08 | Measured whether the XBART grow-from-root warm start should default on; it costs accuracy in noisy/large-n settings, so `n.grow.sweeps` stays opt-in. |
| tree-mixing-proposals.md | COMPLETE (survey), 2026-08-09; ADDENDUM sec 12, 2026-08-10; ADDENDUM sec 13, 2026-09-06; ADDENDUM sec 14, 2026-09-06; ADDENDUM sec 15, 2026-09-07; ADDENDUM sec 16, 2026-09-07 | Surveys tree-space proposal mechanisms beyond the shipped four moves; ranks a same-variable cut move ("perturb") as the best candidate, without shipping one. Section 13 measures the then-shipped mixture against birth/death-only and against change-without-swap and finds no difference on accuracy, coverage, ESS or variable inclusion; section 14 re-runs those arms across a response swap, where the change move does separate and the swap move does not - the finding that dropped swap out of the default. Section 15 is a four-lens proposal brainstorm with a refutation pass: it establishes that the integrated likelihood reads a leaf only through its sufficient statistics (so a partition-preserving move has likelihood ratio exactly 1, the monotone leaf excepted) and that the cut scan gives exact neighbourhood scores where the discrete-MCMC literature pays a Taylor surrogate, then ranks ten mechanisms with a validity, cost, novelty and falsifier column each, recommending none. Section 16 is the first-principles sequel section 15 asks for, run under a novelty gate: it establishes that the tree prior reads no predictor value and that a tree induces the identical row partition in every chain of one sampler, derives from the move census that shipped trees carry two to three leaves, and ranks eight more mechanisms - an exact Gibbs draw on the level-reallocation fibre, lifted cut and birth/death moves, a same-temperature one-tree exchange between chains, a pairwise-collapsed split transfer - with a refutation verdict on each and no recommendation. |
| swap-removal.md | LANDED, 2026-09-07; AMENDED 2026-09-07 (swap restored at default zero, section 9) | Records the decision to remove the swap tree-proposal before 1.0 (its no-op rate, and the response-swap recovery run where change without swap matches the shipped default), the file-by-file removal, the `proposal.probs` surface defaulting to birth_death 0.6 / change 0.4, and the one bundled baseline re-record it takes; section 9 records the partial reversal, the move returning to the kernel at a default of zero because a one-tree fit needs it to cross between rootings. |
| perturb-move.md | PROPOSED, 2026-09-07; AMENDED 2026-09-07 (slices sized, the benefit stage re-primaried on minimum ESS); SLICE 1 LANDED 2026-09-07 (the kernel at weight zero, ab49f83a); SLICE 2 LANDED 2026-09-07 (perturb-balance.R, 30472110); SLICE 3 RUN 2026-09-07: KILL at w = 1, d = 0.16 (d73fb4e0) | Designs the same-variable cut move ("perturb"), the tree-space candidate the mixing survey pre-registers: a fourth kernel that keeps a node's split variable and displaces its cut by one grid position, with a clipped-window proposal correction, a prior-only detailed-balance gate against the occupancy-truncated CGM prior, a twenty-two-file surface enumeration beside the restored swap move, and a benefit study pre-registered on minimum effective sample size at the He and Hahn cell's shipped four-chain configuration; the kernel shipped at default weight 0, the balance gate landed, and the benefit study's confirmatory run killed the move on its pre-registered cell, closing a nonzero default share while leaving the kernel in the tree at weight zero. |
| nog-gibbs.md | PROPOSED, 2026-09-07; AMENDED 2026-09-07 (the veto's real law and the neighbourhood as a rank stratum, the cost table at 1 - stump%, the cost instrument, the balance gate sized, the surface at twenty-four files); SLICE 1 LANDED 2026-09-07 (the kernel at weight zero, 7fb166ca); SLICE 2 LANDED 2026-09-07 (rule-gibbs-balance.R, d888c9f3); SLICE 3 RUN 2026-09-07: NOT KILLED at d = 0.16, the coverage secondary fails (50032833); CUT-ONLY PILOT 2026-09-08: the private cut-only variant keeps about half the Trig+poly gain and all of the Single index one at 1.04 sweep-equivalents against 2.21 (a6f44e12); coverage flag dissolved by the reference arm 2026-09-08 (e002c10d); DOSE RESPONSE 2026-09-08 (17505c50). ADOPTED POST-RELEASE 2026-09-08: the maintainer adopts the cut-only rule draw at `d` = 0.16, to land after the first release. | Designs `rule_gibbs`, an exact collapsed Gibbs draw of the split rule at a nog node (an interior node whose two children are both leaves), replacing the Metropolis change proposal there: the ancestor-determined neighbourhood of every (available ordinal variable, admissible cut) restricted to the branch-rank stratum the empty-leaf veto's lexicographic law makes current, weighted by the cut scan's rank-admitted marginal times the prior factors that do not cancel, acceptance identically one and no reverse count, with a per-side rank-aware scan so the prior-only balance arm can run, a twenty-four-file surface enumeration beside the landed perturb move, and a benefit study pre-registered on minimum effective sample size at the He and Hahn cell; nothing is built. |
| level-fibre.md | PROPOSED, 2026-09-07; AMENDED 2026-09-07 (the linear leaf out of slice 1 and recorded as a door with its `m n` price, the perturbation algebra halved, the empty-leaf reason restated on the zero pin, the pilot as the residual channel with an advisory bar on medians, the backfit-exact gate repaired by profiling, the control slot fixed at creation, the cost against 16.3's own unit); SLICE 1 LANDED 2026-09-07 (the step behind the flag at default off, cbe80534); SLICE 2 PILOT CONFIRMS 2026-09-07 (bf1a4c9e); SLICE 3 RUN 2026-09-07: KILLED, the primary reads -0.9 and -2.5 against a +8 bar (ca92c11f); AUTO UNDER THE FROZEN MIXTURE 2026-09-08 (127f04ee) | Designs an exact Gibbs draw on the level fibre - add a constant to every leaf of a tree with the constants summing to zero across the forest, leaving the fitted function unchanged and the leaf prior as the whole conditional - deriving the m - 1 dimensional draw and its per-leaf-precision form, settling the leaf-model and family reach (the variance forest's multiplicative fibre has no closed form, GP leaves have no level out of sample, the linear leaf's dense fit slab a door priced at a third of a sweep), placing it at the top of the sweep for under 1/5000 of a sweep with no data pass at the constant leaf, and pre-registering a frozen-structure pilot ahead of C1's summed minimum ESS. Shipped as a tri-state `levelGibbs` whose `NA` default takes the step for a forest exactly where that forest's structural mixture is frozen; no baseline moves. |
| benchmark-surfaces.md | COMPLETE (survey), 2026-09-06 | Surveys the test problems the tree, nonparametric-regression, causal-inference and MCMC literatures use, and proposes a measurement battery: an average-case core of four realistic residual surfaces plus eight pathologies covering distinct failure modes, with the acceptance rule that a change must be neutral-or-better on the core and better on at least one pathology. |
| interaction-constraints.md | LANDED, 2026-07-21 | Adds per-forest interaction constraints (`interactions =`): a max order cap and/or named co-occurrence deny/allow rules. |
| threaded-predict.md | LANDED, 2026-08-25 | Wires `predict()`'s `n.threads` argument to a real thread fan-out over (chain, draw), bitwise identical at every thread count. |
