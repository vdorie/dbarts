# Documentation inventory and disagreements, for a later triage

Two parts. Part 1 is the per-file inventory of everything under the documentation tree, with each file's own status line and the content that exists nowhere else; it is the input to a decision about what to keep, fold or delete. Part 2 collects every place a document, a news entry or a code comment says something the code or the git history contradicts, gathered from all five code-side slices, the downstream sweep and the two evidence sweeps.

## Part 1. Inventory

882,316 words across docs/plans, docs/plans/archive, docs/design,
docs/architecture.md and docs/README.md; 179,652 of that is the
review-2026-08-24 subtree. Status is each file's own live line.

The counts and rows below are as the tree stood at b7266da9 and have not been
recounted since. Eleven files have been added: the six plans of the arcs that
landed after it - pure-c-header.md, front-door.md,
interfaces-and-dependencies.md, engine-performance.md,
memory-footprint-audit.md, per-draw-callbacks.md - plus sparse-formula-audit.md
under docs/plans and this directory's own six files, and, under docs/design,
memory-footprint.md, engine-constants.md, per-draw-callbacks.md and
engine-generics-review.md. Every one is indexed and passes the documentation
freshness check; none is inventoried here.

### docs/plans (active, 35 files)

| plan doc | words | status | unique content not recorded elsewhere |
|---|---|---|---|
| README.md | 2127 | process contract | The only statement of roles, the plan template, the RNG gate classes and which gates each requires, and the CI path-filter map |
| INDEX.md | 5028 | manifest, stamped 849f08ea 2026-09-02 | One-line purpose for all 33 active and 132 archived plans; the cheapest replacement for reading the archive |
| adoption-slate.md | 16382 | LANDED (S1-S8), 2026-08-15 | The eight-slice gate log, the "Doors held open" list (cmp-D17) and the arc-closure anchor-refresh audit; the features themselves are in docs/design/r-c-division.md |
| architecture-numerical-review.md | 1472 | reference, findings only | The two readers' verdict lists, including the accepted-risk register (cancellation bounds, drift magnitudes, floor behaviour) that exists nowhere else |
| bartcore-review-tour.md | 2841 | current at 127f04ee | The merge case itself and its reading order; VD's `rc-gate` names this file |
| bcf-cross-host.md | 5145 | LANDED 3f532af2 | The corrected exemption list, the discrimination probes and their exit codes, and the seeds-axis door (cmp-D04) |
| bcf-latent-evidence.md | 5825 | LANDED gate + derivation + SBC; 3 equivalence scenarios PROPOSED | The latent quadrature derivation, the chain-length finding that keeps both latent arms out of the SBC matrix, and the scope note for cmp-U12 |
| capi-shape.md | 5309 | LANDED 9df0cb50 | Slice log for the header's pre-1.0 shape; the shape itself is in inst/include/dbarts/dbarts.h and docs/design/public-surface.md |
| column-kind-consolidation.md | 18092 | LANDED (7 slices) | The "Four residues" section (cmp-U24, cmp-U10), the denseBorrowed naming-inversion argument, and the narrowing-at-the-boundary rule |
| composition-refusals.md | 7113 | LANDED 936825d7 | Slice log; the refusals themselves are in R/spec.R and man/dbarts.Rd |
| correctness-audit.md | 3973 | reference, all findings fixed | Term-by-term re-derivation of every acceptance ratio and conjugate update - the only such derivation in the repo |
| dbarts-h-freeze.md | 2927 | LANDED 6446ddce | Slice log; the header contract is the header |
| gp-followups.md | 186 | research-open | Nothing beyond the TODO entry, plus one stale blocker claim |
| gpu-bart.md | 332 | DONE (memo complete) | Nothing: the verdict is docs/design/gpu-bart.md |
| group-by-exposure.md | 201 | RETIRED 2026-09-06 | Nothing: superseded by docs/design/retire-grouped-random-effects.md |
| latent-subset-mask.md | 11005 | LANDED 2026-08-13, arc complete | Per-family composition rules for the active-row mask and the errata about the empty-leaf veto; the shipped contract is docs/design/active-rows-mask.md |
| mixing-program-report.md | 5618 | RECORD, 2026-09-08 | The maintainer-facing account of the whole mixing program and the adoption decision (cmp-V01, cmp-U22); written for VD, not derivable from the design docs |
| multiforest-veto-rate-falsifier.md | 11743 | RUN AND REPORTED (YELLOW) | The pre-registered falsifier, its acceptance-rate tables and the YELLOW verdict that authorized the multi-forest mutation arc |
| nameable-calibration.md | 7960 | ARC COMPLETE | Slice log; the shipped surface is docs/design/nameable-calibration.md |
| pre-review-cleanup.md | 2801 | LANDED 7cd71f2d | VD's rulings on the four adversarial pre-review reports - the only place those rulings are recorded |
| predict-surface.md | 7778 | LANDED 78f334c1 | Slice log; the shipped signatures are man/ and R/ |
| prerc-surface-freeze.md | 726 | DECIDED 2026-08-25 | The nine pre-release-candidate rulings and the "Post-1.0 by rule" list (cmp-X03); short and load-bearing |
| python-bindings.md | 325 | research-open, no spike run | Nothing beyond the TODO entry |
| rd-records.md | 6911 | LANDED 52c10e02 | Slice log for documentation corrections; the corrections are in man/ |
| release-candidate-review.md | 47895 | SPECCED 2026-08-17, in execution | The per-commit landing notes for the whole release-candidate program, including gate counts per commit - the closest thing to a change log for the last month of the branch |
| repo-modernization.md | 480 | standing hygiene | Nothing beyond the TODO entry |
| retire-grouped-random-effects.md | 2326 | LANDED 1e5f80b2, two prerequisites open | The two release prerequisites and their numeric bars (cmp-L05, cmp-L06); the argument is docs/design/retire-grouped-random-effects.md |
| sbc-calibration.md | 6408 | DONE, all tiers | The running calibration log and the BCF scale-ridge diagnosis that motivated the a-ridge and b-ridge moves |
| sbc-family-tiers.md | 2988 | BUILT | The measured burn ladders per family and the three open arms (cmp-U09); the ladders exist nowhere else |
| simd-survey.md | 2153 | reference, read-only | The arm64 SIMD candidate survey and the invariant that keeps the sufficient-statistic kernel scalar for bitwise reproducibility |
| sparse-extensions.md | 281 | mixed | Nothing beyond the TODO entry |
| surface-refusals.md | 8183 | LANDED d48aef8a | Slice log; the refusals are in R/ and inst/tinytest |
| tau-slice-review.md | 3313 | RETIRED 2026-09-06 | The interweaving derivation in sec 4(c) - now the shape stan4bart's tau-mixing bar (cmp-L05) would need; the rest describes deleted code |
| weighted-binary.md | 325 | parked memo | Nothing beyond the TODO entry |
| x86-simd-plan.md | 2679 | partly landed, three items open | The x86 dispatch reality on the bench box and the three open items (cmp-U23); companion to simd-survey.md |

### docs/plans/review-2026-08-24 (102 files, 179,652 words)

The second whole-branch review's working directory. Only the prose files are
listed individually; the run scripts and captured output are grouped.

| plan doc | words | status | unique content not recorded elsewhere |
|---|---|---|---|
| consolidated-report.md | 8732 | snapshot at b102e17c | The review's numbered findings and its DEFER bin (cmp-D19); the fixes themselves are recorded in release-candidate-review.md |
| decision-brief.md | 3379 | read-only, 2026-08-24 | The three maintainer decisions the review put up, with the probe numbers behind each |
| gate-ledger.md | 3956 | snapshot at b102e17c | Gate and baseline inventory as of 2026-08-24; superseded by the current MANIFEST and the gates ledger |
| gate-ledger-read.md | 2109 | independent read | Its refutation of the gate ledger, drawn from a 3091-run CI dump - the only place that dump is summarized |
| anchor-main.md | 2179 | read-only | The one cross-implementation anchor against main's 0.9-34 engine; every other equivalence record is bartcore-vs-bartcore |
| calibration-sbc.md | 2073 | measurement | The 11-arm ensemble-scale SBC run behind the family tiers; two of its arms describe deleted grouped code |
| matrix-results.md | 3319 | executable matrix, run | The consistency-matrix outcomes; the script (matrix.R) reproduces them |
| matrix-review-entries.md | 2695 | review 1 | Per-entry findings on fitting entries and their documentation |
| matrix-review-generics.md | 2587 | review 2 | Per-cell findings on fit classes x generics |
| generics-survey.md | 4106 | survey | The generics census behind the phase-2 spec, including the ordinal tail-precision door (cmp-U08) |
| generics-phase2-spec.md | 656 | ruling | The orchestrator's rulings on the survey's ambiguities and the doors it left |
| reading-R-list.md | 3739 | candidate list | Per-symbol removal candidates in R/ with cost and confidence; nothing decided |
| reading-engine-list.md | 3070 | candidate list | The same for the engine, bridge and support libraries |
| mutation-A-findings.md | 2514 | leg A | Mutation results over 65 changed test files |
| mutation-B-findings.md | 3085 | leg B | Mutation results over the C++ component suite |
| mutation-C-findings.md | 2537 | leg C | Mutation results over untouched test files, first half |
| mutation-D-findings.md | 2194 | leg D | Mutation results over untouched test files, second half, including the shared-helper gap (cmp-U25) |
| mutation-B-evidence.md | 1730 | evidence | Long-form ladder audit behind leg B |
| wave3-plan.md | 6389 | implementer spec | The wave-3 fix spec; what landed is in release-candidate-review.md |
| xbart-oracle-memo.md | 1544 | memo | What an oracle for the cross-validation entry would be, and what is already covered |
| leaf-scale-pin-memo.md | 1344 | scoping memo | Prices the creation-time leaf-scale pin that would have unblocked the AFT SBC arm; superseded by the status setter that shipped instead |
| memos/backlog-value-scan-2026-08-24.md | 2467 | scan | Per-backlog-item value verdicts (pre-RC / post-1.0 / gated) |
| memos/backlog-value-scan-critique-2026-08-24.md | 1947 | critique | Refutations of that scan, probe-verified |
| memos/prerc-lens1-surface.md | 3277 | audit | The public-surface count and its freeze findings; the rulings are prerc-surface-freeze.md |
| memos/prerc-lens2-backlog.md | 1919 | audit | The doors-and-gates read that fed prerc-surface-freeze.md; names the consumer-absence gating rule violation |
| memos/prerc-lens3-external.md | 2153 | audit | The CRAN-reviewer / first-time-user / linked-consumer read from a built tarball |
| memos/pre-review-audit-review-chain.md | 5886 | adversarial report | Findings on the review chain itself; VD's rulings are pre-review-cleanup.md |
| memos/pre-review-critic.md | 2342 | adversarial report | Completeness critique; rulings as above |
| memos/pre-review-cruft.md | 12610 | adversarial report | The agent-accumulation and cruft inventory - the largest single catalogue of what could be deleted |
| memos/pre-review-yagni.md | 4161 | adversarial report | The YAGNI inventory; rulings as above |
| memos/monotone-branch-fill-bench.md | 2590 | not adopted | The bench arm behind TODO `monotone-leaf-quadrature` (cmp-U03); the variant is on branch archive/monotone-branch-fill |
| memos/predict-replay-slice-spec.md | 1105 | spec | The out-of-sample per-forest replay slice; landed |
| memos/threaded-predict-memo.md | 3486 | memo r1 | Superseded by r2 |
| memos/threaded-predict-memo-r2.md | 3791 | memo r2 | The design that landed; the shipped record is docs/design/threaded-predict.md |
| memos/threaded-predict-critique.md | 3582 | critique | 35 findings against that memo, 33 accepted |
| memos/tree-store-burnin-memo.md | 1476 | scope memo | The saved-tree store's rotation across successive run calls |
| sbc-logs/summary.txt | 1112 | run log | The 11-arm band table; two arms describe deleted grouped code |
| consol/*.R (3 files) | 2997 | reproduction scripts | Independent reproduction of the review's R-surface claims |
| generics/*.R, *.csv (18 files) | 32988 | census scripts + captured grids | The executed generics census; the findings are in generics-survey.md |
| matrix*.R, matrix-grid.csv (3 files) | 11670 | matrix scripts + grid | Reproduce matrix-results.md |
| mutation-{A,C,D}-evidence/* (26 files) | 9159 | drivers, tables, captured runs | Per-mutant kill/survive records behind the four findings files |
| sbc-logs/*.R (5 files) | 1533 | SBC drivers | Reproduce the calibration run |

### docs/plans/archive (132 files)

Every row's status is the file's own; every one-line purpose is already in
docs/plans/INDEX.md's "Archived" table, so "nothing beyond the INDEX row"
below means exactly that: deleting the file loses the slice-by-slice gate log
and nothing a reader of the INDEX plus the paired design doc would miss.

| plan doc | words | status | unique content not recorded elsewhere |
|---|---|---|---|
| architecture-doc.md | 443 | LANDED | Nothing beyond the INDEX row; the product is docs/architecture.md |
| autoconf-dead-code.md | 518 | LANDED | Nothing beyond the INDEX row; the result is configure.ac |
| bart2-argument-consolidation.md | 16953 | COMMITTED SPEC, all 8 forks VD-decided | The eight recorded VD fork resolutions and the argument-semantics rules that recur (subset forwarding, per-forest packaging, the cross-validation grid override); not restated anywhere |
| bcf-b-ridge.md | 3731 | shelved / unimplemented | The full generalized-inverse-Gaussian derivation for the treatment-scale rescale, plus its implementation traps - the whole basis for cmp-X02 |
| bcf-bartcause-relocation.md | 15125 | LANDED, arc closed | The eight recorded VD decisions on relocating the causal fit function to bartCause, and the dbarts-side guards it needed |
| bcf-public-surface.md | 8455 | LANDED (S0-S6) | Slice log; the shipped surface was renamed by multiforest-extension-surface M2 |
| bcf-ridge-interweaving.md | 4128 | LANDED 9617c94 | The prognostic-ridge derivation and the re-scoped acceptance that spawned cmp-X02 |
| bcf-sigma-residual.md | 1494 | RESOLVED | The burn-transient diagnosis and the extreme-tail door (cmp-D01) |
| bcf-testfits-guard.md | 708 | LANDED | Nothing beyond the INDEX row |
| binary-kforest-prior-default.md | 17184 | ARC COMPLETE | The prior-coverage argument for the family-aware amplitude default and the refutation of the mixing hypothesis; also the recorded refusal that became cmp-D02 |
| block-fusion-stage-a.md | 4923 | SUPERSEDED | Nothing beyond the INDEX row; docs/design/block-fusion.md carries the verdict |
| block-fusion-stage-b.md | 4422 | NO-GO | The measured 4-6x slowdown table; summarized in docs/design/block-fusion.md |
| blocked-jacobi-trees.md | 3003 | KILLED as a build target | The head-to-head against within-chain threading that reversed the earlier GO |
| bridge-error-path-leaks.md | 573 | LANDED | Nothing beyond the INDEX row |
| c-api-growth.md | 7620 | LANDED | The additive-growth rules for the results struct and by-name state blocks; the surviving rule is in TODO's release block |
| capi-callbacks.md | 953 | LANDED | Nothing beyond the INDEX row |
| capi-dispatch-table.md | 3897 | DECIDED (VD 2026-07-16), arc complete | VD's three numbered decisions on the ABI mechanism, and the record of the cross-repo sanitizer job landed then removed (cmp-X04) |
| change-move-fix.md | 2520 | LANDED | Nothing beyond docs/design/change-move-balance.md |
| cheap-uniformity.md | 4526 | LANDED (S0-S4) | The measured sparse-predict win (34.5x faster, 9.96x smaller peak memory) and three defects found along the way |
| chi-default-research.md | 1070 | LANDED | The 48-cell study behind the binary k default change - the only record of that simulation |
| chi-hyperprior-df.md | 517 | LANDED | Nothing beyond the INDEX row |
| chi-k-empty-leaf-count.md | 501 | LANDED | Nothing beyond the INDEX row |
| chi-k-runaway.md | 769 | LANDED 4797bc0 | The recorded reasoning for leaving the state restore uncapped |
| collapse-merge.md | 287 | LANDED | Nothing beyond the INDEX row |
| composition-mixing-probe.md | 21517 | KILLED (harm clause fired) | The registered gate architecture, the eleven adopted critique findings, and the harm-clause tables that withdrew the survey's top-ranked candidate |
| constant-leaf-fits.md | 987 | LANDED | The x86 bench discharge numbers (14-18% win, 22-28% mutation regression) |
| constant-leaf-suffstat.md | 561 | LANDED | Nothing beyond the INDEX row |
| consumer-spec-surface.md | 495 | LANDED | Nothing beyond docs/design/consumer-spec-surface.md |
| convergence-diagnostics.md | 458 | LANDED | Nothing beyond the INDEX row |
| cran-readiness.md | 3134 | CLEAN 2026-07-25 | The submission battery TODO's release block points at |
| csc-code-validation.md | 753 | LANDED df79f17 | Nothing beyond the INDEX row |
| cutpoints-shrink-orphan.md | 332 | LANDED | Nothing beyond the INDEX row |
| data-ownership.md | 331 | COMPLETE | Nothing beyond docs/design/data-ownership.md |
| data-ownership-1-container.md | 1661 | LANDED | Nothing beyond docs/design/data-ownership.md and data-store.md |
| data-ownership-2-ingestion.md | 1756 | LANDED | As above |
| data-ownership-3-mutation.md | 3244 | LANDED | As above |
| data-ownership-4-views.md | 2519 | LANDED | As above |
| data-ownership-5-sparse.md | 4086 | LANDED | As above, plus the program-complete declaration |
| data-review-remediation.md | 791 | LANDED | Nothing beyond the INDEX row |
| data-store-consolidation.md | 1164 | LANDED | Nothing beyond docs/design/data-store.md |
| data-store-residuals.md | 748 | LANDED | Nothing beyond docs/design/data-store.md |
| dbarts-h-reshape.md | 16722 | LANDED 2026-08-13 | The four sister-package migration records and the still-owed treatSens call-site note (four sites must pass 1, not 0) |
| empty-leaf-veto.md | 913 | LANDED (NO-GO on removal) | Nothing beyond docs/design/empty-leaf-veto.md |
| equivalence-ci.md | 425 | LANDED | Nothing beyond the workflow and the MANIFEST |
| facade-shape.md | 666 | LANDED 40082c7 | Nothing beyond the INDEX row |
| family-on-model.md | 398 | LANDED | Nothing beyond the INDEX row |
| flat-format-v2.md | 605 | LANDED | Nothing beyond the INDEX row; the format is documented in the state-format code |
| forest-combiner.md | 2972 | LANDED | Nothing beyond docs/design/forest-combiner.md |
| forest-split-bcf.md | 4250 | LANDED (two phases) | Nothing beyond docs/design/bcf.md |
| fuzz-state-roundtrip.md | 409 | LANDED | Nothing beyond the INDEX row |
| gate-blindspot-audit.md | 1407 | LANDED | The 16-poison sweep results and the feature x gate coverage matrix - the origin of most current gates |
| gate-hardening-1.0.md | 2060 | LANDED | Nothing beyond the gates it built |
| gp-cache-test-flake.md | 446 | LANDED | Nothing beyond the INDEX row |
| grouped-equivalence.md | 226 | CLOSED, LANDED | Nothing: describes deleted grouped scenarios |
| grow-from-root-categorical-scan.md | 10279 | LANDED (S0a-S5) | The pre-registered falsifier confirming the shipped-versus-exact draw-law gap, and the missing-code split asymmetry it found |
| grow-from-root-default-study.md | 2959 | KILLED (measured) | Pre-registration; the full data is docs/design/grow-from-root-default.md |
| grow-from-root-warm-start.md | 901 | LANDED | Nothing beyond docs/design/grow-from-root.md |
| grow-from-root.md | 552 | LANDED (memo phase) | Nothing beyond docs/design/grow-from-root.md |
| heteroscedastic.md | 1025 | LANDED (C1-C4) | Nothing beyond docs/design/heteroscedastic.md |
| hot-layer-u8.md | 592 | NO-GO (phase 2) | The measured absence of an arm64 partition-throughput win at 8-bit widths |
| hurdle.md | 1405 | LANDED | Nothing beyond docs/design/hurdle.md |
| interaction-constraints.md | 837 | LANDED f455d7c | Nothing beyond docs/design/interaction-constraints.md |
| interaction-constraints-p4.md | 518 | LANDED | Records that soft path penalties and formal heredity stay deferred (cmp-V07) |
| interface-review.md | 4037 | LANDED | The review-2 audit of the exported R surface and its 11 taste calls |
| kernel-cleanups.md | 327 | LANDED | Nothing beyond the INDEX row |
| linear-leaf-reuse.md | 701 | LANDED | The measured 18-35% linear-leaf speedup |
| monotone-bart.md | 1000 | LANDED | The bias an intermediate design carried and how the exact gate caught it |
| monotone-prior-draw.md | 634 | LANDED 173a710 | Records that this unblocks the monotone SBC arm (cmp-U09) |
| moves-degenerate-root-guard.md | 477 | LANDED | Nothing beyond the INDEX row |
| multi-forest-models.md | 399 | LANDED (tracker) | Nothing beyond the INDEX row |
| multiforest-extension-surface.md | 50694 | ARC COMPLETE (M0-M4.5) | The largest file in the repo: the knob map from the old causal argument names to `forests = list(forest(...))`, the refusal-set census, and the recorded VD fork resolutions. The shipped design is docs/design/multiplier-combiner.md |
| multiforest-mutation-gaps.md | 1339 | LANDED (4 commits) | The two doors that became TODO `multiforest-mutation-gaps` (cmp-D08) |
| multiforest-predictor-mutation.md | 11922 | ARC COMPLETE (SL, S0-S4) | The per-slice stranding hazards and the per-forest column mask as the opt-out |
| multinomial.md | 6037 | ARC CLOSED (C1-C7) | Nothing beyond docs/design/multinomial.md |
| multinomial-counts.md | 2849 | LANDED (C1-C3) | Nothing beyond docs/design/multinomial.md |
| multinomial-counts-mutation.md | 8622 | ARC COMPLETE (S1-S5) | The reopen door that competing risks depends on (cmp-V05) |
| multinomial-formula.md | 1903 | LANDED (C1-C2) | The count-matrix formula pattern TODO `multinomial-doors` says to adopt (cmp-D09) |
| multinomial-level-centering.md | 1421 | LANDED ec2a3d0 | The exact leaf-space conditional that closed the calibration flag |
| multinomial-margins.md | 773 | LANDED | Nothing beyond the INDEX row |
| multinomial-test-surface.md | 3374 | LANDED (C1-C2) | Nothing beyond docs/design/multinomial.md |
| multinomial-varcounts.md | 1518 | LANDED (C1) | Nothing beyond the INDEX row |
| mutation-fuzzing.md | 474 | LANDED | Nothing beyond the INDEX row |
| mutation-journal.md | 611 | LANDED | Nothing beyond the INDEX row |
| negative-binomial.md | 933 | LANDED (C1-C3) | Nothing beyond docs/design/negative-binomial.md |
| ordinal-outcomes.md | 759 | LANDED | Nothing beyond docs/design/ordinal.md |
| package-review-remediation.md | 2535 | LANDED | The seven-reviewer panel's tier-1 through tier-5 findings; tier 5 is the only record of the three untracked performance items (cmp-D18) |
| parallel-falsifiers.md | 948 | LANDED | The three measured falsifier outcomes behind docs/design/parallel-bart-frontier.md sec 5 |
| pointwise-loglik.md | 941 | LANDED | Nothing beyond the INDEX row |
| post-mutation-assertions.md | 341 | LANDED | Nothing beyond the INDEX row |
| ppd-sigma-pairing.md | 1239 | LANDED | Nothing beyond the INDEX row |
| pre-release-surface-fixes.md | 663 | LANDED | Nothing beyond the INDEX row |
| prior-constants.md | 424 | LANDED | Nothing beyond docs/design/prior-defaults.md |
| prior-predictive.md | 1511 | LANDED | Nothing beyond the INDEX row |
| r-ingestion-cleanups.md | 368 | LANDED | Nothing beyond the INDEX row |
| r5-cleanup.md | 382 | LANDED | The only record of why the thread start/stop methods were removed; inst/NEWS.Rd does not mention it |
| range-scaling.md | 779 | LANDED | The decision memo keeping range scaling over standardization, and the origin of the `updateScale` flag |
| rbart-custom-prior-divergence.md | 767 | LANDED | Nothing: describes deleted grouped code |
| rbart-loop-fix.md | 682 | LANDED | Nothing: describes deleted grouped code |
| rbart-loop-profile.md | 915 | LANDED (measurement) | Nothing: describes deleted grouped code |
| rchk-ci.md | 604 | LANDED | Nothing beyond the workflow |
| readability-review.md | 612 | LANDED | Nothing beyond the INDEX row |
| retrospective-reviews.md | 970 | LANDED (program complete) | The six-review status log with per-review totals (cmp-K17) |
| revdep-smoke-ci.md | 349 | LANDED | Nothing beyond the workflow |
| review-perf-followups.md | 1527 | ARC CLOSED | Which tier-5 items landed and which did not (cmp-D18) |
| roadmap-survey.md | 934 | LANDED | The 15-year ecosystem ranking; superseded in breadth by docs/design/bart-landscape.md |
| robust-errors.md | 835 | ARC CLOSED | Nothing beyond docs/design/robust-errors.md |
| runsbcbcf-repair.md | 2766 | LANDED 62caed0 | The "setData door survey" TODO `multiforest-mutation-gaps` cites (cmp-D08) |
| sampletreesfromprior-midchain.md | 655 | LANDED 1947b10 | Nothing beyond the INDEX row |
| sbc-ci-gate.md | 864 | LANDED | Nothing beyond the workflow |
| setpredictor-leafof-rebuild.md | 1816 | CLOSED 2026-08-05 | The measured ceiling (+9.7%) that declined three mechanisms - the reason not to retry them |
| sigma-df-zero-weights.md | 530 | LANDED | Nothing beyond the INDEX row |
| simd-flag-multiversioning.md | 1745 | NO-GO | The measurement that closed flag-built wide-instruction variants |
| small-validation-fixes.md | 1255 | LANDED | Nothing beyond the INDEX row |
| snapshot-tests.md | 456 | LANDED | The rationale for per-family snapshot files, which CLAUDE.local.md's regeneration note depends on |
| state-continuation.md | 476 | LANDED | Nothing beyond the INDEX row |
| state-format-policy.md | 390 | LANDED | Nothing beyond the INDEX row |
| survival-grouped-surface.md | 927 | LANDED | Nothing: the grouped half is deleted |
| survival-models.md | 933 | LANDED except grouped surface | The follow-up doors TODO `survival-followups` cites (cmp-V05) |
| tau-cauchy-exact-ig.md | 924 | LANDED | Nothing: describes deleted grouped code |
| tau-slice-stepout-cap.md | 524 | LANDED | Nothing: describes deleted grouped code |
| test-data-parity.md | 3970 | LANDED (CLOSED) | The measured test-side memory shrink (1.83x-6.98x) |
| test-fit-parallel.md | 356 | LANDED | Nothing beyond the INDEX row |
| test-suite-trim.md | 433 | LANDED | Nothing beyond the INDEX row |
| tests-cpp-split.md | 330 | LANDED | Nothing beyond the INDEX row |
| typed-ingestion.md | 5019 | LANDED (slices 1, 2a, 2b) | The recorded typed-ingestion doors TODO `sparse-extensions` cites (cmp-V04) |
| variance-forest-mutation-routing.md | 4896 | LANDED (S1-S5) | The two doors held open, including the scale-leaf staleness (cmp-D12) and the unverified rescale factor |
| vignette-refresh.md | 316 | LANDED | Nothing beyond the vignettes |
| warm-starts.md | 926 | LANDED | Nothing beyond the INDEX row |
| weighted-binary-ppd.md | 788 | LANDED | Nothing beyond the INDEX row |
| within-chain-threading.md | 803 | NO-GO CLOSED | Nothing beyond docs/design/within-chain-threading.md sec 8 |
| x86-simd.md | 534 | CLOSED/SUPERSEDED | Nothing beyond the INDEX row and x86-simd-plan.md |
| zero-weight-exactness.md | 6690 | ARC COMPLETE (S0-S3) | The exact-zero snap rule and the caller-settable per-forest weight the mask arc later composed with |

### docs/design (55 files) plus the two orientation docs

| plan doc | words | status | unique content not recorded elsewhere |
|---|---|---|---|
| INDEX.md | 2741 | manifest | One-line purpose and status for all 55 design docs |
| active-rows-mask.md | 1668 | reference, current | The shipped per-observation mask contract, per family; the only place it is stated |
| aft-status-setter.md | 3334 | slices 1-2 LANDED; 3-4 PROPOSED | The slice list (cmp-U11, cmp-U11) and the SBC admission measurement; status line is stale, see Disagreements |
| aft-variance-forest.md | 2466 | LANDED 2026-09-06 | The argument that the latent-channel reason for refusing a variance forest does not hold for the survival family |
| bart-as-a-component.md | 2756 | LANDED 2026-08-19 | The internal embedding contract: which mutations are legal between sweeps and what state a mutation does not carry. The user-facing twin is vignettes/dbarts-as-a-component.Rmd |
| bart-landscape.md | 5852 | snapshot 2026-08-12 | A survey of 35 BART implementations; purely descriptive, no in-repo dependency |
| bcf.md | 4020 | LANDED | The causal-forest model, its calibration and its exact-posterior gate |
| benchmark-surfaces.md | 23639 | COMPLETE (survey), 2026-09-06 | The measurement battery's cells, their generators and every recorded arm result; the mixing report is its summary |
| block-fusion.md | 7095 | CLOSED, WONT-DO | The full engineering design plus the measurement that killed it |
| change-move-balance.md | 1108 | LANDED 2026-07-08 | The detailed-balance defect and its repair - a posterior-moving correction against 0.9-x |
| consumer-spec-surface.md | 964 | LANDED | The specification-resolution contract for linked consumers |
| core-generalization.md | 6802 | LANDED, phase 6 open | The founding engine design; its phase 6 is cmp-X01. States an R >= 4.3 floor the shipped DESCRIPTION contradicts (build ledger) |
| correlated-outcomes.md | 1425 | RESOLVED 2026-07-22 | Why the multivariate case needs no engine change and the serial-correlation case does (cmp-D03) |
| data-layout.md | 3617 | CLOSED - SHELVED | The per-node contiguous layout design and the re-measurement (~10%) that shelved it |
| data-ownership.md | 3024 | COMPLETE | The owned-quantized-codes decision and what it rejected |
| data-store.md | 4085 | REFERENCE | The normative predictor-store invariants; required reading before data work |
| empty-leaf-veto.md | 4083 | DECIDED (keep-and-document) | Why the veto stays rather than becoming occupancy-aware; also the source of cmp-D13 |
| error-style.md | 7492 | ADOPTED for new messages | The message-style rule and its evidence base (base, stats, Matrix, survival, lme4, mgcv), plus VD's refinement that tidyverse practice carries no authority |
| feature-matrix.md | 2704 | LIVING REFERENCE | The per-model capability matrix, cite-checked by the freshness guard; the fastest answer to "can family X do Y" |
| forest-combiner.md | 3856 | LANDED | The combiner hierarchy's shape |
| forest-ranef-interweaving.md | 3761 | SUPERSEDED 2026-09-06 | Nothing live: the door it held closed with the grouped path's deletion |
| gp-leaves.md | 3939 | Part 1 LANDED; Part 2 unscheduled | The two-axis separation and the unbuilt non-conjugate move strategy (cmp-X01); also the re-krig nugget note |
| gpu-bart.md | 2559 | NO-GO (survey) | The seven-direction survey and the conditions that would reopen it |
| grouped-random-effects.md | 1797 | RETIRED 2026-09-06 | Nothing live: historical record of deleted code |
| grow-from-root-default.md | 9135 | KILLED (measured) | The full pre-registered study data; the plan file carries only the short form |
| grow-from-root.md | 3962 | MIXED | Why root-down construction ships only as a warm start and not as a sampler |
| heteroscedastic.md | 7174 | LANDED 2026-07-20 | The variance-forest model and its mutation routing |
| hurdle.md | 4537 | LANDED 2026-07-20 | The conditional-independence finding that makes a two-part fit two ordinary fits |
| interaction-constraints.md | 2503 | LANDED 2026-07-21 | The constraint design and the two must-fix critique findings folded in |
| kernel-vocabulary.md | 1179 | REFERENCE | The normative kernel contract; cites a configure flag that no longer exists |
| level-fibre.md | 8859 | slices 1-3 landed, slice 4 killed; auto default 2026-09-08 | The leaf-shift derivation, the frozen-structure pilot and the kill; the shipped default is cmp-K24/K63 |
| linear-leaves.md | 2138 | LANDED 2026-07-04 | The linear leaf model |
| memory-wall-frontier.md | 5563 | CLOSED as idea map | The ranked lever map, the measured fused-pass result (cmp-U04) and which branches closed |
| mia-missingness.md | 1214 | LANDED 2026-07-04 | The missing-direction split rule and the deltas found while landing |
| model-space-survey.md | 6776 | COMPLETE (survey) | The model-class survey behind the multi-forest doors (cmp-D02, cmp-D08) |
| monotone.md | 5562 | LANDED 2026-07-19 | The constrained leaf model; also the record that a monotone constraint silently rewrites the proposal mixture |
| multinomial-mutation-arc.md | 10655 | LANDED 2026-08-24 | The pre-arc inventory and fork prices; sections 1-4 describe code that no longer exists |
| multinomial.md | 3625 | LANDED 2026-07-15 | The softmax model and its level-centering move |
| multiplier-combiner.md | 6017 | LANDED 2026-08-13/14 | The general K-forest basis/amplitude family and its calibration map; records the treatment-ridge door as shut |
| nameable-calibration.md | 2059 | ARC COMPLETE | The named-calibration surface in response units |
| negative-binomial.md | 6741 | LANDED 2026-07-18 | The Polya-Gamma design and the exact-versus-approximate argument gating cmp-V02 |
| nog-gibbs.md | 12469 | PROPOSED; slices 1-3 landed | The exact rule draw's correctness argument, the neighbourhood as a rank stratum, and the cut-only variant's pricing (cmp-V01, cmp-U22) |
| ordinal.md | 5002 | LANDED 2026-07-18 | The cumulative-probit design and its cutpoint block |
| parallel-bart-frontier.md | 3569 | MIXED (survey) | The frontier ranking; three of its candidates are now closed elsewhere |
| perturb-move.md | 10066 | PROPOSED; slices 1-3 landed, slice 4 killed | The cut-move design, its census pricing and the kill (cmp-U05, cmp-K26) |
| pooled-masks.md | 1553 | LANDED 2026-07-04 | The 65535-level categorical design and its inline/pooled boundary |
| prior-defaults.md | 943 | REFERENCE | Every shipped default and its source; the fastest check against a user expectation |
| public-surface.md | 4182 | MIXED | The major-version surface decisions recorded inline as DECIDED; three of its claims are contradicted by the shipped header (capi ledger) |
| r-c-division.md | 3865 | ACCEPTED (VD 2026-08-11) | The standing R-versus-C++ rule in VD's own amended words, and the census that priced it |
| reduced-precision-storage.md | 3854 | LANDED / COMPLETE | The narrowed-storage design, its validation, and VD's build directive |
| retire-grouped-random-effects.md | 4474 | LANDED 1e5f80b2; two prerequisites | The whole argument for the removal, the measured tau-mixing gap, and the two release prerequisites with their numeric bars (cmp-L05, cmp-L06) |
| robust-errors.md | 2044 | LANDED 2026-07-17 | The Student-t augmentation |
| sparse-columns.md | 3898 | LANDED 2026-07-04 | The representation study and the densification threshold |
| survival.md | 6174 | LANDED | Both survival families and the follow-up doors (cmp-V05) |
| swap-removal.md | 4607 | LANDED 2026-09-07; AMENDED | The swap census, VD's option-A choice, and the partial reversal; sections 2-8 describe symbols that exist again (rapi ledger) |
| threaded-predict.md | 2470 | LANDED 2026-08-25 | The threading design and what shipped differently from the proposal |
| tree-mixing-proposals.md | 40731 | COMPLETE (survey) + 5 addenda | The stickiness diagnosis, the move census, the two brainstorm rounds and every candidate's ranking; the mixing report is its executive summary and is 7x shorter |
| weighted-logistic.md | 1285 | LANDED 2026-07-05 | Why logistic weights are tractable and probit weights are not (cmp-V03) |
| within-chain-threading.md | 4700 | CLOSED, NO-GO | The full threading design, its Amdahl model and the measurement that killed it on both architectures |
| ../architecture.md | (orientation) | current-state doc | The engine's layering and mechanism map; the entry point CLAUDE.local.md points at |
| ../README.md | 632 | wayfinding | The only statement of the status-header conventions and the known name collisions between plan and design files |

## Part 2. Disagreements with the documentation

Each bullet gives the path and one line. Grouped by where the drift sits.

### The news file and the shipped man pages

- inst/NEWS.Rd, the 1.0-0 section: no entry records the removal of the sampler's startThreads and stopThreads methods, both public, both documented in main's sampler man page, both announced in the 0.9-31 section.
- inst/NEWS.Rd, upgrading notes: a transactional updatePredictor reference-class method is described that exists on neither branch; the only such name in the tree is a C entry point reached from the bridge.
- inst/NEWS.Rd, the 1.0-0 section: the abbreviation for Bayesian causal forests is used four times before it is introduced anywhere; the phrase is spelled out only in later and unconnected entries.
- inst/NEWS.Rd, upgrading notes: the renaming of the run method's thread-count argument reads as a change of contract, but main's own sampler man page already documented the new spelling, so this is the code catching up to a wrong document, and 0.9-x callers who followed the man page were already broken.
- main's man/dbartsSampler-class.Rd documented a scale-update argument on three test-side setters and documented two plotting entries with formals none of which matched main's code; the new page corrects all five and the news file records none of the corrections.
- inst/NEWS.Rd describes the sum-of-squared-residuals fix two ways, as a division by the data range in one entry and by its square in another; the same fact, two descriptions.
- inst/NEWS.Rd, upgrading notes: the swap proposal is described as keeping its place in the kernel with its mass moving to the birth-death move. Both hold at the tip, but one commit deleted the move outright and a later one restored it, so a reader cannot tell that the equivalence baseline named for the first commit was recorded from a build with no swap move at all.
- inst/NEWS.Rd names the deleted C++ headers accurately but does not say that the tree-counting entry and the two per-observation predictor updates have no flat-C replacement at all; a consumer reading only the news file would expect a port, not a gap.
- R/diagnostics.R and man/summary.bart.Rd keep a random-effect scale in the default variable list of the summary method and of the four draws-conversion method pairs, and the man page itself says no shipped family carries it.

### The design documents

- docs/design/public-surface.md names VD nowhere, in any spelling; every paragraph it marks as decided is agent-made. It is the cited record for most of the flat-C decisions and for the factor and ingestion defaults.
- docs/design/public-surface.md still specifies a single packed API version constant and its accessor; neither exists, both having been replaced by the major and minor pair plus the hash accessor. A later paragraph in the same section says so, and the earlier one was never amended.
- docs/design/public-surface.md reports the audit found 19 callables; main's registration table exposed about 70, and 19 is one consumer's usage, not the surface. The same paragraph says the response setter is looked up but never called, which is true of that consumer and false of the other, whose ported branch does call it.
- docs/design/public-surface.md says the version-one header needs no engine additions; its own later paragraph records two engine and state additions the port did need.
- docs/design/public-surface.md says the cut-point integer width build option stays; it was removed.
- docs/design/core-generalization.md states an R 4.3 toolchain floor; the shipped DESCRIPTION states 4.2.0, deliberately, per the commit that set it.
- docs/design/prior-defaults.md calls itself current and still documents the binary node hyperprior default as the old value; the code ships the new one, and two other design docs carry the current value.
- docs/design/retire-grouped-random-effects.md says the bartCause grouped route remains a release prerequisite and that the reverse-dependency smoke test's bartCause leg is an accepted red; nothing in the news file or DESCRIPTION signals to a user that the one known consumer has no working path at this tip.
- docs/design/swap-removal.md is titled and written as a removal, and its sections 2 to 8 describe deleted symbols that exist again; only its section 9 records the reversal.
- docs/design/aft-status-setter.md says the survival arm's matrix admission and its slices 3 and 4 remain proposed; the admission landed at 2c766437, and the same file's own body 250 lines lower says the arm is in the matrix and cites the workflow row, which exists.
- docs/design/kernel-vocabulary.md cites a configure flag that no longer exists.
- docs/plans/archive/bcf-b-ridge.md is indexed as a no-go while its own later note says the treatment-scale ridge is not implemented and its section 7 recommends landing the move later; a no-go and a land-it-later are not the same verdict, and neither is in the backlog.
- docs/plans/archive/bcf-ridge-interweaving.md files a follow-up under a backlog name that is now resolved by a different change, so the half it actually pointed at was never picked up.

### The backlog and the plan process documents

- TODO's Python-binding entry opens by asserting the engine is R-free below the bridge. It is not: the model header calls three Rmath density functions and the sampler and chain headers reach R's print and error entry points through the external IO unit. The C++ test makefile says the same thing about the RNG object and links R to satisfy it.
- TODO's tree-mixing entry says the nog-node balance script is owed; it landed at d888c9f3, an hour before the backlog text denying it was written at 5f93be0b. Only the equal-cost arm is still owed.
- TODO's tree-mixing entry predates the 2026-09-08 addenda to the design doc: the backlog says the coverage regression blocks a nonzero default share and the arm is owed, while the design doc says the coverage flag dissolved, the arm ran, and the kernel was adopted for the next release.
- docs/plans/weighted-binary.md still carries a mid-July statement that no item is deferred to post-release, which a later backlog stamp reverses.
- docs/plans/gp-followups.md says the item is blocked on a state-format plan that is landed, and the format has since moved on a version.
- docs/plans/adoption-slate.md lists two tickets among its residue that no longer exist; both were deleted with recorded closures and the residue list was never updated.
- docs/plans/prerc-surface-freeze.md closes with a post-1.0-by-rule list of seven additive items, none of which appears in the backlog, so the live backlog does not track them; one is moot since the grouped path was deleted.
- docs/plans/archive/package-review-remediation.md says its tier-5 performance items were recorded as backlog entries; none of the three names is in the current backlog.
- docs/plans/INDEX.md's counts are correct at the tip but its status stamp is a week old, and seven plan status lines have moved since.
- docs/plans/bartcore-review-tour.md lists bartCause's mandatory source edits as none, on the ground that it uses the R API only. That is true of the header and of C symbols, but the compat branch carries 22 commits of mandatory R edits, and the edits main still needs are real: the grouped fit function at three sites, the control argument in the optimizer, and the chain count.
- docs/plans/README.md's CI section says a documentation-only push fires the freshness gate alone and that this is the whole gate. True, and it is why the branch tip has been through no engine gate; the last commit that was is 7ffe9032 for the exact gates and 6dfa3524 for everything else.

### The gate records

- docs/plans/review-2026-08-24/sbc-logs/summary.txt records 11 calibration arms including two grouped ones; the grouped path was deleted at 1e5f80b2 and the workflow's matrix now has seven arms, none grouped. The summary's grouped flag describes code that no longer exists.
- benchmarks/baselines/MANIFEST marks the three equivalence baselines recorded at fbff1989 current, but the engine tip is nine source commits later, and the manifest carries no row saying they were revalidated at the tip.
- .github/workflows/sanitizers.yaml states its per-file floors are set at 60 to 65 percent of measured; the total floor of 5200 is 65 percent of the 8037 the same job reports, so the guard tolerates deleting a third of the suite.
- benchmarks/R/mutation-battery.R's header and docs/plans/release-candidate-review.md both describe an inventory of 23 entries; the file now holds 25.
- The build slice reported that none of the five schedule-and-dispatch workflows has ever executed. The gates slice's run history shows each of the five has run at least once. Verified: all five carry, in addition to a schedule and a dispatch trigger, a push trigger on this branch limited to their own workflow file, and that is how each ran. The accurate statement is that no schedule or dispatch trigger has ever fired, because GitHub binds both to the default branch.

### Code comments and the project instructions

- CLAUDE.local.md says a signature change in the shipped header breaks a linked consumer silently. Verified false in both directions: main registered 71 callable names and this branch registers 48, with zero overlap, so a stale consumer binary fails loudly at lookup with R's function-not-provided error; and a rebuilt consumer fails at compile, because the generated stubs re-derive every signature from the header. The residual silent case is narrow and named in the header itself, a same-width in-place type swap under an unchanged field name, and it is reachable only after the exact-ABI opt-in is dropped from the consumers.
- CLAUDE.local.md's layout section describes the two support libraries as linear algebra and threads, and RNG and IO. Both now carry IO: the first gained a new unit holding host-injected print pointers, while the second keeps the R-bound print entry.
- CLAUDE.local.md's gotchas list a generated configuration header among those the build consumes; the 1377-line flat C API implementation does not include it, and only three other translation units do.
- CLAUDE.local.md says the equivalence comparison should expect identical draws while sampling is untouched. In CI that is never what happens: the workflow runs the statistical mode on one architecture against baselines recorded on another and reports a maximum deviation, not identical streams. Bitwise is local only.
- configure.ac's comment on the cut-point integer width is accurate, but the template it names is unchanged from main, so that file still documents itself as configurable while configure hard-codes the one value.
- src/Makevars.win's comment says the Windows ARM64 architecture string is unverified because no such R build was available to probe; the next relevant commit added a CI job that runs natively on that platform and asserts the string, and the comment was not updated.
- The claim that a pre-enforcement consumer binary fails silently rather than loudly does not hold for either published consumer: both predate the handshake entirely and both fail loudly through lookup on names that no longer exist. Where it surfaces is the consumer's choice: at load for the one whose lookups sit in its package init, at the first analysis call for the one whose lookups sit inside the entry point.
- The engine slice records the causal treatment forest's likelihood-invariant ridge as an implemented move shipped disabled; the archive records the treatment-scale ridge as never implemented. Verified: both are true of different objects. The generic per-forest ridge is implemented in the combiner and applies to any forest whose flag is set; the amplitude specification sets that forest's flag false, and the bridge derives each forest's flag from its amplitude prior scale, so the specification's two flags are unreachable from R. The separate generalized-inverse-Gaussian joint rescale of the archived plan was never built.
- The review claim that embedding consumers called the sampler's thread start and stop methods is unsupported. Verified: on main those methods drove the hierarchical thread manager's lifetime through two C entry points, so they were not no-ops there; on this branch the manager is deleted, so they had become no-ops before removal, and the cleanup plan's revdep sweep found no consumer calling them. The removal is real and the news file does not record it.
- bairrtt's readme and CI give the reason for its dbarts branch pin as an entry point being new in 1.0-0; that entry point is exported from main's namespace and announced in main's news file under 0.9-34. What forces the branch install is bairrtt's own version floor.
- TODO's release block says to push dbarts and stan4bart and submit both to CRAN together. CRAN reviews submissions one at a time; together is not available, and the gap between the two acceptances is the unavoidable window in the merge-order table.
- dbarts's README.md tells users to install from GitHub, which serves the default branch. Today that is 0.9-34, the same as CRAN; at the merge the same line silently starts serving 1.0-0. The readme already describes the removed C++ ABI in the past tense while main still ships it.
