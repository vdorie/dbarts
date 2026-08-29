# Gate and baseline ledger - bartcore @ b102e17c, 2026-08-24

Snapshot at b102e17c, 2026-08-24; counts below are as of that commit.

At HEAD (8e8a63ad): .github/workflows/ now holds 12 files - doc-freshness.yaml
was added after this snapshot and is not counted in section 1 below.
benchmarks/baselines/ has been pruned to 11 files (the MANIFEST's row count
no longer matches file count 1:1; section 3's "64 rows, 64 files" is stale).
The central claim - that equivalence.yaml, sbc.yaml, rchk.yaml, valgrind.yaml
and revdep-smoke.yaml have never run in CI - still holds: `gh run list
--workflow=<file> --limit 2` returns "workflow not found on the default
branch" for all five at HEAD, meaning GitHub has never registered a run for
any of them.

Read-only. Every count below was produced in this pass by direct commands (git, gh api/run list, grep,
wc) against the checkout and the GitHub API for vdorie/dbarts, not copied from
scratch/review-lenses-memo.md; where a number matches that memo it is noted, where it differs or adds
precision that is called out. main's only workflow is check-standard.yaml (git ls-tree main --
.github/workflows); bartcore carries 11 files plus one now-deleted one (abi-contract.yaml, footnote).

## 1. Workflows in .github/workflows/ (11 files, all only on bartcore)

Push-triggered on bartcore (6) - fire on every push; PR trigger targets main/master only (sanitizers:
bartcore+main) and has NEVER fired once: across the repo's full history (3091 runs, `gh api
.../actions/runs --paginate`), every run's event is "push" or "dynamic" (pages deploy) - zero
pull_request, zero schedule, zero workflow_dispatch anywhere. All 6 last ran at the tip push
2026-08-24T22:07 UTC, success.

- check-standard.yaml (R-CMD-check): 5-way OS/R matrix (macos-release, windows-release,
  ubuntu-devel/release/oldrel-1); runs tests/tinytest.R (the whole inst/tinytest/ suite) inside R CMD
  check, plus a standalone re-run of test-simd.R. History: 488 success/25 failure/6 cancelled.
- cpp-tests.yaml: `R CMD INSTALL --preclean .` then builds+runs tests/cpp/test_bartcore on macos-latest.
- exact-gates.yaml: runs 20 benchmarks/R gate scripts (Sec.4) in "quick" mode. History: 200 runs, ALL
  event=push - workflow_dispatch has NEVER been invoked, so mode=full (the long grid) has never run in
  CI, ever.
- lint.yaml: `lintr::lint_package()` + `air format --check .`.
- pkgdown.yaml: builds+deploys the docs site; not a correctness gate.
- sanitizers.yaml: ASAN+UBSAN, clang-asan and gcc-asan r-hub containers, over the tinytest suite;
  `ASAN_OPTIONS=detect_leaks=0` (leak detection explicitly OFF - valgrind's job).

Schedule/dispatch-only (5) - NEVER run, not once, on any branch, by any trigger: absent entirely from
`gh api .../actions/workflows` (GitHub only registers a workflow once it has a run or lives on the
default branch; these do neither) and absent from the full 3091-run history query above. Each file's
own header states GitHub binds `schedule` to the DEFAULT branch (main), which does not carry these
files, so the cron cannot fire; each is workflow_dispatch on bartcore only, and dispatch has never been
used either.

- equivalence.yaml: cron Mon 06:00 UTC. Runs equivalence.R (--strict-coverage), bcf-equivalence.R,
  multinomial-equivalence.R against the 3 "current" MANIFEST baselines (Sec.3).
- sbc.yaml: cron Wed 03:00 UTC. 6-arm matrix (Sec.5).
- rchk.yaml: cron Mon 07:00 UTC. kalibera/rchk PROTECT-balance analyzer over the built tarball; fails
  on any [PB]/[UP]/"Suspicious call" line.
- valgrind.yaml: cron daily 07:00 UTC, comment admits "does not fire until this workflow lands on
  main". r-hub valgrind container, memcheck over the tinytest suite.
- revdep-smoke.yaml: cron monthly (1st, 08:00 UTC). 3-package matrix (Sec.7).

Footnote: abi-contract.yaml existed twice (added 9fe39856/e98d94f7, dropped 1f8a360e/99b356d8, absent
now); its one run FAILED on bartcore push 2026-07-16, then it was removed. Not counted in the 11.

## 2. benchmarks/R (38 files) and tools/ (4 scripts) - CI wiring and orphans

Wired (24 of 38, by grep of exec filename against .github/workflows/*.yaml): 20 to exact-gates.yaml
(Sec.4), 3 to equivalence.yaml (equivalence.R, bcf-equivalence.R, multinomial-equivalence.R), 1 to
sbc.yaml (sbc.R).

Orphan (14 of 38 - referenced by no workflow; matches the memo's count, but MEMBERSHIP moved since the
2026-08-17 audit at .claude/rc-review-orphan-legs-2026-08-17.md: hazard-reduction.R/hurdle-reduction.R
were orphans there and are now wired into exact-gates.yaml; backfit-exact.R, composition-matrix.R,
geweke-mc.R, mutation-battery.R are new files since then and are orphans now):

- bench-sampler.R - NOT a true orphan: documented, deliberately-manual speed gate (CLAUDE.local.md,
  docs/architecture.md, 9 docs/plans files); excluded from CI by design ("never run concurrently with
  other load"). Current baseline bench-sampler-ab1dc52.csv is dated 2026-08-08, 16 days stale at tip.
- backfit-exact.R, geweke-mc.R - genuine correctness oracles (gaussian backfit conditional-exactness;
  Geweke 2004 joint-distribution MC test), each cited as the ORACLE for a specific MANIFEST re-record,
  each run exactly once, manually, 2026-08-17/18, never in CI.
- composition-matrix.R - executable check of docs/design/feature-matrix.md claims across 180 (family x
  capability) cells; run once manually (P5, 7997e50c, 2026-08-18: 175 confirmed/0 disagreements/5
  resolved), never CI.
- mutation-battery.R, grouped-mixing.R - self-described "Repeatable" / "Permanent... a fixed gate to
  score against" in their own headers, never wired to CI. The 2026-08-17 audit already found
  grouped-mixing.R's own re-measured numbers (IACT ratio 3.3x/1.9x) diverge from the header's own cited
  prior figures (~25x/~3x) - undetected because nothing re-runs it.
- bartcore-shim.R - not a harness: a shared helper source()'d by equivalence.R (no top-level statement
  of its own); miscategorized as an orphan by naive basename grep.
- change-fix-instrumentation.R, change-fix-stage2.R, parallel-falsifiers.R - BROKEN, confirmed today:
  each depends on src/-side instrumentation (DBARTS_CHANGE_LOG, DBARTS_CHANGE_STATS, BC_FALSIFIER_*)
  absent from src/ right now (`grep -rl <VAR> src/` = 0 hits for all three) - reverted before commit or
  never landed. Cannot produce output at all.
- forest-ranef-collapse-proto.R, negbin-r-update-mixing.R, ordinal-cutpoint-mixing.R - self-contained
  design-research prototypes that already informed a shipped decision (two ended in NO-GO doors); not
  meant to be regression gates.
- linear-leaf-xcheck.R - true orphan: referenced nowhere in the tree outside itself, not even in a
  landing note.

tools/ (4 scripts, none referenced by any workflow, none in `grep tools/ .github/workflows/*.yaml`):
check-build-freshness.R (stale-install guard, inherently local-only - CI always installs fresh),
check-doc-freshness.R (docs INDEX.md completeness) and check-win-drift.R (*.win vs *.in/DESCRIPTION
drift) - both legitimate CI candidates not wired - regenerate-snapshots.R (a writer, not a gate, Sec.8).

## 3. benchmarks/baselines/ via MANIFEST (64 rows, 64 files, exact 1:1)

Role counts (`awk` over the MANIFEST's role column): 4 current, 59 historical, 1 historical-classic -
matches the memo exactly. By family: equivalence 33 (1 historical-classic, 31 historical, 1 current),
bcf-equivalence 8 (7 historical, 1 current), multinomial-equivalence 12 (11 historical, 1 current),
bench-sampler 11 (10 historical, 1 current). Host/arch column: 64/64 rows say "arm64 (macOS)" - single
host and arch for every recorded baseline ever made; no x86 or Linux row exists.

historical-classic: equivalence-5430fdb.rds, 9 scenarios, recorded 168665a8 (07-02) under the classic
(BayesTree-descended) engine, deleted at b354f3ab (07-03). This is the ONLY row in all 64 anchored to
anything other than bartcore's own prior output; its provenance is PERMANENTLY IRREPRODUCIBLE - the
source that produced it no longer exists at any reachable commit. Quantified cross-engine gap vs its
bartcore-side pair (equivalence-b354f3a.rds): 9 scenarios, max |z| 3.83 over 329 summaries. Two rows
(equivalence-19c499a.rds, equivalence-fbd2168.rds) were mislabeled historical-classic for weeks,
relabeled historical 2026-08-17 after an audit found they were recorded 07-05, 2 days after deletion.

bcf-equivalence and multinomial-equivalence have ZERO historical-classic rows at any point in their
history - both features were built entirely under bartcore, so every row in both families, from their
first recording onward, is bartcore output re-recorded against bartcore output; "the compare's own
statistical check IS the oracle" is the MANIFEST's own phrase for the current equivalence-5a3bc276.rds
row (43 scenarios, max |z| 2.06 across 87 summaries on the one moved scenario).

Current rows: equivalence-5a3bc276.rds (43 scenarios), bcf-equivalence-6e3b9fb8.rds (12), multinomial-
equivalence-4d9a3337.rds (11), bench-sampler-ab1dc52.csv (speed) - the 3 .rds files equivalence.yaml
compares against (never run in CI, Sec.1).

## 4. Exact-posterior gates (exact-gates.yaml, 20 scripts, quick-mode only in CI)

Tolerance is computed in-script per gate (z-score bound or absolute posterior-moment gap), not a shared
constant (e.g. bcf-exact.R: 0.05 quick / 0.015 full, absolute). Scale, by grep of each script's
n.trees/numTrees literal:

- Single-tree (n.trees=1L or numTrees=1L hardcoded, in BOTH quick and full mode - "full" only
  lengthens draws/seeds, never ensemble size): 15 of 20 - bd-balance, change-balance, swap-balance,
  aft-exact, bcf-exact, bcf-exact-weak, bcf-exact-restricted, categorical-exact, linear-exact,
  hazard-exact, hurdle-exact, t-exact, monotone-reference, negbin-exact, ordinal-exact.
- Multi-tree: heteroscedastic-exact.R (mean forest 20L, variance forest 1-2L), multinomial-exact.R
  (mixes 1L and 50L across scenarios), hazard-reduction.R (40L fixed), hurdle-reduction.R (30L fixed),
  logistic-reference.R (ntree = 50L quick / 200L full - the only gate reaching bart()'s classic default
  of 200, but full mode has never run in CI, so no CI-verified run of ANY gate has reached 200 trees,
  and none reach dbarts()/bart2()'s own default of 75 either).

Last CI run: 2026-08-24 (tip push), success, quick mode - as for all 200 runs on record.

## 5. SBC (sbc.yaml + benchmarks/R/sbc.R)

6 arms in the CI matrix: discrete-selfcheck, gaussian, ordinal, nbinom, t, multinom - matches the
memo's "six-arm". Default ensemble size nTrees=50L (one config override at 25L); "ensemble" only
relative to the single-tree exact gates above - neither shipped default (75 or 200). Excluded from the
matrix (workflow's own header): BCF (slow-mixing glue-sigma channel is a recorded finding, not gated)
and aft/hazard/hurdle/heteroscedastic/monotone (sbc-family-tiers.md: ill-posed or not liftable yet).

Never run in CI (Sec.1). One manual local run exists: gaussian arm only, `Rscript benchmarks/R/sbc.R
gaussian 200 200 30`, at a39da5d9 (2026-08-17, docs/plans/release-candidate-review.md "Wave 0c"),
SBC_FAIL_ON_FLAG=1, ALL SEVEN FUNCTIONALS PASS (avg.f chisq p 0.294, f.star1-5 chisq p 0.255-0.849,
sigma chisq 0.962). The other 5 arms (ordinal, nbinom, t, multinom, discrete-selfcheck) have no
recorded run at all, manual or CI, anywhere in docs/plans.

## 6. rchk / valgrind / TSAN / gctorture

rchk.yaml and valgrind.yaml exist (Sec.1) but have 0 CI runs ever. Both were run manually/locally once,
2026-08-19, docs/plans/release-candidate-review.md:
- rchk: initial FAIL (11 [UP], triaged as all false positives - rchk does not model rooting via
  SET_VECTOR_ELT or via attribute access on a rooted argument), fixed with a 10-line defensive patch
  (69de27ac), then a merged-tip run at 8dbc0ce9 CLEAN: 0 [UP], 0 [PB], only the 5 expected
  analyzer-capacity bailouts.
- valgrind: found and fixed 2 real defects - a definite leak (2464 bytes, applyBCFSpec/unwindProtect
  closure body, fixed 7e42dc93) and a pre-bartcore (since 2015, R's own makeModelMatrixFromDataFrame.c)
  OOB heap read reachable via bart()/bart2() extra-factor-level inputs, fixed 7be7a126. Final re-run at
  7be7a126: "All ok, 6419 results", 0 errors, 0 lost bytes.

TSAN: no workflow, script, or reference anywhere in .github/workflows/ or benchmarks/ - does not exist
in this repo in any form. gctorture: no workflow; used exactly once, manually, in the same 2026-08-19
rchk triage, to confirm 4 PROTECT sites empirically - a one-off aid, never a repeatable gate/script.

## 7. revdep / consumer bundle

revdep-smoke.yaml (never run in CI, Sec.1) checks 3 packages against their COMPAT BRANCHES, not CRAN
releases (stan4bart@bartcore, bartCause@dbarts-1.0, treatSens@dbarts-1.0) - CRAN releases predate the
bartcore rework and cannot pair with dev dbarts. One manual run exists, covering 4 packages (P12,
142a50b6, 2026-08-18, tip 8e1e674c, each in a private lib from a git-archive stage): stan4bart clean
(tinytest 531/531, R CMD check 0E/0W/0N), treatSens clean (testthat 186/186, matches --as-cran
baseline), bairrtt clean (206/206, NOT in the CI matrix), bartCause found 2 issues, neither a dbarts
regression (stale RNG-shift snapshots, regenerated; a stan4bart-version mismatch in the tester's lib).

## 8. tinytest hardcoded-value / snapshot pins

6 files carry RNG-locked literals, 34 expect_ assertions total (7+3+14+1+3 in the 5 files below, +6 in
the 6th) - matches the memo. 5 are covered by tools/regenerate-snapshots.R (its own file list, lines
26-30) and share a header declaring themselves a "Seeded-drift tripwire, not a correctness test":
test-reproducibility-{binaryResponse,continuousResponse-multithreaded,continuousResponse-
singleThreaded,rbart,xbart}.R. The 6th, test-rbart-loop-callback.R (15 literals, 6 expect_ sites,
tolerance 1e-8 per its own MANIFEST pin-site note), is explicitly NOT covered by regenerate-snapshots.R
and is hand-regenerated on its own schedule. The tool itself is invoked by no CI workflow (Sec.2) - by
design: a writer, run only after a human confirms a draw-shift is intentional; a CI-invoked regenerator
would silently launder a real regression into a new "expected" value.

## Summary table (gate - proves - anchor - last run - blind to)

- check-standard/cpp-tests.yaml: build + inst/tinytest/tests/cpp pass - NONE (regression only) -
  2026-08-24 push, success - wrong-but-stable computation; PR-only paths (never exercised).
- exact-gates.yaml (quick): draws match an in-script closed-form target at m=1 for 15/20 gates -
  closed-form/analytic - 2026-08-24 push, success (always quick) - defects only visible at ensemble
  scale (75/200 trees; full mode never run).
- sanitizers.yaml: no ASAN/UBSAN finding under tinytest - dynamic instrumentation - 2026-08-24 push,
  success - leaks (disabled), data races (no TSAN).
- lint.yaml/pkgdown.yaml: style/build hygiene, not correctness - n/a - 2026-08-24 push, success -
  everything semantic.
- equivalence.yaml: draws bitwise-identical or posterior unmoved (|z|<4) vs a recording made BY THIS
  ENGINE at a prior sha, one arm64 mac - SELF for 63/64 rows, classic-engine for 1 - NEVER run in CI -
  classic-engine-shared defects, non-macOS divergence, anything since bcf/multinomial inception.
- sbc.yaml: recovered-parameter ranks are uniform, 6 arms at nTrees=50 - SBC self-consistency - NEVER
  run in CI; gaussian arm manual pass, a39da5d9, 2026-08-17 - shipped ensemble scale (75/200); other 5
  arms never actually run at all.
- backfit-exact.R/geweke-mc.R (orphan): a sampler step/whole joint dist. matches a hand-derived
  closed form - closed-form/MC test - each run once manually, 2026-08-17/18, never CI - anything
  outside their one fixed scenario.
- composition-matrix.R (orphan): behavior matches feature-matrix.md's own claim - internal (code vs
  its own docs) - run once, 2026-08-18, never CI - a claim that is itself wrong.
- mutation-battery.R (orphan): suite catches a planted-defect catalog - mutation testing - undated
  last run, never CI - mutations outside its catalog.
- rchk.yaml: no PROTECT-balance defect in the bridge - static analysis (CRAN's tool) - NEVER run in
  CI; manual clean, 8dbc0ce9, 2026-08-19 - logic bugs, leaks, races.
- valgrind.yaml: no invalid access/leak under the full suite - dynamic instrumentation - NEVER run in
  CI; manual full-suite clean, 7be7a126, 2026-08-19 - logic bugs, races, non-arm64/ubuntu platforms.
- TSAN: does not exist - n/a - never - all data races.
- gctorture: no workflow; one manual investigation, 2026-08-19 - dynamic (forced GC) - anything not a
  missing PROTECT.
- revdep-smoke.yaml: 3 consumers build+check against dev dbarts from COMPAT BRANCHES - real code, but
  a pre-release snapshot - NEVER run in CI; 4-package manual run clean, 142a50b6, 2026-08-18 - CRAN
  consumer versions; bairrtt (manual only, not the CI matrix).
- tinytest snapshot pins (6 files, 34 assertions): draws match this code's own prior output - SELF,
  disclaimed as "not a correctness test" - every push - a wrong-but-stable answer by construction.

## Findings

1. 5 of 11 workflows (equivalence, sbc, rchk, valgrind, revdep-smoke) have ZERO runs ever - schedule
   is bound to main (lacks the files), nobody has ever dispatched them manually either. RECOMMEND:
   land these 5 on main now (no compat constraint blocks it), or dispatch each by hand on a cadence.
2. exact-gates.yaml has 200 CI runs, all push-triggered quick mode; full mode has never run once.
   RECOMMEND: dispatch full mode at least once before 1.0 and record the result.
3. 63 of 64 baseline rows, and both the bcf-equivalence and multinomial-equivalence families in their
   entirety, are bartcore recording itself; the sole external anchor (equivalence-5430fdb.rds, classic
   engine) is permanently irreproducible since its source is deleted. RECOMMEND: add one independently-
   implemented reference check for a core gaussian scenario before 1.0, or document why none is planned.
4. Two MANIFEST rows were mislabeled historical-classic for weeks (caught 2026-08-17, both actually
   recorded after the classic engine's deletion). RECOMMEND: add a MANIFEST lint checking
   recording-commit ordering against the classic-engine-deletion commit.
5. 15 of exact-gates.yaml's 20 gates run at n.trees=1 unconditionally; no CI-verified exact-posterior
   run has ever reached either shipped ensemble default (75 or 200 trees). RECOMMEND: add an
   ensemble-scale arm to the gaussian/probit/multinomial gates.
6. 14 of 38 benchmarks/R harnesses and all 4 tools/ scripts are invoked by no CI workflow. Five
   (backfit-exact.R, geweke-mc.R, composition-matrix.R, mutation-battery.R, grouped-mixing.R) call
   themselves gates/oracles/permanent harnesses in their own headers, yet have exactly one recorded run
   each, all manual, 2026-08-16 through 08-18. RECOMMEND: wire these five into a weekly CI workflow.
7. change-fix-instrumentation.R, change-fix-stage2.R, and parallel-falsifiers.R cannot currently
   produce output at all - each depends on src/-side env-var instrumentation confirmed absent from
   src/ today. RECOMMEND: delete; each already served its one-time design decision.
8. linear-leaf-xcheck.R is referenced nowhere in the tree outside itself. RECOMMEND: delete, or add
   one line to benchmarks/README.md.
9. rchk and valgrind both passed clean under manual local runs dated 2026-08-19 (8dbc0ce9, 7be7a126)
   but their CI legs have never executed and cannot until main carries the files. RECOMMEND: see #1.
10. No TSAN leg exists anywhere in the repo. RECOMMEND: no immediate action (ASAN+UBSAN cover per-push
    memory safety); revisit if the engine grows real threading beyond n.threads dispatch.
11. revdep-smoke.yaml's CI matrix (3 packages) is narrower than the one manual run that has ever
    happened (4, adds bairrtt). RECOMMEND: add bairrtt to the CI matrix or document the exclusion.
12. Zero pull_request-triggered runs exist in this repo's history; 6 workflows declare a pull_request
    trigger that has never fired, since all work lands via direct push to bartcore. RECOMMEND: no code
    action; confirm with VD whether PRs are ever planned pre-merge-to-main.
13. All 64 recorded baselines are single-host (arm64 macOS) - none on x86 or Linux (dbarts-x86-bench-
    box note). RECOMMEND: no action while the x86 box is unavailable; CI's ubuntu/windows legs are
    host-portable by design but carry no bitwise baseline comparison.
