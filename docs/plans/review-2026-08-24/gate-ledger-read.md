# Independent read of gate-ledger.md - bartcore @ b102e17c, 2026-08-24

Refute-then-extend. Verdicts are from my own commands (full 3091-run `gh api .../actions/runs --paginate` dump; all 11
trigger blocks read whole; MANIFEST/harness/header greps). Cross-implementation anchoring: see anchor-main.md.

## 1. Verdicts

CONFIRMED:
- 11 workflows on bartcore, 1 on main (`git ls-tree main` = check-standard.yaml; default=main).
- 3091 runs; ZERO pull_request/schedule/workflow_dispatch (`cut -f4|uniq -c` = 2627 push, 464 dynamic). Branches 2621
  bartcore / 464 gh-pages / 6 main.
- equivalence, sbc, rchk, valgrind, revdep-smoke: zero runs. Path x event has only check-standard 519, cpp-tests 496,
  lint 495, pkgdown 491, sanitizers 425, exact-gates 200, abi-contract 1.
- exact-gates 200 runs all event=push, so mode=full (a dispatch input) never ran. check-standard 488/25/6 = 519.
  abi-contract footnote exact. revdep matrix = the 3 compat branches, bairrtt absent. bench-sampler-ab1dc52.csv dated
  2026-08-08.
- benchmarks/R = 38, 14 orphans, membership identical; tools/ = 4 .R, none wired. Broken 3 confirmed: `grep -rl src/`
  = 0 for DBARTS_CHANGE_LOG/_STATS, BC_FALSIFIER, DBARTS_MOVE_SUMMARY, _CHANGE_KERNEL.
- MANIFEST 64 rows / 64 files, exact 1:1 (`comm -3` empty); 4 current / 59 historical / 1 historical-classic; family x
  role exact; 64/64 "arm64 (macOS)". equivalence-5430fdb.rds is the only historical-classic row (9 scenarios, max |z|
  3.83, 329 summaries, the 19c499a/fbd2168 relabel).
- No external reference implementation: BayesTree 0 hits in benchmarks/R, absent from Suggests. NUANCE the ledger
  omits: equivalence-4a42620a.rds "reproduces the implementer's independent own-build recording 42/42" - an
  independent BUILD of the same source, not an implementation, so finding #3 stands. 15 of 20 exact gates at
  n.trees/numTrees=1L unconditionally (re-derived by grep). SBC 6 arms as listed, SBC_FAIL_ON_FLAG set, nTrees=50L,
  one 25L override (sbc.R:2175).
- tinytest pins: 6 files, 7+3+14+1+3=28 plus 6 = 34 expect_ exactly; regenerate-snapshots.R lists 5. I swept all 167
  tinytest files for >=6-decimal literals; the 4 extra hits are re-derived analytic constants, comments, or a string
  check - not RNG pins. "6 files" is right.

WRONG:
- "PR trigger targets main/master only". lint.yaml and pkgdown.yaml declare a BARE `pull_request:`, no branch filter,
  so both accept a PR into any base including bartcore; sanitizers is [bartcore, main]. A PR into bartcore fires 3 of
  the 6 PR-declaring workflows but NOT R-CMD-check, cpp-tests, or exact-gates. Finding #12 is a real asymmetry, not
  "no action".
- "none reach dbarts()/bart2()'s own default of 75". multinomial-exact.R arm 7, line 625, runs `n.trees = 75L`
  UNCONDITIONALLY, quick mode, every push. SUBSTANCE SURVIVES: x is constant there, so no column has a valid cut point
  and every tree is a root - 75 stumps, no ensemble structure. logistic-reference.R is also BOTH classes (an
  n.trees=1L arm beside ntree=50/200). The 200-tree half is CONFIRMED.

OVERSTATED:
- "5-way matrix plus a standalone re-run of test-simd.R". check-standard.yaml has TWO jobs: the 5-config matrix AND
  `windows-arm64-neon` on windows-11-arm - a SIXTH platform class (native aarch64) that asserts NEON then pins
  test-simd.R at an exact count. A platform leg, not a re-run.
- "rchk run once". Two campaigns on 08-19: on MAIN for the CRAN patch (48 [UP] -> zero real findings at cb290550,
  rc-review:731-745) and on bartcore (11 [UP] -> 69de27ac -> 8dbc0ce9 clean, :1138-1172).
- valgrind. The manual run was NOT the workflow's configuration: native arm64 ubuntu:24.04 + stock R 4.3.3, because
  "the r-hub amd64 image crashes under Rosetta emulation on this host" (rc-review:1072), while valgrind.yaml runs the
  r-hub amd64 container. So valgrind has never run in its declared configuration anywhere. "All ok, 6419 results" vs
  6438 on the dev host = 19 skips.
- exact-gates is 193 success / 2 failure / 5 cancelled. linear-leaf-xcheck.R: `grep -rl` returns NOTHING - absent even
  from its own file. bench-sampler is documented in 50 docs/plans files, not 9 (understated).

## 2. Missed

- A. THE FLAT C API HAS THE BRANCH'S STRONGEST GATE AND THE LEDGER HAS NO ROW FOR IT. src/C_interface.cpp:461,:465 are
  two constexpr FNV-1a static_asserts - a signature token, and a full fold over signatures + both ABI enums + the
  callback parameters + the compiler-reported LAYOUT (size, each field name paired with its offset) of the three
  crossing structs - checked against DBARTS_C_API_HASH (dbarts.h:142), beside offsetof asserts at :331. Any covered
  ABI change FAILS THE BUILD until the literal is re-baked, on every push, in every workflow that compiles
  (check-standard's 6 platforms, cpp-tests, exact-gates, both sanitizer legs). Beside it test-capi.R (233 assertions)
  compiles a real consumer.c against the INSTALLED headers, resolves every entry via R_GetCCallable as a LinkingTo
  consumer would, and force-compiles a stale token to require the refusal; its exit_file() paths all hard-`stop()`
  under CI, so it cannot silently skip. Residue: a same-width in-place type swap moves no offset and no name.
- B. RCHK.YAML IS FALSE-PASSING BY CONSTRUCTION AND THE CURE IS ALREADY WRITTEN DOWN. rc-review:1170 records a durable
  GATE LESSON - an empty [UP] grep PASSES on an OOM-killed analyzer (observed: "Killed" left an empty .bcheck), so the
  gate must also require the "Analyzed N functions" line and the absence of "Killed". rchk.yaml does neither: `find
  ... 2>/dev/null | grep -E '...' || true`, failing only on a non-empty match. No .bcheck at all is a green job.
- C. SANITIZERS.YAML CARRIES THE ONLY ANTI-HOLLOWING GATE IN CI, UNCOUNTED: a total tinytest floor (>=5200, baseline
  5799 at f6c8979d), six per-file floors (test-capi.R 150, -bartcore 160, -multinomial-surface 80, three more), and a
  log grep for "runtime error:|AddressSanitizer|SUMMARY: .*Sanitizer", so a non-assertion UBSAN finding still fails.
  It is alone because R CMD check runs tests opaquely.
- D. SANITIZERS HAS NO MAIN BINDING: push branches [bartcore], while every other push workflow lists [main, master,
  bartcore]. On merge to main the per-push memory-safety gate AND floor C stop.
- E. THE EVIDENCE SUBSTRATE IS UNTRACKED. `.claude/` is gitignored (.gitignore:27) and `git ls-files .claude/` = 0,
  yet tracked files cite into it 120 times - including MANIFEST, which names
  .claude/rc-review-p10-baseline-provenance-2026-08-17.md as its authority for the historical-classic relabel. The
  rchk logs, the valgrind final log and the orphan-leg audit all live there; on a fresh clone none of it exists.
- F. PROSE-ONLY GATES (a one-time agent run is not a gate), instrumentation confirmed absent from src/: (1) the
  "instrumented counter build" measuring equal-rank realizations (2 maskprobit, 1 maskordinal, 0 elsewhere) - the
  MANIFEST's CITED ORACLE for the equivalence-d15a2bfb.rds re-record, unreproducible by anything in the tree (the
  qualitative law IS gated, by tests/cpp testEqualRankOneComparison and bd-balance.R's veto arm); (2)
  monotone-prior-draw.md:73-98, instrumentation "removed before landing" giving the P(cap exhausted)~2e-11 that
  justifies the 1e6 retry cap; (3) setpredictor-leafof-rebuild.md:74-98, an rdtsc + md5 probe absent from the tree;
  (4) bcf-bartcause-relocation.md, ~50 MEASURED tags from a .claude/ probe.
- G. THE SPEED GATE IS 1225 COMMITS STALE, 327 IN src/ (`git rev-list --count ab1dc52..HEAD`). "16 days" understates
  it. The equivalence baseline 5a3bc276 is 24 commits / 7 src/ behind - healthy.
- H. MEASURED HIT RATE, NEVER REPORTED. Failures over full history: lint 62, pkgdown 26, check-standard 25, sanitizers
  19, exact-gates 2 (08-06, 08-19), cpp-tests 2. The per-push legs turn red ~4.5-4.8% of runs - proven live; the five
  never-run workflows have had zero opportunities.
- I. THE RECORDED BLAS-SENSITIVITY CLASS IS COVERED AND THE MATRIX CAUGHT IT. cf4a290c went red on all three ubuntu
  R-CMD-check legs and both sanitizer legs while macOS and both Windows legs passed; container-reproduced on amd64 R
  4.6.1 + OpenBLAS, NOT on native arm64 (rc-review:986-1000). Three BLAS classes are covered and the one that bit is
  the one that caught it. The ledger scores it "anchor NONE".
- J. TESTS/CPP IS NOT "ANCHOR NONE" EITHER: main.cpp:1 declares testing "against independently coded reference math",
  and it delivers - test_scan.cpp checks log-likelihoods against an independent brute-force per-cut recompute,
  test_interaction.cpp brute-forces variableAvailable off a hand-built p-bitset. A genuine independent oracle, on
  every push.
- K. THE SOLE EXTERNAL ANCHOR IS A RECORDING OF A KNOWN-DEFECTIVE ENGINE. correctness-audit.md found a MAJOR
  change-move proposal-density omission of CGM lineage, INHERITED FROM THE DELETED CLASSIC ENGINE, biasing tree
  structure toward low-cardinality variables (confirmed by change-balance.R). equivalence-5430fdb.rds records that
  engine: cutover evidence, not a correctness reference.
- L. (fixed-at 74e2e050: both scripts named below are now wired - check-win-drift.R and check-rc-codoc.R from
  lint.yaml's doc-freshness job, check-doc-freshness.R from its own no-paths-ignore doc-freshness.yaml.)
  TWO tools/ SCRIPTS COVERED FAILURE MODES NOTHING ELSE SAW. check-win-drift.R is the only check of the
  hand-maintained src/*.win config headers against their *.in templates and DESCRIPTION (Windows has no configure,
  nothing regenerates them, drift is silent). check-doc-freshness.R resolves every file:line anchor in
  feature-matrix.md against the live tree. Both are fast, deterministic, host-portable, and wired to nothing.

## 3. Ranked remediation

1. agent-fix, Sonnet, ~10 lines. Harden rchk.yaml per its own recorded lesson (require "Analyzed N functions"; fail on
   "Killed" or a missing .bcheck), so a gate that today passes on an OOM-killed analyzer cannot. Miss B; must precede
   trusting any rchk result.
2. agent-fix, Sonnet, 1 line. Add [main, master] to sanitizers.yaml's push branches, so the per-push memory-safety
   gate and the ONLY suite-count floor survive the merge to main. Miss C/D.
3. VD-judgement then agent-fix, Sonnet, ~5 lines + 5 dispatches. Land the five schedule-bound workflows on main
   (main-side action is VD's alone) or dispatch each once now and record it: five declared gates become evidence. Item
   1 first. Findings #1/#9; equivalence-ci.md's own landing note already defers "the workflow's first live run", so
   this is a recorded open item, not an oversight.
4. agent-fix, Opus, one dispatch + a recorded note. Run exact-gates mode=full once. Buys: the only 200-tree
   exact-posterior evidence that will exist, plus the tighter full-mode tolerances. Finding #2.
5. DONE (fixed-at 74e2e050): check-win-drift.R runs from lint.yaml, check-doc-freshness.R from its own
   doc-freshness.yaml. Bought: coverage of silent Windows config drift and of dead doc anchors. Miss L closed.
6. agent-fix, Sonnet, small. Add bartcore to the PR filters of check-standard/cpp-tests/exact-gates, or drop the bare
   `pull_request:` from lint/pkgdown, so a bartcore PR stops getting style checks without correctness checks. Confirm
   PR intent with VD first.
7. agent-fix, Sonnet, 3-4 deletions. Delete change-fix-instrumentation.R, change-fix-stage2.R, parallel-falsifiers.R
   (cannot execute); delete or document linear-leaf-xcheck.R.
8. VD-judgement, sizeable. Move the .claude/ artifacts that TRACKED files cite - the p10 baseline-provenance audit
   MANIFEST names, the rchk .bcheck logs, the valgrind final log - into a tracked tree, or rewrite those citations to
   stand alone, so provenance survives a clone. Miss E.
9. agent-fix, Opus, ~1 day. Add an ensemble-scale arm (75 trees, real cut points) to the gaussian and probit gates:
   the first exact-posterior evidence at shipped scale, since miss A shows the existing 75-tree arm is 75 stumps.
   Finding #5, re-scoped.
10. agent-fix, Sonnet, quiet machine. Re-record bench-sampler at tip; closes a 327-src-commit gap. Miss G. Needs the
    quiet-machine grant, so it schedules around VD.
11. VD-judgement. Whether to keep the 59 historical rows and whether to re-record against an external anchor - defer
    to anchor-main.md, weighing miss K: the classic engine carried a confirmed change-move defect, so it is cutover
    evidence only.
12. agent-fix, Sonnet, ~40 lines. Wire backfit-exact.R, geweke-mc.R, composition-matrix.R, mutation-battery.R and
    grouped-mixing.R into one weekly workflow. Finding #6; ranked here because it needs item 3 to fire at all.
13. defer: TSAN (#10), bairrtt in the revdep matrix (#11), an x86 baseline (#13), a MANIFEST recording-order lint
    (#4). Each is real; none blocks 1.0-0.

## 4. For the maintainer

What the gates prove today: the branch compiles and passes its full 167-file tinytest suite on six platform classes
spanning three BLAS implementations on every push, and it has demonstrably caught real defects (25 red R-CMD-check
runs, 19 red sanitizer runs, one a genuine OpenBLAS-only numeric-sensitivity bug); the engine's component math agrees
with independently coded brute-force oracles in tests/cpp on every push; the shipped flat C API cannot change a
signature, an ABI enum, its callback, or a crossing struct's layout without failing the build, and a real consumer is
compiled and driven against the installed headers on every push; and twenty deterministic exact-posterior gates match
in-script analytic targets. What they do NOT prove: nothing has ever exercised the statistical gates in CI
- equivalence, SBC, rchk, valgrind and revdep-smoke have zero runs between them, so calibration evidence is one manual
  gaussian SBC arm with five never run at all, and rchk/valgrind rest on local runs logged to an untracked directory,
  valgrind's never in its declared configuration; the exact gates are almost entirely single-tree, so the ensemble
  interaction that is the point of BART has no closed-form gate at shipped scale; all 64 baselines are one arm64 mac
  recording itself, the sole cross-engine anchor being an engine later proven to carry a real change-move defect; and
  the speed baseline is 327 engine commits stale. Read the green checks as strong evidence of memory safety,
  portability, API stability and per-step engine math, and as none about calibration at ensemble scale or behavior on
  any machine but this one.
