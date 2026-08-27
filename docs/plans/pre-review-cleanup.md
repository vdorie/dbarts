# Pre-review cleanup: four adversarial reviews and their rulings

Status: LANDED 2026-08-26 at 7cd71f2d (base 3080a9c5; a final anchor pass
follows). Ten commits, 2518a0b5..7cd71f2d.

## 1. Charter

VD, 2026-08-26: before VD's own manual read of the branch
(docs/plans/bartcore-review-tour.md), run four adversarial reviews the tour
itself cannot substitute for - staleness (is the tour and feature-matrix.md
still true of the tree), completeness (defects a pre-review pass should catch
so VD does not have to), YAGNI (speculative generality nobody named a use for),
and agent accumulation (comment/records/test sediment left by many independent
sessions writing to the same files). VD's bar for the tour specifically, as
quoted in the staleness audit: "I'm the consumer, and it needs to be worth my
time to read before we merge bartcore into main" - a merge document, not an
audit syllabus.

Four reviews were commissioned and are tracked verbatim (findings unedited;
paths scrubbed of session-local locations, ASCII, a one-line header naming
this doc as the record of VD's rulings on them):

| review | file | what it covers |
|---|---|---|
| staleness | docs/plans/review-2026-08-24/memos/pre-review-audit-review-chain.md | Every citation and count in bartcore-review-tour.md and feature-matrix.md checked against the live tree at 3080a9c5; a proposed 7-section restructuring of the tour as a merge document. |
| completeness | docs/plans/review-2026-08-24/memos/pre-review-critic.md | Twelve PRE-REVIEW defects, verified live, the first five user-visible; a RESIDUE-OK list of correctly-deferred items. |
| YAGNI | docs/plans/review-2026-08-24/memos/pre-review-yagni.md | 24 numbered candidates for speculative generality, each stated with its strongest value claim before a verdict. |
| accumulation/cruft | docs/plans/review-2026-08-24/memos/pre-review-cruft.md | Sections A-J: R/C++ comment load, helper factoring, multi-agent blame layering, records/CI residue, engine dead-code and duplication census, test-suite redundancy, the "duplicate by design" ruling on the equivalence harnesses. |

## 2. Rulings

Every numbered fork the four reviews put to VD, and the ruling. "Landed at"
is the commit in section 3; a ruling with no commit is recorded here but not
executed by this stack (an open door, section 5, or a decision that changes
no file).

| # | fork | source | ruling | landed at |
|---|---|---|---|---|
| 1 | Archive the 118 landed plan docs with no present-facing reader | cruft H1, H3, H4 | YES | be931ec7 |
| 2 | `src/misc`'s two thread managers + threaded moment wrappers (~2470 lines, zero callers) | yagni #1 | ARCHIVE (branch `archive/misc-thread-managers` at 3080a9c5), then CUT from bartcore | c0ddd93e |
| 3 | Dead harness files (#2), the `engine=new` shim (#3), the `AmplitudeSpec` test fixture (#9), `sliceSample`'s `log = FALSE` arm (#10), `unpack`/named `massign` (#11), `SamplerOptions::forestColumns` (#15), `PredictorUpdateResult::unsupportedSource` (#18), `SamplerShape::varianceLeafPrior` (#19) | yagni | YES, cut | c0ddd93e, 89577d7d |
| 4 | `makeind(x, all = )` - a formal that does nothing, worst of FINISH/CUT/as-is | yagni #12 | KEEP as documented no-op; neither implemented nor dropped this round | (none - residue, sec 6) |
| 5 | Six post-hoc `if (!engine->set*(...)) Rf_error(...)` guards the predicates already make unreachable | cruft E9 | DELETE (the "no unreachable branches" house call; the one unreachable fractional-weight message fixed either way) | 89577d7d |
| 6 | `cutPoints` (split-candidate grid) vs `cutpoints` (ordinal latent thresholds) - one concept, two names one capital letter apart, the second shipped in `dbarts.h` | cruft E7 | YES, rename the internal/C-API one to `ordinalThresholds`; the split grid's `cutPoints` untouched | 89577d7d |
| 7 | `bcf-equivalence.R`/`multinomial-equivalence.R`'s ~250-400 shared lines, previously ruled duplicate-by-design | cruft G1-G3 | OVERTURNED - share via `equivalence-common.R`; `compareCrossHost` differed between the two files by six `sprintf` column-width characters | 754309f1 |
| 8 | Three CI/gate defects: doc-freshness's own `paths-ignore` excludes its declared inputs (H9); the BCF/multinomial cross-host jobs in `equivalence.yaml` duplicate `exact-gates.yaml`'s per-push compares (H10); `valgrind.yaml` installs neither `survival` nor `posterior`, so 37+17 tinytest files silently skip under it and still report green (H10) | cruft H9, H10 | YES, all three, plus the adjacent packaging items (`tools/*.R` shipped in the tarball, `master` referenced where no such branch exists, the kernels Makefile's missing header-dependency tracking) | 2518a0b5 |
| 9 | The R comment corpus's provenance/derivation/defense prose (93 of 177 blocks over six lines) | cruft A1-A2 | YES, cut the offending sentences, keep the constraints | 6fb008ba, 0f5b7b0c |
| 10 | GP leaves (773 engine lines, one `FunctionLeafModel` instantiation, an O(n_leaf^3) Cholesky, zero consumers on any branch) | yagni #4 | KEEP - the value claim holds (a published BART variant, R-reachable and tested); flagged for a deliberate "is this in 1.0" conversation, not cut | (none) |
| 11 | The flat C API's forest/amplitude block, weakest two members `dbarts_sampler_getForestCalibration` and `dbarts_sampler_setForestPriorScale` | yagni #5 | KEEP the whole block, including both weak members - `dbarts-h-reshape.md`'s named enabled models are the gate, and consumer absence is not the gating fact | (none) |
| 12 | `DESCRIPTION`'s `Authors@R` credits code no longer in the tree (the classic engine's radix tree; the `m4` macros `ax_cxx_namespace_std.m4`, `ax_func_posix_memalign.m4`; file lists on live contributors naming absent files) | cruft H12, yagni #23 | REMOVE the credits for code that is gone; trim the file lists that name absent files; keep every credit for code still present | c0ddd93e |
| 13 | `bartcore-review-tour.md`'s reading order and length (975 lines / ~40-50 min, orientation front-loaded, evidence for the merge decision back-loaded) | audit-review-chain sec 5 | The audit's proposed 7-section restructuring (the merge in one page; what a user's script sees; what a linked consumer sees; gates and what they prove; risks a merge carries; doors needing VD's ruling; the code walk) is the going-forward structure. Not executed this round (sec 5, TODO `second-review-followups`). | (none) |
| 14 | Formal heredity (ensemble lattice prior) as an interaction-constraints door | TODO `interaction-constraints`, restated by yagni #? (settled call, not relitigated) | The first post-1.0 arc | (none - already TODO'd; unchanged by this stack) |

Also ruled, as one block rather than one row each: the twelve PRE-REVIEW
completeness defects (critic-pre-review.md, including the gaussian
case-weight bug cruft E10 found independently) - FIX ALL NOW, before the
manual review, not deferred as residue. Landed 74e2e050 (section 3).

## 3. Per-commit landing

| sha | what | gate evidence |
|---|---|---|
| 2518a0b5 | doc-freshness.yaml split out of lint.yaml with no paths-ignore; tools/\*.R .Rbuildignore'd; equivalence.yaml's duplicate BCF/multinomial cross-host jobs dropped (exact-gates.yaml already runs them every push); master -> [main, bartcore] in six workflow files; valgrind.yaml installs survival+posterior; benchmarks/kernels/Makefile gains -MMD -MP dependency tracking; MANIFEST's 7 absent-file baselines annotated `(recording not retained)` | YAML/Makefile/text only, no R or C++ touched; no sampling code, no gate to run beyond the doc-freshness workflow itself firing on the next docs-only push |
| 6fb008ba | Cut provenance/derivation prose from 8 R files (R/bartcore.R, data.R, dbarts.R, mixedMatrix.R, model.R, sliceSample.R, spec.R, xbart.R), -83 comment lines net across this and 0f5b7b0c, keeping every stated constraint | Comment-only; no code line changed; verified by inspection against cruft A2's five worst-block excerpts, all addressed |
| 17b971e2 | docs/design/feature-matrix.md's anchors re-derived by content against 3080a9c5 (36 cells) | docs/design is out of this task's scope; landed by the design-side agent inside this stack |
| 74e2e050 | The twelve critic-pre-review.md defects plus the independently-found gaussian-weight bug (cruft E10/Item 0): n.thin driving printEvery to 0 floored at 1 and dbartsControl's validity refusing printEvery = 0 directly; the extract/plotTree catch-alls fixed to also refuse positional extras; predict's slot-4 offset formal added to bartOrdinal/bartHurdle; fitted's ci.level/sample message now names the moved argument; dbarts_sampler_setVerbose refuses printEvery = 0 instead of dividing by it; enforceBinaryWeightPolicy validates gaussian weights (finite, non-negative) and is now the sole check at every entrance, including the flat dbarts_sampler_setWeights, which previously installed them unchecked; the two cross-host equivalence scripts' skip arm now sets anyFailure; test-plot-generics.R's vacuous plot.bartMultinomial pin replaced with one that fails under the old argument order; six error-style slips fixed; summary's empty-table message interpolates the requested vars; the dead GPGaussianLeaf::maxLeafSize() accessor and the stale sparseFactor.R header comment fixed; two stale doc anchors (docs/design/multinomial-mutation-arc.md, docs/plans/bcf-cross-host.md) corrected; two internal-linkage forwarders marked static | tinytest 7475/0 and tests/cpp all 265 ok (both re-run at the stack tip, section 4); the gaussian-weight fix is exercised by the new refusal, not an equivalence scenario (no scenario supplies a negative gaussian weight) |
| c0ddd93e | -4821 net: src/misc's two thread managers (hierarchicalThreadManager.c, blockingThreadManager.c) and ~1850 lines of threaded moment wrappers in moments.c cut (preserved on branch archive/misc-thread-managers at 3080a9c5); the engine=new shim (bartcore-shim.R, rshim.cpp, the equivalence.R/bench-sampler.R flag); benchmarks/R/linear-leaf-xcheck.R + benchmarks/kernels/linear_leaf.cpp; ordinal-cutpoint-mixing.R, negbin-r-update-mixing.R; R/multipleAssignment.R's unpack; sliceSample.R's log = FALSE arm cut from R/sliceSample.R and pruned from test-slice-sample.R; SamplerOptions::forestColumns removed from chain.hpp/facade.hpp/sampler.hpp; tools/m4/ax_log1p_in_namespace_std.m4 and tools/build-aux/install-sh (0 bytes) deleted; DESCRIPTION's Authors@R credits pruned for code no longer present | tinytest 7475/0 and tests/cpp all 265 ok cover the surviving surface at the stack tip; the deleted tests/cpp instances are accounted for in section 4, not silently lost |
| 754309f1 | bcf-equivalence.R and multinomial-equivalence.R's ~406 shared lines (drawMatrix, drawSummary, lag1Acf, toleranceRatio, essCompare, nonFiniteParts, and compareCrossHost - the 157-line two-tier gate itself, which differed between the two files by six sprintf column-width characters) collapsed into benchmarks/R/equivalence-common.R, sourced by both; runScenarios and settingsList (each file's own scenario grid) kept separate | Both harnesses' compare mode against their pinned baselines was not re-run in this pass (a >=90-minute gate per file, section 4); the merge touches no line any recorded draw passes through - the shared functions read the .rds, they do not produce it |
| 0f5b7b0c | Cut provenance/defense prose from the remaining R comment blocks (R/bart.R, diagnostics.R, generics.R, rbart.R, utility.R) | Comment-only |
| be931ec7 | 132 landed plan docs moved to docs/plans/archive/ verbatim (255 inbound path rewrites); TODO 347 -> 316 lines (closed-entry recaps collapsed to one-line pointers, per cruft H5); four doc-vs-code falsehoods fixed (change-move-fix.md's "kept unchanged" claim about a deleted script; gate-ledger-read.md calling two wired checkers unwired; bcf-bartcause-relocation.md naming two nonexistent fixtures; a fourth, folded into the same commit) | Docs/TODO only; no code touched |
| 89577d7d | The six post-hoc engine-bool guards deleted (rulings table #5); cutpoints -> ordinalThresholds across 164 sites in the engine, both bridges, and dbarts.h's dbarts_drawLatents (the split grid's cutPoints untouched) with DBARTS_C_API_HASH re-baked to 0xb6c0e97dc0688991 and dbarts_apiSignatureToken to 0x0b33edcf638a3cd3, zero consumer hits (stan4bart and treatSens grep clean for cutpoints in a dbarts_drawLatents call); three single-caller forwarders inlined; five test-only functions renamed \*ForTesting; PredictorUpdateResult::unsupportedSource cut (rulings #3) | tinytest 7475/0 and tests/cpp all 265 ok at the stack tip; the hash re-bake is exercised by test-capi.R's stale-token block (section 4) |
| 7cd71f2d | Test suite consolidated: -631 net over five files (test-generics-errors.R, test-multinomial-surface.R, the ten-file rbart preamble family, the xbart trio, test-gp-leaves.R/test-linear-leaves.R); the rbart preamble and leaf-model scaffolding hoisted into inst/common; 20 divergent same-named helpers (seededControl, makeSampler, and others) given distinct per-file names; assertion counts held identical through the consolidation | tinytest 7475/0 (assertion count preserved by design, verified by the commit's own subject and re-confirmed at the tip in this pass) |
| the anchor pass that follows | Re-aligns every docs/design (and any remaining docs/plans) anchor this stack's edits shifted, across all ten commits at once, per house practice (composition-refusals.md sec 12: run once, last, over the stacked tree) | tools/check-doc-freshness.R, re-run after |

## 4. Coverage deltas recorded

- **tests/cpp lost one tested instance of the warm-start column-mask
  containment gate: the mean-forest one mediated by `SamplerOptions::forestColumns`.**
  `testForestColumnRestriction` (test_sampler.cpp) and `testColumnMaskContainment`
  (test_state.cpp, its comment: "the warm-start install path had no feasibility
  gate for a forest's per-forest column mask... a same-mask donor is accepted")
  both used `forestColumns` to confine forest 0 (the mean/prognostic forest) to
  a column subset and assert `setState`/`installForests` refuse an out-of-mask
  donor. `forestColumns` was a dead entrance nothing wrote in R (rulings table
  #3, yagni #15: 0 writers in `src/`, 0 in `R/`, 3 in `tests/cpp` - only the two
  deleted tests themselves) - so cutting it lost a route to the gate, not a
  live one, and no LIVE coverage is lost: the same containment mechanism
  (`columnMaskStateFeasible`, chain.hpp) is exercised through
  `ForestStructureSpec::columns` (R's `vars =` restriction) by
  `inst/tinytest/test-blocks.R`, `test-heteroscedastic-warm-start.R`, and the
  surviving `tests/cpp/test_state.cpp` interaction/block-additive tests.
- **`testVarianceRollback`'s counter genuinely weakened**, and stands as
  recorded (not restorable without `forestColumns`). Before c0ddd93e, the test
  restricted the mean forest away from columns 2-3 (via the now-cut
  `restrictMean`) so a rollback on those columns had the variance forest as
  its *sole possible objector*; the assertion was `varianceObjected > 0`, "at
  least one rollback had the variance forest as sole objector" - a claim that
  isolated the variance-forest-alone case. With `restrictMean` gone, the test
  now asserts the strictly weaker `varianceRerouted > 0`, "at least one
  rollback had variance partitions to restore" (tests/cpp/test_fuzz.cpp;
  re-run at the stack tip, "variance rollback identity (4 rolled-back
  transactions, 3 with variance partitions to restore)"). The underlying
  rollback mechanism is unchanged; the test's discriminating power over it is
  narrower.
- **A stale in-test cite, left as residue**: `inst/tinytest/test-blocks.R:311-313`'s
  comment ("the same mechanism and error text as the BCF columnMask refusal
  test in test-interactions.R") points at a test that
  `inst/tinytest/test-interactions.R` no longer carries under that
  description; comment-only, left for the next test-touching slice (section 6).
- **The kernels header-dependency gotcha is obsolete.** CLAUDE.local.md's
  standing note - "The benchmark kernel shims (benchmarks/kernels) have no
  header dependency tracking - delete their binaries after header edits" - no
  longer describes the tree: 2518a0b5 added `-MMD -MP` plus an `-include`
  of the generated `.d` file to benchmarks/kernels/Makefile (matching
  tests/cpp/Makefile's existing pattern), so a bartcore header edit is now
  picked up automatically. The note lives outside this task's edit scope
  (CLAUDE.local.md is untracked); recorded here so the next reader of that
  file knows the gotcha it names is gone.

## 5. Open doors

- **The R-visible `cutpoints` spelling vs the C++/C-API `ordinalThresholds`
  rename - VD's call.** 89577d7d renamed the internal engine identifier and
  `dbarts_drawLatents`'s parameters, disambiguating them from the split
  grid's `cutPoints`. It did not touch the R-facing vocabulary: `vars =
  c("cutpoints", "sigma", "k", "tau")` (R/diagnostics.R, `summary.bart`'s
  ordinal default) and the corresponding `man/summary.bart.Rd` prose still
  say `"cutpoints"`. Whether the R surface should follow the C-API rename, or
  whether `"cutpoints"` stays as the R-facing spelling (the ordinary
  statistical-English word, distinct from the disambiguation problem the
  C++/C rename solves) is undecided.
- **The seeds-axis re-record.** TODO's `equivalence-harness-statistical-mode`
  entry: a post-RC re-record of bcf-equivalence.R/multinomial-equivalence.R
  replacing the single-chain autocorrelated draws with independent
  per-scenario seeds. Untouched by this stack.
- **`drawShippedGlue` - FINISH-OR-CUT, unruled.** yagni #6: a duplicate
  amplitude conditional (~125 engine lines) whose own comment names its
  deletion trigger ("at which point `AmplitudeSpec::generalAmplitudeDraw`
  becomes the default and this branch is deleted"); kept today because
  cutting it re-records the bcf-equivalence baseline, so this is a baseline
  decision as much as a design one. Not ruled this round.
- **The unwired gates yagni #7 lists as FINISH, beyond the three landed in
  2518a0b5.** `benchmarks/R/backfit-exact.R` is the one `*-exact.R` script
  omitted from exact-gates.yaml's 20-gate list; `composition-matrix.R` (766
  lines, an executable oracle against feature-matrix.md) runs in no
  workflow; `tools/check-doc-freshness.R`'s Part 4 (baseline-hash
  resolvability) self-skips under CI's shallow checkout; `mutation-battery.R`
  pins a MANIFEST role=`historical` baseline; `bench-sampler.R`'s `biggrid`
  is checked before mode dispatch; `sbc.R` accepts 26 `which` values and
  `sbc.yaml` runs 6. None of these was wired or cut this round.
- **The English word "cutpoint"/"cut points" in prose is kept, deliberately
  narrower than the rename.** 89577d7d's scope was the two colliding
  *identifiers*; ordinary descriptive English (docs, Rd prose, comments using
  "cut point" as a term of art rather than as the ambiguous variable name) was
  not swept, and should not be - the collision was between two spellings of
  one C++ variable name, not between two uses of the English word.

## 6. Residue

- `makeind(x, all = )` stays a documented no-op (rulings #4): `man/makeind.Rd`
  says "not currently implemented; retained for signature compatibility",
  which is defensible today only because a `BayesTree`-migrating caller
  passing `all = FALSE` gets an error under `all = TRUE`'s semantics rather
  than silently different output - still the worst of FINISH/CUT/as-is per
  yagni #12, left for a deliberate decision rather than folded into this batch.
- The test-blocks.R stale in-test cite (section 4) and the R-visible
  `cutpoints` spelling (section 5) are the two items a reader of this record
  should not mistake for closed.
- `GP leaves` (yagni #4) and the flat C API's `getForestCalibration`/
  `setForestPriorScale` pair (yagni #5) were both explicitly KEPT rather than
  silently passed over - see rulings #10-11 for the value claims that held.
- The tour restructuring (rulings #13) is adopted as a going-forward
  structure, not executed: `docs/plans/bartcore-review-tour.md` is unchanged
  by this stack (still 975 lines, stamped ae5b91d8) and TODO's
  `second-review-followups` entry now names the regeneration as the step
  before VD's own read.
