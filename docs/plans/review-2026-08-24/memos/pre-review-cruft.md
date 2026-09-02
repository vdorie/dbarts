One of four adversarial pre-review reports read ahead of the manual review; VD's rulings on its findings are recorded in docs/plans/pre-review-cleanup.md.

# Agent accumulation and cruft review - dbarts @ bartcore 3080a9c5

Read-only. No edits made. Adversarial premise applied: every comment over three
lines and every one-caller helper had to justify itself.

## 0. Scale of the thing (context for every number below)

| layer | lines | comment lines | ratio |
|---|---|---|---|
| `R/` (25 files) | 20826 | 3788 | 18% |
| `src/bartcore/*.hpp` (11 files) | 19642 | 6169 | **31%** |
| `src/*.cpp` | 9173 | 2040 | 22% |
| `inst/tinytest/*.R` (169 files) | 41375 | 6646 | 16% |
| `man/*.Rd` (29 files) | 4100 | 99 | 2% |

- Shipped source (R + bartcore + src/*.cpp): **53,097 lines**.
- `docs/` markdown: **94,420 lines** across 194 `docs/plans` + 50 `docs/design`.
  Records are **1.78x the size of the code they describe.**
- `NEWS.Rd`: 2741 lines, **401 `\item`s**, one release.
- `TODO`: 347 lines.
- Branch: 1319 commits off `main`, of which **408 (31%) are `Record ...`
  commits** - pure records-agent output, no code.
- `docs/plans` already holds **18 whole-branch review documents** (incl. a
  `review-2026-08-24/` directory with `memos/prerc-lens{1,2,3}-*.md`, two gate
  ledgers and a consolidated report). This review is at least the third.

The single clearest quantitative signature of the process: **the record layer is
larger than the artifact**, and a third of all commits produced no shipped line.

---

# A. Comments (R layer)

## A1. 177 comment blocks over 6 lines in `R/`, 1866 lines; ~93 carry at least one sentence the comment rule forbids

Census of contiguous `#` blocks in `R/*.R`:

| file | blocks >3 | blocks >6 | prov | defend | narrate | clean |
|---|---|---|---|---|---|---|
| R/bart.R | 68 | 32 | 6 | 0 | 11 | 15 |
| R/generics.R | 61 | 27 | 8 | 4 | 7 | 8 |
| R/model.R | 35 | 20 | 3 | 1 | 4 | 12 |
| R/bartcore.R | 30 | 15 | 5 | 2 | 2 | 6 |
| R/data.R | 33 | 13 | 4 | 2 | 4 | 3 |
| R/dbarts.R | 28 | 13 | 3 | 1 | 4 | 5 |
| R/spec.R | 28 | 11 | 3 | 2 | 2 | 4 |
| R/formulaTerms.R | 12 | 8 | 1 | 1 | 0 | 6 |
| R/utility.R | 25 | 8 | 0 | 1 | 2 | 5 |
| R/mixedMatrix.R | 11 | 7 | 2 | 0 | 1 | 4 |
| others | 70 | 23 | 5 | 3 | 0 | 15 |
| **total** | **401** | **177** | **40** | **16** | **37** | **84** |

Classification is keyword-seeded and hand-verified on the 18 longest blocks;
read it as "93 of 177 long blocks contain at least one offending sentence",
not "93 blocks are wholly cruft". **1101 of the 1866 long-block lines sit in a
block with an offending sentence.**

Note the counter-example, because it matters for calibration: `R/validateComposition.R`
(landed recently) is clean - its comments state real constraints (why the rank
rule needs atom-tagging; why an exported function must restore `.Random.seed`)
and nothing else. The rule is writable. Most of `R/` was written before it, or
by agents that did not hold it.

## A2. The worst five, verbatim

**1. `[[R/model.R:1102-1138@3080a9c5]]` - 37 lines. The longest comment in `R/`. Half of it is a design derivation defending a constant to a reviewer.**

```
## count - so at the literal 1 the all-basis index prior grew as
## 1.04912 sqrt(K) s without bound, making the prior on the combined location
## depend on how the caller DECOMPOSED the mean rather than on what they said
## about it, the same model written as two forests or as four differing by
## sqrt(2). sqrt(2/K) is that invariance, at K = 2 the identity, so bcf's own
## convention is preserved rather than judged. A DEFAULT and not map algebra: a
## declared `sd` keeps its per-forest reading at every K, and the basis-free
## channel is untouched, a Cauchy having no finite variance to enter the budget
## with - so the all-basis index prior is 1.4837 s at every K while the shipped
## shape's fixed-variance channel is only CAPPED by it, rising from 0.699x of
## the classic k = 2 budget toward 0.989x rather than passing 2x at K = 10.
```
The first 7 lines (what the eight transported doubles are, in bridge order) are
a genuine constraint the code cannot show - keep. Everything from "A DEFAULT and
not map algebra" is `docs/design` material: numeric evidence that the chosen
default is right. **Fix: shorten to ~12 lines, move the derivation to the
amplitude design doc. Delta: -25.** Risk: none (comment only).

**2. `[[R/model.R:417-441@3080a9c5]]` - 25 lines, ends on an experiment result.**
```
## tail mass to within 3.4 percent (P(p < 0.01 or p > 0.99) 0.2468 against
## 0.2387, where 2 gave 0.3764).
```
Reporting a calibration experiment's numbers in a shipped file is the purest
form of "why the change is correct". **Fix: shorten to ~10. Delta: -15.**

**3. `[[R/mixedMatrix.R:547-572@3080a9c5]]` - 26 lines narrating `installPredictorColumns`'s
control flow branch by branch** ("rows = NULL replaces columns whole: a
dense-backed slot j is repointed straight at... A CSC-backed column instead has
its entries spliced... a mixed call naming both kinds dispatches per column on
the sign of the map"). This is a prose transcription of the `if` ladder directly
below it. The load-bearing sentence is one: per-observation mutation of a
CSC-backed column is refused upstream, so the partial merge only ever addresses
dense columns. **Fix: shorten to ~6. Delta: -20.**

**4. `[[R/data.R:1704-1726@3080a9c5]]` - 23 lines, justifies a magic constant *and* the
absence of a second check.**
```
# No distinct-value-count check backs this up:
# a rounding collapse to a handful of distinct doubles means the range
# spans at most a few ulps (~1e-15 relative), which the ratio already
# catches five orders of magnitude earlier, while low-cardinality data
# the ratio does NOT flag has values separated by many ulps -
```
Defending a check you did not write, to a reviewer who asked why. The
constraint worth keeping is two lines: doubles near magnitude s are spaced
~2.22e-16*s apart, so 1e-10 is ~1e6x ulp spacing. **Fix: shorten to ~8.
Delta: -15.**

**5. `[[R/bart.R:1436-1454@3080a9c5]]` and `[[R/bartcore.R:907-931@3080a9c5]]` - provenance closers.**
```
# ... so bart2's usual tree.prior/node.prior/resid.prior/control
# machinery resolves n.trees, n.chains, the tree prior and k exactly as it
# would for any other family, and one bartcore_create builds the K-forest
# engine $fit wraps: no abandoned host, no second engine.
# benchmarks/R/multinomial-equivalence.R's own scenarios still exercise the
# lower-level bartcoreMultinomialSampler path directly, which is now a thin
# wrapper over the same factory (see R/bartcore.R).
```
"no abandoned host, no second engine" answers a review comment. "still exercise
... is now a thin wrapper" is provenance. And it cites a path that is not
shipped. **Fix: delete both closers. Delta: -14.**

## A3. 31 cross-file repo-path citations inside shipped `R/` comments; 6 point at non-shipped `benchmarks/`

The rule bans `docs/` paths; the same logic bans `benchmarks/`, which is not in
the built package at all, so the reference is unresolvable for any reader of the
installed source.

- `[[R/bart.R:825@3080a9c5]]`, `[[R/bart.R:1445@3080a9c5]]`, `[[R/bart.R:1540@3080a9c5]]`, `[[R/bart.R:2312@3080a9c5]]`; `[[R/bartcore.R:954@3080a9c5]]`; `[[R/spec.R:476@3080a9c5]]`
  all cite `benchmarks/R/*.R`. **Fix: delete the citation, keep the claim.
  Delta: -8.**
- The other 25 are `R/foo.R`-style sibling citations (`R/dbarts.R` has 7,
  `R/bart.R` 11). These are defensible - they name where a rule actually lives -
  but they are also the thing that rots first. Lower priority; a symbol name
  alone (`resolveFormulaBasisSubset`) is greppable and does not rot.

Per-file counts: bart.R 11, dbarts.R 7, xbart.R 3, bartcore.R 2, rbart.R 2,
utility.R 2, A_class.R 1, data.R 1, generics.R 1, spec.R 1.

---

# B. Helpers and factoring (R layer) - three premises tested and REJECTED

I ran the adversarial premises against `R/` and three of them do not hold. I am
reporting the negative results first because they change how much of the rest
should be believed.

## B1. REJECTED: "duplicated helpers because two agents never saw each other's work"

**Zero** of 379 top-level function names in `R/` is defined twice. Not one
collision across 25 files and 1319 commits. Whatever the process cost, it did
not cost namespace collisions.

## B2. REJECTED: "dead branches left because deletion felt risky"

I found 27 functions with zero direct in-package call sites and chased every
one. **All 27 are live.** The mechanism is indirect dispatch, which defeats a
naive grep:

| function | how it is actually reached |
|---|---|
| `parsePriors` (`[[R/model.R:92@3080a9c5]]`) | `[[R/spec.R:313@3080a9c5]]` rewrites a call onto it via `quoteInNamespace` |
| `processHit` (`[[R/formulaTerms.R:298@3080a9c5]]`) | `[[R/formulaTerms.R:417@3080a9c5]]` `lapply(walked$hits, processHit, ...)` |
| `findTermInFormulaData` (`[[R/data.R:359@3080a9c5]]`) | `[[R/data.R:1133@3080a9c5]], [[R/data.R:1189@3080a9c5]], [[R/data.R:1349@3080a9c5]]` call-rewriting |
| `rethrowValidityError` (`[[R/utility.R:5@3080a9c5]]`) | `[[R/model.R:386@3080a9c5]]` / `[[R/utility.R:15@3080a9c5]]` as a `tryCatch` handler |
| `validateArgumentsInEnvironment` (`[[R/dbarts.R:253@3080a9c5]]`) | `[[R/xbart.R:80@3080a9c5]]`, `[[R/dbarts.R:598@3080a9c5]]` `quoteInNamespace` |
| `validateForestKnobs` (`[[R/model.R:985@3080a9c5]]`) | `[[R/model.R:1043@3080a9c5]]` `lapply(forests, ...)` |
| the 5 `as_draws_*` methods, `.onLoad`/`.onUnload` | registered dynamically / by R |
| `pdbart`, `pd2bart`, `makeind`, `dbartsDrawLatents`, `dbartsWorkingResponse`, `makeTestModelMatrix` | exported in `NAMESPACE`, documented in `man/` |

**No dead code in `R/`. Nothing to delete here.**

The three near-misses are `bartcoreBCFSampler` (127 lines incl. comment),
`bartcoreMultinomialSampler` (52) and `bartcoreMultinomialCountSampler` (48):
227 shipped lines with **zero production callers**, reached only from 20+
tinytest files and 5 `benchmarks/R` harnesses. That is a deliberate internal
harness API - the comment at `[[R/bartcore.R:907@3080a9c5]]` says so in as many words - so
it is not cruft. Noted and cleared.

## B3. REJECTED (mostly): "one-caller helpers are cruft until proven otherwise"

86 non-public helpers have exactly one caller. I sampled the population: the
overwhelming majority are named pipeline stages that a human would also write -
`bart2Ordinal` -> `packageOrdinalResults`, `bart2Negbin` -> `packageNegbinResults`,
`splitHurdleResponse` -> `bart2Hurdle`, the `pdbart.*` quartet, and the
`bartcoreSamplerSet*` bodies that back one R5 method each. Breaking a
2000-line dispatch into named stages is good practice, not residue.

**The `refuse*` family in particular is well-factored, contrary to the premise.**
27 `refuse*` helpers; call-site counts: `refuseUnusedGenericArgs` 27,
`refuseWithoutTrees` 10, `refuseCountsMutation` 10, `refuseResidualsSample` 6,
`refuseAmplitudeMutation` 5, `refuseFlatOffsetOnMultinomial` 5,
`refuseMatrixOffset` 4, `refusePlotTreeMethod` 4,
`refuseSurvivalProbabilitiesMethod` 4, `refuseColliding` 4, `refuseClassCiLevel` 4,
`refuseResponseFreeFormula` 4, and so on. **Only 6 of 27 have a single caller**
(`refuseBartRedirectedFamily`, `refuseHurdlePositiveMissingness`,
`refuseSamplerExtractArgs`, `refuseSparseTestReferenceAgainstTrainTypes`,
`refuseTestMissingness`, `refuseWiderTestColumns`) and four of those six are
30-44 line validators that earn a name. **Do not touch this family.** My own
first pass claimed 44 one-caller guards; that was a broken caller-counter, and
the corrected number kills the finding.

## B4. UPHELD: guard-verb drift - four names for one job

`refuse*` (27), `validate*` (14), `check*` (2), `require*` (1), `enforce*` (1).
Each slice's agent picked a verb. `refuse*` (raise, possibly conditionally) vs
`validate*` (check and return a coerced value) is a coherent split the code
mostly already honours; `checkMissingPolicy`, `checkFamilyUnsupportedArgs`,
`requireCountsCapability`, `enforceWeightPolicy` are the four that should join
one of the two. **Delta: 0 lines, -3 names.** Risk: tinytest only. Low value,
cheap.

## B5. UPHELD: 28 scattered `*UnusedArgs` / `*ForeignReasons` reason tables in `R/generics.R`

`generics.R` carries 28 separate list-valued reason tables at lines 254, 280,
284, 1042, 1051, 1172, 1226, 1412, 1503, 1551, 1710, 1791, 1836, 2087, 2097,
2109, 2117, 2118, 2128, 2129, 2140, 2141, 2152, 2153, 2154, 2287, 2459 - each
added by the slice that added its family. Lines 1987-2170 alone are 184 lines
of nothing but these tables and their accessors.

This is a real design (the surface explains every foreign argument by name), so
it is not cruft per se. The cruft is the *scattering*: three tables sit inside
the multinomial block, three inside ordinal, three inside negbin, one inside
hurdle, and the shared ones sit in a fourth place. **Fix: collect all 28 into one
contiguous table block near `refuseUnusedGenericArgs`. Delta: ~0 lines, large
readability win, removes the "where does my family's table go" question that
produced the scattering.** Risk: tinytest only.

## B6. UPHELD: the four own-class `predict.*` methods share a 20-line prologue and a 12-line epilogue verbatim

`predict.bartOrdinal` (85L) vs `predict.bartNegbin` (70L): **0.66 similarity, 51
identical lines.** The shared skeleton is:

```
type <- validateType(type, eval(formals(predict.bartX)$type))
refuseUnusedGenericArgs(list(...), "predict", "bartX", c(xUnusedArgs, ...,
  foreignArgsFor(predictForeignReasons, names(formals(predict.bartX)))))
if (is.null(object[["<family channel>.raw"]])) refuseWithoutTrees("predict")
... newdata <- validateXTest(...); n.chains <- ...; raw <- ...
```
and the epilogue (`if (length(dim(raw)) == 2L) dim(raw) <- c(dim(raw), 1L)`,
the chains collapse, `if (!is.null(ci.level)) return(posteriorInterval(...))`).

**Honest judgement: the middles are genuinely per-family and must stay parallel.
Only the prologue/epilogue should be factored** into a `predictPrologue(object,
type, dots, ...)` / `predictEpilogue(x, n.chains, combineChains, ci.level)` pair.
**Delta: -60 to -80 across the four families.** Risk: tinytest only, but it
touches every own-class predict path - run the full suite, not one file.

Similarity across the other parallel arms (`fitted` 0.36, `extract` 0.35,
`residuals` 0.51-0.58, `print` 0.54) is low enough that factoring them would be
false economy. **Do not touch those.**

## B7. UPHELD: `R/dbarts.R` 46 R5 methods hand-rolling the same try/store epilogue

`tryCatch` appears 10x, `inherits(tryResult, "error")` 6x, `storeState(ptr)` 13x,
and the normalized-duplicate pass finds the same 8-line
`tryCatch -> inherits -> stop -> updateState -> storeState -> invisible(NULL)`
window at `[[R/dbarts.R:990@3080a9c5]], 1004, 1032, 1269, 1359, 1362, 1448, 1451, 1492, 1495,
1621, 1654`. **Fix: one `withSamplerCall(ptr, updateState, expr)` wrapper.
Delta: -50.** Risk: this is the mutable-method surface `stan4bart`/`bairrtt`
drive; behavioural equivalence must hold exactly. tinytest plus the R5 surface
tests suffice - no C recompile involved.

## B8. Duplicate-window pass: 200 normalized 8-line duplicate windows in `R/`

Most are argument-forwarding lists (`[[R/bart.R:1189@3080a9c5]]/1474/1564/1790/2041` - the
same 8-arg `buildHostSamplerCall` preamble five times) and family vocabulary
vectors (`[[R/bart.R:692@3080a9c5]]`, `[[R/dbarts.R:375@3080a9c5]]`, `[[R/spec.R:799@3080a9c5]]` all spell the same
9-element `family = c("auto", "gaussian", ...)` default). The family vector is
worth centralising - three copies of a controlled vocabulary will drift.
**Fix: one `bartFamilies` constant. Delta: -20.** The argument-forwarding
duplication is inherent to R's calling convention; leave it.

## B9. REJECTED: "three names for one concept"

Checked the obvious candidates. `forestIndex` 96 uses / 0 competitors.
`keepTrees` 89 / 0. `chainNums` 22 / 0. `amplitude` 82 vs `magnitude` 7 (and the
7 are prose, not identifiers). The only real splits are
`n.test`(47)/`numTest`(13)/`nTest`(10) and
`columnNames`(58)/`predictorNames`(44)/`varNames`(29) - and both splits track a
real distinction (dotted user-facing formals vs camelCase internals; container
column names vs the design's predictor names). **Verdict: no meaningful naming
drift in `R/`. The guard-verb split (B2) is the one genuine instance.**

---

# C. Multi-agent layering (blame evidence)

## C1. Function-level layering is severe in the surface files

- `R/dbarts.R`: 172 distinct commits. `R/bart.R`: 144. `R/data.R`: 79.
  `R/generics.R`: 69 (65 distinct commits still visible in blame).
- `predict.bart` (`[[R/generics.R:289-457@3080a9c5]]`, 169 lines): **17 distinct commits**
  still authoring live lines. One line in ten comes from a different session.
- `resolveSamplerSpec` (`[[R/spec.R:102-777@3080a9c5]]`, 675 lines): **23 distinct commits.**

This is the mechanism behind A1 and B3: each session appended its clause to the
prologue and its table next to its own code, because appending is the edit that
cannot break the other slices.

---

# D. Records and process residue

## D1. The docs tree, counted three consistent ways

To avoid confusion between the figures quoted in different sections below, all
verified directly:

| measure | count |
|---|---|
| `docs/plans/*.md` (top level) | **162** |
| all `.md` under `docs/plans` (incl. the `review-2026-08-24/` subdirectory) | **194** |
| `docs/design/*.md` | **50** |
| **all `.md` lines under `docs/`** | **94,420** |
| **all tracked files under `docs/`** (adds `.R`, `.rds`, `.log`, `.txt`, `.py`) | **352** |
| all tracked lines under `docs/` | **108,720** |

Against **53,097 lines of shipped source** (`R/` + `src/bartcore` +
`src/*.cpp`): **1.78 lines of markdown record per line of code**, or 2.05
counting everything tracked under `docs/`. Section H1 works from the 162
top-level plans; section H3 handles the 32 subdirectory docs separately.

Not proposing the record system go away (settled). The ratio is the point.

## D2. 18 whole-branch review documents already exist in `docs/plans`

`architecture-numerical-review.md`, `bartcore-review-tour.md`,
`data-review-remediation.md`, `interface-review.md`,
`package-review-remediation.md`, `readability-review.md`,
`release-candidate-review.md`, `retrospective-reviews.md`,
`review-perf-followups.md`, `tau-slice-review.md`, plus the whole
`docs/plans/review-2026-08-24/` tree (`gate-ledger.md`, `gate-ledger-read.md`,
`matrix-review-entries.md`, `matrix-review-generics.md`, and three
`memos/prerc-lens*.md`).

Ten of these describe reviews whose remediations have landed. They have no
present-facing reader. **Fix: `docs/plans/archive/` (archive, not delete).
Delta: ~10-14 files moved.** Risk: none, but check for inbound links first -
`tools/check-doc-freshness.R` walks these.

## D3. WITHDRAWN: `tools/check-build-freshness.R` is wired, just not via CI

I initially flagged it as the one `tools/` script nothing runs. It is executed
by `[[benchmarks/R/mutation-battery.R:504@3080a9c5]]`. Its three siblings
(`check-doc-freshness`, `check-rc-codoc`, `check-win-drift`) are wired into
`[[.github/workflows/lint.yaml:86-92@3080a9c5]]`, and `regenerate-snapshots.R` is invoked by
five `test-reproducibility-*.R` files. **All five `tools/*.R` are live. No
action.** (The `.R` scripts should still be `.Rbuildignore`d - see H7 - but that
is about what ships, not about whether they are used.)

## D4. TODO contradicts its own stated contract in 79 of its 348 lines

`[[TODO:1@3080a9c5]]` says, verbatim: *"An unordered, forward-facing backlog; completed work
and its rationale [live elsewhere]"*. It then carries **35 entries, 10 of them
flagged CLOSED / LANDED / SHIPPED, consuming 79 lines**:

| entry | line | lines |
|---|---|---|
| `equivalence-harness-statistical-mode` | 71 | 20 |
| `binary-kforest-k1-reachability` | 46 | 13 |
| `multiforest-predictor-mutation` | 153 | 9 |
| `x86-simd` | 289 | 9 |
| `correlated-outcomes` | 63 | 8 |
| `multiforest-extension-surface` | 130 | 7 |
| `bcf-baseline-cross-host` | 33 | 4 |
| `monotone-leaf-branch-fill` | 126 | 4 |
| `r-c-division` | 185 | 3 |
| `rc-gate` | 188 | 2 |

Each keeps a full rationale paragraph explaining why it closed - which is
exactly what the header says lives in the landing note instead. **Fix: collapse
each to one line (`name: CLOSED <hash>, see docs/plans/<doc>.md`) or drop
outright. Delta: -60.** Risk: none. Highest value-per-minute item in the review.

## D5. `inst/common` holds 9 fixtures (723 lines) for 169 test files; 105 files source it

64 test files bypass the shared fixtures entirely. `bartcoreHandle.R` alone is
484 of the 723 lines - the other eight are 12-29 line data generators. The
generator set has not grown since July while the suite doubled, which is why
setup preambles are copy-pasted (see the test-suite section).

---

## D6. CLEARED: generated files, NEWS hygiene

- **No generated file is tracked.** `src/Makevars`, `src/config.hpp`,
  `src/misc/config.h`, `src/external/config.h`, `src/include/misc/types.h` are
  all untracked and gitignored, with `.in` and `.win` counterparts present.
  Clean.
- **NEWS is user-facing, not engine chatter.** Of 297 `\item`s, only 11 mention
  internal machinery, and those 11 are legitimately user-observable (the engine
  replacement, the speed delta, the removed `R_C_interface.hpp` ABI). **Zero
  near-duplicate item openings.** Size is maintainer-timed and out of scope, but
  the content is not the problem.

---

## D7. CLEARED: the mechanical house rules are being followed, without exception

This is worth stating plainly, because it bounds the whole review:

- **Zero non-ASCII characters** across `R/*.R`, `src/bartcore/*.hpp`,
  `src/*.cpp`, `src/*.hpp`, `inst/include/dbarts/*.h`, `man/*.Rd`. No en-dashes,
  no em-dashes, no arrow glyphs. Perfect compliance.
- **Zero `docs/` paths in shipped source.** The one grep hit
  (`[[src/include/external/io.h:16@3080a9c5]]`) is an open-std.org URL for a C standards
  document, i.e. a genuine citation, not a repo path.
- **Only 4 test-path references in shipped C++**, and two of them
  (`[[src/bartcore/chain.hpp:4515@3080a9c5]]`, `[[src/bartcore/chain.hpp:4648@3080a9c5]]` - "defines NDEBUG, so this is live
  only in tests/cpp") are real constraints about where an assert is live.
  The other two (`[[grow.hpp:149@3080a9c5]]`, `[[R_interface_bartcore.cpp:4167@3080a9c5]]`) are
  provenance and should go. **Delta: -2.**
- **No slice codenames in shipped source.**

So the failures found in this review are all *judgement-level* - what a comment
chooses to say, how much record gets written - never mechanical
non-compliance. That is a meaningfully better starting position than the
premise assumed.

## D8. Cross-check: the C/C++ layer carries nearly twice R's long-comment load

Same block census run over the engine and bridges:

| file | lines | blocks >3 | blocks >6 |
|---|---|---|---|
| `src/bartcore/chain.hpp` | 5513 | 178 | 75 |
| `src/R_interface_bartcore.cpp` | 7826 | 180 | 63 |
| **`src/bartcore/combiner.hpp`** | **2086** | **99** | **59** |
| `src/bartcore/model.hpp` | 4924 | 125 | 40 |
| `src/bartcore/data.hpp` | 1858 | 57 | 22 |
| `src/bartcore/sampler.hpp` | 2000 | 55 | 20 |
| `src/bartcore/tree.hpp` | 2159 | 34 | 19 |
| others | - | 93 | 30 |
| **total** | | **821** | **328** |

**328 comment blocks over six lines** against R's 177. `combiner.hpp` is the
densest file in the repo: **one 6+ line comment block every 35 lines.** Longest
individual blocks: `[[grow.hpp:125@3080a9c5]]` (53 lines), `[[scan.hpp:61@3080a9c5]]` (43),
`[[combiner.hpp:1884@3080a9c5]]` (38), `[[chain.hpp:4587@3080a9c5]]` (36).

**But read this number differently from R's.** I read the two longest engine
blocks in full. `[[grow.hpp:125-177@3080a9c5]]` documents the exact RNG draw discipline
(which node draws how many coins, in what order) - load-bearing for bitwise
equivalence and genuinely unshowable in code. `[[combiner.hpp:1884-1921@3080a9c5]]` derives
the level-centering move's exact conditional and warns that
`ConstantGaussianLeaf::scale` already carries the `1/sqrt(m_k)` so dividing
again double-counts it - a real trap, correctly documented.

In both, the rule violations are *surgical sentences inside otherwise earned
comments*:

- `[[grow.hpp:129@3080a9c5]]` "(documented so the equivalence picture stays predictable)" -
  meta-provenance.
- `[[grow.hpp:149@3080a9c5]]` "test_grow.cpp chi-squares the realized root-rule
  frequencies..." - cites a non-shipped test from shipped source.
- `[[combiner.hpp:1895@3080a9c5]]` "(a Jensen bias, confirmed against the exact gate) - so it
  is deliberately NOT used" - defends the choice to a reviewer.
- `[[combiner.hpp:1908@3080a9c5]]` "Uniform absorption is also the better mixing device
  than..." - same.

**Judgement: the engine's comment load is not where to cut.** It is derivation,
not narration, and it is the opposite of R's problem. Trim the four sentences
above and their kin; leave the mathematics alone. Cutting here has a real
downside (losing a documented draw discipline) that cutting in `R/` does not.

---

---

# E. C/C++ engine and bridges

## E1. The two-channel refusal design is the BEST-factored code in the repo, and its "defense in depth" is not defense in depth

I was asked to judge whether the predicate + raiser pairs' post-hoc checks are
cruft. **They are not, and the structure deserves the opposite of criticism.**

`src/R_interface_bartcore_common.hpp` (347 lines) is a single shared header
holding every refusal rule, and the pairing is explicit in the exports:

| predicate | raiser |
|---|---|
| `isMultiForest` | `refuseMultiForestMutation`, `refuseMultiForestResponseMutation` |
| `responseConduitIsFixed` | `refuseVarianceForestScaleUpdate`, `refuseGroupedScaleUpdate` |
| `testFitsAreUndefined` | `refuseUndefinedTestFits` |
| `sigmaIsPinned` | `refusePinnedSigmaChange` |
| `familyCarriesNoWeights` | `refuseBinaryWeightChange`, `enforceBinaryWeightPolicy` |
| - | `refuseEmptyTreeStore`, `refuseNonBinaryMask`, `refuseCscReferenceAgainstStore`, `refuseSparseLeafCovariate` |

`src/C_interface.cpp` carries **22 `using bartcore_bridge::` declarations**
(`[[combiner.hpp:23-44@3080a9c5]]`) pulling exactly these in. **There is one implementation of each rule,
and both entrances - the R bridge and the flat C API - call it.** That is not
duplication and not defense in depth; it is one source of truth with two doors.

The one rule implemented locally in `C_interface.cpp` rather than shared is
`refuseTestMissingness` (`[[combiner.hpp:256-267@3080a9c5]]`), and its comment states precisely why:

```
// A test NA takes a rule's learned missing direction, and a rule learns one
// only where the training column had NAs (ColumnStore::hasMissing gates the
// draw), so on a complete column it would take one fixed branch at every
// split. The R surface refuses first and names the column; this is the flat
// entrances' backstop, which has only the index.
```

**This is a model comment under the maintainer's rule.** It states the
constraint (why an NA cannot be routed), and the "backstop" sentence is a
*design fact* - the flat API is reachable by `LinkingTo: dbarts` consumers
without R ever running, so the C-side check is the only check on that path, and
it differs in what it can report (index, not name). Not provenance, not defense.

**Verdict: no change. The premise's "defensive checks stacked on top of checks
that already guarantee the condition" does not occur here.**

## E2. Bridge structure: 63 entry points, median 40 lines, no copy-pasted preamble pile

`src/R_interface_bartcore.cpp` at 7825 lines sounds like a per-slice dumping
ground. It is not: **63 `SEXP` entry points, median 40 lines each.** The
SEXP-extraction preamble is not duplicated per entry - the `rc.a` helpers do it
(`rc_getLength` 38 uses, `rc_getDouble` 20, `rc_getInt` 12, `rc_getBool` 6
across the whole file), which is far too few uses for 63 hand-rolled preambles.
The 314 `Rf_error` calls are the bridge's validation surface, which is where R-
facing argument checking belongs.

## E3. Comment load is derivation, not narration - see D8

Full census and the four surgical violations are in D8 above. Summary: 328
comment blocks over six lines across the engine and bridges (vs 177 in `R/`),
densest in `combiner.hpp` (one per 35 lines), but the content is exact
conditionals and RNG draw disciplines that are load-bearing for the bitwise
equivalence gate. **This is the one comment corpus in the repo I would not cut.**

## E4. Shipped C++ hygiene: 4 test-path references, 2 of them legitimate

Detailed in D7. `[[src/bartcore/chain.hpp:4515@3080a9c5]]` and `[[src/bartcore/chain.hpp:4648@3080a9c5]]` state a real
constraint about `NDEBUG` making an assert live only under `tests/cpp`.
`[[src/bartcore/grow.hpp:149@3080a9c5]]` ("test_grow.cpp chi-squares the realized root-rule
frequencies...") and `[[src/R_interface_bartcore.cpp:4167@3080a9c5]]` are provenance and
should go. **Delta: -2 lines.**

## E5. CORRECTED: one genuinely dead function in the engine, and my own clearance of it was wrong

**I got this wrong the first time and am correcting it in place.** I originally
cleared all 36 zero-call candidates. A second census, run with clang
(`-Wunused-function -Wunused-member-function -Wunused-private-field
-Wunreachable-code-aggressive`) rather than grep, found one survivor:

**`GPGaussianLeaf::maxLeafSize()` at `[[src/bartcore/model.hpp:1352@3080a9c5]]` is dead.**
```
  std::size_t maxLeafSize() const { return maxLeafSize_; }
```
I had cleared it against `[[src/R_interface_bartcore.cpp:1356@3080a9c5]]`. That line is
```
    model.gpMaxLeafSize = static_cast<size_t>(maxLeafSize);
```
where `maxLeafSize` is a **local `int` declared three lines above** by
`rc_getInt`, not a call to the accessor. Every other hit in `model.hpp` is
either the member field `maxLeafSize_` (`[[src/R_interface_bartcore.cpp:1455@3080a9c5]], [[src/R_interface_bartcore.cpp:1523@3080a9c5]], [[src/R_interface_bartcore.cpp:1597@3080a9c5]], [[src/R_interface_bartcore.cpp:1663@3080a9c5]], [[src/R_interface_bartcore.cpp:2087@3080a9c5]]`), a
constructor parameter (`[[src/R_interface_bartcore.cpp:1360@3080a9c5]], [[src/R_interface_bartcore.cpp:1363@3080a9c5]]`), or a comment (`[[src/R_interface_bartcore.cpp:1324@3080a9c5]], [[src/R_interface_bartcore.cpp:1447@3080a9c5]]`). The
accessor has zero callers. **Fix: delete. Delta: -1.** Risk: none.

**Methodology note that caused my error, and that matters for anyone repeating
this work:** a local, gitignored agent-worktree directory holds a **full stale
copy of the source tree** (complete with `chain.hpp`, `model.hpp`,
`combiner.hpp`) inside a build subdirectory. Any `grep -r` from the repo root
manufactures phantom callers out of it. It is `.gitignore`d and
`.Rbuildignore`d, so it neither ships nor is tracked - but it is a trap for
exactly this kind of audit, and well over a thousand source-like files sit
under that directory in total. **Confirm every "it has a caller" claim against
an explicit path list, never against `.`.**

The rest of the section stands: the other 35 candidates are alive.

## E5b. QUALIFIED: the test-only surface is *mostly* disciplined - five functions break the convention

My "a reader can never mistake it for production surface" was too strong. Beyond
the 16 `*ForTesting` symbols, **five test-only functions carry no marker at all**
and read as production API:

| function | declared | only callers |
|---|---|---|
| `Chain::interweaveGlueRidge` | `[[chain.hpp:1331@3080a9c5]]` | `[[tests/cpp/test_sampler.cpp:2730@3080a9c5]]` |
| `TResponse::estimatesResidualDf` | `[[model.hpp:4167@3080a9c5]]` | `tests/cpp/` |
| `NBResponse::estimatesDispersion` | `[[model.hpp:4463@3080a9c5]]` | `tests/cpp/` |
| `ColumnStore::testColumnIsSparse` | `[[data.hpp:1807@3080a9c5]]` | `[[tests/cpp/test_model.cpp:5767@3080a9c5]]` |
| `ColumnStore::testSparseColumn` | `[[data.hpp:1810@3080a9c5]]` | `[[tests/cpp/test_model.cpp:5768@3080a9c5]]` |

I actually observed three of these myself and recorded them as "alive" without
drawing the conclusion - being reached *only* from `tests/cpp` is the finding,
not the clearance. **Fix: rename to `*ForTesting` so the convention is
exhaustive. Delta: 0 lines, 5 renames.** Risk: `tests/cpp` only.

## E5c. Original census (the 35 that are alive)

I ran a function-definition census over `src/bartcore/*.hpp` plus the bridges:
**638 function-like definitions, 36 with zero apparent call sites.** I chased
the non-obvious ones. **All 36 are alive.** The detector's misses were
template instantiation, member calls, and function-pointer use:

| candidate | actually reached by |
|---|---|
| `bcfGlue` (`[[combiner.hpp:820@3080a9c5]]`) | 23 mentions; live |
| `blendSoftmax` (`[[combiner.hpp:2043@3080a9c5]]`) | called at `[[combiner.hpp:1879@3080a9c5]]` |
| `interweaveGlueRidge` | `[[tests/cpp/test_sampler.cpp:2730@3080a9c5]]` |
| `partitionByPredicate` (`[[tree.hpp:730@3080a9c5]]`) | `[[tree.hpp:761@3080a9c5]]`, `[[tree.hpp:771@3080a9c5]]` (lambda arguments) |
| `mapCutPointsBelow` (`[[tree.hpp:1284@3080a9c5]]`) | recursive self-call at `[[tree.hpp:1294@3080a9c5]]` |
| `runTestFitRange` (`[[chain.hpp:4051@3080a9c5]]`) | function pointer to `misc_mt_runTasks`, `[[chain.hpp:4083@3080a9c5]]` |
| `maxLeafSize` | `[[src/R_interface_bartcore.cpp:1356@3080a9c5]]` |
| `testColumnIsSparse` | `[[tests/cpp/test_model.cpp:5767-5768@3080a9c5]]` |
| `monotoneIntegrate` (`[[model.hpp:366@3080a9c5]]`) | one caller, `[[model.hpp:842@3080a9c5]]` - and it names a mathematical operation, so it earns the name |

**13 of the 36 are `*ForTesting` accessors** - a deliberate, self-documenting
test-only surface. It is **16 distinct symbols, 24 occurrences**, split
`chain.hpp` 11 / `model.hpp` 12 / `sampler.hpp` 1, and **every one of the 16 is
consumed by `tests/cpp`** (set difference is empty). Small and fully used - but
see E5b: the convention is not exhaustive, five more test-only functions carry
no marker.

## E6. Engine clone duplication is near zero - 25 windows in 19,642 lines - but two sites are real

Normalized 10-line duplicate-window pass over `src/bartcore/*.hpp`: **25
windows, collapsing to 9 distinct sites, every one a pair.** For scale, `R/`
produced 200 windows at an 8-line window in a comparable 20,826 lines. **The
engine is the least clone-duplicated code in the repository.** The nine sites:
`[[tree.hpp:1322@3080a9c5]]/1399`, `[[model.hpp:3922@3080a9c5]]/4175`, `[[model.hpp:3117@3080a9c5]]/3256`, `[[model.hpp:1610@3080a9c5]]/2007`,
`[[model.hpp:1175@3080a9c5]]/1192`, `[[model.hpp:1054@3080a9c5]]/1419`, `[[model.hpp:1013@3080a9c5]]/1381`, `[[chain.hpp:5146@3080a9c5]]/5210`, `[[chain.hpp:1418@3080a9c5]]/1990`.

Two deserve names:

**(a) `[[model.hpp:1013@3080a9c5]]` vs `[[model.hpp:1381@3080a9c5]]` - a genuine factorable duplicate, ~13 lines x2.**
The covariate-standardization setup is copied verbatim between two leaf models,
and the *only* difference is which cache the first line clears:
```
    clearStatisticsCache();        <- model.hpp:1013
    clearKernelCaches();           <- model.hpp:1381
    std::size_t numColumns = numCovariates_;
    numObservations_ = data.numObservations;
    means_.assign(numColumns, 0.0);
    sds_.assign(numColumns, 1.0);
    u_.resize(numObservations_ * numColumns);
    for (std::size_t j = 0; j < numColumns; ++j) {
      const double* column = data.rawColumn(columns_[j]);
      double mean, sd;
      if (!data.suppliedStandardization(columns_[j], &mean, &sd))
        standardizationMomentsForColumn(column, numObservations_, &mean, &sd);
      means_[j] = mean;
      sds_[j] = sd;
```
This is the one place in the engine where two agents plainly wrote the same
block twice. **Fix: a shared `standardizeCovariates(data)` in the common base;
each leaf model keeps only its own cache-clear. Delta: -13.**
**Risk: it touches standardization, which feeds every draw. Gate: the bitwise
equivalence harness, not tinytest.** Small payoff against a real gate - do it
only when something else is already opening these files.

**(b) `[[chain.hpp:1418@3080a9c5]]` vs `[[chain.hpp:1990@3080a9c5]]` - parallel backfit skeletons. LEAVE ALONE.**
Both iterate `forests_`, bind `forestY`/`forestWeights` and branch on
`combiner_`, but they diverge immediately: `[[chain.hpp:1990@3080a9c5]]` calls
`combiner_->drawForestGlue(f, rng_, forests_)` before forming the response and
`[[chain.hpp:1418@3080a9c5]]` does not - **a different RNG consumption order.** Factoring these would
put a draw-order difference behind a shared function's parameter, which is
exactly how a bitwise gate breaks quietly. **Legitimately parallel. No change.**

## E7. Engine naming: mostly disciplined, but `cutPoints` and `cutpoints` are TWO DIFFERENT CONCEPTS one capital letter apart

I closed the naming question directly. Most axes are clean singletons:
`numTrees` 188/0, `forestIndex` 78/0, `numCategories` 34/0. Two apparent splits
are legitimate: `numObs` (145) is a local shorthand for the long accessor
(`std::size_t numObs = node.numObservations();`, `[[model.hpp:1257@3080a9c5]]`), not a
competing API name; `workingResponse` (25) vs `resid` (52) are the GLM working
response and the plain residual, genuinely different things.

**But this one is real, and it is worse than drift:**

| identifier | uses | what it is | declared |
|---|---|---|---|
| `cutPoints` | **61** | the per-column **split candidate grid** for tree rules | `[[data.hpp:527@3080a9c5]]` `std::vector<std::vector<double>> cutPoints;` |
| `cutpoints` | **42** | the ordinal cumulative-probit's **K-1 latent thresholds** | `[[combiner.hpp:80@3080a9c5]]` `std::vector<double> cutpoints;` |

These are unrelated quantities. Both are plural, both are `vector<double>` or
`const double*`, they live in files that include one another, and
**`cutpoints` is in the shipped public header** (`[[inst/include/dbarts/dbarts.h:1272@3080a9c5]]`,
"cutpoints the ordinal's numCutpoints strictly increasing interior..."), so a
`LinkingTo: dbarts` consumer meets it. `grep -i cutpoints` returns both, and
nothing but a capital P separates the split grid from an ordered probit's
thresholds.

Distribution is by file, which is what has kept it survivable so far:
`cutPoints` in `data.hpp`/`tree.hpp`/`sampler.hpp`; `cutpoints` in
`combiner.hpp`/`chain.hpp`/`model.hpp`/`facade.hpp` and both bridges. They
overlap in `sampler.hpp`, `R_interface_bartcore.cpp` and
`R_interface_bartcore.hpp`.

**Fix: rename the engine-internal one. `cutpoints` -> `ordinalThresholds`
(42 sites) is the safer direction than touching `cutPoints`' 61, but note the
public header spells it `cutpoints` in `dbarts_*` signatures - so either accept
that the C API name differs from the internal one, or take the ABI change now
while nothing is frozen. Delta: ~0 lines, one concept disambiguated.**
**Risk: mechanical rename, but it crosses the shipped header, so it is an ABI
decision, not a cleanup. Gate: `tests/cpp` plus a `LinkingTo` consumer rebuild
(stan4bart).** This did not appear in any prior review tour.

*(Coverage caveat: a deeper adversarial sweep of engine internals was dispatched
and had not returned when this report was assembled. E1-E7 are my own direct
evidence and now cover every question that was asked, naming included.)*

## E10. NOT CRUFT - A REAL BUG. `dbarts_sampler_setWeights` installs negative gaussian case weights

**This is the only actual defect in the review and it should be fixed before
anything else here is touched.** I verified it independently and it is
unambiguous, because the codebase states the violated contract in its own words
and honours it everywhere else.

**The contract**, `[[src/R_interface_bartcore_common.hpp:182-183@3080a9c5]]`, verbatim:
```
/// gaussian passes through, validated for non-negativity by its callers.
```
So `enforceBinaryWeightPolicy` checks probit (refuse) and logistic (positive
integers) and **deliberately delegates the gaussian non-negativity check to each
caller.** There are four callers. Three honour it. One does not:

| caller | non-negativity check? |
|---|---|
| `[[R_interface_bartcore.cpp:2525@3080a9c5]]` (creation) | via the R surface ahead of it |
| `[[R_interface_bartcore.cpp:4735@3080a9c5]]` (setData) | via the R surface ahead of it |
| `[[R_interface_bartcore.cpp:4950@3080a9c5]]` (R bridge setWeights) | **YES** - `[[R_interface_bartcore.cpp:4954@3080a9c5]]` `if (!(weights[i] >= 0.0)) Rf_error("weights must be non-negative");` |
| **`[[C_interface.cpp:702@3080a9c5]]`** (`dbarts_sampler_setWeights`) | **NO** |

The flat entrance runs the family policy and then installs directly:
```
  enforceBinaryWeightPolicy(weightShape.family, weights,
                            weightShape.numObservations);
  samplerOf(sampler).setWeights(weights);      // C_interface.cpp:704
```

**Why this is a bug and not a judgement call:** the flat C API has no R layer in
front of it - that is the whole point of it, and the file says so at
`[[R_interface_bartcore.cpp:4730-4732@3080a9c5]]` ("the R surface checks it first, so this is
a no-op there and the real gate for the flat C API"). And the *same file* gets
it right 400 lines later for a different weight vector:

```
// C_interface.cpp:1111-1115, dbarts_sampler_setForestWeights
  if (weights != NULL)
    ...
      if (!R_FINITE(weights[i]) || weights[i] < 0.0)
        Rf_error("dbarts_sampler_setForestWeights: weights must be finite and "
                 ...
```
So forest weights are guarded, case weights are not, in one file, by one author.

**Impact.** A gaussian case weight is a precision multiplier (`sd_i = sigma /
sqrt(w_i)`). A negative one is a negative precision contribution to leaf
sufficient statistics: it can yield a negative posterior variance, `NaN` draws,
or - worst - a finite but silently wrong fit with no error raised anywhere. The
exposed consumers are exactly the `LinkingTo: dbarts` packages the flat API
exists for (stan4bart, bairrtt), which reach it without R validation.
`R_FINITE` is unchecked too, so `NA`/`Inf` weights install as well.

**Fix:** three lines at `[[C_interface.cpp:703@3080a9c5]]`, mirroring `[[C_interface.cpp:1111-1115@3080a9c5]]`:
```
  for (size_t i = 0; i < weightShape.numObservations; ++i)
    if (!R_FINITE(weights[i]) || weights[i] < 0.0)
      Rf_error("dbarts_sampler_setWeights: weights must be finite and "
               "non-negative");
```
**Delta: +4. Gate: `tests/cpp` plus a new `test-capi.R` case; no equivalence
run needed (it only adds a refusal on input that was previously undefined
behaviour).** Nothing about the two-channel design changes - this is a missing
caller-side check, not a design question.

**One item I could not fully close:** the flat *creation* path's weight
validation at `[[C_interface.cpp:1249-1254@3080a9c5]]` gates on
`readsWeights = in.weights != NULL && resolved == RF::logistic`, so gaussian
creation weights appear to be unchecked there too. I did not trace every flat
creation entrance, so treat creation as *likely* exposed and confirm it while
fixing the mutation path.

## E8. Eight single-caller trivial forwarders in the engine, -28 lines

Distinct from B3's R-side finding, and this time the premise holds. Of 193
single-caller names in the engine, eight are pure forwarders that name nothing
the caller did not already know. Verified individually:

| forwarder | declared | its one caller | why it is cruft |
|---|---|---|---|
| `ColumnStore::clearTest()` | `[[data.hpp:1789@3080a9c5]]` | `[[sampler.hpp:1366@3080a9c5]]` | a rename of the already-public `resetTestStorage()` |
| `testColumnIsCscBacked()` | `[[data.hpp:1043@3080a9c5]]` | `[[data.hpp:1052@3080a9c5]]` | consumed **nine lines below its own definition** |
| `ColumnSource::isRank()` | `[[data.hpp:158@3080a9c5]]` | `[[data.hpp:469@3080a9c5]]` | a second named layer over one enum compare |
| `treesSplittingOnColumn()` | `[[chain.hpp:2247@3080a9c5]]` | `[[sampler.hpp:1888@3080a9c5]]` | one-line partial application of public `collectSplittingTrees` |
| `varianceTreesSplittingOnColumn()` | `[[chain.hpp:2271@3080a9c5]]` | `[[sampler.hpp:1892@3080a9c5]]` | same, and its call site is **four lines from the previous one** |

The last pair is the clearest case: two wrappers whose only callers sit four
lines apart in the same function, both forwarding to one public method.
**Fix: inline all eight. Delta: -28.** Risk: `tests/cpp` plus a normal build.

**One exclusion I endorse:** `gaussianPdf` (`[[model.hpp:340@3080a9c5]]`) is also a
single-caller helper, but it sits on the monotone-leaf marginal, so inlining it
moves floating-point evaluation order. **Leave it** - or gate the change on the
equivalence harness against a monotone baseline. Not worth it for one line.

## E9. A genuine style fork for the maintainer to settle: six post-hoc `Rf_error` guards

`[[src/R_interface_bartcore.cpp:3874@3080a9c5]], [[src/R_interface_bartcore.cpp:3922@3080a9c5]], [[src/R_interface_bartcore.cpp:3962@3080a9c5]], [[src/R_interface_bartcore.cpp:4010@3080a9c5]], [[src/R_interface_bartcore.cpp:4059@3080a9c5]], [[src/R_interface_bartcore.cpp:4094@3080a9c5]]` each read

```
if (!engine->set*(...)) Rf_error(...);
```

where the engine's only `false` arm is exactly the condition the bridge already
probed immediately above. **The two reviewers of this code disagree, and the
disagreement is legitimate**, so I am presenting it rather than resolving it:

- **Delete them (-42).** They are dead by construction; the file's own comment
  at `[[src/R_interface_bartcore.cpp:4164-4167@3080a9c5]]` argues precisely this way about a sibling case.
- **Keep them (-6, cutting the three verbatim-triplicated comments to one line
  each).** Deleting the test means discarding a `bool` the engine deliberately
  returns, and clang cannot prove the branches unreachable.

The resolution is a house rule, not a review finding: *"never discard a `bool`
returned by the engine"* keeps them; *"no unreachable branches"* deletes them.
**Both positions agree on one real defect regardless of which rule wins:**
`[[src/R_interface_bartcore.cpp:4092-4095@3080a9c5]]` carries three lines of user-facing prose about fractional weights
that **cannot be reached**, because `refuseNonBinaryMask` raises first. That is
a broken error message, not a style question. **Fix it either way. -3.**

Note this is a *different* question from E1. E1 concerns the shared
predicate/raiser pairs in `R_interface_bartcore_common.hpp`, which are one
implementation serving two entrances and are not cruft. E9 concerns bridge-local
re-tests of an engine return value. My E1 verdict is unaffected.

---

# F. Test suite (169 files, 44,963 lines, 5432 assertions)

**Headline: the suite is in far better shape than the premise predicts.**
Corpus-wide clone redundancy is **6.1% of normalized code lines** (2164 of
35,223 at an 8-line window; only 1.3% at a 15-line window, meaning the
duplication is short knob-blocks, not cloned test bodies). Recoverable: ~850
lines, **1.9% of the suite**.

**And the four largest files are the four cleanest**: `test-capi.R` (2012
lines, 270 assertions) has **zero** redundant normalized lines;
`test-bartcore.R` (1105) zero; `test-bcf-creation.R` (1255) 16 of 1020;
`test-sparse-factor.R` (1118) 26 of 886. The adversarial instinct to go after
the big files is wrong here.

## F1. REJECTED: "tests that assert a value no code path can change"

Across all 5432 assertions: **zero** `expect_true(TRUE)`, **zero**
`expect_equal(<literal>, <same literal>)`, **zero** assertions on R-language
guarantees dressed as package behaviour. Total genuinely dead assertions found:
**five**, worth ~8 lines.

1. `[[inst/tinytest/test-generics-multithreaded.R:270-273@3080a9c5]]` - the one true
   self-comparison:
   ```r
   expect_identical(
     predict(rbartFit, xRe, group.by = g, n.threads = 2L),
     predict(rbartFit, xRe, group.by = g, n.threads = 2L)
   )
   ```
   The comment says the intent is "`group.by` follows `...`, so it is matched
   by name only". **That intent is untestable this way** - if `group.by` were
   swallowed by `...`, both sides would be equally wrong and it would still
   pass. The 1-thread-vs-2-thread assertion at `[[inst/tinytest/test-generics-multithreaded.R:260-263@3080a9c5]]` already covers it.
   **Delete. -4.**
2. `[[inst/tinytest/test-argument-surface.R:853@3080a9c5]]`
   `expect_equal(length(formals(dbarts::xbart)), 32L)` - fails on every
   legitimate signature change, names nothing, and the `for (knob in
   xbartKnobs)` loop at `[[inst/tinytest/test-argument-surface.R:845-850@3080a9c5]]` already asserts the surface by name.
   **Delete. -1.**
3. `[[inst/tinytest/test-summary-nondefault-families.R:31@3080a9c5]], 44, 56` - three
   `expect_false(identical(class(...), "summary.default"))` that cannot fail
   while the `expect_true(inherits(..., "summary.bart"))` two lines below
   passes. Also calls `summary()` a second time for nothing. **Delete. -3.**

## F2. REJECTED: forbidden references in tests - zero hits

Grepped the whole suite plus `inst/common` for `docs/` paths, `.md` filenames,
PR numbers, commit hashes and slice codenames: **0 hits.** Every occurrence of
the word "slice" is an R variable (`sweepSlice`, `savedSlice`) or English about
array indexing. One prose wobble only: `[[test-generics-intervals.R:152-154@3080a9c5]]` says
"already had ci.level third before this slice" - a temporal self-reference, the
single provenance leak in 45k lines.

## F3. UPHELD: 811 section banners in three competing styles

| style | count | files |
|---|---|---|
| `# --- text ---` | 483 | 56 |
| `## text` | 269 | 15 |
| `## --- text ---` | 59 | 7 |
| `# ------` (bare rule) | 3 | 1 |

**The drift is between agents, not within them**: only 3 files mix styles
(`test-convergence-diagnostics.R`, `test-heteroscedastic-warm-start.R`,
`test-mutate-sparse-valued.R`). Each slice's agent picked a style and held it.
Sub-drift in numbering: `## --- 1. ...`, `# --- (1) ...`, `# --- family gate 2:`
and roman `# --- (i)/(ii)` (`test-forest-basis-subset.R`) all coexist.
**Fix: normalize on the 60%-share style opportunistically. Delta 0. Risk none.**
Not worth a dedicated pass.

## F4. UPHELD, and the most important test finding: 25 helper names defined in more than one file, 20 with DIFFERENT bodies

This is the real accretion damage, and it costs comprehension rather than disk:

```
seededControl   7 defs / 7 files, 7 DISTINCT shapes (6 functions + 1 plain list)
                test-bcf-creation.R:30           (n.trees=50, seed=17, n.chains=2)
                test-bcf-family.R:32             (25 / 29 / 1)
                test-bcf-r5-surface.R:24         (40 / 23 / 1)
                test-forest-basis-r5.R:31        (25 / 41 / 1)
                test-forest-basis-subset.R:29    (5 / 7 / 1)
                test-forest-weights-r5.R:22      (30 / 71 / 1)
                test-summary-nondefault-families.R:13  <- a list(), not callable
makeSampler     6 defs / 5 files, 6 DISTINCT bodies
logisticSampler 3 defs / 3 files, 3 distinct
earlyRMSE       3 defs / 3 files, 2 distinct
control         3 defs / 3 files, 3 distinct  (and `control` as a function name
                                               is itself a hazard)
recordChannels / buildSampler / parityArm / refusalArm / roundTrip / handleOf
                3 defs each, all distinct
```

tinytest runs each file in a fresh environment, so there is **no runtime
collision** - the cost is that a reader who learns what `seededControl()` means
in one file learns something false about the next six. Five helpers
(`countWarnings`, `splitTrees`, `recombine`, `parseLeaves`,
`leafOfObservations`) have identical bodies and should simply move to
`inst/common` (**-52 lines, zero risk**). The six `seededControl` functions
differ only in three numeric literals and unify trivially - **but the per-file
seed and `n.trees` literals must survive verbatim at each call site or RNG draws
move.**

## F5. UPHELD: the rbart group-by preamble, character-identical in 10 files

`test-rbart-{bartcore,error,example,generics,groupby,multithreaded,options,reproducibility,weights}.R:1-27`
and `[[test-reproducibility-rbart.R:6-27@3080a9c5]]` are byte-for-byte identical (modulo one
extra `source(captureWarnings.R)` in two of them): the Friedman source, the
`RNGkind(sample.kind = "Rounding")` dance, the group draw, `sigma.b <- 1.5`,
the `b[g]` shift and the `rm()`.

**Fix: `inst/common/rbartGroupData.R`. Delta: -190.**
**Risk: this is the one merge needing a gate beyond "tinytest passes."** The
block advances the RNG stream and `test-reproducibility-rbart.R` carries 5
hardcoded snapshot values downstream of it. `source(local = TRUE)` runs the same
calls in the same order so the stream is preserved - but re-run both
`*reproducibility*rbart*` files specifically to confirm.

`inst/common` uptake is otherwise good: `friedmanData.R` 50 files,
`bartcoreHandle.R` 27, `stateContinuation.R` 13, `captureWarnings.R` 12.
Only `pdData.R` (17 lines) has a single consumer.

## F6. CLEARED: hardcoded snapshots are contained and correctly shaped

Only **8 files** carry a 4+ decimal literal, **84 numbers total**, and **zero
lines carry three or more** (no jammed inline vectors). The six
`test-reproducibility-*.R` files each open with a header naming the regeneration
path and the review workflow: *"After an intentional shift, run
tools/regenerate-snapshots.R and eyeball that the new values move by a plausible
magnitude."* A compact digest would be shorter but would destroy exactly that
workflow. **No change recommended.**

## F7. The five worst files, with merge operations

| file(s) | now | after | delta | operation |
|---|---|---|---|---|
| `test-generics-errors.R` | 755 | ~635 | **-120** | Three arms (predict `:1-410`, extract `[[inst/tinytest/test-generics-errors.R:412-581@3080a9c5]]`, fitted `[[inst/tinytest/test-generics-errors.R:583-755@3080a9c5]]`) each rebuild the same six fits with byte-identical knobs. `bart2FitKT:113` and `bartFitKT:415` are identical including `keepTrees=TRUE`. The dataset is built three times unseeded (`xSmall:225`, `xSmall2:482`, `xSmall3:587`) though every assertion is an error-message check that never reads it. Build once at the top. Also 11 fewer sampler runs at test time. |
| `test-multinomial-surface.R` | 983 | ~830 | **-150** | Highest absolute redundancy in the corpus (168/781 = 22%). `internalMultinomialFit:30-59` and `internalMultinomialCountFit:403-432` are 30-line near-clones - **the file's own comment at `[[inst/tinytest/test-multinomial-surface.R:400@3080a9c5]]` admits it**. Plus the same 8-line knob block 10x (`[[inst/tinytest/test-multinomial-surface.R:89@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:117@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:165@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:459@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:498@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:551@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:654@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:680@3080a9c5]], [[inst/tinytest/test-multinomial-surface.R:702@3080a9c5]]`). |
| rbart preamble family (10 files) | ~230 | ~65 | **-190** | See F5. Worst single file `test-rbart-groupby.R` (44% redundant): the preamble plus three near-identical surface arms at `[[inst/tinytest/test-rbart-groupby.R:170-200@3080a9c5]],[[inst/tinytest/test-rbart-groupby.R:235-270@3080a9c5]],[[inst/tinytest/test-rbart-groupby.R:285-315@3080a9c5]]` differing only in `combineChains` and a 0.90-vs-0.80 bound. |
| xbart trio | 437 | ~300 | **-137** | `test-xbart-reproducibility.R` is **68% redundant** - four copies of the same 12-line `xbart(seed=0L)` call varying only `n.threads` (1,1,2,2). Two runs prove determinism; four copies do not. Plus a shared `checkXvalShape()` across `test-xbart-method.R`/`test-xbart-loss.R`. |
| `test-gp-leaves.R` + `test-linear-leaves.R` | 644 | ~560 | **-84** | Three byte-identical blocks modulo `node.prior` and one tolerance: state round-trip (`gp:229-264` = `linear:182-217`), the plotTree/prior-swap refusal (`gp:123-143` = `linear:73-94`), the `rm()` teardown (`gp:290-307` = `linear:263-280`). Per-family assertions stay put - that part is legitimately parallel. |

**Bonus trap found:** `test-reproducibility-xbart.R` and
`test-xbart-reproducibility.R` are different files with near-anagram names (a
snapshot tripwire vs a seed-determinism check). Same for
`test-reproducibility-rbart.R` / `test-rbart-reproducibility.R`. Two agents, two
naming conventions, no coordination.

## F8. Legitimately parallel - do NOT touch

- Cross-family generic-refusal blocks (`[[test-hurdle.R:289-310@3080a9c5]]`,
  `[[test-nbinom.R:476-511@3080a9c5]]`, `[[test-ordinal.R:495-524@3080a9c5]]`): each family must refuse
  independently. The near-duplication *is* the test.
- `[[test-rng.R:84@3080a9c5]], [[test-rng.R:95@3080a9c5]], [[test-rng.R:109@3080a9c5]], [[test-rng.R:120@3080a9c5]]`: four seeded `bart()` calls whose whole point is
  that identical arguments give identical draws.
- The five worst in-file "repeat" offenders (`test-fits-without-offset.R` 28,
  `test-data-testOffset.R` 26, `test-makeModelMatrix.R` 16) all apply the
  repeated predicate to a **rebuilt** object. Real per-construction checks.
- `test-mutate-then-serialize.R` (0 assertions) delegates to `statesAgree()` in
  `inst/common/stateContinuation.R`, which asserts ~30x internally. 13 files
  delegate this way; the density metric, not the file, is misleading.

---

# G. The "duplicate by design" ruling on the equivalence harnesses - OVERTURN IT

I was asked to judge whether today's ruling (that the two equivalence harnesses
duplicate ~250 lines of helpers *by design*) should stand. **It should not, and
the pair identified in the ruling is probably the wrong pair.**

## G1. The duplicating pair is `multinomial-equivalence.R` / `bcf-equivalence.R`, and the overlap is 406 lines, not 250

Normalized-diff over the five largest `benchmarks/R` harnesses:

| pair | shared runs >=5 lines | shared lines |
|---|---|---|
| **`multinomial-equivalence.R` / `bcf-equivalence.R`** | **21** | **406** |
| `equivalence.R` / `multinomial-equivalence.R` | 3 | 26 |
| `equivalence.R` / `bcf-equivalence.R` | 2 | 21 |
| `equivalence.R` / `sbc.R` | 0 | 0 |
| `equivalence.R` / `composition-matrix.R` | 0 | 0 |

`equivalence.R` (1913 lines) is genuinely independent of everything - 26 shared
lines at most. The duplication is entirely between the multinomial and bcf
harnesses, which share **11 of their 11-and-12 top-level function names**:
`compareCrossHost`, `drawMatrix`, `drawSummary`, `essCompare`, `lag1Acf`,
`makeControl`, `nonFiniteParts`, `recordChannels`, `runScenarios`,
`settingsList`. These two files are the same harness with different scenarios.

## G2. The design argument does not survive contact with what is actually shared

The design argument would be one of two things: *(a) a bitwise-equivalence
oracle must not share code with what it validates*, or *(b) the harnesses must
stay independently runnable at different commits*. I read the shared functions.
**Neither holds.**

Every shared helper is a **statistic over draws**, not oracle logic:

- `lag1Acf` (`[[multinomial-equivalence.R:153@3080a9c5]]` = `[[bcf-equivalence.R:170@3080a9c5]]`, 12
  lines): row-wise lag-1 autocorrelation. Pure arithmetic on a matrix.
- `essCompare` (`[[bcf-equivalence.R:183@3080a9c5]]` = `[[bcf-equivalence.R:200@3080a9c5]]`, 13 lines): a two-sample z-statistic with an
  ESS correction. Reads only the two arms' draw matrices.
- `drawSummary` (`[[bcf-equivalence.R:146@3080a9c5]]` = `[[bcf-equivalence.R:163@3080a9c5]]`, 4 lines): `rowMeans` / `var` / `ncol`.
- `nonFiniteParts` (`[[bcf-equivalence.R:367@3080a9c5]]` = `[[bcf-equivalence.R:384@3080a9c5]]`, 24 lines): `sprintf` formatting of
  non-finite cell counts.
- `makeControl` (`[[bcf-equivalence.R:95@3080a9c5]]` = `[[bcf-equivalence.R:113@3080a9c5]]`, 8 lines): a `dbartsControl(...)` builder.

**None of these reimplements sampler behaviour, so none can launder a bug into
the verdict.** They are applied *symmetrically to both arms*: a bug in
`lag1Acf` degrades the test's power identically on both sides whether the
function is shared or copied. A *copied* bug is strictly worse, because it has
to be found and fixed twice, and the second copy is the one everyone forgets.

Argument (b) fails too: both files live in the same repo at the same commit, and
`benchmarks/baselines/*` are keyed by commit hash - checking out an old baseline
checks out the old harness, and would equally check out an old shared helper
file. The independence is not actually purchased by the duplication.

**Fix: `benchmarks/R/equivalence-helpers.R`, sourced by both. Delta: -180 to
-200** (406 shared lines collapse to one ~200-line file plus two `source()`
lines). **Risk: real but bounded.** These harnesses gate bitwise equivalence, so
the merge must be verified by running both against their existing baselines and
confirming identical verdicts - not merely that they run. That is the gate, and
it is a gate beyond tinytest.

## G3. Independent second pass, same verdict, with a sharper fact

A second sweep measured the same pair by whole-file LCS rather than 8-line
windows and found **587 of `bcf-equivalence.R`'s 909 lines matching (65%)**,
498 of them in contiguous runs of 6+ lines (largest: `bcf:143-246`, 104 lines).
Helper-by-helper:

| helper | bcf | mn | identical? |
|---|---|---|---|
| `drawMatrix` | 4 | 4 | YES |
| `drawSummary` | 4 | 4 | YES |
| `lag1Acf` | 12 | 12 | YES |
| `toleranceRatio` | 6 | 6 | YES |
| `essCompare` | 13 | 13 | YES |
| `nonFiniteParts` | 24 | 24 | YES |
| **`compareCrossHost`** | **157** | **157** | **151 of 157** |
| `makeControl` | 8 | 8 | 7 of 8 |
| `settingsList` | 11 | 10 | 9 of 11 |
| `recordChannels` | 13 | 17 | 6 of 13 |
| `runScenarios` | 281 | 377 | 35 of 281 (genuinely per-model) |

**The decisive fact: `compareCrossHost` - the 157-line two-tier cross-host
verdict, which IS the gate logic - differs between the two files in exactly six
lines, and all six are a `sprintf` column width (`"%-14s"` vs `"%-6s"`).** That
is a copy with a formatting edit, not an independent implementation.

And argument (a) is **already violated by these same two files**:
`[[bcf-equivalence.R:56-59@3080a9c5]]` and `[[multinomial-equivalence.R:66-69@3080a9c5]]` both
`source(system.file("common", "bartcoreHandle.R", package = "dbarts"))` - 484
lines of handle API loaded **from the installed package under test**. If that is
acceptable (and it is; 27 tinytest files and 14 benchmark scripts do it), then a
sibling `.R` in `benchmarks/R/` is *strictly safer* - it is not built by the
compiler whose output is being validated. `[[equivalence.R:39-42@3080a9c5]]` already sources
a sibling (`bartcore-shim.R`), so the pattern is established.

The duplication also provides **zero** independent-implementation value: sharing
`compareCrossHost` to 151/157 lines means a defect in the tier-1 bound is
present in both, while guaranteeing a fix to one is silently absent from the
other.

**Fix: `benchmarks/R/equivalence-common.R` holding the six byte-identical
helpers plus `compareCrossHost` (column width as a parameter). Delta: -220 net.
Do NOT merge `runScenarios` (281 vs 377, only 35 shared) or `settingsList` -
those encode each model's scenario grid and must stay separate.**
**Risk: LOW and quantified.** No recorded draw passes through any of the seven
helpers - they read the `.rds`, they do not produce it. Gate: run both compares
against the current pinned baselines before and after; both must stay bitwise.

---

# H. Records, scaffolding and CI

**Scale first:** `docs/` is 352 files and **108,720 lines** - more than
`R/` (22,071) + engine and bridge (29,306) + tinytest (46,469) + `man/` (4,445)
= 102,291, **combined**.

## H1. 118 of 162 plan docs (32,992 lines, 54% of the tree) are residue with no present-facing reader

Criterion applied: a doc has a present-facing reader if its INDEX status is
open/reference/NO-GO/research, or the root TODO names it. Otherwise it is a plan
for landed work that its own landing note and the code now supersede.

Largest: `bart2-argument-consolidation.md` 1965, `dbarts-h-reshape.md` 1863,
`bcf-bartcause-relocation.md` 1693, `latent-subset-mask.md` 1228,
`multiforest-veto-rate-falsifier.md` 1185, `grow-from-root-categorical-scan.md`
1130, `rd-records.md` 967, `c-api-growth.md` 964, `bcf-public-surface.md` 959.

**Fix: `git mv` into `docs/plans/archive/`. Delta: 118 files out of the
top-level listing, 0 lines deleted.** **Risk: link rot, and it is the whole
cost** - 49 of the 118 are cited by a `docs/design/*.md` "Plan:" footer, 10 by
the nav docs, 43 by a still-facing plan: ~102 mechanical path rewrites plus 118
INDEX row moves. Nothing else breaks.

## H2. FALSIFIED: "an INDEX with purpose text longer than some docs"

`docs/plans/INDEX.md` is 298 lines, 159 rows. Total per-entry purpose text is
**29,262 chars against 3,586,566 chars of indexed docs - a ratio of 0.008.**
Worst single row is `grouped-equivalence.md` at 0.109 (180 chars indexing a
1,644-char doc). **No row approaches its target's size. The INDEX is not the
problem; leave it.**

Real INDEX defects instead, both small:
- `docs/plans/gate-hardening-1.0.md` (269 lines) has **no INDEX row**. 162 files,
  159 rows, and the header claims "161 files" - off by two in both directions.
- `docs/plans/review-2026-08-24/` (138 files) is entirely un-indexed, yet 11
  tracked files cite into it by path. **Fix: +6 lines. Risk: none.**

## H3. `docs/plans/review-2026-08-24/` is the sharpest residue in the repo

138 tracked files, 2.2 MB: 37 `.R`, 32 `.md` (8,337 lines), 18 `.txt`, 18
`.log`, 17 `.rds`, 7 `.py`, 3 `.jsonl`, 3 `.csv`, 2 `.out`, 1 `.json`. **The
`.rds`/`.log`/`.out`/`.jsonl` alone are 476 KB of binary and raw run logs checked
into a docs tree.** Its own scripts still point at their pre-promotion home -
`[[matrix.R:4@3080a9c5]]`, `[[matrix-entries.R:7@3080a9c5]], [[matrix-entries.R:18@3080a9c5]]`, `[[matrix-sampler-census.R:18@3080a9c5]]`,
`[[wave3-plan.md:3@3080a9c5]], [[wave3-plan.md:670@3080a9c5]]` all cite their old scratch location, dangling for
anyone who clones.

**Fix: keep the 32 `.md` (TODO cites their ledger); move the 37 `.R` + 7 `.py`
to `benchmarks/R/` or delete; drop the 476 KB of regenerable logs/rds; rewrite
the 5 dangling scratch paths. Delta: -40 files, -476 KB, ~-6,000 lines. Risk: low.**

## H4. Records template drift: five status dialects, 115 heading signatures, 20 landing-note spellings, seven codename alphabets

- **Front matter:** 85 plans carry no status marker at all; 42 use `## Status`,
  22 `Status:`, 12 `status:`, 1 `STATUS:`. The 85 unmarked force
  `tools/check-doc-freshness.R` and the INDEX to infer status from prose -
  which INDEX.md's own header admits.
- **115 distinct H2 heading signatures across 162 docs.** The dominant template
  accounts for only 22 docs across its 4 spellings.
- **20 spellings of the landing-note heading**: `## Landing notes` 21,
  `## Landing note (DATE, HASH)` 20, `## Landing` 14, `## Landings` 12, plus 12
  per-slice variants and singletons like `## Landing record`.
- **Landing notes live in both trees**: 87 of 162 plans and 9 of 50 design docs
  carry one, with no stated rule deciding which.
- **Seven concurrent codename alphabets**: `S0..S14` (29 docs), `C1..C7` (27),
  `D0..D4` (15), `B0..B3` (15), `M0..M4.5` (14), `step 1..n` (15), plus
  `Stage A/B`, `Phase 1/2`, `wave-3`, `FA5`, `FB16`, `P16/P17`. Several docs mix
  two in one file (`bcf-bartcause-relocation.md` uses D0/D2/D3 *and* B0/B0b/B1).

**Fix: one status line on all 162 (+85 lines); one landing-note spelling and a
rule that the note lives in the plan (9 section moves); one line in
`docs/plans/README.md` naming `S<n>` as the going-forward default. No retro-rename
- hashes and NEWS quote the old labels.** Do the landing-note pass *before*
H1's archive or the notes archive with the plans.

## H5. TODO: -55 lines now, -79 once the RC review closes

347 lines, 34 items. 13 carry a closed/landed marker; **6 are closed with no
open door at all** (pure recap the plan doc and NEWS already carry):
`augmentation-helpers` 5, `cheap-uniformity` 4 (*it literally says "No open door
recorded"*), `monotone-leaf-branch-fill` 4, `multiforest-extension-surface` 7,
`multiforest-predictor-mutation` 9, `r-c-division` 3. **Collapse: -26.**

Plus: landed-recap prefixes before the first `Open:` on items that *do* have a
door (`equivalence-harness-statistical-mode` 15 of 20 lines,
`second-review-followups` 9 of 14) **-21**; `bcf-baseline-cross-host` and
`equivalence-harness-statistical-mode` describe the **same** landing (3f532af2)
and name the same remaining door - merge, **-4**; `multiforest-mutation-gaps`
spends 4 of 16 lines on a `monotoneTreeIsFeasible` paragraph concluding
"deliberately left", a settled non-decision that belongs in the code comment,
**-4**. `release-candidate-review` is 38 lines (14% of the item budget) whose
findings ledger already lives in the review directory - **-24 deferred until it
closes.** Risk: none; every collapsed line is reproduced in the named plan doc.

## H6. CLEARED, emphatically: NEWS content

2,741 lines, 400 items, 24 sections; 1.0-0 is 2,182 lines / 222 items (80% of
the file). Size is maintainer-timed and out of scope. **Content is the strongest
part of the record system:**

- **Duplicates: essentially none.** Pairwise Jaccard over 218 parsed items finds
  exactly one pair >= 0.45, and it is by design (the UPGRADING preamble says
  "each is detailed under the other headings").
- **Unobservable internal entries: none.** Every entry is written from an R- or
  C-API-visible surface.
- **No slice codenames anywhere.**

Three real defects only: `[[NEWS.Rd:1925@3f532af2]], [[NEWS.Rd:1941@3f532af2]], [[NEWS.Rd:1956@3f532af2]]` name unexported internals
(`bartcore*` functions, `bartcore_create`) - rewrite to the observable effect;
BUG FIXES item #6 (18 lines) opens "That conditioning is now applied per
forest..." as a continuation of #5 promoted to a standalone `\item`, unreadable
out of order. Flagged only (per the settled call): **43 of the 118 BUG FIXES
entries (368 lines) have as their subject a capability that ships new in 1.0-0**
(multinomial, heteroscedastic, BCF/`forests`, ordinal, nbinom, hazard, hurdle) -
no 0.9-34 user saw those behaviours, so they are not fixes.

## H7. Scaffolding: two dead files, and `tools/` ships to CRAN unnecessarily

**Delete (2 files, -39 lines):**
- `tools/m4/ax_log1p_in_namespace_std.m4` (39 lines) - not `m4_include`d
  (`[[configure.ac:13-16@3f532af2]]` lists exactly four), `AX_CXX_LOG1P_IN_NAMESPACE_STD` has
  0 hits in generated `configure`, `HAVE_LOG1P_IN_NAMESPACE_STD` 0 hits
  anywhere. Last touched 2021-10-04.
- `tools/build-aux/install-sh` - **0 bytes**, `configure:2636` sets
  `ac_aux_files="config.guess config.sub"`, no `AC_PROG_INSTALL`. Last touched
  2014. Risk: low; re-run `autoreconf -i && ./configure` after.

**`tools/` is not `.Rbuildignore`d**, so 2,350 lines of developer R ship in the
CRAN tarball (`check-doc-freshness.R` alone is 1,771 lines / 61 KB).
`build-aux/` and `m4/` must ship; the five `.R` scripts need not.
**Fix: `^tools/[^/]*\.R$` in `.Rbuildignore`. Delta: -2,350 shipped lines.
Risk: none.** Best line-per-keystroke item in the review.

**Everything else in `tools/`, `inst/common`, `tests/cpp` is live.** My earlier
`check-build-freshness.R` suspicion resolves: it is executed by
`[[benchmarks/R/mutation-battery.R:504@3f532af2]]`, so it is wired, just not via CI.
`tests/cpp` is **the most current thing audited** - newest file touched in the
same commit that reshaped `inst/include/dbarts/dbarts.h`; 13 `test_*.cpp` <-> 13
suites in `[[main.cpp:73-90@3f532af2]]` <-> 13 Makefile SOURCES entries, no orphans, zero
references to a removed API.

**`benchmarks/kernels` header-tracking gap confirmed**: its Makefile has no
`-MMD`, no `-MP`, no `-include *.d` (contrast `tests/cpp/Makefile:14` and `[[main.cpp:40@3f532af2]]`,
which have all three), and its only header prerequisite is `partition_u8.h` at
`Makefile:18`. `linear_leaf.cpp` includes the header-only
`<bartcore/bartcore.hpp>`, so an engine header edit changes no prerequisite.
**Fix: 3 lines copied from `tests/cpp`. Removes the documented "delete the
binaries by hand" footgun.**

## H8. FALSE ALARM, corrected on verification: the baseline tag IS pushed

The records sweep flagged as HIGH risk that 9 of 11 baselines - including the
three **current** gate baselines `6e3b9fb8` (bcf), `4d9a3337` (multinomial) and
`ab1dc52` (speed) - are reachable only from the local tag
`bartcore-pre-cran-rebase`, and that a fresh clone therefore could not verify
them.

**I checked `git ls-remote --tags origin` directly. The tag is on origin:**
```
465188178bba641e42589fa4443f460cb2d3c49c  refs/tags/bartcore-pre-cran-rebase
fc2af5ce17f5b3c4858fbb0220c91de9251aac9b  refs/tags/bartcore-pre-cran-rebase^{}
```
`git clone` fetches tags by default, so every baseline hash is reachable from a
fresh clone. **No action. Cleared.** The only residual is that the baselines hang
off a tag rather than a branch, so deleting that tag would orphan them - worth
knowing, not worth doing anything about.

Related, and this one **does** stand: the MANIFEST names **7 baselines with no
file present**
(`equivalence-8b047f8b`, `-4a42620a`, `multinomial-equivalence-1027be5`,
`-b354f3a`, `-19c499a`, `-fbd2168`, `-7903855`), six of them cited as neutrality
evidence for baselines that *are* present - so the evidence chain terminates in
files nobody can open. **Fix: mark each `(recording not retained)`. +7 clauses.**

## H9. NOT CRUFT - a real defect: the doc-freshness CI gate cannot see its own inputs

`.github/workflows/lint.yaml` job `doc-freshness` (lines 74-92) runs
`tools/check-doc-freshness.R`, whose declared inputs are `docs/design/*.md`,
`docs/design/INDEX.md`, `docs/plans/*.md` and `TODO`. **All three are in the
workflow's own `paths-ignore`** (`[[lint.yaml:11@3080a9c5]]` `docs/**`, `[[lint.yaml:12@3080a9c5]]` `TODO`, `[[lint.yaml:13@3080a9c5]]`
`**.md`).

**A push touching only `docs/` or `TODO` never runs the doc-freshness check -
which is every records commit, the common case on this branch** (352 tracked
files under `docs/`, 408 `Record ...` commits).

**Fix: split lines 74-92 into `.github/workflows/doc-freshness.yaml` with no
`paths-ignore`; it needs no package build. Delta: -19 from lint.yaml, +25 new.**
Do **not** instead drop `docs/**` from lint.yaml - that makes every docs push pay
the lint job's full `local::.` compile.

## H10. CI: 54 duplicated lines, half the gate surface automation-dead, 12 stale branch names

- **Duplicate jobs.** `[[equivalence.yaml:69-122@3080a9c5]]` - jobs `bcf-equivalence` (26
  lines) and `multinomial-equivalence` (28 lines) - run the same script, same
  pinned baseline hash and same `--cross-host` flag as the single 24-line step in
  `[[exact-gates.yaml:141-164@3080a9c5]]`, **which fires on every push**. The two are
  byte-identical to each other except 2 of 27 lines. Each baseline hash is pinned
  twice (4 literals for 2 files). **Fix: delete `[[equivalence.yaml:69-122@3080a9c5]]`; keep
  the gaussian job (35-67), which alone uses `--strict-coverage`. Delta: -54
  lines, 3 jobs -> 1, 4 hash pins -> 2. Risk: low - per-push coverage is strictly
  higher than the weekly cron it replaces.**
- **Half the gate surface has never run automatically.** 10 of 11 workflows exist
  only on `bartcore`; `origin/main` carries `check-standard.yaml` alone. GitHub
  binds `schedule` to the **default branch**, so all 5 crons are inert -
  equivalence, rchk, revdep-smoke, sbc, valgrind. **14 of 28 runner jobs (50%),
  including every memory-safety, calibration and revdep gate.** No YAML fix, but
  the green badge covers half of what is declared. Validate all 5 the day
  bartcore reaches main.
- **`master` referenced 12x** across 6 workflows (push + PR triggers) - no such
  branch exists. `[[check-standard.yaml:5@3080a9c5]], [[check-standard.yaml:14@3080a9c5]]`, `[[cpp-tests.yaml:9@3080a9c5]], [[cpp-tests.yaml:18@3080a9c5]]`,
  `[[exact-gates.yaml:38@3080a9c5]], [[exact-gates.yaml:45@3080a9c5]]`, `[[lint.yaml:7@3080a9c5]], [[lint.yaml:16@3080a9c5]]`, `[[pkgdown.yaml:6@3080a9c5]], [[pkgdown.yaml:15@3080a9c5]]`,
  `[[sanitizers.yaml:17@3080a9c5]], [[sanitizers.yaml:26@3080a9c5]]`. **12 edits, zero risk.**
- **`[[valgrind.yaml:41@3080a9c5]]` installs only `c("tinytest","Matrix")`** where
  `[[sanitizers.yaml:64@3080a9c5]]` installs `c("tinytest","Matrix","survival","posterior")`.
  37 tinytest files reference `posterior` and 17 reference `survival`, all behind
  `requireNamespace()` guards - **so under valgrind those blocks skip silently and
  still report green**, and valgrind carries none of `[[sanitizers.yaml:87-125@3080a9c5]]`'s
  count floors to catch the shrinkage. One-line fix.
- **PR/push asymmetry:** 5 of 6 push workflows apply `paths-ignore` to `push` but
  not `pull_request`, so a docs-only PR runs the full 14-job / 12-compile matrix
  while the identical push does not.
- **Two display names are still filenames** (`[[lint.yaml:23@3080a9c5]] name: lint.yaml`,
  `[[pkgdown.yaml:20@3080a9c5]]`). `pkgdown.yaml` is the only push workflow without
  `cancel-in-progress`, and its job-level concurrency groups all non-PR runs under
  the constant `pkgdown-true`, so a push burst **queues** 30-minute site builds.

**CLEARED:** no dead path globs (all 7 `paths-ignore` blocks resolve against
tracked files), **no dead matrix arms** (every matrix variable consumed and
proven), no commented-out YAML, no outdated action pins, and all 20 gate scripts
/ 3 baselines / 3 tools scripts / 6 tinytest files / `tests/cpp/Makefile`
referenced by CI exist - zero broken references.

## H11. CLEARED: generated vs source, and `.win` consistency

All five generated files present, untracked, `.gitignore`d (lines 9, 12-15); all
five `.in` sources tracked. **No generated file is checked in.** The five `.win`
variants are consistent with their `.in` counterparts - the 29 "absent"
`config.hpp` symbols are correct (every consumer uses `#ifdef`, so omission ==
disabled, and the `COMPILER_SUPPORTS_*` set arrives via `src/misc/config.h.win`
plus `-D` flags from `[[src/Makevars.win:13@3080a9c5]]`).

Two gaps worth closing: **`tools/check-win-drift.R` checks the 4 config header
pairs but not `src/Makevars.win`** - and that is the one pair with a real
structural divergence (`AVX_FLAG`/`NEON_FLAG`/`COMPILER` -> `ARCH`).
`PACKAGE_VERSION`/`PACKAGE_STRING` are hardcoded `"1.0-0"` in
`src/config.hpp.win:7-9`, a hand-maintained duplicate of DESCRIPTION - correct
today, silently wrong the first version bump that forgets them.

## H12. Four doc-vs-code falsehoods found in passing

- `[[docs/plans/change-move-fix.md:267@3080a9c5]]` - says `change-fix-stage2.R` "and its
  result files are kept unchanged as the experiment record"; **the `.R` was
  deleted 2026-08-24 in 7318b266.** (Keep the 4 result files in
  `benchmarks/results/` - they are unreproducible evidence for a landed engine
  decision - but fix the sentence.)
- `[[docs/plans/review-2026-08-24/gate-ledger-read.md:103@758bccdd]], [[docs/plans/review-2026-08-24/gate-ledger-read.md:121@758bccdd]]` - calls
  `check-win-drift.R` and `check-doc-freshness.R` **unwired**; both are wired
  (`[[lint.yaml:86@3080a9c5]]`, `[[lint.yaml:89@3080a9c5]]`).
- `[[docs/plans/bcf-bartcause-relocation.md:1224@3080a9c5]]` - names
  `inst/common/linearData.R` and `binaryData.R`; **neither exists.**
- `DESCRIPTION:19-21` - credits `ax_cxx_namespace_std.m4` and
  `ax_func_posix_memalign.m4`; **no such files in the tree.**

## H13. FALSE ALARM, corrected on verification: gh-pages is clean

The records sweep flagged a possible disclosure: `docs-site/` (gitignored
pkgdown output) does contain `CLAUDE.local.html` and `CLAUDE.local.md` - the
private instructions rendered into the site tree by a local pkgdown run.

**I listed `origin/gh-pages` directly: 204 files, zero matching `claude`
case-insensitively. The deployed site does not carry them, and cannot** - CI
builds the site from a clean checkout, and `CLAUDE.local.md` is not checked in,
so the file only exists in a local working tree. **No action. Cleared.**

## H14b. Verified as stated

Two claims I re-checked myself because they drive slate items, both confirmed:

- **H9 (doc-freshness gate)**: `[[lint.yaml:10-14@3080a9c5]]` `paths-ignore` is literally
  `docs/**`, `TODO`, `**.md`, `benchmarks/baselines/**`, and the `doc-freshness`
  job at `[[.github/workflows/lint.yaml:74-92@758bccdd]]` runs `check-win-drift.R`, `check-doc-freshness.R` and
  `check-rc-codoc.R`. A docs-only or TODO-only push skips the whole workflow.
  Confirmed. (The `pull_request` trigger has no `paths-ignore`, so the gate does
  fire on PRs - which softens it only if records changes go through PRs. On a
  branch with 408 direct `Record ...` commits, they do not.)
- **`master` in 6 workflows**: `git branch -a | grep -i master` returns nothing.
  There is no `master` branch, local or remote. Confirmed dead.

---

# I. Prioritized cleanup slate

## Item 0 - NOT a cleanup. Fix this first, on its own.

**`dbarts_sampler_setWeights` (`[[src/C_interface.cpp:702-704@3080a9c5]]`) installs negative,
`NA` and `Inf` gaussian case weights without a check**, violating a contract the
codebase states in its own header (`[[R_interface_bartcore_common.hpp:182@3080a9c5]]`:
"gaussian passes through, validated for non-negativity by its callers") and
honours at its other three callers, including the sibling
`dbarts_sampler_setForestWeights` 400 lines below in the same file. A negative
case weight is a negative precision contribution and can produce a silently
wrong fit. Exposed consumers are the `LinkingTo` packages the flat API exists
for. **+4 lines. Gate: `tests/cpp` + a `test-capi.R` case. No equivalence run
needed.** Full detail and the patch in E10. Also confirm the flat *creation*
path (`[[src/C_interface.cpp:1249-1254@758bccdd]]` gates on logistic only).

## The ten highest-value cleanups

Ordered by value, not by size. "Gate" means what must pass beyond `tinytest`.

| # | item | where | delta | gate |
|---|---|---|---|---|
| 1 | **Split `doc-freshness` into its own workflow with no `paths-ignore`.** The gate's declared inputs (`docs/**`, `TODO`, `**.md`) are all in the workflow's own ignore list, so it never runs on the commits it exists to check - and 408 of 1319 branch commits are exactly those. | `[[.github/workflows/lint.yaml:10-14@3080a9c5]], [[.github/workflows/lint.yaml:74-92@758bccdd]]` | -19 / +25 | CI: confirm the new workflow fires on a docs-only push |
| 2 | **Add `^tools/[^/]*\.R$` to `.Rbuildignore`.** 2350 lines of developer R (incl. `check-doc-freshness.R` at 1771 lines / 61 KB) currently ship in the CRAN tarball. `build-aux/` and `m4/` must stay. | `.Rbuildignore` | **-2350 shipped** | `R CMD build` + `R CMD check` |
| 3 | **Merge the two equivalence harnesses' shared helpers.** `compareCrossHost` - the 157-line gate logic itself - differs between them in six `sprintf` column widths. Both files already source 484 lines from the installed package under test, so the "oracle independence" argument is already spent. | `benchmarks/R/{bcf,multinomial}-equivalence.R` -> new `equivalence-common.R` | **-220** | **run both compares against the pinned baselines before and after; verdicts must stay bitwise** |
| 4 | **Archive 118 landed plan docs (32,992 lines) and strip 476 KB of logs/rds from `docs/plans/review-2026-08-24/`.** `git mv` to `docs/plans/archive/`, nothing deleted. | `docs/plans/` | 118 files moved, -40 files, ~-6000 lines | none; but ~102 inbound path rewrites + 118 INDEX row moves are the whole cost |
| 5 | **Delete the duplicate CI equivalence jobs.** `[[equivalence.yaml:69-122@3080a9c5]]` runs the same script, baseline hash and flag as `[[.github/workflows/exact-gates.yaml:141-164@758bccdd]]`, which already fires on every push. The two jobs differ from each other in 2 of 27 lines. | `.github/workflows/equivalence.yaml` | -54, 3 jobs -> 1, 4 hash pins -> 2 | CI green on next push |
| 6 | **Consolidate the test suite's five worst files.** `test-generics-errors.R` rebuilds the same six fits in three arms (-120, and 11 fewer sampler runs); `test-multinomial-surface.R` (-150); the rbart preamble to `inst/common` (-190); the xbart trio (-137); gp/linear leaf scaffolding (-84). | `inst/tinytest/` | **-680** | **re-run both `*reproducibility*rbart*` files - the shared preamble advances the RNG stream feeding 5 snapshot values** |
| 7 | **Cut the R comment corpus's provenance and design-derivation.** 93 of 177 blocks over six lines carry an offending sentence; 1101 of 1866 long-block lines sit in one. Start with `[[R/model.R:1102@3080a9c5]]` (-25), `[[R/model.R:417@3080a9c5]]` (-15), `[[R/mixedMatrix.R:547@3080a9c5]]` (-20), `[[R/data.R:1704@758bccdd]]` (-15), the 6 `benchmarks/` citations (-8), the 2 shipped-C++ test-path cites (-2). | `R/`, `src/` | -250 to -400 | none (comments only) |
| 8 | **Collapse the TODO's closed entries.** 6 items are closed with no open door at all; another 4 carry landed-recap prefixes before their first `Open:`. Every collapsed line is reproduced verbatim in the named plan doc. The file's own line 1 says completed work lives elsewhere. | `TODO` | -55 (-79 after the RC review closes) | none |
| 9 | **Fix the four doc-vs-code falsehoods and the two dead build files.** `[[docs/plans/change-move-fix.md:267@3080a9c5]]` names a deleted `.R` as "kept"; `[[docs/plans/review-2026-08-24/gate-ledger-read.md:103@758bccdd]], [[docs/plans/review-2026-08-24/gate-ledger-read.md:121@758bccdd]]` calls two wired checkers unwired; `[[docs/plans/bcf-bartcause-relocation.md:1224@3080a9c5]]` names two nonexistent fixtures; `DESCRIPTION:19-21` credits two absent `.m4`. Delete `tools/m4/ax_log1p_in_namespace_std.m4` (39 lines, unreferenced since 2021) and the 0-byte `tools/build-aux/install-sh`. | mixed | -39, 6 edits | `autoreconf -i && ./configure` after the m4 delete |
| 10 | **Rename the 20 divergent test helper names, hoist the 5 identical ones.** `seededControl` means seven different things in seven files (one of them a plain `list()`, not callable); `makeSampler` six. Zero runtime risk - tinytest isolates each file - pure comprehension cost on every future read. | `inst/tinytest/`, `inst/common/` | -52 (+renames) | per-file seed and `n.trees` literals must survive verbatim at each call site |

**Engine items added after the late sweep returned** (all verified
independently): delete the one dead accessor `GPGaussianLeaf::maxLeafSize()`
(`[[model.hpp:1352@3080a9c5]]`, -1); inline the eight single-caller forwarders (-28, E8);
rename the five convention-breaking test-only functions to `*ForTesting`
(E5b, 0 lines); fix the unreachable fractional-weight error prose at
`[[R_interface_bartcore.cpp:4092-4095@3080a9c5]]` (-3). And settle the E9 style fork - six
post-hoc `if (!engine->set*(...)) Rf_error(...)` guards at
`[[R_interface_bartcore.cpp:3874@3080a9c5]], [[R_interface_bartcore.cpp:3922@3080a9c5]], [[R_interface_bartcore.cpp:3962@3080a9c5]], [[R_interface_bartcore.cpp:4010@3080a9c5]], [[R_interface_bartcore.cpp:4059@3080a9c5]], [[R_interface_bartcore.cpp:4094@3080a9c5]]`: **-42 if the
house rule is "no unreachable branches", -6 if it is "never discard a `bool`
from the engine".** Two reviewers split on it; it is a house call, not a finding.

**Not cheap, but decide it before the release candidate:** `cutPoints` (61 uses,
the split candidate grid) and `cutpoints` (42 uses, the ordinal's K-1 latent
thresholds) are two unrelated `vector<double>` concepts separated by one capital
letter, and the second is in the shipped `dbarts.h`. Renaming the internal one
to `ordinalThresholds` costs ~0 lines but crosses the public header, so it is an
ABI call while nothing is frozen - not a cleanup. See E7.

**Runners-up, cheap and worth doing opportunistically:** `master` -> `[main,
bartcore]` in 6 workflows (12 edits, no such branch exists); add `-MMD -MP
-include $(DEPS)` to `benchmarks/kernels/Makefile` (3 lines, kills the documented
"delete the binaries by hand" footgun); add `survival`+`posterior` to
`[[valgrind.yaml:41@3080a9c5]]` (37 tinytest files silently skip under valgrind today and
still report green); mark the MANIFEST's 7 missing baselines
`(recording not retained)`; standardise the four stray guard verbs
(`check*`/`require*`/`enforce*`) onto `refuse*`/`validate*`; factor the
`predict.bart*` prologue/epilogue across the four own-class families (-60 to
-80); one `withSamplerCall()` wrapper for `R/dbarts.R`'s repeated
try/storeState epilogue (-50).

**Aggregate if all ten land: roughly -3,900 lines shipped or tracked (-4,000 with
the engine items), plus 118 files archived and 40 deleted** - one real bug fixed
(item 0), and two gate defects closed (item 1, and via item 5 the double-pinning
of baseline hashes).

---

# J. Judgement: how much accumulation does this branch carry?

**Less than the premise assumes in the code, and far more than it assumes in the
records.** Of the eight failure modes I was asked to hunt, four are simply
absent: there are no duplicated helper names in `R/` (zero collisions in 379
definitions across 1319 commits), essentially no dead code (all 27 zero-caller R
functions turned out to be reached by `lapply`, `tryCatch` or `quoteInNamespace`
rewriting, and 35 of 36 zero-call C++ candidates by templates, member calls or
function pointers - **one real dead accessor in 53,097 lines**, and I had
wrongly cleared even that one), no meaningful naming drift (`forestIndex` 96/0,
`keepTrees` 89/0),
and essentially no tautological tests (five dead assertions in 5,432). The
mechanical house rules are followed without exception - zero non-ASCII
characters and zero `docs/` paths across every shipped file. The two-channel
refusal design is not merely acceptable but the best-factored code in the
repository: one shared header, one implementation per rule, 22 `using`
declarations pulling it into the flat C API, and the single locally-implemented
rule carries a comment explaining exactly why it must be. The test suite's
clone redundancy is 6.1%, and its four largest files carry approximately none of
it. Where the process did leave marks in the code, they are the marks of
*appending*: 28 reason tables scattered through `generics.R` because each slice
put its own next to its own code; a 20-line prologue repeated across four
parallel `predict` methods; `seededControl` meaning seven different things in
seven test files; 483 section banners in three competing styles. That is real,
it is worth roughly 1,500 lines, and none of it is dangerous. The genuine
accumulation is elsewhere and is not code at all: **108,720 lines of `docs/` -
more than `R/`, the engine, the tests and `man/` combined - of which 118 plan
documents totalling 32,992 lines have no present-facing reader; 408 of 1319
commits that produced no shipped line; 18 whole-branch review documents, this
being the third or later; seven concurrent slice-codename alphabets and 20
spellings of "Landing note".** And the one place that sprawl has actually cost
something is the sharpest finding in the review: the CI job that checks
documentation freshness has its own inputs in its own `paths-ignore`, so on a
branch where a third of the commits are records commits, **the records gate has
never run on a records commit.** The code is in better shape than a reader of
the commit log would guess; the record system that describes it has outgrown
both its readers and its own guard rails. One last thing worth saying plainly:
this review was commissioned to hunt cruft, and the single most valuable thing
it turned up is not cruft at all but a **bug** - `dbarts_sampler_setWeights`
installing negative gaussian case weights against a contract the codebase writes
down and honours everywhere else. It survived because it sits on the one path
with no R layer in front of it, which is precisely the path a cruft review is
not looking at. Two of my own confident clearances also turned out to be wrong
(a dead accessor I cleared against a same-named local variable; a test-only
surface I called exhaustive when five functions break its convention), and both
were caught only because a second reviewer ran a compiler where I had run grep.
Take that as the calibration note for the numbers above: the negative results in
this report are strong on the axes that grep can settle, and should be trusted
less on the axes that need a build.
