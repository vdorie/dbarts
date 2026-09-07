# retire-grouped-random-effects

agent: opus (engine + bridge), sonnet (R, docs, tests)
rng: shifting - no surviving fit changes a draw, but the equivalence
     baseline loses two scenarios and re-records
budget: ~3500 deleted / ~200 added, across roughly 100 files

## Goal

Grouped random intercepts are gone from dbarts. `rbart_vi()`, its nine
methods, `R/sliceSample.R`, the `GroupedResponse` decorator and its tau
block, the bridge that reaches them, and the two `dbarts_results` fields no
longer exist; stan4bart is the documented home for multilevel structure. The
per-sweep `setOffset`/`setResponse` conduit, the `ResponseModel` seam and the
`bartcore.*` control-attribute mechanism are untouched.

## Context

- The decision, its rationale, and the exhaustive site list:
  [What goes](../design/retire-grouped-random-effects.md#what-goes) and
  [What stays, and why](../design/retire-grouped-random-effects.md#what-stays-and-why).
- What is being removed:
  [In-core grouped random effects](../design/grouped-random-effects.md#in-core-grouped-random-effects)
  (the landing) and [VERDICT (summary; detail below)](tau-slice-review.md#verdict-summary-detail-below)
  (the tau sampler review).
- The confounding defect that will not be fixed:
  [6. Recommendation: go/no-go](../design/forest-ranef-interweaving.md#6-recommendation-gono-go).
- ABI consequences and the structSize hole:
  [ABI and saved state](../design/retire-grouped-random-effects.md#abi-and-saved-state).
- Consumer migration:
  [Consumers](../design/retire-grouped-random-effects.md#consumers).
- The speed comparison, its firing and the ruling that followed:
  [The decision](../design/retire-grouped-random-effects.md#the-decision).

## Constraints

- **Gate cleared, with two release prerequisites attached.** The speed
  comparison has run. Its 5x tau-ESS-per-second threshold fired on two of the
  three seeds (paired ratios 1.5x / 8.9x / 303x), the decision was re-argued
  and stands, and Steps may start. Two RELEASE PREREQUISITES follow, both
  sister-repo work, filed from those repositories, neither blocking anything
  here:
  1. **stan4bart, tau mixing.** On the K = 20 gaussian design of the
     comparison (n = 2000, Friedman f, tau = 1, one chain, 1000 kept draws
     after 1000 warmup): lag-1 autocorrelation of tau below 0.8 and tau ESS
     at least 100 per 1000 draws. Absolute, so it neither cites `rbart_vi`
     nor needs a script that stops running after this lands; the comparison
     scripts and the three designs are committed to stan4bart's
     `benchmarks/` when the item is filed.
  2. **bartCause, a `group.by` route.** Its stan4bart branch refuses
     `group.by` today, so the one consumer has no working path until that
     changes -
     [Consumers](../design/retire-grouped-random-effects.md#consumers) enumerates
     the four pieces.
  dbarts ships first in the lockstep pair, so the ordering is the risk: once
  it is submitted neither prerequisite can be satisfied by delay.
- The conduit is frozen: `bartcore_setOffset`, `bartcore_setResponse`,
  `bartcore_setPredictor` and the R5 mutation surface change only by losing
  three `refuseGroupedScaleUpdate` call sites.
- `stateFormatVersion` does NOT increment: it is internal encoding, and this
  is pre-release.
- Neither `DBARTS_C_API_MAJOR` nor `DBARTS_C_API_MINOR` moves; only
  `DBARTS_C_API_HASH` re-bakes.
- Out of scope: reverting `refreshLatents`'s sigma parameter (a separate
  optional cleanup); any bartCause or stan4bart edit (VD's side, lockstep
  release); rewriting `docs/plans/release-candidate-review.md`'s ~36 grouped
  paragraph blocks or `docs/plans/review-2026-08-24/anchor-main.md`'s two
  anchor rows, both frozen history.

## Steps

File-disjoint except where noted; 1 must precede 2 only for the build to
stay green at each commit.

1. **Engine and bridge, and the ABI.** Delete `GroupedResponse`, the
   `TauPriorKind`/`logTauPrior`/`logTauPosterior`/`sliceSampleOnce`/
   `drawTauCauchyExactIG`/`drawGroupEffects` block and the four
   `ResponseModel` hooks from `src/bartcore/model.hpp`; the `SamplerOptions`
   group fields, the `Results` members, the decorator construction and the
   four state arms from `src/bartcore/chain.hpp`; `numGroups` from
   `sampler.hpp`, `facade.hpp` and `combiner.hpp`. Delete
   `src/R_interface_rbart.{cpp,hpp}` and its registration; delete
   `applyGroupAttribute`, `refuseGroupedScaleUpdate` (and its declaration),
   the three composition refusals, the two result channels, the two mutation
   guards, the `setData` refusal and the `"ranef"`/`"tau"` state slots from
   `src/R_interface_bartcore.cpp`; drop the two
   `refuseGroupedScaleUpdate` calls from `src/C_interface.cpp`. Remove `tau`
   and `groupEffects` from `dbarts_results`, from `DBARTS_RESULTS_FIELDS` in
   `src/C_interface.cpp` (the list macro driving the alignment asserts and the
   hash fold) and from `DBARTS_RESULTS_INIT`'s twelve positional
   initializers; re-index the offset asserts and take the size assert from 11
   pointer members to 9; amend the header's append-only field-discipline
   paragraph to say what a pre-1.0-0 removal does; re-bake
   `DBARTS_C_API_HASH` from the printed value. Two further C sites compile
   against the header and must move with it: `inst/tinytest/capi/consumer.c`
   (`capi_run_grouped` whole, plus the `results.tau`/`results.groupEffects`
   poison lines in the two structSize probes) and its driver block in
   `inst/tinytest/test-capi.R`, including the grouped-plus-variance-forest
   refusal. Determine and record what restoring a state written by a grouped
   build now does (ignored extra slots, or a shape-check error).
   `R CMD INSTALL . --preclean` (facade virtuals move).
2. **R surface, NAMESPACE, man, pkgdown.** Delete `R/rbart.R` and
   `R/sliceSample.R`; the nine `.rbart` methods from `R/generics.R`,
   `R/plot.R`, `R/bart.R` and `R/diagnostics.R`; the ten NAMESPACE lines; the
   `"rbart"` entry in `R/hooks.R`; the grouped clauses in `R/spec.R`,
   `R/partialDependence.R`, `R/data.R`, `R/utility.R`, `R/xbart.R` and
   `R/dbarts.R`. Also in `R/generics.R`, beyond the methods:
   `rbartUnusedArgs` and its four dispatch references, and the user-facing
   `'group.by' is the grouped (rbart_vi) fit's own predict argument` string
   (its pin in `inst/tinytest/test-generics-errors.R` goes in step 3).
   Delete `man/rbart.Rd`; strip grouped content from the eleven further Rd
   files; drop `- rbart` from `_pkgdown.yml`; drop the
   `Grouped random effects (rbart_vi)` bullet from the shipped `README.md`.
   Replace `vignettes/gibbs_sampler_mixture_model.Rmd`'s two `rbart_vi`
   paragraphs by inlining the deleted loop's core - a ~60-line
   random-intercept Gibbs over `$setOffset` - as the vignette's own worked
   example, rather than pointing at a file that no longer exists.
3. **Tests.** Delete the 17 wholly-grouped tinytest files and
   `inst/common/rbartGroupData.R`; strip the grouped block from the 36
   further tinytest files. Delete the five `testGrouped*` functions and their
   registrations from `tests/cpp/test_model.cpp`, the "grouped intercepts
   delegate" sub-block of `testActiveRows`, the `numGroups` assertion in
   `testConstantGaussian`, and the two `groupEffects`/`groupTau` comparisons
   in `tests/cpp/common.cpp`'s `ChainStateData` equality helper (a compile
   error if missed).
4. **Benchmarks and the baseline.** Delete two whole files:
   `benchmarks/R/grouped-mixing.R` and
   `benchmarks/R/forest-ranef-collapse-proto.R` (the isolation prototype
   behind the NO-GO doc step 5 supersedes). From `benchmarks/R/equivalence.R`
   the `grouped`/`grouped_aft` scenarios and `fitViaRbart`; from
   `benchmarks/R/composition-matrix.R` `probeRbart` and the `"groupedRanef"`
   entry; from `benchmarks/R/mutation-battery.R` the grouped poison test.
   From `benchmarks/R/sbc.R` the four configs, `sbcMakeGroupedFit`,
   `runSbcGrouped`, the `sbcAddGrouping` config extender, the `isGrouped`
   predicate and its four dispatch branches, the module header's usage line
   and four comment blocks - but NOT `sbcMatrixFunctionals`, a literal over
   the five non-grouped matrix arms, so no recorded verdict moves. Add a
   standalone `aft` scenario so AFT keeps equivalence coverage. Record the
   new baseline and update its four places in the recording commit
   (Verification below).
5. **Docs, INDEX and TODO.** Flip
   `docs/design/grouped-random-effects.md` to RETIRED and
   `docs/design/forest-ranef-interweaving.md` to SUPERSEDED, each with a
   one-line pointer to `docs/design/retire-grouped-random-effects.md`; the
   same for `docs/plans/group-by-exposure.md` and
   `docs/plans/tau-slice-review.md`. In `docs/design/feature-matrix.md`,
   remove the `grouped` row from all five tables, the flat-C grouped cell,
   footnotes [f8], [f13], [f14], [f31], [f32], [f37] and [f44], and the
   per-family `grouped` prose block; trim the grouped clauses from [f1],
   [f3], [f27], [f30], [f40] and [f50]. Trim the passing mentions in the 30
   other design docs. Update both INDEX files' status cells to match. In
   `TODO`, drop `group-by-exposure`, drop `negative-binomial`'s "grouped NB"
   clause, the multinomial entry's `rbart_vi` family-token clause and the
   ingestion-mode entry's "grouped frailty", and rewrite
   `sparse-extensions`'s rbart_vi-on-sparse half and
   `correlated-outcomes`'s "rbart_vi pattern" clause (the latter now names
   the `setOffset` conduit). Outside `docs/design/`: `docs/architecture.md`
   (the `R/` layer map and the `GroupedResponse` engine paragraph),
   `docs/plans/bartcore-review-tour.md` (three places, one of them the note
   that `grouped-mixing.R` disagrees with its own header),
   `tools/regenerate-snapshots.R`'s hardcoded
   `"test-reproducibility-rbart.R"`, the MANIFEST pin-site note naming
   `test-rbart-loop-callback.R`, and `.lintr`'s two `rbart_vi` rationale
   comments.
5b. **Cites, as its own pass - the gate in step 5 fails without it.**
   `tools/check-doc-freshness.R` resolves a cite's PATH before it honours a
   `retired:` marker, so `retired:` cannot rescue a cite into a deleted file.
   Two conversions, applied across every document under `docs/`, plus
   `README.md`, `man/*.Rd` and `vignettes/*.Rmd`:
   - a cite whose FILE is deleted becomes a history cite, `path:line`
     pinned to the pre-deletion tip;
   - a cite whose file survives but whose SYMBOL is deleted takes a
     `retired:` marker, with prose that says the construct is gone.
   Exposure measured at the tip: 45 breakable symbol cites in
   `docs/design/feature-matrix.md` (13 of them written against the `rbart.R` alias), 7 in
   `docs/design/forest-ranef-interweaving.md` - flipping that doc's status
   disarms nothing, the guard scans every `docs/**/*.md` and reads no status
   line - 7 in `docs/design/error-style.md`, and 12 across `survival.md`,
   `negative-binomial.md`, `ordinal.md`, `tree-mixing-proposals.md`,
   `threaded-predict.md`, `public-surface.md` and `active-rows-mask.md`. The
   design doc's own cites into deleted files are already in history form
   pinned at 916271f3; its surviving-file cites take `retired:` here.
   Finally drop `"rbart.R"` from `R_ALIAS_FILES` in
   `tools/check-doc-freshness.R` - it resolves unconditionally, with no
   existence test, so leaving it lets a cite into a deleted file resolve to
   nothing forever. Do this LAST, after no cite against the `rbart.R` alias remains.
   Plans under `docs/plans/**`, including `archive/` and
   `review-2026-08-24/`, use sha-pinned history cites already and need no
   change.
6. **NEWS.** 59 `\item` entries in the 1.0-0 section match on grouped text
   and are three classes; do not delete them wholesale, which would destroy
   the release record for surviving features.
   - DELETE the 6 wholly grouped items: `rbart_vi`'s `$fit` list shape,
     `predict.rbart`'s removed `value` argument, the in-core Gibbs item, the
     grouped-sampler response-mutability item, the seeded parallel seed
     restore, and the custom-prior-passed-by-name fix.
   - EDIT the 52 that also name `bart`, `bart2`, `xbart`, `dbarts`,
     `predict.bart` or `plot.bart`, dropping only the `rbart_vi` / grouped
     clause from each.
   - LEAVE the 1 false positive, the "grouped-GAMI decomposition" item.
   Then add the UPGRADING line from
   [The break versus 0.9-34](../design/retire-grouped-random-effects.md#the-break-versus-09-34),
   pointing users at stan4bart.

## Verification

RNG gate class: **shifting**
([RNG classes and their gates](README.md#rng-classes-and-their-gates)). No surviving fit
changes a draw, so the class is earned by the baseline losing scenarios, not
by a moved stream - which is exactly what the neutrality check below proves.

- `cd tests/cpp && make && ./test_bartcore` - all tests ok, five fewer.
- `R CMD INSTALL . --preclean` then `tinytest::test_package("dbarts")` - all
  pass; no RNG-locked snapshot needs regenerating, and any that does is a
  leak, not an expected shift.
- `Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-d4bca4ce.rds`
  (non-strict) BEFORE recording: expect "identical draws (same RNG stream)"
  on all 49 surviving scenarios, 2 skipped. A single non-bitwise scenario
  aborts the step. Then record `equivalence-<sha>.rds` (50 scenarios: the 49
  plus the new `aft`) and update its four places in the same commit - the pin
  in `.github/workflows/equivalence.yaml`, the MANIFEST row carrying that
  neutrality compare as its P17 oracle, feature-matrix [f39]'s count, and
  `equivBaseline` in `benchmarks/R/mutation-battery.R`. There is NO TODO
  ledger line to update: only one of the last three re-records touched TODO
  and no baseline pin lives there today. Demote `d4bca4ce` to historical.
  `equivalence.R compare <new> --strict-coverage` reproduces it 50/50.
- Exact gates, quick mode: unchanged verdicts (none is grouped, and
  `aft-exact.R` is the oracle for the new scenario).
- Sanitizers locally (ASan/UBSan build, tinytest) - the struct-layout change
  and the state-slot removal are the risk.
- `lintr::lint_package()` clean on the touched R files; air format check.
- `R CMD check --as-cran` from a clean `R CMD build` tarball - catches an
  orphaned `\link{rbart_vi}`, a stale alias and an undocumented removal.
- **`revdep-smoke.yaml`'s bartCause leg is EXPECTED RED from this landing**,
  and is an accepted, time-bounded red rather than a failure: bartCause's
  `use.ranef` path calls `dbarts::rbart_vi`, and its stan4bart branch refuses
  `group.by`, so nothing on that side works until release prerequisite 2
  lands. Record the red in the landing note with the prerequisite it waits
  on; it clears when bartCause's dbarts-1.0 branch does, and it must be green
  before the lockstep release. stan4bart's leg stays green (its `rbart_vi`
  call is inside `at_home()`), but `inst/tinytest/test-02-binary.R` is on its
  owed-edit list.
- `Rscript tools/check-doc-freshness.R` ends OK. It is the real gate on step
  5b, not step 5: a cite into a deleted file fails at path resolution, before
  `retired:` is consulted, so the marker cannot rescue one. Run it after 5b,
  not before.
- `Rscript tools/check-rc-codoc.R .` - `sampler.Rd`'s `run()` value docs lose
  `ranef`/`tau`.
- `Rscript -e 'tools::parse_Rd("inst/NEWS.Rd")'` parses.
- No `bench-sampler.R` run owed: nothing on the hot path changes for a fit
  that survives.

## Status

LANDED 2026-09-06 (1e5f80b2); two sister-repo prerequisites gate the release.

## Landing

1e5f80b2, the whole deletion in one commit: 171 files, 887 insertions and
9536 deletions. 25 files deleted outright, as the design enumerated. The
engine keeps the `ResponseModel` seam and every concrete family; the
per-sweep conduit loses only its three `refuseGroupedScaleUpdate` call
sites.

`DBARTS_C_API_HASH` re-bakes to `0x616ffcda8c947777`; neither version
constant moves. Restoring a state written by a grouped build is now
DETERMINED: it succeeds and ignores the two extra slots, because every
per-chain block is read by name and the only length check is on the chain
count.

Gates, all at this tip. tests/cpp 276 ok plain and again under
ASAN/UBSAN with no diagnostic. tinytest 7381/0 - no RNG-locked snapshot
needed regenerating, which is the leak check the Verification named.
Neutrality before recording: `equivalence.R compare` against
`equivalence-d4bca4ce.rds` reports "identical draws (same RNG stream)" for
all 49 survivors, with `grouped` and `grouped_aft` skipped (not produced)
and the new `aft` skipped (not in baseline). The 21 exact gates pass in
quick mode; `bcf-equivalence-3c81d6df` is 12/12 and
`multinomial-equivalence-4d9a3337` 11/11 bitwise. `air format --check`,
`lintr::lint_package()`, `pkgdown::check_pkgdown()`,
`tools/check-rc-codoc.R` and `tools/check-doc-freshness.R` are clean, and
`inst/NEWS.Rd` parses to 291 entries.

NEWS came out three items off the design's census: NINE `\item` entries in
the 1.0-0 section are wholly grouped and were deleted, not six. The three
the census missed - `group.by`/`group.by.test` looked up by name,
`rbart_vi`'s `k` formal defaulting to `NULL`, and the `ranef`-dimnames
crash fix - name a surviving entry point only in passing and have nothing
left once the grouped clause goes. 49 were edited, 1 (grouped-GAMI) left.

Four sites the plan did not name also carried the path and were fixed
here: `DESCRIPTION`'s feature list; two engine comments and one in
`tests/cpp/test_shape.cpp` citing `rbart_vi`'s callback loop as the
single-slab varcount example; and, contrary to step 5b's note that plans
need no change, 36 history cites written with a BARE BASENAME
(`test-grouped-swap.R`, `sliceSample.R`, `rbart.R`, ...) across 15 plan
files. The freshness guard resolves a history cite's path against the
CURRENT tracked inventory and only passes an unresolved token through
verbatim, so a basename that used to resolve stops resolving the moment
its file is deleted; each was spelled out to its full repo path at the
same sha.

`bartcore_runWithCallback` survives with its registration and no caller:
the R Gibbs loop was its only one, and no tinytest reaches it now. It is
left in place as the flat-API callback's R-side twin rather than deleted
in this pass.

revdep-smoke's bartCause leg is EXPECTED RED from here, waiting on release
prerequisite 2 (a `group.by` route that does not go through `rbart_vi`).
