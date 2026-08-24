# bcf-bartcause-relocation

The bcf-public-surface arc's S5 - `bcf()` and the `bartBCF` fit class - RELOCATED
to bartCause under multiforest-extension-surface's fork 4 (RESOLVED VD
2026-08-11), plus the dbarts-side records, guards and engine channel that
relocation needs. Written against dbarts @ 7948e40d (bartcore) and bartCause @
695c603 (dbarts-1.0); every anchor was read before it was cited, and every claim
marked MEASURED was produced by a read-only probe against a private library
carrying dbarts at that tip.

The design-session artifacts - the two-repo census, three blind critiques with
their adjudications, and the B0 suite measurement (2026-08-15/16) - sit in
`.claude/bcf-s5-relocation/`, which is GIT-IGNORED. This file is the only
tracked record and carries the load-bearing facts rather than pointing at them.

agent: D2 sonnet (three R-surface guards). D3 opus (the arc's one ENGINE slice).
  B0 sonnet (one snapshot literal). B0b sonnet (a three-token defect fix, two
  comments, two test oracles). B1 opus (the fit function, its class, the
  response fitter, the moderator exclusion).
sequence: **D2 -> D3 -> B0 -> B0b -> B1.** D0 has LANDED (42846863); its one
  remaining item is a post-B1 records touch.
rng: **D3 is the arc's ONLY engine-touching slice, and it carries the arc's only
  gated neutrality claim and its only baseline re-record.** Everything else is
  bitwise by construction: D0 touched no code; D2 is three R-surface refusal/
  warning paths plus one predicate widening that only ADMITS a previously
  impossible call; B0/B0b/B1 are bartCause-only, so no dbarts baseline can move.
  D3 changes what `storeSample` WRITES (a wider varcount slab), not what it
  draws, so the gaussian (`equivalence-8b047f8b`) and multinomial
  (`multinomial-equivalence-1027be5`) baselines must stay BITWISE and are the
  leak detectors; the bcf baseline needs a SHAPE-ONLY re-record because
  `benchmarks/R/bcf-equivalence.R`:100 records `varcount = result$varcount` and
  compares with `identical()` (:449), which a new dimension fails structurally.
  That re-record is gated on its own falsifier (FB16) proving the added
  dimension is the only difference.
window: inside the pre-release window; sister packages migrate in lockstep.
  Nothing here moves `DBARTS_C_API_HASH`, and that is MADE true rather than
  assumed: D3 gives `storeSample` the caller's declared forest count, so every
  caller-owned varcount buffer - the flat C API's and rbart_vi's - keeps today's
  exact single-slab bytes and `dbarts.h`:147's documented shape stays correct,
  while the R run bridge opts into the widened channel (D3 item 5). The only
  `dbarts.h` edit is one doc COMMENT, which the hash does not cover; no X-list
  entry is added, re-signed or removed, and no other slice adds a C entry.
budget (DENSE-EQUIVALENT lines: non-blank, non-comment, one statement per line;
  bartCause is hand-formatted and carries no air/lintr configuration):
  D2 ~22 R + ~30 test + ~10 Rd/NEWS.
  D3 ~8 engine (the two overrides plus storeSample's count authority) + ~5
  bridge + ~40 R (the packaging reshape, `n.forests`, the `fitSynopsis` arm) +
  ~25 Rd/NEWS/design + ~10 ledger (the re-record's four places) + ~110 test +
  a re-recorded bcf-equivalence baseline.
  B0 1 literal + 2 comment lines.
  B0b ~3 R tokens + ~8 R comment lines + ~6 test-oracle lines + ~10 new test
  + ~4 NEWS.
  B1 ~160 R (bcf.R core) + ~100 R (the response fitter, incl. the two builder
  returns) + ~120 R (five S3 methods) + ~25 R (bartc arm, refit/summary
  widening) + ~10 NAMESPACE + ~195 man + ~10 NEWS + ~310 test (~90 assertions).
  Total ~1215.

## Goal

`bcf()` and the `bartBCF` class ship in bartCause on its dbarts-1.0 branch, per
fork 4 (RESOLVED VD 2026-08-11, `docs/plans/multiforest-extension-surface.md`
:4325-4328). A user writes `bcf(y ~ x1 + x2, data, treatment = z)` - or
`bartc(y, z, x, method.rsp = "bcf")` - and gets a Bayesian causal forest whose
per-draw prognostic surface, treatment surface, glue, per-forest variable counts
and both counterfactual surfaces are reported, computed with NO test matrix, and
whose estimands, common support, plots and summaries are bartCause's existing
ones, unchanged. That last clause is measured rather than assumed: a design-session probe tested the
hypothesis that `sd.cf` collapses toward `sd.obs` under BCF (gutting
`commonSup.rule`) at n = 300, 2 chains, 200 draws against a matched BART fit,
and the two `(sd.cf/sd.obs)^2` distributions were indistinguishable.

dbarts' deliverable is records (D0, LANDED), one small R-surface guard slice
(D2) on the arc's construction seam, and one engine slice (D3) that makes S5's
per-forest varcount contract real. Revision 1's "dbarts gains nothing on any
branch of fork A" is withdrawn.

## Binding decisions inherited (do not reopen)

1. **Home: bartCause, dbarts-1.0.** Fork 4 RESOLVED (VD 2026-08-11).
2. **The moderator EXCLUSION for the propensity-score column is part of that
   cost, not a separate decision** (:4326-4328; bcf-public-surface.md:617-622).
3. **The contract survives verbatim** - per-draw mu, tau, glue, sigma,
   per-forest varcount and both counterfactual surfaces, forest-INDEXED, under
   the option-A element names (bcf-public-surface.md:493-501; the element list
   at :261-262). **As of VD's 2026-08-15 fork C resolution this now holds
   LITERALLY: the per-forest varcount channel is BUILT (D3), not dropped.**
4. **`predict` on new rows is COMPONENTS only; the blended test surface stays
   refused; per-forest saved-tree replay stays a DOOR** (:512-514, :660-662).
5. **`bcf(formula|x, y, treatment, ...)` vocabulary is bartCause-LOCAL.**
6. **No backwards-compatibility constraint.** Nothing is released; costs to
   consumers are enumerated, never design inputs.

## Adjudication: facts this plan is built on

Every item verified live at the anchors given; MEASURED items come from
read-only probes run against a private library at this tip.

- **The engine refuses a test surface on a multi-forest sampler, at creation.**
  `R/spec.R`:509 (`"test predictors" = !is.null(data@x.test)`) inside the
  `!is.null(data@bases)` branch at :413; `Chain::testFitsAreDefined` false for
  the amplitude coupling (`src/bartcore/combiner.hpp`:980; true for softmax at
  :1542); `refuseUndefinedTestFits` (`src/R_interface_bartcore.cpp`:2858; sites
  :4701, :4749, :4779, :5731; `src/C_interface.cpp`:728, :775). So bartCause's
  counterfactual machinery (`R/responseFit.R`:180-193, :245) is UNAVAILABLE, and
  so is the missing-response bookkeeping riding it. Under BCF the counterfactual
  is free instead.
- **The basis must be PRE-SUBSET on the formula path.** `validateForestBases`
  subsets only when handed an index (`R/data.R`:659-661); the x/y branches hand
  it one (dense :1041 under :993-995; sparse :975 under :937), the formula
  branch does not, deliberately (:890-893, comment :891-892), so a full-length
  basis fails the row check at :656-657. MEASURED:

      formula, full-length bases, no subset  -> data 20 x 3, OK
      formula, full-length bases, subset     -> "length of 'bases' must equal length of 'y'"
      formula, PRE-SUBSET bases, subset      -> data 10 x 3, OK
      x/y,     full-length bases, subset     -> data 10 x 3, OK
      dbarts(forests = ) formula + subset    -> "length of 'basis' must equal length of 'y'"

- **The supported construction composition is: bases on the DATA OBJECT,
  knob-only `forests`.** `resolveSamplerSpec` reads `data@bases` when the call
  supplies none (`R/spec.R`:206, branch :413); `dbarts()` passes `bases = NULL`
  from a pre-built object (`R/dbarts.R`:602-604). MEASURED end to end at a
  50-row integer subset of 60 rows, two chains: `bartcore.bcf` present,
  `forestFits` 50x2x12x2, `glue` 3x12x2, live varcounts forest 1 `0 5 3 8`,
  forest 2 `0 2 3 0` over `(z, x1, x2, ps)`. A logical subset works identically.
- **`dbarts()` and `dbartsSpec()` take OPPOSITE sources when both a data object
  and a declaration carry a basis.** MEASURED (2-column basis on the object,
  3-column declaration): `dbartsSpec()` -> ncol 3 (the DECLARATION's, installed
  at `R/spec.R`:405-411 from `bases = basis` at :714, expanded :682-695) and
  `dbarts()` -> ncol 2 (the DATA OBJECT's, its expansion dropped by
  `R/data.R`:742-751 with no warning). The `dbartsSpec()` direction is
  deliberate and pinned by four live assertions
  (`inst/tinytest/test-bcf-creation.R`:312-325 with intent at :310-311,
  :350-357, :751-758, :762-769 with semantics at :759-761). D2 makes the silent
  side loud and states the contract; `dbartsSpec()` is not touched.
- **The counterfactual identity, and the offset.** MEASURED on an
  offset-bearing gaussian BCF:

      max| response.scale*(a*mu + b_z*tau) + response.shift        - train | = 0.7607
      max| response.scale*(a*mu + b_z*tau) + response.shift + off  - train | = 8.882e-16

  so the FULL reconstruction needs the offset, while the difference form this
  plan ships needs neither offset nor shift because both cancel:

      mu.hat.cf = train + response.scale * (b_{1-z} - b_z) * tau_internal

  MEASURED to move the surface by up to 2.772 on that fixture. `b_0` multiplies
  the `(1-z)` indicator and `b_1` the `z` indicator, matching
  `model.matrix(~ factor(z) - 1)` order and the `(a, b_0, b_1)` stacking
  (`man/dbartsSampler-class.Rd`:366).
- **`$getCalibration()` returns one row per chain** (MEASURED 2 x 12;
  `response.scale`/`response.shift` identical across chains).
- **Chain dimensions drop at one chain** (MEASURED: `train` 60x20 vs 60x20x2,
  `forestFits` 60x2x20 vs 60x2x20x2, `glue` 3x20 vs 3x20x2, `sigma` a vector vs
  a matrix), and `run(0L, 0L)` returns NULL outright, so a zero-burn fit has no
  burn sigma.
- **`run()$sigma` exists under probit** (MEASURED, all 1s), so the binary arm's
  deletion is ACTIVE. `run()$tau` is a pre-existing channel unrelated to BCF's
  tau.
- **`crossvalidate` cannot work for BCF on any branch** (`R/model.R`:470;
  `R/spec.R`:486-491; `R/bayesOpt.R`:18, :20, :56; `R/xbart.R`:1-29 has no
  `forests`).
- **Missing responses cannot be served in v1** (`R/responseFit.R`:157-193,
  :262 rides the refused test surface).
- **The moderator exclusion reaches the engine on EVERY forest, including forest
  1** (`R/spec.R`:555-558, masks parked :559-564; MEASURED live varcounts
  `0 5 3 8` / `0 2 3 0`, batched `run()$varcount` sums `0 180 142 233`).
- **The two ingestion paths order `colnames(data@x)` differently** (MEASURED:
  `x1,x2,ps,z` vs `z,x1,x2,ps`) and nothing depends on it - `setdiff` preserves
  order, `resolveModerators`'s `match` (`R/model.R`:705-713) is name-keyed.
- **`blocks` is asymmetric under a column mask; `interactions` is not.**
  MEASURED: the top-level `blocks` (and `forest(1, blocks = )`, which
  `R/spec.R`:224-226 hoists onto it) is resolved at :341 with no
  `availableColumns` and demands a total partition including masked columns;
  forest 2's own `blocks` does not; `interactions(max.order = 1L)` is accepted.
- **The varcount forest axis EXISTS end to end and is exercised by the softmax
  coupling.** `numVariableCountForests()` returns 1 and `variableCountForest(j)`
  returns `reportedForest()` (`combiner.hpp`:651-653, doc :641-650), which is 0
  (:619), unoverridden by `BCFForestCombiner` (:713); the softmax coupling
  overrides both (:1536-1537). `storeSample` loops the writes over the count
  (`chain.hpp`:5005-5019); the bridge sizes the array off
  `shape.numVariableCountForests` (`R_interface_bartcore.cpp`:4219) and ALREADY
  allocates a forest dimension when it exceeds 1 (:4293-4300,
  `allocChannelArray(INTSXP, numPredictors, numVCForests)`), with the seam
  comment at :4216-4218 ("inserting a forest dimension between the predictors
  and the samples exactly as the fits seam inserts locations"). R has the
  matching reshape (`shapeMultinomialChannel`, `R/bart.R`:1363, used for
  varcount at :1440-1446) and the matching documentation (`man/bart.Rd`:445).
  This is what makes D3 six engine lines rather than a channel build.
- **A BCF is creatable through `bart2` today and reports nothing** (MEASURED: no
  `forestFits`, no `glue` in the fit; `packageBartResults`, `R/bart.R`:158-320).
- **`amplitude.prior.variance` is refused on the supported route** (MEASURED);
  `R/model.R`:875-882 tests `spec$basis` where its sibling at :883-891 honours
  `excused`, contradicting the doc at :855-856. D2 item 3.
- **bartCause's flat-sigma reshape has the wrong `byrow`** (MEASURED three ways,
  including a first-principles autocorrelation test), and **the suite MIRRORS
  the defect** (`tests/testthat/test-12-plot.R`:54-55, `test-05-generics.R`:46).
  B0b.
- **The `fit.rsp` element contract is four elements** - `$sigma`, `$first.sigma`,
  `$y`, `$fit` - and the census is COMPLETE (full scan; everything else is
  stan4bart-gated or predict-gated). With `first.sigma` absent `plot_sigma`
  errors at `n.chains > 1` (MEASURED).
- **The propensity-score name cannot be re-derived from `colnames(data@x)`**
  (MEASURED: with a `"psps"` confounder the ladder answers `psps` where the
  score is `ps`), and taking it from the builders is bart-path-neutral because
  `massign`'s unnamed loop (`R/multipleAssignment.R`:84-98) consumes
  positionally and never touches extras.
- **Every subset kind works on the bart path today** (MEASURED: integer,
  logical, negative, character), and a basis built from a bare vector carries
  POSITIONAL row names, so character subsets need normalization.
- **The R/C division puts the outer loop on R's side**
  (`docs/design/r-c-division.md`:84-95, :104-107).
- **`setControl` is not refused on a BCF sampler**, so the driver may mirror
  `bart2`'s two-phase burn (`refuseBCFMutation`'s five sites are
  `R/bartcore.R`:289, :310, :387 and `R/dbarts.R`:1067, :1522;
  `runWithBurnIn`'s test-surface calls are gated at `R/bart.R`:445).
- **The default user library is NOT HEAD**; every bartCause measurement runs
  against a private library built from the dbarts tip under test.

## Slice decomposition

Six slices. D0 has LANDED. D2 -> D3 -> B0 -> B0b -> B1.

### D0 (dbarts, records) - LANDED 42846863

Landed 2026-08-15, TODO-only (52 insertions, 26 deletions), firing no workflow.
It closed `dbartsdata-treatment-slot-debt` against 9c63e9d8 with the four
receipts and the 533-assertion measurement (now `TODO`:259-284, which also
corrects the "tinytest" error at :283-284 and records B0's pending literal at
:275-278); resolved the contradiction with the `bcf-public-surface` entry
(:279-282); repinned `calibrationMapName` to `R_interface_bartcore.cpp`
:4052-4055, literal :4054, call site :4144 (now `TODO`:336-337); completed
`bcf-naming-generalization`'s symbol list with `BCFSpecStorage` and the R-side
`isBCFSampler`/`refuseBCFMutation` pair (:146, :151-159); and split the dead
init-capture into its own entry, `dead-bcf-init-capture` (:285-294), whose text
states it "rides the next engine-adjacent slice rather than a records-only
commit" - **that slice is D3, which carries the deletion.**

Remaining, after B1 lands: update the `bcf-public-surface` entry
(`TODO`:174-192) to record S5 RELOCATED AND LANDED at bartCause dbarts-1.0
@ <hash>. A records touch, not a slice.

### D2 (dbarts, three R-surface guards on the construction seam)

R-only, no src/, no `dbarts.h`, RNG-neutral by construction.

1. **Refuse a discarded basis declaration - in `dbarts()` ONLY.** When
   `forestBasisDeclarations(forests)` yields at least one non-NULL basis AND
   `formula` is already a `dbartsData`, error by name: the declaration cannot
   reach a pre-built data object through `dbarts()`; use `dbartsSpec()`, or put
   the bases on the object with `dbartsData(bases = )`. **`dbartsSpec()` is NOT
   touched** - its first argument must be a `dbartsData` (`R/spec.R`:634-636),
   so the predicate is unconditionally true there, and it does not discard the
   declaration but installs it (:682-695, :714, :405-411), a route four live
   assertions pin. The refusal must also not fire on the supported composition
   (a pre-built object carrying `@bases` plus a knob-only `forests` list), which
   is what B1 uses.
1b. **State the contract in prose**, since the two surfaces genuinely differ:
   `dbartsSpec()` takes the DECLARATION as the source of truth and replaces
   whatever the data object carried (the pinned semantics at
   `test-bcf-creation.R`:759-761); `dbarts()` requires the bases to reach the
   data object it builds, and now refuses rather than dropping the declaration.
   One sentence each in `man/dbarts.Rd`, `man/dbartsSpec.Rd` and
   `man/forest.Rd`'s `basis` item. With item 1 in place no silent combination
   remains on either exported surface, so this closes the asymmetry without a
   door.
2. **Warn on the silently-ignored `bases` argument.** `R/data.R`:743-750's
   ignored-argument warning tests data/test/offset/offset.test only; add
   `bases`. (`dbartsSpec()` does not route through it.)
3. **Honour `hasBasis` in the amplitude-prior check.** `R/model.R`:875-882 tests
   `spec$basis` where its sibling at :883-891 tests `excused`. MOVE line 883
   above 875 - do not duplicate the expression. Forest 1's own refusal survives
   (its `bases[[1]]` is NULL, so `excused` is FALSE there), verified.

Consumer note: item 1 is the dbarts-side twin of stan4bart's
`bart-args-forests-guard` (stan4bart `TODO`:11-18), where a forwarded `forests`
builds a two-forest sampler and dies at the first `setOffset`; here it dies
never. Cross-reference both TODOs.

D2 does not block B1: it makes `b.prior.variance` settable (item 3) and one
wrong call loud (item 1); if it slips, B1 ships with `b.prior.variance` accepted
only at its default and a refusal naming D2.

Gates: the full standing dbarts battery - `R CMD INSTALL --preclean` into a
fresh private library; `tests/cpp` from clean (plain only; no src/ delta); full
tinytest with NO snapshot regenerated; the trio BITWISE
(`equivalence-8b047f8b`, `bcf-equivalence-8b047f8b`,
`multinomial-equivalence-1027be5`, pinned at
`.github/workflows/equivalence.yaml`:61, :87, :113); `air format --check .`;
`lintr::lint`; full local `R CMD check`; CI six-green.
rng: NEUTRAL - two refusal/warning paths and one predicate widening that only
ADMITS a previously-impossible call.

### D3 (dbarts, ENGINE: the per-draw per-forest varcount channel)

VD's fork C resolution (2026-08-15) reverses the drop-and-ticket
recommendation: the channel is the durable solution other multi-forest consumers
will want, so it is built now. This is the arc's ONLY engine-touching slice.

**Shape and naming - DECIDED: WIDEN the existing `varcount` channel's forest
axis; do NOT add a second channel.** The alternatives and the reasons:

- **(chosen) Widen `varcount`.** `run()$varcount` on a K-forest sampler becomes
  `n.predictors x n.forests x n.samples x n.chains` (forest-major within a
  sample, prognostic first) - the RAW channel shape a multinomial run already
  produces. (The packaged FIT shape is a different, draws-first array;
  `shapeMultinomialChannel` maps one to the other, and item 3 is where that
  happens. The two are easy to conflate and are distinguished throughout below.)
  Four independent facts make widening the coherent answer rather than a
  preference: the mechanism is documented as the varcount channel's own forest
  axis (`combiner.hpp`:641-650, "how many forests the per-sample split-usage
  channel records"); the bridge already allocates that dimension when the count
  exceeds 1 (`R_interface_bartcore.cpp`:4293-4300) and its seam comment says so
  (:4216-4218); the softmax coupling already ships the raw shape
  (`combiner.hpp`:1536-1537) with an R reshape helper
  (`shapeMultinomialChannel`, `R/bart.R`:1363, varcount at :1440-1446) and
  documented fit-side dimnames (`man/bart.Rd`:445); and
  `$getForestVariableCounts` (`R/dbarts.R`:1485) already reports per forest,
  1-based, so the per-draw channel and the live read finally address the same
  thing.
- (rejected) **An ADDITIVE per-forest channel** beside `varcount`, on the
  `forestFits`/`glue` precedent. That precedent does not transfer: `forestFits`
  and `glue` were NEW quantities with no existing channel to widen, whereas
  split counts already have a channel whose forest axis is built. Adding a
  second one leaves `varcount` reporting ONE arbitrary forest (the prognostic)
  for every multi-forest model forever - two ways to read the same quantity, the
  primary one silently partial - and diverges from multinomial.
- Cost of the chosen option, stated plainly: `varcount` is an existing
  DOCUMENTED channel, so the R5 run channel's BCF shape changes. **Who is affected is not
  obvious**: the R5 `run()` surface and the equivalence harness are not the
  whole census, and two earlier drafts of this plan stopped there. Items 3 and 5
  enumerate the two consumers that break, both measured.

**Contents.**

1. **Engine, two overrides** (~6 lines) in `BCFForestCombiner`
   (`combiner.hpp`:713; the member `numForests_` is initialized in its
   constructor at :730-731), mirroring the softmax pair at :1536-1537:
   `numVariableCountForests()` returns `numForests_`, `variableCountForest(j)`
   returns `j`. `storeSample` (`chain.hpp`:5005-5019) already loops the writes
   over a count and indexes through `variableCountForest(j)`; item 5 changes
   WHICH count it loops over. **One guard the override owes:** the
   constructor CLAMPS `numForests_(numForests < 2 ? 2 : numForests)`
   (`combiner.hpp`:731) while `Chain` passes `forests_.size()`
   (`chain.hpp`:816-818), and `storeSample` feeds the reported count to
   `forestVariableCounts(f, ...)`, which indexes `forests_[f]`
   (`chain.hpp`:1237-1238). A combiner constructed standalone with fewer than
   two forests would therefore report two and read out of range. Unreachable
   from R (`R/spec.R`:428-443 refuses K < 2) and there is no flat-C
   `createBCFSampler` (`combiner.hpp`:316), but `tests/cpp/test_sampler.cpp`:3724
   already constructs `BCFForestCombiner<...>(data, spec, 3)` standalone, i.e.
   the mismatch is constructible in the very build D3's new assertions live in.
   The new `tests/cpp` arm asserts the reported count against the chain's own
   forest count.
2. **Bridge: nothing structural** (~4 lines, comments plus item 5's one-line
   pin). Re-verified at the tip: `numVCForests` is read from the shape at
   `R_interface_bartcore.cpp`:4219, drives the allocation at :4293-4300, the
   buffer sizing at :4387 and the copy loop at :4431, and is set on the results
   at :4405. Comment edits: :4403-4404 and :4216-4218 say "1 for every additive
   model, K for multinomial" and become "K for multinomial and for a
   multi-forest amplitude model".
3. **R packaging: the break is the RESHAPE, not the naming - REDIAGNOSED.**
   `bart2(<a dbartsData carrying bases>)` runs today (MEASURED) and packages the
   channel through `nameVarcount` (`R/bart.R`:146-156), whose FIRST act is
   `convertSamplesFromDbartsToBart(raw, n.chains, combineChains)`
   (`R/bart.R`:147, body :8-38). That helper dispatches on
   `length(dim(samples)) == 2L` and otherwise does `aperm(samples, c(3L, 2L,
   1L))` (:20, uncombined) or `t(matrix(samples, x[1L], prod(x[-1L])))` (:30,
   combined). Neither is 4-D-aware. MEASURED on a p=3, K=2, S=4, C=2 raw array
   against the 3-D control:

       3-D control  combineChains = FALSE -> 2 x 4 x 3     TRUE -> 8 x 3
       4-D post-D3  combineChains = FALSE -> ERROR "'perm' is of wrong length 3 (!= 4)"
       4-D post-D3  combineChains = TRUE  -> 16 x 3, predictor names correct,
                                             row 2 = forest 2 of sample 1

   So the DEFAULT path (`bart2`'s own default is `combineChains = TRUE`) folds
   the forest axis into the draw axis forest-fastest and returns a matrix
   indistinguishable from a legitimate combined varcount except by row count,
   and the uncombined path errors from inside `aperm`. A "name only the
   predictor margin" alternative is DEAD: naming happens after the damage.

   **The fix, three parts.** (a) Route a multi-forest varcount through
   `shapeMultinomialChannel` (`R/bart.R`:1363), which is already the K-margin
   reshape and handles all three cases - MEASURED on the same raw array: 8 x 3 x
   2 combined, 2 x 4 x 3 x 2 uncombined, 4 x 3 x 2 at one chain, with the
   trailing margin named. The forest count comes from the sampler, not from the
   array's rank (a single-forest multi-chain raw channel is 3-D too):
   `numForests <- if (!is.null(fit$data@bases)) length(fit$data@bases) else 1L`,
   the same `data@bases` probe `isBCFSampler` uses (`R/bartcore.R`:24-26). The
   K-margin names are ENGINE vocabulary - `paste0("forest", seq_len(K))`, never
   causal names; bartCause's `bartBCF` maps them onto `mu`/`tau` in B1.2.
   (b) Record the forest count on the fit as `n.forests` (present only when
   greater than 1, so every existing fit is byte-identical), mirroring the
   multinomial packager's `K` (`R/bart.R`:1433). (c) Teach `fitSynopsis`
   (`R/generics.R`:1746-1757) that count. It derives "kept draws (per chain)"
   from `dim(x[["varcount"]])` whenever no sampler was kept, via
   `length(varcountDims) == 3L -> varcountDims[2L]` else
   `varcountDims[1L] %/% n.chains`. MEASURED today on a `bart2` BCF fit
   (`n.chains = 2L`, `n.samples = 5L`, no kept sampler): `varcount` 10 x 3 and
   "kept draws (per chain): 5" - correct; post-D3 without (c) the row count
   doubles and the synopsis prints 10. After (a) the packaged rank is 3 when
   combined or single-chain and 4 when uncombined, so RANK ALONE IS AMBIGUOUS
   (a single-forest uncombined varcount is also rank 3) - which is why (b)
   exists. The arm mirrors `print.bartMultinomial`'s arithmetic verbatim
   (`R/generics.R`:533): `if (length(d) == 4L) d[2L] else d[1L] %/%
   n.chains`, taken when `n.forests > 1`. `print.bartMultinomial` is why
   multinomial never exposed this: a BCF fit is class `"bart"` (MEASURED) and
   goes through `fitSynopsis`.

   **Rider, stated rather than smuggled:** routing through
   `shapeMultinomialChannel` means D3 ships a public `bart2` output shape
   (draws-first, forests on the trailing margin) for the very configuration the
   K-forest front door is DOORED on. D3 TAKES that decision, for the varcount
   channel only, in the shape multinomial already ships - the door's spelling
   review inherits it as a fact, not a constraint, since the door is about how a
   caller REQUESTS a K-forest fit and how `forestFits`/`glue` reach the fit.
   Scope is exact: D3 makes an EXISTING accidental path CORRECT; it does not
   advertise it. The fit still carries no `forestFits`/`glue`, and no Rd claims
   `bart2` as a multi-forest front door.

   (rejected alternative) **Refuse a multi-forest run result in
   `packageBartResults`** by name, and let the door decide everything. Cheaper
   (~6 R lines), and it pre-empts nothing. Rejected because it REMOVES a
   configuration that runs today for a user who reaches it deliberately through
   `dbartsData(bases = )`, and because a refusal is a worse answer than a
   correct array when the correct array is eleven lines and already has a
   helper. Recorded so it is not re-proposed as a shortcut.
4. **The live read is unchanged.** `$getForestVariableCounts(forest)`
   (`R/dbarts.R`:1485) stays the current-state, 1-based per-forest read; D3 adds
   the DRAW HISTORY it never was, and FB14 pins their relation.
5. **The flat C API keeps TODAY'S EXACT BYTES, and it takes one engine line to
   make that true rather than assumed.** The flat surface is NOT untouched by
   the override alone: `dbarts_sampler_create`
   (`src/C_interface.cpp`:369-376) routes through the same
   `bartcore_bridge::createHolder` the R bridge uses, so the flat API CREATES
   BCF samplers (`dbarts.h`:770-773 says so); `dbarts_sampler_run` sets only
   `engineResults.numReportedLocations` (`src/C_interface.cpp`:386-390) and
   never the forest count, leaving `bartcore::Results::numVariableCountForests`
   at its default 1 (`chain.hpp`:366); and that default is INERT today because
   `storeSample` derives the count from the COMBINER instead
   (`chain.hpp`:5010-5011) - whose own field comment (:363-365) says exactly
   that: "The run bridge sizes variableCounts by it and Sampler strides per
   chain by it; storeSample reads the count from the combiner directly." After
   item 1's override that split becomes a heap overflow: `Sampler::run` strides
   per chain by the CALLER's count (`sampler.hpp`:279, :287, :297-299) while
   `storeSample` writes K slabs per sample into a buffer the caller sized for
   one. The repo ships the consumer and the test: `inst/tinytest/capi/
   consumer.c`:186-190 allocates `p * numSamples * chains`, and
   `inst/tinytest/test-capi.R`:1236-1250 builds a two-forest spec through
   `dbartsSpec`, hands it to `capi_create` and runs it at `n.chains = 2L`,
   `n.samples = 3L` - a 2x overflow of an `R_alloc` block, on the per-push
   sanitizers workflow's own path.

   **The fix (one engine line plus one bridge line).** `storeSample` reads
   `results.numVariableCountForests` instead of asking the combiner. **The clamp
   to the combiner's count happens ONCE, upstream, where the stride is already
   computed**: `Sampler::run` writes the per-chain `Results` at
   `sampler.hpp`:287 (`r.numVariableCountForests = numVarCountForests;`), so
   clamping `numVarCountForests` to the combiner's count there makes the stride
   (`sampler.hpp`:279, :297-299) and the write agree BY CONSTRUCTION and
   `storeSample` needs no clamp of its own. Clamping only in `storeSample` would
   leave a caller who declares MORE than the combiner's count striding per chain
   by the declared count while the writes pack at the clamped one - misaligned
   contents, never an overflow, and nothing in-repo does it, but the one-line
   spelling removes the case. Then: the flat C API, which leaves the field at 1,
   loops `j = 0` only and writes `variableCountForest(0)` = 0 = **the
   PROGNOSTIC forest** - byte-identical to today, and identical before and after
   item 1's override (the base returns `reportedForest()` = 0,
   `combiner.hpp`:619; the override returns `j` = 0) - so `dbarts.h`:147's
   documented `numPredictors x numSamples x numChains` stays TRUE, no X-list
   entry moves and `dbarts_apiHash` does not move. The R run bridge already sets
   the field from the shape (`R_interface_bartcore.cpp`:4405), so it opts in.
   `bartcore_runWithCallback` (:4512) sets it from the shape too and is PINNED
   TO 1 here, which is also MAJOR 3's enforcement (see the gates).

   **Silent prognostic-only is acceptable ONLY documented, so it is documented
   in three places**: `dbarts.h`'s `varcount` doc gains one sentence - the
   channel reports the sampler's REPORTED forest (the prognostic forest of a
   multi-forest model), the flat API declares no forest count, and a caller
   wanting per-forest counts drives the sampler from R; `docs/plans/
   c-api-growth.md` records the field-free contract as a reserved widening (a
   future `dbarts_results` field is a size-guarded append, not a signature
   change); and `test-capi.R` gains an assertion on `length(rBCF$varcount)`
   pinning the contract from the consumer's side. NOTE the doc sentence is a
   COMMENT in `dbarts.h`, not an X-list entry: `DBARTS_C_API_HASH` hashes the
   declaration list, so a comment cannot move it - verified by the same
   mechanism S3's re-bake documented.
6. **The bcf-equivalence baseline is RE-RECORDED, shape-only - and a re-record
   touches FOUR PLACES, all in the SAME COMMIT as the new `.rds`.**
   `benchmarks/R/bcf-equivalence.R`:100 records `varcount = result$varcount` and
   the harness compares with `identical()` (:449), so a new dimension fails
   structurally even though no draw moved. The re-record lands in this slice,
   gated on FB16, and it is the arc's ONLY re-record. The standing obligation is
   `docs/plans/multiforest-extension-surface.md`:3378-3391 ("**A re-record
   touches FOUR places, not one, and all of them land IN THE SAME COMMIT as the
   new `.rds` files** - a workflow pointing at a deleted baseline is a red gate
   that looks like a regression, and a ledger naming a deleted file lies about
   which baseline is current"), whose four legs, re-derived live at this tip:

   (1) **`.github/workflows/equivalence.yaml`:87**, which today names
   `benchmarks/baselines/bcf-equivalence-8b047f8b.rds`; :61 and :113 do NOT move
   (those baselines are unchanged).
   (2) **`benchmarks/baselines/MANIFEST`:42**, today
   `bcf-equivalence-8b047f8b.rds  current  8b047f8b`, which flips to
   `historical` with a new `current` row above it in the existing narrative
   format - scenario count, the neutrality PARTITION (all 12 scenarios' draw
   channels bitwise; `varcount` gains a forest axis whose slab 1 is identical to
   the old array), and the superseding hash. :16 and :48 stay as they are (the
   sibling `current` rows for the two unmoved baselines).
   (3) **The arc's TODO ledger entry**, in the `LANDED <hash> + baselines
   <hash>` idiom the multiforest-predictor-mutation entry uses (`TODO`:535,
   :538, :540). The obligation's own citation (`TODO:258`, `:387`, `:389`) is
   stale at this tip - grep finds no baseline hash in `TODO` at those lines -
   so the live target is this arc's entry: `bcf-public-surface` (`TODO`:174-192)
   today, or the bcf-s5-relocation entry the spec commit creates, whichever is
   current when D3 lands.
   (4) **`docs/design/feature-matrix.md`:719-721**, the `[f39] Current
   baselines` line, which names all three baselines and cites
   `MANIFEST:16, 42, 48`.

   Why this is not optional: the equivalence workflow is `schedule` +
   `workflow_dispatch` only (`.github/workflows/equivalence.yaml`:20-23), so a
   missed yaml bump does not fail the landing commit - it fails the next cron,
   detached from its cause, and reads exactly like the regression FB16 exists to
   rule out. The harness header also needs its own edit: :19-22 says
   `result$varcount` "for a BCF sampler is the PROGNOSTIC forest alone, so
   varcount.tau is the only guard on the treatment forest's counts", which the
   widening falsifies; keep the `varcount.tau` channel (:101) as the live-read
   cross-check against the widened channel.
7. **The dead init-capture rides here.** `dead-bcf-init-capture`
   (`TODO`:285-294) says it "rides the next engine-adjacent slice"; delete
   `treatment = std::vector<double>{}` at `src/R_interface_bartcore.cpp`:2954
   and close the entry at landing. One line, no behavior.
8. **Docs and stale prose.** `man/dbartsSampler-class.Rd`:366's `varcount`
   sentence gains the forest axis, in the language `man/bart.Rd`:445 already
   uses for multinomial (distinguishing the RAW run shape from the packaged fit
   shape); `inst/NEWS.Rd` gains one bullet; `docs/design/bcf.md` and
   `docs/design/multiplier-combiner.md` record that the amplitude coupling now
   reports per-forest split counts; `dbarts.h`'s `varcount` doc gains item 5's
   sentence. **Four in-tree comments are FALSIFIED by this slice and must be
   corrected in the same commit**, three of them in the file the fix lands in:
   `chain.hpp`:316-318 (the `variableCounts` field: "the forest dimension is 1
   for every additive model [...] numCategories for multinomial");
   `chain.hpp`:362-365 ("The run bridge sizes variableCounts by it and Sampler
   strides per chain by it; storeSample reads the count from the combiner
   directly" - the sentence item 5 makes false, and the one this plan quotes as
   its own evidence); `chain.hpp`:5006-5009 ("count 1 (single forest and BCF) is
   the exact current byte layout; a multi-forest combiner (multinomial) records
   each category forest's splits"), falsified twice over and sitting inside the
   block item 5 edits; and `tests/cpp/test_shape.cpp`:293-294 ("Multinomial: the
   only path that widens the reported-location and variable-count-forest
   channels past 1"), in the file the new assertions live in.

**Falsifiers.**

- **FB14 (the per-draw / live-state oracle, both halves).** On a two-forest
  sampler with `n.thin = 1L` and no sweep after the last `storeSample`, the last
  kept draw's slab equals the live read for every chain:
  `run()$varcount[, f, n.samples, chain] == sampler$getForestVariableCounts(f)[, chain]`
  for f = 1, 2 (1-based at the R5 surface; the internal
  `bartcoreForestVariableCounts` is 0-based - the existing single-forest version
  of this oracle is `inst/tinytest/test-multi-forest-seam.R`:365-370). NEGATIVE
  HALF: run one more sweep and the equality must break, so the test proves the
  channel is a draw history rather than a repeated live read.
- **FB15 (mask structural zeros).** On a BCF whose treatment forest is
  restricted with `forest(vars = )`, every masked column's row in the treatment
  forest's slab is zero in EVERY draw, and at least one unmasked column is
  non-zero in some draw. NEGATIVE HALF: without the mask the same column is
  non-zero at the same seed.
- **FB16 (the re-record's own gate).** The widened channel's forest-1 slab is
  `identical()` to the pre-D3 recorded array, element for element, on every
  bcf-equivalence scenario. **Every scenario is SINGLE-CHAIN** (`makeControl()`
  fixes `n.chains = 1L` at `benchmarks/R/bcf-equivalence.R`:85 and all twelve
  scenarios call it), so the recorded array is a 2-D `numPredictors x n.samples`
  matrix and the post-D3 channel is 3-D `p x K x n.samples`: the oracle is
  `identical(new[, 1L, ], old)`, with no chain margin to slice. This is what
  separates "a dimension was added" from "a draw moved"; the re-record does not
  land until it passes.
- **FB17 (softmax and single-forest regression guards).** A multinomial fit's
  `varcount` is unchanged in shape, dimnames and values (bitwise, via
  `multinomial-equivalence-1027be5`, whose harness records both the live
  per-forest reads and `runVarcount = result$varcount`,
  `benchmarks/R/multinomial-equivalence.R`:91, :95 - so it is a real detector);
  and an ordinary sampler's `varcount` keeps its current rank and is bitwise
  (`equivalence-8b047f8b`, whose harness has no multi-forest scenario).
- **FB18 (D3, the `bart2` packaging path - NEW, and nothing else covers it).**
  `bart2(<a dbartsData carrying bases>)` at `combineChains = TRUE` and `FALSE`,
  at one and two chains: the packaged `fit$varcount` has the shape
  `shapeMultinomialChannel` defines (draws x p x K, with a leading chain margin
  when uncombined and multi-chain), carries predictor names on the p margin and
  `forest1..forestK` on the K margin, and its forest-1 slab equals the same
  quantity extracted from the `run()` channel it came from. Plus the synopsis:
  `print()` on that fit reports the TRUE kept-draw count (MEASURED today: 5 at
  `n.samples = 5L, n.chains = 2L`), and `fit$n.forests` is 2. NEGATIVE HALF:
  drop the `n.forests` arm from `fitSynopsis` and the printed count must double,
  so the test proves the arm rather than the arithmetic.
- **FB19 (D3, the flat-C legacy contract - NEW, and it is the BLOCKER-1
  detector).** Three legs. (i) `tests/cpp`: with the caller's declared count at
  1 against a K-forest combiner, exactly `numPredictors` are written per sample
  and the bytes are the prognostic forest's - the buffer is a
  `std::vector<std::uint32_t>` sized `p * numSamples`, the idiom
  `tests/cpp/rshim.cpp`:107 already uses. (ii) `tests/cpp`: the reported count
  never exceeds the chain's own forest count (item 1's clamp guard). (iii)
  `test-capi.R`:1236-1250's existing BCF leg gains a VALUE assertion - the
  returned counts equal the prognostic forest's live
  `$getForestVariableCounts(1L)` at the same seed - beside
  `expect_equal(length(rBCF$varcount), p * n.samples * n.chains)`; a length
  assertion alone cannot move, since the consumer computes the length itself
  from p, samples and chains.
  **NEGATIVE HALF, and it must be leg (i), not leg (iii):** revert `storeSample`
  to consult the combiner and leg (i) must abort under ASAN. Leg (iii) is NOT a
  reliable detector at the shipped fixture size - `consumer.c`:189-190 allocates
  through `R_alloc(p * numSamples * chains, sizeof(uint32_t))`, which at p = 3,
  numSamples = 3, chains = 2 is 72 bytes served from R's own small-vector pages
  rather than an individually malloc'd block, so it carries no ASAN redzone and
  a 2x overrun lands inside a page R already mapped. Leg (i)'s `std::vector` is
  malloc'd, redzoned and container-annotated, so the same overrun aborts
  deterministically. Leg (iii) also has a skip mouth: `test-capi.R`:36-38 calls
  `exit_file("could not compile the C API consumer")` when `R CMD SHLIB` fails,
  so "runs clean" can mean "did not run" - **the capi legs must be confirmed RUN,
  not skipped, for the gate to count.**

**Test budget: ~110 dense-equivalent lines / ~42 assertions.** ~16 tinytest
(FB14, FB15) in a new or extended `inst/tinytest/test-bcf-reporting.R` arm; ~10
`tests/cpp` (the shape predicate asserted on the predicate rather than a forest
count, the S4 pattern; plus FB19 legs (i) and (ii)); ~6 regression (FB17); ~8
tinytest for FB18 (the first tinytest to call `bart2()` on a bases-carrying data
object - grep confirms none does today); ~2 in `test-capi.R` (FB19 leg (iii)).
An earlier ~28 was understated by the two consumer-census falsifiers and by
item 6's ledger work.

**STOP CONDITIONS, read and write sites both** (a read-site-only condition
cannot see the flat-C failure mode item 5 closes):

- **No tinytest assertion is expected to move.** Every BCF varcount assertion in
  the suite is self-relative (`test-bcf-creation.R`:68 compares two BCF
  samplers; `test-multi-forest-seam.R`:207, :220 compares two BCF arms;
  `test-bcf-reporting.R`:135 asserts the run NAME list) and every
  positionally-indexed varcount assertion is on a single-forest or multinomial
  fit (`test-convergence-diagnostics.R`:34, `test-multithreaded.R`:50-51, the
  three `test-reproducibility-*` files, `test-multi-forest-seam.R`:370). The
  `tests/cpp` side is general too (`test_shape.cpp`:86-88 sizes its buffer by
  `numVariableCountForests`, :97 declares it). If one moves anyway, stop.
- **No caller-owned varcount buffer may be reached by a multi-forest sampler.**
  There are exactly two: `inst/tinytest/capi/consumer.c`:186-190 (driven by
  `test-capi.R`:1250) and `R/rbart.R`:599, :618. Both must stay clean with their
  sizing UNCHANGED. The DETECTOR for a botched fix is FB19 leg (i) in
  `tests/cpp`, whose `std::vector` buffer is redzoned; the capi leg is a
  contract pin, not an overflow detector, because `R_alloc` serves its 72 bytes
  from R's small-vector pages where ASAN cannot see a 2x overrun. A `tests/cpp`
  ASAN abort in the varcount arm is THIS slice, never an unrelated flake.

**Gates: the full standing engine battery**, plus three additions:
`R CMD INSTALL --preclean` into a fresh private library; `tests/cpp` from clean,
plain AND ASAN/UBSAN; full tinytest with NO snapshot regenerated; the gaussian
and multinomial baselines BITWISE (the leak detectors) and the bcf baseline
re-recorded only after FB16 passes; `air format --check .`; `lintr::lint`; full
local `R CMD check`; CI six-green.

- **Addition 1: the caller-owned-buffer audit, aimed at the invariant that can
  actually break.** The obvious proposition - that no grouped sampler reaches a
  count above 1, by `R/spec.R`:505 - is VACUOUS at the buffer: `bartcore.groups`
  is set at exactly one place, `R/rbart.R`:367, inside the IN-CORE fast path
  (`builtinTauPrior && is.null(callback)`, :363), which runs through
  `sampler$run` and never touches `bartcore_runWithCallback`; the path that DOES
  own a buffer is the R-loop path (`R/rbart.R`:736-739, buffers at :599 and
  :618), where no `bartcore.groups` attribute is ever set, so `R/spec.R`:505
  cannot fire there. MEASURED: a two-forest `dbartsData` handed to `rbart_vi`
  dies on the R-loop path at `setOffset(updateScale = TRUE) does not support a
  BCF sampler` (`refuseBCFMutation` in `bartcoreSamplerSetOffset`,
  `R/bartcore.R`:308-316, fired from the pre-run rescale at `R/rbart.R`:924-926
  and again at :961-964), and on the in-core path at the grouped refusal. So the
  audit becomes: (i) record BOTH guards by name at the buffer site, the real one
  being `setOffset`'s BCF refusal rather than `R/spec.R`:505; (ii) make the site
  self-enforcing rather than dependent on a guarantee three files away -
  `bartcore_runWithCallback` pins `results.numVariableCountForests = 1`
  (item 5), so the buffer's `p * n.samples` layout is true at the site whatever
  a future slice does upstream. The live door this protects is real:
  `setOffset(updateScale = TRUE)` on a BCF sampler is a recorded door, and the
  next slice that opens it must not inherit a false comment.
- **Addition 2: the re-record's four places, updated IN-COMMIT** - item 6's
  yaml line, MANIFEST row, TODO ledger entry and feature-matrix line, verified
  present in the same commit as the new `.rds` before the slice is called done.
- **Addition 3: the neutrality claim is GATED, not asserted** - D3 changes what
  `storeSample` writes, never what it draws (the varcount block calls only
  `forestVariableCounts` -> `countVariableUses`, a pure tree walk; the bridge's
  extra allocation is outside the `GetRNGstate`/`PutRNGstate` bracket at
  `R_interface_bartcore.cpp`:4421-4424), so the two untouched baselines must be
  bitwise; a divergence on either is a LEAK and aborts the slice.

### B0 (bartCause, restore green against dbarts HEAD) - MEASURED

Measured twice against a private library carrying dbarts @ 231744a0: run 1
(default env) `FAIL 1 | WARN 9 | SKIP 10 | PASS 497`, whose single failure was a
STALE stan4bart binary in the default user library - environment staleness,
resolved by rebuilding into the private lib with zero code changes, no action
owed; run 2 (`NOT_CRAN=true`, stan4bart rebuilt) `FAIL 1 | WARN 9 | SKIP 0 |
PASS 533`.

**Content: exactly one numeric literal** - `tests/testthat/test-06-regression.R`
:52, `0.50402744431802` -> `0.44963899452561451`, inside the
`tmle_version >= "2.1"` branch (:48) behind `skip_on_cran()` (:38) and
`skip_if_not_installed("tmle")`, at `n.chains = 1L` (:45). Lines :10, :20, :30
already carry refreshed values and PASS.

**Gate mechanics.** Plain `Rscript tests/testthat.R` and `R CMD check` do NOT
set `NOT_CRAN`, so the block that consumes the literal is SKIPPED under both;
the gate must set `NOT_CRAN=true` explicitly. Regenerate IN SUITE (whole file),
never from an isolated replication, because the draws depend on the file's
execution history (`set.seed(22)` plus `bartc()`'s own `sample.int` seeding at
`R/bartc.R`:94-95, :167-168).

Gates: full testthat under the private library with `NOT_CRAN=true`;
`R CMD check --as-cran` on a tarball staged from a clean copy outside the tree.
bartCause has no air/lintr configuration and only a manual R-hub workflow.
STOP CONDITION: if the regenerated value differs from `0.44963899452561451`,
stop - the environment moved and the measurement must be re-run first.

Note the sequencing consequence of D3: B0 runs against a dbarts that now
reports a 4-D BCF varcount. bartCause's bart path never sees it (no BCF fit
exists yet), so B0's literal is unaffected; B1's fixtures are the first
consumers.

### B0b (bartCause, the per-chain sigma defect)

Not two tokens: the suite MIRRORS the defect, so the fix moves two existing
oracles and changes one shipped output.

1. Three `byrow = TRUE` tokens dropped: `R/generics.R`:332, `R/plot.R`:19, :20
   (MEASURED: the only three in the package's `R/`).
2. Two comments asserting the false layout rewritten: `R/plot.R`:14-17,
   `R/generics.R`:328-330.
3. Two existing test oracles re-derived: `tests/testthat/test-12-plot.R`:54-55
   (plus its false comment at :31-35) and `test-05-generics.R`:46, the latter
   becoming `expect_equal(sigma, as.vector(t(matrix(fit$fit.rsp$sigma, nrow =
   n.chains))))`.
4. **A user-visible output change with its own NEWS line:** `extract(fit,
   "sigma")` at its default `combineChains = TRUE` currently returns the raw
   flat vector (sample-major) because `byrow = TRUE` makes `combineChains`'s
   `as.vector(t(.))` invert the reshape; after the fix it returns chain-major
   order, matching every other extracted quantity (`R/generics.R`:9). NEWS:
   "`extract(type = \"sigma\")` now returns posterior draws in chain-major
   order, matching every other extracted quantity, and per-chain sigma is now
   correctly assigned to its own chain; previously multi-chain fits interleaved
   chains, which also mis-paired sigma with the posterior predictive draws."
5. A new pin (FB12).

STOP CONDITION: no test outside `test-12-plot.R` and `test-05-generics.R` may
move. B0/B0b non-interference: `test-06-regression.R`:45 runs at
`n.chains = 1L` and both reshape sites are gated on `n.chains > 1L`. B1
non-interference: B1.2a stores matrices, so neither reshape fires for a bcf fit
either way.

Gates: as B0.

### B1 (bartCause, S5 proper)

One commit: `bcf()`, `bartBCF`, five S3 methods, the `method.rsp = "bcf"` arm,
the moderator exclusion, the two builder returns, the `refit`/`summary`
widenings, `man/bcf.Rd`, NEWS, `tests/testthat/test-14-bcf.R`.

#### B1.1 The engine-facing core

One internal fitter, `fitBCF(...)`, is the only place that touches dbarts. Four
calls, not a loop:

```
sampler <- dbarts::dbarts(responseData, control = ctl, forests = list(mu, tau), ...)
sampler$sampleTreesFromPrior(updateState = FALSE)
burn    <- sampler$run(0L, ctl@n.burn,    updateState = FALSE)
samples <- sampler$run(0L, ctl@n.samples, updateState = FALSE)
```

mirroring `bart2`'s standard path (`R/bart.R`:1019, :1022; the two-phase form
:435-475, including the draw-neutral `keepTrainingFits = FALSE`/
`verbose = FALSE` during burn, available on a BCF). `n.chains`, `n.threads`,
`n.burn`, `n.samples`, `keepTrees`, `rngSeed` and `verbose` are `dbartsControl`
slots (`R/dbarts.R`:214-231). The per-draw channels arrive BATCHED from the
second `run()` - `forestFits`, `glue`, and (after D3) the per-forest `varcount`;
the burn run's `$sigma` becomes `first.sigma` (`R/bart.R`:452-455, packaged
:245-251, :268). At `n.burn == 0` there is no burn run at all (`run(0L, 0L)`
returns NULL, MEASURED) - see B1.2a.

Forest declaration, defaults from `forestParams` (`R/model.R`:964-984), with NO
`basis` field on either forest (the bases ride the data object, B1.1a):

- forest 1, prognostic: `forest(vars = muVars, sd = sd.control,
  update.amplitude = update.a, interactions = mu.interactions)`. No basis, so
  its amplitude is the scalar `a` under a half-Cauchy mixture whose median is
  `sd` (`R/model.R`:980). Tree count and structure prior are the fit's own
  (`R/spec.R`:217-226). `mu.blocks` is REFUSED - B1.3.
- forest 2, treatment: `forest(vars = tauVars, n.trees = n.trees.treatment,
  base = treatment.base, power = treatment.power, sd = sd.moderate,
  amplitude.prior.variance = b.prior.variance, update.amplitude = update.b,
  interactions = tau.interactions, blocks = tau.blocks)`. With a basis reaching
  it from `data@bases`, `sd` multiplies the node scale through 0.674
  (`R/model.R`:977-978) and `amplitude.prior.variance` is the fixed-variance
  channel (:979); it needs D2 item 3.

**Response family:** `gaussian`, `probit`, `logistic` (`R/spec.R`:449-471). For
a binary fit the combination is on the LATENT scale and the link is applied LAST
to both surfaces, as `extract.bart` does (`R/generics.R`:303-305); the `sigma`
element is then DELETED, since `run()` supplies 1s (MEASURED) and
`responseIsBinary` (`R/bartc.R`:211-217) keys on its absence.

**The response transform:** read `cal <- sampler$getCalibration(1L)` before the
sampler goes out of scope, take `cal[1L, "response.scale"]`, and assert the rows
agree.

#### B1.1a Construction under `subset`

- **(i) CHOSEN: the fitter normalizes the subset and pre-subsets the basis**,
  then uses the bases-on-the-data-object composition. bartCause already resolves
  the subset expression one line before it builds the data object
  (`R/responseFit.R`:161); the bcf fitter does the same, NORMALIZES to positive
  integer positions, builds the basis from the FULL-length treatment vector,
  subsets it, binds it as a symbol and passes `dbartsDataCall$bases <- <symbol>`
  - the idiom `bartc()` uses for the propensity score (`R/bartc.R`:148-153).
- (ii) Refuse `subset` for bcf. REJECTED: breaks `man/bartc.Rd`:397's canonical
  workflow (a LOGICAL subset) and `tests/testthat/test-04-bartc.R`:252-256 (an
  INTEGER one), for a mechanism measured to work.
- (iii) Make the formula branch subset the bases. REJECTED for this arc: it
  makes a shipped argument's contract ambiguous at equal row counts, and the
  asymmetry is deliberate (`R/data.R`:630-632, :891-892). A door.

**The normalization is TOTAL, because all four kinds reach the bart path today
(MEASURED, `n.obs = 20` each):**

```
s <- if (is.logical(s))        which(s)
     else if (is.character(s)) match(s, rownames(<the frame the subset addresses>))
     else if (any(s < 0))      seq_len(n)[s]
     else                      s
```

logical -> `which()`; an NA-bearing logical subset is refused TODAY on both
paths ("response contains missing values"), so it is stated, not guarded.
character -> `match()` against the frame's row names, because a basis built from
a bare vector carries POSITIONAL names (MEASURED: `1,2,3` vs `r1,r2,r3`);
unmatched names refuse by name. negative -> `seq_len(n)[s]` (MEASURED aligned
either way; written out to make the normalization total). positive integer ->
identity. MEASURED through `dbartsData` with this normalization: all four kinds
give 20 rows of 30.

Two guards:

- Build the basis from the FULL-length z, then subset. Building from the subset
  z ERRORS when an arm is lost (MEASURED: "contrasts can be applied only to
  factors with 2 or more levels"), a confusing failure from inside a fitter.
- **Refuse when either arm is empty after subsetting** - the silent hazard: an
  all-zero basis column is ACCEPTED and fits with a dead amplitude (MEASURED).

#### B1.2 The `bartBCF` structure

| element | shape | content |
|---|---|---|
| `forests` | named list, length K | forest-indexed per-draw INTERNAL-scale fits; names `c("mu", "tau")` at K = 2 |
| `forests$mu`, `forests$tau` | [n.chains, n.samples, n.obs] | one slab of `run()$forestFits` each |
| `glue` | [n.chains, n.samples, 3] | dimnames `c("a", "b.0", "b.1")` on the last axis |
| `mu.hat.obs` | [n.chains, n.samples, n.obs] | `run()$train`; latent scale for a binary fit until the link is applied |
| `mu.hat.cf` | [n.chains, n.samples, n.obs] | `train + response.scale * (b_{1-z} - b_z) * tau`, same scale as `mu.hat.obs`; the LINK is applied LAST to both |
| `sigma` | [n.chains, n.samples] matrix | ABSENT on a binary fit (an active deletion) |
| `first.sigma` | [n.chains, n.burn] matrix | `[n.chains, 0]` when `n.burn == 0`; ABSENT on a binary fit |
| `y` | length n.obs | the training response |
| `varcount` | named list, length K, each [n.chains, n.samples, n.pred] | **the PER-FOREST per-draw split counts from D3's widened channel**, forest-indexed and named exactly as `forests` is, with predictor names on the last margin. The raw channel is `p x K x n.samples x n.chains` and DROPS the chain margin at one chain (`p x K x n.samples`), so the fitter's reshape has two cases, exactly as it does for `forestFits` |
| `trt`, `name.trt`, `name.p.score`, `family`, `call`, `n.chains` | | `name.p.score` is the builders' own resolved name (B1.4), NULL when no score was supplied |
| `fit` | `dbartsSampler` | only when `keepSampler = TRUE` |

`class(result) <- "bartBCF"`. `varcount` mirrors `forests`' container so the two
forest-indexed quantities are addressed identically and a K > 2 model needs no
reshape. **The structural-zero rider stands and `man/bcf.Rd` must carry it:**
under the masks this design installs, the treatment forest's slab is
structurally zero on the treatment column and on the propensity column, and the
prognostic forest's is structurally zero on the treatment column - a zero there
is the mask, not a finding about the data.

#### B1.2a The `fit.rsp` element contract

| element | read by | required shape |
|---|---|---|
| `$sigma` | `R/plot.R`:13, `R/generics.R`:331, `R/bartc.R`:215 | `[n.chains, n.samples]`, absent when binary |
| `$first.sigma` | `R/plot.R`:12 | `[n.chains, n.burn]`, `[n.chains, 0]` at `n.burn == 0`, absent when binary |
| `$y` | `R/generics.R`:371 | length n.obs |
| `$fit` | `R/generics.R`:107, :191 | only under `keepSampler`; predict refuses anyway |

The census is complete (full scan; everything else stan4bart- or
predict-gated). **Why matrices and not `bart2`'s flat vectors:** both bartCause
readers try to un-flatten with the wrong `byrow` (B0b); storing
`[n.chains, n.samples]` means `is.null(dim(s))` is FALSE at `R/plot.R`:18 and
`R/generics.R`:332, neither reshape fires, and every consumer is correct on its
own terms - `plot_sigma` takes `cbind(warmup, sample)` with chain ROWS
(:29-31, :39-41), `extract.bartcFit` returns the matrix and combines through
bartCause's own helper (:334, `R/generics.R`:1-13), and `sampleFromPPD` gets a
chain-fastest flatten aligned with `ev`'s `[chain, sample, obs]` layout. It is
exactly what `bart2(combineChains = FALSE)` produces (MEASURED).

**Accepted inconsistency, stated:** at `n.chains == 1`, `combineChains(m, 1L)`
returns the matrix (`R/generics.R`:4), so `extract(fit, "sigma")` yields a
`1 x n.samples` MATRIX where a bart fit yields a vector (MEASURED). Harmless -
`fitted.bartcFit` has no `"sigma"` type (`R/generics.R`:252-253) and `summary`
consumes it as `mean(sigma^2)` (`R/summary.R`:98-99, :128-129) - and
`man/bcf.Rd` says so.

**`combineChains` and `n.chains`, settled at the `fitBCF` seam:** `fitBCF`
always returns per-chain arrays and takes no `combineChains`;
`getBCFResponseFit` hands the namedList UNCOMBINED 3-D surfaces, matching
`getBartResponseFit`'s `combineChains = FALSE` extracts
(`R/responseFit.R`:243-245); `bcf()` exposes `combineChains = TRUE` and applies
it at the ACCESSORS; `n.chains` defaults to 4 in `bcf()` (bart2's) and 10 on the
bartc arm (`R/responseFit.R`:221), one per front door.

#### B1.3 The signature

```
bcf(formula, data, subset, weights, offset,          # or (x, y) as bart2 does
    treatment,
    moderators = NULL,                               # tau's split set
    p.score = NULL,                                  # pihat, a MU-ONLY column
    n.trees = 200L, base = 0.95, power = 2.0,        # forest 1 (the fit's own)
    n.trees.treatment = 50L,
    treatment.base = 0.25, treatment.power = 3,
    sd.control = NULL, sd.moderate = 1,
    b.prior.variance = 0.5,
    update.a = TRUE, update.b = TRUE,
    mu.interactions = NULL, tau.interactions = NULL,
    tau.blocks = NULL,
    n.samples = 500L, n.burn = 500L, n.chains = 4L,
    n.threads = dbarts::guessNumCores(), combineChains = TRUE,
    keepSampler = FALSE, verbose = TRUE, seed = NA_integer_, ...)
```

Flat arguments built into the spec internally (bcf-public-surface.md:277-280),
in `bartcoreBCFSampler`'s vocabulary (`R/bartcore.R`:637-653).

**`mu.blocks` is deliberately ABSENT.** MEASURED: a first-forest `blocks` - by
either spelling, since `R/spec.R`:224-226 hoists `forest(1, blocks = )` onto the
top-level argument - is resolved at `R/spec.R`:341 with no `availableColumns`
(doc :340: "The partition covers the full design") and demands a total partition
including the masked treatment column. Auto-repair is worse: `blocks()` carries
`trees.per.group` and NULL distributes trees evenly (`R/model.R`:1514-1515,
definition :1518, work in `resolveBlocks`), so a block of only masked columns
produces trees that can never split. v1 refuses by name, points at `tau.blocks`
(MEASURED to work) and at the door. `mu.interactions` is unaffected (MEASURED).

#### B1.4 The moderator EXCLUSION, and z

- **pihat out of TAU.** `p.scoreAsCovariate` (default TRUE) injects the score
  (`R/bartc.R`:146-159; builders splice at `R/argParse.R`:155-171 decisive :162,
  and :326-337 with the RHS at :346-348). Under BCF it belongs to mu only
  (`R/bartcore.R`:628-630). Its own alignment under `subset` is already correct
  on the literal path (`R/argParse.R`:331-336 scatters into subset POSITIONS).
- **z out of BOTH.** The treatment enters through tau's BASIS; a mu forest free
  to split on z could absorb the effect. bartCause's builders always put the
  treatment in the design (`R/argParse.R`:112, :323-324) and the bookkeeping
  reads it there (`R/responseFit.R`:235), so keep z in the design and mask it
  out of both forests' `vars` - MEASURED to work.

**The propensity column's name is RETURNED, never re-derived.** The builders
already compute it (`R/argParse.R`:141-143, :327-329); append it to the lists
they return (:216, :408) and consume it in the fitter. Bart-path-neutral:
`massign`'s unnamed loop (`R/multipleAssignment.R`:84-98) consumes positionally
and never touches extras, and both existing sites take exactly four unnamed
targets (`R/responseFit.R`:144, :152), as does the stan4bart fitter's (:35,
:45). Re-deriving from `colnames(data@x)` is not an option: MEASURED, with a
confounder named `"psps"` the ladder answers `psps` while the score is `ps`, so
`setdiff` would drop the confounder from tau and leave the score in it - the
exact defect this section prevents - while a membership assert passes. The
resolved name rides the fit as `name.p.score`.

Mechanically: `muVars <- setdiff(colnames(data@x), treatmentName)`,
`tauVars <- setdiff(colnames(data@x), c(treatmentName, p.scoreName))`, or the
caller's `moderators` with those names removed and a refusal - not a silent drop
- if either was explicitly requested.

#### B1.5 The response-fitter contract

`getBCFResponseFit` mirrors `getBartResponseFit`'s assembly half and skips its
test half, reusing the data/literal branch and the two argParse calls
(`R/responseFit.R`:140-155) **with `getBCFResponseFit` substituted at :143 and
:148** - those lines are `addCallDefaults(<call>,
eval(quoteInNamespace(getBartResponseFit)))` and `addCallDefaults`
(`R/utility.R`:112-130) fills from `fn`'s formals, so a verbatim copy would take
the bart fitter's defaults. It also reuses the subset resolution (:161), `trt`
(:235), the common-support computation (:279-282) and the `y` backfill (:284),
and DROPS the counterfactual `x.test` construction (:177-193), the test extract
(:245) and the missing-row recombination (:247-277). Its return is the namedList
contract at `R/responseFit.R`:286, consumed by `assignAll` (`R/bartc.R`
:172-174).

**Four sites need attention:**

1. **`R/plot.R`:12** - a missing `first.sigma`, closed by B1.2a.
2. `refit.bartcFit` (`R/generics.R`:506): widen to `%in% c("bart", "bcf")`.
   SEVERITY: `object$est` is read only under the `p.weight`/`tmle` gate
   (`R/generics.R`:339 -> :346; `R/summary.R`:184 -> :231) and the bart path
   returns `est = NULL` (`R/responseFit.R`:286), so post-refit numbers are
   correct either way; the widening stops a user-reachable element being NULL.
3. `summary`'s group-effects total row (`R/summary.R`:301, `%in% "bart"`).
4. `predict.bartcFit`'s refusal message (`R/generics.R`:102-103), false for bcf.

Plus two `bartc()` guards: `R/bartc.R`:107-112 needs a bcf reading of
`p.scoreAsCovariate`/`method.trt`, and :164-165 must refuse `crossvalidate`.
`method.rsp` gains `"bcf"` in the FORMAL at `R/bartc.R`:3 (the validity loop at
:42-47 reads it); the reserved line :123 becomes
`bcf = redirectCall(matchedCall, quoteInNamespace(getBCFResponseFit))`; :128-130
matches `formals(getBartResponseFit)`, so give `getBCFResponseFit` the same
leading formals or widen that line.

#### B1.6 The five S3 methods

- `print.bartBCF` - call, family, K and each forest's tree count, chains x
  samples, whether a sampler was kept.
- `fitted.bartBCF(object, type = c("mu.obs", "mu.cf", "mu.1", "mu.0", "icate",
  "mu", "tau", "glue", "sigma", "varcount"), ...)` - posterior means (for
  `"varcount"`, per-draw means per forest).
- `extract.bartBCF(object, type = <same>, forest = NULL, combineChains = TRUE,
  ...)` - the draws; `forest` selects one element of the forest-indexed
  containers (`forests`, `varcount`) and defaults to all. `"mu"`/`"tau"` are
  INTERNAL scale, named so in the Rd (which also notes `run()$tau` is an
  unrelated dbarts channel); `"icate"` returns
  `response.scale * (b_1 - b_0) * tau`. Calls
  `issueWarningForUnknownArguments()` first (`R/utility.R`:294-319).
- `residuals.bartBCF(object, ...)` - `y - fitted(object, "mu.obs")`;
  `NAMESPACE`:17 imports `fitted, predict, sd` but not `residuals`, so B1 adds
  it.
- `predict.bartBCF(object, newdata, type = c("mu", "tau"), ...)` - **defined and
  REFUSING in v1**, naming the door.

`NAMESPACE`: `export(bcf)`; `S3method` for print/fitted/extract/predict/
residuals on `bartBCF`; `importFrom(stats, residuals)`.

#### B1.7 Documentation and records

`man/bcf.Rd` (usage, arguments, the two-surface value section, the moderator
exclusion, the per-forest `varcount` container WITH the structural-zero rider,
the `n.chains == 1` sigma shape, and every refusal: missing responses,
`crossvalidate`, `group.by` with `use.ranef = TRUE`, `parametric`, `mu.blocks`,
`predict`); `man/bartc.Rd` gains `"bcf"` at :9-10, :57-66, :174-190;
`man/generics.Rd` gains the bartBCF methods if it is the right home - read it
first. `inst/NEWS.Rd` gains a NEW FEATURES subsection folded into the unreleased
1.0-10 section (:4-12), since 1.0-10 awaits submission after dbarts is accepted
(dbarts `TODO`:892-895); if VD has already shipped it, bump to 1.0-11.

Gates: full testthat under the private library with `NOT_CRAN=true`;
`R CMD check --as-cran` on a tarball staged outside the tree (the S5 mandate
names `--as-cran`, bcf-public-surface.md:515);
`tools:::.build_news_db_from_package_NEWS_Rd` parses `inst/NEWS.Rd`.
rng: NEUTRAL for dbarts.

## Fork A. How does `bcf()` reach the engine? - RESOLVED

**RESOLVED (VD 2026-08-15): (a), the driver. As drafted.** bartCause drives
`dbarts()` + `$run()` itself: four calls, with chains/threads/burn as
`dbartsControl` slots and the per-draw channels arriving batched from the second
`run()`; no per-sweep loop, so the S0/S4 1.04x figure does not apply.

The alternative, (b) - dbarts grows `bart2`'s output channels first and bcf()
becomes one more `redirectCall(matchedCall, dbarts::bart2)` arm - does NOT buy
what it appeared to: missing-row handling, keepTrees prediction and xbart
crossvalidation are unavailable for a BCF under BOTH branches (`R/spec.R`:509,
:486-491; `refuseUndefinedTestFits`; `R/xbart.R`:1-29). Revision 1's memory
argument for (a) is withdrawn (the batched array is materialized in R either
way); the recommendation stands on the R/C division putting the outer loop on
R's side (`docs/design/r-c-division.md`:84-95, :104-107) and on (b)'s
inheritance claims being refuted. Both critiques CONCURRED.

**The (b) content is DOORED, not scheduled** - see "The K-forest batched front
door" under Doors held.

## Fork B. Standalone `bcf()`, a `bartc()` arm, or both? - RESOLVED

**RESOLVED (VD 2026-08-15): (a), both front doors. As drafted.** Exported
`bcf()` plus internal `getBCFResponseFit`, with `R/bartc.R`:123 re-spelled, over
ONE engine-facing core (`fitBCF`) and ONE output class; `getBCFResponseFit`
reshapes into the namedList exactly as `getPWeightResponseFit`
(`R/responseFit.R`:399-403) and `getTMLEResponseFit` (:682-686) wrap and reshape
`getBartResponseFit`. (b) - the bartc arm alone - orphans the five S3 methods
the S5 contract names; (c) - the standalone alone - leaves bartCause's headline
estimator unreachable through its own front door. The two riders are settled in
this plan: the `combineChains`/`n.chains` defaults (B1.2a) and the shared core
carrying `first.sigma` and the resolved propensity-score name (B1.2a, B1.4).

## Fork C. Per-forest varcount - RESOLVED, REVERSING THE RECOMMENDATION

**RESOLVED (VD 2026-08-15): BUILD THE CHANNEL NOW.** The draft recommended
(a) drop from the v1 contract and ticket (c); VD judged the per-draw per-forest
varcount channel **the durable solution other multiforest consumers will want**
and took (c) directly. Consequences, all worked through above:

- **D3** is the resulting slice: two engine overrides, a bridge that already
  sizes off the axis, an R packaging guard, an Rd/NEWS pass, four falsifiers,
  the aliasing audit, and the arc's one baseline re-record.
- **S5's contract now holds VERBATIM** - per-draw mu, tau, glue, sigma,
  per-forest varcount and both counterfactual surfaces - so binding decision 3
  is literal rather than amended, and `bartBCF$varcount` is a forest-indexed
  container (B1.2) rather than a single slab with a documentation rider. The
  structural-zero rider survives, because it describes the MASKS, not the
  channel: masked columns are zero by construction and `man/bcf.Rd` says so.
- The rejected options are recorded so they are not re-proposed: (a) dropping it
  leaves `varcount` reporting one arbitrary forest for every multi-forest model;
  (b) accumulating it in the driver would reintroduce the per-sweep loop the
  batched channels removed, and is unavailable under a bart2-based fitter.

## Fork D (NOT a fork - recorded settled)

`bcf-naming-generalization` (dbarts `TODO`:144-173) does not interact: the
rename surface is internal and the one user-visible leg, `calibrationMapName`
(`R_interface_bartcore.cpp`:4052-4055, call site :4144), returns "two-forest
calibration map", correct at K = 2 - the case `bcf()` fits. D0 completed that
ticket's symbol list (`TODO`:146, :151-159); that was the only contact.

## Pre-registered falsifiers

- **FB0 (B1, load-bearing).** With `seed` set, `bcf()` reproduces a hand-written
  `dbarts(<the same data object>, forests = <the same knobs>)` +
  `sampleTreesFromPrior` + `run(0, n.burn)` + `run(0, n.samples)` sequence
  BITWISE on mu fits, tau fits, glue and sigma. NEGATIVE HALF: perturb the
  sequence and the draws must DIFFER.
- **FB1 (B1).** On an OFFSET-BEARING fixture,
  `response.scale * (a*mu + b_z*tau) + response.shift + offset == run()$train`
  to <= 1e-14 (MEASURED 8.9e-16; without the offset it fails by 0.76 on that
  fixture). Second half: flipping the `b` index moves the surface by exactly
  `response.scale * (b_{1-z} - b_z) * tau`, sign-checked separately on treated
  and control rows; `mu.hat.obs` is `identical()` to `run()$train` before the
  link. At `n.chains = 2L`.
- **FB2 (B1, both halves).** With `p.score` supplied and `keepSampler = TRUE`,
  the treatment forest's counts are zero on the pihat and treatment columns and
  the prognostic forest's zero on the treatment column - asserted through BOTH
  the live read and (after D3) the per-draw container, on the column named by
  `fit$name.p.score`. NEGATIVE HALF: with the masks removed at the same seed at
  least one becomes non-zero.
- **FB3 (B1).** `crossvalidate = TRUE` with `method.rsp = "bcf"` refuses, naming
  the pinned k; the identical call at `"bart"` succeeds.
- **FB4 (B1).** An `NA` response refuses under bcf, naming the door; the same
  data under `"bart"` fits (`tests/testthat/test-03-responseFit.R`:127).
- **FB5 (B1).** `predict` refuses on a `bartBCF` and on a bcf `bartcFit`, each
  naming the replay door.
- **FB6 (B1, integration).** `bartc(y, z, x, method.rsp = "bcf")` then
  `summary()`, `fitted(type = "cate")`, `extract(type = "icate")` and
  `plot_indiv()` run and return finite values of the documented shapes, at
  `n.chains = 2L`.
- **FB7 (B1, binary).** A probit bcf fit reports both surfaces in (0, 1), has NO
  `sigma`/`first.sigma`, `responseIsBinary` is TRUE, and
  `extract(type = "sigma")` errors with the existing binary message
  (`R/generics.R`:322-323). At `n.chains = 2L`.
- **FB8 (B1, refit).** `refit(fit, commonSup.rule = "sd")` moves
  `commonSup.sub` and populates `est`; the negative half rides a comment.
- **FB9 - RETIRED.** It gated D1, which is now a door (below). If the door is
  ever taken, its own plan re-registers it.
- **FB10 (B1, three legs).** `plot_sigma(bartc(..., method.rsp = "bcf",
  n.chains = 2L))` runs and draws the burn-in line correctly;
  `$first.sigma` is `[n.chains, n.burn]` and `$sigma` `[n.chains, n.samples]`.
  NEGATIVE HALF: delete `first.sigma` and the call must error. THIRD LEG: at
  `n.burn = 0L`, `first.sigma` is `[n.chains, 0]` and `plot_sigma` still runs.
- **FB11 (B1, four legs plus a negative).** `bartc(..., subset = s,
  method.rsp = "bcf")` fits for `s` as an INTEGER index, a LOGICAL mask, a
  NEGATIVE index and a CHARACTER row-name vector, `n.obs` equal to the subset
  size each time; and `man/bartc.Rd`:397's workflow runs end to end. NEGATIVE
  HALF: a full-length basis produces dbarts' loud row-count refusal, asserted by
  message; plus the empty-arm guard refusing by name.
- **FB12 (B0b, both halves).** `extract(fit, "sigma", combineChains = FALSE)`
  equals the same model's `combineChains = FALSE` fit `$sigma` chain for chain.
  NEGATIVE HALF: restore `byrow = TRUE` and it must fail.
- **FB13 (D2, five legs).** (i) `dbarts(<a pre-built dbartsData>, forests =
  <declaring a basis>)` refuses by name. (ii) `dbartsSpec(<same>)` STILL BUILDS
  and the four shipped assertions pass unchanged. (iii) The supported
  composition still builds a BCF and reports `forestFits`. (iv)
  `dbartsData(<a dbartsData>, bases = ...)` warns. (v)
  `amplitude.prior.variance` is accepted on a forest whose basis arrived via
  `dbartsData(bases = )` and still refused on a forest with no basis anywhere.
- **FB14 (D3, both halves).** With `n.thin = 1L` and no sweep after the last
  `storeSample`, `run()$varcount[, f, n.samples, chain]` equals
  `sampler$getForestVariableCounts(f)[, chain]` for every forest and chain.
  NEGATIVE HALF: one more sweep breaks the equality, proving a draw history
  rather than a repeated live read.
- **FB15 (D3, both halves).** Every masked column's row in the treatment
  forest's slab is zero in EVERY draw and some unmasked column is non-zero in
  some draw; without the mask, at the same seed, that column is non-zero.
- **FB16 (D3, the re-record's gate).** The widened channel's forest-1 slab is
  `identical()` to the pre-D3 recorded array on every bcf-equivalence scenario -
  a 2-D `p x n.samples` matrix, since every scenario is single-chain
  (`bcf-equivalence.R`:85), so the oracle is `identical(new[, 1L, ], old)`. The
  re-record does not land until this passes.
- **FB17 (D3, regression guards).** A multinomial fit's `varcount` is unchanged
  in shape, dimnames and values (`multinomial-equivalence-1027be5` bitwise), and
  a single-forest fit's keeps its current rank and is bitwise
  (`equivalence-8b047f8b`).
- **FB18 (D3, the `bart2` packaging path).** `bart2(<a dbartsData carrying
  bases>)` at `combineChains` TRUE and FALSE, at one and two chains: the
  packaged shape, the predictor and `forest1..forestK` dimnames, the forest-1
  slab against the `run()` channel it came from, `fit$n.forests`, and the
  synopsis's kept-draw count. NEGATIVE HALF: drop the `fitSynopsis` arm and the
  printed count must double. Full statement in D3.
- **FB19 (D3, the flat-C legacy contract - the BLOCKER-1 detector).** A
  caller-declared count of 1 against a K-forest combiner writes exactly
  `numPredictors` per sample, the prognostic forest's; the reported count never
  exceeds the chain's forest count; and `test-capi.R`:1236-1250 runs clean under
  the ASAN/UBSAN battery with `consumer.c`'s buffer sizing UNCHANGED. NEGATIVE
  HALF: consult the combiner again and the sanitizer leg must abort. Full
  statement in D3.

## Test budget and stop conditions

DENSE-EQUIVALENT (one assertion per line).

- **D2: ~30 assertions** in `inst/tinytest/` (FB13's five legs). STOP: if item
  1's refusal fires on the supported composition or on any `dbartsSpec()` call,
  stop - the predicate is wrong and would break B1 and four shipped assertions.
- **D3: ~42 assertions** (~16 tinytest for FB14/FB15, ~10 `tests/cpp` for the
  shape predicate and FB19's first two legs, ~6 regression for FB17, ~8 tinytest
  for FB18, ~2 in `test-capi.R` for FB19's third). STOP, both halves: no
  tinytest assertion is expected to move (every BCF varcount assertion is
  self-relative and every positionally-indexed one is single-forest or
  multinomial), AND no caller-owned varcount buffer may be reached by a
  multi-forest sampler - `consumer.c`:186-190 via `test-capi.R`:1250, and
  `R/rbart.R`:599/:618. A sanitizer abort in `test-capi.R` is THIS slice.
- **B0: one literal.** STOP: if the regenerated value differs, stop.
- **B0b: ~10 new assertions** plus two re-derived oracles. STOP: no test outside
  `test-12-plot.R` and `test-05-generics.R` may move.
- **B1: ~90 assertions**, ~58 in `tests/testthat/test-14-bcf.R` (FB0-FB2, FB4,
  FB5, FB7, FB10, FB11, the `bartBCF` structural shape), ~10 in
  `test-04-bartc.R`, ~14 in `test-05-generics.R`, ~8 in `test-03-responseFit.R`.
  Fixtures: the existing `inst/common/linearData.R` and `binaryData.R`, plus one
  offset-bearing variation inline for FB1 and one row-named frame for FB11.
  **`n.chains = 2L` for FB1, FB6, FB7, FB10 and every structural assertion**;
  keep an explicit 1-chain arm. STOP: (1) stop at the first oracle that cannot
  be an equality against a hand-written dbarts call - "runs without error" is
  capped at ~8 of the 90; (2) no third fixture data set; (3) no test may depend
  on `tmle`, `lme4` or `stan4bart`.

## Consumer costs (enumerated, never constraining)

- **bartCause**: the arc IS its cost - one export, one class, five methods, one
  `method.rsp` value, two builder returns, plus two pre-existing defects fixed
  (B0b and the four contract sites in B1.5).
- **dbarts**: D0 (landed), D2 (three R-surface guards), D3 (the engine channel
  plus one dead-capture deletion). **None touches `dbarts.h`'s DECLARED
  surface** - D3 edits one doc comment inside `dbarts_results_t`, which the
  token cannot see (`dbarts.h`:95-100: it "covers the entry-point SIGNATURES
  [...] and NOT the layout of the structs those signatures name", pinned by
  `test-capi.R`:56-60) - so the hash does not move and no sister package
  rebuilds. D3's only cross-package surface is the BCF `run()$varcount` shape,
  which no sister package reads (stan4bart drives the flat C API; treatSens and
  bairrtt are single-forest) - and it stays that way BECAUSE of item 5's
  caller-authority pin, which is the fact to check before anyone relaxes it.
- **stan4bart / treatSens / bairrtt**: zero. One cross-reference: stan4bart's
  `bart-args-forests-guard` (`TODO`:11-18) is D2 item 1's twin.

## Doors held (recorded, priced, not scheduled)

- **The K-forest batched front door - SPELLING UNDECIDED.** Whether `bart2`
  reaches a multi-forest model through a flat `forests =` formal (the shape rev
  3 drafted as D1: a `bart2` formal, per-draw `forestFits`/`glue` on the fit
  object, the response transform parked there, and an `extract()` per-forest
  selector) or through formula-level term syntax (`y ~ bart(...)`, stan4bart's
  shape) is UNDECIDED (VD 2026-08-15). The question rides the
  `bart2-argument-consolidation` review (dbarts `TODO`:80-143, which names this
  interaction at :118-125 and records that bart2 carries 52 formals with no
  `control =` among them, :83-84, :98-108). Not scheduled here; if taken, it
  needs its own plan and re-registers its own falsifier. Creation through
  `bart2` already works via `dbartsData(bases = )` (MEASURED); what is missing
  is the OUTPUT channel, so the door is an output slice whatever spelling wins.
  Two riders. (i) D3 item 3 takes ONE piece of that output surface - the
  packaged `varcount` shape for a multi-forest fit, in the shape multinomial
  already ships - because the alternative was a silently wrong array; the door
  inherits it as a fact, and `forestFits`/`glue` remain entirely undecided.
  (ii) The reciprocal pointer needs repointing: `TODO`:119-121 still cites
  "bcf-s5-relocation's D1" plus design-draft anchors that no longer exist; they
  are repointed at this file and this door section.
- **Per-forest saved-tree replay** (out-of-sample mu(x), tau(x)), needed by
  `predict.bartBCF` and by bairrtt (bcf-public-surface.md:660-662). Engine
  slice, unpriced here. LANDED 2026-08-20 at dbarts 63df524e - see the
  landing note at the end of this file.
- **A test treatment vector / test basis**, which would retire
  `refuseUndefinedTestFits` (:666) and is the only thing that would restore
  missing-response handling to bcf(). A modelling decision first. The
  PREDICT-TIME half - the combined blend at new rows for an already-fit
  amplitude coupling - LANDED 2026-08-24 at dbarts 139a1976 (`bases =`, or
  the bart2 `forest()` term auto-route); see the landing note at the end of
  this file. What is still doored is the FIT-TIME channel: a resident test
  basis, `yhat.test`, NA-y rows - the modelling decision above.
- **First-forest `blocks` against its own column mask.** `forestColumns` is
  computed at `R/spec.R`:555-558 inside the `data@bases` branch while
  `blockSpec` is resolved at :341 outside it; the door is hoisting the forest-1
  mask above :341 and passing `availableColumns` when the fit is multi-forest,
  plus the doc sentence at :340. ~10 R + ~8 test. Ships `mu.blocks`.
- **Formula-path bases subsetting in dbarts** (B1.1a route iii). ~15 R + ~15
  test, and it makes a shipped argument's contract ambiguous at equal row
  counts.
- **`bcf()` over a data handle / row-subset view** (bcf-public-surface.md:683).
- **Joint x/y/z `setData` on a BCF sampler.** Stays VD-TIMED AND UNDESIGNED per
  `TODO`:497-500; not designed, priced or used here.
- **Grouped x BCF** (refused at `R/spec.R`:505), which is why bcf() refuses
  `group.by` under `use.ranef = TRUE`.
- **`update.a`/`update.b` mutable mid-chain**; bcf() exposes them at creation.

## Resolutions (VD walkthrough, 2026-08-15)

Recorded so nothing is re-litigated. All eight items settled.

1. **Fork A = (a), the driver.** As specced. (Fork A section.)
2. **D1 = DOORED, not scheduled:** "K-forest batched front door - spelling
   undecided (`forests =` argument vs formula-level term syntax a la
   stan4bart's `y ~ bart(...)`)", riding the `bart2-argument-consolidation`
   review (`TODO`:80-143, interaction at :118-125). (Doors held, first entry.)
3. **Fork B = (a), both front doors.** As drafted. (Fork B section.)
4. **Fork C = BUILD THE CHANNEL NOW**, reversing the drop-and-ticket
   recommendation: the per-draw per-forest varcount channel is the durable
   solution other multiforest consumers will want. (Fork C section; slice D3;
   B1.2's `varcount` element.)
5. **D2 before B1.** (Sequence line; D2 section.)
6. **B0b kept in-arc, immediately after B0.** (Sequence line; B0b section.)
7. **The `predict.bartcFit` propensity-column collision is TICKETED and fixed
   right after B1.** MEASURED live: with a confounder named `"psps"` the ladder
   (`R/generics.R`:110-112) resolves the wrong column and :183 overwrites it
   with predicted scores. The fix records the resolved name on the `bartcFit` -
   a namedList widening at `R/responseFit.R`:286 touching all four fitters plus
   `R/bartc.R`:172-181, ~12 lines across five sites - and is mechanical once
   B1's builder-return change has supplied the name; folding it into B1 would
   put a bart-path contract change inside a bcf slice.
8. **`mu.blocks` refused in v1, door held.** (B1.3; Doors held.)

## Departures and corrections (record)

This plan's facts were adjudicated across three blind critiques; what follows is
what CHANGED under them, kept because each correction is a fact a reader would
otherwise re-derive wrongly.

- The bartCause anchors this arc was scoped from were shifted by up to ten lines
  (the response-fitter contract is `R/responseFit.R`:286); the regression
  snapshots are at `test-06-regression.R`:10, :20, :30, :52.
- Creation of a BCF through `bart2` is OPEN (via `dbartsData(bases = )`); what
  is closed is the OUTPUT channel. The reverse was assumed at first.
- Fork A(b)'s "inherits missing-row handling, keepTrees prediction and xbart
  crossvalidation for free" is REFUTED: all three are refused for a BCF under
  both branches. Fork A(a) was correspondingly over-priced (no per-sweep loop is
  needed), and a memory argument once offered for it is WITHDRAWN (the batched
  array is materialized either way).
- D2 must NOT touch `dbartsSpec()`, whose declared-basis route over a pre-built
  `dbartsData` is supported and pinned by four tinytest assertions.
- B0b is three tokens plus two comments plus two test oracles plus a shipped
  output change: the suite MIRRORS the sigma defect rather than missing it.
- The propensity column's name comes from the builders, never from a ladder
  re-derived against `colnames(data@x)`.
- The subset normalization is total over four kinds, not two.
- `first.sigma` needs a zero-column form at `n.burn == 0`, and both
  counterfactual surfaces are linked LAST, never one before the other.
- D3's consumer census twice stopped short: the flat C API creates BCF samplers
  and owns a varcount buffer (item 5), and `bart2`'s packager breaks in its
  RESHAPE rather than its naming (item 3). The rbart_vi audit was first aimed at
  `R/spec.R`:505, which cannot fire on the path that owns the buffer.
- The re-record carries the standing four-place obligation (item 6), FB16's
  baseline is 2-D because every bcf-equivalence scenario is single-chain, and
  FB19's overflow detector is its `tests/cpp` leg, not its capi leg.

Ten hypotheses were hunted and FAILED, and are recorded so they are not re-run:
K = 1 back-compat holds (the bridge branch is `numVCForests > 1`); the two
overrides suffice; masked columns are structurally zero in the per-draw channel;
the glue contributes no counts; both leak-detector baselines are immovable by D3
(`benchmarks/R/equivalence.R` has no multi-forest scenario, and the multinomial
harness records varcount and would detect a leak); the RNG-neutrality argument
is sound (the varcount block is a pure tree walk); FB16 catches a transposed axis
or a swapped forest mapping; no tinytest or `tests/cpp` assertion moves; no
bartCause reader trips (its only `varcount` token is
`tests/testthat/test-04-bartc.R`:232-234, on the single-forest propensity fit);
and common support does NOT degenerate under BCF (measured at n = 300, 2 chains,
200 draws against a matched BART fit: the two `(sd.cf/sd.obs)^2` distributions
are indistinguishable).

## Open verification obligations

- **The B0 literal has not been regenerated in-tree**; the landing regenerates
  it in-suite under `NOT_CRAN=true` and stops if it differs from the measured
  0.44963899452561451.
- **`man/generics.Rd` was not read in full**, so B1.7's "if it is the right
  home" stays an instruction to the implementer.
- **No probe of the end-to-end bcf path exists**, because the code does not
  exist yet; every measurement is of the dbarts primitives the design composes
  plus bartCause's existing paths.
- **D3's widened channel has not been built or run**; its shape claims are
  derived from the softmax coupling's live behavior at the same seam plus the
  bridge's own allocation branch, and FB16 is the falsifier that would catch a
  wrong derivation before the re-record lands. What HAS been measured is every
  consumer-side claim item 3 and item 5 rest on: `nameVarcount` on a 4-D raw
  array (error at `combineChains = FALSE`, silent 16 x 3 fold at TRUE),
  `shapeMultinomialChannel` on the same array (8 x 3 x 2 / 2 x 4 x 3 x 2 /
  4 x 3 x 2 with the trailing margin named), and today's `bart2` BCF fit
  reporting "kept draws (per chain): 5" at `n.samples = 5L, n.chains = 2L`.
- **The flat-C overflow was reasoned from the write path, not executed** -
  building the override and running the varcount arm under ASAN is the
  implementer's first step, and FB19's negative half (leg (i), the `tests/cpp`
  buffer) is written to reproduce it deliberately.

## Landing notes

ARC CLOSED 2026-08-16. All six slices landed: D0 42846863, D2 511ea2b3, D3
6e3b9fb8 + 9a20bb4a + e8404d93 (dbarts), B0 d7c855f4, B0b e33b6f20, B1
34461ce + 8c07947 (bartCause). Final tips: dbarts origin/bartcore e8404d93
(CI six-green, sanitizers included, at both 511ea2b3 and e8404d93);
bartCause origin/dbarts-1.0 8c07947. Walkthrough resolutions discharged:
Fork A = (a), the driver (four calls, no per-sweep loop); Fork B = (a),
both front doors (exported `bcf()` plus the `bartc()` arm over one
`fitBCF` core); Fork C REVERSED the plan's own drop-and-ticket
recommendation and BUILT THE CHANNEL NOW. Because of that reversal the S5
contract ships VERBATIM rather than amended - per-draw mu, tau, glue,
sigma, per-forest varcount and both counterfactual surfaces, forest-
indexed - with per-forest varcount reaching `bartBCF` through D3's widened
`run()$varcount` channel (`shapeMultinomialChannel`, `n.forests` recorded
on the fit) rather than a documentation rider. Ticketed residues, by
name: `bartcause-subset-pscore` (two pre-existing bartCause bugs surfaced
in gate-running, neither caused by `bcf()`),
`predict-bartcfit-propensity-collision` (VD-resolved to land right after
this closure), `getforestvariablecounts-dimnames` (small, dbarts-side).
Doors held, unchanged by this arc: `mu.blocks` (refusal ships, points at
`tau.blocks`); the K-forest batched front door (spelling undecided, rides
`bart2-argument-consolidation` - `TODO`'s own pointer into this arc cited
gitignored design-session artifacts that no longer exist as citable
anchors and is repointed at this plan's Doors held section instead);
per-forest saved-tree replay (both `predict` refusals name it); grouped
BCF (`setOffset`'s BCF refusal is the real guard).

B1 LANDED 34461ce + trailing 8c07947 (bartCause; S5 proper - `bcf()`,
`bartBCF`, five S3 methods, the `method.rsp = "bcf"` arm, the moderator
exclusion; 14 files +1612/-12). `fitBCF` drives `dbarts(forests =
list(forest(vars = muVars), forest(basis = cbind(1-z, z), vars =
tauVars)))` -> `sampleTreesFromPrior` -> two batched `run()` calls, no
per-sweep loop; exported `bcf()` plus x/y doors; print/fitted/extract/
residuals defined and predict DEFINED AND REFUSING, naming the per-forest
saved-tree replay door. `normalizeBCFSubset` handles all four subset
kinds (logical/int-pos/int-neg/character); an empty-arm basis refuses by
name. The moderator exclusion: builders now RETURN the resolved
`pScoreName` (appended positionally, `massign`-safe); `tauVars =
setdiff(design, c(treatment, score))`; a moderator naming the treatment
or the score is REFUSED, not silently trimmed. `fit.rsp` carries
`sigma`/`first.sigma` as `[n.chains, n.samples]` MATRICES at every chain
count (B0b's reshape branches never fire); `missingRows` is length
`n.obs` all-FALSE, a FIX over the bart arm's full-length vector, which
misaligns with `sd.obs` under `subset` (gate-runner measured). `varcount`
is a forest-indexed container off D3's channel. `man/bcf.Rd` (284 + 8
trailing) carries the structural-zero rider; `test-14-bcf.R` has 173
assertions, all falsifiers both halves; the trailing 8c07947 commit adds
seven bookkeeping \value items, FB tags and the FB6 `plot_indiv` leg
(test-05 122 -> 123).

Battery: suite FAIL 0 PASS 711 (+173, exactly the new assertions; the
parent baseline was rebuilt to prove no existing assertion moved). FB0
(the load-bearing falsifier) is bitwise identical to a hand driver on
mu/tau/glue/sigma/first.sigma/mu.hat.obs/varcount.tau; the negative half
(an extra burn sweep) differs. FB1's offset reconstruction is exact to
8.9e-16 (2.27 without the offset - non-vacuous). `check --as-cran`: the
pre-existing Version/Date WARNING only. Gate-runner CONFIRMED with its
OWN hand driver (seed 4242, x/y interface, identical on every channel),
moderator varcount totals (tau z=0 ps=0, mu ps=5410 z=0; negative half
all nonzero; the "psps" name collision resolves correctly - psps
survives in tau at 109 while ps=0), the bartc arm end to end
(summary/estimands finite, refit moves commonSup 120 -> 84, the binary
probit arm clean, group.by's total row correct), all four subset kinds
(n = 76/70/85/74), and reconstruction to 4.4e-16/8.9e-16. Seven
deviations ALL CONCUR, one count corrected: EIGHT extra fit elements,
not three, all load-bearing and documented. Doors this slice leaves
held, unchanged: `mu.blocks` (refuses by name, points at `tau.blocks`
and the door); component-predict replay (both `predict` refusals name
it); grouped BCF (`setOffset`'s BCF refusal is the real guard). Residues
recorded rather than fixed: two "runs without error" oracles in
`test-14-bcf.R`, under the plan's cap of 8; and the
`predict.bartcFit` propensity-collision fix, VD-resolved to land next
(`predict-bartcfit-propensity-collision`, `TODO`).

B0b LANDED e33b6f20 (bartCause, the per-chain sigma defect - not two
tokens: the suite MIRRORED it). Three `byrow = TRUE` drops
(`extract.bartcFit`'s sigma reshape, `R/generics.R`; `plot_sigma`'s
sigma + first.sigma, `R/plot.R`); two comments corrected from
chain-major to sample-major; two shipped test oracles re-derived
(`test-05-generics.R`:46 now `expect_equal(sigma,
as.vector(t(matrix(fit$fit.rsp$sigma, nrow = n.chains))))`;
`test-12-plot.R`:54-55 with its comment at :31-35). FB12 added: a
cross-check against a `bart2(combineChains = FALSE)` refit that never
round-trips combine/uncombine - proven both halves, reverting `byrow`
breaks it. NEWS BUG-FIXES bullet: `extract(type = "sigma")` now returns
chain-major order, matching every other extracted quantity; previously
multi-chain fits interleaved chains, mis-pairing sigma with the
posterior predictive draws. Suite FAIL 0 PASS 538. Gate-runner CONFIRMED
(pair battery):
per-chain bit-identity vs a never-combined fit (max abs diff 0 0 0 0
across 4 chains; the OLD reshape showed 0.191/0.173/0.197/0.267 -
discriminating power); multi-chain `plot_sigma` clean; NEWS parses;
`check --as-cran` 1 pre-existing WARNING (Version/Date, DESCRIPTION
untouched).

B0 LANDED d7c855f4 (bartCause, restore green against dbarts e8404d93).
Pre-edit receipt: [FAIL 1 | WARN 9 | SKIP 10] default env, the NOT_CRAN
set isolating the single failure - `test-06-regression.R`:52 expected
0.50402744431802, actual 0.44963899452561451, identical to the
231744a0-era measurement and proving D2/D3 draw-neutrality end to end (no
leak). One literal changed, regenerated in-suite (a temporary `%.17g`
print, removed after). After: FAIL 0 PASS 534.

D3 LANDED 6e3b9fb8 (engine) + 9a20bb4a (records/baseline) + e8404d93
(comment/pointer correction); CI six-green at e8404d93, sanitizers
included. The arc's one ENGINE slice: `BCFForestCombiner` overrides
`numVariableCountForests()` -> `numForests_` and `variableCountForest(j)`
-> `j`; `Chain::numVariableCountForests` clamps the combiner's report to
`forests_.size()` (the owed guard - the combiner rounds a one-forest
amplitude spec up to two); `storeSample` reads
`results.numVariableCountForests`; a single caller-count clamp in
`Sampler::run`'s stride computation is exported into the shape at
`facade.hpp`:382 so the R bridge sizes off it. Flat C NEVER sets the
field (default 1), so it keeps byte-identical legacy prognostic-only
bytes; `dbarts.h` gains one doc COMMENT only (token-stream proven, 2870
identical tokens); `apiHash` UNMOVED at 0x85bd1ef04beb3848
(static_assert + a compiled consumer + an independent `R_GetCCallable`
probe). R packaging: a K>1 fit routes `varcount` through
`shapeMultinomialChannel` (`forest1..forestK` on the trailing margin,
predictor names leading), records `n.forests` when K>1, and the
`fitSynopsis` arm reports the true kept-draw count; the K=1 path is
UNTOUCHED (legacy `nameVarcount`). Riders: the dead init-capture
(`treatment = std::vector<double>{}`) is deleted, closing
`dead-bcf-init-capture`; `bartcore_runWithCallback` is pinned to 1 with
an audit comment (the real guard is `setOffset`'s BCF refusal -
`R/spec.R`:505 never reaches the buffer path, measured); six falsified
comments corrected (the plan named four; the implementer found
`chain.hpp`:305-312, the gate-runner found `combiner.hpp`:641-651, the
latter fixed in e8404d93 along with two stale doc baseline pointers,
line counts preserved - `feature-matrix.md` 882,
`multiforest-extension-surface.md` 5324).

Battery: 5318 tinytest FAILURES 0 (+34, additions only); `tests/cpp` 245
ok from clean; ASAN negative half = the designed FB19 abort
(heap-buffer-overflow WRITE of 16 at 0 bytes after a 64-byte region,
`Chain::forestVariableCounts` `chain.hpp`:1256 via `storeSample` :5036);
ASAN positive clean. equivalence 37/37 + multinomial 10/10 BITWISE (the
leak detectors, unmoved). FB16: the old 2-D `varcount` (4x40) vs the new
4x2x40 - forest-1 slab `identical()` on all 12 scenarios, slab 2's last
draw equals the recorded live `varcount.tau`; the harness-vs-old-baseline
mismatch set is EXACTLY the `varcount` channel, 12/12. SHAPE-ONLY
re-record: `bcf-equivalence-6e3b9fb8.rds` (766813 bytes) with the
four-place same-commit bookkeeping - `equivalence.yaml`:87; MANIFEST new
current row :42, old row historical :43, multinomial row shifted :48 ->
:49; TODO ledger "D3 LANDED 6e3b9fb8 + baselines 6e3b9fb8";
`feature-matrix.md` [f39] renamed and repinned to MANIFEST:16,42,49.
Gate-runner CONFIRMED with two independent probes: K=1 A/B
BYTE-IDENTITY vs a base-build library (`varcount`, dimnames, sigma,
`names(fit)` all identical; `n.forests` absent); K=2 raw 4x2x6x2 with
live-read cross-checks and masked-forest structural zeros, packaged
shapes/dimnames/synopsis correct both `combineChains` ways; `R CMD check
--as-cran` Status OK.

Deviations (adjudicated, adopted): two commits recorded, not one (the
baseline names the commit whose behavior it records, repo precedent);
the clamp is spelled in `Chain::numVariableCountForests`; FB18's negative
half pins p (the unarmed 3-D branch reads the predictor count, not a
doubled draw count); a fifth and sixth falsified comment found beyond
the plan's four; the two stale doc pointers fixed in e8404d93. Residue
found in gate-running, ticketed rather than folded in:
`$getForestVariableCounts` (`R/dbarts.R`:1508) returns an UNNAMED matrix
while the packaged channel it mirrors carries predictor names -
`getforestvariablecounts-dimnames` (`TODO`).

D2 LANDED 511ea2b3 (dbarts, three R-surface guards on the construction
seam; CI six-green). `dbarts()` (only) refuses a basis-carrying forest
declaration over a pre-built `dbartsData`, which was a silent
single-forest fit before; `dbartsSpec()` is untouched, since it installs
the declaration rather than discarding it. `dbartsData`'s ignored-args
warning gains `bases`. `amplitude.prior.variance` is excused for
`hasBasis` forests (`R/model.R`'s `excused` hoist moved above the basis
test). Three Rd sentences added: `dbarts.Rd`'s forests item,
`dbartsSpec.Rd`'s shared item stating declaration-REPLACES-bases,
`forest.Rd`'s basis item. 28 assertions over FB13's five legs appended to
`test-bcf-creation.R`. Battery: 5284 tinytest FAILURES 0; the trio all
bitwise; air/lintr clean; `R CMD check --as-cran` OK (two pre-existing
INFOs); `tests/cpp` green. Gate-runner CONFIRMED with two independent
stop-condition probes: the supported composition is unaffected, and
`dbartsSpec()`'s declaration-wins semantics measured (ncol 2 -> 3
replacement).

D0 LANDED 42846863 (dbarts, records/TODO only, no CI fired). Rewrote the
treatment-slot ticket CLOSED against `dbartsdata-treatment-slot-debt`'s
9c63e9d8 with live receipts plus the measured bartCause facts (a
533-assertion NOT_CRAN run, zero argument-passing errors); added the
`dead-bcf-init-capture` ticket (closed by D3 above); completed
`bcf-naming-generalization`'s symbol list (`BCFSpecStorage` plus uses,
`isBCFSampler`/`refuseBCFMutation` plus five call sites); re-pinned
`calibrationMapName` (definition `R_interface_bartcore.cpp`:4052-4055,
the literal at :4054, the call site at :4144); corrected a
tinytest-vs-testthat mislabel.

2026-08-16, post-arc: the bartcause-subset-pscore residue landed at
bartCause 45c7397 (dbarts-1.0). The data.frame + subset misassignment
was fixed by placing the subset-fitted score at the resolved subset
positions of a full-length column (verification found the worse
variant: when the subset length divides nrow(data), the old assignment
silently RECYCLED and the response model trained on a misaligned ps
column). The p.score = <vector> failure was re-diagnosed: bartc has no
p.score formal, so the vector partial-matched p.scoreAsCovariate -
subset-independent; it now fails as a validated length-1 logical with
a method.trt = <vector> hint. FB11's bartc-arm caveat ("holds only via
the literal path") is LIFTED for the data.frame path; supplying scores
as a vector remains spelled method.trt = <vector> by design.

2026-08-16, post-arc: the K-forest batched front door's spelling is
DECIDED by the bart2-argument-consolidation review
(docs/plans/bart2-argument-consolidation.md, section 5): formula-level
terms with colon sugar - y ~ x1 + x2 + z:forest(x1 + x2), general
named form forest(x1 + x2, basis = ~z) - chosen over a flat forests=
formal (declined, kept as a door). The in-sample per-forest output
channel (forestFits/glue packaging, extract(type = "forest")) builds
in that arc as slices S11-S13; measurement there confirmed run()
already emits both channels, so the R-side packaging is the whole
slice. Out-of-sample per-forest replay (predict) remains doored on
the engine-side saved-tree work recorded above.

2026-08-20: the per-forest replay door DISCHARGED at dbarts 63df524e
("Replay each forest separately at new rows"). Engine:
`Chain::predictPerForestFromSavedSample` plus a live-tree twin sum one
forest's trees at caller rows into a forest-major slab, RAW - no
fitScale, no fitShift, no offset - with each forest's own numTrees
driving its loop (ragged counts safe); `Sampler::predictPerForest`
fans over chains and saved slots through one new `SamplerBase`
virtual. Bridge: `bartcore_predictPerForest` gates on
`forestReportingIsDefined` and receives-and-refuses an offset by name
(a per-forest total takes no offset for the same reason it takes no
shift). Surface: `predict(fit, newdata, type = "forest", forest = )`
reusing extract's forest vocabulary on the response scale, and the
sampler's `$predictForests` on the internal scale. The combined
predict on an amplitude coupling stays refused, message now pointing
at replay-and-recombine; `inst/include/dbarts/dbarts.h` untouched; no
state-format bump. 15 files, +651/-18; new test file
test-predict-forest.R so no existing file's RNG history moved.
Gating: independent Opus diff review (verdict LAND-WITH-FIXES; its
one bug - the forest-name margin read a hardcoded axis 3 and returned
a zero-width margin under n.chains > 1 with combineChains = FALSE -
fixed pre-landing with a regression case verified to fail on the
unfixed code, and its sparse-newdata coverage gap closed with a
dgCMatrix case pinned to 1e-12 dense agreement) plus an independent
gate battery: tinytest 6494/0; trio bitwise 42/12/11 identical draws,
zero max|z| lines; tests/cpp full + sampler green; ASAN/UBSAN zero
findings; mutation proofs M1 and M4 re-run independently and
discriminating; R CMD check --as-cran Status OK; air/lintr/NEWS/
doc-freshness clean. The compiled-code gates ran on a pre-amend tree
whose src/ is byte-identical to the landed commit (the fix was
R-and-test only); the R-side gates were re-run on the amended tree.
Known accepted narrowings, recorded not fixed: `predict(type =
"forest")` silently ignores ci.level/weights/n.threads (the spec's
sanctioned bypass), and the `$bc`-handle R wrapper for the bridge
entry existed to plant mutation M3 and was unreachable from the public
families until the multinomial mutation arc replaced the host shell:
`$bc` and the host shell are deleted (multinomial-mutation-arc.md S4).
The binding spec lived at gitignored scratch/
predict-replay-slice-spec.md (critique-amended); this note is the
durable record.

2026-08-24: the combined-predict door DISCHARGED at dbarts 139a1976 ("Blend the
per-forest replay into combined predictions at new rows"). R-only, 792 raw added /
26 removed over 12 files: R/generics.R 243, test-predict-blend.R 371 new,
R/model.R 52, R/formulaTerms.R 38, feature-matrix.md 32, NEWS.Rd 18, R/dbarts.R
14, man/bart.Rd 7, R/bart.R 7, R/data.R 5, bart-as-a-component.md 4,
dbartsSampler-class.Rd 1 - over the 630 stop, accepted on the dense-equivalent
discount: air's one-fragment-per-line formatting of five refusal messages costs
~90 lines for ~25 of message, and the oracle needed six fixture configs, not one.
Interface: `predict(type = "ev"/"ppd"/"bart")` on an amplitude-coupled fit blends
`eta = shift + sum_k (glue_k %*% t(B_k)) * F_k + offset`, link after. New `bases`
formal, after `forest`: length-K list, or a bare value when one forest carries a
basis; refused on single-forest fits and `type = "forest"`; overrides the bart2
auto-route, which derives each basis from `newdata` via `basis.terms` (formula +
xlev + evaluated levels), stashed on `attr(spec$control, "bartcore.forests")`
after `resolveSamplerSpec` and copied to the fit. `expandForestBasis(atPrediction
= TRUE)` skips the all-zero-column/single-level-factor checks - unamended, they'd
refuse the z = 0 counterfactual. `validateForestBases` gains a `rows` noun:
predict names 'newdata', not 'y'. Both amplitude arms now gate on `keepTrees`;
`type = "forest"` lacked it, silently returning one draw before. `offset` is
accepted on the blend, validated after the `type = "forest"` early return, keeping
the C-side `refuseUndefinedTestFits` pin reachable. man/bart.Rd:297's "term's own
text" labels claim is fixed: labels are NULL on the term route. Critique blocker:
bart2's `z:forest(x)` yields a ONE-column basis, not (1-z, z); tests use each
route's own shape. Oracle: in-sample identity vs `yhat.train` at 1e-12, gaussian
(shift/scale non-unit), probit, logistic, two-chain uncombined, q = 3,
forest-1-with-basis; bitwise identity with the documented recombination on both
counterfactual arms, plus the exact arm difference (forest 2's total moved between
amplitudes). Mutations: M1 drop shift, 7 fails (binary green, as predicted); M2
glue positional, 8; M3 basis ignored, 13; M4 double-scale, 9. Gates: tinytest
6721/0; trio bitwise 43/12/11, no max|z|; --as-cran 1 NOTE (days since update);
lintr/air clean; NEWS 294 rows; freshness OK. Anchor pass re-derived ~60
docs/design citations by content; multinomial-mutation-arc.md left as a dated
transcript. Side finding, ticketed, not fixed here: the saved-tree store's
sampler-level write cursor carries across recorded `run()` calls while
`predict`/`predictPerForest`/`predictVariance`/`getTrees` walk slots
0..capacity-1, so a hand-driven `keepTrees` sampler with two recorded runs gets
rotated draws, fixed at 124259d0 (see release-candidate-review.md's burn-down
note). Still open, unchanged: the FIT-TIME test-basis channel (a
`dbartsData` slot, a `forest()` formal, `$setTestBasis`, `yhat.test`, NA-y rows,
flat C entries) is the maintainer-held modelling decision; the sampler-level
`refuseUndefinedTestFits` stands.
