# multinomial-varcounts

agent: opus (a REPORTING arc, not a posterior arc: the per-sample varcount
  channel reads tree state the sampler already produced, so no draw moves. The
  care is entirely in keeping the varcount write/stride byte-identical on the
  single-forest and BCF paths while widening it to K category channels for
  multinomial - the multinomial-test-surface arc's discipline, verbatim).
rng: ONE CLASS: RNG-NEUTRAL. forestVariableCounts counts existing tree splits
  (no draw); widening the storeSample varcount loop, the bridge allocation, and
  the chain stride moves no rng. Every commit bitwise on ALL THREE standing
  anchors: single-forest equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds,
  BCF 5 scenarios x 6 channels IDENTICAL vs bcf-equivalence-99205ee.rds, and the
  multinomial fixture's FOUR existing channels IDENTICAL vs
  multinomial-equivalence-bcefa63.rds (the new per-sample varcount channel is
  additive). A moved draw on any anchor means the varcount seam touched a draw
  path it must not.
window: NONE. Internal engine + bartcore .Call + the bart2 R surface; dbarts.h
  FROZEN (multinomial is internal, stan4bart never runs it); family =
  "multinomial" already ships (bb29d00), so no new family slot.
budget: ~180-230 lines. C1 ~170: engine ~70 (two combiner virtuals + the
  storeSample varcount loop + the chain/sampler stride + a SamplerBase accessor
  + a Results field), bridge ~30 (varcount alloc widened like the C2 fits seam +
  scratch/copy sizing), R surface ~20 (packageMultinomialResults varcount
  reshape), tests ~50 (tinytest shape/names/reproduction + the fixture channel +
  MANIFEST + one tests/cpp per-category storeSample check). C2 ~40 docs. Chiefly
  src/bartcore/combiner.hpp, chain.hpp, sampler.hpp, facade.hpp,
  R_interface_bartcore.cpp, R/bart.R, benchmarks/R/multinomial-equivalence.R + a
  re-recorded baseline, benchmarks/baselines/MANIFEST,
  inst/tinytest/test-multinomial-surface.R, tests/cpp/test_sampler.cpp.
  Header edits -> --preclean; delete stale tests/cpp binaries.

## Goal

The shipped multinomial fit gains a per-sample per-CATEGORY variable-count
channel: bart2(family = "multinomial")$varcount reports each category's forest's
split usage per draw, the symmetric per-category counts mbart2 exposes (its
varcount is a length-K list of ndpost x p matrices). RNG-neutral: the counts
read existing tree state, so single-forest and BCF stay byte-identical. Closes
the follow-up the multinomial C7 landing note filed as omitted (multinomial.md
"The surface"; docs/plans/multi-forest-models.md).

## Context (ground truth, read in code)

- The loss point: storeSample writes results.variableCounts (numPredictors x
  numSamples) for the ONE reportedIndex forest only -
  forestVariableCounts(reportedIndex, results.variableCounts + sampleNum *
  numPredictors) (chain.hpp:2282-2284); reportedIndex = reportedForest()
  (chain.hpp:2211), which is 0 for both BCF and multinomial (combiner.hpp:287,
  not overridden). So a multinomial run's per-sample varcount channel is
  CATEGORY 0's forest only - the "addresses category 0 only" the design note
  cites (multinomial.md). Unlike trainingFits/testFits, the varcount
  channel was NOT widened by the C2 location seam.
- Where per-category counts live TODAY: forestVariableCounts(f, out)
  (chain.hpp:469-474) -> SamplerBase (facade.hpp:160,359) ->
  bartcoreForestVariableCounts (R/bartcore.R:604) reports forest f's CURRENT
  (final-state) cumulative counts - a query, not a per-sample channel; the fit
  keeps it only under keepTrees. The multinomial fixture records this final-state
  query as its "varcount" channel (multinomial-equivalence.R:81-85), NOT the
  per-sample run channel.
- The stride/alloc plumbing to widen, mirroring numReportedLocations exactly:
  Results.variableCounts doc + a new count field (chain.hpp:199,214); the
  per-chain varcount stride (sampler.hpp:314-316, no location factor today); the
  bridge varcountExpr alloc as p x numSamples (x numChains) + scratch + uint32
  copy (R_interface_bartcore.cpp:2177-2182,2216-2217,2242-2244); the callback
  path aliases a per-sweep p buffer (single-forest, 2300-2306).
- Combiner reporting map: numReportedLocations() = K for multinomial (=
  numCategories_, combiner.hpp:589), 1 for BCF/base (297). packageMultinomialResults
  OMITS varcount by name (R/bart.R:743-778); shapeMultinomialChannel reshapes a
  p-leading K-trailing raw array to (n.chains x) n.samples x lead x K (713-732) -
  reusable verbatim for varcount (predictors in the lead slot).

## Design

The widening is per-FOREST, keyed on a combiner-owned forest set, NOT on
numReportedLocations. Add two base virtuals: numVariableCountForests() (default
1) and variableCountForest(std::size_t j) (default reportedForest()).
MultinomialForestCombiner overrides them to numCategories_ and to j (category
k = forest k). storeSample loops j, writing forestVariableCounts(f, base +
(sampleNum * numVCForests + j) * numPredictors); a new Results.numVariableCount-
Forests (default 1) drives the bridge alloc (p x numVCForests x numSamples x
chains, an allocFitsArray-style INTSXP) and the chain stride. At numVCForests =
1 every offset collapses to sampleNum * numPredictors - single-forest and BCF
byte-identical, BCF staying prognostic-only by intent (chain.hpp:190-191), not
coincidence.

ALTERNATIVE (reuse numReportedLocations for the count): saves the one count
virtual since it equals K for multinomial and 1 for BCF/single-forest TODAY, but
conflates "softmax probability output channels" with "forests to report", still
needs a location->forest map, and mis-shapes a future model (heteroscedastic/
hurdle, flagged to re-touch these seams) whose fits-location count diverges from
its forest count. Rejected: the varcount axis is per-forest, mirroring the
existing forestVariableCounts query; the extra virtual is 3 lines.

R surface: packageMultinomialResults reshapes samples$varcount (p x K x
n.samples x chains) through shapeMultinomialChannel to (n.chains x) n.samples x
p x K, levels(y) on the K margin, colnames(x) threaded onto the p margin (pass
predictor names in, as standard bart varcount does, R/bart.R:184-189); drop the
"deliberately omitted" comment. No extract type, no fitted change (varcount is a
fit field). Single-forest/BCF fits UNCHANGED (numVCForests = 1: the R varcount
path is untouched).

## Commits

C1. The per-sample per-category varcount channel, one gated commit (RNG-neutral;
   the engine seam is exercised end to end only through the surface + fixture,
   the multinomial-arc discipline). Engine (combiner.hpp two virtuals + override,
   chain.hpp storeSample loop + accessor + Results field, sampler.hpp stride,
   facade.hpp SamplerBase accessor) + bridge (R_interface_bartcore.cpp alloc +
   scratch/copy + numVariableCountForests set on both run and callback paths) +
   R (packageMultinomialResults varcount) + the fixture channel (record
   res$varcount in recordChannels, re-record multinomial-equivalence-<hash>.rds,
   demote bcefa63 to historical, MANIFEST neutrality note) + tinytest (a
   K = 3 and a K = 2 fit produce a levels+predictor-named n.samples x p x K
   varcount; a public fit's varcount reproduces the internal bartcoreRun channel
   bit-for-bit at a fixed seed - the C7 reproduction pattern) + tests/cpp (a
   per-category storeSample varcount case: channel k equals
   forestVariableCounts(k), the K = 2 reduction). Files above. Gate: all THREE
   anchors IDENTICAL (single-forest 22/22, BCF 5x6, multinomial's four existing
   channels) + the new varcount channel identical to itself + full tinytest (no
   regen - RNG-neutral) + tests/cpp from make clean + air + rchk note
   (bartcore_run alloc change on the bridge). Size: L. --preclean; delete
   tests/cpp binaries. Abort: any anchor or existing multinomial channel moves =
   the seam touched a draw path.

C2. Docs. Correct multinomial.md "The surface" (varcount now DEFINED per-sample
   per-category) and the C7 landing follow-up list; mark the item landed in
   docs/plans/multi-forest-models.md and the TODO; this plan's landing notes.
   Files: docs/design/*, docs/plans/*, TODO. Gate: full tinytest; R CMD check man
   unaffected (no man/ topic added - varcount is a documented bartMultinomial
   field). Size: S.

## Verification

- Every commit: the three standing anchors bitwise (a moved draw = a
  reporting path touched a sampled quantity, stop). No bench-sampler: the
  widening adds work only at storeSample (K forestVariableCounts calls vs 1,
  off the per-sweep partition/suffstat hot path) and only for multinomial;
  single-forest/BCF storeSample is byte-identical. Note it, skip the quiet-window
  bench unless a reviewer disputes the off-hot-path claim.
- C1 additionally: the multinomial fixture's NEW varcount channel identical to
  itself, and the tinytest public==internal reproduction at K = 2 and K = 3.
- dbarts.h unchanged -> no stan4bart lockstep; C1's bridge alloc earns "rchk on
  next scheduled run".

## Open questions for VD

- Q1 RESOLVED (VD 2026-07-16, "internal consistency is fine"): the single
  array (n.chains x) n.samples x p x K, no varcount.mean. The original
  fork, for the record:
- Q1 (varcount fit shape/naming: single array vs mbart2's list-of-K).
  RECOMMEND a single array (n.chains x) n.samples x p x K, levels on the K
  margin, predictors on the p margin - the convention this fit's yhat.train
  already uses (K on the trailing margin), reusing shapeMultinomialChannel, and
  the least-surprise shape for a dbarts user who meets every other K-output as an
  array. AGAINST: mbart2, the closest precedent, returns varcount as a length-K
  LIST of n.samples x p matrices plus a K x p varcount.mean; a user porting from
  BART meets a list. COST of matching mbart2: a list where every other
  bartMultinomial K-output (yhat.train/test) is an array - an internal
  inconsistency a dbarts user meets directly. I weight dbarts-internal
  consistency over cross-package parity here (the array carries the same
  information, sliceable to per-category). Sub-point: standard dbarts bart ships
  no varcount.mean, so I omit it (a user applies mean over the draw margin);
  say if mbart2 parity should override and add $varcount.mean (K x p).

## Landing notes

C1 LANDED (the commit carrying this note): numVariableCountForests /
variableCountForest(j) on ForestCombiner (defaults 1 / reportedForest;
multinomial overrides K / j); storeSample loops the varcount write
forest-major within a sample; Results/Chain/Sampler/facade thread the
count; the bridge allocates p x K x samples (x chains) INTSXP through
an allocVarcountArray lambda mirroring the fits seam, byte-identical
matrix/3D shapes at count 1; packageMultinomialResults reports
(n.chains x) n.samples x p x K with levels on K and predictor names on
p (Q1 resolution). Gates (implementer + independent orchestrator
re-run): 22/22 identical draws, BCF 5x6 bitwise, multinomial neutrality
4/4 existing channels bitwise vs bcefa63; tinytest 2938/0; tests/cpp
clean rebuild; air clean; no new .Call (rchk note: bartcore_run's
varcount allocation widened, same PROTECT pattern as allocFitsArray).
The fixture gains the per-sample runVarcount channel; recorded and
registered in the follow-up commit named by this landing's hash.
