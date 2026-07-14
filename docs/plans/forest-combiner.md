# forest-combiner

agent: opus (the combiner is the multi-forest coupling seam - one owner keeps
  per-forest residual formation, the coupling draw, combined fits, and per-forest
  reporting coherent, and every relocation bitwise-neutral on both the
  single-forest and the BCF path).
rng: neutral, EVERY step. Pure relocation of existing math; no draw order moves.
  Single-forest neutrality is guarded by equivalence 22/22 IDENTICAL vs
  equivalence-ac6ec2c.rds. The BCF path is NOT in that anchor (equivalence fits
  through dbarts(), which never builds a BCF sampler - finding below), so it is
  guarded by a NEW bitwise BCF fixture (step 1) at every BCF-touching step, with
  bcf-exact{,-restricted,-weak}.R as the posterior-correctness backstop.
window: none. Internal engine reorg; dbarts.h and the R/bridge BCF surface
  (bartcoreBCFSampler, setTreatment, getForestFits, getBCFGlue) are UNCHANGED.
  Lands BEFORE any multi-forest model (multi-forest-models.md); unblocks
  multinomial/heteroscedastic/hurdle and grouped-x-multi-forest composition.
budget: ~700-1000 lines, chiefly src/bartcore/chain.hpp (the combiner is defined
  beside Forest<L>), plus benchmarks/R/bcf-equivalence.R + baseline,
  tests/cpp/test_sampler.cpp, docs/design/forest-combiner.md, and doc updates.
  Header edits -> --preclean; delete stale tests/cpp binaries.

## Goal

The BCF glue leaves Chain for a polymorphic ForestCombiner<L> (virtual base,
matching the ResponseModel side's dispatch style). The combiner owns per-forest
residual formation, the coupling draw, combined fits, and per-forest reporting
addressing; BCF is its first instance (BCFForestCombiner<L>). Chain keeps the
sweep skeleton and the single-forest fast path: combiner_ is null for one forest,
every touchpoint stays behind `if (combiner_)` (today's `if (bcf_)`), so the
single-forest sweep is byte-identical and pays no new indirection. Chain's sweep
loops forests and calls combiner_->... generically; the two-forest-literal math
(mu/tau, a/b0/b1) is confined to the BCF subclass, so a future combiner subclass
lands without touching Chain.

## Binding contracts inherited (do not reopen)

- dbarts.h FROZEN; the internal .Call + facade BCF surface UNCHANGED. This is an
  engine-header reorg only: SamplerBase's numForests/setTreatment/bcfGlue/
  forestTotalFits (facade.hpp:139-147) and Sampler's fan-out (sampler.hpp:1000-
  1010) forward to Chain exactly as today; Chain forwards to combiner_.
- The Forest split already shipped (forest-split-bcf.md): Chain holds
  std::vector<Forest<L>>, size 1 off BCF; Forest owns trees/fits/residual/leaf/
  selector/keepTrees/columnMask (chain.hpp:271-326). The combiner does NOT touch
  the Forest split; it only relocates the coupling that combines forests.
- Per-forest tree-query addressing is already plumbed: savedTree/savedTreeSlopes/
  savedTreeMasks/flattenTree carry a defaulted forestIndex, numTreesInForest is
  virtual (forest-split-bcf Q-C, f39f335). That plumbing stays; the combiner owns
  the SEMANTIC reporting map (which channel addresses which forest), not the
  tree-query threading.
- BCF's calibration map, glue priors, updateA/updateB switches, and the ridge-
  interweaving move (bcf.md Calibration; bcf-ridge-interweaving.md) are the exact
  math to relocate, not to change. The rng draw order within a sweep is load-
  bearing (below) and is preserved verbatim.

## Context (seams, read in code)

- Chain members to reparent (chain.hpp:2472-2476): forests_ stays; bcf_
  (std::unique_ptr<BCFState>) becomes combiner_ (std::unique_ptr<ForestCombiner
  <L>>). BCFState (chain.hpp:361-369: z, a/b0/b1, aVariance, priors, updateA/B,
  the combined/forestResponse/forestWeights scratch) moves into BCFForestCombiner.
- The combining math, all Chain private methods today, all keyed on forest 0 = mu
  / forest 1 = tau (chain.hpp): forestMultiplier (2283), formForestResponse (2293
  - per-forest residual r/m, weight w*m^2), combinedFits (2308 - null-short-
  circuits to forests_[0].totalFits.data()), drawGlue (2324 - a, then aVariance
  IG, then b0, then b1), interweaveGlueRidge (594 - the post-glue GIG rescale
  reaching forest 0's treeFits/totalFits/test fits/saved-tree FlatNodes).
- Two sweep call sites duplicate the BCF branch and MUST both route through the
  combiner: run() (chain.hpp:733-847: per-forest formForestResponse at 740,
  combinedFits at 831, drawGlue+interweave at 844-847) and growForestFromRoot()
  (chain.hpp:995-1064, the same branch). The default single-forest body of each
  is the byte-identity anchor.
- The sweep's per-sweep draw order to preserve exactly: per-forest tree draws;
  refreshLatents; drawSigma; drawGlue (a, aVariance, b0, b1); interweaveGlueRidge
  (v); per-forest k draw; DART. combiner_->drawGlue and ->afterCombine are called
  at the drawGlue/interweave points, same sequence.
- Per-forest reporting addressing (storeSample, chain.hpp:2366-2456): trainingFits
  is the a*mu+b_z*tau blend under BCF (2379) else forest 0 (2387); testFits (2398)
  and logLikelihood (2445) NaN-flag under BCF (no test treatment vector);
  variableCounts/k/splitProbabilities read forest 0 (2375, 2415-2428). getState/
  setState glue (2/(de)serialize a/aVariance/b0/b1, chain.hpp:1695-1701, 1863-
  1867, 1915-1919; ChainStateData.hasBCF at 261-262).
- BCF-construction helpers stay in the Chain BCF constructor (not the combiner):
  scaledResponseSd (2226), buildBCFForest (2244), the calibration map
  (chain.hpp:528-546). Construction wires combiner_ = a BCFForestCombiner built
  from BCFSpec; the combiner holds the borrowed z and the glue state thereafter.
- GroupedResponse (model.hpp:2639) is a ResponseModel decorator: it subtracts
  group effects from workingResponse() BELOW the combining step. The sweep feeds
  the combiner y = response_->workingResponse() (group-adjusted already), so a
  combiner and a grouped decorator compose without either knowing the other -
  the composition the architecture review flagged as blocked (architecture-
  numerical-review.md:130-133) becomes expressible. This refactor makes it
  expressible; it does NOT build grouped-BCF.

## Constraints

- Draw-neutral on BOTH paths at EVERY step. Single-forest: equivalence 22/22
  IDENTICAL. BCF: the step-1 bitwise fixture IDENTICAL. bcf-exact.R quick/full is
  the correctness backstop but is STATISTICAL (tol 0.015-0.05), not bitwise - it
  cannot alone prove neutrality; the fixture is the bitwise guard. If a step
  cannot be made bitwise on either path, the relocation leaked - stop.
- Single-forest fast path unchanged: combiner_ == nullptr for one forest; every
  touchpoint is a single null test (identical to `if (bcf_)` today); combinedFits
  returns the direct forests_[0].totalFits.data() pointer with NO virtual call and
  NO copy; NO NullCombiner sentinel object (a real object would force per-sweep
  virtual dispatch the single-forest chain must not pay). One confirmatory
  bench-sampler compare vs bench-sampler-4008675.csv at close (quiet window).
- Combiner methods are per-SWEEP granularity, never per-observation - the design's
  "never per-observation virtual dispatch" (core-generalization.md:225) holds; the
  combiner's O(n) loops run inside one virtual call per sweep, as ResponseModel's
  refreshLatents/drawSigma already do.
- fast over safe in C/C++; header-only C++20; Doxygen/LLVM/ASCII; comment/code
  delta bounded.
- OUT of scope: any multi-forest MODEL (multinomial/heteroscedastic/hurdle);
  moving GroupedResponse; a per-observation sigma vector; the BCF-construction
  helpers (they stay Chain-side); dbarts.h / R / bridge changes.

## Scope decisions (recorded)

1. VIRTUAL DISPATCH, templated on the leaf. ForestCombiner<L> is a virtual base
   with BCFForestCombiner<L> the first subclass, selected at CREATION (single vs
   BCF vs future) - a runtime choice at per-sweep granularity, exactly
   ResponseModel's shape. Templated on L because the ridge move reaches Forest
   <L>'s buffers and saved-tree FlatNodes; BCF instantiates only at
   ConstantGaussianLeaf (static_assert as today). Compile-time selection is
   rejected: it would template Chain on the combiner and multiply instantiations.
2. NULL-POINTER SHORT-CIRCUIT, not a NullCombiner. combiner_ == nullptr is the
   single-forest fast path; the sweep branches `if (combiner_)` where it branches
   `if (bcf_)` today. This is what guarantees zero new hot-path indirection.
3. DEFINE IN chain.hpp beside Forest<L>. The combiner needs Forest<L> (defined in
   chain.hpp) for the ridge move, so a separate combiner.hpp would fight the
   include order or force a Forest extraction - churn this neutral refactor should
   not carry. Extract to src/bartcore/combiner.hpp when the SECOND combiner
   (multinomial) lands and gives it a second consumer. "ResponseModel side" is a
   RESPONSIBILITY statement (combining is a response-side concern owned by a
   dedicated object), not a header location.
4. THE COMBINER API GENERALIZES ITS INPUT SIDE BEYOND TWO FORESTS without
   building past BCF: (a) combinedFits(forests) and formForestResponse(f,
   forests, ...) take the whole forest vector, not a hardcoded mu/tau pair;
   (b) formForestResponse already produces a (response, weights) pair per
   forest, so a heteroscedastic variance forest can later route into the
   WEIGHT channel rather than the additive location - the seam exists,
   unused; (c) drawGlue and the post-combine afterCombine hook are virtual,
   so each future model owns its coupling and its move; (d) storeSample asks
   the combiner which forest each reported channel addresses - BCF's
   forest-0 map is the only one built. WHAT STILL RE-CARVES when the queued
   models land (review finding, recorded so multi-forest-models.md plans
   against reality): the combined-fit OUTPUT is a single n-vector
   (combinedFits -> const double*; refreshLatents/drawSigma take one
   location; results.trainingFits is one channel) - multinomial's n x K
   combined object forces signature changes there; Chain holds ONE
   response_/sigma_ and a vector<Forest<L>> with a single leaf type L -
   heteroscedastic's positive-leaf variance forest and hurdle's per-forest
   response families break those CHAIN-level invariants, not the combiner
   API; and the state wire format (ChainStateData's BCF-shaped
   a/aVariance/b0/b1) needs a format bump for any non-BCF combiner. This
   plan freezes the input-side shape and the virtual coupling hooks; the
   named Chain-level invariants are the honest remaining work.

## The BCF path's bitwise gate (finding, load-bearing)

equivalence.R's scenarios (friedman, probit, weighted, splitprobs, chik, ...) all
fit through dbarts(), which builds only single-forest samplers; BCF is reachable
only via the internal bartcoreBCFSampler. So equivalence-ac6ec2c.rds guards the
single-forest paths and says NOTHING about BCF. bcf-exact.R runs the BCF sampler
but matches a closed-form posterior to MC error - statistical, not bitwise. A
neutral refactor of the BCF path therefore has no existing bitwise guard.
test-bcf.R's "default neutrality" echo (test-bcf.R:222) is a within-session pair,
not a before/after-refactor gate. Step 1 establishes the missing guard.

## Steps

1. BCF bitwise fixture + cpp gate tightenings (gate enablers; NO engine
   change). Add benchmarks/R/bcf-equivalence.R in equivalence.R's
   record/compare idiom, mirroring its settings guard (pin seeds, draw
   counts, tree count, n.threads = 1 in the baseline meta; error on
   mismatch). Scenarios: default two-forest; one restricted-moderator; one
   with updateA/updateB toggled; one WEIGHTED (weights ~ U(0.5, 2) - the
   w*m^2 and wi weight channels of formForestResponse/drawGlue are
   otherwise unobserved by any gate, review finding); one setTreatment
   scenario (run, setTreatment(new z), run, record - the only bitwise guard
   on the step-4 setTreatment routing). EVERY scenario records the RAW
   per-forest fits (bartcoreForestFits, both forests), the glue
   (bartcoreBCFGlue), sigma, AND the reported result$train and
   result$varcount channels: train/varcount are the ONLY bitwise guard on
   storeSample's BCF branches, whose glue reads get mechanically rewritten
   in step 2 when combiner_ replaces bcf_ (the live fits/glue accessors
   bypass storeSample entirely, review finding); train also catches a
   pre/post-interweave a-vs-mu mismatch the live state is self-consistently
   blind to. Compare asserts BITWISE identity (identical(), not tolerance).
   Do NOT bitwise-compare a getState -> setState continuation: restore is
   structural by contract (test-bcf.R), and a bitwise assertion there is a
   false gate. Freeze the recorded-column set BEFORE capturing the
   baseline, and confirm equivalence 22/22 still holds at the recording
   commit so the BCF baseline does not bless a drifted engine. Also close
   the two component-gate holes the fixture cannot reach: tests/cpp gains a
   BCF growForestFromRoot case (that branch is unreachable from R - without
   it the plan relocates code no gate observes) and
   testBCFInterweaveKeepTrees is tightened to numThin > 1 (the sampleNum ->
   saved-slot addressing of the saved-tree rescale is pinned only at
   numThin = 1 today). Record the baseline from the pre-refactor HEAD - it
   is the BCF analog of equivalence-ac6ec2c.rds and the per-step guard for
   steps 2-4. Files: benchmarks/R/bcf-equivalence.R,
   benchmarks/baselines/bcf-equivalence-<hash>.rds,
   tests/cpp/test_sampler.cpp. Tests: the script's own record-then-compare
   round-trips identical; the new/tightened cpp cases pass pre-refactor.
   Gate class: RNG-neutral (test-only) - equivalence 22/22 + tinytest
   unchanged + tests/cpp (grown) + the new fixture identical to itself.
   Size: M (grown from S by review).

2. ForestCombiner<L> base + BCFForestCombiner<L>; combining relocated (engine,
   BYTE-IDENTICAL). Introduce the virtual ForestCombiner<L> (formForestResponse,
   combinedFits, plus the reporting/serialization virtuals stubbed for steps 3-4)
   and BCFForestCombiner<L> holding BCFState's fields. combiner_ replaces bcf_
   (scope decision 2). Move forestMultiplier/formForestResponse/combinedFits into
   BCFForestCombiner; Chain's combinedFits becomes `combiner_ ? combiner_->
   combinedFits(forests_) : forests_[0].totalFits.data()`; both sweep loops
   (run() and growForestFromRoot()) call combiner_->formForestResponse behind
   `if (combiner_)`. drawGlue/interweave stay Chain methods reading combiner state
   for now. NOTE (review finding): renaming the owner forces EVERY glue
   reader - storeSample, getState/setState/installForest, setTreatment,
   bcfGlue - to be mechanically rewritten to read through combiner_ in THIS
   step; their semantics move in steps 3-4 but their reads change here, and
   the fixture's train/varcount columns are what makes that rewrite honest.
   The BCF Chain ctor builds combiner_. Files: chain.hpp. Tests:
   tests/cpp BCF sanity unchanged (testBCFTwoForest and kin, test_sampler.cpp:
   1445+). Gate class: RNG-neutral - equivalence 22/22 + BCF fixture identical +
   tests/cpp + bcf-exact.R quick + full tinytest. Size: L. Abort: any single-
   forest fit moves (relocation leaked into the fast path) or the BCF fixture
   moves.

3. Coupling draw + post-combine move into the combiner (engine, BYTE-IDENTICAL).
   drawGlue and interweaveGlueRidge become virtual ForestCombiner methods
   (base: drawGlue no-op, afterCombine no-op); Chain calls combiner_->drawGlue
   (rng_, sigma_, ...) then combiner_->afterCombine(forests_, record, sampleNum,
   rng_) at the exact points and in the exact rng order it calls drawGlue/
   interweave today (a, aVariance, b0, b1, then v). The ridge move's Forest and
   saved-tree writes move with it, taking forests_ by reference. Files: chain.hpp.
   Tests: tests/cpp testBCFInterweave{,KeepTrees} unchanged (test_sampler.cpp:
   1588, 1665). Gate class: RNG-neutral - equivalence 22/22 + BCF fixture
   identical + tests/cpp + bcf-exact.R quick+full + bcf-exact-restricted.R +
   bcf-exact-weak.R + tinytest. Size: M. (Foldable into step 2 if the combined
   diff stays reviewable - the forest-split precedent.) Abort: any bcf-exact gap
   over tolerance or the fixture moves = the draw order or the ridge write shifted.

4. Reporting, serialization, and mutation through the combiner (engine, BYTE-
   IDENTICAL). storeSample asks combiner_ for the combined training fit and for
   channel definedness (BCF: testFits/logLikelihood undefined -> NaN; forest-0
   varcount/k/splitProbs the combiner's reported-forest map); getState/setState
   AND installForest (the warm-start glue restore, chain.hpp:1863-1867 - the
   third glue site, review finding) delegate glue (de)serialization to
   combiner_->serialize/restore (ChainStateData glue fields stay the wire
   format, filled by the combiner); setTreatment and bcfGlue route through
   combiner_. Single-forest storeSample/getState/setState
   keep their direct forest-0 path behind `if (combiner_)`. Files: chain.hpp.
   Tests: tests/cpp state round-trip / fuzz unchanged; test-bcf.R unchanged. Gate
   class: RNG-neutral - equivalence 22/22 + BCF fixture identical + tests/cpp
   (incl. the state fuzzer) + bcf-exact.R quick + full tinytest 2832 no regen.
   Size: M. Abort: a BCF state round-trip or the fixture moves.

5. Design note + anticipation record + confirmatory bench (docs + gate).
   docs/design/forest-combiner.md: the ForestCombiner responsibilities, the null-
   short-circuit fast-path guarantee, and - WITHOUT building - exactly what the API
   anticipates for multinomial (K forests, PG coupling), heteroscedastic (variance
   forest into the weight channel; the deferred per-observation-sigma decision),
   hurdle (two response-linked forests), and grouped-x-multi-forest (GroupedResponse
   composes below the combiner). Record honestly (review findings): the
   Chain-level re-carves scope decision 4 names (combined-output shape,
   single response_/sigma_, single-L forest vector, BCF-shaped wire format);
   that post-mutation fit rebuild (applyNewData/revalidateTrees/
   rebuildFitsFromParameters) remains forest-0-only, a pre-existing
   forest-split interim this plan does not close; and that the combiner is a
   standalone ForestCombiner<L> hierarchy rather than a ResponseModel
   subclass, the accepted reading of the architecture review's debt #1.
   Record that BCF is the first instance and that
   Forest<L>/combiner extraction to combiner.hpp is deferred to the second
   combiner. Update multi-forest-models.md (blocker discharged) and architecture-
   numerical-review.md (debt #1 closed). Files: docs/design/forest-combiner.md,
   docs/plans/*. Gate class: R CMD check man unaffected; one confirmatory
   bench-sampler compare vs bench-sampler-4008675.csv (single-forest hot path;
   quiet window) - no metric regresses. Size: S.

## Verification

- Every commit: equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds (single-
  forest neutrality) + the bcf-equivalence fixture IDENTICAL (BCF neutrality) +
  full tinytest 2832 from a preclean install, no snapshot regen + tests/cpp
  (delete stale binaries after the header edits).
- BCF-touching commits (2-4): bcf-exact.R quick+full UNMOVED (max gaps as recorded
  in forest-split-bcf's step-5/C2 landing notes); bcf-exact-restricted.R and
  bcf-exact-weak.R pass where the step touches the restriction/prior-glue path.
- Close: bench-sampler compare vs bench-sampler-4008675.csv, single-forest arms at
  parity (the fast-path guarantee); dbarts.h unchanged, no stan4bart lockstep, no
  rchk (the bridge is untouched - this is engine-header-only).

## Open questions for VD

RESOLVED (VD, 2026-07-14): Q1-Q4 as recommended. Two independent design
reviews (architecture fit; gate sufficiency) then amended the plan, both
returning "implement with named amendments/additions": scope decision 4
softened to the input-side generalization with the future Chain-level
re-carves named; step 1 grew the fixture (weighted and setTreatment
scenarios, train/varcount reported channels, the settings guard, the
baseline pitfalls) plus the two cpp gate tightenings (growForestFromRoot
BCF case, keepTrees numThin > 1); step 2 records the mechanical
glue-reader rewrite it forces; step 4 names installForest's glue restore.
Construction detail from review: the combiner needs the observation count
(and the ridge move data_) - hold a ColumnStore reference or receive n at
construction. One reading AWAITING VD confirmation (recommended: accept): the combiner
is a standalone ForestCombiner<L> hierarchy, not a ResponseModel subclass
- architecture-numerical-review.md:130-133's "ResponseModel side" is read
as a responsibility statement (combining belongs to a dedicated
response-side object), since ResponseModel's per-observation-location
interface does not fit per-forest residual formation. Steps 2-4 depend on
this reading; step 1 does not.

- Q1 (combiner dispatch). RECOMMEND virtual ForestCombiner<L>, null-pointer
  short-circuit for single forest. It matches ResponseModel (creation-selected,
  per-sweep, runtime-virtual), and the null pointer keeps the single-forest fast
  path byte-identical with no virtual call. What would change it: nothing short of
  measuring a per-sweep virtual call as non-negligible, which it is not (the
  ResponseModel calls already dominate).
- Q2 (GroupedResponse placement). RECOMMEND leave it a ResponseModel decorator,
  unmoved. The combiner feeds the response chain the combined fit; the grouped
  decorator sits below and subtracts group effects from the working response the
  combiner then reads - the two compose already, which is precisely what making
  combining a response-side object unblocks. What would change it: if a future
  model needs group effects computed AFTER combining rather than before, revisit;
  BCF and the queued models do not.
- Q3 (does the coupling draw's rng order constrain the combiner call order). YES,
  and the plan honors it: combiner_->drawGlue/afterCombine are called at the same
  sweep points and draw in the same sequence (a, aVariance, b0, b1, v) as today,
  so BCF stays bitwise. RECOMMEND the combiner own its internal draw order and
  Chain fix only WHERE in the sweep drawGlue/afterCombine fire (unchanged); a
  future model defines its own internal order freely. What would change it: a
  model needing a draw interleaved among the tree sweeps (not after) would need a
  finer combiner hook - none queued does.
- Q4 (combiner header location). RECOMMEND define in chain.hpp beside Forest<L>
  now; extract Forest<L>+ForestCombiner<L> to combiner.hpp when multinomial (the
  second combiner) lands. Header extraction is pure motion best done with a second
  consumer to shape it. What would change it: if chain.hpp's size becomes a
  review burden before multinomial, extract earlier as its own neutral commit.
