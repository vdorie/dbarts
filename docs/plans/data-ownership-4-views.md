# data-ownership-4-views

agent: opus (the split loop, the facade/chain seam, and the bridge parse
  path are engine-internal; one owner keeps the view contract coherent).
rng: neutral - every new capability defaults OFF. The per-forest column
  restriction is nullptr-guarded (no restriction = the current
  collectAvailableVariables/grow path, byte-for-byte); the handle
  column-subset is a read-only view over shared cuts+codes, identical
  binning to today's full-column view. A consumer that actually restricts
  columns (BCF's tau moderators) changes its own posterior, but that is
  forest-split-bcf's gate, not this plan's. Gate: equivalence 22/22
  IDENTICAL vs equivalence-ac6ec2c.rds at EVERY commit; full tinytest 2771
  from a preclean install, no snapshot regen.
window: plan 4 of docs/design/data-ownership.md (FROZEN 2026-07-11); builds
  on plan-1 container, plan-2 ingestion (6c90507), plan-3 mutation
  (a141338). Plan 5 ships the sparse-categorical kernel + test-side sparse.
budget: ~450-750 lines across ~11 files (src/bartcore/chain.hpp, grow.hpp,
  tree.hpp, facade.hpp, data.hpp, src/R_interface_bartcore.cpp,
  R/bartcore.R, a tests/cpp component test, a mutation/handle tinytest,
  man/*.Rd, the design landing note). Header edits -> --preclean; delete
  stale tests/cpp binaries.
bench: the split loop is hot, but the mask is nullptr-guarded so the
  default path is unchanged. One confirmatory bench-sampler compare vs
  bench-sampler-32fc7c8.csv (single-forest and full-store BCF; no
  regression). Sub-ms metrics carry 6-8% noise; re-run a single flag.

## Goal

A single owned container is shared by several BART models through
column-subset VIEWS: a view is a column-index list over the shared cuts +
codes (kernels read one column at a time - codes.data() + j*n - so no
per-model contiguous block is needed). Two objects share: (a) the forests
of one multi-forest sampler (BCF's prognostic mu reads every column, its
treatment tau reads only the moderator subset), and (b) separate samplers
attached to one standalone handle (CV folds today; hurdle/IRT-style
embeddings next). Sharing is READ-ONLY under the single-writer rule:
predictors are immutable through a shared view (the view holds no
re-quantizable raw, so the mutation and re-quantize surface refuses it -
data.hpp:133 isView, R_interface_bartcore.cpp:1134 refuseViewSampler);
response-side mutation stays per-sampler. Lifetime: a view is a
self-contained copy of the cuts it needs (no back-pointer to the parent
after construction), GC-anchored by the handle's data object via PROT_DATA.

## Binding contracts inherited (plans 1-3, do not reopen)

- data@x is the COLLECTED RAW SOURCE: per-column values R mapped at
  creation, an accepted mutation swapping columns in by R reference at the
  mutating call. NEVER a view of engine quantized state; never written
  during sampling. R reference counting is the ownership tracker.
- The engine keeps NO predictor raw except the leaf-covariate gather
  (data.hpp:193 gatheredRawColumns). The mutable-raw flag was KILLED, never
  built - do not re-add engine-owned mutable columns (that is the ~400 MB
  borrow re-entering by the side door).
- extract(sampler, "predictors") is the on-demand numeric-code
  materializer over data@x; no engine reroute.
- Off-hot-path raw consumers take raw at call time (getTrees replay gained
  its trailing trainingData arg at plan 1; setState/setCutPoints source the
  live R container). dbarts.h is unchanged here; no stan4bart lockstep.
- sparseFactor is exported but refused before the engine (S-CAT, plan 5);
  its CSC-categorical kernel is NOT this plan's.

## Context (what already holds, read in code)

- The standalone handle EXISTS: DataHandle wraps a ColumnStore
  (R_interface_bartcore.cpp:1170); bartcore_createDataHandle:1375 gathers
  every column raw; R helper bartcoreDataHandle (R/bartcore.R:403).
- Row-subset views EXIST: ColumnStore::buildFromParent (data.hpp:742)
  densifies and COPIES the subset's codes over the parent's cut grid,
  gathering leaf-covariate raw with parent standardization; the view is
  self-contained (data.hpp:740). bartcore_createFromHandle:1424 +
  bartcoreSamplerFromHandle (R/bartcore.R:416) build a sampler over one;
  xbart is the consumer (R/xbart.R:552). Views refuse the raw-x mutation +
  re-quantize surface (refuseViewSampler). buildFromParent has NO column
  axis today: it always spans parent.numPredictors.
- Forest is the composable unit (chain.hpp:257); a Chain holds
  std::vector<Forest>; BCF is two forests over ONE shared store. The BCF
  ctor comment records the exact blocker: "both forests read the full
  store (the moderator subset arrives with data ownership)"
  (chain.hpp:471; forest-split-bcf.md step-3 landing).
- Split-variable selection loops over data.numPredictors gated by
  scratch.available[j] (grow.hpp:88-105); availability is filled from
  data + node ancestry by Tree::collectAvailableVariables (tree.hpp:434) -
  there is NO per-forest column axis. Per-forest fixedSplitProbabilities
  already exist and install into treePrior.splitProbabilities
  (chain.hpp:272,442-446), so a per-forest column axis is a natural sibling.
- Design sources: Sharing (data-ownership.md:147-156), open consideration 2
  (:58-71); standalone handle (core-generalization.md:188-192);
  public-surface.md:288 section 5 (DECIDED internal-first; serialization +
  exposure explicitly still open).

## Scope drift to reconcile in the landing note

1. The design's Sharing section says "codes are shared ... one update
   visible to every model." Landed reality: views are self-contained code
   COPIES (buildFromParent, public-surface.md:319 "self-contained copies,
   no lifetime coupling") and REFUSE mutation. This plan reconciles to
   READ-ONLY single-writer sharing (design's own fallback,
   public-surface.md:306-308); shared MUTABLE codes across samplers is the
   open Q3, not built here.
2. The design frames the column subset as a "view" (an index list). Within
   ONE BCF sampler the two forests share ONE store, so the faithful
   realization of "no contiguous block per model" is a per-forest
   AVAILABILITY MASK over the shared codes, NOT a second ColumnStore (a
   second store would duplicate codes and defeat the sharing intent). The
   handle path (separate samplers) keeps the copy-view mechanic.
3. The standalone handle is listed under the design's plan-4 line but
   landed earlier (public-surface.md:314). This plan EXTENDS it (a column
   axis), it does not introduce it.

## Constraints

- Draw-neutral: every new field/argument defaults to "all columns," and
  the default code path is bitwise-identical (nullptr guard on the hot
  loop, not a defaulted-true mask that reorders draws). Equivalence 22/22
  is the gate; any deviation on the default path = defect, stop.
- Single-writer, read-only: a shared view refuses the raw-x mutation and
  re-quantize surface exactly as today (refuseViewSampler). Do not weaken
  that; do not add engine-owned mutable raw (contract above).
- OUT of scope: the sparse-categorical kernel (plan 5); multi-forest
  SAMPLING LOGIC (forest-split-bcf owns BCF's tau update, its R
  `moderators` argument, and the two-forest exact-posterior gate - this
  plan supplies the mechanism the tau forest consumes, see Q2); shared
  MUTABLE predictors across samplers (Q3); handle serialization and public
  exposure (Q4). dbarts.h unchanged.

## Steps

1. Per-forest column-availability mask (engine + facade). Forest
   (chain.hpp:257) gains an optional column list (empty = all). grow.hpp
   and collectAvailableVariables consult it via a nullptr/empty guard so a
   restriction only clears bits, never reweights the default draw. Carry it
   on SamplerOptions per-forest and through createSampler/
   createSamplerOverStore (facade.hpp). Gate: tests/cpp - a forest
   restricted to a subset never proposes a split outside it, AND an
   unrestricted forest's draws are bitwise-identical to a build without the
   field; equivalence 22/22; full tinytest 2771 no regen. Abort: any
   default-path draw moves.
2. Carry the per-forest column list to BCFSpec.tau (bridge + chain).
   BCFForestSpec/BCFSpec (chain.hpp:314,322) gain tau's moderator column
   list; bartcore_createBCF parses it; the BCF ctor (chain.hpp:472) installs
   it on the tau forest. Default (no list) = both forests read the full
   store, byte-for-byte today. The R-facing `moderators` argument + the
   two-forest exact-posterior gate are forest-split-bcf's (Q2). Gate:
   equivalence 22/22 (default BCF unchanged); the existing BCF tinytest +
   bcf-exact.R quick unmoved; full tinytest 2771. Abort: default BCF draws
   move.
3. Column-subset on the handle view path (data.hpp + bridge + R).
   buildFromParent (data.hpp:742) gains an optional column-index list
   (default = all parent columns) so a view spans a subset; the view's
   numPredictors, cut/code copy, leaf gather, and test codes index through
   it. bartcore_createFromHandle + bartcoreSamplerFromHandle
   (R/bartcore.R:416) gain a `columns` argument; xbart passes all columns
   (unchanged). Gate: tests/cpp - a column-subset view bins its columns
   identically to the parent, and an all-columns view is bitwise-identical
   to today's view; equivalence 22/22; full tinytest (xbart snapshots
   UNREGENERATED - all-columns default). Abort: xbart draws move.
4. Docs + landing note. Record the READ-ONLY single-writer resolution, the
   per-forest-mask-not-a-second-store realization (drift 2), the three
   scope-drift reconciliations, and the deferrals (Q3, Q4). Note the
   forest-split-bcf seam: this plan removes the moderator blocker; that plan
   wires `moderators` and gates the posterior. Gate: R CMD check man; full
   tinytest 2771.

## Verification

- Full tinytest 2771 per commit from a preclean install (--preclean;
  header edits in commits 1-3).
- Equivalence 22/22 IDENTICAL vs equivalence-ac6ec2c.rds at every commit -
  the frozen-default gate; any deviation is a defect, stop.
- tests/cpp: restricted-forest containment + unrestricted bitwise identity
  (commit 1); column-subset-view binning + all-columns identity (commit 3).
  Delete stale binaries after header edits.
- xbart, BCF, and mutation snapshots pass UNREGENERATED (every default is
  all-columns; a forced regen means the default path leaked a change).
- One bench-sampler compare vs bench-sampler-32fc7c8.csv at the end
  (single-forest + full-store BCF; no regression - the mask is
  nullptr-guarded). Sub-ms noise caveat.
- rchk on the next scheduled run (commits 2-3 touch the bridge).
- dbarts.h unchanged; no stan4bart lockstep delta (bridge growth is
  internal .Call only, the plan-1/2 precedent).

## Open questions for VD

- Q1 (column-restriction mechanism). RECOMMEND a per-forest availability
  mask cleared in collectAvailableVariables, NOT reusing splitProbabilities
  with zero entries. A hard mask gives exact "this forest may only split on
  these columns" semantics and is trivially neutral when absent; zero
  probabilities entangle the restriction with the DART/fixed-weight draw
  math and the -inf log terms, risking the neutrality gate. What would
  change it: if forest-split-bcf wants moderators expressed as a prior
  weight rather than a hard set.
- Q2 (BCF `moderators` argument + exact-posterior gate: this plan vs
  forest-split-bcf). RECOMMEND forest-split-bcf. That plan owns BCF
  sampling logic and already carries the two-forest gate; this plan should
  stop at the mechanism (steps 1-2 leave tau's list defaulted to all). What
  would change it: if VD wants the moderator surface shipped in lockstep
  with the mechanism rather than as a forest-split-bcf follow-up.
- Q3 (shared MUTABLE codes across samplers - bairrtt's setPredictorJointly
  collapse). RECOMMEND keep read-only single-writer; do NOT build shared
  mutable engine codes now. The design promises "one update visible to
  every model," but that reintroduces shared mutable engine state plans 1-3
  deliberately removed, and bairrtt is an R-API-only consumer whose joint
  update already works per-sampler. What would change it: a measured
  workload where the per-sampler re-quantize is the bottleneck.
- Q4 (handle serialization + public exposure). RECOMMEND defer both, per
  the design's DECIDED internal-first (public-surface.md:310-312,328); keep
  the handle unserialized and unexported. What would change it: a scheduled
  consumer (hurdle, or a bairrtt migration) that needs a persisted or
  user-visible handle.
