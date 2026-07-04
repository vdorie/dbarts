# Linear leaves

Wave-2 model, phase 5 of core-generalization.md. The first leaf model
beyond the constant leaf: each bottom node holds a small Bayesian linear
regression, so a forest fits smoothly-varying coefficients
(Chipman-George-McCulloch treed regression; the varying-slope structure
of Bayesian Causal Forests when the leaf covariate is a treatment
indicator).

ALL FOUR STAGES LANDED 2026-07-04; xbart view support remains the one
deliberate deferral. The four open decisions at the end of this doc were
implemented per their recommendations (explicit column designation, slope
sd = intercept sd, chi-k over all coordinates, beta.<name> reporting) -
any can be revisited before release. Stage-4 (R surface) notes:

- node.prior = linear(columns, k) on dbarts(): a dbartsLinearPrior
  whose raw designation (names or indices) resolves against the model
  matrix in parsePriors (resolveLeafCovariates: exact column names or
  1-based indices, duplicates and factor columns rejected via
  data@varTypes, at most 8). An unresolved designation cannot enter
  dbartsModel directly, mirroring the split.probs guard. bart2 and
  rbart_vi construct node.prior = normal(k) internally and do not
  reach linear; xbart has no node.prior at all, and the data-handle
  creation path refuses designations outright.
- The bridge dispatches through bartcore::createSampler (both
  instantiations now ship; speed compare stayed at-or-below baseline,
  .so 442 -> 529KB total including the format code). parseModel reads
  the resolved columns off the prior; setModel refuses a designation
  change (the leaf model is a template instantiation, fixed at
  creation).
- State objects gain tree.params/saved.params slots (REALSXP, each
  tree's slopes concatenated, pre-order by leaf; split by the (m+1)/2
  leaf counts of tree.sizes). Linear-leaf states require them.
- getTrees emits one slope column per covariate (NA on internal
  nodes), generically named C-side and renamed beta.<column name> by
  the R wrapper (open decision 4, implemented as recommended);
  plotTree labels linear leaves "b0 + b1 name ...". printTrees prints
  "b:" slope lists at leaves.
- Gates: tinytest test-linear-leaves.R (17 asserts: designation
  validation, slope recovery, saved-tree replay vs recorded test
  fits, beta columns, bitwise state round trip through the R slots,
  fixed-designation refusal, forced predictor mutation, probit
  smoke); full suite 2209/60; equivalence identical draws; speed
  clean; R CMD check --as-cran.

Stage-3 deltas:

- Each live tree persists its parameter blocks (Chain::paramsByTree_,
  arena-indexed, stride numParams) - the source of truth for flatten,
  state, prediction, and the mutation flows, replacing recovery from
  fits. Draws land there directly; sampleTreesFromPrior zeroes them.
- FlatNode stays {variable, value, flags} with the intercept in value;
  slopes ride side arrays, numParams - 1 doubles per leaf in pre-order:
  Tree::flatten/buildFromFlat/collapseEmptyNodes/mapOldCutPointsOntoNew
  and the print helpers take a defaulted paramStride (and slope array),
  so every scalar call site compiles the identical arithmetic (verified:
  equivalence identical draws). ChainStateData gains treeParams and
  savedTreeParams; stateIsValid requires (m + 1) / 2 slope blocks for an
  m-record tree.
- Prediction replays flat trees through addFlatLinearPredictionsBelow,
  standardizing designated columns on the fly with the training
  constants - bit-identical to the engine's own uTest math, so saved-
  tree predictions match recorded test fits.
- Mutation semantics: in-place predictor changes (setPredictor,
  updatePredictor, per-observation sweeps, setCutPoints) regather the
  covariates under the EXISTING standardization constants (calibration
  is sticky, like refreshCutsForColumn keeping the cut count); rollback
  regathers the restored values exactly. Whole-data setData
  re-initializes constants the way it rebuilds the cut grid, carrying
  the persisted parameters across as the same approximate continuation
  the split remap embodies. Collapse merges average per coordinate.
- Every stage-2 refusal is lifted; xbart views remain the stage-4 item
  (ColumnStore views hold no raw covariates yet), and the bridge still
  instantiates only the constant leaf.

Stage-2 notes follow.
The leaf concept split into LeafModelCore plus Scalar/VectorLeafModel
shapes: the moves see only logIntegratedLikelihoodForNode, the constant
leaf keeps its scalar draw interface and code paths textually unchanged
(gates: equivalence identical draws, suite 2192/0, speed within noise),
and vector models write numParams() doubles per leaf with fits evaluated
per observation; Chain branches on L::hasVectorParams at compile time.
Deltas from the proposal discovered while landing:

- No arena-indexed sufficient-statistic blocks: (U'WU, U'Wz, z'Wz)
  accumulate per call over the node's index segment instead, the same
  O(n_leaf (q+1)^2) the Costs section already budgeted per score. Cached
  blocks would need snapshot/restore coupling with the moves' subtree
  rollbacks (Node's two scalars ride the existing snapshots; side
  storage would not) for a constant-factor win; revisit only if
  profiling demands.
- Missing leaf-covariate values enter at the standardized mean (zero)
  rather than erroring: composes with MIA, whose rules still route the
  missingness itself. A constant or all-missing designated column keeps
  sd 1 and degrades to an extra intercept the ridge absorbs.
- Posterior and prior draws fill coordinates in order, intercept first;
  empty leaves zero their block without consuming generator draws,
  matching the scalar path.
- Each chain owns its standardized covariate copy (n x q doubles);
  share through the sampler later if multi-chain memory matters.
- Stage-2 refusals (until the stage-3 formats): keepTrees, getState/
  setState, flattenTree/printTrees, predict from live or saved trees,
  and the whole raw-x mutation surface (setPredictor, updatePredictor,
  the per-observation sessions, setData, setCutPoints) refuse
  gracefully - rolled back, false, empty, or zeroed, never wrong.
- createSampler (facade.hpp) dispatches on options.leafCovariateColumns
  and validates the designation (at most 8 columns, in range, not
  categorical; null on failure). The bridge stays on
  createClassicSampler until stage 4: linking the second instantiation
  costs +65KB .so (442 -> 506KB, +15%), +0.5s on the bridge TU (2.9 ->
  3.8s), and its code-layout shift alone moved hot microbenchmarks 3-9%
  with no behavioral change, so the cost should land with the feature.
- Component tests check the marginal against an independently coded R
  reference (weighted/unweighted, children of a split, q = 2, and q = 0
  equality with the constant leaf), posterior draw moments against R,
  and end-to-end varying-slope recovery (unit slope on the active side,
  zero on the flat side) plus a sampled-k smoke test.

Original proposal follows.

## Model

A designated, small set of q predictor columns u(x) enters every leaf:
leaf fit = b0 + b'u(x), with conjugate normal priors on (b0, b), so the
marginal over the parameters stays closed-form and the existing
conjugate MH moves apply unchanged. The constant leaf is the q = 0
special case and keeps its own code path untouched.

- Leaf covariates are ordinal (continuous) columns only; a categorical
  column in the set is an error (its codes are unordered - interact via
  splits instead). Leaf models consume raw column values, never codes.
- Internally the leaf covariates are standardized (training mean and sd,
  stored on the store) so the prior is expressed on a comparable scale;
  reported parameters stay on the internal standardized scale like leaf
  values today.
- Priors: b0 ~ N(0, (scale / k)^2) exactly as the constant leaf;
  b_j ~ N(0, (scale / k)^2) on the standardized covariates. scale is
  node.scale / sqrt(numTrees) as today. No cross-coordinate prior
  correlation; the posterior is the usual ridge normal with
  V = (U'WU / sigma^2 + P)^-1 solved by Cholesky of a (q+1) x (q+1)
  block.

## Engine

- IntegrableLeafModel widens from scalar arguments to node context: the
  leaf model scores a node from (data, tree, nodeIndex, y, weights, k,
  sigma^2) and computes/refreshes its own sufficient statistics.
  ConstantGaussianLeaf's implementation reproduces the current
  average/numEffectiveObservations/variance sequence exactly - the
  bitwise gate on constant-leaf fits.
- Node keeps its two cached doubles (the constant leaf's stats, also
  used by the empty-leaf veto and orphanChildren merging). The linear
  leaf maintains arena-indexed (U'WU, U'Wz, z'Wz) blocks in its own
  side storage, rebuilt where setNodeAverages runs today and refreshed
  by birth/refreshSubtree. Accumulation is one pass over the leaf's
  index segment reading q raw columns.
- Draws stay at tree granularity (sampleParametersAndSetFits): per leaf,
  posterior mean/Cholesky solve plus q+1 standard normals. paramByNode
  becomes (q+1) doubles per node (flat vector, stride q+1). Fits change
  from set-to-constant to a dot product per member observation; test
  fits and predict route as today then evaluate b0 + b'u.
- The parameter-recovery-from-fits flows (predictor mutation, setData
  remapping, state restore) assume fits are constant within a leaf.
  Linear-leaf chains instead persist their parameters explicitly
  (per-tree paramByNode retained between sweeps, T x nodes x (q+1)
  doubles - small) and the recovery helpers dispatch through the leaf
  model.
- The k (chi) hyperprior update accumulates sum-of-squared parameters
  over every coordinate, with the leaf count scaled by q+1 - the same
  scaled-chi posterior since all coordinates share the scale/k prior sd.
  (Alternative if unwanted: refuse updateK with linear leaves at first.)
- Moves, tree priors, DART, and split machinery are untouched: they are
  already templated on the leaf model and only see branch scores.

## Formats

- FlatNode stays {variable, value, flags}. A leaf keeps its intercept in
  value (constant-leaf consumers degrade gracefully); the slopes live in
  a per-tree params array parallel to the pre-order leaves, appended to
  state as tree.params/saved.params (REALSXP, q doubles per leaf in
  pre-order, required for linear-leaf states, absent otherwise). Replay
  and prediction helpers take the params alongside the flat nodes.
- getTrees: leaf rows keep the intercept in value and gain one numeric
  column per leaf covariate, named beta.<column name>, NA on internal
  nodes (and on everything for constant-leaf samplers - the columns
  appear only when the model has linear leaves). plotTree prints
  "b0 + b1 u1 + ..." with signif()-rounded coefficients.
- Views (the xbart data handle) hold no raw x, but linear leaves need
  the raw leaf covariates: buildFromParent additionally gathers those q
  columns (plus their standardization constants), keeping views
  self-contained. Everything else about views is unchanged.

## R surface

- The leaf model is the node prior: dbartsPriors gains
  linear(columns, k = 2), accepted by node.prior everywhere normal(k)
  is; bart2 and rbart_vi reach it through their existing prior-object
  plumbing. columns is a character or integer vector naming predictor
  columns; with factors as categorical a factor column is a single
  (rejected) column, while under indicators its dummies are ordinary
  ordinal columns and legal.
- xbart: linear leaves arrive through node.prior once views gather raw
  leaf columns; until then xbart rejects them with a clear error.
- Weights, probit, and logistic compose automatically: the working
  response and Polya-Gamma weights flow through U'WU like any weighted
  regression.
- bart keeps the frozen shim (constant leaves only).

## Costs

Per proposed move the branch score costs O(n_leaf q^2 + q^3) against
O(n_leaf) today; q is expected to be 1-4 (a treatment arm, a running
variable). Sampler<LinearGaussianLeaf> is a second template
instantiation; compile time and .so size are the bounded risk from
core-generalization.md and get measured at landing.

## Gates

Constant-leaf behavior must be bitwise unchanged: full suite, the
equivalence compare (identical draws), speed benchmarks. New coverage:
a fixed-tree marginal likelihood and posterior draw checked against
independently coded R math; end-to-end recovery of a varying slope
(y = x2 * step(x1) + noise, leaf covariate x2); state round trip with
params; predict/getTrees consistency; probit and logistic smoke tests;
the chi-k accumulation if kept.

## Open decisions (Vincent)

1. Leaf covariate designation: this proposal says an explicit column
   set on the prior (linear(columns, k)), no all-columns default -
   right call? (All columns is reachable by listing them; with p large
   it is the expensive, original CGM-2002 formulation.)
2. Slope prior scale: standardized covariates with slope sd = intercept
   sd (scale/k). A slope.scale multiplier could be exposed later
   without breaking anything - start without it?
3. k hyperprior over all coordinates versus fixed-k-only at first.
4. getTrees reporting as beta.<name> columns - or would a single
   list-column or a separate attribute be preferable?
