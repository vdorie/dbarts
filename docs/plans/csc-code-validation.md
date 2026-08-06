# csc-code-validation

status: LANDED df79f17 2026-08-05, all CI green incl. sanitizers
agent: opus
rng: neutral
budget: ~400 lines (creation validator rewrite, one new container
  validator, mutation guards, tinytest)

## Goal

Every categorical ingestion entrance - training and test, dense
slice, CSC slice, reference code - bounds codes against the
TRAINING-side K: declared (cscCategoryCounts) for CSC-trained,
inferred max+1 for dense-trained, store numCuts once a store exists.
Closes the verified mis-bin/shift-UB/heap-overread class
(typed-ingestion.md Survey); refusals match the dense error texts.

## Context

The hole is SPECIFICALLY creation-time validation plus the container
parse; predict, getTrees(newdata), setData, setPredictor value
bounds, setCutPoints (categorical refused), and setState/
installForests (buildFromFlatBelow bounds masks) are verified safe.
- validateCategoricalPredictors (R_interface_bartcore.cpp) is broken
  four ways: skips CSC-trained columns; bails on CSC-backed TEST
  columns whatever the training source (rawParsedTestColumn NULL);
  rawTrainingColumn is NULL for CSC-trained (naive skip-removal
  derefs it - scan the CSC slice instead); and it runs pre-store, so
  the bound must be reconstructed, not read off numCuts.
- resolveCscCategoricalReferences bounds the reference code against
  the CONTAINER's declared count, not the training K (test-side call
  passes categoryCountsOut null) - and refCode is what every implicit
  row reads, the higher-volume variant.
- Training side: a container whose stored codes reach its own
  declared K is accepted (numCuts = cscCategoryCount; codeFor is a
  bare cast) - the same class lands on the training store.
- Critique record (2026-08-05, refuting): BLOCKER x5, all folded in
  here; its counterfactual validators passed full tinytest 3484,
  bitwise-identical draws on 4 designs, zero false refusals across
  the valid battery (in-range declared-unobserved, larger-K and
  smaller-K containers, NA, xbart). Its baseline alarm was moot
  (equivalence-fbd2168 is historical-classic, not a gate).

## Constraints

- Refuse, never clamp; validation keyed on the VIEW being ingested
  (dense slice | CSC slice | reference code), bound from the training
  side. In-range declared-but-unobserved stays accepted.
- Residual, recorded not fixed: a level-order-permuted foreign
  container with in-range codes passes (semantic alignment is
  typed-ingestion's job); tests/cpp builds CSC stores directly
  (in-tree, ours - not a counterexample to bridge validation).
- rng neutral; no engine code change (one stale doc comment fix at
  mutateCscColumnFromDense allowed).
- Out of scope: dense declared-K asymmetry (typed-ingestion slice 1),
  factors="indicators" collapse, the mutation-shape lift.

## Steps

1. Creation: rewrite validateCategoricalPredictors - per-column
   training bound (declared K via data.cscCategoryCounts, else
   inferred max+1); validate the training CSC slice against its
   declared K; validate the test view whatever its backing (dense
   matrix, dense slice, CSC slice, reference code) against the bound.
2. Mutation: new validateTestContainerAgainstStore(const ColumnStore&,
   const ParsedTestContainer&) bounding stored values, dense-backed
   container columns, and the reference code against store numCuts;
   call it after parseTestContainer at BOTH entrances
   (setTestPredictor and setTestPredictorAndOffset).
3. Guard the corrupting mutation entrance: bridge setPredictor/
   updatePredictor entries and dbarts_sampler_setPredictor
   (C_interface.cpp) refuse on a design with any CSC-backed column,
   mirroring the R5 refusal - mutateCscColumnFromDense rebuilds the
   nonzero pattern from value != 0, silently merging the reference
   level (verified by direct .Call); the real fix is the sparse
   mutation-shape lift. Fix its stale "ordinal only" comment.
4. tinytest: out-of-range refused at creation and both container
   entrances x {stored codes, reference code} x {CSC-trained,
   dense-trained} x {inline, pooled}; training-side codes >= declared
   K; lock-ins that predict/getTrees still refuse; valid battery
   stays accepted (expect_error targets the sampler constructor, not
   dbartsData).
5. Acceptance: rerun the probe and critique scripts (scratchpad
   typed-ingestion-probe/ and critic/) - every invalid case refuses,
   every valid case fits identically.

## Verification

- R CMD INSTALL .; full tinytest (3484 + new); air format --check .;
  local R CMD check (tinytest is R code).
- Equivalence trio bitwise via the dedicated harnesses (MANIFEST
  current: equivalence-7903855, bcf-equivalence-99205ee,
  multinomial-equivalence-ec2a3d0).
- tests/cpp make && ./test_bartcore unchanged-green; CI incl.
  sanitizers to green before landing more on top.

## Landing

- df79f17 (396+/37-): validateCategoricalPredictors rewritten per
  step 1; validateTestContainerAgainstStore at both container
  entrances; refuseCscCategoricalMutation at bartcore_setPredictor/
  updatePredictor and both flat C API predictor entries;
  mutateCscColumnFromDense comment corrected; 25 tinytest assertions
  (suite total 3509).
- Two deviations, both accepted at review: the mutation guard keys on
  CSC-backed CATEGORICAL columns only (the literal any-CSC predicate
  would break the documented dgCMatrix ordinal per-column mutation
  path, test-data-sparse.R); dbarts_sampler_updatePredictor is
  guarded too (identical hazard, unnamed in the plan).
- Gates run independently at review on top of the implementer's own
  battery: INSTALL --preclean; tests/cpp clean-build all pass;
  tinytest 3509/0; trio bitwise (27/27, 5 BCF scenarios x 6 channels,
  3 multinomial x 5); air clean; R CMD check Status OK; all six CI
  legs green incl. sanitizers.
- Residual as recorded in Constraints: a level-order-permuted
  in-range container is still accepted; semantic alignment belongs to
  typed-ingestion.
