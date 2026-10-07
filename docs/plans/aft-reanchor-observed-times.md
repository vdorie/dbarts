# aft-reanchor-observed-times: a survival sampler's response range is read from the observed times, always

Status: LANDED 2026-10-07, 21bbbfd7 to 11377982 (dec-B296; slice A of the four that work was cut into, see
[response-scale-rows.md](response-scale-rows.md)).

agent: one opus implementer (engine and its tests); opus reviewer.
rng: by call.
- NEUTRAL for every sampler that is not aft; for an aft sampler that never calls `setOffset` with
  `updateScale = TRUE`; and for that call on an aft sampler with no censored row, or one whose censored
  rows still sit at their censoring times, which is straight after creation only: `setResponse`, a status
  change and a `setActiveRows` reactivation each redraw censored rows without a sweep, and the call after
  any of them is in the changed class below (corrected 2026-10-07; see Calls made in planning).
- POSTERIOR-CHANGING for one sequence: an aft sampler with a censored row, `setOffset(updateScale = TRUE)`
  through the R method or the flat C entry, after a sweep or after a restore of drawn times. Today each
  chain takes its range from its own drawn times; afterwards every chain takes the range of the observed
  times less the new offset.
window: straight after [leaf-conversions.md](leaf-conversions.md), which rewrites what the same call does
to kept draws; before slice B of response-scale-rows, which makes the same routine read the rows in the
likelihood. Serial with any other work in model.hpp.
budget: ~240 lines (engine ~40, tests/cpp ~100, tinytest ~70, help, design note and comments ~30), upper
figure 450. It rests on one routine read in the leaf-conversions tree and on there being no bridge, R or
header change; plans have run 1.5 to 2 times low, and tests are where they ran over.

## Goal

After a re-derivation of the response range, every chain of a survival (aft) sampler is in one transform,
the one a sampler created on the same observed times and offset would have. A drawn time never enters a
range.

## Context

Marks: (ran) on the tip's build and again on a library built from the leaf-conversions branch after its
last engine commit (not built by the planner), the two outputs equal line for line; (read) in the source
of the leaf-conversions tree.

- (ran) 200 rows, about half censored, two chains, 20 sweeps. At creation both chains hold shift 1.3856
  and scale 5.5858, the midpoint and width of the log observed times, [-1.4073, 4.1785]. After the
  sweeps the log times in force span [-1.2354, 4.9796] in chain 1 and [-1.2354, 5.4308] in chain 2.
  `setOffset(rep(0, n), updateScale = TRUE)` then leaves chain 1 at shift 1.8721, scale 6.2150 and chain
  2 at 2.0977, 6.6662: two chains, two leaf priors, and the model's `response.range` records chain 1's.
- (ran) With an offset in force and re-derived at that offset after 20 sweeps: 1.9072, 6.6168 and 1.9020,
  6.6065, where the observed times less the offset give 1.3410, 5.9580. A copy of that sampler holds
  chain 1's transform in both chains.
- (ran) The same call straight after creation gives both chains 1.3410, 5.9580 exactly; with no censored
  row both chains agree after any number of sweeps; and `setResponse(updateScale = TRUE)` gives one
  transform, the observed times', after sweeps.
- (read) Why. [`AFTResponse::setOffset`](../../src/bartcore/model.hpp) forwards to the Gaussian response
  it holds, whose [`GaussianResponse::rescale`](../../src/bartcore/model.hpp) takes the minimum and
  maximum of the vector it is fitted to; for aft that vector holds a drawn time at each censored row.
  [`AFTResponse::setResponse`](../../src/bartcore/model.hpp) and
  [`AFTResponse::setData`](../../src/bartcore/model.hpp) write every censored row back to its censoring
  time first, so they read observed times, as creation does. `setOffset` is the one site.
- (read) leaf-conversions has [`Chain::setOffset`](../../src/bartcore/chain.hpp) read the transform before
  and after the call and restate that chain's kept draws; it is right for whatever transform the response
  installs, so this slice changes no line of it. Its tests re-anchor an aft sampler through `setResponse`
  only.
- (read) [aft-status-setter.md](../design/aft-status-setter.md) records the behaviour as a property
  ("re-anchors on latents"), not as a defect; that sentence goes.
- (ran) Released 0.9-34 has no aft family: no `family` argument and no survival code in its namespace.
- (read, by search) No scenario of the three equivalence corpora, no exact gate and no snapshot file
  re-derives an aft range. One tinytest does
  (["expect_silent(sampler$setOffset(rep(0.1, n), updateScale = TRUE))"](../../inst/tinytest/test-aft-heteroscedastic.R))
  and reads no draw afterwards. No consumer package fits aft.

## The rule

The range of an aft sampler is the minimum and maximum of the log observed time less the offset: an
event's time, a censored row's censoring time. It is read that way at creation, at `setResponse` and
`setData` (as now) and at `setOffset` with `updateScale = TRUE` (new). What the trees are fitted to is
still built from the times in force, drawn ones included: the call draws nothing and moves no drawn time.
The residual sd and its prior are restated in the new units as now.

## Constraints

- Every call named NEUTRAL above has the bits it has today. In particular `setOffset` with
  `updateScale = FALSE`, and with `TRUE` where no row is censored, runs the code it runs now.
- No bridge, R, flat-header or facade change. `--preclean` on the engine commit, for the edited header.
- The routine is left in the shape slice B extends: one place computes a range, from a vector its caller
  names, and aft's caller names the observed times. Slice B adds which rows that place reads.
- Kept draws: leaf-conversions' restatement is not edited and must still hold `predict` to 1e-12 across
  the call.

## Steps

1. Engine. [`GaussianResponse::rescale`](../../src/bartcore/model.hpp) takes the vector the range is read
   from, its own response where none is given, and builds the working response from its own response as
   now. [`GaussianResponse::setOffset`](../../src/bartcore/model.hpp) hands such a vector through.
   [`AFTResponse::setOffset`](../../src/bartcore/model.hpp), at `updateScale` true and with a censored
   row, hands the observed log times: the times in force with each censored row at its bound
   ([`censorBound_`](../../src/bartcore/model.hpp)). Its comment and the class comment say that a range
   is never read from a draw.
2. tests/cpp, beside [`testAFTStatusSetter`](../../tests/cpp/test_model.cpp), one case, two chains with
   their own generators, censored rows among them, 20 sweeps, then `setOffset(offset, true)`:
   - one transform: both chains report the same scale and shift, bit for bit, equal to the literal from
     the observed log times less the offset written in the test (fails today: two transforms, neither
     the literal);
   - nothing drawn: each chain's log times are bit for bit what they were, its generator has consumed
     nothing, and its working response is those times under the new transform to 1e-12;
   - the residual sd in response units is unchanged to 1e-12;
   - neutral: before any sweep the call gives, bit for bit, the scale, shift and working response that
     `setResponse` with `updateScale` true onto the response in force gives under the same offset, which
     reads observed times on both builds;
   - kept draws on: the chain's replay at fixed points is unchanged to 1e-12 across the call.
3. tinytest, in test-aft.R: two chains, a run, `setOffset(offset, updateScale = TRUE)` with the offset in
   force and with a new one: the per-chain leaf prior (`getLeafPrior`) is one row repeated, its
   `response.shift` and `response.scale` equal the literal, the model's `response.range` equals it, a
   `copy()` holds it, and five more sweeps are finite. Once through the flat entry
   ([`dbarts_sampler_setOffset`](../../src/C_interface.cpp)), by the existing C-API test shim.
4. Mutations (Verification).
5. Records. The `updateScale` item of [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) gains one
   sentence: for an aft sampler the range is always that of the observed times, a censored row's
   censoring time among them. [aft-status-setter.md](../design/aft-status-setter.md): the sentence that
   the transform "re-anchors on latents" is replaced by the rule, dated. No NEWS item: 0.9-34 had no aft
   family, so no released behaviour changes. The landing note here.

## Verification

- `R CMD INSTALL --preclean` into the slice's own library; tests/cpp builds and passes, clean under ASan
  and UBSan; the full tinytest suite.
- On a reference build: the four `test-reproducibility-*.R` files pass unchanged, and the three compares
  are bitwise, every scenario reporting identical draws, counted, with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against `bcf-equivalence-1b7d730c.rds`, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. Nothing is re-recorded. A scenario that is not identical is a
  finding.
- Every gate [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) lists, in `quick` mode,
  unchanged; `aft-exact.R` and `aft-hetero-pit.R` among them. No new exact arm: after the fix a
  re-derived sampler is in the transform a sampler created on the same observed times and offset holds,
  which those gates cover, and the fix is held to that identity (step 2, "one transform").
- Mutations, each to fail the named check:
  - the observed times are not handed on (today's code): step 2 "one transform", step 3's one row
    repeated;
  - the working response is built from the observed times too (drawn times lost): step 2 "nothing drawn";
  - the range's arithmetic is reordered on the way (a bit moves where nothing was drawn): step 2
    "neutral";
  - the range is read from the observed times but before the new offset is installed: step 2's literal,
    step 3's new offset.
- `lintr::lint_package()`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`, each on its own exit
  status; `R CMD check --as-cran` on a tarball from a clean copy (man/ is touched).
- Cost: one pass over the rows and one temporary vector per re-derivation on an aft sampler with a
  censored row. Not a sweep path; no bench compare.

## What a consumer sees

Nothing. stan4bart, bartCause, treatSens and bairrtt fit no aft model (searched their R and src
directories for the family; read, not run). stan4bart and treatSens do call the flat `setOffset` entry
with `updateScale` true, on gaussian and probit samplers, where this slice changes no bit.

## Out of scope, and where it goes

- Which rows a range reads (the rows in the likelihood): slice B, [response-scale-rows.md](response-scale-rows.md).
- A re-derivation that keeps the fitted function when no response changed: slice C, stated there.
- The record keeping one chain's transform when chains disagree: after this slice no call left makes them
  disagree; nothing is added to detect it.

## Calls made in planning

Agent-made, for the maintainer's later mark.

- The aft range is always read from the observed times, never from a chain's drawn times, at creation and
  at every re-derivation (the coordinator's ruling on the critique of the response-scale-rows design,
  2026-10-07). The alternative, today's behaviour, gives each chain its own prior.
- This is its own slice and lands straight after leaf-conversions, ahead of the change of rows (the same
  ruling): it is a defect of today's call whatever the rows rule is.
- The range routine takes the vector to read as an argument, and aft builds that vector at the call. Not
  taken: aft computing the minimum and maximum itself and installing them, which would be a second place
  that knows how a range is taken, to be kept in step with slice B's rows.
- No exact-gate arm and no NEWS item, for the reasons under Verification and step 5.
- The help gains a sentence although no released behaviour changes: `updateScale` on an aft sampler is
  new in 1.0-0 and the sentence says what it reads.
- Rechecked against the tip before building, 2026-10-07, after leaf-conversions had landed. Every symbol
  cited is there under its name and the routines read as Context says. A probe on a build of the tip gave
  the first measurement to the digit (1.3856 and 5.5858 at creation; 1.8721, 6.2150 and 2.0977, 6.6662
  after the call; the record chain 1's) and, on an offset of its own, the facts of the second and third:
  two transforms after sweeps, a copy holding chain 1's in both chains, one transform straight after
  creation and with no censored row.
- What the recheck moved: the `rng:` line counts "straight after `setResponse`" with the calls whose
  censored rows still sit at their censoring times. They do not: `setResponse` redraws every censored row
  above its bound before it returns (ran: all of them above, on the tip's build). So `setOffset` with
  `updateScale = TRUE` straight after a `setResponse` belongs with "after a sweep", in the changed class;
  straight after creation is the one neutral case with a censored row. The rule is as written.
- Step 2's "neutral" check is built to that. Before any sweep the call's scale and shift are, bit for bit,
  those of a response swap that re-derives the range under the same offset; its working response is that
  swap's on every row the swap does not redraw (the events), and on every row it is that of a sampler
  created at the offset, the other path that reads observed times on both builds.
- Calls made in building. [`GaussianResponse::setOffset`](../../src/bartcore/model.hpp) gains a second,
  non-virtual form that carries the vector, and the virtual forwards to it with none, so no virtual
  changes. The range is taken in the working buffer before the working response is built there, so the
  one temporary vector is aft's. The engine test censors the row whose observed time less the offset is
  the largest, so that every chain holds a drawn time above the observed range. The R test reads each
  chain's transform from the bridge's per-chain leaf-prior matrix, as test-calibration-midchain.R does
  (`getLeafPrior` reports the first chain's), and reaches the flat entry through the shared consumer in
  test-aft.R, skipped where there is no compiler, as test-monotone.R does. No landing note is written
  here before the push.
- Found by the mutations: from R the fourth one (the range read before the new offset is installed) is
  seen only on a sampler made without an offset. The bridge and the flat entry copy a new offset over
  the one in force before the engine is reached, so by then the old offset is gone and "a new offset"
  cannot tell the two orders apart. Step 3 therefore has a third case, an offset where there was none,
  and the flat entry's check uses it. The third mutation moves a bit only where the reordered sum does
  not round back; adding and subtracting 1 left the engine fixture's two extremes as they were (R caught
  it once, through the recorded range), 3.3 did not.
- After the review, 2026-10-07: both tests gained a sampler with one censored row, the one whose observed
  time less the offset is the largest. Until then nothing failed when the last censored row was left at
  its drawn time or when a sampler with a single censored row took the old path (the reviewer's two
  surviving mutations). The `rng:` line is corrected in place.

## Landing note

Landed 2026-10-07 in 21bbbfd7 to 11377982, one push with step 5 of
[leaf-conversions.md](leaf-conversions.md). The independent review, told to refute, found nothing
blocking. It reproduced the defect on the build before (two chains at (2.003, 5.536) and (2.032, 5.594)
where the observed times give (1.383, 4.401)) and found 23 of 23 re-derivation routes on the observed
times' range afterwards, against 7 of 23 before, with the drawn times unmoved and kept draws holding to
9e-16, a variance forest included. Of 33 seeded fits compared between the two builds 31 are identical,
the two that differ being the two built for the changed class; the identical ones cover every family
built on the gaussian response and every aft call this plan names neutral. Two of the review's ten
mutations survived the first build's tests, a skipped last censored row and one censored row left on the
old path; a case with a single censored row holding the largest time, swept and re-anchored, was added
and fails under both. The `rng:` line's "straight after `setResponse`" was wrong and is corrected:
`setResponse` redraws every censored row, as a status change and a row brought back by the mask do, so
straight after creation is the one neutral case with a censored row. Gates: as in the landing note of
leaf-conversions.md, the same push.
