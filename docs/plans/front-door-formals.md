# front-door-formals

agent: opus (S1 the fold, the control formal and the proposal-mixture
  move); sonnet (S2 manual, NEWS, pkgdown, consumer check). Serial.
rng: neutral for both slices: the same prior objects, proposal mixture
  and control reach `dbarts()`, so the three bitwise equivalence
  baselines (gaussian, BCF, multinomial) are expected IDENTICAL and no
  snapshot moves.
window: R plus one bridge read (the proposal mixture moves from the
  model object to the control). Pre-release: names lock at the CRAN push
  (dec-B79), so this arc changes the surface once, and its retirements
  expire with the existing tombstones at 1.1-0 rather than opening a
  second round.
budget: S1 ~250 R (bart ~120, xbart ~40, tombstones ~40, class and bridge
  ~50) + ~250 tests; S2 ~300 Rd/NEWS.

Decision: dec-B116 in [docs/decisions.md](../decisions.md).

## Goal

`bart()` reads as a model specification plus the run knobs everyone
touches, the way `lmer()` and `gam()` do: priors are objects with one
promoted scalar (`k`), and every other engine or sampler setting is
reached through `control = dbartsControl(...)`. No user of `bart()` or
`xbart()` has to drop to `dbarts()` and a sampler loop to change a
setting the control carries.

## Context

- The count today: `bart()` has 51 formals, `xbart()` 31, `dbarts()`
  about 35, `dbartsControl()` 23. Nine `bart()` scalars mirror slots on
  the three prior objects the August consolidation admitted
  ([front-door.md](front-door.md), S2): `power`, `base`, `split.probs`
  on [`cgm`](../../R/model.R) and [`dart`](../../R/model.R); `k`,
  `prior.scale` on [`normal`](../../R/model.R); `sigdf`, `sigquant` on
  [`chisq`](../../R/model.R). `sigest` is a data-derived scale anchor
  with no object, and `proposal.probs` is the tree-move mixture, a
  sampler setting held on [`dbartsModel`](../../R/A_class.R) as its
  `p.birth_death` and sibling slots and read there by
  [`parseModel`](../../src/R_interface_bartcore.cpp).
- Neither `bart()` nor `xbart()` accepts a control. The August ruling
  (archive/bart2-argument-consolidation.md, fork 3.a) kept `control =`
  off because a control carries fit state: [`dbartsSpec`](../../R/spec.R)
  parks the variance forest, dispersion, survival status and forest map
  on `bartcore.*` attributes, and `setControl` copies them forward, so a
  control taken from an earlier fit would smuggle that state into a new
  one. The four engine settings of dec-B91 (`categoricalExhaustiveCap`,
  `testFitParallelCutoff`, `predictParallelCutoff`,
  `sparseDensityThreshold`) are therefore reachable from `dbarts()`
  only.
- The retirement mechanism exists: a retired name arrives through the
  dots, is read by name from the matched call by
  [`resolveConsolidatedArgs`](../../R/tombstones.R), warned once per
  session, applied and cleared; any other name in the dots is refused by
  [`refuseForeignFrontDoorArgs`](../../R/tombstones.R). The registry
  [`dbartsTombstones`](../../R/tombstones.R) carries one expiry, and
  [test-tombstones.R](../../inst/tinytest/test-tombstones.R) asserts it
  against NAMESPACE and NEWS.
- `xbart()` sweeps `k`, `power`, `base` and `n.trees` as grid axes, so
  those stay its formals; only `control` changes there. Its retired
  `control` (the `xbartControlReason` tombstone) is reversed, which
  dec-B79 permits before the CRAN push.

## Constraints

- Gates: tinytest whole suite; tests/cpp (the bridge read moves);
  equivalence trio IDENTICAL; `R CMD check --as-cran` from a clean
  tarball; `lintr::lint_package()`; pkgdown check; NEWS parses;
  `tools/check-doc-freshness.R`; `--preclean` on the bridge commit.
- Draws unchanged: every default value is preserved, only its home
  moves.
- Out of scope: folding the residual prior onto the family object
  (`gaussian(sigma = chisq(3, 0.9))`, the completion of dec-B98's rule;
  a follow-on recorded in the TODO); any engine change beyond the bridge
  read; `bartBT`, which keeps 0.9-34's 31 formals; `dbarts()`'s own
  formals beyond `proposal.probs`.

## Steps

S1, the fold, the control formal, the proposal mixture:

1. `proposal.probs` becomes a `dbartsControl` slot and formal (default
   unchanged, validity rule moves with it, refused by `setControl` after
   creation like the four engine settings unless the sampler already
   accepts a mid-run change through `setModel`, in which case keep that
   route working through `setControl` instead). The bridge reads it from
   the control expression. `dbartsModel` loses the `p.*` slots. The
   saved-state reader tolerates a state that carries the old model
   slots. `dbarts()` loses the formal; its dots carry the retired name.
2. `bart()` loses `power`, `base`, `split.probs`, `prior.scale`,
   `sigdf`, `sigquant` and `proposal.probs`. Each joins
   `consolidatedArgsFor$bart` with a reason naming its object (`cgm()`,
   `normal()`, `chisq()`, `dbartsControl()`), warned once per session,
   applied with the old semantics (the existing scalar-to-object
   resolver keeps working on the consolidated values), refused by name
   when the object is also supplied. `k` and `sigest` stay. The legacy
   door `bartBT` is untouched.
3. `bart()` and `xbart()` gain `control = dbarts::dbartsControl()`
   immediately after `keepFits`/`callback` (bart) and after `tree.prior`
   (xbart). Precedence, one rule: a flat formal the caller supplied wins
   over the control's slot; a slot the caller did not name flat is taken
   from the control; the fields a door forces (xbart's n.chains,
   keepTrees, keepTrainingFits, updateState, verbose; bart's n.samples
   handling) are forced after both. A control carrying any `bartcore.*`
   attribute is refused with a message naming a fresh `dbartsControl()`.
   The `xbartControlReason` tombstone and `refuseRetiredXbartControl`
   are deleted; the registry and its test follow.
4. Tests: the seven retirements warn once and produce identical draws to
   the object spelling; the collision refusals; `control` reaches the
   four engine settings from both doors; the flat-wins rule; the
   fit-state refusal; xbart's grid still runs; `dbartsControl(proposal.probs = )`
   validity and the bridge printout under verbose.

S2, manual, NEWS, pkgdown, consumers:

5. `man/bart.Rd` and `man/xbart.Rd`: Usage and Arguments shrink; the
   retired spellings are one paragraph under `...`; `control` is
   documented with the precedence rule and the refusal; the Arguments
   block is grouped under headings (data and formula; family and priors;
   sampling run; output; engine) in the order brms and mgcv use.
   `man/dbartsControl.Rd` gains `proposal.probs`. NEWS records the
   retirements and the control formal in the existing 1.0-0 tombstone
   list. pkgdown check.
6. Consumer check: bartCause (dbarts-1.0) forwards dots to `bart()`, so
   any of the seven names its tests pass now warn; stan4bart (bartcore)
   filters on `dbartsControl` and `dbartsSpec` formals, so
   `proposal.probs` moving onto the control changes which filter admits
   it. Run both suites against the slice library; fix and push to the
   sister branches where needed.

## Verification

    R_LIBS=<lib> Rscript -e 'tinytest::test_package("dbarts")'
    cd tests/cpp && make && ./test_bartcore
    R_LIBS=<lib> Rscript benchmarks/R/equivalence.R compare benchmarks/baselines/equivalence-f0236082.rds   # 52 identical
    (BCF and multinomial compares per MANIFEST; identical)
    R CMD check --as-cran <tarball>
    R_LIBS=<lib> Rscript -e 'length(formals(dbarts::bart))'   # 45

## Landing note, S1 (2026-09-10)

LANDED at 1a6da4d8313b6a086939ad7947161c8d320711f4, four commits (the fold and the control formal, then
one refinement and two reviewer fixes):

- babe360c3a7b9244ca8933c8ce8f78dc0f0e95ff Fold bart's prior scalars onto their objects and open both front doors to a control
- 64278e15ed7f67eb47c206808b9eed50389927cf Let a control speak for a slot it named, not only for one that differs
- 7fb7c498929f6cd55c7fbea56387d948e30ca9d3 Let a control speak for the settings the doors read as locals
- 1a6da4d8313b6a086939ad7947161c8d320711f4 Refuse the retired mixture beside a control that named the same slot

`bart` has 45 formals (was 51), `xbart` 32, `dbarts` 26, `dbartsControl`
24. The seven retired names ride [`consolidatedArgsFor`](../../R/tombstones.R)
with reasons naming their objects; `k` and `sigest` stay. The tree-move
mixture is a `dbartsControl` slot and formal read by
[`parseProposalProbs`](../../src/R_interface_bartcore.cpp) from the
control; a mid-run change still reaches the engine, through
`$setControl` re-installing the model, and a refused install rolls the
stored control back. Both doors merge a supplied control through
[`mergeFrontDoorControl`](../../R/dbarts.R): a flat name in the door's
matched call wins; otherwise the control's slot is taken when the
constructor's `dbarts.supplied` record names it or the slot differs from
a fresh control, so an explicit default and a post-construction slot
edit both speak. A control carrying a `bartcore.*` attribute is refused
by [`refuseFitStateControl`](../../R/dbarts.R). Review found and fixed:
`keepTrees` and `seed` re-read as locals after the merge; xbart's cell
samplers all seeded from a control-carried seed; the stale control
after a refused mixture install; the retired `proposal.probs` beside a
control naming the slot not refused like the other six. Gates at the
tip: tinytest 8539/0; tests/cpp all passed; equivalence 52/52, BCF
12/12, multinomial 11/11 identical; `R CMD check --as-cran` OK; air,
lintr, doc-freshness, rc-codoc, pkgdown, NEWS parse clean. S2 (manual
prose, argument grouping, NEWS, consumers) open.

## Landing note, S2 (2026-09-10)

LANDED at b1d2ba6b0c80b1fdd25b4a96049bfd3418e297f1, three commits (the manual and NEWS, then one reviewer
fix):

- ee8a2ae23e689bff87e2c8310fe1a37d8900a033 Group bart's Arguments block and finish the front-door manual prose
- 219aa0fe50b38db1f672113d36a04f7ffd88e623 Record the control formal and proposal.probs move in NEWS
- b1d2ba6b0c80b1fdd25b4a96049bfd3418e297f1 Split split.probs out of the power/base tombstone entry

The `bart` Arguments block is grouped by bold lead-ins inside
`\arguments` (data and formula; family and priors; sampling run;
output; engine), the same 51 items reordered and none dropped; the
seven retired spellings are one paragraph under `...`, each naming its
object, and `split.probs` names `cgm()` alone since `dart()` has no
such formal. `control` is documented on both doors with the precedence
rule and the fit-state refusal; `proposal.probs` is documented on
[dbartsControl](../../man/dbartsControl.Rd) with its six elements and
the frozen mixture; the stale xbart-control text in
man/dbarts-deprecated.Rd and NEWS is corrected. Consumers, each
installed against this tip and run in full: bartCause dbarts-1.0 764
tests, 0 failures, no retired name in its sources; stan4bart bartcore
542 tests, 0 failures, its formals filter admitting `proposal.probs`
through both constructors harmlessly. Gates: `R CMD check --as-cran`
OK; doc-freshness, rc-codoc, pkgdown, air, NEWS parse clean;
test-tombstones.R 139/0, test-front-door-control.R 100/0. The arc is
complete; the residual prior onto the family object stays the recorded
follow-on in the TODO.
