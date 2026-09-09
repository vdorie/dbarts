# Removing the swap tree-proposal

Status: LANDED, 2026-09-07; AMENDED 2026-09-07 (swap restored at default zero, section 9).

Amended by [pure-c-header](../plans/pure-c-header.md#pure-c-header): the flat C header creates no sampler and
no longer declares the predictor, test-data, weight, active-row, per-forest, state,
tree-extraction or augmentation entries - each is a method on the R sampler object the
handle is now read from. The `retired:` cites below name constructs that are gone; what
this record says about the R and engine sides still holds.

The swap move is removed from the MCMC kernel before 1.0 (VD, 2026-09-07, option A of the fork the mixing survey left open). This
record states the decision, its evidence and exactly what the slice did; it does not reargue the call. Sections 2 to 5 are written
in the present tense of the proposal and describe work that has since landed. Section 9 records the partial reversal: the move is
back in the kernel at a default of zero, so every symbol sections 2 to 8 call deleted exists again and their cites are live.

## 1. The decision, and its evidence

Swap is nearly all no-op: 70.14 to 76.97 percent of its proposals never reach a score in every census cell, so it accepts 1.71 to
4.30 percent of the proposals it MAKES even though per SCORED proposal it accepts 5.73 to 15.06 percent, and 23.55 in the earlier
200-tree count - the weakness is the no-op rate, not the acceptance ratio
([6.1 Stage 0 - the move census (pilot; no kill criterion)](tree-mixing-proposals.md#61-stage-0---the-move-census-pilot-no-kill-criterion)).
On the one criterion where the shipped mixtures separate at all - re-adaptation after a response swap, this package's distinguishing
use - the arm carrying change with swap at zero matches the shipped default on all four primary contrasts and every secondary one,
and the arm dropping both moves is the only one that loses ([14.6 Reading](tree-mixing-proposals.md#146-reading)). The one measured
cost is confined to a single tree: on P2's duplicate-column null at `m = 1` the no-swap arm parks 5 of 40 chains on an x3 root,
change being vetoed once a child splits on x1 and swap the only rule-rotating move - but at 50 and 200 trees the trap is gone, zero
stuck trees in either arm and x3 root shares of 0.281 against 0.283 and 0.323 against 0.323
([10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)).

## 2. What is removed

**[`swapMove`](../../src/bartcore/moves.hpp)**, deleted, 115 lines with its template header, and everything only it uses: the swappable-node
collection into `MoveScratch::nodeScratch`, the `swapIsSensible` pair of `ruleIsValid` calls and the whole-subtree
`interactionSubtreeIsValid` check behind them, the `applySwap`/`undoSwap` lambdas, the snapshot-refresh-restore path, and its two
census macro calls. Also in that file: [`MoveContext`](../../src/bartcore/moves.hpp)'s `swapProbability` field, the `swap`
enumerator of [`StepType`](../../src/bartcore/moves.hpp) (whose only branches outside the file are on birth and death), the header
comment naming three moves, and the census legend's two swap clauses.
[`metropolisJumpForTree`](../../src/bartcore/moves.hpp)'s dispatch chain loses its middle branch:

    was  bd = ctx.birthOrDeathProbability;  if (u < bd) birthOrDeathMove;  else if (u < bd + ctx.swapProbability) swapMove;
                                                                            else changeMove;
    now  bd = ctx.birthOrDeathProbability;  if (u < bd) birthOrDeathMove;   else changeMove;

**[`fillSwappable`, `rulesAreEqual`](../../src/bartcore/tree.hpp)**, both deleted, 18 lines, both with `swapMove` as their only caller;
`maskEquals` stays, the categorical change path and `Tree` itself reading it.
**Ten sites in [`SamplerOptions`, `ModelParameters`, `VarianceForest`](../../src/bartcore/chain.hpp)** - three struct fields, four
copies into a forest (creation, `setModel`, `buildSpecifiedForest`, `buildMultinomialForest`), the variance forest's copy and the
two `MoveContext` initializers - and **three fields in
[`Forest`, `ForestStructureSpec`, `MultinomialForestSpec`](../../src/bartcore/combiner.hpp)**, `Forest` being what `MoveContext` is
built from and the two specs what BCF and multinomial fits read. With [`ParsedModel`](../../src/R_interface_bartcore.cpp) that is
SEVEN initializers carrying the old literals; `ParsedModel`'s are dead on the create path, `parseModel` overwriting all four, but
the field still goes.

**Eight sites in [`parseModel`, `refuseUnsupportedAmplitudeComposition`, `buildMultinomialSampler`, `bartcore_setModel`](../../src/R_interface_bartcore.cpp)**:
`ParsedModel`'s field, the `p.swap` slot read, the sum-to-one check's third term, the creation printout's format string and its
argument, the `SamplerOptions` copy, the two-forest refusal's hard-coded mixture, the multinomial spec copy and `setModel`'s
`ModelParameters` copy. The refusal's literal is load-bearing: the two-forest path never reads a proposal probability off the model
at all - `ForestStructureSpec`'s own defaults supply every BCF fit's mixture - so it is both the gate and the record of what BCF
runs, and it moves to 0.6 / 0.4 with the struct defaults or refuses the new default outright. **Nothing in src/C_interface.cpp or
the shipped header**: retired: [`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) - the model crosses as a
`SEXP` and no flat struct carries a proposal probability, so the hash and the major/minor pair stand and nothing recompiles.

**tests/cpp: fourteen positional `MoveContext` initializers, plus five assignments.** The swap probability is the FOURTH element of
a positional aggregate initializer, so dropping the field binds a `double` to `const double* weights` - a hard compile error, not a
silent one. Twelve are in tests/cpp/test_moves.cpp, which the slice must edit though it names no swap move; two survive in
tests/cpp/test_interaction.cpp and a third goes with **[`testSwapSiblingStrand`](../../tests/cpp/test_interaction.cpp), deleted**,
81 lines, whose shape - a swap co-occurring a forbidden pair with neither swapped variable the stranded one - is unreachable without
the move, `testChangeStrandInvariant` keeping the walk gated. The five assignments are two in test_sampler.cpp, three in
test_model.cpp.

**benchmarks/R/swap-balance.R deleted**, 407 lines, the move's own exact-posterior detailed-balance gate. That takes TWO edits in
.github/workflows/exact-gates.yaml - the gate-list token and the header comment's "(birth/death, change, swap detailed balance)" -
plus benchmarks/README.md's balance-gate roster, which spells "swap-balance = swap" in prose, and strips poison 5, whose `find`
text was verbatim `swapMove`'s `swapIsSensible` block
(["m05"](../../benchmarks/R/mutation-battery.R)).

## 3. The surface

`proposal.probs` goes from three structural names to two, `birth_death` and `change`; the default becomes `birth_death 0.6, change
0.4`, `birth 0.5` unchanged - arm C of [14.2 Design](tree-mixing-proposals.md#142-design), the mixture 10.1 and 14.4 have already
measured. SIX R default VECTORS spell it: [`defaultProposalProbs`](../../R/model.R); the [`dbarts`](../../R/dbarts.R), `bart2`
(R/bart.R) and [`dbartsSpec`](../../R/spec.R) formals; and two literals in the monotone branch of
[`resolveSamplerSpec`](../../R/spec.R) - the comparison default, which does NOT read `defaultProposalProbs`, and the
birth/death-only rewrite, not a mixture and landing as `c(birth_death = 1, change = 0, birth = 0.5)`. FOUR further `"swap"` literals
sit in the [`dbartsModel`](../../R/A_class.R) initializer in R/model.R: the three-name subset, the `names(probs) <-` assignment in
the one-NA branch, the all-NA fallback subset and the slot write. The class loses `p.swap` from slots, prototype and validity's sum.

**A, refuse a stale `swap` by name - and the check runs FIRST.** One helper, at the head of
`resolveSamplerSpec` - ahead of the monotone branch and the BCF `unsupported` block, both downstream inside it - the funnel
`dbarts()`, `bart2()` (through its `dbarts()` call) and `dbartsSpec()` all reach; the same helper guards the `dbartsModel`
initializer, so a direct `new("dbartsModel", proposal.probs = )` gets the identical message. Later does not work: the monotone stop
preempts it for a swap-carrying non-default, and the monotone rewrite DISCARDS a swap-carrying vector matching the new default.

| the caller's `proposal.probs` | unconstrained | with `monotone` |
|---|---|---|
| `birth_death 0.5, swap 0.1, change 0.4, birth 0.5` (the old default) | refused, naming the removal | refused, same message |
| `birth_death 0.6, swap 0, change 0.4, birth 0.5` | refused | refused |
| `swap` alone | refused | refused |
| `birth_death 0.6, change 0.4, birth 0.5` (the new default) | proceeds; sum enforced at `setValidity` | rewritten birth/death-only, silently, as today |
| `birth_death 0.7` alone | proceeds; `change` filled with 0.3 | monotone stop: not the default vector |

The cost is row two: a vector semantically identical to the new default is an error, the price of one unambiguous message for every
spelling naming a move the kernel does not have. **B, accept and ignore a zero**: `swap = 0` is silent and right, but the old
documented default reaches `setValidity` as 0.5 + 0.4 and fails with "rule proposal probabilities must sum to 1", naming nothing
about the removal, and a vector naming ONLY `swap` presents two NAs, takes the all-NA branch and silently gets the default. **C,
fold its mass into birth_death**: runs a kernel the caller did not ask for and blinds the sum check to a typo. **RECOMMEND A.** The
sum-to-one tolerance path then drops a term on both sides - R validity at `sqrt(.Machine$double.eps)`, the bridge's
`sumToOneTolerance`, over `(p.birth_death, p.change)` - and the one-NA fill sits at two names: one NA takes the residual, naming
`birth_death` alone resolving `change`, and all-NA falls back to the default. Under A the removed name never reaches the fill.

**Editing the tests is not "drop the `swap` element".** An element may be dropped only where the remainder still sums to one;
otherwise the vector is RE-VALUED, and one that is no longer a reachable kernel is DELETED. FIVE tinytest files spell the mixture.
test-sum-to-one-tolerance.R drops cleanly, its `makeModel` naming both survivors, so the 1e-9-accept / 1e-7-refuse pair holds.
test-monotone.R drops cleanly at its non-default vector, which still trips the monotone stop on `change`; it separately loses two
`p.swap` reads and MOVES its default assertion from 0.5 to 0.6. The monotone check precedes the model funnel, so it is only
because the helper sits ahead of BOTH that a swap-carrying vector under `monotone` meets the removal's message rather than the
monotone stop. test-spec.R loses a `p.swap` read. test-bcf-creation.R must be
RE-VALUED: trimmed, its vector sums to 0.9 and fails `setValidity` with a message its `expect_error` pattern does not match, so it
needs a non-default summing to one. test-argument-surface.R pins the OLD default in a defaulted-versus-explicit equivalence check
and is re-valued. In benchmarks both P2 scripts' `default` arm is DELETED - no longer a reachable kernel - while `birthdeath`,
`noswap` and bd-balance.R's arm drop their element cleanly.

**Four Rd files**: the usage lines of man/dbarts.Rd, man/bart2.Rd and man/dbartsSpec.Rd, and the `proposal.probs` argument text of
those plus man/bart.Rd (`bart()` itself takes `NULL`, so the Rd alone). **NEWS**: a new item in inst/NEWS.Rd's 1.0-0 UPGRADING
block, and a CORRECTION to the existing 1.0-0 entry spelling the old vector verbatim (section 6). **Consumer exposure: none by
name** - `git -C ../stan4bart grep -n proposal bartcore -- R` and `git -C ../bartCause grep -n proposal dbarts-1.0 -- R` are both
empty, stan4bart forwarding `bart_args` through `formals(dbarts::dbartsSpec)` generically and bairrtt's `proposal` hits being its
own theta step. An unnamed caller inherits the new default; a stan4bart user whose `bart_args` carries `swap` is refused.

## 4. RNG, and the one re-record

**Bitwise neutrality against the current baselines is impossible.** The default moves from 0.5 / 0.1 / 0.4 to 0.6 / 0.4, so a tree
whose move-type uniform lands in [0.5, 0.6) takes birth/death where it took swap; at 50 trees or more that happens inside the first
sweep with probability at least `1 - 0.9^50`, above 0.994, and the stream desynchronizes from there. **But the OLD kernel at the NEW
mixture is bitwise the new kernel.** With `swapProbability = 0.0` the middle test is `u < bd + 0.0`, which is `u < bd` exactly in
IEEE for any finite `bd` and is already known false, so control reaches `changeMove` at the same stream position and `swapMove`
draws nothing. The selection consumes one uniform per tree per sweep either way, and nothing anywhere reads the mixture for a
proposal count or a veto budget. That identity is the oracle section 5 rests on.

**Every RNG-locked artifact in the repo, and its obligation.**

1. `benchmarks/baselines/equivalence-1e5f80b2.rds` (50 scenarios) - re-recorded, a new file plus a MANIFEST row marked `current`
   with its P17 oracle, 1e5f80b2 demoted; the .github/workflows/equivalence.yaml pin and mutation-battery.R's path both move.
2. `bcf-equivalence-3c81d6df.rds` (12) - same, plus the `--cross-host` pin in .github/workflows/exact-gates.yaml.
3. `multinomial-equivalence-4d9a3337.rds` (11) - same, plus that workflow's second pin.
4. feature-matrix.md's Evidence paragraph, naming all three with their counts; check-doc-freshness recomputes those from
   [`makeScenarios`](../../benchmarks/R/equivalence.R) and its two siblings, and refuses a cited file that is gone.
5. The gate ledger, docs/plans/review-2026-08-24/gate-ledger.md: its [f39] baseline footnote, and section 4's counts, 20 scripts
   falling to 19 and 15 single-tree gates of 20 to 14 of 19.
6. The four seeded-drift tripwires, the only tinytest files pinning exact draws:
   test-reproducibility-continuousResponse-singleThreaded.R, -multithreaded.R, -binaryResponse.R and -xbart.R, regenerated by
   ["Regenerates the reference values"](../../tools/regenerate-snapshots.R) - it replays each file top to bottom and rewrites every
   `referenceX` list, values depending on the file's full execution history rather than the preceding seed. Every other numeric
   literal in inst/tinytest is an analytic constant or a tolerance band.
7. `bench-sampler-ab1dc52.csv` - a COMPARE run, not a re-record. Birth/death takes swap's 0.1 and is a scored move where 70 to 77
   percent of swaps were no-ops, so a small slowdown is expected; only an arm past the harness's 1.05 ratio makes this a re-record
   with its own MANIFEST row.

**SBC pins no ranks.** .github/workflows/sbc.yaml fixes seeds and records verdicts in docs/plans, and says outright that a
draw-shifting commit reshuffles the same fixed-seed stream - which is why the gate is non-blocking. Nothing to re-record; run the
matrix once after the landing and record a verdict only where an arm flags. MANIFEST's P17 rule wants an oracle for a draw-changing
re-record; here it is a bitwise identity, not an adjudication. **P5, P6 and C1 are NOT re-run**: their pilot numbers were measured
at the former default and stay as recorded, section 10 of benchmark-surfaces.md saying so once at the top of its records while the
per-cell sentences stop calling that mixture "shipped" (section 6). Re-running three cells to re-attribute numbers no verdict rests
on is not worth the compute.

## 5. Gates

- **tests/cpp**: `make && ./test_bartcore` exits 0 once the fourteen `MoveContext` initializers are shortened and
  `testSwapSiblingStrand` is deleted. **tinytest**: 0 fail, after the five mixture-spelling files are edited and the four tripwires
  regenerated.
- **The bitwise oracle, and the P17 row.** Build a THROWAWAY private library at the parent commit with only the mixture DEFAULTS
  moved: the FIVE default R vectors to `birth_death 0.6, swap 0, change 0.4`, and the seven C++ initializers and
  `refuseUnsupportedAmplitudeComposition`'s literal to the same, leaving `swapMove` and the dispatch standing. The monotone
  birth/death-only REWRITE is not a default and must stay at `birth_death = 1`, or every monotone fit gets a change move the
  constrained leaf cannot score. Record all three equivalence baselines there, then `compare` from the landed tip: every scenario
  must report "identical draws (same RNG stream)", so the removal provably moved no draw. The baselines that SHIP are recorded at
  the landed tip, and a second `--preclean` install must replay them 50/50, 12/12 and 11/11 bitwise, gaussian `--strict-coverage`.
- **P2, an end-to-end path check, not independent evidence.** ["noswap"](../../benchmarks/R/surfaces/P2-confounded-step.R) is
  already the arm; at the new kernel it IS the default arm. Run the script at the parent commit and at the landed tip on the same
  five seeds: 10.1's no-swap row - pooled 0.457, 70.8 switches per chain, minimum 0, between-chain sd 0.149, 5 of 40 chains parked -
  must come back exactly, and ["noswap"](../../benchmarks/R/surfaces/P2-null-at-scale.R)'s rows at `m = 50` and `m = 200` must stay
  at zero stuck trees. By section 4's identity this can only fail if the identity fails, which the three oracles test over 73
  scenarios; its value is exercising the whole R-to-engine path. Deleting the `default` arm costs 10.1's default rows their
  reproducing script, which the record must say.
- **Speed**: `bench-sampler.R compare` on a quiet machine, no arm past 1.05. **The move census, re-run at the new default** on a
  `-DBARTCORE_MOVE_CENSUS` private build: [`moveTable`](../../benchmarks/R/move-census.R)'s pooled and per-move acceptance for
  birth, death and change over the four cells, recorded in the 6.1 addendum beside the mixture that produced it.
  [5. Benefit, pre-registered](perturb-move.md#5-benefit-pre-registered) recomposes those rates arithmetically from the
  swap-carrying census; this replaces the arithmetic with a measurement.

**Poisons.** For the bitwise oracle: leave the throwaway build's swap probability at 0.1 and the compare must FAIL - every scenario
falling back to a statistical comparison rather than an identical stream - showing the oracle tests the mixture, not one binary
re-run. For the refusal: a call passing the old documented default must error naming `swap`, under `monotone` as well as without it,
pinned by a tinytest. The other exact gates take no poison - bd-balance.R and change-balance.R are untouched beyond the former's arm
losing its `swap = 0`.

## 6. Docs updated in the same landing

- **The six live `swapMove` symbol cites fail the guard the moment it goes** - three in tree-mixing-proposals.md, one in
  empty-leaf-veto.md, one in benchmark-surfaces.md, one here - as do this document's `fillSwappable` / `rulesAreEqual` and
  `testSwapSiblingStrand` cites. Each becomes a `retired:` marker before the link, with prose saying the symbol is gone, or is
  rewritten per the cite grammar in tools/check-doc-freshness.R's header.
- **tree-mixing-proposals.md section 2** is the largest stale block in docs/design: the three-arm dispatch fence, the shipped-mixture
  sentence with its three cites, the `StepType` enumerator list and the Swap bullet. Also its swap-balance line-count sentence, its
  three "arm A is the shipped default" descriptions, and 14.6, whose "This does not move the shipped default" is restated to record
  the decision and point here. 6.1's two addenda keep their swap rows as measurements of a removed move and say so.
- **benchmark-surfaces.md**: one statement at the top of section 10 that its records were measured at the former default, after
  which the per-cell "shipped default"/"shipped mixture" phrasings in the 10 preamble, 10.1's table header and arm-contrast and
  no-swap paragraphs, 10.2, 10.3, 10.4 and 10.5's second, third, fifth, sixth and last bullets stop claiming a current kernel; 6.3's
  three-mixture note and 6.5's likewise. **docs/architecture.md's "Tree moves"**: the "birth/death, change, and swap" sentence and
  the `StepType` list. **monotone.md**: its unconstrained-default mix and the default vector it quotes.
- **feature-matrix.md** and **the review tour** state only the current kernel - the Evidence paragraph's baseline names, and for the
  tour the stamp alone, its tree-moves section handing over a reading order and deferring the move set to docs/architecture.md.
  **perturb-move.md's** [8. Slices](perturb-move.md#8-slices) slice 0, which says the decision is not yet written down in
  docs/design, and its section 3 sentence describing this removal's inverse.
- **inst/NEWS.Rd**: the new 1.0-0 UPGRADING item, AND the existing 1.0-0 entry that spells
  `c(birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5)` verbatim as `bart2`'s newly live default - rewritten to the vector
  that ships. **benchmarks/README.md**'s balance-gate roster and **benchmarks/R/surfaces/README.md**'s three-arm description.
  **docs/design/INDEX.md**'s row for this document, whose status phrase check-doc-freshness pairs against the `Status:` line above.

## 7. The slice, and where it lands

**One slice, not staged.** The kernel and its surface are one commit: a build with the R names removed and `swapMove` standing is a
build nothing can reach the move from, and one with the move removed and the names standing does not compile. Deletions: moves.hpp
-119, tree.hpp -18, chain.hpp -10, combiner.hpp -3, R_interface_bartcore.cpp -12, tests/cpp -86, swap-balance.R -407,
mutation-battery.R -20, two edits in exact-gates.yaml. Edits: R/ about 16 lines net, man/ 8, five tinytest files about 25, fourteen
tests/cpp initializers, the benchmark arm lists and rosters 10, plus section 6's documents. Additions: inst/NEWS.Rd +8 and a
tinytest pinning the refusal, +12. Roughly 720 lines removed and 60 added across about 35 files, plus four regenerated snapshots and
three re-recorded baselines. **It lands BEFORE the perturb kernel and before any SBC or benefit arm**, being
[8. Slices](perturb-move.md#8-slices)' slice 0 and taking the one bundled re-record, so everything downstream records against the
new kernel once: perturb's default-weight-0 neutrality claim in
[6. RNG and baselines](perturb-move.md#6-rng-and-baselines) is a claim about THIS dispatch chain, and Stage 2's control arm is this
default.

## 8. Landing

Three commits on 2026-09-07: fbff1989 the kernel, the surface and the tests; ff1d18ee the bundled re-record; the third this
document and the rest of section 6.

**The oracle held.** A throwaway library at the parent commit 39692087, with only the mixture defaults moved and `swapMove` and
its dispatch left standing, recorded all three equivalence baselines; `compare` from fbff1989 reported "identical draws (same RNG
stream)" on every scenario of all three - gaussian 50 compared / 0 skipped with zero `max |z|` lines, BCF 12 compared / 0 skipped
every channel bitwise, multinomial 11 compared / 0 skipped every channel bitwise. So the deletion moved no draw the mixture change
did not. The poison fired as designed: the same compare against the outgoing baselines, recorded at swap 0.1, fell back to the
statistical mode on every scenario of all three with zero identical streams. The gaussian partition against 1e5f80b2 is max
|z| = 3.73 over 3687 summaries, 17 with |z| > 3 and none at |z| > 4.

**The shipped baselines** are `equivalence-fbff1989.rds`, `bcf-equivalence-fbff1989.rds` and
`multinomial-equivalence-fbff1989.rds`, recorded at fbff1989 and replayed 50/50, 12/12 and 11/11 from a second `--preclean`
install. The BCF and multinomial harnesses replicate no seed, so their draws-axis fallback against the outgoing baselines reports
large |z| (up to 20.82 and 9.31) that is not a calibrated posterior comparison; their MANIFEST rows say so and rest on the bitwise
identity plus the four exact-posterior gates, all of which PASS at this tip in quick mode - bcf-exact E[mu] gap 0.0005,
bcf-exact-weak E[tau] 0.0012, bcf-exact-restricted E[mu] 0.0007, multinomial-exact all arms.

**Gates.** tests/cpp 277 ok, all tests passed, and the same under `-fsanitize=address,undefined`; the count is 278 less
`testSwapSiblingStrand`. tinytest 7449 tests, 0 failures. `air format --check`, `lintr::lint_package()`,
`tools/check-doc-freshness.R` and `tools/check-rc-codoc.R` all clean. The refusal poison: deleting
retired: [`refuseRemovedProposalNames`](../../R/model.R) - the helper is itself gone, deleted by section 9's reversal - failed
exactly the tinytest that pinned it, and restoring it passed.
`DBARTS_C_API_HASH` is unchanged.

**Two tests were pinned on the old stream rather than on semantics** and were rewritten to their intent rather than re-valued:
tests/cpp's `installForests` precondition swept a fixed two sweeps to move the twin's tree off the donor's and now sweeps until it
moves, and test-weighted-binary-ppd.R required every `w = 5` column to draw an intermediate count where one column's posterior mean
probability is 0.0005, so it now reads only the columns whose posterior mean probability is interior. Both keep the property they
were written for.

**The `default` arm of both P2 scripts is deleted**, being a mixture the kernel can no longer run. That costs
[10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)'s `default` rows their
reproducing script: those numbers stand as a record of the former kernel and cannot be re-derived from the tree. The `noswap` and
`birthdeath` arms still run, and `noswap` is now the shipped mixture.

**The P2 path check reproduces.** `P2-confounded-step.R` at the landed tip returns the duplicate-column null's no-swap row
verbatim - pooled 0.457 (0.359-0.555), 70.8 (57.9-81.2) switches per chain, minimum 0, between-chain sd 0.149 (0.036-0.251), 5 of
40 chains parked - and its birth/death rows likewise, on both designs. `P2-null-at-scale.R` returns 0.283 (0.268-0.305) at
`m = 50` and 0.323 (0.317-0.330) at `m = 200`, between-chain sd 0.009 (0.007-0.012) and 0.003 (0.002-0.004), ZERO stuck-on-x3
trees at either count. Every figure is 10.1's, to the digit. So the whole R-to-engine path carries the identity the three oracles
test at the harness level.

**Not done here, and owed.** Three things. The SPEED compare was not run at the landing: `bench-sampler.R` needs a quiet machine
and the landing box was not one (1-minute load above 3, a virtual machine resident). It RAN 2026-09-08 on a quiet-machine grant at
the level-auto tip, as a same-machine A/B against a rebuild of `bench-sampler-ab1dc52.csv`'s commit (6 alternating rounds, per-round
minima): no arm past 1.05, the birth/death-heavy run arms 0.985 to 0.996, so swap's 0.1 moving to birth/death cost nothing
measurable; the one real cost, embedded-offset at 1.032, is on the offset-mutate loop and is recorded with the re-recorded baseline
`bench-sampler-127f04ee.csv` in the MANIFEST. The MOVE CENSUS has not been re-run at the new default: section 5 asks for
[`moveTable`](../../benchmarks/R/move-census.R)'s per-move acceptance on a `-DBARTCORE_MOVE_CENSUS` build, recorded in 6.1's
addendum, and [5. Benefit, pre-registered](perturb-move.md#5-benefit-pre-registered) still recomposes those rates arithmetically
from the swap-carrying census. The SBC matrix has not been run since the landing; section 4 says it pins no ranks, so nothing is
owed there beyond one run of the matrix and a verdict only where an arm flags.

## 9. Reversal: the move returns at default zero

The removal is PARTIALLY reversed (VD, 2026-09-07). The swap move is back in the kernel and `swap` is again a legal name in
`proposal.probs`; the SHIPPED DEFAULT is unchanged, `birth_death 0.6, swap 0, change 0.4, birth 0.5`. Sections 1 to 8 stand as the
record of the removal and of the evidence for zeroing the move at production forest sizes, which is not revisited.

**The evidence is a one-tree exact-posterior gate.** `benchmarks/R/hazard-exact.R` fits ONE tree on two live columns and compares
the sampler's subject-level hazard against a brute-force enumeration over the 62 reachable trees. At the two-move kernel it FAILS,
max hazard gap 0.0120 against a full-mode tolerance of 0.004; at the three-move kernel the same gate reads 0.0008. Across the full
sweep that is 20 of 21 gates passing at the two-move kernel against 21 of 21 at the three-move one. The mechanism is specific to a
single tree: [`changeMove`](../../src/bartcore/moves.hpp) redraws a node's split variable and then its cut from the
descendant-valid set, so once the splits BENEATH the root depend on the root's own variable no new variable has a valid cut there
and the proposal no-ops; [`swapMove`](../../src/bartcore/moves.hpp) is the only move that rotates a child's rule UP, which is how a
one-tree chain crosses between rootings. Section 1's own measured cost - 5 of 40 chains parked on an x3 root at `m = 1`
([10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)) - is the same failure seen
through a different statistic. At 50 and 200 trees the ensemble self-averages it away, so nothing about the shipped default changes.

**Bitwise neutrality holds in the other direction.** Section 4's identity is symmetric: with `swapProbability = 0.0` the dispatch
test `u < bd + 0.0` is `u < bd` exactly, so restoring the branch moves no draw at the shipped default. The three canonical
baselines recorded at fbff1989 replay 50/50, 12/12 and 11/11 bitwise from the restored build, gaussian with
`--strict-coverage` and zero `max |z|` lines, and the four seeded-drift tripwires pass UNCHANGED. Nothing is re-recorded.

**What the slice did.** Restored from 39692087: `swapMove` with its helpers and dispatch branch, `MoveContext::swapProbability`,
`StepType::swap` and the two census macro calls; `fillSwappable` and `rulesAreEqual`; the swap-probability field on
[`SamplerOptions`, `ModelParameters`, `VarianceForest`](../../src/bartcore/chain.hpp) and on
[`Forest`, `ForestStructureSpec`, `MultinomialForestSpec`](../../src/bartcore/combiner.hpp); the bridge's parsed field, slot read,
three-term sum check, creation printout and four copies, with
[`refuseUnsupportedAmplitudeComposition`](../../src/R_interface_bartcore.cpp)'s literal at 0.6 / 0.0 / 0.4; the fourteen positional
`MoveContext` initializers, the five assignments and
[`testSwapSiblingStrand`](../../tests/cpp/test_interaction.cpp); `benchmarks/R/swap-balance.R` with its exact-gates token and header
clause and benchmarks/README.md's roster line; and poison 5
(["m05"](../../benchmarks/R/mutation-battery.R)). The surface keeps the new numbers everywhere: six R default vectors and the
monotone comparison default at `birth_death 0.6, swap 0, change 0.4, birth 0.5`, the monotone rewrite at `birth_death 1, swap 0,
change 0, birth 0.5`, [`dbartsModel`](../../R/A_class.R)'s `p.swap` slot back with prototype 0 and a three-term validity sum, and
the one-NA fill and both sum-to-one tolerance paths over three names. The fill resolves two unnamed elements as well when one of
them is swap, which takes its zero so the other can take the residual - `c(birth_death = 0.7)` is 0.7 / 0 / 0.3 - while `swap`
named alone leaves the birth/death-versus-change split undetermined and is an error. `refuseRemovedProposalNames` is deleted, and
["a caller-supplied three-move mixture"](../../inst/tinytest/test-proposal-probs.R) replaces the tinytest that pinned it: the
three-name surface, a one-tree fit created and run at `swap = 0.1`, the printout, and the round trip.

**The one-tree gates now ask for the move.** Every exact-posterior gate that fits a single tree and whose creation path accepts a
caller mixture passes `birth_death 0.5, swap 0.1, change 0.4, birth 0.5` explicitly, with one header sentence saying why: aft-,
categorical-, linear-, hazard-, hurdle-, t-, negbin-, ordinal- and heteroscedastic-exact, logistic-reference, and
multinomial-exact's one-tree arm. Four cannot and stay at the shipped default: the BCF two-forest path reads
[`ForestStructureSpec`](../../src/bartcore/combiner.hpp)'s own struct defaults and refuses a non-default mixture (bcf-exact,
bcf-exact-weak, bcf-exact-restricted, bcf-latent-exact), and monotone-reference's constraint rewrites the mixture birth/death-only.
Those five remain the standing one-tree exposure. bd-balance and change-balance set their own mixtures and are untouched; the P2
surface scripts regain the swap-carrying arm, named `swap`, so
[10.1 P2, the confounded step function](benchmark-surfaces.md#101-p2-the-confounded-step-function)'s rows reproduce again.
