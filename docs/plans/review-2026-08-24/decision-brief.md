# Three maintainer decisions - substance

Read-only against bartcore tip 0045507c. Every number is from a tool run in this pass:
file:line from the working tree, probe output from a staged `git archive HEAD` built with
`R CMD INSTALL --preclean` into a private library. Probes are in `scratchpad/brief/`.

## DECISION 1 - three places where one modelling rule has two implementations

### (i) The amplitude-forest combiner routes on basis shape, not on the prior it carries

WHAT AN AMPLITUDE SPEC IS. A multi-forest fit writes its mean as a weighted sum of forests:
`sum_f dot(a_f, B_f(i,.)) * f_f(x_i)`, where `B_f` is forest f's basis (an n x q_f matrix)
and `a_f` is that forest's amplitude vector - the glue scaling the forest's contribution.
This is the BCF shape: forest 0 prognostic, implicit all-ones basis, one scalar amplitude
`a`; forest 1 the treatment forest, two-column 0/1 indicator pair, amplitude pair
`(b0, b1)`. Each block's prior has two spellings: a fixed-variance normal
(`ForestAmplitudePrior.variance`, `[[combiner.hpp:440@0045507c]]`) or a half-Cauchy scale mixture, in
which the variance is a live inverse-gamma auxiliary redrawn after the block
(`halfCauchyScale`, `[[combiner.hpp:441@0045507c]]`).

THE DIVERGENCE. `drawGlue` (`[[combiner.hpp:986@0045507c]]`) picks between two implementations of the
same conditional: `drawShippedGlue` (`[[combiner.hpp:1142@0045507c]]`) when `forests.size() == 2 &&
shippedShape()`, otherwise the general q-variate sweep `drawAmplitudes` (`[[combiner.hpp:1205@0045507c]]`).
`shippedShape()` (`[[combiner.hpp:1474@0045507c]]`) tests only basis widths and canonicality - it never reads the
prior. But the two bodies do not agree about the prior: the shipped one refreshes
`prior[0].variance` alone (`[[combiner.hpp:1171@0045507c]]`), the general one loops every forest with a positive
scale (`[[combiner.hpp:1213@0045507c]]`). So a K=2 spec with canonical bases and a positive half-Cauchy scale on
forest 1 gets a FIXED-variance prior on `(b0, b1)` where the general path would sample it.
Those are two different models, not two roundings of one.

WHAT "TWO ADMISSIBLE SPECS" MEANS CONCRETELY - and the correction the review needs. There
is no `bart2`/`dbarts(forests = ...)` call that reaches it. `forestParams`
(`[[R/model.R:1157@0045507c]]`) emits the half-Cauchy scale as `if (withBasis) 0 else declared(spec$sd,
default)`: any forest carrying a basis gets a hard 0, and `resolveForests`
(`[[R/model.R:1061-1067@0045507c]]`) refuses a basis-free forest past the first. The internal creator
agrees - `bartcoreBCFSampler` (`[[R/bartcore.R:792@0045507c]]`) hardcodes 0 in the treatment forest's
slot. PROBE (`brief/probe-bcf.R`): a stock two-forest fit reports forest 1 params
`50 0.25 3 1 1 1 2 1` and forest 2 `50 0.25 3 1 0.674 0.5 0 1` - the scale is 2 on the
basis-free forest and 0 on the basis forest, exactly as the source says. The divergent
region is reachable only by hand-setting the undocumented `bartcore.forests` control
attribute and calling `new("dbartsSampler", ...)` (a spelling `[[dbartsSpec.Rd:83@0045507c]]` shows).
PROBE (`brief/probe-bcf2.R`) does that: with a declared scale of 2 on forest 2, its
amplitude prior variance reads 0.5 before and 0.5 after 40 sweeps, while forest 1's moves
1 -> 2.658. The declared scale-mixture is silently inert. Note the side effect: setting
that slot also flips `forest.ridge` (`[[R_interface_bartcore.cpp:2180@0045507c]]`), which consumes a
GIG draw per sweep - forest 1's variance differs between the two arms (1.859 vs 2.658)
for that reason, not because of the routing.
No C++ fixture reaches it either: `[[tests/cpp/test_sampler.cpp:6449@0045507c]]` puts a scale on
forest 0 only, at K = 3, where `shippedShape()` is already false.

CANDIDATE RULINGS.
(a) Make the predicate honest: add "and no forest past 0 declares a positive scale" to
    `shippedShape()`. ~3 engine lines. Both admissible spellings then get the general
    model. Draw-neutral on every R-reachable spec, so all three canonical baselines stay
    bitwise; nothing to re-record.
(b) Fix the body instead: loop f >= 1 in `drawShippedGlue` and refresh each scale-mixture
    variance after its block. ~12 lines. Same reachability, so also bitwise; but it
    duplicates the general path's refresh a second time, which is what caused this.
(c) Refuse the spec at the bridge: a basis forest may not declare a positive amplitude
    prior scale. ~4 lines in `applyAmplitudeSpec`, draw-neutral; shuts the door rather
    than deciding what is behind it.
(d) Delete the branch: make `generalAmplitudeDraw` the default. The comment at
    `[[combiner.hpp:966-975@0045507c]]` records the measurement: the general path CANNOT reproduce bcf
    bitwise (the two moment accumulators fuse differently; 21 accumulation variants tried).
    Re-records `bcf-equivalence-6e3b9fb8.rds` (all 12 scenarios) and the `bart2twoforest`
    scenario of `equivalence-5a3bc276.rds`. `multinomial-equivalence-4d9a3337.rds` is
    untouched (the multinomial combiner has no amplitudes).
GATE COVERAGE: `bcf-exact.R` enumerates the two-forest posterior with a free `a` under
`Cauchy(0, aPriorScale)` and free `(b0, b1)` under `N(0, bPriorVariance)` - it gates the
shipped model, and would gate a corrected one only if a new arm gave forest 1 a scale
mixture. The SBC lane's BCF arm covers the identified effect and prognostic functions,
not the amplitude prior spelling.

RECOMMENDATION: (a). It repairs the actual defect - a selector that reads half the state
it selects on - at three lines, is bitwise on all three baselines, and leaves the
performance branch that (d) would delete standing until a bcf re-record is wanted for its
own reasons. Priority: low. Nothing on the documented surface can express the offending
spec, so this is a latent-consistency fix, not a live-model bug.

### (ii) Variance-leaf positivity is enforced on three of four state paths

WHAT THE RULE IS. Under a heteroscedastic fit the variance forest models `s^2(x)` as a
PRODUCT of per-tree leaf factors. `rebuildVarianceForest` (`[[chain.hpp:4331-4368@658869ac]]`) scatters
each leaf factor `h` into `combinedVariance[i] *= h`, and the sweep then forms
`divisor[i] = combinedVariance[i] / factor[i]` and `treeResidual[i] = meanResidual[i] /
sqrt(divisor[i])` (`[[chain.hpp:463-464@658869ac]]`) and weights the mean side by `1/s^2(x)`. A
non-positive leaf is therefore not a bad value, it is a broken model.

THE FOUR PATHS AND THE GAP. Construction and sampling draw leaves from an inverse-chi-
square, so they are positive by construction and carry no check. Two validation sites do
check: `[[chain.hpp:3301@658869ac]]` (setState's live variance trees) and `[[chain.hpp:3334@658869ac]]` (setState's saved
buffer). The warm-start install has two arms and only one checks: the slot-sourced arm at
`[[sampler.hpp:919-923@658869ac]]` applies the law with a comment naming exactly why ("the buffer is
hand-buildable and a rebuild scatters the leaf straight into a divisor"); the live-sourced
arm at `[[sampler.hpp:891@658869ac]]` is a bare `dst.varianceTrees = src.varianceTrees;` with nothing.
Neither donor parser checks either (`[[R_interface_bartcore.cpp:7237-7251@658869ac]]`).

WHAT A USER CAN DO. The review recorded this as "reachable from a .Call with a hand-built
state, not from R's own state objects". That is too generous. `warmStartState`
(`[[R/dbarts.R:750-753@658869ac]]`) accepts a raw `bartcoreState` and returns it unchanged, so
`sampler$installTrees(state)` takes a hand-edited state directly. PROBE
(`brief/probe-var3.R`, `probe-var4.R`), on a heteroscedastic sampler with `keepTrees =
FALSE` so the donor has no saved slots and the install takes the live arm: read
`storeState()`, rewrite one leaf's 8-byte payload to -2.5, then
  - `setState(poisoned)`  -> REFUSED, "state is not consistent with this sampler"
  - `installTrees(poisoned)` -> ACCEPTED, no error
and the subsequent run, at the same seed as the clean install, produces fitted values
ranging to 24039.18 where the clean install ranges to 2.117. Silent numerical corruption,
not a crash: the negative factor makes `combinedVariance` negative on that leaf's rows, so
the mean-side precision goes negative.

CANDIDATE RULINGS.
(a) Apply the same loop at `[[sampler.hpp:891@658869ac]]` that `[[sampler.hpp:919-923@658869ac]]` already carries. ~5 engine
    lines, and the existing `WarmStartResult::varianceMismatch` code is already wired to a
    named R error. All four paths then answer the same way.
(b) Check in the donor parser (`R_interface_bartcore.cpp` around [[sampler.hpp:7237@658869ac]]) so both arms and
    any future reader inherit it. ~8 lines, better message, but puts a model law in the
    bridge rather than the engine.
(c) Leave it and document that `installTrees` trusts its donor. Free, and wrong: the
    sibling arm already decided this and `setState` refuses the identical object.
COST TO BASELINES: none under any option. The check only converts an install that today
corrupts into a refusal; every legitimate donor's leaves are positive. `equivalence`,
`bcf-equivalence` and `multinomial-equivalence` are all bitwise-unchanged. TESTS THAT
MOVE: none break - `[[test-heteroscedastic-warm-start.R:91-92@658869ac]]` exercises the live-sourced
arm with a legitimate donor and keeps passing. That file is the natural home for one new
`expect_error` beside its existing slot-mismatch refusals at `[[test-heteroscedastic-warm-start.R:149@658869ac]]` and `[[test-heteroscedastic-warm-start.R:169@658869ac]]`. No exact
gate or SBC arm covers this and none should: it is a validation law, not a posterior.
RECOMMENDATION: (a), and land it independently of the other two. It is the only one of the
three that is live on the documented R surface, it costs five lines, it moves no draws,
and the fix is a copy of the sibling arm's own body. Priority: highest of the three.

### (iii) forests_[0] is hardcoded in applyNewData and recoverTreeParameters

WHAT THE TWO ROUTINES DO. They are the two phases of a whole-data replacement.
`recoverTreeParameters` (`[[chain.hpp:2420@658869ac]]`) reads every tree's leaf parameters against the
CURRENT fits and partitions, before the store moves. `applyNewData` (`[[chain.hpp:2448@658869ac]]`)
then swaps the response, resizes the per-observation storage, remaps split indices onto
the rebuilt cut grid, and collapses whatever is left invalid. Both open with
`Forest<L, ResidT>& forest = forests_[0];` (`[[chain.hpp:2421@658869ac]]`, `[[chain.hpp:2453@658869ac]]`) where the siblings around
them loop - `repartitionTrees` (`[[chain.hpp:2320@658869ac]]`), `revalidateTrees` (`[[chain.hpp:2269@658869ac]]`),
`rebuildFitsFromParameters` (`[[chain.hpp:2328@658869ac]]`), `dropStaleMissingDirections` (`[[chain.hpp:2526@658869ac]]`),
`forceRefreshTrees` (`[[chain.hpp:2540@658869ac]]`). On a K-forest chain, forests 1..K-1 would keep fits against
the old grid.

WHEN A K-FOREST FIT WOULD REACH THEM - and why it does not. The only caller is
`SamplerBase::setData` (`[[sampler.hpp:1163@658869ac]]`, `[[sampler.hpp:1178@658869ac]]`), whose only entry is
`bartcore_setData`, guarded by `refuseMultiForestMutation`
(`[[R_interface_bartcore.cpp:4632@658869ac]]`, guard at `[[R_interface_bartcore.cpp:2610@658869ac]]`). There is no `dbarts_sampler_setData`
in the shipped flat header, so the C API cannot reach it either. PROBE
(`brief/probe-forest0.R`): on a two-forest fit `setData` refuses with "setData does not
support a sampler that carries forest amplitudes"; on a K = 3 multinomial fit it refuses
with "$setData is not available on a multinomial sampler". `setPredictor` is ACCEPTED on
both, and correctly - it routes through `revalidateTrees`/`rebuildFitsFromParameters`,
which loop `forests_`. Save/reload does not reach these two routines at all.

CANDIDATE RULINGS.
(a) Loop `forests_` in both, matching the siblings. ~15 engine lines. But the guard's own
    comment (`[[R_interface_bartcore.cpp:2591-2600@658869ac]]`) explains that looping is NOT sufficient:
    every per-forest amplitude basis is copied at the observation count it was installed
    at, so a replacement that grows n would over-read every basis and a fixed-n one would
    silently re-pair old bases with new rows. A real lift must take the bases in the same
    call. So (a) as written buys correctness for the multinomial shape only.
(b) Leave the bodies and add an engine-level assert or early return - the guard is two
    layers away in the bridge and nothing in the engine states the precondition. ~4 lines
    plus one C++ test.
(c) Do nothing.
COST TO BASELINES: zero under all three - no shipped path executes the divergent code.

RECOMMENDATION: (b). This is the one of the three that is NOT draw-moving anywhere - the
only route in is refused two layers up, by two different messages, on both multi-forest
shapes. Correcting the bodies would be writing untested code for a door that is shut and
whose full opening needs a basis-carrying setData that nobody has asked for. Record the
precondition where it is violated and move on. Priority: lowest.

NOTE ON THE REVIEW'S "ALL THREE MOVE DRAWS": true of none of them on the documented
surface. (ii) reaches documented R (`installTrees`) but only with a hand-edited state;
(i) needs an undocumented control attribute; (iii) is unreachable. Under the recommended
options the slice moves no draws, needs no re-record, and need not land last.

## DECISION 2 - the docs/ citation policy for shipped comments

THE COUNT, at the tip. 316 lines under `R/`, `src/`, `inst/`, `man/` cite a `docs/` path,
across 80 files. Split: `R/` 109 (13 files), `src/` 151 (`src/bartcore` 102 in 11 files,
top-level `src/*.cpp|hpp` 42, `src/include` 3, `src/misc` 1), `inst/include` 2,
`inst/tinytest` 54 (46 files), `man/` 0. The review's "262 tree-wide" is exactly
109 + 151 + 2 - shipped code excluding the test suite; adding `inst/tinytest` gives 316.
ALL 316 cite `docs/design`. ZERO cite `docs/plans`. Verified by grep for both strings; the
standing rule against plan references in shipped comments is already fully satisfied, so
alternative (a)'s "strip only the plans ones" is a no-op today.
40 distinct design documents are cited, all 40 exist and all 49 files in `docs/design` are
git-tracked; 74 of the citations carry a section or number anchor. Heaviest: `bcf.md` 49,
`multinomial.md` 38, `negative-binomial.md` 28, `ordinal.md` 25, `survival.md` 19.

THE DANGLING-POINTER FACT. `.Rbuildignore` line 17 is `^docs$`. Built the tarball to
confirm: `dbarts_1.0-0.tar.gz` (1,310,085 bytes, no built vignettes) contains 0 paths
under `dbarts/docs`. So a CRAN user reading the installed source follows the pointer to
nothing, while a reader on github.com/vdorie/dbarts - the canonical source, where all 49
files are tracked - follows it to the document.

THE C2 NARRATIVE LINES. One site, not 28 sites: the block comment at
`[[src/R_interface_bartcore.cpp:6272-6311@658869ac]]` (40 lines) above `stateFormatVersion = 3` at
`[[src/R_interface_bartcore.cpp:6312@658869ac]]`. Load-bearing core is the registry rule at `[[src/R_interface_bartcore.cpp:6285-6296@658869ac]]` (append-only block names,
frozen encodings, when the version and the read floor bump) plus the attributes clause at
`[[src/R_interface_bartcore.cpp:6305-6307@658869ac]]`. Narrative: `[[src/R_interface_bartcore.cpp:6272-6284@658869ac]]` (13 lines) recounting three pre-release format
iterations, opening "The shipped format (version 2)" 40 lines above the literal 3, and
`[[src/R_interface_bartcore.cpp:6296-6303@658869ac]]` (8 lines) on what versions 2 and 3 changed. Both are provenance for a package
that has not shipped - the block says so itself at `[[src/R_interface_bartcore.cpp:6310@658869ac]]`, "no release ever shipped
format 3". About 21 to 26 lines depending on where the sentence boundaries are cut; the
review's ~28 is right to within the arithmetic of splitting shared lines.

ALTERNATIVES.
(a) Keep the design citations, strip nothing else. Zero work - there are no plan
    citations to strip. Buys: the comment stays a pointer to the reasoning for a reader
    on GitHub, where anyone modifying the engine is. Costs: 262 shipped-code lines a CRAN
    reader cannot resolve, with nothing in the comment saying the target was stripped.
(b) Strip every `docs/` citation and inline the constraint. Buys: every shipped comment
    is self-contained. Costs: 316 edits over 80 files, each a judgement about what in the
    design document is the CONSTRAINT and what is derivation; done badly it turns a
    one-line pointer into a paragraph.
(c) Ship the design documents. Measured: `docs/design` is 1.5 MB over 49 files, 594 KB as
    a gzipped tar member, against a current tarball of 1.25 MB - roughly 1.85 MB total,
    comfortably inside CRAN's 5 MB guidance. But `docs/` cannot ship where it stands: R's
    own check carries "Non-standard file/directory found at top level" (confirmed in the
    installed `tools` package's message database), so it would have to move to
    `inst/docs/design`, at which point all 316 citation strings are wrong by one path
    component in the source tree and right in the installed tree, or vice versa.

RECOMMENDATION: (a), with one qualifier and one carve-out.
Qualifier: every citation is a design citation and every target exists, so the policy
question is only "may a shipped comment cite a path the tarball strips". It may - the
audience for an engine comment is someone reading the engine, and that person is on
GitHub. Say so once, in a conspicuous header, rather than 262 times.
Carve-out: the C2 block is a separate problem from the citation policy and should be fixed
regardless. Cut `[[src/R_interface_bartcore.cpp:6272-6284@658869ac]]` and `[[src/R_interface_bartcore.cpp:6296-6303@658869ac]]`, keep the registry rule and the attributes
clause, and change "The shipped format (version 2)" to name 3 as the first shipped format.
That is a ~26-line deletion in one file, no draws, no tests.

## DECISION 3 - restore xbart's k-grid sort?

WHAT 0.9-34 DID. `git show main:[[R/xbart.R:67-69@658869ac]]` computes `kOrder <- order(k, decreasing =
TRUE)`, its inverse, and sorts `k` before the `.Call`; `[[R/xbart.R:127-132@658869ac]]` un-permutes the result
array and `k` afterwards, so the reported array is in the user's order. WHAT THE TIP DOES:
`[[R/xbart.R:424-428@658869ac]]` builds `cells <- expand.grid(iBase, iPower, iK, iTrees)` and sweeps in
that order; `[[R/xbart.R:689@658869ac]]` gives the first cell of each tree count a fresh sampler burned
`n.burn[1]` and `[[R/xbart.R:692@658869ac]]` gives every later cell a warm start burned `n.burn[2]`.
Grepped `docs/`, `inst/NEWS.Rd` and the landing notes for any record of the removal: none.
Grepped 0.9-34's `man/xbart.Rd` and `src/crossvalidate.cpp` for a stated reason: none.
The nearest thing is `[[xbart.Rd:5@658869ac]]`, "sharing burn-in between parameter settings", and the
`n.burn` argument's second element - the warm start is documented, the sort never was.

IS THE ORDER DEPENDENCE WARM-START BIAS OR STREAM NOISE? Measured
(`brief/probe-xbart2.R`): Friedman data, n = 120, 4 folds, 25 trees, n.reps = 12,
24 seeds, `n.burn = c(20, 5)`, comparing k = c(2,8), k = c(8,2), and each k run alone.
The control validates the pairing exactly: the first-in-sweep cell reproduces the run-alone
cell to 0.0000 with sd 0, since it is the same fresh sampler off the same stream. Against
that zero:
  k = 8 warm after k = 2, minus fresh:  -0.0370  (se 0.0056, t = -6.59)
  k = 2 warm after k = 8, minus fresh:  +0.0204  (se 0.0101, t = +2.03)
across-seed sd of the fresh cells: 0.021 (k=8), 0.038 (k=2). Both signs point the same
way - the warm cell is dragged toward its predecessor's fit - so this is warm-start bias,
not stream noise, and the memo's probe C (8.4% at n.reps = 1) was mostly single-draw
noise. It is a transient: at `n.burn[2]` of 5 / 20 / 100 the k = 8 shift is
-0.0370 / -0.0178 / -0.0095, and at the SHIPPED default `n.burn = c(200, 150)` with
n.samples 50 it is -0.0087 (0.27% of the loss level, t = -2.19) and +0.0082 (t = 0.55, not
significant). So on defaults it is real, systematic, and about a tenth of one seed's noise.

THE GATE. `[[benchmarks/R/equivalence.R:1387-1400@658869ac]]`, scenario `xbart`, runs
`k = c(1, 3)` - ASCENDING. A decreasing sort reorders it to `c(3, 1)`, changing which cell
is fresh and which is warm and shifting the RNG stream, so the scenario's recorded channel
(20 seeds x 8 loss cells, `loss.1..loss.8`) MOVES. Restoring the sort is therefore NOT
bitwise for the gate: `equivalence-5a3bc276.rds` needs a re-record, with the other 42
scenarios expected identical (`EQUIVALENCE_SCENARIOS` allows a filtered check to prove
that before the full record). `bcf-equivalence` and `multinomial-equivalence` are
untouched. No exact gate and no SBC arm covers xbart at all.

ALTERNATIVES.
(a) Restore the decreasing sort. Buys: results invariant to the order the user typed the
    grid in - two callers who list the same k values differently get the same array - and
    0.9-34 parity. Costs: ~8 R lines; one baseline re-record; and it does NOT remove the
    bias, it only fixes its direction. Under decreasing k the sweep runs most-shrunk to
    least-shrunk, so the flattered cells are the small-k ones; under the current ascending
    convention the flattered cells are the large-k ones. Worth saying plainly in the
    decision: the sort buys determinism, not unbiasedness.
(b) Keep `expand.grid` order and document the order dependence in `xbart.Rd`. Costs ~4
    doc lines, no code, no re-record, nothing moves. Buys honesty and nothing else; the
    user is told their answer depends on how they typed the grid, and told to sort it
    themselves if they care.
(c) Drop warm-starting between cells - every cell gets a fresh sampler at `n.burn[1]`.
    Removes the bias entirely and makes the answer order-invariant for free. Costs: the
    whole point of `n.burn`'s two-element form, documented at `[[xbart.Rd:54-55@5a3bc276]]`, plus a
    large runtime increase (every cell pays 200 burn sweeps rather than 150); re-records
    the `xbart` scenario; and `n.burn[2]` becomes a dead argument to deprecate.

RECOMMENDATION: (a), with (b)'s documentation sentence added rather than instead. The
order-invariance is a genuine surface property - the same call written two ways should
give one answer - and this is the pre-release window in which the public face is settled.
The measured bias at defaults (0.27%) is too small to justify (c)'s cost, and too
systematic to leave silently order-dependent under (b) alone. Cost is ~8 R lines, one
sentence in `xbart.Rd` naming the sweep order and the warm start's residual bias, and one
re-record of `equivalence-5a3bc276.rds` in which only the `xbart` row may move.
