# forest-defaults-by-kind: a forest's defaults go by its kind, and a forest is selected by its label

Status: PUSH 1 LANDED 2026-10-07 (dec-B274, dec-B276; dec-B241 and dec-B246 as dec-B274 restates them):
8ae00bdb to 16069daf. Push 2, selecting a forest by its label, is not built. Followed push 3 of
[written-surface.md](written-surface.md). Amended 2026-10-07 after the critique of
the multiplier law and dec-B281 and dec-B282: push 1 gains the interim refusal of the held shapes the tip
gets wrong (step 1.10).

agent: two pushes. Push 1 (the defaults): opus implementer for the R code, sonnet for the respelled tests
and the help once the code is fixed, opus reviewer told to refute. Push 2 (selection by label): sonnet
implementer, opus reviewer. The reason for opus on push 1: it is about 250 lines in the one block that
decides what model a call fits, and every slip there is a fit with another tree count or tree prior and no
message. The three reviews of written-surface's second push each found such a fit one step to the side of
the last; this plan lists the same steps to the side for this slice and tests each (Verification).
rng: stated per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- POSTERIOR-CHANGING for one sequence: a model in which no forest is without a basis, created with no tree
  count and no tree prior given to the fitting function, whose first forest leaves any of `n.trees`, `base`
  and `power` unstated. That forest took the fitting function's default for each (75, 0.95, 2) and now
  takes the default of its kind (50, 0.25, 3). The doors: a formula whose forests all have a basis, a
  `forests` list of such forests, a data object whose `bases` are all given, at `bart()`, `dbarts()` and
  `dbartsSpec()`. The oracle is the same call with the three numbers stated on that forest, fitted on the
  base build: identical draws, five pairs run for this plan. No such call is in a test's pinned draw, a
  baseline, an exact gate, a help example, a vignette or a consumer (searched, and run where it could be).
- NEUTRAL for every other fit that both builds accept: one forest, and several with the forest that has no
  basis first. 17 seeded pairs, the 55 and 15 scenarios of the two equivalence baselines that reach this
  code, and bartCause's `bcf()` were identical with a prototype of the rule in the slice's place.
- Refused where accepted: a tree count, tree prior, leaf prior, `interactions` or `blocks` given to the
  fitting function of a model in which every forest has a basis. 44 creations in 6 test files do it, each
  with a count and nothing else; none elsewhere. Respelled onto the first forest each keeps its draws.
- Accepted where refused: a `forests` list, and a data object's `bases`, whose forest with no basis is not
  the first. No draw existed before; the evidence is in Verification.
- Refused where accepted: `amplitude = fixed()` on a forest the tip holds at a value the help does not
  state (a basis of two columns that is not the second forest, a basis of three or more columns anywhere),
  and a swap that changes the width of a held forest's basis. No creation in the suite, the benchmarks or
  a consumer does the first (traced on push 3's build: 49 held forests, all of the two shapes kept); three
  swaps in one test file do the second.
- Push 2 moves no draw: it adds a way to name a forest.
window: pre-release, directly after written-surface push 3 and before the sd unit, the kind by class and
the multiplier law, so that the tips on which written order decides a model are written-surface's alone.
Serial with any other work in [`resolveForests`](../../R/model.R), [`forestParams`](../../R/model.R), the
multi-forest block of [`resolveSamplerSpec`](../../R/spec.R),
[`buildHostSamplerCall`](../../R/bart.R), [`resolveForestIndex`](../../R/bartcore.R),
[`resolveForestSelection`](../../R/generics.R) or the per-forest methods of the sampler. No engine, bridge
or header file changes, so it shares only R/dbarts.R and man/dbartsSampler-class.Rd with
[leaf-conversions.md](leaf-conversions.md) and
[cross-family-state-install.md](cross-family-state-install.md), other functions and items, and may land
before, between or after them.
budget: ~1620 lines changed (R ~430; tinytest ~880, of which ~680 in two new files and ~200 respelled or
changed; help ~180; design note, architecture, public-surface, TODO and the two indexes ~130), upper
figure 3050 (R 850, tinytest 1600, help 350, records 250). By push: 1 ~1070 (upper 2050), of which the
interim refusal of step 1.10 is ~120; 2 ~550 (upper 1000). The design estimated 600 and planned for 1100;
written-surface's pushes landed at 1.8 and 2.7 times their plans.

## Goal

A forest's default tree count and tree prior depend on whether it has a basis and on nothing else: the
forest with no basis takes the fitting function's, a forest with a basis takes the multiplied forest's (50
trees, `base = 0.25`, `power = 3`), wherever either is written, in a formula, in a `forests` list and in a
data object's `bases` alike. Where no forest is without a basis, what the fitting function states for
that forest (a tree count, a tree prior, a leaf prior, `interactions`, `blocks`) is refused by name and
the model is accepted. Every method that takes a forest takes its label as well as its position. Until
the multiplier law lands, a coefficient is held only where the tip holds it at the value the help
states.

## Context

Measured at the tip (82b00b7d, written-surface pushes 1 and 2 landed) on the shipped build, R 4.6.1: 150
rows; x1, x2, x3, dose, age, a 0/1 z; one chain. "Count C" is a count named in the call, 15 in the probes.
The spellings below are the tip's, with a tilde on a basis; push 3 drops it.

- The suite. 227 files run on a shipped build, 16143 results at push 2's landing. Rerun for this plan
  without the three files that ask for more than two threads: 224 files, 16065 results, none failing, 100
  seconds in one process.
- Defaults today. Sources: F the fitting function's default (75 trees; `cgm(2, 0.95)`), S a value stated
  to the fitting function, O the forest's own argument, M the multiplied forest's default (50; 0.25, 3).
  A forest's own `n.trees`, `base` and `power` each win where given, at every door.

  | model, as written | door | forest 1 | later forests |
  |---|---|---|---|
  | one forest: plain terms, `forest(x1 + x2)`, `list(forest())` | `bart`, `dbarts`, matrix, `dbartsSpec` | F or S | |
  | no-basis forest first, then one or two with a basis | formula at `bart` and `dbarts`; list on a formula, on a matrix, in `dbartsSpec`; `bases = list(NULL, w)` | F or S | M |
  | no-basis forest written second or last | formula | it is forest 1: F or S | M |
  | the same | list, matrix, `dbartsSpec`, `bases = list(w, NULL)` | refused: "forest 2 needs a 'basis'" | |
  | every forest with a basis | every door, `bases = list(w1, w2)` included | F or S: 75, 0.95, 2 with nothing stated; C with a count named | M |
  | two forests with no basis; one forest with a basis | list, data object | refused | |

  So `y ~ forest(x1, basis = ~ dose) + forest(x2, basis = ~ age)` has 75 and 50 trees, 15 and 50 under a
  control that names 15, and the two terms swapped are 75 and 50 again with the forests exchanged.
  `tree.prior = cgm(3, 0.8)` beside it gives the first forest 0.8 and 3; a first forest stating
  `base = 0.4` under that prior has 0.4 and 3.
- Where the numbers live. The bridge takes forest 1's count from the control and its tree prior from the
  model, and every later forest's from that forest's own eight numbers
  ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)); forest 1's first three of the eight are
  carried and unread, and say 50, 0.25, 3 whatever forest 1 runs under. Which of a forest's two size
  channels its `sd` reaches already goes by whether it has a basis, at any position
  ([`forestParams`](../../R/model.R)). The engine runs a forest with no basis at any position: a
  prototype of the rule fitted, ran, restored, copied, swapped bases, predicted and extracted on one.
- A count given twice. `bart(y ~ forest(x1 + x2, n.trees = 7) + ..., n.trees = 15)` is refused
  ([`refuseTreeCountGivenTwice`](../../R/formulaTerms.R)), by the formula's text: with `basis = NULL`
  written on that forest it is accepted and has 7 trees. At `dbarts()` a control naming 15 beside a
  no-basis forest stating 7 gives 7, as man/forest.Rd says.
- What is stated. `bart()` builds its own control, which always names a count, and always hands
  `dbarts()` a tree prior, so `dbarts()` cannot tell what `bart()`'s caller wrote. A control records the
  arguments its constructor was given ([`controlSuppliedSlots`](../../R/dbarts.R)):
  `dbartsControl(n.trees = 75L)` records `n.trees`, through `do.call` and through a wrapper that forwards
  its own default too. So a named 75 is told from an untouched control at the tip. Not told apart: a slot
  edited back to 75, and `new("dbartsControl")`.
- The leaf prior. On any model of several forests the fitting function's `leaf.prior` is accepted only
  where it changes nothing: `normal`, `normal(k = 2)`, `bart`'s `k = 2`. `normal(k = 3)`, `normal(sd = 2)`
  and a linear leaf are refused whatever the forests are.
- `interactions` and `blocks` given to the fitting function are forest 1's, whatever it is; stated there
  and on forest 1 both, refused ("declared both at the top level and on the first forest").
- A control taken from a sampler and given to another fit carries the first fit's count:
  `dbarts(y ~ x1 + x2, d, control = s$control)` has 7 trees when the no-basis forest of `s` stated 7.
- `print` of a fit shows one line, `n.trees:`, the control's, for any number of forests.
- Labels. `names(forests)` where the list has any, an unnamed entry being `""`; a formula's forests and an
  unnamed list have none. [The label rule](written-surface.md#the-label-rule) is push 3's.
- Selecting a forest. Nine methods of the sampler take one (`getLeafPrior`, `getK`, `getForestFits`,
  `getForestAmplitudes`, `getForestVariableCounts`, `getTrees`, `plotTree`, `setForestWeights`,
  `setForestBasis`), all through [`resolveForestIndex`](../../R/bartcore.R); `getTrees` alone takes
  several. On a fit, `extract` with `type = "forest"`, `"k"` and `"leaf.prior.sd"` and `predict` with
  `type = "forest"` take one or several ([`resolveForestSelection`](../../R/generics.R)), and
  `extract(type = "trees")` and `plotTree` hand theirs to the sampler. Run on three forests, the second
  named `dose` in the list:

  | `forest =` | the nine methods, and `extract(type = "trees")` | `extract` and `predict`, the other types |
  |---|---|---|
  | `2L`, `2` | forest 2 | forest 2 |
  | `"dose"` | refused: "must be coercible to type: integer" | refused: "must name one of 'forest1', 'forest2', 'forest3'" |
  | `"forest2"` | refused, the same | forest 2 |
  | `"2"` | forest 2 | refused, the same |
  | `TRUE`; `factor("dose")`; `list(2)` | forest 1; forest 1, by the factor's code; forest 2 | the same |
  | `NA`, `""`, `0L`; `4L`; `2.5` | refused, three texts | refused |

- Lists given by position. `$setLeafPrior(forests = )` refuses a name that is not its position's label.
  `predict(bases = )` ignores names: `list("factor(z)" = b1, dose = b2, forest1 = b3)` is read as forests
  1, 2, 3.
- What states an argument beside forests that all have a basis, counted by running the suite on a
  prototype that logs it: 44 creations in 6 test files (32 in test-bcf-family.R, 5 in test-bcf-creation.R,
  4 in test-formula-terms.R, one each in test-forest-arguments.R, test-forest-basis-r5.R and
  test-predict-blend.R), every one a count and nothing else, named on the control in 43 and as `bart()`'s
  `n.trees` in one. None in benchmarks/R:
  its 11 call sites in 10 scripts each put the no-basis forest first (read). None in a help example (the
  eight pages with one, run) or a vignette (none fits such a model).
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) builds a list of two forests, the first with no
  basis, states `n.trees` on both and hands its tree prior to `dbarts()`; stan4bart (bartcore, a9d081b),
  treatSens (dbarts-1.0, aecec71) and bairrtt (main, 3f57f61) declare one forest. None passes a string as
  a forest.
- Retaken on push 3's build (4eeaf03d, 2026-10-07), the suite run with every creation traced: 230 files
  without the three that ask for more than two threads, 17121 results, none failing; 1186 models of
  several forests; 66 creations whose first forest has a basis, in 6 files (test-bcf-family.R 54,
  test-bcf-creation.R 5, test-formula-terms.R 4, one each in test-forest-basis-r5.R,
  test-forest-labels.R and test-predict-blend.R). Which of the 66 state a count was not retraced; step
  1.9 counts them on the landed tip.
- A held coefficient goes by position. `amplitude = fixed()` holds a forest's coefficients where the
  engine starts them ([`AmplitudeState`](../../src/bartcore/combiner.hpp),
  [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp)): forest 1's at 1, forest 2's at 0 for its
  first column and 1 for the rest, every later forest's at 1. Read off push 3's build: a two-level
  factor held as the second forest is (0, 1), the help's value; held as the third, or as the first where
  no forest is without a basis, it is (1, 1), a function added on every row; three levels held second
  are (0, 1, 1); `cbind(1 - z, z)` handed through `dbartsData(bases = )` and held second is (0, 1). A
  forest with no basis held as the second forest is held at 0 and leaves the model: R cannot write that
  at the tip and can once step 1.1 lands. One numeric column is refused already
  ([`refuseHeldOneColumn`](../../R/model.R), dec-A171). `$setForestBasis` on a held two-level forest takes
  a three-level factor and holds it at (0, 1, 1).
- What holds a coefficient, traced on push 3's build: the suite creates 49 held forests, 17 with no basis
  as forest 1 of 2 and 32 with two columns as forest 2 of 2 (one of the 32 two numeric columns, in
  test-forest-basis-terms.R). bartCause's `bcf()` holds its first forest, which has no basis, under
  `update.a = FALSE` and its second under `update.b = FALSE`, and that second forest's basis is
  `cbind(1 - z, z)`, two NUMERIC columns through `dbartsData(bases = )`; its test-14-bcf.R fits both.
  The benchmarks hold `factor(z)` as the second forest and a forest with no basis as the first.

## The rule

A forest is of one of two kinds: it has a basis or it has none. The kind is read in one place, from the
bases the model ends with, after a formula's forests, a list's declarations and a data object's `bases`
have been put together ([`resolveSamplerSpec`](../../R/spec.R)). The plain forest is the one forest of a
model of several that has no basis. A model has at most one, and may have none.

| | the plain forest | a forest with a basis |
|---|---|---|
| tree count | the fitting function's: `bart`'s `n.trees`, the control's at `dbarts` and `dbartsSpec`; 75 | 50 |
| tree prior, `base` and `power` | the fitting function's `tree.prior`; 0.95, 2 | 0.25, 3 |
| `interactions`, `blocks` | the fitting function's or its own; both given is refused | its own |
| leaf prior `sd` | by kind already; unchanged | by kind already; unchanged |

- A forest's own `n.trees`, `base` and `power` each govern that forest, as today. Each left out takes the
  default of the forest's kind, one by one: `forest(basis = dose, base = 0.5)` has 50 trees and
  `power = 3`.
- Position decides nothing. In a formula the plain forest is forest 1 wherever written and the others
  follow as written, as today; a list and a data object keep their own order, and the plain forest may
  stand anywhere in them. Two forests of one kind written in the other order are the same model and not
  the same draws: forests are swept in order.
- Stated. At `bart()`: `n.trees` named in the call, whatever its value, or a `control` whose constructor
  named `n.trees` or whose count is not that constructor's default; `tree.prior`, or the retired `power`,
  `base` or `split.probs`,
  named; `leaf.prior` or `k` named. At `dbarts()` and `dbartsSpec()`: a control as above; `tree.prior`
  named; `leaf.prior` or the retired `node.prior` named. `interactions` and `blocks`: not `NULL`. A count
  edited on a control back to 75 and a control made by `new()` read as not stated.
- With no plain forest, each of the five that is stated is refused, in the order count, tree prior, leaf
  prior, `interactions`, `blocks`, before anything else is read from them; with it moved to a forest or
  left out the model is created.
- The plain forest's own count beside the fitting function's: at `bart()` refused, as today, and judged
  against the model's plain forest, not the formula's text; at `dbarts()` the forest's own governs, as
  today. The refusal there belongs to the control-migration arc (dec-B241).
- The control's count and the model's tree prior hold forest 1's, as the bridge reads them: 50 and
  `cgm(3, 0.25)` for a first forest with a basis that states none. Each forest's eight numbers start with
  the count, base and power that forest runs under, forest 1's included. A control handed from one fit
  to another carries the first fit's plain forest's count, as today, or the fitting function's own where
  no forest was plain; never a multiplied forest's.

The same shapes after the slice, in the table of Context:

| model, as written | door | forest 1 | later forests |
|---|---|---|---|
| one forest | every door | F or S | |
| plain forest first | every door | F or S | M |
| plain forest written second or last | formula | it is forest 1: F or S | M |
| the same | list, matrix, `dbartsSpec`, `bases = list(w, NULL)` | M | the plain forest F or S, the others M |
| every forest with a basis, nothing stated | every door | M | M |
| every forest with a basis, something stated to the fitting function | every door | refused, naming the argument | |
| two forests with no basis; one forest with a basis | list, data object | refused | |

## A held coefficient, until the multiplier law

Interim, by the width of the forest's basis and the forest's position, which is all every door knows at
this tip: `dbartsData(bases = )` has no record of a basis's class until
[forest-kind-by-class.md](forest-kind-by-class.md). Checked where one numeric column is refused today.

| held forest | position | the tip holds it at | this slice | afterwards |
|---|---|---|---|---|
| no basis | first; third or later | 1 | accepted | the multiplier law makes its `sd` exact |
| no basis | second | 0: the forest leaves the model | refused | accepted at 1, by the multiplier law's first push |
| a basis of two columns | second | (0, 1) | accepted | a two-level factor: unchanged. Two numeric columns: the basis itself is refused for good (dec-B282, forest-kind-by-class push 2) |
| a basis of two columns | any other | (1, 1): no contrast | refused | a two-level factor: accepted at (0, 1), by the multiplier law's first push. Two numeric columns: refused for good |
| a basis of three or more columns | any | (0, 1, 1, ...) or all 1 | refused | refused for good: the basis itself is, drawn or held (dec-B281, dec-B282, forest-kind-by-class push 2) |
| a basis of one column | any | refused already (dec-A171) | unchanged | accepted at 1, by the multiplier law's second push |
| any, at `$setForestBasis` | | the new width's values by position | another width refused | for good: after forest-kind-by-class push 2 no forest changes width |

So two of these refusals are interim (a forest with no basis second; a two-level factor away from the
second place), one is interim already (one column), and two are permanent in effect, their texts
replaced when the basis itself is refused (three or more columns; a width changed by a swap). What stays
accepted and is not a model the help describes: two numeric columns held as the second forest, at (0, 1),
the first column's term zero. It cannot be told from a factor's indicator columns at the data door, it
is how bartCause's `bcf()` holds its treatment forest, and dec-B282 closes it after bartCause has moved
to `factor(z)`.

## Selecting a forest

Written against [The label rule](written-surface.md#the-label-rule) as push 3 lands it: every forest of a
model of several has a label, labels are unique, and a list name of the form `forest<digits>` is its own
position's.

1. `forest` is a number or a string, or `NULL` where a method takes every forest. A vector where the
   method takes several today.
2. A number is a position, read by today's code with today's refusals.
3. A string is never a position. It is the forest with that label; failing that, the forest whose label
   is the same code (`"I(dose / 30)"` finds `I(dose/30)`); and `forest<i>` is forest i, the name position
   i has on every per-forest margin, so `forest = "forest2"` selects what it selects today.
4. A string that is the label of one forest and the name of another's position, or that is the same code
   as two labels, is refused, naming both.
5. `"2"` is the forest labelled `"2"`, or no forest; the refusal says a position is given as a number.
   With a list named `"2"`, `"1"`, the string `"2"` is forest 1 and the number 2 is forest 2, as base R
   reads a list by name and by position.
6. `NA`, `""`, a factor, a logical and a list are refused by name.
7. A sampler of one forest and a multinomial one have no labels: `forest<i>` is taken where the number
   is, and any other string is refused.
8. What comes back is what the position gives, to the bit. Nothing is renamed by label in this slice.

A list given forest by forest (`$setLeafPrior(forests = )`, `predict(bases = )`) is read by position, and
a name on an entry that is neither its position's label nor `forest<i>` is refused.

## Refused forms, with their texts

Base R's style, as in written-surface. `<arg>` is the argument as the caller wrote it; the push that adds
each is in brackets.

    [1] 'n.trees' given to the fitting function is the tree count of the forest with no basis, and every forest of this model has a basis; state a count on a forest, as forest(x1, basis = a, n.trees = 100)
    [1] the control's 'n.trees' is the tree count of the forest with no basis, and every forest of this model has a basis; state a count on a forest, as forest(x1, basis = a, n.trees = 100), and leave it out of dbartsControl()
    [1] '<arg>' given to the fitting function is the tree prior of the forest with no basis, and every forest of this model has a basis; state it on a forest, as forest(x1, basis = a, base = 0.25, power = 3)            (tree.prior, power, base, split.probs)
    [1] '<arg>' given to the fitting function is the leaf prior of the forest with no basis, and every forest of this model has a basis; state a forest's size on the forest, as forest(x1, basis = a, sd = 2)            (leaf.prior, k, node.prior)
    [1] 'interactions' given to the fitting function is a constraint on the forest with no basis, and every forest of this model has a basis; state it on a forest, as forest(x1, basis = a, interactions = interactions(max.order = 1))            ('blocks' the same)
    [1] forests 1 and 3 have no 'basis': a model has one forest with no multiplier, and every other forest states a 'basis'
    [1] 'interactions' is given to the fitting function and to the forest with no basis, which are the same constraint; give one            ('blocks' the same)
    [1] 'n.trees' is given to the fitting function and to the forest with no basis, which are the same count; give one
    [1] forest 2: amplitude = fixed() on a forest with no basis is not supported yet where it is the second forest; it would hold the forest at zero. Put the forest with no basis first, or let the coefficient be drawn
    [1] forest 3: amplitude = fixed() on a basis of two columns is not supported yet unless the forest is the second; it would hold both coefficients at 1. Put the forest second, or let the coefficients be drawn
    [1] forest 2: amplitude = fixed() on a basis of 3 columns is not supported; it would hold every column but the first at 1. Let the coefficients be drawn
    [1] $setForestBasis cannot change the width of forest 2's basis (2 to 3): its coefficients are held (amplitude = fixed()), and the held value is defined for that width only; make a new sampler
    [2] 'forest' names no forest of this model: "age"; its forests are "forest1", "scale(age)", "I(dose/30)"
    [2] 'forest' names no forest of this model: "2"; its forests are "forest1", "scale(age)", "I(dose/30)"; a position is given as a number, forest = 2
    [2] 'forest' names no forest of this model: "dose"; this sampler's forests have no labels, so select one by position
    [2] 'forest' ("forest2") is the label of forest 3 and the name of position 2; select by position, as forest = 3
    [2] 'forest' ("a +b") is the label of forests 2 and 3 ("a+b", "a + b"); give one exactly, or select by position
    [2] 'forest' must not be NA or an empty string
    [2] 'forest' must be a number, the forest's position, or a string, its label; not a factor            (a logical, a list)
    [2] 'bases' names forest 2 "age", and its label is "dose"; give the bases in the forests' order

Replaced: "forest 2 needs a 'basis': the amplitudes multiplying it are what distinguishes it from the
first", in its two variants, by the sixth; "is declared both at the top level and on the first forest",
by the seventh. The eighth is written-surface's text without the forest quoted, for a plain forest the
formula's text does not show to be one; where it shows, the text and its quote are unchanged. Kept as
they are: the at-least-two refusal, what a model of several forests refuses whatever its forests (a DART
tree prior, `split.probs`, a `k` other than 2, a named `sd`, other leaves), and a number's refusals.
The tree prior's text gives `base` and `power` because `forest()` has those two until `tree.prior`
arrives on it; that slice changes the example.

## Constraints

- Every fit with one forest, and every fit of several whose plain forest is first, draws what it draws on
  the base build. No baseline is re-recorded and no snapshot file regenerated.
- No change to src/, to the flat C header, to the stored state, or to the order and meaning of the eight
  numbers per forest. The bridge keeps reading forest 1's count from the control and its tree prior from
  the model.
- A forest's default `sd`, the unit of `sd`, and the two size channels are untouched. No held value
  and no law of a held forest changes: step 1.10 refuses, and builds nothing.
- No test this slice adds fits a factor of three or more levels or a basis of several numeric columns,
  the refused rows of step 1.10 apart: both are refused two slices on (dec-B281, dec-B282), and a test
  written on one now is deleted then.
- The order of a model's forests is untouched at every door: forest i is the forest it is today wherever
  both builds accept the model.
- `forest()`'s formals are untouched. Nothing of `tree.prior` or `leaf.prior` on `forest()`, of
  `updateBasisScale`, or of `dbarts()`'s own `n.trees` rides along.
- No reader, `extract` or margin is renamed by label; dec-B275's list is the multiplier-law slice's.
- A number given as `forest` does what it does today, in every method, with today's texts.
- The two mutation-battery anchors in [`writeForestSpreads`](../../R/dbarts.R) are not edited.
- No consumer is edited. bartCause's `bcf()` holds the two shapes step 1.10 keeps.
- Each push leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Pushes

Two, each gated on its own and each a coherent tip.

1. The defaults (dec-B274). POSTERIOR-CHANGING for the one sequence of the `rng:` line; the design note
   lands here. With it, because this is the push that first lets a forest with no basis stand second:
   the interim refusal of the held shapes the tip gets wrong.
2. Selection by label (dec-B276). No draw moves.

What waits on written-surface push 3. Push 2 whole: before it a formula's forests have no label, two
forests may share one, and a second forest may be named `forest1`. Of push 1: its tests and help are
written in push 3's spelling of a basis, it edits test files push 3 edits, and the count of calls to
respell is retaken on push 3's tip. What can be written before: push 1's R code, which reads no label and
nothing of how a basis is written, and push 2's one function with its table of cases over a vector of
labels.

## Steps

"Fails today" is what the tip does where the test expects otherwise. New functions are named for the
reader's sake. Calls are in push 3's spelling; a fixture is 150 rows with x1, x2, x3, dose, age and a 0/1
z, and "the table" is each forest's count from the engine and its base and power from its record.

### Push 1: the defaults

1.1 One finder of the plain forest, and a plain forest anywhere. A function of the resolved bases
    (`plainForest`) gives the position of the forest with no basis, none, or the refusal of two; it is the
    only place a forest's kind is asked, by [`resolveForests`](../../R/model.R) for a list and by
    [`resolveSamplerSpec`](../../R/spec.R) for a data object. The refusal of a forest past the first with
    no basis goes from both. Tests, a new file test-forest-defaults.R:
    `forests = list(forest(basis = dose), forest())` on a formula, on a matrix and in `dbartsSpec`, and
    `dbartsData(bases = list(dose, NULL))`, are created (fail today: refused) and each is the reversed
    list's model forest for forest: the same eight numbers per forest, the same columns, the half-Cauchy
    channel on the plain forest wherever it stands (`getLeafPrior()$prior.sd.of`); it runs, restores,
    copies, takes `$setForestBasis` on either forest and predicts, `predict` at the training rows
    reproducing the training fit; `list(forest(), forest())`, `list(forest(), forest(basis = z),
    forest())` and `bases = list(NULL, NULL)` are refused with the two-forests text naming the positions;
    one forest with a basis keeps the at-least-two text.
1.2 Defaults by kind. In [`resolveSamplerSpec`](../../R/spec.R): the fitting function's count and tree
    prior are read once, before any forest's own statement is applied; forest 1 takes its own statement,
    else its kind's default, onto the control and the model; [`forestParams`](../../R/model.R) gives
    every forest, forest 1 included, the count, base and power it runs under. Tests: the table for the
    seven shapes of "The rule" at `bart`, `dbarts` on a formula, the list on a formula and on a matrix,
    the data object and `dbartsSpec`, with nothing stated and with a count and a tree prior stated where
    a plain forest takes them (fails today for every row whose first forest has a basis); an all-basis
    model with nothing stated is, draw for draw, the same model with `n.trees = 50, base = 0.25,
    power = 3` written on every forest, at the three doors and under probit (fails today: 75, 0.95, 2 on
    the first); `forest(basis = dose, base = 0.5)` first has 50 trees and power 3;
    `list(forest(basis = dose, n.trees = 9), forest())` under a control naming 15 has 9 and 15, and under
    `tree.prior = cgm(3, 0.8)` the plain forest has 0.8 and 3 and the first 0.25 and 3;
    `list(forest(basis = dose), forest())` under a control naming 15 is, draw for draw, the list with
    `n.trees = 50, base = 0.25, power = 3` on the first forest and `n.trees = 15, base = 0.95, power = 2`
    on the second under an untouched control; a formula's two multiplied terms swapped have the same
    table with the rows exchanged, and other draws.
1.3 What is stated. [`bart`](../../R/bart.R) records, for the sampler it builds, which of its own
    arguments its caller stated (the count, flat or through a control that speaks for it as
    [`mergeFrontDoorControl`](../../R/dbarts.R) judges; the tree prior or a retired shorthand of it; the
    leaf prior or `k`), on the control it hands over, and [`resolveSamplerSpec`](../../R/spec.R) reads the
    record and removes it; without one it reads the control ([`controlSuppliedSlots`](../../R/dbarts.R),
    or a count that is not the constructor's default) and the call. Tests, on an all-basis model, each
    refused or created as "The rule" says: `bart()` with nothing, `n.trees = 75L`, `n.trees = 100L`, a
    partial name, `do.call`, a wrapper forwarding its own default, `control = dbartsControl()`,
    `control = dbartsControl(n.trees = 75L)`, a control edited to 30; `dbarts()` and `dbartsSpec()` with
    no control, an untouched one, a named 75, a named 20, `do.call(dbartsControl, )`, an edit to 30, an
    edit back to 75 and `new("dbartsControl")` (the last two created, pinned as the limit of what can be
    told); `tree.prior = cgm` and `leaf.prior = normal`, the defaults named. A fit's control carries no
    such record afterwards.
1.4 No plain forest. The five refusals of "Refused forms", made in one block where the plain forest is
    found, before the fitting function's `interactions` and `blocks` are read against any forest. Tests:
    each of the five at `bart`, `dbarts` on a formula, the list, the data object and `dbartsSpec`, with its
    text, the retired names each named as written (fail today: created); two given, the first in order
    refused; each refused call with the argument moved to a forest, or dropped, is created, and with the
    count and the tip's tree prior moved to the first forest it is the tip's model (the pair script);
    beside a plain forest at any position none is refused.
1.5 `interactions` and `blocks` follow the plain forest. Given to the fitting function they are resolved
    against the plain forest's columns and tree count at its position; forest 1's own stay forest 1's.
    Both given is refused against the plain forest wherever it stands. Tests: with the plain forest
    second, the fitting function's `interactions(max.order = 1)` is on forest 2 and not on forest 1 (the
    two records, and over 300 sweeps no tree of forest 2 splitting on two predictors along one path
    while trees of forest 1 do); the same for `blocks` with the fitting function's count;
    `interactions` on the first, multiplied forest and at the top are two constraints on two forests; on
    the plain forest and at the top, refused with the reworded text at both positions.
1.6 The count given twice at `bart()`. Where the model's plain forest states `n.trees` and `bart()`'s own
    was named, the refusal is made on the model, so `forest(x1 + x2, n.trees = 7, basis = NULL)` and a
    `basis` holding `NULL` are refused as `forest(x1 + x2, n.trees = 7)` is (fails today: 7 trees, the 15
    dropped). The early check on the formula's text stays for its quoted text. Test: the three spellings;
    a multiplied forest's count beside `bart()`'s stays accepted.
1.7 A control carried to another fit. Where forest 1 has a basis the control's count is a multiplied
    forest's, so the forests' record keeps beside it the count the next fit is to inherit: the plain
    forest's as it ran, or the fitting function's own where no forest was plain. Where
    [`resolveSamplerSpec`](../../R/spec.R) clears a carried control's fit attributes it first puts that
    count back on the control. A control from a fit whose plain forest is first is not touched. Tests:
    `s` an all-basis sampler; `dbarts(y ~ x1 + x2, d, control = s$control)` has 75 trees (fails without
    the step: 50), and the same all-basis model is created again (fails without it: refused for a count
    nobody wrote); a plain-second fit under a control naming 15 carries 15, and one whose plain forest
    stated 7 carries 7, as the plain-first fit of the same forests does at the tip.
1.8 `print` of a fit of several forests gives each forest's count from the engine, in order
    ([`fitSynopsis`](../../R/generics.R)): `n.trees: 50, 75`; one forest prints as today. Test: the
    line for a plain-second fit and for one forest.
1.9 Respell and repair. The 44 creations: each call drops the count from its control, or from its
    `bart()` call, and states on its first forest the three numbers the tip gave it, `n.trees = <that
    count>, base = 0.95, power = 2`, in a one-entry `forests` list where the call gave none; identical
    draws, by the pair script. A control helper that plain-first fits share keeps its count for them.

    | file | creations as run | where the count is named | respelled |
    |---|---|---|---|
    | test-bcf-family.R | 32 | `seededControlBcfFamily()`, 25 | a second control helper with no count for the three all-basis sites (`basisSampler`, `transportParams`, the loop over the gaussian fixtures); each gives forest 1 the three numbers, added to the first entry where a list is passed |
    | test-bcf-creation.R | 5 | `control`, 50 | the two lists of two factor-basis forests, the `dbartsSpec` of the same, the two fits on `dataBasisOnFirst` |
    | test-formula-terms.R | 4 | Block F's own control, 15, three; the `fit` helper's `n.trees = 3` given to `bart()`, one | the block is rewritten (below); the all-basis formula of the predict loop states the three numbers on its first term and hands the helper `n.trees = NULL` |
    | test-forest-arguments.R, test-forest-basis-r5.R, test-predict-blend.R | 1 each | each file's control helper: 10, 25, 8 | the one call |

    Block F of test-formula-terms.R, ["every forest multiplied"](../../inst/tinytest/test-formula-terms.R),
    is the pin written-surface left for this slice (its
    [Push 2: the forests of a formula](written-surface.md#push-2-the-forests-of-a-formula), step 2.6): 12
    assertions that the first written forest takes the control's 15. Turned around: nothing stated, both
    forests 50, 0.25, 3; the formula and the list of the two identical in draws; the terms swapped the
    same table and other draws; the control's 15 refused. 11 further assertions take a new expectation:
    four pins of the text "forest 2 needs a 'basis'", which is gone, and its like take
    the two-forests text and two become creations; three pins of "on the first forest" take the reworded
    text (two in test-bcf-creation.R, one in test-formula-terms.R); two pins of the eight numbers in
    test-bcf-creation.R expect forest 1's own count and prior in the first three. Run the suite first and
    repair what it shows.
1.10 The held shapes the tip gets wrong, interim. Where [`resolveSamplerSpec`](../../R/spec.R) refuses
    a held single column today ([`refuseHeldOneColumn`](../../R/model.R)), one function
    (`refuseHeldShape`) applies the table of "A held coefficient, until the multiplier law" to every
    held forest, by `NCOL()` of its final basis and its position; and where
    [`setForestBasis`](../../R/dbarts.R) refuses one column on a held forest it refuses any width but
    the forest's own. Nothing is held anew and no value changes. Tests, a new block of
    test-forest-defaults.R: every row of the table at a `forests` list on a formula and on a matrix, in
    `dbartsSpec()`, at `dbartsData(bases = )` and, where a formula can write it, as a formula term; an
    accepted row reads `getForestAmplitudes()` as 1, or 0 and 1, before and after 20 sweeps; a refused
    row has its text, and with the hold dropped is created (fail today: a plain forest second held at
    0 once step 1.1 lets it be written; a two-level factor third at (1, 1); three levels at (0, 1, 1)).
    `cbind(1 - z, z)` through the data door, held second, is created and held at (0, 1), with a comment
    that forest-kind-by-class push 2 turns the assertion around. The three swaps of
    test-forest-arguments.R that change a held forest's width become pins of the swap's text, with the
    sampler `identical()` to an untouched twin afterwards. dec-A171's pins stand.
1.11 Help and records. man/forest.Rd: the `n.trees, base, power` item (the two kinds; a forest's own
    governs; where none is plain, state them on a forest; the sentence about the forest written first
    goes), `interactions, blocks` (the plain forest's), the paragraph on terms. man/bart.Rd: the section
    "Formula Terms" and the `n.trees`, `tree.prior`, `leaf.prior`, `interactions` and `blocks` items, one
    sentence each. man/dbarts.Rd: `forests` (at most one forest with no basis, anywhere), `tree.prior`.
    man/dbartsControl.Rd: `n.trees` is the plain forest's count in a model of several. man/dbartsData.Rd:
    `bases` (one `NULL` at most, anywhere). docs/design/forest-defaults-by-kind.md with its index row: the
    rule, what is stated at each door, the one changed sequence with its oracle, what a control carries,
  the swap law's figures.
    docs/design/public-surface.md and docs/architecture.md where they give forest 1 the fit's count. TODO:
    this entry closes; `forest-prior-args` names it landed. man/forest.Rd's `amplitude` item says which
    held forests are taken until the multiplier law lands, in the table's words; the design note
    records the table with which rows are interim.
1.12 Mutations (Verification): apply each, install, run the named test, record the failing count, revert,
    `touch` the file.

### Push 2: selection by label

2.1 One reader of a selection (`selectForest`): given `forest`, the labels and the number of forests it
    returns positions, by "Selecting a forest". A number goes to today's code unchanged.
    [`resolveForestIndex`](../../R/bartcore.R) and [`resolveForestSelection`](../../R/generics.R) both
    call it, the sampler reading its labels from the forests' record and a fit from its `forest.labels`
    attribute. Tests, a new file test-forest-selection.R, on the function alone over five label sets:
    every case of the rule and every text marked [2].
2.2 The sampler. The nine methods take a label. Fixture: four forests in a list, the plain forest second,
    one named in the list, one labelled by its basis text, one on a value, with a different `n.trees` and
    `sd` each, so that every reader's answer differs by forest. Tests: for each method and each forest,
    the label, `forest<i>` and the position give identical results, and the two writers change the forest
    named and no other (`data@bases` and the forest weights before and after); `getTrees` with a vector of
    labels; `"2"`, `TRUE`, `factor("dose")` and `list(2)` refused (fail today: forests 2, 1, 1, 2); on a
    sampler of one forest and on a multinomial one, `forest1` taken and a label refused.
2.3 A fit. `extract` with `type = "forest"` (and `contribution`), `"k"`, `"leaf.prior.sd"` and `"trees"`,
    `predict(type = "forest")` and `plotTree`, on a `bart` fit of three `forest()` terms and on a
    multinomial fit's `extract(type = "trees")`. [`selectForests`](../../R/generics.R) maps a label to its
    position and never matches it against the margin's own names. Tests: label, `forest<i>` and position
    identical for each, margins named `forest<i>` as today; a vector of labels in another order returns
    that order; the suite's five uses of a `forest<i>` string keep their results; the arms that refuse
    `forest` keep their texts.
2.4 Lists by position. [`resolveForestBases`](../../R/generics.R) refuses a named `bases` entry whose name
    is neither its position's label nor `forest<i>` (fails today: read by position). Test: named in order
    accepted and identical to unnamed; two names exchanged refused; `$setLeafPrior(forests = )` with a
    list in label order accepted and in another order refused, as today.
2.5 Help. [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) and the eight other methods'
    `forest` item and docstrings, which `tools/check-rc-codoc.R` holds to each other; man/bartBT.Rd's
    `forest` and `bases` items; man/plotTree.Rd; man/forest.Rd where push 3 describes labels: a string is
    a label, a number a position, `forest<i>` names position i. docs/design/forest-defaults-by-kind.md
    gains the selection rule.
2.6 Mutations, as 1.12.

## Verification

Every push, against the slice's own library (`R CMD INSTALL -l <lib> .`, `R_LIBS=<lib>` on every call;
check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series, at most two
cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`: unchanged and passing (nothing under src/ moves).
- The full tinytest suite on the shipped build, in one process, counted per file: no failure, no file
  stopping, and at least the base build's count plus the new files' assertions less Block F's 12; the
  landing note gives the figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted per scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds` (about twelve minutes at two cores), 15 against
  `bcf-equivalence-1b7d730c.rds`, 11 against `multinomial-equivalence-80b1c8d4.rds`. Nothing is
  re-recorded: no scenario is in the changed sequence.
- The pair script (below), old side on the base build, new side on the slice's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on
  a tarball from a clean copy.
- The consumers, each suite whole against a private install of the slice, none failing and none edited:
  bartCause on dbarts-1.0 (1412 expectations at its last run; its test-14-bcf.R fits `update.a = FALSE`
  and `update.b = FALSE` together, the two held shapes step 1.10 keeps), stan4bart on bartcore (582),
  treatSens on dbarts-1.0 (306), bairrtt on main (207). The last three declare one forest with
  `forest(n.trees = )` and nothing else of this surface (searched again 2026-10-07).

Push 1 in addition: every gate `.github/workflows/exact-gates.yaml` lists, in `quick`, unchanged (the
class is posterior-changing and a fit's control carries another record); the design note; the swap law,
run once and recorded in the note: two multiplied forests in either order, and the plain forest first and
second, 16 seeds each way with 400 sweeps discarded and 800 kept, the difference in mean sigma within 2
standard errors and the z of the difference in posterior mean fit, over the rows, with a standard
deviation near 1 (run for this plan: z 0.90 and -0.90 for sigma, 1.11 and 1.04 over the rows, against 0.55
and 0.95 for one order under other seeds).

The pair script. Each row is one model, fitted on the base build as written there and on the slice's as
written here with the same seed: sampler fits compare the train draws, sigma and the amplitudes, `bart`
fits `yhat.train`, sigma and the glue. `identical()` on every row but the last. Run for this plan with a
prototype of the rule: 30 identical, the last three differing.

| | base build | slice |
|---|---|---|
| 01 to 04 | one forest under a named count; plain terms and a term; three `forest()` terms under `cgm(3, 0.8)`; the plain forest written second with its own count and base | the same text |
| 05 to 07 | a list, plain first, on a formula; on a matrix with counts and an `sd`; `bases = list(NULL, w)` under `cgm(2.5, 0.7)` | the same text |
| 08, 09 | probit, a list; logistic, a formula term | the same text |
| 10 to 12 | `bart` with `n.trees = 15`; with counts on both forests, a tree prior and two chains; with `control = dbartsControl(n.trees = 15)` | the same text |
| 13 to 15 | plain first: `interactions` on the forest; at the top; `blocks` at the top | the same text |
| 16, 17 | one listed forest stating its count; both coefficients held, counts and a tree prior | the same text |
| 18 to 20 | all `bases` given, a control naming 25: no list; a list of two `sd`s; a list of two factor-basis forests under 50 | no count on the control, the first forest stating the count, `base = 0.95, power = 2` |
| 21 | a formula of two multiplied terms under a control naming 15 | the three numbers on the first term |
| 22, 23 | probit, a list of two multiplied forests under 25; logistic, three `bases` under 25 | the same respelling |
| 24 | a held factor basis first, under 10 | the same respelling |
| 25 | `bart(..., n.trees = 15, tree.prior = cgm(3, 0.8))`, two multiplied terms | `n.trees = 15, base = 0.8, power = 3` on the first term |
| 26 to 30 | two multiplied forests with `n.trees = 50, base = 0.25, power = 3` stated on the first: a list, a formula, a data object, `bart` under probit, and with `power = 2.5` beside the other two | nothing stated on it (30: `power = 2.5` alone) |
| 31 to 33 | two multiplied forests, nothing stated: a list, a formula, `bart` | the same text: NOT identical, the changed sequence |

Mutations, each expected to fail the named test and no gate before it:

- push 1, forest 1 keeps the fitting function's count when it has a basis: step 1.2's table and its
  draw-for-draw identity;
- push 1, forest 1 keeps the fitting function's tree prior when it has a basis: the same two;
- push 1, the default by position at the data object alone (a first forest whose basis came from `bases`
  is treated as plain): step 1.2's table at that door;
- push 1, the plain forest is taken to be forest 1: step 1.4's refusals on an all-basis model, and step
  1.5's "on forest 2 and not on forest 1";
- push 1, a plain forest past the first takes 50 trees, and in a second mutation `cgm(3, 0.25)`: step
  1.2's plain-second rows;
- push 1, the fitting function's count is read after forest 1's own is written: step 1.2's 9 and 15;
- push 1, each of the five refusals dropped, one mutation each: step 1.4's test of that argument;
- push 1, `bart()`'s record dropped: step 1.3's `bart(n.trees = 75L)`; the control's record ignored: step
  1.3's `dbartsControl(n.trees = 75L)`;
- push 1, both-given checked against forest 1: step 1.5's plain-second refusal;
- push 1, the count is not put back on a carried control: step 1.7; `print` reads the control's count:
  step 1.8;
- push 1, two forests with no basis accepted when the first has one: step 1.1's `list(forest(basis = z),
  forest(), forest())`, added to its refusals;
- push 1, the held refusal decided by class and not by width (two numeric columns held second refused):
  step 1.10's `cbind(1 - z, z)` through the data door, and bartCause's suite;
- push 1, a held forest with no basis refused at every place but the first: step 1.10's third place;
  accepted second: its refusal;
- push 1, a held forest's swap to another width let through: step 1.10's three swaps;
- push 2, a label resolved to the next forest; to its rank among the sorted labels: step 2.2's identity
  for each method;
- push 2, the first of two matching labels taken: step 2.1's two labels of one code, and the label that
  is another position's name;
- push 2, a string of digits coerced to a position: step 2.1's `"2"`, and the list named `"2"`, `"1"`;
- push 2, a factor coerced: step 2.2's `factor("dose")`;
- push 2, `selectForests` matching a label against the margin's names: step 2.3's `"leaf.prior.sd"` by
  label;
- push 2, one method left reading `forest` by the old path: step 2.2's row for that method;
- push 2, the check of `bases` names dropped: step 2.4.

Not a hot-path change: nothing a sweep runs is touched.

## NEWS

No new item: the forests of a model, `forest()` and every `forest` argument are new in 1.0-0 and nothing
released changes. The "Multi-forest models" item is reread and changed only if it says which forest
takes the fitting function's count.

## Out of scope, and where it goes

- `dbarts()`'s own `n.trees`, and the refusal of a count on the control beside one on the plain forest:
  the control-migration arc (dec-B241). The tree prior and the leaf prior given at both places: the
  slices that put `tree.prior` and `leaf.prior` on `forest()` (dec-B246).
- The unit of `sd` with the reader's and `extract`'s lists named by label, the kind by class at the data
  door with the refusal of a basis of three or more levels or of several numeric columns, the multiplier
  law with the held values by kind and the printed block: their slices, in TODO `forest-prior-args`.
- What a held coefficient is where step 1.10 refuses it: the multiplier law's first push, which lifts
  the two interim rows.
- `updateBasisScale` and `leaf.prior = normal(sd = )` on `forest()`: after the merge to main (dec-B276).
- To TODO as new entries: a named `bases` list read by name and not only checked; the names of a data
  object's `bases` as labels; a multinomial fit's categories as labels for `getTrees`; and, if push 3
  leaves it possible, a basis whose text is `forest<digits>` of another position, which this slice
  refuses at selection.

## Calls made in planning

- A `forests` list, and a data object's `bases`, whose plain forest is not first are accepted in their
  own order. dec-B274 says a forest's defaults go by kind "in a formula and in a list alike" and that
  addition commutes; it does not name this list. Refusing it would keep a rule by position, and the
  refusal's own words ("what distinguishes it from the first"). The cost: two shapes with no draw on the
  base build to compare with, so their evidence is the forest-for-forest records, the surfaces of step 1.1
  and the swap law, not a bitwise pair.
- Stated means named, the default's own value included: 75 named and 75 untouched are told apart at every
  door. A wrapper that always names a count, tree prior or leaf prior is therefore refused on a model
  with no plain forest and must leave it out. The alternative, refusing only a value that differs from
  the default, would let `n.trees = 75` beside two multiplied forests mean 50 in silence.
- The leaf prior, `interactions` and `blocks` are refused with no plain forest, though dec-B274 names the
  count and the tree prior: the design's reading of dec-B241 and dec-B246. A leaf prior is refused though
  none changes a model of several forests today, so that giving the fitting function's `leaf.prior` to
  the plain forest later changes no accepted call. With the maintainer (NOTES); planned as written.
- The five refusals come in one block, first, ahead of what any model of several forests refuses, so a
  DART tree prior beside two multiplied forests is told that a tree prior belongs on a forest and not
  that DART is unsupported. Placing them after would leave `blocks` read against forest 1's columns
  first and fail there in other words.
- The count's text has a second form for a control, which says to leave it out of `dbartsControl()`: at
  `dbarts()` the caller wrote no `n.trees` of the fitting function's.
- Each forest's first three numbers say what it runs under, forest 1's too, where the design left forest
  1's at 50, 0.25, 3, unread. Cost: two pins. A record that is true at every position is what the tests
  of step 1.2 read, and the control-migration arc builds its per-forest record from it.
- "On the first forest" becomes "the forest with no basis" in the both-given text, three pins: the first
  forest is no longer the one meant.
- The count given twice at `bart()` is judged on the model (step 1.6). Not in the design; the text-only
  check takes `basis = NULL` for a multiplied forest, a second definition of the plain forest.
- A control never carries a multiplied forest's count to its next fit (step 1.7). Not in the design.
  Found on the prototype: a control from an all-basis fit gave a later one-forest fit 50 trees in
  silence, and the same all-basis model the refusal of a count its caller never wrote. What a control
  from a plain-first fit carries is the tip's and is not changed.
- `print` gives a count per forest (step 1.8). The one line would print a multiplied forest's 50 as the
  fit's. The designed block with labels and sizes stays the multiplier-law slice's.
- A string is never a position. The nine sampler methods take `"2"` as forest 2 at the tip, and
  `TRUE`, a factor and a list by coercion, `factor("dose")` reaching forest 1 where `dose` is forest 2;
  all four are refused. None is in a test, a help page or a consumer (searched).
- `forest<i>` stays a name for position i in every method, beside the label. The design had a string be a
  label only; five test lines and man/bartBT.Rd select `"forest2"`, which under push 3 is the position of
  a forest labelled by its basis. A string that is one forest's label and another position's name is
  refused.
- A label matches exactly, then as code. The alternative, exact only, refuses `"I(dose / 30)"` as the
  caller wrote it.
- `predict(bases = )` checks names and does not reorder by them. Reordering is an addition; reading a
  named list by position in silence, once forests have names users see, is the fault the reviews found
  most often.
- The numbers the design gave, retaken at the tip: 44 creations in 6 files where it counted 39 in 4, 23
  changed assertions with Block F where it counted 6, 1500 lines where it estimated 600.
- Two pushes: the defaults change a model and the selection changes none, and they are reviewed against
  different risks.
- The held shapes the tip gets wrong are refused here, by width and position (step 1.10). The critique of
  the multiplier law found them accepted on every tip until that law's held push, four slices on; the
  coordinator took the finding and placed the refusal in this push, the one that first lets a forest
  with no basis stand second. By width and position and not by class: at this tip the data door cannot
  tell a factor's two indicator columns from two numbers, and a rule by class would refuse bartCause's
  `bcf(update.b = FALSE)`. Narrower than the critique's "every held shape but two": a forest with no
  basis held as the third or a later forest is held at 1, the help's value, and is kept.
- Two numeric columns held as the second forest stay accepted until dec-B282's refusal lands in
  forest-kind-by-class push 2. Refusing them here would need the class at the data door.
- The pair script and the swap law are run at landing and not tracked: one needs the base build, the
  other is a statistical run. What stays in the suite is step 1.2's identity between two spellings the
  slice accepts.
- Rechecked on the landed tip of written-surface (2aabf1ef) before push 1 was built. Every symbol the steps
  name stands and every rule holds; what moved is counts and three texts.
  - The 44 creations that state a count are 44 still, every one a count and nothing else, in six files,
    but test-forest-labels.R has one (a control naming 5) where test-forest-arguments.R had it: that
    file's call sits inside a refusal and creates nothing. 49 held forests of the two shapes kept (17 with
    no basis first, 31 of two indicator columns and one of two numeric columns second) and three swaps
    that change a held forest's width, all in test-forest-arguments.R, as counted.
  - `bart()` takes no `forests` list, so its doors are a formula and a data object.
  - Pair 24 of the pair script, a held factor basis first, is a shape step 1.10 refuses; the pair was run
    with the factor drawn.
- Calls made in building push 1, beyond the steps.
  - The refusal of a held basis of three or more columns says "every column but the first" at the second
    place only and "every column" elsewhere, which is what the engine would hold there.
  - The text for `blocks` gives its own example, `forest(x1 + x2, basis = a, blocks = blocks(list("x1",
    "x2")))`: the one for `interactions` does not fit it.
  - Two tree priors or two leaf priors named at once (`power` and `base`) are refused under the first.
  - `$setControl` takes the control a sampler was created under where the first forest has a basis. The
    stored control holds the first forest's count there, and a control with the caller's own count, or
    none, was refused as changing `n.trees`, which the base build took.
  - `print` reads the counts from the engine while the sampler's pointer is live and from the forests'
    record after a reload.
  - The count given twice at `bart()` is refused for a model of one forest too, where that forest is
    written with `basis = NULL`.
  - No test creates a model with the hold dropped for a basis of three columns: such a basis is refused
    two slices on.
- Rechecked on the landed tip of push 1 (c7fab5bc, 2026-10-07) before push 2 was built, by runs on a
  four-forest sampler (the plain forest second, one list name, two basis texts) and a three-term `bart`
  fit. Every function the steps name stands, and the table of "Selecting a forest" holds as written for
  the nine methods and for `extract` and `predict`. What the steps did not say:
  - The labels live in the forests' record only where a model has several forests: a sampler of one
    forest, one declared as `forests = list(a = forest())`, and a multinomial one record none, so the
    third rule's "no labels" is the record's absence and `forest<i>` needs the forest count from the
    engine. A fit's `attr(fit, "forest.labels")` is there for every fit of several forests.
  - `extract(type = "k")` and `"leaf.prior.sd"` on a fit of one forest are refused as a model parameter
    before any selection is read, and `getTrees` takes its vector through `vapply` over
    `resolveForestIndex`, as `plotTree` takes its one forest; those are the only callers besides the
    nine methods, and a vector of labels needs a branch in `getTrees`.
  - `$setLeafPrior(forests = )` accepts a name only equal to its position's label, so `forest1` on a
    forest labelled `dose` is refused today; the rule of "Lists given by position" accepts `forest<i>`
    at position i, which is built, a loosening of that one call. `predict(bases = )` reads a named list
    by position and ignores the names, as the plan says.
  - Two existing pins change their text: test-forest-labels.R ("must be coercible to type: integer", a
    string given to a sampler) and test-predict-forest.R ("must name one of", a name given to
    `predict`).
  - The suite at the tip is 18412 results over 232 files, as push 1's landing note says.
- The refusal of a stated argument on a model with no plain forest ends, at `bartBT()` and at `bart()`
  given a data object, with "leave it out, or fit with dbarts()" (for the control's count, "leave it out
  of dbartsControl(), or fit with dbarts()"), as the coordinator ruled after the review of push 1, and
  keeps the remedy of "Refused forms" at the formula doors and at `dbarts()` and `dbartsSpec()`. Built as
  its own commit, ahead of push 2; the first two lines of "Refused forms" are the text at the doors that
  keep it.

## Landing note

Push 1 landed 2026-10-07 as 8ae00bdb to 16069daf. Two independent reviews, each told to refute. The first
found the rule sound at every door and two faults at the edges of what a control carries: a count a
caller wrote into a control taken from a sampler was put back to the carried count without a message
(written 20, held 75), and `bartBT` could no longer create a model whose every forest has a basis,
naming a count, tree prior and leaf prior its own caller had not stated. The fix round puts the carried
count back only while the slot still holds the first forest's own count, and has `bartBT` record what its
caller stated, as `bart` does. With every forest holding a basis and nothing stated `bartBT` now gives 50
and 50 trees, the plan's rule, where the build before the slice gave 200 and 50. The second review said
land: 72 fits in sequence were right wherever the caller had not written the first forest's own count,
nine mutations were caught, the three bitwise compares were identical (55, 15, 11) and the four bcf
exact gates passed. Left, and recorded as dec-A179: a count written into a carried control that equals
the count it already shows cannot be told from no edit; a control taken from a `bart()` or `bartBT()` fit
and given to `dbarts()` for a model whose every forest has a basis, and `setControl` given an all-basis
sampler's control on another sampler, are refused where the build before took them, the refusal naming
a count nobody stated (root TODO, control-carried-count-edges, for the control migration); the refusal's
remedy, to state the count on a forest, cannot be written at `bartBT` or at `bart()` given a data object
(push 2 rewords it). Gates on a clean copy of the rebased head: install, tests/cpp 357 lines passing, the
suite 18412 results and none failed over 232 files, lint, format, the document checks, the mutation
anchors, build, and check with the one standing NOTE; bartCause's suite on the landed build is in the
records commit's message.
