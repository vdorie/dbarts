# A forest's defaults go by its kind

Status: LANDED 2026-10-07 (dec-B274, dec-B276): the defaults as 8ae00bdb to 16069daf, the selection of a forest
by its label as 13829fb9 to 3ef31f97. Plan: [forest-defaults-by-kind.md](../plans/forest-defaults-by-kind.md). Rulings: dec-B274, dec-B241 and
dec-B246 as dec-B274 restates them, dec-B281 and dec-B282 for the held shapes, in docs/decisions.md.

## The rule

A forest has a basis or it has none. A forest with no basis, a plain forest, may stand at any place in a
`forests` list or in a data object's `bases`; in a formula the plain forests come first, in the order
written. One that states no sd takes the fitting function's leaf prior as the one forest of a fit does, held at
1, whatever stands beside it (dec-B509). A model may hold several: two or more that state no sd share that leaf
prior as one total and split its tree count ([forests-shared-total.md](forests-shared-total.md)). The kind is asked in
one place, [`plainForest`](../../R/model.R), of the bases the model ends with.

| | the plain forest | a forest with a basis |
|---|---|---|
| tree count | the fitting function's `n.trees`, at `bart`, `dbarts` and `dbartsSpec`; 75 | 50 |
| tree prior | the fitting function's `power` and `base`; `base = 0.95`, `power = 2` | `base = 0.25`, `power = 3` |
| `interactions`, `blocks` | the fitting function's or its own; both given is refused; beside a second plain forest, its own | its own |
| leaf prior and size | the fitting function's leaf prior, held; `sd` states its own, held; `normal()` or `student()` draws | its `sd`, drawn |

A forest's own `n.trees`, `base` and `power` each govern that forest, and each one left out takes the
default of the forest's kind. Position decides nothing: two forests of one kind written in the other order
are the same model, and other draws, since forests are swept in order.

Where no forest is plain, what the fitting function states for the plain forest has no forest to belong
to. Each of the tree count, the tree prior (`power` and `base`, named one by one), the leaf prior, `interactions` and `blocks` that is stated is
refused by name, in that order, before anything else is read from it, and the model is created once it is
moved to a forest or left out ([`refuseStatedWithNoPlainForest`](../../R/model.R)).

## What counts as stated

Stated means named with a value, the default's own value included: `n.trees = 75` beside two multiplied
forests is refused, so that it cannot come to mean 50 in silence ([`plainForestStated`](../../R/model.R)).

A written `NULL` follows the argument's default in the signature. Where the default is `NULL` (`k` and
`leaf.prior` on `bart()` and `xbart()`, every argument of `forest()` and `varianceForest()`) it is the name
left out and states nothing: `bart(y ~ forest(x1 + x2, sd = 2), k = NULL)` fits as it does without `k`, and
`forest(x1 + x2, power = NULL)` beside `power = 3` on the fitting function is no collision. Anywhere else
(`n.trees`, `power` and `base` on every door, `k` on `rbart_vi()` and `bartBT()`, `leaf.prior` on `dbarts()`
and `dbartsSpec()`) it is refused, the message naming the argument and saying to leave it out. Names are judged
by the value they hold, not by the code that wrote them: a wrapper's `power = power` is read as what `power`
holds in the wrapper, and a `NULL` passed through `...` is the same as one written in place ([`nullWritten`](../../R/model.R),
[`refuseNullArgument`](../../R/model.R)). `split.probs`, `sigest`, `seed` and the other arguments that default
to `NULL` are absent when `NULL`, as before.

| door | the count | the tree prior | the leaf prior |
|---|---|---|---|
| `bart()` | `n.trees` named, or a `control` that writes one | `power` or `base` | `leaf.prior` or `k` |
| `dbarts()`, `dbartsSpec()` | `n.trees` named, or a `control` that writes one | `power` or `base`, or the retired `tree.prior` | `leaf.prior`, or the retired `node.prior` |
| `bartBT()` | `ntree` | `power` or `base` | `k` |

`split.probs` is the plain forest's alone and is not among the stated names: a forest takes none, so it
reaches no forest that could refuse it by name. Beside forests that all have a basis a uniform one states the
default and is taken, and one that is not uniform, or `dart()`, meets the several-forest refusal. `bartBT()`
forwards `power` and `base` always, and `splitprobs` only where its caller named it.

Two levels, one rule. A `power` or `base` named on the fitting function and on the plain forest is one
setting given twice and is refused by name, as `n.trees` is; named on only one of the two it governs the
plain forest, the forest's own value winning where the fitting function's is the default. A variance forest
that states none reads the plain forest's resolved values, and never rescues a value: a lone mean forest
stating `power` beside the fitting function's `power` is refused whatever the variance forest states. Where
no forest is plain the refusal names the argument given and the call to write, as
`forest(x1, basis = a, power = 3)`.

A control writes a count only when its constructor's call named `n.trees` or its slot differs from the one the
sampler left on it ([`controlValueKind`](../../R/records.R)); a slot edited back to the sampler's reads as
carried. `bart()` and `bartBT()` build their own control and pass the count flat, so each hands `dbarts()` a
control naming none and a tree prior, and what the caller stated rides the specification until it is resolved.
`bartBT()` reaches a model of several forests through a data object's `bases` alone. A count named flat beats a
written one, which waits with the carried ones for the not-used warning.

`bart()` refuses its own `n.trees` beside one on the plain forest, judged on the model, so a forest
written with `basis = NULL` is the plain forest too. At `dbarts()` a count written on the control beside
one a forest states is refused.

## Where the numbers live

Each forest's record holds the count, base and power that forest runs under, the first forest's included
([`forestParams`](../../R/model.R)), and the bridge reads them from there; the model's tree prior holds the first
forest's, so a first forest with a basis that states none holds 50 and `base = 0.25, power = 3`. The control
holds none of it: its `n.trees` is a copy the sampler writes, the count at the fitting function's level
([`fittingLevelCount`](../../R/records.R)), the forest with no basis's or else the first forest's, by
[`mirrorControl`](../../R/records.R).

A control taken from a fit and given to another sets no count: the copy is read as carried, not written
([`controlValueKind`](../../R/records.R)), so the next fit takes the count it is given and the forests' own,
and a warning says when the control's differs. A count a caller writes on a control is a stated one, the
plain forest's, and is refused where every forest has a basis or the forest states its own.
`$setControl` keeps the sampler's count. `print` of a fit gives a count for each forest,
`n.trees: 50, 75`.

## What changes

One sequence of calls: a model in which every forest has a basis, with no tree count and no tree prior
given to the fitting function, whose first forest leaves any of the three unstated. That forest ran under
75, 0.95 and 2 and now runs under 50, 0.25 and 3; at `bartBT()`, whose count is 200, under 200. The same
call with the three numbers written on that forest, fitted on the build before this one, gives identical
draws: six such pairs, across a list, a formula, a data object, `bart` under probit and `dbartsSpec` under
logistic. Every fit of one forest, and every fit of several whose plain forest is first, draws what it
drew: 33 pairs identical. Twelve calls that state a count, a tree prior, `interactions` or `blocks` beside
forests that all have a basis are refused now and, respelled onto the first forest, draw what they drew.

The same model in the other order has the same distribution. Sixteen seeds each way, 400 sweeps discarded
and 800 kept, 150 rows: two multiplied forests in either order differ in mean sigma by 0.57 standard
errors, and the z of the difference in posterior mean fit has standard deviation 0.96 over the rows; the
plain forest first against second, -0.25 and 1.07. One order under other seeds gives -0.84 and 1.17, and
1.45 and 0.94.

## A held coefficient, for now

The engine holds a coefficient where it starts it, by the forest's position: the second forest's first
column at 0 and its others at 1, every other forest's at 1. Until a held coefficient has a value of its
own for every shape, `amplitude = fixed()` is taken only where that is the value the help states
(retired: [`refuseHeldShape`](../../R/model.R), gone with forest-multiplier's push 2, when every
shape held has its value).

| held forest | position | held at | taken | lasts |
|---|---|---|---|---|
| no basis | first; third or later | 1 | yes | |
| no basis, the only one | second | 0, so the forest would leave the model | refused | interim |
| no basis, sharing a total | any | 1 | yes | |
| a basis of two columns | second | (0, 1) | yes | two numeric columns are refused later with the basis itself (dec-B282) |
| a basis of two columns | any other | (1, 1), no contrast | refused | interim for a factor of two levels |
| a basis of three or more columns | any | (0, 1, 1, ...) or all 1 | refused | for good (dec-B281, dec-B282) |
| a basis of one column | any | | refused | interim (dec-A171) |

`$setForestBasis` refuses a held forest a replacement of another width, a forest created with no basis
taking none. Nothing is held anew and no held value changes. Since forest-kind-by-class's second push a basis
is one numeric column or two levels, refused otherwise at creation before a hold is read, so the rows for two
numeric columns and for three or more are that refusal's, and no forest changes width, held or drawn
([forest-kind-by-class.md](forest-kind-by-class.md)).

## Selecting a forest

Every place that takes a forest takes its label as well as its position: the sampler's nine methods
(`getLeafPrior`, `getK`, `getForestFits`, `getForestAmplitudes`, `getForestVariableCounts`, `getTrees`,
`plotTree`, `setForestWeights`, `setForestBasis`), `extract` with `type = "forest"`, `"k"`, `"leaf.prior.sd"`
and `"trees"`, and `predict(type = "forest")`. One function reads the argument,
[`selectForest`](../../R/bartcore.R), and both the sampler and a fit go through it.

- A number is a position, as base R reads a list element by number, and goes to the code that read it before, with its texts.
- A string is never a position. It is the forest with that label; failing that, the forest whose label is
  the same code, so `"I(dose / 30)"` finds `I(dose/30)`; and `forest<i>` is forest i, the name position i
  has on every per-forest margin, so `"forest2"` selects what it did.
- A string that is one forest's label and another position's name, or the same code as two labels, is
  refused naming both. `"2"` is the forest labeled `"2"` or no forest, and the refusal says a position is
  given as a number.
- `NA`, an empty string, a factor, a logical and a list are refused by name. `TRUE`, a factor and a list
  were read as numbers before.
- A sampler of one forest, one declared with a single `forests` entry, and a multinomial one record no
  labels: `forest<i>` is taken where the number is, and any other string is refused.
- The margins of a fit keep their names, `forest1`, `forest2`; a label is mapped to its position and never
  matched against them. What comes back is what the position gives, to the bit, and nothing is renamed by
  label. `extract(type = "trees")` and `plotTree` hand their forest to the sampler.

A list given forest by forest is read by position. A name on an entry of `predict(bases = )` or of
`$setLeafPrior(forests = )` must be that position's label or `forest<i>`, and any other is refused, so that a
list in another order is not read in silence as the order of its positions.
