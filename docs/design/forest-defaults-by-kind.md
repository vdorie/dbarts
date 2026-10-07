# A forest's defaults go by its kind

Status: PROPOSED 2026-10-07 (dec-B274); the defaults are built, the selection of a forest by its label is
not. Plan: [forest-defaults-by-kind.md](../plans/forest-defaults-by-kind.md). Rulings: dec-B274, dec-B241 and
dec-B246 as dec-B274 restates them, dec-B281 and dec-B282 for the held shapes, in docs/decisions.md.

## The rule

A forest has a basis or it has none. A model of several forests has at most one forest with no basis, the
plain forest, and it may stand at any place in a `forests` list or in a data object's `bases`; in a formula
it is the first forest wherever it is written. The kind is asked in one place,
[`plainForest`](../../R/model.R), of the bases the model ends with.

| | the plain forest | a forest with a basis |
|---|---|---|
| tree count | the fitting function's: `bart`'s `n.trees`, the control's at `dbarts` and `dbartsSpec`; 75 | 50 |
| tree prior | the fitting function's `tree.prior`; `base = 0.95`, `power = 2` | `base = 0.25`, `power = 3` |
| `interactions`, `blocks` | the fitting function's or its own; both given is refused | its own |

A forest's own `n.trees`, `base` and `power` each govern that forest, and each one left out takes the
default of the forest's kind. Position decides nothing: two forests of one kind written in the other order
are the same model, and other draws, since forests are swept in order.

Where no forest is plain, what the fitting function states for the plain forest has no forest to belong
to. Each of the tree count, the tree prior, the leaf prior, `interactions` and `blocks` that is stated is
refused by name, in that order, before anything else is read from it, and the model is created once it is
moved to a forest or left out ([`refuseStatedWithNoPlainForest`](../../R/model.R)).

## What counts as stated

Stated means named, the default's own value included: `n.trees = 75` beside two multiplied forests is
refused, so that it cannot come to mean 50 in silence ([`plainForestStated`](../../R/model.R)).

| door | the count | the tree prior | the leaf prior |
|---|---|---|---|
| `bart()` | `n.trees` named, or a `control` that states one | `tree.prior`, or the retired `power`, `base` or `split.probs` | `leaf.prior` or `k` |
| `dbarts()`, `dbartsSpec()` | a control that states one | `tree.prior` | `leaf.prior`, or the retired `node.prior` |
| `bartBT()` | `ntree` | `power`, `base` or `splitprobs` | `k` |

A control states a count when its constructor's call named `n.trees` or its slot is not the constructor's
default. A slot edited back to 75 and a control made by `new()` read as not stated: that is the limit of
what can be told. `bart()` and `bartBT()` each build their own control and always hand `dbarts()` a tree
prior, so each records what its own caller stated, under its own names, on the control it hands over, and
the record is taken off when the sampler's specification is resolved. `bartBT()` reaches a model of
several forests through a data object's `bases` alone.

`bart()` refuses its own `n.trees` beside one on the plain forest, judged on the model, so a forest
written with `basis = NULL` is the plain forest too. At `dbarts()` the forest's own count governs over the
control's.

## Where the numbers live

The bridge reads the first forest's tree count from the control and its tree prior from the model, so
those two hold the first forest's, whatever its kind: 50 and `cgm(3, 0.25)` for a first forest with a basis
that states none. Each forest's record starts with the count, base and power that forest runs under, the
first forest's included ([`forestParams`](../../R/model.R)).

A control taken from a fit and given to another never carries a multiplied forest's count. Where the first
forest has a basis the forests' record keeps the count the next fit is to inherit, the plain forest's as
it ran or the fitting function's own where no forest was plain, and it goes back on the control before
that fit reads it. It goes back only while the control still holds the first forest's count: any other
count there is the caller's edit and stands as a stated count does, the plain forest's in the next fit
and refused where every forest of it has a basis. An edit to the first forest's own count cannot be told
from none. `$setControl` takes the control a sampler was created under as well as the sampler's
own. `print` of a fit gives a count for each forest, `n.trees: 50, 75`.

## What changes

One sequence of calls: a model in which every forest has a basis, with no tree count and no tree prior
given to the fitting function, whose first forest leaves any of the three unstated. That forest ran under
75, 0.95 and 2 and now runs under 50, 0.25 and 3; at `bartBT()`, whose count is 200, under 200. The same
call with the three numbers written on that forest, fitted on the build before this one, gives identical
draws: six such pairs, across a list, a formula, a data object, `bart` under probit and `dbartsSpec` under
logistic. Every fit of one forest, and every fit of several whose plain forest is first, draws what it
drew: 33 pairs identical. Twelve calls that state a count, a tree prior, `interactions` or `blocks` beside
forests that all have a basis are refused now and, respelled onto the first forest, draw what they drew.

The same model in the other order has the same law. Sixteen seeds each way, 400 sweeps discarded and 800
kept, 150 rows: two multiplied forests in either order differ in mean sigma by 0.57 standard errors, and
the z of the difference in posterior mean fit has standard deviation 0.96 over the rows; the plain forest
first against second, -0.25 and 1.07. One order under other seeds gives -0.84 and 1.17, and 1.45 and 0.94.

## A held coefficient, for now

The engine holds a coefficient where it starts it, by the forest's position: the second forest's first
column at 0 and its others at 1, every other forest's at 1. Until a held coefficient has a value of its
own for every shape, `amplitude = fixed()` is taken only where that is the value the help states
([`refuseHeldShape`](../../R/model.R)).

| held forest | position | held at | taken | lasts |
|---|---|---|---|---|
| no basis | first; third or later | 1 | yes | |
| no basis | second | 0, so the forest would leave the model | refused | interim |
| a basis of two columns | second | (0, 1) | yes | two numeric columns are refused later with the basis itself (dec-B282) |
| a basis of two columns | any other | (1, 1), no contrast | refused | interim for a factor of two levels |
| a basis of three or more columns | any | (0, 1, 1, ...) or all 1 | refused | for good (dec-B281, dec-B282) |
| a basis of one column | any | | refused | interim (dec-A171) |

`$setForestBasis` refuses a held forest a replacement of another width, a forest created with no basis
taking none. Nothing is held anew and no held value changes.
