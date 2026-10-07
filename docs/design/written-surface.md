# A model of several forests is written one way, forest by forest

Status: BUILT; pushes 1 and 2 LANDED 2026-10-06 and 2026-10-07, push 3 not yet landed. Plan:
[written-surface.md](../plans/written-surface.md). Rulings: dec-B266 to dec-B275, dec-A171 and dec-A173 in
docs/decisions.md.

## What was wrong

A model of several forests could be written in two ways that did not agree. In a formula a multiplier was a
colon, `z:forest(x1)`, and the first forest could not be written at all. In a `forests` list a multiplier was
`basis = ~ z`. A basis was R code, so `~ dose + age` was one column, the sum, two columns needed `cbind()`, and
the package carried code of its own to rebuild `scale()` and `poly()` at new rows. The two doors read the same
basis on different rows under `subset`, found its names in different places, and only one of them could
predict. A held coefficient was a flag, `update.amplitude`.

## The grammar

    forest(vars = NULL, basis = NULL, sd = NULL, n.trees = NULL, base = NULL, power = NULL,
           amplitude = NULL, interactions = NULL, blocks = NULL, amplitude.prior.variance = NULL)

Only the first argument may be given without its name, so that any later argument is an addition.

| written | means |
|---|---|
| `y ~ forest(x1 + x2) + forest(x1, basis = dose)` | two forests; the second is multiplied by `dose` |
| `y ~ x1 + x2 + forest(x1, basis = dose)` | the same model: the plain terms are the forest with no basis |
| `forests = list(forest(), forest(x1, basis = dose))` | the same model, in a list |
| `basis = dose + age` | two columns, a coefficient each |
| `basis = I(dose + age)` | one column, the sum |
| `basis = factor(z)`, a character or a logical column | one column for each level that a kept row has |
| `basis = scale(age)`, `poly(dose, 2)`, `log(dose)` | a term as in `lm`, rebuilt at new rows as `lm` rebuilds it |
| `basis = dose:age` | one column, the product |
| `basis = 1 + dose`; `0 + dose`, `dose - 1` | a column of ones and `dose`; `dose` |
| `basis = ~ dose + age` | the same as without the tilde, in every respect |
| `b <- ~ dose + age` and then `basis = b` | the same columns; read late, as R reads a formula |
| `amplitude = fixed()` | the forest's coefficient is held, not drawn |

The forests of a formula, a forest's predictors and the two arguments `sd` and `amplitude` are pushes 1 and 2
and are described in the plan. This note is about the basis, which is push 3.

## A basis is the right-hand side of a model formula

A basis written as code is read by [`basisCode`](../../R/forestBasis.R) and built by R's own model frame and
model matrix on `~ 0 + <basis>`. So it means what it means in `lm`, at the fit and at new rows, and the package
no longer carries a reading of its own.

Before anything is evaluated, [`parseBasisGrammar`](../../R/forestBasis.R) walks the top of the code. `+`
separates terms, parentheses group, `:` between terms is their product, and `1`, `0 +` and `- 1` are read as in
`lm`. Everything else at the top is a term and is left to R. Refused there, by name and with the form to write:

- `*` between terms (dec-B273). In `lm` it is both columns and their product, and before this change it was the
  product alone, so either reading would have surprised someone.
- A term multiplied or divided by a number, as `dose / 30`. A size for each column may later be written this way
  (dec-B272), so the form is kept free. `I(dose / 30)` is accepted and is the rescaled column.
- `-` other than `- 1`, `^`, `/`, `%in%`, `|`, `.` and `offset()`: in a model formula they remove, cross or nest
  terms, and read as arithmetic they would be another basis.
- `cbind()`: `+` is the one way to write several columns. The refusal says what to write, and says nothing where
  no sum of terms gives the same columns in the same order.
- A term that calls `normal`, `fixed`, `student`, `cauchy`, `linear`, `gp`, `cgm`, `dart`, `chisq`, `chi`,
  `invchi`, `forest` or `varianceForest`, whoever defines the function: a prior on a term may later be written
  so (dec-B272).

A basis is one kind: a single factor, character column or logical column, which gives a column for every level
with no level dropped, or numeric terms ([`basisColumns`](../../R/forestBasis.R)). A factor beside other terms
is refused, the engine having one law for a forest's coefficients.

## One reading at both doors

A `forest()` term of a formula and a `forest()` in a `forests` list reduce a basis to the same thing, code with
the environment its names are looked up in, and from there one function reads it
([`readForestBasis`](../../R/forestBasis.R)) and one builds it ([`buildCodeBasis`](../../R/forestBasis.R)). The
two doors cannot disagree about a basis because neither has code of its own for one.

Rows are decided once, and in one place. The data object evaluates `subset` once, where it is built, and
applies the na.action; nothing else reads the caller's `subset`. Every basis rides the data object's `bases`
argument and is cut there. A value rides as itself. For a basis written as code the numbers of the data's rows
ride in its place ([`basisRowNumbers`](../../R/forestBasis.R)) and come back as the rows the fit kept, in the
fit's order ([`buildFitBases`](../../R/forestBasis.R)). So a subset that draws its rows, `subset = sample(n,
100)`, gives every basis the fit's own rows; read a second time it gave them another draw's, without a message,
which the first review of this push found at every door.

A basis written as code is read over every row of the data the fit was given, and only on the rows kept is it
decided which levels have a row, whether a value is missing and whether a column is all zeros. What a term
computes across rows, `scale()`, `poly()` and `I(age - mean(age))` alike, comes from every row, as in `lm`. A
level that the kept rows leave empty is no column, whatever emptied it.

A basis handed over as a value has no code, and is cut by the data object, which refuses a value that was
already cut. A value is not built again on the kept rows: a level that `subset` empties keeps a column of
zeros, where code drops it. Giving a value and code one builder is the kind-by-class slice's.

At new rows the record a fit keeps is R's `terms` object with its environment, and
[`replayForestBasis`](../../R/model.R) rebuilds the basis through `model.frame()`. This works whichever door
wrote the basis; a basis written in a list could not be rebuilt before. A column of the fit's data must be among
the new rows. A number found where the basis was written is used. A vector with a value for every fitted row is
refused, since its values belong to the fitted rows.

## Where a name is found

It depends on who wrote the code, and the rule is dec-A173's with the rulings made on the review of this push.

A call of `forest()` is made by the caller, in a `forests` list or ahead of the fit
([`captureForestBasis`](../../R/forestBasis.R)):

- A name of a column of `data` is that column. A column hides a caller's variable of the same name.
- Anything else is what it was where `forest()` was called, at that moment. The value is taken once, quietly;
  what it warns of is kept and raised by a fit that uses it.
- No variable of the caller's is looked up later than the call. Every name the code uses that the caller binds,
  in the frames from the call's own up to the workspace, is copied when `forest()` is called
  ([`bindBasisAtCall`](../../R/forestBasis.R)), a number written beside a column included, so
  `for (k in c(10, 30)) forest(basis = I(dose / k))` gives each forest its own `k`, at the fit and at
  `predict`. What R or an attached package supplies, `scale` or `pi`, is not copied and is looked up at the
  fit.
- A tilde written in place changes none of this: `forest(basis = ~ I(dose / k))` is read exactly as the same
  code without the tilde.

A formula made elsewhere is R's own object and keeps R's own reading: a formula held in a variable or handed
over, a `forest()` term inside a fit's formula, and the formula given to `$setForestBasis` are read against the
data and then in the formula's environment, when the model is built and again by `predict`. A formula made in a
loop and held in a variable therefore sees the loop's last value in every forest.

Two cases are not the caller's code at all. An argument that a function of base R writes, `lapply()`'s
`X[[i]]` and `Map()`'s `dots[[2L]][[1L]]`, is a value handed over: its text is no label and its names are not
looked up in the data ([`isLoopMachinery`](../../R/model.R) tells it by the frame the argument stands in). And a
basis forwarded through dots by a call that has returned cannot be read where it was written: its value is used
when its code names no column of the data, and it is refused by name when it names one
([`forwardedBasis`](../../R/forestBasis.R)).

So a forest built in a loop, by `lapply()` or by `Map()`, a variable changed or removed before the fit, and a
forest saved and read back each fit the basis the call was given. A forest carries the names its basis uses and
no frame of the function that built it, except through a held formula, which carries its environment as any
formula does.

## Names and labels

The columns of a basis are named as `coef(lm(y ~ 0 + <basis>))` names them: `dose`, `age`; `poly(dose, 2)1`,
`poly(dose, 2)2`; `factor(z)0`, `factor(z)1`; `(Intercept)`. A value keeps its names when every column has one
and no two are alike, and otherwise has none. Two columns alike are refused. `$setForestBasis` takes columns by
position and the recorded names stay; a replacement with the recorded names in another order is refused
([`refuseReorderedBasisNames`](../../R/forestBasis.R)).

Every forest has a label, fixed at creation ([`forestLabels`](../../R/forestBasis.R)): the list name where one
is written, else the text of the basis, else `forest<i>`. Labels are unique. They are recorded and checked and
select nothing yet; selecting a forest by its label is a later slice.

Both sets of names are fixed now because the per-column reader and `extract` will be named by them for good
(dec-B272, dec-B275).

## The texts whose meaning changed

Each was new in 1.0-0 and none had been released. Each is now a model that another text fitted before, so that
text's draws are the oracle. Measured with seeded fits, the old text on the build before this push and the new
text on this one.

| text | before | after | identical to, on the build before |
|---|---|---|---|
| `basis = ~ dose + age` | one column, the sum | two columns | `~ cbind(dose, age)` |
| `basis = ~ dose - 1` | `dose` minus one | `dose` | `~ dose` |
| `basis = ~ 1 + dose` | `dose` plus one | a column of ones and `dose` | `~ cbind(1, dose)` |
| `basis = dose` in a term, beside a caller's `dose` | the caller's | the data's | `~ dose` |
| `basis = dose` in a list, beside a caller's `dose` | the caller's | the data's | `~ dose` |
| a term computed across rows (`scale()`, `poly()`, `I(age - mean(age))`, `cut(age, 3)`) in a formula's term under `subset` | computed on the rows kept | on every row | the constants written out, to 1e-10 |
| `~ factor(g)` in a list, a level emptied by `subset` | a column of zeros for that level | no column | the same text in a formula's term |

The first three rows hold with no tilde too, for a basis written beside the caller's own vectors, which the
build before took as a value: `basis = dv + av` was one column, the sum, and is two; `dv - 1` and `1 + dv`
likewise. The last row was found while building. Before, the two doors differed: the term dropped the level and
the list kept an all-zero column, whose coefficient the data never touched. The old meaning of each of the
first three is still written, inside `I()`: `I(dose + age)`, `I(dose - 1)`, `I(1 + dose)`, each identical to the
old text.

A model whose text is in none of these families, and whose `subset` is the same each time it is read, draws
what it drew: 218 comparisons of seeded fits over three families, one and two chains, both doors, held
and drawn coefficients, factor, two-column and value bases, `subset`, weights, an offset and a test set. A
`subset` that draws its rows is the exception by design: the build before gave its bases the rows of a second
draw.

## What is fixed for later

dec-B272 defers how a size for each column is stated, on the condition that adding it changes no call. What
this slice does for that:

- One sd on a basis of several columns states that size for each column.
- A longer or named sd, a term over or times a number, and a prior's name on a term are refused by name.
- The names of a basis's columns and of the forests follow a complete rule.
- A star in a basis is refused (dec-B273).

## Left open

A missing value in a basis is refused where `lm` would drop the row. A hazard fit of several forests does not
check that a basis covers the data under `subset`, its rows being expanded. Defaults by a forest's kind, the
unit of `sd`, the kind of a value by its class, the per-column reader and selecting a forest by label are later
slices (TODO `forest-prior-args`).
