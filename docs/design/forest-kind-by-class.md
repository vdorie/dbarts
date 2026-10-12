# What a forest multiplies is decided by the class of its basis

Status: PUSH 1 BUILT on wt/forest-kind-recheck (e3aac0ea), not landed; push 2 (one size a forest) is planned. The
plan is [forest-kind-by-class.md](../plans/forest-kind-by-class.md), which carries the steps, the gates and the
calls made in planning; this note is the standing statement of the rule push 1 builds. Rulings: dec-B260 and
dec-B254 for the kind as structure, dec-B289 and dec-B308 for a factor's levels, dec-B346 for a column of zeros,
dec-B452 for a swap that never changes the kind, dec-B453 for numbers refused for a factor forest at predict,
dec-B231 and dec-B418 for a state that carries no data, in docs/decisions.md.

A forest multiplies levels or numbers, and the class of the basis handed over says which, at every door that
takes one: a `forest()` term, a `forests` list, `dbartsSpec()`, `dbartsData(bases = )`, `$setForestBasis` and
`predict(bases = )`. The content is never looked at.

## 1. The rule

1. Class. A factor, ordered or not, a character vector and a logical vector are levels, with a coefficient for
   each: the factor's own order, the sorted values of a character vector, `FALSE` then `TRUE`. An integer or a
   double, a vector or a matrix, is numbers, whatever its values: a 0/1 number is one column and one
   coefficient, and `factor(z)` is how two levels are written. A matrix of logicals or of characters, of one
   column or several, a data frame and any other class are refused by name. One function answers
   ([`basisClass`](../../R/forestBasis.R)); it refuses nothing, and the reader that called it words the refusal.
2. Rows. A basis is judged on the rows the fit keeps, after `subset` and the na.action, in every spelling:
   a level no kept row has gives no column, with no message; fewer than two left is refused; a missing value on
   a kept row is refused and on a dropped row not looked at; a factor with NA as a level and a kept row at it is
   refused, and with none is dropped as any level no row has; a numeric column of zeros on every kept row is
   refused. A row of weight zero or masked out is a row the fit keeps. This is done once, at a creation, by
   [`finishBasisRows`](../../R/forestBasis.R), which `dbartsData()` and `dbartsSpec()` call.
3. The record. A factor's levels ride the data object, one entry for each forest, in the slot `basis.levels`
   (`NULL` for a forest of numbers or with none, and for a data object with no bases). The block in `bases`
   stays the engine's form, a factor's indicator columns among them. A forest's kind is read from the block
   and the levels by [`forestBasisKind`](../../R/forestBasis.R) and from nowhere else. A plain function
   ([`checkBasisRecord`](../../R/data.R)) holds that a forest's recorded levels are as many as its block has
   columns; the sampler's `initialize`, `dbarts()` and `dbartsSpec()` call it, so `copy()` does too, and a
   reload's re-creation does not.
4. After creation ([`conformBasis`](../../R/forestBasis.R)), by one function that every spelling reaches:

   | forest created as | takes | refuses |
   |---|---|---|
   | none | nothing | any basis |
   | levels | a factor, a character vector or a logical vector, matched to the forest's levels by name; a level of the forest no row has is a column of zeros | numbers, the sampler's own stored block among them; a row at a level the forest lacks, or at NA |
   | numbers | a numeric vector or a matrix of one column, a column of zeros among them | a factor, a character vector, a logical vector |

   A level the value declares and no row has is passed by, so the factor a forest was created from is taken
   again. A name matches character for character. A levels forest has the same columns for life: coefficient j
   multiplies level j's rows after any swap, and the recorded levels are creation's whatever a swap declares.
   A numeric forest that is drawn may still take another width until push 2.
5. What a fit stores. A fit's `bases` holds, for a forest of levels, the factor, one value for each row the fit
   kept, its levels those recorded at creation ([`storedForestBases`](../../R/forestBasis.R)); for a numeric
   forest the block; `NULL` for none. Each entry is what `predict(bases = )` takes for its forest, so
   `predict(fit, x, bases = fit$bases)` at the training rows is the training fit. The prediction turns the
   factor back to its block with [`storedBasisBlock`](../../R/forestBasis.R).
6. What the sampler keeps. The block and the levels, not the factor: a caller that will hand a basis back keeps
   the factor it handed over, and reads the forest's levels at `data@basis.levels[[f]]`.
7. State. A state carries no data, and the levels are data. Installed in a sampler whose factor has the same
   levels in another order, it goes by position, with no message; the recipient's levels and block are its own.
8. `$setData` on a sampler of one forest refuses a data object that carries bases, before anything is stored
   ([`bartcoreSamplerSetData`](../../R/bartcore.R)); the sampler has no forest to multiply them.
9. A saved sampler. A data object saved without the slot has no levels
   ([`dataBasisLevels`](../../R/data.R) answers `NULL`), so its forests read as multiplying numbers and its
   record passes the check; a saved factor forest then takes numbers at a swap and refuses a factor.

## 2. The changed sequences and their oracles

Four call sequences fit another model than before and appear in no test, benchmark, help example or consumer:

1. A logical vector handed to `dbartsData(bases = )`: one 0/1 numeric column before, two levels now. Oracle: the
   base build's fit of it through a `forests` list, identical draws.
2. `$setForestBasis` on a factor forest, given a factor whose levels stand in another order or lack one: matched
   by position before, by name now. Oracle: the base build given the same vector as
   `factor(v, levels = <the levels at creation>)`.
3. `predict(bases = )` with such a vector: the same, with the same oracle.
4. A factor of three or more levels handed over as a value, where `subset` or the na.action leaves one level
   with no row and two or more with one: a column of zeros with a coefficient of its own before, the levels that
   have a row now. Oracle: the base build's fit of the same factor written as code, identical draws.

A fit that both builds accept, these four apart, is the same to the bit, which the equivalence harnesses and the
reference-build snapshot files check.

## 3. What it leaves

Push 2 refuses, at every door, a factor, character or logical basis of three or more levels and a basis of
several numeric columns, with the form to write, and respells the tests' blocks of indicator columns as the
factors they stand for. The multiplier slice gives a factor and a number different priors. The data object's
column naming at the data door is as it was.
