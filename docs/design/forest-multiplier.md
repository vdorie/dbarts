# What a forest's coefficient is, held or drawn, and what a numeric multiplier's sd is the size of

Status: LANDED: pushes 1 and 2 2026-10-17 (5a62933e, 5a1e9dca), push 3 2026-10-23 (4c65a3c0). The plan is
[forest-multiplier.md](../plans/forest-multiplier.md), which carries the steps, the gates and the calls.
This note states the rule the first two pushes build: what a forest's coefficient is, held or drawn,
the default of each kind of forest, a numeric column's scale and its record, the bound under which a row
leaves a forest's update, and the unit on a response with no spread. Rulings: dec-B290, dec-B291 and
dec-B293 for the defaults, dec-B303 with dec-B302 and dec-B309 for the rows a column's scale reads,
dec-B307 for a constant column, dec-B353 and dec-B354 for a held coefficient, dec-B376 for what a fit
records of a hold, dec-B346 and dec-B452 for a swap and a kind that no swap changes, in
docs/decisions.md.

## 1. The rule

A forest of a model of several multiplies nothing, the two levels of a factor, or a numeric column,
by the class of its basis ([forest-kind-by-class.md](forest-kind-by-class.md)). Its term at row i is
the basis row times the forest's coefficients, times the forest's own fit F. s is the forest's `sd`
in the response's units; U is the unit ([forest-sd-unit.md](forest-sd-unit.md)).

| the forest multiplies | drawn | held, `amplitude = fixed()` |
|---|---|---|
| nothing | F has sd U, the coefficient's size half-Cauchy with median s / U | the coefficient is 1; F has sd s, exactly |
| two levels | each level's coefficient N(0, 1/2); F has sd s / 0.674 | (0, 1); F has sd s, exactly, the second level against the first |
| a numeric column w | the coefficient N(0, 1), starting at 1; F has sd s / 0.674, per unit of w | the coefficient is 1; F has sd s, exactly, per unit of w |

- A held coefficient is 1, and what it multiplies is the basis (dec-B354): 1 with no basis, 0 and 1
  for a factor's two levels, the column as given for a number. It goes by the forest's kind and never
  by where the forest stands. `fixed()` takes no value but 1; another coding is another basis, as
  `basis = z - 0.5`, and another size is `sd`.
- A numeric column is used as given: never centered and never divided, at the fit, after a swap and at
  new rows. A stated sd is per unit of the column, and nothing is read from the data for it.
- Defaults. A factor takes U, and a number U for one standard deviation of the column, s = U / sd(w),
  held or drawn and whatever the number of forests (dec-B291). A drawn forest with no basis takes 2 U
  under gaussian and U under probit and logistic. A held forest with no basis takes, under gaussian, the
  size the one forest of a single-forest chain has at its defaults, a quarter of the response's range as
  the engine holds it (dec-B290), and U under probit and logistic (dec-B293). A held numeric coefficient
  that states no sd takes the drawn default: right at ordinary slopes, short under probit when the slope
  is two or three times larger, which the help says.
- The variance of a coefficient is not an argument: `amplitude.prior.variance` is gone.
- A gaussian response with no spread over the rows its scale is read from is scaled by a window of
  width one; the anchor is then a quarter of that window, which is what the one forest of a
  single-forest chain gets there, and U is a quarter of the window too.

## 2. A numeric column's scale and its record

sd(w) is the sample standard deviation of the column, n - 1 in the divisor, the mean taken first and
corrected and the deviations from it squared, as R's `sd()` does, unweighted, over the rows in the
likelihood when the sampler is first created (dec-B303). Those are the rows the response's scale was
read from: under gaussian the rows left by `subset` and the na.action with a positive weight, and every
row where none is in or the rows in hold fewer than two distinct values of the response (dec-B302,
dec-B309); under probit and logistic, where nothing is read from the response, the rows the 0/1 weights'
active-row mask keeps, and every row where it keeps none. Either fallback is reported by the engine as a
status and raised by R as one warning, which under gaussian names the response beside the column. A
column constant over those rows, by its values, or whose standard deviation there is at most 2^-26 of
its largest absolute value, has no scale: a forest on it that states no sd is refused by name, the
refusal giving the constant or the ratio, and with an sd stated it is created (dec-B307).

The engine takes it, once, for a numeric forest that states no sd, in
[`basisColumnScale`](../../src/bartcore/combiner.hpp), handed the rows by the K-forest constructor of
[`Chain`](../../src/bartcore/chain.hpp). R records it on the data object, `data@basis.scale`, in the
sampler's `initialize` where the anchor is recorded, for a forest whose entry the object lacks
([`recordBasisScales`](../../R/dbarts.R)), and the bridge hands it back to every later creation
([`applyForestBases`](../../src/R_interface_bartcore.cpp)), checked there and only there: one positive
finite number for the forest's one column. On every path:

| path | the scale |
|---|---|
| first creation, at every door | taken over the rows in and recorded |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | handed from the record: the creation's prior, bit for bit |
| `$setForestBasis` | nothing moves: a swap installs the basis and derives nothing (dec-B452), a column of zeros among them (dec-B346) |
| `$setLeafPrior(forests = )` | the forest is stated at s; the record stays, so a default put back later has the creation's scale |
| `setResponse`, `setOffset`, `setWeights`, `setPredictor`, `setActiveRows` | nothing moves |
| `setState`, a state in a reload | nothing moves: a state holds no scale and installs under the recipient's own record |
| a declaration that gives a forest of a data object another basis | that forest's entry is dropped and taken from the new column |
| a control or a `forests` list carried to another fit | the new fit takes its own: a scale is on the data, never on the control |
| a record cleared by hand | the next creation takes the column then in force |

The mask a latent family's 0/1 weights resolve to is handed to the creation
([`dbartsSampler$initialize`](../../R/dbarts.R)), which installs it as it did before and reads the rows
of a column's scale from it; a mask installed after creation moves no prior.

## 3. The row bound

A row leaves a forest's update, with zero weight and zero response, when its multiplier is at most
2^-26 of the largest absolute multiplier that forest has in the sweep, and for no smaller reason
([`AmplitudeForestCombiner::formForestResponse`](../../src/bartcore/combiner.hpp)). The division then
amplifies by at most 2^26 relative to the forest's largest multiplier. An absolute bound, which the bound
was before, left a column of values below about 1e-8 with no row at all; relative, a default fit on w
and on c w is the same fit for any c that leaves both representable, and to the bit for a power of two.
The first commit of push 2 is this change alone, gated bitwise on every recorded scenario.

## 4. Where it lives

The kind is stated to the engine at every construction, from the data object: the bridge's
[`parseData`](../../src/R_interface_bartcore.cpp) reads it as R's
[`forestBasisKind`](../../R/forestBasis.R) does, and [`ForestSpec`](../../src/bartcore/combiner.hpp)
carries it as a [`ForestKind`](../../src/bartcore/combiner.hpp), with the scale handed back as
`basisScale`. [`forestMapEntry`](../../src/bartcore/combiner.hpp) resolves the rule from the statement,
the kind, the hold and the scale, and [`defaultForestSd`](../../src/bartcore/combiner.hpp) holds the
defaults; the held default of a forest with no basis is
[`singleForestDefaultSize`](../../src/bartcore/chain.hpp), called and not copied, and so is the anchor on
a response with no spread. One function sets a forest's leaf scale from the entry in force,
[`mapLeafScale`](../../src/bartcore/chain.hpp), at construction and at a write. The combiner sets a held
block to its kind's value, and a drawn numeric coefficient to 1, when it is built
([`AmplitudeForestCombiner::holdAtKindValue`](../../src/bartcore/combiner.hpp)), and throws for a numeric
basis of several columns, for a held block of levels of other than two columns, and for a kind that does
not fit its basis. retired: [`basisRowNorm`](../../src/bartcore/chain.hpp), the median row norm the
leaf scale was divided by, is gone.

The reader ([`reportLeafPrior`](../../R/dbarts.R)) gives `sd` (a numeric forest's per unit of its column
and named by it), `sd.stated`, `multiplier` ("none", "factor", "numeric"), `amplitude` (`fixed()` where
held, absent where drawn), `prior.sd.of` ("forest", "level difference", "forest, per unit of basis",
each with " (prior median)" where drawn) and, on a numeric forest whose data records one, `basis.scale`.
The channels the map carried its numbers in are no longer reported. `extract(type = "leaf.prior.sd")`
returns the entries by label, and one forest's entry itself. `samplePriorPredictive` refuses a sampler
of several forests by name, their coefficients not being drawn from their prior there.

## 5. What a state and a saved sampler hold

A state holds the coefficients and one variance a forest, no sd, no scale, no kind; every install puts
them in as stored (dec-B418), and a held block in the state is passed by. A saved sampler holds the
scale on its data object and the anchor and the unit on its control, so a reload builds the creation's
prior. Neither the stored state's format nor the shipped header changes. The header has no create
entry and cannot describe a sampler of several forests, so every such sampler is created from R, which
keeps the record.

## 6. What moved

- A model of three or more forests in which a forest with a basis states no sd (push 1): its default is
  one unit where it was sqrt(2 / K) of one. Bit for bit at two forests. Oracle: the base build stating
  the unit on each such forest gives the same draws, to the bit.
- A held two-level factor (push 1): its own sd is s where it was s / 0.674. Oracle: the base build
  stating 0.674 s, equal to 1e-10; and the held arms of benchmarks/R/forest-prior-exact.R.
- A held forest with no basis (push 1): its own sd is s where it was the unit whatever s, and under
  gaussian a quarter of the range by default. Under probit and logistic at the default no bit moves.
  Oracle: the exact gate's arms and the equality with `bart()`'s own forest.
- A gaussian model of several forests on a constant response (push 1): its forests had size 0.
- Every model with a numeric multiplier, drawn, at the default or a stated sd (push 2): the column is
  no longer divided by its median row norm, its coefficient's variance is 1 where it was 1/2, a default
  is per standard deviation of the column over the rows in, and a drawn coefficient starts at 1 where,
  standing second, it started at 0. The base build has no spelling of it. Oracles: the numeric arms of
  forest-prior-exact.R (a 0/1 column and one of four values, default and stated, drawn and held, rows of
  weight zero, and an arm of three trees, the coefficient's probability of lying within 0.5, 1 and 2
  among the functionals), the literals of tests/cpp's `testNumericMultiplierLaw`, and the units twin.
- Accepted where refused (push 2): a held numeric column, at 1. Refused where accepted: a numeric column
  constant, or constant to rounding, over the rows its scale reads, on a forest that states no sd.
- Bit for bit (push 2's first commit): the relative bound, on every recorded scenario, snapshot and pin.
- The BCF baseline: push 1 re-recorded `glue_toggle`; push 2 adds four scenarios, a numeric column at
  its default, through a swap, a write and a copy, held, and a held forest with no basis under gaussian,
  so that a later slice that moves the record is gated against recorded draws.

## 7. What the engine keeps general, and what is left for later

A coefficient for each column of any block, drawn as one block under one variance, a state's glue
block for any widths, and a drawn block of levels of any width under the factor's rule, which R never
reaches. A held block has a value for nothing, two levels and one numeric column. A kind is a statement
and never read off a block's width, so a numeric block of several columns, whichever meaning it is given
after the release, is a new row and no change to how a 1.0 model is read; today it throws.

Each later addition is an addition and not a change of a 1.0 model:

- `$setForestBasis(updateBasisScale = TRUE)`: the argument does not exist, so no call states it; with
  none a swap keeps the scale; the record is a slot a call can rewrite and `basisColumnScale` takes the
  rows it is handed; no state or shipped signature is involved.
- `leaf.prior = normal(sd = )` on `forest()`, and `normal()` for the default back: the reader already
  says stated or not, and a stated forest keeps its record, so a default put back has the creation's
  scale.
- A heavy-tailed coefficient: `amplitude` takes `fixed()` alone today; it is one row of
  `forestMapEntry`, and the state already carries one variance a forest.
- A forest for each level of a factor of three or more, and several numeric columns with a size for
  each: refused today; the reader's `sd` and `basis.scale` and `extract`'s entries are vectors by
  column, named by the column, and `data@basis.scale` holds a vector for each forest.
- `samplePriorPredictive` on several forests, which needs the coefficients drawn from their prior.
- A drawn coefficient's start by kind for every kind, which would move the draws of every model of three
  forests; today a forest with no basis or a factor starts by position and a number at 1.
- `setData` on a sampler of several forests, refused: when it opens, a new object's record is used and
  what it lacks is taken.

## 8. The printed block

A model of several forests prints a block of one line for each forest, built by
[`forestPriorLines`](../../R/forestPrint.R) from the reader's entries and, for a numeric forest, the mean and
standard deviation of its column as the object holds it: a sampler's data object, a fit's stored `bases`, over
the rows in the likelihood when the block is printed (a fit's stored weights and mask, a sampler's own).
`print` of a fit and `show` of a sampler give the short block, a line's label, its kind, a numeric column's
mean and sd and the forest's sd, a numeric forest's per sd of its column as it stands (per unit where the
column is constant), under a plain header and over a last line that points to the full block. `summary` of a
fit, after its table, and a sampler created with `verbose = TRUE` give the full block: under a header that
names the unit, a numeric forest's sd per unit and per sd of its column, the column's range, where the sd came
from and whether the coefficient is held. Nothing of the column is stored: the mean and sd are those of the
column in force, so after a swap or a change of weights they are the new ones, and a default's full line ends
with the sd of the column it was taken from where that differs. The one number a fit stores for the block is
the unit the sds are stated in, the engine's own, `attr(fit, "forest.unit")`, which a sampler already holds on
its data object; R takes no second unit from the response. The header names it as the help does: sd(y) as the `sd`
item takes it, or the link's error sd. A fit made with `keepFits = FALSE` keeps no bases, labels or unit and
prints no block, and its summary keeps each forest's sd on the line under the table. The line under the table
of a fit that does print the block leaves its sds out. A copy is not a creation and prints nothing. A fit of
one forest prints and summarizes as before. The format is a constraint on nothing: the help names the two
things to read in the full block, a column's mean and range and the sd per sd of the column.
