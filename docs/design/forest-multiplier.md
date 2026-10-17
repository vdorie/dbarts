# What a forest's coefficient is, held or drawn

Status: PUSH 1 LANDED 2026-10-17 (5a62933e); pushes 2 and 3 planned. The plan is
[forest-multiplier.md](../plans/forest-multiplier.md), which carries the steps, the gates and the calls.
This note states the rule push 1 builds: the held half of it, the default of a forest with a basis, and
the unit on a response with no spread. Rulings: dec-B290, dec-B291 and dec-B293 for the defaults,
dec-B353 and dec-B354 for a held coefficient, dec-B376 for what a fit records of a hold, dec-B452 for
a kind that no swap changes, in docs/decisions.md. Push 2 completes the note with the numeric multiplier.

## 1. The rule, push 1's half

A forest of a model of several multiplies nothing, the two levels of a factor, or a numeric column,
by the class of its basis ([forest-kind-by-class.md](forest-kind-by-class.md)). Its term at row i is
the basis row times the forest's coefficients, times the forest's own fit F. s is the forest's `sd`
in the response's units; U is the unit ([forest-sd-unit.md](forest-sd-unit.md)).

| the forest multiplies | drawn | held, `amplitude = fixed()` |
|---|---|---|
| nothing | F has sd U, the coefficient's size half-Cauchy with median s / U | the coefficient is 1; F has sd s, exactly |
| two levels | each level's coefficient N(0, 1/2); F has sd s / 0.674 | (0, 1); F has sd s, exactly, the second level against the first |
| a numeric column | as before this slice, until push 2 | refused (dec-B353) until push 2 |

- A held coefficient is 1, and what it multiplies is the basis (dec-B354): 1 with no basis, 0 and 1
  for a factor's two levels. It goes by the forest's kind and never by where the forest stands.
  `fixed()` takes no value but 1; another coding is another basis, and another size is `sd`.
- Defaults. A forest with a basis takes one unit, U, held or drawn and whatever the number of forests
  (dec-B291). A drawn forest with no basis takes 2 U under gaussian and U under probit and logistic. A
  held forest with no basis takes, under gaussian, the size the one forest of a single-forest chain
  has at its defaults, a quarter of the response's range as the engine holds it (dec-B290), and U
  under probit and logistic (dec-B293).
- The variance of a coefficient is not an argument: `amplitude.prior.variance` is gone, and a
  forest's coefficients with a basis are drawn under variance 1/2.
- A gaussian response with no spread over the rows its scale is read from is scaled by a window of
  width one; the anchor is then a quarter of that window, which is what the one forest of a
  single-forest chain gets there, and U is a quarter of the window too.

## 2. Where it lives

The kind is stated to the engine at every construction, from the data object: the bridge's
[`parseData`](../../src/R_interface_bartcore.cpp) reads it as R's
[`forestBasisKind`](../../R/forestBasis.R) does, a basis with recorded levels being levels and any other
basis numbers, and [`ForestSpec`](../../src/bartcore/combiner.hpp) carries it as a
[`ForestKind`](../../src/bartcore/combiner.hpp). [`forestMapEntry`](../../src/bartcore/combiner.hpp)
resolves the rule from the statement, the kind and the hold, and
[`defaultForestSd`](../../src/bartcore/combiner.hpp) holds the defaults; the held default of a forest
with no basis is [`singleForestDefaultSize`](../../src/bartcore/chain.hpp), called and not copied, and so
is the anchor on a response with no spread ([`Chain::latentScaleAnchor`](../../src/bartcore/chain.hpp)).
The combiner sets a held block to its kind's value when it is built, after the layout that keeps every
drawn block's start ([`AmplitudeForestCombiner::holdAtKindValue`](../../src/bartcore/combiner.hpp)), and
throws for a held block with no value (levels of other than two columns, a numeric column) and for a
kind that does not fit its basis. A held block is passed by at a state install, as before.

R refuses a held numeric column by its kind ([`refuseHeldShape`](../../R/model.R)). The reader
([`reportLeafPrior`](../../R/dbarts.R)) gives `multiplier` ("none", "factor", "numeric"), `amplitude`
(`fixed()` where held, absent where drawn) and `prior.sd.of`: "forest" or "level difference", with
" (prior median)" where the coefficient is drawn; "forest total" for a numeric column until push 2.

## 3. What moved

- A model of three or more forests in which a forest with a basis states no sd: its default is one unit
  where it was sqrt(2 / K) of one. Bit for bit at two forests. Oracle: the base build stating the unit
  on each such forest gives the same draws, to the bit.
- A held two-level factor: its own sd is s where it was s / 0.674. Oracle: the base build stating 0.674
  s, equal to 1e-10; and the held arms of benchmarks/R/forest-prior-exact.R.
- A held forest with no basis: its own sd is s where it was the unit whatever s, and under gaussian a
  quarter of the range by default. Under probit and logistic at the default no bit moves. Oracle: the
  exact gate's arms and the equality with `bart()`'s own forest.
- A gaussian model of several forests on a constant response: its forests had size 0.
- The BCF baseline's `glue_toggle`, which holds the treatment forest's coefficients, re-recorded.

## 4. What the engine keeps general

A coefficient for each column of any block, drawn as one block under one variance, a state's glue
block for any widths, and a drawn block of levels of any width under the factor's rule, which R never
reaches. A held block has a value only for nothing and two levels. A kind is a statement and never read
off a block's width, so a numeric block of several columns, whichever meaning it is given after the
release, is a new row and no change to how a 1.0 model is read.
