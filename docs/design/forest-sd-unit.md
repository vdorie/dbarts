# A forest's sd is in the response's units

Status: LANDED 2026-10-11 (e84c0dd0 to 37b59ea0, calls dec-A208). The plan is
[forest-sd-unit.md](../plans/forest-sd-unit.md), which carries the steps, the gates and the calls made in
planning; this note is the standing statement of the rule. Rulings: dec-B246 as revised, dec-B253, dec-B266,
dec-B269 as dec-B280 and dec-B282 revise it, and dec-B296, dec-B302 and dec-B309 for the rows the unit is read
from, in docs/decisions.md.

A forest's `sd` is a number of the response's units wherever it is written or read: `forest(sd = )` at
creation, `$setLeafPrior(forests = )`, what `$getLeafPrior()` returns, and `extract(type = "leaf.prior.sd")`.
The units are those of y under a gaussian response and those of the latent index under probit and logistic.

## 1. The rule

1. The unit. One number of the sampler converts a stated sd: the response's standard deviation under
   gaussian, `sd(y - offset)` over the rows in the likelihood when the sampler is created
   ([response-scale-rows.md](response-scale-rows.md)); 1 under probit and pi / sqrt(3) under logistic, the
   standard deviation of the link's own error distribution.
2. Stated. The engine is handed the number as written and divides it by the unit. What it then does with the
   quotient is what the calibration map did with a number in units of the response's scale.
3. Not stated. The engine is told so and takes the default itself, a multiple of the unit that never passes
   through the unit on the way in: 2 under gaussian and 1 under a latent family for a forest with no basis,
   sqrt(2 / K) for a forest with a basis in a model of K forests.
4. Who holds the map. The engine, in one function
   ([`forestMapEntry`](../../src/bartcore/combiner.hpp)), with the defaults in another
   ([`defaultForestSd`](../../src/bartcore/combiner.hpp)). A forest's statement
   ([`ForestStatement`](../../src/bartcore/combiner.hpp)) is an sd or none and the hold; whether the forest
   has a basis was the one thing of its kind the map read, until forest-multiplier's push 1 stated the kind
   and dropped `amplitude.prior.variance` ([forest-multiplier.md](forest-multiplier.md)). R derives none of it and divides nothing.
5. Read. `$getLeafPrior(f)` returns `leaf.prior`, the statement (`forest(sd = v)`, v the number as written,
   where one was stated, and `forest()` where none was); `sd`, the number in force in the response's units,
   which is v bit for bit or the default's multiple times the unit; and `sd.stated`. With no forest it returns
   those lists named by the forests' labels.
6. Write. `$setLeafPrior(forests = )` with `forest(sd = v)` states v: the forest is a stated one from then on,
   also where v is the number in force, in which case no bit of the prior moves
   ([`Chain::setForestMapSd`](../../src/bartcore/chain.hpp)). An entry `forest()` leaves its forest as it is.
7. `extract(type = "leaf.prior.sd")` on a fit of several forests returns the reader's `sd` entries as a list
   named by label, and the entry itself where one forest is selected.
8. The unit is derived once, by the engine, at a first creation, and recorded beside the anchor
   ([`derivedForestSdUnit`](../../src/bartcore/chain.hpp)). Every re-creation is handed the record: a copy, a
   reload, a sampler made again from its own control, model and data. No later change of the response, the
   weights or a mask moves it. A control or a `forests` list carried to other data states the same numbers of
   that data's units: the unit is not carried over, and the new sampler derives its own.
9. One forest. `forest(sd = s)` on the one forest of a single-forest model is that model's
   `leaf.prior = normal(sd = s)`, in every family, with that family's own verdict on a stated sd; stated at both
   places it is refused.
10. A name on a single `sd` is dropped, at `forest()` and in the writer. A longer `sd` is refused.

## 2. What the number is the size of

Measured with the engine's own prior draws under a flat likelihood; each entry is the prior median of the
absolute size of the forest's contribution over sd times the unit.

| forest | gaussian | probit | logistic |
|---|---|---|---|
| no basis, coefficient drawn, default | 0.996 | 0.977 | 1.002 |
| no basis, coefficient drawn, `sd = 0.7` | 0.967 | 0.987 | 0.955 |
| no basis, coefficient held, default | 0.504 | 1.002 | 1.002 |
| no basis, coefficient held, `sd = 0.7` | 1.423 | 1.415 | 1.430 |
| factor, coefficient drawn, default; `sd = 0.7` | 1.014; 0.948 | 1.012; 1.006 | 1.023; 0.990 |
| factor, coefficient held, `sd = 0.7`; default | 1.483; 1.494 | 1.494; 1.498 | 1.474; 1.477 |
| one numeric column, coefficient drawn, default; `sd = 0.7` | 0.713; 0.691 | 0.701; 0.712 | 0.736; 0.693 |

So `sd` is the size of the forest's contribution where the forest has no basis or a factor and its
coefficient is drawn. Where the coefficient is held, or the basis is a numeric column, it is not: a held
forest with no basis has the size of the unit whatever its sd is, a held factor 1.48 times it, and a numeric
column 0.71 times it at a row whose basis has the median nonzero norm. The help says so until the multiplier
changes it.

## 3. Constraints

- A fit that states no forest sd draws what it drew before the unit moved, to the bit, in every family.
- The calibration map is not touched: which channel a forest's number reaches, the 0.674, the coefficient's
  variance, the row norm, the held values and the starts by position, the defaults and the K in them. Every
  expression keeps its order of operations.
- The conversion is made once for each construction and once for each restatement, in
  [`forestMapEntry`](../../src/bartcore/combiner.hpp), against the recorded unit. A stated sd equal to the
  number its default comes to ([`forestSdInForce`](../../src/bartcore/combiner.hpp)) resolves to the default's
  multiple and not to the quotient, which can differ from it in the last place: a default read back and stated
  gives the prior it already had, at the write and at every later construction. A caller who states a number
  that happens to equal the default's gets the default's bits.
- A stated sd on a response with no spread (a unit that is not positive and finite) is refused, at creation
  and at a write, with a text that names nothing of the engine. A model of several forests on a constant
  response is created when no sd is stated.
- The forests' record on the control holds four numbers a forest (tree count, base, power, the hold), `sd`
  and `variance` (one number or NA a forest), `hasBasis`, the anchor and the unit, with the record's other
  entries as they were. `hasBasis` stays: the four numbers do not record which channel the forest's sd
  reached, which the map reads from the forest's first kind for the life of the sampler, and
  `$setForestBasis` refuses a basis for a forest created with none
  (retired: [`ForestSpec::createdWithoutBasis`](../../src/bartcore/combiner.hpp), gone with
  forest-multiplier's push 1, which states each forest's kind from the data object and drops
  `hasBasis` from the record). The record can mark a forest as created without a basis; it cannot
  give a forest with no basis in the data one.
- The flat C header, the stored state and every state block are untouched.

## 4. The changed sequence and its oracle

A model of several forests that states `sd` on a forest under a gaussian or a logistic response and is left
as written is another prior: 0.7 was 0.7 standard deviations of the response and is 0.7 units of it. The
oracle is the same model on the base build with the number divided by the unit; respelled as `r * sd(y)` the
prior is the same to rounding, and draws agree to 1e-14 to 3e-14 where they are not identical bits. No
recorded scenario and no snapshot file states a forest sd. The four BCF exact gates and the BCF arm of
`sbc.R` state theirs as multiples of the response's standard deviation or of the link's scale.

## 5. What it leaves

The reader's `amplitude.prior.variance`, `amplitude.prior.scale`, `leaf.scale.factor`, `leaf.scale.divisor`,
`basis.row.norm` and `prior.sd.of` stay the engine's channels, in multiples of the unit, and are said to be so
in the help. A numeric forest's `sd` is not named by its column, and the printed block of a fit is unchanged.
The multiplier removes or renames these.
