# leaf-prior-reader-shape: the leaf prior alone, and k as chain state

agent: Opus (R, man and tests; no engine or bridge code)
rng: neutral
budget: about 400 lines over R/, man/, inst/tinytest, inst/NEWS.Rd, benchmarks

## Goal

`$getLeafPrior()` returns the leaf prior alone, as a list, rather than a
chains x 14 matrix. The per-chain values of a drawn k are chain state and move
to a new reader, `$getK()`. Ruling: dec-B141, with dec-A108's calibration
entries and dec-A123's current-value reader as precedent.

## Reader shape

One forest's prior is a named list:

- `leaf.prior`: the specification in the terms it was named in, a
  `normal()`, `linear()` or `gp()` object carrying exactly one of `k` (a
  number or a `chi()` law) or `sd` (a number or an `invchi()` law). It goes
  back into `$setLeafPrior()` or a fitting function's `leaf.prior` as is. An
  unnamed prior is stated as the family default it resolved to. A fixed value
  is read off the engine, so it is what is in force; a law comes from the
  model, with an `invchi()` scale read off the anchor in force. On a forest
  whose scale a calibration map sets, k is pinned at 1 and the entry is
  `normal(sd = )` at the map's leaf scale.
- `leaf.model` and `prior.sd.of`: the former attributes, as elements.
- `prior.mean`, `anchor`, `response.scale`, `response.shift`. `anchor` is the
  one k is relative to, so the spread in force on each chain is
  `anchor / $getK()`: the data's anchor under a k-named prior and under
  `sd = invchi(df, 0)`; otherwise, under an sd-named prior, twice the sd or
  `invchi()` scale in force (the reference k of 2), not the data's anchor,
  which the list does not carry (an accepted cost); the map's leaf scale on a
  map forest.
- On a forest whose scale a calibration map sets, and absent (so `NULL`)
  elsewhere: `basis.row.norm`, `leaf.scale.factor` and `leaf.scale.divisor`
  (`NA` after a state install brings a calibration the map did not derive),
  and one of `amplitude.prior.variance` or `amplitude.prior.scale`.

These are sampler-level: every chain carries the same values, except that
`$setState` accepts chains saved from different samplers, each bringing its
own response transform, leaf scale and k (measured: a two-chain sampler whose
second chain came from a sampler fit to 10y + 3 reports two response scales).
A quantity the chains disagree on reads `NA`, the sampler holding no single
value. An `NA` spread in `leaf.prior` is refused on write by the class
validity, `$setLeafPrior`, `$setModel` and the fitting functions.

At `forest = NULL` a multi-forest sampler returns an unnamed list, one prior
per forest; a single-forest sampler's `NULL` read is its forest-1 read.

## getK

`$getK(forest = NULL)` returns each chain's current k, the value `run()$k`
records per draw, read without running: bitwise the last draw after a run.
A fixed k repeats per chain; a map forest's is 1. It is k whatever terms the
prior was named in, relative to the prior's `anchor`. One forest gives a
length n.chains vector; `NULL` on a multi-forest sampler a forests x chains
matrix, the forest margin before the chains as on the sibling readers.

## Callers

bart's and predict's per-forest channels read `response.scale` and
`response.shift` from the list; rbart asks whether the spec's k or sd is a
law. bartCause's bcf changes in lockstep (orchestrator).

## Landing

Pending.
