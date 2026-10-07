# A saved state holds the chain, not the model

Status: LANDED 2026-10-02 (3f2c46fc to 8b5191d0). Plan: docs/plans/state-not-model.md. Rulings: dec-B195, dec-B196,
dec-B197 and dec-B200 in docs/decisions.md; the calls made under them, dec-A146 and dec-A148.

The package's original division was four kinds of thing. Data is what the caller supplies. The model is the
priors and the values held fixed. State is where the sampler is as it walks the posterior that model and data
define. Scratch is whatever can be worked out again from the other three. A saved state should carry state and
nothing else, and installing one should leave the model alone.

## What a state carries

| block | what it is | kind | on install |
|---|---|---|---|
| trees, leaf values, a linear leaf's slopes, a gp leaf's fits | the forest | state | installed |
| saved draws (the six `saved.*` blocks) | the kept trees `predict` replays | state | installed |
| latents, ordinal thresholds | the augmentation variables | state | installed |
| DART split weights and delay counter | the split prior's current draw | state | installed |
| the generator | the random stream | state | installed |
| variance-forest trees, live and saved | the scale surface | state | installed |
| amplitudes, and a scale-mixture forest's amplitude variance | drawn coupling values | state | installed |
| sigma, k, the Student-t df, the negative-binomial shape, the DART concentration | state where drawn, model where fixed | written only where drawn | installed only where the recipient draws it; absent keeps the recipient's |
| `fit.scale` | the response units the stored numbers are in | the units of the state | compared with the sampler's; converted when they differ |
| `cutPoints`, a leaf's covariate standardization, a heuristic gp lengthscale | the frame the trees and slopes are read through | scratch, frozen | installed as is |
| a supplied gp lengthscale | the kernel | model | kept; saved draws under another are refused |
| weights and censoring digests | which data the latents were drawn against | data, by digest | compared; mismatches are reconciled |
| `family` | which response family's chain this is | model, by name | compared with the sampler's; never installed |

Not carried: the leaf scale, a fixed amplitude prior variance or fixed amplitudes, the tree prior, the move
probabilities, the sigma prior, the variance forest's leaf prior, the monotone directions and the bases. A state
written before these were dropped still installs: a block the reader no longer wants is ignored.

## The response transform

Each stored leaf value is a number on an internal scale. What it means on the response scale depends on the
transform - a multiplier and a shift - in force when it was stored. A `k`-named leaf prior is stated against
that same transform: its centre is the transform's shift and its width a constant times the range. So the
transform is both the units of the state and part of the model.

The sampler now holds one transform. It is set when the sampler is created, and again when `setResponse` or
`setOffset` with `updateScale = TRUE`, or `setData`, re-anchors it; nothing else moves it. The R object records
it on the model, so a copy or a reload is re-created in it. A state stored in other units has its numbers
rewritten into the sampler's as it is installed: every leaf value is multiplied by the ratio of the two ranges
and shifted by the difference of the shifts, split evenly over the trees; slopes and gp fits take the ratio
alone, and variance factors its square, split over the variance trees. Amplitudes are multipliers and stay
as they are: the leaf values beneath them carry the ratio. A constant response's transform spans 1
upward from its value, the window from c to c + 1, and converts like any other. The replayed function agrees
with the stored one to rounding. A state in the sampler's own units is not touched and installs bit for bit.

Two kinds of chain cannot take a different shift - a gp leaf's saved draws carry no mean term, and with
amplitudes no single forest owns the location - and a state stored under another shift is refused there by
name, a gp state whether or not it holds saved draws. Another range at the same shift converts.

A re-anchor is a model change, so restoring a state saved before one does not undo it: the state is converted
into the new units. To roll a re-anchoring proposal back, re-anchor to the old response and then restore.

## Where the division holds, and where it does not

It holds for every prior quantity and every value held fixed: an install never changes what `getLeafPrior`,
`getSigmas` on a fixed sigma, `getShape` on a fixed shape or the fixed df report, and every chain of a sampler
runs under one prior. Chains spliced from several samplers, as stan4bart's kept-tree replay does, are converted
into one set of units and run as one posterior.

It does not hold in four places.

- The frame. The cut grid, a linear or gp leaf's covariate standardization and a heuristic lengthscale are
  derived from data once and then kept while the data moves, so only the state still has them. They shape a
  prior too - the split rule, the slope prior, the kernel - and they install with the state as they are. This
  is the TODO item state-frame-prior.
- Writes through the C API. A host that re-anchors or writes a fixed sigma through `dbarts.h` changes the
  engine and not the R object, so a sampler re-created from that object holds the last values R saw.
- A sampler made from a model and data as a first creation - `new("dbartsSampler", control, model, data)` -
  anchors to that data, whatever the model it was handed recorded. That keeps a model reused on other data
  anchored to its own data, as before; a host that wants the saver's anchor re-creates through `copy` or a
  reload. stan4bart's restored samplers are first creations, so the prior they report is anchored to the
  data, not to the per-chain anchors its chains drew under; its replay reads no prior.
- A value written to a sampler that holds it fixed is recorded on the model for sigma (`setSigma`) and for the
  leaf prior (`setLeafPrior`), but there is no writer for a fixed df, shape or concentration after creation,
  so the question does not arise for them.

## The family a state was stored under

Added 2026-10-07 ([cross-family-state-install.md](../plans/cross-family-state-install.md)).

The blocks of a state do not say which response family drew them. A probit state, a logistic one and a hazard
one hold the same blocks, and the latent block is a latent response under probit, ordinal and aft and a
precision - a Student-t scale, a Polya-Gamma variate - under Student-t, logistic and negative binomial. Before
the record, `setState` took most states of another family. Measured on 60 rows, 5 trees and one chain, each
ordered pair in a process of its own, an install and then 23 sweeps (rows: the state; columns: the sampler;
"refused" is `state is not consistent with this sampler`):

| state | gaussian | Student-t | probit | logistic | ordinal | nbinom | aft | multinomial |
|---|---|---|---|---|---|---|---|---|
| gaussian | - | refused | runs | runs | refused | refused | runs | refused |
| Student-t | refused | - | runs | runs | refused | runs | runs | refused |
| probit | refused | not finite | - | no return | refused | refused | runs | refused |
| logistic | refused | runs | runs | - | refused | refused | runs | refused |
| ordinal | refused | not finite | runs | no return | - | refused | runs | refused |
| nbinom | refused | runs | runs | runs | refused | - | runs | refused |
| aft | refused | not finite | runs | no return | refused | no return | - | refused |
| multinomial | refused | refused | refused | refused | refused | refused | refused | - |

Of the 56 pairs 31 were refused, 18 installed and ran, 3 left fits that were not finite and 4 a sweep that did
not return. Every breaking pair put latent responses, many of them negative, where the sampler holds
precisions; the hazard and two-forest forms of probit broke a logistic sampler as the plain pair did.

The rule. Every state carries a top-level attribute `family`, one string: the family the sampler runs, as the
`family` argument spells it. A hazard sampler writes its link's family, the model it runs on its expanded
rows. The leaf model, the forests and a variance forest are not in it. The record is model by name: it is
compared with the sampler's own and never installed, as the two digests are, so a state still carries no
model. On an install, by `setState`, by `copy` and by a reload alike:

1. The record is the sampler's family: the state is installed by the other rules.
2. The record is another family's: the state is refused, naming both families, before anything else of it is
   read. A state of another family is named as that whatever else about it differs.
3. The record is present and is not one string: the state is refused as malformed.
4. There is no record, as on a state stored before it existed: the state is installed by the other rules.
5. Whatever the record says, a Student-t, logistic or negative-binomial sampler refuses a latent block holding
   a value that is not positive and finite. The other families' latents are real numbers and are not judged.

After any of these refusals the sampler, its stored state and its generators are as they were. Every pair of
different families is refused, the 18 that ran included. A gaussian state holds no latents and used to install
into a probit, logistic or aft sampler; it is refused with the rest, because a chain's trees mean something
only under the family they were drawn in, which is the line dec-B254 draws for `setModel`. Pairs of one family
that differ in a value held fixed - the Student-t df or the negative-binomial shape, fixed in one and drawn in
the other - or in a leaf constraint install as before. The warm start, `installTrees`, reads neither the
record nor the latents, and is the way to start a sampler from a fit of another family.

## Where it lives

The writer and the install rule: [`Chain::getState`](../../src/bartcore/chain.hpp),
[`Chain::setState`](../../src/bartcore/chain.hpp), [`installForest`](../../src/bartcore/chain.hpp),
[`installDrawnScalars`](../../src/bartcore/chain.hpp) and
[`AmplitudeForestCombiner::restoreGlue`](../../src/bartcore/combiner.hpp). The units:
[`Sampler::setAnchor`](../../src/bartcore/sampler.hpp), [`unitsDiffer`](../../src/bartcore/sampler.hpp),
[`convertStateUnits`](../../src/bartcore/chain.hpp) and [`moveScale`](../../src/bartcore/chain.hpp). The record:
[`recordAnchor`](../../R/dbarts.R), [`applyAnchor`](../../R/dbarts.R) and
[`recreatePointer`](../../R/dbarts.R). The family record:
[`stateFamilyName`](../../src/R_interface_bartcore.cpp), written by
[`storeState`](../../src/R_interface_bartcore.cpp) and compared first in
[`setState`](../../src/R_interface_bartcore.cpp); the floor on precisions:
[`ResponseModel::canHoldLatents`](../../src/bartcore/model.hpp), asked by
[`Chain::stateIsValid`](../../src/bartcore/chain.hpp).
