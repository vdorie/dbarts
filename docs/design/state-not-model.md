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
| `fit.scale` | the response mapping the chain was under when the state was read | a record | recorded, not compared: no install reads it (dec-B418) |
| `cutPoints`, a leaf's covariate standardization, a heuristic gp lengthscale | the frame the trees and slopes are read through | scratch, frozen | installed as is |
| a supplied gp lengthscale | the kernel | model | kept; saved draws under another are refused |
| weights and censoring digests | which data the latents were drawn against | data, by digest | compared; mismatches are reconciled |

A state whose `cutPoints` hold a value twice in one column is refused, naming the column: no grid repeats a
point ([cut-grid.md](cut-grid.md)), so a stored split's value names one position.

Not carried: the leaf scale, a fixed amplitude prior variance or fixed amplitudes, the tree prior, the move
probabilities, the sigma prior, the variance forest's leaf prior, the monotone directions and the bases. A state
written before these were dropped still installs: a block the reader no longer wants is ignored.

## The response transform

Each stored leaf value is a number on an internal scale. What it means on the response scale depends on the
transform - a multiplier and a shift - of the sampler that reads it. The leaf prior's `k` is stated against
that same transform: its centre is the transform's shift and its width, `k.scale`, a constant of the family
times the range, whichever of `k` and `sd` the prior is written with
([leaf-scale-rules.md](leaf-scale-rules.md)). So the transform is part of the model.

The sampler holds one transform. It is set when the sampler is created, and again when `setResponse` or
`setOffset` with `updateScale = TRUE`, or `setData`, re-anchors it; nothing else moves it. The R object records
it on the model, and a copy or a reload is re-created in it: the re-created chains are moved to the record at
once, before any state goes in, since no install moves them.

Every install - `setState`, `copy()`, a reload, a warm start - puts the state in as stored (dec-B418): trees,
leaf values, slopes, gp fits, variance factors, the kept draws, `k` and sigma are numbers on the internal scale
and are read against the recipient's transform, never converted. Sigma is carried the same way, the state
holding the chain's internal value, so an install writes it back with no pass through response units. A state
stored under the transform in force puts the chain back where it was stored, value for value. A state stored
before a re-anchor, or by a sampler on another response, is legal and has no special meaning: the same
internal numbers are a position the chain could hold, and the next draws move them where the data say. A
response that was only rescaled is then no change at all on the internal scale; one whose centre moved too is
fitted wrongly until the draws catch up. No leaf model refuses such a state: a gp leaf's saved draws and
forests coupled through amplitudes, which had no mean term to carry a converted shift, take it as any other
does. A constant response's transform is the window of width 1 centred on its value, c - 0.5 to c + 0.5
(dec-B386; it was c to c + 1 until 2026-10-08), recorded as (c, c).

A re-anchor and an install differ on the draws the sampler has kept. A re-anchor rewrites them into the new
transform, so they stay the functions they were ([leaf-conversions.md](leaf-conversions.md)); an install
brings them in as stored, so after a restore across transforms `predict` reads them on the new scale.

A re-anchor is a model change, so restoring a state saved before one does not undo it. A re-anchoring
proposal's rollback is a restore of the response with the state: re-anchor to the old response, then restore.

## Where the division holds, and where it does not

It holds for every prior quantity and every value held fixed: an install never changes what `getLeafPrior`,
`getSigmas` on a fixed sigma, `getShape` on a fixed shape or the fixed df report, and every chain of a sampler
runs under one prior. Chains spliced from several samplers go in as stored and run as one posterior under the
sampler's one transform; a chain drawn under another transform is then read on a scale it was not drawn on, so
a host that splices chains across transforms restores each into a sampler at that chain's own.

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

## A state of another family

Measured for [cross-family-state-install.md](../plans/cross-family-state-install.md); the rule is dec-B283's.

A state does not say which response family stored it, and its blocks do not show it. A probit state, a logistic
one and a hazard one hold the same blocks, and the latent block is a latent response under probit, ordinal and
aft and a precision - a Student-t scale, a Polya-Gamma variate - under Student-t, logistic and negative
binomial. So `setState`, `copy` and a reload do not ask the family: a state whose blocks fit the sampler is
installed, and what the sampler then holds of another family's state is not promised.

One thing is refused for its values. A Student-t, logistic or negative-binomial sampler refuses a latent block
holding a value that is not positive and finite, in the words of every other state that does not fit, `state is
not consistent with this sampler`. Its sweep divides by those values and draws against them, and latent
responses in their place, many of them negative, leave fits that are not finite or a sweep that does not
return. The other families' latents are real numbers and are not judged. After the refusal the sampler, its
stored state and its generators are as they were.

On 60 rows, 5 trees and one chain, each ordered pair in a process of its own, an install and then 23 sweeps
(rows: the state; columns: the sampler; "refused" is that message):

| state | gaussian | Student-t | probit | logistic | ordinal | nbinom | aft | multinomial |
|---|---|---|---|---|---|---|---|---|
| gaussian | - | refused | runs | runs | refused | refused | runs | refused |
| Student-t | refused | - | runs | runs | refused | runs | runs | refused |
| probit | refused | refused | - | refused | refused | refused | runs | refused |
| logistic | refused | runs | runs | - | refused | refused | runs | refused |
| ordinal | refused | refused | runs | refused | - | refused | runs | refused |
| nbinom | refused | runs | runs | runs | refused | - | runs | refused |
| aft | refused | refused | runs | refused | refused | refused | - | refused |
| multinomial | refused | refused | refused | refused | refused | refused | refused | - |

Of the 56 pairs 38 are refused and 18 install and run. Seven are refused for the precisions alone: a probit,
ordinal or aft state by a Student-t and by a logistic sampler, and an aft state by a negative-binomial one. The
hazard and two-forest forms of probit are refused by a logistic sampler as the plain pair is. A gaussian state
holds no latents and installs into a probit, logistic or aft sampler. Among the pairs of different families
that install, most leave the sampler holding the other family's latent block as stored until its next sweep.
Pairs of one family that differ in a value held fixed - the Student-t df or the negative-binomial shape, fixed
in one and drawn in the other - install, as does a monotone sampler's state into a plain one; a plain state
goes into a monotone sampler only when its leaf values are in order. The warm start, `installTrees`, reads no
latents and is the way to start a sampler from a fit of another family: it takes the donor's trees and, where
the sampler draws them, its sigma, k and DART state. A linear leaf's coefficients are restated from the
donor's covariate standardization into the sampler's own ([leaf-conversions.md](leaf-conversions.md)).

## Where it lives

The writer and the install rule: [`Chain::getState`](../../src/bartcore/chain.hpp),
[`Chain::setState`](../../src/bartcore/chain.hpp), [`installForest`](../../src/bartcore/chain.hpp),
[`installDrawnScalars`](../../src/bartcore/chain.hpp) and
[`AmplitudeForestCombiner::restoreGlue`](../../src/bartcore/combiner.hpp). The transform:
[`Sampler::setAnchor`](../../src/bartcore/sampler.hpp) and [`moveScale`](../../src/bartcore/chain.hpp); a
re-anchor's rewrite of the kept draws, [`restateSavedDraws`](../../src/bartcore/chain.hpp). The record:
[`recordAnchor`](../../R/dbarts.R), [`applyAnchor`](../../R/dbarts.R) and
[`recreatePointer`](../../R/dbarts.R). The floor on precisions:
[`ResponseModel::canHoldLatents`](../../src/bartcore/model.hpp), asked by
[`Chain::stateIsValid`](../../src/bartcore/chain.hpp).
