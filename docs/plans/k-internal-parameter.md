# k-internal-parameter: k is the chain's parameter on the internal scale, and the sd spelling is translated

Status: PLANNED 2026-10-09 for dec-B414, for one ruling on the whole. Nothing is built. A blind critique
follows; then the ruling table goes to the maintainer.

agent: opus implementer, one; one opus reviewer who runs the mutants under Tests.
rng: SHIFTING, as [RNG classes and their gates](README.md#rng-classes-and-their-gates) defines it: the
posterior never changes. NEUTRAL, bit for bit, for every fit whose leaf prior is written with k, in every
operation but a state install or warm start across response units with a drawn k (dec-B384, planned
already). Draws move for fits written with an sd, by rounding only (ran: within 2.3e-14 over 1000
draws); and for the operations whose rule changes: `setModel` with another leaf prior, an install
across a change of prior into or out of the sd spelling, a re-anchor under a drawn sd, and an install or
warm start across response units with a drawn k.
window: before the merge to main; engine slices stay serial. Recommended first of the three plans that
touch the leaf scale (Open call 4).
budget: planned ~950 lines (engine ~150, bridge ~40, R ~60 net of removals, tests/cpp ~200, tinytest
~350, help ~60, docs ~90). Forecast and stops: [Budget and stops](#budget-and-stops).

## Summary for the maintainer

Every chain holds one number for the leaf spread, k, and the spread is a reference scale divided by k.
Today that reference scale depends on how the prior was written: with k it is the data's own scale
(half the response range for a continuous response, 3 for probit, about 5.44 for logistic); with an sd it
is twice the named sd, so that k sits at 2. Because the yardstick moves with the spelling, keeping k and
keeping the spread are the same only until the spelling changes, and the rulings so far keep one in some
operations and the other in others.

Under the rule you asked for, the yardstick never moves: it is always the data's scale. k is the
chain's parameter against it. Writing the prior with an sd becomes a translation: an sd of s is k equal
to the data's scale over s, and a scaled inverse chi prior on the sd with scale s is a chi prior on k
with scale the data's scale over s, the reciprocal pair the help already states. Since nothing else
redefines the yardstick, keeping k is keeping the spread, and the two rulings that seemed to disagree
(a change of prior keeps the spread; a restored state keeps its k) now say the same thing.

What changes for a user:

- A prior written with k: nothing. Every such fit is bit for bit what it is today, except a state or
  warm start moved onto a response of another spread, which already had a ruling (keep the spread)
  waiting to be built.
- A prior written with an sd: the same model and the same posterior. Each k is today's times one
  constant, and the fits agree with today's to about 1e-14. The k a fit reports is a
  different number (the data's scale over the spread, not twice the named sd over it), and the reported
  scale k is measured against becomes the data's.
- setModel given another leaf prior now keeps the spread, as setLeafPrior does. Today it keeps k,
  and the spread jumps when the spelling or a named sd changes.
- A state stored before a change of prior, restored after it, keeps its spread. Today it jumps when
  the change crossed between the k and sd spellings or changed a named sd.
- After the response is re-derived (a re-anchor), a drawn sd's current value stretches with the fit
  until the next draw, as a drawn k's does; its prior stays in response units. Today the current value
  stays put while the fit stretches.
- A warm start onto a response of another spread keeps the spread, as a restored state will.

The engine changes, by a small amount: it holds a named sd and translates it whenever the data's scale
moves, which also fixes a gap where the C interface's re-anchor (stan4bart's warm-up loop) did not keep
a named sd. The C interface's signatures do not change. No recorded baseline or snapshot is expected to
move. The work is about 950 lines planned, 1400 to 1900 at the usual overrun. Five calls are open below,
each with a recommendation.

## Goal

On the internal scale the leaf prior's reference scale, `k.scale`, is a constant of the family and the
response transform, never redefined by a prior's spelling, and k is each chain's parameter against it.
`normal(sd = s)` is `k = k.scale / s`, and `normal(sd = invchi(df, s))` is `k = chi(df, k.scale / s)`,
retranslated whenever the transform moves so the named sd keeps its value and distribution in response
units. Only an operation that rescales the leaves rescales k; a re-anchor leaves the leaves and k as they
were; a change of prior touches neither.

## Context

The note this answers: [leaf-scale-rules.md](../design/leaf-scale-rules.md), its
[2. The rules today](../design/leaf-scale-rules.md#2-the-rules-today) and
[4. Two consistent rules](../design/leaf-scale-rules.md#4-two-consistent-rules). The direction, in the
maintainer's words (dec-B414): "What I'm getting at is, what if we kept k.scale fixed and didn't redefine
it, but translated sds and the scale of prior on sds so that they had the intended value and
distribution?" and "But does k.scale change when the response changes? 2 is 2, right? The response is
still mapped to -0.5 to 0.5. The parameter itself of course lags, but that's to be expected."

How the engine holds the leaf scale today (read, on 4cf51b11):
- Each forest holds an internal per-tree leaf scale and k; the spread in force is the leaf scale over k.
  [`resolvedNodeScale`](../../src/bartcore/chain.hpp) sets the leaf scale from the family's internal
  value ([`defaultLeafScale`](../../R/model.R): 0.5 gaussian, 3 probit, pi sqrt(3) logistic, 3 nbinom)
  unless a named prior scale is finite, in which case it is that scale over the transform's multiplier.
  [`priorScaleFactor`](../../src/bartcore/chain.hpp) converts the internal scale to the response-unit
  `k.scale` that [`forestCalibration`](../../src/bartcore/chain.hpp) reports.
- The R model carries a named sd as `prior.scale` = 2 s with k fixed at 2, or with `chi(df, 2)` for
  `invchi(df, s)` ([`resolveLeafPrior`](../../R/model.R), the [`dbartsModel`](../../R/A_class.R) slot).
  `invchi(df, 0)` is `chi(df, Inf)` with no named scale.
- After every re-anchor, re-creation, copy, install and warm start, R writes the named scale back
  ([`reissueNamedLeafSd`](../../R/dbarts.R)) through
  [`bartcore_setLeafPrior`](../../src/R_interface_bartcore.cpp) to
  [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp), which rewrites the internal leaf scale.
  The flat C entries [`dbarts_sampler_setOffset`](../../src/C_interface.cpp) and
  [`dbarts_sampler_setResponse`](../../src/C_interface.cpp) re-anchor in the engine and write nothing
  back (read), so on that path a named sd is not kept in response units today.
- `$setLeafPrior` keeps the drawn spread by multiplying k by the ratio of the two `k.scale` values
  ([`keepDrawnSpread`](../../R/dbarts.R), [`leafKScale`](../../R/dbarts.R),
  [`Chain::scaleDrawnK`](../../src/bartcore/chain.hpp)). [`Chain::setModel`](../../src/bartcore/chain.hpp)
  re-derives the leaf scale from the model and keeps a drawn k.
- A state's k is installed as stored by [`Chain::setState`](../../src/bartcore/chain.hpp) and the warm
  start's [`Chain::installForest`](../../src/bartcore/chain.hpp); both first pass a state in other units
  through [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp)
  ([`Sampler::setState`](../../src/bartcore/sampler.hpp),
  [`Sampler::installForests`](../../src/bartcore/sampler.hpp)), which rescales the leaves by the units
  ratio r and leaves k.
- The k draw is a scaled gamma ([`ChiKHyperprior`](../../src/bartcore/model.hpp)): k squared is gamma
  with rate half the leaves' sum of squares over the leaf scale squared plus half over the chi scale
  squared. Multiplying the leaf scale and the chi scale by one factor multiplies the drawn k by it, for the
  same generator draws.

### Verified on the current build (ran)

Probes scratch/kint/01-today.R, 02-today-more.R and 03-emulate.R on a private install of 4cf51b11
(library scratch/libs/kint, shipped build): gaussian, 200 rows, half range 1.9694, 20 trees, one chain,
200 sweeps under `normal(k = chi(1.5, 2))` leaving k 2.1901 and the spread 0.8992. Each line is what
the build does today; the columns are k, `k.scale`, spread.

| operation | today (ran) |
|---|---|
| creation, `k = chi(1.5, 2)` or `chi(1.5, 4)` | 2, 1.9694, 0.9847 |
| creation, `k = 3` | 3, 1.9694, 0.6565 |
| creation, `sd = 0.5` or `sd = invchi(3, 0.5)` | 2, 1.0000, 0.5000 |
| creation, `sd = invchi(3, 0)` | 2, 1.9694, 0.9847 |
| creation, probit `k = chi(1.5, 2)`; `sd = 0.5`; logistic `k = chi(1.5, 2)` | 2, 3, 1.5; 2, 1, 0.5; 2, 5.4414, 2.7207 |
| `setLeafPrior` into `sd = invchi(3, 1)`, then `invchi(3, 2)` | 2.2241, 2, 0.8992; 4.4482, 4, 0.8992 |
| then into `sd = 0.5`, then `k = chi(1.5, 2)`, then `k = 3` | 2, 1, 0.5; 3.9388, 1.9694, 0.5; 3, 1.9694, 0.6565 |
| `setModel` with the model of an `sd = invchi(3, 1)` sampler | 2.1901, 2, 0.9132 |
| store, `setLeafPrior(sd = invchi(3, 1))`, restore the stored state | 2.1901, 2, 0.9132 |
| that state into a fresh `sd = invchi(3, 1)` sampler; into a fixed `k = 3` one | 2.1901, 2, 0.9132; 3, 1.9694, 0.6565 |
| `setState` onto a response 3 times as wide, k spelling | 2.1901, 5.9082, 2.6977 |
| the same, sd spelling (donor at 2.4954, 0.8015) | 2.4954, 2, 0.8015 |
| `installTrees` into `sd = invchi(3, 1)`, same response | 2.1901, 2, 0.9132 |
| `installTrees` onto a response 3 times as wide, k spelling | 2.1901, 5.9082, 2.6977 |
| `setResponse(3 y, updateScale = TRUE)`, k spelling, drawn; fixed `k = 3` | 2.1901, 5.9082, 2.6977; 3, 5.9082, 1.9694 |
| the same, `sd = 0.5`; `sd = invchi(3, 1)` (at 2.4954) | 2, 1, 0.5; 2.4954, 2, 0.8015 |
| `setOffset(5, updateScale = TRUE)` (shift only) | nothing moves, either spelling |
| one sweep after that re-anchor, `sd = invchi(3, 1)` | 0.7982, 2, 2.5057 |
| fit `sd = invchi(3, 0.5)`: fit$k mean, `k.scale`, `leaf.prior.sd` mean | 1.1836, 1, 0.8685 |
| fit `sd = 0.5`: fit$k, fit$fixed$k, `leaf.prior.sd` | NULL, 2, 0.5 |

The emulation (03-emulate.R): today's k spelling at k = `k.scale` / s computes what the new sd spelling
would. A fixed `sd = 0.5` and `sd = 0.9847` are bitwise identical to `k = k.scale / s` over 1000 draws;
`sd = 1.3` differs by at most 2.26e-14. A drawn `sd = invchi(3, s)` against `k = chi(3, k.scale / s)`
started at `k.scale / s`: every k draw is today's times exactly `k.scale / (2 s)` (to 12 digits, both s),
training fits within 2.04e-14 and 1.84e-14 over 1000 draws, no draw past 1e-12.

## The rulings under this rule

For one ruling on the whole. "Unchanged": the rule and what it builds stand. "Restated": the outcome
stands and the mechanism or wording changes. "Reversed": part of the ruling no longer holds; the part is
named. Built behaviour that changes without reversing a ruling is in the next section.

| ruling | what it says | under this rule |
|---|---|---|
| dec-B201 | `k.scale` is the value k is relative to; under an sd spelling twice the named sd, "so that the engine's k sits at 2"; "I want k to be the engine's k always."; "users have a way to convert between k and standard deviation." | Reversed in its definition: `k.scale` is the data's scale under every spelling. The name, k as the engine's k, and the conversion (spread = `k.scale` / k) stand. |
| dec-A124 | (the reader slice) `k.scale` the data's under k, twice the named sd under sd, "so the data's anchor is no longer readable under an sd-named prior" | Reversed with dec-B201: the data's scale is readable on every fit. |
| dec-B141 | `getK` reports the engine's k whatever the spelling | Unchanged; under an sd spelling the number is `k.scale` / spread. |
| dec-B192 | a fit carries the k the sampler recorded; on an sd-named fit "relative to twice the named scale" | Restated: relative to the data's scale. |
| dec-B193, dec-B194 | extract answers k and `leaf.prior.sd` on every fit; one number for a held value | Unchanged; `leaf.prior.sd` values unchanged, k values move under an sd spelling. |
| dec-B376 | no fit stores the leaf sd; `fit$fixed$k` for a held k | Unchanged; a held sd's `fit$fixed$k` becomes `k.scale` / s. |
| dec-A105 | named by k (relative) or sd (absolute); "A named sd stays in response units across every call that re-anchors the response scale" | Restated: the named sd, fixed or the scale of its `invchi()`, is retranslated at every re-anchor and keeps its value and distribution. A drawn sd's current value lags with the leaves until the next draw ("The parameter itself of course lags", dec-B414); today it stays. |
| dec-B356 | a switch between spellings of a drawn prior keeps the spread, "k becoming k_old x k.scale_new / k.scale_old" ("B. Keep the spread.") | Restated: `k.scale` never differs, so k is untouched and the spread is kept. |
| dec-B369 | a changed `invchi()` scale keeps the spread at the call ("OK, then A. Keep the spread at the call.") | Restated as dec-B356. |
| dec-B392 | a fixed sd turned drawn keeps the sd in force ("They're trying to actively set the sd of the leaf prior, so both of their values should be interpretted as such.") | Restated: the fixed sd s is k = `k.scale` / s, and the drawn prior keeps that k. |
| dec-B393 | every switch keeps the sd into a drawn prior; a stated value is literal; "k.scale always belongs to the new prior's spelling"; "Changing the prior never moves the state; no setter for k exists." | Rule stands; its `k.scale` clause is reversed: `k.scale` belongs to the data. "I don't see changing the prior as changing the state / the active parameter" now holds literally: neither k nor the spread moves. |
| dec-B401 | "The k itself should literally transfer - a parameter is a parameter, regardless of the prior. It may be a bad fit, but so be it." | Unchanged; the spread is now kept too, since no prior moves `k.scale`. The setState plan's H4 closes. |
| dec-B384, dec-B407 | a state installed onto another response spread keeps its spread, k re-expressed; "Keep it for scale changes. If it helps at all, I guess we can think of `k` as just the internal representation and the parameter as the sd." | Restated: k is divided by the units ratio with which the leaves are converted, under every spelling (today the sd spelling's `k.scale` ignores the response, so its formula left k alone). The two readings in dec-B407 now agree. |
| dec-A146 | (agents' call) a warm start follows the setState rule; a drawn k installs as stored | Its install-as-stored part, already revised by dec-B384 for setState, is revised for the warm start too (Open call 3). |
| dec-B343 | "A. The donor gives its trees, and sigma and k only where the new fit draws them (as built)." | Unchanged in what transfers; k is converted with the leaves when the donor is in other units (Open call 3). |
| dec-B200 | a state in other units is converted on install ("Convert on install.") | Unchanged; k joins the leaves in the conversion. |
| dec-B195, dec-B196 | a state carries no model; "Model when fixed, state when drawn." | Unchanged. A fixed named sd's k is model, retranslated with the data's scale. |
| dec-B254 | `setModel` changes parameters, "the leaf scale as k or sd" among them | Unchanged in scope; it now keeps the spread across a change of leaf prior (the note's defect a). |
| dec-B396 | `copy()` stores the current state when none is stored | Unchanged. |
| dec-A107 | `setLeafPrior` changes only the spread or its prior | Unchanged. |
| dec-A121 | (agents' call) the engine's k carried across a switch of spelling | Already revised by dec-B356; moot. |
| dec-A187 | (agents' call) `setLeafPrior` keeps the spread through `scaleDrawnK`, which state-install-keeps-spread would reuse | Superseded: no write moves `k.scale`, and the install's conversion divides k itself; the primitive goes. |
| dec-B371, dec-A190 | the probit rescaling step divides k by the drawn factor | Unchanged: it acts on the internal k and the leaves together. |
| dec-A13 | no cap on a drawn k; `chi(df, Inf)` improper | Unchanged; `invchi(df, 0)` is still `chi(df, Inf)`. |
| dec-B331, dec-B330 | a re-derived scale converts saved draws, the live leaves kept; a gp sampler holding saved draws refuses | Unchanged; "both lag" is how the live leaves already behave. |
| dec-B361 | prior draws under the prior the sampler runs under | Unchanged; the private sampler's named sd is translated by the engine. |
| dec-B362, dec-B364, dec-B405 | several forests refuse a re-derivation; a kept scale (fewer than two values); the status through the C API | Unchanged. A kept scale moves no transform, so nothing is retranslated. |
| dec-B122 | a variance forest recalibrates on a re-anchor | Unchanged (it has no k). |
| dec-B142, dec-B253, dec-B275 | the multi-forest writer and a forest's sd | Unchanged: a named leaf-prior sd is refused on several forests, multinomial and hurdle fits. |

Reversed, in sum: dec-B201's definition of `k.scale` under the sd spelling (with dec-A124's), dec-B393's
`k.scale` clause, and the as-stored warm start across units of dec-A146 (Open call 3). Superseded: dec-A187's
primitive. Every other ruling stands, several with their mechanism restated.

## Operation by operation

Today is the run above. "New" uses the same fixture: data's `k.scale` 1.9694, rho = `k.scale` / (2 s)
the factor between today's and the new k under an sd spelling. New values are arithmetic on the run
values, not runs. "Draws" says which fits move.

| operation | today (ran) | new | draws |
|---|---|---|---|
| creation, k spelling | as above | unchanged | none |
| creation, `sd = 0.5` | k 2, `k.scale` 1, spread 0.5 | k 3.9388, `k.scale` 1.9694, spread 0.5 | rounding (bitwise in the run) |
| creation, `sd = invchi(3, 0.5)` | k starts 2 against 1 | k ~ chi(3, 3.9388), starts 3.9388 (the start translated: spread s, as today) | rounding; k draws times rho |
| creation, `sd = invchi(df, 0)` | `chi(df, Inf)` | unchanged | none |
| creation, probit `sd = invchi(1.5, 1.5)` | k 2 against 3, the k default's chain | rho = 1: identical | none |
| `setLeafPrior` into a fixed value | the stated k; k 2 for a stated sd | the stated k; k = `k.scale` / s for a stated sd | rounding, sd spelling |
| `setLeafPrior` into a drawn prior | k times `k.scale` ratio; spread kept | k untouched; spread kept (0.8992 throughout) | rounding where an sd spelling is in force |
| `setModel` with another leaf prior | k kept; spread 0.8992 to 0.9132 | k kept; spread kept, 0.8992 | the chain after the call, where the spelling or a named sd changes |
| install of a state from the same sampler and prior | as stored, bitwise | unchanged | none |
| install across a change of prior (store, setLeafPrior, restore; or another sampler's state) | k as stored; spread 0.9132 | k as stored; spread 0.8992 | where the k.scale differed: one side sd-spelled |
| install into a sampler holding k fixed | the sampler's | unchanged | none |
| `setState`, `copy()`, reload onto another response spread | k as stored; k spelling 0.8992 to 2.6977; sd spelling kept | k divided by the units ratio: 6.5703, spread 0.8992, either spelling (dec-B384) | every drawn k installed across units |
| warm start across a change of prior | k as stored; 0.8992 to 0.9132 | k as stored; spread kept | where one side is sd-spelled |
| warm start onto another response spread | k as stored; 0.8992 to 2.6977 | k divided by the units ratio; spread kept (Open call 3) | every drawn k warm-started across units |
| re-anchor, k spelling | k kept; spread stretches with the fit | unchanged | none |
| re-anchor, fixed sd `0.5` on 3 y | k 2, `k.scale` 1, spread 0.5 | k retranslated 11.8164, `k.scale` 5.9082, spread 0.5 | rounding |
| re-anchor, drawn `sd = invchi(3, 1)` on 3 y | k kept, spread 0.8015 held; next draw from invchi(3, 1) | k kept (2.4572 at the new numbers), spread lags to 2.4044 with the fit; prior retranslated to chi(3, 5.9082), the next draw from invchi(3, 1) in response units | the sweep after each such re-anchor |
| re-anchor through the flat C entries, named sd | not retranslated (read): a fixed sd stretches with the response | retranslated, as through R | stan4bart with a named sd in `bart_args` |
| `setResponse`/`setOffset` at `updateScale = FALSE`; a kept scale (dec-B364) | nothing moves | unchanged | none |
| probit rescaling step | k divided by the factor | unchanged | none |

The setLeafPrior chain, step by step under the new rule (same fixture): before k 2.1901, spread 0.8992;
into `sd = invchi(3, 1)` k 2.1901, 0.8992; into `invchi(3, 2)` the same; into `sd = 0.5` k 3.9388, 0.5;
into `k = chi(1.5, 2)` k 3.9388, 0.5 (as today); into `k = 3` k 3, 0.6565 (as today).

## What a fit and the readers report

For a k-spelled fit, nothing changes. For an sd-spelled fit (fixture values):

- `fit$k`, `fit$first.k`, `extract(fit, "k")`, `summary(fit, vars = "k")`, `getK`, `run()$k`: k times
  rho, the data's scale over the spread (fit mean 1.1836 to 2.331).
- `fit$fixed$k` for a held sd: `k.scale` / s (2 to 3.9388).
- `fit$leaf.prior$k.scale` and `getLeafPrior()$k.scale`: the data's scale (1 to 1.9694; probit 3).
- `fit$leaf.prior$leaf.prior` and `getLeafPrior()$leaf.prior`: unchanged, the prior as named, written
  back through `setLeafPrior` moving no bit.
- `extract(fit, "leaf.prior.sd")`, summary's default leaf-scale row on an sd-named fit, print: unchanged
  values (to rounding). The extractor already computes `k.scale` / k from the fit's own pair
  (["leaf.prior.sd"](../../R/diagnostics.R)), so a fit saved before the change still reads correctly.
- The C API ([dbarts.h](../../inst/include/dbarts/dbarts.h)): [`dbarts_results`](../../inst/include/dbarts/dbarts.h)
  and [`dbarts_draw`](../../inst/include/dbarts/dbarts.h) carry k per draw; for an sd-named sampler the
  numbers move by rho. No entry, struct or enum changes, so
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move (it covers signatures, enums
  and layouts, read): not an ABI event. The header states nothing about what k is relative to (read);
  one sentence is added on the `k` field.
- stan4bart (bartcore 963956b, read): copies the k draws into its fit's `k` and reads
  `dbarts_sampler_kIsSampled`; it forwards `bart_args$leaf.prior` to dbarts, so a user may name an sd
  there, and then calls `dbarts_sampler_setOffset` with `updateScale` true during warm-up. Its k numbers
  move for such a fit and its named sd is now kept across the warm-up re-anchors. Its tests write k
  (`normal(k = 3)`, test-09-bartArgs.R, read). bartCause reads `fit$k` from a k-spelled fit; treatSens
  writes `chi()` and fixed k; bairrtt names no leaf prior (read, their dbarts-1.0 and main branches).

## Where it is built

The engine changes; the translation cannot live in R alone. Creation needs the data's scale, which the
engine computes from the rows, and the flat C re-anchor never returns to R (Open call 1).

Engine ([chain.hpp](../../src/bartcore/chain.hpp), [sampler.hpp](../../src/bartcore/sampler.hpp),
[facade.hpp](../../src/bartcore/facade.hpp)):
1. A forest holds `namedSd`, NaN when the prior is written with k. The `priorScale` fields of
   [`ModelParameters`](../../src/bartcore/chain.hpp) and the creation options become `namedSd`.
2. The internal leaf scale is always the family's: [`resolvedNodeScale`](../../src/bartcore/chain.hpp)
   loses its named branch, at creation and in [`Chain::setModel`](../../src/bartcore/chain.hpp).
3. One translation routine: where `namedSd` is finite, a fixed k is `k.scale / namedSd` and a drawn
   prior's chi scale is `k.scale / namedSd`, with `k.scale` the reader's own expression
   ([`priorScaleFactor`](../../src/bartcore/chain.hpp) times the leaf scale), so the reader and the
   write agree bit for bit. At creation a drawn k starts at the same value. A drawn k is never written by
   it.
4. It runs wherever a chain's transform moves: [`Chain::setResponse`](../../src/bartcore/chain.hpp) and
   [`Chain::setOffset`](../../src/bartcore/chain.hpp) at `updateScale`,
   [`applyNewData`](../../src/bartcore/chain.hpp), [`moveScale`](../../src/bartcore/chain.hpp) (under
   [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp), both arms), and
   [`Chain::installForest`](../../src/bartcore/chain.hpp)'s scale restore. tests/cpp holds the invariant
   after every mutation the fuzz and mutation harnesses make (step 9).
5. [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp) becomes a named-sd writer: it sets
   `namedSd` and translates, leaving a drawn k; a write of the value in force is skipped. Its facade
   virtual on [`SamplerBase`](../../src/bartcore/facade.hpp) is renamed with it.
6. [`Chain::setModel`](../../src/bartcore/chain.hpp): a fixed k is the model's or the translation; a
   drawn prior is translated; a drawn k is kept, as today.
7. [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp) divides each forest's state k by the units
   ratio, so [`Chain::setState`](../../src/bartcore/chain.hpp) and the warm start install it with the
   converted leaves (dec-B384; Open call 3). An absent k (NaN) stays absent; a recipient holding k fixed
   ignores it, as today. A ratio of 1 moves no bit.
8. [`Chain::scaleDrawnK`](../../src/bartcore/chain.hpp) and its facade virtual go: no caller remains.
   `--preclean` on every install (facade virtuals change).

Bridge ([R_interface_bartcore.cpp](../../src/R_interface_bartcore.cpp),
[R_interface.cpp](../../src/R_interface.cpp)): the model parse computes the named sd once from the R
model, `prior.scale` over the fixed k or the chi scale (2 s / 2, exact), for creation
([`optionsFromParsed`](../../src/R_interface_bartcore.cpp)) and `setModel` alike;
[`bartcore_setLeafPrior`](../../src/R_interface_bartcore.cpp) takes the named sd;
[`bartcore_scaleDrawnK`](../../src/R_interface_bartcore.cpp) and its registration go. The flat C entries
do not change: the engine retranslates under them.

R ([dbarts.R](../../R/dbarts.R)):
- [`reissueNamedLeafSd`](../../R/dbarts.R) and its six callers go (re-anchors, setData, the install
  shared by re-creation and setState, copy, installTrees, samplePriorPredictive): the engine keeps the
  invariant.
- [`keepDrawnSpread`](../../R/dbarts.R) and [`leafKScale`](../../R/dbarts.R) go; `$setLeafPrior`'s
  same-prior write passes the named sd ([`writeLeafPrior`](../../R/dbarts.R)).
- [`reportLeafPrior`](../../R/dbarts.R) takes an sd-named specification from the model's own
  `leaf.prior`, not from `k.scale` / k, so a fixed sd reads back exactly.
- [`resolveLeafPrior`](../../R/model.R) and the model's encoding are unchanged (settled below); their
  comments and the [`dbartsModel`](../../R/A_class.R) slot comment say what the encoding now means.

## State format and saved objects

- A stored k means the chain's k against the data's scale. For a k-spelled sampler that is what it
  meant. For an sd-spelled sampler the number changes meaning; no format version moves. The registry
  rule at [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) holds both numbers at 1 until the
  first release, no serialized format having shipped (read), and 0.9-34 states are refused already.
- Consequence: a sampler built with a drawn sd on a development build, saved and reloaded on the new
  build, restarts at its spread times rho. A held sd is model and is retranslated on reload, so it is
  unaffected. No migration (Open call 5).
- A saved fit carries k and `k.scale` as a pair; the extractor divides one by the other, so old and new
  fits both report the right `leaf.prior.sd` (read). The R model's encoding does not change, so a saved
  model or sampler object reads as before.
- No `.rds` fixture in the tree holds an sd-spelled state (read: the tracked `.rds` files are the
  equivalence baselines and the 0.9-34 classic compare).

## Composition with the other plans

The setState slice (setstate-force-update.md on wt/install-surface-plan, 26406631, read):
- Its B5 builds dec-B384 as `k = fs.k / r` in `Chain::setState` behind a `kFollowsUnits` flag that R
  sets only under the k spelling, because today's sd `k.scale` ignores the response. Under this rule the
  division is unconditional and sits in `convertStateUnits` (step 7), covering the warm start too: B5,
  its flag, its bridge argument (arity 4 to 5, not 6) and its R argument go.
- Its spread test "under an sd-named drawn prior ... getK unchanged" reverses (k divided under every
  spelling), and its mutant "k divided under an sd-named prior" becomes "k not divided".
- Its H4 (the spread of a state stored under another prior) closes: dec-B401 keeps k, and k now keeps
  the spread.
- Part A (missingness first seen) does not touch the leaf scale.

response-scale-rows (wt/rsr-recheck, faf8ca7d, read): its step 5's re-derivation after a count fit's
creation under a mask "restates a named sd" from R; under this plan the engine does, and the R write
goes. Its tinytest's named-sd arm and its `k.scale` literal then read the data's scale over the rows in,
under both spellings. Its flat-entry status (dec-B405) and dec-B364's keep need nothing here: a kept
scale moves no transform.

Order (Open call 4): this plan first, after gp-copy-continuation; then setState Parts A and B with B5
struck; then response-scale-rows B. Whichever lands second of this and response-scale-rows rereads the
other's named-sd arms. Serial with the logistic scale move if one is built (dec-B404), which would write
k in the same file. The wording sweep queued last stays last.

## Constraints

- Frozen: every k-spelled fit's bits, in every operation but an install or warm start across response
  units with a drawn k; the R model's encoding; the C API's signatures, structs and enums.
- Out of scope: several forests, multinomial and hurdle fits (a named sd is refused there, unchanged);
  the variance forest's leaf prior; the probit rescaling step; a setter for k (dec-B393: "For that, you'd
  need a setK() function.", not asked for); what a re-anchor does to the live leaves (they lag, as
  built).
- The setState slice's B5 is not built twice: whichever lands second takes the other's form
  (Composition).

## Steps

1. Engine steps 1 to 3 and 6 (creation and setModel), tests/cpp for them; the emulation identity:
   `sd = s` against `k = k.scale / s` within 1e-12 over 1000 draws, and bitwise where rho is 1.
2. Engine step 4 (the transform-moving sites) and the invariant check.
3. Engine step 5 and the bridge's `setLeafPrior`; engine step 8 and the bridge removals.
4. Engine step 7 (k with the units), its tests/cpp check; the warm start's per Open call 3.
5. R: the removals, `writeLeafPrior`, `reportLeafPrior`, comments.
6. tinytest: the new file and the changed pins (Tests).
7. Help: [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) (`getLeafPrior`'s `k.scale`; the
   `setLeafPrior` paragraph, k untouched because `k.scale` never moves; the named-sd paragraph, a drawn
   value lagging; `setState`, copy and `installTrees` keeping the spread across units);
   [dbartsPriors.Rd](../../man/dbartsPriors.Rd) (the named sd translated to k against the table's
   `k.scale`; its "absolute" paragraph); [bartBT.Rd](../../man/bartBT.Rd) and [bart.Rd](../../man/bart.Rd)
   where they say what a fit's k is relative to; the `k` field in dbarts.h. NEWS: nothing (the sd
   spelling is new in 1.0-0; dec-B384's line rides the slice that carries it).
8. Docs: [leaf-scale-rules.md](../design/leaf-scale-rules.md) Status (decided, with this plan) and a
   short closing section saying which view was taken; [state-not-model.md](../design/state-not-model.md)
   where it says what k a state carries.
9. Gates (below), then the records at landing: the ledger entry quoting the ruling, this plan's Status
   and Landing, TODO.

## Tests

tests/cpp ([test_sampler.cpp](../../tests/cpp/test_sampler.cpp), [test_state.cpp](../../tests/cpp/test_state.cpp),
[test_facade.cpp](../../tests/cpp/test_facade.cpp)):
- Creation: a fixed named sd's k equals the reader's `k.scale` over it bitwise, the internal leaf scale
  equals a k-spelled sampler's; a drawn one's chi scale and start likewise.
- Invariant: after every transform move (both re-anchor setters, applyNewData, setAnchor both arms, an
  install that restores a scale), a fixed k equals `k.scale` over the named sd, a drawn chi scale
  likewise, and a drawn k is the value before the move.
- setModel keeps a drawn k across a change of named sd and of spelling.
- The named-sd writer: fixed and drawn, a drawn k untouched, an equal write moving no bit.
- convertStateUnits: k divided by the ratio; NaN stays NaN; ratio 1 bitwise.
- The facade list: the renamed virtual, scaleDrawnK gone.

tinytest, new `test-k-internal.R`: one block per row of the operation table, asserting k,
`k.scale` and the spread (to 1e-12) on the probe's fixture; the binary identity (probit `sd =
invchi(1.5, 1.5)` and `k = chi(1.5, 2)` identical); the fit readers listed above; a write-back of
`getLeafPrior()$leaf.prior` moving no bit under both sd forms.

tinytest, changed (counts of k readers in sd-spelled files, ran: grep):
[test-calibration-midchain.R](../../inst/tinytest/test-calibration-midchain.R) (40, the spelling-switch
and named-reading pins among them), [test-fit-stores-k.R](../../inst/tinytest/test-fit-stores-k.R) (16),
[test-calibration-creation.R](../../inst/tinytest/test-calibration-creation.R) (11),
[test-state-not-model.R](../../inst/tinytest/test-state-not-model.R),
[test-leaf-prior-k-or-sd.R](../../inst/tinytest/test-leaf-prior-k-or-sd.R) (its header comment; its
translation pins hold, the R encoding being unchanged),
[test-nbinom.R](../../inst/tinytest/test-nbinom.R),
[test-embedding-recipes.R](../../inst/tinytest/test-embedding-recipes.R); the full suite finds the rest.

Mutants the reviewer runs, each failing a test: the translation skipped at each transform-moving site in
turn; a drawn k retranslated at a re-anchor (it must lag); a drawn sd started at k 2; setModel
re-expressing k against a ratio (today's setLeafPrior arithmetic); convertStateUnits not dividing k, or
dividing a recipient's fixed k; `k.scale` reported as twice the named sd; reportLeafPrior reading a fixed
sd as `k.scale` / k; the named-sd writer moving a drawn k.

## Gates

On the slice's tip and its own library, independently of the implementer, shifting class:
- tests/cpp plain and under `-fsanitize=address,undefined`; R-loaded ASAN over the touched test files.
- The full tinytest suite (`at_home = TRUE`); stan4bart's suite against the build.
- Reference build, `--preclean`: the equivalence trio `compare --bitwise` (equivalence.R also
  `--strict-coverage`) against [MANIFEST](../../benchmarks/baselines/MANIFEST)'s current files, and the
  four arm64 snapshot files. Expected: all identical. No equivalence scenario or snapshot file writes an
  sd spelling, a leaf-prior writer, `setModel`, a state install or a warm start; the only re-anchors are
  `setData` scenarios under the k spelling (ran: grep). Any move is a stop.
- Every exact gate in quick mode ([exact-gates.yaml](../../.github/workflows/exact-gates.yaml)): what a
  fit carries changes value. [aft-exact.R](../../benchmarks/R/aft-exact.R),
  [t-exact.R](../../benchmarks/R/t-exact.R) and [logistic-reference.R](../../benchmarks/R/logistic-reference.R)
  write a fixed sd and move by rounding; each must pass. No baseline moves, so no oracle is needed; if a
  gate fails, the emulation identity (step 1) is the oracle: today's build at `k = k.scale / s` against
  the new build at `sd = s`.
- `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.
- Not hot-path: the translation runs at a mutation, never per sweep; no bench-sampler compare.

## Budget and stops

Planned ~950 lines. Engine slices on this surface ran 1.1 to 2.9 times their plans; the last ran 1.76.
Forecast at 1.5 to 2 times: 1425 to 1900. Stop at 2400 and report, without working around, when:
- the diff passes 2400 lines;
- any equivalence scenario or snapshot moves;
- a k-spelled fit moves a bit in any operation but an install or warm start across units;
- an sd-spelled fit departs from the emulation by more than 1e-10 within 1000 draws;
- a transform-moving site is found that the invariant check cannot reach;
- a reading needs a state-format change, or a test needs a call no ruling, settled call or open call
  here makes.

## Settled in planning

- The R model keeps its encoding (`prior.scale` 2 s with k at 2 or `chi(df, 2)`); the bridge divides it
  out. Saved models read as before, the translation pins hold, and nothing R-side moves but comments.
  Re-encoding the slot as s would touch every reader of it for no change a user sees.
- A drawn sd's start is the translation of today's: the chain starts at the named sd. A k-spelled
  `chi()` still starts at 2.
- `getLeafPrior` reports a named sd from the model, exactly; the engine's fixed k is `k.scale` / s,
  which need not divide back to s's bits.

## Open calls

Each with a recommendation, for the maintainer unless marked.

1. Where the translation lives (the orchestrator's). The engine holds the named sd and retranslates
   wherever the transform moves (recommended): the flat C re-anchor, which stan4bart calls every warm-up
   step, keeps a named sd as R's re-anchor does, and one invariant replaces six R write-backs. The
   alternative, R and the bridge translating after the engine computes the scale, leaves the engine
   untouched but needs a two-step creation and leaves the flat path's gap (a named sd stretching with
   each warm-up re-anchor), against dec-A105's every call that re-anchors.
2. A drawn sd's current value at a re-anchor. Recommended: it lags with the leaves until the next draw,
   the prior retranslated at once, as dec-B414 states ("The parameter itself of course lags, but that's
   to be expected."). It is what a drawn k does and did on 0.9-34, where every `setResponse` re-anchored
   and k stayed; the re-anchor changes the target and the next draw comes from the named prior in
   response units. The alternative re-expresses k at each re-anchor so the value stays, which is today's
   built behaviour and needs a ratio write at every re-anchor, the one move this rule removes.
3. The warm start onto another response spread. Recommended: k is converted with the leaves, as a
   restored state's is (dec-B384); one line serves both, since they share the conversion. A user
   warm-starting onto a rescaled response means the same function and the same spread in response
   units; as built, the spread restarts as many times larger as the response is wider (0.8992 to 2.6977
   on 3 y, ran). dec-B343's "as built" answered which values transfer, not their units. 0.9-34 had no
   warm start. The alternative keeps k as stored, with the help saying so.
4. Order. Recommended: this plan first, after gp-copy-continuation; the setState slice's B5 then shrinks
   to nothing and H4 closes, and response-scale-rows writes its named-sd tests once against the data's
   `k.scale`. The alternative, setState Part B first, builds the `kFollowsUnits` flag this plan then
   deletes (about 30 more lines here).
5. Saved development-build states of an sd-spelled sampler with a drawn sd (the orchestrator's).
   Recommended: no migration, per the state registry's pre-release rule; such a state restarts at its
   spread times rho. The alternative, a state attribute recording the reference scale k was drawn
   against, costs about 60 lines and serves only states written by builds that never shipped.

## Evidence

Ran on 4cf51b11 (private library scratch/libs/kint; probes and outputs in scratch/kint/): today's
behaviour per operation (01-today.R, 02-today-more.R), the emulation of the new sd spelling
(03-emulate.R), and the greps for sd-spelled tests, baselines, snapshot files and exact gates. Read:
every engine, bridge, R and header claim cited by symbol; stan4bart's bartcore branch at 963956b, its
bart_args forwarding and its flat-API calls; bartCause, treatSens and bairrtt for their k use; the setState
plan at 26406631 and the response-scale-rows plan at faf8ca7d; the rulings in the table. New values in
the operation table are arithmetic on the run values, not runs.
