# k-internal-parameter: k is the chain's parameter on the internal scale, and the sd spelling is translated

Status: PLANNED 2026-10-09 for dec-B414, revised from its blind critique (scratch/kcrit2/critique.md,
build after corrections); for one ruling on the whole. Nothing is built.

agent: opus implementer, one; one opus reviewer who runs the mutants under Tests.
rng: SHIFTING for every fit made and moved through R: the posterior does not change. NEUTRAL, bit for
bit, for every fit whose leaf prior is written with k, in every operation but a state install or warm
start across response units with a drawn k. Draws move by rounding only for fits written with an sd
(ran: within 2.3e-14 over 1000 draws on gaussian; the critique ran every family and leaf model that
takes an sd, within 4.4e-13 over 500 draws, the k ratio exact). Draws move in trajectory, the posterior
unchanged, after `setModel` with another leaf prior, an install across a change of prior into or out
of the sd spelling, a re-anchor under a drawn sd, an `xbart` sd grid's warm sweep into a drawn cell, and
an install or warm start across response units with a drawn k. POSTERIOR-CHANGING on one path: a
sampler with a named sd re-anchored through the flat C entries, which stan4bart does every warm-up step
when `bart_args` names an sd. Today nothing restores the sd there, so sampling runs under the sd
stretched by the last warm-up re-anchor; after this plan it runs under the named sd, a fix toward
dec-A105. No exact gate reaches that path; it is gated by a new arm of the C API test, whose compiled
consumer calls the flat `setOffset`, and by stan4bart's suite (Gates).
window: before the merge to main; engine slices stay serial. Its place in the engine queue is Open
call 4.
budget: planned ~1020 lines (engine ~160, bridge ~40, R ~80 net of removals, tests/cpp ~210, tinytest
~390, help ~70, docs and a benchmark comment ~70). Forecast and stops:
[Budget and stops](#budget-and-stops).

## Summary for the maintainer

Every chain holds one number for the leaf spread, k, and the spread is a reference scale divided by k.
Today that reference scale depends on how the prior was written: with k it is the data's own scale
(half the response range for a continuous response, 3 for probit, about 5.44 for logistic); with an sd it
is twice the named sd, so that k sits at 2. Because the yardstick moves with the spelling, keeping k and
keeping the spread are the same only until the spelling changes, and the rulings so far keep one in some
operations and the other in others.

Under the rule you asked for, the yardstick never moves: it is always the data's scale, and k is the
chain's parameter against it. Writing the prior with an sd becomes a translation: an sd of s is k equal
to the data's scale over s, and a scaled inverse chi prior on the sd with scale s is a chi prior on k
with scale the data's scale over s, the reciprocal pair the help already states. Since nothing else
redefines the yardstick, keeping k is keeping the spread, and the two rulings that seemed to disagree
(a change of prior keeps the spread; a restored state keeps its k) now say the same thing.

What changes for a user:

- A prior written with k: nothing. Every such fit is bit for bit what it is today, except a state
  restored onto a response of another spread, which your earlier ruling to keep the spread already
  covers and which is not built yet, and the warm start onto such a response (open call 3).
- A prior written with an sd: the same model and the same posterior. Each k is today's times one
  constant, and the fits agree with today's to about 1e-13 or better on every family and leaf model.
  The k a fit reports is a different number (the data's scale over the spread, not twice the named sd
  over it), and the reported scale k is measured against becomes the data's. A held sd is still
  reported exactly as you named it.
- setModel given another leaf prior now keeps the spread, as setLeafPrior does. Today it keeps k, and
  the spread jumps when the spelling or a named sd changes. xbart's sweep over an sd grid goes through
  the same path, so a drawn cell now starts where the previous cell left the spread (as a k grid's
  cells already do), not at its own named scale.
- A state stored before a change of prior, restored after it, keeps its spread. Today it jumps when the
  change crossed between the k and sd spellings or changed a named sd.
- After the response is re-derived (a re-anchor), a drawn sd's current value stretches with the fit
  until the next draw, as a drawn k's does, while its prior stays in response units. Today the current
  value stays put while the fit stretches (open call 2).
- stan4bart with an sd named in its BART arguments: its warm-up re-derives the scale through dbarts's C
  interface, which today lets a named sd stretch with the response, so its sampling ran under a
  different prior than the one named. After this change it runs under the named sd: a posterior change,
  and a fix.

Rulings reversed, in part:

- dec-B201's definition of the reference scale under the sd spelling (twice the named sd), with the
  agents' matching call dec-A124. Your words there ("I want k to be the engine's k always", a way to
  convert between k and sd) still hold.
- In dec-B393, the register's clause that the reference scale is twice the named value under the sd
  spelling. Your words there ("reset the k.scale", "changing the prior" is not "changing the state")
  still hold, more simply than before.
- dec-A105, for a drawn sd's current value at a re-anchor: it lags with the fit instead of staying in
  response units, if you take open call 2's recommendation. The named sd itself still stays in response
  units.
- dec-A146's agents' call that a warm start installs k as stored across a change of response spread, if
  you take open call 3's recommendation.

Open calls for you (numbered as in the plan), each with a recommendation:

2. A drawn sd's current value at a re-anchor. Recommended: it lags with the fit until the next draw, as
   your words in dec-B414 have it ("The parameter itself of course lags") and as a drawn k does. The
   cost: a drawn residual sd stays put in response units across the same call (measured: 0.3933 before
   and after), so the drawn leaf sd and the drawn residual sd part there; the leaf sd follows the leaves,
   which are kept on the internal scale, and the residual sd follows the residuals, which are not. The
   note showed you the k view keeping this value, as built; the alternative does that, and then a drawn
   k and a drawn sd behave differently at a re-anchor.
3. A warm start onto a response of another spread. Recommended: k is converted with the leaves, so the
   spread is kept, as you ruled for a restored state. As built the spread restarts as many times larger
   as the response is wider (0.90 to 2.70 on a response three times as wide).
4. Order among the engine slices. Recommended: gp-copy (in review), then the width-weighted cut rule
   (ruled, ready), then this plan, then setState's two parts, then response-scale-rows, then the
   logistic scale move and whatever the sampler survey brings. This plan goes before setState and
   response-scale-rows so neither builds or tests what it removes, and before any new k sampler so that
   sampler is written once against k's final meaning.

Two further calls, 1 and 5, are implementation (where the translation lives: the engine; saved development-build
states: no migration). The work is about 1020 lines planned, 1530 to 2040 at the usual overrun. No
recorded baseline or snapshot is expected to move.

## Goal

On the internal scale the leaf prior's reference scale, `k.scale`, is a constant of the family and the
response transform, never redefined by a prior's spelling, and k is each chain's parameter against it.
`normal(sd = s)` is `k = k.scale / s`, and `normal(sd = invchi(df, s))` is `k = chi(df, k.scale / s)`,
retranslated whenever the transform moves so the named sd keeps its value and distribution in response
units. Only an operation that rescales the leaves rescales k; a re-anchor leaves the leaves and a drawn k
as they were, a held k being retranslated; a change of prior touches neither.

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
- A re-creation (a reload, a copy, `setState` on a dead pointer) creates the sampler with the install
  to follow, so [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp) leaves the chains and
  [`Chain::setState`](../../src/bartcore/chain.hpp)'s own scale restore moves them (the critique, read).
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

The critique widened it (scratch/kcrit2/01-families.R, 03-monotone.R, ran): fixed and drawn pairs over
500 draws on gaussian with two chains, probit, logistic, nbinom, ordinal, aft, Student-t, linear and gp
leaves, a variance forest and monotone (fixed): every pair within 4.4e-13 (logistic 1.4e-13, gp
4.4e-13), the k ratio exact to 10 digits, so the probit rescaling step is invariant too; 5000 draws on
gaussian, 3.9e-14. It also ran (02-ops.R) a drawn sigma across `setResponse(3 y, updateScale = TRUE)`
and across `setState` onto 3 y: 0.3933 before and after, in response units, while the leaf spread
stretched (0.8992 to 2.6977); and xbart's warm path, `setModel` from a fixed `sd = 0.25` into
`invchi(3, 1)` (spread 0.25 to 1.0) and from a chain at 0.9292 under `invchi(3, 0.5)` into
`invchi(3, 2)` (to 3.7169).

## The rulings under this rule

For one ruling on the whole. "Unchanged": the rule and what it builds stand. "Restated": the outcome
stands and the mechanism or wording changes. "Reversed": part of the ruling no longer holds; the part is
named. Built behaviour that changes without reversing a ruling is in the next section. Words in
quotation marks are the maintainer's; "Register:" marks the decision register's own wording, which is
not.

| ruling | what it says | under this rule |
|---|---|---|
| dec-B201 | Register: `k.scale` is the value k is relative to; under an sd-named prior twice the sd or `invchi()` scale named, set so that the engine's k sits at 2. The maintainer: "I want k to be the engine's k always."; "Ultimately, the exact value isn't itself important. The only thing that really matter is that users have a way to convert between k and standard deviation." | Reversed in the register's definition: `k.scale` is the data's scale under every spelling. The maintainer's words stand: k is the engine's k, and spread = `k.scale` / k converts. |
| dec-A124 | (agents' call) Register: `k.scale` the data's under k and twice the named sd under sd, so the data's scale is not readable under an sd-named prior | Reversed with dec-B201: the data's scale is readable on every fit. |
| dec-B141 | `getK` reports the engine's k whatever the spelling | Unchanged; under an sd spelling the number is `k.scale` / spread. |
| dec-B192 | Register: a fit carries the k the sampler recorded, on an sd-named fit relative to twice the named scale. The maintainer: "The samples of `k` (if there are any) should be stored as `k`. End-users should use `extract`" | Restated: relative to the data's scale. The maintainer's words stand. |
| dec-B193, dec-B194 | extract answers k and `leaf.prior.sd` on every fit; one number for a held value | Unchanged; `leaf.prior.sd` values unchanged, k values move under an sd spelling. |
| dec-B376 | Register: no fit stores the leaf sd; `fit$fixed$k` for a held k | Unchanged; a held sd's `fit$fixed$k` becomes `k.scale` / s, and extract and summary read the held sd itself from the recorded prior, exactly. |
| dec-A105 | Register: named by k (relative) or sd (absolute); a named sd stays in response units across every call that re-anchors the response scale, re-applied after it, as a fixed residual sd and the variance forest's prior already do. The maintainer, on re-anchoring: "OK, proceed using it." | Reversed for a drawn sd's value in force at a re-anchor, pending Open call 2: it lags with the leaves until the next draw ("The parameter itself of course lags", dec-B414), where it stays today; the leaf-scale note's k-view column showed it kept, as built. Restated for the named sd, fixed or the scale of its `invchi()`: retranslated at every re-anchor, keeping its value and distribution in response units, now on the flat C path too. |
| dec-B356 | A switch between spellings of a drawn prior keeps the spread. Register: k becoming k_old x k.scale_new / k.scale_old. The maintainer: "B. Keep the spread." | Restated: `k.scale` never differs, so k is untouched and the spread is kept. |
| dec-B369 | a changed `invchi()` scale keeps the spread at the call ("OK, then A. Keep the spread at the call.") | Restated as dec-B356. |
| dec-B392 | a fixed sd turned drawn keeps the sd in force ("They're trying to actively set the sd of the leaf prior, so both of their values should be interpretted as such.") | Restated: the fixed sd s is k = `k.scale` / s, and the drawn prior keeps that k. |
| dec-B393 | Every switch keeps the sd into a drawn prior; a stated value is literal. Register: `k.scale` always belongs to the new prior's spelling (the data's scale under k, twice the named value under sd); changing the prior never moves the state; no setter for k exists. The maintainer: "if a person installs a fixed value of k, they clearly mean it to be interpretted literally so that should also reset the k.scale", "Treat drawn k priors and fixed k priors the same - reset the k scales. The person is trying to put a distribution on k itself and otherwise the interpretation is off.", "I don't see changing the prior as changing the state / the active parameter. For that, you'd need a setK() function.", "Yes, keep until the next draw." | Rule stands. Reversed: only the register's "twice the named value under sd". The maintainer's words hold: every prior's `k.scale` is the data's, so a k, fixed or drawn, is read literally against it, and neither k nor the spread moves at a change of prior. |
| dec-B401 | "The k itself should literally transfer - a parameter is a parameter, regardless of the prior. It may be a bad fit, but so be it." | Unchanged; the spread is now kept too, since no prior moves `k.scale`. The setState plan's H4 closes. |
| dec-B384, dec-B407 | A state installed onto another response spread keeps its spread. Register: k becoming k_state x k.scale_recipient / k.scale_state. The maintainer: "3. Use your recommendation." (dec-B384); "Keep it for scale changes. If it helps at all, I guess we can think of `k` as just the internal representation and the parameter as the sd." (dec-B407) | Restated: k is divided by the units ratio with which the leaves are converted, under every spelling (today the sd spelling's `k.scale` ignores the response, so the register's formula left k alone there). The two readings in dec-B407 now agree. |
| dec-A146 | (agents' call, one of those "not put" to the maintainer) Register: a warm start follows the setState rule; a drawn k installs as stored | Its install-as-stored part, revised by dec-B384 for setState, is reversed for the warm start too, pending Open call 3. |
| dec-B343 | "A. The donor gives its trees, and sigma and k only where the new fit draws them (as built)." | Unchanged in what transfers; k is converted with the leaves when the donor is in other units (Open call 3). |
| dec-B200 | a state in other units is converted on install ("Convert on install.") | Unchanged; k joins the leaves in the conversion. |
| dec-B195, dec-B196 | A state carries no model. The maintainer: "Model when fixed, state when drawn." | Unchanged. A held named sd's k is model, retranslated with the data's scale. |
| dec-B254 | Register: `setModel` changes parameters, the leaf scale as k or sd among them, never structure. The maintainer: "Parameters yes, structure no." | Unchanged in scope; it now keeps the spread across a change of leaf prior (the note's defect a), xbart's sd grid included. |
| dec-B396 | Register: `copy()` stores the current state when none is stored | Unchanged. |
| dec-A107 | Register: `setLeafPrior` changes only the spread or its prior | Unchanged. |
| dec-A121 | (agents' call) the engine's k carried across a switch of spelling | Already revised by dec-B356; moot. |
| dec-A187 | (agents' call) `setLeafPrior` keeps the spread through `scaleDrawnK`, which state-install-keeps-spread would reuse | Superseded: no write moves `k.scale`, and the install's conversion divides k itself; the primitive goes. |
| dec-B371, dec-A190 | the probit rescaling step divides k by the drawn factor | Unchanged: it acts on the internal k and the leaves together. |
| dec-A13 | no cap on a drawn k; `chi(df, Inf)` improper | Unchanged; `invchi(df, 0)` is still `chi(df, Inf)`. |
| dec-B331, dec-B330 | a re-derived scale converts saved draws, the live leaves kept; a gp sampler holding saved draws refuses | Unchanged; "both lag" is how the live leaves already behave. |
| dec-B361 | prior draws under the prior the sampler runs under | Unchanged; the private sampler's named sd is translated by the engine. |
| dec-B362, dec-B364, dec-B405 | several forests refuse a re-derivation; a kept scale (fewer than two values); the status through the C API | Unchanged. A kept scale moves no transform, so nothing is retranslated. |
| dec-B122 | a variance forest recalibrates on a re-anchor | Unchanged (it has no k). |
| dec-B142, dec-B253, dec-B275 | the multi-forest writer and a forest's sd | Unchanged: a named leaf-prior sd is refused on several forests, multinomial and hurdle fits. |

Reversed, in sum: the register's definition of `k.scale` under the sd spelling in dec-B201 (with
dec-A124's); the register's "twice the named value under sd" in dec-B393; dec-A105 for a drawn sd's
value in force at a re-anchor (pending Open call 2); and dec-A146's as-stored warm start across response
units (pending Open call 3). Superseded: dec-A187's primitive. Every other ruling stands, several with
their mechanism restated; every quoted word of the maintainer's stands.

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
| `xbart` sd grid, the warm sweep into a drawn cell (through `setModel`) | fixed `sd = 0.25` into `invchi(3, 1)`: spread 0.25 to 1.0; a chain at 0.9292 under `invchi(3, 0.5)` into `invchi(3, 2)`: to 3.7169 (ran by the critique) | the previous cell's spread kept: 0.25; 0.9292, as a k grid's `chi()` cells already keep it | every grid with a drawn sd cell after the first; the first cell still starts at its named scale |
| install of a state from the same sampler and prior | as stored, bitwise | unchanged | none |
| install across a change of prior (store, setLeafPrior, restore; or another sampler's state) | k as stored; spread 0.9132 | k as stored; spread 0.8992 | where the k.scale differed: one side sd-spelled |
| install into a sampler holding k fixed | the sampler's | unchanged | none |
| `setState`, `copy()`, reload onto another response spread | k as stored; k spelling 0.8992 to 2.6977; sd spelling kept | k divided by the units ratio: 6.5703, spread 0.8992, either spelling (dec-B384) | every drawn k installed across units |
| warm start across a change of prior | k as stored; 0.8992 to 0.9132 | k as stored; spread kept | where one side is sd-spelled |
| warm start onto another response spread | k as stored; 0.8992 to 2.6977 | k divided by the units ratio; spread kept (Open call 3) | every drawn k warm-started across units |
| re-anchor, k spelling | k kept; spread stretches with the fit | unchanged | none |
| re-anchor, fixed sd `0.5` on 3 y | k 2, `k.scale` 1, spread 0.5 | k retranslated 11.8164, `k.scale` 5.9082, spread 0.5 | rounding |
| re-anchor, drawn `sd = invchi(3, 1)` on 3 y | k kept, spread 0.8015 held; next draw from invchi(3, 1) | k kept (2.4572 at the new numbers), spread lags to 2.4044 with the fit; prior retranslated to chi(3, 5.9082), the next draw from invchi(3, 1) in response units | the sweep after each such re-anchor |
| re-anchor through the flat C entries, named sd | not retranslated (read): the sd stretches with the response, and sampling after warm-up runs under the last stretch | retranslated, as through R | POSTERIOR: stan4bart with a named sd in `bart_args` now samples under the named sd (a fix toward dec-A105) |
| reload, `copy()` or `setState` on a dead pointer after a re-anchor, named sd | R writes the named sd back after the install (read) | the engine retranslates at the install's scale restore | rounding |
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
  values. A drawn spread is `k.scale` / k from the fit's own pair, so a fit saved before the change still
  reads correctly (read: ["leaf.prior.sd"](../../R/diagnostics.R)), equal to today's to rounding. A held
  sd is read from the recorded prior, `fit$leaf.prior$leaf.prior`, where a single forest names one
  ([`extractParameter`](../../R/generics.R), which summary's line of held values,
  [`fixedSummaryValues`](../../R/diagnostics.R), calls): `k.scale` / (`k.scale` / s) misses s by one
  ulp for 5 to 14 percent of scales (the critique ran 10000 random scales at s = 0.7, 1.3, 0.1), where
  today 2 s / 2 is exact.
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
4. It runs wherever a chain's response transform changes, beside the variance forest's recalibration
   that already sits at each such site: the `updateScale` arms of
   [`Chain::setResponse`](../../src/bartcore/chain.hpp) and [`Chain::setOffset`](../../src/bartcore/chain.hpp);
   [`applyNewData`](../../src/bartcore/chain.hpp) (setData); and every restore of the response's scale,
   in [`moveScale`](../../src/bartcore/chain.hpp) (the moving arm of
   [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp)),
   [`Chain::installForest`](../../src/bartcore/chain.hpp) and
   [`Chain::setState`](../../src/bartcore/chain.hpp). The last is the re-creation path: a reload, a copy
   and `setState` on a dead pointer install with `setAnchor` leaving the chains, so the install's own
   restore moves them. One Chain routine serves all six; the reviewer greps chain.hpp for every call that
   changes the response's transform (its setResponse and setOffset at `updateScale`, setData,
   restoreScale) and finds the routine after each. tests/cpp holds the invariant after each site (Tests).
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
- [`extractParameter`](../../R/generics.R): a held leaf-prior sd on a single-forest fit whose recorded
  prior names a numeric sd is that sd; summary's held line follows through it.
- [`resolveLeafPrior`](../../R/model.R) and the model's encoding are unchanged (settled below); their
  comments and the [`dbartsModel`](../../R/A_class.R) slot comment say what the encoding now means.
- xbart needs no code change: its warm sweep calls [`bartcoreSetModel`](../../R/xbart.R) per cell, which
  takes the engine's new `setModel` rule.

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

Order: Open call 4. The serial engine queue as it stands (the orchestrator's brief, read):
gp-copy-continuation in review; default-rule-per-gap, width-weighted gaps, planned and ruled to build
(dec-B406, dec-B409); this plan; setstate-force-update, Part A then Part B, waiting on its held questions;
response-scale-rows slice B; and, only if their studies call for them, dec-B404's logistic scale move
(if a drawn k is kept) and samplers from dec-B415's survey of k's slow tail, either of which would write
k's update in chain.hpp. The wording sweep queued last stays last.

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
5. R: the removals, `writeLeafPrior`, `reportLeafPrior`, `extractParameter`'s held sd, comments.
6. tinytest: the new file and the changed pins (Tests).
7. Help: [xbart.Rd](../../man/xbart.Rd) (its `sd` item: "invchi(df, c) starts its chain at the spread
   c" holds for the first cell; each later cell of the warm sweep starts at the previous cell's spread,
   whichever spelling); [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd) (`getLeafPrior`'s `k.scale`; the
   `setLeafPrior` paragraph, k untouched because `k.scale` never moves; the named-sd paragraph, a drawn
   value lagging; `setState`, copy and `installTrees` keeping the spread across units);
   [dbartsPriors.Rd](../../man/dbartsPriors.Rd) (the named sd translated to k against the table's
   `k.scale`; its "absolute" paragraph); [bartBT.Rd](../../man/bartBT.Rd) and [bart.Rd](../../man/bart.Rd)
   where they say what a fit's k is relative to; the `k` field in dbarts.h. NEWS: nothing (the sd
   spelling is new in 1.0-0; dec-B384's line rides the slice that carries it).
8. Docs and a benchmark comment: [composition-matrix.R](../../benchmarks/R/composition-matrix.R)'s
   "calibration" probe comment, which describes the reference k of 2, restated (not a gate);
   [leaf-scale-rules.md](../design/leaf-scale-rules.md) Status (decided, with this plan) and a
   short closing section saying which view was taken; [state-not-model.md](../design/state-not-model.md)
   where it says what k a state carries.
9. Gates (below), then the records at landing: the ledger entry quoting the ruling, this plan's Status
   and Landing, TODO.

## Tests

tests/cpp ([test_sampler.cpp](../../tests/cpp/test_sampler.cpp), [test_state.cpp](../../tests/cpp/test_state.cpp),
[test_facade.cpp](../../tests/cpp/test_facade.cpp)):
- Creation: a fixed named sd's k equals the reader's `k.scale` over it bitwise, the internal leaf scale
  equals a k-spelled sampler's; a drawn one's chi scale and start likewise.
- Invariant, one check per site of step 4: after each re-anchor setter at `updateScale`, applyNewData,
  setAnchor's moving arm, installForest's restore and Chain::setState's restore (a state in other units
  into a sampler created with the install to follow), a fixed k equals `k.scale` over the named sd, a
  drawn chi scale likewise, and a drawn k is the value before the move.
- setModel keeps a drawn k across a change of named sd and of spelling.
- The named-sd writer: fixed and drawn, a drawn k untouched, an equal write moving no bit.
- convertStateUnits: k divided by the ratio; NaN stays NaN; ratio 1 bitwise.
- The facade list: the renamed virtual, scaleDrawnK gone.

tinytest, new `test-k-internal.R`: one block per row of the operation table, asserting k,
`k.scale` and the spread (to 1e-12) on the probe's fixture; the binary identity (probit `sd =
invchi(1.5, 1.5)` and `k = chi(1.5, 2)` identical); the fit readers listed above; a write-back of
`getLeafPrior()$leaf.prior` moving no bit under both sd forms. And:
- Re-creation after a re-anchor: a sampler with a held and with a drawn named sd, re-anchored, then
  copied, saved and reloaded with `saveRDS`, and given `setState` on a dead pointer: the held spread is s
  to 1e-12 and the drawn chi scale is `k.scale` / s, each way.
- xbart's path: [`bartcoreSetModel`](../../R/xbart.R) from a fixed `sd = 0.25` into `invchi(3, 1)` keeps
  the spread at 0.25, and from a drawn chain into `invchi(3, 2)` keeps its spread; an `xbart` grid
  mixing a fixed and a drawn cell of different scales runs and repeats under a seed. (The existing grid
  tests pair each drawn cell with a fixed cell of the same s and cannot see the change.)
- A held sd read exactly: for s in 0.7, 1.3 and 0.1 on several responses, `extract(fit,
  "leaf.prior.sd")` and summary's held line are identical to s.

tinytest, [test-capi.R](../../inst/tinytest/test-capi.R), a new arm: a sampler with a held and with a
drawn named sd, re-anchored through the compiled consumer's flat `setOffset` with `updateScale` on an
offset that moves the range; the held spread stays s to 1e-12, the drawn chi scale is `k.scale` / s, a
drawn k is unchanged. It fails today. The file skips where the consumer cannot be compiled, as now.

tinytest, changed (counts of k readers in sd-spelled files, ran: grep):
[test-calibration-midchain.R](../../inst/tinytest/test-calibration-midchain.R) (40, the spelling-switch
and named-reading pins among them), [test-fit-stores-k.R](../../inst/tinytest/test-fit-stores-k.R) (16),
[test-calibration-creation.R](../../inst/tinytest/test-calibration-creation.R) (11),
[test-state-not-model.R](../../inst/tinytest/test-state-not-model.R),
[test-leaf-prior-k-or-sd.R](../../inst/tinytest/test-leaf-prior-k-or-sd.R) (its header comment; its
translation pins hold, the R encoding being unchanged),
[test-nbinom.R](../../inst/tinytest/test-nbinom.R),
[test-embedding-recipes.R](../../inst/tinytest/test-embedding-recipes.R); the full suite finds the rest.

Mutants the reviewer runs, each failing a test: the translation skipped at each of step 4's six sites in
turn, Chain::setState's restore among them; a drawn k retranslated at a re-anchor (it must lag); a drawn sd started at k 2; setModel
re-expressing k against a ratio (today's setLeafPrior arithmetic); convertStateUnits not dividing k, or
dividing a recipient's fixed k; `k.scale` reported as twice the named sd; reportLeafPrior reading a fixed
sd as `k.scale` / k; `extractParameter` computing a held sd as `k.scale` over the held k; the named-sd
writer moving a drawn k.

## Gates

On the slice's tip and its own library, independently of the implementer, shifting class:
- tests/cpp plain and under `-fsanitize=address,undefined`; R-loaded ASAN over the touched test files.
- The full tinytest suite (`at_home = TRUE`), test-capi.R's new arm compiled and run (a skip there is a
  stop: it is the only check of the flat path in this package); stan4bart's suite against the build.
  The flat-path posterior change has no exact gate: stan4bart's tests write k (read), so its suite
  checks that nothing else moved, and test-capi.R checks the fix.
- Reference build, `--preclean`: the equivalence trio `compare --bitwise` (equivalence.R also
  `--strict-coverage`) against [MANIFEST](../../benchmarks/baselines/MANIFEST)'s current files, and the
  four arm64 snapshot files. Expected: all identical. No equivalence scenario or snapshot file writes an
  sd spelling, a leaf-prior writer, `setModel`, a state install or a warm start; the only re-anchors are
  `setData` scenarios under the k spelling (ran: grep). Any move is a stop.
- [composition-matrix.R](../../benchmarks/R/composition-matrix.R) writes a fixed sd through
  `setLeafPrior` but is not a CI gate; it is run once, all cells as recorded.
- Every exact gate in quick mode ([exact-gates.yaml](../../.github/workflows/exact-gates.yaml)): what a
  fit carries changes value. [aft-exact.R](../../benchmarks/R/aft-exact.R),
  [t-exact.R](../../benchmarks/R/t-exact.R) and [logistic-reference.R](../../benchmarks/R/logistic-reference.R)
  write a fixed sd and move by rounding; each must pass. No baseline moves, so no oracle is needed; if a
  gate fails, the emulation identity (step 1) is the oracle: today's build at `k = k.scale / s` against
  the new build at `sd = s`.
- `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.
- Not hot-path: the translation runs at a mutation, never per sweep; no bench-sampler compare.

## Budget and stops

Planned ~1020 lines (the first draft's ~950 and the critique's additions: Chain::setState's site, xbart's
help and tests, the C API arm, the held sd's exact read). Engine slices on this surface ran 1.1 to 2.9
times their plans; the last ran 1.76. Forecast at 1.5 to 2 times: 1530 to 2040. Stop at 2550 and report,
without working around, when:
- the diff passes 2550 lines;
- any equivalence scenario or snapshot moves;
- a k-spelled fit moves a bit in any operation but an install or warm start across units;
- an sd-spelled fit departs from the emulation by more than 1e-10 within 1000 draws;
- a site that changes the response's transform is found that step 4's routine cannot be put after;
- test-capi.R's consumer cannot be compiled on the gate machine;
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
   to be expected."). It is what a drawn k does, and what k did on 0.9-34, where every `setResponse`
   re-anchored and k stayed; the leaf sd is the scale of the leaves, which a re-anchor keeps on the
   internal scale, so the two stay consistent (the note's disagreement 4 goes), and the next draw comes
   from the named prior in response units. Against it: a drawn sigma stays in response units across
   the same call (0.3933 before and after `setResponse(3 y, updateScale = TRUE)`, ran by the critique),
   and dec-A105's register text tied a named sd to "a fixed residual sd", so the drawn leaf sd and the
   drawn sigma part at a re-anchor; sigma is held in response units and the leaves are not. The
   leaf-scale note's k-view column showed this value kept, as built; this reading of dec-B414 does not
   keep it. The alternative re-expresses an sd-spelled drawn k at each re-anchor so the value stays
   (about 20 lines, a ratio write the rule otherwise removes), and then a drawn k and a drawn sd behave
   differently at a re-anchor.
3. The warm start onto another response spread. Recommended: k is converted with the leaves, as a
   restored state's is (dec-B384); one line serves both, since they share the conversion, which has
   exactly those two callers and leaves sigma in response units (the critique, read and ran). A user
   warm-starting onto a rescaled response means the same function and the same spread in response
   units; as built, the spread restarts as many times larger as the response is wider (0.8992 to 2.6977
   on 3 y, ran). dec-A146's install-as-stored was a call not put to the maintainer, and dec-B343
   confirmed which values transfer, not their units. 0.9-34 had no warm start. The alternative keeps k
   as stored, with the help saying so.
4. Order in the serial engine queue. Recommended: gp-copy-continuation (in review), then
   default-rule-per-gap, then this plan, then setstate-force-update Parts A and B, then
   response-scale-rows B, then dec-B404's logistic move and dec-B415's samplers if their studies call
   for them. Why: default-rule-per-gap is ruled and ready while this plan waits on a ruling, and it
   touches the cut grids, not the leaf scale, so neither changes the other's tests; if this plan is
   ruled first it may go ahead of it. This plan before setState Part B, so B5's `kFollowsUnits` flag and
   H4 are never built (about 30 lines saved here and more there); Part A is independent and may go
   either side. Before response-scale-rows, so its named-sd arms and `k.scale` literals are written once
   against the data's scale. Before any new k sampler, so that sampler's code and exact gate are written
   once against k's final meaning (under the sd spelling the chi scale it sees becomes `k.scale` / s).
   The alternative, this plan after setState Part B, deletes the flag that slice builds.
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
plan at 26406631 and the response-scale-rows plan at faf8ca7d; the rulings in the table; the
orchestrator's queue for Open call 4. New values in the operation table are arithmetic on the run
values, not runs. The blind critique's runs (scratch/kcrit2/, private library scratch/libs/kcrit2 on the
same code): the probes re-run identically, the emulation over every family and leaf model, sigma across a
re-anchor and an install, xbart's warm path, and the held sd's one-ulp miss; and read: the re-creation
path through Chain::setState's restore.
