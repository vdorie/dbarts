# k-internal-parameter: k is the chain's parameter on the internal scale, the sd spelling is sugar, and installs convert nothing

Status: PLANNED 2026-10-09; RULED as a whole (dec-B416), its open calls ruled (dec-B417, dec-B418,
dec-A191). Revised from its blind critique (scratch/kcrit2/critique.md) and then to the rulings, which
widen it by dec-B418's removal of every install conversion. One question is held for the maintainer
([Held for the maintainer](#held-for-the-maintainer)). Nothing is built.

agent: opus implementer, one; one opus reviewer who runs the mutants under Tests.
rng: SHIFTING: no fit's posterior changes through R. NEUTRAL, bit for bit, for every fit whose leaf
prior is written with k in every operation but an install (`setState`, `copy()`, a reload, a warm start)
whose stored mapping differs from the recipient's; every install within one mapping stays bitwise. Draws
move by rounding only for fits written with an sd (ran: within 2.3e-14 over 1000 draws on gaussian; the
critique ran every family and leaf model that takes an sd, within 4.4e-13 over 500 draws, the k ratio
exact). Draws move in trajectory after `setModel` with another leaf prior, an install across a change of
prior into or out of the sd spelling, a re-anchor under a drawn sd, an `xbart` sd grid's warm sweep into
a drawn cell, and every install across a change of mapping, where leaves, k and sigma now go in as stored
on the internal scale where today the leaves are converted and sigma held in response units (dec-B418).
POSTERIOR-CHANGING on one path: a sampler with a named sd re-anchored through the flat C entries, which
stan4bart does every warm-up step when `bart_args` names an sd; today sampling runs under the sd
stretched by the last warm-up re-anchor, after this plan under the named sd, a fix toward dec-A105. No
exact gate reaches it; a new arm of the C API test and stan4bart's suite gate it (Gates). And stan4bart's
restore of a continuous fit across per-chain mappings changes its predictions unless stan4bart changes
in lockstep ([Held for the maintainer](#held-for-the-maintainer)).
window: before the merge to main; engine slices stay serial; after default-rule-per-gap and before
setstate-force-update (dec-A191).
budget: planned ~1670 lines (engine ~330, of which ~170 removed; bridge ~80; R ~110; tests/cpp ~330;
tinytest ~560; help ~110; docs and a benchmark comment ~150); stan4bart's restore (~120) is outside it.
Forecast and stops: [Budget and stops](#budget-and-stops).

## Summary for the maintainer

You approved the plan as a whole (dec-B416) and ruled its three open calls. This is the plan revised to
those rulings, with one question that the install ruling raises.

What the rule is. Every chain holds one number for the leaf spread, k, measured against the data's own
scale: the response mapped to [-0.5, 0.5], or the fixed 3 of probit and about 5.44 of logistic. That
scale is never redefined by how a prior is written. A prior written with an sd is sugar: an sd of s is k
equal to the data's scale over s, and a scaled inverse chi prior on the sd is the matching chi prior on
k. When the response's mapping is re-derived, a named sd is retranslated so it keeps its meaning in
response units, while the chain's current k, and so its current sd, stays put and lags with the leaves
until the next draw (dec-B417). Every install (a warm start, a restored state, a copy, a reload) puts
the chain in as it was stored on the internal scale, read against the receiving sampler's mapping, and
converts nothing (dec-B418).

What changes for a user:

- A prior written with k: draws are bit for bit today's, except an install whose stored mapping differs
  from the sampler's. There the leaves, k and the residual sd now go in as stored instead of being
  converted, so the restored fit is the stored one read on the new scale; the next draws move it where
  the data say. A pure rescaling of the response is then no change at all on the internal scale.
- A prior written with an sd: the same model and the same posterior (fits within 4e-13 of today's on
  every family and leaf model run). The k a fit reports is the real k against the data's scale, and the
  reported reference scale is the data's. A held sd is reported exactly as named.
- setModel, and xbart's sweep over an sd grid, keep the spread across a change of prior, as setLeafPrior
  does. A state restored after a change of prior keeps its spread.
- A drawn sd lags at a re-anchor, as a drawn k does.
- Restoring within one mapping continues the chain bit for bit, as today. Restoring across a change of
  mapping is legal and has no special meaning; the help says so. A Gaussian-process leaf or a model with
  amplitude forests no longer refuses such a state.
- stan4bart with an sd named in its BART arguments now samples under the named sd: a posterior change,
  and a fix.

Rulings reversed, in part: dec-B201's definition of the reference scale under the sd spelling (with the
agents' dec-A124) and the register's matching clause in dec-B393, your words in both standing; dec-A105
for a drawn sd's current value at a re-anchor (dec-B417); dec-B200's conversion on install, dec-B384 and
dec-B407 (dec-B418).

One question for you, raised by the install ruling. Its premise was that a copy or reload shares the
response's mapping. stan4bart's restore of a saved fit does not: each of its chains re-derives the
mapping during its own warm-up, so the chains end on different mappings, and the restore installs all of
them into one sampler rebuilt from the data. Today each chain is converted into that sampler's mapping
and predictions from the kept draws match to 1e-15. Without conversion, in a run of that pattern, the
restored fit was off by up to 0.97 where the fit itself has a spread of 0.70: a restored continuous
stan4bart fit would predict wrongly. Recommended: stan4bart restores one sampler per chain, each created
at that chain's own final mapping, so every install stays within one mapping and is exact, and combines
their predictions; about 120 lines in stan4bart, landing with this slice. The alternatives are a sampler
that holds a mapping per chain, stan4bart no longer re-deriving per chain, or keeping the conversion for
this one case.

The work is about 1670 lines planned, most of it removal and test changes, 2500 to 3340 at the usual
overrun; stan4bart's part is separate. No recorded baseline or snapshot is expected to move.

## Goal

On the internal scale the leaf prior's reference scale, `k.scale`, is a constant of the family and the
response mapping, never redefined by a prior's spelling, and k is each chain's parameter against it.
`normal(sd = s)` is `k = k.scale / s`, and `normal(sd = invchi(df, s))` is `k = chi(df, k.scale / s)`,
retranslated whenever the mapping is re-derived so the named sd keeps its value and distribution in
response units. A re-anchor leaves the leaves and a drawn k as they were, a held k being retranslated; a
change of prior touches neither. Every install puts the chain's parameters (tree structure, leaves, k,
sigma and the rest) in as stored on the internal scale against the recipient's mapping and converts
nothing; a restore within one mapping is bitwise, one across mappings legal and of no special meaning.

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
  retired: [`resolvedNodeScale`](../../src/bartcore/chain.hpp) sets the leaf scale from the family's internal
  value ([`defaultLeafScale`](../../R/model.R): 0.5 gaussian, 3 probit, pi sqrt(3) logistic, 3 nbinom)
  unless a named prior scale is finite, in which case it is that scale over the transform's multiplier.
  [`priorScaleFactor`](../../src/bartcore/chain.hpp) converts the internal scale to the response-unit
  `k.scale` that [`forestCalibration`](../../src/bartcore/chain.hpp) reports.
- The R model carries a named sd as `prior.scale` = 2 s with k fixed at 2, or with `chi(df, 2)` for
  `invchi(df, s)` ([`resolveLeafPrior`](../../R/model.R), the [`dbartsModel`](../../R/A_class.R) slot).
  `invchi(df, 0)` is `chi(df, Inf)` with no named scale.
- After every re-anchor, re-creation, copy, install and warm start, R writes the named scale back
  (retired: [`reissueNamedLeafSd`](../../R/dbarts.R)) through
  [`bartcore_setLeafPrior`](../../src/R_interface_bartcore.cpp) to
  retired: [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp), which rewrites the internal leaf scale.
  The flat C entries [`dbarts_sampler_setOffset`](../../src/C_interface.cpp) and
  [`dbarts_sampler_setResponse`](../../src/C_interface.cpp) re-anchor in the engine and write nothing
  back (read), so on that path a named sd is not kept in response units today.
- A re-creation (a reload, a copy, `setState` on a dead pointer) creates the sampler with the install
  to follow, so [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp) leaves the chains and
  [`Chain::setState`](../../src/bartcore/chain.hpp)'s own scale restore
  (retired: [`installsScale`](../../src/bartcore/chain.hpp)) moves them (the critique, read).
- `$setLeafPrior` keeps the drawn spread by multiplying k by the ratio of the two `k.scale` values
  (retired: [`keepDrawnSpread`](../../R/dbarts.R), retired: [`leafKScale`](../../R/dbarts.R),
  retired: [`Chain::scaleDrawnK`](../../src/bartcore/chain.hpp)). [`Chain::setModel`](../../src/bartcore/chain.hpp)
  re-derives the leaf scale from the model and keeps a drawn k.
- Every install converts units today (read). [`Sampler::setState`](../../src/bartcore/sampler.hpp) and
  the warm start's [`Sampler::installForests`](../../src/bartcore/sampler.hpp) pass a state whose
  mapping differs (retired: [`Sampler::unitsDiffer`](../../src/bartcore/sampler.hpp)) through
  retired: [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp): live and saved leaves times the ratio r of
  the two ranges plus the shift split over the trees, slopes and gp fits times r, variance factors by
  r^(2/m'); a gp leaf or forests with amplitudes under another shift are refused
  (retired: [`unitsRefusalMessage`](../../src/R_interface_bartcore.cpp)). k installs as stored where the
  recipient draws it. sigma is stored in response units ([`Chain::getState`](../../src/bartcore/chain.hpp)
  writes `sigma()`, the internal value times the mapping's scale) and installed in response units
  ([`installDrawnScalars`](../../src/bartcore/chain.hpp)), so it too stays where it was in response
  units. The warm start (dec-B343) carries the same: trees and leaves (converted), k and sigma where the
  new fit draws them (k as stored, sigma in response units), the DART split weights, and the variance
  trees (converted). [`Chain::installForest`](../../src/bartcore/chain.hpp) then restores the scale the
  conversion already put it at. Latents and ordinal thresholds are installed as stored either way.
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

dec-B418 emulated on today's build (scratch/kint/04-no-convert.R, ran): a chain run 100 sweeps,
re-anchored through `setOffset(0.8 x3 + 0.5, updateScale = TRUE)` as stan4bart's warm-up does, run 100
more, its state installed into a sampler created from the same data. Converted, as built: the
recipient's fit equals the donor's to 6.7e-16. Unconverted (the state's `fit.scale` relabelled as the
recipient's, so nothing converts): the fit is off by up to 0.97, against a spread of the donor's fit of
0.70 sd, the mapping's centre having moved as well as its width; k and sigma are as stored.

## The rulings under this rule

"Unchanged": the rule and what it builds stand. "Restated": the outcome stands and the mechanism or
wording changes. "Reversed": part of the ruling no longer holds; the part is named. Words in quotation
marks are the maintainer's; "Register:" marks the decision register's own wording, which is not.

| ruling | what it says | under this rule |
|---|---|---|
| dec-B416 | The plan approved as a whole. The maintainer: "Shouldn't `k` be `k`, the data scale be fixed, and `sd` syntactic sugar?", then "Ah. OK, proceed." | Recorded: this plan. |
| dec-B417 | Open call 2: "It lags, like k. It too is just sugar over the actual parameter." | Recorded: a drawn sd's current value at a re-anchor is k's, which stays; the named prior is retranslated. |
| dec-B418 | Open call 3, widened: every install (warm start, `setState`, `copy()`, reload) puts the chain's parameters in as stored on the internal scale against the recipient's mapping, never converted. The maintainer: "The goal of a warm start is to take tree structure (for the most part), and I guess `k` can ride along for free."; "Does restoring across a re-anchor even make sense? The mental model of the sampler is a linear one through time, with each sample being a snapshot."; "OK, let's do it." | Recorded: step 7 and step 8. Its premise that a copy or reload shares the mapping fails for stan4bart's restore (Held for the maintainer). |
| dec-A191 | (orchestrator's call; "Up to you. I don't care how this lands.") The engine order: gp-copy, the width-weighted default cut rule, this plan carrying dec-B418, setstate-force-update, response-scale-rows, then dec-B404's move and dec-B415's survey | Recorded: Composition. |
| dec-B201 | Register: `k.scale` is the value k is relative to; under an sd-named prior twice the sd or `invchi()` scale named, set so that the engine's k sits at 2. The maintainer: "I want k to be the engine's k always."; "Ultimately, the exact value isn't itself important. The only thing that really matter is that users have a way to convert between k and standard deviation." | Reversed in the register's definition: `k.scale` is the data's scale under every spelling. The maintainer's words stand: k is the engine's k, and spread = `k.scale` / k converts. |
| dec-A124 | (agents' call) Register: `k.scale` the data's under k and twice the named sd under sd, so the data's scale is not readable under an sd-named prior | Reversed with dec-B201: the data's scale is readable on every fit. |
| dec-B141 | `getK` reports the engine's k whatever the spelling | Unchanged; under an sd spelling the number is `k.scale` / spread. |
| dec-B192 | Register: a fit carries the k the sampler recorded, on an sd-named fit relative to twice the named scale. The maintainer: "The samples of `k` (if there are any) should be stored as `k`. End-users should use `extract`" | Restated: relative to the data's scale. The maintainer's words stand. |
| dec-B193, dec-B194 | extract answers k and `leaf.prior.sd` on every fit; one number for a held value | Unchanged; `leaf.prior.sd` values unchanged, k values move under an sd spelling. |
| dec-B376 | Register: no fit stores the leaf sd; `fit$fixed$k` for a held k | Unchanged; a held sd's `fit$fixed$k` becomes `k.scale` / s, and extract and summary read the held sd itself from the recorded prior, exactly. |
| dec-A105 | Register: named by k (relative) or sd (absolute); a named sd stays in response units across every call that re-anchors the response scale, re-applied after it, as a fixed residual sd and the variance forest's prior already do. The maintainer, on re-anchoring: "OK, proceed using it." | Reversed by dec-B417 for a drawn sd's value in force at a re-anchor: it lags with the leaves until the next draw ("It lags, like k. It too is just sugar over the actual parameter."). Restated for the named sd, fixed or the scale of its `invchi()`: retranslated at every re-anchor, keeping its value and distribution in response units, now on the flat C path too. |
| dec-B356 | A switch between spellings of a drawn prior keeps the spread. Register: k becoming k_old x k.scale_new / k.scale_old. The maintainer: "B. Keep the spread." | Restated: `k.scale` never differs, so k is untouched and the spread is kept. |
| dec-B369 | a changed `invchi()` scale keeps the spread at the call ("OK, then A. Keep the spread at the call.") | Restated as dec-B356. |
| dec-B392 | a fixed sd turned drawn keeps the sd in force ("They're trying to actively set the sd of the leaf prior, so both of their values should be interpretted as such.") | Restated: the fixed sd s is k = `k.scale` / s, and the drawn prior keeps that k. |
| dec-B393 | Every switch keeps the sd into a drawn prior; a stated value is literal. Register: `k.scale` always belongs to the new prior's spelling (the data's scale under k, twice the named value under sd); changing the prior never moves the state; no setter for k exists. The maintainer: "if a person installs a fixed value of k, they clearly mean it to be interpretted literally so that should also reset the k.scale", "Treat drawn k priors and fixed k priors the same - reset the k scales. The person is trying to put a distribution on k itself and otherwise the interpretation is off.", "I don't see changing the prior as changing the state / the active parameter. For that, you'd need a setK() function.", "Yes, keep until the next draw." | Rule stands. Reversed: only the register's "twice the named value under sd". The maintainer's words hold: every prior's `k.scale` is the data's, so a k, fixed or drawn, is read literally against it, and neither k nor the spread moves at a change of prior. |
| dec-B401 | "The k itself should literally transfer - a parameter is a parameter, regardless of the prior. It may be a bad fit, but so be it." | Unchanged; the spread is now kept too, since no prior moves `k.scale`. The setState plan's H4 closes. |
| dec-B384, dec-B407 | A state installed onto another response spread keeps its spread. Register: k becoming k_state x k.scale_recipient / k.scale_state. The maintainer: "3. Use your recommendation." (dec-B384); "Keep it for scale changes. If it helps at all, I guess we can think of `k` as just the internal representation and the parameter as the sd." (dec-B407) | Reversed by dec-B418: k installs as stored on the internal scale; against a wider mapping the spread widens with it, as the leaves do. Never built. |
| dec-A146 | (agents' call, one of those "not put" to the maintainer) Register: a warm start follows the setState rule; a drawn k installs as stored | Superseded by dec-B418, one rule for every install. |
| dec-B343 | "A. The donor gives its trees, and sigma and k only where the new fit draws them (as built)." | Unchanged in what transfers; restated by dec-B418 in its units: leaves, k and sigma as stored on the internal scale, never converted (today the leaves are converted and sigma held in response units). |
| dec-B200 | A state stored in other units is converted on install; the mapping is model, recorded on the R object and moved only at creation or a re-anchor. The maintainer: "Convert on install." | Reversed by dec-B418 in its conversion: every install puts the state in as stored on the internal scale against the recipient's mapping ("OK, let's do it."), and the refusal of a gp leaf or amplitude forests under another shift goes. Stands: the mapping is model, recorded, and moved by nothing but creation and a re-anchor; an install no longer moves it at all. |
| dec-B195, dec-B196 | A state carries no model. The maintainer: "Model when fixed, state when drawn." | Unchanged. A held named sd's k is model, retranslated with the data's scale. |
| dec-B254 | Register: `setModel` changes parameters, the leaf scale as k or sd among them, never structure. The maintainer: "Parameters yes, structure no." | Unchanged in scope; it now keeps the spread across a change of leaf prior (the note's defect a), xbart's sd grid included. |
| dec-B396 | Register: `copy()` stores the current state when none is stored | Unchanged. |
| dec-A107 | Register: `setLeafPrior` changes only the spread or its prior | Unchanged. |
| dec-A121 | (agents' call) the engine's k carried across a switch of spelling | Already revised by dec-B356; moot. |
| dec-A187 | (agents' call) `setLeafPrior` keeps the spread through `scaleDrawnK`, which state-install-keeps-spread would reuse | Superseded: no write moves `k.scale`, and the install's conversion divides k itself; the primitive goes. |
| dec-B371, dec-A190 | the probit rescaling step divides k by the drawn factor | Unchanged: it acts on the internal k and the leaves together. |
| dec-A13 | no cap on a drawn k; `chi(df, Inf)` improper | Unchanged; `invchi(df, 0)` is still `chi(df, Inf)`. |
| dec-B331, dec-B330 | a re-anchor converts saved draws, the live leaves kept; a gp sampler holding saved draws refuses a re-anchor ("Eh, you can keep converting. Let's not be too paternalistic.") | Unchanged: they govern a re-anchor, not an install. An install no longer converts saved draws (dec-B418), so a re-anchor and an install now differ there; the help says so. |
| dec-B361 | prior draws under the prior the sampler runs under | Unchanged; the private sampler's named sd is translated by the engine. |
| dec-B362, dec-B364, dec-B405 | several forests refuse a re-derivation; a kept scale (fewer than two values); the status through the C API | Unchanged. A kept scale moves no transform, so nothing is retranslated. |
| dec-B122 | a variance forest recalibrates on a re-anchor | Unchanged (it has no k). |
| dec-B142, dec-B253, dec-B275 | the multi-forest writer and a forest's sd | Unchanged: a named leaf-prior sd is refused on several forests, multinomial and hurdle fits. |

Reversed, in sum: the register's definition of `k.scale` under the sd spelling in dec-B201 (with
dec-A124's); the register's "twice the named value under sd" in dec-B393; dec-A105 for a drawn sd's value
in force at a re-anchor (dec-B417); dec-B200's conversion on install, dec-B384 and dec-B407 (dec-B418).
Superseded: dec-A146's warm-start call and dec-A187's primitive. Every other ruling stands, several with
their mechanism restated; every quoted word of the maintainer's stands.

## Operation by operation

Today is the run above. "New" uses the same fixture: data's `k.scale` 1.9694, rho = `k.scale` / (2 s)
the factor between today's and the new k under an sd spelling. New values are arithmetic on the run
values, not runs, but for the unconverted install, emulated (04-no-convert.R). "Draws" says which fits
move.

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
| install within one mapping: a sampler's own restore, `copy()`, a reload of a state stored since the last re-anchor | as stored, bitwise | unchanged, bitwise | none |
| install across a change of prior, same mapping | k as stored; spread 0.9132 where one side is sd-spelled | k as stored; spread 0.8992 | where one side is sd-spelled |
| install into a sampler holding k or sigma fixed | the sampler's | unchanged | none |
| restore after a re-anchor (store, re-anchor, restore) | leaves and saved draws converted into the new mapping; k as stored; sigma held in response units | leaves, saved draws, k and sigma as stored on the internal scale, read against the new mapping, as the live chain held them at the re-anchor | every such restore |
| `setState`, `copy()`, reload onto a response 3 times as wide, k spelling | leaves converted (the function kept); k 2.1901, spread 2.6977; sigma 0.3933 | nothing converted: the function, the spread and sigma all 3 times wider (k 2.1901, spread 2.6977, sigma 1.1799), a pure rescale being no change on the internal scale | every such install |
| the same, sd spelling (donor at k 2.4572 in the new numbers, spread 0.8015) | k as stored, spread kept 0.8015 | k 2.4572 as stored, spread 2.4044 | every such install |
| install from a mapping with another centre (another offset; stan4bart's restore) | converted: the donor's fit to 6.7e-16 (ran) | as stored: fit off by up to 0.97 against a fit sd of 0.70 (ran, emulated) | every such install; stan4bart: Held |
| install of a gp leaf, or amplitude forests, stored under another shift | refused | installed as stored | was an error |
| warm start onto another mapping (`bart(warm.start = )`, `installTrees`) | leaves converted; k as stored (spread 0.8992 to 2.6977 on 3 y); sigma in response units | leaves, k and sigma as stored on the internal scale; the next draws move them | every warm start across mappings |
| warm start, same mapping, another prior | k as stored; 0.8992 to 0.9132 | k as stored; spread kept | where one side is sd-spelled |
| re-anchor, k spelling | k kept; spread stretches with the fit | unchanged | none |
| re-anchor, fixed sd `0.5` on 3 y | k 2, `k.scale` 1, spread 0.5 | k retranslated 11.8164, `k.scale` 5.9082, spread 0.5 | rounding |
| re-anchor, drawn `sd = invchi(3, 1)` on 3 y | k kept, spread 0.8015 held; next draw from invchi(3, 1) | k kept (2.4572), spread lags to 2.4044 with the fit; prior retranslated to chi(3, 5.9082), the next draw from invchi(3, 1) in response units (dec-B417) | the sweep after each such re-anchor |
| re-anchor through the flat C entries, named sd | not retranslated (read): the sd stretches with the response, and sampling after warm-up runs under the last stretch | retranslated, as through R | POSTERIOR: stan4bart with a named sd in `bart_args` now samples under the named sd (a fix toward dec-A105) |
| `setResponse`/`setOffset` at `updateScale = FALSE`; a kept scale (dec-B364) | nothing moves | unchanged | none |
| probit rescaling step | k divided by the factor | unchanged | none |

The setLeafPrior chain, step by step under the new rule (same fixture): before k 2.1901, spread 0.8992;
into `sd = invchi(3, 1)` k 2.1901, 0.8992; into `invchi(3, 2)` the same; into `sd = 0.5` k 3.9388, 0.5;
into `k = chi(1.5, 2)` k 3.9388, 0.5 (as today); into `k = 3` k 3, 0.6565 (as today).

A re-anchor and an install now differ on the kept draws: a re-anchor still converts them so `predict`
on them does not move (dec-B331), while an install brings them in as stored, so after a restore across
mappings `predict` on the kept draws reads them on the new scale.

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
- A state's `sigma` entry holds the internal value (step 8); no R code reads it (read: grep).
- The C API ([dbarts.h](../../inst/include/dbarts/dbarts.h)): [`dbarts_results`](../../inst/include/dbarts/dbarts.h)
  and [`dbarts_draw`](../../inst/include/dbarts/dbarts.h) carry k per draw; for an sd-named sampler the
  numbers move by rho. No entry, struct or enum changes, so
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move (it covers signatures, enums
  and layouts, read): not an ABI event. The header states nothing about what k is relative to (read);
  one sentence is added on the `k` field. No flat entry installs a state or converts units (read: the
  entry list), so dec-B418 changes nothing there.
- stan4bart (bartcore 963956b, read): copies the k draws into its fit's `k` and reads
  `dbarts_sampler_kIsSampled`; it forwards `bart_args$leaf.prior` to dbarts, so a user may name an sd
  there, and then calls `dbarts_sampler_setOffset` with `updateScale` true during warm-up. Its k numbers
  move for such a fit and its named sd is now kept across the warm-up re-anchors. Its tests write k
  (`normal(k = 3)`, test-09-bartArgs.R, read). Its restore of a saved fit goes through R's `setState`
  and is changed by dec-B418 (Held for the maintainer). bartCause reads `fit$k` from a k-spelled fit;
  treatSens writes `chi()` and fixed k; bairrtt names no leaf prior; none installs a state (read, their
  dbarts-1.0 and main branches).

## Where it is built

The engine changes; the translation cannot live in R alone. Creation needs the data's scale, which the
engine computes from the rows, and the flat C re-anchor never returns to R (Open call 1). dec-B418
removes more engine code than the translation adds.

Engine ([chain.hpp](../../src/bartcore/chain.hpp), [sampler.hpp](../../src/bartcore/sampler.hpp),
[facade.hpp](../../src/bartcore/facade.hpp)):
1. A forest holds `namedSd`, NaN when the prior is written with k. The `priorScale` fields of
   [`ModelParameters`](../../src/bartcore/chain.hpp) and the creation options become `namedSd`.
2. The internal leaf scale is always the family's: retired: [`resolvedNodeScale`](../../src/bartcore/chain.hpp)
   loses its named branch, at creation and in [`Chain::setModel`](../../src/bartcore/chain.hpp).
3. One translation routine: where `namedSd` is finite, a fixed k is `k.scale / namedSd` and a drawn
   prior's chi scale is `k.scale / namedSd`, with `k.scale` the reader's own expression
   ([`priorScaleFactor`](../../src/bartcore/chain.hpp) times the leaf scale), so the reader and the
   write agree bit for bit. At creation a drawn k starts at the same value. A drawn k is never written by
   it.
4. It runs wherever a chain's response mapping changes, beside the variance forest's recalibration that
   already sits at each such site: the `updateScale` arms of
   [`Chain::setResponse`](../../src/bartcore/chain.hpp) and [`Chain::setOffset`](../../src/bartcore/chain.hpp);
   [`applyNewData`](../../src/bartcore/chain.hpp) (setData); and
   [`moveScale`](../../src/bartcore/chain.hpp), under [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp).
   After step 7 no install moves the mapping, so these four are all; the reviewer greps chain.hpp for
   every call that changes the response's mapping (setResponse and setOffset at `updateScale`, setData,
   restoreScale) and finds the routine after each. tests/cpp holds the invariant after each site.
5. retired: [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp) becomes a named-sd writer: it sets
   `namedSd` and translates, leaving a drawn k; a write of the value in force is skipped. Its facade
   virtual on [`SamplerBase`](../../src/bartcore/facade.hpp) is renamed with it.
6. [`Chain::setModel`](../../src/bartcore/chain.hpp): a fixed k is the model's or the translation; a
   drawn prior is translated; a drawn k is kept, as today.
7. No install converts (dec-B418). Removed: the units pass in
   [`Sampler::setState`](../../src/bartcore/sampler.hpp) (the `unitsRefused` out-parameter, and the
   `valuesMoved` contribution to its `altered` report) and in
   [`Sampler::installForests`](../../src/bartcore/sampler.hpp) (`WarmStartResult::unitsMismatch`),
   retired: [`Sampler::unitsDiffer`](../../src/bartcore/sampler.hpp),
   retired: [`Chain::convertStateUnits`](../../src/bartcore/chain.hpp), and the scale restore at the head of
   [`Chain::setState`](../../src/bartcore/chain.hpp) and
   [`Chain::installForest`](../../src/bartcore/chain.hpp) (retired: [`installsScale`](../../src/bartcore/chain.hpp)
   and the variance leaf's recalibration after it). The leaf, slope, gp and variance-factor helpers stay
   for [`restateSavedDraws`](../../src/bartcore/chain.hpp), which a re-anchor still runs (dec-B331).
   [`Sampler::setAnchor`](../../src/bartcore/sampler.hpp) moves the chains on every call: a re-creation
   no longer leaves the move to the install it is followed by, so a copy or reload is at the recorded
   mapping before its state goes in. A state's `fit.scale` is still written, recording the mapping it
   was drawn under, and no install reads it.
8. sigma on the internal scale: [`Chain::getState`](../../src/bartcore/chain.hpp) stores the internal
   value where the chain draws sigma, and [`installDrawnScalars`](../../src/bartcore/chain.hpp) writes it
   back unscaled, for setState and the warm start alike. A held sigma is model and is untouched (dec-B196).
   Latents and ordinal thresholds install as stored, as today.
9. retired: [`Chain::scaleDrawnK`](../../src/bartcore/chain.hpp) and its facade virtual go: no caller remains.
   `--preclean` on every install (facade virtuals change).

Bridge ([R_interface_bartcore.cpp](../../src/R_interface_bartcore.cpp),
[R_interface.cpp](../../src/R_interface.cpp)): the model parse computes the named sd once from the R
model, `prior.scale` over the fixed k or the chi scale (2 s / 2, exact), for creation
([`optionsFromParsed`](../../src/R_interface_bartcore.cpp)) and `setModel` alike;
[`bartcore_setLeafPrior`](../../src/R_interface_bartcore.cpp) takes the named sd;
retired: [`bartcore_scaleDrawnK`](../../src/R_interface_bartcore.cpp) and its registration go;
retired: [`unitsRefusalMessage`](../../src/R_interface_bartcore.cpp) and its two uses (setState's and the warm
start's) go; [`bartcore_anchor`](../../src/R_interface_bartcore.cpp) loses its install-follows argument
(arity 3 to 2). The flat C entries do not change: the engine retranslates under them.

R ([dbarts.R](../../R/dbarts.R)):
- retired: [`reissueNamedLeafSd`](../../R/dbarts.R) and its six callers go (re-anchors, setData, the install
  shared by re-creation and setState, copy, installTrees, samplePriorPredictive): the engine keeps the
  invariant.
- retired: [`keepDrawnSpread`](../../R/dbarts.R) and retired: [`leafKScale`](../../R/dbarts.R) go; `$setLeafPrior`'s
  same-prior write passes the named sd ([`writeLeafPrior`](../../R/dbarts.R)).
- [`reportLeafPrior`](../../R/dbarts.R) takes an sd-named specification from the model's own
  `leaf.prior`, not from `k.scale` / k, so a fixed sd reads back exactly.
- [`applyAnchor`](../../R/dbarts.R) and [`recreatePointer`](../../R/dbarts.R) lose `installFollows`.
- The docstrings of `$setState`, `$copy` and `$installTrees` (rc-codoc) say what dec-B418 requires
  (step 7 of Steps).
- [`extractParameter`](../../R/generics.R): a held leaf-prior sd on a single-forest fit whose recorded
  prior names a numeric sd is that sd; summary's held line follows through it.
- [`resolveLeafPrior`](../../R/model.R) and the model's encoding are unchanged (settled below); their
  comments and the [`dbartsModel`](../../R/A_class.R) slot comment say what the encoding now means.
- xbart needs no code change: its warm sweep calls [`bartcoreSetModel`](../../R/xbart.R) per cell, which
  takes the engine's new `setModel` rule.

## State format and saved objects

- A stored k means the chain's k against the data's scale; a stored sigma, after step 8, the internal
  value. For a k-spelled sampler k means what it meant; sigma and an sd-spelled sampler's k change
  meaning. No format version moves: the registry rule at
  [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) holds both numbers at 1 until the first
  release, no serialized format having shipped (read), and 0.9-34 states are refused already.
- Consequence: a sampler that draws sigma, or draws an sd-spelled k, saved on a development build and
  reloaded on the new build, restarts with sigma read as internal (its response-unit value divided by
  the scale, wrong by that factor) and such a k times rho. Held values are model and are unaffected. No
  migration (Open call 5).
- A saved fit carries k and `k.scale` as a pair; the extractor divides one by the other, so old and new
  fits both report the right `leaf.prior.sd` (read). The R model's encoding does not change, so a saved
  model object reads as before.
- No `.rds` fixture in the tree holds a state (read: the tracked `.rds` files are the equivalence
  baselines and the 0.9-34 classic compare).

## Composition with the other plans

The setState slice (setstate-force-update.md on wt/install-surface-plan, 20e08954, read): it is written
against this plan's dec-B418 form. It expects a `Sampler::setState` with no units pass and no
`unitsRefused`, which step 7 delivers; its B5 (`kFollowsUnits`) and its H4 are gone; it owns the forced
and unforced forms, the kept store's size and the missing-value record, none of which this plan touches.
Its `altered` report loses the units pass's contribution here, so a cross-mapping install is clean.

response-scale-rows (wt/rsr-recheck, faf8ca7d, read): its step 5's re-derivation after a count fit's
creation under a mask "restates a named sd" from R; under this plan the engine does, and the R write
goes. Its tinytest's named-sd arm and its `k.scale` literal then read the data's scale over the rows in,
under both spellings. Its "Held" arm (readers unchanged across a copy and a reload) holds: a copy and a
reload share the mapping. Its flat-entry status (dec-B405) and dec-B364's keep need nothing here: a kept
scale moves no mapping.

Order (dec-A191): gp-copy-continuation (landed, dec-A192), default-rule-per-gap, this plan,
setstate-force-update Parts A and B, response-scale-rows B, then dec-B404's logistic move and dec-B415's
samplers if their studies call for them. The wording sweep queued last stays last. stan4bart's restore
change, if the held question is answered as recommended, lands in lockstep with this slice on its
bartcore branch.

## Constraints

- Frozen: every k-spelled fit's bits in every operation but an install across a change of mapping;
  every install within one mapping, bitwise; the R model's encoding; the C API's signatures, structs and
  enums.
- Out of scope: several forests, multinomial and hurdle fits (a named sd is refused there, unchanged);
  the variance forest's leaf prior; the probit rescaling step; a setter for k (dec-B393: "For that, you'd
  need a setK() function.", not asked for); what a re-anchor does to the live leaves and to kept draws
  (dec-B331, unchanged); a warm start's conversion of a linear leaf's slopes between two covariate
  standardizations ([leaf-conversions.md](../design/leaf-conversions.md)), which is about the predictors'
  frame, not the response's mapping, and which dec-B418 does not name (Open call 6).

## Steps

1. Engine steps 1 to 3 and 6 (creation and setModel), tests/cpp for them; the emulation identity:
   `sd = s` against `k = k.scale / s` within 1e-12 over 1000 draws, and bitwise where rho is 1.
2. Engine step 4 (the mapping-changing sites) and the invariant check.
3. Engine step 5 and the bridge's `setLeafPrior`; engine step 9 and the bridge removal.
4. Engine steps 7 and 8 (no install converts; sigma internal), the bridge's refusal message and
   `bartcore_anchor`'s argument, their tests/cpp checks.
5. R: the removals, `writeLeafPrior`, `reportLeafPrior`, `applyAnchor`, `extractParameter`'s held sd,
   docstrings, comments.
6. tinytest: the new file and the changed pins (Tests).
7. Help. [dbartsSampler-class.Rd](../../man/dbartsSampler-class.Rd): `setState`, `copy` and the Saving
   section say, for dec-B418, that a state goes in as stored - its trees, leaf values, k, sigma and the
   kept draws are numbers on the sampler's internal scale and are read against this sampler's response
   mapping, never converted; a restore within one mapping continues the chain exactly; a restore after
   a re-anchor, or a state from another sampler, is legal and has no special meaning, the next draws
   moving the parameters where the data say; a kept draw restored across mappings predicts on the new
   scale, where a re-anchor converts the kept draws. `installTrees` and [bart.Rd](../../man/bart.Rd)'s
   `warm.start` say the same of the donor: its trees, and k and sigma where this fit draws them, as stored
   on the internal scale, not converted to this fit's response. Also `getLeafPrior`'s `k.scale`; the
   `setLeafPrior` paragraph (k untouched because `k.scale` never moves); the named-sd paragraph (a drawn
   value lagging). [xbart.Rd](../../man/xbart.Rd): its `sd` item, "invchi(df, c) starts its chain at the
   spread c" holding for the first cell, each later cell of the warm sweep starting at the previous
   cell's spread, whichever spelling. [dbartsPriors.Rd](../../man/dbartsPriors.Rd): the named sd
   translated to k against the table's `k.scale`; its "absolute" paragraph.
   [bartBT.Rd](../../man/bartBT.Rd) where it says what a fit's k is relative to; the `k` field in
   dbarts.h. NEWS: nothing for the sd spelling (new in 1.0-0); the install rule is new in 1.0-0 too, 0.9-34
   having installed a state as stored with no record of its units, so NEWS states none.
8. Docs and a benchmark comment: [state-not-model.md](../design/state-not-model.md)'s first table
   (`fit.scale`: recorded, not compared) and its section
   [The response transform](../design/state-not-model.md#the-response-transform) (installs as stored; the
   gp and amplitude refusal gone; a re-anchor's rollback is a restore of the response with the state);
   [leaf-conversions.md](../design/leaf-conversions.md) where it says an install's conversion;
   [leaf-scale-rules.md](../design/leaf-scale-rules.md) Status (decided, with this plan) and a short
   closing section; [composition-matrix.R](../../benchmarks/R/composition-matrix.R)'s "calibration"
   probe comment, which describes the reference k of 2 (not a gate).
9. Gates (below), then the records at landing: the ledger entry for the calls made, this plan's Status
   and Landing, TODO (state-install-keeps-spread closes, reversed).

## Tests

tests/cpp ([test_sampler.cpp](../../tests/cpp/test_sampler.cpp), [test_state.cpp](../../tests/cpp/test_state.cpp),
[test_facade.cpp](../../tests/cpp/test_facade.cpp)):
- Creation: a fixed named sd's k equals the reader's `k.scale` over it bitwise, the internal leaf scale
  equals a k-spelled sampler's; a drawn one's chi scale and start likewise.
- Invariant, one check per site of step 4: after each re-anchor setter at `updateScale`, applyNewData and
  setAnchor, a fixed k equals `k.scale` over the named sd, a drawn chi scale likewise, and a drawn k is
  the value before the move.
- setModel keeps a drawn k across a change of named sd and of spelling.
- The named-sd writer: fixed and drawn, a drawn k untouched, an equal write moving no bit.
- Installs as stored: a state from another mapping, into setState and into installForests, leaves the
  recipient's mapping unmoved, and the chain's state read back (getState) equals the installed one in
  every tree, leaf, slope, gp fit, variance factor, saved draw, k and sigma, bitwise; a gp leaf and an
  amplitude-coupled state under another shift install. The 19 conversion checks in test_state.cpp
  (ran: grep) become these.
- sigma round trip: getState then setState of a chain drawing sigma, within one mapping and across one,
  leaves the internal sigma bitwise.
- setAnchor moves the chains on a re-creation; the chain's mapping after a re-creation and install is the
  record's.
- The facade list: the renamed virtual, scaleDrawnK gone.

tinytest, new `test-k-internal.R`: one block per row of the operation table, asserting k, `k.scale`, the
spread and sigma (to 1e-12) on the probe's fixture; the binary identity (probit `sd = invchi(1.5, 1.5)`
and `k = chi(1.5, 2)` identical); the fit readers listed above; a write-back of
`getLeafPrior()$leaf.prior` moving no bit under both sd forms. And:
- Within one mapping: a sampler's own restore, a copy and a reload (`saveRDS` and `readRDS`) continue
  bit for bit, with a held and a drawn named sd, before and after a re-anchor whose state was stored
  after it.
- Across mappings: setState and installTrees onto 3 y put in the stored internal numbers (the
  recipient's next stored state equals the installed one apart from `fit.scale`, before any sweep), sigma
  3 times its response-unit value, the spread 3 times; no conversion, no refusal for a gp leaf or a bcf
  fit under another shift.
- A re-anchor still converts kept draws (`predict` on them unchanged), and an install across mappings
  does not (`predict` on them reads the new scale): the two pinned side by side.
- Re-creation after a re-anchor, named sd held and drawn: copy, reload and setState on a dead pointer
  put the chains at the recorded mapping; the held spread is s to 1e-12, the drawn chi scale
  `k.scale` / s.
- xbart's path: [`bartcoreSetModel`](../../R/xbart.R) from a fixed `sd = 0.25` into `invchi(3, 1)` keeps
  the spread at 0.25, and from a drawn chain into `invchi(3, 2)` keeps its spread; an `xbart` grid mixing
  a fixed and a drawn cell of different scales runs and repeats under a seed. (The existing grid tests
  pair each drawn cell with a fixed cell of the same s and cannot see the change.)
- A held sd read exactly: for s in 0.7, 1.3 and 0.1 on several responses, `extract(fit,
  "leaf.prior.sd")` and summary's held line are identical to s.

tinytest, [test-capi.R](../../inst/tinytest/test-capi.R), a new arm: a sampler with a held and with a
drawn named sd, re-anchored through the compiled consumer's flat `setOffset` with `updateScale` on an
offset that moves the range; the held spread stays s to 1e-12, the drawn chi scale is `k.scale` / s, a
drawn k is unchanged. It fails today. The file skips where the consumer cannot be compiled, as now.

tinytest, changed (ran: grep): the sd-spelling pins in
[test-calibration-midchain.R](../../inst/tinytest/test-calibration-midchain.R) (40 k readers),
[test-fit-stores-k.R](../../inst/tinytest/test-fit-stores-k.R) (16),
[test-calibration-creation.R](../../inst/tinytest/test-calibration-creation.R) (11),
[test-leaf-prior-k-or-sd.R](../../inst/tinytest/test-leaf-prior-k-or-sd.R) (its header comment; its
translation pins hold), [test-nbinom.R](../../inst/tinytest/test-nbinom.R),
[test-embedding-recipes.R](../../inst/tinytest/test-embedding-recipes.R); and the conversion and
units-refusal pins, [test-state-not-model.R](../../inst/tinytest/test-state-not-model.R) (10 matching
lines), the bcf, multi-forest and gp files that pin "another response shift" refused
([test-bcf.R](../../inst/tinytest/test-bcf.R), [test-multi-forest-seam.R](../../inst/tinytest/test-multi-forest-seam.R),
[test-forest-basis-r5.R](../../inst/tinytest/test-forest-basis-r5.R) among them). The full suite finds the
rest.

Mutants the reviewer runs, each failing a test: the translation skipped at each of step 4's four sites in
turn; a drawn k retranslated at a re-anchor (it must lag); a drawn sd started at k 2; setModel
re-expressing k against a ratio (today's setLeafPrior arithmetic); a units pass left in setState, or in
installForests; k divided by the units ratio on install (dec-B384's build); sigma stored or installed in
response units; setAnchor not moving the chains on a re-creation; the install restoring the state's
`fit.scale`; the gp or amplitude refusal kept; `k.scale` reported as twice the named sd; reportLeafPrior
reading a fixed sd as `k.scale` / k; `extractParameter` computing a held sd as `k.scale` over the held k;
the named-sd writer moving a drawn k.

## Gates

On the slice's tip and its own library, independently of the implementer, shifting class:
- tests/cpp plain and under `-fsanitize=address,undefined`; R-loaded ASAN over the touched test files.
- The full tinytest suite (`at_home = TRUE`), test-capi.R's new arm compiled and run (a skip there is a
  stop: it is the only check of the flat path in this package); stan4bart's suite against the build,
  with its restore change if the held question is answered as recommended. The flat-path posterior
  change has no exact gate: stan4bart's tests write k (read), so its suite checks that nothing else
  moved, and test-capi.R checks the fix.
- Reference build, `--preclean`: the equivalence trio `compare --bitwise` (equivalence.R also
  `--strict-coverage`) against [MANIFEST](../../benchmarks/baselines/MANIFEST)'s current files, and the
  four arm64 snapshot files. Expected: all identical. No equivalence scenario or snapshot file writes an
  sd spelling, a leaf-prior writer, `setModel`, a state install or a warm start; the only re-anchors are
  `setData` scenarios under the k spelling (ran: grep). Any move is a stop.
- Every exact gate in quick mode ([exact-gates.yaml](../../.github/workflows/exact-gates.yaml)): what a
  fit carries changes value. [aft-exact.R](../../benchmarks/R/aft-exact.R),
  [t-exact.R](../../benchmarks/R/t-exact.R) and [logistic-reference.R](../../benchmarks/R/logistic-reference.R)
  write a fixed sd and move by rounding; negbin-mixing.R installs a sampler's own state within one
  mapping, bitwise (read, the setState plan's callers). Each must pass. No baseline moves, so no oracle is
  needed; if a gate fails, the emulation identity (step 1) is the oracle.
- [composition-matrix.R](../../benchmarks/R/composition-matrix.R) writes a fixed sd through
  `setLeafPrior` but is not a CI gate; it is run once, all cells as recorded.
- `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.
- Not hot-path: the translation runs at a mutation, never per sweep, and the install loses work; no
  bench-sampler compare.

## Budget and stops

Planned ~1670 lines: the reviewed plan's ~1020, less dec-B384's division and its tests (~40), plus
dec-B418's removal of every install conversion (engine ~170 removed and ~40 added, bridge ~40, R ~30,
tests/cpp ~120, tinytest ~170, help ~40, design docs ~80). stan4bart's restore change (~120) is its own.
Engine slices on this surface ran 1.1 to 2.9 times their plans; the last ran 1.76. Forecast at 1.5 to 2
times: 2500 to 3340. Stop at 4200 and report, without working around, when:
- the diff passes 4200 lines;
- any equivalence scenario or snapshot moves;
- a k-spelled fit moves a bit in any operation but an install across a change of mapping;
- an install within one mapping moves a bit;
- an sd-spelled fit departs from the emulation by more than 1e-10 within 1000 draws;
- a site that changes the response's mapping is found that step 4's routine cannot be put after;
- test-capi.R's consumer cannot be compiled on the gate machine;
- stan4bart's suite fails on its restore and the held question is not answered;
- a reading needs a state-format change beyond sigma's meaning, or a test needs a call no ruling,
  settled call or open call here makes.

## Settled in planning

- The R model keeps its encoding (`prior.scale` 2 s with k at 2 or `chi(df, 2)`); the bridge divides it
  out. Saved models read as before, the translation pins hold, and nothing R-side moves but comments.
- A drawn sd's start is the translation of today's: the chain starts at the named sd. A k-spelled
  `chi()` still starts at 2.
- `getLeafPrior` reports a named sd from the model, exactly; the engine's fixed k is `k.scale` / s,
  which need not divide back to s's bits.
- sigma is carried on the internal scale by storing the internal value in the state (step 8), not by
  reading the state's `fit.scale` at install: dec-B418 has the install read no mapping but the
  recipient's.

## Held for the maintainer

stan4bart's restore across per-chain mappings. dec-B418 was put with "a copy or reload shares the
mapping". stan4bart's restore of a saved fit (restoreBartSampler, reached after a reload and at fit
end; read on bartcore 963956b) does not: each chain runs its own dbarts sampler, re-deriving the mapping
through `dbarts_sampler_setOffset(..., updateScale)` during its own warm-up, so the chains end on
different mappings (dec-B191 records the same finding), and the restore installs every chain into one
sampler built from the fit's control, model and data. Today each chain is converted into that sampler's
mapping, and kept-draw predictions match to about 1e-15 (dec-B200's measurement; 6.7e-16 in
04-no-convert.R). Unconverted, in that pattern, the restored fit is off by up to 0.97 against a fit sd of
0.70 (ran, emulated): a restored continuous stan4bart fit would predict from the wrong mapping. Binary
stan4bart fits are unaffected, their mapping being fixed.
- Recommended: stan4bart restores one dbarts sampler per chain, each created at its chain's final
  mapping (the model's recorded mapping set from that chain's stored `fit.scale`), so every install is
  within one mapping and bitwise, and combines the chains' predictions itself. It keeps dec-B418 whole
  and its linear model of a chain, each chain continuing its own history. About 120 lines in stan4bart's
  R and its prediction path, landing in lockstep with this slice; its suite's restore tests are the gate.
- A sampler holding one mapping per chain: an engine change against dec-B200's one-mapping sampler, and
  chains of one sampler then run under k-spelled priors of different widths in response units, which
  dec-B191's discussion set aside.
- stan4bart re-deriving no mapping per chain (or one shared across chains): changes stan4bart's
  posterior and the warm-up practice dec-B330 describes.
- Keeping the conversion for an install into a sampler re-created to receive the state: a partial
  reversal of dec-B418.

## Ruled

Open call 2 (a drawn sd at a re-anchor): dec-B417, it lags. Open call 3 (the warm start across a response
spread): dec-B418, widened to every install, nothing converted. Open call 4 (order): dec-A191.

## Open calls

For the orchestrator, each with a recommendation.

1. Where the translation lives. The engine holds the named sd and retranslates wherever the mapping
   moves (recommended): the flat C re-anchor, which stan4bart calls every warm-up step, keeps a named sd
   as R's re-anchor does, and one invariant replaces six R write-backs. The alternative, R and the
   bridge translating after the engine computes the scale, leaves the engine untouched but needs a
   two-step creation and leaves the flat path's gap, against dec-A105's every call that re-anchors.
5. Saved development-build states of a sampler drawing sigma or an sd-spelled k. Recommended: no
   migration, per the state registry's pre-release rule; such a state restarts with sigma and k misread
   by the factors above. The alternative, a state attribute marking the old meanings, costs about 60
   lines and serves only states written by builds that never shipped.
6. A warm start's conversion of a linear leaf's slopes between two covariate standardizations (the
   leaf-conversions rule). dec-B418 names the response's mapping; the covariate standardization is the
   predictors' frame, like the cut grid a warm start already remaps. Recommended: unchanged, outside
   this slice, and named to the maintainer only if a reading of dec-B418 as "every frame" is wanted.
   The alternative removes it too (about 60 lines and the warm-start half of test-leaf-conversions.R).

## Evidence

Ran on 4cf51b11 and on the rebased tips (private library scratch/libs/kint; probes and outputs in
scratch/kint/): today's behaviour per operation (01-today.R, 02-today-more.R), the emulation of the new
sd spelling (03-emulate.R), an install across a re-derived mapping with and without conversion
(04-no-convert.R), and the greps for sd-spelled tests, conversion tests, baselines, snapshot files and
exact gates. Read: every engine, bridge, R and header claim cited by symbol; stan4bart's bartcore branch
at 963956b, its bart_args forwarding, its flat-API calls, its per-chain warm-up re-derivation and its
restore; bartCause, treatSens and bairrtt for their k use and installs; the setState plan at 20e08954 and
the response-scale-rows plan at faf8ca7d; the rulings in the table. New values in the operation table
are arithmetic on the run values, not runs, except the emulated install. The blind critique's runs
(scratch/kcrit2/, private library scratch/libs/kcrit2 on the same code): the probes re-run identically,
the emulation over every family and leaf model, sigma across a re-anchor and an install, xbart's warm
path, and the held sd's one-ulp miss; and read: the re-creation path through Chain::setState's restore.
