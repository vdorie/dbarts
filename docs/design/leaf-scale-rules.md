# Leaf scale rules: k or the sd

Status: FOR DECISION, 2026-10-09. Nothing here is built or changed; it lays out
the rules as they stand and two ways to make them one rule.

## 1. The question

The leaf sd is the prior sd of the forest's total at a point; each leaf's own
prior sd is that divided by the square root of the number of trees. The
spelling of a leaf prior is whether it is written with k, `normal(k = )`, or
with the sd, `normal(sd = )`; either is fixed or drawn (`chi()` on k,
`invchi()` on the sd). Whatever the spelling, each chain holds one number, k,
and

    leaf sd = k.scale / k

`k.scale` depends on the spelling. With k it is the data's scale: half the
training range of a continuous response, 3 for probit, pi sqrt(3) for
logistic. With `sd = s` or `sd = invchi(df, s)` it is 2 s, so k sits at 2 when
the leaf sd is s; with `sd = invchi(df, 0)` it is the data's scale.

A state is a chain's current trees, sigma and k. Installing one means
`setState`, `copy()`, a reload with `readRDS`, or a warm start
(`bart(warm.start = )`, `installTrees`). A re-anchor is a call that re-derives
the data's scale from a new response: `setResponse` or `setOffset` with
`updateScale = TRUE`, or `setData`.

While `k.scale` stays put, keeping k and keeping the leaf sd are the same. They
part when an operation changes `k.scale` under a chain: a new spelling or
`invchi()` scale, a state installed under another prior, a re-anchor, or a
state installed into a sampler whose response has another spread. The rulings
so far keep k in some of these and the sd in others. This note sets them out
and asks which number is the chain's parameter.

Keeping the sd means, here, keeping it relative to the leaves it governs. Where
an operation leaves the leaves as they are, the sd stays the same number.
Where it rescales them - a state converted into a sampler on another response,
or a re-anchor, which stretches the live fit with the response - the sd is
rescaled with them.

## 2. The rules today

Measured on the current build with a gaussian sampler whose response has half
range 3.71, 20 trees, after 200 sweeps under `normal(k = chi(1.5, 2))`: one
chain at k = 2.86, leaf sd 1.30. Built: the current code does this (run).

| operation | k | leaf sd | built |
|---|---|---|---|
| creation, fixed `k = c` | c | `k.scale` / c | yes |
| creation, drawn `k = chi()` | 2, whatever the chi scale | `k.scale` / 2 | yes |
| creation, `sd = s` or `sd = invchi(df, s)`, s > 0 | 2 | s | yes |
| `setLeafPrior` into a fixed value, either spelling | the stated k; 2 for a stated sd | as stated | yes |
| `setLeafPrior` into a drawn prior, from any prior | old k x new `k.scale` / old `k.scale` | kept (1.30 stays 1.30); the new prior acts from the next draw | yes |
| `setModel` with another leaf prior | kept | 1.30 to 0.70 into `sd = invchi(3, 1)` | yes, note a |
| install of a state from the same sampler and prior | as stored | as stored | yes, bit for bit |
| install of a state drawn under another prior, or before a `setLeafPrior` | as stored | the recipient's `k.scale` / k: 1.30 to 0.70 | yes |
| any install into a sampler that holds k or the sd fixed | the sampler's | the sampler's | yes |
| `copy()` of a sampler with no stored state | stores the current state, then copies it | as stored | yes |
| `setState`, `copy()` or a reload into a sampler whose response has another spread | re-expressed against the recipient's `k.scale` | kept | no, note b |
| warm start, across a change of prior or of response spread | as stored | the recipient's `k.scale` / k: 1.30 to 0.70 across a prior, 1.30 to 3.90 across a response three times as wide | yes, note c |
| re-anchor, k spelling | kept | stretches with the live fit (x3 on 3 y) | yes |
| re-anchor, sd spelling | kept | the named sd, fixed or drawn, stays in response units while the live fit stretches x3 | yes |
| `setResponse` or `setOffset` at the default `updateScale = FALSE` | unchanged | unchanged | yes, note d |
| probit rescaling step | divided by a drawn factor each sweep | multiplied by it | yes |

The probit rescaling step is a move, taken each sweep by a probit fit with one
forest, constant leaves and a drawn k, that multiplies the latent responses and
the leaves by one drawn factor and divides k by it.

- a. No ruling names `setModel`. The rule of the spelling-switch rulings,
  "Changing the prior never moves the state", may already cover it, in which
  case keeping k there is a defect against that ruling.
- b. Ruled, not built. Today k installs as stored: under the k spelling the sd
  triples (1.30 to 3.90) on a response three times as wide, while the converted
  leaves keep the donor's function. Under the sd spelling `k.scale` does not
  depend on the response, so the sd is already kept.
- c. The warm-start ruling carries the donor's k where the new fit draws k. An
  earlier call put the warm start under the same install rule as `setState`,
  and the response-spread ruling revised that call, so whether its planned
  build covers the warm start is open.
- d. 0.9-34's `setResponse` re-derived the scale on every call.

A fit stores k and its draws in `fit$k`, with `k.scale` in `fit$leaf.prior`;
`extract(fit, "leaf.prior.sd")` computes `k.scale` / k; `getK` returns k.

The rules rest on these words of yours:

- Offered keeping k or keeping the sd across a switch between the k and sd
  spellings of a drawn prior: "B. Keep the spread."
- Offered keeping the sd at the call or letting it jump, for a changed
  `invchi()` scale: "OK, then A. Keep the spread at the call."
- On a fixed sd made drawn: "Changing the named value within the sd spelling
  should be solved by the user's intent. They're trying to actively set the sd
  of the leaf prior, so both of their values should be interpretted as such."
- On a switch out of a fixed prior: "if a person installs a fixed value of k,
  they clearly mean it to be interpretted literally so that should also reset
  the k.scale" and "I don't see changing the prior as changing the state / the
  active parameter. For that, you'd need a setK() function."; then, offered
  keeping the sd until the next draw or keeping k's number: "Yes, keep until the
  next draw."
- Offered keeping k on a restore across a change of prior, or the state
  carrying the sd so a restore puts it back: "The k itself should literally
  transfer - a parameter is a parameter, regardless of the prior. It may be a
  bad fit, but so be it."
- Offered keeping the response-spread rule for a change of response scale, or
  k as stored there too: "Keep it for scale changes. If it helps at all, I
  guess we can think of `k` as just the internal representation and the
  parameter as the sd. This whole thing has been and continues to be a mess."
- Offered whether a restore reinstates a value the sampler holds fixed: "Model
  when fixed, state when drawn."
- Offered the warm start adopting the donor's prior or keeping the new fit's:
  "A. The donor gives its trees, and sigma and k only where the new fit draws
  them (as built)."
- On naming `k.scale`: "I want k to be the engine's k always." and
  "Ultimately, the exact value isn't itself important. The only thing that
  really matter is that users have a way to convert between k and standard
  deviation."

## 3. Disagreements

Two rulings say a change of prior leaves the parameter alone, and mean
different numbers by it. Of `setLeafPrior` you said "I don't see changing the
prior as changing the state / the active parameter", and the rule keeps the
leaf sd. Of a restored state you said "a parameter is a parameter, regardless
of the prior", and the rule keeps k.

Sorted by what each rule keeps, keeping the sd meant as in section 1:

| keeps the leaf sd | keeps k |
|---|---|
| `setLeafPrior` into a drawn prior (spelling-switch rulings; built) | an install across a change of prior (the ruling that a parameter is a parameter; built) |
| an install into a sampler on another response spread (response-spread ruling; not built) | `setModel` with another leaf prior (built; may be a defect, note a) |
| | the warm start, across a prior and across a response spread (built; open, note c) |
| | a re-anchor with a drawn sd: `k.scale` and k do not move, so the sd stays one number while the leaves stretch (named-sd ruling; built) |

A re-anchor under the k spelling keeps both: k stays and the sd stretches with
the leaves.

Where a user meets the difference:

1. Store, change the prior, restore. At sd 1.30 under `normal(k = chi(1.5,
   2))`, `setLeafPrior(normal(sd = invchi(3, 1)))` keeps 1.30. Restoring the
   state stored a moment before, at 1.30, gives 0.70.
2. Two writers of one prior. `setLeafPrior` and `setModel` given the same new
   leaf prior leave the chain at 1.30 and at 0.70.
3. Two installers of one state. Once the response-spread ruling is built, if it
   does not reach the warm start, `setState` and a warm start put the same state
   into the same sampler on a response three times as wide at 1.30 and at 3.90.
4. A drawn sd across a re-anchor. The leaves stretch with the response and the
   sd they were drawn under does not, so the two differ by a factor equal to the
   ratio of the spreads, 3 here. Under the k spelling they move together.

## 4. Two consistent rules

Both views take a fixed value as stated, agree wherever `k.scale` does not
move, and leave the probit rescaling step as it is, since it moves the leaves
and the sd together. They differ in the rows below. Each cell says what a user
sees and which ruling it reverses, if any.

| operation | the sd is the parameter | k is the parameter |
|---|---|---|
| `setLeafPrior` into a drawn prior | sd kept, as built | k kept: the sd jumps at the call by the ratio of the two `k.scale` values (1.30 to 0.70 in the example). Reverses the spelling-switch rulings |
| `setModel` with another leaf prior | sd kept: 1.30 stays 1.30. Changes built behaviour, which may be a defect today (note a) | k kept, as built |
| an install across a change of prior | the state's sd: a state saved at 1.30 restarts at 1.30 under any prior, so a saved state records its sd. Reverses the ruling that a parameter is a parameter | k as stored, as built and ruled |
| an install into a sampler on another response spread | sd kept, as ruled | k as stored: under the k spelling the sd restarts as many times larger as the response is wider, while the leaves keep the donor's function. Reverses the response-spread ruling and its confirmation |
| warm start | the donor's sd, across a prior and a response spread. Reverses the warm-start ruling's as-built carrying of k | the donor's k, as built |
| re-anchor with a drawn sd | the drawn sd stretches with the leaves, the named `invchi()` scale staying in response units. Reverses the named-sd ruling for the value in force | the sd stays in response units, as built |
| what a fit reports | either. Reporting the sd as the parameter puts its draws in the fit where `fit$k` is now, with `extract(fit, "k")` computing k, and reverses the ruling that no fit stores the leaf sd | `fit$k` holds k's draws, `extract(fit, "leaf.prior.sd")` computes the sd, as built |

The sd view keeps the spelling-switch rulings and the response-spread ruling.
It reverses the ruling that a parameter is a parameter, the warm start's
carrying of k, and the named-sd ruling for a drawn sd at a re-anchor; it
changes `setModel`; and, if the report follows, it reverses the ruling that no
fit stores the leaf sd.

The k view keeps every built behaviour but `setLeafPrior`. It reverses the
spelling-switch rulings and the response-spread ruling with its confirmation.
A rule relative to the data throughout would also leave a moved state's leaves
unconverted, as 0.9-34 did, against your "Convert on install."

## 5. Other packages

The parameter column is the number each package holds for the leaf scale, the
one it would carry if a chain were moved.

| package | leaf scale written as | when the prior, a stored state or the data's scale changes | the parameter |
|---|---|---|---|
| BART 2.9 (`wbart`, `pbart`, `gbart`) | k, a fixed number, against half the range, or against 3 for probit and 6 for logistic; `wbart` also takes `sigmaf`, an sd of f | set once per call; no sampler object, no restart, no hyperprior; `sigest` sets the sigma prior once | the per-tree sd, from k or `sigmaf` |
| BayesTree | k, after rescaling y to [-0.5, 0.5] by the training range | one call; nothing changes | k |
| bartMachine 1.4 | k (cross-validated over 2, 3, 5) on the same rescaling | serialized only to predict; no continuation | k |
| bcf 2.0 | `sd_control`, `sd_moderate` in response units (2 sd(y), sd(y)) | no warm start; scale multipliers drawn under half-Cauchy and half-normal priors | the sd |
| stochtree 0.4 | a leaf variance per tree on the standardized outcome, `sigma2_leaf`, drawn under an inverse gamma prior | a warm start installs the earlier fit's trees and, where the new fit draws it, its leaf variance as stored, whatever the new prior; values are carried in standardized units, so across data of another spread they are relative to each fit's sd(y) | the leaf variance, relative to sd(y) |
| SoftBart 1.0 | k at the call, turned at once into an sd, `sigma_mu`, drawn under a half-Cauchy | its forest object for loops holds `sigma_mu`, takes the response in whatever units it is handed and never rescales; the prior has no setter; k is reported as a function of `sigma_mu` | the sd |
| dbarts 0.9-34 (run) | k, fixed or `chi()` | `setResponse` re-derived the scale on every call, keeping k and stretching the live fit; a state carried k always and installed it as stored, even into a fixed-k sampler; a state moved to a response three times as wide kept k and its leaves were read in the new units, the function tripling | k, relative to the current response |

No other package's fitting functions change a leaf prior on a live chain;
stochtree's low-level interface lets a hand-written loop set the leaf variance
directly. The two that restart a chain, stochtree and 0.9-34, carry the stored
value as stored, relative to the data's scale. bcf names the sd and stochtree a
variance. SoftBart's user writes k, but the sampler holds and reports the sd,
deriving k. BART, BayesTree and bartMachine name k, BART's `wbart` also taking
an sd.

## 6. Uses in the wild

GitHub code search, the dbarts issue tracker and a web search found about a
dozen uses of the sampler in loops outside CRAN packages, listed in Appendix B.
What they assume:

- All use a fixed k, most through 0.9-x's `node.prior = normal(k)`. None uses a
  drawn k in a loop, `setLeafPrior`, a state across a change of prior, or a
  state moved between responses. No current or ruled leaf rule changes any of
  them.
- Most call `setResponse` every iteration on a latent or working response (an
  ordinal latent, a Polya-Gamma working response, a SUR or VAR residual). Under
  0.9-34 each call re-derived the scale, so the leaf sd tracked each draw's
  range and the live fit was stretched each time; under 1.0 the scale stays as
  created. One user asked for the 1.0 behaviour in a dbarts issue of July 2026
  ("editing the prior every MCMC iteration (which is incorrect)"), and code in
  a related repository routes the response through the offset, citing that
  issue.
- Users reason in sds. The same user holds a fixed reference response at
  creation, passing the real one through `setOffset(updateScale = FALSE)` so
  the leaf prior never moves; with the reference spanning plus and minus
  sqrt(G0), its half range is sqrt(G0), so k = 1 gives a forest sd of sqrt(G0).
  Five related VAR and quantile repositories pass k under the name `sd.mu`.
  Under 1.0 both would write `normal(sd = )`.
- Two restore stored states, each into a sampler on the same model and data
  with k fixed, so the leaf rules leave them alone. One, across processes
  through `saveRDS`, checks predictions bit for bit. The other edits a chain
  state with `chain_state@savedTrees <- integer(0)`, which errors under 1.0,
  where a chain state is a list.

Most at risk is the negative-binomial loop that passes a Polya-Gamma working
response, (y - xi) / (2 omega), to `setResponse` each iteration under a fixed
k. As written it does not run under 1.0: its creation call gives
`resid.prior = fixed(1)` with `sigma = 1`, which now stops with "'sigma' has no
effect under a fixed residual scale"; the reference-response loop above makes
the same call. Once the user drops `sigma = 1`, its posterior changes. Under
0.9-34 its leaf sd followed each draw's range, which the smallest omega
dominates; under 1.0 it is fixed at the range of the one draw made before the
sampler was created. Neither is the log-mean scale of the 1.0
negative-binomial family (`k.scale` 3). This comes from 1.0 no longer
re-deriving the scale on `setResponse`, not from the leaf rules weighed here.

On k in general, Stack Overflow and Cross Validated turned up nothing. What
turned up is the naming above (k passed as `sd.mu`, k worked back from a wanted
sd) and the issue on `setResponse` moving the prior.

## Appendix A. Limits of the evidence

- stochtree, SoftBart, bcf and BayesTree are not installed here; their rows
  come from reading their CRAN sources, not from runs. BART, bartMachine and
  dbarts 0.9-34 were run or their installed code read.
- Code search reaches only public, indexed repositories; replication archives
  outside GitHub and private code are not covered.
- The uses in the wild were read, not run, apart from the two 1.0 failures
  quoted, which were run.

## Appendix B. Evidence

Rulings, in the decision register:

- Spelling-switch rulings (keep the sd; "Changing the prior never moves the
  state"): dec-B356, dec-B369, dec-B392, dec-B393.
- The ruling that a parameter is a parameter: dec-B401.
- The response-spread ruling: dec-B384, confirmed by dec-B407; planned, not
  built.
- The earlier call putting the warm start under the install rule: dec-A146.
- `k.scale` named, k reported as is: dec-B201; the reader: dec-B141.
- The named-sd ruling (k or sd spelling, a named sd absolute across a
  re-anchor): dec-A105.
- Fixed values are model, drawn values state: dec-B195, dec-B196.
- Convert on install: dec-B200.
- The warm-start ruling: dec-B343.
- `copy()` stores the current state: dec-B396.
- `setModel` changes parameters, never structure: dec-B254.
- No fit stores the leaf sd: dec-B376.
- The probit rescaling step: dec-B371, dec-A190.

Uses found in the wild:

- Reference response, k from a target sd, and the `sigma = 1` creation call:
  https://github.com/EoghanONeill/TobitBART
- The same author: https://github.com/EoghanONeill/SURBART,
  https://github.com/EoghanONeill/BAVART
- dbarts issue "setResponse updates the prior scale":
  https://github.com/vdorie/dbarts/issues/80
- Mixed-frequency and quantile VARs (`sd.mu` as k; the offset comment):
  https://github.com/mpfarrho/mf-bavart, https://github.com/mpfarrho/qf-bart,
  https://github.com/mpfarrho/gp-mf; a copy in
  https://github.com/danielvitonet/TESI
- Negative-binomial BART: https://github.com/jacobenglert/nbbart
- Principal stratification with BART (arXiv 2408.03777):
  https://github.com/AlkemaLab/prince_BART, https://github.com/AlkemaLab/nkids
- https://github.com/bilalafzalshafi/sebart
- Sensitivity analysis restoring a reference chain:
  https://github.com/cochran4/wadsi
- State round trip across processes from Python:
  https://github.com/harry1310/WeatherProbabilistic
- Python wrappers for benchmarks, no leaf-scale assumptions:
  https://github.com/Gattocrucco/bart-gp-article

Other packages' sources: BART 2.9.10 and bartMachine 1.4.2 (installed); CRAN
sources of stochtree 0.4.5, SoftBart 1.0.3, bcf 2.0.2 and BayesTree 0.3-1.5.
