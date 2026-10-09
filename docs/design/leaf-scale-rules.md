# Leaf scale rules: k or the sd

Status: FOR DECISION, 2026-10-09. Nothing here is built or changed; it lays out
the rules as they stand and two ways to make them one rule.

## 1. The question

A leaf prior is written as `normal(k = )` or `normal(sd = )`, fixed or drawn
(`chi()` on k, `invchi()` on the sd). Either way the sampler holds one number
per chain, k, and the leaf sd in force is

    sd = k.scale / k

where `k.scale` depends on the spelling: under `k =` it is the data's scale
(half the training range for a continuous response, 3 for probit, pi sqrt(3)
for logistic); under `sd = s` or `sd = invchi(df, s)` it is 2 s, so that k
sits at 2 when the sd is s.

So long as `k.scale` stays put, keeping k and keeping the sd are the same
thing. They part only when an operation changes `k.scale` under a chain: a
change of spelling or of the `invchi()` scale, a state installed under another
prior, or a change of the response's spread under the `k` spelling. The
rulings so far answer those cases differently, some keeping the sd and some
keeping k. This note sets them side by side and asks which quantity is the
chain's parameter.

## 2. The rules today

Measured on the current build with a gaussian sampler whose response has half
range 3.71, 20 trees, after 200 sweeps under `normal(k = chi(1.5, 2))`: one
chain at k = 2.86, sd = 1.30. "Built" means the behaviour was run and matches.

| operation | k | leaf sd | built |
|---|---|---|---|
| creation, fixed `k = c` | c | `k.scale` / c | yes |
| creation, drawn `k = chi()` | starts at 2, whatever the chi scale | `k.scale` / 2 | yes |
| creation, `sd = s` or `sd = invchi(df, s)` | 2 | s | yes |
| `setLeafPrior` into a fixed value, either spelling | the stated k, or 2 for a stated sd | as stated | yes |
| `setLeafPrior` into a drawn prior, from anything (spelling switch, new `chi()` or `invchi()` scale, fixed to drawn) | re-expressed: old k times new `k.scale` / old `k.scale` | kept at the call (1.30 stays 1.30); the new prior acts from the next draw | yes |
| `setModel` with another leaf prior | kept as is | moves by the ratio of the two `k.scale` values (1.30 to 0.70 into `sd = invchi(3, 1)`) | yes; not ruled |
| `setState`, `copy()`, reload of a state from the same sampler and prior | as stored | as stored | yes, bit for bit |
| `setState` of a state drawn under another prior (or before a `setLeafPrior`) | installed as stored | recipient's `k.scale` / k (1.30 becomes 0.70) | yes |
| any install into a sampler that holds k or the sd fixed | the sampler's fixed value; the state's k is ignored | the fixed value | yes |
| `copy()` of a sampler with no stored state | stores the current state first, then as above | as above | yes |
| a state installed into a sampler whose response has another spread | ruled: re-expressed so the sd is kept, as the leaves are converted to keep the function | ruled: kept | no: today k installs as stored, so under the `k` spelling the sd triples (1.30 to 3.90) on a response three times as wide while the converted leaves keep the donor's function; under the `sd` spelling `k.scale` does not depend on the response and the sd is already kept |
| warm start (`bart(warm.start = )`, `installTrees`), where the new fit draws k | donor's k as stored, across a change of prior or of response spread | follows the recipient's `k.scale` (1.30 to 0.70 across a prior; 1.30 to 3.90 across a threefold response) | yes; the response-spread case is not named by any ruling |
| re-anchor (`setResponse` or `setOffset` with `updateScale = TRUE`, `setData`), `k` spelling | kept | stretches with the response (x3), as the live fit does | yes |
| re-anchor, `sd` spelling | kept | a named sd stays in response units, fixed or drawn, while the live fit stretches x3 | yes |
| `setResponse` or `setOffset` at the default `updateScale = FALSE` | unchanged | unchanged | yes; 0.9-34's `setResponse` re-anchored on every call |
| probit rescaling step (probit, one forest, constant leaves, k drawn) | divided by a drawn factor each sweep | multiplied by it, with the leaves and latents | yes |
| what a fit reports | `k` and its draws; `getK` the engine's k | computed: `extract(fit, "leaf.prior.sd")` is `k.scale` / k | yes |

The words the rules rest on, in the maintainer's own:

- On the spelling switches (four rulings, 2026-10-08): "B. Keep the spread.";
  "OK, then A. Keep the spread at the call."; "Changing the named value within
  the sd spelling should be solved by the user's intent. They're trying to
  actively set the sd of the leaf prior, so both of their values should be
  interpretted as such."; "if a person installs a fixed value of k, they
  clearly mean it to be interpretted literally"; "I don't see changing the
  prior as changing the state / the active parameter. For that, you'd need a
  setK() function."; "Yes, keep until the next draw."
- On a state across a change of prior (2026-10-09): "The k itself should
  literally transfer - a parameter is a parameter, regardless of the prior. It
  may be a bad fit, but so be it."
- On a state across a change of response spread (2026-10-09, confirming the
  ruling of 2026-10-08): "Keep it for scale changes. If it helps at all, I guess
  we can think of `k` as just the internal representation and the parameter as
  the sd. This whole thing has been and continues to be a mess."
- On fixed values: "Model when fixed, state when drawn."
- On the warm start: "A. The donor gives its trees, and sigma and k only where
  the new fit draws them (as built)."
- On what k is: "I want k to be the engine's k always." and "Ultimately, the
  exact value isn't itself important. The only thing that really matter is that
  users have a way to convert between k and standard deviation."

## 3. Where the rules disagree

Two rulings both say a change of prior leaves the parameter alone ("I don't see
changing the prior as changing the state / the active parameter" and "a
parameter is a parameter, regardless of the prior"). They differ on which
number the parameter is: `setLeafPrior` keeps the sd, `setState` keeps k.

| rules that keep the sd | rules that keep k |
|---|---|
| `setLeafPrior` into any drawn prior (the four spelling rulings) | `setState`, `copy()`, reload across a change of prior |
| a state moved to another response spread (ruled, not built) | `setModel` with another leaf prior (built, not ruled) |
| a named sd held in response units across a re-anchor | the warm start, across both a prior and a response spread (built) |
| | the re-anchor under the `k` spelling (k kept, sd follows the data) |
| | reporting: draws of k stored, the sd computed |

Where that shows:

1. Store, change the prior, restore. Under `normal(k = chi(1.5, 2))` at sd 1.30,
   `setLeafPrior(normal(sd = invchi(3, 1)))` keeps 1.30; restoring the state
   stored a moment before gives 0.70, though the state was taken at 1.30.
2. Two writers of one prior. `setLeafPrior` and `setModel` given the same new
   leaf prior leave the chain at 1.30 and 0.70.
3. Two installers of one state. Once the response-spread ruling is built,
   `setState` and the warm start will place the same state into the same
   sampler on a wider response at different sds (1.30 against 3.90 in the
   example).
4. A drawn sd across a re-anchor. Under `sd = invchi()`, a re-anchor stretches
   the live leaves with the response but holds the drawn sd in response units,
   so the leaves and the sd they were drawn under part by the ratio of the
   spreads. Under the `k` spelling the two move together.

## 4. Two consistent rules

The views agree wherever `k.scale` is unchanged, and on fixed values (taken
as stated) and the probit step (which moves the leaves and the sd together).
They differ only in the rows below.

| operation | sd is the parameter | k is the parameter |
|---|---|---|
| `setLeafPrior` into a drawn prior | keep the sd (as built) | keep k; the sd moves by the ratio of `k.scale` values. Reverses the four built spelling rulings |
| `setModel` with another leaf prior | keep the sd. Changes built behaviour | keep k (as built) |
| `setState`, `copy()`, reload across a change of prior | install the state's sd; the state must record the sd, or the `k.scale` its k was drawn against. Reverses the "a parameter is a parameter" ruling; the option shown then was about 90 lines against 15 | k as stored (as ruled) |
| a state into a sampler on another response spread | keep the sd (as ruled, not built) | k as stored; under the `k` spelling the sd jumps by the ratio of spreads while the leaves keep the function. Reverses the response-spread ruling and its confirmation; drops unbuilt work. A fully k-relative rule would also leave the leaves unconverted, as 0.9-34 did, against "Convert on install." |
| warm start | donor's sd. Changes built behaviour | donor's k (as built) |
| re-anchor under a drawn `sd` | the drawn sd stretches with the leaves, k re-expressed. Changes built behaviour | k kept, the sd held in response units (as built) |
| what a fit reports | the sd and its draws first; k as the internal number behind `getK` and `k.scale`. Changes what a fit stores | as built |

In short: the sd view keeps every ruling of 2026-10-08 and the response-spread
confirmation, and moves `setModel`, the warm start, the state record, the drawn
sd at a re-anchor, and perhaps the report. The k view keeps every built
behaviour except `setLeafPrior`, which returns to keeping k, and drops the
unbuilt response-spread ruling.

## 5. What other packages do

| package | leaf scale named as | when the prior, a stored state or the data's scale changes | the parameter |
|---|---|---|---|
| BART 2.9 (`wbart`, `pbart`, `gbart`) | k (fixed number) against half the range, or 3 or 6 on a latent scale; `wbart` also takes `sigmaf`, an sd of f | computed once per call; no sampler object, no restart, no hyperprior; `sigest` sets the sigma prior once | the per-tree sd, from k or `sigmaf` |
| BayesTree | k, after rescaling y to [-0.5, 0.5] by the training range | nothing to change: one call | k, relative to the training range |
| bartMachine 1.4 | k (cross-validated over 2, 3, 5) on the same rescaling | serialized only to predict; no continuation | k |
| bcf 2.0 | `sd_control`, `sd_moderate` in response units (2 sd(y), sd(y)) | no warm start; scale multipliers drawn under half-Cauchy and half-normal priors | the sd |
| stochtree 0.4 | a leaf variance per tree on the standardized outcome, `sigma2_leaf`, drawn under an inverse gamma prior | a warm start installs the previous fit's trees and, when the new fit draws it, its leaf variance as stored, whatever the new prior; values are carried in standardized units, so across data of another spread they are relative to each fit's sd(y) | the leaf variance, relative to sd(y) |
| SoftBart 1.0 | k at the call, turned at once into `sigma_mu`; drawn under a half-Cauchy | its forest object for embedding in a loop holds `sigma_mu` as state, takes the response in whatever units it is handed and never rescales; no setter for the prior; k reported as a function of `sigma_mu` | the sd |
| dbarts 0.9-34 (run) | k, fixed or `chi()` | `setResponse` re-anchored on every call, k kept and the live fit stretched; a state carried k always and installed it as stored, even into a fixed-k sampler; a state moved to a response three times as wide kept k and its leaves were read in the new units, the function tripling | k, relative to the current response |

No other package's fitting functions change a leaf prior on a live chain;
stochtree's low-level interface lets a hand-written loop set the leaf variance
directly, the variance being then the number held. The two that restart a
chain, stochtree and 0.9-34, carry the stored value as stored, relative to the
data's scale. The packages written after BART name the scale
as an sd or variance (bcf, stochtree, SoftBart); SoftBart reports k only as a
function of the sd.

## 6. What users expect

GitHub code search, the dbarts issue tracker and a web search found about a
dozen uses of the sampler in loops outside CRAN packages (listed in Appendix
B). What they assume:

- All use a fixed k, most through 0.9-x's `node.prior = normal(k)`. None uses a
  drawn k in a loop, `setLeafPrior`, a state across a change of prior, or a
  state moved between responses. No current or ruled leaf rule changes any of
  them.
- Most call `setResponse` every iteration on a latent or working response (an
  ordinal latent, a Polya-Gamma working response, a SUR or VAR residual). Under
  0.9-34 each call re-anchored, so the leaf sd tracked each draw's range and
  the live fit was stretched each time; under 1.0 the scale stays as created.
  Two users say the 1.0 behaviour is the one they expected: the dbarts issue
  of July 2026 asking that `setResponse` not re-derive the prior ("editing the
  prior every MCMC iteration (which is incorrect)"), and a code comment that
  routes the response through the offset "since setResponse rescales the
  prior".
- Users reason in sds. One author holds a fixed reference response at
  creation and passes the real response through `setOffset(updateScale =
  FALSE)` so the leaf prior never moves, and computes k from a target sd (k =
  1 against a reference spanning plus and minus sqrt(G0), for a prior variance
  G0). Five related VAR and quantile codes pass k under the name `sd.mu`. Under 1.0 both
  would write `normal(sd = )`.
- Two restore stored states: one into a hand-built one-chain sampler on the
  same model and data, one across processes via `saveRDS` and a scaffold fit on
  the same data, checked bit for bit. Both hold k fixed and stay in one set of
  units, so they install bit for bit under every rule here.

Most at risk is the negative-binomial loop that passes a Polya-Gamma working
response, (y - xi) / (2 omega), to `setResponse` each iteration under a fixed
k. Under 0.9-34 its leaf sd followed each draw's range, which the smallest
omega dominates; under 1.0 it is fixed at the range of the single draw made
before the sampler was created. Its posterior changes, and neither scale is the
log-mean scale the 1.0 negative-binomial family uses (`k.scale` 3). This rests
on the response-scale lock, not on the leaf rules this note weighs.

Confusion about k in general: none found on Stack Overflow or Cross Validated.
What turned up is the naming above (k passed as `sd.mu`, k back-solved from a
wanted sd) and the issue on `setResponse` moving the prior.

## Appendix A. Not verified

- stochtree, SoftBart, bcf and BayesTree are not installed here; their rows come
  from reading their CRAN sources, not from runs. BART, bartMachine and dbarts
  0.9-34 were run or their installed code read.
- Code search reaches only public, indexed repositories; replication archives
  outside GitHub and private code are not covered.
- The in-the-wild uses were read, not run.

## Appendix B. Evidence

Rulings, in the decision register:

- Spelling switches keep the sd: dec-B356, dec-B369, dec-B392, dec-B393.
- k transfers literally across a change of prior: dec-B401.
- A state across a change of response spread keeps the sd: dec-B384, confirmed
  by dec-B407; its build is the TODO item state-install-keeps-spread.
- k is the engine's k, `k.scale` named: dec-B201; reader shape dec-B141.
- The k or sd spelling, a named sd absolute across a re-anchor: dec-A105.
- Fixed values are model, drawn values state: dec-B195, dec-B196.
- Convert on install: dec-B200; the install calls dec-A146.
- Warm start keeps the new fit's prior: dec-B343.
- `copy()` stores the current state: dec-B396.
- `setModel` changes parameters, never structure: dec-B254.
- What a fit stores: dec-B376.
- The probit rescaling step: dec-B371, dec-A190.

Uses found in the wild:

- TobitBART (fixed reference response, k from a target sd):
  https://github.com/EoghanONeill/TobitBART
- SURBART and BAVART, same author: https://github.com/EoghanONeill/SURBART,
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

Other packages' sources: BART 2.9.10 and bartMachine (installed); CRAN sources
of stochtree 0.4.5, SoftBart 1.0.3, bcf 2.0.2 and BayesTree 0.3-1.5.
