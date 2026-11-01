# forests-shared-total: several forests with no basis share one prior total

Status: BUILT on this branch 2026-10-31, not landed. Slices B, 1 (the engine and bridge), 2 (the R surface), 3 (the
constant basis column, dec-B528, dec-B543) and F (the default flip) are built, each reviewed and its findings fixed;
no call is held. Slice F's fix round changes engine behavior in one place, the median of `student()` with no sd under
`normal(sd = invchi())` ([Calls made](#calls-made) 32), so a fresh review of that round comes before the landing; the
[Landing note](#landing-note) holds the rest. Planned 2026-10-29; amended the same day after its blind critique
(build after amendments; scratch/crit-2plans/report.md) and the maintainer's rulings dec-B503 (the tree count),
dec-B505 (an equal part per forest) and dec-B506 (`blocks()` changed to match, slice B); amended again the same day
for dec-B508 (`forest()` takes no `leaf.prior`), dec-B509 (a forest with no basis takes the fitting function's leaf
prior, held, in every model: D1 revised, the default flip in [Slice F: the default flip](#slice-f-the-default-flip))
and dec-B510. `amplitude = fixed()` on a lone forest with no basis is ruled, taken in slice F (dec-B513). Rulings:
dec-B292, dec-B503, dec-B505, dec-B506, dec-B509, dec-B513, dec-B528 and dec-B543 (with dec-B291, dec-B274, dec-B307,
dec-B508, dec-B510, dec-B525); window "Post-merge, pre-1.0" (VD 2026-10-07).

agent: one opus implementer for slice B, then slices 1 and 2, then slice F; one opus review of slice B, one over
slices 1 and 2 and one of slice F, each running the mutants below, the same reviewer checking the fix round unless
the fixes change engine behavior. Tier "Changes draws" ([Process by risk](README.md#process-by-risk)): slice B
changes the `blocks()` prior; slice 1 edits the engine's calibration map, k draw and rescaling step; slice 2
opens a new prior; slice F changes a default. Slice 3 is R surface, planned for one sonnet review; the second opus
review of slices 1 and 2 read it instead, its two findings fixed there.
rng: slice B POSTERIOR-CHANGING for a `blocks()` fit whose groups hold unequal tree counts, the default 75 trees
over two groups included, and bitwise neutral for every other fit (its gates in [Slice B](#slice-b)).
Slices 1 and 2 NEUTRAL, bit for bit, for every call main accepts: no such call has two forests with no basis, and
every new branch is taken only where two or more state no sd; a lone sharer beside a forest with a basis keeps main's
path until slice F. Gate: a
bitwise A/B against an e2f37f6d build, the equivalence trio 55, 19 and 11 identical (equivalence-32418b6e,
bcf-equivalence-20eed7ed, multinomial-equivalence-80b1c8d4), the four snapshot files unchanged, no re-record.
New behavior, with nothing on main to compare against: a model of several forests with no basis, a new prior
gated by new arms of two exact-posterior gates, a statistical equivalence against the one-forest `blocks()` model
it equals, and a design note. Slice F POSTERIOR-CHANGING (dec-B509) for every model of two or more forests in which a
forest with no basis states no sd or states sd with no `amplitude`: the default drawn coefficient becomes the fitting
function's leaf prior, held; under probit and logistic `amplitude = fixed()` with no sd moves from 1 unit held to the
drawn k. Its gates and the baselines it moves are in [Slice F: the default flip](#slice-f-the-default-flip).
Slice 3 draws nothing differently: a call it refuses ran on main, and is the ruling's change (dec-B528).
window: after the merge to main, before 1.0-0 (dec-B292, dec-B509). Slice F lands only after
coefficient-prior-opt-in and bartCause's `bcf()` change ([Interaction with the other pending plans](#interaction-with-the-other-pending-plans)).
budget: ~3570 lines in the tree. Slice B ~580, built (engine ~110, tests/cpp ~150, a new exact gate ~200,
tinytest ~70, help ~10, design note ~40). Slices 1 and 2 ~2340: engine ~400, bridge ~70, R ~280, help and NEWS
~130, tests/cpp ~400, tinytest ~500, exact-gate arms and the equivalence script ~400, docs and records ~160
(slices: ~870, ~1470). Slice F ~590: engine ~40, tests/cpp ~90, R ~70, tinytest ~260 (re-pins in about twelve
files ~200), gate scripts ~30, help and NEWS ~40, design notes and records ~60. Slice 3 ~60. The pooled k over
per-forest leaf scales, the rescaling step's K-forest reading and the K-forest chain's single-forest leaf prior
are each a new path, not a parameter; the shared coefficient's draw is gone (dec-B509, D3). Out of the tree,
~450 lines of scratch: the A/B scripts. Agent time, estimated from comparable landings and not measured: slice B
~7 agent-hours (spent); slices 1 and 2 ~24 (implementers ~15, one review and its fix round ~5, gates ~4); slice
F ~8 (implementer ~4, review with mutants ~1.5, gates and re-records ~2.5); slice 3 ~2. bartCause's change, ~1.5,
is its own.

## Summary for the maintainer

Today a model takes one forest with no basis. `y ~ forest(x1) + forest(x2) + forest(x3)` is refused, and the
way to write forests that simply add is a basis column of ones with a stated size. After this item that formula
fits. The forests with no basis share one prior size the way the trees of one forest do: together they have
the size a single such forest would have, so three forests of equal tree counts are each about 0.58 of it. A
forest with a basis keeps its own size, as you ruled.

"The size one such forest would have" is `bart()`'s own forest in every model, as you ruled (dec-B509): its `k`,
a stated `leaf.prior`, and under probit and logistic its drawn `k`. Without a forest that has a basis the model is
then exactly `bart()`'s one forest with its trees divided among the formula's terms, which is what `blocks()`
already does inside one forest. A drawn `k` is drawn once for the group, so the forests grow and shrink together,
as trees do. A forest that states its own `sd` has that size and leaves the share.

Beside a forest with a basis this changes today's default: the forest with no basis no longer draws a coefficient
under a half-Cauchy with median twice the response's sd. That prior becomes `amplitude = student(1)` with `sd`
the median, which bartCause's `bcf()` will state. The check you asked for ran first: on the bcf paper's
simulation design, 8 scenarios and 800 paired fits, holding the size against drawing it changed the 90 percent
coverage of the individual effects by +0.004 [0.001, 0.008] and of the average effect by -0.003 [-0.019, 0.014],
with a 3.7 percent smaller error in the effects, at equal time. Coverage of the prognostic function was about half
a point lower, and 0.7 to 1.5 points lower with about 3 percent more error where that function is nonlinear. The
flip is its own slice, about 590 lines and 8 agent-hours, and it lands after the `student()` prior and
bartCause's change.

The default sizes were checked by simulating the prior. Each forest at the full size instead makes the sum
sqrt(K) times as wide: at three forests under probit at a fixed size, 38 percent of the prior mass goes to
probabilities past 0.01 or 0.99, where one forest puts 14 percent.

You ruled the tree count: the forests with no basis split `bart()`'s count among themselves as `blocks()`
splits it, and a forest that states its own count keeps it. Print will show each forest's count and the total.
You ruled that each such forest takes an equal part of the shared size whatever its tree count, a count saying
how complex a term may be and not how much it explains; and that `blocks()` changes to match, each group an
equal part of its forest's size. The formula and the matching `blocks()` call are then one model at every count.

The `blocks()` change on its own: about 580 lines and 7 agent-hours, built first. It changes results only for
`blocks()` fits whose groups hold different numbers of trees. At the default split the change is small: 75 trees
over two groups, 38 and 37, move each group's prior size by under 1 percent, and over four groups by at most 2
percent. Where a caller set uneven counts it is large: 5 and 15 trees out of 20 move the groups from 0.50 and
0.87 of the forest's size to 0.71 each. `blocks()` was not in 0.9-34, so no released fit changes. The whole
item is about 3,570 lines and 41 agent-hours, these estimates from comparable work.

You ruled that a constant basis column, today's workaround, is now refused in favor of the plain spelling, sd stated
or not (dec-B528). You ruled that `amplitude = fixed()` on a lone forest with no basis, which after the flip only
restates its default, is taken there too, where it is refused today (dec-B513).

## Goal

A model may hold any number of forests with no basis, beside forests with a basis or alone, at every door that
takes several forests: a formula's `forest()` terms, a `forests` list and a data object's `bases`. The forests
with no basis that state no `sd` share one prior total, each an equal part of its variance whatever its tree
count (dec-B505); the total is the fitting function's leaf prior, held, in every model (dec-B509); a drawn k is one
draw for the group; the fitting function's tree count is the group's (dec-B503). A forest partitioned by `blocks()`
gives each group an equal part of its prior variance (dec-B506), so the formula and the matching `blocks()` call
are one model at every count. Every model main accepts draws what it drew, except `blocks()` fits with unequal
groups and, after slice F, models with a forest with a basis whose forest with no basis states no sd.

## What exists today

Read, and run on main at e2f37f6d (shipped build, arm64 macOS; the probes are in scratch/fst-plan, untracked):

- The refusal, gone with slice 2. [`plainForest`](../../R/model.R) refused a second forest with no basis ("a
  model has one forest with no multiplier"), and the formula reader refused it first (run: "the formula has 2
  forest() terms with no basis"), pinned in inst/tinytest/test-forest-defaults.R and
  inst/tinytest/test-formula-terms.R, whose pins are now creations. Plain terms beside a `forest()` with no
  basis were refused as two such forests; they are refused by their two readings now.
- The engine already takes any number of forests of kind none: [`ForestKind`](../../src/bartcore/combiner.hpp)
  is stated per forest and nothing in [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) counts them.
  Each forest has its own coefficient block, a drawn one with no basis under a half-Cauchy mixture.
  The cap is R's alone.
- Sizes. [`forestMapEntry`](../../src/bartcore/combiner.hpp) and
  [`defaultForestSd`](../../src/bartcore/combiner.hpp) give a forest with no basis, beside one with a basis, a
  drawn coefficient whose half-Cauchy median is 2 units under gaussian and 1 under probit and logistic, the
  forest's own fit at one unit (U, sd(y) or the link's error sd); held, a quarter of the response's range under
  gaussian ([`singleForestDefaultSize`](../../src/bartcore/chain.hpp), dec-B290) and U otherwise (dec-B293).
  Run: the reader gives the drawn one `sd` 1.67 (2U), `k.scale` 0.833 (U), k 1; held, `sd` and `k.scale` 0.827,
  a quarter of the range.
- A lone forest is the single-forest chain: `y ~ forest(x1)` reports k 2 and `k.scale` 1.65, half the range;
  `forest(x1, sd = 2)` is held as k = `k.scale` / 2. Under probit and logistic `bart()` draws k under
  chi(1.5, 2), `k.scale` 3 and 5.44 (run). dec-B420 keeps the drawn k: a fixed one was not good enough
  ([binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)).
- A model of several forests refuses at [`resolveSamplerSpec`](../../R/spec.R) every family but gaussian, probit
  and logistic, a k hyperprior, a k other than 2, test predictors, a variance forest, linear and gp leaves,
  DART, monotone and more (run: `bart(..., k = 3)` and `leaf.prior = normal(k = 3)` give "a treatment forest
  does not support a non-default 'k'"; `dbarts()` has no `k`, `dbarts(k = 3)` being an unused argument;
  `leaf.prior = normal(sd = 1)` is refused with a redirect, "a model of several forests states each forest's size
  on the forest: write the size of the forest with no basis as forest(x1 + x2, sd = 2) ... and leave 'sd' out of
  'leaf.prior'", re-run on 1ffd6f92's code). `leaf.prior = normal(2)`, the default restated, is taken: the exact
  gates bcf-exact.R, bcf-exact-weak.R and bcf-exact-restricted.R pass it beside forests that all state sd.
  [`resolveLeafHyperprior`](../../R/model.R) forces the binary default to a fixed 2 on a model of several forests,
  the map pinning k at 1.
- Who draws the default half-Cauchy beside a forest with a basis (dec-B509's audit): bartCause's `bcf()` states
  no `amplitude` and no sd on its prognostic forest where it draws it; bcf-equivalence.R's scenarios, the
  bart2twoforest scenario of equivalence.R, bcf-exact.R, bcf-exact-weak.R, bcf-latent-exact.R and sbc.R's bcf arm
  leave `amplitude` out on the forest with no basis (all but the first two state its sd); about twelve tinytest
  files pin that default or its refusals. R's `2 * sd(y)` and the engine's 2U differ in the last place
  (1.2166873523600565 against 1.2166873523600563, run): a forest stating `sd = 2 * sd(y)` under gaussian draws
  differently from the default, bit for bit; under probit `sd = 1` and under logistic `sd = pi / sqrt(3)` draw
  `identical()` to it.
- The check dec-B509 asked for, run before this amendment: a held prognostic forest against the drawn default on
  Hahn, Murray and Carvalho (2020) section 6.1's design, 8 scenarios, 800 paired fits, held minus drawn pooled:
  tau(x_i) 90 percent coverage +0.004 [0.001, 0.008], ATE coverage -0.003 [-0.019, 0.014], tau RMSE 3.7 percent
  lower, at equal time (ratio 1.011). Not clearly worse. f's coverage about 0.5 points lower pooled, 0.7 to 1.5
  points with about 3 percent higher RMSE where the prognostic function is nonlinear.
- The workaround: `forest(x2, basis = one, sd = 2)`, `one` a column of ones, is a numeric forest with a drawn
  N(0, 1) coefficient starting at 1 and its own size (dec-B307; run: created, `sd` 2, `k.scale` 2.97).
- `blocks()` ([`blocks`](../../R/model.R)) confines each tree of one forest to a group of predictors, the trees
  split evenly among the groups by default, the first `n %% G` groups taking one more; every tree keeps the
  forest's leaf sd. That is an additive model whose groups share the forest's total in proportion to their trees.
- Tree counts in print: a fit's per-forest lines ([`forestPriorLines`](../../R/forestPrint.R)) show no count;
  only the synopsis a fit without its call prints has them ([`fitTreeCounts`](../../R/generics.R)).
- Defaults by kind ([The rule](../design/forest-defaults-by-kind.md#the-rule)): the forest with no basis takes
  the fitting function's `n.trees` (75 at `bart()` and `dbarts()`), tree prior, `interactions` and `blocks`; a
  forest with a basis takes 50 trees, base 0.25, power 3. [`forestParams`](../../R/model.R) writes each forest's
  count; forest 1's goes over the control's slot, where the bridge reads it.
- Starts. A drawn coefficient on a forest with no basis starts by position: 1 at the first, 0 at the second
  ([`AmplitudeState`](../../src/bartcore/combiner.hpp)'s seed layout; TODO drawn-coefficient-start-by-kind).
- k and `k.scale` on several forests: k is 1 and `k.scale` the forest's own prior sd before its coefficient
  (forest-multiplier's call (2), held for this item).

## Design

### Statistics: one total split, against each at the full total

At a point x a forest's fit is the sum of one leaf from each tree, the leaves independent N(0, s^2 / m) over m
trees, so f(x) ~ N(0, s^2) exactly, whatever the trees. BART's k / sqrt(m) rule sets s = `k.scale` / k for the
forest and splits it evenly over its trees. With K forests with no basis:

- Each at the full total, as K separate fits would be: each component has sd T, the sum sqrt(K) T. The sum
  widens with every forest added, 1.41 T at two and 1.73 T at three.
- One total split (the ruling): forest j's component has variance s_j T^2, its shares s_j summing to 1, so the
  sum has sd T exactly, at every x. A part per tree, s_j = m_j / M over the group's M trees, gives every leaf
  N(0, T^2 / M): one forest of M trees with its trees in groups, the `blocks()` model as main builds it. An equal
  part per forest, s_j = 1 / K, gives forest j's leaves N(0, T^2 / (K m_j)). The two agree where the counts are
  equal: each component has sd T / sqrt(K), 0.71 T at two, 0.58 T at three, 0.5 T at four. Ruled: an equal part
  per forest (dec-B505), and `blocks()` changed to give each of its G groups 1 / G (dec-B506), so the two
  spellings stay one model at every count.
- A drawn total. One coefficient a ~ Cauchy(0, c) on the sum, f = a (F_1 + ... + F_K), each F_j at U sqrt(s_j):
  the sum has the one forest's prior exactly. Each forest its own a_j ~ Cauchy(0, c_K) instead: given the F's, the
  sum at x is Cauchy with scale c_K (|F_1| + ... + |F_K|), which no c_K makes the one forest's: at K = 3,
  c_K = c / sqrt(K), the variance split, leaves it about twice as wide, and c / K matches its upper tail but not
  its median (the critique's check at 4e5 draws: 2.1 times, and 1.24 times at the median).

The probe (scratch/fst-plan/prior.R; pure R, trees from the CGM prior, base 0.95, power 2, 100 points in
[0, 1]^3, 400 draws, three one-variable forests of 75 trees against one forest on all three, T = 1):

| model | sd f(x) | sd of one component | median spread over x | probit tail, T = 1.5 |
|---|---|---|---|---|
| one forest | 1.04 | - | 0.58 | 0.138 |
| three, each at the full total | 1.78 | 1.04 | 1.02 | 0.383 |
| three, equal split | 1.01 | 0.59 | 0.58 | 0.126 |
| three of 25 trees, leaf sd of one forest (`blocks()`) | 1.03 | 0.58 | 0.59 | 0.129 |

The probit tail is the probability that Phi(f(x)) is below 0.01 or above 0.99 at a fixed total of 1.5, the k = 2
binary size. With the drawn coefficient at c = 2 (median |f(x)|, 90th percentile, median spread over x): one
forest 1.26, 9.6, 1.21; one shared coefficient 1.22, 8.9, 1.21; each its own at c / sqrt(3) 2.41, 15.1, 2.03; at
c / 3 1.30, 9.1, 1.09; at c 4.60, 28.7, 4.06.

Practice. BART's own rule is the even split over trees. Fitters of several forests or terms give each its own
prior, whatever the count: flexBART 0.5 for every forest of a model of several, bcf and stochtree each forest
its own, mgcv and rstanarm each term its own (dec-B293, dec-B291). The shared total is this package's rule, from
the ruling, and it is the trees' rule carried to forests.

### The questions

D1. Which total. RULED, dec-B509: the fitting function's leaf prior in every model, whatever other forests it has:
`bart()`'s `k`, `leaf.prior = normal(k = )`, a drawn k under `chi()` or `normal(sd = )` at either function, with
the coefficient held at 1. The defaults are `bart()`'s: k 2 under gaussian, the drawn k under chi(1.5, 2) under
probit and logistic. Under gaussian at the defaults that is the held size main already gives a held forest with
no basis (dec-B290, a quarter of the range), so a held gaussian fit keeps its value; a drawn one does not. Under
probit and logistic it supersedes dec-B293's held size of 1 unit for a forest stating no sd, which dec-B509 does
not name (its "whatever other forests the model has" and its "a drawn k under chi()"); the records should mark
dec-B293 revised in part at slice F's landing. Built in two steps: slices 1 and 2 for models with two or more
forests with no basis that state no sd, with or without a forest with a basis (new models, neutral), and
[Slice F: the default flip](#slice-f-the-default-flip) for the lone such forest beside any other, main's model among
them.

D2. How the engine tells them apart. The sharers are the forests of kind none that state no sd, read from what
every forest already states ([`ForestStatement`](../../src/bartcore/combiner.hpp), its kind); whether the model
has a forest with a basis is read the same way. No new field, no flag from R, one engine function that answers
both. They stay separate engine forests, so every per-forest surface (labels, fits, counts, trees, weights, move
settings, tree priors) works as it does. MY CALL. Alternatives: an explicit group handed by R, which the ruling's
"no switch" leaves nothing to set; and one engine forest whose trees are the forests (a generalized `blocks()`),
which would need every per-forest knob and reader made per block.

D3. A drawn total is one draw: a drawn k is one k, drawn from the pooled leaves of every sharer, in every model.
MY CALL: it is what makes the sum exactly one forest's prior, and "as the trees of a forest do" is one scale.
Alternative: each its own k, which matches the one forest only in scale (the table above) and has no total that
is shared. The shared drawn coefficient this plan had beside a forest with a basis is gone with dec-B509: no
sharer draws a coefficient by default. `amplitude = student(df)` (coefficient-prior-opt-in under dec-B509, the
opt-in drawn coefficient) is taken on a forest with no basis that states its own `sd`, which is outside the share
and draws its own coefficient at median sd, and on the one forest with no basis that states none where it shares
with no other; it is refused by name where two or more forests share the total ("amplitude = student() draws one
forest's coefficient; these forests share one total: state each forest's 'sd' to give it its own"). MY CALL, the
smallest correct rule: it builds no shared coefficient. Alternative: one coefficient for the group under the t,
every sharer stating the same df, the shared draw this plan had (engine ~60 lines, tests ~130, a bcf-exact mode),
built for an opt-in nobody has asked to share; accepting it later changes no call.

D4. A stated `sd` on one of them. It is that forest's own size, as `sd` is everywhere (dec-B253): exact, held,
or the median where it states `amplitude = student(df)`, with its own coefficient or fixed k. It leaves the share, and the forests that state none
share the total among themselves. MY CALL. Alternatives: the stated number as the total the forest takes 1 /
sqrt(K) of, which no other `sd` means; and keeping K at every forest with no basis, so stating one forest's size
shrinks the others' total in silence.

D5. Probit and logistic. In every model the sharers take `bart()`'s drawn k under chi(1.5, 2), one k for the
group (dec-B509; dec-B420 found a fixed k not good enough for binary fits). The probit rescaling step
([`Chain::drawForestRescale`](../../src/bartcore/chain.hpp)), which mixes that k, is extended to the group and
to every other forest the model has: it scales the latents, every sharer's leaves and k by one factor, the leaf
count M the group's, and reads the current fit of every forest that is not a sharer (a forest with a basis,
drawn or held, and a forest stating sd) as part of the offset. The step's conditional already carries an offset
the group leaves unscaled, the term Q = sum o_i (z_i - f_i): with o_i the user's offset plus the non-sharers'
fit and f_i the sharers', (alpha z - o - alpha f)^2 expands as it does today, the truncation of z is invariant for
alpha > 0, and the Jacobian and the leaf and k priors are the sharers' alone. So the step is exact in every
model whose sharers draw k, with no change to its density, and the earlier decline where a forest states sd goes.
What it does not do is move a forest with a basis: near separation, where the fit's whole scale is pinned by the
latents, the part a forest with a basis carries does not move with them. MY CALL. Alternatives: (a) scale every
forest's contribution too, a drawn coefficient block a by alpha (adding q to the power and |a|^2 / (2V) to the
e^(2v) coefficient) and a held forest's leaves (adding their count and their standardized squares), exact and a
fuller move of the scale, ~60 more engine lines through the combiner's coefficient blocks, unmeasured; (b)
decline beside a forest with a basis, which leaves the default probit bcf model's k mixing as a lone forest's did
before the step (a quarter of the SBC census datasets under 50 effective draws of 99 at thin 1000, against none
with it; [probit-k-scale-move.md](../design/probit-k-scale-move.md)). (a) is the step to add if a binary causal
model's k is found to mix slowly. The logistic collapsed scale move (dec-B404), when built, covers the group
the same way. A forest with no basis stating `amplitude = student(df)` has no k (D3, coefficient-prior-opt-in),
so a model whose one forest with no basis draws a coefficient has no step, as today.

D6. The tree count. Ruled by dec-B503: the fitting function's `n.trees` is the group's, split as evenly as
`blocks()` splits it, and a forest that states `n.trees` keeps its own. A model of forests with no basis costs
one forest. MY CALLS under the ruling, [Calls made](#calls-made) 11 to 14: a stated count stands outside the
total; the forests stating none split all of `n.trees` among themselves; `bart(n.trees = )` beside a stated
count is taken, refused only where every forest with no basis states one; print shows each forest's count and
the total. So `forest(x1, n.trees = 10) + forest(x2) + forest(x3)` at the default is 10, 38 and 37 trees.

D7. The tree prior, `interactions` and `blocks`. Each forest with no basis takes the fitting function's tree
prior where it states none, as the one does today. `interactions` and `blocks` given to the fitting function
beside two or more forests with no basis are refused by name, telling the caller to state them on a forest, as
they are beside none ([`refuseStatedWithNoPlainForest`](../../R/model.R)): a partition of one forest's
predictors cannot be each forest's. MY CALL. Alternative: give each forest the fitting function's, which works
for `interactions` and never for `blocks`.

D8. k and `k.scale` (forest-multiplier's call (2)). The rule that covers every forest: `k.scale` / k is the
forest's own prior sd before its coefficient, in the response's units. A forest whose size a calibration map
sets, a forest with a basis or one stating sd, reports k 1 and `k.scale` that sd, as built. A sharer, in every
model (dec-B509), reports the group's k, one value reported on each, and `k.scale` its share times k, so `k.scale`
/ k is its share, and `getK` and `extract(type = "k")` give the drawn k on each. `bart(k = )`, `leaf.prior =
normal(k = )` and a `chi()` beside a forest with a basis are the group's, no longer refused. MY CALL. Alternative:
k 1 everywhere, which would hide the drawn k.

The item after this one, forest-amplitude-normal (dec-B542), changes the first sentence's scope. "A forest with a
basis or one stating sd reports k 1 and `k.scale` that sd" is true for slice F; after that item it holds only where
the coefficient is held: a drawn forest of any kind reports k 1 and `k.scale` the response's unit, and `$sd` keeps
the sd. So slice F's new pins read `$sd`, and pin a drawn `k.scale` only on a forest with no basis stating
`student(df)`.

D9. The constant basis column (dec-B307's revisit). RULED, dec-B528, [Open calls](#open-calls) 2: refused at creation,
sd stated or not, naming the plain spelling, in slice 3. Constant is what
[`basisColumnScale`](../../src/bartcore/combiner.hpp) finds no scale in, over the same rows: by its values, or to
rounding, over the rows in the likelihood, and over every row where too few are in. It is one refusal at the one
function every fitting function's door goes through ([`resolveSamplerSpec`](../../R/spec.R)): a formula's term, a
`forests` list, a data object's `bases`. The spelling it names is the same model: a stated sd is per unit of the
column, so it offers `sd` as the column's absolute value times the sd stated (`basis = rep(3, n), sd = 2` is
`forest(x1, sd = 6)`), and it keeps `amplitude = fixed()` or `student()` as written. A later swap to a constant column
stays taken (dec-B346), and a re-creation reads what the swap left, so a copy and a reload are not creations. `new()`
over a live sampler's control and model and a data object whose column was edited to a constant is taken too, where
an sd is stated and, with none, where the data's `basis.scale` record stands; the engine refuses only a column with
neither (run). It cannot be told from `new()` over a swapped sampler's parts, and the route is nearly closed, a plain
control refused and a fresh `dbartsData` failing on `n.cuts`, so no refusal is added there. A factor of one level and
a logical column at one value stay refused by their own messages at every door, sd stated or not; a column of zeros
by its own. A factor, character or logical basis left with one level in the likelihood by the weights is refused the
same way, naming no spelling, a level with no row leaving its coefficient at its prior (dec-B543, at creation only,
as dec-B528).

D10. The rest of the surface, each MY CALL:

- Holds. One rule for `amplitude` on a forest with no basis, in every model, after slice F: `fixed()` restates
  the default, held, and is taken, on a sharer, on a forest stating sd and on a lone forest (the last is
  [Open calls](#open-calls) 4, ruled by dec-B513: refused today, "this model has one forest, which has no
  coefficient to hold"; taken in slice F, with no message); `student(df)` is the opt-in drawn coefficient, taken as D3 says and refused
  on a lone forest, whose single-forest chain has no coefficient (the refusal reworded: "a model of one forest has
  no coefficient to draw: for a drawn size give leaf.prior = normal(k = chi())"). Sharers holding and not holding
  cannot mix: none draws. Slices 1 and 2 build this for models with two or more forests with no basis, all new;
  a lone forest with no basis beside any other forest keeps main's drawn default and main's `amplitude` rules
  until slice F, and a lone forest keeps main's refusal of `fixed()` until slice F takes it (dec-B513).
- Starts. Where two or more forests have no basis, every drawn coefficient on such a forest (one stating
  `student()`) starts at 1, so none starts out of the model; a model with one keeps its start by position. No
  model main accepts moves; TODO drawn-coefficient-start-by-kind is noted.
- A restated sd. dec-B525 replaces this paragraph's refusal: a sharer may state its own size live, leaving the group
  with the others' shares unchanged. The refusal goes to forest-door's slice 3, a sharer leaves the group, which
  removes it; this branch keeps it, among two or more sharers, until then, and a lone sharer takes the sd already
  ([Calls made](#calls-made) 24). `$setLeafPrior(forests = )` restating a sharer's sd is refused by name: membership is fixed at
  creation, and a sharer leaving would move the others' share. A forest that stated sd restates it as today, and
  `forest()` written back stays the no-op it is (dec-B508: `forest()` states a size as `sd` alone). On a live
  sampler with sharers, in every model, the group's leaf prior (its k, drawn or held, or a named sd) is a
  parameter (dec-B254) and is set through `$setModel`'s leaf prior, built in setmodel-parameters-not-structure's
  slice 4 after this item lands; until then `leaf.prior =` on such a sampler stays refused, as on every sampler of
  several forests.
- Families and compositions. Gaussian, probit and logistic, and every refusal a model of several forests has
  today stands, test predictors and the variance forest among them. Each is an addition later.
- The formula. Plain terms beside a `forest()` with no basis stay refused: mgcv reads `y ~ x1 + s(x2)` with x1 a
  linear term, so `y ~ x1 + forest(x2)` has two readings. The message now says to write each as a `forest()`.
  A `forest()` with no predictors keeps its reading, every predictor of the fit.
- Reading. `getForestAmplitudes` on a sharer gives its held 1; `print` and `summary` mark a sharer's line as a
  share of the total and show each forest's tree count and, after the forests, the total and the group's k; the
  verbose creation summary prints the group once. A forest stating `student()` reads as coefficient-prior-opt-in
  says.

D11. How the total divides among sharers of unequal tree counts. Ruled by dec-B505: an equal part per forest,
s_j = 1 / K, whatever its count, its leaf scale dividing that part among its own trees, so each sharer has its own
leaf scale. dec-B506 changes `blocks()` to match ([Slice B](#slice-b)), so
the model is `blocks()` with `trees.per.group` at every count.

### R surface

1. [`plainForest`](../../R/model.R) returns every position with no basis and refuses nothing; its callers take a
   vector. The formula reader's count refusal goes; its plain-terms refusal is reworded (D10).
2. [`forestParams`](../../R/model.R), by dec-B503 and D6: the forests with no basis that state no `n.trees`
   split the fitting function's whole count evenly among themselves, the first `n %% K` taking one more, as
   `blocks()` splits; a forest's own count governs it and is outside the split. The fitting function's
   `n.trees` beside a stated count is taken, and refused by name only where every forest with no basis states
   one, as on one forest today. A split leaving a forest no tree is refused, as `blocks()` refuses it. Where
   forest 1's count is not the fitting function's, the
   forests' record keeps the count the next fit inherits, as it does where forest 1 has a basis
   ([Where the numbers live](../design/forest-defaults-by-kind.md#where-the-numbers-live)). Each sharer takes the
   fitting function's base and power where it states none.
3. `interactions` and `blocks` given to the fitting function beside two or more forests with no basis are refused
   (D7), in [`refuseStatedWithNoPlainForest`](../../R/model.R)'s style.
4. [`resolveSamplerSpec`](../../R/spec.R): where a model has sharers (two or more in slices 1 and 2; any in
   slice F), the k hyperprior and a k other than 2 are no longer refused, and a fitting function's `leaf.prior =
   normal(sd = )` states the group's total instead of being refused with its redirect;
   [`resolveLeafHyperprior`](../../R/model.R) gives such a model the binary default chi(1.5, 2) instead of a
   fixed 2. Where no forest shares (every forest with no basis states sd, or none exists), a leaf prior given to
   the fitting function has no forest to belong to: the default restated (`normal(2)`, `k = 2`) is taken
   silently, as main takes it, and anything else is refused by name. Everything else it refuses for several
   forests stands.
5. Holds (D10): `fixed()` on a forest with no basis taken as the default restated; `student()` on a sharer refused
   where two or more share (D3).
6. `$setLeafPrior(forests = )` refuses a sharer's sd ([`resolveForestSpreads`](../../R/dbarts.R)). dec-B525 replaces
   this: the refusal goes to forest-door's slice 3, a sharer leaves the group, and this branch keeps it until then.
7. [`reportLeafPrior`](../../R/dbarts.R) and the k readers follow D8; `prior.sd.of` gains ", shared" on a sharer;
   [`forestPriorLines`](../../R/forestPrint.R) marks the share and gains each forest's tree count, and the fit's
   print the total, on every model of several forests. A fit of sharers with a drawn k stores it as a one-forest
   fit stores its k, wherever the sharers stand. `extract(type = "k")` reads it per forest as
   `type = "leaf.prior.sd"` does: a sharer the draws, a forest outside the group its own number, a name or index that
   is no forest refused, `forest = NULL` the named list. `summary` tabulates the drawn k, and `vars = "leaf.prior.sd"`
   the row `extract` has; a named drawn total, `normal(sd = invchi())`, reads `leaf.prior.sd`, as one forest does.

### Bridge

1. Where a model has sharers, the K-forest creation hands the chain the leaf prior the single-forest path reads
   (k, whether it is drawn and its chi prior, a named sd), which today it drops behind R's refusal.
2. `bartcore_getLeafPrior` reports each forest's k and the share its map entry carries; the verbose summary
   ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)) prints the group.
3. Nothing in [`applyForestBases`](../../src/R_interface_bartcore.cpp) or the data record changes: a sharer has
   no basis and no scale.

### Engine

1. Sharing. One function reads the specs: the sharers (kind none, no sd stated, no df stated) and their tree
   counts. [`forestMapEntry`](../../src/bartcore/combiner.hpp) takes the sharer's share s_j = 1 / K (dec-B505)
   and multiplies a sharer's node scale by sqrt(s_j). A lone sharer takes no branch: the expression as it stands,
   operand for operand. Sharers are held (dec-B509): no coefficient is drawn for one.
2. (Gone with dec-B509: the shared coefficient. [`drawAmplitudes`](../../src/bartcore/combiner.hpp) keeps its
   per-forest loop.)
3. Starts (D10) at construction, after the layout: every drawn kind-none block at 1 where two or more forests
   are of kind none.
4. Sharers' leaf prior. The K-forest [`Chain`](../../src/bartcore/chain.hpp) constructor takes the options' k, its
   draw and chi prior, and the single-forest chain's node scale (the expression
   [`singleForestDefaultSize`](../../src/bartcore/chain.hpp) and [`translateNamedSd`](../../src/bartcore/chain.hpp)
   rest on, called and not copied) times sqrt(s_j), in place of the map's k = 1, for every sharer, whatever other
   forests the model has; a forest stating sd and a forest with a basis stay on the map. One exception, MY CALL:
   a lone sharer beside a forest with a basis whose k is held at the default 2 keeps dec-B290's held map entry,
   operand for operand, so every held gaussian fit main accepts keeps its draws bit for bit (slice F's gates read
   it); its reader still reports k 2 by D8. Alternative: the k path for it too, the same prior with draws that
   differ in rounding, which would move bcf-equivalence's held scenario and every `fixed()` pin for no gain.
5. The pooled k. The k draw after each sweep draws one k for the group from the sum over sharers of each forest's
   sum of squared leaves over its own leaf scale squared, and the sum of their leaf counts, on unit scale
   ([`ChiKHyperprior`](../../src/bartcore/chain.hpp)), and writes it to every sharer. Each sharer has its own leaf
   scale (dec-B505), so each forest's sum enters standardized by its own scale, never by one shared scale; with
   equal counts the scales are equal and this is the single-forest draw over M trees. Any other forest draws as
   now.
6. The probit rescaling step (D5): [`Chain::forestRescaleApplies`](../../src/bartcore/chain.hpp) admits a K-forest
   chain whose sharers draw k on constant leaves, today refused by its `combiner_` test, whatever other forests it
   has; the step reads R and Q with f the sharers' summed fit and the offset the user's plus every other forest's
   current contribution (the combined fit less the sharers'), and scales each sharer's leaves, `totalFits` and
   `totalFitsScale` under the existing bound and rule. The combiner's combined-fit cache (`glue_.combined`) needs no
   refill: it is read only inside [`AmplitudeForestCombiner::combinedFits`](../../src/bartcore/combiner.hpp), which
   fills it, so a refill planned here was dropped, the draws bitwise without it. M, the group's leaf count, cancels
   between the Jacobian and the leaf prior, so the step reads none. No other forest's leaves, coefficient or fits are touched. The step's conditional reads no leaf scale: (mu, k) -> (alpha mu,
   k / alpha) leaves every standardized leaf mu k / s unchanged whatever s is, so the leaf prior's Jacobian is
   alpha^-M with one scale per forest or per group as with one in all. Per-forest scales need no weighting there;
   a tests/cpp pin holds the step's density at unequal scales to a direct evaluation.
7. State. The coefficient layout and [`glueIsValid`](../../src/bartcore/combiner.hpp) do not change: each
   sharer's block holds its held 1. `stateFormatVersion` does not move.
8. [`ForestCalibration`](../../src/bartcore/chain.hpp) carries whether a forest shares, its share s_j and, for a
   sharer, k in force.

No facade virtual changes. The shipped header ([`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)) has no
create entry and does not change; stan4bart builds one forest.

### Slice B

`blocks()` gives each group an equal part. dec-B506: a forest partitioned by `blocks()` gives each of its G groups 1 / G of the forest's prior variance,
whatever `trees.per.group` says, each group's leaf scale dividing its part among its own trees; with equal counts
nothing changes. Today every tree has the forest's leaf scale s / sqrt(m) ([`installBlockMasks`](../../src/bartcore/chain.hpp)
assigns only columns), so group g's part is m_g / m. After it, a tree of group g has leaf scale
s / sqrt(G m_g), the forest's times f_g = sqrt(m / (G m_g)).

Engine, on every chain that takes `blocks()`: the single-forest chain, a forest of a K-forest model and the
category forests of a multinomial model.

1. `installBlockMasks` computes f_g from the group counts and stores it on the forest only where the counts are
   unequal; a forest with equal groups or no groups stores none and takes today's code, operand for operand.
2. Each tree's update: before tree t's moves and parameter draw, the leaf's scale is the forest's times
   f_{g(t)} (from `blockOfTree`), restored after; every read inside a tree's update (the leaf model's draws and
   integrated likelihood, the rule scan's prior variance in [scan.hpp](../../src/bartcore/scan.hpp)) sees the
   group's scale with no change of its own. The same guard is required, not optional, on the grow-from-root
   sweep and the prior leaf draw: both draw leaves the factored k statistic then reads, and a prior draw must
   be the group's prior. The monotone joint-prior acceptance in `sampleTreesFromPrior` stays unfactored: it
   orders leaves within one tree, and the order does not depend on the scale.
3. The reads outside a tree's update take the tree's factor: the k statistic accumulates (mu / f_g)^2, so the k
   draw on the forest's scale is exact; [`Chain::drawLevelShift`](../../src/bartcore/chain.hpp)'s per-leaf
   precision; and the multinomial grand-level shift in
   [`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp), which reads one s per category
   forest. The probit rescaling step reads no leaf scale (engine item 6).
4. Nothing in R or the bridge changes: the group counts already reach the engine. No facade virtual, no state
   block; a state, a copy and a warm start place trees by slot and the slot's group gives the scale.

Help: man/blocks.Rd says each group takes an equal part of the forest's prior size, `trees.per.group` setting
how finely each divides it. Design: the design note of this item carries the rule;
[P4 landing (2026-07-21)](../design/interaction-constraints.md#p4-landing-2026-07-21) gains a pointer. NEWS: none
of its own, `blocks()` being new in 1.0-0; its existing entry states the rule.

What moves. A `blocks()` fit whose groups hold unequal counts draws differently from main; each group's prior sd
as a fraction of the forest's moves from sqrt(m_g / m) to sqrt(1 / G):

| groups (trees) | today | after |
|---|---|---|
| 2 of 75 (38, 37) | 0.712, 0.702 | 0.707 each |
| 4 of 75 (19, 19, 19, 18) | 0.503, 0.490 | 0.5 each |
| 2 of 20 stated (5, 15) | 0.5, 0.866 | 0.707 each |

At the default splits the posterior change is well inside Monte Carlo error; at stated uneven counts it is the
prior the ruling asks for. Equal splits (75 over 3 or 5) and every fit without `blocks()` are bitwise as main.

Baselines and snapshots, by name, with the verdict expected: no scenario of equivalence-32418b6e,
bcf-equivalence-20eed7ed or multinomial-equivalence-80b1c8d4 states `blocks()`, so the trio stays 55, 19 and 11
identical under `compare --bitwise --strict-coverage`, no re-record and no z mode needed; none of the four
seeded-drift snapshot files (test-reproducibility-binaryResponse.R,
test-reproducibility-continuousResponse-singleThreaded.R, test-reproducibility-continuousResponse-multithreaded.R,
test-reproducibility-xbart.R) uses `blocks()`, so all pass unchanged; bench-sampler-127f04ee is not re-recorded.
A difference in any of them is a stop: it would mean the factor reached a forest with equal groups.

Tinytest pins: none in test-blocks.R or the other files that use `blocks()` (test-bcf-creation.R,
test-constructor-vocabulary.R, test-forest-defaults.R, test-setstate-force-update.R, test-single-forest-vars.R,
test-warm-start-forms.R, test-warm-start-kinds.R) reads a draw's value: they pin confinement, counts, refusals,
warm-start placement and run-to-run identity, all of which hold. Unequal-count fits there (5 and 15, 15 and 5,
two groups of 75's default) draw differently and keep passing; the implementer re-runs each file and reports any
pin that moves as a finding.

Tests and gates, slice B:

- tests/cpp: the factor stored only for unequal counts, every count compared ({7, 7, 6} stores three factors,
  {5, 5, 5, 5} none); each tree's leaf scale s / sqrt(G m_g) during its update and the forest's after, in a
  sweep, a grow-from-root sweep and a prior leaf draw; the k statistic against a hand sum on the constant leaf
  ((mu / f_g)^2), the linear leaf ((param / f_g)^2 per coordinate) and the GP leaf (each leaf's sum of squares
  over f_g^2), each apart from the unfactored sum; the level shift's precision with factors; the multinomial
  shift with factors; equal counts bitwise the old sweep.
- Exact posterior, a new benchmarks/R/blocks-exact.R joining exact-gates.yaml: gaussian, sigma fixed, one forest
  of three trees in two groups, x1 (two trees) and x2 (one tree), each a 0/1 column with one cut; the 8 tree
  structures (each tree a root or its one split) enumerated, leaves integrated in closed form under the group
  variances T^2 / 4 and T^2 / 2; the sampler's posterior split probability per tree and E f(cell) held to the
  enumeration. A second arm with k drawn under chi(1.5, 2), an outer quadrature over k. The per-tree prior
  (T^2 / 3 on every leaf) fails it, which the gate is run against once to prove it discriminates.
- tinytest, in test-blocks.R: an unequal fit (5 and 15) created, run, copied and reloaded after `storeState`,
  the copy continuing `identical()` to the reload and within 1e-12 of the live sampler (re-creation re-derives
  the cached fits, as on main); a warm start between two unequal samplers.
- The A/B against an e2f37f6d library, every channel `identical()`: no `blocks()`; `blocks()` with equal groups
  (75 over 3) on the single-forest chain, a forest of a bcf model and a multinomial model; a state, copy, reload
  and warm start of each. The unequal arm is excluded and covered by the exact gate.
- Gates: tests/cpp plain and under ASan/UBSan, the R-loaded sanitizer run over test-blocks.R, the full suite,
  the trio and snapshots above, `R CMD check --as-cran`, `bench-sampler.R` compare (maintainer-run): the sweep
  gains one branch per tree.
- Mutants, each to fail a named expectation: the factor sqrt(G m_g / m) inverted; the factor stored for equal
  counts (fails its tests/cpp pin; a factor of 1 is exact, so the A/B cannot see it); the equal-counts check
  reading only the first two groups (fails the {7, 7, 6} pin); the linear or GP leaf's k statistic unfactored
  (fails its hand-sum pin); the scale not restored after a tree (fails the next tree's tests/cpp pin); the k
  statistic unfactored (fails the k arm); `drawLevelShift` or the multinomial shift reading the forest's scale
  (fails its tests/cpp pin).

### Slice F: the default flip

dec-B509 for every model of two or more forests, by its "whatever other forests the model has": a lone forest with no
basis that states no sd, beside anything, takes the fitting function's leaf prior, held; a forest with no basis
stating sd and no `amplitude` is held at it. So it moves main's model, a forest beside a forest with a basis, and a
model slices 1 and 2 newly take with main's drawn defaults, `forest(x1) + forest(x2, sd = 2)`, refused on main, whose
unstated forest shares with no other and whose stated forest draws a coefficient at median 2. The drawn coefficient is
opt-in, `amplitude = student(df)`. Built after a rebase onto a main that carries coefficient-prior-opt-in, and merged only after
bartCause's `bcf()` states its prior ([Interaction with the other pending plans](#interaction-with-the-other-pending-plans)).
dec-B539's forest-amplitude-normal, the drawn `amplitude = normal()`, builds after this slice, taken where
`student(df)` is. dec-B540's door, `$setLeafPrior` on the one forest that draws a coefficient, reads the scale creation
set, which engine item 3's median is taken from and the record carries.

Engine:

1. [`forestMapEntry`](../../src/bartcore/combiner.hpp): a forest of kind none draws a coefficient only where it
   states a df (coefficient-prior-opt-in's `amplitudeDf`); otherwise it is held, whatever `updateAmplitude` says.
   The hold keys on `isnan(df)`, the statement, not on a local that excludes +inf, and `dfStated = !isnan(df)` stays
   in `forestMapEntry`: forest-amplitude-normal (dec-B539) states a normal coefficient as an internal df of +inf on
   `ForestStatement::amplitudeDf`, and a hold reading anything else would hold a `normal()` forest silently once it
   lands.
2. The lone sharer takes engine item 4's path (the held k = 2 default on dec-B290's entry, its exception) and
   item 6's rescaling step.
3. A forest of kind none stating a df and no sd, `!isnan(df)` by the same reading as item 1 so that dec-B539's +inf
   takes it too: its median, in place of
   [`defaultForestSd`](../../src/bartcore/combiner.hpp)'s 2 or 1, is the prior median of the size it would
   otherwise have: `k.scale` / k at a held k, `k.scale` over the chi prior's median k at a drawn one (the size is
   monotone in k; scale times the square root of the chi-square median at its df, through the support library),
   a named sd as given, and under `normal(sd = invchi(df, s))` that prior's median, s over the chi's, so that it and
   `normal(k = chi(df, k.scale / s))`, one prior, give one median. At `bart()`'s binary default that is 3 / 1.906 = 1.574 under probit and 5.44 / 1.906 =
   2.85 under logistic, in the link's units; the forest then has no k (coefficient-prior-opt-in).
4. [`expandForestSpecs`](../../src/bartcore/combiner.hpp): the K = 2 mu/tau shape states df 1 on mu where
   `updateA`, so the tests/cpp bcf fixtures keep bcf's prior and their draws, the shape staying a pure re-spelling
   of bcf.

R: item 4 of [R surface](#r-surface) for a lone sharer beside any other forest (k, `chi()`,
`normal(sd = )` taken; the binary default drawn); item 5's holds, with `fixed()` on a lone forest taken as
the default restated, with no message (dec-B513, [Open calls](#open-calls) 4), `student()` there still refused; the record unchanged, the hold of a forest with no basis being the engine's
default and `heldAmplitude` the statement.

What moves, a forest with no basis in a model of two or more forests. The column main is its reading beside a forest
with a basis; `forest(x1) + forest(x2, sd = 2)`, refused on main, reads the same after slices 1 and 2 and moves the same
way, `forest(x1)` by the first rows and `forest(x2, sd = 2)` by the `sd` row:

| statement | main | after |
|---|---|---|
| none, gaussian | drawn, half-Cauchy median 2U | held at k 2, a quarter of the range |
| none, probit or logistic | drawn, median U | held, k drawn under chi(1.5, 2) |
| `fixed()`, gaussian | held, a quarter of the range | the same, bit for bit |
| `fixed()`, probit or logistic | held at U | held, k drawn under chi(1.5, 2) |
| `sd` alone | drawn, median sd | held at sd |
| `sd` and `fixed()` | held at sd | the same, bit for bit |
| `student(df)`, no sd | median 2U or U | median of the size the forest would hold |
| `bart(k = )`, `chi()`, `normal(sd = )` | refused | the forest's leaf prior |

Baselines and snapshots, by name, with the verdict expected:

- equivalence-32418b6e: bart2twoforest (`y ~ x1 + x2 + forest(x1 + x2, basis = ~zf)`, gaussian, default) moves.
  `compare --bitwise --strict-coverage` with it excluded: 54 identical; then re-recorded, the one baseline
  dec-B509 named.
- bcf-equivalence-20eed7ed: its `bcfForests` helper restates bcf's prior where it draws, `amplitude = student(1)`
  and `sd` 2 sd(y) under gaussian, 1 under probit, pi / sqrt(3) under logistic. Expected: latent_probit identical
  (its unit is exactly 1, run), the held gaussian scenario identical (engine item 4's exception), the logistic
  scenarios identical (R's pi / sqrt(3) is the engine's unit bit for bit, its draws `identical()`, run), and the
  gaussian drawn scenarios changed in rounding only, R's sd(y) and the engine's unit differing in the last place:
  `bcf-equivalence.R compare <baseline> --cross-host` against 20eed7ed, tier 1 with varcount identical, then
  re-recorded. z mode was the first gate named here and fails main against itself under a changed seed (17 of 19
  scenarios |z| 2.2 to 15.3, run). MY CALL. Alternative: state sd as twice the engine's unit read from a
  probe sampler (a factor forest's default sd is exactly one unit), bitwise with no re-record, at a probe per
  scenario in a gate script. As built: 16 of 19 bitwise, glue_toggle, set_treatment and per_observation_partial
  moved by rounding, tier 1's largest deviation ratio 9.3e-07, varcount identical on all 19.
- multinomial-equivalence-80b1c8d4: 11 identical. The four snapshot files use no forest with a basis: unchanged.
- Exact gates: bcf-exact.R, bcf-exact-weak.R and bcf-latent-exact.R add `amplitude = student(1)` to their drawn
  modes, which already state sd: the same chains, passing unchanged; their held modes and bcf-exact-restricted.R
  and forest-prior-exact.R state sd and `fixed()` and do not move, save forest-prior-exact.R's arm of a forest with
  no basis under `student(1)` at no sd, whose oracle's median moves from 2 units to engine item 3's size. probit-k-scale-exact.R's arm from slice 1
  covers the step beside a forest with a basis, and gains a lone-factor arm: one frozen forest with no basis under
  one drawn k beside a held factor, the lone sharer's own path. sbc.R's bcf arm adds `student(1)` where it draws (its next full
  run, not CI). bench-sampler-127f04ee is not re-recorded.
- tinytest: about twelve files pin the drawn default beside a forest with a basis or the k refusal
  (test-bcf-creation.R, test-bcf.R, test-bcf-family.R, test-calibration-midchain.R, test-forest-arguments.R,
  test-forest-basis-r5.R, test-forest-defaults.R, test-forest-held.R, test-forest-leaf-prior-reader.R,
  test-forest-sd-unit.R, test-formula-terms.R, test-multiforest-leaf-prior-writer.R): a pin of the drawn
  coefficient states `student(1)` and its sd; a pin of the default is re-derived; each reported per file.

Tests and gates, slice F:

- tests/cpp: kind none held without a df and drawn with one; `expandForestSpecs`' df 1, the bcf fixtures bitwise;
  the lone sharer at a held k of 2 bitwise dec-B290's entry; the no-sd median at a held k, a drawn k (1.574 under
  probit), a named sd and a named sd under a drawn k; `forestsTakeLeafPrior` directly.
- tinytest: the default beside a forest with a basis and beside a forest with no basis stating sd, under gaussian,
  probit and logistic (the reader's k, `k.scale` and no `amplitude`; the coefficient 1 at every draw); `bart(k = 3)`,
  `chi()` and `normal(sd = )` taken; `fixed()` on a lone forest taken (dec-B513); `student(1)` with `sd = 1` under probit
  `identical()` to the pre-flip default; the two spellings of one drawn prior one median; the sd in force written on a
  lone forest the same prior, to rounding under a named sd.
- An A/B (scratch) against a library of main after coefficient-prior-opt-in: every fit of one forest, every model in
  which each forest with no basis states sd with `fixed()` or `student(1)`, and the held gaussian default identical;
  the moved rows of the table changed, as listed. The same rows beside a forest with no basis stating sd are held
  against this branch before the flip, main refusing that model: `forest(x1) + forest(x2, sd = 2)` changes, both
  forests moving, and `forest(x1, sd = 1, amplitude = fixed()) + forest(x2, sd = 2, amplitude = fixed())` does not.
- exact-gates.yaml's list in quick mode; the trio as above; the full suite; tests/cpp and the R-loaded sanitizer
  run (the chi median is new numerics); `R CMD check --as-cran`; lint, air, the codoc and freshness checks.
- Mutants: kind none drawn without a df (the default pins); the dec-B290 exception dropped (the held A/B);
  `expandForestSpecs` without df 1 (the tests/cpp bcf fixtures); the no-sd median left at 2U (its pin); the chi
  median at the wrong df (the 1.574 pin); the k refusal left in place (`bart(k = 3)`).

Help, design and records: man/forest.Rd's `amplitude` and sizes by kind (bcf's prognostic prior written
`amplitude = student(1), sd = 2 * sd(y)`); man/bart.Rd and man/dbarts.Rd on `k` and `leaf.prior` beside a forest
with a basis; NEWS edits 1.0-0's forest item, no entry of its own (nothing released moves);
[The rule](../design/forest-defaults-by-kind.md#the-rule) and [1. The rule](../design/forest-multiplier.md#1-the-rule)
amended; docs/design/bcf.md's dec-B18 record noted revised for a formula model; this item's design note gains the
flip and the check; dec-B293 marked revised in part; TODO no-basis-forest-prior removed at the merge.

## Fit and state

- No new record: the group is read from the forests' statements and kinds at every creation, so `copy()`, a
  reload and a dead-pointer `setState` rebuild it.
- A stored state holds a sharer's held 1 on its block; a state installs across samplers as today.
- A control taken from such a fit carries the forests' record as today, with the inheritable count of R item 2.

## Help, NEWS and docs

- man/forest.Rd: several forests with no basis, what they share and how; `sd` on one takes it out of the share;
  `amplitude` and the holds; the count split, that a stated count is outside it and sets how finely that
  forest divides its equal part, not its size. man/dbarts.Rd and man/bart.Rd: the `forests` and formula items,
  `n.trees` (the group's, taken beside a stated count) and `k` beside several forests with no basis; the
  fitting function's `leaf.prior = normal(sd = )`, the group's total, beside `forest(sd = )`, one forest's own
  size (dec-B508: `forest()` takes no `leaf.prior`). man/dbartsSampler-class.Rd: `getLeafPrior` and `getK`
  per D8, `setLeafPrior(forests = )`'s refusal. man/blocks.Rd: one sentence relating the two.
- NEWS, in 1.0-0's forest item: a model may hold several forests with no basis, sharing one total.
- A design note, docs/design/forests-shared-total.md, the rule and its statistics, landing with slice 2. Amended
  with it: [The rule](../design/forest-defaults-by-kind.md#the-rule) (at most one forest with no basis);
  [The calibration map, general in K](../design/multiplier-combiner.md#the-calibration-map-general-in-k)
  ("Adaptivity is capped at one forest" and the dispersion paragraph); [1. The rule](../design/forest-multiplier.md#1-the-rule)
  (defaults); TODO forests-shared-total removed, drawn-coefficient-start-by-kind and forest-multiplier's call (2)
  noted.

## Tests and gates

tests/cpp:

- The map: a sharer's node scale is the lone forest's times sqrt(1 / K), K = 2 and 3, equal and unequal counts,
  held and drawn, gaussian and probit; a forest stating sd beside sharers keeps its own; K = 1 bitwise the old
  expression.
- The pooled k: one k on every sharer; the draw's inputs against a hand sum, each forest standardized by its own
  leaf scale, at unequal counts; with slice B, a K-forest chain of sharers and the single-forest `blocks()` chain
  of the same counts give equal leaf scales tree for tree and equal k statistics on equal leaves.
- The rescaling step on a two-sharer probit chain, alone, beside a forest stating sd and beside a factor forest:
  admitted, M the group's, R and Q against a hand sum with the other forests' fit in the offset; after a step,
  each sharer's `totalFits` and the combiner's combined fits equal a re-derivation from the leaves, and no other
  forest's leaves, coefficient or fits move.
- Starts at 1 with two kind-none forests; by position with one.
- Fixtures that build two forests of kind none today are re-read: each states its kind or moves to the rule with
  its pins re-derived, reported.

tinytest, a new test-forests-shared-total.R, small n:

- Created at every door: a formula, a `forests` list, a data object's `bases`; gaussian, probit, logistic.
- The reader: each sharer's `sd` and `k.scale` the lone forest's times sqrt(s_j), K = 2 and 3, drawn and held; a
  stated forest's its own; D8's k on each.
- The parts, per forest (dec-B505): with no forest with a basis and k fixed, each sharer's `sd` is T / sqrt(K)
  and its leaf scale T / sqrt(K m_j) (to 1e-12), at K = 3 of 75 trees, at K = 2 of 75 trees (38 and 37: leaf
  variances T^2 / 76 and T^2 / 74, not 75's), and with a stated count (`forest(x1, n.trees = 10) + forest(x2) +
  forest(x3)`: each forest T / sqrt(3), its leaves over 10, 38 and 37 trees). A part per tree fails the last two.
  The equality with `blocks()` at those counts is a tests/cpp pin (each sharer's leaf scale equal to the matching
  group's after slice B) and the statistical equivalence below.
- Sharers held: `getForestAmplitudes` 1 on every sharer at every draw, beside a forest with a basis too.
- Counts (D6): 75 split 25, 25, 25 and 38, 37; a stated count beside unstated ones, 10, 38, 37, with and
  without `bart(n.trees = )`, which is then the unstated forests' count; `bart(n.trees = )` refused where every
  forest with no basis states one, and on one forest as today; `print`'s per-forest counts and total.
- Holds and opt-ins (D3, D10): `fixed()` on a sharer taken, the same draws as without it; `student(1)` on a
  sharer among two refused; `student(1)` with `sd` on a forest with no basis beside sharers taken, its own
  coefficient drawn; a leaf prior at the fitting function where no forest shares, `normal(2)` taken and
  `normal(k = 3)` refused.
- Refusals, each by name: `interactions` and `blocks` on the fitting function; a sharer's restated sd; plain terms beside a `forest()` with no basis (new text); a family or
  composition several forests refuse.
- Copy, `saveRDS` and `readRDS` after `storeState`: the next draws `identical()` to the live sampler's.
- Pins updated: the two refusals named in [What exists today](#what-exists-today) become creations.

Mutants, each to fail a named expectation: the scale times s_j rather than sqrt(s_j); the share m_j / M where
1 / K is ruled; the share counting forests that state sd; the pooled k standardized by the first sharer's scale
for all; the pooled k from the first forest's leaves; the rescaling step's leaf write stopping at the first sharer,
the step reading no M for the M of one forest planned here to break; the rescaling step scaling a non-sharer's fit, or leaving it out of the
offset (each fails the probit arm beside a frozen basis tree); `student()` taken on a sharer among two; the start rule applied at K = 1; the count split counting a stated forest, or subtracting it from the
total; `bart(n.trees = )` refused beside one stated count; `leaf.prior = normal(sd = )` still redirected with no
basis forest.

Gates, on the branch's own libraries:

- `R CMD INSTALL --preclean`, plain and `--configure-args=--enable-reference-build`, `buildInfo()$mode`
  checked; tests/cpp plain and under address and undefined sanitizers; the R-loaded sanitizer run on slice 2's
  tip, new numerics being reachable from R.
- The A/B, on both slices: a scratch script run on the reference build against an e2f37f6d library and the
  branch's, every recorded channel compared with `identical()`, two chains: one forest under gaussian, probit
  (the rescaling step on) and logistic; bcf drawn and held, under gaussian, probit and logistic; a forest with no
  basis second (its start at 0); three forests (none, a factor, a numeric column) under gaussian, probit and
  logistic, these binary K-forest chains being where the rescaling admission and the pooled k change code; a
  lone forest with no basis beside a forest with a basis, drawn by default and with `sd`, which keep main's path
  until slice F; the constant-column spelling with sd; a held forest with no basis; `$setLeafPrior(forests = )`; a warm start
  (`installTrees`) from one K-forest sampler into another; xbart; a copy and a `saveRDS`/`readRDS`
  continuation of each. All identical; a difference is a stop.
- The equivalence trio, `compare --bitwise --strict-coverage`: 55, 19 and 11 identical, 0 skipped, no
  `max |z|` line; the four snapshot files unchanged on the reference build.
- Exact posterior, new arms in quick mode, joining exact-gates.yaml: benchmarks/R/forest-prior-exact.R, held arms
  in closed form: two sharers beside a held factor, three forests with no basis one stating sd, two sharers
  alone at k 2 and at `normal(k = 3)`, and two sharers of unequal counts, each at T / sqrt(2);
  benchmarks/R/probit-k-scale-exact.R, an arm of two frozen trees sharing one drawn k, the second a stump
  integrated in an outer quadrature, an arm of the same two sharers beside a frozen single-tree factor forest
  held at a stated sd (the offset reading of D5), and an arm with a user offset, which alone sees the step drop it,
  each with the rescaling step on; the review ran them with the step removed too, from a library of its own, the arms
  passing, so the step keeps the target. The existing arms of both pass unchanged. bcf-exact.R gains no mode: the
  shared coefficient is gone.
- Statistical equivalence (a new script in benchmarks/R, z mode over seeds): `y ~ forest(x1) +
  forest(x2) + forest(x3)` against `y ~ forest(x1 + x2 + x3, blocks = )` with the same trees per group, gaussian
  at k 2 and probit under the drawn k; two forests at 75 trees (38 and 37) against `blocks()` at 75; and a
  stated count, `forest(x1, n.trees = 10) + forest(x2) + forest(x3)` against `trees.per.group = c(10, 38, 37)`:
  posterior mean fit by row, sigma, k, with two scenarios that see a k error, gaussian under a drawn k and unequal
  counts under one. With slice B in, the two spellings are one model at every count (dec-B506); a |z| past 4 is a
  stop. It joins exact-gates.yaml.
- The full tinytest suite, one process; `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`,
  `tools/check-win-drift.R`, `tools/check-doc-freshness.R`, each on its own exit status; `inst/NEWS.Rd` parses;
  `R CMD check --as-cran` from a tarball of a clean copy staged outside the tree; exact-gates.yaml's list in quick
  mode (a fit object carries new k on sharers).
- No `bench-sampler.R` compare: the new code runs once per sweep, never per node or row.

## Interaction with the other pending plans

forest-after-merge (8adc5d1a): dec-B508 removed its long form, `forest(leaf.prior = )`, and its first slice;
`forest()` states a size as `sd` alone, so this plan's long-form notes and the side-by-side pin they asked for are
gone. The fitting function's `leaf.prior = normal(sd = )` is the group's total (R item 4) and `forest(sd = )` one
forest's own size (D4); both help pages say so. Shared file: [`forestMapEntry`](../../src/bartcore/combiner.hpp)
if its remaining slices add a call site; whichever lands second passes the share there. This item adds no facade
virtual. Landing: independent; whichever is green first lands first.

coefficient-prior-opt-in (its branch, amended with this plan for dec-B509 and dec-B510): it lands first. It
accepts `amplitude = student(df)` on a forest with no basis, `student(1)` there being main's drawn default operand
for operand, so bartCause can state bcf's prior before the default it restates goes away. This plan then carries
the flip ([Slice F: the default flip](#slice-f-the-default-flip)), the median of a forest with no basis stating
`student()` and no sd (its engine item 3), and the refusal of `student()` on a sharer among two or more (D3).
Shared code: [`forestMapEntry`](../../src/bartcore/combiner.hpp) (its df and median, this plan's share and the
kind-none hold), [`validateForestAmplitude`](../../R/model.R) and benchmarks/R/forest-prior-exact.R; this plan
rebases onto it. Its pins of `student(1)` on a forest with no basis against the pre-flip default are re-pinned by
slice F against the stated prior.

bartCause: its `bcf()` builds `forest(vars = , n.trees = , sd = sd.control, amplitude = NULL or fixed())` for the
prognostic forest, drawing under dbarts's default and stating sd only where it holds it under gaussian. Before
slice F lands it states bcf's prior: where it draws, `amplitude = student(1)` and `sd` 2 sd(y) under gaussian or
the link's unit under probit and logistic (1, pi / sqrt(3)), the default it has today; where it holds, the
link's unit under probit and logistic as well, the held default moving to the drawn k there. Under probit and
logistic its draws stay bit for bit; under gaussian they differ in rounding, R's sd(y) and dbarts's unit differing in the last place
(run). That change is a commit on bartCause's own branch, after coefficient-prior-opt-in is on dbarts's main;
bartCause's tests read shapes and pass either way. Order of landing: coefficient-prior-opt-in; bartCause's
`bcf()`; this item, slices B, 1, 2 and F, merged once.

setmodel-parameters-not-structure (its slices 1 to 3 independent of this item): its slice 4 lands after this
item and is planned against D7 and D10, `$setModel`'s tree prior reaching every forest with no basis that states
none, and the group's leaf prior set through `$setModel`. Engine features merge serially.

## Order of work

0. Slice B, `blocks()` per group (~580 lines: engine ~110, tests/cpp ~150, blocks-exact.R ~200, tinytest ~70,
   help and design ~50; ~7 agent-hours). Tier: changes draws (a prior; posterior-changing for `blocks()` fits
   with unequal groups, neutral bitwise for all else). Built and reviewed first, on its own: it is independent of
   sharing, smaller, and slice 2's equivalence against `blocks()` at unequal counts needs it. Gates as
   [Slice B](#slice-b) lists.
1. Slice 1, engine and bridge (~870 lines, tests/cpp ~400 of them). Tier: changes draws (engine). RNG: neutral;
   nothing in R reaches a sharer yet. Gates: tests/cpp with sanitizers, the A/B, the equivalence trio, the
   snapshots, the full suite.
2. Slice 2, R surface and the new prior (~1470 lines: R ~280, help and NEWS ~130, tinytest ~500, gate arms and
   the equivalence script ~400, docs ~160). Tier: changes draws (a new prior). RNG: neutral for every call main
   accepts (the A/B again on the final tip); posterior-defining for the new models, gated by the exact arms, the
   equivalence against `blocks()` and the design note. Then the one review over slices 1 and 2, the fix round and
   the whole gate list. Slices 1 and 2: ~24 agent-hours. They may be built while coefficient-prior-opt-in is
   built; neither waits on the other until slice F.
F. Slice F, the default flip (~590 lines; ~8 agent-hours). Tier: changes draws (a default). RNG:
   posterior-changing for main's model with a forest with a basis. Built after a rebase onto a main carrying
   coefficient-prior-opt-in; its own opus review with mutants; gates and re-records as
   [Slice F: the default flip](#slice-f-the-default-flip) lists. The branch merges only after bartCause's `bcf()`
   change is on bartCause's main.
3. Slice 3, the constant column (~60 lines: R ~15, help ~10, tinytest ~30, records). Tier: R surface, no draw
   change; sonnet review. Built under dec-B528. Independent of slices B, 1 and 2. ~2 agent-hours.
4. Records: the ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, the design
   note's Status, TODO, and dec-B506's record line pointing at slice B's landing; dec-B509's at slice F's, and
   dec-B293 marked revised in part. One merge at the end.

## Cost

- The engine gains a second meaning of k on the K-forest chain (the group's, on every sharer), and the rescaling
  step reads other forests' fits as an offset.
- A forest's count and size now depend on how many forests with no basis the model has, by design; adding one
  takes trees from the others (dec-B503's cost), ten at the default having 7 or 8 each.
- A leaf scale per sharer (dec-B505) and per `blocks()` group (dec-B506): the pooled k and the level shifts read a
  scale per forest or per tree, and a `blocks()` fit with unequal groups draws differently from main.
- A model of forests with no basis inherits every refusal of a model of several forests: no test predictors at
  fit time, no variance forest, gaussian, probit and logistic only, constant leaves only. `bart(x.test = )` users
  of an additive model must predict after.
- Slice F changes the default posterior of every model of two or more forests in which a forest with no basis states
  no sd or states sd with no `amplitude`, none released (dec-B509's cost): one equivalence baseline re-recorded by the ruling's count and a
  second, bcf-equivalence, re-recorded for rounding; about twelve tinytest files re-pinned. The check found a
  held prognostic forest not clearly worse for the effects, and f's coverage about half a point lower.
- Under probit and logistic, `amplitude = fixed()` on a forest with no basis now means bart()'s drawn k, not 1
  unit held (dec-B293 superseded in part); a caller wanting that size states `sd`.
- The probit rescaling step does not move a forest with a basis (D5): a binary causal model's k mixes better than
  without the step, by an amount unmeasured.

## Landing note

Owed before the merge: the fresh review of slice F's fix round; slice B's bench-sampler compare on a quiet machine;
the rebase, the x86 leg and CI on the tip; bartCause's branch bcf-tests-state-prior, its two tests' hand-written
drivers and a comment, merged with this item; the records of Order of work's item 4. Current baselines:
equivalence-1380accf and bcf-equivalence-1380accf, multinomial-equivalence-80b1c8d4 unmoved.

For the next item: slice F's re-recorded bcf-equivalence baseline holds until forest-amplitude-normal (dec-B542), which re-records
again, judged by the same tier 1 of `--cross-host`. That item also re-pins what test-bcf-family.R,
test-forest-sd-unit.R and test-forest-leaf-prior-reader.R already pin of a drawn basis forest's units; slice F added
none. forest-door's slice 3 takes the refusal of a sharer's sd (dec-B525).

## Calls made

1. RULED, dec-B509: the total is the fitting function's leaf prior, held, in every model (D1).
2. The engine reads the sharers from the statements it has; forests stay engine forests (D2).
3. A drawn k is one k for the group (D3). `student()` on a sharer among two or more is refused; no shared
   coefficient is built (D3).
4. A stated sd is the forest's own size and takes it out of the share; the rest share among themselves (D4).
5. Probit and logistic: the drawn k shared in every model (dec-B509); the rescaling step extended to every model
   whose sharers draw k, reading the other forests' fit as offset (D5).
6. The fitting function's tree prior is each sharer's; its `interactions` and `blocks` refused beside two or
   more (D7).
7. k and `k.scale`: `k.scale` / k is the forest's own sd; sharers report the group's k in every model (D8).
8. Holds, starts, restatement, families, the formula and reading as D10 lists them.
9. No new record and no state format change; a stored state's sharers hold their held 1.
10. No bench compare for slices 1 and 2: their new code runs once per sweep. Slice B takes one, its branch being
    per tree.
11. Under dec-B503, a stated count stands outside the total. Alternative: inside it, the stated count
    subtracted from the fitting function's, which stated counts can use up. Three of the four panelists named
    this reading unprompted.
12. The forests stating no count split all of `n.trees` among themselves, so stating a count on one forest never
    changes the others'. Alternative: they split what is left, the reading 11 rejects.
13. `bart(n.trees = )` beside a stated count is taken, governing the unstated forests, and refused by name only
    where every forest with no basis states one, which keeps today's one-forest refusal. Alternative: refused
    beside any stated count, as on one forest today, which leaves the unstated forests' count settable only by
    stating each.
14. Print shows each forest's tree count and the total, on every model of several forests. Alternative: the
    counts in the synopsis only, as today; all four panelists asked for this.
15. The rescaling step, admitted on a K-forest chain, reads every other forest's fit as offset, the combined fit less
    the sharers' (D5); the combined-fit cache has no reader but the fill, so it is neither scaled nor refilled. Alternatives:
    scale every forest's contribution, ~60 more engine lines, a fuller move, unmeasured; or decline beside a
    forest with a basis, which leaves the default probit bcf model's k mixing as a lone forest's did before the
    step.
16. Slice B stores a group factor only where the group counts are unequal, so every other fit takes today's code
    operand for operand and the baselines stay bitwise. Alternative: a factor on every `blocks()` forest, 1 where
    the counts are equal, which moves no posterior but rounds differently and would move every `blocks()` fit's
    draws.
17. Slice B sets a tree's scale for the length of its update, so the leaf models and the rule scan read it unchanged,
    and fixes the three reads outside a tree's update by name. Alternative: a scale per tree threaded through every
    leaf-model call, a wider diff through model.hpp for the same result.
18. Slice B lands first, reviewed on its own, ahead of slices 1 and 2. Alternative: folded into slice 1, which
    makes slice 1 posterior-changing for main's `blocks()` fits and puts both under one review.
19. The flip is its own slice, F, after slices 1 and 2, so their bitwise A/B against main stays valid and the one
    posterior change of slices 1 to 3 has its own review and gates. Alternative: in slice 2, which mixes a
    neutral slice with a default change under one A/B.
20. A lone sharer beside a forest with a basis whose k is held at 2 keeps dec-B290's held map entry, so held
    gaussian fits stay bitwise (engine item 4).
21. A leaf prior given to the fitting function where no forest shares is taken silently when it restates the
    default (`normal(2)`, `k = 2`), as main takes it, and refused by name otherwise (R item 4). Alternative:
    refuse any, which breaks three exact gates and refuses a harmless restatement.
22. bcf-equivalence's helper restates bcf's prior with R's sd(y), its gaussian drawn scenarios compared with
    `--cross-host`, tier 1 with varcount identical, and re-recorded (Slice F). Alternative: twice the engine's unit
    read from a probe sampler, bitwise.
23. A forest with no basis stating `student()` and no sd takes as median the prior median of the size it would
    otherwise have, the chi prior's median k where k is drawn (Slice F engine item 3); set with
    coefficient-prior-opt-in, which records the alternatives.
24. A lone sharer's `$setLeafPrior(forests = list(forest(sd = )))` is taken, leaving the leaf prior with k pinned at
    1 (dec-B525's end state); the refusal's reason, moving others' shares, does not apply to one. The model's leaf
    prior goes back to a fixed k of 2 when no forest takes it, so a copy passes the bridge. Alternative: refuse it
    as among two or more.
25. `forestMapEntry` and `defaultForestSd` lose their family argument: after the flip nothing reads it.
26. Engine item 4's exception keys on a held k equal to the default and no named sd, with no branch on the family:
    probit at a stated k of 2 also takes the held map entry, the same prior.
27. A lone sharer reads `prior.sd.of` "forest", not "forest, shared", and has no group line in print or verbose; its
    print line is "no basis; sd X (default)", "(leaf prior)" where k is not the default, or "no basis; k chi(1.5, 2)".
    Its reader's amplitude is `fixed()`, as a sharer's among two or more (the plan's test wording said none; slice 2
    built it so).
28. `student()` with no sd under a flat chi prior or an infinite k is refused by name, having no median, the engine's
    throw as a backstop.
29. The reader's `leaf.prior` written back through `forests =` is nothing written; `$setLeafPrior(leaf.prior)` on a
    lone sharer takes only a restatement of the one in force, so `normal(k = 2)` is refused under probit, naming
    `forests =`.
30. `fit$fixed$amplitude` is read from the reader, so a forest stating sd, now held, is listed.
31. A leaf prior with no forest to take it is refused in new words, "no forest of this model takes a leaf prior's
    'k'" and "... 'sd'": the old ones, "a model of several forests does not support a non-default 'k'" and "...
    states each forest's size on the forest", are false once such a model takes both where a forest has no sd.
    A multinomial model, which declares no bases, keeps the old 'k' words.
32. `student()` with no sd under `normal(sd = invchi(df, s))` has as median that prior's, s over the chi's median, as
    the same prior written `normal(k = chi(df, k.scale / s))` has. As first built it read s, which gave one prior two
    medians (0.976 and 1.501 on the review's call).
33. A lone forest's print line says "(default)" only where its family holds k at 2 by default, gaussian, and "(leaf
    prior)" at `normal(k = 2)` under probit and logistic, whose default is drawn. A forest with no basis stating
    `student()` and no sd says "(leaf prior)" at every leaf prior: that is true of each, and neither a fit's record nor
    the reader says whether a leaf prior was stated. Alternative: "(leaf prior)" on every lone forest, or a new field
    on the reader to tell a student() forest's default apart.
34. The sd in force written on a lone sharer under `normal(sd = )` re-derives its entry, so it keeps the prior to
    rounding, 2.4e-15 gaussian and 3.1e-15 probit, not bitwise; the help says so. At a held k of 2 and 3 it is identical.
    Alternative: a branch keeping the entry, for no difference a posterior can see.

## Open calls

1. RULED, dec-B503 (VD 2026-10-29: "Ok, use option 2."). The tree count of forests with no basis is the group's,
   split as evenly as `blocks()` splits it; a forest that states `n.trees` keeps its own. The alternative, the
   full count on each forest, was not taken. The calls made under it are [Calls made](#calls-made) 11 to 14.
2. RULED, dec-B528 (VD 2026-10-30: "Use your recommendation for this item."): option (a), built in slice 3 (D9). The
   alternatives, (b) and a warning, were not taken.
   The constant basis column, `forest(x2, basis = one, sd = 2)` (dec-B307 said to decide here):
   - (a) Recommended and taken. Refused at creation, sd stated or not, the message naming `forest(x2)` with no basis.
     A treatment column left constant by a subset is caught, and one model has one spelling. Cost: a forest with
     its own sign-symmetric coefficient on a constant can no longer be written; a later swap to a constant
     column is still taken (dec-B346).
   - (b) Kept as built, created where sd is stated. It is a different model from the plain spelling (its own
     N(0, 1) coefficient, outside the share). Cost: the subset case fits with no message.
3. RULED, dec-B505 (VD 2026-10-29: "The motivation behind setting different numbers of trees is to express that
   different forests are likely to have to represent different levels of complexity in the relationships among
   their inputs. That doesn't mean that each forest shouldn't attempt to capture the same amount of residual
   variation."). How the shared total divides among sharers of unequal counts: an equal part per forest, 1 / K,
   each forest's leaf scale dividing its part among its own trees. The alternative, a part per tree, which this
   plan had recommended, was not taken. dec-B506 (VD: "Change blocks too.") changes `blocks()` to match, each
   group 1 / G of its forest's prior: [Slice B](#slice-b).
4. RULED, dec-B513 (VD 2026-10-29: "Use your recommendation for this item."): option (a), built in slice F. The
   alternative, (b), was not taken.
   `amplitude = fixed()` on a lone forest with no basis (`y ~ forest(x1, amplitude = fixed())`, or `bart()` with
   one `forest()` term). After slice F every forest with no basis is held by default, so `fixed()` there restates
   the default. Today it is refused: "this model has one forest, which has no coefficient to hold; 'amplitude'
   needs a model of several forests" (run). Beside other forests it is taken (D10).
   - (a) Recommended and taken. Accepted on a lone forest, the default restated, with no message: one rule for
     `fixed()` on every forest with no basis, and a call moved from a model of several forests to one keeps
     working. Cost: a word that does nothing is accepted; ~10 lines and a pin.
   - (b) Refused as today. Cost: `fixed()` is taken on a forest with no basis in one model and refused in another
     that differs by a forest.
   `student()` on a lone forest stays refused under both: the single-forest chain has no coefficient to draw.

## Not verified

- The equality of slice 2's models with `blocks()` after slice B is argued from the leaf scales and the sweep order;
  forests-blocks-equivalence.R gates it by z over seeds, which sees a difference and does not prove none.
- The probe is pure R at 400 draws (the drawn-coefficient quantiles carry Cauchy noise of several percent); it
  does not use the engine.
- How well the shared k mixes under logistic before dec-B404's move is built; coverage at equal time is the
  measure, as for one forest.
- How well the drawn k mixes beside a forest with a basis under probit, with the step reading that forest as
  offset, and under logistic with no step; the check measured gaussian alone.
- The effect of slice B on main's `blocks()` fits at the default splits is computed on the prior (the table in
  [Slice B](#slice-b)), not measured on a posterior.
- The agent-hours are estimates from comparable landings.
