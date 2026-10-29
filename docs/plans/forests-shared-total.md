# forests-shared-total: several forests with no basis share one prior total

Status: PLANNED 2026-10-29; two calls held for the maintainer ([Open calls](#open-calls)). Ruling: dec-B292
(with dec-B291, dec-B274, dec-B307); window "Post-merge, pre-1.0" (VD 2026-10-07).

agent: one opus implementer for slices 1 and 2; one opus review over both that runs the mutants below, the same
reviewer checking the fix round unless the fixes change engine behavior. Tier "Changes draws"
([Process by risk](README.md#process-by-risk)): slice 1 edits the engine's calibration map, coefficient draw and
k draw; slice 2 opens a new prior. Slice 3 is R surface, one sonnet review.
rng: NEUTRAL, bit for bit, for every call main accepts: no such call has two forests with no basis, and every
new branch is taken only where two or more do. Gate: a bitwise A/B against an e2f37f6d build, the equivalence
trio 55, 19 and 11 identical (equivalence-32418b6e, bcf-equivalence-20eed7ed, multinomial-equivalence-80b1c8d4),
the four snapshot files unchanged, no re-record. New behavior, with nothing on main to compare against: a model
of several forests with no basis. That is a new prior (posterior-changing class for those models): new arms of
three exact-posterior gates, a statistical equivalence against the one-forest `blocks()` model it equals, and a
design note.
window: after the merge to main, before 1.0-0 (dec-B292). Additive: until it lands the model is refused.
budget: ~2150 lines in the tree: engine ~260, bridge ~70, R ~260, help and NEWS ~130, tests/cpp ~350, tinytest
~480, exact-gate arms and the equivalence script ~380, docs and records ~220 (slices: ~680, ~1410, ~60). Out of the tree, ~220 lines of
scratch: the A/B script.

## Summary for the maintainer

Today a model takes one forest with no basis. `y ~ forest(x1) + forest(x2) + forest(x3)` is refused, and the
way to write forests that simply add is a basis column of ones with a stated size. After this item that formula
fits. The forests with no basis share one prior size the way the trees of one forest do: together they have
the size a single such forest would have, and each gets an equal part of its variance, so with three forests
each is about 0.58 of the total. A forest with a basis keeps its own size, as you ruled.

What "the size one such forest would have" means depends on the model. With no forest that has a basis, it is
`bart()`'s own forest: its `k`, a stated `leaf.prior`, and under probit and logistic its drawn `k`. The model is
then exactly `bart()`'s one forest with its trees divided among the formula's terms, which is what `blocks()`
already does inside one forest. Beside a forest with a basis, it is what the forest with no basis has there
today: a drawn size whose median is twice the response's sd, or the link's unit under probit and logistic. A
size that is drawn is drawn once for the group, so the forests grow and shrink together, as trees do. A forest
that states its own `sd` has that size and leaves the share.

The default sizes were checked by simulating the prior. Each forest at the full size instead makes the sum
sqrt(K) times as wide: at three forests under probit at a fixed size, 38 percent of the prior mass goes to
probabilities past 0.01 or 0.99, where one forest puts 14 percent.

Two questions for you, below: whether `bart()`'s tree count goes to each forest or is split among them (I
recommend split, the cost of one forest), and whether a constant basis column, today's workaround, should now
be refused in favor of the plain spelling (I recommend refusing it).

## Goal

A model may hold any number of forests with no basis, beside forests with a basis or alone, at every door that
takes several forests: a formula's `forest()` terms, a `forests` list and a data object's `bases`. The forests
with no basis that state no `sd` share one prior total, each an equal part of its variance; the total is the size
one forest with no basis has in the same model; a drawn total is one draw for the group. Every model main
accepts draws what it drew.

## What exists today

Read, and run on main at e2f37f6d (shipped build, arm64 macOS; the probes are in scratch/fst-plan, untracked):

- The refusal. [`plainForest`](../../R/model.R) refuses a second forest with no basis ("a model has one forest
  with no multiplier"), and the formula reader refuses it first
  (["forest() terms with no basis"](../../R/formulaTerms.R); run: "the formula has 2 forest() terms with no
  basis"). Pinned in ["twoPlain"](../../inst/tinytest/test-forest-defaults.R) and
  ["the formula has 2 forest() terms with no basis"](../../inst/tinytest/test-formula-terms.R). Plain terms
  beside a `forest()` with no basis are refused as two such forests.
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
  does not support a non-default 'k'"; `leaf.prior = normal(sd = 1)` is sent to the forest's `sd`).
- The workaround: `forest(x2, basis = one, sd = 2)`, `one` a column of ones, is a numeric forest with a drawn
  N(0, 1) coefficient starting at 1 and its own size (dec-B307; run: created, `sd` 2, `k.scale` 2.97).
- `blocks()` ([`blocks`](../../R/model.R)) confines each tree of one forest to a group of predictors, the trees
  split evenly among the groups by default; every tree keeps the forest's leaf sd. That is an additive model
  whose groups share the forest's total in proportion to their trees.
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
- One total split evenly (the ruling): forest j's leaves N(0, T^2 / (K m_j)), so each component has sd
  T / sqrt(K) (0.71 T at two, 0.58 T at three, 0.5 T at four) and the sum T exactly, at every x. With equal tree
  counts this is one forest of K m trees with its trees in groups: the `blocks()` model.
- A drawn total. One coefficient a ~ Cauchy(0, c) on the sum, f = a (F_1 + ... + F_K), each F_j at U / sqrt(K):
  the sum has the one forest's prior exactly. Each forest its own a_j ~ Cauchy(0, c_K) instead: given the F's, the
  sum at x is Cauchy with scale c_K (|F_1| + ... + |F_K|), which no c_K makes the one forest's; c_K = c / sqrt(K),
  the variance split, leaves it sqrt(K) too wide, and c / K matches it in scale only.

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

D1. Which total. "The size one forest with no basis has in the same model" (dec-B292: "whose sum has the prior
one such forest has"). Where no forest has a basis that forest is `bart()`'s one forest, so the sum takes the
fitting function's leaf prior: `k`, `leaf.prior = normal(k = )` or `normal(sd = )`, and under probit and
logistic the drawn k under chi(1.5, 2). Beside a forest with a basis it is the forest with no basis there today:
drawn, median 2U under gaussian and U under probit and logistic; held, a quarter of the range or U.
MY CALL. Alternative: the second for every model, which makes `y ~ forest(x1) + forest(x2)` a heavier-tailed,
wider prior than `y ~ forest(x1)` (half-Cauchy median 2 sd(y) against a fixed quarter of the range, about 1.4
sd(y) at n = 300) and leaves `bart(k = )` refused for it. The cost of the call: the two kinds of model differ in
what their forests with no basis share, as one such forest already differs today between a lone forest and one
beside a forest with a basis.

D2. How the engine tells them apart. The sharers are the forests of kind none that state no sd, read from what
every forest already states ([`ForestStatement`](../../src/bartcore/combiner.hpp), its kind); whether the model
has a forest with a basis is read the same way. No new field, no flag from R, one engine function that answers
both. They stay separate engine forests, so every per-forest surface (labels, fits, counts, trees, weights, move
settings, tree priors) works as it does. MY CALL. Alternatives: an explicit group handed by R, which the ruling's
"no switch" leaves nothing to set; and one engine forest whose trees are the forests (a generalized `blocks()`),
which would need every per-forest knob and reader made per block.

D3. A drawn total is one draw. Beside a forest with a basis, the sharers take one coefficient: their blocks hold
one value, drawn once from its full conditional, whose design row is the sum of the sharers' fits, with one
half-Cauchy variance refreshed after it. With no forest with a basis, a drawn k is one k, drawn from the pooled
leaves of every sharer. MY CALL: it is what makes the sum exactly one forest's prior, and "as the trees of a
forest do" is one scale. Alternative: each its own coefficient at c / K, which keeps the engine as it is and
gives each component its own adaptive size, but matches the one forest only in scale (the table above) and has
no total that is shared.

D4. A stated `sd` on one of them. It is that forest's own size, as `sd` is everywhere (dec-B253): the size a lone
forest with no basis stating that sd has in the same model, exact where held or fixed, the median where its
coefficient is drawn, with its own coefficient or fixed k. It leaves the share, and the forests that state none
share the total among themselves. MY CALL. Alternatives: the stated number as the total the forest takes 1 /
sqrt(K) of, which no other `sd` means; and keeping K at every forest with no basis, so stating one forest's size
shrinks the others' total in silence.

D5. Probit and logistic. With no forest with a basis, the sharers take `bart()`'s drawn k under chi(1.5, 2),
one k for the group: dec-B420 found a fixed k not good enough for binary fits. The probit rescaling step
([`Chain::drawForestRescale`](../../src/bartcore/chain.hpp)), which mixes that k, is extended to the group: it
scales the latents, every sharer's leaves and k by one factor, the leaf count M the group's. It declines, as it
does today where it does not apply, when any forest states an sd, a fixed-size forest's leaves having no k to
absorb the factor. The logistic collapsed scale move (dec-B404), when built, covers the group. Beside a forest
with a basis the sharers take the forest with no basis's binary default, U drawn or held (dec-B293). MY CALL.

D6. The tree count. Held for the maintainer, [Open calls](#open-calls) 1. The plan is written for the
recommendation, the fitting function's count split evenly among the sharers, so a model of forests with no basis
costs one forest and at defaults is `blocks()` exactly.

D7. The tree prior, `interactions` and `blocks`. Each forest with no basis takes the fitting function's tree
prior where it states none, as the one does today. `interactions` and `blocks` given to the fitting function
beside two or more forests with no basis are refused by name, telling the caller to state them on a forest, as
they are beside none ([`refuseStatedWithNoPlainForest`](../../R/model.R)): a partition of one forest's
predictors cannot be each forest's. MY CALL. Alternative: give each forest the fitting function's, which works
for `interactions` and never for `blocks`.

D8. k and `k.scale` (forest-multiplier's call (2)). The rule that covers every forest: `k.scale` / k is the
forest's own prior sd before its coefficient, in the response's units. A forest whose size a calibration map
sets reports k 1 and `k.scale` that sd, as built; for a sharer beside a forest with a basis that is its share,
U / sqrt(K) by default. A sharer in a model with no forest with a basis reports the group's k, one value
reported on each, and `k.scale` its share times k, so `k.scale` / k is its share, and `getK` and `extract(type =
"k")` give the drawn k on each. MY CALL. Alternative: k 1 everywhere, which would hide the drawn k.

D9. The constant basis column (dec-B307's revisit). Held for the maintainer, [Open calls](#open-calls) 2. The
plan is written for the recommendation, refusing it at creation and naming the plain spelling, in slice 3.

D10. The rest of the surface, each MY CALL:

- Holds. With no forest with a basis there is no coefficient: every forest is held at 1, a forest stating sd
  at exactly that size, and `amplitude` is refused there as on a lone forest ("this model has one forest, which
  has no coefficient to hold" reworded for several). Beside one, the
  sharers are all held or all drawn, a mix refused by name; a forest stating sd holds or draws on its own.
- Starts. Where two or more forests have no basis, every drawn coefficient on such a forest starts at 1, the
  shared one included, so none starts out of the model; a model with one keeps its start by position. No model
  main accepts moves; TODO drawn-coefficient-start-by-kind is noted.
- A restated sd. `$setLeafPrior(forests = )` restating a sharer's sd is refused by name: membership is fixed at
  creation, and a sharer leaving would move the others' share. A forest that stated sd restates it as today, and
  `forest()` written back stays the no-op it is. `leaf.prior =` on a sampler of sharers stays refused as on every
  sampler of several forests (an addition later).
- Families and compositions. Gaussian, probit and logistic, and every refusal a model of several forests has
  today stands, test predictors and the variance forest among them. Each is an addition later.
- The formula. Plain terms beside a `forest()` with no basis stay refused: mgcv reads `y ~ x1 + s(x2)` with x1 a
  linear term, so `y ~ x1 + forest(x2)` has two readings. The message now says to write each as a `forest()`.
  A `forest()` with no predictors keeps its reading, every predictor of the fit.
- Reading. `getForestAmplitudes` on a sharer gives the shared value; `print` and `summary` mark a sharer's line
  as a share of the total, and the verbose creation summary prints the group once.

### R surface

1. [`plainForest`](../../R/model.R) returns every position with no basis and refuses nothing; its callers take a
   vector. The formula reader's count refusal goes; its plain-terms refusal is reworded (D10).
2. [`forestParams`](../../R/model.R): under OC1's recommendation the sharers that state no `n.trees` split the
   fitting function's count evenly, the first `n %% K` taking one more, as `blocks()` splits; a forest's own
   count governs it and is not counted in the split. Where forest 1's count is not the fitting function's, the
   forests' record keeps the count the next fit inherits, as it does where forest 1 has a basis
   ([Where the numbers live](../design/forest-defaults-by-kind.md#where-the-numbers-live)). Each sharer takes the
   fitting function's base and power where it states none.
3. `interactions` and `blocks` given to the fitting function beside two or more forests with no basis are refused
   (D7), in [`refuseStatedWithNoPlainForest`](../../R/model.R)'s style.
4. [`resolveSamplerSpec`](../../R/spec.R): where no forest has a basis and two or more are given, the k
   hyperprior and a k other than 2 are no longer refused, and a fitting function's `leaf.prior = normal(sd = )`
   states the group's total instead of being sent to a forest; where every such forest states sd, a leaf prior
   given to the fitting function has no forest to belong to and is refused by name. Everything else it refuses
   for several forests stands. Beside a forest with a basis nothing changes, and the redirect message names "each forest" where
   several have no basis.
5. Holds (D10): `amplitude` refused where no forest has a basis; a mixed hold among the sharers refused.
6. `$setLeafPrior(forests = )` refuses a sharer's sd ([`resolveForestSpreads`](../../R/dbarts.R)).
7. [`reportLeafPrior`](../../R/dbarts.R) and the k readers follow D8; `prior.sd.of` gains ", shared" on a sharer;
   [`forestPriorLines`](../../R/forestPrint.R) marks the share. A fit of sharers with a drawn k stores it as a
   one-forest fit stores its k.

### Bridge

1. Where no forest has a basis, the K-forest creation hands the chain the leaf prior the single-forest path
   reads (k, whether it is drawn and its chi prior, a named sd), which today it drops behind R's refusal.
2. `bartcore_getLeafPrior` reports each forest's k and the share its map entry carries; the verbose summary
   ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)) prints the group.
3. Nothing in [`applyForestBases`](../../src/R_interface_bartcore.cpp) or the data record changes: a sharer has
   no basis and no scale.

### Engine

1. Sharing. One function reads the specs: the sharers (kind none, no sd stated) and whether any forest has a
   basis. [`forestMapEntry`](../../src/bartcore/combiner.hpp) takes the sharer count K and divides a sharer's
   node scale by sqrt(K), held or drawn; the drawn coefficient's median is the group's, c, unchanged. K = 1 is the
   expression as it stands, operand for operand.
2. The shared coefficient. [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) keeps one block per
   forest and gains, for each forest, the first sharer it reads. [`drawAmplitudes`](../../src/bartcore/combiner.hpp)
   draws the group's value at the first sharer's place in index order, from the conditional whose design row is
   the sum of the sharers' fits against the residual net of every other forest, refreshes the one variance, and
   writes the value to every sharer's block; later sharers are passed by. One sharer or none is today's loop.
3. Starts (D10) at construction, after the layout: every drawn kind-none block at 1 where two or more forests
   are of kind none.
4. A model with no forest with a basis. The K-forest [`Chain`](../../src/bartcore/chain.hpp) constructor takes
   the options' k, its draw and chi prior, and the single-forest chain's node scale (the expression
   [`singleForestDefaultSize`](../../src/bartcore/chain.hpp) and [`translateNamedSd`](../../src/bartcore/chain.hpp)
   rest on, called and not copied), in place of the map's k = 1, for every sharer; a forest stating sd stays on
   the map.
5. The pooled k. The k draw after each sweep draws one k for the group from the sum over sharers of each forest's
   sum of squared leaves over its own leaf scale squared, and the sum of their leaf counts, on unit scale
   ([`ChiKHyperprior`](../../src/bartcore/chain.hpp)), and writes it to every sharer. Any other forest draws as now.
6. The probit rescaling step (D5): [`Chain::forestRescaleApplies`](../../src/bartcore/chain.hpp) admits a K-forest
   chain whose forests are all sharers with a drawn k and constant leaves, and the step scales each sharer's
   leaves and cache under the existing bound and rule; M is the group's leaf count.
7. State. The coefficient layout and [`glueIsValid`](../../src/bartcore/combiner.hpp) do not change: each
   sharer's block holds the shared value. An install takes the first sharer's value for all (no sampler writes
   others). `stateFormatVersion` does not move.
8. [`ForestCalibration`](../../src/bartcore/chain.hpp) carries whether a forest shares, the count K and, for a
   sharer of a model with no forest with a basis, k in force.

No facade virtual changes. The shipped header ([`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h)) has no
create entry and does not change; stan4bart builds one forest.

## Fit and state

- No new record: the group is read from the forests' statements and kinds at every creation, so `copy()`, a
  reload and a dead-pointer `setState` rebuild it.
- A stored state holds a value on every sharer's block, equal; a state installs across samplers as today.
- A control taken from such a fit carries the forests' record as today, with the inheritable count of R item 2.

## Help, NEWS and docs

- man/forest.Rd: several forests with no basis, what they share and how; `sd` on one takes it out of the share;
  `amplitude` and the holds; the count split. man/dbarts.Rd and man/bart.Rd: the `forests` and formula items,
  `n.trees` and `k` beside several forests with no basis. man/dbartsSampler-class.Rd: `getLeafPrior` and `getK`
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

- The map: a sharer's node scale is the lone forest's over sqrt(K), K = 2 and 3, held and drawn, gaussian and
  probit; a forest stating sd beside sharers keeps its own; K = 1 bitwise the old expression.
- The shared coefficient: every sharer's block equal after each sweep; the draw against a closed-form conditional
  on a fixture; one sharer draws bitwise what the old loop drew.
- The pooled k: one k on every sharer; the draw's inputs against a hand sum.
- The rescaling step on a two-sharer probit chain: admitted, declined with a stated sd, M the group's.
- Starts at 1 with two kind-none forests; by position with one.
- Fixtures that build two forests of kind none today are re-read: each states its kind or moves to the rule with
  its pins re-derived, reported.

tinytest, a new test-forests-shared-total.R, small n:

- Created at every door: a formula, a `forests` list, a data object's `bases`; gaussian, probit, logistic.
- The reader: each sharer's `sd` and `k.scale` the lone forest's over sqrt(K), K = 2 and 3, drawn and held; a
  stated forest's its own; D8's k on each.
- Equal by construction: with no forest with a basis and k fixed, the sharers' leaf scales equal `blocks()`'s one
  forest of the same trees (to 1e-12).
- The shared coefficient: `getForestAmplitudes` equal on every sharer at every draw.
- Counts (OC1): 75 split 25, 25, 25 and 38, 37; a stated count beside unstated ones; `print`'s counts.
- Refusals, each by name: `interactions` and `blocks` on the fitting function; `amplitude` with no basis forest;
  a mixed hold; a sharer's restated sd; plain terms beside a `forest()` with no basis (new text); a family or
  composition several forests refuse.
- Copy, `saveRDS` and `readRDS` after `storeState`: the next draws `identical()` to the live sampler's.
- Pins updated: the two refusals named in [What exists today](#what-exists-today) become creations.

Mutants, each to fail a named expectation: the share over K rather than sqrt(K); K counting forests that state
sd; the shared draw writing the first sharer only; the design row of one sharer; the variance refreshed per
sharer; the pooled k from the first forest's leaves; the rescaling step's M of one forest; the start rule applied
at K = 1; the count split counting a stated forest; `leaf.prior = normal(sd = )` still redirected with no basis
forest.

Gates, on the branch's own libraries:

- `R CMD INSTALL --preclean`, plain and `--configure-args=--enable-reference-build`, `buildInfo()$mode`
  checked; tests/cpp plain and under address and undefined sanitizers; the R-loaded sanitizer run on slice 2's
  tip, new numerics being reachable from R.
- The A/B, on both slices: a scratch script run on the reference build against an e2f37f6d library and the
  branch's, every recorded channel compared with `identical()`, two chains: one forest under gaussian, probit
  (the rescaling step on) and logistic; bcf drawn and held; a forest with no basis second (its start at 0);
  three forests (none, a factor, a numeric column); the constant-column spelling with sd; a held forest with no
  basis; `$setLeafPrior(forests = )`; xbart; a copy and a `saveRDS`/`readRDS` continuation of each. All identical;
  a difference is a stop.
- The equivalence trio, `compare --bitwise --strict-coverage`: 55, 19 and 11 identical, 0 skipped, no
  `max |z|` line; the four snapshot files unchanged on the reference build.
- Exact posterior, new arms in quick mode, joining exact-gates.yaml: benchmarks/R/forest-prior-exact.R, held arms
  in closed form: two sharers beside a held factor, three forests with no basis one stating sd, and two sharers
  alone at k 2 and at `normal(k = 3)`; benchmarks/R/bcf-exact.R, a mode with two single-tree sharers under one
  drawn coefficient beside the treatment forest, the quadrature over a unchanged; benchmarks/R/probit-k-scale-exact.R,
  an arm of two frozen trees sharing one drawn k, the second a stump integrated in an outer quadrature, with the
  rescaling step on and off. The existing arms of all three pass unchanged.
- Statistical equivalence (a new script in benchmarks/R, z mode over seeds): `y ~ forest(x1) +
  forest(x2) + forest(x3)` against `y ~ forest(x1 + x2 + x3, blocks = )` with the same trees per group, gaussian
  at k 2 and probit under the drawn k: posterior mean fit by row, sigma, k. The two are one model; a |z| past
  4 is a stop.
- The full tinytest suite, one process; `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`,
  `tools/check-win-drift.R`, `tools/check-doc-freshness.R`, each on its own exit status; `inst/NEWS.Rd` parses;
  `R CMD check --as-cran` from a tarball of a clean copy staged outside the tree; exact-gates.yaml's list in quick
  mode (a fit object carries new k on sharers).
- No `bench-sampler.R` compare: the new code runs once per sweep, never per node or row.

## Order of work

1. Slice 1, engine and bridge (~680 lines, tests/cpp ~350 of them). Tier: changes draws (engine). RNG: neutral;
   nothing in R reaches a sharer yet. Gates: tests/cpp with sanitizers, the A/B, the equivalence trio, the
   snapshots, the full suite.
2. Slice 2, R surface and the new prior (~1410 lines: R ~245, help and NEWS ~120, tinytest ~450, gate arms and
   the equivalence script ~380, docs ~215). Tier: changes draws (a new prior). RNG: neutral for every call main
   accepts (the A/B again on the final tip); posterior-defining for the new models, gated by the exact arms, the
   equivalence against `blocks()` and the design note. Waits on OC1 for item 2's count. Then the one review over
   both slices, the fix round and the whole gate list.
3. Slice 3, the constant column (~60 lines: R ~15, help ~10, tinytest ~30, records). Tier: R surface, no draw
   change; sonnet review. Waits on OC2; under (b) it is the help sentence alone. Independent of slices 1 and 2.
4. Records: the ledger entry for the calls made, this plan's Status and Landing note, the INDEX row, the design
   note's Status, TODO. One merge at the end.

## Cost

- The engine gains a second meaning of k on the K-forest chain (the group's, where no forest has a basis) and a
  coefficient read from another forest's block.
- A forest's count and size now depend on how many forests with no basis the model has, by design.
- A model of forests with no basis inherits every refusal of a model of several forests: no test predictors at
  fit time, no variance forest, gaussian, probit and logistic only, constant leaves only. `bart(x.test = )` users
  of an additive model must predict after.
- Two models one formula apart differ in what is shared: adding a forest with a basis turns the sharers' total
  from `bart()`'s k into the drawn coefficient (D1), as it does for a lone forest today.

## Calls made

1. The total is the size one forest with no basis has in the same model: `bart()`'s leaf prior where no forest
   has a basis, the forest with no basis's default of a model of several beside one (D1).
2. The engine reads the sharers from the statements it has; forests stay engine forests (D2).
3. A drawn total is one coefficient or one k for the group (D3).
4. A stated sd is the forest's own size and takes it out of the share; the rest share among themselves (D4).
5. Probit and logistic: the drawn k shared, the rescaling step extended, the binary defaults kept beside a basis
   forest (D5).
6. The fitting function's tree prior is each sharer's; its `interactions` and `blocks` refused beside two or
   more (D7).
7. k and `k.scale`: `k.scale` / k is the forest's own sd; sharers with no basis forest report the group's k (D8).
8. Holds, starts, restatement, families, the formula and reading as D10 lists them.
9. No new record and no state format change; a stored state's sharers hold equal values.
10. No bench compare.

## Open calls

1. The tree count of forests with no basis. `bart()`'s `n.trees` (75) is today the one forest with no basis's.
   With three of them:
   - (a) Each takes the full count: 75 each, 225 in all. A forest's count does not depend on its neighbors, as
     each mgcv smooth and each stochtree forest has its own. Cost: K times the time of one forest at defaults.
   - (b) Recommended. The count is the group's, split evenly, 25 each, as `blocks()` splits a forest's trees; a
     forest's own `n.trees` is its own. The model costs one forest, and at defaults it is exactly `bart()`'s
     forest with its trees divided among the terms. Cost: adding a forest takes trees from the others, and at
     ten forests each has seven or eight.
2. The constant basis column, `forest(x2, basis = one, sd = 2)` (dec-B307 said to decide here):
   - (a) Recommended. Refused at creation, sd stated or not, the message naming `forest(x2)` with no basis.
     A treatment column left constant by a subset is caught, and one model has one spelling. Cost: a forest with
     its own sign-symmetric coefficient on a constant can no longer be written; a later swap to a constant
     column is still taken (dec-B346).
   - (b) Kept as built, created where sd is stated. It is a different model from the plain spelling (its own
     N(0, 1) coefficient, outside the share). Cost: the subset case fits with no message.

## Not verified

- Nothing is built. The equality of slice 2's default with the `blocks()` model is argued from the leaf scales
  and the sweep order and is a gate, not a measurement.
- The probe is pure R at 400 draws (the drawn-coefficient quantiles carry Cauchy noise of several percent); it
  does not use the engine.
- Whether the K-forest chain's gaussian sigma prior and starting sigma equal the single-forest chain's; the
  `blocks()` equivalence gate checks it.
- How well the shared k mixes under logistic before dec-B404's move is built; coverage at equal time is the
  measure, as for one forest.
- The cost of the shared coefficient's draw is not measured.
