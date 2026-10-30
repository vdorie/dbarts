# tree-prior-flat: a tree's depth is plain base and power at every door, and DART a value of split.probs

Status: LANDED 2026-10-30 (dec-A246), built in three slices, slice 1 with one review and its fix round, slices 2
and 3 with one review; no call was held. Rulings: dec-B512, with dec-B116, dec-B246 and dec-B252 (each revised
by it), dec-B241 and dec-B503 (the collision rule it extends), dec-B254, dec-B266, dec-B274, dec-B508, dec-B514,
dec-B515 and dec-B518 (the live door).

agent: slice 1 (the R that decides which model a call fits) one opus implementer (effort high) and one opus
review (effort high) that runs the mutants below and checks the fix round. Slices 2 and 3 (help, docs, benchmark
scripts, the two consumer branches) one sonnet implementer and one sonnet review (effort medium). Tier "R
surface, no draw change" ([Process by risk](README.md#process-by-risk)): no engine, bridge or C header line
moves; slice 1 takes the opus pair because a slip there is another fitted model with no message, as for
forest-defaults-by-kind's first push. The blind critique has run.
rng: NEUTRAL, bit for bit, for every call main accepts with the same values: the model's tree-prior record keeps
its classes and slots, and every door writes the same numbers into it. Calls main accepts and this item refuses:
a setting stated on the fitting function and on the forest with no basis, which main resolves silently in the
forest's favor; `tree.prior` on `bart()` and `xbart()`, where it never shipped; and `xbart(tree.prior = cgm())`,
which on main fits the default depth with no message. Calls main refuses and this item takes: a uniform
`split.probs` named beside forests that all have a basis, and a 0.9-34-saved `cgm` object given as `tree.prior`.
No call of dbarts changes its draws; nothing is re-recorded. In bartCause, `bartc(crossvalidate = TRUE, power = )`
changes its chosen k and its draws, the xbart defect being fixed under it (slice 3).
window: for 1.0-0; on main before setmodel-parameters' slice 4 and before forests-shared-total's slice 2
([Order of landing](#order-of-landing)).
budget: ~1580 lines in the tree. Slice 1 ~1140 (R ~600, tinytest ~540: a new file ~300, re-pins in about 39
files ~240). Slice 2 ~440 (help ~210, NEWS ~25, vignettes and README ~10, design docs ~80, benchmark scripts ~70,
the scenario respelled in equivalence.R included, records ~45). Slice 3, outside this tree: bartCause ~150,
stan4bart ~50. Out of the tree, ~400 lines of scratch A/B script.

## Summary for the maintainer

You ruled that a tree's depth is written as plain `power` and `base` everywhere, that split probabilities are
`split.probs`, a vector or `dart()`, and that `tree.prior` leaves `bart()` and is retired with a warning on
`dbarts()` and `dbartsSpec()`. This plan builds it.

Nothing a user fits with dbarts changes its draws. The tree prior a sampler holds inside stays the same object it
is today, so `$setModel`, copies, saved fits and the consumers' reads of `model@tree.prior` keep working. What
changes is what the doors take and what they refuse. A 0.9-34 `bart2(power = 3)` call works as it did in 0.9-34,
with no warning; a 0.9-34 `dbarts(tree.prior = cgm(power = 3))` call, or one handing back a tree prior saved by
0.9-34, gets one warning a session and fits the same model. `split.probs = dart()` turns on DART at `bart()`,
`dbarts()`, `dbartsSpec()` and `xbart()`.

Two defects on main go: `xbart(tree.prior = cgm(power = 3))` fitting the default depth, and a forest's `power`
silently beating the fitting function's. The second becomes a refusal by name, as `n.trees` given at both levels
already is. In bartCause, a cross-validated fit given `power` now cross-validates at that power, so its chosen k
and draws change; its NEWS says so.

The live door. setmodel-parameters' slice 4 had planned `$setTreePrior(tree.prior = cgm(...))`; under your
ruling it takes `power`, `base` and `split.probs` plainly, and only the fitting function's level (dec-B518). A
forest's own depth changes through `$setForest(forest, power = , base = )` and the variance forest's through
`$setVarianceForest`, both built by the forest-door item. Where a statement at the fitting function's level
reaches no forest, because every forest with no basis stated its own, it is refused with the exact call to
write, as creation refuses the same statement (dec-B515). Whether an unstated variance forest follows a live
change of the tree prior is that plan's open question with you (held decision #6); slice 4 waits on it.

bartCause and stan4bart drop their translations into `cgm()` and pass the plain arguments through; treatSens
does not change. Their branches merge right after this one.

About 1580 lines with tests and docs, and about 200 in the two consumers; about 24 agent hours. Plans here run
1.5 to 2 times low, so budget 36 to 48.

## Goal

Every door takes a tree's depth as `power` and `base` and its split probabilities as `split.probs`, a vector or
`dart()`; `tree.prior` is refused on `bart()` and `xbart()` with the plain spelling to write, and accepted on
`dbarts()` and `dbartsSpec()` with a once-per-session warning until 1.1-0; a setting stated on the fitting
function and on every forest with no basis is refused by name; and every call main accepts and this item keeps
draws what it drew.

## What exists today

Run on main at 83cefe7d (shipped build, arm64 macOS, private library scratch/libs/tpf-plan; probes
scratch/tpf-plan/q0.R to q6.R and the critic's scratch/tpf-crit/c1.R to c6.R, untracked, and the design agent's
scratch/tree-prior-design/p1.R to p10.R rerun there as q2, q5 and q8). Data: 120 rows, `x1` to `x3` uniform, a
0/1 `z`, a gaussian `y`.

- 0.9-34 (scratch/libs/cran): `bart2` takes `power`, `base` and `split.probs` plainly; `dbarts` takes
  `tree.prior` (default `cgm`, resolved by bare name, `cgm(power = 2, base = 0.95, split.probs = 1/num.vars)`) as
  its tenth formal, and no `power`; `xbart` takes `power` and `base` grids and no `tree.prior` or `split.probs`;
  `rbart_vi` takes `power` and `base`; `bart` is BayesTree's (`power`, `base`, `splitprobs`). No `dart` exists.
  Nothing exports `cgm`. A `cgm` object saved by 0.9-34 has slots `power`, `base`, `splitProbabilities` and no
  `splitProbabilitiesSpec`; main given it as `tree.prior` fails "no slot of name 'splitProbabilitiesSpec'"
  (c3a, c3b).
- `bart()` takes `tree.prior` as a formal; `power`, `base` and `split.probs` arrive through `...` and warn once
  ("'power' has left 'bart': ... write tree.prior = cgm(power = ) or tree.prior = dart(power = )"), the values
  used ([`resolveConsolidatedArgs`](../../R/tombstones.R), [`consolidatedArgReasons`](../../R/tombstones.R));
  [`buildSamplerPriors`](../../R/bart.R) rebuilds them into a `cgm()` call for `dbarts()` and refuses
  `tree.prior` beside any of them ([`refuseColliding`](../../R/model.R)).
- `dbarts()` and `dbartsSpec()` take `tree.prior = cgm`; `dbarts(power = 3)` is "unused argument 'power'
  passed to 'dbarts'"; `dbarts(tree.prior = NULL)` is "'tree.prior' must be a tree prior specification".
  [`parsePriors`](../../R/model.R) resolves the object; [`resolveSplitProbabilities`](../../R/model.R)
  normalizes a vector against the columns, with `num.vars` and `numvars` in its vocabulary, and a uniform vector
  to none. A prior read off a sampler has `splitProbabilitiesSpec` NULL and its resolved probabilities and
  `update.delay`, and handed back to `dbarts()` keeps both, the delay too under another burn-in (c3b).
- `dart(power = 2, base = 0.95, a, b, rho, alpha, update.alpha, update.delay)`
  ([`dart`](../../R/model.R)) builds a [`dbartsDartPrior`](../../R/A_class.R), which extends
  [`dbartsCGMPrior`](../../R/A_class.R). `cgm(power = 3)` stores base 0.95: an unstated half cannot be told from
  the default.
- A model's `tree.prior` slot holds the first forest's prior; the bridge reads `power`, `base`,
  `splitProbabilities` and, by class, DART's slots from it, and the verbose summary prints "power and base for
  tree prior" ([`printInitialSummary`](../../src/R_interface_bartcore.cpp)). `print(sampler)`, `print(fit)` and
  `summary(fit)` show no tree prior; `s$model` shows the slots. A DART sampler made with `dart(power = 3)` copies
  and reloads with the same record (q4).
- `forest()` and `varianceForest()` take `base` and `power` plainly ([`forest`](../../R/model.R),
  [`varianceForest`](../../R/model.R)) and no `tree.prior` or `split.probs` (R's unused-argument error).
- Defect 1, xbart. `xbart(tree.prior = cgm(power = 3, base = 0.5))` returns errors identical to `xbart()`'s at
  the same seed (1.17472, 1.126016 both), while `xbart(power = 3, base = 0.5)` returns 1.12657, 1.120871 (q5):
  [`cellModel`](../../R/xbart.R) writes the grid's `power[1]` and `base[1]`, at their defaults, over the
  object's.
- Defect 2, two levels. `bart(y ~ forest(x1 + x2 + x3, power = 1) + forest(x2, basis = factor(z)),
  tree.prior = cgm(power = 3))` fits power 1, no message; the same at `dbarts(forests = list(forest(power = 1),
  forest(basis = factor(z))), tree.prior = cgm(power = 3))`; with the forest with no basis written second, power 1
  and base 0.5 from `cgm(power = 3, base = 0.5)` (q3); a lone term, `bart(y ~ forest(x1 + x2, power = 1),
  tree.prior = cgm(power = 3, base = 0.5))`, fits power 1, base 0.5 (c1). The writer is
  [`resolveSamplerSpec`](../../R/spec.R), which writes the first forest's knobs over the model's prior. `n.trees`
  given to `bart()` and to that forest is refused: "'n.trees' is given to the fitting function and to the
  forest with no basis ('forest(x1 + x2, n.trees = 10)'), which are the same count; give one"
  ([`refuseGivenTwiceInFormula`](../../R/formulaTerms.R)); at `dbarts()` a control's count beside the forest's
  leaves the forest's to govern.
- No forest with no basis. `dbarts(y ~ forest(x1, basis = factor(z)) + forest(x2, basis = ~ x3), tree.prior =
  cgm(power = 3))`: "'tree.prior' given to the fitting function is the tree prior of the forest with no basis,
  and every forest of this model has a basis; state it on a forest, as forest(x1, basis = a, base = 0.25, power =
  3)" ([`refuseStatedWithNoPlainForest`](../../R/model.R), dec-B274). A `split.probs` named there is refused the
  same way, a uniform `c(x1 = 1, x2 = 1)` included (pinned in test-forest-defaults.R), with a remedy no forest can
  take.
- Several forests refuse DART and a resolved non-uniform `split.probs` ("a treatment forest does not support a
  DART tree prior"; "... 'split.probs'") and take a uniform one: `bart(dbartsData(..., bases = list(NULL,
  factor(z))), split.probs = c(1, 1, 1))`, `bart(y ~ forest(x1 + x2) + forest(x3, basis = factor(z)),
  split.probs = 1/num.vars)` and `bartBT()` on such a data object at its default `splitprobs = 1/numvars` all fit;
  `bartBT(splitprobs = c(1, 2, 1))` there is refused (c5, c6). `family = "multinomial"` refuses `split.probs` by
  naming at `bart()` and DART in the specification. A variance forest cannot join several forests ("a treatment
  forest does not support 'variance'").
- The variance forest. An unstated one takes the resolved record's base and power at creation
  ([`resolveSamplerSpec`](../../R/spec.R)): 0.5 and 3 under `cgm(power = 3, base = 0.5)`, 0.95 and 3 under
  `dart(power = 3)`, and a lone forest term's own statement, 0.5 and 1 under `forest(x1 + x2 + x3, power = 1,
  base = 0.5)` with no fitting-function prior (c2); split probabilities never reach it (q3).
- `xbart()` takes `power` and `base` grids, `split.probs`, and `tree.prior` (`cgm()` or `dart()`);
  `xbart(split.probs = dart())` is "could not find function 'dart'" and `xbart(split.probs = c(2, 1, 1)/num.vars)`
  "object 'num.vars' not found" (c2); it refuses forest terms and data objects with bases.
- `pdbart(y ~ ., power = 3)` warns "'power' is BayesTree's spelling; 'pdbart' now fits through 'bart', which takes
  'tree.prior = cgm(power = )'" ([`pdbartBayesTreeNames`](../../R/tombstones.R),
  [`translatePdbartCall`](../../R/tombstones.R)), and rewrites into `tree.prior`; neither reached a release.
- `bartBT()` takes `power`, `base`, `splitprobs = 1/numvars` and always builds a `cgm()` call with all three;
  `rbart_vi()` (deprecated, removed in 1.1-0) takes `power`, `base`, builds `cgm(power, base)`, and refuses
  `split.probs` and `tree.prior`.
- The live sampler. `$setModel` with `m@tree.prior@power <- 1` is honored on a one-forest sampler, refused on a
  DART one and on one of several forests; `$setLeafPrior(forests = list(forest(power = 1)))` is "'power' is fixed
  at creation" (q8). No `$setTreePrior` exists, and no `$setForest`.
- Consumers. bartCause's `bcf()` builds `quote(cgm(power, base))` for `dbarts(tree.prior = )` and states
  `treatment.base`, `treatment.power` on the treatment forest; its `translateRetiredBartArguments` turns `bartc`'s
  `power`, `base`, `split.probs` into `tree.prior = cgm(...)` for `bart()` and, under `crossvalidate = TRUE`, for
  the `xbart()` call its bayesOpt.R builds from the same call, where the defect drops them; `retiredSpellings`
  lists them as `tree.prior`'s; man/bartc.Rd documents that; test-17, test-21 read `model@tree.prior@power`;
  test-03 and test-14 pass `dbartsPriors$cgm(2.0, 0.95)` to `dbarts()`. Released bartCause 1.0-10 forwarded to
  0.9-34's `bart2`, which refused `tree.prior`, so `bartc(tree.prior = )` never shipped. stan4bart forwards every
  `bart_args` name that is a `dbartsSpec()` formal, then rebuilds `power`, `base`, `split.probs` into `cgm()` and
  refuses them beside `bart_args$tree.prior`; its test greps the verbose "power and base for tree prior" line.
  treatSens calls `dbartsSpec()` with one `forest(n.trees = )` and `pdbart()` with neither; nothing to change.
  All three install dbarts from GitHub main in CI.

## Inventory

Every place `base`, `power`, `split.probs`, `tree.prior`, `cgm()` and `dart()` are read, stored, printed or
documented, by grep of R/, man/, src/R_interface*.cpp, inst/, benchmarks/R, vignettes, README.md and docs/design,
and what changes. "Respell" means the call is rewritten in the plain spelling with the same values.

R (slice 1):

- R/bart.R. [`bart`](../../R/bart.R) gains `power = 2, base = 0.95, split.probs = NULL` after `k`, where
  0.9-34's `bart2` had them, and loses `tree.prior`; [`buildSamplerPriors`](../../R/bart.R) stops building a
  `cgm()` call and forwards the three plain (`split.probs` unevaluated, as now); the consolidated-name reading of
  the three goes. The `family = "multinomial"` refusal tests the resolved value (D6). [`bart2Hurdle`](../../R/bart.R)
  hands the three to both components where it handed `tree.prior`. [`bartBT`](../../R/bart.R) forwards `power`
  and `base` plain and `split.probs` only where its caller named `splitprobs`. The stated record
  (`plainStatedAttr`) names each of the three separately.
- R/tombstones.R. `power`, `base`, `split.probs` leave `consolidatedArgReasons`, `consolidatedPriorScalars`,
  `consolidatedArgsFor$bart` and the registry rows owned by `bart`. A `tree.prior` row for `dbarts` and
  `dbartsSpec` joins, with its warning (D3). The refusal of `tree.prior` at `bart()` and `xbart()` names the
  plain spelling (D4). `power` and `base` leave [`pdbartBayesTreeNames`](../../R/tombstones.R) and its spelling
  table: they are `bart()`'s own names, passed through with no warning; `splitprobs` maps to `split.probs` and
  still warns. The rewrite into `cgm()` in [`translatePdbartCall`](../../R/tombstones.R) goes.
- R/model.R. [`parsePriors`](../../R/model.R) builds the record from `power`, `base` and `split.probs` (D1, D2).
  [`dart`](../../R/model.R) loses `power` and `base`. [`cgm`](../../R/model.R) stays in
  [`dbartsPriors`](../../R/model.R) as `tree.prior`'s retired value until 1.1-0. `plainForestArgNames` and
  [`plainForestStated`](../../R/model.R) record `power` and `base` by name; `split.probs` leaves the stated set
  (D6). The collision check (D5) is one helper beside [`refuseStatedWithNoPlainForest`](../../R/model.R), whose
  remedy names the argument given. [`forest`](../../R/model.R)'s and [`varianceForest`](../../R/model.R)'s
  comments, and the `dbartsVarianceForest` format's "<mean forest's>", say what D5 says.
- R/formulaTerms.R. [`refuseGivenTwiceInFormula`](../../R/formulaTerms.R) becomes the formula half of D5 for any
  of `n.trees`, `power`, `base`.
- R/spec.R. [`resolveSamplerSpec`](../../R/spec.R): the fitting function's values come from the plain
  arguments; the first forest's knobs are written over the record only where the fitting function's did not
  reach it (D5 refuses the rest); the unstated variance forest reads the resolved record, as now (D5). The
  several-forest and multinomial refusals keep main's test on the resolved value (D6).
  [`dbartsSpec`](../../R/spec.R) gains the three plain arguments and keeps `tree.prior = NULL` as a trailing,
  warned formal (D3).
- R/dbarts.R. [`dbarts`](../../R/dbarts.R) the same, `tree.prior = NULL` trailing beside `sigma` and
  `node.prior`. `$setModel`'s DART refusal text (replaced on setmodel-parameters' branch) and every other read of
  `model@tree.prior` stay.
- R/xbart.R. [`xbart`](../../R/xbart.R) loses `tree.prior` (D4); `split.probs` takes `dart()` and `num.vars`,
  resolved in the split-probability vocabulary; [`cellModel`](../../R/xbart.R) writes each cell's grid values onto
  the record of either class, as now. The comment defending the grid's override of an object's power goes.
- R/rbart.R. [`rbart_vi`](../../R/rbart.R) forwards `power` and `base` plain; nothing else (deprecated).
- R/A_class.R, R/bartcore.R ([`setModelRefused`](../../R/bartcore.R)), R/forestBasis.R (`cgm` and `dart` among the
  names reserved in a basis term): unchanged.
- src/R_interface_bartcore.cpp ([`bartcore_setModel`](../../src/R_interface_bartcore.cpp), the creation reads,
  the verbose line) and inst/include/dbarts/dbarts.h: unchanged. tests/cpp: one comment names `cgm`; unchanged.

Help, NEWS, docs (slice 2):

- man/bart.Rd: usage, `power`, `base`, `split.probs` items (a vector, `NULL` uniform, `num.vars`, or `dart()`),
  `tree.prior` gone from the `tree.prior, leaf.prior` item, the `cgm(3, 0.25)` of the forest paragraph written
  `base = 0.25, power = 3`, examples. man/dbarts.Rd and man/dbartsSpec.Rd: the same three, the same
  `cgm(3, 0.25)`; `tree.prior` gets its own item as a deprecated alias, as `sigma` and `node.prior` have, not a
  live one. man/xbart.Rd: `tree.prior` item gone, `split.probs` takes `dart()`, the `...` item's "dart, now
  tree.prior = dart()" respelled. man/dbartsPriors.Rd: `dart(a, b, rho, alpha, update.alpha, update.delay)`,
  given as `split.probs`; `cgm` marked as `tree.prior`'s retired spelling; examples. man/dbartsFamilies.Rd: its
  "inside tree.prior" sentence names `split.probs`. man/rbart.Rd: the `power`, `base` item stops saying `bart` no
  longer takes them. man/forest.Rd, man/varianceForest.Rd: the fitting function's `power` and `base`, the
  governing paragraph rewritten for D5. man/monotone.Rd: the mBART prior as `power = 0.8, base = 0.25`.
  man/dbarts-deprecated.Rd: the `power`, `base`, `split.probs` items go; a `tree.prior` item on `dbarts` and
  `dbartsSpec`. man/pdbart.Rd, man/dbarts-package.Rd, man/dbartsSampler-class.Rd (the `setModel` DART sentence):
  one line each. man/bartBT.Rd: unchanged.
- inst/NEWS.Rd: [NEWS](#news).
- vignettes/gibbs_sampler_mixture_model.Rmd (`tree.prior = dart()`) and vignettes/dbarts-as-a-component.Rmd
  (`tree.prior = dbartsPriors$cgm(power = 0.5)`), README.md: respell.
- docs/design: [3a. Prior specification](../design/public-surface.md#3a-prior-specification) and
  [4. DART exposure](../design/public-surface.md#4-dart-exposure) gain a dated revision paragraph;
  [The rule](../design/forest-defaults-by-kind.md#the-rule) and
  [What counts as stated](../design/forest-defaults-by-kind.md#what-counts-as-stated) are rewritten for the
  plain names, D5 and D6; monotone.md, prior-defaults.md, consumer-spec-surface.md, pdbart-on-bart.md and
  bart-landscape.md respell sentences stating present behavior. The rest (nog-gibbs.md, core-generalization.md,
  level-fiber.md, robust-errors.md, written-surface.md, interaction-constraints.md, multinomial.md,
  feature-matrix.md) are records or name the prior in passing: a sentence stating present behavior is respelled,
  a dated record is not.

Tests and benchmarks:

- inst/tinytest: about 200 call sites in 39 files state `tree.prior = cgm(...)`, `dart(...)` or the retired
  `power`/`base`/`split.probs` on `bart()`; each is respelled, and the pins of the retired spellings' warnings
  (test-tombstones.R, test-front-door-control.R, test-argument-surface.R, test-pdbart.R, test-model-priors.R,
  test-forest-defaults.R, test-xbart-model.R) are inverted or replaced. test-forest-defaults.R's uniform
  `split.probs` refusal beside forests that all have a basis becomes an acceptance (D6). The four snapshot files
  state none.
- benchmarks/R: 28 scripts, about 60 sites, respelled; equivalence.R's `dart` scenario,
  `tree.prior = dbarts:::dart(2, 0.95)`, becomes `split.probs = dbarts:::dart()` (positional `dart(2, 0.95)`
  would now be `a = 2, b = 0.95`). bcf-equivalence.R and multinomial-equivalence.R state none.

## Design

D1. The record does not change. A model's `tree.prior` slot keeps holding a `dbartsCGMPrior` (power, base, the
normalized `splitProbabilities`) or a `dbartsDartPrior` (those and DART's settings), the first forest's, as now.
Every door resolves its plain arguments into that record: `power` and `base` into the slots, a vector through
[`resolveSplitProbabilities`](../../R/model.R) as now, `dart()` as the class with its settings and the resolved
`update.delay`. So `$setModel` on model objects (`m@tree.prior@power <- 1`), `copy()`, `storeState` and a reload,
setmodel-parameters' slices 1 to 3 (which read DART's settings and the class off the record), the bridge, and
bartCause's tests reading `model@tree.prior@power` keep working, and the record main writes for a call is the
record this item writes for its respelling, slot for slot. Getters: `s$model` shows the same slots; the verbose
summary's "power and base for tree prior" line stays; `print` and `summary` show no tree prior, as now. No reader
is added here (a forest's reader is `$getForest(forest)`, built by the forest-door item).

D2. `dart()` takes `a`, `b`, `rho`, `alpha`, `update.alpha`, `update.delay` and returns the record class with
`power` and `base` at 2 and 0.95, which the door overwrites with its own. `dart(power = 3)` is R's
unused-argument error, as `forest(leaf.prior = )` is under dec-B508; `dart` is new in 1.0-0. A bare `dart` is
`dart()`, as `tree.prior = cgm` meant `cgm()`. `split.probs` resolves in the vocabulary `dart`, `num.vars`,
`numvars`, layered over the caller's frame, at every door, `xbart()` included, and a forwarded `..N` is recovered
as now. A `cgm()` object or call given as `split.probs` is refused by name: "'split.probs' takes a vector of
split probabilities or dart(); cgm() is the retired tree.prior's value: write power, base and split.probs".

D3. `dbarts()` and `dbartsSpec()` take `power = 2`, `base = 0.95`, `split.probs = NULL` where `tree.prior` stood
(after `n.samples`, and after `control`), and keep `tree.prior = NULL` as a trailing formal, as `sigma` and
`node.prior` are kept, so a caller keyed on the formals (stan4bart's `bart_args`) still reaches it. `NULL` states
nothing and does not warn. Any other `tree.prior` warns once a session, class `dbartsDeprecatedWarning`:
"'tree.prior' has left 'dbarts': write the tree prior as power, base and split.probs (split.probs = dart() for
DART), as bart() takes them; 'tree.prior' is removed in dbarts 1.1-0. The value was used." It is read as the plain
arguments it states, in four cases:

1. A written call whose head is `cgm` or `dart` in any spelling: bare, `dbartsPriors$cgm`,
   `dbarts::dbartsPriors$cgm`, `dbarts:::cgm`, and the same for `dart`. It is matched against the constructor's
   signature on main (`cgm(power, base, split.probs)`, `dart(power, base, a, b, rho, alpha, update.alpha,
   update.delay)`) and states the arguments it names: `cgm(power = 3)` states `power` alone, `cgm(3, 0.5)` power
   and base; `dart(power = 3, a = 1)` is `power = 3, split.probs = dart(a = 1)`; a `split.probs` it names is
   stated as written, so a uniform one restates the default (D6).
2. A bare name: `cgm` states nothing; `dart` is `split.probs = dart()`.
3. A value (any other expression), an object built by the constructor or a record read off a sampler: it states
   `power` and `base`, and `split.probs` only where it is DART or its probabilities, written
   (`splitProbabilitiesSpec`) or resolved (`splitProbabilities`), are not uniform; a DART record's settings are
   taken as they stand, its resolved `update.delay` included, as main keeps them (c3b).
4. A `cgm` object saved by 0.9-34, which has no `splitProbabilitiesSpec` slot: read as case 3 from `power`,
   `base` and `splitProbabilities`, where main fails "no slot of name".

`tree.prior` beside any of `power`, `base` or `split.probs` is refused: "'tree.prior' and 'power' both set the
tree prior of 'dbarts'; give 'power' alone". A tree-prior object given positionally as the tenth argument of
`dbarts()`, where 0.9-34's `tree.prior` stood and `power` now stands, is read as `tree.prior`, warned the same.

D4. `bart()` and `xbart()` take no `tree.prior`, which never shipped on either. Given one, through `...`, each
refuses it by name with the spelling D3's four cases translate it to: "'tree.prior' is not an argument of 'bart':
write the tree prior as plain arguments, power = 3, split.probs = dart(a = 1)". On `xbart()` the message adds
that `power` and `base` are grid axes. So `xbart(tree.prior = cgm(power = 3))` can no longer fit the default
depth; `xbart(power = 3)` fits power 3, as on main.

D5. Two levels, one rule. The fitting function's `power` and `base` reach every mean forest with no basis that
does not state its own (dec-B512); a forest stating its own keeps it. A variance forest that states none reads
the resolved record of the mean forest, as main: the fitting function's values, or a lone forest term's own where
it states them (c2); a variance forest stands only beside one forest, so the two agree wherever both are stated
and D5's refusal applies. A name given to the fitting function that reaches no mean forest is refused by name,
before anything is built:

- where every mean forest has a basis (dec-B274, now per name): "'power' given to the fitting function is the
  tree prior of the forest with no basis, and every forest of this model has a basis; state it on a forest, as
  forest(x1, basis = a, power = 3)", the remedy naming the argument given;
- where every forest with no basis states it (dec-B512, as for `n.trees`): "'power' is given to the fitting
  function and to the forest with no basis ('forest(x1 + x2, power = 1)'), which are the same setting; give
  one" (the code where a formula wrote it, "(forest 1)" in a list or a data object's `bases`).

Stated means named, the default's own value included, as for `n.trees`: `bart(power = 2)` beside
`forest(x, power = 1)` is refused, and `bart()` forwards which of the two its caller named. Names are checked one
by one: `power` on the fitting function and `base` on the forest are both taken. A variance forest never rescues
a value: `bart(y ~ forest(x, power = 1), power = 3, variance = ~ x3)` is refused, the lone mean forest stating
`power` twice. The helper takes the positions of the forests with no basis, so forests-shared-total, which makes
them several, only widens `plainForest`'s answer; its own `n.trees` rule (refused only where every sharer states
one) is this rule. A lone forest written as a term, `y ~ forest(x, power = 1)`, is the model's forest, and
`power = 1` there is the model's tree prior, as `sd` there is its leaf prior (dec-B266).

D6. `split.probs` is the mean forest's alone, as on main: the variance forest never takes split probabilities,
and `forest()` and `varianceForest()` take no `split.probs`. Where no forest can take split probabilities (a model
of several forests, `family = "multinomial"`) main's test stands: a resolved non-uniform vector or `dart()` is
refused, by its own message, "a model of several forests does not support 'split.probs' (a vector that is not
uniform, or dart()); drop it or fit a single-forest model"; a uniform one (`c(1, 1, 1)`, `1/num.vars`,
`bartBT()`'s default) states the default and is taken. `split.probs` is therefore not in D5's naming rule: beside
forests that all have a basis a uniform one, which main refuses by naming, is taken, and a non-uniform one meets
D6's refusal, not D5's remedy, which no forest could follow. `bartBT()` forwards `splitprobs` only where its
caller named it.

D7. The live door: setmodel-parameters' slice 4 ([The live door](#the-live-door)).

D8. `cgm` stays in `dbartsPriors`, documented as the value of the retired `tree.prior`, and leaves with it at
1.1-0, its class staying as the record's. A TODO entry with the 1.1-0 removals records it.

### The live door

Slice 4 of setmodel-parameters (on its branch, planned at a628e48e) is to give a live sampler the tree priors of
every forest, through `$setTreePrior(tree.prior = , forests = , variance = )`, its open call 1 (a). Under dec-B512
and dec-B518 the door becomes:

    $setTreePrior(power = NULL, base = NULL, split.probs = NULL)

- One level: `power` and `base` are the fitting function's, and `$setTreePrior` takes no `forests` or `variance`.
  A forest's own depth changes through `$setForest(forest, power = , base = )`, by label or position, and the
  variance forest's through `$setVarianceForest(power = , base = )`; both belong to the forest-door item, whose
  name, labels for a multinomial sampler's category forests and spelling for returning a forest to the fitting
  function's level are open with the maintainer. Statedness is per forest and per name, as creation recorded it
  (a lone forest term stating `power` states it at the forest's level, so slice 4 records it on the one-forest
  path as well) and as `$setForest` sets it from then on.
- `power`, `base`: each reaches every forest with no basis that does not state its own; NULL leaves it as it is.
  Where, for a name, every forest with no basis stated its own, the call is refused, because creation refuses the
  same statement (dec-B512's collision rule; dec-B515: a setter refuses only where creation refuses the same
  statement). The check is per name, and the refusal names the forests that stated it and the call to write:
  "'power' given to $setTreePrior reaches no forest: forest 1, the forest with no basis, states its own; change
  it with $setForest(1, power = 3)". The same rule covers `$setModel` with a changed `power` or `base` on a
  sampler of several forests, where slice 4 makes its tree prior the fitting function's level, and
  `$setLeafPrior(normal(sd = ))` where every forest with no basis stated its own `sd` (dec-B515's writer, its own
  item). On a sampler of one forest `$setModel` keeps main's behavior: its model's tree prior is that forest's
  record (D1), which it replaces whole.
- `split.probs`: a vector changes the values (a parameter, dec-B254); `dart(...)` on a DART sampler changes
  DART's settings, a setting it leaves out taking creation's default for it, as creation with `dart(...)` would
  (dec-B514's principle; slice 3's install); `1` restores equal probabilities, as creation with none has; a
  non-uniform vector on a DART sampler or `dart()` on an ordinary one is structure and is refused: "split.probs =
  dart() makes a DART sampler of an ordinary one, which changes the model's structure; make a new sampler
  instead".
- `$setForest` takes `sd`, `amplitude` and the move settings beside `power` and `base`, refuses structure
  (vars, n.trees, basis, interactions, blocks) by name, and refuses a statement where creation refuses the
  same; it is not built here.
- No `tree.prior` argument: it never shipped on a live sampler. `$setModel` keeps taking model objects
  (dec-B512), and slice 4's install from the record (its item 3) is unchanged.

Whether a variance forest that stated no tree prior follows a live change of the fitting function's (slice 4's
open call 2, held decision #6) is pending with the maintainer; slice 4 waits on it. This item does not decide it.

This item builds none of the door: slice 4 builds `$setTreePrior` after forests-shared-total, and the forest-door
item builds `$setForest` and `$setVarianceForest`. The amendment to the setmodel-parameters plan (its door, item
4, its open call 1, and the respelling below) is a records step for that branch's writer, before slice 4 is
built.

The setmodel-parameters branch's slices 1 to 3 add about 36 lines in the old spellings, each respelled when the
branch rebases onto this item (counted at a628e48e): inst/tinytest/test-setmodel-parameters.R 13 (`tree.prior =
dart(...)` and `dart(power = 3, base = 0.8)`, the last as `power = 3, base = 0.8, split.probs = dart()`),
inst/tinytest/test-split-probs-zero.R 8 and benchmarks/R/split-probs-exact.R 3 (`cgm(split.probs = )`),
inst/tinytest/test-setmodel-structure.R 4, inst/tinytest/test-single-forest-vars.R 1, inst/NEWS.Rd 1,
man/dbarts.Rd 1 (its `tree.prior` item, merged with this item's), and R/bartcore.R 1. The last is slice 1's
structural refusal, which names "the tree prior's kind" with the values "dart()" and "cgm()"; it is respelled to
name the setting the caller now writes: "'split.probs' differs in kind: the sampler's is dart(), the model's a
vector; make a new sampler instead", the pins in test-setmodel-structure.R following it.

## Collisions

Every pair of levels at which one setting can be stated twice, and what happens.

| levels | example | result |
|---|---|---|
| fitting function, forest with no basis | `bart(y ~ forest(x, power = 1) + forest(x, basis = z), power = 3)` | refused, D5 |
| the same, other names | `bart(... forest(x, base = 0.5) ..., power = 3)` | power 3, base 0.5 |
| fitting function, forest with a basis | `bart(... forest(x, basis = z, power = 4), power = 3)` | taken: 3 on the forest with no basis, 4 on the other |
| fitting function, no forest with no basis | `dbarts(y ~ forest(x1, basis = a) + forest(x2, basis = b), power = 3)` | refused, D5 (dec-B274) |
| fitting function, variance forest | `bart(y ~ x, power = 3, variance = varianceForest(x3, power = 1))` | taken: 3 and 1 |
| lone forest term, unstated variance forest | `bart(y ~ forest(x, power = 1), variance = ~ x3)` | taken: 1 on both, as main |
| fitting function, lone forest term | `bart(y ~ forest(x, power = 1), power = 3)` | refused, D5 |
| uniform `split.probs`, several forests | `bart(y ~ forest(x1) + forest(x2, basis = z), split.probs = 1/num.vars)` | taken, as main (D6) |
| `tree.prior` (retired) and the plain names | `dbarts(..., tree.prior = cgm(), power = 3)` | refused, D3 |
| `tree.prior` (retired), forest | `dbarts(y ~ forest(x, power = 1) + ..., tree.prior = cgm(power = 3))` | refused, D5, after D3's warning |
| `tree.prior` object, forest, other name | `tp <- cgm(power = 3); dbarts(... forest(x, base = 0.5) ..., tree.prior = tp)` | refused: a value states base too |
| `tree.prior` at `bart()`, `xbart()` | `bart(y ~ x, tree.prior = cgm(power = 3))` | refused, D4, naming `power = 3` |
| pdbart's BayesTree name and bart's | `pdbart(y ~ ., splitprobs = p, split.probs = p)` | refused, as now ("both set one setting") |
| bcf()'s and bart()'s, bartCause | `bcf(..., power =, treatment.power =)` | two forests; no collision |
| stan4bart `bart_args` | `bart_args = list(power = 2.5, tree.prior = cgm())` | refused by dbartsSpec, D3 |
| live door, forest stated at creation | `$setTreePrior(power = 3)`, every forest with no basis having stated power | refused naming those forests and `$setForest` |
| live door, one forest | `$setForest(1, power = 2)` | changes that forest only; the fitting function's level is `$setTreePrior`'s |
| live door and `$setModel` | separate calls | each installs; the later stands |

A value no forest can take, each refused by name: the fitting function's `power` or `base` where no mean forest
has no basis or every one states it (D5); a non-uniform `split.probs` or `dart()` in a model of several forests or
under multinomial (D6); `cgm()` as `split.probs` (D2); `split.probs = dart()` swapped for a non-uniform vector on
a live sampler or the reverse (structure); `power` not a single positive number and `base` not a single number in
(0, 1), at every door with the messages `forest()` gives ("'power' must be a single positive number", "'base' must
be a single number in (0, 1)"), `xbart()` keeping its grid's "must contain only" wording. `forest(split.probs = )`,
`varianceForest(split.probs = )`, `forest(tree.prior = )` and `dart(power = )` are R's unused-argument error.

## NEWS

Scope is what differs from 0.9-34; nothing that never reached a release.

- The tombstone item loses "On bart, power, base and split.probs become tree.prior = cgm(power, base,
  split.probs)": `bart()` takes them as 0.9-34's `bart2` did. It gains: on `dbarts`, `tree.prior` becomes
  `power`, `base` and `split.probs`, as on `bart`; `tree.prior = cgm(power, base, split.probs)`, or a tree prior
  saved by 0.9-34, is accepted with a once-per-session warning, the value used, and is refused from 1.1-0. A
  0.9-34 `dbarts(tree.prior = cgm(power = 3))` caller sees that warning once and fits the same model.
- The "Tree priors" item becomes: `split.probs = dart()` runs DART (Linero 2018) on `bart`, `dbarts`,
  `dbartsSpec` and `xbart`; "`bart` gains `tree.prior`" goes; `leaf.prior` stays. The `xbart` item says
  `split.probs` is new there, a vector or `dart()`, and names no `tree.prior`.
- Nothing about `pdbart`'s `power` (its warning never reached a release), the xbart defect, or the two-level
  refusal: none reached a release (0.9-34 had no forests and no `tree.prior` on `xbart`).
- bartCause's own NEWS (slice 3): `bartc(crossvalidate = TRUE)` given `power`, `base` or `split.probs` now
  cross-validates under them, where the development builds dropped them; its chosen k and draws change.

## Order of landing

1. This item, slices 1 and 2 on branch tree-prior-flat, their reviews and the gates; then bartCause's and
   stan4bart's branches (slice 3), built and green against this branch's tip installed from the branch; then
   dbarts's merge, followed in the same sitting by the two consumer merges, and their CI re-run on main. Between
   the dbarts merge and each consumer's: `bartc(power = )` is refused (D4, `bart()` given `tree.prior`); stan4bart's
   `bart_args$power`, `base` or `split.probs` is refused (D3, its formals loop now forwarding them beside the
   `tree.prior` it still builds), while `bart_args$tree.prior` keeps working with the warning. The consumers'
   branches are ready first so the window is minutes. A compatibility shim in the consumers (reading
   `formals(dbarts::bart)`) would remove it at the cost of dead code; not taken (Calls made 14).
2. coefficient-prior-opt-in (being built on its branch): independent; textual conflicts in R/model.R
   (`dbartsPriors` gains `student`, this item edits `dart` and `cgm` there), man/dbartsPriors.Rd, man/forest.Rd,
   inst/NEWS.Rd and tinytest files. Whichever is green first lands first; the second rebases. Pinned by whichever
   lands second: `split.probs = student(3)` is refused as not a split-probability value, and `amplitude = dart()`
   as not a coefficient prior.
3. forests-shared-total: this item lands before its slice 2 (the R surface) rebases. Its D7 and R item 2 (each
   sharer takes the fitting function's base and power where it states none) then read the plain arguments, and
   D5's helper takes its vector of plain positions. Its slice B (built) and slice 1 (engine) do not meet this
   item. Shared: [`resolveSamplerSpec`](../../R/spec.R), [`forestParams`](../../R/model.R),
   [`plainForestStated`](../../R/model.R), test-forest-defaults.R.
4. setmodel-parameters: slices 1 to 3 meet this item only textually (the respelling listed in
   [The live door](#the-live-door); slice 3's DART install reads the unchanged record). This item lands before
   slice 4, whose door is amended as that section says; the branch, which merges once after slice 4, rebases onto
   it. Slice 4 also waits on its open call 2.
5. control-migration (ruled, unplanned): it puts `n.trees`, `n.cuts`, `useQuantiles` on `dbarts()` beside this
   item's three; D5's helper covers `dbarts(n.trees = )` beside a stated forest count. Whichever is second
   rebases; this item is expected first.

## Slices

### Slice 1: the doors and the record (R)

Steps: the R inventory above, D1 to D6 and D8; the tinytest re-pins and the new test file.

1. [`parsePriors`](../../R/model.R) takes `power`, `base`, `split.probs` (and the retired `tree.prior`, read by
   D3's four cases first) and builds the record; `dart()` per D2.
2. `dbarts()`, `dbartsSpec()`: formals, the trailing `tree.prior` and its warning, the positional tenth argument,
   the both-given refusal.
3. `bart()`, `bartBT()` (its `splitprobs` forwarded only where named), `bart2Hurdle()`, the multinomial path,
   `rbart_vi()`: forward the plain three; the `bart()` tombstones for them go; D4's refusal.
4. `xbart()`: D4, `split.probs = dart()` and `num.vars`.
5. D5: the stated record per name, the helper for formula, list and data-object doors, the per-name refusal
   where no forest has no basis; D6's test on the resolved value.
6. pdbart's table.
7. Tests ([Tests and gates](#tests-and-gates)).

RNG: neutral for every call main accepts; refusals added, two refusals lifted (D3 case 4, D6). Cost: ~1140
lines; 8 agent hours.

### Slice 2: help, NEWS, docs, benchmark scripts

Steps: the help, NEWS and docs inventory; respell the 28 benchmark scripts and the scenario in equivalence.R; the
records' TODO entry for 1.1-0 (D8). Cost: ~440 lines; 3 agent hours.

### Slice 3: consumers, each on its own branch

- bartCause (~150 lines, 2.5 agent hours): `bcf()` sets `samplerCall$power` and `$base` instead of `tree.prior`;
  `translateRetiredBartArguments` stops building `cgm()`, so `power`, `base`, `split.probs` pass through to
  `bart()` and to the cross-validation's `xbart()` call; `retiredSpellings` and `dropOverriddenSpellings` lose the
  `tree.prior` entry; `bartc(tree.prior = )`, which never shipped (1.0-10 forwarded to 0.9-34's `bart2`, which
  refused it), is refused by bartCause with the translated spelling, D4's text under `bartc`'s name. Tests:
  test-17 and test-21 keep reading `model@tree.prior@power` (D1) and lose the both-given refusal of `tree.prior`;
  test-03 and test-14 respell `dbartsPriors$cgm(2.0, 0.95)`; a new pin that `bartc(crossvalidate = TRUE, power =
  1.5)` cross-validates at power 1.5. man/bartc.Rd: `power`, `base`, `split.probs` are passed through, not
  collected into `tree.prior`. NEWS: the crossvalidate line ([NEWS](#news)).
- stan4bart (~50 lines, 1.5 agent hours): the `cgm()` rebuild and its refusal go, the formals loop forwarding
  `power`, `base`, `split.probs` as `dbartsSpec()` formals; `bart_args$tree.prior` reaches `dbartsSpec()`'s kept
  formal, which warns (D3); test-09-bartArgs.R's "both 'tree.prior' and 'power'" pin takes dbarts's message, and
  its `tree.prior = cgm(power = 2.5, base = 0.9)` arm counts the one warning. The verbose line is unchanged.
- treatSens: nothing (searched: one `dbartsSpec(forests = list(forest(n.trees = )))`, `pdbart()` with neither).

Then: the opus review of slice 1, 2.5 hours; the sonnet review of slices 2 and 3, 1.5; the gates, 3 hours (the
A/B grows); the fix rounds, 2; records, 0.5. About 24 agent hours; plans here run 1.5 to 2 times low, so budget
36 to 48.

## Tests and gates

tinytest, a new test-tree-prior-flat.R, small n, two chains where draws are compared:

- Every door takes the plain three and writes main's record: `bart()`, `dbarts()`, `dbartsSpec()`, `xbart()`
  (per cell, through the fit's model, CGM and DART records both), `bartBT()`, `rbart_vi()`, `pdbart()`, the hurdle
  family's two components, monotone; `model@tree.prior` `identical()` to the record a constructed `cgm()` or
  `dart()` with the same values gives after resolution.
- `split.probs`: a named vector, `c(2, 1, 1) / num.vars` (at `xbart()` too), `dart()`, bare `dart`,
  `dbartsPriors$dart(a = 1)`, through a wrapper's `...`; DART's `update.delay` default half the burn-in;
  `dart(power = 3)` an unused-argument error; `dart(0.3)` is `a = 0.3`; `cgm()` as `split.probs` refused.
- D3, each case: written `cgm(power = 3)`, `cgm(3, 0.5)`, `dbartsPriors$cgm(power = 3)`, `dbarts:::cgm(...)`,
  `dart(power = 3, a = 1)`; bare `cgm` and `dart`; a constructed object; a record read off a sampler (CGM with
  non-uniform probabilities; DART with its resolved `update.delay` under another burn-in); a 0.9-34-saved object
  (a fixture saved by 0.9-34, as scratch/tpf-crit/c3a.R makes it, read with `power`, `base` and its
  probabilities); `tree.prior = NULL` silent; the positional tenth argument; each warning exactly once a session
  (counted with `withCallingHandlers` over two calls, after resetting the once-per-session key); `tree.prior`
  beside `power` refused.
- D4: `bart(tree.prior = cgm(power = 3, base = 0.5))` and `xbart(tree.prior = dart(power = 3, a = 1))` refused,
  each message naming the translation.
- D5, by door (formula, `forests` list, a data object's `bases`) and by name: the refusal where the forest with no
  basis states it, the default's own value included (`bart(power = 2)`); other names taken; a forest with a basis
  taken; no forest with no basis, each name refused with its own remedy; the variance-forest case refused; a lone
  forest term alone taken as the model's.
- D6: beside two forests and under multinomial, a non-uniform vector and `dart()` refused with the new message;
  `c(1, 1, 1)`, `1/num.vars` taken; beside forests that all have a basis a uniform one taken (the moved pin);
  `bartBT()` on a data object with bases at its default `splitprobs` taken, and with `splitprobs = c(1, 2, 1)`
  refused.
- The variance forest: unstated takes the fitting function's `power` and `base` (its record 3 and 0.5) and a lone
  forest term's own (1 and 0.5, c2's call), stated keeps its own, `split.probs` never reaches it.
- pdbart: `power`, `base` pass with no warning; `splitprobs` still warns and maps.
- Continuation: a DART sampler made with `split.probs = dart()`, run, `copy()`, `storeState` with `saveRDS` and
  `readRDS`: next draws `identical()`; `$setModel` with `m@tree.prior@power <- 1` on a one-forest sampler taken,
  as on main.
- Re-pins: the 39 files respelled; every pin of a retired-spelling warning on `bart()` inverted to no warning.

Gates, on the branch's own library:

- `R CMD INSTALL -l <lib>`; `buildInfo()$mode` checked.
- A neutrality A/B (scratch, ~400 lines): an 83cefe7d library against the branch's, two chains, every recorded
  channel compared with `identical()`: each door in each spelling main accepts against its plain respelling on the
  branch, for every family path [`buildSamplerPriors`](../../R/bart.R) feeds (gaussian, student, probit,
  logistic, aft, ordinal, nbinom, multinomial, hazard and its probit and logistic forms, hurdle), each with a
  stated `power` and `base` and, where the family takes them, `split.probs` and `dart()`; `xbart`'s grid with and
  without DART; `bartBT` (its default `splitprobs`, a stated one, a data object with bases), `rbart_vi`, `pdbart`,
  monotone; heteroscedastic with and without a stated variance tree prior and with a lone forest term stating
  `power` and `base` (H2's call); bcf-shaped forests with the forest with no basis first and second, with a
  uniform `split.probs` and with a `cgm` object (bartCause's test shapes); `split.probs = 1/num.vars` beside two
  forests; the default call with no tree prior stated, which must not move. A difference is a stop.
- The equivalence trio on the reference build, `compare --bitwise --strict-coverage`: 55 of 55
  (equivalence-32418b6e, its `dart` scenario respelled), 19 of 19 (bcf-equivalence-20eed7ed), 11 of 11
  (multinomial-equivalence-80b1c8d4), identical, no `max |z|` line; the four snapshot files unchanged.
- The full tinytest suite in one process; `lintr::lint_package()` (names move), `air format --check .`,
  `tools/check-rc-codoc.R`, `tools/check-win-drift.R`, `tools/check-doc-freshness.R`, each on its own exit
  status; `inst/NEWS.Rd` parses; `R CMD check --as-cran` from a tarball of a clean copy staged outside the tree.
- exact-gates.yaml's list in quick mode once on the branch tip (its scripts are respelled; no number may move).
  CI's exact-gates runs it too: benchmarks are not ignored there.
- Consumers: bartCause's and stan4bart's suites against the branch's library, each branch green.
- No sanitizer run beyond CI's (no compiled code moves); no `bench-sampler.R` compare.

Mutants, each to fail a named expectation:

- [`parsePriors`](../../R/model.R) drops a stated `base` for the default (the record pins; the A/B).
- D5 check removed, the forest's value winning (the two-level refusal pin); D5 checking names jointly, refusing
  `power` beside a forest's `base` (the other-names pin); `bart()` not forwarding which names its caller stated,
  so `bart(power = 2)` beside `forest(power = 1)` is taken (the default-restated pin); a missing `power` read as
  stated, refusing `bart(y ~ forest(x, power = 1))` (the lone-term pin); the variance forest counted as reached
  (the variance-case pin).
- The unstated variance forest reading the fitting function's arguments instead of the resolved record (the
  lone-term variance pin; the A/B's H2 arm).
- D6 testing by naming instead of the resolved value (the uniform acceptance pins; the A/B's bcf-shaped `cgm`
  arm); `bartBT()` forwarding its default `splitprobs` (the bartBT data-object pin).
- D3: a value read as stating only its non-default slots, or `split.probs` where uniform (the value pins); a
  written `cgm(power = 3)` read as stating `base` too (the written cross-name pin); a `dbarts:::cgm` head not
  recognized (the head pin); `dart(power = 3, a = 1)` translated with `power` inside `dart()` (the translation
  pin); a 0.9-34 object read through the missing slot (the fixture pin); the warning every call, none, or on
  `NULL` (the count pins).
- D6's refusal after D5's, giving "state it on a forest" for `split.probs` (the message pin).
- `split.probs = dart()` reaching the variance forest's record (the variance record pin).
- `dart()` keeping `power` as its first formal (the `dart(0.3)` pin; equivalence's `dart` scenario).
- `split.probs` evaluated in the caller's frame (the `num.vars` pins at `bart()` and `xbart()`; the wrapper pin).
- [`rbart_vi`](../../R/rbart.R) dropping `power` (its record pin); [`cellModel`](../../R/xbart.R) not writing the
  grid's `power` onto a DART record, or dropping DART's settings (the xbart DART per-cell pin).
- pdbart still warning for `power` (the no-warning pin).
- equivalence.R's scenario left as `dart(2, 0.95)` (the bitwise compare: `a = 2`).

## Calls made

1. The record does not change (D1): `model@tree.prior` keeps its classes and slots, so neutrality is structural
   and nothing below R moves. Alternative: plain `power`, `base`, `split.probs` slots on the model, which matches
   the surface but moves the bridge, every reader, setmodel-parameters' slices 1 to 3 and bartCause's tests, for no
   difference a user sees.
2. `dart()` returns the record class with its depth slots at the defaults, overwritten by the door (D2).
   Alternative: a class of its own with no depth slots, ~40 lines and a change at every class test, so that
   `print(dart())` shows no `power`.
3. The formals show the defaults, `power = 2, base = 0.95, split.probs = NULL`, as 0.9-34's `bart2` and main's
   `xbart` do, and "stated" is "named" (`missing()`), as for `n.trees`. Alternative: `NULL` defaults, the house's
   "absent", which hide the numbers from the help's usage line.
4. The retired `tree.prior` is read by D3's four cases: a written call by the names it states, a value by its
   depth and its non-uniform or DART split probabilities, a 0.9-34 object from its three slots. Alternative: all
   slots always, which refuses `tree.prior = cgm(power = 3)` beside a forest's `base` and a uniform object beside
   several forests, calls main fits as meant.
5. `tree.prior` on `bart()` and `xbart()` is refused with the translation, not accepted with a warning (D4):
   dec-B512 retires it with a warning only where it shipped. Alternative: a plain unused-argument error.
6. `split.probs` stays the mean forest's: no `split.probs` on `forest()` or `varianceForest()`, none reaching the
   variance forest (D6), as main. Alternative: `forest(split.probs = )` on a lone forest, ~40 lines, refused in
   every model of several forests, where the engine has none.
7. One rule for two levels (D5): a fitting function's value that reaches no mean forest is refused, which is
   dec-B274's refusal and dec-B512's collision at once, and forests-shared-total's `n.trees` rule for several
   forests. The variance forest does not rescue a value, and an unstated one reads the resolved record, as main.
   Alternative: the variance forest counting as reached, which takes `bart(y ~ forest(x, power = 1), power = 3,
   variance = ~ x3)` with power 1 on the mean forest.
8. `split.probs` is judged by its resolved value, not by naming (D6): a uniform one states the default and is
   taken wherever main takes it, and beside forests that all have a basis too, where main refused it by naming;
   a non-uniform one beside several forests meets its own message, not "state it on a forest".
9. A positional `dart(2, 0.95)` now means `a = 2, b = 0.95`, unguarded: `dart` never shipped, and the one such
   call in the tree (equivalence.R) is respelled.
10. `cgm` stays until 1.1-0 and leaves with `tree.prior` (D8); its class stays as the record.
11. pdbart's `power` and `base` leave the BayesTree table, so `pdbart(..., power = 3)` warns no more.
12. `tree.prior` stays a trailing formal of `dbarts()` and `dbartsSpec()`, warned, documented as a deprecated
    alias (D3), so stan4bart's `bart_args$tree.prior`, keyed on the formals, keeps working; a positional tenth
    argument holding a tree prior is read as it. Alternative: `tree.prior` on `...`, which turns stan4bart's
    `bart_args$tree.prior` into an unknown name.
13. The live door takes no `tree.prior`, only the fitting function's level (dec-B518), and `split.probs = 1`
    restores equal probabilities there. A statement that reaches no forest is refused per name, naming the forests
    that stated it and `$setForest`, by dec-B512's collision rule and dec-B515; the same rule covers `$setModel` on several forests and `$setLeafPrior(normal(sd = ))`
    (the coordinator's call, formerly this plan's open call 1). Alternatives: the statement overwriting the
    forests' own, or taken and changing nothing.
14. Consumers merge right after dbarts, from branches made green against it first. Alternative: a shim in each
    consumer that reads dbarts's formals, which removes a minutes-long window at the cost of dead code.
15. Slice 1 takes an opus implementer and review (the coordinator's call); slices 2 and 3 sonnet, the bitwise A/B
    and equivalence trio being the oracle that no draw moves.
16. `bartc(tree.prior = )` is refused by bartCause with the translated spelling: it never shipped there.

## Open calls

None for the maintainer. The question this plan had held, `$setTreePrior(power = 3)` where every forest with no
basis stated its own, is the coordinator's call 13. setmodel-parameters' open call 2, whether an unstated variance
forest follows a live change of the tree prior, remains that plan's, pending.

## Not verified

- Nothing is built; the sizes are from the inventory, not a prototype.
- That every one of the ~200 tinytest and ~60 benchmark sites respells without a change of value is assumed from
  the grep; the A/B, the equivalence trio and the exact gates are the check.
- 0.9-34's positional `dbarts(..., tree.prior = cgm(3, 0.5, c(0.5, 0.25, 0.25)))` and a saved 0.9-34 prior were
  run on 0.9-34 and main only; their reading under D3 is planned, not run.
- setmodel-parameters' slice 4 and forests-shared-total were read on their branches at a628e48e and 0f16f8d5;
  their calls may move. The control migration has no plan.
- The consumers' changes are from reading their code at bartCause f7edf5b and stan4bart 566a1d9; the bartCause
  crossvalidate change of k is from reading bayesOpt.R, not run.
- On a one-forest sampler `$setModel` keeps main's behavior where the critic's third condition would refuse a
  changed power after a lone forest term stated its own; the live door section says why.

## Landing note

Landed 2026-10-30 (dec-A246) on coefficient-prior-opt-in: the doors of slice 1, the help, NEWS, docs and benchmark
scripts of slice 2, and the consumer branches of slice 3, merged into bartCause's and stan4bart's main in the same
sitting. Neutral for every call main accepts, by the A/B on both builds and the equivalence trio; x86-64 leg and
CI green on the branch tip. The live door stays with setmodel-parameters' slice 4.
