# tree-prior-flat: a tree's depth is plain base and power at every door, and DART a value of split.probs

Status: PLANNED 2026-10-29 on 83cefe7d (dec-B512). One call is open for the maintainer, and it concerns
setmodel-parameters' slice 4 door, not this item ([Open calls](#open-calls)). Rulings: dec-B512, with dec-B116,
dec-B246 and dec-B252 (each revised by it), dec-B241 and dec-B503 (the collision rule it extends), dec-B254,
dec-B266, dec-B274, dec-B508 and dec-B514.

agent: one sonnet implementer (effort medium) for slices 1 and 2; one sonnet review over both (effort medium),
told to refute, that runs the mutants below and checks the fix round. Tier "R surface, no draw change"
([Process by risk](README.md#process-by-risk)): no engine, bridge or C header line moves. Slice 3 is two
consumer branches, one sonnet agent each, reviewed by the same reviewer. No blind critique (not that tier).
rng: NEUTRAL, bit for bit, for every call main accepts with the same values: the model's tree-prior record keeps
its classes and slots, and every door writes the same numbers into it. Calls main accepts and this item refuses:
a setting stated on the fitting function and on the forest with no basis, which main resolves silently in the
forest's favor; `tree.prior` on `bart()` and `xbart()`, where it never shipped; and `xbart(tree.prior = cgm())`,
which on main fits the default depth with no message. No call changes its draws; nothing is re-recorded.
window: for 1.0-0; on main before setmodel-parameters' slice 4 and before forests-shared-total's slice 2
([Order of landing](#order-of-landing)).
budget: ~1450 lines in the tree. Slice 1 ~990 (R ~520, tinytest ~470: a new file ~240, re-pins in about 39 files
~230). Slice 2 ~420 (help ~200, NEWS ~25, vignettes and README ~10, design docs ~80, benchmark scripts ~70,
the scenario respelled in equivalence.R included, records ~35). Slice 3, outside this tree: bartCause ~70,
stan4bart ~45. Out of the tree, ~250 lines of scratch A/B script.

## Summary for the maintainer

You ruled that a tree's depth is written as plain `power` and `base` everywhere, that split probabilities are
`split.probs`, a vector or `dart()`, and that `tree.prior` leaves `bart()` and is retired with a warning on
`dbarts()` and `dbartsSpec()`. This plan builds it.

Nothing a user fits changes its draws. The tree prior a sampler holds inside stays the same object it is today,
so `$setModel`, copies, saved fits and the consumers' reads of `model@tree.prior` keep working. What changes is
what the doors take and what they refuse. A 0.9-34 `bart2(power = 3)` call works as it did in 0.9-34, with no
warning; a 0.9-34 `dbarts(tree.prior = cgm(power = 3))` call gets one warning a session and fits the same
model. `split.probs = dart()` turns on DART at `bart()`, `dbarts()`, `dbartsSpec()` and `xbart()`.

Two defects on main go: `xbart(tree.prior = cgm(power = 3))` fitting the default depth, and a forest's
`power` silently beating the fitting function's. The second becomes a refusal by name, as `n.trees` given at both
levels already is.

The live door. setmodel-parameters' slice 4 had planned `$setTreePrior(tree.prior = cgm(...))`; under your
ruling it takes `power`, `base` and `split.probs` plainly, with `forests` and `variance` as planned. One choice
there is not settled by a ruling: what `$setTreePrior(power = 3)` does when the only forest with no basis stated
its own power at creation (open call 1). I recommend refusing it with the exact call to write.

bartCause and stan4bart drop their translations into `cgm()` and pass the plain arguments through; treatSens
does not change. Their branches merge right after this one.

About 1450 lines with tests and docs, and about 120 in the two consumers; about 19 agent hours. Plans here run
1.5 to 2 times low, so budget 28 to 38.

## Goal

Every door takes a tree's depth as `power` and `base` and its split probabilities as `split.probs`, a vector or
`dart()`; `tree.prior` is refused on `bart()` and `xbart()` with the plain spelling to write, and accepted on
`dbarts()` and `dbartsSpec()` with a once-per-session warning until 1.1-0; a setting stated on the fitting
function and on every forest with no basis is refused by name; and every call main accepts and this item keeps
draws what it drew.

## What exists today

Run on main at 83cefe7d (shipped build, arm64 macOS, private library scratch/libs/tpf-plan; probes
scratch/tpf-plan/q0.R to q6.R, untracked, and the design agent's scratch/tree-prior-design/p1.R to p10.R rerun
there as q2, q5 and q8). Data: 120 rows, `x1` to `x3` uniform, a 0/1 `z`, a gaussian `y`.

- 0.9-34 (scratch/libs/cran): `bart2` takes `power`, `base` and `split.probs` plainly; `dbarts` takes
  `tree.prior` (default `cgm`, resolved by bare name, `cgm(power = 2, base = 0.95, split.probs = 1/num.vars)`) and
  no `power`; `xbart` takes `power` and `base` grids and no `tree.prior` or `split.probs`; `rbart_vi` takes
  `power` and `base`; `bart` is BayesTree's (`power`, `base`, `splitprobs`). No `dart` exists. Nothing exports
  `cgm`.
- `bart()` takes `tree.prior` as a formal; `power`, `base` and `split.probs` arrive through `...` and warn once
  ("'power' has left 'bart': ... write tree.prior = cgm(power = ) or tree.prior = dart(power = )"), the values
  used ([`resolveConsolidatedArgs`](../../R/tombstones.R), [`consolidatedArgReasons`](../../R/tombstones.R));
  [`buildSamplerPriors`](../../R/bart.R) rebuilds them into a `cgm()` call for `dbarts()` and refuses
  `tree.prior` beside any of them ([`refuseColliding`](../../R/model.R)).
- `dbarts()` and `dbartsSpec()` take `tree.prior = cgm`; `dbarts(power = 3)` is "unused argument 'power'
  passed to 'dbarts'". [`parsePriors`](../../R/model.R) resolves the object; [`resolveSplitProbabilities`](../../R/model.R)
  normalizes a vector against the columns, with `num.vars` and `numvars` in its vocabulary.
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
  and base 0.5 from `cgm(power = 3, base = 0.5)` (q3). The writer is
  [`resolveSamplerSpec`](../../R/spec.R), which writes the first forest's knobs over the model's prior. `n.trees`
  given to `bart()` and to that forest is refused: "'n.trees' is given to the fitting function and to the
  forest with no basis ('forest(x1 + x2, n.trees = 10)'), which are the same count; give one"
  ([`refuseTreeCountGivenTwice`](../../R/formulaTerms.R)); at `dbarts()` a control's count beside the forest's
  leaves the forest's to govern.
- No forest with no basis. `dbarts(y ~ forest(x1, basis = factor(z)) + forest(x2, basis = ~ x3), tree.prior =
  cgm(power = 3))`: "'tree.prior' given to the fitting function is the tree prior of the forest with no basis,
  and every forest of this model has a basis; state it on a forest, as forest(x1, basis = a, base = 0.25, power =
  3)" ([`refuseStatedWithNoPlainForest`](../../R/model.R), dec-B274). `split.probs` there gets the same remedy,
  which no forest can take.
- Several forests refuse DART and split probabilities ("a treatment forest does not support a DART tree prior";
  "... 'split.probs'"); `family = "multinomial"` refuses both at `bart()` and in the specification. A variance
  forest cannot join several forests ("a treatment forest does not support 'variance'").
- The variance forest. An unstated one takes the model's base and power at creation
  ([`resolveSamplerSpec`](../../R/spec.R)): 0.5 and 3 under `cgm(power = 3, base = 0.5)`, 0.95 and 3 under
  `dart(power = 3)`; split probabilities never reach it (q3).
- `xbart()` takes `power` and `base` grids, `split.probs`, and `tree.prior` (`cgm()` or `dart()`); `xbart(split.probs
  = dart())` is "could not find function 'dart'".
- `pdbart(y ~ ., power = 3)` warns "'power' is BayesTree's spelling; 'pdbart' now fits through 'bart', which takes
  'tree.prior = cgm(power = )'" ([`pdbartBayesTreeNames`](../../R/tombstones.R),
  [`translatePdbartCall`](../../R/tombstones.R)), and rewrites into `tree.prior`.
- `bartBT()` takes `power`, `base`, `splitprobs` and builds a `cgm()` call; `rbart_vi()` (deprecated, removed in
  1.1-0) takes `power`, `base`, builds `cgm(power, base)`, and refuses `split.probs` and `tree.prior`.
- The live sampler. `$setModel` with `m@tree.prior@power <- 1` is honored on a one-forest sampler, refused on a
  DART one and on one of several forests; `$setLeafPrior(forests = list(forest(power = 1)))` is "'power' is fixed
  at creation" (q8). No `$setTreePrior` exists.
- Consumers. bartCause's `bcf()` builds `quote(cgm(power, base))` for `dbarts(tree.prior = )` and states
  `treatment.base`, `treatment.power` on the treatment forest; its `translateRetiredBartArguments` turns `bartc`'s
  `power`, `base`, `split.probs` into `tree.prior = cgm(...)` for `bart()`, and its `retiredSpellings` lists
  them as `tree.prior`'s; its tests read `model@tree.prior@power` and pass `dbartsPriors$cgm(2.0, 0.95)` to
  `dbarts()`. stan4bart rebuilds `bart_args`' `power`, `base`, `split.probs` into `cgm()` for `dbartsSpec()` and
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
  the three goes. The `family = "multinomial"` refusal reads `split.probs` (a vector or `dart()`) instead of the
  object's class. [`bart2Hurdle`](../../R/bart.R) hands the three to both components where it handed
  `tree.prior`. [`bartBT`](../../R/bart.R) forwards `power`, `base` and `split.probs` plain. The stated record
  (`plainStatedAttr`) names each of the three separately.
- R/tombstones.R. `power`, `base`, `split.probs` leave `consolidatedArgReasons`, `consolidatedPriorScalars`,
  `consolidatedArgsFor$bart` and the registry rows owned by `bart`. A `tree.prior` row for `dbarts` and
  `dbartsSpec` joins, carried on `...` as retired names are, with its warning ([Design](#design) D3). The
  refusal of `tree.prior` at `bart()` and `xbart()` names the plain spelling (D4). `power` and `base` leave
  [`pdbartBayesTreeNames`](../../R/tombstones.R) and its spelling table: they are `bart()`'s own names, passed
  through with no warning; `splitprobs` maps to `split.probs` and still warns. The rewrite into `cgm()` in
  [`translatePdbartCall`](../../R/tombstones.R) goes.
- R/model.R. [`parsePriors`](../../R/model.R) builds the record from `power`, `base` and `split.probs` (D1, D2).
  [`dart`](../../R/model.R) loses `power` and `base`. [`cgm`](../../R/model.R) stays in
  [`dbartsPriors`](../../R/model.R) as `tree.prior`'s retired value until 1.1-0. `plainForestArgNames` and
  [`plainForestStated`](../../R/model.R) record `power`, `base`, `split.probs` by name. The collision check (D5)
  is one helper beside [`refuseStatedWithNoPlainForest`](../../R/model.R), whose remedy names the argument given.
  [`forest`](../../R/model.R)'s and [`varianceForest`](../../R/model.R)'s comments, and the
  `dbartsVarianceForest` format's "<mean forest's>", say "the fitting function's".
- R/formulaTerms.R. [`refuseTreeCountGivenTwice`](../../R/formulaTerms.R) becomes the formula half of D5 for any
  of `n.trees`, `power`, `base`.
- R/spec.R. [`resolveSamplerSpec`](../../R/spec.R): the fitting function's values come from the plain
  arguments; the first forest's knobs are written over the record only where the fitting function's did not
  reach it (D5 refuses the rest); the variance forest reads the fitting function's base and power, as now
  (which equals the model's, a variance forest standing only beside one forest). The several-forest and
  multinomial refusals name `split.probs`. [`dbartsSpec`](../../R/spec.R) gains the three plain arguments and
  takes `tree.prior` through `...`.
- R/dbarts.R. [`dbarts`](../../R/dbarts.R) the same. `$setModel`'s DART refusal text (replaced on
  setmodel-parameters' branch) and every other read of `model@tree.prior` stay.
- R/xbart.R. [`xbart`](../../R/xbart.R) loses `tree.prior` (D4); `split.probs` takes `dart()`, resolved in the
  split-probability vocabulary; [`cellModel`](../../R/xbart.R) as now. The comment defending the grid's
  override of an object's power goes.
- R/rbart.R. [`rbart_vi`](../../R/rbart.R) forwards `power` and `base` plain; nothing else (deprecated).
- R/A_class.R, R/bartcore.R ([`setModelRefused`](../../R/bartcore.R)), R/forestBasis.R (`cgm` and `dart` among the
  names reserved in a basis term): unchanged.
- src/R_interface_bartcore.cpp ([`bartcore_setModel`](../../src/R_interface_bartcore.cpp), the creation reads,
  the verbose line) and inst/include/dbarts/dbarts.h: unchanged. tests/cpp: one comment names `cgm`; unchanged.

Help, NEWS, docs (slice 2):

- man/bart.Rd: usage, `power`, `base`, `split.probs` items (a vector, `NULL` uniform, `num.vars`, or `dart()`),
  `tree.prior` gone from the `tree.prior, leaf.prior` item, examples. man/dbarts.Rd and man/dbartsSpec.Rd: the same
  three; `tree.prior` described as retired. man/xbart.Rd: `tree.prior` item gone, `split.probs` takes `dart()`, the
  `...` item's "dart, now tree.prior = dart()" respelled. man/dbartsPriors.Rd: `dart(a, b, rho, alpha,
  update.alpha, update.delay)`, given as `split.probs`; `cgm` marked as `tree.prior`'s retired spelling; examples.
  man/forest.Rd, man/varianceForest.Rd: "the fitting function's `power` and `base`", the governing paragraph
  rewritten for D5. man/monotone.Rd: the mBART prior as `power = 0.8, base = 0.25`. man/dbarts-deprecated.Rd: the
  `power`, `base`, `split.probs` items go; a `tree.prior` item on `dbarts` and `dbartsSpec`. man/pdbart.Rd,
  man/dbarts-package.Rd, man/dbartsSampler-class.Rd (the `setModel` DART sentence): one line each. man/bartBT.Rd:
  unchanged.
- inst/NEWS.Rd: [NEWS](#news).
- vignettes/gibbs_sampler_mixture_model.Rmd (`tree.prior = dart()`) and vignettes/dbarts-as-a-component.Rmd
  (`tree.prior = dbartsPriors$cgm(power = 0.5)`), README.md: respell.
- docs/design: [3a. Prior specification](../design/public-surface.md#3a-prior-specification) and
  [4. DART exposure](../design/public-surface.md#4-dart-exposure) gain a dated revision paragraph;
  [The rule](../design/forest-defaults-by-kind.md#the-rule) and
  [What counts as stated](../design/forest-defaults-by-kind.md#what-counts-as-stated) are rewritten for the
  plain names and D5; monotone.md, prior-defaults.md, consumer-spec-surface.md, pdbart-on-bart.md and
  bart-landscape.md respell sentences stating present behavior. The rest (nog-gibbs.md, core-generalization.md,
  level-fiber.md, robust-errors.md, written-surface.md, interaction-constraints.md, multinomial.md,
  feature-matrix.md) are records or name the prior in passing: a sentence stating present behavior is respelled,
  a dated record is not.

Tests and benchmarks:

- inst/tinytest: about 200 call sites in 39 files state `tree.prior = cgm(...)`, `dart(...)` or the retired
  `power`/`base`/`split.probs` on `bart()`; each is respelled, and the pins of the retired spellings' warnings
  (test-tombstones.R, test-front-door-control.R, test-argument-surface.R, test-pdbart.R, test-model-priors.R,
  test-forest-defaults.R, test-xbart-model.R) are inverted or replaced. The four snapshot files state none.
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
is added here (a per-forest tree-prior reader belongs with slice 4's door).

D2. `dart()` takes `a`, `b`, `rho`, `alpha`, `update.alpha`, `update.delay` and returns the record class with
`power` and `base` at 2 and 0.95, which the door overwrites with its own. `dart(power = 3)` is R's
unused-argument error, as `forest(leaf.prior = )` is under dec-B508; `dart` is new in 1.0-0. A bare `dart` is
`dart()`, as `tree.prior = cgm` meant `cgm()`. `split.probs` resolves in the vocabulary `dart`, `num.vars`,
`numvars`, layered over the caller's frame, at every door, and a forwarded `..N` is recovered as now.

D3. `dbarts()` and `dbartsSpec()` take `power = 2`, `base = 0.95`, `split.probs = NULL` where `tree.prior` stood
(after `n.samples` and after `control`), and `tree.prior` on `...`. A `tree.prior` given there warns once a
session, class `dbartsDeprecatedWarning`: "'tree.prior' has left 'dbarts': write the tree prior as power, base and
split.probs (split.probs = dart() for DART), as bart() takes them; 'tree.prior' is removed in dbarts 1.1-0. The
value was used." It is then read as the plain arguments it states: a written `cgm(...)` or `dart(...)` call
states the arguments it names (positional ones matched to `cgm`'s formals: power, base, split.probs), so
`cgm(power = 3)` states `power` alone; a value, an object built elsewhere, states all of its slots; `dart(...)`
is `split.probs = dart(...)`; a bare `cgm` states nothing. `tree.prior` beside any of `power`, `base` or
`split.probs` is refused: "'tree.prior' and 'power' both set the tree prior of 'dbarts'; give 'power' alone".

D4. `bart()` and `xbart()` take no `tree.prior`, which never shipped on either. Given one, through `...`, each
refuses it by name with the spelling D3 would translate it to: "'tree.prior' is not an argument of 'bart': write
the tree prior as plain arguments, power = 3, base = 0.5" (for `dart(...)`, "split.probs = dart(...)"). On
`xbart()` the message adds that `power` and `base` are grid axes. So `xbart(tree.prior = cgm(power = 3))` can no
longer fit the default depth; `xbart(power = 3)` fits power 3, as on main.

D5. Two levels, one rule. The fitting function's `power` and `base` reach every mean forest with no basis that
does not state its own, and a variance forest that states none (dec-B512); a forest stating its own keeps it. A
name given to the fitting function that reaches no mean forest is refused by name, before anything is built:

- where every mean forest has a basis (dec-B274, now per name): "'power' given to the fitting function is the
  tree prior of the forest with no basis, and every forest of this model has a basis; state it on a forest, as
  forest(x1, basis = a, power = 3)", the remedy naming the argument given;
- where every forest with no basis states it (dec-B512, as for `n.trees`): "'power' is given to the fitting
  function and to the forest with no basis ('forest(x1 + x2, power = 1)'), which are the same setting; give
  one" (the code where a formula wrote it, "(forest 1)" in a list or a data object's `bases`).

Stated means named, the default's own value included, as for `n.trees`: `bart(power = 2)` beside
`forest(x, power = 1)` is refused, and `bart()` forwards which of the three its caller named. Names are checked
one by one: `power` on the fitting function and `base` on the forest are both taken. A variance forest never
rescues a value: `bart(y ~ forest(x, power = 1), power = 3, variance = ~ x3)` is refused, the lone mean forest
stating `power` twice. The helper takes the positions of the forests with no basis, so forests-shared-total, which
makes them several, only widens `plainForest`'s answer; its own `n.trees` rule (refused only where every sharer
states one) is this rule. A lone forest written as a term, `y ~ forest(x, power = 1)`, is the model's forest, and
`power = 1` there is the model's tree prior, as `sd` there is its leaf prior (dec-B266).

D6. `split.probs` is the mean forest's alone, as on main: the variance forest takes `base` and `power` and never
split probabilities, and `forest()` and `varianceForest()` take no `split.probs`. In a model of several forests,
and under `family = "multinomial"`, `split.probs` (a vector or `dart()`) is refused first, by its own message,
not by D5's "state it on a forest": "a model of several forests does not support 'split.probs' (a vector or
dart()); drop it or fit a single-forest model".

D7. The live door: setmodel-parameters' slice 4 ([The live door](#the-live-door)).

D8. `cgm` stays in `dbartsPriors`, documented as the value of the retired `tree.prior`, and leaves with it at
1.1-0, its class staying as the record's. A TODO entry with the 1.1-0 removals records it.

### The live door

Slice 4 of setmodel-parameters (on its branch, planned at a628e48e) is to give a live sampler the tree priors of
every forest, through `$setTreePrior(tree.prior = , forests = , variance = )`, its open call 1 (a). Under dec-B512
the door becomes:

    $setTreePrior(power = NULL, base = NULL, split.probs = NULL, forests = NULL, variance = NULL)

- `power`, `base`: the fitting function's. Each reaches every forest with no basis that does not state its own,
  and a variance forest that stated none, as at creation; NULL leaves it as it is. On a sampler of one forest the
  forest's statement is the model's, so it always takes them. Where the forests with no basis all state the name
  it is refused, naming the call to write (open call 1): "'power' given to $setTreePrior reaches no forest: the
  forest with no basis (forest 1) states its own; change it with forests = list(forest(power = 3), forest())".
- `split.probs`: a vector changes the values (a parameter, dec-B254); `dart(...)` on a DART sampler changes
  DART's settings, a setting it leaves out taking creation's default for it, as creation with `dart(...)` would
  (dec-B514's principle; slice 3's install); `1` restores equal probabilities, as creation with none has; a
  vector on a DART sampler or `dart()` on an ordinary one is structure and is refused: "split.probs = dart()
  makes a DART sampler of an ordinary one, which changes the model's structure; make a new sampler instead".
- `forests`: a list of `forest(power = , base = )`, by position or label, a forest set there stating its own from
  then on (slice 4 item 4); any other knob refused as fixed at creation, `sd` naming `$setLeafPrior`. Stated at
  the top and on that forest in one call, D5 refuses it.
- `variance`: `varianceForest(power = , base = )`.
- No `tree.prior` argument: it never shipped on a live sampler. `$setModel` keeps taking model objects
  (dec-B512), and slice 4's install from the record (its item 3) is unchanged.

Under dec-B512 an unstated variance forest takes the fitting function's tree prior, and under dec-B514 a setter
gives what creation with its statement gives; together they settle slice 4's open call 2 as its (a): an unstated
variance forest follows `$setTreePrior`'s and `$setModel`'s `power` and `base`. That moves one call main honors (a
heteroscedastic `$setModel` with a new tree prior), gated in slice 4 by the exact arm its plan names. Calls made 12.

This item builds none of the door: slice 4 builds it after forests-shared-total. The amendment to that plan
(its door, item 4, its open calls 1 and 2, its tinytest pins of `tree.prior`) is a records step for that branch's
writer, before slice 4 is built.

## Collisions

Every pair of levels at which one setting can be stated twice, and what happens.

| levels | example | result |
|---|---|---|
| fitting function, forest with no basis | `bart(y ~ forest(x, power = 1) + forest(x, basis = z), power = 3)` | refused, D5 |
| the same, other names | `bart(... forest(x, base = 0.5) ..., power = 3)` | power 3, base 0.5 |
| fitting function, forest with a basis | `bart(... forest(x, basis = z, power = 4), power = 3)` | taken: 3 on the forest with no basis, 4 on the other |
| fitting function, no forest with no basis | `dbarts(y ~ forest(x1, basis = a) + forest(x2, basis = b), power = 3)` | refused, D5 (dec-B274) |
| fitting function, variance forest | `bart(y ~ x, power = 3, variance = varianceForest(x3, power = 1))` | taken: 3 and 1 |
| fitting function, lone forest term | `bart(y ~ forest(x, power = 1), power = 3)` | refused, D5 |
| `tree.prior` (retired) and the plain names | `dbarts(..., tree.prior = cgm(), power = 3)` | refused, D3 |
| `tree.prior` (retired), forest | `dbarts(y ~ forest(x, power = 1) + ..., tree.prior = cgm(power = 3))` | refused, D5, after D3's warning |
| `tree.prior` object, forest, other name | `tp <- cgm(power = 3); dbarts(... forest(x, base = 0.5) ..., tree.prior = tp)` | refused: an object states base too |
| `tree.prior` at `bart()`, `xbart()` | `bart(y ~ x, tree.prior = cgm(power = 3))` | refused, D4, naming `power = 3` |
| pdbart's BayesTree name and bart's | `pdbart(y ~ ., splitprobs = p, split.probs = p)` | refused, as now ("both set one setting") |
| bcf()'s and bart()'s, bartCause | `bcf(..., power =, treatment.power =)` | two forests; no collision |
| stan4bart `bart_args` | `bart_args = list(power = 2.5, tree.prior = cgm())` | refused by dbartsSpec, D3 |
| live door, top and forest | `$setTreePrior(power = 3, forests = list(forest(power = 2), forest()))` | refused, D5 |
| live door, top, forest stated at creation | `$setTreePrior(power = 3)`, the forest with no basis having stated power | refused naming `forests =` (open call 1) |
| live door and `$setModel` | separate calls | each installs; the later stands |

A value no forest can take, each refused by name: the fitting function's `power` or `base` where no mean forest
has no basis or every one states it (D5); `split.probs` in a model of several forests or under multinomial (D6);
`split.probs = dart()` swapped for a vector on a live sampler or the reverse (structure); `power` not a single
positive number and `base` not a single number in (0, 1), at every door with the messages `forest()` gives
("'power' must be a single positive number", "'base' must be a single number in (0, 1)"), `xbart()` keeping its
grid's "must contain only" wording. `forest(split.probs = )`, `varianceForest(split.probs = )`,
`forest(tree.prior = )` and `dart(power = )` are R's unused-argument error.

## NEWS

Scope is what differs from 0.9-34; nothing that never reached a release.

- The tombstone item loses "On bart, power, base and split.probs become tree.prior = cgm(power, base,
  split.probs)": `bart()` takes them as 0.9-34's `bart2` did. It gains: on `dbarts`, `tree.prior` becomes `power`, `base` and `split.probs`, as on `bart`; `tree.prior = cgm(power, base,
  split.probs)` is accepted with a once-per-session warning, the value used, and is refused from 1.1-0. A 0.9-34
  `dbarts(tree.prior = cgm(power = 3))` caller sees that warning once and fits the same model.
- The "Tree priors" item becomes: `split.probs = dart()` runs DART (Linero 2018) on `bart`, `dbarts`,
  `dbartsSpec` and `xbart`; `split.probs` on `xbart` is new, a vector or `dart()`. "`bart` gains `tree.prior`"
  goes; `leaf.prior` stays.
- `pdbart` and `pd2bart`: `power` and `base` are no longer among the BayesTree spellings it warns of.
- Nothing about the xbart defect or the two-level refusal: neither reached a release (0.9-34 had no forests and
  no `tree.prior` on `xbart`).

## Order of landing

1. This item, slices 1 and 2 on branch tree-prior-flat, one review, the gates; then bartCause's and stan4bart's
   branches (slice 3), built and green against this branch's tip installed from the branch; then dbarts's merge,
   followed in the same sitting by the two consumer merges, and their CI re-run on main. Between the dbarts merge
   and bartCause's, `bartc(power = )` is refused (D4); the consumers' branches are ready first so the window is
   minutes. A compatibility shim in the consumers (reading `formals(dbarts::bart)`) would remove the window at the
   cost of dead code; not taken (Calls made 14).
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
4. setmodel-parameters: slices 1 to 3 do not meet this item but textually (slice 1's `$setModel` refusals in
   R/dbarts.R and R/bartcore.R; slice 3's DART install reads the unchanged record). This item lands before slice
   4, whose door is amended as [The live door](#the-live-door) says; the branch, which merges once after slice 4,
   rebases onto it.
5. control-migration (ruled, unplanned): it puts `n.trees`, `n.cuts`, `useQuantiles` on `dbarts()` beside this
   item's three; D5's helper covers `dbarts(n.trees = )` beside a stated forest count. Whichever is second
   rebases; this item is expected first.

## Slices

### Slice 1: the doors and the record (R)

Steps: the R inventory above, D1 to D6 and D8; the tinytest re-pins and the new test file.

1. [`parsePriors`](../../R/model.R) takes `power`, `base`, `split.probs` (and the retired `tree.prior`, translated
   by D3 first) and builds the record; `dart()` per D2.
2. `dbarts()`, `dbartsSpec()`: formals, the `tree.prior` tombstone and translation, the both-given refusal.
3. `bart()`, `bartBT()`, `bart2Hurdle()`, the multinomial path, `rbart_vi()`: forward the plain three; the
   `bart()` tombstones for them go; D4's refusal.
4. `xbart()`: D4, `split.probs = dart()`.
5. D5: the stated record per name, the helper for formula, list and data-object doors, the per-name refusal
   where no forest has no basis; D6's refusal first.
6. pdbart's table.
7. Tests ([Tests and gates](#tests-and-gates)).

RNG: neutral for every call main accepts; refusals added. Cost: ~990 lines; 7 agent hours.

### Slice 2: help, NEWS, docs, benchmark scripts

Steps: the help, NEWS and docs inventory; respell the 28 benchmark scripts and the scenario in equivalence.R; the
records' TODO entry for 1.1-0 (D8). Cost: ~420 lines; 3 agent hours.

### Slice 3: consumers, each on its own branch

- bartCause: `bcf()` sets `samplerCall$power` and `$base` instead of `tree.prior`; `translateRetiredBartArguments`
  stops building `cgm()` (the three pass through to `bart()`), and `retiredSpellings` loses its `tree.prior`
  entry; test-14-bcf.R's `dbarts(tree.prior = dbartsPriors$cgm(2.0, 0.95))` respelled `power = 2, base = 0.95`.
  Its tests reading `model@tree.prior@power` stand (D1). ~70 lines, 1.5 agent hours.
- stan4bart: `bart_args`' `power`, `base`, `split.probs` pass to `dbartsSpec()` as its own formals (the
  `known_bart_args` list already takes every `dbartsSpec` formal); the `cgm()` rebuild and its refusal go;
  `bart_args$tree.prior` reaches `dbartsSpec()`, which warns (D3); test-09-bartArgs.R's "both 'tree.prior' and
  'power'" pin takes dbarts's message, and its `tree.prior = cgm(power = 2.5, base = 0.9)` arm counts the one
  warning. The verbose line is unchanged. ~45 lines, 1.5 agent hours.
- treatSens: nothing (searched: one `dbartsSpec(forests = list(forest(n.trees = )))`, `pdbart()` with neither).

Both: the review over slices 1 and 2, 2 hours; the gates, 2 hours; the fix round, 1.5; records, 0.5. About 19
agent hours; plans here run 1.5 to 2 times low, so budget 28 to 38.

## Tests and gates

tinytest, a new test-tree-prior-flat.R, small n, two chains where draws are compared:

- Every door takes the plain three and writes main's record: `bart()`, `dbarts()`, `dbartsSpec()`, `xbart()`
  (per cell, through the fit's model), `bartBT()`, `rbart_vi()`, `pdbart()`, the hurdle family's two components,
  monotone; `model@tree.prior` `identical()` to the record a constructed `cgm()` or `dart()` with the same values
  gives after resolution.
- `split.probs`: a named vector, `c(2, 1, 1) / num.vars`, `dart()`, bare `dart`, `dbartsPriors$dart(a = 1)`, through a
  wrapper's `...`; DART's `update.delay` default half the burn-in; `dart(power = 3)` an unused-argument error;
  `dart(0.3)` is `a = 0.3`.
- D3: `dbarts(tree.prior = cgm(power = 3))` and `dbartsSpec()` warn exactly once a session (counted with
  `withCallingHandlers` over two calls, after resetting the once-per-session key), fit power 3, base 0.95; `dart()`
  there is DART; a bare `cgm`; an object; `tree.prior` beside `power` refused.
- D4: `bart(tree.prior = cgm(power = 3, base = 0.5))` and `xbart(tree.prior = dart(a = 1))` refused, each message
  naming the translation.
- D5, by door (formula, `forests` list, a data object's `bases`) and by name: the refusal where the forest with no
  basis states it, the default's own value included (`bart(power = 2)`); other names taken; a forest with a basis
  taken; no forest with no basis, each name refused with its own remedy; the variance-forest case refused; a lone
  forest term alone taken as the model's.
- D6: `split.probs` vector and `dart()` beside two forests and under multinomial, with the new message.
- The variance forest: unstated takes the fitting function's `power` and `base` (its record 3 and 0.5), stated
  keeps its own, `split.probs` never reaches it.
- pdbart: `power`, `base` pass with no warning; `splitprobs` still warns and maps.
- Continuation: a DART sampler made with `split.probs = dart()`, run, `copy()`, `storeState` with `saveRDS` and
  `readRDS`: next draws `identical()`; `$setModel` with `m@tree.prior@power <- 1` on a one-forest sampler taken,
  as on main.
- Re-pins: the 39 files respelled; every pin of a retired-spelling warning on `bart()` inverted to no warning.

Gates, on the branch's own library:

- `R CMD INSTALL -l <lib>`; `buildInfo()$mode` checked.
- A neutrality A/B (scratch, ~250 lines): an 83cefe7d library against the branch's, two chains, every recorded
  channel compared with `identical()`: each door in each spelling main accepts against its plain respelling on the
  branch (`bart(tree.prior = cgm(...))` against `bart(power = ...)`, `dbarts(tree.prior = dart(...))` against
  `split.probs = dart(...)`, xbart's grid with and without DART, `bartBT`, `rbart_vi`, `pdbart`, hurdle, monotone,
  heteroscedastic with and without a stated variance tree prior, bcf-shaped forests with the forest with no basis
  first and second, the default call), and the default call with no tree prior stated, which must not move. A
  difference is a stop.
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
- D3: an object read as stating only its non-default slots (the object pin); a written `cgm(power = 3)` read as
  stating `base` too (the written cross-name pin); the warning every call, or none (the count pin).
- D6's refusal after D5's, giving "state it on a forest" for `split.probs` (the message pin).
- `split.probs = dart()` reaching the variance forest's record, or the variance forest losing the fitting
  function's `base` (the variance record pins).
- `dart()` keeping `power` as its first formal (the `dart(0.3)` pin; equivalence's `dart` scenario).
- `split.probs` evaluated in the caller's frame (the `num.vars` pin; the wrapper pin).
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
4. The retired `tree.prior` is read by the names a written call states, and as all its slots where it is a value
   (D3). Alternative: all slots always, which refuses `tree.prior = cgm(power = 3)` beside a forest's `base`, a
   call main fits as it was meant.
5. `tree.prior` on `bart()` and `xbart()` is refused with the translation, not accepted with a warning (D4):
   dec-B512 retires it with a warning only where it shipped. The translated message is the cheap part of the
   hand-holding. Alternative: a plain unused-argument error.
6. `split.probs` stays the mean forest's: no `split.probs` on `forest()` or `varianceForest()`, none reaching the
   variance forest (D6), as main. Alternative: `forest(split.probs = )` on a lone forest, ~40 lines, refused in
   every model of several forests, where the engine has none.
7. One rule for two levels (D5): a fitting function's value that reaches no mean forest is refused, which is
   dec-B274's refusal and dec-B512's collision at once, and forests-shared-total's `n.trees` rule for several
   forests. The variance forest does not rescue a value. Alternative: the variance forest counting as reached,
   which takes `bart(y ~ forest(x, power = 1), power = 3, variance = ~ x3)` with power 1 on the mean forest.
8. `split.probs` beside several forests is refused by its own message before D5 (D6); main's "state it on a
   forest" names a remedy no forest has.
9. A positional `dart(2, 0.95)` now means `a = 2, b = 0.95`, unguarded: `dart` never shipped, and the one such
   call in the tree (equivalence.R) is respelled.
10. `cgm` stays until 1.1-0 and leaves with `tree.prior` (D8); its class stays as the record.
11. pdbart's `power` and `base` leave the BayesTree table, so 0.9-34's `pdbart(..., power = 3)` warns no more.
12. Slice 4's open call 2 is settled by dec-B512 with dec-B514, as its (a) ([The live door](#the-live-door)).
    Alternative: its (b), an unstated variance forest keeping creation's values after a live change, against
    "gives what creation gives".
13. The live door takes no `tree.prior`, and `split.probs = 1` restores equal probabilities there.
14. Consumers merge right after dbarts, from branches made green against it first. Alternative: a shim in each
    consumer that reads dbarts's formals, which removes a minutes-long window at the cost of dead code.
15. Tier "R surface, no draw change", one sonnet review, with the bitwise A/B and equivalence trio as the
    oracle that no draw moves. Alternative: an opus implementer, as forest-defaults-by-kind took for the block
    that decides what model a call fits; the A/B is what catches that class here.

## Open calls

1. For setmodel-parameters' slice 4, not for this item: `$setTreePrior(power = 3)` on a sampler of several forests
   whose forests with no basis all stated their own `power` at creation.
   - (a) Recommended. Refused, naming the call to write: "'power' given to $setTreePrior reaches no forest: the
     forest with no basis (forest 1) states its own; change it with forests = list(forest(power = 3), forest())".
     The same rule as creation, where the two statements in one call are refused; a sampler of one forest always
     takes it. Cost: one refused call with its remedy.
   - (b) Taken onto every forest with no basis, clearing their own statements, the latest word winning. Cost: a
     forest's own statement means "kept" at creation and under a group change elsewhere, and "overwritten" here;
     with several forests with no basis, one call can silently undo a statement the caller made per forest.
   - (c) Taken and changes nothing, the forests keeping their own, as dec-B514 keeps a stated sd under a bare
     `normal()`. Cost: a call that does nothing with no message.

## Not verified

- Nothing is built; the sizes are from the inventory, not a prototype.
- That every one of the ~200 tinytest and ~60 benchmark sites respells without a change of value is assumed from
  the grep; the A/B, the equivalence trio and the exact gates are the check.
- 0.9-34's `dbarts(tree.prior = cgm(3, 0.5, c(0.5, 0.25, 0.25)))` positional form was run on 0.9-34 only; its
  translation under D3 is planned, not run.
- setmodel-parameters' slice 4 and forests-shared-total were read on their branches at a628e48e and 0f16f8d5;
  their calls may move. The control migration has no plan.
- The consumers' changes are from reading their code at bartCause f7edf5b and stan4bart 566a1d9.
