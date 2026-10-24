# leaf-one-row-fallback: fewer than two rows in read a leaf covariate from every row, under every family

Status: BUILT 2026-10-24 on wt/leaf-one-row off bartcore c41f25c0, not landed; implements dec-B485.

agent: one opus implementer (engine, R text, tests, help, notes); one opus review that runs mutants.
rng: draw-changing only where a case refused today now runs. NEUTRAL, bit for bit, for every other call, with
one corner that moves (below, [Built](#built)). No recorded scenario has one row in with a leaf covariate:
expected 55, 19 and 11 identical against equivalence-32418b6e, bcf-equivalence-20eed7ed and
multinomial-equivalence-80b1c8d4, no re-record.
window: before the merge (TODO `leaf-one-row-fallback`).
budget: ~300 lines of package code plus tests; built at 40 lines added and 33 removed, most of them comments.

## Goal

Where fewer than two rows are in the likelihood (positive weight, not masked out) at a derivation of a linear or
gp leaf's standardization, the covariates' center and scale and a gp leaf's default lengthscale are taken over
every row, with one warning that counts the rows in, under every family, at creation, at a `setData` that derives
and in an `xbart` fold. Gaussian, Student-t and nbinom keep dec-B479's trigger and text; two or more rows in at one
value of a covariate stay refused (dec-B465); a re-derivation by the row setters keeps (dec-B364).

## Context

The ruling is dec-B485, revising that call of dec-A226; it rests on dec-B302 (no row in reads every row),
dec-B479 (the response's fallback carries the leaf) and dec-B465. Push 1 of
[response-scale-rows-d.md](response-scale-rows-d.md) built the rows-in rule and the response fallback
([Push 1 as built](response-scale-rows-d.md#push-1-as-built),
[Push 1 fix round](response-scale-rows-d.md#push-1-fix-round)); the standing statement is rule 9 of
[response-scale-rows.md](../design/response-scale-rows.md). The engine decides the rows in
[`Sampler::leafRowsFrom`](../../src/bartcore/sampler.hpp), every door reading them through it: creation and an
`xbart` fold through the constructor tail ([`Sampler::initializeChains`](../../src/bartcore/sampler.hpp)), a
replacement through [`Sampler::setData`](../../src/bartcore/sampler.hpp), a re-derivation through
[`Sampler::rederiveStandardization`](../../src/bartcore/sampler.hpp). R words the warning
([`warnRowsFallBack`](../../R/bartcore.R)).

## Built

[`Sampler::leafRowsFrom`](../../src/bartcore/sampler.hpp) marks one row in out of more than one as it marked no row
in, so every family and every door that derives reads every row there and flags the covariates
(`DerivedRecords::leafFromAllRows`, which the creation and `setData` warnings read). The flag of
[`LeafRows`](../../src/bartcore/model.hpp) is renamed `noneIn` to `tooFewIn`, which it already meant under
dec-B479. [`Sampler::withResponseFallback`](../../src/bartcore/sampler.hpp) is unchanged: one row in now reaches it
already marked, the same `LeafRows` it returned, so gaussian, Student-t and nbinom take the same bits and the same
text. A re-derivation at one row in keeps every marked covariate, as it did: the covariates held one value there
before, and are kept as too few now. [`warnRowsFallBack`](../../R/bartcore.R) gains the one-row lead,
"one observation is in the likelihood (1 of 200 has a positive weight), so the center and scale of "x1" are
taken from all 200"; the bridge's refusal needs no change, the count it is handed being 0 where every row was read
("over all 200 observations, read because too few are in the likelihood").

The other families. Today probit, ordinal and hazard.probit refuse one row in through `dbarts()`, `bart()` and a
deriving `setData`; logistic and aft reach the same refusal only by hand, a mask given to `new()` (logistic's
weights are counts and aft takes none), aft's data object without the record its fit wrote. All five now read
every row with the one-row text. No `xbart` fold reaches one row in outside gaussian: `xbart` refuses probit's 0/1
weights, logistic's weights are positive counts and it fits no other family, and a fold's sampler takes no mask;
a gaussian fold with one row in falls back with its response, as before. The engine path a fold takes, the
constructor tail over a store view, is pinned in tests/cpp.

Calls made in building, each with its alternative:

- The trigger is in `leafRowsFrom`, by count, ahead of the response's test. Alternative: a count test in
  `withResponseFallback`, which `setData` calls only where it derives the range, so a gaussian `setData` keeping the
  range at one row in would still refuse; under the ruling it reads every row with the leaf's text.
- One row in out of one row is every row in, refused as one value "(1 of 1 observations are in it)", as before.
  Alternative: read as too few, the same one row read under the words "read because too few are in the
  likelihood", which says something false of a one-row store.
- The text mirrors the no-row one: "(1 of n has a positive weight)". Alternative: "(1 of n observations is in it)",
  true also where a mask given by hand, not a weight, leaves the rest out (a logistic or aft sampler through
  `new(active = )`); the no-row text says "have a positive weight" there too.
- A gp leaf whose lengthscale is stated reads every row at one row in, as at no row in and as gaussian's one row
  in under dec-B479. This is the one call the ruling's "no fit that is created today moves" does not cover: such a
  leaf under probit, ordinal or hazard.probit at creation or a deriving `setData`, or under gaussian at a `setData`
  that keeps or takes the range, was created on the placeholder (the one value, divided by 1, silently) and is now
  centered and scaled over every row with the warning. Alternative: keep the placeholder there, so no accepted call
  moves, at the cost of a stated-lengthscale gp reading every row at no row in and its one row at one row in under
  these families, and every row at one row in under gaussian's creation.
- NEWS keeps its entry, restated to stay true: the covariates are centered and scaled over the rows in, over every
  row where fewer than two are in, and a covariate the rows read hold at one value is refused. xbart's help needs
  no change, its folds reaching one row in only under gaussian, whose response falls back.

Tests. test-leaf-standardization-rows.R: the refusal pin ("a family that reads no response ... refuses one row in")
restated to the one-row text; gaussian and Student-t keep dec-B479's text, unchanged, through `dbarts()` and gaussian
through `bart()`; a gaussian `setData` keeping the range reads every row in the leaf's words, the range kept; for
probit, ordinal and hazard.probit, linear and gp, at creation and at a deriving `setData`, the text, the record and
the reader over every row, the all-rows flag and finite sweeps; probit through `bart()` and `new()` by hand; logistic
and aft by hand with a mask; a stated-lengthscale gp; the one row in missing the covariate; two covariates named; a
column constant over every row refused; two rows in read, and refused at one value; a one-row store refused as one
value. tests/cpp `testLeafStandardizationRows`: one row in by the creation mask under gaussian, probit, nbinom and
logistic, linear and gp; a store view's sampler under probit, as a fold builds it; a probit replacement by its mask
and a gaussian one by weight keeping the range; two rows in read; a re-derivation at one row in keeping every
covariate.

What moves, by a two-library compare against c41f25c0 over 292 fixtures (every family, linear, two covariates,
gp and a stated-lengthscale gp, every row, half, two, one and no row in, at creation, a deriving `setData`, a
gaussian `setData` keeping the range, the row setters' re-derivations and `xbart`): 264 identical in record,
reader, warnings and draws, and the 28 that move all have one row in: 22 refused and now created, and 6 the
stated-lengthscale corner above. Gaussian, Student-t and nbinom at one row in, every re-derivation and every
`xbart` fit are identical.

Size: package code 40 lines added and 33 removed (engine 30 and 27, most of them comments and the rename; R 10 and
6); tinytest 202 and 22; tests/cpp 104 and 1; help, NEWS and rule 9 16 and 13.

Gates, on the package code of 0b940941 from `--preclean` installs in a private library: tests/cpp 412 ok plain and
under ASAN and UBSAN, no diagnostic; `tinytest::test_package` in one process, 268 files, 30762 results, none
failed (test-leaf-standardization-rows.R 566); lint (0 lints), air, check-rc-codoc (48 methods), check-win-drift,
check-doc-freshness and verify-anchors clean, NEWS parsing at 227 entries; on a `--enable-reference-build` install
the four snapshot files (14, 3, 7 and 3, none failed) and the bitwise compares, 55 of 55 identical against
equivalence-32418b6e (`--strict-coverage`, 0 skipped), 19 of 19 against bcf-equivalence-20eed7ed and 11 of 11
against multinomial-equivalence-80b1c8d4, no "max |z|" line, and both cross-host compares; the 43 commands of
exact-gates.yaml at `quick` on the shipped build, all passing (linear-exact.R among them); `R CMD check --as-cran
--no-manual` from a tarball built from a clean export, the Date field's NOTE alone. Not run: the R-loaded path
under ASAN (no new numerics: one row in reaches the every-row derivation no row in already runs), stan4bart, the
x86-64 leg.

Mutants, each alone, an engine one installed with `--preclean` and run with tests/cpp; failures in
test-leaf-standardization-rows.R of 566: one row in not marked (the refusal as before), stops at probit's creation,
tests/cpp aborts on the uncaught refusal; the one-row-store guard dropped, 1; two rows in marked too, 3, tests/cpp
1; the re-derivation deriving where too few are in, 4, tests/cpp 2; the linear leaf's all-rows flag dropped, 40,
tests/cpp 8; the refusal's count not 0 where every row was read, 2; the creation mask ignored, 47, tests/cpp 14;
`setData`'s weights ignored, 26, tests/cpp 1; the one-row text dropped, 19; the one-row text at no row in, 5; the
joint text under every family, 20; `setData` counting the rows in without its weights, 3. All twelve killed,
two by tinytest alone (the guard and the count, which tests/cpp does not reach).
