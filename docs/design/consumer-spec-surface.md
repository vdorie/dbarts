# Consumer specification surface (dbartsSpec): design

Status: PROPOSED 2026-07-25. Plan: docs/plans/consumer-spec-surface.md. Adds one
exported function, `dbartsSpec()`, that resolves a `(control, model, data)` triple
plus its family token exactly as `dbarts()` does, without constructing a sampler.
It closes the last gap in the public path to the flat C API: a `LinkingTo: dbarts`
consumer can already build a `dbartsData` and call `dbarts_sampler_create`, but has
no supported way to produce the `dbartsModel` and the resolved `dbartsControl` that
call expects. rng class: BYTE-neutral (a pure factoring of `dbarts()`'s existing
resolution; no draw-law change, no engine code).

## 1. The problem

Every feature the bartcore program added is configured through attributes parked on
the control and model SEXPs plus a family string: `bartcore.n.categories` (ordinal),
`bartcore.dispersion` (nbinom), `bartcore.survival` (aft), `bartcore.variance`
(heteroscedastic), `bartcore.hazard.periods`, `monotone`, `resid.df`,
`interaction.max.order`, `interaction.forbidden`, `block.of.column`,
`block.tree.counts`. All three SEXPs cross `dbarts_sampler_create(control, model,
data, family)` untouched, so the flat C API can in principle reach the whole
single-forest feature set.

In practice it cannot, because ATTACHING those attributes is unexported logic living
in the body of `dbarts()`. The only consumer to try (stan4bart) resolved it by
hand-constructing `new("dbartsModel", ...)` positionally and calling
`dbarts:::parsePriors` through a `quoteInNamespace` shim. That arrangement has
already broken once - the 1.0 port had to chase a `parsePriors` signature change -
and it reaches nothing beyond the priors: family, monotone, interactions, blocks,
variance, and missingness are all unavailable to it.

The cost is paid twice. Consumers cannot use features that are otherwise complete,
and dbarts cannot change an internal signature without breaking a package that had
no supported alternative.

## 2. What already resolves where

`dbarts()` is the single choke point: `bart2()` and `rbart_vi()` both reach the
engine by `redirectCall`ing into it, and only the multi-forest hosts (BCF,
multinomial, hurdle) build samplers another way. Its body divides cleanly at the
point where the response is materialized:

- INGESTION (formula/matrix dispatch, `Surv` parsing, discrete-time hazard
  expansion, the `dbartsData` call) - formula-interface concerns, upstream of the
  triple.
- RESOLUTION (family from the response shape, ordinal recoding, weight policy,
  the sigma estimate, offset normalization, prior parsing, model assembly, every
  attribute above) - the part a consumer with its own design matrix needs and
  cannot reach.

`dbartsData()` is already exported, so ingestion is already public for consumers
that want it. Only resolution is missing.

## 3. Design

Move resolution out of `dbarts()` into ONE internal, `resolveSamplerSpec()` (R/spec.R),
and export a wrapper over it. No logic is duplicated: `dbarts()` and `dbartsSpec()`
call the same function, so a family can never resolve two ways.

One internal rather than three (response / priors / model assembly, the obvious
split) because the three stages share half a dozen locals - `fixedUnitScale` feeds
the resid.prior override, which feeds `node.scale`; `monotoneDirections` feeds both
the prior parse and the model attributes - and threading those through three
signatures buys nothing but an opportunity to thread one wrong. The move is then a
single contiguous block of code, which is what makes it verifiable against the
bitwise equivalence gate.

The prior parse stays call-shaped inside it: the prior vocabulary is NSE
(`normal(chi(1.5))` must resolve in dbarts's vocabulary regardless of what the
caller has attached), so each entry point hands over its own `match.call()` and
formals and `parsePriors` evaluates the argument expressions in the supplied
environment. That is exactly the shim consumers were writing by hand.

`dbartsSpec()` takes a ready `dbartsData` and returns
`list(control, model, data, family)`.

## 4. Why the family token is returned rather than left on the model

`model@family` always carries the resolved token, and the R-side sampler constructor
reads it there. A flat-C-API consumer must pass it as the fourth argument to
`dbarts_sampler_create` instead, where the empty string means "dispatch on the
response shape". That default is right for gaussian, probit, ordinal, and nbinom -
each is inferable from the response coding or a control attribute the bridge reads -
but WRONG for `aft` and `logistic`, which are indistinguishable from gaussian and
probit by shape and are silently mis-fit if the string is dropped. Returning the
resolved token in the list makes the correct call the obvious one, and the
documentation names the trap.

## 5. Alternatives rejected

- A NARROWER export (parsePriors alone, made public). Fixes stan4bart's immediate
  breakage and nothing else: the attribute wiring stays unreachable, so every future
  family costs a consumer-side commit that duplicates a dbarts internal.
- A HELPER THAT DUPLICATES the attribute wiring rather than sharing it. Drifts by
  construction; the ordinal/nbinom/aft resolution order is subtle (fixedUnitScale
  feeds the resid.prior override, which feeds node.scale) and would be re-derived
  wrong.
- EXPORTING A SAMPLER-BUILDING function instead. Does not serve the motivating
  consumer, which holds its sampler C-side through the flat API and needs the
  specification without an R-side instantiation.

## 6. Timing

The export is the only piece of this that is time-sensitive. Appending a
`dbarts_results` field or adding a BCF entry point are both additive-MINOR by the
ABI's own rules and cost a lockstep release whenever they happen. An export that
consumers should be coding against, by contrast, wants to exist BEFORE those
consumers ship against a frozen 1.0 surface with the `:::` workaround baked in.
