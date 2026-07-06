# state-format-policy

agent: sonnet
rng: neutral
window: pre-release (the first release ships the first states users keep)
budget: ~120 lines + docs

## Goal

Saved sampler states carry an explicit format version and package
provenance; loading an incompatible state refuses cleanly with a
message saying what to do; the guarantee level is documented where
users will find it.

## Context

- States are "opaque and engine-specific"
  (docs/design/public-surface.md section 2) - deliberately unspecified,
  which silently means saved fits may not survive upgrades and the
  failure mode today would be a validation error at best.
- Three in-flight items bump the format (flat-format-v2,
  state-continuation, forest-split-bcf); they need ONE shared version
  scheme, defined here first (cheap to do ahead of them).
- getPointer re-creates samplers from stored state after save/load;
  that path is where the check lives.

## Constraints

- Policy recommendation: states load within a release series
  (same format version); cross-version loads refuse with the version
  pair named and "re-fit or use the release that wrote it". No
  migration shims unless a future item buys one deliberately.
- The version check must be cheap and precede any structural
  validation (refuse before touching content).
- Out of scope: the format changes themselves (owned by the items
  above); serializing the standalone data handle (open per
  public-surface.md section 5).

## Steps

1. Add formatVersion + writing package version to the state object;
   setState/re-creation checks version first and errors with the
   documented message.
2. Document the guarantee: man page for the sampler state surface +
   one paragraph in public-surface.md section 2's landing record.
3. tinytest: a state with a doctored version refuses cleanly and
   mutates nothing; a current state round-trips.

## Verification

- Full tinytest incl. the refusal test; component tests unchanged;
  equivalence exact.
