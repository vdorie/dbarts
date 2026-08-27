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

- Scoped down (VD, 2026-07-06): nothing from the 1.0 branch is
  released, so no cross-version guarantee apparatus and nothing for
  binary formats - at most a read helper for serialized R state
  objects. Realization: version-tag new states; setState/re-creation
  refuses unversioned or mismatched states cleanly ("re-fit or use the
  release that wrote it"). No migration shims. The shared version
  scheme for the format-bumping items stands.
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

## Landing note (2026-07-06, d3f992a)

Landed at the scoped-down level: storeState stamps formatVersion (1)
and the writing package version; setState refuses mismatched or
unversioned states before structural parsing, naming both versions
(unversioned reads as version 0, "written by dbarts unknown"); the
guarantee paragraph sits in man/dbartsSampler-class.Rd and
public-surface.md section 2. Refusal-without-mutation and
bitwise-restore covered in test-bartcore.R. Gates: tinytest 2464 ok,
component tests pass, equivalence exact 18/18 vs af04d0c, lint clean.
No deviations; bench-sampler not applicable (cold path).
