# setstate-mixed-chains: setState refuses chains saved under different leaf priors

Status: PLANNED 2026-10-02 under dec-B191 in [decisions.md](../decisions.md). Starts after
[nbinom-dispersion-name.md](nbinom-dispersion-name.md) lands: both edit the bridge's state code.

agent: opus implementer, one; opus reviewer.
rng: NEUTRAL. A state that installed and ran correctly before installs and runs bitwise as before; only a state
whose chains disagree is newly refused.
window: pre-release, before the 1.0-0 merge.
budget: ~250 lines (bridge ~60, R ~30, manual ~20, tinytest ~120, records ~20). Plans have run 1.5-2x low.

## Goal

A sampler's chains always share one leaf prior per forest and one response transform. `setState` refuses a state
whose chains disagree on any of them and names the one that differs, so
[`getLeafPrior`](../../man/dbartsSampler-class.Rd) reports real values for them and never an NA that means "the
chains disagree".

## Context

- The reader collapses per-chain values in [`reportLeafPrior`](../../R/dbarts.R): its `shared` helper returns NA
  when the chains' values differ. The per-chain table comes from
  [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp).
- A saved chain carries its own response transform (the `fit.scale` block) and, per forest, its k and leaf scale;
  the install paths that read them are in the same bridge file.
- Probe, 2026-10-02: a chain saved from a sampler made with `normal(sd = 1)` spliced into one made with
  `normal(sd = 0.3)` installs, the sampler runs, the reader reports the sd and anchor as NA, and writing the read
  value back is refused with a message about disagreeing chains.
- The reader has a second NA, on a forest whose scale a calibration map sets, after a state install brings a
  calibration the map did not derive. There the chains agree with each other. It is out of scope.

## Constraints

- The check covers exactly what the reader reports once for all chains: per forest the anchor, the fixed spread
  and the prior mean, and the response scale and shift. Use the quantities the reader's `shared` helper collapses,
  so the two cannot drift apart.
- A drawn k differs by chain and must still install. So must chains spliced from separate samplers built with
  the same prior on the same data, and every state the test suite and the consumers install today.
- Equality is exact. Chains under one prior and one response compute these values by the same arithmetic.
- The refusal happens before any chain is installed: a refused `setState` leaves the sampler as it was.
- Every path that installs chains from a saved state is covered, not only the `setState` method; name each one
  in the Landing note.
- The message names the quantity and the forest, and says the chains were saved under different leaf priors or
  responses. One message shape, pinned by pattern in the tests.
- Afterwards the reader's NA-on-disagreement branch is unreachable for these quantities: remove it, and narrow
  the wording of the write refusal and of the manual so they describe only the calibration case that remains.
- No change to what `storeState` writes, to either state-format constant, or to the shipped C header.
- Out of scope: the foreign-calibration NA, the anchor's meaning under an sd-named prior (dec-A124, not ruled).

## Steps

1. Bridge: compare the incoming chains before installing; refuse with the message. Find every install path.
2. R: drop the reader's disagreement branch; narrow the write refusal's message.
3. Manual: `getLeafPrior`'s and `setState`'s entries in `man/dbartsSampler-class.Rd`, and the reference-class
   docstrings they mirror.
4. Tests: a mixed spread, a mixed response, a mixed forest on a multi-forest sampler, each refused with the
   sampler unchanged afterwards (draws bitwise equal to an untouched twin); a same-prior splice across two
   samplers, a drawn-k state and a full install from a sampler with a different prior, each accepted.
5. Records: the Landing note here, the index row, the TODO item.

## Verification

Against a private library, installed with `--preclean`:

- `cd tests/cpp && make && ./test_bartcore` passes.
- `tinytest::test_package("dbarts")` passes with no new warning; the new tests fail when the check is removed
  (run that mutation once and say so).
- stan4bart's and bartCause's suites pass against a fresh private-library chain built on this tip.
- `lintr::lint_package()`, `air format --check .`, `tools/check-rc-codoc.R`, `tools/check-win-drift.R` and
  `tools/check-doc-freshness.R` pass, each on its own exit status.
