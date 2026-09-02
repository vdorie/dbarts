# family-on-model

agent: sonnet
rng: neutral
window: pre-release (the S4 slot layout freezes at 1.0-0)
budget: ~200 lines

## Goal

The response family lives with the rest of the statistical
specification, or a recorded decision says why it stays on the control
object. Today it sits on dbartsControl purely by analogy to classic's
Control::responseIsBinary.

## Decision

Question: move family from dbartsControl ([[R/A_class.R:189@07230c35]]) to
dbartsModel? Recommendation: move. dbartsModel holds tree.prior,
node.prior, node.hyperprior, resid.prior - everything else that defines
the model - and the engine's own decomposition makes ResponseModel the
top-level statistical concept (src/bartcore/model.hpp). The slot was
placed by analogy during the rewrite (commit 0a02263 "the way it
preserves binary"), not by design. Counter-evidence that would change
the call: family also drives pointer re-creation after save/load and
setControl preservation - if the migration of those paths turns out to
touch the state format, keep it and record why. Signed off (VD,
2026-07-06): move.

## Constraints

- bart/bart2/dbarts signatures do not change; only where the resolved
  value is stored.
- Saved fits from pre-release dev versions are not a compat concern
  (nothing is released); do not add migration shims.
- Out of scope: adding families; rbart_vi's probit-only refusal.

## Steps

1. Move the slot and its validity checks to dbartsModel; family
   resolution in dbarts()/bart2 targets the model object.
2. Chase consumers: the bridge's family argument, pointer re-creation
   after save/load, setControl's preservation (moves to setModel),
   control prototype, show methods.
3. Update man/dbarts.Rd, man/bart.Rd, man/dbartsControl.Rd,
   man/dbartsModel.Rd.
4. If the decision lands "keep": instead add one paragraph to
   docs/design/public-surface.md section 3 recording the reason.

## Verification

- Full tinytest suite (test-bartcore.R exercises family dispatch and
  save/load re-creation).
- R CMD check codoc clean.
- Equivalence exact (neutral).

## Landing note (2026-07-06, 3c02121)

Landed: slot, prototype, and validity moved to dbartsModel; family
resolution targets the model in dbarts()/xbart(); creation and both
re-creation paths pass model@family to the bridge; setControl no
longer carries it and setModel preserves it (tested both ways);
consumers chased through bart/rbart/xbart result packaging and the
xbart loss dispatch. Plan's man-page step proved a no-op: no Rd file
documents slots and dbartsControl.Rd never mentioned family; the
dbarts.Rd family item describes the argument, not storage. No show
methods exist. Implementer hit a usage limit after the substance;
reviewer ran the gates and landed: tinytest 2465 ok, codoc clean,
equivalence exact vs af04d0c, lint zero; component tests vacuously
unchanged (no C/C++ delta). State format untouched.
