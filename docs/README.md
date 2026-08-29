# docs/ orientation

This is the wayfinding map for `docs/`. It does not itself record any
decision - it explains where decisions live and how to find one.

## Layout

- `docs/architecture.md` - the current-state orientation doc for the
  bartcore engine (layering, concepts, where to look). Start here for "how
  does the engine work today."
- `docs/design/` - design proposals, rationale, and landing/decision
  records. One doc per feature or investigation: the problem, the forks
  considered, the decision, and (for landed work) the landing notes.
  Manifest: [`docs/design/INDEX.md`](design/INDEX.md) - every design doc,
  grouped by theme, with current status and a one-liner. Use it to answer
  "does X already have a design doc" or "what happened to Y."
- `docs/plans/` - implementation plans: one file per TODO backlog item or
  landed feature, with Goal/Context/Steps/Verification and a closing
  Status or Landing note. Manifest:
  [`docs/plans/INDEX.md`](plans/INDEX.md) - every plan doc, grouped by
  cluster, with current status and a one-liner. Same use as the design
  index, for `docs/plans/`.
- `docs/plans/README.md` - the process/contract doc: roles (who plans, who
  implements, who reviews), the plan-file template, RNG gate classes and
  the gates each requires, the brevity rubric, the review checklist. HOW
  plans get written, gated, and reviewed - not a listing of the
  directory's contents; that is `docs/plans/INDEX.md`'s job.
- repo-root `TODO` - the live, unordered backlog of OPEN work, forward-
  facing only. Most entries name an implementation plan in `docs/plans/`;
  a few record instead in a design doc or a differently-named plan.
  Completed work is trimmed out once landed; its record lives in git
  history and in the design/plan docs themselves, not here.

If you're deciding what to work on next, read the TODO. If you're writing
a plan, read `docs/plans/README.md`. If you're trying to find or reconcile
a doc, start at the two INDEX files above.

## Doc conventions

- **Status header.** Most docs carry a status marker near the top, in one
  of two forms - both forms are in use:
  - `Status: <WORD>, <date>...` as a top-of-file line (common in
    `docs/design/`).
  - `## Status` as a section heading, with the verdict in the prose below
    it (common in `docs/plans/`).
  A few standing-reference docs (not proposals - e.g. `docs/design/data-store.md`,
  `docs/design/kernel-vocabulary.md`) carry no status line by design; their
  purpose is self-evident from the opening paragraph
  (`docs/design/prior-defaults.md` instead carries `Status: reference,
  current`, updated in place when a default changes). Where a doc has
  neither form and isn't a standing reference, the two INDEX files fall
  back to a status drawn from the doc's own opening text rather than
  flagging a gap.
- **Authoritative resolutions.** Long docs close with a dated,
  explicitly-"authoritative" section that overrides earlier ones (example:
  `hurdle.md` section 13, `forest-ranef-interweaving.md` section 9); read
  that section first.
- **Design <-> plan pairing.** The common case is `docs/design/X.md` paired
  1:1 with `docs/plans/X.md` (same basename). Several features instead pair
  a design doc with a cluster of differently-named plan files (see
  collisions below), and a few design docs (standing references, or
  NO-GO/parked investigations) have no plan file at all - that's expected,
  not a gap.

## Known name collisions and non-obvious pairings

These are intentional, not renamed for consistency - navigate them:

- `docs/plans/archive/c-api-growth.md` vs. `docs/plans/archive/capi-callbacks.md` /
  `docs/plans/archive/capi-dispatch-table.md` - same ABI-growth subsystem, two
  spellings of the same prefix (`c-api` vs `capi`). Three distinct
  stages of one arc, not redundant (only capi-dispatch-table.md
  explicitly cites c-api-growth.md).
- `docs/plans/archive/x86-simd.md` vs. `docs/plans/x86-simd-plan.md` - NOT
  duplicates. `x86-simd-plan.md` is a READ-ONLY measurement memo (the
  x86-measured follow-up to `simd-survey.md`); `x86-simd.md` is the
  action-plan-shaped file that resulted from it, now CLOSED/SUPERSEDED.
- `docs/design/monotone.md` <-> `docs/plans/archive/monotone-bart.md` - design doc
  and plan file do not share a basename.
- `docs/design/ordinal.md` <-> `docs/plans/archive/ordinal-outcomes.md` - same
  mismatch pattern.
- `docs/plans/archive/within-chain-threading.md` is a full plan file (goal,
  steps, Landing notes summarizing the NO-GO), but the deep analysis
  lives in `docs/design/within-chain-threading.md` section 8, which its
  status line points at.
