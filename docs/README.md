# docs/ orientation

This is the wayfinding map for `docs/`. It does not itself record any
decision - it explains where decisions live and how to find one.

## Layout

- `docs/architecture.md` - the current-state orientation doc for the
  bartcore engine (layering, concepts, where to look). Start here for "how
  does the engine work today."
- `docs/design/` - design proposals, rationale, and landing/decision
  records. One doc per feature or investigation: the problem, the forks
  considered, the decision, and (for landed work) the landing notes. ~37
  files. Manifest: `docs/design/INDEX.md`.
- `docs/plans/` - implementation plans: one file per TODO backlog item or
  landed feature, with Goal/Context/Steps/Verification and a closing
  Status or Landing note. ~120 files (plus the process doc below).
  Manifest: `docs/plans/INDEX.md`.
- `docs/plans/README.md` - the process/contract doc: roles (who plans, who
  implements, who reviews), the plan-file template, RNG gate classes and
  the gates each requires, the brevity rubric, the review checklist. It is
  NOT an index and does not enumerate the directory's contents - that job
  belongs to `docs/plans/INDEX.md`.
- repo-root `TODO` - the live, unordered backlog of OPEN work. Every entry
  names a plan file in `docs/plans/`. Completed work is deliberately
  trimmed out of the TODO once landed; its record lives in git history and
  in the design/plan docs themselves, not here.

## Three navigation surfaces, three different jobs

Do not expect any one of these to do another's job:

1. **The root TODO** - what's still OPEN, right now. Forward-facing only;
   landed items are removed once closed.
2. **`docs/plans/README.md`** - HOW plans get written, gated, and
   reviewed. Process, not content. Has no file-by-file listing.
3. **`docs/design/INDEX.md` and `docs/plans/INDEX.md`** (this pair) - the
   COMPLETE map of every doc in each directory, landed and closed included,
   with a current status and a one-line purpose. Use these to answer "does
   X already have a design doc / plan" or "what happened to Y."

If you're deciding what to work on next, read the TODO. If you're writing
a plan, read `docs/plans/README.md`. If you're trying to find or reconcile
a doc, start at the two INDEX files.

## Doc conventions

- **Status header.** Most docs carry a status marker near the top, in one
  of two forms - both in active use, normalize your reading of either:
  - `Status: <WORD>, <date>...` as a top-of-file line (common in
    `docs/design/`).
  - `## Status` as a section heading, with the verdict in the prose below
    it (common in `docs/plans/`).
  A few standing-reference docs (not proposals - e.g. `docs/design/data-store.md`,
  `docs/design/kernel-vocabulary.md`) carry no status line by design; their
  purpose is self-evident from the opening paragraph
  (`docs/design/prior-defaults.md` instead carries `Status: reference,
  current`, updated in place when a default changes). Where a doc has neither form and isn't a standing
  reference, that's a real wayfinding gap - the two INDEX files flag it.
- **Authoritative resolutions / dated supersedes.** The best docs in this
  corpus close with a dated, explicitly-"authoritative" section that
  states which earlier sections it overrides (e.g. `hurdle.md` section 13,
  `forest-ranef-interweaving.md` section 9). This is what keeps a dense
  document navigable without a rewrite: state the final answer once, dated,
  scoped to what it overrides - don't require the reader to un-believe
  earlier sections themselves.
- **Design <-> plan pairing.** The common case is `docs/design/X.md` paired
  1:1 with `docs/plans/X.md` (same basename). Several features instead pair
  a design doc with a cluster of differently-named plan files (see
  collisions below), and a few design docs (standing references, or
  NO-GO/parked investigations) have no plan file at all - that's expected,
  not a gap.

## Known name collisions and non-obvious pairings

These are real, intentional, and not to be renamed - just navigate them:

- `docs/plans/c-api-growth.md` vs. `docs/plans/capi-callbacks.md` /
  `docs/plans/capi-dispatch-table.md` - same ABI-growth subsystem, two
  spellings of the same prefix (`c-api` vs `capi`). Three distinct
  stages of one arc, not redundant (only capi-dispatch-table.md
  explicitly cites c-api-growth.md).
- `docs/plans/x86-simd.md` vs. `docs/plans/x86-simd-plan.md` - NOT
  duplicates. `x86-simd-plan.md` is a READ-ONLY measurement memo (the
  x86-measured follow-up to `simd-survey.md`); `x86-simd.md` is the
  action-plan-shaped file that resulted from it, now CLOSED/SUPERSEDED.
- `docs/design/monotone.md` <-> `docs/plans/monotone-bart.md` - design doc
  and plan file do not share a basename.
- `docs/design/ordinal.md` <-> `docs/plans/ordinal-outcomes.md` - same
  mismatch pattern.
- `docs/plans/within-chain-threading.md` is a full plan file (goal,
  steps, Landing notes summarizing the NO-GO), but the deep analysis
  lives in `docs/design/within-chain-threading.md` section 8, which its
  status line points at.
- Multi-part programs, three different suffix conventions for "part N of a
  bigger arc": `data-ownership-{1-container,2-ingestion,3-mutation,
  4-views,5-sparse}.md` (numbered-slug), `block-fusion-stage-{a,b}.md`
  (lettered-stage), and `forest-split-bcf.md` (no suffix at all - a second
  "Phase 2" is bolted onto the same file instead of getting its
  own name).

## Index files

- [`docs/design/INDEX.md`](design/INDEX.md) - every design doc, grouped by
  theme, with current status and a one-liner.
- [`docs/plans/INDEX.md`](plans/INDEX.md) - every plan doc, grouped by
  cluster, with current status and a one-liner.
