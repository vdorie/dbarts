# docs/ orientation

A map of `docs/`: where the design records, plans and process live. It
records no decision itself.

## Layout

- `docs/architecture.md` - how the engine and its R and C layers work
  today. Start here.
- `docs/design/` - one doc per feature or investigation: the problem, the
  options considered, the decision and, for landed work, the landing notes.
  [`docs/design/INDEX.md`](design/INDEX.md) lists every design doc by theme
  with its status and a one-line summary.
- `docs/plans/` - implementation plans, one per backlog item or feature,
  each closing with a Status or Landing note.
  [`docs/plans/INDEX.md`](plans/INDEX.md) lists them the same way.
  `docs/plans/archive/` holds landed or closed plans that only record how
  the branch got here; the index lists them in its last section.
- `docs/plans/README.md` - the plan process: the plan template, the
  citation grammar, the RNG gate classes and the gates each requires,
  landing, the reviewer checklist, gate hygiene, and what CI runs.
- `docs/decisions.md` - the decision register: one prose entry per design
  decision the bartcore branch carries, with who made it and on what
  evidence, and a Marked line the maintainer writes on. A claim that the
  maintainer decided something cites an entry id here; a claim with no
  entry has no attribution.
- `docs/plans/bartcore-review-tour.md` - the case for merging bartcore into
  main, with a reading order for the code.
- `docs/plans/bartcore-landing/` - the landing memo (`memo.md`) and the
  registers behind it: changes against main (`changes.md`), completion
  status (`completion.md`), a doc inventory for triage (`doc-inventory.md`),
  neutral evaluations of agent-made decisions (`evaluations.md`), the rubric
  for documents the maintainer reads (`rubric.md`), and the check sheets for
  the memo and the register (`memo-check.md`, `ledger-check.md`).
- `docs/plans/review-2026-08-24/` - the working records of the 2026-08-24
  whole-branch review: findings, evidence, and the scripts and logs behind
  them. Historical; TODO cites a few of its memos.

Outside `docs/`: the repo-root `TODO` is the open backlog, forward-facing
only, and most entries name a plan; `benchmarks/README.md` describes the
timing, equivalence and exact-posterior harnesses.

To choose what to work on, read `TODO`. To write a plan, read
`docs/plans/README.md`. To find a doc, start at the two INDEX files.

## Conventions

- **Status.** Most docs state their status near the top, either as a
  `Status: <word>, <date>` line (usual in `docs/design/`) or as a
  `## Status` section (usual in `docs/plans/`). Standing references such as
  `docs/design/data-store.md` and `docs/design/kernel-vocabulary.md` carry
  none; `docs/design/prior-defaults.md` carries `Status: reference,
  current` and is updated in place when a default changes. For a doc with
  neither, the INDEX files take the status from its opening text.
- **Authoritative sections.** Some long design docs end with a dated section
  marked authoritative that overrides the earlier ones (`hurdle.md`
  section 13, `forest-ranef-interweaving.md` section 9). Read it first.
- **Design and plan pairs.** `docs/design/X.md` usually pairs with
  `docs/plans/X.md`. Standing references and NO-GO or parked
  investigations may have no plan. The pairs that do not share a name:
  - `docs/design/monotone.md` with `docs/plans/archive/monotone-bart.md`,
    and `docs/design/ordinal.md` with
    `docs/plans/archive/ordinal-outcomes.md`.
  - `docs/design/memory-footprint.md` with
    `docs/plans/memory-footprint-audit.md`: the model is in the design
    doc, the measurements that check it in the plan.
  - `docs/plans/archive/within-chain-threading.md` records the NO-GO; the
    analysis is section 8 of `docs/design/within-chain-threading.md`.
  - `docs/plans/x86-simd-plan.md` is a measurement memo, the x86 companion
    to `docs/plans/simd-survey.md`; `docs/plans/archive/x86-simd.md` is the
    closed action plan on the same gaps.
  - The C API growth work is spread over three distinct docs,
    `docs/plans/archive/c-api-growth.md`, `capi-callbacks.md` and
    `capi-dispatch-table.md`.
