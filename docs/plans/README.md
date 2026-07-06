# Plan process

The TODO at the repo root is an unordered backlog; every item names a
plan file here. This document is the contract between the planning
agent, the implementing agents, and the reviewer.

## Roles

- Fable (orchestrator): writes and owns plans, routes items, reviews
  diffs, updates the TODO. Plans are written against the code, not from
  memory: every anchor in a plan was read before it was cited.
- Opus (implementer): engine internals (src/bartcore/), numerics,
  statistics, MCMC-kernel changes, C API design.
- Sonnet (implementer): mechanical edits, deletions, build system,
  R surface, tests, documentation.

Routing is the plan's `agent:` field. When in doubt: could a mistake
change the posterior? Opus. Otherwise Sonnet.

## Plan format

One file per TODO item, `<item>.md`, hard cap 80 lines. Front block:
`agent:`, `rng:`, `window:` (if any), `budget:` (expected diff size).
Sections:

- Goal: two or three sentences; what is true after the item lands.
- Context: file:line anchors and design-doc pointers; no narrative a
  pointer can replace.
- Decision (decision-gated items only): the question, a recommendation,
  and what evidence would change it. VD signs off before implementation.
- Constraints: gates, contract freezes, explicit out-of-scope list.
- Steps: numbered; each independently verifiable.
- Verification: exact commands and expected outcomes.

## RNG classes and their gates

- neutral: draws unchanged.
  Gates: tests/cpp component tests; full tinytest suite.
- shifting: draws change, the posterior does not.
  Gates: the above, plus regenerate RNG-locked snapshots by replaying
  whole test files, re-record the equivalence baseline, and pass the
  statistical (z) mode against the previous baseline.
- posterior-changing: the stationary distribution or a default changes.
  Gates: all of the above, plus the exact-posterior gates
  (benchmarks/R/logistic-reference.R, categorical-exact.R) and a design
  note in docs/design/.

Hot-path changes of any class additionally need bench-sampler.R compare
on a quiet machine (maintainer-run; never concurrent with other load).

## Implementation protocol

- One item per agent, in an isolated worktree. The prompt is the plan
  file path, CLAUDE.local.md, and the report format below. The plan is
  the spec; do not restate it.
- After any C/C++ change: R CMD INSTALL . before tinytest (--preclean
  after header or virtual changes).
- Stop conditions: a step fails twice; the diff exceeds 1.5x budget; a
  needed change is out of scope. Report and stop; do not improvise.

## Brevity (binding; Opus especially)

- Final report <= 20 lines: files touched, gate results (pass/fail
  lines verbatim), deviations from plan. Nothing else - no methodology,
  no plan recap, no prose about what the code now does.
- New comments only for constraints the code cannot show. The diff's
  comment-line delta must not exceed its code-line delta; the reviewer
  checks this mechanically.
- No references in code or comments to this process, the plan, or the
  conversation. "The classic engine did X" only when compatibility is
  itself the constraint.
- Match the surrounding file's comment density and idiom (Doxygen, LLVM
  house style, ASCII).

## Review (Fable)

1. Diff vs plan: every step landed; nothing outside scope.
2. Gates: re-run or verify transcripts; RNG class honored.
3. Brevity: comment/code delta, report length, no narration.
4. Readability (standing): could a maintainer with no access to any
   agent conversation understand the diff from code and comments alone?
   Reject change-log comments, restated behavior, hedging prose.
5. At most two fix rounds by message to the same agent; after that, fix
   directly or withdraw and replan.
