Third whole-branch review of dbarts bartcore (R package, C++ BART engine), shared brief.

Setting. Release 1.0-0 is near. Two earlier whole-branch reviews ran (the second ended at 7ad0bbea,
2026-08-25); since then 1,294 commits landed (~27k lines in src/ R/ man/), reviewed slice by slice, and
slice reviews have found real defects in every slice. Your job is to find the defects that slipped through,
in your lens only. Assume they exist; try to refute correctness, not to confirm it.

Ground rules.
- READ-ONLY. Pinned code: .claude/worktrees/review3 (detached at 01dee4b4).
  Read its CLAUDE.local.md first (architecture pointers, build/test notes). Do not edit, commit, push,
  or install into it; do not touch any other checkout. Spawn no sub-agents.
- A library built from that tree: R_LIBS=<scratchpad>/r3-lib
  (prefix every R call). Released 0.9-x is on branch main (git -C <tree> show main:<path>); to probe old
  behaviour, CRAN's dbarts may be installed into your own scratch lib if useful.
- Scratch files only under <scratchpad>/r3-<lens>-*.
  Run long commands in the FOREGROUND; never a pgrep wait loop; kill any background shell before you finish.
  Threads: n.threads <= 2 in any fit you run.
- Known items are not findings: check root TODO and docs/decisions.md (ledger) before reporting; if a
  finding is already filed or ruled, drop it or mark it KNOWN with the TODO/ledger id only if you have new
  evidence it is worse than recorded. docs/ explains intended behaviour: check it before calling something
  accidental.
- Evidence. Every finding carries a reproducing probe (the R or C++ snippet and its actual output) or, when a
  probe is infeasible, an exact code argument naming file and symbol. Say why the existing gates (tinytest,
  tests/cpp, the exact gates, equivalence) did not catch it. No speculation-only findings.
- Severity: BLOCKER (silently wrong results, crash, memory corruption, data loss); MAJOR (wrong or surprising
  for a user, a regression against 0.9-x not documented as intended, a documented claim that is false);
  MINOR (edge-case inconsistency, misleading message). Skip style.

Output. Write your findings to scratch/review-3/<lens>.md (only that file):
a header with what you covered and what you did not, then one entry per finding:
  ID (<lens>-NN), severity, location (file + symbol), claim (one sentence), probe + output, why gates missed it,
  suggested fix (one or two sentences). Then a short "checked and found correct" list.
Final report to the orchestrator: <= 20 lines - counts by severity, the BLOCKER/MAJOR one-liners, the file path.
Budget: about 2 hours of work; at 3 hours stop and write up what you have.
