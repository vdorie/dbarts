# rchk-ci

agent: sonnet
rng: neutral
budget: one workflow + fixes for findings

## Goal

rchk (the PROTECT-balance static analyzer CRAN itself runs) checks the
bridge on a schedule, and its current findings are fixed or annotated.

## Context

- The bridge hand-manages a protection scheme (the PROT_* fixed slots
  in src/R_interface_bartcore.cpp:35-50, retain/release around
  transactional setters) - exactly rchk's bug class, and
  data-ownership will churn this code.
- Available harness: the kalibera/rchk docker image; community actions
  exist, or a container job running rchk's rcheck against the built
  package.

## Constraints

- Scheduled + manual dispatch (rchk is slow; not per-push).
- False positives get a tracked annotation (rchk supports suppression
  comments sparingly), not silence.
- Out of scope: rewriting protection style wholesale (data-ownership
  owns the structural change).

## Steps

1. .github/workflows/rchk.yaml: container job, build, run rchk,
   fail on new findings (baseline file if the first run is noisy).
2. Triage the first run: fix real imbalances; record false positives.
3. One line in docs/plans/README.md's gate list: bridge-touching items
   note "rchk on next scheduled run" in review.

## Verification

- The workflow runs green via workflow_dispatch; a deliberately
  unprotected allocation on a scratch branch is caught (then
  reverted).

## Landing note (2026-07-07, 3b30cf5)

Landed the workflow (weekly + dispatch; builds the tarball, runs
kalibera/rchk in a container, fails on [PB]/[UP]/maacheck suspicious-
call tags, tolerates the benign too-many-states bailout) and the
README review-gate line. Step-2 triage DEFERRED: Docker Desktop's
hub proxy on this machine blackholes pulls of non-cached images
(hello-world included; registry auth itself answers), so no local
rchk output was obtainable - first triage rides the first dispatched
run after a push. One unverified candidate noted for that run:
bartcore_run's namesExpr passed unprotected to Rf_setAttrib
(R_interface_bartcore.cpp:~1434), matching rchk's documented
setAttrib pattern; deliberately not fixed without tool confirmation.
