# revdep-smoke-ci

agent: sonnet
rng: neutral
budget: one workflow

## Goal

Breakage in the two tightest reverse dependencies surfaces monthly,
not at release-sweep time: CI installs dev dbarts and runs R CMD check
on stan4bart and bartCause.

## Context

- The release procedure runs a full manual 23-package sweep; between
  releases nothing watches. stan4bart is the only compiled-boundary
  consumer (dbarts.h; lockstep releases, VD maintains both) and
  bartCause is VD's own downstream - both break first and matter most.
- Sources: stan4bart and bartCause from CRAN (released) plus,
  when their dev branches exist on GitHub, those too - the dev-vs-dev
  pairing is the one that predicts the next lockstep release.

## Constraints

- Monthly cron + manual dispatch; runtime dominated by stan4bart
  compilation (~tens of minutes) - acceptable on a schedule.
- Failures notify (issue creation or workflow failure is enough; VD
  watches the repo).
- Out of scope: the full 23-package sweep (stays manual at release);
  fixing what it finds.

## Steps

1. .github/workflows/revdep-smoke.yaml: install dbarts from checkout,
   fetch stan4bart + bartCause (CRAN tarballs; dev remotes where
   configured), R CMD check both with _R_CHECK_FORCE_SUGGESTS_=false.
2. Matrix note: run on the platform the packages' users have (ubuntu
   is fine; macOS optional later).

## Verification

- workflow_dispatch run is green against current CRAN versions;
  a deliberate dbarts.h signature change on a scratch branch makes
  stan4bart's leg fail (then reverted).
