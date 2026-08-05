# python-bindings

status: RESEARCH-OPEN, decision-gated (as of 2026-08-04) - no spike run
        recorded; horizon item, likely out-of-repo (TODO python-bindings
        carries the two non-R-host contracts a binding must take over)
agent: opus (prototype); decision-gated, likely out-of-repo
rng: neutral (no engine change intended)
budget: memo + spike; product work lives in its own repository

## Goal

A recorded feasibility answer on exposing the bartcore engine to
Python, with a spike proving or refuting "prototype-scale after
kernel-cleanups".

## Context

- What blocks it today: misc.a prints via Rprintf (kernel-cleanups
  removes this); external.a's RNG seeds from R in the bridge but the
  engine takes ext_rng handles (standalone seeding exists - the
  equivalence rshim and tests/cpp already run engine code without a
  live R session, via rshim stubs); the C API takes R spec objects at
  creation, but the facade underneath takes plain C++ structs - a
  binding would construct those directly and own its ingestion
  (no dbartsData; NumPy arrays in, quantization via ColumnStore's
  existing build path).
- Audience evidence: pymc-bart's traction; no fast embeddable BART
  core exists in Python with a mutation surface (the IRT/BCF embedding
  story is unique).
- Cost honesty: a binding forks the surface-maintenance burden
  (docs, tests, releases) - the memo must weigh that against demand;
  out-of-repo (a separate package vendoring or submoduling the
  headers) keeps CRAN dbarts unaffected.

## Constraints

- Blocked on kernel-cleanups.
- The spike must not modify dbarts sources beyond what kernel-cleanups
  already planned; if it needs more, that IS the finding - record it.
- Out of scope: packaging, wheels, API design polish.

## Steps

1. Spike (scratch repo): pybind11 wrapper over SamplerBase facade -
   construct from arrays, run, return train fits; measure against
   pymc-bart on one synthetic problem.
2. Memo: what it took, what dbarts-side seams would need formalizing
   (a stable facade header contract), maintenance model; VD decides.

## Verification

- The spike fits and matches dbarts-in-R draws under a shared seed
  protocol (same ext_rng streams); memo recorded; TODO updated.
