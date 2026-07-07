# bridge-error-path-leaks

agent: opus
rng: neutral (error paths only; success-path draws untouched)
budget: ~150 lines across the bridge entry points

## Goal

No bridge entry point leaks C++ heap on an Rf_error path: the
sanitizers workflow's valgrind leg (its first run, 2026-07-07, run id
28865156756) reports zero definitely-lost blocks.

## Context

- Rf_error longjmps past C++ destructors. The creation parsers fill
  std::vector members of stack-local ParsedControl/ParsedData/
  ParsedModel and then call Rf_error (directly and through the rc_*
  validators) - every deliberate error-path test in the suite leaks
  those buffers. setState already handles this with an error
  accumulator (see its comment); the parsers and setters do not.
- Valgrind-confirmed sites (first leg run): parseData columnTypes
  (R_interface_bartcore.cpp:368) and maxNumCuts (:440) reached from
  createHolder (:1087), createFromHandle (:1208), and setData (:1514);
  scratch vectors in createFromHandle (:1224) and setCutPoints
  (:1865). 15 error contexts, 31 blocks, 1457 bytes definitely lost.
  Treat the listed sites as examples of the class, not the whole of
  it: audit every extern "C" entry point that mixes owning C++ locals
  with Rf_error-capable calls.

## Constraints

- Success paths bitwise-unchanged (equivalence exact gates).
- Two acceptable fix shapes; pick per site, favoring the first:
  1. setState's accumulator: validate into an error message, destroy
     the C++ state, Rf_error at the end (works where control flow is
     simple).
  2. R_UnwindProtect around the owning scope so destructors run
     before the jump continues (works where rc_* validators error
     deep inside; keep the wrapper small and single-purpose).
- Do not weaken any validation or reword existing error messages
  (tests match on them).
- Out of scope: leaks R itself reports (none observed); the PROT_
  retain scheme (data-ownership owns its future).

## Steps

1. Enumerate the class: every bridge entry point whose stack owns
   C++ containers across a possible Rf_error (creation, setData,
   setCutPoints, setModel, setState relatives, predict paths).
2. Fix per the constraint above.
3. Local verification: macOS has no valgrind; verify mechanically
   (e.g. a scoped canary asserting destruction on the error path in
   the component tests, or targeted valgrind in the r-hub container
   if the image is available locally).
4. The authoritative gate is the CI valgrind leg on the next push.

## Verification

- Component tests and full tinytest green; equivalence exact 18/18.
- Next push's valgrind leg: ERROR SUMMARY: 0 errors, definitely lost:
  0 bytes.

## Landing note (2026-07-07, a41bd22)

One file, src/R_interface_bartcore.cpp, net +56 (raw +499/-443, the
churn is closure reindentation). An 18-line unwindProtect helper runs
each owning scope as a heap-held closure under R_UnwindProtect; the
cleanup deletes the closure on both the normal return and the error
jump, so by-value-captured containers destruct instead of leaking.
Ten entry points wrapped (creation, setData, setModel, the predictor
updaters, setCutPoints, getTrees/printTrees); the hot paths
(bartcore_run, per-observation update) free-then-error instead
(shape 1) and stay wrapper-free. Audited-and-excluded: setControl
(scalars only), pointer-only setters/getters/predict, setState
(already accumulator-form). Destructor-on-jump was demonstrated with
a temporary canary (removed). Gates: install, component 87/87,
tinytest 2472/0, equivalence exact 18/18. Pending maintainer gates:
the CI valgrind leg on this push is the authoritative zero-leak
check; bench-sampler's setPredictor scenarios cover the wrapper's
per-call allocation at the next quiet-machine run (expected
negligible: one R alloc per mutation call).

Confirmed: the first post-fix valgrind leg (run 28871332412,
2026-07-07) passed - zero findings; all three sanitizer legs green.
