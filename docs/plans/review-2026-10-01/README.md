# Third whole-branch review, 2026-10-01

The working records of the third whole-branch review of bartcore, run on
2026-10-01 against code pinned at 01dee4b4, 1,294 commits after the second
review ended at 7ad0bbea. Ten reviewers each took one lens: nine read one
area of the package and were told to try to break it, and the tenth compared
posteriors against released 0.9-34. Every finding carries a probe that
reproduces it. Independent verifiers then re-derived each finding of the
nine with their own probes. Alongside ran a mutation leg: planted defects,
to test the tests.

The files are as the review left them, except that line citations are now
links to the code at 01dee4b4 and session paths are gone. Probe scripts
named "scratchpad r3-..." lived in the review session's scratch space and
were not kept; neither was the pinned review3 worktree.

What the review found, what the maintainer ruled and how the fixes landed is
in the landing note
[Third whole-branch review: lenses, findings and fixes (01dee4b4..94b8fd40, 2026-10-01 to 2026-10-02)](../release-candidate-review.md#third-whole-branch-review-lenses-findings-and-fixes-01dee4b494b8fd40-2026-10-01-to-2026-10-02).

First wave:

- BRIEF.md - the brief every lens reviewer shared: ground rules, severity
  scale, output format.
- engine.md - the engine: monotone leaves, successive-conditional checks,
  threading, mutation invariants, edge cases.
- engine-verify.md - verification of engine.md.
- bridge.md - the R-to-C++ bridge and the flat C API.
- bridge-verify.md - verification of bridge.md.
- rfit.md - the R fitting functions against base R conventions and 0.9-34.
- rgen.md - methods on fits and samplers: predict and the other generics,
  save, load and copy, the deprecation stubs.
- r-verify.md - verification of rfit.md and rgen.md.
- docs.md - help pages and NEWS checked against the code.
- docs-verify.md - verification of docs.md.

Second wave:

- ingest.md - data ingestion: factors, missing values, sparse input, dates,
  the cut grid.
- auxiliary.md - xbart, partial dependence, rbart_vi and diagnostics.
  Written as aux.md, a name Windows cannot check out, and cited under
  that name by the files here.
- wave2a-verify.md - verification of ingest.md and auxiliary.md.
- families.md - every response family end to end.
- families-verify.md - verification of families.md.
- multiforest.md - models with several forests: BCF, variance forests,
  multinomial, interactions and blocks, linear and GP leaves.
- multiforest-verify.md - verification of multiforest.md.
- anchor.md - 0.9-34 against 1.0-0 on every model both can fit, a re-run of
  the second review's [anchor-main.md](../review-2026-08-24/anchor-main.md).
- mutation.md - the mutation leg: 62 one-token defects planted in code
  changed since 7ad0bbea, and which tests caught each.
