# Plan process

The TODO at the repo root is an unordered backlog; most items name a
plan file here, a few record instead in a design doc or a
differently-named plan. This document covers the plan-file format, the
citation rule, the RNG-class gates a change needs, how a landing is
recorded and carried out, the reviewer's checklist, gate hygiene, and
what CI runs on a push.

## Plan format

One file per TODO item, `<item>.md`. A small decision-gated stub
targets roughly 80 lines; a landed or multi-slice item's file commonly
runs to several hundred lines (occasionally thousands) once Landing
notes, amendments, and session logs accumulate - there is no enforced
cap on the final file. Front block, should carry:
`agent:`, `rng:`, `window:` (if any), `budget:` (expected diff size).
Sections, should carry:

- Goal: two or three sentences; what is true after the item lands.
- Context: durable references (a code symbol name, a doc section title
  or item label) and design-doc pointers; no narrative a pointer can
  replace. Never a bare line number - see Cross-references.
- Decision (decision-gated items only): the question, a recommendation,
  and what evidence would change it. VD signs off before implementation.
- Constraints: gates, contract freezes, explicit out-of-scope list.
- Steps: numbered; each independently verifiable.
- Verification: exact commands and expected outcomes.

## Cross-references

Documentation cites code by SYMBOL. Every citation - of code or of
another document - is written in one machine-checked grammar: an
ordinary markdown link, so it renders and resolves on GitHub. The target
is the file, written relative to the citing document's own directory;
the text is what in that file is meant:

    [`rollTreeResidual`](../../src/bartcore/chain.hpp)  current state, by symbol
    [`GaussianResponse`](../../src/bartcore/model.hpp), [`ProbitResponse`](../../src/bartcore/model.hpp)  several symbols, one cite
    [`bart2()`](../../R/bart.R)  a trailing () is stripped
    [`dbartsSampler$setResponse`](../../man/dbartsSampler-class.Rd)  split on $ as well as ::
    ["hurdle.lognormal"](../../inst/tinytest/test-argument-surface.R)  verbatim fragment
    [Implementation record](../design/data-ownership.md#implementation-record)  doc to doc, by heading slug
    [src/bartcore/chain.hpp:2581](https://github.com/vdorie/dbarts/blob/d477a46be658d885b769dc186c2832991348187d/src/bartcore/chain.hpp#L2581)  history, a line at a commit
    retired: [`refuseHostMutation`](../../R/spec.R)  the construct is gone
    unresolved: [inst/include/dbarts/R_C_interface.hpp:40](https://github.com/vdorie/dbarts/blob/d477a46be658d885b769dc186c2832991348187d/inst/include/dbarts/R_C_interface.hpp#L40)  cannot be placed

- A target is a path relative to the citing document, and it must
  resolve to a file that exists: a citation into another repository is
  prose, not a cite. A link whose text is the target's own path or
  basename is a plain file link, not a cite; any OTHER link into a
  tracked file fails, because it reads as a cite and is checked by
  nothing.
- A symbol is checked component-wise - split on `::` and `$`, each part a
  whole-word token somewhere in that file. Adjacency is not required and
  no line number is involved, so a definition that moves needs no re-pin.
  Several symbols of one file are several links separated by `, `; that
  run is one cite, and a marker in front of it covers all of it.
- A `.md` target carrying a `#` fragment is a heading cite: the fragment
  must be GitHub's slug of one of the target document's headings -
  lowercased, everything but letters, digits, spaces, hyphens and
  underscores dropped, spaces to hyphens - and the link text is that
  heading, so a renamed section breaks the cite instead of drifting
  quietly. A heading of the citing document itself is `(#slug)`, with no
  path.
- A quoted fragment is checked as a literal substring. That is the form
  for a file with no symbols to name - a tinytest script; a tests/cpp
  citation names its test function instead. A bracket inside the fragment
  is backslash-escaped so the link still parses.
- A line number appears ONLY in the history form, a link to
  `https://github.com/vdorie/dbarts/blob/<sha>/<path>#L<a>-L<b>` whose
  sha is the full 40 hex digits and must be an ancestor of HEAD, and
  whose text repeats `path:a-b`, the claim the URL makes. The claim is
  read AT THAT COMMIT: the path must exist in that commit's tree and the
  cited line must be within the file's length there. A file deleted since
  therefore still cites cleanly, while a line that never existed at the
  named commit does not - which is what catches a bare `:NNN` whose file
  was guessed wrong. Landing notes and a plan's Landing section cite this
  way, pinned to the commit the note describes.
- `retired:` before a cite says the named CONSTRUCT is gone, so its
  content is not checked; the location must still be real - the target
  resolves, and a history cite's sha, path and line are checked as usual.
  The prose around it must say the thing is gone.
- `unresolved:` before a HISTORY cite says the LOCATION itself cannot be
  established: the target of a frozen record's line reference could not
  be placed at any candidate commit. The link must still parse and the
  sha must still be an ancestor; the path and the line are not checked.
  It is the honest marker for a record that cannot be repaired, never a
  way to quiet a cite that can be.
- A marker binds only when it sits immediately before the link, so
  `un-retired:` and any other word ending in `retired:` disarms nothing.
- A cite parked inside a fenced code block is not checked at all: a fence
  is a code sample or a transcript. Do not put a live cite in one.
- A symbol is found by token search over the whole file, so a name
  occurring only in a comment or a string satisfies its cite. Qualify a
  name the target file defines more than once - `ProbitResponse::computeLogLikelihood`,
  not `computeLogLikelihood` - so the token that answers is the one the
  sentence means.

Backticks belong inside a symbol cite's link text and nowhere else in
the grammar. `tools/check-doc-freshness.R`
checks every cite in everything under docs/, the top-level README.md,
man/*.Rd and vignettes/*.Rmd, and fails on any bare line reference left
outside the history form, so an unconverted citation cannot pass
unnoticed. State each fact in one home doc; elsewhere link to it by
title, do not restate it (a copied fact is a second thing to keep in
sync, and the one that rots).

## RNG classes and their gates

- neutral: draws unchanged.
  Gates: tests/cpp component tests; full tinytest suite.
- shifting: draws change, the posterior does not.
  Gates: the above, plus regenerate RNG-locked snapshots by replaying
  whole test files, re-record the equivalence baseline, and pass the
  statistical (z) mode against the previous baseline.
- posterior-changing: the stationary distribution or a default changes.
  Gates: all of the above, plus the exact-posterior gates
  (.github/workflows/exact-gates.yaml's per-family list) and a design
  note in docs/design/.

Hot-path changes of any class additionally need bench-sampler.R compare
on a quiet machine (maintainer-run; never concurrent with other load).

## Landing

The landing record is two edits: append the plan's `## Landing` (or
`## Landing note`) note, and bump the matching `docs/design/<x>.md`
`Status:` line to `LANDED <date> (<commit>)`. The design Status line is
the record most often missed; check it explicitly at every landing.
release-candidate-review.md inserts its landing notes newest-first
under its "## Landing notes" heading; every other plan file appends at
the end. Never insert into the middle of a note sequence another doc
cites.

The procedure around those edits, as practiced:

1. Implement in a linked worktree branched off `origin/bartcore`, with
   its own private R library (`R CMD INSTALL -l <lib> .` and `R_LIBS=<lib>`
   on every R call; `~/.Renviron` overrides `R_LIBS_USER`, so the prefix
   is not optional). One writer per worktree.
2. Diff review by a second reader, then the gate battery for the change's
   RNG class (above) run independently of the implementer, against the
   slice's own library. `--preclean` on every engine commit: a stale
   object silently fails the bitwise gates.
3. Push the reviewed sha, fast-forward the main checkout, then a separate
   records commit carrying the real landed hash. Cherry-picking into an
   integration worktree renames every commit, so records cite only hashes
   that are ancestors of `origin/bartcore` (`git merge-base --is-ancestor`).
4. Clean up the worktree and its library only after the push succeeded,
   chaining with `&&`, and never through a pipe: `git merge --ff-only X |
   tail -1` reports tail's exit status, not the merge's.
5. Watch CI to green (below) before the next slice branches off the tip.

Independent file-disjoint slices may run in parallel worktrees off one
base, each with its own implementer and gate run; they stack by rebase in
slice order and the batch is gated once with one merged-tree battery plus
cross-slice probes. Hard chains stay serial.

## Reviewer checklist

Cheap gates are re-run by the reviewer, not trusted from the report.

- `air format --check .` tree-wide, and `lintr::lint()` on every touched R
  file against the slice's own library (a stale library manufactures
  `object_usage_linter` false positives). `lintr::lint_package()` whenever
  a slice moves names, touches many R files or edits NAMESPACE; per-file
  lintr is necessary, not sufficient. `air.toml` excludes docs/, whose
  verbatim records must not be reformatted.
- A `_pkgdown.yml` entry plus `pkgdown::check_pkgdown(".")` for every new
  exported Rd topic.
- `inst/NEWS.Rd` must parse when touched: gate on a non-NULL result from
  `tools:::.build_news_db_from_package_NEWS_Rd("inst/NEWS.Rd")` and its
  entry count. `tools:::.build_news_db(".")` is a silent no-op.
- `R CMD check --as-cran` from a tarball built from a clean copy staged
  outside the tree, for any commit touching R/ or man/. CRAN's core limit
  fires parallel's own error before the package's validation does.
- `tools/check-rc-codoc.R` (reference-class methods, which `codoc` cannot
  see) and `tools/check-doc-freshness.R`, each gated on its OWN exit
  status. A `Rscript ... | tail -1` chain masks a failure behind tail's
  exit 0.
- Build provenance before trusting a suite: the installed package's mtime
  must postdate the source it claims to test.

## Gate hygiene

- The equivalence harnesses' terminal summary line is not evidence.
  Count the per-scenario "identical draws (same RNG stream)" lines and
  require the full scenario count with no "max |z|" line. Any skipped
  scenario or any |z| means a wrong baseline file or a real change.
  benchmarks/baselines/MANIFEST names the current baselines.
- Count warnings with `withCallingHandlers`; a bare `expect_warning`
  is blind to extras. Muffle a third-party warning at its source only when
  it is meaningless to the user, and pin new leaks by pattern per class
  rather than retreating to a pattern-only expectation.
- Prove a gate discriminates by running the mutation it exists to catch.
  After reverting a mutation, `touch` the file: an `mv` back preserves the
  mtime and the next install keeps the mutated object.
- `expect_equal` on an always-NULL field is a no-op, and a fixture whose
  factor has one level vacates its pin. A constant "filler" predictor
  column is not inert either; never use one to disable splitting.
- A snapshot-valued test is replayed as a whole file when its draws move;
  the values depend on the file's full execution history.
- Sanitizers, locally, before pushing any engine commit that makes new
  numerics reachable: build tests/cpp with
  `OPT="-O2 -g -fsanitize=address,undefined"` and run with
  `ASAN_OPTIONS=detect_container_overflow=0` (container-overflow is a
  false positive from mixing instrumented objects with the
  uninstrumented static libraries). On macOS symbolization spawns `atos`,
  which raises the Developer Tool Access prompt; a prompt means a
  diagnostic fired, so read the count.
- That bullet covers tests/cpp; the R-loaded path (the bridge, the flat C
  entry file, anything only a `.Call` reaches) needs its own run, and it
  IS reachable on macOS. SIP strips `DYLD_INSERT_LIBRARIES` from the
  `bin/R` shell wrapper, not from `$(R RHOME)/bin/exec/R`. Build with
  `R_MAKEVARS_USER=<file adding -fsanitize=address to CFLAGS, CXXFLAGS
  and LDFLAGS> R CMD INSTALL --preclean --no-test-load -l <lib> .` (the
  load test runs under the wrapper and would fail on uninstalled
  interceptors), then run the suite as
  `R_HOME=$(R RHOME) R_LIBS=<lib>
  DYLD_INSERT_LIBRARIES=<.../libclang_rt.asan_osx_dynamic.dylib>
  ASAN_OPTIONS=detect_container_overflow=0 $(R RHOME)/bin/exec/R --vanilla
  --no-echo -f <driver.R>`; `R_HOME` is not optional, the exec binary not
  setting it itself. Skip any test file that spawns a worker R process -
  the child goes through the wrapper and aborts on "Interceptors are not
  working".
- A renamed OPTIONAL state block is a silent misread, not an error, unless
  the format floor moves with it; see the registry rule at
  [`stateFormatVersion`](../../src/R_interface_bartcore.cpp).

## CI

Per-push workflows and what they ignore (read the yaml for the current
lists; this is the shape):

- check-standard, cpp-tests, sanitizers: `docs/**`, `TODO`, `**.md`,
  `benchmarks/**`. cpp-tests and sanitizers also ignore `inst/NEWS.Rd`;
  cpp-tests ignores `man/**`.
- exact-gates: `docs/**`, `TODO`, `**.md`, the MANIFEST, `inst/NEWS.Rd`.
  It does not ignore the rest of benchmarks/ and cancels an in-progress
  run on the same branch, so a records push that touches benchmarks/
  right after a slice push cancels the slice's run; space the pushes or
  rerun.
- lint: `docs/**`, `TODO`, `**.md`, `benchmarks/baselines/**`.
- pkgdown: `docs/**`, `TODO`, `benchmarks/**`, but not `**.md`
  (README is an input).
- doc-freshness: nothing ignored, by design.
- check-api-hash (part of check-standard.yaml): a dedicated job, full
  checkout history, running
  [`tools/check-api-hash.sh`](../../tools/check-api-hash.sh). It reads
  the newest tag matching `v1.*` or `1.*` - `v<major>.<minor>-<patch>`,
  matching R's own version string, or the same without the `v`, this
  repo's existing spelling (`0.8-7`) - picking the newest by version
  across either spelling, and fails when
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) moved since
  that tag but the major/minor pair did not. No `1.*`-or-later tag
  exists yet, so it prints "no release tag, skipped" and passes; this
  repo's pre-1.0-0 tags (`0.8-7`, `bartcore-pre-cran-rebase`) predate
  the major-1 line and do not match either pattern.

So a docs-only or TODO-only push fires doc-freshness alone; run the
freshness guard locally and that is the whole gate. Any `.github/` touch
fires everything. equivalence, rchk, revdep-smoke, sbc and valgrind are
schedule and dispatch only, and GitHub binds those triggers to the
default branch, so they stay dormant until bartcore reaches main.

Reading a red run: "cancelled" usually means a step hit its
`timeout-minutes`; compare the duration to the limit before rerunning. A
red that goes green on a same-commit rerun is the runner-hardware class,
and the same-commit rerun is the first probe. A Linux-only red with macOS
and Windows green is the BLAS/platform class. Either way the fix is
host-independent code, never a retry, and the red run stays red in
history. Poll with `gh run list --commit <full sha>` (a short sha
silently matches nothing) and, after a rerun, by run id. A `run: |` step
wrapping `Rscript -e '...'` is a single-quoted shell string, so no
apostrophe anywhere inside it, comments included.
