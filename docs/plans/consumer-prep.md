# consumer-prep: the consumer packages stop depending on what the control migration moves

Status: PLANNED.

agent: sonnet implementer, one per package; opus reviewer, one per package.
rng: NEUTRAL. R only, in four other repositories; no dbarts file changes. Every seeded fit is bit for bit
what it is now.
window: any time; before the dbarts change that removes the control's cut-count slot, which breaks stan4bart
as it stands.
budget: ~130 lines across the four packages (stan4bart ~55, treatSens ~45, bartCause ~15, bairrtt ~15),
upper figure 200. Plans have run 1.5-2x low.

## Goal

stan4bart, bartCause, treatSens and bairrtt state a tree count on the forest it belongs to, a cut count on
the data object, and read a family from the model, in spellings that build the same sampler on today's
dbarts and after the control migration (dec-B239, dec-B240, dec-B241). The later dbarts changes then need
no consumer commit on the same day for anything listed here.

## Context

- The migration takes the tree count, the cut count, the quantile rule and the binary flag off the control.
  Along the way the control loses its `n.cuts` slot outright (dec-B241: no mirrored copy), then
  `dbartsControl` stops taking `n.trees`, `n.cuts` and `useQuantiles` as formals and warns when given them
  (dec-A02, dec-B250). The slots `n.trees`, `useQuantiles` and `binary` stay readable as mirrors.
- Three spellings exist today and are the end state's own: a tree count as the first forest's,
  `forests = list(forest(n.trees = ))` (dec-B241; man/forest.Rd: the first forest's count governs); a cut
  count written to the data object's `n.cuts`, which [`dbartsSpec`](../../R/spec.R) keeps when it holds one
  count per column; the family read from `model@family`.
- Three spellings do not exist yet, so nothing here uses them: `dbarts(n.trees = )` (refused: unused
  argument), a bare `forest()` outside a list (refused), and any home for the quantile rule other than the
  control (`dbartsData` has no such slot or argument).
- A single forest declared through `forests =` gets a `forest` column on `getTrees()` today, which dec-B249
  removes. The spelling is therefore used only where no tree table is read: treatSens, bairrtt, bartCause's
  `bcf` (two forests, the column is there either way) and stan4bart's `mvbart`. stan4bart's main fit hands
  out tree tables and keeps its count where it is (below).
- Run on the tip's build, old spelling against new, per candidate: the engine's tree count per forest, the
  data object's cut counts, the split points held, the quantile rule, the binary flag, the family and the
  draws of a seeded short run are identical for stan4bart's cut count (gaussian, probit), `bcf` (gaussian,
  probit, logistic), bairrtt's two surfaces after `setCutPoints` and `sampleTreesFromPrior`, and `mvbart`'s
  per-equation sampler.
- Run on each package itself, tip against tip plus the edit, each installed privately against the tip's
  dbarts: stan4bart, six seeded fits (default, 20 cuts, 100 as a double, no trees kept, quantile cuts,
  binary) identical; `mvbart`, six seeded fits identical; treatSens's `cibart` under one seed with three
  treatment models identical; bairrtt, seeded two-chain fits at 75, 7 and 8 trees identical; bartCause's
  `bcf` suite, whose hand-written dbarts drivers still name the count on the control, passes with the same
  counts (145 expectations in its bcf file).
- treatSens's hand-built triple and the one `dbartsSpec` returns have the same control, model and data
  slots for the outcome model. For the propensity model one slot differs and no draw does: the hand-built
  model stores a chi-squared residual prior on a probit model, which has no residual scale, and `dbartsSpec`
  stores the fixed one.

## What each package needs, and what is left out

| site | today | change here | why it is forward-compatible |
|---|---|---|---|
| stan4bart, `stan4bart_fit`: `bart_args$n.cuts` | routed to `dbartsControl` by its formals, then read back from `control@n.cuts` onto the data object | validated in stan4bart and written to the data object; never given to the control; named in the list of known `bart_args` | the slot read stops working when the slot goes; the data object is the count's home before and after; the name stays accepted when it leaves the control's formals |
| stan4bart, `mvbart`: `n.trees` | on the shared control | on each equation's `dbarts` call as its one forest's count | the forest is the count's home; `mvbart` reads no tree table |
| bartCause, `bcf`: prognostic count; binary; the reported counts | `dbartsControl(n.trees = )`; `sampler$control@binary`; `control@n.trees` | the first `forest(n.trees = )`; the model's family; `bcf`'s own argument | all three are the end state's spellings |
| treatSens, `makeBartSpecs` | builds the triple by hand: count and cut count on the control, `control@binary` written, two unexported dbarts functions called | `dbartsSpec(data, control, leaf.prior, forests = list(forest(n.trees = )), family = )`; the cut count on the data object as now | one exported entry; nothing stated on the control; no slot written |
| bairrtt, `irt_causal_bart`: `n_trees` | `dbartsControl(n.trees = )` shared by two `dbarts` calls | `forests = list(dbartsForests$forest(n.trees = ))` on both | as `mvbart` |

Left out, each because it would need a second edit later:
- The quantile rule in stan4bart's and `mvbart`'s `bart_args`. Its only home today is the control. Both keep
  working unedited until `dbartsControl` drops the name; the edit lands once the data object can hold it.
- `mvbart`'s `bart_args$n.cuts`, which belongs on `dbarts(n.cuts = )`, an argument that does not exist yet.
  Same timing.
- Anything that picks a spelling by inspecting dbarts's formals at run time: one of its two paths is dead on
  release day.

Needs no edit at all, checked by reading the routing against each stage of the migration:
- stan4bart's `bart_args$n.trees`. It goes to the control while that is the only formal of that name, to
  both the control and `dbartsSpec` (same value; the flat argument wins, dec-B116) once `dbartsSpec` has
  one, and to `dbartsSpec` alone after the control drops it. The two existing loops do this.
- `bcf`'s extra arguments: a name the control takes goes to the control, a name `dbarts` takes goes to
  `dbarts`; a cut count moves from the first loop to the second by itself.
- bartCause's `bart` argument filter and its `bart_args` forwarding: `bart` keeps all three names.
- stan4bart's C read of `keepTrees`; treatSens's benchmarking control, which names none of the three.

## Constraints

- No dbarts change. Each package is installed into its own private library with the dbarts tip's library
  ahead of it on `R_LIBS`; nothing goes into the user library.
- Draws do not move: each package's seeded results are identical before and after, compared by the
  reviewer on two private installs, not taken from the implementer.
- No new user-visible behaviour but two: stan4bart's messages for a bad `n.cuts` are its own (same values
  refused: zero, negative, fractional, NA, non-numeric); `mvbart` refuses `bart_args$n.trees` by name where
  today it is ignored in favour of `mvbart`'s own `n.trees` (run: 99 there beside 7 gives the 7-tree fit).
- Each repository's own conventions: no attribution, no markdown in commit messages. bairrtt is committed
  and pushed on main; the others on their compat branches, pushed after review.
- Out of scope: the per-forest prior arguments of dec-B246, which rework `bcf`'s `forest()` calls in their
  own arc; respelling tests inside dbarts.

## Steps

One implementer per package, each step verifiable alone. The four are independent.

1. bairrtt (branch main). In `irt_causal_bart`'s chain setup, drop `n.trees` from the `dbartsControl`
   call and pass `forests = list(dbarts::dbartsForests$forest(n.trees = ))` to the response and the
   assignment `dbarts` calls. Test, in test-causal.R: seeded fits at 7 and 8 trees differ. Gates:
   `R CMD INSTALL .`, `tinytest::test_package("bairrtt")`, `air format --check .`,
   `lintr::lint_package()`.
2. stan4bart (branch bartcore), two commits.
   a. `stan4bart_fit`: `n.cuts` joins the explicit known names and is excluded from the control loop; when
      given, a small resolver returns one whole positive count per BART column (length one or one per
      column) and its result is written to the data object; when not, the data object is left as built and
      dbarts supplies its default. Test, in test-09-bartArgs.R: the check that
      reads `control@n.cuts` reads the stored data object's `n.cuts`, at 20, not at the default it cannot
      tell from; each bad value is refused naming `n.cuts`; a 20-cut fit differs from the default under one
      seed.
   b. `mvbart`: `n.trees` leaves the shared control and goes on each equation's `dbarts` call as
      `forests = list(dbarts::dbartsForests$forest(n.trees = ))`; `bart_args$n.trees` is refused beside the
      seed. Test, in test-16-mvbart.R: fits at 7 and 8 trees differ; the refusal; the file's independent
      BART comparator is respelled the same way.
   Gates: `R CMD INSTALL --preclean .`; `NOT_CRAN=true tinytest::test_package("stan4bart", at_home = TRUE)`;
   `Rscript benchmarks/R/compare-posterior.R compare-all` on the baseline its MANIFEST marks current.
3. treatSens (branch dbarts-1.0, the linked worktree). `makeBartSpecs` keeps its data object and the write
   of 100 cuts per column, builds a control that names neither count, and returns the control, model and
   data of one `dbarts::dbartsSpec` call given the leaf prior expression, the family ("probit" or
   "gaussian") and the forest's count. The sigma estimate, the `control@binary` write, the `parsePriors`
   call and the hand-made model go. Test, a new testthat file: for an outcome and a propensity triple, the
   sampler built from it holds the asked number of trees (distinct tree ids in `getTrees()`), 100 cuts per
   column, the family asked for, and a finite data sigma for the outcome. Gates:
   `R CMD INSTALL --preclean .`; the testthat suite with `NOT_CRAN=true`; `R CMD check --as-cran` on a
   built tarball.
4. bartCause (branch dbarts-1.0), two commits, after the implementer now working there has committed
   (the working tree holds uncommitted edits to other files; start from a clean HEAD).
   a. `bcf`'s sampler construction: the coerced `n.trees` goes on the prognostic `forest()` call and leaves
      the control's arguments; the binary flag is the family being probit or logistic; the result's
      `n.trees` reads the argument. The suite is run now, before 4b: its hand-written drivers still name
      the count on the control, so their bitwise comparisons prove the two spellings equal.
   b. The two hand-written drivers (test-03-responseFit.R, test-14-bcf.R) state the count on the first
      forest. The file's binary test (no `sigma`, the link applied last) pins the family read for probit;
      a logistic fit is added to it.
   Gates: private install; the testthat suite with `NOT_CRAN=true`, expectation and warning counts as
   before; `R CMD check --as-cran` on a built tarball.

## Verification

- Per package, by the reviewer: two private installs, the commit before and the commit after, against the
  dbarts tip; the seeded fits listed in Context identical; the package's gates above.
- `git grep` over each package's R, src and tests for `control@n.cuts`, `control@binary`, `control@n.trees`
  and for `n.trees` or `n.cuts` inside a `dbartsControl(` call: the only hits left are the two routing
  loops that take `bart_args` or extra arguments by the control's formals, and tests that pass a count
  through `bart_args`.
- Mutation, one per package: drop the new statement (the forest's count; the write to the data object) so
  the default applies: the package's test fails. The tests pin that the value reaches the sampler; the
  search above is what checks the spelling, both spellings building the same sampler today.
- No dbarts gate applies: no dbarts file changes.

## Calls made in planning

- `mvbart` refuses `bart_args$n.trees` by name, where the value was silently ignored: a stated value is
  used or refused, never dropped.
- Found on the way, not planned here: treatSens's hand-built probit model stores a residual prior the
  family cannot use (it disappears with step 3); stan4bart's test of `n.cuts` asserts the default value, so
  it passes whether or not the argument is honoured (fixed in step 2a).
