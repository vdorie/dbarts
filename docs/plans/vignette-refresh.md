# vignette-refresh

agent: sonnet
rng: neutral (docs; vignette code must run)
budget: ~300 lines of Rmd

## Goal

The vignettes describe the 1.0-0 surface: family selection, DART,
sparse inputs, missingness, and the new leaf models appear where users
would look for them, and the saved-trees vignette handles non-constant
leaves.

## Context

- Neither vignette mentions family, dart, sparse, missing, or the
  linear/gp leaf models (TODO note).
- vignettes/working_with_saved_trees.Rmd's traversal code assumes
  constant leaves; linear/gp change the getTrees columns.
- If flat-format-v2's reporting decision lands (value = NA +
  directions for all categorical rules), the saved-trees examples must
  match - sequence after that decision or write against the current
  output and note the dependency.

## Constraints

- Vignette chunks execute at build time: keep added fits tiny (small n,
  few iterations, fixed seeds) so R CMD build stays fast.
- Tone and depth match the existing vignettes; new sections, not a
  rewrite.
- Out of scope: pkgdown site structure; man pages (prior-constants
  covers their gaps).

## Steps

1. Main vignette: a short section each for family (probit vs logistic,
   the latent scale), dart (when and why, varprobs), data.frame
   ingestion with factors/sparse/missing (the categorical-vs-indicators
   story and the missing = "incorporate" default).
2. Saved-trees vignette: branch the traversal example on the leaf
   model; show the directions column for a categorical rule.
3. Cross-reference docs/design/ for readers who want the internals.

## Verification

- R CMD build (vignettes execute); R CMD check clean.
- Build time increase measured and kept small.

## Landing note (2026-07-07, see merge commit)

Landed: the Gibbs vignette gains Response Families, DART, and data
frame/factor/missing sections; the saved-trees vignette documents the
full getTrees column set (directions, missing, beta.*), adds a
categorical-splits section, and rebuildTree handles linear leaves -
every example verified against installed behavior before the prose
was written. Build times grew < 0.1s per vignette. Gates: install,
clean vignette builds, tinytest 2470 ok. 115 insertions.
