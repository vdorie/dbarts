# Splits a getTrees() data frame into one group per (chain, sample, tree)
# triple present in its columns.
splitTrees <- function(trees) {
  cols <- intersect(c("chain", "sample", "tree"), names(trees))
  split(trees, trees[cols], drop = TRUE)
}
