# Splits a getTrees() data frame into one group per (forest, chain, sample,
# tree) tuple present in its columns.
splitTrees <- function(trees) {
  cols <- intersect(c("forest", "chain", "sample", "tree"), names(trees))
  split(trees, trees[cols], drop = TRUE)
}
