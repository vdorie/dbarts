## A forests list whose forest k multiplies bases[[k]], handed over as a value,
## a forest with no basis where that entry is NULL; `forests` gives each
## forest's other settings, a forest() stating none where it is shorter.
basisForests <- function(bases, forests = NULL) {
  forest <- dbarts::dbartsForests$forest
  lapply(seq_len(max(length(bases), length(forests))), function(k) {
    spec <- if (k <= length(forests)) forests[[k]] else forest()
    spec["basis"] <- list(if (k <= length(bases)) bases[[k]])
    spec
  })
}
