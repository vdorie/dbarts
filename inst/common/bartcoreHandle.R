# A public dbartsSampler's own $getTrees always reads forest 1 (the first;
# see docs/plans/test-scaffolding-consolidation.md, decision D1), so a
# Bayesian causal forest's treatment forest, or a multinomial fit's category
# past the first, has no route back into R through the R5 surface. This is
# that route for tests and benchmarks: a thin .Call directly on
# sampler$getPointer(), mirroring $getTrees's own arguments (and none of its
# defaulting, categorical-split decoding or linear-leaf column renaming) for
# an arbitrary forest. forest indexes from 1, as with
# $getForestFits/$getForestAmplitudes/$getLeafPrior - not the 0-based
# convention the bridge's C_dbarts_bartcore_getTrees entry takes.
forestTrees <- function(
  sampler,
  forest,
  treeNums,
  chainNums,
  sampleNums = NULL,
  current = FALSE,
  newdata = NULL
) {
  if (missing(chainNums)) {
    chainNums <- seq_len(sampler$control@n.chains)
  }
  if (!is.null(newdata)) {
    newdata <- as.matrix(newdata)
    storage.mode(newdata) <- "double"
  }
  .Call(
    dbarts:::C_dbarts_bartcore_getTrees,
    sampler$getPointer(),
    as.integer(chainNums),
    if (is.null(sampleNums)) NULL else as.integer(sampleNums),
    as.integer(treeNums),
    as.logical(current),
    newdata,
    dbarts:::rawPredictorMatrix(sampler$data@x),
    dbarts:::resolveForestIndex(forest)
  )
}
