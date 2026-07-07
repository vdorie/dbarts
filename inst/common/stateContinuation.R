# Structural round-trip predicate for the state tests after bitwise
# continuation was dropped: a restored sampler must reconstruct the model -
# tree structure and leaf parameters - exactly, and sigma to within the
# original-scale round trip. Returns a logical so callers wrap it in a
# top-level expect_true (tinytest does not collect expect_ calls nested in a
# helper). A function leaf's flattened value is a reporting mean whose sum
# order the canonical rebuild does not preserve, so tree.values is not
# compared; the constant-leaf value round trip is gated by test-capi.R's
# predictions-identical check.
statesAgree <- function(reState, saved) {
  agree <- length(reState) == length(saved)
  for (ci in seq_along(saved)) {
    reForests <- reState[[ci]]$forests
    savedForests <- saved[[ci]]$forests
    agree <- agree && length(reForests) == length(savedForests)
    for (fi in seq_along(savedForests)) {
      agree <- agree &&
        identical(reForests[[fi]]$tree.vars, savedForests[[fi]]$tree.vars) &&
        identical(reForests[[fi]]$tree.sizes, savedForests[[fi]]$tree.sizes) &&
        identical(reForests[[fi]]$tree.flags, savedForests[[fi]]$tree.flags) &&
        identical(reForests[[fi]]$tree.params, savedForests[[fi]]$tree.params)
    }
    agree <- agree &&
      isTRUE(all.equal(
        reState[[ci]]$sigma,
        saved[[ci]]$sigma,
        tolerance = 1e-8
      ))
  }
  agree
}
