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
    agree <- agree &&
      identical(reState[[ci]]$tree.vars, saved[[ci]]$tree.vars) &&
      identical(reState[[ci]]$tree.sizes, saved[[ci]]$tree.sizes) &&
      identical(reState[[ci]]$tree.flags, saved[[ci]]$tree.flags) &&
      identical(reState[[ci]]$tree.params, saved[[ci]]$tree.params) &&
      isTRUE(all.equal(
        reState[[ci]]$sigma,
        saved[[ci]]$sigma,
        tolerance = 1e-8
      ))
  }
  agree
}
