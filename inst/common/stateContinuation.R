# Structural round-trip predicate for the state tests after bitwise
# continuation was dropped: a restored sampler must reconstruct the model -
# tree structure and leaf parameters - exactly, and sigma to within the
# original-scale round trip. A function leaf's flattened value is a
# reporting mean whose sum order the canonical rebuild does not preserve, so
# tree.values is not compared; the constant-leaf value round trip is gated
# by test-capi.R's predictions-identical check.
#
# Asserts each compared field directly, so a mismatch attributes to that
# field at the caller's line; callers no longer wrap the call in
# expect_true. Every source() of this file must pass local = TRUE for those
# assertions to reach the run's masked expect_* bindings. Pass
# expect = FALSE for a silent boolean predicate, e.g. to assert two states
# DIFFER.
statesAgree <- function(reState, saved, expect = TRUE) {
  agree <- TRUE
  isField <- function(current, target, label) {
    agree <<- agree && identical(current, target)
    # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
    if (expect) expect_identical(current, target, info = label)
  }
  if (length(saved) == 0L) {
    if (expect) {
      # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
      expect_true(
        FALSE,
        info = "statesAgree: zero-length state, nothing was compared"
      )
    }
    return(invisible(FALSE))
  }
  isField(length(reState), length(saved), "number of chains")
  # the leaf prior's two halves (k, leaf.scale): a restore that reconstructed
  # the trees but not the calibration they were drawn under would leave a
  # hybrid, and a comparison of the trees alone cannot see it
  forestFields <- c(
    "tree.vars",
    "tree.sizes",
    "tree.flags",
    "tree.params",
    "k",
    "leaf.scale"
  )
  # chain-level blocks the per-forest loop cannot reach: the variance forest
  # sits outside the forests list, so its LIVE and SAVED trees are siblings of
  # it. The saved buffer is the only record of the kept samples' scale surface,
  # and a restore that rebuilt the live trees alone would leave predict
  # replaying the destination's identity fill unnoticed. variance.values IS
  # compared: a scale leaf's value is the surface itself, not a reporting mean.
  # A homoscedastic pair reads NULL on both sides and compares vacuously, so
  # the loop only bites for a heteroscedastic caller.
  chainFields <- c(
    "variance.vars",
    "variance.values",
    "variance.sizes",
    "variance.flags",
    "variance.masks",
    "variance.saved.vars",
    "variance.saved.values",
    "variance.saved.sizes",
    "variance.saved.flags",
    "variance.saved.masks"
  )
  for (ci in seq_along(saved)) {
    reForests <- reState[[ci]]$forests
    savedForests <- saved[[ci]]$forests
    isField(
      length(reForests),
      length(savedForests),
      sprintf("chain %d forest count", ci)
    )
    for (fi in seq_along(savedForests)) {
      label <- sprintf("chain %d forest %d", ci, fi)
      for (field in forestFields) {
        isField(
          reForests[[fi]][[field]],
          savedForests[[fi]][[field]],
          paste(label, field)
        )
      }
    }
    for (field in chainFields) {
      isField(
        reState[[ci]][[field]],
        saved[[ci]][[field]],
        sprintf("chain %d %s", ci, field)
      )
    }
    sigmaLabel <- sprintf("chain %d sigma", ci)
    agree <- agree &&
      isTRUE(all.equal(
        reState[[ci]]$sigma,
        saved[[ci]]$sigma,
        tolerance = 1e-8
      ))
    if (expect) {
      # nolint next: object_usage_linter. tinytest attaches expect_* at run time.
      expect_equal(
        reState[[ci]]$sigma,
        saved[[ci]]$sigma,
        tolerance = 1e-8,
        info = sigmaLabel
      )
    }
  }
  invisible(agree)
}
