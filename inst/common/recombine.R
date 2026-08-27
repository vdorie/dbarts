# Recombines a per-forest prediction array back to the single-forest scale
# using the caller's own 'fit' (glue, n.forests, glueForest) and 'shift'.
# Every source() of this file must pass local = TRUE so recombine() resolves
# 'fit' and 'shift' from the caller's own environment at call time.
recombine <- function(perForest, bases, nRows) {
  # nolint next: object_usage_linter. resolves from the caller (file header).
  out <- matrix(shift, nrow(fit$glue), nRows)
  # nolint next: object_usage_linter. resolves from the caller (file header).
  for (k in seq_len(fit$n.forests)) {
    basis <- if (is.null(bases[[k]])) matrix(1, nRows, 1L) else bases[[k]]
    # nolint next: object_usage_linter. resolves from the caller (file header).
    g <- fit$glue[, glueForest == dimnames(perForest)[[3L]][k], drop = FALSE]
    out <- out + (g %*% t(basis)) * perForest[,, k]
  }
  out
}
