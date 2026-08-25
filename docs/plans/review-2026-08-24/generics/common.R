suppressPackageStartupMessages({
  library(dbarts); library(survival)
  if (requireNamespace("posterior", quietly = TRUE)) library(posterior)
})
grDevices::pdf(NULL)

OUT <- Sys.getenv("MXR2_OUT")
TAG <- Sys.getenv("MXR2_TAG")
csvPath <- file.path(OUT, "fill-grid.csv")

escapeCsv <- function(s) {
  s <- as.character(s); s[is.na(s)] <- ""
  s <- gsub('"', '""', s, fixed = TRUE); s <- gsub("[\r\n]+", " ", s)
  paste0('"', s, '"')
}
conn <- file(csvPath, open = if (identical(TAG, "live")) "wt" else "at")
if (identical(TAG, "live"))
  writeLines("tag,fit,condition,generic,type,sample,combine,outcome,detail", conn)
row <- function(fit, condition, generic, type, sample, combine, outcome, detail) {
  writeLines(paste(escapeCsv(TAG), escapeCsv(fit), escapeCsv(condition), escapeCsv(generic),
                   escapeCsv(type), escapeCsv(sample), escapeCsv(combine),
                   escapeCsv(outcome), escapeCsv(detail), sep = ","), conn)
  flush(conn)
}

rawErrorPatterns <- c("^'arg' should be one of", "unused argument", "subscript out of bounds",
  "could not find function", "is missing, with no default", "missing value where TRUE/FALSE needed",
  "non-numeric argument", "invalid subscript", "object '.*' not found",
  "attempt to apply non-function", "[$] operator is invalid for atomic vectors", "invalid 'type'",
  "NA/NaN/Inf", "unable to find an inherited method", "no applicable method",
  "incorrect number of dimensions", "non-conformable", "cannot coerce", "undefined columns selected",
  "matched by multiple actual arguments", "argument .* is missing",
  "arguments imply differing number of rows", "invalid first argument",
  "argument is of length zero", "is not a function", "must be a", "invalid graphics")
classifyError <- function(msg) {
  for (p in rawErrorPatterns) if (grepl(p, msg, perl = TRUE)) return("error-without-reason")
  "refused"
}
describeValue <- function(v) {
  if (is.null(v)) return("NULL")
  cls <- paste(class(v), collapse = "/")
  d <- tryCatch(dim(v), error = function(e) NULL)
  shape <- if (!is.null(d)) paste(d, collapse = "x") else paste0("len", length(v))
  dn <- tryCatch(dimnames(v), error = function(e) NULL)
  dnStr <- if (!is.null(dn)) paste0("dn=[", paste(vapply(dn, function(z)
      if (is.null(z)) "NULL" else paste0("{", paste(head(z,3), collapse=","),
        if (length(z)>3) ",..." else "", "}"), character(1)), collapse=","), "]") else ""
  extra <- setdiff(names(attributes(v)), c("dim","dimnames","class","names","row.names"))
  aStr <- if (length(extra)) paste0("attrs=[", paste(extra, collapse=","), "]") else ""
  trimws(paste(cls, shape, dnStr, aStr))
}
evalCell <- function(fit, condition, generic, type, sample, combine, expr) {
  out <- withCallingHandlers(
    tryCatch({ v <- force(expr); list(ok=TRUE, val=v, msg=NA_character_) },
             error = function(e) list(ok=FALSE, val=NULL, msg=conditionMessage(e))),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
  if (out$ok) row(fit, condition, generic, type, sample, combine, "accepted", describeValue(out$val))
  else row(fit, condition, generic, type, sample, combine, classifyError(out$msg), out$msg)
  invisible(out)
}
