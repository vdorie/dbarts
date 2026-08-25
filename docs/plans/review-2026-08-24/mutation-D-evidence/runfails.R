suppressPackageStartupMessages(library(dbarts))
args <- commandArgs(trailingOnly = TRUE)
listFile <- args[1]; outFile <- args[2]
files <- readLines(listFile)
root <- Sys.getenv("SRCROOT")
res <- character(0)
for (f in files) {
  out <- tryCatch({
    r <- suppressWarnings(tinytest::run_test_file(file.path(root, f), verbose = 0))
    lg <- as.logical(r)
    ok <- sum(lg, na.rm=TRUE); bad <- sum(!lg | is.na(lg))
    det <- ""
    if (bad > 0) {
      idx <- which(!lg | is.na(lg))
      det <- paste(vapply(idx, function(i) {
        a <- r[[i]]
        paste0(attr(a,"file")," :", attr(a,"fst"), "-", attr(a,"lst"))
      }, ""), collapse=";")
    }
    sprintf("%s\tPASS=%d\tFAIL=%d\t%s", basename(f), ok, bad, det)
  }, error = function(e) sprintf("%s\tPASS=NA\tFAIL=ERR\t%s", basename(f), gsub("[\t\n]"," ",conditionMessage(e))))
  res <- c(res, out)
  cat(out, "\n", sep="")
}
writeLines(res, outFile)
