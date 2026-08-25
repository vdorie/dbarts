source(file.path(Sys.getenv("MXR2_FILL"), "common.R"))
source(file.path(Sys.getenv("MXR2_FILL"), "fixture.R"))
source(file.path(Sys.getenv("MXR2_FILL"), "grid.R"))
rdsDir <- file.path(OUT, "rds")
for (nm in names(builders)) {
  p <- file.path(rdsDir, paste0(nm, ".rds")); if (!file.exists(p)) next
  obj <- tryCatch(readRDS(p), error = function(e) e)
  if (inherits(obj, "error")) { row(nm, "reload", NA,NA,NA,NA, "RELOAD-FAILED", conditionMessage(obj)); next }
  row(nm, "reload", NA,NA,NA,NA, "reloaded", paste0("class=", paste(class(obj), collapse="/")))
  fieldCensus(nm, obj)
  runGenericGrid(nm, obj, extraCombine = FALSE)
}
close(conn); cat("rds phase done\n")
