source(file.path(Sys.getenv("MXR2_FILL"), "common.R"))
source(file.path(Sys.getenv("MXR2_FILL"), "fixture.R"))
source(file.path(Sys.getenv("MXR2_FILL"), "grid.R"))
rdsDir <- file.path(OUT, "rds"); dir.create(rdsDir, showWarnings=FALSE, recursive=TRUE)
for (nm in names(builders)) {
  obj <- tryCatch(builders[[nm]](), error = function(e) e)
  if (inherits(obj, "error")) { row(nm, "build", NA,NA,NA,NA, "BUILD-FAILED", conditionMessage(obj)); next }
  row(nm, "build", NA,NA,NA,NA, "built", paste0("class=", paste(class(obj), collapse="/")))
  fieldCensus(nm, obj)
  runGenericGrid(nm, obj, extraCombine = TRUE)
  saveRDS(obj, file.path(rdsDir, paste0(nm, ".rds")))
}
close(conn); cat("live phase done\n")
