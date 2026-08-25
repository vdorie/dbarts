# ------- documented-field census (man/bart.Rd Value section) -------
bartRdFields <- c("call","yhat.train","yhat.test","yhat.train.mean","yhat.test.mean",
  "sigma","first.sigma","varcount","sigest","resid.dist","y","fit","n.chains","family")
fieldCensus <- function(nm, obj) {
  if (!is.list(obj)) { row(nm, "fields", "names", NA, NA, NA, "n/a", paste0("class=", paste(class(obj), collapse="/"))); return(invisible()) }
  have <- names(obj)
  present <- bartRdFields[bartRdFields %in% have & !vapply(bartRdFields, function(f) is.null(obj[[f]]), logical(1))]
  absent <- setdiff(bartRdFields, present)
  row(nm, "fields", "documented-present", NA, NA, NA, "info", paste(present, collapse=";"))
  row(nm, "fields", "documented-absent", NA, NA, NA, "info", paste(absent, collapse=";"))
  row(nm, "fields", "all-names", NA, NA, NA, "info", paste(have, collapse=";"))
}

genericNames <- c("predict","extract","fitted","residuals","summary","print","plot",
                  "plotTree","survivalProbabilities","as_draws_array","as_draws_df")
classesOfInterest <- c("bart","rbart","bartMultinomial","bartOrdinal","bartNegbin",
                       "bartHurdle","pdbart","pd2bart","dbartsSampler")

methodFor <- function(generic, obj) {
  cl <- class(obj)
  for (c1 in cl) { m <- tryCatch(getS3method(generic, c1, optional=TRUE), error=function(e) NULL)
                   if (!is.null(m)) return(list(m=m, cls=c1)) }
  list(m = NULL, cls = NA_character_)
}
formalVals <- function(fn, nm) {
  if (is.null(fn)) return(NA_character_)
  f <- formals(fn); if (!(nm %in% names(f))) return(NA_character_)
  v <- tryCatch(eval(f[[nm]]), error=function(e) NA_character_)
  if (!is.character(v)) NA_character_ else v
}

runGenericGrid <- function(nm, obj, extraCombine) {
  for (g in genericNames) {
    gf <- tryCatch(get(g, envir=globalenv()), error=function(e) NULL)
    if (is.null(gf)) { row(nm, "generic", g, NA,NA,NA, "generic-absent", "generic not visible"); next }
    mi <- methodFor(g, obj)
    dispatch <- if (is.null(mi$m)) "DEFAULT" else paste0("method:", g, ".", mi$cls)
    types <- formalVals(mi$m, "type"); samples <- formalVals(mi$m, "sample")
    types   <- if (identical(types, NA_character_)) NA_character_ else c(types, "zzz-bad-type")
    samples <- if (identical(samples, NA_character_)) NA_character_ else c(samples, "zzz-bad-sample")
    fn <- names(formals(mi$m))
    combines <- if (!is.null(mi$m) && "combineChains" %in% fn && extraCombine) c(TRUE, FALSE) else NA
    for (ty in types) for (sm in samples) for (cb in combines) {
      a <- list(obj)
      if (!is.null(mi$m)) {
        if ("newdata" %in% fn) a$newdata <- xTest
        if ("group.by" %in% fn) a$group.by <- groupTest
        if (identical(g, "plotTree")) a$treeNum <- 1L
        if (identical(g, "survivalProbabilities")) a$times <- c(0.5,1,1.5)
        if (!identical(ty, NA_character_)) a$type <- ty
        if (!identical(sm, NA_character_)) a$sample <- sm
        if (!identical(cb, NA)) a$combineChains <- cb
      }
      evalCell(nm, dispatch, g, ty, sm, cb,
        expr = { z <- utils::capture.output(v <- do.call(gf, a)); v })
    }
  }
}
