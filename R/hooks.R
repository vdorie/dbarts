.onUnload <- function(libpath)
{
  ## gc is necessary to collect external pointers who have not yet been collected
  ## that have finalizers pointing to the soon-to-unloaded dll
  gc(FALSE)
  .Call("cbart_finalize")
  library.dynam.unload("cbart", libpath)
}
