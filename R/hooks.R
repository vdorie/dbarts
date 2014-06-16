.onUnload <- function(libpath)
{
  .Call("cbart_finalize")
  library.dynam.unload("cbart", libpath)
}
