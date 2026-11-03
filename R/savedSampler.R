## A sampler saved by 0.9-x carries its own method closures, which R reaches
## through `$` on the object, and a method it did not carry is this version's;
## both fail on 0.9-x's missing slots and symbols with errors that do not say
## why. Every method access goes through here, so the refusal sits where any
## use of such a sampler starts. The five fields 0.9-x had read as before. The
## test is one attribute lookup on the data: 0.9-x's has no 'counts'.
savedSamplerFields <- c("pointer", "control", "model", "data", "state")

## The inherited method reads its name argument from the call as written, so
## it is called with the name as a value and not through callNextMethod.
refObjectDollar <- methods::getMethod("$", "envRefClass")

setMethod("$", "dbartsSampler", function(x, name) {
  if (!(name %in% savedSamplerFields)) {
    data <- get0("data", envir = as.environment(x), inherits = FALSE)
    if (
      inherits(data, "dbartsData") &&
        !("counts" %in% names(attributes(data)))
    ) {
      refuseSavedBy09x(x, "sampler")
    }
  }
  do.call(refObjectDollar, list(x, name))
})
