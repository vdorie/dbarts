# A dbartsData object whose missing-value policy is "error": the slot no
# longer has an entry point of its own, missing predictors being modelled
# for whichever rows the fit's na.action kept, but a host driving a sampler
# directly can still ask for new predictors to be refused rather than
# routed. Built here so the tests that exercise that refusal say what they
# are doing.
strictData <- function(...) {
  data <- dbarts::dbartsData(...)
  data@missing <- "error"
  data
}
