# Each forest record's count, base, power and drawn flag, one numeric a forest.
forestParamsOf <- function(model) {
  lapply(model@forests, function(record) {
    as.double(c(
      record@n.trees,
      record@base,
      record@power,
      if (record@held.amplitude) 0 else 1
    ))
  })
}
