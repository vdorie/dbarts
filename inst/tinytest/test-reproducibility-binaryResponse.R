# Seeded-drift tripwire, not a correctness test: pins draws from a fixed
# seed so unintended RNG-affecting changes are caught. After an intentional
# shift, run tools/regenerate-snapshots.R and eyeball that the new values
# move by a plausible magnitude.

# The pinned values only mean anything on the reference build: its draw path
# is scalar and fixed-order, where the shipped build's is free to reassociate
# and lands elsewhere from the same seed. Regenerate on the reference build
# too - tools/regenerate-snapshots.R evaluates this file top to bottom, and
# the exit below is not a function it defines.
if (!identical(dbarts:::buildInfo()$mode, "reference")) {
  exit_file("seeded-drift snapshots are pinned to the reference build")
}

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# basic probit example
n.burn <- 10L
n.sims <- 100L

set.seed(99)
bartFit <- dbarts::bart(
  y.train = testData$Z,
  x.train = testData$X,
  ntree = 50L,
  ndpost = n.sims,
  nskip = n.burn,
  k = 4.5,
  verbose = FALSE
)

referenceBase <- list(
  yhatTrain = c(
    -0.0644582062346995,
    0.198293939988424,
    0.947723301886893,
    0.438016094709192,
    -0.410785945473055
  ),
  varcount = c(24L, 24L, 24L)
)

expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceBase$yhatTrain)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims, ], referenceBase$varcount)
expect_equal(extract(bartFit), pnorm(bartFit$yhat.train))
rm(bartFit, n.sims, n.burn)

# basic probit example with a non-zero binary offset
n.burn <- 10L
n.sims <- 100L

set.seed(99)
bartFit <- dbarts::bart(
  y.train = testData$Z,
  x.train = testData$X,
  ntree = 50L,
  ndpost = n.sims,
  nskip = n.burn,
  k = 4.5,
  binaryOffset = 0.1,
  verbose = FALSE
)

n.sims <- nrow(bartFit$yhat.train)

referenceOffset <- list(
  yhatTrain = c(
    0.0401438230497843,
    0.520943287890571,
    0.188529541776154,
    0.303160260166396,
    0.182047452939775
  ),
  varcount = c(19L, 16L, 23L)
)

expect_equal(bartFit$yhat.train[n.sims, 1L:5L], referenceOffset$yhatTrain)
expect_null(bartFit$yhat.test)
expect_equal(bartFit$varcount[n.sims, ], referenceOffset$varcount)
rm(bartFit, n.sims, n.burn)

rm(testData)
