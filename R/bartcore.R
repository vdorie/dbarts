# The bartcore engine behind dbartsSampler. The dbartsSampler methods
# delegate to the bartcoreSampler* functions below; the C side borrows
# vectors and pins them in the external pointer's protection slot.
#
# Not supported (methods error): weights with binary responses (weighted
# probit has no coherent latent-variable form), and setControl changes to
# anything fixed at creation (chain/tree counts, generators, and the cut
# grid).
#
# keepTrees/getTrees/predict use test-pinned formats. State serialization
# (storeState/setState, or runs with updateState) produces an engine-specific
# opaque object; restoring it into a sampler over the same data - including
# the transparent re-creation getPointer performs after save/load - continues
# the chains bitwise identically.

# A public dbartsSampler created through the treatment = spec branch
# (docs/design/bcf.md) carries its 0/1 treatment vector on data@treatment,
# mirroring data@weights; this is the R5-layer capability probe, cheaper than
# a round trip through the bridge's own (Chain::bcfGlue-based) one.
isBCFSampler <- function(sampler) {
  !is.null(sampler$data@treatment)
}

# BCF-specific wording for a mutation the bridge refuses through a guard
# shared with every multi-forest model (refuseMultiForestMutation and its
# siblings in R_interface_bartcore.cpp also cover a future multinomial
# creation route, so their own message cannot name BCF by itself). Raised
# R-side, before the .Call, so a BCF sampler never reaches the bridge's
# generic "multi-forest" phrasing.
refuseBCFMutation <- function(sampler, what, ...) {
  if (isBCFSampler(sampler)) {
    stop(what, " does not support a BCF sampler: ", ...)
  }
}

# Drives a dbartsSampler (the R5 sampler layer), reading its control defaults
# and delegating through its external pointer; cf. bartcoreRun, which drives a
# low-level bartcore handle directly.
bartcoreSamplerRun <- function(sampler, numBurnIn, numSamples) {
  control <- sampler$control
  if (is.na(numBurnIn)) {
    numBurnIn <- control@n.burn
  }
  if (is.na(numSamples)) {
    numSamples <- control@n.samples
  }
  if (is.na(numSamples)) {
    stop("bartcore engine samplers require 'numSamples' to be specified")
  }

  result <- .Call(
    C_dbarts_bartcore_run,
    sampler$getPointer(),
    as.integer(numBurnIn),
    as.integer(numSamples)
  )
  if (is.null(result)) {
    return(invisible(NULL))
  }
  result
}

bartcoreSamplerSetPredictor <- function(
  sampler,
  x,
  column,
  forceUpdate,
  updateCutPoints
) {
  # A sparse design - a pure dgCMatrix or a mixed dense/sparse container -
  # accepts column-granular and whole-matrix mutation
  # (docs/design/sparse-columns.md), maintained R-side by
  # installPredictorColumns rather than by the dense branch's pointer swap;
  # only per-observation replacement of a sparse-backed column stays fixed at
  # creation. Read before data@x is swapped.
  sparseSource <- predictorSourceIsSparse(sampler$data@x)

  # no BCF pre-check on the partial path either: the session's cell guard
  # caches every forest, pruned to the trees the column can move, so a row
  # installs only if it empties no leaf anywhere and a two-forest sampler takes
  # it (docs/plans/multiforest-predictor-mutation.md)
  partialUpdate <- !is.null(forceUpdate) &&
    is.character(forceUpdate) &&
    length(forceUpdate) == 1L &&
    !is.na(forceUpdate) &&
    forceUpdate == "partial"

  if (!is.null(column) && is.character(column)) {
    if (is.null(colnames(sampler$data@x))) {
      stop(
        "column names not specified at initialization, so cannot be ",
        "replaced by name"
      )
    }
    column <- match(column, colnames(sampler$data@x))
    if (anyNA(column)) stop("column name not found in names of current X")
  }

  ptr <- sampler$getPointer()

  if (partialUpdate) {
    if (is.null(column)) {
      stop("partial updates require a single 'column' to be specified")
    }
    if (length(column) != 1L) {
      stop("partial updates can only be applied to a single column")
    }
    column <- as.integer(column)
    # a CSC-backed column's rank storage cannot take a cell-at-a-time write
    # without an O(nnz) shift per cell; a DENSE-backed column of a mixed design
    # can, and is the motivating IRT latent case (docs/design/sparse-columns.md
    # extension (i)), so the refusal is per column rather than per design
    if (predictorColumnIsSparseBacked(sampler$data@x, column)) {
      stop(
        "per-observation updates require a dense-backed column; replace a ",
        "sparse column wholesale with a non-partial update"
      )
    }
    if (isTRUE(coerceOrError(updateCutPoints, "logical"))) {
      stop("partial updates cannot also update cut points")
    }

    x <- as.double(x)
    installed <- .Call(
      C_dbarts_bartcore_updatePredictorPerObservation,
      ptr,
      x,
      as.integer(column)
    )
    # the engine keeps no predictor matrix, so maintain data@x R-side for the
    # observations the scan installed; install by reference - the merge
    # starts from the old column and overwrites only the installed rows
    sampler$data@x <- installPredictorColumns(
      sampler$data@x,
      installed,
      column,
      x[installed]
    )
    return(installed)
  }

  forceUpdate <- if (is.null(forceUpdate)) {
    is.null(column)
  } else {
    coerceOrError(forceUpdate, "logical")
  }
  updateCutPoints <- coerceOrError(updateCutPoints, "logical")

  # no BCF pre-check here: a transactional whole-matrix or column update
  # revalidates every forest and rolls the whole change back if any leaf of any
  # tree of any forest would empty, so a two-forest sampler takes it
  # (docs/plans/multiforest-predictor-mutation.md). The per-observation session
  # above is the one that still refuses.

  # dim(), not is.matrix(): the latter is FALSE for every Matrix class, so a
  # transposed dgCMatrix argument (same total length, wrong shape) fell
  # through to the length-only check below and was silently reinterpreted
  # column-major
  xDim <- dim(x)
  if (is.null(column)) {
    if (!is.null(xDim)) {
      if (xDim[2L] != ncol(sampler$data@x)) {
        stop("dimension of x must be equal to ", ncol(sampler$data@x))
      }
      if (xDim[1L] != nrow(sampler$data@x)) {
        stop("dimension of x must be equal to ", nrow(sampler$data@x))
      }
    } else if (length(x) != prod(dim(sampler$data@x))) {
      stop("length of new x does not match old")
    }
    # a sparse-valued argument onto a sparse-backed design rides to the bridge
    # as supplied: it materializes there, under the store's own implicit rule,
    # rather than being densified here. Every other argument - a plain vector,
    # a sparseVector, any Matrix class the bridge does not ingest - keeps the
    # as.double path, as does a plain-matrix design.
    if (!(sparseSource && predictorSourceIsSparse(x))) {
      x <- if (!is.null(xDim)) {
        matrix(as.double(x), xDim[1L])
      } else {
        matrix(as.double(x), nrow(sampler$data@x))
      }
    }
    if (!sparseSource) {
      # a pointer swap: the engine borrows data@x, so install there first and
      # revert if the transaction rolls back
      oldX <- sampler$data@x
      sampler$data@x <- x
      tryResult <- tryCatch(
        updateSuccessful <- .Call(
          C_dbarts_bartcore_setPredictor,
          ptr,
          sampler$data@x,
          forceUpdate,
          updateCutPoints
        ),
        error = function(e) {
          sampler$data@x <- oldX
          e
        }
      )
      if (inherits(tryResult, "error")) {
        stop(tryResult)
      }
      if (!forceUpdate && !updateSuccessful) sampler$data@x <- oldX
    } else {
      # a sparse-bearing source: the engine borrows the argument matrix rather
      # than data@x, so nothing needs installing until it accepts. Splice the
      # replacement columns into the container BEFORE the call, so a throw
      # there cannot leave data@x describing the old design (sampler
      # re-creation after save/load reads it). Replacing a sparse column
      # densifies its storage - every row now differs from the implicit value.
      newX <- installPredictorColumns(
        sampler$data@x,
        NULL,
        seq_len(ncol(sampler$data@x)),
        x
      )
      updateSuccessful <- .Call(
        C_dbarts_bartcore_setPredictor,
        ptr,
        x,
        forceUpdate,
        updateCutPoints
      )
      if (isTRUE(updateSuccessful)) sampler$data@x <- newX
    }
  } else {
    column <- as.integer(column)
    if (any(column < 1L | column > ncol(sampler$data@x))) {
      stop(
        "column '",
        column[which(column < 1L | column > ncol(sampler$data@x))[1L]],
        "' is out of range"
      )
    }
    # length() counts a container's fields rather than its cells, so the shape
    # check reads dim() wherever the argument carries one; only a dimensionless
    # argument falls back to the total-length check
    if (!is.null(xDim)) {
      if (xDim[2L] != length(column)) {
        stop(
          "number of columns of new x does not match length of columns to ",
          "replace"
        )
      }
      if (xDim[1L] != nrow(sampler$data@x)) {
        stop("length of new x does not match y")
      }
    } else if (length(x) != nrow(sampler$data@x) * length(column)) {
      stop("length of new x does not match y")
    }
    if (!(sparseSource && predictorSourceIsSparse(x))) {
      x <- as.double(x)
    }
    # the engine keeps no predictor matrix, so maintain data@x R-side when the
    # update is applied (forceUpdate, or a non-rolled-back transaction);
    # install by reference - only the addressed columns change, the rest of the
    # container is shared. Build the replacement BEFORE the engine commits: a
    # CSC-backed column's install rewrites the container's sparse slots, and a
    # throw there after acceptance would leave data@x describing the old design
    # (sampler re-creation after save/load reads it).
    newX <- installPredictorColumns(sampler$data@x, NULL, column, x)
    updateSuccessful <- .Call(
      C_dbarts_bartcore_updatePredictor,
      ptr,
      x,
      column,
      forceUpdate,
      updateCutPoints
    )
    if (isTRUE(updateSuccessful)) {
      sampler$data@x <- newX
    }
  }

  if (!forceUpdate) updateSuccessful else invisible(NULL)
}

bartcoreSamplerSetResponse <- function(sampler, y, updateScale = FALSE) {
  y <- as.double(y)
  if (anyNA(y)) {
    stop("response contains missing values")
  }
  if (isTRUE(updateScale)) {
    refuseBCFMutation(
      sampler,
      "setResponse(updateScale = TRUE)",
      "both forests keep leaf calibrations stated against the response ",
      "transform fixed at creation; use updateScale = FALSE instead"
    )
  }
  # validate (the C length check) before installing, so a rejected y never
  # leaves data@y holding the bad replacement
  .Call(
    C_dbarts_bartcore_setResponse,
    sampler$getPointer(),
    y,
    updateScale
  )
  sampler$data@y <- y
  invisible(NULL)
}

bartcoreSamplerSetOffset <- function(sampler, offset, updateScale) {
  if (isTRUE(updateScale)) {
    refuseBCFMutation(
      sampler,
      "setOffset(updateScale = TRUE)",
      "both forests keep leaf calibrations stated against the response ",
      "transform fixed at creation; use updateScale = FALSE instead"
    )
  }
  # a synced test offset follows the regular one;
  # NA marks "leave the test offset alone"
  offset.test <- NA
  if (is.null(offset)) {
    if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
      offset.test <- NULL
    }
  } else {
    offset <- as.double(offset)
    if (anyNA(offset)) {
      stop("'offset' contains missing values")
    }
    if (length(offset) == 1L) {
      if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
        offset.test <- if (!is.null(sampler$data@x.test)) {
          rep_len(offset, nrow(sampler$data@x.test))
        } else {
          NULL
        }
      }
      offset <- rep_len(offset, length(sampler$data@y))
    } else {
      if (length(offset) != length(sampler$data@y)) {
        stop(
          "length of replacement offset is not equal to number of observations"
        )
      }
      if (identical(sampler$data@testUsesRegularOffset, TRUE)) {
        offset.test <- if (
          !is.null(sampler$data@x.test) &&
            length(offset) == nrow(sampler$data@x.test)
        ) {
          offset
        } else {
          NULL
        }
      }
    }
  }

  ptr <- sampler$getPointer()

  sampler$data@offset <- offset
  .Call(
    C_dbarts_bartcore_setOffset,
    ptr,
    sampler$data@offset,
    as.logical(updateScale)
  )

  if (!identical(offset.test, NA)) {
    oldOffset.test <- sampler$data@offset.test
    sampler$data@offset.test <- offset.test
    tryResult <- tryCatch(
      .Call(C_dbarts_bartcore_setTestOffset, ptr, sampler$data@offset.test),
      error = function(e) {
        sampler$data@offset.test <- oldOffset.test
        e
      }
    )
    if (inherits(tryResult, "error")) stop(tryResult)
  }

  invisible(NULL)
}

bartcoreSamplerSetData <- function(sampler, newData) {
  if (!inherits(newData, "dbartsData")) {
    stop("'data' must inherit from dbartsData")
  }
  refuseBCFMutation(
    sampler,
    "setData",
    "both forests are calibrated against the data at creation; make a new ",
    "sampler instead"
  )
  if (ncol(newData@x) != ncol(sampler$data@x)) {
    stop("bartcore setData requires the same predictors")
  }

  newData@n.cuts <- sampler$data@n.cuts
  newData@sigma <- sampler$data@sigma

  ptr <- sampler$getPointer()

  oldData <- sampler$data
  sampler$data <- newData
  tryResult <- tryCatch(
    .Call(C_dbarts_bartcore_setData, ptr, sampler$data),
    error = function(e) {
      sampler$data <- oldData
      e
    }
  )
  if (inherits(tryResult, "error")) {
    stop(tryResult)
  }

  invisible(NULL)
}

bartcoreSamplerSetCutPoints <- function(sampler, cuts, column) {
  if (!is.null(column) && is.character(column)) {
    if (is.null(colnames(sampler$data@x))) {
      stop(
        "column names not specified at initialization, so cannot be ",
        "replaced by name"
      )
    }
    column <- match(column, colnames(sampler$data@x))
    if (anyNA(column)) stop("column name not found in names of current X")
  }
  if (is.null(column)) {
    column <- seq_len(ncol(sampler$data@x))
  }

  if (!is.list(cuts)) {
    cuts <- list(cuts)
  }
  cuts <- lapply(cuts, as.double)

  # the engine re-quantizes the transient borrow of the current predictors
  .Call(
    C_dbarts_bartcore_setCutPoints,
    sampler$getPointer(),
    cuts,
    as.integer(column),
    rawPredictorMatrix(sampler$data@x)
  )
  invisible(NULL)
}

bartcoreSamplerSetTestPredictor <- function(sampler, x.test, column) {
  if (!is.null(column) && is.character(column)) {
    if (is.null(colnames(sampler$data@x.test))) {
      stop(
        "column names not specified at initialization, so cannot be ",
        "replaced by name"
      )
    }
    column <- match(column, colnames(sampler$data@x.test))
    if (anyNA(column)) {
      stop("column name not found in names of current test predictor matrix")
    }
  }

  if (is.null(column)) {
    # NULL removes the test data; a frame/sparse input becomes a container the
    # bridge codes against the training cuts. The bridge clears any test offset
    # with a NULL removal.
    x.test <- validateXTest(x.test, sampler$data@x)
  } else {
    column <- as.integer(column)
    if (any(column < 1L | column > ncol(sampler$data@x.test))) {
      stop(
        "column '",
        column[which(column < 1L | column > ncol(sampler$data@x.test))[1L]],
        "' is out of range"
      )
    }
    xTestDim <- dim(x.test)
    if (!is.null(xTestDim) && xTestDim[2L] != length(column)) {
      stop(
        "number of columns of new x does not match length of columns to replace"
      )
    }
    if (length(x.test) != nrow(sampler$data@x.test) * length(column)) {
      stop("length of new x does not match old x.test")
    }
    if (inherits(sampler$data@x.test, "dbartsMixedMatrix")) {
      # a container's per-column storage decision (dense vs CSC-backed) is
      # preserved; installPredictorColumns splices the replacement in place,
      # canonicalizing a CSC-backed target column against its implicit
      x.test <- installPredictorColumns(
        sampler$data@x.test,
        NULL,
        column,
        x.test
      )
    } else {
      # the engine replaces the whole matrix; column updates copy-modify it
      new.x.test <- sampler$data@x.test
      new.x.test[, column] <- as.double(x.test)
      x.test <- new.x.test
    }
  }

  # install the new test set R-side, then roll back if the bridge refuses it
  # (a container whose leaf-covariate column is CSC-backed), keeping the R5
  # object and the engine's prior test store consistent
  oldX.test <- sampler$data@x.test
  oldOffset.test <- sampler$data@offset.test
  sampler$data@x.test <- x.test
  if (is.null(x.test)) {
    sampler$data@offset.test <- NULL
  }
  tryResult <- tryCatch(
    .Call(
      C_dbarts_bartcore_setTestPredictor,
      sampler$getPointer(),
      sampler$data@x.test
    ),
    error = function(e) {
      sampler$data@x.test <- oldX.test
      sampler$data@offset.test <- oldOffset.test
      e
    }
  )
  if (inherits(tryResult, "error")) {
    stop(tryResult)
  }
  invisible(NULL)
}

# Internal interface used by the tests and the equivalence harness; a
# bartcore handle is constructed from the validated control/model/data of an
# existing dbartsSampler regardless of that sampler's own engine. family
# selects the response model for binary responses ("probit", the default,
# or "logistic"); it is the only access to the logistic sampler until the
# new public surface lands. Likewise, columns marked CATEGORICAL (1L) in
# data@varTypes split by category subset rather than by threshold - a
# bartcore-only capability with no public ingestion path yet, so callers
# flip the type by hand; values must be integer codes 0..K-1, K <= 32.

bartcoreSampler <- function(sampler, family = "") {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_create,
    sampler$control,
    sampler$model,
    sampler$data,
    as.character(family)
  )
  # the engine keeps no predictor matrix, so the low-level wrappers track the
  # current predictors R-side to feed the re-quantize entry points
  # (setCutPoints, setState, saved-tree getTrees); null for sparse designs
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# A built predictor store (cuts + codes) shared across row-subset samplers;
# the shared-handle design (public-surface.md section 5), internal and
# unserializable. control
# contributes useQuantiles; data contributes x, the column types, and
# n.cuts. leafCovariateColumns names (1-based) the columns whose raw values
# a view's leaf model will read; the handle owns raw only for those, so a
# constant-leaf caller passes none. A view designating an undeclared column
# is refused when it is built.
bartcoreDataHandle <- function(control, data, leafCovariateColumns = NULL) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createDataHandle,
    control,
    data,
    if (!is.null(leafCovariateColumns)) as.integer(leafCovariateColumns)
  )
  result
}

# A sampler over a row subset of a handle: it copies the handle's cut grid
# and gathers its rows' codes, so folds bin identically to the full data.
# data is the full data object the handle was built from; trainRows and
# testRows index its rows, y/weights/offset are sliced by trainRows, and a
# test offset comes from offset[testRows] (xbart's fold semantics). columns,
# when given, are 1-based indices restricting the view to a column subset (NULL
# spans every column). The result refuses raw-predictor mutation (setPredictor
# and friends, setData, setCutPoints, setState); family is as bartcoreSampler's.
bartcoreSamplerFromHandle <- function(
  handle,
  control,
  model,
  data,
  trainRows,
  testRows = NULL,
  family = "",
  columns = NULL
) {
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createFromHandle,
    control,
    model,
    data,
    handle$ptr,
    as.integer(trainRows),
    if (!is.null(testRows)) as.integer(testRows),
    as.character(family),
    if (!is.null(columns)) as.integer(columns)
  )
  # a view refuses the raw-predictor and re-quantize surface, and its live-tree
  # getTrees needs no training replay, so it carries no predictor matrix
  result$x <- NULL
  result
}

# A BCF two-forest sampler (docs/design/bcf.md), internal and gaussian only:
# the model spec is the prognostic forest mu(x, pihat) - the caller supplies
# pihat as an ordinary predictor column - the arguments below the treatment
# forest tau(x) and the glue. z is the 0/1 treatment. sd.control and
# sd.moderate are the prognostic and effect magnitudes in sd(y) units; the
# calibration map converts them to the per-forest leaf scales at creation,
# overriding the host model's node prior and k for the mu forest. update.a /
# update.b hold the matching glue block fixed when FALSE. moderators restricts
# the treatment forest to a subset of the shared design's columns (a character
# vector matched to colnames(data@x), or a 1-based numeric index vector; NULL =
# all columns); the prognostic forest always reads the full store. In the BCF
# of Hahn, Murray and Carvalho (2020) the estimated propensity score enters
# mu's design but not tau's - with a shared design that means carrying pihat as
# a data column and leaving it out of moderators.
bartcoreBCFSampler <- function(
  sampler,
  z,
  n.trees.treatment = 50L,
  treatment.base = 0.25,
  treatment.power = 3,
  sd.control = 2,
  sd.moderate = 1,
  b.prior.variance = 0.5,
  update.a = TRUE,
  update.b = TRUE,
  moderators = NULL,
  mu.interactions = NULL,
  tau.interactions = NULL,
  mu.blocks = NULL,
  tau.blocks = NULL
) {
  moderators <- resolveModerators(moderators, sampler$data)
  bcfParams <- as.double(c(
    n.trees.treatment,
    treatment.base,
    treatment.power,
    sd.control,
    sd.moderate,
    b.prior.variance,
    update.a,
    update.b
  ))
  # per-forest interaction constraints (docs/design/interaction-constraints.md):
  # mu and tau each resolve their own interactions() prior against the shared
  # design, so the treatment forest can be capped additive-or-low-order while
  # the prognostic forest stays free (the calibrated-additivity causal use)
  muInteractions <- resolveInteractions(mu.interactions, sampler$data)
  tauInteractions <- resolveInteractions(tau.interactions, sampler$data)

  # per-forest block-additive constraints (each whole tree confined to one
  # declared column group): mu partitions the full
  # design over its own tree count; tau partitions its available columns (the
  # moderator subset if restricted, else the full design) over the treatment
  # tree count, and the engine intersects tau's block rows with the moderator
  # mask at install. The deterministic capacity consumes no rng.
  muBlocks <- resolveBlocks(mu.blocks, sampler$data, sampler$control@n.trees)
  tauBlocks <- resolveBlocks(
    tau.blocks,
    sampler$data,
    as.integer(n.trees.treatment),
    availableColumns = moderators
  )

  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createBCF,
    sampler$control,
    sampler$model,
    sampler$data,
    as.double(z),
    bcfParams,
    # resolved 1-based moderator indices, or NULL for an unrestricted forest
    moderators,
    # resolved interactions() lists (max.order + 0-based forbidden pairs), or
    # NULL; mu is forest 0 (prognostic), tau is forest 1 (treatment)
    muInteractions,
    tauInteractions,
    # resolved blocks() lists (0-based per-column group + per-group tree
    # capacity), or NULL; per-forest as with the interactions above
    muBlocks,
    tauBlocks
  )
  # BCF requires dense predictors; track them R-side for the re-quantize surface
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# A K-forest multinomial (softmax) sampler (docs/design/multinomial.md),
# internal and single-trial: the K symmetric category forests couple through a
# softmax likelihood with an interleaved one-vs-rest Polya-Gamma augmentation
# and a level-centering move. labels are the 0-based category codes (0..K-1),
# one per observation - the response, so the host sampler's own response is
# ignored. K defaults to one past the largest code. The per-forest leaf scale
# follows the multinomial calibration (pi*sqrt(3)/sqrt(2)), not the host node
# prior. The run's train channel carries the K softmax probabilities
# (n.observations x K x n.samples); per-category fits and split counts read
# through bartcoreForestFits / bartcoreForestVariableCounts (0-based forest =
# category). When the data object carries x.test, the run's test channel
# reports the same K softmax probabilities on the held-out rows (the K forests'
# totalTestFits blended by softmax). Under keepTrees, out-of-sample predict()
# replays all K forests' saved trees and softmaxes them (bartcorePredict).
bartcoreMultinomialSampler <- function(sampler, labels, K = NULL) {
  labels <- as.integer(labels)
  if (anyNA(labels)) {
    stop("multinomial labels must be integer category codes 0..K-1")
  }
  if (any(labels < 0L)) {
    stop("multinomial labels must be nonnegative category codes 0..K-1")
  }
  if (is.null(K)) {
    K <- max(labels) + 1L
  }
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createMultinomial,
    sampler$control,
    sampler$model,
    sampler$data,
    labels,
    as.integer(K)
  )
  # the engine keeps no predictor matrix; track it R-side for the re-quantize
  # surface, as the other dense-predictor wrappers do
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# The grouped-count analog of bartcoreMultinomialSampler (docs/design/
# multinomial.md): the response is an n x K matrix of nonnegative integer
# counts, column k holding category k's success counts per observation, with
# trials n_i = sum_k counts[i, k] (>= 1). K defaults to the column count. Same
# K-forest softmax engine; the single-trial label path is the special case of a
# one-hot matrix with every trial 1. Validated R-side (safe over fast); the
# engine re-derives the trials and re-checks the invariants.
bartcoreMultinomialCountSampler <- function(sampler, counts, K = NULL) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) {
    stop("multinomial counts must be a numeric matrix of nonnegative integers")
  }
  if (anyNA(counts)) {
    stop("multinomial counts must not contain missing values")
  }
  if (any(counts < 0)) {
    stop("multinomial counts must be nonnegative")
  }
  if (any(counts != round(counts))) {
    stop("multinomial counts must be whole numbers")
  }
  if (is.null(K)) {
    K <- ncol(counts)
  }
  if (ncol(counts) != K) {
    stop("multinomial counts must have K columns")
  }
  if (K < 2L) {
    stop("multinomial requires at least two categories")
  }
  if (any(rowSums(counts) < 1)) {
    stop("every multinomial count row must have at least one trial (n_i >= 1)")
  }
  storage.mode(counts) <- "integer"
  result <- new.env(parent = emptyenv())
  # counts is column-major (category-major), the layout the combiner reads
  result$ptr <- .Call(
    C_dbarts_bartcore_createMultinomialCounts,
    sampler$control,
    sampler$model,
    sampler$data,
    counts,
    as.integer(K)
  )
  # the engine keeps no predictor matrix; track it R-side for the re-quantize
  # surface, as the other dense-predictor wrappers do
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# Replaces a multinomial sampler's response: an n x K matrix of nonnegative
# integer counts in the same layout the count creation entry takes (column k is
# category k), with trials n_i = sum_k counts[i, k] >= 1. n and K are fixed at
# creation - every combiner buffer is sized by n, and K is the forest count - so
# only the values may change. The trees carry over, fitted to the previous
# counts, exactly as a single-forest setResponse leaves them; the next run forms
# every category's working response against the new counts.
#
# The sweep draws n_i Polya-Gamma variates per observation per category, so
# replacing single-trial labels with grouped counts multiplies sweep cost by
# mean(n_i).
#
# Validated R-side (safe over fast); the engine re-derives the trials and
# re-checks the invariants, and refuses everything before it installs anything,
# so a rejected matrix leaves the sampler untouched.
bartcoreSetCounts <- function(bcSampler, counts) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) {
    stop("multinomial counts must be a numeric matrix of nonnegative integers")
  }
  if (anyNA(counts)) {
    stop("multinomial counts must not contain missing values")
  }
  if (any(counts < 0)) {
    stop("multinomial counts must be nonnegative")
  }
  if (any(counts != round(counts))) {
    stop("multinomial counts must be whole numbers")
  }
  if (any(rowSums(counts) < 1)) {
    stop("every multinomial count row must have at least one trial (n_i >= 1)")
  }
  storage.mode(counts) <- "integer"
  # counts is column-major (category-major), the layout the combiner reads
  invisible(.Call(C_dbarts_bartcore_setCounts, bcSampler$ptr, counts))
}

# The 0/1 treatment the treatment forest contrasts on; re-forms b_{z_i} and
# both residuals on the next run.
bartcoreSetTreatment <- function(bcSampler, z) {
  invisible(.Call(C_dbarts_bartcore_setTreatment, bcSampler$ptr, as.double(z)))
}

# A per-forest, per-observation weight s on forest `forest` (0 prognostic, 1
# treatment): a multiplicative precision factor on that forest's own leaf
# conditionals, composing with the case weights so the forest's draws see
# w_i * m_f^2 * s_i. It does not remove the row from occupancy, from the
# empty-leaf veto, from the combination (the row still receives m_f f_f(x_i)),
# or from the residual sigma degrees of freedom; s_i = 0 says only that row i
# carries no information about forest f, and that forest's leaves over such
# rows stay well-defined prior draws.
#
# Two edges, or a consumer is misled. At s_i = 0 with a nonzero multiplier only
# the WEIGHT is zeroed - the response stays the reparameterized residual - so
# the reported-fit exactness an exactly zero multiplier buys does not follow
# this channel. And the weight lives on the sampler, not in its state, so a
# pipeline that REBUILDS a sampler and restores a stored state silently drops
# the weight and fits a different model while the states still agree.
bartcoreSetForestWeights <- function(bcSampler, forest, weights) {
  weights <- as.double(weights)
  if (!all(is.finite(weights)) || any(weights < 0)) {
    stop("forest weights must be finite and non-negative")
  }
  invisible(.Call(
    C_dbarts_bartcore_setForestWeights,
    bcSampler$ptr,
    as.integer(forest),
    weights
  ))
}

# The glue on the combining response, a 3 x n.chains matrix of (a, b0, b1).
bartcoreBCFGlue <- function(bcSampler) {
  .Call(C_dbarts_bartcore_getBCFGlue, bcSampler$ptr)
}

# A forest's internal-scale function values (0 prognostic, 1 treatment),
# n.observations x n.chains.
bartcoreForestFits <- function(bcSampler, forest) {
  .Call(C_dbarts_bartcore_getForestFits, bcSampler$ptr, as.integer(forest))
}

# A forest's current per-predictor split counts (0 prognostic, 1 treatment),
# n.predictors x n.chains; the per-forest analog of the run's varcount channel.
bartcoreForestVariableCounts <- function(bcSampler, forest) {
  .Call(
    C_dbarts_bartcore_getForestVariableCounts,
    bcSampler$ptr,
    as.integer(forest)
  )
}

bartcoreSetModel <- function(bcSampler, model, control, data) {
  invisible(.Call(
    C_dbarts_bartcore_setModel,
    bcSampler$ptr,
    model,
    control,
    data
  ))
}

# Drives a low-level bartcore handle (a bcSampler env holding $ptr) directly,
# not a dbartsSampler; cf. bartcoreSamplerRun, the R5 sampler-layer entry.
bartcoreRun <- function(bcSampler, numBurnIn = 0L, numSamples = 1L) {
  .Call(
    C_dbarts_bartcore_run,
    bcSampler$ptr,
    as.integer(numBurnIn),
    as.integer(numSamples)
  )
}

bartcoreSetOffset <- function(bcSampler, offset, updateScale = FALSE) {
  if (!is.null(offset)) {
    offset <- as.double(offset)
  }
  invisible(.Call(
    C_dbarts_bartcore_setOffset,
    bcSampler$ptr,
    offset,
    as.logical(updateScale)
  ))
}

bartcoreSetResponse <- function(bcSampler, y, updateScale = FALSE) {
  invisible(.Call(
    C_dbarts_bartcore_setResponse,
    bcSampler$ptr,
    as.double(y),
    as.logical(updateScale)
  ))
}

bartcoreSetWeights <- function(bcSampler, weights) {
  invisible(.Call(
    C_dbarts_bartcore_setWeights,
    bcSampler$ptr,
    as.double(weights)
  ))
}

bartcoreSetTestOffset <- function(bcSampler, offset.test) {
  if (!is.null(offset.test)) {
    offset.test <- as.double(offset.test)
  }
  invisible(.Call(C_dbarts_bartcore_setTestOffset, bcSampler$ptr, offset.test))
}

bartcoreSetData <- function(bcSampler, data) {
  invisible(.Call(C_dbarts_bartcore_setData, bcSampler$ptr, data))
}

bartcoreSetTestPredictor <- function(bcSampler, x.test) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  invisible(.Call(C_dbarts_bartcore_setTestPredictor, bcSampler$ptr, x.test))
}

bartcoreSetPredictor <- function(
  bcSampler,
  x,
  forceUpdate = FALSE,
  updateCutPoints = FALSE
) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  result <- .Call(
    C_dbarts_bartcore_setPredictor,
    bcSampler$ptr,
    x,
    as.logical(forceUpdate),
    as.logical(updateCutPoints)
  )
  # track the current predictors R-side (the engine keeps none) when installed
  if (isTRUE(result)) {
    bcSampler$x <- x
  }
  result
}

# Overwrite columns of the current predictors; the engine keeps no matrix, so
# the wrapper maintains the R-side copy the re-quantize entry points read.
bartcoreUpdatePredictor <- function(
  bcSampler,
  x,
  columns,
  forceUpdate = FALSE,
  updateCutPoints = FALSE
) {
  columns <- as.integer(columns)
  result <- .Call(
    C_dbarts_bartcore_updatePredictor,
    bcSampler$ptr,
    as.double(x),
    columns,
    as.logical(forceUpdate),
    as.logical(updateCutPoints)
  )
  if (isTRUE(result) && !is.null(bcSampler$x)) {
    bcSampler$x[, columns] <- matrix(as.double(x), nrow(bcSampler$x))
  }
  result
}

bartcoreUpdatePredictorPerObservation <- function(bcSampler, x, column) {
  x <- as.double(x)
  column <- as.integer(column)
  installed <- .Call(
    C_dbarts_bartcore_updatePredictorPerObservation,
    bcSampler$ptr,
    x,
    column
  )
  if (!is.null(bcSampler$x)) {
    bcSampler$x[installed, column] <- x[installed]
  }
  installed
}

bartcoreUpdatePredictorPerObservationJointly <- function(
  bcSamplers,
  x,
  columns
) {
  x <- as.double(x)
  columns <- as.integer(columns)
  installed <- .Call(
    C_dbarts_bartcore_updatePredictorPerObservationJointly,
    lapply(bcSamplers, function(s) s$ptr),
    x,
    columns
  )
  for (k in seq_along(bcSamplers)) {
    if (!is.null(bcSamplers[[k]]$x)) {
      bcSamplers[[k]]$x[installed, columns[k]] <- x[installed]
    }
  }
  installed
}

# cutPoints is a list of strictly increasing numeric vectors, one per column;
# trees are refreshed unconditionally, collapsing splits the new cuts orphan
bartcoreSetCutPoints <- function(bcSampler, cutPoints, columns) {
  if (!is.list(cutPoints)) {
    cutPoints <- list(cutPoints)
  }
  cutPoints <- lapply(cutPoints, as.double)
  invisible(.Call(
    C_dbarts_bartcore_setCutPoints,
    bcSampler$ptr,
    cutPoints,
    as.integer(columns),
    bcSampler$x
  ))
}

bartcoreGetLatents <- function(bcSampler) {
  .Call(C_dbarts_bartcore_getLatents, bcSampler$ptr, NULL)
}

bartcorePredict <- function(bcSampler, x.test, offset.test = NULL) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  if (!is.null(offset.test)) {
    offset.test <- as.double(offset.test)
  }
  .Call(C_dbarts_bartcore_predict, bcSampler$ptr, x.test, offset.test)
}

bartcoreGetTrees <- function(
  bcSampler,
  chainNums,
  sampleNums = NULL,
  treeNums,
  current = FALSE,
  newdata = NULL,
  forest = 0L
) {
  if (!is.null(newdata)) {
    newdata <- as.matrix(newdata)
    storage.mode(newdata) <- "double"
  }
  # forest is 0-based, as for bartcoreForestFits (0 prognostic, 1 treatment)
  .Call(
    C_dbarts_bartcore_getTrees,
    bcSampler$ptr,
    as.integer(chainNums),
    if (is.null(sampleNums)) NULL else as.integer(sampleNums),
    as.integer(treeNums),
    as.logical(current),
    newdata,
    bcSampler$x,
    as.integer(forest)
  )
}

bartcoreStoreState <- function(bcSampler) {
  .Call(C_dbarts_bartcore_storeState, bcSampler$ptr)
}

bartcoreSetState <- function(bcSampler, state) {
  invisible(.Call(
    C_dbarts_bartcore_setState,
    bcSampler$ptr,
    state,
    bcSampler$x
  ))
}
