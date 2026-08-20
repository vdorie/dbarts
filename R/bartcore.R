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

# A public dbartsSampler created through the forests = spec branch
# (docs/design/bcf.md) carries its forests' amplitude bases on data@bases,
# mirroring data@weights; this is the R-level capability probe, cheaper than
# a round trip through the bridge's own (totalAmplitudes-based) one. It is a
# CAPABILITY test - "carries amplitudes" - and deliberately not a forest count:
# a K-forest multinomial carries several forests and no amplitudes at all, so a
# numForests probe would misfire on it, as the bridge and the flat C entry each
# record independently.
isBCFSampler <- function(sampler) {
  !is.null(sampler$data@bases)
}

# BCF-specific wording for a mutation the bridge refuses through a guard
# shared with every multi-forest model (refuseMultiForestMutation and its
# siblings in R_interface_bartcore.cpp also cover the multinomial creation
# route, so their own message cannot name BCF by itself). Raised R-side,
# before the .Call, so a BCF sampler never reaches the bridge's generic
# "multi-forest" phrasing.
refuseBCFMutation <- function(sampler, what, ...) {
  if (isBCFSampler(sampler)) {
    stop(what, " does not support a BCF sampler: ", ...)
  }
}

# Drives a dbartsSampler (the R-level sampler layer), reading its control
# defaults and delegating through its external pointer; cf. bartcoreRun,
# which drives a low-level bartcore handle directly.
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

# Resolves a character 'column' against source's colnames into a 1-based
# integer index (or indices); NULL or an already-numeric 'column' passes
# through unchanged. 'what' names source for the not-found message.
resolveColumnIndex <- function(source, column, what) {
  if (is.null(column) || !is.character(column)) {
    return(column)
  }
  if (is.null(colnames(source))) {
    stop(
      "column names not specified at initialization, so cannot be ",
      "replaced by name"
    )
  }
  column <- match(column, colnames(source))
  if (anyNA(column)) {
    stop("column name not found in names of ", what)
  }
  column
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
  # installs only if it empties no leaf anywhere and a two-forest sampler
  # takes it
  partialUpdate <- !is.null(forceUpdate) &&
    is.character(forceUpdate) &&
    length(forceUpdate) == 1L &&
    !is.na(forceUpdate) &&
    forceUpdate == "partial"

  column <- resolveColumnIndex(sampler$data@x, column, "current X")

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

  # no BCF pre-check here either: a transactional whole-matrix or column
  # update revalidates every forest and rolls the whole change back if any
  # leaf of any tree of any forest would empty, so a two-forest sampler takes
  # it - the same as the per-observation session above.

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
      stop("'x' must have length ", prod(dim(sampler$data@x)))
    }
    # a sparse-valued argument onto a sparse-backed design rides to the bridge
    # as supplied: it materializes there, under the store's own implicit rule,
    # rather than being densified here. Every other argument - a plain vector,
    # a sparseVector, any Matrix class the bridge does not ingest - keeps the
    # as.double path, as does a plain-matrix design.
    if (!(sparseSource && predictorSourceIsSparse(x))) {
      # matrix(as.double(x), ...) strips every attribute, so the incoming
      # dimnames would otherwise vanish; carry them onto the replacement,
      # falling back to the sampler's current names when x supplies none -
      # the shapes already agree by the checks above
      xDimnames <- dimnames(x)
      x <- if (!is.null(xDim)) {
        matrix(as.double(x), xDim[1L])
      } else {
        matrix(as.double(x), nrow(sampler$data@x))
      }
      dimnames(x) <- if (!is.null(xDimnames)) {
        xDimnames
      } else {
        dimnames(sampler$data@x)
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
        stop("'x' must have ", length(column), " column(s)")
      }
      if (xDim[1L] != nrow(sampler$data@x)) {
        stop("'x' must have ", nrow(sampler$data@x), " row(s)")
      }
    } else if (length(x) != nrow(sampler$data@x) * length(column)) {
      stop("'x' must have length ", nrow(sampler$data@x) * length(column))
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
      "every forest keeps its leaf calibration stated against the anchor ",
      "fixed at creation; use updateScale = FALSE instead"
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
      "every forest keeps its leaf calibration stated against the anchor ",
      "fixed at creation; use updateScale = FALSE instead"
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
    "every forest is calibrated against the data at creation; make a new ",
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
  column <- resolveColumnIndex(sampler$data@x, column, "current X")
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
  column <- resolveColumnIndex(
    sampler$data@x.test,
    column,
    "current test predictor matrix"
  )

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
      stop("'x.test' must have ", length(column), " column(s)")
    }
    if (length(x.test) != nrow(sampler$data@x.test) * length(column)) {
      stop(
        "'x.test' must have length ",
        nrow(sampler$data@x.test) * length(column)
      )
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
  # (a container whose leaf-covariate column is CSC-backed), keeping the
  # R-level object and the engine's prior test store consistent
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
# or "logistic"), the same choice bart2's own family argument exposes
# publicly. Columns marked CATEGORICAL (1L) in data@varTypes split by
# category subset rather than by threshold; factors = "categorical", the
# default on the public data constructor, sets this directly, so callers
# no longer need to flip the type by hand. Values must be integer codes
# 0..K-1, K <= 32.

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

# A BCF two-forest sampler (docs/design/bcf.md), internal:
# the model spec is the prognostic forest mu(x, pihat) - the caller supplies
# pihat as an ordinary predictor column - the arguments below the treatment
# forest tau(x) and the glue. z is the 0/1 treatment. sd.control and
# sd.moderate are the prognostic and effect magnitudes in units of the
# family's own latent scale (sd(y) under gaussian, 1 under probit); the
# calibration map converts them to the per-forest leaf scales at creation,
# overriding the host model's node prior and k for the mu forest. A NULL
# sd.control takes the family's own default median, 2 under gaussian and 1
# under the latent families, exactly as the public route does; sd.moderate's 1
# is the K-aware default's own K = 2 value, sqrt(2/2), so this two-forest
# spelling states it as the literal it has always been. update.a /
# update.b hold the matching glue block fixed when FALSE. moderators restricts
# the treatment forest to a subset of the shared design's columns (a character
# vector matched to colnames(data@x), or a 1-based numeric index vector; NULL =
# all columns); the prognostic forest always reads the full store. In the BCF
# of Hahn, Murray and Carvalho (2020) the estimated propensity score enters
# mu's design but not tau's - with a shared design that means carrying pihat as
# a data column and leaving it out of moderators.
#
# family names the response family the two forests are combined under
# (gaussian, probit or logistic), and is written into the model this call
# hands the bridge, which is the ONE place the family is read from. Supplying
# it here rather than only on the model keeps the family visible at the call
# site, as the other internal creators do, while leaving one source of truth.
bartcoreBCFSampler <- function(
  sampler,
  z,
  family = NULL,
  n.trees.treatment = 50L,
  treatment.base = 0.25,
  treatment.power = 3,
  sd.control = NULL,
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
  # the bridge derives the family from this model's own slot, so a supplied name
  # is written into the copy this call hands it; an unknown one is refused
  # there, against the response shape, as every other route's is. Resolved
  # BEFORE the transport below rather than beside the .Call, because
  # sd.control's default is the family's: a NULL left inside as.double(c(...))
  # collapses that vector to length 7 and surfaces as the bridge's length-8
  # refusal, naming a shape rather than the family.
  model <- sampler$model
  if (!is.null(family)) {
    model@family <- family
  }
  if (is.null(sd.control)) {
    # the same helper the public route takes (R/model.R), so the two routes
    # stay on one set of defaults - which is what the creation oracle's
    # draw-for-draw comparison exists to hold - and a family this route cannot
    # build still meets the bridge's own refusal rather than one from here
    sd.control <- defaultAmplitudePriorScale(model@family)
  }
  # the K-length per-forest transport, at K = 2: the prognostic forest carries
  # a plain scalar amplitude under the half-Cauchy scale mixture sd.control is
  # the median of, and the treatment forest a fixed-variance amplitude pair
  # whose magnitude rides its node scale through the half-normal median
  bcfParams <- list(
    as.double(c(0, 0, 0, 1, 1, 1, sd.control, update.a)),
    as.double(c(
      n.trees.treatment,
      treatment.base,
      treatment.power,
      sd.moderate,
      0.674,
      b.prior.variance,
      0,
      update.b
    ))
  )
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

  # the treatment forest's basis is the two-level factor basis of z, whose
  # amplitudes are exactly (b0, b1); the prognostic forest declares none, so
  # the engine gives it the implicit intercept its single amplitude scales
  z <- as.double(z != 0)
  bases <- list(NULL, cbind(1 - z, z))

  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createBCF,
    sampler$control,
    model,
    sampler$data,
    bases,
    bcfParams,
    # resolved 1-based column indices per forest, or NULL for unrestricted
    list(NULL, moderators),
    # resolved interactions() lists (max.order + 0-based forbidden pairs), or
    # NULL; forest 1 is prognostic, forest 2 the treatment forest
    list(muInteractions, tauInteractions),
    # resolved blocks() lists (0-based per-column group + per-group tree
    # capacity), or NULL; per-forest as with the interactions above
    list(muBlocks, tauBlocks)
  )
  # BCF requires dense predictors; track them R-side for the re-quantize surface
  result$x <- rawPredictorMatrix(sampler$data@x)
  result
}

# Integer coercion for an already-validated count matrix. The RANGE check comes
# BEFORE the coercion: storage.mode() turns a value past .Machine$integer.max
# into NA with a warning, which the engine would then report as a negative
# count - a true refusal naming the wrong reason.
asCountMatrix <- function(counts) {
  if (any(counts > .Machine$integer.max)) {
    stop("multinomial counts must be representable as integers")
  }
  storage.mode(counts) <- "integer"
  counts
}

# The n x K category offset shared by the two multinomial creators and the
# setter: NULL means none, anything else must be a numeric n x K matrix with
# every entry finite (an infinite entry propagates through the log-sum-exp
# margin into a NaN for every category of that row). Only the row-centred part
# of the offset is identified - adding a constant to every entry of a row
# leaves the softmax unchanged - so the input is passed through as given rather
# than silently re-centred. `what` names the offset in the refusals, since the
# train rows and the test rows each carry their own.
validateCategoryOffset <- function(offset, n, K, what = "category offset") {
  if (is.null(offset)) {
    return(NULL)
  }
  offset <- as.matrix(offset)
  if (!is.numeric(offset)) {
    stop(sprintf("multinomial %s must be a numeric matrix", what))
  }
  if (nrow(offset) != n || ncol(offset) != K) {
    stop(sprintf(
      "multinomial %s must be a %d x %d matrix",
      what,
      n,
      K
    ))
  }
  if (!all(is.finite(offset))) {
    stop(sprintf("multinomial %s must be finite", what))
  }
  storage.mode(offset) <- "double"
  offset
}

# The creation-time test twin, shaped against the host data object's test rows:
# without them there is nothing for a per-test-row offset to describe, and
# accepting one would leave it silently unread.
validateCategoryTestOffset <- function(offset.test, sampler, K) {
  if (is.null(offset.test)) {
    return(NULL)
  }
  if (is.null(sampler$data@x.test)) {
    stop(
      "multinomial category test offset requires test data: it carries one ",
      "row per test row"
    )
  }
  validateCategoryOffset(
    offset.test,
    nrow(sampler$data@x.test),
    K,
    "category test offset"
  )
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
# offset, when given, is the n x K category offset bartcoreSetCategoryOffset
# installs; offset.test the nTest x K one bartcoreSetCategoryTestOffset
# installs over the host data object's x.test. See those two for the semantics
# and for what they refuse.
bartcoreMultinomialSampler <- function(
  sampler,
  labels,
  K = NULL,
  offset = NULL,
  offset.test = NULL
) {
  labels <- as.integer(labels)
  if (anyNA(labels)) {
    stop("multinomial labels must be integer category codes 0..K-1")
  }
  if (any(labels < 0L)) {
    stop("multinomial labels must be non-negative category codes 0..K-1")
  }
  if (is.null(K)) {
    K <- max(labels) + 1L
  }
  offset <- validateCategoryOffset(offset, length(labels), K)
  offset.test <- validateCategoryTestOffset(offset.test, sampler, K)
  result <- new.env(parent = emptyenv())
  result$ptr <- .Call(
    C_dbarts_bartcore_createMultinomial,
    sampler$control,
    sampler$model,
    sampler$data,
    labels,
    as.integer(K),
    offset,
    offset.test
  )
  # the engine keeps no predictor matrix; track it R-side for the re-quantize
  # surface, as the other dense-predictor wrappers do
  result$x <- rawPredictorMatrix(sampler$data@x)
  # K, so the response-side entries can state the expected shape R-side; the
  # engine reads its own (numReportedLocations) and re-checks
  result$K <- as.integer(K)
  result
}

# The grouped-count analog of bartcoreMultinomialSampler (docs/design/
# multinomial.md): the response is an n x K matrix of nonnegative integer
# counts, column k holding category k's success counts per observation, with
# trials n_i = sum_k counts[i, k] (>= 1). K defaults to the column count. Same
# K-forest softmax engine; the single-trial label path is the special case of a
# one-hot matrix with every trial 1. Validated R-side (safe over fast); the
# engine re-derives the trials and re-checks the invariants.
bartcoreMultinomialCountSampler <- function(
  sampler,
  counts,
  K = NULL,
  offset = NULL,
  offset.test = NULL
) {
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) {
    stop("multinomial counts must be a numeric matrix of non-negative integers")
  }
  if (anyNA(counts)) {
    stop("multinomial counts must not contain missing values")
  }
  if (any(counts < 0)) {
    stop("multinomial counts must be non-negative")
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
  counts <- asCountMatrix(counts)
  offset <- validateCategoryOffset(offset, nrow(counts), K)
  offset.test <- validateCategoryTestOffset(offset.test, sampler, K)
  result <- new.env(parent = emptyenv())
  # counts is column-major (category-major), the layout the combiner reads
  result$ptr <- .Call(
    C_dbarts_bartcore_createMultinomialCounts,
    sampler$control,
    sampler$model,
    sampler$data,
    counts,
    as.integer(K),
    offset,
    offset.test
  )
  # the engine keeps no predictor matrix; track it R-side for the re-quantize
  # surface, as the other dense-predictor wrappers do
  result$x <- rawPredictorMatrix(sampler$data@x)
  # K, so the response-side entries can state the expected shape R-side; the
  # engine reads its own (numReportedLocations) and re-checks
  result$K <- as.integer(K)
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
    stop("multinomial counts must be a numeric matrix of non-negative integers")
  }
  if (anyNA(counts)) {
    stop("multinomial counts must not contain missing values")
  }
  if (any(counts < 0)) {
    stop("multinomial counts must be non-negative")
  }
  if (any(counts != round(counts))) {
    stop("multinomial counts must be whole numbers")
  }
  if (any(rowSums(counts) < 1)) {
    stop("every multinomial count row must have at least one trial (n_i >= 1)")
  }
  counts <- asCountMatrix(counts)
  # counts is column-major (category-major), the layout the combiner reads
  invisible(.Call(C_dbarts_bartcore_setCounts, bcSampler$ptr, counts))
}

# Installs (or clears, at NULL) a multinomial sampler's n x K category offset:
# the latent becomes f_ik + o_ik, so the offset enters the log-sum-exp margins,
# each category's working response and the reported softmax probabilities, and
# never a leaf value. n and K are fixed at creation, as they are for the counts.
#
# This is the response-side counterpart of bartcoreSetCounts, not
# bartcoreSetOffset: the response model's offset is added to every reported
# channel AFTER the K forests are blended, which for a softmax is the wrong side
# of the nonlinearity, and a flat per-observation offset is the softmax's own
# null direction in any case. Only the row-centred part is identified - adding a
# constant to a whole row of the matrix leaves every reported probability
# unchanged - and the entrance leaves the input as given rather than re-centring
# it.
#
# This one shifts the TRAIN latent only. The test rows are other rows, so they
# carry their own offset (bartcoreSetCategoryTestOffset), and predict takes its
# own matrix per call; neither is derived from this one.
bartcoreSetCategoryOffset <- function(bcSampler, offset) {
  # the shape check needs the handle's own K, which only a multinomial handle
  # carries; off one, the entry's capability probe is the refusal, and it names
  # the family situation rather than a shape
  if (!is.null(bcSampler$K)) {
    offset <- validateCategoryOffset(offset, nrow(bcSampler$x), bcSampler$K)
  }
  invisible(.Call(
    C_dbarts_bartcore_setCategoryOffset,
    bcSampler$ptr,
    offset
  ))
}

# Installs (or clears, at NULL) a multinomial sampler's nTest x K category test
# offset: the recorded test channel becomes softmax(f_test + o_test), formed
# where the train blend forms softmax(f + o). The test fits enter no
# likelihood, so this moves the reported test probabilities and nothing else -
# no draw, no working response, no train channel.
#
# Its rows are the CURRENT test rows. Replacing those rows while it is
# installed is refused rather than silently reinterpreted (clear it first);
# out-of-sample predict does not read it at all, taking its own matrix for the
# rows it is given.
bartcoreSetCategoryTestOffset <- function(bcSampler, offset.test) {
  # as for the train offset, the shape check needs the handle's own K; off a
  # multinomial handle the entry's capability probe is the refusal. The ROW
  # count belongs to the sampler's test store, which the handle does not track,
  # so the engine is what pins it - here only K and finiteness are checked
  if (!is.null(bcSampler$K) && !is.null(offset.test)) {
    offset.test <- as.matrix(offset.test)
    offset.test <- validateCategoryOffset(
      offset.test,
      nrow(offset.test),
      bcSampler$K,
      "category test offset"
    )
  }
  invisible(.Call(
    C_dbarts_bartcore_setCategoryTestOffset,
    bcSampler$ptr,
    offset.test
  ))
}

# Forest `forest`'s (0-based) amplitude basis, an n x q numeric matrix; re-forms
# every multiplier and both residuals on the next run. The sole basis-mutation
# route, at any forest and any width: the amplitudes are preserved and remapped,
# and a width-preserving install is the bitwise identity on all of them.
bartcoreSetForestBasis <- function(bcSampler, forest, basis) {
  basis <- as.matrix(basis)
  storage.mode(basis) <- "double"
  invisible(.Call(
    C_dbarts_bartcore_setForestBasis,
    bcSampler$ptr,
    as.integer(forest),
    basis
  ))
}

# A per-forest, per-observation weight s on forest `forest` (0 prognostic, 1
# treatment): a multiplicative precision factor on that forest's own leaf
# conditionals, composing with the case weights and the active-row mask so
# the forest's draws see (w_i * a_i) * m_f^2 * s_i. It does not remove the
# row from occupancy, from the empty-leaf veto, from the combination (the
# row still receives m_f f_f(x_i)), or from the residual sigma degrees of
# freedom; s_i = 0 says only that row i carries no information about forest
# f, and that forest's leaves over such rows stay well-defined prior draws.
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

# Forest `forest`'s (0-based) amplitudes on the combining response, a
# q_forest x n.chains matrix, or - at forest = NULL - the whole vector stacked
# forest-major, sum q_f x n.chains, which is the shape the run's own `glue`
# channel carries. Ragged forest by forest, which is why a forest is named:
# bcf's forest 0 carries a, its forest 1 the pair (b0, b1).
bartcoreForestAmplitudes <- function(bcSampler, forest = NULL) {
  .Call(
    C_dbarts_bartcore_getForestAmplitudes,
    bcSampler$ptr,
    if (is.null(forest)) NULL else as.integer(forest)
  )
}

# A forest's internal-scale function values (0 prognostic, 1 treatment),
# n.observations x n.chains. For a multinomial handle (forest = category,
# 0-based) this is the OFFSET-FREE totalFits: with a category offset
# installed, softmax_k(bartcoreForestFits(bc, k) + offset[, k]) reproduces
# the reported train channel, not softmax_k(bartcoreForestFits(bc, k)) alone.
bartcoreForestFits <- function(bcSampler, forest) {
  .Call(C_dbarts_bartcore_getForestFits, bcSampler$ptr, as.integer(forest))
}

# The combined per-observation location on the RESPONSE scale and without the
# offset, n.observations x n.chains - the quantity the run's train channel
# carries with the offset folded in. Refused on a multinomial handle, whose
# reported channels are per-category softmax probabilities and not one additive
# location; bartcorePredict(bc, x) serves that read instead, and reports the
# saved samples rather than the current state under keepTrees.
bartcoreFitsWithoutOffset <- function(bcSampler) {
  .Call(C_dbarts_bartcore_getFitsWithoutOffset, bcSampler$ptr)
}

# A forest's leaf-prior calibration in RESPONSE units, one row per chain
# (columns prior.scale, prior.sd, prior.mean, k, k.has.hyperprior,
# response.scale, response.shift, then the calibration map's own
# amplitude.prior.variance, amplitude.prior.scale, node.scale.factor,
# node.scale.divisor and basis.row.norm, NaN off the map) with the leaf model
# on a "leaf.model" attribute. Total over forests: a combiner's forests report
# the calibration its own map fixed. Forest is 0-based here, as everywhere on
# this layer; the R-level $getCalibration indexes from 1.
bartcoreForestCalibration <- function(bcSampler, forest) {
  .Call(C_dbarts_bartcore_getCalibration, bcSampler$ptr, as.integer(forest))
}

# Restates a forest's leaf prior on every chain so the forest total's prior sd
# at k = 1 is prior.scale, response units. Refused when a combiner owns the
# calibration; nothing else moves, and a write reproducing what is in force is
# skipped bitwise.
bartcoreSetForestPriorScale <- function(bcSampler, forest, prior.scale) {
  invisible(.Call(
    C_dbarts_bartcore_setCalibration,
    bcSampler$ptr,
    as.integer(forest),
    validateLiveScale(prior.scale, "prior.scale")
  ))
}

# The R-level calibration surface indexes forests from 1, as R indexes; the
# bridge and the flat C entries count from 0, as the engine does.
resolveForestIndex <- function(forest) {
  forest <- as.integer(forest)
  if (length(forest) != 1L || is.na(forest) || forest < 1L) {
    stop("'forest' must be a single positive integer (1 selects the first)")
  }
  forest - 1L
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

bartcoreSetModel <- function(bcSampler, model, data) {
  invisible(.Call(C_dbarts_bartcore_setModel, bcSampler$ptr, model, data))
}

# Drives a low-level bartcore handle (a bcSampler env holding $ptr) directly,
# not a dbartsSampler; cf. bartcoreSamplerRun, the R-level sampler-layer entry.
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

# A per-observation 0/1 mask saying which rows are in the data set for this
# sampler this sweep. An inactive row leaves every sufficient statistic, every
# family-level parameter update and its own latent draw, but keeps its leaf
# occupancy and still receives a fitted value. Absolute and independent of the
# case weights: the family serves w * a in either call order.
#
# NULL clears. An all-ones mask reports success and installs nothing - the
# opposite of setWeights, where all-ones installs - because one channel is
# membership and the other precision. An all-zeros mask is accepted and runs,
# with every forest at its prior. Values other than 0 and 1 refuse the whole
# call and install nothing.
bartcoreSetActiveRows <- function(bcSampler, active) {
  if (!is.null(active)) {
    active <- as.double(active)
  }
  invisible(.Call(C_dbarts_bartcore_setActiveRows, bcSampler$ptr, active))
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

# offset.test takes the shape of the surface: a flat per-row vector for every
# additive family, and for a multinomial handle an nNew x K matrix, one row per
# predicted row, entering the raw fits before the softmax. It is never taken
# from the sampler - these rows are the caller's - so a handle carrying a
# category offset must be given one here, an all-zero matrix being how the
# offset-free surface is asked for.
bartcorePredict <- function(bcSampler, x.test, offset.test = NULL) {
  x.test <- as.matrix(x.test)
  storage.mode(x.test) <- "double"
  if (!is.null(offset.test)) {
    if (!is.null(bcSampler$K)) {
      offset.test <- validateCategoryOffset(
        offset.test,
        nrow(x.test),
        bcSampler$K,
        "predict category offset"
      )
    } else {
      offset.test <- as.double(offset.test)
    }
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
