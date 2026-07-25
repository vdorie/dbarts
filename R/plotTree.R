## Pad the categorical direction strings a getTrees data.frame carries to the
## declared level count. C reports every categorical rule already decoded over
## its observed categories (character k sends the variable's k'th level down
## the "L"eft or "R"ight child; the rule's value is NA); this right-pads with
## "L" to the level count, exactly what unobserved trailing levels decode to
## under the canonical gauge. Ordinal rules and leaves stay NA.
decodeCategoricalSplits <- function(trees, x, varTypes) {
  factorLevels <- attr(x, "factor.levels")
  numLevels <- integer(ncol(x))
  for (j in which(varTypes == CATEGORICAL_VARIABLE)) {
    numLevels[j] <- if (!is.null(factorLevels) && !is.null(factorLevels[[j]])) {
      length(factorLevels[[j]])
    } else {
      as.integer(max(x[, j])) + 1L
    } # hand-typed matrix: the levels are the codes
  }

  directions <- if (!is.null(trees$directions)) {
    trees$directions
  } else {
    rep.int(NA_character_, nrow(trees))
  }
  isCategoricalRule <- trees$var > 0L
  isCategoricalRule[isCategoricalRule] <-
    varTypes[trees$var[isCategoricalRule]] == CATEGORICAL_VARIABLE
  for (i in which(isCategoricalRule)) {
    padding <- numLevels[trees$var[i]] - nchar(directions[i])
    if (padding > 0L) {
      directions[i] <- paste0(
        directions[i],
        paste(rep.int("L", padding), collapse = "")
      )
    }
  }
  trees$directions <- directions
  trees
}

getTreeDepthAndSize <- function(node) {
  if (node$var[1L] == -1) {
    return(c(depth = 1, size = 1))
  }

  left <- getTreeDepthAndSize(node[-1L, , drop = FALSE])
  right <- getTreeDepthAndSize(node[
    seq.int(2L + left[["size"]], nrow(node)),
    ,
    drop = FALSE
  ])

  c(
    depth = 1 + max(left[["depth"]], right[["depth"]]),
    size = 1 + left[["size"]] + right[["size"]]
  )
}

fillPlotCoordinatesForNode <- function(node, maxDepth, currDepth, index) {
  if (node$var[1L] == -1) {
    node$y[1L] <- 1 # maxDepth
    node$x[1L] <- index
    node$index[1L] <- index + 1L
    return(node)
  }

  leftSize <- getTreeDepthAndSize(node[-1L, , drop = FALSE])

  leftNodes <- seq.int(2L, 1L + leftSize[["size"]])
  rightNodes <- seq.int(2L + leftSize[["size"]], nrow(node))

  left <- node[leftNodes, , drop = FALSE]
  left <- fillPlotCoordinatesForNode(left, maxDepth, currDepth + 1L, index)
  index <- left$index[1L]
  left$index[1L] <- NA_integer_

  right <- node[rightNodes, , drop = FALSE]
  right <- fillPlotCoordinatesForNode(right, maxDepth, currDepth + 1L, index)
  node$index[1L] <- right$index[1L]
  right$index[1L] <- NA_integer_

  node$y[1L] <- maxDepth - currDepth + 1L
  node$x[1L] <- (left$x[1L] + right$x[1L]) / 2

  node[leftNodes, ] <- left
  node[rightNodes, ] <- right
  node
}

plotNode <- function(node, sampler, plotPars) {
  nodeHeight <- plotPars[["nodeHeight"]]
  nodeWidth <- plotPars[["nodeWidth"]]
  nodeGap <- plotPars[["nodeGap"]]
  if (node$var[1L] != -1) {
    directions <- if (!is.null(node$directions)) {
      node$directions[1L]
    } else {
      NA_character_
    }
    varName <- if (!is.null(colnames(sampler$data@x))) {
      colnames(sampler$data@x)[node$var[1L]]
    } else {
      NULL
    }
    if (!is.na(directions)) {
      # a categorical rule; like the ordinal '<=', the label states the
      # condition for going left
      factorLevels <- attr(sampler$data@x, "factor.levels")
      levelNames <-
        if (!is.null(factorLevels) && !is.null(factorLevels[[node$var[1L]]])) {
          factorLevels[[node$var[1L]]]
        } else {
          as.character(seq_len(nchar(directions)) - 1L)
        }
      expr1 <- paste0(
        if (!is.null(varName)) varName else paste0("x[", node$var[1L], "]"),
        " in {",
        paste(
          levelNames[strsplit(directions, "")[[1L]] == "L"],
          collapse = ","
        ),
        "}"
      )
    } else {
      expr1 <- expression(a <= b)
      if (!is.null(varName)) {
        expr1[[1]][[2]] <- varName
      } else {
        expr1[[1]][[2]] <- quote(x[a])
        expr1[[1]][[2]][[3]] <- node$var[1L]
      }
      expr1[[1]][[3]] <- signif(node$value[1L], 3)
    }
    naRoute <- if (!is.null(node$missing)) node$missing[1L] else NA_character_
    if (!is.na(naRoute)) {
      # the rule also routes missing values; the label becomes plain text
      if (is.expression(expr1)) {
        expr1 <- paste0(
          if (!is.null(varName)) varName else paste0("x[", node$var[1L], "]"),
          " <= ",
          signif(node$value[1L], 3)
        )
      }
      expr1 <- paste0(expr1, ", NA->", naRoute)
    }
  } else {
    slopeColumns <- startsWith(names(node), "beta.")
    if (any(slopeColumns)) {
      # a linear leaf: label with its fitted plane's coefficients
      slopes <- unlist(node[1L, slopeColumns])
      covariateNames <- sub("^beta\\.", "", names(node)[slopeColumns])
      expr1 <- paste0(
        signif(node$value[1L], 3),
        paste0(
          ifelse(slopes < 0, " - ", " + "),
          signif(abs(slopes), 3),
          " ",
          covariateNames,
          collapse = ""
        )
      )
    } else if (is.na(node$value[1L])) {
      # a function-valued (gp) leaf: the fit varies within the node and
      # rides prediction only
      expr1 <- "gp"
    } else {
      expr1 <- expression(mu == b)
      expr1[[1]][[3]] <- signif(node$value[1L], 3)
    }
  }
  expr2 <- expression(n == a)
  expr2[[1]][[3]] <- node$n[1L]

  plotExpr <- expression(atop(a, b))
  plotExpr[[1]][[2]] <- expr1[[1]]
  plotExpr[[1]][[3]] <- expr2[[1]]
  y <- node$y[1L] * nodeHeight - nodeHeight / 2 + (node$y[1L] - 1) * nodeGap
  x <- node$x[1L] * nodeWidth - nodeWidth / 2

  cex <- par("cex")
  verticalOffset <- graphics::strheight("\n", cex = cex) / 2
  graphics::text(x, y + verticalOffset, expr1, adj = c(0.5, 0.5), cex = cex)
  graphics::text(
    x,
    y - verticalOffset,
    expr2,
    adj = c(0.5, 0.5),
    cex = cex
  )

  if (node$var[1L] != -1) {
    leftSize <- getTreeDepthAndSize(node[-1L, , drop = FALSE])

    leftNodes <- seq.int(2L, 1L + leftSize[["size"]])
    rightNodes <- seq.int(2L + leftSize[["size"]], nrow(node))

    left <- node[leftNodes, , drop = FALSE]
    right <- node[rightNodes, , drop = FALSE]

    plotNode(left, sampler, plotPars)
    plotNode(right, sampler, plotPars)

    y.l <- left$y[1L] * nodeHeight - nodeHeight / 2 + (left$y[1L] - 1) * nodeGap
    x.l <- left$x[1L] * nodeWidth - nodeWidth / 2
    y.r <- right$y[1L] *
      nodeHeight -
      nodeHeight / 2 +
      (right$y[1L] - 1) * nodeGap
    x.r <- right$x[1L] * nodeWidth - nodeWidth / 2

    skippedSpace <- (node$y[1L] - left$y[1L] - 1) * (nodeHeight + nodeGap)

    y.m <- (y + y.l) / 2
    x.m <- (x + x.l) / 2
    theta <- atan2(y - y.m, x - x.m)
    segmentLength <- (nodeGap + skippedSpace) / 2
    y.1 <- segmentLength * sin(theta) + y.m
    x.1 <- segmentLength * cos(theta) + x.m
    y.2 <- y.m - segmentLength * sin(theta)
    x.2 <- x.m - segmentLength * cos(theta)
    lines(c(x.1, x.2), c(y.1, y.2))

    skippedSpace <- (node$y[1L] - right$y[1L] - 1) * (nodeHeight + nodeGap)

    y.m <- (y.r + y) / 2
    x.m <- (x.r + x) / 2
    theta <- atan2(y.m - y, x.m - x)
    segmentLength <- (nodeGap + skippedSpace) / 2
    y.1 <- segmentLength * sin(theta) + y.m
    x.1 <- segmentLength * cos(theta) + x.m
    y.2 <- y.m - segmentLength * sin(theta)
    x.2 <- x.m - segmentLength * cos(theta)
    lines(c(x.1, x.2), c(y.1, y.2))
  }
}
