buildTree <- function(treeChars)
{
  if (treeChars[1] == ".") return(list(remainder = treeChars[-1]))
  
  splitVar <- as.integer(treeChars[1]) + 1L
  splitIndex <- as.integer(treeChars[2]) + 1L
  
  leftChild <- buildTree(treeChars[-c(1, 2)])
  rightChild <- buildTree(leftChild$remainder)
  leftChild$remainder <- NULL
  remainder <- rightChild$remainder
  rightChild$remainder <- NULL
  
  result <- namedList(splitVar, splitIndex, leftChild, rightChild, remainder)
  leftChild$parent <- result
  rightChild$parent <- result
  
  result
}

fillObservationsForNode <- function(node, sampler, cutPoints)
{
  if (!is.null(node$leftChild)) {
    goesLeft <- sampler$data@x[node$indices, node$splitVar] <= cutPoints[[node$splitVar]][node$splitIndex]
    node$leftChild$indices  <- node$indices[goesLeft]
    node$rightChild$indices <- node$indices[!goesLeft]
    
    node$leftChild  <- fillObservationsForNode(node$leftChild,  sampler, cutPoints)
    node$rightChild <- fillObservationsForNode(node$rightChild, sampler, cutPoints)
  }
  node
}

fillPlotInfoForNode <- function(node, sampler, treeFits)
{
  if (!is.null(node$leftChild)) {
    node$leftChild  <- fillPlotInfoForNode(node$leftChild,  sampler, treeFits)
    node$rightChild <- fillPlotInfoForNode(node$rightChild, sampler, treeFits)
  } else {
    node$mu <- treeFits[node$indices[1]]
  }
  node
}

getNumEndNodes <- function(node)
{
  if (!is.null(node$leftChild))
    return(getNumEndNodes(node$leftChild) + getNumEndNodes(node$rightChild))
  1
}

getMaxDepth <- function(node)
{
  if (!is.null(node$leftChild))
    return(1 + max(getMaxDepth(node$leftChild), getMaxDepth(node$rightChild)))
  1
}

fillPlotCoordinatesForNode <- function(node, maxDepth, currDepth, index)
{
  if (is.null(node$leftChild)) {
    node$y <- 1L # maxDepth
    node$x <- index
    node$index <- index + 1L
    return(node)
  }
  
  node$leftChild <- fillPlotCoordinatesForNode(node$leftChild, maxDepth, currDepth + 1L, index)
  index <- node$leftChild$index
  node$leftChild$index <- NULL
  node$rightChild <- fillPlotCoordinatesForNode(node$rightChild, maxDepth, currDepth + 1L, index)
  node$index <- node$rightChild$index
  node$rightChild$index <- NULL
  
  node$y <- maxDepth - currDepth + 1L
  node$x <- (node$leftChild$x + node$rightChild$x) / 2
  node
}

plotNode <- function(node, sampler, cutPoints, plotPars)
{
  compress <- FALSE
  if (!is.null(node$leftChild)) {
    expr1 <- expression(a <= b)
    if (!is.null(colnames(sampler$data@x))) {
      expr1[[1]][[2]] <- colnames(sampler$data@x)[node$splitVar]
    } else {
      expr1[[1]][[2]] <- quote(x[a])
      expr1[[1]][[2]][[3]] <- node$splitVar
    }
    expr1[[1]][[3]] <- signif(cutPoints[[node$splitVar]][node$splitIndex], 2)
  } else {
    expr1 <- expression(mu == b)
    expr1[[1]][[3]] <- signif(node$mu, 2)
  }
  if (length(node$indices) <= 4) {
    expr2 <- expression(i == group("(", list(a), ")"))
    for (i in seq_along(node$indices))
      expr2[[1]][[3]][[3]][[1 + i]] <- node$indices[i]
    compress <- TRUE
  } else {
    expr2 <- expression(n == a)
    expr2[[1]][[3]] <- length(node$indices)
  }
  plotExpr <- expression(atop(a, b))
  plotExpr[[1]][[2]] <- expr1[[1]]
  plotExpr[[1]][[3]] <- expr2[[1]]
  y <- node$y * plotPars$nodeHeight - plotPars$nodeHeight / 2 + (node$y - 1) * plotPars$nodeGap
  x <- node$x * plotPars$nodeWidth - plotPars$nodeWidth / 2
  
  cex <- par("cex")
  verticalOffset <- strheight("\n", cex = cex) / 2
  text(x, y + verticalOffset, expr1, adj = c(0.5, 0.5), cex = cex)
  text(x, y - verticalOffset, expr2, adj = c(0.5, 0.5), cex = cex * if (compress) 0.8 else 1)
  
  if (!is.null(node$leftChild)) {
    plotNode(node$leftChild, sampler, cutPoints, plotPars)
    plotNode(node$rightChild, sampler, cutPoints, plotPars)
    
    y.l <- node$leftChild$y * plotPars$nodeHeight - plotPars$nodeHeight / 2 + (node$leftChild$y - 1) * plotPars$nodeGap
    x.l <- node$leftChild$x * plotPars$nodeWidth - plotPars$nodeWidth / 2
    y.r <- node$rightChild$y * plotPars$nodeHeight - plotPars$nodeHeight / 2 + (node$rightChild$y - 1) * plotPars$nodeGap
    x.r <- node$rightChild$x * plotPars$nodeWidth - plotPars$nodeWidth / 2
    
    skippedSpace <- (node$y - node$leftChild$y - 1) * (plotPars$nodeHeight + plotPars$nodeGap)
    
    y.m <- (y + y.l) / 2
    x.m <- (x + x.l) / 2
    theta <- atan2(y - y.m, x - x.m)
    segmentLength <- (plotPars$nodeGap + skippedSpace) / 2
    y.1 <- segmentLength * sin(theta) + y.m
    x.1 <- segmentLength * cos(theta) + x.m
    y.2 <- y.m - segmentLength * sin(theta)
    x.2 <- x.m - segmentLength * cos(theta)
    lines(c(x.1, x.2), c(y.1, y.2))
    
    skippedSpace <- (node$y - node$rightChild$y - 1) * (plotPars$nodeHeight + plotPars$nodeGap)
    
    y.m <- (y.r + y) / 2
    x.m <- (x.r + x) / 2
    theta <- atan2(y.m - y, x.m - x)
    segmentLength <- (plotPars$nodeGap + skippedSpace) / 2
    y.1 <- segmentLength * sin(theta) + y.m
    x.1 <- segmentLength * cos(theta) + x.m
    y.2 <- y.m - segmentLength * sin(theta)
    x.2 <- x.m - segmentLength * cos(theta)
    lines(c(x.1, x.2), c(y.1, y.2))
  }
}

createCutPoints <- function(sampler)
{
  if (sampler$control@useQuantiles) {
    cutter <- function(j) {
      uniqueElements <- unique(sampler$data@x[,j])
      numUnique <- length(uniqueElements)
      if (numUnique <= sampler$data@n.cuts[j] + 1L) {
        numCuts <- numUnique - 1L
        step <- 1L
        offset <- 0L
      } else {
        numCuts <- numUnique
        step <- numCuts %/% numUnique
        offset <- step %/% 2
      }
      indices <- sapply(seq.int(numCuts) * step + offset, function(x) min(x, numUnique - 1L))
      sortedElements <- sort(uniqueElements)
      (sortedElements[indices] + sortedElements[indices + 1L]) / 2
    }
  } else {
    cutter <- function(j) {
      m <- min(sampler$data@x[,j]); M <- max(sampler$data@x[,j])
      seq(m, M, length.out = sampler$data@n.cuts[j])
    }
  }
  return(lapply(seq_len(ncol(sampler$data@x)), cutter))
}
