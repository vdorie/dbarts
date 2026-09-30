#!/usr/bin/env Rscript

# Size of the monotone leaf order in fitted trees, the input to the count
# limit in docs/plans/monotone-exact-birth-death.md (Decision). The exact
# normalizer Z_T = e(P_T) / L! is counted by a DP over the down-sets of the
# leaf order P_T, one component at a time; its work is about (down-sets) x
# (leaves) per component. A fit's trees are read with getTrees(); each tree's
# order is built with the engine's adjacency rule (j < k when j's box ends
# where k's starts on a constrained axis and the boxes overlap on every other
# axis), closed transitively, split into components, and its down-sets are
# counted by memoized recursion.
#
# Per fit it prints the sampler's time per kept sweep, the most leaves any
# kept tree has, and, as median and maximum over kept sweeps: the largest
# whole-tree down-set count (the product over components) in the sweep, with
# the share of sweeps where it exceeds 2^22; the largest component count; and
# the sweep's count work, summed over its trees, of a move whose touched leaf
# is uniform: sum over components of (|C| / L) D(C) |C|, in down-set x leaf
# units. On every fifth kept sweep it also applies every death and 10 random
# births to each tree and prints the largest component a move creates, which
# is what the count sees, the share of trees with a death whose component
# passes 2^24, and how often per sweep such a death is proposed. masks=<file> writes every distinct component a state
# or one of those moves holds for benchmarks/kernels/monotone_count, which
# times the engine's layered count on them.
#
# closure mode checks which measures deaths never raise, over every death of
# random guillotine trees, and that no birth lowers the number of linear
# extensions e (so a birth's Z_T0 / Z_T* is at most L0 + 1). It then prints the
# two-by-two tree whose every death raises the largest component's count (a
# per-component budget of 4 traps it away from the root).
#
# logz mode prints, for a fit, how much -log Z_T its trees carry: the per-tree
# weight the documented prior removes relative to the unnormalized one.
#
# pairs mode writes, on every fifth kept sweep, 10 random births and every
# distinct death of each tree, as the pair of orders a move's ratio needs
# (the components of the finer tree holding the two children, and the
# component of the coarser tree holding their merged leaf), for
# benchmarks/kernels/monotone_ratio.
#
# Usage: Rscript monotone-order-size.R nTrees n p nConstrained [seed] [masks=f]
#        Rscript monotone-order-size.R closure
#        Rscript monotone-order-size.R logz nTrees n p nConstrained [seed]
#        Rscript monotone-order-size.R pairs nTrees n p nConstrained seed file
#   x is uniform on [0, 1]^p; y increases in the first nConstrained axes and
#   has 2 sin(2 pi x_free) (1 + x1) in the first free one; the installed
#   engine runs 1000 burn-in and 200 kept sweeps (logz: 500 and 20).

args <- commandArgs(trailingOnly = TRUE)
budget <- 2^22
guard <- 2^24
# a component whose recursion outgrows this many states is reported as Inf
memoCap <- 2^23

# ---- trees, orders and counts ------------------------------------------------

# a tree is list() for a leaf or list(var, cut, left, right); left is x <= cut
fromPreorder <- function(var, value) {
  i <- 0L
  build <- function() {
    i <<- i + 1L
    if (var[i] == -1L) {
      return(list())
    }
    node <- list(var[i], value[i])
    node[[3L]] <- build()
    node[[4L]] <- build()
    node
  }
  build()
}

# each leaf's box and its path (the list indices that reach it)
leafBoxes <- function(tree, p) {
  out <- list()
  walk <- function(node, lo, hi, path) {
    if (!length(node)) {
      out[[length(out) + 1L]] <<- list(lo = lo, hi = hi, path = path)
      return(invisible())
    }
    v <- node[[1L]]
    hiLeft <- hi
    hiLeft[v] <- node[[2L]]
    walk(node[[3L]], lo, hiLeft, c(path, 3L))
    loRight <- lo
    loRight[v] <- node[[2L]]
    walk(node[[4L]], loRight, hi, c(path, 4L))
  }
  walk(tree, rep(-Inf, p), rep(Inf, p), integer(0L))
  out
}

setNode <- function(tree, path, node) {
  if (!length(path)) {
    return(node)
  }
  tree[[path]] <- node
  tree
}

# R[j, k] when leaf j lies below leaf k, transitively closed; every
# constrained axis increasing
leafOrder <- function(boxes, nc) {
  n <- length(boxes)
  p <- length(boxes[[1L]]$lo)
  rel <- matrix(FALSE, n, n)
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      if (j == k) {
        next
      }
      a <- boxes[[j]]
      b <- boxes[[k]]
      for (i in seq_len(nc)) {
        if (a$hi[i] != b$lo[i]) {
          next
        }
        other <- setdiff(seq_len(p), i)
        if (
          all(pmax(a$lo[other], b$lo[other]) < pmin(a$hi[other], b$hi[other]))
        ) {
          rel[j, k] <- TRUE
        }
      }
    }
  }
  repeat {
    closed <- rel | (rel %*% rel > 0)
    if (all(closed == rel)) {
      return(rel)
    }
    rel <- closed
  }
}

components <- function(rel) {
  adj <- rel | t(rel)
  n <- nrow(adj)
  label <- seq_len(n)
  repeat {
    nextLabel <- vapply(
      seq_len(n),
      function(j) min(label[adj[j, ] | seq_len(n) == j]),
      integer(1L)
    )
    if (all(nextLabel == label)) {
      return(unname(split(seq_len(n), label)))
    }
    label <- nextLabel
  }
}

# down-sets of S: those without x (nothing above x) plus those with x (all of
# its down-closure), each a down-set of what is left
countDownSets <- function(rel) {
  memo <- new.env(hash = TRUE)
  states <- 0L
  count <- function(set) {
    if (!length(set)) {
      return(1)
    }
    key <- paste(set, collapse = ".")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    states <<- states + 1L
    if (states > memoCap) {
      stop("memo cap")
    }
    x <- set[1L]
    above <- set[rel[x, set]]
    below <- set[rel[set, x]]
    result <- count(setdiff(set, c(x, above))) +
      count(setdiff(set, c(x, below)))
    assign(key, result, envir = memo)
    result
  }
  tryCatch(count(seq_len(nrow(rel))), error = function(e) Inf)
}

# linear extensions: e(S) sums e(S - x) over the maximal elements x of S
countExtensions <- function(rel) {
  memo <- new.env(hash = TRUE)
  count <- function(set) {
    if (length(set) <= 1L) {
      return(1)
    }
    key <- paste(set, collapse = ".")
    hit <- memo[[key]]
    if (!is.null(hit)) {
      return(hit)
    }
    top <- set[!vapply(set, function(x) any(rel[x, set]), TRUE)]
    result <- sum(vapply(top, function(x) count(setdiff(set, x)), 0))
    assign(key, result, envir = memo)
    result
  }
  count(seq_len(nrow(rel)))
}

logNormalizer <- function(rel) {
  sum(vapply(
    components(rel),
    function(ix) {
      log(countExtensions(rel[ix, ix, drop = FALSE])) - lfactorial(length(ix))
    },
    0
  ))
}

orderSize <- function(tree, p, nc) {
  rel <- leafOrder(leafBoxes(tree, p), nc)
  parts <- components(rel)
  counts <- vapply(
    parts,
    function(ix) countDownSets(rel[ix, ix, drop = FALSE]),
    numeric(1L)
  )
  sizes <- lengths(parts)
  c(
    leaves = nrow(rel),
    whole = prod(counts),
    largest = max(counts),
    work = sum(sizes / nrow(rel) * counts * sizes)
  )
}

# ---- fits --------------------------------------------------------------------

runFit <- function(nTrees, n, p, nc, seed, nBurn, nKept) {
  suppressPackageStartupMessages(library(dbarts))
  set.seed(seed)
  x <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  f <- 2 * rowSums(x[, seq_len(nc), drop = FALSE])
  if (p > nc) {
    f <- f + 2 * sin(2 * pi * x[, nc + 1L]) * (1 + x[, 1L])
  }
  y <- f + rnorm(n, 0, 0.3)
  control <- dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = nTrees,
    n.samples = nKept,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = seed + 1L
  )
  monotone <- setNames(rep(1L, nc), colnames(x)[seq_len(nc)])
  sampler <- dbarts(x, y, control = control, monotone = monotone)
  invisible(sampler$run(nBurn, 0L))
  seconds <- system.time(invisible(sampler$run(0L, nKept)))[["elapsed"]]
  trees <- sampler$getTrees()
  keys <- unique(trees[c("sample", "tree")])
  list(
    trees = lapply(seq_len(nrow(keys)), function(r) {
      trees[trees$sample == keys$sample[r] & trees$tree == keys$tree[r], ]
    }),
    keys = keys,
    secondsPerSweep = seconds / nKept
  )
}

# components written for the C++ counter, keyed by their relation
maskStore <- new.env(hash = TRUE)
# each element's predecessor mask as two 64-bit hex words (high, low)
maskWords <- function(sub) {
  hexWord <- function(bits) {
    nibbles <- matrix(as.integer(bits), 4L)
    paste(
      rev(sprintf("%x", colSums(nibbles * c(1L, 2L, 4L, 8L)))),
      collapse = ""
    )
  }
  words <- vapply(
    seq_len(nrow(sub)),
    function(k) {
      bits <- c(sub[, k], rep(FALSE, 128L - nrow(sub)))
      paste(hexWord(bits[65:128]), hexWord(bits[1:64]))
    },
    ""
  )
  paste(words, collapse = " ")
}

storeComponent <- function(sub) {
  key <- paste(nrow(sub), paste(which(sub), collapse = ","))
  if (is.null(maskStore[[key]]) && nrow(sub) <= 128L) {
    id <- sprintf("c%d", length(ls(maskStore)) + 1L)
    assign(key, paste(id, nrow(sub), maskWords(sub)), maskStore)
  }
}

# the largest down-set count among the components holding the given leaves
touchedCount <- function(tree, p, nc, paths, store) {
  boxes <- leafBoxes(tree, p)
  rel <- leafOrder(boxes, nc)
  leafPaths <- vapply(boxes, function(b) paste(b$path, collapse = "."), "")
  touched <- which(leafPaths %in% paths)
  best <- 0
  for (ix in components(rel)) {
    if (any(touched %in% ix)) {
      sub <- rel[ix, ix, drop = FALSE]
      if (store) {
        storeComponent(sub)
      }
      best <- max(best, countDownSets(sub))
    }
  }
  best
}

# largest component created by every death and by nBirths random births
moveSizes <- function(tree, p, nc, nBirths, store) {
  if (store) {
    rel <- leafOrder(leafBoxes(tree, p), nc)
    for (ix in components(rel)) {
      storeComponent(rel[ix, ix, drop = FALSE])
    }
  }
  death <- 0
  overGuard <- 0
  nogs <- nogPaths(tree)
  for (path in nogs) {
    if (length(path)) {
      dead <- setNode(tree, path, list())
      created <- touchedCount(dead, p, nc, paste(path, collapse = "."), store)
      death <- max(death, created)
      overGuard <- overGuard + (created > guard) / length(nogs)
    }
  }
  birth <- 0
  boxes <- leafBoxes(tree, p)
  for (b in seq_len(nBirths)) {
    box <- boxes[[sample.int(length(boxes), 1L)]]
    v <- sample.int(p, 1L)
    grid <- 1:99 / 100
    cuts <- grid[grid > box$lo[v] & grid < box$hi[v]]
    if (!length(cuts)) {
      next
    }
    cut <- cuts[sample.int(length(cuts), 1L)]
    grown <- setNode(tree, box$path, list(v, cut, list(), list()))
    children <- vapply(
      c(3L, 4L),
      function(side) paste(c(box$path, side), collapse = "."),
      ""
    )
    birth <- max(birth, touchedCount(grown, p, nc, children, store))
  }
  c(death = death, birth = birth, overGuard = overGuard)
}

fitSizes <- function(nTrees, n, p, nc, seed, masks) {
  fit <- runFit(nTrees, n, p, nc, seed, 1000L, 200L)
  trees <- lapply(fit$trees, function(d) fromPreorder(d$var, d$value))
  sizes <- t(vapply(trees, orderSize, numeric(4L), p = p, nc = nc))
  sample <- fit$keys$sample
  perSweep <- function(column, f) tapply(sizes[, column], sample, f)
  whole <- perSweep("whole", max)
  largest <- perSweep("largest", max)
  moved <- which(sample %% 5L == 0L)
  moves <- vapply(
    trees[moved],
    moveSizes,
    numeric(3L),
    p = p,
    nc = nc,
    nBirths = 10L,
    store = !is.null(masks)
  )
  q <- function(v) {
    sprintf("%.2g / %.2g", stats::median(v), max(v))
  }
  cat(sprintf(
    paste0(
      "trees %d n %d p %d nc %d seed %d: %.2f ms/sweep, leaves max %d; ",
      "per sweep med/max: whole-tree %s (over 2^22 in %.0f%%), largest ",
      "component %s, count work %s; largest component created: by a death ",
      "%.2g, by a birth %.2g (state %.2g); trees with a death past 2^24 ",
      "%.0f%%, one proposed per sweep %.3f\n"
    ),
    nTrees,
    n,
    p,
    nc,
    seed,
    1000 * fit$secondsPerSweep,
    max(sizes[, "leaves"]),
    q(whole),
    100 * mean(whole > budget),
    q(largest),
    q(perSweep("work", sum)),
    max(moves["death", ]),
    max(moves["birth", ]),
    max(sizes[moved, "largest"]),
    100 * mean(moves["overGuard", ] > 0),
    # a death is proposed with probability 1/2, at a uniform nog node
    mean(tapply(moves["overGuard", ], sample[moved], function(f) {
      1 - prod(1 - f / 2)
    }))
  ))
  if (!is.null(masks)) {
    writeLines(unlist(mget(ls(maskStore), maskStore)), masks)
  }
}

# -log Z_T per tree and per constrained split in a fit
fitNormalizers <- function(nTrees, n, p, nc, seed) {
  fit <- runFit(nTrees, n, p, nc, seed, 500L, 20L)
  res <- t(vapply(
    fit$trees,
    function(d) {
      rel <- leafOrder(leafBoxes(fromPreorder(d$var, d$value), p), nc)
      c(
        leaves = nrow(rel),
        splits = sum(d$var >= 1L & d$var <= nc),
        logZ = logNormalizer(rel)
      )
    },
    numeric(3L)
  ))
  cat(sprintf(
    paste0(
      "trees %d n %d p %d nc %d seed %d: mean leaves %.2f, trees with a ",
      "constrained split %.0f%%, -log Z_T per tree %.3f, per constrained ",
      "split %.3f\n"
    ),
    nTrees,
    n,
    p,
    nc,
    seed,
    mean(res[, "leaves"]),
    100 * mean(res[, "splits"] > 0),
    -mean(res[, "logZ"]),
    -sum(res[, "logZ"]) / sum(res[, "splits"])
  ))
}

# ---- move pairs for the ratio kernel ----------------------------------------

# one move between tree (T0, the leaf at path) and grown (T*, its children at
# path 3 and 4): the components of T* holding the children, with the lower
# child c1 and the other c2, and the component of T0 holding the leaf
pairLine <- function(id, tree, grown, path, p, nc) {
  indexOf <- function(boxes, target) {
    which(vapply(boxes, function(b) identical(b$path, target), TRUE))
  }
  boxes <- leafBoxes(grown, p)
  rel <- leafOrder(boxes, nc)
  children <- c(indexOf(boxes, c(path, 3L)), indexOf(boxes, c(path, 4L)))
  parts <- components(rel)
  held <- sort(unlist(parts[vapply(
    parts,
    function(ix) any(children %in% ix),
    TRUE
  )]))
  boxes0 <- leafBoxes(tree, p)
  rel0 <- leafOrder(boxes0, nc)
  leaf <- indexOf(boxes0, path)
  parts0 <- components(rel0)
  c0 <- parts0[[which(vapply(parts0, function(ix) leaf %in% ix, TRUE))]]
  if (length(held) > 128L) {
    return(NULL)
  }
  sprintf(
    "pair %s %d %d %d %d %s | %d %s",
    id,
    (if (length(path)) grown[[path]] else grown)[[1L]],
    length(held),
    match(children[1L], held) - 1L,
    match(children[2L], held) - 1L,
    maskWords(rel[held, held, drop = FALSE]),
    length(c0),
    maskWords(rel0[c0, c0, drop = FALSE])
  )
}

# on every fifth kept sweep, 10 random births of each tree and all of its
# deaths (each distinct death once), for benchmarks/kernels/monotone_ratio
fitPairs <- function(nTrees, n, p, nc, seed, file) {
  fit <- runFit(nTrees, n, p, nc, seed, 1000L, 200L)
  con <- file(file, "w")
  on.exit(close(con))
  seen <- new.env(hash = TRUE)
  for (r in seq_along(fit$trees)) {
    if (fit$keys$sample[r] %% 5L) {
      next
    }
    tree <- fromPreorder(fit$trees[[r]]$var, fit$trees[[r]]$value)
    tag <- sprintf("s%d_t%d", fit$keys$sample[r], fit$keys$tree[r])
    boxes <- leafBoxes(tree, p)
    for (b in 1:10) {
      box <- boxes[[sample.int(length(boxes), 1L)]]
      v <- sample.int(p, 1L)
      grid <- 1:99 / 100
      cuts <- grid[grid > box$lo[v] & grid < box$hi[v]]
      if (!length(cuts)) {
        next
      }
      cut <- cuts[sample.int(length(cuts), 1L)]
      grown <- setNode(tree, box$path, list(v, cut, list(), list()))
      line <- pairLine(sprintf("B_%s_%d", tag, b), tree, grown, box$path, p, nc)
      writeLines(line, con)
    }
    for (path in nogPaths(tree)) {
      if (!length(path)) {
        next
      }
      dead <- setNode(tree, path, list())
      line <- pairLine(
        sprintf("D_%s_%s", tag, paste(path, collapse = ".")),
        dead,
        tree,
        path,
        p,
        nc
      )
      key <- sub("^pair \\S+ ", "", line)
      if (length(line) && is.null(seen[[key]])) {
        assign(key, TRUE, seen)
        writeLines(line, con)
      }
    }
  }
  invisible()
}

# ---- closure under deaths ----------------------------------------------------

randomTree <- function(nLeaves, p) {
  tree <- list()
  paths <- list(integer(0L))
  boxes <- list(list(lo = rep(0, p), hi = rep(1, p)))
  while (length(paths) < nLeaves) {
    j <- sample.int(length(paths), 1L)
    box <- boxes[[j]]
    v <- sample.int(p, 1L)
    cuts <- (1:15 / 16)[1:15 / 16 > box$lo[v] & 1:15 / 16 < box$hi[v]]
    if (!length(cuts)) {
      next
    }
    cut <- cuts[sample.int(length(cuts), 1L)]
    node <- list(v, cut, list(), list())
    if (length(paths[[j]])) {
      tree[[paths[[j]]]] <- node
    } else {
      tree <- node
    }
    left <- box
    left$hi[v] <- cut
    right <- box
    right$lo[v] <- cut
    paths <- c(paths[-j], list(c(paths[[j]], 3L), c(paths[[j]], 4L)))
    boxes <- c(boxes[-j], list(left, right))
  }
  tree
}

# paths to the internal nodes whose children are both leaves
nogPaths <- function(node, path = integer(0L)) {
  if (!length(node)) {
    return(list())
  }
  if (!length(node[[3L]]) && !length(node[[4L]])) {
    return(list(path))
  }
  c(nogPaths(node[[3L]], c(path, 3L)), nogPaths(node[[4L]], c(path, 4L)))
}

closureCheck <- function(nTrees = 300L, seed = 5L) {
  set.seed(seed)
  raised <- c(whole = 0L, largest = 0L, work = 0L)
  deaths <- 0L
  fewerExtensions <- 0L
  worstRatio <- 0
  for (r in seq_len(nTrees)) {
    p <- sample(2:3, 1L)
    nc <- sample.int(p, 1L)
    tree <- randomTree(sample(6:16, 1L), p)
    before <- orderSize(tree, p, nc)
    relBefore <- leafOrder(leafBoxes(tree, p), nc)
    logEBefore <- logNormalizer(relBefore) + lfactorial(nrow(relBefore))
    for (path in nogPaths(tree)) {
      dead <- tree
      if (length(path)) {
        dead[[path]] <- list()
      } else {
        dead <- list()
      }
      after <- orderSize(dead, p, nc)
      deaths <- deaths + 1L
      # the death's reverse is a birth from dead to tree: e(tree) >= e(dead)
      relAfter <- leafOrder(leafBoxes(dead, p), nc)
      logEAfter <- logNormalizer(relAfter) + lfactorial(nrow(relAfter))
      fewerExtensions <- fewerExtensions + (logEBefore < logEAfter - 1e-9)
      worstRatio <- max(worstRatio, exp(logEAfter - logEBefore))
      for (m in names(raised)) {
        raised[[m]] <- raised[[m]] + (after[[m]] > before[[m]] * (1 + 1e-12))
      }
    }
  }
  cat(sprintf(
    "%d deaths of %d random trees raise: whole-tree %d, largest component %d, move work %d\n",
    deaths,
    nTrees,
    raised[["whole"]],
    raised[["largest"]],
    raised[["work"]]
  ))
  cat(sprintf(
    "births that lower e: %d of %d; largest e(T0) / e(T*) %.3g\n",
    fewerExtensions,
    deaths,
    worstRatio
  ))

  # x1 constrained, x2 free: x1 split, then x2 on both sides; A < R1, B < R2
  trap <- list(
    1L,
    0.5,
    list(2L, 0.5, list(), list()),
    list(2L, 0.5, list(), list())
  )
  show <- function(label, tree) {
    s <- orderSize(tree, 2L, 1L)
    cat(sprintf(
      "%s: whole-tree %g, largest component %g\n",
      label,
      s[["whole"]],
      s[["largest"]]
    ))
  }
  show("two-by-two tree", trap)
  lower <- trap
  lower[[3L]] <- list()
  show("  after the lower death", lower)
  upper <- trap
  upper[[4L]] <- list()
  show("  after the upper death", upper)
}

masks <- sub("^masks=", "", grep("^masks=", args, value = TRUE))
args <- grep("^masks=", args, value = TRUE, invert = TRUE)
if (identical(args[1L], "pairs")) {
  a <- as.integer(args[2:6])
  fitPairs(a[1L], a[2L], a[3L], a[4L], a[5L], args[7L])
} else if (identical(args[1L], "closure")) {
  closureCheck()
} else if (identical(args[1L], "logz")) {
  a <- as.integer(args[-1L])
  fitNormalizers(a[1L], a[2L], a[3L], a[4L], if (length(a) > 4L) a[5L] else 1L)
} else {
  a <- as.integer(args)
  fitSizes(
    a[1L],
    a[2L],
    a[3L],
    a[4L],
    if (length(a) > 4L) a[5L] else 1L,
    if (length(masks)) masks else NULL
  )
}
