rejectionSample <- function(target, dgenerator, rgenerator, constant, boundary, log = TRUE, maxIter = 100L)
{
  useLog <- log; rm(log)
  numIters <- 0
  if (useLog) {
    while (TRUE) {
      u <- -rexp(1)
      x <- rgenerator()
      numIters <- numIters + 1
      if (x <= boundary[1] || x >= boundary[2]) next
      if (u < target(x) - constant - dgenerator(x)) return(x)
      if (numIters == maxIter) stop("unable to obtain rejection sample after ", maxIter)
    }
  } else {
    while (TRUE) {
      u <- runif(1)
      x <- rgenerator()
      numIters <- numIters + 1
      if (x <= boundary[1] || x >= boundary[2]) next
      if (u < (target(x) / (constant * dgenerator(x)))) return(x)
      if (numIters == maxIter) stop("unable to obtain rejection sample after ", maxIter)
    }
  }
}

sliceSample <- function(target, start, numSamples = 100L, width = NA, maxIter = 100L, boundary = c(-Inf, Inf), log = TRUE) {
  useLog <- log; rm(log)
  
  findMode <- function(target, start, boundary) {
    optimResult <- tryCatch(optim(start, target, method = "L-BFGS-B", lower = boundary[1L], upper = boundary[2L], hessian = TRUE, control = list(fnscale = -1)),
                            error = function(e) e)
    if (inherits(optimResult, "error")) ## ||
        ##(is.finite(boundary[1]) && abs(optimResult$par - boundary[1]) < 1e-6) ||
        ##(is.finite(boundary[2]) && abs(optimResult$par - boundary[2]) < 1e-6))
    {
      ## if optim fails, do own gradient ascent
      delta <- 1e-6
      while (start - 2 * delta <= boundary[1L] && start + 2 * delta >= boundary[2L]) delta <- delta / 2
      if (start - delta <= boundary[1L]) {
        lh <- target(start); mh <- target(start + delta); rh <- target(start + 2 * delta)
        deriv <- -(3 * lh - 4 * mh + rh) / (2 * delta)
        hess <- (lh - 2 * mh + rh) / (delta^2)
      } else if (start + delta >= boundary[2L]) {
        lh <- target(start - 2 * delta); mh <- target(start - delta); rh <- target(start)
        deriv <-  (3 * rh - 4 * mh + lh) / (2 * delta)
        hess <- (rh - 2 * mh + lh) / (delta^2)
      } else {
        rh <- target(start + delta); mh <- target(start); lh <- target(start - delta)
        deriv <- (rh - lh) / (2 * delta)
        hess <- (rh - 2 * mh + lh) / (delta^2)
      }
      step <- abs(deriv / hess) / 5
      
      lh <- start - step; mh <- start; rh <- start + step
      if (lh <= boundary[1L]) lh <- if (is.finite(boundary[1L])) boundary[1L] + delta else start - delta
      if (rh >= boundary[2L]) rh <- if (is.finite(boundary[2L])) boundary[2L] - delta else start + delta
      lf <- target(lh); mf <- target(mh); rf <- target(rh)
      
      ## keep going until we've got a middle that is higher than one side or the other
      while ((lf < mf && mf < rf) || (lf > mf && mf > rf)) {
        if (lf >= mf) {
          mh <- lh; mf <- lf
          lh <- lh - step
          if (lh <= boundary[1L]) lh <- if (is.finite(boundary[1L])) boundary[1L] + delta else mh - delta
          lf <- target(lh)
        } else {
          mh <- rh; mf <- rf
          rh <- rh + step
          if (rh >= boundary[2L]) rh <- if (is.finite(boundary[2L])) boundary[2L] - delta else mh + delta
          rf <- target(rh)
        }
      }
            
      optimResult <- tryCatch(optim(mh, target, method ="L-BFGS-B", lower = boundary[1L], upper = boundary[2L],
                                    hessian = TRUE, control = list(fnscale = -1)),
                              error = function(e) e)
    }
    optimResult
  }
  
  getInterval <- function(f, x, width, height, boundary) {
    r <- runif(1L)
    x.l <- x - r * width
    x.r <- x + (1 - r) * width

    if (is.finite(boundary[1L])) {
      while (x.l > boundary[1L] && f(x.l) > height) x.l <- x.l - width
      if (x.l < boundary[1L]) x.l <- boundary[1L]
    } else {
      while(f(x.l) > height) x.l <-  x.l - width
    }
    if (is.finite(boundary[2L])) {
      while (x.r < boundary[2L] && f(x.r) > height) x.r <- x.r + width
      if (x.r > boundary[2L]) x.r <- boundary[2L]
    } else {
      while(f(x.r) > height) x.r <-  x.r + width
    }
    
    return(c(x.l, x.r))
  }
  shrinkInterval <- function(x, x.p, int) {
    if (x.p > x) int[2L] <- x.p else int[1L] <- x.p
    return(int)
  }
  
  f <- target
  if (is.na(width)) {
    optimResult <- findMode(target, start, boundary)
    
    if (!inherits(optimResult, "error")) {
      if (useLog == TRUE) {
        normalizingConstant <- NULL ## for R CMD check
        evalEnv <- list2env(list(target = target, normalizingConstant = optimResult$value))
        f <- function(x) exp(target(x) - normalizingConstant)
        environment(f) <- evalEnv
        optimResult$value <- 1
        ## optimResult$hessian is this is theoretically unchanged by the transformation, since f'(x_0) = 0 && (h(x_0) - normConst) = 0, however it can be inaccurate numerically
        optimResult$hessian <- optimHess(optimResult$par, f, control = list(fnscale = -1))
      }
      # width is derived from a normal approximation at the mode based on the equality of second derivatives
      # going out two standard deviations
      #width <- 2 * abs(sqrt(2 * pi) * optimResult$value / optimResult$hessian[1L])^(1/3)
      width <- 2 * abs(optimResult$hessian[1L] * sqrt(2 * pi))^(-1/3)
      if (is.nan(width) || is.infinite(width)) width <- 1000
    } else {
      width <- 1000
    }
  }
  
  result <- rep(NA_real_, numSamples)
  x <- start
  f.x <- f(x)
  if (f.x <= 1e-2) {
    ## if the starting point has a really low density, we grab a different one using
    ## rejection sampling
    mu <- NULL; sigma <- NULL ## for R CMD check
    # us normal approximation with a slightly inflated standard deviation
    evalEnv <- list2env(list(mu = start, sigma = 1.15 * width / 2))
    r <- function() rnorm(1, mu, sigma)
    d <- function(x) dnorm(x, mu, sigma, log = TRUE)
    environment(r) <- evalEnv; environment(d) <- evalEnv
    if (exists("optimResult") && useLog == TRUE) {
      if (inherits(optimResult, "error")) stop("slice sampler failed: unable to determine initial curvature")
      evalEnv$mu <- optimResult$par
      # special case for variance components
      if (is.finite(boundary[1L]) && optimResult$par - boundary[1L] < evalEnv$sigma)
        evalEnv$sigma <- optimResult$par - boundary[1L]
      if (is.finite(boundary[2L]) && boundary[2L] - optimResult$par < evalEnv$sigma)
        evalEnv$sigma <- boundary[2L] - optimResult$par
      c <- target(optimResult$par) - d(optimResult$par)
    } else {
      stop("rejection start for case without optimization and/or not on log scale not yet implemented")
    }
      
    tryResult <- tryCatch(x <- rejectionSample(target, d, r, c, boundary, maxIter = maxIter), error = function(e) e)
    if (inherits(tryResult, "error")) {
      warning("rejection sample failed after ", maxIter, " iterations; dominating function may require hand-tuning")
      x <- start
    } else {
      f.x <- f(x)
    }
  }
  for (i in seq_len(numSamples)) {
    u.p <- runif(1L, 0, f.x)
    int <- getInterval(f, x, width, u.p, boundary)
    for (j in seq_len(maxIter)) {
      x.p <- runif(1, int[1L], int[2L])
      f.x <- f(x.p)
      if (is.nan(f.x) || is.infinite(f.x)) { stop("slice sampler failed: likely due to underflow") }
      if (f.x > u.p) break
      int <- shrinkInterval(x, x.p, int)
    }
    if (j == maxIter) { stop("slice sampler failed: maxIter reached") }
    x <- x.p
    result[i] <- x
  }
  result
}

