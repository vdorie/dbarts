# test that rbart compares favorably to lmer for nonlinear models
if (
  length(find.package("lme4", quiet = TRUE)) > 0L
  && tryCatch(
    is.environment(lme4 <- asNamespace("lme4")), error = function(e) FALSE
  )
) {
  f <- function(x) {
      10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
  }
  
  set.seed(99L)
  sigma <- 1.0
  n     <- 500
  
  x  <- matrix(runif(n * 10), n, 10)
  Ey <- f(x)
  y  <- rnorm(n, Ey, sigma)
  
  n.g <- 15
  g <- sample(n.g, length(y), replace = TRUE)
  sigma.b <- 1.5
  b <- rnorm(n.g, 0, sigma.b)
  
  y <- y + b[g]
  
  df <- as.data.frame(x)
  colnames(df) <- paste0("x_", seq_len(ncol(x)))
  df$y <- y
  df$g <- g
  
  
  rbartFit <- dbarts::rbart_vi(
    y ~ . - g, df, group.by = g,
    n.samples = 200L, n.burn = 100L, n.thin = 2L, n.chains = 2L,
    n.trees = 50L, n.threads = 1L,
    verbose = FALSE
  )
  ranef.rbart <- rbartFit$ranef.mean
  
  lmerFit <- suppressWarnings(lme4$lmer(y ~ . - g + (1 | g), df))
  ranef.lmer <- lme4$ranef.merMod(lmerFit)[[1L]][[1L]]
  
  b_rmse.rbart <- sqrt(mean((b - ranef.rbart)^2))
  b_rmse.lmer <- sqrt(mean((b - ranef.lmer)^2))
  expect_true(b_rmse.rbart < b_rmse.lmer)
  
  
  rho <- 0.4
  p.y <- pnorm((Ey - mean(Ey)) / sd(Ey) + rho * .75 * b[g])
  set.seed(99L)
  y <- rbinom(n, 1L, p.y)
  df <- as.data.frame(x)
  colnames(df) <- paste0("x_", seq_len(ncol(x)))
  df$y <- y
  df$g <- g
  
  rbartFit <- dbarts::rbart_vi(
    y ~ . - g, df, group.by = g,
    n.samples = 240L, n.burn = 120L, n.thin = 3L, n.chains = 2L,
    n.trees = 50L, n.threads = 1L,
    verbose = FALSE
  )
  ranef.rbart <- rbartFit$ranef.mean
  
  glmerFit <- lme4$glmer(
    y ~ . - g + (1 | g),
    df,
    family = binomial(link = "probit")
  )
  
  rbart.mu.hat <- fitted(rbartFit)
  expect_equal(apply(dbarts:::extract(rbartFit), 2L, mean), rbart.mu.hat)
  glmer.mu.hat <- fitted(glmerFit, type = "response")
  
  dev.rbart <- -2 * mean(log(ifelse(df$y == 1, rbart.mu.hat, 1 - rbart.mu.hat)))
  dev.glmer <- -2 * mean(log(ifelse(df$y == 1, glmer.mu.hat, 1 - glmer.mu.hat)))
  expect_true(dev.rbart < dev.glmer)

  rm(dev.glmer, dev.rbart, glmer.mu.hat, rbart.mu.hat, glmerFit)
  rm(ranef.rbart, rbartFit, df, y, p.y, rho)
  rm(b_rmse.lmer, b_rmse.rbart, ranef.lmer, lmerFit)
  rm(b, g, Ey, x, n, sigma, f, lme4)
}

