# negbin-r-update-mixing.R
#
# No-trees prototype for docs/design/negative-binomial.md, decision 2 (the r
# update). Isolates the dispersion (r) full-conditional update from the mean
# model: the linear predictor psi_i = logit(p_i) is HELD FIXED (the logit-p
# parameterization, where p_i does NOT depend on r), so both samplers and a
# fine grid target the identical 1-D posterior p(r | y, {p_i}) under a shared
# Gamma(a0, b0) prior. Two schemes for a REAL-valued r are compared:
#
#   (A) CRT-Gamma  : Zhou-Carin augmentation. L_i ~ CRT(y_i, r), then
#                    r | . ~ Gamma(a0 + sum L_i, rate b0 - sum log(1 - p_i)).
#                    Exact Gibbs, conjugate, integer table draws only, real r.
#   (B) MH-log-r   : random-walk Metropolis on xi = log r against the closed-
#                    form NB log-likelihood (no augmentation), symmetric
#                    proposal, target includes the exp(xi) Jacobian.
#
# WHY logit-p AND NOT log-mean: the CRT-Gamma conjugacy needs (1-p_i)^r to be
# const^r, i.e. p_i INDEPENDENT of r. That holds under the logit-p param
# (psi = f + offset, p_i fixed by the fit) but NOT under log-mean (log mu =
# f + offset => p_i = mu_i/(mu_i+r) moves with r, and (r/(mu_i+r))^r is not a
# Gamma kernel). So the two schemes only target the SAME posterior under
# logit-p; that is the setting compared here. Under log-mean only MH-log-r is
# valid for the r-update. This dependence is decision-1<->decision-2 coupling,
# documented in the note.
#
# ROLE IN THE NOTE (post-review): the note's fork (A) recommendation (integer
# r on a capped grid, exact) supersedes both schemes for v1; this study is the
# evidence base behind the fork (B) real-r door and for the CRT > MH ranking
# that fixes (B)'s r update. It compares r-update KERNELS in isolation.
#
# HONEST CAVEATS: (1) holding psi fixed removes the r-vs-mean-level confounding
# a real BART fit induces (mean = r exp(psi): r and the fit's level both scale
# the mean; analogous to robust-errors' lambda-sigma coupling and ordinal's
# f-vs-cutpoint ridge). The ESS numbers here are therefore an OPTIMISTIC bound
# on in-sampler r mixing, not a dbarts prediction; the transferable findings
# are the RANKING and the correctness agreement with the grid. (2) This
# isolation never draws omega at all, so it can witness NEITHER the
# fractional-PG truncation bias NOR sweep-ordering (r-vs-omega) bugs - it
# validates the r kernels, nothing about their composition into the sweep.
# (3) Each (n, r) cell is a single replicate.
#
# Self-contained, seed-fixed. Run: Rscript negbin-r-update-mixing.R

set.seed(20260718)

## ---- prior (shared by both samplers and the grid) -------------------------
## Gamma(shape a0, rate b0); the robust-errors tail-bounding precedent
## (gamma(2, 0.1)) adapted to a dispersion. Weakly informative, proper.
a0 <- 2.0
b0 <- 0.1

## ---- ESS via Geyer initial-positive-sequence ------------------------------
ess_ips <- function(x) {
  x <- x - mean(x)
  n <- length(x)
  v0 <- sum(x * x) / n
  if (v0 <= 0) return(NA_real_)
  acf_full <- as.numeric(acf(x, lag.max = n - 1, plot = FALSE, demean = FALSE)$acf)
  ## pair up (Gamma_k = rho_{2k} + rho_{2k+1}); stop when a pair goes <= 0
  gsum <- 0
  k <- 1
  while (2 * k + 1 < n) {
    g <- acf_full[2 * k] + acf_full[2 * k + 1]  # acf_full[1] is lag 0
    if (g <= 0) break
    gsum <- gsum + g
    k <- k + 1
  }
  tau <- 1 + 2 * gsum          # integrated autocorrelation time (rho_0 = 1)
  n / tau
}

## ---- CRT draw: L_i ~ CRT(y_i, r), vectorized over observations by table j --
draw_CRT_total <- function(y, r) {
  L <- numeric(length(y))
  maxy <- max(y)
  if (maxy < 1) return(0)
  for (j in 1:maxy) {
    active <- y >= j
    pj <- r / (r + (j - 1))
    L[active] <- L[active] + (runif(sum(active)) < pj)
  }
  sum(L)
}

## ---- closed-form NB log-posterior in r (p fixed, logit-p param) -----------
## y_i ~ NB(r, p_i); log lik = sum lgamma(y+r) - lgamma(r) - lgamma(y+1)
##  + r log(1-p_i) + y log p_i, with p_i FIXED (independent of r).
logpost_r <- function(r, y, p) {
  if (r <= 0) return(-Inf)
  ll <- sum(lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
            r * log1p(-p) + y * log(p))
  ll + (a0 - 1) * log(r) - b0 * r          # + Gamma(a0,b0) log-prior kernel
}

## ---- one experimental cell ------------------------------------------------
run_cell <- function(n, r_true, n_iter = 6000L, burn = 1500L) {
  ## logit-p model: psi_i = logit(p_i) fixed (does NOT depend on r); the mean
  ## r*exp(psi_i) varies with r, which is precisely the r-vs-level confounding
  ## a real fit would carry and this isolation deliberately removes.
  x   <- rnorm(n)
  psi <- -0.4 + 0.5 * x                    # spread of p_i in ~(0.2, 0.8)
  p   <- plogis(psi)                       # p_i fixed
  y   <- rnbinom(n, size = r_true, prob = 1 - p)  # R prob = 1 - p (our p)
  sumLog1mp <- sum(log1p(-p))              # sum log(1-p_i) (< 0), fixed

  ## ---- (A) CRT-Gamma Gibbs ----
  rA <- numeric(n_iter); r_cur <- 1.0
  for (t in 1:n_iter) {
    Ltot  <- draw_CRT_total(y, r_cur)
    shape <- a0 + Ltot
    rate  <- b0 - sumLog1mp                # b0 + sum(-log(1-p_i)) > 0
    r_cur <- rgamma(1, shape = shape, rate = rate)
    rA[t] <- r_cur
  }

  ## ---- (B) MH on xi = log r ----
  rB <- numeric(n_iter); xi <- log(1.0)
  lp <- logpost_r(exp(xi), y, p) + xi      # target in xi carries the Jacobian
  prop_sd <- 0.4; acc <- 0
  adapt_to <- burn
  for (t in 1:n_iter) {
    xi_p <- xi + prop_sd * rnorm(1)
    lp_p <- logpost_r(exp(xi_p), y, p) + xi_p
    if (log(runif(1)) < lp_p - lp) { xi <- xi_p; lp <- lp_p; acc <- acc + 1 }
    rB[t] <- exp(xi)
    ## diminishing adaptation of the proposal sd during burn-in only (frozen after)
    if (t <= adapt_to) {
      target <- 0.44
      prop_sd <- exp(log(prop_sd) + (1 / sqrt(t)) *
                     ((acc / t) - target))
    }
  }

  ## ---- fine-grid exact posterior (same prior) ----
  grid <- seq(0.02, 60, length.out = 8000)
  lg <- vapply(grid, logpost_r, numeric(1), y = y, p = p)
  w  <- exp(lg - max(lg)); w <- w / sum(w)
  grid_mean <- sum(grid * w)
  grid_sd   <- sqrt(sum((grid - grid_mean)^2 * w))

  keep <- (burn + 1):n_iter
  list(n = n, r_true = r_true,
       ess_CRT = ess_ips(rA[keep]), ess_MH = ess_ips(rB[keep]),
       mh_accept = acc / n_iter,
       mean_CRT = mean(rA[keep]), mean_MH = mean(rB[keep]),
       mean_grid = grid_mean, sd_grid = grid_sd)
}

## ---- run the grid of cells ------------------------------------------------
cells <- expand.grid(n = c(200L, 2000L), r_true = c(0.5, 2, 10))
res <- lapply(seq_len(nrow(cells)),
              function(i) run_cell(cells$n[i], cells$r_true[i]))

cat(sprintf("prior: Gamma(shape=%.2f, rate=%.2f); 6000 iter, 1500 burn; seed 20260718\n\n",
            a0, b0))
cat(sprintf("%5s %6s | %8s %8s %7s | %8s %8s %8s %8s\n",
            "n", "r", "ESS_CRT", "ESS_MH", "MHacc",
            "post_CRT", "post_MH", "grid", "grid_sd"))
for (r in res)
  cat(sprintf("%5d %6.1f | %8.0f %8.0f %7.2f | %8.3f %8.3f %8.3f %8.3f\n",
              r$n, r$r_true, r$ess_CRT, r$ess_MH, r$mh_accept,
              r$mean_CRT, r$mean_MH, r$mean_grid, r$sd_grid))

## ---- Part B: real-shape PG(b, z) truncation-error characterization --------
## With real r the NB-PG mean update needs omega_i ~ PG(y_i + r, psi_i) with
## non-integer shape. The shipped Devroye sampler draws PG(1, z) only; integer
## b is a sum of b PG(1, z) draws (model.hpp:2544-2546). The candidate
## fractional method is the Devroye/Zhou gamma-sum TRUNCATION:
##   PG(a, z) = (1/(2 pi^2)) sum_{k=1}^{K} g_k / ((k-1/2)^2 + z^2/(4 pi^2)),
##   g_k ~ Gamma(a, 1),
## which is APPROXIMATE: truncation drops a nonnegative tail, so draws are
## biased LOW by < a/(2 pi^2 (K - 1/2)) in the mean (relative ~ 2/(pi^2 K) at
## small z, ~1e-3 at K = 200). This check CHARACTERIZES that error against the
## exact mean E[PG(b, z)] = (b/(2 z)) tanh(z/2): observed relative deviations
## (<~1.2% at 4000 draws) are consistent with the ~0.1% one-sided bias bound
## plus MC noise. It does NOT certify an exact sampler - no exact real-shape
## PG sampler is established (note section 2B) - and here even the integer
## part uses the truncated sum (the shipped design would use exact Devroye
## for it). Evidence for the note's fork (B) error budget.
pg_gamma_sum <- function(a, z, K = 200L) {
  k <- 1:K
  denom <- (k - 0.5)^2 + (z^2) / (4 * pi^2)
  g <- rgamma(K, shape = a, rate = 1)
  sum(g / denom) / (2 * pi^2)
}
pg_real <- function(b, z, K = 200L) {
  ib <- floor(b); fb <- b - ib
  s <- 0
  if (ib >= 1) for (i in 1:ib) s <- s + pg_gamma_sum(1, z, K)
  if (fb > 0) s <- s + pg_gamma_sum(fb, z, K)
  s
}
pg_mean_exact <- function(b, z) (b / (2 * z)) * tanh(z / 2)

cat("\nPart B: real-shape PG(b,z) gamma-sum composition, mean check (K=200 terms)\n")
cat(sprintf("%6s %5s | %10s %10s %8s\n", "b", "z", "emp_mean", "exact", "relerr"))
for (bz in list(c(2.5, 1.0), c(0.5, 0.8), c(10.7, 1.5), c(3.3, 0.3))) {
  b <- bz[1]; z <- bz[2]
  draws <- replicate(4000, pg_real(b, z))
  em <- mean(draws); ex <- pg_mean_exact(b, z)
  cat(sprintf("%6.1f %5.1f | %10.4f %10.4f %8.4f\n",
              b, z, em, ex, abs(em - ex) / ex))
}
