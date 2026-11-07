## Exact toy: theta in {1,2} (trees), a in {0,1} (row j active), omega in {1,2} (row j's latent).
## Target: pi(theta, a=0) = p0(theta); pi(theta, a=1, omega) = p1(theta) c(omega | theta).
## Sweep: latent step (IMH for c(.|theta) if a = 1; frozen if a = 0), tree step
## theta ~ pi(theta | a, omega) (inactive row ignored), then the caller's activity update.
set.seed(3)
p0 <- c(0.15, 0.25); p1 <- c(0.35, 0.25)               # sums to 1
cc <- rbind(c(0.7, 0.3), c(0.2, 0.8))                  # c(omega | theta), rows theta
g  <- c(0.45, 0.55)                                    # pseudo-prior, no theta
qp <- rbind(c(0.5, 0.5), c(0.6, 0.4))                  # IMH proposal given theta
S <- expand.grid(th = 1:2, a = 0:1, om = 1:2); n <- nrow(S)
idx <- function(th, a, om) which(S$th == th & S$a == a & S$om == om)
imh <- function(th) {                                  # 2x2 IMH kernel on omega for c(.|th)
  K <- matrix(0, 2, 2)
  for (x in 1:2) for (y in 1:2) if (x != y)
    K[x, y] <- qp[th, y] * min(1, cc[th, y] * qp[th, x] / (cc[th, x] * qp[th, y]))
  diag(K) <- 1 - rowSums(K); K
}
Klat <- matrix(0, n, n); Ktree <- Klat
for (s in 1:n) {
  th <- S$th[s]; a <- S$a[s]; om <- S$om[s]
  if (a == 0) Klat[s, s] <- 1 else for (y in 1:2) Klat[s, idx(th, 1, y)] <- imh(th)[om, y]
  w <- if (a == 0) p0 else p1 * cc[, om]; w <- w / sum(w)
  for (t in 1:2) Ktree[s, idx(t, a, om)] <- w[t]
}
act <- function(rule) {
  K <- matrix(0, n, n)
  for (s in 1:n) {
    th <- S$th[s]; om <- S$om[s]
    w <- if (rule == "ext") c(p0[th] * g[om], p1[th] * cc[th, om])  # given theta and every omega
         else c(p0[th], p1[th])                                    # omega integrated out
    w <- w / sum(w)
    K[s, idx(th, 0, om)] <- w[1]; K[s, idx(th, 1, om)] <- w[2]
  }
  K
}
stat <- function(P) { e <- eigen(t(P)); v <- Re(e$vectors[, 1]); v / sum(v) }
pext <- ifelse(S$a == 0, p0[S$th] * g[S$om], p1[S$th] * cc[cbind(S$th, S$om)])
for (rule in c("ext", "int")) {
  P <- Klat %*% Ktree %*% act(rule); st <- stat(P)
  m <- tapply(st, list(S$th, S$a), sum); mt <- cbind(p0, p1)
  cat(rule, ": max |pi_ext P - pi_ext| =", signif(max(abs(pext %*% P - pext)), 3),
      " TV of (theta, a) marginal from target =", signif(sum(abs(m - mt)) / 2, 3), "\n")
}
## the engine's other steps leave pi_ext invariant on their own
cat("latent step:", signif(max(abs(pext %*% Klat - pext)), 3),
    " tree step:", signif(max(abs(pext %*% Ktree - pext)), 3), "\n")
