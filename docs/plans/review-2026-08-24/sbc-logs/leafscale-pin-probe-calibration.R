suppressMessages(library(dbarts))
set.seed(11L); n <- 150L; p <- 3L
x <- matrix(runif(n * p), n, p)
ctl <- dbartsControl(n.trees = 50L, n.chains = 1L, n.threads = 1L,
                     updateState = FALSE, verbose = FALSE)
yA <- seq(-2.5, 2.5, length.out = n)          # anchor build vector, range 5
yB <- 4 * yA + 7                              # a "replication" y, range 20
cal <- function(s) s$getCalibration(1L)[1L, c("prior.scale","prior.sd","prior.mean","response.scale","response.shift")]

cat("--- gaussian, DEFAULT node prior (inherited calibration) ---\n")
a0 <- dbarts(x, yA, control = ctl); b0 <- dbarts(x, yB, control = ctl)
print(rbind(anchor = cal(a0), rebuild = cal(b0)))

S <- unname(cal(a0)["prior.scale"])
cat("\n--- gaussian, NAMED node.prior = normal(2, scale =", S, ") ---\n")
a1 <- dbarts(x, yA, node.prior = normal(2, scale = S), control = ctl)
b1 <- dbarts(x, yB, node.prior = normal(2, scale = S), control = ctl)
print(rbind(anchor = cal(a1), rebuild = cal(b1)))
cat("identical prior.scale across rebuilds: ",
    identical(cal(a1)[["prior.scale"]], cal(b1)[["prior.scale"]]), "\n")
cat("bart2 formal: ")
print(names(formals(bart2))[grepl("prior", names(formals(bart2)))])
cat("model slot after creation: a1@model@prior.scale = ", a1$model@prior.scale, "\n")

cat("\n--- setResponse(updateScale = FALSE) vs TRUE on the anchor ---\n")
a2 <- dbarts(x, yA, control = ctl)
before <- cal(a2); a2$setResponse(yB); mid <- cal(a2)
a3 <- dbarts(x, yA, control = ctl); a3$setResponse(yB, updateScale = TRUE); aft3 <- cal(a3)
print(rbind(created = before, `setResponse FALSE` = mid, `setResponse TRUE` = aft3))

cat("\n--- setData (whole-data swap) ---\n")
a4 <- dbarts(x, yA, control = ctl)
a4$setData(dbartsData(x, yB))
print(rbind(created = before, `after setData` = cal(a4)))

cat("\n--- setData under a NAMED calibration ---\n")
a5 <- dbarts(x, yA, node.prior = normal(2, scale = S), control = ctl)
a5$setData(dbartsData(x, yB))
print(rbind(`named, created` = cal(a1), `named, after setData` = cal(a5)))
a5$setModel(a5$model)
cat("after setModel(model) re-issue: ", cal(a5)[["prior.scale"]], "\n")
