set.seed(20260824)
n <- 40L; nTest <- 12L
x <- matrix(runif(n * 3L), n, 3L); colnames(x) <- c("x1","x2","x3")
xTest <- matrix(runif(nTest * 3L), nTest, 3L); colnames(xTest) <- colnames(x)
group <- factor(sample(letters[1:4], n, TRUE))
groupTest <- factor(sample(letters[1:4], nTest, TRUE))
yGaussian <- as.numeric(2*x[,1] - x[,2] + rnorm(n, 0, 0.3))
pBin <- plogis(2*(x[,1]-0.5))
yBinary <- factor(rbinom(n, 1L, pBin), levels = c(0,1))
etaM <- cbind(2*(x[,1]-0.5), x[,2]-x[,3], 0); pM <- exp(etaM); pM <- pM/rowSums(pM)
labM <- vapply(seq_len(n), function(i) sample.int(3L,1L,prob=pM[i,]), integer(1))
yMultinom <- factor(c("a","b","c")[labM], levels=c("a","b","c"))
yOrdinal <- ordered(c("lo","mid","hi")[1L+(x[,1]>0.33)+(x[,1]>0.66)], levels=c("lo","mid","hi"))
yNbinom <- rnbinom(n, size=5L, mu=exp(0.4*(x[,1]-0.5)))
survTime <- rexp(n, rate=exp(-0.4*(x[,1]-0.5))); survStatus <- rbinom(n,1L,0.85)
ySurv <- survival::Surv(survTime, survStatus)
yHurdle <- ifelse(runif(n) < 0.3, 0, rlnorm(n, meanlog = 0.2*(x[,1]-0.5)))
zAmp <- rbinom(n, 1L, 0.5)
ampDf <- data.frame(y = yGaussian, x1 = x[,1], x2 = x[,2], x3 = x[,3], z = zAmp)

b2 <- function(y, ..., test = xTest) {
  args <- list(x, y, n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
               n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE,
               verbose = FALSE, seed = 1L)
  if (!is.null(test)) args$test <- test
  do.call(bart2, utils::modifyList(args, list(...)))
}
rb <- function(y, ...) {
  args <- list(x, y, group.by = group, group.by.test = groupTest, test = xTest,
               n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
               n.threads = 1L, keepTrees = TRUE, verbose = FALSE)
  do.call(rbart_vi, utils::modifyList(args, list(...)))
}

# name -> zero-arg builder.  conditions: base / kT0 (keepTrees+keepSampler FALSE)
# / nc2 (2 chains, combined) / nc2u (2 chains, uncombined)
builders <- list(
  gaussian.base      = function() b2(yGaussian),
  gaussian.kT0       = function() b2(yGaussian, keepTrees = FALSE, keepSampler = FALSE),
  gaussian.nc2       = function() b2(yGaussian, n.chains = 2L),
  gaussian.nc2u      = function() b2(yGaussian, n.chains = 2L, combineChains = FALSE),
  gaussian.nc2.kT0   = function() b2(yGaussian, n.chains = 2L, keepTrees = FALSE, keepSampler = FALSE),
  probit.base        = function() b2(yBinary, family = "probit"),
  probit.kT0         = function() b2(yBinary, family = "probit", keepTrees = FALSE, keepSampler = FALSE),
  probit.nc2u        = function() b2(yBinary, family = "probit", n.chains = 2L, combineChains = FALSE),
  logistic.base      = function() b2(yBinary, family = "logistic"),
  logistic.kT0       = function() b2(yBinary, family = "logistic", keepTrees = FALSE, keepSampler = FALSE),
  aft.base           = function() b2(ySurv, family = "aft"),
  aft.kT0            = function() b2(ySurv, family = "aft", keepTrees = FALSE, keepSampler = FALSE),
  hazard.base        = function() b2(ySurv, family = "hazard", test = NULL),
  hazard.kT0         = function() b2(ySurv, family = "hazard", test = NULL, keepTrees = FALSE, keepSampler = FALSE),
  hetero.base        = function() b2(yGaussian, variance = ~x1),
  hetero.kT0         = function() b2(yGaussian, variance = ~x1, keepTrees = FALSE, keepSampler = FALSE),
  hetero.nc2u        = function() b2(yGaussian, variance = ~x1, n.chains = 2L, combineChains = FALSE),
  student.base       = function() bart2(x, yGaussian, n.trees = 4L, n.samples = 6L,
                          n.burn = 4L, n.chains = 1L, n.threads = 1L, keepTrees = TRUE,
                          keepTrainingFits = TRUE, test = xTest, resid.dist = student(),
                          verbose = FALSE, seed = 1L),
  amplitude.base     = function() bart2(y ~ x1 + x2 + z:forest(x1 + x2), ampDf,
                          n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
                          n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE,
                          verbose = FALSE, seed = 1L),
  amplitude.kT0      = function() bart2(y ~ x1 + x2 + z:forest(x1 + x2), ampDf,
                          n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 1L,
                          n.threads = 1L, keepTrees = FALSE, keepSampler = FALSE,
                          keepTrainingFits = TRUE, verbose = FALSE, seed = 1L),
  amplitude.nc2u     = function() bart2(y ~ x1 + x2 + z:forest(x1 + x2), ampDf,
                          n.trees = 4L, n.samples = 6L, n.burn = 4L, n.chains = 2L,
                          n.threads = 1L, keepTrees = TRUE, keepTrainingFits = TRUE,
                          combineChains = FALSE, verbose = FALSE, seed = 1L),
  multinomial.base   = function() b2(yMultinom, family = "multinomial"),
  multinomial.kT0    = function() b2(yMultinom, family = "multinomial", keepTrees = FALSE, keepSampler = FALSE),
  multinomial.nc2    = function() b2(yMultinom, family = "multinomial", n.chains = 2L),
  multinomial.nc2u   = function() b2(yMultinom, family = "multinomial", n.chains = 2L, combineChains = FALSE),
  ordinal.base       = function() b2(yOrdinal, family = "ordinal"),
  ordinal.kT0        = function() b2(yOrdinal, family = "ordinal", keepTrees = FALSE, keepSampler = FALSE),
  ordinal.nc2u       = function() b2(yOrdinal, family = "ordinal", n.chains = 2L, combineChains = FALSE),
  nbinom.base        = function() b2(yNbinom, family = "nbinom"),
  nbinom.kT0         = function() b2(yNbinom, family = "nbinom", keepTrees = FALSE, keepSampler = FALSE),
  nbinom.nc2u        = function() b2(yNbinom, family = "nbinom", n.chains = 2L, combineChains = FALSE),
  hurdle.base        = function() b2(yHurdle, family = "hurdle.lognormal", test = NULL),
  hurdle.kT0         = function() b2(yHurdle, family = "hurdle.lognormal", test = NULL, keepTrees = FALSE, keepSampler = FALSE),
  hurdle.nc2u        = function() b2(yHurdle, family = "hurdle.lognormal", test = NULL, n.chains = 2L, combineChains = FALSE),
  rbart.base         = function() rb(yGaussian),
  rbart.kT0          = function() rb(yGaussian, keepTrees = FALSE),
  rbart.nc2          = function() rb(yGaussian, n.chains = 2L),
  rbart.nc2u         = function() rb(yGaussian, n.chains = 2L, combineChains = FALSE),
  rbart.aft          = function() rbart_vi(ySurv ~ x1 + x2, data.frame(x1=x[,1], x2=x[,2]),
                          group.by = group, n.trees = 4L, n.samples = 6L, n.burn = 4L,
                          n.chains = 1L, n.threads = 1L, keepTrees = TRUE, verbose = FALSE),
  xbart.result       = function() xbart(x, yGaussian, n.samples = 5L, n.reps = 2L, n.test = 5L,
                          n.burn = c(4L,3L), n.trees = 4L, n.threads = 1L,
                          method = "k-fold", verbose = FALSE)
)
