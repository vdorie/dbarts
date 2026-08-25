suppressPackageStartupMessages(library(dbarts))
set.seed(1); n<-40L; nTest<-12L
x <- matrix(runif(n*3),n,3); colnames(x)<-c("x1","x2","x3")
xTest <- matrix(runif(nTest*3),nTest,3); colnames(xTest)<-colnames(x)
y <- as.numeric(2*x[,1]-x[,2]+rnorm(n,0,0.3))
f <- bart2(x,y,n.trees=3L,n.samples=4L,n.burn=3L,n.chains=1L,n.threads=1L,keepTrees=TRUE,keepTrainingFits=TRUE,test=xTest,verbose=FALSE,seed=1L)
saveRDS(f, file.path(Sys.getenv("D"),"f-fresh.rds"))          # nothing touched first
g <- bart2(x,y,n.trees=3L,n.samples=4L,n.burn=3L,n.chains=2L,n.threads=1L,keepTrees=TRUE,keepTrainingFits=TRUE,test=xTest,verbose=FALSE,seed=1L)
saveRDS(g, file.path(Sys.getenv("D"),"g-fresh.rds"))
s <- dbarts(x,y,control=dbartsControl(n.chains=1L,n.threads=1L,n.trees=3L,n.burn=3L,n.samples=4L,keepTrees=TRUE))
invisible(s$run())
saveRDS(s, file.path(Sys.getenv("D"),"s-fresh.rds"))
cat("saved\n")
