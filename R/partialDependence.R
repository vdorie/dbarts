## create the contents to be used in partial dependence plots
pdbart <- function (
   x.train, y.train,
   xind = seq_len(ncol(x.train)),
   levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
   pl = TRUE, plquants = c(0.05, 0.95),
   ...
)
{
   n = nrow(x.train)
   nvar = length(xind)
   nlevels = rep(0,nvar)
   if(is.null(levs)) {
      levs = list()
      for(i in 1:nvar) {
         ux = unique(x.train[,xind[i]])
	 if(length(ux) < length(levquants)) levs[[i]] = sort(ux)
	 else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
      }
   } 
   nlevels = unlist(lapply(levs,length))
   x.test=NULL
   for(i in 1:nvar) {
      for(v in levs[[i]]) {
         temp = x.train
         temp[,xind[i]] = v
         x.test = rbind(x.test,temp)
      }
   }
   pdbrt = bart(x.train,y.train,x.test,...)
   fdr = list() 
   cnt=0
   for(j in 1:nvar) {
      fdrtemp=NULL
      for(i in 1:nlevels[j]) {
         cind = cnt + ((i-1)*n+1):(i*n)
         fdrtemp = cbind(fdrtemp,(apply(pdbrt$yhat.test[,cind],1,mean)))
      }
      fdr[[j]] = fdrtemp
      cnt = cnt + n*nlevels[j]
   }
   if(is.null(colnames(x.train))) xlbs = paste('x',xind,sep='')
   else xlbs = colnames(x.train)[xind]
   if('sigma' %in% names(pdbrt)) {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      first.sigma=pdbrt$first.sigma,sigma=pdbrt$sigma,
      yhat.train.mean=pdbrt$yhat.train.mean,sigest=pdbrt$sigest,y=pdbrt$y)
   } else {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      y=pdbrt$y)
   }
   class(retval) = 'pdbart'
   if(pl) plot(retval,plquants=plquants)
   return(retval)
}

pd2bart <- function (
   x.train, y.train,
   xind = c(1, 2),
   levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
   pl = TRUE, plquants = c(0.05, 0.95), 
   ...
)
{
   n = nrow(x.train)
   nlevels = rep(0,2)
   if(is.null(levs)) {
      levs = list()
      for(i in 1:2) {
         ux = unique(x.train[,xind[i]])
	 if(length(ux) <= length(levquants)) levs[[i]] = sort(ux)
	 else levs[[i]] = unique(quantile(x.train[,xind[i]],probs=levquants))
      }
   } 
   nlevels = unlist(lapply(levs,length))
   xvals <- as.matrix(expand.grid(levs[[1]],levs[[2]]))
   nxvals <- nrow(xvals)
   if (ncol(x.train)==2){
      cat('special case: only 2 xs\n')
      x.test = xvals
   } else {
      x.test=NULL
      for(v in 1:nxvals) {
         temp = x.train
         temp[,xind[1]] = xvals[v,1]
         temp[,xind[2]] = xvals[v,2]
         x.test = rbind(x.test,temp)
      }
   }
   pdbrt = bart(x.train,y.train,x.test,...)
   if (ncol(x.train)==2) {
      fdr = pdbrt$yhat.test
   } else {
      fdr = NULL 
      for(i in 1:nxvals) {
         cind =  ((i-1)*n+1):(i*n)
         fdr = cbind(fdr,(apply(pdbrt$yhat.test[,cind],1,mean)))
      }
   }
   if(is.null(colnames(x.train))) xlbs = paste('x',xind,sep='')
   else xlbs = colnames(x.train)[xind]
   if('sigma' %in% names(pdbrt)) {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      first.sigma=pdbrt$first.sigma,sigma=pdbrt$sigma,
      yhat.train.mean=pdbrt$yhat.train.mean,sigest=pdbrt$sigest,y=pdbrt$y)
   } else {
   retval = list(fd = fdr,levs = levs,xlbs=xlbs,
      bartcall=pdbrt$call,yhat.train=pdbrt$yhat.train,
      y=pdbrt$y)
   }
   class(retval) = 'pd2bart'
   if(pl) plot(retval,plquants=plquants)
   return(retval)
}
