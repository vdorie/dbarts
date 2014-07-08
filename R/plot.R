plot.bart = function(
   x,
   plquants = c(0.05, 0.95), cols = c('blue','black'),
   ...
)
{
   if('sigma' %in% names(x)) {
   par(mfrow=c(1,2))
   plot(c(x$first.sigma,x$sigma),col=rep(c('red','black'),
     c(length(x$first.sigma),length(x$sigma))),ylab='sigma',...)
   } else {
     par(mfrow=c(1,1))
   }

  if('sigma' %in% names(x)) {
     ql <- apply(x$yhat.train,2,quantile,probs=plquants[1])
     qm <- apply(x$yhat.train,2,quantile,probs=.5)
     qu <- apply(x$yhat.train,2,quantile,probs=plquants[2])
     plot(x$y,qm,ylim=range(ql,qu),xlab='y',ylab=
      'posterior interval for E(Y|x)',...)
     for (i in 1:length(qm))
       lines(rep(x$y[i],2),c(ql[i],qu[i]),col=cols[1])
     abline(0,1,lty=2,col=cols[2])
  } else {
     pdrs = pnorm(x$yhat.train) #draws of p(Y=1 | x)
     ql <- apply(pdrs,2,quantile,probs=plquants[1])
     qm <- apply(pdrs,2,quantile,probs=.5)
     qu <- apply(pdrs,2,quantile,probs=plquants[2])
     plot(qm,qm,ylim=range(ql,qu),xlab='meadian of p',ylab=
      'posterior interval for P(Y=1|x)',...)
     for (i in 1:length(qm))
       lines(rep(qm[i],2),c(ql[i],qu[i]),col=cols[1])
     abline(0,1,lty=2,col=cols[2])
  }
}

plot.pdbart = function(
   x,
   xind = seq_len(length(x$fd)),
   plquants = c(0.05, 0.95), cols = c('black','blue'),
   ...
)
{
   rgy = range(x$fd)
   for(i in xind) {
         tsum = apply(x$fd[[i]],2,quantile,probs=c(plquants[1],.5,plquants[2]))
         plot(range(x$levs[[i]]),rgy,type='n',xlab=x$xlbs[i],ylab='partial-dependence',...)
         lines(x$levs[[i]],tsum[2,],col=cols[1],type='b')
         lines(x$levs[[i]],tsum[1,],col=cols[2],type='b')
         lines(x$levs[[i]],tsum[3,],col=cols[2],type='b')
   }
}

plot.pd2bart = function(
   x,
   plquants = c(0.05, 0.95), contour.color='white',
   justmedian = TRUE,
   ...
)
{
   pdquants = apply(x$fd,2,quantile,probs=c(plquants[1],.5,plquants[2]))
   qq <- vector('list',3)
   for (i in 1:3) 
      qq[[i]]  <- matrix(pdquants[i,],nrow=length(x$levs[[1]]))
   if(justmedian) {
     zlim = range(qq[[2]])
     vind = c(2)
   } else {
      par(mfrow=c(1,3))
      zlim = range(qq)
      vind = 1:3
   }
   for (i in vind) {
     image(x=x$levs[[1]],y=x$levs[[2]],qq[[i]],zlim=zlim,
       xlab=x$xlbs[1],ylab=x$xlbs[2],...)
     contour(x=x$levs[[1]],y=x$levs[[2]],qq[[i]],zlim=zlim,
       ,add=TRUE,method='edge',col=contour.color)
     title(main=c('Lower quantile','Median','Upper quantile')[i])
   }
}
