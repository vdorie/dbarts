plot.bart <- function(
   x,
   plquants = c(0.05, 0.95), cols = c('blue','black'),
   ...
)
{
  if ("sigma" %in% names(x)) {
    par(mfrow = c(1L, 2L))
    if (!is.null(dim(x$sigma))) {
      plot(NULL, type = "n", ylab = "sigma",
           xlim = c(1, ncol(x$first.sigma) + ncol(x$sigma)), ylim = range(x$first.sigma, x$sigma))
      for (i in seq_len(nrow(x$sigma))) {
        lines(c(seq_len(ncol(x$first.sigma)), ncol(x$first.sigma) + 0.5),
              c(x$first.sigma[i,], 0.5 * (x$first.sigma[i, ncol(x$first.sigma)] + x$sigma[i,1L])), col = "red", lty = i)
        lines(c(ncol(x$first.sigma) + 0.5, seq.int(ncol(x$first.sigma) + 1, length.out = ncol(x$sigma))),
              c(0.5 * (x$first.sigma[i, ncol(x$first.sigma)] + x$sigma[i,1L]), x$sigma[i,]), lty = i)
      }
    } else {
      plot(c(x$first.sigma, x$sigma),
           col = rep(c("red", "black"), c(length(x$first.sigma), length(x$sigma))),
           ylab = "sigma", ...)
    }
  }

  if('sigma' %in% names(x)) {
     ql <- apply(x$yhat.train, length(dim(x$yhat.train)), quantile,probs=plquants[1])
     qm <- apply(x$yhat.train, length(dim(x$yhat.train)), quantile,probs=.5)
     qu <- apply(x$yhat.train, length(dim(x$yhat.train)), quantile,probs=plquants[2])
     plot(x$y,qm,ylim=range(ql,qu),xlab='y',ylab=
      'posterior interval for E(Y|x)',...)
     for (i in 1:length(qm))
       lines(rep(x$y[i],2),c(ql[i],qu[i]),col=cols[1])
     abline(0,1,lty=2,col=cols[2])
  } else {
     pdrs = pnorm(x$yhat.train) #draws of p(Y=1 | x)
     ql <- apply(pdrs, length(dim(pdrs)), quantile,probs=plquants[1])
     qm <- apply(pdrs, length(dim(pdrs)), quantile,probs=.5)
     qu <- apply(pdrs, length(dim(pdrs)), quantile,probs=plquants[2])
     plot(qm,qm,ylim=range(ql,qu),xlab='meadian of p',ylab=
      'posterior interval for P(Y=1|x)',...)
     for (i in 1:length(qm))
       lines(rep(qm[i],2),c(ql[i],qu[i]),col=cols[1])
     abline(0,1,lty=2,col=cols[2])
  }
}

plot.rbart <- function(x, plquants = c(0.05, 0.95), cols = c('blue','black'), ...)
{
  if ("sigma" %in% names(x)) {
    par(mfrow = c(1L, 2L))
    if (!is.null(dim(x$sigma))) {
      plot(NULL, type = "n", ylab = "sigma",
           xlim = c(1, ncol(x$first.sigma) + ncol(x$sigma)), ylim = range(x$first.sigma, x$sigma))
      for (i in seq_len(nrow(x$sigma))) {
        lines(c(seq_len(ncol(x$first.sigma)), ncol(x$first.sigma) + 0.5),
              c(x$first.sigma[i,], 0.5 * (x$first.sigma[i, ncol(x$first.sigma)] + x$sigma[i,1L])), col = "red", lty = i)
        lines(c(ncol(x$first.sigma) + 0.5, seq.int(ncol(x$first.sigma) + 1, length.out = ncol(x$sigma))),
              c(0.5 * (x$first.sigma[i, ncol(x$first.sigma)] + x$sigma[i,1L]), x$sigma[i,]), lty = i)
      }
    } else {
      plot(c(x$first.sigma, x$sigma),
           col = rep(c("red", "black"), c(length(x$first.sigma), length(x$sigma))),
           ylab = "sigma", ...)
    }
  }
  
  if (length(dim(x$ranef)) > 2L) 
    ranef <- x$ranef[,,as.integer(x$group.by)]
  else
    ranef <- x$ranef[,as.integer(x$group.by)]
  yhat.train <- x$yhat.train + ranef
  
  if ("sigma" %in% names(x)) {
    ql <- apply(yhat.train, length(dim(yhat.train)), quantile, probs = plquants[1L])
    qm <- apply(yhat.train, length(dim(yhat.train)), quantile, probs = .5)
    qu <- apply(yhat.train, length(dim(yhat.train)), quantile, probs = plquants[2L])
    plot(x$y, qm, ylim = range(ql, qu), xlab = "y", ylab= "posterior interval for E(Y | x)", ...)
    
    for (i in seq_along(qm)) lines(rep(x$y[i], 2L), c(ql[i], qu[i]), col = cols[1L])
    abline(0, 1, lty = 2L, col = cols[2L])
  } else {
    ## shouldn't happen for now
    pdrs <- pnorm(yhat.train) #draws of p(Y=1 | x)
    ql <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[1L])
    qm <- apply(pdrs, length(dim(pdrs)), quantile, probs = .5)
    qu <- apply(pdrs, length(dim(pdrs)), quantile, probs = plquants[2L])
    plot(qm, qm, ylim = range(ql, qu), xlab = "median of p",
         ylab= "posterior interval for P(Y = 1|  x)", ...)
    for (i in seq_along(qm)) lines(rep(qm[i], 2L), c(ql[i], qu[i]), col = cols[1L])
    abline(0, 1, lty = 2L, col = cols[2L])
  }
}


plot.pdbart <- function(
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

plot.pd2bart <- function(
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

