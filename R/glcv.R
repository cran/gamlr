cvselect <- function(fit, nfold=10, npen=20, verb=TRUE){
  if(is.null(fit$path)) stop("you can only cross validate a path-fit object")

  ## build the grid
  n <- nrow(fit$X)
  mugrid <- exp(seq(max(log(fit$penalty)),min(log(fit$penalty)),length=npen))
  rord <- sample(1:n)
  cvpve <- matrix(0, nrow=nfold, ncol=npen)
  if(nfold>=n) nfold <- n
  ogrid <- c(ceiling(seq(1,n+1,length=nfold+1)))

  ## out-of-sample loop
  for(f in 1:nfold){
    test <- rord[((1:n) >= ogrid[f]) & ((1:n) < ogrid[f+1])]
    fitfold <- do.call(gamlr,
                       c(list(covars=fit$X[-test,],
                              response=fit$response[-test],
                              free=fit$free,
                              fixe=fit$fixe[-test],
                              family=fit$family,
                              scale=FALSE, mugrid=mugrid), fit$path))
    predfold <- predict(fitfold, newdata=fit$X[test,])
    pvefold <- apply(predfold, 2, function(f) pve(fit$response[test],f,fit$family))
    cvpve[f,match(paste(fitfold$penalty)[-1],mugrid)] <- pvefold[-1]
    if(verb) cat(sprintf("fold %d\n",f))
  }
  ## format and plot
  fit$cvpve <- as.data.frame(cvpve)
  dimnames(fit$cvpve) <- list(fold=1:nfold, pen.seg=1:ncol(cvpve))
  if(verb){
    boxplot(fit$cvpve, xlab="-log E[penalty]",ylab="PVE", xaxt="n",
            main=sprintf("%d-fold cross validation",nfold))
    tck <- floor(seq(2,ncol(fit$cvpve)-1,length=4))
    axis(side=1,at=tck,labels=round(-log(mugrid),1)[tck]) }

  ## select
  mu <- mugrid[which.max(colMeans(fit$cvpve))]
  fit$cvindex <- which.min(abs(fit$penalty - mu) + 1e5*(fit$penalty<mu))
  fit$cvgrid <- mugrid

  return(fit)
}


cvplot <- function(fit, ...){
  argl = list(...)
  if(is.null(argl$ylab)) argl$xlab="-log E[lambda]"
  if(is.null(argl$xlab)) argl$ylab=sprintf("%d-fold CV PVE",nrow(fit$cvpve))
  
  suppressWarnings(do.call(boxplot, c(list(x=fit$cvpve, xaxt="n"), argl)))
  tck <- floor(seq(2,ncol(fit$cvpve)-1,length=5))
  axis(side=1,at=tck,labels=round(-log(fit$cvgrid[tck]),1))
}
