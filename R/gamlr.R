##### estimation for gamma-lasso penalized regression  ######

## Wrapper function; most happens in .c and glpath
gamlr <- function(covars, response, family="linear",
                  free=NULL, scale=TRUE, fixe=NULL, store=TRUE,
                  penvar=1, ortho=TRUE, step=NULL,
                  stoparg=list(rule="BF",val=10),
                  cdpar=list(tol=1e-5), verb=FALSE, ...)
{
    on.exit(.C("gamlr_cleanup", PACKAGE = "gamlr"))

    ## check and clean all arguments
    chk <- glcheck(covars=covars, response=response, family=family,
                   fixe=fixe, scale=scale, free=free,
                   cdpar=cdpar, verb=verb, store=store)

    path <- c(list(...),
              penvar=penvar, ortho=ortho, step=step, list(stoparg=stoparg))
    out <- glpath(chk$xlist, chk$ylist, chk$algo, path)

    ## add scaling info
    out$scaled=scale
    out$covarMean=chk$xlist$covarMean
    out$covarSD=chk$xlist$covarSD

    ## class and return
    class(out) <- "gamlr"
    invisible(out)
}


###### implemented methods for gamlr objects ######

## s3 plot function
plot.gamlr<- function(x, against="logpen", select=NULL, ...)
{
  if(!is.null(x$path) && ncol(x$load)>1){
    p <- ncol(x$load)
    nzr <- unique(x$load$i)
    nzr <- nzr[!(nzr%in%x$free)]
    cols <- rainbow(length(nzr))
    names(cols) <- nzr
    if(against=="segment"){
      xv <- 1:p
      xvn <- "path segment"
    } else if(against=="pve"){
      xv <- x$pve
      xvn <- "pve"
    } else if(against=="pen"){
      xv <- 1/x$penalty
      xvn <- "1/E[lambda]"
    } else if(against=="logpen"){
      xv <- log(1/x$penalty)
      xvn <- "-log(E[lambda])"
    } else
    stop("unrecognized 'against' argument.  options are pve,pen,logpen.")

    argl = list(...)
    if(is.null(argl$ylim)) argl$ylim=range(x$load[nzr,])
    if(is.null(argl$ylab)) argl$ylab="loading"
    if(is.null(argl$xlab)) argl$xlab=xvn
    if(!is.null(argl$col)){
      cols[1:length(cols)] <- argl$col
      argl$col <- NULL }
    do.call(plot, c(list(x=xv, y=rep(0,p), col="grey70", type="l"), argl))
    for(i in nzr) lines(xv, c(as.matrix(x$load[i,])), col=cols[paste(i)])

    if(!is.null(select)){
      seli <- glselect(x, select)
      abline(v=xv[seli],col="grey20",lty=2) }
  } else
  plot(x$fitted ~ x$response, ...)
}

## S3 method coef function
coef.gamlr <- function(object, origscale=TRUE, select=NULL, ...){
  if(!is.null(select) && !is.null(object$path)){
    seli <- glselect(object,select)
    loads <- object$loadings[,seli]
    ind <- loads$i
    coef <- matrix(c(object$intercept[seli],loads$v),
                   dimnames=list(c("intercept", rownames(loads)[ind])))
  } else{
    ind <- 1:(ncol(object$X))
    coef <- rbind(matrix(object$intercept,nrow=1,dimnames=list("intercept")),
                  as.matrix(object$loadings)) }

  if(origscale && !is.null(object$covarSD)){
    coef[-1,] <-  coef[-1,]/object$covarSD[ind]
    coef[1,] <- coef[1,] - col_sums(object$covarMean[ind]*coef[-1,,drop=FALSE])
  }
  return( coef )
}

## S3 method predict function
predict.gamlr <- function(object, newdata, select=NULL, ...)
{
  p <- nrow(object$loadings)
  if(is.vector(newdata))
    newdata <- matrix(newdata, nrow=1)

  if(ncol(newdata)!=p)
    stop(sprintf("newdata must be a matrix with %d columns",p))

  newdata <- as.simple_triplet_matrix(cbind(rep(1,nrow(newdata)),newdata))

  B <- coef(object, select=select, ...)

  if(!is.null(select)){
    BS <- rep(0,p+1)
    names(BS) <- c("intercept",rownames(object$loadings))
    BS[rownames(B)] <- B
    B <- matrix(BS)
  } else{ colnames(B) <- object$penalty }

  return(tcrossprod_simple_triplet_matrix(newdata, t(B)))
}


## S3 method summary function
summary.gamlr <- function(object,   top=20, ...){
  print(object)

  if(!is.null(object$path) && length(object$load$v)>0){
    bnames <- names(object$order)[1:min(length(object$order),top)]
    bsign <- rep(" ",length(bnames))
    bsign[as.matrix(object$load[bnames,ncol(object$load)]>0)] <- "+"
    bsign[as.matrix(object$load[bnames,ncol(object$load)]<0)] <- "-"

    cat(" Top variables by order of entry: \n  ")
    cat(paste(bsign,bnames,'\n '))
    cat("\n")
  }
}

print.gamlr <- function(x, ...){
  cat(sprintf("\n %s gamlr object with ", x$family))
  if(ncol(x$load)==1)
    cat(sprintf("%g%% nonzero coefficients. PVE = %g \n\n",
                round(100*length(x$loadings$v)/nrow(x$loadings),2),
                round(x$pve,2)))
  else cat(sprintf("%d path segments. \n",  ncol(x$load)))
}


## AIC
AIC.gamlr <- function(object, ..., k=2)
  return(-2*object$llhd + k*object$df)

