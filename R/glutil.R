#######  Undocumented gamlr utilities #########

## utility for penalty selection
glselect <- function(object, select="adjBF"){
  if(select=="CV"){
    if(is.null(object$cvindex)) stop("you need to run cvselect before select='CV'")
    seli <- object$cvindex }
  else if(select=="AIC") seli <- which.min(AIC(object))
  else if(select=="BIC") seli <- which.min(AIC(object, k=log(object$nobs)))
  else if(select=="BF") seli <- which.max(object$BF)
  else seli <- which.max(object$BF-object$zerodf)
  return(seli) }

## wrappers for c
glnllhd <- function(famid, eta, Y)
{
    cobj <- .C("R_nllhd",
               famid = as.integer(famid),
               l = double(1),
               n = as.integer(length(Y)),
               e = as.double(eta),
               y = as.double(Y),
               PACKAGE="gamlr",
               dup=FALSE)
    return(cobj$l)
}

glgrad <- function(famid, X, xcnt, eta, Y)
{
    cobj <- .C("R_grad",
               famid = as.integer(famid),
               G = double(ncol(X)),
               d = as.integer(ncol(X)),
               X = as.double(X$v),
               xr = as.integer(X$i-1),
               xcnt = as.integer(xcnt),
               eta = as.double(eta),
               Y = as.double(Y),
               PACKAGE="gamlr",
               dup=FALSE)
    names(cobj$G) <- colnames(X)
    return(cobj$G)
}

glcurve <- function(famid, X, xcnt=rep(0,ncol(X)+1), eta=NULL)
{
    if(is.null(X$cs)) X$cs <- rep(0,ncol(X)+1)
    cobj <- .C("R_curve",
               famid = as.integer(famid),
               H = double(ncol(X)),
               d = as.integer(ncol(X)),
               X = as.double(X$v),
               xr = as.integer(X$i-1),
               xcnt = as.integer(xcnt),
               eta = as.double(eta),
               PACKAGE="gamlr",
               dup=FALSE)
    names(cobj$H) <- colnames(X)
    return(cobj$H)
}
