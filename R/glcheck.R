#######  Undocumented gamlr argument  manipulation #########

## check the inputs and bin
glcheck <- function(covars, response, family,
                    fixe=NULL, free=NULL, scale=TRUE,
                    cdpar=NULL, verb=FALSE, store=FALSE)
{

    ## data checks and list building
    ylist <- glresponse(response, fixe, family)
    xlist <- glcovars(covars, scale, free)
    algo <- do.call(algochk, c(cdpar,verb=verb,store=store))
    
    ## names and check
    if(length(ylist$Y) != nrow(xlist$X))
        stop("different number of predictor and response observations")
    if(is.null(names(ylist$Y))) names(ylist$Y) <- rownames(xlist$X)
    else if(is.null(rownames(xlist$X))) rownames(xlist$X) <- names(ylist$Y)

    return( list(ylist=ylist, xlist=xlist, algo=algo) )
}

## check the coordinate descent par
algochk <- function(tol=1e-5, qn=10^5, trust=1, verb=FALSE, store=FALSE){
  stopifnot(tol>0)
  return(list(tol=as.double(tol),
              qn=as.double(qn),
              trust=as.double(trust),
              verb=verb,
              store=store)) }

## check the response
glresponse <- function(response, fixe, family){

    if(family=="linear") famid <- 1
    else if(family=="poisson") famid <- 2
    else if(family=="binomial") famid <- 3
    else stop("Bad family argument")

    if(!is.vector(response)){
        if(is.matrix(response) && ncol(response)==1) response <- response[,1]
        else stop("gamlr only takes vector response arguments") }
    if(famid==3){
        if(is.factor(response)) response <- unclass(response)-1
        response <- (response>0) - (response <=0) }

    n <- length(response)
    if(is.null(fixe))
        fixe <- rep(0, n)
    else if(length(fixe)!=n)
        stop("fixed effects length must match response length.")

    ## add intercept
    Y <- response
    ysum <- sum(Y)

    return(list(Y=as.double(Y),
                ysum=as.double(ysum),
                fixe=as.double(fixe),
                famname=family,
                famid=as.integer(famid)))
}

## check the covariates
glcovars <- function(covars, scale, free=NULL){

    if(is.null(dim(covars)))
        covars <- matrix(covars, dimnames=list(names(covars),
                                 deparse(substitute(covars))))
    if(is.data.frame(covars))
        covars <- sapply(covars,as.numeric)

    ## standardize for mean 0 sample variance 1.
    if(scale){
        covarSD <- sdev(covars)*sqrt(1-1/nrow(covars))
        if(is.simple_triplet_matrix(covars)) covarMean <- rep(0,ncol(covars))
        else{ covars <- as.matrix(covars)
              covarMean <- colMeans(covars) }
        covars <- normalize(covars, m=covarMean, s=covarSD) }
    else{ covarMean <- covarSD <- NULL }

    ##  and check names
    if(is.null(colnames(covars)))
        colnames(covars) <- paste('v',1:ncol(covars),sep="")
    X <- as.simple_triplet_matrix(covars)
    X <- cbind(intercept=rep(1,nrow(covars)), X)

    ## format X with col major index ordering for passing to c
    o <- order(X$j)
    X$i <- X$i[o]
    X$j <- X$j[o]
    X$v <- as.double(X$v[o])
    xcnt <- as.integer(c(0,cumsum(col_sums(X!=0))))

    ## unpenalized var (always incl. intercept)
    free <- c(1,free+1)
    if(any(free>ncol(X) | free<1)) stop("bad `free' argument")

    return(list(X=X, xcnt=xcnt, free=free,
                scale=scale, covarMean=covarMean, covarSD=covarSD))
}

glchkpos <- function(arg, d=1){
    if(d>1 && length(arg)==1)
        arg <- rep(arg,d)
    if(any(arg<=0) || length(arg)!=d)
        stop(paste("bad argument: ", deparse(substitute(arg))))
    return(as.double(arg))
}
