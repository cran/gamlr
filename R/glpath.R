
## undocumented path estimation function
glpath <- function(xlist, ylist, algo, path)
{
    X <- xlist$X
    n <- nrow(xlist$X)
    p <- ncol(xlist$X)
    xr <- as.integer(X$i-1)
    xc <- as.integer(X$j-1)
    xcnt <- xlist$xcnt

    ## penalized variable indexing
    isfree <- sort(xlist$free)
    notfree <- (1:p)[-isfree]
    pfree <- length(isfree)
    pnotfree <- length(notfree)
    neverz <- notfree

############# build path & define stopping rule
    path <- chkpath(path,p=p,isfree=isfree, famname=ylist$famname)
    if(path$stoparg$rule=="BF")
      stopcnd <- function(){
        if(seg<(path$stoparg$val+1)) return(TRUE)
        bfd <- BF[(seg-path$stoparg$val):seg]
        if(any(is.nan(bfd))) return(FALSE)
        return(!all( diff(bfd) < 0 )) }
    else if(path$stoparg$rule=="maxvar")
        stopcnd <- function() return(nz < path$stoparg$val)
    else if(path$stoparg$rule=="maxpve")
        stopcnd <- function() return(pve[seg] < path$stoparg$val)
    else if(path$stoparg$rule=="minpen")
        stopcnd <- function() return(mu > path$stoparg$val)
    else if(path$stoparg$rule=="maxseg")
        stopcnd <- function() return(seg <= path$stoparg$val)
    D <- as.double(rep(algo$trust,p))

########## initial free-variable fit
    fit <- .C("R_gamlr",
               famid = ylist$famid,
               n = n,
               d = pfree,
               Y = ylist$Y,
               ysum = ylist$ysum,
               Nx = length(xlist$X[,xlist$free]$v),
               X = xlist$X[,xlist$free]$v,
               xr = as.integer(xlist$X[,xlist$free]$i-1),
               xc = as.integer(xlist$X[,xlist$free]$j-1),
               xcnt = as.integer(c(0,cumsum(col_sums(xlist$X[,xlist$free]!=0)))),
               B = as.double(rep(0,pfree)),
               E = ylist$fixe,
               mod = as.integer(rep(0,pfree)),
               par = double(0),
               rate = double(0),
               D = D[1:pfree],
               G = double(pfree),
               H = double(pfree),
               tol = algo$tol,
               qn = algo$qn,
               verb = as.integer(algo$verb > 1),
               PACKAGE="gamlr",
               dup=FALSE)


    B <- simple_triplet_matrix(i=xlist$free,
                               j=rep(1,pfree),
                               v=fit$B,
                               nrow=p, ncol=1)
    fit$B <- as.matrix(B)
    fit$G <- glgrad(ylist$famid, X, xcnt, fit$E, ylist$Y)
    fit$H <- glcurve(ylist$famid, X, xcnt, fit$E)

###  initial log marginal lhd and stats
    gradcut = abs(fit$G)
    gradcut[isfree] <- 0
    pve <- pve(ylist$Y, fit$E, fam=ylist$famname)
    active <- 0
    glselect <- .C("R_glselect",
                estlam = as.integer(path$estlam),
                famid = ylist$famid,
                n = n,
                Y = ylist$Y,
                eta = as.double(fit$E),
                p = as.integer(0),
                df = as.double(pfree),
                beta = 0,hess = 0,rate = 0,mu = 0,
                LPY = double(1),
                LLHD = double(1),
                PACKAGE="gamlr",
                dup=FALSE)
    lpy0 <- glselect$LPY
    llhd <- glselect$LLHD
    BF <- 0.0
    bfdrops <- 0
    dof <- pfree
    zerodf <- 0
    s2 <- disp(ylist$Y, fit$E,
               family=ylist$famname, df=pfree)

######### likelihood curvature calculations

    if(path$estlam){
      if(!is.null(path$fixpenvar)) lamvar <- path$fixpenvar
      else{
        ## calculate baseline "predictions"
        if(ylist$famname=="binomial") h <- 1/4
        else if(ylist$famname=="poisson") h <- ylist$Y
        else h <- 1
        ## fill in curvature values
        designcurve <- rep(0.0,p)
        if(path$ortho){
          if(algo$verb) cat("calculating orthogonal curvature ...\n")
          phis <- ## near correlations
            tcrossprod_simple_triplet_matrix(t(normalize(X[,notfree],m=0)))/n
          diag(phis) <- 0
          ## pls directions
          Z <- tcrossprod_simple_triplet_matrix(X[,notfree],phis)
          ## utility matrix
          Xhat <- X[,notfree]
          Xindex <- cbind(Xhat$i,Xhat$j)
          ## project
          Xhat$v <- Xhat$v*Z[Xindex]
          BZ <- col_sums(Xhat)/colSums(Z^2)
          Xhat$v <- Z[Xindex]*BZ[Xhat$j]
          ## curvature based on residuals
          Xhat <- X[,notfree]-Xhat
          designcurve[notfree] <- col_sums(h*Xhat^2)
          rm(Z,Xhat,BZ,phis) }
        else
          designcurve[notfree] <- col_sums(h*X[,notfree]^2)
        designcurve[notfree][designcurve[notfree]==0] <- 1e-8
        lamvar <- designcurve*path$penvar
        rm(designcurve)
      }
    } else lamvar <- 0


##############  main function

    ## initialize
    nz <- 0
    ord <- c()

    if(path$estlam){ # gamma lasso
      ## set the initial penalty
      absG <- abs(fit$G)[notfree]
      mu <- as.double(absG[which.max(absG)])
      rate <- as.double(rep(0,p))
      rate[notfree] <- as.double(mu/lamvar[notfree])
    } else{ # lasso
      rate <- as.double(path$rate)
      absG <- abs(fit$G/rate)[notfree]
      mu <- as.double(absG[which.max(absG)])
    }

    mupath <- mu
    seg <- 1
    if(algo$store){
      R <- list('1'=rate)
      G <- list('1'=fit$G)
      H <- list('1'=fit$H)
      E <- list('1'=fit$E) }

    if(algo$verb)
      cat(sprintf("t=1, pen=%g at null model.\n",round(mu,2)))

    ## outer path loop
    while(stopcnd()){

      ## update path segment penalty
      pnz <- length(neverz)
      if(!is.null(path$mugrid))
        mu <- max(path$mugrid[path$mugrid<mu])
      else if(path$fixstep) mu <- mupath[seg]*(1-path$step)
      else{
        if(pnz>1)
          mu <- as.double(sort(abs(fit$G)[neverz],partial=pnz-1)[pnz-1])
        else mu <- 0.01*mu
        if(mu >= mupath[seg]-0.001) mu <- 0.95*mupath[seg]
      }


      if(path$estlam)
        rate[notfree] <- as.double(mu/lamvar[notfree])

      ## inner coordinate descent
      fit <- .C("R_gamlr",
                famid = ylist$famid,
                n = n,
                d = p,
                Y = ylist$Y,
                ysum = ylist$ysum,
                Nx = length(X$v),
                X = X$v,
                xr = xr,
                xc = xc,
                xcnt = xcnt,
                B = as.double(fit$B),
                E = as.double(fit$E),
                mod = path$mod,
                par = c(mu,mu),
                rate = rate,
                D = D,
                G = as.double(fit$G),
                H = as.double(fit$H),
                tol = algo$tol,
                qn = algo$qn,
                verb = as.integer((algo$verb)>1),
                PACKAGE="gamlr",
                dup=FALSE)


      ## pve
      pvenew = pve(ylist$Y, fit$E, fam=ylist$famname)
      if(pvenew < pve[seg]){
        warning("stopped due to a drop in PVE")
        break  }
      else pve <- c(pve, pvenew)
      seg <- seg+1

      ## track variable entry
      isnz <- notfree[fit$B[notfree]!=0]
      nz <- length(isnz)
      active <- c(active, nz)
      newb <- setdiff(isnz,ord)
      ord <- c(ord,newb)
      neverz <- setdiff(notfree,ord)

      ########  model selection ##########

      ## get the degrees of freedom
      if(path$estlam){
        gradcut[neverz] <- abs(fit$G[neverz])
        Edf <- pgamma(gradcut, rate[notfree]*mu, rate[notfree])
        dfmu <- pfree + sum(Edf)
        zerodf <- c(zerodf, sum(Edf[!(notfree%in%isnz)]))
      } else dfmu <- pfree + 1.2*nz

      ## call to c for marginal likelihood
      glselect <-  .C("R_glselect",
                      estlam = as.integer(path$estlam),
                      famid = ylist$famid,
                      n = n,
                      Y = ylist$Y,
                      eta = as.double(fit$E),
                      p = nz,
                      df = as.double(dfmu),
                      beta = as.double(fit$B[isnz]),
                      hess = as.double(fit$H[isnz]),
                      rate = c(rate[isnz]),
                      mu = mu,
                      LPY = double(1),
                      LLHD = double(1),
                      PACKAGE="gamlr",
                      dup=FALSE)

      ## track the results
      if(!is.finite(glselect$LPY)) glselect$LPY <- NaN
      BF <- c(BF, glselect$LPY - lpy0)
      llhd <- c(llhd, glselect$LLHD)
      dof <- c(dof, dfmu)
      s2 <- c(s2, disp(ylist$Y, fit$E,
                       family=ylist$famname, df=dfmu))

      ## update fits
      B <- cbind(B,fit$B)
      mupath <- c(mupath, mu)
      if(algo$store){
        R[[paste(seg)]] <- rate
        G[[paste(seg)]] <- fit$G
        H[[paste(seg)]] <- fit$H
        E[[paste(seg)]] <- fit$E }

      if(algo$verb)
        cat(sprintf("t=%d, pen=%g, pve=%g%% with %d nonzero coefficients\n",
                    seg, round(mu,2),round(pve[seg]*100,2), nz))
    }

    ## names, formatting, and adjustments
    rownames(B) <- colnames(X)
    if(length(ord)>0) names(ord) <- colnames(X)[ord]
    names(mupath) <- names(llhd) <- names(BF) <- 1:seg
    s2[s2<0] <- NaN
    if(algo$store){
      E <- as.matrix(as.data.frame(E))
      G <- as.matrix(as.data.frame(G))
      H <- as.matrix(as.data.frame(H))
      R <- as.matrix(as.data.frame(R))
      rownames(G) <- rownames(H) <- colnames(X)
      rownames(R) <- colnames(X)
      rownames(E) <- rownames(X)
      G <- G[notfree,]
      H <- H[notfree,]
      R <- R[notfree,]
    } else{ E <- G <- H <- R <- NULL }

    ## construct the return object
    return( list(response=ylist$Y,  X=xlist$X[,-1], fixe=ylist$fixe,
                 family=ylist$famname, free=xlist$free[-1]-1, nobs=n,
                 intercept=as.matrix(B[1,]),
                 loadings=B[-1,], order=ord,
                 pve=pve, dispersion=s2, llhd=llhd,
                 df=dof, BF=BF, zerodf=zerodf, active=active,
                 penalty=mupath, lamvar=lamvar, path=path,
                 fitted.values=E, gradient=G, rate=R, curvature=H) )

  }


chkpath <- function(path, p, isfree, famname){

    ## partially orthogonalize the path?
    if(is.null(path$ortho)) path$ortho <- TRUE

    ## [undocumentd] use the guaranteed convex algo?
    if(is.null(path$stable)) path$stable <- 1
    if(!path$stable) path$stable <- -1

    ## penalty model and variance
    if(is.null(path$penvar) && !is.null(path$penvar)) path$penvar <- path$penvar
    if(is.null(path$penvar)) path$penvar <- 1

    if(path$penvar > 0) path$estlam <- TRUE
    else{
      path$estlam <- FALSE
      path$rate <- rep(1,p) }

    if(path$estlam) path$mod <- rep(1,p)*path$stable
    else path$mod <- rep(2,p)

    path$mod[isfree] <- 0
    path$mod <- as.integer(path$mod)

    ## step size
    if(is.null(path$step))
        path$fixstep <- FALSE
    else{
        path$fixstep <- TRUE
        stopifnot(all(c(path$step>0,path$step<1))) }

    ## [undocumented] fixed variance
    if(!is.null(path$fixpenvar)){
      if(!path$estlam) stop("no variance to fix under the lasso")
      path$fixpenvar <- glchkpos(path$fixpenvar,p)
      path$fixpenvar[isfree] <- 0 }

    ## [undocumented] adaptive lasso weights
    if(!is.null(path$adapt)){
      if(path$estlam) stop("no 'adaptive' gamma lasso")
      path$rate <- path$adapt
      path$adapt <- TRUE
      if(length(path$rate)==(p-1)) path$rate <- c(1,path$rate)
      path$rate <- glchkpos(path$rate,p)
    }

    ## check stopping conditions
    if(is.null(path$stoparg))
      path$stoparg <-list(rule="BF",val=10)
    if(is.character(path$stoparg))
      path$stoparg <- list(rule=path$stoparg)
    if(!is.list(path$stoparg)) path$stoparg <- as.list(path$stoparg)
    if(is.null(names(path$stoparg)) && length(path$stoparg)==2)
      names(path$stoparg) <- c("rule","val")

    if(!is.null(path$mugrid)){
      path$stoparg$rule <- "minpen"
      path$stoparg$val <- min(path$mugrid)
    } else if(path$stoparg$rule=="BF"){
        if(is.null(path$stoparg$val)) path$stoparg$val <- 10
    } else if(path$stoparg$rule=="maxvar"){
        if(is.null(path$stoparg$val)) path$stoparg$val <- p-length(isfree)
        if(path$stoparg$val > (p-length(isfree)))
          stop("bad 'stoparg$val' maxvar")
    } else if(path$stoparg$rule=="maxpve"){
        if(is.null(path$stoparg$val)) path$stoparg$val <- 0.9
        if(path$stoparg$val >= 1 || path$stoparg$val <= 0)
          stop("bad 'stoparg$val' maxpve")
    } else if(path$stoparg$rule=="minpen"){
        if(is.null(path$stoparg$val)) path$stoparg$val <- exp(1)
        if(path$stoparg$val <= 0)
          stop("bad 'stoparg$val' minpen")
    } else if(path$stoparg$rule=="maxseg"){
        if(is.null(path$stoparg$val)) path$stoparg$val <- 500
        if(path$stoparg$val <= 0)
          stop("bad 'stoparg$val' maxseg")
    } else stop("bad stopping rule specification")

    return(path)
}

