\name{cvselect}
\alias{cvselect}
\alias{cvplot}
\title{ Cross Validation for gamlr
}
\description{Cross Validated gamma lasso penalty selection
  }
\usage{
cvselect(fit, nfold=10, npen=20, verb=TRUE)
cvplot(fit, ...)
}
\arguments{
  \item{fit}{ A fitted gamlr object.}
  \item{nfold}{ The number of cross validation folds. }
  \item{npen}{ The size of the penalty grid to consider. }
  \item{verb}{ Whether to print and plot progress. }
  \item{...}{ Extra plot arguments. }
}
\details{ Performs \code{nfold} out-of-sample prediction exercises 
on \code{n/nfold} sized test sets, and uses resulting PVE to choose the penalty size
(i.e. the size of the expected L1 regularization weight).  See \code{gamlr} for details. \code{cvplot} can be used to boxplot the results. }
\value{
  The returned \code{gamlr} object is the same as \code{fit}, but augmented with 
  \item{cvpve}{Out-of-sample proportion of variance explained.}
  \item{cvindex}{The index of the CV optimal penalty (within \code{fit$penalty}).}
  \item{cvgrid}{The grid of actual cross validated penalties.}
}
\references{
  Taddy 2013, "The Gamma-Lasso".
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{

n <- 100
p <- 100

xvar <- matrix(ncol=p,nrow=p)
for(i in 1:p) for(j in i:p) xvar[i,j] <- 0.5^{abs(i-j)}
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)
beta <- matrix( (-1)^(1:p)*exp(-(1:p)/10) )
mu = x\%*\%beta
y <- rbinom(n, prob=1/(1+exp(-2*mu)), size=1)

fit <- cvselect(gamlr(x, y, verb=TRUE, family="binomial"), npen=5, nfold=5)

## plot some model selection criterion	

par(mfrow=c(1,3))
## again, the CV plot
cvplot(fit)
## BF
plot(-log(fit$penalty),fit$BF,ylab="log BF", pch=21,bg="darkorange")
## adj BF
points(-log(fit$penalty),fit$BF-fit$zerodf,pch=21,bg="dodgerblue")
## fitted values and all crit
plot(fit, select="CV", col=8)
abline(v=-log(fit$penalty)[which.min(AIC(fit))], col="green", lty=2) # AIC
abline(v=-log(fit$penalty)[which.min(AIC(fit, k=log(fit$n)))], col="red", lty=2) # BIC
abline(v=-log(fit$penalty)[which.max(fit$BF)], col="darkorange", lty=2) # BF
abline(v=-log(fit$penalty)[which.max(fit$BF-fit$zerodf)], col="dodgerblue", lty=2) # adjBF
 }
\seealso{gamlr}