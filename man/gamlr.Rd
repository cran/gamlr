\name{gamlr}
\alias{gamlr}
\alias{print.gamlr}
\alias{summary.gamlr}
\alias{predict.gamlr}
\alias{plot.gamlr}
\alias{coef.gamlr}
\alias{AIC.gamlr}
\title{ Gamma-Lasso regression }
\description{
MAP estimation of regression models.
  }
\usage{
gamlr(covars, response, family="linear",
      free=NULL, scale=TRUE, fixe=NULL, store=TRUE,
      penvar=1, ortho=TRUE, step=NULL,
      stoparg=list(rule="BF",val=10),
      cdpar=list(tol=1e-5), verb=FALSE, ...)
}
\arguments{
  \item{covars}{
    A \code{matrix} or \code{simple_triplet_matrix} of
    \code{ncol(covars)} covariate values
    for each of the \code{length(response)} observations.
    This does not include the intercept, which is ALWAYS added in the design matrix.
  }
  \item{response}{ A vector of response values. }
  \item{family}{ Model type; either "linear", "poisson", or
	"binomial".}
  \item{free}{ Columns of \code{covars} that are to be fit without
	penalty.  This over-rides any other relevant arguments.}
  \item{scale}{
    Whether or not to scale the covariates. Default is \code{TRUE}.
    If \code{covars} is a \code{matrix}, this will scale the inputs to
    have mean zero and standard deviation of one.
    If \code{covars} is a \code{simple_triplet_matrix}, we assume that
    you want to stay in sparse format; hence the inputs are scaled to
    have sd = 1 but left unshifted.
  }
  \item{fixe}{ Possible fixed effects to be added to the linear equation
	  of your model.}
  \item{store}{ Whether to store fitted, gradient, and curvature
    values.}
  \item{penvar}{ The "prior variance" scaling on your L1 penalty. The default
  				 is 1, and higher values correspond to a MORE concave penalty.
				 \code{penvar=0} is the lasso.

  				 Actual lambda variance for each coefficient is \code{penvar}
				 times a measure of negative log likelihood curvature for that variable
				 (see \code{ortho} below and details). }
  \item{ortho}{ If TRUE, we calculate likelihood curvature (which is used to set
  				penalty variance; see \code{penvar}) based on
  				residuals for each covariate regressed onto all others.
				\code{ortho=FALSE} just uses raw covariate values.

				Note that if this is too slow you can just set a smaller \code{penvar}
				and achieve similar path continuity results. }
  \item{step}{ Possible fixed size of segment steps as a
			   proportion of previous penalty (must be between zero and one).  If
			   \code{step=NULL} the path jumps according to gradient information. }
  \item{stoparg}{ A list with items "rule" and "val" giving the rule and threshold
  				  for terminating the regularization path. Possible settings include:
				  \code{list("BF",val)}: stop after
				  \code{val} drops in estimated Bayes Factor,
				  \code{list("maxpve",val)}: stop at \code{val} PVE,
				  \code{list("maxvar",val)}: stop when the model
				  includes \code{val} nonzero penalized coefficients,
				  \code{list("minpen",val)}: stop at \code{val}	penalty expectation.
				  If only the rule is set then defaults are used for val.

				  You can also specify \code{mugrid} to give a
				  fixed grid of penalties.  This over-rides \code{stoparg}.

				  Note that the path \emph{always} stops if there is a drop in PVE.}
  \item{cdpar}{ List of coordinate descent algorithm parameters,
  				including \code{tol}: the optimization convergence tolerance for
				absolute change over an update to the full parameter set.  }
  \item{verb}{ Control for print-statement output.  0 prints nothing,
  			   higher integers print more.  }
  \item{...}{ Extra undocumented arguments. }
}
\details{ Finds the posterior modes along a regularization path
		   for regression parameters using coordinate descent.

  \emph{Coefficient penalization} is based upon the precision parameters
  \eqn{\lambda} of independent Laplace priors on each penalized
  	   regression coefficient.  Here, the Laplace density is \eqn{p(z) =
       (\lambda/2)exp[-\lambda|z|]}, with variance \eqn{2/\lambda}. Via the
       \code{penvar} argument, this precision is either fixed (\code{penvar=0}), which
       corresponds to the L1 penalty \eqn{\lambda|z|}, or it is assigned a
       \eqn{Gamma(s, r)} prior, with expectation \eqn{s/r} and rate
       \eqn{r}, which corresponds to the log penalty \eqn{s*log[r + |z|]}.
	   In the latter case, coefficient penalties are specified such that
	   \eqn{var(\lambda) = s/r^2 =} \code{penvar*curve}, where curve is
	   a measure of negative log likelihood curvature for that covariate.


  \emph{S3 Method functions} (\code{plot,coef,summary,print}) take standard
  arguments. In addition:

  \code{plot},\code{coef}, and \code{predict} take the argument
  \code{select} to return results for a single penalty selected through
  some model evaluation criterion.  Options for \code{select} are "BF" for
  the active-set Bayes Factor, "adjBF" for the BF adjusted to account for the number
  of coefficients set to zero, "AIC","BIC", and "CV".  Note that for "CV" you need
  to have first run \code{cvselect}.  \code{select=TRUE} defaults to "adjBF".

  The \code{plot.gamlr} function takes the argument
  \code{against} to indicate if you'd like to show
  solution paths against either "pen", "logpen", or "pve".

  The \code{summary.glm} function takes argument \code{top}
  to set how many "most important" coefficients.

  The \code{coef.gamlr} function takes \code{origscale=TRUE} which
  indicates whether coefficients should be returned on the original
  covariate scale (as opposed to on the scale of normalized covariates).

  The \code{predict.gamlr} function returns new values on the  natural linear
  scale  \eqn{X'\beta} (same as \code{fitted} above).
}
\value{
  The returned \code{gamlr} object list including some of the
  input variables (possibly cleaned or updated), as well as
  \item{X}{A simple triplet matrix of the (scaled) input covariates.}
  \item{nobs}{ The number of observations.}
  \item{intercept}{The fitted intercept.}
  \item{loadings}{Simple triplet matrix of regression coefficients.}
  \item{pve}{The in-sample proportion of variance explained.}
  \item{dispersion}{Fitted residual dispersion.}
  \item{llhd}{Fitted log likelihood.}
  \item{df}{Approximate degrees of freedom.}
  \item{BF}{log Bayes factors with respect to the model with only
	unpenalized coefficients.}
  \item{zerodf}{Approximate degrees of freedom for the in-active set.  The
  							adjusted Bayes factor is then \code{BF - zerodf}. }
  \item{active}{Number of nonzero penalized covariates.}
  \item{penalty}{Path of either the expected or fixed L1 penalty.}
  \item{lamvar}{Variance of the L1 penalties for each coefficient.}
  \item{covarMean}{If \code{scale}, the amount covariates were
    shifted (original means for matrix \code{covars}, 0 for sparse stm
    \code{covars}).  Otherwise empty.}
  \item{covarSD}{If \code{scale}, the original
    covariate standard deviations.  Otherwise empty.}
  \item{rate}{If \code{store}, By-coefficient rate parameters along the path.}
  \item{fitted}{If \code{store}, Fitted natural parameters \eqn{X'\beta}.}
  \item{grad}{If \code{store}, negative llhd radient at solutions. }
  \item{curve}{If \code{store}, negative llhd curvature at solutions. }
}
\references{
  Taddy 2013, "The Gamma-Lasso".
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{


### a low-D test (super multi-collinear)

n <- 1000
p <- 3
xvar <- matrix(0.9, nrow=p,ncol=p)
diag(xvar) <- 1
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)
y <- 4 + 3*x[,1] + -1*x[,2] + rnorm(n)

## stoparg="maxvar" defaults to list(rule="maxvar",val=ncol(covars))
fitl <- gamlr(x, y, step=.05, penvar=0, stoparg="maxvar") # lasso
fitgl <- gamlr(x, y, step=.05, stoparg="maxvar") # default log penalty
fitglno <- gamlr(x, y, step=.05, ortho=FALSE, stoparg="maxvar") # without orthogonalized curvature

par(mfrow=c(1,3))
ylim = range(c(fitgl$load$v,fitglno$load$v))
plot(fitl, ylim=ylim)
plot(fitgl, ylim=ylim)
plot(fitglno, ylim=ylim)


### a larger more adversarial example

n <- 100
p <- 200

xvar <- matrix(ncol=p,nrow=p)
for(i in 1:p) for(j in i:p) xvar[i,j] <- 0.5^{abs(i-j)}
x <- matrix(rnorm(p*n), nrow=n)\%*\%chol(xvar)

beta <- matrix( (-1)^(1:p)*exp(-(1:p)/10) )
mu = x\%*\%beta
y = mu + rnorm(n,sd=sd(c(mu))/2)

fitl <- gamlr(x, y, penvar=0)
fitgl <- gamlr(x, y)
fitglno <- gamlr(x, y, ortho=FALSE)

par(mfrow=c(1,3))
ylim = range(c(fitgl$load$v,fitglno$load$v))
plot(fitl, ylim=ylim, select="BF")
plot(fitgl, ylim=ylim, select="BF")
plot(fitglno, ylim=ylim, select="BF")

 }
\seealso{cvselect,hockey}
