\name{pve}
\alias{pve}
\alias{disp}
\title{Proportion of Variance Explained and model dispersion
}
\description{Proportion of Variance Explained, and estimated dispersion parameters, 
for arbitrary fit and response in linear, binomial (logistic) and poisson (log) linear models. }
\usage{
pve(y, f, family="linear")	
disp(y, f, family="linear", df=1)	
}
\arguments{
  \item{y}{ A vector of response values.}
  \item{f}{ A vector of linear model fits (i.e. on logit scale for "binomial", on log scale for "poisson").}
  \item{family}{Regression family name.}
  \item{df}{For dispersion, the degrees of freedom used in estimating \code{y}.}
}
\details{ \code{pve} returns \eqn{1 - var(y-E[y|f])/var(y)}, 
and \code{disp} returns \eqn{[var(y-E[y|f])/var(y|f)]*[n/(n-df)]}.  
See the functions for formulas. }
\value{
  PVE on a [0,1] scale or dispersion. }
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
x <- rnorm(1000)
y <- x + rnorm(1000,sd=1/2)
f <- lm(y~x)$fitted
pve(y,f)
disp(y,f)
 }
\seealso{gamlr}