\name{normalize}
\alias{normalize}
\title{
Normalize
}
\description{
Normalize matrix columns. }
\usage{
normalize(x, m=NULL, s=sdev(x))
}
\arguments{
  \item{x}{
    A matrix.
}
\item{m}{Optional column shifts.}
\item{s}{Column scaling factor (we divide by \code{s}).}
}
\value{
 If \code{x} is a matrix, the default returns mean-zero and variance-one columns.
 If \code{x} is a \code{simple_triplet_matrix}, columns are scaled by \code{sdev(x)} but left unshifted; the function returns a \code{simple_triplet_matrix}.   If shift and scale are specified AND \code{x} is a \code{matrix}, columns are shifted by \code{-m}  and divided by \code{s}.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
x <- matrix(1:9, ncol=3)
normalize( x )
as.matrix(normalize(as.simple_triplet_matrix(x)))
 }
