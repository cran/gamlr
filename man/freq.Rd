\name{freq}
\alias{freq}
\alias{tfidf}
\title{
Frequency Matrix Conversion
}
\description{
Convert a count matrix to the corresponding frequency or term frequency
- inverse document frequency matrix.}
\usage{
freq(x, byrow=TRUE)
tfidf(x, freq=FALSE, offset=1)
}
\arguments{
  \item{x}{A matrix or \code{simple_triplet_matrix} with count entries. }
  \item{byrow}{An indicator for whether you have observation-rows and
	category-columns, or vice versa.}
  \item{freq}{An indicator for whether \code{x} is already a frequency
	matrix.}
  \item{offset}{An offset added to the total number of documents. }
}
\value{
 \code{freq} returns a matrix with row (\code{byrow=TRUE}) or column
 (\code{byrow=FALSE}) sums of one.  \code{tfidf} returns a matrix
 with entries \eqn{f_{ij} \log[ (n+o)/d_j ]}, where \eqn{f_{ij}} is term-j frequency in document-i, 
  \eqn{d_j} is the number of documents containing term-j, and \eqn{o} is
  your offset.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
F <- freq( t(rmultinom(10, 20, c(1/2,1/4,1/8,1/8))) )
tfidf(F, freq=TRUE)
 }
