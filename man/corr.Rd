\name{corr}
\alias{corr}
\alias{sdev}
\title{Sparse Matrix Correlation}
\description{
Correlation calculations for a \code{simple_triplet_matrix} and a \code{matrix}.}
\usage{
corr(x, y)
sdev(x)
}
\arguments{
  \item{x}{A \code{simple_triplet_matrix} (or a \code{matrix}, in which case the 
  function returns \code{cor(x,y)}). }
  \item{y}{A \code{matrix} with \code{nrow(y)=nrow(x)}. }
}
\value{
 An \code{ncol(x)} by \code{ncol(y)} matrix containing correlation
 between x and y, or just the vector of column standard
 deviations for \code{sdev}.
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\examples{
data(hockey)

## some good players
r <- corr(Xplayer, Ygoals)
sort(r[,1], decreasing=TRUE)[1:20]

## players with lots of ice-time
s <- sdev(Xplayer)
sort(s, decreasing=TRUE)[1:20]
}
