## correlation for slam simple_triplet_matrix and regular matrix
corr <- function(x, y){
  if(!inherits(x, "simple_triplet_matrix")){ return(cor(x,y) ) }

  n <- nrow(x)
  v <- t(normalize(y))

  r <- tcrossprod_simple_triplet_matrix(t(x)/sdev(x), v)/(nrow(x)-1)
  dimnames(r) <- list(dimnames(x)[[2]], dimnames(y)[[2]])
  return( r ) }

## column standard deviation for a simple_triplet_matrix
sdev <- function(x){
  if(!inherits(x, "simple_triplet_matrix")){ return(apply(as.matrix(x),2,sd)) }
  n <- nrow(x)
  sqrt(col_sums(x^2)/(n-1) - col_sums(x)^2/(n^2 - n))
  return( sqrt(col_sums(x^2)/(n-1) - col_sums(x)^2/(n^2 - n)) ) }

##  normalizing design matrices
normalize <- function(x, m=NULL, s=sdev(x)){

  if(!is.null(ncol(x)))
    if(length(s)!=ncol(x)) stop("length(s)!=ncol(x)")
  s[s==0] <- 1

  if(is.simple_triplet_matrix(x)){
    x$v <- x$v/s[x$j]
    return(x) }
  
  x <- as.matrix(x)
  if(is.null(m)) m <- col_means(x)
  return( t((t(x) - m)/s) )
}

## converting count to frequency matrix
freq <- function(x, byrow=TRUE){
    if(byrow){ s <- row_sums(x)
               s[s==0] <- 1
               return( x/s ) }
    else{
      s <- col_sums(x)
      s[s==0] <- 1
      return(t(t(x)/s)) }
}

## converting a count/freq matrix to tf-idf
tfidf <- function(x, freq=FALSE, offset=1){

  if(!freq){tf <- freq(x)}
  else{tf <- x}

  idf <- log( nrow(x)+offset ) - log(col_sums(x>0))
  x <- t( t(tf) * idf )

  return( x ) }
