## calculate PVE
pve <- function(y, f, family="linear"){
  if(family=="binomial")
    if(is.factor(y)) y <- as.numeric(y)>1
    else y <- as.numeric(y>0)
   
  vy <- var(y)
  if(vy == 0) return(0)
  if(family=="linear")
    return( 1 - var(y-f)/vy)
  else if(family=="poisson")
    return( 1 - var(y-exp(f))/vy )
  else if(family=="binomial")
    return( 1 - var( (y>0) - 1/(1+exp(-f)) )/vy ) 
  else stop("bad gamlr object in glpve")
}

## calculate dispersion
disp <- function(y, f, family="linear", df=1){
  if(family=="binomial")
    if(is.factor(y)) y <- as.numeric(y)>1
    else y <- as.numeric(y>0)
  
  adj <- 1-df/length(y)
  if(family=="linear")
    return( mean( (y-f)^2 )/adj )
  else if(family=="poisson")
    return( mean( (y-exp(f))^2/exp(f) )/adj )
  else if(family=="binomial"){
    q <- 1/(1+exp(-f))
    return( mean( (y-q)^2/(q*(1-q)) )/adj ) }
  else stop("bad gamlr object in glpve")
}
