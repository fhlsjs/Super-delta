normalize.rank <- function(X, na.rm=TRUE){
  ngenes <- dim(X)[1]
  nsamp <- dim(X)[2]
  Xrank <- matrix(0,nrow=ngenes, ncol=nsamp)
  for (i in 1:nsamp){
    Xrank[,i] <- rank(X[,i])/ngenes
  }
  return(Xrank)
}
## Test for this function
## x<-matrix(rnorm(100,0,1),10,10)
## x
## normalize.rank(x)