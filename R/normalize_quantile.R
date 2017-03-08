## Self-written quantile normalization
normalize.quantile<-function(X,na.rm=TRUE){
  ngene<-nrow(X)
  nsamp<-ncol(X)
  if (is.numeric(X)==0) stop("Error: Input matrix must be numeric!")
  Xquant<-Xsort<-Xmean<-matrix(0,ngene,nsamp)
  for (j in 1:nsamp){
    Xsort[,j]<-sort(X[,j])
  }
  for (i in 1:ngene){
    for (j in 1:nsamp){
    Xmean[i,j]<-mean(Xsort[i,])
    }
  }
  for (j in 1:nsamp){
    Xquant[,j]<-Xmean[,j][rank(X[,j],ties.method = "random")]
  }
  return(Xquant)
}
## Test for this function
## x<-matrix(c(8,7,3,1,9,15,2,6,5,13,9,7,5,2,6,13,15,8,9,11),nrow=5,ncol=4)
## x
## normalize.quantile(x)
