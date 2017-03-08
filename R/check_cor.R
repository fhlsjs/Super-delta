########################################
### This function is for randomly    ###
### selecting genes from a data set  ###
### and check correlations,          ###
### mainly for verifying raw data    ###
########################################
check_cor <- function(X, n=1000, na.rm=TRUE){
  if (nrow(as.matrix(X))<=n){
    COR <- cor(t(X))
    COR1 <- COR[1,2:nrow(X)]
  }
  else{
    X_ran <- as.matrix(X)[sample(nrow(as.matrix(X)),n),]
    COR <- cor(t(X_ran))
    COR1 <- COR[1,2:n]
  }
  hist(COR1, main="Histogram of correlations between \n 1000 randomly selected genes",
       xlim=c(-1,1), xlab="Correlations",ylab="Frequency")
}
## Test for the check_cor function
## check_cor(Total_Raw,n=100)
## Most of the correlations are huge (>0.8)
## Based on Qiu et. al. (2005), these data are really raw data
## This indicates strong requirement of normalization procedures