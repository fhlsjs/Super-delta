##########################################
###  Write a function to filter genes  ###
###  with very small variability       ###
##########################################
## Take IQR as a robust measurement of variability
## Default trimming proportion is 20%
## We may want to work on different proportions, 10%, 15%, 20%, 25%, and 30%
filt_IQR <- function(X, prop=0.2, na.rm=TRUE)
  {
    rownum <- nrow(as.matrix(X))
    filt <- round(rownum*prop)
    X_filt <- as.matrix(X)[-order(apply(as.matrix(X),1,IQR))[1:filt],]
    return(X_filt)
}

## Test this function
## X<-matrix(rnorm(5380,0,1),nrow=538,ncol=10)
## Y<-filt_IQR(X,prop=0.18)
#538*0.18=96.84~97
#Don't use too "good" numbers
## View(X)
## View(Y)
## dim(Y)
## order(apply(as.matrix(X),1,IQR))[1:97]