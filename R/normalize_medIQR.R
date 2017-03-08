## X is raw data, fin_prod is the normalized data
## set.seed(1234)
## iqr_norm_yuhang<-function(ngene,narray,mug,sigmag,mun=0,sigman){
##     X<-matrix(data=NA,ngene,narray)
##     Y<-matrix(data=NA,ngene,narray)
##     Z<-matrix(data=NA,ngene,narray)
##     fin_prod<-matrix(data=NA,ngene,narray)

##     ##narray = number of arrays
##     ##ngene = number of genes
##     ##mug = mean vector of gene expressions, to be estimated from real data
##     ##mun = mean of noise, noise is array-wise, WLOG assume mun=0
##     ##sigmag = var-cov matrix of gene expressions, to be estimated from real data
##     ##sigman = standard deviation of noise, to be estimated from real data
##     ##noise is like random effect in ANOVA model
##     for (i in 1:narray){
##         X[,i]<-mvrnorm(n=1,mug,sigmag)+rnorm(1,0,sigman)
##     }

##     ##Ignore boxplots for large number of arrays
##     ##Parallel boxplot for raw data
##     ##boxplot(as.data.frame(X),xlab="Channel (Array)",ylab="Log Signal",axes=F)
##     ##axis(1,labels=1:narray,at=1:narray)
##     ##axis(2)
##     ##box()
##                                         #
##     med<-apply(X,2,median)
##     MED<-median(med)
##     iqr<-apply(X,2,IQR)
##     IQRR<-median(iqr)
##     for (i in 1:narray){
##         Y[,i]<-X[,i]-med[i]
##     }
##     for (i in 1:narray){
##         Z[,i]<-Y[,i]*IQRR/iqr[i]
##     }
##     for (i in 1:narray){
##         fin_prod[,i]<-Z[,i]+MED
##     }

##     ##Ignore boxplots for large number of arrays
##     ##Parallel boxplot for normalized data
##     ##boxplot(as.data.frame(fin_prod),xlab="Channel (Array)",
##     ##ylab="Normalized Data",axes=F)
##     ##axis(1,labels=1:narray,at=1:narray)
##     ##axis(2)
##     ##box()
##     return(fin_prod)
## }

normalize.medIQR <- function(X, na.rm=TRUE, TOL=.001){
    med<-apply(X,2,median, na.rm=na.rm)
    iqr<-apply(X,2,IQR, na.rm=na.rm)
    for (i in 1:ncol(X)){
        ## we only need to generate warnings when iqr are too small
        ## (and may result in numerical instability.
        if (iqr[i] <= TOL)
            warning(paste0("Column ", i, ": IQR is too small (<", TOL, "). Median/IQR normalization may be numerically unstable. Removal of this array is recommended."))
    }
    ## Y is the centered and scaled matrix. Needs to revise the code
    ## below so that it generates warnings for slides with (almost)
    ## constant expressions (IQR==0 or very small).
    Y <- sweep(sweep(X, 2, med, "-"), 2, iqr, "/")
    ## Xnorm is the normalized X
    Xnorm <- Y*median(iqr) + median(med)
    return(Xnorm)
}
