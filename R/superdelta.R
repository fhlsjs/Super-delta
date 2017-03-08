is.wholenumber<-function(x, tol = .Machine$double.eps^0.5){
  res<-numeric(length(x))
  for (i in 1:length(x)){
    res[i]<-ifelse(!is.numeric(x[i]), 0, abs(x[i]-round(x[i]))<tol)
  }
  return(as.logical(res))
}
### turn a vector of t-statistics into two sided pvalues
#Define a function to handle 2-sided and/or 1-sided test#
.t2p <- function(x, df, side=c("two.sided", "less", "greater")){
    if (side=="two.sided"){
        return(2-2*pt(abs(x),df))
    } else if (side=="less"){
        return(pt(x, df))
    } else if (side=="greater"){
        return(1-pt(x, df))
    } else {
        stop(paste("side must be one of the following three options: two.sided (default), less, or greater."))
    }
}

## Median; fold; trim; median
#Trim the outliers on both sides by MAD (Median Absolute Distance)#
#YL made this note on Jul 13, 2015#
.mftm <- function(tvec, trim=0.2){
    tvec.no.na <- na.omit(tvec)
    ## tvec.folded is firsrt centered around median and then folded
    ## around zero.
    tvec.folded <- abs(tvec.no.na-median(tvec.no.na))
    ## Remove trim amount of extreme observations from tvec.folded
    t.est <- median(tvec.no.na[tvec.folded<=quantile(tvec.folded, 1-trim)])
    ##t.est is the median after trimming
    return(t.est)
}

## cluster tvec by mclust
.tclust1 <- function(tvec, ...) {
    ## get rid of the NAs first
    tvec <- na.omit(tvec)
    mm.i <- summary(mclustBIC(tvec, G=3, modelName="V", ...), tvec)
    NC1 <- sum(mm.i$classification==1)
    NC2 <- sum(mm.i$classification==2)
    NC3 <- sum(mm.i$classification==3)
    largest.cluster <- which.max(c(NC1, NC2, NC3))
    t.est <- as.double(mm.i$parameters$mean[largest.cluster])
    return(t.est)
}

.tclust2 <- function(tvec, trim=0.2, ...){
    ## mclust is not that robust, so use the following wrapper as a backup.
    .substitute <- function(x) {
        warning("Mclust does not converge, use the robust method instead.")
        .mftm(tvec, trim=trim)
    }
    t.est <- tryCatch(.tclust1(tvec, ...), warning=.substitute, error=.substitute)
    #If no condition is signaled when evaluating expr#
    #then tryCatch returns the value of the expression#
    return(t.est)
}

## this function takes tmat, which is a matrix of
## t-statistics. tmat[i,j] is the tstat associated with the ith gene
## normalized by the jth gene.  It estimates one best t-stat according
## to different methods.
.est.t <- function(tmat, method=c("robust", "mean", "median", "mclust"), trim=0.2, ...){
    method=match.arg(method)
    if (method=="robust"){
        t.est <- apply(tmat, 1, .mftm, trim=trim)
    } else if (method=="median"){
        t.est <- apply(tmat, 1, median, na.rm=TRUE)
    } else if (method=="mean"){
        t.est <- rowMeans(tmat, na.rm=TRUE)
    } else if (method=="mclust"){
        t.est <- apply(tmat, 1, .tclust2)
    } else {
        stop(paste("Method", method, "is not implemented!"))
    }
    ## 7/30/2011.  Try t-stat with discounted variance. This discount
    ## is useful because the variance of delta is twice as much as the
    ## variance of the original expression.
    return(sqrt(2)*t.est)
}

superdelta <- function(X, classlabel, test="t", side=c("two.sided", "less", "greater"), methods="robust", trim=0.2, baseline="auto", ...)
{
    test <- match.arg(test)
    side <- match.arg(side)
    ngenes <- dim(X)[1]  #nrow#
    nslides <- dim(X)[2]  #ncol#
    teststat <- matrix(0, nrow=ngenes, ncol=length(methods))
    ind <- matrix(0, nrow=ngenes, ncol=length(methods))
    colnames(teststat) <- colnames(ind) <- methods
    rownames(teststat) <- rownames(ind) <- rownames(X)
    ## Let's turn the classlabel into a binary vector that is required
    ## by two_sample_tstat() in C.
    classlabel <- as.integer(as.factor(classlabel)) - 1
    Xa <- X[, classlabel==0]; Xb <- X[, classlabel==1]
    n1 <- ncol(Xa); n2 <- ncol(Xb)
    ## gene-order is calculated via pooled sample variance/STD
    Va <- apply(Xa, 1, var); Vb <- apply(Xb, 1, var)
    DFnum <- (Va/n1 + Vb/n2)^2
    DFdenom <- Va^2/n1^2/(n1-1) + Vb^2/n2^2/(n2-1)
    DF <- round(DFnum/DFdenom)
    ## Now deal with the baseline genes.  First, if baseline is a
    ## vector of genes/names, we assume it is a pre-specified list of
    ## baseline genes.
    if (length(baseline)>1){
        ##assume this is a vector of pre-selected baseline genes (heuristic)
        if (!all(is.wholenumber(baseline))){
            ## try to translate the list of gene symbols into indices
            baseline.genes <- which(rownames(X) %in% baseline)
        } else {
            ## assume that this is a list of integer-valued indices
            baseline.genes <- baseline
        }
    } else if (baseline=="auto") {
        ##if ngenes is greater than 1000, we only use 1000
        ##genes. Otherwise, we use all genes.
        if (ngenes>1000){
            baseline.genes <- sample(ngenes, 1000)
        } else {
            ##use all genes
            baseline.genes <- 1:ngenes
        }
    } else if (baseline=="all"){
        baseline.genes <- 1:ngenes
    } else {
        stop("Variable baseline must be either: 1. A vector of pre-specified baseline genes (length>1). 2. String value 'auto' or 'all'. ")
    }
    ## We have to reshape the data matrix in order to use C function
    ## superdelta_tstats().
    glist <- c(baseline.genes, setdiff(1:ngenes, baseline.genes))
    ## Y must be passed to C program as a numerical VECTOR!
    Y <- as.double(as.matrix(X)[glist, ])
    ## The real calculation starts
    Tmat <- rep(0, ngenes*length(baseline.genes))
    res<-.C("superdelta_tstats",as.double(Y), as.integer(ngenes),
            as.integer(nslides),as.integer(length(baseline.genes)),
            as.integer(classlabel),
            Tmat=as.double(Tmat), PACKAGE="superdelta")$Tmat
    ## Now make a matrix out of the res and re-organize it in the
    ## original ordering.
    tmat <- matrix(res, nrow=ngenes)[order(glist),]
    
    ## Now estimate the centers
    for (mm in methods){
        teststat[,mm] <- .est.t(tmat, method=mm, trim=trim, ...)
    }
    ## Now find the indices of matching pair gene of each gene
    for (mm in methods){
      for (i in 1:ngenes){
        ind[i,mm] <- baseline.genes[which.min(abs(sqrt(2)*tmat[i,]-teststat[i,mm]))]
      }
    }
    ## two sided or one sided test; different test statistics, etc ##
    if (test=="t"){
        rawp=.t2p(teststat, df=DF, side=side)
    } else {
        stop(paste("Test statistic", test, "is not implemented yet!"))
    }

    ## two representations ##
    if (ncol(teststat)==1){
        rr <- data.frame(teststat=teststat, rawp=rawp, ind=ind)
        #rr <- list(teststat=teststat, rawp=rawp, ind=ind, tmat=tmat)
        colnames(rr) <- c("teststat", "rawp", "pairing gene index")
        return(rr)
    } else {
        return(list(teststat=teststat, rawp=rawp, ind=ind))
    }
}