t_test_genes<-function(X,classlabel,na.rm=TRUE){
  ngenes <- nrow(X)
  CL <- as.numeric(as.factor(classlabel)) - 1
  if (is.numeric(X)==0 | is.matrix(X)==0)
    stop("Error: Input data must be a numeric matrix!")
  if (length(unique(CL))!=2)
    stop("Error: Input data must have 2 phenotypic groups!")
  Xa <- X[,CL==0]; Xb <- X[,CL==1]
  teststat <- rawp <- numeric(ngenes)
  for (i in 1:ngenes){
    myttest <- t.test(Xa[i,], Xb[i,], alternative = "two.sided", mu = 0, 
                      paired = FALSE, var.equal = FALSE, conf.level = 0.95)
    teststat[i] <- myttest$statistic
    rawp[i] <- myttest$p.value
  }
  rr <- data.frame(teststat=teststat, rawp=rawp)
  colnames(rr) <- c("teststat", "rawp")
  return(rr)
  #return(rawp)
}

## rawp<-numeric(17826)
## for (i in 1:17826){
  ##rawp<-t.test_genes(normalize.rank(TA_Raw_filt_pCR),normalize.rank(TA_Raw_filt_RD))
  ##sum(p.adjust(rawp,method="BH",n=length(rawp))<=0.05)
## }