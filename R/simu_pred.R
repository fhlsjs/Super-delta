simu_pred <- function(NGENES=10000, n1=50, n2=50, bl=NULL, fc=2^rnorm(NGENES,0,0.5), SIM=c("SIM1", "SIM2", "SIM3"), TopGenes=50, eta=1.5, sigma=1){
  if (length(fc)!=NGENES){
    stop("Input signal length must match total number of genes!")
  }
  if (!is.null(bl) & length(bl)!=NGENES){
    stop("Input baseline mean expression length must match total number of genes!")
  }
  ## bl is user input baseline mean log2 expression of one group
  ## fc is fold-change on raw scale of all input genes, Tusher et al. (2001)
  ## m is number of genes for SIM1
  ## n1, n2 are the two group sample sizes
  ## eta and sigma are model parameters, default signal-to-noise ratio 1.5
  ## eff_size is the vector of true signal before log2 transformation
  m <- ifelse(SIM=="SIM1", NGENES, round(NGENES/2))
  #print(NGENES); print(m)
  ind <- order(abs(log2(fc)), decreasing=T)[1:min(1000,round(0.4*m))]
  eff_size <- fc[ind]
  eff_size1 <- eff_size[eff_size<1]; eff_size2 <- setdiff(eff_size, eff_size1)
  #print(eff_size[1:50]); print(eff_size1[1:30]); print(eff_size2[1:30])
  eff_size_use <- {if(SIM=="SIM3") eff_size1 else eff_size}
  #print(length(eff_size_use))
  cutoff <- length(eff_size_use)
  #print(cutoff)
  alpha1 <- rnorm(n1+n2, 0, eta); alpha2 <- rnorm(n1+n2, 0, eta)
  epsilon1 <- matrix(rnorm(m*(n1+n2), 0, sigma), nrow = m, ncol = n1+n2, byrow = T)
  epsilon2 <- matrix(rnorm(m*(n1+n2), 0, sigma), nrow = m, ncol = n1+n2, byrow = T)
  if (is.null(bl)){
    baseline <- rnorm(m, 7, 2)
  } else{
    baseline <- bl[sample(1:NGENES, m, replace=F)]
  }
  X <- matrix(data = rep(0, m*(n1+n2)), nrow = m, ncol = n1+n2, byrow = T)
  for (i in 1:cutoff){
    X[i,1:n1] <- baseline[i]
    X[i,(n1+1):(n1+n2)] <- log2(2^X[i,1]*eff_size_use[i])
  }
  for (i in (cutoff+1):m){
    X[i,] <- baseline[i]
  }
  ## Xoracle are the oracle expression matrices
  Xoracle1 <- X + epsilon1; Xoracle2 <- X + epsilon2
  ## Y is after adding slide-specific noise (random factor)
  Y1 <- Y2 <- matrix(data = rep(0, m*(n1+n2)), nrow = m, ncol = n1+n2, byrow = T)
  for (i in 1:(n1+n2)){
    Y1[,i] <- Xoracle1[,i] + alpha1[i]
    Y2[,i] <- Xoracle2[,i] + alpha2[i]
  }
  rownames(Y1) <- rownames(Y2) <- paste("Gene", as.character(1:m), sep = "")
  colnames(Y1) <- colnames(Y2) <- paste("Sample", as.character(1:(n1+n2)), sep = "")
  ## lab are the group labels
  lab <- c(rep(1, n1), rep(-1, n2))
  tstat <- matrix(0, nrow = m, ncol = 2, byrow = T)
  #rawp <- adjp_BH <- adjp_Bonf <- matrix(0, nrow = m, ncol = 8, byrow = T)
  Y1_quant <- normalize.quantile(Y1)
  rownames(Y1_quant) <- paste("Gene", as.character(1:m), sep = "")
  colnames(Y1_quant) <- paste("Sample", as.character(1:(n1+n2)), sep = "")
  #print(sum(colnames(Y1_quant)==colnames(Y1))); print(sum(rownames(Y1_quant)==rownames(Y1)))
  ddd <- t_test_genes(Y1_quant, classlabel = lab, na.rm = TRUE)
  tstat[,1] <- ddd$teststat
  hhh <- superdelta(Y1, classlabel = lab, test = "t", side = "two.sided", methods = "robust", trim = 0.2, baseline = "all")
  tstat[,2] <- hhh$teststat
  #ind1 <- order(abs(ddd$teststat),decreasing = T)[1:TopGenes]
  ind1 <- sample(order(abs(ddd$teststat),decreasing = T)[1:min(1000,round(0.2*m))], min(TopGenes,round(0.2*m)), replace = F)
  #print(ind1)
  #ind2 <- order(abs(hhh$teststat),decreasing = T)[1:round(TopGenes/2)]
  ind2 <- sample(order(abs(hhh$teststat),decreasing = T)[1:min(1000,round(0.2*m))], min(round(TopGenes/2),round(0.2*m)), replace = F)
  #print(ind2)
  ## Next is quantile, easier:
  XY1_train <- as.data.frame(t(rbind(Y1_quant[ind1,], lab)))
  XY1_test <- as.data.frame(t(rbind(Y2[ind1,], lab)))
  colnames(XY1_train)[TopGenes+1] <- colnames(XY1_test)[TopGenes+1] <- "Label"
  #print(dim(XY1_train)); print(dim(XY1_test))
  #View(XY1_train); View(XY1_test)
  mod1.svm <- svm(as.factor(Label) ~ ., data = XY1_train, type = "C-classification")
  pred1.svm <- predict(mod1.svm, newdata = XY1_test, type = "class")
  acry1.svm <- sum(pred1.svm==XY1_test$Label)/(n1+n2)
  ## Next is super-delta, much harder:
  XY2_train <- as.data.frame(t(rbind(Y1[ind2,] - Y1[hhh$`pairing gene index`[ind2],], lab)))
  XY2_test <- as.data.frame(t(rbind(Y2[ind2,] - Y2[hhh$`pairing gene index`[ind2],], lab)))
  for (i in 1:round(TopGenes/2)){
    colnames(XY2_train)[i] <- colnames(XY2_test)[i] <- paste("Difference", as.character(i), sep = "")
  }
  colnames(XY2_train)[round(TopGenes/2)+1] <- colnames(XY2_test)[round(TopGenes/2)+1] <- "Label"
  #View(XY2_train); View(XY2_test)
  mod2.svm <- svm(as.factor(Label) ~ ., data = XY2_train, type = "C-classification")
  pred2.svm <- predict(mod2.svm, newdata = XY2_test, type = "class")
  acry2.svm <- sum(pred2.svm==XY2_test$Label)/(n1+n2)
  return(c(acry1.svm, acry2.svm))
  #return(acry1.svm)
}