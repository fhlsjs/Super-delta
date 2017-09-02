simu <- function(NGENES = 10000, n1 = 50, n2 = 50, bl = NULL, fc = 2^rnorm(NGENES,0,0.5), SIM = c("SIM1", "SIM2", "SIM3"), eta = 1.5, sigma = 1){
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
  #nonDE_ind <- sample(setdiff(1:NGENES, ind), size = m-cutoff, replace = FALSE)
  alpha <- rnorm(n1+n2, 0, eta)
  epsilon <- matrix(rnorm(m*(n1+n2), 0, sigma), nrow = m, ncol = n1+n2, byrow = T)
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
  ## Xoracle is the oracle expression matrix
  Xoracle <- X + epsilon
  ## Y is after adding slide-specific noise (random factor)
  Y <- matrix(data = rep(0, m*(n1+n2)), nrow = m, ncol = n1+n2, byrow = T)
  for (i in 1:(n1+n2)){
    Y[,i] <- Xoracle[,i] + alpha[i]
  }
  ## lab are the group labels
  lab <- c(rep(1, n1), rep(-1, n2))
  tstat <- rawp <- adjp_BH <- adjp_Bonf <- matrix(0, nrow = m, ncol = 8, byrow = T)
  aaa <- t_test_genes(Xoracle, classlabel = lab, na.rm = TRUE)
  rawp[,1] <- aaa$rawp; tstat[,1] <- aaa$teststat
  bbb <- t_test_genes(normalize.global(Y), classlabel = lab, na.rm = TRUE)
  rawp[,2] <- bbb$rawp; tstat[,2] <- bbb$teststat
  ccc <- t_test_genes(normalize.medIQR(Y), classlabel = lab, na.rm = TRUE)
  rawp[,3] <- ccc$rawp; tstat[,3] <- ccc$teststat
  ptm <- proc.time()
  ddd <- t_test_genes(normalize.quantile(Y), classlabel = lab, na.rm = TRUE)
  rawp[,4] <- ddd$rawp; tstat[,4] <- ddd$teststat
  TM <- proc.time() - ptm  ##Record time for quantile normalization
  eee <- t_test_genes(normalizeCyclicLoess(Y, method = "fast"), classlabel = lab, na.rm = TRUE)
  rawp[,5] <- eee$rawp; tstat[,5] <- eee$teststat
  fff <- superdelta(Y, classlabel = lab, test = "t", side = "two.sided", methods = "mean", baseline = "auto")
  rawp[,6] <- fff$rawp; tstat[,6] <- fff$teststat
  ggg <- superdelta(Y, classlabel = lab, test = "t", side = "two.sided", methods = "median", baseline = "auto")
  rawp[,7] <- ggg$rawp; tstat[,7] <- ggg$teststat
  ptm <- proc.time()
  hhh <- superdelta(Y, classlabel = lab, test="t",side="two.sided",methods="robust",trim=0.2,baseline="auto")
  rawp[,8] <- hhh$rawp; tstat[,8] <- hhh$teststat
  TM2 <- proc.time() - ptm  ##Record time for super-delta + MFTM
  simu_res <- matrix(0, nrow = 6, ncol = 8, byrow = T)
  ## Display results in true/false positive rates!
  rownames(simu_res) <- c("none(power)", "none(type I error)", "BH(power)", "BH(type I error)", "Bonf(power)", "Bonf(type I error)")
  colnames(simu_res) <- colnames(rawp) <- c("oracle", "global", "medIQR", "quantile", "cyc_loess", "mean", "median", "robust")
  for (i in 1:8){
    adjp_BH[,i] <- p.adjust(rawp[,i], method = "BH", n = m)
    adjp_Bonf[,i] <- p.adjust(rawp[,i], method = "bonferroni", n = m)
  }
  for (i in 1:8){
    simu_res[1,i] <- round(100*sum(rawp[1:cutoff, i]<=0.05)/cutoff, 2)
    simu_res[2,i] <- round(100*sum(rawp[(cutoff+1):m, i]<=0.05)/(m-cutoff), 2)
    simu_res[3,i] <- round(100*sum(adjp_BH[1:cutoff, i]<=0.05)/cutoff, 2)
    simu_res[4,i] <- round(100*sum(adjp_BH[(cutoff+1):m, i]<=0.05)/(m-cutoff), 2)
    simu_res[5,i] <- round(100*sum(adjp_Bonf[1:cutoff, i]<=0.05)/cutoff, 2)
    simu_res[6,i] <- round(100*sum(adjp_Bonf[(cutoff+1):m, i]<=0.05)/(m-cutoff), 2)
  }
  return(list(simu_res, tstat, TM, TM2))
}