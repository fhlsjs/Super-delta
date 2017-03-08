## This is the original delta-sequence function.  For simplicity, I only implemented t.test with two.sided hypothesis
deltaseq_YL <- function(X, classlabel){
  ngenes <- nrow(X)
  ## Let's turn the classlabel into a binary vector
  classlabel <- as.integer(as.factor(classlabel)) - 1
  ## define Xa, Xb
  Xa <- X[, classlabel==0]; Xb <- X[, classlabel==1]
  n1 <- ncol(Xa); n2 <- ncol(Xb)
  ## gene-order is calculated via pooled sample variance/STD
  Va <- apply(Xa, 1, var); Vb <- apply(Xb, 1, var)
  DFnum <- (Va/n1 + Vb/n2)^2
  DFdenom <- Va^2/n1^2/(n1-1) + Vb^2/n2^2/(n2-1)
  DF <- round(DFnum/DFdenom)
  Sp <- (n1-1)*Va + (n2-1)*Vb
  ## order and the reverse order
  o <- order(Sp); ro <- order(o)
  Xord <- X[o,]
  ndeltas <- ngenes %/%2
  if (ngenes %% 2==1){                #odd number of genes
    ## First kind of pairing.
    ## We take the row-differences and append the last one to it
    deltas1 <- rbind(Xord[seq(2, ndeltas*2, 2),] - Xord[seq(1, ndeltas*2-1, 2),], Xord[ngenes, ])
    ## compute the statistics at the deltas level
    delta.stats1 <- fastt(deltas1, classlabel)
    ## now translate them back to genes
    deltas2genes1 <- c(rep(1:ndeltas, each=2), ndeltas+1)
    ## multiply (-1,1, ..., -1, 1, 1) to ensure the sign of the
    ## t-stats are consistent.
    stats1 <- delta.stats1[deltas2genes1] * c(rep(c(-1,1), ndeltas), 1)
    ## the second kind of pairing
    deltas2 <- rbind(Xord[1,], Xord[seq(3, ndeltas*2+1, 2),] - Xord[seq(2, ndeltas*2, 2),])
    delta.stats2 <- fastt(deltas2, classlabel)
    deltas2genes2 <- c(1, rep(2:(ndeltas+1), each=2))
    stats2 <- delta.stats2[deltas2genes2] * c(1, rep(c(-1,1), ndeltas))
  } else {                            #even number of genes
    ## First kind of pairing
    deltas1 <- Xord[seq(2, ndeltas*2, 2),] - Xord[seq(1, ndeltas*2-1, 2),]
    delta.stats1 <- fastt(deltas1, classlabel)
    deltas2genes1 <- rep(1:ndeltas, each=2)
    stats1 <- delta.stats1[deltas2genes1] * rep(c(-1,1), ndeltas)
    ## Second kind of pairing. Leave the top and bottom alone
    deltas2 <- rbind(Xord[1,], Xord[seq(3, ndeltas*2-1, 2),] - Xord[seq(2, (ndeltas-1)*2, 2),], Xord[ngenes,])
    delta.stats2 <- fastt(deltas2, classlabel)
    deltas2genes2 <- c(1, rep(2:ndeltas, each=2), ndeltas+1)
    stats2 <- delta.stats2[deltas2genes2] * c(1, rep(c(-1,1), ndeltas-1), 1)
  }
  ## The less extreme of the two tstats gets reported. Note that I
  ## have to reverse-order these stats so that they are consistent
  ## with the original gene list.
  ## 11/18/2015.  Multiply the test stats by sqrt(2) to boost power.    
  teststat <- sqrt(2)*apply(cbind(stats1, stats2), 1, function(tt) {at <- abs(tt); ifelse(at[1] <= at[2], tt[1], tt[2])})[ro]
  ## The raw pvalue. Only two-sided test is implemented at this moment.
  rawp=.t2p(teststat, df=DF, side="two.sided")
  rr <- data.frame(teststat=teststat, rawp=rawp)
  colnames(rr) <- c("teststat", "rawp")
  return(rr)
}