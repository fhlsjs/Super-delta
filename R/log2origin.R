log2origin<-function(X){
  if (is.numeric(X)==0)
    stop("Error: Input data must be a numeric matrix or vector!")
  Y<-2^X
  return(Y)
}

#View(TA.probe)
#log2origin(unlist(TA.probe[1:3,1:3]))
#View(unlist(TA.probe[1:3,1:3]))
#View(as.matrix(TA.probe)[1:3,1:3])
#log2origin(as.matrix(TA.probe)[1:3,1:3])
