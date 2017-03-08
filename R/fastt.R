## A convenient function to replace fastT() or rowttest() from
## genefilter package.

fastt <- function(X, classlabel){
    ngenes <- nrow(X); nslides <- ncol(X)
    Tvec <- rep(0, ngenes)

    res<-.C("fastt",as.double(X), as.integer(ngenes),
            as.integer(nslides),
            as.integer(classlabel), Tvec=as.double(Tvec),
            PACKAGE="superdelta")$Tvec
    return(res)
}
