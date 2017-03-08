\name{fastt}
\alias{fastt}
\title{A C-implementation of two-sample t-test for matrix data.}

\description{ A convenient function to replace fastT() or rowttest() from
 genefilter package.  Only Student's two sample t-test is
  implemented in the current version.  }

\usage{ fastt(X, classlabel) }

\arguments{
  \item{X}{
    A data frame or matrix, with \eqn{p} rows corresponding to
    variables (hypotheses) and \eqn{n} columns to observations. In the
    case of gene expression data, rows correspond to genes and columns to
    mRNA samples. The data can be read using \code{\link{read.table}}.
  }

  \item{classlabel}{ A vector of integers corresponding to observation
    (column) class labels. For \eqn{k} classes, the labels must be
    integers between 0 and \eqn{k-1}.  We only implemented two-sample
    test as of 08/01/2015.
  }

}

\details{ Some details. }

\value{A vector of two-sample t-statistics
}

\author{Xing Qiu, Yuhang Liu, Jinfeng Zhang}
\seealso{\code{\link{t.test}}}
\examples{
## The Golub data.  On my computer, it takes about 3 seconds to run, YMMV.
## The result is a data.frame with two columns: teststat and rawp.

data(golub)
classlabel<-golub.cl
system.time(tvec1 <- fastt(golub, classlabel))

## You can try the standard t.test(). It is much slower.
\dontrun{system.time(tvec2 <- apply(golub, 1, function(x) t.test(x~classlabel, var.equal=TRUE)$statistics))}

}

\keyword{htest}% __ONLY ONE__ keyword per line
