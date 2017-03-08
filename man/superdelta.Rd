\name{superdelta}
\alias{superdelta}
\title{Super-delta sequence method for testing gene significance.}

\description{ The Super-delta sequence method is a combination of a)
  gene normalization procedure based on pairwise normalization; b) and a
  gene selection procedure based on a suitable test statistic.  It must
  be applied to log-transformed, but \emph{not normalized} gene
  expression level data.  Only Student's two sample t-test is
  implemented in the current version.  }

\usage{ superdelta(X, classlabel, test="t", side=c("two.sided", "less", "greater"), methods="robust", trim=0.2, baseline="auto", ...) }

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

  \item{test}{A character string specifying the statistic to be used to
    test the null hypothesis of no association between the variables and
    the class labels.\cr If \code{test="t"}, the tests are based on
    Student's two-sample t-statistics (equal variances). \cr If
    \code{test="f"}, the tests are based on F-statistics.\cr If
    \code{test="pairt"}, the tests are based on paired t-statistics. We
    only implemented two-sample test as of 08/01/2015.
  }

  \item{side}{A character string specifying the type of rejection region.\cr
    If \code{side="two.sided"}, two-tailed tests, the null hypothesis is rejected for large absolute values of the test statistic.\cr
    If \code{side="more"}, one-tailed tests, the null hypothesis is
    rejected for large values of the test statistic.\cr
    If \code{side="less"}, one-tailed tests,  the null hypothesis is rejected for small values of the test statistic.
  }

  \item{methods}{A list of estimating methods to be used with the
    super-delta sequence method.  The default value "robust" should be
    the best method for most applications.  Other options include
    "mean", "median", and "mclust" (takes about 100 times more computing
    time). These advanced options are implemented for specialized
    situations only.
  }

  \item{trim}{ A real number between 0 and 1.  It controls the
    proportion of "extreme" test statistics to be exluded from the test
    statistic estimation step.  This parameter is similiar in idea, but
    not identical to, the same parameter in function \code{mean()} for
    computing the robust sample mean.  
  }

  \item{baseline}{ See Section details. }

  \item{...}{Other advanced options. Currently only used for mclust method. }
}

\details{ The super-delta sequence method is a combination of a) gene
  normalization procedure based on pairwise normalization; b) and a gene
  selection procedure based on a suitable test statistic.  It must be
  applied to background corrected, log-transformed, but \emph{not
  normalized} gene expression level data.  Only Student's two sample
  t-test is implemented in the current version.  Future plan includes
  paired t-test and ANOVA F-test.  For more technical details of the
  super-delta sequence method please see the reference.

  The \code{trim} parameter which controls how robust the super-delta
  sequence method is.  In theory, a safe \code{trim} value should be no
  less than the proportion of true differentially expressed genes.  The
  default value 0.2 should work for most practical applications for two
  reasons: a) it is not that common to have more than 20\% genes
  selected as DEG; b) even if more than 20\% of genes are DEG, chances
  are some of the "outliers" caused by pairing with the up-regulated
  genes get cancelled with those caused by the down-regulated genes; c)
  trimming is by nature a conservative operation, so some residual bias
  will not cause noticeable inflation of type I error.  Please see the
  reference for more details.

  The baseline paramter defines the set of baseline genes.  If set to
  "all" or 0, every gene is normalized to every other gene.  If set to a
  vector of integers, genes with these indices will be used as baseline
  genes, so every gene is normalized by these baseline genes only.  This
  can save computation time and memory drastically for large data.  If
  it is a positive integer, a random subset of genes of this size is
  created by sample(ngenes, baseline) and used as baseline genes.  The
  default ("auto") is to use 1000 genes when the number of genes is larger
  than 1000 and all genes for smaller data. }

\value{If one method is specified (the default behavior), it produce a
  data frame with two components: 
  \item{teststat}{Vector (or matrix, if more than one method is
    specified) of test statistics.}
  \item{rawp}{Vector (or matrix) of raw (unadjusted) \eqn{p}-values.}
  If more than one method is specified, \code{teststat} and \code{rawp}
  are organized as data frames and put together as a list.
}

\author{Xing Qiu, Yuhang Liu, Jinfeng Zhang}
\seealso{\code{\link{mean}}, \code{\link{t.test}}, \code{\link{p.adjust}}}
\examples{
## The Golub data.  On my computer, it takes about 3 seconds to run, YMMV.
## The result is a data.frame with two columns: teststat and rawp.

data(golub)
classlabel<-golub.cl
system.time(res1 <- superdelta(golub, classlabel))

## Golub data has 3,051 genes. You can use all genes as baseline genes,
## in which case it may take about 10 seconds to run.
\dontrun{system.time(res2 <- superdelta(golub, classlabel, baseline="all"))}

## use both robust and the mclust methods together with a lower-side
## test to find up-regulated genes.  A random set of 500 genes are used.
## The results are organized in a list of two data.frames instead. This
## command takes several minutes to run because Mclust is
## computationally very expensive.

\dontrun{library(mclust); res3 <- superdelta(golub, classlabel,
methods=c("robust", "mclust"), side="lower", baseline=500)}
}

\keyword{htest}% __ONLY ONE__ keyword per line
