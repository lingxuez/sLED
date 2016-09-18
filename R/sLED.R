#' A Function that performs the sparse leading eigenvalue driven (sLED) Test for 
#' tow-sample high-dimensional covariance and relationship matrices.
#' 
#' @description 
#' Suppose X, Y are p-dimensional random vectors independently coming from two populations.
#' Let \eqn{D} be the differential matrix given by
#' \deqn{D = A(Y) - A(X)} 
#' sLED tests the following hypothesis:
#' \deqn{H_0: D=0 versus H_1: D != 0}
#' where A() represents some p-by-p relationship matrix among features, e.g., covariance matrices, 
#' correlation matrices, as well as the weighted adjacency matrices defined as
#' \deqn{A_{ij} = |corr(i, j)|^b}
#' for some constant b > 0, 1 <= i, j <= p.
#' 
#' @param X n1-by-p matrix for samples from the first population. 
#'        Rows are samples/observations, while columns are the features.
#' @param Y n2-by-p matrix for samples from the second population. 
#'        Rows are samples/observations, while columns are the features.
#' @param adj.beta a positive number representing the power to transform correlation matrices 
#'        to weighted adjacency matrices by \eqn{A_{ij} = |r_ij|^adj.beta}, 
#'        where \eqn{r_ij} represents the Pearson correlation.
#'        When adj.beta=0, the regular correlation marix \eqn{(r_ij)} is used.
#' @param rho a large positive constant such that \eqn{A(X)-A(Y)+diag(rep(rho, p))} is positive definite.
#' @param sumabs.seq a numeric vector specifing the sequence of sparsity parameters to use, 
#'        each between \eqn{1/sqrt(p)} and 1.
#' @param npermute number of permutations to use, default is 100
#' @param seeds a numeric vector with the length equals to npermute, specifying the  
#'        seed for randomly permutating samples in each permutation. 
#'        If set tu NULL, then the default seeds are used.
#' @param verbose whether to print the progress during permutation tests
#' @param niter the number of iterations to use in the PMD algorithm (see \code{fastPMD()})
#' @param trace whether to trace the progress of PMD algorithm (see \code{fastPMD()})
#'
#' @return Tn the test statistic
#' @return Tn.perm the test statistic for permuted samples
#' @return Tn.perm.sign the sign for permuted samples: "pos" if the permuted test statistic is given by sEig(D), 
#'          and "neg" if is given by sEig(-D)
#' @return pVal the p-value of sLED test  
#' @return sumabs.seq  the sequence of sparsity parameters
#' @return rho the positive constant used
#' @return stats the sequence of test statistics
#' @return sign signs when using different sparsity parameters,
#'          "pos" if the test statistic is given by sEig(D), and "neg" if is given by sEig(-D).
#' @return v the sequence of sparse leading eigenvectors
#' @return leverage the leverage score for genes (defined as \eqn{v^2} element-wise) using different sparsity parameters
#' 
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Differential Matrices, 
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.
#' 
#' @seealso \code{getTestStat()}, \code{fastPMD()}.
#' @export
sLED <- function(X, Y, adj.beta=0, rho=1000, 
                 sumabs.seq=0.2, npermute=100, seeds=NULL, verbose=TRUE, niter=20, trace=FALSE,
                 useMC=TRUE, ncore=2) {
  
  n1 <- nrow(X)
  n2 <- nrow(Y)
  # center observations
  X <- scale(X, center=TRUE, scale=FALSE)
  Y <- scale(Y, center=TRUE, scale=FALSE)
  Z <- rbind(X, Y)
  
  # test statistic
  D.hat <- getDiffMatrix(X, Y, adj.beta)  
  pma.results <- sLEDTestStat(Dmat=D.hat, rho=rho,
                            sumabs.seq=sumabs.seq, niter=niter, trace=trace)  
  Tn <- pma.results$stats
 
  ## permutation test
  hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
  if (useMC) {
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.")
    }
    ## use multi-core parallelization
    permute.results <- sLEDpermuteMC(Z=Z, n1=n1, n2=n2, adj.beta=adj.beta, rho=rho, ncore=ncore,
                                   sumabs.seq=sumabs.seq, npermute=npermute, seeds=seeds, 
                                   verbose=verbose, niter=niter, trace=trace) 
  } else {
    permute.results <- sLEDpermute(Z=Z, n1=n1, n2=n2, adj.beta=adj.beta, rho=rho,
                                  sumabs.seq=sumabs.seq, npermute=npermute, seeds=seeds, 
                                  verbose=verbose, niter=niter, trace=trace)
  }

  ## p-value
  pVal = rowSums(permute.results$Tn.permute > Tn) / npermute
  
  return(c(pma.results, permute.results, list(Tn = Tn, pVal = pVal)))
}
