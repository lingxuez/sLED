#' Calculate the sLED test statistic for differential matrices.
#' 
#' @description 
#' A differential matrix is the difference between two symmetric relationship matrices.
#' For any symmetric differential matrix \eqn{D}, sLED test statistic is defined as
#' \deqn{max{ sEig(D), sEig(-D) }}
#' where sEig() is the sparse leading eigenvalue, defined as
#' \deqn{$max_{v} v^T A v$}{max_{v} t(v)*A*v}
#' subject to
#' \deqn{$||v||_2 \leq 1, ||v||_1 \leq s$}{||v||_2 <= 1, ||v||_1 <= s} 
#' 
#' @param Dmat p-by-p numeric matrix, the differential matrix 
#' @param rho a large positive constant such that \eqn{D+diag(rep(rho, p))} and \eqn{-D+diag(rep(rho, p))} 
#'        are positive definite.
#' @param sumabs.seq a numeric vector specifing the sequence of sparsity parameters, each between \eqn{1/sqrt(p)} and 1.
#'        Each sumabs*\eqn{$sqrt(p)$}{sqrt(p)} is the upperbound of the L_1 norm of leading sparse eigenvector \eqn{v}.
#' @param niter the number of iterations to use in the PMD algorithm (see \code{fastPMD()})
#' @param trace whether to trace the progress of PMD algorithm (see \code{fastPMD()})
#' 
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
#' @seealso \code{sLED()}, \code{fastPMD()}.
#' @export
sLEDTestStat <- function(Dmat, rho=1000, sumabs.seq=0.2,
                        niter=20, trace=FALSE) {
  ndim <- 1 ## only consider the first sparse eigenvector
  p <- ncol(Dmat)  
  ntest <- length(sumabs.seq)
  
  results <- list()
  results$sumabs.seq <- sumabs.seq
  results$rho <- rho
  
  results$stats <- rep(NA, ntest)
  results$sign <- rep(NA, ntest)
  results$v <- matrix(NA, nrow=ntest, ncol=p)
  results$leverage <- matrix(NA, nrow=ntest, ncol=p)
  
  ## for each sparsity parameter
  for (i in 1:ntest) {
    sumabs <- sumabs.seq[i]
    pos.out <- fastPMD(Dmat + rho * diag(p), 
                       sumabs=sumabs, trace=trace, niter=niter)    
    neg.out <- fastPMD(- Dmat + rho * diag(p), 
                       sumabs=sumabs, trace=trace, niter=niter)
    
    if (pos.out$d >= neg.out$d) {
      results$sign[i] <- "pos"
      results$stats[i] <- pos.out$d - rho 
      results$v[i, ] <- pos.out$v
      results$leverage[i, ] <- (pos.out$v)^2
    } else {
      results$sign[i] <- "neg"
      results$stats[i] <- neg.out$d - rho 
      results$v[i, ] <- neg.out$v
      results$leverage[i, ] <- (neg.out$v)^2
    }
  }
  
  return(results)
}
