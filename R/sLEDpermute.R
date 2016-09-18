#' Get the sLED test statistic on permuted samples.
#' 
#' @param Z (n1+n2)-by-p matrix, containing the samples from two populations with p features
#' @param n1 the first n1 rows in Z represent the first population
#' @param n2 the (n1+1):(n1+n2) rows in Z represent the second population
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
#' @return Tn.perm the test statistic for permuted samples
#' @return Tn.perm.sign the sign for permuted samples: "pos" if the permuted test statistic is given by sEig(D), 
#'          and "neg" if is given by sEig(-D)
#'          
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Differential Matrices, 
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.
#' 
#' @seealso \code{getTestStat()}, \code{sLED()}.
sLEDpermute <- function(Z, n1, n2, adj.beta=0, rho=1000,
                        sumabs.seq=0.2, npermute=100, seeds=NULL, verbose=TRUE, niter=20, trace=FALSE) {
  
  ntest <- length(sumabs.seq)
  Tn.permute <- matrix(NA, nrow=ntest, ncol=npermute)
  Tn.permute.sign <- matrix(NA, nrow=ntest, ncol=npermute)
  
  ## permutation
  if (verbose) {
    cat(npermute, "permutations started:\n")
  }
  
  for (i in 1:npermute) {
    permuteResult <- sLEDOnePermute(i=i, Z=Z, n1=n1, n2=n2, seeds=seeds,  
                                    sumabs.seq=sumabs.seq, 
                                    adj.beta=adj.beta, rho=rho,
                                    verbose=verbose, niter=niter, trace=trace)   
    Tn.permute[, i] <- permuteResult$Tn.permute
    Tn.permute.sign[, i] <- permuteResult$Tn.permute.sign
  }
  
  if (verbose) {
    cat(npermute, "permutations finished.", fill=TRUE)
  }
  
  return(list(Tn.permute = Tn.permute, Tn.permute.sign = Tn.permute.sign))
}


