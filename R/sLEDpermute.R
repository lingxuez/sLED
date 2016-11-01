#' sLED permutation
#' 
#' @description Get the sLED test statistic on permuted samples.
#' 
#' @param Z (n1+n2)-by-p matrix, containing the samples from two populations with p features
#' @param n1 the first n1 rows in Z represent the first population
#' @param n2 the (n1+1):(n1+n2) rows in Z represent the second population
#' @param adj.beta a positive number representing the power to transform correlation matrices 
#'        to weighted adjacency matrices by \eqn{A_{ij} = |r_ij|^adj.beta}, 
#'        where \eqn{r_ij} represents the Pearson correlation.
#'        When adj.beta=0, the regular covariance marix is used.
#' @param rho a large positive constant such that \eqn{A(X)-A(Y)+diag(rep(rho, p))} is positive definite.
#' @param ncore the number of cores to use
#' @param sumabs.seq a numeric vector specifing the sequence of sparsity parameters to use, 
#'        each between \eqn{1/sqrt(p)} and 1.
#' @param npermute number of permutations to use, default is 100
#' @param useMC logical, whether to use multi-core version
#' @param mc.cores a number indicating how many cores to use in parallelization
#' @param seeds a numeric vector with the length equals to \code{npermute}, 
#'        where \code{seeds[i]} specifies the seeding for the i-th permutation. 
#'        Set to \code{NULL} if do not want to specify.
#' @param verbose whether to print the progress during permutation tests
#' @param niter the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param trace whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#' 
#' @return A list containing the following components:
#'  \item{Tn.perm}{a numeric vector with length \code{npermute}, the test statistic in permutations.}
#'  \item{Tn.perm.sign}{a vector of characters with length \code{npermute}, the sign in permutations: 
#'      "pos" if the permuted test statistic is given by sEig(D), and "neg" if is given by sEig(-D),
#'      where \code{sEig} denotes the sparse leading eigenvalue.}
#'          
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Differential Matrices, 
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.
#' 
#' @seealso \code{sLED()}.
sLEDpermute <- function(Z, n1, n2, adj.beta=0, rho=1000,
                        sumabs.seq=0.2, npermute=100, 
                        useMC=FALSE, mc.cores=1, seeds=NULL, verbose=TRUE, niter=20, trace=FALSE) {
  
  ntest <- length(sumabs.seq)
  Tn.permute <- matrix(NA, nrow=ntest, ncol=npermute)
  Tn.permute.sign <- matrix(NA, nrow=ntest, ncol=npermute)
  
  ## permutation
  if (verbose) {
    cat(npermute, "permutation started:\n")
  }
  
  if (useMC) {
    ## try to use multi-core parallelization
    hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.
           Otherwise set useMC=FALSE for non-parallel computation.")
    }
    perm.results <- parallel::mclapply(1:npermute, sLEDOnePermute, 
                                       Z=Z, n1=n1, n2=n2, seeds=seeds,  
                                       sumabs.seq=sumabs.seq, 
                                       adj.beta=adj.beta, rho=rho,
                                       verbose=FALSE, niter=niter, trace=trace,
                                       mc.cores=mc.cores, mc.preschedule=TRUE)
  } else {
    ## without parallelization
    perm.results <- lapply(1:npermute, sLEDOnePermute, 
                                       Z=Z, n1=n1, n2=n2, seeds=seeds,  
                                       sumabs.seq=sumabs.seq, 
                                       adj.beta=adj.beta, rho=rho,
                                       verbose=verbose, niter=niter, trace=trace)
  }
  
  ## extract test statistics and signs
  for (i in 1:npermute) {
    Tn.permute[, i] <- perm.results[[i]]$Tn.permute
    Tn.permute.sign[, i] <- perm.results[[i]]$Tn.permute.sign
  } 
  
  if (verbose) {
    cat("permutations finished.", fill=TRUE)
  } 
  
  return(list(Tn.permute = Tn.permute, Tn.permute.sign = Tn.permute.sign))
}


#' Permute indices
#' 
#' @param n1 number of observations from first population
#' @param n2 number of observations from second population
#' 
#' @return A list with the following components:
#'  \item{i1}{a numeric vector with length \code{n1}, the indices for permuted first population}
#'  \item{i2}{a numeric vector with length \code{n2}, the indices for permuted second population}
#'  
#' @seealso \code{sLEDpermute()}.
permuteIndex <- function(n1, n2){
  i.sample <- sample(1:(n1+n2), replace=FALSE)  
  return(list(i1=i.sample[1:n1], i2=i.sample[(n1+1):(n1+n2)]))  
}


#' One permutation of sLED
#' 
#' @param i the i-th permutation
#' @param Z (n1+n2)-by-p matrix, containing the pooled samples from two populations. Columns are features.
#' @param n1 the first n1 rows in Z represent the first population
#' @param n2 the (n1+1):(n1+n2) rows in Z represent the second population
#' @param adj.beta a positive number representing the power to transform correlation matrices 
#'        to weighted adjacency matrices by \eqn{A_{ij} = |r_ij|^adj.beta}, 
#'        where \eqn{r_ij} represents the Pearson correlation.
#'        When adj.beta=0, the regular correlation marix \eqn{(r_ij)} is used.
#' @param rho a large positive constant such that \eqn{A(X)-A(Y)+diag(rep(rho, p))} is positive definite.
#' @param sumabs.seq a numeric vector specifing the sequence of sparsity parameters to use, 
#'        each between \eqn{1/sqrt(p)} and 1.
#' @param seeds a numeric vector, where \code{seeds[i]} specifies the seeding for 
#'        the i-th permutation. Set to \code{NULL} if do not want to specify.
#'
#' @return A list containing the following components:
#'  \item{Tn.permute}{A number represents the test statistic in this permutation.}
#'  \item{Tn.permute.sign}{A string, "pos" if the test statistic is given by sEig(D), and "neg" if is given by sEig(-D),
#'        where \code{sEig} denotes the sparse leading eigenvalue.}
#'        
#' @seealso \code{sLEDpermute()}.
sLEDOnePermute <- function(i, Z, n1, n2, 
                           seeds=NULL, sumabs.seq=0.2, adj.beta=0, rho=1000,
                           verbose=TRUE, niter=20, trace=FALSE) {
  if (!is.null(seeds)){
    set.seed(seeds[i])
  }
  
  i.permute <- permuteIndex(n1, n2)
  D.permute <- getDiffMatrix(X=Z[i.permute$i1, ], Y=Z[i.permute$i2, ], adj.beta)
  test.permute <- sLEDTestStat(Dmat=D.permute, rho=rho, 
                               sumabs.seq=sumabs.seq, niter=niter, trace=trace)  
  
  if (verbose  && (i %% 10)==0){
    cat(i, ",")
  }
  return(list(Tn.permute=test.permute$stats, 
              Tn.permute.sign=test.permute$sign))
}






