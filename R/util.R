#' A helper function to permute sample indices
permuteIndex <- function(n1, n2){
  i.sample <- sample(1:(n1+n2), replace=FALSE)  
  
  return(list(i1=i.sample[1:n1], i2=i.sample[(n1+1):(n1+n2)]))  
}


#' A helper function to get the differential matrix
#' \deqn{D = A(Y) - A(X)}
#' where A() represents the weighted adjacency matrices defined as
#' \deqn{A_{ij} = |corr(i, j)|^b}
#' for some constant b > 0, 1 <= i, j <= p.
#' We let A represent the regular correlation matrix when b=0.
#'  
#' @param X n1-by-p matrix for samples from the first population. 
#'        Rows are samples/observations, while columns are the features.
#' @param Y n2-by-p matrix for samples from the second population. 
#'        Rows are samples/observations, while columns are the features.
#' @param adj.beta Power to transform correlation matrices to weighted adjacency matrice
#'        by \eqn{A_{ij} = |r_ij|^adj.beta} where \eqn{r_ij} represents the Pearson correlation.
#'        adj.beta=0 indicates using the regular correlation marix \eqn{(r_ij)}.
#'      
#' @return Dmat the differential matrix \eqn{D = A(Y) - A(X)}
getDiffMatrix <- function(X, Y, adj.beta) {
  if (adj.beta == 0) {
    Dmat <- cor(Y) - cor(X)
  } else {
    Dmat <- abs(cor(Y))^adj.beta - abs(cor(X))^adj.beta
  }
  return(Dmat)
}

#' A function to perform 1 permutation of sLED
#' 
#' @param i the i-th permutation
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
#' @param seeds a numeric vector with the length equals to npermute, specifying the  
#'        seed for randomly permutating samples in each permutation. 
#'        If set tu NULL, then the default seeds are used.
#' @param verbose whether to print the progress during permutation tests
#' @param niter the number of iterations to use in the PMD algorithm (see \code{fastPMD()})
#' @param trace whether to trace the progress of PMD algorithm (see \code{fastPMD()})
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


