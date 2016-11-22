#' The sparse leading eigenvalue driven (sLED) test 
#' 
#' @description 
#' The sLED test for two-sample high-dimensional covariance and relationship matrices.
#' Suppose X, Y are p-dimensional random vectors independently coming from two populations.
#' Let \eqn{D} be the differential matrix given by
#' \deqn{D = A(Y) - A(X)} 
#' sLED tests the following hypothesis:
#' \deqn{H_0: D=0 versus H_1: D != 0}
#' where A() represents some p-by-p relationship matrix among features, including covariance matrices, 
#' correlation matrices, or the weighted adjacency matrices defined as
#' \deqn{A_{ij} = |corr(i, j)|^b}
#' for some constant b > 0, 1 <= i, j <= p. 
#' Let A represent the regular correlation matrix when b=0, and covariance matrix when b<0.
#' 
#' @param X n1-by-p matrix for samples from the first population. 
#'        Rows are samples/observations, while columns are the features.
#' @param Y n2-by-p matrix for samples from the second population. 
#'        Rows are samples/observations, while columns are the features.
#' @param adj.beta a positive number representing the power to transform correlation matrices 
#'        to weighted adjacency matrices by \eqn{A_{ij} = |r_ij|^adj.beta}, 
#'        where \eqn{r_ij} represents the Pearson correlation.
#'        When \code{adj.beta=0}, the correlation marix is used. 
#'        When \code{adj.beta<0}, the covariance matrix is used.
#'        The default value is \code{adj.beta=-1}.
#' @param rho a large positive constant such that \eqn{A(X)-A(Y)+diag(rep(rho, p))} is positive definite.
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
#' @param trace logical, whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#'
#' @return A list containing the following components:
#'  \item{Tn}{the test statistic}
#'  \item{Tn.perm}{the test statistic for permuted samples}
#'  \item{Tn.perm.sign}{the sign for permuted samples: 
#'          "pos" if the permuted test statistic is given by sEig(D), 
#'          and "neg" if is given by sEig(-D),
#'          where \code{sEig} denotes the sparse leading eigenvalue.}
#'  \item{pVal}{the p-value of sLED test}  
#'  \item{sumabs.seq}{a numeric vector for a sequence of sparsity parameters. Default is 0.2.
#'              The numbers must be between \eqn{1/sqrt{p}} and 1.}
#'  \item{rho}{a positive constant to augment the diagonal of the differential matrix \eqn{D}
#'            such that \eqn{D + rho*I} becomes positive definite.}
#'  \item{stats}{a numeric vector of test statistics when using different sparsity parameters
#'            (corresponding to \code{sumabs.seq}).}
#'  \item{sign}{a vector of signs when using different sparsity parameters (corresponding to \code{sumabs.seq}).
#'          Sign is "pos" if the test statistic is given by sEig(D), and "neg" if is given by sEig(-D),
#'          where \code{sEig} denotes the sparse leading eigenvalue.}
#'  \item{v}{the sequence of sparse leading eigenvectors, each row corresponds to one sparsity 
#'          parameter given by \code{sumabs.seq}.}
#'  \item{leverage}{the leverage of genes (defined as \eqn{v^2} element-wise) using 
#'          different sparsity parameters. Each row corresponds to one sparsity 
#'          parameter given by \code{sumabs.seq}.}
#' 
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Covariance Matrices, 
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.
#' 
#' @details For large data sets, the multi-core version is recommended:
#' \code{useMC=TRUE} and \code{mc.cores=n}, where \code{n} is the number of cores to use.
#' 
#' @seealso \code{symmPMD()}.
#' 
#' @export
#' 
#' @examples 
#' # Run sLED on a synthetic dataset under the null hypothesis
#' # where cov(X) = cov(Y)
#' n <- 50
#' p <- 100
#' set.seed(99)
#' X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
#' set.seed(42)
#' Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
#' 
#' # run sLED and check the p-value
#' result <- sLED(X=X, Y=Y, npermute=50)
#' result$pVal
#' 
#' 
#' # Run sLED on a synthetic dataset under the alternative hypothesis
#' # where cov(X) != cov(Y), and the difference occur at the first 10 coordinates
#' n <- 50
#' p <- 100
#' set.seed(99)
#' X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
#' s <- 10 ## signals
#' sigma.2 <- diag(p)
#' sigma.2[1:s, 1:s] <- sigma.2[1:s, 1:s] + 0.2
#' set.seed(42)
#' Y2 <- MASS::mvrnorm(n, mu=rep(0, p), Sigma=sigma.2)
#' 
#' # run sLED and check the p-value
#' result <- sLED(X=X, Y=Y2, sumabs.seq=0.25, npermute=100, seeds = c(1:100))
#' result$pVal
#' 
#' # the signalling coordinates detected by sLED
#' which(result$leverage != 0)
#' 
sLED <- function(X, Y, adj.beta=-1, rho=1000, sumabs.seq=0.2, npermute=100, 
                 useMC=FALSE, mc.cores=1, seeds=NULL, verbose=TRUE, niter=20, trace=FALSE) {

  # test statistic
  D.hat <- getDiffMatrix(X, Y, adj.beta)  
  pma.results <- sLEDTestStat(Dmat=D.hat, rho=rho,
                            sumabs.seq=sumabs.seq, niter=niter, trace=trace)  
  Tn <- pma.results$stats
 
  ## permutation test
  n1 <- nrow(X)
  n2 <- nrow(Y)
  Z <- rbind(X, Y)
  permute.results <- sLEDpermute(Z=Z, n1=n1, n2=n2, adj.beta=adj.beta, rho=rho,
                                 sumabs.seq=sumabs.seq, npermute=npermute, 
                                 useMC=useMC, mc.cores=mc.cores, seeds=seeds, 
                                 verbose=verbose, niter=niter, trace=trace)

  ## p-value
  pVal = rowSums(permute.results$Tn.permute > Tn) / npermute
  
  return(c(pma.results, permute.results, list(Tn = Tn, pVal = pVal)))
}


#' Test statistic of sLED
#' 
#' @description 
#' Calculate the sLED test statistic given a differential matrix \eqn{D}.
#' A differential matrix is the difference between two symmetric relationship matrices.
#' For any symmetric differential matrix \eqn{D}, sLED test statistic is defined as
#' \deqn{max{ sEig(D), sEig(-D) }}
#' where \code{sEig()} is the sparse leading eigenvalue, defined as
#' \deqn{$max_{v} v^T A v$}{max_{v} t(v)*A*v}
#' subject to
#' \eqn{$||v||_2 \leq 1, ||v||_1 \leq s$}{||v||_2 <= 1, ||v||_1 <= s}. 
#' 
#' @param Dmat p-by-p numeric matrix, the differential matrix 
#' @param rho a large positive constant such that \eqn{D+diag(rep(rho, p))} and \eqn{-D+diag(rep(rho, p))} 
#'        are positive definite.
#' @param sumabs.seq a numeric vector specifing the sequence of sparsity parameters, each between \eqn{1/sqrt(p)} and 1.
#'        Each sumabs*\eqn{$sqrt(p)$}{sqrt(p)} is the upperbound of the L_1 norm of leading sparse eigenvector \eqn{v}.
#' @param niter the number of iterations to use in the PMD algorithm (see \code{symmPMD()})
#' @param trace whether to trace the progress of PMD algorithm (see \code{symmPMD()})
#' 
#' @return A list containing the following components:
#'  \item{sumabs.seq}{the sequence of sparsity parameters}
#'  \item{rho}{a positive constant to augment the diagonal of the differential matrix 
#'            such that \eqn{D + rho*I} becomes positive definite.}
#'  \item{stats}{a numeric vector of test statistics when using different sparsity parameters
#'            (corresponding to \code{sumabs.seq}).}
#'  \item{sign}{a vector of signs when using different sparsity parameters (corresponding to \code{sumabs.seq}).
#'          Sign is "pos" if the test statistic is given by sEig(D), and "neg" if is given by sEig(-D),
#'          where \code{sEig} denotes the sparse leading eigenvalue.}
#'  \item{v}{the sequence of sparse leading eigenvectors, each row corresponds to one sparsity 
#'          parameter given by \code{sumabs.seq}.}
#'  \item{leverage}{the leverage score for genes (defined as \eqn{v^2} element-wise) using 
#'          different sparsity parameters. Each row corresponds to one sparsity 
#'          parameter given by \code{sumabs.seq}.}
#' 
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Covariance Matrices, 
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.
#' 
#' @seealso \code{sLED()}.
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
    pos.out <- symmPMD(Dmat + rho * diag(p), 
                       sumabs=sumabs, trace=trace, niter=niter)    
    neg.out <- symmPMD(- Dmat + rho * diag(p), 
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


#' The differential matrix
#' 
#' @description Given observations from two populations X and Y, 
#' compute the differential matrix
#' \deqn{D = A(Y) - A(X)}
#' where A() is the covariance matrix, or the weighted adjacency matrices defined as
#' \deqn{A_{ij} = |corr(i, j)|^b}
#' for some constant b > 0, 1 <= i, j <= p.
#' Let A represent the regular correlation matrix when b=0, and covariance matrix when b<0.
#'
#' @param X n1-by-p matrix for samples from the first population. 
#'        Rows are samples/observations, while columns are the features.
#' @param Y n2-by-p matrix for samples from the second population. 
#'        Rows are samples/observations, while columns are the features.
#' @param adj.beta Power to transform correlation matrices to weighted adjacency matrice
#'        by \eqn{A_{ij} = |r_ij|^adj.beta} where \eqn{r_ij} represents the Pearson correlation.
#'        When \code{adj.beta=0}, the correlation marix is used. 
#'        When \code{adj.beta<0}, the covariance matrix is used.
#'        The default value is \code{adj.beta=-1}.
#'      
#' @return The p-by-p differential matrix \eqn{D = A(Y) - A(X)}
getDiffMatrix <- function(X, Y, adj.beta=-1) {
  if (adj.beta < 0) {
    Dmat <- cov(Y) - cov(X)
  } else if (adj.beta == 0) {
    Dmat <- cor(Y) - cor(X)
  } else {
    Dmat <- abs(cor(Y))^adj.beta - abs(cor(X))^adj.beta
  }
  return(Dmat)
}

