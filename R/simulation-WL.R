#' Two-sample covariance test (Wu and Li 2015)
#' 
#' @description The two-sample covariance test using random projections 
#' proposed in Wu and Li (2015) 
#' "Tests for High-Dimensional Covariance Matrices Using Random Matrix Projection".
#' 
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' @param nproj number of random projections to use
#' @param useMC logical variable indicating whether to use multicore parallelization.
#'        R packages \code{parallel} and \code{doParallel} are required if set to \code{TRUE}.
#' @param mc.cores decide the number of cores to use when \code{useMC} is set to \code{TRUE}.
#' 
#' @return A list containing the following components:
#'  \item{test.stat}{test statistic}
#'  \item{pVal}{the p-value calculated using the limiting distribution 
#'            (max of independent standard normal)}
#' @seealso \code{Cai.max.test()}, \code{Chang.maxBoot.test()}, \code{LC.U.test()},
#'  \code{Schott.Frob.test()}.
#'  
#' @references 
#' Wu and Li (2015) 
#' "Tests for High-Dimensional Covariance Matrices Using Random Matrix Projection",
#' arXiv preprint arXiv:1511.01611.
#'  
#' @export
WL.randProj.test <- function(X, Y, nproj=100, useMC=FALSE, mc.cores=1) {
  p <- ncol(X)
  n1 <- nrow(X)
  n2 <- nrow(Y)
  
  ## random projections
  if (useMC) {
    ## try to use multi-core parallelization
    hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.
           Otherwise set useMC=FALSE for non-parallel computation.")
    }
    Fstats <- parallel::mclapply(1:nproj, oneProj, Xsample=X, Ysample=Y, n1=n1, n2=n2, p=p,
                                mc.cores=mc.cores, mc.preschedule=TRUE)
    Fstats = unlist(Fstats)
    
  } else {
    ## without parallelization
    Fstats <- sapply(1:nproj, oneProj, Xsample=X, Ysample=Y, n1=n1, n2=n2, p=p)
  }

  ## test stat
  test.stat <- max(Fstats)
  pVal <- 1 - (pnorm(test.stat, mean=0, sd=1, lower.tail=TRUE))^nproj
  
  return(list(test.stat=test.stat, pVal=pVal, Fstats=Fstats))
}


#' F-statistic in one random projection
#' 
#' @param i i-th repetition
#' @param Xsample n1 by p matrix, observation of the first population, columns are features
#' @param Ysample n2 by p matrix, observation of the second population, columns are features
#' @param n1 number of rows in Xsample
#' @param n2 number of rows in Ysample
#' @param p number of columns
#' 
#' @seealso \code{WL.randProj.test()}.
oneProj <- function(i, Xsample, Ysample, n1, n2, p) {
  Rvec <- rnorm(p, mean=0, sd=1)
  Sx <- t(Rvec) %*% var(Xsample) %*% Rvec * ((n1-1)/n1)
  Sy <- t(Rvec) %*% var(Ysample) %*% Rvec * ((n2-1)/n2)
  Fstat <- log(Sy / Sx) / sqrt(2/n1 + 2/n2)
  
  return (Fstat[1, 1])
}