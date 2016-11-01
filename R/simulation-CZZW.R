#' Two-sample covariance test (Chang et al. 2016)
#' 
#' @description The two-sample test for high-dimensional covariance matrices in
#'  Chang, Zhou, Zhou, and Wang (2016)
#' "Comparing Large Covariance Matrices under Weak Conditions on the Dependence Structure".
#' Using a wild boostrap procedure, the test statistic is essentially the square root of
#'  the statistic proposed in Cai, Liu and Xia (2013).
#' 
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' @param nresample the number of bootstraps to perform
#' @param useMC logical variable indicating whether to use multicore parallelization.
#'        R packages \code{parallel} and \code{doParallel} are required if set to \code{TRUE}.
#' @param mc.cores decide the number of cores to use when \code{useMC} is set to \code{TRUE}.
#'        
#' @return A list with the following components:
#'  \item{test.stat}{test statistic}
#'  \item{test.stat.boot}{bootstrap test statistics, a numeric vector with length nresample}
#'  \item{pVal}{bootstrap p-value}
#' 
#' @seealso \code{Cai.max.test()}, \code{LC.U.test()}, \code{WL.randProj.test()}, \code{Schott.Frob.test()}.
#' 
#' @references 
#' Chang, Zhou, Zhou, and Wang (2016)
#' "Comparing Large Covariance Matrices under Weak Conditions on the Dependence Structure",
#' arXiv preprint arXiv:1505.04493.
#' 
#' @export
Chang.maxBoot.test <- function(X, Y, nresample=1000, useMC=TRUE, mc.cores=1) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  
  # statistics for estimating cov(X), cov(Y) and their variance (equation (2.2)).
  # note that these are the same as in Cai, Liu and Xia (2013).
  stats.x <- Cai.theta.sigma(X)
  stats.y <- Cai.theta.sigma(Y)

  # the test statistic; the denominator keeps unchanged across permutations
  Tdenominator <- sqrt(stats.x$theta/n1 + stats.y$theta/n2)
  test.stat <- max(abs(stats.x$sigma - stats.y$sigma) / Tdenominator)

  # wild bootstrap
  test.stat.boot <- Chang.wildBootstrap(X=X, Y=Y, sigma.x=stats.x$sigma, sigma.y=stats.y$sigma, 
                                        Tdenominator=Tdenominator, nresample=nresample,
                                        useMC=useMC, mc.cores=mc.cores)
    
  pVal <- sum(test.stat.boot > test.stat) / nresample
  return(list(test.stat=test.stat, test.stat.boot=test.stat.boot, pVal=pVal))
}

#' Wild bootstrap
#'
#' @description Perform the wild boostrap procedure for testing two-sample covariance
#' matrices in Chang et al. (2016).
#' 
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' @param sigma.x p-by-p matrix, sample covariance of X, pre-calculated
#' @param sigma.y p-by-p matrix, sample covariance of X, pre-calculated
#' @param Tdenominator the denominator of the test statistics, which is pre-calculated 
#'        and remains unchanged in bootstrap.
#' @param nresample the number of bootstraps to perform
#' @param useMC logical variable indicating whether to use multicore parallelization.
#'        R packages \code{parallel} and \code{doParallel} are required if set to \code{TRUE}.
#' @param mc.cores decide the number of cores to use when \code{useMC} is set to \code{TRUE}.
#'        
#' @return A numeric vector with length \code{nresample}, containing the test statistics
#' in wild boostrap repetitions. 
#' 
#' @seealso \code{Chang.maxBoot.test()}.
Chang.wildBootstrap <- function(X, Y, sigma.x, sigma.y, Tdenominator, nresample=1000,
                                useMC=TRUE, mc.cores=1){
  n1 <- nrow(X)
  n2 <- nrow(Y)
  X.center <- scale(X,center=TRUE,scale=FALSE)
  Y.center <- scale(Y,center=TRUE,scale=FALSE)
  
  ## bootstrap
  if (useMC) {
    ## try to use multi-core parallelization
    hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.
           Otherwise set useMC=FALSE for non-parallel computation.")
    }
    test.stat.boot <- parallel::mclapply(1:nresample, oneBootstrapStat, 
                                         X.center=X.center, Y.center=Y.center, n1=n1, n2=n2, 
                                         sigma.x=sigma.x, sigma.y=sigma.y, Tdenominator=Tdenominator,
                                         mc.cores=mc.cores, mc.preschedule=TRUE)
    test.stat.boot = unlist(test.stat.boot)
    
  } else {
      ## without parallelization
      test.stat.boot <- sapply(1:nresample, oneBootstrapStat, 
                               X.center=X.center, Y.center=Y.center, n1=n1, n2=n2, 
                               sigma.x=sigma.x, sigma.y=sigma.y, Tdenominator=Tdenominator)
  }

  return(test.stat.boot)
}

#' One repetition in wild bootstrap
#'
#' @description Helper function for \code{Chang.wildBootstrap()}.
#' 
#' @param X.center n1 by p matrix, observation of the first population, 
#'        columns are centered features
#' @param Y.center n2 by p matrix, observation of the second population, 
#'        columns are centered features
#' @param n1 number of rows of \code{X.center}
#' @param n2 number of rows of \code{Y.center}
#' @param sigma.x p-by-p matrix, sample covariance of X, pre-calculated
#' @param sigma.y p-by-p matrix, sample covariance of X, pre-calculated
#' @param Tdenominator the denominator of the test statistics, pre-calculated
#'        
#' @return The test statistic in this boostrap repetition. 
#' 
#' @seealso \code{Chang.maxBoot.test()}.
oneBootstrapStat <- function(i, X.center, Y.center, n1, n2, sigma.x, sigma.y, Tdenominator) {

  ## coefficient for wild boostrap
  g.boot.x <- rnorm(n1)
  g.boot.y <- rnorm(n2)
  
  ## bootstrap test statistic
  boot.cov.x <- t(X.center) %*% diag(g.boot.x) %*% X.center - sum(g.boot.x)*sigma.x
  boot.cov.y <- t(Y.center) %*% diag(g.boot.y) %*% Y.center - sum(g.boot.y)*sigma.y
  
  test.stat <- max(abs(boot.cov.x / n1 - boot.cov.y / n2) / Tdenominator)
  
  return(test.stat)
}
