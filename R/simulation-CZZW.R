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
Chang.maxBoot.test <- function(X, Y, nresample=1000) {
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
                                        Tdenominator=Tdenominator, nresample=nresample)
    
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
#' @return A numeric vector with length \code{nresample}, containing the test statistics
#' in wild boostrap repetitions. 
#' 
#' @seealso \code{Chang.maxBoot.test()}.
Chang.wildBootstrap <- function(X, Y, sigma.x, sigma.y, Tdenominator, nresample=1000){
  n1 <- nrow(X)
  n2 <- nrow(Y)
  X.center <- scale(X,center=TRUE,scale=FALSE)
  Y.center <- scale(Y,center=TRUE,scale=FALSE)
  
  test.stat.boot <- rep(NA, nresample)
  for (rep in c(1:nresample)) {
    ## coefficient for wild boostrap
    g.boot.x <- rnorm(n1)
    g.boot.y <- rnorm(n2)
    
    ## bootstrap test statistic
    boot.cov.x <- t(X.center) %*% diag(g.boot.x) %*% X.center - sum(g.boot.x)*sigma.x
    boot.cov.y <- t(Y.center) %*% diag(g.boot.y) %*% Y.center - sum(g.boot.y)*sigma.y
    test.stat.boot[rep] <- max(abs(boot.cov.x / n1 - boot.cov.y / n2) / Tdenominator)
  }

  return(test.stat.boot)
}
