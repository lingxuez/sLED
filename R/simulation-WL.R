#' Two-sample covariance test (Wu and Li 2015)
#' 
#' @description The two-sample covariance test using random projections 
#' proposed in Wu and Li (2015) 
#' "Tests for High-Dimensional Covariance Matrices Using Random Matrix Projection".
#' 
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' @param nproj number of random projections to use
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
WL.randProj.test <- function(X, Y, nproj=100) {
  p <- ncol(X)
  n1 <- nrow(X)
  n2 <- nrow(Y)
  
  ## random projections
  Fstat <- rep(NA, nproj)
  for (i in 1:nproj) {
    Rvec <- rnorm(p, mean=0, sd=1)
    Sx <- t(Rvec) %*% var(X) %*% Rvec * ((n1-1)/n1)
    Sy <- t(Rvec) %*% var(Y) %*% Rvec * ((n2-1)/n2)
    
    ## variance stabilizing transformation of F-stat=Sx/Sy
    Fstat[i] <- log(Sy / Sx) / sqrt(2/n1 + 2/n2)
  }
  
  ## test stat
  test.stat <- max(Fstat)
  pVal <- pnorm(test.stat, mean=0, sd=1, lower.tail=FALSE)
  pVal <- 1 - (pnorm(test.stat, mean=0, sd=1, lower.tail=TRUE))^nproj
  
  return(list(test.stat=test.stat, pVal=pVal))
}