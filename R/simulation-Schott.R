#' Two-sample covariance test (Schott 2007)
#' 
#' @description 
#' The two-sample covariance test proposed in Schott (2007) 
#' "A test for the equality of covariance matrices when
#'  the dimension is large relative to the sample sizes".
#'  The original test procedure handles multiple multivariate Normal distributions, 
#'  but for the purpose of our simulation, the implementation is for two populations only. 
#'  
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' 
#' @return A list containing the following components:
#'  \item{test.stat}{test statistic}
#'  \item{pVal}{the p-value calculated using the limiting distribution (standrad normal)}
#'  
#' @seealso \code{Cai.max.test()}, \code{Chang.maxBoot.test()}, 
#' \code{LC.U.test()}, \code{WL.randProj.test()}.
Schott.Frob.test <- function(X, Y) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  n <- n1 + n2
  
  SX <- cov(X) * (n1-1)/n1
  SY <- cov(Y) * (n2-1)/n2
  Spool <- cov(rbind(X, Y)) * (n-1)/n
  
  ## an estimation of tr[(cov(X)-cov(Y))^2] (equation (1))
  frob.est <- (1 - (n1-2) / ((n1-1)*(n1+2))) * sum(SX^2) +
              (1 - (n2-2) / ((n2-1)*(n2+2))) * sum(SY^2) -
    2 * tr(SX %*% SY) - n1 / ((n1-1)*(n1+2)) * (tr(SX))^2 - n2 / ((n2-1)*(n2+2)) * (tr(SY))^2
  
  ## an estimation of sd(frob.est) (equation (4))
  a <- (sum(Spool^2) - (tr(Spool))^2/n) * n^2 / ((n-1)*(n+2))
  frob.sd <- 2 * (n1+n2)/(n1*n2) * a
  
  ## test statistic
  test.stat <- frob.est / frob.sd
  pVal <- pnorm(test.stat, lower.tail=FALSE)
  
  return (list(test.stat=test.stat, pVal=pVal))
}

#' Trace of a matrix
tr <- function(A) {
  return (sum(diag(A)))
}
