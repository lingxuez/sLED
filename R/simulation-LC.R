#' Two-sample covariance test (Li and Chen 2012)
#' 
#' @description Testing the equality of two high-dimensional covariance matrices
#' based on the \eqn{L_2} norm, proposed in Li and Chen (2012)
#' "Two Sample Tests for High-dimensional Covariance Matrices"
#'
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' @return A list containing the following components:
#'  \item{Tn}{the U statistic for ||cov(X)-cov(Y)||_F^2}
#'  \item{Tn.sd}{the estimated standard deviation of \code{Tn}}
#'  \item{test.stat}{test statistic, \code{Tn}/\code{Tn.sd}}
#'  \item{pVal}{the p-value calculated using the limiting distribution (standard normal)}
#' 
#' @seealso \code{Cai.max.test()}, \code{Chang.maxBoot.test()}, 
#' \code{WL.randProj.test()}, \code{Schott.Frob.test()}.
LC.U.test <- function(X, Y){
  n1 = nrow(X)
  n2 = nrow(Y)
  X.center = scale(X, center=TRUE, scale=FALSE)
  Y.center = scale(Y, center=TRUE, scale=FALSE)
  
  ## the A matrix in equation (2.1) 
  A.x = LC.Avalue(X.center)
  A.y = LC.Avalue(Y.center)
  ## the C matrix in equation (2.2)
  C.xy = LC.Cvalue(X.center, Y.center)
  
  ## test statistic and estimated sd
  Tn = A.x + A.y - 2*C.xy 
  Tn.sd = A.x*2/n1 + A.y*2/n2
  
  ## limiting distribution is standard normal
  pVal = pnorm(Tn/Tn.sd, lower.tail=FALSE)
  
  return(list(Tn = Tn, Tn.sd = Tn.sd, test.stat=Tn/Tn.sd, pVal = pVal))
}

#' U-statistic of \eqn{tr(cov(X)^2)}
#' 
#' @description Calculate the A value as in Li and Chen (2012) equation (2.1), 
#' for estimating \eqn{tr(cov(X)^2)}.
#' As suggested in the paper, only the first term is used after centering data.
#' 
#' @param X n by p matrix, columns are features, centered such that colMean=0
#' @return A number, the U statistic for \eqn{tr(cov(X)^2)}.
LC.Avalue <- function(X){
  n = nrow(X)
  
  ## for each pair of p-dimensional X1, X2 i.i.d., 
  ## t(X1) %*% X2 (scalar) is an unbiased estimation of tr(cov(X)^2)
  sample.est <- (X %*% t(X))^2
  diag(sample.est) <- 0 ## exclude self-self products
  Avalue <- sum(sample.est) / (n * (n-1))
  
  return (Avalue)
}

#' U-statistic of \eqn{tr(cov(X)cov(Y))}.
#' @description Calculate the C value as in Li and Chen (2012) equation (2.2), 
#' for estimating \eqn{tr(cov(X)cov(Y))}.
#' As suggested in the paper, only the first term is used after centering data.
#' 
#' @param X n1 by p matrix of the first population, centered with colMean=0
#' @param Y n2 by p matrix of the second population, centered with colMean=0
#' @return A number, the U statistic for \eqn{tr(cov(X)cov(Y))}.
LC.Cvalue <- function(X, Y){
  n1 = nrow(X)
  n2 = nrow(Y)
  
  ## for each p-dimensional X1, Y1 independent,
  ## t(X1) %*% Y1 (scalar) is an unbiased estimation of tr(cov(X)cov(Y))
  sample.est <- (X %*% t(Y))^2 
  Cvalue <- sum(sample.est) / (n1 * n2)
  
  return (Cvalue)
}
