#' Two-sample covariance test (Cai et al. 2013)
#' 
#' @description Testing the equality of two high-dimensional covariance matrices
#' based on the \eqn{L_\infinity} norm, proposed in Cai, Liu and Xia (2013)
#' "Two-Sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings".
#' 
#' @param X n1 by p matrix, observation of the first population, columns are features
#' @param Y n2 by p matrix, observation of the second population, columns are features
#' @return A list with the following components:
#'  \item{Mn}{the largest M_ij as defined in Cai (2013) equation (2)}
#'  \item{test.stat}{test statistic (calculated as Mn - 4*log p + log log p)}
#'  \item{pVal}{p-value given by the limiting distribution (Gumbol distribution)}
#'  
#' @seealso \code{Chang.maxBoot.test()}, \code{LC.U.test()}, \code{WL.randProj.test()}, 
#' \code{Schott.Frob.test()}.
Cai.max.test <- function(X, Y){
  p <- ncol(X)
  M <- Cai.Mmatrix(X, Y)
  Mn <- max(as.numeric(M)) 
  test.stat <- Mn - 4*log(p) + log(log(p)) # test statistic
  
  ## p-value given by limiting distribution
  pVal <- 1 - exp(- exp(-test.stat/2) / sqrt(8*pi))
  
  return (list(Mn = Mn, test.stat = test.stat, pVal = pVal))
}

#' The standardized statistics for (cov(Y)-cov(X))^2
#' 
#' @description A helper function to calculate the M-matrix in Cai et al. (2013), equation (2)
#' 
#' @param X n1-by-p matrix, observation of the first population, columns are features
#' @param Y n2-by-p matrix, observation of the second population, columns are features
#' @return A p-by-p numeric matrix, the estimation of (cov(Y)-cov(X))^2 
#' 
#' @seealso \code{Cai.max.test()}.
Cai.Mmatrix <- function(X, Y){
  n1 <- nrow(X)
  n2 <- nrow(Y)
  stats.x <- Cai.theta.sigma(X)
  stats.y <- Cai.theta.sigma(Y)
  M <- (stats.x$sigma - stats.y$sigma)^2 / (stats.x$theta/n1 + stats.y$theta/n2)
  return (M)
}

#' Estimation of cov(X) and its variance
#' 
#' @description A helper function to calculate an estimation of cov(X) (\code{sigma})
#'  and its variance (\code{theta}). See equation (2) in Cai (2013)
#' 
#' @param X n by p observation matrix, columns are features
#' @return A list with the following components:
#'  \item{sigma}{A p-by-p matrix, estimation of cov(X)}
#'  \item{theta}{A p-by-p matrix, estimation of the variance of cov(X)}
#' 
#' @seealso \code{Cai.max.test()}.
Cai.theta.sigma <- function(X){
  X.center <- scale(X, center=TRUE, scale=FALSE)
  n <- nrow(X)
  
  # ## naive implementation
  # p <- ncol(X)
  # sigma = cov(X.center) * (n-1) / n
  # theta = matrix(NA, nrow=p, ncol=p)
  # for (i in 1:p){
  #   for (j in i:p){
  #     Zij = X.center[,i] * X.center[,j]
  #     theta[i, j] = (n-1)/n * var(Zij)
  #   }
  # }
  # theta[lower.tri(theta)] <- t(theta)[lower.tri(t(theta))]
  
  sigma <- t(X.center) %*% X.center / n
  theta <- t(X.center^2) %*% (X.center^2) / n - sigma^2
  
  ## suffices to use only upper triangle elements due to symmetry
  return (list(theta = theta, sigma = sigma))
}
