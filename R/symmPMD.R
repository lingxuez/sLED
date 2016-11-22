#' Symmetric Penalized Matrix Decomposition.
#' 
#' @description 
#' This function solves for the Sparse Principal Component Analysis given 
#' a positive definite matrix A:
#'   \deqn{$max_{v} v^T A v$}{max_{v} t(v)*A*v}
#' subject to
#'   \deqn{$||v||_2 \leq 1, ||v||_1 \leq s$}{||v||_2 <= 1, ||v||_1 <= s} 
#' The solution v is the sparse leading eigenvector, and the corresponding objective 
#' \eqn{$v^T A v$}{t(v)*A*v} is the sparse leading engenvalue.
#' 
#' The algorithm uses an iterative procedure similar to the R Package "PMA", but speeds up the computation
#' using the extra constraint that the decomposition is symmetric. 
#' 
#' @param x p-by-p matrix, symmetric and positive definite
#' @param sumabs sumabs*\eqn{$sqrt(p)$}{sqrt(p)} is the upperbound of the L_1 norm of \eqn{$v$}{v},
#'       controling the sparsity of solution. Must be between \eqn{$1/sqrt(p)$}{1/sqrt(p)} and 1.
#' @param v the starting value of the algorithm, either a pre-calculated first singular vector of x, or NULL.
#' @param niter number of iterations to perform the iterative optimizations
#' @param trace whether to print tracing info during optimization
#' 
#' @return A list containing the following components:
#'  \item{v}{the sparse leading eigenvector v}
#'  \item{d}{the sparse leading eigenvalue \eqn{$d=v^T A v$}{d=t(v)*A*v }}
#'  \item{sumabs}{sumabs*\eqn{$sqrt(p)$}{sqrt(p)} is the upperbound of the L_1 norm of \eqn{$v$}{v}}
#' 
#' @references Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Covariance Matrices, 
#' with Application to Detecting Schizophrenia Risk Genes", arXiv:1606.00252.
#' @references Witten, Tibshirani and Hastie (2009), "A penalized matrix decomposition, 
#' with applications to sparse principal components and canonical correlation analysis", Biostatistics 10(3):515-534.
#' 
#' @export
symmPMD <- function(x, sumabs=0.3, niter=50, v=NULL, trace=TRUE) {
  if (!isSymmetric(x)) {
    stop("x must be a symmetric matrix.")
  }
  
  x[is.na(x)] <- mean(x[!is.na(x)])
  p <- nrow(x)
  sumabsv <- sqrt(p) * sumabs
  K <- 1 
  
  # initial value: first singular vector
  if (is.null(v)) {
    v <- rARPACK::eigs_sym(x, K, "LA")$vectors
  }
  
  # main algorithm (only works for K=1!)
  out <- solvePMD(x, sumabsv=sumabsv, v=v, niter=niter, trace=trace)
  
  return(list(v=out$v, d=out$d, v.init=out$v.init, sumabs=sumabs))
}


#' Solving symmetric Penalized Matrix Decomposition
#'
#' @description 
#' An iterative algorithm that solves the Sparse Principal Component Analysis problem: given 
#' a positive definite matrix A:
#'   \deqn{$max_{v} v^T A v$}{max_{v} t(v)*A*v }
#' subject to
#'   \deqn{$||v||_2 \leq 1, ||v||_1 \leq s$}{||v||_2 <= 1, ||v||_1 <= s} 
#' The solution v is the sparse leading eigenvector, and the corresponding objective 
#' \eqn{$v^T A v$}{t(v)*A*v} is the sparse leading engenvalue.
#' 
#' @param x p-by-p matrix, symmetric and positive definite
#' @param sumabsv the upperbound of the L_1 norm of \eqn{$v$}{v}, controling the sparsity of solution. 
#'       Must be between 1 and \eqn{$sqrt(p)$}{sqrt(p)}.
#' @param v the starting value of the algorithm.
#' @param niter number of iterations to perform the iterative optimizations
#' @param trace whether to print tracing info during optimization
#' 
#' @return A list containing the following components:
#'  \item{v}{the sparse leading eigenvector v}
#'  \item{d}{the sparse leading eigenvalue \eqn{$d=v^T A v$}{d=t(v)*A*v }}
#'  \item{v.init}{the initial value of v}
#'
#' @seealso \code{symmPMD()}.
solvePMD <- function(x, sumabsv, v, niter=50, trace=TRUE) {
  if (!isSymmetric(x)) {
    stop("x must be a symmetric matrix.")
  }
  
  ## initialize
  v.init <- v
  oldv <- rnorm(ncol(x))
  
  ## iterative updates
  for (iter in 1:niter) {   
    if (sum(abs(oldv - v)) <= 1e-07) {
      break
    }
    if (trace) {
      cat(iter, fill=FALSE)
    }
    
    oldv <- v
    
    ## v <- S(X*oldv, lamv) / ||S(X*oldv, lamv)||_2
    ## where S() is soft threshold, lamv >= 0 is such that ||v||_1=sumabsv
    argv <- x %*% v
    lamv <- BinarySearch(argv, sumabsv)
    sv <- soft(argv, lamv)
    v <- matrix(sv / l2n(sv), ncol=1)     
  }
  
  ## optimal objective criteria
  d <- as.numeric(t(v) %*% (x %*% v))
  if (trace) {
    cat(fill=TRUE)
  }
  
  return(list(d=d, v=v, v.init=v.init))
}


#' Search soft threshold
#' @description
#' A binary search to find proper soft threshold \code{lamv} such that
#' \deqn{sv = soft(argv, lamv) / ||soft(argv, lamv)||_2, ||sv||_1 = sumabsv}
#' 
#' @param argv the vector to be soft thresholded
#' @param sumabsv upperbound of the L_1 norm of sv
#' @param maxiter max iteration to perform binary search
#' @return the proper threshold level \code{lamv}.
#' 
#' @seealso \code{symmPMD()}.
BinarySearch <- function (argv, sumabsv, maxiter=150) {
  ## no thresholding
  if (l2n(argv) == 0 || sum(abs(argv/l2n(argv))) <= sumabsv) {
    return(0)
  } 
  
  ## binary search
  lam1 <- 0
  lam2 <- max(abs(argv)) - 1e-05
  iter <- 1
  while (iter < maxiter) {
    sv <- soft(argv, (lam1 + lam2)/2)
    if (sum(abs(sv/l2n(sv))) < sumabsv) {
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    if ((lam2 - lam1) < 1e-06) 
      return((lam1 + lam2)/2)
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}

#' L2 norm for vector
#' 
#' @param vec a numeric vector
#' @return the L2 norm of \code{vec}.
l2n <- function (vec) {
  a <- sqrt(sum(vec^2))
  if (a == 0) 
    a <- 0.05
  return(a)
}

#' Soft threshold
#' 
#' @param x a numeric vector
#' @param d the soft threshold level
#' @return the vector after soft thresholding \code{x} at level \code{d}.
#' 
#' @seealso \code{symmPMD()}.
soft <- function (x, d) {
  return(sign(x) * pmax(0, abs(x) - d))
}
