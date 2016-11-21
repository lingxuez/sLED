#' Helper functions to simulate data
library(rARPACK)

#' Simulate samples
#' 
#' @param n1 number of samples in first population
#' @param n2 number of samples in second population
#' @param p feature dimension
#' @param gamma.1 square root of the 1st covariance matrix
#' @param gamma.2 square root of the 1st covariance matrix
#' @param sample.dist base distributions. See \strong{Details} below.
#' 
#' @details The function generates the random variable Z
#' with independent features coming from 
#' the given base distribution \code{sample.dist}. 
#' Then let
#' X = gamma.1 %*% Z1,
#' Y.null = gamma.1 %*% Z2,
#' Y.alt = gamma.2 %*% Z3.
#' we will have cov(X)=cov(Y.null)=sigma.1, and cov(Y.alt)=sigma.2.
#' 
#' The possible values for \code{sample.dist} include
#' \tabular{ll}{
#' "normal" \tab Standard normal distribution. \cr
#' "nb" \tab Negative binomial with mu=2 (mean), size=2 (dispersion),
#'            which implies variance = mu+mu^2/size = 4.
#'            Substract the samples by mean.\cr
#' "t" \tab t-distribution with df=12.\cr
#' "gamma" \tab centralized Gamma(4, 0.5) (i.e., substract the expectation 2)
#' }
#' 
#' @return A list containing
#' \item{X}{Independent samples from 1st population, \eqn{n1} by \eqn{p} matrix.}
#' \item{Y.null}{Independent samples from 1st population, \eqn{n2} by \eqn{p} matrix.}
#' \item{Y.alt}{Independent samples from 2nd population, \eqn{n2} by \eqn{p} matrix.}
simulate.samples <- function(n1, n2, p, gamma.1, gamma.2, sample.dist="normal") {
  if (sample.dist=="normal") {
    Z1 <- matrix(rnorm(n1*p), nrow=n1, ncol=p)
    Z2 <- matrix(rnorm(n2*p), nrow=n2, ncol=p)
    Z3 <- matrix(rnorm(n2*p), nrow=n2, ncol=p)
  } else if (sample.dist == "nb") {
    mu <- 2
    size <- 2
    Z1 <- matrix(rnbinom(n1*p, mu=mu, size=size), nrow=n1, ncol=p) - mu
    Z2 <- matrix(rnbinom(n2*p, mu=mu, size=size), nrow=n2, ncol=p) - mu
    Z3 <- matrix(rnbinom(n2*p, mu=mu, size=size), nrow=n2, ncol=p) - mu
  } else if (sample.dist == "t") {
    df <- 12
    Z1 <- matrix(rt(n1*p, df=df), nrow=n1, ncol=p)
    Z2 <- matrix(rt(n2*p, df=df), nrow=n2, ncol=p)
    Z3 <- matrix(rt(n2*p, df=df), nrow=n2, ncol=p)
  } else if (sample.dist == "gamma") {
    shape <- 4
    scale <- 0.5
    Z1 <- matrix(rgamma(n1*p, shape=shape, scale=scale), nrow=n1, ncol=p) - shape*scale
    Z2 <- matrix(rgamma(n2*p, shape=shape, scale=scale), nrow=n2, ncol=p) - shape*scale
    Z3 <- matrix(rgamma(n2*p, shape=shape, scale=scale), nrow=n2, ncol=p) - shape*scale
  } else {
    stop("Unknown values of sample.dist.")
  }
  
  return(list(X = Z1 %*% gamma.1, Y.null = Z2 %*% gamma.1, Y.alt = Z3 %*% gamma.2))
}


#' Generate covariance matrix
#'
#' Generate two covariance matrices \code{sigma.1} and \code{sigma.2}
#' for null and alternative hypothesis, respectively.
#' For the convenience of later simulations, the square root version
#' \code{gamma.1} and \code{gamm.2} are provided, such that
#' sigma.1 = gamma.1 %*% t(gamma.1),
#' sigma.2 = gamma.2 %*% t(gamma.2).
#' 
#' @param p feature dimension
#' @param cov.method method for generating the covariance matrix under null.
#'            See \strong{Details} below.
#' @param diff.method the method for generating the differential matrix. 
#'            See \strong{Details} below.
#'            
#' @details The \code{cov.method} argument is a character string
#' that specifies the type of covariance matrix to be generated
#' for null hypothesis.
#' Possible values are:
#' 
#' \tabular{ll}{
#' "block.diag" \tab block diagonal with 10 blocks, as in 
#'                  Cai et al. (2013) and Chang et al. (2016).\cr
#' "exp.decay" \tab exponential decay as in 
#'                  Cai et al. (2013) and Chang et al. (2016).\cr
#' "noisy.diag" \tab diagonal matrix with sparse nose,
#'                   as in Cai et al. (2013).                                 
#' }
#' 
#' The \code{diff.method} argument is a character string
#' that specifies the type of differential matrix to be generated,
#' such that \code{sigma.2} = \code{sigma.1} + \code{diff.matrix}
#' for alternative hypothesis.
#' Possible values are:
#' \tabular{ll}{
#' "block" \tab \code{diff.matrix} is non-zero in an s-by-s sub-block,
#'              with entries generated from uniform distribution.\cr
#' }
#'    
#' @return A list of two square root of covariance matrix:
#' \item{gamma.1}{The square root of covariance matrix under null hypothesis.}
#' \item{gamma.2}{The square root of covariance matrix under alternative hypothesis.}
get.cov.matrix <- function(p, cov.method="block.diag", diff.method="block",
                           theta=5, s=NULL,
                           eigengenes=NULL, maxCor=0.9, modProportions=NULL ## for wgcna
                           ) {
  ## base matrix
  sigma.1 <- get.base.matrix(p, cov.method, eigengenes, maxCor, modProportions)
  diff.matrix <- get.diff.matrix(p=p, diff.method=diff.method, theta=theta, 
                                 max.var=max(abs(diag(sigma.1))), s=s)
  sigma.2 <- sigma.1 + diff.matrix 
  if ((!isSymmetric(sigma.1)) || (!isSymmetric(sigma.2))) {
    stop("Something is wrong: sigmas are not symmetric.")
  }
  
  ## make both sigma matrices positive definite
  adj.constant <- adj.base.matrix(sigma.1, sigma.2)
  diag(sigma.1) <- diag(sigma.1)  + adj.constant
  diag(sigma.2) <- diag(sigma.2) + adj.constant
  
  return(list(gamma.1=get.matrix.sqrt(sigma.1),
              gamma.2=get.matrix.sqrt(sigma.2),
              diff.matrix=diff.matrix))
}

#' The square root of a positive definite matrix
#' @param A a symmetric positive definite matrix
#' @return the square root of \eqn{A}.
get.matrix.sqrt <- function(A) {
  A.eig <- eigen(A)
  A.sqrt <- A.eig$vectors %*% diag(sqrt(A.eig$values)) %*% t(A.eig$vectors)
  return (A.sqrt)
}

#' Generate the base matrix 
#' 
#' Note: the base matrix may not be positive definite.
#' @param p feature dimension
#' @param cov.method method for generating the base matrix.
#'            see \strong{Details} below.
#'            
#' @details The \code{cov.method} argument is a character string
#' that specifies the type of base matrix to be generated.
#' Possible values are:
#' 
#' \tabular{ll}{
#' "block.diag" \tab block diagonal with 10 blocks, as in 
#'                  Cai et al. (2013) and Chang et al. (2016).\cr
#' "exp.decay" \tab exponential decay as in 
#'                  Cai et al. (2013) and Chang et al. (2016).\cr
#' "noisy.diag" \tab diagonal matrix with sparse nose,
#'                   as in Cai et al. (2013).                                
#' }
#' @return A \eqn{p} by \eqn{p} base matrix, which will be used
#' to generate the covariance matrix. 
get.base.matrix <- function(p, cov.method="block.diag",
                            eigengenes=NULL, maxCor=0.9, modProportions=NULL ## for wgcna
                            ) {
  
  if (cov.method == "block.diag") {
    ## parameter: block size and entries in block
    block.num <- 10
    off.value <- 0.55
    ## covariance
    A <- matrix(0, nrow=p, ncol=p)
    block.size <- ceiling(p/block.num)
    for (k in c(1:block.num)) {
      i.start <- (k-1)*block.size + 1
      i.end <- min(k*block.size, p)
      A[i.start:i.end, i.start:i.end] <- off.value
    }
    diag(A) <- 1

  } else if (cov.method == "exp.decay") {
    ## parameter: sigma_ij = delta^{|i-j|^rho}
    delta <- 0.5
    rho <- 1
    ## covariance
    A <- matrix(nrow=p, ncol=p)
    A <- delta^(abs(row(A) - col(A))^rho)

  } else if (cov.method == "noisy.diag") {
    A <- diag(p)
    A[upper.tri(A, diag=FALSE)] <- rbinom(n=p*(p-1)/2, size=1, prob=0.05)
    A[lower.tri(A, diag=FALSE)] <- t(A)[lower.tri(A, diag=FALSE)] 
    
  } else {
    stop("Unknown value of cov.method.")
  }
  
  ## add variance on diagonal
  var.range <- c(0.5, 2.5)
  diag.sd <- sqrt(runif(p, min=var.range[1], max=var.range[2]))
  base.matrix <- diag(diag.sd) %*% A %*% diag(diag.sd)
  if (!isSymmetric(base.matrix)) {
    stop("Something is wrong: base.matrix is not symmetric!")
  }
  return (base.matrix)
}

#' Adjust base matrix to be positive definite
#' 
#' @param base.matrix.1 a \eqn{p} by \eqn{p} base matrix
#' @param base.matrix.2 optional, another \eqn{p} by \eqn{p} base matrix.
#' 
#' @details
#' If both \code{base.matrix.1} and \code{base.matrix.2} are provided, then 
#' the function finds a postive constant \code{lambda}, such that 
#' \code{base.matrix.1} + \code{lambda} I 
#' and 
#' \code{base.matrix.2} + \code{lambda} I 
#' are both positive definite, with minimum eigenvalues being at least 0.05.
#' If only \code{base.matrix.1} is provided, the function finds
#' \code{base.matrix.1} + \code{lambda} I
#' such that its minimum eigenvalue is at least 0.05.
#' 
#' @return \code{lambda}, a positive constant to adjust given base matrix(es).
adj.base.matrix <- function(base.matrix.1, base.matrix.2=NULL) {
  min.eigen <- rARPACK::eigs_sym(base.matrix.1, k=1, which="SA")$values
  if (! is.null(base.matrix.2)) {
    min.eigen.2 <- rARPACK::eigs_sym(base.matrix.2, k=1, which="SA")$values
    min.eigen <- min(min.eigen, min.eigen.2)
  }
  return( abs(min.eigen) + 0.05)
}


#' Generate the differential matrix
#' 
#' @param p feature dimension
#' @param diff.method the method for generating the differential matrix. 
#'            See \strong{Details} below.
#'
#' @details The \code{diff.method} argument is a character string
#' that specifies the type of differential matrix to be generated,
#' such that \code{sigma.2} = \code{sigma.1} + \code{diff.matrix}
#' for alternative hypothesis.
#' Possible values are:
#' \tabular{ll}{
#' "block" \tab \code{diff.matrix} is non-zero in an s-by-s sub-block,
#'              with entries generated from uniform distribution.\cr
#' "spike" \tab \code{diff.matrix} is non-zero in (2s)-by-(2s) sub-block,
#'              with larger entries in an s-by-s sub-block.
#' }
#' 
#' @return A \eqn{p} by \eqn{p} symmetric differential matrix.
get.diff.matrix <- function(p, diff.method, theta=1, max.var=1, s=2) {
  if (diff.method == "block") {
    ## signal level
    delta <- theta * sqrt(max.var * log(p))
    signal <- matrix(runif(s*s, min=delta/2, max=2*delta), nrow=s, ncol=s)
    ## make it symmetric
    signal[lower.tri(signal, diag=FALSE)] <- t(signal)[lower.tri(signal, diag=FALSE)]
    
    diff.matrix <- matrix(0, nrow=p, ncol=p)
    diff.matrix[1:s, 1:s] <- signal

  } else if (diff.method == "spike") {
    ## unit vector
    v <- rep(0, p)
    v[1:(2*s)] <- rnorm(2*s, mean=0.1, sd=0.1) ## weak signals
    v[1:s] <- rnorm(s, mean=1, sd=0.1) ## strong signals
    v <- v / sqrt(sum(v^2)) ## normalize
    
    ## signal level
    delta <- theta * sqrt(max.var * log(p))
    diff.matrix <- delta * v %*% t(v)
    
  } else {
    stop("Unknown value of diff.method.")
  }
  
  return (diff.matrix)
}



