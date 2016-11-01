#' Permutation p-value
#' 
#' @description Obtain the p-value using permutation test.
#' 
#' @param X n1 by p matrix of the first population
#' @param Y n2 by p matrix of the second population
#' @param testFunction a function that takes in two arguments, X and Y,
#'        and returns a list including the field \code{test.stat} for test statistic
#' @param npermute the number of permutations to conduct
#' @param useMC logical variable indicating whether to use multicore parallelization.
#'        R packages \code{parallel} and \code{doParallel} are required if set to \code{TRUE}.
#' @param mc.cores decide the number of cores to use when \code{useMC} is set to \code{TRUE}.
#' @param seeds a numeric vector with length \code{npermute}, 
#'        where \code{seeds[i]} specifies the seeding for the i-th permutation. 
#'        Set to \code{NULL} if do not want to specify.
#' 
#' @return A list containing the following components:
#'  \item{test.res}{the results of the original test, including test statistic,
#'        sometimes also include the p-values calculated using the limiting distribution}
#'  \item{test.stat.permute}{a numeric vector with length \code{permute}, 
#'          for the test statsitics in permutations}
#'  \item{pVal}{the permutation pvalue, calculated as the fraction of \code{test.stat.permute}
#'      being larger than \code{test.res$test.stat}}
#'      
#' @seealso \code{Cai.max.test()}, \code{Chang.maxBoot.test()}, 
#'        \code{LC.U.test()}, \code{WL.randProj.test()}, \code{Schott.Frob.test()}.
#'        
#' @export
getPermutePval <- function(X, Y, testFunction, npermute=100, useMC=FALSE, mc.cores=1, seeds=NULL) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  Z <- rbind(X, Y)
  
  ## test statistics on original data
  test.res <- testFunction(X=X, Y=Y)
  
  ## permutation
  if (useMC) {
    ## try to use multi-core parallelization
    hasPackage <- requireNamespace("doParallel", quietly = TRUE) && requireNamespace("parallel", quietly = TRUE)
    if (!hasPackage) {
      stop("Please install packages 'doParallel' and 'parallel' for multi-core parallelization.
           Otherwise set useMC=FALSE for non-parallel computation.")
    }
    test.stat.permute <- parallel::mclapply(1:npermute, onePermuteStat, 
                                       Z=Z, n1=n1, n2=n2, testFunction=testFunction, seeds=seeds,  
                                       mc.cores=mc.cores, mc.preschedule=TRUE)
    test.stat.permute = unlist(test.stat.permute)
    
  } else {
    ## without parallelization
    test.stat.permute <- sapply(1:npermute, onePermuteStat, 
                                Z=Z, n1=n1, n2=n2, testFunction=testFunction, seeds=seeds)
  }
  
  pVal <- sum(test.stat.permute > test.res$test.stat) / npermute
  
  return(list(test.res=test.res, test.stat.permute=test.stat.permute,
              pVal=pVal))
}


#' Test statistic in one permutation
#' 
#' @description 
#' Helper function for \code{getPermutePval()}.
#' 
#' @param i the i-th permutation, in the range of {1, ..., \code{npermute}}.
#' @param Z (n1+n2)-by-p matrix, containing the pooled samples from two populations. Columns are features.
#' @param n1 the first n1 rows in Z represent the first population
#' @param n2 the (n1+1):(n1+n2) rows in Z represent the second population
#' @param testFunction a function that takes in two arguments, X and Y,
#'        and returns a list including the field \code{test.stat} for test statistic
#' @param seeds a numeric vector with length \code{npermute}, 
#'        where \code{seeds[i]} specifies the seeding for the i-th permutation. 
#'        Set to \code{NULL} if do not want to specify.
#' 
#' @return Return the test statistic in the i-th permutation.
#' 
#' @seealso \code{getPermutePval()}.
onePermuteStat <- function(i, Z, n1, n2, testFunction, seeds=NULL) {
  
  if (!is.null(seeds)){
    set.seed(seeds[i])
  }
  i.permute <- permuteIndex(n1, n2)
  
  test.stat <- testFunction(X=Z[i.permute$i1, ], Y=Z[i.permute$i2, ])$test.stat
  return(test.stat)
}

