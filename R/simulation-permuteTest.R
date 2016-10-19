#' Permutation p-value
#' 
#' @description Obtain the p-value using permutation test.
#' 
#' @param X n1 by p matrix of the first population
#' @param Y n2 by p matrix of the second population
#' @param testFunction a function that takes in two arguments, X and Y,
#'        and returns a list including the field \code{test.stat} for test statistic
#' @param npermute the number of permutations to conduct 
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
getPermutePval <- function(X, Y, testFunction, npermute=100) {
  test.res <- testFunction(X=X, Y=Y)
  test.stat.permute <- getPermuteStat(X=X, Y=Y, testFunction=testFunction, npermute=npermute)
  
  pVal <- sum(test.stat.permute > test.res$test.stat) / npermute
  
  return(list(test.res=test.res, test.stat.permute=test.stat.permute,
              pVal=pVal))
}

#' Permute test statistics
#' 
#' @description Obtain the null distribution of test statistics in permutation.
#' 
#' @param X n1 by p matrix of the first population
#' @param Y n2 by p matrix of the second population
#' @param testFunction a function that takes in two arguments, X and Y,
#'        and returns a list including the field \code{test.stat} for test statistic
#' @param npermute the number of permutations to conduct 
#' @return a numeric vector with length \code{npermute} for the test statsitics in permutations.
getPermuteStat <- function(X, Y, testFunction, npermute=100) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  Z <- rbind(X, Y)
  
  test.stat.permute <- rep(NA, npermute)
  for (rep in c(1:npermute)) {
    ## permute labels in two groups
    i.permute <- sample(n1+n2, size=n1+n2, replace=FALSE)
    test.stat.permute[rep] <- testFunction(X=Z[i.permute[1:n1], ], 
                                           Y=Z[i.permute[(n1+1):(n1+n2)], ])$test.stat
  }
  
  return (test.stat.permute)
}

