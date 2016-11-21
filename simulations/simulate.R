## 
## This is the code to re-produce the simulation results in
## Zhu, Lei, Devlin and Roeder (2016), 
## "Testing High Dimensional Covariance Matrices, 
## with Application to Detecting Schizophrenia Risk Genes", 
## (arXiv:1606.00252)
##
## For more detailed description of the simulation scenarios,
## please refer to Section 3 of the paper.
##
## Please note that the simulation can take hours.
## 

rm(list=ls())
library(sLED)
library(MASS)
source("simulation_utils.R")

#######################
## Parameters
#######################

## The distribution of samples.
## Possible values are: 
## 1. "normal": Z_ij are sampled from standard normal
## 2. "nb": Z_ij are sampled from centralized NegativeBinomial(2, 2) 
## 3. "t": Z_ij are sampled from t-distribution with df=12
## 4. "gamma": Z_ij are sampled from centralized Gamma(4, 0.5)
sample.dist <- "normal"

## Covariance structure for $\Sigma_1$. 
## Possible values are: 
## 1. "noisy.diag": noisy diagonal
## 2. "block.diag": block-diagonal
## 3. "exp.decay": exponential decay
## The method "WGCNA" in the paper is not provided here 
## because of the restricted access of the CMC data.
cov.method <- "noisy.diag"

## Differential matrix structure. 
## Possible values are: 
## 1. "block": sparse-block difference
## 2. "spike": soft-sparse spiked difference
diff.method <- "block"

## Signal strength:
## the magnitude of entries in the differential matrix is 
## d = \theta * \sqrt{ max_j \Sigma_{jj}^* \log(p) } 
if (diff.method == "block") {
  theta <- 0.5
} else if (diff.method == "spike") {
  theta <- 4
}

## Other parameters
n1 <- 100 ## sample size in the first population
n2 <- 100 ## sample size in the second population 
p <- 100  ## feature dimension
s <- floor(p*0.1) ## signal sparsity

nrep <- 100 ## numer of repetitions to estimate empirical power
npermute <- 100 ## number of permutations to compute each p-value
nproj <- 100 ## number of projections to use for "RProj" (Wu and Li (2015))

useMC <- TRUE ## use multi-core parallelization to speed up computation
mc.cores <- 2 ## number of cores to use for parallelization

## Specify a sequence of smoothing parameters for sLED.
## This is the parameter $c$ in equation (2.16).
if (p >= 100) {
  sumabs.seq <- seq(0.1, 0.3, by=0.02)
} else {
  ## Note that the parameter $c$ must satisfy
  ## $c \sqrt{p} >= 1$, so $c=0.1$ is not valid when $p<100$.
  stop("p is too small for c=0.1.")
}


#######################
## Generate covariance matrices
#######################
cov.matrices <- get.cov.matrix(p=p, cov.method=cov.method,
            diff.method=diff.method, theta=theta, s=s)

#######################
## Simulations
#######################
## store the p-values in each repetition
Cai.alt <- rep(NA, nrep)
Chang.alt <- rep(NA, nrep)
LC.alt <- rep(NA, nrep)
Schott.alt <- rep(NA, nrep)
WL.alt <- rep(NA, nrep)
## for sLED, a sequence of smoothing parameters is used
sLED.alt <- matrix(NA, nrow=nrep, ncol=length(sumabs.seq))

## obtain p-values for all tests
for (rep in c(1:nrep)) {
  ## track progress
  cat(rep, "/", nrep, "\n")
  
  ## simulate samples
  samples <- simulate.samples(n1=n1, n2=n2, p=p, gamma.1=cov.matrices$gamma.1,
                              gamma.2=cov.matrices$gamma.2, sample.dist=sample.dist)
  X <- samples$X
  Y.alt <- samples$Y.alt

  ## sLED
  sLED.alt[rep, ] <- sLED(X, Y.alt, sumabs.seq=sumabs.seq, npermute=npermute, useMC=useMC, mc.cores=mc.cores)$pVal

  ## Cai et al. (2013)
  Cai.alt[rep] <- getPermutePval(X, Y.alt, Cai.max.test, npermute=npermute, useMC=useMC, mc.cores=mc.cores)$pVal

  ## Chang et al. (2016)
  Chang.alt[rep] <- Chang.maxBoot.test(X, Y.alt, nresample=npermute, useMC=useMC, mc.cores=mc.cores)$pVal

  ## Li and Chen (2012)
  LC.alt[rep] <- getPermutePval(X, Y.alt, LC.U.test, npermute=npermute, useMC=useMC, mc.cores=mc.cores)$pVal

  ## Schott (2007)
  Schott.alt[rep] <- getPermutePval(X, Y.alt, Schott.Frob.test, npermute=npermute, useMC=useMC, mc.cores=mc.cores)$pVal

  ## Wu and Li (2015)
  ## Note that this method can be very slow if \code{nproj} is larger.
  WL.alt[rep] <- getPermutePval(X, Y.alt, 
                                    testFunction=function(x, y){
                                      WL.randProj.test(x, y, nproj=nproj)
                                    },
                                    npermute=npermute, useMC=useMC, mc.cores=mc.cores)$pVal
}

testing.results <- list(sLED.alt=sLED.alt, Cai.alt=Cai.alt,
                        Chang.alt=Chang.alt, LC.alt=LC.alt,
                        Schott.alt=Schott.alt, WL.alt=WL.alt)


########################
## Compute empirical power
########################
level <- 0.05 ## level of test
i.sled.pick <- 10 ## use c=0.3 for sLED
summary.results <- data.frame(method=names(testing.results),
                              power=sapply(testing.results,
                                          function(x){
                                            if (!is.null(dim(x))) {
                                              return(sum(x[, i.sled.pick] <= level) / nrep)
                                            } else {
                                              return(sum(x <= level) / nrep)
                                            }
                                          }))
print(summary.results)
