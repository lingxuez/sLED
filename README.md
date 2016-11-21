# sLED: A two-sample test for high-dimensional covariance matrices

This is the R Code for
> Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Covariance Matrices, with Application to Detecting Schizophrenia Risk Genes", [arXiv:1606.00252](https://arxiv.org/abs/1606.00252).

Pease cite sLED in your publication if it helps your research:
```
@article{zhu2016testing,
    title={Testing High Dimensional Covariance Matrices, with Application to Detecting Schizophrenia Risk Genes},
    author={Zhu, Lingxue and Lei, Jing and Devlin, Bernie and Roeder, Kathryn},
    journal={arXiv preprint arXiv:1606.00252},
    year={2016}
}
```

## A short introduction to sLED
Suppose X, Y are p-dimensional random vectors independently coming from two populations.
Let `D` be the differential matrix

> D = cov(Y) - cov(X)

The goal for sLED is to test the following hypothesis:

> H_0: D=0 versus H_1: D \neq 0

and to identify the non-zero entries in `D` if the null hypothesis is rejected. sLED is more powerful than many existing two-sample testing procedures for high-dimensional covariance matrices (that is, when the dimension of features `p` is larger than the sample sizes), even when the signal is both weak and sparse.

sLED can also be used to compare other p-by-p relationship matrices, including correlation matrices and weighted adjacency matrices. 


## Installation
This package can be installed through `devtools` in R:
```{r}
install.packages("devtools") ## if not installed
library("devtools")
devtools::install_github("lingxuez/sLED")
```


## Examples
### Size of sLED
First, let's try sLED under the null hypothesis. We generate 100 samples from standard Normal distributions with p=100:
```{r}
n <- 50
p <- 100
set.seed(99)
X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
set.seed(42)
Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
```

Now we apply sLED. For illustration, we use 50 permutations here and leave all other arguments as default:
```{r}
library("sLED")
result <- sLED(X=X, Y=Y, npermute=50)
```

Let's check the p-value of the test, which hopefully is not too small (since `X` and `Y` are identically distributed):
```{r}
result$pVal
## [1] 0.88
```

### Power of sLED
Now let's try a more interesting example to show the power of sLED. Let's generate another 50 Gaussian samples with a different covariance structure:
```{r}
n <- 50
p <- 100
  
## The first population is still standard normal
set.seed(99)
X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  
## For the second population, the first 10 genes have different correlation structure
s <- 10
sigma.2 <- diag(p)
sigma.2[1:s, 1:s] <- sigma.2[1:s, 1:s] + 0.2
  
set.seed(42)
Y2 <- MASS::mvrnorm(n, mu=rep(0, p), Sigma=sigma.2)
```

Now we run sLED. Note that the changes in covariance matrices happen at 10% of the genes, so the ideal sparsity parameter `sumabs` should be around (usually slightly less than) `sqrt(0.1)=0.32`. Here, we pick `sumabs=0.25`, and use 100 permutations:
```{r}
## for reproducibility, let's also set the seeds for permutation
result <- sLED(X=X, Y=Y2, sumabs.seq=0.25, npermute=100, seeds = c(1:100))
result$pVal
## [1] 0
```
The p-value is near zero. Further more, let's check which genes are detected by sLED, that is, have non-zero leverage:
```{r}
which(result$leverage != 0)
## [1]  1  2  4  5  7  8  9 10 30
```
We see that sLED correctly identifies most of the first 10 signaling genes!

We can also run sLED across a range of sparsity parameters `sumabs` at once:
```{r}
## here we let sumabs.seq to be a vector of 3 different sparsity parameters
result <- sLED(X=X, Y=Y2, sumabs.seq=c(0.2, 0.25, 0.3), npermute=100, seeds = c(1:100))
                 
## we can check the 3 p-values, which are all zero
result$pVal
## [1] 0 0 0

## let's also look at which genes have non-zero leverage
detected.genes <- apply(result$leverage, 1, function(x){which(x!=0)})
names(detected.genes) <- paste0("sumabs=",result$sumabs.seq)
detected.genes
## $`sumabs=0.2`
## [1]  2  4  5  7  9 10
##
## $`sumabs=0.25`
## [1]  1  2  4  5  7  8  9 10 30
##
## $`sumabs=0.3`
## [1]  1  2  3  4  5  6  7  8  9 10 20 24 30 39 50 56 60 88 99
```
As shown above, the sLED p-value is usually pretty stable for a reasonable range of sparsity parameters. With larger `sumabs`, the solution becomes denser, i.e. more genes will have non-zero leverage.

## Parallelization

The test can get computationally expensive with large number of permutations. For big data sets, the multi-core version is recommended:
```{r}
result_multicore <- sLED(X=X, Y=Y, npermute=1000, useMC=TRUE, mc.cores=2)
```
Please note that you need to have R packages `doParallel` and `parallel` installed.


## Simulations
We provide the code to re-produce the simulation results in Section 3 of the paper. Please refer to the R script

> simulations/simulate.R


## Tests
This package is still under developement, and has only been tested on Mac OS X 10.11.6.
