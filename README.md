# sLED: A two-sample test for high-dimensional covariance matrices

This is the R Code for
> Zhu, Lei, Devlin and Roeder (2017) "Testing high-dimensional covariance matrices, with application to detecting schizophrenia risk genes", *Annals of Applied Statistics*, **11**(3):1810-1831.  ([arxiv](https://arxiv.org/abs/1606.00252))

Pease cite sLED in your publication if it helps your research:
```
@ARTICLE{zhu2017testing,
    AUTHOR = {Lingxue Zhu and Jing Lei and Bernie Devlin and Kathryn Roeder},
     TITLE = {Testing high-dimensional covariance matrices, with application to detecting schizophrenia risk genes},
   JOURNAL = {Ann. Appl. Statist.},
  FJOURNAL = {Annals of Applied Statistics},
      YEAR = {2017},
    VOLUME = {11},
    NUMBER = {3},
     PAGES = {1810-1831},
}
```

## A short introduction to sLED
Suppose `X`, `Y` are `p`-dimensional random vectors independently coming from two populations.
Let `D` be the differential matrix

> D = cov(Y) - cov(X)

The goal for sLED is to test the following hypothesis:

> H_0: D = 0 versus H_1: D != 0

and to identify the non-zero entries in `D` if the null hypothesis is rejected. sLED is more powerful than many existing two-sample testing procedures for high-dimensional covariance matrices (that is, when the dimension of features `p` is larger than the sample sizes), even when the signal is both weak and sparse.

sLED can also be used to compare other `p`-by-`p` relationship matrices, including correlation matrices and weighted adjacency matrices. 


## Installation
This package can be installed through `devtools` in R:
```{r}
install.packages("devtools") ## if not installed
library("devtools")
devtools::install_github("lingxuez/sLED")
```

## Examples
### Size of sLED
First, let's try sLED under the null hypothesis. We generate 100 samples from standard Normal distributions with `p=100`:
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

Let's check the p-value of sLED, which hopefully is not too small (since `X` and `Y` are identically distributed):
```{r}
result$pVal
## [1] 0.86
```

### Power of sLED
Now let's try a more interesting example to show the power of sLED. Let's generate another 50 Gaussian samples with a different covariance structure:
```{r}
n <- 50
p <- 100
  
## The first population is still standard normal
set.seed(99)
X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  
## For the second population, the first 10 features have different correlation structure
s <- 10
sigma.2 <- diag(p)
sigma.2[1:s, 1:s] <- sigma.2[1:s, 1:s] + 0.2
  
set.seed(42)
Y2 <- MASS::mvrnorm(n, mu=rep(0, p), Sigma=sigma.2)
```

Now we run sLED. Note that the changes in covariance matrices occur at 10% of the features, so the ideal sparsity parameter `sumabs` should be around (usually slightly less than) `sqrt(0.1)=0.32`. Here, we pick `sumabs=0.25`, and use 100 permutations:
```{r}
## for reproducibility, let's also set the seeds for permutation
result <- sLED(X=X, Y=Y2, sumabs.seq=0.25, npermute=100, seeds = c(1:100))
result$pVal
## [1] 0.03
```
The p-value is near zero. Further more, let's check which features are detected by sLED, that is, have non-zero leverage:
```{r}
which(result$leverage != 0)
## [1]  1  2  3  4  6  7  9 10 16 26 46 50 52 83
```
We see that sLED correctly identifies most of the first 10 signaling features.
We can also run sLED across a range of sparsity parameters `sumabs` at once:
```{r}
## here we let sumabs.seq to be a vector of 3 different sparsity parameters
result <- sLED(X=X, Y=Y2, sumabs.seq=c(0.2, 0.25, 0.3), npermute=100, seeds = c(1:100))
                 
## we can check the 3 p-values, which are all < 0.05
result$pVal
## [1] 0.04 0.03 0.02

## let's also look at which features have non-zero leverage
detected.genes <- apply(result$leverage, 1, function(x){which(x!=0)})
names(detected.genes) <- paste0("sumabs=",result$sumabs.seq)
detected.genes
## $`sumabs=0.2`
## [1]  2  6  7  9 16 26 50 52
## 
## $`sumabs=0.25`
## [1]  1  2  3  4  6  7  9 10 16 26 46 50 52 83
## 
## $`sumabs=0.3`
## [1]  1  2  3  4  6  7  9 10 16 18 19 20 26 35 45 46 50 52 60 63 83
```
As shown above, the sLED p-value is usually pretty stable for a reasonable range of sparsity parameters. With larger `sumabs`, the solution becomes denser, i.e. more features will have non-zero leverage.


## Sparsity parameter

A key tuning parameter for sLED is the sparsity parameter, `sumabs`. This corresponds to the parameter `c` in equation (2.16) in our paper. Larger values of `sumabs` correspond to denser solutions. Roughly speaking, `sumabs^2` provides a loose lowerbound on the proportion of features to be detected. `sumabs` can take any value between `1/\sqrt{p}` and `1`. In practice, a smaller value is usually preferred for better interpretability.  


## Parallelization

The test can get computationally expensive with large number of permutations. For larger data sets, the multi-core version is recommended:
```{r}
result_multicore <- sLED(X=X, Y=Y, npermute=1000, useMC=TRUE, mc.cores=2)
```
Please note that this requires the R packages `doParallel` and `parallel`.


## Simulations

We provide the code to re-produce the simulation results in Section 3 of the paper. Please refer to the R script

> simulations/simulate.R

Note that the simulation can take hours to finish. Multi-core parallelization is highly recommended here.


## Tests
This package is still under developement, and has only been tested on Mac OS X 10.11.6.
