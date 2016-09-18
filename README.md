# sLED: A two-sample test for high-dimensional differential matrices

## sLED
Suppose X, Y are p-dimensional random vectors independently coming from two populations.
Let D be the differential matrix

D = A(Y) - A(X)

where A() is some p-by-p relationship matrix among features in the two populations, including the covariance matrices and correlation matrices. 

sLED is a powerful procedure to test the following hypothesis:

H_0: D=0 versus H_1: D \neq 0

and to identify the non-zero entries if the null hypothesis is rejected.

See more details in
> Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Differential Matrices, with Application to Detecting Schizophrenia Risk Genes", [arXiv:1606.00252](https://arxiv.org/abs/1606.00252).

## Installation
This package can be installed through `devtools` in R as follows:
```{r}
install.packages("devtools")
library("devtools")
devtools::install("/path/to/sLEDpackage")
```
Alternatively, it can be installed from command line by 
> R CMD INSTALL sLED_0.0.0.9000.tar.gz

## Example
Here is a simple example to apply sLED test, when the data is generated from null hypothesis.

First, we generate the data X, Y, each containing 50 samples independently coming from standard Normal distributions with p=100:
```{r}
  n <- 50
  p <- 100
  set.seed(99)
  X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  set.seed(42)
  Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
```

Now we apply sLED. For illustration, we use 10 permutations here and leave all other arguments as default:
```{r}
library("sLED")
result <- sLED(X=X, Y=Y, npermute=10)
```

Now you can check the p-value of the test, which hopefully is not too small:
```{r}
result$pVal
## 0.6
```

Finally, note that the test can get computationally expensive with more permutations. For large data sets, the multi-core version is recommended. This requires the R packages `doParallel` and `parallel`. Once they are installed, a multi-core sLED can be run by
```{r}
result_multicore <- sLED(X=X, Y=Y, useMC=TRUE, ncore=2)
```

