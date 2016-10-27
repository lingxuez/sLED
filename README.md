# sLED: A two-sample test for high-dimensional differential matrices

This is the R Code for
> Zhu, Lei, Devlin and Roeder (2016), "Testing High Dimensional Differential Matrices, with Application to Detecting Schizophrenia Risk Genes", [arXiv:1606.00252](https://arxiv.org/abs/1606.00252).

Pease cite sLED in your publication if it helps your research:
```
@article{zhu2016testing,
    title={Testing High Dimensional Differential Matrices, with Application to Detecting Schizophrenia Risk Genes},
    author={Zhu, Lingxue and Lei, Jing and Devlin, Bernie and Roeder, Kathryn},
    journal={arXiv preprint arXiv:1606.00252},
    year={2016}
}
```

## A short introduction to sLED
Suppose X, Y are p-dimensional random vectors independently coming from two populations.
Let D be the differential matrix

D = A(Y) - A(X)

where A() is some p-by-p relationship matrix among features in the two populations, including the covariance matrices and correlation matrices. 

The goal for sLED is to test the following hypothesis:

H_0: D=0 versus H_1: D \neq 0

and to identify the non-zero entries if the null hypothesis is rejected. sLED is more powerful than many existing two-sample testing procedures for high-dimensional covariance matrices (that is, when p is larger than the sample sizes), even when the signal is both weak and sparse.


## Installation
This package can be installed through `devtools` in R:
```{r}
install.packages("devtools") ## if not installed
library("devtools")
devtools::install_github("lingxuez/sLED")
```

Alternatively, you can download the files, open the terminal, go to the directory that contains the package directory, and use
```
$ R CMD build sLED
* checking for file ‘sLED/DESCRIPTION’ ... OK
* preparing ‘sLED’:
* checking DESCRIPTION meta-information ... OK
* checking for LF line-endings in source and make files
* checking for empty or unneeded directories
* building ‘sLED_0.0.0.9000.tar.gz’
```
and you should have the file `sLED_0.0.0.9000.tar.gz` in the directory. Now run
```
R CMD INSTALL sLED_0.0.0.9000.tar.gz
```


## Example
First, let's try sLED under the null hypothesis. We generate 100 samples from standard Normal distributions with p=100:
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

Let's check the p-value of the test, which hopefully is not too small (since `X` and `Y` are identically distributed):
```{r}
result$pVal
## 0.8
```

... more examples to come ...


## Parallelization

The test can get computationally expensive with large number of permutations. For big data sets, the multi-core version is recommended:
```{r}
result_multicore <- sLED(X=X, Y=Y, npermute=1000, useMC=TRUE, ncore=2)
```
Please note that you need to have R packages `doParallel` and `parallel` installed.

## Tests
This package is still under developement, and has only been tested on Mac OS X 10.11.6.
