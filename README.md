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

and to identify the non-zero entries if the null hypothesis is rejected. sLED is more powerfull than many existing two-sample testing procedures for high-dimensional covariance matrices, even when the signal is both weak and sparse.


## Installation
This package can be installed through `devtools` in R as follows:
```{r}
install.packages("devtools")
library("devtools")
devtools::install("path/to/this/directory")
```
where `"path/to/this/directory"` should be repalced by the path to this directory in your computer.

Alternatively, in terminal, go to the directory that contains the package directory, and use
```
$ R CMD build sLED
* checking for file ‘sLED/DESCRIPTION’ ... OK
* preparing ‘sLED’:
* checking DESCRIPTION meta-information ... OK
* checking for LF line-endings in source and make files
* checking for empty or unneeded directories
* building ‘sLED_0.0.0.9000.tar.gz’
```
and you should have the file `sLED_0.0.0.9000.tar.gz` in the directory. Now use
```
R CMD INSTALL sLED_0.0.0.9000.tar.gz
```


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

## Test
This package is still under developement, and has only been tested on Mac OS X 10.11.6.
