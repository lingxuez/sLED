
test_that("sLED() under null hypothesis when n1=n2=50, p=100", {
  n <- 50
  p <- 100
  npermute <- 40
  seeds <- c(1:npermute)
  
  ## simulate samples: X, Y ~ Normal(0, I)  
  set.seed(99)
  X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  set.seed(42)
  Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  
  ## sLED
  result <- sLED(X=X, Y=Y, adj.beta=0, rho=1000, 
                sumabs.seq=0.2, npermute=npermute, seeds=seeds, 
                verbose=TRUE, niter=20, trace=FALSE, useMC=FALSE)
  
  ## check results
  expect_that(result$pVal, equals(0.6), info=info)
  i.supp <- which(result$leverage != 0)
  i.expect_supp <- c(27, 31, 32, 38, 57, 59, 67)
  expect_that(i.supp, equals(i.expect_supp), info=info)
})


test_that("sLED() under null hypothesis when n1=n2=50, p=100, using multi-core", {
  n <- 50
  p <- 100
  npermute <- 40
  seeds <- c(1:npermute)
  
  ## simulate samples: X, Y ~ Normal(0, I)  
  set.seed(99)
  X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  set.seed(42)
  Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  
  ## sLED
  result <- sLED(X=X, Y=Y, adj.beta=0, rho=1000, 
                 sumabs.seq=c(0.15, 0.2), npermute=npermute, seeds=seeds, 
                 verbose=FALSE, niter=20, trace=FALSE, 
                 useMC=TRUE, mc.cores=2)
  
  ## check results
  expect_that(result$pVal, equals(c(0.9, 0.6)), info=info)
  num.supp <- rowSums(result$leverage != 0)
  num.expect_supp <- c(11, 7)
  expect_that(num.supp, equals(num.expect_supp), info=info)
})