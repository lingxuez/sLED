
test_that("sLED() under null hypothesis when n1=n2=50, p=100", {
  n <- 50
  p <- 100
  
  ## simulate samples: X, Y ~ Normal(0, I)  
  set.seed(99)
  X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  set.seed(42)
  Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  
  ## sLED
  seeds <- c(1:10)
  result <- sLED(X=X, Y=Y, adj.beta=0, rho=1000, 
                sumabs.seq=0.2, npermute=10, seeds=seeds, 
                verbose=TRUE, niter=20, trace=FALSE, useMC=FALSE)
  
  ## check results
  expect_that(result$pVal, equals(0.8), info=info)
  i.supp <- which(result$leverage != 0)
  i.expect_supp <- c(8, 27, 46, 67, 77, 94)
  expect_that(i.supp, equals(i.expect_supp), info=info)
})


test_that("sLED() under null hypothesis when n1=n2=50, p=100, using multi-core", {
  n <- 50
  p <- 100
  
  ## simulate samples: X, Y ~ Normal(0, I)  
  set.seed(99)
  X <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  set.seed(42)
  Y <- matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)
  
  ## sLED
  seeds <- c(1:10)
  result <- sLED(X=X, Y=Y, adj.beta=0, rho=1000, 
                 sumabs.seq=0.2, npermute=10, seeds=seeds, 
                 verbose=TRUE, niter=20, trace=FALSE, 
                 useMC=TRUE, ncore=2)
  
  ## check results
  expect_that(result$pVal, equals(0.8), info=info)
  i.supp <- which(result$leverage != 0)
  i.expect_supp <- c(8, 27, 46, 67, 77, 94)
  expect_that(i.supp, equals(i.expect_supp), info=info)
})