% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation-CLX.R
\name{Cai.max.test}
\alias{Cai.max.test}
\title{Two-sample covariance test (Cai et al. 2013)}
\usage{
Cai.max.test(X, Y)
}
\arguments{
\item{X}{n1 by p matrix, observation of the first population, columns are features}

\item{Y}{n2 by p matrix, observation of the second population, columns are features}
}
\value{
A list with the following components:
 \item{Mn}{the largest M_ij as defined in Cai (2013) equation (2)}
 \item{test.stat}{test statistic (calculated as Mn - 4*log p + log log p)}
 \item{pVal}{p-value given by the limiting distribution (Gumbol distribution)}
}
\description{
Testing the equality of two high-dimensional covariance matrices
based on the \eqn{L_\infinity} norm, proposed in Cai, Liu and Xia (2013)
"Two-Sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings".
}
\references{
Cai, Liu and Xia (2013)
"Two-Sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings",
Journal of the American Statistical Association.
}
\seealso{
\code{Chang.maxBoot.test()}, \code{LC.U.test()}, \code{WL.randProj.test()}, 
\code{Schott.Frob.test()}.
}

