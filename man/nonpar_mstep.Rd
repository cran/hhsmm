% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonpar_mstep.R
\name{nonpar_mstep}
\alias{nonpar_mstep}
\title{the M step function of the EM algorithm}
\usage{
nonpar_mstep(x, wt, control = list(K = 5, lambda0 = 0.5))
}
\arguments{
\item{x}{the observation matrix}

\item{wt}{the state probabilities matrix (number of observations
times number of states)}

\item{control}{the parameters to control the M-step function.
The simillar name is chosen with that of \code{\link{dnonpar}},
to be used in \code{...} argument of the \code{\link{hhsmmfit}} function.
Here, it contains the following items:
\itemize{
\item \code{K} the degrees of freedom for the B-spline, default is \code{K=5}
\item \code{lambda0} the initial value of the smoothing parameter, default is \code{lambda0=0.5}
}}
}
\value{
list of emission (nonparametric mixture of splines) parameters:
(\code{coef})
}
\description{
The M step function of the EM algorithm for the mixture
of splines nonparametric density estimator
}
\examples{
x <- rmvnorm(100, rep(0, 2), matrix(c(4, 2, 2, 3), 2, 2))
wt <- matrix(rep(1, 100), 100, 1)
emission = nonpar_mstep(x, wt)
coef <- emission$coef[[1]]
x_axis <- seq(min(x[, 1]), max(x[, 1]), length.out = 100)
y_axis <- seq(min(x[, 2]), max(x[, 2]), length.out = 100)
f1 <- function(x, y) { 
  data = matrix(c(x, y), ncol = 2)
  tmpmodel = list(parms.emission = emission)
	 dnonpar(data, 1, tmpmodel)
}
z1 <- outer(x_axis, y_axis, f1)
f2 <- function(x, y) { 
  data = matrix(c(x, y), ncol = 2)
  dmvnorm(data, rep(0, 2), matrix(c(4, 2, 2, 3), 2, 2))
}
z2 <- outer(x_axis, y_axis, f2)
par(mfrow = c(1, 2))
persp(x_axis, y_axis, z1, theta = -60, phi = 45, col = rainbow(50))
persp(x_axis, y_axis, z2, theta = -60, phi = 45, col = rainbow(50))

}
\references{
Langrock, R., Kneib, T., Sohn, A., & DeRuiter, S. L. (2015).
Nonparametric inference in hidden Markov
models using P-splines. \emph{Biometrics}, 71(2), 520-528.
}
\author{
Morteza Amini, \email{morteza.amini@ut.ac.ir},
Reza Salehian,  \email{reza.salehian@ut.ac.ir}
}
