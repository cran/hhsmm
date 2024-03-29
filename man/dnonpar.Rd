% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnonpar.R
\name{dnonpar}
\alias{dnonpar}
\title{pdf of the mixture of B-splines for hhsmm}
\usage{
dnonpar(x, j, model, control = list(K = 5))
}
\arguments{
\item{x}{an observation vector or matrix}

\item{j}{a specified state between 1 to nstate}

\item{model}{a hhsmmspec model}

\item{control}{the parameters to control the density function.
The simillar name is chosen with that of \code{\link{nonpar_mstep}},
to be used in \code{...} argument of the \code{\link{hhsmmfit}} function.
Here, it contains only the parameter \code{K} which is the degrees of freedom for
the B-spline, default is \code{K=5}}
}
\value{
the probability density function value
}
\description{
The probability density function of a mixture of B-splines
for a specified observation vector, a specified state and a specified
model's parameters
}
\examples{
J <- 3
initial <- c(1, 0, 0)
semi <- c(FALSE, TRUE, FALSE)
P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), 
nrow = J, byrow = TRUE)
par <- list(mu = list(list(7, 8), list(10, 9, 11), list(12, 14)),
sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
remission = rmixmvnorm)
clus = initial_cluster(train, nstate = 3, nmix = NULL, ltr = FALSE,
final.absorb = FALSE, verbose = TRUE)
semi <- c(FALSE, TRUE, FALSE)
initmodel = initialize_model(clus = clus, mstep = nonpar_mstep,
	sojourn = "gamma", M = max(train$N), semi = semi)
p = dnonpar(train$x, 1, initmodel)

}
\author{
Morteza Amini, \email{morteza.amini@ut.ac.ir},
Reza Salehian,  \email{reza.salehian@ut.ac.ir}
}
