% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmultinomial-hhsmm.R
\name{dmultinomial.hhsmm}
\alias{dmultinomial.hhsmm}
\title{pdf of the multinomial emission distribution for hhsmm}
\usage{
dmultinomial.hhsmm(x, j, model, n)
}
\arguments{
\item{x}{the observation vector}

\item{j}{a specified state between 1 to nstate}

\item{model}{a hhsmmspec model}

\item{n}{the maximum possible level of the multinomial vector (i.e. from 1 to n)}
}
\value{
the probability density function value
}
\description{
The probability density function of a multinomial emission distribution
for a specified observation vector, a specified state and a specified
model's parameters
}
\examples{
J <- 2
initial <- c(1, 0)
semi <- rep(TRUE, 2)
P <- matrix(c(0, 1, 1, 0), 
nrow = J, byrow = TRUE)
par <- list(prob = list(c(0.6,  0.2, 0.2),
                           c(0.2, 0.6,  0.2)))
sojourn <- list(shape = c(1, 3), scale = c(2, 10), type = "gamma")
model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
dens.emis = dmultinomial.hhsmm, remission = rmultinomial.hhsmm,
 mstep = mstep.multinomial,sojourn = sojourn, semi = semi)
train <- simulate(model, nsim = c(20, 30, 42, 50), seed = 1234, 
remission = rmultinomial.hhsmm)
clus = initial_cluster(train = train, nstate = 2, nmix = NULL,
ltr = FALSE, final.absorb = FALSE, verbose = TRUE)
initmodel = initialize_model(clus = clus, mstep = mstep.multinomial, n = 3,
dens.emission = dmultinomial.hhsmm, sojourn = "gamma", semi = rep(TRUE, 2), 
M = max(train$N),verbose = TRUE)
fit1 = hhsmmfit(x = train, model = initmodel, mstep = mstep.multinomial, n = 3, 
M = max(train$N))
homogeneity(fit1$yhat,train$s)

}
\author{
Morteza Amini, \email{morteza.amini@ut.ac.ir}
}
