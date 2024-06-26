% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hhsmmfit.R
\name{hhsmmfit}
\alias{hhsmmfit}
\title{hhsmm model fit}
\usage{
hhsmmfit(
  x,
  model,
  mstep = NULL,
  ...,
  M = NA,
  par = list(maxit = 100, lock.transition = FALSE, lock.d = FALSE, lock.init = FALSE,
    graphical = FALSE, verbose = TRUE)
)
}
\arguments{
\item{x}{a data of class \code{"hhsmmdata"}, which can also contain missing values (NA or NaN)}

\item{model}{an initial model created  by \code{hhsmm.spec} or \code{initialize_model}}

\item{mstep}{the M step function for the EM algorithm, which also can be given in the model}

\item{...}{additional parameters for the dens.emission and mstep functions}

\item{M}{the maximum duration in each state}

\item{par}{additional list of control parameters of the \code{hhsmmfit} function
including the following items:
\itemize{
\item \code{maxit} the maximum number of iterations for the EM algorithm
\item \code{lock.transition} logical. if TRUE the transition matrix will not be updated through the EM algorithm
\item \code{lock.d} logical. if TRUE the sojourn probability matrix d will not be updated through the EM algorithm
\item \code{lock.init} logical. if TRUE the initial probability vector will not be updated through the EM algorithm
\item \code{graphical} logical. if TRUE a plot of the sojourn probabilities will be plotted through the EM algorithm
\item \code{verbose} logical. if TRUE the outputs will be printed
}}
}
\value{
a list of class \code{"hhsmm"} containing the following items:
\itemize{
\item \code{loglike} the log-likelihood of the fitted model
\item \code{AIC} the Akaike information criterion of the fitted model
\item \code{BIC} the Bayesian information criterion of the fitted model
\item \code{model} the fitted model
\item \code{estep_variables} the E step (forward-backward) probabilities of the final iteration of the EM algorithm
\item \code{M} the maximum duration in each state
\item \code{J} the number of states
\item \code{NN} the vector of sequence lengths
\item \code{f} the emission probability density function
\item \code{mstep} the M step function of the EM algorithm
\item \code{yhat} the estimated sequence of states
}
}
\description{
Fits a hidden hybrid Markov-semi-Markov model to a data of class \code{"hhsmmdata"} and using an initial
model created  by \code{\link{hhsmmspec}} or \code{\link{initialize_model}}
}
\examples{
J <- 3
initial <- c(1, 0, 0)
semi <- c(FALSE, TRUE, FALSE)
P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), nrow = J, 
byrow = TRUE)
par <- list(mu = list(list(7, 8), list(10, 9, 11), list(12, 14)),
sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
remission = rmixmvnorm)
clus = initial_cluster(train, nstate = 3, nmix = c(2 ,2, 2),ltr = FALSE,
final.absorb = FALSE, verbose = TRUE)
initmodel1 = initialize_model(clus = clus, sojourn = "gamma", 
M = max(train$N), semi = semi)
fit1 = hhsmmfit(x = train, model = initmodel1, M = max(train$N))

}
\references{
Guedon, Y. (2005). Hidden hybrid Markov/semi-Markov chains.
\emph{Computational statistics and Data analysis}, 49(3), 663-688.

OConnell, J., & Hojsgaard, S. (2011). Hidden semi Markov
models for multiple observation sequences: The mhsmm package
for R. \emph{Journal of Statistical Software}, 39(4), 1-22.
}
\author{
Morteza Amini, \email{morteza.amini@ut.ac.ir}, Afarin Bayat, \email{aftbayat@gmail.com}
}
