#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for estimation of the 
#' parameters of the multinomial emission distribution
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, 
#'
#' @param x the observation matrix
#' @param wt the state probabilities matrix (number of observations 
#' times number of states)
#' @param n the maximum possible level of the multinomial vector (i.e. from 1 to n) 
#' 
#' @return list of multinomial emission parameters:
#' (\code{prob})
#'
#' @examples
#'J <- 2
#'initial <- c(1, 0)
#'semi <- rep(TRUE, 2)
#'P <- matrix(c(0, 1, 1, 0), 
#' nrow = J, byrow = TRUE)
#' par <- list(prob = list(c(0.6,  0.2, 0.2),
#'                            c(0.2, 0.6,  0.2)))
#' sojourn <- list(shape = c(1, 3), scale = c(2, 10), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmultinomial.hhsmm, remission = rmultinomial.hhsmm,
#'  mstep = mstep.multinomial,sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(20, 30, 42, 50), seed = 1234, 
#' remission = rmultinomial.hhsmm)
#' clus = initial_cluster(train = train, nstate = 2, nmix = NULL,
#' ltr = FALSE, final.absorb = FALSE, verbose = TRUE)
#' initmodel = initialize_model(clus = clus, mstep = mstep.multinomial, n = 3,
#' dens.emission = dmultinomial.hhsmm, sojourn = "gamma", semi = rep(TRUE, 2), 
#' M = max(train$N),verbose = TRUE)
#' fit1 = hhsmmfit(x = train, model = initmodel, mstep = mstep.multinomial, n = 3, 
#' M = max(train$N))
#' homogeneity(fit1$yhat,train$s)
#' 
#' @export
#'
mstep.multinomial <- function (x, wt, n){ 
	k = ncol(wt)
	prob = list()
	for(i in 1:k) prob[[i]] = numeric(n)
	y = t(sapply(1:n, function(i) 1*(x == i))) 
	for (i in 1:k) prob[[i]]=sapply(1:n, function(j) weighted.mean(y[j,],wt[,i]))
	list(prob = prob)
}