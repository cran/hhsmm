#' the M step function of the EM algorithm
#'
#' The M step function of the EM algorithm for the robust emission proposed by 
#' Qin et al. (2024) using the 
#' observation matrix and the estimated weight vectors
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param x the observation matrix
#' @param wt the state probabilities matrix (number of observations 
#' times number of states)
#' @param control a list containing the control parameter k with the default 
#' value equal to 1.345
#'
#' @return list of emission parameters:
#' (\code{mu} and \code{sigma})
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- c(FALSE, TRUE, FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), nrow = J, 
#' byrow = TRUE)
#' par <- list(mu = list(list(7, 8), list(10, 9, 11), list(12, 14)),
#' sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
#' mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
#' sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
#' remission = rmixmvnorm)
#' clus = initial_cluster(train, nstate = 3, nmix = NULL ,ltr = FALSE,
#' final.absorb = FALSE, verbose = TRUE)
#' initmodel1 = initialize_model(clus = clus, sojourn = "gamma", 
#' M = max(train$N), semi = semi, dens.emission = drobust, mstep = robust_mstep)
#' # not test
#' # fit1 = hhsmmfit(x = train, model = initmodel1, M = max(train$N), 
#' # mstep = robust_mstep)
#'
#'
#' @references
#' Qin, S., Tan, Z., & Wu, Y. (2024). On robust estimation of hidden semi-Markov 
#' regime-switching models. Annals of Operations Research, 1-33.
#' 
#' @importFrom stats optim
#' 
#' @export
robust_mstep <- function(x, wt, control = list(k = 1.345)) 
{
	k = control$k
  	J = ncol(wt)
  	wt <- wt / rowSums(wt)
  	n = nrow(x)
	d = ncol(x)
	emission = list(mu = list(), sigma =  list())
  	for (j in 1:J) {
    	mloglike = function(theta){
      		mu = theta[1:d]
	  		sig = exp(theta[(d+1):(2*d)])
	  		tmu = tsigma = list()
	  		for(jj in 1:J){
	      		tmu[[jj]] = mu
    	  		tsigma[[jj]] = sig
	  		}
      		tmodel = list(parms.emission = 
                      list(mu = tmu,
						sigma = tsigma))
			z <- drobust(x, j, 
				model = tmodel, 
				control = control)
      		loglike = t(wt[, j]) %*% log(z)
      		return(-loglike) 
    	}	
   		start = rep(0,2*d)
    	fit = suppressWarnings(optim(start,mloglike))
		thetahat = fit$par 
		emission$mu[[j]] = thetahat[1:d]
		emission$sigma[[j]] = exp(thetahat[(d+1):(2*d)])
  	}# for j
  	emission
}# end of function
