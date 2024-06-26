#' initialize the hhsmmspec model for a specified emission distribution
#'
#' Initialize the \code{\link{hhsmmspec}} model by using an initial clustering 
#' obtained by \code{\link{initial_cluster}} and the emission distribution 
#' characterized by mstep and dens.emission
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param clus initial clustering obtained by \code{initial_cluster}
#' @param mstep the mstep function of the EM algorithm with an style 
#' simillar to that of \code{\link{mixmvnorm_mstep}}. If NULL, the 
#' \code{\link{mixmvnorm_mstep}} is considered for the complete data set and 
#' \code{\link{miss_mixmvnorm_mstep}} is considered for the data with missing 
#' values (NA or NaN)
#' @param dens.emission the density of the emission distribution with an style simillar to that of \code{\link{dmixmvnorm}}
#' @param sojourn one of the following cases:
#' \itemize{
#' \item \code{"nonparametric"} non-parametric sojourn distribution
#' \item \code{"nbinom"} negative binomial sojourn distribution
#' \item \code{"logarithmic"} logarithmic sojourn distribution
#' \item \code{"poisson"} poisson sojourn distribution
#' \item \code{"gamma"} gamma sojourn distribution
#' \item \code{"weibull"} weibull sojourn distribution
#' \item \code{"lnorm"} log-normal sojourn distribution
#' \item \code{"auto"} automatic determination of the sojourn distribution using the chi-square test
#' }
#' @param semi logical and of one of the following forms:
#' \itemize{
#' \item a logical value: if TRUE all states are considered as semi-Markovian else Markovian
#' \item a logical vector of length nstate: the TRUE associated states are considered as semi-Markovian
#' and FALSE associated states are considered as Markovian
#' \item \code{NULL} if \code{ltr}=TRUE then \code{semi = c(rep(TRUE,nstate-1),FALSE)}, else 
#' \code{semi = rep(TRUE,nstate)}
#' }
#' @param M maximum number of waiting times in each state
#' @param verbose logical. if TRUE the outputs will be printed
#' the normal distributions will be estimated 
#' @param ... additional parameters of the \code{mstep} function
#'
#' @return a \code{\link{hhsmmspec}} model containing the following items:
#' \itemize{
#' \item \code{init} initial probabilities of states
#' \item \code{transition} transition matrix
#' \item \code{parms.emission} parameters of the mixture normal emission (\code{mu}, \code{sigma}, \code{mix.p})
#' \item \code{sojourn} list of sojourn time distribution parameters and its \code{type}
#' \item \code{dens.emission} the emission probability density function
#' \item \code{mstep} the M step function of the EM algorithm
#' \item \code{semi} a logical vector of length nstate with the TRUE associated states are considered as semi-Markovian
#' }
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
#' clus = initial_cluster(train, nstate = 3, nmix = c(2 ,2, 2),ltr = FALSE,
#' final.absorb = FALSE, verbose = TRUE)
#' initmodel = initialize_model(clus = clus, sojourn = "gamma", 
#' M = max(train$N))
#'
#' @export
#'
initialize_model <- function(clus, mstep = NULL, dens.emission = dmixmvnorm, sojourn = NULL, semi = NULL, 
	M, verbose = FALSE, ...)
{
	if(is.null(mstep)){ 
		if(clus$miss){
			mstep=miss_mixmvnorm_mstep
			par  = initial_estimate(clus,mstep,verbose=verbose,par=NULL)
		}else{
 			mstep = mixmvnorm_mstep
			par  = initial_estimate(clus,mstep,verbose=verbose)
		}
	}else{
		par  = initial_estimate(clus,mstep,verbose=verbose,...)
	}
	if(verbose) cat("Initializing model ... \n")
	init_model= make_model(par,mstep,dens.emission,semi=semi,M=M,sojourn=sojourn)
	return(init_model)
}
