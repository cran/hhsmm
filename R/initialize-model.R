#' initialize the hhsmmspec model
#'
#' Initialize the \code{\link{hhsmmspec}} model by using an initial clustering of class \code{"hhsmm.clust"}  
#' obtained by \code{\link{initial_cluster}}
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param clus initial clustering of class \code{"hhsmm.clust"} obtained by \code{initial_cluster}
#' @param sojourn one of the following cases:
#' \itemize{
#' \item \code{"nonparametric"}{ non-parametric sojourn distribution}
#' \item \code{"nbinom"}{ negative binomial sojourn distribution}
#' \item \code{"logarithmic"}{ logarithmic sojourn distribution}
#' \item \code{"poisson"}{ poisson sojourn distribution}
#' \item \code{"gamma"}{ gamma sojourn distribution}
#' \item \code{"weibull"}{ weibull sojourn distribution}
#' \item \code{"lnorm"}{ log-normal sojourn distribution}
#' \item \code{"auto"}{automatic determination of the sojourn distribution using the chi-square test}
#' }
#' @param semi logical and of one of the following forms:
#' \itemize{
#' \item a logical value: if TRUE all states are considered as semi-Markovian else Markovian
#' \item a logical vector of length nstate: the TRUE associated states are considered as semi-Markovian
#' and FALSE associated states are considered as Markovian
#' \item \code{NULL}{ if \code{ltr}=TRUE then \code{semi = c(rep(TRUE,nstate-1),FALSE)}, else 
#' \code{semi = rep(TRUE,nstate)}}
#' }
#' @param M maximum number of waiting times in each state
#' @param verbose logical. if TRUE the outputs will be printed
#'
#' @return a \code{\link{hhsmmspec}} model containing the following items:
#' \itemize{
#' \item \code{init}{ initial probabilities of states}
#' \item \code{transition}{ transition matrix}
#' \item \code{parms.emission}{ parameters of the mixture normal emission (\code{mu}, \code{sigma}, \code{mix.p})}
#' \item \code{sojourn}{ list of sojourn time distribution parameters and its \code{type}}
#' \item \code{dens.emission}{ =\code{dmixmvnorm} function as emission probability density function}
#' \item \code{mstep}{ \code{= mixmvnorm_mstep} as the M step function of the EM algorithm}
#' \item \code{semi}{ a logical vector of length nstate with the TRUE associated states are considered as semi-Markovian}
#' }
#'
#' @examples
#' J <- 3
#' initial <- c(1,0,0)
#' semi <- c(FALSE,TRUE,FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), nrow = J, byrow=TRUE)
#' par <- list(mu = list(list(7,8),list(10,9,11),list(12,14)),
#' sigma = list(list(3.8,4.9),list(4.3,4.2,5.4),list(4.5,6.1)),
#' mix.p = list(c(0.3,0.7),c(0.2,0.3,0.5),c(0.5,0.5)))
#' sojourn <- list(shape = c(0,3,0), scale = c(0,10,0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(10,8,8,18), seed = 1234, remission = rmixmvnorm)
#' clus = initial_cluster(train,nstate=3,nmix=c(2,2,2),ltr=FALSE,
#' final.absorb=FALSE,verbose=TRUE)
#' initmodel = initialize_model(clus=clus,sojourn="gamma",M=max(train$N))
#'
#' @export
#'
initialize_model<-function(clus,sojourn=NULL,semi=NULL,M,verbose=FALSE){
	par  = initial_estimate(clus,verbose=verbose)
	if(verbose) cat("Initializing model ... \n")
	init_model= make_model(par,semi=semi,M=M,sojourn=sojourn)
	return(init_model)
}
