#' pdf of the mixture of the robust emission proposed by 
#' Qin et al. (2024)
#'
#' The probability density function of the robust emission proposed by 
#' Qin et al. (2024) for a specified observation vector, a 
#' specified state and a specified model's parameters
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param x an observation vector or matrix
#' @param j a specified state between 1 to nstate
#' @param model a hhsmmspec model
#' @param control a list containing the control parameter k with the default 
#' value equal to 1.345
#' 
#' @return the probability density function value
#'
#' @examples
#' J <- 3
#' initial <- c(1, 0, 0)
#' semi <- c(FALSE, TRUE, FALSE)
#' P <- matrix(c(0.8, 0.1, 0.1, 0.5, 0, 0.5, 0.1, 0.2, 0.7), 
#' nrow = J, byrow = TRUE)
#' par <- list(mu = list(list(7, 8),list(10, 9, 11), list(12, 14)),
#' sigma = list(list(3.8, 4.9), list(4.3, 4.2, 5.4), list(4.5, 6.1)),
#' mix.p = list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.5, 0.5)))
#' sojourn <- list(shape = c(0, 3, 0), scale = c(0, 10, 0), type = "gamma")
#' model <- hhsmmspec(init = initial, transition = P, parms.emis = par,
#' dens.emis = dmixmvnorm, sojourn = sojourn, semi = semi)
#' train <- simulate(model, nsim = c(10, 8, 8, 18), seed = 1234, 
#' remission = rmixmvnorm)
#' mu = list(0,1)
#' sigma = list(1,2)
#' robustmodel = list(parms.emission = list(mu = mu,sigma = sigma))
#' p = drobust(train$x, 1, robustmodel)
#'
#' @references
#' Qin, S., Tan, Z., and Wu, Y. (2024). On robust estimation of hidden semi-Markov 
#' regime-switching models. Annals of Operations Research, 1-33.
#' 
#' @export
#' 
drobust <- function(x, j, model, control = list(k = 1.345))
{
	k = control$k
	n = nrow(x)
	d = ncol(x)
	mu = model$parms.emission$mu[[j]]
	sigma = model$parms.emission$sigma[[j]]
	dens <- numeric(n)
  	for(i in 1:n){
    	index <- (abs((x[i,] - mu)/sigma) <= k) * 1
    	dens[i] = prod(1/sigma * exp(-sum(0.5*((x[i,] - mu)/sigma)^2 * index - (k * (abs((x[i,] - mu)/sigma) - 0.5*k)) * (index - 1))))
  	}
	dens
}
