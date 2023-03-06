#' Random data generation from the multinomial emission distribution for hhsmm model
#'
#' Generates a vector of observations from multinomial emission distribution 
#' in a specified state and using the parameters of a specified model
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}
#'
#' @param j a specified state
#' @param model a \code{\link{hhsmmspec}} model
#'
#' @return a random vector of observations from multinomial emission distribution 
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
#' plot(train)
#' 
#'
#' @importFrom stats rmultinom
#'
#' @export
#'
rmultinomial.hhsmm <- function (j, model){ 
  	n = length(model$parms.emission$prob[[j]])
  	(1:n) %*% rmultinom(1, 1, model$parms.emission$prob[[j]]) 
}
