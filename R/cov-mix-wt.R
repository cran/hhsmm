#' weighted covariance 
#'
#' The weighted means and variances using the 
#' observation matrix and the estimated weight vectors
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat, \email{aftbayat@@gmail.com}
#'
#' @param x the observation matrix
#' @param wt1 the state probabilities matrix (number of observations 
#' times number of states)
#' @param wt2 the mixture components probabilities list (of length 
#' nstate) of matrices (number of observations times number of 
#' mixture components)
#' @param cor logical. if TRUE the weighted correlation is also given
#' @param center logical. if TRUE the weighted mean is also given
#' @param method with two possible entries:
#' \itemize{
#' \item \code{"unbiased"}{ the unbiased estimator is given}
#' \item \code{"ML"}{ the maximum likelihood estimator is given}
#' }
#'
#' @return list of emission (mixture multivariate normal) parameters:
#' (\code{mu}, \code{sigma} and \code{mix.p})
#'
#' @examples
#' data(CMAPSS)
#' n = nrow(CMAPSS$train$x)
#' wt1 = matrix(runif(3*n),nrow=n,ncol=3)
#' wt2 = list()
#' for(j in 1:3) wt2[[j]] = matrix(runif(5*n),nrow=n,ncol=5)
#' emission = mixmvnorm_mstep(CMAPSS$train$x, wt1, wt2)
#' 
#' @export
#'
cov.mix.wt<- function (x, wt1 = rep(1/nrow(x), nrow(x)), wt2 = rep(1/nrow(x), nrow(x)) , 
 cor = FALSE, center = TRUE, method = c("unbiased", "ML")) 
{
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    else if (!is.matrix(x)) 
        stop("'x' must be a matrix or a data frame")
    if (!all(is.finite(x))) 
        stop("'x' must contain finite values only")
    n <- nrow(x)
    if (length(wt1) != n || length(wt2) != n) 
        stop("length of 'wt's must equal the number of rows in 'x'")
    if (any(wt1 < 0) || any(wt2 < 0)) 
        stop("weights must be non-negative!")
    if ((s1 <- sum(wt1)) == 0) 
        stop("state weights must be not all zero")
    if ((s2 <- sum(wt2)) == 0) 
        warning("for some mixture components weights are all zero! The components are not used!")
	wt <- wt1 * wt2 / sum(wt1 * wt2)
    if (is.logical(center)) {
        center <- if (center){ 
            colSums(wt * x) 
        }else 0
    }
    else {
        if (length(center) != ncol(x)) 
            stop("length of 'center' must equal the number of columns in 'x'")
    }
    x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)
    cov <- switch(match.arg(method), unbiased = crossprod(x)/(1 - 
        sum(wt^2)), ML = crossprod(x))
    	y <- list(cov = cov, center = center, n.obs = n)
    	y$wt1 <- wt1
    	y$wt2 <- wt2
	y$pmix = mean(wt2)
    if (cor) {
        Is <- 1/sqrt(diag(cov))
        R <- cov
        R[] <- Is * cov * rep(Is, each = nrow(cov))
        y$cor <- R
    }
    y
}
