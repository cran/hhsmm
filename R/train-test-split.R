#' Splitting the data sets to train and test 
#'
#' A function to split the train data of class \code{"hhsmmdata"}
#' to train and test subsets with an option to right trim the sequences
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat, \email{aftbayat@@gmail.com}
#'
#' @param train the train data of class \code{"hhsmmdata"}
#' @param train.ratio a number in (0,1) which determines the ratio of the train subset
#' @param trim logical. if TRUE the sequences will be right trimmed with random lengths
#'
#' @return a list containing:
#' \itemize{
#' \item\code{train}{ the randomly selected subset of train data of class \code{"hhsmmdata"}}
#' \item\code{test}{ the randomly selected subset of test data of class \code{"hhsmmdata"}. 
#' If \code{trim}=TRUE, the test list also contains the RUL (residual useful lifetime) values}
#' }
#' 
#' @details In reliability applications, the hhsmm models are often left-to-right
#' and the modeling aims to predict the future states. In such cases, the 
#' test sets are right trimmed and the prediction aims to predict the 
#' residual useful lifetime (RUL) of a new sequence. 
#'
#' @examples
#'\donttest{
#' data(CMAPSS)
#' tt = train_test_split(CMAPSS$train,train.ratio = 0.7, trim = TRUE)
#'}
#' 
#' @export
#'
train_test_split <- function(train,train.ratio = 0.7, trim = FALSE){
	if(class(train)!="hhsmmdata") stop("train must be of class hhsmmdata")
	if(train.ratio<=0 | train.ratio>=1) stop("train.ratio must be in (0,1)")
	x = train$x
	N = train$N
	if(!is.null(train$s)) s = train$s else s = NULL
	n = length(N)
	ntest = trunc((1-train.ratio) * n)
	sam = sample(1:n, ntest)
	Ns = cumsum(c(0,N))
	xtest = xtrain = rep(0,ncol(x))
	Ntest = Ntrain = c()
	if(!is.null(s)) stest = strain = c()
	for(i in 1:n){
		if(i %in% sam){
			xtest = rbind(xtest, x[(Ns[i]+1):Ns[i+1],])
			Ntest = c(Ntest,N[i])
			if(!is.null(s)) stest = c(stest,s[(Ns[i]+1):Ns[i+1]])
		} else {
			xtrain = rbind(xtrain, x[(Ns[i]+1):Ns[i+1],])
			Ntrain = c(Ntrain,N[i])
			if(!is.null(s)) strain = c(strain,s[(Ns[i]+1):Ns[i+1]])
		}
	}
	xtest = xtest[-1,]
	xtrain = xtrain[-1,]
	if(trim){
		xtest.trim = rep(0,ncol(x))
		u = runif(ntest)
		Ntrim = trunc(Ntest * u)
		Nts = cumsum(c(0,Ntest))
		for(i in 1:ntest){
			xtest.trim = rbind(xtest.trim,xtest[(Nts[i]+1):(Nts[i]+Ntrim[i]),])
		}
		xtest = xtest.trim[-1,]
	}
	if(!is.null(s)){
		train = list(x = xtrain, N = Ntrain , s = strain)
			if(trim) test = list(x = xtest, N = Ntrim , s = stest, RUL = Ntest - Ntrim)
			else test = list(x = xtest, N = Ntest , s = stest)
	} else {
		train = list(x = xtrain, N = Ntrain)
			if(trim) test = list(x = xtest, N = Ntrim , RUL = Ntest - Ntrim)
			else test = list(x = xtest, N = Ntest)
	}
	class(train) <- class(test) <- "hhsmmdata"
	list(train = train , test = test)
}
