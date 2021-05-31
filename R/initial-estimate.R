#' initial estimation of the model parameters
#'
#' Provides the initial estimates of the model parameters (mixture normal emission 
#' parameters and waiting times in each state) for an initial clustering of the class 
#' \code{"hhsmm.clust"} obtained by \code{\link{initial_cluster}}
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param clus an initial clustering of the class \code{"hhsmm.clust"} obtained by \code{initial_cluster}
#' @param verbose logical. if TRUE the outputs will be printed 
#'
#' @return a list of class \code{"hhsmm.par"} containing the following items:
#' \itemize{
#' \item \code{mu}{ list of mean vectors of mixture normal for each state and each mixture component}
#' \item \code{sigma}{ list of covariance matrices of mixture normal for each state and each mixture component}
#' \item \code{mix.p}{ list of mixture probabilities for each state}
#' \item \code{leng}{ list of waiting times in each state for each sequence}
#' \item \code{clusters}{ the exact clusters of each observation (available if \code{ltr}=FALSE)}
#' \item \code{nmix}{ the number of mixture components (a vector of positive (non-zero) integers of length \code{nstate})}
#' \item \code{ltr}{ logical. if TRUE a left to right hidden hybrid Markovian/semi-Markovianmodel is assumed}
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
#' par = initial_estimate(clus,verbose=TRUE)
#'
#' @export
#'
initial_estimate<-function(clus,verbose=FALSE){
		if(class(clus)!="hhsmm.clust") stop("an object of class hhsmm.clust is required !")
		mix.clus = clus$mix.clus
		state.clus = clus$state.clus
		ltr = clus$ltr
		final.absorb = clus$final.absorb
		nmix = clus$nmix
		clust.X = clus$clust.X
		num.units = length(clust.X)
		p = ncol(clust.X[[1]][[1]])
		nstate = length(nmix)
		if(verbose) cat("Intitial estimation .... \n")
		mu = sig = leng = list()
		mix.p = list()
		Tx = list()
		for(j in 1:nstate){
			Tx[[j]]= clust.X[[1]][[j]]
			if(num.units>1) for(m in 2:num.units) Tx[[j]]= rbind(Tx[[j]],clust.X[[m]][[j]])
			Tx[[j]] = as.matrix(Tx[[j]][apply(Tx[[j]],1,function(xx) !any(is.na(xx))),])
			if(verbose) cat("State ",j," estimation \n")
			if(nmix[j]==1){
				mu[[j]]=rep(0,p)
				sig[[j]]=matrix(0,p,p)
				mix.p[[j]] = 1
				leng[[j]]= 0
			}else{ 
				mix.p[[j]] = rep(0,nmix[j])
				mu[[j]] = sig[[j]] = list()
				leng[[j]]= 0
				for(k in 1:nmix[j]){
					mu[[j]][[k]]=rep(0,p)
					sig[[j]][[k]]=matrix(0,p,p)
				}# for k
			}# if else 
			for(m in 1:num.units){
				if(j < nstate | !final.absorb)
					if(ltr)
						if(!all(is.na(clust.X[[m]][[j]])))	
							leng[[j]][m]=nrow(clust.X[[m]][[j]])
				if(verbose) .progress(x=m,max=num.units)
			}# for m
			leng[[j]][is.na(leng[[j]])]=0
 			if(nmix[j]>1){
				for(k in 1:nmix[j]){
					mu[[j]][[k]] = colMeans(as.matrix(Tx[[j]][mix.clus[[j]] == k,]))
	  				sig[[j]][[k]] = cov(as.matrix(Tx[[j]][mix.clus[[j]] == k,]))
					if(verbose) cat("Mixture component ",k," estimation\n")
					mix.p[[j]][k] = sum(mix.clus[[j]] == k)/length(mix.clus[[j]])
					names(mu[[j]][[k]])<-	colnames(sig[[j]][[k]])<-rownames(sig[[j]][[k]])<-colnames(clust.X[[1]][[1]])
				}# for k 
			}else{
				mu[[j]] = colMeans(Tx[[j]])
	  			sig[[j]] = cov(Tx[[j]])
				names(mu[[j]])<-	colnames(sig[[j]])<-rownames(sig[[j]])<-colnames(clust.X[[1]][[1]])
			}#if else		
		}# for j
		ret = list(mu=mu, sig=sig, mix.p=mix.p, leng=leng, state.clus=state.clus, ltr=ltr, nmix=nmix)
		class(ret) <- "hhsmm.par"
	ret
}
