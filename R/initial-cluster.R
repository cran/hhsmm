#' initial clustering of the data set 
#'
#' Provides an initial clustering for a data of class \code{"hhsmmdata"} which 
#' determines the initial states and mixture components (if necessary) 
#' to be used for initial parameter and model estimation
#'
#' @author Morteza Amini, \email{morteza.amini@@ut.ac.ir}, Afarin Bayat,  \email{aftbayat@@gmail.com}
#'
#' @param train the train data set of class \code{"hhsmmdata"}
#' @param nstate number of states 
#' @param nmix number of mixture components which is of one of the following forms:
#' \itemize{
#' \item a vector of positive (non-zero) integers of length \code{nstate}
#' \item a positive (non-zero) integer 
#' \item the text \code{"auto"}: the number of mixture components will be determined 
#' automatically based on the within cluster sum of squares 
#' }
#' @param ltr logical. if TRUE a left to right hidden hybrid Markov/semi-Markov model is assumed
#' @param final.absorb logical. if TRUE the final state of the sequence is assumed to be the absorbance state
#' @param verbose logical. if TRUE the outputs will be printed 
#'
#' @return a list of class \code{"hhsmm.clust"} containing the following items:
#' \itemize{
#' \item \code{clust.X}{ a list of clustered observations for each sequence and state}
#' \item \code{mix.clus}{ a list of the clusters for the mixtures for each state}
#' \item \code{state.clus}{ the exact state clusters of each observation (available if \code{ltr}=FALSE)}
#' \item \code{nmix}{ the number of mixture components (a vector of positive (non-zero) integers of length \code{nstate})}
#' \item \code{ltr}{ logical. if TRUE a left to right hidden hybrid Markov/semi-Markov model is assumed}
#' \item \code{final.absorb}{ logical. if TRUE the final state of the sequence is assumed to be the absorbance state}
#' }
#'
#' @details In reliability applications, the hhsmm models are often left-to-right
#' and the modeling aims to predict the future states. In such cases, the
#' \code{ltr}=TRUE and \code{final.absorb}=TRUE should be set. 
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
#' test <-  simulate(model, nsim = c(7,3,3,8), seed = 1234, remission = rmixmvnorm)
#' clus = initial_cluster(train,nstate=3,nmix=c(2,2,2),ltr=FALSE,
#' final.absorb=FALSE,verbose=TRUE)
#'
#' @export
#'
initial_cluster<-function(train,nstate,nmix,ltr=FALSE,final.absorb=FALSE,verbose=FALSE){
		if(length(nmix)==1 & mode(nmix)=="numeric") nmix = rep(nmix,nstate)
		if(length(nmix)!=nstate & mode(nmix)=="numeric") stop("length of nmix must be 1 or equal the number of states.")
		if(class(train)!="hhsmmdata") stop("class of train data must be hhsmmdata !")
		if(!ltr & final.absorb){
			final.absorb = FALSE
			warning("The final.absorb is for left to right case only! Changed to FALSE...")
		}
		Tx = list()
		data = as.matrix(train$x)
		num.units= length(train$N)
		for(j in 1:nstate){
			Tx[[j]] = matrix(0,nrow = 1,ncol=ncol(data))
		}# for j
		Ns = c(0,cumsum(train$N))
		xt = list()
		if(verbose) cat("Within sequence clustering ... \n")
		if(ltr){
			clusters = NULL
		} else {
			clusters = list()
		}
		if(ltr){
			for(m in 1:num.units){
				xt[[m]] = list()
				if(verbose) .progress(x=m,max=num.units)
				C=as.matrix(data[(Ns[m]+1):Ns[m+1],])
				if(final.absorb) D= as.matrix(C[-nrow(C),]) else D = C
				clus = ltr_clus(D,nstate-final.absorb)
				for(j in 1:(nstate-final.absorb)){
  					if(sum(clus==j)>0){
						xt[[m]][[j]]=matrix(D[clus==j,],
							nrow=sum(clus==j),ncol=ncol(D))
						colnames(xt[[m]][[j]]) <- colnames(train$x)
						Tx[[j]]= rbind(Tx[[j]],xt[[m]][[j]])
						colnames(Tx[[j]]) <- colnames(train$x)
					} else{
						xt[[m]][[j]] = NA
					}#if else
				}# for j
				if(final.absorb){
					Tx[[nstate]]=rbind(Tx[[nstate]],C[nrow(C),])
					colnames(Tx[[nstate]]) <- colnames(train$x)
					xt[[m]][[nstate]]=matrix(C[nrow(C),],
							nrow=1,ncol=ncol(D))	
					colnames(xt[[m]][[nstate]]) <- colnames(train$x)
				}
			}# for m 
		} else{
				clus = kmeans(data,nstate,nstart=10)$cluster
				for(m in 1:num.units){
					clusters[[m]] = clus[(Ns[m]+1):Ns[m+1]]
					xt[[m]] = list()
					D = as.matrix(data[(Ns[m]+1):Ns[m+1],])
					for(j in 1:nstate){
						xt[[m]][[j]]=matrix(D[clusters[[m]]==j,],
							nrow=sum(clusters[[m]]==j),ncol=ncol(data))
						Tx[[j]]= rbind(Tx[[j]],xt[[m]][[j]])
						colnames(Tx[[j]]) <- colnames(train$x)
					}# for j
				}# for m 
		}
		for(j in 1:nstate) Tx[[j]] = as.matrix(Tx[[j]][-1,])
		anmix = c()
		mix.clus = list()
		for(j in 1:nstate){
			if(verbose) cat("State ",j,"\n")
			if(verbose) cat("Between sequence clustering ... \n")
			if(length(nmix)==1){ if(nmix=="auto"){
				if(verbose) cat("Automatic determination of the number of mixture components ... \n")
				continue = TRUE
				DW = Inf
				oldW = (nrow(Tx[[j]])-1)*sum(apply(Tx[[j]],2,var))
				cntr = 0
				eps = 1e-2
				anmix[j] = 1
				while(continue & (cntr+1) < (nrow(Tx[[j]])*0.5) ){
					cntr = cntr + 1
					tmpclus = kmeans(Tx[[j]],cntr+1,nstart=10)
					newW = sum(tmpclus$withinss)
					DW = c(DW,oldW - newW)
					oldW = newW
					if(cntr>2){
						DDW = -diff(DW)/DW[-1]
						DDDW = -diff(DDW)
						if(any(DDDW<=0) | cntr > 10){
							anmix[j] = which.max(DDW[-1]) + 2
							continue = FALSE
						}# if 
					}# if 
				}# while
				mix.clus[[j]] = kmeans(Tx[[j]],anmix[j],nstart=10)$cluster	
			}# if 
			} else {
				mix.clus[[j]] = kmeans(Tx[[j]],nmix[j],nstart=10)$cluster	
				anmix[j] = nmix[j]
			}# if else 
		}# for j
	out = list(clust.X=xt, mix.clus=mix.clus, state.clus = clusters, 
		nmix=anmix, ltr = ltr, final.absorb = final.absorb)
	class(out) <- "hhsmm.clust"
	return(out)
}