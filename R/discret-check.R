#'
.discret.check <- function(data,nstate){
	a <- unique(data)
	1 * (all(trunc(a) == a) & length(a) == nstate) + 
	2 * (all(trunc(a) == a) & length(a) < nstate)
}