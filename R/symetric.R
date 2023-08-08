.symetric <- function(M){
	p = nrow(M)
	for(i in 1:(p-1)){
		for(j in (i+1):p){
			M[j,i] = M[i,j]
		}
	}
	M
}