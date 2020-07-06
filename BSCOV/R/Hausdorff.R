#' Calculate the Hausdorff distance
#''
#'@param P input set
#'@param Q infput set
#'@param MV the distance defined for empty set
#'@return Hausdorff distance
Hausdorff <-
function(P,Q,MV=Inf){
	m <- length(P)
    n <- length(Q)
    if(!m) P=NA
    if(!n) Q=NA
	if(is.na(P)&&!is.na(Q)) return(c(0,MV))
	if(!is.na(P)&&is.na(Q)) return(c(MV,0))
	if(is.na(P)&&is.na(Q)) return(c(MV,MV))

    XX <- matrix(rep(P,n), m, n, byrow = F)
    YY <- matrix(rep(Q,m), m, n, byrow = T)
    dist <- abs(XX-YY)
    return (c(max(apply(dist,1,min)),max(apply(dist,2,min))))
}
