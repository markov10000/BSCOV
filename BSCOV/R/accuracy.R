#' Calculate the estimated error of the estimate
#'
#'@param P input set
#'@param hatP estiamtes of set P
#'@param T a positive number to calculate thr
#'@param thr if \eqn{thr=NULL}, let \eqn{thr=log(T)}
#'@return the estimated error of the estimate
#' @export
accuracy<-
function(P, hatP, T, thr=NULL){
	m <- length(P)
    n <- length(hatP)
    ACU <- rep(0,m)
    if(!m|!n) return(ACU)
    if(is.null(thr)) thr <- log(T)
    XX <- matrix(rep(P,n), m, n, byrow = FALSE)
    YY <- matrix(rep(hatP,m), m, n, byrow = TRUE)
    ACU <- rowSums(abs(XX-YY)<=thr)>0
    return (ACU)
}
