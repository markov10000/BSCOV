#' Select the thresholding parameter using the  quantile of CUSUM stiatistics
#'
#' @param x input time series matrix, with each row representing a time series
#' @param dw trims off the interval of consideration in the binary segmentation algorithm and determines the minimum length of a stationary segment; if dw=NULL, a default value is chosen as described in the Appendix of Barigozzi, Cho & Fryzlewicz (2016)
#' @param SN.op if \eqn{SN.op==1} Self Normalized; \eqn{SN.op==2} mean of absolute; \eqn{SN.op==3} \eqn{mad(abs(diff(x)-mad(diff(x))))} ...
#' @param M number of re-sampling
#' @param select quantile position
#' @return quantile of CUSUM stiatistics
#' @export

select.thr <-
function(x,dw=1, SN.op=0,M=1,select=seq(0.99,1,0.001)){
	len <- dim(x)[2]
	cdist <- c()
	for (Mi in 1:M){
		dist<- abs(calcul.cusum(x[,sample(len)],dw=dw,SN.op=SN.op))
		cdist <- rbind(cdist, dist)
		}
	return(quantile(apply(cdist[,-c((1:dw), (len-dw+1):len),drop=F],1,max),select))
}
