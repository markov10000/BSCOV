#' Using Wild Sparsified Binary Segmentation to find chang points in the second order structure in both factors and idiosyncratic errors.
#'
#' @param y input time series matrix, with each row representing a time series
#' @param dw trims off the interval of consideration in the binary segmentation algorithm and determines the minimum length of a stationary segment; if dw=NULL, a default value is chosen as described in the Appendix of Barigozzi, Cho & Fryzlewicz (2016)
#' @param screen.id id given by screening
#' @param mathrnd random intervals
#' @param SN.op if \eqn{SN.op==1} Self Normalized; \eqn{SN.op==2} mean of absolute; \eqn{SN.op==3} \eqn{mad(abs(diff(x)-mad(diff(x))))} ...
#' @param do.parallel if \eqn{do.parallel=TRUE}, a set of copies of R running in parallel are created and used for bootstrap procedure
#' @return matrix with first row MaxCUSUM, second row location of MaxCUSUM
#' @author Yu-Ning Li \email{yl3409@york.ac.uk}

calcul.matrnd <- function(y,dw,screen.id,matrnd,normf,SN.op,do.parallel){
	n <- dim(y)[1]
    T <- dim(y)[2]
    M <- dim(matrnd)[2]


	if(do.parallel&&length(screen.id)>500){#do.parallel only if dimension is big
		calcul.cusum=calcul.cusum
		mat34 <- foreach::foreach(MM=(1:M), .combine=cbind, .packages=c()) %dopar% {
			len <- matrnd[2,MM]-matrnd[1,MM] + 1
			stat<-calcul.cusum(y[screen.id, matrnd[1,MM]:matrnd[2,MM],drop=F],dw=dw,SN.op=SN.op)
			norm.stat <- normf(stat[,-c((1:dw), (len-dw+1):len,drop=F)])
			max.stat <- max(norm.stat)
			hat.chp <- matrnd[1,MM] + dw + min(which(norm.stat==max.stat)) -1
			c(max.stat,hat.chp)
		}
	}else{
		mat34 <- matrix(0,2,M)
		for (MM in 1:M){
			len <- matrnd[2,MM]-matrnd[1,MM] + 1
			stat<- calcul.cusum(y[screen.id, matrnd[1,MM]:matrnd[2,MM],drop=F],dw=dw,SN.op=SN.op)
			norm.stat <- normf(stat[,-c((1:dw), (len-dw+1):len),drop=F])
			mat34[1,MM] <- max(norm.stat)
			mat34[2,MM] <- matrnd[1,MM] + dw + min(which(norm.stat==mat34[1,MM])) -1   #hat.chp
		}
	}
	return(mat34)
}
