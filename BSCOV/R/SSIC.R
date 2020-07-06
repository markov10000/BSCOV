#' Using SSIC selects the number of breaks
#'
#' @param x input time series matrix, with each row representing a time series
#' @param mat information of change points
#' @param SSIC.slack two positive constants as the thresholding parameters
#' @param SIC.const pennalty function
#' @return information of breaks with one of the thresholding parameters \itemize{
#'  \item{"pen"}{pennalty function}
#'  \item{"ic.curve"}{information criteria}
#'  \item{"dim.cpt.ic"}{dimension that has the breaks}
#'  \item{"no.cpt.ic "}{number of breaks}
#'  \item{"cpt.ic"}{location of breaks}
#' }
#' @export

SSIC <- function(x,mat,SSIC.slack=c(1/2,1/2),SIC.const=NULL){
	w.cpt1 <- c()
	w.cpt2 <- c()
	x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
	n <- nrow(x)
    T <- ncol(x)

if(is.null(SIC.const)) {
    SIC.const <- min((T)^0.5,log(T)^4.01)
}
if(is.null(mat)) {
		w.cpt1$cpt.ic<- as.numeric()
		w.cpt1$no.cpt.ic <- 0
		w.cpt1$max.diff.pen <- rep(0,n)
		return(list(w.cpt1=w.cpt1,w.cpt2=w.cpt1))
	}else{
		if(SIC.const<=0){#all candidates are return
			w.cpt$cpt.ic<- cpt.trylist
			w.cpt$no.cpt.ic <- len.cpt
			w.cpt$max.diff.pen <- rep(0,n)
			return(list(w.cpt1=w.cpt1,w.cpt2=w.cpt1))
		}else{
			SIC.const1 <- SSIC.slack[1] * SIC.const
			SIC.const2 <- SSIC.slack[2] * SIC.const
		}
	}


	cpt.trylist <- mat[3,order(mat[5,],decreasing = TRUE)]
	len.cpt <- length(cpt.trylist)
	w.cpt1$pen <- SIC.const1
	w.cpt2$pen <- SIC.const2
	#ic part

		w.cpt1$ic.curve <- w.cpt2$ic.curve <- matrix(0,n,len.cpt+1)

		if(len.cpt) for(i in len.cpt:1){
			min.log.lik <- T/2 * log(rowSums((x -  means.between.cpt(x,cpt.trylist[1:i]))^2)/T) # n*1 matrix
			w.cpt1$ic.curve[,i+1] <- min.log.lik +  i * w.cpt1$pen
			w.cpt2$ic.curve[,i+1] <- min.log.lik +  i * w.cpt2$pen
		}
		w.cpt1$ic.curve[,1] <- w.cpt2$ic.curve[,1] <- T/2 * log(apply(x,1,var))
		# if(len.cpt) w.cpt$hatpen <- max(w.cpt$ic.curve[,len.cpt]-w.cpt$ic.curve[,len.cpt+1] + w.cpt$pen)

		w.cpt1$cpt.ic <- w.cpt2$cpt.ic <- list()
	    diff.penalty1 <- t(diff(t(w.cpt1$ic.curve))<0)  ##decreasing
		diff.penalty2 <- t(diff(t(w.cpt2$ic.curve))<0)  ##decreasing
		w.cpt1$dim.cpt.ic <- colSums(matrix(diff.penalty1,nrow(diff.penalty1)))
		w.cpt2$dim.cpt.ic <- colSums(matrix(diff.penalty1,nrow(diff.penalty2)))
		w.cpt1$no.cpt.ic <- which(c(w.cpt1$dim.cpt.ic,0)==0)[1] - 1
		w.cpt2$no.cpt.ic <- which(c(w.cpt2$dim.cpt.ic,0)==0)[1] - 1
		if(w.cpt1$no.cpt.ic>0){w.cpt1$cpt.ic <- cpt.trylist[1:w.cpt1$no.cpt.ic]}else{w.cpt1$cpt.ic <- numeric(0)}
		if(w.cpt2$no.cpt.ic>0){w.cpt2$cpt.ic <- cpt.trylist[1:w.cpt2$no.cpt.ic]}else{w.cpt2$cpt.ic <- numeric(0)}

#w.cpt$ic.curve= NULL
return(list(w.cpt1=w.cpt1,w.cpt2=w.cpt2))
}
