#' Data transformataion before construction of CUSUM stiatistics
#'
#' @param raw input time series matrix, with each row representing a time series
#' @param idio.diag if(\eqn{idio.diag=TRUE}) only diagnal elements (variance) are used
#' @param input.op options of data tramsformation. if \eqn{input.op=1},  second moments, if \eqn{input.op=2}, square of difference and square of sum, if \eqn{input.op=3}, wavelet transformation
#' @return the transformed data
calcul.input <- function(raw, idio.diag=FALSE, input.op=1)
{
	n <- dim(raw)[1]
	T <- dim(raw)[2]

if(input.op==1){#covariance
	if(idio.diag){
 		vech.e <- raw^2
 	}else{
		vech.e <- matrix(0, n*(n+1)/2,T)
		qq.ind <- 1
		for (qq1 in 1:n){
			for (qq2 in 1:qq1){
				vech.e[qq.ind,] <- raw[qq1,] * raw[qq2,]
				qq.ind <- qq.ind + 1
			}
		}

	}
}
if(input.op==2){#square of plus and minus
	if(idio.diag){
 		vech.e <- raw^2 * 2
 	}else{
		vech.e <- matrix(0, n*n,T)
		vech.e[1:n,] <- raw^2 * 2
		qq.ind <- 1
		for (qq1 in 2:n){
			for (qq2 in 1:(qq1-1)){
				vech.e[n+qq.ind,] <- (raw[qq1,] + raw[qq2,])^2
				vech.e[n+qq.ind+1,] <- (raw[qq1,] - raw[qq2,])^2
				qq.ind <- qq.ind + 2
			}
		}
	}
}
if(input.op==3){#wavelet
	scales <- -(1:floor(log(log(T,2),2)))
	vech.e <- c()
	for(sc in scales){
		cc <- func_coef(raw, sc)
		if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
		if(idio.diag){
 			vech.e <- rbind(vech.e, t(func_input_on(cc)))
 		}else{
	  		sgn <- sign(cc%*%t(cc))
			vech.e <- rbind(vech.e, t(func_input(cc,sgn)))
		}
	}
}
return(vech.e)
}
