#' Construction of CUSUM stiatistics
#'
#' @param x input time series matrix, with each row representing a time series
#' @param dw trims off the interval of consideration in the binary segmentation algorithm and determines the minimum length of a stationary segment; if dw=NULL, a default value is chosen as described in the Appendix of Barigozzi, Cho & Fryzlewicz (2016)
#' @param alpha CUSUM parameter
#' @param SN.op if \eqn{SN.op==1} Self Normalized; \eqn{SN.op==2} mean of absolute; \eqn{SN.op==3} \eqn{mad(abs(diff(x)-mad(diff(x))))} ...
#' @return CUSUM stiatistics




calcul.cusum <-
function(x,dw=0, alpha=0.5,SN.op=0){
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation

    leftsums <- apply(x,1,cumsum)


    t <- 1:(n-1)
    tnt <- t*(n-t)
    cusum<-c()
    # constructing CUSUM matrix
    #rightsums <- t(leftsums[n,]-t(leftsums))
    #cusum$cusum0 <- t((rightsums[t,]/(n-t) - leftsums[t,]/t)*sqrt(n)*(tnt/n^2)^alpha)
    cusum$cusum <- t((leftsums[t,]*n - t%o%leftsums[n,])*(tnt)^(alpha-1) * n^(-2*alpha+0.5))
    cusum$SN <- rep(1,p)
    if(SN.op==0){#return(cusum)
    }else if (SN.op==1){#Self Normalize
        if(dw*2>n) {dw = 0; warning("dw*2>n; set dw=0")}
        rightsums <- t(leftsums[n,]-t(leftsums))
        #cpt at cusum$SNk+1
        cusum$SNk <- dw + apply(cusum$cusum[,-c(1:dw,(n-dw+1):n)],1,which.max)
        #        std<-function(x)mad(diff(x)/sqrt(2))
        for (pii in 1:p){
            SNn <- cusum$SNk[pii]
            SNt <- 1:(SNn-1)
            SNcumsum1 <- (leftsums[SNt,pii]*SNn - SNt%o%leftsums[SNn,pii]) * SNn^(-1)
            SNt2 <- 1:(n-SNn-1)
            SNcumsum2 <- (rightsums[(n-1):(SNn+1),pii]*(n-SNn) - SNt2%o%rightsums[SNn,pii]) * (n-SNn)^(-1)
            cusum$SN[pii] <- sum(rbind(SNcumsum1^2,SNcumsum2^2))^0.5/n

        }
    }else if (SN.op==2){#chi-square
        cusum$SN <- rowMeans(abs(x))
    }else if (SN.op==3){#normal
        std<-function(x) mad(abs(diff(x)-mad(diff(x))))
        cusum$SN <- apply(x,1,std)
    }else if (SN.op==4){#normal
        std<-function(x) mad(abs(x-mad(x)))
        cusum$SN <- apply(x,1,std)
    }else if (SN.op==5){#sd
        cusum$SN <- apply(x,1,sd)
        }
    return(cusum$cusum/cusum$SN)
}
