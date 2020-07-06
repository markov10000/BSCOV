
get.fac.mod <- function(x, max.q=NULL, q=NULL, bn.op=2, normalisation=FALSE){
	#library(rARPACK)
	# library(RSpectra)
	T <- ncol(x); n <- nrow(x)
	cnt <- min(n, T)
	if(is.null(max.q)) max.q <- round(sqrt(cnt))

	if(normalisation){
		mx <- matrix(rep(apply(x, 1, mean), each=T), byrow=TRUE, nrow=n)
		x <- x-mx
		sdx <- apply(x, 1, sd)
		x <- x/sdx
	} else{
		mx <- rep(0, n); sdx <- rep(1, n)
	}

	#if(!max.q) return()
	xx <- t(x)%*%x/n
	eig <- rARPACK::eigs_sym(A=xx,k=max.q) # only calculate first k eigen, using (rARPACK)
	max.q <- min(max.q,length(eig$value))
	f <- t(eig$vectors[, 1:max.q, drop=FALSE])*sqrt(T)
	lam <- x%*%t(f)/T

	if(bn.op==10){
		ic <- exp(-diff(log(eig$values)))
		q.hat <- which.max(ic)
		return(list(lam = lam, f = f, eig.val=eig$values[1:max.q], norm.x=x, mean.x=mx, sd.x=sdx, q.hat=q.hat, max.q=max.q, ic=ic))
	}

	if(bn.op){
		ic <- rep(0, 1+max.q)
		ic[1] <- (bn.op <= 4)*log(mean(x^2)) + (bn.op==5)*mean(x^2) + (bn.op >= 6)*log(mean(x^2))
		l <- 1
		while(l<=max.q){
			hchi <- lam[, 1:l, drop=FALSE]%*%f[1:l, , drop=FALSE]

			ic[l+1] <- (bn.op <= 4)*log(mean((x-hchi)^2)) + (bn.op >= 6)*log(mean((x-hchi)^2)) +
				(bn.op==1)*l*(n+T)/(n*T)*log(n*T/(n+T)) +  ##Bai IC1
				(bn.op==2)*l*(n+T)/(n*T)*log(cnt) +		   ##Bai IC2
				(bn.op==3)*l*log(cnt)/cnt + 			   ##Bai IC3
				(bn.op==4)*l*log(n*T/(n+T))/cnt +		   ##LI  IC4
				(bn.op==5)*(mean((x-hchi)^2)+l*mean((x-hchi)^2)*(n+T-l)*log(n*T)/(n*T))+
				(bn.op==6)*l*log(mean(eig$values[1:l]))/mean(eig$values[1:l]) + ##L1
				(bn.op==7)*l*log(n*sum(eig$values[1:l]))/(n*eig$values[l])   ##L2
			l <- l+1
		}
		q.hat <- which(ic==min(ic))-1
	} else{
		ic <- rep(0, max.q)
		q.hat <- q
	}

	return(list(lam = lam, f = f, eig.val=eig$values[1:max.q], norm.x=x, mean.x=mx, sd.x=sdx, q.hat=q.hat, max.q=max.q, ic=ic))
}
