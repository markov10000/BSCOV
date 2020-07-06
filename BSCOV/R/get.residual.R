get.residual<-
function(gfm0=gfm0,common.est.cps=common.est.cps,refact=1,bn.op = 2){
	x<-gfm0$norm.x	
	T <- ncol(x)
	no.comm.cpt <- length(common.est.cps)
	common.est.cps <- c(0,common.est.cps,T)
	q.hat <- gfm0$q.hat

    if(no.comm.cpt&&refact==1){
        hat.vep <- x*0   
        for(ccpti in 1:(no.comm.cpt+1)){
            gfm <- get.fac.mod(x[,(1+common.est.cps[ccpti]):(common.est.cps[ccpti+1])], bn.op = bn.op, max.q = q.hat, normalisation=F)
            q <- gfm$q.hat
            hat.lam <- gfm$lam[, 1:q, drop=FALSE]; hat.f <- gfm$f[1:q, , drop=FALSE]; hat.nx <- gfm$norm.x
            hat.vep[,(1+common.est.cps[ccpti]):(common.est.cps[ccpti+1])] <- hat.nx - hat.lam%*%hat.f   
        } 
    }else{   	
        hat.lam <- gfm0$lam[, 1:q.hat, drop=FALSE]; hat.f <- gfm0$f[1:q.hat, , drop=FALSE]; hat.nx <- gfm0$norm.x
        hat.vep <- hat.nx - hat.lam%*%hat.f
    }
    return(hat.vep)
}