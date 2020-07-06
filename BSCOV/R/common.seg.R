common.seg <- function(gfm, q, dw=NULL, rule=NULL, SN.op=0, sig.lev=1,
	common.norm=2, common.norm.thr=0, norm.op=FALSE, SSIC.slack=c(1/2,1/2),
	ssic=0, WBS=0, M=0, input.op=1, do.parallel=FLASE,setseed=0){
	cusum.op=1 ##fixed 
	lam <- gfm$lam; f <- gfm$f; nx <- gfm$norm.x
	f <- f[1:q, , drop=FALSE] 
	# hat.chi <- lam[, 1:q, drop=FALSE]%*%f[1:q, , drop=FALSE]
	T <- dim(nx)[2]
	if (is.null(dw)) dw <- round(min(log(T)^2, T^(6/7)/4))

	
	vech.f <- calcul.input(f, idio.diag=FALSE, input.op=input.op)
	d <- nrow(vech.f); len <- ncol(vech.f)
	if(is.null(rule)) rule <- round(log(len, 2)/2)

	mt <- make.tree(vech.f, dw, rule, SN.op=SN.op, WBS=WBS,SBS=0, M=M, norm0=common.norm, norm.thr=common.norm.thr, norm.op=norm.op,
		do.parallel=do.parallel, setseed=setseed)
	mat <- mt$mat
	matrnd <- mt$matrnd
	ct.pre <- SSIC(vech.f,mat,SSIC.slack=SSIC.slack,SIC.const=NULL)
	ct <- c()
	ct <- ct.pre$w.cpt1
	cq <- select.thr(vech.f - means.between.cpt(vech.f, ct$cpt.ic), dw, SN.op, select= c(sig.lev,seq(0.99,1,0.001)))

	ls <- list( BSmat=mat, matrnd=matrnd, ct=ct.pre, est.cps=ct$cpt.ic, norm=common.norm,cus.qtl=cq[-1], cus.thr=cq[1])

	return(ls)
}