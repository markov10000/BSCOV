idio.seg <- function(hat.vep, dw=NULL, rule=NULL, SN.op=0, sig.lev=1,
idio.norm=2, idio.norm.thr=Inf, norm.op=FALSE, SSIC.slack=c(1/2,1/2),
	SIC.op=1, WBS=0, SBS=1, M=0, idio.diag=FALSE, input.op=1, do.parallel=TRUE,setseed=0){
    T <- dim(hat.vep)[2]
    if (is.null(dw))  dw <- round(min(log(T)^2, T^(6/7)/4))
	
	vech.e <- calcul.input(hat.vep, idio.diag=idio.diag, input.op=input.op)
	d <- nrow(vech.e); len <- ncol(vech.e)
	if(is.null(rule)) rule <- round(log(len, 2)/2)
	mt <- make.tree(vech.e, dw, rule, SN.op=SN.op,
	 WBS=WBS, SBS=SBS, M=M, norm0=idio.norm, norm.thr=idio.norm.thr, norm.op=norm.op,
	 do.parallel=do.parallel, setseed=setseed)
	mat <- mt$mat
	matrnd <- mt$matrnd
	ct.pre <- SSIC(vech.e,mat,SSIC.slack=SSIC.slack,SIC.const=NULL)
	ct <- c()
    ct <- ct.pre$w.cpt2
	cq <- select.thr(vech.e - means.between.cpt(vech.e, ct$cpt.ic), dw, SN.op, select= c(sig.lev,seq(0.99,1,0.001)))

	ls <- list( BSmat=mat, matrnd=matrnd, ct=ct.pre, est.cps=ct$cpt.ic, norm=idio.norm,idio.norm.thr=idio.norm.thr, cus.qtl=cq[-1], cus.thr=cq[1])
		
	return(ls)
}