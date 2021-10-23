#' Using Wild Sparsified Binary Segmentation to find chang points in the second order structure in both factors and idiosyncratic errors.
#'
#' @param x input time series matrix, with each row representing a time series
#' @param normalisation if \eqn{normalisation=TRUE}, normalized the input data
#' @param r the number of factors, if \eqn{r=NULL}, screening over a range of factor number candidates is performed as described in the paper
#' @param bn.op an index number for the information criterion-based estimator of Bai and Ng (2002) for the number of factors
#' @param sig.lev
#' @param dw trims off the interval of consideration in the binary segmentation algorithm and determines the minimum length of a stationary segment; if dw=NULL, a default value is chosen as described in the Appendix of Barigozzi, Cho & Fryzlewicz (2016)
#' @param rule the depth of a binary tree for change-point analysis, see the Appendix of Barigozzi, Cho & Fryzlewicz (2016)
#' @param SN.op if \eqn{SN.op==1} Self Normalized; \eqn{SN.op==2} mean of absolute; \eqn{SN.op==3} \eqn{mad(abs(diff(x)-mad(diff(x))))} ...
#' @param refact.op if \eqn{refact.op=TRUE}, redo the PCA on each segment after the breaks in the factors are found
#' @param common.norm norm used in the  CUSUM statistics aggregation of factors 
#' @param idio.norm norm used in the  CUSUM statistics aggregation of idiosyncratic errors
#' @param common.norm.thr threshold used in the CUSUM statistics of factors 
#' @param idio.norm.thr threshold used in the CUSUM statistics of idiosyncratic errors
#' @param norm.op if \eqn{idio.diag=TRUE} the threshold is multiplied by \eqn{((e-s+1)/T)^0.5}
#' @param SSIC.slack positive constant used in the SSIC 
#' @param max.q the maximum number of factors, if \eqn{max.q=NULL}, a default value is chosen as described in the paper
#' @param WBS Wild Binary Segmentation
#' @param SBS Sparsified Binary Segmentation
#' @param M the number of random intervals
#' @param input.op options of data tramsformation. if \eqn{input.op=1},  second moments, if \eqn{input.op=2}, square of difference and square of sum, if \eqn{input.op=3}, wavelet transformation
##' @param idio.diag if \eqn{idio.diag=TRUE}, only the diagonal wavelet-transform is employed in order to generate the panel of statistics from the idiosyncratic components
#' @param do.parallel if \eqn{do.parallel=TRUE}, a set of copies of R running in parallel are created and used for bootstrap procedure
#' @param no.proc sets the number of processes to be spawned when do.parallel=TRUE
#' @param setseed sets the ramdon seed 
#' @return change points
#' \itemize{
#'  \item{"gfm"}{information of factor model}
#'  \item{"para"}{Parameters used in the algotrithm}
#'  \item{"hat.vep"}{estimate of idiosyncratic components}
#'  \item{"common.seg.res"}{information of common component breaks}
#'  \item{"idio.seg.res"}{information of idiosyncratic component breaks}
#'  \item{"common.est.cps"}{location of common component breaks}
#'  \item{"idio.est.cps"}{location idiosyncratic component breaks}
#' }
#' @author Yu-Ning Li \email{yl3409@york.ac.uk}
BSCOV <-
function (x,normalisation=FALSE, r = NULL, bn.op = 2, sig.lev=1, dw = NULL, rule = NULL, SN.op=0, refact.op=0,
      common.norm=2, idio.norm=2, common.norm.thr=0, idio.norm.thr=Inf, norm.op=FALSE, SSIC.slack=c(1/2, 1/2),
     max.q = NULL, WBS=1, SBS=1, M=500,  input.op=1,
     idio.diag = FALSE, 
    do.parallel = TRUE, no.proc = 2, setseed=0) 
{

if(input.op==3){#wavelet
    require(Rcpp)
    require(RcppArmadillo)
     sourceCpp("input.cpp")
}
    n <- nrow(x)
    T <- ncol(x)
    if (length(SSIC.slack)==1) SSIC.slack = rep(SSIC.slack,2)
    if (is.null(max.q)) 
        max.q <- max(round(20, sqrt(min(n, T))))
    if (is.null(dw)) 
        dw <- round(min(log(T)^2, T^(6/7)/4))

    if(is.null(rule)) rule <- round(log(T, 2)/2)
    if(idio.norm.thr==0)SBS=0


    para<-list(M=M,r=r, SN.op=SN.op, is.NoF.fixed=r, max.q=max.q, sig.lev=sig.lev, dw=dw, rule=rule, normalisation=normalisation, refact.op=refact.op,
        common.norm=common.norm, idio.norm=idio.norm,idio.norm.thr=idio.norm.thr,norm.op=norm.op, SSIC.slack=SSIC.slack, 
        idio.diag =idio.diag, WBS=WBS,SBS=SBS,input.op=input.op, 
        do.parallel=do.parallel, no.proc = no.proc, setseed=setseed)

    gfm0 <- get.fac.mod(x, bn.op = bn.op, max.q = max.q, normalisation=normalisation)
    x <- gfm0$norm.x
    if (is.null(r)) 
        q.hat <- gfm0$q.hat
    else q.hat <- r # given NoF 

    if (do.parallel) {
        cl.cores <- parallel::detectCores()
        cl <- parallel::makeCluster(no.proc)
        doParallel::registerDoParallel(cl)
    }

    if(q.hat){#NoF>0
        cs <- common.seg(gfm0, q = q.hat, dw = dw, rule = rule, SN.op=0, sig.lev=sig.lev,
            common.norm=common.norm, common.norm.thr=0, norm.op=norm.op, SSIC.slack=SSIC.slack,
            WBS=WBS, M=M, input.op=input.op,
            do.parallel = do.parallel, setseed=setseed)
        common.est.cps <- sort(cs$est.cps)
        hat.vep<-get.residual(gfm0=gfm0, r=q.hat, common.est.cps=common.est.cps, refact=refact.op,bn.op = bn.op)
   }else{#NoF=0
        cs <- NULL
        common.est.cps <- as.numeric()
        hat.vep <- x
   } 
    is <- idio.seg(hat.vep, dw = dw, rule = rule,  SN.op=SN.op, sig.lev=sig.lev,
        idio.norm=idio.norm, idio.norm.thr=idio.norm.thr, norm.op=norm.op, SSIC.slack=SSIC.slack,
         WBS=WBS, SBS=SBS, M=M, idio.diag=idio.diag, input.op=input.op,
    do.parallel = do.parallel, setseed=setseed)
    para$idio.norm.thr <-is$idio.norm.thr

     idio.est.cps <- is$est.cps
    if (do.parallel) parallel::stopCluster(cl)
return(list(gfm = gfm0, para=para, hat.vep=hat.vep, common.seg.res=cs, common.est.cps = common.est.cps, 
        idio.seg.res = is, idio.est.cps = idio.est.cps))
}