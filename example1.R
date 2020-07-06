## This example create an approximate factor model with 5 factors.
## Breaks for common components and idiosyncratic components are  

###Data Generation
n <- 200; T <- 200
r <- 5 # factor number
scenario = 1 # 0 for no change point
eta1 = round(T/3) #cpt at T=67 in 200
eta2 = round(2*T/3) #cpt at T=133 in 200
beta = 0.2
beta_i = matrix(runif(n,min=-beta,max=beta), nrow=n)
rho_f = 0.4
rho_fi = rho_f-0.05*(0:(r-1))


u = matrix(rnorm(r*T), nrow=r)

phi = 1

theta = phi # * r /(1-rho_f^2) 

# idiosyncratic components
rho_epsilon <- -0.5
Sigma0 <- matrix(runif(n,min=0.5,max=1.5), nrow=n) ##variance
Sigma <- rho_epsilon^ abs(matrix(rep(1:n,times=n),n,n) - matrix(rep(1:n,each=n),n,n) )
Sigma <- Sigma0%*%t(Sigma0)*Sigma
epsilon <- t(MASS::mvrnorm(n=T, rep(0, n), Sigma))
varpho = 0.1 ##0.1,0.5,1
cardS = floor(varpho*n/2)
ranS1 = sample(1:n, cardS, replace = F)
ranS2 = sample(1:(n-cardS), cardS, replace = F)
ranS2 = (1:n)[-ranS1][ranS2]
epsilon0 = epsilon[ranS1,(eta2+1):T]
epsilon[ranS1,(eta2+1):T] = epsilon[ranS2,(eta2+1):T]
epsilon[ranS2,(eta2+1):T] = epsilon0


f <- matrix(0, nrow=r, ncol=T) # factors
f[,1] = u[,1]
if (scenario==1){
	for (TT in 2:T){
		f[,TT] <- rho_fi*f[,TT-1]+u[,TT]
	}
}


Lam <- matrix(rnorm(n*r), nrow=n, ncol=r) # loadings
if (scenario==1){
	varpho = 1 # NOT CHANGE
	cardS = round(varpho*n)
	ranS = sample(1:n, cardS, replace = F)
  	deltasigma = sqrt(4)
	delta = matrix(rnorm(cardS*r) * deltasigma, nrow=cardS)

	chi <- epsilon*0 # common component
	chi[, 1:eta1] <- Lam%*%f[, 1:eta1] # eta # change-point
	Lam[ranS,] <- Lam[ranS,] +  delta 
	chi[, (eta1+1):T] <- Lam%*%f[, (eta1+1):T]
}
#epsilon <- matrix(rnorm(n*T), nrow=n)
x <- chi + sqrt(theta) * epsilon


###Detection
###BSCOV
res1 = BSCOV(x,do.parallel=0,WBS=0,SBS=0, M=400,no.proc=3,SN.op=3, norm.op=1,idio.norm.thr=0, idio.diag=T,bn.op=2)
###SBSCOV
res1 = BSCOV(x,do.parallel=0,WBS=0,SBS=1, M=400,no.proc=3,SN.op=3, norm.op=1,idio.norm.thr=0, idio.diag=T,bn.op=2)
###WBSCOV
res1 = BSCOV(x,do.parallel=0,WBS=1,SBS=0, M=400,no.proc=3,SN.op=3, norm.op=1,idio.norm.thr=0, idio.diag=T,bn.op=2)
###WSBSCOV
res1 = BSCOV(x,do.parallel=0,WBS=1,SBS=1, M=400,no.proc=3,SN.op=3, norm.op=1,idio.norm.thr=0, idio.diag=T,bn.op=2)




