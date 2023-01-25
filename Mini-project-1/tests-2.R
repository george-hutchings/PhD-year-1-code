#rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm)


####FAVB with ARD continous variable
E_step_omega<-function(Y,B_k,B_sigma,sigma_k,A_0){
  p = ncol(Y)
  tmp = 0
  for (i in 1:p) {
    tmp = tmp + sigma_k[i]*(B_k[i,]%*%t(B_k[i,])+B_sigma[,,i])
  }
  Matrix     <-tmp+A_0
  inv_Matrix <-solve(Matrix)
  omega      <- t((inv_Matrix%*%t(B_k)*sigma_k)%*%t(Y))
  list(omega=omega,Matrix=Matrix,inv_Matrix=inv_Matrix)
}

E_step_lambda<-function(Y,M,omega,Ealpha,sigma_k){
  n=nrow(Y)
  p=ncol(Y)
  k=ncol(omega)
  tmp=0
  for (i in 1:n) {
    tmp=tmp+(M+omega[i,]%*%t(omega[i,]))
  }
  B_sigma   = array(0,c(k,k,p))
  B_new     = matrix(0,p,k)
  for (i in 1:p) {
    B_sigma[,,i] <- solve(sigma_k[i]*tmp+Ealpha)
    B_new[i,]    <- ((t(Y[,i]))%*%omega)%*%B_sigma[,,i]*sigma_k[i]
  }
  return(list(B_new=B_new,B_sigma=B_sigma))
}





n=1000L
k      = 4
p      = 8# = d
#B_kT   = matrix(c(0.9078809,0,0, 0, 0.6828301, 0.7861421, 0, 0, 0, 0,0.8, 0.0, 0.0, 0.0 ,0.0,0.0, 0,1,0,0,0,0,0,0,0.7,0, 0, 0.7, 0, 0, 0.651872, 0.85 ),p,k)
B_kT   = matrix(c(0.9078809,0,0, 0, 0.6828301, 0.7861421, 0, 0, 0, 0,0.8, 0.5, 0.0, 0.0 ,0.0,0.0, .5,.7,0,0,0,0,0,0,0.7,0, 0, 0.7, 0, 0, 0.651872, 0.85 ),p,k)
omegaT = matrix(rnorm(n*k),n,k) #true etas
muT    = omegaT%*%t(B_kT)
Y      = muT+matrix(rnorm(n*p,0,.5),n,p)


# initialise eta, Lambda
omega   = matrix(rnorm(n*k),n,k)
B_k     = matrix(rnorm(p*k),p,k)
conv=TRUE #whether or not converged 
iter=0
Ealpha = diag(k)
B_k     = matrix(rnorm(p*k),p,k)
B_sigma = array(0,c(k,k,p))
sigma_k = rep(1,p)
A_0     = diag(k) #identity matrix
b0 =.001
alpha0=.001

for(i in 1:p){B_sigma[,,i]=diag(k)} 
while (conv) {
  iter=iter+1
  
  E_step1 <- E_step_omega(Y,B_k,B_sigma,sigma_k,A_0)
  omega 	<-E_step1$omega
  M			  <-E_step1$inv_Matrix
  #M-Step
  #B_new     <- (t(Z)%*%omega)%*%solve(n*M+t(omega)%*%omega)
  E_step2 <-E_step_lambda(Y,M,omega,Ealpha,sigma_k)
  B_new   <- E_step2$B_new
  B_sigma <- E_step2$B_sigma
  
  Alpha     <- .5*(diag(t(B_new)%*%B_k+apply(B_sigma,1:2,sum)))+b0
  Ealpha    <-diag((alpha0+p/2)/Alpha)
  
  mm        <-omega%*%t(B_new)
  for (i in 1:p) {
    tmp = 0
    E_BB = B_sigma[,,i]+B_new[i,]%*%t(B_new[i,])
    for (j in 1:n) {
      E_FF = omega[j,]%*%t(omega[j,])+M
      tmp = tmp + Y[j,i]^2-2*Y[j,i]*mm[j,i]+sum(diag(E_BB%*%E_FF))
    }
    sigma_k[i] = (alpha0+n/2)/(.5*tmp+b0)
  }
  
  #sigma_k   <- diag(t(Z)%*%Z+E_stepZ$Vz-B_new%*%t(omega)%*%Z)/n
  
  
  thr = max((B_new-B_k))
  if(thr<.01){conv=FALSE}
  print(iter)
  print(B_new)
  print(sigma_k)
  print(thr)
  B_k = B_new
}

