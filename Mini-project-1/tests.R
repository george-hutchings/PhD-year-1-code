#rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm, microbenchmark)

## Generate Data
N = 100L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:9, 2] = 2L
Lambda[10, 3] = 3L #TODO doesnt work if this is small (eg 1)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + diag(1L,D))


parameters$beta = numeric(parameters$D)
trASigma = sum(diag(EEtaTEta)%*%parameters$S.Eta)
for (j in 1:parameters$D){
  parameters$beta[j] =  sum(Y[,j]**2) - 
    2*(Y[,j]%*%parameters$M.Eta)%*%parameters$M.Lambda[j,]
+ trASigma + crossprod(parameters$M.Lambda[j,], EEtaTEta%*%parameters$M.Lambda[j,])
}
parameters$beta = parameters$beta0 + 0.5*parameters$beta
