#rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm, microbenchmark)

# Generate Data
N = 2000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:8, 2] = 1L
Lambda[9:10, 3] = 1L #doesnt work if small??
Sigma = diag(1L, D)
Eta = matrix(rnorm(N*K), nrow=N)
Y = matrix(NA, nrow=N, ncol = D)
for (i in 1:N) {Y[i,] = rnorm(Lambda%*%Eta[i,], diag(Sigma))}

parameters = c()
parameters$N = nrow(Y)
parameters$D = ncol(Y)
parameters$K = K
parameters$a0 = parameters$b0 = 1e-3
parameters$alpha0 = parameters$beta0 = 1e-3


## Initialise Lambda parameters
# Initialise mean of Lambda by PCA 
pY = prcomp(Y, scale=TRUE, rank=parameters$K)
parameters$M.Lambda = pY$rotation%*%diag(pY$sdev[1:parameters$K]**2) 
#Initialised by Wishart distribution, last dimension refers to the vector the covariance matrix corresponds to 
parameters$S.Lambda =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K)) 

## Initialise Eta parameters by std normal and Wishart
parameters$M.Eta = matrix(rnorm(parameters$N*parameters$K), nrow=parameters$N, ncol=parameters$K)
parameters$S.Eta =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K))[,,1] #the same for each observation

## Initialise b parameters (a does not change)
parameters$b = rgamma(parameters$K, shape=2, rate= 1/parameters$b0)
parameters$a = parameters$a0 + (0.5*parameters$D)

## Initialise b parameters (a does not change)
parameters$beta = rgamma(parameters$D, shape=2, rate= 1/parameters$beta0)
parameters$alpha = parameters$alpha0 + (0.5*parameters$N)


E.LambdaTHLambda <- function(parameters){
  
}


  
  
  
  
  
  
  
  