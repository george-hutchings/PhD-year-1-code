rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS)
a0 = b0 = 1e-3

N = 100
D = 10
K = 3

# Generate Data
Lambda = matrix(0, D, K)#DxK
Lambda[1:4, 1] = 1
Lambda[5:9, 2] = 2
Lambda[10, 3] = -1

#Sigma = diag(runif(D))
Sigma = diag(,D)
Y = mvrnorm(N, rep(0, D) , tcrossprod(Lambda) + Sigma)

parameters = c()
parameters$M.Lambda = matrix(1, nrow=D, ncol=K)
parameters$S.Lambda =  matrix(1, nrow =D, ncol=K)

parameters$M.Eta = matrix(0, nrow=N, ncol=K)
parameters$S.Eta =  matrix(1, nrow =K, ncol=K)

#parameters$M.Tau = matrix(1, nrow=D, ncol=K)
parameters$b = b0
parameters$a = a0

E.LambdaTLambda <- function(parameters) {
  crossprod(parameters$M.Lambda) + diag(colSums(parameters$S.Lambda))
}

E.Lambda <- function(parameters) {
  parameters$M.Lambda
}

E.Tau <- function(parameters) {
  parameters$a / parameters$b
}

Sum.n.E.Eta2 <- function(parameters) {
  N * diag(parameters$S.Lambda) + apply(parameters$M.Lambda, 2, function(x)
    sum(x ** 2)) #v
}

E.Eta <- function(parameters) {
  parameters$M.Eta
}

parameters_old = parameters

iterations = 10000
for (iter in 1:iterations){
parameters_old = parameters
parameters$a = a0 + 1

parameters$S.Eta = solve(E.LambdaTLambda(parameters) + diag(1,K,K)) #TODO dims to make it work, nb the same for each observation
parameters$M.Eta = tcrossprod(Y, tcrossprod(parameters$S.Eta, E.Lambda(parameters))) #TODO check correct mult with Y
parameters$S.Lambda = 1 / (E.Tau(parameters) + matrix(rep(Sum.n.E.Eta2(parameters), D), D, byrow=TRUE))


#paramters$M.Lambda
for (i in 1:D) {
  Yi = Y[, i]
  Mi.Lambda = parameters$M.Lambda[i,]
  for (j in 1:K) {
    Mi.0j = Mi.Lambda
    Mi.0j[j] = 0
    YiMj =sum(Y[, i] * parameters$M.Eta[,j])
    parameters$M.Lambda[i, j] = (YiMj - crossprod(Mi.0j, crossprod(parameters$M.Eta,
                                                    parameters$M.Eta[, j]))) * parameters$S.Lambda[i, j] #TODO this formula is wrong specifically the YiMi bit (M.Eta only goes up to K not D...)
  }
}
parameters$b = b0 + 0.5 * (parameters$S.Lambda + parameters$M.Lambda ** 2)

if (iter%%10000){print(iter/iterations)}

}