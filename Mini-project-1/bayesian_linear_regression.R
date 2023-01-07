rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS)

set.seed(10)
lambda = 0.5
alpha0 = 1e-1
beta0 = 1e-1
I = 100000


set.seed(12345)       # For reproducibility
D     <- 1            # For simplicity we have only one covariate
coefs <- c(-0.5, 2)   # Generate y with these coefficients
# Precision parameter

X <- cbind(1, replicate(D, rnorm(I)))             # Generate X data
y <- X %*% coefs  + rnorm(I, sd = sqrt(1/lambda)) # Generate y data



## simulate data 


E.wTw = 0
E.wTw.old = 0
E.Tau.old = 0

for (i in 1:100){
Alpha = alpha0 + D/2
Beta = beta0 + E.wTw/2
E.Tau = Alpha/Beta

Sinv =  lambda*t(X)%*%X
diag(Sinv)= diag(Sinv) + E.Tau


m = lambda*solve(Sinv, t(X)%*%y)

mTm = norm(m, type="2")**2
E.wTw = sum(1/eigen(Sinv, symmetric =TRUE, only.values=TRUE)$values) + mTm
E.wTw.old = E.wTw
E.Tau.old = E.Tau
}
print(abs(E.wTw-E.wTw.old))
print(abs(E.Tau - E.Tau.old))

Xstar = cbind(1, c(-1500:1500)/1000)
ystar = Xstar%*%m
plot(X[,2],y)
lines(Xstar[,2],ystar)
S = solve(Sinv)
var = diag((Xstar%*%S%*%t(Xstar)))+1/lambda
lines(Xstar[,2], ystar + 1.96*sqrt(var), col='red')
lines(Xstar[,2], ystar - 1.96*sqrt(var), col='red')
print(mean(var))

tau = c(0:1000)/100
plot(tau, dgamma(tau, shape=Alpha, rate=Beta),'l')

print(mean(diag((Xstar%*%S%*%t(Xstar)))))




