rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm)

## Generate Data
N = 10000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:3, 1] = 1L
Lambda[4:7, 2] = 1L
Lambda[8:10, 3] = 1L #TODO doesnt work if this is small (eg 1)
Sigma = diag(1L,D)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + Sigma)

# Convert dimensions to discrete
discreteDims = c(1:3)
groups = list(c(0.1,0.1,0.5,0.3), c(0.1,0.4,0.2,0.2,0.1), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5))# ith element is probability of being in group i 
MakeDiscrete <- function(Y, discreteDims, groups, doplot=TRUE){
  discreteY = sapply(1:length(discreteDims), function(i) stepfun(cumsum(groups[[i]]), 1:(length(groups[[i]])+1L))( pnorm(scale(Y[, discreteDims[i]])) ) )
  Y[,discreteDims] = discreteY
  Y
}

i = 1

boundaries = c(0, cumsum(groups[[1]]))
boundaries = qnorm(boundaries)
boundaries[1] = -3
boundaries[length(boundaries)] = 3
pseudoBoundaries = boundaries

n = length(boundaries) -1L
temp = brewer.pal(n = n, name = "Dark2")
plot(pseudoBoundaries, dnorm(pseudoBoundaries), ylim=c(0, 0.5), xlab = 'z_n', ylab='')
for (j in 1:n){
  xx = linspace(pseudoBoundaries[j], pseudoBoundaries[j+1])
  lines(xx, dnorm(xx), col=temp[j])
}
