rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer)

## Generate Data
N = 1000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:9, 2] = 2L
Lambda[10, 3] = 1L
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + diag(1L,D))


discreteDims = c(2, ncol(Y))
boundaries = list(c(-2, -1, 1, 2 ), c(-2, 0, 1 ))
# Convert dimension to discrete
MakeDiscrete <- function(Y, discreteDims, boundaries, doplot=TRUE){
discreteParams = vector(mode = "list", length = length(discreteDims))
for (i in 1:length(discreteDims)){
  discreteParams[[i]]$dimension = discreteDims[i]
  discreteParams[[i]]$boundary = boundaries[[i]]
  discreteParams[[i]]$Ngroups = length( discreteParams[[i]]$boundary)+1
  discreteParams[[i]]$scaledY = scale(Y[,discreteParams[[i]]$dimension])
  discreteParams[[i]]$discreteY = as.integer(findInterval(discreteParams[[i]]$scaledY, c(-Inf, discreteParams[[i]]$boundary, Inf)))
  if (doplot){
    temp = brewer.pal(n = discreteParams[[i]]$Ngroups, name = "Dark2")
    pseudoBoundaries = c(discreteParams[[i]]$boundary[1]-3,
                         discreteParams[[i]]$boundary, 
                         discreteParams[[i]]$boundary[length(discreteParams[[i]]$boundary)] +3)
    plot(pseudoBoundaries, dnorm(pseudoBoundaries), ylim=c(0, 0.5), xlab = 'x', ylab='')
    for (j in 1:discreteParams[[i]]$Ngroups){
      xx = linspace(pseudoBoundaries[j], pseudoBoundaries[j+1])
      lines(xx, dnorm(xx), col=temp[j])
    }
    plot(Y[,discreteParams[[i]]$dimension], discreteParams[[i]]$discreteY)
  }
  Y[, discreteParams[[i]]$dimensio] = discreteParams[[i]]$discreteY
}
return(Y)
}
Y = MakeDiscrete(Y, discreteDims, boundaries)


y = Y[,1] # Y[,discreteDims[1]]
x = ecdf(y)
plot(x)

set.seed = 1908
z = rnorm(N, 2, 5)

z=y
zhat = scale(z)[,1]
foo  = pnorm(zhat)
yhat = quantile(x, pnorm(zhat))
#plot(z, col = yhat)
plot(yhat)
plot(y)
plot(y, yhat)
abline(a =0, b=1)
              




