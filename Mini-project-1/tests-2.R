#rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, tmvtnorm, truncnorm)


E.Z <- function(parameters){
  sapply(1:parameters$NdiscreteDims, function(i) etruncnorm(a=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,1], b=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,2], mean=parameters$M.Z[,i]))
}
foo = E.Z(parameters)

i = 1

mu=parameters$M.Z[,i]
LowerUpper = parameters$Zboundaries[[i]][Y[,discreteDims[i]], ]
lower = temp[,1]
upper = temp[,2]
#mtmvnorm(mean = mu, lower=lower, upper=upper, doComputeVariance = FALSE)

temp = 

print(temp)
print(foo[[i]])
temp - foo[[i]]

# update Z (pseudoY with discrete dimensions)
parameters$pseudoY[,discreteDims] = tcrossprod(parameters$M.Eta, parameters$M.Lambda[discreteDims,])
i = 1
foo = parameters$M.Eta %*% parameters$M.Lambda[discreteDims[i],]


