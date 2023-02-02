rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm)

# Generate Data
N = 2000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:8, 2] = 1L
Lambda[9:10, 3] = 1L #doesnt work if small??
set.seed(1234)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + diag(1L,D))

# Convert dimensions to discrete
discreteDims = c(1:10)
groups = list(c(0.1,0.1,0.5,0.3), c(0.1,0.4,0.2,0.2,0.1), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5))# ith element is probability of being in group i 
MakeDiscrete <- function(Y, discreteDims, groups, doplot=TRUE){
  discreteY = sapply(1:length(discreteDims), function(i) stepfun(cumsum(groups[[i]]), 1:(length(groups[[i]])+1L))( pnorm(scale(Y[, discreteDims[i]])) ) )
  Y[,discreteDims] = discreteY
  Y
}
Y = MakeDiscrete(Y, discreteDims, groups)


CAVI <- function(Y, discreteDims=c(), maxiterations=1000L, tol = 0.1, seed=NULL){
  if (!is.null(seed)){set.seed(seed)}
  parameters = c()
  parameters$N = nrow(Y)
  parameters$D = ncol(Y)
  parameters$K = K
  parameters$a0 = parameters$b0 = 1e-3
  
  parameters$NdiscreteDims = length(discreteDims)
  #parameters$empericalCDFs = apply(Y[,discreteDims], 2L, ecdf)
  
  boundaries = function(x){
    toFind = 1L
    counts = c(0)
    count = sum(x == toFind)
    while (count) {
      counts = append(counts, count)
      toFind = toFind+1L
      count = sum(x == toFind)
    }
    counts = cumsum(counts)
    counts = qnorm(counts/counts[length(counts)])
    
    m = matrix(NA, nrow=length(counts)-1L, ncol=2)
    for (i in 1:nrow(m)){
      m[i,] = counts[i:(i+1L)]
    }
    m
  }
  
  Y[,discreteDims] = Y[,discreteDims] - (apply(Y[,discreteDims],2,min)-1L) #ensure starting 1L TODO: groups must be non empty and be 1,2,3,4....
  
  parameters$Zboundaries =  apply(Y[,discreteDims], 2L, boundaries, simplify=FALSE)
  parameters$M.Z = matrix(rnorm(parameters$N*parameters$NdiscreteDims), nrow=parameters$N, ncol=parameters$NdiscreteDims) #not the actual mean! only the mean in the truncated normal
  #parameters$Zcov = diag(1L, nrow=parameters$N) #for each discrete dimension
  
  parameters$pseudoY = Y
  parameters$pseudoY[, discreteDims] = sapply(1:parameters$NdiscreteDims, 
                                      function(i) rtruncnorm(1, a=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,1],
                                      b=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,2], 
                                      mean=parameters$M.Z[,i]))
  
  
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
  parameters$a = parameters$a0 + 0.5
  parameters$b = matrix(rgamma(parameters$D*parameters$K, shape=2, rate= 1/parameters$a), nrow=parameters$D, ncol=parameters$K) # has mean parameters$a
  
  
  # Functions for calculating expectations
  E.Lambda <- function(parameters){
    parameters$M.Lambda
  }
  E.LambdaTLambda <- function(parameters) {
    crossprod(parameters$M.Lambda) + rowSums(parameters$S.Lambda,dims=2L)
  }
  E.Lambda2 <-function(parameters){
    (parameters$M.Lambda**2) + t(apply(parameters$S.Lambda, 3L, diag))
  }
  E.Eta<- function(parameters){
    parameters$M.Eta
  }
  E.EtaTEta <- function(parameters) {
    crossprod(parameters$M.Eta) + parameters$N*parameters$S.Eta
  }
  E.Tau <- function(parameters) {
    parameters$a / parameters$b
  }
  
  E.Z <- function(parameters){
    sapply(1:parameters$NdiscreteDims, 
           function(i) etruncnorm(a=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,1],
                                  b=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,2], 
                                  mean=parameters$M.Z[,i]))
  }
 
  
  
  oldLambda = parameters$M.Lambda*Inf
  for (iter in 1:maxiterations){
    parameters$S.Eta = solve(E.LambdaTLambda(parameters) + diag(1,parameters$K,parameters$K)) # The same for each observation
    parameters$M.Eta = tcrossprod(parameters$pseudoY, tcrossprod(parameters$S.Eta, E.Lambda(parameters)))
    
    
    ETau = E.Tau(parameters)
    EEtaTEta = E.EtaTEta(parameters)
    for (i in 1:parameters$D){
      parameters$S.Lambda[,,i] = solve(EEtaTEta + diag(ETau[i,]))
      parameters$M.Lambda[i,] = parameters$S.Lambda[,,i]%*%(crossprod(parameters$M.Eta,parameters$pseudoY[,i]))
      if (i%in%discreteDims){
      }
    }
    
    # update Z (pseudoY with discrete dimensions)
    parameters$M.Z = tcrossprod(parameters$M.Eta, parameters$M.Lambda[discreteDims,])
    parameters$pseudoY[,discreteDims] = E.Z(parameters)
    
    parameters$b = parameters$b0 + 0.5 *E.Lambda2(parameters)
    
    cond = max(abs(oldLambda - parameters$M.Lambda))<tol
    if (cond){
      print(paste('Iterations', iter))
      break }
    oldLambda = parameters$M.Lambda
  }
  
  varimaxLambda = varimax(parameters$M.Lambda)$loadings[1:parameters$D,]
  return( list(Lambda=varimaxLambda) )}


# Perform CAVI on data repeats times (with random initialisations)
repeats=1
results = vector(mode = "list", length = repeats)
for (i in 1:repeats){
  results[[i]] = CAVI(Y, discreteDims = discreteDims, tol=0.1, seed = 1908+i)
  print(i/repeats)
}


VILambda = results[[i]]$Lambda
R = diag(1,K)
R[1,1]=-1
VILambda=VILambda%*%R
# pdf(file= "figures/heatmap1disconly.pdf")
# heatmap(Lambda, scale='none', Rowv = NA, Colv =NA, main = 'True Lambda discrete only')
# dev.off()
pdf(file= "figures/heatmap2disconly.pdf")
heatmap(VILambda, scale = 'none', Rowv = NA, Colv =NA, main='VI Lambda discrete only')
dev.off()