rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm)
#setwd("/home/hutchings/Documents/PhD-year-1/PhD-year-1-code/Mini-project-1")
# Generate Data
N = 2000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:8, 2] = 1L
Lambda[9:10, 3] = 1L #doesnt work if small??
set.seed(1234)
Eta = matrix(NA, N, K)
Yfull = matrix(NA, N, D)
Sigma = diag(1L,D)
tmp = diag(1L,K)
tmp2 = rep(0L, K)
for (n in 1:N){
  Eta[n,] = mvrnorm(1, tmp2, tmp)
  Yfull[n,] = mvrnorm(1, Lambda%*%Eta[n,], Sigma)
}
#Yfull = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + diag(1L,D))

# Convert dimensions to discrete
discreteDims = c(1,7,9)
groups = list(c(0.1,0.1,0.5,0.3), c(0.1,0.4,0.2,0.2,0.1), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5))# ith element is probability of being in group i 
MakeDiscrete <- function(Y, discreteDims, groups, doplot=TRUE){
  discreteY = sapply(1:length(discreteDims), function(i) stepfun(cumsum(groups[[i]]), 1:(length(groups[[i]])+1L))( pnorm(scale(Y[, discreteDims[i]])) ) )
  Y[,discreteDims] = discreteY
  Y
}
Yfull = MakeDiscrete(Yfull, discreteDims, groups)
prob = 0.1 # probability of data begin missing


## set msk=1,2,3 dependign on how much missing data required
msk=2
if (msk==1){
  missingMask = matrix(0, N,D)
}else if (msk==2){
  prob=0.1
  missingMask = matrix(rbinom(N*D,1,prob=prob), nrow=N)
}else if (msk==3){
  missingMask = matrix(rbinom(N*D,1,prob=prob), nrow=N)
  missingMask[1:floor(N*0.75),c(1,3,5,6)] = 1
}
prob = mean(missingMask)
pdf(file= paste0("figures/missing-", msk, ".pdf"))
heatmap(missingMask, scale = "none", Rowv = NA, Colv = NA, col = c('white', 'black'), labRow = '', xlab = 'Dimension', ylab='Individual', main = paste('Missing Data Proportion: ', prob))
dev.off()
Y=Yfull
Y[as.logical(missingMask)] = NA 



CAVI <- function(Y, discreteDims=c(), maxiterations=1000L, tol = 0.001, seed=NULL){
  if (!is.null(seed)){set.seed(seed)}
  parameters = c()
  parameters$N = nrow(Y)
  parameters$D = ncol(Y)
  parameters$K = K
  parameters$a0 = parameters$b0 = 1e-3
  parameters$NdiscreteDims = length(discreteDims)
  
  boundaries = function(x){
  #assumes that classes are 1,...,n and all non empty
    toFind = 1L
    counts = c(0)
    count = sum(x == toFind, na.rm=TRUE)
    while (count) {
      counts = append(counts, count)
      toFind = toFind+1L
      count = sum(x == toFind, na.rm=TRUE)
    }
    counts = cumsum(counts)
    counts = qnorm(counts/counts[length(counts)])
    
    m = matrix(NA, nrow=length(counts)-1L, ncol=2)
    for (i in 1:nrow(m)){
      m[i,] = counts[i:(i+1L)]
    }
    m
    rbind(m, c(-Inf, Inf)) #adds a final boundary for unbounded truncated normal
  }
  
  
  parameters$Zboundaries =  apply(Y[,discreteDims], 2L, boundaries, simplify=FALSE)
  parameters$M.Z = matrix(rnorm(parameters$N*parameters$NdiscreteDims), nrow=parameters$N, ncol=parameters$NdiscreteDims) #not the actual mean! only the mean in the truncated normal
  #parameters$Zcov = diag(1L, nrow=parameters$N) #for each discrete dimension
  
  parameters$missingMask = is.na(Y)
  
  #calculate avg value for of each column, for initialisation
  Ymu = colMeans(Y, na.rm=TRUE)
  Ymu[discreteDims] = round(Ymu[discreteDims])
  # creates another category for missing data in the discrete variables
  for (i in discreteDims){
    Y[parameters$missingMask[,i], i] =  max(Y[,i], na.rm=TRUE)+ 1L
  }
  
  
  parameters$pseudoY = Y
  #initialise missing values of pseudoY with their mean
  for (d in 1:parameters$D){parameters$pseudoY[parameters$missingMask[,d],d] = Ymu[d]}
  
  
  # making all continuous
  parameters$pseudoY[, discreteDims] = sapply(1:parameters$NdiscreteDims, 
                                      function(i) rtruncnorm(1, a=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,1],
                                      b=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,2], 
                                      mean=parameters$M.Z[,i]))
  
  
  ## Initialise Lambda parameters
  # Initialise mean of Lambda by PCA 
  parameters$M.Lambda = unname(factanal(parameters$pseudoY, parameters$K)$loadings[1:parameters$D,])
  #Initialised by Wishart distribution, last dimension refers to the vector the covariance matrix corresponds to 
  parameters$S.Lambda =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K)) 
  
  ## Initialise Eta parameters by std normal and Wishart
  parameters$M.Eta = matrix(rnorm(parameters$N*parameters$K), nrow=parameters$N, ncol=parameters$K)
  parameters$S.Eta =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K))[,,1] #the same for each observation
  
  ## Initialise b parameters (a does not change)
  parameters$a = parameters$a0 + 0.5
  parameters$b = matrix(parameters$a, nrow=parameters$D, ncol=parameters$K) 
  
  
  
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
  
  #check
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
    
    #update the continuous missing ones:
    for (d in 1:parameters$D){
      if (!any(d==discreteDims)){
        tmp = (parameters$M.Eta%*%parameters$M.Lambda[d,])[,1]
        parameters$pseudoY[parameters$missingMask[,d],d] = tmp[parameters$missingMask[,d]]
      }
    }
    
    parameters$b = parameters$b0 + 0.5*E.Lambda2(parameters)
    
    cond = max(abs(oldLambda - parameters$M.Lambda))<tol
    if (cond){
      print(paste('Iterations', iter))
        break }
    oldLambda = parameters$M.Lambda
    #print(parameters$M.Lambda)
  }
  
  tmp = varimax(parameters$M.Lambda)
  R = tmp$rotmat
  Eta = parameters$M.Eta%*%R
  varimaxLambda = tmp$loadings[1:parameters$D,]
  
  for (i in 1:parameters$NdiscreteDims){
    d = discreteDims[i]
    #convert boundaries to vector
    boundaries = parameters$Zboundaries[[i]][,1]
    boundaries[length(boundaries)] = Inf
    print(boundaries)
    parameters$pseudoY[,d] = findInterval(parameters$pseudoY[,d] - (parameters$M.Eta%*%parameters$M.Lambda[d,]), boundaries)
  }
  print(parameters$M.Eta%*%parameters$M.Lambda[d,])
  return( list(Lambda=varimaxLambda, Eta=Eta, Yhat=parameters$pseudoY) )}


# Perform CAVI on data repeats times (with random initialisations)
repeats=1
results = vector(mode = "list", length = repeats)
for (i in 1:repeats){
  results[[i]] = CAVI(Y, discreteDims = discreteDims, seed = 1908+i)
  print(i/repeats)
}


print(Lambda)
R = matrix(0,3,3)
R[1,2]=R[2,1]=R[3,3]=1
predLambda = results[[i]]$Lambda%*%R
predEta = results[[i]]$Eta%*%R
round(predLambda,2)
# pdf(file= paste0("figures/heatmap1-", msk, ".pdf"))
# heatmap(Lambda, scale='none', Rowv = NA, Colv =NA, main = 'True Lambda')
# dev.off()
pdf(file= paste0("figures/heatmap2-", msk, ".pdf"))
heatmap(predLambda, scale = 'none', Rowv = NA, Colv =NA, main='VI Lambda')
dev.off()

#scaling Y then looking at Its errors
if (any((msk==c(2,3)))){
  Yhat = results[[i]]$Yhat
  tmp= as.logical(missingMask)
  Yhatscaled =Yfullscaled = Yhat*0
  for (d in 1:D){
    mu = mean(Yfull[, d])
    sdY = sd(Yfull[, d])
    Yhatscaled[,d] = (Yhat[,d]-mu)/sdY
    Yfullscaled[,d] = (Yfull[,d]-mu)/sdY
  }
  RMSerror = sqrt(mean((Yhatscaled[tmp] - Yfullscaled[tmp])**2))
  print(RMSerror)
  Errors = Yhatscaled[tmp] - Yfullscaled[tmp]
  pdf(file= paste0("figures/errorshist-", msk, ".pdf"))
  hist(Errors, main=paste('Histogram of errors, RMSE=', round(RMSerror,3)))
  #plot(sort(abs(Errors)), main=paste('Sorted error of Imputed values RMSE=', round(RMSerror,3)), lty=1)
  dev.off()
}

RMSerror = sqrt(mean((predEta- Eta)**2))
print(RMSerror)
Errors = as.vector(predEta - Eta)
pdf(file= paste0("figures/errorsetahist-", msk, ".pdf"))
hist(Errors, main=paste('Histogram of Eta errors, RMSE=', round(RMSerror,3)))
#plot(sort(abs(Errors)), main=paste('Sorted error of Imputed values RMSE=', round(RMSerror,3)), lty=1)
dev.off()


