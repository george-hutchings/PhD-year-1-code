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
Lambda[1:2, 1] = 1L
Lambda[3:7, 2] = 1L
Lambda[8:10, 3] = 1L #TODO doesnt work if this is small (eg 1)
Sigma = diag(1L,D)
Yfull = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + Sigma)

# Convert dimensions to discrete
discreteDims = c(1:10)
groups = list(c(0.1,0.1,0.5,0.3), c(0.1,0.4,0.2,0.2,0.1), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5), c(0.5,0.5))# ith element is probability of being in group i 
MakeDiscrete <- function(Y, discreteDims, groups, doplot=TRUE){
  discreteY = sapply(1:length(discreteDims), function(i) stepfun(cumsum(groups[[i]]), 1:(length(groups[[i]])+1L))( pnorm(scale(Y[, discreteDims[i]])) ) )
  Y[,discreteDims] = discreteY
  Y
}
Yfull = MakeDiscrete(Yfull, discreteDims, groups)
prob = 0.05 # probability of data begin missing
missingMask = matrix(rbinom(N*D,1,prob=prob), nrow=N)
heatmap(missingMask, scale = "none", Rowv = NA, Colv = NA, col = c('white', 'black'), main = paste('p = ', prob))
Y=Yfull
Y[as.logical(missingMask)] = NA 



CAVI <- function(Y, discreteDims=c(), maxiterations=1000L, tol = 0.1, seed=NULL){
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
    rbind(m, c(-Inf, Inf))
  }
  
  
  parameters$Zboundaries =  apply(Y[,discreteDims], 2L, boundaries, simplify=FALSE)
  parameters$M.Z = matrix(rnorm(parameters$N*parameters$NdiscreteDims), nrow=parameters$N, ncol=parameters$NdiscreteDims) #not the actual mean! only the mean in the truncated normal
  #parameters$Zcov = diag(1L, nrow=parameters$N) #for each discrete dimension
  
  parameters$missingMask = is.na(Y)
  #calculate avg value for of each column, for initialisation
  Ymu = colMeans(Y, na.rm=TRUE)
  Ymu[discreteDims] = round(Ymu[discreteDims])
  # creates another category for missing data in the discrete variables
  for (i in 1:parameters$NdiscreteDims){
    d = discreteDims[i]
    Y[parameters$missingMask[,i], discreteDims[i]] =  max(Y[,discreteDims[i]], na.rm=TRUE)+ 1L
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
  parameters$b = matrix(rgamma(parameters$D*parameters$K, shape=2, rate= 1/parameters$a), nrow=parameters$D, ncol=parameters$K) 
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
 
  
  ## Function for calculating ELBO TODO: check it
  ELBO <- function(parameters){
    pt1 = numeric(4)
    pt1[1] = -0.5*(parameters$N*sum(diag(E.LambdaTLambda(parameters)%*%parameters$S.Eta)) + 
                     -2*sum((Y%*%parameters$M.Lambda)*parameters$M.Eta) + 
                     sum((parameters$M.Eta%*%E.LambdaTLambda(parameters))*parameters$M.Eta))
    pt1[2] = -0.5*(sum(diag(parameters$S.Eta)) + sum(parameters$M.Eta**2))
    pt1[3] = -0.5*( sum(log(parameters$b)) + 
                      parameters$a*sum(t(apply(parameters$S.Lambda, 3L, diag))/parameters$b) +
                      parameters$a*sum(parameters$M.Lambda*parameters$M.Lambda/parameters$b) )
    pt1[4] = (1-parameters$a0)*sum(log(parameters$b)) -parameters$b0*parameters$a*sum(1/parameters$b)
    
    pt2 = numeric(3)
    pt2[1] = N*0.5*determinant(parameters$S.Eta)$modulus
    pt2[2] = 0.5*sum(apply(parameters$S.Lambda, 3L, function(x) determinant(x)$modulus))
    pt2[3]= -sum(log(parameters$b))
    sum(pt1, pt2)
  }
  
  
  
  # Loop for performing CAVI, until ELBO - mean(ELBO) is < tol for 5 iterations
  ELBOvec = c(1:maxiterations) #numeric(maxiterations)
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
    
    parameters$b = parameters$b0 + 0.5*E.Lambda2(parameters)
    
    #ELBOvec[iter] = ELBO(parameters)
    if (iter>5){
      val = ELBOvec[(iter-5):iter]
      if (all(abs(val - mean(val))<tol)){ 
        ELBOvec = ELBOvec[1:iter]
        break }
    }
    #print(parameters$M.Lambda)
  }
  
  varimaxLambda = varimax(parameters$M.Lambda)$loadings[1:parameters$D,]
  return( list(Lambda=varimaxLambda, ELBO=ELBOvec) )}


# Perform CAVI on data repeats times (with random initialisations)
repeats=5
results = vector(mode = "list", length = repeats)
ELBOlast = numeric(repeats)
for (i in 1:repeats){
  results[[i]] = CAVI(Y,maxiterations = 200, discreteDims = discreteDims, tol=0.1, seed = 1908+i)
  ELBOlast[i] = results[[i]]$ELBO[length(results[[i]]$ELBO)]
  print(i/repeats)
}

idx = which.max(ELBOlast)
plot(results[[idx]]$ELBO, type='o', col='red', xlab = 'Iteration', ylab= 'ELBO')
temp = c(1:repeats)
for (j in sample(temp[-idx], 4)){
  lines(results[[j]]$ELBO)
}
print(ELBOlast[idx])
VILambda = results[[idx]]$Lambda
print(Lambda)
print(round(VILambda,3))

for (i in 1:repeats){print(round(results[[i]]$Lambda,2))}
