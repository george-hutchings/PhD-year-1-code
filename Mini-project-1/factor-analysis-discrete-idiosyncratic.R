rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer, truncnorm)

## Generate Data
N = 8000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:8, 2] = 1L
Lambda[9:10, 3] = 1L #TODO doesnt work if this is small (eg 1)
Sigma = diag(0.5,D)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + Sigma)

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
  parameters$alpha0 = parameters$beta0 = 1e-3
  
  parameters$NdiscreteDims = length(discreteDims)
  
  parameters$empericalCDFs = apply(Y[,discreteDims], 2L, ecdf)
  
  # # seems no point to do this?
  # foo = sapply(seq_len(parameters$NdiscreteDims), function(i) qnorm(parameters$empericalCDFs[[i]]((Y[,discreteDims[i]]))))
  # i=1
  # plot(Y[,discreteDims[i]])
  # plot(foo[,i]) #inf not plotted
  # plot(scale(Y[,discreteDims[i]])[,1])
  
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
  
  parameters$pseudoY = Y
  parameters$pseudoY[, discreteDims] = sapply(1:parameters$NdiscreteDims, 
                                      function(i) rtruncnorm(1, a=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,1],
                                      b=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,2], 
                                      mean=parameters$M.Z[,i]))
  
  ## Initialise Lambda parameters
  # Initialise mean of Lambda by PCA 
  tmp = factanal(parameters$pseudoY, parameters$K, scores = 'reg')
  parameters$M.Lambda = unname(tmp$loadings[1:parameters$D,])
  tmp = tmp$scores
  #Initialised by Wishart distribution, last dimension refers to the vector the covariance matrix corresponds to 
  parameters$S.Lambda =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K)) 
  
  ## Initialise Eta parameters by std normal and Wishart
  parameters$M.Eta = matrix(rnorm(parameters$N*parameters$K), nrow=parameters$N, ncol=parameters$K)
  parameters$S.Eta =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K))[,,1] #the same for each observation
  
  ## Initialise b parameters (a does not change)
  parameters$b = matrix(rgamma(parameters$D*parameters$K, shape=2, rate= 1/parameters$b0), nrow=parameters$D, ncol=parameters$K) 
  parameters$a = parameters$a0 + 1/2
  
  ## Initialise beta parameters (alpha does not change)
  parameters$alpha = parameters$alpha0 + 0.5*parameters$N
  parameters$beta = rgamma(parameters$D, shape=parameters$alpha, rate= 1/diag(var(parameters$pseudoY)))
  parameters$beta= parameters$alpha/diag(cov(parameters$pseudoY -tcrossprod(tmp, parameters$M.Lambda)))
  parameters$beta = rep(parameters$alpha, parameters$D)
  
  
  # Functions for calculating expectations
  E.Lambda <- function(parameters){
    parameters$M.Lambda
  }
  E.LambdaTHLambda <-function(parameters){
    m = parameters$S.Lambda
    tempE.H = E.H(parameters)
    for (i in 1:parameters$D){
      m[,,i] = parameters$S.Lambda[,,i]*tempE.H[i]
    }
    crossprod(parameters$M.Lambda*tempE.H, parameters$M.Lambda) + rowSums(m, dims=2L)
  }
  
  E.H <- function(parameters) {
    parameters$alpha / parameters$beta
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
                                  mean=parameters$M.Z[,i], sd = 1))
  }
  
  S.Z <- function(parameters){
    sapply(1:parameters$NdiscreteDims, 
           function(i) vtruncnorm(a=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,1],
                                  b=parameters$Zboundaries[[i]][Y[,discreteDims[i]], ][,2], 
                                  mean=parameters$M.Z[,i], sd = 1))
  }
 
  
  ## Function for calculating ELBO TODO: check it
  ELBO <- function(parameters){
    pt1 = numeric(4)
    pt1[1] = -0.5*(parameters$N*sum(diag(E.LambdaTHLambda(parameters)%*%parameters$S.Eta)) + 
                     -2*sum((Y%*%parameters$M.Lambda)*parameters$M.Eta) + 
                     sum((parameters$M.Eta%*%E.LambdaTHLambda(parameters))*parameters$M.Eta))
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
    EH = E.H(parameters)
    
    parameters$S.Eta = solve(E.LambdaTHLambda(parameters) + diag(1,parameters$K)) # The same for each observation
    parameters$M.Eta = tcrossprod(parameters$pseudoY, tcrossprod(parameters$S.Eta, E.Lambda(parameters)*EH))
    
    
    ETau = E.Tau(parameters)
    EEtaTEta = E.EtaTEta(parameters)
    for (i in 1:parameters$D){
      parameters$S.Lambda[,,i] = solve((EH[i]*EEtaTEta) + diag(ETau[i,]))
      parameters$M.Lambda[i,] = parameters$S.Lambda[,,i]%*%(crossprod(parameters$M.Eta, parameters$pseudoY[,i]))*EH[i]
    }
    
    # update Z (pseudoY with discrete dimensions)
    parameters$M.Z = tcrossprod(parameters$M.Eta, parameters$M.Lambda[discreteDims,])
    parameters$pseudoY[,discreteDims]   = E.Z(parameters) #=  t(t(E.Z(parameters))/sqrt(EH))
    
    parameters$b = parameters$b0 + 0.5 *E.Lambda2(parameters)
    
    varZ = S.Z(parameters) #= t(t(S.Z(parameters))/EH)
    for (j in 1:parameters$D){
      tmpcond = (j==discreteDims)
      if (any(tmpcond)) {
        tmp = sum(varZ[, which(tmpcond)])
      } else{
       tmp = 0
      }
      parameters$beta[j] = tmp + sum(parameters$pseudoY[, j]**2)
      - 2*crossprod(parameters$pseudoY[,j],parameters$M.Eta)%*%parameters$M.Lambda[j,]
      + sum(diag(EEtaTEta%*%parameters$S.Lambda[,,j])) + crossprod(parameters$M.Lambda[j,], EEtaTEta%*%parameters$M.Lambda[j,])
    }
    
    parameters$beta = parameters$beta0 + 0.5*parameters$beta
    
    #ELBOvec[iter] = ELBO(parameters)
    if (iter>5){
      val = ELBOvec[(iter-5):iter]
      if (all(abs(val - mean(val))<tol)){ 
        ELBOvec = ELBOvec[1:iter]
        break }
    }
    #print(parameters$M.Lambda)
    #print(EH)
  }
  print(parameters$M.Lambda)
  varimaxLambda = varimax(parameters$M.Lambda)$loadings[1:parameters$D,]
  return( list(Lambda=varimaxLambda, ELBO=ELBOvec) )}


# Perform CAVI on data repeats times (with random initialisations)
repeats=5
results = vector(mode = "list", length = repeats)
ELBOlast = numeric(repeats)
for (i in 1:repeats){
  results[[i]] = CAVI(Y, maxiterations= 50, discreteDims = discreteDims, tol=0.1, seed = 1908+i)
  ELBOlast[i] = results[[i]]$ELBO[length(results[[i]]$ELBO)]
  print(i/repeats)
}

idx = which.max(ELBOlast)
plot(results[[idx]]$ELBO, type='o', col='red', xlab = 'Iteration', ylab= 'ELBO')
temp = c(1:repeats)
for (j in sample(temp[-idx], 4)){
  lines(results[[j]]$ELBO)
}

for (i in 1:repeats){print(round(results[[i]]$Lambda,1))}
