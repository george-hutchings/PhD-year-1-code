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


# Generate Data
N = 2000L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:8, 2] = 1L
Lambda[9:10, 3] = 1L #doesnt work if small??
Sigma = diag(0.5, D)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + Sigma)



CAVIfewTau <- function(Y, maxiterations=1000L, tol = 0.1, seed=NULL){
  if (!is.null(seed)){set.seed(seed)}
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
  
  ## Initialise beta parameters (alpha does not change)
  parameters$alpha = parameters$a0 + 0.5*parameters$N
  parameters$beta = rgamma(parameters$D, shape=parameters$alpha, rate= 1/diag(var(Y)))
  #parameters$beta = rep(parameters$alpha *0.5, parameters$D)
  
  
  # Functions for calculating expectations
  E.Lambda <- function(parameters){
    parameters$M.Lambda
  }
  E.LambdaTLambda <- function(parameters) {
    crossprod(parameters$M.Lambda) + rowSums(parameters$S.Lambda,dims=2L)
  }
  
  #checked 
  E.H <- function(parameters) {
    parameters$alpha / parameters$beta
  }
  
  #checked
  E.LambdaTHLambda <-function(parameters){
    m = parameters$S.Lambda
    tempE.H = E.H(parameters)
    for (i in 1:parameters$D){
      m[,,i] = parameters$S.Lambda[,,i]*tempE.H[i]
    }
    crossprod(parameters$M.Lambda*tempE.H, parameters$M.Lambda) + rowSums(m, dims=2L)
  }
  
  ElambdaTlambda <-function(parameters){
    tmp = colSums(parameters$M.Lambda**2)
    for (i in 1:parameters$K)
      tmp[i] = tmp[i] + sum(parameters$S.Lambda[i,i,])
    tmp
  }
  #checked
  E.Eta<- function(parameters){
    parameters$M.Eta
  }
  #checked
  E.EtaTEta <- function(parameters) {
    crossprod(parameters$M.Eta) + parameters$N*parameters$S.Eta
  }
  E.Tau <- function(parameters) {
    parameters$a / parameters$b
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
    EH = E.H(parameters)
    
    parameters$S.Eta = solve(E.LambdaTHLambda(parameters) + diag(1,parameters$K)) # The same for each observation #checked
    parameters$M.Eta = tcrossprod(Y, tcrossprod(parameters$S.Eta, E.Lambda(parameters)*EH)) #checked
    
    
    
    ETau = E.Tau(parameters)
    EEtaTEta = E.EtaTEta(parameters)
    for (i in 1:parameters$D){
      parameters$S.Lambda[,,i] = solve((EH[i]*EEtaTEta) + diag(ETau))
      parameters$M.Lambda[i,] = parameters$S.Lambda[,,i]%*%(crossprod(parameters$M.Eta,Y[,i]))*EH[i]
    }
    
    parameters$b = parameters$b0 + 0.5 *ElambdaTlambda(parameters)
    
    for (j in 1:parameters$D){
      parameters$beta[j] =  sum(Y[,j]**2) 
      - 2*(crossprod(Y[,j],parameters$M.Eta)%*%parameters$M.Lambda[j,])
      + sum(diag(EEtaTEta%*%parameters$S.Lambda[,,j])) + crossprod(parameters$M.Lambda[j,], EEtaTEta%*%parameters$M.Lambda[j,])
    }
    parameters$beta = parameters$beta0 + 0.5*parameters$beta
    
    # ELBOvec[iter] = ELBO(parameters) #UNCOMMENT TO IMPLEMENT ELBO
    if (iter>5){
      val = ELBOvec[(iter-5):iter]
      if (all(abs(val - mean(val))<tol)){ 
        ELBOvec = ELBOvec[1:iter]
        break }
    }
  }
  
  #print(1/sqrt(E.H(parameters)))
  
  varimaxLambda = varimax(parameters$M.Lambda)$loadings[1:parameters$D,]
  return( list(Lambda=varimaxLambda, ELBO=ELBOvec) )}


CAVIlotsTau <- function(Y, maxiterations=1000L, tol = 0.1, seed=NULL){
  if (!is.null(seed)){set.seed(seed)}
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
  parameters$a = parameters$a0 + 0.5
  parameters$b = matrix(rgamma(parameters$D*parameters$K, shape=2, rate= 1/parameters$alpha), nrow=parameters$D, ncol=parameters$K) 
  
  
  ## Initialise beta parameters (alpha does not change)
  parameters$alpha = parameters$a0 + 0.5*parameters$N
  parameters$beta = rgamma(parameters$D, shape=parameters$alpha, rate= 1/diag(var(Y)))
  
  
  # Functions for calculating expectations
  E.Lambda <- function(parameters){
    parameters$M.Lambda
  }
  E.LambdaTLambda <- function(parameters) {
    crossprod(parameters$M.Lambda) + rowSums(parameters$S.Lambda,dims=2L)
  }
  
  #checked 
  E.A <- function(parameters) {
    parameters$alpha / parameters$beta
  }
  
  #checked
  E.LambdaTALambda <-function(parameters){
    m = parameters$S.Lambda
    tempE.A = E.A(parameters)
    for (i in 1:parameters$D){
      m[,,i] = parameters$S.Lambda[,,i]*tempE.A[i]
    }
    crossprod(parameters$M.Lambda*tempE.A, parameters$M.Lambda) + rowSums(m, dims=2L)
  }
  
  E.Lambda2 <-function(parameters){
    (parameters$M.Lambda**2) + t(apply(parameters$S.Lambda, 3L, diag))
  }
  #checked
  E.Eta<- function(parameters){
    parameters$M.Eta
  }
  #checked
  E.EtaTEta <- function(parameters) {
    crossprod(parameters$M.Eta) + parameters$N*parameters$S.Eta
  }
  E.Tau <- function(parameters) {
    parameters$a / parameters$b
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
    EA = E.A(parameters)
    
    parameters$S.Eta = solve(E.LambdaTALambda(parameters) + diag(1,parameters$K,parameters$K)) # The same for each observation #checked
    parameters$M.Eta = tcrossprod(Y, tcrossprod(parameters$S.Eta, E.Lambda(parameters)*EA)) #checked
    
    
    
    ETau = E.Tau(parameters)
    EEtaTEta = E.EtaTEta(parameters)
    for (i in 1:parameters$D){
      parameters$S.Lambda[,,i] = solve((EA[i]*EEtaTEta) + diag(ETau[i,]))
      parameters$M.Lambda[i,] = parameters$S.Lambda[,,i]%*%(crossprod(parameters$M.Eta,Y[,i]))*EA[i]
    }
    
    parameters$b = parameters$b0 + 0.5*E.Lambda2(parameters) #todo: this is producing negative betas
    stopifnot(all(E.Lambda2(parameters)>0))
    
    #trEEtaTEtaSigma = sum(diag(EEtaTEta%*%parameters$S.Lambda[,,j]))
    for (j in 1:parameters$D){
      parameters$beta[j] =  sum(Y[,j]**2) 
      - 2*(Y[,j]%*%parameters$M.Eta)%*%parameters$M.Lambda[j,]
      + sum(diag(EEtaTEta%*%parameters$S.Lambda[,,j])) + crossprod(parameters$M.Lambda[j,], EEtaTEta%*%parameters$M.Lambda[j,])
    }
    parameters$beta = parameters$beta0 + 0.5*parameters$beta
    
    
    # ELBOvec[iter] = ELBO(parameters) #UNCOMMENT TO IMPLEMENT ELBO
    if (iter>5){
      val = ELBOvec[(iter-5):iter]
      if (all(abs(val - mean(val))<tol)){ 
        ELBOvec = ELBOvec[1:iter]
        break }
    }
  }
  
  
  varimaxLambda = varimax(parameters$M.Lambda)$loadings[1:parameters$D,]
  return( list(Lambda=varimaxLambda, ELBO=ELBOvec) )}



# Generate Data
N = 100L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:8, 2] = 1L
Lambda[9:10, 3] = 1L #doesnt work if small??
Sigma = diag(0.5, D)
iN=10
countsmat = matrix(NA, nrow=iN, ncol=2)
for (j in 1:10){
counts = c(0,0)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + Sigma)
for (i in 1:iN){
  tmp  = round(CAVIlotsTau(Y)$Lambda,1)
  print(tmp)
  if (any(colSums(tmp) == 0)){
    counts[2] = counts[2] +1
  } 
  tmp  = round(CAVIfewTau(Y)$Lambda,1)
  print(tmp)
  if (any(colSums(tmp) == 0)){
    counts[1] = counts[1] +1 
  } 
  print(strrep('=',30))
}
countsmat[j,] = counts
print(j)}
countsmat = countsmat/iN
colnames(countsmat) = c("Few Tau", "Lots Tau")
print(countsmat)
