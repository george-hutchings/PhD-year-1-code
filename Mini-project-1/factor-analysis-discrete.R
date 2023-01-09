rm(list = ls()) # Remove variables
# Install a package manager and packages
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(MASS, pracma, RColorBrewer)

## Generate Data
N = 100L
D = 10L
K = 3L
Lambda = matrix(0L, D, K)#DxK
Lambda[1:4, 1] = 1L
Lambda[5:9, 2] = 2L
Lambda[10, 3] = 3L #TODO doesnt work if this is small (eg 1)
Y = mvrnorm(N, rep(0L, D) , tcrossprod(Lambda) + diag(1L,D))

# Convert dimensions to discrete
discreteDims = c(2L, ncol(Y))
boundaries = list(c(-2, -1, 1, 2 ), c(-2, 0, 1 ))
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
#Y = MakeDiscrete(Y, discreteDims, boundaries)



CAVI <- function(Y, maxiterations=1000L, tol = 0.1, seed=NULL){
  if (!is.null(seed)){set.seed(seed)}
  parameters = c()
  parameters$N = nrow(Y)
  parameters$D = ncol(Y)
  parameters$K = K
  parameters$a0 = parameters$b0 = 1e-3
  
  ## Initialise Lambda parameters
  # Initialise mean of Lambda by PCA (and varimax)
  pY = prcomp(Y, scale=TRUE, rank=parameters$K)
  temp = varimax(pY$rotation%*%diag(pY$sdev[1:parameters$K]**2), normalize=FALSE) #TODO normalize??
  parameters$M.Lambda = temp$loadings[1:parameters$D,]
  #Initialised by Wishart distribution, last dimension refers to the vector the covariance matrix corresponds to 
  parameters$S.Lambda =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K)) 
  
  ## Initialise Eta parameters by std normal and Wishart
  parameters$M.Eta = matrix(rnorm(parameters$N*parameters$K), nrow=parameters$N, ncol=parameters$K)
  parameters$S.Eta =  rWishart(parameters$D, parameters$K, diag(1/parameters$K,parameters$K))[,,1] #the same for each observation
  
  ## Initialise b parameters (a does not change)
  parameters$b = matrix(rgamma(parameters$D*parameters$K, shape=2, rate= 1/parameters$b0), nrow=parameters$D, ncol=parameters$K) 
  parameters$a = parameters$a0 + 1
  
  
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
  ELBOvec = numeric(maxiterations)
  for (iter in 1:maxiterations){
    parameters$S.Eta = solve(E.LambdaTLambda(parameters) + diag(1,parameters$K,parameters$K)) # The same for each observation
    parameters$M.Eta = tcrossprod(Y, tcrossprod(parameters$S.Eta, E.Lambda(parameters)))
    
    
    ETau = E.Tau(parameters)
    EEtaTEta = E.EtaTEta(parameters)
    for (i in 1:parameters$D){
      parameters$S.Lambda[,,i] = solve(EEtaTEta + diag(ETau[i,]))
      parameters$M.Lambda[i,] = parameters$S.Lambda[,,i]%*%(crossprod(parameters$M.Eta,Y[,i]))
    }
    
    parameters$b = parameters$b0 + 0.5 *E.Lambda2(parameters)
    
    ELBOvec[iter] = ELBO(parameters)
    if (iter>5){
      val = ELBOvec[(iter-5):iter]
      if (all(abs(val - mean(val))<tol)){ 
        ELBOvec = ELBOvec[1:iter]
        break }
    }
  }
  
  
  varimaxLambda = varimax(parameters$M.Lambda)$loadings[1:parameters$D,]
  return( list(Lambda=varimaxLambda, ELBO=ELBOvec) )}


# Perform CAVI on data repeats times (with random initialisations)
repeats=100L
results = vector(mode = "list", length = repeats)
ELBOlast = numeric(repeats)
for (i in 1:repeats){
  results[[i]] = CAVI(Y, tol=0.01, seed = 1908+i)
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
print(round(VILambda,1))


nonIncreasingELBOidx = which(sapply(results, function(x) is.unsorted(x$ELBO)))
plot(results[[nonIncreasingELBOidx[1]]]$ELBO)
plot(diff(results[[nonIncreasingELBOidx[1]]]$ELBO)>0)




