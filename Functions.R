#--------------------------------------------------------------------------------------------------------------------------------------------
#----                        ----------------------------------------------------------------------------------------------------------------
#----  Simulation Functions  ----------------------------------------------------------------------------------------------------------------
#----                        ----------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Directed ZIPLPCM Simulation
Simulation_Directed_ZIPLPCM <- function(beta,P,mu,tau,d=3,z=NULL,seed=NULL){
  # beta: LPM intercept parameter
  # P: K X K matrix with each entry denoting p_{gh} where K is the number of non-empty groups in z
  # mu: d X K matrix with each column denoting a mu_k
  # tau: 1 X K matrix with each entry denoting tau_k
  # d: dimension of latent positions; d=3 by default
  # z: 1 X N vector; the specified memberships; note that z here should not contain any empty clusters
  
  set.seed(seed)
  
  if (!is.null(z)) {
    N <- length(z)
    K <- max(z)
    Z <- t(t(matrix(z,N,K))==(1:K))*1 # transform membership vector to a N X K matrix
    if (length(which(apply(Z,2,sum)==0))>0){
      stop("Empty groups exist in z!")
    }
  }else{
    stop("The clustering z is not specified!")
  }
  
  # Sample the latent positions
  library(mvtnorm)
  n_k <- table(z) # obtain n_k for each group
  U <- matrix(0,N,d)
  for (k in 1:K){
    U[which(z==k),] <- mvtnorm::rmvnorm(n_k[k],mu[,k],(1/tau[k])*diag(d)) # Rfast also has the Rfast::rmvnorm function which is faster, but it will reset the RNG seed automatically by set.seed(NULL)
  }
  library(Rfast) # for Rfast::Dist()
  Dist_U <- Rfast::Dist(U)
  
  # Sample the augmented unusual zero variable; 1 corresponds to a unusual zero; 0 indicates Poisson distributed
  P_ij <- Z%*%P%*%t(Z)
  nu <- matrix(rbinom(N^2,1,P_ij),N,N)
  diag(nu) <- 0
  
  # Directly sample the realization Y without sampling X
  Y <- matrix(rpois(N^2, exp(beta-Dist_U)),N,N)*(nu==0)
  diag(Y) <- 0
  
  return(list(U=U, Z=Z, z=z, nu=nu, Y=Y))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Directed Poisson-LPCM Simulation
Simulation_Directed_PoissonLPCM <- function(beta,mu,tau,d=3,z=NULL,seed=NULL){
  # beta: LPM intercept parameter
  # P: K X K matrix with each entry denoting p_{gh} where K is the number of non-empty groups in z
  # mu: d X K matrix with each column denoting a mu_k
  # tau: 1 X K matrix with each entry denoting tau_k
  # d: dimension of latent positions; d=2 by default
  # z: 1 X N vector; the specified memberships; note that z here should not contain any empty clusters
  
  set.seed(seed)
  
  if (!is.null(z)) {
    N <- length(z)
    K <- max(z)
    Z <- t(t(matrix(z,N,K))==(1:K))*1 # transform membership vector to a N X K matrix
    if (length(which(apply(Z,2,sum)==0))>0){
      stop("Empty groups exists in z!")
    }
  }else{
    stop("The clustering z is not specified!")
  }
  
  # Sample the latent positions
  library(mvtnorm)
  n_k <- table(z) # obtain n_k for each group
  U <- matrix(0,N,d)
  for (k in 1:K){
    U[which(z==k),] <- mvtnorm::rmvnorm(n_k[k],mu[,k],(1/tau[k])*diag(d))
  }
  library(Rfast) # for Rfast::Dist()
  Dist_U <- Rfast::Dist(U)
  
  # Directly sample the realization Y without sampling X
  Y <- matrix(rpois(N^2, exp(beta-Dist_U)),N,N)
  diag(Y) <- 0
  
  return(list(U=U, Z=Z, z=z, Y=Y))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Directed ZIPLPCM Simulation Process
Simulation_Directed_ZIPSBM <- function(P,Lambda,z=NULL,seed=NULL){
  # beta: LPM intercept parameter
  # P: K X K matrix with each entry denoting p_{gh} where K is the number of non-empty groups in z
  # Lambda: K X K matrix with each entry denoting lambda_{gh} where K is the number of non-empty groups in z
  # z: 1 X N vector; the specified memberships; note that z here should not contain any empty clusters
  
  set.seed(seed)
  
  if (!is.null(z)) {
    N <- length(z)
    K <- max(z)
    Z <- t(t(matrix(z,N,K))==(1:K))*1 # transform membership vector to a N X K matrix
    if (length(which(apply(Z,2,sum)==0))>0){
      stop("Empty groups exists in z!")
    }
  }else{
    stop("The clustering z is not specified!")
  }
  
  # Sample the augmented unusual zero variable; 1 corresponds to a unusual zero; 0 indicates Poisson distributed
  P_ij <- Z%*%P%*%t(Z)
  nu <- matrix(rbinom(N^2,1,P_ij),N,N)
  diag(nu) <- 0
  
  # Construct the matrix of lambda_ij
  Lambda_ij <- Z%*%Lambda%*%t(Z)
  
  # Directly sample the realization Y without sampling X
  Y <- matrix(rpois(N^2, Lambda_ij),N,N)*(nu==0)
  diag(Y) <- 0
  
  return(list(Z=Z, z=z, nu=nu, Y=Y))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

### The following simulation function is for undirected ZIP-LPCM networks

# Undirected ZIPLPCM Simulation
Simulation_UnDirected_ZIPLPCM <- function(beta,P,mu,tau,d=3,z=NULL,seed=NULL){
  # beta: LPM intercept parameter
  # P: K X K matrix with each entry denoting p_{gh} where K is the number of non-empty groups in z
  # mu: d X K matrix with each column denoting a mu_k
  # tau: 1 X K matrix with each entry denoting tau_k
  # d: dimension of latent positions; d=2 by default
  # z: 1 X N vector; the specified memberships; note that z here should not contain any empty clusters
  
  library("Rfast")
  library("gdata") # for upperTriangle() and lowerTriangle() functions
  set.seed(seed)
  
  if (!is.null(z)) {
    N <- length(z)
    K <- max(z)
    Z <- t(t(matrix(z,N,K))==(1:K))*1 # transform membership vector to a N X K matrix
    if (length(which(apply(Z,2,sum)==0))>0){
      stop("Empty groups exist in z!")
    }
  }else{
    stop("The clustering z is not specified!")
  }
  
  # Sample the latent positions
  library(mvtnorm)
  n_k <- Rfast::Table(z) # obtain n_k for each group
  U <- matrix(0,N,d)
  for (k in 1:K){
    U[which(z==k),] <- mvtnorm::rmvnorm(n_k[k],mu[,k],(1/tau[k])*diag(d)) # Rfast also has the Rfast::rmvnorm which is much faster, but it will set.seed(NULL)
  }
  Dist_U <- Rfast::Dist(U)
  
  # Sample the augmented missing zero variable; 1 means missing zero; 0 means Poisson distributed
  P_ij <- Z%*%P%*%t(Z)
  nu <- matrix(0,N,N)
  upperTriangle(nu) <- rbinom(N*(N-1)/2,1,upperTriangle(P_ij))
  lowerTriangle(nu,byrow=TRUE)<-upperTriangle(nu) # or nu <- nu+t(nu)
  
  # Directly sample the realization Y without sampling X
  Y <- matrix(0,N,N)
  upperTriangle(Y) <- rpois(N*(N-1)/2, exp(beta-upperTriangle(Dist_U)))*(upperTriangle(nu)==0)
  lowerTriangle(Y,byrow=TRUE)<-upperTriangle(Y)
  
  return(list(U=U, Z=Z, z=z, nu=nu, Y=Y))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#----                                    ----------------------------------------------------------------------------------------------------
#----  Pre- and Post- Processing Functions  ----------------------------------------------------------------------------------------------------
#----                                    ----------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Label switching function based on  Strehl and J. Ghosh (2003) and Rastelli and Friel (2018)
LabelSwitching_SG2003 <- 
  function(z=NULL, Z=NULL, vector_1byK=NULL,matrix_dbyK=NULL,matrix_KbyK=NULL,matrix2_KbyK=NULL){
    # z: the membership vector to be label-switched
    # Z: the membership matrix to be label-switched
    # vector_1byK: clustering dependent 1 X K vector to be label-switched
    # matrix_dbyK: clustering dependent d X K matrix to be label-switched
    # matrix_KbyK: clustering dependent K X K matrix to be label-switched
    # matrix2_KbyK: the 2nd clustering dependent K X K matrix to be label-switched
    # Note that no empty groups are allowed to exist in this function
    
    if (!is.null(z)&is.null(Z)){
      N <- length(z); K <- max(z)
      Z <- t(t(matrix(z,N,K))==(1:K))*1 # matrix transformation
    }else if(is.null(z)&!is.null(Z)){
      N <- nrow(Z); K <- ncol(Z)
      z <- c(Z%*%1:K)
    }else if(!is.null(z)&!is.null(Z)){
      N <- length(z); K <- ncol(Z)
    }else{
      stop("Either Z or z should be provided!")
    }
    
    if (length(which(apply(Z,2,sum)==0))>0){ # if empty groups exists
      Z <- Z[,-which(apply(Z,2,sum)==0)] # remove the empty groups
      K <- ncol(Z)
      z <- c(Z%*%1:K)
    }
    
    # Construct the reassignment rule which is a K X 2 matrix; example: c(3,1) in a row means the initial group 3 is relabeled as group 1
    reassignment_rule <- cbind(order(apply(Z,2,which.max)),1:K)
    # First create a K X N matrix where each column is reassignment_rule[,1] and we apply t() to obtain the transport N X K matrix.
    # Then we compare each column of the transport with z and replace the output "true" by the corresponding element in the reassignment_rule[,2] N X K transport matrix
    # Finally, we apply row sums to obtain the relabeled clustering
    z <- apply((t(matrix(reassignment_rule[,1],K,N))==z)*t(matrix(reassignment_rule[,2],K,N)),1,sum)
    
    if (!is.null(vector_1byK)){
      vector_1byK <- vector_1byK[reassignment_rule[,1]]
    }
    if (!is.null(matrix_dbyK)){
      matrix_dbyK <- matrix_dbyK[,reassignment_rule[,1]]
    }
    if (!is.null(matrix_KbyK)){
      matrix_KbyK <- matrix_KbyK[reassignment_rule[,1],reassignment_rule[,1]]
    }
    if (!is.null(matrix2_KbyK)){
      matrix2_KbyK <- matrix2_KbyK[reassignment_rule[,1],reassignment_rule[,1]]
    }
    return(list(z=z, Z=Z, vector_1byK=vector_1byK,matrix_dbyK=matrix_dbyK,matrix_KbyK=matrix_KbyK,matrix2_KbyK=matrix2_KbyK))
  }

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Define the non-negative weight, W_{N,K} function in the MFM
MFM_WNK_TruncPois1 <- function(N,K,alpha){
  res <- 0
  k <- K
  repeat{
    Add <- exp(sum(log(k-0:(K-1)))-sum(log(alpha*k+0:(N-1)))-log(factorial(k)*(exp(1)-1)))
    res <- res + Add
    if (Add>.Machine$double.xmin){
      k <- k + 1
    }else{
      break
    }
  }
  return(res)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Evaluating the complete likelihood of (23) in the paper for the directed ZIP-LPCM-MFM
Directed_ZIPLPCM_MFM_CompleteLikelihood <- function(X,U,beta,nu,P,z, alpha1,alpha2,omega, alpha, A=NULL,omega_c=NULL){
  # X: N X N matrix, missing data imputed adjacency matrix
  # U: N X d matrix, latent positions
  # beta: latent position model intercept parameter
  # nu: N X N matrix, unusual zero indicator
  # P: K X K matrix, unusual zero probability for each pair of clusters
  # z: 1 X N vector, clustering without empty clusters
  
  # alpha_1,alpha_2,omega: prior parameters of MVN(u_i|...)
  # alpha: MFM parameter
  
  # A: 1 X N vector, categorical node attributes, this correspond to the node attributes \boldsymbol{c} in the paper
  # omega_c: node attributes cohesion parameter
  
  library(Rfast) # for Rfast::Dist()
  
  K <- max(z) # extract the number of non-empty clusters
  N <- length(z) # extract the network size
  d <- ncol(U) # extract the dimension of the latent positions
  res <- 0 # initialize the output
  
  # Evaluate the X|U,\beta term
  Dist_U <- Rfast::Dist(U)
  res_temp <- dpois(X,exp(beta-Dist_U),log=TRUE)
  diag(res_temp) <- 0
  res <- res + sum(res_temp)
  
  # Evaluate the \nu|P,z term
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
  res_temp <- dbinom(nu,1,P_ij,log=TRUE)
  diag(res_temp) <- 0
  res <- res + sum(res_temp)
  
  # Evaluate the U|z term
  res <- res + (alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega))*K
  n_k_table <- table(z)
  for (k in 1:K){
    U_k <- matrix(U[z==k,],ncol=d) # matrix() function here ensures matrix form of one single vector case
    n_k <- n_k_table[[k]]
    res <- res + lgamma(alpha1+n_k*d/2)-n_k*d/2*log(pi)-d/2*log(omega+n_k)-
      (alpha1+n_k*d/2)*log(sum(U_k^2)-sum(colSums(U_k)^2)/(n_k+omega)+alpha2)
  }
  
  if (!is.null(A)){ # if node attributes are provided
    # Evaluate the A|z term
    C <- max(A)
    table.z.A <- table(z,A)
    sum_omega_c <- omega_c*C # note that no empty level is allowed in A for c = 1,2,...,C
    res <- res + sum(lgamma(table.z.A+omega_c)) - sum(lgamma(n_k_table+sum_omega_c)) + (lgamma(sum_omega_c)-C*lgamma(omega_c))*K
  }
  
  # Evaluate the z term
  res <- res + log(MFM_WNK_TruncPois1(N,K,alpha))
  for (k in 1:K){
    res <- res + sum(log(alpha+0:(n_k_table[[k]]-1)))
  }
  
  return(res)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Evaluating the complete likelihood of the directed Poisson-LPCM-MFM
# This is similar to the ZIP-LPCM-MFM, but excluding the \nu|P,z term

Directed_PoissonLPCM_MFM_CompleteLikelihood <- function(X,U,beta,z, alpha1,alpha2,omega, alpha, A=NULL,omega_c=NULL){
  # X: N X N matrix, observed adjacency matrix
  # U: N X d matrix, latent positions
  # beta: latent position model intercept parameter
  # z: 1 X N vector, clustering without empty clusters
  
  # alpha_1,alpha_2,omega: prior parameters of MVN(u_i|...)
  # alpha: MFM parameter
  
  # A: 1 X N vector, categorical node attributes, this correspond to the node attributes \boldsymbol{c} in the paper
  # omega_c: node attributes cohesion parameter
  
  library(Rfast) # for Rfast::Dist()
  
  K <- max(z) # extract the number of non-empty clusters
  N <- length(z) # extract the network size
  d <- ncol(U) # extract the dimension of the latent positions
  res <- 0 # initialize the output
  
  # Evaluate the X|U,\beta term
  Dist_U <- Rfast::Dist(U)
  res_temp <- dpois(X,exp(beta-Dist_U),log=TRUE)
  diag(res_temp) <- 0
  res <- res + sum(res_temp)
  
  # Evaluate the U|z term
  res <- res + (alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega))*K
  n_k_table <- table(z)
  for (k in 1:K){
    U_k <- matrix(U[z==k,],ncol=d) # matrix() function here ensures matrix form of one single vector case
    n_k <- n_k_table[[k]]
    res <- res + lgamma(alpha1+n_k*d/2)-n_k*d/2*log(pi)-d/2*log(omega+n_k)-
      (alpha1+n_k*d/2)*log(sum(U_k^2)-sum(colSums(U_k)^2)/(n_k+omega)+alpha2)
  }
  
  if (!is.null(A)){ # if node attributes are provided
    # Evaluate the A|z term
    C <- max(A)
    table.z.A <- table(z,A)
    sum_omega_c <- omega_c*C # note that no empty level is allowed in A for c = 1,2,...,C
    res <- res + sum(lgamma(table.z.A+omega_c)) - sum(lgamma(n_k_table+sum_omega_c)) + (lgamma(sum_omega_c)-C*lgamma(omega_c))*K
  }
  
  # Evaluate the z term
  res <- res + log(MFM_WNK_TruncPois1(N,K,alpha))
  for (k in 1:K){
    res <- res + sum(log(alpha+0:(n_k_table[[k]]-1)))
  }
  
  return(res)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Evaluating the complete likelihood of (23) in the paper for the undirected ZIP-LPCM-MFM

UnDirected_ZIPLPCM_MFM_CompleteLikelihood <- function(X,U,beta,nu,P,z, alpha1,alpha2,omega, alpha, A=NULL,omega_c=NULL){
  # X: N X N matrix, missing data imputed adjacency matrix
  # U: N X d matrix, latent positions
  # beta: latent position model intercept parameter
  # nu: N X N matrix, unusual zero indicator
  # P: K X K matrix, unusual zero probability for each pair of clusters
  # z: 1 X N vector, clustering without empty clusters
  
  # alpha_1,alpha_2,omega: prior parameters of MVN(u_i|...)
  # alpha: MFM parameter
  
  # A: 1 X N vector, categorical node attributes, this correspond to the node attributes \boldsymbol{c} in the paper
  # omega_c: node attributes cohesion parameter
  
  library("gdata")
  library("Rfast")
  
  K <- max(z) # extract the number of non-empty clusters
  N <- length(z) # extract the network size
  d <- ncol(U) # extract the dimension of the latent positions
  res <- 0 # initialize the output
  
  # Evaluate the X|U,\beta term
  Dist_U <- Rfast::Dist(U)
  res_temp <- dpois(upperTriangle(X),exp(beta-upperTriangle(Dist_U)),log=TRUE)
  res <- res + sum(res_temp)
  
  # Evaluate the \nu|P,z term
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
  res_temp <- dbinom(upperTriangle(nu),1,upperTriangle(P_ij),log=TRUE)
  res <- res + sum(res_temp)
  
  # Evaluate the U|z term
  res <- res + (alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega))*K
  n_k_table <- Rfast::Table(z)
  U_special_matrix <- matrix(U,nrow=N,ncol=d*K)*matrix(do.call("rbind", replicate(d,Z,simplify = FALSE)),nrow=N,ncol=d*K)
  res <- res + sum(lgamma(alpha1+n_k_table*d/2)-n_k_table*d/2*log(pi)-d/2*log(omega+n_k_table)-
                     (n_k_table*d/2+alpha1)*log(rowSums(t(Z)%*%(U_special_matrix^2))-
                                                  rowSums((t(Z)%*%U_special_matrix)^2)/(n_k_table+omega)+alpha2))
  
  if (!is.null(A)){ # if node attributes are provided
    # Evaluate the A|z term
    C <- max(A)
    table.z.A <- table(z,A) # Rfast::Table() doesn't work here
    sum_omega_c <- omega_c*C # Note that no empty level is allowed in A for c = 1,2,...,C
    res <- res + sum(lgamma(table.z.A+omega_c)) - sum(lgamma(n_k_table+sum_omega_c)) + (lgamma(sum_omega_c)-C*lgamma(omega_c))*K
  }
  
  # Evaluate the z term
  res <- res + log(MFM_WNK_TruncPois1(N,K,alpha))
  for (k in 1:K){
    res <- res + sum(log(alpha+0:(n_k_table[[k]]-1)))
  }
  
  return(res)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#----                            ------------------------------------------------------------------------------------------------------------
#----  Implementation Functions  ------------------------------------------------------------------------------------------------------------
#----                            ------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

### We start from directed cases

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# The partially collapsed Metropolis-within-Gibbs algorithm for directed ZIP-LPCM-MFM 
MwG_Directed_ZIPLPCM <- function(Y,T, omega,alpha1,alpha2,alpha, beta1,beta2,sigma2prop_beta,sigma2prop_U,d=3,z=NULL,p_eject=0.5,A=NULL,omega_c=NULL){
  library(mvtnorm) # for rmvnorm
  library(ergmito) # for geodesic
  library(Rfast) # for Rfast::Dist()
  
  N <- nrow(Y) # Extract the number of nodes
  log_p_0 <- function(a,n_k){lgamma(2*a)-lgamma(a)+lgamma(a+n_k)-lgamma(2*a+n_k)} # build the log look-up table for choosing Beta(a,a) in the truncated AE step
  Set_a_for_P_E <- Vectorize(log_p_0,c("a","n_k")) # vectorize the variables a and n_k
  lookuptable_for_a <- outer(seq(0.01,100,0.01),1:N,Set_a_for_P_E) # evaluate for each pair of a and n_k where a=seq(0.01,100,0.01) and n_k=0:N
  # lookuptable_for_a[is.nan(lookuptable_for_a)]<-0 # transform NaN to 0 # no NaN in log form
  # Note that we set logP(n_g=0)=logP(n_{K+1}=0)=log(0.01) in this algorithm
  
  # Y: The N X N observed adjacency matrix
  # A: The 1 X N categorical exogenous node attributes which correspond to \boldsymbol{c} in the paper
  # omega_c: The parameter of DM cohesion for A where omega_c is assumed to be the same for all c = 1,2,...,C
  if (!is.null(A)){ # if A is provided
    C <- max(A)
    if (!is.null(omega_c)){
      sum_omega_c <- omega_c*C # Note that no empty level is allowed in A for c = 1,2,...,C
    }else{
      stop("omega_c should be provided if node attributes is included!")
    }
  }
  
  # Recall here that
  # beta: a LPM intercept parameter
  # P: a K X K matrix of clustering dependent unusual zero probability
  # U: N X d matrix, latent positions
  # nu: N X N matrix, unusual zero indicator
  # z: 1 X N vector, memberships
  
  # mu: d X K matrix MVN mean (collapsed, may be inferred if required)
  # tau: 1 X K vector MVN precision, tau = 1/sigma2 (collapsed, may be inferred if required)
  # Pi: 1 X K vector membership probability (collapsed, may be inferred if required)
  
  # omega, (alpha1,alpha2): conjugate prior parameters for mu_k and tau_k: k=1,2,...,K
  # alpha: prior parameter for Pi; (beta1,beta2): prior parameter for each p_gh: g,h=1,2,...,K
  # sigma2prop_U,sigma2prop_beta: the proposal variance proposed, respectively, for U and beta
  # d: the dimension of the latent positions
  # z: initial clustering specified
  
  # Note that K is the number of non-empty groups, and no empty groups is allowed in z within this function
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize the clustering
  if (is.null(z)){ # if initial clustering is not provided
    z <- 1:N # set the clustering where each node is in a singleton group as the initial clustering
  }
  K <- max(z)
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize beta and P
  beta <- 0
  P <- matrix(beta1/(beta1+beta2),K,K) # set p_gh prior mean for initial state of each p_gh: g,h=1,2,...,K
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize latent positions by Classical (Metric) Multidimensional Scaling (MDS) on asymmetric distance matrix
  Y_GM <- geodesic(Y)[[1]] # geodesic matrix (GM) of Y, the length/the number of edges of the shortest paths between each pair of nodes
  Y_GM[is.na(Y_GM)] <- 0 # remove NA as 0
  U <- cmdscale(Y_GM,d) # obtain the initial latent positions by MDS
  # Note here that MDS may bring the same position for multiple nodes
  Dist_U <- Rfast::Dist(U)
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize the missing zero indicator matrix nu and the missing data imputed adj X
  nu <- (Y==0)*1
  diag(nu) <- 0
  Imputation_position <- which(nu==1) # extract the position of those y_{ij}=0 except the diagonal elements
  Imputation_n <- length(Imputation_position) # the number of possible missing data imputed elements
  
  P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
  Imp_P_ij <- P_ij[Imputation_position] # extract the p_ij for those possibly imputed elements
  nu[Imputation_position] <- rbinom(Imputation_n,1,Imp_P_ij/(Imp_P_ij + (1-Imp_P_ij)*dpois(0,exp(beta-Dist_U[Imputation_position])))) # Initialize nu
  
  X <- Y
  nu_1_position <- which(nu==1) # extract the position of those updated nu_ij=1
  X[nu_1_position] <- rpois(length(nu_1_position),exp(beta-Dist_U[nu_1_position])) # Initialize X
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Posterior chains
  z_list <- matrix(z,nrow=1) # the posterior chain will be a (T+1) X N matrix whose (t+1)th row is the t'th clustering z^{(t)}
  U_list <- list(U) # the posterior chain will be a list()
  nu_list <- list(nu) # the posterior chain will be a list()
  beta_list <- c(beta) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th beta^{(t)}
  P_list <- list(P) # the posterior chain will be a list()
  K_list <- c(K) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th K^{(t)}
  # Note that K here is the number of non-empty clusters
  X_list <- list(X) # the posterior chain will be a list()
  acceptance_count_U <- matrix(0,T,N) # Count the number of acceptance for the M-H step of each latent position u_i
  acceptance_count_beta <- matrix(0,1,T) # Count the number of acceptance for the M-H step of beta
  acceptance_count_TAE <- matrix(0,1,T) # Count the number of acceptance for the M-H step of AE move
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Main algorithm
  for (t in 1:T){
    if ((t%%100) == 0){
      cat("t=",t,", K=",K,"\n") # monitor the process
    }
    # Note that (.)^{(t)} is the current state and (.)^{(t+1)}is the new state within this function
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sample nu
    P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
    Imp_P_ij <- P_ij[Imputation_position] # extract the p_ij for those possibly imputed elements, i.e. all the y_ij=0 except diagonal
    nu[Imputation_position] <- rbinom(Imputation_n,1,Imp_P_ij/(Imp_P_ij + (1-Imp_P_ij)*dpois(0,exp(beta-Dist_U[Imputation_position]))))
    nu_list[[t+1]] <- nu
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sample X
    X <- Y
    nu_1_position <- which(nu==1) # extract the position of those updated nu_ij=1
    X[nu_1_position] <- rpois(length(nu_1_position),exp(beta-Dist_U[nu_1_position]))
    X_list[[t+1]] <- X
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update beta by M-H step
    beta_c <- beta # reserve the current beta
    beta <- rnorm(1,beta_c,sd=sqrt(sigma2prop_beta)) # Note here that sigma2prop_beta is the proposal variance
    log_alpha_beta_right <- dpois(X,exp(beta-Dist_U),log=TRUE)-dpois(X,exp(beta_c-Dist_U),log=TRUE)
    diag(log_alpha_beta_right) <- 0
    log_alpha_beta_right <- sum(log_alpha_beta_right)
    
    if (log(runif(1)) <= min(0,log_alpha_beta_right)){ # accept or not
      beta_list[t+1] <- beta
      acceptance_count_beta[t] <- 1
    }else{
      beta_list[t+1] <- beta <- beta_c
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each u_i by M-H step
    n_k <- table(z) # n_k of z
    position_matrix <- Z==1 # return TRUE if z_ik = 1 and FALSE otherwise in Z matrix
    for (i in sample.int(N,N)){ # update in random order
      U_c <- U # Update the current latent positions
      Dist_U_c <- Dist_U # Update the Dist_U_c based on U_c
      n_zi <- n_k[[z[i]]] # extract the current number of nodes in group z_i;
      position_zi <- position_matrix[,z[i]] # z==z[i] # extract the positions within 1:N of those nodes from group z_i
      U_c_zi <- matrix(U_c[position_zi,],ncol=d) # extract the current latent positions of nodes in group z_i
      # Note that the matrix(,ncol=d) transformation here aims to deal with the single vector case which would not be in the form of a 1 X 2 matrix
      U[i,] <- mvtnorm::rmvnorm(1,U_c[i,],sigma2prop_U*diag(d)) # proposal distribution # Rfast also has the Rfast::rmvnorm which is faster, but it will reset the RNG seed by set.seed(NULL)
      Dist_U <- Rfast::Dist(U) # Update the Dist_U based on proposed U
      U_p_zi <- matrix(U[position_zi,],ncol=d) # extract the proposed latent positions of nodes in group z_i
      log_alpha_U_right <-  (n_zi*d/2 + alpha1)*(log(sum(U_c_zi^2) - sum(colSums(U_c_zi)^2)/(n_zi + omega) + alpha2)-
                                                   log(sum(U_p_zi^2) - sum(colSums(U_p_zi)^2)/(n_zi + omega) + alpha2))+ # f(U|z) terms
        sum(dpois(X[i,-i],exp(beta-Dist_U[i,-i]),log=TRUE)+dpois(X[-i,i],exp(beta-Dist_U[-i,i]),log=TRUE)-
              dpois(X[i,-i],exp(beta-Dist_U_c[i,-i]),log=TRUE)-dpois(X[-i,i],exp(beta-Dist_U_c[-i,i]),log=TRUE)) # f(X|beta,U) terms
      
      if (log(runif(1)) <= min(0,log_alpha_U_right)){ # accept or not
        acceptance_count_U[t,i] <- 1
      }else{
        U <- U_c
        Dist_U <- Dist_U_c
      }# end if
    }# end i
    U_list[[t+1]] <- U
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each z_i via MFM
    # Note that this MFM inference step has no empty clusters
    nu_gh <- t(Z)%*%nu%*%Z # {nu_gh} of z
    if (length(n_k)>1){
      n_gh <- outer(n_k,n_k)-diag(n_k) # {n_gh} of z; n_k already obtained above
    }else{ # if z only has one cluster
      n_gh <- matrix(n_k^2-n_k)
    }
    
    for (i in sample.int(N,N)){
      nonempty <- which(colSums(matrix(Z[-i,],nrow=N-1)) > 0) # extract non-empty groups after removing node i
      Z <- Z[,nonempty] # remove the group if removing node i makes it empty
      if (length(nonempty)==1){Z <- matrix(Z,N,1)} # ensure the matrix form of Z when K = 1 after removing the possible empty group
      K <- ncol(Z) # update K after removing node i
      Z_Ri <- Z[-i,] # remove node i, i.e., the i'th row from Z
      if (K==1){Z_Ri <- matrix(Z_Ri,nrow=N-1,ncol=1)} # Z_Ri becomes vector form if k=1, so we transform it to matrix to avoid crashing
      z_Ri <- c(Z_Ri%*%1:K) # z after removing node i
      U_Ri <- U[-i,] # U after removing node i
      n_k_Ri_full <- colSums(Z_Ri) # an 1 X K matrix of {n_k} after removing node i
      
      # 1. First obtain f(U|z) terms of the full-conditional distribution by assuming that node i is removed and then is reassigned
      # Construct a N-1 X d*K matrix of latent positions where 
      # the rows which correspond to cluster 1 individuals in the 1st three columns store the corresponding latent positions, and other elements in the 1st three columns are zeros
      # The rows which correspond to cluster 2 individuals in the 2nd three columns store the corresponding latent positions, and other elements in the 2nd three columns are zeros
      # Similarly for all the rest rows and columns
      U_special_matrix_Ri <- matrix(U_Ri,nrow=N-1,ncol=d*K)*matrix(do.call("rbind", replicate(d,Z_Ri,simplify = FALSE)),nrow=N-1,ncol=d*K)
      # Use rbind() function to combine the above matrix with a new K X d*K matrix, where matrix[1,1:3], matrix[2,3:6],matrix[3,7:9],... of this new matrix are U[i,] if d=3 while all other elements are zeros
      U_special_matrix_Ai <- rbind(U_special_matrix_Ri,t(t(matrix(do.call("rbind", replicate(d,diag(K),simplify = FALSE)),nrow=K,ncol=d*K))*U[i,]))
      Z_Ri_spcecial_Ai <- rbind(Z_Ri,diag(K)) # create a special clustering based on Z_Ri, where each cluster is added by one more node whose latent position is U[i,]
      n_k_Ri_full_Ai <- n_k_Ri_full+1 # an 1 X K matrix of {n_k} after assuming that the node i is added to each cluster
      
      # Ri means removing node i; evaluate each non-empty k'th component of the full-conditional function U\{u_i} | z\{z_i}
      # Ai means adding node i; evaluate each non-empty k'th component of the full-conditional function U|z by assuming that node i is assigned in each non-empty cluster
      Uk_terms_Ai <- lgamma(alpha1+n_k_Ri_full_Ai*d/2)-n_k_Ri_full_Ai*d/2*log(pi)-d/2*log(omega+n_k_Ri_full_Ai)-
        (n_k_Ri_full_Ai*d/2+alpha1)*log(rowSums(t(Z_Ri_spcecial_Ai)%*%(U_special_matrix_Ai^2))-
                                          rowSums((t(Z_Ri_spcecial_Ai)%*%U_special_matrix_Ai)^2)/(n_k_Ri_full_Ai+omega)+alpha2)
      Uk_terms_Ri <- lgamma(alpha1+n_k_Ri_full*d/2)-n_k_Ri_full*d/2*log(pi)-d/2*log(omega+n_k_Ri_full)-
        (n_k_Ri_full*d/2+alpha1)*log(rowSums(t(Z_Ri)%*%(U_special_matrix_Ri^2))-
                                       rowSums((t(Z_Ri)%*%U_special_matrix_Ri)^2)/(n_k_Ri_full+omega)+alpha2)
      log_update_prob_U <- c(Uk_terms_Ai-Uk_terms_Ri,lgamma(alpha1+d/2)-d/2*log(pi)-d/2*log(omega+1) + alpha1*log(alpha2)+d/2*log(omega)-lgamma(alpha1)-
                               (d/2+alpha1)*log(sum(U[i,]^2)*(1-1/(1+omega))+alpha2)) # here we evaluate the f(U|z\{z_i},z_i)/f(U\{u_i}|z\{z_i}) as well as the case which assigns node i to the new group K+1
      # If we let U_k be the set of latent positions of nodes from cluster k, then some of the terms above correspond to:
      # rowSums(t(Z_Ri)%*%(U_special_matrix_Ri^2)) is {sum(U_k^2)} before adding node i to each cluster
      # rowSums((t(Z_Ri)%*%U_special_matrix_Ri)^2) is {sum(colSums(U_k)^2)} before adding node i to each cluster
      # rowSums(t(Z_Ri_spcecial_Ai)%*%(U_special_matrix_Ai^2)) is {sum(U_k^2)} after adding node i to each cluster
      # rowSums((t(Z_Ri_spcecial_Ai)%*%U_special_matrix_Ai)^2) is {sum(colSums(U_k)^2)} after adding node i to each cluster
      
      # 2. Then we obtain the f(nu|z) term
      nu_gh_Ri <- matrix(nu_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      n_gh_Ri <- matrix(n_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      nu_ik <- crossprod(Z_Ri, nu[i,-i]) # {nu_ik}
      nu_ki <- crossprod(Z_Ri, nu[-i,i]) # {nu_ki}
      z_i <- which(Z[i,] > 0) # this is z_i, e.g., = k; or = integer(0) with length==0 if removing node i leading to an empty cluster
      if (length(z_i)==1){ # if removing node i doesn't bring an empty cluster
        nu_gh_resid <- matrix(0,K,K) # calculate the difference in nu_gh after removing node i
        nu_gh_resid[z_i,] <- nu_ik; nu_gh_resid[,z_i] <- nu_gh_resid[,z_i]+nu_ki
        nu_gh_Ri <- nu_gh_Ri-nu_gh_resid
        
        n_gh_resid <- matrix(0,K,K) # calculate the difference in n_gh after removing node i
        n_gh_resid[z_i,] <- n_k_Ri_full; n_gh_resid[,z_i] <- n_gh_resid[,z_i]+n_k_Ri_full
        n_gh_Ri <- n_gh_Ri-n_gh_resid
      }
      
      if (K>1){
        nu_ik_row_Matrix <- matrix(1,K,1)%*%matrix(nu_ik,1,K) # construct a K X K matrix where each row is nu_ik
        nu_ki_column_Matrix <- matrix(nu_ki,K,K) # construct a K X K matrix where each column is nu_ki
        n_k_Ri_full_row_Matrix <- matrix(1,K,1)%*%matrix(n_k_Ri_full,1,K) # construct a K X K matrix where each row is n_k_Ri_full; thus the transport is a K X K matrix where each column is n_k_Ri_full
        
        changed_nu_gh_row <- nu_gh_Ri+nu_ik_row_Matrix+diag(diag(nu_ki_column_Matrix)) # construct a K X K matrix where each row is the row of corresponding nu_gh after adding node i
        changed_nu_gh_column <- nu_gh_Ri+nu_ki_column_Matrix # construct a K X K matrix where each column is the column of corresponding nu_gh after adding node i but excluing diagonal
        diag(changed_nu_gh_column) <- NA # excluding diagonal
        changed_nu_gh_column <- matrix(changed_nu_gh_column[!is.na(changed_nu_gh_column)],nrow=K-1,ncol=K) # excluding diagonal leading to a K-1 X K matrix
        
        changed_n_gh_row <- n_gh_Ri+n_k_Ri_full_row_Matrix+diag(n_k_Ri_full) # construct a K X K matrix where each row is the row of corresponding n_gh after adding node i
        changed_n_gh_column <- n_gh_Ri+t(n_k_Ri_full_row_Matrix) # construct a K X K matrix where each column is the column of corresponding n_gh after adding node i but excluing diagonal
        diag(changed_n_gh_column) <- NA # excluding diagonal
        changed_n_gh_column <- matrix(changed_n_gh_column[!is.na(changed_n_gh_column)],nrow=K-1,ncol=K) # excluding diagonal leading to a K-1 X K matrix
        
        log_update_prob_nu <- c(rowSums(lbeta(changed_nu_gh_row+beta1,changed_n_gh_row-changed_nu_gh_row+beta2))+ # add the row terms including diag for each z_i=k case
                                  colSums(lbeta(changed_nu_gh_column+beta1,changed_n_gh_column-changed_nu_gh_column+beta2))- # add the column terms excluding diag for each z_i=k case
                                  rowSums(lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2))-colSums(lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2))+ # evaluate the demominator of the f(\nu|z_{-i},z_i=k)/f(\nu_{-i}|z_{-i}) for each k = 1,2,...,K
                                  lbeta(diag(nu_gh_Ri)+beta1,diag(n_gh_Ri-nu_gh_Ri)+beta2),# remove the repeat diag term
                                sum(lbeta(c(nu_ik,nu_ki)+beta1,c(n_k_Ri_full,n_k_Ri_full)-c(nu_ik,nu_ki)+beta2)-lbeta(beta1,beta2))) # add the term relevant to the new group K+1
      }else{ # if removing node i leads to only 1 group
        log_update_prob_nu <- c(lbeta(nu_gh_Ri+nu_ki+nu_ik+beta1,n_gh_Ri+2*n_k_Ri_full-(nu_gh_Ri+nu_ki+nu_ik)+beta2)-
                                  lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2),
                                sum(lbeta(c(nu_ik,nu_ki)+beta1,c(n_k_Ri_full,n_k_Ri_full)-c(nu_ik,nu_ki)+beta2)-lbeta(beta1,beta2))) # add the term relevant to the new group K+1
      }
      
      # 3. Obtain the f(z) term
      log_update_prob_z <- c(log(n_k_Ri_full+alpha),log(alpha)+log(MFM_WNK_TruncPois1(N,K+1,alpha))-
                               log(MFM_WNK_TruncPois1(N,K,alpha)))
      
      # 4. Obtain the assignment probability
      log_update_prob <- log_update_prob_U+log_update_prob_nu+log_update_prob_z
      if (!is.null(A)){ # if A is provided
        A_Ri <- A[-i] # remove c_i from node attributes A
        table.z_Ri.A_Ri <- table(z_Ri,factor(A_Ri,levels=1:C)) # it's possible that an empty level might exist in A_Ri, so we use factor(,levels) to avoid such a problem
        log_update_prob <- log_update_prob + c(log(table.z_Ri.A_Ri[,A[i]]+omega_c)-log(n_k_Ri_full+sum_omega_c),log(omega_c)-log(sum_omega_c))
      }
      
      update_prob <- exp(log_update_prob - max(log_update_prob))
      
      # 5. Determine the reassignment and update Z,K,nu_gh,n_gh,n_k
      z_i_new <- which(rmultinom(1,1,update_prob) > 0)
      if(length(z_i) == 1){Z[i,z_i]<-0} # if removing node i doesn't bring an empty group, let Z[i,] be zeros; otherwise, Z[i,] are already zeros
      if (z_i_new==K+1){
        Z <- cbind(Z,0)
        K <- K+1
        Z[i,z_i_new] <- 1
        nu_gh_Ri <- rbind(cbind(nu_gh_Ri,0),0)
        n_gh_Ri <- rbind(cbind(n_gh_Ri,0),0)
        nu_ik <- c(nu_ik,0); nu_ki <- c(nu_ki,0)
        n_k_Ri_full <- c(n_k_Ri_full,0)
      }else{
        Z[i,z_i_new] <- 1
      }
      
      nu_gh_resid <- matrix(0,K,K) # nu_gh difference by adding node i to group k
      nu_gh_resid[z_i_new,] <- nu_ik; nu_gh_resid[,z_i_new] <- nu_gh_resid[,z_i_new]+nu_ki
      nu_gh <- nu_gh_Ri + nu_gh_resid
      
      n_gh_resid <- matrix(0,K,K) # n_gh difference by adding node i to group k
      n_gh_resid[z_i_new,] <- n_k_Ri_full; n_gh_resid[,z_i_new] <- n_gh_resid[,z_i_new]+n_k_Ri_full
      n_gh <- n_gh_Ri + n_gh_resid
    } # end for i
    
    z <- c(Z%*%c(1:ncol(Z)))
    z_list <- rbind(z_list,z)
    K_list[t+1] <- K
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # TAE step
    if (K == 1){ # if only one cluster, eject
      AE_decision <- 0
      p_eject_p <- 1 # proposal probability of the eject move
      p_absorb_p <- 1-p_eject # proposal probability of the reverse absorb move
    }else if (K == N){ # if N clusters, absorb
      AE_decision <- 1
      p_eject_p <- p_eject
      p_absorb_p <- 1
    }else{ # with probability p_eject to eject or absorb
      AE_decision <- runif(1)
      p_eject_p <- p_eject # P(eject) in the proposal
      p_absorb_p <- 1-p_eject
    }
    
    if (AE_decision<=p_eject){ # eject move
      g <- sample.int(K,1) # randomly choose one cluster
      n_k <- table(z) # n_k of z
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      U_g_c <- matrix(U[z==g,],ncol=d) # extract all the u_j for current clustering z_j=g
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_c]-log(0.01)))] # Recall: check "a" ranging from 0.01 to 100 and check n from 1 to N; 0.01 here is set as the Pr(n_{K+1}=0);
      p_E <- rbeta(1,a,a) # sample the probability of the nodes being assigned to cluster K+1
      K_p <- K+1 # proposed number of clusters
      z_p <- z # initialize proposed clustering
      z_p[z_p==g]<-rbinom(n_g_c,1,p_E)*(K_p-g)+g # apply the eject move; for each node j in cluster g, if rbinom() returns 1, assign (K_p-g)+g to z_j; if rbinom() returns 0, assign g back to z_j
      n_k_p <- table(factor(z_p,levels=1:K_p)) # {n_k} of z_p
      n_g_p <- n_k_p[[g]] # extract n_g for proposed clustering after eject move
      # Note that it's possible to have empty cluster in the proposed z and table() can not directly return 0 for empty clusters;
      # Thus we use table(factor(z_p,levels=1:K_p)) to let the empty cluster appear.
      if (n_g_p!=0 & n_g_p!=n_g_c){ # if any of the ejected cluster is not empty
        if (K_p==N){ # If the proposed state has K == N
          p_absorb_p <- 1 # the probability of the backward absorb move for the proposed state is 1
        }
        U_g_p <- matrix(U[z_p==g,],ncol=d) # extract all the u_j for proposed clustering z_j=g after eject move
        n_Kp_p <- n_k_p[[K_p]] # extract n_Kp for proposed clustering after eject move
        U_Kp_p <- matrix(U[z_p==K_p,],ncol=d) # extract all the u_j for proposed clustering z_j=K_p after eject move
        log_alpha_E_right <-  alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega)+log(p_absorb_p)-log(p_eject_p)+ # constant terms
          lgamma(alpha1+n_g_p*d/2) - n_g_p*d/2*log(pi) - d/2*log(omega+n_g_p) - (alpha1+n_g_p*d/2)*log(sum(U_g_p^2)-sum(colSums(U_g_p)^2)/(n_g_p+omega)+alpha2)+ # proposed lambda'_g term
          lgamma(alpha1+n_Kp_p*d/2) - n_Kp_p*d/2*log(pi) - d/2*log(omega+n_Kp_p) - (alpha1+n_Kp_p*d/2)*log(sum(U_Kp_p^2)-sum(colSums(U_Kp_p)^2)/(n_Kp_p+omega)+alpha2)- # proposed lambda'_Kp term
          (lgamma(alpha1+n_g_c*d/2) - n_g_c*d/2*log(pi) - d/2*log(omega+n_g_c) - (alpha1+n_g_c*d/2)*log(sum(U_g_c^2)-sum(colSums(U_g_c)^2)/(n_g_c+omega)+alpha2))+ # current lambda^{(t+1)]}_g term; note that z^{(t+1)} was updated from z^{(t)} in the previous step, so current state is z^{(t+1)}
          log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))+sum(log(alpha+0:(n_Kp_p-1)))- # z' MFM prior term
          log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))+ # z MFM prior term
          2*lgamma(a)-lgamma(2*a)+lgamma(2*a+n_g_c)-lgamma(a+n_g_p)-lgamma(a+n_Kp_p)+ # proposal terms
          log(1-2*exp(log_p_0(a,n_g_c))) # normalizing term of truncating at the empty cluster cases
        
        if (!is.null(A)){ # if A is provided
          table.z_c.A <- table(z,A)
          table.z_p.A <- table(z_p,A)
          log_alpha_E_right <- log_alpha_E_right +
            sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)+
            sum(lgamma(table.z_p.A[K_p,]+omega_c))-lgamma(n_Kp_p+sum_omega_c)+lgamma(sum_omega_c)-C*lgamma(omega_c)-
            sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c) # the ratio of A|z term
        }
        
        Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
        nu_gh_p <- t(Z_p)%*%nu%*%Z_p # {nu_gh} of z_p
        n_gh_p <- outer(n_k_p,n_k_p)-diag(n_k_p) # {n_gh} of z_p
        
        changed_gh_position <- matrix(FALSE,K,K) # extract those changed blocks
        changed_gh_position[g,] <- TRUE; changed_gh_position[,g] <- TRUE
        changed_gh_position_p <- matrix(FALSE,K_p,K_p) # extract those changed blocks
        changed_gh_position_p[c(g,K_p),] <- TRUE; changed_gh_position_p[,c(g,K_p)] <- TRUE
        
        changed_nu_gh <- nu_gh[changed_gh_position]; changed_n_gh <- n_gh[changed_gh_position]
        changed_nu_gh_p <- nu_gh_p[changed_gh_position_p]; changed_n_gh_p <- n_gh_p[changed_gh_position_p]
        
        log_alpha_E_right <- log_alpha_E_right + sum(lbeta(changed_nu_gh_p+beta1,changed_n_gh_p-changed_nu_gh_p+beta2))-
          sum(lbeta(changed_nu_gh+beta1,changed_n_gh-changed_nu_gh+beta2))-lbeta(beta1,beta2)*(sum(changed_gh_position_p)-sum(changed_gh_position))
        
        if (log(runif(1)) <= min(0,log_alpha_E_right)){ # accept or not
          z_list[t+1,] <- z <- z_p # update the new clustering
          Z <- Z_p
          K_list[t+1] <- K <- K_p
          acceptance_count_TAE[t] <- 1 # 1 means accept, 0 means not accept
          nu_gh <- nu_gh_p
          n_gh <- n_gh_p
        }# else: stay at the current state and do nothing
      }# else: stay at the current state and do nothing
    }else{ # absorb step
      random_two_picks <- sort(sample.int(K,2))
      g <- random_two_picks[1] # g < h, we let g absorb h
      h <- random_two_picks[2]
      n_k <- table(z)
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      n_h_c <- n_k[[h]] # extract n_h for current clustering
      U_g_c <- matrix(U[z==g,],ncol=d) # extract all the u_j for current clustering z_j=g
      U_h_c <- matrix(U[z==h,],ncol=d) # extract all the u_j for current clustering z_j=h
      K_p <- K-1
      if (K_p==1){ # If the proposed state has K == 1
        p_eject_p <- 1 # the probability of the backward eject move for the proposed state is 1
      }
      z_p <- z # initialize proposed clustering
      z_p[z_p==h] <- g # assign all nodes from cluster h to cluster g
      z_p[z_p>h] <- z_p[z_p>h]-1 # relabel the cluster labels which are greater than h
      n_g_p <- n_g_c+n_h_c # table(z_p)[[g]] # extract n_g for proposed clustering after absorb move
      U_g_p <- matrix(U[z_p==g,],ncol=d) # extract all the u_j for proposed clustering z_j=g after absorb move
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_p]-log(0.01)))] # Recall: check a ranging from 0.01 to 100 and check n from 0 to N
      log_alpha_A_right <-  -(alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega))-log(p_absorb_p)+log(p_eject_p)+ #constant terms
        lgamma(alpha1+n_g_p*d/2) - n_g_p*d/2*log(pi) - d/2*log(omega+n_g_p) - (alpha1+n_g_p*d/2)*log(sum(U_g_p^2)-sum(colSums(U_g_p)^2)/(n_g_p+omega)+alpha2)- # proposed lambda'_g term
        (lgamma(alpha1+n_h_c*d/2) - n_h_c*d/2*log(pi) - d/2*log(omega+n_h_c) - (alpha1+n_h_c*d/2)*log(sum(U_h_c^2)-sum(colSums(U_h_c)^2)/(n_h_c+omega)+alpha2)+ # current lambda^{(t+1)}_h term
           lgamma(alpha1+n_g_c*d/2) - n_g_c*d/2*log(pi) - d/2*log(omega+n_g_c) - (alpha1+n_g_c*d/2)*log(sum(U_g_c^2)-sum(colSums(U_g_c)^2)/(n_g_c+omega)+alpha2))+ # current lambda^{(t+1)}_g term
        log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))-
        log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))-sum(log(alpha+0:(n_h_c-1)))+
        lgamma(2*a)-2*lgamma(a)-lgamma(2*a+n_g_p)+lgamma(a+n_g_c)+lgamma(a+n_h_c)-
        log(1-2*exp(log_p_0(a,n_g_p))) # normalizing term of truncating at the empty cluster cases
      
      if (!is.null(A)){ # if A is provided
        table.z_c.A <- table(z,A)
        table.z_p.A <- table(z_p,A)
        log_alpha_A_right <- log_alpha_A_right +
          sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)-lgamma(sum_omega_c)+C*lgamma(omega_c)-
          sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c)-
          sum(lgamma(table.z_c.A[h,]+omega_c))+lgamma(n_h_c+sum_omega_c) # the ratio of A|z term
      }
      
      # Evaluate the ratio of the nu terms next
      n_k_p <- table(z_p) # n_k of z_p
      Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
      nu_gh_p <- t(Z_p)%*%nu%*%Z_p # {nu_gh} of z_p
      if (length(n_k_p)>1){ # if absorbing leading to K_p>1
        n_gh_p <- outer(n_k_p,n_k_p)-diag(n_k_p) # {n_gh} of z_p
      }else{ # if absorbing leading to K_p=1
        n_gh_p <- matrix(n_k_p^2-n_k_p)
      }
      changed_gh_position <- matrix(FALSE,K,K) # extract those changed blocks
      changed_gh_position[c(g,h),] <- TRUE; changed_gh_position[,c(g,h)] <- TRUE
      changed_gh_position_p <- matrix(FALSE,K_p,K_p) # extract those changed blocks
      changed_gh_position_p[g,] <- TRUE; changed_gh_position_p[,g] <- TRUE
      
      changed_nu_gh <- nu_gh[changed_gh_position]; changed_n_gh <- n_gh[changed_gh_position]
      changed_nu_gh_p <- nu_gh_p[changed_gh_position_p]; changed_n_gh_p <- n_gh_p[changed_gh_position_p]
      
      log_alpha_A_right <- log_alpha_A_right + sum(lbeta(changed_nu_gh_p+beta1,changed_n_gh_p-changed_nu_gh_p+beta2))-
        sum(lbeta(changed_nu_gh+beta1,changed_n_gh-changed_nu_gh+beta2))-lbeta(beta1,beta2)*(sum(changed_gh_position_p)-sum(changed_gh_position))
      
      if (log(runif(1)) <= min(0,log_alpha_A_right)){ # accept or not
        z_list[t+1,] <- z <- z_p
        Z <- Z_p
        K_list[t+1] <- K <- K_p
        acceptance_count_TAE[t] <- 1
        nu_gh <- nu_gh_p
        n_gh <- n_gh_p
      }# else: stay at the current state and do nothing
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update P by Gibbs
    P_list[[t+1]] <- P <- matrix(rbeta(K^2,nu_gh+beta1,n_gh-nu_gh+beta2),K,K)
    #--------------------------------------------------------------------------------------------------------------------------------------------
  } # end T
  return(list(z = z_list,U = U_list,nu = nu_list,K = K_list,X = X_list,P = P_list,beta = beta_list,
              acceptance_count_U = acceptance_count_U, acceptance_count_beta = acceptance_count_beta, acceptance_count_TAE = acceptance_count_TAE))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# The partially collapsed Metropolis-within-Gibbs algorithm for directed Poisson-LPCM-MFM 
MwG_Directed_PoissonLPCM <- function(Y,T, omega,alpha1,alpha2,alpha, sigma2prop_beta,sigma2prop_U,d=3,z=NULL,p_eject=0.5,A=NULL,omega_c=NULL){
  library(mvtnorm) # for rmvnorm
  library(ergmito) # for geodesic
  library(Rfast) # for Rfast::Dist()
  
  N <- nrow(Y) # Extract the number of nodes
  log_p_0 <- function(a,n_k){lgamma(2*a)-lgamma(a)+lgamma(a+n_k)-lgamma(2*a+n_k)} # build the log look-up table for choosing Beta(a,a) in the truncated AE step
  Set_a_for_P_E <- Vectorize(log_p_0,c("a","n_k")) # vectorize the variables a and n_k
  lookuptable_for_a <- outer(seq(0.01,100,0.01),1:N,Set_a_for_P_E) # evaluate for each pair of a and n_k where a=seq(0.01,100,0.01) and n_k=0:N
  # lookuptable_for_a[is.nan(lookuptable_for_a)]<-0 # transform NaN to 0 # no NaN in log form
  # Note that we set logP(n_g=0)=logP(n_{K+1}=0)=log(0.01) in this algorithm
  
  # Y: The N X N observed adjacency matrix
  # A: The 1 X N categorical exogenous node attributes which correspond to \boldsymbol{c} in the paper
  # omega_c: The parameter of DM cohesion for A where omega_c is assumed to be the same for all c = 1,2,...,C
  if (!is.null(A)){ # if A is provided
    C <- max(A)
    if (!is.null(omega_c)){
      sum_omega_c <- omega_c*C # Note that no empty level is allowed in A for c = 1,2,...,C
    }else{
      stop("omega_c should be provided if node attributes A is included!")
    }
  }
  
  # Recall here that
  # beta: a LPM intercept parameter
  # U: N X d matrix, latent positions
  # z: 1 X N vector, memberships
  
  # mu: 2 X K matrix MVN mean (collapsed, may be inferred if required)
  # tau: 1 X K vector MVN precision, tau = 1/sigma2 (collapsed, may be inferred if required)
  # Pi: 1 X K vector membership probability (collapsed, may be inferred if required)
  
  # omega, (alpha1,alpha2): Conjugate prior parameters for mu_k and tau_k, respectively
  # alpha: prior parameter for Pi; 
  # sigma2prop_U,sigma2prop_beta: the proposal variance proposed, respectively, for U and beta
  # d: the dimension of the latent positions
  # z: initial clustering specified
  
  # Note that K is the number of non-empty groups and no empty groups is allowed in z within this function
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize clustering
  if (is.null(z)){ # if initial clustering is not provided
    z <- 1:N # set the clustering where each node is in a singleton group as the initial clustering
  }
  K <- max(z)
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize beta and P
  beta <- 0
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize latent positions by Classical (Metric) Multidimensional Scaling (MDS) on asymmetric distance matrix
  Y_GM <- geodesic(Y)[[1]] # geodesic matrix (GM) of Y, the length/the number of edges of the shortest paths between each pair of nodes
  Y_GM[is.na(Y_GM)] <- 0 # remove NA as 0
  U <- cmdscale(Y_GM,d) # obtain the initial latent positions by MDS
  # Note here that MDS may bring the same position for multiple nodes
  Dist_U <- Rfast::Dist(U)
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  X <- Y
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Posterior chains
  z_list <- matrix(z,nrow=1) # the posterior chain will be a (T+1) X N matrix whose (t+1)th row is the t'th clustering z^{(t)}
  U_list <- list(U) # the posterior chain will be a list()
  beta_list <- c(beta) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th beta^{(t)}
  K_list <- c(K) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th K^{(t)}
  # Note that K here is the number of non-empty clusters
  acceptance_count_U <- matrix(0,T,N) # Count the number of acceptance for the M-H step of each latent position u_i
  acceptance_count_beta <- matrix(0,1,T) # Count the number of acceptance for the M-H step of beta
  acceptance_count_TAE <- matrix(0,1,T) # Count the number of acceptance for the M-H step of AE move
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Main algorithm
  for (t in 1:T){
    if ((t%%100) == 0){
      cat("t=",t,", K=",K,"\n") # monitor the process
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update beta by M-H step
    beta_c <- beta
    beta <- rnorm(1,beta_c,sd=sqrt(sigma2prop_beta)) # Note here that sigma2prop_beta is the proposal variance
    log_alpha_beta_right <- dpois(X,exp(beta-Dist_U),log=TRUE)-dpois(X,exp(beta_c-Dist_U),log=TRUE)
    diag(log_alpha_beta_right) <- 0
    log_alpha_beta_right <- sum(log_alpha_beta_right)
    
    if (log(runif(1)) <= min(0,log_alpha_beta_right)){ # accept or not
      beta_list[t+1] <- beta
      acceptance_count_beta[t] <- 1
    }else{
      beta_list[t+1] <- beta <- beta_c
    }
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each u_i by M-H step
    n_k <- table(z) # n_k of z
    position_matrix <- Z==1 # return TRUE if z_ik = 1 and FALSE otherwise in Z matrix
    for (i in sample.int(N,N)){ # update in random order
      U_c <- U # Update the current latent positions
      Dist_U_c <- Dist_U # Update the Dist_U_c based on U_c
      n_zi <- n_k[[z[i]]] # extract the current number of nodes in group z_i; 
      position_zi <- position_matrix[,z[i]] # z==z[i] # extract the positions within 1:N of those nodes from group z_i 
      U_c_zi <- matrix(U_c[position_zi,],ncol=d) # extract the current latent positions of nodes in group z_i
      # Note that the matrix(,ncol=d) transformation here aims to deal with the single vector case which would not be in the form of a 1 X 2 matrix
      U[i,] <- mvtnorm::rmvnorm(1,U_c[i,],sigma2prop_U*diag(d)) # proposal distribution # Rfast also has the Rfast::rmvnorm which is faster, but it will reset the RNG seed by set.seed(NULL)
      Dist_U <- Rfast::Dist(U) # Update the Dist_U based on proposed U
      U_p_zi <- matrix(U[position_zi,],ncol=d) # extract the proposed latent positions of nodes in group z_i
      log_alpha_U_right <-  (n_zi*d/2 + alpha1)*(log(sum(U_c_zi^2) - sum(colSums(U_c_zi)^2)/(n_zi + omega) + alpha2)-
                                                   log(sum(U_p_zi^2) - sum(colSums(U_p_zi)^2)/(n_zi + omega) + alpha2))+ # f(U|z) terms
        sum(dpois(X[i,-i],exp(beta-Dist_U[i,-i]),log=TRUE)+dpois(X[-i,i],exp(beta-Dist_U[-i,i]),log=TRUE)-
              dpois(X[i,-i],exp(beta-Dist_U_c[i,-i]),log=TRUE)-dpois(X[-i,i],exp(beta-Dist_U_c[-i,i]),log=TRUE)) # f(X|beta,U) terms
      
      if (log(runif(1)) <= min(0,log_alpha_U_right)){ # accept or not
        acceptance_count_U[t,i] <- 1
      }else{
        U <- U_c
        Dist_U <- Dist_U_c
      }# end if
    }# end i
    U_list[[t+1]] <- U
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each z_i via MFM
    # Note that this MFM inference step has no empty clusters
    
    for (i in sample.int(N,N)){
      nonempty <- which(colSums(matrix(Z[-i,],nrow=N-1)) > 0) # extract non-empty groups after removing node i
      Z <- Z[,nonempty] # remove the group if removing node i makes it empty
      if (length(nonempty)==1){Z <- matrix(Z,N,1)} # ensure the matrix form of Z when K = 1 after removing the possible empty group
      K <- ncol(Z) # update K after removing node i
      Z_Ri <- Z[-i,] # remove node i from Z
      if (K==1){Z_Ri <- matrix(Z_Ri,nrow=N-1,ncol=1)} # Z_Ri becomes vector form if k=1
      z_Ri <- c(Z_Ri%*%1:K) # z after removing node i
      U_Ri <- U[-i,] # U after removing node i
      n_k_Ri_full <- colSums(Z_Ri) # an 1 X K matrix of {n_k} after removing node i
      
      z_i <- which(Z[i,] > 0) # this is z_i, e.g., = k; or = integer(0) with length==0 if removing node i leading to an empty cluster
      
      # 1. First obtain f(U|z) terms of the full-conditional distribution by assuming that node i is removed and then is reassigned
      # Construct a N-1 X d*K matrix of latent positions where 
      # the rows which correspond to cluster 1 individuals in the 1st three columns store the corresponding latent positions, and other elements in the 1st three columns are zeros
      # The rows which correspond to cluster 2 individuals in the 2nd three columns store the corresponding latent positions, and other elements in the 2nd three columns are zeros
      # Similarly for all the rest rows and columns
      U_special_matrix_Ri <- matrix(U_Ri,nrow=N-1,ncol=d*K)*matrix(do.call("rbind", replicate(d,Z_Ri,simplify = FALSE)),nrow=N-1,ncol=d*K)
      # Use rbind() function to combine the above matrix with a new K X d*K matrix, where matrix[1,1:3], matrix[2,3:6],matrix[3,7:9],... of this new matrix are U[i,] if d=3 while all other elements are zeros
      U_special_matrix_Ai <- rbind(U_special_matrix_Ri,t(t(matrix(do.call("rbind", replicate(d,diag(K),simplify = FALSE)),nrow=K,ncol=d*K))*U[i,]))
      Z_Ri_spcecial_Ai <- rbind(Z_Ri,diag(K)) # create a special clustering based on Z_Ri, where each cluster is added by one more node whose latent position is U[i,]
      n_k_Ri_full_Ai <- n_k_Ri_full+1 # an 1 X K matrix of {n_k} after assuming that the node i is added to each cluster
      
      # Ri means removing node i; evaluate each non-empty k'th component of the full-conditional function U\{u_i} | z\{z_i}
      # Ai means adding node i; evaluate each non-empty k'th component of the full-conditional function U|z by assuming that node i is assigned in each non-empty cluster
      Uk_terms_Ai <- lgamma(alpha1+n_k_Ri_full_Ai*d/2)-n_k_Ri_full_Ai*d/2*log(pi)-d/2*log(omega+n_k_Ri_full_Ai)-
        (n_k_Ri_full_Ai*d/2+alpha1)*log(rowSums(t(Z_Ri_spcecial_Ai)%*%(U_special_matrix_Ai^2))-
                                          rowSums((t(Z_Ri_spcecial_Ai)%*%U_special_matrix_Ai)^2)/(n_k_Ri_full_Ai+omega)+alpha2)
      Uk_terms_Ri <- lgamma(alpha1+n_k_Ri_full*d/2)-n_k_Ri_full*d/2*log(pi)-d/2*log(omega+n_k_Ri_full)-
        (n_k_Ri_full*d/2+alpha1)*log(rowSums(t(Z_Ri)%*%(U_special_matrix_Ri^2))-
                                       rowSums((t(Z_Ri)%*%U_special_matrix_Ri)^2)/(n_k_Ri_full+omega)+alpha2)
      log_update_prob_U <- c(Uk_terms_Ai-Uk_terms_Ri,lgamma(alpha1+d/2)-d/2*log(pi)-d/2*log(omega+1) + alpha1*log(alpha2)+d/2*log(omega)-lgamma(alpha1)-
                               (d/2+alpha1)*log(sum(U[i,]^2)*(1-1/(1+omega))+alpha2)) # here we evaluate the f(U|z\{z_i},z_i)/f(U\{u_i}|z\{z_i}) as well as the case which assigns node i to the new group K+1
      # If we let U_k be the set of latent positions of nodes from cluster k, then some of the terms above correspond to:
      # rowSums(t(Z_Ri)%*%(U_special_matrix_Ri^2)) is {sum(U_k^2)} before adding node i to each cluster
      # rowSums((t(Z_Ri)%*%U_special_matrix_Ri)^2) is {sum(colSums(U_k)^2)} before adding node i to each cluster
      # rowSums(t(Z_Ri_spcecial_Ai)%*%(U_special_matrix_Ai^2)) is {sum(U_k^2)} after adding node i to each cluster
      # rowSums((t(Z_Ri_spcecial_Ai)%*%U_special_matrix_Ai)^2) is {sum(colSums(U_k)^2)} after adding node i to each cluster
      
      # 3. Obtain the f(z) term
      log_update_prob_z <- c(log(n_k_Ri_full+alpha),log(alpha)+log(MFM_WNK_TruncPois1(N,K+1,alpha))-
                               log(MFM_WNK_TruncPois1(N,K,alpha)))
      
      # 4. Obtain the assignment probability
      log_update_prob <- log_update_prob_U+log_update_prob_z
      if (!is.null(A)){ # if A is provided
        A_Ri <- A[-i] # remove c_i from node attributes A
        table.z_Ri.A_Ri <- table(z_Ri,factor(A_Ri,levels=1:C)) # it's possible that an empty level might exist in A_Ri, so we use factor(,levels) to avoid such a problem
        log_update_prob <- log_update_prob + c(log(table.z_Ri.A_Ri[,A[i]]+omega_c)-log(n_k_Ri_full+sum_omega_c),log(omega_c)-log(sum_omega_c))
      }
      
      update_prob <- exp(log_update_prob - max(log_update_prob))
      
      # 5. Determine the reassignment and update Z,K,nu_gh,n_gh,n_k
      z_i_new <- which(rmultinom(1,1,update_prob) > 0)
      if(length(z_i) == 1){Z[i,z_i]<-0} # if removing node i doesn't bring an empty group, let Z[i,] be zeros; otherwise, Z[i,] are already zeros
      if (z_i_new==K+1){
        Z <- cbind(Z,0)
        K <- K+1
        Z[i,z_i_new] <- 1
      }else{
        Z[i,z_i_new] <- 1
      }
    } # end for i
    
    z <- c(Z%*%c(1:ncol(Z)))
    z_list <- rbind(z_list,z)
    K_list[t+1] <- K
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # TAE step
    if (K == 1){ # if only one cluster, eject
      AE_decision <- 0
      p_eject_p <- 1 # proposal probability of the eject move
      p_absorb_p <- 1-p_eject # proposal probability of the reverse absorb move
    }else if (K == N){ # if N clusters, absorb
      AE_decision <- 1
      p_eject_p <- p_eject
      p_absorb_p <- 1
    }else{ # with probability p_eject to eject or absorb
      AE_decision <- runif(1)
      p_eject_p <- p_eject # P(eject) in the proposal
      p_absorb_p <- 1-p_eject
    }
    
    if (AE_decision<=p_eject){ # eject move
      g <- sample.int(K,1) # randomly choose one cluster
      n_k <- table(z) # n_k of z
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      U_g_c <- matrix(U[z==g,],ncol=d) # extract all the u_j for current clustering z_j=g
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_c]-log(0.01)))] # Recall: check "a" ranging from 0.01 to 100 and check "n" from 1 to N; 0.01 here is set as the Pr(n_{K+1}=0);
      p_E <- rbeta(1,a,a) # sample the probability of the nodes being assigned to cluster K+1
      K_p <- K+1 # proposed number of clusters
      z_p <- z # initialize proposed clustering
      z_p[z_p==g]<-rbinom(n_g_c,1,p_E)*(K_p-g)+g # apply the eject move; for each node j in cluster g, if rbinom() returns 1, assign (K_p-g)+g to z_j; if rbinom() returns 0, assign g back to z_j
      n_k_p <- table(factor(z_p,levels=1:K_p)) # {n_k} of z_p
      n_g_p <- n_k_p[[g]] # extract n_g for proposed clustering after eject move
      # Note that it's possible to have empty cluster in the proposed z and table() can not directly return 0 for empty clusters;
      # Thus we use table(factor(z_p,levels=1:K_p)) to let the empty cluster appear.
      if (n_g_p!=0 & n_g_p!=n_g_c){ # if any of the ejected cluster is not empty
        if (K_p==N){ # If the proposed state has K == N
          p_absorb_p <- 1 # the probability of the backward absorb move for the proposed state is 1
        }
        U_g_p <- matrix(U[z_p==g,],ncol=d) # extract all the u_j for proposed clustering z_j=g after eject move 
        n_Kp_p <- n_k_p[[K_p]] # extract n_Kp for proposed clustering after eject move
        U_Kp_p <- matrix(U[z_p==K_p,],ncol=d) # extract all the u_j for proposed clustering z_j=K_p after eject move 
        log_alpha_E_right <-  alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega)+log(p_absorb_p)-log(p_eject_p)+ # constant terms
          lgamma(alpha1+n_g_p*d/2) - n_g_p*d/2*log(pi) - d/2*log(omega+n_g_p) - (alpha1+n_g_p*d/2)*log(sum(U_g_p^2)-sum(colSums(U_g_p)^2)/(n_g_p+omega)+alpha2)+ # proposed lambda'_g term
          lgamma(alpha1+n_Kp_p*d/2) - n_Kp_p*d/2*log(pi) - d/2*log(omega+n_Kp_p) - (alpha1+n_Kp_p*d/2)*log(sum(U_Kp_p^2)-sum(colSums(U_Kp_p)^2)/(n_Kp_p+omega)+alpha2)- # proposed lambda'_Kp term
          (lgamma(alpha1+n_g_c*d/2) - n_g_c*d/2*log(pi) - d/2*log(omega+n_g_c) - (alpha1+n_g_c*d/2)*log(sum(U_g_c^2)-sum(colSums(U_g_c)^2)/(n_g_c+omega)+alpha2))+ # current lambda^{(t+1)]}_g term; note that z^{(t+1)} was updated from z^{(t)} in the previous step, so current state is z^{(t+1)}
          log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))+sum(log(alpha+0:(n_Kp_p-1)))- # z' MFM prior term
          log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))+ # z MFM prior term
          2*lgamma(a)-lgamma(2*a)+lgamma(2*a+n_g_c)-lgamma(a+n_g_p)-lgamma(a+n_Kp_p)+ # proposal terms
          log(1-2*exp(log_p_0(a,n_g_c))) # normalizing term of truncating at the empty cluster cases
        
        if (!is.null(A)){ # if A is provided
          table.z_c.A <- table(z,A)
          table.z_p.A <- table(z_p,A)
          log_alpha_E_right <- log_alpha_E_right +
            sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)+
            sum(lgamma(table.z_p.A[K_p,]+omega_c))-lgamma(n_Kp_p+sum_omega_c)+lgamma(sum_omega_c)-C*lgamma(omega_c)-
            sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c) # the ratio of A|z term
        }
        
        Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
        
        if (log(runif(1)) <= min(0,log_alpha_E_right)){ # accept or not
          z_list[t+1,] <- z <- z_p # update the new clustering
          Z <- Z_p
          K_list[t+1] <- K <- K_p
          acceptance_count_TAE[t] <- 1 # 1 means accept, 0 means not accept
        }# else: stay at the current state and do nothing
      }# else: stay at the current state and do nothing
    }else{ # absorb step
      random_two_picks <- sort(sample.int(K,2))
      g <- random_two_picks[1] # g < h, we let g absorb h
      h <- random_two_picks[2]
      n_k <- table(z)
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      n_h_c <- n_k[[h]] # extract n_h for current clustering
      U_g_c <- matrix(U[z==g,],ncol=d) # extract all the u_j for current clustering z_j=g
      U_h_c <- matrix(U[z==h,],ncol=d) # extract all the u_j for current clustering z_j=h
      K_p <- K-1
      if (K_p==1){ # If the proposed state has K == 1
        p_eject_p <- 1 # the probability of the backward eject move for the proposed state is 1
      }
      z_p <- z # initialize proposed clustering
      z_p[z_p==h] <- g # assign all nodes from cluster h to cluster g
      z_p[z_p>h] <- z_p[z_p>h]-1 # relabel the cluster labels which are greater than h
      n_g_p <- n_g_c+n_h_c # table(z_p)[[g]] # extract n_g for proposed clustering after absorb move
      U_g_p <- matrix(U[z_p==g,],ncol=d) # extract all the u_j for proposed clustering z_j=g after absorb move 
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_p]-log(0.01)))] # Recall: check a ranging from 0.01 to 100 and check n from 0 to N 
      log_alpha_A_right <-  -(alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega))-log(p_absorb_p)+log(p_eject_p)+ #constant terms
        lgamma(alpha1+n_g_p*d/2) - n_g_p*d/2*log(pi) - d/2*log(omega+n_g_p) - (alpha1+n_g_p*d/2)*log(sum(U_g_p^2)-sum(colSums(U_g_p)^2)/(n_g_p+omega)+alpha2)- # proposed lambda'_g term
        (lgamma(alpha1+n_h_c*d/2) - n_h_c*d/2*log(pi) - d/2*log(omega+n_h_c) - (alpha1+n_h_c*d/2)*log(sum(U_h_c^2)-sum(colSums(U_h_c)^2)/(n_h_c+omega)+alpha2)+ # current lambda^{(t+1)}_h term
           lgamma(alpha1+n_g_c*d/2) - n_g_c*d/2*log(pi) - d/2*log(omega+n_g_c) - (alpha1+n_g_c*d/2)*log(sum(U_g_c^2)-sum(colSums(U_g_c)^2)/(n_g_c+omega)+alpha2))+ # current lambda^{(t+1)}_g term
        log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))-
        log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))-sum(log(alpha+0:(n_h_c-1)))+
        lgamma(2*a)-2*lgamma(a)-lgamma(2*a+n_g_p)+lgamma(a+n_g_c)+lgamma(a+n_h_c)-
        log(1-2*exp(log_p_0(a,n_g_p))) # normalizing term of truncating at the empty cluster cases
      
      if (!is.null(A)){ # if A is provided
        table.z_c.A <- table(z,A)
        table.z_p.A <- table(z_p,A)
        log_alpha_A_right <- log_alpha_A_right +
          sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)-lgamma(sum_omega_c)+C*lgamma(omega_c)-
          sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c)-
          sum(lgamma(table.z_c.A[h,]+omega_c))+lgamma(n_h_c+sum_omega_c) # the ratio of A|z term
      }
      
      # Evaluate the ratio of the nu terms next
      Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
      
      if (log(runif(1)) <= min(0,log_alpha_A_right)){ # accept or not
        z_list[t+1,] <- z <- z_p
        Z <- Z_p
        K_list[t+1] <- K <- K_p
        acceptance_count_TAE[t] <- 1
      }# else: stay at the current state and do nothing
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
  } # end T
  return(list(z = z_list,U = U_list,K = K_list,beta = beta_list,
              acceptance_count_U = acceptance_count_U, acceptance_count_beta = acceptance_count_beta, acceptance_count_TAE = acceptance_count_TAE))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# The partially collapsed Gibbs algorithm for directed ZIP-SBM-MFM 
Gibbs_Directed_ZIPSBM <- function(Y,T, alpha1,alpha2,alpha, beta1,beta2,z=NULL,p_eject=0.5,A=NULL,omega_c=NULL){
  
  if (sum(diag(Y))!=0){
    stop("diag(Y) are not zeros! Self-loops are not allowed!")
  }
  N <- nrow(Y) # Extract the number of nodes
  log_p_0 <- function(a,n_k){lgamma(2*a)-lgamma(a)+lgamma(a+n_k)-lgamma(2*a+n_k)} # build the lookup table for choosing Beta(a,a) in the TAE step
  Set_a_for_P_E <- Vectorize(log_p_0,c("a","n_k")) # Vectorize the variables a and n_k
  lookuptable_for_a <- outer(seq(0.01,100,0.01),1:N,Set_a_for_P_E) # evaluate for each pair of a and n_k where a=seq(0.01,100,0.01) and n_k=0:N
  # lookuptable_for_a[is.nan(lookuptable_for_a)]<-0 # transform NaN to 0 # no NaN in log form
  # Note that we set P(n_g=0)=P(n_{K+1}=0)=0.01 in this algorithm
  
  # Y: The N X N observed adjacency matrix
  # A: The 1 X N categorical exogenous node attributes which correspond to \boldsymbol{c} in the paper
  # omega_c: The parameter of DM cohesion for A where omega_c is assumed to be the same for all c = 1,2,...,C
  if (!is.null(A)){ # if A is provided
    C <- max(A)
    if (!is.null(omega_c)){
      sum_omega_c <- omega_c*C # Note that no empty level is allowed in A for c = 1,2,...,C
    }else{
      stop("omega_c should be provided if node attributes A is included!")
    }
  }
  
  # P: a K X K matrix of clustering dependent missing zero probability
  # Lambda: a K X K matrix of clustering dependent Poisson rate
  # nu: N X N matrix, missing zero indicator
  # z: 1 X N vector, memberships
  # Pi: 1 X K vector membership probability (collapsed, may be inferred if required)
  
  # (alpha1,alpha2): Conjugate prior parameters for each lambda_gh: g,h=1,2,...,K
  # alpha: prior parameter for Pi embeded inside the MFM; 
  # (beta1,beta2): prior parameter for each p_gh: g,h=1,2,...,K
  # d: the dimension of the latent positions
  # z: initial clustering specified
  
  # Note that K is the number of non-empty groups and no empty groups is allowed in z within this function
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize clustering
  if (is.null(z)){ # if initial clustering is not provided
    z <- 1:N # set the clustering where each node is in a singleton group as the initial clustering
  }
  K <- max(z)
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize Lambda and P
  Lambda <- matrix(alpha1/alpha2,K,K) # set lambda_gh prior mean for initial state of each lambda_gh: g,h=1,2,...,K
  P <- matrix(beta1/(beta1+beta2),K,K) # set p_gh prior mean for initial state of each p_gh: g,h=1,2,...,K
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize the missing zero indicator matrix nu and the missing data imputed adj X
  nu <- (Y==0)*1
  diag(nu) <- 0
  Imputation_position <- which(nu==1) # extract the position of those y_{ij}=0 except the diagonal elements
  Imputation_n <- length(Imputation_position) # the number of possible missing data imputed elements
  
  Lambda_ij <- Z%*%Lambda%*%t(Z) # find the lambda_{zizj} for each \nu_ij
  Imp_Lambda_ij <- Lambda_ij[Imputation_position] # extract the lambda_ij for those possibly imputed elements
  P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
  Imp_P_ij <- P_ij[Imputation_position] # extract the p_ij for those possibly imputed elements
  nu[Imputation_position] <- rbinom(Imputation_n,1,Imp_P_ij/(Imp_P_ij + (1-Imp_P_ij)*dpois(0,Imp_Lambda_ij))) # Initialize nu
  
  X <- Y
  nu_1_position <- which(nu==1) # extract the position of those updated nu_ij=1
  X[nu_1_position] <- rpois(length(nu_1_position),Lambda_ij[nu_1_position]) # Initialize X
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Posterior chains
  z_list <- matrix(z,nrow=1) # the posterior chain will be a (T+1) X N matrix whose (t+1)th row is the t'th clustering z^{(t)}
  nu_list <- list(nu) # the posterior chain will be a list()
  Lambda_list <- list(Lambda) # the posterior chain will be a list()
  P_list <- list(P) # the posterior chain will be a list()
  K_list <- c(K) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th K^{(t)}
  # Note that K here is the number of non-empty clusters
  X_list <- list(X) # the posterior chain will be a list()
  acceptance_count_TAE <- matrix(0,1,T) # Count the number of acceptance for the M-H step of AE move
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Main algorithm
  for (t in 1:T){
    if ((t%%100) == 0){
      cat("t=",t,", K=",K,"\n") # monitor the process
    }
    # Note that (.)^{(t)} is the current state and (.)^{(t+1)}is the new state within this function
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sample nu
    Lambda_ij <- Z%*%Lambda%*%t(Z) # find the lambda_{zizj} for each \nu_ij
    Imp_Lambda_ij <- Lambda_ij[Imputation_position] # extract the lambda_ij for those possibly imputed elements
    P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
    Imp_P_ij <- P_ij[Imputation_position] # extract the p_ij for those possibly imputed elements, i.e. all the y_ij=0 except diagonal
    nu[Imputation_position] <- rbinom(Imputation_n,1,Imp_P_ij/(Imp_P_ij + (1-Imp_P_ij)*dpois(0,Imp_Lambda_ij)))
    nu_list[[t+1]] <- nu
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sample X
    X <- Y
    nu_1_position <- which(nu==1) # extract the position of those updated nu_ij=1
    X[nu_1_position] <- rpois(length(nu_1_position),Lambda_ij[nu_1_position])
    X_list[[t+1]] <- X
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each z_i via MFM
    # Note that this MFM inference step has no empty clusters
    X_gh <- t(Z)%*%X%*%Z # {X_gh} of z
    nu_gh <- t(Z)%*%nu%*%Z # {nu_gh} of z
    n_k <- table(z) # n_k of z
    if (length(n_k)>1){
      n_gh <- outer(n_k,n_k)-diag(n_k) # {n_gh} of z; n_k already obtained above
    }else{
      n_gh <- matrix(n_k^2-n_k)
    }
    
    for (i in sample.int(N,N)){
      nonempty <- which(colSums(matrix(Z[-i,],nrow=N-1)) > 0) # extract non-empty groups after removing node i
      Z <- Z[,nonempty] # remove the group if removing node i makes it empty
      if (length(nonempty)==1){Z <- matrix(Z,N,1)} # ensure the matrix form of Z when K = 1 after removing the possible empty group
      K <- ncol(Z) # update K after removing node i
      Z_Ri <- Z[-i,] # remove node i from Z
      if (K==1){Z_Ri <- matrix(Z_Ri,nrow=N-1,ncol=1)} # Z_Ri becomes vector form if k=1
      z_Ri <- c(Z_Ri%*%1:K) # z after removing node i
      n_k_Ri_full <- colSums(Z_Ri) # an 1 X K matrix of {n_k} after removing node i
      
      # 1. and 2. First obtain f(X|z) and f(nu|z) terms of the full-conditional distribution by assuming that node i is removed and then is reassigned
      X_gh_Ri <- matrix(X_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      nu_gh_Ri <- matrix(nu_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      n_gh_Ri <- matrix(n_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      X_ik <- crossprod(Z_Ri, X[i,-i]) # {X_ik}
      X_ki <- crossprod(Z_Ri, X[-i,i]) # {X_ki}
      nu_ik <- crossprod(Z_Ri, nu[i,-i]) # {nu_ik}
      nu_ki <- crossprod(Z_Ri, nu[-i,i]) # {nu_ki}
      z_i <- which(Z[i,] > 0) # this is z_i, e.g., = k; or = integer(0) with length==0 if removing node i leading to an empty cluster
      if (length(z_i)==1){ # if removing node i doesn't bring an empty cluster
        X_gh_resid <- matrix(0,K,K) # calculate the difference in nu_gh after removing node i
        X_gh_resid[z_i,] <- X_ik; X_gh_resid[,z_i] <- X_gh_resid[,z_i]+X_ki
        X_gh_Ri <- X_gh_Ri-X_gh_resid
        
        nu_gh_resid <- matrix(0,K,K) # calculate the difference in nu_gh after removing node i
        nu_gh_resid[z_i,] <- nu_ik; nu_gh_resid[,z_i] <- nu_gh_resid[,z_i]+nu_ki
        nu_gh_Ri <- nu_gh_Ri-nu_gh_resid
        
        n_gh_resid <- matrix(0,K,K) # calculate the difference in n_gh after removing node i
        n_gh_resid[z_i,] <- n_k_Ri_full; n_gh_resid[,z_i] <- n_gh_resid[,z_i]+n_k_Ri_full
        n_gh_Ri <- n_gh_Ri-n_gh_resid
      }#else{X_gh_Ri <- X_gh_Ri;nu_gh_Ri <- nu_gh_Ri;n_gh_Ri <- n_gh_Ri} # if removing node i brings an empty cluster which is already removed
      
      if (K>1){
        X_ik_row_Matrix <- matrix(1,K,1)%*%matrix(X_ik,1,K) # construct a K X K matrix where each row is X_ik
        X_ki_column_Matrix <- matrix(X_ki,K,K) # construct a K X K matrix where each column is X_ki
        nu_ik_row_Matrix <- matrix(1,K,1)%*%matrix(nu_ik,1,K) # construct a K X K matrix where each row is nu_ik
        nu_ki_column_Matrix <- matrix(nu_ki,K,K) # construct a K X K matrix where each column is nu_ki
        n_k_Ri_full_row_Matrix <- matrix(1,K,1)%*%matrix(n_k_Ri_full,1,K) # construct a K X K matrix where each row is n_k_Ri_full; thus the transport is a K X K matrix where each column is n_k_Ri_full
        
        changed_X_gh_row <- X_gh_Ri+X_ik_row_Matrix+diag(diag(X_ki_column_Matrix)) # construct a K X K matrix where each row is the row of corresponding nu_gh after adding node i
        changed_X_gh_column <- X_gh_Ri+X_ki_column_Matrix # construct a K X K matrix where each column is the column of corresponding nu_gh after adding node i but excluding diagonal
        diag(changed_X_gh_column) <- NA # excluding diagonal
        changed_X_gh_column <- matrix(changed_X_gh_column[!is.na(changed_X_gh_column)],nrow=K-1,ncol=K) # excluding diagonal leading to a K-1 X K matrix
        
        changed_nu_gh_row <- nu_gh_Ri+nu_ik_row_Matrix+diag(diag(nu_ki_column_Matrix)) # construct a K X K matrix where each row is the row of corresponding nu_gh after adding node i
        changed_nu_gh_column <- nu_gh_Ri+nu_ki_column_Matrix # construct a K X K matrix where each column is the column of corresponding nu_gh after adding node i but excluding diagonal
        diag(changed_nu_gh_column) <- NA # excluding diagonal
        changed_nu_gh_column <- matrix(changed_nu_gh_column[!is.na(changed_nu_gh_column)],nrow=K-1,ncol=K) # excluding diagonal leading to a K-1 X K matrix
        
        changed_n_gh_row <- n_gh_Ri+n_k_Ri_full_row_Matrix+diag(n_k_Ri_full) # construct a K X K matrix where each row is the row of corresponding n_gh after adding node i
        changed_n_gh_column <- n_gh_Ri+t(n_k_Ri_full_row_Matrix) # construct a K X K matrix where each column is the column of corresponding n_gh after adding node i but excluing diagonal
        diag(changed_n_gh_column) <- NA # excluding diagonal
        changed_n_gh_column <- matrix(changed_n_gh_column[!is.na(changed_n_gh_column)],nrow=K-1,ncol=K) # excluding diagonal leading to a K-1 X K matrix
        
        log_update_prob_Lambda <- c(rowSums(lgamma(changed_X_gh_row+alpha1)-(changed_X_gh_row+alpha1)*log(changed_n_gh_row+alpha2))+ # add the row terms including diag for each z_i=k case
                                      colSums(lgamma(changed_X_gh_column+alpha1)-(changed_X_gh_column+alpha1)*log(changed_n_gh_column+alpha2))- # add the column terms excluding diag for each z_i=k case
                                      rowSums(lgamma(X_gh_Ri+alpha1)-(X_gh_Ri+alpha1)*log(n_gh_Ri+alpha2))-colSums(lgamma(X_gh_Ri+alpha1)-(X_gh_Ri+alpha1)*log(n_gh_Ri+alpha2))+ # evaluate the demominator of the f(\nu|z_{-i},z_i=k)/f(\nu_{-i}|z_{-i}) for each k = 1,2,...,K
                                      lgamma(diag(X_gh_Ri)+alpha1)-(diag(X_gh_Ri)+alpha1)*log(diag(n_gh_Ri)+alpha2),# remove the repeat diag term
                                    sum(lgamma(c(X_ik,X_ki)+alpha1)-(c(X_ik,X_ki)+alpha1)*log(c(n_k_Ri_full,n_k_Ri_full)+alpha2)+alpha1*log(alpha2)-lgamma(alpha1))) # add the term relevant to the new group K+1
        
        log_update_prob_nu <- c(rowSums(lbeta(changed_nu_gh_row+beta1,changed_n_gh_row-changed_nu_gh_row+beta2))+ # add the row terms including diag for each z_i=k case
                                  colSums(lbeta(changed_nu_gh_column+beta1,changed_n_gh_column-changed_nu_gh_column+beta2))- # add the column terms excluding diag for each z_i=k case
                                  rowSums(lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2))-colSums(lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2))+ # evaluate the demominator of the f(\nu|z_{-i},z_i=k)/f(\nu_{-i}|z_{-i}) for each k = 1,2,...,K
                                  lbeta(diag(nu_gh_Ri)+beta1,diag(n_gh_Ri-nu_gh_Ri)+beta2),# remove the repeat diag term
                                sum(lbeta(c(nu_ik,nu_ki)+beta1,c(n_k_Ri_full,n_k_Ri_full)-c(nu_ik,nu_ki)+beta2)-lbeta(beta1,beta2))) # add the term relevant to the new group K+1
      }else{ # if removing node i leads to only 1 group
        log_update_prob_Lambda <- c(lgamma(X_gh_Ri+X_ki+X_ik+alpha1)-(X_gh_Ri+X_ki+X_ik+alpha1)*log(n_gh_Ri+2*n_k_Ri_full+alpha2)-
                                      lgamma(X_gh_Ri+alpha1)+(X_gh_Ri+alpha1)*log(n_gh_Ri+alpha2),
                                    sum(lgamma(c(X_ik,X_ki)+alpha1)-(c(X_ik,X_ki)+alpha1)*log(c(n_k_Ri_full,n_k_Ri_full)+alpha2)+alpha1*log(alpha2)-lgamma(alpha1))) # add the term relevant to the new group K+1
        
        log_update_prob_nu <- c(lbeta(nu_gh_Ri+nu_ki+nu_ik+beta1,n_gh_Ri+2*n_k_Ri_full-(nu_gh_Ri+nu_ki+nu_ik)+beta2)-
                                  lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2),
                                sum(lbeta(c(nu_ik,nu_ki)+beta1,c(n_k_Ri_full,n_k_Ri_full)-c(nu_ik,nu_ki)+beta2)-lbeta(beta1,beta2))) # add the term relevant to the new group K+1
      }
      
      # 3. Obtain the f(z) term
      log_update_prob_z <- c(log(n_k_Ri_full+alpha),log(alpha)+log(MFM_WNK_TruncPois1(N,K+1,alpha))-
                               log(MFM_WNK_TruncPois1(N,K,alpha)))
      
      # 4. Obtain the assignment probability
      log_update_prob <- log_update_prob_Lambda+log_update_prob_nu+log_update_prob_z
      if (!is.null(A)){ # if A is provided
        A_Ri <- A[-i] # remove c_i from node attributes A
        table.z_Ri.A_Ri <- table(z_Ri,factor(A_Ri,levels=1:C)) # it's possible that an empty level might exist in A_Ri, so we use factor(,levels) to avoid such a problem
        log_update_prob <- log_update_prob + c(log(table.z_Ri.A_Ri[,A[i]]+omega_c)-log(n_k_Ri_full+sum_omega_c),log(omega_c)-log(sum_omega_c))
      }
      update_prob <- exp(log_update_prob - max(log_update_prob))
      
      # 5. Determine the reassignment and update Z,K,X_gh,nu_gh,n_gh,n_k
      z_i_new <- which(rmultinom(1,1,update_prob) > 0)
      if(length(z_i) == 1){Z[i,z_i]<-0} # if removing node i doesn't bring an empty group, let Z[i,] be zeros; otherwise, Z[i,] are already zeros
      if (z_i_new==K+1){
        Z <- cbind(Z,0)
        K <- K+1
        Z[i,z_i_new] <- 1
        X_gh_Ri <- rbind(cbind(X_gh_Ri,0),0)
        nu_gh_Ri <- rbind(cbind(nu_gh_Ri,0),0)
        n_gh_Ri <- rbind(cbind(n_gh_Ri,0),0)
        X_ik <- c(X_ik,0); X_ki <- c(X_ki,0)
        nu_ik <- c(nu_ik,0); nu_ki <- c(nu_ki,0)
        n_k_Ri_full <- c(n_k_Ri_full,0)
      }else{
        Z[i,z_i_new] <- 1
      }
      
      X_gh_resid <- matrix(0,K,K) # X_gh difference by adding node i to group k
      X_gh_resid[z_i_new,] <- X_ik; X_gh_resid[,z_i_new] <- X_gh_resid[,z_i_new]+X_ki
      X_gh <- X_gh_Ri + X_gh_resid
      
      nu_gh_resid <- matrix(0,K,K) # nu_gh difference by adding node i to group k
      nu_gh_resid[z_i_new,] <- nu_ik; nu_gh_resid[,z_i_new] <- nu_gh_resid[,z_i_new]+nu_ki
      nu_gh <- nu_gh_Ri + nu_gh_resid
      
      n_gh_resid <- matrix(0,K,K) # n_gh difference by adding node i to group k
      n_gh_resid[z_i_new,] <- n_k_Ri_full; n_gh_resid[,z_i_new] <- n_gh_resid[,z_i_new]+n_k_Ri_full
      n_gh <- n_gh_Ri + n_gh_resid
    } # end for i
    
    z <- c(Z%*%c(1:ncol(Z)))
    z_list <- rbind(z_list,z)
    K_list[t+1] <- K
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # TAE step
    if (K == 1){ # if only one cluster, eject
      AE_decision <- 0
      p_eject_p <- 1
      p_absorb_p <- 1-p_eject
    }else if (K == N){ # if N cluster, absorb
      AE_decision <- 1
      p_eject_p <- p_eject
      p_absorb_p <- 1
    }else{ # with probability p_eject to eject or absorb
      AE_decision <- runif(1)
      p_eject_p <- p_eject # P(eject) in the proposal
      p_absorb_p <- 1-p_eject
    }
    if (AE_decision<=p_eject){ # eject move
      g <- sample.int(K,1) # randomly choose one cluster
      n_k <- table(z) # n_k of z
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_c]-log(0.01)))] # Recall: check a ranging from 0.01 to 100 and check n from 1 to N; 0.01 here is set as the Pr(n_{K+1}=0);
      p_E <- rbeta(1,a,a) # sample the probability of the nodes being assigned to cluster K+1
      K_p <- K+1 # proposed number of clusters
      z_p <- z # initialize proposed clustering
      z_p[z_p==g]<-rbinom(n_g_c,1,p_E)*(K_p-g)+g # apply the eject move; for each node j in cluster g, if rbinom() returns 1, assign (K_p-g)+g to z_j; if rbinom() returns 0, assign g back to z_j
      n_k_p <- table(factor(z_p,levels=1:K_p)) # {n_k} of z_p
      n_g_p <- n_k_p[[g]] # extract n_g for proposed clustering after eject move
      # note that it's possible to have empty cluster in the proposed z and table() can not directly return 0 for empty clusters;
      # Thus we use table(factor(z_p,levels=1:K_p)) to let the empty cluster appear.
      if (n_g_p!=0 & n_g_p!=n_g_c){ # if any of the ejected cluster is not empty
        if (K_p==N){ # If the proposed state has K == N
          p_absorb_p <- 1 # the probability of the backward absorb move for the proposed state is 1
        }
        n_Kp_p <- n_k_p[[K_p]] # extract n_Kp for proposed clustering after eject move
        log_alpha_E_right <-  log(p_absorb_p)-log(p_eject_p)+
          log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))+sum(log(alpha+0:(n_Kp_p-1)))- # z' MFM prior term
          log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))+ # z MFM prior term
          2*lgamma(a)-lgamma(2*a)+lgamma(2*a+n_g_c)-lgamma(a+n_g_p)-lgamma(a+n_Kp_p)+ # proposal terms
          log(1-2*exp(log_p_0(a,n_g_c))) # normalizing term of truncating at the empty cluster cases
        
        if (!is.null(A)){ # if A is provided
          table.z_c.A <- table(z,A)
          table.z_p.A <- table(z_p,A)
          log_alpha_E_right <- log_alpha_E_right +
            sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)+
            sum(lgamma(table.z_p.A[K_p,]+omega_c))-lgamma(n_Kp_p+sum_omega_c)+lgamma(sum_omega_c)-C*lgamma(omega_c)-
            sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c) # the ratio of A|z term
        }
        
        Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
        X_gh_p <- t(Z_p)%*%X%*%Z_p # {nu_gh} of z_p
        nu_gh_p <- t(Z_p)%*%nu%*%Z_p # {nu_gh} of z_p
        n_gh_p <- outer(n_k_p,n_k_p)-diag(n_k_p) # {n_gh} of z_p
        
        changed_gh_position <- matrix(FALSE,K,K) # extract those changed blocks
        changed_gh_position[g,] <- TRUE; changed_gh_position[,g] <- TRUE
        changed_gh_position_p <- matrix(FALSE,K_p,K_p) # extract those changed blocks
        changed_gh_position_p[c(g,K_p),] <- TRUE; changed_gh_position_p[,c(g,K_p)] <- TRUE
        
        changed_X_gh <- X_gh[changed_gh_position];changed_nu_gh <- nu_gh[changed_gh_position]; changed_n_gh <- n_gh[changed_gh_position]
        changed_X_gh_p <- X_gh_p[changed_gh_position_p];changed_nu_gh_p <- nu_gh_p[changed_gh_position_p]; changed_n_gh_p <- n_gh_p[changed_gh_position_p]
        
        log_alpha_E_right <- log_alpha_E_right + 
          sum(lgamma(changed_X_gh_p+alpha1)-(changed_X_gh_p+alpha1)*log(changed_n_gh_p+alpha2))-
          sum(lgamma(changed_X_gh+alpha1)-(changed_X_gh+alpha1)*log(changed_n_gh+alpha2))+(alpha1*log(alpha2)-lgamma(alpha1))*(sum(changed_gh_position_p)-sum(changed_gh_position))+
          sum(lbeta(changed_nu_gh_p+beta1,changed_n_gh_p-changed_nu_gh_p+beta2))-
          sum(lbeta(changed_nu_gh+beta1,changed_n_gh-changed_nu_gh+beta2))-lbeta(beta1,beta2)*(sum(changed_gh_position_p)-sum(changed_gh_position))
        
        if (log(runif(1)) <= min(0,log_alpha_E_right)){ # accept or not
          z_list[t+1,] <- z <- z_p # update the new clustering
          Z <- Z_p
          K_list[t+1] <- K <- K_p
          acceptance_count_TAE[t] <- 1 # 1 means accept, 0 means not accept
          X_gh <- X_gh_p
          nu_gh <- nu_gh_p
          n_gh <- n_gh_p
        }# else: stay at the current state and do nothing
      }# else: stay at the current state and do nothing
    }else{ # absorb step
      random_two_picks <- sort(sample.int(K,2))
      g <- random_two_picks[1] # g < h, we let g absorb h
      h <- random_two_picks[2]
      n_k <- table(z)
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      n_h_c <- n_k[[h]] # extract n_h for current clustering
      K_p <- K-1
      if (K_p==1){ # If the proposed state has K == 1
        p_eject_p <- 1 # the probability of the backward eject move for the proposed state is 1
      }
      z_p <- z # initialize proposed clustering
      z_p[z_p==h] <- g # assign all nodes from cluster h to cluster g
      z_p[z_p>h] <- z_p[z_p>h]-1 # relabel the cluster labels which are greater than h
      n_g_p <- n_g_c+n_h_c # table(z_p)[[g]] # extract n_g for proposed clustering after absorb move
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_p]-log(0.01)))] # Recall: check a ranging from 0.01 to 100 and check n from 0 to N
      log_alpha_A_right <-  -log(p_absorb_p)+log(p_eject_p)+ #constant terms
        log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))-
        log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))-sum(log(alpha+0:(n_h_c-1)))+
        lgamma(2*a)-2*lgamma(a)-lgamma(2*a+n_g_p)+lgamma(a+n_g_c)+lgamma(a+n_h_c)-
        log(1-2*exp(log_p_0(a,n_g_p))) # normalizing term of truncating at the empty cluster cases
      
      if (!is.null(A)){ # if A is provided
        table.z_c.A <- table(z,A)
        table.z_p.A <- table(z_p,A)
        log_alpha_A_right <- log_alpha_A_right +
          sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)-lgamma(sum_omega_c)+C*lgamma(omega_c)-
          sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c)-
          sum(lgamma(table.z_c.A[h,]+omega_c))+lgamma(n_h_c+sum_omega_c) # the ratio of A|z term
      }
      
      # Evaluate the ratio of the Lambda and nu terms next
      n_k_p <- table(z_p) # n_k of z_p
      Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
      X_gh_p <- t(Z_p)%*%X%*%Z_p # {nu_gh} of z_p
      nu_gh_p <- t(Z_p)%*%nu%*%Z_p # {nu_gh} of z_p
      if (length(n_k_p)>1){ # if absorbing leading to K_p>1
        n_gh_p <- outer(n_k_p,n_k_p)-diag(n_k_p) # {n_gh} of z_p
      }else{ # if absorbing leading to K_p=1
        n_gh_p <- matrix(n_k_p^2-n_k_p)
      }
      changed_gh_position <- matrix(FALSE,K,K) # extract those changed blocks
      changed_gh_position[c(g,h),] <- TRUE; changed_gh_position[,c(g,h)] <- TRUE
      changed_gh_position_p <- matrix(FALSE,K_p,K_p) # extract those changed blocks
      changed_gh_position_p[g,] <- TRUE; changed_gh_position_p[,g] <- TRUE
      
      changed_X_gh <- X_gh[changed_gh_position];changed_nu_gh <- nu_gh[changed_gh_position]; changed_n_gh <- n_gh[changed_gh_position]
      changed_X_gh_p <- X_gh_p[changed_gh_position_p];changed_nu_gh_p <- nu_gh_p[changed_gh_position_p]; changed_n_gh_p <- n_gh_p[changed_gh_position_p]
      
      log_alpha_A_right <- log_alpha_A_right + 
        sum(lgamma(changed_X_gh_p+alpha1)-(changed_X_gh_p+alpha1)*log(changed_n_gh_p+alpha2))-
        sum(lgamma(changed_X_gh+alpha1)-(changed_X_gh+alpha1)*log(changed_n_gh+alpha2))+(alpha1*log(alpha2)-lgamma(alpha1))*(sum(changed_gh_position_p)-sum(changed_gh_position))+
        sum(lbeta(changed_nu_gh_p+beta1,changed_n_gh_p-changed_nu_gh_p+beta2))-
        sum(lbeta(changed_nu_gh+beta1,changed_n_gh-changed_nu_gh+beta2))-lbeta(beta1,beta2)*(sum(changed_gh_position_p)-sum(changed_gh_position))
      
      if (log(runif(1)) <= min(0,log_alpha_A_right)){ # accept or not
        z_list[t+1,] <- z <- z_p
        Z <- Z_p
        K_list[t+1] <- K <- K_p
        acceptance_count_TAE[t] <- 1
        X_gh <- X_gh_p
        nu_gh <- nu_gh_p
        n_gh <- n_gh_p
      }# else: stay at the current state and do nothing
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update Lambda and P by Gibbs
    Lambda_list[[t+1]] <- Lambda <- matrix(rgamma(K^2,X_gh+alpha1,n_gh+alpha2),K,K)
    P_list[[t+1]] <- P <- matrix(rbeta(K^2,nu_gh+beta1,n_gh-nu_gh+beta2),K,K)
    #--------------------------------------------------------------------------------------------------------------------------------------------
  } # end T
  return(list(z = z_list,Lambda = Lambda_list,nu = nu_list,K = K_list,X = X_list,P = P_list,acceptance_count_TAE = acceptance_count_TAE))
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Evaluating the likelihood of ZIP-LPCM-MFM without collapsing Lambda and P
Directed_ZIPSBM_MFM_CompleteLikelihood <- function(X,Lambda,nu,P,z, alpha, A=NULL,omega_c=NULL){
  # X: N X N matrix, missing data imputed adjacency matrix
  # Lambda: K X K matrix, unusual zero probability for each pair of clusters
  # nu: N X N matrix, unusual zero indicator
  # P: K X K matrix, unusual zero probability for each pair of clusters
  # z: 1 X N vector, clustering without empty clusters
  # alpha: MFM parameter
  # A: 1 X N vector, categorical node attributes, this correspond to the node attributes \boldsymbol{c} in the paper
  # omega_c: node attributes cohesion parameter
  
  K <- max(z) # extract the number of non-empty clusters
  N <- length(z) # extract the network size
  n_k_table <- table(z)
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  res <- 0 # initialize the posterior
  
  # Evaluate the X|Lambda,z term
  Lambda_ij <- Z%*%Lambda%*%t(Z) # find the p_{zizj} for each \nu_ij
  res_temp <- dpois(X,Lambda_ij,log=TRUE)
  diag(res_temp) <- 0
  res <- res + sum(res_temp)
  
  # Evaluate the \nu|P,z term
  P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
  res_temp <- dbinom(nu,1,P_ij,log=TRUE)
  diag(res_temp) <- 0
  res <- res + sum(res_temp)
  
  if (!is.null(A)){
    # Evaluate the A|z term
    C <- max(A)
    table.z.A <- table(z,A)
    sum_omega_c <- omega_c*C # Note that no empty level is allowed in A for c = 1,2,...,C
    res <- res + sum(lgamma(table.z.A+omega_c)) - sum(lgamma(n_k_table+sum_omega_c)) + (lgamma(sum_omega_c)-C*lgamma(omega_c))*K
  }
  
  # Evaluate the z term
  res <- res + log(MFM_WNK_TruncPois1(N,K,alpha))
  for (k in 1:K){
    res <- res + sum(log(alpha+0:(n_k_table[[k]]-1)))
  }
  
  return(res)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

### The above PCMwG ZIP-LPCM algorithm for directed networks can be easily reduced to the undirected case shown below

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# The partially collapsed Metropolis-within-Gibbs algorithm for undirected ZIP-LPCM-MFM 
MwG_UnDirected_ZIPLPCM <- function(Y,T, omega,alpha1,alpha2,alpha, beta1,beta2,sigma2prop_beta,sigma2prop_U,d=3,z=NULL,p_eject=0.5,A=NULL,omega_c=NULL){
  library(mvtnorm) # for rmvnorm
  library(ergmito) # for geodesic
  library(Rfast) # for Rfast::Dist()
  library(gdata) # for upperTriangle() and lowerTriangle()
  
  N <- nrow(Y) # Extract the number of nodes
  log_p_0 <- function(a,n_k){lgamma(2*a)-lgamma(a)+lgamma(a+n_k)-lgamma(2*a+n_k)} # build the look-up table for choosing Beta(a,a) in the TAE step
  Set_a_for_P_E <- Vectorize(log_p_0,c("a","n_k")) # vectorize the variables a and n_k
  lookuptable_for_a <- outer(seq(0.01,100,0.01),1:N,Set_a_for_P_E) # evaluate for each pair of a and n_k where a=seq(0.01,100,0.01) and n_k=0:N
  # lookuptable_for_a[is.nan(lookuptable_for_a)]<-0 # transform NaN to 0 # no NaN in log form
  # Note that we set logP(n_g=0)=logP(n_{K+1}=0)=log(0.01) in this algorithm
  
  # Y: The N X N observed adjacency matrix
  # A: The 1 X N categorical exogenous node attributes which correspond to \boldsymbol{c} in the paper
  # omega_c: The parameter of DM cohesion for A where omega_c is assumed to be the same for all c = 1,2,...,C
  if (!is.null(A)){ # if A is provided
    C <- max(A)
    if (!is.null(omega_c)){
      sum_omega_c <- omega_c*C # Note that no empty level is allowed in A for c = 1,2,...,C
    }else{
      stop("omega_c should be provided if node attributes A is included!")
    }
  }
  
  # Recall here that
  # beta: a LPM intercept parameter
  # P: a K X K matrix of clustering dependent unusual zero probability
  # U: N X d matrix, latent positions
  # nu: N X N matrix, unusual zero indicator
  # z: 1 X N vector, memberships
  
  # mu: d X K matrix MVN mean (collapsed, may be inferred if required)
  # tau: 1 X K vector MVN precision, tau = 1/sigma2 (collapsed, may be inferred if required)
  # Pi: 1 X K vector membership probability (collapsed, may be inferred if required)
  
  # omega, (alpha1,alpha2): Conjugate prior parameters for mu_k and tau_k, respectively
  # alpha: prior parameter for Pi; (beta1,beta2): prior parameter for each p_gh: g,h=1,2,...,K
  # sigma2prop_U,sigma2prop_beta: the proposal variance proposed, respectively, for U and beta
  # d: the dimension of the latent positions
  # z: initial clustering specified
  
  # Note that K is the number of non-empty groups and no empty groups is allowed in z within this function
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize clustering
  if (is.null(z)){ # if initial clustering is not provided
    z <- 1:N # set the clustering where each node is in a singleton group as the initial clustering
  }
  K <- max(z)
  Z <- t(t(matrix(z,N,K))==(1:K))*1 # transfer vector z to matrix Z
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize beta and P
  beta <- 0
  P <- matrix(beta1/(beta1+beta2),K,K) # set p_gh prior mean for initial state of each p_gh: g,h=1,2,...,K
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize latent positions by Classical (Metric) Multidimensional Scaling (MDS) on asymmetric distance matrix
  Y_GM <- geodesic(Y)[[1]] # geodesic matrix (GM) of Y, the length/the number of edges of the shortest paths between each pair of nodes
  Y_GM[is.na(Y_GM)] <- 0 # remove NA as 0
  U <- cmdscale(Y_GM,d) # obtain the initial latent positions by MDS
  # Note here that MDS may bring the same position for multiple nodes
  Dist_U <- Rfast::Dist(U)
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize the missing zero indicator matrix nu and the missing data imputed adj X
  nu <- (Y==0)*1 # note that Y is undirected and thus is symmetric
  diag(nu) <- 0
  Imputation_position <- which(upperTriangle(nu)==1) # extract the position of those y_{ij}=0 for upper triangle elements
  Imputation_n <- length(Imputation_position) # the number of possible missing data imputed elements
  
  P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij # note that P is also symmetric 
  Imp_P_ij <- upperTriangle(P_ij)[Imputation_position] # extract the p_ij for those possibly imputed elements
  upperTriangle(nu)[Imputation_position] <- rbinom(Imputation_n,1,Imp_P_ij/(Imp_P_ij + (1-Imp_P_ij)*dpois(0,exp(beta-upperTriangle(Dist_U)[Imputation_position])))) # Initialize nu
  lowerTriangle(nu,byrow=TRUE)<-upperTriangle(nu)
  
  X <- Y
  nu_1_position <- which(upperTriangle(nu)==1) # extract the position of those updated nu_ij=1
  upperTriangle(X)[nu_1_position] <- rpois(length(nu_1_position),exp(beta-upperTriangle(Dist_U)[nu_1_position])) # Initialize X
  lowerTriangle(X,byrow=TRUE)<-upperTriangle(X)
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Posterior chains
  z_list <- matrix(z,nrow=1) # the posterior chain will be a (T+1) X N matrix whose (t+1)th row is the t'th clustering z^{(t)}
  U_list <- list(U) # the posterior chain will be a list()
  nu_list <- list(nu) # the posterior chain will be a list()
  beta_list <- c(beta) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th beta^{(t)}
  P_list <- list(P) # the posterior chain will be a list()
  K_list <- c(K) # the posterior chain will be a 1 X (T+1) matrix whose (t+1)'th entry is the t'th K^{(t)}
  # Note that K here is the number of non-empty clusters
  X_list <- list(X) # the posterior chain will be a list()
  acceptance_count_U <- matrix(0,T,N) # Count the number of acceptance for the M-H step of each latent position u_i
  acceptance_count_beta <- matrix(0,1,T) # Count the number of acceptance for the M-H step of beta
  acceptance_count_TAE <- matrix(0,1,T) # Count the number of acceptance for the M-H step of AE move
  
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Main algorithm
  for (t in 1:T){
    if ((t%%100) == 0){
      cat("t=",t,", K=",K,"\n") # monitor the process
    }
    # Note that (.)^{(t)} is the current state and (.)^{(t+1)}is the new state within this function
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sample nu
    P_ij <- Z%*%P%*%t(Z) # find the p_{zizj} for each \nu_ij
    Imp_P_ij <- upperTriangle(P_ij)[Imputation_position] # extract the p_ij for those possibly imputed elements, i.e. all the y_ij=0 except diagonal
    upperTriangle(nu)[Imputation_position] <- rbinom(Imputation_n,1,Imp_P_ij/(Imp_P_ij + (1-Imp_P_ij)*dpois(0,exp(beta-upperTriangle(Dist_U)[Imputation_position]))))
    lowerTriangle(nu,byrow=TRUE)<-upperTriangle(nu)
    nu_list[[t+1]] <- nu
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sample X
    X <- Y
    nu_1_position <- which(upperTriangle(nu)==1) # extract the position of those updated nu_ij=1
    upperTriangle(X)[nu_1_position] <- rpois(length(nu_1_position),exp(beta-upperTriangle(Dist_U)[nu_1_position]))
    lowerTriangle(X,byrow=TRUE)<-upperTriangle(X)
    X_list[[t+1]] <- X
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update beta by M-H step
    beta_c <- beta
    beta <- rnorm(1,beta_c,sd=sqrt(sigma2prop_beta)) # Note here that sigma2prop_beta is the proposal variance
    log_alpha_beta_right <- sum(dpois(upperTriangle(X),exp(beta-upperTriangle(Dist_U)),log=TRUE)-dpois(upperTriangle(X),exp(beta_c-upperTriangle(Dist_U)),log=TRUE))
    
    if (log(runif(1)) <= min(0,log_alpha_beta_right)){ # accept or not
      beta_list[t+1] <- beta
      acceptance_count_beta[t] <- 1
    }else{
      beta_list[t+1] <- beta <- beta_c
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each u_i by M-H step
    n_k <- Rfast::Table(z) # n_k of z
    position_matrix <- Z==1 # return TRUE if z_ik = 1 and FALSE otherwise in Z matrix
    for (i in sample.int(N,N)){ # update in random order
      U_c <- U # Update the current latent positions
      Dist_U_c <- Dist_U # Update the Dist_U_c based on U_c
      n_zi <- n_k[[z[i]]] # extract the current number of nodes in group z_i; 
      position_zi <- position_matrix[,z[i]] # z==z[i] # extract the positions within 1:N of those nodes from group z_i 
      U_c_zi <- matrix(U_c[position_zi,],ncol=d) # extract the current latent positions of nodes in group z_i
      # Note that the matrix(,ncol=d) transformation here aims to deal with the single vector case which would not be in the form of a 1 X 2 matrix
      U[i,] <- mvtnorm::rmvnorm(1,U_c[i,],sigma2prop_U*diag(d)) # proposal distribution # Rfast also has the Rfast::rmvnorm which is much faster, but it will set.seed(NULL)
      Dist_U <- Rfast::Dist(U) # Update the Dist_U based on proposed U
      U_p_zi <- matrix(U[position_zi,],ncol=d) # extract the proposed latent positions of nodes in group z_i
      log_alpha_U_right <-  (n_zi*d/2 + alpha1)*(log(sum(U_c_zi^2) - sum(colSums(U_c_zi)^2)/(n_zi + omega) + alpha2)-
                                                   log(sum(U_p_zi^2) - sum(colSums(U_p_zi)^2)/(n_zi + omega) + alpha2))+ # f(U|z) terms
        sum(dpois(X[i,-i],exp(beta-Dist_U[i,-i]),log=TRUE)-dpois(X[i,-i],exp(beta-Dist_U_c[i,-i]),log=TRUE)) # f(X|beta,U) terms
      
      if (log(runif(1)) <= min(0,log_alpha_U_right)){ # accept or not
        acceptance_count_U[t,i] <- 1
      }else{
        U <- U_c
        Dist_U <- Dist_U_c
      }# end if
    }# end i
    U_list[[t+1]] <- U
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update each z_i via MFM
    # Note that this MFM inference step has no empty clusters
    nu_gh <- t(Z)%*%nu%*%Z # {nu_gh} of z
    diag(nu_gh) <- diag(nu_gh)/2 # diag/2
    if (length(n_k)>1){
      n_gh <- outer(n_k,n_k)-diag(n_k) # {n_gh} of z; n_k already obtained above
      diag(n_gh) <- diag(n_gh)/2 # diag/2
    }else{ # if only one cluster
      n_gh <- matrix((n_k^2-n_k)/2)
    }
    
    for (i in sample.int(N,N)){
      nonempty <- which(colSums(matrix(Z[-i,],nrow=N-1)) > 0) # extract non-empty groups after removing node i
      Z <- Z[,nonempty] # remove the group if removing node i makes it empty
      if (length(nonempty)==1){Z <- matrix(Z,N,1)} # ensure the matrix form of Z when K = 1 after removing the possible empty group
      K <- ncol(Z) # update K after removing node i
      Z_Ri <- Z[-i,] # remove node i from Z
      if (K==1){Z_Ri <- matrix(Z_Ri,nrow=N-1,ncol=1)} # Z_Ri becomes vector form if k=1
      z_Ri <- c(Z_Ri%*%1:K) # z after removing node i
      U_Ri <- U[-i,] # U after removing node i
      n_k_Ri_full <- colSums(Z_Ri) # an 1 X K matrix of {n_k} after removing node i
      
      # 1. First obtain f(U|z) terms of the full-conditional distribution by assuming that node i is removed and then is reassigned
      # Construct a N-1 X d*K matrix of latent positions where 
      # the rows which correspond to cluster 1 individuals in the 1st three columns store the corresponding latent positions, and other elements in the 1st three columns are zeros
      # The rows which correspond to cluster 2 individuals in the 2nd three columns store the corresponding latent positions, and other elements in the 2nd three columns are zeros
      # Similarly for all the rest rows and columns
      U_special_matrix_Ri <- matrix(U_Ri,nrow=N-1,ncol=d*K)*matrix(do.call("rbind", replicate(d,Z_Ri,simplify = FALSE)),nrow=N-1,ncol=d*K)
      # Use rbind() function to combine the above matrix with a new K X d*K matrix, where matrix[1,1:3], matrix[2,3:6],matrix[3,7:9],... of this new matrix are U[i,] if d=3 while all other elements are zeros
      U_special_matrix_Ai <- rbind(U_special_matrix_Ri,t(t(matrix(do.call("rbind", replicate(d,diag(K),simplify = FALSE)),nrow=K,ncol=d*K))*U[i,]))
      Z_Ri_spcecial_Ai <- rbind(Z_Ri,diag(K))
      n_k_Ri_full_Ai <- n_k_Ri_full+1
      
      # Ri means remove node i; evaluate each non-empty k'th component of the full-conditional function U\{u_i} | z\{z_i}
      # Ai means add node i; evaluate each non-empty k'th component of the full-conditional function U|z by assuming that z_i is in a non-empty cluster k
      Uk_terms_Ai <- lgamma(alpha1+n_k_Ri_full_Ai*d/2)-n_k_Ri_full_Ai*d/2*log(pi)-d/2*log(omega+n_k_Ri_full_Ai)-
        (n_k_Ri_full_Ai*d/2+alpha1)*log(rowSums(t(Z_Ri_spcecial_Ai)%*%(U_special_matrix_Ai^2))-
                                          rowSums((t(Z_Ri_spcecial_Ai)%*%U_special_matrix_Ai)^2)/(n_k_Ri_full_Ai+omega)+alpha2)
      Uk_terms_Ri <- lgamma(alpha1+n_k_Ri_full*d/2)-n_k_Ri_full*d/2*log(pi)-d/2*log(omega+n_k_Ri_full)-
        (n_k_Ri_full*d/2+alpha1)*log(rowSums(t(Z_Ri)%*%(U_special_matrix_Ri^2))-
                                       rowSums((t(Z_Ri)%*%U_special_matrix_Ri)^2)/(n_k_Ri_full+omega)+alpha2)
      log_update_prob_U <- c(Uk_terms_Ai-Uk_terms_Ri,lgamma(alpha1+d/2)-d/2*log(pi)-d/2*log(omega+1) + alpha1*log(alpha2)+d/2*log(omega)-lgamma(alpha1)-
                               (d/2+alpha1)*log(sum(U[i,]^2)*(1-1/(1+omega))+alpha2)) # here we evaluate the f(U|z\{z_i},z_i)/f(U\{u_i}|z\{z_i}) as well as the case which assigns node i to the new group K+1
      # If we let U_k be the set of latent positions of nodes from cluster k, then some of the terms above correspond to:
      # rowSums(t(Z_Ri)%*%(U_special_matrix_Ri^2)) is {sum(U_k^2)} before adding node i to each cluster
      # rowSums((t(Z_Ri)%*%U_special_matrix_Ri)^2) is {sum(colSums(U_k)^2)} before adding node i to each cluster
      # rowSums(t(Z_Ri_spcecial_Ai)%*%(U_special_matrix_Ai^2)) is {sum(U_k^2)} after adding node i to each cluster
      # rowSums((t(Z_Ri_spcecial_Ai)%*%U_special_matrix_Ai)^2) is {sum(colSums(U_k)^2)} after adding node i to each cluster
      
      # 2. Then we obtain the f(nu|z) term
      nu_gh_Ri <- matrix(nu_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      n_gh_Ri <- matrix(n_gh[nonempty,nonempty],K,K) # remove the group if removing i brings an empty group
      nu_ik <- crossprod(Z_Ri, nu[i,-i]) # {nu_ik}
      z_i <- which(Z[i,] > 0) # this is z_i, e.g., = k; or = integer(0) with length==0 if removing node i leading to an empty cluster
      if (length(z_i)==1){ # if removing node i doesn't bring an empty cluster
        nu_gh_resid <- matrix(0,K,K) # calculate the difference in nu_gh after removing node i
        nu_gh_resid[z_i,] <- nu_ik; nu_gh_resid[,z_i] <- nu_ik
        nu_gh_Ri <- nu_gh_Ri-nu_gh_resid
        
        n_gh_resid <- matrix(0,K,K) # calculate the difference in n_gh after removing node i
        n_gh_resid[z_i,] <- n_k_Ri_full; n_gh_resid[,z_i] <- n_k_Ri_full
        n_gh_Ri <- n_gh_Ri-n_gh_resid
      }
      nu_ik_row_Matrix <- matrix(1,K,1)%*%matrix(nu_ik,1,K) # construct a K X K matrix where each row is nu_ik
      n_k_Ri_full_row_Matrix <- matrix(1,K,1)%*%matrix(n_k_Ri_full,1,K) # construct a K X K matrix where each row is n_k_Ri_full; thus the transport is a K X K matrix where each column is n_k_Ri_full
      log_update_prob_nu <- c(rowSums(lbeta(nu_gh_Ri+nu_ik_row_Matrix+beta1,n_gh_Ri+n_k_Ri_full_row_Matrix-(nu_gh_Ri+nu_ik_row_Matrix)+beta2)-
                                        lbeta(nu_gh_Ri+beta1,n_gh_Ri-nu_gh_Ri+beta2)), # log nu terms for k = 1,2,...,K
                              sum(lbeta(nu_ik+beta1,n_k_Ri_full-nu_ik+beta2)-lbeta(beta1,beta2))) # add the term relevant to the new group K+1
      
      # 3. Obtain the f(z) term
      log_update_prob_z <- c(log(n_k_Ri_full+alpha),log(alpha)+log(MFM_WNK_TruncPois1(N,K+1,alpha))-
                               log(MFM_WNK_TruncPois1(N,K,alpha)))
      
      # 4. Obtain the assignment probability
      log_update_prob <- log_update_prob_U+log_update_prob_nu+log_update_prob_z
      if (!is.null(A)){ # if A is provided
        A_Ri <- A[-i] # remove c_i from node attributes A
        table.z_Ri.A_Ri <- table(z_Ri,factor(A_Ri,levels=1:C)) # it's possible that an empty level might exist in A_Ri, so we use factor(,levels) to avoid such a problem
        log_update_prob <- log_update_prob + c(log(table.z_Ri.A_Ri[,A[i]]+omega_c)-log(n_k_Ri_full+sum_omega_c),log(omega_c)-log(sum_omega_c))
      }
      
      update_prob <- exp(log_update_prob - max(log_update_prob))
      
      # 5. Determine the reassignment and update Z,K,nu_gh,n_gh,n_k
      z_i_new <- which(rmultinom(1,1,update_prob) > 0)
      if(length(z_i) == 1){Z[i,z_i]<-0} # if removing node i doesn't bring an empty group, let Z[i,] be zeros; otherwise, Z[i,] are already zeros
      if (z_i_new==K+1){
        Z <- cbind(Z,0)
        K <- K+1
        Z[i,z_i_new] <- 1
        nu_gh_Ri <- rbind(cbind(nu_gh_Ri,0),0)
        n_gh_Ri <- rbind(cbind(n_gh_Ri,0),0)
        nu_ik <- c(nu_ik,0)
        n_k_Ri_full <- c(n_k_Ri_full,0)
      }else{
        Z[i,z_i_new] <- 1
      }
      
      nu_gh_resid <- matrix(0,K,K) # nu_gh difference by adding node i to group k
      nu_gh_resid[z_i_new,] <- nu_ik; nu_gh_resid[,z_i_new] <- nu_ik
      nu_gh <- nu_gh_Ri + nu_gh_resid
      
      n_gh_resid <- matrix(0,K,K) # n_gh difference by adding node i to group k
      n_gh_resid[z_i_new,] <- n_k_Ri_full; n_gh_resid[,z_i_new] <- n_k_Ri_full
      n_gh <- n_gh_Ri + n_gh_resid
    } # end for i
    
    z <- c(Z%*%c(1:ncol(Z)))
    z_list <- rbind(z_list,z)
    K_list[t+1] <- K
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # TAE step
    if (K == 1){ # if only one cluster, eject
      AE_decision <- 0
      p_eject_p <- 1 # proposal probability of the eject move
      p_absorb_p <- 1-p_eject # proposal probability of the reverse absorb move
    }else if (K == N){ # if N clusters, absorb
      AE_decision <- 1
      p_eject_p <- p_eject 
      p_absorb_p <- 1
    }else{ # with probability p_eject to eject or absorb
      AE_decision <- runif(1)
      p_eject_p <- p_eject # P(eject) in the proposal
      p_absorb_p <- 1-p_eject
    }
    
    if (AE_decision<=p_eject){ # eject move
      g <- sample.int(K,1) # randomly choose one cluster
      n_k <- table(z) # n_k of z
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      U_g_c <- matrix(U[z==g,],ncol=d) # extract all the u_j for current clustering z_j=g
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_c]-log(0.01)))] # Recall: check a ranging from 0.01 to 100 and check n from 1 to N; 0.01 here is set as the Pr(n_{K+1}=0);
      p_E <- rbeta(1,a,a) # sample the probability of the nodes being assigned to cluster K+1
      K_p <- K+1 # proposed number of clusters
      z_p <- z # initialize proposed clustering
      z_p[z_p==g]<-rbinom(n_g_c,1,p_E)*(K_p-g)+g # apply the eject move; for each node j in cluster g, if rbinom() returns 1, assign (K_p-g)+g to z_j; if rbinom() returns 0, assign g back to z_j
      n_k_p <- table(factor(z_p,levels=1:K_p)) # {n_k} of z_p # cannot use Rfast::table here which doesn't return 0 when a group is empty
      n_g_p <- n_k_p[[g]] # extract n_g for proposed clustering after eject move
      # Note that it's possible to have empty cluster in the proposed z and table() can not directly return 0 for empty clusters;
      # Thus we use table(factor(z_p,levels=1:K_p)) to let the empty cluster appear.
      if (n_g_p!=0 & n_g_p!=n_g_c){ # if any of the ejected cluster is not empty
        if (K_p==N){ # If the proposed state has K == N
          p_absorb_p <- 1 # the probability of the backward absorb move for the proposed state is 1
        }
        U_g_p <- matrix(U[z_p==g,],ncol=d) # extract all the u_j for proposed clustering z_j=g after eject move 
        n_Kp_p <- n_k_p[[K_p]] # extract n_Kp for proposed clustering after eject move
        U_Kp_p <- matrix(U[z_p==K_p,],ncol=d) # extract all the u_j for proposed clustering z_j=K_p after eject move 
        log_alpha_E_right <-  alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega)+log(p_absorb_p)-log(p_eject_p)+ # constant terms
          lgamma(alpha1+n_g_p*d/2) - n_g_p*d/2*log(pi) - d/2*log(omega+n_g_p) - (alpha1+n_g_p*d/2)*log(sum(U_g_p^2)-sum(colSums(U_g_p)^2)/(n_g_p+omega)+alpha2)+ # proposed lambda'_g term
          lgamma(alpha1+n_Kp_p*d/2) - n_Kp_p*d/2*log(pi) - d/2*log(omega+n_Kp_p) - (alpha1+n_Kp_p*d/2)*log(sum(U_Kp_p^2)-sum(colSums(U_Kp_p)^2)/(n_Kp_p+omega)+alpha2)- # proposed lambda'_Kp term
          (lgamma(alpha1+n_g_c*d/2) - n_g_c*d/2*log(pi) - d/2*log(omega+n_g_c) - (alpha1+n_g_c*d/2)*log(sum(U_g_c^2)-sum(colSums(U_g_c)^2)/(n_g_c+omega)+alpha2))+ # current lambda^{(t+1)]}_g term; note that z^{(t+1)} was updated from z^{(t)} in the previous step, so current state is z^{(t+1)}
          log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))+sum(log(alpha+0:(n_Kp_p-1)))- # z' MFM prior term
          log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))+ # z MFM prior term
          2*lgamma(a)-lgamma(2*a)+lgamma(2*a+n_g_c)-lgamma(a+n_g_p)-lgamma(a+n_Kp_p)+ # proposal terms
          log(1-2*exp(log_p_0(a,n_g_c))) # normalizing term of truncating at the empty cluster cases
        
        if (!is.null(A)){ # if A is provided
          table.z_c.A <- table(z,A)
          table.z_p.A <- table(z_p,A)
          log_alpha_E_right <- log_alpha_E_right +
            sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)+
            sum(lgamma(table.z_p.A[K_p,]+omega_c))-lgamma(n_Kp_p+sum_omega_c)+lgamma(sum_omega_c)-C*lgamma(omega_c)-
            sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c) # the ratio of A|z term
        }
        
        Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
        nu_gh_p <- t(Z_p)%*%nu%*%Z_p # {nu_gh} of z_p
        diag(nu_gh_p) <- diag(nu_gh_p)/2
        n_gh_p <- outer(n_k_p,n_k_p)-diag(n_k_p) # {n_gh} of z_p
        diag(n_gh_p) <- diag(n_gh_p)/2
        
        changed_gh_position_p <- matrix(FALSE,K_p,K_p) # extract those changed blocks
        changed_gh_position_p[g,] <- TRUE; changed_gh_position_p[,K_p] <- TRUE
        
        changed_nu_gh <- nu_gh[g,]; changed_n_gh <- n_gh[g,]
        changed_nu_gh_p <- nu_gh_p[changed_gh_position_p]; changed_n_gh_p <- n_gh_p[changed_gh_position_p]
        
        log_alpha_E_right <- log_alpha_E_right + sum(lbeta(changed_nu_gh_p+beta1,changed_n_gh_p-changed_nu_gh_p+beta2))-
          sum(lbeta(changed_nu_gh+beta1,changed_n_gh-changed_nu_gh+beta2))-lbeta(beta1,beta2)*K_p
        
        if (log(runif(1)) <= min(0,log_alpha_E_right)){ # accept or not
          z_list[t+1,] <- z <- z_p # update the new clustering
          Z <- Z_p
          K_list[t+1] <- K <- K_p
          acceptance_count_TAE[t] <- 1 # 1 means accept, 0 means not accept
          nu_gh <- nu_gh_p
          n_gh <- n_gh_p
        }# else: stay at the current state and do nothing
      }# else: stay at the current state and do nothing
    }else{ # absorb step
      random_two_picks <- sort(sample.int(K,2))
      g <- random_two_picks[1] # g < h, we let g absorb h
      h <- random_two_picks[2]
      n_k <- table(z)
      n_g_c <- n_k[[g]] # extract n_g for current clustering
      n_h_c <- n_k[[h]] # extract n_h for current clustering
      U_g_c <- matrix(U[z==g,],ncol=d) # extract all the u_j for current clustering z_j=g
      U_h_c <- matrix(U[z==h,],ncol=d) # extract all the u_j for current clustering z_j=h
      K_p <- K-1
      if (K_p==1){ # If the proposed state has K == 1
        p_eject_p <- 1 # the probability of the backward eject move for the proposed state is 1
      }
      z_p <- z # initialize proposed clustering
      z_p[z_p==h] <- g # assign all nodes from cluster h to cluster g
      z_p[z_p>h] <- z_p[z_p>h]-1 # relabel the cluster labels which are greater than h
      n_g_p <- n_g_c+n_h_c # table(z_p)[[g]] # extract n_g for proposed clustering after absorb move
      U_g_p <- matrix(U[z_p==g,],ncol=d) # extract all the u_j for proposed clustering z_j=g after absorb move 
      a <- seq(0.01,100,0.01)[which.min(abs(lookuptable_for_a[,n_g_p]-log(0.01)))] # Recall: check a ranging from 0.01 to 100 and check n from 0 to N 
      log_alpha_A_right <-  -(alpha1*log(alpha2)-lgamma(alpha1)+d/2*log(omega))-log(p_absorb_p)+log(p_eject_p)+ #constant terms
        lgamma(alpha1+n_g_p*d/2) - n_g_p*d/2*log(pi) - d/2*log(omega+n_g_p) - (alpha1+n_g_p*d/2)*log(sum(U_g_p^2)-sum(colSums(U_g_p)^2)/(n_g_p+omega)+alpha2)- # proposed lambda'_g term
        (lgamma(alpha1+n_h_c*d/2) - n_h_c*d/2*log(pi) - d/2*log(omega+n_h_c) - (alpha1+n_h_c*d/2)*log(sum(U_h_c^2)-sum(colSums(U_h_c)^2)/(n_h_c+omega)+alpha2)+ # current lambda^{(t+1)}_h term
           lgamma(alpha1+n_g_c*d/2) - n_g_c*d/2*log(pi) - d/2*log(omega+n_g_c) - (alpha1+n_g_c*d/2)*log(sum(U_g_c^2)-sum(colSums(U_g_c)^2)/(n_g_c+omega)+alpha2))+ # current lambda^{(t+1)}_g term
        log(MFM_WNK_TruncPois1(N,K_p,alpha))+sum(log(alpha+0:(n_g_p-1)))-
        log(MFM_WNK_TruncPois1(N,K,alpha))-sum(log(alpha+0:(n_g_c-1)))-sum(log(alpha+0:(n_h_c-1)))+
        lgamma(2*a)-2*lgamma(a)-lgamma(2*a+n_g_p)+lgamma(a+n_g_c)+lgamma(a+n_h_c)-
        log(1-2*exp(log_p_0(a,n_g_p))) # normalizing term of truncating at the empty cluster cases
      
      if (!is.null(A)){ # if A is provided
        table.z_c.A <- table(z,A)
        table.z_p.A <- table(z_p,A)
        log_alpha_A_right <- log_alpha_A_right +
          sum(lgamma(table.z_p.A[g,]+omega_c))-lgamma(n_g_p+sum_omega_c)-lgamma(sum_omega_c)+C*lgamma(omega_c)-
          sum(lgamma(table.z_c.A[g,]+omega_c))+lgamma(n_g_c+sum_omega_c)-
          sum(lgamma(table.z_c.A[h,]+omega_c))+lgamma(n_h_c+sum_omega_c) # the ratio of A|z term
      }
      
      # Evaluate the ratio of the nu terms next
      n_k_p <- table(z_p) # n_k of z_p
      Z_p <- t(t(matrix(z_p,N,K_p))==(1:K_p))*1 # transfer vector z to matrix Z
      nu_gh_p <- t(Z_p)%*%nu%*%Z_p # {nu_gh} of z_p
      diag(nu_gh_p) <- diag(nu_gh_p)/2
      if (length(n_k_p)>1){ # if absorbing leading to K_p>1
        n_gh_p <- outer(n_k_p,n_k_p)-diag(n_k_p) # {n_gh} of z_p
        diag(n_gh_p) <- diag(n_gh_p)/2
      }else{ # if absorbing leading to K_p=1
        n_gh_p <- matrix((n_k_p^2-n_k_p)/2)
      }
      changed_gh_position <- matrix(FALSE,K,K) # extract those changed blocks
      changed_gh_position[g,] <- TRUE; changed_gh_position[,h] <- TRUE
      
      changed_nu_gh <- nu_gh[changed_gh_position]; changed_n_gh <- n_gh[changed_gh_position]
      changed_nu_gh_p <- nu_gh_p[g,]; changed_n_gh_p <- n_gh_p[g,]
      
      log_alpha_A_right <- log_alpha_A_right + sum(lbeta(changed_nu_gh_p+beta1,changed_n_gh_p-changed_nu_gh_p+beta2))-
        sum(lbeta(changed_nu_gh+beta1,changed_n_gh-changed_nu_gh+beta2))+lbeta(beta1,beta2)*K
      
      if (log(runif(1)) <= min(0,log_alpha_A_right)){ # accept or not
        z_list[t+1,] <- z <- z_p
        Z <- Z_p
        K_list[t+1] <- K <- K_p
        acceptance_count_TAE[t] <- 1
        nu_gh <- nu_gh_p
        n_gh <- n_gh_p
      }# else: stay at the current state and do nothing
    }
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update P by Gibbs
    P <- matrix(0,K,K)
    nu_gh_upperdiag <- upperTriangle(nu_gh,diag=TRUE)
    upperTriangle(P,diag=TRUE) <- rbeta(K*(K+1)/2,nu_gh_upperdiag+beta1,upperTriangle(n_gh,diag=TRUE)-nu_gh_upperdiag+beta2)
    lowerTriangle(P,byrow=TRUE)<-upperTriangle(P)
    P_list[[t+1]] <- P
    #--------------------------------------------------------------------------------------------------------------------------------------------
  } # end T
  return(list(z = z_list,U = U_list,nu = nu_list,K = K_list,X = X_list,P = P_list,beta = beta_list,
              acceptance_count_U = acceptance_count_U, acceptance_count_beta = acceptance_count_beta, acceptance_count_TAE = acceptance_count_TAE))
}

