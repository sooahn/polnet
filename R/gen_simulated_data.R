#'@name gen_simulated_data
#'@param m 
#'@param n
#'@param k
#'@param m.center.cor
#'@param n.center.cor
#'@return A list of simulated data

#'@import MASS
#'@useDynLib polnet, .registration = TRUE
#'@export


# D=2
# n.cluster=4
# clients.center=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5
# legislators.center=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5
# clients=100
# legislators=50
# client_space = NULL
# legislator_space = NULL
# mu = NULL
# Sigma = NULL
# tau = 1
# alpha_popularity = NULL
# Beta_popularity = NULL
# v = 0
# sigma_sq_L = 2
# sigma_sq_P = 3


# D=2
# n.cluster=1
# clients.center=c(0,0)
# legislators.center=c(0,0)
# clients=100
# legislators=50
# client_space = NULL
# legislator_space = NULL
# mu = NULL
# Sigma = NULL
# tau = 1
# alpha_popularity = NULL
# Beta_popularity = NULL
# v = 0
# sigma_sq_L = 2
# sigma_sq_P = 3
# 
# D=1
# n.cluster=1
# clients.center=0
# legislators.center=1
# sigma_sq_L = 1
# sigma_sq_P = 2
# v=3
# tau = 0.5

# Generate clustered LSNM data
random_LSNM_data_cluster <- function(n.cluster=4, clients.center=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5,
                                     legislators.center=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5, D=2, 
                                     clients=100, legislators=50,
                                     Sigma = NULL, tau = 1, v = 0, sigma_sq_L = 2, sigma_sq_P = 3){
  
  # Generate Positive Definite Matrices
  Posdef <- function (n, ev = runif(n, 0, 1)) {
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp) 
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
  }
  
  if (is.null(Sigma)){
    if (D>=2){
      Sigma <- Posdef(n=D)  
    }
    else Sigma <- runif(1, 1, 2)
  }
  
  # Initialize Actor Coordinates
  clients.cor <- matrix(0, nrow=0, ncol=D+1)
  legislators.cor <- matrix(0, nrow=0, ncol=D+1)
  
  # Generate number of units in each clusters randomly
  clients.division <- rmultinom(n=1, size=clients, prob=rep(1/n.cluster,n.cluster))
  legislators.division <- rmultinom(n=1, size=legislators, prob=rep(1/n.cluster,n.cluster))
  
  # Simulated Coordinates Generation
  if (n.cluster>=2){
    for (i in 1:n.cluster){
      # Generate Coordinates for Actor 1
      if (D>=2){
        clients.sigma <- diag(D)*tau
      } else clients.sigma <- tau
      temp.new.clients.cor <- mvrnorm(n=clients.division[i], mu=clients.center[i,], 
                                      Sigma=clients.sigma)
  
      temp.new.clients.cor <- cbind(temp.new.clients.cor, i)
      clients.cor <- rbind(clients.cor, temp.new.clients.cor)
      
      # Generate Coordinates for Actor 2
      legislators.sigma <- Sigma
      temp.new.legislators.cor <- mvrnorm(n=legislators.division[i], mu=legislators.center[i,], 
                                Sigma=legislators.sigma)
      temp.new.legislators.cor <- cbind(temp.new.legislators.cor, i)
      legislators.cor <- rbind(legislators.cor, temp.new.legislators.cor)
    }
  } else {
    
    # Generate Coordinates for Actor 1
    if (D>=2){
      clients.sigma <- diag(D)*tau
    } else clients.sigma <- tau
    temp.new.clients.cor <- mvrnorm(n=clients.division, mu=clients.center, 
                                    Sigma=clients.sigma)
    
    temp.new.clients.cor <- cbind(temp.new.clients.cor, 1)
    clients.cor <- rbind(clients.cor, temp.new.clients.cor)
    
    # Generate Coordinates for Actor 2
    legislators.sigma <- Sigma
    temp.new.legislators.cor <- mvrnorm(n=legislators.division, mu=legislators.center, 
                                        Sigma=legislators.sigma)
    temp.new.legislators.cor <- cbind(temp.new.legislators.cor, 1)
    legislators.cor <- rbind(legislators.cor, temp.new.legislators.cor)
    
  }
    
  # Randomly Generate Actors Popularity
  clients.popularity <- rnorm(clients, v, sqrt(sigma_sq_L))
  legislators.popularity <- rnorm(legislators, 0, sqrt(sigma_sq_P))

  
  # Call existing function
  LSNM_data <- random_LSNM_data(D=D, clients=clients, legislators = legislators, client_space = clients.cor[,-(D+1)],
                        legislator_space = legislators.cor[,-(D+1)], alpha_popularity = clients.popularity,
                        Beta_popularity = legislators.popularity)
  
  # return
  ret.list <- list(LSNM_data=LSNM_data, clients.popularity=clients.popularity, legislators.popularity=legislators.popularity, clients.cor=clients.cor,
                   legislators.cor=legislators.cor, D=D, clients.center=clients.center,
                   legislators.center=legislators.center, clients=clients,
                   legislators=legislators, Sigma=Sigma, tau=tau, v=v, sigma_sq_L = sigma_sq_L, 
                   sigma_sq_P = sigma_sq_P)
}
  

# m <- clients
# n <- legislators
# m.cor <- clients.cor
# n.cor <- legislators.cor 
# m.popularity <- clients.popularity
# n.popularity <- legislators.popularity
# 
# # Generate Simulated Data
# gen_simulated_data <- function(m=100, n=50, k=4, dimension=2,
#                                m.center.cor=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5,
#                                n.center.cor=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5){
#   # Number of Actor 1: m
#   # Number of Actor 2: n
#   # Number of Clusters: k
#   # Center Coordinates of Clusters are m.center.cor and n.center.cor
# 
#   # Initialize Actor Coordinates
#   m.cor <- matrix(0, nrow=0, ncol=dimension+1)
#   n.cor <- matrix(0, nrow=0, ncol=dimension+1)
#   
#   # Generate randomly number of units in each clusters 
#   m.division <- rmultinom(n=1, size=m, prob=rep(1/k,k))
#   n.division <- rmultinom(n=1, size=n, prob=rep(1/k,k))
#   
#   # Simulated Coordinates Generation
#   for (i in 1:k){
#     # Generate Coordinates for Actor 1
#     if (dimension>=2){
#       temp.sigma <- Posdef(dimension)
#     } else temp.sigma <- runif(1,1,2)
#     temp.new.m.cor <- mvrnorm(n=m.division[i], mu=m.center.cor[i,], 
#             Sigma=temp.sigma)
#     
#     temp.new.m.cor <- cbind(temp.new.m.cor, i)
#     m.cor <- rbind(m.cor, temp.new.m.cor)
#     
#     # Generate Coordinates for Actor 2
#     if (dimension>=2){
#       temp.sigma <- Posdef(dimension)
#     } else temp.sigma <- runif(1,1,2)
#     temp.new.n.cor <- mvrnorm(n=n.division[i], mu=n.center.cor[i,], 
#                               Sigma=temp.sigma)
#     temp.new.n.cor <- cbind(temp.new.n.cor, i)
#     n.cor <- rbind(n.cor, temp.new.n.cor)
#   }
#   
#   # Randomly Generate Actors Popularity
#   m.popularity <- runif(m, 2, 5)
#   n.popularity <- runif(n, 0, 0)
#   
#   # Initialize Adjacency Matrix
#   A <- matrix(0, nrow=m, ncol=n)
#   
#   # Define Pairwise Distrance Function
#   vectorized_pdist <- function(A,B) {
#     an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
#     bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
#   
#     m = nrow(A)
#     n = nrow(B)
#   
#     tmp = matrix(rep(an, n), nrow=m) 
#     tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
#     sqrt( tmp - 2 * tcrossprod(A,B) )
#   }
#   
#   # Compute Distance Matrix
#   distance <- vectorized_pdist(cbind(m.cor[,1], m.cor[,2]), 
#                         cbind(n.cor[,1], n.cor[,2]))
#   
#   # Compute Mean Matrix
#   mean <- -distance+matrix(rep(m.popularity, n), nrow=m)+
#     matrix(rep(n.popularity, each=m), ncol=n)
#   
#   mean <- exp(mean)
#   
#   # Generate Adjacency Matrix
#   mean.vec <- as.vector(mean)
#   A <- matrix(rpois(m*n, mean), nrow=m)
#   
#   colnames(m.cor) <- c(paste0("d",1:dimension),"group")
#   colnames(n.cor) <- c(paste0("d",1:dimension),"group")
#   
#   ret <- list(A=A,
#               m.cor=m.cor,
#               n.cor=n.cor)
#   
#   return(ret)
# }
