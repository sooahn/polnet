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
  return(Z)}

# Generate Simulated Data
gen_simulated_data <- function(m=100, n=50, k=4,
                               m.center.cor=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5,
                               n.center.cor=rbind(c(-1,-1), c(-1, 1), c(1, 1), c(1, -1))*5){
  # Number of Actor 1: m
  # Number of Actor 2: n
  # Number of Clusters: k
  # Center Coordinates of Clusters are m.center.cor and n.center.cor

  # Initialize Actor Coordinates
  m.cor <- matrix(0, nrow=0, ncol=3)
  n.cor <- matrix(0, nrow=0, ncol=3)
  
  # Generate randomly number of units in each clusters 
  m.division <- rmultinom(n=1, size=m, prob=rep(1/k,k))
  n.division <- rmultinom(n=1, size=n, prob=rep(1/k,k))
  
  # Simulated Coordinates Generation
  for (i in 1:k){
    # Generate Coordinates for Actor 1
    temp.new.m.cor <- mvrnorm(n=m.division[i], mu=m.center.cor[i,], 
            Sigma=Posdef(2))
    temp.new.m.cor <- cbind(temp.new.m.cor, i)
    m.cor <- rbind(m.cor, temp.new.m.cor)
    
    # Generate Coordinates for Actor 2
    temp.new.n.cor <- mvrnorm(n=n.division[i], mu=n.center.cor[i,], 
                              Sigma=Posdef(2))
    temp.new.n.cor <- cbind(temp.new.n.cor, i)
    n.cor <- rbind(n.cor, temp.new.n.cor)
  }
  
  # Randomly Generate Actors Popularity
  m.popularity <- runif(m, 0, 5)
  n.popularity <- runif(n, 0, 5)
  
  # Initialize Adjacency Matrix
  A <- matrix(0, nrow=m, ncol=n)
  
  # Define Pairwise Distrance Function
  vectorized_pdist <- function(A,B) {
    an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
    bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
    m = nrow(A)
    n = nrow(B)
  
    tmp = matrix(rep(an, n), nrow=m) 
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    sqrt( tmp - 2 * tcrossprod(A,B) )
  }
  
  # Compute Distance Matrix
  distance <- vectorized_pdist(cbind(m.cor[,1], m.cor[,2]), 
                        cbind(n.cor[,1], n.cor[,2]))
  
  # Compute Mean Matrix
  mean <- -distance+matrix(rep(m.popularity, n), nrow=m)+
    matrix(rep(n.popularity, each=m), ncol=n)
  
  mean <- exp(mean)
  
  # Generate Adjacency Matrix
  mean.vec <- as.vector(mean)
  A <- matrix(rpois(m*n, mean), nrow=m)
  
  # Plot M and N
  # plot(m.cor[,1], m.cor[,2], col=c("red", "blue", "green", "black")[m.cor[,3]], main="Actor 1")
  # plot(n.cor[,1], n.cor[,2], col=c("red", "blue", "green", "black")[n.cor[,3]], main="Actor 2")
  
  colnames(m.cor) <- c("x","y","group")
  colnames(n.cor) <- c("x","y","group")
  
  ret <- list(A=A,
              m.cor=m.cor,
              n.cor=n.cor)
  
  return(ret)
}
