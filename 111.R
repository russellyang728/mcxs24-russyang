set.seed(2048)
RW11 <- arima.sim(model= list(order = c(0, 1, 0)), n=1000, mean=0, sd=1)
RW22 <- arima.sim(model= list(order = c(0, 1, 0)), n=1000, mean=0, sd=1)

RW.1 <- cbind(RW11,RW22)

Y <- RW.1[2:nrow(RW.1),]

Y            = RW.1[2:nrow(RW.1),]
X            = matrix(1,nrow(Y),1)
X            = cbind(X,RW.1[2: nrow(RW.1)-1,])
N            = 2            # number of variables
p            = 1            # number of lags
K            = 1+p*N
S            = 5000         # sample size

sign.restrictions = c(1,1)

A.hat = solve(t(X)%*%X)%*%t(X)%*%Y                
Sigma.hat = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/nrow(Y)

kappa.1     <- 0.02^2
kappa.2     <- 100

A.prior     <- matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] <- diag(N)
V.prior     <- diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N)))
S.prior     <- diag(diag(Sigma.hat))
nu.prior    <- N+1
m_a =1

Sigma.posterior.store    <- array(NA, c(N,N,S))
A.posterior.store        <- array(NA, c((1+p*N),N,S))
B0.posterior.store       <- array(NA,c(N,N,S))
B1.posterior.store       <- array(NA,c(N,K,S))

compute_sbar2 = function(s2, Sigma, V) {
  sbar2 = 1/s2
  for (i in 1:nrow(Sigma)) {
    sbar2 = sbar2 + 1/(Sigma[i,i]*V[i+1,i+1])
  }
  return(1/sbar2)
}


compute_a_posterior = function(sbar2, Sigma,V,A){
  a_posterior = m_a/sbar2
  n=nrow(Sigma)
  A_1 = A[2:(n-1),]
  for(i in 1:n) {
    a_posterior = a_posterior + A_1[i,i]/(Sigma[i,i]*V[i+1,i+1])
  }
  return(sbar2*a_posterior)
}

s2 = 1
a.posterior.store = numeric(S)
Sigma.posterior   <- rWishart(1, df=nu.prior, Sigma=S.prior)[,,1]
V.bar = V.prior
A.posterior = A.prior

for (s in 1:S){
  sbar2 = compute_sbar2(s2, Sigma.posterior,V.bar)
  a_posterior = compute_a_posterior(sbar2, Sigma.posterior,V.bar, A.posterior)
  a.posterior = rnorm(1, a_posterior, sqrt(sbar2))
  a.posterior.store[s] = a.posterior
  
  # Matrix normal-inverse Wishart posterior parameters
  V.bar.inv <- t(X)%*%X + diag(1/diag(V.prior))
  V.bar     <- solve(V.bar.inv)
  A.bar     <- V.bar%*%(t(X)%*%Y + diag(1/diag(V.prior))%*%(A.prior*a.posterior))
  nu.bar    <- nrow(Y) + nu.prior
  S.bar     <- S.prior + t(Y)%*%Y + t((A.prior*a.posterior))%*%diag(1/diag(V.prior))%*%(A.prior*a.posterior) - t(A.bar)%*%V.bar.inv%*%A.bar
  S.bar.inv <- solve(S.bar)
  
  
  # Draw Posterior distribution
  Sigma.posterior   <- rWishart(1, df=nu.bar, Sigma=S.bar.inv)[,,1]
  Sigma.posterior.store[,,s]   <- solve(Sigma.posterior)
  
  # Draw A from matrix-variate normal distribution
  A.posterior = matrix(mvtnorm::rmvnorm(1, mean=as.vector(A.bar), sigma=Sigma.posterior.store[,,s]%x%V.bar), ncol=N)
  A.posterior.store[,,s] = A.posterior
  
  ## Draw from the Structural Form
  cholSigma.s        <- chol(Sigma.posterior.store[,,s])
  B0.posterior.store[,,s]  <- solve(t(cholSigma.s)) 
  B1.posterior.store[,,s]  <- B0.posterior.store[,,s]%*%t(A.posterior.store [,,s])
}

# Identification via sign restrictions 
R1 <- diag(sign.restrictions)

# Storage matrices for Q identified estimates
i.vec <- c()
Q.store      <- array(NA,c(N,N,(S)))
B0.store    <- array(NA,c(N,N,(S)))
B1.store     <- array(NA,c(N,K,(S)))
A.store     <- array(NA,c(K,N,S))
Sigma.store  <- array(NA,c(N,N,S))

for (s in 1:S){
  A             <- A.posterior.store[,,s]
  Sigma         <- Sigma.posterior.store[,,s]
  B0.tilde      <- B0.posterior.store[,,s]
  B1.tilde      <- B1.posterior.store[,,s]
  
  sign.restrictions.do.not.hold = TRUE
  i=1
  while (sign.restrictions.do.not.hold){
    X           <- matrix(rnorm(N*N),N,N)         
    QR          <- qr(X, tol = 1e-10)
    Q           <- qr.Q(QR,complete=TRUE)
    R           <- qr.R(QR,complete=TRUE)
    Q           <- t(Q %*% diag(sign(diag(R))))
    B0          <- Q%*%B0.tilde                    
    B1          <- Q%*%B1.tilde                   
    B0.inv      <- solve(B0)      
    check       <- all(c(B0[1,1], B0[2,2]) > 0)
    
    if (check){sign.restrictions.do.not.hold=FALSE}
    i=i+1 
  }
  i.vec <- c(i.vec,i) 
  
  Q.store[,,s] <- Q
  B0.store[,,s] <- B0
  B0.mean <- apply(B0.store,1:2,mean)
  B1.store[,,s] <- B1
  B1.mean <- apply(B1.store,1:2,mean)
  A.store[,,s]       <- A
  A.mean             <- apply(A,1:2,mean)
  Sigma.store[,,s]   <- Sigma
  Sigma.mean         <- apply(Sigma,1:2,mean)
}
