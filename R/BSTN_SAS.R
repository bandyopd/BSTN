#' Fit Bayesian Skewed Tensor Normal (BSTN) model with tensor spike-and-slab prior to GAAD data
#' Lee et al. (2021+) Bayesian Regression Analysis of Skewed Tensor Response
#'
#' @param Y t x s x b x n array of skewed tensor response
#' @param X p x n matrix of covariates
#' @param vecy tsb x n matrix: vectorized tensor
#' @param n.burn burn-in period
#' @param n.save number of posterior samples
#' @param thin thinning size
#'
#' @return Returns a list with the following components:
#' \item rho: posterior samples of correlation of each mode of tensor response (3 x n.save)
#' \item sigma.sq: posterior samples of variance parameter (b x n.save)
#' \item lam.est: posterior samples of skewness parameters (b x n.save)
#' \item est.est: posterior samples of common effects of covariates (p x n.save)
#' \item omega: posterior samples of (zeros/ones) that particular element is
#'              included in the gamma.est (t x s x b x p x n.save)
#' \item gamma.est: posterior samples of sparsity elements that
#'                  identify different effects of each tooth-sites (t x s x b x p x n.save)
#'
#' @export
#'

BSTN_SAS <- function(Y,X,vecy, n.burn = 10, n.save = 100, thin = 1){

  t = dim(Y)[1]; s = dim(Y)[2]; b = dim(Y)[3]; p = dim(X)[1]; n = dim(X)[2];

  # Store MCMC out

  B.est.save <- array(NA, c(p, t*s*b, n.save))
  rho.save <- matrix(NA, 3, n.save)
  sigma.sq.save <- matrix(NA, 2, n.save)
  lam.est.save <- matrix(NA, 2, n.save)
  omega.save <- array(NA, dim = c(t,s,b,p,n.save))
  eta.save <- matrix(NA, p, n.save)
  gamma.est.save <- array(NA, dim = c(t,s,b,p,n.save))

  source("./functions_tensor.R")

  #load required packages
  library(expm); library(Matrix); library(matrixcalc); library(LaplacesDemon)
  library(MASS); library(msm); library(truncnorm); library(abind)
  library(magrittr); library(doParallel); library(TruncatedNormal)
  registerDoParallel(cores=2)

  #---------------------------------------------------------------------------------

  # Update rho1, rho2 & rho3 (M-H algorithm) using closed form of equicorrelation assumption for R1 & R2 & R3

  rho.update <- function(t,s,rho,vecy,X,B.est,lam.est,W,sigma.sq){

    rho.prop = c(rbeta(1,2.5,2), rbeta(1,2.5,2), rbeta(1,2.5,2)) #(rho1.prop, rho2.prop, rho3.prop)

    det.curr = ((sigma.sq[1]*sigma.sq[2]*((1-rho[3])^(b-1))*(1 + (b-1)*rho[3]))^(t*s))*((((1-rho[2])^(s-1))*(1 + (s-1)*rho[2]))^(t*b))*((((1-rho[1])^(t-1))*(1 + (t-1)*rho[1]))^(s*b))
    det.prop = ((sigma.sq[1]*sigma.sq[2]*((1-rho.prop[3])^(b-1))*(1 + (b-1)*rho.prop[3]))^(t*s))*((((1-rho.prop[2])^(s-1))*(1 + (s-1)*rho.prop[2]))^(t*b))*((((1-rho.prop[1])^(t-1))*(1 + (t-1)*rho.prop[1]))^(s*b))

    inv.curr = kron(diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2])))%*%((1/(1-rho[3]))*(diag(b) - (rho[3]/ (1 + (b-1)*rho[3]))*matrix(1,b,b)))%*%diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2]))), ((1/(1-rho[2]))*(diag(s) - (rho[2]/ (1 + (s-1)*rho[2]))*matrix(1,s,s))), ((1/(1-rho[1]))*(diag(t) - (rho[1]/ (1 + (t-1)*rho[1]))*matrix(1,t,t)))) # \{D_{\sigma}^{-1}R_{\rho_3}^{-1}D_{\sigma}^{-1} \otimes R_{\rho_2}^{-1} \otimes R_{\rho_1}^{-1}\}
    inv.prop = kron(diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2])))%*%((1/(1-rho.prop[3]))*(diag(b) - (rho.prop[3]/ (1 + (b-1)*rho.prop[3]))*matrix(1,b,b)))%*%diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2]))), ((1/(1-rho.prop[2]))*(diag(s) - (rho.prop[2]/ (1 + (s-1)*rho.prop[2]))*matrix(1,s,s))), ((1/(1-rho.prop[1]))*(diag(t) - (rho.prop[1]/ (1 + (t-1)*rho.prop[1]))*matrix(1,t,t))))

    inv.R21.curr = kron( ((1/(1-rho[2]))*(diag(s) - (rho[2]/ (1 + (s-1)*rho[2]))*matrix(1,s,s))), ((1/(1-rho[1]))*(diag(t) - (rho[1]/ (1 + (t-1)*rho[1]))*matrix(1,t,t))) )
    inv.R21.prop = kron( ((1/(1-rho.prop[2]))*(diag(s) - (rho.prop[2]/ (1 + (s-1)*rho.prop[2]))*matrix(1,s,s))), ((1/(1-rho.prop[1]))*(diag(t) - (rho.prop[1]/ (1 + (t-1)*rho.prop[1]))*matrix(1,t,t))) )

    inv.R3.curr = ((1/(1-rho[3]))*(diag(b) - (rho[3]/ (1 + (b-1)*rho[3]))*matrix(1,b,b)))
    inv.R3.prop = ((1/(1-rho.prop[3]))*(diag(b) - (rho.prop[3]/ (1 + (b-1)*rho.prop[3]))*matrix(1,b,b)))

    logdens.curr = -0.5*n*log(det.curr) - 0.5*sum( (t(vecy) - t(X)%*%B.est - t(W)%*%(kron(diag(lam.est,2,2),diag(t*s))) )%*%inv.curr%*%(vecy - t(B.est)%*%X - (kron(diag(lam.est,2,2),diag(t*s))%*%W) ) )
    logdens.prop = -0.5*n*log(det.prop) - 0.5*sum( (t(vecy) - t(X)%*%B.est - t(W)%*%(kron(diag(lam.est,2,2),diag(t*s))) )%*%inv.prop%*%(vecy - t(B.est)%*%X - (kron(diag(lam.est,2,2),diag(t*s))%*%W) ) )

    logratio = logdens.prop - logdens.curr
    if(log(runif(1)) > logratio) {rho = rho} else {rho = rho.prop}
    if(log(runif(1)) > logratio) {inv.Sigma = inv.curr} else {inv.Sigma = inv.prop}
    if(log(runif(1)) > logratio) {inv.R21 = inv.R21.curr} else {inv.R21 = inv.R21.prop}
    if(log(runif(1)) > logratio) {inv.R3 = inv.R3.curr} else {inv.R3 = inv.R3.prop}
    if(log(runif(1)) > logratio) {det = det.curr} else {det = det.prop}


    return(list(rho, inv.Sigma, inv.R21, inv.R3, det))

  }

  # Update variance parameter sigma.sq.est = c(sigma1.est, sigma2.est)
  # Apply Theorem 1 (Marginalization)

  sigma.sq.update <- function (t,s,b,n,Y,X,B.est,lam.est,W,inv.R21,inv.R3,sigma.sq,rho){

    g1 = 2; g2 = 2; #prior distribution for inv.sigma.sq ~ Ga(2,2)
    Sww <-  foreach(l = 1:n,.combine='+') %dopar% {crossprod(W[1:(t*s),l])}
    #S <- sum( ( mat(Y[,,1,],3) - t(X)%*%B.est[,1:(t*s)] - lam.est[1]*t(W[1:(t*s),]) )%*%kron(inv.R3,inv.R21)[1:(t*s),1:(t*s)]%*%t( mat(Y[,,1,],3) - t(X)%*%B.est[,1:(t*s)] - lam.est[1]*t(W[1:(t*s),]) ) )
    S <- foreach(l = 1:n,.combine='+') %dopar% { ( c(Y[,,1,l]) - t(X)[l,]%*%B.est[,1:(t*s)] - lam.est[1]*t(W[1:(t*s),l]) )%*%kron(inv.R3,inv.R21)[1:(t*s),1:(t*s)]%*%t( c(Y[,,1,l]) - t(X)[l,]%*%B.est[,1:(t*s)] - lam.est[1]*t(W[1:(t*s),l]) ) }
    inv.sigma.sq <- rgamma(1, g1 + n*t*s, 0.5*S + 0.5*Sww + g2)
    sigma1.sq <- 1/inv.sigma.sq

    Sww2 <- foreach(l = 1:n,.combine='+') %dopar% {crossprod(W[(t*s + 1):(2*t*s),l])}
    #S2 <- sum( ( mat(Y[,,2,],3) - t(X)%*%B.est[,(t*s + 1):(2*t*s)] - lam.est[2]*t(W[(t*s + 1):(2*t*s),]) )%*%kron(inv.R3,inv.R21)[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%t( mat(Y[,,2,],3) - t(X)%*%B.est[,(t*s + 1):(2*t*s)] - lam.est[2]*t(W[(t*s + 1):(2*t*s),]) ) )
    S2 <- as.numeric(foreach(l = 1:n,.combine='+') %dopar% { ( c(Y[,,2,l]) - t(X)[l,]%*%B.est[,(t*s + 1):(2*t*s)] - lam.est[2]*t(W[(t*s + 1):(2*t*s),l]) )%*%kron(inv.R3,inv.R21)[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%t( c(Y[,,2,l]) - t(X)[l,]%*%B.est[,(t*s + 1):(2*t*s)] - lam.est[2]*t(W[(t*s + 1):(2*t*s),l]) ) })
    inv.sigma.sq2 <- rgamma(1, g1 + n*t*s, 0.5*S2 + 0.5*Sww2 + g2)
    sigma2.sq <- 1/inv.sigma.sq2

    return(c(sigma1.sq,sigma2.sq))
  }

  #---------------------------------------------------------------------------------
  # Update skewness parameter lam.est = c(lam1.est, lam2.est)
  # Apply Theorem 1 (Marginalization)

  lam.est.update <- function (W,inv.Sigma,B.est,X,Y,t,s,b,n){

    Swinvw <- foreach(l = 1:n,.combine='+') %dopar% {t(W[1:(t*s),l])%*%inv.Sigma[1:(t*s),1:(t*s)]%*%W[1:(t*s),l]}
    A.lam <- (b^2)*Swinvw + 1
    B.lam <- foreach(l = 1:n,.combine='+') %dopar% {( c(Y[,,1,l]) - t(X)[l,]%*%B.est[,1:(t*s)])%*%inv.Sigma[1:(t*s),1:(t*s)]%*%W[1:(t*s),l]*(b^2) + 4 }
    lam1.est <- rnorm(1, mean = B.lam/A.lam, sd = 1/(sqrt(A.lam)))

    Swinvw2 <- foreach(l = 1:n,.combine='+') %dopar% {t(W[(t*s + 1):(2*t*s),l])%*%inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%W[(t*s + 1):(2*t*s),l]}
    A.lam2 <- (b^2)*Swinvw2 + 1
    B.lam2 <- foreach(l = 1:n,.combine='+') %dopar% {( c(Y[,,2,l]) - t(X)[l,]%*%B.est[,(t*s + 1):(2*t*s)])%*%inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%W[(t*s + 1):(2*t*s),l]*(b^2) + 4}
    lam2.est <- rnorm(1, mean = B.lam2/A.lam2, sd = 1/(sqrt(A.lam2)))

    return(c(lam1.est,lam2.est))
  }

  #---------------------------------------------------------------------------------
  # Update W = abs(Z_2) : Note that the each component is sampled from univariate truncated normal distribution

  W.update <- function(t,s,b,n, lam.est, inv.Sigma, B.est, Y, X, sigma.sq, W){

      for (jkl in 1:(t*s*b)){
        for (N in 1:n){

          D <-  kron(diag(c(lam.est^2)), diag(t*s))%*%inv.Sigma + diag(t*s*b); D <- solve(D)
          E <- kron(diag(c(lam.est^2)), diag(t*s))%*%inv.Sigma%*%( t(mat(Y,4)) - t(B.est)%*%X)
          W[jkl,N] <- truncnorm::rtruncnorm(1, a=0, b=Inf, mean = c(D%*%E[,N])[jkl], sd = sqrt(D[jkl,jkl]) )

        }
      }

    return(W)

  }

  #---------------------------------------------------------------------------------
  # Update eta

  eta.est.update <- function(X,vecy,inv.Sigma,W,lam.est, t,s,p){

    eta.xy <- X%*%t(vecy)%*%inv.Sigma
    eta.xw <- X%*%t(W)%*%inv.Sigma*lam.est
    eta.xx <- X%*%t(X)
    eta.A.inv <- solve(kron(eta.xx,inv.Sigma) + 10*diag(t*s*b*p) )
    vec.eta <- MASS::mvrnorm(1, mu = eta.A.inv%*%(as.vector(t(eta.xy)) - as.vector(t(eta.xw))), Sigma = eta.A.inv, tol = 1e-3)
    eta.est <- apply(mat(array(vec.eta, dim = c(t,s,p)),3),1,mean)

    return(eta.est)

  }


  #---------------------------------------------------------------------------------
  # Update gamma and beta
  # Apply Theorem 1 (Marginalization)

  B.est.update <- function(X,Y,W,lam.est, inv.Sigma, inv.R1, inv.R2, t,s,p, omega, eta.est){

    psi <- rbeta(1, 0.1 + sum(omega), 0.1 + t*s*b*p - sum(omega) )
    B <- array(NA, dim = c(t,s,2,p)) ;

    Sxx <- matrix(NA,p,1); Sxy <- matrix(NA,p,t*s); Sxw <- matrix(NA,p,t*s); piup <- matrix(NA,p,1); gamma.est <- array(NA, dim = c(t,s,p))
    Sxx2 <- matrix(NA,p,1); Sxy2 <- matrix(NA,p,t*s); Sxw2 <- matrix(NA,p,t*s); piup2 <- matrix(NA,p,1); gamma.est2 <- array(NA, dim = c(t,s,p))
    A.inv <- array(NA, dim = c(p,t*s,t*s)); A.inv0 <- array(NA, dim = c(p,t*s,t*s)); A.inv2 <- array(NA, dim = c(p,t*s,t*s)); A.inv02 <- array(NA, dim = c(p,t*s,t*s));
    exponent0 <- matrix(NA,p,1); log.det0 <- matrix(NA,p,1);  l0 <- matrix(NA,p,1); exponent1 <- matrix(NA,p,1); log.det1 <- matrix(NA,p,1);  l1 <- matrix(NA,p,1);
    exponent02 <- matrix(NA,p,1); log.det02 <- matrix(NA,p,1);  l02 <- matrix(NA,p,1); exponent12 <- matrix(NA,p,1); log.det12 <- matrix(NA,p,1);  l12 <- matrix(NA,p,1);

    for (i1 in 1:t){
      for (i2 in 1:s){
        for (j in 1:p){

          Sxx[j,] <- X[j,]%*%t(X)[,j]
          Sxy[j,] <- X[j,]%*%mat(Y[,,1,],3)%*%inv.Sigma[1:(t*s),1:(t*s)]
          Sxw[j,] <- X[j,]%*%t(W[1:(t*s),])%*%inv.Sigma[1:(t*s),1:(t*s)]*lam.est[1]
          A.inv[j,,] <- solve(kron(Sxx[j,],inv.Sigma[1:(t*s),1:(t*s)]) + 100*diag(t*s*1))
          A.inv0[j,,] <- solve(kron(Sxx[j,],inv.Sigma[1:(t*s),1:(t*s)]))

          log.det0[j,] <- (-(t*s*1)/2)*log(2*pi) + 0.5*logdet(A.inv0[j,,]) + (-(t*s*1)/2)*log(2)
          exponent0[j,] <- (-0.5)*t(as.vector(t(Sxy[j,])) - as.vector(t(Sxw[j,])))%*%A.inv0[j,,]%*%(as.vector(t(Sxy[j,])) - as.vector(t(Sxw[j,])))

          log.det1[j,] <- (-(t*s*1)/2)*log(2*pi) + 0.5*logdet(A.inv[j,,]) + (-(t*s*1)/2)*log(2)
          exponent1[j,] <- (-0.5)*t(as.vector(t(Sxy[j,])) - as.vector(t(Sxw[j,])))%*%A.inv[j,,]%*%(as.vector(t(Sxy[j,])) - as.vector(t(Sxw[j,])))

          l1[j,] <- (log.det1[j,] + exponent1[j,]); l0[j,] <- (log.det0[j,] + exponent0[j,])
          piup[j,] <- exp(l0[j,] - l1[j,])

          omega[i1,i2,1,j] <- rbinom(1,1, psi /  (psi + (1-psi)*piup[j,]) )
          gamma.est[,,j] <- array(MASS::mvrnorm(1, mu = A.inv[j,,]%*%(as.vector(t(Sxy[j,])) - as.vector(t(Sxw[j,]))), Sigma = A.inv[j,,], tol = 1e-3) ,dim = c(t,s))

          Sxy2[j,] <- X[j,]%*%mat(Y[,,2,],3)%*%inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]
          Sxw2[j,] <- X[j,]%*%t(W[(t*s + 1):(2*t*s),])%*%inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]*lam.est[2]
          A.inv2[j,,] <- solve(kron(Sxx[j,],inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]) + 100*diag(t*s*1))
          A.inv02[j,,] <- solve(kron(Sxx[j,],inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]))

          log.det02[j,] <- (-(t*s*1)/2)*log(2*pi) + 0.5*logdet(A.inv02[j,,]) + (-(t*s*1)/2)*log(2)
          exponent02[j,] <- (-0.5)*t(as.vector(t(Sxy2[j,])) - as.vector(t(Sxw2[j,])))%*%A.inv02[j,,]%*%(as.vector(t(Sxy2[j,])) - as.vector(t(Sxw2[j,])))

          log.det12[j,] <- (-(t*s*1)/2)*log(2*pi) + 0.5*logdet(A.inv2[j,,]) + (-(t*s*1)/2)*log(2)
          exponent12[j,] <- (-0.5)*t(as.vector(t(Sxy2[j,])) - as.vector(t(Sxw2[j,])))%*%A.inv2[j,,]%*%(as.vector(t(Sxy2[j,])) - as.vector(t(Sxw2[j,])))

          l12[j,] <- (log.det12[j,] + exponent12[j,]); l02[j,] <- (log.det02[j,] + exponent02[j,])
          piup2[j,] <- exp(l02[j,] - l12[j,])

          omega[i1,i2,2,j] <- rbinom(1,1, psi /  (psi + (1-psi)*piup2[j,]) )
          gamma.est2[,,j] <- array(MASS::mvrnorm(1, mu = A.inv2[j,,]%*%(as.vector(t(Sxy2[j,])) - as.vector(t(Sxw2[j,]))), Sigma = A.inv2[j,,], tol = 1e-3) ,dim = c(t,s))

        }
      }
    }

    for (i1 in 1:t){
      for (i2 in 1:s){
        for (j in 1:p){
          if(omega[i1,i2,1,j]=="0"){gamma.est[i1,i2,j] <- 0}
          B[i1,i2,1,j] <- eta.est[j] + gamma.est[i1,i2,j]
        }}}



    for (i1 in 1:t){
      for (i2 in 1:s){
        for (j in 1:p){
          if(omega[i1,i2,2,j]=="0"){gamma.est2[i1,i2,j] <- 0}
          B[i1,i2,2,j] <- eta.est[j] + gamma.est2[i1,i2,j]
        }}}

    B.est <- mat(B,4)

    return(list(B.est, omega, gamma.est))

  }

  #initial values

  B.est <- solve(X%*%t(X))%*%X%*%t(vecy) # OLS
  W <- matrix(0.2, t*s*b,n) # vecotrize W = |Z_{2i}|
  rho <- matrix(c(0.6,0.6,0.6), 3, 1)
  sigma.sq <- c(1,1)
  lam.est <- c(1,1)
  omega <- array(rbinom(t*s*b*p,1,0.05),dim = c(t,s,b,p))

  R1 <- matrix(0,t,t)
  for(j in 1:t){
    for(k in 1:t){
      if (j != k ){
        R1[j,k] = 0.6
      }
      else{
        R1[j,k] <- 1
      }
    }
  }

  R2 <- matrix(0,s,s)
  for(j in 1:s){
    for(k in 1:s){
      if (j != k ){
        R2[j,k] = 0.6
      }
      else{
        R2[j,k] <- 1
      }
    }
  }

  R3 <- matrix(0,b,b)
  for(j in 1:b){
    for(k in 1:b){
      if (j != k ){
        R3[j,k] = 0.6
      }
      else{
        R3[j,k] <- 1
      }
    }
  }

  # fill missing responses with the new missing values

  for (N in 1:n){
    vecy[,N] = ifelse (is.na(vecy[,N]), MASS::mvrnorm(length(delta_p[,N][(delta_p[,N]==1)]), mu = t(B.est)%*%X[,N] + kron(diag(lam.est),diag(t*s))%*%W[,N], Sigma = kron(R3,R2,R1)) , vecy[,N])
  }
  Y <- array(vecy, dim = c(t,s,b,n))

  begin_sampler <- proc.time()[3]

  # MCMC iterations
  for (i in 1:(n.burn + n.save*thin)) { #}

    rho.result = rho.update(t,s,rho,vecy,X,B.est,lam.est,W,sigma.sq)
    rho = rho.result[[1]]; inv.Sigma = rho.result[[2]]; inv.R21 = rho.result[[3]]; inv.R3 = rho.result[[4]]
    sigma.sq <- sigma.sq.update(t,s,b,n,Y,X,B.est,lam.est,W,inv.R21,inv.R3,sigma.sq,rho)
    lam.est = lam.est.update(W,inv.Sigma,B.est,X,Y,t,s,b,n)

    #W =  W.update(t,s,b,n, lam.est, inv.Sigma, B.est, Y, X, sigma.sq, W)
    eta.est = eta.est.update(X,vecy,inv.Sigma,W,lam.est, t,s,p)
    B.result =  B.est.update(X,Y,W,lam.est,inv.Sigma, inv.R1, inv.R2, t,s,p, omega, eta.est)
    B.est = B.result[[1]]; omega = B.result[[2]]; gamma.est = B.result[[3]]

    ##Keep track of MCMC output:
    if(i > n.burn & (i - n.burn)%%thin==0){
      ii = (i - n.burn)/thin

      rho.save[,ii] <- rho
      sigma.sq.save[,ii] <- sigma.sq
      lam.est.save[,ii] <- lam.est
      eta.save[,ii] <- eta.est
      B.est.save[,,ii] <- B.est
      omega.save[,,,,ii] <- omega
      gamma.est.save[,,,,ii] <- gamma.est

    }

    if(i %% ceiling((n.save)/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time(),
          " Iteration # ",i," of ", (n.save*thin + n.burn),
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - begin_sampler, " seconds #####\n",
          "##### Time per each iteration: ",(proc.time()[3] - begin_sampler)/i, " seconds #####\n"
        )
      )
    }
    #print(c("Iteration",i, round(rho,4), round(sigma.sq,4), round(lam.est,4)))

  } # end sampler

  result = list(rho.save, sigma.sq.save, lam.est.save, B.est.save, eta.save, omega.save, gamma.est.save)
  names(result) = c('rho','sigma.sq','lam.est','B.est','eta.est', 'omega', 'gamma.est')

  return(result)

} # end BSTN function




