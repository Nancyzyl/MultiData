########################################################
################## Generate data
############################################################################################
library(MASS)
library(Matrix)
library(mvtnorm)
library(ncvreg)
library(grpreg)
library(splines)
library(plyr)

##### n=250, p=50, S1
## Generate nonparametric functions
genr.np.func <- function(U, np.perc){
  n <- nrow(U)
  q <- ncol(U)
  f <- vector(mode = "list", 3)
  f[[1]] <- function(x)  - 3 * sin(pi * x - pi)       # integrate(f,-3,3)=0
  f[[2]] <- function(x)  8.85 * (x / 3 + 0.2) ^ 2 - exp( - 2 * x / 3 + 0.6)
  f[[3]] <- function(x) - 3 * sin(pi * x - pi)
  fU <- matrix(0, n, q)
  for (i in 1:q)
  {
    if (i%%3 == 1) fU[, i] <- f[[1]](U[, i])
    else if (i%%3 == 2) fU[, i] <- f[[2]](U[, i])
    else fU[, i] <- f[[3]](U[, i])
  }
  fU <- fU * np.perc
  return(fU)
}
my.index <- function(num, p){
  a <- (p - 1):1
  a <- cumsum(a)
  a <- c(0, a)
  for(j in 1:(p - 1))
  {
    if(num > a[j] & num <= a[j + 1])
    {
      r1 <- j
      r2 <- num - a[j] + j
    }
  }
  result <- c(r1, r2)
  return(result)
}

cresigma <- function(p){
  sigma <- matrix(nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i, j] <- 0.5 ^ (abs(i - j))
    }
  }
  return(sigma)
}

## get interactions of X and X or X and fU
get.interac <- function(X, fU = NULL){
  p <- ncol(X)
  if (is.null(fU))
  {
    IXX <- NULL
    for (k in 1:(p - 1)) IXX <- cbind(IXX, X[, (k + 1):p] * X[, k])
    return(IXX)
  }
  else
  {
    if (nrow(fU) != nrow(X)) stop("nrows of fU is not equal to nrows of X.")
    q <- ncol(fU)
    IXU <- NULL
    for (k in 1:q) IXU <- cbind(IXU, X * fU[, k])
    return(IXU)
  }
}
## ????????ʽ?任
changeNotes <- function(x, n, totaln, m, p){
  ind <- m * seq(0, p - 1)
  ind3 <- cumsum(n)
  ind2 <- c(0, ind3[-length(ind3)]) + 1
  xx <- c()
  if(m > 1){
    for(i in 1:m){
      xx0 <- matrix(0, n[i], m * p)
      xx0[, i + ind] <- x[(ind2[i]:ind3[i]), ]
      xx <- rbind(xx, xx0)
    }
  } else xx <- x
  return(xx)
}

##### get the subscript of the interaction
index.interac <- function(index.main, p, m){
  if (!all(index.main %in% 1:p)) stop("elements of index.main is wrong.")
  index.main <- sort(index.main)
  low <- as.vector(lower.tri(matrix(0, p, p)))
  low[low] <- 1:(p * (p - 1) / 2)
  Matr <- matrix(low, p, p, byrow = T)
  index.int <- as.vector(t(Matr[index.main, index.main]))
  index.int <- index.int[index.int != 0]
  index.int <- ((index.int - 1) * m + 1) : (index.int * m) 
  return(index.int)
}

##### Generate the data for semiparametric interaction model with continuous response #####
data.gener <- function(n, p, q, q1, np.ind, np.perc, beta.true, alpha.true){
  # n: number of subjects, it's a vector
  # p: number of covariables in parametric components X
  # q: number of covariables in nonparametric components U
  # beta.true: the coefficients of X: p-dimensional
  # alpha.true: the coefficients of interactions of X: p*(p-1)/2-dimensional
  # m: number of datasets
  # totaln: number of total subjects
  m <- length(n); totaln <- sum(n)
  X <- matrix(0, totaln, p); U <- matrix(0, totaln, q); err <- rep(0, totaln)
  X <- rmvnorm(totaln, mean = rep(0, p), sigma = 3 * diag(p)) #covariates are independent
  # X = rmvnorm(n,mean=rep(0,p),sigma=3*corr.M(p,0.5,"AR1")) #covariates are correlated
  U <- matrix(runif(totaln * q, min = -3, max = 3), nr = totaln)   #var(U)=(max-min)^2/12=3
  err <- rnorm(totaln, mean = 0, sd = 1)
  fU <- genr.np.func(U[, 1:q1, drop = FALSE], np.perc = np.perc)   # nonparametric function of U
  IX <- get.interac(X)          # interactions of X and X
  IXX <- changeNotes(IX, n, totaln, m, p*(p-1)/2)
  fU <- fU%*%diag(np.ind)
  # XX: changed X
  # X: original X
  XX <- changeNotes(X, n, totaln, m, p)
  Y <- drop(XX %*% beta.true + rowSums(fU) + IXX %*% alpha.true + err)
  data <- list(Y = Y, X = X, XX = XX, IX = IX, IXX = IXX, U = U, fU = fU)
  return(list(data = data, Y = Y, X = X, XX = XX, U = U, IX = IX, IXX = IXX, fU = fU, n = n,
              p = p, q = q, q1 = q1, beta.true = beta.true, alpha.true = alpha.true))
}
data.regener <- function(data, n, p){
  Y <- data$Y; X <- data$X;  U <- data$U; fU <- data$fU;  IX <- data$IX;
  n <- cumsum(c(0, n))
  newdata <- list(Y = Y[(n[p] + 1):n[p + 1]], XX = X[(n[p] + 1):n[p + 1], ], 
                  IXX = IX[(n[p] + 1):n[p + 1], ], U = U[(n[p] + 1):n[p + 1], ], 
                  fU = fU[(n[p] + 1):n[p + 1], ])
  return(newdata)
}

utmatr.to.vec <- function(Matr, p){
  if (nrow(Matr) != p | ncol(Matr) != p) stop("dimensions of the matrix are wrong.")
  vec <- as.vector(t(Matr))
  low <- as.vector(lower.tri(t(Matr)))
  vec <- vec[low]
  return(vec)
}

vec.to.utmatr <- function(vec, p){
  if (length(vec) != p * (p - 1)/2) stop("length of the vector is wrong.")
  low <- as.vector(lower.tri(matrix(0, p, p)))
  low[low] <- vec
  Matr <- matrix(low, p, p, byrow = T)
  return(Matr)
}
## get alpha (interactions of beta) from beta (and gamma)
beta.to.alpha2 <- function(beta, gamma = NULL){
  p <- length(beta)
  Ibeta <- NULL
  for (k in 1:(p - 1)) Ibeta <- c(Ibeta, beta[k] * beta[(k + 1):p])
  if (is.null(gamma)) alpha <- Ibeta
  else
  {
    if (length(gamma) != p * (p - 1) / 2) stop("length of gamma and beta must match.")
    else alpha <- Ibeta*gamma
  }
  return(alpha)
}
beta.to.alpha <- function(beta, m, gamma = NULL){
  beta <- matrix(beta, ncol = m, byrow = T)
  Ibeta <- NULL
  for (i in 1:ncol(beta)) Ibeta = c(Ibeta, beta.to.alpha2(beta[, i]))
  Ibeta <- matrix(Ibeta, nrow = m, byrow = T)
  Ibeta <- as.vector(Ibeta)
  if (is.null(gamma)) alpha <- Ibeta
  else {alpha <- Ibeta * gamma}
  return(alpha)
}

## get nonparametric function from B-spline basis matrix B and coeffieients phi
get.np.func <- function(B, phi, q, df){
  if (ncol(B) != q * df | length(phi) != q * df) stop("dimensions of B or phi are wrong.")
  n <- nrow(B)
  fU <- matrix(0, n, q)
  for (i in 1:q) fU[, i] <- B[, ((i - 1) * df + 1):(i * df)] %*% phi[((i - 1) * df + 1):(i * df)]
  return(fU)
}
#######################################################
# algorithm
checkConvergence <- function(beta, beta_old, eps = .001) {
  converged <- 1
  J <- length(beta)
  for (j in 1:J) {
    if (beta[j] != 0 & beta_old[j] != 0) {
      if (abs(beta[j] - beta_old[j]) > eps) {
        converged <- 0
        break
      }
    } else if (beta[j] == 0 & beta_old[j] != 0) {
      converged <- 0
      break
    } else if (beta[j] != 0 & beta_old[j] == 0) {
      converged <- 0
      break
    }
  }
  return(converged)
}

est.perf <- function(data, n, N.bs = seq(4, 12, 4), beta.true = NA, alpha.true = NA, case = 1){
  Y <- data$Y; X <- data$XX;  U <- data$U; fU <- data$fU;  IXX <- data$IXX;
  #Y<-scale(Y);X<-scale(X);IXX<-scale(IXX);U<-scale(U);fU<-scale(fU) 
  m <- length(n); p <- ncol(X)/m; q <- ncol(U); totaln <- sum(n); q1 <- ncol(fU)
  fU <- cbind(fU, matrix(0, totaln, q - q1))
  Y.obs <- Y; X.obs <- X; fU.obs <- fU; IXX.obs <- IXX
  Esti_sh_sp <- vector(mode = "list", length(N.bs))
  Esti_sp <- vector(mode = "list", length(N.bs))
  #PE = rep(0,length(N.bs))    # PE: prediction error
  RASE_sh_sp <- vector(length = length(N.bs))
  RASE_sp <- vector(length = length(N.bs))
  for (i in 1:length(N.bs)){
    print(N.bs[i])
    Esti_sh_sp[[i]] <- est.block(n, Y, X, U, IXX, p, N.bs[i], eps1 = .001)
    beta.est <- Esti_sh_sp[[i]]$beta
    alph.est <- Esti_sh_sp[[i]]$alpha
    fU.est <- Esti_sh_sp[[i]]$fU
    lam.best1 <- Esti_sh_sp[[i]]$lam1
    lam.best2 <- Esti_sh_sp[[i]]$lam2
    RASE_sh_sp[i] <- sqrt(sum(apply(fU.est - fU.obs, 2, crossprod)) / length(Y.obs))
  }
  N.bs_sh_sp <- N.bs[which.min(RASE_sh_sp)]
  print(paste0("N.bs_sh_sp=",N.bs_sh_sp))
  Esti_sh_sp <- Esti_sh_sp[[which.min(RASE_sh_sp)]]
  beta.est <- Esti_sh_sp$beta
  alpha.est <- Esti_sh_sp$alpha
  fU.est <- Esti_sh_sp$fU
  lam.best1 <- Esti_sh_sp$lam1
  lam.best2 <- Esti_sh_sp$lam2
  PE <- crossprod(Y.obs - X.obs %*% beta.est - IXX.obs %*% alpha.est - rowSums(fU.est)) / length(Y.obs)
  RASE <- min(RASE_sh_sp)
  beta.adjust <- alpha.adjust <- NULL
  print(case)
  if (case != 1){
    for (i in 1:p){
      beta.adjust <- c(beta.adjust, rep(beta.est[i], 3)) 
    }
    for(i in 1:(p * (p - 1) / 2)){
      alpha.adjust <- c(alpha.adjust, rep(alpha.est[i], 3))
    }
    MSE <- crossprod(beta.adjust - beta.true)/length(beta.true)
    MSE2 <- crossprod(alpha.adjust - alpha.true)/length(alpha.true)
  } else {
    MSE <- crossprod(beta.est - beta.true)/length(beta.true)
    MSE2 <- crossprod(alpha.est - alpha.true)/length(alpha.true)
  }
  return(list(beta = beta.est, alpha = alpha.est, fU.est = fU.est, PE = PE, MSE = MSE, MSE2 = MSE2,
              RASE = RASE, N.bs = N.bs_sh_sp, lam.best1 = lam.best1, lam.best2 = lam.best2))
  #return(list(beta=beta.est,alpha=alpha.est,PE=PE,MSE=MSE,RASE=RASE,Fit=Fit,Sel=Sel.times))
}
est.block <- function(n, Y, X, U, IXX, p, df_bs = 4, eps1 = .001){
  totaln <- length(Y); m <- length(n); q <- ncol(U); pp <- p * m
  B <- matrix(0, totaln, q * df_bs)   # B-spline basis matrix, df_bs is number of the B-spline basis
  for (i in 1:q) 
    B[, ((i - 1) * df_bs + 1):(i * df_bs)] <- sqrt(3) * scale(as.matrix(ns(U[, i], df = df_bs)))  
  #ns() Generate the B-spline basis matrix for a natural cubic spline.
  #Q = chol(solve(V))  # Q:upper triangular matrix, t(Q)%*%Q=solve(V)
  XB <- cbind(X, IXX, B) 
  #OLS estimates with Generalized Inverse or ridge regression as initial estimates
  sol <- lm.ridge(Y ~ XB, lambda = seq(0, 1, 0.001))
  #print(sol$lambda[which.min(sol$GCV)])
  ini_ols <- as.vector(sol$coef[, which.min(sol$GCV)])
  phi_old <- ini_ols[-(1:(m * (p + p * (p - 1) / 2)))]   # coefficients of the B-spline basis matrix: B
  fU_old <- get.np.func(B, phi_old, q, df_bs)
  #fU = genr.np.func(U[,1:q1,drop=FALSE],np.perc=np.perc)
  XB <- cbind(X, IXX) 
  YB <- Y - rowSums(fU_old)
  sol <- lm.ridge(YB ~ XB, lambda = seq(0, 1, 0.001))
  #sol = lm.ridge(Y~XB,lambda = seq(0,1,0.001))
  #print(sol$lambda[which.min(sol$GCV)])
  ini_ols <- as.vector(sol$coef[, which.min(sol$GCV)])
  beta_old <- ini_ols[1:pp]
  #print(beta_old)
  gamma_old <- ini_ols[(pp + 1):(m * (p + p * (p - 1) / 2))] / beta.to.alpha(beta_old, m)
  
  beta0 <- beta_old
  gamma0 <- gamma_old
  #write.table(beta0,paste0(df_bs,"beta0.txt"))
  #write.table(gamma0,paste0(df_bs,"gamma0.txt"))
  t <- 1      # t indexes the inner iterations
  repeat
  {
    gamma0 <- gamma_old
    beta0 <- beta_old
    # update gamma and eta from beta and phi
    #fU_old <- get.np.func(B, phi_old, q, df_bs)
    Y_tilde <- Y - X %*% beta_old - rowSums(fU_old) 
    Ibeta <- beta.to.alpha(beta_old, m)
    #print(Ibeta)
    XX_tilde <- IXX %*% diag(Ibeta) 
    # Note that maybe all elements in one column of XXU_tilde are zero
    index <- 1 - apply(XX_tilde, 2, function(x) all(x == 0))
    ind <- matrix(index, ncol = m, byrow = T)
    ind[apply(ind, 1, sum)!=0, ] <- rep(1, m)
    index <- as.vector(t(ind))
    if (any(index == 0)){
      if (all(index == 0)) {
        gamma_new <- gamma_old
      } else {
        index0 <- which(index == 0)
        XX_tilde <- as.matrix(XX_tilde[, -index0])
        #a=ncol(XX_tilde)
        #if(a<2) print(a)
        gamma_new <- rep(0, length(gamma_old))
        pq <- ncol(XX_tilde) / m
        groups <- NULL
        for(i in 1:pq){
          groups <- c(groups, rep(i, m))
        }
        alasso_res <- grpreg(XX_tilde, Y_tilde, group = groups, penalty = "grLasso", family = "gaussian",
                             lambda.min = 1e-1, nlambda = 500)
        #alasso_res = adaptive.lasso_ncvreg(XX_tilde, Y_tilde, weight=gamma0[-index0], lam.vect=lambda)
        gamma_new[-index0] <- drop(select(alasso_res, criterion = "BIC")$beta)[-1]
        lambda_gamma <- drop(select(alasso_res, criterion = "BIC")$lambda)
        #gamma_matrix = matrix(rep(gamma_old, length(lambda_gamma_vect)), nrow=length(gamma_old))
        #gamma_matrix[-index0, ] = alasso_res$coeff 
      }
    } else { 
      groups <- NULL
      for(i in 1:(p * (p - 1) / 2)){
        groups <- c(groups, rep(i, m))
      }
      alasso_res <- grpreg(XX_tilde, Y_tilde, group = groups, penalty = "grLasso", family = "gaussian",
                           lambda.min = 1e-1, nlambda = 500)
      #alasso_res = adaptive.lasso_ncvreg(XX_tilde, Y_tilde, weight=gamma0, lam.vect=lambda)
      gamma_new <- drop(select(alasso_res, criterion = "BIC")$beta)[-1]
      lambda_gamma <- drop(select(alasso_res, criterion="BIC")$lambda)
      #gamma_matrix = alasso_res$coeff        
    }
    ## update phi (or fU) from beta,  gamma 
    alpha_new <- beta.to.alpha(beta_old, m, gamma_new)    
    Y_tilde <- Y - X %*% beta_old - IXX %*% alpha_new 
    XXU_tilde <- B
    if(ncol(B) >= nrow(B)){
      sol <- lm.ridge(Y_tilde ~ XXU_tilde, lambda = seq(0, 0.1, 0.0001))
      #print(which.min(sol$GCV))
      phi_new <- as.vector(sol$coef[, which.min(sol$GCV)])
    } else {
      sol <- lm(Y_tilde~XXU_tilde-1)
      phi_new <- sol$coef
    }
    #phi_new<-phi_old  
    fU_new <- get.np.func(B, phi_new, q, df_bs)
    
    ## update beta from gamma,  phi (or fU) and eta
    beta_new <- beta_old
    # M_eta <- matrix(eta_new, p, q)           # p*q matrix of eta
    beta_new <- beta_old
    Y_tilde <- Y - IXX %*% alpha_new - rowSums(fU_new) 
    XXU_tilde <- X
    groups <- NULL
    for(i in 1:p){
      groups <- c(groups, rep(i, m))
    }
    block_res <- grpreg(XXU_tilde, Y_tilde, group = groups, penalty = "grLasso", family = "gaussian",
                        lambda.min = 1e-1, nlambda = 500)
    beta_new <- drop(select(block_res, criterion = "BIC")$beta)[-1]
    lambda_beta <- drop(select(block_res, criterion = "BIC")$lambda)
    #print(beta_matrix)
    
    if ((crossprod(beta_new - beta_old) <= eps1 &
         crossprod(gamma_new - gamma_old) <= eps1) | t > 30)
      break
    t <- t + 1
    beta_old <- beta_new
    #print(beta_old)
    gamma_old <- gamma_new
    phi_old <- phi_new
    fU_old <- fU_new
  }
  print(paste0("t=", t))
  #print(paste0("lambda_gamma=",lambda_gamma))
  #print(paste0("lambda_beta=",lambda_beta))
  #write.table(lambda_gamma,paste0(df_bs,"lambda_gamma.txt"))
  #write.table(lambda_beta,paste0(df_bs,"lambda_beta.txt"))
  beta <- beta_new
  gamma <- gamma_new
  phi <- phi_new
  fU <- get.np.func(B, phi, q, df_bs)  
  alpha <- beta.to.alpha(beta, m, gamma) 
  #write.table(beta,paste0(df_bs,"beta.est.txt"))
  #write.table(gamma,paste0(df_bs,"gamma.est.txt"))
  #write.table(alpha,paste0(df_bs,"alpha.est.txt"))
  #write.table(fU,paste0(df_bs,"fU.txt"))
  return(list(beta = beta, alpha = alpha, fU = fU, B = B, phi = phi, 
              lam1 = lambda_beta, lam2 = lambda_gamma))
}

PerformEva <- function(true, est){
  tp <- sum((colSums(true) == 0) * (colSums(est) == 0)) / sum((colSums(true) == 0))
  tn <- sum((colSums(true) != 0) * (colSums(est) != 0)) / sum((colSums(true) != 0))
  fp <- sum((colSums(true) == 0) * (colSums(est) != 0)) / sum((colSums(true) == 0))
  fn <- sum((colSums(true) != 0) * (colSums(est) == 0)) / sum((colSums(true) != 0))
  ans <- c(tp, tn, fp, fn)
  return(ans)
}
set.seed(525)
run <- 2
n <- c(100, 80, 70)
p <- 10
q <- 4
q1 <- 2
m <- 3
beta.true <- rep(0, m * p)
beta.true[1:15] <- c(2, 0.5, 2, 1.5, 1, 1, 1, 1, 0.8, 1, 1.5, 1, 0.8, 2, 1.5)
alpha.true <- rep(0, m * p * (p - 1) / 2)
alpha.true[index.interac(1:2, p, m)] <- c(2, 2, 1.5)
alpha.true[index.interac(c(1, 3), p, m)] <- c(1.5, 1, 1)
alpha.true[index.interac(2:3, p, m)] <- c(1, 0.5, 0.5)
alpha.true[index.interac(4:5, p, m)] <- c(0.5, 1.5, 2)

dir.create(paste("./sim1", run, "_p", p, sep = ''), showWarnings = FALSE)
setwd(paste("./sim1", run, "_p", p, sep = ''))

sim.run <- function(run, n, p, q, q1, beta.true, alpha.true){
  #make arrays to save the results
  m <- length(n); pp <- p * m
  beta.trues <- matrix(beta.true, ncol = p)
  alpha.trues <- matrix(alpha.true, ncol = p * (p - 1) / 2)
  beta.list <- beta.temp <- array(0, dim = c(p * m, run), dimnames = list(paste0("X", 1:(p * m), NULL)))
  beta.list3 <- array(0, dim = c(p, run), dimnames = list(paste0("X", 1:p, NULL)))
  beta.lists <- array(0, dim = c(p, run, m))
  fU.list <- fU.list3 <- fU.lists <- array(0, dim = c(sum(n), run))
  name.X <- name.XX <- NULL
  for (k in 1:(p - 1)) name.X <- c(name.X, (paste0(paste0("X", k), paste0("X", (k + 1):p))))
  for (k in 1:length(name.X)) name.XX <- c(name.XX, rep(name.X[k], m))
  alpha.list <- alpha.temp <- array(0, dim = c(m * p * (p - 1) / 2, run), dimnames = list(name.XX, NULL))
  alpha.list3 <- array(0, dim = c(p * (p - 1) / 2, run), dimnames = list(name.X, NULL))
  alpha.lists <- array(0, dim = c(p * (p - 1) / 2, run, m))
  ##################################################
  Err.list <- Err.temp <- matrix(0, run, 4, dimnames = list(NULL, c("PE", "MSE.beta", "MSE.alpha", "RASE")))
  Err.list3 <- matrix(0, run, 4, dimnames = list(NULL, c("PE", "MSE.beta", "MSE.alpha", "RASE")))
  Err.lists <- array(0, dim = c(run, 4, m), dimnames = list(NULL, c("PE", "MSE.beta", "MSE.alpha", "RASE")))
  perf.beta.m3 <- perf.alpha.m3 <- matrix(0, nrow = 4, ncol = run)
  perf.beta.m2 <- perf.alpha.m2 <- matrix(0, nrow = 4, ncol = run)
  perf.beta.m1 <- perf.alpha.m1 <- matrix(0, nrow = 4, ncol = run)
  N.bs <- N.bs3 <- rep(0, run)
  N.bss <- matrix(0, nrow = m, ncol = run)
  #repeat simulations for run times 
  for (i in 1:run){
    print(paste("run", i, " "))
    my <- data.gener(n, p, q = 4, q1 = 2, np.ind = c(1, 1), np.perc = 1, beta.true, alpha.true)
    result2 <- est.perf(my$data, n, N.bs = seq(4, 24, 4), beta.true, alpha.true) #M2
    perf.beta.m2[, i] <- PerformEva(beta.trues, matrix(result2$beta, ncol = p))
    perf.alpha.m2[, i] <- PerformEva(alpha.trues, matrix(result2$alpha, ncol = p * (p - 1) / 2))
    beta.list[, i] <- result2$beta
    alpha.list[, i] <- result2$alpha
    fU.list[, i] <- rowSums(result2$fU)
    Err.list[i, ] <- c(result2$PE, result2$MSE, result2$MSE2, result2$RASE)
    N.bs[i] <- result2$N.bs
    sepdata <- list(Y = my$Y, X = my$X, IX = my$IX, U = my$U, fU = my$fU)
    sepdata3 <- list(Y = my$Y, XX = my$X, IXX = my$IX, U = my$U, fU = my$fU)
    result3 <- est.perf(sepdata3, sum(n), N.bs = seq(4, 24, 4), beta.true = beta.true, alpha.true = alpha.true,
                        case = 3) #M3
    perf.beta.m3[, i] <- PerformEva(beta.trues, matrix(result3$beta, ncol = p))
    perf.alpha.m3[, i] <- PerformEva(alpha.trues, matrix(result3$alpha, ncol = p * (p - 1) / 2))
    beta.list3[, i] <- result3$beta
    alpha.list3[, i] <- result3$alpha
    fU.list3[, i] <- rowSums(result3$fU)
    Err.list3[i, ] <- c(result3$PE, result3$MSE, result3$MSE2, result3$RASE)
    N.bs3[i] <- result3$N.bs
    fU.temp <- NULL
    for(t in 1:length(n)){    #M1
      print(paste("lay", t, " "))
      newdata <- data.regener(sepdata, n, t)
      result <- est.perf(newdata, n[t], N.bs = seq(4, 24, 4), beta.trues[t, ], alpha.trues[t, ])
      beta.lists[, i, t] <- result$beta
      alpha.lists[, i, t] <- result$alpha
      fU.temp <- c(fU.temp, rowSums(result$fU))
      Err.lists[i, , t] <- c(result$PE, result$MSE, result$MSE2, result$RASE)
      N.bss[t, i] <- result$N.bs
    }
    beta.temp[, i] <- as.vector(t(beta.lists[, i, ]))
    alpha.temp[, i] <- as.vector(t(alpha.lists[, i, ]))
    fU.lists[, i] <- fU.temp
    Err.temp[i, ] <- apply(Err.lists[i, , ], 1, mean)
    perf.beta.m1[, i] <- PerformEva(beta.trues, matrix(beta.temp[, i], ncol = p))
    perf.alpha.m1[, i] <- PerformEva(alpha.trues, matrix(alpha.temp[, i], ncol = p * (p - 1) / 2))
  }
  #M2
  beta <- cbind(apply(beta.list, 1, mean), apply(beta.list, 1, sd))
  alpha <- cbind(apply(alpha.list, 1, mean), apply(alpha.list, 1, sd))
  fU <- cbind(apply(fU.list, 1, mean), apply(fU.list, 1, sd))
  Err <- rbind(colMeans(Err.list), apply(Err.list, 2, sd))
  perf.beta2 <- apply(perf.beta.m2, 1, mean) 
  perf.alpha2 <- apply(perf.alpha.m2, 1, mean)
  #M1
  betas <- cbind(apply(beta.temp, 1, mean), apply(beta.temp, 1, sd))
  alphas <- cbind(apply(alpha.temp, 1, mean), apply(alpha.temp, 1, sd))
  fUs <- cbind(apply(fU.lists, 1, mean), apply(fU.lists, 1, sd))
  Errs <- rbind(colMeans(Err.temp), apply(Err.temp, 2, sd))
  perf.beta1 <- apply(perf.beta.m1, 1, mean)
  perf.alpha1 <- apply(perf.alpha.m1, 1, mean)
  #M3
  beta3 <- cbind(apply(beta.list3, 1, mean), apply(beta.list3, 1, sd))
  alpha3 <- cbind(apply(alpha.list3, 1, mean), apply(alpha.list3, 1, sd))
  fU3 <- cbind(apply(fU.list3, 1,mean), apply(fU.list3, 1, sd))
  Err3 <- rbind(colMeans(Err.list3), apply(Err.list3, 2, sd))
  perf.beta3 <- apply(perf.beta.m3, 1, mean) 
  perf.alpha3 <- apply(perf.alpha.m3, 1, mean)
  #
  perf.beta <- rbind(perf.beta1, perf.beta2, perf.beta3)
  perf.alpha <- rbind(perf.alpha1, perf.alpha2, perf.alpha3)
  return(list(beta.list = beta.list, beta = beta, alpha.list = alpha.list, alpha = alpha,
              fU = fU, Err.list = Err.list, Err = Err, N.bs = N.bs, perf.beta = perf.beta,
              beta.temp = beta.temp, betas = betas, alpha.temp = alpha.temp, alphas = alphas,
              Err.temp = Err.temp, Errs = Errs, N.bss = N.bss, perf.alpha = perf.alpha,
              beta.list3 = beta.list3, beta3 = beta3, alpha.list3 = alpha.list3, alpha3 = alpha3,
              fU3 = fU3, Err.list3 = Err.list3, Err3 = Err3, N.bs3 = N.bs3))
}

est <- sim.run(run, n, p, q, q1, beta.true, alpha.true)

est.beta <- est$beta.list
est.beta[est.beta != 0] = 1
timb <- 0
for(i in 1:run){
  temp <- matrix(est.beta[, i], ncol=p)
  res_sum <- apply(temp, 2, sum)
  res_sum[res_sum != 0] = 1
  timb <- timb + res_sum
}
est.alpha <- est$alpha.list
est.alpha[est.alpha != 0] = 1
tima <- 0
for(i in 1:run){
  temp <- matrix(est.alpha[, i], nrow=m)
  res_sum <- apply(temp, 2, sum)
  res_sum[res_sum != 0] = 1
  tima <- tima + res_sum
}


est.betas <- est$beta.temp
est.betas[est.betas != 0] = 1
timbs <- 0
for(i in 1:run){
  temp <- matrix(est.betas[, i], ncol=p)
  res_sum <- apply(temp, 2, sum)
  res_sum[res_sum != 0] = 1
  timbs <- timbs + res_sum
}
est.alphas <- est$alpha.temp
est.alphas[est.alphas != 0] = 1
timas <- 0
for(i in 1:run){
  temp <- matrix(est.alphas[, i], nrow=m)
  res_sum <- apply(temp, 2, sum)
  res_sum[res_sum != 0] = 1
  timas <- timas + res_sum
}

est.beta3 <- est$beta.list3
est.beta3[est.beta3 != 0] = 1
timb3 <- apply(est.beta3, 1, sum)
est.alpha3 <- est$alpha.list3
est.alpha3[est.alpha3 != 0] = 1
tima3 <- apply(est.alpha3, 1, sum)

write.csv(rbind(timbs, timb, timb3), "timb.csv")
write.csv(rbind(timas, tima, tima3), "tima.csv")
write.csv(round(rbind.fill(data.frame(est$Errs), data.frame(est$Err), data.frame(est$Err3)), 3), "err.csv")
write.csv(round(est$perf.beta, 3), 'perfBeta.csv')
write.csv(round(est$perf.alpha, 3), 'perfAlpha.csv')
