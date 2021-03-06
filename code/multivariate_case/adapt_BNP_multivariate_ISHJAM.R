Rcpp::sourceCpp("code/multivariate_case/adapt_BNP_multivariate_helpers.cpp")
library(LaplacesDemon)
library(MCMCpack)

# CAP model  --  ISHWARAN JAMES

Update.pi <- function(alphabeta, aDir, J) {
  alpha        <- alphabeta[, 1]
  ns           <- table(factor(alpha, levels = 0:J))
  nj           <- c(ns[-1], ns[1])
  aDirtilde    <- nj + aDir
  return(MCMCpack::rdirichlet(1, aDirtilde))
}

################################ robust estimates extractor

Omegatilde.maker <- function(pi, omega, J) {
  return(c(pi[-(J + 1)], exp(log(pi[J + 1]) + log(omega))))
}

log_omegatilde.maker <- function(pi, omega, J) {
  return(c(log(pi[-(J + 1)]), log(pi[J + 1]) + log(omega)))
}


StartingEst <- function(X, categ, raw_MCD, h_MCD, J) {
  
  if (h_MCD == 1) {
    #standard non-robust methods are computed
    estimates_from_train <-
      mclust::covw(X = X, Z = mclust::unmap(categ))
    xbar_j               <- estimates_from_train$mean
    S2_j                 <- estimates_from_train$S
  } else {
    robust_estimates_from_train <-
      lapply(unique(sort(categ)), function(group)
        rrcov::CovMcd(x = X[categ == group, ], alpha = h_MCD))
    
    if (raw_MCD) {
      xbar_j <-
        sapply(1:J, function(g)
          robust_estimates_from_train[[g]]$center)
      S2_j   <-
        array(
          sapply(1:J, function(g)
            robust_estimates_from_train[[g]]$cov) / nrow(X) * (nrow(X) - 1),
          dim = c(p, p, J)
        )
    } else{
      xbar_j <-
        sapply(1:J, function(g)
          robust_estimates_from_train[[g]]$raw.center)
      S2_j   <-
        array(
          sapply(1:J, function(g)
            robust_estimates_from_train[[g]]$raw.cov) / nrow(X) * (nrow(X) -
                                                                     1),
          dim = c(p, p, J)
        )
    }
  }
  return(list(xbar_j, S2_j))
}

# alphaneta_to_zeta <- function(AB, n, J) {
#   Z           <- AB[, 1]
#   Z[AB[, 2] > 0] <- AB[AB[, 2] > 0, 2] + J
#   return(Z)
# }

# zeta_to_alphabeta <- function(Z, n, J) {
#   AB   <- matrix(0, n, 2)
#   ind1 <- which(Z <= J)
#   AB[ind1, 1]  <- Z[ind1]
#   AB[-ind1, 2] <- Z[-ind1] - J
#   return(AB)
# }


###################################################################################

MCMC_adaptBNP_multivariate_ISHJAMES <- function(Y,
                                             X,
                                             categ,
                                             prior,
                                             L,
                                             h_MCD = .75,
                                             # percentage of untrimmed observations for computing robust prior mean vector and scatter from train set
                                             raw_MCD = F,
                                             NSIM,
                                             thinning,
                                             burn_in,
                                             verbose = T,
                                             fixed_alphaDP = F,
                                             light = F,
                                             learning_type = c("transductive")){#, "inductive")) {
  #####################################################################################
  Y <-
    data.matrix(Y) # for avoiding problems in dealing with data.frame objects
  X <- data.matrix(X)
  #####################################################################################
  n      <- nrow(Y) # sample size test data
  p      <- ncol(Y)
  J      <- length(unique(categ)) # # known categories
  #####################################################################################
  SERob  <-
    StartingEst(
      X = X,
      categ = categ,
      raw_MCD = raw_MCD,
      h_MCD = h_MCD,
      J=J
    )
  xbar_j <- SERob[[1]]
  S2_j   <- SERob[[2]]
  #####################################################################################
  ### Hyperparameters
  aDir     <- prior$aDir
  aDP      <- prior$aDP
  m_H      <- prior$m_H
  k_H      <- prior$k_H
  v_H      <- prior$v_H
  S_H      <- prior$S_H
  k_g      <- prior$k_g
  v_g      <- prior$v_g
  a_alpha  <- prior$a_alphaDP
  b_alpha  <- prior$b_alphaDP
  #####################################################################################
  ### Containers
  Alpha_Beta   <- array(NA, dim = c(n, 2, NSIM))
  MU_new       <- vector("list", length = NSIM)
  SIGMA_new    <- vector("list", length = NSIM)
  
  if (learning_type %in% c("transductive", "training_stochastic")) {
    MU_train     <- array(NA, dim = c(J, p, NSIM))
    SIGMA_train  <- array(NA, dim = c(p, p, J, NSIM))
  }
  
  PiDir        <- matrix(NA, NSIM, J + 1) # J+1-th ? pi0
  Omega        <- vector("list", length = NSIM)
  AlphaDP      <- numeric(NSIM)
  U_Slice      <- matrix(NA, n, NSIM)
  #####################################################################################
  # inizialization
  pidir         <- c(rdirichlet(1, aDir))
  alphabeta     <- matrix(0, n, 2)
  alphabeta[, 1] <- sample(0:J, n, replace = T)
  u             <- rbeta(n = L,
                         shape1 = 1,
                         shape2 = aDP)
  omega         <- StickBreaker_cpp(V = u)
  alphabeta[alphabeta[, 1] == 0, 2] <-
    as.numeric(factor(sample(
      1:L,
      size =
        sum(alphabeta[, 1] ==
              0),
      replace = T,
      prob = omega
    )))
  Sigma         <- replicate(L, riwish(v = v_H, S = S_H))
  mu            <-
    sapply(1:L, function(ind)
      mvtnorm::rmvnorm(
        n = 1,
        mean = m_H,
        sigma = Sigma[, , ind] / k_H
      ))
  mu.train      <- xbar_j
  Sigma.train   <- S2_j

  if (verbose) {
    cat("MCMC progress:\n")
    flush.console()
    pbar <- txtProgressBar(min = 1,
                           max = NSIM * thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  
  
  for (sim in 1:(NSIM * thinning + burn_in)) {
    ##################################################################################
    pidir      <-
      c(Update.pi(
        alphabeta = alphabeta,
        aDir = aDir,
        J = J
      ))
    ###################################################################################
    u          <-
      UPD_Sticks_Beta_cpp(AB = alphabeta,
                          L = L,
                          alphaDP = aDP)
    if (L == 1) {
      omega <- 1
    } else{
      omega <- StickBreaker_cpp(u)
    }
    
    log_omegatilde <- log_omegatilde.maker(pidir, omega, J)
    ###################################################################################
    ####################################################################################
    TH         <- Update_Theta_cpp(
      Y = Y,
      AB = alphabeta,
      m_H = m_H,
      k_H = k_H,
      v_H = v_H,
      S_H = S_H,
      L = L,
      p = p
    )
    mu         <- TH[[1]]
    Sigma      <- TH[[2]]
    #####################################################################################
    if (learning_type %in% c("transductive", "training_stochastic")) {
      THETA_g    <- Update_Theta_cpp_TRAIN(
        Y = Y,
        AB = alphabeta,
        MU_g = xbar_j,
        k_g = k_g,
        v_g = v_g,
        SIGMA_g =  S2_j * (v_g - p - 1),
        G = J,
        p = p
      )
      
      mu.train    <- THETA_g[[1]]
      Sigma.train <- THETA_g[[2]]
    }
    ############################################################################################
    alphabeta  <- Upd_alphabeta_cpp(Y=Y,
                                    mu =  mu, Sigma = Sigma,
                                    pidir =  pidir,omega =  omega,
                                    J =  J,
                                    x_bar =  mu.train, 
                                    S2_j =  Sigma.train,
                                    L=L,
                                    n=n, poss_lab = 1:(L+J))

    index <- alphabeta[alphabeta[,2]>0,2]
    u.ind <- unique(sort(index))
    mu    <- mu[,u.ind]
    Sigma <- Sigma[,,u.ind]
    alphabeta[alphabeta[,2]>0,2] <- as.numeric(factor(alphabeta[alphabeta[,2]>0,2]))
    #########################################################################################################
    if(fixed_alphaDP  & sim == 1){
      aDP <- aDP
    }else if (fixed_alphaDP == F){
      aDP  <- rgamma(1, shape =  a_alpha + (L-1),  rate =    b_alpha -  sum(log(1-u[-L])) )
    }
    ############################################################################################
    
    if (sim > burn_in &&
        ((sim - burn_in) %% thinning == 0)) {
      rr                 <- floor((sim - burn_in) / thinning)
      
      PiDir[rr, ]         <- pidir
      Alpha_Beta[, , rr]   <- alphabeta
      Omega[[rr]]        <- omega
      MU_new[[rr]]       <- t(mu)
      SIGMA_new[[rr]]    <- Sigma
      if (learning_type %in% c("transductive", "training_stochastic")) {
        MU_train[, , rr]     <- t(mu.train)
        SIGMA_train[, , , rr] <- Sigma.train
      }
      AlphaDP[rr]        <- aDP
      }
    ################################################################
    ################################################################
    
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }
  
  # If training fixed I only report the output from robust information extraction
  if (learning_type %in% c("inductive", "training_fixed")) {
    MU_train     <- t(mu.train)
    SIGMA_train <- Sigma.train
  }
  
  if (light) {
    out <-  list(
      P         = PiDir,
      AB        = Alpha_Beta,
      O         = Omega,
      Mu        = MU_new,
      Sigma     = SIGMA_new,
      Mu.train  = MU_train,
      Sig.train = SIGMA_train,
      x_b       = xbar_j,
      s_b       = S2_j,
      alphaDP   = AlphaDP
    )
    
  } else{
    out <-  list(
      P         = PiDir,
      AB        = Alpha_Beta,
      O         = Omega,
      Mu        = MU_new,
      Sigma     = SIGMA_new,
      Mu.train  = MU_train,
      Sig.train = SIGMA_train,
      x_b       = xbar_j,
      s_b       = S2_j,
      alphaDP   = AlphaDP
    )
  }
  
  
  return(out)
}
