library(LaplacesDemon)
library(MCMCpack)
#Rcpp::sourceCpp(here::here("code/functional_extension_code/SLICE_functional_helpers.cpp"))

Update.pi <- function(alphabeta,aDir,J){
  alpha        <- alphabeta[,1]
  ns           <- table(factor(alpha,levels = 0:J))
  # tabulate??
  nj           <- c(ns[-1],ns[1])
  aDirtilde    <- nj+aDir
  return(rdirichlet(1,aDirtilde))
}

a.AND.b_sigmaT_train <- function(sigma.est.train, # matrix TxG
                                 vg,t,G){
  ab      <- array(NA,c(t,G,2))
  ab[,,1] <- 2 + (sigma.est.train)^2/vg
  ab[,,2] <- ((sigma.est.train)^2/vg+1)*sigma.est.train
  
  return(ab)
}



Omegatilde.maker <- function(pi,omega,J){
  
  return(c(pi[-(J+1)],pi[J+1]*omega))
  
}

log_omegatilde.maker <- function(pi, omega, J) {
  return(c(log(pi[-(J + 1)]), log(pi[J + 1]) + log(omega)))
  
}

alphaneta_to_zeta <- function(AB,n,J){
  Z           <- AB[,1]
  Z[AB[,2]>0] <- AB[AB[,2]>0,2] + J
  return(Z)
}

zeta_to_alphabeta <- function(Z,n,J){
  AB   <- matrix(0,n,2)
  ind1 <- which(Z<=J)
  AB[ind1,1]  <- Z[ind1]
  AB[-ind1,2] <- Z[-ind1] - J
  return(AB)
}



MCMCap.funct.SLICE <- function(Y,
                               prior,
                               L,
                               # L numero possibili gruppi
                               nsim,
                               thinning,
                               burn_in,
                               fixed_alphaDP = F,
                               kappa = .5,
                               verbose = 1,
                               learning_type = c("transductive", "inductive")) {
  
  #####################################################################
  ############
  basis       <- prior$basis
  beta.train  <- prior$beta.train
  sigma.train0<- prior$sigma.train # now this is a matrix of functions
  effe.train0 <- basis%*%beta.train
  vg          <- prior$vg
  KappaG      <- prior$KappaG   
  n           <- ncol(Y)
  t           <- nrow(Y)
  J           <- ncol(beta.train)
  K           <- nrow(beta.train) # numero basi
  ############
  R           <- rowSums(basis^2)
  alphaD      <- rep(1,J)
  OneOverR    <- 1/R
  # HYPERPRIOR PARAM for STOCHASTIC TRAINING
  A_B_IG_Train_Stochastic <- a.AND.b_sigmaT_train(sigma.est.train = sigma.train0,vg = vg,
                                                  t = t,G = J)
  ###############################################
  # Containers
  EFFE.TEST   <- vector("list",length = nsim)#array( NA,dim = c(t,L,nsim))
  SIGMA.TEST  <- vector("list",length = nsim)
  ###
  PROB        <- matrix(NA,nsim, J+1)
  OMEGA       <- vector("list",length = nsim)
  TAU         <- vector("list",length = nsim)
  AK          <- vector("list",length = nsim)
  ALPHA.BETA  <- array( NA,dim = c(n,2,nsim))
  AlphaDP     <- numeric(nsim)
  USLICE      <- matrix(NA,n,nsim)
  ###
  if (learning_type %in% c("transductive", "training_stochastic")) {
    EFFE.TRAIN  <- array(NA, dim = c(t, J, nsim))
    SIGMA.TRAIN <- array(NA, dim = c(t, J, nsim))
  }
  ###
  ##############################################
  ### Hyperparameters
  aDir     <- prior$aDir
  aDP      <- prior$aDP
  a_H      <- prior$a_H  # prior varianza gruppi nuovi 
  b_H      <- prior$b_H      
  a_tau    <- prior$a_tau  # varianza dei betini
  b_tau    <- prior$b_tau  # varianza dei betini
  a_alpha  <- prior$a_alphaDP
  b_alpha  <- prior$b_alphaDP
  s_tau    <- prior$s_tau
  # inizialization
  pidir         <- c(rdirichlet(1,aDir))
  alphabeta     <- matrix(NA,n,2)   
  alphabeta[,1] <- sample(0:J,n,replace = T)
  alphabeta[,2] <- ifelse(alphabeta[,1]>0,0,sample(L))
  u             <- rbeta(n = L, shape1 = 1, shape2 = aDP)
  omega         <- StickBreaker_cpp(V = u)
  effe.test     <- matrix(NA,t,L)
  sigma.test    <- matrix(NA,t,L)
  effe.train    <- effe.train0 #matrix(NA,t,J)
  sigma.train   <- sigma.train0 #matrix(NA,t,J)
  tau           <- numeric(L)
  ak            <- numeric(L)
  for(k in 1:L){ 
    sigma.test[,k] <- rinvgamma(t,a_H,b_H)
    tau[k]         <- rinvgamma(1,a_tau,b_tau)
    ak[k]          <- rnorm(1,0,sqrt(s_tau))
    effe.test[,k]  <- rnorm(t,0,sd=sqrt(tau[k]*R))
  }
  ZETA <- alphaneta_to_zeta(alphabeta,n = n,J = J)
  #####################################################################################
  #g <- function(j) (1-kappa)*kappa^(j-1) * (j>0)
  # magari 
  #g2 <- function(j) ifelse(j<=J, (1-kappa)/J,  kappa*(1-kappa)^(j-1)  )
  ######################################################################################
  #oppure
  # g2 <- function(j,J) ifelse(j<=J, (1-kappa)/(J+1), 
  #                            (J*kappa+1)/(J+1) * (1-kappa)/(J*kappa+1) * ((J*kappa+kappa)/(J*kappa+1)) ^ (j-J-1)  )
  g2 <- function(j,J) ifelse(j<=J, (1-kappa)/(J+1), 
                             (1-kappa)/(J+1) *
                               (((J+1)*kappa)/(J*kappa+1)) ^ (j-J-1)  )
  
  # oppure, primo stick delle novelty can be higher  delle vecchie
  # cat("look at me before playing around!")
  # g2 <- function(j,J,theta,kappa=.9) ifelse(j<=J, (theta)/(J+1), 
  #                                     (J*theta+1)/(J+1) * (1-kappa) * (kappa) ^ (j-J-1)  )
  # # Best configuaration so far, primo stick novelty minore delle old, ma progressione a zero molto migliore
  #plot(g2(1:20,J = 4,theta = .90,kappa = .85))
  
  if (verbose) {
    cat("MCMC progress:\n")
    flush.console()
    pbar <- txtProgressBar(min = 1,
                           max = nsim*thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  
  for(sim in 1:(nsim*thinning + burn_in)){
    
    ############################################################################
    Ui    <- runif(n, 0, 100*g2(ZETA,J = J))/100
    #L.z   <- 1 + floor((log(Ui) - log(1 - kappa)) / log(kappa))
    L.z   <- J + 1 + floor(
      (log(Ui) - (log( (J * kappa + 1) / (J + 1) ) + log( (1 - kappa) / (J * kappa + 1) )))/log(((J * kappa + kappa) / (J * kappa + 1)) ))
    if(length(L.z)==0){ L.z <- 1 }
    J.z  <- max(alphabeta[,2])
    L    <- max(c(L.z, J.z))
    xi   <- g2(1:L, J = J)
    PL   <- c(1:(J+L))
    #############################################################################
    pidir   <- c(Update.pi(alphabeta = alphabeta, aDir = aDir, J = J))
    #######################################################################
    u       <- UPD_Sticks_Beta_cpp(AB = alphabeta,L = L,alphaDP = aDP)
    if(L==1){
      omega <- 1  
    }else{
      omega <- StickBreaker_cpp(u)
    }
    # omegatilde <- Omegatilde.maker(pidir,omega,J)
    log_omegatilde <- log_omegatilde.maker(pidir, omega, J)
    #######################################################################
    tau        <- Update_tau_l(a_tau = a_tau,b_tau = b_tau,
                               effe_test = as.matrix(effe.test),
                               L = L,
                               t = t,
                               OOR = OneOverR,
                               a_l = ak)
    ak         <- Update_a_l(effe_test = as.matrix(effe.test),L = L,tau_l = tau,
                             R = R,s = s_tau,t = t,OneOverR = OneOverR)
    #######################################################################
    effe.test  <- Update_effe_test_t_cpp(Y = Y,
                                         alphabeta =  alphabeta, L =  L,
                                         R =  R, 
                                         ak = ak,
                                         t = t,
                                         tau_k = tau,
                                         sigma_test = as.matrix(sigma.test))
    
    sigma.test <- Update_sigma_test_t_cpp(Y = Y, alphabeta = alphabeta, t = t,
                                          asig = a_H,bsig =  b_H, L = L,
                                          effe_test = effe.test)
    #######################################################################
    if (learning_type %in% c("transductive", "training_stochastic")) {
      
      effe.train  <- Update_effe_train_t_cpp(
        Y = Y,
        G = J,
        alphabeta =  alphabeta,
        t = t,
        f_bar_g = effe.train0,
        sigma_train = sigma.train,
        KAPPAG = KappaG
      )
      
      sigma.train <- Update_sigma_train_t_cpp(
        Y = Y,
        alphabeta = alphabeta,
        t = t,
        effe_train = effe.train,
        G = J,
        a_priorG = A_B_IG_Train_Stochastic[, , 1],
        b_priorG = A_B_IG_Train_Stochastic[, , 2]
      )
    }
    #######################################################################
    effe_all  <- cbind(effe.train, effe.test)
    sigma_all <- cbind(sigma.train,sigma.test)
    ###############################################################################
    ZETA <- Upd_ZETA_t_SLICE(Y = Y,effe_ALL = effe_all,
                         sigma_ALL = sigma_all,
                         Uitilde = Ui,
                         xitilde = xi,
                         log_omegatilde = log_omegatilde,
                         J = J,n = n,L = L,poss_lab = PL)
    ########################################################################sicuri che funzioni?
    ind2         <- ZETA[ZETA>J]
    u.ind        <- unique(sort(ind2))-J
    effe.test    <- effe.test[,u.ind]
    sigma.test   <- sigma.test[,u.ind]
    ZETA[ZETA>J] <- as.numeric(factor(ZETA[ZETA>J]))+J
    
    alphabeta <- zeta_to_alphabeta(ZETA,n,J)
    #############################################################################################
    if(fixed_alphaDP  & sim == 1){
      aDP      <- aDP
    }else if (fixed_alphaDP == F){
      beta0  <- alphabeta[,2]
      uz     <- length(unique(beta0[beta0>0]))
      eta    <- rbeta(1,aDP+1,n)
      Q      <- (a_alpha+uz-1)/(n*(b_alpha-log(eta)))
      pi_eta <- Q/(1+Q)
      aDP    <- ifelse(runif(1)<pi_eta,  rgamma(1,a_alpha+uz,   b_alpha-log(eta)),
                       rgamma(1, a_alpha+uz-1,b_alpha-log(eta))  )
    }
    ############################################################################################
    
    
    if (sim > burn_in && ((sim - burn_in) %% thinning == 0)) {
      rr                <- floor((sim - burn_in)/thinning);
      PROB[rr,]         <- pidir
      ALPHA.BETA[,,rr]  <- alphabeta
      if (learning_type %in% c("transductive", "training_stochastic")) {
        EFFE.TRAIN[, , rr]  <- effe.train
        SIGMA.TRAIN[, , rr] <- sigma.train
      }
      OMEGA[[rr]]       <- omega
      EFFE.TEST[[rr]]   <- effe.test
      SIGMA.TEST[[rr]]  <- sigma.test
      TAU[[rr]]         <- tau
      AK[[rr]]          <- ak 
      AlphaDP[rr]       <- aDP
      USLICE[,rr]       <- Ui
    }
    ################################################################   
    ################################################################   
    #   if (sim%%(verbose.step*thinning) == 0) {
    #     cat(paste("Sampling iteration: ", sim, " out of ",nsim*thinning + burn_in,"\n"))}
    # }
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }
  
  # If training fixed I only report the output from robust information extraction
  if (learning_type %in% c("inductive", "training_fixed")) {
    EFFE.TRAIN  <- effe.train
    SIGMA.TRAIN <- sigma.train
  }
  out <-  list(
    P   =  PROB,
    AB  =  ALPHA.BETA,
    O   =  OMEGA,
    TAU =  TAU,
    AK  =  AK,
    aDP =  AlphaDP,
    FTr =  EFFE.TRAIN,
    STr =  SIGMA.TRAIN,
    FTr0 =  effe.train0,
    STr0 =  sigma.train0,
    FTe =  EFFE.TEST,
    STe =  SIGMA.TEST,
    USL = USLICE)
  return(out)
}

