
library(tidyverse)
library(here)
library(rsample)
library(mclust)
library(mcclust.ext)
library(caret)
library(fda)
library(patchwork)

#library(BNPadaptRDA)

library(brand)

# spectra
#NIR:780-2500
#IR:2500-250000
load("../adaptdaBNP/data/meatdata.Rdata")

set.seed(33)
t <- length(wavelength)
x <- wavelength
df <- spectra
# df <- apply(df,2,function(x)x-mean(x))
cl <- factor(type)

extr_sup <- 2500
#extr_sup <- 1100
matplot(x = wavelength[wavelength<extr_sup], y=spectra[which(wavelength<extr_sup),], col = factor(cl), type = "l")
df <- spectra[which(wavelength<extr_sup),]
x <- wavelength[wavelength<extr_sup]
# Define the basis
n_basis <- 100
# n_basis <- 800

basis_order <- 5 # spline
basis <-
  create.bspline.basis(rangeval = range(x),
                       nbasis = n_basis,
                       norder = basis_order)

basis_evaluated <- eval.basis(evalarg = x, basisobj = basis)

# func_df <- Data2fd(argvals=x, 
#         y=as.matrix(spectra), basisobj=basis) # create functional data
# df <- eval.fd(evalarg = x, fdobj = func_df, Lfdobj = 1)
# plot.fd(func_df, col=factor(cl), Lfdobj = 1)
# matplot(df, col=factor(cl), type="l")
in_train <- createDataPartition(y = cl, times = 1, p = .5)[[1]]

extra_group_in_train <-  in_train %in% which(cl == "Beef")
in_train_without_extra_group <-
  in_train[!extra_group_in_train]
in_test_with_extra_group <-
  setdiff(1:length(cl), in_train_without_extra_group)
cl_train <- cl[in_train_without_extra_group]
cl_test <- cl[in_test_with_extra_group]

###################################################################################
df <- apply(df,2,scale) # FIXME shall we think about robust standardization ?
###################################################################################


X <- df[, in_train_without_extra_group]
Y <- df[, in_test_with_extra_group]

# I add 4 outliyng curves, as done in \cite{FernandezPierna2007a}
# matplot(X, col = cl_train, type = "l",ylim=c(0.5,1.5))
# matplot(Y, col = cl_test, type = "l",ylim=c(0.5,1.5))
matplot(X, col = cl_train, type = "l")
matplot(Y, col = cl_test, type = "l")
cl_test
# obs_to_be_modified <- c(35, 41, 43, 45)
obs_to_be_modified <- c(73,74,41, 76)
cl_test[obs_to_be_modified]
cl_no_modified_obs <- cl_test[-obs_to_be_modified]
# Shift
# out_1 <- c(Y[-c(1:6),obs_to_be_modified[1]],rowMeans(Y[1045:1050, which(cl_test == "Turkey")]))
out_1 <- c(Y[-c(1:15),obs_to_be_modified[1]],rowMeans(Y[1036:1050, which(cl_test == "Pork")]))
matplot(cbind(Y[,obs_to_be_modified[1]], out_1), type="l")
# Random noise
# out_2 <- Y[,obs_to_be_modified[2]] + rnorm(n = nrow(Y),sd = .1)
out_2 <- Y[,obs_to_be_modified[2]] + rnorm(n = nrow(Y),sd = .25)
matplot(cbind(Y[,obs_to_be_modified[2]], out_2), type="l")
# Spike at the end of the spectrum
out_3 <- Y[,obs_to_be_modified[3]]
# out_3[800] <- .5+out_3[800]
out_3[800] <- 3+out_3[800]
matplot(cbind(Y[,obs_to_be_modified[3]], out_3), type="l")
# adding slope
# out_4 <- Y[,obs_to_be_modified[4]]*1.2
out_4 <- Y[,obs_to_be_modified[4]]*1.2
matplot(cbind(Y[,obs_to_be_modified[4]], out_4), type="l")
Y_out <- cbind(Y[,-obs_to_be_modified], out_1,out_2,out_3, out_4)
cl_test_out <-
  factor(c(cl_no_modified_obs, rep("out", 4)),
         labels = c("Beef", "Chicken", "Lamb", "Pork", "Turkey", "out"))
matplot(Y_out, col = cl_test_out, type = "l")
# A <- X[, which(cl_train == "Poultry")] # Poultry: known groups tr set
A <- X[, which(cl_train == "Chicken")] 
B <- X[, which(cl_train == "Lamb")]
C <- X[, which(cl_train == "Pork")]
D <- X[, which(cl_train == "Turkey")] 
# D <- X[, which(cl_train == "Beef")]
# AY <- Y[, which(cl_test == "Chicken")] # Poultry: known groups tr set
# BY <- Y[, which(cl_test == "Turkey")] # Poultry: known groups tr set
# CY <- Y[, which(cl_test == "Pork")]
# DY <- Y[, which(cl_test == "Lamb")]

# beta_coefs_tr <-
#   extract_tr_coef_splines(single_func_set = list(A, C, D), x = x)
beta_coefs_tr <-
  extract_robust_tr_v2(
    single_func_set = list(A, B, C, D),
    x = x,basis_evaluated = basis_evaluated, 
    h_MCD = .75
  )

matplot(D,
        lty = 2,
        type = "l", col = 1)

# matplot(
#   DY,
#   lty = 2,
#   type = "l",
#   col = 2,
#   add = T
# )


effe.train  <- cbind(beta_coefs_tr[1,][[1]],
                     beta_coefs_tr[1,][[2]],
                     beta_coefs_tr[1,][[3]],
                     beta_coefs_tr[1,][[4]])
beta.train  <- cbind(beta_coefs_tr[3,][[1]],
                     beta_coefs_tr[3,][[2]],
                     beta_coefs_tr[3,][[3]],
                     beta_coefs_tr[3,][[4]])

matplot(effe.train)

sigma.trainA <-
  apply((apply(A, 2, function(z)
    z -   effe.train[, 1])), 1, function(k)
      sum(k ^ 2) / (length(k) - 1))
sigma.trainB <-
  apply((apply(B, 2, function(z)
    z -   effe.train[, 2])), 1, function(k)
      sum(k ^ 2) / (length(k) - 1))
sigma.trainC <-
  apply((apply(C, 2, function(z)
    z -   effe.train[, 3])), 1, function(k)
      sum(k ^ 2) / (length(k) - 1))
sigma.trainD <-
  apply((apply(D, 2, function(z)
    z -   effe.train[, 4])), 1, function(k)
      sum(k ^ 2) / (length(k) - 1))

sigma.train <- cbind(beta_coefs_tr[2,][[1]],
                     beta_coefs_tr[2,][[2]],
                     beta_coefs_tr[2,][[3]],
                     beta_coefs_tr[2,][[4]])

# sigma.train_old <- cbind(sigma.trainA,
#                          sigma.trainB,
#                          sigma.trainC,
#                          sigma.trainD)



plot(ts(sigma.train),type="l",lty=1) 
#plot(ts(sigma.train_old),type="l",lty=1) 
 
 
# var cresce a manetta at the end of the wavelength
# sigma.train <- c(sigma.trainA, sigma.trainC,sigma.trainD)
J <- length(unique(cl_train))





prior <- list(
  
  basis       = basis_evaluated,
  beta.train  = beta.train,
  aDir     = c(rep(1,J),.5), 
  
  sigma.train = sigma.train,
  
  # prior varianza gruppi nuovi
  a_H      = 5,
  b_H      = 1,
  
  # varianza dei betini
  a_tau    = 3,
  b_tau    = 1,
  s_tau    = 1,
  aDP      =  1,
  a_alphaDP = 3,
  b_alphaDP = 3,
  
  # Necessary for stochastic train
  vg = 1e-5,
  KappaG = 1e-5
)


set.seed(3111)
RES <- Brand_fct(Y = as.matrix(Y_out),
                          prior = prior,
                          L=20, # L numero possibili gruppi
                          nsim=10000, 
                          thinning=1, 
                          burn_in=10000,
                          verbose=1,
                          fixed_alphaDP = 0,
                          kappa = .5,
                          learning_type = "inductive",light.output = T)

#saveRDS(RES,"Codice_review/Diagram//new_functional_meat.RDS")

RES

#RES <- saveRDS(RES,"../adaptdaBNP/Codice_review/Diagram/new_functional_meat_26Jan.RDS")
RES <- readRDS("../adaptdaBNP/Codice_review/Diagram/new_functional_meat_26Jan.RDS")

RES$prior
plot(as.ts(RES$aDP))


major.vote <- function(x){
  as.numeric(names(sort(table(x),decreasing = TRUE))[1]) 
}


cl_alpha <- apply(RES$AB[,1,],1,major.vote)

table(cl_test_out, cl_alpha)
table(cl_alpha, cl_test_out)

classError(classification = cl_alpha, class = cl_test_out)

xtable::xtable(table(cl_alpha, cl_test_out))

PPN <- rowMeans(RES$AB[,1,]==0)

new_obs <- which(cl_alpha==0)


BET    <- RES$AB[new_obs,2,]
psmBET <- comp.psm(t(BET)+1)
image(psmBET)
cl_beta_binder <- minbinder(psmBET)$cl
cl_beta_VI <- minVI(psmBET)$cl
cl_beta <- cl_alpha
cl_beta[cl_alpha==0] <- cl_beta_VI+J # I add J to separate from the original tr labels
matplot(Y,type="l",col=cl_beta+1) 
table(cl_beta_VI, cl_test_out[cl_alpha==0]) # cool: it identifies two classes:
# one beef and the other in which units from the previous classes
# (that were wrongly labeled as novelty) are gathered
cl_test_out[which(cl_beta==5)]
table(cl_beta, cl_test_out)

# Plotting

DF_4_PLOT <-
  t(Y_out) %>%
  as_tibble() %>%
  set_names(nm = paste0("WL", x)) %>%
  rowid_to_column() %>%
  mutate(
    PPN = PPN,
    cl_alpha = factor(
      cl_alpha,
      labels = c("Novelty", "Chicken", "Lamb", "Pork", "Turkey")
    ),
    cl_true = cl_test_out,
    cl_beta = factor(
      cl_beta,
      labels = c(
        "Chicken",
        "Lamb",
        "Pork",
        "Turkey",
        "Novelty 1",
        "Novelty 2"
      )
    )
  ) %>%
  dplyr::select(rowid, PPN, everything()) %>%
  pivot_longer(
    -c(rowid, PPN, cl_alpha, cl_beta, cl_true),
    names_to = "WL",
    values_to = "Absorbance"
  ) %>%
  mutate(
    WL = parse_integer(str_extract(WL, pattern = "\\d+")),
    PPN_title = "Posterior probability of being a novelty",
    CL_title = "A-posteriori classification (majority vote)",
    CL2_title = "Classified as:",
    NOV_title = "DP clustering of novelties"
  ) %>%
  arrange(WL)

# Avg of novel components
# df_novel_groups <-
#   lapply(1:10000, function(nsim)
#     RES$FTe[[nsim]][, c(1, 2)]) %>%
#   Reduce(f = "+") / 1000

# DF_EFFE_TEST <-
#   t(df_novel_groups) %>%
#   as_tibble() %>%
#   set_names(nm = paste0("WL", x)) %>%
#   rowid_to_column() %>%
#   mutate(
#     cl_novelty = factor(
#       c(1,2),
#       labels = c("Novelty 2", "Novelty 1")
#     )
#   ) %>%
#   pivot_longer(-c(rowid, cl_novelty),
#                names_to = "WL",
#                values_to = "Absorbance") %>%
#   mutate(WL = parse_integer(str_extract(WL, pattern = "\\d+"))) %>%
#   arrange(WL)

# Prob of being a novelty

plot_A <- ggplot(DF_4_PLOT, aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=PPN, alpha=.1),
            show.legend = c(alpha=FALSE)) +
  # xlab("Wavenumber (cm-1)") +
  labs(x="Wavelength (nm)",x="", col=quote(hat(PPN))) +
  scale_color_viridis_c()+
  theme_bw() +
  facet_wrap(~PPN_title)


# Classification via majority vote

plot_B <- ggplot(DF_4_PLOT, aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=cl_alpha, lty=cl_alpha), alpha=.6) +
  labs(x="Wavelength (nm)", color="Meat type", lty="Meat type") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~CL_title)

# plot_A/plot_B
plot_A+theme(strip.text.x = element_text(size=12))
# ggsave("ppn_clas_meat_data_with_outliers.pdf",width = 12.1/3*2, height = 7.19/3*2)
ggsave("Codice_review/Diagram/ppn_clas_meat_data_with_outliers.pdf",width = 10.7, height = 4.22)

# Classification via Variation of Information loss for novelty component (avg EFFE TEST superimposed)
plot_C <- ggplot(filter(DF_4_PLOT, cl_alpha == "Novelty"), aes(WL, Absorbance)) +
  geom_path(aes(group = rowid, col = cl_beta), alpha = .3) +
  # geom_path(
  #   data = DF_EFFE_TEST,
  #   mapping = aes(WL, Absorbance, col = cl_novelty),
  #   cex = .8
  # ) +
  labs(x = "", color = "Meat type", lty = "Meat type") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~NOV_title)

# Focus on the outlying turkey


plot_D <-
  ggplot(mutate(filter(DF_4_PLOT, cl_true == "Turkey"), tt = "Focus on Turkey sub-population"),
         aes(WL, Absorbance)) +
  geom_path(aes(group = rowid, col = "Turkey"),
            alpha = .3) +
  geom_path(
    data = filter(DF_4_PLOT, cl_true == "Turkey", cl_alpha == "Novelty"),
    aes(group = rowid, col = "Outlying Turkey"),
    alpha = .7
  ) + labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  theme_bw() +
  facet_wrap( ~ tt) +
  # theme(legend.position = "bottom") +
  scale_color_manual(values = c(Turkey = "#66A61E", # same color as in the first plot
                                "Outlying Turkey"= "darkred")) +
  scale_linetype_manual(values = c(Turkey = 4, # same color as in the first plot
                                   "Outlying Turkey"= 1))


plot_C/plot_D
ggsave("Codice_review/Diagram/clas_novelty_meat_data_with_outliers.pdf",width = 12.1/3*2, height = 7.19/3*2)


(plot_C+theme(legend.position = "bottom",strip.text.x = element_text(size=12)))+
  (plot_D+theme(legend.position = "bottom",strip.text.x = element_text(size=12)))
ggsave("Codice_review/Diagram/clas_novelty_meat_data_with_outliers_long.pdf",width = 12.1/3*3, height = 7.19/3*2)
# ggsave("clas_novelty_meat_data.pdf",width = 12.1/2, height = 7.19/2)

# Faceting plot B

ggplot(DF_4_PLOT %>% mutate(cl_alpha=paste("Classified as:",cl_alpha)), aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=cl_true), alpha=.7) +
  labs(x="Wavelength (nm)", color="Meat type", lty="Meat type") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_manual(values = c(1,"royalblue",4,"orange","forestgreen",2)) +
  facet_wrap(~cl_alpha,ncol = 2)+ theme(legend.position = c(0.75, 0.15),
                                        strip.text.x = element_text(size=12))
ggsave("Codice_review/Diagram/allclassfunctional2.pdf",width = 12.1/3*2, height = 7.19/3*3)

# Mean value plot B

effe.train


EFFE_TRAIN_AVG_4_PLOT <-
  t(effe.train) %>%
  as_tibble() %>%
  set_names(nm = paste0("WL", x)) %>%
  rowid_to_column() %>%
  mutate(
    cl_alpha = factor(
      c("Chicken", "Lamb", "Pork", "Turkey"),
      levels = c("Novelty","Chicken", "Lamb", "Pork", "Turkey"),
    )) %>%
  dplyr::select(rowid, everything()) %>%
  pivot_longer(
    -c(rowid, cl_alpha),
    names_to = "WL",
    values_to = "Absorbance"
  ) %>%
  mutate(
    WL = parse_integer(str_extract(WL, pattern = "\\d+"))
    # AVG_title = ""
  ) %>%
  arrange(WL)

# DF_4_PLOT %>% 
#   group_by(WL,cl_alpha) %>% 
#   mutate(Absorbance=mean(Absorbance)) %>% 
#   ungroup() %>% 
#   distinct(WL,Absorbance, .keep_all = TRUE) %>% 
filter(DF_4_PLOT, cl_alpha == "Novelty") %>% 
  ggplot(aes(WL, Absorbance)) +
  geom_path(
    aes(group = rowid, col = cl_alpha, lty = cl_alpha),
    alpha = .2
  ) +
  geom_path(data=EFFE_TRAIN_AVG_4_PLOT,
            aes(group = rowid, col = cl_alpha, lty = cl_alpha), alpha =
              1) +
  labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")













# Why Chicken and Turkey are the same? ------------------------------------










filter(DF_4_PLOT, cl_alpha == "Turkey"| cl_alpha =="Chicken") %>% 
  ggplot(aes(WL, Absorbance)) +
  geom_path(
    aes(group = rowid, col = cl_alpha),
    alpha = 1
  ) +theme_bw()+scale_color_viridis_d()+
  labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  theme_bw() 

filter(DF_4_PLOT, cl_alpha == "Turkey"| cl_alpha =="Chicken") %>% 
  ggplot(aes(WL, Absorbance)) +
  geom_path(
    aes(group = rowid, col = cl_alpha),
    alpha = 1
  ) +theme_bw()+scale_color_viridis_d()+
  labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  theme_bw() 





effe.train




DF_TRAIN_PLOT <-
  reshape2::melt(effe.train) %>%
  as_tibble() %>%
  mutate(
    WL=rep(x,4),
    Absorbance=value,
    cl_alpha = factor(
      Var2,
      labels = c("Chicken", "Lamb", "Pork", "Turkey"),
    ),
    title = "Robustly Extracted Mean Functions")



convince1 <- DF_TRAIN_PLOT %>% 
  ggplot(aes(WL, Absorbance)) +
  geom_path(
    aes(group = cl_alpha, col = cl_alpha),
    alpha = 1,lwd=1
  ) +theme_bw()+scale_color_viridis_d(end=.95)+
  labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  facet_wrap(~title)+theme(legend.position = "bottom",strip.text.x = element_text(size=12))


convince2 <- filter(DF_4_PLOT %>% mutate(title="Poultry: Test Set"), cl_alpha == "Turkey"| cl_alpha =="Chicken") %>% 
  ggplot(aes(WL, Absorbance)) +
  geom_path(
    aes(group = rowid, col = cl_alpha),
    alpha = .8
  ) +theme_bw()+scale_color_manual(values=c(2,4))+
  labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  theme_bw() +  facet_wrap(~title)+
  theme(legend.position = "bottom",strip.text.x = element_text(size=12))


library(patchwork)
convince1+convince2
ggsave("Codice_review/Diagram/convinceme.pdf",width = 12.1/3*3, height = 7.19/3*2)




a <- ggplot(DF_4_PLOT %>% filter(cl_true=="Turkey") %>% mutate(cl_tit="All Turkey units") , aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=cl_alpha), alpha=.8,lwd=.6) +
  labs(x="Wavelength (nm)", color="Turkey classfied as") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_manual(values=1:4) +
  facet_wrap(~cl_tit)+
  theme(legend.position = "none")+theme(strip.text.x = element_text(size=12))
# what about within the true turkeys?
b <- ggplot(DF_4_PLOT %>% filter(cl_true=="Turkey") , aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=cl_alpha), alpha=1,lwd=.6) +
  labs(x="Wavelength (nm)", color="Turkey classfied as") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_manual(values=1:4) +
  facet_wrap(~cl_alpha)+
  theme(legend.position = "bottom")+theme(strip.text.x = element_text(size=12))

a/b
ggsave("Codice_review/Diagram/allturkeyclassfunctional.pdf",width = 12.1/3*2, height = 7.19/3*5)

















##

convince3 <- filter(DF_4_PLOT %>% 
                      mutate(title="Beef and Outliers"), cl_alpha=="Novelty") %>% 
  ggplot(aes(WL, Absorbance)) +
  geom_path(data=filter(DF_4_PLOT %>% 
                          mutate(title="Beef and Outliers"), cl_true=="out" & cl_beta == "Novelty 1"),
    aes(group = rowid,col="Outliers"),
    alpha = .7
  ) +
  geom_path(data=filter(DF_4_PLOT %>% 
                          mutate(title="Beef and Outliers"), cl_true=="Beef" & cl_beta == "Novelty 2"),
            aes(group = rowid,col="Beef"),
            alpha = .7
  )+
   geom_path(data=filter(DF_4_PLOT %>% 
              mutate(title="Beef and Outliers"), cl_true=="out" & cl_beta == "Novelty 2"),
    aes(group = rowid,col="Out w/ Beef"),
    alpha = 1
  )+
  theme_bw()+
  scale_color_manual(values=c("blue","red","gray"))+
  labs(x = "Wavelength (nm)", color = "Meat type", lty = "Meat type") +
  theme_bw() +  facet_wrap(~title)+
  theme(legend.position = "bottom",strip.text.x = element_text(size=12))
convince3
ggsave("Codice_review/Diagram/BeefandOutliers.pdf",width = 12.1/3*2, height = 7.19/3*3)
