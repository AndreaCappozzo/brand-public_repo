
library(tidyverse)
library(here)
library(rsample)
library(mclust)
library(mcclust.ext)
library(caret)
library(fda)
library(patchwork)

#library(BNPadaptRDA)

Rcpp::sourceCpp("code/functional_case/SLICE_functional.cpp")
source(here("code/functional_case/SLICE_functional.R"))

# spectra
#NIR:780-2500
#IR:2500-250000
load("data/meatdata.Rdata")

set.seed(33)
t <- length(wavelength)
x <- wavelength
df <- spectra
cl <- factor(type)

extr_sup <- 2500
matplot(x = wavelength[wavelength<extr_sup], y=spectra[which(wavelength<extr_sup),], col = factor(cl), type = "l")
df <- spectra[which(wavelength<extr_sup),]
x <- wavelength[wavelength<extr_sup]
# Define the basis
n_basis <- 100

basis_order <- 5 # spline
basis <-
  create.bspline.basis(rangeval = range(x),
                       nbasis = n_basis,
                       norder = basis_order)

basis_evaluated <- eval.basis(evalarg = x, basisobj = basis)

in_train <- createDataPartition(y = cl, times = 1, p = .5)[[1]]

extra_group_in_train <-  in_train %in% which(cl == "Beef")
in_train_without_extra_group <-
  in_train[!extra_group_in_train]
in_test_with_extra_group <-
  setdiff(1:length(cl), in_train_without_extra_group)
cl_train <- cl[in_train_without_extra_group]
cl_test <- cl[in_test_with_extra_group]

df <- apply(df,2,scale) 


X <- df[, in_train_without_extra_group]
Y <- df[, in_test_with_extra_group]

# I add 4 outliyng curves, as done in \cite{FernandezPierna2007a}
matplot(X, col = cl_train, type = "l")
matplot(Y, col = cl_test, type = "l")
cl_test
obs_to_be_modified <- c(73,74,41, 76)
cl_test[obs_to_be_modified]
cl_no_modified_obs <- cl_test[-obs_to_be_modified]
# Shift
out_1 <- c(Y[-c(1:15),obs_to_be_modified[1]],rowMeans(Y[1036:1050, which(cl_test == "Pork")]))
matplot(cbind(Y[,obs_to_be_modified[1]], out_1), type="l")
# Random noise
out_2 <- Y[,obs_to_be_modified[2]] + rnorm(n = nrow(Y),sd = .25)
matplot(cbind(Y[,obs_to_be_modified[2]], out_2), type="l")
# Spike at the end of the spectrum
out_3 <- Y[,obs_to_be_modified[3]]
out_3[800] <- 3+out_3[800]
matplot(cbind(Y[,obs_to_be_modified[3]], out_3), type="l")
# adding slope
out_4 <- Y[,obs_to_be_modified[4]]*1.2
matplot(cbind(Y[,obs_to_be_modified[4]], out_4), type="l")
Y_out <- cbind(Y[,-obs_to_be_modified], out_1,out_2,out_3, out_4)
cl_test_out <-
  factor(c(cl_no_modified_obs, rep("out", 4)),
         labels = c("Beef", "Chicken", "Lamb", "Pork", "Turkey", "out"))
matplot(Y_out, col = cl_test_out, type = "l")
A <- X[, which(cl_train == "Chicken")] 
D <- X[, which(cl_train == "Turkey")] 
C <- X[, which(cl_train == "Pork")]
B <- X[, which(cl_train == "Lamb")]

single_func_set <- A
extract_robust_tr_coef_splines <-
  Vectorize(function(single_func_set,
                     x = x,
                     h_MCD = .75) {
    smooth_representation  <-   smooth.basis(
      argvals = x,
      y = as.matrix(single_func_set),
      fdParobj = basis
    )
    basis_coef_matrix <- t(smooth_representation$fd$coefs)
    robust_procedure <- rrcov::CovMrcd(x = basis_coef_matrix, alpha = h_MCD)
    # regularized version, suitable for high dimensional data
    robust_procedure$center
  }, vectorize.args = "single_func_set")

beta_coefs_tr <-
  extract_robust_tr_coef_splines(
    single_func_set = list(A, B, C, D),
    x = x,
    h_MCD = .5
  )

matplot(D,
        lty = 2,
        type = "l", col = 1)

matplot(
  basis_evaluated %*% beta_coefs_tr[, 4],
  lty = 2,
  type = "l",
  col = 3,
  add = T,
  lwd = 2
)


effe.train  <- basis_evaluated%*%beta_coefs_tr

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

sigma.train <- cbind(sigma.trainA,sigma.trainB,sigma.trainC,sigma.trainD)

J <- length(unique(cl_train))

prior <- list(
  
  basis       = basis_evaluated,
  beta.train  = beta_coefs_tr,
  aDir     = c(rep(1,J),.5), 
  
  sigma.train = sigma.train,
  
  # prior varianza gruppi nuovi
  a_H      = 10,
  b_H      = 1,
  
  # varianza dei betini
  a_tau    = 1,
  b_tau    = 2,
  
  aDP      = .1,
  a_alphaDP = 1,
  b_alphaDP = 1,
  
  s_tau = .1,
  # Necessary for stochastic train
  vg = 1e-10,
  KappaG = 1e-10
)


set.seed(33)
RES <- MCMCap.funct.SLICE(Y = as.matrix(Y_out),
                          prior = prior,
                          L=10, # L numero possibili gruppi
                          nsim=5000, 
                          thinning=1, 
                          burn_in=5000,
                          verbose=1,
                          fixed_alphaDP = FALSE,
                          kappa = .5,
                          learning_type = "inductive")

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
    NOV_title = "DP clustering of novelties"
  ) %>%
  arrange(WL)

# Avg of novel components
df_novel_groups <-
  lapply(1:1000, function(nsim)
    RES$FT[[nsim]][, c(5, 6)]) %>%
  Reduce(f = "+") / 1000

DF_EFFE_TEST <-
  t(df_novel_groups) %>%
  as_tibble() %>%
  set_names(nm = paste0("WL", x)) %>%
  rowid_to_column() %>%
  mutate(
    cl_novelty = factor(
      c(1,2),
      labels = c("Novelty 2", "Novelty 1")
    )
  ) %>%
  pivot_longer(-c(rowid, cl_novelty),
               names_to = "WL",
               values_to = "Absorbance") %>%
  mutate(WL = parse_integer(str_extract(WL, pattern = "\\d+"))) %>%
  arrange(WL)

# Prob of being a novelty

plot_A <- ggplot(DF_4_PLOT, aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=PPN, alpha=.1),
            show.legend = c(alpha=FALSE)) +
  # xlab("Wavenumber (cm-1)") +
  labs(x="", col=quote(hat(PPN))) +
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
plot_A
# ggsave("ppn_clas_meat_data_with_outliers.pdf",width = 12.1/3*2, height = 7.19/3*2)
ggsave("ppn_clas_meat_data_with_outliers.pdf",width = 10.7, height = 4.22)

# Classification via Variation of Information loss for novelty component (avg EFFE TEST superimposed)
plot_C <- ggplot(filter(DF_4_PLOT, cl_alpha == "Novelty"), aes(WL, Absorbance)) +
  geom_path(aes(group = rowid, col = cl_beta, lty = cl_beta), alpha = .3) +
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
  geom_path(aes(group = rowid, col = "Turkey", lty="Turkey"),
            alpha = .3) +
  geom_path(
    data = filter(DF_4_PLOT, cl_true == "Turkey", cl_alpha == "Novelty"),
    aes(group = rowid, col = "Outlying Turkey", lty = "Outlying Turkey"),
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
ggsave("clas_novelty_meat_data_with_outliers.pdf",width = 12.1/3*2, height = 7.19/3*2)
  # ggsave("clas_novelty_meat_data.pdf",width = 12.1/2, height = 7.19/2)

# Faceting plot B

ggplot(DF_4_PLOT, aes(WL, Absorbance)) +
  geom_path(aes(group=rowid, col=cl_alpha, lty=cl_alpha), alpha=.6) +
  labs(x="Wavelength (nm)", color="Meat type", lty="Meat type") +
  theme_bw() +
  # theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(CL_title~cl_alpha)

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

