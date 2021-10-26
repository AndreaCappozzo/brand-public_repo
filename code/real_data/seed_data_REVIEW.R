
library(tidyverse)
library(here)
library(rsample)
library(mclust)
library(mcclust.ext)
library(GGally)
library(brand)
#library(BNPadaptRDA)

# Rcpp::sourceCpp(here('code/January2021/MLVT_Slice.cpp'))
# source(here('code/January2021/MLVT_Slice.R'))

seed_col_names <- c(
  "area",
  "perimeter",
  "compactness",
  "length",
  "width",
  "asymmetry",
  "length groove",
  "seed_type"
)
seeds_dataset <- read_table2(here("data/seeds_dataset.txt"),
                             col_names = seed_col_names)
set.seed(1245)

initial_tr_te <-
  initial_split(data = seeds_dataset,
                prop = .5,
                strata = "seed_type")
X_tr <- training(initial_tr_te)
Y_te <- testing(initial_tr_te)

disc_variable <- "seed_type"


# X_tr_with_LB <- X_tr
#
# obs_label_noise <- sample(x = which(x =X_tr$seed_type==3), size = 10, replace = FALSE)
# X_tr_with_LB$seed_type[obs_label_noise] <- 4

# X <- filter(X_tr_with_LB, get(disc_variable) != 3)
X <- filter(X_tr, get(disc_variable) != 3)
# Y <- bind_rows(Y_te, filter(X_tr, get(disc_variable) == 3))
Y <- Y_te

# obs_LB <- which(X$seed_type==4)
# X$seed_type[obs_LB] <- 1

color_values <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73")
shape_values <- c("1" = 5, "2" = 20, "3" = 18)

# Learning scenario plot --------------------------------------------------

df_4_plot <- bind_rows(mutate(X, set = "Training set"),
                       mutate(Y, set = "Test set")) %>%
  mutate(set = factor(
    set,
    levels = c("Training set", "Test set"),
    ordered = T
  ),
  seed_type = factor(seed_type))

ggplot(df_4_plot) +
  geom_point(aes(compactness, perimeter, shape = seed_type, col = seed_type)) +
  facet_wrap( ~ set) + theme_bw() +
  scale_color_manual("Seed Type", values = color_values) +
  scale_shape_manual("Seed Type", values = shape_values) +
  theme(legend.position = "bottom")#, panel.spacing = unit(5, "lines"))

# ggsave(filename = "seed_learning_scenario.pdf",
#        width = 9.61/2,
#        height =  5.78/2)

pairs(X[, -8], col = X$seed_type)
pairs(Y[, -8], col = Y$seed_type)
categ <- X$seed_type
p <- ncol(X) - 1
J <- length(unique(categ))

prior <- list(
  aDir = c(rep(1, J), .1),
  aDP = 1,
  #(.5?)
  m_H = rep(0, p),
  k_H = .01,
  v_H = 10,
  S_H = diag(p) ,
  k_g = 25,
  v_g = 1000,
  a_alphaDP = 1,
  b_alphaDP = 1
)


set.seed(12345)

mod <-
  Brand_mlvt(
    Y = as.matrix(Y[, -8]),
    X = as.matrix(X[, -8]),
    categ = categ,
    prior = prior,
    L = 50,
    burn_in = 20000,
    thinning = 1,
    NSIM = 10000,
    fixed_alphaDP = F,
    h_MCD = .95,
    raw_MCD = F,
    learning_type = "inductive"
  )



x <- 2
y <- 3
plot(Y[,c(x,y)])
# Train fixed
lines(ellipse::ellipse(mod$Sig.train[c(x,y),c(x,y),1],
                       centre = mod$Mu.train[1,c(x,y)]))
lines(ellipse::ellipse(mod$Sig.train[c(x,y),c(x,y),2],
                       centre = mod$Mu.train[2,c(x,y)]))
# Train sthocastic
# lines(ellipse::ellipse(mod$Sig.train[c(x,y),c(x,y),1,1000],
#                        centre = mod$Mu.train[1,c(x,y),1000]))
# lines(ellipse::ellipse(mod$Sig.train[c(x,y),c(x,y),2,1000],
#                        centre = mod$Mu.train[2,c(x,y),1000]))
lines(ellipse::ellipse(mod$s_b[c(x,y),c(x,y),1],
                       centre = mod$x_b[c(x,y),1]),col=2,lwd=2,lty=2)
lines(ellipse::ellipse(mod$s_b[c(x,y),c(x,y),2],
                       centre = mod$x_b[c(x,y),2]),col=2,lwd=2,lty=2)



hist(mod$alphaDP)
ALP    <- mod$AB[, 1, ]

major.vote <- function(x) {
  as.numeric(names(sort(table(x),decreasing = T))[1])
}

dim(mod$Mu.train)

cl_pred <- apply(mod$AB[, 1, ], 1, major.vote)
class_y <- Y$seed_type
table(cl_pred, Y$seed_type)

plot(Y[,c(x,y)],pch=21, col = cl_pred + 1)
plot(Y[,c(x,y)],pch=21, bg = class_y)
points(X[,c(x,y)],pch=21, bg = categ)


plot(Y[,c(x,y)],pch=20+class_y, bg = cl_pred + 1)


ALP <- mod$AB[, 1,]
psmALP <- comp.psm(t(ALP) + 1)
image(psmALP)
cl_pred_binder <- minbinder(psmALP)$cl
cl_pred_VI <- minVI(psmALP)$cl
table(cl_pred_binder, cl_pred)
table(cl_pred_binder, cl_pred_VI)
table(Y$seed_type, cl_pred_VI)
table(Y$seed_type, cl_pred_binder)

plot(Y, col = cl_pred_binder + 1)
# Evaluate class accuracy
class_y <- Y$seed_type
mclust::classError(classification = cl_pred_binder, class = class_y)
mclust::classError(classification = cl_pred, class = class_y)
table(classification = cl_pred, class = class_y)
table(classification = cl_pred_binder, class = class_y)

newobs <- cl_pred == 0
#newobs <- cl_pred_binder == 2

pnew <- (apply(mod$AB[, 1,] == 0, 1, mean))
ggpairs_plot_pnew <-
  ggpairs(
    Y[, -8],
    aes(col = pnew),
    upper = list(continuous = "points"),
    lower = list(
      continuous = "points",
      combo =
        "facetdensity",
      discrete = "ratio",
      na = "na"
    ),
    diag = list(continuous = "blankDiag")
  )
for (i in 1:ggpairs_plot_pnew$nrow) {
  for (j in 1:ggpairs_plot_pnew$ncol) {
    ggpairs_plot_pnew[i, j] <- ggpairs_plot_pnew[i, j] +
      scale_color_viridis_c() + theme_bw()
  }
}

color_values <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73")
# color_values <- c("1" = "#E69F00", "2" = "#56B4E9", "0" = "#009E73")
shape_values <- c("1" = 5, "2" = 20, "3" = 18)
# shape_values <- c("1" = 5, "2" = 20, "0" = 18)

ggpairs_plot_cl_pred <-
  ggpairs(
    Y[, -8],
    aes(col = as.character(cl_pred_binder), shape=as.character(cl_pred_binder)),
    # aes(col = as.character(cl_pred), shape=as.character(cl_pred)),
    upper = list(continuous = "points"),
    lower = list(
      continuous = "points",
      combo =
        "facetdensity",
      discrete = "ratio",
      na = "na"
    ),
    diag = list(continuous = "blankDiag")
  )


for (i in 1:ggpairs_plot_pnew$nrow) {
  for (j in 1:ggpairs_plot_pnew$ncol) {
    if (j > i) {
      ggpairs_plot_pnew[i, j] <- ggpairs_plot_cl_pred[i, j] +
        scale_fill_manual(values = color_values) +
        scale_color_manual(values = color_values) +
        scale_shape_manual(values = shape_values) +
        theme_bw()
    } else{
      next
    }
  }
}
ggpairs_plot_pnew
ggsave(filename = "seed_BNP_ggplots_AF.pdf")

# How many new groups? FIXME ----------------------------------------------------

BET <- mod$AB[newobs, 2, ]
psmBET <- comp.psm(t(BET) + 1)
image(psmBET)
clB <- minbinder(psmBET)$cl
clV <- minVI(psmBET)$cl

plot(Y[, 1:2], pch = ".", cex = 2)
points(Y[newobs,1:2 ], col = clB)

plot(Y, pch = ".", cex = 2)
points(Y[newobsV, ], col = clV)



# plot(Y[newobs,],col=minbinder(psmBET)$cl)
# abline(h=-3.5)

pnew <- (apply(mod$AB[, 1, ] == 0, 1, mean))
qplot(pull(Y[, 1]), pull(Y[, 2]), col = pnew) + scale_color_viridis_c() + theme_bw()
