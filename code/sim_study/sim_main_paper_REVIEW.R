

library(mvnfast)
library(mcclust)
library(mclust)
library(ggplot2)
library(mcclust.ext)
library(tidyverse)
library(forecast)
library(adaptDA) # install.packages("~/Downloads/adaptDA_1.0.tar",repos = NULL, type = "source")
library(raedda) # remotes::install_github("AndreaCappozzo/raedda")
library(foreach)
library(parallel)
library(patchwork)
library(xtable)
library(doRNG)
library(brand)

#library(BNPadaptRDA)

N <- 1000
M <- 1000
N_vec <- list(300, 300, 400)

G <- 3
H <- 7
p <- 2
cl_train <- rep(1:G, N_vec)

MU <-
  list(c(-5, 5),
       c(-4, -4),
       c(4, 4),
       c(-0, 0),
       c(5, -10),
       c(5, -10),
       c(-10, -10))

SIGMA <-
  list(
    matrix(c(1, .9, .9, 1), 2, 2),
    diag(2),
    diag(2),
    matrix(c(1, -.75, -.75, 1), 2, 2),
    matrix(c(1, .9, .9, 1), 2, 2),
    matrix(c(1, -.9, -.9, 1), 2, 2),
    diag(.01, 2)
  )

BURN_IN <- 20000
length_chain <- 20000

novelty_size <- "not_small"
label_noise <- TRUE

sim_study_adaptBNP <-
  function(label_noise = c(TRUE, FALSE),
           # v_tilde_train = c(10, 1000),
           # alpha_MCD = c(.75, 1),
           novelty_size = c("small", "not_small")) {
    M_vec <- switch(
      novelty_size,
      "small" = list(350, 250, 250, 49, 50, 50, 1),
      "not_small" = list(200, 200, 250, 90, 100, 100, 10)
    )
    
    cl_test <- rep(1:7, M_vec)
    
    cl_train_label_noise <- cl_train
    if (label_noise == TRUE) {
      obs_LB <-
        sample(which(cl_train_label_noise %in% c(2, 3)), size = 120)
      cl_train_label_noise[obs_LB] <-
        ifelse(test = cl_train_label_noise[obs_LB] == 2, 3, 2) # add label noise
    }
    
    X <-
      purrr::map_dfr(1:G,  ~ as.data.frame(mvnfast::rmvn(
        n = N_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]]
      )))
    
    Y <-
      purrr::map_dfr(1:H,  ~ as.data.frame(mvnfast::rmvn(
        n = M_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]]
      )))
    
    # plot(X, col=cl_train_label_noise)
    # plot(Y, col=cl_test)
    
    J <- length(unique(cl_train_label_noise))
    
    adapt_BNP_res_collection <- function(alpha_MCD = c(.75, 1),
                                         k_tilde_train = c(10, 1000)) {
      prior <- list(
        aDir = c(rep(1, J), .1),
        aDP = .1,
        #(.5?)
        m_H     = rep(0, p),
        k_H = .01,
        v_H = 10,
        S_H = diag(p) * 10,
        k_g = k_tilde_train,
        v_g = 5,
        a_alphaDP = 1,
        b_alphaDP = 1
      )
      
      t_adapt_BNP <- system.time(
        fit_adapt_RDA_bnp <-
          Brand_mlvt(
            Y = Y,
            X = X,
            categ = cl_train_label_noise,
            prior = prior,
            L = 20,
            burn_in = BURN_IN,
            thinning = 1,
            nsim = length_chain,
            fixed_alphaDP = FALSE,
            h_MCD = alpha_MCD,
            raw_MCD = FALSE,
            kappa = .25,
            learning_type = "transductive",
            light = TRUE
          )
      )[3]
      
      cl_adapt_BNP <-
        apply(fit_adapt_RDA_bnp$AB[, 1, ], 1, major.vote)
      novelty_adapt_BNP <- which(cl_adapt_BNP == 0)
      BET    <- fit_adapt_RDA_bnp$AB[novelty_adapt_BNP, 2,]
      psmBET <- comp.psm(t(BET) + 1)
      cl_beta_VI <- minVI(psmBET)$cl
      cl_beta_adapt_BNP <- cl_adapt_BNP
      cl_beta_adapt_BNP[cl_adapt_BNP == 0] <-
        cl_beta_VI + G # I add G to separate from the original tr labels
      a_posteriori_prob_novelty <-
        apply(fit_adapt_RDA_bnp$AB[, 1, ] == 0, 1, mean)
      results <-
        list(
          cl_alpha = cl_adapt_BNP,
          cl_beta = cl_beta_adapt_BNP,
          comp_time = t_adapt_BNP,
          a_posteriori_prob_novelty = a_posteriori_prob_novelty
        )
      results
    }
    
    # I collect results for 4 different situations of adaptBNP model
    
    adapt_BNP_cases <- purrr::cross2(c(1, .75), # alpha_MCD
                                     c(10, 1000)) # k_tilde_train
    
    fit_adapt_BNP_list <-
      purrr::map(.x = adapt_BNP_cases,
                 ~ adapt_BNP_res_collection(alpha_MCD = .x[[1]], k_tilde_train = .x[[2]]))
    
    # COMPETITORS: AMDA and RAEDDA inductive -------------------------------------------------------------
    
    # AMDAt
    # t_amda_t <- system.time(fit_amda_t_list <-
    #                           lapply(
    #                             G:H,
    #                             FUN = function(h)
    #                               amdat(
    #                                 X = X,
    #                                 cls = cl_train_label_noise,
    #                                 Y = Y,
    #                                 K = h,
    #                                 model = "qda"
    #                               )
    #                           ))[3]
    #
    # best_BIC_t <-
    #   which.max(sapply(
    #     1:length(fit_amda_t_list),
    #     FUN = function(g)
    #       fit_amda_t_list[[g]]$crit$bic
    #   ))
    #
    # fit_amda_t <- fit_amda_t_list[[best_BIC_t]]
    
    t_amda_i_1 <- Sys.time()
    fit_amda_learning <-
      adaptDA::amdai(X = X,
                     cls = cl_train_label_noise,
                     model = "qda")
    
    fit_amda_i_list <-
      lapply(
        G:H,
        FUN = function(h)
          tryCatch(
            predict.amdai(fit_amda_learning, Y = Y, K = h),
            error = function(e)
              list(crit = list(bic = NA))
          )
      )
    
    t_amda_i_2 <- Sys.time()
    t_amda_i <- difftime(t_amda_i_2, t_amda_i_1)
    best_BIC_i <-
      which.max(sapply(
        1:length(fit_amda_i_list),
        FUN = function(g)
          fit_amda_i_list[[g]]$crit$bic
      ))
    fit_amda_i <- fit_amda_i_list[[best_BIC_i]]
    
    # # RAEDDAt
    #
    # t_raedda_t <- system.time(
    #   fit_raedda_t <- raedda_t(
    #     X_train = X,
    #     class_train = cl_train_label_noise,
    #     X_test = Y,
    #     G = NULL,
    #     alpha_train = .1,
    #     alpha_test = 0.05,
    #     model_names_l = NULL,
    #     model_names_d = NULL,
    #     restr_factor_d = NULL,
    #     ctrl_init = control_init(n_samp = 20, n_start_extra_classes = 50),
    #     ctrl_EM = control_EM(nstart_EM = 100)
    #   )
    # )[3]
    
    # RAEDDAi
    
    t_raedda_i <- system.time(
      fit_raedda_i <- raedda_i(
        X_train = X,
        class_train = cl_train_label_noise,
        X_test = Y,
        G = NULL,
        alpha_train = .12,
        alpha_discovery = 0.05,
        model_names_l = NULL,
        model_names_d = NULL,
        restr_factor_d = NULL,
        ctrl_init = control_init(n_samp = 20, n_start_extra_classes = 50),
        ctrl_EM = control_EM(nstart_EM = 100)
      )
    )[3]
    
    cl_raedda_i <-
      fit_raedda_i$discovery_phase$Best$test$cl_after_trimming
    
    raedda_post_processing <- function(data = Y, fit_raedda) {
      # I compute the comp density and I reassigned those units that have comp density
      # higher than the quantile used to trim. Still experimental, need to think more about it
      # FIXME needs to be modified to take into account both transductive and inductive approaches
      # NOW IT ONLY WORKS WITH INDUCTIVE
      fit_raedda <- fit_raedda$discovery_phase
      test_density <- do.call(mclust::dens, c(
        list(
          data = data,
          logarithm = TRUE,
          modelName = fit_raedda$Best$model_name
        ),
        fit_raedda$Best
      ))
      
      test_comp_density <- do.call(mclust::cdens, c(
        list(
          data = data,
          logarithm = TRUE,
          modelName = fit_raedda$Best$model_name
        ),
        fit_raedda$Best
      ))
      
      test_assigned_comp_density <-
        sapply(1:nrow(data), function(x)
          test_comp_density[x, fit_raedda$Best$test$cl[x]])
      
      alpha_u_quantile <-
        quantile(test_density, probs = fit_raedda$Best$test$alpha_discovery)
      cl_map <-
        factor(fit_raedda$Best$test$cl, levels = c(levels(fit_raedda$Best$test$cl), "0"))
      cl_map[test_assigned_comp_density < alpha_u_quantile] <- 0
      cl_map
    }
    
    cl_raedda_i_map <-
      raedda_post_processing(data = Y, fit_raedda = fit_raedda_i)
    
    # METRIC COLLECTION -------------------------------------------------------
    
    # % hidden groups detected
    obs_hidden_groups <-
      which(cl_test > 3)
    
    # adaptBNP
    
    # novelty_adapt_BNP <- which(cl_adapt_BNP == 0)
    novelty_adapt_BNP <-
      purrr::map(fit_adapt_BNP_list, ~ which(.x$cl_alpha == 0))
    
    # prop_hidden_groups_adaptBNP <-
    #   mean(novelty_adapt_BNP %in% obs_hidden_groups)
    prop_hidden_groups_adaptBNP <-
      map_dbl(novelty_adapt_BNP, ~ mean(.x %in% obs_hidden_groups))
    
    # amda
    novelty_adma_i <- which(fit_amda_i$cls > 3)
    prop_hidden_groups_amda_i <-
      mean(novelty_adma_i %in% obs_hidden_groups)
    
    # raedda
    
    obs_novelty_raedda_i <- which(
      cl_raedda_i == "HIDDEN_GROUP_1" |
        cl_raedda_i == "HIDDEN_GROUP_2" |
        cl_raedda_i == "HIDDEN_GROUP_3" |
        cl_raedda_i == 0
    )
    
    obs_novelty_raedda_i_map <- which(
      cl_raedda_i_map == "HIDDEN_GROUP_1" |
        cl_raedda_i_map == "HIDDEN_GROUP_2" |
        cl_raedda_i_map == "HIDDEN_GROUP_3" |
        cl_raedda_i_map == 0
    )
    
    prop_hidden_groups_raedda_i <-
      mean(obs_novelty_raedda_i %in% obs_hidden_groups)
    
    prop_hidden_groups_raedda_i_map <-
      mean(obs_novelty_raedda_i_map %in% obs_hidden_groups)
    
    # METRIC 1: precision
    precision_hidden_group <-
      c(
        prop_hidden_groups_adaptBNP,
        prop_hidden_groups_raedda_i,
        prop_hidden_groups_raedda_i_map,
        prop_hidden_groups_amda_i
      )
    
    # ARI
    cl_hat <- vector("list", length = 7)
    cl_hat[1:4] <-
      purrr::map(fit_adapt_BNP_list,  ~ .x$cl_beta)
    
    cl_hat[[5]]  <- cl_raedda_i
    cl_hat[[6]]  <- cl_raedda_i_map
    cl_hat[[7]]  <- fit_amda_i$cls
    
    
    # METRIC 2: Adjusted Rand Index
    ARI <-
      map_dbl(cl_hat, ~
                mclust::adjustedRandIndex(x = cl_test, y = .x))
    
    # Accuracy on the training group
    # (for these classes I have the label, I can compute the accuracy)
    
    # METRIC 3: Accuracy on the known subset
    
    accuracy_on_the_known_subset <-
      map_dbl(purrr::map(cl_hat, .f = ~ .x[cl_test < 4]), ~
                mean(cl_test[cl_test < 4] == .x))
    
    # # ARI on the novelties only
    #
    # ARI_on_the_novelty_subset <-
    #   map_dbl(purrr::map(cl_hat, .f = ~ .x[cl_test >= 4]),
    #           ~
    #             mclust::adjustedRandIndex(x = cl_test[cl_test >= 4], y = .x))
    
    # Results collection ------------------------------------------------------
    
    model_names <-
      c(
        "MCD_1_lambda_tr_10",
        "MCD_75_lambda_tr_10",
        "MCD_1_lambda_tr_1000",
        "MCD_75_lambda_tr_1000",
        "RAEDDA",
        "RAEDDA_with_MAP",
        "AMDA"
      )
    
    comp_times <-
      c(map_dbl(fit_adapt_BNP_list,  ~ .x$comp_time),
        t_raedda_i,
        t_amda_i)
    
    names(comp_times) <- model_names[-6]
    names(ARI) <- names(accuracy_on_the_known_subset) <-
      names(precision_hidden_group) <-
      model_names
    
    res <-
      list(
        ARI = ARI,
        # ARI_on_the_novelty_subset = ARI_on_the_novelty_subset,
        accuracy_on_the_known_subset = accuracy_on_the_known_subset,
        precision_hidden_group = precision_hidden_group,
        comp_times = comp_times,
        cl = cbind(
          set_names(map_dfc(
            fit_adapt_BNP_list, .f = ~ cbind(.x$cl_alpha, .x$cl_beta)
          ), model_names[1:4]),
          cl_raedda = cl_raedda_i,
          cl_raedda_map = cl_raedda_i_map,
          cl_amdai = fit_amda_i$cls,
          cl_true = cl_test
        ),
        Y_test = bind_cols(Y, set_names(
          map_dfc(fit_adapt_BNP_list, ~ .x$a_posteriori_prob_novelty),
          model_names[1:4]
        ))
      ) # a posteriori classification
    res
  }

i <- NULL # dummy to trick R CMD check
# NSIM <- 100
MC_SIM <- 100

#THE FOLLOWING LINE TRIGGERS THE SIMULATION STUDY

clustvarsel::startParallel(TRUE)


for (perc_novelty in c("small", "not_small")) {
  for (lb in c(TRUE, FALSE)) {
    out <- foreach(
      i = 1:MC_SIM,
      .packages = c("Rcpp", "RcppArmadillo","mvnfast", "purrr", "mclust", "brand", "adaptDA", "raedda")
    ) %dopar% {
      set.seed(i*10) # reproducibility
      out <- sim_study_adaptBNP(label_noise = lb,
                                novelty_size = perc_novelty)
    }
    
    saveRDS(
      out,
      file = paste0(
        "Codice_review/multivariate_case/results/REVIEWsim_study_adaptDA_BNP_label_noise_",
        lb,
        "_novelty_",
        perc_novelty,
        ".Rds"
      )
    )
    
    cat(
      paste0(
        "sim_study_adaptDA_BNP_label_noise_",
        lb,
        "_novelty_",
        perc_novelty,
        ".Rds\n"
      )
    )
  }
}


# Load results ------------------------------------------------------------

theme_set(theme_bw())

collected_results_names <-
  str_subset(list.files("Codice_review/multivariate_case/results/"),pattern = "^RE") %>%
  str_split(pattern = "_",
            n = 5,
            simplify = TRUE) %>%
  magrittr::extract(, 5) %>%
  str_remove(pattern = "\\.Rds")

sim_results <-
  lapply(str_subset(list.files("Codice_review/multivariate_case/results/"),pattern = "^RE"), function(file)
    readRDS(paste0("Codice_review/multivariate_case/results/", file)))

names(sim_results) <- collected_results_names

sim_names_split <-
  str_split(
    collected_results_names,
    pattern = "_",
    n = 4,
    simplify = TRUE
  )

sim_names_split[, 3] <- str_to_title(sim_names_split[, 3])
sim_names_split[, 4] <-
  rep(c("Novelty~size==~Not~small", "Novelty~size==~Small"), 2)

metrics_name <- # Metrics I want to extract from simulations
  c("ARI",
    # "ARI_on_the_novelty_subset",
    "accuracy_on_the_known_subset",
    "precision_hidden_group")



tidy_BNP_name <- function(text) {
  values_extraction <-
    str_extract_all(text, pattern = "\\d+", simplify = TRUE)
  if (values_extraction[, 1] == "75") {
    values_extraction[, 1] <- "0.75"
  }
  paste0(
    "Brand~",
    "(",
    "list(eta[MCD]==~",
    values_extraction[, 1],
    ",~lambda[Tr]==",
    values_extraction[, 2],
    ")",
    ")"
  )
}


metric_collector <- function(n_sim_scenario, a_post_prob = FALSE) {
  if (a_post_prob == TRUE) {
    sim_metric_subset <-
      map_dfr(sim_results[[n_sim_scenario]], ~ .x$Y_test) %>%
      pivot_longer(cols = starts_with("MCD"),
                   names_to = "method",
                   values_to = "a_posteriori_prob_novelty")
    
  } else {
    sim_metric_subset <-
      purrr::map(sim_results[[n_sim_scenario]], ~ .x[metrics_name]) %>%
      unlist %>%
      enframe() %>%
      separate(col = name,
               into = c("metric", "method"),
               sep = "\\.")
  }
  
  sim_metric_subset  %>%
    mutate(
      method = fct_relevel(factor(
        map_chr(
          sim_metric_subset$method,
          .f = ~ ifelse(
            str_detect(.x, pattern = "^MCD"),
            yes = tidy_BNP_name(.x),
            no = .x
          )
        )
      ), "AMDA", after = Inf),
      label_noise = paste0("list(Label~Noise==~", sim_names_split[n_sim_scenario, 3], ")"),
      novelty_size = sim_names_split[n_sim_scenario, 4]
    )
  
}

df_metric <-
  map_dfr(.x = 1:length(sim_results),
          metric_collector,
          a_post_prob = FALSE)

df_metric <- df_metric %>% 
  filter(method!="RAEDDA") %>% # I keep only the RAEDDA with MAP and I call it RAEDDA
  dplyr::mutate(method=fct_recode(.f = method, "RAEDDA"="RAEDDA_with_MAP"))

boxplot_creator <- function(df, selected_metric, y_ticks = TRUE) {
  metric_name_4_plot <-
    switch (
      selected_metric,
      "ARI" = "Adjusted Rand Index",
      "ARI_on_the_novelty_subset" = "ARI on the novelty subset",
      "accuracy_on_the_known_subset" = "Accuracy on the observed classes subset",
      "perc_hidden_group" = "% Novelty Detected",
      "precision_hidden_group" = "Novelty predictive value"
    )
  
  df %>%
    filter(metric == selected_metric) %>%
    ggplot() +
    geom_boxplot(aes(y = method, x = value)) +
    labs(x = metric_name_4_plot, y = "") +
    facet_grid(label_noise ~ novelty_size,
               scales = "free",
               labeller = label_parsed) +
    scale_y_discrete(labels = ifelse(y_ticks, scales::label_parse(), list(NULL)))
  
  # gg_p <- df_metric %>%
  #   filter(metric == "ARI", novelty_size == "novelty_not_small") %>%
  #   ggplot() +
  #   labs(x = "", y = selected_metric) +
  #   geom_boxplot(aes(method, value)) +
  #   facet_grid(label_noise ~ alpha_MCD,
  #              scales = "free",
  #              labeller = label_parsed)
  # gg_p + gg_p_small
}

boxplot_creator(df = df_metric,
                selected_metric = "ARI",
                y_ticks = TRUE)

boxplot_creator(df = df_metric, selected_metric = "precision_hidden_group")

boxplot_creator(df = df_metric, selected_metric = "accuracy_on_the_known_subset")

walk(
  metrics_name,
  .f = ~ ggsave(
    filename = paste0(
      "~/Google Drive/university/research/publications/BRAND/after_review/figures/boxplot_",
      # change it accordingly
      .x,
      "_REVIEW.pdf"
    ),
    plot = boxplot_creator(df = df_metric, selected_metric = .x),
    width = 6.97,
    height = 3.83
  )
)

# Hex plot PPN ------------------------------------------------------------


hex_df <- # takes ~ 1 min, probably it is shitty implemented
  map_dfr(.x = 1:length(sim_results),
          metric_collector,
          a_post_prob = TRUE)
hex_df <- hex_df %>% 
  filter(!(V1>8 & V2>2))
ggplot(
  filter(hex_df, novelty_size == "Novelty~size==~Not~small"),
  mapping = aes(V1, V2, z = a_posteriori_prob_novelty)
) +
  stat_summary_hex(
    fun = function(x)
      mean(x), bins=50
  ) +
  facet_grid(label_noise ~ method, labeller = label_parsed) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(x = "", y = "", fill = quote(hat("PPN")))

ggsave(
  filename = paste0(
    "~/Google Drive/university/research/publications/BRAND/after_review/figures/hex_plot",
    # change it accordingly
    "_REVIEW.pdf"
  ),
  plot = last_plot(),
  width = 13,
  height = 6.39
)

# Sample scenario plot ----------------------------------------------------

novelty_size <- "not_small"
label_noise <- TRUE

M_vec <- switch(
  novelty_size,
  "small" = list(350, 250, 250, 49, 50, 50, 1),
  "not_small" = list(200, 200, 250, 90, 100, 100, 10)
)

cl_test <- rep(1:7, M_vec)

cl_train_label_noise <- cl_train
if (label_noise == TRUE) {
  obs_LB <-
    sample(which(cl_train_label_noise %in% c(2, 3)), size = 100)
  cl_train_label_noise[obs_LB] <-
    ifelse(test = cl_train_label_noise[obs_LB] == 2, 3, 2) # add label noise
}
X <-
  purrr::map_dfr(1:G,  ~ as.data.frame(mvnfast::rmvn(
    n = N_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]]
  )))

Y <-
  purrr::map_dfr(1:H,  ~ as.data.frame(mvnfast::rmvn(
    n = M_vec[[.x]], mu = MU[[.x]], sigma = SIGMA[[.x]]
  )))

plot(X, col = cl_train_label_noise)
plot(Y, col = cl_test)

sample_plot_df <-
  mutate(X, cl = cl_train_label_noise, set = "Training set") %>%
  bind_rows(mutate(Y, cl = cl_test, set = "Test set")) %>%
  mutate(cl = factor(cl), set = fct_rev(set))

ggplot(sample_plot_df) +
  geom_point(aes(V1, V2, col = cl, shape = cl),
             size = .7,
             alpha = .6) +
  facet_grid( ~ set) +
  scale_shape_manual(values = c(
    "1" = 1,
    "2" = 2,
    "3" = 3,
    "4" = 4,
    "5" = 5,
    "6" = 6,
    "7" = 7
  )) +
  scale_color_brewer(palette = "Dark2") +
  # coord_fixed() +
  labs(x = "",
       y = "",
       col = "Class",
       shape = "Class") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 2)))
ggsave("example_sim_study.pdf", width = 7.96, height = 3.76)


# Summary tables creation -------------------------------------------------

table_creator <-
  function(lb = "list(Label~Noise==~False)", ns = "Novelty~size==~Not~small") {
    df_4_table_avg <- df_metric %>%
      filter(label_noise == lb,
             novelty_size == ns,
             metric != "ARI_on_the_novelty_subset") %>%
      group_by(metric, method, label_noise, novelty_size) %>%
      summarise(avg = mean(value), sd = sd(value)) %>%
      ungroup %>%
      dplyr::select(metric, method, avg) %>%
      pivot_wider(names_from = metric, values_from = avg) %>%
      arrange(desc(method)) %>%
      dplyr::select(2:4)
    
    df_4_table_sd <- df_metric %>%
      filter(label_noise == lb,
             novelty_size == ns,
             metric != "ARI_on_the_novelty_subset") %>%
      group_by(metric, method, label_noise, novelty_size) %>%
      summarise(avg = mean(value), sd = sd(value)) %>%
      ungroup %>%
      dplyr::select(metric, method, sd) %>%
      pivot_wider(names_from = metric, values_from = sd) %>%
      arrange(desc(method)) %>%
      dplyr::select(2:4)
    
    summary_table <-
      matrix(ncol = (n_distinct(df_metric$metric) ) ,
             nrow = n_distinct(df_metric$method) * 2)
    
    summary_table[seq(from = 1,
                      to = nrow(summary_table),
                      by = 2), ] <- as.matrix(round(df_4_table_avg, 3))
    summary_table[seq(from = 2,
                      to = nrow(summary_table),
                      by = 2), ] <- paste0("(", round(as.matrix(df_4_table_sd), 3), ")")
    summary_table
  }


table_cases <- purrr::cross2(unique(df_metric$label_noise),
                             unique(df_metric$novelty_size))

tables_list <-
  purrr::map(.x = table_cases,
             ~ table_creator(lb = .x[[1]], ns = .x[[2]]))

tb_final <-
  (rbind(
    cbind(tables_list[[1]], tables_list[[2]]),
    cbind(tables_list[[3]], tables_list[[4]])
  ))

row_nm_tb_final <-
  paste0(
    "$",
    c(
      "AMDA",
      "1",
      "RAEDDA",
      "2",
      "BRAND(\\alpha_{MCD}=1,\\lambda_{Tr}=1000)",
      "3",
      "BRAND(\\alpha_{MCD}=1,\\lambda_{Tr}=10)",
      "4",
      "BRAND(\\alpha_{MCD}=0.75,\\lambda_{Tr}=1000)",
      "5",
      "BRAND(\\alpha_{MCD}=0.75,\\lambda_{Tr}=10)",
      "6"
    ),
    "$"
  )

col_nm_tb_final <- c("Accuracy on the known subset", "ARI",
                     "\\% Novelty detected")
colnames(tb_final) <- rep(col_nm_tb_final, 2)
rownames(tb_final) <- rep(row_nm_tb_final, 2)
library(xtable)

print.xtable(xtable(tb_final[1:12, ]),
             sanitize.rownames.function = identity)

print.xtable(xtable(tb_final[13:24, ]),
             sanitize.rownames.function = identity)
