# use the code in figure.R to produce plots for cf comparison

# Part 1: compare causal forest performances  using super learner with big sample size (n = 10000)
# step 1: compare different parameter options for cf with n = 10000, nsims = 100, confounder_number = 5-15
# step 2: create the table
# read in the output files using regular expression
pacman::p_load(tidyverse, ggplot2, looplot)
               
plot_data <- function(metric){
  setwd("H:/RA/HTE/output/sim100_cf/min.node.size_30")
  results_min.node <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv) %>%
    filter(grepl("^causal_forest", estimator)) %>%
    mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
                             ifelse(grepl("_cate_m0$", estimator), "cate_m0",
                                    "cate_m1"))) %>%
    select(group, confounder_number, estimate, all_of(metric)) %>%
    dplyr::rename(min.node.size30 = all_of(metric))
  
  
  setwd("H:/RA/HTE/output/sim100_cf/num.trees_10000")
  results_num.tree <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv) %>%
    filter(grepl("^causal_forest", estimator)) %>%
    mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
                             ifelse(grepl("_cate_m0$", estimator), "cate_m0",
                                    "cate_m1"))) %>%
    select(group, confounder_number, estimate, all_of(metric)) %>%
    dplyr::rename(num.tree.10k = all_of(metric))
  
  setwd("H:/RA/HTE/output/sim100_cf/tune.parameter=size&prune")
  results_tune.node.prune <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv) %>%
    filter(grepl("^causal_forest", estimator)) %>%
    mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
                             ifelse(grepl("_cate_m0$", estimator), "cate_m0",
                                    "cate_m1"))) %>%
    select(group, confounder_number, estimate, all_of(metric)) %>%
    dplyr::rename(tune.node.prune = all_of(metric))
  
  
  results <- results_min.node %>%
    left_join(results_num.tree, by = c("group","confounder_number", "estimate")) %>%
    left_join(results_tune.node.prune, by = c("group", "confounder_number", "estimate"))
  
  # read in the truth file
  truth_data <- read_csv(
    here::here("data", "truth_data_simulation_all.csv")
  )
  
  # get the truth data for each group
  truth <- truth_data %>%
    filter(c_dim == 15) %>%
    mutate(group = row_number(), RD_m0 = round(risk_diff_m0, 3), RD_m1 = round(risk_diff_m1, 3)) %>%
    select(group, RD_m0, RD_m1)
  
  # merge the truth data with the result
  plot_result <- results %>%
    left_join(truth, by = "group") %>%
    filter(confounder_number != 15)%>% # we lack some results when confounder_num = 15
  select(-group)
  
  return(plot_result)
}

plot_result <- plot_data(metric = "coverage_prob")
p = nested_loop_plot(resdf = plot_result,
                     x = "confounder_number", steps = c("estimate","RD_m1","RD_m0"),
                     methods = c("min.node.size30", "num.tree.10k", "tune.node.prune"),
                     steps_y_base = -0.1, steps_y_height = 0.05, steps_y_shift = 0.05,
                     x_name = "Confounder number", y_name = "Coverage",
                     y_breaks = seq(0, 1, by = 0.25),
                     spu_x_shift = 10,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                     hline_intercept = 0,
                     y_expand_add = c(0.1, NULL),
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 5)
                       )
                     )) +
  labs(title = "Causal forest performances using super learner with big sample size",
       subtitle = "Sample size = 10,000, simulation number = 100")
ggsave("cf_n10K_coverage.png", p, path = here::here("output"), width = 12, height = 7)


plot_result <-plot_data(metric = "bias_mean")
p = nested_loop_plot(resdf = plot_result,
                     x = "confounder_number", steps = c("estimate","RD_m1","RD_m0"),
                     steps_y_base = -0.045, steps_y_height = 0.003, steps_y_shift = 0.003,
                     x_name = "Confounder number", y_name = "Bias",
                     y_breaks = seq(-0.02, 0.02, by = 0.01),
                     spu_x_shift = 10,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                     hline_intercept = -0.04,
                     y_expand_add = c(0.005, NULL),
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 5)
                       )
                     ))+
  labs(title = "Causal forest performances using super learner with big sample size",
       subtitle = "Sample size = 10,000, simulation number = 100")

ggsave("cf_n10k_bias.png", p, path = here::here("output"), width = 12, height = 7)

plot_result <-plot_data(metric = "MSE")
p = nested_loop_plot(resdf = plot_result,
                     x = "confounder_number", steps = c("estimate","RD_m1","RD_m0"),
                     steps_y_base = -0.0002, steps_y_height = 0.0001, steps_y_shift = 0.0002,
                     x_name = "Confounder number", y_name = "MSE",
                     y_breaks = seq(0, 0.005, by = 0.001),
                     spu_x_shift = 10,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                     hline_intercept = 0,
                     y_expand_add = c(0.0001, NULL),
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 5)
                       )
                     ))+
  labs(title = "Causal forest performances using super learner with big sample size",
       subtitle = "Sample size = 10,000, simulation number = 100")
ggsave("cf_n10k_MSE.png", p, path = here::here("output"), width = 12, height = 7)

# Part 1: compare causal forest performances using super learner with small sample size (n = 1000): tune vs. untuned

plot_data2 <- function(metric){
  setwd("H:/RA/HTE/output/experiment output/super learner/tuned")
  results_tuned <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv) %>%
    filter(grepl("^causal_forest", estimator)) %>%
    mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
                             ifelse(grepl("_cate_m0$", estimator), "cate_m0",
                                    "cate_m1"))) %>%
    select(group, estimate, all_of(metric)) %>%
    dplyr::rename(causal_forest_tuned = all_of(metric))
  
  
  setwd("H:/RA/HTE/output/experiment output/super learner/untuned")
  results_untuned <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv) %>%
    filter(grepl("^causal_forest", estimator)) %>%
    mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
                             ifelse(grepl("_cate_m0$", estimator), "cate_m0",
                                    "cate_m1"))) %>%
    select(group, estimate, all_of(metric)) %>%
    dplyr::rename(causal_forest_untuned = all_of(metric))
  
  
  results <- results_untuned %>%
    left_join(results_tuned, by = c("group", "estimate"))
  
  # read in the truth file
  truth_data <- read_csv(
    here::here("data", "truth_data_simulation_all.csv")
  )
  
  # get the truth data for each group
  truth <- truth_data %>%
    filter(c_dim == 15) %>%
    mutate(group = row_number(), RD_m0 = round(risk_diff_m0, 3), RD_m1 = round(risk_diff_m1, 3)) %>%
    select(group, RD_m0, RD_m1)
  
  # merge the truth data with the result
  plot_result <- results %>%
    left_join(truth, by = "group") %>%
    select(-group)
  
  return(plot_result)
}

plot_result <- plot_data2(metric = "coverage_prob")
p = nested_loop_plot(resdf = plot_result,
                     x = "RD_m0", steps = c("estimate", "RD_m1"),
                     steps_y_base = -0.1, steps_y_height = 0.05, steps_y_shift = 0.05,
                     x_name = "Estimate", y_name = "Coverage",
                     y_breaks = seq(0, 1, by = 0.25),
                     spu_x_shift = 0.025,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                     hline_intercept = 0,
                     y_expand_add = c(0.1, NULL),
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 5)
                       )
                     )) +
  labs(title = "Causal forest performances using super learner: tuned vs. untuned",
       subtitle = "Sample size = 1,000, confounder number = 15, simulation number = 100")
ggsave("cf_n1000_coverage.png", p, path = here::here("output"), width = 12, height = 7)


plot_result <-plot_data2(metric = "bias_mean")
p = nested_loop_plot(resdf = plot_result,
                     x = "RD_m0", steps = c("estimate", "RD_m1"),
                     steps_y_base = -0.045, steps_y_height = 0.003, steps_y_shift = 0.003,
                     x_name = "Confounder number", y_name = "Bias",
                     y_breaks = seq(-0.02, 0.02, by = 0.01),
                     spu_x_shift =0.025,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                     hline_intercept = -0.04,
                     y_expand_add = c(0.005, NULL),
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 5)
                       )
                     )) +
  labs(title = "Causal forest performances using super learner: tuned vs. untuned",
       subtitle = "Sample size = 1,000, confounder number = 15, simulation number = 100")

ggsave("cf_n1000_bias.png", p, path = here::here("output"), width = 12, height = 7)

plot_result <-plot_data2(metric = "MSE")
p = nested_loop_plot(resdf = plot_result,
                     x = "RD_m0", steps = c("estimate", "RD_m1"),
                     steps_y_base = -0.0002, steps_y_height = 0.0001, steps_y_shift = 0.0002,
                     x_name = "Confounder number", y_name = "MSE",
                     y_breaks = seq(0, 0.005, by = 0.001),
                     spu_x_shift =0.025,
                     steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                     hline_intercept = 0,
                     y_expand_add = c(0.0001, NULL),
                     post_processing = list(
                       add_custom_theme = list(
                         axis.text.x = element_text(angle = -90,
                                                    vjust = 0.5,
                                                    size = 5)
                       )
                     )) +
  labs(title = "Causal forest performances using super learner: tuned vs. untuned",
       subtitle = "Sample size = 1,000, confounder number = 15, simulation number = 100")
ggsave("cf_n1000_MSE.png", p, path = here::here("output"), width = 12, height = 7)