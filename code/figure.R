# purpose: create nested loop plots for the HTE simulation project
# by Echo Jan 25, 2024
# plan: 
# step 1:methods: compare DR learner and causal forest, both with super learner
# step 2:create the table for simulation results with confounder number, sample size and true effects
#        for bias, MSE and coverage 
# step 3:use the code in the link https://matherealize.github.io/looplot_demo.html
#                                 https://matherealize.github.io/looplot_gallery.html


# install.packages("devtools")
# devtools::install_github("matherealize/looplot", force = TRUE)

pacman::p_load(tidyverse, ggplot2, looplot,reshape, data.table) #other package needed: here

nested_plot <- function(metric = "coverage_prob", y_name = "Coverage",estimate = "cate_m1"){
  # step 2: create the table
  # read in the output files using regular expression
  # setwd("H:/RA/HTE/output/sim1000_n10000")
  # results_big <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
  #   map_dfr(read_csv)
  # setwd("H:/RA/HTE/output/sim1000_n1000")
  # results_small <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
  #   map_dfr(read_csv)
  # results <- rbind(results_big, results_small)
  setwd("H:/RA/HTE/output/sim1000_n10000 April2024")
  results <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv)
  
  # read in the truth file
  truth_data <- read_csv(
    here::here("data", "truth_data_simulation_all_update.csv")
  )
  
  # get the truth data for each group
  truth <- truth_data %>%
    filter(c_dim == 15) %>%
    mutate(group = row_number(), RD_m0 = round(risk_diff_m0, 3), RD_m1 = round(risk_diff_m1, 3)) %>%
    select(group, RD_m0, RD_m1)
  
  # merge the truth data with the result
  plot_result <- results %>%
    mutate(truth = round(truth, 3)) %>%
    mutate(RMSE = sqrt(MSE)) %>%
    left_join(truth, by = "group")
  
  # double check if the truth is the same as we expected for each group
  #table(plot_result[(grep("_cate_m0$", plot_result$estimator)),]$truth) # not exactly same
  #table(plot_result[(grep("_cate_m1$", plot_result$estimator)),]$truth) # exactly same
  
  # subset data to only the targeted estimate (ate vs. cate) and evaluation metrics 
  plot_result2 <- plot_result[(grep(paste0(estimate, "$"), plot_result$estimator)),] %>%
    select(estimator, group, RD_m0, RD_m1, sample_size, confounder_number, all_of(metric)) %>%
    #filter(sample_size %in% c(500, 5000, 10000)) %>% # decrease the number of data points (sample size values) on x axis to make it look more clear
    mutate(sample_size = factor(sample_size))%>%
    filter(complete.cases(.)) # oracle doesn;t have value for coverage 
  
  # transform data from long to wide for plots
  plot_data <- data.table::dcast(setDT(plot_result2), RD_m0+RD_m1+sample_size+confounder_number~estimator, value.var = metric)
  
  # Define custom names for the steps (without underscores)
  plot_names <- c(
    "confounder_number" = "Confounder Number",
    "RD_m1" = "RD m1",
    "RD_m0" = "RD m0",
    "DRLearner_cate_m0" = "DR Learner",
    "causal_forest_SL_cate_m0" = "Causal Forest with SL",
    "causal_forest_cate_m0" = "Causal Forest",
    "gComp_cate_m0" = "g Computation",
    "oracle_cate_m0" = "Oracle",
    "DRLearner_cate_m1" = "DR Learner",
    "causal_forest_SL_cate_m1" = "Causal Forest with SL",
    "causal_forest_cate_m1" = "Causal Forest",
    "gComp_cate_m1" = "g Computation",
    "oracle_cate_m1" = "Oracle",
    "DRLearner_ate" = "DR Learner",
    "causal_forest_SL_ate" = "Causal Forest with SL",
    "causal_forest_ate" = "Causal Forest",
    "gComp_ate" = "g Computation",
    "oracle_ate" = "Oracle"
  )
  
  # Rename the columns in your data before plotting
  plot_data_renamed <- plot_data
  for (old_name in names(plot_names)) {
    if (old_name %in% names(plot_data_renamed)) {
      names(plot_data_renamed)[names(plot_data_renamed) == old_name] <- plot_names[old_name]
    }
  }
  
  # plot_all <- plot_result %>%
  #   mutate(estimate = ifelse(grepl("_ate$", estimator), "ate",
  #          ifelse(grepl("_cate_m0$", estimator), "cate_m0",
  #          "cate_m1"))) 
  
  
  # step 3:
  if(metric=="coverage_prob"){
    p = nested_loop_plot(resdf = plot_data_renamed,
                         x = "sample_size", steps = c("Confounder Number","RD m1","RD m0"),
                         steps_y_base = -0.1, steps_y_height = 0.05, steps_y_shift = 0.08,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(0, 1, by = 0.25),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 6,
                         hline_intercept = 0,
                         y_expand_add = c(0.1, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3","darkorange"), # https://r-graph-gallery.com/ggplot2-color.html
                         line_linetypes = c(1, 2, 3,4),
                         point_shapes = c(1,2,3,4),
                         base_size = 25,
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 10)
                           )
                         ))
  }else if(metric=="bias_mean"){
    p = nested_loop_plot(resdf = plot_data_renamed,
                         x = "sample_size", steps = c("Confounder Number","RD m1","RD m0"),
                         steps_y_base = -0.06, steps_y_height = 0.008, steps_y_shift = 0.015,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(-0.02, 0.10, by = 0.02),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 6,
                         hline_intercept = -0.04,
                         y_expand_add = c(0.02, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3","darkorange", "yellow"),
                         line_linetypes = c(1,2,3,4,5),
                         point_shapes = c(1,2,3,4,5),
                         base_size = 25,
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 10)
                           )
                         ))
  } else{
    p = nested_loop_plot(resdf = plot_data_renamed,
                         x = "sample_size", steps = c("Confounder Number","RD m1","RD m0"),
                         steps_y_base = -0.01, steps_y_height = 0.005, steps_y_shift = 0.012,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(0, 0.1, by = 0.02),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 6,
                         hline_intercept = 0,
                         y_expand_add = c(0.01, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3","darkorange", "yellow"),
                         line_linetypes = c(1,2,3,4,5),
                         point_shapes = c(1,2,3,4,5),
                         base_size = 25,
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 10)
                           )
                         ))
  }
  
  
  # figure in a panel
  # p2 = nested_loop_plot(resdf = plot_data, 
  #                      x = "sample_size", steps = c("confounder_number"),
  #                      grid_rows = "RD_m1", grid_cols = "RD_m0",
  #                      steps_y_base = -0.1, steps_y_height = 0.1, 
  #                      x_name = "Sample Size", y_name = y_name, ylim = c(-0.5,1),
  #                      spu_x_shift = 1000,
  #                      steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
  #                      hline_intercept = 0, 
  #                      y_expand_add = c(0.2, NULL), 
  #                      post_processing = list(
  #                        add_custom_theme = list(
  #                          axis.text.x = element_text(angle = -90, 
  #                                                     vjust = 0.5, 
  #                                                     size = 5) 
  #                        )
  #                      ))
  # p3 = nested_loop_plot(resdf = plot_data,
  #                       x = "sample_size", steps = c("RD_m1","RD_m0"),
  #                       grid_cols = "confounder_number",
  #                       steps_y_base = -0.1, steps_y_height = 0.1,
  #                       x_name = "Sample Size", y_name = y_name,
  #                       spu_x_shift = 1000,
  #                       steps_values_annotate = TRUE, steps_annotation_size = 2.5,
  #                       hline_intercept = 0,
  #                       y_expand_add = c(0.2, NULL),
  #                       post_processing = list(
  #                         add_custom_theme = list(
  #                           axis.text.x = element_text(angle = -90,
  #                                                      vjust = 0.5,
  #                                                      size = 5)
  #                         )
  #                       ))
  ggsave(paste0(y_name, "_", estimate, ".png"), p, path = here::here("output"), width = 20, height = 8)
}

nested_plot(metric = "coverage_prob", y_name = "Coverage",estimate = "ate")
nested_plot(metric = "coverage_prob", y_name = "Coverage",estimate = "cate_m0")
nested_plot(metric = "coverage_prob", y_name = "Coverage",estimate = "cate_m1")

nested_plot(metric = "bias_mean", y_name = "Bias",estimate = "ate")
nested_plot(metric = "bias_mean", y_name = "Bias",estimate = "cate_m0")
nested_plot(metric = "bias_mean", y_name = "Bias",estimate = "cate_m1")

nested_plot(metric = "RMSE", y_name = "RMSE",estimate = "ate")
nested_plot(metric = "RMSE", y_name = "RMSE",estimate = "cate_m0")
nested_plot(metric = "RMSE", y_name = "RMSE",estimate = "cate_m1")

# update Oct 8, 2025: panel plot
panel_plot <- function(metric = "coverage_prob", y_name = "Coverage"){
  # step 2: create the table
  setwd("H:/RA/HTE/output/sim1000_n10000 April2024")
  results <-list.files(path = getwd(), pattern = "Summary_truth") %>% 
    map_dfr(read_csv)
  
  # read in the truth file
  truth_data <- read_csv(
    here::here("data", "truth_data_simulation_all_update.csv")
  )
  
  # get the truth data for each group
  truth <- truth_data %>%
    filter(c_dim == 15) %>%
    mutate(group = row_number(), RD_m0 = round(risk_diff_m0, 3), RD_m1 = round(risk_diff_m1, 3)) %>%
    select(group, RD_m0, RD_m1)
  
  # merge the truth data with the result
  plot_result <- results %>%
    mutate(truth = round(truth, 3)) %>%
    mutate(RMSE = sqrt(MSE)) %>%
    left_join(truth, by = "group")
  
  # double check if the truth is the same as we expected for each group
  #table(plot_result[(grep("_cate_m0$", plot_result$estimator)),]$truth) # not exactly same
  #table(plot_result[(grep("_cate_m1$", plot_result$estimator)),]$truth) # exactly same
  
  # subset data to only the targeted estimate (cate) and evaluation metrics 
  plot_result2 <- plot_result[(grep("_cate_", plot_result$estimator)),] %>%
    select(estimator, group, RD_m0, RD_m1, sample_size, confounder_number, all_of(metric)) %>%
    #filter(sample_size %in% c(500, 5000, 10000)) %>% # decrease the number of data points (sample size values) on x axis to make it look more clear
    mutate(sample_size = factor(sample_size))%>%
    mutate(RD = case_when(
      RD_m0 == 0.025 & RD_m1 == 0.025 ~ 0.025,
      RD_m0 == -0.025 & RD_m1 == -0.025 ~ -0.025,
      TRUE ~ NA_real_    # fallback for all other cases
    ))%>%
    mutate(M = as.integer(sub(".*cate_m", "", estimator))) %>%
    filter(complete.cases(.)) %>% # oracle doesn;t have value for coverage 
    mutate(estimator_base = sub("_m[01]$", "", estimator))
  
  # transform data from long to wide for plots
  plot_data <- data.table::dcast(setDT(plot_result2), RD+sample_size+confounder_number+M~estimator_base, value.var = metric)
  
  # Define custom names for the steps (without underscores)
  plot_names <- c(
    "confounder_number" = "Confounder Number",
    #"RD_m1" = "RD m1",
    #"RD_m0" = "RD m0",
    "DRLearner_cate" = "DR Learner",
    "causal_forest_SL_cate" = "Causal Forest with SL",
    "causal_forest_cate" = "Causal Forest",
    "gComp_cate" = "g Computation",
    "oracle_cate" = "Oracle"
  )
  
  # Rename the columns in your data before plotting
  plot_data_renamed <- plot_data
  for (old_name in names(plot_names)) {
    if (old_name %in% names(plot_data_renamed)) {
      names(plot_data_renamed)[names(plot_data_renamed) == old_name] <- plot_names[old_name]
    }
  }
  
  
  # step 3:
  if(metric=="coverage_prob"){
    p = nested_loop_plot(resdf = plot_data_renamed,
                         x = "sample_size", steps = c("Confounder Number"),
                         grid_rows = "M", grid_cols = "RD", 
                         steps_y_base = -0.2, steps_y_height = 0.05, steps_y_shift = 0.08,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(0, 1, by = 0.25),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 6,
                         hline_intercept = 0,
                         y_expand_add = c(0.2, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3","darkorange"), # https://r-graph-gallery.com/ggplot2-color.html
                         line_linetypes = c(1, 2, 3,4),
                         point_shapes = c(1,2,3,4),
                         base_size = 25,
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 10)
                           )
                         ))
  }else if(metric=="bias_mean"){
    p = nested_loop_plot(resdf = plot_data_renamed,
                         x = "sample_size", steps = c("Confounder Number"),
                         grid_rows = "M", grid_cols = "RD", 
                         steps_y_base = -0.08, steps_y_height = 0.01, steps_y_shift = 0.015,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(-0.02, 0.10, by = 0.02),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 6,
                         hline_intercept = -0.04,
                         y_expand_add = c(0.04, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3","darkorange", "yellow"),
                         line_linetypes = c(1,2,3,4,5),
                         point_shapes = c(1,2,3,4,5),
                         base_size = 25,
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 10)
                           )
                         ))
  } else{
    p = nested_loop_plot(resdf = plot_data_renamed,
                         x = "sample_size", steps = c("Confounder Number"),
                         grid_rows = "M", grid_cols = "RD", 
                         steps_y_base = -0.03, steps_y_height = 0.01, steps_y_shift = 0.012,
                         x_name = "Sample Size", y_name = y_name,
                         y_breaks = seq(0, 0.1, by = 0.02),
                         #spu_x_shift = 2000,
                         steps_values_annotate = TRUE, steps_annotation_size = 6,
                         hline_intercept = 0,
                         y_expand_add = c(0.04, NULL),
                         colors = c("blue4", "darkcyan", "chartreuse3","darkorange", "yellow"),
                         line_linetypes = c(1,2,3,4,5),
                         point_shapes = c(1,2,3,4,5),
                         base_size = 25,
                         post_processing = list(
                           add_custom_theme = list(
                             axis.text.x = element_text(angle = -90,
                                                        vjust = 0.5,
                                                        size = 10)
                           )
                         ))
  }
  
  
  # figure in a panel
  # p2 = nested_loop_plot(resdf = plot_data, 
  #                      x = "sample_size", steps = c("confounder_number"),
  #                      grid_rows = "RD_m1", grid_cols = "RD_m0",
  #                      steps_y_base = -0.1, steps_y_height = 0.1, 
  #                      x_name = "Sample Size", y_name = y_name, ylim = c(-0.5,1),
  #                      spu_x_shift = 1000,
  #                      steps_values_annotate = TRUE, steps_annotation_size = 2.5, 
  #                      hline_intercept = 0, 
  #                      y_expand_add = c(0.2, NULL), 
  #                      post_processing = list(
  #                        add_custom_theme = list(
  #                          axis.text.x = element_text(angle = -90, 
  #                                                     vjust = 0.5, 
  #                                                     size = 5) 
  #                        )
  #                      ))
  # p3 = nested_loop_plot(resdf = plot_data,
  #                       x = "sample_size", steps = c("RD_m1","RD_m0"),
  #                       grid_cols = "confounder_number",
  #                       steps_y_base = -0.1, steps_y_height = 0.1,
  #                       x_name = "Sample Size", y_name = y_name,
  #                       spu_x_shift = 1000,
  #                       steps_values_annotate = TRUE, steps_annotation_size = 2.5,
  #                       hline_intercept = 0,
  #                       y_expand_add = c(0.2, NULL),
  #                       post_processing = list(
  #                         add_custom_theme = list(
  #                           axis.text.x = element_text(angle = -90,
  #                                                      vjust = 0.5,
  #                                                      size = 5)
  #                         )
  #                       ))
  ggsave(paste0(y_name, "_panel.png"), p, path = here::here("output"), width = 20, height = 8)
}

panel_plot(metric = "coverage_prob", y_name = "Coverage")
panel_plot(metric = "bias_mean", y_name = "Bias")
panel_plot(metric = "RMSE", y_name = "RMSE")
