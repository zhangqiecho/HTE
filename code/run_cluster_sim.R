# update Feb 2024: 
# 1) update the truth file
# 2) add the causal forest without super learner

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# run the function of cluster_sim using foreach with doParallel in HPC

packages <- c("data.table","tidyverse","SuperLearner","boot",
              "here","mvtnorm","foreach","doParallel","broom",
              "ggplot2","grf", "lmtest", "sandwich", "glmnet") # added a package

for (package in packages) {
  library(package, character.only=T)
}

func.packages <- c("data.table","tidyverse","SuperLearner","boot",
                   "here","mvtnorm","broom", "grf", "lmtest", "sandwich", "glmnet")

c_number <- as.numeric(args[1])
cat(paste("Number of Confounders for Each Sim:", c_number), "\n")
sample_size <- as.numeric(args[2])
cat(paste("Sample Size for Each Sim:", sample_size), "\n")
number_sims <- as.numeric(args[3])
cat(paste("Number of Simulations:", number_sims, "\n"))

# read in the seed data
seed_data <- read_csv(
  here("data", "random_seed_values.csv")
)

random_seeds <- seed_data$random_seeds

# read in the truth file
truth_data <- read_csv(
  here("data", "truth_data_simulation_all_update.csv")
)

# use c_number to extract the parameters needed from the truth data
truth_data <- truth_data %>%
  filter(c_dim == c_number) %>%
  mutate(group = row_number())

# param_dat1 <- list(param1 = 1.29,
#                   param2 = c(1.05, 1.29, 1.6))
# 
# data_grid1 <- expand.grid(param_dat1, KEEP.OUT.ATTRS = FALSE)
# 
# param_dat2 <- list(param1 = 1,
#                    param2 = c(.815, 1, 1.24))
# 
# data_grid2 <- expand.grid(param_dat2, KEEP.OUT.ATTRS = FALSE)
# 
# param_dat3 <- list(param1 = 0.765,
#                    param2 = c(.625, .765, .945))
# 
# data_grid3 <- expand.grid(param_dat3, KEEP.OUT.ATTRS = FALSE)
# 
# data_grid <- rbind(data_grid1, data_grid2, data_grid3)

source("cluster_sim_func.R")

# save the random seed file
# set.seed(123)
# random_seeds <- runif(1000, min = 1, max = 249999)
# random_seeds <- unique(round(random_seeds, digits = 0))
# write_csv(data.frame(random_seeds), here("data","random_seed_values.csv"))

#parallel::detectCores()

#n.cores <- parallel::detectCores() - 2 # the detected number of cores might not be active on HPC

n.cores = 20

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK" # qz: FORK doesn't work on windows,doesn't work for Rollins cluster
)

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()


results <- foreach(i = 1:nrow(truth_data), .packages = func.packages) %dopar% {
  
  cat(paste("for loop iteration",i), "\n")
  
  result_list <- lapply(random_seeds[1:number_sims], function(x) cluster_sim(seed_value = x,
                                                                        number_sims = number_sims,
                                                                        sample_size = sample_size, 
                                                                        c_number = c_number,
                                                                        param1 = log(truth_data[i,]$param1), 
                                                                        param2 = log(truth_data[i,]$param2),
                                                                        group = i))
  
}

parallel::stopCluster(cl = my.cluster)

# extract the blp output and result table to separate files
result_all <- blp_cf_all <- blp_cf_sl_all <- blp_drl_all <- NULL
for (i in seq(nrow(truth_data))) {
  for (j in seq(number_sims)) {
    result_all <- rbind(result_all, results[[i]][[j]][[1]])
    blp_cf_all <- rbind(blp_cf_all, results[[i]][[j]][[2]])
    blp_cf_sl_all <- rbind(blp_cf_sl_all, results[[i]][[j]][[3]])
    blp_drl_all <- rbind(blp_drl_all, results[[i]][[j]][[4]])
  }
}


# combine the truth with the result
truth <- gather(truth_data, estimand, truth, risk_diff, risk_diff_m0, risk_diff_m1) %>%
  mutate(type = ifelse(estimand == "risk_diff", "ate", ifelse(estimand == "risk_diff_m0", "_m0", "_m1"))) %>%
  select(group, truth, type)
result_all_truth <- result_all %>%
  mutate(type = substr(estimator, nchar(estimator)-2, nchar(estimator))) %>%
  left_join(truth, by = c("type", "group"))

# calculate the bias, MSE and coverage probability compared with truth
summary_truth <- result_all_truth %>%
  mutate(bias = estimate - truth, 
         coverage = ifelse(estimate - 1.96*std.err <= truth & truth <= estimate + 1.96*std.err, 1, 0)) %>%
  group_by(group, estimator) %>%
  mutate(p_exposure_mean = mean(p_exposure), bias_mean = mean(bias), MSE = mean((bias)^2), coverage_prob = mean(coverage)) %>%
  select(group, estimator, sample_size, confounder_number, or_exposure, or_interaction, p_exposure_mean, truth, bias_mean, MSE, coverage_prob) %>%
  unique()
 
# combine the oracle truth with the result
oracle <- result_all_truth %>%
  filter(substr(estimator, 1, 6) == "oracle") %>%
  select(group,  seed_number, estimator, estimate, type) %>% 
  dplyr::rename(oracle = estimate) %>%
  select(group, seed_number, oracle, type)
result_all_oracle <- result_all %>%
  filter(substr(estimator, 1, 6) != "oracle") %>%
  mutate(type = substr(estimator, nchar(estimator)-2, nchar(estimator))) %>%
  left_join(oracle, by = c("type", "group", "seed_number"))
# calculate the bias, MSE and coverage probability compared with oracle
summary_oracle <- result_all_oracle %>%
  mutate(bias = estimate - oracle, 
         coverage = ifelse(estimate - 1.96*std.err <= oracle & oracle <= estimate + 1.96*std.err, 1, 0)) %>%
  group_by(group, estimator) %>%
  mutate(p_exposure_mean = mean(p_exposure), oracle_mean = mean(oracle), bias_mean = mean(bias), MSE = mean((bias)^2), coverage_prob = mean(coverage))%>%
  select(group, estimator, sample_size, confounder_number, or_exposure, or_interaction, p_exposure_mean, oracle_mean, bias_mean, MSE, coverage_prob) %>%
  unique()


# save the files
## all results
write_csv(result_all, here("output", paste0("All_results_c", c_number, "_n", sample_size, "_sims", number_sims, "_", Sys.Date(), ".csv")))
## summary results
write_csv(summary_truth, here("output", paste0("Summary_truth_c", c_number, "_n", sample_size, "_sims", number_sims, "_", Sys.Date(), ".csv")))
write_csv(summary_oracle, here("output", paste0("Summary_oracle_c", c_number, "_n", sample_size, "_sims", number_sims, "_", Sys.Date(), ".csv")))
## blp results
write_csv(blp_cf_all, here("output", paste0("Blp_cf_c", c_number, "_n", sample_size, "_sims", number_sims, "_", Sys.Date(), ".csv")))
write_csv(blp_cf_sl_all, here("output", paste0("Blp_cf_sl_c", c_number, "_n", sample_size, "_sims", number_sims, "_", Sys.Date(), ".csv")))
write_csv(blp_drl_all, here("output", paste0("Blp_drl_c", c_number, "_n", sample_size, "_sims", number_sims, "_", Sys.Date(), ".csv")))

