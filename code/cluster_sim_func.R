## FOCUS ON PROBLEM OF USING CF AND DR LEARNER TO FIND THE RIGHT MODIFIER
# function of the simulation
# including data generation and analysis
## update in Jan 2024: 
# 1) round down the continuous confounders
## update in Feb 2024:
# 1) update the calculation of oracle truth
# 2) add the causal forest without super learner

# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

cluster_sim <- function(seed_value, number_sims, sample_size, c_number, param1, param2, group){
  
  set.seed(seed_value)
  param1 <- exp(param1) 
  param2 <- exp(param2)
  n = sample_size
  p = c_number
  
  ## CONFOUNDERS
  sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
  # round down the continuous confounders to reduce the number of split points to search over
  c     <- round(mvtnorm::rmvnorm(n, mean=rep(0,p), sigma=sigma), 3) 
  
  # DESIGN MATRIX FOR C 
  dMat <- model.matrix(
    as.formula(
      paste("~(",
            paste("c[,",1:ncol(c),"]", collapse="+")
            ,")")
    )
  )
  
  beta <- rep(log(1.75),c_number) # outcome-confounder ORs:1.75
  
  theta <- c(-.5, rep(log(1.5),c_number)) #exposure-confounder ORs:1.5, Q: why intercept -0.5?
  pi_ <- expit(dMat%*%theta) # predicted exposure and modifier probability (they are the same)
  
  # EXPOSURE AND MODIFIER
  x    <- rbinom(n,1,pi_)
  m    <- rbinom(n,1,pi_)
  
  balancing_intercept <- -log((1/.1)-1) - # outcome mu of .1
    log(param1)*pi_ -       # offsetting exposure effect
    log(2)*pi_ +            # offsetting modifier effect: outcome-modifier OR: 2
    log(param2)*pi_*pi_   # offsetting interaction term
  
  # PREDICTED OUTCOME PROBABILITIES
  mu_y <- expit(balancing_intercept +
                  log(param1)*x+ log(2)*m - log(param2)*x*m + # this is the exposure and modifier effect
                  dMat[,-1]%*%beta # this is the confounder-outcome part
  )
  
  y <- rbinom(n, 1, mu_y)
  
  dat <- data.frame(y,m,c,x)
  names(dat) <- c("y","m",paste0("c",1:c_number),"x")
  
  # mu_exposure_modifier: Q: what if modifier is also a confounder?
  mu11 <- expit(balancing_intercept + 
                  log(param1)*1 + log(2)*1 - 
                  log(param2)*1*1 + dMat[,-1]%*%beta)
  
  mu10 <- expit(balancing_intercept + 
                  log(param1)*1 + log(2)*0 - 
                  log(param2)*1*0 + dMat[,-1]%*%beta)
  
  mu01 <- expit(balancing_intercept + 
                  log(param1)*0 + log(2)*1 - 
                  log(param2)*0*1 + dMat[,-1]%*%beta)
  
  mu00 <- expit(balancing_intercept + 
                  log(param1)*0 + log(2)*0 - 
                  log(param2)*0*0 + dMat[,-1]%*%beta)
  # oracle truth in the data generating mechanism
  rd_m0 <- mean(mu10 - mu00)
  rd_m1 <- mean(mu11 - mu01)
  # rd_ate <- mean(rd_m0*(1-pi_) + rd_m1*pi_)
  
  mu1 <-  expit(balancing_intercept +
                  log(param1)*1 + log(2)*pi_ -
                  log(param2)*1*pi_ + dMat[,-1]%*%beta)
  mu0 <-  expit(balancing_intercept +
                  log(param1)*0 + log(2)*pi_ -
                  log(param2)*0*pi_ + dMat[,-1]%*%beta)
  rd_ate <- mean(mu1-mu0) # udpated calculation of oracle ATE
  
  
  
  ## DATA GENERATION FINISHED; DATA ANALYSIS START
  
  # ORACLE MODEL, G COMP
  mod_form <- as.formula(paste("y ~ x + m + x*m +", paste0("c", 1:c_number, collapse = "+")))
  mod_true <- glm(mod_form, data = dat, family = binomial(link = "logit"))
  
  # c_matrix1 <- c(0, (as.matrix(attr(summary(mod_true)$terms,"term.labels")) == c("x")))
  # c_matrix2 <- c_matrix1 + c(0, (as.matrix(attr(summary(mod_true)$terms,"term.labels")) == c("x:m")))
  # 
  # ate_m0 <- c(coef(mod_true) %*% c_matrix1, 
  #             sqrt(t(c_matrix1) %*% vcov(mod_true) %*% c_matrix1))
  # 
  # ate_m1 <- c(coef(mod_true) %*% c_matrix2, 
  #             sqrt(t(c_matrix2) %*% vcov(mod_true) %*% c_matrix2))
  
  mu1_gcomp <- mean(predict(mod_true, newdata = transform(dat, x = 1), type = "response"))
  mu0_gcomp <- mean(predict(mod_true, newdata = transform(dat, x = 0), type = "response"))
  
  ate_gcomp <- mu1_gcomp - mu0_gcomp
  
  mu11_gcomp <- mean(predict(mod_true, newdata = transform(dat, x = 1, m = 1), type = "response"))
  mu01_gcomp <- mean(predict(mod_true, newdata = transform(dat, x = 0, m = 1), type = "response"))
  
  ate_m1_gcomp <- mu11_gcomp - mu01_gcomp
  
  mu10_gcomp <- mean(predict(mod_true, newdata = transform(dat, x = 1, m = 0), type = "response"))
  mu00_gcomp <- mean(predict(mod_true, newdata = transform(dat, x = 0, m = 0), type = "response"))
  
  ate_m0_gcomp <- mu10_gcomp - mu00_gcomp
  
  boot_func <- function(data, index){
    
    boot_dat <- data[index, ]
    
    mod_true_ <- glm(mod_form, data = boot_dat, family = binomial(link = "logit"))
    
    mu1_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 1), type = "response"))
    mu0_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 0), type = "response"))
    
    ate_gcomp_ <- mu1_ - mu0_
    
    mu11_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 1, m = 1), type = "response"))
    mu01_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 0, m = 1), type = "response"))
    
    ate_m1_gcomp_ <- mu11_ - mu01_
    
    mu10_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 1, m = 0), type = "response"))
    mu00_ <- mean(predict(mod_true_, newdata = transform(boot_dat, x = 0, m = 0), type = "response"))
    
    ate_m0_gcomp_ <- mu10_ - mu00_
    
    return(c(ate_gcomp_,
             ate_m0_gcomp_,
             ate_m1_gcomp_))
    
  }
  
  boot_res <- boot(boot_func, data = dat, R = 500)
  
  ate_gcomp <- c(ate_gcomp, sd(boot_res$t[,1]))
  names(ate_gcomp) <- c("estimate", "std.err")
  
  ate_m0_gcomp <- c(ate_m0_gcomp, sd(boot_res$t[,2]))
  names(ate_m0_gcomp) <- c("estimate", "std.err")
  
  ate_m1_gcomp <- c(ate_m1_gcomp, sd(boot_res$t[,3]))
  names(ate_m1_gcomp) <- c("estimate", "std.err")
  
  ## use super learner to predict propensity score and outcome for both causal forest and DR learner
  ## 
  covariates_matrix <- dat[,c("m", paste0("c", 1:c_number))]
  covariates_matrix_x <- dat[,c("x", "m", paste0("c", 1:c_number))]
  outcome <-  as.matrix(dat$y)
  exposure <- as.matrix(dat$x)
  
  # construct a sample size variable 
  n <- nrow(dat)
  
  # We use 10-fold cross fitting, and determine precisely which of the observations get assigned to each fold:
  num.folds <- 10
  folds <- sort(seq(n) %% num.folds) + 1
  
  mean_learner <- "SL.mean"
  glm_learner <- "SL.glm"
  
  # ranger learner
  ranger_learner <- "SL.ranger"
  
  # glmnet learner
  glmnet_learner <- "SL.glmnet"
  
  # xgboost learner
  xgboost_learner <- "SL.xgboost"
  
  # earth learner
  #earth_learner <- "SL.earth"
  
  sl_lib <- c(ranger_learner,
              glmnet_learner, 
              xgboost_learner, 
              mean_learner,
              glm_learner)
  
  # Specify the number of folds for V-fold cross-validation
  # Use same folds for causal forest and DR learner
  # Doing cross-validation this way automatically deploys cross-fitting
  fold_dat <- tibble(id = 1:n, folds)
  fold_index <- split(fold_dat$id, fold_dat$folds)
  
  fit_mu <- CV.SuperLearner(Y = outcome,
                            X = covariates_matrix_x, 
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index),
                            control = list(saveCVFitLibrary = T),
                            parallel = "seq",
                            verbose = T)
  
  fit_pi <- CV.SuperLearner(Y = exposure,
                            X = covariates_matrix,
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index),#, stratifyCV = TRUE),
                            control = list(saveCVFitLibrary = T),
                            parallel = "seq",
                            verbose = T)
  
  pscore <- as.matrix(fit_pi$SL.predict)
  
  mu_hat <- as.matrix(fit_mu$SL.predict)
  
  ## HERE WE HAVE CAUSAL FOREST CODE TO IDENTIFY APPROPORATE MODIFERS
  # we construct a causal forest using the `causal_forest` function from the grf package:
  # causal forest without using prediction from super learner
  forest <- causal_forest(X = covariates_matrix, 
                             Y = outcome, 
                             W = exposure, 
                             num.trees = 2000,
                             honesty = TRUE,
                             #min.node.size = 10, # might need increase this, tunable
                             #alpha = .05, #  tunable
                             #imbalance.penalty = 0, #  tunable
                             stabilize.splits = TRUE,
                             tune.parameters = "all", # might need change this
                             tune.num.trees = 200,
                             tune.num.reps = 50,
                             tune.num.draws = 1000,
                             compute.oob.predictions = TRUE,
                             num.threads = 10,
                             seed = 123,
                             clusters = folds) # https://bit.ly/42aBuPP
  # causal forest using prediction from super learner
  forest_sl <- causal_forest(X = covariates_matrix, 
                          Y = outcome, 
                          W = exposure, 
                          Y.hat = mu_hat, # use super learner for w hat and Y hat
                          W.hat = pscore,
                          num.trees = 2000,
                          honesty = TRUE,
                          #min.node.size = 10, # might need increase this, tunable
                          #alpha = .05, #  tunable
                          #imbalance.penalty = 0, #  tunable
                          stabilize.splits = TRUE,
                          tune.parameters = "all", # might need change this
                          tune.num.trees = 200,
                          tune.num.reps = 50,
                          tune.num.draws = 1000,
                          compute.oob.predictions = TRUE,
                          num.threads = 10,
                          seed = 123,
                          clusters = folds) # https://bit.ly/42aBuPP
  
  #tau.hat = predict(forest)$predictions
  # start modifying algorithms form here by adding another cf..............
  ## COMPUTE THE ATE FROM THE CAUSAL FOREST ALGORITHM
  # we then use the `average_treatment_effect()` function to obtain an ATE estimate
  cf_ate <- average_treatment_effect(forest)
  cf_sl_ate <- average_treatment_effect(forest_sl)
  
  blp_cf <- tidy(best_linear_projection(forest, covariates_matrix)) %>% 
    filter(term != "(Intercept)") %>% 
    arrange(desc(abs(statistic))) %>% 
    mutate(parameter1 = param1,
           parameter2 = param2,
           sample_size = n,
           confounder_number = c_number,
           seed_number = seed_value,
           method = "BLP_causal_forest",
           group = group)
  
  blp_cf_sl <- tidy(best_linear_projection(forest_sl, covariates_matrix)) %>% 
    filter(term != "(Intercept)") %>% 
    arrange(desc(abs(statistic))) %>% 
    mutate(parameter1 = param1,
           parameter2 = param2,
           sample_size = n,
           confounder_number = c_number,
           seed_number = seed_value,
           method = "BLP_causal_forest",
           group = group)
  ## need to pick the "strong" ones from this, and store the results
  ## what we should do is rank the |test statistic| and pick the first one to check if it's m
  ## maybe keep this whole thing??
  # if(first(blp_cf$term)=="m"){
  #   m_cf = 1
  # } else {
  #   m_cf = 0
  # }
  aipw_scores <- get_scores(forest,
                            subset = NULL,
                            debiasing.weights = NULL,
                            num.trees.for.weights = 2000
  )
  
  aipw_scores_sl <- get_scores(forest_sl,
                            subset = NULL,
                            debiasing.weights = NULL,
                            num.trees.for.weights = 2000
  )
  
  # function to estimate cate based on aipw_scores
  cate_predict <- function(aipw_scores){
    fmla <- as.formula(aipw_scores ~ factor(m)) # Q: why not use tau.hat?
    score_reg <- lm(fmla, 
                    data=transform(dat, aipw_scores=aipw_scores))
    
    c_matrix1 <- c(1, 0)
    c_matrix2 <- c(1, 1)
    
    cate_m0 <- data.frame(t(c(coef(score_reg) %*% c_matrix1, 
                              sqrt(t(c_matrix1) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix1))
    ))
    names(cate_m0) <- c("estimate", "std.err")
    
    cate_m1 <- data.frame(t(c(coef(score_reg) %*% c_matrix2, 
                              sqrt(t(c_matrix2) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix2))
    ))
    names(cate_m1) <- c("estimate", "std.err")
    return(list(cate_m0, cate_m1))
  }
  
  cf_cate <- cate_predict(aipw_scores = aipw_scores)
  cf_cate_m0 <- cf_cate[[1]]
  cf_cate_m1 <- cf_cate[[2]]
  cf_sl_cate <- cate_predict(aipw_scores = aipw_scores_sl)
  cf_sl_cate_m0 <- cf_sl_cate[[1]]
  cf_sl_cate_m1 <- cf_sl_cate[[2]]
  # DO THE EXACT SAME FOR DRLEARNER.
  
  mu_hat1 <- NULL
  for(i in 1:num.folds){
    mu_hat1 <- rbind(mu_hat1, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = base::transform(
                               covariates_matrix_x[fold_index[[i]],], x = 1), 
                             onlySL=T)$pred)
  }
  
  mu_hat0 <- NULL
  for(i in 1:num.folds){
    mu_hat0 <- rbind(mu_hat0, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = base::transform(
                               covariates_matrix_x[fold_index[[i]],], x = 0), 
                             onlySL=T)$pred)
  }
  
  ## aipw
  aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1){
    aipw_score <- ((2*exposure - 1)*(outcome - mu_hat))/((2*exposure - 1)*pscore + (1 - exposure)) + (mu_hat1 - mu_hat0)
    return(aipw_score)
  }
  
  aipw_score <- aipw_func(exposure, 
                          outcome, 
                          pscore, 
                          mu_hat, 
                          mu_hat0, 
                          mu_hat1)
  colnames(aipw_score) <- NULL
  
  aipw_psi <- mean(aipw_score)
  
  aipw_se <- sd(aipw_score)/sqrt(n)
  
  aipw_ate <- data.frame(t(c(aipw_psi, aipw_se)))
  names(aipw_ate) <- c("estimate", "std.err")
  
  ## best linear projection
  score_dat <- data.frame(aipw_score, covariates_matrix)
  score_fit <- lm(aipw_score ~ ., data = score_dat)
  
  blp_aipw <- coeftest(score_fit, vcov = vcovHC(score_fit, type = "HC3"))
  
  blp_drl <- tidy(blp_aipw) %>% 
    filter(term != "(Intercept)") %>% 
    arrange(desc(abs(statistic))) %>% 
    mutate(parameter1 = param1,
           parameter2 = param2,
           sample_size = n,
           confounder_number = c_number,
           seed_number = seed_value,
           method = "BLP_dr_learner",
           group = group) 
  
  
  # if(first(blp_drl$term)=="m"){
  #   m_drl = 1
  # } else {
  #   m_drl = 0
  # }
  
  drl_cate <- cate_predict(aipw_scores = aipw_score)
  drl_cate_m0 <- drl_cate[[1]]
  drl_cate_m1 <- drl_cate[[2]]
  
  # fmla <- as.formula(aipw_score ~ factor(m))
  # score_reg <- lm(fmla, 
  #                 data=transform(dat, aipw_scores=aipw_score))
  # 
  # c_matrix1 <- c(1, 0)
  # c_matrix2 <- c(1, 1)
  # 
  # drl_cate_m0 <- data.frame(t(c(coef(score_reg) %*% c_matrix1, 
  #                               sqrt(t(c_matrix1) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix1))
  # ))
  # names(drl_cate_m0) <- c("estimate", "std.err")
  # 
  # drl_cate_m1 <- data.frame(t(c(coef(score_reg) %*% c_matrix2, 
  #                               sqrt(t(c_matrix2) %*% vcovHC(score_reg, type = "HC3") %*% c_matrix2))
  # ))
  # names(drl_cate_m1) <- c("estimate", "std.err")
  
  # objects to save:
  result <- rbind(
    data.frame(estimator = "oracle_ate",            estimate=rd_ate, std.err=NA),
    data.frame(estimator = "oracle_cate_m0",        estimate=rd_m0, std.err=NA),
    data.frame(estimator = "oracle_cate_m1",        estimate=rd_m1, std.err=NA),
    data.frame(estimator = "gComp_ate",             t(ate_gcomp)),
    data.frame(estimator = "gComp_cate_m0",         t(ate_m0_gcomp)),
    data.frame(estimator = "gComp_cate_m1",         t(ate_m1_gcomp)),
    data.frame(estimator = "causal_forest_ate",     t(cf_ate)),
    data.frame(estimator = "causal_forest_cate_m0", cf_cate_m0),
    data.frame(estimator = "causal_forest_cate_m1", cf_cate_m1),
    data.frame(estimator = "causal_forest_SL_ate",     t(cf_sl_ate)),
    data.frame(estimator = "causal_forest_SL_cate_m0", cf_sl_cate_m0),
    data.frame(estimator = "causal_forest_SL_cate_m1", cf_sl_cate_m1),
    data.frame(estimator = "DRLearner_ate",         aipw_ate),
    data.frame(estimator = "DRLearner_cate_m0",     drl_cate_m0),
    data.frame(estimator = "DRLearner_cate_m1",     drl_cate_m1)
    ) %>%
    mutate(group = group, seed_number = seed_value, sample_size = n, confounder_number = c_number, or_exposure = param1, or_interaction = param2, p_exposure = mean(pi_)) %>%
    select(group, estimator, seed_number, sample_size, confounder_number, or_exposure, or_interaction, p_exposure, estimate, std.err)
  ## writing files for each group of parameters:
  write_csv(blp_cf, here("output",paste0("Blp_cf_c", c_number, "_n", n, "_sims", number_sims, ".csv")), append = T)
  write_csv(blp_cf_sl, here("output",paste0("Blp_cf_sl_c", c_number, "_n", n, "_sims", number_sims, ".csv")), append = T)
  write_csv(blp_drl, here("output",paste0("Blp_drl_c", c_number, "_n", n, "_sims", number_sims, ".csv")), append = T)
  write_csv(result, here("output",paste0("All_results_c", c_number, "_n", n, "_sims", number_sims, ".csv")), append = T)
  return(list(result, blp_cf, blp_cf_sl, blp_drl))
  
}   

