pacman::p_load(# install if needed and load packages 
  here,         
  tidyverse,     
  lmtest,
  sandwich,
  broom,
  mvtnorm  
)

# 	different effect parametrizations
# •	risk differences
# •	prevalence of Y = .1
# •	RD = -0.025, 0, 0.025 cross-combined
# 	sample sizes
# •	500, 1000
# 	and number of confounding variables
# •	5, 10, 15


# CREATE EXPIT AND LOGIT FUNCTIONS
expit <- function(x){ exp(x)/(1+exp(x)) }
logit <- function(x){ log(x/(1-x)) }

data_sim <- function(sample_size, c_number, param1, param2, exposure, referent, modifier1, modifier0){
  
    param1 = exp(param1)
    param2 = exp(param2)
    n = sample_size
    p = c_number
    set.seed(123)
    
    ## CONFOUNDERS
    sigma <- matrix(0,nrow=p,ncol=p); diag(sigma) <- 1
    c     <- rmvnorm(n, mean=rep(0,p), sigma=sigma)
    
    # DESIGN MATRIX FOR c
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

    balancing_intercept <- -log((1/.1)-1) - # outcome mu of .1
      log(param1)*pi_ -       # offsetting exposure effect
      log(2)*pi_ +            # offsetting modifier effect: outcome-modifier OR: 2
      log(param2)*pi_*pi_   # offsetting interaction term

    # mu_exposure_modifier: Q: what if modifier is not a confounder?
    mu11 <- expit(balancing_intercept + 
                  log(param1)*exposure + log(2)*modifier1 - 
                  log(param2)*exposure*modifier1 + dMat[,-1]%*%beta)
    
    mu10 <- expit(balancing_intercept + 
                    log(param1)*exposure + log(2)*modifier0 - 
                    log(param2)*exposure*modifier0 + dMat[,-1]%*%beta)
    
    mu01 <- expit(balancing_intercept + 
                    log(param1)*referent + log(2)*modifier1 - 
                    log(param2)*referent*modifier1 + dMat[,-1]%*%beta)
    
    mu00 <- expit(balancing_intercept + 
                    log(param1)*referent + log(2)*modifier0 - 
                    log(param2)*referent*modifier0 + dMat[,-1]%*%beta)
    
    rd_m0 <- mean(mu10 - mu00)
    rd_m1 <- mean(mu11 - mu01)
    #rd <- mean(rd_m0*(1-pi_) + rd_m1*pi_)
    # update:
    mu1 <-  expit(balancing_intercept +
                    log(param1)*1 + log(2)*pi_ -
                    log(param2)*1*pi_ + dMat[,-1]%*%beta)
    mu0 <-  expit(balancing_intercept +
                    log(param1)*0 + log(2)*pi_ -
                    log(param2)*0*pi_ + dMat[,-1]%*%beta)
    rd <- mean(mu1-mu0)
    
    
    res <- data.frame(
      sample_size = sample_size,
      mean_pi = mean(pi_),
      c_dim = c_number,
      param1 = param1, 
      param2 = param2,
      risk_diff_m0 = rd_m0,
      risk_diff_m1 = rd_m1,
      risk_diff = rd
      )
    
    return(res)
  
}

# rd_m0 = 0.025
#param1_list <- c(1.31) 
#param2_list <- c(1.065, 1.31, 1.625)
## c_number = 15
param1_list <- c(1.29) 
param2_list <- c(1.05, 1.29, 1.6)
## c_number = 10
param1_list <- c(1.291) 
param2_list <- c(1.06, 1.291, 1.59)
## c_number = 5
param1_list <- c(1.31) 
param2_list <- c(1.08, 1.31, 1.61)

# rd_m0 = 0
## c_number = 15
param1_list <- c(1) 
param2_list <- c(.815, 1, 1.24)
## c_number = 10
param1_list <- c(1) 
param2_list <- c(.825, 1, 1.23)
## c_number = 5
param1_list <- c(1) 
param2_list <- c(.83, 1, 1.22)

# rd_m0 = -0.025
## c_number = 15
param1_list <- c(.765) 
param2_list <- c(.625, .765, .945)
## c_number = 10
param1_list <- c(.761) 
param2_list <- c(.63, .761, .93)
## c_number = 5
param1_list <- c(.745) 
param2_list <- c(.62, .745, .90)

true_dat <- NULL
for(i in param1_list){
  for(j in param2_list){
    true_dat <- rbind(true_dat, 
                      data_sim(sample_size = 5e5,
                                  c_number = 10, 
                                  param1 = log(i), 
                                  param2 = log(j), 
                                  exposure = 1, 
                                  referent = 0, 
                                  modifier1 = 1, 
                                  modifier0 = 0)
    ) 
  }
}

true_dat

# save the truth to a file
## c_number = 15
param_dat1 <- list(param1 = 1.29,
                   param2 = c(1.05, 1.29, 1.6))

data_grid1 <- expand.grid(param_dat1, KEEP.OUT.ATTRS = FALSE)

param_dat2 <- list(param1 = 1,
                   param2 = c(.815, 1, 1.24))

data_grid2 <- expand.grid(param_dat2, KEEP.OUT.ATTRS = FALSE)

param_dat3 <- list(param1 = 0.765,
                   param2 = c(.625, .765, .945))

data_grid3 <- expand.grid(param_dat3, KEEP.OUT.ATTRS = FALSE)

data_grid <- rbind(data_grid1, data_grid2, data_grid3)
table <- lapply(1:nrow(data_grid), 
                function(x) data_sim(sample_size = 5e6,
                                    c_number = 15, 
                                    param1 = log(data_grid[x,]$param1), 
                                    param2 = log(data_grid[x,]$param2), 
                                    exposure = 1, 
                                    referent = 0, 
                                    modifier1 = 1, 
                                    modifier0 = 0))


truth <- do.call(rbind, table)
#write_csv(truth, "H:/RA/HTE/Data/truth_data_simulation.csv")

## c_number= 10
param_dat1 <- list(param1 = 1.291,
                   param2 = c(1.06, 1.291, 1.59))

data_grid1 <- expand.grid(param_dat1, KEEP.OUT.ATTRS = FALSE)

param_dat2 <- list(param1 = 1,
                   param2 = c(.825, 1, 1.23))

data_grid2 <- expand.grid(param_dat2, KEEP.OUT.ATTRS = FALSE)

param_dat3 <- list(param1 = 0.761,
                   param2 = c(.63, .761, .93))

data_grid3 <- expand.grid(param_dat3, KEEP.OUT.ATTRS = FALSE)

data_grid <- rbind(data_grid1, data_grid2, data_grid3)
table2 <- lapply(1:nrow(data_grid), 
                function(x) data_sim(sample_size = 5e6,
                                     c_number = 10, 
                                     param1 = log(data_grid[x,]$param1), 
                                     param2 = log(data_grid[x,]$param2), 
                                     exposure = 1, 
                                     referent = 0, 
                                     modifier1 = 1, 
                                     modifier0 = 0))


truth2 <- do.call(rbind, table2)


## c_number= 5
param_dat1 <- list(param1 = 1.31,
                   param2 = c(1.08, 1.31, 1.61))

data_grid1 <- expand.grid(param_dat1, KEEP.OUT.ATTRS = FALSE)

param_dat2 <- list(param1 = 1,
                   param2 = c(.83, 1, 1.22))

data_grid2 <- expand.grid(param_dat2, KEEP.OUT.ATTRS = FALSE)

param_dat3 <- list(param1 = .745,
                   param2 = c(.62, .745, .90))

data_grid3 <- expand.grid(param_dat3, KEEP.OUT.ATTRS = FALSE)

data_grid <- rbind(data_grid1, data_grid2, data_grid3)
table3 <- lapply(1:nrow(data_grid), 
                 function(x) data_sim(sample_size = 5e6,
                                      c_number = 5, 
                                      param1 = log(data_grid[x,]$param1), 
                                      param2 = log(data_grid[x,]$param2), 
                                      exposure = 1, 
                                      referent = 0, 
                                      modifier1 = 1, 
                                      modifier0 = 0))


truth3 <- do.call(rbind, table3)

truth_all <- rbind(truth, truth2, truth3)
write_csv(truth_all, "H:/RA/HTE/Data/truth_data_simulation_all_update.csv")

# Sigma = matrix(c(1,0.5,0.5,1), ncol=2)
# R = chol(Sigma) # Sigma == t(R)%*%  R
# n = 5
# X = t(R) %*% matrix(rnorm(n*2), 2)
# 
# X %*% t(X)/n # test