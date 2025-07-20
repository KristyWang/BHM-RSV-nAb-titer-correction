# This function simulates dilution series data. 
# The generated datasets include virus control (VC) data, internal positive control (PC) data,
# International Standard (IS) control data at different concentrations (denoted as IS1 and IS2),
# and serum sample data.

# Load required libraries
library(MASS)
library(dplyr)

simu_func <- function(seed,                                                     # Random seed for reproducibility
                      num_sto,                                                  # Number of working virus stocks
                      num_bat_per_sto,                                          # Number of batches per working virus stocks
                      num_vc_obs_per_bat,                                       # Number of VC observations per batch
                      num_sam_per_bat,                                          # Number of serum samples per batch
                      
                      num_replicate_per_sam = 2,                                # Number of replicates per serum sample
                      dilution_per_sam = c(40, 120, 360, 1080, 3240,            # Serial dilutions
                                           9720, 29160, 87480),
                      
                      sto_effect = c(50, 47, 38, 28),                           # Effects of each working virus stock
                      f_PC = c(0.945, 0.004, 0.977, 6.524),                     # 4PL parameters of PC
                      f_IS1 = c(0.984, -0.019, 0.965, 6.357),                   # 4PL parameters of IS1
                      f_IS2 = c(1.068, -0.018, 0.981, 6.808),                   # 4PL parameters of IS2
                      f_pop = c(0.923, -0.028, 0.972, 5.545),                   # 4PL parameters at the population level
                      
                      Rho_f_bat = matrix(c(1, 0, 0, 0, 0,                       # Correlation matrix of batch-level random effects on 4PL parameters and VC
                                           0, 1, 0, 0, -0.7,
                                           0, 0, 1, 0.6, 0,
                                           0, 0, 0.6, 1, 0,
                                           0, -0.7, 0, 0, 1), 
                                         nrow = 5, ncol = 5, byrow = T),
                      
                      sd_f_bat = c(0.028, 0.211, 0.009, 0.25, 3.67),            # Standard deviations for batch-level effects
                      
                      Rho_f_sam = matrix(c( 1.000,  0.000,  0.000,  0.000,      # Correlation matrix of sample-level random effects
                                            0.000,  1.000,  0.000,  0.000,
                                            0.000,  0.000,  1.000,  0.600,
                                            0.000,  0.000,  0.600,  1.000), 
                                         nrow = 4, ncol = 4, byrow = T),
                      
                      sd_f_sam = c(0.1, 0.1, 0.006, 1.00)){                     # Standard deviations for sample-level effects
  
## Set seed ----
  set.seed(seed)
  
  ## Simulate batch-level effects ----
  bat_effect <- mvrnorm(n = sum(num_bat_per_sto),
                        mu = c(0, 0, 0, 0, 0), 
                        Sigma = diag(sd_f_bat) %*% Rho_f_bat %*% diag(sd_f_bat))
  bat_effect_vc <- bat_effect[, 5]
  bat_effect <- bat_effect[, 1:4]
  
  ## Simulate sample-level effects ----
  n_total <- sum(num_bat_per_sto) * num_sam_per_bat
  mu = c(0, 0, 0, 0)
  Sigma = diag(sd_f_sam) %*% Rho_f_sam %*% diag(sd_f_sam)
  sam_effect <- mvrnorm(n = n_total, mu = mu, Sigma = Sigma)
  
  ## Generate VC data ----
  df_vc <- data.frame(sto_id = rep(1:num_sto, times = num_bat_per_sto),
                      bat_id = 1:sum(num_bat_per_sto),
                      sto_effect = rep(sto_effect, times = num_bat_per_sto),
                      bat_effect = bat_effect_vc)
  df_vc <- df_vc %>% slice(rep(1:n(), each = num_vc_obs_per_bat)) %>%
    mutate(mu_vc_per_bat = sto_effect + bat_effect) 
  df_vc$Count <- rpois(nrow(df_vc), df_vc$mu_vc_per_bat)
  
  ## Generate PC data ----
  df_PC <- data.frame(sto_id = rep(1:num_sto, times = num_bat_per_sto),
                      bat_id = 1:sum(num_bat_per_sto),
                      sto_effect = rep(sto_effect, times = num_bat_per_sto))
  df_PC[, c("b", "c", "d", "e")] <- matrix(rep(f_PC, times = sum(num_bat_per_sto)), ncol = 4, byrow = T) + bat_effect
  df_PC <- df_PC %>% slice(rep(1:n(), each = num_replicate_per_sam)) %>%
    mutate(rpt = rep(1:num_replicate_per_sam, times = sum(num_bat_per_sto))) %>%
    slice(rep(1:n(), each = length(dilution_per_sam))) %>%
    mutate(dilution = rep(dilution_per_sam, times = num_replicate_per_sam * sum(num_bat_per_sto))) %>%
    mutate(mu_y =  c + (d - c) / (1 + exp(b * (log(dilution) - e)))) %>% 
    mutate(mu_Count = sto_effect * (1 - mu_y))
  df_PC$Count <- rpois(nrow(df_PC), df_PC$mu_Count)
  
  ## Generate IS1 data ----
  df_IS1 <- data.frame(sto_id = rep(1:num_sto, times = num_bat_per_sto),
                       bat_id = 1:sum(num_bat_per_sto),
                       sto_effect = rep(sto_effect, times = num_bat_per_sto))
  df_IS1[, c("b", "c", "d", "e")] <- matrix(rep(f_IS1, times = sum(num_bat_per_sto)), ncol = 4, byrow = T) + bat_effect
  df_IS1 <- df_IS1 %>% slice(rep(1:n(), each = num_replicate_per_sam)) %>%
    mutate(rpt = rep(1:num_replicate_per_sam, times = sum(num_bat_per_sto))) %>%
    slice(rep(1:n(), each = length(dilution_per_sam))) %>%
    mutate(dilution = rep(dilution_per_sam, times = num_replicate_per_sam * sum(num_bat_per_sto))) %>%
    mutate(mu_y =  c + (d - c) / (1 + exp(b * (log(dilution) - e)))) %>% 
    mutate(mu_Count = sto_effect * (1 - mu_y))
  df_IS1$Count <- rpois(nrow(df_IS1), df_IS1$mu_Count)
  
  ## Generate IS2 data ----
  df_IS2 <- data.frame(sto_id = rep(1:num_sto, times = num_bat_per_sto),
                       bat_id = 1:sum(num_bat_per_sto),
                       sto_effect = rep(sto_effect, times = num_bat_per_sto))
  df_IS2[, c("b", "c", "d", "e")] <- matrix(rep(f_IS2, times = sum(num_bat_per_sto)), ncol = 4, byrow = T) + bat_effect
  df_IS2 <- df_IS2 %>% slice(rep(1:n(), each = num_replicate_per_sam)) %>%
    mutate(rpt = rep(1:num_replicate_per_sam, times = sum(num_bat_per_sto))) %>%
    slice(rep(1:n(), each = length(dilution_per_sam))) %>%
    mutate(dilution = rep(dilution_per_sam, times = num_replicate_per_sam * sum(num_bat_per_sto))) %>%
    mutate(mu_y =  c + (d - c) / (1 + exp(b * (log(dilution) - e)))) %>% 
    mutate(mu_Count = sto_effect * (1 - mu_y))
  df_IS2$Count <- rpois(nrow(df_IS2), df_IS2$mu_Count)
  
  ## Generate serum sample data ----
  df_sam <- data.frame(sto_id = rep(1:num_sto, times = num_bat_per_sto),
                       bat_id = 1:sum(num_bat_per_sto),
                       sto_effect = rep(sto_effect, times = num_bat_per_sto),
                       b_bat_effect = bat_effect[, 1],
                       c_bat_effect = bat_effect[, 2],
                       d_bat_effect = bat_effect[, 3],
                       e_bat_effect = bat_effect[, 4])
  
  df_sam <- df_sam %>% slice(rep(1:n(), each = num_sam_per_bat)) %>% 
    mutate(simu_sam_id = 1:(sum(num_bat_per_sto) * num_sam_per_bat),
           sam_bat_id = rep(1:num_sam_per_bat, times = sum(num_bat_per_sto))) %>% 
    mutate(b = f_pop[1] + b_bat_effect + sam_effect[, 1],
           c = f_pop[2] + c_bat_effect + sam_effect[, 2],
           d = f_pop[3] + d_bat_effect + sam_effect[, 3],
           e = f_pop[4] + e_bat_effect + sam_effect[, 4]) %>%
    slice(rep(1:n(), each = num_replicate_per_sam)) %>%
    mutate(rpt = rep(1:num_replicate_per_sam, times = (sum(num_bat_per_sto) * num_sam_per_bat))) %>%
    slice(rep(1:n(), each = length(dilution_per_sam))) %>%
    mutate(dilution = rep(dilution_per_sam, times = (num_replicate_per_sam * sum(num_bat_per_sto) * num_sam_per_bat))) %>%
    mutate(mu_y = c + (d - c) / (1 + exp(b * (log(dilution) - e)))) %>% 
    mutate(mu_Count = sto_effect * (1 - mu_y))
  
  df_sam$mu_Count <- ifelse(df_sam$mu_Count < 0, 0, df_sam$mu_Count)
  df_sam$Count <- rpois(nrow(df_sam), df_sam$mu_Count)
  
  ## Calculate true titers ----
  which_bat_per_sam = df_sam %>% group_by(simu_sam_id) %>% slice(1) %>%
    dplyr::select(sto_id, bat_id, simu_sam_id, sam_bat_id)
  
  which_bat_per_sam$Nt_true = NA
  for (i in 1:nrow(which_bat_per_sam)){
    bat_id_i <- which_bat_per_sam[i, ]$bat_id
    bat_effect_i <- bat_effect[bat_id_i, ]
    f_i <- f_pop + bat_effect_i + sam_effect[i, ]
    f_adj_i <-f_pop + sam_effect[i, ]
    
    if (f_adj_i[2] > 0.5){
      which_bat_per_sam[i, "Nt_true"] = log(87480);
    } else if (f_adj_i[3] < 0.5){
      which_bat_per_sam[i, "Nt_true"] = log(40) - log(3);
    } else {
      which_bat_per_sam[i, "Nt_true"] = log((f_adj_i[3]-f_adj_i[2])/(0.5-f_adj_i[2])-1)/(f_adj_i[1]) + f_adj_i[4];
    }
  }
  
  which_bat_per_sam$Nt_true = log2(exp(which_bat_per_sam$Nt_true))
  Nt_PC_true <- log2(exp(log((f_PC[3]-f_PC[2])/(0.5-f_PC[2])-1)/(f_PC[1]) + f_PC[4]))
  Nt_IS1_true <- log2(exp(log((f_IS1[3]-f_IS1[2])/(0.5-f_IS1[2])-1)/(f_IS1[1]) + f_IS1[4]))
  Nt_IS2_true <- log2(exp(log((f_IS2[3]-f_IS2[2])/(0.5-f_IS2[2])-1)/(f_IS2[1]) + f_IS2[4]))
  
  ## Return all simulated data ----
  return(list(df_vc = df_vc,
              df_PC = df_PC,
              df_IS1 = df_IS1,
              df_IS2 = df_IS2,
              df_sam = df_sam,
              Nt_PC_true = Nt_PC_true,
              Nt_IS1_true = Nt_IS1_true, 
              Nt_IS2_true = Nt_IS2_true,
              sto_effect = sto_effect,
              bat_effect_vc = bat_effect_vc,
              f_PC = f_PC,
              f_IS1 = f_IS1,
              f_IS2 = f_IS2,
              f_pop = f_pop,
              bat_effect = bat_effect,
              sam_effect = sam_effect,
              which_bat_per_sam = which_bat_per_sam))
  
}

# An example ----
# simu_data <- simu_func(seed = 123,
#                        num_sto = 4,
#                        num_bat_per_sto = c(7, 7, 7, 7),
#                        num_vc_obs_per_bat = 22,
#                        num_sam_per_bat = 20)
