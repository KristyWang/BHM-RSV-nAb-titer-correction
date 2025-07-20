# Load required libraries
library(rstan)
library(bayesplot)

# ---------------------------------------------------
# Model 1: Effect of Each Working Virus Stock
# ---------------------------------------------------
# This model estimates the effect (delta) of each working virus stock.
# Only virus control (VC) data are used for model fitting.
# The goal is to isolate how different working virus stocks contribute to 
# observed foci counts in the VC condition (i.e., background signal),
# while accounting for batch-to-batch variability through a random effect.

model1_code <- "

data {
  
  int<lower=1> num_sto_vc;                                            // Number of working virus stocks
  int<lower=1> num_bat_vc;                                            // Number of batches
  int<lower=1> num_obs_vc;                                            // Number of observations 
  int<lower=1, upper=num_sto_vc> which_sto_vc[num_obs_vc];            // Index mapping each observation to a working virus stock
  int<lower=1, upper=num_bat_vc> which_bat_vc[num_obs_vc];            // Index mapping each observation to a batch
  int<lower=1, upper=num_bat_vc> which_sto_per_bat[num_bat_vc];       // Index mapping each batch to a working virus stock
  int<lower=0> count_vc[num_obs_vc];                                  // Foci count for each observation
  
  vector<lower=0>[num_sto_vc] delta_prior;                            // Prior means for the foci count of each working virus stock
  real<lower=0, upper=1> sample_from_prior;                           // Flag to control whether to sample from priors only (1=Yes, 0=No)
  
}

parameters {

  vector[num_sto_vc] z_delta;                                          // Standardized working virus stock effects (to be scaled later)
  vector[num_bat_vc] z_bat_vc;                                         // Standardized batch effects
  real<lower=0> sd_bat_vc;                                             // Standard deviation for batch effects
 
}

transformed parameters{

  vector[num_sto_vc] delta = delta_prior + 10 * z_delta;                // Working virus stock effects: prior + scaled z_delta
  vector[num_bat_vc] bat_vc = z_bat_vc * sd_bat_vc;                     // Batch effects: scaled by sd_bat_vc
  vector[num_obs_vc] mu_vc;                                             // Linear predictor for each observation
  
  for (i in 1:num_obs_vc){                                              // Construct the linear predictor 
      mu_vc[i] = delta[which_sto_vc[i]] + bat_vc[which_bat_vc[i]];
  }
  
}

model{

  // priors
  z_delta ~ std_normal();
  z_bat_vc ~ std_normal();
  sd_bat_vc ~ std_normal();
  
  // Likelihood
  if (sample_from_prior == 0){
    target += poisson_lpmf(count_vc | mu_vc);
  }
 
}

generated quantities{

  // Posterior mean foci counts for each batch 
  vector[num_bat_vc] mu_vc_per_bat;
  for (i in 1:num_bat_vc){
    mu_vc_per_bat[i] = delta[which_sto_per_bat[i]] + bat_vc[i];
  }
  
}
"

# ---------------------------------------------------
# Model 2: Estimating nAb Titers with Batch Effect Correction
# ---------------------------------------------------
# This model estimates neutralizing antibody (nAb) titers for each serum sample.
# The relationship between dilutions and foci reduction is modeled using
# a Bayesian hierarchical 4-parameter logistic (4PL) model,
# with random effects included for both serum samples and assay batches.
#
# Foci reduction is calculated relative to the modeled mean foci counts 
# from virus control (VC) data for each working virus stock, as estimated in Model 1.
#
# Dilution series data from assay controls—including the internal positive control (PC)
# and the International Standard control (IS)—are used to estimate and correct for batch effects.
# Dilution series from serum samples are modeled to generate accurate nAb titer estimates,
# appropriately adjusted for batch-to-batch variation.

model2_code <- "

 // 4-parameter logistic (4PL) function
functions {
  real LL4(real x_log, real b, real c, real d, real e) {              // x_log: log-transformed dilution
    return c + (d - c) / (1 + exp(b * (x_log - e)));                  // c: lower asymptote; d: upper asymptote; 
  }                                                                   // b: slope; e: log-transformed midrange dilution
}

data {

  int<lower=1> num_bat_total;                                         // Total number of batches
  
  // Assay control: internal positive control (PC)
  int<lower=1> num_obs_pc;                                            // number of observations
  int<lower=1, upper=num_bat_total> which_bat_pc[num_obs_pc];         // batch index
  vector<lower=log(40), upper=log(87480)>[num_obs_pc] x_log_pc;       // log dilution
  vector<upper=1>[num_obs_pc] y_pc;                                   // Foci reduction
  
  // Assay control: International Standard Control (IS500)
  int<lower=1> num_obs_IS1;                                           // number of observations
  int<lower=1, upper=num_bat_total> which_bat_IS1[num_obs_IS1];       // batch index
  vector<lower=log(40), upper=log(87480)>[num_obs_IS1] x_log_IS1;     // log dilution
  vector<upper=1>[num_obs_IS1] y_IS1;                                 // Foci reduction
  
   // Assay control: International Standard Control (IS1000)
  int<lower=1> num_obs_IS2;                                           // number of observations
  int<lower=1, upper=num_bat_total> which_bat_IS2[num_obs_IS2];       // batch index
  vector<lower=log(40), upper=log(87480)>[num_obs_IS2] x_log_IS2;     // log dilution
  vector<upper=1>[num_obs_IS2] y_IS2;                                 // Foci reduction
  
  // Serum sample data
  int<lower=1> num_sam;                                               // number of serum samples
  int<lower=1> num_obs;                                               // number of observations
  int<lower=1, upper=num_bat_total> which_bat[num_obs];               // batch index
  int<lower=1, upper=num_sam> which_sam[num_obs];                     // sample index
  vector<lower=log(40), upper=log(87480)>[num_obs] x_log;             // log dilution
  vector<upper=1>[num_obs] y;                                         // Foci reduction
  
  // Prior means
  vector[4] f_pc_prior;                                               // 4PL parameters of PC
  vector[4] f_IS1_prior;                                              // 4PL parameters of IS1
  vector[4] f_IS2_prior;                                              // 4PL parameters of IS2
  vector[4] f_pop_prior;                                              // 4PL parameters at the population level
  
  // Lower and upper bounds
  vector[4] f_lower;                                                  // Lower bounds for 4PL parameters (population-level and assay controls)
  vector[4] f_upper;                                                  // Upper bounds for 4PL parameters (population-level and assay controls)
  
  int<lower=1, upper=num_bat_total> which_bat_per_sam[num_sam];       // Index mapping each sample to a batch
  real<lower=0, upper=1> sample_from_prior;                           // Flag to control whether to sample from priors only (1=Yes, 0=No)
  
}

parameters {

  // Standardized 4PL parameters of PC (to be scaled later)
  real<lower=(f_lower[1]-f_pc_prior[1])/0.25, upper=(f_upper[1]-f_pc_prior[1])/0.25> b_pc_raw;
  real<lower=(f_lower[2]-f_pc_prior[2])/0.05, upper=(f_upper[2]-f_pc_prior[2])/0.05> c_pc_raw;
  real<lower=(f_lower[3]-f_pc_prior[3])/0.05, upper=(f_upper[3]-f_pc_prior[3])/0.05> d_pc_raw;
  real<lower=(f_lower[4]-f_pc_prior[4])/1.00, upper=(f_upper[4]-f_pc_prior[4])/1.00> e_pc_raw;
  
  // Standardized 4PL parameters of IS500
  real<lower=(f_lower[1]-f_IS1_prior[1])/0.25, upper=(f_upper[1]-f_IS1_prior[1])/0.25> b_IS1_raw;
  real<lower=(f_lower[2]-f_IS1_prior[2])/0.05, upper=(f_upper[2]-f_IS1_prior[2])/0.05> c_IS1_raw;
  real<lower=(f_lower[3]-f_IS1_prior[3])/0.05, upper=(f_upper[3]-f_IS1_prior[3])/0.05> d_IS1_raw;
  real<lower=(f_lower[4]-f_IS1_prior[4])/1.00, upper=(f_upper[4]-f_IS1_prior[4])/1.00> e_IS1_raw;
  
  // Standardized 4PL parameters of IS1000
  real<lower=(f_lower[1]-f_IS2_prior[1])/0.25, upper=(f_upper[1]-f_IS2_prior[1])/0.25> b_IS2_raw;
  real<lower=(f_lower[2]-f_IS2_prior[2])/0.05, upper=(f_upper[2]-f_IS2_prior[2])/0.05> c_IS2_raw;
  real<lower=(f_lower[3]-f_IS2_prior[3])/0.05, upper=(f_upper[3]-f_IS2_prior[3])/0.05> d_IS2_raw;
  real<lower=(f_lower[4]-f_IS2_prior[4])/1.00, upper=(f_upper[4]-f_IS2_prior[4])/1.00> e_IS2_raw;
  
  // Standardized 4PL parameters for serum samples at the population level
  real<lower=(f_lower[1]-f_pop_prior[1])/0.25, upper=(f_upper[1]-f_pop_prior[1])/0.25> b_pop_raw;
  real<lower=(f_lower[2]-f_pop_prior[2])/0.05, upper=(f_upper[2]-f_pop_prior[2])/0.05> c_pop_raw;
  real<lower=(f_lower[3]-f_pop_prior[3])/0.05, upper=(f_upper[3]-f_pop_prior[3])/0.05> d_pop_raw;
  real<lower=(f_lower[4]-f_pop_prior[4])/1.00, upper=(f_upper[4]-f_pop_prior[4])/1.00> e_pop_raw;
  
  // Standardized batch-level random effects
  matrix[num_bat_total, 4] z_f_bat;
  vector<lower=0>[4] z_sd_f_bat;
  
  // Standardized sample-level random effects (including batch effects)
  matrix[num_sam, 4] z_f_sam_bat;
  vector<lower=0>[4] z_sd_f_sam_bat;
  
  // Standardized observation-level noise
  real<lower=0> z_sd_y;
  
}

transformed parameters{

  // Rescaled 4PL parameters of PC
  vector[4] f_pc;
  f_pc[1] = f_pc_prior[1] + 0.25 * b_pc_raw;
  f_pc[2] = f_pc_prior[2] + 0.05 * c_pc_raw;
  f_pc[3] = f_pc_prior[3] + 0.05 * d_pc_raw;
  f_pc[4] = f_pc_prior[4] + e_pc_raw;
  
  // Rescaled 4PL parameters of IS500
  vector[4] f_IS1;
  f_IS1[1] = f_IS1_prior[1] + 0.25 * b_IS1_raw;
  f_IS1[2] = f_IS1_prior[2] + 0.05 * c_IS1_raw;
  f_IS1[3] = f_IS1_prior[3] + 0.05 * d_IS1_raw;
  f_IS1[4] = f_IS1_prior[4] + e_IS1_raw;
  
  // Rescaled 4PL parameters of IS1000
  vector[4] f_IS2;
  f_IS2[1] = f_IS2_prior[1] + 0.25 * b_IS2_raw;
  f_IS2[2] = f_IS2_prior[2] + 0.05 * c_IS2_raw;
  f_IS2[3] = f_IS2_prior[3] + 0.05 * d_IS2_raw;
  f_IS2[4] = f_IS2_prior[4] + e_IS2_raw;
  
  // Rescaled population-level 4PL parameters of serum samples
  vector[4] f_pop;
  f_pop[1] = f_pop_prior[1] + 0.25 * b_pop_raw;
  f_pop[2] = f_pop_prior[2] + 0.05 * c_pop_raw;
  f_pop[3] = f_pop_prior[3] + 0.05 * d_pop_raw;
  f_pop[4] = f_pop_prior[4] + e_pop_raw;
  
  // Rescaled batch-level random effects 
  vector[4] sd_f_bat;
  sd_f_bat[1] = 0.5 * z_sd_f_bat[1];
  sd_f_bat[2] = 0.5 * z_sd_f_bat[2];
  sd_f_bat[3] = 0.5 * z_sd_f_bat[3];
  sd_f_bat[4] = 1.0 * z_sd_f_bat[4];
  matrix[num_bat_total, 4] f_bat = z_f_bat * diag_matrix(sd_f_bat);
  
  // Rescaled sample-level random effects (including batch effects)
  vector[4] sd_f_sam_bat;
  sd_f_sam_bat[1] = 0.5 * z_sd_f_sam_bat[1];
  sd_f_sam_bat[2] = 1.0 * z_sd_f_sam_bat[2];
  sd_f_sam_bat[3] = 0.5 * z_sd_f_sam_bat[3];
  sd_f_sam_bat[4] = 2.0 * z_sd_f_sam_bat[4];
  matrix[num_sam, 4] f_sam_bat = z_f_sam_bat * diag_matrix(sd_f_sam_bat);
  
  // Predicted foci reduction for PC
  vector[num_obs_pc] mu_y_pc;
  for (i in 1:num_obs_pc){
    real b = f_pc[1] + f_bat[which_bat_pc[i], 1];
    real c = f_pc[2] + f_bat[which_bat_pc[i], 2];
    real d = f_pc[3] + f_bat[which_bat_pc[i], 3];
    real e = f_pc[4] + f_bat[which_bat_pc[i], 4];
    mu_y_pc[i] = LL4(x_log_pc[i], b, c, d, e);
  }
  
  // Predicted foci reduction for IS500
  vector[num_obs_IS1] mu_y_IS1;
  for (i in 1:num_obs_IS1){
    real b = f_IS1[1] + f_bat[which_bat_IS1[i], 1];
    real c = f_IS1[2] + f_bat[which_bat_IS1[i], 2];
    real d = f_IS1[3] + f_bat[which_bat_IS1[i], 3];
    real e = f_IS1[4] + f_bat[which_bat_IS1[i], 4];
    mu_y_IS1[i] = LL4(x_log_IS1[i], b, c, d, e);
  }
  
  // Predicted foci reduction for IS1000
  vector[num_obs_IS2] mu_y_IS2;
  for (i in 1:num_obs_IS2){
    real b = f_IS2[1] + f_bat[which_bat_IS2[i], 1];
    real c = f_IS2[2] + f_bat[which_bat_IS2[i], 2];
    real d = f_IS2[3] + f_bat[which_bat_IS2[i], 3];
    real e = f_IS2[4] + f_bat[which_bat_IS2[i], 4];
    mu_y_IS2[i] = LL4(x_log_IS2[i], b, c, d, e);
  }
  
  // Predicted foci reduction for serum samples
  vector[num_obs] mu_y;
  for (i in 1:num_obs){
    real b = f_pop[1] + f_sam_bat[which_sam[i], 1];
    real c = f_pop[2] + f_sam_bat[which_sam[i], 2];
    real d = f_pop[3] + f_sam_bat[which_sam[i], 3];
    real e = f_pop[4] + f_sam_bat[which_sam[i], 4];
    mu_y[i] = LL4(x_log[i], b, c, d, e);
  }
  
  // Observation-level noise
  real sd_y = 0.15 * z_sd_y;
  
}

model{

  // Priors
  b_pc_raw ~ normal(0, 1); 
  c_pc_raw ~ normal(0, 1); 
  d_pc_raw ~ normal(0, 1); 
  e_pc_raw ~ normal(0, 1); 
  
  b_IS1_raw ~ normal(0, 1); 
  c_IS1_raw ~ normal(0, 1); 
  d_IS1_raw ~ normal(0, 1); 
  e_IS1_raw ~ normal(0, 1); 
  
  b_IS2_raw ~ normal(0, 1); 
  c_IS2_raw ~ normal(0, 1); 
  d_IS2_raw ~ normal(0, 1); 
  e_IS2_raw ~ normal(0, 1); 
  
  b_pop_raw ~ normal(0, 1); 
  c_pop_raw ~ normal(0, 1); 
  d_pop_raw ~ normal(0, 1); 
  e_pop_raw ~ normal(0, 1); 
  
  to_vector(z_f_bat) ~ std_normal();
  z_sd_f_bat ~ normal(0, 1); 
  
  to_vector(z_f_sam_bat) ~ std_normal();
  z_sd_f_sam_bat ~ normal(0, 1); 
  
  z_sd_y ~ normal(0, 1); 

  // Likelihood
  if (sample_from_prior == 0){
    target += normal_lpdf(y_pc | mu_y_pc, sd_y);
    target += normal_lpdf(y_IS1 | mu_y_IS1, sd_y);
    target += normal_lpdf(y_IS2 | mu_y_IS2, sd_y);
    target += normal_lpdf(y | mu_y, sd_y);
  }

}

generated quantities{

  // Calculate nAb titers for assay controls
  real Nt_PC;                                                         // nAb titer for PC
  real Nt_IS1;                                                        // nAb titer for IS500
  real Nt_IS2;                                                        // nAb titer for IS1000
  
  if (f_pc[2] > 0.5){
    Nt_PC = log(87480);
  } else if (f_pc[3] < 0.5){
    Nt_PC = log(40) - log(3);
  } else{
    Nt_PC = log((f_pc[3]-f_pc[2])/(0.5-f_pc[2])-1)/(f_pc[1]) + f_pc[4];
  }
  
  if (f_IS1[2] > 0.5){
    Nt_IS1 = log(87480);
  } else if (f_IS1[3] < 0.5){
    Nt_IS1 = log(40) - log(3);
  } else{
    Nt_IS1 = log((f_IS1[3]-f_IS1[2])/(0.5-f_IS1[2])-1)/(f_IS1[1]) + f_IS1[4];
  }
  
  if (f_IS2[2] > 0.5){
    Nt_IS2 = log(87480);
  } else if (f_IS2[3] < 0.5){
    Nt_IS2 = log(40) - log(3);
  } else{
    Nt_IS2 = log((f_IS2[3]-f_IS2[2])/(0.5-f_IS2[2])-1)/(f_IS2[1]) + f_IS2[4];
  }
  
  // Fitted foci reduction data for assay controls
  vector[8] y_pc_adj;
  vector[8] y_IS1_adj;
  vector[8] y_IS2_adj;
  
  for (i in 1:8){
    y_pc_adj[i] = LL4(log(40 * 3 ^ (i-1)), f_pc[1], f_pc[2], f_pc[3], f_pc[4]);
  }
  
  for (i in 1:8){
    y_IS1_adj[i] = LL4(log(40 * 3 ^ (i-1)), f_IS1[1], f_IS1[2], f_IS1[3], f_IS1[4]);
  }
  
  for (i in 1:8){
    y_IS2_adj[i] = LL4(log(40 * 3 ^ (i-1)), f_IS2[1], f_IS2[2], f_IS2[3], f_IS2[4]);
  }
  
  // Calculate nAb titers for serum samples
  vector[num_sam] Nt_sam;                                               // BHM-unadjusted nAb titers
  vector[num_sam] Nt_sam_adj;                                           // BHM-adjusted nAb titers
  matrix[8, num_sam] y_sam_adj;                                         // BHM-adjusted foci reduction data 
  
  for (i in 1:num_sam){
    real b = f_pop[1] + f_sam_bat[i, 1];
    real c = f_pop[2] + f_sam_bat[i, 2];
    real d = f_pop[3] + f_sam_bat[i, 3];
    real e = f_pop[4] + f_sam_bat[i, 4];
    
    real b_adj = b - f_bat[which_bat_per_sam[i], 1];
    real c_adj = c - f_bat[which_bat_per_sam[i], 2];
    real d_adj = d - f_bat[which_bat_per_sam[i], 3];
    real e_adj = e - f_bat[which_bat_per_sam[i], 4];
    
    for (l in 1:8){
      y_sam_adj[l, i] = LL4(log(40 * 3 ^ (l-1)), b_adj, c_adj, d_adj, e_adj);
    }
    
    if (b < 0){
      Nt_sam[i] = log(40) - log(3);
    } else if (c > 0.5){
      Nt_sam[i] = log(87480);
    } else if (d < 0.5){
      Nt_sam[i] = log(40) - log(3);
    } else {
      Nt_sam[i] = log((d-c)/(0.5-c)-1)/(b) + e;
    }
    
    if (b_adj < 0){
      Nt_sam_adj[i] = log(40) - log(3);
    } else if (c_adj > 0.5){
      Nt_sam_adj[i] = log(87480);
    } else if (d_adj < 0.5){
      Nt_sam_adj[i] = log(40) - log(3);
    } else {
      Nt_sam_adj[i] = log((d_adj-c_adj)/(0.5-c_adj)-1)/(b_adj) + e_adj;
    }
  }
  
}

"