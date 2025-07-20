# This script is designed for the simulation study presented in the manuscript.
# It consists of the following four parts:
# 1) Generating simulated data that capture key features of the real RSV FRNT experimental data,
# 2) Applying standard methods — including the Karber method and the 4PL model — to the simulated data,
# 3) Fitting Bayesian hierarchical model (BHM) to the simulated data,
# 4) Creating figures to visualize and interpret the simulation results.

# Clear workspace
rm(list = ls())

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Source custom functions and model scripts
source("1_Func_Karber_and_4PL.R")
source("2_Model_Stan.R")
source("3_Simu_Func.R")

# 1. Generate simulated data ----
simu_data <- simu_func(seed = 20250428,
                       num_sto = 4,
                       num_bat_per_sto = c(7, 7, 7, 7),
                       num_vc_obs_per_bat = 22,
                       num_sam_per_bat = 20)
rm(simu_func)
save(simu_data, file = "../1. Data/simu_data.RData")

# Check simulated data
selected_samples <- sample(simu_data$which_bat_per_sam$simu_sam_id, 40)
simu_data$df_sam %>% filter(simu_sam_id %in% selected_samples) %>%
  ggplot(aes(x = log(dilution), y = mu_y)) +
  geom_point() +
  facet_wrap(~ simu_sam_id, nrow = 5) 
rm(selected_samples)

# 2. Apply standard methods to the simulated data ----
mean_vc <- simu_data$df_vc %>% 
  group_by(bat_id) %>% 
  summarise(mean_vc = mean(Count))

df_sam <- left_join(simu_data$df_sam, mean_vc)

# Karber method
Karber_data <- Karber_function(data = df_sam, 
                               sam_id_var_name = "simu_sam_id")
# 4PL Model
FRNT50 <- DRM_function(data = df_sam, 
                       sam_id_var_name = "simu_sam_id")
# Merge Karber and 4PL results
Karber_and_4PL <- left_join(Karber_data, FRNT50[, c(1,2,7)])
Karber_and_4PL <- left_join(simu_data$which_bat_per_sam, Karber_and_4PL)

Karber_and_4PL[Karber_and_4PL$Nt_true <= log2(40/3), ]$Nt_true <-  log2(40/3)
Karber_and_4PL[Karber_and_4PL$Nt.Karber <= log2(40/3), ]$Nt.Karber <-  log2(40/3)

# Check results
round(cor(Karber_and_4PL$Nt_true, Karber_and_4PL$Nt.Karber, method = "spearman"), 2)
round(cor(Karber_and_4PL$Nt_true, Karber_and_4PL$Nt.4pl, method = "spearman"), 2)

plot(Karber_and_4PL$Nt_true, Karber_and_4PL$Nt.Karber)
abline(0, 1)
plot(Karber_and_4PL$Nt_true, Karber_and_4PL$Nt.4pl)
abline(0, 1)

rm(Karber_data, FRNT50, mean_vc)
rm(Karber_function, DRM_function)

# 3. Fit BHM Model ----
#### Estimate effect of each working virus stock -----
df_vc <- simu_data$df_vc
which_sto_per_bat <- df_vc %>% group_by(bat_id) %>% slice(1)

# Fit simulated data (virus control only)
data1_stan <- list(num_sto_vc = 4,
                   num_bat_vc = length(unique(df_vc$bat_id)),
                   num_obs_vc = nrow(df_vc),
                   which_sto_vc = df_vc$sto_id,
                   which_bat_vc = df_vc$bat_id,
                   count_vc = df_vc$Count,
                   
                   which_sto_per_bat = which_sto_per_bat$sto_id,
                   # Priors
                   delta_prior = c(48, 44, 40, 25),
                   sample_from_prior = 0)

# model <- stan_model(model_code=model1_code)
# fit1 <- stan(model_code=model1_code, data=data1_stan,
#              chains=4, cores=4, iter=3000,
#              control = list(adapt_delta = 0.95))
# saveRDS(fit1, file = "../3. Output/Simu_fit1.rds")
fit1 <- readRDS("../3. Output/Simu_fit1.rds")
posterior1 <- extract(fit1)

# Check results
plot(simu_data$sto_effect, apply(posterior1$delta, 2, median))
abline(0, 1)
plot(simu_data$bat_effect_vc, apply(posterior1$bat_vc, 2, median))
abline(0, 1)

# Sample from priors only 
data1_stan <- list(num_sto_vc = 4,
                   num_bat_vc = length(unique(df_vc$bat_id)),
                   num_obs_vc = nrow(df_vc),
                   which_sto_vc = df_vc$sto_id,
                   which_bat_vc = df_vc$bat_id,
                   count_vc = df_vc$Count,
                   
                   which_sto_per_bat = which_sto_per_bat$sto_id,
                   # Priors
                   delta_prior = c(48, 44, 40, 25),
                   sample_from_prior = 1)

# fit1_prior <- stan(model_code=model1_code, data=data1_stan,
#                    chains=4, cores=4, iter=3000,
#                    control = list(adapt_delta = 0.95))
# saveRDS(fit1_prior, file = "../3. Output/Simu_fit1_prior.rds")
fit1_prior <- readRDS("../3. Output/Simu_fit1_prior.rds")
posterior1_prior <- extract(fit1_prior)

## Calculate Foci reduction data (FR) ----
sto_param <- data.frame(sto_id = 1:4,
                        post_sto_effect = apply(posterior1$delta, 2, median))

df_PC <- left_join(simu_data$df_PC, sto_param)
df_PC$FR <- 1 - df_PC$Count/df_PC$post_sto_effect

df_IS1 <- left_join(simu_data$df_IS1, sto_param)
df_IS1$FR <- 1 - df_IS1$Count/df_IS1$post_sto_effect

df_IS2 <- left_join(simu_data$df_IS2, sto_param)
df_IS2$FR <- 1 - df_IS2$Count/df_IS2$post_sto_effect

df_sam <- left_join(df_sam, sto_param)
df_sam$FR <- 1 - df_sam$Count/df_sam$post_sto_effect

rm(sto_param)

## Remove 'clear negative' samples ----
pos_sam <- df_sam %>% filter(dilution == 40) %>% 
  group_by(simu_sam_id, bat_id) %>% 
  summarise(max_y = mean(FR)) %>% 
  filter(max_y > 0.3) %>% 
  arrange(simu_sam_id) 
pos_sam$sam_id <- 1:nrow(pos_sam)

df_sam <- df_sam %>% filter(simu_sam_id %in% pos_sam$simu_sam_id)
df_sam <- left_join(df_sam, pos_sam[, c("simu_sam_id", "sam_id")])

## Estimate batch effects and nAb titers ----
# Fit simulated data 
data2_stan <- list(num_bat_total = 28,
                   
                   num_obs_pc = nrow(df_PC),
                   which_bat_pc = df_PC$bat_id,
                   x_log_pc = log(df_PC$dilution),
                   y_pc = df_PC$FR,
                   
                   num_obs_IS1 = nrow(df_IS1),
                   which_bat_IS1 = df_IS1$bat_id,
                   x_log_IS1 = log(df_IS1$dilution),
                   y_IS1 = df_IS1$FR,
                   
                   num_obs_IS2 = nrow(df_IS2),
                   which_bat_IS2 = df_IS2$bat_id,
                   x_log_IS2 = log(df_IS2$dilution),
                   y_IS2 = df_IS2$FR,
                   
                   num_sam = length(unique(df_sam$sam_id)),
                   num_obs = nrow(df_sam),
                   which_bat = df_sam$bat_id,
                   which_sam = df_sam$sam_id,
                   x_log = log(df_sam$dilution),
                   y = df_sam$FR,
                   
                   f_pc_prior = c(1, 0, 1, 6),
                   f_IS1_prior = c(1, 0, 1, 6),   
                   f_IS2_prior = c(1, 0, 1, 6),   
                   f_pop_prior = c(1, 0, 1, 6),   
                   
                   f_lower = c(0, -1, 0, log(40)),
                   f_upper = c(2,  1, 1, log(87480)),
                   
                   which_bat_per_sam = pos_sam$bat_id,
                   
                   sample_from_prior = 0)

# model <- stan_model(model_code=model2_code)
# fit2 <- stan(model_code=model2_code, data=data2_stan,
#              chains=4, cores=4, iter=3000,
#              control = list(adapt_delta = 0.95))
# saveRDS(fit2, file = "../3. Output/Simu_fit2.rds")
fit2 <- readRDS("../3. Output/Simu_fit2.rds")
posterior2 <- extract(fit2)

# Sample from priors only 
data2_stan <- list(num_bat_total = 28,
                   
                   num_obs_pc = nrow(df_PC),
                   which_bat_pc = df_PC$bat_id,
                   x_log_pc = log(df_PC$dilution),
                   y_pc = df_PC$FR,
                   
                   num_obs_IS1 = nrow(df_IS1),
                   which_bat_IS1 = df_IS1$bat_id,
                   x_log_IS1 = log(df_IS1$dilution),
                   y_IS1 = df_IS1$FR,
                   
                   num_obs_IS2 = nrow(df_IS2),
                   which_bat_IS2 = df_IS2$bat_id,
                   x_log_IS2 = log(df_IS2$dilution),
                   y_IS2 = df_IS2$FR,
                   
                   num_sam = length(unique(df_sam$sam_id)),
                   num_obs = nrow(df_sam),
                   which_bat = df_sam$bat_id,
                   which_sam = df_sam$sam_id,
                   x_log = log(df_sam$dilution),
                   y = df_sam$FR,
                   
                   f_pc_prior = c(1, 0, 1, 6),
                   f_IS1_prior = c(1, 0, 1, 6),   
                   f_IS2_prior = c(1, 0, 1, 6),   
                   f_pop_prior = c(1, 0, 1, 6),   
                   
                   f_lower = c(0, -1, 0, log(40)),
                   f_upper = c(2,  1, 1, log(87480)),
                   which_bat_per_sam = pos_sam$bat_id,
                   
                   sample_from_prior = 1)
# 
# fit2_prior <- stan(model_code=model2_code, data=data2_stan,
#                    chains=4, cores=4, iter=3000,
#                    control = list(adapt_delta = 0.95))
# saveRDS(fit2_prior, file = "../3. Output/Simu_fit2_prior.rds")
fit2_prior <- readRDS("../3. Output/Simu_fit2_prior.rds")
posterior2_prior <- extract(fit2_prior)

# rm(model)
rm(data1_stan, data2_stan, model1_code, model2_code)

# 4. Figures ----
## Dist plot 1: effect of each working virus stock ----
delta <- extract(fit1, "delta")$delta
dimnames(delta) <- list(iteration = 1:6000, sto = 1:4)
delta <- reshape2::melt(delta, varnames = c("Iteration", "sto"))

delta_prior <- extract(fit1_prior, "delta")$delta
dimnames(delta_prior) <- list(iteration = 1:6000, sto = 1:4)
delta_prior <- reshape2::melt(delta_prior, varnames = c("Iteration", "sto"))

tiff("../3. Output/Dist plot 1_sto_effects.tiff", width = 8, height = 2.5, units = "in", res = 300)
ggplot() +
  facet_wrap(~ sto, nrow = 1, scales = "free_y",
             labeller = as_labeller(c("1" = "Virus working dilution = 200",
                                      "2" = "Virus working dilution = 300",
                                      "3" = "Virus working dilution = 330",
                                      "4" = "Virus working dilution = 400"))) +
  geom_density(data = delta %>% filter(Iteration > 3000),
               aes(x = value, color = "Posterior", fill =  "Posterior"), alpha = 0.5) +
  geom_density(data = delta_prior %>% filter(Iteration > 3000),
               aes(x = value, color = "Prior", fill = "Prior"), alpha = 0.5) +
  geom_vline(data = delta %>% filter(Iteration > 3000) %>%
               group_by(sto) %>% summarise(median_value = median(value)),
             aes(xintercept = median_value, color = "Posterior"),
             linetype = "dashed") +
  geom_vline(data = data.frame(sto = 1:4, delta = simu_data$sto_effect),
             aes(xintercept = delta, color = "True", fill = "True"),
             linetype = "dashed") +
  scale_x_continuous("Effects of each working virus stock") +
  scale_y_continuous("Density") +
  scale_color_manual("",values = c("Prior" = "#6BAED6", "Posterior" = "#D53E4F", "True" = "black")) +
  scale_fill_manual("", values = c("Prior" = "#9ECAE1", "Posterior" = "#FDAE61", "True" = "black")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = 'top')
dev.off()
rm(delta_prior, delta)

## Dist plot 2: population-level parameters ----
f_pop <- extract(fit2, "f_pop")$f_pop
dimnames(f_pop) <- list(iteration = 1:6000,  Pars = c("b", "c", "d", "e"))
f_pop <- reshape2::melt(f_pop, varnames = c("Iteration", "Pars"))

f_pop_prior <- extract(fit2_prior, "f_pop")$f_pop
dimnames(f_pop_prior) <- list(iteration = 1:6000,  Pars = c("b", "c", "d", "e"))
f_pop_prior <- reshape2::melt(f_pop_prior, varnames = c("Iteration", "Pars"))

tiff("../3. Output/Dist plot 2_pop_param.tiff", width = 8, height = 2.5, units = "in", res = 300)
ggplot() +
  scale_color_manual("",values = c("Prior" = "#6BAED6", "Posterior" = "#D53E4F", "True" = "black")) +
  scale_fill_manual("", values = c("Prior" = "#9ECAE1", "Posterior" = "#FDAE61", "True" = "black")) +
  facet_wrap(~Pars, nrow = 1, scales = "free",
             labeller = as_labeller(c(
               "b" = "italic(b^pop)",
               "c" = "italic(c^pop)",
               "d" = "italic(d^pop)",
               "e" = "italic(hat(e)^pop)"),
               label_parsed)) +
  geom_density(data = f_pop %>% filter(Iteration > 3000),
               aes(x = value, color =  "Posterior", fill =  "Posterior"), alpha = 0.5) +
  geom_density(data = f_pop_prior %>% filter(Iteration > 3000),
               aes(x = value, color =  "Prior", fill =  "Prior"), alpha = 0.5) +
  geom_vline(data = data.frame(Pars = c("b", "c", "d", "e"),
                               true = simu_data$f_pop),
             aes(xintercept = true, color = "True", fill = "True"),
             linetype = "dashed") +
  geom_vline(data = data.frame(f_pop %>% filter(Iteration > 3000) %>% 
                                 group_by(Pars) %>% summarise(post = median(value))),
             aes(xintercept = post, color = "Posterior", fill = "Posterior"),
             linetype = "dashed") +
  scale_x_continuous("Value") +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "top")
dev.off()
rm(f_pop_prior, f_pop)

## Dist plot 3: Batch level parameters ----
f_bat <- extract(fit2, "f_bat")$f_bat 
dimnames(f_bat) <- list(iteration = 1:6000, bat_id = 1:28, Pars = c("b", "c", "d", "e"))
f_bat <- reshape2::melt(f_bat, varnames = c("Iteration", "bat_id", "Pars"))

f_bat_prior <- extract(fit2_prior, "f_bat")$f_bat 
dimnames(f_bat_prior) <- list(iteration = 1:6000, bat_id = 1:28, Pars = c("b", "c", "d", "e"))
f_bat_prior <- reshape2::melt(f_bat_prior, varnames = c("Iteration", "bat_id", "Pars"))

f_bat_true <- simu_data$bat_effect
dimnames(f_bat_true) <- list(bat_id = 1:28, Pars = c("b", "c", "d", "e"))
f_bat_true <- reshape2::melt(f_bat_true, varnames = c("bat_id", "Pars"))

# select 6 batches for plot 
set.seed(123)
selected_batches <- sample(unique(f_bat_true$bat_id), 6)

tiff("../3. Output/Dist plot 3_batch_param.tiff", width = 12, height = 6, units = "in", res = 300)
ggplot() +
  scale_color_manual("",values = c("Prior" = "#6BAED6", "Posterior" = "#D53E4F", "True" = "black")) +
  scale_fill_manual("", values = c("Prior" = "#9ECAE1", "Posterior" = "#FDAE61", "True" = "black")) +
  facet_grid(
    Pars ~ bat_id, scales = "free",
    labeller = labeller(
      Pars = as_labeller(c(
        "b" = "italic(gamma[b])",
        "c" = "italic(gamma[c])",
        "d" = "italic(gamma[d])",
        "e" = "italic(gamma[hat(e)])"
      ), label_parsed),
      bat_id = label_value)) +
  geom_density(data = f_bat %>% filter(bat_id %in% selected_batches & Iteration > 3000),
               aes(x = value, color =  "Posterior", fill =  "Posterior"), alpha = 0.5) +
  geom_density(data = f_bat_prior %>% filter(bat_id %in% selected_batches & Iteration > 3000), 
               aes(x = value, color =  "Prior", fill =  "Prior"), alpha = 0.5) +
  geom_vline(data = f_bat_true %>% filter(bat_id %in% selected_batches),
             aes(xintercept = value, color = "True", fill = "True"),
             linetype = "dashed") +
  geom_vline(data = f_bat %>% filter(bat_id %in% selected_batches & Iteration > 3000) %>%
               group_by(Pars, bat_id) %>% summarise(post = median(value)),
             aes(xintercept = post, color = "Posterior", fill = "Posterior"),
             linetype = "dashed") +
  scale_x_continuous("Values", limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5)) +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "top")
dev.off()

## Dist plot 4: Sample level parameters ----
f_sam_bat <- extract(fit2, "f_sam_bat")$f_sam_bat
dimnames(f_sam_bat) <- list(iteration = 1:6000, 
                            sam_id = 1:nrow(pos_sam), 
                            Pars = c("b", "c", "d", "e"))
f_sam_bat <- reshape2::melt(f_sam_bat, 
                            varnames = c("Iteration", "sam_id", "Pars"))

f_sam_bat_prior <- extract(fit2_prior, "f_sam_bat")$f_sam_bat 
dimnames(f_sam_bat_prior) <- list(iteration = 1:6000, 
                                  sam_id = 1:nrow(pos_sam), 
                                  Pars = c("b", "c", "d", "e"))
f_sam_bat_prior <- reshape2::melt(f_sam_bat_prior, 
                                  varnames = c("Iteration", "sam_id", "Pars"))

f_sam_true <- simu_data$sam_effect
dimnames(f_sam_true) <- list(simu_sam_id = 1:length(unique(simu_data$df_sam$simu_sam_id)), 
                             Pars = c("b", "c", "d", "e"))
f_sam_true <- reshape2::melt(f_sam_true, varnames = c("simu_sam_id", "Pars"))

# select 6 samples for plot 
set.seed(123)
selected_samples <- sample(unique(pos_sam$sam_id), 6)
df <- pos_sam %>% filter(sam_id %in% selected_samples) %>%
  dplyr::select(simu_sam_id, sam_id, bat_id)
colnames(f_sam_bat)[4] <- "value_sam_bat"
colnames(f_bat)[4] <- "value_bat"
f_sam_bat <- left_join(df, f_sam_bat)
f_sam_bat <- left_join(f_sam_bat, f_bat)
f_sam_bat <- f_sam_bat %>% mutate(value = value_sam_bat - value_bat)

colnames(f_sam_bat_prior)[4] <- "value_sam_bat"
colnames(f_bat_prior)[4] <- "value_bat"
f_sam_bat_prior <- left_join(df, f_sam_bat_prior)
f_sam_bat_prior <- left_join(f_sam_bat_prior, f_bat_prior)
f_sam_bat_prior <- f_sam_bat_prior %>% 
  mutate(value = value_sam_bat - value_bat)

f_sam_true <- left_join(df, f_sam_true)
f_sam_true <- f_sam_true %>% filter(sam_id %in% selected_samples)

tiff("../3. Output/Dist plot 4_sample_param.tiff", width = 12, height = 6, units = "in", res = 300)
ggplot() +
  scale_color_manual("",values = c("Prior" = "#6BAED6", "Posterior" = "#D53E4F", "True" = "black")) +
  scale_fill_manual("", values = c("Prior" = "#9ECAE1", "Posterior" = "#FDAE61", "True" = "black")) +
  facet_grid(
    Pars ~ sam_id, scales = "free",
    labeller = labeller(
      Pars = as_labeller(c(
        "b" = "italic(lambda[b])",
        "c" = "italic(lambda[c])",
        "d" = "italic(lambda[d])",
        "e" = "italic(lambda[hat(e)])"
      ), label_parsed),
      bat_id = label_value)) +
  geom_density(data = f_sam_bat %>% filter(Iteration > 3000),
               aes(x = value, color =  "Posterior", fill =  "Posterior"), alpha = 0.5) +
  geom_density(data = f_sam_bat_prior %>% filter(Iteration > 3000), 
               aes(x = value, color =  "Prior", fill =  "Prior"), alpha = 0.5) +
  geom_vline(data = f_sam_true, aes(xintercept = value, color = "True", fill = "True"),
             linetype = "dashed") +
  geom_vline(data = f_sam_bat %>% filter(Iteration > 3000) %>%
               group_by(Pars, sam_id) %>% summarise(post = median(value)),
             aes(xintercept = post, color = "Posterior", fill = "Posterior"),
             linetype = "dashed") +
  scale_x_continuous("Values", limits = c(-5, 5), breaks = seq(-5, 5, 2.5)) +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "top")
dev.off()
rm(f_bat_prior, f_bat, f_bat_true, selected_batches)
rm(df, f_sam_bat_prior, f_sam_bat, f_sam_true, selected_samples)

## Main Fig 1: Agreement in titer estimates ----
# all nAb titer estimates 
pos_sam$Nt.unadj <- log2(exp(apply(posterior2$Nt_sam, 2, median)))
pos_sam$Nt.adj <- log2(exp(apply(posterior2$Nt_sam_adj, 2, median)))
df <- left_join(Karber_and_4PL, pos_sam)
df <- df %>% dplyr::select(simu_sam_id, sam_id, bat_id, Nt_true, 
                    Nt.Karber, Nt.4pl, Nt.unadj, Nt.adj)
df <- data.table(df)
df[Nt.unadj <= log2(40/3) | is.na(Nt.unadj), Nt.unadj := log2(40/3)]
df[Nt.adj <= log2(40/3) | is.na(Nt.adj), Nt.adj := log2(40/3)]

# Spearman correlation
rho1 <- round(cor(df$Nt_true, df$Nt.Karber, method = "spearman"), 2)
rho1_text <- bquote(Spearman~italic(rho)~"="~.(rho1))
rho2 <- round(cor(df$Nt_true, df$Nt.4pl, method = "spearman"), 2)
rho2_text <- bquote(Spearman~italic(rho)~"="~.(rho2))
rho3 <- round(cor(df$Nt_true, df$Nt.unadj, method = "spearman"), 2)
rho3_text <- bquote(Spearman~italic(rho)~"="~.(rho3))
rho4 <- round(cor(df$Nt_true, df$Nt.adj, method = "spearman"), 2)
rho4_text <- bquote(Spearman~italic(rho)~"="~.(rho4))

# Karber method
df %>% ggplot(aes(x = Nt_true)) + 
  geom_point(aes(y = Nt.Karber), color = "#2B83BA", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  scale_x_continuous("True titers", limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  scale_y_continuous("NAb titer estimates",limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  
  annotate("text",  x = 3, y = 15.8, hjust = 0, vjust = 1, label = "Karber Method", 
           size = 2.5, fontface = "bold") +
  annotate("text",  x = 11, y = 4, hjust = 0, vjust = 1, label = rho1_text, size = 2.5) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p1

# 4PL Model
df %>% ggplot(aes(x = Nt_true)) + 
  geom_point(aes(y = Nt.4pl), color = "#41AB5D", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  scale_x_continuous("True titers", limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  scale_y_continuous("NAb titer estimates",limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  
  annotate("text",  x = 3, y = 15.8, hjust = 0, vjust = 1, label = "4PL Model", 
           size = 2.5, fontface = "bold") +
  annotate("text",  x = 11, y = 4, hjust = 0, vjust = 1, label = rho2_text, size = 2.5) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p2

# BHM-Unadjusted
df %>% ggplot(aes(x = Nt_true)) + 
  geom_point(aes(y = Nt.unadj), color = "#FDAE61", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  scale_x_continuous("True titers", limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  scale_y_continuous("NAb titer estimates",limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  
  annotate("text",  x = 3, y = 15.8, hjust = 0, vjust = 1, label = "BHM-Unadjusted", 
           size = 2.5, fontface = "bold") +
  annotate("text",  x = 11, y = 4, hjust = 0, vjust = 1, label = rho3_text, size = 2.5) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p3

# BHM-Adjusted
df %>% ggplot(aes(x = Nt_true)) + 
  geom_point(aes(y = Nt.adj), color = "#D7191C", alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = log2(40), linetype = 2, color = "grey") +
  geom_vline(xintercept = log2(40), linetype = 2, color = "grey") +
  scale_x_continuous("True titers", limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  scale_y_continuous("NAb titer estimates",limits = c(2.5,16.5), breaks = seq(4, 16, 2),
                     labels = c(expression(2^4), expression(2^6), expression(2^8), 
                                expression(2^10), expression(2^12), expression(2^14),
                                expression(2^16))) +
  
  annotate("text",  x = 3, y = 15.8, hjust = 0, vjust = 1, label = "BHM-Adjusted", 
           size = 2.5, fontface = "bold") +
  annotate("text",  x = 11, y = 4, hjust = 0, vjust = 1, label = rho4_text, size = 2.5) +

  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p4

# Difference
df_diff <- df
df_diff[, 5:8] <- df[, 5:8] - df[, c(4,4,4,4)]
df_diff <- df_diff %>% 
  tidyr::pivot_longer(cols = -c("simu_sam_id", "sam_id", "bat_id", "Nt_true"), 
                      names_to = "Method", values_to = "diff")
df_diff$Method <- factor(df_diff$Method,
                         levels = c("Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj"))

ggplot(df_diff, aes(x = Method, y = diff, color = Method, fill = Method)) +
  scale_color_manual("",
                     breaks = c("Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj"),
                     labels = c("Karber Method", "4PL Model", 
                                "BHM-Unadjusted", "BHM-Adjusted"),
                     values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("",
                    breaks = c("Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj"),
                    labels = c("Karber Method", "4PL Model", 
                               "BHM-Unadjusted", "BHM-Adjusted"),
                    values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  geom_violin(alpha = 0.3, trim = FALSE, linewidth = 0.25) +
  geom_boxplot(width = 0.1, fill = "white", linewidth = 0.25, outlier.shape = NA) + 
  geom_hline(aes(yintercept = 0), color = "grey", linetype = "dashed") +
  scale_x_discrete("Methods", labels = c("Karber\nMethod",
                                         "4PL\nModel",
                                         "BHM\nUnadjusted",
                                         "BHM\nAdjusted")) +
  scale_y_continuous("Difference from true titers", limits = c(-6, 12),
                     breaks = seq(-6, 12, 2)) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7)) -> p5

# Output
tiff("../3. Output/Plot_simu_agreement.tiff", width = 9, height = 5 , units = "in", res = 600)
ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3, labels = c("A. ", "B. ", "C. ", "D. ", "E."),
          font.label = list(size = 10))
dev.off()
rm(p1, p2, p3, p4, p5, df_diff)
rm(rho1, rho2, rho3, rho4)
rm(rho1_text, rho2_text, rho3_text, rho4_text)
rm(Karber_and_4PL, pos_sam)

## Main Fig 2: illustrate simulated data ----
vc_param <- which_sto_per_bat[, c("sto_id", "bat_id", "mu_vc_per_bat")]
vc_param$median <- apply(posterior1$mu_vc_per_bat, 2, median)
vc_param$lower <- apply(posterior1$mu_vc_per_bat, 2, quantile, 0.025)
vc_param$upper <- apply(posterior1$mu_vc_per_bat, 2, quantile, 0.975)
colnames(vc_param)[3] <- "True" 

#### VC ----
ggplot(data = vc_param) + 
  scale_color_manual("Working virus stock",
                     breaks = 1:4, 
                     labels = c("A", "B", "C", "D"),
                     values = c("#FC8D62","#66C2A5","#999999","#8DA0CB")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_point(aes(x = True, y = median, color = factor(sto_id))) +
  geom_errorbar(data = vc_param, aes(x = True, ymin = lower, ymax =upper,
                                     color = factor(sto_id)), width = 0) +
  scale_x_continuous("True foci counts",
                     limits = c(25, 56),
                     breaks = seq(25, 55, 5)) +
  scale_y_continuous("Estimated foci counts\n", 
                     limits = c(23, 57),
                     breaks = seq(25, 55, 5)) +
  annotate("text", x = 52, y = 25, fontface = "bold",
           label = "Virus Control", size = 3) +
  theme_bw() + theme(legend.position = c(0.2, 0.7),
                     legend.background = element_blank(),
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 8),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7),
                     plot.margin = unit(c(5, 5, 5, 10), "point"),
                     panel.grid = element_blank()) -> p1

#### PC ----
# Unadjusted curve
df_PC$Method <- "Batch effects included"
df_PC$x_log <- log(df_PC$dilution)
df_PC$mu_y_true <- df_PC$mu_y
df_PC$mu_y_BHM <- apply(posterior2$mu_y_pc, 2, median)
df_PC$mu_y_BHM_lower <- apply(posterior2$mu_y_pc, 2,quantile, 0.025)
df_PC$mu_y_BHM_upper <- apply(posterior2$mu_y_pc, 2,quantile, 0.975)
unadj_curve <- df_PC %>% 
  dplyr::select(Method, bat_id, x_log, mu_y_true, mu_y_BHM, mu_y_BHM_lower, mu_y_BHM_upper)

# Adjusted curve
f_PC <- simu_data$f_PC
x_log <- unique(df_PC$x_log)
adj_curve <- rbind(data.frame(Method = "Batch effects excluded",
                              bat_id = 0,
                              x_log = x_log, 
                              mu_y_true = f_PC[2] + (f_PC[3] - f_PC[2]) / 
                                       (1 + exp(f_PC[1] * (x_log - f_PC[4]))),
                              mu_y_BHM = apply(posterior2$y_pc_adj, 2, median),
                              mu_y_BHM_lower = apply(posterior2$y_pc_adj, 2, quantile, 0.025),
                              mu_y_BHM_upper = apply(posterior2$y_pc_adj, 2, quantile, 0.975)))
PC_curve <- rbind(adj_curve, unadj_curve)

bat_sel <- c(0, 20, 21)
ggplot(PC_curve %>% filter(bat_id %in% bat_sel)) +
  geom_ribbon(aes(x = x_log, fill = Method, 
                  ymin = mu_y_BHM_lower, ymax = mu_y_BHM_upper,
                  group = factor(bat_id)), alpha = 0.2) +
  geom_line(aes(x = x_log, color = Method, y = mu_y_true, 
                linetype = "True", group = bat_id)) +
  geom_line(aes(x = x_log, color = Method, y = mu_y_BHM, 
                linetype = "Estimated", group = bat_id)) +
  scale_color_manual("", values = c("Batch effects excluded" = "#D7191C", 
                                    "Batch effects included" =  "#FDAE61")) +
  scale_fill_manual("", values = c("Batch effects excluded" = "#D7191C", 
                                   "Batch effects included" =  "#FDAE61")) +
  scale_linetype_manual("", values = c("True" = "solid", "Estimated" = "dashed")) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey", linewidth = 0.25) +
  geom_point(data = df_PC %>% filter(bat_id %in% bat_sel),
             aes(x = log(dilution), y = FR), 
             shape = 1, color =  "#FDAE61", size = 0.8, alpha = 0.8) +
  scale_x_continuous("Serial dilutions",
                      breaks = x_log, labels = unique(df_PC$dilution)) +
  scale_y_continuous("Foci reduction",
                      limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.25)) +
  annotate("text", x = log(50000), y = 0.63, label = "Batch 1", size = 2.5) +
  annotate("segment", x = log(50000), xend = log(50000), y = 0.56, yend = 0.47, 
           linewidth = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
  annotate("text", x = log(50000), y = -0.40, label = "Batch 2", size = 2.5) +
  annotate("segment", x = log(50000), xend = log(50000), y = -0.33, yend = -0.25, 
           linewidth = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
  annotate("text", x = log(20000), y = 0.94, label = "Positive Control", size = 3,
           fontface = "bold") +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = c(0.3, 0.27),
                     legend.key.height = unit(11, "point"),
                     legend.background = element_blank(),
                     legend.text = element_text(size = 7),
                     legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
                     legend.spacing.y = unit(-0.3, "cm"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 7)) +
  guides(linetype = guide_legend(order = 1),
         color = guide_legend(order = 2),
         fill = guide_legend(order = 2)) -> p2

#### Sample ----
x_log1 <- seq(x_log[1], x_log[8], 0.1)
true_f_pop <- simu_data$f_pop
true_sam_effect <- simu_data$sam_effect
df_sam$mu_y_BHM_unadj <- apply(posterior2$mu_y, 2, median)
mu_y_BHM_adj <- apply(posterior2$y_sam_adj, c(2, 3), median)

# select some representative samples for plot
selected_cases <- c(387, 407, 404) 
df1 <- data.frame(sam_id = selected_cases, id = 1:3)
df1 <- left_join(df1, df)
df2 <- data.frame(bat_id = unique(df1$bat_id), batch = 1:length(unique(df1$bat_id)))
df1 <- left_join(df1, df2)
df1$legend.true <- "True"
df1$legend.adj <- "BHM-Adjusted"
df1$legend.unadj <- "BHM-Unadjusted"
df1$legend.4pl <- "4PL Model"
df1$legend.karber <- "Karber Method"
df1$strip <- paste("Sample ", df1$id, ", Batch ", df1$batch, sep="")

# data for plot 
plot_data <- df_sam %>% filter(sam_id %in% selected_cases)
plot_data <- left_join(plot_data, df1)
mean_obs_vc <- df_vc %>% 
  filter(bat_id %in% unique(plot_data$bat_id)) %>% 
  group_by(bat_id) %>%
  summarise(mean_obs_vc = mean(Count))
plot_data <- left_join(plot_data, mean_obs_vc)
plot_data <- plot_data %>% mutate(FR1 = 1 - Count/mean_obs_vc)

for (j in selected_cases){
  k = df[df$sam_id == j, ]$bat_id
  tmp <- plot_data[plot_data$sam_id == j, ]
  
  # True
  j1 = df1[df1$sam_id == j, ]$simu_sam_id
  b <- as.numeric(true_f_pop[1] + true_sam_effect[j1, 1])
  c <- as.numeric(true_f_pop[2] + true_sam_effect[j1, 2])
  d <- as.numeric(true_f_pop[3] + true_sam_effect[j1, 3])
  e <- as.numeric(true_f_pop[4] + true_sam_effect[j1, 4])
  y_true.j <- c + (d - c) / (1 + exp(b*(x_log-e)))
  
  # 4PL model
  model <- drm(data = tmp, FR1 ~ dilution, fct = LL.4())
  FR_4pl.j <- predict(model, newdata = data.frame(dilution = exp(x_log1)))
  FR_4pl.j <- data.frame(sam_id = j, 
                         bat_id = k,
                         x_log = x_log1,
                         y_m = FR_4pl.j)
  # BHM-Adjusted
  BHM_adj.j <- data.frame(sam_id = j, 
                          bat_id = k,
                          x_log = x_log,
                          y_true = y_true.j,
                          y_m = mu_y_BHM_adj[, j])
  if (j == selected_cases[1]){
    FR_4pl <- FR_4pl.j
    BHM_adj <- BHM_adj.j
  } else{
    FR_4pl <- rbind(FR_4pl, FR_4pl.j)
    BHM_adj <- rbind(BHM_adj, BHM_adj.j)
  }
}

FR_4pl <- left_join(FR_4pl, df1)
BHM_adj <- left_join(BHM_adj, df1)

ggplot() + 
  facet_wrap(~strip, ncol = 3) +
  scale_color_manual("", values = c("True" = "black",
                                    "Karber Method" = "#2B83BA",
                                    "4PL Model" =  "#41AB5D",
                                    "BHM-Unadjusted" = "#FDAE61",
                                    "BHM-Adjusted" = "#D7191C")) +
  # True 
  geom_line(data = BHM_adj,
            aes(x = x_log, y = y_true, color = legend.true), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt_true), y = -0.32, color = legend.true),
             shape = 15, size = 1) +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt_true), xend = log(2^Nt_true),
                   y = -0.32, yend = 0.5, color = legend.true)) +
  
  # Karber Method
  geom_point(data = df1, aes(x = log(2^Nt.Karber), y = -0.32, color = legend.karber),
             shape = 6, size = 1) +
  
  # BHM-Unadjusted
  geom_point(data = plot_data %>% filter(sam_id %in% df1$sam_id),
             aes(x = log(dilution), y = FR), color = "#FDAE61",
             size = 0.8, alpha = 0.8, shape = 1) +
  geom_line(data = plot_data %>% filter(sam_id %in% df1$sam_id),
            aes(x = log(dilution), y = mu_y, color = legend.unadj), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt.unadj), y = -0.32, color = legend.unadj),
             shape = 7, size = 1) +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt.unadj), xend = log(2^Nt.unadj),
                   y = -0.32, yend = 0.5, color = legend.unadj)) +
  
  # 4PL Model
  geom_point(data = plot_data %>% filter(sam_id %in% df1$sam_id),
             aes(x = log(dilution), y = FR1), color = "#41AB5D",
             size = 0.8, alpha = 0.8, shape = 1) +
  geom_line(data = FR_4pl %>% filter(sam_id %in% df1$sam_id),
            aes(x = x_log, y = y_m, color = legend.4pl), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt.4pl), y = -0.32, color = legend.4pl),
             shape = 4, size = 1) +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt.4pl), xend = log(2^Nt.4pl),
                   y = -0.25, yend = 0.5, color = legend.4pl)) +
  
  # BHM-Adjusted
  geom_line(data = BHM_adj %>% filter(sam_id %in% df1$sam_id),
            aes(x = x_log, y = y_m, color = legend.adj), linewidth = 0.5) +
  geom_point(data = df1, aes(x = log(2^Nt.adj), y = -0.32, color = legend.adj),
             shape = 3, size = 1) +
  annotate("segment", x = log(40/3), xend = log(87480), y = 0.5, yend = 0.5, 
           linetype = "dashed", linewidth = 0.25, color = "grey") +
  geom_segment(data = df1, linetype = "dashed", linewidth = 0.25,
               aes(x = log(2^Nt.adj), xend = log(2^Nt.adj),
                   y = -0.32, yend = 0.5, color = legend.adj)) +
  
  scale_x_continuous("Serial dilutions",
                     breaks = x_log, labels = unique(plot_data$dilution)) +
  scale_y_continuous("Foci reduction", limits = c(-0.5, 1),
                     breaks = seq(-0.5, 1, 0.25)) + 
  
  theme_bw() + theme(strip.background = element_blank(),
                     panel.grid = element_blank(),
                     legend.position = "bottom",
                     legend.text = element_text(size = 10),
                     axis.text.x = element_text(angle = 25, size = 7, hjust = 1)) -> p3
#### Output ----
tiff("../3. Output/Plot_simu_curve.tiff", width = 7.5, height = 5.5 , units = "in", res = 600)
l1 <- ggarrange(p1, p2, nrow = 1, widths = c(1.35, 1), labels = c("A. ", "B. "),
                font.label = list(size = 12))
ggarrange(l1, p3, nrow = 2, labels = c("", "C. "), heights = c(1, 1), font.label = list(size = 12))
dev.off()


