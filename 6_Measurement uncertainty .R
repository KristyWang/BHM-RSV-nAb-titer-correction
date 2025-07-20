# This script is used to present the variability in RSV FRNT measurements 
# based on experimental data from assay controls, 
# including a virus control (VC), an internal standard control (PC), 
# and an International Standard control (IS, concentration = 1000). 
# The objectives are:
# 1) To characterize the variability in the raw experimental readouts (i.e., foci counts)
# and in foci reduction, which is calculated from the raw counts;
# 2) To quantify the contributions of within-batch and between-batch variability to titer estimates 
# obtained using standard methods, including the Karber formula and the 4-parameter logistic (4PL) model.

# Clear workspace
rm(list = ls())

# Load required libraries
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(lme4)
library(lmerTest)

# Source custom functions
source("1_Func_Karber_and_4PL.R")

# Import raw experimental data ----
source("5_Clean_experimental_data.R")
data$stock <- as.numeric(factor(data$vd, levels = c(200, 300, 330, 400), 1:4))
data_VC$stock <- as.numeric(factor(data_VC$vd, levels = c(200, 300, 330, 400), 1:4))
VC_data = data_VC[, c("batch", "stock", "Count")]
mean_obs_vc <- VC_data %>% group_by(batch, stock) %>% summarise(mean_obs_vc = mean(Count))
data <- left_join(data, mean_obs_vc)
data <- data %>% mutate(FR1 = 1 - Count/mean_obs_vc)
PC_data <- data[data$Serum_id == "PC", ]
IS_data <- data[data$Serum_id == "IS1000", ]
rm(data, data_VC, mean_obs_vc)

# Foci count: VC ----
ggplot() + 
  geom_boxplot(data = VC_data, linewidth = 0.2, width = 0.7,
               aes(x = factor(batch), y = Count, color = factor(stock)), 
               alpha = 0.5, fill = "white", outlier.shape = 1) +
  geom_point(data = VC_data %>% group_by(stock, batch) %>% summarise(mean= mean(Count)),
             aes(x = batch, y = mean, color = factor(stock))) +
  scale_color_manual("Virus working dilutions",
                     breaks = 1:4, 
                     labels = c("1:200", "1:300", "1:330", "1:400"),
                     values = c("#FC8D62","#66C2A5","#999999","#8DA0CB")) +
  scale_x_discrete("Batch id") +
  scale_y_continuous("Count of foci", limits = c(2.5, 72.5),
                     breaks = seq(10, 70, 10)) +
  annotate("text", x = 2.5, y = 5.5, label = "Virus Control", size = 2.5,
           fontface = "bold") +
  theme_bw() + theme(legend.position = c(0.5, 0.9),
                     legend.direction = "horizontal",
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 8),
                     plot.margin = unit(c(5, 10, 0, 12.5), "point"),
                     panel.grid = element_blank()) -> p0

# Foci count: PC ----
batch_sel <- c(20, 22)

PC_data %>% filter(batch %in% batch_sel) %>% mutate(group = paste0(batch, "-",rpt)) %>%
  ggplot(aes(x = log(dilution), y = Count, colour = factor(batch), group = group)) +
  geom_point() +
  geom_line(linetype = "dashed", alpha = 0.5, linewidth = 0.35) +
  scale_x_continuous("Serial dilutions",
                     breaks = log(unique(PC_data$dilution)), 
                     labels = unique(PC_data$dilution)) +
  scale_y_continuous("Count of foci") + 
  scale_color_manual("Batch id", values = c("#A65628", "#4DAF4A")) +
  annotate("text", x = log(27000), y = 2, label = "Positive Control", fontface = "bold",
           hjust = 0.5, vjust = 0.5, size = 2.5) +
  theme_bw() + theme(strip.background = element_blank(),
                     legend.position = c(0.15, 0.78),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     legend.key.height = unit(13, "point"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 8),
                     panel.grid = element_blank(),
                     plot.margin = unit(c(10, 10, 5, 10), "point"),
                     axis.text.x = element_text(angle = 0, hjust = 0.5)) -> p1

# Foci count: IS ----
IS_data %>% filter(batch %in% batch_sel) %>% mutate(group = paste0(batch, "-",rpt)) %>%
  ggplot(aes(x = log(dilution), y = Count, colour = factor(batch), group = group)) +
  geom_point() +
  geom_line(linetype = "dashed", alpha = 0.5, linewidth = 0.35) +
  scale_x_continuous("Serial dilutions",
                     breaks = log(unique(PC_data$dilution)), 
                     labels = unique(PC_data$dilution)) +
  scale_y_continuous("Count of foci") + 
  scale_color_manual("Batch id", values = c("#A65628", "#4DAF4A")) +
  annotate("text", x = log(17000), y = 3, label = "International Standard\n(Concentration = 1000)", 
           fontface = "bold", hjust = 0.5, vjust = 0.5, size = 2.5) +
  theme_bw() + theme(strip.background = element_blank(),
                     legend.position = c(0.15, 0.78),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     legend.key.height = unit(13, "point"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 8),
                     panel.grid = element_blank(),
                     plot.margin = unit(c(10, 10, 5, 10), "point"),
                     axis.text.x = element_text(angle = 0, hjust = 0.5)) -> p2

# Foci reduction: PC ----
model1 <- drm(data = PC_data %>% filter(batch == batch_sel[1]) %>%
                group_by(dilution) %>%
                summarise(mean_FR1 = mean(FR1)), mean_FR1 ~ dilution, fct = LL.4())
model1.coef <- data.frame(summary(model1)$coefficients)[, 1]
titer1 <- log2(exp((1/model1.coef[1]) * log((model1.coef[3] - model1.coef[2])/(0.5 - model1.coef[2]) - 1) + log(model1.coef[4])))
model1.1 <- drm(data = PC_data %>% filter(batch == batch_sel[1] & rpt == 1), FR1 ~ dilution, fct = LL.4())
model1.2 <- drm(data = PC_data %>% filter(batch == batch_sel[1] & rpt == 2), FR1 ~ dilution, fct = LL.4())

model2 <- drm(data = PC_data %>% filter(batch == batch_sel[2]) %>%
                group_by(dilution) %>%
                summarise(mean_FR1 = mean(FR1)), mean_FR1 ~ dilution, fct = LL.4())
model2.coef <- data.frame(summary(model2)$coefficients)[, 1]
titer2 <- log2(exp((1/model2.coef[1]) * log((model2.coef[3] - model2.coef[2])/(0.5 - model2.coef[2]) - 1) + log(model2.coef[4])))
model2.1 <- drm(data = PC_data %>% filter(batch == batch_sel[2] & rpt == 1), FR1 ~ dilution, fct = LL.4())
model2.2 <- drm(data = PC_data %>% filter(batch == batch_sel[2] & rpt == 2), FR1 ~ dilution, fct = LL.4())

x <- exp(seq(log(40), log(87580), 0.2))
newdata <- data.frame(dilution = x)
pred.data.all <- rbind(data.frame(batch = rep(batch_sel[1], times = 39), 
                                  dilution = x,
                                  fit = predict(model1, newdata = newdata)),
                       data.frame(batch = rep(batch_sel[2], times = 39), 
                                  dilution = x,
                                  fit = predict(model2, newdata = newdata)))

pred.data <- rbind(data.frame(batch = rep(batch_sel[1], times = 39), 
                              rpt = rep(1, times = 39),
                              dilution = x,
                              fit = predict(model1.1, newdata = newdata)),
                   data.frame(batch = rep(batch_sel[1], times = 39), 
                              rpt = rep(2, times = 39),
                              dilution = x,
                              fit = predict(model1.2, newdata = newdata)),
                   data.frame(batch = rep(batch_sel[2], times = 39), 
                              rpt = rep(1, times = 39),
                              dilution = x,
                              fit = predict(model2.1, newdata = newdata)),
                   data.frame(batch = rep(batch_sel[2], times = 39), 
                              rpt = rep(2, times = 39),
                              dilution = x,
                              fit = predict(model2.2, newdata = newdata)))
pred.data <- pred.data %>% mutate(group = paste0(batch, "-",rpt))

PC_data %>% filter(batch %in% batch_sel) %>% mutate(group = paste0(batch, "-",rpt)) %>%
  ggplot(aes(x = log(dilution), y = FR1, colour = factor(batch))) +
  geom_point(alpha = 0.3) +
  geom_line(data = pred.data, aes(group = group, y = fit), 
            linetype = "dashed", alpha = 0.5, linewidth = 0.35) +
  geom_line(data = pred.data.all, aes(x = log(dilution), y = fit, 
                                  group = factor(batch))) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  scale_x_continuous("Serial dilutions",
                     breaks = log(unique(PC_data$dilution)), 
                     labels = unique(PC_data$dilution)) +
  scale_y_continuous("Foci reduction", limits = c(-0.4, 1.1), 
                     breaks = seq(-0.4, 1, 0.2)) + 
  annotate("point", x = log(2^titer1), y = 0.5, color = "#A65628", size = 2) +
  annotate("segment", x = log(2^titer1), xend = log(2^titer1), y = -0.4, yend = 0.5,
           color = "#A65628", linewidth = 0.2, linetype = "dashed") +
  annotate("point", x = log(2^titer2), y = 0.5, color = "#4DAF4A", size = 2) +
  annotate("segment", x = log(2^titer2), xend = log(2^titer2), y = -0.4, yend = 0.5,
           color = "#4DAF4A", linewidth = 0.2, linetype = "dashed") +
  scale_color_manual("Batch id", values = c("#A65628", "#4DAF4A")) +
  scale_linetype_manual("Batch id", values = c(1,2)) +
  annotate("text", x = log(120), y = -0.36, label = "Positive Control", fontface = "bold",
           hjust = 0.5, vjust = 0.5, size = 2.5) +
  theme_bw() + theme(strip.background = element_blank(),
                     panel.grid = element_blank(),
                     legend.position = c(0.85, 0.78),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     legend.key.height = unit(13, "point"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 8),
                     plot.margin = unit(c(10, 10, 5, 5), "point"),
                     axis.text.x = element_text(angle = 0, hjust = 0.5)) -> p3

# Foci reduction: IS ----
model1 <- drm(data = IS_data %>% filter(batch == batch_sel[1]) %>%
                group_by(dilution) %>%
                summarise(mean_FR1 = mean(FR1)), mean_FR1 ~ dilution, fct = LL.4())
model1.coef <- data.frame(summary(model1)$coefficients)[, 1]
titer1 <- log2(exp((1/model1.coef[1]) * log((model1.coef[3] - model1.coef[2])/(0.5 - model1.coef[2]) - 1) + log(model1.coef[4])))
model1.1 <- drm(data = IS_data %>% filter(batch == batch_sel[1] & rpt == 1), FR1 ~ dilution, fct = LL.4())
model1.2 <- drm(data = IS_data %>% filter(batch == batch_sel[1] & rpt == 2), FR1 ~ dilution, fct = LL.4())

model2 <- drm(data = IS_data %>% filter(batch == batch_sel[2]) %>%
                group_by(dilution) %>%
                summarise(mean_FR1 = mean(FR1)), mean_FR1 ~ dilution, fct = LL.4())
model2.coef <- data.frame(summary(model2)$coefficients)[, 1]
titer2 <- log2(exp((1/model2.coef[1]) * log((model2.coef[3] - model2.coef[2])/(0.5 - model2.coef[2]) - 1) + log(model2.coef[4])))
model2.1 <- drm(data = IS_data %>% filter(batch == batch_sel[2] & rpt == 1), FR1 ~ dilution, fct = LL.4())
model2.2 <- drm(data = IS_data %>% filter(batch == batch_sel[2] & rpt == 2), FR1 ~ dilution, fct = LL.4())

x <- exp(seq(log(40), log(87580), 0.2))
newdata <- data.frame(dilution = x)
pred.data.all <- rbind(data.frame(batch = rep(batch_sel[1], times = 39), 
                                  dilution = x,
                                  fit = predict(model1, newdata = newdata)),
                       data.frame(batch = rep(batch_sel[2], times = 39), 
                                  dilution = x,
                                  fit = predict(model2, newdata = newdata)))

pred.data <- rbind(data.frame(batch = rep(batch_sel[1], times = 39), 
                              rpt = rep(1, times = 39),
                              dilution = x,
                              fit = predict(model1.1, newdata = newdata)),
                   data.frame(batch = rep(batch_sel[1], times = 39), 
                              rpt = rep(2, times = 39),
                              dilution = x,
                              fit = predict(model1.2, newdata = newdata)),
                   data.frame(batch = rep(batch_sel[2], times = 39), 
                              rpt = rep(1, times = 39),
                              dilution = x,
                              fit = predict(model2.1, newdata = newdata)),
                   data.frame(batch = rep(batch_sel[2], times = 39), 
                              rpt = rep(2, times = 39),
                              dilution = x,
                              fit = predict(model2.2, newdata = newdata)))
pred.data <- pred.data %>% mutate(group = paste0(batch, "-",rpt))

IS_data %>% filter(batch %in% batch_sel) %>% mutate(group = paste0(batch, "-",rpt)) %>%
  ggplot(aes(x = log(dilution), y = FR1, colour = factor(batch))) +
  geom_point(alpha = 0.3) +
  geom_line(data = pred.data, aes(group = group, y = fit), 
            linetype = "dashed", alpha = 0.5, linewidth = 0.35) +
  geom_line(data = pred.data.all, aes(x = log(dilution), y = fit, 
                                      group = factor(batch))) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  scale_x_continuous("Serial dilutions",
                     breaks = log(unique(IS_data$dilution)), 
                     labels = unique(IS_data$dilution)) +
  scale_y_continuous("Foci reduction", limits = c(-0.4, 1.1), 
                     breaks = seq(-0.4, 1, 0.2)) + 
  annotate("point", x = log(2^titer1), y = 0.5, color = "#A65628", size = 2) +
  annotate("segment", x = log(2^titer1), xend = log(2^titer1), y = -0.4, yend = 0.5,
           color = "#A65628", linewidth = 0.2, linetype = "dashed") +
  annotate("point", x = log(2^titer2), y = 0.5, color = "#4DAF4A", size = 2) +
  annotate("segment", x = log(2^titer2), xend = log(2^titer2), y = -0.4, yend = 0.5,
           color = "#4DAF4A", linewidth = 0.2, linetype = "dashed") +
  scale_color_manual("Batch id", values = c("#A65628", "#4DAF4A")) +
  scale_linetype_manual("Batch id", values = c(1,2)) +
  annotate("text", x = log(200), y = -0.305, label = "International Standard\n(Concentration = 1000)", 
           fontface = "bold", hjust = 0.5, vjust = 0.5, size = 2.5) +
  theme_bw() + theme(strip.background = element_blank(),
                     panel.grid = element_blank(),
                     legend.position = c(0.85, 0.78),
                     legend.text = element_text(size = 8),
                     legend.title = element_text(size = 10),
                     legend.key.height = unit(13, "point"),
                     axis.title = element_text(size = 10),
                     axis.text = element_text(size = 8),
                     plot.margin = unit(c(10, 10, 5, 5), "point"),
                     axis.text.x = element_text(angle = 0, hjust = 0.5)) -> p4

# Estimate titers for each assay replicate ----
PC_data <- PC_data %>% mutate(id = paste(batch, rpt, sep = "-"))
PC_Karber <- Karber_function(data = PC_data,
                             sam_id_var_name = "id",
                             bat_id_var_name = "batch",
                             dilution_var_name = "dilution", 
                             Count_var_name = "Count", 
                             mean_vc_var_name = "mean_obs_vc",
                             max_dilution = 87480,
                             dilution_interval = 3)

PC_4PL <- DRM_function(data = PC_data,
                       sam_id_var_name = "id",
                       bat_id_var_name = "batch",
                       dilution_var_name = "dilution", 
                       Count_var_name = "Count", 
                       mean_vc_var_name = "mean_obs_vc",
                       max_dilution = 87480,
                       min_dilution = 40,
                       dilution_interval = 3)

Karber_4PL_data <- left_join(PC_Karber, PC_4PL)
plot_data <- Karber_4PL_data %>%
  pivot_longer(cols = c(Nt.Karber, Nt.4pl), 
               names_to = "Method", values_to = "Value")
plot_data[plot_data$Method == "Nt.Karber", ]$Method <- "Karber Method"
plot_data[plot_data$Method == "Nt.4pl", ]$Method <- "4PL Model"
plot_data$Method <- factor(plot_data$Method, levels = c("Karber Method", "4PL Model"))

## Analysis of variance ----
#### karber method ----
model <- lmer(Nt.Karber ~ (1 | batch), data = Karber_4PL_data, REML = TRUE)
summary(model)
var_components <- as.data.frame(VarCorr(model))
var_between_batch <- var_components$vcov[1]  # Between-batch variance
var_within_batch <- var_components$vcov[2]  # within-batch variance
var_total <- var_between_batch  + var_within_batch
ratio_var_between_batch <- var_between_batch / var_total
ratio_var_within_batch <- var_within_batch / var_total

variance_data <- data.frame(
  Source = c("Between-batch", "Within-batch"),
  Variance = c(ratio_var_between_batch, ratio_var_within_batch)) %>% 
  mutate(label = paste0(Source, " (", round(100 * Variance, 1), "%)"))

ggplot(variance_data, aes(x = 1.5, y = Variance, fill = label)) + 
  geom_bar(stat = "identity", width = 1, color = "white") +  
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual("Source of variance", values = c("#FC9272","#FCBBA1")) +  
  xlim(c(0.1, 2)) +
  theme_void() +
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "point"),
        plot.margin = unit(c(5, 5, 0, 0), "point"),
        legend.direction = "vertical") -> p5.1

plot_data %>% filter(Method == "Karber Method") %>% ggplot() +
  geom_density(aes(x = Value, color = Method, fill = Method), alpha = 0.2) +
  scale_color_manual("Methods", values = c( "#FC9272", "#9ECAE1")) +   
  scale_fill_manual("Methods", values = c( "#FC9272", "#9ECAE1")) +
  scale_x_continuous(name = "NAb titers estimates",
                     limits = c(6.5, 14), 
                     breaks = seq(7, 14, 1),
                     labels = c(expression(2^7), expression(2^8), expression(2^9), 
                                expression(2^10), expression(2^11), expression(2^12), 
                                expression(2^13), expression(2^14))) +
  scale_y_continuous(name = "Density", limits = c(0, 0.65)) +
  annotate("text", x = 9.4, y = 0.162, size = 3,
           label = paste0("Total variance = ", round(var_total, 2))) +
  annotate("text", x = 6.5, y = 0.64, label = "Karber Method", size = 3, vjust = 1, hjust = 0, fontface = "bold") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(10, "point"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(5, 5, 0, 5), "point"),
        legend.position = "none") -> p5

p5 <- p5 + inset_element(p5.1, left = 0.53, bottom = 0.4, right = 0.99, top = 0.99)

#### 4PL Model ----
model <- lmer(Nt.4pl ~ (1 | batch), data = Karber_4PL_data, REML = TRUE)
summary(model)
var_components <- as.data.frame(VarCorr(model))
var_between_batch <- var_components$vcov[1]  # Between-batch variance
var_within_batch <- var_components$vcov[2]  # within-batch variance
var_total <- var_between_batch  + var_within_batch
ratio_var_between_batch <- var_between_batch / var_total
ratio_var_within_batch <- var_within_batch / var_total

variance_data <- data.frame(
  Source = c("Between-batch", "Within-batch"),
  Variance = c(ratio_var_between_batch, ratio_var_within_batch)) %>% 
  mutate(label = paste0(Source, " (", round(100 * Variance, 1), "%)"))

ggplot(variance_data, aes(x = 1.5, y = Variance, fill = label)) + 
  geom_bar(stat = "identity", width = 1, color = "white") +  
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual("Source of variance", values = c("#9ECAE1", "#C6DBEF")) +  
  xlim(c(0.1, 2)) +
  theme_void() +
  theme(legend.title = element_text(size = 8, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(10, "point"),
        plot.margin = unit(c(0,5,5,0), "point"),
        legend.direction = "vertical") -> p6.1

plot_data %>% filter(Method == "4PL Model") %>% ggplot() +
  geom_density(aes(x = Value, color = Method, fill = Method), alpha = 0.5) +
  scale_color_manual("Methods", values = c("#9ECAE1", "#FC9272")) +   
  scale_fill_manual("Methods", values = c( "#9ECAE1", "#FC9272")) +
  scale_x_continuous(name = "NAb titer estimates",
                     limits = c(6.5, 14), 
                     breaks = seq(7, 14, 1),
                     labels = c(expression(2^7), expression(2^8), expression(2^9), 
                                expression(2^10), expression(2^11), expression(2^12), 
                                expression(2^13), expression(2^14))) +
  scale_y_continuous(name = "Density", limits = c(0, 0.65)) +
  annotate("text", x = 9.3, y = 0.162, size = 3,
           label = paste0("Total variance = ", round(var_total, 2))) +
  annotate("text", x = 6.5, y = 0.64, label = "4PL Model", size = 3, vjust = 1, hjust = 0, fontface = "bold") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(10, "point"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(5, 5, 0, 5), "point"),
        legend.position = "none") -> p6

p6 <- p6 + inset_element(p6.1, left = 0.53, bottom = 0.37, right = 0.99, top = 0.96)

# Output ----
tiff("../3. Output/Plot_uncertainty in raw data.tiff", width = 7, height = 8, units = "in", res = 300)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, 
          labels = c("B. ", "C. ", "D.", "E."), 
          font.label = list(size = 10)) -> l2
ggarrange(p0, l2, ncol = 1, nrow = 2, heights = c(1, 2),
          labels = c("A. ", ""), 
          font.label = list(size = 10)) 

dev.off()

tiff("../3. Output/Plot_uncertainty in titer estimates.tiff", width = 7, height = 3, units = "in", res = 300)
ggarrange(p6, p5, ncol = 2, nrow = 1, 
          labels = c("A. ", "B. "), 
          font.label = list(size = 10))
dev.off()



