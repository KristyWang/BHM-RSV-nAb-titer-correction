# This script was used to assess how different titer estimation methods
# influence RSV serological measures at the population level, 
# including age-stratified seroprevalence, geometric mean titers (GMTs), fold-rises and seroconversion.

# Clear workspace
rm(list = ls())

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# 1. Import and merge data ----
## Import nAb titer estimates for each sample
nAb_titer <- data.table(read.csv("../3. Output/Final_Nt.csv"))[, -1]

## Import basic information for each sample
Basic_info <- data.table(read.csv("../1. Data/Basic_info.csv"))[, c(1,2,3,4,5,6,10)]

## Merge data 
nAb_titer[Serum_id == 60538, Serum_id := "6053C8"]
nAb_titer[Serum_id == 62878, Serum_id := "6287C8"]
nAb_titer[Serum_id == 62128, Serum_id := "6212C8"]
data_all <- data.table(left_join(nAb_titer, Basic_info))
data_all <- data_all[!is.na(age), ]
rm(nAb_titer, Basic_info)

## Define age groups
data_all <- data_all[, -c(2, 8)]
data_all[age < 7/12, age_grp := "0-6m"]
data_all[age >= 7/12 & age < 13/12, age_grp := "7-12m"]
data_all[age >= 13/12 & age < 3, age_grp := "1-2y"]
data_all[age >= 3 & age < 6, age_grp := "3-5y"]
data_all[age >= 6 & age < 18, age_grp := "6-17y"]
data_all[age >= 18 & age < 65, age_grp := "18-64y"]
data_all[age >= 65, age_grp := "65+y"]

data_all$age_grp <- factor(data_all$age_grp, 
                           levels = c("0-6m", "7-12m", "1-2y", 
                                      "3-5y", "6-17y","18-64y", "65+y"))
# Sample size by group
table(data_all$age_grp)

# 2. GMTs and seroprevalence ---- 
## Transform data
data_all_long <- data_all %>% 
  tidyr::pivot_longer(cols = -c("Serum_id", "batch", "kid", "mid", "age", 
                                "Visit", "sex", "source", "age_grp"),
                      names_to = "method", values_to = "value")

data_all_long$age_grp <- factor(data_all_long$age_grp, 
                                levels = c("0-6m", "7-12m", "1-2y", 
                                           "3-5y", "6-17y","18-64y", "65+y"))
data_all_long$method <- factor(data_all_long$method, 
                               levels = c("Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj"))

## GMTs 
data_all_long %>% group_by(age_grp, method) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  mutate(lower = mean - 1.96 * sd /sqrt(n),
         upper = mean + 1.96 * sd /sqrt(n)) %>% 
  ggplot(aes(x = age_grp, color = method)) +
  geom_point(aes(y = mean), position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(0.4), width = 0) + 
  geom_hline(yintercept = log2(40), linetype = "dashed", color = "grey") +
  scale_color_manual("",
                     breaks = c("Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj"),
                     labels = c("Karber Method", "4PL Model", 
                                "BHM-Unadjusted", "BHM-Adjusted"),
                     values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age groups") +
  scale_y_continuous("\nGeometric Mean Titers", limits = c(3.5,9.2), breaks = seq(4, 9, 1),
                     labels = c(expression(2^4), expression(2^5), expression(2^6), 
                                expression(2^7), expression(2^8), expression(2^9))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 7)) -> p1

## Seroprevalence
data_all_long$pos <- 0
data_all_long[data_all_long$value >= log2(40), ]$pos <- 1
summary_df <- data_all_long %>% group_by(method, age_grp) %>%
  summarise(n = n(), pos_count = sum(pos), .groups = "drop") %>%
  mutate(pos_rate = pos_count/n) %>% 
  mutate(A = 2 * n * pos_rate + 1.96 ^ 2,
         B = 1.96 * sqrt(1.96 ^ 2 + 4 * n * pos_rate * (1 - pos_rate)), 
         C = 2 * (n + 1.96 ^ 2)) %>% 
  mutate(ci_lower = (A - B) / C,
         ci_upper = (A + B) / C)

ggplot(summary_df, aes(x = age_grp, y = pos_rate, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(color = method), position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = method),
                position = position_dodge(width = 0.9), width = 0.2) +
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
  scale_x_discrete("Age groups") +
  scale_y_continuous("\nSeroprevalence (%)", 
                     limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 7)) -> p2

## Output
tiff("../3. Output/Plot_GMTs and seroprevalence.tiff", width = 6, height = 5, units = "in", res = 600)
ggarrange(p1, p2, nrow = 2, align = "hv",
          labels = c(" ", " "),
          font.label = list(size = 8), 
          common.legend = T, legend = "top")
dev.off()

rm(p1, p2, summary_df, data_all_long)

# 3. New infection ----
## Transform data
table(data_all$Visit)
colnames(data_all)
df1 <- data_all[data_all$Visit %in% c("C1", "C2", "C3", "C4", "C5", "C6"),
                c("kid", "Visit", "age", "Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj")]
colnames(df1) <- c("kid","Visit1", "age1", "Nt.Karber1", "Nt.4PL1", "Nt.unadj1", "Nt.adj1")

df2 <- data_all[data_all$Visit %in% c("C2", "C3", "C4", "C5", "C6", "C7"),
                c("kid", "Visit", "age", "Nt.Karber", "Nt.4pl", "Nt.unadj", "Nt.adj")]
colnames(df2) <- c("kid","Visit2", "age2", "Nt.Karber2", "Nt.4PL2", "Nt.unadj2", "Nt.adj2")

df1 <- data.table(df1)
df1[Visit1 == "C1", Visit2 := "C2"]
df1[Visit1 == "C2", Visit2 := "C3"]
df1[Visit1 == "C3", Visit2 := "C4"]
df1[Visit1 == "C4", Visit2 := "C5"]
df1[Visit1 == "C5", Visit2 := "C6"]
df1[Visit1 == "C6", Visit2 := "C7"]
df <- left_join(df1, df2)

## Define new infection outcomes ----
df <- df %>% mutate(sero.conv.Karber = (Nt.Karber2 >= log2(40) & Nt.Karber1 < log2(40)),
                    sero.2fld.Karber = ((sero.conv.Karber == F) & ((2^Nt.Karber2) / (2^Nt.Karber1) >= 2)),
                    sero.4fld.Karber = ((sero.conv.Karber == F) & ((2^Nt.Karber2) / (2^Nt.Karber1) >= 4)),
                    sero.conv.4PL = (Nt.4PL2 >= log2(40) & Nt.4PL1 < log2(40)),
                    sero.2fld.4PL = ((sero.conv.4PL == F) & ((2^Nt.4PL2) / (2^Nt.4PL1) >= 2)),
                    sero.4fld.4PL = ((sero.conv.4PL == F) & ((2^Nt.4PL2) / (2^Nt.4PL1) >= 4)),
                    sero.conv.unadj = (Nt.unadj2 >= log2(40) & Nt.unadj1 < log2(40)),
                    sero.2fld.unadj = ((sero.conv.unadj == F) & ((2^Nt.unadj2) / (2^Nt.unadj1) >= 2)),
                    sero.4fld.unadj = ((sero.conv.unadj == F) & ((2^Nt.unadj2) / (2^Nt.unadj1) >= 4)),
                    sero.conv.adj = (Nt.adj2 >= log2(40) & Nt.adj1 < log2(40)),
                    sero.2fld.adj = ((sero.conv.adj == F) & ((2^Nt.adj2) / (2^Nt.adj1) >= 2)),
                    sero.4fld.adj = ((sero.conv.adj == F) & ((2^Nt.adj2) / (2^Nt.adj1) >= 4)))
df <- df %>% filter(!is.na(age2))

## Population-level agreements ----
df_total <- df
df_total$Visit2 <- "Overall"
df <- rbind(df, df_total)
df.Karber <- df %>% group_by(Visit2) %>% 
  summarise(n = n(),
            n.conv = sum(sero.conv.Karber),
            n.2fld = sum(sero.2fld.Karber),
            n.4fld = sum(sero.4fld.Karber)) %>% 
  mutate(Method = "Karber Method")

df.4PL <- df %>% group_by(Visit2) %>% 
  summarise(n = n(),
            n.conv = sum(sero.conv.4PL),
            n.2fld = sum(sero.2fld.4PL),
            n.4fld = sum(sero.4fld.4PL)) %>% 
  mutate(Method = "4PL Model")

df.unadj <- df %>% group_by(Visit2) %>% 
  summarise(n = n(),
            n.conv = sum(sero.conv.unadj),
            n.2fld = sum(sero.2fld.unadj),
            n.4fld = sum(sero.4fld.unadj)) %>% 
  mutate(Method = "BHM-Unadjusted")

df.adj <- df %>% group_by(Visit2) %>% 
  summarise(n = n(),
            n.conv = sum(sero.conv.adj),
            n.2fld = sum(sero.2fld.adj),
            n.4fld = sum(sero.4fld.adj)) %>% 
  mutate(Method = "BHM-Adjusted")
df_sum <- rbind(df.Karber, df.4PL, df.unadj, df.adj)

df_sum <- df_sum %>%
  mutate(p.conv = n.conv / n,
         p.conv.low = p.conv - 1.96 * sqrt(p.conv * (1 - p.conv) / n),
         p.conv.high = p.conv + 1.96 * sqrt(p.conv * (1 - p.conv) / n),
         p.2fld = n.2fld / n, 
         p.2fld.low = p.2fld - 1.96 * sqrt(p.2fld * (1 - p.2fld) / n),
         p.2fld.high = p.2fld + 1.96 * sqrt(p.2fld * (1 - p.2fld) / n),
         p.4fld = n.4fld / n,
         p.4fld.low = p.4fld - 1.96 * sqrt(p.4fld * (1 - p.4fld) / n),
         p.4fld.high = p.4fld + 1.96 * sqrt(p.4fld * (1 - p.4fld) / n))

df_sum$Method <- factor(df_sum$Method, 
                        levels = c("Karber Method", 
                                   "4PL Model",
                                   "BHM-Unadjusted",
                                   "BHM-Adjusted"))
#### Seroconvension ----
df_sum %>% ggplot(aes(x = Visit2, y = p.conv, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(colour = Method), position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = p.conv.low, ymax = p.conv.high, colour = Method), 
                position = position_dodge(0.9), width = 0) +
  scale_color_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age at follow-up visits", labels = c("2m", "4m", "6m", "1y", "2y", "3y", "Overall")) +
  scale_y_continuous("\nSeroconversion (%)", 
                     limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1),
                     labels = seq(0, 60, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7)) -> p3

#### 2-Fold titer rise ----
df_sum %>% ggplot() +
  geom_bar(aes(x = Visit2, y = p.2fld, fill = Method),
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(x = Visit2, y = p.2fld, colour = Method), 
             position = position_dodge(0.9)) +
  geom_errorbar(aes(x = Visit2, ymin = p.2fld.low, ymax = p.2fld.high, colour = Method), 
                position = position_dodge(0.9), width = 0) +
  scale_color_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age at follow-up visits", labels = c("2m", "4m", "6m", "1y", "2y", "3y", "Overall")) +
  scale_y_continuous("\nTwo-fold titer rise (%)", 
                     limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1),
                     labels = seq(0, 60, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7)) -> p4

#### 4-Fold titer rise ----
df_sum %>% ggplot() +
  geom_bar(aes(x = Visit2, y = p.4fld, fill = Method),
           stat = "identity", position = position_dodge(width = 0.9), alpha = 0.2) +
  geom_point(aes(x = Visit2, y = p.4fld, colour = Method), 
             position = position_dodge(0.9)) +
  geom_errorbar(aes(x = Visit2, ymin = p.4fld.low, ymax = p.4fld.high, colour = Method), 
                position = position_dodge(0.9), width = 0) +
  scale_color_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_fill_manual("", values = c("#2B83BA", "#41AB5D", "#FDAE61", "#D7191C")) +
  scale_x_discrete("Age at follow-up visits", labels = c("2m", "4m", "6m", "1y", "2y", "3y", "Overall")) +
  scale_y_continuous("\nFour-fold titer rise (%)", 
                     limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1),
                     labels = seq(0, 60, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 7)) -> p5

#### Output ----
tiff("../3. Output/Plot_new_inf_pop_agreements.tiff", width = 10, height = 2.5, units = "in", res = 600)
ggarrange(p3, p4, p5, nrow = 1, ncol = 3, common.legend = T, 
          labels = c("A. ", "B. ", "C. "),
          font.label = list(size = 8))
dev.off()

## Individual-level agreements ----
#### Sample size  ----
df <- df %>% filter(Visit2 != "Overall")
df %>% summarise(n1 = sum(sero.conv.adj), n2 = n() - n1, n = n()) -> n.conv
df %>% summarise(n1 = sum(sero.2fld.adj), n2 = n() - n1, n = n()) -> n.2fld
df %>% summarise(n1 = sum(sero.4fld.adj), n2 = n() - n1, n = n()) -> n.4fld
n.samplesize <- c(t(n.conv), t(n.2fld), t(n.4fld))
#### Karber Method ----
n1 <- nrow(df[(df$sero.conv.adj == T & df$sero.conv.Karber == T), ])
n2 <- nrow(df[(df$sero.conv.adj == F & df$sero.conv.Karber == F), ])
n = n1 + n2
n.conv <- c(n1, n2, n)

n1 <- nrow(df[(df$sero.2fld.adj == T & df$sero.2fld.Karber == T), ])
n2 <- nrow(df[(df$sero.2fld.adj == F & df$sero.2fld.Karber == F), ])
n = n1 + n2
n.2fld <- c(n1, n2, n)

n1 <- nrow(df[(df$sero.4fld.adj == T & df$sero.4fld.Karber == T), ])
n2 <- nrow(df[(df$sero.4fld.adj == F & df$sero.4fld.Karber == F), ])
n = n1 + n2
n.4fld <- c(n1, n2, n)
n.Karber <- c(t(n.conv), t(n.2fld), t(n.4fld))
Agree.Karber <- n.Karber/n.samplesize

#### 4PL Model ----
n1 <- nrow(df[(df$sero.conv.adj == T & df$sero.conv.4PL == T), ])
n2 <- nrow(df[(df$sero.conv.adj == F & df$sero.conv.4PL == F), ])
n = n1 + n2
n.conv <- c(n1, n2, n)

n1 <- nrow(df[(df$sero.2fld.adj == T & df$sero.2fld.4PL == T), ])
n2 <- nrow(df[(df$sero.2fld.adj == F & df$sero.2fld.4PL == F), ])
n = n1 + n2
n.2fld <- c(n1, n2, n)

n1 <- nrow(df[(df$sero.4fld.adj == T & df$sero.4fld.4PL == T), ])
n2 <- nrow(df[(df$sero.4fld.adj == F & df$sero.4fld.4PL == F), ])
n = n1 + n2
n.4fld <- c(n1, n2, n)
n.4PL <- c(t(n.conv), t(n.2fld), t(n.4fld))
Agree.4PL <- n.4PL/n.samplesize

#### BHM-unadjusted ----
n1 <- nrow(df[(df$sero.conv.adj == T & df$sero.conv.unadj == T), ])
n2 <- nrow(df[(df$sero.conv.adj == F & df$sero.conv.unadj == F), ])
n = n1 + n2
n.conv <- c(n1, n2, n)

n1 <- nrow(df[(df$sero.2fld.adj == T & df$sero.2fld.unadj == T), ])
n2 <- nrow(df[(df$sero.2fld.adj == F & df$sero.2fld.unadj == F), ])
n = n1 + n2
n.2fld <- c(n1, n2, n)

n1 <- nrow(df[(df$sero.4fld.adj == T & df$sero.4fld.unadj == T), ])
n2 <- nrow(df[(df$sero.4fld.adj == F & df$sero.4fld.unadj == F), ])
n = n1 + n2
n.4fld <- c(n1, n2, n)
n.unadj <- c(t(n.conv), t(n.2fld), t(n.4fld))
Agree.unadj <- n.unadj/n.samplesize

Agreement <- data.frame(n = n.samplesize, Agree.Karber = Agree.Karber, 
                        Agree.4PL = Agree.4PL, Agree.unadj = Agree.unadj)
Sen <- Agreement[c(1,4,7), ] %>% 
  mutate(outcome = c("sero.conv", "sero.2fld", "sero.4fld")) %>% 
  tidyr::pivot_longer(cols = -c("n", "outcome"), 
                      names_to = "Method", 
                      values_to = "Sen") %>% 
  mutate(lci = Sen - 1.96 * sqrt(Sen * (1 - Sen) / n),
         uci = Sen + 1.96 * sqrt(Sen * (1 - Sen) / n))

Spe <- Agreement[c(2,5,8), ] %>% 
  mutate(outcome = c("sero.conv", "sero.2fld", "sero.4fld")) %>% 
  tidyr::pivot_longer(cols = -c("n", "outcome"), 
                      names_to = "Method", 
                      values_to = "Spe") %>%
  mutate(lci = Spe - 1.96 * sqrt(Spe * (1 - Spe) / n),
         uci = Spe + 1.96 * sqrt(Spe * (1 - Spe) / n))

overall.agreement <- Agreement[c(3,6,9), ] %>% 
  mutate(outcome = c("sero.conv", "sero.2fld", "sero.4fld")) %>% 
  tidyr::pivot_longer(cols = -c("n", "outcome"), 
                      names_to = "Method", 
                      values_to = "Agr") %>%
  mutate(lci = Agr - 1.96 * sqrt(Agr * (1 - Agr) / n),
         uci = Agr + 1.96 * sqrt(Agr * (1 - Agr) / n))

#### Output ----
write.csv(Agreement, file = "../3. Output/Table_new_inf_ind_agreements.csv")