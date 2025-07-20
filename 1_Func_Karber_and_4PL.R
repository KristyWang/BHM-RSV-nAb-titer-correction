# Load required libraries
library(drc)
library(dplyr)

# Function for calculating individual titers using Karber method ---- 
Karber_function <- function(data, 
                            sam_id_var_name = "sam_id", 
                            bat_id_var_name = "bat_id",
                            dilution_var_name = "dilution", 
                            Count_var_name = "Count", 
                            mean_vc_var_name = "mean_vc",
                            max_dilution = 87480,
                            dilution_interval = 3){
  
  Karber_data <- data.frame(sam_id = character(),
                            bat_id = numeric(),
                            Nt.Karber = numeric())

  data$sam_id <- data[, which(colnames(data) == sam_id_var_name)]
  data$bat_id <- data[, which(colnames(data) == bat_id_var_name)]
  data$dilution <- data[, which(colnames(data) == dilution_var_name)]
  data$Count <- data[, which(colnames(data) == Count_var_name)]
  data$mean_vc <- data[, which(colnames(data) == mean_vc_var_name)]
  
  for (j in unique(data$sam_id)){
    sample_j <- data %>% filter(sam_id == j)
    
    for (k in unique(sample_j$bat_id)){
      bat_k <- sample_j %>% filter(bat_id == k)
      
      # Other fixed variables not specified: mean_vc, dilution, Count
      p_k <-bat_k %>% group_by(dilution) %>% summarise(p = mean(Count/mean_vc))
      log10_ND50 <- log10(max_dilution) - log10(dilution_interval) * (sum(p_k$p) - 0.5)
      ND50 <- log2(10 ^ log10_ND50)
      tmp <- data.frame(sam_id = j,
                        bat_id = k,
                        Nt.Karber = ND50)
      Karber_data <- rbind(Karber_data, tmp)
    }
  }
  
  colnames(Karber_data)[1] <- sam_id_var_name
  colnames(Karber_data)[2] <- bat_id_var_name
  return(Karber_data)
  
}

# Function for calculating individual titers using 4PL model ---- 
DRM_function <- function(data, 
                         sam_id_var_name = "sam_id", 
                         bat_id_var_name = "bat_id",
                         dilution_var_name = "dilution", 
                         Count_var_name = "Count", 
                         mean_vc_var_name = "mean_vc",
                         max_dilution = 87480,
                         min_dilution = 40,
                         dilution_interval = 3){

  DRM_data <- data.frame(sam_id = character(),
                         bat_id = numeric(),
                         b = numeric(),
                         c = numeric(),
                         d = numeric(),
                         e = numeric(),
                         Nt.4pl = numeric())
  
  data$sam_id <- data[, which(colnames(data) == sam_id_var_name)]
  data$bat_id <- data[, which(colnames(data) == bat_id_var_name)]
  data$dilution <- data[, which(colnames(data) == dilution_var_name)]
  data$Count <- data[, which(colnames(data) == Count_var_name)]
  data$mean_vc <- data[, which(colnames(data) == mean_vc_var_name)]
  data$FR <- 1 - data$Count/data$mean_vc
  
  max_y <- data %>% filter(dilution == min_dilution) %>% 
    group_by(sam_id, bat_id) %>%
    summarise(max_y = mean(FR))
  pos_id <- max_y %>% filter(max_y >= 0.3)
  neg_id <- max_y %>% filter(max_y < 0.3)
  for (j in unique(pos_id$sam_id)){
    sample_j <- data %>% filter(sam_id == j)
    for (k in unique(sample_j$bat_id)){
      bat_k <- sample_j %>% filter(bat_id == k) %>% 
        group_by(dilution) %>%
        summarise(mean_FR = mean(FR))
      bat_k$dilution <- as.numeric(as.character(bat_k$dilution))
      tmp <- tryCatch({
        model <- drm(data = bat_k, mean_FR ~ dilution, 
                     fct = LL.4())
        summary(model)
        result <- summary(model)$coefficients[, 1]
        b <- result[1]
        c <- result[2]
        d <- result[3]
        e <- log(result[4])
        if ((d <= 0.5) | (b <= 0)){
          FRNT = log2(min_dilution) - log2(dilution_interval)
        } else if (c > 0.5){
          FRNT = log2(max_dilution)
        } else{
          FRNT = pmax(log2(exp(log((d-c)/(0.5-c)-1)/(b) + e)), log2(min_dilution) - log2(dilution_interval))
          FRNT = pmin(FRNT, log2(max_dilution))
        }
        
        tmp <- data.frame(sam_id = j,
                          bat_id = k,
                          b = b, c = c, d = d, e = e,
                          Nt.4pl = FRNT)
        tmp
      },
      error = function(cond){
        tmp <- data.frame(sam_id = j,
                          bat_id = k,
                          b = NA, c = NA, d = NA, e = NA,
                          Nt.4pl = NA)
        tmp
      })
      DRM_data <- rbind(DRM_data, tmp)
    }
  }
  
  tmp <- neg_id[, c(1,2)]
  tmp$b <- NA
  tmp$c <- NA
  tmp$d <- NA
  tmp$e <- NA
  tmp$Nt.4pl <- log2(min_dilution) - log2(dilution_interval)
  colnames(tmp) <- colnames(DRM_data)
  DRM_data <- rbind(DRM_data, tmp)
  
  colnames(DRM_data)[1] <- sam_id_var_name
  colnames(DRM_data)[2] <- bat_id_var_name
  return(DRM_data)
}