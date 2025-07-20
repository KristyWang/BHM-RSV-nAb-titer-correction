# Load required libraries
library(readxl)
library(data.table)
library(dplyr)

# Assay control ####
load("../1. Data/VC_raw.dat")
load("../1. Data/PC_raw.dat")
load("../1. Data/IS_raw.dat")

lab_rcd <- data.table(read_excel("../1. data/stock per batch.xlsx", sheet = 'Sheet2'))
lab_rcd$plate <- paste("plate", lab_rcd$plate_pc, sep=" ")

# VC 
colnames(VC_out) 
VC_out$Count <- as.numeric(VC_out$VC1)
VC_out$vd <- as.factor(VC_out$Virus_working_dilution)
VC_out <- left_join(VC_out, lab_rcd)
data_VC <- VC_out[, c("x", "vd", "plate","Count")]

# PC
PC_out$Serum_id <- "PC"
PC_out$Count <- as.numeric(PC_out$PC1)
PC_out$rpt <- NA
for (i in unique(PC_out$x)){
  num <- nrow(PC_out[PC_out$x == i, ])/8
  PC_out[PC_out$x == i, ]$rpt <- rep(1:num, each = 8)
}
PC_out$vd <- as.factor(PC_out$Virus_working_dilution)
PC_out <- left_join(PC_out, lab_rcd)
PC_out <- PC_out[, c("x", "vd", "Serum_id", 
                     "rpt", "dilution", "Count")]

# IS
IS_out$Serum_id <- paste("IS", IS_out$IS_concentration, sep="")
IS_out$Count <- as.numeric(IS_out$IS1)
IS_out$x1 <- paste(IS_out$x, IS_out$IS_concentration, sep="_")
IS_out$rpt <- NA
for (i in unique(IS_out$x1)){
  num <- nrow(IS_out[IS_out$x1 == i, ])/8
  IS_out[IS_out$x1 == i, ]$rpt <- rep(1:num, each = 8)
}
IS_out$vd <- as.factor(IS_out$Virus_working_dilution)
IS_out <- left_join(IS_out, lab_rcd)
IS_out <- IS_out[, c("x", "vd", "Serum_id", 
                     "rpt", "dilution", "Count")]
rm(i, num, VC_out, lab_rcd)


# Serum sample ----
data <- data.table(read_excel("../1. Data/Nt-raw-data-MN.xlsx", sheet = 'data'))

for (i in seq(1,(nrow(data)+1),12)){
  tmp <- data[seq(i,i+10,1),]
  test <- rep(as.character(tmp[1,1]),times=8)
  plate <- rep(as.character(tmp[1,2]),times=8)
  dilution <- tmp[4:11,1]
  
  id <- rep(as.character(tmp[2,3]),times=8)
  loc1 <- tmp[4:11,3]
  loc2 <- tmp[4:11,4]
  tmp1 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,5]),times=8)
  loc1 <- tmp[4:11,5]
  loc2 <- tmp[4:11,6]
  tmp2 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,7]),times=8)
  loc1 <- tmp[4:11,7]
  loc2 <- tmp[4:11,8]
  tmp3 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,9]),times=8)
  loc1 <- tmp[4:11,9]
  loc2 <- tmp[4:11,10]
  tmp4 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,11]),times=8)
  loc1 <- tmp[4:11,11]
  loc2 <- tmp[4:11,12]
  tmp5 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,13]),times=8)
  loc1 <- tmp[4:11,13]
  loc2 <- tmp[4:11,14]
  tmp6 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  if (i == 1) {data_out <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, use.names=F)}
  if (i != 1) {data_out <- rbind(data_out, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, use.names=F)}
  
  rm(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
}

data_out <- data_out[!is.na(id),]
colnames(data_out) <- c("test_date", "plate", "dilution", "Serum_id", "loc1", "loc2")
data_MN <- data_out
rm(data, data_out, dilution, loc1, loc2, tmp, i, id, plate, test)

#crs
raw <- data.table(read_excel("../1. Data/Nt-raw-data-crs.xlsx", sheet = "data"))

for (i in seq(1,(nrow(raw)+1),12)){
  tmp <- raw[seq(i,i+10,1),]
  test <- rep(as.character(tmp[1,1]),times=8)
  plate <- rep(as.character(tmp[1,2]),times=8)
  dilution <- tmp[4:11,1]
  
  id <- rep(as.character(tmp[2,3]),times=8)
  loc1 <- tmp[4:11,3]
  loc2 <- tmp[4:11,4]
  tmp1 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,5]),times=8)
  loc1 <- tmp[4:11,5]
  loc2 <- tmp[4:11,6]
  tmp2 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,7]),times=8)
  loc1 <- tmp[4:11,7]
  loc2 <- tmp[4:11,8]
  tmp3 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,9]),times=8)
  loc1 <- tmp[4:11,9]
  loc2 <- tmp[4:11,10]
  tmp4 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,11]),times=8)
  loc1 <- tmp[4:11,11]
  loc2 <- tmp[4:11,12]
  tmp5 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  id <- rep(as.character(tmp[2,13]),times=8)
  loc1 <- tmp[4:11,13]
  loc2 <- tmp[4:11,14]
  tmp6 <- cbind(test,plate,dilution, id, loc1, loc2)
  
  if (i == 1) {data_out <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, use.names=F)}
  if (i != 1) {data_out <- rbind(data_out, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, use.names=F)}
  
  rm(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
}

data_out <- data_out[!is.na(id),]
colnames(data_out) <- c("test_date", "plate", "dilution", "Serum_id", "loc1", "loc2")
data_cross <- data_out
data <- data.table(rbind(data_MN, data_cross))
rm(raw, data_out, dilution, loc1, loc2, tmp, i, id, plate, test)
rm(data_MN, data_cross)

data$x <- substr(data$test_date, nchar(data$test_date)-7, nchar(data$test_date))
data[data$x == "240411-1", ]$x <- "20240411-1"
data[data$x == "240411-2", ]$x <- "20240411-2"
data[data$test_date == "23_1_20240130", ]$x <- "20240130-1"
data[data$test_date == "23_2_20240130", ]$x <- "20240130-2"
data[data$test_date == "20240320", ]$x <- "20240314"

data <- left_join(data, data_VC %>% group_by(x) %>% slice(1) %>% dplyr::select(x, vd))

data1 <- data[, c("x", "vd", "Serum_id", "dilution", "loc1")]
data2 <- data[, c("x", "vd", "Serum_id", "dilution", "loc2")]
data1$rpt <- 1
data2$rpt <- 2
data <- rbind(data1, data2, use.names = F)
data$Count <- as.numeric(data$loc1)
data <- data[, c("x", "vd", "Serum_id", "rpt", "dilution", "Count")]
rm(data1, data2)

# Merge data ----
data <- rbind(PC_out, IS_out, data)
list <- data.frame(batch = 1:28, x = unique(data$x))
data <- left_join(data, list)
data$dilution <- as.numeric(data$dilution)
data_VC <- left_join(data_VC, list)
rm(IS_out, PC_out, list)

# Final dataset ---- 
data_VC <- data_VC[, c("batch", "vd", "Count")]
data <- data[, c("Serum_id", "batch", "rpt", "vd", "dilution", "Count")]
