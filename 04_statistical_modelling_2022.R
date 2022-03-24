# Statistical analysis script
# This script uses: 
# (1) the dataset of variables produced in "phylo_know_analysis_2021.R"
# (2) socioeconomic variables fitted to botanical countries
# (3) shapefile of botanical countries
#
# It conducts a statistical analysis of the dataset through spatial autoregressive modelling

# total runtime: 1 hour (due to AIC Weight calculations)




# ---- Assemble and check the data ---- #




# Load in required packages
library(spdep)
library(ncf)
library(spdplyr)
library(rgdal)
library(RColorBrewer)
library(tidyverse)
library(ggExtra)
library(spgwr)
library(e1071)
library(olsrr)
library(MASS)
library(lmtest)
library(car)
library(ggplot2)
library(data.table)


# Read in the dataset and the shapefile containing spatial references for the botanical countries (Without Bouvet Island) and merge
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

bio_var <- fread("data\\variables_2022.csv", quote = ",")
socioeco_var <- fread("data\\socioeco_var.csv", quote = ",")
data_var <- merge(bio_var, socioeco_var, by.x = "LEVEL_3_CO", by.y = "LEVEL_3_CO", all.x = T)
l3_shape <- readOGR("data\\level3.shp")
l3_shape <- l3_shape[-c(40),] # Excluding Bouvet Island
map01 <- merge(l3_shape, data_var, by.x = "LEVEL_3_CO", by.y = "LEVEL_3_CO", all.x = T)


# Make a correlation matrix to detect possible correlation between variables
cordf <- as.data.frame(data_var[,c(2:21)])
dfnum <- data.matrix(cordf)
resCOM <- cor(dfnum, use = "complete.obs", method = "spearman")
round(resCOM, 2)

# By looking at the correlation matrix, range_mean, range_median, endemism_total, gdp_sum and pop_count is excluded 

# Correlation matrix of included explanatory variables
cordf <- as.data.frame(data_var[,c(4,7,9,15,16,18:21)])
dfnum <- data.matrix(cordf)
resCOM <- cor(dfnum, use = "complete.obs", method = "kendall")
round(resCOM, 2) # Population Density and Road Density = 0.53, range and endemism = -0.41, range and sr = -0.36

# Check for heteroscedasticity
full_model_IC <- lm(RELATIVE ~ RICHNESS + RANGE_AREA + ENDEMISM_R + GDP_CAPITA + ROAD_DENSITY + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP, data = map01)
par(mfrow = c(2,2))
plot(full_model_IC)

# There seems to be homoscedasticity in the residuals, so the variables will not be transformed




# ---- Define neigborhood structure to control for SAC after testing multiple options ---- #



# define coordinates of botanical countries
coords <- coordinates(l3_shape)

# define distance-based neighboorhood structure
dist_nb_2000 <- dnearneigh(coords, 0, 2000, longlat=TRUE)

# Assign weights to the best neighborhood structure ---- #
dn2000w <- nb2listw(dist_nb_2000, zero.policy = T)





# Standardize the variables to enable comparison between variable slopes
map01$TOTAL_std <- scale(map01$TOTAL)
map01$RELATIVE_std <- scale(map01$RELATIVE)
map01$GEN_DUP_std <- scale(map01$GEN_DUP)
map01$RICHNESS_std <- scale(map01$RICHNESS)
map01$RANGE_MEDIAN_std <- scale(map01$RANGE_MEDIAN)
map01$RANGE_AREA_std <- scale(map01$RANGE_AREA)
map01$ENDEMISM_R_std <- scale(map01$ENDEMISM_R)
map01$BIEN_OCCUR_std <- scale(map01$BIEN_OCCUR)
map01$GDP_CAPITA_std <- scale(map01$GDP_CAPITA)
map01$ROAD_DENSITY_std <- scale(map01$ROAD_DENSITY)
map01$POP_DENSITY_std <- scale(map01$POP_DENSITY)
map01$SECURITY_std <- scale(map01$SECURITY)
map01$RESEARCH_EXP_std <- scale(map01$RESEARCH_EXP)
map01$EDUCATION_EXP_std <- scale(map01$EDUCATION_EXP)
map01$RICHNESS2_std <- map01$RICHNESS_std^2




# ---- Pick the best models by AIC using DNN = 2000 km as the neighborhood structure ---- #




# Create vectors of predictor variables (preds_qua includes quadratic term for sr)
preds <- c("RICHNESS_std" , "RANGE_AREA_std" , "ENDEMISM_R_std" , "GDP_CAPITA_std" , "ROAD_DENSITY_std" , "POP_DENSITY_std" , "SECURITY_std" , "RESEARCH_EXP_std" , "EDUCATION_EXP_std")
preds_qua <- c("RICHNESS_std" , "RICHNESS2_std" , "RANGE_AREA_std" , "ENDEMISM_R_std" , "GDP_CAPITA_std" , "ROAD_DENSITY_std" , "POP_DENSITY_std" , "SECURITY_std" , "RESEARCH_EXP_std" , "EDUCATION_EXP_std")

# Create lists that contains all possible combinations of the 9 predictors
allmodels <- list()
allmodels_qua <- list()

for(m in 1:9){
  combinations <- combn(preds, m = m)
  for(c in 1:ncol(combinations)){
    allmodels[[length(allmodels)+1]] <- combinations[,c]
  }
}
for(m in 1:10){
  combinations_qua <- combn(preds_qua, m = m)
  for(c in 1:ncol(combinations_qua)){
    allmodels_qua[[length(allmodels_qua)+1]] <- combinations_qua[,c]
  }
}



# -- Inventory completeness of phylogenetic knowledge as response variable -- #


# Make spatial error models for all combinations of predictors and extract the AIC of each model
AICs_IC <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_IC <- c(AICs_IC, AIC(model))
}

# Determine the variable combination that has the lowest AIC
sort(AICs_IC, index.return=TRUE, decreasing=F) # Lowest AIC belongs to IC variable combination 266

# dAIC > 2 = 406.1484 407.0816 407.2772 407.5421 407.5735 408.0124
# dAIC > 2 idx = 266 400 393 263 398 387

# Print variable combinations
allmodels[[266]]

# Based on AIC, create the best IC model
IC_model <- errorsarlm(RELATIVE_std ~ RICHNESS_std + RANGE_AREA_std + ENDEMISM_R_std + POP_DENSITY_std + RESEARCH_EXP_std , data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)



# -- BIEN as response variable -- #



# Make spatial error models for all combinations of predictors and extract the AIC of each model
AICs_BIEN <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_BIEN <- c(AICs_BIEN, AIC(model))
}

# Determine the variable combination that has the lowest AIC
sort(AICs_BIEN, index.return=TRUE, decreasing=F) # Lowest AIC belongs to BIEN variable combination 150

# dAIC > 2 = 788.9082 789.0484 789.3227 789.3926 789.9776 790.0059 790.4324 790.4735 790.6426 790.8164 790.9014 790.9067
# dAIC > 2 idx = 150 286 410 289 483 407 290 415 414 280 270 416

# Print variable combination 150
allmodels[[150]]

# Based on AIC, create the best BIEN model
BIEN_model <- errorsarlm(BIEN_OCCUR_std ~ RICHNESS_std + RANGE_AREA_std + RESEARCH_EXP_std + EDUCATION_EXP_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)



# -- Total count of sequences as response variable -- #



# Make spatial error models for all combinations of predictors and extract the AIC of each model
AICs_TOTAL <- c()

for(i in 1:length(allmodels_qua)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels_qua[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_TOTAL <- c(AICs_TOTAL, AIC(model))
}

# Determine the variable combination that has the lowest AIC
sort(AICs_TOTAL, index.return=TRUE, decreasing=F) # Lowest AIC belongs to TC variable combination 390

# dAIC > 2 = -116.958876 -116.133910 -115.749704 -115.210768 -115.174308 -115.130895 -115.062235
# dAIC > 2 idx = 390  648  853  650  645  652  641

# Print variable combination 390
allmodels_qua[[390]]

# Based on AIC, create the best TC model
TOTAL_model <- errorsarlm(TOTAL_std ~ RICHNESS_std + RICHNESS2_std + RANGE_AREA_std + ENDEMISM_R_std + RESEARCH_EXP, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)



# -- Inventory completeness with all markers as response variable -- #



AICs_markers <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_markers <- c(AICs_markers, AIC(model))
}

# Determine the variable combination that has the lowest AIC
sort(AICs_markers, index.return=TRUE, decreasing=F) # Lowest AIC belongs to IC variable combination 266

# dAIC > 2 = 338.1415 338.4114 338.4575 338.5203 339.4227 339.6523 339.7089 339.7932 339.8911 339.8946 340.0752
# dAIC > 2 idx = 398 266 476 393 480 263 507 395 472 502 400

# Print variable combinations
allmodels[[398]]

# Based on AIC, create the best IC model
markers_model <- errorsarlm(GEN_DUP_std ~ RICHNESS_std + RANGE_AREA_std + ENDEMISM_R_std + POP_DENSITY_std + RESEARCH_EXP_std + SECURITY_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)





# ---- Summaries of the best model for each response variable ---- #




# summaries including Nagelkerke Pseudo R-squared
summary(IC_model, Nagelkerke = T)
summary(BIEN_model, Nagelkerke = T)
summary(TOTAL_model, Nagelkerke = T)
summary(markers_model, Nagelkerke = T)




# ---- Check if residuals are still spatially correlated ---- #




# Add residuals to dataset
map01$RESID_IC_SAR <- residuals(IC_model)
map01$RESID_BIEN_SAR <- residuals(BIEN_model)

# detect spatial autocorrelation with a Monto Carlo simulation of Moran's I
moran.mc(map01$RESID_IC_SAR, dn2000w, nsim = 999, zero.policy = T)
moran.mc(na.omit(map01$RESID_BIEN_SAR), dn2000w, nsim = 999, zero.policy = T)

# Plot the residuals
par(mfrow = c(1,2))
spplot(map01, "RESID_IC_SAR")
spplot(map01, "RESID_BIEN_SAR")

# Residuals are no longer spatially correlated





# ---- Apply variance partitioning for species richness and mean range size for the best model ---- #


# Best model
summary(IC_model <- errorsarlm(RELATIVE_std ~ RICHNESS_std + RANGE_AREA_std + ENDEMISM_R_std + POP_DENSITY_std + RESEARCH_EXP_std , data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20), Nagelkerke = T)

# Best model without SR
summary(IC_model_noSR <- errorsarlm(RELATIVE_std ~ RANGE_AREA_std + ENDEMISM_R_std + POP_DENSITY_std + RESEARCH_EXP_std , data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20), Nagelkerke = T)

# Best model without Range
summary(IC_model_noRANGE <- errorsarlm(RELATIVE_std ~ RICHNESS_std + ENDEMISM_R_std + POP_DENSITY_std + RESEARCH_EXP_std , data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20), Nagelkerke = T)

Rt = 0.83051
Rrange = 0.80664
Rsr = 0.77345
RPrange = Rt - Rsr
RPsr = Rt - Rrange
Rmx = Rt - (Rrange + Rsr)

print(RPrange)
print(RPsr)
print(Rmx)




# ---- Determine thw weighted average coefficient for all variables for each model ---- #

 


# Summed weights
sum_weigth_IC <- c()
sum_weigth_BIEN <- c()
sum_weigth_TOTAL <- c()
sum_weigth_markers <- c()

for(i in AICs_IC){
  sum_weigth_IC <- sum(exp(-0.5*((i-406.1484)*-1)))
}

for(i in AICs_BIEN){
  sum_weigth_BIEN <- sum(exp(-0.5*((i-788.9082)*-1)))
}

for(i in AICs_TOTAL){
  sum_weigth_TOTAL <- sum(exp(-0.5*((i--116.958876))))
}

for(i in AICs_markers){
  sum_weigth_markers <- sum(exp(-0.5*((i-338.1415))))
}



# Akaike weights



AKAIKEW_IC <- c()
AKAIKEW_BIEN <- c()
AKAIKEW_TOTAL <- c()
AKAIKEW_markers <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AKAIKEW_IC[[i]] <- exp(-0.5*(AIC(model) - 406.1484))/sum_weigth_IC
}
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AKAIKEW_BIEN[[i]] <- exp(-0.5*(AIC(model) - 788.9082))/sum_weigth_BIEN
}
for(i in 1:length(allmodels_qua)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels_qua[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AKAIKEW_TOTAL[[i]] <- exp(-0.5*(AIC(model) - -116.958876))/sum_weigth_TOTAL
}
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AKAIKEW_markers[[i]] <- exp(-0.5*(AIC(model) - 338.1415))/sum_weigth_markers
}

AKAIKEW_IC[[266]]
AKAIKEW_BIEN[[150]]
AKAIKEW_TOTAL[[390]]
AKAIKEW_markers[[398]]



# -- Weighted mean slope for all variables -- #

# Inventory completeness weighted mean
weight_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  weight_IC[[i]] <- mean(as.list(coef(summary(model)))$RICHNESS_std)
}
weight_IC_0 <- unlist(weight_IC)
weight_IC_0[is.na(weight_IC_0)] <- 0
ic.mean.richness <- weighted.mean(weight_IC_0, unlist(AKAIKEW_IC))


range_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  range_IC[[i]] <- mean(as.list(coef(summary(model)))$RANGE_AREA_std)
}
range_IC_0 <- unlist(range_IC)
range_IC_0[is.na(range_IC_0)] <- 0
ic.mean.range <- weighted.mean(range_IC_0, unlist(AKAIKEW_IC))


endemism_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  endemism_IC[[i]] <- mean(as.list(coef(summary(model)))$ENDEMISM_R_std)
}
endemism_IC_0 <- unlist(endemism_IC)
endemism_IC_0[is.na(endemism_IC_0)] <- 0
ic.mean.endemism <- weighted.mean(endemism_IC_0, unlist(AKAIKEW_IC))


GDP_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  GDP_IC[[i]] <- mean(as.list(coef(summary(model)))$GDP_CAPITA_std)
}
GDP_IC_0 <- unlist(GDP_IC)
GDP_IC_0[is.na(GDP_IC_0)] <- 0
ic.mean.GDP <- weighted.mean(GDP_IC_0, unlist(AKAIKEW_IC))


road_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  road_IC[[i]] <- mean(as.list(coef(summary(model)))$ROAD_DENSITY_std)
}
road_IC_0 <- unlist(road_IC)
road_IC_0[is.na(road_IC_0)] <- 0
ic.mean.road <- weighted.mean(road_IC_0, unlist(AKAIKEW_IC))


pop_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  pop_IC[[i]] <- mean(as.list(coef(summary(model)))$POP_DENSITY_std)
}
pop_IC_0 <- unlist(pop_IC)
pop_IC_0[is.na(pop_IC_0)] <- 0
ic.mean.pop <- weighted.mean(pop_IC_0, unlist(AKAIKEW_IC))


sec_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  sec_IC[[i]] <- mean(as.list(coef(summary(model)))$SECURITY_std)
}
sec_IC_0 <- unlist(sec_IC)
sec_IC_0[is.na(sec_IC_0)] <- 0
ic.mean.sec <- weighted.mean(sec_IC_0, unlist(AKAIKEW_IC))


res_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  res_IC[[i]] <- mean(as.list(coef(summary(model)))$RESEARCH_EXP_std)
}
res_IC_0 <- unlist(res_IC)
res_IC_0[is.na(res_IC_0)] <- 0
ic.mean.res <- weighted.mean(res_IC_0, unlist(AKAIKEW_IC))


edu_IC <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  edu_IC[[i]] <- mean(as.list(coef(summary(model)))$EDUCATION_EXP_std)
}
edu_IC_0 <- unlist(edu_IC)
edu_IC_0[is.na(edu_IC_0)] <- 0
ic.mean.edu <- weighted.mean(edu_IC_0, unlist(AKAIKEW_IC))




# BIEN weighted mean




richness_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  richness_BIEN[[i]] <- mean(as.list(coef(summary(model)))$RICHNESS_std)
}
richness_BIEN_0 <- unlist(richness_BIEN)
richness_BIEN_0[is.na(richness_BIEN_0)] <- 0
BIEN.mean.richness <- weighted.mean(richness_BIEN_0, unlist(AKAIKEW_BIEN))


range_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  range_BIEN[[i]] <- mean(as.list(coef(summary(model)))$RANGE_AREA_std)
}
range_BIEN_0 <- unlist(range_BIEN)
range_BIEN_0[is.na(range_BIEN_0)] <- 0
BIEN.mean.range <- weighted.mean(range_BIEN_0, unlist(AKAIKEW_BIEN))


endemism_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  endemism_BIEN[[i]] <- mean(as.list(coef(summary(model)))$ENDEMISM_R_std)
}
endemism_BIEN_0 <- unlist(endemism_BIEN)
endemism_BIEN_0[is.na(endemism_BIEN_0)] <- 0
BIEN.mean.endemism <- weighted.mean(endemism_BIEN_0, unlist(AKAIKEW_BIEN))


GDP_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  GDP_BIEN[[i]] <- mean(as.list(coef(summary(model)))$GDP_CAPITA_std)
}
GDP_BIEN_0 <- unlist(GDP_BIEN)
GDP_BIEN_0[is.na(GDP_BIEN_0)] <- 0
BIEN.mean.GDP <- weighted.mean(GDP_BIEN_0, unlist(AKAIKEW_BIEN))


road_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  road_BIEN[[i]] <- mean(as.list(coef(summary(model)))$ROAD_DENSITY_std)
}
road_BIEN_0 <- unlist(road_BIEN)
road_BIEN_0[is.na(road_BIEN_0)] <- 0
BIEN.mean.road <- weighted.mean(road_BIEN_0, unlist(AKAIKEW_BIEN))


pop_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  pop_BIEN[[i]] <- mean(as.list(coef(summary(model)))$POP_DENSITY_std)
}
pop_BIEN_0 <- unlist(pop_BIEN)
pop_BIEN_0[is.na(pop_BIEN_0)] <- 0
BIEN.mean.pop <- weighted.mean(pop_BIEN_0, unlist(AKAIKEW_BIEN))


sec_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  sec_BIEN[[i]] <- mean(as.list(coef(summary(model)))$SECURITY_std)
}
sec_BIEN_0 <- unlist(sec_BIEN)
sec_BIEN_0[is.na(sec_BIEN_0)] <- 0
BIEN.mean.sec <- weighted.mean(sec_BIEN_0, unlist(AKAIKEW_BIEN))


res_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  res_BIEN[[i]] <- mean(as.list(coef(summary(model)))$RESEARCH_EXP_std)
}
res_BIEN_0 <- unlist(res_BIEN)
res_BIEN_0[is.na(res_BIEN_0)] <- 0
BIEN.mean.res <- weighted.mean(res_BIEN_0, unlist(AKAIKEW_BIEN))


edu_BIEN <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  edu_BIEN[[i]] <- mean(as.list(coef(summary(model)))$EDUCATION_EXP_std)
}
edu_BIEN_0 <- unlist(edu_BIEN)
edu_BIEN_0[is.na(edu_BIEN_0)] <- 0
BIEN.mean.edu <- weighted.mean(edu_BIEN_0, unlist(AKAIKEW_BIEN))



# Phylogenetic effort weighted mean



richness_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  richness_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$RICHNESS_std)
}
richness_TOTAL_0 <- unlist(richness_TOTAL)
richness_TOTAL_0[is.na(richness_TOTAL_0)] <- 0
TOTAL.mean.richness <- weighted.mean(richness_TOTAL_0, unlist(AKAIKEW_TOTAL))


richness2_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  richness2_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$RICHNESS2_std)
}
richness2_TOTAL_0 <- unlist(richness2_TOTAL)
richness2_TOTAL_0[is.na(richness2_TOTAL_0)] <- 0
TOTAL.mean.richness2 <- weighted.mean(richness2_TOTAL_0, unlist(AKAIKEW_TOTAL))


range_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  range_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$RANGE_AREA_std)
}
range_TOTAL_0 <- unlist(range_TOTAL)
range_TOTAL_0[is.na(range_TOTAL_0)] <- 0
TOTAL.mean.range <- weighted.mean(range_TOTAL_0, unlist(AKAIKEW_TOTAL))


endemism_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  endemism_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$ENDEMISM_R_std)
}
endemism_TOTAL_0 <- unlist(endemism_TOTAL)
endemism_TOTAL_0[is.na(endemism_TOTAL_0)] <- 0
TOTAL.mean.endemism <- weighted.mean(endemism_TOTAL_0, unlist(AKAIKEW_TOTAL))


GDP_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  GDP_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$GDP_CAPITA_std)
}
GDP_TOTAL_0 <- unlist(GDP_TOTAL)
GDP_TOTAL_0[is.na(GDP_TOTAL_0)] <- 0
TOTAL.mean.GDP <- weighted.mean(GDP_TOTAL_0, unlist(AKAIKEW_TOTAL))


road_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  road_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$ROAD_DENSITY_std)
}
road_TOTAL_0 <- unlist(road_TOTAL)
road_TOTAL_0[is.na(road_TOTAL_0)] <- 0
TOTAL.mean.road <- weighted.mean(road_TOTAL_0, unlist(AKAIKEW_TOTAL))


pop_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  pop_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$POP_DENSITY_std)
}
pop_TOTAL_0 <- unlist(pop_TOTAL)
pop_TOTAL_0[is.na(pop_TOTAL_0)] <- 0
TOTAL.mean.pop <- weighted.mean(pop_TOTAL_0, unlist(AKAIKEW_TOTAL))


sec_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  sec_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$SECURITY_std)
}
sec_TOTAL_0 <- unlist(sec_TOTAL)
sec_TOTAL_0[is.na(sec_TOTAL_0)] <- 0
TOTAL.mean.sec <- weighted.mean(sec_TOTAL_0, unlist(AKAIKEW_TOTAL))


res_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  res_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$RESEARCH_EXP_std)
}
res_TOTAL_0 <- unlist(res_TOTAL)
res_TOTAL_0[is.na(res_TOTAL_0)] <- 0
TOTAL.mean.res <- weighted.mean(res_TOTAL_0, unlist(AKAIKEW_TOTAL))


edu_TOTAL <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  edu_TOTAL[[i]] <- mean(as.list(coef(summary(model)))$EDUCATION_EXP_std)
}
edu_TOTAL_0 <- unlist(edu_TOTAL)
edu_TOTAL_0[is.na(edu_TOTAL_0)] <- 0
TOTAL.mean.edu <- weighted.mean(edu_TOTAL_0, unlist(AKAIKEW_TOTAL))



# Phylogenetic completeness of all markers weighted mean



richness_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  richness_markers[[i]] <- mean(as.list(coef(summary(model)))$RICHNESS_std)
}
richness_markers_0 <- unlist(richness_markers)
richness_markers_0[is.na(richness_markers_0)] <- 0
markers.mean.richness <- weighted.mean(richness_markers_0, unlist(AKAIKEW_markers))


range_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  range_markers[[i]] <- mean(as.list(coef(summary(model)))$RANGE_AREA_std)
}
range_markers_0 <- unlist(range_markers)
range_markers_0[is.na(range_markers_0)] <- 0
markers.mean.range <- weighted.mean(range_markers_0, unlist(AKAIKEW_markers))


endemism_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  endemism_markers[[i]] <- mean(as.list(coef(summary(model)))$ENDEMISM_R_std)
}
endemism_markers_0 <- unlist(endemism_markers)
endemism_markers_0[is.na(endemism_markers_0)] <- 0
markers.mean.endemism <- weighted.mean(endemism_markers_0, unlist(AKAIKEW_markers))


GDP_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  GDP_markers[[i]] <- mean(as.list(coef(summary(model)))$GDP_CAPITA_std)
}
GDP_markers_0 <- unlist(GDP_markers)
GDP_markers_0[is.na(GDP_markers_0)] <- 0
markers.mean.GDP <- weighted.mean(GDP_markers_0, unlist(AKAIKEW_markers))


road_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  road_markers[[i]] <- mean(as.list(coef(summary(model)))$ROAD_DENSITY_std)
}
road_markers_0 <- unlist(road_markers)
road_markers_0[is.na(road_markers_0)] <- 0
markers.mean.road <- weighted.mean(road_markers_0, unlist(AKAIKEW_markers))


pop_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  pop_markers[[i]] <- mean(as.list(coef(summary(model)))$POP_DENSITY_std)
}
pop_markers_0 <- unlist(pop_markers)
pop_markers_0[is.na(pop_markers_0)] <- 0
markers.mean.pop <- weighted.mean(pop_markers_0, unlist(AKAIKEW_markers))


sec_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  sec_markers[[i]] <- mean(as.list(coef(summary(model)))$SECURITY_std)
}
sec_markers_0 <- unlist(sec_markers)
sec_markers_0[is.na(sec_markers_0)] <- 0
markers.mean.sec <- weighted.mean(sec_markers_0, unlist(AKAIKEW_markers))


res_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  res_markers[[i]] <- mean(as.list(coef(summary(model)))$RESEARCH_EXP_std)
}
res_markers_0 <- unlist(res_markers)
res_markers_0[is.na(res_markers_0)] <- 0
markers.mean.res <- weighted.mean(res_markers_0, unlist(AKAIKEW_markers))


edu_markers <- c()
for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("GEN_DUP_std ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  edu_markers[[i]] <- mean(as.list(coef(summary(model)))$EDUCATION_EXP_std)
}
edu_markers_0 <- unlist(edu_markers)
edu_markers_0[is.na(edu_markers_0)] <- 0
markers.mean.edu <- weighted.mean(edu_markers_0, unlist(AKAIKEW_markers))



# ---- Print the results of the AIC weighted coefficient calculations ---- #


print(ic.mean.richness)
print(ic.mean.range)
print(ic.mean.endemism)
print(ic.mean.GDP)
print(ic.mean.road)
print(ic.mean.pop)
print(ic.mean.sec)
print(ic.mean.res)
print(ic.mean.edu)

print(BIEN.mean.richness)
print(BIEN.mean.range)
print(BIEN.mean.endemism)
print(BIEN.mean.GDP)
print(BIEN.mean.road)
print(BIEN.mean.pop)
print(BIEN.mean.sec)
print(BIEN.mean.res)
print(BIEN.mean.edu)

print(TOTAL.mean.richness)
print(TOTAL.mean.richness2)
print(TOTAL.mean.range)
print(TOTAL.mean.endemism)
print(TOTAL.mean.GDP)
print(TOTAL.mean.road)
print(TOTAL.mean.pop)
print(TOTAL.mean.sec)
print(TOTAL.mean.res)
print(TOTAL.mean.edu)

print(markers.mean.richness)
print(markers.mean.range)
print(markers.mean.endemism)
print(markers.mean.GDP)
print(markers.mean.road)
print(markers.mean.pop)
print(markers.mean.sec)
print(markers.mean.res)
print(markers.mean.edu)



# ---- Establish relationship between the two definitions of phylogenetic knowledge:



PK_model <- errorsarlm(RELATIVE_std ~ GEN_DUP_std , data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
summary(PK_model, Nagelkerke = T)

