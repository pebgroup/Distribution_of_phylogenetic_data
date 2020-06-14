# Statistical analysis script
# This script uses a dataset ("all_variables_2303.csv") compiled of biodiversity variables produced in the Biodiversity Variables script and socioeconomic variables produced in ArcMap
# It conducts a statistical analysis of the dataset through spatial autoregressive modelling




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


# Read in the dataset and the shapefile containing spatial references for the botanical countries (Without Bouvet Island)
data_var <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\results\\VAR_230320\\all_variables_2303.csv", header = T)
l3_shape <- readOGR("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\shapefile\\level3.shp")
l3_shape <- l3_shape[-c(40),] # Excluding Bouvet Island


# Merge with shapefile and visualize
map01 <- merge(l3_shape, data_var, by.x = "LEVEL_3_CO", by.y = "LEVEL_3_CO", all.x = T)
color_palette <- brewer.pal(n = 9, name = "YlOrRd") 
spplot(map01, "RELATIVE", col.regions = color_palette, cuts = 8, col = "gray", main = "Inventory completeness", xlab = "", ylab ="")


# Make a correlation matrix to detect possible correlation between variables
cordf <- as.data.frame(data_var[,c(4:19)])
dfnum <- data.matrix(cordf)
resCOM <- cor(dfnum, use = "complete.obs", method = "kendall")
round(resCOM, 2)
# By looking at the correlation matrix, range_mean, endemism_total, gdp_sum and pop_count is excluded 


# Check for heteroscedasticity
full_model_IC <- lm(RELATIVE ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R + GDP_CAPITA + ROAD_DENSITY + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP, data = map01)
par(mfrow = c(1,1))
plot(full_model_IC)

# There seems to be homoscedasticity in the residuals, so the variables will not be transformed




# ---- Define potential neigborhood structure ---- #




# Define a variety of different neighborhood structures to control for SAC
coords <- coordinates(l3_shape)

queen <- poly2nb(l3_shape)
rook <- poly2nb(l3_shape, queen = F)
rows <- row.names(as(l3_shape, "data.frame"))
kn1 <- knn2nb(knearneigh(coords, k = 1), row.names = rows)
kn2 <- knn2nb(knearneigh(coords, k = 2), row.names = rows)
kn3 <- knn2nb(knearneigh(coords, k = 3), row.names = rows)
dist_nb_1000 <- dnearneigh(coords, 0, 1000, longlat=TRUE)
dist_nb_1500 <- dnearneigh(coords, 0, 1500, longlat=TRUE)
dist_nb_2000 <- dnearneigh(coords, 0, 2000, longlat=TRUE)
dist_nb_5000 <- dnearneigh(coords, 0, 5000, longlat=TRUE)

# Repeat the process with a subset shapefile that excludes polygons for which there is missing data. This enables the later creation of correlograms for the residuals of the full linear model. 

# Subset the shapefile
incomplete_cases <- data_var[!complete.cases(data_var),]
incomplete_cases <- c(incomplete_cases$LEVEL_3_CO)
alt_l3 <- l3_shape[-c(incomplete_cases),]
alt_coords <- coordinates(alt_l3)

# Neighborhood structures only containing polygons with complete rows of data
queen_alt <- poly2nb(alt_l3)
rook_alt <- poly2nb(alt_l3, queen = F)
rows_alt <- row.names(as(alt_l3, "data.frame"))
kn1_alt <- knn2nb(knearneigh(alt_coords, k = 1), row.names = rows_alt)
kn2_alt <- knn2nb(knearneigh(alt_coords, k = 2), row.names = rows_alt)
kn3_alt <- knn2nb(knearneigh(alt_coords, k = 3), row.names = rows_alt)
dist_nb_1000_alt <- dnearneigh(alt_coords, 0, 1000, longlat=TRUE)
dist_nb_1500_alt <- dnearneigh(alt_coords, 0, 1500, longlat=TRUE)
dist_nb_2000_alt <- dnearneigh(alt_coords, 0, 2000, longlat=TRUE)
dist_nb_5000_alt <- dnearneigh(alt_coords, 0, 5000, longlat=TRUE)


# Plot and visualize the different neighborhoods
par(mfrow = c(2,4))
plot(l3_shape)
plot(queen, coords, add = T, main = "Queen")
plot(l3_shape)
plot(rook, coords, add = T, main = "Rook")
plot(l3_shape)
plot(kn1, coords, add = T, main = "KNN = 1")
plot(l3_shape)
plot(kn2, coords, add = T, main = "KNN = 2")
plot(l3_shape)
plot(kn3, coords, add = T, main = "KNN = 3")
plot(l3_shape)
plot(dist_nb_1000, coords, add = T, main = "1000km")
plot(l3_shape)
plot(dist_nb_1500, coords, add = T, main = "1500km")
plot(l3_shape)
plot(dist_nb_2000, coords, add = T, main = "2000km")




# ---- Produce correlograms of the different neighborhoods ---- #




# Full Inventory Completeness (IC) model and IC model with only complete cases 
full_model_IC <- lm(RELATIVE ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R + GDP_CAPITA + ROAD_DENSITY + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP, data = map01)
CC_model_IC <- lm(RELATIVE ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R + GDP_CAPITA + ROAD_DENSITY + POP_DENSITY, data = map01)

# Full IC model correlograms
par(mfrow = c(2,4))
cgram_DN1000f <- sp.correlogram(dist_nb_1000_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN1000f, main = "Correlogram DNN 1000 km")
cgram_DN1500f <- sp.correlogram(dist_nb_1500_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN1500f, main = "Correlogram DNN 1500 km")
cgram_DN2000f <- sp.correlogram(dist_nb_2000_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN2000f, main = "Correlogram DNN 2000 km")
cgram_DN5000f <- sp.correlogram(dist_nb_5000_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN5000f, main = "Correlogram DNN 5000 km")
cgram_rookf <- sp.correlogram(rook_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_rookf, main = "Correlogram rook")
cgram_queenf <- sp.correlogram(queen_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_queenf, main = "Correlogram queen")
cgram_KN2f <- sp.correlogram(kn2_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_KN2f, main = "Correlogram KNN 2")
cgram_KN3f <- sp.correlogram(kn3_alt, full_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_KN3f, main = "Correlogram KNN 3")
# DNN = 2000 km is best

# Subset complete cases IC model correlograms
cgram_DN1000f <- sp.correlogram(dist_nb_1000, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN1000f, main = "Correlogram DNN 1000 km")
cgram_DN1500f <- sp.correlogram(dist_nb_1500, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN1500f, main = "Correlogram DNN 1500 km")
cgram_DN2000f <- sp.correlogram(dist_nb_2000, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN2000f, main = "Correlogram DNN 2000 km")
cgram_DN5000f <- sp.correlogram(dist_nb_5000, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN5000f, main = "Correlogram DNN 5000 km")
cgram_rookf <- sp.correlogram(rook, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_rookf, main = "Correlogram rook")
cgram_queenf <- sp.correlogram(queen, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_queenf, main = "Correlogram queen")
cgram_KN2f <- sp.correlogram(kn2, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_KN2f, main = "Correlogram KNN 2")
cgram_KN3f <- sp.correlogram(kn3, CC_model_IC$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_KN3f, main = "Correlogram KNN 3")
# DNN = 2000 km is still best

# Full BIEN model and BIEN model with only complete cases 
full_model_BIEN <- lm(BIEN_OCCUR ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R + GDP_CAPITA + ROAD_DENSITY + POP_DENSITY + SECURITY + RESEARCH_EXP + EDUCATION_EXP, data = map01)
CC_model_BIEN <- lm(BIEN_OCCUR ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R + GDP_CAPITA + ROAD_DENSITY + POP_DENSITY, data = map01)

# Full BIEN model correlograms
cgram_DN1000f <- sp.correlogram(dist_nb_1000_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_DN1000f, main = "Correlogram DNN 1000 km")
cgram_DN1500f <- sp.correlogram(dist_nb_1500_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_DN1500f, main = "Correlogram DNN 1500 km")
cgram_DN2000f <- sp.correlogram(dist_nb_2000_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_DN2000f, main = "Correlogram DNN 2000 km")
cgram_DN5000f <- sp.correlogram(dist_nb_5000_alt, full_model_BIEN$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN5000f, main = "Correlogram DNN 5000 km")
cgram_rookf <- sp.correlogram(rook_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_rookf, main = "Correlogram rook")
cgram_queenf <- sp.correlogram(queen_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_queenf, main = "Correlogram queen")
cgram_KN2f <- sp.correlogram(kn2_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_KN2f, main = "Correlogram KNN 2")
cgram_KN3f <- sp.correlogram(kn3_alt, full_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_KN3f, main = "Correlogram KNN 3")
# DNN = 2000 is best

# Subset complete cases BIEN model correlograms
par(mfrow = c(2,4))
cgram_DN1000f <- sp.correlogram(dist_nb_1000, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_DN1000f, main = "Correlogram DNN 1000 km")
cgram_DN1500f <- sp.correlogram(dist_nb_1500, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_DN1500f, main = "Correlogram DNN 1500 km")
cgram_DN2000f <- sp.correlogram(dist_nb_2000, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_DN2000f, main = "Correlogram DNN 2000 km")
cgram_DN5000f <- sp.correlogram(dist_nb_5000, CC_model_BIEN$residuals, order = 5, method = "I", zero.policy = T)
plot(cgram_DN5000f, main = "Correlogram DNN 5000 km")
cgram_rookf <- sp.correlogram(rook, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_rookf, main = "Correlogram rook")
cgram_queenf <- sp.correlogram(queen, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_queenf, main = "Correlogram queen")
cgram_KN2f <- sp.correlogram(kn2, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_KN2f, main = "Correlogram KNN 2")
cgram_KN3f <- sp.correlogram(kn3, CC_model_BIEN$residuals, order = 10, method = "I", zero.policy = T)
plot(cgram_KN3f, main = "Correlogram KNN 3")
# DNN = 2000 is best 




# ---- Assign weights to the best neighborhood structure and decide on the best spatial autoregressive models ---- #




# Assign weight
dn2000w <- nb2listw(dist_nb_2000, zero.policy = T)

# Test for spatial autocorrelation (SAC)
lm.morantest(full_model_IC, dn2000w, zero.policy = T) # Plenty of SAC
lm.morantest(full_model_BIEN, dn2000w, zero.policy = T) # Plenty of SAC

# Lagrange multiplier test
LRM_IC <-  lm.LMtests(full_model_IC, dn2000w, test = "all", zero.policy = T)
LRM_IC
LRM_BIEN <-  lm.LMtests(full_model_IC, dn2000w, test = "all", zero.policy = T)
LRM_BIEN

# a Spatial Error model is preferred




# ---- Pick the best models by AIC using DNN = 2000 km as the neighborhood structure ---- #




# Create vector of predictor variables
preds <- c("RICHNESS" , "RANGE_MEDIAN" , "ENDEMISM_R" , "GDP_CAPITA" , "ROAD_DENSITY" , "POP_DENSITY" , "SECURITY" , "RESEARCH_EXP" , "EDUCATION_EXP")

# Create a list that contains all possible combinations of the 9 predictors
allmodels <- list()

for(m in 1:9){
  combinations <- combn(preds, m = m)
  for(c in 1:ncol(combinations)){
    allmodels[[length(allmodels)+1]] <- combinations[,c]
  }
}


# -- IC as response variable -- #


# Make spatial error models for all combinations of predictors and extract the AIC of each model
AICs_IC <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("RELATIVE ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_IC <- c(AICs_IC, AIC(model))
  }

# Determine the variable combination that has the lowest AIC
sort(AICs_IC, index.return=TRUE, decreasing=F) # Lowest AIC belongs to IC variable combination 46

# Print variable combination 46
allmodels[[46]]

# Based on AIC, create the best IC model
IC_model <- errorsarlm(RELATIVE ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)


# -- BIEN as response variable -- #


# Make spatial error models for all combinations of predictors and extract the AIC of each model
AICs_BIEN <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("BIEN_OCCUR ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_BIEN <- c(AICs_BIEN, AIC(model))
}

# Determine the variable combination that has the lowest AIC
sort(AICs_BIEN, index.return=TRUE, decreasing=F) # Lowest AIC belongs to BIEN variable combination 277

# Print variable combination 277
allmodels[[277]]

# Based on AIC, create the best BIEN model
BIEN_model <- errorsarlm(BIEN_OCCUR ~ RICHNESS + RANGE_MEDIAN + GDP_CAPITA + POP_DENSITY + EDUCATION_EXP, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)


# -- Total count of sequences as response variable -- #


# Make spatial error models for all combinations of predictors and extract the AIC of each model
AICs_TOTAL <- c()

for(i in 1:length(allmodels)){
  model <- errorsarlm(paste("TOTAL ~", paste(allmodels[[i]], collapse = " + ")), data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
  AICs_TOTAL <- c(AICs_TOTAL, AIC(model))
}

# Determine the variable combination that has the lowest AIC
sort(AICs_TOTAL, index.return=TRUE, decreasing=F) # Lowest AIC belongs to TC variable combination 401

# Print variable combination 401
allmodels[[401]]

# Based on AIC, create the best TC model
TOTAL_model <- errorsarlm(TOTAL ~ RICHNESS + RANGE_MEDIAN + ENDEMISM_R + SECURITY + RESEARCH_EXP + EDUCATION_EXP, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)





# ---- Summaries of the best model for each response variable ---- #



# summaries including Nagelkerke Pseudo R-squared
summary(IC_model, Nagelkerke = T)
summary(BIEN_model, Nagelkerke = T)
summary(TOTAL_model, Nagelkerke = T)
summary(BIEN_model_TOTAL, Nagelkerke = T)

# Standardize the variables to enable comparison between variable slopes
map01$TOTAL_std <- scale(map01$TOTAL)
map01$RELATIVE_std <- scale(map01$RELATIVE)
map01$RICHNESS_std <- scale(map01$RICHNESS)
map01$RANGE_MEDIAN_std <- scale(map01$RANGE_MEDIAN)
map01$ENDEMISM_R_std <- scale(map01$ENDEMISM_R)
map01$BIEN_OCCUR_std <- scale(map01$BIEN_OCCUR)
map01$GDP_CAPITA_std <- scale(map01$GDP_CAPITA)
map01$ROAD_DENSITY_std <- scale(map01$ROAD_DENSITY)
map01$POP_DENSITY_std <- scale(map01$POP_DENSITY)
map01$SECURITY_std <- scale(map01$SECURITY)
map01$RESEARCH_EXP_std <- scale(map01$RESEARCH_EXP)
map01$EDUCATION_EXP_std <- scale(map01$EDUCATION_EXP)


# Best models with standardized variables
IC_model_std <- errorsarlm(RELATIVE_std ~ RICHNESS_std + RANGE_MEDIAN_std + ENDEMISM_R_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
BIEN_model_std <- errorsarlm(BIEN_OCCUR_std ~ RICHNESS_std + RANGE_MEDIAN_std + GDP_CAPITA_std + POP_DENSITY_std + EDUCATION_EXP_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
TOTAL_model_std <- errorsarlm(TOTAL_std ~ RICHNESS_std + RANGE_MEDIAN_std + ENDEMISM_R_std + SECURITY_std + RESEARCH_EXP_std + EDUCATION_EXP_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)

# Second best models
IC_model_std2 <- errorsarlm(RELATIVE_std ~ RICHNESS_std + RANGE_MEDIAN_std + ENDEMISM_R_std + GDP_CAPITA_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
BIEN_model_std2 <- errorsarlm(BIEN_OCCUR_std ~ RICHNESS_std + RANGE_MEDIAN_std + EDUCATION_EXP_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
TOTAL_model_std2 <- errorsarlm(TOTAL_std ~ RICHNESS_std + RANGE_MEDIAN_std + ENDEMISM_R_std + SECURITY_std + RESEARCH_EXP_std + EDUCATION_EXP_std + POP_DENSITY_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)


# Third best models
IC_model_std3 <- errorsarlm(RELATIVE_std ~ RICHNESS_std + RANGE_MEDIAN_std + ENDEMISM_R_std + ROAD_DENSITY_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
BIEN_model_std3 <- errorsarlm(BIEN_OCCUR_std ~ RICHNESS_std + RANGE_MEDIAN_std + GDP_CAPITA_std + POP_DENSITY_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)
TOTAL_model_std3 <- errorsarlm(TOTAL_std ~ RICHNESS_std + RANGE_MEDIAN_std + ENDEMISM_R_std + SECURITY_std + RESEARCH_EXP_std + EDUCATION_EXP_std + POP_DENSITY_std + GDP_CAPITA_std, data = map01, listw = dn2000w, quiet = T, zero.policy = T, tol.solve = 1e-20)


# summaries of the best standardized models
summary(IC_model_std, Nagelkerke = T)
summary(BIEN_model_std, Nagelkerke = T)
summary(TOTAL_model_std, Nagelkerke = T)

# Summaries of second and third best models
summary(IC_model_std2, Nagelkerke = T)
summary(BIEN_model_std2, Nagelkerke = T)
summary(TOTAL_model_std2, Nagelkerke = T)
summary(IC_model_std3, Nagelkerke = T)
summary(BIEN_model_std3, Nagelkerke = T)
summary(TOTAL_model_std3, Nagelkerke = T)




# ---- Check if residuals are still spatially correlated ---- #



# Add residuals to dataset
map01$RESID_IC_SAR <- residuals(IC_model)
map01$RESID_BIEN_SAR <- residuals(BIEN_model_wNA)

# detect spatial autocorrelation with a Monto Carlo simulation of Moran's I
moran.mc(map01$RESID_IC_SAR, dn2000w, nsim = 999, zero.policy = T)
moran.mc(na.omit(map01$RESID_BIEN_SAR), dn2000w, nsim = 999, zero.policy = T)

# Plot the residuals
par(mfrow = c(1,2))
spplot(map01, "RESID_IC_SAR")
spplot(map01, "RESID_BIEN_SAR")

# Residuals are no longer spatially correlated