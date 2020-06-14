# Biodiversity variable creation script
# This script uses datasets produced in the "Data manipulation" script, including:
# 1) ID's and Names
# 2) WCSP distribution data in wide format
# 3) The edited GenBank download,
# As well as a list of species with BIEN occurrences for each botanical country, translated into WCSP ID's (produced by Melanie Tietje),
# it produces multiple variables related to biodiversity and two response variables at the scale of botanical contries, which is gathered in one table. These include:
# 1) Total number of species with relevent molecular data on GenBank
# 2) Inventory completeness of flora with relevent molecular data on GenBank (response)
# 3) Inventory completeness of flora with a record occurrence in BIEN (response)
# 4) Species richness
# 5) Mean number of countries a species of a country is present in
# 6) Median number of countries a species of a country is present in
# 7) Total number of endemic species 
# 8) Endemics relative to species richness




# ---- Read data from Genbank, Checklist, names and id table and the BIEN .Rdata file ---- #




gen.data <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Genetic_data\\GB_230320\\genbank_entries_2303.csv", header = T) # As produced in the genbank_downlad.R script
wcsp.data <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\230320\\wcsp_wide_2303.csv", header = T) # As produced in the data_manipulation.R script
name.id <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\ID_and_Names.csv", header = T) # As produced in the data_manipulation.R script
load("~/SPECIALE/DATA/BIEN/BIEN_in_WCSP_regions.RData") # R workspace with objects "spec.list" and "res". Spec.list is a list of all L3 regions with the species recorded there in BIEN, translated to WCSP ID's. Res is a df with the lengths of the elements of spec.list.  




# ---- Prepare vectors and dataframes for the loops ---- #




# Create vector of area codes and taxa with genetic data
area.vector <- c(colnames(wcsp.data))
genbank.vector <- c(as.character(gen.data$species))

# Create dataframe with data for distribution range
dist.range <- table(unlist(wcsp.data))
dist.range.df <- data.frame(dist.range)

# Convert spec.list into df, as the loop can't handle it as a list
spec.list.equal.length <- lapply(spec.list, 'length<-', max(lengths(spec.list)))
spec.df <- as.data.frame(spec.list.equal.length)

# Create vectors to store results
result.list.total <- c()
result.list.relative <- c()
result.list.diversity <- c()
result.list.distrange.mean <- c()
result.list.distrange.median <- c()
result.list.endemism.total <- c()
result.list.endemism.relative <- c()
bien_results <- c()




# ---- Fill vectors with data for desired variables ---- #




# Number of species in each region with genetic data available (Phylogenetic effort)
for (i in area.vector){
  vector01 <- c(as.character(wcsp.data[,i]))
  result.list.total[[i]] <- length(intersect(genbank.vector, vector01))
}

# Number of species with genetic data relative to total number of species in each region (Phylogenetic knowledge, main response variable)
for (i2 in area.vector){
  vector02 <- c(as.character(wcsp.data[,i2]))
  result.list.relative[[i2]] <- length(intersect(genbank.vector, vector02))/length(vector02[!is.na(vector02)])
}

# Total species richness
for (i3 in area.vector){
  vector03 <- c(as.character(wcsp.data[,i3]))
  result.list.diversity[[i3]] <- length(vector03[!is.na(vector03)])
}

# Mean number of regions in which a species occur for each region (mean range)
for (i4 in area.vector){
  vector04 <- c(na.omit(as.character(wcsp.data[,i4])))
  result.list.distrange.mean[[i4]] <- mean(dist.range.df$Freq[match(vector04, dist.range.df$Var1)])
}

# Median number of regions in which a species occur for each region (median range)
for (i5 in area.vector){
  vector05 <- c(na.omit(as.character(wcsp.data[,i5])))
  result.list.distrange.median[[i5]] <- median(dist.range.df$Freq[match(vector05, dist.range.df$Var1)])
}

# Number of species that only occur in the specific region (total endemics)
for (i6 in area.vector){
  vector06 <- c(na.omit(as.character(wcsp.data[,i6])))
  result.list.endemism.total[[i6]] <- sum((dist.range.df$Freq[match(vector06, dist.range.df$Var1)]) == 1)
}

# number of endemics relative to total number of species in the region (endemics relative to SR)
for (i7 in area.vector){
  vector07 <- c(na.omit(as.character(wcsp.data[,i7])))
  result.list.endemism.relative[[i7]] <- (sum((dist.range.df$Freq[match(vector07, dist.range.df$Var1)]) == 1))/length(dist.range.df$Freq[match(vector07, dist.range.df$Var1)])
}

# Number of species with available BIEN data in a region relative to SR (Distribution Knowledge)
for (i in area.vector)
{
  wcsp.acc.name <- name.id$accepted_species[match(spec.df[,i], name.id$alt_id)]
  wcsp.total <- wcsp.data[,i]
  bien_results[[i]] <- length(intersect(wcsp.acc.name, wcsp.total))/length(wcsp.total[!is.na(wcsp.total)])
}

# if desired, print results to check the results
#print(result.list.total)
#print(result.list.relative)
#print(result.list.diversity)
#print(result.list.distrange.mean)
#print(result.list.distrange.median)
#print(result.list.endemism.total)
#print(result.list.endemism.relative)
#print(bien_results)




# ---- Combine results into single dataframe ---- #




# convert vectors to df's
total.df <- data.frame(result.list.total)
relative.df <- data.frame(result.list.relative)
diversity.df <- data.frame(result.list.diversity)
dist.range.mean.df <- data.frame(result.list.distrange.mean)
dist.range.median.df <- data.frame(result.list.distrange.median)
endemism.total.df <- data.frame(result.list.endemism.total)
endemism.relative.df <- data.frame(result.list.endemism.relative)
bien.results.df <- data.frame(bien_results)

# combine df's into one df
col.headers <- c("TOTAL")
names(total.df) <- col.headers
total.df$RELATIVE <- relative.df$result.list.relative
total.df$RICHNESS <- diversity.df$result.list.diversity
total.df$RANGE_MEAN <- dist.range.mean.df$result.list.distrange.mean
total.df$RANGE_MEDIAN <- dist.range.median.df$result.list.distrange.median
total.df$ENDEMISM_T <- endemism.total.df$result.list.endemism.total
total.df$ENDEMISM_R <- endemism.relative.df$result.list.endemism.relative
total.df$BIEN_OCCUR <- bien.results.df$bien_results
total.df <- cbind(L3 = rownames(total.df), total.df)
total.df$RANGE_MEAN <- round(total.df$DISTRANGEMEAN, digits = 2)
total.df$RANGE_MEDIAN <- round(total.df$DISTRANGEMEDIAN, digits = 2)
total.df$ENDEMISM_R <- round(total.df$ENDEMISMREL, digits = 3)




# ---- Write results to be mapped and used in the statistical analysis ---- #
write.table(total.df,"C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\results\\VAR_230320\\variables_2303.csv",sep = ",", col.names = T, row.names = F, quote = F)

