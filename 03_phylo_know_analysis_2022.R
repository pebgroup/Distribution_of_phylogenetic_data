# Biodiversity variable creation script
#
# This script uses three datasets produced in the "wcvp_subset_2021" script, two subsets of the GenBank download produced in "genbank_download.R, and BIEN data:
# WCVP:
# 1) ID's and Names
# 2) WCVP distribution data in wide format with countries as columns
# 3) WCVP distribution data in a long format with two columns: species and botanical country
# GenBank:
# 4) a list of species that have at least ONE relevant sequence available
# 5) a list of species that preserves duplicate species if they appear with different molecular markers
# BIEN:
# 6) a list of species with BIEN occurrences for each botanical country, translated into WCSP ID's (BIEN_in_WCSP_regions_Sept21_2020download.RData)
# Botanical countries:
# 7) a table of Level 3 botanical counrty codes and their area size
#
# This script produces multiple explanatory variables related to biodiversity and four response variables at the scale of botanical contries, which is all gathered in one table. These variables include:
# 1) Total number of species with relevent molecular data on GenBank (response)
# 2) Inventory completeness of flora with at least one entry with relevant sequence data on GenBank (response)
# 3) Inventory completeness of flora including all relevant molecular data on GenBank (response)
# 4) Inventory completeness of flora with a record occurrence in BIEN (response)
# 5) Species richness
# 6) Mean number of countries a species of a country is present in
# 7) Median number of countries a species of a country is present in
# 8) Mean range size of species
# 9) Total number of endemic species 
# 10) Endemics relative to species richness
# 11) Proporption of endemics sequenced

# Total runtime: ~5 minutes


# ---- Read data from Genbank, WCVP and BIEN ---- #

library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))


library(data.table)
gen.data <- fread("data\\genbank_entries_2022.csv", header = T, quote = "", sep = NULL) # As produced in the genbank_download.R script
wcsp.data <- fread("data\\wcsp_wide_2022.csv", header = T) # As produced in the wcvp_subset_2021.R script
name.id <- fread("data\\ID_and_Names_2022.csv", header = T) # As produced in the wcvp_subset_2021.R script
load("BIEN_in_WCSP_regions_sep2021.RData") # R workspace with objects "spec.list" and "res". Spec.list is a list of all L3 regions with the species recorded there in BIEN, translated to WCSP ID's. Res is a df with the lengths of the elements of spec.list.
wcsp.long <- fread("data\\wcsp_long_2022.csv", quote = "") # As produced in the wcvp_subset_2021.R script
L3.area <- fread("data\\L3_and_area.csv", header = T) # Botanical countries and their area size
gen.data.duplicate <- fread("data\\genbank_entries_w_duplicates_2022.csv", quote = "", header = T, sep = NULL) # As produced in the genbank_download.R script




# ---- Prepare vectors and dataframes for the variable creation loops ---- #




# Create vector of area codes and taxa with genetic data
area.vector <- c(colnames(wcsp.data))
genbank.vector <- c(as.character(gen.data$species))
genbank.vector.duplicate <- c(as.character(gen.data.duplicate$species))

# Create dataframe with data for distribution range
dist.range <- table(unlist(wcsp.data))
dist.range.df <- data.frame(dist.range)
dist.range.df.endemics <- dist.range.df[dist.range.df$Freq == 1,]

# Create dataframe of species and area, aggregate and sum it, and create a wide-style dataframe with sqkm replacing species.
wcsp.long$area_sq <- L3.area$AREA_SQKM[match(wcsp.long$area_code_l3, L3.area$LEVEL_3_CO)]
L3.area.only <- subset(wcsp.long, select = c("accepted_plant_name" , "area_sq"))
L3.area.sum <- aggregate(. ~ accepted_plant_name, L3.area.only, sum)
L3.area.sum <- as.data.frame(L3.area.sum)
wide.sq <- wcsp.data
wide.sq[] <- lapply(wcsp.data, function(x) L3.area.sum$area_sq[match(x, L3.area.sum$accepted_plant_name)])

# Convert spec.list and wcsp.data into df, as the loop can't handle it as a list
spec.list.equal.length <- lapply(spec.list, 'length<-', max(lengths(spec.list)))
spec.df <- as.data.frame(spec.list.equal.length)
wcsp.data <- as.data.frame(wcsp.data)
wide.sq <- as.data.frame(wide.sq)


# Create vectors to store results
result.list.total <- c()
result.list.relative <- c()
result.list.diversity <- c()
result.list.distrange.mean <- c()
result.list.distrange.median <- c()
result.list.range.area <- c()
result.list.endemism.total <- c()
result.list.endemism.relative <- c()
result.list.bien <- c()
result.list.endemics.list <- c()
result.list.endemics.seq <- c()
result.list.gen.dup <- c()




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

# Total area of the regions in which a species occur (range area)
for (i8 in area.vector){
  vector08 <- c(na.omit(as.numeric(wide.sq[,i8])))
  result.list.range.area[[i8]] <- mean(vector08)
}

# Number of species that only occur in the specific region (total endemics)
for (i6 in area.vector){
  vector06 <- c(na.omit(as.character(wcsp.data[,i6])))
  result.list.endemism.total[[i6]] <- sum((dist.range.df$Freq[match(vector06, dist.range.df$Var1)]) == 1)
}

# Number of endemics relative to total number of species in the region (endemics relative to SR)
for (i7 in area.vector){
  vector07 <- c(na.omit(as.character(wcsp.data[,i7])))
  result.list.endemism.relative[[i7]] <- (sum((dist.range.df$Freq[match(vector07, dist.range.df$Var1)]) == 1))/length(dist.range.df$Freq[match(vector07, dist.range.df$Var1)])
}


# Number of species with available BIEN data in a region relative to SR (Distribution Knowledge)
for (i in area.vector)
{
  wcsp.acc.name <- name.id$taxon_name[match(spec.df[,i], name.id$plant_name_id)]
  wcsp.total <- wcsp.data[,i]
  result.list.bien[[i]] <- length(intersect(wcsp.acc.name, wcsp.total))/length(wcsp.total[!is.na(wcsp.total)])
}


# List of all species that are endemic to a single botanical country
for (i9 in area.vector){
  vector09 <- c(na.omit(as.character(wcsp.data[,i9])))
  result.list.endemics.list[[i9]] <- dist.range.df.endemics$Var1[match(vector09, dist.range.df.endemics$Var1)]
}

# Format the list created above to become fit for looping
endemism.equal.length <- lapply(result.list.endemics.list, 'length<-', max(lengths(result.list.endemics.list)))
endemism.as.df <- as.data.frame(endemism.equal.length)

# % of endemic species sequenced
for (i10 in area.vector){
  vector10 <- c(as.character(endemism.as.df[,i10]))
  result.list.endemics.seq[[i10]] <- length(intersect(genbank.vector, vector10))/length(vector10[!is.na(vector10)])
}


# Inventory completeness of phylogenetic knowledge, when completeness is achieved by having available sequences for all 128 relevant markers for each species
for (i11 in area.vector){
  vector11 <- c(as.character(wcsp.data[,i11]))
  result.list.gen.dup[[i11]] <- length(genbank.vector.duplicate[genbank.vector.duplicate %in% intersect(genbank.vector.duplicate, vector11)])/((length(vector11[!is.na(vector11)]))*128)
}


# Aggregate results into table

# convert vectors to df's
total.df <- data.frame(TOTAL=unlist(result.list.total),
                       RELATIVE=unlist(result.list.relative),
                       RICHNESS=unlist(result.list.diversity),
                       RANGE_MEAN=unlist(result.list.distrange.mean),
                       RANGE_MEDIAN=unlist(result.list.distrange.median),
                       RANGE_AREA=unlist(result.list.range.area),
                       ENDEMISM_T=unlist(result.list.endemism.total),
                       ENDEMISM_R=unlist(result.list.endemism.relative),
                       BIEN_OCCUR=unlist(result.list.bien),
                       ENDEMIC_SEQ=unlist(result.list.endemics.seq),
                       GEN_DUP=unlist(result.list.gen.dup)
                       )
total.df <- data.frame(LEVEL_3_CO = row.names(total.df), total.df)



# ---- Create table to be used for logistic regression of the sequencing status of species (1 = sequenced, 0 = not sequenced) on their range size, measured as the sum of the areas (in km^2) of all the botanical countries a species occurs in (logit.Rmd) 



L3.area.sum$sequenced <- ifelse(L3.area.sum$accepted_plant_name %in% gen.data$species, 1, 0)



# ---- Write results to be mapped and used in the statistical analysis ---- #




write.table(total.df,"data\\variables_2022.csv",sep = ",", col.names = T, row.names = F, quote = F)

write.table(endemics.seq.df,"data\\endemic_completeness_2022.csv",sep = ",", col.names = F, row.names = T, quote = F)

write.table(L3.area.sum,"data\\seq_area_2022.csv", sep = ",", col.names = T, row.names = F, quote = F)



# ---- Additional information of interest included in the paper ---- #


# % of species with phylogentically relevant data

length(intersect(unique(wcsp.long$accepted_plant_name) , unique(gen.data$species))) # total number

length(intersect(unique(wcsp.long$accepted_plant_name) , unique(gen.data$species)))/length(unique(wcsp.long$accepted_plant_name)) # %
