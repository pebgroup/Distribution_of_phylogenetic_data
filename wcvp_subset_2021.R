# WCVP data script
# This script uses the WCVP datasets for plant names and their distribution
# It produduces 4 subset datasets:
# 1) a table of all WCVP ID's and the corresponding name (names and id's)
# 2) a table with two columns: accepted species name and distribution
# 3) a table with columns for all botanical countries (L3 regions), and the species present in thoses countries
# 4) a table with all WCVP synonyms and the corresponding accepted name (synonym translator)

library(data.table)
library(rgdal)
library(tidyr)

# ---- Read in the WCVP data and subset for relevent entries ---- #




wcsp.names <- fread("C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\WCSP2021\\checklist_names.txt", quote = "") # WCVP names file
wcsp.distribution <- fread("C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\WCSP2021\\checklist_distribution.txt", quote = "") # WCVP distribution file

# Exclude fern and moss families
ferns <- fread("C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\WCSP\\fern_list.txt", quote = "", header = F, col.names = "family") # list of ferns
moss <- fread("C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\WCSP\\bryophyta.csv", quote = "", header = F, col.names = "family") # list of mosses
wcsp.names <- subset(wcsp.names, !family %in% ferns$family)
wcsp.names <- subset(wcsp.names, !family %in% moss$family)

# Include only accepted entries
wcsp.names.species <- wcsp.names[grepl("Accepted", wcsp.names$taxon_status),]

# Synonymize plant name IDs in distribution data to use accepted plant name IDs, and subset to taxa with accepted plant name IDs

# Merge with accepted names:
wcsp.distribution <- merge(wcsp.distribution, wcsp.names.species[,c("accepted_plant_name_id", "plant_name_id", "taxon_rank")],
                           by="plant_name_id", all.x=TRUE)

# Remove all entries with no accepted plant name ID
wcsp.dist.subset <- wcsp.distribution[!is.na(wcsp.distribution$accepted_plant_name_id),]

# Remove subsp., var., f. etc that was introduced to the data after synonymizing from distribution data
wcsp.dist.subset <- wcsp.dist.subset[wcsp.dist.subset$taxon_rank == "Species", ]

# Remove extinct and introduced species
wcsp.dist.subset <- subset(wcsp.dist.subset, wcsp.dist.subset$extinct == 0)
wcsp.dist.subset <- subset(wcsp.dist.subset, wcsp.dist.subset$introduced == 0)

# Remove or correct mistaken Area codes
shape <- readOGR("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\shapefile\\level3.shp")
wcsp.dist.subset$area_code_l3 = toupper(wcsp.dist.subset$area_code_l3)
wcsp.dist.subset <- wcsp.dist.subset[wcsp.dist.subset$area_code_l3 %in% shape$LEVEL_3_CO,]

# produce table of WCSP ID's and corresponding names
name.id <- wcsp.names[grepl("Accepted", wcsp.names$taxon_status),]
name.id <- name.id[!(name.id$species == ""), ]
name.id <- subset(name.id, select = c(plant_name_id, taxon_name))

# Replace WCSP ID's with accepted names using the ID and Names object
wcsp.dist.subset$accepted_plant_name <- name.id$taxon_name[match(wcsp.dist.subset$accepted_plant_name_id, name.id$plant_name_id)]
wcsp.dist.subset <- na.omit(wcsp.dist.subset)

# Replace "Ã-" with "x" to indicate accepted hybrids, as this is the format used in GenBank
wcsp.dist.subset <- as.data.frame(sapply(wcsp.dist.subset, gsub, pattern = "Ã-", replacement = "x"))




# ---- Format the WCSP subset into a dataframe that includes L3 areas as columns ---- #




# Create vector of unique L3 area codes to use in the loop
l3.vector <- c(as.character(unique(wcsp.dist.subset$area_code_l3)))

# Create an empty list to fill
area.list <- list()

# Loop and search for all species associated with an l3 region and arrange them in a list
for (i in l3.vector){
  inventory <- wcsp.dist.subset$accepted_plant_name[wcsp.dist.subset$area_code_l3 == i]
  area.list[[i]] <- inventory
}

# Convert the list into the final dataframe format, which will be used in the "diversity_variables" script, where the response variable and biodiversity-related explanatory variables will be produced.
list.equal.length <- lapply(area.list, 'length<-', max(lengths(area.list)))
wcsp.df <- as.data.frame(list.equal.length) # The table is written in the end of the script




# ---- Create a table of synonyms and their accepted names to be used to translate synonyms from Genbank ---- #




# Retrive records for all not-accepted names and remove entries at genus level
wcsp.synonyms <- wcsp.names[!grepl("Accepted", wcsp.names$taxon_status),]
wcsp.synonyms <- wcsp.synonyms[!(wcsp.synonyms$species == ""), ]

# Combine columns: Genus, Genus hybrid marker, Species hybrid marker, Species, Infraspecific rank and infraspecific epithet to producea new column with the full synonym
wcsp.synonyms <- replace(wcsp.synonyms, wcsp.synonyms == "", NA)

# Remove irrelvant columns
wcsp.synonyms <- subset(wcsp.synonyms, select = c(taxon_name,accepted_plant_name_id))

# Substitute plant ID with plant name from "Names and Id's"

wcsp.synonyms$accepted_plant_name_id <- name.id$taxon_name[match(wcsp.synonyms$accepted_plant_name_id, name.id$plant_name_id)]
wcsp.synonyms <- na.omit(wcsp.synonyms)

# Remove cases where synonym = accepted name to reduce computation time during translation
# Note: These cases occur when the only difference between synonym and accepted name is the "author" columns. This is never included in the name on GenBank, so it is safe to remove these)
wcsp.synonyms <- wcsp.synonyms[wcsp.synonyms$taxon_name != wcsp.synonyms$accepted_plant_name_id,]

# Remove duplicated synonyms
# Note: Again the "author" columns is the culprit. These always refer to the same accepted name, so one is expandable to reduce computation time.
wcsp.synonyms <- wcsp.synonyms[!duplicated(wcsp.synonyms$taxon_name),]

# Remove cases where a synonym is the accepted name of another species. 
# Note: This one is complicated. Several species names exists as both synonym and accepted name, but for different species. This is, again, due to differences in the author extension to the synonym, which is not included here, as it does not occur on Genbank. -
# - Including these cases would result in the possible translation of an accepted name into another accepted name, which is undesirable
wcsp.synonyms <- wcsp.synonyms[!(wcsp.synonyms$taxon_name %in% intersect(wcsp.synonyms$accepted_plant_name_id, wcsp.synonyms$taxon_name)),]

# Translate "Ã-" to "x" to adapt to the GenBank download
wcsp.synonyms <- as.data.frame(sapply(wcsp.synonyms, gsub, pattern = "Ã-", replacement = "x"))

# Remove cases where the accepted name is "Unplaced Unplaced"
wcsp.synonyms <- wcsp.synonyms[wcsp.synonyms$accepted_plant_name_id != "Unplaced Unplaced",]


# The final synonym translator file can be written in the end of the script




# ---- Write tables for further use ---- #




# Write a table that includes all accepted species and their ID's as metadata (ID's and Names)
write.table(name.id,"C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\results_2021\\ID_and_Names_220921.csv", sep = ",", row.names = F, quote = F)

# Write a table with the relevant subset of the WCSP distribution data as useful metadata (WCSP long format)
write.table(wcsp.dist.subset,"C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\results_2021\\wcsp_long_220921.csv", sep = ",", row.names = F, col.names = T, quote = F)

# Write a table with the relevant subset of the WCSP distribution data in the format necessary for further analysis (WCSP wide format)
write.table(wcsp.df,"C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\results_2021\\wcsp_wide_220921.csv", sep = ",", row.names = F, col.names = T, quote = F)

# Write a table with synonyms and their accepted names to be used to translate the GenBank download into only accepted names (Synonym translator)
write.table(wcsp.synonyms,"C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\results_2021\\synonym_translator_220921.csv", sep = ",", row.names = F, col.names = T, quote = F)




# ---- Additional information of interest included in the paper ---- #



# Number of species with distribution data
length(unique(wcsp.dist.subset$accepted_plant_name))

# number of species in total
wcsp.total <- wcsp.names.species[wcsp.names.species$taxon_rank == "Species", ]
length(unique(wcsp.total$accepted_plant_name_id))

# Percentage of species with distribution data
length(unique(wcsp.dist.subset$accepted_plant_name))/length(unique(wcsp.total$accepted_plant_name_id))


