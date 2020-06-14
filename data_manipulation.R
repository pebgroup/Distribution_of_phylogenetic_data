# Initial data manipulation script
# This script uses WCSP data for plant names and distribution
# It produduces:
# 1) a table of all WCSP ID's and the corresponding name
# 2) a table with columns for all botanical countries (L3 regions), and the species present in thoses countries
# 3) a table with all WCSP synonyms and the corresponding accepted name




# ---- Read in the WCSP data in and combine them into two large dataset ---- #




published.names <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\published_names_19_10_2018.csv", sep = "|", stringsAsFactors = F)
unpublished.names <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\unpublished_names_19_10_2018.csv", sep = "|", stringsAsFactors = F)
wcsp.names <- rbind(published.names, unpublished.names)

published.distribution <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\published_distribution_19_10_2018.csv", sep = "|", stringsAsFactors = F)
unpublished.distribution <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\unpublished_distribution_19_10_2018.csv", sep = "|", stringsAsFactors = F)
wcsp.distribution <- rbind(published.distribution, unpublished.distribution)




# ---- Produce a table that includes all accepted species and their WCSP ID, which can be used in other analyses as well ---- #




# Subset the WCSP name dataset for only accepted species, and removing entries that are only at genus level
wcsp.names.subset <- wcsp.names[grepl("Accepted", wcsp.names$taxon_status_description),]
wcsp.names.subset <- wcsp.names.subset[!(wcsp.names.subset$species == ""), ]
name.id <- subset(wcsp.names.subset, select = -c(checklist,ipni_id,family,genus,genus_hybrid_marker,species_hybrid_marker,species,infraspecific_rank,infraspecific_epithet,primary_author,parenthetical_author,taxon_status_description,taxon_name,lifeform,climate,redlist_category,accepted_name_id,accepted_family) )

# A useful table of ID's and names can be written in the end of the script




# ---- Subset the WCSP distribution data by removing entries that are synonyms, extinct, introduced or sub- species ---- #




# Remove synonyms by matching the WCSP distribution data to the ID and Names object created above, which only contains accepted species
wcsp.dist.subset <- wcsp.distribution[(wcsp.distribution$checklist_id %in% name.id$checklist_id),]

# Remove extinct and introduced species
wcsp.dist.subset <- subset(wcsp.dist.subset, wcsp.dist.subset$extinct == 0)
wcsp.dist.subset <- subset(wcsp.dist.subset, wcsp.dist.subset$introduced == 0)

# Remove all irrelevant columns from the subset and entries with no distribution data
wcsp.dist.subset <- subset (wcsp.dist.subset, select = c("checklist_id","area_code_l3"))
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == ""), ]

# Replace WCSP ID's with accepted names using the ID and Names object
wcsp.dist.subset$checklist_id <- name.id$accepted_name[match(wcsp.dist.subset$checklist_id, name.id$checklist_id)]
wcsp.dist.subset <- na.omit(wcsp.dist.subset)

# Remove subspecies and missing data
wcsp.dist.subset <- wcsp.dist.subset[!grepl("subsp.", wcsp.dist.subset$checklist_id),]
wcsp.dist.subset <- wcsp.dist.subset[!grepl("var.", wcsp.dist.subset$checklist_id),]
wcsp.dist.subset <- wcsp.dist.subset[!grepl("f.", wcsp.dist.subset$checklist_id),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$checklist_id == ""), ]

# Remove entries that had misspelled area codes
# Note: these were detected after the first run, as they would show up as an area with each only one species. It is not always obvious which letter is mistyped, and they were therefor removed, instead of corrected
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "NHA"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "SRI"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "ALL"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "BKF"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "BUT"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "CCP"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "MIS"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "QLS"),]
wcsp.dist.subset <- wcsp.dist.subset[!(wcsp.dist.subset$area_code_l3 == "YRM"),]

# Convert area codes that have been mistankenly entered in lowercase to uppercase
wcsp.dist.subset[,2] <- toupper(wcsp.dist.subset[,2])

# Replace "Ã-" with "x" to indicate accepted hybrids, as this is the format used in GenBank
wcsp.dist.subset <- as.data.frame(sapply(wcsp.dist.subset, gsub, pattern = "Ã-", replacement = "x"))




# ---- Format the subset from having to columns with "Accepted name" and "Area code" to a df where unique area codes (368) are columns containg their present species. 




# Create vector of unique L3 area codes to use in the loop
l3.vector <- c(as.character(unique(wcsp.dist.subset$area_code_l3)))

# Create an empty list to fill
area.list <- list()

# Loop and search for all species associated with an l3 region and arrange them in a list
for (i in l3.vector){
  inventory <- wcsp.dist.subset$checklist_id[wcsp.dist.subset$area_code_l3 == i]
  area.list[[i]] <- inventory
}

# Convert the list into the final dataframe format, which will be used in the "diversity_variables" script, where the response variable and biodiversity-related explanatory variables will be produced.
list.equal.length <- lapply(area.list, 'length<-', max(lengths(area.list)))
wcsp.df <- as.data.frame(list.equal.length) # The table is written in the end of the script




# ---- Create a table of synonyms and their accepted names to be used to translate synonyms from Genbank ---- #




# Retrive records for all not-accepted names and remove entries at genus level
wcsp.synonyms <- wcsp.names[!grepl("Accepted", wcsp.names$taxon_status_description),]
wcsp.synonyms <- wcsp.synonyms[!(wcsp.synonyms$species == ""), ]

# Combine columns: Genus, Genus hybrid marker, Species hybrid marker, Species, Infraspecific rank and infraspecific epithet to producea new column with the full synonym
library(tidyr)
wcsp.synonyms <- replace(wcsp.synonyms, wcsp.synonyms == "", NA)
wcsp.synonyms <- unite(wcsp.synonyms, synonym, genus:infraspecific_epithet, sep = " ", na.rm = TRUE)

# Remove irrelvant columns
wcsp.synonyms <- subset(wcsp.synonyms, select = c(synonym,accepted_name))

# Remove cases where synonym = accepted name to reduce computation time during translation
# Note: These cases occur when the only difference between synonym and accepted name is the "author" columns. This is never included in the name on GenBank, so it is safe to remove these)
wcsp.synonyms <- wcsp.synonyms[wcsp.synonyms$synonym != wcsp.synonyms$accepted_name,]

# Remove duplicated synonyms
# Note: Again the "author" columns is the culprit. These always refer to the same accepted name, so one is expandable to reduce computation time.
wcsp.synonyms <- wcsp.synonyms[!duplicated(wcsp.synonyms$synonym),]

# Remove cases where a synonym is the accepted name of another species. 
# Note: This one is complicated. Several species names exists as both synonym and accepted name, but for different species. This is, again, due to differences in the author extension to the synonym, which is not included here, as it does not occur on Genbank. -
# - Including these cases would result in the possible translation of an accepted name into another accepted name, which is undesirable
wcsp.synonyms <- wcsp.synonyms[!(wcsp.synonyms$synonym %in% intersect(wcsp.synonyms$accepted_name, wcsp.synonyms$synonym)),]

# Translate "Ã-" to "x" to adapt to the GenBank download
wcsp.synonyms <- as.data.frame(sapply(wcsp.synonyms, gsub, pattern = "Ã-", replacement = "x"))

# Remove cases where the accepted name is "Unplaced Unplaced"
wcsp.synonyms <- wcsp.synonyms[wcsp.synonyms$accepted_name != "Unplaced Unplaced",]

# The finalized translator can be written in the end of the script




# ---- Write the desired tables for further use ("ID's and Names" and "WCSP wide format" are necessary for variable creation) ---- #




# Write a table that includes all accepted species and their ID's as metadata (ID's and Names)
#write.table(name.id,"C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\ID_and_Names.csv", sep = ",", row.names = F, quote = F)

# Write a table with the relevant subset of the WCSP distribution data as useful metadata (WCSP long format)
write.table(wcsp.dist.subset,"C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\230320\\wcsp_long_2303.csv", sep = ",", row.names = F, col.names = T, quote = F)

# Write a table with the relevant subset of the WCSP distribution data in the format necessary for further analysis (WCSP wide format)
write.table(wcsp.df,"C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\230320\\wcsp_wide_2303.csv", sep = ",", row.names = F, col.names = T, quote = F)

# Write a table with synonyms and their accepted names to be used to translate the GenBank download into only accepted names (Synonym translator)
#write.table(wcsp.synonyms,"C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\synonyms\\synonym_translator_2303.csv", sep = ",", row.names = F, col.names = F, quote = F)
