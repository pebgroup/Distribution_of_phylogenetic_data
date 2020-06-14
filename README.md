# Mapping and explaining our phylogenetic knowledge
The respository includes scripts used to map phylogenetic knowledge of plants based on the World Checklist of Selected Plant Families (WCSP) and GenBank. It also includes a statistical scripts that creates a model to determine predictors of phylogenetic knowledge, as well as a .csv with all response- and explanatory variables.



## R scripts:

**data_manipulation.R:** The input of the script is published and unpublished WCSP names and distributions. The script manipulates this data to produce three files: 
* (1) a list of all WCSP names and their correspond WCSP ID, 
* (2) a table with columns for all botanical countries and the species present in thoses countries 
* (3) A list of all WCSP synonyms and their corresponding accepted names


**genbank_download.R:** The input to the script is a subset of a GenBank download, which contains only species name and description. It produces a list of species that have data relevant to phylogenetics available in GenBank


**data_analysis_01.R:** The input to the script are the products of "data_manipulation.R" and "genbank_download.R", as decribed above. It produces a table, which includes values for eight variables for all 368 botanical countries.


**statistics.R:** The input to the script is the table "all_variables_2303.csv", and the shapefile "level3.shp". It produces a statistical analysis of which variables can predict phylogenetic knowledge, phylogenetic effort and distribution knowledge with a spatial error model.



## Metadata:

**all_variables_2303.csv:** A table based on the product of "data_analysis_01.R", with the additional inclusion of socioeconomic variables adjusted to botanical countries in ArcMap (ESRI, 2011). This table is the input to "statistics.R"


**level3.shp:** a shapefile of botanical countries (L3 regions)
