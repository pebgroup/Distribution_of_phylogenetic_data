# The Darwinian shortfall in plants: phylogenetic knowledge is driven by range size
The respository includes scripts used to map and analyze phylogenetic knowledge of plants based on the World Checklist of Vascular Plants (WCVP) and GenBank. 



## R scripts:

The data that is required to run the scripts is available for download at https://doi.org/10.5061/dryad.2547d7wrz. Unzip the data in a folder named 'data' within the folder containing these scripts. Run scripts in numerical order. 

**NOTE: "02_genbank_download.R" has a runtime that exceeds 3 hours. The script can be skipped, however, as the results are included in the Dryad repository.**

**01_wcvp_subset_2022:** The input of the script is WCVP datasets with names and distributions. The script produces four files: 
* 1) a table of all WCVP ID's and the corresponding name (names and id's)
* 2) a table with two columns: accepted species name and distribution
* 3) a table with columns for all botanical countries (L3 regions), and the species present in thoses countries
* 4) a table with all WCVP synonyms and the corresponding accepted name (synonym translator)

**02_genbank_download.R_2022:** The input to the script is a subset of a GenBank download, which contains only species name and description. It produces two lists of species that have data relevant to phylogeny available in GenBank:
* 1) a list of species that have at least ONE relevant sequence available
* 2) a list of species that preserves duplicate species if they appear with different molecular markers
* NOTE: LONG RUNTIME

**03_phylo_know_analysis_2022.R:** This script uses three datasets produced in the "wcvp_subset_2021" script, two subsets of the GenBank download produced in "genbank_download.R, and BIEN data:  
WCVP:
* ID's and Names
* WCVP distribution data in wide format with countries as columns
* WCVP distribution data in a long format with two columns: species and botanical country  
GenBank:
* a list of species that have at least ONE relevant sequence available
* a list of species that preserves duplicate species if they appear with different molecular markers  
BIEN:
* a list of species with BIEN occurrences for each botanical country, translated into WCSP ID's (BIEN_in_WCSP_regions_Sept21_2020download.RData)  
Botanical countries:
* a table of Level 3 botanical counrty codes and their area size  

The script produces multiple explanatory variables related to biodiversity and four response variables at the scale of botanical contries, which is all gathered in one table. These variables include:
* Total number of species with relevent molecular data on GenBank (response)
* Inventory completeness of flora with at least one entry with relevant sequence data on GenBank (response)
* Inventory completeness of flora including all relevant molecular data on GenBank (response)
* Inventory completeness of flora with a record occurrence in BIEN (response)
* Species richness
* Mean number of countries a species of a country is present in
* Median number of countries a species of a country is present in
* Mean range size of species
* Total number of endemic species 
* Endemics relative to species richness
* Proporption of endemics sequenced

**04_statistical_modelling_2022.R:** The input to the script is: 
* the dataset of variables produced in "phylo_know_analysis_2021.R"
* socioeconomic variables fitted to botanical countries
* shapefile of botanical countries

It conducts a statistical analysis of the dataset through spatial autoregressive modelling.

**05_logit_2022.RMD**
* the script uses "seq_area_2022.csv" as produced in "03_phylo_know_analysis_2021.R"
* Logistic regression of the sequencing status of species (1 = sequenced, 0 = not sequenced) on their range size, measured as the sum of the areas (in km^2) of all the botanical countries a species occurs in. 

**bien_TDWG_intersection_Sep2021_SERVER.R**
* script to run on GenomeDK due to memomry requirements (min 16GB, 2h runtime). Assigns BIEN occurrences to botanical countries
