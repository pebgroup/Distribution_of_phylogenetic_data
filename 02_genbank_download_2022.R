# GenBank download editing script
# This script uses a subset of a download of the GenBank database
# It produces two lists of species that have data relevant to phylogeny available in GenBank:
# 1) a list of species that have at least ONE relevant sequence available
# 2) a list of species that preserves duplicate species if they appear with different molecular markers

#### TOTAL RUNTIME: ~ 3-4 HOURS ####

# Due to the long runtime, the results of the script are already uploaded to the Dryad repository


library(data.table)
library(dplyr)


# ---- Read in the subset of the GenBank download, which contains species name and description (subset was made in SQLite prior to local download to limit file size) ---- #


library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))


genbank.data <- fread("data\\gb_280321(4).txt", quote = "")
gdt <- data.frame(genbank.data, stringsAsFactors = F)
colnames(gdt) <- c("desc","species")




# ---- Create objects for each of the 128 most used markers in phylogeny according to Hinchliff & Smith (2014) ---- #




# Hinchliff & Smith (2004) describes which keywords that must be present in the GenBank description, and which keywords that cannot be present
accD <- gdt[grepl("accD", gdt$desc), ]
accD_psaI <- gdt[grepl("accD", gdt$desc) & grepl("psaI", gdt$desc) & !grepl("ycf4|rbcL", gdt$desc), ]
atp1 <- gdt[grepl("atp1", gdt$desc), ]
atp9 <- gdt[grepl("atp9|atp synthase subunit 9", gdt$desc) & !grepl("rps12", gdt$desc), ]
atpA <- gdt[grepl("atpA", gdt$desc), ]
atpA_mt <- gdt[grepl("atpA|atp1|atpase alpha subunit", gdt$desc) & grepl("mitochondria|mitochondrion|mitochondrial", gdt$desc) & !grepl("chloroplast|trnR|atp6|atpF|Arabidopsis|plasma membrane", gdt$desc), ]
atpB <- gdt[grepl("atpB", gdt$desc), ]
atpB_rbcL <- gdt[grepl("atpB", gdt$desc) & grepl("rbcL", gdt$desc) & !grepl("accd|atpe", gdt$desc), ]
atpE <- gdt[grepl("atpE", gdt$desc), ]
atpF <- gdt[grepl("atpF", gdt$desc), ]
atpH <- gdt[grepl("atpH", gdt$desc), ]
atpI <- gdt[grepl("atpI", gdt$desc), ]
atpI_atpH <- gdt[grepl("atpI", gdt$desc) & grepl("atpH", gdt$desc) & !grepl("rps2|atpF", gdt$desc), ]
ccsA <- gdt[grepl("ccsA", gdt$desc), ]
cemA <- gdt[grepl("cemA", gdt$desc), ]
clpP <- gdt[grepl("clpP", gdt$desc), ]
cob <- gdt[grepl("cob", gdt$desc) & !grepl("trnL|matK|rbcL|rpl16|ndhF|trnD|trnH|ribosulation|ribosylation|deoxyhypusine|rps14|cobi420|COBRA|Arabidopsis thaliana", gdt$desc), ]
cobi420_mosses <- gdt[grepl("cob", gdt$desc) & !grepl("trnL|matK|rbcL|rpl16|ndhF|trnD|trnH|ribosulation|deoxyhypusine|rps14", gdt$desc), ]
cox1 <- gdt[grepl("cox1", gdt$desc) & !grepl("cox17|cox2|coxii|cox3|coxiii cox19|cox1i624", gdt$desc), ]
cox1i624 <- gdt[grepl("cox1", gdt$desc), ]
cox2 <- gdt[grepl("cox2|coxii|cytochrome oxidase subunit 2|cytochrome oxidase subunit ii", gdt$desc) & !grepl("cox3|coxiii|nad3|cox1|coxi|psba|trnl|transcribed spacer|trnS|chloroplast", gdt$desc), ]
cox3 <- gdt[grepl("cox2|coxii|cytochrome oxidase subunit 2|cytochrome oxidase subunit ii", gdt$desc) & !grepl("cox2|coxii|nad3|cox1|coxi|psba|trnl|transcribed spacer|trnS|chloroplast", gdt$desc), ]
ETS <- gdt[grepl("external transcribed spacer|ets", gdt$desc), ]
infA <- gdt[grepl("infA", gdt$desc), ]
ITS <- gdt[grepl("its|its2|58S|internal transcribed spacer", gdt$desc) & !grepl("45S|5S|23S", gdt$desc), ]
ITS_23S_5S <- gdt[grepl("45S|5S|23S|internal transcribed spacer|its", gdt$desc) & !grepl("18S|58S|26S", gdt$desc), ]
matK <- gdt[grepl("matK", gdt$desc), ]
matR <- gdt[grepl("matR", gdt$desc), ]
nad2 <- gdt[grepl("nad2", gdt$desc), ]
nad5 <- gdt[grepl("nad5", gdt$desc), ]
nad5_nad4 <- gdt[grepl("nad5", gdt$desc) & grepl("nad4", gdt$desc), ]
ndhA <- gdt[grepl("ndhA", gdt$desc), ]
ndhB <- gdt[grepl("ndhB", gdt$desc), ]
ndhC <- gdt[grepl("ndhC", gdt$desc), ]
ndhD <- gdt[grepl("ndhD", gdt$desc), ]
ndhE <- gdt[grepl("ndhE", gdt$desc), ]
ndhF <- gdt[grepl("ndhF", gdt$desc), ]
ndhF_rpl32_trnL <- gdt[grepl("ndhF", gdt$desc) & grepl("rpl32", gdt$desc) & grepl("trnL", gdt$desc) & !grepl("ccsA|ycf1", gdt$desc), ]
ndhG <- gdt[grepl("ndhG", gdt$desc), ]
ndhH <- gdt[grepl("ndhH", gdt$desc), ]
ndhI <- gdt[grepl("ndhI", gdt$desc), ]
ndhJ <- gdt[grepl("ndhJ", gdt$desc), ]
ndhK <- gdt[grepl("ndhK", gdt$desc), ]
petA <- gdt[grepl("petA", gdt$desc), ]
petA_psbJ <- gdt[grepl("petA", gdt$desc) & grepl("psbJ", gdt$desc) & !grepl("cemA|psbL|psbF|psbE", gdt$desc), ]
petB <- gdt[grepl("petB", gdt$desc), ]
petD <- gdt[grepl("petD", gdt$desc), ]
petG <- gdt[grepl("petG", gdt$desc), ]
petL <- gdt[grepl("petL", gdt$desc), ]
petN <- gdt[grepl("petN", gdt$desc), ]
phyA <- gdt[grepl("phyA|phytochrome A", gdt$desc) & !grepl("Arabidposis thaliana|Pennisetum glaucum", gdt$desc), ]
phyB_monocots <- gdt[grepl("phyB|phytochrome B", gdt$desc) & !grepl("phyB1|phyB2", gdt$desc), ]
psaA <- gdt[grepl("psaA", gdt$desc), ]
psaB <- gdt[grepl("psaB", gdt$desc), ]
psaC <- gdt[grepl("psaC", gdt$desc), ]
psaI <- gdt[grepl("psaI", gdt$desc), ]
psaJ <- gdt[grepl("psaJ", gdt$desc), ]
psbA <- gdt[grepl("psbA", gdt$desc), ]
psbA_trnH <- gdt[grepl("trnH", gdt$desc) & !grepl("trnK|matK|rpl2|rpl32", gdt$desc), ]
psbA_trnK <- gdt[grepl("psbA", gdt$desc) & grepl("trnK", gdt$desc) & !grepl("rps16|trnH", gdt$desc), ]
psbB <- gdt[grepl("psbB", gdt$desc), ]
psbC <- gdt[grepl("psbC", gdt$desc), ]
psbD <- gdt[grepl("psbD", gdt$desc), ]
psbE <- gdt[grepl("psbE", gdt$desc), ]
psbE_petL <- gdt[grepl("psbE", gdt$desc) & grepl("petL", gdt$desc) & !grepl("psbF|psbL|petG", gdt$desc), ]
psbF <- gdt[grepl("psbF", gdt$desc), ]
psbH <- gdt[grepl("psbH", gdt$desc), ]
psbI <- gdt[grepl("psbI", gdt$desc), ]
psbJ <- gdt[grepl("psbJ", gdt$desc), ]
psbK <- gdt[grepl("psbK", gdt$desc), ]
psbL <- gdt[grepl("psbL", gdt$desc), ]
psbM <- gdt[grepl("psbM", gdt$desc), ]
psbM_trnD <- gdt[grepl("psbM", gdt$desc) & grepl("trnD", gdt$desc) & !grepl("trnY|trnE|petN|ycf6|trnC", gdt$desc), ]
psbN <- gdt[grepl("psbN", gdt$desc), ]
psbT <- gdt[grepl("psbT", gdt$desc), ]
psbT_psbN <- gdt[grepl("psbT", gdt$desc) & grepl("psbN", gdt$desc) & !grepl("psbB|psbH|Development of species_specific molecular markers for eight Salix species based on cpDNA polymorphism", gdt$desc), ]
psbZ <- gdt[grepl("psbZ", gdt$desc), ]
rbcL <- gdt[grepl("rbcL", gdt$desc), ]
rpl14 <- gdt[grepl("rpl14", gdt$desc), ]
rpL14_rpS8_rpL36 <- gdt[grepl("rpL14", gdt$desc) & grepl("rpS8", gdt$desc) & grepl("rpL36", gdt$desc) & !grepl("rps3", gdt$desc), ]
rpl16 <- gdt[grepl("rpl16", gdt$desc), ]
rpl2 <- gdt[grepl("rpl2", gdt$desc), ]
rpl20 <- gdt[grepl("rpl20", gdt$desc), ]
rpl22 <- gdt[grepl("rpl22", gdt$desc), ]
rpl23 <- gdt[grepl("rpl23", gdt$desc), ]
rpl32 <- gdt[grepl("rpl32", gdt$desc), ]
rpl33 <- gdt[grepl("rpl33", gdt$desc), ]
rpl36 <- gdt[grepl("rpl36", gdt$desc), ]
rpoA <- gdt[grepl("rpoA", gdt$desc), ]
rpoB <- gdt[grepl("rpoB", gdt$desc), ]
rpoC1 <- gdt[grepl("rpoC1", gdt$desc), ]
rpoC2 <- gdt[grepl("rpoC2", gdt$desc), ]
rps11 <- gdt[grepl("rps11", gdt$desc), ]
rps12 <- gdt[grepl("rps12", gdt$desc), ]
rps12_rpl20 <- gdt[grepl("rps12", gdt$desc) & grepl("rpl20", gdt$desc) & !grepl("clpP|rps18", gdt$desc), ]
rps14 <- gdt[grepl("rps14", gdt$desc), ]
rps15 <- gdt[grepl("rps15", gdt$desc), ]
rps16 <- gdt[grepl("rps16", gdt$desc), ]
rps16_trnK <- gdt[grepl("rps16", gdt$desc) & grepl("trnK", gdt$desc) & !grepl("trnQ|matK|Veronica", gdt$desc), ]
rps18 <- gdt[grepl("rps18", gdt$desc), ]
rps19 <- gdt[grepl("rps19", gdt$desc), ]
rps2 <- gdt[grepl("rps2", gdt$desc), ]
rps3 <- gdt[grepl("rps3", gdt$desc), ]
rps4 <- gdt[grepl("rps4", gdt$desc), ]
rps4_trnT_trnL <- gdt[grepl("rpS4", gdt$desc) & grepl("trnT", gdt$desc) & grepl("trnL", gdt$desc) & !grepl("trnS|trnF", gdt$desc), ]
rps7 <- gdt[grepl("rps7", gdt$desc), ]
rps8 <- gdt[grepl("rps8", gdt$desc), ]
rrn16 <- gdt[grepl("rrn16", gdt$desc), ]
rrn23 <- gdt[grepl("rrn23", gdt$desc), ]
rrn4.5 <- gdt[grepl("rrn4.5", gdt$desc), ]
rrn5 <- gdt[grepl("rrn5", gdt$desc), ]
trnC_rpoB <- gdt[grepl("trnC", gdt$desc) & grepl("rpoB", gdt$desc) & !grepl("petN|ycf6|rpoC1", gdt$desc), ]
trnC_ycf6_psbM <- gdt[grepl("ycf6|petN", gdt$desc) & grepl("psbM", gdt$desc) & !grepl("rpoB|trnD", gdt$desc), ]
trnF_ndhJ <- gdt[grepl("trnF", gdt$desc) & grepl("ndhJ", gdt$desc) & !grepl("trnL|ndhK", gdt$desc), ]
trnG_intron_bryos <- gdt[grepl("trnG", gdt$desc) & !grepl("trnR|trnS", gdt$desc), ]
trnG_intron_tracheophyta <- gdt[grepl("trnG", gdt$desc) & !grepl("trnR|trnS", gdt$desc), ]
trnK_intron <- gdt[grepl("trnK", gdt$desc) & !grepl("rps16|psbA", gdt$desc), ]
trnQ_rps16 <- gdt[grepl("trnQ", gdt$desc) & grepl("rps16", gdt$desc) & !grepl("trnK|psbK", gdt$desc), ]
trnS_UGA_trnfM <- gdt[grepl("trnS|tRNA_Ser", gdt$desc) & grepl("trnfM|trnM", gdt$desc) & !grepl("rps4|psbI|psbK|psbC|trnG|trnQ|ycf3|trnS_GCU|trnR|psbD", gdt$desc), ]
trnS_rpS4 <- gdt[grepl("trnS", gdt$desc) & grepl("rps4", gdt$desc) & !grepl("trnT|ycf3", gdt$desc), ]
trnT_psbD <- gdt[grepl("trnT", gdt$desc) & grepl("psbD", gdt$desc) & !grepl("psbC|trnE|trnY", gdt$desc), ]
trnT_trnL_trnF <- gdt[grepl("trnT|trnL|trnF", gdt$desc) & !grepl("trnE|nad7|rps4|ndhJ|ndhB|psbM|psbA|trnS|trnfM|trnD", gdt$desc), ]
trnV_ndhC <- gdt[grepl("trnV", gdt$desc) & grepl("ndhC", gdt$desc), ]
trnV_trnM <- gdt[grepl("trnV", gdt$desc) & grepl("trnM", gdt$desc), ]
trnY_trnE <- gdt[grepl("trnY", gdt$desc) & grepl("trnE", gdt$desc) & !grepl("trnT", gdt$desc), ]
ycf1 <- gdt[grepl("ycf1", gdt$desc), ]
ycf2 <- gdt[grepl("ycf2", gdt$desc), ]
ycf3 <- gdt[grepl("ycf3", gdt$desc), ]



# 1)
# ---- Produce a table that contains the names of all species with a GenBank entry relevant to phylogeny ---- #




# Move all species from the above produced objects into a single list
list_of_species <- list(accD = accD$species, 
                        accD_psaI = accD_psaI$species, 
                        atp1 = atp1$species, 
                        atp9 = atp9$species, 
                        atpA = atpA$species, 
                        atpA_mt = atpA_mt$species, 
                        atpB = atpB$species, 
                        atpB_rbcL = atpB_rbcL$species, 
                        atpE = atpE$species, 
                        atpF = atpF$species, 
                        atpH = atpH$species, 
                        atpI = atpI$species, 
                        atpI_atpH = atpI_atpH$species, 
                        ccsA = ccsA$species, 
                        cemA = cemA$species, 
                        clpP = clpP$species, 
                        cob = cob$species, 
                        cobi420_mosses = cobi420_mosses$species, 
                        cox1 = cox1$species, 
                        cox1i624 = cox1i624$species, 
                        cox2 = cox2$species, 
                        cox3 = cox3$species, 
                        ETS = ETS$species, 
                        infA = infA$species, 
                        ITS = ITS$species, 
                        ITS_23S_5S = ITS_23S_5S$species, 
                        matK = matK$species, 
                        matR = matR$species, 
                        nad2 = nad2$species, 
                        nad5 = nad5$species, 
                        nad5_nad4 = nad5_nad4$species, 
                        ndhA = ndhA$species, 
                        ndhB = ndhB$species, 
                        ndhC = ndhC$species, 
                        ndhD = ndhD$species, 
                        ndhE = ndhE$species, 
                        ndhF = ndhF$species, 
                        ndhF_rpl32_trnL = ndhF_rpl32_trnL$species, 
                        ndhG = ndhG$species, 
                        ndhH = ndhH$species, 
                        ndhI = ndhI$species, 
                        ndhJ = ndhJ$species, 
                        ndhK = ndhK$species, 
                        petA = petA$species, 
                        petA_psbJ = petA_psbJ$species, 
                        petB = petB$species, 
                        petD = petD$species, 
                        petG = petG$species, 
                        petL = petL$species, 
                        petN = petN$species, 
                        phyA = phyA$species, 
                        phyB_monocots = phyB_monocots$species, 
                        psaA = psaA$species, 
                        psaB = psaB$species, 
                        psaC = psaC$species, 
                        psaI = psaI$species, 
                        psaJ = psaJ$species, 
                        psbA = psbA$species, 
                        psbA_trnH = psbA_trnH$species, 
                        psbA_trnK = psbA_trnK$species, 
                        psbB = psbB$species, 
                        psbC = psbC$species, 
                        psbD = psbD$species, 
                        psbE = psbE$species, 
                        psbE_petL = psbE_petL$species, 
                        psbF = psbF$species, 
                        psbH = psbH$species, 
                        psbI = psbI$species, 
                        psbJ = psbJ$species, 
                        psbK = psbK$species, 
                        psbL = psbL$species, 
                        psbM = psbM$species, 
                        psbM_trnD = psbM_trnD$species, 
                        psbN = psbN$species, 
                        psbT = psbT$species, 
                        psbT_psbN = psbT_psbN$species, 
                        psbZ = psbZ$species, 
                        rbcL = rbcL$species, 
                        rpl14 = rpl14$species, 
                        rpL14_rpS8_rpL36 = rpL14_rpS8_rpL36$species, 
                        rpl16 = rpl16$species, 
                        rpl2 = rpl2$species, 
                        rpl20 = rpl20$species, 
                        rpl22 = rpl22$species, 
                        rpl23 = rpl23$species, 
                        rpl32 = rpl32$species, 
                        rpl33 = rpl33$species, 
                        rpl36 = rpl36$species, 
                        rpoA = rpoA$species, 
                        rpoB = rpoB$species, 
                        rpoC1 = rpoC1$species, 
                        rpoC2 = rpoC2$species, 
                        rps11 = rps11$species, 
                        rps12 = rps12$species, 
                        rps12_rpl20 = rps12_rpl20$species, 
                        rps14 = rps14$species, 
                        rps15 = rps15$species, 
                        rps16 = rps16$species, 
                        rps16_trnK = rps16_trnK$species, 
                        rps18 = rps18$species, 
                        rps19 = rps19$species, 
                        rps2 = rps2$species, 
                        rps3 = rps3$species, 
                        rps4 = rps4$species, 
                        rps4_trnT_trnL = rps4_trnT_trnL$species, 
                        rps7 = rps7$species, 
                        rps8 = rps8$species, 
                        rrn16 = rrn16$species, 
                        rrn23 = rrn23$species, 
                        rrn4.5 = rrn4.5$species, 
                        rrn5 = rrn5$species, 
                        trnC_rpoB = trnC_rpoB$species, 
                        trnC_ycf6_psbM = trnC_ycf6_psbM$species, 
                        trnF_ndhJ = trnF_ndhJ$species, 
                        trnG_intron_bryos = trnG_intron_bryos$species, 
                        trnG_intron_tracheophyta = trnG_intron_tracheophyta$species, 
                        trnK_intron = trnK_intron$species, 
                        trnQ_rps16 = trnQ_rps16$species, 
                        trnS_UGA_trnfM = trnS_UGA_trnfM$species, 
                        trnS_rpS4 = trnS_rpS4$species, 
                        trnT_psbD = trnT_psbD$species, 
                        trnT_trnL_trnF = trnT_trnL_trnF$species, 
                        trnV_ndhC = trnV_ndhC$species, 
                        trnV_trnM = trnV_trnM$species, 
                        trnY_trnE = trnY_trnE$species, 
                        ycf1 = ycf1$species, 
                        ycf2 = ycf2$species, 
                        ycf3 = ycf3$species
)

# Convert into dataframe
list_as_df <- data.frame(unlist(list_of_species))
list_as_df2 <- list_as_df$unlist.list_of_species.
df1 <- data.frame(list_as_df2)



# ---- Edit the list to limit length and better match the WCSP names ---- #




# Remove duplicates and add a header
unique.sp <- unique(df1)
unique.sp <- as.data.frame(unique.sp)
colnames(unique.sp) <- c("species")

# Remove entries called "Hybrid cultivar", not specified species (sp), and uncertain entries (aff and cf)
unique.sp <- unique.sp[!grepl(" hybrid | sp | aff | cf ", unique.sp$species),]

# Edit names to match WCSP, e.g. "Sommera donnell_smithii" to "Sommera donnell-smithii"
unique.sp <- gsub("_", "-", unique.sp) 

# include "." to match WCSP names
unique.sp <- gsub(" subsp ", " subsp. ", unique.sp)
unique.sp <- gsub(" var ", " var. ", unique.sp)
unique.sp <- gsub(" f ", " f. ", unique.sp)

unique.sp <- as.data.frame(unique.sp)
colnames(unique.sp) <- c("species")




# ---- Translate any synonyms included in the Genbank Download to their accepted name using a list of synonyms and accepted names produced from WCSP data.




# Read in the synonym translator and list of accepted WCVP names produced in the wcvp_subset_2021.R script
synonym.translator <- fread("synonym_translator_2022.csv", quote = ",", header = T)
names.id <- fread("ID_and_Names_2022.csv", sep = ",", header = T)

# Translate the synonyms using the synonym translator table and remove new duplicates 
indx <- match(unlist(unique.sp), synonym.translator$taxon_name, nomatch = 0)
unique.sp$species[indx !=0] <- synonym.translator$accepted_plant_name_id[indx]
translated.gb.list <- unique(unique.sp$species)




# ---- Remove any subspecies extension from names. The binomial name is preserved, as the subspecies represents an entry for the species, yet it will not be recognized as a match with WCVP data, if it includes an extension ---- #




fin.gb.list <- gsub(" subsp. .*$", "", translated.gb.list)
fin.gb.list <- gsub(" var. .*$", "", fin.gb.list)
fin.gb.list <- gsub(" f. .*$", "", fin.gb.list)

fin.gb.list <- as.data.frame(fin.gb.list)
colnames(fin.gb.list) <- c("species")




# ---- Write the finished table of species with phylogenetically relevant GenBank entries (1) ---- #




write.table(fin.gb.list, "genbank_entries_2022.csv", sep = ",", row.names = F, col.names = T, quote = F)



# 2)
# ---- Additional analysis of GenBank data, where all relevant sequences are preserved. This is to account for the depth of phylogenetic knowledge for each species ---- #




# Create data frame with two colomns: species and marker(desc)
accD$desc <- "accD"
accD_psaI$desc <- "accD_psaI"
atp1$desc <- "atp1"
atp9$desc <- "atp9"
atpA$desc <- "atpA"
atpA_mt$desc <- "atpA_mt"
atpB$desc <- "atpB"
atpB_rbcL$desc <- "atpB_rbcL"
atpE$desc <- "atpE"
atpF$desc <- "atpF"
atpH$desc <- "atpH"
atpI$desc <- "atpI"
atpI_atpH$desc <- "atpI_atpH"
ccsA$desc <- "ccsA"
cemA$desc <- "cemA"
clpP$desc <- "clpP"
cob$desc <- "cob"
cobi420_mosses$desc <- "cobi420_mosses"
cox1$desc <- "cox1"
cox1i624$desc <- "cox1i624"
cox2$desc <- "cox2"
cox3$desc <- "cox3"
ETS$desc <- "ETS"
infA$desc <- "infA"
ITS$desc <- "ITS"
ITS_23S_5S$desc <- "ITS_23S_5S"
matK$desc <- "matK"
matR$desc <- "matR"
nad2$desc <- "nad2"
nad5$desc <- "nad5"
nad5_nad4$desc <- "nad5_nad4"
ndhA$desc <- "ndhA"
ndhB$desc <- "ndhB"
ndhC$desc <- "ndhC"
ndhD$desc <- "ndhD"
ndhE$desc <- "ndhE"
ndhF$desc <- "ndhF"
ndhF_rpl32_trnL$desc <- "ndhF_rpl32_trnL"
ndhG$desc <- "ndhG"
ndhH$desc <- "ndhH"
ndhI$desc <- "ndhI"
ndhJ$desc <- "ndhJ"
ndhK$desc <- "ndhK"
petA$desc <- "petA"
petA_psbJ$desc <- "petA_psbJ"
petB$desc <- "petB"
petD$desc <- "petD"
petG$desc <- "petG"
petL$desc <- "petL"
petN$desc <- "petN"
phyA$desc <- "phyA"
phyB_monocots$desc <- "phyB_monocots"
psaA$desc <- "psaA"
psaB$desc <- "psaB"
psaC$desc <- "psaC"
psaI$desc <- "psaI"
psaJ$desc <- "psaJ"
psbA$desc <- "psbA"
psbA_trnH$desc <- "psbA_trnH"
psbA_trnK$desc <- "psbA_trnK"
psbB$desc <- "psbB"
psbC$desc <- "psbC"
psbD$desc <- "psbD"
psbE$desc <- "psbE"
psbE_petL$desc <- "psbE_petL"
psbF$desc <- "psbF"
psbH$desc <- "psbH"
psbI$desc <- "psbI"
psbJ$desc <- "psbJ"
psbK$desc <- "psbK"
psbL$desc <- "psbL"
psbM$desc <- "psbM"
psbM_trnD$desc <- "psbM_trnD"
psbN$desc <- "psbN"
psbT$desc <- "psbT"
psbT_psbN$desc <- "psbT_psbN"
psbZ$desc <- "psbZ"
rbcL$desc <- "rbcL"
rpl14$desc <- "rpl14"
rpL14_rpS8_rpL36$desc <- "rpL14_rpS8_rpL36"
rpl16$desc <- "rpl16"
rpl2$desc <- "rpl2"
rpl20$desc <- "rpl20"
rpl22$desc <- "rpl22"
rpl23$desc <- "rpl23"
rpl32$desc <- "rpl32"
rpl33$desc <- "rpl33"
rpl36$desc <- "rpl36"
rpoA$desc <- "rpoA"
rpoB$desc <- "rpoB"
rpoC1$desc <- "rpoC1"
rpoC2$desc <- "rpoC2"
rps11$desc <- "rps11"
rps12$desc <- "rps12"
rps12_rpl20$desc <- "rps12_rpl20"
rps14$desc <- "rps14"
rps15$desc <- "rps15"
rps16$desc <- "rps16"
rps16_trnK$desc <- "rps16_trnK"
rps18$desc <- "rps18"
rps19$desc <- "rps19"
rps2$desc <- "rps2"
rps3$desc <- "rps3"
rps4$desc <- "rps4"
rps4_trnT_trnL$desc <- "rps4_trnT_trnL"
rps7$desc <- "rps7"
rps8$desc <- "rps8"
rrn16$desc <- "rrn16"
rrn23$desc <- "rrn23"
rrn4.5$desc <- "rrn4.5"
rrn5$desc <- "rrn5"
trnC_rpoB$desc <- "trnC_rpoB"
trnC_ycf6_psbM$desc <- "trnC_ycf6_psbM"
trnF_ndhJ$desc <- "trnF_ndhJ"
trnG_intron_bryos$desc <- "trnG_intron_bryos"
trnG_intron_tracheophyta$desc <- "trnG_intron_tracheophyta"
trnK_intron$desc <- "trnK_intron"
trnQ_rps16$desc <- "trnQ_rps16"
trnS_UGA_trnfM$desc <- "trnS_UGA_trnfM"
trnS_rpS4$desc <- "trnS_rpS4"
trnT_psbD$desc <- "trnT_psbD"
trnT_trnL_trnF$desc <- "trnT_trnL_trnF"
trnV_ndhC$desc <- "trnV_ndhC"
trnV_trnM$desc <- "trnV_trnM"
trnY_trnE$desc <- "trnY_trnE"
ycf1$desc <- "ycf1"
ycf2$desc <- "ycf2"
ycf3$desc <- "ycf3"


spec.mark <- rbind(accD,
                   accD_psaI,
                   atp1,
                   atp9,
                   atpA,
                   atpA_mt,
                   atpB,
                   atpB_rbcL,
                   atpE,
                   atpF,
                   atpH,
                   atpI,
                   atpI_atpH,
                   ccsA,
                   cemA,
                   clpP,
                   cob,
                   cobi420_mosses,
                   cox1,
                   cox1i624,
                   cox2,
                   cox3,
                   ETS,
                   infA,
                   ITS,
                   ITS_23S_5S,
                   matK,
                   matR,
                   nad2,
                   nad5,
                   nad5_nad4,
                   ndhA,
                   ndhB,
                   ndhC,
                   ndhD,
                   ndhE,
                   ndhF,
                   ndhF_rpl32_trnL,
                   ndhG,
                   ndhH,
                   ndhI,
                   ndhJ,
                   ndhK,
                   petA,
                   petA_psbJ,
                   petB,
                   petD,
                   petG,
                   petL,
                   petN,
                   phyA,
                   phyB_monocots,
                   psaA,
                   psaB,
                   psaC,
                   psaI,
                   psaJ,
                   psbA,
                   psbA_trnH,
                   psbA_trnK,
                   psbB,
                   psbC,
                   psbD,
                   psbE,
                   psbE_petL,
                   psbF,
                   psbH,
                   psbI,
                   psbJ,
                   psbK,
                   psbL,
                   psbM,
                   psbM_trnD,
                   psbN,
                   psbT,
                   psbT_psbN,
                   psbZ,
                   rbcL,
                   rpl14,
                   rpL14_rpS8_rpL36,
                   rpl16,
                   rpl2,
                   rpl20,
                   rpl22,
                   rpl23,
                   rpl32,
                   rpl33,
                   rpl36,
                   rpoA,
                   rpoB,
                   rpoC1,
                   rpoC2,
                   rps11,
                   rps12,
                   rps12_rpl20,
                   rps14,
                   rps15,
                   rps16,
                   rps16_trnK,
                   rps18,
                   rps19,
                   rps2,
                   rps3,
                   rps4,
                   rps4_trnT_trnL,
                   rps7,
                   rps8,
                   rrn16,
                   rrn23,
                   rrn4.5,
                   rrn5,
                   trnC_rpoB,
                   trnC_ycf6_psbM,
                   trnF_ndhJ,
                   trnG_intron_bryos,
                   trnG_intron_tracheophyta,
                   trnK_intron,
                   trnQ_rps16,
                   trnS_UGA_trnfM,
                   trnS_rpS4,
                   trnT_psbD,
                   trnT_trnL_trnF,
                   trnV_ndhC,
                   trnV_trnM,
                   trnY_trnE,
                   ycf1,
                   ycf2,
                   ycf3
)

spec.mark <- data.frame(spec.mark, stringsAsFactors = F)

# Repeat the cleaning process and synchronize with WSVP names, WITHOUT removing duplicate entries yet 

# Remove entries called "Hybrid cultivar", not specified species (sp), and uncertain entries (aff and cf)
spec.mark <- spec.mark[!grepl(" hybrid | sp | aff | cf ", spec.mark$species),]

# Edit names to match WCSP, e.g. "Sommera donnell_smithii" to "Sommera donnell-smithii"
spec.mark$species <- gsub("_", "-", spec.mark$species)

# include "." to match WCSP names
spec.mark$species <- gsub(" subsp ", " subsp. ", spec.mark$species)
spec.mark$species <- gsub(" var ", " var. ", spec.mark$species)
spec.mark$species <- gsub(" f ", " f. ", spec.mark$species)

# Translate the synonyms using the synonym translator table and remove new duplicates 
synonym.translator <- fread("data\\synonym_translator_2022.csv", quote = ",", header = T)
names.id <- fread("data\\ID_and_Names_2022.csv", sep = ",", header = T)
index <- match(spec.mark$species, synonym.translator$taxon_name, nomatch = 0)
spec.mark$species[index !=0] <- synonym.translator$accepted_plant_name_id[index]

# Remove any subspecies extension from names. The binomial name is preserved, as the subspecies represents an entry for the species, yet it will not be recognized as a match with WCSP data, if it includes an extension
spec.mark$species <- gsub(" subsp. .*$", "", spec.mark$species)
spec.mark$species <- gsub(" var. .*$", "", spec.mark$species)
spec.mark$species <- gsub(" f. .*$", "", spec.mark$species)

#save(spec.mark, file = "C:\\Users\\alexa\\Documents\\SPECIALE\\publication\\220921\\results_2021\\specmark.RData")


# ---- Remove duplicate species entries within each marker to preserve maximum 1 entry per species per marker.  


# To avoid memory exhaustion, divide the df into smaller sizes and merge them after the process
spec.mark1 <- unique(spec.mark[1:3000000,c("desc" , "species")])
spec.mark2 <- unique(spec.mark[3000001:6000000,c("desc" , "species")])
spec.mark3 <- unique(spec.mark[6000001:9000000,c("desc" , "species")])
spec.mark4 <- unique(spec.mark[9000001:12000000,c("desc" , "species")])
spec.mark5 <- unique(spec.mark[12000001:(length(spec.mark$species)),c("desc" , "species")])

# Merge the divided df's and run unique function again to account for duplicates that has been preserved due to the division. This can now be done on the entire df, because of the smaller size.
spec.mark.uni <- rbind(spec.mark1, spec.mark2, spec.mark3, spec.mark4, spec.mark5)
spec.mark.unique <- unique(spec.mark.uni[c("desc" , "species")])


# Remove irrelevant columns by excluding entries that will not match the WSVP to reduce computation time an later analysis
species.all.markers <- as.data.frame(spec.mark.unique$species[spec.mark.unique$species %in% intersect(spec.mark.unique$species, names.id$taxon_name)], stringsAsFactors = F)
colnames(species.all.markers) <- c("species")

# Write table of genbank entries that preserve duplicate species entries from different molecular markers (2)
write.table(species.all.markers, "data\\genbank_entries_w_duplicates_2022.csv", sep = ",", row.names = F, col.names = T, quote = F)

