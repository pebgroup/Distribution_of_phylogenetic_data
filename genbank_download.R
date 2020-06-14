# GenBank download editing script
# This script uses a subset of a download of the GenBank database
# It produces a list of species that have data relevant to phylogeny available in GenBank




# ---- Read in the subset of the GenBank download, which contains species name and description (subset was made in SQLite prior to local download to limit file size) ---- #




genbank.data <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Genetic_data\\sqlite_download\\ascsv\\desc_data2.csv", sep = ",", header = T)
gdt <- data.frame(genbank.data, stringsAsFactors = F)




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




# Read in the translator produced in the data_manipulation.R script
synonym.translator <- read.csv("C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Checklist\\synonyms\\synonym_translator_2303.csv", sep = ",", header = T)

# Translate the synonyms using the synonym translator table and remove new duplicates 
indx <- match(unlist(unique.sp), synonym.translator$synonym, nomatch = 0)
unique.sp$species[indx !=0] <- synonym.translator$accepted_name[indx]
translated.gb.list <- unique(unique.sp$species)




# ---- Remove any subspecies extension from names. The binomial name is preserved, as the subspecies represents an entry for the species, yet it will not be recognized as a match with WCSP data, if it includes an extension ---- #




fin.gb.list <- gsub(" subsp. .*$", "", translated.gb.list)
fin.gb.list <- gsub(" var. .*$", "", fin.gb.list)
fin.gb.list <- gsub(" f. .*$", "", fin.gb.list)

fin.gb.list <- as.data.frame(fin.gb.list)
colnames(fin.gb.list) <- c("species")




# ---- Write the finished table of species with phylogenetically relevant GenBank entries ---- #



write.table(fin.gb.list, "C:\\Users\\alexa\\Documents\\SPECIALE\\DATA\\Genetic_data\\GB_230320\\genbank_entries_2303.csv", sep = ",", row.names = F, col.names = T, quote = F)

