

#### GEOGRAPHY ##########
library(data.table)
library(rgdal)

dat <- readRDS("data/BIEN_WCSP_merged_Sep21_no_centroids.rds")
dat <- dat[-which(is.na(dat$accepted_plant_name_id)),]
# 
# # Cleaning
# ## coordinates
# dat <- dat[!abs(dat$lat)>90,]
# dat <- dat[!abs(dat$lng)>180,]


# sort data
dat <- setorder(dat, accepted_plant_name_id)

# Read botanical countries shapefile
shp <- readOGR("data/shapefile/level3.shp")
proj4string(shp)

# loop over regions
coord <- dat[,c("lat", "lng")]
coordinates(coord) <- ~  lng + lat  # convert dataframe to spatial points object
proj4string(coord) <- proj4string(shp) # match projection attributes of both objects

# see rgdal's transition to PROJ6.
# throws an error when trying to assign projections attributes


# BIEN species counts & list per WCSP region #####################
# record all occurring species and calculate the intersection percentage with the WCSP species in that region

regions <- shp$LEVEL_3_CO
res <- data.frame(region=rep(NA,length(regions)), sr=rep(0, length(regions)))
spec.list <- list()
for(i in 1:length(regions)){
  temp <- shp[i,]
  # get species that occur in this region
  region_occurrences <-  over(temp, geometry(coord), returnList = TRUE)
  region_species <- dat$accepted_plant_name_id[as.numeric(unlist(region_occurrences))]
  # accumulate species occurrences to presence data
  res[i, 1] <- as.character(temp$LEVEL_3_CO)
  res[i, 2] <- length(unique(region_species))
  # write species names into the list for each region
  #names(spec.list[i]) <- as.character(temp$LEVEL_3_CO)
  spec.list[[i]] <- unique(region_species)
  if(!i%%1)cat(i,"\r")
}
names(spec.list) <- shp$LEVEL_3_CO


save(res, spec.list, file="BIEN_in_WCSP_regions_Sept21.RData")

#if(getwd()=="/data_vol/melanie"){
#  save(res, spec.list, file="BIEN_in_WCSP_regions.RData")
#}else{
#  save(res, spec.list, file="BIEN_in_WCSP_regions_local_subset.RData")
#}


# add notification if run is finished



################### END SERVER OUTSOURCING #############################

