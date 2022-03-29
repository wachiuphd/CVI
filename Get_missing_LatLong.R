library(data.table)
library(stringr)
library(dplyr)
library(naniar)
library(readxl)
datafolder <- "Data"

# Census tracts from 22-02_CVI_state_county_tract.xlsx
tractsraw <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/Other/22-02_CVI_state_county_tract.xlsx",
                       sheet="Tract")
tracts <- tractsraw[,c("STATE","County_Name")]
tracts$GEOID.State <- tractsraw$STATEFP10
tracts$GEOID.County <- tractsraw$FIPS
tracts$GEOID.Tract <- tractsraw$GEOID10
tracts <- tracts[order(tracts$GEOID.Tract),]

tracts_latlong2014 <- fread(file.path(datafolder,"2014_Gaz_tracts_national.txt"),
                            keepLeadingZeros = TRUE)
tracts_latlong2014[, c("USPS","ALAND","AWATER","ALAND_SQMI","AWATER_SQMI"):=NULL]
setnames(tracts_latlong2014,"GEOID","GEOID.Tract")
tracts_latlong2014$LatLong <- paste0(tracts_latlong2014$INTPTLAT,",",tracts_latlong2014$INTPTLONG)
tracts_latlong2014[, c("INTPTLAT","INTPTLONG"):=NULL]

# Join lat long to tracts
tracts2014 <- left_join(tracts,tracts_latlong2014)

toxpiboundaries2014 <- fread(file.path(datafolder,"ToxPiBoundaries.csv"),keepLeadingZeros = TRUE)
toxpiboundaries2014$GEOID.Tract <- toxpiboundaries2014$GEOID10
tracts2014_boundaries <- left_join(tracts2014,toxpiboundaries2014)

tracts2014_latlong <- tracts2014_boundaries[,c("STATE","County_Name","GEOID.State","GEOID.County","GEOID.Tract","INTPTLAT10","INTPTLON10")]

fwrite(tracts2014_latlong,file.path(datafolder,"Tracts_LatLong.csv"),quote=TRUE)
fwrite(subset(tracts2014_latlong,is.na(INTPTLAT10)),file.path(datafolder,"Tracts_LatLong_missing.csv"),quote=TRUE)

