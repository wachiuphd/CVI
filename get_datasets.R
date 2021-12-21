library(data.table)
library(openxlsx)
datafolder <- "Data"

### Get census tracts from 2019 Census Tract Gazetteer File
### https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2019_Gazetteer/2019_Gaz_tracts_national.zip
tracts <- fread(file.path(datafolder,"2019_Gaz_tracts_national.txt"),
                keepLeadingZeros = TRUE)
tracts[, c("ALAND","AWATER","ALAND_SQMI","AWATER_SQMI","INTPTLAT","INTPTLONG"):=NULL]
tracts$GEOID.State <- substring(tracts$GEOID,1,2)
tracts$GEOID.County <- substring(tracts$GEOID,1,5)

### Get master sheet

cvi.master <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/CurrentCVIIndicatorsDoc - 12.20.2021.xlsx",
                        sheet="Subcategories (In Progress)")
for (j in 1:41) {
  fname<-(file.path("~/Dropbox/Climate Health Vulnerability Index",
                            cvi.master$`Path To File`[j]))
  fname.exists <- file.exists(fname)
  cat(fname," exists? ",fname.exists,"\n")
  if (fname.exists) {
    tmp <- fread(fname,
                 keepLeadingZeros = TRUE)
    cat(names(tmp),sep="\n")
    cat("\n")
    # cols<-c("GEOID",cvi.master$`Indicator Name`[j])
    # setnames(tmp, cvi.master$`GEOID Column Name`[j], cols[1])
    # setnames(tmp, cvi.master$`Data Column Name`[j], cols[2])
    # tmp <- tmp[,..cols]
  }
}
