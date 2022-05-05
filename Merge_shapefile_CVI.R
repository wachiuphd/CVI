# Merge ToxPi Scores and Data into 2010 shapefile tracts
library(data.table)
library(sp)
library(raster)
library(rgdal) ## use install.packages(c("rgdal"),repos="https://mac.R-project.org")
p <- shapefile(file.path("Data","2010 Tracts","2010tracts.shp"))

# merge overall toxpi
shpdir <- file.path("Data","CVI ToxPi Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.toxpi.df <- fread("CVI-pct-comb.csv",integer64 = "double",
                      keepLeadingZeros = TRUE)
m <- merge(p, cvi.toxpi.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIToxPiTracts.shp"),overwrite=TRUE)
zip(paste0(shpdir,".zip"),shpdir) 
unlink(shpdir,recursive=TRUE) # clean up

# merge all pct data
shpdir <- file.path("Data","CVI Pct Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.pct.df <- fread("CVI_data_pct.csv",integer64 = "double",
                    keepLeadingZeros = TRUE)
m <- merge(p, cvi.pct.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIPctTracts.shp"),overwrite=TRUE)
zip(paste0(shpdir,".zip"),shpdir) 
unlink(shpdir,recursive=TRUE) # clean up

# merge all original data
shpdir <- file.path("Data","CVI Data Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.df<-fread("CVI_data_current.csv",integer64 = "double",
              keepLeadingZeros = TRUE)
m <- merge(p, cvi.df, by.x = "GEOID10", by.y = "GEOID.Tract")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIDataTracts.shp"),overwrite=TRUE)
zip(paste0(shpdir,".zip"),shpdir) 
unlink(shpdir,recursive=TRUE) # clean up

# merge each category
indicators.df<-fread("CVI_indicators_current.csv")
categories <- unique(indicators.df$`Baseline Vulnerability`)

for (i in 1:length(categories)) { 
  onecat <- categories[i]
  print(onecat)
  # Create directory for shapefile
  catdir <- file.path("Data",
                      paste0("CVI ToxPi.",onecat))
  if (!dir.exists(catdir)) dir.create(catdir)
  # read in one category ToxPi
  cvi.pct.toxpi.cat <- fread(paste0("CVI-pct-cat-",
                                    gsub(": ","-",onecat),".csv"),
                             integer64 = "double",
                             keepLeadingZeros = TRUE)
  m <- merge(p, cvi.pct.toxpi.cat, by.x = "GEOID10", by.y = "FIPS")
  # save as shapefile again
  shapefile(m, file.path(catdir,
                         paste0("CVIToxPi",onecat,".shp")),
                         overwrite=TRUE)
  # create zip file for uploading
  zip(paste0(catdir,".zip"),catdir) 
  unlink(catdir,recursive=TRUE) # clean up
}


