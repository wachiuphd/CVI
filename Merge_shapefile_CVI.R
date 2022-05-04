# Merge ToxPi Scores and Data into 2010 shapefile tracts
library(sp)
library(raster)
library(rgdal) ## use install.packages(c("rgdal"),repos="https://mac.R-project.org")
p <- shapefile(file.path("Data","2010 Tracts","2010tracts.shp"))

# # Load CSV files
cvi.toxpi.df <- fread("CVI-pct-comb.csv",integer64 = "double",
                      keepLeadingZeros = TRUE)
cvi.df<-fread("CVI_data_current.csv",integer64 = "double",
              keepLeadingZeros = TRUE)

# merge on common variable
m <- merge(p, cvi.toxpi.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path("Data","CVI ToxPi Tracts","CVIToxPiTracts.shp"),overwrite=TRUE)
save(m,file=file.path("Data","CVI ToxPi Tracts.RData"))

# merge on common variable
m <- merge(p, cvi.df, by.x = "GEOID10", by.y = "GEOID.Tract")
# save as shapefile again
shapefile(m, file.path("Data","CVI Tracts","CVItracts.shp"),overwrite=TRUE)
save(m,file=file.path("Data","CVI Tracts.RData"))

