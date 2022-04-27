# Merge ToxPi Scores and Data into 2010 shapefile tracts
library(sp)
library(raster)
library(rgdal) ## use install.packages(c("rgdal"),repos="https://mac.R-project.org")
p <- shapefile(file.path("Data","2010 Tracts","2010tracts.shp"))

# # Load CSV file
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)
# # merge on common variable, here called 'key'
# m <- merge(p, d, by.x = "GEOID10", by.y = "GEOID.Tract")
# # save as shapefile again
# shapefile(m, "path/merged.shp")
