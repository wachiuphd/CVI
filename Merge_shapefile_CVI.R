# Merge ToxPi Scores and Data into 2010 shapefile tracts
library(sp)
library(raster)
library(rgdal) ## Problem installing
p <- shapefile(file.path("Data","2010 Tracts","2010tracts.shp"))
