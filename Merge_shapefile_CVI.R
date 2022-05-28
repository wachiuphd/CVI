# Merge ToxPi Scores and Data into 2010 shapefile tracts
library(data.table)
library(sp)
library(raster)
library(rgdal) ## use install.packages(c("rgdal"),repos="https://mac.R-project.org")
p <- shapefile(file.path("Data","2010 Tracts","2010tracts.shp"))

# merge overall toxpi
shpdir <- file.path("Shapefiles","CVI ToxPi Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.toxpi.df <- fread(file.path("CVI-pct","CVI-pct-comb.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
names(cvi.toxpi.df) <- gsub("[[:punct:]]", "", names(cvi.toxpi.df))
m <- merge(p, cvi.toxpi.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIToxPiTracts.shp"),overwrite=TRUE)

# merge baseline toxpi
shpdir <- file.path("Shapefiles","CVI Baseline ToxPi Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.toxpi.df <- fread(file.path("CVI-pct","CVI-pct-comb-baseline.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
names(cvi.toxpi.df) <- gsub("[[:punct:]]", "", names(cvi.toxpi.df))
m <- merge(p, cvi.toxpi.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIBaselineToxPiTracts.shp"),overwrite=TRUE)

# merge climate toxpi
shpdir <- file.path("Shapefiles","CVI Climate ToxPi Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.toxpi.df <- fread(file.path("CVI-pct","CVI-pct-comb-baseline.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
names(cvi.toxpi.df) <- gsub("[[:punct:]]", "", names(cvi.toxpi.df))
m <- merge(p, cvi.toxpi.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIClimateToxPiTracts.shp"),overwrite=TRUE)

# merge all pct data
shpdir <- file.path("Shapefiles","CVI Pct Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.pct.df <- fread(file.path("CVI-pct","CVI_data_pct.csv"),integer64 = "double",
                    keepLeadingZeros = TRUE)
names(cvi.pct.df) <- gsub("[[:punct:]]", "", names(cvi.pct.df))
m <- merge(p, cvi.pct.df, by.x = "GEOID10", by.y = "FIPS")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIPctTracts.shp"),overwrite=TRUE)

# merge all original data
shpdir <- file.path("Shapefiles","CVI Data Tracts")
if (!dir.exists(shpdir)) dir.create(shpdir)
cvi.df<-fread("CVI_data_current.csv",integer64 = "double",
              keepLeadingZeros = TRUE)
names(cvi.df) <- gsub("[[:punct:]]", "", names(cvi.df))
m <- merge(p, cvi.df, by.x = "GEOID10", by.y = "GEOIDTract")
# save as shapefile again
shapefile(m, file.path(shpdir,"CVIDataTracts.shp"),overwrite=TRUE)

# merge each category
indicators.df<-fread("CVI_indicators_current.csv")
categories <- unique(indicators.df$`Baseline Vulnerability`)

for (i in 1:length(categories)) { 
  onecat <- categories[i]
  onecatname <- gsub(" ","",gsub("[[:punct:]]", ".", onecat))
  print(onecat)
  # Create directory for shapefile
  catdir <- file.path("Shapefiles",
                      paste0("CVI ToxPi.",onecatname))
  if (!dir.exists(catdir)) dir.create(catdir)
  # read in one category ToxPi
  cvi.pct.toxpi.cat <- fread(file.path("CVI-pct",paste0("CVI-pct-cat-",
                                    gsub(": ","-",onecat),".csv")),
                             integer64 = "double",
                             keepLeadingZeros = TRUE)
  names(cvi.pct.toxpi.cat) <- gsub("[[:punct:]]", "", names(cvi.pct.toxpi.cat))
  m <- merge(p, cvi.pct.toxpi.cat, by.x = "GEOID10", by.y = "FIPS")
  # save as shapefile again
  shapefile(m, file.path(catdir,
                         paste0("CVIToxPi",onecatname,".shp")),
                         overwrite=TRUE)
}


