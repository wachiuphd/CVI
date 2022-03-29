library(data.table)
indicators.df<-fread("CVI_indicators_current.csv")
categories <- unique(indicators.df$`Baseline Vulnerability`)
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)
tx.part.fips.counties <- fread("TX.part.FIPS.counties.csv",
                               keepLeadingZeros = TRUE)
indx.TX.part <- cvi.df$STATE=="TX" & (cvi.df$GEOID.County %in% tx.part.fips.counties$FIPS.County)
for (i in 1:length(categories)) {
  onecat <- categories[i]
  cvi.pct.toxpi.cat<-fread(paste0("CVI-pct-cat-",gsub(": ","-",onecat),".csv"),
                           keepLeadingZeros = TRUE)
  TX.part.cvi.pct.toxpi.cat <- cvi.pct.toxpi.cat[indx.TX.part,]
  fwrite(TX.part.cvi.pct.toxpi.cat,
         paste0("TX.part.","CVI-pct-cat-",gsub(": ","-",onecat),".csv"))
  
  cvi.pct.toxpi.cat.gis<-fread(paste0("CVI-pct-cat-",gsub(": ","-",onecat),".gis.csv"),
                           keepLeadingZeros = TRUE)
  TX.part.cvi.pct.toxpi.cat.gis <- cvi.pct.toxpi.cat.gis[indx.TX.part,]
  fwrite(TX.part.cvi.pct.toxpi.cat.gis,
         paste0("TX.part.","CVI-pct-cat-",gsub(": ","-",onecat),".gis.csv"))
  
}

cvi.pct.toxpi <- fread("CVI-pct-comb.csv",
                       keepLeadingZeros = TRUE)
TX.part.cvi.pct.toxpi <- cvi.pct.toxpi[indx.TX.part,]
fwrite(TX.part.cvi.pct.toxpi,
       paste0("TX.part.","CVI-pct-comb.csv"))

cvi.pct.toxpi.gis <- fread("CVI-pct-comb.gis.csv",
                       keepLeadingZeros = TRUE)
TX.part.cvi.pct.toxpi.gis <- cvi.pct.toxpi.gis[indx.TX.part,]
fwrite(TX.part.cvi.pct.toxpi.gis,
       paste0("TX.part.","CVI-pct-comb.gis.csv"))
