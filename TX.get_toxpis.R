library(data.table)
indicators.df<-fread("CVI_indicators_current.csv")
categories <- unique(indicators.df$`Baseline Vulnerability`)
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)
indx.TX <- cvi.df$STATE=="TX"
for (i in 1:length(categories)) {
  onecat <- categories[i]
  cvi.pct.toxpi.cat<-fread(paste0("CVI-pct-cat-",gsub(": ","-",onecat),".csv"),
                           keepLeadingZeros = TRUE)
  TX.cvi.pct.toxpi.cat <- cvi.pct.toxpi.cat[indx.TX,]
  fwrite(TX.cvi.pct.toxpi.cat,
         paste0("TX.","CVI-pct-cat-",gsub(": ","-",onecat),".csv"))
  
  cvi.pct.toxpi.cat.gis<-fread(paste0("CVI-pct-cat-",gsub(": ","-",onecat),".gis.csv"),
                           keepLeadingZeros = TRUE)
  TX.cvi.pct.toxpi.cat.gis <- cvi.pct.toxpi.cat.gis[indx.TX,]
  fwrite(TX.cvi.pct.toxpi.cat.gis,
         paste0("TX.","CVI-pct-cat-",gsub(": ","-",onecat),".gis.csv"))
  
}

cvi.pct.toxpi <- fread("CVI-pct-comb.csv",
                       keepLeadingZeros = TRUE)
TX.cvi.pct.toxpi <- cvi.pct.toxpi[indx.TX,]
fwrite(TX.cvi.pct.toxpi,
       paste0("TX.","CVI-pct-comb.csv"))

cvi.pct.toxpi.gis <- fread("CVI-pct-comb.gis.csv",
                       keepLeadingZeros = TRUE)
TX.cvi.pct.toxpi.gis <- cvi.pct.toxpi.gis[indx.TX,]
fwrite(TX.cvi.pct.toxpi.gis,
       paste0("TX.","CVI-pct-comb.gis.csv"))
