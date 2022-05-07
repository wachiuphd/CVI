library(data.table)
Will.IL.dir <- "Will.IL"
if (!dir.exists(Will.IL.dir)) dir.create(Will.IL.dir)
indicators.df<-fread("CVI_indicators_current.csv")
categories <- unique(indicators.df$`Baseline Vulnerability`)
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)
indx.will.il <- cvi.df$STATE=="IL" & cvi.df$County_Name=="Will"

for (i in 1:length(categories)) {
  onecat <- categories[i]
  cvi.pct.toxpi.cat<-fread(paste0("CVI-pct-cat-",gsub(": ","-",onecat),".csv"),
                           keepLeadingZeros = TRUE)
  will.il.cvi.pct.toxpi.cat <- cvi.pct.toxpi.cat[indx.will.il,]
  fwrite(will.il.cvi.pct.toxpi.cat,
         file.path(Will.IL.dir,
                   paste0("Will.IL.","CVI-pct-cat-",
                          gsub(": ","-",onecat),".csv")))
  
  cvi.pct.toxpi.cat.gis<-fread(paste0("CVI-pct-cat-",gsub(": ","-",onecat),".gis.csv"),
                           keepLeadingZeros = TRUE)
  will.il.cvi.pct.toxpi.cat.gis <- cvi.pct.toxpi.cat.gis[indx.will.il,]
  fwrite(will.il.cvi.pct.toxpi.cat.gis,
         file.path(Will.IL.dir,
                   paste0("Will.IL.","CVI-pct-cat-",
                          gsub(": ","-",onecat),".gis.csv")))
  
}

cvi.pct.toxpi <- fread("CVI-pct-comb.csv",
                       keepLeadingZeros = TRUE)
will.il.cvi.pct.toxpi <- cvi.pct.toxpi[indx.will.il,]
fwrite(will.il.cvi.pct.toxpi,
       file.path(Will.IL.dir,
                 paste0("Will.IL.","CVI-pct-comb.csv")))

cvi.pct.toxpi.gis <- fread("CVI-pct-comb.gis.csv",
                       keepLeadingZeros = TRUE)
will.il.cvi.pct.toxpi.gis <- cvi.pct.toxpi.gis[indx.will.il,]
fwrite(will.il.cvi.pct.toxpi.gis,
       file.path(Will.IL.dir,
                 paste0("Will.IL.","CVI-pct-comb.gis.csv")))
