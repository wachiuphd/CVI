library(data.table)
library(toxpiR)
library(dplyr)
library(ggplot2)
library(GGally)
library(grid)
library(moments)

######## Data post-processing
diagdir <- "Diagnostics"
if (!dir.exists(diagdir)) dir.create(diagdir)

# 10 colors, color-blind friendly (# removed, all lower case)
Tol_muted <- tolower(c('88CCEE', '44AA99', '117733', '332288', 'DDCC77', '999933','CC6677', '882255', 'AA4499', 'DDDDDD'))

indicators.df<-fread("CVI_indicators_current.csv")
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)
pdf(file.path(diagdir,"CheckDist.pdf"))
for (j in 1:nrow(indicators.df)) {
  par(mfcol=c(2,2))
  y <- cvi.df[[j+6]]
  y <- y[!is.na(y)]
  try({qqnorm(y,main="",pch=15,cex=0.2); qqline(y);})
  try(hist(y,main="",xlab="Value"))
  try({qqnorm(log(y[y>0]),main="",pch=15,cex=0.2); qqline(log(y[y>0]));})
  try(hist(log(y[y>0]),main="",xlab="Log Value [positive only]"))
  mtext(paste("Row",j,"\n",indicators.df$Parameters[j]),outer=TRUE,line=-2,cex=0.75)
}
dev.off()

nareplcols <- indicators.df$Parameters[indicators.df$`Replace NA with median`==1]
# View((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols])
print(as.numeric((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))
# Replace NA by county median
cvi.df[, (nareplcols) := lapply(.SD, function(x) nafill(x, type = "const", fill = median(x, na.rm = TRUE)))
    , by = GEOID.County
    , .SDcols = nareplcols]
print(as.numeric((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))
# if still NA replace remaining by state median
cvi.df[, (nareplcols) := lapply(.SD, function(x) nafill(x, type = "const", fill = median(x, na.rm = TRUE)))
       , by = GEOID.State
       , .SDcols = nareplcols]
print(as.numeric((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))
# if still NA replace remaining by overall median
cvi.df[, (nareplcols) := lapply(.SD, function(x) nafill(x, type = "const", fill = median(x, na.rm = TRUE)))
       , .SDcols = nareplcols]
print(as.numeric((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))

na0cols <- indicators.df$Parameters[indicators.df$`Replace NA with median`==0]

# other columns replace NA with 0
cvi.df[, (na0cols) := lapply(.SD, function(x) nafill(x, type = "const", fill = 0))
       , .SDcols = na0cols]
print(as.numeric((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))))

cvi.dat.df <- cvi.df[,-(1:6)]

## Rank correlations
pdf(file.path(diagdir,"CVI-corr.pdf"))
cvi.abbr <- cvi.df[,-(1:6)]
names(cvi.abbr) <- abbreviate(names(cvi.abbr),minlength=8)
pp<-ggcorr(data=slice_sample(cvi.abbr,n=7000),
           method=c("pairwise.complete.obs","spearman"))
print(pp)
dev.off()

# Absolute value when adverse direction is absolute value
na_adverse <- which(is.na(indicators.df$`Adverse Direction`))
if (length(na_adverse) > 0) {
  cvi.dat.df <- as.data.frame(cvi.dat.df)
  cvi.dat.df[,na_adverse] <- abs(cvi.dat.df[,na_adverse])
  cvi.dat.df <- as.data.table(cvi.dat.df)
  indicators.df$`Adverse Direction`[na_adverse] <- 1
}

# idcols for use in GUI
idcols_gui.df <- data.table(
  `row#` = 1:nrow(cvi.df),
  Name=paste0(cvi.df$STATE,", ",cvi.df$County_Name,", ",cvi.df$GEOID.Tract),
  FIPS=cvi.df$GEOID.Tract,
  Source=cvi.df$LatLong
)

idcols_gis.df <- data.table(
  FIPS=cvi.df$GEOID.Tract,
  Name=paste0(cvi.df$STATE,", ",cvi.df$County_Name,", ",cvi.df$GEOID.Tract),
  Source=cvi.df$LatLong
)

############ Percentiles
pctdir <- "CVI-pct"
if (!dir.exists(pctdir)) dir.create(pctdir)

# Standardize to percentile from 0 to 1
# ToxPi will treat NA as zero by default
cvi.pct.df<-sweep(cvi.dat.df,2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.pct.df<-as.data.frame(apply(cvi.pct.df,2,rank,ties.method="min",na.last="keep")) # rank
cvi.pct.df<-sweep(cvi.pct.df-1,2,apply(cvi.pct.df,2,max,na.rm=T)-1,"/") # turn into percentile 0-1
pdf(file.path(pctdir,"CVI-pct.pdf"),height=10,width=7)
boxplot(as.list(cvi.pct.df),
        horizontal = TRUE,pars=list(outpch=15,cex=0.3))
dev.off()
cvi.pct.df<-cbind(idcols_gui.df, cvi.pct.df) # save for use by ToxPi GUI
fwrite(cvi.pct.df,file.path(pctdir,"CVI_data_pct.csv"),quote=TRUE)

# Simple ToxPi - by Category only - all indicators equal weights, each category equal weight
categories <- unique(indicators.df$`Baseline Vulnerability`)
indicators.bycat <- list()
for (i in 1:length(categories)) {
  onecat <- categories[i]
  indicators.bycat[[i]] <- indicators.df$Parameters[
    indicators.df$`Baseline Vulnerability`==onecat
  ]
}

f.slices <- TxpSliceList(Baseline.Health=
                           TxpSlice(txpValueNames = indicators.bycat[[1]]),
                         Baseline.SocialEconomic=
                           TxpSlice(txpValueNames = indicators.bycat[[2]]),
                         Baseline.Infrastructure=
                           TxpSlice(txpValueNames = indicators.bycat[[3]]),
                         Baseline.Environment=
                           TxpSlice(txpValueNames = indicators.bycat[[4]]),
                         ClimateChange.Health=
                           TxpSlice(txpValueNames = indicators.bycat[[5]]),
                         ClimateChange.SocialEconomic=
                           TxpSlice(txpValueNames = indicators.bycat[[6]]),
                         ClimateChange.ExtremeEvents=
                           TxpSlice(txpValueNames = indicators.bycat[[7]]))
f.model <-TxpModel(txpSlices = f.slices)

## Pct results
f.pct.results <- txpCalculateScores(model=f.model,
                                      input=cvi.pct.df,
                                      id.var="FIPS")
indx <- txpRanks(f.pct.results)<=10 |
  txpRanks(f.pct.results)>=(max(txpRanks(f.pct.results))-9)

# save for use by ToxPi GUI
cvi.pct.toxpi <- cbind(idcols_gui.df, 
                         data.table(`ToxPi Score`=f.pct.results@txpScores),
                         data.table(f.pct.results@txpSliceScores)) 
fwrite(cvi.pct.toxpi,file.path(pctdir,"CVI-pct-allinone.csv"))
# save for use by ToxPi GIS
cvi.pct.toxpi.gis <- cbind(data.table(`ToxPi Score`=f.pct.results@txpScores),
                             idcols_gis.df, 
                             data.table(f.pct.results@txpSliceScores)) 
slicenames <- names(f.slices@listData)
newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
setnames(cvi.pct.toxpi.gis,slicenames,newslicenames)
fwrite(cvi.pct.toxpi.gis,file.path(pctdir,"CVI-pct-allinone.gis.csv"))

pdf(file.path(pctdir,"ToxPi-pct-allinone.pdf"),height=8,width=10)
plot(f.pct.results[indx,],fills=paste0("#",Tol_muted))
plot(f.pct.results,y=txpRanks(f.pct.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(f.pct.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()

# Second level ToxPi - Subcategory ToxPis first

indicators.bysubcat.list <- list()
subcategories.list <- list()
f.slices.list <- list()
f.model.list <- list()
f.pct.results.list <- list()

# Each category separately
pdf(file.path(pctdir,"ToxPi-pct-subcat.pdf"),height=8,width=10)
for (i in 1:length(categories)) {
  onecat <- categories[i]
  subcategories <- unique(indicators.df$Subcategory[
    indicators.df$`Baseline Vulnerability`==onecat])
  indicators.bysubcat <- list()
  for (j in 1:length(subcategories)) {
    onesubcat <- subcategories[j]
    indicators.bysubcat[[j]] <- indicators.df$Parameters[
      indicators.df$`Baseline Vulnerability`==onecat &
        indicators.df$Subcategory==onesubcat
    ]
  }
  fcat.slices <- list()
  for (j in 1:length(subcategories)) {
    fcat.slices[[j]] <- TxpSlice(txpValueNames = indicators.bysubcat[[j]])
  }
  names(fcat.slices)<-subcategories
  fcat.slices<-as.TxpSliceList(fcat.slices)
  fcat.model <-TxpModel(txpSlices = fcat.slices)
  ## pct
  fcat.pct.results <- txpCalculateScores(model=fcat.model,
                                           input=cvi.pct.df,
                                           id.var="FIPS")
  indx.pct <- txpRanks(fcat.pct.results)<=10 |
    txpRanks(fcat.pct.results)>=(max(txpRanks(fcat.pct.results))-9)
  plot(fcat.pct.results[indx.pct,],name=onecat,fills=paste0("#",Tol_muted))
  grid.text("Pct",y=0.05,x=0.5,just="top")
  plot(fcat.pct.results,y=txpRanks(fcat.pct.results),name=onecat)
  grid.text(paste("Pct",onecat),y=0.95,x=0.5,just="top")
  if (ncol(fcat.pct.results@txpSliceScores)>1) {
    pp<-ggpairs(slice_sample(as.data.table(fcat.pct.results@txpSliceScores),
                             n=7000),
                lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
    print(pp)
  }
  subcategories.list[[i]]<-subcategories
  indicators.bysubcat.list[[i]]<-indicators.bysubcat
  f.slices.list[[i]]<-fcat.slices
  f.model.list[[i]]<-fcat.model
  f.pct.results.list[[i]]<-fcat.pct.results
  # save for use by ToxPi GUI
  cvi.pct.toxpi.cat <- cbind(idcols_gui.df, 
                               data.table(`ToxPi Score`=fcat.pct.results@txpScores),
                               data.table(fcat.pct.results@txpSliceScores)) 
  fwrite(cvi.pct.toxpi.cat,file.path(pctdir,
                                     paste0("CVI-pct-cat-",
                                    gsub(": ","-",onecat),".csv")))
  # save for use by ToxPi GIS
  cvi.pct.toxpi.cat.gis <- cbind(data.table(`ToxPi Score`=fcat.pct.results@txpScores),
                                   idcols_gis.df, 
                                   data.table(fcat.pct.results@txpSliceScores)) 
  slicenames <- names(fcat.slices@listData)
  newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
  setnames(cvi.pct.toxpi.cat.gis,slicenames,newslicenames)
  fwrite(cvi.pct.toxpi.cat.gis,file.path(pctdir,
                                         paste0("CVI-pct-cat-",
                                        gsub(": ","-",onecat),".gis.csv")))
}
dev.off()

## Combined 
fcomb.slices <- TxpSliceList(Baseline.Health=
                               TxpSlice(txpValueNames = categories[1]),
                             Baseline.SocialEconomic=
                               TxpSlice(txpValueNames = categories[2]),
                             Baseline.Infrastructure=
                               TxpSlice(txpValueNames = categories[3]),
                             Baseline.Environment=
                               TxpSlice(txpValueNames = categories[4]),
                             ClimateChange.Health=
                               TxpSlice(txpValueNames = categories[5]),
                             ClimateChange.SocialEconomic=
                               TxpSlice(txpValueNames = categories[6]),
                             ClimateChange.ExtremeEvents=
                               TxpSlice(txpValueNames = categories[7]))
fcomb.model <-TxpModel(txpSlices = fcomb.slices)

pdf(file.path(pctdir,"ToxPi-pct-subcat-comb.pdf"),height=8,width=10)

## Pct
cvi.pct.cat.scores <- idcols_gui.df
for (i in 1:length(categories)) {
  cvi.pct.cat.scores <- cbind(cvi.pct.cat.scores,f.pct.results.list[[i]]@txpScores)
}
names(cvi.pct.cat.scores) <- c(names(idcols_gui.df),categories)
fcomb.pct.results <- txpCalculateScores(model=fcomb.model,
                                input=cvi.pct.cat.scores,
                                id.var="FIPS")
indx.pct <- txpRanks(fcomb.pct.results)<=10 |
  txpRanks(fcomb.pct.results)>=(max(txpRanks(fcomb.pct.results))-9)

# save for use by ToxPi GUI
cvi.pct.toxpi.comb <- cbind(idcols_gui.df, 
                              data.table(`ToxPi Score`=fcomb.pct.results@txpScores),
                              data.table(fcomb.pct.results@txpSliceScores)) 
fwrite(cvi.pct.toxpi.comb,file.path(pctdir,paste0("CVI-pct-comb.csv")))
# save for use by ToxPi GIS
cvi.pct.toxpi.comb.gis <- cbind(data.table(`ToxPi Score`=fcomb.pct.results@txpScores),
                                  idcols_gis.df, 
                                  data.table(fcomb.pct.results@txpSliceScores)) 
slicenames <- names(fcomb.slices@listData)
newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
setnames(cvi.pct.toxpi.comb.gis,slicenames,newslicenames)
fwrite(cvi.pct.toxpi.comb.gis,file.path(pctdir,paste0("CVI-pct-comb.gis.csv")))

plot(fcomb.pct.results[indx.pct,],fills=paste0("#",Tol_muted))
plot(fcomb.pct.results,y=txpRanks(fcomb.pct.results))
grid.text(paste("Pct","Overall"),y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(fcomb.pct.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)

dev.off()

## Percentile file for GIS
cvi.pct.df.namesfixed <- cvi.pct.df
names(cvi.pct.df.namesfixed)<-gsub(",","",names(cvi.pct.df.namesfixed))
names(cvi.pct.df.namesfixed)<-gsub("\"","",names(cvi.pct.df.namesfixed))
fwrite(cvi.pct.df.namesfixed,file.path(pctdir,"CVI_data_pct.gis.csv"))

########## Scale 0-1 without transformations
scaledir <- "CVI-scale"
if (!dir.exists(scaledir)) dir.create(scaledir)
# 
# # Standardize to 0-1 without transformations
cvi.scale.df<-sweep(cvi.dat.df,2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.scale.df<-as.data.frame(scale(cvi.scale.df,center=apply(cvi.scale.df,2,min),
                                  scale=apply(cvi.scale.df,2,FUN=function(x){diff(range(x))})))
pdf(file.path(scaledir,"CVI-scale.pdf"),height=10,width=7)
boxplot(as.list(cvi.scale.df),
        horizontal = TRUE,pars=list(outpch=15,cex=0.3))
dev.off()
cvi.scale.df<-cbind(idcols_gui.df,cvi.scale.df) # save for use by ToxPi GUI
fwrite(cvi.scale.df,file.path(scaledir,"CVI_data_scale.csv"),quote=TRUE)


## Scale results
f.scale.results <- txpCalculateScores(model=f.model,
                                      input=cvi.scale.df,
                                      id.var="FIPS")
indx <- txpRanks(f.scale.results)<=10 |
  txpRanks(f.scale.results)>=(max(txpRanks(f.scale.results))-9)

# save for use by ToxPi GUI
cvi.scale.toxpi <- cbind(idcols_gui.df, 
                         data.table(`ToxPi Score`=f.scale.results@txpScores),
                         data.table(f.scale.results@txpSliceScores)) 
fwrite(cvi.scale.toxpi,file.path(scaledir,"CVI-scale-allinone.csv"))
# save for use by ToxPi GIS
cvi.scale.toxpi.gis <- cbind(data.table(`ToxPi Score`=f.scale.results@txpScores),
                             idcols_gis.df, 
                             data.table(f.scale.results@txpSliceScores)) 
slicenames <- names(f.slices@listData)
newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
setnames(cvi.scale.toxpi.gis,slicenames,newslicenames)
fwrite(cvi.scale.toxpi.gis,file.path(scaledir,"CVI-scale-allinone.gis.csv"))

pdf(file.path(scaledir,"ToxPi-scale-allinone.pdf"),height=8,width=10)
plot(f.scale.results[indx,],fills=paste0("#",Tol_muted))
plot(f.scale.results,y=txpRanks(f.scale.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(f.scale.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()



f.scale.results.list <- list()
# Each category separately
pdf(file.path(scaledir,"ToxPi-scale-subcat.pdf"),height=8,width=10)
for (i in 1:length(categories)) {
  onecat <- categories[i]
  subcategories <- unique(indicators.df$Subcategory[
    indicators.df$`Baseline Vulnerability`==onecat])
  indicators.bysubcat <- list()
  for (j in 1:length(subcategories)) {
    onesubcat <- subcategories[j]
    indicators.bysubcat[[j]] <- indicators.df$Parameters[
      indicators.df$`Baseline Vulnerability`==onecat &
        indicators.df$Subcategory==onesubcat
    ]
  }
  fcat.slices <- list()
  for (j in 1:length(subcategories)) {
    fcat.slices[[j]] <- TxpSlice(txpValueNames = indicators.bysubcat[[j]])
  }
  names(fcat.slices)<-subcategories
  fcat.slices<-as.TxpSliceList(fcat.slices)
  fcat.model <-TxpModel(txpSlices = fcat.slices)
  ## Scale
  fcat.scale.results <- txpCalculateScores(model=fcat.model,
                                           input=cvi.scale.df,
                                           id.var="FIPS")
  indx.scale <- txpRanks(fcat.scale.results)<=10 |
    txpRanks(fcat.scale.results)>=(max(txpRanks(fcat.scale.results))-9)
  plot(fcat.scale.results[indx.scale,],name=onecat,fills=paste0("#",Tol_muted))  
  grid.text("Scale",y=0.05,x=0.5,just="top")
  plot(fcat.scale.results,y=txpRanks(fcat.scale.results),name=onecat)
  grid.text(paste("Scale",onecat),y=0.95,x=0.5,just="top")
  if (ncol(fcat.scale.results@txpSliceScores)>1) {
    pp<-ggpairs(slice_sample(as.data.table(fcat.scale.results@txpSliceScores),
                             n=7000),
                lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
    print(pp)
  }
  subcategories.list[[i]]<-subcategories
  indicators.bysubcat.list[[i]]<-indicators.bysubcat
  f.slices.list[[i]]<-fcat.slices
  f.model.list[[i]]<-fcat.model
  f.scale.results.list[[i]]<-fcat.scale.results
  # save for use by ToxPi GUI
  cvi.scale.toxpi.cat <- cbind(idcols_gui.df, 
                               data.table(`ToxPi Score`=fcat.scale.results@txpScores),
                               data.table(fcat.scale.results@txpSliceScores)) 
  fwrite(cvi.scale.toxpi.cat,file.path(scaledir,
                                       paste0("CVI-scale-cat-",
                                    gsub(": ","-",onecat),".csv")))
  # save for use by ToxPi GIS
  cvi.scale.toxpi.cat.gis <- cbind(data.table(`ToxPi Score`=fcat.scale.results@txpScores),
                                   idcols_gis.df, 
                                   data.table(fcat.scale.results@txpSliceScores)) 
  slicenames <- names(fcat.slices@listData)
  newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
  setnames(cvi.scale.toxpi.cat.gis,slicenames,newslicenames)
  fwrite(cvi.scale.toxpi.cat.gis,file.path(scaledir,
                                           paste0("CVI-scale-cat-",
                                        gsub(": ","-",onecat),".gis.csv")))
  
}
dev.off()

pdf(file.path(scaledir,"ToxPi-scale-subcat-comb.pdf"),height=8,width=10)

## Scale
cvi.scale.cat.scores <- idcols_gui.df
for (i in 1:length(categories)) {
  cvi.scale.cat.scores <- cbind(cvi.scale.cat.scores,f.scale.results.list[[i]]@txpScores)
}
names(cvi.scale.cat.scores) <- c(names(idcols_gui.df),categories)
fcomb.scale.results <- txpCalculateScores(model=fcomb.model,
                                          input=cvi.scale.cat.scores,
                                          id.var="FIPS")
indx.scale <- txpRanks(fcomb.scale.results)<=10 |
  txpRanks(fcomb.scale.results)>=(max(txpRanks(fcomb.scale.results))-9)


# save for use by ToxPi GUI
cvi.scale.toxpi.comb <- cbind(idcols_gui.df, 
                              data.table(`ToxPi Score`=fcomb.scale.results@txpScores),
                              data.table(fcomb.scale.results@txpSliceScores)) 
fwrite(cvi.scale.toxpi.comb,file.path(scaledir,paste0("CVI-scale-comb.csv")))
# save for use by ToxPi GIS
cvi.scale.toxpi.comb.gis <- cbind(data.table(`ToxPi Score`=fcomb.scale.results@txpScores),
                                  idcols_gis.df, 
                                  data.table(fcomb.scale.results@txpSliceScores)) 
slicenames <- names(fcomb.slices@listData)
newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
setnames(cvi.scale.toxpi.comb.gis,slicenames,newslicenames)
fwrite(cvi.scale.toxpi.comb.gis,file.path(scaledir,paste0("CVI-scale-comb.gis.csv")))


plot(fcomb.scale.results[indx.scale,],fills=paste0("#",Tol_muted))
plot(fcomb.scale.results,y=txpRanks(fcomb.scale.results))
grid.text(paste("Scale","Overall"),y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(fcomb.scale.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()

##########

plot(f.pct.results@txpScores,f.scale.results@txpScores)
plot(f.pct.results@txpRanks,f.scale.results@txpRanks)
cor(f.pct.results@txpScores,f.scale.results@txpScores)
cor(f.pct.results@txpScores,f.scale.results@txpScores,method = "spearman")

plot(fcomb.pct.results@txpScores,fcomb.scale.results@txpScores)
plot(fcomb.pct.results@txpRanks,fcomb.scale.results@txpRanks)
cor(fcomb.pct.results@txpScores,fcomb.scale.results@txpScores)
cor(fcomb.pct.results@txpScores,fcomb.scale.results@txpScores,method = "spearman")

