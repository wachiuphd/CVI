library(data.table)
library(toxpiR)
library(dplyr)
library(ggplot2)
library(GGally)
library(grid)
library(moments)

indicators.df<-fread("CVI_indicators_current.csv")
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)

nareplcols <- indicators.df$`Indicator Name`[indicators.df$`Replace NA with median`==1]
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

na0cols <- indicators.df$`Indicator Name`[indicators.df$`Replace NA with median`==0]

# other columns replace NA with 0
cvi.df[, (na0cols) := lapply(.SD, function(x) nafill(x, type = "const", fill = 0))
       , .SDcols = na0cols]
print(as.numeric((apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))))


## Rank correlations
pdf("CVI-corr.pdf")
cvi.abbr <- cvi.df[,-(1:6)]
names(cvi.abbr) <- abbreviate(names(cvi.abbr),minlength=8)
pp<-ggcorr(data=slice_sample(cvi.abbr,n=7000),
           method=c("pairwise.complete.obs","spearman"))
print(pp)
dev.off()

# 
# # Standardize to 0-1 without transformations
cvi.scale.df<-sweep(cvi.df[,-(1:6)],2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.scale.df<-scale(cvi.scale.df,center=apply(cvi.scale.df,2,min),
                    scale=apply(cvi.scale.df,2,FUN=function(x){diff(range(x))})) 
cvi.scale.df<-cbind(cvi.df[,1:6],cvi.scale.df) # reattach geoids
fwrite(cvi.scale.df,"CVI_data_scale.csv")
pdf("CVI-scale.pdf",height=10,width=7)
boxplot(as.list(cvi.scale.df[,-(1:6)]),
        horizontal = TRUE,pars=list(outpch=15,cex=0.3))
dev.off()

# Standardize to percentile from 0 to 1
# ToxPi will treat NA as zero by default
cvi.pct.df<-sweep(cvi.df[,-(1:6)],2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.pct.df<-apply(cvi.pct.df,2,rank,ties.method="min",na.last="keep") # rank
cvi.pct.df<-sweep(cvi.pct.df-1,2,apply(cvi.pct.df,2,max,na.rm=T)-1,"/") # turn into percentile 0-1
cvi.pct.df<-cbind(cvi.df[,1:6], cvi.pct.df) #reattach geoids
fwrite(cvi.pct.df,"CVI_data_pct.csv")
pdf("CVI-pct.pdf",height=10,width=7)
boxplot(as.list(cvi.pct.df[,-(1:6)]),
        horizontal = TRUE,pars=list(outpch=15,cex=0.3))
dev.off()

# Simple ToxPi - by Category only - all indicators equal weights, each category equal weight
categories <- unique(indicators.df$`Baseline Vulnerability`)
indicators.bycat <- list()
for (i in 1:length(categories)) {
  onecat <- categories[i]
  indicators.bycat[[i]] <- indicators.df$`Indicator Name`[
    indicators.df$`Baseline Vulnerability`==onecat
  ]
}

f.slices <- TxpSliceList(Baseline.Health=
                           TxpSlice(txpValueNames = indicators.bycat[[1]]),
                         Baseline.Social=
                           TxpSlice(txpValueNames = indicators.bycat[[2]]),
                         Baseline.Infrastructure=
                           TxpSlice(txpValueNames = indicators.bycat[[3]]),
                         Baseline.Environment=
                           TxpSlice(txpValueNames = indicators.bycat[[4]]),
                         ClimateChange.Health=
                           TxpSlice(txpValueNames = indicators.bycat[[5]]),
                         ClimateChange.Economic=
                           TxpSlice(txpValueNames = indicators.bycat[[6]]),
                         ClimateChange.ExtremeEvents=
                           TxpSlice(txpValueNames = indicators.bycat[[7]]))
f.model <-TxpModel(txpSlices = f.slices)

## Scale results
f.scale.results <- txpCalculateScores(model=f.model,
                                input=cvi.scale.df,
                                id.var="GEOID.Tract")
indx <- txpRanks(f.scale.results)<=10 |
  txpRanks(f.scale.results)>=(max(txpRanks(f.scale.results))-9)
cvi.scale.toxpi <- cbind(cvi.scale.df[,1:6],
                       data.table(ToxPiScore=f.scale.results@txpScores),
                       data.table(f.scale.results@txpSliceScores))
fwrite(cvi.scale.toxpi,"CVI-scale-allinone.csv")

pdf("ToxPi-scale-allinone.pdf",height=8,width=10)
plot(f.scale.results[indx,])
plot(f.scale.results,y=txpRanks(f.scale.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(f.scale.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()

## Pct results
f.pct.results <- txpCalculateScores(model=f.model,
                                      input=cvi.pct.df,
                                      id.var="GEOID.Tract")
indx <- txpRanks(f.pct.results)<=10 |
  txpRanks(f.pct.results)>=(max(txpRanks(f.pct.results))-9)
cvi.pct.toxpi <- cbind(cvi.pct.df[,1:6],
                         data.table(ToxPiScore=f.pct.results@txpScores),
                         data.table(f.pct.results@txpSliceScores))
fwrite(cvi.pct.toxpi,"CVI-pct-allinone.csv")

pdf("ToxPi-pct-allinone.pdf",height=8,width=10)
plot(f.pct.results[indx,])
plot(f.pct.results,y=txpRanks(f.pct.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(f.pct.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()

plot(f.pct.results@txpScores,f.scale.results@txpScores)
plot(f.pct.results@txpRanks,f.scale.results@txpRanks)
cor(f.pct.results@txpScores,f.scale.results@txpScores)
cor(f.pct.results@txpScores,f.scale.results@txpScores,method = "spearman")

# Second level ToxPi - Subcategory ToxPis first

indicators.bysubcat.list <- list()
subcategories.list <- list()
f.slices.list <- list()
f.model.list <- list()
f.scale.results.list <- list()
f.pct.results.list <- list()

# Each category separately
pdf("ToxPi-subcat.pdf",height=8,width=10)
for (i in 1:length(categories)) {
  onecat <- categories[i]
  subcategories <- unique(indicators.df$Subcategory[
    indicators.df$`Baseline Vulnerability`==onecat])
  indicators.bysubcat <- list()
  for (j in 1:length(subcategories)) {
    onesubcat <- subcategories[j]
    indicators.bysubcat[[j]] <- indicators.df$`Indicator Name`[
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
                                         id.var="GEOID.Tract")
  indx.scale <- txpRanks(fcat.scale.results)<=10 |
    txpRanks(fcat.scale.results)>=(max(txpRanks(fcat.scale.results))-9)
  plot(fcat.scale.results[indx.scale,],name=onecat)
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
  cvi.scale.toxpi.cat <- cbind(cvi.scale.df[,1:6],
                         data.table(ToxPiScore=fcat.scale.results@txpScores),
                         data.table(fcat.scale.results@txpSliceScores))
  fwrite(cvi.scale.toxpi.cat,paste0("CVI-scale-cat-",
                                    gsub(": ","-",onecat),".csv"))
  ## pct
  fcat.pct.results <- txpCalculateScores(model=fcat.model,
                                           input=cvi.pct.df,
                                           id.var="GEOID.Tract")
  indx.pct <- txpRanks(fcat.pct.results)<=10 |
    txpRanks(fcat.pct.results)>=(max(txpRanks(fcat.pct.results))-9)
  plot(fcat.pct.results[indx.pct,],name=onecat)
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
  cvi.pct.toxpi.cat <- cbind(cvi.pct.df[,1:6],
                               data.table(ToxPiScore=fcat.pct.results@txpScores),
                               data.table(fcat.pct.results@txpSliceScores))
  fwrite(cvi.pct.toxpi.cat,paste0("CVI-pct-cat-",
                                  gsub(": ","-",onecat),".csv"))
}
dev.off()

## Combined 
fcomb.slices <- TxpSliceList(Baseline.Health=
                               TxpSlice(txpValueNames = categories[1]),
                             Baseline.Social=
                               TxpSlice(txpValueNames = categories[2]),
                             Baseline.Infrastructure=
                               TxpSlice(txpValueNames = categories[3]),
                             Baseline.Environment=
                               TxpSlice(txpValueNames = categories[4]),
                             ClimateChange.Health=
                               TxpSlice(txpValueNames = categories[5]),
                             ClimateChange.Economic=
                               TxpSlice(txpValueNames = categories[6]),
                             ClimateChange.ExtremeEvents=
                               TxpSlice(txpValueNames = categories[7]))
fcomb.model <-TxpModel(txpSlices = fcomb.slices)

pdf("ToxPi-subcat-comb.pdf",height=8,width=10)
## Scale
cvi.scale.cat.scores <- cvi.scale.df[,1:6]
for (i in 1:length(categories)) {
  cvi.scale.cat.scores <- cbind(cvi.scale.cat.scores,f.scale.results.list[[i]]@txpScores)
}
names(cvi.scale.cat.scores) <- c(names(cvi.scale.cat.scores)[1:6],categories)
fcomb.scale.results <- txpCalculateScores(model=fcomb.model,
                                        input=cvi.scale.cat.scores,
                                        id.var="GEOID.Tract")
indx.scale <- txpRanks(fcomb.scale.results)<=10 |
  txpRanks(fcomb.scale.results)>=(max(txpRanks(fcomb.scale.results))-9)

cvi.scale.comb.toxpi <- cbind(cvi.scale.df[,1:6],
                            data.table(ToxPiScore=fcomb.scale.results@txpScores),
                            data.table(fcomb.scale.results@txpSliceScores))
fwrite(cvi.scale.comb.toxpi,paste0("CVI-scale-comb.csv"))

plot(fcomb.scale.results[indx.scale,])
plot(fcomb.scale.results,y=txpRanks(fcomb.scale.results))
grid.text(paste("Scale","Overall"),y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(fcomb.scale.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)

## Pct
cvi.pct.cat.scores <- cvi.pct.df[,1:6]
for (i in 1:length(categories)) {
  cvi.pct.cat.scores <- cbind(cvi.pct.cat.scores,f.pct.results.list[[i]]@txpScores)
}
names(cvi.pct.cat.scores) <- c(names(cvi.pct.cat.scores)[1:6],categories)
fcomb.pct.results <- txpCalculateScores(model=fcomb.model,
                                input=cvi.pct.cat.scores,
                                id.var="GEOID.Tract")
indx.pct <- txpRanks(fcomb.pct.results)<=10 |
  txpRanks(fcomb.pct.results)>=(max(txpRanks(fcomb.pct.results))-9)

cvi.pct.comb.toxpi <- cbind(cvi.pct.df[,1:6],
                       data.table(ToxPiScore=fcomb.pct.results@txpScores),
                       data.table(fcomb.pct.results@txpSliceScores))
fwrite(cvi.pct.comb.toxpi,paste0("CVI-pct-comb.csv"))
plot(fcomb.pct.results[indx.pct,])
plot(fcomb.pct.results,y=txpRanks(fcomb.pct.results))
grid.text(paste("Pct","Overall"),y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(fcomb.pct.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)

dev.off()

plot(fcomb.pct.results@txpScores,fcomb.scale.results@txpScores)
plot(fcomb.pct.results@txpRanks,fcomb.scale.results@txpRanks)
cor(fcomb.pct.results@txpScores,fcomb.scale.results@txpScores)
cor(fcomb.pct.results@txpScores,fcomb.scale.results@txpScores,method = "spearman")
