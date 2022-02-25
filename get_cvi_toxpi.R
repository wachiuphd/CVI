library(data.table)
library(toxpiR)
library(dplyr)
library(ggplot2)
library(GGally)
library(grid)

indicators.df<-fread("CVI_indicators_current.csv")
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE)

## Rank correlations
pdf("CVI-corr.pdf")
cvi.abbr <- cvi.df[,-(1:5)]
names(cvi.abbr) <- abbreviate(names(cvi.abbr),minlength=8)
pp<-ggcorr(data=slice_sample(cvi.abbr,n=7000),
           method=c("pairwise.complete.obs","spearman"))
print(pp)
dev.off()

# 
# # Standardize to z-score
# cvi.scale.df<-scale(cvi.df[,-(1:5)]) # scale
# cvi.scale.df<-sweep(cvi.scale.df,2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
# cvi.scale.df<-cbind(cvi.df[,1:5],cvi.scale.df) # reattach geoids
# boxplot(as.list(cvi.scale.df[,-(1:5)]),
#         horizontal = TRUE,pars=list(outpch=15,cex=0.3))

# Standardize to percentile from 0 to 1 - keep NA values for now as NA
# ToxPi will treat NA as zero by default.
# Some indicators we will replace NA with 0.5 (=median)
cvi.pct.df<-sweep(cvi.df[,-(1:5)],2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.pct.df<-apply(cvi.pct.df,2,rank,ties.method="min",na.last="keep") # rank
cvi.pct.df<-sweep(cvi.pct.df-1,2,apply(cvi.pct.df,2,max,na.rm=T)-1,"/") # turn into percentile 0-1
cvi.pct.df<-cbind(cvi.df[,1:5], cvi.pct.df) #reattach geoids
pdf("CVI-pct.pdf")
boxplot(as.list(cvi.pct.df[,-(1:5)]),
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
f.results <- txpCalculateScores(model=f.model,
                                input=cvi.pct.df,
                                id.var="GEOID.Tract")
indx <- txpRanks(f.results)<=20
cvi.pct.toxpi <- cbind(cvi.pct.df[,1:5],
                       data.table(ToxPiScore=f.results@txpScores),
                       data.table(f.results@txpSliceScores))


pdf("ToxPi-allinone.pdf")
plot(f.results[indx,])
plot(f.results,y=txpRanks(f.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(f.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()

# Second level ToxPi - Subcategory ToxPis first

indicators.bysubcat.list <- list()
subcategories.list <- list()
f.slices.list <- list()
f.model.list <- list()
f.results.list <- list()

# Each category separately
pdf("ToxPi-subcat.pdf")
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
  fcat.results <- txpCalculateScores(model=fcat.model,
                                  input=cvi.pct.df,
                                  id.var="GEOID.Tract")
  indx <- txpRanks(fcat.results)<=sort(txpRanks(fcat.results))[20]
  plot(fcat.results[indx,],name=onecat)
  plot(fcat.results,y=txpRanks(fcat.results),name=onecat)
  grid.text(onecat,y=0.95,x=0.5,just="top")
  if (ncol(fcat.results@txpSliceScores)>1) {
    pp<-ggpairs(slice_sample(as.data.table(fcat.results@txpSliceScores),
                             n=7000),
                lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
    print(pp)
  }
  subcategories.list[[i]]<-subcategories
  indicators.bysubcat.list[[i]]<-indicators.bysubcat
  f.slices.list[[i]]<-fcat.slices
  f.model.list[[i]]<-fcat.model
  f.results.list[[i]]<-fcat.results
}

cvi.cat.scores <- cvi.pct.df[,1:5]
for (i in 1:length(categories)) {
  cvi.cat.scores <- cbind(cvi.cat.scores,f.results.list[[i]]@txpScores)
}
names(cvi.cat.scores) <- c(names(cvi.cat.scores)[1:5],categories)

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
fcomb.results <- txpCalculateScores(model=fcomb.model,
                                input=cvi.cat.scores,
                                id.var="GEOID.Tract")
indx <- txpRanks(fcomb.results)<=20

cvi.pct.comb.toxpi <- cbind(cvi.pct.df[,1:5],
                       data.table(ToxPiScore=fcomb.results@txpScores),
                       data.table(fcomb.results@txpSliceScores))
plot(fcomb.results[indx,])
plot(fcomb.results,y=txpRanks(fcomb.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(slice_sample(as.data.table(fcomb.results@txpSliceScores),
                         n=7000),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)
dev.off()
