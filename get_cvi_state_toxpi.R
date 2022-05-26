library(data.table)
library(toxpiR)
library(dplyr)
library(ggplot2)
library(GGally)
library(grid)
library(moments)
library(tigris)

######## Data post-processing

# 10 colors, color-blind friendly (# removed, all lower case)
Tol_muted <- tolower(c('88CCEE', '44AA99', '117733', '332288', 'DDCC77', '999933','CC6677', '882255', 'AA4499', 'DDDDDD'))

indicators.df<-fread("CVI_indicators_current.csv")
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE,integer64 = "numeric")

st.df10 <- as.data.frame(states(year=2010))
st.df10 <- st.df10[,names(st.df10)%in%c("GEOID10","INTPTLAT10","INTPTLON10")]
st.df10$LatLong <- paste(st.df10$INTPTLAT10,st.df10$INTPTLON10,sep=",")
rownames(st.df10)<-st.df10$GEOID10

############ State 
cvi.state.df <- cvi.df[,c(1:6)]
cvi.state.df$LatLong <- st.df10[cvi.state.df$GEOID.State,"LatLong"]
cvi.state.df$County_Name <- ""
cvi.state.df$GEOID.County <- ""
cvi.state.df$GEOID.Tract <- ""
cvi.state.df <- cvi.state.df[!duplicated(cvi.state.df),]
for (j in 7:(ncol(cvi.df))) {
  cat(paste0(j,"...",names(cvi.df)[j],"..."))
  tmp <- aggregate(. ~ STATE+GEOID.State,
                   data=as.data.frame(cvi.df)[,c(1,3,j)],
                   FUN = median, na.rm=T)
  cvi.state.df <- left_join(cvi.state.df,tmp)
}

# Replace main data frame with state data frame
cvi.df <- cvi.state.df

###### 

nareplcols <- indicators.df$Parameters[indicators.df$`Replace NA with median`==1]
# View((base::apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols])
print(as.numeric((base::apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))
# if still NA replace remaining by state median
cvi.df[, (nareplcols) := lapply(.SD, function(x) nafill(x, type = "const", fill = median(x, na.rm = TRUE)))
       , by = GEOID.State
       , .SDcols = nareplcols]
print(as.numeric((base::apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))
# if still NA replace remaining by overall median
cvi.df[, (nareplcols) := lapply(.SD, function(x) nafill(x, type = "const", fill = median(x, na.rm = TRUE)))
       , .SDcols = nareplcols]
print(as.numeric((base::apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))[nareplcols]))

na0cols <- indicators.df$Parameters[indicators.df$`Replace NA with median`==0]

# other columns replace NA with 0
cvi.df[, (na0cols) := lapply(.SD, function(x) nafill(x, type = "const", fill = 0))
       , .SDcols = na0cols]
print(as.numeric((base::apply(cvi.df,2,FUN=function(x) {sum(is.na(x))}))))

cvi.dat.df <- cvi.df[,-(1:6)]

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
  Name=paste0(cvi.df$STATE),
  FIPS=cvi.df$GEOID.State,
  Source=cvi.df$LatLong
)

idcols_gis.df <- data.table(
  FIPS=cvi.df$GEOID.State,
  Name=paste0(cvi.df$STATE),
  Source=cvi.df$LatLong
)

############ Percentiles
pctdir <- "CVI-state-pct"
if (!dir.exists(pctdir)) dir.create(pctdir)

# Standardize to percentile from 0 to 1
# ToxPi will treat NA as zero by default
cvi.pct.df<-sweep(cvi.dat.df,2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.pct.df<-as.data.frame(base::apply(cvi.pct.df,2,rank,ties.method="min",na.last="keep")) # rank
cvi.pct.df<-sweep(cvi.pct.df-1,2,base::apply(cvi.pct.df,2,max,na.rm=T)-1,"/") # turn into percentile 0-1
pdf(file.path(pctdir,"CVI-state-pct.pdf"),height=10,width=7)
boxplot(as.list(cvi.pct.df),
        horizontal = TRUE,pars=list(outpch=15,cex=0.3))
dev.off()
cvi.pct.df<-cbind(idcols_gui.df, cvi.pct.df) # save for use by ToxPi GUI
fwrite(cvi.pct.df,file.path(pctdir,"CVI-state_data_pct.csv"),quote=TRUE)

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
fwrite(cvi.pct.toxpi,file.path(pctdir,"CVI-state-pct-allinone.csv"))
# save for use by ToxPi GIS
cvi.pct.toxpi.gis <- cbind(data.table(`ToxPi Score`=f.pct.results@txpScores),
                           idcols_gis.df, 
                           data.table(f.pct.results@txpSliceScores)) 
slicenames <- names(f.slices@listData)
newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
setnames(cvi.pct.toxpi.gis,slicenames,newslicenames)
fwrite(cvi.pct.toxpi.gis,file.path(pctdir,"CVI-state-pct-allinone.gis.csv"))

pdf(file.path(pctdir,"ToxPi-state-pct-allinone.pdf"),height=8,width=10)
plot(f.pct.results[indx,],fills=paste0("#",Tol_muted))
plot(f.pct.results,y=txpRanks(f.pct.results))
grid.text("Overall",y=0.95,x=0.5,just="top")
pp<-ggpairs(as.data.table(f.pct.results@txpSliceScores),
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
pdf(file.path(pctdir,"ToxPi-state-pct-subcat.pdf"),height=8,width=10)
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
    pp<-ggpairs(as.data.table(fcat.pct.results@txpSliceScores),
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
                                     paste0("CVI-state-pct-cat-",
                                            gsub(": ","-",onecat),".csv")))
  # save for use by ToxPi GIS
  cvi.pct.toxpi.cat.gis <- cbind(data.table(`ToxPi Score`=fcat.pct.results@txpScores),
                                 idcols_gis.df, 
                                 data.table(fcat.pct.results@txpSliceScores)) 
  slicenames <- names(fcat.slices@listData)
  newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
  setnames(cvi.pct.toxpi.cat.gis,slicenames,newslicenames)
  fwrite(cvi.pct.toxpi.cat.gis,file.path(pctdir,
                                         paste0("CVI-state-pct-cat-",
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

pdf(file.path(pctdir,"ToxPi-state-pct-subcat-comb.pdf"),height=8,width=10)

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
fwrite(cvi.pct.toxpi.comb,file.path(pctdir,paste0("CVI-state-pct-comb.csv")))
# save for use by ToxPi GIS
cvi.pct.toxpi.comb.gis <- cbind(data.table(`ToxPi Score`=fcomb.pct.results@txpScores),
                                idcols_gis.df, 
                                data.table(fcomb.pct.results@txpSliceScores)) 
slicenames <- names(fcomb.slices@listData)
newslicenames <- paste0(slicenames,"!1","!0x",Tol_muted[1:length(slicenames)],"ff")
setnames(cvi.pct.toxpi.comb.gis,slicenames,newslicenames)
fwrite(cvi.pct.toxpi.comb.gis,file.path(pctdir,paste0("CVI-state-pct-comb.gis.csv")))

plot(fcomb.pct.results[indx.pct,],fills=paste0("#",Tol_muted))
plot(fcomb.pct.results,y=txpRanks(fcomb.pct.results))
grid.text(paste("Pct","Overall"),y=0.95,x=0.5,just="top")
pp<-ggpairs(as.data.table(fcomb.pct.results@txpSliceScores),
            lower = list(continuous = wrap("points", alpha=0.1, size=0.1)))
print(pp)

dev.off()

## Percentile file for GIS
cvi.pct.df.namesfixed <- cvi.pct.df
names(cvi.pct.df.namesfixed)<-gsub(",","",names(cvi.pct.df.namesfixed))
names(cvi.pct.df.namesfixed)<-gsub("\"","",names(cvi.pct.df.namesfixed))
fwrite(cvi.pct.df.namesfixed,file.path(pctdir,"CVI-state_data_pct.gis.csv"))

