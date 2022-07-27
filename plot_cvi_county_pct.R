library(data.table)
library(choroplethr)
library(choroplethrMaps)
library(ggplot2)
library(tigris)
library(dplyr)

data(df_pop_county)
pctdat <- fread(file.path("CVI-county-pct","CVI-county_data_pct.csv"),
                keepLeadingZeros = TRUE)

pdf(file.path("CVI-county-pct","CVI-county_data_pct.pdf"),height=6,width=10)
for (j in 5:(ncol(pctdat))) {
  dat.df <- data.frame(region=as.numeric(pctdat$FIPS),
                       value=100*pctdat[[j]])
  plt<-CountyChoropleth$new(dat.df)
  plt$set_num_colors(1)
  plt$set_zoom(NULL)
  #plt$ggplot_scale <- scale_fill_viridis_c("County\nVulnerability\nPercentile",option="magma")
  plt$title<-paste0(names(pctdat)[j],"\n(median of census tracts in county)")
  plt$ggplot_polygon <- geom_polygon(aes(fill = value),color=NA)
  p<-plt$render()
  print(p)
}
dev.off()

############ Max census tract map

indicators.df<-fread("CVI_indicators_current.csv")
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE,integer64 = "numeric")

ct.df10 <- as.data.frame(counties(year=2010))
ct.df10 <- ct.df10[,names(ct.df10)%in%c("GEOID10","INTPTLAT10","INTPTLON10")]
ct.df10$LatLong <- paste(ct.df10$INTPTLAT10,ct.df10$INTPTLON10,sep=",")
rownames(ct.df10)<-ct.df10$GEOID10
############ County - MAX
cvi.county.df <- cvi.df[,c(1:6)]
cvi.county.df$LatLong <- ct.df10[cvi.county.df$GEOID.County,"LatLong"]
cvi.county.df$GEOID.Tract <- ""
cvi.county.df <- cvi.county.df[!duplicated(cvi.county.df),]
for (j in 7:(ncol(cvi.df))) {
  cat(paste0(j,"...")) 
  if (is.na(indicators.df$`Adverse Direction`[j-6])) {
    cat("||...")
    tmp <- aggregate(. ~ STATE+County_Name+GEOID.State+GEOID.County,
                     data=as.data.frame(cvi.df)[,c(1,2,3,4,j)],
                     FUN = function (x) {max(abs(x), na.rm=T)})
  } else if (indicators.df$`Adverse Direction`[j-6] == 1) {
    cat("+...")
    tmp <- aggregate(. ~ STATE+County_Name+GEOID.State+GEOID.County,
                     data=as.data.frame(cvi.df)[,c(1,2,3,4,j)],
                     FUN = max, na.rm=T)
  } else if (indicators.df$`Adverse Direction`[j-6] == -1) {
    cat("-...")
    tmp <- aggregate(. ~ STATE+County_Name+GEOID.State+GEOID.County,
                     data=as.data.frame(cvi.df)[,c(1,2,3,4,j)],
                     FUN = min, na.rm=T)
  } 
  cat(paste0(names(cvi.df)[j],"..."))
  cvi.county.df <- left_join(cvi.county.df,tmp)
}

# Replace main data frame with county data frame
cvi.df <- cvi.county.df


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
  Name=paste0(cvi.df$STATE,", ",cvi.df$County_Name),
  FIPS=cvi.df$GEOID.County,
  Source=cvi.df$LatLong
)

idcols_gis.df <- data.table(
  FIPS=cvi.df$GEOID.County,
  Name=paste0(cvi.df$STATE,", ",cvi.df$County_Name),
  Source=cvi.df$LatLong
)

############ Percentiles
pctdir <- "CVI-county-pct"
if (!dir.exists(pctdir)) dir.create(pctdir)

# Standardize to percentile from 0 to 1
# ToxPi will treat NA as zero by default
cvi.pct.df<-sweep(cvi.dat.df,2,indicators.df$`Adverse Direction`,"*") # multiple by adverse direction
cvi.pct.df<-as.data.frame(base::apply(cvi.pct.df,2,rank,ties.method="min",na.last="keep")) # rank
cvi.pct.df<-sweep(cvi.pct.df-1,2,base::apply(cvi.pct.df,2,max,na.rm=T)-1,"/") # turn into percentile 0-1
cvi.pct.df<-cbind(idcols_gui.df, cvi.pct.df) # save for use by ToxPi GUI
fwrite(cvi.pct.df,file.path(pctdir,"CVI-county_data_maxpct.csv"),quote=TRUE)

##

pdf(file.path("CVI-county-pct","CVI-county_data_maxpct.pdf"),height=6,width=10)
for (j in 5:(ncol(cvi.pct.df))) {
  dat.df <- data.frame(region=as.numeric(cvi.pct.df$FIPS),
                       value=100*cvi.pct.df[[j]])
  plt<-CountyChoropleth$new(dat.df)
  plt$set_num_colors(1)
  plt$set_zoom(NULL)
  # plt$ggplot_scale <- scale_fill_viridis_c("County\nVulnerability\nPercentile",option="magma")
  plt$title<-paste0(names(cvi.pct.df)[j],"\n(MAX of census tracts in county)")
  plt$ggplot_polygon <- geom_polygon(aes(fill = value),color=NA)
  p<-plt$render()
  print(p)
}
dev.off()


############ NA Census tract county map
pctdir <- "CVI-county-pct"
indicators.df<-fread("CVI_indicators_current.csv")
cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE,integer64 = "numeric")

ct.df10 <- as.data.frame(counties(year=2010))
ct.df10 <- ct.df10[,names(ct.df10)%in%c("GEOID10","INTPTLAT10","INTPTLON10")]
ct.df10$LatLong <- paste(ct.df10$INTPTLAT10,ct.df10$INTPTLON10,sep=",")
rownames(ct.df10)<-ct.df10$GEOID10
############ County - fraction NA
cvi.county.df <- cvi.df[,c(1:6)]
cvi.county.df$LatLong <- ct.df10[cvi.county.df$GEOID.County,"LatLong"]
cvi.county.df$GEOID.Tract <- ""
cvi.county.df <- cvi.county.df[!duplicated(cvi.county.df),]
for (j in 7:(ncol(cvi.df))) {
  cat(paste0(j,"...")) 
  tmpdat <- as.data.frame(cvi.df)[,c(1,2,3,4,j)]
  tmpdat[[5]] <- is.na(tmpdat[[5]])
#  if (indicators.df$`Replace NA with median`[j-6]==0) tmpdat[is.na(tmpdat[[5]]),5] <- 0
  tmp <- aggregate(. ~ STATE+County_Name+GEOID.State+GEOID.County,
                   data=tmpdat,
                   FUN = function (x) {sum(x)/length(x)})
  cat(paste0(names(cvi.df)[j],"..."))
  cvi.county.df <- left_join(cvi.county.df,tmp)
}

# Replace main data frame with county data frame
cvi.df <- cvi.county.df

cvi.dat.df <- cvi.df[,-(1:6)]

# idcols for use in GUI
idcols_gui.df <- data.table(
  `row#` = 1:nrow(cvi.df),
  Name=paste0(cvi.df$STATE,", ",cvi.df$County_Name),
  FIPS=cvi.df$GEOID.County,
  Source=cvi.df$LatLong
)
cvi.na.df<-cbind(idcols_gui.df, cvi.dat.df) # save for use by ToxPi GUI
fwrite(cvi.na.df,file.path(pctdir,"CVI-county_data_fracNA.csv"),quote=TRUE)

##

pdf(file.path("CVI-county-pct","CVI-county_data_fracNA.pdf"),height=6,width=10)
for (j in 5:(ncol(cvi.na.df))) {
  cat(j,"...")
  dat.df <- data.frame(region=as.numeric(cvi.na.df$FIPS),
                       value=cvi.na.df[[j]])
  if (length(unique(dat.df$value))==1) ncolors <- 0 else ncolors=1
  mtitle <- paste0(names(cvi.na.df)[j],"\n(Fraction of NA census tracts in county)")
  if (indicators.df$`Replace NA with median`[j-4]==0) {
    mtitle <- paste(mtitle, "[NA means 0]")
  }
  plt<-CountyChoropleth$new(dat.df)
  plt$set_num_colors(ncolors)
  plt$set_zoom(NULL)
  # plt$ggplot_scale <- scale_fill_viridis_c("County\nFraction",option="magma",limits=c(0,1))
  plt$title<-mtitle
  plt$ggplot_polygon <- geom_polygon(aes(fill = value),color=NA)
  p<-plt$render()
  print(p)
}
dev.off()

