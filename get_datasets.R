library(data.table)
library(stringr)
library(dplyr)
library(naniar)
library(readxl)
datafolder <- "Data"

# Census tracts from 22-02_CVI_state_county_tract.xlsx
tractsraw <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/Other/22-02_CVI_state_county_tract.xlsx",
  sheet="Tract")
tracts <- tractsraw[,c("STATE","County_Name")]
tracts$GEOID.State <- tractsraw$STATEFP10
tracts$GEOID.County <- tractsraw$FIPS
tracts$GEOID.Tract <- tractsraw$GEOID10
tracts <- tracts[order(tracts$GEOID.Tract),]

# ### Get census tracts from 2019 Census Tract Gazetteer File - has lat long
# ### https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2019_Gazetteer/2019_Gaz_tracts_national.zip
tracts_latlong <- fread(file.path(datafolder,"2019_Gaz_tracts_national.txt"),
                        keepLeadingZeros = TRUE)
tracts_latlong[, c("USPS","ALAND","AWATER","ALAND_SQMI","AWATER_SQMI"):=NULL]
setnames(tracts_latlong,"GEOID","GEOID.Tract")
tracts_latlong$LatLong <- paste0(tracts_latlong$INTPTLAT,",",tracts_latlong$INTPTLONG)
tracts_latlong[, c("INTPTLAT","INTPTLONG"):=NULL]

# Join lat long to tracts
tracts <- left_join(tracts,tracts_latlong)

### Get master sheet

cvi.master <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/CurrentCVIIndicatorsDoc - 02.23.22.xlsx",
                        sheet="Subcategories (In Progress)",trim_ws = FALSE)

checkdatrow <- function(jrow) {
  j<-jrow
  fname<-(file.path("~/Dropbox/Climate Health Vulnerability Index",
                    cvi.master$`Path To File`[j]))
  fname.exists <- file.exists(fname)
  cat(fname," exists? ",fname.exists,"\n")
  tmp <- data.frame()
  if (fname.exists) {
    tmp <- fread(fname,
                 keepLeadingZeros = TRUE)
    print(names(tmp))
    cat("****GEOID Column Name ",cvi.master$`GEOID Column Name`[j],
        cvi.master$`GEOID Column Name`[j] %in% names(tmp),"\n")
    if (!is.na(cvi.master$`Subset Column Name`[j])) {
      subsetcolnum <- which(names(tmp)==cvi.master$`Subset Column Name`[j])
      tmp.subsets<-unique(tmp[[subsetcolnum]])
      cvi.master$`Indicator Name`[j] %in% tmp.subsets
      cat("****Indicator subset ",cvi.master$`Indicator Name`[j],
          cvi.master$`Indicator Name`[j] %in% tmp.subsets,"\n")
      tmp <- tmp[tmp[[subsetcolnum]]==cvi.master$`Indicator Name`[j],]
    }
    cat("****Data Column Name ",cvi.master$`Data Column Name`[j],
        cvi.master$`Data Column Name`[j] %in% names(tmp),"\n")
    cols<-c("GEOID",
            cvi.master$`Indicator Name`[j]
            )
    setnames(tmp, cvi.master$`GEOID Column Name`[j], cols[1])
    setnames(tmp, cvi.master$`Data Column Name`[j], cols[2])
    tmp <- tmp[,..cols]
    if (class(tmp[[2]]) == "character") {
      tmp[[2]][tmp[[2]]=="N/A"] <- NA
      tmp[[2]][tmp[[2]]=="Missing"] <- NA
      tmp[[2]] <- gsub("%","",tmp[[2]]) # get rid of % sign
      tmp[[2]] <- as.numeric(tmp[[2]])
    }
    tmp[[2]][tmp[[2]]==-999] <- NA
    tmp <- subset(tmp,!is.na(GEOID)) # Remove row without GEOID
    tmp <- subset(tmp,GEOID != "") # Remove row without GEOID
    cat("number of rows: ",nrow(tmp),"\n")
    cat("number of na rows: ",sum(is.na(tmp[,2]))," or ",sum(tmp[,2]<0),"\n")
#    tmp$Category <- cvi.master$`Baseline Vulnerability`[j]
#    tmp$Subcategory <- cvi.master$Subcategory[j]
    print(head(tmp))
    cat("-------------------------------\n\n")
  }
  if (nrow(tmp)>0 & ncol(tmp)==2 & sum(!is.na(tmp[[2]])>0)) {
    return(tmp)
  } else {
    return(NA)
  }
}

pdf("CheckDist.pdf")
options(width=200)
capture.output( {
for (j in 1:nrow(cvi.master)) {
  print(paste("-----",j))
  tmp.df <- try(checkdatrow(j))
  if (is.data.frame(tmp.df)) {
    tmp.df$GEOID <- as.character(tmp.df$GEOID)
    tmp.df <- unique(tmp.df) # removed duplicated rows
    if (length(unique(str_length(tmp.df$GEOID)))==1) {
      if (str_length(tmp.df$GEOID[1]) == 12) { 
        # Block Groups - convert to median of tract
        tmp.df$GEOID.Tract <- substr(tmp.df$GEOID,1,11)
        catname <- names(tmp.df)[2]
        names(tmp.df)[2] <- "value"
        tmp.df <- tmp.df[,.(y = median(value,na.rm=TRUE)), by = GEOID.Tract]
        names(tmp.df)[2] <- catname  
        tmp.df$GEOID.County <- substr(tmp.df$GEOID.Tract,1,5)
        tmp.df$GEOID.State <- substr(tmp.df$GEOID.Tract,1,2)
        tracts <- left_join(tracts,tmp.df)
      } else if (str_length(tmp.df$GEOID[1]) == 11) {
        tmp.df$GEOID.County <- substr(tmp.df$GEOID,1,5)
        tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
        names(tmp.df)[1] <- "GEOID.Tract"
        tracts <- left_join(tracts,tmp.df)
      } else if (str_length(tmp.df$GEOID[1]) == 5) {
        tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
        names(tmp.df)[1] <- "GEOID.County"
        tracts <- left_join(tracts,tmp.df)
      } else if (str_length(tmp.df$GEOID[1]) == 2) {
        names(tmp.df)[1] <- "GEOID.State"
        tracts <- left_join(tracts,tmp.df)
      } 
    } else if (length(unique(str_length(tmp.df$GEOID)))==2) {
      geoidlengths <- sort(unique(str_length(tmp.df$GEOID)))
      if (sum(geoidlengths==c(10,11)) == 2) { 
        # tracts, missing leading zeros
        indx10 <- str_length(tmp.df$GEOID)==10
        tmp.df$GEOID[indx10] <- paste0("0",tmp.df$GEOID[indx10])
        tmp.df$GEOID.County <- substr(tmp.df$GEOID,1,5)
        tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
        names(tmp.df)[1] <- "GEOID.Tract"
        tracts <- left_join(tracts,tmp.df)
      } else if (sum(geoidlengths==c(4,5)) == 2) {
        # counties, missing leading zeros
        indx4 <- str_length(tmp.df$GEOID)==4
        tmp.df$GEOID[indx4] <- paste0("0",tmp.df$GEOID[indx4])
        tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
        names(tmp.df)[1] <- "GEOID.County"
        tracts <- left_join(tracts,tmp.df)
      } else if (sum(geoidlengths==c(1,2)) == 2) {
        # states, missing leading zeros
        indx1 <- str_length(tmp.df$GEOID)==1
        tmp.df$GEOID[indx1] <- paste0("0",tmp.df$GEOID[indx1])
        names(tmp.df)[1] <- "GEOID.State"
        tracts <- left_join(tracts,tmp.df)
      }
    }
    print(head(tmp.df))
  }
  if ("GEOID.State" %in% names(tmp.df)) {
    print(paste("******",j,"Verified ******"))
    vartext <- abbreviate(names(tmp.df)[2],minlength = 20)
    par(mfrow=c(2,2))
    y <- tmp.df[[2]]
    y <- y[!is.na(y)]
    if (length(y)>0) {
      if (length(grep("change",cvi.master$`Data Column Name`[j],ignore.case = TRUE)) < 1 &
          length(grep("additional",cvi.master$`Data Column Name`[j],ignore.case = TRUE)) < 1 &
          cvi.master$`Data Column Name`[j] != "CCI_EE_FDmean_days") {
        y <- y[y>=0]
      }
      if (max(y)<=0) {
        vartext <- paste0("-",vartext)
        y <- -y
      }
      qqnorm(y,main="",pch=15,cex=0.2); qqline(y);
      qqnorm(log(y[y>0]),main="",pch=15,cex=0.2); qqline(log(y[y>0]));
      hist(y,main="",xlab=vartext);
      hist(log(y[y>0]),main="",xlab=paste("Log",vartext));
      mtext(paste("Row",j,"\n",names(tmp.df)[2]),outer=TRUE,line=-2,cex=0.75)
    }
  } else {
    print(paste("!!!!!!",j,"Not processed !!!!!!"))
  }
  cat("\n\n")
}
},file="Checkoutput.txt")
options(width=80)
dev.off()
icols <- c("Indicator Name","Adverse Direction","Replace NA with median","Baseline Vulnerability ","Subcategory","Parameters","Agency or data source","Year of data release","Geographic Level")
indicators.df <- data.table(`Indicator Name`=names(tracts)[-(1:6)])
indicators.df <- left_join(indicators.df,as.data.table(cvi.master)[,..icols])
indicators.df$`Adverse Direction`<-as.numeric(indicators.df$`Adverse Direction`)
# replace "n/a" with 0
indicators.df$`Replace NA with median`<-as.numeric(indicators.df$`Replace NA with median`)
indicators.df$`Replace NA with median`[is.na(indicators.df$`Replace NA with median`)]<-0
fwrite(tracts,"CVI_data_current.csv")
fwrite(indicators.df,"CVI_indicators_current.csv")
print(as.numeric((apply(tracts,2,FUN=function(x) {sum(!is.na(x))}))))
