library(data.table)
library(stringr)
library(dplyr)
library(naniar)
library(readxl)
datafolder <- "Data"

# Census tracts from 22-02_CVI_state_county_tract_updated.xlsx
# Nine census tracts internal points moved to be internal to tract boundary
tractsraw <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/Other/22-02_CVI_state_county_tract_updated.xlsx",
  sheet="Tract")
tracts <- tractsraw[,c("STATE","County_Name")]
tracts$GEOID.State <- tractsraw$STATEFP10
tracts$GEOID.County <- tractsraw$FIPS
tracts$GEOID.Tract <- tractsraw$GEOID10
tracts$LatLong <- paste0(tractsraw$INTPTLAT10,",",tractsraw$INTPTLON10)
tracts <- tracts[order(tracts$GEOID.Tract),]

### Get master sheet

cvi.master <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/CVI Indicators_ForBeta_052222xlsx.xlsx",
                        sheet="Alpha Indicators",trim_ws = FALSE)
indicator.verified <- rep(FALSE,nrow(cvi.master))
indicator.geo <- rep("",nrow(cvi.master))

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
    ## Special case treat 999 as zero
    if (cvi.master$`Data Column Name`[j] %in% c("WalkScore","TransitScore","BikeScore")) {
      tmp[[2]][tmp[[2]]==999] <- 0
    }
    tmp <- subset(tmp,!is.na(GEOID)) # Remove row without GEOID
    tmp <- subset(tmp,GEOID != "") # Remove row without GEOID
    cat("number of rows: ",nrow(tmp),"\n")
    cat("number of na rows: ",sum(is.na(tmp[,2]))," or ",sum(tmp[,2]<0),"\n")
    names(tmp)[2] <- cvi.master$Parameters[j]
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

options(width=200)
capture.output( {
for (j in 1:nrow(cvi.master)) {
  print(paste("-----",j))
  tmp.df <- try(checkdatrow(j))
  if (is.data.frame(tmp.df)) {
    tmp.df$GEOID <- as.character(tmp.df$GEOID)
    tmp.df <- unique(tmp.df) # removed duplicated rows
    if (sum(duplicated(tmp.df$GEOID))==0) { # No duplicated GEOIDs
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
          indicator.geo[j] <- "Tract"
        } else if (str_length(tmp.df$GEOID[1]) == 11) {
          tmp.df$GEOID.County <- substr(tmp.df$GEOID,1,5)
          tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
          names(tmp.df)[1] <- "GEOID.Tract"
          tracts <- left_join(tracts,tmp.df)
          indicator.geo[j] <- "Tract"
        } else if (str_length(tmp.df$GEOID[1]) == 5) {
          tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
          names(tmp.df)[1] <- "GEOID.County"
          tracts <- left_join(tracts,tmp.df)
          indicator.geo[j] <- "County"
        } else if (str_length(tmp.df$GEOID[1]) == 2) {
          names(tmp.df)[1] <- "GEOID.State"
          tracts <- left_join(tracts,tmp.df)
          indicator.geo[j] <- "State"
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
          indicator.geo[j] <- "Tract"
        } else if (sum(geoidlengths==c(4,5)) == 2) {
          # counties, missing leading zeros
          indx4 <- str_length(tmp.df$GEOID)==4
          tmp.df$GEOID[indx4] <- paste0("0",tmp.df$GEOID[indx4])
          tmp.df$GEOID.State <- substr(tmp.df$GEOID,1,2)
          names(tmp.df)[1] <- "GEOID.County"
          tracts <- left_join(tracts,tmp.df)
          indicator.geo[j] <- "County"
        } else if (sum(geoidlengths==c(1,2)) == 2) {
          # states, missing leading zeros
          indx1 <- str_length(tmp.df$GEOID)==1
          tmp.df$GEOID[indx1] <- paste0("0",tmp.df$GEOID[indx1])
          names(tmp.df)[1] <- "GEOID.State"
          tracts <- left_join(tracts,tmp.df)
          indicator.geo[j] <- "State"
        }
      }
    } else {
      print("!!!!!!!!!Non-unique GEOIDs")
    }
    print(head(tmp.df))
  }
  if ("GEOID.State" %in% names(tmp.df)) {
    print(paste("******",j,"Verified ******"))
    indicator.verified[j] <- TRUE
  } else {
    print(paste("!!!!!!",j,"Not processed !!!!!!"))
  }
  cat("\n\n")
}
},file="Checkoutput.txt")
options(width=80)
cvi.master$Verified <- indicator.verified
cvi.master$GeographicScale <- indicator.geo
fwrite(cvi.master,"CVI_master.csv")

icols <- c("Indicator Name","Adverse Direction","Replace NA with median","Baseline Vulnerability","Subcategory","Parameters","Agency or data source","Year of data release","Geographic Level","GeographicScale")
indicators.df <- as.data.table(subset(cvi.master,Verified==TRUE))[,..icols]
indicators.df$`Adverse Direction`<-as.numeric(indicators.df$`Adverse Direction`)
# replace "n/a" with 0
indicators.df$`Replace NA with median`<-as.numeric(indicators.df$`Replace NA with median`)
indicators.df$`Replace NA with median`[is.na(indicators.df$`Replace NA with median`)]<-0
fwrite(tracts,"CVI_data_current.csv")
fwrite(indicators.df,"CVI_indicators_current.csv")
print(as.numeric((base::apply(tracts,2,FUN=function(x) {sum(!is.na(x))}))))
