library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

pctdir <- "CVI-pct"
datafolder <- "Data"
supfigfolder <- "SuppFigures"

cvi.toxpi.df <- fread(file.path(pctdir,"CVI-pct-comb.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
names(cvi.toxpi.df)[5] <- "CVI ToxPi Score"
cvi.toxpi.baseline.df <- fread(file.path(pctdir,"CVI-pct-comb-baseline.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
names(cvi.toxpi.baseline.df)[5] <- "Baseline ToxPi Score"
cvi.toxpi.climate.df <- fread(file.path(pctdir,"CVI-pct-comb-climate.csv"),integer64 = "double",
                               keepLeadingZeros = TRUE)
names(cvi.toxpi.climate.df)[5] <- "Climate Change ToxPi Score"

cejst.dat <- fread(file.path(datafolder,"CEJST.csv"),integer64 = "double",
                   keepLeadingZeros = TRUE)
names(cejst.dat)[1]<-"FIPS"
ntracts <- nrow(cvi.toxpi.df)
cvi.cejst <- left_join(cvi.toxpi.df,cvi.toxpi.baseline.df)
cvi.cejst <- left_join(cvi.cejst,cvi.toxpi.climate.df)
cvi.cejst <- left_join(cvi.cejst,cejst.dat[,c("FIPS","Identified as disadvantaged")])
cvi.cejst$CEJST <- "Not disadvantaged"
cvi.cejst$CEJST[cvi.cejst$`Identified as disadvantaged`] <- "Disadvantaged"

cvi.cejst.df <- pivot_longer(cvi.cejst,cols=5:14)

cvi.cejst.df$name <- factor(cvi.cejst.df$name,levels=
                              c("CVI ToxPi Score",
                                "Baseline ToxPi Score",
                                "Baseline.Health",
                                "Baseline.SocialEconomic",
                                "Baseline.Infrastructure",
                                "Baseline.Environment",
                                "Climate Change ToxPi Score",
                                "ClimateChange.Health",
                                "ClimateChange.SocialEconomic",
                                "ClimateChange.ExtremeEvents"))

p<-ggplot(subset(cvi.cejst.df,name %in% c("CVI ToxPi Score","Baseline ToxPi Score","Climate Change ToxPi Score")))+
  geom_boxplot(aes(x=value,y=CEJST,
                   fill=CEJST))+
  theme_bw()+facet_wrap(~name,ncol=1)+theme(axis.title.y = element_blank(),
                                            axis.text.y = element_blank())+
  xlim(0,1)+
  scale_fill_discrete("CEJST",breaks=c("Not disadvantaged","Disadvantaged"))

ggsave(file.path(supfigfolder,"CVI_CEJST_boxplot.pdf"),p,height=3.4,width=5.6)

indxorder <- order(cvi.cejst$`CVI ToxPi Score`,decreasing = TRUE)

cat("top 100 in CEJST disadvantaged:",100*sum(cvi.cejst[indxorder[1:100],c("Identified as disadvantaged")])/100,"%\n")
cat("top 1% in CEJST disadvantaged:",100*sum(cvi.cejst[indxorder[1:round(ntracts/100)],c("Identified as disadvantaged")])/round(ntracts/100),"%\n")
cat("top 10% in CEJST disadvantaged:",100*sum(cvi.cejst[indxorder[1:round(ntracts/10)],c("Identified as disadvantaged")])/round(ntracts/10),"%\n")

