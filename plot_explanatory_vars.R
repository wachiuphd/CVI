rm(list = ls())
library(readxl)
library(PerformanceAnalytics)
library(dplyr)
library(purrr)
library(ggplot2)

working <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working)
datarepo=paste0(working,"/Data/")
cvidatarepo=paste0(working,"/CVI-pct/")
tracts=paste0(datarepo,"tracts2010.xlsx")
nlcd=paste0(datarepo,"nlcd.xlsx")
other=paste0(datarepo,"other.xlsx")
income=paste0(datarepo,"income.csv")
baseline=paste0(cvidatarepo,"CVI-pct-comb-baseline.csv")
climate=paste0(cvidatarepo,"CVI-pct-comb-climate.csv")
overall=paste0(cvidatarepo,"CVI-pct-comb.csv")

tractsdf=read_xlsx(tracts)
nlcddf=read_xlsx(nlcd)
otherdf=read_xlsx(other)

tractsdf[,5]=as.numeric(unlist(tractsdf[,5]))
nlcddf[,1]=as.numeric(unlist(nlcddf[,1]))
otherdf[,1]=as.numeric(unlist(otherdf[,1]))

baselinedf=read.csv(baseline)
baselinedf=baselinedf[,c(3,5)]
colnames(baselinedf)<-c("GEOID10","BaselineScore")

climatedf=read.csv(climate)
climatedf=climatedf[,c(3,5)]
colnames(climatedf)<-c("GEOID10","ClimateScore")

overalldf=read.csv(overall)
overalldf=overalldf[,c(3,5)]
colnames(overalldf)<-c("GEOID10","Score")

alldomains=inner_join(overalldf,baselinedf,"GEOID10")
alldomains=inner_join(alldomains,climatedf,"GEOID10")

listraw=list(overalldf,baselinedf,climatedf)

nlcddomains=merge(alldomains, nlcddf,by="GEOID10")
otherdomains=merge(alldomains,otherdf,by="GEOID10")
otherdomains[,6]=as.numeric(otherdomains[,6])

incomedf=read.csv(income)
incomedomains=merge(alldomains,incomedf,by="GEOID10")
incomedomains[,7]=as.numeric(incomedomains[,7])


#Charts
chart.Correlation(nlcddomains[,c(2:4,5:11)],pch=".") #correlation of domains with land use (NLCD 2019)
chart.Correlation(incomedomains[,c(2:4,7)],pch=".") #correlation of top domains plus income per capita (2012-2016, ACS)
chart.Correlation(otherdomains[,c(2:4,6)],pch=".") #correlation of top domains plus population density

labs=c("","Overall CVI Score", "Baseline CVI Score", "Climate Change Risk")

#HOLC Grade charts
xvar=otherdomains$HOLC_GRADE
for(i in 2:4){
  yvar=otherdomains[,i]
  
  g <- ggplot(otherdomains, aes(x=xvar, y=yvar, group=xvar, fill=as.factor(xvar))) + geom_violin(trim=FALSE) + stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="gray")+geom_text(aes(label=format(..count..,big.mark = ",",trim = TRUE)), y=0.10, stat='count', colour="black", size=4)+geom_hline(yintercept=mean(yvar),linetype="dashed") + theme(legend.position="none",axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"))+xlab("HOLC Grade")+ylab(labs[i])+ylim(0,1)
  
  plot(g)
  # ggsave(filename = paste0(working,"/Exports/HOLC_",labs[i],".pdf"),dpi=300)  #uncomment to export PDFs of plots
}

#EPA region charts
xvar=otherdomains$EPA_REG

for(i in 2:4){
  yvar=otherdomains[,i]
  g <- ggplot(otherdomains, aes(x=xvar, y=yvar, group=xvar, fill=as.factor(xvar))) + geom_violin(trim=FALSE) + scale_x_continuous(breaks=seq(1,10,1)) + stat_summary(fun=mean, geom="point", shape=20, size=5, color="black", fill="gray")+ geom_hline(yintercept=mean(yvar),linetype="dashed")+geom_text(aes(label=format(..count..,big.mark = ",",trim = TRUE)), y=0.10, stat='count', colour="black",size=3)+theme(legend.position="none",axis.text=element_text(size=10),axis.title=element_text(size=12,face="bold"))+xlab("EPA Region")+ylab(labs[i])+ylim(0,1)
  
  plot(g)
  # ggsave(filename = paste0(working,"/Exports/EPAReg_",labs[i],".pdf"),dpi=300)  #uncomment to export PDFs of plots
}