library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(choroplethr)
library(choroplethrMaps)
library(tidyr)
library(ggpubr)

data(state.regions)
state.regions$`Census Region` <- "South"
state.regions$`Census Region`[state.regions$abb %in% 
                             c("CT","ME","MA","NH","RI","VT",
                               "NY","NJ","PA")] <- "Northeast"
state.regions$`Census Region`[state.regions$abb %in% 
                             c("IL","IN","MI","OH","WI",
                               "IA","KS","MN","MO","NE",
                               "ND","SD")] <- "Midwest"
state.regions$`Census Region`[state.regions$abb %in% 
                             c("AZ","CO","ID","MT","NV","NM",
                               "UT","WY","AK","CA","HI","OR",
                               "WA")] <- "West"
rownames(state.regions)<-state.regions$abb

cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE,integer64 = "numeric")
cvi.toxpi.df <- fread(file.path("CVI-pct","CVI-pct-comb.csv"),
                      integer64 = "double",keepLeadingZeros = TRUE)

scores.df <- cvi.toxpi.df[,c(3,5:12)]
names(scores.df)[1] <- "GEOID.Tract"
scores.df <- left_join(scores.df,cvi.df[,1:5])
scores.df$state <- factor(scores.df$STATE)
scores.df$county <- factor(paste(scores.df$County_Name,scores.df$STATE,sep=", "))
scores.df$`Census Region` <- state.regions[scores.df$STATE,"Census Region"]
scores.df$`Census Region` <- factor(scores.df$`Census Region`,levels=
                                      c("Northeast","Midwest","South","West"))

scores.df.med <- aggregate(.~state,FUN=median,data=scores.df[,c(2:9,14)])
scores.df.max <- aggregate(.~state,FUN=max,data=scores.df[,c(2:9,14)])

names(scores.df)[2:9] <- c(
  "CVI Score",
  "Baseline Vulnerability: Health",
  "Baseline Vulnerability: Social and Economic",
  "Baseline Vulnerability: Infrastructure",
  "Baseline Vulnerability: Enviroment",
  "Climate Change Risk: Health",
  "Climate Change Risk: Social and Economic",
  "Climate Change Risk: Extreme Events"
)


############# County level maps

scores.df.med.county <- aggregate(.~GEOID.County,
                                  FUN=median,data=scores.df[,c(2:9,13)])
mapplt.list <-list()

scores.df.mean.county <- aggregate(.~GEOID.County,
                                  FUN=mean,data=scores.df[,c(2:9,13)])

scores.df.var.county <- aggregate(.~GEOID.County,
                                   FUN=function(x) {if (length(x)>1) var(x) else 0},data=scores.df[,c(2:9,13)])

scores.df.relvar.county <- cbind(GEOID.County=scores.df.var.county[,1],
                                 as.data.frame(sweep(data.matrix(scores.df.var.county[,2:9]),
                                       2,base::apply(as.matrix(scores.df[,2:9]),2,FUN=var),
                                       FUN = "/")))

scores.df.max.county <- aggregate(.~GEOID.County,
                                   FUN=max,data=scores.df[,c(2:9,13)])

viridisopts <- c("A","B","C","D","E","A","B","C")

for (j in 1:8) {
  dat.df <- data.frame(region=as.numeric(scores.df.med.county$GEOID.County),
                       value=scores.df.med.county[[j+1]])
  plt<-CountyChoropleth$new(dat.df)
  plt$set_num_colors(1)
  plt$set_zoom(NULL)
  plt$ggplot_scale <- scale_fill_viridis_c("",option=viridisopts[j],limits=c(0,1))
  plt$title<-paste(letters[j],names(scores.df.med.county)[j+1])
  plt$ggplot_polygon <- geom_polygon(aes(fill = value),color=NA)
  mapplt.list[[j]]<-plt$render()
}
figmaps <- ggarrange(plotlist=mapplt.list,nrow=4,ncol=2)
ggsave("CVI_maps.pdf",figmaps,height=7,width=6.5,scale=2)


for (j in 1:8) {
  dat.df <- data.frame(region=as.numeric(scores.df.max.county$GEOID.County),
                       value=scores.df.max.county[[j+1]])
  plt<-CountyChoropleth$new(dat.df)
  plt$set_num_colors(1)
  plt$set_zoom(NULL)
  plt$ggplot_scale <- scale_fill_viridis_c("",option=viridisopts[j],limits=c(0,1))
  plt$title<-paste(letters[j],names(scores.df.max.county)[j+1])
  plt$ggplot_polygon <- geom_polygon(aes(fill = value),color=NA)
  mapplt.list[[j]]<-plt$render()
}
figmaps <- ggarrange(plotlist=mapplt.list,nrow=4,ncol=2)
ggsave("CVI_maps_max.pdf",figmaps,height=7,width=6.5,scale=2)

for (j in 1:8) {
  dat.df <- data.frame(region=as.numeric(scores.df.relvar.county$GEOID.County),
                       value=scores.df.relvar.county[[j+1]])
  plt<-CountyChoropleth$new(dat.df)
  plt$set_num_colors(1)
  plt$set_zoom(NULL)
  plt$ggplot_scale <- scale_fill_viridis_c("",option=viridisopts[j])
  plt$title<-paste(letters[j],names(scores.df.relvar.county)[j+1])
  plt$ggplot_polygon <- geom_polygon(aes(fill = value),color=NA)
  mapplt.list[[j]]<-plt$render()
}
figmaps <- ggarrange(plotlist=mapplt.list,nrow=4,ncol=2)
ggsave("CVI_maps_relvar.pdf",figmaps,height=7,width=6.5,scale=2)

############# Source of variability

scores.df.long <- pivot_longer(scores.df,cols = 2:9)
scores.df.long$name <- factor(scores.df.long$name,
                              levels=unique(scores.df.long$name)
)
scores.df.long$state <- factor(scores.df.long$state,
                               levels=as.character(
                                 scores.df.med$state[order(scores.df.med$`ToxPi Score`)]))
scores.df.long$label <- as.character(gsub(": ",":\n",scores.df.long$name))
scores.df.long$label <- factor(scores.df.long$label,levels=
                                 as.character(gsub(": ",":\n",levels(scores.df.long$name))))

tmp.topcat <- 
  aggregate(value ~ GEOID.Tract,
            data = subset(scores.df.long,name != "CVI Score"),
            FUN = which.max)
tmp.topcat$TopCat <- factor(levels(scores.df.long$name)[tmp.topcat$value+1],
                            levels=rev(levels(scores.df.long$name)))
scores.df <- left_join(scores.df,tmp.topcat[,c(1,3)])

tmp.topcat <- 
  aggregate(value ~ GEOID.Tract,
            data = subset(scores.df.long,name %in% c(
              "Baseline Vulnerability: Health",
              "Baseline Vulnerability: Social and Economic",
              "Baseline Vulnerability: Infrastructure",
              "Baseline Vulnerability: Enviroment"
            )),
            FUN = which.max)
tmp.topcat$TopCatVulnerability <- factor(levels(scores.df.long$name)[tmp.topcat$value+1],
                                         levels=rev(levels(scores.df.long$name)))

scores.df <- left_join(scores.df,tmp.topcat[,c(1,3)])

tmp.topcat <- 
  aggregate(value ~ GEOID.Tract,
            data = subset(scores.df.long,name %in% c(
              "Climate Change Risk: Health",
              "Climate Change Risk: Social and Economic",
              "Climate Change Risk: Extreme Events"
            )),
            FUN = which.max)
tmp.topcat$TopCatRisk <- factor(levels(scores.df.long$name)[tmp.topcat$value+5],
                                levels=rev(levels(scores.df.long$name)))
scores.df <- left_join(scores.df,tmp.topcat[,c(1,3)])

fwrite(scores.df,"Scores_summary.csv")

##### Summary boxplots

scoresboxplt<-
  ggplot(scores.df.long)+geom_boxplot(aes(x=value,y=state,fill=`Census Region`),outlier.size = 0.2)+
  scale_fill_viridis_d(begin=0.3)+
  facet_wrap(~label,nrow=1)+xlab("")+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
print(scoresboxplt)
ggsave("StateSummaryBoxplots.pdf",scoresboxplt,height=3,width=6.5,scale=2)

#####

indicators.df <- fread("CVI_indicators_current.csv")
indicators.df$`Baseline Vulnerability` <- 
  factor(indicators.df$`Baseline Vulnerability`,
         levels=unique(indicators.df$`Baseline Vulnerability`))

indicators.geo <- 
  rbind(
    as.matrix(t(table(indicators.df$GeographicScale))),
    as.matrix(ftable(indicators.df$`Baseline Vulnerability`,
                     indicators.df$GeographicScale)))
rownames(indicators.geo)<- c(
  "CVI Score",
  "Baseline Vulnerability: Health",
  "Baseline Vulnerability: Social and Economic",
  "Baseline Vulnerability: Infrastructure",
  "Baseline Vulnerability: Enviroment",
  "Climate Change Risk: Health",
  "Climate Change Risk: Social and Economic",
  "Climate Change Risk: Extreme Events"
)
indicators.geo <- as.data.frame(indicators.geo)
indicators.geo$label <- gsub(": ",":\n",rownames(indicators.geo))
indicators.geo$label <- factor(indicators.geo$label,levels=
                                 rev(indicators.geo$label))
indicators.geo.df <- pivot_longer(indicators.geo,1:4)
indicators.geo.df$name <- factor(indicators.geo.df$name,
                                 levels=rev(c("State","County","Tract","Tract (raster)")))

pgeo <- ggplot(indicators.geo.df)+
  geom_col(aes(x=label,y=value,fill=name),position="fill")+
  scale_y_continuous(label = scales::percent)+
  scale_fill_viridis_d(begin=0,end=0.9,option="C")+
  labs(x="",y="Geographic Scale (% of indicators)")+theme_bw()+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.text.y=element_text(hjust=0))+
  guides(fill = guide_legend(reverse = TRUE))
print(pgeo)
ggsave("GeoScale.pdf",pgeo,height=3,width=2.33,scale=2)


##### Heterogeneity - state and county level
scores.df.state.SSres <- aggregate(.~STATE,data=scores.df[,2:10],
                                 FUN=function(x) {var(x)*length(x)})
state.r2 <- 1 - (base::apply(scores.df.state.SSres[,-1],2,sum))/
  (base::apply(scores.df[,2:9],2,var)*nrow(scores.df))

scores.df.county.SSres <- aggregate(.~GEOID.County,data=scores.df[,c(2:9,13)],
                                   FUN=function(x) {if (length(x)>1) (var(x)*length(x)) else 0})
county.r2 <- 1 - (base::apply(scores.df.county.SSres[,-1],2,sum))/
  (base::apply(scores.df[,2:9],2,var)*nrow(scores.df))

r2.df <- rbind(data.frame(Variable=rep("State",8),
                          Category=gsub(": ",":\n",names(state.r2)),
                          R2=state.r2),
               data.frame(Variable=rep("County",8),
                          Category=gsub(": ",":\n",names(county.r2)),
                          R2=county.r2)
)
r2.df$Category<-factor(r2.df$Category,levels=
                         rev(gsub(": ",":\n",names(state.r2))))
pr2<-ggplot(r2.df)+
  geom_col(aes(x=Category,y=R2,
               group=Variable,fill=Variable),position="dodge")+
  scale_y_continuous(label = scales::percent)+
  scale_fill_viridis_d(begin=0.6,end=0.9,option="C")+
  labs(x="",y=bquote('% census tract variance due to state or county '(R^2)))+
  coord_flip()+theme_bw()+theme(legend.title = element_blank(),
                                axis.text.y=element_text(hjust=0))+
  guides(fill = guide_legend(reverse = TRUE))

print(pr2)
ggsave("R2.pdf",pr2,height=3,width=2.33,scale=2)

###### Bar graph of top categories

scores.df.cat <- cbind(data.frame(Sample=rep("All",8)),
                       as.data.frame(prop.table(table(scores.df$TopCat))))
scores.df.cat <- rbind(scores.df.cat,
                       cbind(data.frame(Sample=rep("Top quartile",8)),
                             as.data.frame(prop.table(table(
                               subset(scores.df,`CVI Score` >= quantile(scores.df$`CVI Score`,prob=0.75))$TopCat)))
                       )
)
scores.df.cat <- rbind(scores.df.cat,
                       cbind(data.frame(Sample=rep("Top decile",8)),
                             as.data.frame(prop.table(table(
                               subset(scores.df,`CVI Score` >= quantile(scores.df$`CVI Score`,prob=0.9))$TopCat)))
                       )
)
names(scores.df.cat)[3] <- "Fraction of census tracts"
scores.df.cat$`Dominant Category` <-
  as.character(gsub(": ",":\n",scores.df.cat$Var1))
scores.df.cat$`Dominant Category` <- factor(scores.df.cat$`Dominant Category`,
                                            levels=
                                              as.character(gsub(": ",":\n",levels(scores.df.cat$Var1))))
scores.df.cat <- subset(scores.df.cat,`Dominant Category` != "CVI Score")
scores.df.cat$Sample <- factor(scores.df.cat$Sample,
                               levels=c("All","Top quartile","Top decile"))
pcat <- ggplot(scores.df.cat)+
  geom_col(aes(x=`Dominant Category`,y=`Fraction of census tracts`,
               fill=Sample,group=Sample),position="dodge")+
  scale_y_continuous(label = scales::percent)+
  scale_fill_viridis_d(begin=0.2,end=0.8,option="inferno")+xlab("")+
  labs(subtitle = "Dominant Category",y="% of census tracts")+
  coord_flip()+theme_bw()+theme(legend.title = element_blank(),
                                axis.text.y=element_text(hjust=0))+
  guides(fill = guide_legend(reverse = TRUE))
print(pcat)
ggsave("Top.Categories.pdf",pcat,height=3,width=2.33,scale=2)
###

pfig <- ggarrange(scoresboxplt,
                  ggarrange(pgeo,pr2,pcat,labels=c("b","c","d"),nrow=1),
                  labels=c("a",""),ncol=1,heights=c(2,1))
ggsave("State_County_summary_fig.pdf",pfig,height=4.5,width=6.5,scale=2)

### Not run
# 
# toxpi.r2 <- data.frame(name=names(scores.df)[2:9])
# toxpi.r2$R2.state <- NA
# toxpi.r2$R2.county <- NA
# pdf("R2_diagnostics_state.pdf")
# par(mfrow=c(2,2))
# for (j in 1:8) {
#   print(j)
#   cols <- c(12,13,j+1)
#   tmp.df <- scores.df[,..cols]
#   names(tmp.df)[3] <- "y"
#   lm.state <- lm(y ~ GEOID.State,data=tmp.df)
#   plot(lm.state,pch=15,cex=0.1)
#   toxpi.r2$R2.state[j]<-summary(lm.state)$adj.r.squared
# }
# dev.off()
# fwrite(toxpi.r2,"R2state.county.csv")
# 
# ## County
# pdf("R2_diagnostics_county.pdf")
# par(mfrow=c(2,2))
# for (j in 1:8) {
#   print(j)
#   cols <- c(12,13,j+1)
#   tmp.df <- scores.df[,..cols]
#   names(tmp.df)[3] <- "y"
#   lm.county <- lm(y ~ GEOID.County,data=tmp.df)
#   plot(lm.county,pch=15,cex=0.1)
#   toxpi.r2$R2.county[j]<-summary(lm.county)$adj.r.squared
# }
# dev.off()
# fwrite(toxpi.r2,"R2state.county.csv")

###########
# 
# tmp.df.topdecile <- data.frame(state=state.regions$abb)
# pdf("State Summaries.pdf",height=3.5,width=6)
# for (j in 1:8) {
#   ytitle <- names(scores.df)[j+1]
#   cols <- c(14,j+1)
#   tmp.df <- scores.df[,..cols]
#   cols.m <- c(1,j+1)
#   tmp.df.med <- scores.df.med[,cols.m]
#   names(tmp.df)[2] <- "y"
#   names(tmp.df.med)[2] <- "y"
#   tmp.df$state <- factor(tmp.df$state,
#                          levels=as.character(
#                            tmp.df.med$state[order(tmp.df.med$y)]))
#   r2state <- toxpi.r2$R2.state[j]
#   plt <- ggplot(tmp.df)+geom_boxplot(aes(y=y,x=state),outlier.size = 0.5)+
#     ylab(ytitle)+theme_bw()+ggtitle(ytitle)+ylim(0,1)+
#     annotate("text",x=1,y=1,
#                   label=paste0("R2.state=",signif(r2state,2)),
#              vjust=1,hjust=0)+
#     theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
#   print(plt)
#   
#   tmp.df.max <- scores.df.max[,cols.m]
#   names(tmp.df.max)[2] <- "y"
#   tmp.df.max$y75 <- tmp.df.max$y >= quantile(tmp.df$y,prob=0.75)
#   numstates<-sum(tmp.df.max$y75)
#   cat(ytitle,": ",numstates,"/51 states with census tracts > 75th %ile\n",
#       sep="")
#   tmp.df.max$region <- state.regions[as.character(tmp.df.max$state),"region"]
#   tmp.df.max$value <- as.numeric(tmp.df.max$y75)
#   pst75 <- StateChoropleth$new(tmp.df.max)
#   pst75$title = paste0(ytitle,"\nStates with census tracts > 75th %ile")
#   pst75$set_num_colors(1)
#   pst75$set_zoom(NULL)
#   pst75$show_labels = FALSE
#   print(pst75$render())
#   
#   tmp.df.max$y90 <- tmp.df.max$y >= quantile(tmp.df$y,prob=0.90)
#   numstates<-sum(tmp.df.max$y90)
#   cat(ytitle,": ",numstates,"/51 states with census tracts > 90th %ile\n",
#       sep="")
#   tmp.df.max$value <- as.numeric(tmp.df.max$y90)
#   pst90 <- StateChoropleth$new(tmp.df.max)
#   pst90$title = paste0(ytitle,"\nStates with census tracts > 90th %ile")
#   pst90$set_num_colors(1)
#   pst90$set_zoom(NULL)
#   pst90$show_labels = FALSE
#   print(pst90$render())
#   
#   names(tmp.df.max)[6] <- ytitle
#   tmp.df.max$state <- as.character(tmp.df.max$state)
#   tmp.df.topdecile <- left_join(tmp.df.topdecile,tmp.df.max[,c(1,6)])
# }
# dev.off()
# 
# tmp.df.topdecile$Baseline <- tmp.df.topdecile$Baseline.Health+tmp.df.topdecile$Baseline.SocialEconomic+tmp.df.topdecile$Baseline.Infrastructure+tmp.df.topdecile$Baseline.Environment > 0
# tmp.df.topdecile$ClimateChange <- tmp.df.topdecile$ClimateChange.Health+tmp.df.topdecile$ClimateChange.SocialEconomic+tmp.df.topdecile$ClimateChange.ExtremeEvents > 0



# library(dbscan)
# library(GGally)
# set.seed(3.14159)
# samp.indx<-sample.int(nrow(scores.df),7000)
# ggpairs(scores.df[samp.indx,3:9],
#         lower = list(continuous = wrap("points", alpha=0.2,size=0.1)))
# 
# TX.indx <- scores.df$STATE=="TX" 
# ggpairs(scores.df[TX.indx,3:9],
#         lower = list(continuous = wrap("points", alpha=0.2,size=0.1)))
# 
# 
# harris.indx <- scores.df$STATE=="TX" & scores.df$County_Name=="Harris"
# ggpairs(scores.df[harris.indx,3:9],
#         lower = list(continuous = wrap("points", alpha=0.2,size=0.1)))
# 
# 
# scores.mat <- as.matrix(scores.df[samp.indx,3:9])
# h <- hdbscan(scores.mat,minPts=100)
