library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(naniar)
library(readxl)
datafolder <- "Data"
pctdir <- "CVI-pct"
# 10 colors, color-blind friendly (# removed, all lower case)
Tol_muted <- tolower(c('88CCEE', '44AA99', '117733', '332288', 'DDCC77', '999933','CC6677', '882255', 'AA4499', 'DDDDDD'))

# Census tracts from 22-02_CVI_state_county_tract_updated.xlsx
# Nine census tracts internal points moved to be internal to tract boundary
tractsraw <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/Other/22-02_CVI_state_county_tract_updated.xlsx",
                       sheet="Tract")
tractsraw$FIPS <- tractsraw$GEOID10

tractspop <- read_xlsx("~/Dropbox/Climate Health Vulnerability Index/Other/TotalPopulationbyTract2010.xlsx",
                       sheet="TotalPopulationbyTract2010")
tractspop$Population <- as.numeric(tractspop$P001001)
tractspop <- rename(tractspop,FIPS = `Formatted GEOID`)

tractsdat <- left_join(tractsraw,tractspop)
tractsdat$PopDens <- tractsdat$Population/(tractsdat$ALAND10/1e6) # Pop per sq km land

cvi.toxpi.df <- fread(file.path(pctdir,"CVI-pct-comb.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
cvi.toxpi.df <- left_join(cvi.toxpi.df,tractsdat[,c("STATE","FIPS","ALAND10","Population","PopDens","INTPTLAT10","INTPTLON10")])

x90<-base::apply(cvi.toxpi.df[,6:12],2,quantile,prob=0.9)
xge90<-scale(cvi.toxpi.df[,6:12],center=x90,scale=FALSE)
xge90.df<-data.frame(STATE=cvi.toxpi.df$STATE,
                     Baseline=base::apply(xge90[,1:4],1,max),
                     Climate=base::apply(xge90[,5:7],1,max))
xge90.df$both <- xge90.df$Baseline>=0 & xge90.df$Climate>=0
xge90.df.state<-aggregate(both~STATE,xge90.df,max)

### K-means clustering 

x <- cvi.toxpi.df[,6:12]

# Gap statistic to determine optimal number of clusters = 6
library(factoextra)
set.seed(3.14159)
indx <- sample(nrow(x),nrow(x)*0.1) # 10% sample for determining # of clusters
indx2 <- sample(nrow(x),nrow(x)*0.1) # 10% sample for plotting
xsamp <- x[indx,]
# numclus.gap <- fviz_nbclust(xsamp,kmeans,method="gap_stat",iter.max=20)
# print(numclus.gap)

## For each cluster, plot cluster center

library(gridExtra)
x.kmeans <- kmeans(x,6)
catnames <- gsub("\\.","\n",colnames(x.kmeans$centers))
centers <- x.kmeans$centers
centers.cvi <- base::apply(centers,1,mean)
clusorder <- order(base::apply(x.kmeans$centers,1,mean),decreasing=TRUE)
centers.list <- list()
centers.list[[1]] <- pieGridGrob(centers[clusorder,],
                                 labels=paste0("Cluster ",LETTERS[1:6],"\nCVI Score: ",
                                               round(centers.cvi[clusorder],2)),
                                 fills=paste0("#",Tol_muted),vp=viewport(width=0.9, height=0.9))
#for (i in 1:6) centers.list[[i]] <- pieGrob(centers[clusorder[i],],fills=paste0("#",Tol_muted),labels="Cluster")
#centers.list[[7]] <- legendGrob(catnames,pch=15,gp=gpar(cex=0.7,col=paste0("#",Tol_muted)))
#grid.arrange(grobs=centers.list,layout_matrix = matrix(c(1,4,2,5,3,6,7,7),nrow=2,ncol=4))
centers.list[[2]] <- legendGrob(catnames,pch=15,gp=gpar(col=paste0("#",Tol_muted)))
g <- arrangeGrob(grobs=centers.list,layout_matrix = matrix(c(1,1,1,1,1,1,2,2),nrow=2,ncol=4))
grid.newpage()
grid.draw(g)
ggsave("k-means cluster-centers.pdf",g,height=3,width=6,scale=1.5)

## Plot distribution of overall scores by cluster
cvi <- cvi.toxpi.df[,5]
clusletter <- LETTERS[1:6]
names(clusletter) <- clusorder
cvi$cluster <- clusletter[paste(x.kmeans$cluster)]
cvi <- cbind(cvi,x)

cvi.latlong <- cbind(cvi,data.frame(lat=as.numeric(cvi.toxpi.df[[17]]),
                                    long=as.numeric(cvi.toxpi.df[[18]])))
pdf("k-means cluster-points.pdf",height=5,width=10)
print(ggplot(cvi.latlong) + geom_point(aes(x=long,y=lat,color=cluster))+
  xlim(-124,-67)+ylim(24,50)+facet_wrap(~cluster))

print(ggplot(cvi.latlong) + geom_point(aes(x=long,y=lat,color=cluster))+
  xlim(-175,-124)+ylim(53,72)+facet_wrap(~cluster))

print(ggplot(cvi.latlong) + geom_point(aes(x=long,y=lat,color=cluster))+
  xlim(-160,-155)+ylim(18,23)+facet_wrap(~cluster))
dev.off()

pdf("k-means cluster-means.pdf",height=6,width=6)
par(mar=c(4,15,2,2),mfrow=c(3,2))
for (i in 1:6) {
  clusnow <- LETTERS[i]
  barplot(rev(base::apply(subset(cvi,cluster==clusnow)[,-2],2,mean)),
          las=1,horiz = TRUE,main=paste("Cluster",clusnow),
          xlab="Mean",xlim=c(0,1))
}
dev.off()

cvi.clusters.diff <- cbind(data.frame(cluster=LETTERS[1:6]),sweep(aggregate(.~cluster,data=cvi[,-1],mean)[,-1],2,
                                                                  base::apply(cvi[,-(1:2)],2,mean)))
cvi.clusters.diff.df <- pivot_longer(cvi.clusters.diff,cols=2:8)
ggplot(cvi.clusters.diff.df)+geom_col(aes(x=value,y=name))+facet_wrap(~cluster)
pdf("k-means cluster-variances.pdf",height=6,width=6)
par(mar=c(4,15,2,2),mfrow=c(3,2))
for (i in 1:6) {
  clusnow <- LETTERS[i]
  barplot(rev(base::apply(subset(cvi,cluster==clusnow)[,-2],2,var)),
          las=1,horiz = TRUE,main=paste("Cluster",clusnow),
          xlab="Variance",xlim=c(0,0.02))
}
dev.off()

pdf("k-means cluster-CV.pdf",height=6,width=6)
par(mar=c(4,15,2,2),mfrow=c(3,2))
for (i in 1:6) {
  clusnow <- LETTERS[i]
  barplot(rev(base::apply(subset(cvi,cluster==clusnow)[,-2],2,function(y) {sd(y)/mean(y)})),
          las=1,horiz = TRUE,main=paste("Cluster",clusnow),
          xlab="CV",xlim=c(0,0.4))
}
dev.off()

library(lsr)

cvi.df <- pivot_longer(cvi,col=c(1,3:9))
pdf("k-means cluster-dist.pdf",height=6,width=8)
print(ggplot(cvi.df)+geom_histogram(aes(value))+facet_grid(cluster~name))
print(ggplot(cvi.df)+geom_histogram(aes(value))+facet_grid(name~cluster))
dev.off()


## For each cluster, plot representative ToxPis (highest ranked, 75%, median, 25%, lowest ranked)
for (i in 1:6) {
  clusnum <- clusorder[i]
  indx.clus <- x.kmeans$cluster == clusnum
  scores.sum <- quantile(cvi$`ToxPi Score`[indx.clus],type=1,
                         prob=c(0,0.01,0.05,0.25,0.5,0.75,0.95,0.99,1))
  indx.clus.quants <- indx.clus & (cvi$`ToxPi Score` %in% scores.sum)
  indx.clus.quants <- (1:length(indx.clus.quants))[indx.clus.quants]
  indx.clus.quants <- indx.clus.quants[order(
    cvi$`ToxPi Score`[indx.clus.quants],decreasing=TRUE)]
  pi.list <- list()
  pi.list[[1]] <- pieGridGrob(as.matrix(x[indx.clus.quants,]),
                              labels=paste0(cvi.toxpi.df$Name[indx.clus.quants],
                                            "\nCVI Score: ",round(cvi$`ToxPi Score`[indx.clus.quants],2)),
                              fills=paste0("#",Tol_muted),vp=viewport(width=0.9, height=0.9))
  pi.list[[2]] <- legendGrob(c(paste("Cluster",LETTERS[i]),catnames),
                             pch=15,
                             gp=gpar(col=c(NA,paste0("#",Tol_muted))))
  g <- arrangeGrob(grobs=pi.list,layout_matrix = matrix(c(1,1,1,1,1,1,1,1,1,10,10,10),nrow=3,ncol=4))
  grid.newpage()
  grid.draw(g)
  ggsave(paste0("k-means cluster ",LETTERS[i],".pdf"),g,height=4.5,width=6,scale=2)
  
  #for (j in 1:9) pi.list[[j]] <- pieGrob(as.numeric(x[indx.clus.quants[j],]),fills=paste0("#",Tol_muted))
  #pi.list[[10]] <- legendGrob(catnames,pch=15,gp=gpar(cex=0.7,col=paste0("#",Tol_muted)))
  #grid.arrange(grobs=pi.list,layout_matrix = matrix(c(1,4,7,2,5,8,3,6,9,10,10,10),nrow=3,ncol=4))
  #mtext(i,side=3)
}


cvi$State <- cvi.toxpi.df$STATE
cvi.state.cluster.tab <- table(cvi[,c("State","cluster")])
cvi.state.cluster.tab.prop <- sweep(cvi.state.cluster.tab,1,base::apply(cvi.state.cluster.tab,1,sum),FUN="/")

stateorder <- base::order(cvi.state.cluster.tab.prop[,"A"]+
                          cvi.state.cluster.tab.prop[,"C"]+
                          cvi.state.cluster.tab.prop[,"B"],
                          cvi.state.cluster.tab.prop[,"D"]+
                          cvi.state.cluster.tab.prop[,"E"]+
                          cvi.state.cluster.tab.prop[,"F"],
                          cvi.state.cluster.tab.prop[,"F"],
                          decreasing=TRUE)

cvi.state.cluster.tab.df <- as.data.frame(cvi.state.cluster.tab)
cvi.state.cluster.tab.df$State <- factor(cvi.state.cluster.tab.df$State,
                                         levels=rownames(cvi.state.cluster.tab)[stateorder])

pdf("k-means cluster-bystate.pdf",height=6,width=8)

print(ggplot(cvi.state.cluster.tab.df)+
  geom_col(aes(y=State,x=Freq,fill=cluster),position="fill")+
  scale_fill_viridis_d())

cvi.state.cluster.tab.prop.df <- as.data.frame(cvi.state.cluster.tab.prop)
cvi.state.cluster.tab.prop.df$State <- factor(cvi.state.cluster.tab.prop.df$State,
                                         levels=rownames(cvi.state.cluster.tab.prop)[stateorder])
print(ggplot(cvi.state.cluster.tab.prop.df)+
  geom_col(aes(x=cluster,y=Freq))+facet_wrap(~State))
dev.off()



# library(pca3d)
# pca.x <- princomp(x)
# pca2d(pca.x,group=x.kmeans$cluster)
# pca2d(pca.x,group=x.kmeans$cluster,components = 2:3)
# 
# pca.xsamp <- princomp(x[indx2,])
# pca3d(pca.xsamp,group=x.kmeans$cluster[indx2])
# 
# load(file.path(pctdir,paste0("CVI-pct-comb.Rdata")))
# 
# pca3d(as.matrix(xsamp),group=x.kmeans$cluster[indx2],legend="topleft")
# 
# pca2d(pca.x,group=x.kmeans$cluster,show.ellipses = TRUE,legend = "topleft")
# 


# plot(fcomb.pct.results[indx.clus,],y=txpRanks(fcomb.pct.results[indx.clus,]),pch=15)


# 
# baseline.toxpi.df <- fread(file.path("CVI-pct","CVI-pct-comb-baseline.csv"),integer64 = "double",
#                       keepLeadingZeros = TRUE)
# 
# climate.toxpi.df <- fread(file.path("CVI-pct","CVI-pct-comb-climate.csv"),integer64 = "double",
#                       keepLeadingZeros = TRUE)
