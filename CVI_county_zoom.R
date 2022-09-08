library(data.table)
library(stringr)
library(dplyr)
library(tidyr)
library(naniar)
library(readxl)
library(ggplot2)
library(choroplethr)
library(choroplethrMaps)
library(ggpubr)
library(lsr)
library(toxpiR)
library(grid)
library(tigris)


figdir <- "Figures"
supfigdir <- "SuppFigures"
datafolder <- "Data"
pctdir <- "CVI-pct"
catnames <-   c(  "Baseline Vulnerability:\nHealth",
                  "Baseline Vulnerability:\nSocial and Economic",
                  "Baseline Vulnerability:\nInfrastructure",
                  "Baseline Vulnerability:\nEnvironment",
                  "Climate Change Risk:\nHealth",
                  "Climate Change Risk:\nSocial and Economic",
                  "Climate Change Risk:\nExtreme Events"
)

options(tigris_year=2019)
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

# 10 colors, color-blind friendly (# removed, all lower case)
Tol_muted <- tolower(c('88CCEE', '44AA99', '117733', '332288', 'DDCC77', '999933','CC6677', '882255', 'AA4499', 'DDDDDD'))

Tol_max <- tolower(c('009FEE','00AA8E','007728','170088','DDB800','999900','CC0022','880044','AA008E','010101'))

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

cvi.toxpi.df <- fread(file.path(pctdir,"CVI-pct-comb-clusters.csv"),integer64 = "double",
                      keepLeadingZeros = TRUE)
cvi.toxpi.df <- left_join(cvi.toxpi.df,tractsdat[,c("County_Name","FIPS")])
cvi.toxpi.df$pctrank <- round(100*(rank(cvi.toxpi.df$`ToxPi Score`)-1)/(nrow(cvi.toxpi.df)-1),3)

indicators.df<-fread("CVI_indicators_current.csv")
indicators.df$`Category` <- 
  factor(indicators.df$`Category`,
         levels=unique(indicators.df$`Category`))

cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE,integer64 = "numeric")

cvi.pct.df<-fread(file.path(pctdir,"CVI_data_pct.csv"),
              keepLeadingZeros = TRUE,integer64 = "numeric")
cvi.pct.df$GEOID.Tract <- cvi.pct.df$FIPS

categories <- unique(indicators.df$`Category`)

cvi.toxpi.cat.list <- list()

for (i in 1:length(categories)) {
  onecat <- categories[i]
  cvi.toxpi.cat.list[[i]]  <- fread(file.path(pctdir,
                                        paste0("CVI-pct-cat-",
                                               gsub(": ","-",onecat),".csv")),
                              keepLeadingZeros = TRUE,integer64 = "numeric")
}

onecounty.df <- subset(cvi.toxpi.df,STATE=="TX" & County_Name=="Harris")
county_zoom <- as.numeric(onecounty.df$GEOID.County[1])
dat.df <- data.frame(region=as.numeric(onecounty.df$FIPS),
                     value=onecounty.df$`ToxPi Score`)
plt<-TractChoropleth$new("texas",dat.df)
plt$set_num_colors(1)
plt$set_zoom_tract(county_zoom=county_zoom,tract_zoom=NULL)
plt$ggplot_scale <- list(scale_fill_viridis_c("",option="A",limits=range(cvi.toxpi.df$`ToxPi Score`)),
                         scale_color_viridis_c("",option="A",limits=range(cvi.toxpi.df$`ToxPi Score`)))
plt$title<-""
plt$ggplot_polygon <- geom_polygon(aes(fill = value,color=value))
p <- plt$render()
p

toptract <- onecounty.df[order(onecounty.df$pctrank,decreasing = TRUE)[1:3],]
toptract$label <- paste0("Tract ",toptract$FIPS,"\nCVI Score: ",
                         round(toptract$`ToxPi Score`,4),
                         " (",toptract$pctrank,"%)")
tractorder<-order(toptract$`ToxPi Score`,decreasing = FALSE)
toptract <- toptract[tractorder,]
plt<-TractChoropleth$new("texas",dat.df)
plt$set_num_colors(1)
tract_zoom <- as.numeric(toptract$FIPS)
plt$set_zoom_tract(county_zoom=county_zoom,tract_zoom=tract_zoom)
plt$title<-""
plt$ggplot_scale <- list(scale_fill_viridis_c("",option="A",limits=range(cvi.toxpi.df$`ToxPi Score`)),
                         scale_color_viridis_c("",option="A",limits=range(cvi.toxpi.df$`ToxPi Score`)))
plt$ggplot_polygon <- geom_polygon(aes(fill = value,color=value))
p2 <- plt$render()+geom_text(data=toptract,aes(INTPTLON10,INTPTLAT10,label=label),size=3)+
  theme(legend.position = "none")

ppi_fill <- pieGridGrob(1+0*as.matrix(toptract[,6:12]),
                        labels=NULL,
                        fills=paste0("#",Tol_muted[10]),vp=viewport(width=0.9, height=2, angle=90),
                        #gp=gpar(cex=0.8),
                        ncol=1)#nrow=1)

ppi <- pieGridGrob(as.matrix(toptract[3:1,6:12]),
                               labels=NULL,
                               fills=paste0("#",Tol_muted),vp=viewport(width=0.9, height=2, angle=90),
                               #gp=gpar(cex=0.8),
                               ncol=1)#nrow=1)
onetract <- toptract[order(toptract$`ToxPi Score`,decreasing = TRUE)[1],]


ppi_cat.list <- list()
ppi_leg.list <- list()
for (i in 1:length(categories)) {
  onecat <- categories[i]
  cvi.toxpi.cat<-cvi.toxpi.cat.list[[i]]
  onetract.cat <- subset(cvi.toxpi.cat,FIPS==onetract$FIPS)
  nsubcat <- ncol(onetract.cat)-5
  fillcols <- paste0("#",paste0(Tol_max[i],as.hexmode(round((1+(nsubcat:1))*255/(1+nsubcat)))))
  ppi_onefill <- pieGridGrob(as.matrix(1),
                             labels=NULL,
                             fills="#EEEEEE",vp=viewport(width=0.75, height=0.75,angle=90),
                             gp=gpar(cex=1),
                             ncol=1)
  ppi_onelab <- pieGridGrob(as.matrix(1),
                             labels=catnames[i],
                             fills="#FFFFFF00",vp=viewport(width=0.9, height=0.9),
                             gp=gpar(cex=1),
                             ncol=1)
  ppi_cat.list[[i]]<-gList(ppi_onefill,
    pieGridGrob(as.matrix(onetract.cat[,6:ncol(onetract.cat)]),
                labels = NULL,
          fills=fillcols,
          vp=viewport(width=0.75, height=0.75,angle=90),
          gp=gpar(cex=1)),
    ppi_onelab
    )
  ppi_leg.list[[i]] <- legendGrob(labels=names(onetract.cat)[6:ncol(onetract.cat)],
                                  pch=15,gp=gpar(cex=0.7,
                                                 col=fillcols),
                                  vp=viewport(width=0.9, height=0.9))
}

psubcat <- ggarrange(plotlist = c(ppi_cat.list,ppi_leg.list),nrow=2,ncol=7,
                     heights=c(1,1))

pzoom <- ggarrange(ggarrange(p,ggarrange(p2,gList(ppi_fill,ppi),ncol=1,
                                         heights=c(2,1),
                             labels=c("Top Three Ranking Tracts","")),
                             ncol=2,
                             widths = c(3,2.5)),
                   psubcat,nrow=2,heights=c(2,1.25),
                   labels=c("Harris County, TX",paste("Top Ranked Tract:",last(onetract$FIPS))),
                   label.y=c(1,1.1),label.x=c(0,-0.03)
            )
ggsave(file.path(figdir,"Harris County Zoom.pdf"),pzoom,height=5,width=6,scale=2.5)


pscores.list <- list()

for (i in 1:length(categories)) {
  onecat <- categories[i]
  cvi.toxpi.cat<-cvi.toxpi.cat.list[[i]]
  onetract.cat <- subset(cvi.toxpi.cat,FIPS==onetract$FIPS)
  
  parameters <- indicators.df$Parameters[
    indicators.df$Category==onecat]
  numparms <- length(parameters)
  subcategories <- indicators.df$Subcategory[
    indicators.df$Category==onecat]
  indx <- c("FIPS",parameters)
  onetract_scores.df <- pivot_longer(subset(cvi.pct.df,
                                            FIPS==onetract$FIPS)[,..indx],
                                     cols = 1+1:numparms)
  onetract_scores.df$name <- factor(onetract_scores.df$name,
                                    levels=rev(parameters))
  onetract_scores.df$Category <- onecat
  onetract_scores.df$Subcategory <- factor(subcategories,
                                           levels=unique(subcategories))
  pscores.list[[i]]<-
    ggplot(onetract_scores.df,aes(x=value,y=name,fill=Subcategory))+geom_col()+
    scale_y_discrete(position="right",
                     labels = function(x) str_wrap(x, width = 40))+
    scale_x_continuous(label = scales::percent,limits=c(0,1))+
    scale_fill_manual(values=paste0("#",Tol_muted))+
    theme_bw()+
    theme(legend.position ="none",axis.title.y = element_blank())+
    ggtitle(paste0(catnames[i]))+
    xlab("Percentile")+
    facet_wrap(~Subcategory,scale="free_y",ncol=1)
    #facet_grid(Subcategory~.,scales="free_y",space="free_y")

  # ggplot(onetract_scores.df,aes(x=value,y=name,fill=Subcategory))+geom_col()+
  #   scale_y_discrete(position="right")+
  #   scale_x_continuous(label = scales::percent,limits=c(0,1))+
  #   scale_fill_manual(values=paste0("#",Tol_muted))+
  #   theme_bw()+
  #   ggtitle(paste0(catnames[i]))+
  #   xlab("Percentile")+ylab("Parameter")
  # 
}

pscores.baseline1 <- ggarrange(plotlist = pscores.list[1:3],nrow=1,ncol=3)
pscores.climate <- ggarrange(plotlist = pscores.list[5:7],nrow=1,ncol=3)
pscores <- ggarrange(ggarrange(pscores.baseline1,pscores.climate,
                               nrow=2,ncol=1,
                               heights=c(4,3)),
                     pscores.list[[4]],widths=c(3,1))
ggsave(file.path(supfigdir,"Harris County TopTract scores.pdf"),pscores,height=5.5,width=6,scale=4)
