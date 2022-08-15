library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(choroplethr)
library(choroplethrMaps)
library(tidyr)
library(ggpubr)

pctdir <- "CVI-pct"
figdir <- "SuppFigures"

# 10 colors, color-blind friendly (# removed, all lower case)
Tol_muted <- tolower(c('88CCEE', '44AA99', '117733', '332288', 'DDCC77', '999933','CC6677', '882255', 'AA4499', 'DDDDDD'))

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

indicators.df<-fread("CVI_indicators_current.csv")
indicators.df$Subcategory <- 
  factor(indicators.df$Subcategory,
         levels=unique(indicators.df$Subcategory))

cvi.df<-fread("CVI_data_current.csv",
              keepLeadingZeros = TRUE,integer64 = "numeric")
cvi.pct.df <- fread(file.path(pctdir,"CVI_data_pct.csv"),
                    keepLeadingZeros = TRUE,integer64 = "numeric")
cvi.pct.df$GEOID.Tract <- cvi.pct.df$FIPS
subcategories <- unique(indicators.df$Subcategory)

for (i in 1:length(subcategories)) {
  onesubcat <- subcategories[i]
  subcategoryname <- onesubcat
  print(subcategoryname)
  
  parameters <- unique(indicators.df$Parameters[
    indicators.df$Subcategory==onesubcat])
  numparm <-length(parameters)

  onecat <- indicators.df$Category[indicators.df$Subcategory==onesubcat &
                                     indicators.df$Parameters==parameters[1]]
  
  cvi.pct.toxpi.cat  <- fread(file.path(pctdir,
                                        paste0("CVI-pct-cat-",
                                               gsub(": ","-",onecat),".csv")),
                              keepLeadingZeros = TRUE,integer64 = "numeric")
  indx <- c(3,grep(onesubcat,names(cvi.pct.toxpi.cat)))  
  scores.df <- cvi.pct.toxpi.cat[,..indx]
  names(scores.df)[1] <- "GEOID.Tract"
  indx <- c("GEOID.Tract",parameters)
  scores.df <- left_join(scores.df,cvi.pct.df[,..indx])
  scores.df <- left_join(scores.df,cvi.df[,1:5])
  scores.df$state <- factor(scores.df$STATE)
  scores.df$county <- factor(paste(scores.df$County_Name,scores.df$STATE,sep=", "))
  scores.df$`Census Region` <- state.regions[scores.df$STATE,"Census Region"]
  scores.df$`Census Region` <- factor(scores.df$`Census Region`,levels=
                                        c("Northeast","Midwest","South","West"))
  names(scores.df)[2] <- "Subcategory Score"
  
  ## State aggregates
  indx<-c(2+0:numparm,grep("state",names(scores.df)))
  scores.df.med <- aggregate(.~state,FUN=median,data=scores.df[,..indx])
  scores.df.max <- aggregate(.~state,FUN=max,data=scores.df[,..indx])
  
  ############# County level maps
  indx<-c(2+0:numparm,grep("GEOID.County",names(scores.df)))
  scores.df.med.county <- aggregate(.~GEOID.County,
                                    FUN=median,data=scores.df[,..indx])
  ## Alternative metrics to map
  scores.df.mean.county <- aggregate(.~GEOID.County,
                                     FUN=mean,data=scores.df[,..indx])
  
  scores.df.var.county <- aggregate(.~GEOID.County,
                                    FUN=function(x) {if (length(x)>1) var(x) else 0},data=scores.df[,..indx])
  
  scores.df.max.county <- aggregate(.~GEOID.County,
                                    FUN=max,data=scores.df[,..indx])
  indx<-2+0:numparm
  scores.df.relvar.county <- cbind(GEOID.County=scores.df.var.county[,1],
                                   as.data.frame(sweep(data.matrix(scores.df.var.county[,indx]),
                                                       2,base::apply(as.matrix(scores.df[,..indx]),2,FUN=var),
                                                       FUN = "/")))
  ## Category-Subcategory Maps
  subcatmap.list <-list()
  parmmap.list <- list()
  viridisopts <- c("A","B","C","D","E","A","B","C","B","C")
  for (j in 1:(numparm+1)) {
    dat.df <- data.frame(region=as.numeric(scores.df.med.county$GEOID.County),
                         value=scores.df.med.county[[j+1]])
    plt<-CountyChoropleth$new(dat.df)
    plt$set_num_colors(1)
    plt$set_zoom(NULL)
    plt$ggplot_scale <- list(scale_fill_viridis_c("",option=viridisopts[j],limits=c(0,1)),
                             scale_color_viridis_c("",option=viridisopts[j],limits=c(0,1)))
    if (j == 1) {
      plt$title<-paste0("     ",subcategoryname)
    } else {
      plt$title<-paste0("     ",names(scores.df.med.county)[j+1])
    }
    plt$ggplot_polygon <- geom_polygon(aes(fill = value,color=value))
    if (j == 1) subcatmap.list[[1]] <- plt$render() else parmmap.list[[j-1]] <- plt$render()
  }
  indx <- 1+1:numparm
  numrow <- ceiling(numparm/2)
  figheight <- 1.75*(numrow+2)
  figmap <- ggarrange(subcatmap.list[[1]],
                      ggarrange(plotlist=parmmap.list,
                                nrow=numrow,ncol=2,labels=letters[indx]),
                      labels=c(letters[1],""),
                      nrow=2,ncol=1,heights = c(2,numrow))
  ggsave(file.path(figdir,
                   paste0("CVI_subcategorymap-",i,"-",
                          gsub(": ","-",as.character(onesubcat)),"-map.pdf")),
         figmap,height=figheight,width=6.5,scale=2)
  
  ## 
  
  ############# Source of variability
  indx <- 2+0:numparm
  scores.df.long <- pivot_longer(scores.df,cols = all_of(indx))
  scores.df.long$name <- factor(scores.df.long$name,
                                levels=unique(scores.df.long$name)
  )
  scores.df.long$state <- factor(scores.df.long$state,
                                 levels=as.character(
                                   scores.df.med$state[order(scores.df.med[[2]])]))
  scores.df.long$label <- factor(scores.df.long$name,levels=unique(scores.df.long$name))
  
  tmp.topcat <- 
    aggregate(value ~ GEOID.Tract,
              data = subset(scores.df.long,as.numeric(name) != 1),
              FUN = which.max)
  tmp.topcat$TopCat <- factor(levels(scores.df.long$name)[tmp.topcat$value+1],
                              levels=rev(levels(scores.df.long$name)))
  scores.df <- left_join(scores.df,tmp.topcat[,c(1,3)])
  
  ##### Summary boxplots
  scoresboxplt<-
    ggplot(scores.df.long)+geom_boxplot(aes(x=value,y=state,fill=`Census Region`),outlier.size = 0.2)+
    scale_fill_viridis_d(begin=0.3)+
    facet_wrap(~label,nrow=1)+xlab("")+theme_bw()+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
    ggtitle(subcategoryname)
  print(scoresboxplt)
  
  ##### Geographic scale
  subcat.indicators <- subset(indicators.df,Subcategory == onesubcat)
  subcat.indicators$Parameters <- factor(subcat.indicators$Parameters,
                                       levels=unique(subcat.indicators$Parameters))
  indicators.geo <- 
    rbind(
      as.matrix(t(table(subcat.indicators$GeographicScale))),
      as.matrix(ftable(subcat.indicators$Parameters,
                       subcat.indicators$GeographicScale)))
  rownames(indicators.geo)[1] <- "Subcategory Score"
  # rownames(indicators.geo) <- gsub("& ","&\n",rownames(indicators.geo))
  # rownames(indicators.geo) <- gsub("and ","and\n",rownames(indicators.geo))
  indicators.geo <- as.data.frame(indicators.geo)
  indicators.geo$label <- factor(rownames(indicators.geo),levels=
                                   rev(rownames(indicators.geo)))
  indicators.geo.df <- pivot_longer(indicators.geo,1:(ncol(indicators.geo)-1))
  indicators.geo.df$name <- factor(indicators.geo.df$name,
                                   levels=rev(c("State","County","Tract","Tract (raster)")))
  pgeo <- ggplot(indicators.geo.df)+
    geom_col(aes(x=label,y=value,fill=name),position="fill")+
    scale_y_continuous(label = scales::percent)+
    scale_fill_viridis_d(begin=0,end=0.9,option="C",drop=FALSE)+
    labs(x="",y="Geographic Scale (% of indicators)")+theme_bw()+
    coord_flip()+
    theme(legend.title = element_blank(),
          axis.text.y=element_text(hjust=0))+
    guides(fill = guide_legend(reverse = TRUE))
  print(pgeo)
  
  ##### Heterogeneity - state and county level
  indx<-c(2+0:numparm,grep("STATE",names(scores.df)))
  scores.df.state.SSres <- aggregate(.~STATE,data=scores.df[,..indx],
                                     FUN=function(x) {var(x)*length(x)})
  indx<-2+0:numparm
  state.r2 <- 1 - (base::apply(scores.df.state.SSres[,-1],2,sum))/
    (base::apply(scores.df[,..indx],2,var)*nrow(scores.df))
  # names(state.r2) <- gsub("& ","&\n",names(state.r2))
  # names(state.r2) <- gsub("and ","and\n",names(state.r2))
  
  indx<-c(2+0:numparm,grep("GEOID.County",names(scores.df)))
  scores.df.county.SSres <- aggregate(.~GEOID.County,data=scores.df[,..indx],
                                      FUN=function(x) {if (length(x)>1) (var(x)*length(x)) else 0})
  indx<-2+0:numparm
  county.r2 <- 1 - (base::apply(scores.df.county.SSres[,-1],2,sum))/
    (base::apply(scores.df[,..indx],2,var)*nrow(scores.df))
  # names(county.r2) <- gsub("& ","&\n",names(county.r2))
  # names(county.r2) <- gsub("and ","and\n",names(county.r2))
  
  r2.df <- rbind(data.frame(Variable=rep("State",1+numparm),
                            Parameter=names(state.r2),
                            R2=state.r2),
                 data.frame(Variable=rep("County",1+numparm),
                            Parameter=names(county.r2),
                            R2=county.r2)
  )
  r2.df$Parameter<-factor(r2.df$Parameter,levels=
                              rev(names(state.r2)))
  pr2<-ggplot(r2.df)+
    geom_col(aes(x=Parameter,y=R2,
                 group=Variable,fill=Variable),position="dodge")+
    scale_y_continuous(label = scales::percent,limits=c(0,1))+
    scale_fill_viridis_d(begin=0.6,end=0.9,option="C")+
    labs(x="",y=bquote('% variance due to state or county '(R^2)))+
    coord_flip()+theme_bw()+theme(legend.title = element_blank(),
                                  axis.text.y=element_text(hjust=0))+
    guides(fill = guide_legend(reverse = TRUE))
  
  print(pr2)
  ###### Bar graph of top subcategories
  scores.df.cat <- cbind(data.frame(Sample=rep("All",1+numparm)),
                         as.data.frame(prop.table(table(scores.df$TopCat))))
  scores.df.cat <- rbind(scores.df.cat,
                         cbind(data.frame(Sample=rep("Top quartile",1+numparm)),
                               as.data.frame(prop.table(table(
                                 (scores.df[scores.df[[2]]  >= quantile(scores.df[[2]],prob=0.75),])$TopCat)))
                         )
  )
  scores.df.cat <- rbind(scores.df.cat,
                         cbind(data.frame(Sample=rep("Top decile",1+numparm)),
                               as.data.frame(prop.table(table(
                                 (scores.df[scores.df[[2]]  >= quantile(scores.df[[2]],prob=0.9),])$TopCat)))
                         )
  )
  names(scores.df.cat)[3] <- "Fraction of census tracts"
  # scores.df.cat$Var1 <- gsub("& ","&\n",scores.df.cat$Var1)
  # scores.df.cat$Var1 <- gsub("and ","and\n",scores.df.cat$Var1)
  
  scores.df.cat$`Dominant Parameter` <- factor(
    as.character(scores.df.cat$Var1),
    levels=unique(as.character(scores.df.cat$Var1)))
  scores.df.cat <- subset(scores.df.cat,(`Dominant Parameter` %in% parameters))
  scores.df.cat$Sample <- factor(scores.df.cat$Sample,
                                 levels=c("All","Top quartile","Top decile"))
  pcat <- ggplot(scores.df.cat)+
    geom_col(aes(x=`Dominant Parameter`,y=`Fraction of census tracts`,
                 fill=Sample,group=Sample),position="dodge")+
    scale_y_continuous(label = scales::percent,limits=c(0,1))+
    scale_fill_viridis_d(begin=0.2,end=0.8,option="inferno")+xlab("")+
    labs(subtitle = "Dominant Parameter",y="% of census tracts")+
    coord_flip()+theme_bw()+theme(legend.title = element_blank(),
                                  axis.text.y=element_text(hjust=0))+
    guides(fill = guide_legend(reverse = TRUE))
  print(pcat)

  indx<-2+1:numparm
  scores.cor <- cor(scores.df[,..indx])
  # rownames(scores.cor) <- gsub("& ","&\n",rownames(scores.cor))
  # rownames(scores.cor) <- gsub("and ","and\n",rownames(scores.cor))
  # colnames(scores.cor) <- gsub("& ","&\n",colnames(scores.cor))
  # colnames(scores.cor) <- gsub("and ","and\n",colnames(scores.cor))
  pcor <- ggcorrplot(scores.cor, show.diag = TRUE,
                     type ="lower", lab =TRUE,tl.cex=10)+
    theme(legend.position = "none",axis.text.x=element_blank())
    
  pfig <- ggarrange(scoresboxplt,
                    ggarrange(pgeo,pr2,pcat,pcor,labels=c("b","c","d","e"),ncol=2,nrow=2),
                    labels=c("a",""),ncol=1,heights=c(2,2))
  ggsave(file.path(figdir,
                     paste0("CVI_subcategorymap-",i,"-",
                            gsub(": ","-",as.character(onesubcat)),"-summary.pdf")),
         pfig,height=6,width=6.5,scale=2)
}
