#############################################################################
####Species NMDS Nech Sar - Oct 2012 - Nov 2013 ####
#Stuart Smith
#12/10/2018 
############################################################################
# Libraries
rm(list=ls())
library(Matrix)
library(lattice)
library(MASS)
library(tidyr)
library(lubridate)
library(vegan)
library(ggplot2)
library(nlme)
library(dplyr)
############################################################################

nsSpp<-read.csv("VegCompSeason3.csv",header=T,sep=",")

names(nsSpp)
head(nsSpp)
str(nsSpp)

nsreharvest3<-read.table("ProductivitySeason3.txt",header=T,sep="\t")
names(nsreharvest3)

# Only exclosure and open treatments - not encloure/rodent work etc.
nsSpp2<-droplevels(nsSpp[nsSpp$Treatment=="Control" | nsSpp$Treatment=="Exclosure",])
nsreharvest3<-droplevels(nsreharvest3[nsreharvest3$Treatment=="Control" | nsreharvest3$Treatment=="Exclosure",])

# Add regrowth data to species dataset
nsreharvest3<-droplevels(nsreharvest3[!is.na(nsreharvest3$Harvest.date),])
nsreharvest3<-droplevels(nsreharvest3[nsreharvest3$Harvest=="double",])
nsreharvest3$Season<-nsreharvest3$Reharvest.date
levels(nsreharvest3$Season)<-c("Short I" ,  "Long" , "Short II")
levels(nsSpp2$Season)<-c("Short I" ,  "Long" , "Short II")
nsSpp2$Trt.name<-as.factor(nsSpp2$Trt.name)
nsSpp2$Season<-as.factor(nsSpp2$Season)
nsreharvest3$Season <- factor(nsreharvest3$Season, levels(nsreharvest3$Season)[c(3,2,1)])
nsSpp2$plot_code<-as.factor(with(nsSpp2, paste(Transect,Block,Trt.name,Season,Replicate, sep="_")))
nsreharvest3$plot_code<-as.factor(with(nsreharvest3, paste(Transect,Block,Trt.name,Season,Replicate, sep="_")))

dim(nsSpp2)
dim(nsreharvest3)

nsreharvest3b<-nsreharvest3[c("Transect","Block","Treatment","Trt.name","Season","Replicate","plot_code","Livestock.density","GrassNetReharvestBiomass1","DwarfShrubNetReharvestBiomass1",
                              "HerbNetReharvestBiomass1","ClimberNetReharvestBiomass1","TotalBiomass1")]
nsSpp3<-left_join(nsSpp2,nsreharvest3b, by=c("Transect","Treatment","Block","Trt.name","Season","Replicate","Livestock.density","plot_code"),drop=F)

# Block as a factor
nsSpp3$fBlock<-as.factor(nsSpp3$Block)

#### SEPERATE PLANT FUNCTIONAL GROUP DATASETS ####

# Combine herbivores and climbers
levels(nsSppFx$Fx.group)<-c("Herbs","Grass","Dwarf shrub","Grass","Herb","Dwarf shrub")

# Repeat analysis - but just for grasses, woody, herbs and climbers
nsSppFx<-read.csv("Spp.Fx.group.csv",header=T,sep=",")
nsSppFxG<-droplevels(nsSppFx[nsSppFx$Fx.group!="Grass",])
nsSppFxW<-droplevels(nsSppFx[nsSppFx$Fx.group!="Dwarf shrub",])
nsSppFxC<-droplevels(nsSppFx[nsSppFx$Fx.group!="Herb",])

# Remove . and spaces = same names..Grasses
names(nsSpp3)  <- gsub("\\.", "",  names(nsSpp3))
to.remove  <- gsub(" ", "", levels(nsSppFxG$Species))
to.remove  <- gsub("\\.", "", to.remove )
to.remove  <- gsub("[()]", "",to.remove)

# Remove . and spaces = same names..Woody
to.removeW  <- gsub(" ", "", levels(nsSppFxW$Species))
to.removeW  <- gsub("\\.", "", to.removeW )
to.removeW  <- gsub("[()]", "",to.removeW)

# Remove . and spaces = same names..Herbs and Climbers
to.removeC  <- gsub(" ", "", levels(nsSppFxC$Species))
to.removeC  <- gsub("\\.", "", to.removeC )
to.removeC  <- gsub("[()]", "",to.removeC)

`%ni%` <- Negate(`%in%`)
nsSpp3G<-nsSpp3[ , !(names(nsSpp3) %in% to.remove)]
#nsSpp2G<-subset(nsSpp2,select = names(nsSpp2) %ni% to.remove)
names(nsSpp3G) # Subset is now just grasses

# Subset woody - dwarf shrubs
nsSpp3W<-nsSpp3[ , !(names(nsSpp3) %in% to.removeW)]
nsSpp3W<-subset(nsSpp3,select = names(nsSpp3) %ni% to.removeW)
names(nsSpp3W) # Subset is now just dwarf shrubs

# Subset Herbs and climbers
nsSpp3C<-nsSpp3[ , !(names(nsSpp3) %in% to.removeC)]
nsSpp3C<-subset(nsSpp3,select = names(nsSpp3) %ni% to.removeC)
names(nsSpp3C) # Subset is now herbs and climber

################################################################################
#### Non multi-dimensional scaling - explore spp data ####
################################################################################

# USE NMDS
names(nsSpp3[,10:74]) # 65 species # Rosett?
mdsNS<-metaMDS(nsSpp3[,10:74], trace =F) 
mdsNS
#Stress:    0.2910474    # Not great stress - but under 0.3

# Stressplot
vare.dis<-vegdist(nsSpp3[,10:74]) 
vare.dis2<-as.matrix(vare.dis)
vare.dis
vare.mds0<-isoMDS(vare.dis) 
stressplot(vare.mds0,vare.dis) # final  value 24.971572
# Metri = 0.94, LinearR2 =0.73

# Basic NMDS plot
par(mfrow=c(1,1))
plot(mdsNS, type="p")
plot(mdsNS$points[,1]~mdsNS$points[,2]) 
# No major outliers - strong central clustering

# Enviro fit - fit season, treatment, livestock density 
# Create a factor to combining livestock density, treatment and date
nsSpp3$harvest_code<-as.factor(with(nsSpp3, paste(Livestockdensity,Treatment,Date, sep="-")))

# Envfit
NS.ev <- envfit(mdsNS~Livestockdensity+Treatment+ Season,data = nsSpp3, 
                  perm=999,strata=as.numeric(nsSpp3$Block))
NS.ev 
# Weak difference in season - biggest difference is livestock density - r2=0.26

NS.har <- envfit(mdsNS~ harvest_code,data = nsSpp3, 
                perm=999,strata=as.numeric(nsSpp3$Block))
NS.har

# Extract centroids
vec.sp.df<-as.data.frame(cbind(NS.har$factors$centroids*sqrt(NS.har$factors$r)))
vec.sp.df$Livestockdensity<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                               "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.df$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                       "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.df$Season<-rep(c(1,3,2))

names(nsSpp3)
TotBio<-aggregate(TotalBiomass1~harvest_code,nsSpp3,mean)
vec.sp.df$TotalBiomass1<-TotBio[,2]
vec.sp.df$TotalBiomass1<-as.numeric(vec.sp.df$TotalBiomass1)

# Reorder by season and Harvest code
vec.sp.df$Season2<- factor(vec.sp.df$Season, levels = c("1","2","3"))
vec.sp.df<- vec.sp.df[order(vec.sp.df$Season2),] 
vec.sp.df$harvest_code<-as.factor(with(vec.sp.df, paste(Livestockdensity,Treatment, sep="-")))

# Draw ordiellipse in ggplot2
# NEEDS TO CONNECT TO SPP MULTIVARIATE SCRIPT...
#NMDS = data.frame(NMDS1 = mdsNS$points[,1], NMDS2 = mdsNS$points[,2],group=nsSpp3$Livestockdensity)
#NMDS.mean=aggregate(cbind(NMDS$NMDS1,NMDS$NMDS2)~group,data=NMDS,mean)
#NMDS = data.frame(NMDS1 = SerEbio5$NMDS1, NMDS2 = SerEbio5$NMDS2,group=SerEbio5$fn.non.N)
#NMDS.mean=aggregate(cbind(NMDS1,NMDS2)~fn.non.N,SerSppfullwide,mean)

#NMDS<-na.omit(NMDS)
#plot.new()
#ord<-ordiellipse(mdsNS,nsSpp3$Livestockdensity, conf=0.95)
#ord<-ordiellipse(mdsSero,SerSppfullwide$Trt_code, conf=0.95)

#Sero.evHR <- envfit(mdsSero ~ Trt_code,data = SerSppfullwide)
#plot(Sero.evHR$factors$centroids[,1]~Sero.evHR$factors$centroids[,2])

# Extracting the path for the ordiellipse - to be drawn in ggplot2
#veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
#{
#  theta <- (0:npoints) * 2 * pi/npoints
#  Circle <- cbind(cos(theta), sin(theta))
#  t(center + scale * t(Circle %*% chol(cov)))
#}

#df_ell <- data.frame()
#for(g in levels(NMDS$group)){
#  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
#                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
#                                ,group=g))
#}

#df_ell <- data.frame()
#for(g in levels(NMDS$group)){
#  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
#                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
#                                ,group=g))
#}#

#colnames(df_ell)[3]<-"Livestockdensity"
#df_ell$NMDS1<-df_ell$NMDS1/4
#df_ell$NMDS2<-df_ell$NMDS2/4

# Plot centroids
CenPlot<-ggplot(vec.sp.df[order(vec.sp.df$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity))
#CenPlot<-CenPlot+geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=Livestockdensity), size=1,show.legend=T)
CenPlot<-CenPlot+geom_point(aes(shape=Treatment,size=TotalBiomass1), stroke=1)
CenPlot<-CenPlot+geom_path(aes(group=harvest_code),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
CenPlot<-CenPlot +geom_text(aes(label=Season),hjust=0, vjust=-.95, show.legend = F)
CenPlot<-CenPlot +scale_colour_manual(values=c("black","grey80","grey50"))
CenPlot<-CenPlot +scale_fill_manual(values=c("black","grey80","grey50"))
CenPlot<-CenPlot +scale_radius(range=c(1,8))
CenPlot<-CenPlot +scale_shape_manual(values=c(21,22))
CenPlot<-CenPlot +ggtitle("Total biomass")
CenPlot<-CenPlot + #theme_bw() +
  theme(rect = element_rect(fill ="transparent")
        ,panel.background=element_rect(fill="transparent")
        ,plot.background=element_rect(fill="transparent",colour=NA)
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,panel.grid.major.x = element_blank()
        ,panel.grid.major.y = element_blank()
        ,axis.text=element_text(size=12,color="black")
        ,axis.title.y=element_text(size=12,color="black")
        ,axis.title.x=element_text(size=12,color="black")
        ,axis.text.x=element_text(size=11,color="black",
                                  margin=margin(2.5,2.5,2.5,2.5,"mm"))
        ,axis.ticks.length=unit(-1.5, "mm")
        ,axis.text.y = element_text(margin=margin(2.5,2.5,2.5,2.5,"mm"))
        ,axis.text.y.right =element_text(margin=margin(2.5,2.5,2.5,2.5,"mm"))
        ,axis.line.y = element_line(color="black", size = .5)
        ,axis.line.x = element_line(color="black", size = .5)
        ,plot.margin = unit(c(2.5,2.5,2.5,2.5), "mm")
        ,strip.background = element_rect(fill="transparent",colour=NA)
        ,strip.text.x = element_text(size=12,margin = margin(.5,.5,.5,.5, "mm"),hjust = .02)
        ,strip.text.y = element_blank()
        ,panel.spacing = unit(.1, "lines")
        ,legend.text=element_text(size=12)
        ,legend.title=element_text(size=12)
        ,legend.position = "right"
        ,legend.justification = "top"
        ,legend.direction="vertical"
        ,legend.key.width = unit(1.2,"cm"))
CenPlot

# PERMANOVA/ADONIS - Total 
# Distance matrix
vare.dis<-vegdist(nsSpp3[,10:74]) 
vare.dis2<-as.matrix(vare.dis)
vare.dis

# PERMANOVA grasses
nsSpp3$time_code<-as.factor(with(nsSpp3, paste(Transect,Block,Treatment,Replicate, sep="-")))
PermT<-adonis(vare.dis2 ~ Livestockdensity+Treatment+Date+
                Date:Treatment+Livestockdensity:Date+
                Livestockdensity:Treatment+
                Treatment:Livestockdensity:Date, 
              strata=as.numeric(nsSpp3$Block)/as.numeric(nsSpp3$time_code),
              method = "bray",perm=999, data=nsSpp3)
PermT$aov.tab
PermT$aov.tab$R2
str(PermT)
install.packages("Rcmdr")
library(BiodiversityR)
Ordination.model1 <- CAPdiscrim(nsSpp3[,10:74]~Treatment:Livestockdensity:Date, data=nsSpp,
                                dist="bray", axes=2, m=0, add=FALSE)
Ordination.model1
plot1 <- ordiplot(Ordination.model1, type="none")
ordisymbol(plot1, dune.env, "Management", legend=TRUE)


#Contrast: High_Low                      average
#BothriochloainsculptaHochstExARichACamus 1.440e-01 # 0.144
#ChrysopogonplumulosusHochst              1.281e-01
#CynodonnlemfuensisVanderyst              1.153e-01
#HeteropogoncontortusLRoemSchult          7.716e-02
# Strongest species differences in grasses
boxplot(BothriochloainsculptaHochstExARichACamus~Livestockdensity, nsSpp3) # High in high
boxplot(ChrysopogonplumulosusHochst~Livestockdensity, nsSpp3) # High in high
boxplot(CynodonnlemfuensisVanderyst~Livestockdensity, nsSpp3) # High in Low livestock
boxplot(HeteropogoncontortusLRoemSchult~Livestockdensity, nsSpp3) # Higher in Low

#### Grasses ####
# USE NMDS on Grass only
names(nsSpp3G[,10:28]) # 
mdsNSG<-metaMDS(nsSpp3G[,10:28], trace =F) 
mdsNSG # 0.1550549

NS.evG <- envfit(mdsNSG~Livestockdensity+Treatment+ Season,data = nsSpp3G, 
                 perm=999,strata=as.numeric(nsSpp3G$Block))
NS.evG # Livestock density and treatment significant - not season

nsSpp3G$harvest_code<-as.factor(with(nsSpp3G, paste(Livestockdensity,Treatment,Date, sep="-")))
NS.harG <- envfit(mdsNSG~ harvest_code,data = nsSpp3G, 
                  perm=999,strata=as.numeric(nsSpp3G$Block))
NS.harG

# Extract centroids
#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258

vec.sp.dfG<-as.data.frame(cbind(NS.harG$factors$centroids*sqrt(NS.harG$factors$r)))
#vec.sp.dfG<-as.data.frame(scores(mdsNSG, display = "sites"))
vec.sp.dfG$Livestockdensity<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                                "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.dfG$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                        "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.dfG$Season<-rep(c(1,3,2))

GReharvestBio<-aggregate(GrassNetReharvestBiomass1~harvest_code,nsSpp3G,mean)
vec.sp.dfG$GrassNetReharvestBiomass1<-GReharvestBio[,2]
vec.sp.dfG$GrassNetReharvestBiomass1<-as.numeric(vec.sp.dfG$GrassNetReharvestBiomass1)
#NS.harG$scores
#Grass species scores
spp.scrs <- as.data.frame(scores(mdsNSG, display = "species"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Livestockdensity<-c("Low")
spp.scrs$Season<-c("Short I")
spp.scrs$Treatment<-c("Control")

# Reorder by season and Harvest code
vec.sp.dfG$Season2<- factor(vec.sp.dfG$Season, levels = c("1","2","3"))
vec.sp.dfG<- vec.sp.dfG[order(vec.sp.dfG$Season2),] 
vec.sp.dfG$harvest_code<-as.factor(with(vec.sp.dfG, paste(Livestockdensity,Treatment, sep="-")))

vec.sp.dfG$Livestockdensity<-as.factor(vec.sp.dfG$Livestockdensity)
levels(vec.sp.dfG$Livestockdensity)

# Draw ordiellipse in ggplot2
# NEEDS TO CONNECT TO SPP MULTIVARIATE SCRIPT...
#NMDS = data.frame(NMDS1 = mdsNSG$points[,1], NMDS2 = mdsNSG$points[,2],group=nsSpp3G$Livestockdensity)
#NMDS.mean=aggregate(cbind(NMDS$NMDS1,NMDS$NMDS2)~group,data=NMDS,mean)
#NMDS = data.frame(NMDS1 = SerEbio5$NMDS1, NMDS2 = SerEbio5$NMDS2,group=SerEbio5$fn.non.N)
#NMDS.mean=aggregate(cbind(NMDS1,NMDS2)~fn.non.N,SerSppfullwide,mean)

#NMDS<-na.omit(NMDS)
#plot.new()
#ord<-ordiellipse(mdsNS,nsSpp3$Livestockdensity, conf=0.95)
#ord<-ordiellipse(mdsSero,SerSppfullwide$Trt_code, conf=0.95)

#Sero.evHR <- envfit(mdsSero ~ Trt_code,data = SerSppfullwide)
#plot(Sero.evHR$factors$centroids[,1]~Sero.evHR$factors$centroids[,2])

# Extracting the path for the ordiellipse - to be drawn in ggplot2
#veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
#{
#  theta <- (0:npoints) * 2 * pi/npoints
#  Circle <- cbind(cos(theta), sin(theta))
#  t(center + scale * t(Circle %*% chol(cov)))
#}

#df_ell <- data.frame()
#for(g in levels(NMDS$group)){
#  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
#                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
#                                ,group=g))
#}

#df_ell <- data.frame()
#for(g in levels(NMDS$group)){
#  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
#                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
#                                ,group=g))
#}
#colnames(df_ell)[3]<-"Livestockdensity"
#df_ell$NMDS1<-df_ell$NMDS1/6
#df_ell$NMDS2<-df_ell$NMDS2/6

# Grass NMDS movement plot
CenPlotG<-ggplot(vec.sp.dfG[order(vec.sp.df$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity))
#CenPlotG<-CenPlotG+geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=Livestockdensity), size=1,show.legend=T)
CenPlotG<-CenPlotG+geom_point(aes(size=GrassNetReharvestBiomass1,shape=Treatment), stroke=1)
CenPlotG<-CenPlotG+geom_path(aes(group=harvest_code),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
#CenPlotG<-CenPlotG+geom_text(data=spp.scrs,aes(label=Species),colour="light grey")
CenPlotG<-CenPlotG+geom_text(aes(label=Season),hjust=0, vjust=-.95, show.legend = F)
CenPlotG<-CenPlotG+scale_colour_manual(values=c("black","grey80","grey50"))
CenPlotG<-CenPlotG+scale_fill_manual(values=c("black","grey80","grey50"))
CenPlotG<-CenPlotG+scale_radius(range=c(1,8))
CenPlotG<-CenPlotG+scale_shape_manual(values=c(21,22))
CenPlotG<-CenPlotG+ggtitle("Grass")
CenPlotG<-CenPlotG+ #theme_bw() +
  theme(rect = element_rect(fill ="transparent")
        ,panel.background=element_rect(fill="transparent")
        ,plot.background=element_rect(fill="transparent",colour=NA)
        #,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,panel.grid.major.x = element_blank()
        ,panel.grid.major.y = element_blank()
        ,axis.text=element_text(size=12,color="black")
        ,axis.title.y=element_text(size=12,color="black")
        ,axis.title.x=element_text(size=12,color="black")
        ,axis.text.x=element_text(size=11,color="black",
                                  margin=margin(2.5,2.5,2.5,2.5,"mm"))
        ,axis.ticks.length=unit(-1.5, "mm")
        ,axis.text.y = element_text(margin=margin(2.5,2.5,2.5,2.5,"mm"))
        ,axis.text.y.right =element_text(margin=margin(2.5,2.5,2.5,2.5,"mm"))
        ,axis.line.y = element_line(color="black", size = .5)
        ,axis.line.x = element_line(color="black", size = .5)
        ,plot.margin = unit(c(2.5,2.5,2.5,2.5), "mm")
        ,strip.background = element_rect(fill="transparent",colour=NA)
        ,strip.text.x = element_text(size=12,margin = margin(.5,.5,.5,.5, "mm"),hjust = .02)
        ,strip.text.y = element_blank()
        ,panel.spacing = unit(.1, "lines")
        ,legend.text=element_text(size=12)
        ,legend.title=element_text(size=12)
        ,legend.position = "right"
        ,legend.justification = "top"
        ,legend.direction="vertical"
        ,legend.key.width = unit(1.2,"cm"))
CenPlotG

# PERMANOVA/ADONIS - GRASSES
# Distance matrix
#nsSpp3GSI<-droplevels(nsSpp3G[nsSpp3G$Date!="21.10.2012",])
vare.disG<-vegdist(nsSpp3G[,10:28]) 
vare.dis2G<-as.matrix(vare.disG)

# PERMANOVA grasses
nsSpp3G$time_code<-as.factor(with(nsSpp3G, paste(Transect,Block,Treatment,Replicate, sep="-")))
PermG<-adonis(vare.dis2G ~ Livestockdensity+Treatment+Date+
                Date:Treatment+Livestockdensity:Date+
                  Livestockdensity:Treatment+
                  Treatment:Livestockdensity:Date, 
                strata=as.numeric(nsSpp3G$Block)/as.numeric(nsSpp3G$time_code),
              method = "bray",perm=999, data=nsSpp3G)
PermG$aov.tab
PermG$aov.tab$R2

## Repeated measure
library(permute)
library(vegan)

sp <- as.matrix(nsSpp3G[,10:28])
env<-(nsSpp3G$Livestockdensity)
time<-nsSpp3G$Date
nsSpp3G$time_code<-as.factor(with(nsSpp3G, paste(Transect,Block,Treatment,Replicate, sep="-")))
Block<-as.numeric(nsSpp3G$Block)
rep.mes <- as.numeric(nsSpp3G$time_code)

#### NMDS ####
#sol <- metaMDS(sp, trymax = 5)
#fit <- envfit(sol~nsSpp3G$Livestockdensity+nsSpp3G$Treatment+nsSpp3G$Livestockdensity:nsSpp3G$Treatment, permutations = 1) ## perms now won't work!
print(fit <- adonis(vare.dis2G ~ time,strata =Block/rep.mes,method = "bray",perm=1, data=nsSpp3G))
#print(fit <- adonis(vare.dis2G ~ Livestockdensity+Treatment+Livestockdensity:Treatment,
#                    strata=as.numeric(nsSpp3G$Block),method = "bray",perm=0, data=nsSpp3G))
fit$aov.tab$R2

### setting up frame for population of r2 values:
#3^60 = 4.239116e+28
B <- 999 ## number of perms
pop <- rep(NA, B + 1)
pop[1] <-fit$aov.tab[1,5]
#pop[1] <- fit$aov.tab[1, 5]

## set-up a Control object:
ctrl <-how(plots=Plots(strata =Block/rep.mes),within = Within(type = "series",,mirror = TRUE))

## Number of observations
nobs <- nrow(sp)

## check it works
matrix(shuffle(nobs, control = ctrl), ncol = 3, byrow = TRUE)

### loop:
set.seed(1)
for(i in 2:(B+1)){
idx <- shuffle(nobs, control = ctrl)
fit.rand <- adonis(vare.dis2G ~ time[idx],method = "bray",perm=1)
pop[i] <- fit.rand$aov.tab$R2[1]
}

### get the p-value:
print(pval <- sum(pop >= pop[1]) / (B + 1))
### 0.001 = time = significant - species differ through time points

hist(pop, xlab = "Population R2")
abline(v = pop[1], col = 2, lty = 3) 
# Significant variation through time

# Grass pecies most asscoicated with specific livestock x treatment x season combination
nsSpp3G.env<-nsSpp3G[,c("fBlock","Livestockdensity","Season","Treatment")]
nsSpp3GI<-droplevels(nsSpp3G[nsSpp3G$Date!="21.10.2012",])
nsSpp3G.envI<-droplevels(nsSpp3G.env[nsSpp3G.env$Date!="21.10.2012",])
simG <- with(nsSpp3G.envI, simper(nsSpp3GI[,10:28],nsSpp3GI$Livestockdensity),permutations=999)
summary(simG,ordered = TRUE)

# Plots specific species by livestockdensity, treatment and date
#average       sd  ratio
#BothriochloainsculptaHochstExARichACamus 0.181135 0.16779 1.0795
#ChrysopogonplumulosusHochst              0.179700 0.17066 1.0530
#CynodonnlemfuensisVanderyst              0.140528 0.24070 0.5838
#HeteropogoncontortusLRoemSchult          0.095644 0.18057 0.5297

# High and medium livestock - lower in low livestock
BotIns<-ggplot(nsSpp3G,aes(x=Livestockdensity,y=BothriochloainsculptaHochstExARichACamus,shape=Treatment,colour=Livestockdensity, linetype=Date))
BotIns<-BotIns+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
BotIns

HetCon<-ggplot(nsSpp3G,aes(x=Livestockdensity,y=HeteropogoncontortusLRoemSchult,shape=Treatment,colour=Livestockdensity, linetype=Date))
HetCon<-HetCon+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
HetCon # Huge increase in low livesock exclosures - season II

# Medium and livestock higher
ChrPlu<-ggplot(nsSpp3G,aes(x=Livestockdensity,y=ChrysopogonplumulosusHochst,shape=Treatment,colour=Livestockdensity, linetype=Date))
ChrPlu<-ChrPlu+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
ChrPlu

# Low livestock - higher - in exclosures
CynNle<-ggplot(nsSpp3G,aes(x=Livestockdensity,y=CynodonnlemfuensisVanderyst,shape=Treatment,colour=Livestockdensity, linetype=Date))
CynNle<-CynNle+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
CynNle # Huge increase in low livestock exclosures- season III


##########################################################################################
#### Species richness ####
# Paritioning the beta - alpha and gamma
library(betapart)
##########################################################################################
names(nsSpp3G)

# Create a factor to define each plot shared across seasons
nsSpp3G$plot_codeII<-as.factor(with(nsSpp3G, paste(Transect,Block,Trtname,Replicate, sep="_")))

# Compare diversity through time
levels(nsSpp3G$Date)
nsSpp3GSeasonI<-droplevels(nsSpp3G[nsSpp3G$Date=="21.10.2012",])
nsSpp3GSeasonII<-droplevels(nsSpp3G[nsSpp3G$Date=="21.11.2013",])
nsSpp3GSeasonIII<-droplevels(nsSpp3G[nsSpp3G$Date=="21.6.2013",])

# Presence and absence of species
names(nsSpp3GSeasonI[,10:28])
names(nsSpp3GSeasonII[,10:28])
names(nsSpp3GSeasonIII[,10:28])
GpresabsI<-ifelse(nsSpp3GSeasonI[,10:28]>0,1,0)
GpresabsII<-ifelse(nsSpp3GSeasonII[,10:28]>0,1,0)
GpresabsIII<-ifelse(nsSpp3GSeasonIII[,10:28]>0,1,0)

# Convert to betapart objects
GpresabsI.core <- betapart.core(GpresabsI)
GpresabsII.core <- betapart.core(GpresabsII)
GpresabsIII.core <- betapart.core(GpresabsIII)

# Assign plot (without season) to each row
row.names(GpresabsI) <- paste(nsSpp3GSeasonI$plot_codeII, 1:nrow(GpresabsI), sep="")
row.names(GpresabsII) <- paste(nsSpp3GSeasonII$plot_codeII, 1:nrow(GpresabsII), sep="")
row.names(GpresabsIII) <- paste(nsSpp3GSeasonIII$plot_codeII, 1:nrow(GpresabsIII), sep="")


dimnames(GpresabsI) 
dimnames(GpresabsII)
dimnames(GpresabsIII)

# Beta.temp
btI<-beta.temp(GpresabsI,GpresabsII,index.family="sor") # this is soreson
btII<-beta.temp(GpresabsII,GpresabsIII,index.family="sor") 

nsSpp3GS2_3<-rbind(nsSpp3GSeasonII,nsSpp3GSeasonIII)
betaSeason<-rbind(btI,btII)
nsSpp3Gbeta<-cbind(nsSpp3GS2_3,betaSeason)
names(nsSpp3Gbeta)

# Analyse turnover component (sim) - GRASSES
turnGbt<-lmer(beta.sim~Livestockdensity+Treatment+Date+
             # Date:Treatment+Livestockdensity:Date+
#              Livestockdensity:Treatment+
              #Treatment:Livestockdensity:Date+#Simplified via LRT 
               (1 |Block), data=nsSpp3Gbeta)
summary(turnGbt)
anova(turnGbt) #
AIC(turnGbt) # -37.78991
drop1(turnGbt, test="Chi")

# Analyse nested component (sne) - GRASSES
nestGbt<-lmer(beta.sne~Livestockdensity+Treatment+#Date+
                 #Date:Treatment+Livestockdensity:Date+
                              Livestockdensity:Treatment+
               # Treatment:Livestockdensity:Date+#Simplified via LRT 
                (1 |Block), data=nsSpp3Gbeta)
summary(nestGbt)
anova(nestGbt) #
AIC(nestGbt) # -190.4819
drop1(nestGbt, test="Chi") # 
#Livestockdensity:Treatment 0.02193 *

sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
nsSpp3GbetaMean<-aggregate(beta.sne~Livestockdensity+Treatment,nsSpp3Gbeta,mean)
nsSpp3GbetaSD<-aggregate(beta.sne~Livestockdensity+Treatment,nsSpp3Gbeta,sem)
nsSpp3GbetaMean<-cbind(nsSpp3GbetaMean,nsSpp3GbetaSD[3])
colnames(nsSpp3GbetaMean)[4]<-"sd"
pd <- position_dodge(0.5)
ggplot(nsSpp3GbetaMean,aes(x=Livestockdensity,y=beta.sne,colour=Treatment))+
  geom_point(size=4.5, position=pd)+geom_errorbar(aes(x = Livestockdensity, ymin=beta.sne-sd,ymax=beta.sne+sd),
      position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)

# Plot Beta by livestock density and treatment
ggplot(nsSpp3Gbeta,aes(x=beta.sor,colour=Livestockdensity, linetype=Treatment))+geom_density()
ggplot(nsSpp3Gbeta,aes(x=beta.sim,colour=Livestockdensity, linetype=Treatment))+geom_density()
ggplot(nsSpp3Gbeta,aes(x=beta.sne,colour=Livestockdensity, linetype=Treatment))+geom_density()

# Comparison of the square root transformed bsim and bsne components of bsor b
with(nsSpp3Gbeta, plot(sqrt(beta.sim) ~ sqrt(beta.sne), 
                 type='n', ylab=expression(sqrt(beta[sim])), 
                 xlab=expression(sqrt(beta[sne]))))
with(nsSpp3Gbeta, text(y= sqrt(beta.sim), x=sqrt(beta.sne), col=c(Livestockdensity)))

# Nested - seems to be two plots that are outliers - are these low livestock?
plot(hclust(dist(nsSpp3Gbeta$beta.sim), method="single"),col=c(nsSpp3Gbeta$Livestockdensity),
     hang=-1, main='', sub='', xlab='')

##########################################################################################
#### Species richness #### 
##########################################################################################
#### Overall summary - species richness
nsSpp3$Richness<-apply(nsSpp3[,10:74]>0,1,sum)
aggregate(Richness~Livestockdensity+Treatment+Season,nsSpp3,mean)
nsSpp3$Shannon<-diversity(nsSpp3[,10:74], index="shannon")
ShanD1<-lme(Shannon~Treatment+Livestockdensity+Season+
              Treatment:Livestockdensity+Livestockdensity:Season+
              Treatment:Season+
              Treatment:Livestockdensity:Season,
            random= ~ 1|fBlock,data=nsSpp3)
summary(ShanD1)
anova(ShanD1)
AIC(ShanD1) #1136.222    
aggregate(Shannon~Livestockdensity+Season,nsSpp3,mean)
# Diversity higher in low livestock and third season

# Grasses
nsSpp3G$Richness<-apply(nsSpp3G[,10:28]>0,1,sum)
aggregate(Richness~Season+Livestockdensity,nsSpp3G,mean)






#### Biplot #####
# Biplot # Species - just grasses 
par(mfrow=c(1,1))
plot(mdsNSG, type="n",xlim=c(-2,2), ylim=c(-1,1),
     ylab="NMDS 2", xlab="NMDS 1",mgp=c(1.75,.45,0), 
     tck=.02, las=1, lwd=1.75, bty='l', main="Species")
with(nsSpp3G,text(mdsNSG,display="species",
                  col="grey", cex=.75))#Add spp
with(nsSpp3G,ordiellipse(mdsNSG,Livestockdensity:Treatment, conf=0.95,
                         cex=1.5, col="black", lwd=2, lty=c(1,2,3),label=F))

# USE NMDS on Dwarf shrub only
names(nsSpp3W0)
nsSpp3W$woodySUM<-rowSums(nsSpp3W[,c(10:19)]) # Need to remove rows with sum zero
nsSpp3W0<-droplevels(nsSpp3W[nsSpp3W$woodySUM!=0, ])
mdsNSW<-metaMDS(nsSpp3W0[,c(10:19)], trace =F) 
mdsNSW #  0.09276031

NS.evW <- envfit(mdsNSW~Livestockdensity+Treatment+ Season,data = nsSpp3W0, 
                 perm=999,strata=as.numeric(nsSpp3W0$Block))
NS.evW # Livestock density and treatment significant - not season

nsSpp3W0$harvest_code<-as.factor(with(nsSpp3W0, paste(Livestockdensity,Treatment,Date, sep="-")))

NS.harW <- envfit(mdsNSW~ harvest_code,data = nsSpp3W0, 
                  perm=999,strata=as.numeric(nsSpp3W0$Block))
NS.harW

# Extract centroids
#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258

vec.sp.dfW<-as.data.frame(cbind(NS.harW$factors$centroids*sqrt(NS.harW$factors$r)))
#vec.sp.dfG<-as.data.frame(scores(mdsNSG, display = "sites"))
vec.sp.dfW$Livestock.density<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                                "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.dfW$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                        "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.dfW$Season<-rep(c(1,3,2))

WReharvestBio<-aggregate(DwarfShrubNetReharvestBiomass1~harvest_code,nsSpp3W0,mean)
vec.sp.dfW$DwarfShrubNetReharvestBiomass1<-WReharvestBio[,2]
vec.sp.dfW$DwarfShrubNetReharvestBiomass1<-as.numeric(vec.sp.dfW$DwarfShrubNetReharvestBiomass1)
#NS.harG$scores
#spp.scrs <- as.data.frame(scores(mdsNSG, display = "species"))
#spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

# Plot woody
CenPlotW<-ggplot(vec.sp.dfW,aes(x=NMDS1,y=NMDS2))
#CenPlotG<-CenPlotG+geom_text(data=spp.scrs,aes(label=Species),colour="light grey")
CenPlotW<-CenPlotW+geom_point(aes( colour=Livestock.density,fill=Livestock.density,shape=Treatment,size=DwarfShrubNetReharvestBiomass1))
CenPlotW<-CenPlotW+geom_text(aes(label=Season),hjust=0, vjust=-.5)
CenPlotW<-CenPlotW+theme_classic()
CenPlotW




# Need to edit below...

# Exclosure species move away from open plots - just explore low intensity livestock

# Focus in on low intensity and season effects - grasses
nsSpp2GL<-nsSpp2G[nsSpp2G$Livestockdensity=="Low",]
nsSpp2GL$Trt_Liv<-as.factor(with(nsSpp2GL, paste(Livestockdensity,Treatment, sep="-")))
nsSpp2G$Trt_Liv<-as.factor(with(nsSpp2G, paste(Livestockdensity,Treatment, sep="-")))
nsSpp2$Trt_Liv<-as.factor(with(nsSpp2, paste(Livestockdensity,Treatment, sep="-")))

nsSpp2GLex<-nsSpp2GL[nsSpp2GL$Treatment=="Exclosure",]

nsSpp2GLex1<-nsSpp2GLex[nsSpp2GLex$Season=="Season I",]
nsSpp2GLex23<-nsSpp2GLex[nsSpp2GLex$Season!="Season I",]

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMean<- function(data) sapply(data, mean, na.rm = TRUE)
Mean1<-colMean(nsSpp2GLex1[,10:28])
Mean23<-colMean(nsSpp2GLex23[,10:28])
PlantMeans<-as.data.frame(cbind(Mean1,Mean23))
names(PlantMeans)
plot(Mean23~Mean1,PlantMeans)
abline(0, 1) # Loss of grass species - replacement by herbs?


# Focus in on low intensity and season effects - climbers

nsSpp2CL<-nsSpp2C[nsSpp2C$Livestockdensity=="Low",]
nsSpp2CL$Trt_Liv<-as.factor(with(nsSpp2CL, paste(Livestockdensity,Treatment, sep="-")))
nsSpp2C$Trt_Liv<-as.factor(with(nsSpp2C, paste(Livestockdensity,Treatment, sep="-")))

nsSpp2CLex<-nsSpp2CL[nsSpp2CL$Treatment=="Exclosure",]

nsSpp2CLex1<-nsSpp2CLex[nsSpp2CLex$Season=="Season I",]
nsSpp2CLex23<-nsSpp2CLex[nsSpp2CLex$Season!="Season I",]
names(nsSpp2CLex)

Mean1C<-colMean(nsSpp2CLex1[,10:13])
Mean23C<-colMean(nsSpp2CLex23[,10:13])
PlantMeansC<-as.data.frame(cbind(Mean1C,Mean23C))
plot(Mean23C~Mean1C,PlantMeansC)
abline(0, 1) # Loss of grass species - replacement by herbs?



#### Biplot #####
# Biplot # Species - just grasses 
par(mfrow=c(1,1))
plot(mdsNSG, type="n",xlim=c(-2,2), ylim=c(-2,2),
     ylab="NMDS 2", xlab="NMDS 1",mgp=c(1.75,.45,0), 
     tck=.02, las=1, lwd=1.75, bty='l', main="Species")
with(nsSpp2G,text(mdsNSG,display="species",
                  col="grey", cex=.75))#Add spp
with(nsSpp2G,ordiellipse(mdsNSG,Livestockdensity:Treatment, conf=0.95,
                         cex=1.5, col="black", lwd=2, lty=c(1,2,3),label=F))
names(nsSpp2G)

#CynodonnlemfuensisVanderyst
CynNle<-ggplot(nsSpp2GL,aes(x=Trt_Liv,y=CynodonnlemfuensisVanderyst,shape=Treatment,colour=Season))
CynNle<-CynNle+geom_boxplot(outlier.shape=NA,show.legend=F)+geom_jitter(size=2,stroke=1,show.legend=F)
CynNle

#IschaemumafrumJFGmeLDandy
HetCon<-ggplot(nsSpp2G,aes(x=Trt_Liv,y=IschaemumafrumJFGmeLDandy,shape=Treatment,colour=Season))
HetCon<-HetCon+geom_boxplot(outlier.shape=NA,show.legend=F)+geom_jitter(size=2,stroke=1,show.legend=F)
HetCon

# Heteropogoncontortus
names(nsSpp2GL)
HetCon<-ggplot(nsSpp2GL,aes(x=Trt_Liv,y=HeteropogoncontortusLRoemSchult,shape=Treatment,colour=Season))
HetCon<-HetCon+geom_boxplot(outlier.shape=NA,show.legend=F)+geom_jitter(size=2,stroke=1,show.legend=F)
HetCon



# Herb also increases...
#RhynchosiaminimaLDC
names(nsSpp2G)
RhyMin<-ggplot(nsSpp2,aes(x=Trt_Liv,y=RhynchosiaminimaLDC  ,shape=Treatment,colour=Livestockdensity))
RhyMin<-RhyMin+geom_boxplot(outlier.shape=NA,show.legend=F)+geom_jitter(size=2,stroke=1,show.legend=F)
RhyMin


#### END ####