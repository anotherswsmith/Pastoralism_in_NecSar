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
levels(nsreharvest3$Season)<-c("Season III" ,  "Season II" , "Season I")
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
dim(nsSpp3)

################################################################################
#### Non multi-dimensional scaling - explore spp data ####
################################################################################

# USE NMDS
names(nsSpp3)
names(nsSpp3[,10:74]) # 65 species # Rosett?
mdsNS<-metaMDS(nsSpp3[,10:74], trace =F) 
mdsNS
#Stress:    0.2910275   # Not great stress - but under 0.3

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
nsSpp3$harvest_code<-as.factor(with(nsSpp3, paste(Livestock.density,Treatment,Date, sep="-")))

# Envfit
NS.ev <- envfit(mdsNS~Livestock.density+Treatment+ Season,data = nsSpp3, 
                  perm=999,strata=as.numeric(nsSpp3$Block))
NS.ev 
# Weak difference in season - biggest difference is livestock density - r2=0.26

NS.har <- envfit(mdsNS~ harvest_code,data = nsSpp3, 
                perm=999,strata=as.numeric(nsSpp3$Block))
NS.har

# Extract centroids
vec.sp.df<-as.data.frame(cbind(NS.har$factors$centroids*sqrt(NS.har$factors$r)))
vec.sp.df$Livestock.density<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                               "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.df$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                       "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.df$Season<-rep(c(1,3,2))

names(nsSpp3)
TotBio<-aggregate(TotalBiomass1~harvest_code,nsSpp3,mean)
vec.sp.df$TotalBiomass1<-TotBio[,2]
vec.sp.df$TotalBiomass1<-as.numeric(vec.sp.df$TotalBiomass1)


# Plot centroids
CenPlot<-ggplot(vec.sp.df,aes(x=NMDS1,y=NMDS2, colour=Livestock.density,fill=Livestock.density,shape=Treatment))
CenPlot<-CenPlot+geom_point(aes(size=TotalBiomass1))
CenPlot<-CenPlot +geom_text(aes(label=Season),hjust=0, vjust=-.5)
CenPlot<-CenPlot +ggtitle("Total biomass")
CenPlot<-CenPlot +theme_classic()
CenPlot

# Repeat analysis - but just for grasses, woody, herbs and climbers
nsSppFx<-read.csv("Spp.Fx.group.csv",header=T,sep=",")
nsSppFxG<-droplevels(nsSppFx[nsSppFx$Fx.group!="Grass",])
nsSppFxW<-droplevels(nsSppFx[nsSppFx$Fx.group!="Dwarf shrub",])
nsSppFxC<-droplevels(nsSppFx[nsSppFx$Fx.group!="Climber",])

# Remove . and spaces = same names..Grasses
names(nsSpp3)  <- gsub("\\.", "",  names(nsSpp3))
to.remove  <- gsub(" ", "", levels(nsSppFxG$Species))
to.remove  <- gsub("\\.", "", to.remove )
to.remove  <- gsub("[()]", "",to.remove)

# Remove . and spaces = same names..Woody
to.removeW  <- gsub(" ", "", levels(nsSppFxW$Species))
to.removeW  <- gsub("\\.", "", to.removeW )
to.removeW  <- gsub("[()]", "",to.removeW)

# Remove . and spaces = same names..Climbers
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

# Subset climbers
nsSpp3C<-nsSpp3[ , !(names(nsSpp3) %in% to.removeC)]
nsSpp3C<-subset(nsSpp3,select = names(nsSpp3) %ni% to.removeC)
names(nsSpp3C) # Subset is now just climber

# USE NMDS on Grass only
names(nsSpp3G[,10:28]) # 
mdsNSG<-metaMDS(nsSpp3G[,10:28], trace =F) 
mdsNSG # 0.1550549

names(nsSpp3G)
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
vec.sp.dfG$Livestock.density<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                               "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.dfG$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                       "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.dfG$Season<-rep(c(1,3,2))

GReharvestBio<-aggregate(GrassNetReharvestBiomass1~harvest_code,nsSpp3G,mean)
vec.sp.dfG$GrassNetReharvestBiomass1<-GReharvestBio[,2]
vec.sp.dfG$GrassNetReharvestBiomass1<-as.numeric(vec.sp.dfG$GrassNetReharvestBiomass1)
#NS.harG$scores
#spp.scrs <- as.data.frame(scores(mdsNSG, display = "species"))
#spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

# Plot
CenPlotG<-ggplot(vec.sp.dfG,aes(x=NMDS1,y=NMDS2))
#CenPlotG<-CenPlotG+geom_text(data=spp.scrs,aes(label=Species),colour="light grey")
CenPlotG<-CenPlotG+geom_point(aes( colour=Livestock.density,fill=Livestock.density,shape=Treatment,size=GrassNetReharvestBiomass1))
CenPlotG<-CenPlotG+geom_text(aes(label=Season),hjust=0, vjust=-.5)
CenPlotG<-CenPlotG+ggtitle("Grass")
CenPlotG<-CenPlotG+theme_classic()
CenPlotG

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