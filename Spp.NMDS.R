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
############################################################################

nsSpp<-read.csv("VegCompSeason3.csv",header=T,sep=",")

names(nsSpp)
head(nsSpp)
str(nsSpp)

################################################################################
#### Non multi-dimensional scaling - explore spp data ####
################################################################################

# USE NMDS
names(nsSpp)
names(nsSpp[,8:72]) # 65 species # Rosett?
mdsNS<-metaMDS(nsSpp[,8:72], trace =F) 
mdsNS
#Stress:     0.2925064  # Not great stress

# Stressplot
vare.dis<-vegdist(nsSpp[,8:72]) 
vare.dis2<-as.matrix(vare.dis)
vare.dis
vare.mds0<-isoMDS(vare.dis) 
stressplot(vare.mds0,vare.dis) # final  value 23.725340 
# Metri = 0.94, LinearR2 =0.77

# Basic NMDS plot
par(mfrow=c(1,1))
plot(mdsNS, type="p")
plot(mdsNS$points[,1]~mdsNS$points[,2]) 
# No major outliers - strong central clustering

# Enviro fit - fit season, treatment, livestock density 
names(nsSpp)
NS.ev <- envfit(mdsNS~ Season.Year+Treatment+Livestock.density,data = nsSpp, 
                  perm=999) #strata=as.numeric(nsSpp$Transect)/as.numeric(nsSpp$Quadrat))
NS.ev # All seemily important for spp composition - low R2

# Issue transect (experimental design) corresponds to dung

# Biplot # Species
par(mfrow=c(1,1))
plot(mdsNS, type="n",xlim=c(-2,2), ylim=c(-2,2),
     ylab="NMDS 2", xlab="NMDS 1",mgp=c(1.75,.45,0), 
     tck=.02, las=1, lwd=1.75, bty='l', main="Species")
with(nsSpp,text(mdsNS,display="species",
                         col="grey", cex=.75))#Add spp
with(nsSpp,ordiellipse(mdsNS,Livestock.density, conf=0.95,
                                cex=1.5, col="black", lwd=2, lty=c(1,2,3),label=T))

# Wider range of species in low livestock category - but smaller when wide
# Categorical division of livestock, but was it the relationship as raw dung/metabolic bio