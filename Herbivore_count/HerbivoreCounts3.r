#####################################################
# Nech Sar - Herbivore observation data
rm(list=ls())
library(MASS)
library(vegan)
library(ggplot2)
library(plyr)
library(rgeos)
library(sp)
library(raster)
library(rgdal)
library(dplyr)
library(Matrix)
library(ggplot2)
library(glmmADMB)
library(glmmTMB)
library(multcomp)
library(lme4)
library(INLA)

#####################################################
# All herbivore count data
#####################################################

# Herbivore counts and metabolic biomass
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/Herbivore_count")
HerbProb<-read.csv(file="AllHerbivoreCounts_3prob.csv", sep=",",header=TRUE)
names(HerbProb)
str(HerbProb)
head(HerbProb)
tail(HerbProb)

levels(HerbProb$Species)
#[1] "Burchells_Zebra"    "Cattle"             "Grants_Gazelle"     "Greater_Kudu"       "Pastoralist"       
#[6] "Swaynes_Hartebeest"

# Remove pastoralist
HerbProbTp<-droplevels(HerbProb[!HerbProb$Species=="Pastoralist",])

###############################################################################################
#### Number of animals #####
###############################################################################################
## Calculate metabolic biomass per herbivore ###
### Mean metabilic biomass per species ###

# Aggregate herbivore count per Zone1 per herbivore/pastoralist
names(HerbProbTp)

# Sum all counts per species
HerbProb2<-HerbProbTp%>% 
  group_by(Species,date) %>% 
  summarise(Total = sum(Total))

# Divide by total area of grassland km2
HerbProb2$No.km2<-HerbProb2$Total/44.5

aggregate(No.km2~Species, HerbProb2,mean)
#            Species     No.km2
#1    Burchells_Zebra 8.06741573
#2             Cattle 1.36329588
#3     Grants_Gazelle 1.00374532
#4       Greater_Kudu 0.32958801
#5 Swaynes_Hartebeest 0.05992509


#### Metabolic biomass ###
HerbProbTbz<-HerbProbTp[HerbProbTp$Species=="Burchells_Zebra",]
HerbProbTbz$Bio<-HerbProbTbz$Total*215
HerbProbTbz$MetBio<-HerbProbTbz$Bio^.75

HerbProbTc<-HerbProbT[HerbProbT$Species=="Cattle",]
HerbProbTc$Bio<-HerbProbTc$Total*250
HerbProbTc$MetBio<-HerbProbTc$Bio^.75

HerbProbTgg<-HerbProbT[HerbProbT$Species=="Grants_Gazelle",]
HerbProbTgg$Bio<-HerbProbTgg$Total*50
HerbProbTgg$MetBio<-HerbProbTgg$Bio^.75

HerbProbTgk<-HerbProbT[HerbProbT$Species=="Greater_Kudu",]
HerbProbTgk$Bio<-HerbProbTgk$Total*200
HerbProbTgk$MetBio<-HerbProbTgk$Bio^.75

HerbProbTsh<-HerbProbT[HerbProbT$Species=="Swaynes_Hartebeest",]
HerbProbTsh$Bio<-HerbProbTsh$Total*171
HerbProbTsh$MetBio<-HerbProbTsh$Bio^.75


####################################################################################

# Join to distance from boma
#ZoneDist<-read.csv(file="BomaZoneDistance.csv", sep=",",header=TRUE)
#ZoneDist2<-ZoneDist[,c("Zone","NEAR_DIST")]
#colnames(ZoneDist2)[1]<-"Zone1"
#HerbProbTp2<-left_join(HerbProbTp,ZoneDist2,by=c("Zone1"))

# Join zone area
#zonearea<-read.table("CountZonePolys3.txt",header=T,sep="\t")
#colnames(zonearea)[2]<-"Zone1"
#HerbProbTp3<-left_join(HerbProbTp2,zonearea,by=c("Zone1"))
#HerbProbTp3$Density<-HerbProbTp3$Total/HerbProbTp3$Area_km2
#HerbProbTp3$MetBioDensity<-HerbProbTp3$MetBio/HerbProbTp3$Area_km2
#HerbProbTp3$Dist.km<-HerbProbTp3$NEAR_DIST/1000
#ggplot(HerbProbTp3, aes(y=Total,x=NEAR_DIST, colour=Species, fill=Species))+geom_point()
#ggplot(HerbProbTp3, aes(y=Total,x=NEAR_DIST))+geom_point()+facet_wrap(~Species, ncol=5)
#ggplot(HerbProbTp3, aes(y=Density,x=NEAR_DIST))+geom_point()+facet_wrap(~Species, ncol=5)
#ggplot(HerbProbTp3, aes(y=MetBio,x=NEAR_DIST))+geom_point()+facet_wrap(~Species, ncol=5)
#ggplot(HerbProbTp3, aes(y=MetBioDensity,x=NEAR_DIST))+geom_point()+facet_wrap(~Species, ncol=5)

# Aggregate herbivore count per Zone1 per herbivore/pastoralist
#HerbProbTp5<-HerbProbTp3%>% 
#  group_by(Species,Zone1,NEAR_DIST, Area_km2) %>% 
#  summarise(Total = sum(Total))

#names(HerbProbTp5)
#HerbProbTp5$Density<-HerbProbTp5$Total/HerbProbTp5$Area_km2
#ggplot(HerbProbTp5, aes(y=Total,x=NEAR_DIST))+geom_point()
#BomaDist<-read.table(file="DistanceToBomaCell1.txt",header=T,sep="\t")
#ggplot(data=BomaDist, aes(x=XCoord, y=YCoord, colour=CORE))+geom_point()+theme_bw()

###############################################################################################
#### Number of individuals ####
###############################################################################################

# Remove Swaynes - only one occurence
HerbProbTp4<-droplevels(HerbProbTp3[!HerbProbTp3$Species=="Swaynes_Hartebeest",])
HerbProbZ<-droplevels(HerbProbTp3[HerbProbTp3$Species=="Burchells_Zebra",])


1/HerbProbTp4$Area_km2
M1<-glm.nb(Total~NEAR_DIST+Species+NEAR_DIST:Species+offset(1/Area_km2),data=HerbProbTp4)
M1b<-glm.nb(Total~NEAR_DIST+offset(1/Area_km2),data=HerbProbTp5)
summary(M1)
anova(M1)
anova(M1b)

# Dispersion statistic:
E2 <- resid(M1, type = "pearson")
N  <- nrow(HerbProbTp4)
p  <- length(coef(M1)) + 1  # '+1' is due to theta # Need to add plus extra parameter that is theta

sum(E2^2) / (N - p) # 1.35 # Overdispersed

library(pscl)
library(lmtest)

M4 <- zeroinfl(Total ~ NEAR_DIST+Species+NEAR_DIST:Species+offset(1/Area_km2) | 
                 NEAR_DIST+Species+NEAR_DIST:Species+offset(1/Area_km2), 
               dist = 'negbin',data = HerbProbTp4)

M4b <- zeroinfl(Total ~ NEAR_DIST+Species+offset(1/Area_km2) | 
                 NEAR_DIST+Species+offset(1/Area_km2), 
               dist = 'negbin',data = HerbProbTp4)

M4c <- zeroinfl(Total ~ Species+offset(1/Area_km2) | 
                  Species+offset(1/Area_km2), 
                dist = 'negbin',data = HerbProbTp4)

summary(M4)
lrtest(M4,M4b)
lrtest(M4b,M4c)

# Dispersion statistic:
E2 <- resid(M4, type = "pearson")
N  <- nrow(HerbProbTp4)
p  <- length(coef(M4)) + 1  # '+1' is due to theta # Need to add plus extra parameter that is theta

sum(E2^2) / (N - p) # 1.099759 # Very slight overdisp

drop1(M4, test="Chisq")



#Model interpretation of the NB GLM
#Sketch the fitted values

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = HerbProbTp4$NEAR_DIST, 
     y = HerbProbTp4$Total,
     xlab = "Distance",
     ylab = "Number of herbivores",
     cex.lab = 2,
     pch = 16, 
     type = "p")

# Create an artifical grid of covariate values
M1<-glm.nb(Total~NEAR_DIST+Species+NEAR_DIST:Species+offset(1/Area_km2),data=HerbProbTp4)
MyData <- expand.grid(NEAR_DIST = seq(from  = min(HerbProbTp4$NEAR_DIST),to= max(HerbProbTp4$NEAR_DIST),length = 25),
                      Area_km2 = seq(from  = min(HerbProbTp4$Area_km2),to= max(HerbProbTp4$Area_km2),length = 25),
                                     Species=levels(HerbProbTp4$Species))
M1<-glm.nb(Total~NEAR_DIST,data=HerbProbTp4)
#MyData <- expand.grid(NEAR_DIST = seq(from  = min(HerbProbTp4$NEAR_DIST),to= max(HerbProbTp4$NEAR_DIST),length = 25))
# Predict the expected squirrel values                                       
P1 <- predict(M1, newdata = MyData, type = "link", se = TRUE)
lines(x = MyData$NEAR_DIST, 
      y = exp(P1$fit), 
      lwd = 3)
lines(x = MyData$NEAR_DIST, 
      y = exp(P1$fit + 1.96 * P1$se.fit), 
      lwd = 3, 
      lty = 2)
lines(x = MyData$NEAR_DIST, 
      y = exp(P1$fit - 1.96 * P1$se.fit), 
      lwd = 3, 
      lty = 2)



HerbProbTp4$ZoneID<-as.numeric(HerbProbTp4$Zone1)
f1A <- Density ~ f(NEAR_DIST,model="rw2")+Species+NEAR_DIST:Species +f(Date, model="ar1", replicate=ZoneID)
f1B <- MetBioDensity~ NEAR_DIST+Species+NEAR_DIST:Species +f(Date, model="ar1", replicate=ZoneID)
RM1<-inla(f1A, data=HerbProbTp4, family="gaussian", # zeroinflatednbinomial0
        control.compute=list(dic=T)) #quantiles = c(0.25,0975))
RM2<-inla(f1B, data=HerbProbTp4, family="gaussian", # zeroinflatednbinomial0
          control.compute=list(dic=T)) #quantiles = c(0.25,0975))


summary(RM1)
bri.hyperpar.plot(RM1)
bri.fixed.plot(RM1)
bri.fixed.plot(RM2)

# Neg Bin dist on count data
M1 <- gamm(Density ~ NEAR_DIST+Species+NEAR_DIST:Species+Date,random=list(Zone1=~1),data = HerbProbTp4)
summary(M1)

(inla.models()$likelihood)

M1s <- gam(Density ~ NEAR_DIST,data = HerbProbTp4)
summary(M1s)



par(mfrow = c(1,2),mar = c(5,5,2,2)) # Nonlinear pattern distance
plot(M1, cex.lab = 1.5,  shade = TRUE)

names(HerbProbTp4)
HerbProbTp4$Dist.km<-HerbProbTp4$NEAR_DIST/1000
M2<-lmer(Density ~Dist.km+Species+#Dist.km:Species+
           (1|Zone1), #          (1|Zone1),
          HerbProbTp4)
summary(M2)
anova(M2)
drop1(M2, test="Chisq")
?lmer
# Dispersion statistic:
E2 <- resid(M2, type = "pearson")
N  <- nrow(HerbProbTp2)
p  <- length(coef(M2)) + 1  # '+1' is due to theta # Need to add plus extra parameter that is theta

sum(E2^2) / (N - p) # 0.894216
summary(M2)
anova(M2)

# Plot residuals vs fitted values
F3 <- fitted(M2)
E3 <- resid(M2, type = "pearson") # person residuals similar to 
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(x = F3, 
     y = E3,
     xlab = "Fitted values",
     ylab = "Pearson residuals",
     cex.lab = 1.5)
abline(h = 0, lty = 2)


###############################################################################################
#### Metabolic biomass #####
###############################################################################################

####################################################################################
## Calculate metabolic biomass per herbivore per zone per time ###
### Sum total metabilic biomass per zone ###
### Divide species metabolic biomass per zone by total metabolic biomass per zone ###

# Aggregate herbivore count per Zone1 per herbivore/pastoralist
HerbProbTp4<-HerbProbTp2%>% 
  group_by(Species,Zone1) %>% 
  summarise(Total = sum(Total))

HerbProbTbz<-HerbProbT[HerbProbT$Species=="Burchells_Zebra",]
HerbProbTbz$Bio<-HerbProbTbz$Total*215
HerbProbTbz$MetBio<-HerbProbTbz$Bio^.75

HerbProbTc<-HerbProbT[HerbProbT$Species=="Cattle",]
HerbProbTc$Bio<-HerbProbTc$Total*250
HerbProbTc$MetBio<-HerbProbTc$Bio^.75

HerbProbTgg<-HerbProbT[HerbProbT$Species=="Grants_Gazelle",]
HerbProbTgg$Bio<-HerbProbTgg$Total*50
HerbProbTgg$MetBio<-HerbProbTgg$Bio^.75

HerbProbTgk<-HerbProbT[HerbProbT$Species=="Greater_Kudu",]
HerbProbTgk$Bio<-HerbProbTgk$Total*200
HerbProbTgk$MetBio<-HerbProbTgk$Bio^.75

HerbProbTsh<-HerbProbT[HerbProbT$Species=="Swaynes_Hartebeest",]
HerbProbTsh$Bio<-HerbProbTsh$Total*171
HerbProbTsh$MetBio<-HerbProbTsh$Bio^.75

# Neg Bin dist on count data
M3 <- lm(MetBio ~ NEAR_DIST+Species+NEAR_DIST:Species+offset(1/Area_km2),data = HerbProbTp4)
summary(M3)
anova(M3)
drop1(M3,test="Chisq")

# Convert to binary data - presence vs absence and probability of detection
HerbProbTp$PresAb<-HerbProbTp$Total
HerbProbTp$PresAb[HerbProbTp$PresAb>0]<-1

# Binominal probablity
HPmod<- glm(PresAb~Species+Zone1+Species:Zone1,data=HerbProbTp,family=binomial())
summary(HPmod)
anova(HPmod)

# Binominal probability with mixed model with date 
HPmod<-glmmadmb(PresAb~Species+Zone1+Species:Zone1+
                  (1|date), 
                family="binomial",#zeroInflation = TRUE,
                #admb.opts=admbControl(shess=FALSE,noinit=FALSE),
                admb.opts=admbControl(shess=FALSE,noinit=FALSE, impSamp=200,maxfn=500,imaxfn=500,maxph=5),
                data=HerbProbTp)
summary(HPmod)

#Model matrix
MyData <- expand.grid(Species = levels(HerbProbTp$Species),
                      Zone1 = levels(HerbProbTp$Zone1))
head(MyData)

#Convert the covariate values into an X matrix
Xp <- model.matrix(~ Species+Zone1+Species:Zone1, data = MyData)

#Extract parameters and parameter covariance matrix
betas    <- fixef(HPmod)
Covbetas <- vcov(HPmod)

#Calculate the fitted values in the predictor scale
MyData$eta <- Xp %*% betas
MyData$Pi  <- exp(MyData$eta) / (1 + exp(MyData$eta))

#Calculate the SEs on the scale of the predictor function
MyData$se    <- sqrt(diag(Xp %*% Covbetas %*% t(Xp)))
MyData$SeUp  <- exp(MyData$eta + 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  + 1.96 *MyData$se))
MyData$SeLo  <- exp(MyData$eta - 1.96 *MyData$se) / 
  (1 + exp(MyData$eta  - 1.96 *MyData$se))

MyData


ggplot(MyData, aes(y=Pi, x=Zone1, colour =Species, fill=Species))+geom_point()+theme_classic()
MyData$WildLiv<-MyData$Species
levels(MyData$WildLiv)<-c("Wildlife","Cattle","Wildlife","Wildlife","Wildlife")
Mydata<-aggregate(Pi~WildLiv+Zone1,MyData, mean)
ggplot(Mydata, aes(y=V1, x=Zone1, colour=WildLiv, fill=WildLiv))+geom_hline(yintercept=.5)+
  geom_point(size=4)+theme_classic()

####################################################################################
## Calculate metabolic biomass per herbivore per zone per time ###
### Sum total metabilic biomass per zone ###
### Divide species metabolic biomass per zone by total metabolic biomass per zone ###

# Aggregate herbivore count per Zone1 per herbivore/pastoralist
HerbProbT<-HerbProb %>% 
  group_by(Species,Zone1) %>% 
  summarise(Total = sum(Total))

HerbProbTbz<-HerbProbT[HerbProbT$Species=="Burchells_Zebra",]
HerbProbTbz$Bio<-HerbProbTbz$Total*215
HerbProbTbz$MetBio<-HerbProbTbz$Bio^.75

HerbProbTc<-HerbProbT[HerbProbT$Species=="Cattle",]
HerbProbTc$Bio<-HerbProbTc$Total*250
HerbProbTc$MetBio<-HerbProbTc$Bio^.75

HerbProbTgg<-HerbProbT[HerbProbT$Species=="Grants_Gazelle",]
HerbProbTgg$Bio<-HerbProbTgg$Total*50
HerbProbTgg$MetBio<-HerbProbTgg$Bio^.75

HerbProbTgk<-HerbProbT[HerbProbT$Species=="Greater_Kudu",]
HerbProbTgk$Bio<-HerbProbTgk$Total*200
HerbProbTgk$MetBio<-HerbProbTgk$Bio^.75

HerbProbTsh<-HerbProbT[HerbProbT$Species=="Swaynes_Hartebeest",]
HerbProbTsh$Bio<-HerbProbTsh$Total*171
HerbProbTsh$MetBio<-HerbProbTsh$Bio^.75

TotalWildMetBio<-cbind(HerbProbTbz$MetBio,HerbProbTgg$MetBio,
      HerbProbTgk$MetBio,HerbProbTsh$MetBio)

HerbProbTbz$TotalWildMetBio<-rowSums(TotalWildMetBio)
HerbProbTc$TotalWildMetBio<-rowSums(TotalWildMetBio)
HerbProbTgg$TotalWildMetBio<-rowSums(TotalWildMetBio)
HerbProbTgk$TotalWildMetBio<-rowSums(TotalWildMetBio)
HerbProbTsh$TotalWildMetBio<-rowSums(TotalWildMetBio)

HerbProbMetBio<-rbind(HerbProbTbz,HerbProbTc,HerbProbTgg,HerbProbTgk,HerbProbTsh)
#HerbProbMetBio$PerBio<-(HerbProbMetBio$MetBio/HerbProbMetBio$TotalMetBio)
#HerbProbMetBio$PerBio[is.na(HerbProbMetBio$PerBio)]<-0

# Plot herbivore biomass across zones through time
names(HerbProbMetBio)
HBioT<-ggplot(HerbProbMetBio, aes(y=MetBio, x=Zone1))
HBioT<-HBioT+facet_wrap(~Species)
HBioT<-HBioT+geom_point()
HBioT<-HBioT+theme_classic()
HBioT

####################################################################################
dim(MyData) # 85
dim(HerbProbMetBio) # 85
head(MyData)
head(HerbProbMetBio)
HPMB<-left_join(HerbProbMetBio, MyData, by=c("Species","Zone1"))

names(HPMB)
HPMB$PiBio<-HPMB$Pi*HPMB$PerBio

Mydata2<-aggregate(TotalWildMetBio~WildLiv+Zone1,HPMB, mean)
Mydata$TotalWildMetBio<-Mydata2$TotalWildMetBio
ggplot(Mydata, aes(y=V1, x=Zone1, colour=WildLiv, fill=WildLiv, size=TotalWildMetBio))+geom_hline(yintercept=.5)+
  geom_point()+theme_classic()


#####################################################################################################
#### James old script #####
#####################################################################################################
#Herbivore count analysis after 3rd count

herbcounts3<-read.table("AllHerbivoreCounts_3.txt",header=T,sep="\t")

zonearea<-read.table("CountZonePolys3.txt",header=T,sep="\t")


herbcounts3<-merge(herbcounts3,zonearea,by.x="Zone1",by.y="Zone")
Density<-with(herbcounts3,Total/Area_km2)
herbcounts3<-cbind(herbcounts3,Density)
names(herbcounts3)

# Zone + date
herbcounts3$zone_date<-as.factor(with(herbcounts3, paste(Zone1,date, sep="-")))


with(herbcounts3,tapply(Density,list(Zone1,Species,Counting),mean,na.rm=T))
meandens<-with(herbcounts3,tapply(Density,list(zone_date,Species),sum))/4
meandens<-data.frame(cbind(rownames(meandens),meandens))
names(meandens)[1]<-"Zone"
zonearea<-zonearea[order(zonearea$Zone),]

meandensdf<-merge(zonearea,meandens,by.x="Zone",by.y="Zone",all.x=T,all.y=T)

dim(meandensdf)
head(meandensdf)
meandensdf[,2:12] <- as.numeric(as.matrix(meandensdf[,2:12]))
meandensdf[is.na(meandensdf)]<-0
#write.table(meandensdf,"M:\\Africa\\Herbivore Counts\\CountAverageFeb13.txt")
#write.table(meandensdf,"CountAverageFeb13_Date.txt")



herbcounts3a<-merge(aggregate(Density~(Zone1:Counting:Species),sum,data=herbcounts3),herbcounts3)
herbcounts3b<-reshape(herbcounts3a,direction="wide",timevar="Species",v.names=c("Total","Density"),idvar=c("Zone1","Counting"))
herbcounts3b[,20:37][is.na(herbcounts3b[,20:37])]<-0

with(herbcounts3b,plot(Xcent,Ycent,cex=Density.Burchells_Zebra))

(averageherbdens<-merge(aggregate(Density~(Zone1:Species),mean,data=herbcounts3a),herbcounts3a,all.x=T))


 with(meandensdf,plot(Burchells_Zebra,Cattle))
with(meandensdf,points(Burchells_Zebra,Greater_Kudu,col=grey(0.8),pch=16)) 
with(meandensdf,points(Burchells_Zebra,Grants_Gazelle,col=grey(0.4),pch=16)) 
 
   with(meandensdf, plot( Xcent,Ycent,cex=Burchells_Zebra,pch=1))
  with(meandensdf, points( Xcent,Ycent,cex=Cattle,col=2,pch=2))
    with(meandensdf, points( Xcent,Ycent,cex=Grants_Gazelle,col=3,pch=3))
      with(meandensdf, points( Xcent,Ycent,cex=Greater_Kudu,col=4,pch=4))
 
 barplot(as.matrix(t(meandensdf[,c(6:7,9:10,12)])),legend=T,names.arg=meandensdf$Zone)
 
 require(mapplots)
with(meandensdf,basemap(xlim=c(349500,354000),ylim=c(650000,660000),bg="white"))
with(meandensdf, draw.pie(Xcent,Ycent,as.matrix(meandensdf[,c(6:7,9:10,12)]),radius=100*rowSums(meandensdf[,c(6:7,9:10,12)])))
legend.pie("topr",labels=c("Zebra","Cattle","Gazelle","Kudu","Hartebeast"),radius=1000)

with(meandensdf,basemap(xlim=c(349500,354000),ylim=c(650000,660000),bg="white"))
with(meandensdf[8,], draw.pie(Xcent,Ycent,as.matrix(meandensdf[8,c(6:7,9:10,12)]),radius=50*rowSums(meandensdf[8,c(6:7,9:10,12)]),col=c("white",1,grey(0.6),grey(0.4))))
with(meandensdf[c(1:7,9:17),], draw.pie(Xcent,Ycent,as.matrix(meandensdf[c(1:7,9:17),c(6:7,9:10,12)]),radius=50*rowSums(meandensdf[c(1:7,9:17),c(6:7,9:10,12)]),col=c("white",1,grey(0.8),grey(0.6),grey(0.4))))
legend.pie("bottomr",labels=c("Zebra","Cattle","Gazelle","Kudu","Hartebeest"),radius=600,mab=3,init.angle=225,inset=-0.01,col=c("white",1,grey(0.8),grey(0.6),grey(0.4)))
text(353500,652700,"Legend",cex=1.2)
  
colSums(meandensdf[,5:12])/sum(meandensdf[,2])
colSums(meandensdf[,c(14,17,20,23,26)])/sum(meandensdf[,2])

#Zebra density  

library(maptools)
  myshp<-readShapePoly("M:\\Africa\\Plots\\CountZones.shp")
  plot(myshp,add=T)
  
  
#Met bio and TLUs

meandensdf$Burchells_ZebraBio<-meandensdf$Burchells_Zebra*216            #Biomass
meandensdf$Burchells_ZebraMetBio<-meandensdf$Burchells_ZebraBio^0.75   #Metabolic biomass =Biomass^0.75
meandensdf$Burchells_ZebraTLU<-meandensdf$Burchells_ZebraMetBio/250^0.75  #TLU = tropical livestock units (1 250kg Cow)

meandensdf$CattleBio<-meandensdf$Cattle*250            #Biomass
meandensdf$CattleMetBio<-meandensdf$CattleBio^0.75   #Metabolic biomass =Biomass^0.75
meandensdf$CattleTLU<-meandensdf$CattleMetBio/250^0.75  #TLU = tropical livestock units (1 250kg Cow)

meandensdf$Grants_GazelleBio<-meandensdf$Grants_Gazelle*50            #Biomass
meandensdf$Grants_GazelleMetBio<-meandensdf$Grants_GazelleBio^0.75   #Metabolic biomass =Biomass^0.75
meandensdf$Grants_GazelleTLU<-meandensdf$Grants_GazelleMetBio/250^0.75  #TLU = tropical livestock units (1 250kg Cow)

meandensdf$Greater_KuduBio<-meandensdf$Greater_Kudu*200            #Biomass
meandensdf$Greater_KuduMetBio<-meandensdf$Greater_KuduBio^0.75   #Metabolic biomass =Biomass^0.75
meandensdf$Greater_KuduTLU<-meandensdf$Greater_KuduMetBio/250^0.75  #TLU = tropical livestock units (1 250kg Cow)

meandensdf$Swaynes_HartebeestBio<-meandensdf$Swaynes_Hartebeest*171            #Biomass
meandensdf$Swaynes_HartebeestMetBio<-meandensdf$Swaynes_HartebeestBio^0.75   #Metabolic biomass =Biomass^0.75
meandensdf$Swaynes_HartebeestTLU<-meandensdf$Swaynes_HartebeestMetBio/250^0.75  #TLU = tropical livestock units (1 250kg Cow)

write.table(meandensdf,"M:\\Africa\\Herbivore Counts\\CountAverageFeb13.txt")


#MetBio Pie
with(meandensdf,basemap(xlim=c(349000,355500),ylim=c(649500,662000),bg="white"))
 plot(myshp,add=T)
with(meandensdf[8,], draw.pie(Xcent,Ycent,as.matrix(meandensdf[8,c(14,17,20,23,26)]),radius=rowSums(meandensdf[8,c(14,17,20,23,26)]),col=c("white",1,grey(0.6),grey(0.4))))
with(meandensdf[c(1:7,9:17),], draw.pie(Xcent,Ycent,as.matrix(meandensdf[c(1:7,9:17),c(14,17,20,23,26)]),radius=rowSums(meandensdf[c(1:7,9:17),c(14,17,20,23,26)]),col=c("white",1,grey(0.8),grey(0.6),grey(0.4))))
legend.pie("bottomr",labels=c("Zebra","Cattle","Gazelle","Kudu","Hartebeest"),radius=600,mab=3,init.angle=225,inset=-0.01,col=c("white",1,grey(0.8),grey(0.6),grey(0.4)))
text(353500,652700,"Legend",cex=1.2)
rowMeans(meandensdf[,c(14,17,20,23,26)])
  
#TLU Pie
with(meandensdf,basemap(xlim=c(340000,360000),ylim=c(640000,670000),bg="white"))
with(meandensdf[8,], draw.pie(Xcent,Ycent,as.matrix(meandensdf[8,c(15,18,21,22,27)]),radius=5*rowSums(meandensdf[8,c(15,18,21,22,27)]),col=c("white",1,grey(0.6),grey(0.4))))
with(meandensdf[c(1:7,9:17),], draw.pie(Xcent,Ycent,as.matrix(meandensdf[c(1:7,9:17),c(15,18,21,22,27)]),radius=5*rowSums(meandensdf[c(1:7,9:17),c(15,18,21,22,27)]),col=c("white",1,grey(0.8),grey(0.6),grey(0.4))))
legend.pie("bottomr",labels=c("Zebra","Cattle","Gazelle","Kudu","Hartebeest"),radius=600,mab=3,init.angle=225,inset=-0.01,col=c("white",1,grey(0.8),grey(0.6),grey(0.4)))
text(353500,652700,"Legend",cex=1.2)

#Mean species densities
 apply(meandensdf[,c(6:7,9:10,12)],2,mean)
 apply(meandensdf[,c(6:7,9:10,12)],2,sem)
#Mean species Met bio
 apply(meandensdf[,c(14,17,20,23,26)],2,mean)
 apply(meandensdf[,c(14,17,20,23,26)],2,sem)

