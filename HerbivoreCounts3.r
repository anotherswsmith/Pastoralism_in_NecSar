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

#####################################################
# All herbivore count dat
#####################################################

# Herbivore probability
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

