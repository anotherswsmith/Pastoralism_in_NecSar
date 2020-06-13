########################################################################
#### Nech Sar - Grassland Biomass regrowth and pastoralism ####
rm(list=ls())
library(MASS)
library(vegan)
library(ggplot2)
library(plyr)
library(rgeos)
library(sp)
library(raster)
library(rgdal)
library(spatstat)
library(maptools)
########################################################################
#### Exclosure biomass ####
#########################################################################
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/Exclosures_spp_regrowth/")
nsbiomass3<-read.table("BiomassSeason3a.txt",header=T,sep="\t")

# Structure of data
names(nsbiomass3)
head(nsbiomass3)
tail(nsbiomass3)
str(nsbiomass3)

# Standard error function
sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# Factors 
#nsbiomass3$Livestock.density<-factor(nsbiomass3$Livestock.density,levels=c("Low","Medium","High"),ordered=T)
nsbiomass3$Boma.density<-factor(nsbiomass3$Boma.density,levels=c("High","Low"),ordered=T)
nsbiomass3$HerbClimb1<-nsbiomass3$HerbNetBiomassSeason1+nsbiomass3$ClimberNetBiomassSeason1
nsbiomass3$HerbClimb2<-nsbiomass3$HerbNetBiomassSeason2+nsbiomass3$ClimberNetBiomassSeason2
nsbiomass3$HerbClimb3<-nsbiomass3$HerbNetBiomassSeason3+nsbiomass3$ClimberNetBiomassSeason3
nsbiomass3$plotid<-paste(nsbiomass3$Treatment,nsbiomass3$Livestock.density)

#Convert to gm^-2
#nsbiomass3[,6:23]<-nsbiomass3[,6:23]*4 # Original - changed order?
names(nsbiomass3)
names(nsbiomass3[,c(8:22,25:27)])
nsbiomass3[,c(8:22,25:27)]<-nsbiomass3[,c(8:22,25:27)]*4
levels(nsbiomass3$Treatment)
nsbiomass3orig<- droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",])
with(nsbiomass3orig,tapply(TotalSeason3,list(Treatment,Boma.density),mean))

########################################################################
##### Combine dung and biomass regrowth exclosure xy coordinates ####
########################################################################
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/Exclosures_spp_regrowth")
nsreharvest3<-read.table("ProductivitySeason3.txt",header=T,sep="\t")
names(nsreharvest3)

setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/GIS")
Ex_location<-read.csv(file="PlotNames.csv", sep=",",header=TRUE)
names(Ex_location)
levels(Ex_location$OtherName)
colnames(Ex_location)[7]<-"Trt.name"

with(plot(Y~X,Ex_location))

#### Join regrowth and dung data ####
# Regrowth biomass
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/Herbivore_count")
nsherb3<-read.table("CountAverageFeb13_Date.txt",header=T,sep="\t")
names(nsherb3)

myvars <- c("Trt.name", "X", "Y")
Ex_location1<-Ex_location[myvars]
nsreharvest3loc<-merge(nsreharvest3,Ex_location1, by="Trt.name",drop=F)

# Spatial join - herbivore observations and regrowth
utmproj<-"+proj=utm +north +zone=37 +init=EPSG:32637" # central-ethiopia-37n
wp222loc<- cbind(nsreharvest3loc$X,nsreharvest3loc$Y)# get geographical location
Biosp<-SpatialPointsDataFrame(wp222loc,nsreharvest3loc, proj4string=CRS(utmproj),match.ID = TRUE)
Herbloc<-cbind(nsherb3$Xcent,nsherb3$Ycent)# get geographical location
Herbsp<-SpatialPointsDataFrame(Herbloc,nsherb3, proj4string=CRS(utmproj),match.ID = TRUE)

# Nearest number code
nsreharvest3loc$nearest_in_set2 <- apply(gDistance(Herbsp,Biosp, byid=TRUE), 1, which.min)
nsherb3$nearest_in_set2<-seq(1:17)
dim(Herbsp) # 68
dim(Biosp) # 630
dim(nsreharvest3loc) #630

#### Pastoral settlements - surveyed by Kjirsten ####
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/GIS")
bomas<-readOGR(dsn=".", "boma_points")
bomas_proj <-bomas
latlongproj<-("+proj=longlat +datum=WGS84")
utmproj<-"+proj=utm +north +zone=37 +init=EPSG:32637"
coordinates(bomas_proj)<- ~coords.x1 + coords.x2
proj4string(bomas_proj)<-latlongproj
bomas_projutm<-spTransform(bomas_proj,utmproj)
extent(bomas_projutm)
dim(bomas_proj)

# Exclosure locations
excl_proj <-nsreharvest3loc
coordinates(excl_proj)<- ~ X + Y
proj4string(excl_proj)<-utmproj
excl_projLAT<-spTransform(excl_proj,latlongproj)
extent(excl_projLAT)
par(mfrow=c(1,1))
plot(bomas_projutm,pch=16)
plot(excl_proj, add=T,pch=16, col="red")
BOMAScord<-coordinates(bomas_projutm)
BOMAScord<-BOMAScord[,1:2]
EXCLcord<-coordinates(excl_proj)
EXCLcord<-EXCLcord[,1:2]

#### Nech Sar outline ####
NechSar<-readOGR(dsn=".", "NSNPboundary")
NechSar_proj <-NechSar
proj4string(NechSar_proj)<-utmproj
area(NechSar) #391202776 = 514 km2
514/391202776
391202776*1.313897e-06

#### Grassland outline ####
NechSarGrassland<-readOGR(dsn=".", "HerbivoreCountZonest")
plot(NechSarGrassland)

#### Smooth buffer around grassland ####
NechSarBuff<-buffer(NechSarGrassland, width=175, dissolve=T)
NechSarBuff2<-gBuffer(NechSarGrassland, width = 174, byid = T)
BuffDiff<- gDifference(NechSarBuff, NechSarBuff2)
plot(BuffDiff)
plot(NechSarGrassland, add=T, col="red") # Adds 175 buffer
NechSarGrassland<-BuffDiff
NechSarGrassland_proj<-NechSarGrassland
proj4string(NechSarGrassland_proj)<-utmproj
plot(NechSarGrassland)
area(NechSarGrassland) #33914.14
33914.14*1.313897e-03#grassland is 44.5 km2

### Water features ####
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Mapping/GLWD-level1")
Lakes2<- readOGR(dsn=".", "glwd_1")
Lakes_proj<-Lakes2
proj4string(Lakes_proj)<-latlongproj
LakesLAT<-spTransform(NechSar_proj,latlongproj)
bb<-extent(LakesLAT)
Lake<-crop(Lakes_proj,LakesLAT)
Lakes<-spTransform(Lake,utmproj)

### Rivers ####
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/GIS")
Rivers<-readOGR(dsn=".","chamo_rivers") 
Rivers_proj<-Rivers
proj4string(Rivers_proj)<-utmproj
River<-crop(Rivers_proj,NechSar_proj)
plot(River)

#### Interplote  ####
NechSarOWIN<-as.owin(NechSar)
NechSarGrasslandOWIN<-as.owin(NechSarGrassland)
class(NechSarOWIN)
plot(NechSarOWIN)
plot(NechSarGrassland,col="green",add=T)
plot(Lakes,col="dodgerblue",add=T)
plot(River,col="dodgerblue",lwd=2,add=T)
plot(bomas_projutm, pch=21,add=T)
plot(excl_proj, add=T,pch=16, col="red")
pts<-coordinates(bomas_projutm)[,1:2]
head(pts)
extent(NechSar)

#Bounding box
coords = matrix(c(34.07309 , -3.426145 ,
                  34.07309 , -1.629324,
                  35.6555, -1.629324,
                  35.6555, -3.426145 ,
                  34.07309 , -3.426145 ), 
                ncol = 2, byrow = TRUE)

coords = matrix(c(338260.1 , 645081.5,
                  338260.1 , 673426.7,
                  365575.2, 673426.7,
                  365575.2, 645081.5,
                  338260.1 , 645081.5), 
                ncol = 2, byrow = TRUE)
P1 <- Polygon(coords)
bb <- SpatialPolygons(list(Polygons(list(P1), ID = "a")), proj4string=CRS("+proj=utm +north +zone=37 +init=EPSG:32637"))
NechSarBB<-as.owin(bb)

####ppp - point pattern ####
p<-ppp(pts[,1],pts[,2], window=NechSarBB) # 30 points outside boundary
plot(p)

#### Kernel density pastoral settlements - Kmeans ####
ds<-density(p,kernel = c("gaussian"))
plot(ds, main ="boma density")
plot(NechSarOWIN, add=T)

plot(excl_proj, add=T)
plot(bomas_projutm, add=T, pch=16)
format(4*10^-6,scientific = FALSE)

# Nearest distance settlements to exclosures
dist<-gDistance(excl_proj,bomas_projutm,byid=T)
min_distExclosures<-apply(dist, 2, min)
nsreharvest3loc$min_distExclosures<-min_distExclosures
# Distance not great / does not account for density of points

#### Intersect density plot ####
pts2<-coordinates(excl_proj)[,1:2]
p2<-ppp(pts2[,1],pts2[,2], window=NechSarOWIN)
nsreharvest3loc$boma_density<-ds[p2,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(p2)]

# Create random points and display blocks in probability distibution of points on map
pts <- spsample(NechSar, 9999, type = 'random')
pRand<-coordinates(pts)[,1:2]
pRand2<-ppp(pRand[,1],pRand[,2], window=NechSarOWIN)
BomaRandom<-ds[pRand2,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(pRand2)]

ptsG <- spsample(NechSarGrassland, 9999, type = 'random')
pRandG<-coordinates(ptsG)[,1:2]
pRandG2<-ppp(pRandG[,1],pRandG[,2], window=NechSarGrasslandOWIN)
BomaRandomGrassland<-ds[pRandG2,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(pRandG2)]

plot(pRandG2)
plot(NechSarGrassland_proj,add=T)
plot(bomas_projutm, pch=21,add=T)
plot(excl_proj, add=T,pch=16, col="red")

#### Settlement density - Nech Sar ####
x <- BomaRandom
df <- approxfun(density(x))
plot(density(BomaRandom))
xnew <- aggregate(boma_density~Plot.pair,nsreharvest3loc, mean)
#xnew<-c(0.00001,0.00002,0.00003)
points(xnew$boma_density,df(xnew$boma_density),pch=21,bg="red",col="dark red", cex=1.5)

#### Settlement density - Nech Sar Grassland ####
xG <- BomaRandomGrassland
dfG <- approxfun(density(xG))
plot(density(BomaRandomGrassland))
xnewG <- aggregate(boma_density~Plot.pair,nsreharvest3loc, mean)
#xnew<-c(0.00001,0.00002,0.00003)
points(xnewG$boma_density,dfG(xnewG$boma_density),pch=21,bg="red",col="dark red", cex=1.5)

# Testing alternative kernel densities
#ds<-density(p,kernel = c("gaussian"))
#dsE<-density(p,kernel = c("epanechnikov"))
#BomaRandom<-ds[pRand2,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(pRand2)]
#BomaRandomE<-dsE[pRand2,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(pRand2)]
#plot(BomaRandom~BomaRandomE)

#### Density based scanning ####
#library(factoextra)
#library(fpc)
#library(dbscan)
#Tukulpts<-coordinates(bomas_proj)[,1:2]
#km.res <- kmeans(Tukulpts, 2, nstart = 25)
#fviz_cluster(km.res, Tukulpts, frame = FALSE, geom = "point")
#db <- fpc::dbscan(Tukulpts, eps = .15, MinPts = 5)

# Plot DBSCAN results
#plot(db, Tukulpts, main = "DBSCAN", frame = FALSE) # One big cluster...
#fviz_cluster(db, Tukulpts, ellipse = 0.05, frame = FALSE, geom = "point")
#dbscan::kNNdistplot(Tukulpts, k =  5)

#distH<-gDistance(excl_proj,bomas_projutm,byid=T)
#min_distHerbs<-apply(dist, 2, min)

nsreharvest3loc$nearest_in_set2<-min_constructionDistance
plot(TotalBiomass1~nearest_in_set2,nsreharvest3loc)
plot(TotalBiomass1~boma_density,nsreharvest3loc)

#### Join herbivore biomass ####
MyHerb<-c("nearest_in_set2","Date","Burchells_ZebraMetBio","CattleMetBio","Grants_GazelleMetBio","Greater_KuduMetBio","Swaynes_HartebeestMetBio")
nsherb3sub<-nsherb3[MyHerb]

nsreharvest3b<-droplevels(nsreharvest3[nsreharvest3$Treatment=="Control" | nsreharvest3$Treatment=="Exclosure",])
dim(nsreharvest3)

#### Recheck plot pairs ####
nsbiomass3herbio<-left_join(nsreharvest3b,nsreharvest3loc, by=c("Plot.name","Reharvest.date","Harvest","rain.mm"),drop=T)
 #merge(nsreharvest3,nsreharvest3loc,by.x='Plot.name',by.y='Plot.name',all.x=T,all.y=F)
nsbiomass3herbio$fX<-as.factor(nsbiomass3herbio$X)
levels(nsbiomass3herbio$fX)<-c("B1","B2","B3","B4","B5","B6","B7","B8","B9")
#write.table(nsbiomass3herbio, "nsbiomass3herbio.txt", row.name=F, sep="\t")
dim(nsbiomass3herbio)

# Regrowth and herbivore metabolic biomass
levels(nsreharvest3loc$Reharvest.date) #"21.11.2013" "21.6.2013"  "25.11.2012"
levels(nsherb3sub$Date)#25.3.2012 26.10.2012 8.2.2013 9.2.2013
levels(nsherb3sub$Date)<-c("21.6.2013","21.6.2013" ,"21.11.2013", "21.11.2013")
nsherb3sub$Reharvest.date<-nsherb3sub$Date
nsherb3sub<-aggregate(.~Reharvest.date+nearest_in_set2,nsherb3sub,mean)

library(dplyr)
nsreharvest3locHerb<-left_join(nsreharvest3loc,nsherb3sub, by=c("Reharvest.date","nearest_in_set2"),drop=F)
nsReharvest<-nsreharvest3locHerb

#### Alternative - create krigging maps of herbivores ####
#### Create krigging maps of herbivore metabolic biomass ####
coordinates(nsherb3) <- ~ Xcent + Ycent
nsherb3.grid<-nsreharvest3locHerb
coordinates(nsherb3.grid) <- ~ X + Y

table(nsReharvestH1$Block,nsReharvestH1$Plot.name)
table(nsReharvestH1$Block,nsReharvestH1$Treatment,nsReharvestH1$Reharvest.date) #,nsReharvestH1$Transect)
write.table(nsReharvestH1, "nsReharvestH1.txt", row.name=F, sep="\t")
boxplot(boma_density~Block,nsReharvestb)

# Cattle variogram
#cattle.vgm <- variogram(CattleMetBio~1, nsherb3)
#cattlef = function(x) attr(cattle.fit <<- fit.variogram(cattle.vgm, vgm(,"Mat",nugget=NA,kappa=x)),"SSErr")
#optimize(cattlef, c(0.1, 5))
#plot(cattle.vgm, cattle.fit)
#cattle.kriged <- krige(CattleMetBio ~ 1, nsherb3, nsherb3.grid, model=cattle.fit)

#cattle.kriged %>% as.data.frame %>%
 #ggplot(aes(x=X, y=Y)) + geom_tile(aes(fill=var1.pred), height=500,width=500) + coord_equal() +
 # scale_fill_gradient(low = "yellow", high="red") +
#  theme_bw()
# Largerly under-estimates biomass values of cattle - 
# Max value ~250...due to plateau in semi-variogram

#nsreharvest3locHerb$cat.krig<-cattle.kriged$var1.pred
#names(nsreharvest3locHerb)
#plot(nsreharvest3locHerb$cat.krig,nsreharvest3locHerb$CattleMetBio)
# No much difference here # Krig = lower estimate of cattle 

#### Boma density - plot pair ####
# 9 plot pairs in experiment 
ggplot(nsreharvest3b,aes(y=boma_density,x=Plot.pair, colour=Block))+
  geom_boxplot(size=2)+theme_classic()#+coord_flip()


#### Grass and tin roof settlements - surveyed by Aramde Fetene et al. 2019 ####
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/GIS/Fetene_settlements")
tinroof<-readOGR(dsn=".", "Tin_roof")
grassroof<-readOGR(dsn=".", "Grass_roof")
tinroof_proj <-tinroof
grassroof_proj <-grassroof
latlongproj<-("+proj=longlat +datum=WGS84")
utmproj<-"+proj=utm +north +zone=37 +init=EPSG:32637"
coordinates(tinroof_proj)<- ~coords.x1 + coords.x2
coordinates(grassroof_proj)<- ~coords.x1 + coords.x2
proj4string(tinroof_proj)<-latlongproj
proj4string(grassroof_proj)<-latlongproj
tinroof_projutm<-spTransform(tinroof_proj,utmproj)
grassroof_projutm<-spTransform(grassroof_proj,utmproj)


plot(bomas_projutm,pch=16)
plot(tinroof_projutm, add=T,pch=16, col="blue")
plot(grassroof_projutm, add=T,pch=16, col="green")
plot(excl_proj, add=T,pch=16, col="red")

colnames(grassroof_projutm@data)[4]<-"settlements"
colnames(tinroof_projutm@data)[4]<-"settlements"

FenteSettlements<-rbind(grassroof_projutm,tinroof_projutm)

#### Kernel density pastoral settlements - Kmeans ####
ptsF<-coordinates(FenteSettlements)[,1:2]
pF<-ppp(ptsF[,1],ptsF[,2], window=NechSarBB) # 30 points outside boundary
plot(pF)

dsF<-density(pF,kernel = c("gaussian"))
plot(dsF, main ="boma density")
plot(NechSarOWIN, add=T)
plot(NechSarGrassland_proj,add=T)
plot(FenteSettlements, pch=21,add=T)
plot(excl_proj, add=T,pch=16, col="red")

pts2<-coordinates(excl_proj)[,1:2]
p2<-ppp(pts2[,1],pts2[,2], window=NechSarOWIN)
nsreharvest3loc$boma_densityFetene<-dsF[p2,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(p2)]

names(nsreharvest3loc)
distinct_ns <-distinct(nsreharvest3loc , Y, .keep_all = TRUE)

# TEST RANDOM SETTLEMENT POINTS IN GRASSLAND
dim(bomas_proj) # 679
ptsG679 <- spsample(NechSarGrassland, 9999, type = 'random')
pRandG679<-coordinates(ptsG679)[,1:2]
pRandG679b<-ppp(pRandG679[,1],pRandG679[,2], window=NechSarGrasslandOWIN)
BomaRandomGrassland679<-ds[pRandG679b,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(pRandG679b)]
BomaRandomGrassland679F<-dsF[pRandG679b,drop=TRUE, tight=FALSE, raster=NULL, rescue=is.owin(pRandG679b)]

str(BomaRandomGrassland679)
plot(BomaRandomGrassland679~BomaRandomGrassland679F)
summary(lm(BomaRandomGrassland679~BomaRandomGrassland679F))
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             2.653e-07  1.540e-08   17.23   <2e-16 ***
#  BomaRandomGrassland679F 7.480e+00  5.278e-02  141.72   <2e-16 ***
#Residual standard error: 2.776e-07 on 677 degrees of freedom
#Multiple R-squared:  0.9674,	Adjusted R-squared:  0.9673 
#F-statistic: 2.008e+04 on 1 and 677 DF,  p-value: < 2.2e-16


# TEST EXCLOSURE LOCATIONS
plot(distinct_ns$boma_densityFetene~distinct_ns$boma_density)
abline(lm(boma_densityFetene~boma_density,distinct_ns))
summary(lm(boma_densityFetene~boma_density,distinct_ns))
# F 528.4, df=1,13, r2=0.97, p<0.001
RankSettlement<-wilcox.test(distinct_ns$boma_densityFetene,distinct_ns$boma_density, paired = TRUE)
#RankSettlement

######################################################################
#### Biomass + Regrowth - H1 and H2 ####
########################################################################
names(nsReharvest)
# Set up as date
nsReharvest$Harvest.date<-as.Date(nsReharvest$Harvest.date,"%d.%m.%Y")
nsReharvest$Reharvest.date<-as.Date(nsReharvest$Reharvest.date,"%d.%m.%Y")
levels(as.factor(nsReharvest$Reharvest.date))
# Only need exclosure and control
levels(nsReharvest$Treatment)
nsReharvestb<-droplevels(nsReharvest[nsReharvest$Treatment=="Control" | nsReharvest$Treatment=="Exclosure",])

# Remove single reharvested quadrats from sampling
dim(nsReharvestb)
nsReharvestb2<- droplevels(nsReharvestb[nsReharvestb$Harvest!="single",])
nsReharvestSingle0<- droplevels(nsReharvestb[nsReharvestb$Harvest!="single" | nsReharvestb$Reharvest.date!="2013-06-21",])
nsReharvestSingle1<- droplevels(nsReharvestSingle0[nsReharvestSingle0$Harvest!="double" | nsReharvestSingle0$Reharvest.date!="2013-11-21",])
nsReharvestSingle<- droplevels(nsReharvestSingle1[nsReharvestSingle1$Reharvest.date!="2012-11-25",])
nsReharvestDouble<- droplevels(nsReharvestb[nsReharvestb$Harvest!="double" & nsReharvestb$Reharvest.date=="2013-11-21",])

# Restack 
# Total biomass is T0 = double = "2012-11-25" "2013-06-21" "2013-11-21"
# Regrowth H1 is T1 = single = (double) "2013-06-21" and (single) "2013-11-21"
# Regrowth H2 is T1 = double = "2013-11-21"
names(nsReharvestb2)
myvars<-c("X","Y","Harvest", "Harvest.date","Reharvest.date","Plot.name" ,"Trt.name","Treatment","Boma.density","Plot.pair",
         "GrassNetReharvestBiomass0","DwarfShrubNetReharvestBiomass0", "HerbNetReharvestBiomass0", "ClimberNetReharvestBiomass0",
          "TotalBiomass0", "min_distExclosures","boma_density")

myvars2<-c("X","Y","Harvest", "Harvest.date","Reharvest.date","Plot.name","Trt.name","Treatment","Boma.density","Plot.pair",
          "GrassNetReharvestBiomass1",  "DwarfShrubNetReharvestBiomass1", "HerbNetReharvestBiomass1", "ClimberNetReharvestBiomass1",
          "TotalBiomass1", "min_distExclosures","boma_density")

TotalBiomass<-nsReharvestb2[,myvars]
SingleBiomass<-nsReharvestSingle[,myvars2]
DoubleBiomass<-nsReharvestDouble[,myvars2]

names(TotalBiomass)
colnames(TotalBiomass)[11:15]<-c( "GrassNetReharvestBiomass", "DwarfShrubNetReharvestBiomass","HerbNetReharvestBiomass",
                               "ClimberNetReharvestBiomass", "TotalBiomass")
colnames(SingleBiomass)[11:15]<-c( "GrassNetReharvestBiomass", "DwarfShrubNetReharvestBiomass","HerbNetReharvestBiomass",
                                  "ClimberNetReharvestBiomass", "TotalBiomass")
colnames(DoubleBiomass)[11:15]<-c( "GrassNetReharvestBiomass", "DwarfShrubNetReharvestBiomass","HerbNetReharvestBiomass",
                                  "ClimberNetReharvestBiomass", "TotalBiomass")

dim(TotalBiomass)
dim(SingleBiomass)
dim(DoubleBiomass)

# Total biomass - change harvest name
levels(TotalBiomass$Harvest)<-"biomass"
levels(SingleBiomass$Harvest)<-c("regrowth-single","regrowth-single")
levels(DoubleBiomass$Harvest)<-"regrowth-double"

# Combined datasets
NechSarBiomass<-rbind(TotalBiomass,SingleBiomass,DoubleBiomass)

# Graph biomass
BiomassMeans<-aggregate(TotalBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassSD<-aggregate(TotalBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeans$sd<-BiomassSD$TotalBiomass

ggplot(BiomassMeans, aes(y=TotalBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
  geom_errorbar(aes(ymin=TotalBiomass-sd, ymax=TotalBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
  geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

#### Total -Mixed model with season as fixed factor ####
NechSarBiomass$fPlot.pair<-as.factor(NechSarBiomass$Plot.pair)
NechSarBiomass$fPlot.name<-as.factor(NechSarBiomass$Plot.name)
NechSarBiomass$fBoma.density<-as.factor(NechSarBiomass$Boma.density)
NechSarBiomass$Reharvest.date2<-as.factor(NechSarBiomass$Reharvest.date)
levels(NechSarBiomass$Reharvest.date2)<-c("Short I", "Long", "Short II")

NechSarBiomass$TotalBiomass[NechSarBiomass$TotalBiomass==0]<-0.1 # One zero!

table(NechSarBiomass$fPlot.name,NechSarBiomass$fPlot.pair) # 6 in each plot

# Total biomass model
BioRegrow<-glmmadmb(TotalBiomass~fBoma.density+Harvest+Treatment+
                fBoma.density:Harvest+Harvest:Treatment+
                fBoma.density:Treatment+
                fBoma.density:Harvest:Treatment+
                (1|fPlot.pair),family="Gamma", data=NechSarBiomass)
summary(BioRegrow)
anova(BioRegrow) 
AIC(BioRegrow) # 5435.36

# Residual vs fitted values
E0 <- resid(BioRegrow, type ="pearson")
F0 <- fitted(BioRegrow)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F0, 
     y = E0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Conical with gaussian - better gamma - non conical
#Improved with GLMM fitted gamma?

drop1(BioRegrow, test="Chisq")

# Generate pvalues
BioRegrow1<-update(BioRegrow,~.-fBoma.density:Harvest:Treatment)
BioRegrow2<-update(BioRegrow1,~.-Harvest:Treatment)
BioRegrow3<-update(BioRegrow1,~.-Treatment:fBoma.density)
BioRegrow4<-update(BioRegrow1,~.-fBoma.density:Harvest)
BioRegrow5<-update(BioRegrow2,~.-Treatment:fBoma.density)
BioRegrow6<-update(BioRegrow5,~.-fBoma.density:Harvest)

BioRegrow7<-update(BioRegrow6,~.-Harvest)
BioRegrow8<-update(BioRegrow6,~.-fBoma.density)
BioRegrow9<-update(BioRegrow6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrow,BioRegrow1) #
anova(BioRegrow1,BioRegrow2) # 
anova(BioRegrow1,BioRegrow3) #
anova(BioRegrow1,BioRegrow4) # 
anova(BioRegrow6,BioRegrow7) #
anova(BioRegrow6,BioRegrow8) #
anova(BioRegrow6,BioRegrow9) #

#  NoPar  LogLik Df Deviance Pr(>Chi)
#2    14 -2707.1  2     0.42   0.8106 #fBoma.density:Harvest:Treatment
#2    12 -2707.3  2      0.5   0.7788 #Harvest:Treatment
#2    12 -2707.3  1     1.14   0.2857 #Treatment:fBoma.density
#2    12 -2707.3  2     7.94  0.01887 * #fBoma.density:Harvest
#2     7 -2712.0  2   206.42 < 2.2e-16 *** #Harvest
#2     7 -2712.0  1     1.58   0.2088 #fBoma.density
#2     7 -2712.0  1     68.4 < 2.2e-16 ***#Treatment

aggregate(TotalBiomass~Treatment+Harvest,NechSarBiomass,mean)
#  Treatment         Harvest TotalBiomass
#1   Control         biomass     87.26889
#2 Exclosure         biomass    128.89259
#3   Control regrowth-single     37.03667
#4 Exclosure regrowth-single     53.48444
#5   Control regrowth-double     37.57556
#6 Exclosure regrowth-double     56.73333
(128.89259-87.26889)/87.26889 #0.4769592
(mean(53.48444,56.73333)-mean(37.03667,37.57556))/mean(37.03667,37.57556) # 0.4440942

aggregate(TotalBiomass~fBoma.density+Harvest,NechSarBiomass,mean)
aggregate(TotalBiomass~fBoma.density+Harvest,NechSarBiomass,sd)
mean(53.48444,56.73333)
mean(37.03667,37.57556)
mean(34.56563,38.66825)
mean(22.16488,25.24784)

(57.64200-52.79200)/52.79200*100
aggregate(TotalBiomass~Harvest,NechSarBiomass,mean)
57.64200-52.79200

# Interaction / tukuls density and season
with(NechSarBiomass, {interaction.plot(fBoma.density,Harvest,TotalBiomass,
                                       xlab = "Settlement",
                                       ylab = "Biomass",
                                       fun=mean)})

with(NechSarBiomass, {interaction.plot(Treatment,Harvest,TotalBiomass,
                                      xlab = "Treatment",
                                      ylab = "Biomass",
                                      fun=mean)})
# Driven by the slight increase in regrowth biomass 

#### Difference total biomass and regrowth ####
myvars3<-c("X","Y","Harvest", "Harvest.date","Reharvest.date","Trt.name","Plot.name","Treatment","Boma.density","Plot.pair",
           "min_distExclosures","boma_density","Regrow","RegrowGrass","RegrowWoody","RegrowForb")

# Calculate differences total biomass - regrowth
# Total
nsReharvestSingle$Regrow<-nsReharvestSingle$TotalBiomass1-nsReharvestSingle$TotalBiomass0
nsReharvestDouble$Regrow<-nsReharvestDouble$TotalBiomass1-nsReharvestDouble$TotalBiomass0
# Grass
nsReharvestSingle$RegrowGrass<-nsReharvestSingle$GrassNetReharvestBiomass1-nsReharvestSingle$GrassNetReharvestBiomass0
nsReharvestDouble$RegrowGrass<-nsReharvestDouble$GrassNetReharvestBiomass1-nsReharvestDouble$GrassNetReharvestBiomass0
#Woody
nsReharvestSingle$RegrowWoody<-nsReharvestSingle$DwarfShrubNetReharvestBiomass1-nsReharvestSingle$DwarfShrubNetReharvestBiomass0
nsReharvestDouble$RegrowWoody<-nsReharvestDouble$DwarfShrubNetReharvestBiomass1-nsReharvestDouble$DwarfShrubNetReharvestBiomass0
# Forbs
nsReharvestSingle$RegrowForb<-(nsReharvestSingle$ClimberNetReharvestBiomass1+nsReharvestSingle$HerbNetReharvestBiomass1)-(nsReharvestSingle$ClimberNetReharvestBiomass0+nsReharvestSingle$HerbNetReharvestBiomass0)
nsReharvestDouble$RegrowForb<-(nsReharvestDouble$ClimberNetReharvestBiomass1+nsReharvestDouble$HerbNetReharvestBiomass1)-(nsReharvestDouble$ClimberNetReharvestBiomass0+nsReharvestDouble$HerbNetReharvestBiomass0)

SingleBiomass2<-nsReharvestSingle[,myvars3]
DoubleBiomass2<-nsReharvestDouble[,myvars3]

levels(SingleBiomass2$Harvest)<-c("regrowth-single","regrowth-single")
levels(DoubleBiomass2$Harvest)<-"regrowth-double"

# Combined datasets
NechSarRegrow<-rbind(SingleBiomass2,DoubleBiomass2)

# Graph regrowth
RegrowMeans<-aggregate(Regrow~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowSD<-aggregate(Regrow~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeans$sd<-RegrowSD$Regrow
BiomassMean1<-BiomassMeans[BiomassMeans$Harvest!="biomass",]
RegrowMeans$TotalBiomass<-BiomassMean1$TotalBiomass

ggplot(RegrowMeans, aes(y=Regrow, x=Boma.density, col=Harvest,shape=Treatment, size=TotalBiomass))+
  geom_hline(yintercept =0, col="grey", linetype="dashed")+
  geom_errorbar(aes(ymin=Regrow-sd, ymax=Regrow+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
  geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
  theme_classic()

#### Total regrowth -Mixed model with season as fixed factor ####
NechSarRegrow$fPlot.pair<-as.factor(NechSarRegrow$Plot.pair)
NechSarRegrow$fPlot.name<-as.factor(NechSarRegrow$Plot.name)
NechSarRegrow$fBoma.density<-as.factor(NechSarRegrow$Boma.density)
NechSarRegrow$Reharvest.date2<-as.factor(NechSarRegrow$Reharvest.date)
levels(NechSarRegrow$Reharvest.date2)<-c("Long", "Short II")

BioRegrowB<-lmer(Regrow~fBoma.density+Treatment+Harvest+
                   Harvest:Treatment+ fBoma.density:Harvest+  
                  fBoma.density:Treatment+
                   fBoma.density:Treatment:Harvest+
                  (1|fPlot.pair),
                  data=NechSarRegrow)
summary(BioRegrowB)
anova(BioRegrowB) # Highly signficant
AIC(BioRegrowB) # 2826.343

# Residual vs fitted values
E02<- resid(BioRegrowB, type ="pearson")
F02 <- fitted(BioRegrowB)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F02, 
     y = E02,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Good

drop1(BioRegrowB, test="Chisq")

# Generate pvalues
BioRegrowB1<-update(BioRegrowB,~.-fBoma.density:Harvest:Treatment)
BioRegrowB2<-update(BioRegrowB1,~.-Harvest:Treatment)
BioRegrowB3<-update(BioRegrowB1,~.-Treatment:fBoma.density)
BioRegrowB4<-update(BioRegrowB1,~.-fBoma.density:Harvest)
BioRegrowB5<-update(BioRegrowB2,~.-Treatment:fBoma.density)
BioRegrowB6<-update(BioRegrowB5,~.-fBoma.density:Harvest)

BioRegrowB7<-update(BioRegrowB6,~.-Harvest)
BioRegrowB8<-update(BioRegrowB6,~.-fBoma.density)
BioRegrowB9<-update(BioRegrowB6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowB,BioRegrowB1) #
anova(BioRegrowB1,BioRegrowB2) # 
anova(BioRegrowB1,BioRegrowB3) #
anova(BioRegrowB1,BioRegrowB4) # 
anova(BioRegrowB6,BioRegrowB7) #
anova(BioRegrowB6,BioRegrowB8) #
anova(BioRegrowB6,BioRegrowB9) #

#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#BioRegrowB  10 2878.6 2914.6 -1429.3   2858.6 3.0025      1    0.08314 . #fBoma.density:Harvest:Treatment
#BioRegrowB  10 2878.6 2914.6 -1429.3   2858.6 3.9806      2     0.1367 #Harvest:Treatment
#BioRegrowB1  9 2879.6 2912.0 -1430.8   2861.6 2.3153      1     0.1281 #Treatment:fBoma.density
#BioRegrowB1  9 2879.6 2912.0 -1430.8   2861.6 0.9624      1     0.3266 #fBoma.density:Harvest
#BioRegrowB6  6 2877.8 2899.4 -1432.9   2865.8 0.379      1     0.5381 #Harvest
#BioRegrowB6  6 2877.8 2899.4 -1432.9   2865.8 0.206      1     0.6499 #fBoma.density
#BioRegrowB6  6 2877.8 2899.4 -1432.9   2865.8 12.264      1  0.0004616 *** #Treatment

#### Grass -Mixed model with season as fixed factor ####
names(NechSarBiomass)
#"GrassNetReharvestBiomass""DwarfShrubNetReharvestBiomass" "HerbNetReharvestBiomass" "ClimberNetReharvestBiomass" 

#### TEST DISTANCE TO SETTLEMENTS - NOT DENSITY #####
names(nsReharvestH1) #min_distExclosures
plot(nsReharvestH1$min_distExclosures)

#### Total -Mixed model with season as fixed factor - NEAREST Settlements ####
aggregate(TotalBiomass~min_distExclosures+Treatment, NechSarBiomass,mean)
aggregate(TotalBiomass~min_distExclosures+Treatment, NechSarBiomass,sd)

Tot0modDIS<-glmmadmb(TotalBiomass~min_distExclosures+Harvest+Treatment+
                     Treatment:min_distExclosures+Treatment:Harvest+
                     Harvest:min_distExclosures+
                     Harvest:min_distExclosures:Treatment+
                   (1|fPlot.pair),family="Gamma", data=NechSarBiomass)
summary(Tot0modDIS)
AIC(Tot0modDIS) #5422.1

# Residual vs fitted values
E0 <- resid(Tot0modDIS, type ="pearson")
F0 <- fitted(Tot0modDIS)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F0, 
     y = E0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK

#drop1(Tot0modDIS, test="Chisq")

# Generate pvalues
Tot0modDIS1<-update(Tot0modDIS,~.-Treatment:min_distExclosures:Harvest)
Tot0modDIS2<-update(Tot0modDIS1,~.-Treatment:Harvest)
Tot0modDIS3<-update(Tot0modDIS1,~.-Treatment:min_distExclosures)
Tot0modDIS4<-update(Tot0modDIS1,~.-Harvest:min_distExclosures)
Tot0modDIS5<-update(Tot0modDIS2,~.-Treatment:min_distExclosures)
Tot0modDIS6<-update(Tot0modDIS5,~.-Harvest:min_distExclosures)

Tot0modDIS7<-update(Tot0modDIS6,~.-Harvest)
Tot0modDIS8<-update(Tot0modDIS6,~.-min_distExclosures)
Tot0modDIS9<-update(Tot0modDIS6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Tot0modDIS,Tot0modDIS1) #
anova(Tot0modDIS1,Tot0modDIS2) # 
anova(Tot0modDIS1,Tot0modDIS3) #
anova(Tot0modDIS1,Tot0modDIS4) # 
anova(Tot0modDIS6,Tot0modDIS7) #
anova(Tot0modDIS6,Tot0modDIS8) #
anova(Tot0modDIS6,Tot0modDIS9) #

#NoPar  LogLik Df Deviance Pr(>Chi)
#2    14 -2697.1  2     1.54    0.463 # Treatment:min_distExclosures:Harvest
#2    12 -2697.8  2     0.64   0.7261 #Treatment:Harvest
#2    12 -2697.8  1    20.82 5.045e-06 *** #Treatment:min_distExclosures
#2    12 -2697.8  2      1.3    0.522 #Harvest:min_distExclosures
#2     7 -2709.2  2    206.3 < 2.2e-16 *** #Harvest
#2     7 -2709.2  1     7.14 0.007538 ** #min_distExclosures
#2     7 -2709.2  1    69.02 < 2.2e-16 *** #Treatment

# Interaction / settlement density and season
xyplot(TotalBiomass~min_distExclosures|Treatment*Harvest, NechSarBiomass)
names(NechSarRegrow)
# Distance to settlement regrowth difference
BioRegrowDist<-lmer(Regrow~min_distExclosures+Harvest+Treatment+
                      Treatment:min_distExclosures+Treatment:Harvest+
                      Harvest:min_distExclosures+
                      Harvest:min_distExclosures:Treatment+
                   (1|fPlot.pair),
                 data=NechSarRegrow)
summary(BioRegrowDist)
anova(BioRegrowDist) # Highly signficant
AIC(BioRegrowDist) # 2890.265

# Residual vs fitted values
E02d<- resid(BioRegrowDist, type ="pearson")
F02d <- fitted(BioRegrowDist)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F02d, 
     y = E02d,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Good

#drop1(BioRegrowDist, test="Chisq")

# Generate pvalues
BioRegrowDist1<-update(BioRegrowDist,~.-min_distExclosures:Harvest:Treatment)
BioRegrowDist2<-update(BioRegrowDist1,~.-Harvest:Treatment)
BioRegrowDist3<-update(BioRegrowDist1,~.-Treatment:min_distExclosures)
BioRegrowDist4<-update(BioRegrowDist1,~.-min_distExclosures:Harvest)
BioRegrowDist5<-update(BioRegrowDist2,~.-Treatment:min_distExclosures)
BioRegrowDist6<-update(BioRegrowDist5,~.-min_distExclosures:Harvest)

BioRegrowDist7<-update(BioRegrowDist6,~.-Harvest)
BioRegrowDist8<-update(BioRegrowDist6,~.-min_distExclosures)
BioRegrowDist9<-update(BioRegrowDist6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowDist,BioRegrowDist1) #
anova(BioRegrowDist1,BioRegrowDist2) # 
anova(BioRegrowDist1,BioRegrowDist3) #
anova(BioRegrowDist1,BioRegrowDist4) # 
anova(BioRegrowDist6,BioRegrowDist7) #
anova(BioRegrowDist6,BioRegrowDist8) #
anova(BioRegrowDist6,BioRegrowDist9) #

#NoPar  LogLik Df Deviance Pr(>Chi)
#BioRegrowDist  10 2880.6 2916.6 -1430.3   2860.6 0.8472      1     0.3574 # Treatment:min_distExclosures:Harvest
#BioRegrowDist1  9 2879.4 2911.8 -1430.7   2861.4 0.969      1     0.3249  #Treatment:Harvest
#BioRegrowDist1  9 2879.4 2911.8 -1430.7   2861.4 0.5492      1     0.4587 #Treatment:min_distExclosures
#BioRegrowDist1  9 2879.4 2911.8 -1430.7   2861.4 0.2859      1     0.5928 #Harvest:min_distExclosures
#BioRegrowDist6  6 2875.2 2896.8 -1431.6   2863.2 0.379      1     0.5381 #Harvest
#BioRegrowDist6  6 2875.2 2896.8 -1431.6   2863.2 2.7848      1    0.09516 . #min_distExclosures
#BioRegrowDist6  6 2875.2 2896.8 -1431.6   2863.2 12.264      1  0.0004616 *** #Treatment


##### Graph grass biomass ####
BiomassMeansG<-aggregate(GrassNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassGSD<-aggregate(GrassNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeansG$sd<-BiomassGSD$GrassNetReharvestBiomass

ggplot(BiomassMeansG, aes(y=GrassNetReharvestBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=GrassNetReharvestBiomass-sd, ymax=GrassNetReharvestBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

NechSarBiomass$GrassNetReharvestBiomass[NechSarBiomass$GrassNetReharvestBiomass==0]<-0.1
#NechSarBiomassG<-NechSarBiomass[!NechSarBiomass$GrassNetReharvestBiomass==0,]
NechSarBiomassG<-NechSarBiomassG[!is.na(NechSarBiomassG$GrassNetReharvestBiomass),]

BioRegrowG<-glmmadmb(GrassNetReharvestBiomass~fBoma.density+Harvest+Treatment+
                  fBoma.density:Harvest+Harvest:Treatment+
                  fBoma.density:Treatment+
                  fBoma.density:Harvest:Treatment+
                  (1|fPlot.pair),family="Gamma", data=NechSarBiomassG)
summary(BioRegrowG)
#anova(BioRegrowG) 
AIC(BioRegrowG) # 5240.86

# Residual vs fitted values
EG0 <- resid(BioRegrowG, type ="pearson")
FG0 <- fitted(BioRegrowG)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = FG0, 
     y = EG0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Conical...

#drop1(BioRegrowG, test="Chisq")

# Generate pvalues
BioRegrowG1<-update(BioRegrowG,~.-fBoma.density:Harvest:Treatment)
BioRegrowG2<-update(BioRegrowG1,~.-Harvest:Treatment)
BioRegrowG3<-update(BioRegrowG1,~.-Treatment:fBoma.density)
BioRegrowG4<-update(BioRegrowG1,~.-fBoma.density:Harvest)
BioRegrowG5<-update(BioRegrowG2,~.-Treatment:fBoma.density)
BioRegrowG6<-update(BioRegrowG5,~.-fBoma.density:Harvest)

BioRegrowG7<-update(BioRegrowG6,~.-Harvest)
BioRegrowG8<-update(BioRegrowG6,~.-fBoma.density)
BioRegrowG9<-update(BioRegrowG6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowG,BioRegrowG1) #
anova(BioRegrowG1,BioRegrowG2) # 
anova(BioRegrowG1,BioRegrowG3) #
anova(BioRegrowG1,BioRegrowG4) # 
anova(BioRegrowG6,BioRegrowG7) #
anova(BioRegrowG6,BioRegrowG8) #
anova(BioRegrowG6,BioRegrowG9) #

#  NoPar  LogLik Df Deviance Pr(>Chi)
#  2    14 -2606.4  2     0.24   0.8869 #fBoma.density:Harvest:Treatment
#  2    12 -2606.6  2     0.98   0.6126 #Harvest:Treatment
#  2    12 -2606.6  1     0.22    0.639 #Treatment:fBoma.density
#  2    12 -2606.6  2     3.58    0.167 #fBoma.density:Harvest
#  2     7 -2608.9  2   219.26 < 2.2e-16 *** #Harvest
#  2     7 -2608.9  1     0.38   0.5376 #fBoma.density
#  2     7 -2608.9  1    87.18 < 2.2e-16 *** #Treatment

#### Grass regrowth model ####
# Graph Grass regrowth
RegrowMeansG<-aggregate(RegrowGrass~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowGSD<-aggregate(RegrowGrass~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeansG$sd<-RegrowGSD$RegrowGrass
BiomassMeanG1<-BiomassMeansG[BiomassMeans$Harvest!="biomass",]
RegrowMeansG$GrassNetReharvestBiomass<-BiomassMeanG1$GrassNetReharvestBiomass

ggplot(RegrowMeansG, aes(y=RegrowGrass, x=Boma.density, col=Harvest,shape=Treatment, size=GrassNetReharvestBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=RegrowGrass-sd, ymax=RegrowGrass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
theme_classic()

#### Grass regrowth -Mixed model with season as fixed factor ####
BioRegrowBG<-lmer(RegrowGrass~fBoma.density+Treatment+Harvest+
                   Harvest:Treatment+ fBoma.density:Harvest+  
                   fBoma.density:Treatment+
                   fBoma.density:Treatment:Harvest+
                   (1|fPlot.pair), data=NechSarRegrow)
summary(BioRegrowBG)
anova(BioRegrowBG) # Highly signficant
AIC(BioRegrowBG) # 2806.08

# Residual vs fitted values
E02G<- resid(BioRegrowBG, type ="pearson")
F02G <- fitted(BioRegrowBG)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F02G, 
     y = E02G,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Good

drop1(BioRegrowBG, test="Chisq")

# Generate pvalues
BioRegrowBG1<-update(BioRegrowBG,~.-fBoma.density:Harvest:Treatment)
BioRegrowBG2<-update(BioRegrowBG1,~.-Harvest:Treatment)
BioRegrowBG3<-update(BioRegrowBG1,~.-Treatment:fBoma.density)
BioRegrowBG4<-update(BioRegrowBG1,~.-fBoma.density:Harvest)
BioRegrowBG5<-update(BioRegrowBG2,~.-Treatment:fBoma.density)
BioRegrowBG6<-update(BioRegrowBG5,~.-fBoma.density:Harvest)

BioRegrowBG7<-update(BioRegrowBG6,~.-Harvest)
BioRegrowBG8<-update(BioRegrowBG6,~.-fBoma.density)
BioRegrowBG9<-update(BioRegrowBG6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowBG,BioRegrowBG1) #
anova(BioRegrowBG1,BioRegrowBG2) # 
anova(BioRegrowBG1,BioRegrowBG3) #
anova(BioRegrowBG1,BioRegrowBG4) # 
anova(BioRegrowBG6,BioRegrowBG7) #
anova(BioRegrowBG6,BioRegrowBG8) #
anova(BioRegrowBG6,BioRegrowBG9) #

#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#BioRegrowBG  10 2857.4 2893.4 -1418.7   2837.4 2.0435      1     0.1529 #fBoma.density:Harvest:Treatment
#BioRegrowBG1  9 2857.4 2889.8 -1419.7   2839.4 1.3234      1       0.25 #Harvest:Treatment
#BioRegrowBG1  9 2857.4 2889.8 -1419.7   2839.4 2.4034      1     0.1211 #Treatment:fBoma.density
#BioRegrowBG1  9 2857.4 2889.8 -1419.7   2839.4 0.7802      1     0.3771 #fBoma.density:Harvest
#BioRegrowBG6  6 2855.9 2877.5 -1422.0   2843.9 0.4814      1     0.4878 #Harvest
#BioRegrowBG6  6 2855.9 2877.5  -1422   2843.9 0.102      1     0.7495 #fBoma.density
#BioRegrowBG6  6 2855.9 2877.5 -1422.0   2843.9 12.649      1  0.0003758 *** #Treatment

#### Forb -Mixed model with season as fixed factor ####
names(NechSarBiomass)
#"GrassNetReharvestBiomass""DwarfShrubNetReharvestBiomass" "HerbNetReharvestBiomass" "ClimberNetReharvestBiomass" 
NechSarBiomass$ForbNetReharvestBiomass<-NechSarBiomass$HerbNetReharvestBiomass+NechSarBiomass$ClimberNetReharvestBiomass

# Graph forb biomass
BiomassMeansF<-aggregate(ForbNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassFSD<-aggregate(ForbNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeansF$sd<-BiomassFSD$ForbNetReharvestBiomass

ggplot(BiomassMeansF, aes(y=ForbNetReharvestBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=ForbNetReharvestBiomass-sd, ymax=ForbNetReharvestBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

length(NechSarBiomass$ForbNetReharvestBiomass[NechSarBiomass$ForbNetReharvestBiomass==0])
length(NechSarBiomass$ForbNetReharvestBiomass)
283/540 # 53% = zero
NechSarBiomass$ForbNetReharvestBiomass[NechSarBiomass$ForbNetReharvestBiomass==0]<-0.01

# Forb biomass model 
BioRegrowF<-glmmadmb(ForbNetReharvestBiomass~fBoma.density+Harvest+Treatment+
                   fBoma.density:Harvest+Harvest:Treatment+
                   fBoma.density:Treatment+
                   fBoma.density:Harvest:Treatment+
                   (1|fPlot.pair), family="Gamma",#zeroInflation=TRUE,
                   data=NechSarBiomass) #zeroInflation=TRUE
summary(BioRegrowF)
anova(BioRegrowF) 
AIC(BioRegrowF) # 1287.426

# Residual vs fitted values
EF0 <- resid(BioRegrowF, type ="pearson")
FF0 <- fitted(BioRegrowF)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = FF0, 
     y = EF0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # strong edge zeros - negative zero fit with guassian
# Need zero inflation term = non negative predictions

#drop1(BioRegrowF, test="Chisq")

# Generate pvalues
BioRegrowF1<-update(BioRegrowF,~.-fBoma.density:Harvest:Treatment)
BioRegrowF2<-update(BioRegrowF1,~.-Harvest:Treatment)
BioRegrowF3<-update(BioRegrowF1,~.-Treatment:fBoma.density)
BioRegrowF4<-update(BioRegrowF1,~.-fBoma.density:Harvest)
BioRegrowF5<-update(BioRegrowF2,~.-Treatment:fBoma.density)
BioRegrowF6<-update(BioRegrowF5,~.-fBoma.density:Harvest)

BioRegrowF7<-update(BioRegrowF6,~.-Harvest)
BioRegrowF8<-update(BioRegrowF6,~.-fBoma.density)
BioRegrowF9<-update(BioRegrowF6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowF,BioRegrowF1) #
anova(BioRegrowF1,BioRegrowF2) # 
anova(BioRegrowF1,BioRegrowF3) #
anova(BioRegrowF1,BioRegrowF4) # 
anova(BioRegrowF6,BioRegrowF7) #
anova(BioRegrowF6,BioRegrowF8) #
anova(BioRegrowF6,BioRegrowF9) #

#  NoPar  LogLik Df Deviance  Pr(>Chi)   
# 2    14 -629.71  2    1.994    0.369 #fBoma.density:Harvest:Treatment
# 2    12 -630.71  2    1.766   0.4135 #Harvest:Treatment
# 2    12 -630.71  1    30.87 2.759e-08 *** #Treatment:fBoma.density
# 2    12 -630.71  2    3.136   0.2085 #fBoma.density:Harvest
# 2     7 -646.76  2    4.268   0.1184 #Harvest
# 2     7 -646.76  1    6.116   0.0134 *#fBoma.density
# 2     7 -646.76  1     5.38  0.02037 *#Treatment


# Graph Forb regrowth
names(NechSarRegrow)
RegrowMeansF<-aggregate(RegrowForb~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowFSD<-aggregate(RegrowForb~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeansF$sd<-RegrowFSD$RegrowForb
BiomassMeanF1<-BiomassMeansF[BiomassMeans$Harvest!="biomass",]
RegrowMeansF$ForbNetReharvestBiomass<-BiomassMeanF1$ForbNetReharvestBiomass

ggplot(RegrowMeansF, aes(y=RegrowForb, x=Boma.density, col=Harvest,shape=Treatment, size=ForbNetReharvestBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=RegrowForb-sd, ymax=RegrowForb+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
theme_classic()

#### Forb regrowth -Mixed model with season as fixed factor ####
BioRegrowBF<-lmer(RegrowForb~fBoma.density+Treatment+Harvest+
                    Harvest:Treatment+ fBoma.density:Harvest+  
                    fBoma.density:Treatment+
                    fBoma.density:Treatment:Harvest+
                    (1|fPlot.pair), data=NechSarRegrow)
summary(BioRegrowBF)
anova(BioRegrowBF) # Highly signficant
AIC(BioRegrowBF) # 2224.878

# Residual vs fitted values
E02F<- resid(BioRegrowBF, type ="pearson")
F02F <- fitted(BioRegrowBF)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F02F, 
     y = E02F,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Good

drop1(BioRegrowBF, test="Chisq")

# Generate pvalues
BioRegrowBF1<-update(BioRegrowBF,~.-fBoma.density:Harvest:Treatment)
BioRegrowBF2<-update(BioRegrowBF1,~.-Harvest:Treatment)
BioRegrowBF3<-update(BioRegrowBF1,~.-Treatment:fBoma.density)
BioRegrowBF4<-update(BioRegrowBF1,~.-fBoma.density:Harvest)
BioRegrowBF5<-update(BioRegrowBF2,~.-Treatment:fBoma.density)
BioRegrowBF6<-update(BioRegrowBF5,~.-fBoma.density:Harvest)

BioRegrowBF7<-update(BioRegrowBF6,~.-Harvest)
BioRegrowBF8<-update(BioRegrowBF6,~.-fBoma.density)
BioRegrowBF9<-update(BioRegrowBF6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowBF,BioRegrowBF1) #
anova(BioRegrowBF1,BioRegrowBF2) # 
anova(BioRegrowBF1,BioRegrowBF3) #
anova(BioRegrowBF1,BioRegrowBF4) # 
anova(BioRegrowBF6,BioRegrowBF7) #
anova(BioRegrowBF6,BioRegrowBF8) #
anova(BioRegrowBF6,BioRegrowBF9) #

#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#BioRegrowBF  10 2256.5 2292.5 -1118.3   2236.5 0.0601      1     0.8064 #fBoma.density:Harvest:Treatment
#BioRegrowBF1  9 2254.6 2287.0 -1118.3   2236.6 0.1327      1     0.7157 #Harvest:Treatment
#BioRegrowBF1  9 2254.6 2287.0 -1118.3   2236.6 1.2895      1     0.2561 #Treatment:fBoma.density
#BioRegrowBF1  9 2254.6 2287.0 -1118.3   2236.6 0.0883      1     0.7664 # fBoma.density:Harvest
#BioRegrowBF6  6 2250.1 2271.7 -1119.0   2238.1 0.1214      1     0.7275 #Harvest
#BioRegrowBF6  6 2250.1 2271.7 -1119.0   2238.1 0.1363      1      0.712 #fBoma.density
#BioRegrowBF6  6 2250.1 2271.7 -1119.0   2238.1 0.668       1     0.4137 #Treatment


#### Dwarf -Mixed model with season as fixed factor ####
names(NechSarBiomass)
#"GrassNetReharvestBiomass""DwarfShrubNetReharvestBiomass" "HerbNetReharvestBiomass" "ClimberNetReharvestBiomass" 

#### Woody biomass ####
BiomassMeansW<-aggregate(DwarfShrubNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassWSD<-aggregate(DwarfShrubNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeansW$sd<-BiomassWSD$DwarfShrubNetReharvestBiomass

ggplot(BiomassMeansW, aes(y=DwarfShrubNetReharvestBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=DwarfShrubNetReharvestBiomass-sd, ymax=DwarfShrubNetReharvestBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

# Dwarfshrubs
length(NechSarBiomass$DwarfShrubNetReharvestBiomass[NechSarBiomass$DwarfShrubNetReharvestBiomass==0])
length(NechSarBiomass$DwarfShrubNetReharvestBiomass)
342/540 # 63% = zero
NechSarBiomass$DwarfShrubNetReharvestBiomass[NechSarBiomass$DwarfShrubNetReharvestBiomass==0]<-0.01
NechSarBiomass$DwarfShrubNetReharvestBiomass[NechSarBiomass$DwarfShrubNetReharvestBiomass<0]<-1.5  # One negative value -1.5??

BioRegrowW<-glmmadmb(DwarfShrubNetReharvestBiomass~fBoma.density+Harvest+Treatment+
                   fBoma.density:Harvest+Harvest:Treatment+
                   fBoma.density:Treatment+
                   fBoma.density:Harvest:Treatment+
                   (1|fPlot.pair), family="gamma",
#admb.opts=admbControl(shess=FALSE,noinit=FALSE),
#admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                   data=NechSarBiomass)

summary(BioRegrowW)
AIC(BioRegrowW) # 679.61

# Residual vs fitted values
EW0 <- resid(BioRegrowW, type ="pearson")
FW0 <- fitted(BioRegrowW)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = FW0, 
     y = EW0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # strong edge zeros - negative zero fit with guassian
# Need zero inflation term = non negative predictions

#drop1(BioRegrowW, test="Chisq")

# Generate pvalues
BioRegrowW1<-update(BioRegrowW,~.-fBoma.density:Harvest:Treatment)
BioRegrowW2<-update(BioRegrowW1,~.-Harvest:Treatment)
BioRegrowW3<-update(BioRegrowW1,~.-Treatment:fBoma.density)
BioRegrowW4<-update(BioRegrowW1,~.-fBoma.density:Harvest)
BioRegrowW5<-update(BioRegrowW2,~.-Treatment:fBoma.density)
BioRegrowW6<-update(BioRegrowW5,~.-fBoma.density:Harvest)

BioRegrowW7<-update(BioRegrowW6,~.-Harvest)
BioRegrowW8<-update(BioRegrowW6,~.-fBoma.density)
BioRegrowW9<-update(BioRegrowW6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowW,BioRegrowW1) #
anova(BioRegrowW1,BioRegrowW2) # 
anova(BioRegrowW1,BioRegrowW3) #
anova(BioRegrowW1,BioRegrowW4) # 
anova(BioRegrowW6,BioRegrowW7) #
anova(BioRegrowW6,BioRegrowW8) #
anova(BioRegrowW6,BioRegrowW9) #

#  NoPar  LogLik Df Deviance  Pr(>Chi)   
#  2    14 -325.81  2    0.864   0.6492 #fBoma.density:Harvest:Treatment
#  2    12 -326.24  2    0.818   0.6643 #Harvest:Treatment
#  2    12 -326.24  1    0.114   0.7356 #Treatment:fBoma.density
#  2    12 -326.24  2    2.696   0.2598 #fBoma.density:Harvest
#  2     7 -328.11  2    7.796  0.02028 * #Harvest
#  2     7 -328.11  1    0.994   0.3188 #fBoma.density
#  2     7 -328.11  1    0.434     0.51 #Treatment

aggregate(DwarfShrubNetReharvestBiomass~Harvest,NechSarBiomass, mean)


### Woody regrowth ####
RegrowMeansW<-aggregate(RegrowWoody~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowWSD<-aggregate(RegrowWoody~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeansW$sd<-RegrowWSD$RegrowWoody
BiomassMeanW1<-BiomassMeansW[BiomassMeans$Harvest!="biomass",]
RegrowMeansW$DwarfShrubNetReharvestBiomass<-BiomassMeanW1$DwarfShrubNetReharvestBiomass

ggplot(RegrowMeansW, aes(y=RegrowWoody, x=Boma.density, col=Harvest,shape=Treatment, size=DwarfShrubNetReharvestBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=RegrowWoody-sd, ymax=RegrowWoody+sd),width=.1,lwd=1,position=position_dodge(width=.5),show.legend=F)+
geom_point(position=position_dodge(width=.5))+scale_radius(range=c(2,8))+
theme_classic()

plot(NechSarRegrow$RegrowWoody)

#### Woody regrowth -Mixed model with season as fixed factor ####
BioRegrowBW<-lmer(RegrowWoody~fBoma.density+Treatment+Harvest+
                    Harvest:Treatment+ fBoma.density:Harvest+  
                    fBoma.density:Treatment+
                    fBoma.density:Treatment:Harvest+
                    (1|fPlot.pair), data=NechSarRegrow)
summary(BioRegrowBW)
anova(BioRegrowBW) # Highly signficant
AIC(BioRegrowBW) # 2148.093

# Residual vs fitted values
E02W<- resid(BioRegrowBW, type ="pearson")
F02W <- fitted(BioRegrowBW)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F02W, 
     y = E02W,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Clusters?

drop1(BioRegrowBW, test="Chisq")

# Generate pvalues
BioRegrowBW1<-update(BioRegrowBW,~.-fBoma.density:Harvest:Treatment)
BioRegrowBW2<-update(BioRegrowBW1,~.-Harvest:Treatment)
BioRegrowBW3<-update(BioRegrowBW1,~.-Treatment:fBoma.density)
BioRegrowBW4<-update(BioRegrowBW1,~.-fBoma.density:Harvest)
BioRegrowBW5<-update(BioRegrowBW2,~.-Treatment:fBoma.density)
BioRegrowBW6<-update(BioRegrowBW5,~.-fBoma.density:Harvest)

BioRegrowBW7<-update(BioRegrowBW6,~.-Harvest)
BioRegrowBW8<-update(BioRegrowBW6,~.-fBoma.density)
BioRegrowBW9<-update(BioRegrowBW6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(BioRegrowBW,BioRegrowBW1) #
anova(BioRegrowBW1,BioRegrowBW2) # 
anova(BioRegrowBW1,BioRegrowBW3) #
anova(BioRegrowBW1,BioRegrowBW4) # 
anova(BioRegrowBW6,BioRegrowBW7) #
anova(BioRegrowBW6,BioRegrowBW8) #
anova(BioRegrowBW6,BioRegrowBW9) #

#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#BioRegrowBW  10 2176.5 2212.5 -1078.2   2156.5 0.9808      1      0.322 #fBoma.density:Harvest:Treatment
#BioRegrowBW1  9 2175.5 2207.9 -1078.7   2157.5 0.712      1     0.3988  #Harvest:Treatment
#BioRegrowBW1  9 2175.5 2207.9 -1078.7   2157.5 1.4405      1     0.2301  #Treatment:fBoma.density
#BioRegrowBW1  9 2175.5 2207.9 -1078.7   2157.5 0.0165      1     0.8978  #fBoma.density:Harvest
#BioRegrowBW6  6 2171.6 2193.2 -1079.8   2159.6 0.0446      1     0.8327  #Harvest
#BioRegrowBW6  6 2171.6 2193.2 -1079.8   2159.6 3.4215      1    0.06435 .  #fBoma.density
#BioRegrowBW6  6 2171.6 2193.2 -1079.8   2159.6 1.4769      1     0.2243  #Treatment

######################################################################
#### Biomass and regrowth - plotting data ####
########################################################################

##### Grassland biomass #####
BiomassMeans<-aggregate(TotalBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassSD<-aggregate(TotalBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeans$sd<-BiomassSD$TotalBiomass

ggplot(BiomassMeans, aes(y=TotalBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=TotalBiomass-sd, ymax=TotalBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()


##### Grassland regrowth #####
RegrowMeans<-aggregate(Regrow~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowSD<-aggregate(Regrow~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeans$sd<-RegrowSD$Regrow
BiomassMean1<-BiomassMeans[BiomassMeans$Harvest!="biomass",]
RegrowMeans$TotalBiomass<-BiomassMean1$TotalBiomass

ggplot(RegrowMeans, aes(y=Regrow, x=Boma.density, col=Harvest,shape=Treatment, size=TotalBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=Regrow-sd, ymax=Regrow+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
theme_classic()

##### Grass biomass ####
BiomassMeansG<-aggregate(GrassNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassGSD<-aggregate(GrassNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeansG$sd<-BiomassGSD$GrassNetReharvestBiomass

ggplot(BiomassMeansG, aes(y=GrassNetReharvestBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=GrassNetReharvestBiomass-sd, ymax=GrassNetReharvestBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

##### Grass regrowth ####
RegrowMeansG<-aggregate(RegrowGrass~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowGSD<-aggregate(RegrowGrass~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeansG$sd<-RegrowGSD$RegrowGrass
BiomassMeanG1<-BiomassMeansG[BiomassMeans$Harvest!="biomass",]
RegrowMeansG$GrassNetReharvestBiomass<-BiomassMeanG1$GrassNetReharvestBiomass

ggplot(RegrowMeansG, aes(y=RegrowGrass, x=Boma.density, col=Harvest,shape=Treatment, size=GrassNetReharvestBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=RegrowGrass-sd, ymax=RegrowGrass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
theme_classic()

##### Forb biomass #####
BiomassMeansF<-aggregate(ForbNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassFSD<-aggregate(ForbNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeansF$sd<-BiomassFSD$ForbNetReharvestBiomass

ggplot(BiomassMeansF, aes(y=ForbNetReharvestBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=ForbNetReharvestBiomass-sd, ymax=ForbNetReharvestBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

#### Forb regrowth ####
RegrowMeansF<-aggregate(RegrowForb~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowFSD<-aggregate(RegrowForb~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeansF$sd<-RegrowFSD$RegrowForb
BiomassMeanF1<-BiomassMeansF[BiomassMeans$Harvest!="biomass",]
RegrowMeansF$ForbNetReharvestBiomass<-BiomassMeanF1$ForbNetReharvestBiomass

ggplot(RegrowMeansF, aes(y=RegrowForb, x=Boma.density, col=Harvest,shape=Treatment, size=ForbNetReharvestBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=RegrowForb-sd, ymax=RegrowForb+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
theme_classic()

#### Woody biomass ####
BiomassMeansW<-aggregate(DwarfShrubNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,mean)
BiomassWSD<-aggregate(DwarfShrubNetReharvestBiomass~Treatment+Boma.density+Harvest,NechSarBiomass,sd)
BiomassMeansW$sd<-BiomassWSD$DwarfShrubNetReharvestBiomass

ggplot(BiomassMeansW, aes(y=DwarfShrubNetReharvestBiomass, x=Boma.density, colour=Harvest, shape=Treatment))+
geom_errorbar(aes(ymin=DwarfShrubNetReharvestBiomass-sd, ymax=DwarfShrubNetReharvestBiomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

#### Woody regrowth ####
RegrowMeansW<-aggregate(RegrowWoody~Treatment+Boma.density+Harvest,NechSarRegrow,mean)
RegrowWSD<-aggregate(RegrowWoody~Treatment+Boma.density+Harvest,NechSarRegrow,sd)
RegrowMeansW$sd<-RegrowWSD$RegrowWoody
BiomassMeanW1<-BiomassMeansW[BiomassMeans$Harvest!="biomass",]
RegrowMeansW$DwarfShrubNetReharvestBiomass<-BiomassMeanW1$DwarfShrubNetReharvestBiomass

ggplot(RegrowMeansW, aes(y=RegrowWoody, x=Boma.density, col=Harvest,shape=Treatment, size=DwarfShrubNetReharvestBiomass))+
geom_hline(yintercept =0, col="grey", linetype="dashed")+
geom_errorbar(aes(ymin=RegrowWoody-sd, ymax=RegrowWoody+sd),width=.1,lwd=1,position=position_dodge(width=.5),show.legend=F)+
geom_point(position=position_dodge(width=.5))+scale_radius(range=c(2,8))+
theme_classic()


#### Combine datasets ####

colnames(BiomassMeans)[4]<-"Biomass"
colnames(BiomassMeansG)[4]<-"Biomass"
colnames(BiomassMeansF)[4]<-"Biomass"
colnames(BiomassMeansW)[4]<-"Biomass"
BiomassMeans$FxGroup<-"Total"
BiomassMeansG$FxGroup<-"Grasses"
BiomassMeansF$FxGroup<-"Forbs"
BiomassMeansW$FxGroup<-"Woody"

BiomassMeansALL<-rbind(BiomassMeans,BiomassMeansG,BiomassMeansF,BiomassMeansW)
BiomassMeansALL$FxGroup<-as.factor(BiomassMeansALL$FxGroup)
levels(BiomassMeansALL$FxGroup)<-c("Forbs","Grasses","Total","Woody")
BiomassMeansALL$FxGroup<- factor(BiomassMeansALL$FxGroup, levels = c("Total", "Grasses","Forbs","Woody"))

BioGraph<-ggplot(BiomassMeansALL, aes(y=Biomass, x=Boma.density, colour=Harvest, shape=Treatment))+
 geom_errorbar(aes(ymin=Biomass-sd, ymax=Biomass+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
  facet_wrap(~FxGroup, ncol=1, scale="free")+
  geom_point(size=4.5,position=position_dodge(width=.65))+theme_classic()

# Regrowth
colnames(RegrowMeans)[4]<-"Regrow"
colnames(RegrowMeansG)[4]<-"Regrow"
colnames(RegrowMeansF)[4]<-"Regrow"
colnames(RegrowMeansW)[4]<-"Regrow"
colnames(RegrowMeans)[6]<-"Biomass"
colnames(RegrowMeansG)[6]<-"Biomass"
colnames(RegrowMeansF)[6]<-"Biomass"
colnames(RegrowMeansW)[6]<-"Biomass"

RegrowMeans$FxGroup<-"Total"
RegrowMeansG$FxGroup<-"Grasses"
RegrowMeansF$FxGroup<-"Forbs"
RegrowMeansW$FxGroup<-"Woody"


RegrowMeansALL<-rbind(RegrowMeans,RegrowMeansG,RegrowMeansF,RegrowMeansW)
RegrowMeansALL$FxGroup<-as.factor(RegrowMeansALL$FxGroup)
levels(RegrowMeansALL$FxGroup)<-c("Forbs","Grasses","Total","Woody")
RegrowMeansALL$FxGroup<- factor(RegrowMeansALL$FxGroup, levels = c("Total", "Grasses","Forbs","Woody"))

RegrowGraph<-ggplot(RegrowMeansALL, aes(y=Regrow, x=Boma.density, col=Harvest,shape=Treatment, size=Biomass))+
  geom_hline(yintercept =0, col="grey", linetype="dashed")+
  geom_errorbar(aes(ymin=Regrow-sd, ymax=Regrow+sd),width=.1,lwd=1,position=position_dodge(width=.65),show.legend=F)+
  geom_point(position=position_dodge(width=.65))+scale_radius(range=c(2,8))+
  facet_wrap(~FxGroup, ncol=1, scale="free")+
  theme_classic()

library(grid)
library(gridExtra)
library(egg)

egg::ggarrange(BioGraph,RegrowGraph, ncol = 2) # common.legend = T,legend="right")

######################################################################
#### Biomass regrowth - plotting data ####
########################################################################

# What is the spatial distribution of regrowth measurements
# Regrowth biomass
#nsreharvest3<-read.table("ProductivitySeason3.txt",header=T,sep="\t")

names(nsReharvest)
#plot(nsReharvest$boma_density)
#nsReharvest$Boma.density<-nsReharvest$boma_density
#nsReharvest$Boma.density[nsReharvest$Boma.density>1.75e-06]<-"High"
#nsReharvest$Boma.density[nsReharvest$Boma.density<1.75e-06]<-"Low"
#nsReharvest$Boma.density[nsReharvest$Boma.density==8.36311009231797e-07]<-"High"

# Set up as date
nsReharvest$Harvest.date<-as.Date(nsReharvest$Harvest.date,"%d.%m.%Y")
nsReharvest$Reharvest.date<-as.Date(nsReharvest$Reharvest.date,"%d.%m.%Y")

# Only need exclosure and control
levels(nsReharvest$Treatment)
nsReharvestb<-droplevels(nsReharvest[nsReharvest$Treatment=="Control" | nsReharvest$Treatment=="Exclosure",])
nsReharvestDUP<-nsReharvestb

# Livestock density and regrowth only - not the original biomass - average seperately
nsReharvestb<- droplevels(nsReharvestb[nsReharvestb$Harvest!="original",])
nsReharvestb$harvest_code<-as.factor(with(nsReharvestb, paste(Boma.density,Treatment, sep="-")))
levels(nsReharvestb$harvest_code) # 18

# Remove single reharvested quadrats from sampling
nsReharvestb<- droplevels(nsReharvestb[nsReharvestb$Harvest!="single",])

# Regrowth - Total Biomass
nsReharvestavg<-aggregate(TotalBiomass1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb, mean)
nsReharvestsem<-aggregate(TotalBiomass1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb,sd)
nsReharvestavg<-cbind(nsReharvestavg,nsReharvestsem[6])
colnames(nsReharvestavg)[6]<-"Biomass"
colnames(nsReharvestavg)[7]<-"se"
nsReharvestavg$fxgroup<-"Total"

# Regrowth - Grass Biomass
nsreharvest3Gavg<-aggregate(GrassNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb, mean)
nsreharvest3Gsem<-aggregate(GrassNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb,sd)
nsReharvestGavg<-cbind(nsreharvest3Gavg,nsreharvest3Gsem[6])
colnames(nsReharvestGavg)[6]<-"Biomass"
colnames(nsReharvestGavg)[7]<-"se"
nsReharvestGavg$fxgroup<-"Grass"

# Regrowth - Woody Biomass
nsreharvest3Wavg<-aggregate(DwarfShrubNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb, mean)
nsreharvest3Wsem<-aggregate(DwarfShrubNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb,sd)
nsReharvestWavg<-cbind(nsreharvest3Wavg,nsreharvest3Wsem[6])
colnames(nsReharvestWavg)[6]<-"Biomass"
colnames(nsReharvestWavg)[7]<-"se"
nsReharvestWavg$fxgroup<-"Woody"

# Regrowth - Herb + Climbers
nsReharvestb$HerbClimber1<-nsReharvestb$HerbNetReharvestBiomass1+nsReharvestb$ClimberNetReharvestBiomass1
nsreharvest3Havg<-aggregate(HerbClimber1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb, mean)
nsreharvest3Hsem<-aggregate(HerbClimber1~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb,sd)
nsReharvestHavg<-cbind(nsreharvest3Havg,nsreharvest3Hsem[6])
colnames(nsReharvestHavg)[6]<-"Biomass"
colnames(nsReharvestHavg)[7]<-"se"
nsReharvestHavg$fxgroup<-"Herb & Climber"

# Combine Total, Grass, Woody and Herb layers
nsReharvestAll<-rbind(nsReharvestavg,nsReharvestGavg,nsReharvestWavg,nsReharvestHavg)

# Original biomass - Total biomass
#nsReharvestb0<- droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
nsReharvestb0<- nsReharvestb
nsReharvestavg0<-aggregate(TotalBiomass0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0, mean)
nsReharvestsem0<-aggregate(TotalBiomass0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0,sd)
nsReharvestavg0<-cbind(nsReharvestavg0,nsReharvestsem0[6])
colnames(nsReharvestavg0)[6]<-"Biomass"
colnames(nsReharvestavg0)[7]<-"se"
nsReharvestavg0$fxgroup<-"Total"

# Original biomass - Grass biomass
nsReharvestavgG0<-aggregate(GrassNetReharvestBiomass0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0, mean)
nsReharvestsemG0<-aggregate(GrassNetReharvestBiomass0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0,sd)
nsReharvestavgG0<-cbind(nsReharvestavgG0,nsReharvestsemG0[6])
colnames(nsReharvestavgG0)[6]<-"Biomass"
colnames(nsReharvestavgG0)[7]<-"se"
nsReharvestavgG0$fxgroup<-"Grass"

# Original biomass - Woody biomass
nsReharvestavgW0<-aggregate(DwarfShrubNetReharvestBiomass0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0, mean)
nsReharvestsemW0<-aggregate(DwarfShrubNetReharvestBiomass0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0,sem)
nsReharvestavgW0<-cbind(nsReharvestavgW0,nsReharvestsemW0[6])
colnames(nsReharvestavgW0)[6]<-"Biomass"
colnames(nsReharvestavgW0)[7]<-"se"
nsReharvestavgW0$fxgroup<-"Woody"

# Original biomass - Herb biomass
nsReharvestb0$HerbClimber0<-nsReharvestb0$HerbNetReharvestBiomass0+nsReharvestb0$ClimberNetReharvestBiomass0
nsReharvestavgH0<-aggregate(HerbClimber0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0, mean)
nsReharvestsemH0<-aggregate(HerbClimber0~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb0,sd)
nsReharvestavgH0<-cbind(nsReharvestavgH0,nsReharvestsemH0[6])
colnames(nsReharvestavgH0)[6]<-"Biomass"
colnames(nsReharvestavgH0)[7]<-"se"
nsReharvestavgH0$fxgroup<-"Herb & Climber"

# Combine Total, Grass, Woody and Herb layers and 
nsReharvestAll0<-rbind(nsReharvestavg0,nsReharvestavgG0,nsReharvestavgW0,nsReharvestavgH0)

# Relevel factors
# Plant groups
nsReharvestAll$fxgroup<-as.factor(nsReharvestAll$fxgroup)
nsReharvestAll0$fxgroup<-as.factor(nsReharvestAll0$fxgroup)
nsReharvestAll$fxgroup<- factor(nsReharvestAll$fxgroup, levels = c("Total","Grass","Herb & Climber","Woody"))
nsReharvestAll0$fxgroup<- factor(nsReharvestAll0$fxgroup, levels = c("Total","Grass","Herb & Climber","Woody"))

# Boma.density
nsReharvestAll$Boma.density<-as.factor(nsReharvestAll$Boma.density)
nsReharvestAll0$Boma.density<-as.factor(nsReharvestAll0$Boma.density)
nsReharvestAll$Boma.density<- factor(nsReharvestAll$Boma.density, levels = c("Low","High"))
nsReharvestAll0$Boma.density<- factor(nsReharvestAll0$Boma.density, levels = c("Low","High"))
aggregate(Biomass~Reharvest.date+Boma.density+Treatment+fxgroup,nsReharvestAll,mean)

# Convert date to season
nsReharvestAll$Reharvest.date2<-as.factor(nsReharvestAll$Reharvest.date)
nsReharvestAll0$Reharvest.date2<-as.factor(nsReharvestAll0$Reharvest.date)
levels(nsReharvestAll$Reharvest.date2)<-c("Short I", "Long", "Short II")
levels(nsReharvestAll0$Reharvest.date2)<-c("Short I", "Long", "Short II")

aggregate(Biomass~Reharvest.date+Boma.density+Treatment+fxgroup,nsReharvestAll,mean)

# Harvest treatment code
nsReharvestAll$harvest_trt<-as.factor(with(nsReharvestAll, paste(Harvest,Treatment,Boma.density, sep="-")))
nsReharvestAll0$harvest_trt<-as.factor(with(nsReharvestAll0, paste(Harvest,Treatment,Boma.density, sep="-")))

# Total Biomass
ReHavTot<-nsReharvestAll[nsReharvestAll$fxgroup=="Total",]
ReHavTot0<-nsReharvestAll0[nsReharvestAll0$fxgroup=="Total",]
levels(ReHavTot$Boma.density)<-c("Far away","Close")  
levels(ReHavTot0$Boma.density)<-c("Far away","Close")  
levels(ReHavTot$Treatment)<-c("Open","Exclosed")  
levels(ReHavTot0$Treatment)<-c("Open","Exclosed") 

# Adding second axis for rainfall
rainavg<-aggregate(rain.mm~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb, mean)
rainsem<-aggregate(rain.mm~Harvest+Reharvest.date+Treatment+Boma.density+harvest_code,nsReharvestb,sd)
rainavg2<-cbind(rainavg,rainsem[6])
colnames(rainavg2)[6]<-"rain.mm"
colnames(rainavg2)[7]<-"se"
rainavg2$fxgroup<-"Total"

rainavg2$Reharvest.date2<-as.factor(rainavg2$Reharvest.date)
levels(rainavg2$Reharvest.date2)<-c("Short I", "Long", "Short II")
rainavg2$Boma.density<- factor(rainavg2$Boma.density, levels = c("Low","High"))
levels(rainavg2$Boma.density)<-c("Far away","Close")  
levels(rainavg2$Treatment)<-c("Open","Exclosed")  

scaleFactor <- mean(rainavg2$rain.mm,na.rm=T)/mean(ReHavTot$Biomass,na.rm=T)

# Position dodge
pd <- position_dodge(0.5)

# Filling code
ReHavTot$LivTrt<-as.factor(with(ReHavTot, paste(Boma.density , Treatment, sep="")))
ReHavTot0$LivTrt<-as.factor(with(ReHavTot0, paste(Boma.density , Treatment, sep="")))
rainavg2$LivTrt<-as.factor(with(rainavg2, paste(Boma.density , Treatment, sep="")))

# Total Biomass Only 
Regrow<-ggplot(ReHavTot, aes(x=Reharvest.date2, y=Biomass, group=Harvest,linetype=Harvest)) #shape=Treatment,colour=Boma.density,fill=LivTrt 

#Regrow<-Regrow+geom_smooth(data=rainavg2,aes(y=(rain.mm+100)/1.4),method="loess",span=.9, se=F,linetype="dotted",colour="dodgerblue1", size=1,show.legend=F, alpha=.76)
#Regrow<-Regrow+geom_point(data=rainavg2,aes(y=(rain.mm+100)/1.4),fill="dodgerblue1",colour="dodgerblue1", size=3.5,show.legend=F, size=2,alpha=.65)
Regrow<-Regrow+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
Regrow<-Regrow+geom_errorbar(aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),position=pd,stat = "identity",linetype="solid",width=.1,show.legend=F)
Regrow<-Regrow+geom_point(position=pd,stat = "identity",colour="grey50",fill="white",size=3.5, shape=21,stroke=1)
Regrow<-Regrow+facet_wrap(~Treatment+Boma.density, scale="fixed", ncol=2)
#Regrow<-Regrow+scale_colour_manual(values=c("grey50")) #"grey35"
#Regrow<-Regrow+scale_fill_manual(values=c("white"))#"black","white","grey70","white","grey35","white"#"grey60",
Regrow<-Regrow+scale_shape_manual(values=c(21,21))
Regrow<-Regrow+geom_errorbar(data=ReHavTot0,aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),width=.1, linetype="solid", colour="black",alpha=.99,show.legend=F)
Regrow<-Regrow+geom_point(data=ReHavTot0,size=3.5,stroke=1,alpha=.99,colour="black", fill="black",shape=22,show.legend=F)
Regrow<-Regrow+scale_linetype_manual(values =c("double" ="solid", single="dotted"))
Regrow<-Regrow+xlab("Season") + ylab(expression(paste("Biomass (g ",m^-2,")")))
#Regrow<-Regrow+scale_y_continuous(limits=c(0,275),sec.axis = sec_axis(~ . *1.4, breaks = c(100,200,300,400), labels = c(0,100,200,300), name = "Rainfall (mm)"),expand=c(0,0))
Regrow<-Regrow+ 
  theme(rect = element_rect(fill ="transparent")
        ,panel.background=element_rect(fill="transparent")
        ,plot.background=element_rect(fill="transparent",colour=NA)
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
        ,legend.key=element_rect(colour = NA, fill = NA)
        ,legend.key.width = unit(1.2,"cm"))
Regrow<- Regrow+annotate(geom = 'segment', y =-Inf, yend =-Inf, color = 'black', x =  -Inf, xend = Inf, size = .75) 
Regrow<- Regrow+annotate(geom = 'segment', x =-Inf, xend =-Inf, color = 'black', y =  -Inf, yend = Inf, size = .75) 

Regrow

#Regrow<- Regrow+guides(colour=F,size=F,alpha=F, fill=F, linetype=F,
#                       shape=guide_legend(order=2,"Treatment",override.aes = list(shape=c(21,21), alpha=0.99, size=4,linetype=NA,fill=c("grey30","white"),col=c("grey30"),stroke=1)))
#Regrow<- Regrow+guides(linetype=F, fill=F,#size=T,linetype=F,
#                         colour = guide_legend("Livestock density",override.aes = list(shape=c(21,21,21), size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1, linetype=NA)),
#                         shape = guide_legend("Treatment",override.aes = list(shape=c(21,21), size=3.5,fill=c("white","grey50"),col=c("grey30","grey30"), stroke=1, linetype=NA)))

#ReHavTot2<-ReHavTot
#ReHavTot2$Sampling<-as.factor(ReHavTot2$Treatment)
#levels(ReHavTot2$Sampling)<-c("Standing biomass","Regrowth") #"One harvest","Two Harvests","Three Harvests"
#Regrow2 <-  Regrow+ geom_point(data = ReHavTot2, aes(size=Sampling, shape = NA), colour = "grey50")
#Regrow2 <-  Regrow2 + guides(size=guide_legend("Harvests", override.aes=list(shape=c(22,21),linetype=c("blank","solid"), size=1,fill=c("orangered3","white"),col=c("black","grey20"), stroke=1))) # "blank","dotted","solid" # "black","white","grey60" #"black","grey30","grey20"
#Regrow2

# Enlarge symbol around line manually
#library(grid)
#library(ggpubr)
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-1.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-1.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-1.6-2-6-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-2.6-2-6-2", size = unit(4, "mm")) 

#Regrow2b <- grid.grab()
#is.grob(Regrow2b)

# Export graph
#filename <- paste0("TotBioRegrowth_NecSar", "_",Sys.Date(), ".jpeg" )
#jpeg (filename, width=22, height=12, res=400, unit="cm")
#Regrow2
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#Regrow2b <- grid.grab()
#dev.off()

# All other functional groups - Grass, Woody and Herbs - supporting info

#### Grasses only - Regrowth graph ####
ReHavG<-nsReharvestAll[nsReharvestAll$fxgroup=="Grass",]
ReHavG0<-nsReharvestAll0[nsReharvestAll0$fxgroup=="Grass",]
levels(ReHavG$Boma.density)<-c("Low density","High density")  
levels(ReHavG0$Boma.density)<-c("Low density","High density")  
levels(ReHavG$Treatment)<-c("Open","Exclosed")  
levels(ReHavG0$Treatment)<-c("Open","Exclosed")

# Filling code
ReHavG$LivTrt<-as.factor(with(ReHavG, paste(Boma.density , Treatment, sep="")))
ReHavG0$LivTrt<-as.factor(with(ReHavG0, paste(Boma.density , Treatment, sep="")))

#### Grass regrowth ####
RegrowG<-ggplot(ReHavG, aes(x=Reharvest.date2, y=Biomass,
                             group=Harvest,linetype=Harvest,shape=Treatment,colour=Boma.density,fill=LivTrt)) 
RegrowG<-RegrowG+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
RegrowG<-RegrowG+geom_errorbar(aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
RegrowG<-RegrowG+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
RegrowG<-RegrowG+facet_wrap(~Treatment+Boma.density, scale="fixed", ncol=2)
RegrowG<-RegrowG+scale_colour_manual(values=c("grey70","black")) #"grey35"
RegrowG<-RegrowG+scale_fill_manual(values=c("black","white","grey70","white","grey35","white"))#"grey60",
RegrowG<-RegrowG+scale_shape_manual(values=c(21,21))
RegrowG<-RegrowG+geom_errorbar(data=ReHavG0,aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),width=.2, linetype="solid", colour="black",alpha=.99,show.legend=F)
RegrowG<-RegrowG+geom_point(data=ReHavG0,size=3.5,stroke=1,alpha=.99,colour="black", fill="orangered3",shape=22,show.legend=F)
RegrowG<-RegrowG+scale_linetype_manual(values =c("double" ="solid", single="dotted"))
RegrowG<-RegrowG+xlab("Rainy season") + ylab(expression(paste("Biomass (g ",m^-2,")")))
RegrowG<-RegrowG+ 
  theme(rect = element_rect(fill ="transparent")
        ,panel.background=element_rect(fill="transparent")
        ,plot.background=element_rect(fill="transparent",colour=NA)
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
        ,legend.key=element_rect(colour = NA, fill = NA)
        ,legend.key.width = unit(1.2,"cm"))
RegrowG<-RegrowG+annotate(geom = 'segment', y =-Inf, yend =-Inf, color = 'black', x =  -Inf, xend = Inf, size = .75) 
RegrowG<-RegrowG+annotate(geom = 'segment', x =-Inf, xend =-Inf, color = 'black', y =  -Inf, yend = Inf, size = .75) 

#Regrow<- Regrow+guides(colour=F,size=F,alpha=F, fill=F, linetype=F,
#                       shape=guide_legend(order=2,"Treatment",override.aes = list(shape=c(21,21), alpha=0.99, size=4,linetype=NA,fill=c("grey30","white"),col=c("grey30"),stroke=1)))
#RegrowG<-RegrowG+guides(linetype=F, fill=F,#size=T,linetype=F,
#                       colour = guide_legend("Livestock density",override.aes = list(shape=c(21,21,21), size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1, linetype=NA)),
#                       shape = guide_legend("Treatment",override.aes = list(shape=c(21,21), size=3.5,fill=c("white","grey50"),col=c("grey30","grey30"), stroke=1, linetype=NA)))
#RegrowG

#ReHavG2<-ReHavG
#ReHavG2$Sampling<-as.factor(ReHavG2$Treatment)
#levels(ReHavG2$Sampling)<-c("Standing biomass","Regrowth") #"One harvest","Two Harvests","Three Harvests"
#RegrowG2 <-  RegrowG+ geom_point(data = ReHavG2, aes(size=Sampling, shape = NA), colour = "grey50")
#RegrowG2 <-  RegrowG2 + guides(size=guide_legend("Harvests", override.aes=list(shape=c(22,21),linetype=c("blank","solid"), size=1,fill=c("orangered3","white"),col=c("black","grey20"), stroke=1))) # "blank","dotted","solid" # "black","white","grey60" #"black","grey30","grey20"
#RegrowG2

# Enlarge symbol around line manually
#library(grid)
#library(ggpubr)
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-1.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-1.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-1.6-2-6-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-2.6-2-6-2", size = unit(4, "mm")) 

#RegrowG2b <- grid.grab()
#is.grob(RegrowG2b)

# Export graph
#filename <- paste0("GrassBioRegrowth_NecSar", "_",Sys.Date(), ".jpeg" )
#jpeg (filename, width=22, height=12, res=400, unit="cm")
#RegrowG2
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#RegrowG2b <- grid.grab()
#dev.off()

#### Herbs and climbers only - Regrowth graph ####
ReHavH<-nsReharvestAll[nsReharvestAll$fxgroup=="Herb & Climber",]
ReHavH0<-nsReharvestAll0[nsReharvestAll0$fxgroup=="Herb & Climber",]
levels(ReHavH$Boma.density)<-c("Low density","High density")  
levels(ReHavH0$Boma.density)<-c("Low density","High density")  
levels(ReHavH$Treatment)<-c("Open","Exclosed")  
levels(ReHavH0$Treatment)<-c("Open","Exclosed") 

# Filling code
ReHavH$LivTrt<-as.factor(with(ReHavH, paste(Boma.density , Treatment, sep="")))
ReHavH0$LivTrt<-as.factor(with(ReHavH0, paste(Boma.density , Treatment, sep="")))

# Herb regrowth graph
RegrowH<-ggplot(ReHavH, aes(x=Reharvest.date2, y=Biomass,
                            group=Harvest,linetype=Harvest,shape=Treatment,colour=Boma.density,fill=LivTrt)) 
RegrowH<-RegrowH+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
RegrowH<-RegrowH+geom_errorbar(aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
RegrowH<-RegrowH+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
RegrowH<-RegrowH+facet_wrap(~Treatment+Boma.density, scale="fixed", ncol=2)
RegrowH<-RegrowH+scale_colour_manual(values=c("grey70","black")) #"grey35"
RegrowH<-RegrowH+scale_fill_manual(values=c("black","white","grey70","white","grey35","white"))#"grey60",
RegrowH<-RegrowH+scale_shape_manual(values=c(21,21))
RegrowH<-RegrowH+geom_errorbar(data=ReHavH0,aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),width=.2, linetype="solid", colour="black",alpha=.99,show.legend=F)
RegrowH<-RegrowH+geom_point(data=ReHavH0,size=3.5,stroke=1,alpha=.99,colour="black", fill="orangered3",shape=22,show.legend=F)
RegrowH<-RegrowH+scale_linetype_manual(values =c("double" ="solid", single="dotted"))
RegrowH<-RegrowH+xlab("Rainy season") + ylab(expression(paste("Biomass (g ",m^-2,")")))
RegrowH<-RegrowH+
  theme(rect = element_rect(fill ="transparent")
        ,panel.background=element_rect(fill="transparent")
        ,plot.background=element_rect(fill="transparent",colour=NA)
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
        ,legend.key=element_rect(colour = NA, fill = NA)
        ,legend.key.width = unit(1.2,"cm"))
RegrowH<-RegrowH+annotate(geom = 'segment', y =-Inf, yend =-Inf, color = 'black', x =  -Inf, xend = Inf, size = .75) 
RegrowH<-RegrowH+annotate(geom = 'segment', x =-Inf, xend =-Inf, color = 'black', y =  -Inf, yend = Inf, size = .75) 

#Regrow<- Regrow+guides(colour=F,size=F,alpha=F, fill=F, linetype=F,
#                       shape=guide_legend(order=2,"Treatment",override.aes = list(shape=c(21,21), alpha=0.99, size=4,linetype=NA,fill=c("grey30","white"),col=c("grey30"),stroke=1)))
#RegrowH<-RegrowH+guides(linetype=F, fill=F,#size=T,linetype=F,
#                        colour = guide_legend("Livestock density",override.aes = list(shape=c(21,21,21), size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1, linetype=NA)),
#                        shape = guide_legend("Treatment",override.aes = list(shape=c(21,21), size=3.5,fill=c("white","grey50"),col=c("grey30","grey30"), stroke=1, linetype=NA)))

#ReHavH2<-ReHavH
#ReHavH2$Sampling<-as.factor(ReHavH2$Treatment)
#levels(ReHavH2$Sampling)<-c("Standing biomass","Regrowth") #"One harvest","Two Harvests","Three Harvests"
#RegrowH2 <-  RegrowH+ geom_point(data = ReHavH2, aes(size=Sampling, shape = NA), colour = "grey50")
#RegrowH2 <-  RegrowH2 + guides(size=guide_legend("Harvests", override.aes=list(shape=c(22,21),linetype=c("blank","solid"), size=1,fill=c("orangered3","white"),col=c("black","grey20"), stroke=1))) # "blank","dotted","solid" # "black","white","grey60" #"black","grey30","grey20"
#RegrowH2

# Enlarge symbol around line manually
#library(grid)
#library(ggpubr)
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-1.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-1.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-1.6-2-6-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-2.6-2-6-2", size = unit(4, "mm")) 

#RegrowG2b <- grid.grab()
#is.grob(RegrowG2b)

# Export graph
#filename <- paste0("HerbsBioRegrowth_NecSar", "_",Sys.Date(), ".jpeg" )
#jpeg (filename, width=22, height=12, res=400, unit="cm")
#RegrowH2
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#RegrowH2b <- grid.grab()
#dev.off()

#### Woody only - Regrowth graph ####
ReHavW<-nsReharvestAll[nsReharvestAll$fxgroup=="Woody",]
ReHavW0<-nsReharvestAll0[nsReharvestAll0$fxgroup=="Woody",]
levels(ReHavW$Boma.density)<-c("Low density","High density")  
levels(ReHavW0$Boma.density)<-c("Low density","High density")  
levels(ReHavW$Treatment)<-c("Open","Exclosed")  
levels(ReHavW0$Treatment)<-c("Open","Exclosed") 

# Filling code
ReHavW$LivTrt<-as.factor(with(ReHavW, paste(Boma.density , Treatment, sep="")))
ReHavW0$LivTrt<-as.factor(with(ReHavW0, paste(Boma.density , Treatment, sep="")))

# Woody regrowth
RegrowW<-ggplot(ReHavW, aes(x=Reharvest.date2, y=Biomass,
                            group=Harvest,linetype=Harvest,shape=Treatment,colour=Boma.density,fill=LivTrt)) 
RegrowW<-RegrowW+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
RegrowW<-RegrowW+geom_errorbar(aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
RegrowW<-RegrowW+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
RegrowW<-RegrowW+facet_wrap(~Treatment+Boma.density, scale="fixed", ncol=2)
RegrowW<-RegrowW+scale_colour_manual(values=c("grey70","black")) #"grey35"
RegrowW<-RegrowW+scale_fill_manual(values=c("black","white","grey70","white","grey35","white"))#"grey60",
RegrowW<-RegrowW+scale_shape_manual(values=c(21,21))
RegrowW<-RegrowW+geom_errorbar(data=ReHavW0,aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),width=.2, linetype="solid", colour="black",alpha=.99,show.legend=F)
RegrowW<-RegrowW+geom_point(data=ReHavW0,size=3.5,stroke=1,alpha=.99,colour="black", fill="orangered3",shape=22,show.legend=F)
RegrowW<-RegrowW+scale_linetype_manual(values =c("double" ="solid", single="dotted"))
RegrowW<-RegrowW+xlab("Rainy season") + ylab(expression(paste("Biomass (g ",m^-2,")")))
RegrowW<-RegrowW+
  theme(rect = element_rect(fill ="transparent")
        ,panel.background=element_rect(fill="transparent")
        ,plot.background=element_rect(fill="transparent",colour=NA)
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
        ,legend.key=element_rect(colour = NA, fill = NA)
        ,legend.key.width = unit(1.2,"cm"))
RegrowW<-RegrowW+annotate(geom = 'segment', y =-Inf, yend =-Inf, color = 'black', x =  -Inf, xend = Inf, size = .75) 
RegrowW<-RegrowW+annotate(geom = 'segment', x =-Inf, xend =-Inf, color = 'black', y =  -Inf, yend = Inf, size = .75) 

#Regrow<- Regrow+guides(colour=F,size=F,alpha=F, fill=F, linetype=F,
#                       shape=guide_legend(order=2,"Treatment",override.aes = list(shape=c(21,21), alpha=0.99, size=4,linetype=NA,fill=c("grey30","white"),col=c("grey30"),stroke=1)))
#RegrowW<-RegrowW+guides(linetype=F, fill=F,#size=T,linetype=F,
#                        colour = guide_legend("Livestock density",override.aes = list(shape=c(21,21,21), size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1, linetype=NA)),
#                        shape = guide_legend("Treatment",override.aes = list(shape=c(21,21), size=3.5,fill=c("white","grey50"),col=c("grey30","grey30"), stroke=1, linetype=NA)))

#ReHavW2<-ReHavW
#ReHavW2$Sampling<-as.factor(ReHavW2$Treatment)
#levels(ReHavW2$Sampling)<-c("Standing biomass","Regrowth") #"One harvest","Two Harvests","Three Harvests"
#RegrowW2 <-  RegrowW+ geom_point(data = ReHavW2, aes(size=Sampling, shape = NA), colour = "grey50")
#RegrowW2 <-  RegrowW2 + guides(size=guide_legend("Harvests", override.aes=list(shape=c(22,21),linetype=c("blank","solid"), size=1,fill=c("orangered3","white"),col=c("black","grey20"), stroke=1))) # "blank","dotted","solid" # "black","white","grey60" #"black","grey30","grey20"
#RegrowW2

# Enlarge symbol around line manually
#library(grid)
#library(ggpubr)
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-1.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-1.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-1.6-2-6-2", size = unit(4, "mm")) 
#grid.gedit("key-5-1-2.6-2-6-2", size = unit(4, "mm")) 

#RegrowW2b <- grid.grab()
#is.grob(RegrowW2b)

# Export graph
#filename <- paste0("WoodyBioRegrowth_NecSar", "_",Sys.Date(), ".jpeg" )
#jpeg (filename, width=22, height=12, res=400, unit="cm")
#RegrowW2
#grid.ls(grid.force()) 
#grid.gedit("key-3-1-2.4-2-4-2", size = unit(4, "mm")) 
#grid.gedit("key-4-1-2.5-2-5-2", size = unit(4, "mm")) 
#RegrowW2b <- grid.grab()
#dev.off()

RegrowW<-ggplot(ReHavW, aes(x=Reharvest.date2, y=Biomass,
                            group=Harvest,linetype=Harvest,shape=Treatment,colour=Harvest,fill=Harvest)) 
RegrowW<-RegrowW+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
RegrowW<-RegrowW+geom_errorbar(aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
RegrowW<-RegrowW+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
RegrowW<-RegrowW+facet_wrap(~Treatment+Livestock.density, scale="fixed", ncol=3)
RegrowW<-RegrowW+scale_colour_manual(values=c("grey20","grey30"))
RegrowW<-RegrowW+scale_fill_manual(values=c("grey60","white"))
RegrowW<-RegrowW+scale_shape_manual(values=c(21,22))
RegrowW<-RegrowW+geom_errorbar(data=ReHavW0,aes(x = Reharvest.date2, ymin=Biomass-se,ymax=Biomass+se),width=.2, linetype="solid", colour="black",alpha=.99,show.legend=F)
RegrowW<-RegrowW+geom_point(data=ReHavW0,size=3.5,stroke=1,alpha=.99,colour="black", fill="black",show.legend=F)
RegrowW<-RegrowW+scale_linetype_manual(values =c("double" ="solid", "single"="dotted"))
RegrowW<-RegrowW+xlab("Rainy season") + ylab(expression(paste("Biomass (g ",m^-2,")"))) + ggtitle("Woody")
RegrowW<-RegrowW+#theme_bw() +
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
        ,legend.key=element_rect(colour = NA, fill = NA)
        ,legend.key.width = unit(1.2,"cm"))
RegrowW<-RegrowW+annotate(geom = 'segment', y =-Inf, yend =-Inf, color = 'black', x =  -Inf, xend = Inf, size = .75) 
RegrowW<-RegrowW+annotate(geom = 'segment', x =-Inf, xend =-Inf, color = 'black', y =  -Inf, yend = Inf, size = .75) 

#RegrowW<-RegrowW+guides(colour=F,size=F,alpha=F, fill=F, linetype=F,
#                        shape=guide_legend(order=2,"Treatment",override.aes = list(shape=c(21,22), alpha=0.99, size=4,linetype=NA,fill=NA,col=c("grey30"),stroke=1)))
#ReHavW2<-ReHavW
#ReHavW2$Sampling<-as.factor(ReHavW2$Livestock.density)
#levels(ReHavW2$Sampling)<-c("One harvest","Two Harvests","Three Harvests")
#RegrowW2 <-  RegrowW+ geom_point(data = ReHavW2, aes(size=Sampling, shape = NA, linetype=NA), colour = "grey50")
#RegrowW2 <-  RegrowW2 + guides(size=guide_legend("Harvests", override.aes=list(shape=c(21),linetype=c("blank","dotted","solid"), size=1,fill=c("black","white","grey60"),col=c("black","grey30","grey20"), stroke=1)))
#RegrowW2

# Combine graphs
#library(grid)
#library(gridExtra)
#library(egg)
#library(ggpubr)
#egg::ggarrange(Regrow+ theme(legend.position="none"),Regrow2+ theme(legend.position="none"), ncol=1) #common.legend = T)

#### Functional group split ####
mean(nsReharvestb$GrassNetReharvestBiomass1/nsReharvestb$TotalBiomass1*100,na.rm=T)
sd(nsReharvestb$GrassNetReharvestBiomass1/nsReharvestb$TotalBiomass1*100,na.rm=T)
# Grasses 82.86805 + 23
nsReharvestb$HerbClimberBiomass1<-nsReharvestb$ClimberNetReharvestBiomass1 + nsReharvestb$HerbNetReharvestBiomass1

mean(nsReharvestb$HerbClimberBiomass1/nsReharvestb$TotalBiomass1*100,na.rm=T)
sd(nsReharvestb$HerbClimberBiomass1/nsReharvestb$TotalBiomass1*100,na.rm=T)
#Herb and climber 9.214427 + 17

mean(nsReharvestb$DwarfShrubNetReharvestBiomass1/nsReharvestb$TotalBiomass1*100,na.rm=T)
sd(nsReharvestb$DwarfShrubNetReharvestBiomass1/nsReharvestb$TotalBiomass1*100,na.rm=T)
#Herb and climber 7.917522 + 16

#### Functional group averages ####
# Combine Total, Grass, Woody and Herb layers and 
nsReharvestAll0 # Total biomass
nsReharvestAll # Regrowth
nsReharvestAllsub<-droplevels(nsReharvestAll[nsReharvestAll$Reharvest.date!="2012-11-25",])

nsReharvestAllsub$Harvest<-"regrowth"
nsReharvestAll0$Harvest<-"total"

nsReharvestFxGrp<-rbind(nsReharvestAllsub,nsReharvestAll0)

nsReharvestAllavg<-aggregate(Biomass~Harvest+Boma.density+Treatment+fxgroup,nsReharvestFxGrp, mean)
nsReharvestAllsd<-aggregate(Biomass~Harvest+Boma.density+Treatment+fxgroup,nsReharvestFxGrp, sd)

nsReharvestAllavg$sd<-nsReharvestAllsd$Biomass

ggplot(nsReharvestAllavg, aes(y=Biomass, x=Boma.density, colour=fxgroup, shape=Treatment))+
  geom_errorbar(aes(x=Boma.density,ymin=Biomass-sd,ymax=Biomass+sd),width=.1,position=pd)+
  geom_point(position=pd,size=3.5)+facet_wrap(~Harvest)+theme_bw()

########################################################################
#### Percentage difference ####
########################################################################

#### Percent difference Original biomass ####
nsReharvestbO<-droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
nsReharvestbH1_2<-droplevels(nsReharvestb[!is.na(nsReharvestb$Harvest.date),])
#nsReharvestbH1_2$TotalBiomass0<-nsReharvestbH1_2[nsReharvestbH1_2$TotalBiomass0<0.1,]<-.1

# Percent of biomass at each time interval - treatment
names(nsReharvestbH1_2)
nsReharvestbH1_2$Tot.Per.diff<-(nsReharvestbH1_2$TotalBiomass1-nsReharvestbH1_2$TotalBiomass0)/sd(nsReharvestbH1_2$TotalBiomass0)
nsReharvestbH1_2$Woody.Per.diff<-(nsReharvestbH1_2$DwarfShrubNetReharvestBiomass1-nsReharvestbH1_2$DwarfShrubNetReharvestBiomass0)/sd(nsReharvestbH1_2$DwarfShrubNetReharvestBiomass0)

#nsReharvestbH1_2$G.Per.diff<-((nsReharvestbH1_2$GrassNetReharvestBiomass1-nsReharvestbO$GrassNetReharvestBiomass1)/nsReharvestbO$GrassNetReharvestBiomass1)*100
#nsReharvestbH1_2$W.Per.diff<-((nsReharvestbH1_2$DwarfShrubNetReharvestBiomass1-nsReharvestbO$DwarfShrubNetReharvestBiomass1)/nsReharvestbO$DwarfShrubNetReharvestBiomass1)*100

# Percent diff total biomass compared to the original plot
nsReharvestavg<-aggregate(Tot.Per.diff~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestbH1_2, mean)
nsReharvestsem<-aggregate(Tot.Per.diff~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestbH1_2,sd)
nsReharvestavg<-cbind(nsReharvestavg,nsReharvestsem[6])
colnames(nsReharvestavg)[7]<-"se"

# Total biomassregrowth
Regrow2<-ggplot(nsReharvestavg, aes(x=Reharvest.date, y=Tot.Per.diff,linetype=harvest_code,shape=Treatment,colour=Harvest,fill=Harvest)) 
Regrow2<-Regrow2+facet_wrap(~Livestock.density, scale="fixed")
Regrow2<-Regrow2+geom_line(aes(linetype=harvest_code), show.legend = F) +theme_classic()
Regrow2<-Regrow2+geom_point(size=4)
Regrow2<-Regrow2+ggtitle("Total Biomass percentage")
Regrow2

# Percent diff woody biomass # compared to the original plot
nsReharvestavgW<-aggregate(Woody.Per.diff~Harvest+Reharvest.date+Treatment+Livestock.density,nsReharvestbH1_2, mean)
nsReharvestsemW<-aggregate(Woody.Per.diff~Harvest+Reharvest.date+Treatment+Livestock.density,nsReharvestbH1_2,sd)

# Difference between exclosed and open, livestock density for each season
aggregate(DwarfShrubNetReharvestBiomass1~Livestock.density+Treatment+Season,nsReharvestb,mean)
((8.520000-2.100000)/2.100000)*100

# Difference between exclosed and open by livestock density
aggregate(TotalBiomass1~Livestock.density+Treatment,nsReharvestb,mean)
((127.89556-108.96000)/108.96000)*100 # 19 % difference # LOW
((93.22889-55.90667)/55.90667)*100 # 66% medium 
((68.86444-39.11333)/39.11333)*100 # 76 % high 

########################################################################
#### Data exploration ####
# A Missing values?
# B Outliers in Y / Outliers in X
# C Collinearity X
# D Relationships Y vs X
# E Spatial/temporal aspects of sampling design (not relevant here)
# F Interactions (is the quality of the data good enough to include them?)
# G Zero inflation Y
# H Are categorical covariates balanced?
########################################################################

# Only need exclosure and control
nsReharvestb<-droplevels(nsReharvest[nsReharvest$Treatment=="Control" | nsReharvest$Treatment=="Exclosure",])

# Livestock density and regrowth
nsReharvestb$harvest_code<-as.factor(with(nsReharvestb, paste(Boma.density,Treatment, sep="-")))
levels(nsReharvestb$harvest_code) # 4 codes Boma-exclosure treatment

#### Only reharvest ####
nsReharvestbO<-droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
nsReharvestbH1_2<-droplevels(nsReharvestb[!is.na(nsReharvestb$Harvest.date),])

nsReharvestbH1_2$Plot.ID<-as.numeric(nsReharvestbH1_2$Plot.name)
Rdate2<-strptime(as.character(nsReharvestbH1_2$Reharvest.date),format="%d.%m.%Y",tz="Africa/Nairobi" )
nsReharvestbH1_2$YrMonth<-format(Rdate2, "%Y-%m")

# ONLY double harvest...
nsReharvestH1<- droplevels(nsReharvestbH1_2[nsReharvestbH1_2$Harvest=="double",])
#nsReharvestH1<-nsReharvestbH1_2
nsReharvestH1$Harvest

# Block as factor
nsReharvestH1$fBlock<-as.factor(nsReharvestH1$Block)
nsReharvestH1$fPlot.ID<-as.factor(nsReharvestH1$Plot.ID)

#B Collinearity X

# Herbivore covariates
#TotalHerb<-cbind(nsReharvestb$Burchells_ZebraMetBio,
#nsReharvestb$CattleMetBio,nsReharvestb$Grants_GazelleMetBio,
#nsReharvestb$Greater_KuduMetBio,nsReharvestb$Swaynes_HartebeestMetBio)
#WildHerb<-cbind(nsReharvestb$Burchells_ZebraMetBio,nsReharvestb$Grants_GazelleMetBio,
#                 nsReharvestb$Greater_KuduMetBio,nsReharvestb$Swaynes_HartebeestMetBio)
#nsReharvestb$TotalHerb<-rowSums(TotalHerb, na.rm=T)
#nsReharvestb$WildHerb<-rowSums(WildHerb, na.rm=T)

# Boxplot for factors
names(nsReharvestb)
par(mfrow = c(1, 1), mar = c(4, 3, 3, 2))
boxplot(rain.mm ~ Season, 
        xlab = "Season",
        ylab = "Rain.mm",
        data = nsReharvestb) # Collinear...

par(mfrow = c(1, 1), mar = c(4, 3, 3, 2))
boxplot(rain.mm ~ Livestock.density, 
        xlab = "Livestock.density",
        ylab = "Rain.mm",
        data = nsReharvestb) # Not colinear - low area not drier

bwplot(rain.mm~Season|Livestock.density,nsReharvestb)
# Second short season - medium livestock gets a bit more rainfall...

#corvif
library(car)
names(nsReharvestb)
GRAZE.lm1<-lm(TotalBiomass1~Treatment+Harvest.date+rain.mm+
                Burchells_ZebraMetBio+CattleMetBio+Grants_GazelleMetBio+
                Greater_KuduMetBio+Swaynes_HartebeestMetBio, data = nsReharvestb)
vif(GRAZE.lm1)
# Cannot use rainfall and season in the same model
#Treatment             Harvest.date                  rain.mm    Burchells_ZebraMetBio 
#1.000000              1727.356013              1702.018066                 1.859672 
#CattleMetBio     Grants_GazelleMetBio       Greater_KuduMetBio Swaynes_HartebeestMetBio 
#2.073992                 1.523146                 2.387375                 1.486234 mm           

# C Collinearity Y
names(nsReharvestb)
MyVar<-c("GrassNetReharvestBiomass1","DwarfShrubNetReharvestBiomass1",
         "HerbNetReharvestBiomass1","ClimberNetReharvestBiomass1",
         "Burchells_ZebraMetBio","CattleMetBio","Grants_GazelleMetBio","Greater_KuduMetBio",
         "Swaynes_HartebeestMetBio" )
pairs(nsReharvestb[,MyVar],lower.panel = panel.cor)
# Poor correlation between functional group regrowth biomass

########################################################################
#### Biomass regrowth - data analysis ####
library(nlme)
library(lme4)
########################################################################

# Only need exclosure and control
nsReharvestb<-droplevels(nsReharvest[nsReharvest$Treatment=="Control" | nsReharvest$Treatment=="Exclosure",])

# Livestock density and regrowth
nsReharvestb$harvest_code<-as.factor(with(nsReharvestb, paste(Boma.density,Treatment, sep="-")))
levels(nsReharvestb$harvest_code) # 4

#### Only reharvest ####
nsReharvestbO<-droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
nsReharvestbH1_2<-droplevels(nsReharvestb[!is.na(nsReharvestb$Harvest.date),])

nsReharvestbH1_2$Plot.ID<-as.numeric(nsReharvestbH1_2$Plot.name)
Rdate2<-strptime(as.character(nsReharvestbH1_2$Reharvest.date),format="%d.%m.%Y",tz="Africa/Nairobi" )
nsReharvestbH1_2$YrMonth<-format(Rdate2, "%Y-%m")

# ONLY double harvest...
nsReharvestH1<- droplevels(nsReharvestbH1_2[nsReharvestbH1_2$Harvest=="double",])
#nsReharvestH1<-nsReharvestbH1_2
nsReharvestH1$Harvest

# Block as factor
nsReharvestH1$fBlock<-as.factor(nsReharvestH1$Block)
nsReharvestH1$fPlot.ID<-as.factor(nsReharvestH1$Plot.ID)

nsReharvestbH1_2$fBlock<-as.factor(nsReharvestbH1_2$Block)
nsReharvestbH1_2$fPlot.ID<-as.factor(nsReharvestbH1_2$Plot.ID)

# Testing whether auto-correlation structure is necessary
# Auto correlation structure - unneccessary
#cs1AR1 <- corAR1(0.2, form = ~Reharvest.date|fBlock/fPlot.ID)
#cs1AR1. <- Initialize(cs1AR1, data =nsReharvestH1) 
#corMatrix(cs1AR1.)
nsReharvestH1T<- nsReharvestH1[is.finite(nsReharvestH1$Tot.Per.diff),]

# Herbivore covariates
names(nsReharvestH1)
TotalHerb<-cbind(nsReharvestH1$Burchells_ZebraMetBio,
                 nsReharvestH1$CattleMetBio,nsReharvestH1$Grants_GazelleMetBio,
                 nsReharvestH1$Greater_KuduMetBio,nsReharvestH1$Swaynes_HartebeestMetBio)
WildHerb<-cbind(nsReharvestH1$Burchells_ZebraMetBio,nsReharvestH1$Grants_GazelleMetBio,
                nsReharvestH1$Greater_KuduMetBio,nsReharvestH1$Swaynes_HartebeestMetBio)
nsReharvestH1$TotalHerb<-rowSums(TotalHerb, na.rm=T)
nsReharvestH1$WildHerb<-rowSums(WildHerb, na.rm=T)

########################################################################################################
#### TOTAL BIOMASS REGROWTH MODEL ####
########################################################################################################

# Only need exclosure and control
nsReharvestb<-droplevels(nsReharvest[nsReharvest$Treatment=="Control" | nsReharvest$Treatment=="Exclosure",])
nsReharvestb$Harvest # Single and Double harvest

# Livestock density and regrowth
nsReharvestb$harvest_code<-as.factor(with(nsReharvestb, paste(Boma.density,Treatment, sep="-")))
levels(nsReharvestb$harvest_code) # 4

#### Only reharvest - REMOVES first harvest (not reharvest) ####
#nsReharvestbO<-droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
#nsReharvestbH1_2<-droplevels(nsReharvestb[!is.na(nsReharvestb$Harvest.date),])

#nsReharvestbH1_2$Plot.ID<-as.numeric(nsReharvestbH1_2$Plot.name)
#Rdate2<-strptime(as.character(nsReharvestbH1_2$Reharvest.date),format="%d.%m.%Y",tz="Africa/Nairobi" )
#nsReharvestbH1_2$YrMonth<-format(Rdate2, "%Y-%m")

#nsReharvestb$Plot.ID<-as.numeric(nsReharvestb$Plot.name)
#Rdate3<-strptime(as.character(nsReharvestb$Reharvest.date),format="%Y-%m-%d",tz="Africa/Nairobi" )
#nsReharvestb$YrMonth<-format(Rdate3, "%Y-%m")

# ONLY double harvest...REMOVES FIRST HARVEST
nsReharvestH1<- droplevels(nsReharvestb[nsReharvestb$Harvest=="double",])
#nsReharvestH1<-nsReharvestbH1_2
nsReharvestH1$Reharvest.date
str(nsReharvestH1)
#nsReharvestbO<-droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
#nsReharvestbH1_2<-droplevels(nsReharvestb[!is.na(nsReharvestb$Harvest.date),])

# Duplicates of harvest 1 and 2 in time = issues with temporal model
nsReharvestH1$Plot.ID<-as.numeric(nsReharvestH1$Plot.name)
Rdate2<-strptime(as.character(nsReharvestH1$Reharvest.date),format="%d.%m.%Y",tz="Africa/Nairobi" )
nsReharvestH1$YrMonth<-format(Rdate2, "%Y-%m")
#nsReharvestH1$YrMonth<-format(as.Date(nsReharvestH1$Reharvest.date), "%Y.%m")
nsReharvestH1$Reharvest.date2<-as.factor(nsReharvestH1$Reharvest.date)
levels(nsReharvestH1$Reharvest.date2)<-c("Short II", "Long", "Short I")
nsReharvestH1$Reharvest.date2 <- ordered(nsReharvestH1$Reharvest.date2, levels=c("Short I", "Long", "Short II"))

#nsReharvestb$Plot.ID<-as.numeric(nsReharvestb$Plot.name)
#nsReharvestb$YrMonth<-format(as.Date(nsReharvestb$Reharvest.date), "%Y-%m")
#nsReharvestb$Reharvest.date2<-as.factor(nsReharvestb$Reharvest.date)
#levels(nsReharvestb$Reharvest.date2)<-c("Short I", "Long", "Short II")

#nsReharvestH1<- droplevels(nsReharvestbH1_2[nsReharvestbH1_2$Harvest=="double",])
#table(nsReharvestH1$YrMonth,nsReharvestH1$Plot.ID)
#table(nsReharvestH1$Block,nsReharvestH1$Plot.ID)

#### Housekeeping ####
nsReharvestH1$fBlock<-as.factor(nsReharvestH1$Block)
nsReharvestH1$fPlot.ID<-as.factor(nsReharvestH1$Plot.ID)
nsReharvestH1$fBoma.density<-as.factor(nsReharvestH1$Boma.density)

nsReharvestb$fBlock<-as.factor(nsReharvestb$Block)
nsReharvestb$fPlot.ID<-as.factor(nsReharvestb$Plot.ID)
nsReharvestb$fBoma.density<-as.factor(nsReharvestb$Boma.density)
nsReharvestb$Plotcode<-as.factor(with(nsReharvestb, paste(fBlock,Treatment, sep="_")))

#### Regrowth Model  - only double harvest ####

# Lme auto-correlation
# Frequentist - time series
# turn a date into a 'monthnumber' relative to an origin
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01"));  lt$year*12 + lt$mon } 
# compute a month difference as a difference between two monnb's
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }
# take it for a spin
mondf(as.Date("2008-01-01"), Sys.Date())
#[1] 24
mondf(c("2002-03-31", "2002-04-30", "2002-05-31"), "2002-06-30")
nsReharvestH1$Reharvest.date
nsReharvestH1$YrMonthNumber<-mondf(c(as.POSIXlt(as.Date(nsReharvestH1$Reharvest.date,format="%Y-%m-%d",tz="Africa/Nairobi" ))), "2012-11-25")*-1
class(nsReharvestH1$YrMonthNumber)

#### Total biomass temporal autocorrelation model  ####
cs1AR1 <- corAR1(0.2, form = ~YrMonthNumber|fBlock/fPlot.ID) # AR matrix needs to be unique
cs1AR1. <- Initialize(cs1AR1, data = nsReharvestH1) 
corMatrix(cs1AR1.)

# Standardize covariate
MyStd <- function(x) {  (x - mean(x))  /  sd(x) } # Standardize covariates....

nsReharvestH1$TotalBiomass0.c     <- MyStd(nsReharvestH1$TotalBiomass0)

plot(TotalBiomass0~TotalBiomass0.c,nsReharvestH1)

Tot0<-lme(TotalBiomass0~Treatment+fBoma.density+Treatment:fBoma.density,
        random= ~ 1|fBlock/fPlot.ID, na.action=na.pass, method="ML",
        correlation=corAR1(0.2, form=~YrMonthNumber|fBlock/fPlot.ID),data=nsReharvestH1)
summary(Tot0)
anova(Tot0) # Highly signficant
AIC(Tot0) #1852.415

# Residual vs fitted values
E0 <- resid(Tot0, type ="pearson")
F0 <- fitted(Tot0)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F0, 
     y = E0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK - becoming conical - not too bad

##### Regrowth autocorrelation model ####
Tot1<-lme(TotalBiomass1~Treatment+fBoma.density+Treatment:fBoma.density,
          random= ~ 1|fBlock/fPlot.ID, na.action=na.pass, method="REML",
          correlation=corAR1(0.2, form=~YrMonthNumber|fBlock/fPlot.ID),data=nsReharvestH1SL)
summary(Tot1)
anova(Tot1) # Highly signficant
AIC(Tot0) #1852.415

# Residual vs fitted values
E1 <- resid(Tot1, type ="pearson")
F1 <- fitted(Tot1)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # Bad = conical

# Issues residual vs fit - conical 

#### Total -Mixed model with season as fixed factor ####
Tot0mod<-lmer(TotalBiomass0~Treatment+fBoma.density+Reharvest.date2+ 
                   Treatment:fBoma.density+Treatment:Reharvest.date2+
                   Reharvest.date2:fBoma.density+
                   Reharvest.date2:fBoma.density:Treatment+
                   (1|fBlock), data=nsReharvest)
summary(Tot0mod)
anova(Tot0mod) # Highly signficant
AIC(Tot0mod) #1852.415

# Residual vs fitted values
E0 <- resid(Tot0mod, type ="pearson")
F0 <- fitted(Tot0mod)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F0, 
     y = E0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK

drop1(Tot0mod, test="Chisq")

# Generate pvalues
Tot0mod1<-update(Tot0mod,~.-Treatment:fBoma.density:Reharvest.date2)
Tot0mod2<-update(Tot0mod1,~.-Treatment:Reharvest.date2)
Tot0mod3<-update(Tot0mod1,~.-Treatment:fBoma.density)
Tot0mod4<-update(Tot0mod1,~.-Reharvest.date2:fBoma.density)
Tot0mod5<-update(Tot0mod2,~.-Treatment:fBoma.density)
Tot0mod6<-update(Tot0mod5,~.-Reharvest.date2:fBoma.density)

Tot0mod7<-update(Tot0mod6,~.-Reharvest.date2)
Tot0mod8<-update(Tot0mod6,~.-fBoma.density)
Tot0mod9<-update(Tot0mod6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Tot0mod,Tot0mod1) #
anova(Tot0mod1,Tot0mod2) # 
anova(Tot0mod1,Tot0mod3) #
anova(Tot0mod1,Tot0mod4) # 
anova(Tot0mod6,Tot0mod7) #
anova(Tot0mod6,Tot0mod8) #
anova(Tot0mod6,Tot0mod9) #

#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#14 2926.5 2976.8 -1449.2   2898.5 5.7604      2    0.05612 . #Treatment:fBoma.density:Reharvest.date2
#12 2928.2 2971.4 -1452.1   2904.2 1.5376      2     0.4636 # Treatment:Reharvest.date2
#12 2928.2 2971.4 -1452.1   2904.2 0.4669      1     0.4944 # Treatment:fBoma.density
#12 2928.2 2971.4 -1452.1   2904.2 12.952      2    0.00154 ** # Reharvest.date2:fBoma.density
#7 2933.1 2958.3 -1459.5   2919.1 70.742      2   4.35e-16 *** #Reharvest.date2
#7 2933.1 2958.3 -1459.5   2919.1 1.0656      1     0.3019 # fBoma.density
#7 2933.1 2958.3 -1459.5   2919.1 39.664      1  3.016e-10 *** #Treatment

aggregate(TotalBiomass0~fBoma.density,nsReharvestH1,mean)
aggregate(TotalBiomass0~fBoma.density,nsReharvestH1,sd)

#     fBoma.density TotalBiomass0
#1          High       97.6125 + 55.5
#2           Low      116.4553 + 79.6

aggregate(TotalBiomass0~fBoma.density+Treatment,nsReharvestH1,mean)
#fBoma.density Treatment TotalBiomass1
#1          High   Control      74.43667
#2           Low   Control      97.53467
#3          High Exclosure     120.78833
#4           Low Exclosure     135.37600
(135.37600-97.53467)/97.53467*100 # 38.79782 % Low
(120.78833-74.43667)/74.43667*100 # 62.26993 % High 

#### Regrowth - Mixed model with season as fixed factor ####
# N.B. Only two seasons for regrowth model ####
Tot1mod<-lmer(TotalBiomass1~Treatment+fBoma.density+Reharvest.date2+ 
                Treatment:fBoma.density+Treatment:Reharvest.date2+
                Reharvest.date2:fBoma.density+
                Reharvest.date2:fBoma.density:Treatment+
                (1|fBlock), data=nsReharvestH1SL)
summary(Tot1mod)
anova(Tot1mod) # Highly signficant
AIC(Tot1mod) #1852.415

# Residual vs fitted values
E1 <- resid(Tot1mod, type ="pearson")
F1 <- fitted(Tot1mod)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK - largest values played - but not second largest

# Generate pvalues
Tot1mod1<-update(Tot1mod,~.-Treatment:fBoma.density:Reharvest.date2)
Tot1mod2<-update(Tot1mod1,~.-Treatment:Reharvest.date2)
Tot1mod3<-update(Tot1mod1,~.-Treatment:fBoma.density)
Tot1mod4<-update(Tot1mod1,~.-Reharvest.date2:fBoma.density)
Tot1mod5<-update(Tot1mod2,~.-Treatment:fBoma.density)
Tot1mod6<-update(Tot1mod5,~.-Reharvest.date2:fBoma.density)

Tot1mod7<-update(Tot1mod6,~.-Reharvest.date2)
Tot1mod8<-update(Tot1mod6,~.-fBoma.density)
Tot1mod9<-update(Tot1mod6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Tot1mod,Tot1mod1) #
anova(Tot1mod,Tot1mod2) # 
anova(Tot1mod1,Tot1mod3) #
anova(Tot1mod1,Tot1mod4) # 
anova(Tot1mod6,Tot1mod7) #
anova(Tot1mod6,Tot1mod8) #
anova(Tot1mod6,Tot1mod9) #

# NoPar  LogLik Df Deviance Pr(>Chi)
# 10 1759.3 1791.2 -869.63   1739.3 2.3555      1     0.1248  # Treatment:fBoma.density:Reharvest.date2
# 10 1759.3 1791.2 -869.63   1739.3 5.4948      2    0.06409 . #Treatment:Reharvest.date2
# 9 1759.6 1788.4 -870.81   1741.6 0.2073      1     0.6489  #Treatment:fBoma.density
# 9 1759.6 1788.4 -870.81   1741.6 6.9232      1   0.008508 ** #Reharvest.date2:fBoma.density
# 6 1763.8 1782.9 -875.88   1751.8 6.6198      1    0.01008 * #Reharvest.date2
# 6 1763.8 1782.9 -875.88   1751.8 2.6963      1     0.1006 # fBoma.density
# 6 1763.8 1782.9 -875.88   1751.8 21.443      1  3.645e-06 *** #Treatment

aggregate(TotalBiomass1~fBoma.density,nsReharvestH1,mean)
aggregate(TotalBiomass1~fBoma.density,nsReharvestH1,sd)

#fBoma.density TotalBiomass1
#1          High      66.29417 + 57.92966
#2           Low      95.15533 + 78.77888

aggregate(TotalBiomass1~fBoma.density+Treatment,nsReharvestH1,mean)
#fBoma.density Treatment TotalBiomass1
#1          High   Control      51.62333
#2           Low   Control      81.08933
#3          High Exclosure      80.96500
#4           Low Exclosure     109.22133
(109.22133-81.08933)/81.08933*100 # 34.6926 % Low
(80.96500-51.62333)/51.62333*100 # 56.838 % High 

# Total - regrowth difference
nsReharvestH1SL$TotBioDiff<-nsReharvestH1SL$TotalBiomass0-nsReharvestH1SL$TotalBiomass1
TotDmod<-lmer(TotBioDiff~Treatment+fBoma.density+Reharvest.date2+ 
                Treatment:fBoma.density+Treatment:Reharvest.date2+
                Reharvest.date2:fBoma.density+
                Reharvest.date2:fBoma.density:Treatment+
                (1|fBlock), data=nsReharvestH1SL)
summary(TotDmod)
anova(TotDmod) # Highly signficant
AIC(TotDmod) #1852.415

# Residual vs fitted values
Ed <- resid(TotDmod, type ="pearson")
Fd <- fitted(TotDmod)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = Fd, 
     y = Ed,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK - largest values played - but not second largest

# Generate pvalues
TotDmod1<-update(TotDmod,~.-Treatment:fBoma.density:Reharvest.date2)
TotDmod2<-update(TotDmod1,~.-Treatment:Reharvest.date2)
TotDmod3<-update(TotDmod1,~.-Treatment:fBoma.density)
TotDmod4<-update(TotDmod1,~.-Reharvest.date2:fBoma.density)
TotDmod5<-update(TotDmod2,~.-Treatment:fBoma.density)
TotDmod6<-update(TotDmod5,~.-Reharvest.date2:fBoma.density)

TotDmod7<-update(TotDmod6,~.-Reharvest.date2)
TotDmod8<-update(TotDmod6,~.-fBoma.density)
TotDmod9<-update(TotDmod6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(TotDmod,TotDmod1) #
anova(TotDmod,TotDmod2) # 
anova(TotDmod1,TotDmod3) #
anova(TotDmod1,TotDmod4) # 
anova(TotDmod6,TotDmod7) #
anova(TotDmod6,TotDmod8) #
anova(TotDmod6,TotDmod9) #

# NoPar  LogLik Df Deviance Pr(>Chi)
# TotDmod  10 1947.0 1978.9 -963.49   1927.0 10.426      1   0.001243 **  # Treatment:fBoma.density:Reharvest.date2
# TotDmod  10 1947.0 1978.9 -963.49   1927.0 15.212      2  0.0004974 *** #Treatment:Reharvest.date2
# TotDmod1  9 1955.4 1984.1 -968.70   1937.4 0.4891      1     0.4843  #Treatment:fBoma.density
# TotDmod1  9 1955.4 1984.1 -968.7   1937.4 0.7919      1     0.3735 #Reharvest.date2:fBoma.density
# TotDmod6  6 1955.4 1974.6 -971.72   1943.4 3.7837      1    0.05175 . #Reharvest.date2
# TotDmod6  6 1955.4 1974.6 -971.72   1943.4 1.3686      1     0.2421 # fBoma.density
# TotDmod6  6 1955.4 1974.6 -971.72   1943.4 5.9279      1     0.0149 * #Treatment

# Within factors - Season x tukuls density x treatment
ref.grid.TotDmodFINAL <- ref_grid(TotDmod) #at = list(Region = c("Dry", "Wet"))) #Want to remove Intermadiate in the grid as it can't be compared to the other landuses, except wild.
ref.grid.TotDmodFINAL
#Check threeway first:
emmip(ref.grid.TotDmodFINAL, Reharvest.date2~Treatment|fBoma.density, type="response")#
emmeans.TotDmodFINAL <- emmeans(ref.grid.TotDmodFINAL, pairwise~Treatment*fBoma.density*Reharvest.date2|fBoma.density,type="response") #
emmeans.TotDmodFINAL$contrasts
emmeans.TotDmodFINAL$emmeans
emmeans.TotDmodFINAL.pairs <- pairs(emmeans.TotDmodFINAL,simple = "each", combine =TRUE)
plot(emmeans.TotDmodFINAL, comparisons = FALSE)


#Check dates...
#nsReharvestH1SL<-nsReharvestH1[nsReharvestH1$Reharvest.date2!="Short I",]

#### Temporal autocorrelation model with gamma ####
#Tot1TMB<-glmmTMB(TotalBiomass1~Treatment+fBoma.density+TotalBiomass0+ 
#                   Treatment:fBoma.density+Treatment:TotalBiomass0+
#                   TotalBiomass0:fBoma.density+
#                   TotalBiomass0:fBoma.density:Treatment+
#          (1|fBlock/fPlot.ID) + 
#          ar1(YrMonth-1|fBlock/fPlot.ID),
#        data=nsReharvestH1SL,family=Gamma)
#          #list(family="Gamma",link="inverse"))
#### Testing INLA ####
#library(INLA)
#library(brinla)
#library(nlme)
#library(mgcv)

# Formula
#fBioRainMod0<- TotalBiomass1~Treatment+fBoma.density+TotalBiomass0.c+
#  Treatment:fBoma.density+Treatment:TotalBiomass0.c+
 # TotalBiomass0.c:fBoma.density+
#  TotalBiomass0.c:fBoma.density:Treatment+ f(YrMonth, model="ar1", replicate=Plot.ID)+
 # f(fBlock, model = "iid")
#nsReharvestH1$TotalBiomass1[nsReharvestH1$TotalBiomass1==0]<-0.01
#BioRainMod <- inla(fBioRainMod0, family = "gamma",
#                   control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE),
#                   data = nsReharvestH1)
#BioRainMod$waic$waic #2714.748
#pvalue histogram....
#par(mfrow = c(1, 1), mar = c(4, 3, 3, 2))
#pval<-rep(NA, nrow=(nsReharvestH1))
#for(i in 1:nrow(nsReharvestH1)){
#  pval[i]<-inla.pmarginal(q=nsReharvestH1$TotalBiomass1[i],
#                          marginal=BioRainMod$marginals.fitted.values[[i]])
#}
#hist(pval) # Bimodial - not good fit...
#bri.fixed.plot(BioRainMod)

#### Rereharvested model ####
#nsReharvestH1$TotalBiomass1[nsReharvestH1$TotalBiomass1==0]<-0.01
#nsReharvestH1b<-droplevels(nsReharvestH1[nsReharvestH1$Reharvest.date2!="Short I",])
#nsReharvestH1b$Reharvest.date2

#### GLMM with Gamma ####
#Tot1gam<-glmmadmb(TotalBiomass1~Treatment+fBoma.density+TotalBiomass0+ 
#                    Treatment:fBoma.density+Treatment:TotalBiomass0+
#                    TotalBiomass0:fBoma.density+
#                    TotalBiomass0:fBoma.density:Treatment+
#                    (1|fBlock)+(1|Reharvest.date2),family="Gamma",
#                  # admb.opts=admbControl(shess=FALSE,noinit=FALSE),
#                  data=nsReharvestH1SL) #+(1|rain.mm)
#summary(Tot1gam)
#AIC(Tot1gam) #2671.18
#drop1(Tot1gam,test="Chisq")

# Residual vs fitted values
#E1 <- resid(Tot1gam, type ="pearson")
#F1 <- fitted(Tot1gam)
#par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
#plot(x = F1, 
#     y = E1,
#     xlab = "Fitted values",
#     ylab = "Residuals")
#abline(v = 0, lwd = 2, col = 2)
#abline(h = 0, lty = 2, col = 1) # Good spread - not skewed

# Check deviance
#sum(resid(Tot1gam))
#-2*logLik(Tot1gam) # 'log Lik.'
#(2647.18--0.5372177)/2647.18 #  1.000203 # Excellent

# Interaction / boma density and season
with(nsReharvestH1, {interaction.plot(Reharvest.date2,fBoma.density,TotalBiomass1,
                                     xlab = "Season",
                                     ylab = "Biomass Reharvest",
                                     fun=mean)})

aggregate(TotalBiomass1~fBoma.density+Reharvest.date2,nsReharvestH1,mean)

aggregate(TotalBiomass1~Boma.density+Reharvest.date+Harvest,nsReharvest,mean)
nsReharvestH1$Harvest

ggplot(nsReharvestH1SL,aes(y=TotalBiomass1,x=Reharvest.date2, col=fBoma.density))+facet_wrap(~fBoma.density+Treatment)+geom_point()
table(nsReharvestH1$fBoma.density,nsReharvestH1$fBoma.density)

#### Non-reharvested model ####
#nsReharvestH1$TotalBiomass0[nsReharvestH1$TotalBiomass0==0]<-0.01
#Tot0gam<-glmmadmb(TotalBiomass0~Treatment+fBoma.density+Reharvest.date2+
#            Reharvest.date2:fBoma.density+ Treatment:fBoma.density+
#               Treatment:Reharvest.date2+
#             Treatment:fBoma.density:Reharvest.date2+
#                 (1|fBlock),family="Gamma",
#              # admb.opts=admbControl(shess=FALSE,noinit=FALSE),
#               data=nsReharvestH1) #+(1|rain.mm)
#correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1)
#summary(Tot0gam)
#AIC(Tot0gam)
#drop1(Tot0gam,test="Chisq")

# Residual vs fitted values
#E0 <- resid(Tot0gam, type ="pearson")
#F0 <- fitted(Tot0gam)
#par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
#plot(x = F0, 
#     y = E0,
#     xlab = "Fitted values",
#     ylab = "Residuals")
#abline(v = 0, lwd = 2, col = 2)
#abline(h = 0, lty = 2, col = 1) # Looks good

# Check deviance
#sum(resid(Tot0gam))
#-2*logLik(Tot0gam) # 'log Lik.'
#(2899.06-0.007022121)/2899.06 #  0.9999976 # Excellent

# Interaction / boma density and season
with(nsReharvestH1, {interaction.plot(Reharvest.date2,fBoma.density,TotalBiomass0,
                                      xlab = "Season",
                                      ylab = "Biomass Harvest",
                                      fun=mean)})

nsReharvestH1Op<-nsReharvestH1[nsReharvestH1$Treatment=="Control",]
nsReharvestH1Ex<-nsReharvestH1[nsReharvestH1$Treatment=="Exclosure",]
with(nsReharvestH1Op, {interaction.plot(Reharvest.date2,fBoma.density,TotalBiomass0,
                                      xlab = "Season",
                                      ylab = "Biomass Harvest",
                                      fun=mean)})
with(nsReharvestH1Ex, {interaction.plot(Reharvest.date2,fBoma.density,TotalBiomass0,
                                        xlab = "Season",
                                        ylab = "Biomass Harvest",
                                        fun=mean)})

#Zeb<-ggplot(nsReharvestH1,aes(y=TotalBiomass1, x=Burchells_ZebraMetBio))
#Zeb<-Zeb+geom_point(data=nsReharvestbH1_2,aes(y=TotalBiomass0, x=Burchells_ZebraMetBio),size=3, colour="light grey", alpha=.5)
#Zeb<-Zeb+geom_point(size =3,col="dark grey")
#Zeb<-Zeb+facet_wrap(~Treatment)
#Zeb<-Zeb+theme_classic()
#Zeb

# Standing biomass - i.e. no regrowth
#Tot1b<-lmer(TotalBiomass0~Treatment+
#             +Burchells_ZebraMetBio+#CattleMetBio+                 
#             +Grants_GazelleMetBio+Swaynes_HartebeestMetBio+ #Greater_KuduMetBio
#             +Burchells_ZebraMetBio:Treatment+#CattleMetBio:Treatment+                 
             # +Grants_GazelleMetBio:Treatment+Greater_KuduMetBio:Treatment+
             #Swaynes_HartebeestMetBio:Treatment+ 
#             (1|fBlock)+(1|rain.mm),data=nsReharvestbH1_2)
#correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1)

#### TEST DISTANCE TO SETTLEMENTS - NOT DENSITY #####

names(nsReharvestH1) #min_distExclosures
plot(nsReharvestH1$min_distExclosures)

#### Total -Mixed model with season as fixed factor - NEAREST TUKULS ####
Tot0modDIS<-lmer(TotalBiomass0~Treatment+min_distExclosures+Reharvest.date2+ 
                Treatment:min_distExclosures+Treatment:Reharvest.date2+
                Reharvest.date2:min_distExclosures+
                Reharvest.date2:min_distExclosures:Treatment+
                (1|fBlock), data=nsReharvestH1)
summary(Tot0modDIS)
anova(Tot0modDIS) # Highly signficant
AIC(Tot0modDIS) #1852.415

# Residual vs fitted values
E0 <- resid(Tot0modDIS, type ="pearson")
F0 <- fitted(Tot0modDIS)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F0, 
     y = E0,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK

drop1(Tot0modDIS, test="Chisq")

# Generate pvalues
Tot0modDIS1<-update(Tot0modDIS,~.-Treatment:min_distExclosures:Reharvest.date2)
Tot0modDIS2<-update(Tot0modDIS1,~.-Treatment:Reharvest.date2)
Tot0modDIS3<-update(Tot0modDIS1,~.-Treatment:min_distExclosures)
Tot0modDIS4<-update(Tot0modDIS1,~.-Reharvest.date2:min_distExclosures)
Tot0modDIS5<-update(Tot0modDIS2,~.-Treatment:min_distExclosures)
Tot0modDIS6<-update(Tot0modDIS5,~.-Reharvest.date2:min_distExclosures)

Tot0modDIS7<-update(Tot0modDIS6,~.-Reharvest.date2)
Tot0modDIS8<-update(Tot0modDIS6,~.-min_distExclosures)
Tot0modDIS9<-update(Tot0modDIS6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Tot0modDIS,Tot0modDIS1) #
anova(Tot0modDIS1,Tot0modDIS2) # 
anova(Tot0modDIS1,Tot0modDIS3) #
anova(Tot0modDIS1,Tot0modDIS4) # 
anova(Tot0modDIS6,Tot0modDIS7) #
anova(Tot0modDIS6,Tot0modDIS8) #
anova(Tot0modDIS6,Tot0modDIS9) #

#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#Tot0modDIS  14 2928.4 2978.8 -1450.2   2900.4 4.4809      2     0.1064 #Treatment:min_distExclosures:Reharvest.date2
#Tot0modDIS1 12 2928.9 2972.1 -1452.4   2904.9 1.5342      2     0.4644 #Treatment:Reharvest.date2
#Tot0modDIS1 12 2928.9 2972.1 -1452.4   2904.9 1.8284      1     0.1763 #Treatment:min_distExclosures
#Tot0modDIS1 12 2928.9 2972.1 -1452.4   2904.9 2.1245      2     0.3457 #Reharvest.date2:min_distExclosures
#Tot0modDIS6  7 2924.3 2949.5 -1455.2   2910.3 72.74      2  < 2.2e-16 *** #Reharvest.date2
#Tot0modDIS6  7 2924.3 2949.5 -1455.2   2910.3 9.8183      1   0.001728 ** #min_distExclosures
#Tot0modDIS6  7 2924.3 2949.5 -1455.2   2910.3 40.874      1  1.624e-10 *** #Treatment

plot(TotalBiomass0~min_distExclosures,nsReharvestH1)
abline(lm(TotalBiomass0~min_distExclosures,nsReharvestH1), lwd=2, col="red")

#### Regrowth - Mixed model with season as fixed factor - min_distExclosures ####
# N.B. Only two seasons for regrowth model ####
Tot1modDIS<-lmer(TotalBiomass1~Treatment+min_distExclosures+Reharvest.date2+ 
                Treatment:min_distExclosures+Treatment:Reharvest.date2+
                Reharvest.date2:min_distExclosures+
                Reharvest.date2:min_distExclosures:Treatment+
                (1|fBlock), data=nsReharvestH1SL)
summary(Tot1modDIS)
anova(Tot1modDIS) # Highly signficant
AIC(Tot1modDIS) #1852.415

# Residual vs fitted values
E1 <- resid(Tot1modDIS, type ="pearson")
F1 <- fitted(Tot1modDIS)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK - largest values played - but not second largest

# Generate pvalues
Tot1modDIS1<-update(Tot1modDIS,~.-Treatment:min_distExclosures:Reharvest.date2)
Tot1modDIS2<-update(Tot1modDIS1,~.-Treatment:Reharvest.date2)
Tot1modDIS3<-update(Tot1modDIS1,~.-Treatment:min_distExclosures)
Tot1modDIS4<-update(Tot1modDIS1,~.-Reharvest.date2:min_distExclosures)
Tot1modDIS5<-update(Tot1modDIS2,~.-Treatment:min_distExclosures)
Tot1modDIS6<-update(Tot1modDIS5,~.-Reharvest.date2:min_distExclosures)

Tot1modDIS7<-update(Tot1modDIS6,~.-Reharvest.date2)
Tot1modDIS8<-update(Tot1modDIS6,~.-min_distExclosures)
Tot1modDIS9<-update(Tot1modDIS6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Tot1modDIS,Tot1modDIS1) #
anova(Tot1modDIS,Tot1modDIS2) # 
anova(Tot1modDIS1,Tot1modDIS3) #
anova(Tot1modDIS1,Tot1modDIS4) # 
anova(Tot1modDIS6,Tot1modDIS7) #
anova(Tot1modDIS6,Tot1modDIS8) #
anova(Tot1modDIS6,Tot1modDIS9) #

# NoPar  LogLik Df Deviance Pr(>Chi)
# Tot1modDIS  10 1737.4 1769.3 -858.71   1717.4 8.7092      1   0.003166 ** # Treatment:min_distExclosures:Reharvest.date2
# Tot1modDIS  10 1737.4 1769.3 -858.71   1717.4 12.05      2   0.002418 ** #Treatment:Reharvest.date2
# Tot1modDIS1  9 1744.1 1772.9 -863.06   1726.1 0.3528      1     0.5525  #Treatment:min_distExclosures
# Tot1modDIS1  9 1744.1 1772.9 -863.06   1726.1 14.265      1  0.0001588 *** #Reharvest.date2:min_distExclosures
# Tot1modDIS6  6 1755.8 1775.0 -871.90   1743.8 6.7457      1   0.009397 ** #Reharvest.date2
# Tot1modDIS6  6 1755.8 1775.0 -871.90   1743.8 10.664      1   0.001093 ** # min_distExclosures
# Tot1modDIS6  6 1755.8 1775.0 -871.90   1743.8 21.831      1  2.978e-06 *** #Treatment

xyplot(TotalBiomass1~min_distExclosures|Treatment*R,nsReharvestH1SL)

#### Total - regrowth difference - Nearest Tukuls####
nsReharvestH1SL$TotBioDiff<-nsReharvestH1SL$TotalBiomass0-nsReharvestH1SL$TotalBiomass1
TotDmodDIS<-lmer(TotBioDiff~Treatment+min_distExclosures+Reharvest.date2+ 
                Treatment:min_distExclosures+Treatment:Reharvest.date2+
                Reharvest.date2:min_distExclosures+
                Reharvest.date2:min_distExclosures:Treatment+
                (1|fBlock), data=nsReharvestH1SL)
summary(TotDmodDIS)
anova(TotDmodDIS) # Highly signficant
AIC(TotDmodDIS) #1852.415

# Residual vs fitted values
Ed <- resid(TotDmodDIS, type ="pearson")
Fd <- fitted(TotDmodDIS)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = Fd, 
     y = Ed,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1) # OK - largest values played - but not second largest

# Generate pvalues
TotDmodDIS1<-update(TotDmodDIS,~.-Treatment:min_distExclosures:Reharvest.date2)
TotDmodDIS2<-update(TotDmodDIS1,~.-Treatment:Reharvest.date2)
TotDmodDIS3<-update(TotDmodDIS1,~.-Treatment:min_distExclosures)
TotDmodDIS4<-update(TotDmodDIS1,~.-Reharvest.date2:min_distExclosures)
TotDmodDIS5<-update(TotDmodDIS2,~.-Treatment:min_distExclosures)
TotDmodDIS6<-update(TotDmodDIS5,~.-Reharvest.date2:min_distExclosures)

TotDmodDIS7<-update(TotDmodDIS6,~.-Reharvest.date2)
TotDmodDIS8<-update(TotDmodDIS6,~.-min_distExclosures)
TotDmodDIS9<-update(TotDmodDIS6,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(TotDmodDIS,TotDmodDIS1) #
anova(TotDmodDIS,TotDmodDIS2) # 
anova(TotDmodDIS1,TotDmodDIS3) #
anova(TotDmodDIS1,TotDmodDIS4) # 
anova(TotDmodDIS6,TotDmodDIS7) #
anova(TotDmodDIS6,TotDmodDIS8) #
anova(TotDmodDIS6,TotDmodDIS9) #

# NoPar  LogLik Df Deviance Pr(>Chi)
# TotDmodDIS  10 1949.5 1981.4 -964.75   1929.5 9.5003      1   0.002054 ** #Treatment:min_distExclosures:Reharvest.date2
# TotDmodDIS  10 1949.5 1981.4 -964.75   1929.5 14.293      2  0.0007877 *** #Treatment:Reharvest.date2
# TotDmodDIS1  9 1957.0 1985.7 -969.50   1939.0 0.3145      1     0.5749 #Treatment:min_distExclosures
# TotDmodDIS1  9 1957.0 1985.7 -969.50   1939.0 0.5333      1     0.4652 #Reharvest.date2:min_distExclosures
# TotDmodDIS6  6 1956.6 1975.8 -972.31   1944.6 3.7985      1     0.0513 . #Reharvest.date2
# TotDmodDIS6  6 1956.6 1975.8 -972.31   1944.6 0.1958      1     0.6581  #min_distExclosures
# TotDmodDIS6  6 1956.6 1975.8 -972.31   1944.6 5.9513      1    0.01471 * #Treatment

#####################################################################################
#### GRASS BIOMASS Regrowth only double harvest ####
#####################################################################################

#nsReharvestHmeanSL<-nsReharvestHmean[nsReharvestHmean$Reharvest.date2!="Short I",]

Grass1<-lmer(GrassNetReharvestBiomass1~Livestock.density+Treatment+Season+
             Season:Treatment+Livestock.density:Season+
             Livestock.density:Treatment+ # Simplified via LRT
             Treatment:Livestock.density:Season+
             (1|fBlock), REML=F,data=nsReharvestH1)
summary(Grass1)
anova(Grass1)
AIC(Grass1) #1745.778

# Simplfy model using LRT
drop1(Grass1,test="Chisq")

# Update and remove factors # issues with interactions
Grass2<-lmer(GrassNetReharvestBiomass1~Livestock.density+Treatment+Season+(1|fBlock), REML=F,data=nsReharvestH1)
Grass1a <- update(Grass1, .~. -Treatment:Livestock.density:Season)
Grass1a2 <- update(Grass1a, .~. -Livestock.density:Treatment)
Grass1a3 <- update(Grass1a, .~. -Livestock.density:Season)
Grass1a4 <- update(Grass1a, .~. -Season:Treatment)
Grass2a <- update(Grass2, .~. -Season)
Grass2b <- update(Grass2, .~. -Treatment)
Grass2c <- update(Grass2, .~. -Livestock.density)

anova(Grass1,Grass1a)
anova(Grass1a,Grass1a2)
anova(Grass1a,Grass1a3)
anova(Grass1a,Grass1a4)
anova(Grass2,Grass2a)
anova(Grass2,Grass2b)
anova(Grass2,Grass2c)

#      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#Grass1  14 1745.8 1790.5 -858.89   1717.8 7.1326      2    0.02826 * #Treatment:Livestock.density:Harvest.date
#Grass1a  12 1748.9 1787.2 -862.46   1724.9 2.8283      2     0.2431 # Livestock.density:Treatment
#Grass1a  12 1748.9 1787.2 -862.46   1724.9 19.745      2  5.157e-05 *** #Livestock.density:Harvest.date
#Grass1a  12 1748.9 1787.2 -862.46   1724.9 1.6315      1     0.2015 # Harvest.date:Treatment
#Grass2   7 1762.6 1785.0 -874.31   1748.6 3.8582      1     0.0495 * # Harvest.date
#Grass2   7 1762.6 1785 -874.31   1748.6 21.256      1  4.019e-06 *** # Treatment
#Grass2   7 1762.6 1785.0 -874.31   1748.6 9.2415      2   0.009846 ** # Livestock.density

####################################################################################
#### Herbs and Climbers BIOMASS Regrowth only double harvest####
#####################################################################################
nsReharvestH1$HerbClimber0<-nsReharvestH1$HerbNetReharvestBiomass0+nsReharvestH1$ClimberNetReharvestBiomass0
nsReharvestH1SL$HerbClimber0<-nsReharvestH1SL$HerbNetReharvestBiomass0+nsReharvestH1SL$ClimberNetReharvestBiomass0
nsReharvestH1SL$HerbClimber1<-nsReharvestH1SL$HerbNetReharvestBiomass1+nsReharvestH1SL$ClimberNetReharvestBiomass1

names(nsReharvestH1)

#### Herb total biomass #####
Herb0<-lmer(HerbClimber0~Treatment+fBoma.density+Reharvest.date2+ 
              Treatment:fBoma.density+Treatment:Reharvest.date2+
              Reharvest.date2:fBoma.density+
              Reharvest.date2:fBoma.density:Treatment+
               (1|fBlock), REML=F,data=nsReharvestH1)
summary(Herb0)
anova(Herb0)
AIC(Herb0) #1357.897

#Inspect chosen model for homogeneity:
E1 <- resid(Herb0, type ="pearson")
F1 <- fitted(Herb0)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1)

#drop1(Herb0, test="Chisq")

# Generate pvalues
Herb01<-update(Herb0,~.-Treatment:fBoma.density:Reharvest.date2)
Herb02<-update(Herb01,~.-Treatment:Reharvest.date2)
Herb03<-update(Herb01,~.-Treatment:fBoma.density)
Herb04<-update(Herb01,~.-Reharvest.date2:fBoma.density)
Herb05<-update(Herb02,~.-Treatment:fBoma.density)
Herb06<-update(Herb05,~.-Reharvest.date2:fBoma.density)

Herb07<-update(Herb06,~.-Reharvest.date2)
Herb08<-update(Herb06,~.-fBoma.density)
Herb09<-update(Herb06,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Herb0,Herb01) #
anova(Herb01,Herb02) # 
anova(Herb01,Herb03) #
anova(Herb01,Herb04) # 
anova(Herb06,Herb07) #
anova(Herb06,Herb08) #
anova(Herb06,Herb09) #

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#Herb0  14 2166.1 2216.5 -1069.0   2138.1 2.1975      2     0.3333 #Treatment:fBoma.density:Reharvest.date2
#Herb01 12 2164.3 2207.5 -1070.1   2140.3 1.325      2     0.5156 # Treatment:Reharvest.date2
#Herb01 12 2164.3 2207.5 -1070.1   2140.3 9.0708      1   0.002597 ** # Treatment:fBoma.density
#Herb01 12 2164.3 2207.5 -1070.1   2140.3 2.403      2     0.3007 # Reharvest.date2:fBoma.density
#Herb06  7 2166.9 2192.1 -1076.5   2152.9 7.5914      2    0.02247 * # Reharvest.date2
#Herb06  7 2166.9 2192.1 -1076.5   2152.9 2.7375      1    0.09801 . #fBoma.density
#Herb06  7 2166.9 2192.1 -1076.5   2152.9 1.0959      1     0.2952 # Treatment

#### Herb regrowth ####
Herb1<-lmer(HerbClimber1~Treatment+fBoma.density+Reharvest.date2+ 
                Treatment:fBoma.density+Treatment:Reharvest.date2+
                Reharvest.date2:fBoma.density+
                Reharvest.date2:fBoma.density:Treatment+
                (1|fBlock), data=nsReharvestH1SL)
summary(Herb1)
anova(Herb1)
AIC(Herb1) #1357.897

#Inspect chosen model for homogeneity:
E1 <- resid(Herb1, type ="pearson")
F1 <- fitted(Herb1)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1)

# Simplfy model using LRT
#drop1(Herb1,test="Chisq") # Nothing singificant....

# Generate pvalues
Herb11<-update(Herb1,~.-Treatment:fBoma.density:Reharvest.date2)
Herb12<-update(Herb11,~.-Treatment:Reharvest.date2)
Herb13<-update(Herb11,~.-Treatment:fBoma.density)
Herb14<-update(Herb11,~.-Reharvest.date2:fBoma.density)
Herb15<-update(Herb12,~.-Treatment:fBoma.density)
Herb16<-update(Herb15,~.-Reharvest.date2:fBoma.density)

Herb17<-update(Herb16,~.-Reharvest.date2)
Herb18<-update(Herb16,~.-fBoma.density)
Herb19<-update(Herb16,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Herb1,Herb11) #
anova(Herb11,Herb12) # 
anova(Herb11,Herb13) #
anova(Herb11,Herb14) # 
anova(Herb16,Herb17) #
anova(Herb16,Herb18) #
anova(Herb16,Herb19) #

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#Herb1  10 1353.3 1385.3 -666.67   1333.3 0.3991      1     0.5276 #Treatment:fBoma.density:Reharvest.date2
#Herb11  9 1351.7 1380.5 -666.87   1333.7 0.8398      1     0.3594 #Treatment:Reharvest.date2
#Herb11  9 1351.7 1380.5 -666.87   1333.7 1.3451      1     0.2461 #Treatment:fBoma.density
#Herb11  9 1351.7 1380.5 -666.87   1333.7 3e-04      1     0.9854 # Reharvest.date2:fBoma.density
#Herb16  6 1347.9 1367.1 -667.96   1335.9 4.5544      1    0.03283 *# Reharvest.date2
#Herb16  6 1347.9 1367.1 -667.96   1335.9 2.6103      1     0.1062 #fBoma.density
#Herb16  6 1347.9 1367.1 -667.96   1335.9 0.048      1     0.8266# Treatment

#### Herb regrowth ####
nsReharvestH1SL$HerbDiff<-nsReharvestH1SL$HerbClimber0-nsReharvestH1SL$HerbClimber1

HerbDifmod<-lmer(HerbDiff~Treatment+fBoma.density+Reharvest.date2+ 
              Treatment:fBoma.density+Treatment:Reharvest.date2+
              Reharvest.date2:fBoma.density+
              Reharvest.date2:fBoma.density:Treatment+
              (1|fBlock), data=nsReharvestH1SL)

summary(HerbDifmod)
anova(HerbDifmod)
AIC(HerbDifmod) #1357.897

drop1(HerbDifmod, test="Chisq")


########################################################################################
#### Woody BIOMASS Regrowth only double harvest####
Woody0<-lmer(DwarfShrubNetReharvestBiomass0~Treatment+fBoma.density+Reharvest.date2+ 
               Treatment:fBoma.density+Treatment:Reharvest.date2+
               Reharvest.date2:fBoma.density+
               Reharvest.date2:fBoma.density:Treatment+
               (1|fBlock), REML=F,data=nsReharvestH1)
summary(Woody0)
anova(Woody0)
AIC(Woody0) #1259.642

#Inspect chosen model for homogeneity:
E1 <- resid(Woody0, type ="pearson")
F1 <- fitted(Woody0)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 1)

# Generate pvalues
Woody01<-update(Woody0,~.-Treatment:fBoma.density:Reharvest.date2)
Woody02<-update(Woody01,~.-Treatment:Reharvest.date2)
Woody03<-update(Woody01,~.-Treatment:fBoma.density)
Woody04<-update(Woody01,~.-Reharvest.date2:fBoma.density)
Woody05<-update(Woody02,~.-Treatment:fBoma.density)
Woody06<-update(Woody05,~.-Reharvest.date2:fBoma.density)

Woody07<-update(Woody06,~.-Reharvest.date2)
Woody08<-update(Woody06,~.-fBoma.density)
Woody09<-update(Woody06,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Woody0,Woody01) #
anova(Woody01,Woody02) # 
anova(Woody01,Woody03) #
anova(Woody01,Woody04) # 
anova(Woody06,Woody07) #
anova(Woody06,Woody08) #
anova(Woody06,Woody09) #

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
#Woody0  14 2136.4 2186.8 -1054.2   2108.4 3.0666      2     0.2158 #Treatment:fBoma.density:Reharvest.date2
#Woody01 12 2135.4 2178.6 -1055.7   2111.4 6.4849      2    0.03907 * #Treatment:Reharvest.date2
#Woody01 12 2135.4 2178.6 -1055.7   2111.4 0.7046      1     0.4013 #Treatment:fBoma.density
#Woody01 12 2135.4 2178.6 -1055.7   2111.4 2.2002      2     0.3328 # Reharvest.date2:fBoma.density
#Woody06  7 2134.8 2159.9 -1060.4   2120.8 0.6339      2     0.7284 # Reharvest.date2
#Woody06  7 2134.8 2159.9 -1060.4   2120.8 0.1687      1     0.6812 #fBoma.density
#Woody06  7 2134.8 2159.9 -1060.4   2120.8 0.0056      1     0.9404 # Treatment

Woody1<-lmer(DwarfShrubNetReharvestBiomass1~Treatment+fBoma.density+Reharvest.date2+ 
               Treatment:fBoma.density+Treatment:Reharvest.date2+
               Reharvest.date2:fBoma.density+
               Reharvest.date2:fBoma.density:Treatment+
               (1|fBlock), REML=F,data=nsReharvestH1SL)
summary(Woody1)
anova(Woody1)
AIC(Woody1) #1259.642

# Simplfy model using LRT
drop1(Woody1,test="Chisq")

# Generate pvalues
Woody11<-update(Woody1,~.-Treatment:fBoma.density:Reharvest.date2)
Woody12<-update(Woody11,~.-Treatment:Reharvest.date2)
Woody13<-update(Woody11,~.-Treatment:fBoma.density)
Woody14<-update(Woody11,~.-Reharvest.date2:fBoma.density)
Woody15<-update(Woody12,~.-Treatment:fBoma.density)
Woody16<-update(Woody15,~.-Reharvest.date2:fBoma.density)

Woody17<-update(Woody16,~.-Reharvest.date2)
Woody18<-update(Woody16,~.-fBoma.density)
Woody19<-update(Woody16,~.-Treatment)

# Pvalues
#NoPar  LogLik Df Deviance Pr(>Chi)
anova(Woody1,Woody11) #
anova(Woody11,Woody12) # 
anova(Woody11,Woody13) #
anova(Woody11,Woody14) # 
anova(Woody16,Woody17) #
anova(Woody16,Woody18) #
anova(Woody16,Woody19) #


# Medium - woody biomass increases in exclosures
# Low and High livestock - exclosures woody biomass lower in exclosures

###########################
# Explore different grazers
Tot1<-lmer(Tot.Per.diff~Treatment+Burchells_ZebraMetBio,
           #Burchells_ZebraMetBio:Reharvest.date,#CattleMetBio+Grants_GazelleMetBio+Treatment:Burchells_ZebraMetBio
        random= ~ 1|fBlock,na.action=na.pass, method="ML",data=nsReharvestH1)
        #correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1)
summary(Tot1)
anova(Tot1)
AIC(Tot1) #2622.329

plot(Tot1) #OK
plot(ACF(Tot1),alpha=0.05) 
drop1(Tot1,test="Chisq") # Burchells_ZebraMetBio  1 1750.1 8.1234  0.00437 **

# 1. Get the betas and gammas
#beta.mcmc <- outB$sims.list$beta 
#dim(beta.mcmc) #3000    6
Tot2<-lme(TotalBiomass1~Burchells_ZebraMetBio:Reharvest.date,
          random= ~ 1|fBlock,na.action=na.pass, method="ML",data=nsReharvestH1)
Betas5 <- fixef(Tot2)
Covbetas5 <- vcov(Tot2)
Betas5

#2. Define a grid of covariate values
#   without extrapolation
MyData5 <- expand.grid(#Treatment=levels(nsReharvestH1$Treatment),
                       Burchells_ZebraMetBio=seq(min(nsReharvestH1$Burchells_ZebraMetBio),max(nsReharvestH1$Burchells_ZebraMetBio), length=25),
                       Reharvest.date=c(min(as.Date(nsReharvestH1$Reharvest.date)),max(as.Date(nsReharvestH1$Reharvest.date))))

head(MyData5)

#Convert the covariate values into an X matrix
Xp <- model.matrix(~Burchells_ZebraMetBio:Reharvest.date, data = MyData5)
dim(Xp) #50 3

#C. Calculate the predicted MCMC values
MyData5$Pred<- Xp %*% as.matrix(Betas5)
dim(MyData5$Pred) # 1250    1

min(MyData5$Pred) # -1.112951

# Upper and lower error
MyData5$SE <- sqrt(  diag(Xp %*% vcov(Tot2) %*% t(Xp))  )

#And using the Pred and SE values, we can calculate
#a 95% confidence interval
MyData5$SeUp <- MyData5$Pred + 1.96 * MyData5$SE
MyData5$SeLo <- MyData5$Pred - 1.96 * MyData5$SE
names(MyData5)
colnames(MyData5)[3]<-"TotalBiomass1"

# Averages for difference herbivores
aggregate(TotalBiomass1~ Burchells_ZebraMetBio + CattleMetBio + Grants_GazelleMetBio + Greater_KuduMetBio+Treatment,nsReharvestH1,mean)
ZebraPer<-aggregate(TotalBiomass1~ Burchells_ZebraMetBio+Reharvest.date+Treatment,nsReharvestH1,mean)
CattlePer<-aggregate(TotalBiomass1~ CattleMetBio+Reharvest.date+Treatment,nsReharvestH1,mean)
GrantPer<-aggregate(TotalBiomass1~ Grants_GazelleMetBio+Reharvest.date+Treatment,nsReharvestH1,mean)
KuduPer<-aggregate(TotalBiomass1~ Greater_KuduMetBio+Reharvest.date+Treatment,nsReharvestH1,mean)

ZebraPerS<-aggregate(TotalBiomass1~ Burchells_ZebraMetBio+Reharvest.date+Treatment,nsReharvestH1,sd)
CattlePerS<-aggregate(TotalBiomass1~ CattleMetBio+Reharvest.date+Treatment,nsReharvestH1,sd)
GrantPerS<-aggregate(TotalBiomass1~ Grants_GazelleMetBio+Reharvest.date+Treatment,nsReharvestH1,sd)

ZebraPer<-cbind(ZebraPer,ZebraPerS[4])
colnames(ZebraPer)[5]<-"se"

# Regrowth different herbivore densities
HerbivoreRegrowth<-ggplot(ZebraPer,aes(x=Burchells_ZebraMetBio,y=TotalBiomass1))
HerbivoreRegrowth<-HerbivoreRegrowth+geom_hline(yintercept =0, linetype="dashed")
HerbivoreRegrowth<-HerbivoreRegrowth+ geom_ribbon(data =MyData5,  aes(x =Burchells_ZebraMetBio,
                                              ymax = SeUp, ymin = SeLo),alpha = 0.5, show.legend = F)
HerbivoreRegrowth<-HerbivoreRegrowth+geom_line(data = MyData5, aes(x = Burchells_ZebraMetBio),size=1,colour="grey30",show.legend=F)
#HerbivoreRegrowth<-HerbivoreRegrowth+geom_errorbar(data=CattlePerS,aes(x=CattleMetBio,ymin=CattlePer$TotalBiomass1-TotalBiomass1,ymax=CattlePer$TotalBiomass1+TotalBiomass1),colour="grey", width=.2, alpha=.6)
HerbivoreRegrowth<-HerbivoreRegrowth+geom_point(data=CattlePer,aes(x=CattleMetBio,y=TotalBiomass1, shape=Treatment),colour="grey", alpha=.6)
HerbivoreRegrowth<-HerbivoreRegrowth+geom_point(data=GrantPer,aes(x=Grants_GazelleMetBio,y=TotalBiomass1, shape=Treatment),colour="grey", alpha=.6)
HerbivoreRegrowth<-HerbivoreRegrowth+geom_point(data=KuduPer,aes(x=Greater_KuduMetBio,y=TotalBiomass1, shape=Treatment),colour="grey")
HerbivoreRegrowth<-HerbivoreRegrowth+geom_errorbar(aes(x=Burchells_ZebraMetBio,ymin=TotalBiomass1-se,ymax=TotalBiomass1+se),colour="black", width=.2)
HerbivoreRegrowth<-HerbivoreRegrowth+geom_point(data=ZebraPer,aes(x=Burchells_ZebraMetBio,y=TotalBiomass1, shape=Treatment), stroke=1,fill="white",size=3,colour="black")
HerbivoreRegrowth<-HerbivoreRegrowth+facet_wrap(~Reharvest.date)
HerbivoreRegrowth<-HerbivoreRegrowth+theme_classic()
HerbivoreRegrowth

# Plot interaction...
par(mfrow=c(1,1))
with(nsReharvestH1 , {interaction.plot(Burchells_ZebraMetBio,Reharvest.date,Tot.Per.diff,
                                  xlab = "tree canopy",
                                  ylab = "Biomass Rainfall correlation",
                                  fun=mean)})




# Grass biomass
G1<-lme(G.Per.diff~Livestock.density,#Treatment+Treatment:Livestock.density,
          random= ~ 1|fBlock,na.action=na.pass, method="ML",data=nsReharvestH1)
#correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1)
summary(G1)
anova(G1)
AIC(G1) #1728

plot(G1) #OK
drop1(G1,test="Chisq") # Burchells_ZebraMetBio  1 1750.1 8.1234  0.00437 **

names(nsReharvestH1)
G1<-lme(G.Per.diff ~ Treatment+Treatment:Burchells_ZebraMetBio + 
          CattleMetBio + Grants_GazelleMetBio + Greater_KuduMetBio,
        random= ~ 1|fBlock,na.action=na.pass, method="ML",
        correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1)
summary(G1)
anova(G1)
AIC(G1) #5

plot(G1) #OK
drop1(G1,test="Chisq") # 0.08311 . marginal

nsReharvestH1W<- nsReharvestH1[is.finite(nsReharvestH1$W.Per.diff),]
W1<-lme(W.Per.diff ~ Treatment+Treatment:Burchells_ZebraMetBio + 
          CattleMetBio + Grants_GazelleMetBio + Greater_KuduMetBio,
        random= ~ 1|fBlock,na.action=na.pass, method="ML",
        correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1W)
summary(W1)
anova(W1) # Interaction Treatment:Livestock.density
AIC(W1) #5

plot(W1) #OK
drop1(W1,test="Chisq") # 0.08311 . marginal


###############################################################################
#### OLD SCRIPT ####
###############################################################################

########################################################################
#### Herbivore count observations ####
########################################################################

nsherb3<-read.table("CountAverageFeb13_Date.txt",header=T,sep="\t")

names(nsherb3)
nsherb3$Zone # Each herbivore transect corresponds to a unique zone unique zone? 

# Spatial correlation
par(mfrow = c(2, 2), mar = c(4, 3, 3, 2))
par(pty = "s", mar = c(5,5,2,2), cex.lab = 1.5)       
plot(x =  nsherb3$Xcent,
     y =  nsherb3$Ycent,
     type = "p",
     main="GazelleMetBio",
     cex=(nsherb3$Grants_GazelleMetBio/40), #nsdung3$Cattle #Burchells_Zebra # Grants_GazelleMetBio
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

plot(x =  nsherb3$Xcent,
     y =  nsherb3$Ycent,
     type = "p",
     main="Burchells_Zebra",
     cex=(nsherb3$Burchells_ZebraMetBio/125), 
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

plot(x =  nsherb3$Xcent,
     y =  nsherb3$Ycent,
     type = "p",
     main="Greater_KuduMetBio",
     cex=(nsherb3$Greater_KuduMetBio/40), 
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

plot(x =  nsherb3$Xcent,
     y =  nsherb3$Ycent,
     type = "p",
     main="CattleMetBio",
     cex=(nsherb3$CattleMetBio/40), 
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

plot(x =  nsherb3$Xcent,
     y =  nsherb3$Ycent,
     type = "p",
     main="CattleMetBio",
     cex=(nsherb3$Swaynes_HartebeestMetBio/40), 
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

plot(x =  nsherb3$Xcent,
     y =  nsherb3$Ycent,
     type = "p",
     main="CattleMetBio",
     cex=(nsherb3$Total_Met/130), 
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

# Explore structure of data - but this is transect...
#ordihull(nsherb3[, c("Xcent", "Ycent")],
#         draw = "polygon",
#         groups = nsherb3[, "Zone"],
#         label = F,
#         col = "red")     
# Zone = each unique herbivore transect

# Need Zuur's functional crib sheet for this
# Collinearity amongst MEtBio
#names(nsherb3)
#MyVar<-c("Burchells_ZebraMetBio","CattleMetBio","Grants_GazelleMetBio","Greater_KuduMetBio","Swaynes_HartebeestMetBio")
#pairs(nsherb3[,MyVar],lower.panel = panel.cor)
# Zebra and gazelle are positively related...0.7
###############################################################################



###############################################################################
##### Spatial model - Biomass regrowth in relation to herbivore metabolic weights##### 
#library(devtools)
#devtools::install_url("https://www.math.ntnu.no/inla/R/stable/bin/macosx/contrib/3.3/INLA_0.0-1483775362.tgz")
#INLA_0.0-1485844051.tgz	# More recent
library(INLA)
library(mgcv)
library(gstat)
library(sp)
library(brinla)
###############################################################

# INLA # Bayesian approach
# All poissible model interactions...from singl eterms in model

# Duplicates of harvest 1 and 2 in time = issues with temporal model
nsReharvestbH1_2$Plot.ID<-as.numeric(nsReharvestbH1_2$Plot.name)
nsReharvestbH1_2$YrMonth<-format(as.Date(nsReharvestbH1_2$Reharvest.date), "%Y-%m")

nsReharvestH1<- droplevels(nsReharvestbH1_2[nsReharvestbH1_2$Harvest=="single",])
table(nsReharvestH1$YrMonth,nsReharvestH1$Plot.ID)
table(nsReharvestH1$Block,nsReharvestH1$Plot.ID)

f1 <-formula(y ~ Treatment*Burchells_ZebraMetBio*CattleMetBio+
               Grants_GazelleMetBio)#*Greater_KuduMetBio)
terms2 <-attr(terms.formula(f1), "term.labels")
terms2 # 18 total - 16 two-way interactions

# Remove Swaynes_HartebeestMetBio - one point measured

# Subset two way
termssub<-terms2[1:8]
termssub<-unlist(termssub)

# Model formulas
fRegrowOrig<-Tot.Per.diff~f(Block, model = "iid")+#f(YrMonth,model="ar1", replicate=Plot.ID)+
  Livestock.density+Treatment+Treatment:Livestock.density:Harvest

fRegrowModT<-as.formula(paste('Tot.Per.diff~f(YrMonth,model="ar1", replicate=Plot.ID)+f(Block, model = "iid")+', paste(termssub, collapse=" + ")))
fRegrowMod<-as.formula(paste('Tot.Per.diff~f(Block, model = "iid")+', paste(termssub, collapse=" + ")))

nsReharvestb$Plot.ID<-as.numeric(nsReharvestb$Plot.name)
nsReharvestb$YrMonth<-format(as.Date(nsReharvestb$Reharvest.date), "%Y-%m")

#### Temporal model ####
RegrowModOT <- inla(fRegrowOrig,
                    family = "gaussian",
                    control.compute = list(waic=TRUE,dic=TRUE),
                    data = nsReharvestH1)


#### Temporal model ####
RegrowModT <- inla(fRegrowModT,
                   family = "gaussian",
                   control.compute = list(waic=TRUE,dic=TRUE),
                   data = nsReharvestH1)

#### No temporal structure ####
RegrowMod <- inla(fRegrowMod,
                  family = "gaussian",
                  control.compute = list(waic=TRUE,dic=TRUE),
                  data = nsReharvestH1)

#RegrowMod <- inla(fRegrowMod,
#                  family = "gamma",
#data=inla.stack.data(stk.e1),
#control.compute = list(waic=T, dic = TRUE),
#control.predictor = list(A = inla.stack.A(stk.e1)))

c(RegrowModOT$waic$waic,RegrowModT$waic$waic,RegrowMod$waic$waic) # No difference
#Gaussian no difference - poor fit
# Gamma better with temporal strcture

#pvalue histogram....
par(mfrow = c(1, 1), mar = c(4, 3, 3, 2))
pval<-rep(NA, nrow=(nsReharvestH1))
for(i in 1:nrow(nsReharvestH1)){
  pval[i]<-inla.pmarginal(q=nsReharvestH1$Tot.Per.diff[i],
                          marginal=RegrowModOT$marginals.fitted.values[[i]])
}
hist(pval) # Bad

bri.fixed.plot(RegrowModOT) # Some herbivore interactions      



# Frequent - use dredge to get best model
nsReharvestH1$fBlock<-as.factor(nsReharvestH1$Block)
nsReharvestH1$fPlot.ID<-as.factor(nsReharvestH1$Plot.ID)

cs1AR1 <- corAR1(0.2, form = ~Reharvest.date|fBlock/fPlot.ID)
cs1AR1. <- Initialize(cs1AR1, data =nsReharvestH1) 
corMatrix(cs1AR1.)


# Total biomass - Create a plot code + harvest - single or double code
nsReharvestb$harvest_code<-as.factor(with(nsReharvestb, paste(Trt.name,Harvest, sep="-")))
levels(nsReharvestb$harvest_code) # 18

nsReharvestavg<-aggregate(TotalBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsReharvest, mean)
nsReharvestsem<-aggregate(TotalBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsReharvest,sd)

# Spatial regrowth double and single harvest
RegrowSp<-ggplot(nsReharvestavg, aes(x=X, y=Y, size=TotalBiomass1,group=Trt.name,shape=Harvest)) 
RegrowSp<-RegrowSp+geom_point()
RegrowSp<-RegrowSp+facet_wrap(~Treatment+Harvest, scale="fixed")
RegrowSp<-RegrowSp+theme_classic()
RegrowSp

# Plot over time
Regrow<-ggplot(nsReharvestavg, aes(x=Harvest.date, y=TotalBiomass1, group=Trt.name,shape=Harvest)) 
Regrow<-Regrow+geom_point(size=3)
Regrow<-Regrow+facet_wrap(~Treatment+Livestock.density, scale="fixed")
Regrow<-Regrow+geom_line(aes(linetype=Trt.name))
Regrow
# Not much change

# Grass biomass - Create a plot code + harvest - single or double code
nsreharvest3Gavg<-aggregate(GrassNetReharvestBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsreharvest3locb, mean)
nsreharvest3Gsem<-aggregate(GrassNetReharvestBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsreharvest3locb,sd)

# Spatial plot grass biomass
# Spatial regrowth double and single harvest - grass biomass
RegrowSpG<-ggplot(nsreharvest3Gavg, aes(x=X, y=Y, size=GrassNetReharvestBiomass1,group=Trt.name,shape=Harvest)) 
RegrowSpG<-RegrowSpG+geom_point()
RegrowSpG<-RegrowSpG+facet_wrap(~Treatment+Harvest, scale="fixed")
RegrowSpG<-RegrowSpG+theme_classic()
RegrowSpG

# Grass biomass regrowth through time
RegrowG<-ggplot(nsreharvest3Gavg, aes(x=Harvest.date, y=GrassNetReharvestBiomass1, group=Trt.name,shape=Harvest)) 
RegrowG<-RegrowG+geom_point(size=3)
RegrowG<-RegrowG+facet_wrap(~Treatment+Livestock.density, scale="fixed")
RegrowG<-RegrowG+geom_line(aes(linetype=Trt.name))
RegrowG

names(nsreharvest3b)
nsreharvest3Wavg<-aggregate(DwarfShrubNetReharvestBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsreharvest3locb, mean)
nsreharvest3Wsem<-aggregate(DwarfShrubNetReharvestBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsreharvest3locb,sd)

# Spatial regrowth double and single harvest - woody biomass
RegrowSpW<-ggplot(nsreharvest3Wavg, aes(x=X, y=Y, size=DwarfShrubNetReharvestBiomass1,group=Trt.name,shape=Harvest)) 
RegrowSpW<-RegrowSpW+geom_point()
RegrowSpW<-RegrowSpW+facet_wrap(~Treatment+Harvest, scale="fixed")
RegrowSpW<-RegrowSpW+theme_classic()
RegrowSpW

# Woody biomass regrowth through time
RegrowW<-ggplot(nsreharvest3Wavg, aes(x=Harvest.date, y=DwarfShrubNetReharvestBiomass1, group=Trt.name,shape=Harvest)) 
RegrowW<-RegrowW+geom_point(size=3)
RegrowW<-RegrowW+facet_wrap(~Treatment+Livestock.density, scale="fixed")
RegrowW<-RegrowW+geom_line(aes(linetype=Trt.name))
RegrowW

#### OLD SCRIPT #### 

#### SPATIAL INLA MODEL ####
# Did not work - insufficient spatial points
#Spatial position of the sampling locations:

# Reproject
utmproj<-"+proj=utm +north +zone=37 +init=EPSG:32637" # central-ethiopia-37n
latlongproj<-("+proj=longlat +datum=WGS84")
rl_proj <-nsReharvestb
coordinates(rl_proj)<- ~X + Y #define which columns correspond to x's and y's
proj4string(rl_proj)<-utmproj #define projection (originally defined above)
extent(rl_proj) 
sme<-spTransform(rl_proj,latlongproj)
extent(sme) # Now in WGS84
class(sme)
nsReharvestb$X<-sme@coords[,1]
nsReharvestb$Y<-sme@coords[,2]

library(lattice)
par(mfrow=c(1,1), mar=c(1,1,1,1))
#coordinates(SppSNPs) <- c("Lat", "Lon")
loc2 <- cbind(nsReharvestb$X, nsReharvestb$Y)

#Not in the book:
xyplot(loc2[,2] ~ loc2[,1],
       aspext = "iso")

# Create a grid mesh
# Issues running mesh with UTM?
mesh5a <- inla.mesh.2d(loc2, max.edge=c(20, 20),cutoff =.001) # Try this

#Install package splancs # Shape files to define boundaries
#library(splancs)
zzdomain2 <- inla.nonconvex.hull(loc2)
mesh6a <- inla.mesh.2d(boundary = zzdomain2, max.edge=c(20, 20), cutoff = 0.1)
# NEed to know what the value of the cut-off equates to

# Plot mesh
plot(mesh5a, asp=1)
points(loc2,col=2,pch=16, cex=.5)

#And mesh is mesh from here onwards
mesh6a$n 
mesh5a$n 
# WGS 84 = 0.01 = 30 need toand 0.001 = 65 - FAR TOO SMALL...~300 - 500
# issues with UTM...huge numbers
# UTM
# 216250 # 0.001
# 844171 # cut off = 1000! - Still orders of mangitude too high
# 864027 # cut off= 0.001 - LOADS!

# Tell INLA which sampling locations match the points
# on the mesh
A1      <- inla.spde.make.A(mesh5a, loc = loc2) # Tells INLA where the sampling locations
# Values 0 and 1 #one equals sampling location
# Also need to do for the covariates - but covariates are not always at the same position

#Define the Matern correlation on the mesh
spde2   <- inla.spde2.matern(mesh5a, alpha = 2) # Quantify distance between points
# Will be explained on the next Powerpoint slide

# Create YrMonth- needs to be factor
nsReharvestb$YrMonth<-format(as.Date(nsReharvestb$Reharvest.date), "%Y-%m")
#nsReharvestb$YrMonth<-as.factor(nsReharvestb$YrMonth)
nsReharvestb$Plot.ID<-as.numeric(nsReharvestb$Plot.name)

# Set up the model. 
# Create a data frame with an intercept and the covariate
names(nsReharvestb)
N2 <- nrow(nsReharvestb)
X2 <- data.frame(Intercept = rep(1,N2), 
                 Plot.ID=nsReharvestb$Plot.ID,
                 Treatment= nsReharvestb$Treatment,
                 Harvest= nsReharvestb$Harvest,
                 YrMonth= nsReharvestb$YrMonth,
                 Burchells_ZebraMetBio = nsReharvestb$Burchells_ZebraMetBio,
                 CattleMetBio=nsReharvestb$CattleMetBio,
                 Grants_GazelleMetBio=nsReharvestb$Grants_GazelleMetBio,
                 Greater_KuduMetBio = nsReharvestb$Greater_KuduMetBio,
                 Swaynes_HartebeestMetBio=nsReharvestb$Swaynes_HartebeestMetBio
                 
) # Covariates
str(X2)
#Tell INLA that the covariates are sampled at the same
#sampling locations.
# MODEL for total biomass - can repeat for woody, grass etc.
stk.e1 <- inla.stack(
  tag = "est",
  data = list(y = nsReharvestb$TotalBiomass1),  # Y vaiable
  A = list(A1,1),      #This is the confusing bit # Sampling locations      
  effects = list(                 
    s = 1:spde2$n.spde,       #Spatial field  
    X2))                      #Covariates
dim(inla.stack.A(stk.e1)) #270  279
########################################

fRegrowMod<-y ~  -1 +f(YrMonth,model='rw2')+ f(s, model=spde2) +
  Treatment+Harvest+ Burchells_ZebraMetBio+CattleMetBio+
  Grants_GazelleMetBio+Greater_KuduMetBio+
  Swaynes_HartebeestMetBio

RegrowMod <- inla(fRegrowMod,
                  family = "gamma",
                  data=inla.stack.data(stk.e1),
                  control.compute = list(waic=T, dic = TRUE),
                  control.predictor = list(A = inla.stack.A(stk.e1)))
########################################




#####################################################
#### END OF RE-EDITS ####
#####################################################


# Raw biomass for each harvest period
tiff("BiomassJune2017.tif",width=12,height=8,units="in",res=100)

par(mfrow=c(1,3))
meanbio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason1,list(Treatment,Livestock.density),mean,na.rm=T))
sembio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason1,list(Treatment,Livestock.density),sd))
b1<-barplot(meanbio,beside=T,main="November 2012",xlab="Livestock density",ylab="",ylim=c(0,1000),las=1)
title(ylab=expression("Biomass gm"^"-2"),line=2)
arrows(b1,meanbio+sembio,b1,meanbio-sembio,code=3,length=0.05,angle=90)

meanbio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason2,list(Treatment,Livestock.density),mean,na.rm=T))
sembio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason2,list(Treatment,Livestock.density),sd))
b1<-barplot(meanbio,beside=T,main="June 2013",xlab="Livestock density",ylab="",ylim=c(0,1000),las=1)
arrows(b1,meanbio+sembio,b1,meanbio-sembio,code=3,length=0.05,angle=90)
title(ylab=expression("Biomass gm"^"-2"),line=2)

meanbio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason3,list(Treatment,Livestock.density),mean,na.rm=T))
sembio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason3,list(Treatment,Livestock.density),sem))
b1<-barplot(meanbio,beside=T,legend.text=c("Grazed","Exclosure"),main="November 2013",xlab="Livestock density",ylab="",ylim=c(0,1000),las=1)
arrows(b1,meanbio+sembio,b1,meanbio-sembio,code=3,length=0.05,angle=90)
title(ylab=expression("Biomass g m"^"-2"),line=2)
dev.off()



lm1<-with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),lm(TotalSeason1~Treatment*Livestock.density))
anova(lm1)
lm2<-update(lm1,.~.-Treatment:Livestock.density)
anova(lm2)
lm1<-with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),lm(TotalSeason2~Treatment*Livestock.density))
anova(lm1)
lm2<-update(lm1,.~.-Treatment:Livestock.density)
anova(lm2)
lm1<-with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),lm(TotalSeason3~Treatment*Livestock.density))
anova(lm1)

require(nlme)
require(plotrix)

nsbiomass3orig$Livestock.density1<-factor(nsbiomass3orig$Livestock.density,ordered=F)
nsbiomass3orig$plottran<-paste(nsbiomass3orig$Transect,nsbiomass3orig$plotid)
#Stacked barplots
{
#Biomass
#Season1
tiff("M:\\Africa\\VegetationDataSeason3\\Figs\\BiomassJuly2017.tif",width=12,height=8,units="in",res=100)
lmeBioGramS1<-with(nsbiomass3orig,lme(GrassNetBiomassSeason1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioGramS1)
lmeBioDwarfshrubS1<-with(nsbiomass3orig,lme(DwarfShrubNetBiomassSeason1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioDwarfshrubS1)
lmeBioHerbS1<-with(nsbiomass3orig,lme(HerbClimb1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioHerbS1)

lmeBioTotallS1<-with(nsbiomass3orig,lme(TotalSeason1~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
anova(lmeBioTotallS1)

exBioS1<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
pBioGramS1<-predict(lmeBioGramS1,exBioS1,level=0)
Designmat <- model.matrix(eval(eval(lmeBioGramS1$call$fixed)[-2]), exBioS1[-4])
predvar <- diag(Designmat %*% lmeBioGramS1$varFix %*% t(Designmat))
exBioS1$GramS1<-pBioGramS1
exBioS1$GramS1SE <- sqrt(predvar)

pBioDwarfshrubS1<-predict(lmeBioDwarfshrubS1,exBioS1,level=0)
Designmat <- model.matrix(eval(eval(lmeBioDwarfshrubS1$call$fixed)[-2]), exBioS1[-4])
predvar <- diag(Designmat %*% lmeBioDwarfshrubS1$varFix %*% t(Designmat))
exBioS1$DwarfshrubS1<-pBioDwarfshrubS1
exBioS1$DwarfshrubS1SE <- sqrt(predvar)

pBioHerbS1<-predict(lmeBioHerbS1,exBioS1,level=0)
Designmat <- model.matrix(eval(eval(lmeBioHerbS1$call$fixed)[-2]), exBioS1[-4])
predvar <- diag(Designmat %*% lmeBioHerbS1$varFix %*% t(Designmat))
exBioS1$HerbS1<-pBioHerbS1
exBioS1$HerbS1SE <- sqrt(predvar)

pBioTotallS1<-predict(lmeBioTotallS1,exBioS1,level=0)
Designmat <- model.matrix(eval(eval(lmeBioTotallS1$call$fixed)[-2]), exBioS1[-4])
predvar <- diag(Designmat %*% lmeBioTotallS1$varFix %*% t(Designmat))
exBioS1$TotallS1<-pBioTotallS1
exBioS1$TotallS1SE <- sqrt(predvar)
exBioS1

par(mfrow=c(1,3))
par(mar=c(12,4,3,1))
b1<-barplot(t(exBioS1[,c(3,5,7)]),space=c(1,0),main="Biomass \n October 2012",ylim=c(0,1000),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
title(ylab=expression("Biomass g m"^-2),line=2.5)
#arrows(b1-0.4,exBioS1$GramS1+exBioS1$GramS1SE,b1-0.4,exBioS1$GramS1-exBioS1$GramS1SE,code=3,length=0.05,angle=90)
#arrows(b1-0.2,exBioS1$GramS1+exBioS1$DwarfshrubS1+exBioS1$DwarfshrubS1SE,b1-0.2,exBioS1$GramS1+exBioS1$DwarfshrubS1-exBioS1$DwarfshrubS1SE,code=3,length=0.05,angle=90)
#arrows(b1+0.2,exBioS1$GramS1+exBioS1$DwarfshrubS1+exBioS1$HerbS1+exBioS1$HerbS1SE,b1+0.2,exBioS1$GramS1+exBioS1$DwarfshrubS1+exBioS1$HerbS1-exBioS1$HerbS1SE,code=3,length=0.05,angle=90)
arrows(b1,exBioS1$TotallS1+exBioS1$TotallS1SE,b1,exBioS1$TotallS1-exBioS1$TotallS1SE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
lmeBioTotallS1_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeBioTotallS1),3))
lmeBioTotallS1_Tab[which(lmeBioTotallS1_Tab[,5]=="0",),5] <- "<0.001"
#addtable2plot(5,800,lmeBioTotallS1_Tab,cex=0.8)
mtext(side=3,adj=0,"(a)",line=1)

#Season2
lmeBioGramS2<-with(nsbiomass3orig,lme(GrassNetBiomassSeason2~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioGramS2)
lmeBioDwarfshrubS2<-with(nsbiomass3orig,lme(DwarfShrubNetBiomassSeason2~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioDwarfshrubS2)
lmeBioHerbS2<-with(nsbiomass3orig,lme(HerbClimb2~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioHerbS2)
lmeBioTotallS2<-with(nsbiomass3orig,lme(TotalSeason2~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
anova(lmeBioTotallS2)

exBioS2<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
pBioGramS2<-predict(lmeBioGramS2,exBioS2,level=0)
Designmat <- model.matrix(eval(eval(lmeBioGramS2$call$fixed)[-2]), exBioS2[-4])
predvar <- diag(Designmat %*% lmeBioGramS2$varFix %*% t(Designmat))
exBioS2$GramS2<-pBioGramS2
exBioS2$GramS2SE <- sqrt(predvar)

pBioDwarfshrubS2<-predict(lmeBioDwarfshrubS2,exBioS2,level=0)
Designmat <- model.matrix(eval(eval(lmeBioDwarfshrubS2$call$fixed)[-2]), exBioS2[-4])
predvar <- diag(Designmat %*% lmeBioDwarfshrubS2$varFix %*% t(Designmat))
exBioS2$DwarfshrubS2<-pBioDwarfshrubS2
exBioS2$DwarfshrubS2SE <- sqrt(predvar)

pBioHerbS2<-predict(lmeBioHerbS2,exBioS2,level=0)
Designmat <- model.matrix(eval(eval(lmeBioHerbS2$call$fixed)[-2]), exBioS2[-4])
predvar <- diag(Designmat %*% lmeBioHerbS2$varFix %*% t(Designmat))
exBioS2$HerbS2<-pBioHerbS2
exBioS2$HerbS2SE <- sqrt(predvar)

pBioTotallS2<-predict(lmeBioTotallS2,exBioS2,level=0)
Designmat <- model.matrix(eval(eval(lmeBioTotallS2$call$fixed)[-2]), exBioS2[-4])
predvar <- diag(Designmat %*% lmeBioTotallS2$varFix %*% t(Designmat))
exBioS2$TotallS2<-pBioTotallS2
exBioS2$TotallS2SE <- sqrt(predvar)
exBioS2

par(mar=c(12,4,3,1))
b1<-barplot(t(exBioS2[,c(3,5,7)]),space=c(1,0),main="Biomass \n June 2013",ylim=c(0,1000),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
#b1<-barplot(t(exBioS2[,c(9)]),space=c(1,0),main="Biomass \n June 2013",ylim=c(0,900),las=1,xlab="",col=exBioS2$Treatment
title(ylab=expression("Biomass g m"^-2),line=2.5)
#arrows(b1-0.4,exBioS2$GramS2+exBioS2$GramS2SE,b1-0.4,exBioS2$GramS2-exBioS2$GramS2SE,code=3,length=0.05,angle=90)
#arrows(b1-0.2,exBioS2$GramS2+exBioS2$DwarfshrubS2+exBioS2$DwarfshrubS2SE,b1-0.2,exBioS2$GramS2+exBioS2$DwarfshrubS2-exBioS2$DwarfshrubS2SE,code=3,length=0.05,angle=90)
#arrows(b1+0.2,exBioS2$GramS2+exBioS2$DwarfshrubS2+exBioS2$HerbS2+exBioS2$HerbS2SE,b1+0.2,exBioS2$GramS2+exBioS2$DwarfshrubS2+exBioS2$HerbS2-exBioS2$HerbS2SE,code=3,length=0.05,angle=90)
arrows(b1,exBioS2$TotallS2+exBioS2$TotallS2SE,b1,exBioS2$TotallS2-exBioS2$TotallS2SE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
#lmeBioTotallS2_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeBioTotallS2),3))
#lmeBioTotallS2_Tab[which(lmeBioTotallS2_Tab[,5]=="0",),5] <- "<0.001"
#addtable2plot(5,800,lmeBioTotallS2_Tab,cex=0.8)
mtext(side=3,adj=0,"(b)",line=1)
#Season3
lmeBioGramS3<-with(nsbiomass3orig,lme(GrassNetBiomassSeason3~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioGramS3)
lmeBioDwarfshrubS3<-with(nsbiomass3orig,lme(DwarfShrubNetBiomassSeason3~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioDwarfshrubS3)
lmeBioHerbS3<-with(nsbiomass3orig,lme(HerbClimb3~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeBioHerbS3)
lmeBioTotallS3<-with(nsbiomass3orig,lme(TotalSeason3~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
anova(lmeBioTotallS3)

exBioS3<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
pBioGramS3<-predict(lmeBioGramS3,exBioS3,level=0)
Designmat <- model.matrix(eval(eval(lmeBioGramS3$call$fixed)[-2]), exBioS3[-4])
predvar <- diag(Designmat %*% lmeBioGramS3$varFix %*% t(Designmat))
exBioS3$GramS3<-pBioGramS3
exBioS3$GramS3SE <- sqrt(predvar)

pBioDwarfshrubS3<-predict(lmeBioDwarfshrubS3,exBioS3,level=0)
Designmat <- model.matrix(eval(eval(lmeBioDwarfshrubS3$call$fixed)[-2]), exBioS3[-4])
predvar <- diag(Designmat %*% lmeBioDwarfshrubS3$varFix %*% t(Designmat))
exBioS3$DwarfshrubS3<-pBioDwarfshrubS3
exBioS3$DwarfshrubS3SE <- sqrt(predvar)

pBioHerbS3<-predict(lmeBioHerbS3,exBioS3,level=0)
Designmat <- model.matrix(eval(eval(lmeBioHerbS3$call$fixed)[-2]), exBioS3[-4])
predvar <- diag(Designmat %*% lmeBioHerbS3$varFix %*% t(Designmat))
exBioS3$HerbS3<-pBioHerbS3
exBioS3$HerbS3SE <- sqrt(predvar)

pBioTotallS3<-predict(lmeBioTotallS3,exBioS3,level=0)
Designmat <- model.matrix(eval(eval(lmeBioTotallS3$call$fixed)[-2]), exBioS3[-4])
predvar <- diag(Designmat %*% lmeBioTotallS3$varFix %*% t(Designmat))
exBioS3$TotallS3<-pBioTotallS3
exBioS3$TotallS3SE <- sqrt(predvar)
exBioS3

par(mar=c(12,4,3,1))
b1<-barplot(t(exBioS3[,c(3,5,7)]),space=c(1,0),main="Biomass \n November 2013",ylim=c(0,1000),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
title(ylab=expression("Biomass g m"^-2),line=2.5)
#arrows(b1-0.4,exBioS3$GramS3+exBioS3$GramS3SE,b1-0.4,exBioS3$GramS3-exBioS3$GramS3SE,code=3,length=0.05,angle=90)
#arrows(b1-0.2,exBioS3$GramS3+exBioS3$DwarfshrubS3+exBioS3$DwarfshrubS3SE,b1-0.2,exBioS3$GramS3+exBioS3$DwarfshrubS3-exBioS3$DwarfshrubS3SE,code=3,length=0.05,angle=90)
#arrows(b1+0.2,exBioS3$GramS3+exBioS3$DwarfshrubS3+exBioS3$HerbS3+exBioS3$HerbS3SE,b1+0.2,exBioS3$GramS3+exBioS3$DwarfshrubS3+exBioS3$HerbS3-exBioS3$HerbS3SE,code=3,length=0.05,angle=90)
arrows(b1,exBioS3$TotallS3+exBioS3$TotallS3SE,b1,exBioS3$TotallS3-exBioS3$TotallS3SE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
lmeBioTotallS3_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeBioTotallS3),3))
lmeBioTotallS3_Tab[which(lmeBioTotallS3_Tab[,5]=="0",),5] <- "<0.001"
#addtable2plot(5,800,lmeBioTotallS3_Tab,cex=0.8)
mtext(side=3,adj=0,"(c)",line=1)
legend(5,600,fill=c(grey(0.9),grey(0.6),grey(0.3)),c("Herbs","Dwarf shrubs","Graminoids"))
dev.off()
}


lmeBioTotallS1<-with(nsbiomass3orig,lme(TotalSeason1~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
anova(lmeBioTotallS1)
lmeBioTotallS2<-with(nsbiomass3orig,lme(TotalSeason2~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
anova(lmeBioTotallS2)
lmeBioTotallS3<-with(nsbiomass3orig,lme(TotalSeason3~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
anova(lmeBioTotallS3)



#################################################################################
nsreharvest3<-read.table("M:\\Africa\\VegetationDataSeason3\\ProductivitySeason3.txt",header=T,sep="\t")

head(nsreharvest3)
#Convert to gm^-2
nsreharvest3$HerbClimbReharv1<-nsreharvest3$HerbNetReharvestBiomass1+nsreharvest3$ClimberNetReharvestBiomass1
nsreharvest3$HerbClimbReharv2<-nsreharvest3$HerbNetReharvestBiomass2+nsreharvest3$ClimberNetReharvestBiomass2
nsreharvest3[,9:18]<-nsreharvest3[,9:18]*4
nsreharvest3$TotalProductivityS1<-rowSums(nsreharvest3[,9:12])
nsreharvest3$TotalProductivityS2<-rowSums(nsreharvest3[,13:16])
nsreharvest3$Livestock.density1<-factor(nsreharvest3$Livestock.density,levels=c("Low","Medium","High"),ordered=T)
nsreharvest3orig<-droplevels(nsreharvest3[nsreharvest3$Treatment=="Control" | nsreharvest3$Treatment=="Exclosure",])
nsreharvest3orig$plottran<-paste(nsreharvest3orig$Transect,nsreharvest3orig$Treatment,nsreharvest3orig$Livestock.density)



par(mfrow=c(1,2))
reharvmean<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],tapply(TotalProductivityS1,list(Treatment,Livestock.density1),mean))
reharvsem<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],tapply(TotalProductivityS1,list(Treatment,Livestock.density1),sem))
b1<-barplot(reharvmean,beside=T,ylim=c(0,550),ylab="",xlab="Livestock density",legend.text=c("Grazed","Exclosure"),main="Productivity \n Oct 2012 - April 2013" ,las=1 )
arrows(b1,reharvmean+reharvsem,b1,reharvmean-reharvsem,length=0.05,angle=90,code=3)
title(ylab=(expression(paste("Productivity gm"^"-2", "season"^"-1"))),line=2.5)

#reharvmean<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],tapply(TotalProductivityS2,list(Treatment,Livestock.density1),mean))
#reharvsem<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],tapply(TotalProductivityS2,list(Treatment,Livestock.density1),sem))
#b1<-barplot(reharvmean,beside=T,ylim=c(0,550),ylab="",xlab="Livestock density",legend.text=c("Grazed","Exclosure"),main="Second Reharvest \n Nov 2012 - April 2013 - Nov 2013"  )
#arrows(b1,reharvmean+reharvsem,b1,reharvmean-reharvsem,length=0.05,angle=90,code=3)
#title(ylab=(expression(paste("Productivity gm"^"-2", "season"^"-1"))),line=2.5)

reharvmean<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],tapply(TotalProductivityS1,list(Treatment,Livestock.density1),mean,na.rm=T))
reharvsem<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],tapply(TotalProductivityS1,list(Treatment,Livestock.density1),sem))
b1<-barplot(reharvmean,beside=T,ylim=c(0,550),ylab="",xlab="Livestock density",legend.text=c("Grazed","Exclosure"),main="Productivity \n June 2013 - Nov 2013",las=1 )
arrows(b1,reharvmean+reharvsem,b1,reharvmean-reharvsem,length=0.05,angle=90,code=3)
title(ylab=(expression(paste("Productivity gm"^"-2", "season"^"-1"))),line=2.5)

lm1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],lm(TotalProductivityS1~Treatment*Livestock.density))
anova(lm1)
lm1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],lm(TotalProductivityS2~Treatment*Livestock.density))
anova(lm1)
lm1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],lm(TotalProductivityS1~Treatment*Livestock.density))
summary(lm1)
anova(lm1)


#Stacked
nsreharvest3orig$Livestock.density1<-factor(nsreharvest3orig$Livestock.density,ordered=F)
{
tiff("M:\\Africa\\VegetationDataSeason3\\Figs\\Prod_July2017.tif",width=12,height=8,units="in",res=100)
#First Reharvest \n Nov 2012 - April 2013
  lmeProdGramS1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],lme(GrassNetReharvestBiomass1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
  anova(lmeProdGramS1)
  lmeProdDwarfshrubS1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],lme(DwarfShrubNetReharvestBiomass1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
  anova(lmeProdDwarfshrubS1)
  lmeProdHerbS1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],lme(HerbClimbReharv1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
  anova(lmeProdHerbS1)
  
  lmeProdTotallS1<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="Nov.25-27,2012",],lme(TotalProductivityS1~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
  anova(lmeProdTotallS1)
  
  exProdS1<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
  pProdGramS1<-predict(lmeProdGramS1,exProdS1,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdGramS1$call$fixed)[-2]), exProdS1[-4])
  predvar <- diag(Designmat %*% lmeProdGramS1$varFix %*% t(Designmat))
  exProdS1$GramS1<-pProdGramS1
  exProdS1$GramS1SE <- sqrt(predvar)
  
  pProdDwarfshrubS1<-predict(lmeProdDwarfshrubS1,exProdS1,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdDwarfshrubS1$call$fixed)[-2]), exProdS1[-4])
  predvar <- diag(Designmat %*% lmeProdDwarfshrubS1$varFix %*% t(Designmat))
  exProdS1$DwarfshrubS1<-pProdDwarfshrubS1
  exProdS1$DwarfshrubS1SE <- sqrt(predvar)
  
  pProdHerbS1<-predict(lmeProdHerbS1,exProdS1,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdHerbS1$call$fixed)[-2]), exProdS1[-4])
  predvar <- diag(Designmat %*% lmeProdHerbS1$varFix %*% t(Designmat))
  exProdS1$HerbS1<-pProdHerbS1
  exProdS1$HerbS1SE <- sqrt(predvar)
  
  pProdTotallS1<-predict(lmeProdTotallS1,exProdS1,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdTotallS1$call$fixed)[-2]), exProdS1[-4])
  predvar <- diag(Designmat %*% lmeProdTotallS1$varFix %*% t(Designmat))
  exProdS1$TotallS1<-pProdTotallS1
  exProdS1$TotallS1SE <- sqrt(predvar)
  exProdS1
  
  par(mfrow=c(1,2))
  par(mar=c(12,4,3,1))
  b1<-barplot(t(exProdS1[,c(3,5,7)]),space=c(1,0),main="Regrowth \n Oct 2012 to April 2013",ylim=c(0,550),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
  title(ylab=expression(paste("Productivity g m"^-2," season"^-1)),line=2.5)
# arrows(b1-0.4,exProdS1$GramS1+exProdS1$GramS1SE,b1-0.4,exProdS1$GramS1-exProdS1$GramS1SE,code=3,length=0.05,angle=90)
#  arrows(b1-0.2,exProdS1$GramS1+exProdS1$DwarfshrubS1+exProdS1$DwarfshrubS1SE,b1-0.2,exProdS1$GramS1+exProdS1$DwarfshrubS1-exProdS1$DwarfshrubS1SE,code=3,length=0.05,angle=90)
#  arrows(b1+0.2,exProdS1$GramS1+exProdS1$DwarfshrubS1+exProdS1$HerbS1+exProdS1$HerbS1SE,b1+0.2,exProdS1$GramS1+exProdS1$DwarfshrubS1+exProdS1$HerbS1-exProdS1$HerbS1SE,code=3,length=0.05,angle=90)
  arrows(b1,exProdS1$TotallS1+exProdS1$TotallS1SE,b1,exProdS1$TotallS1-exProdS1$TotallS1SE,code=3,length=0.05,angle=90,lwd=2)
  mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
  lmeProdTotallS1_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeProdTotallS1),3))
  lmeProdTotallS1_Tab[which(lmeProdTotallS1_Tab[,5]=="0",),5] <- "<0.001"
#  addtable2plot(5,450,lmeProdTotallS1_Tab,cex=0.8)
  mtext(side=3,adj=0,"(a)",line=1)
  
   #First Reharvest \n April 2013 - Nov 2013"
  lmeProdGramS1a<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],lme(GrassNetReharvestBiomass1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
  anova(lmeProdGramS1a)
  lmeProdDwarfshrubS1a<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],lme(DwarfShrubNetReharvestBiomass1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
  anova(lmeProdDwarfshrubS1a)
  lmeProdHerbS1a<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],lme(HerbClimbReharv1~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
  anova(lmeProdHerbS1a)
  
  lmeProdTotallS1a<-with(nsreharvest3orig[nsreharvest3orig$Date.First.Harvested=="June21-25,2013",],lme(TotalProductivityS1~Treatment*Livestock.density1,random=~1|plottran,na.action=na.omit))
  anova(lmeProdTotallS1a)
  
  exProdS1a<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
  pProdGramS1a<-predict(lmeProdGramS1a,exProdS1a,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdGramS1a$call$fixed)[-2]), exProdS1a[-4])
  predvar <- diag(Designmat %*% lmeProdGramS1a$varFix %*% t(Designmat))
  exProdS1a$GramS1a<-pProdGramS1a
  exProdS1a$GramS1aSE <- sqrt(predvar)
  
  pProdDwarfshrubS1a<-predict(lmeProdDwarfshrubS1a,exProdS1a,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdDwarfshrubS1a$call$fixed)[-2]), exProdS1a[-4])
  predvar <- diag(Designmat %*% lmeProdDwarfshrubS1a$varFix %*% t(Designmat))
  exProdS1a$DwarfshrubS1a<-pProdDwarfshrubS1a
  exProdS1a$DwarfshrubS1aSE <- sqrt(predvar)
  
  pProdHerbS1a<-predict(lmeProdHerbS1a,exProdS1a,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdHerbS1a$call$fixed)[-2]), exProdS1a[-4])
  predvar <- diag(Designmat %*% lmeProdHerbS1a$varFix %*% t(Designmat))
  exProdS1a$HerbS1a<-pProdHerbS1a
  exProdS1a$HerbS1aSE <- sqrt(predvar)
  
  pProdTotallS1a<-predict(lmeProdTotallS1a,exProdS1a,level=0)
  Designmat <- model.matrix(eval(eval(lmeProdTotallS1a$call$fixed)[-2]), exProdS1a[-4])
  predvar <- diag(Designmat %*% lmeProdTotallS1a$varFix %*% t(Designmat))
  exProdS1a$TotallS1a<-pProdTotallS1a
  exProdS1a$TotallS1aSE <- sqrt(predvar)
  exProdS1a

par(mar=c(12,4,3,1))
b1<-barplot(t(exProdS1a[,c(3,5,7)]),space=c(1,0),main="Regrowth \n June 2013 - Nov 2013",ylim=c(0,550),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
title(ylab=expression(paste("Productivity g m"^-2," season"^-1)),line=2.5)
#arrows(b1-0.4,exProdS1a$GramS1a+exProdS1a$GramS1aSE,b1-0.4,exProdS1a$GramS1a-exProdS1a$GramS1aSE,code=3,length=0.05,angle=90)
#arrows(b1-0.2,exProdS1a$GramS1a+exProdS1a$DwarfshrubS1a+exProdS1a$DwarfshrubS1aSE,b1-0.2,exProdS1a$GramS1a+exProdS1a$DwarfshrubS1a-exProdS1a$DwarfshrubS1aSE,code=3,length=0.05,angle=90)
#arrows(b1+0.2,exProdS1a$GramS1a+exProdS1a$DwarfshrubS1a+exProdS1a$HerbS1a+exProdS1a$HerbS1aSE,b1+0.2,exProdS1a$GramS1a+exProdS1a$DwarfshrubS1a+exProdS1a$HerbS1a-exProdS1a$HerbS1aSE,code=3,length=0.05,angle=90)
arrows(b1,exProdS1a$TotallS1a+exProdS1a$TotallS1aSE,b1,exProdS1a$TotallS1a-exProdS1a$TotallS1aSE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
#lmeProdTotallS1a_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeProdTotallS1a),3))
#lmeProdTotallS1a_Tab[which(lmeProdTotallS1a_Tab[,5]=="0",),5] <- "<0.001"
#addtable2plot(5,450,lmeProdTotallS1a_Tab,cex=0.8)
mtext(side=3,adj=0,"(b)",line=1)
legend(5,400,fill=c(grey(0.9),grey(0.6),grey(0.3)),c("Herbs","Dwarf shrubs","Graminoids"))
dev.off()
}


anova(lmeProdTotallS1)
anova(lmeProdTotallS1a)

#################################################################################
#Composition
require(vegan)

nsspprich<-read.table("M:\\Africa\\VegetationDataSeason3\\spprich.txt",header=T,sep="\t")  
head(nsspprich)

nsspprich$ShrubRich<-nsspprich$Shrub_spprich+nsspprich$Dwarfshrub_spprich
nsspprich$GramRich<-nsspprich$Grass_spprich+nsspprich$Sedge_spprich
nsspprich$HerbClimbRich<-nsspprich$Herb_spprich+nsspprich$Climber_spprich

nsspprich$Livestock.density<-factor(nsspprich$Livestock.density,levels=c("Low","Medium","High"),ordered=T)
nsspprich$Livestock.density1<-factor(nsspprich$Livestock.density,levels=c("Low","Medium","High"),ordered=F)
nsspprich3orig<-droplevels(nsspprich[nsspprich$Treatment=="Control" | nsspprich$Treatment=="Exclosure",])

#Stacked bar
tiff("M:\\Africa\\VegetationDataSeason3\\Figs\\SppRich.tif",width=12,height=8,units="in",res=100)
#Season1
lmeSppGramS1<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season I",],lme(GramRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppGramS1)
lmeSppShrubS1<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season I",],lme(ShrubRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppShrubS1)
lmeSppHerbS1<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season I",],lme(HerbClimbRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppHerbS1)
lmeSppTotallS1<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season I",],lme(Total_spprich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppTotallS1)

exSppS1<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
pSppGramS1<-predict(lmeSppGramS1,exSppS1,level=0)
Designmat <- model.matrix(eval(eval(lmeSppGramS1$call$fixed)[-2]), exSppS1[-4])
predvar <- diag(Designmat %*% lmeSppGramS1$varFix %*% t(Designmat))
exSppS1$GramS1<-pSppGramS1
exSppS1$GramS1SE <- sqrt(predvar)

pSppShrubsS1<-predict(lmeSppShrubS1,exSppS1,level=0)
Designmat <- model.matrix(eval(eval(lmeSppShrubS1$call$fixed)[-2]), exSppS1[-4])
predvar <- diag(Designmat %*% lmeSppShrubS1$varFix %*% t(Designmat))
exSppS1$ShrubsS1<-pSppShrubsS1
exSppS1$ShrubsS1SE <- sqrt(predvar)

pSppHerbS1<-predict(lmeSppHerbS1,exSppS1,level=0)
Designmat <- model.matrix(eval(eval(lmeSppHerbS1$call$fixed)[-2]), exSppS1[-4])
predvar <- diag(Designmat %*% lmeSppHerbS1$varFix %*% t(Designmat))
exSppS1$HerbS1<-pSppHerbS1
exSppS1$HerbS1SE <- sqrt(predvar)

pSppTotallS1<-predict(lmeSppTotallS1,exSppS1,level=0)
Designmat <- model.matrix(eval(eval(lmeSppTotallS1$call$fixed)[-2]), exSppS1[-4])
predvar <- diag(Designmat %*% lmeSppTotallS1$varFix %*% t(Designmat))
exSppS1$TotallS1<-pSppTotallS1
exSppS1$TotallS1SE <- sqrt(predvar)
exSppS1

par(mfrow=c(1,3))
par(mar=c(12,4,3,1))
b1<-barplot(t(exSppS1[,c(3,5,7)]),space=c(1,0),main="Species richness \n October 2012",ylab="Species richness",ylim=c(0,12),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
arrows(b1-0.4,exSppS1$GramS1+exSppS1$GramS1SE,b1-0.4,exSppS1$GramS1-exSppS1$GramS1SE,code=3,length=0.05,angle=90)
arrows(b1-0.2,exSppS1$GramS1+exSppS1$ShrubsS1+exSppS1$ShrubsS1SE,b1-0.2,exSppS1$GramS1+exSppS1$ShrubsS1-exSppS1$ShrubsS1SE,code=3,length=0.05,angle=90)
arrows(b1+0.2,exSppS1$GramS1+exSppS1$ShrubsS1+exSppS1$HerbS1+exSppS1$HerbS1SE,b1+0.2,exSppS1$GramS1+exSppS1$ShrubsS1+exSppS1$HerbS1-exSppS1$HerbS1SE,code=3,length=0.05,angle=90)
arrows(b1,exSppS1$TotallS1+exSppS1$TotallS1SE,b1,exSppS1$TotallS1-exSppS1$TotallS1SE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
lmeSppTotallS1_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeSppTotallS1),3))
lmeSppTotallS1_Tab[which(lmeSppTotallS1_Tab[,5]=="0",),5] <- "<0.001"
addtable2plot(5,10,lmeSppTotallS1_Tab,cex=0.8)
mtext(side=3,adj=0,"(a)",line=1)

#Season 2
lmeSppGramS2<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season II",],lme(GramRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppGramS2)
lmeSppShrubS2<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season II",],lme(ShrubRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppShrubS2)
lmeSppHerbS2<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season II",],lme(HerbClimbRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppHerbS2)
lmeSppTotallS2<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season II",],lme(Total_spprich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppTotallS2)

exSppS2<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
pSppGramS2<-predict(lmeSppGramS2,exSppS2,level=0)
Designmat <- model.matrix(eval(eval(lmeSppGramS2$call$fixed)[-2]), exSppS2[-4])
predvar <- diag(Designmat %*% lmeSppGramS2$varFix %*% t(Designmat))
exSppS2$GramS2<-pSppGramS2
exSppS2$GramS2SE <- sqrt(predvar)

pSppShrubsS2<-predict(lmeSppShrubS2,exSppS2,level=0)
Designmat <- model.matrix(eval(eval(lmeSppShrubS2$call$fixed)[-2]), exSppS2[-4])
predvar <- diag(Designmat %*% lmeSppShrubS2$varFix %*% t(Designmat))
exSppS2$ShrubsS2<-pSppShrubsS2
exSppS2$ShrubsS2SE <- sqrt(predvar)

pSppHerbS2<-predict(lmeSppHerbS2,exSppS2,level=0)
Designmat <- model.matrix(eval(eval(lmeSppHerbS2$call$fixed)[-2]), exSppS2[-4])
predvar <- diag(Designmat %*% lmeSppHerbS2$varFix %*% t(Designmat))
exSppS2$HerbS2<-pSppHerbS2
exSppS2$HerbS2SE <- sqrt(predvar)

pSppTotallS2<-predict(lmeSppTotallS2,exSppS2,level=0)
Designmat <- model.matrix(eval(eval(lmeSppTotallS2$call$fixed)[-2]), exSppS2[-4])
predvar <- diag(Designmat %*% lmeSppTotallS2$varFix %*% t(Designmat))
exSppS2$TotallS2<-pSppTotallS2
exSppS2$TotallS2SE <- sqrt(predvar)
exSppS2

par(mar=c(12,4,3,1))
b1<-barplot(t(exSppS2[,c(3,5,7)]),space=c(1,0),main="Species richness \n June 2013",ylab="Species richness",ylim=c(0,12),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
arrows(b1-0.4,exSppS2$GramS2+exSppS2$GramS2SE,b1-0.4,exSppS2$GramS2-exSppS2$GramS2SE,code=3,length=0.05,angle=90)
arrows(b1-0.2,exSppS2$GramS2+exSppS2$ShrubsS2+exSppS2$ShrubsS2SE,b1-0.2,exSppS2$GramS2+exSppS2$ShrubsS2-exSppS2$ShrubsS2SE,code=3,length=0.05,angle=90)
arrows(b1+0.2,exSppS2$GramS2+exSppS2$ShrubsS2+exSppS2$HerbS2+exSppS2$HerbS2SE,b1+0.2,exSppS2$GramS2+exSppS2$ShrubsS2+exSppS2$HerbS2-exSppS2$HerbS2SE,code=3,length=0.05,angle=90)
arrows(b1,exSppS2$TotallS2+exSppS2$TotallS2SE,b1,exSppS2$TotallS2-exSppS2$TotallS2SE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
lmeSppTotallS2_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeSppTotallS2),3))
lmeSppTotallS2_Tab[which(lmeSppTotallS2_Tab[,5]=="0",),5] <- "<0.001"
addtable2plot(5,10,lmeSppTotallS2_Tab,cex=0.8)
mtext(side=3,adj=0,"(a)",line=1)

#Season 3
lmeSppGramS3<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season III",],lme(GramRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppGramS3)
lmeSppShrubS3<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season III",],lme(ShrubRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppShrubS3)
lmeSppHerbS3<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season III",],lme(HerbClimbRich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppHerbS3)
lmeSppTotallS3<-with(nsspprich3orig[nsspprich3orig$Season.Year=="Season III",],lme(Total_spprich~Treatment*Livestock.density1,random=~1|Transect,na.action=na.omit))
anova(lmeSppTotallS3)

exSppS3<-expand.grid(Treatment=c("Control","Exclosure"),Livestock.density1=c("Low","Medium","High"))
pSppGramS3<-predict(lmeSppGramS3,exSppS3,level=0)
Designmat <- model.matrix(eval(eval(lmeSppGramS3$call$fixed)[-2]), exSppS3[-4])
predvar <- diag(Designmat %*% lmeSppGramS3$varFix %*% t(Designmat))
exSppS3$GramS3<-pSppGramS3
exSppS3$GramS3SE <- sqrt(predvar)

pSppShrubsS3<-predict(lmeSppShrubS3,exSppS3,level=0)
Designmat <- model.matrix(eval(eval(lmeSppShrubS3$call$fixed)[-2]), exSppS3[-4])
predvar <- diag(Designmat %*% lmeSppShrubS3$varFix %*% t(Designmat))
exSppS3$ShrubsS3<-pSppShrubsS3
exSppS3$ShrubsS3SE <- sqrt(predvar)

pSppHerbS3<-predict(lmeSppHerbS3,exSppS3,level=0)
Designmat <- model.matrix(eval(eval(lmeSppHerbS3$call$fixed)[-2]), exSppS3[-4])
predvar <- diag(Designmat %*% lmeSppHerbS3$varFix %*% t(Designmat))
exSppS3$HerbS3<-pSppHerbS3
exSppS3$HerbS3SE <- sqrt(predvar)

pSppTotallS3<-predict(lmeSppTotallS3,exSppS3,level=0)
Designmat <- model.matrix(eval(eval(lmeSppTotallS3$call$fixed)[-2]), exSppS3[-4])
predvar <- diag(Designmat %*% lmeSppTotallS3$varFix %*% t(Designmat))
exSppS3$TotallS3<-pSppTotallS3
exSppS3$TotallS3SE <- sqrt(predvar)
exSppS3

par(mar=c(12,4,3,1))
b1<-barplot(t(exSppS3[,c(3,5,7)]),space=c(1,0),main="Species richness \n November 2013",ylab="Species richness",ylim=c(0,12),las=1,xlab="",col=c(grey(0.3),grey(0.6),grey(0.9)))
arrows(b1-0.4,exSppS3$GramS3+exSppS3$GramS3SE,b1-0.4,exSppS3$GramS3-exSppS3$GramS3SE,code=3,length=0.05,angle=90)
arrows(b1-0.2,exSppS3$GramS3+exSppS3$ShrubsS3+exSppS3$ShrubsS3SE,b1-0.2,exSppS3$GramS3+exSppS3$ShrubsS3-exSppS3$ShrubsS3SE,code=3,length=0.05,angle=90)
arrows(b1+0.2,exSppS3$GramS3+exSppS3$ShrubsS3+exSppS3$HerbS3+exSppS3$HerbS3SE,b1+0.2,exSppS3$GramS3+exSppS3$ShrubsS3+exSppS3$HerbS3-exSppS3$HerbS3SE,code=3,length=0.05,angle=90)
arrows(b1,exSppS3$TotallS3+exSppS3$TotallS3SE,b1,exSppS3$TotallS3-exSppS3$TotallS3SE,code=3,length=0.05,angle=90,lwd=2)
mtext(1,text=c(paste("Low livestock\nGrazed"),paste("Low livestock\nExclosed"),paste("Medium livestock\nGrazed"),paste("Medium livestock\nExclosed"),paste("High livestock\nGrazed"),paste("High livestock\nExclosed")),srt=90,line=1,at=b1,las=2)
lmeSppTotallS3_Tab<-cbind(Factor=c("I","T","L","T:L"),round(anova(lmeSppTotallS3),3))
lmeSppTotallS3_Tab[which(lmeSppTotallS3_Tab[,5]=="0",),5] <- "<0.001"
addtable2plot(5,10,lmeSppTotallS3_Tab,cex=0.8)
mtext(side=3,adj=0,"(c)",line=1)
legend(1,10,fill=c(grey(0.9),grey(0.6),grey(0.3)),c("Herbs","Dwarf shrubs","Graminoids"))
dev.off()

###############
#Relating to herbivore density

herbdens<-read.table("M:\\Africa\\Herbivore Counts\\CountAverageFeb13.txt",header=T)

head(nsbiomass3)
nsbiomass3$pairID<-factor(with(nsbiomass3,paste(Transect,Livestock.density)))
#nsbiomass3$Zone<-nsbiomass3$pairID

#levels(nsbiomass3$Zone)<-c('D2','A1','C2','D3','B2','C3','D4','A3','C4')
#write.table(nsbiomass3,"M:\\Africa\\VegetationDataSeason3\\BiomassSeason3a.txt")

nsbiomass3_herbio<-merge(nsbiomass3,herbdens,by.x='Zone',by.y='Zone',all.x=T,all.y=F)

dfSite<-droplevels(nsbiomass3_herbio[nsbiomass3_herbio$Treatment=="Control" & nsbiomass3_herbio$Replicate==1,])
meanherbs<-aggregate.data.frame(dfSite[,30:52],list(Livestock.density=dfSite$Livestock.density),mean)
seherbs<-aggregate.data.frame(dfSite[,30:52],list(Livestock.density=dfSite$Livestock.density),sem)

tiff("M:\\Africa\\VegetationDataSeason3\\Figs\\GrazerDensity.tif",width=8,height=6,units="in",res=100)
par(mar=c(5,5,1,1))
bh<-barplot(t(meanherbs[,c(14,11,17,23)]),beside=T,names.arg=meanherbs$Livestock.density
        ,xlab='Livestock density',las=1,ylim=c(0,350)
        ,ylab=expression('Grazer density (metabolic biomass) kg km   '^-2),
        col=c('white',grey(0.3),grey(0.5),(grey(0.8))),
                legend.text=c('Cattle','Zebra','Gazelle','Hartebeest')
        ,args.legend=list(x='top'))
arrows(bh,t(meanherbs[,c(14,11,17,23)])+t(seherbs[,c(14,11,17,23)]),bh,t(meanherbs[,c(14,11,17,23)])-t(seherbs[,c(14,11,17,23)]),code=3,length=0.05,angle=90)
dev.off()


dfComp<-droplevels(nsbiomass3_herbio[nsbiomass3_herbio$Treatment=="Control" | nsbiomass3_herbio$Treatment=="Exclosure",])
dfWide<- aggregate.data.frame(dfComp[,c(7:24,30:52)],list(Treatment=dfComp$Treatment,Livestock.density=dfComp$Livestock.density,Transect=dfComp$Transect),mean,na.rm=T)

dfDiff<-dfWide[dfWide$Treatment=='Exclosure',4:21]-dfWide[dfWide$Treatment=='Control',4:21]
dfDiff1<-cbind(dfDiff,dfWide[dfWide$Treatment=='Exclosure',c(2:3,22:44)])

dfDiff1$GrazersTLU<-dfDiff1$CattleTLU+dfDiff1$Burchells_ZebraTLU+dfDiff1$Grants_GazelleTLU+dfDiff1$Swaynes_HartebeestTLU
dfDiff1$GrazersMetBio<-dfDiff1$CattleMetBio+dfDiff1$Burchells_ZebraMetBio+dfDiff1$Grants_GazelleMetBio+dfDiff1$Swaynes_HartebeestMetBio
dfDiff1$TotalTLU<-dfDiff1$CattleTLU+dfDiff1$Burchells_ZebraTLU+dfDiff1$Grants_GazelleTLU+dfDiff1$Swaynes_HartebeestTLU+dfDiff1$Greater_KuduTLU
dfDiff1$TotalMetBio<-dfDiff1$CattleMetBio+dfDiff1$Burchells_ZebraMetBio+dfDiff1$Grants_GazelleMetBio+dfDiff1$Swaynes_HartebeestMetBio+dfDiff1$Greater_KuduMetBio


#cols<-c(grey(0),grey(0.3),grey(0.5),grey(0.8))
dfDiff1[,21:47]<-dfDiff1[,21:47]+0.0000000001

#Figure
{
  tiff("M:\\Africa\\VegetationDataSeason3\\Figs\\TreatmentDifference.tif",width=12,height=8,units="in",res=100)
  
  #  x11(12,12)
par(oma=c(3,3,0,0))
par(mfcol=c(4,2))
par(mar=c(3,3,1,2))
#cols<-c('blue','yellow','orange','brown')
cols<-c('white',grey(0.3),grey(0.5),(grey(0.8)))

with(dfDiff1,plot(GrazersMetBio,TotalSeason3,xlim=c(50,750),ylim=c(-70,300),las=1,type='n',main='Total',
                   xlab='',
                   ylab=''))

for(i in 1:nrow(dfDiff1)){
  with(dfDiff1[i,],floating.pie(GrazersMetBio,TotalSeason3,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
  #  with(dfDiff1[i],draw.circle(GrazersMetBio,TotalProductivityS2,radius=30,col=(Livestock.density)))
}
with(dfDiff1[dfDiff1$Livestock.density=='Low',],points(GrazersMetBio,TotalSeason3,pch=1,col=grey(0.8),cex=1.3))
with(dfDiff1[dfDiff1$Livestock.density=='Medium',],points(GrazersMetBio,TotalSeason3,pch=3,col=grey(0.5),cex=1.3))
with(dfDiff1[dfDiff1$Livestock.density=='High',],points(GrazersMetBio,TotalSeason3,pch=4,col=grey(0),cex=1.3))
#legend('bottomr',fill=cols,c('Cattle','Zebra','Gazelle','Hartebeest'),title='Metabolic biomass')
legend('bottomr',pch=c(1,3,4),col=c(grey(0.8),grey(0.5),grey(0.3)),c('Low','Medium','High'),title='Livestock density')

summary(lm(TotalSeason3~GrazersMetBio,data=dfDiff1))

par(mar=c(3,3,1,2))
with(dfDiff1,plot(GrazersMetBio,GrassNetBiomassSeason3,xlim=c(50,750),ylim=c(-20,280),las=1,type='n',ylab='',main='Grass',xlab=""))
for(i in 1:nrow(dfDiff1)){
  with(dfDiff1[i,],floating.pie(GrazersMetBio,GrassNetBiomassSeason3,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
  #  with(dfDiff1[i],draw.circle(GrazersMetBio,GrassNetReharvestBiomass2,radius=30,col=(Livestock.density)))
}
with(dfDiff1[dfDiff1$Livestock.density=='Low',],points(GrazersMetBio,GrassNetBiomassSeason3,pch=1,col=grey(0.8),cex=1.3))
with(dfDiff1[dfDiff1$Livestock.density=='Medium',],points(GrazersMetBio,GrassNetBiomassSeason3,pch=3,col=grey(0.5),cex=1.3))
with(dfDiff1[dfDiff1$Livestock.density=='High',],points(GrazersMetBio,GrassNetBiomassSeason3,pch=4,col=grey(0),cex=1.3))
#legend('bottomr',fill=cols,c('Cattle','Zebra','Gazelle','Hartebeest'),title='Metabolic biomass')
summary(lm(GrassNetBiomassSeason3~GrazersMetBio,data=dfDiff1))

par(mar=c(3,3,1,2))
with(dfDiff1,plot(GrazersMetBio,HerbClimb3,xlim=c(50,750),ylim=c(-70,70),las=1,type='n',ylab='',main='Herbs',xlab= ""))
for(i in 1:nrow(dfDiff1)){
  with(dfDiff1[i,],floating.pie(GrazersMetBio,HerbClimb3,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
}
with(dfDiff1[dfDiff1$Livestock.density=='Low',],points(GrazersMetBio,HerbClimb3,pch=1,col=grey(0.8)))
with(dfDiff1[dfDiff1$Livestock.density=='Medium',],points(GrazersMetBio,HerbClimb3,pch=3,col=grey(0.5)))
with(dfDiff1[dfDiff1$Livestock.density=='High',],points(GrazersMetBio,HerbClimb3,pch=4,col=grey(0)))
summary(lm(HerbClimb3~GrazersMetBio,data=dfDiff1))

par(mar=c(3,3,1,2))
with(dfDiff1,plot(GrazersMetBio,DwarfShrubNetBiomassSeason3,xlim=c(50,750),ylim=c(-60,80),las=1,type='n',
                  ylab='',main='Shrubs',
                  xlab=''))
for(i in 1:nrow(dfDiff1)){
  with(dfDiff1[i,],floating.pie(GrazersMetBio,DwarfShrubNetBiomassSeason3,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
  #  with(dfDiff1[i],draw.circle(GrazersMetBio,GrassNetReharvestBiomass2,radius=30,col=(Livestock.density)))
}
with(dfDiff1[dfDiff1$Livestock.density=='Low',],points(GrazersMetBio,DwarfShrubNetBiomassSeason3,pch=1,col=grey(0.8)))
with(dfDiff1[dfDiff1$Livestock.density=='Medium',],points(GrazersMetBio,DwarfShrubNetBiomassSeason3,pch=3,col=grey(0.5)))
with(dfDiff1[dfDiff1$Livestock.density=='High',],points(GrazersMetBio,DwarfShrubNetBiomassSeason3,pch=4,col=grey(0)))
summary(lm(DwarfShrubNetBiomassSeason3~GrazersMetBio,data=dfDiff1))
anova(lm(DwarfShrubNetBiomassSeason3~GrazersMetBio,data=dfDiff1))
abline(lm(DwarfShrubNetBiomassSeason3~GrazersMetBio,data=dfDiff1))

mtext(side=1,expression('Grazer density (metabolic biomass) kg km'^-2),line=1,outer=T)
mtext(side=2,expression("Grazing treatment reduction in biomass (Exclosure - Grazed) g m    "^-2),line=0,outer=T)


#Productivity
nsreharvest3$pairID<-factor(with(nsreharvest3,paste(Transect,Livestock.density)))
nsreharvest3$Zone<-nsreharvest3$pairID

levels(nsreharvest3$Zone)<-c('D2','A1','C2','D3','B2','C3','D4','A3','C4')

nsreharvest3_herbio<-merge(nsreharvest3,herbdens,by.x='Zone',by.y='Zone',all.x=T,all.y=F)

dfHComp<-droplevels(nsreharvest3_herbio[nsreharvest3_herbio$Treatment=="Control" | nsreharvest3_herbio$Treatment=="Exclosure",])
dfHWide<- aggregate.data.frame(dfHComp[,c(10:21,27:49)],list(Treatment=dfHComp$Treatment,Livestock.density=dfHComp$Livestock.density,Transect=dfHComp$Transect),mean,na.rm=T)

dfHDiff<-dfHWide[dfHWide$Treatment=='Exclosure',4:15]-dfHWide[dfHWide$Treatment=='Control',4:15]
dfHDiff1<-cbind(dfHDiff,dfHWide[dfHWide$Treatment=='Exclosure',c(2:3,16:38)])

dfHDiff1$GrazersTLU<-dfHDiff1$CattleTLU+dfHDiff1$Burchells_ZebraTLU+dfHDiff1$Grants_GazelleTLU+dfHDiff1$Swaynes_HartebeestTLU
dfHDiff1$GrazersMetBio<-dfHDiff1$CattleMetBio+dfHDiff1$Burchells_ZebraMetBio+dfHDiff1$Grants_GazelleMetBio+dfHDiff1$Swaynes_HartebeestMetBio
dfHDiff1$TotalTLU<-dfHDiff1$CattleTLU+dfHDiff1$Burchells_ZebraTLU+dfHDiff1$Grants_GazelleTLU+dfHDiff1$Swaynes_HartebeestTLU+dfHDiff1$Greater_KuduTLU
dfHDiff1$TotalMetBio<-dfHDiff1$CattleMetBio+dfHDiff1$Burchells_ZebraMetBio+dfHDiff1$Grants_GazelleMetBio+dfHDiff1$Swaynes_HartebeestMetBio+dfHDiff1$Greater_KuduMetBio



#x11(7,12)
#par(mfrow=c(4,1))
#with(dfHDiff1,plot(Cattle,TotalProductivityS2))
#cols<-c(grey(0),grey(0.3),grey(0.5),grey(0.8))
#cols<-c('blue','yellow','orange','brown')
dfHDiff1[,15:41]<-dfHDiff1[,15:41]+0.0000000001

par(mar=c(3,4,1,1))
with(dfHDiff1,plot(GrazersMetBio,TotalProductivityS2,xlim=c(50,750),ylim=c(-30,230),las=1, type='n',main='Total',
                   xlab='',
                   ylab=''))

for(i in 1:nrow(dfHDiff1)){
  with(dfHDiff1[i,],floating.pie(GrazersMetBio,TotalProductivityS2,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
#  with(dfHDiff1[i],draw.circle(GrazersMetBio,TotalProductivityS2,radius=30,col=(Livestock.density)))
}
with(dfHDiff1[dfHDiff1$Livestock.density=='Low',],points(GrazersMetBio,TotalProductivityS2,pch=1,col=grey(0.8),cex=1.3))
with(dfHDiff1[dfHDiff1$Livestock.density=='Medium',],points(GrazersMetBio,TotalProductivityS2,pch=3,col=grey(0.5),cex=1.3))
with(dfHDiff1[dfHDiff1$Livestock.density=='High',],points(GrazersMetBio,TotalProductivityS2,pch=4,col=grey(0),cex=1.3))
legend('bottomr',fill=cols,c('Cattle','Zebra','Gazelle','Hartebeest'),title='Metabolic biomass')
summary(lm(TotalProductivityS2~GrazersMetBio,data=dfHDiff1))


par(mar=c(3,4,1,1))
with(dfHDiff1,plot(GrazersMetBio,GrassNetReharvestBiomass2,xlim=c(50,750),ylim=c(-40,210)
                   ,las=1,type='n',ylab='',main='Grass'))
for(i in 1:nrow(dfHDiff1)){
  with(dfHDiff1[i,],floating.pie(GrazersMetBio,GrassNetReharvestBiomass2,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
  #  with(dfHDiff1[i],draw.circle(GrazersMetBio,GrassNetReharvestBiomass2,radius=30,col=(Livestock.density)))
}
with(dfHDiff1[dfHDiff1$Livestock.density=='Low',],points(GrazersMetBio,GrassNetReharvestBiomass2,pch=1,col=grey(0.8),cex=1.3))
with(dfHDiff1[dfHDiff1$Livestock.density=='Medium',],points(GrazersMetBio,GrassNetReharvestBiomass2,pch=3,col=grey(0.5),cex=1.3))
with(dfHDiff1[dfHDiff1$Livestock.density=='High',],points(GrazersMetBio,GrassNetReharvestBiomass2,pch=4,col=grey(0),cex=1.3))
#legend('bottomr',fill=cols,c('Cattle','Zebra','Gazelle','Hartebeest'),title='Metabolic biomass')
summary(lm(GrassNetReharvestBiomass2~GrazersMetBio,data=dfHDiff1))
anova(lm(GrassNetReharvestBiomass2~GrazersMetBio,data=dfHDiff1))
abline(lm(GrassNetReharvestBiomass2~GrazersMetBio,data=dfHDiff1))


par(mar=c(3,4,1,1))
with(dfHDiff1,plot(GrazersMetBio,HerbClimbReharv2,xlim=c(50,750),ylim=c(-40,70),las=1,type='n'
                   ,ylab='',main='Herbs'))
for(i in 1:nrow(dfHDiff1)){
  with(dfHDiff1[i,],floating.pie(GrazersMetBio,HerbClimbReharv2,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
  #  with(dfHDiff1[i],draw.circle(GrazersMetBio,GrassNetReharvestBiomass2,radius=30,col=(Livestock.density)))
}
with(dfHDiff1[dfHDiff1$Livestock.density=='Low',],points(GrazersMetBio,HerbClimbReharv2,pch=1,col=grey(0.8)))
with(dfHDiff1[dfHDiff1$Livestock.density=='Medium',],points(GrazersMetBio,HerbClimbReharv2,pch=3,col=grey(0.5)))
with(dfHDiff1[dfHDiff1$Livestock.density=='High',],points(GrazersMetBio,HerbClimbReharv2,pch=4,col=grey(0)))
summary(lm(HerbClimbReharv2~GrazersMetBio,data=dfHDiff1))

par(mar=c(3,4,1,1))
with(dfHDiff1,plot(GrazersMetBio,DwarfShrubNetReharvestBiomass2,xlim=c(50,750),ylim=c(-40,70),
                   las=1,type='n',ylab='',main='Shrubs'))
for(i in 1:nrow(dfHDiff1)){
  with(dfHDiff1[i,],floating.pie(GrazersMetBio,DwarfShrubNetReharvestBiomass2,x=c(CattleMetBio,Burchells_ZebraMetBio,Grants_GazelleMetBio,Swaynes_HartebeestMetBio)
                                 ,radius=30,col=cols,edges=500))
  #  with(dfHDiff1[i],draw.circle(GrazersMetBio,GrassNetReharvestBiomass2,radius=30,col=(Livestock.density)))
}
par(mar=c(3,4,1,1))
with(dfHDiff1[dfHDiff1$Livestock.density=='Low',],points(GrazersMetBio,DwarfShrubNetReharvestBiomass2,pch=1,col=grey(0.8)))
with(dfHDiff1[dfHDiff1$Livestock.density=='Medium',],points(GrazersMetBio,DwarfShrubNetReharvestBiomass2,pch=3,col=grey(0.5)))
with(dfHDiff1[dfHDiff1$Livestock.density=='High',],points(GrazersMetBio,DwarfShrubNetReharvestBiomass2,pch=4,col=grey(0)))
summary(lm(DwarfShrubNetReharvestBiomass2~GrazersMetBio,data=dfHDiff1))

mtext(side=1,expression('Grazer density (metabolic biomass) kg km'^-2),line=1,outer=T)
mtext(side=2,expression("Grazing treatment reduction in productivity (Exclosure - Grazed) g m    "^ -2),line=-47,outer=T)
dev.off()
}









#Species richness
nsspprich0<-read.csv("M:\\Africa\\VegetationDataSeason3\\SppRich_1.csv",header=T)
nsspprich1<-nsspprich0[nsspprich0$Season_Year=='Season III',]
nsspprich1$pairID<-factor(with(nsspprich1,paste(Transect,Livestock.density)))
nsspprich1$Zone<-nsspprich1$pairID

levels(nsspprich1$Zone)<-c('D2','A1','C2','D3','B2','C3','D4','A3','C4')
nsspprich_herbio<-merge(nsspprich1,herbdens,by.x='Zone',by.y='Zone',all.x=T,all.y=F)

dfSpComp<-droplevels(nsspprich_herbio[nsspprich_herbio$Treatment=="Control" | nsspprich_herbio$Treatment=="Exclosure",])
dfSpWide<- aggregate.data.frame(dfSpComp[,c(9:80,85:107)],list(Treatment=dfSpComp$Treatment,Livestock.density=dfSpComp$Livestock.density,Transect=dfSpComp$Transect),mean,na.rm=T)

dfSpDiff<-dfSpWide[dfSpWide$Treatment=='Exclosure',4:75]-dfSpWide[dfSpWide$Treatment=='Control',4:75]
dfSpDiff1<-cbind(dfSpDiff,dfSpWide[dfSpWide$Treatment=='Exclosure',c(2:3,76:98)])

dfSpDiff1$GrazersTLU<-dfSpDiff1$CattleTLU+dfSpDiff1$Burchells_ZebraTLU+dfSpDiff1$Grants_GazelleTLU+dfSpDiff1$Swaynes_HartebeestTLU
dfSpDiff1$GrazersMetBio<-dfSpDiff1$CattleMetBio+dfSpDiff1$Burchells_ZebraMetBio+dfSpDiff1$Grants_GazelleMetBio+dfSpDiff1$Swaynes_HartebeestMetBio
dfSpDiff1$TotalTLU<-dfSpDiff1$CattleTLU+dfSpDiff1$Burchells_ZebraTLU+dfSpDiff1$Grants_GazelleTLU+dfSpDiff1$Swaynes_HartebeestTLU+dfSpDiff1$Greater_KuduTLU
dfSpDiff1$TotalMetBio<-dfSpDiff1$CattleMetBio+dfSpDiff1$Burchells_ZebraMetBio+dfSpDiff1$Grants_GazelleMetBio+dfSpDiff1$Swaynes_HartebeestMetBio+dfSpDiff1$Greater_KuduMetBio

with(dfSpDiff1,plot(Cattle,Total_spprich))
with(dfSpDiff1,plot(GrazersMetBio,Total_spprich))
summary(lm(Total_spprich~GrazersMetBio,data=dfSpDiff1))
summary(lm(Total_spprich~Cattle,data=dfSpDiff1))

require(vegan)
dfSpOrd<-dfSpComp
dfSpOrd[dfSpOrd$Treatment=='Exclosure',85:107]<-0

mds1<-metaMDS(dfSpOrd[,9:73])
plot(mds1)
points(scores(mds1)[dfSpOrd$Treatment=='Exclosure',],col='red',pch=16)
points(scores(mds1)[dfSpOrd$Treatment=='Control',],col=grey(0.5),pch=16)
ef1<-envfit(mds1,dfSpOrd[,c(94,97,100,103,106)])
plot(ef1)


pr1<-prcomp(dfSpDiff1[,1:65])
ef2<-envfit(pr1,dfSpDiff1[,c(84,87,90,93,96)])
plot(ef2)
##################
spdata<-nsspprich[,8:72]
mds1<-metaMDS(spdata)
plot(mds1,type='n')
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium',],pch=16,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium',],pch=1,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High',],pch=1,col=grey(0))
ordihull(mds1,paste(nsspprich$Treatment,nsspprich$Livestock.density,nsspprich$Date) )

par(mfrow=c(1,3))
plot(mds1,type='n',main="Low herbivore density")
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season II',],pch=16,col=grey(0.5))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season II',],pch=1,col=grey(0.5))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season III',],pch=16,col=grey(0.3))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season III',],pch=1,col=grey(0.3))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season IV',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season IV',],pch=1,col=grey(0))


plot(mds1,type='n',main="Medium herbivore density")
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season II',],pch=16,col=grey(0.5))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season II',],pch=1,col=grey(0.5))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season III',],pch=16,col=grey(0.3))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season III',],pch=1,col=grey(0.3))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season IV',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season IV',],pch=1,col=grey(0))

plot(mds1,type='n',main="High herbivore density")
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season II',],pch=16,col=grey(0.5))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season II',],pch=1,col=grey(0.5))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season III',],pch=16,col=grey(0.3))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season III',],pch=1,col=grey(0.3))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season IV',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season IV',],pch=1,col=grey(0))


par(mfrow=c(1,3))
plot(mds1,type='n',main="Season I")
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season I',],pch=16,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season I',],pch=1,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season I',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season I',],pch=1,col=grey(0))

plot(mds1,type='n',main="Season II")
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season II',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season II',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season II',],pch=16,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season II',],pch=1,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season II',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season II',],pch=1,col=grey(0))

plot(mds1,type='n',main="Season III")
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season III',],pch=16,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Low'& nsspprich$Season=='Season III',],pch=1,col=grey(0.8))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season III',],pch=16,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='Medium'& nsspprich$Season=='Season III',],pch=1,col=grey(0.4))
points(scores(mds1)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High'& nsspprich$Season=='Season III',],pch=16,col=grey(0))
points(scores(mds1)[nsspprich$Treatment=='Control'& nsspprich$Livestock.density=='High'& nsspprich$Season=='Season III',],pch=1,col=grey(0))



mds2<-metaMDS(spdata,k=3)
require(scatterplot3d)

x11(18,12)
par(mfrow=c(1,3))
s3d<-scatterplot3d(scores(mds2),type='n',main='Low livestock density')
mtext(side=3,adj=0,"a")
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season II',],pch=16,col=grey(0.5))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season III',],pch=16,col=grey(0.3))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season IV',],pch=16,col=grey(0))

s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season II',],pch=1,col=grey(0.5))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season III',],pch=1,col=grey(0.3))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season IV',],pch=1,col=grey(0))

s3d<-scatterplot3d(scores(mds2),type='n',main='Medium livestock density')
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season II',],pch=16,col=grey(0.5))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season III',],pch=16,col=grey(0.3))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season IV',],pch=16,col=grey(0))

s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season II',],pch=1,col=grey(0.5))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season III',],pch=1,col=grey(0.3))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season IV',],pch=1,col=grey(0))

s3d<-scatterplot3d(scores(mds2),type='n',main='High livestock density')
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season II',],pch=16,col=grey(0.5))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season III',],pch=16,col=grey(0.3))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season IV',],pch=16,col=grey(0))

s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season II',],pch=1,col=grey(0.5))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season III',],pch=1,col=grey(0.3))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season IV',],pch=1,col=grey(0))
legend('topl',pch=c(1,1,1,1,16,16,16,16),col=c(grey(0.8),grey(0.5),grey(0.3),grey(0)),c("Grazed - Season 1","Grazed - Season 2", "Grazed - Season 3", "Grazed - Season 4","Exclosed - Season 1","Exclosed - Season 2", "Exclosed - Season 3", "Exclosed - Season 4"))


x11(18,12)
par(mfrow=c(1,3))
s3d<-scatterplot3d(scores(mds2),type='n',main='Season 1')
mtext(side=3,adj=0,"a")
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season I',],pch=1,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season I',],pch=1,col=grey(0.4))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season I',],pch=1,col=grey(0))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season I',],pch=16,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season I',],pch=16,col=grey(0.4))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season I',],pch=16,col=grey(0))

s3d<-scatterplot3d(scores(mds2),type='n',main='Season 2')
mtext(side=3,adj=0,"b")
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season II',],pch=1,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season II',],pch=1,col=grey(0.4))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season II',],pch=1,col=grey(0))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season II',],pch=16,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season II',],pch=16,col=grey(0.4))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season II',],pch=16,col=grey(0))

s3d<-scatterplot3d(scores(mds2),type='n',main='Season 3')
mtext(side=3,adj=0,"c")
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season III',],pch=1,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season III',],pch=1,col=grey(0.4))
s3d$points(scores(mds2)[nsspprich$Treatment=='Exclosure' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season III',],pch=1,col=grey(0))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Low' & nsspprich$Season=='Season III',],pch=16,col=grey(0.8))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='Medium' & nsspprich$Season=='Season III',],pch=16,col=grey(0.4))
s3d$points(scores(mds2)[nsspprich$Treatment=='Control' & nsspprich$Livestock.density=='High' & nsspprich$Season=='Season III',],pch=16,col=grey(0))

legend('topl',pch=c(1,1,1,16,16,16),col=c(grey(0.8),grey(0.4),grey(0)),c("Exclosed - Low livestock","Exclosed - Medium livestock","Exclosed - High livestock","Grazed - Low livestock","Grazed - Medium livestock","Grazed - High livestock"),ncol=2,cex=0.6,bg='white')


permutest(cca(nsspprich[,8:72]),)
betadisper(vegdist(nsspprich[,8:72]),paste(nsspprich$Treatment,nsspprich$Livestock.density))
