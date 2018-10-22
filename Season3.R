#####################################################
# Nech Sar - Biomass regrowth and pastoralism 
rm(list=ls())
library(MASS)
library(vegan)
library(ggplot2)
library(plyr)
library(rgeos)
library(sp)
library(raster)
library(rgdal)

#####################################################
#### Exclosure biomass ####
#####################################################
nsbiomass3<-read.table("BiomassSeason3a.txt",header=T,sep="\t")

# Structure of data
names(nsbiomass3)
head(nsbiomass3)
tail(nsbiomass3)
str(nsbiomass3)

# Standard error function
sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# Factors 
nsbiomass3$Livestock.density<-factor(nsbiomass3$Livestock.density,levels=c("Low","Medium","High"),ordered=T)
nsbiomass3$HerbClimb1<-nsbiomass3$HerbNetBiomassSeason1+nsbiomass3$ClimberNetBiomassSeason1
nsbiomass3$HerbClimb2<-nsbiomass3$HerbNetBiomassSeason2+nsbiomass3$ClimberNetBiomassSeason2
nsbiomass3$HerbClimb3<-nsbiomass3$HerbNetBiomassSeason3+nsbiomass3$ClimberNetBiomassSeason3
nsbiomass3$plotid<-paste(nsbiomass3$Treatment,nsbiomass3$Livestock.density)

#Convert to gm^-2
#nsbiomass3[,6:23]<-nsbiomass3[,6:23]*4 # Original - changed order?
nsbiomass3[,c(6:20,23:25)]<-nsbiomass3[,c(6:20,23:25)]*4

levels(nsbiomass3$Treatment)
nsbiomass3orig<- droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",])
with(nsbiomass3,tapply(TotalSeason3,list(Treatment,Livestock.density),mean))

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
     cex=(nsherb3$Burchells_Zebra/40), 
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

##### Combine dung and productivity exclosure xy coordinates ####

Ex_location<-read.csv(file="PlotNames.csv", sep=",",header=TRUE)
names(Ex_location)
levels(Ex_location$OtherName)
colnames(Ex_location)[7]<-"Trt.name"

#### Join regrowth and dung data ####
# Regrowth biomass
nsreharvest3<-read.table("ProductivitySeason3.txt",header=T,sep="\t")
names(nsreharvest3)

# Check level names
levels(nsreharvest3$Trt.name)
levels(Ex_location$Trt.name)

dim(nsreharvest3)
myvars <- c("Trt.name", "X", "Y")
Ex_location1<-Ex_location[myvars]
nsreharvest3loc<-left_join(nsreharvest3,Ex_location1, by=c("Trt.name"),drop=F)
dim(nsreharvest3loc)

# Spatial join - herbivore observations and regrowth
utmproj<-"+proj=utm +north +zone=37 +init=EPSG:32637" # central-ethiopia-37n

wp222loc<- cbind(nsreharvest3loc$X,nsreharvest3loc$Y)# get geographical location
Biosp<-SpatialPointsDataFrame(wp222loc,nsreharvest3loc, proj4string=CRS(utmproj),match.ID = TRUE)

Herbloc<-cbind(nsherb3$Xcent,nsherb3$Ycent)# get geographical location
Herbsp<-SpatialPointsDataFrame(Herbloc,nsherb3, proj4string=CRS(utmproj),match.ID = TRUE)


# Nearest number code
nsreharvest3loc$nearest_in_set2 <- apply(gDistance(Herbsp,Biosp, byid=TRUE), 1, which.min)
nsherb3$nearest_in_set2<-seq(1:17)

dim(nsherb3)
# Join herbivore biomass
MyHerb<-c("nearest_in_set2","Date","Burchells_ZebraMetBio","CattleMetBio","Grants_GazelleMetBio","Greater_KuduMetBio","Swaynes_HartebeestMetBio")
nsherb3sub<-nsherb3[MyHerb]

# Regrowth and herbivore metabolic biomass
levels(nsreharvest3loc$Reharvest.date) #"21.11.2013" "21.6.2013"  "25.11.2012"
#25.3.2012 26.10.2012 8.2.2013 9.2.2013
levels(nsherb3sub$Date)<-c("21.6.2013","21.6.2013" ,"21.11.2013", "21.11.2013")
nsherb3sub$Reharvest.date<-nsherb3sub$Date
nsherb3sub<-aggregate(.~Reharvest.date+nearest_in_set2,nsherb3sub,mean)

nsreharvest3locHerb<-left_join(nsreharvest3loc,nsherb3sub, by=c("Reharvest.date","nearest_in_set2"),drop=F)
nsReharvest<-nsreharvest3locHerb

nsreharvest3locHerb$Reharvest.date
nsReharvest$Reharvest.date

nsherb3sub$Burchells_ZebraMetBio
nsreharvest3locHerb$Burchells_ZebraMetBio
#### Alternative - create krigging maps of herbivores ####
#### Create krigging maps of herbivore metabolic biomass ####
#coordinates(nsherb3) <- ~ Xcent + Ycent
#nsherb3.grid<-nsreharvest3locHerb
#coordinates(nsherb3.grid) <- ~ X + Y

# Cattle variogram
#cattle.vgm <- variogram(CattleMetBio~1, nsherb3)
#cattlef = function(x) attr(cattle.fit <<- fit.variogram(cattle.vgm, vgm(,"Mat",nugget=NA,kappa=x)),"SSErr")
#optimize(cattlef, c(0.1, 5))
#plot(cattle.vgm, cattle.fit)
#cattle.kriged <- krige(CattleMetBio ~ 1, nsherb3, nsherb3.grid, model=cattle.fit)

#cattle.kriged %>% as.data.frame %>%
#  ggplot(aes(x=X, y=Y)) + geom_tile(aes(fill=var1.pred), height=500,width=500) + coord_equal() +
#  scale_fill_gradient(low = "yellow", high="red") +
#  theme_bw()
# Largerly under-estimates biomass values of cattle - 
# Max value ~250...due to plateau in semi-variogram

#nsreharvest3locHerb$cat.krig<-cattle.kriged$var1.pred
#names(nsreharvest3locHerb)
#plot(nsreharvest3locHerb$cat.krig,nsreharvest3locHerb$CattleMetBio)
# No much difference here # Krig = lower estimate of cattle 

########################################################################
#### Productivity - regrowth ####
########################################################################

# What is the spatial distribution of regrowth measurements
# Regrowth biomass
#nsreharvest3<-read.table("ProductivitySeason3.txt",header=T,sep="\t")
#names(nsreharvest3)

# Set up as date
nsReharvest$Harvest.date<-as.Date(nsReharvest$Harvest.date,"%d.%m.%Y")
nsReharvest$Reharvest.date<-as.Date(nsReharvest$Reharvest.date,"%d.%m.%Y")

# Only need exclosure and control
nsReharvestb<-droplevels(nsReharvest[nsReharvest$Treatment=="Control" | nsReharvest$Treatment=="Exclosure",])

# Livestock density and regrowth
nsReharvestb<- droplevels(nsReharvestb[nsReharvestb$Harvest!="original",])
nsReharvestb$harvest_code<-as.factor(with(nsReharvestb, paste(Livestock.density,Treatment, sep="-")))
levels(nsReharvestb$harvest_code) # 18

nsReharvestavg<-aggregate(TotalBiomass1~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb, mean)
nsReharvestsem<-aggregate(TotalBiomass1~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb,sem)
nsReharvestavg<-cbind(nsReharvestavg,nsReharvestsem[6])
colnames(nsReharvestavg)[7]<-"se"

nsReharvestavg0<-aggregate(TotalBiomass0~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb, mean)
nsReharvestsem0<-aggregate(TotalBiomass0~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb,sem)
nsReharvestavg0<-cbind(nsReharvestavg0,nsReharvestsem0[6])
colnames(nsReharvestavg0)[7]<-"se"
colnames(nsReharvestavg0)[6]<-"TotalBiomass1"

# Total biomassregrowth
Regrow<-ggplot(nsReharvestavg, aes(x=Reharvest.date, y=TotalBiomass1,linetype=harvest_code,shape=Treatment,colour=Harvest,fill=Harvest)) 
Regrow<-Regrow+geom_errorbar(aes(x = Reharvest.date, ymin=TotalBiomass1-se,ymax=TotalBiomass1+se),width=.2,show.legend=F)
Regrow<-Regrow+geom_point(size=3)
Regrow<-Regrow+facet_wrap(~Treatment+Livestock.density, scale="fixed")
Regrow<-Regrow+geom_line(aes(linetype=harvest_code), show.legend = F) +theme_classic()
Regrow<-Regrow+geom_errorbar(data=nsReharvestavg0,aes(x = Reharvest.date, ymin=TotalBiomass1-se,ymax=TotalBiomass1+se),width=.2, colour="green4",show.legend=F)
Regrow<-Regrow+geom_point(data=nsReharvestavg0,size=3, colour="green4",show.legend=F)
Regrow<-Regrow+ggtitle("Total Biomass")
Regrow

# Grass biomass - Create a plot code + harvest - single or double code
nsreharvest3Gavg<-aggregate(GrassNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb, mean)
nsreharvest3Gsem<-aggregate(GrassNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb,sem)
nsReharvestGavg<-cbind(nsreharvest3Gavg,nsreharvest3Gsem[6])
colnames(nsReharvestGavg)[7]<-"se"

# Grass biomass regrowth through time
RegrowG<-ggplot(nsReharvestGavg, aes(x=Reharvest.date, y=GrassNetReharvestBiomass1,shape=Treatment,colour=Harvest,fill=Harvest)) 
RegrowG<-RegrowG+geom_errorbar(aes(x = Reharvest.date, ymin=GrassNetReharvestBiomass1-se,ymax=GrassNetReharvestBiomass1+se),width=.2,show.legend=F)
RegrowG<-RegrowG+geom_point(size=3)
RegrowG<-RegrowG+facet_wrap(~Livestock.density, scale="fixed")
RegrowG<-RegrowG+geom_line(aes(linetype=harvest_code), show.legend = F) +theme_classic()
RegrowG<-RegrowG+ggtitle("Grass Biomass")
RegrowG

# Woody biomass - Create a plot code + harvest - single or double code
nsreharvest3Wavg<-aggregate(DwarfShrubNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb, mean)
nsreharvest3Wsem<-aggregate(DwarfShrubNetReharvestBiomass1~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestb,sem)
nsReharvestWavg<-cbind(nsreharvest3Wavg,nsreharvest3Wsem[6])
colnames(nsReharvestWavg)[7]<-"se"

# Woody biomass regrowth through time
RegrowW<-ggplot(nsReharvestGavg, aes(x=Reharvest.date, y=DwarfShrubNetReharvestBiomass1,shape=Treatment,colour=Harvest,fill=Harvest)) 
RegrowW<-RegrowW+geom_errorbar(aes(x = Reharvest.date, ymin=DwarfShrubNetReharvestBiomass1-se,ymax=DwarfShrubNetReharvestBiomass1+se),width=.2,show.legend=F)
RegrowW<-RegrowW+geom_point(size=3)
RegrowW<-RegrowW+facet_wrap(~Livestock.density, scale="fixed")
RegrowW<-RegrowW+geom_line(aes(linetype=harvest_code), show.legend = F) +theme_classic()
RegrowW<-RegrowW+ggtitle("Woody Biomass")
RegrowW


#### Percent difference Original biomass ####
nsReharvestbO<-droplevels(nsReharvestb[is.na(nsReharvestb$Harvest.date),])
nsReharvestbH1_2<-droplevels(nsReharvestb[!is.na(nsReharvestb$Harvest.date),])
#nsReharvestbH1_2$TotalBiomass0<-nsReharvestbH1_2[nsReharvestbH1_2$TotalBiomass0<0.1,]<-.1

# Percent of biomass at each time interval - treatment
nsReharvestbH1_2$Tot.Per.diff<-(nsReharvestbH1_2$TotalBiomass1-nsReharvestbH1_2$TotalBiomass0)/sd(nsReharvestbH1_2$TotalBiomass0)

#nsReharvestbH1_2$G.Per.diff<-((nsReharvestbH1_2$GrassNetReharvestBiomass1-nsReharvestbO$GrassNetReharvestBiomass1)/nsReharvestbO$GrassNetReharvestBiomass1)*100
#nsReharvestbH1_2$W.Per.diff<-((nsReharvestbH1_2$DwarfShrubNetReharvestBiomass1-nsReharvestbO$DwarfShrubNetReharvestBiomass1)/nsReharvestbO$DwarfShrubNetReharvestBiomass1)*100

# Percent diff
nsReharvestavg<-aggregate(Tot.Per.diff~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestbH1_2, mean)
nsReharvestsem<-aggregate(Tot.Per.diff~Harvest+Reharvest.date+Treatment+Livestock.density+harvest_code,nsReharvestbH1_2,sem)
nsReharvestavg<-cbind(nsReharvestavg,nsReharvestsem[6])
colnames(nsReharvestavg)[7]<-"se"

# Total biomassregrowth
Regrow<-ggplot(nsReharvestavg, aes(x=Reharvest.date, y=Tot.Per.diff,linetype=harvest_code,shape=Treatment,colour=Harvest,fill=Harvest)) 
Regrow<-Regrow+facet_wrap(~Livestock.density, scale="fixed")
Regrow<-Regrow+geom_line(aes(linetype=harvest_code), show.legend = F) +theme_classic()
Regrow<-Regrow+geom_point(size=4)
Regrow<-Regrow+ggtitle("Total Biomass")
Regrow

###############################################################################
#### Mixed linear model (nlme) - auto-correlation? ####

# Only one harvest per time
nsReharvestbH1_2$Plot.ID<-as.numeric(nsReharvestbH1_2$Plot.name)
nsReharvestbH1_2$YrMonth<-format(as.Date(nsReharvestbH1_2$Reharvest.date), "%Y-%m")
nsReharvestH1<- droplevels(nsReharvestbH1_2[nsReharvestbH1_2$Harvest=="double",])
#nsReharvestH1<-nsReharvestbH1_2

nsReharvestH1$fBlock<-as.factor(nsReharvestH1$Block)
nsReharvestH1$fPlot.ID<-as.factor(nsReharvestH1$Plot.ID)

# Auto correlation structure - unneccessary
#cs1AR1 <- corAR1(0.2, form = ~Reharvest.date|fBlock/fPlot.ID)
#cs1AR1. <- Initialize(cs1AR1, data =nsReharvestH1) 
#corMatrix(cs1AR1.)
nsReharvestH1T<- nsReharvestH1[is.finite(nsReharvestH1$Tot.Per.diff),]

Tot1<-lme(Tot.Per.diff~Livestock.density+Treatment+
            Treatment:Livestock.density:Harvest.date,
         # Livestock.density:Harvest.date+
         # Livestock.density:Harvest+
         #Livestock.density:Harvest.date:Treatment, #TotalBiomass0,
          random= ~ 1|fBlock, method="ML",data=nsReharvestH1T)
          #correlation=corAR1(0.2, form=~Harvest.date|fBlock/fPlot.ID),data=nsReharvestH1)
summary(Tot1)
anova(Tot1)
AIC(Tot1) #2610.436
plot(ACF(Tot1),alpha=0.05) # Nothing systematic # But inital break of 1.0?

plot(Tot1) #OK
drop1(Tot1,test="Chisq") # 0.05717 . marginal
#Livestock.density:Treatment  2 727.16 5.6331  0.05981. marginal

# Explore different grazers

Tot1<-lme(Tot.Per.diff~Treatment+Burchells_ZebraMetBio,
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

ZebraPerS<-aggregate(TotalBiomass1~ Burchells_ZebraMetBio+Reharvest.date+Treatment,nsReharvestH1,sem)
CattlePerS<-aggregate(TotalBiomass1~ CattleMetBio+Reharvest.date+Treatment,nsReharvestH1,sem)
GrantPerS<-aggregate(TotalBiomass1~ Grants_GazelleMetBio+Reharvest.date+Treatment,nsReharvestH1,sem)

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
nsReharvestsem<-aggregate(TotalBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsReharvest,sem)

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
nsreharvest3Gsem<-aggregate(GrassNetReharvestBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsreharvest3locb,sem)

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
nsreharvest3Wsem<-aggregate(DwarfShrubNetReharvestBiomass1~Trt.name+Harvest.date+Harvest+Reharvest.date+Treatment+Livestock.density+X+Y,nsreharvest3locb,sem)

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
#### END OF REDITS ####
#####################################################


# Raw biomass for each harvest period
tiff("BiomassJune2017.tif",width=12,height=8,units="in",res=100)

par(mfrow=c(1,3))
meanbio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason1,list(Treatment,Livestock.density),mean,na.rm=T))
sembio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason1,list(Treatment,Livestock.density),sem))
b1<-barplot(meanbio,beside=T,main="November 2012",xlab="Livestock density",ylab="",ylim=c(0,1000),las=1)
title(ylab=expression("Biomass gm"^"-2"),line=2)
arrows(b1,meanbio+sembio,b1,meanbio-sembio,code=3,length=0.05,angle=90)

meanbio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason2,list(Treatment,Livestock.density),mean,na.rm=T))
sembio<- with(droplevels(nsbiomass3[nsbiomass3$Treatment=="Control" | nsbiomass3$Treatment=="Exclosure",]),tapply(TotalSeason2,list(Treatment,Livestock.density),sem))
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
