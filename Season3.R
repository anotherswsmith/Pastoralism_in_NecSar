#####################################################
# Nech Sar - Biomass regrowth and pastoralism 
rm(list=ls())
library(MASS)
library(vegan)
#####################################################

nsbiomass3<-read.table("BiomassSeason3a.txt",header=T,sep="")

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

nsdung3<-read.table("CountAverageFeb13.txt",header=T,sep=" ")

names(nsdung3)
nsdung3$Zone # Each herbivore transect corresponds to a unique zone unique zone? 

nsdung3$Cattle

# Spatial correlation
par(mfrow = c(1, 1), mar = c(4, 3, 3, 2))
par(pty = "s", mar = c(5,5,2,2), cex.lab = 1.5)       
plot(x =  nsdung3$Ycent,
     y =  nsdung3$Xcent,
     type = "p",
     cex=(nsdung3$Burchells_Zebra/10), #nsdung3$Cattle
     pch = 21,
     xlab = "X-coordinates",
     ylab = "Y-coordinates")

ordihull(nsdung3[, c("Ycent", "Xcent")],
         draw = "polygon",
         groups = nsdung3[, "Zone"],
         label = F,
         col = "red")     
# Zone = each unique herbivore transect





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
