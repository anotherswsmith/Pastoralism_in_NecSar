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
names(nsSpp3)
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

# Plot centroids
CenPlot<-ggplot(vec.sp.df[order(vec.sp.df$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity,shape=Treatment))
CenPlot<-CenPlot+geom_point(aes(size=TotalBiomass1), stroke=1)
CenPlot<-CenPlot+geom_path(aes(group=harvest_code),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
CenPlot<-CenPlot +geom_text(aes(label=Season),hjust=0, vjust=-.5, show.legend = F)
CenPlot<-CenPlot +scale_colour_manual(values=c("grey60","gray60","grey20"))
CenPlot<-CenPlot +scale_fill_manual(values=c("grey80","black","grey40"))
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

# Mean distance moved by each plot - overall
nsSpp3$NMDS1<-mdsNS$points[,1]
nsSpp3$NMDS2<-mdsNS$points[,2]
nsSpp31<-droplevels(nsSpp3[nsSpp3$Season=="Short I",])
nsSpp32<-droplevels(nsSpp3[nsSpp3$Season=="Long",])
nsSpp33<-droplevels(nsSpp3[nsSpp3$Season=="Short II",])

nsSpp31sub<-nsSpp31[,c("Quadrats","Season","NMDS1","NMDS2")]
nsSpp32sub<-nsSpp32[,c("Quadrats","Season","NMDS1","NMDS2")]
nsSpp33sub<-nsSpp33[,c("Quadrats","Season","NMDS1","NMDS2")]
dim(nsSpp32sub)

# First season NMDS scores
cnt<-cbind(nsSpp31sub$NMDS1,nsSpp31sub$NMDS2)

#apply(m,1,function(x,cnt) {(sqrt((x[1] - cnt[1])^2+(x[2]-cnt[2])^2))},cnt)
library(ecodist)
euc.dist <- function(x1) sqrt(sum((x1 - cnt) ^ 2)) # bcdist= bray curtis
NMDSdist2<-apply(cbind(nsSpp32sub$NMDS1,nsSpp32sub$NMDS2), 1, bcdist)
NMDSdist3<-apply(cbind(nsSpp33sub$NMDS1,nsSpp33sub$NMDS2), 1, bcdist)

nsSpp32<-cbind(nsSpp32,NMDSdist2)
nsSpp33<-cbind(nsSpp33,NMDSdist3)
colnames(nsSpp32)[85]<-"Bcdist"
colnames(nsSpp33)[85]<-"Bcdist"
reharv23<-rbind(nsSpp32,nsSpp33)
names(reharv23)

# Large outlier - test with and without  - low livestock
#nsSpp3G23b<-nsSpp3G23[-c(90,180),]

# Transform bray curtis distance to positive values - so difference from zero..
abs(reharv23$Bcdist)
reharv23$Bcdist <- abs(reharv23$Bcdist)

# Mixed linear model - grass biomass and Bray curtis dist
library(glmmADMB)
names(reharv23)
min(reharv23$TotalBiomass1)
Bcdist<-glmmadmb(Bcdist~Livestockdensity+Treatment+Season+
                       # Livestockdensity:Treatment+#Season:Treatment+
                       # Livestockdensity:Season+Livestockdensity:Treatment:Season
                      +(1|fBlock), 
                      #admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                      family="gamma",data=reharv23)
summary(Bcdist)

E1 <- resid(Bcdist, type = "pearson") 
F1 <- fitted(Bcdist)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals", 
     xlim = c(min(F1), max(F1)))
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 2) #

# Drop any interactions?
drop1(Bcdist, test="Chi")
#Livestockdensity:Treatment:Season  2 728.08 0.574   0.7505
#Treatment:Season            1 732.77 3.212   0.0731 .
#Livestockdensity:Treatment  2 734.13  5.360 0.068563 . 

# Plot Euclidean distance vs livestock and treatment
sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
x<-aggregate(Bcdist~Livestockdensity+Treatment,reharv23,mean)
xsd<-aggregate(Bcdist~Livestockdensity+Treatment,reharv23,sem)
XEU<-cbind(x,xsd[,3])
colnames(XEU)[4]<-"sd"

p<-ggplot(XEU,aes(x=Livestockdensity,y=Bcdist,colour=Livestockdensity,shape=Treatment))
p<-p+geom_point(stat = "identity", size=3,position=position_dodge(width=.65))   
p<-p+geom_errorbar(aes(ymin=Bcdist-sd, ymax=Bcdist+sd),position=position_dodge(width=.65),stat = "identity", width=.2)
p<-p+ylab("Euclidean distance")
p<-p+theme_classic()
p

Ep<-ggplot(reharv23,aes(x=Eudist,y=TotalBiomass1,colour=Season,shape=Treatment))
Ep<-Ep+facet_wrap(~Livestockdensity)
Ep<-Ep+geom_point(stat = "identity", size=3,position=position_dodge(width=.65))   
Ep<-Ep+xlab("Euclidean distance")
Ep<-Ep+theme_classic()
Ep

# Species most asscoicated with specific livestock x treatment x season combination
nsSpp3.env<-nsSpp3[,c("fBlock","Livestockdensity","Season","Treatment")]
sim <- with(nsSpp3.env, simper(nsSpp3[,10:74],nsSpp3$Livestockdensity),permutations=999)
summary(sim,ordered = TRUE)

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

# Grass NMDS movement plot
CenPlotG<-ggplot(vec.sp.dfG[order(vec.sp.df$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity,shape=Treatment))
CenPlotG<-CenPlotG+geom_point(aes(size=GrassNetReharvestBiomass1), stroke=1)
CenPlotG<-CenPlotG+geom_path(aes(group=harvest_code),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
#CenPlotG<-CenPlotG+geom_text(data=spp.scrs,aes(label=Species),colour="light grey")
CenPlotG<-CenPlotG+geom_text(aes(label=Season),hjust=0, vjust=-.5, show.legend = F)
CenPlotG<-CenPlotG+scale_colour_manual(values=c("grey60","gray60","grey20"))
CenPlotG<-CenPlotG+scale_fill_manual(values=c("grey80","black","grey40"))
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

# Mean distance moved by each plot - grasses
nsSpp3G$NMDS1<-mdsNSG$points[,1]
nsSpp3G$NMDS2<-mdsNSG$points[,2]
nsSpp3G1<-droplevels(nsSpp3G[nsSpp3G$Season=="Short I",])
nsSpp3G2<-droplevels(nsSpp3G[nsSpp3G$Season=="Long",])
nsSpp3G3<-droplevels(nsSpp3G[nsSpp3G$Season=="Short II",])

nsSpp3G1sub<-nsSpp3G1[,c("Quadrats","Season","NMDS1","NMDS2")]
nsSpp3G2sub<-nsSpp3G2[,c("Quadrats","Season","NMDS1","NMDS2")]
nsSpp3G3sub<-nsSpp3G3[,c("Quadrats","Season","NMDS1","NMDS2")]
dim(nsSpp3G2sub)

#cnt = c(mean(m[,1]),mean(m[,2]))
cnt<-cbind(nsSpp3G1sub$NMDS1,nsSpp3G1sub$NMDS2)

# Bray curtis distance movement - Grasses
library(ecodist)
#apply(m,1,function(x,cnt) {(sqrt((x[1] - cnt[1])^2+(x[2]-cnt[2])^2))},cnt)
euc.dist <- function(x1) sqrt(sum((x1 - cnt) ^ 2))
test.NMDSdist2<-apply(cbind(nsSpp3G2sub$NMDS1,nsSpp3G2sub$NMDS2), 1, bcdist)
test.NMDSdist3<-apply(cbind(nsSpp3G3sub$NMDS1,nsSpp3G3sub$NMDS2), 1,bcdist)

#test.NMDSdist2<-as.data.frame(apply(cnt,1,function(x,cnt) {(sqrt((nsSpp3G2sub$NMDS1 - cnt[1])^2+(nsSpp3G2sub$NMDS2-cnt[2])^2))},cnt))
#test.NMDSdist3<-as.data.frame(apply(cnt,1,function(x,cnt) {(sqrt((nsSpp3G3sub$NMDS1 - cnt[1])^2+(nsSpp3G3sub$NMDS2-cnt[2])^2))},cnt))

nsSpp3G2<-cbind(nsSpp3G2,test.NMDSdist2)
nsSpp3G3<-cbind(nsSpp3G3,test.NMDSdist3)
colnames(nsSpp3G2)[38]<-"Bcdist"
colnames(nsSpp3G3)[38]<-"Bcdist"
nsSpp3G23<-rbind(nsSpp3G2,nsSpp3G3)
names(nsSpp3G23)

nsSpp3G1<-cbind(nsSpp3G1,test.NMDSdist2,test.NMDSdist3)
names(nsSpp3G1)
#colnames(nsSpp3G1)[38]<-"V1"
#colnames(nsSpp3G1)[39]<-"V2"
nsSpp3G1$Eudist<-rowSums(nsSpp3G1[,c(38:39)])
nsSpp3G1b<-nsSpp3G1[-90,]

# Large outlier - test with and without  - low livestock
nsSpp3G23b<-nsSpp3G23[-c(90,180),]

# Mixed linear model - grass biomass and Eudist
library(glmmADMB)
names(nsSpp3G23)
nsSpp3G23$fBlock<-as.factor(nsSpp3G23$Block)
nsSpp3G23$GrassNetReharvestBiomass1[nsSpp3G23$GrassNetReharvestBiomass1==0]<-1

# Bray curtis distance as positive
nsSpp3G23$Bcdist <- abs(nsSpp3G23$Bcdist)
nsSpp3G23b$Bcdist <- abs(nsSpp3G23b$Bcdist)

BcdistLMG<-glmmadmb(Bcdist~Livestockdensity+Treatment+Season+
                     Livestockdensity:Treatment+Season:Treatment+
                     Livestockdensity:Season+Livestockdensity:Treatment:Season
                   +(1|fBlock), 
                   #admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                   family="gamma",data=nsSpp3G23)
summary(BcdistLMG)
AIC(BcdistLMG) # Imrpoves AIC 754.038

#Checking assumptions
E1G <- resid(BcdistLMG, type = "pearson") 
F1G <- fitted(BcdistLMG)

par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1G, 
     y = E1G,
     xlab = "Fitted values",
     ylab = "Residuals", 
     xlim = c(min(F1), max(F1)))
abline(v = 0, lwd = 2, col = 2)
abline(h = 0, lty = 2, col = 2) # Better spread no transformation

drop1(BcdistLMG, test="Chi")
#Livestockdensity:Treatment:Season  2 758.06 8.022  0.01812 *

sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
x<-aggregate(Eudist~Livestockdensity+Treatment,nsSpp3G23b,mean)
xsd<-aggregate(Eudist~Livestockdensity+Treatment,nsSpp3G23b,sem)
XEU<-cbind(x,xsd[,3])
colnames(XEU)[4]<-"sd"

nsSpp3G23$Eudist
p<-ggplot(XEU,aes(x=Livestockdensity,y=Eudist,colour=Treatment))
#p<-p+ geom_boxplot()#outlier.shape = NA)
#p<-p+ scale_y_continuous(limits = c(4,11))
p<-p+geom_point(stat = "identity", size=3,position=position_dodge(width=.65))   
p<-p+geom_errorbar(aes(ymin=Eudist-sd, ymax=Eudist+sd),position=position_dodge(width=.65),stat = "identity", width=.2)
p<-p+theme_classic()
p

p1<-ggplot(nsSpp3G1b,aes(x=Eudist,y=GrassNetReharvestBiomass1,colour=Livestockdensity,shape=Treatment))
#p<-p+ geom_errorbar(aes(ymin=Eudist-sd,ymax=Eudist+sd), width=.2)
p1<-p1+geom_point(stat = "identity")   
p1

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