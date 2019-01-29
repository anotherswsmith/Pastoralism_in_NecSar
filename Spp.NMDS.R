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
setwd("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar")
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
levels(nsreharvest3$Season)<-c("Short II" ,  "Long" , "Short I")
levels(nsSpp2$Season)<-c("Short I" ,  "Long" , "Short II")
nsSpp2$Trt.name<-as.factor(nsSpp2$Trt.name)
nsSpp2$Season<-as.factor(nsSpp2$Season)
nsreharvest3$Season <- factor(nsreharvest3$Season, levels(nsreharvest3$Season)[c(3,2,1)])
nsSpp2$plot_code<-as.factor(with(nsSpp2, paste(Transect,Block,Trt.name,Season,Replicate, sep="_")))
nsreharvest3$plot_code<-as.factor(with(nsreharvest3, paste(Transect,Block,Trt.name,Season,Replicate, sep="_")))

dim(nsSpp2)
dim(nsreharvest3)

nsreharvest3b<-nsreharvest3[c("Transect","Block","Treatment","Trt.name","Season","Replicate","plot_code","Livestock.density","GrassNetReharvestBiomass1","DwarfShrubNetReharvestBiomass1",
                              "HerbNetReharvestBiomass1","ClimberNetReharvestBiomass1","TotalBiomass1","rain.mm")]
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
names(nsSpp3G[10:28]) 

# Subset woody - dwarf shrubs
nsSpp3W<-nsSpp3[ , !(names(nsSpp3) %in% to.removeW)]
nsSpp3W<-subset(nsSpp3,select = names(nsSpp3) %ni% to.removeW)
names(nsSpp3W) # Subset is now just dwarf shrubs
names(nsSpp3W[10:19]) 

# Subset Herbs and climbers
nsSpp3C<-nsSpp3[ , !(names(nsSpp3) %in% to.removeC)]
nsSpp3C<-subset(nsSpp3,select = names(nsSpp3) %ni% to.removeC)
names(nsSpp3C) # Subset is now herbs and climber
names(nsSpp3C[10:38]) 

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
nsSpp3$harvest_code<-as.factor(with(nsSpp3, paste(Livestockdensity,Treatment, sep="-")))

# Envfit
NS.ev <- envfit(mdsNS~Livestockdensity+Treatment+ Season,data = nsSpp3, 
                  perm=999,strata=as.numeric(nsSpp3$Block))
NS.ev 
# Weak difference in season - biggest difference is livestock density - r2=0.26

NS.har <- envfit(mdsNS~ harvest_code,data = nsSpp3, 
                perm=999,strata=as.numeric(nsSpp3$Block))
NS.har

# PERMANOVA/ADONIS - Total 
# Distance matrix
vare.dis<-vegdist(nsSpp3[,10:74],"bray") 
vare.dis2<-as.matrix(vare.dis)

# Beta Dispersion
nsSpp3$harvest_code<-as.factor(with(nsSpp3, paste(Livestockdensity,Treatment, sep="-")))
modTLiv <- betadisper(vare.dis, nsSpp3$harvest_code,type = c("centroid"))
modTLivM <- betadisper(vare.dis, nsSpp3$harvest_code,type = c("median"))

# Dispersion sd
SE<- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
CentSD<-tapply(modTLiv$distances, nsSpp3$harvest_code, SE)



# PERMANOVA total biomass
nsSpp3$time_code<-as.factor(with(nsSpp3, paste(Transect,Block,Treatment,Replicate, sep="-")))
PermT<-adonis(vare.dis2 ~ Livestockdensity+Treatment+Season+#TotalBiomass1+
                Season:Treatment+Livestockdensity:Season+
                Livestockdensity:Treatment+#TotalBiomass1:Season+
                #TotalBiomass1:Livestockdensity+#TotalBiomass1:Treatment+
                Treatment:Livestockdensity:Season#+
              #Treatment:Livestockdensity:TotalBiomass1,
              strata=as.numeric(nsSpp3$Block),#/as.numeric(nsSpp3$time_code),
              method = "bray",perm=999, data=nsSpp3)
PermT$aov.tab
#PermT$aov.tab$R2
mean(modTLiv$centroids[,1])
modTLiv$centroids[1]
# Extract centroids 
#vec.sp.df<-as.data.frame(cbind(NS.har$factors$centroids*sqrt(NS.har$factors$r)))
vec.sp.df<-data.frame(cbind(modTLiv$centroids[,1],modTLiv$centroids[,2],CentSD))
vec.sp.df
colnames(vec.sp.df)[1]<-"NMDS1"
colnames(vec.sp.df)[2]<-"NMDS2"
colnames(vec.sp.df)[3]<-"CenSd"
data.frame(cbind(modTLiv$centroids[,1],modTLiv$centroids[,2]))
centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
vec.sp.df$Livestockdensity<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                               "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.df$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                       "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.df$Season<-rep(c("Long","Short I","Short II"))
vec.sp.df$harvest_code<-as.factor(with(vec.sp.df, paste(Livestockdensity,Treatment,Season, sep="-")))

TotBio<-aggregate(TotalBiomass1~harvest_code,nsSpp3,mean)
vec.sp.df<-merge(vec.sp.df, TotBio, by.x = "harvest_code")
vec.sp.df$TotalBiomass1<-as.numeric(vec.sp.df$TotalBiomass1)

# Reorder by season and Harvest code
vec.sp.df$Season<- factor(vec.sp.df$Season, levels = c("Short I","Long","Short II"))
vec.sp.df$Season2<-c(1,2,3)
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

# Importance of rainfall - use ordisurf
#https://oliviarata.wordpress.com/2014/07/17/ordinations-in-ggplot2-v2-ordisurf/
names(nsSpp3)
ordi <- ordisurf(mdsNS  ~ TotalBiomass1, data =  nsSpp3, plot = FALSE, scaling = 3,
                 method = "REML", select = TRUE)
summary(ordi) # Deviance explained 7.53%
anova(ordi)
#s(x1,x2) 3.474  9.000 2.014 0.000118

#plot(ordi)
ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
ordi.mite.na #looks ready for plotting!

colnames(ordi.mite.na)<-c("NMDS1","NMDS2","rain.mm")
#ordi.mite.na$Livestockdensity<-c("Low")
class(ordi.mite.na$rain.mm)
ordi.mite.na$TotalBiomass1<-ordi.mite.na$rain.mm

# Rename seasons, exclosure, - short hand
levels(vec.sp.df$Season)<-c(" ","L","SII")
vec.sp.df$Treatment<-as.factor(vec.sp.df$Treatment)
vec.sp.df$Livestockdensity<-as.factor(vec.sp.df$Livestockdensity)
levels(vec.sp.df$Treatment)<-c("Open","Exclosed")
vec.sp.df$Livestockdensity<- factor(vec.sp.df$Livestockdensity, levels(vec.sp.df$Livestockdensity)[c(2,3,1)])

sizeLegend<-expression(paste("Biomass regrowth (g ",m^-2,")"))

# Plot centroids
CenPlot<-ggplot(vec.sp.df[order(vec.sp.df$Season),],aes(x=NMDS1,y=NMDS2))
#CenPlot<-CenPlot+geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=Livestockdensity), size=1,show.legend=T)
CenPlot<-CenPlot +scale_x_continuous(limits = c(-0.2,.2), expand = c(0,0))
CenPlot<-CenPlot +scale_y_continuous(limits = c(-0.2,.25), expand = c(0,0))
CenPlot<-CenPlot+stat_contour(data = ordi.mite.na, aes(x = NMDS1, y = NMDS2, z = rain.mm, size=TotalBiomass1,alpha= (..level..)),colour = "green4", #colour = rev(..level..),
             binwidth = 2, show.legend = F)
CenPlot<-CenPlot+annotate(geom="text",x=-0.0, y=-0.16, label=expression(paste("70 g ",m^-2,"")), colour = "green4", size=4)
CenPlot<-CenPlot+annotate(geom="text",x=-0.0, y=0.047, label=expression(paste("     80 g ",m^-2,"")), colour = "green4", size=4)
CenPlot<-CenPlot+annotate(geom="text",x=-0.0, y=0.21, label=expression(paste("88 g ",m^-2,"")), colour = "green4", size=4)
CenPlot<-CenPlot+geom_errorbar(aes(colour=Livestockdensity,ymin=NMDS2-CenSd, ymax=NMDS2+CenSd),show.legend=F)
CenPlot<-CenPlot+geom_errorbarh(aes(colour=Livestockdensity,xmin = NMDS1-CenSd,xmax = NMDS1+CenSd),show.legend=F)
CenPlot<-CenPlot+geom_point(aes(shape=Treatment,size=TotalBiomass1,colour=Livestockdensity,fill=Livestockdensity), stroke=1)
CenPlot<-CenPlot+geom_path(aes(group=harvest_code,colour=Livestockdensity),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
CenPlot<-CenPlot +geom_text(aes(label=Season),hjust=0, vjust=-.95, show.legend = F)
CenPlot<-CenPlot +scale_colour_manual(values=c("grey50","grey80","black"))#"black","grey50","grey80"
CenPlot<-CenPlot +scale_fill_manual(values=c("grey50","grey80","black")) #"black","grey50","grey80"
CenPlot<-CenPlot +scale_radius(sizeLegend,range=c(1,8))
#CenPlot<-CenPlot +scale_alpha_manual(values=c(.75,.95))
CenPlot<-CenPlot +scale_shape_manual(values=c(21,22))
#CenPlot<-CenPlot +ggtitle("Total biomass")
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
        ,legend.background=element_blank()
        ,panel.spacing = unit(.1, "lines")
        ,legend.text=element_text(size=12)
        ,legend.title=element_text(size=12)
        ,legend.position = "right"
        ,legend.justification = "top"
        ,legend.direction="vertical"
        ,legend.key.width = unit(1.2,"cm"))
CenPlot

CenPlot<- CenPlot+guides(linetype=F, fill=F,#size=T,
                         colour = guide_legend("Livestock density",override.aes = list(shape=c(21), size=3.5,fill=c("black","grey50","grey80"),col=c("black","grey50","grey80"), stroke=1)),
                         shape = guide_legend("Treatment",override.aes = list(shape=c(21,22), size=3.5,fill=NA,col="grey30", stroke=1)))
      # size = guide_legend((expression(paste("Biomass ( kg ",m^-2,")"))),override.aes = list(shape=c(21),fill=c("grey30"),col=c("grey30"), stroke=1)))
CenPlot

# Export ggplot
ggsave("NMDS2.ordsurf.jpeg",width= 15, height = 12,units ="cm", bg ="transparent",
       dpi = 600, limitsize = TRUE)

# CHECKING NUMBERS FOR ordisurf line
library(directlabels)
#CenPlot2<-ggplot()
#CenPlot2<-CenPlot2 +scale_x_continuous(limits = c(-0.2,.2), expand = c(0,0))
#CenPlot2<-CenPlot2 +scale_y_continuous(limits = c(-0.2,.2), expand = c(0,0))
#CenPlot2<-CenPlot2+stat_contour(data = ordi.mite.na, aes(x = NMDS1, y = NMDS2, z = rain.mm,colour= (..level..)),colour = "green4", #colour = rev(..level..),
#                              binwidth = 2)
#direct.label(CenPlot2, list("visualcenter", colour='green4')) # dodgerblue1

# Predicted versus observed biomass based on changes in sepcies composition
str(ordi)
#y ~ s(x1, x2, k = 10, bs = "tp", fx = FALSE)
#A. Create a grid of covariate values
range(Crayfish2$TotalBiomass1)
MyData <- expand.grid(x1 = vec.sp.df$NMDS1,
                      x2= vec.sp.df$NMDS2,
                      k=10)
                      #TotalBiomass1= seq(min(nsSpp3$TotalBiomass1),max(nsSpp3$TotalBiomass1), length = 25))

predict(ordi,MyData,type="terms",se=TRUE)

#B. Predict Weight values for these artificial covariate values
P2 <- predict(ordi, newdata = MyData, se = TRUE)

#Put everything in MyData (for ggplot2)
MyData$mu   <- P2$fit
MyData$seup <- P2$fit + 1.96 * P2$se.fit
MyData$selo <- P2$fit - 1.96 * P2$se.fit
MyData
colnames(MyData)<-c("NMDS1","NMDS2","k","TotalBiomass1","seup","selo")
names(vec.sp.df)
PredGAM<-ggplot(vec.sp.df,aes(x=NMDS2, y=TotalBiomass1))+geom_point(size=2)
PredGAM<-PredGAM+geom_point(data=MyData,aes(x=NMDS2, y=TotalBiomass1), colour="red")
PredGAM

vec.sp.df$TotalBiomass1
PredGAMnmds2<-aggregate(TotalBiomass1~NMDS2,MyData,mean)
vec.sp.df$PredBio<-PredGAMnmds2$TotalBiomass1
vec.sp.df$PredDiff<-vec.sp.df$TotalBiomass1-PredGAMnmds2$TotalBiomass1

PredDiff2<-ggplot(vec.sp.df,aes(x=TotalBiomass1, y=PredBio, colour=Livestockdensity,shape=Treatment))
PredDiff2<-PredDiff2+geom_abline(slope=1, intercept=0, size =.95) # Reference line
PredDiff2<-PredDiff2+geom_point(size=2)
PredDiff2

# CONTRAST WITHIN PERMANOVA - OVERALL COMMUNITY
#https://thebiobucket.blogspot.com/2011/08/two-way-permanova-adonis-with-custom.html#more
# 1st factor = treatment:
nsSpp3$Livestockdensity<- factor(nsSpp3$Livestockdensity, levels(nsSpp3$Livestockdensity)[c(2,3,1)])
nsSpp3$Treatment<- factor(nsSpp3$Treatment, levels(nsSpp3$Treatment)[c(2,1)])
treat <- nsSpp3$Livestockdensity
levels(nsSpp3$Livestockdensity)
# 2nd factor = impact:
imp <- nsSpp3$Treatment

# simulating effect -
# simulation will add similar effects

## create a design matrix of the contrasts for "imp"
contrasts(imp) <- c(-1, 1)
Imp <- model.matrix(~ imp)[, -1]

## create a design matrix of the contrasts for "treat"
contrasts(treat) <- cbind(c(0,1,0),c(0,0,1))
Treat <- model.matrix(~ treat)[, -1]

imp.in.t1 <- Imp * ifelse(treat == "t1", 1, 0)
imp.in.t2 <- Imp * Treat[, 1]
imp.in.t3 <- Imp * Treat[, 2]

## specify the orthogonal contrasts for "treat"
contrasts(treat) <- cbind(c(1, -1, 0), c(1, 0, -1))

## specify the design matrix of the orthogonal
## contrasts for "treat"
Treat.ortho <- model.matrix(~ treat)[, -1]

## create a factor for each of the orthogonal "treat" contrasts
treat1vs2 <- Treat.ortho[, 1]
treat1vs3 <- Treat.ortho[, 2]

## do the pm-manova with the full model
fm1 <- adonis(vare.dis2~ treat * imp, method = "bray", perm = 999)

## do the pm-manova with the orthogonal contrasts for imp and treat'
## and the interaction contrasts of interest
fm2 <- adonis(vare.dis2~ treat1vs2 + treat1vs3 +
                imp.in.t1 + imp.in.t2 + imp.in.t3,
              method = "bray", perm = 999)
fm1; fm2

# High versus low and medium - (High first) + treatment interactions
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat1vs2   1     3.369  3.3692 12.9304 0.04340  0.001 ***# High vs low
#treat1vs3   1     2.774  2.7738 10.6453 0.03573  0.001 *** # High vs medium
#imp.in.t2   1     2.033  2.0330  7.8023 0.02619  0.001 *** # High vs low excl x open
#imp.in.t3   1     0.401  0.4006  1.5373 0.00516  0.121    # High vs medium excl x open
#Residuals 265    69.050  0.2606         0.88951           
#Total     269    77.627                 1.00000   

# Medium versus high and low - (Medium first) + treatment interactions
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat1vs2   1     1.692  1.6921  6.5455 0.02180  0.001 ***
#treat1vs3   1     4.451  4.4509 17.2171 0.05734  0.001 ***
#imp.in.t2   1     0.944  0.9439  3.6511 0.01216  0.001 ***
#imp.in.t3   1     2.033  2.0330  7.8642 0.02619  0.001 ***
#Residuals 265    68.507  0.2585         0.88252           
#Total     269    77.627                 1.00000     

# Low versus medium and high - (Low first) + treatment interactions
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat1vs2   1     4.153  4.1532 15.6916 0.05350  0.001 *** # Low vs medium
#treat1vs3   1     1.990  1.9898  7.5180 0.02563  0.001 *** # Low vs high
#imp.in.t2   1     0.401  0.4006  1.5134 0.00516  0.130     # Low vs medium exclosure
#imp.in.t3   1     0.944  0.9439  3.5661 0.01216  0.001 *** # Low vs high exclosure
#Residuals 265    70.139  0.2647         0.90355           

### check model with changed effects
eff <- sort(rep(1:6, 10))
eff[treat == "t1" | treat == "t2"] <- 3

# add noise:
dim(vare.dis2)
spdf <- matrix(NA, 270, 10, dimnames =
                 list(1:270, c("sp1", "sp2", "sp3", "sp4","sp5","sp6","sp7","sp8","sp9","sp10")))
spdf$sp1 = eff + rnorm(270, 0, 0.25)
spdf$sp2 = eff + rnorm(270, 0, 0.25)
spdf$sp3 = eff + rnorm(270, 0, 0.25)
spdf$sp4 = eff + rnorm(270, 0, 0.25)
spdf$sp5 = eff + rnorm(270, 0, 0.25)
spdf$sp6 = eff + rnorm(270, 0, 0.25)
spdf$sp7 = eff + rnorm(270, 0, 0.25)
spdf$sp8 = eff + rnorm(270, 0, 0.25)
spdf$sp9 = eff + rnorm(270, 0, 0.25)
spdf$sp10 = eff + rnorm(270, 0, 0.25)

# interaction plot for all species:
par(mfrow=c(2, 2), mar = c(2.5, 2.5, 0.5, 0.5))
for (i in 1:10) {interaction.plot(treat, imp, spdf[, i], lty = c(1, 2), legend = F);
  legend("topleft", bty = "n", cex = 1.2,
         paste("Species", i, sep = " "));
  legend("bottomright", c("no", "yes"),
         bty = "n", lty = c(2, 1))}

vare.disNULL<-vegdist(spdf,"bray") 
vare.dis2<-as.matrix(vare.dis)

fm1 <- adonis(spdf ~ treat * imp, method = "bray", perm = 999)
fm2 <- adonis(spdf ~ treat1vs2 + treat1vs3 +
                imp.in.t1 + imp.in.t2 + imp.in.t3,
              method = "euclidean", perm = 999)
fm1; fm2

#### Mean distance moved by each plot - overall ####
nsSpp3$NMDS1<-mdsNS$points[,1]
nsSpp3$NMDS2<-mdsNS$points[,2]
nsSpp31<-droplevels(nsSpp3[nsSpp3$Season=="Short I",])
nsSpp32<-droplevels(nsSpp3[nsSpp3$Season=="Long",])
nsSpp33<-droplevels(nsSpp3[nsSpp3$Season=="Short II",])

nsSpp31sub<-nsSpp31[,c("Quadrats","Season","NMDS1","NMDS2")]
nsSpp32sub<-nsSpp32[,c("Quadrats","Season","NMDS1","NMDS2")]
nsSpp33sub<-nsSpp33[,c("Quadrats","Season","NMDS1","NMDS2")]
dim(nsSpp32sub)

# Centroids harvest code - livestock and treatment 
modTLiv$centroids[,2] # Axis 2 centroids
modTLivM$centroids[,2]
modTLivT<-as.data.frame(modTLiv$centroids[,2] )
modTLivD<-as.data.frame(modTLivM$centroids[,2] )
modTLivD$harvest_code<- rownames(modTLivD)
colnames(modTLivD)[1]<-"Centroids"
nsSpp3merge<-merge(nsSpp3,modTLivD, by="harvest_code")


nsSpp31<-droplevels(nsSpp3merge[nsSpp3merge$Season=="Short I",])
nsSpp32<-droplevels(nsSpp3merge[nsSpp3merge$Season=="Long",])
nsSpp33<-droplevels(nsSpp3merge[nsSpp3merge$Season=="Short II",])
names(nsSpp31)
nsSpp31sub<-nsSpp31[,c("Quadrats","Season","NMDS1","NMDS2","Centroids")]
nsSpp32sub<-nsSpp32[,c("Quadrats","Season","NMDS1","NMDS2","Centroids")]
nsSpp33sub<-nsSpp33[,c("Quadrats","Season","NMDS1","NMDS2","Centroids")]
dim(nsSpp32sub)

# First season NMDS scores
cnt<-cbind(nsSpp31sub$NMDS1,nsSpp31sub$NMDS2)

# BCdist is row based? i.e. from row to row - therefore not necessary to transpose the data?? rbind? not cbind?
library(ecodist)
#euc.dist <- function(x1) sqrt(sum((x1 - cnt) ^ 2)) # bcdist= bray curtis
#NMDSdist2<-as.data.frame(cbind(apply(cbind(nsSpp31sub$NMDS1,nsSpp32sub$NMDS1),1, bcdist),apply(cbind(nsSpp31sub$NMDS2,nsSpp32sub$NMDS2),1, bcdist)))
#NMDSdist3<-as.data.frame(cbind(apply(cbind(nsSpp32sub$NMDS1,nsSpp33sub$NMDS1),1, bcdist),apply(cbind(nsSpp32sub$NMDS2,nsSpp33sub$NMDS2),1, bcdist)))
tnsSpp32sub<-t(nsSpp32sub$NMDS1)
tnsSpp32subC<-t(nsSpp32sub$Centroids)
tnsSpp32NMDS1<-rbind(tnsSpp32sub,tnsSpp32subC)
bcdist(tnsSpp32NMDS1)
NMDSdist2<-as.data.frame(cbind(apply(cbind(nsSpp32sub$NMDS1,nsSpp32sub$Centroids),1, bcdist),apply(cbind(nsSpp32sub$NMDS2,nsSpp32sub$Centroids),1, bcdist)))
NMDSdist3<-as.data.frame(cbind(apply(cbind(nsSpp33sub$NMDS1,nsSpp33sub$Centroids),1, bcdist),apply(cbind(nsSpp33sub$NMDS2,nsSpp33sub$Centroids),1, bcdist)))

nsSpp32<-cbind(nsSpp32,NMDSdist2)
nsSpp33<-cbind(nsSpp33,NMDSdist3)

colnames(nsSpp32)[92]<-"Bcdist1"
colnames(nsSpp32)[93]<-"Bcdist2"
colnames(nsSpp33)[92]<-"Bcdist1"
colnames(nsSpp33)[93]<-"Bcdist2"
reharv23<-rbind(nsSpp32,nsSpp33)

clnames <- cbind(reharv23$Bcdist1,reharv23$Bcdist2)
mean(reharv23$Bcdist1,reharv23$Bcdist2, data=reharv23)
reharv23$MeanBcdist<-rowMeans(clnames, na.rm = FALSE, dims = 1)
ggplot(reharv23,aes(x=Livestockdensity,y=Bcdist1, colour=Treatment))+geom_point()
ggplot(reharv23,aes(x=Livestockdensity,y=Bcdist2, colour=Treatment))+geom_point()
ggplot(reharv23,aes(x=Livestockdensity,y=MeanBcdist, colour=Treatment))+geom_point()

aggregate(abs(MeanBcdist)~Livestockdensity+Treatment+Season,reharv23,mean)
aggregate(abs(Bcdist2)~Livestockdensity+Treatment+Season,reharv23,mean)
aggregate(abs(Bcdist2)~Livestockdensity+Treatment+Season,reharv23,sem)

# Not working?

# Raw difference from centroid and quadrat NMDS2 score
nsSpp3mergeii<-droplevels(nsSpp3merge[nsSpp3merge$Season!="Short I",])
nsSpp3mergeii$distBC<-nsSpp3mergeii$NMDS2-nsSpp3mergeii$Centroids
ggplot(nsSpp3mergeii,aes(x=Livestockdensity,y=distBC, colour=Treatment, shape=Season))+geom_jitter()

nsSpp3mergeii$Bcdist2 <- abs(nsSpp3mergeii$distBC)
ggplot(nsSpp3mergeii,aes(x=Livestockdensity,y=Bcdist2, colour=Treatment, shape=Season))+geom_jitter()

aggregate(distBC~Livestockdensity+Treatment,nsSpp3mergeii,mean)

library(lme4)
BcDIST<-lmer(distBC~Livestockdensity+Season+Treatment+
               Livestockdensity:Season+Season:Treatment+
               Livestockdensity:Treatment+
               Treatment:Livestockdensity:Season+
               (1 |Block), data=nsSpp3mergeii)
summary(BcDIST)
anova(BcDIST)
AIC(BcDIST)
drop1(BcDIST)
hist(nsSpp3mergeii$distBC) #Analysis not good values bunched around zero - negative and positive
# Large outlier - test with and without  - low livestock
#nsSpp3G23b<-nsSpp3G23[-c(90,180),]


# Transform bray curtis distance to positive values - so difference from zero..
abs(reharv23$Bcdist)
reharv23$Bcdist1 <- abs(reharv23$Bcdist1)
reharv23$Bcdist2 <- abs(reharv23$Bcdist2)
ggplot(reharv23,aes(x=Livestockdensity,y=Bcdist2, colour=Treatment, shape =Season))+geom_jitter()

BetaSXsum<-aggregate(Bcdist2~Livestockdensity+Treatment+Season,reharv23,mean)
BetaSXsum

# Sum distance by quadrat
detach(package:plyr)
library(dplyr)
names(reharv23)
reharv23$time_code<-as.factor(with(reharv23, paste(Transect,Block,Treatment,Replicate, sep="-")))
QuadSUM1<-reharv23 %>% group_by(time_code,Livestockdensity,Treatment,Block) %>% summarise(Bcdist1 = sum(Bcdist1))
QuadSUM2<-reharv23 %>% group_by(time_code,Livestockdensity,Treatment,Block) %>% summarise(Bcdist2 = sum(Bcdist2))
ggplot(QuadSUM1,aes(x=Livestockdensity,y=Bcdist1, colour=Treatment))+geom_point()
ggplot(QuadSUM2,aes(x=Livestockdensity,y=Bcdist2, colour=Treatment))+geom_point()

# Mixed linear model - grass biomass and Bray curtis dist
library(glmmADMB)
Bcdist1<-glmmadmb(Bcdist1~Livestockdensity+Treatment+
                    Livestockdensity:Treatment+
                   +(1|Block), 
                 #admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                 family="gamma",data=QuadSUM1)
summary(Bcdist1)

Bcdist2<-glmmadmb(Bcdist2~Livestockdensity+Season+Treatment+
                    Livestockdensity:Season+Season:Treatment+
                    Livestockdensity:Treatment+
                    Treatment:Livestockdensity:Season+
                    +(1|Block), 
                  #admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                  family="gamma",data=reharv23)
summary(Bcdist2)
drop1(Bcdist2)




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

# Importance of rainfall - use ord.surf
par(mfrow=c(1,1))
ordi <- ordisurf(mdsNS  ~ rainmm, data =  nsSpp3, plot = FALSE, scaling = 3,
                           method = "ML", select = TRUE)
plot(mdsNS, type="n",xlim=c(-1,1), ylim=c(-1,1),
     ylab="NMDS 2", xlab="NMDS 1",mgp=c(1.75,.45,0), 
     tck=.02, las=1, lwd=1.75, bty='l', main="Rain (mm)")
plot(ordi, col = "dodgerblue1",lwd=1.5,npoints=6, labcex=.75, add = TRUE)


# Overall spp most asscoicated with specific livestock x treatment x season combination
nsSpp3.env<-nsSpp3[,c("fBlock","Livestockdensity","Season","Treatment")]
nsSpp3I<-droplevels(nsSpp3[nsSpp3$Date!="21.10.2012",])
nsSpp3.envI<-droplevels(nsSpp3.env[nsSpp3.env$Date!="21.10.2012",])

# Livestock density x treatment as a single factor
nsSpp3I$Liv_Trt<-as.factor(with(nsSpp3I, paste(Livestockdensity, sep="_")))

simT <- with(nsSpp3.envI, simper(nsSpp3I[,10:74],nsSpp3I$Liv_Trt),permutations=999)
SimpSum<-summary(simT,ordered = TRUE)
#sink("SIMPER.summary.txt")
#print(SimpSum)


#Contrast: High_Medium 
#                                           average        sd  ratio       ava       avb cumsum
#BothriochloainsculptaHochstExARichACamus 1.622e-01 0.1425744 1.1376 15.266667 19.366667 0.2442
#ChrysopogonplumulosusHochst              1.473e-01 0.1161515 1.2684 14.766667 18.383333 0.4660
#TriumfettaflavescensHochst               6.768e-02 0.0949657 0.7127  3.283333  6.075000 0.5679
#HeteropogoncontortusLRoemSchult          5.011e-02 0.0920001 0.5447  3.950000  2.666667 0.6433
#TetrapogonvillosusDesj                   2.991e-02 0.0819890 0.3648  0.000000  3.250000 0.6883
#DigitariamacroblepharaHackStapf          2.877e-02 0.0322177 0.8929  3.033333  1.333333 0.7316
#LintonianutansStapf                      2.082e-02 0.0383072 0.5434  2.000000  0.016667 0.7630 
#RhynchosiaminimaLDC                      1.866e-02 0.0486582 0.3835  0.716667  1.766667 0.7910  # 80%
#BarleriahomoiotrichaCBClarke             1.759e-02 0.0399106 0.4408  1.716667  0.250000 0.8175
#CynodonnlemfuensisVanderyst              1.600e-02 0.0573391 0.2791  0.816667  1.083333 0.8416

#Contrast: High_Low 

#Contrast: High_Low 
#average        sd  ratio       ava       avb cumsum
#BothriochloainsculptaHochstExARichACamus 1.365e-01 0.1298052 1.0519 15.266667 10.333333 0.1718
#ChrysopogonplumulosusHochst              1.320e-01 0.1221284 1.0806 14.766667 11.133333 0.3379
#CynodonnlemfuensisVanderyst              1.099e-01 0.1932534 0.5686  0.816667 14.216667 0.4762
#HeteropogoncontortusLRoemSchult          7.837e-02 0.1527033 0.5132  3.950000  7.341667 0.5748
#TriumfettaflavescensHochst               5.658e-02 0.0862453 0.6560  3.283333  5.316667 0.6460
#RhynchosiaminimaLDC                      3.165e-02 0.0404309 0.7829  0.716667  3.808333 0.6858
#DigitariamacroblepharaHackStapf          2.706e-02 0.0327213 0.8271  3.033333  0.416667 0.7199
#DyschoristeradicansNees                  2.400e-02 0.0451341 0.5317  1.425000  1.683333 0.7501 
#AristidakenyensisHenr                    2.274e-02 0.0724896 0.3137  0.416667  2.166667 0.7787
#LintonianutansStapf                      2.269e-02 0.0346046 0.6558  2.000000  0.900000 0.8073 # 80%

#Contrast: Medium_Low 
#average       sd  ratio       ava       avb cumsum
#BothriochloainsculptaHochstExARichACamus 0.1466790 0.123793 1.1849 19.366667 10.333333 0.1944
#ChrysopogonplumulosusHochst              0.1317967 0.112944 1.1669 18.383333 11.133333 0.3691
#CynodonnlemfuensisVanderyst              0.1060495 0.186373 0.5690  1.083333 14.216667 0.5097
#HeteropogoncontortusLRoemSchult          0.0697258 0.136581 0.5105  2.666667  7.341667 0.6021
#TriumfettaflavescensHochst               0.0588615 0.069250 0.8500  6.075000  5.316667 0.6802
#TetrapogonvillosusDesj                   0.0350416 0.089707 0.3906  3.250000  1.666667 0.7266
#RhynchosiaminimaLDC                      0.0339460 0.048535 0.6994  1.766667  3.808333 0.7716
#AristidakenyensisHenr                    0.0201779 0.067430 0.2992  0.166667  2.166667 0.7984 # 80% here
#SetariaincrassataHochstWild              0.0201162 0.091928 0.2188  0.000000  2.833333 0.8250 
#DyschoristeradicansNees                  0.0152932 0.032773 0.4666  0.116667  1.683333 0.8453

nsSpp3$plot_codeX<-as.factor(with(nsSpp3, paste(Livestockdensity,Treatment,Season, sep="_")))
nsSpp3I$plot_codeX<-as.factor(with(nsSpp3I, paste(Livestockdensity,Treatment,Season, sep="_")))

#### Key species ####

# Chr Plu
ChrPlu<-ggplot(nsSpp3I,aes(x=plot_codeX,y=ChrysopogonplumulosusHochst,shape=Treatment,colour=Livestockdensity, linetype=Season))
ChrPlu<-ChrPlu+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
ChrPlu # Higher in high and moderate liverstock - lower in lo

plot(GrassNetReharvestBiomass1~ChrysopogonplumulosusHochst,col=c(nsSpp3G$Livestockdensity),nsSpp3G)
summary(lm(GrassNetReharvestBiomass1~ChrysopogonplumulosusHochst,nsSpp3G))

# ChrPlu Mean cover and occurrence
nsSpp3Cp<-nsSpp3[nsSpp3$ChrysopogonplumulosusHochs!=0,]
mean(nsSpp3Cp$ChrysopogonplumulosusHochs) # 19.0% 
sd(nsSpp3Cp$ChrysopogonplumulosusHochs) # 16.3% 
nrow(nsSpp3[nsSpp3$ChrysopogonplumulosusHochs>0.01,])/270*100 # 80 %

# Bot ins
BotIns<-ggplot(nsSpp3I,aes(x=plot_codeX,y=BothriochloainsculptaHochstExARichACamus,shape=Treatment,colour=Livestockdensity, linetype=Season))
BotIns<-BotIns+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
BotIns # Higher in high and moderate liverstock - lower in low
plot(GrassNetReharvestBiomass1~BothriochloainsculptaHochstExARichACamus,nsSpp3G)
summary(lm(GrassNetReharvestBiomass1~BothriochloainsculptaHochstExARichACamus,nsSpp3G))

# BotIns Mean cover and occurrence
nsSpp3Bi<-nsSpp3[nsSpp3$BothriochloainsculptaHochstExARichACamus!=0,]
mean(nsSpp3Bi$BothriochloainsculptaHochstExARichACamus) # 24.7%
sd(nsSpp3Bi$BothriochloainsculptaHochstExARichACamus) # 22.5% 
nrow(nsSpp3[nsSpp3$BothriochloainsculptaHochstExARichACamus>0.01,])/270*100 # 74.1 %

# Dig mac
Digmac<-ggplot(nsSpp3I,aes(x=plot_codeX,y=DigitariamacroblepharaHackStapf  ,shape=Treatment,colour=Livestockdensity, linetype=Season))
Digmac<-Digmac+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
Digmac # High livestock and moderae - absent in low 

aggregate(DigitariamacroblepharaHackStapf~Liv_Trt,nsSpp3I,mean)
aggregate(DigitariamacroblepharaHackStapf~Liv_Trt,nsSpp3I,sd)

nsSpp3Dm<-nsSpp3[nsSpp3$DigitariamacroblepharaHackStapf!=0,]
mean(nsSpp3Dm$DigitariamacroblepharaHackStapf) # 4.9%
sd(nsSpp3Dm$DigitariamacroblepharaHackStapf) # 4.0% 
nrow(nsSpp3[nsSpp3$DigitariamacroblepharaHackStapf>0.01,])/270*100 # 33.3 %

# HetCon
HetCon<-ggplot(nsSpp3I,aes(x=plot_codeX,y=HeteropogoncontortusLRoemSchult,shape=Treatment,colour=Livestockdensity, linetype=Season))
HetCon<-HetCon+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
HetCon # Contrast low exclosure none - open high 
# High in high exclosures 
plot(GrassNetReharvestBiomass1~HeteropogoncontortusLRoemSchult,nsSpp3G)

nsSpp3Hc<-nsSpp3[nsSpp3$HeteropogoncontortusLRoemSchult!=0,]
mean(nsSpp3Hc$HeteropogoncontortusLRoemSchult) # 19.7%
sd(nsSpp3Hc$HeteropogoncontortusLRoemSchult) # 30.0% 
nrow(nsSpp3[nsSpp3$HeteropogoncontortusLRoemSchult>0.01,])/270*100 # 25.2 %

#CynNle
CynNle<-ggplot(nsSpp3I,aes(x=plot_codeX,y=CynodonnlemfuensisVanderyst,shape=Treatment,colour=Livestockdensity, linetype=Season))
CynNle<-CynNle+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
CynNle # Huge increase in low livestock exclosures
# Some occurence of this grass in moderate at last season - though also present in second season...
plot(GrassNetReharvestBiomass1~CynodonnlemfuensisVanderyst,nsSpp3G)
abline(lm(GrassNetReharvestBiomass1~CynodonnlemfuensisVanderyst,nsSpp3G))
summary(lm(GrassNetReharvestBiomass1~CynodonnlemfuensisVanderyst,nsSpp3G))

aggregate(CynodonnlemfuensisVanderyst~Liv_Trt,nsSpp3I,mean)
aggregate(CynodonnlemfuensisVanderyst~Liv_Trt,nsSpp3I,sd)

nsSpp3Cn<-nsSpp3[nsSpp3$CynodonnlemfuensisVanderyst!=0,]
mean(nsSpp3Cn$CynodonnlemfuensisVanderyst) # 27.4%
sd(nsSpp3Cn$CynodonnlemfuensisVanderyst) # 29.5% 
nrow(nsSpp3[nsSpp3$CynodonnlemfuensisVanderyst>0.01,])/270*100 # 21.1 %

#Lin ton
Linton<-ggplot(nsSpp3I,aes(x=plot_codeX,y=LintonianutansStapf,shape=Treatment,colour=Livestockdensity, linetype=Season))
Linton<-Linton+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
Linton # High in high livestock and outside in low livestock - none in moderate

nsSpp3Lt<-nsSpp3[nsSpp3$LintonianutansStapf!=0,]
mean(nsSpp3Lt$LintonianutansStapf) #6.2%
sd(nsSpp3Lt$LintonianutansStapf) # 10.5% 
nrow(nsSpp3[nsSpp3$LintonianutansStapf>0.01,])/270*100 # 21.5 %

#AriKen
AriKen<-ggplot(nsSpp3I,aes(x=plot_codeX,y=AristidakenyensisHenr,shape=Treatment,colour=Livestockdensity, linetype=Season))
AriKen<-AriKen+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
AriKen # Contrast low exclosure none - open high 
# High in high exclosures 

nsSpp3Ak<-nsSpp3[nsSpp3$AristidakenyensisHenr!=0,]
mean(nsSpp3Ak$AristidakenyensisHenr) # 13.3%
sd(nsSpp3Ak$AristidakenyensisHenr) # 17.0% 
nrow(nsSpp3[nsSpp3$AristidakenyensisHenr>0.01,])/270*100 # 7.8 %

#Tetvil
Tetvil<-ggplot(nsSpp3I,aes(x=plot_codeX,y=TetrapogonvillosusDesj,shape=Treatment,colour=Livestockdensity, linetype=Season))
Tetvil<-Tetvil+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
Tetvil

nsSpp3Tv<-nsSpp3[nsSpp3$TetrapogonvillosusDesj!=0,]
mean(nsSpp3Tv$TetrapogonvillosusDesj) # 19.5%
sd(nsSpp3Tv$TetrapogonvillosusDesj) # 23.0% 
nrow(nsSpp3[nsSpp3$TetrapogonvillosusDesj>0.01,])/270*100 # 10.7 %

aggregate(TetrapogonvillosusDesj~Livestockdensity+Treatment+Season,nsSpp3Tv,mean)

#Setinc
SetariaincrassataHochstWild
nsSpp3Si<-nsSpp3[nsSpp3$SetariaincrassataHochstWild!=0,]
mean(nsSpp3Si$SetariaincrassataHochstWild) # 22.6%
sd(nsSpp3Si$SetariaincrassataHochstWild) # 30.9% 
nrow(nsSpp3[nsSpp3$SetariaincrassataHochstWild>0.01,])/270*100 # 3.3 %

#### Forbs ####
#RhyMin
RhyMin<-ggplot(nsSpp3I,aes(x=plot_codeX,y=RhynchosiaminimaLDC,shape=Treatment,colour=Livestockdensity, linetype=Season))
RhyMin<-RhyMin+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
RhyMin # Higher in low than moderate or high livestock - particularly exclosures

nsSpp3Rm<-nsSpp3[nsSpp3$RhynchosiaminimaLDC!=0,]
mean(nsSpp3Rm$RhynchosiaminimaLDC) # 3.9%
sd(nsSpp3Rm$RhynchosiaminimaLDC) # 6.0% 
nrow(nsSpp3[nsSpp3$RhynchosiaminimaLDC>0.01,])/270*100 # 45.6 %

#DysRad
DysRad<-ggplot(nsSpp3I,aes(x=plot_codeX,y=DyschoristeradicansNees,shape=Treatment,colour=Livestockdensity, linetype=Season))
DysRad<-DysRad+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
DysRad # Higher in low than moderate or high livestock - particularly exclosures

nsSpp3Dr<-nsSpp3[nsSpp3$DyschoristeradicansNees!=0,]
mean(nsSpp3Dr$DyschoristeradicansNees) # 4.0%
sd(nsSpp3Dr$DyschoristeradicansNees) # 4.8% 
nrow(nsSpp3[nsSpp3$DyschoristeradicansNees>0.01,])/270*100 # 18.1 %

#Barhom
nsSpp3Bh<-nsSpp3[nsSpp3$BarleriahomoiotrichaCBClarke!=0,]
mean(nsSpp3Bh$BarleriahomoiotrichaCBClarke) # 3.7%
sd(nsSpp3Bh$BarleriahomoiotrichaCBClarke) # 5.0% 
nrow(nsSpp3[nsSpp3$BarleriahomoiotrichaCBClarke>0.01,])/270*100 # 21.5 %

# Dwarf shrub species
#Trifla
Trifla<-ggplot(nsSpp3I,aes(x=plot_codeX,y=TriumfettaflavescensHochst,shape=Treatment,colour=Livestockdensity, linetype=Season))
Trifla<-Trifla+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
Trifla # Low livestock exclosures and moderate livestock

nsSpp3Tf<-nsSpp3[nsSpp3$TriumfettaflavescensHochst!=0,]
mean(nsSpp3Tf$TriumfettaflavescensHochst) # 7.5%
sd(nsSpp3Tf$TriumfettaflavescensHochst) # 22.5% 
nrow(nsSpp3I[nsSpp3I$TriumfettaflavescensHochst>0.01,])/270*100 # 57.8 %


### Species occurence
library(dplyr)
library(tidyr)
library(reshape)

nsSpp3Long<-melt(nsSpp3I, id.vars=c("Season","Transect","Treatment","Livestockdensity","Quadrats", "Block","Trtname", "Replicate" , "Date",
 "plot_code","GrassNetReharvestBiomass1","DwarfShrubNetReharvestBiomass1"          
,"HerbNetReharvestBiomass1","ClimberNetReharvestBiomass1", "TotalBiomass1", "rainmm", "fBlock", "plot_codeX" )) #,"Liv_Trt"

nsSpp3Long<-droplevels(nsSpp3Long[nsSpp3Long$value>0.01,])
colnames(nsSpp3Long)[19:20]<-c("spp","per_cover")
names(nsSpp3Long)

nsSpp3Long$plot_id<-as.factor(with(nsSpp3Long, paste(Livestockdensity,Treatment,Replicate, sep="_")))
levels(nsSpp3Long$plot_id)

nsSpp3Long$preabs<-ifelse(nsSpp3Long$per_cover>0.01,1,0)

table(nsSpp3Long$spp)

SppOccur<-nsSpp3Long %>%
  group_by(spp) %>% 
  # Count occurrences per species
  summarise(n = n()) %>%
  # Get percent per occurence (number quadrats over experiment = 270)
  mutate(percent = n / 270 * 100) %>%
  # Select results per species
  select(spp, percent) 
SppOccur<-as.data.frame(SppOccur)
SppOccur

# Overall spp most asscoicated with specific livestock x treatment x season combination
nsSpp3.env<-nsSpp3[,c("fBlock","Livestockdensity","Season","Treatment")]
nsSpp3I<-droplevels(nsSpp3[nsSpp3$Date!="21.10.2012",])
nsSpp3.envI<-droplevels(nsSpp3.env[nsSpp3.env$Date!="21.10.2012",])

# Livestock density x treatment as a single factor
nsSpp3I$Liv_Trt<-as.factor(with(nsSpp3I, paste(Livestockdensity, sep="_")))

simT <- with(nsSpp3.envI, simper(nsSpp3I[,10:74],nsSpp3I$Liv_Trt),permutations=999)
SimpSum<-summary(simT,ordered = TRUE)
#sink("SIMPER.summary.txt")
print(SimpSum)

class(SimpSum)
SimpSumDATA<-as.data.frame(SimpSum)
SimpSumDATA<-data.frame(unclass(SimpSum), check.names = FALSE, stringsAsFactors = FALSE)

names(SimpSumDATA)
SimpSumDATA$spp <- rownames(SimpSumDATA)
SimpSumDATA$spp<-as.factor(SimpSumDATA$spp)
class(SppOccur$spp)

SimpOccur<-merge(SimpSumDATA,SppOccur, by="spp")
names(SimpOccur)
SimpOccur$sppAB<-abbreviate(SimpOccur$spp)
plot(High_Low.average~percent,SimpOccur)
abline(lm(High_Low.average~percent,SimpOccur))

ggplot(SimpOccur, aes(x=percent, y= High_Low.average*100, label=sppAB))+
geom_text(aes(label=sppAB),hjust=0, vjust=0)
##########################################################################
#### Grasses ####
##########################################################################
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
# PERMANOVA/ADONIS - GRASSES
# Distance matrix
#nsSpp3GSI<-droplevels(nsSpp3G[nsSpp3G$Date!="21.10.2012",])
vare.disG<-vegdist(nsSpp3G[,10:28], "bray") 
vare.dis2G<-as.matrix(vare.disG)

# Beta Dispersion
nsSpp3G$harvest_code<-as.factor(with(nsSpp3G, paste(Livestockdensity,Treatment,Season, sep="-")))
modTLivG <- betadisper(vare.disG, nsSpp3G$harvest_code,type = c("centroid"))
modTLivG

# PERMANOVA grasses
nsSpp3G$time_code<-as.factor(with(nsSpp3G, paste(Transect,Block,Treatment,Replicate, sep="-")))
PermG<-adonis(vare.disG~Livestockdensity+Treatment+Season+rainmm+
                Season:Treatment+Livestockdensity:Season+
                Livestockdensity:Treatment+rainmm:Season+
                rainmm:Livestockdensity+rainmm:Treatment+
                Treatment:Livestockdensity:Season+
                Treatment:Livestockdensity:rainmm, 
              strata=as.numeric(nsSpp3G$Block),#/as.numeric(nsSpp3G$time_code),
              method = "bray",perm=999, data=nsSpp3G)
PermG$aov.tab
PermG$aov.tab$R2

# Extract centroids
#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258
#vec.sp.dfG<-as.data.frame(cbind(NS.harG$factors$centroids*sqrt(NS.harG$factors$r)))
vec.sp.dfG<-data.frame(cbind(modTLivG$centroids[,1],modTLivG$centroids[,2]))
colnames(vec.sp.dfG)[1]<-"NMDS1"
colnames(vec.sp.dfG)[2]<-"NMDS2"
#vec.sp.dfG<-as.data.frame(scores(mdsNSG, display = "sites"))
vec.sp.dfG$Livestockdensity<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                                "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.dfG$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                        "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.dfG$Season<-rep(c("Long","Short I","Short II"))
vec.sp.dfG$harvest_code<-as.factor(with(vec.sp.dfG, paste(Livestockdensity,Treatment,Season, sep="-")))

TotBio<-aggregate(TotalBiomass1~harvest_code,nsSpp3,mean)
vec.sp.df<-merge(vec.sp.df, TotBio, by.x = "harvest_code")
vec.sp.df$TotalBiomass1<-as.numeric(vec.sp.df$TotalBiomass1)

GReharvestBio<-aggregate(GrassNetReharvestBiomass1~harvest_code,nsSpp3G,mean)
vec.sp.dfG<-merge(vec.sp.dfG, GReharvestBio, by.x = "harvest_code")
vec.sp.dfG$GrassNetReharvestBiomass1<-as.numeric(vec.sp.dfG$GrassNetReharvestBiomass1)

#NS.harG$scores
#Grass species scores
spp.scrs <- as.data.frame(scores(mdsNSG, display = "species"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Livestockdensity<-c("Low")
spp.scrs$Season<-c("Short I")
spp.scrs$Treatment<-c("Control")

# Reorder by season and Harvest code
vec.sp.dfG$Season<- factor(vec.sp.dfG$Season, levels = c("Short I","Long","Short II"))
vec.sp.dfG$Season2<-c(1,2,3)
vvec.sp.dfG<- vec.sp.dfG[order(vec.sp.dfG$Season2),] 
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
CenPlotG<-ggplot(vec.sp.dfG[order(vec.sp.dfG$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity))
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

#### CONTRASTS WITHIN GRASS PERMANOVA ####
# 1st factor = treatment:
treat <- nsSpp3G$Livestockdensity

# 2nd factor = impact:
imp <- nsSpp3G$Treatment

# simulating effect -
# simulation will add similar effects

## create a design matrix of the contrasts for "imp"
contrasts(imp) <- c(-1, 1)
Imp <- model.matrix(~ imp)[, -1]

## create a design matrix of the contrasts for "treat"
contrasts(treat) <- cbind(c(0,1,0),c(0,0,1))
Treat <- model.matrix(~ treat)[, -1]

imp.in.t1 <- Imp * ifelse(treat == "t1", 1, 0)
imp.in.t2 <- Imp * Treat[, 1]
imp.in.t3 <- Imp * Treat[, 2]

## specify the orthogonal contrasts for "treat"
contrasts(treat) <- cbind(c(1, -1, 0), c(1, 0, -1))

## specify the design matrix of the orthogonal
## contrasts for "treat"
Treat.ortho <- model.matrix(~ treat)[, -1]

## create a factor for each of the orthogonal "treat" contrasts
treat1vs2 <- Treat.ortho[, 1]
treat1vs3 <- Treat.ortho[, 2]

## do the pm-manova with the full model
fm1 <- adonis(vare.dis2G~ treat * imp, method = "bray", perm = 999)

## do the pm-manova with the orthogonal contrasts for imp and treat'
## and the interaction contrasts of interest
fm2 <- adonis(vare.dis2G~ treat1vs2 + treat1vs3 +
                imp.in.t1 + imp.in.t2 + imp.in.t3,
              method = "bray", perm = 999)
fm1; fm2

# Contrast of high livestock versus low and medium + treatment interactions
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat1vs2   1     3.082 3.08222 12.0766 0.04072  0.001 ***# High vs low
#treat1vs3   1     2.691 2.69139 10.5452 0.03556  0.001 ***# High vs medium
#imp.in.t2   1     2.055 2.05539  8.0533 0.02716  0.001 ***# High vs low excl x open
#imp.in.t3   1     0.223 0.22279  0.8729 0.00294  0.537# High vs low excl x open    
#Residuals 265    67.634 0.25522         0.89362           
#Total     269    75.686                 1.00000    

#  Contrast with high livestock
# interaction between high and  low livestock and treatment - 
# Issues - this is the same as total biomass - I do not believe the output...
# no interaction with high and medium and treatment
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat1vs2   1     3.082 3.08222 12.0766 0.04072  0.001 ***
#treat1vs3   1     2.691 2.69139 10.5452 0.03556  0.001 ***
#imp.in.t2   1     2.055 2.05539  8.0533 0.02716  0.001 ***
#imp.in.t3   1     0.223 0.22279  0.8729 0.00294  0.495    
#####

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

test.NMDSdist2<-as.data.frame(cbind(apply(cbind(nsSpp3G1sub$NMDS1,nsSpp3G2sub$NMDS1),1, bcdist),apply(cbind(nsSpp3G1sub$NMDS2,nsSpp3G2sub$NMDS2),1, bcdist)))
test.NMDSdist3<-as.data.frame(cbind(apply(cbind(nsSpp3G2sub$NMDS1,nsSpp3G3sub$NMDS1),1, bcdist),apply(cbind(nsSpp3G2sub$NMDS2,nsSpp3G3sub$NMDS2),1, bcdist)))

nsSpp3G2<-cbind(nsSpp3G2,test.NMDSdist2)
nsSpp3G3<-cbind(nsSpp3G3,test.NMDSdist3)
colnames(nsSpp3G2)[39]<-"Bcdist1"
colnames(nsSpp3G2)[40]<-"Bcdist2"
colnames(nsSpp3G3)[39]<-"Bcdist1"
colnames(nsSpp3G3)[40]<-"Bcdist2"
nsSpp3G23<-rbind(nsSpp3G2,nsSpp3G3)
names(nsSpp3G23)

# Transform bray curtis distance to positive values - so difference from zero..
abs(nsSpp3G23$Bcdist1)
nsSpp3G23$Bcdist1 <- abs(nsSpp3G23$Bcdist1)
nsSpp3G23$Bcdist2 <- abs(nsSpp3G23$Bcdist2)

ggplot(nsSpp3G23,aes(x=Livestockdensity,y=Bcdist1, colour=Treatment,shape=Season))+geom_boxplot()+scale_y_continuous(limits=c(0,10))
ggplot(nsSpp3G23,aes(x=Livestockdensity,y=Bcdist2, colour=Treatment, shape=Season))+geom_boxplot()+scale_y_continuous(limits=c(0,10))

nsSpp3G23b<-nsSpp3G23[-c(90,180),]

# Sum distance by quadrat
library(dplyr)
nsSpp3G23$time_code<-as.factor(with(nsSpp3G23, paste(Transect,Block,Treatment,Replicate, sep="-")))
QuadSUMG1<-nsSpp3G23 %>% group_by(time_code,Livestockdensity,Treatment,Block) %>% summarise(Bcdist1 = sum(Bcdist1))
QuadSUMG2<-nsSpp3G23 %>% group_by(time_code,Livestockdensity,Treatment,Block) %>% summarise(Bcdist2 = sum(Bcdist2))
ggplot(QuadSUMG1,aes(x=Livestockdensity,y=Bcdist1, colour=Treatment))+geom_boxplot()+scale_y_continuous(limits=c(0,10))
ggplot(QuadSUMG2,aes(x=Livestockdensity,y=Bcdist2, colour=Treatment))+geom_boxplot()+scale_y_continuous(limits=c(0,10))

# Mixed linear model - grass biomass and Eudist
library(glmmADMB)
#nsSpp3G23$fBlock<-as.factor(nsSpp3G23$Block)
#nsSpp3G23$Bcdist1[nsSpp3G23$Bcdist1==0]<-.01
#nsSpp3G23$Bcdist2[nsSpp3G23$Bcdist2==0]<-.01

BcdistG1<-glmmadmb(Bcdist1~Livestockdensity+Treatment+
                    Livestockdensity:Treatment+
                    +(1|Block), 
                  #admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                  family="gamma",data=QuadSUMG1)
summary(BcdistG1)
drop1(BcdistG1, test="Chi")
BcdistG2<-glmmadmb(Bcdist2~Livestockdensity+Treatment+
                    Livestockdensity:Treatment+
                    +(1|Block), 
                  #admb.opts=admbControl(shess=FALSE,noinit=FALSE,impSamp=200,maxfn=1000,imaxfn=500,maxph=5),
                  family="gamma",data=QuadSUMG2)
summary(BcdistG2)
drop1(BcdistG2, test="Chi")

summary(BcdistLMG)
summary(BcdistLMG2)

drop1(BcdistLMG, test="Chi")
#Livestockdensity:Treatment  2 171.51  3.5748  0.167395  
#Livestockdensity:Treatment:Season  2 171.93 4.358   0.1132

drop1(BcdistLMG2, test="Chi")

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

nsSpp3G$plot_codeX<-as.factor(with(nsSpp3G, paste(Livestockdensity,Treatment,Season, sep="_")))

# High and medium livestock - lower in low livestock
Digmac<-ggplot(nsSpp3,aes(x=plot_codeX,y=DigitariamacroblepharaHackStapf  ,shape=Treatment,colour=Livestockdensity, linetype=Season))
Digmac<-Digmac+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
Digmac # High livestock

BotIns<-ggplot(nsSpp3G,aes(x=plot_codeX,y=BothriochloainsculptaHochstExARichACamus,shape=Treatment,colour=Livestockdensity, linetype=Season))
BotIns<-BotIns+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
BotIns
# Decline in Bot ins during third season moderate exclosures...
plot(GrassNetReharvestBiomass1~BothriochloainsculptaHochstExARichACamus,nsSpp3G)

HetCon<-ggplot(nsSpp3G,aes(x=plot_codeX,y=HeteropogoncontortusLRoemSchult,shape=Treatment,colour=Livestockdensity, linetype=Season))
HetCon<-HetCon+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
HetCon # Huge increase in low livesock exclosures - season II
# High in high exclosures and low open
plot(GrassNetReharvestBiomass1~HeteropogoncontortusLRoemSchult,nsSpp3G)

# Medium and livestock higher
ChrPlu<-ggplot(nsSpp3G,aes(x=plot_codeX,y=ChrysopogonplumulosusHochst,shape=Treatment,colour=Livestockdensity, linetype=Season))
ChrPlu<-ChrPlu+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
ChrPlu
plot(GrassNetReharvestBiomass1~ChrysopogonplumulosusHochst,col=c(nsSpp3G$Livestockdensity),nsSpp3G)
# Decline in Chr Plu during thir season moderate exclosures...
summary(lm(GrassNetReharvestBiomass1~ChrysopogonplumulosusHochst,nsSpp3G))

# Low livestock - higher - in exclosures
CynNle<-ggplot(nsSpp3G,aes(x=plot_codeX,y=CynodonnlemfuensisVanderyst,shape=Treatment,colour=Livestockdensity, linetype=Season))
CynNle<-CynNle+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
CynNle # Huge increase in low livestock exclosures- season III
# Some occurence of this grass in moderate at last season - though also present in second season...
plot(GrassNetReharvestBiomass1~CynodonnlemfuensisVanderyst,nsSpp3G)
abline(lm(GrassNetReharvestBiomass1~CynodonnlemfuensisVanderyst,nsSpp3G))
summary(lm(GrassNetReharvestBiomass1~CynodonnlemfuensisVanderyst,nsSpp3G))

CynNle<-ggplot(nsSpp3G,aes(x=plot_codeX,y=CynodonnlemfuensisVanderyst,shape=Treatment,colour=Livestockdensity, linetype=Season))
CynNle<-CynNle+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
CynNle


# Others
nsSpp3$plot_codeX<-as.factor(with(nsSpp3, paste(Livestockdensity,Treatment,Season, sep="_")))

RhyMin<-ggplot(nsSpp3,aes(x=plot_codeX,y=RhynchosiaminimaLDC,shape=Treatment,colour=Livestockdensity, linetype=Season))
RhyMin<-RhyMin+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
RhyMin 

Trifla<-ggplot(nsSpp3,aes(x=plot_codeX,y=TriumfettaflavescensHochst,shape=Treatment,colour=Livestockdensity, linetype=Season))
Trifla<-Trifla+geom_boxplot(outlier.shape=NA,fill=NA,show.legend=F)+geom_jitter(size=2.5,stroke=1,show.legend=T)
Trifla # Low livestock and moderate # Dwardf shrub

##########################################################################
#### Woody ####
##########################################################################
# USE NMDS on Woody only

# Remove rows with sum = zero
nsSpp3Wb<- nsSpp3W[ rowSums(nsSpp3W[,10:19])!=0, ] 

# NMDS woody biomass
mdsNSW<-metaMDS(nsSpp3Wb[,10:19], trace =F) 
mdsNSW # 0.09492855 

NS.evW <- envfit(mdsNSW~Livestockdensity+Treatment+ Season,data = nsSpp3Wb, 
                 perm=999,strata=as.numeric(nsSpp3Wb$Block))
NS.evW # Livestock density and treatment significant - not season

nsSpp3Wb$harvest_code<-as.factor(with(nsSpp3Wb, paste(Livestockdensity,Treatment,Date, sep="-")))
NS.harW <- envfit(mdsNSW~ harvest_code,data = nsSpp3Wb, 
                  perm=999,strata=as.numeric(nsSpp3Wb$Block))
NS.harW
# PERMANOVA/ADONIS - WOODY PLANTS
# Distance matrix
#nsSpp3GSI<-droplevels(nsSpp3G[nsSpp3G$Date!="21.10.2012",])
vare.disW<-vegdist(nsSpp3Wb[,10:19], "bray") 
vare.dis2W<-as.matrix(vare.disW)

# Beta Dispersion
nsSpp3Wb$harvest_code<-as.factor(with(nsSpp3Wb, paste(Livestockdensity,Treatment,Season, sep="-")))
modTLivW <- betadisper(vare.disW, nsSpp3Wb$harvest_code,type = c("centroid"))
modTLivW

# PERMANOVA grasses
nsSpp3Wb$time_code<-as.factor(with(nsSpp3Wb, paste(Transect,Block,Treatment,Replicate, sep="-")))
PermW<-adonis(vare.disW~Livestockdensity+Treatment+Season+rainmm+
                Season:Treatment+Livestockdensity:Season+
                Livestockdensity:Treatment+rainmm:Season+
                rainmm:Livestockdensity+rainmm:Treatment+
                Treatment:Livestockdensity:Season+
                Treatment:Livestockdensity:rainmm, 
              strata=as.numeric(nsSpp3Wb$Block),#/as.numeric(nsSpp3G$time_code),
              method = "bray",perm=999, data=nsSpp3Wb)
PermW$aov.tab
PermW$aov.tab$R2

# Extract centroids
#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258
#vec.sp.dfG<-as.data.frame(cbind(NS.harG$factors$centroids*sqrt(NS.harG$factors$r)))
vec.sp.dfW<-data.frame(cbind(modTLivW$centroids[,1],modTLivW$centroids[,2]))
colnames(vec.sp.dfW)[1]<-"NMDS1"
colnames(vec.sp.dfW)[2]<-"NMDS2"
#vec.sp.dfG<-as.data.frame(scores(mdsNSG, display = "sites"))
vec.sp.dfW$Livestockdensity<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                               "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.dfW$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                        "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.dfW$Season<-rep(c("Long","Short I","Short II"))
vec.sp.dfW$harvest_code<-as.factor(with(vec.sp.dfW, paste(Livestockdensity,Treatment,Season, sep="-")))

WReharvestBio<-aggregate(DwarfShrubNetReharvestBiomass1~harvest_code,nsSpp3Wb,mean)
vec.sp.dfW<-merge(vec.sp.dfW, WReharvestBio, by.x = "harvest_code")
vec.sp.dfW$DwarfShrubNetReharvestBiomass1<-as.numeric(vec.sp.dfW$DwarfShrubNetReharvestBiomass1)

#NS.harW$scores
#Woody species scores
spp.scrsW <- as.data.frame(scores(mdsNSW, display = "species"))
spp.scrsW <- cbind(spp.scrsW, Species = rownames(spp.scrsW))
spp.scrsW$Livestockdensity<-c("Low")
spp.scrsW$Season<-c("Short I")
spp.scrsW$Treatment<-c("Control")

# Reorder by season and Harvest code
vec.sp.dfW$Season<- factor(vec.sp.dfW$Season, levels = c("Short I","Long","Short II"))
vec.sp.dfW$Season2<-c(1,2,3)
vvec.sp.dfW<- vec.sp.dfW[order(vec.sp.dfW$Season2),] 
vec.sp.dfW$harvest_code<-as.factor(with(vec.sp.dfW, paste(Livestockdensity,Treatment, sep="-")))

vec.sp.dfW$Livestockdensity<-as.factor(vec.sp.dfW$Livestockdensity)
levels(vec.sp.dfW$Livestockdensity)

# Woody plants multivariate movement plot
CenPlotW<-ggplot(vec.sp.dfW[order(vec.sp.dfW$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity))
#CenPlotG<-CenPlotG+geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=Livestockdensity), size=1,show.legend=T)
CenPlotW<-CenPlotW+geom_point(aes(size=DwarfShrubNetReharvestBiomass1,shape=Treatment), stroke=1)
CenPlotW<-CenPlotW+geom_path(aes(group=harvest_code),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
#CenPlotW<-CenPlotW+geom_text(data=spp.scrs,aes(label=Species),colour="light grey")
CenPlotW<-CenPlotW+geom_text(aes(label=Season),hjust=0, vjust=-.95, show.legend = F)
CenPlotW<-CenPlotW+scale_colour_manual(values=c("black","grey80","grey50"))
CenPlotW<-CenPlotW+scale_fill_manual(values=c("black","grey80","grey50"))
CenPlotW<-CenPlotW+scale_radius(range=c(1,8))
CenPlotW<-CenPlotW+scale_shape_manual(values=c(21,22))
CenPlotW<-CenPlotW+ggtitle("Woody")
CenPlotW<-CenPlotW+#theme_bw() +
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
CenPlotW

##########################################################################
#### Herbs and climbers ####
##########################################################################
# USE NMDS on Herbs and climbers only
names(nsSpp3C) 
# Remove rows with sum = zero
nsSpp3Cb<- nsSpp3C[ rowSums(nsSpp3C[,10:38])!=0, ] 
names(nsSpp3Cb[,10:38])
# NMDS woody biomass
mdsNSC<-metaMDS(nsSpp3Cb[,10:38], trace =F) 
mdsNSC # 0.1615137  

NS.evC <- envfit(mdsNSC~Livestockdensity+Treatment+Season+rainmm,data = nsSpp3Cb, 
                 perm=999,strata=as.numeric(nsSpp3Cb$Block))
NS.evC # Only treatment marginally significant

nsSpp3Cb$harvest_code<-as.factor(with(nsSpp3Cb, paste(Livestockdensity,Treatment,Date, sep="-")))
NS.harC <- envfit(mdsNSC~ harvest_code,data = nsSpp3Cb, 
                  perm=999,strata=as.numeric(nsSpp3Cb$Block))
NS.harC # NS

# PERMANOVA/ADONIS - HERBS AND CLIMBERS PLANTS
# Distance matrix
vare.disC<-vegdist(nsSpp3Cb[,10:38], "bray") 
vare.dis2C<-as.matrix(vare.disC)

# Beta Dispersion
nsSpp3Cb$harvest_code<-as.factor(with(nsSpp3Cb, paste(Livestockdensity,Treatment,Season, sep="-")))
modTLivC <- betadisper(vare.disC, nsSpp3Cb$harvest_code,type = c("centroid"))
modTLivC

# PERMANOVA grasses
nsSpp3Cb$time_code<-as.factor(with(nsSpp3Cb, paste(Transect,Block,Treatment,Replicate, sep="-")))
PermC<-adonis(vare.disC~Livestockdensity+Treatment+Season+rainmm+
                Season:Treatment+Livestockdensity:Season+
                Livestockdensity:Treatment+rainmm:Season+
                rainmm:Livestockdensity+rainmm:Treatment+
                Treatment:Livestockdensity:Season+
                Treatment:Livestockdensity:rainmm, 
              strata=as.numeric(nsSpp3Cb$Block),#/as.numeric(nsSpp3G$time_code),
              method = "bray",perm=999, data=nsSpp3Cb)
PermC$aov.tab 
PermC$aov.tab$R2

# Extract centroids
#https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258
#vec.sp.dfC<-as.data.frame(cbind(NS.harC$factors$centroids*sqrt(NS.harC$factors$r)))
vec.sp.dfC<-data.frame(cbind(modTLivC$centroids[,1],modTLivC$centroids[,2]))
colnames(vec.sp.dfC)[1]<-"NMDS1"
colnames(vec.sp.dfC)[2]<-"NMDS2"
#vec.sp.dfC<-as.data.frame(scores(mdsNSC, display = "sites"))
vec.sp.dfC$Livestockdensity<-c("High","High","High","High","High","High","Low","Low","Low","Low","Low","Low",
                               "Medium","Medium","Medium","Medium","Medium","Medium")
vec.sp.dfC$Treatment<-c("Control","Control","Control","Exclosure","Exclosure","Exclosure","Control","Control","Control","Exclosure","Exclosure","Exclosure",
                        "Control","Control","Control","Exclosure","Exclosure","Exclosure")
vec.sp.dfC$Season<-rep(c("Long","Short I","Short II"))
vec.sp.dfC$harvest_code<-as.factor(with(vec.sp.dfC, paste(Livestockdensity,Treatment,Season, sep="-")))

nsSpp3Cb$HerbClimberBiomass1<-nsSpp3Cb$ClimberNetReharvestBiomass1+nsSpp3Cb$HerbNetReharvestBiomass1
CReharvestBio<-aggregate(HerbClimberBiomass1~harvest_code,nsSpp3Cb,mean)
vec.sp.dfC<-merge(vec.sp.dfC, CReharvestBio, by.x = "harvest_code")
vec.sp.dfC$HerbClimberBiomass1<-as.numeric(vec.sp.dfC$HerbClimberBiomass1)

#NS.harC$scores
#Herb and climber spp scores
spp.scrsC <- as.data.frame(scores(mdsNSC, display = "species"))
spp.scrsC <- cbind(spp.scrsC, Species = rownames(spp.scrsC))
spp.scrsC$Livestockdensity<-c("Low")
spp.scrsC$Season<-c("Short I")
spp.scrsC$Treatment<-c("Control")

# Reorder by season and Harvest code
vec.sp.dfC$Season<- factor(vec.sp.dfC$Season, levels = c("Short I","Long","Short II"))
vec.sp.dfC$Season2<-c(1,2,3)
vvec.sp.dfC<- vec.sp.dfW[order(vec.sp.dfC$Season2),] 
vec.sp.dfC$harvest_code<-as.factor(with(vec.sp.dfC, paste(Livestockdensity,Treatment, sep="-")))

vec.sp.dfC$Livestockdensity<-as.factor(vec.sp.dfC$Livestockdensity)
levels(vec.sp.dfC$Livestockdensity)

# Herb and climber multivariate movement plot
CenPlotC<-ggplot(vec.sp.dfC[order(vec.sp.dfC$Season),],aes(x=NMDS1,y=NMDS2, colour=Livestockdensity,fill=Livestockdensity))
#CenPlotG<-CenPlotG+geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=Livestockdensity), size=1,show.legend=T)
CenPlotC<-CenPlotC+geom_point(aes(size=HerbClimberBiomass1,shape=Treatment), stroke=1)
CenPlotC<-CenPlotC+geom_path(aes(group=harvest_code),size=1,arrow = arrow(angle=25,length = unit(3.5, "mm")), show.legend = F)
#CenPlotW<-CenPlotW+geom_text(data=spp.scrs,aes(label=Species),colour="light grey")
CenPlotC<-CenPlotC+geom_text(aes(label=Season),hjust=0, vjust=-.95, show.legend = F)
CenPlotC<-CenPlotC+scale_colour_manual(values=c("black","grey80","grey50"))
CenPlotC<-CenPlotC+scale_fill_manual(values=c("black","grey80","grey50"))
CenPlotC<-CenPlotC+scale_radius(range=c(1,8))
CenPlotC<-CenPlotC+scale_shape_manual(values=c(21,22))
CenPlotC<-CenPlotC+ggtitle("Herbs & Climbers")
CenPlotC<-CenPlotC+#theme_bw() +
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
CenPlotC

# Herb also increases...
#RhynchosiaminimaLDC
names(nsSpp2G)
RhyMin<-ggplot(nsSpp2,aes(x=Trt_Liv,y=RhynchosiaminimaLDC  ,shape=Treatment,colour=Livestockdensity))
RhyMin<-RhyMin+geom_boxplot(outlier.shape=NA,show.legend=F)+geom_jitter(size=2,stroke=1,show.legend=F)
RhyMin

# Combine multivariate plots
library(grid)
library(gridExtra)
library(egg)
library(ggpubr)

egg::ggarrange(CenPlotG+ theme(legend.position="none"),
               CenPlotC+ theme(legend.position="none"), 
               CenPlotW+ theme(legend.position="none"), ncol=3) #common.legend = T)

##########################################################################################
#### Species richness and evenness ####
# Paritioning the beta - alpha and gamma
library(lme4)
library(betapart)
##########################################################################################

# OVERALL COMMUNITY
# Shannon index, species richness and Eveness
nsSpp3$Shannon<-diversity(nsSpp3[,10:74], index = "shannon")
nsSpp3$richness<-specnumber(nsSpp3[,10:74], MARGIN = 1)
nsSpp3$eveness<-nsSpp3$Shannon/log(nsSpp3$richness)

# Species richness and biomass
xyplot(Shannon~TotalBiomass1|Livestockdensity,nsSpp3)
abline(lm(Shannon~TotalBiomass1,nsSpp3))
summary(lm(Shannon~TotalBiomass1,nsSpp3)) # Highly significant
xyplot(richness~TotalBiomass1|Livestockdensity,nsSpp3)
abline(lm(richness~TotalBiomass1,nsSpp3))
summary(lm(richness~TotalBiomass1,nsSpp3))
xyplot(eveness~TotalBiomass1|Livestockdensity,nsSpp3)
abline(lm(eveness~TotalBiomass1,nsSpp3))
summary(lm(eveness~TotalBiomass1,nsSpp3)) # NS

MyVar<-c("Shannon","richness","eveness")
pairs(nsSpp3[,MyVar],lower.panel = panel.cor) # All highligh correlated
#  Richness and eveness less so...Higher diversity = higher even
# Less diverse, less even = more productive

# Drop Short I
nsSpp3i<-droplevels(nsSpp3[nsSpp3$Date!="21.10.2012",])

# Shannon diversity
ShanTot<-lmer(Shannon~Livestockdensity+#Treatment+
               #Livestockdensity:Season+Season:Treatment+
               #Livestockdensity:Treatment+ 
               #Treatment:Livestockdensity:Season+
               (1 |Block), data=nsSpp3i)
summary(ShanTot)
anova(ShanTot) 
AIC(ShanTot) # 184.9415
drop1(ShanTot, test="Chi") # Only livestock

#                 Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#Livestockdensity 1.0024 0.50122     2   177   3.415 0.03506 *

# Update model
ShanTota <- update(ShanTot, .~. -Livestockdensity)
anova(ShanTot,ShanTota)
#         Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#ShanTot   5 172.39 188.35 -81.195   162.39 6.8152      2    0.03312 *

# Resid versus fit
E1 <- resid(ShanTot,type="pearson") 
F1 <- fitted(ShanTot)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1,  y = E1, xlab = "Fitted values",ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

# lsmeans
library(multcomp)
library(multcompView)
library(lsmeans)
library(lmerTest)
library(Hmisc)
library(pbkrtest)
boxplot(Shannon~Livestockdensity,nsSpp3i) # High has higher diversity
summary(glht(ShanTot, mcp(Livestockdensity="Tukey")))  # High has higher diversity - marginally

# Eveness
EveTot<-lmer(eveness~Livestockdensity+Season+rainmm+#Treatment+
                #Livestockdensity:Season+Season:Treatment+
                #Livestockdensity:Treatment+   rainmm:Season+ 
                #rainmm:Livestockdensity+rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3i)
summary(EveTot)
anova(EveTot) 
AIC(EveTot) # -123.8387
drop1(EveTot, test="Chi") # ALL NS
bwplot(eveness~Livestockdensity|Season,nsSpp3)

# Update model
EveTota <- update(EveTot, .~. -Livestockdensity)
EveTotb <- update(EveTot, .~. -Season)
EveTotc <- update(EveTot, .~. -rainmm)
anova(EveTot,EveTota)
anova(EveTot,EveTotb)
anova(EveTot,EveTotc)
#         Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#EveTot   7 -157.48 -135.13 85.740  -171.48 7.5263      2    0.02321 * # Livestockdensity
#EveTot   7 -157.48 -135.13 85.740  -171.48 5.171      1    0.02297 * # Season
#EveTot   7 -157.48 -135.13 85.740  -171.48 5.3941      1     0.0202 * # Rain

# Contrasts within the model
boxplot(eveness~Livestockdensity,nsSpp3i) # Low livestock has lower eveness...
boxplot(eveness~Season,nsSpp3i) # Long is more even than the rest
summary(glht(EveTot, mcp(Livestockdensity="Tukey")))  # Lower evenenss with low livestock 
summary(glht(EveTot, mcp(Season="Tukey"))) 

#### OVERALL COMMUNITY SPECIES CHANGE ####
# Create a factor to define each plot shared across seasons
nsSpp3$plot_codeII<-as.factor(with(nsSpp3, paste(Transect,Block,Trtname,Replicate, sep="_")))

#### Community diversity based on abundance ####
nsSpp3Low<-droplevels(nsSpp3[nsSpp3$Livestockdensity=="Low",])
nsSpp3Med<-droplevels(nsSpp3[nsSpp3$Livestockdensity=="Medium",])
nsSpp3High<-droplevels(nsSpp3[nsSpp3$Livestockdensity=="High",])

OverallPresabsL<-ifelse(nsSpp3Low[,10:74]>0,1,0)
OverallPresabsM<-ifelse(nsSpp3Med[,10:74]>0,1,0)
OverallPresabsH<-ifelse(nsSpp3High[,10:74]>0,1,0)

beta.sample(OverallPresabsL, index.family="sorensen", sites=nrow(nsSpp3Low), samples = 1)
beta.sample(OverallPresabsM, index.family="sorensen", sites=nrow(nsSpp3Med), samples = 1)
beta.sample(OverallPresabsH, index.family="sorensen", sites=nrow(nsSpp3High), samples = 1)

#### Compare diversity through time of overall community ####
nsSpp3SeasonI<-droplevels(nsSpp3[nsSpp3$Date=="21.10.2012",])
nsSpp3SeasonII<-droplevels(nsSpp3[nsSpp3$Date=="21.11.2013",])
nsSpp3SeasonIII<-droplevels(nsSpp3[nsSpp3$Date=="21.6.2013",])

# Presence and absence of species
names(nsSpp3SeasonI[,10:74])
names(nsSpp3SeasonII[,10:74])
names(nsSpp3SeasonIII[,10:74])
PresabsI<-ifelse(nsSpp3SeasonI[,10:74]>0,1,0)
PresabsII<-ifelse(nsSpp3SeasonII[,10:74]>0,1,0)
PresabsIII<-ifelse(nsSpp3SeasonIII[,10:74]>0,1,0)

# Convert to betapart objects
PresabsI.core <- betapart.core(PresabsI)
PresabsII.core <- betapart.core(PresabsII)
PresabsIII.core <- betapart.core(PresabsIII)

# Assign plot (without season) to each row
row.names(PresabsI) <- paste(nsSpp3SeasonI$plot_codeII, 1:nrow(PresabsI), sep="")
row.names(PresabsII) <- paste(nsSpp3SeasonII$plot_codeII, 1:nrow(PresabsII), sep="")
row.names(PresabsIII) <- paste(nsSpp3SeasonIII$plot_codeII, 1:nrow(PresabsIII), sep="")

dimnames(PresabsI) 
dimnames(PresabsII)
dimnames(PresabsIII)

#### Beta.temp - Overall community ####
ObtI<-beta.temp(PresabsI,PresabsII,index.family="sor") # this is soreson
ObtII<-beta.temp(PresabsII,PresabsIII,index.family="sor") 

nsSpp3S2_3<-rbind(nsSpp3SeasonII,nsSpp3SeasonIII)
betaSeason2_3<-rbind(ObtI,ObtII)
nsSpp3beta<-cbind(nsSpp3S2_3,betaSeason2_3)

# Analyse Beta diversiy (turnovr) - OVERALL COMMUNITY
turnSor<-lmer(beta.sor~Livestockdensity+Treatment+Season+
               Livestockdensity:Season+#Season:Treatment+
               Livestockdensity:Treatment+
               #Treatment:Livestockdensity:Season+
               (1 |Block), data=nsSpp3beta)

summary(turnSor)
anova(turnSor) 
AIC(turnSor) #-148.959
drop1(turnSor, test="Chi") 

# Resid versus fit
E1 <- resid(turnSor,type="pearson") 
F1 <- fitted(turnSor)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1,  y = E1, xlab = "Fitted values",ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

# Update and remove factors # issues with interactions
turnSor3<-lmer(beta.sor~Livestockdensity+Treatment+Season+ (1 |Block), data=nsSpp3beta)
turnSora <- update(turnSor, .~. -Livestockdensity:Treatment)
turnSorb <- update(turnSor, .~. -Livestockdensity:Season)
turnSor3b<- update(turnSor3, .~. -Season)
turnSor3c <- update(turnSor3, .~. -Treatment)
turnSor3d <- update(turnSor3, .~. -Livestockdensity)


anova(turnSor,turnSora)
anova(turnSor,turnSorb)
anova(turnSor3,turnSor3b)
anova(turnSor3,turnSor3c)
anova(turnSor3,turnSor3d)

#Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
#turnSor  11 -195.66 -160.54 108.83  -217.66 15.579      2  0.0004141 *** #Livestockdensity:Treatment
#turnSor  11 -195.66 -160.54 108.83  -217.66 6.6669      2    0.03567 * #Livestockdensity:Season
#turnSor3   7 -181.97 -159.62 97.985  -195.97 5.8149      1    0.01589 * #Season
#turnSor3   7 -181.97 -159.62 97.985  -195.97 7.5699      1   0.005935 ** #Treatment
#turnSor3   7 -181.97 -159.62 97.985  -195.97 5.3486      2    0.06896 . #Livestockdensity

# Analyse turnover component (sim) - OVERALL COMMUNITY
turnbt<-lmer(beta.sim~Livestockdensity+Treatment+Season+rainmm+
                Livestockdensity:Season+Season:Treatment+
                Livestockdensity:Treatment+ #  rainmm:Season+ 
                rainmm:Livestockdensity+rainmm:Treatment+
                Treatment:Livestockdensity:Season+
                Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3beta)
summary(turnbt)
anova(turnbt) 
AIC(turnbt) # -55.17815
drop1(turnbt, test="Chi") # ALL NS

# Resid versus fit
E1 <- resid(turnbt,type="pearson") 
F1 <- fitted(turnbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1,  y = E1, xlab = "Fitted values",ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

# Update and remove factors # issues with interactions
turnbt3<-lmer(beta.sim~Livestockdensity+Treatment+Season+rainmm+ (1 |Block), data=nsSpp3beta)
turnbta <- update(turnbt, .~. -Treatment:Livestockdensity:rainmm)
turnbtb <- update(turnbt, .~. -Treatment:Livestockdensity:Season)
turnbt2 <- update(turnbta, .~. -Treatment:Livestockdensity:Season) # No three way interactions..
turnbt2a <- update(turnbt2, .~. -rainmm:Treatment)
turnbt2b <- update(turnbt2, .~. -rainmm:Livestockdensity)
turnbt2c <- update(turnbt2, .~. -Livestockdensity:Treatment)
turnbt2d <- update(turnbt2, .~. -Season:Treatment)
turnbt2e <- update(turnbt2, .~. -Livestockdensity:Season)
turnbt3a <- update(turnbt3, .~. -rainmm)
turnbt3b<- update(turnbt3, .~. -Season)
turnbt3c <- update(turnbt3, .~. -Treatment)
turnbt3d <- update(turnbt3, .~. -Livestockdensity)

anova(turnbt,turnbta)
anova(turnbt,turnbtb)
anova(turnbt2,turnbt2a)
anova(turnbt2,turnbt2b)
anova(turnbt2,turnbt2c)
anova(turnbt2,turnbt2d)
anova(turnbt2,turnbt2e)
anova(turnbt3,turnbt3a)
anova(turnbt3,turnbt3b)
anova(turnbt3,turnbt3c)
anova(turnbt3,turnbt3d)

#Df    AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)   
#turnbta 18 -148.5 -91.023 92.248   -184.5 
#turnbt  20 -155.9 -92.040 97.949   -195.9 11.402      2   0.003342 ** #Treatment:Livestockdensity:rainmm
#turnbt  20 -155.90 -92.040 97.949  -195.90 11.532      2   0.003132 ** #Treatment:Livestockdensity:Season
#turnbt2  16 -149.34  -98.253 90.670  -181.34 0.7724      1     0.3795 #rainmm:Treatment
#turnbt2  16 -149.34  -98.253  90.67  -181.34 6.3596      2    0.04159 * #rainmm:Livestockdensity
#turnbt2  16 -149.34 -98.253 90.670  -181.34 15.957      2  0.0003428 ***#Livestockdensity:Treatment
#turnbt2  16 -149.34  -98.253 90.670  -181.34 0.7924      1     0.3734 #Season:Treatment
#turnbt2  16 -149.34  -98.253 90.670  -181.34 6.5259      2    0.03828 *#Livestockdensity:Season
#turnbt3   8 -141.61 -116.07 78.807  -157.61 6.4437      1    0.01113 * # rain
#turnbt3   8 -141.61 -116.07 78.807  -157.61 6.2094      1    0.01271 * # Season
#turnbt3   8 -141.61 -116.07 78.807  -157.61 10.039      1   0.001533 ** # Treatment
#turnbt3   8 -141.61 -116.07 78.807  -157.61 7.4406      2    0.02423 * # Livestock density

nsSpp3betaMean<-aggregate(beta.sim~Livestockdensity+Treatment+Season,nsSpp3beta,mean)
ggplot(nsSpp3betaMean,aes(x=Livestockdensity,y=beta.sim,colour=Treatment, shape=Season))+
  geom_point(size=4.5, position=pd)
nsSpp3betaRain<-aggregate(beta.sim~Livestockdensity+Treatment+Season+rainmm,nsSpp3beta,mean)
ggplot(nsSpp3betaRain,aes(x=rainmm,y=beta.sim,colour=Treatment, shape=Livestockdensity))+
  geom_point(size=4.5)

# Analyse nested component (sne) - OVERALL
nestbt<-lmer(beta.sne~Livestockdensity+Treatment+Season+rainmm+
               #Livestockdensity:Season+Season:Treatment+
               #Livestockdensity:Treatment+rainmm:Season+ 
               #rainmm:Livestockdensity+rainmm:Treatment+
               #Treatment:Livestockdensity:Season+
               #Treatment:Livestockdensity:rainmm+
               (1 |Block), data=nsSpp3beta)
summary(nestbt)
anova(nestbt) 
AIC(nestbt) # -235.313
drop1(nestbt, test="Chi") # ALL NS

# Resid versus fit
E1 <- resid(turnbt,type="pearson") 
F1 <- fitted(turnbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1,  y = E1, xlab = "Fitted values",ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

#### Graphing Shannon, evenness and species turnover for the overall community

library(betapart)

# Create matrics
nsSpp3$Shannon<-diversity(nsSpp3[,10:74], index = "shannon")
nsSpp3$richness<-specnumber(nsSpp3[,10:74], MARGIN = 1)
nsSpp3$eveness<-nsSpp3$Shannon/log(nsSpp3$richness)

#### Compare diversity through time of overall community ####
nsSpp3SeasonI<-droplevels(nsSpp3[nsSpp3$Date=="21.10.2012",])
nsSpp3SeasonII<-droplevels(nsSpp3[nsSpp3$Date=="21.11.2013",])
nsSpp3SeasonIII<-droplevels(nsSpp3[nsSpp3$Date=="21.6.2013",])

# Presence and absence of species
PresabsI<-ifelse(nsSpp3SeasonI[,10:74]>0,1,0)
PresabsII<-ifelse(nsSpp3SeasonII[,10:74]>0,1,0)
PresabsIII<-ifelse(nsSpp3SeasonIII[,10:74]>0,1,0)

# Assign plot (without season) to each row
row.names(PresabsI) <- paste(nsSpp3SeasonI$plot_codeII, 1:nrow(PresabsI), sep="")
row.names(PresabsII) <- paste(nsSpp3SeasonII$plot_codeII, 1:nrow(PresabsII), sep="")
row.names(PresabsIII) <- paste(nsSpp3SeasonIII$plot_codeII, 1:nrow(PresabsIII), sep="")

#### Beta.temp - Overall community ####
ObtI<-beta.temp(PresabsI,PresabsII,index.family="sor") # this is soreson
ObtII<-beta.temp(PresabsII,PresabsIII,index.family="sor") 

nsSpp3S2_3<-rbind(nsSpp3SeasonII,nsSpp3SeasonIII)
betaSeason2_3<-rbind(ObtI,ObtII)
nsSpp3beta<-cbind(nsSpp3S2_3,betaSeason2_3)

# Remove short season
nsSpp3ii<-droplevels(nsSpp3[nsSpp3$Season!="Short I",])

# Summary
ShanXsum<-aggregate(Shannon~Livestockdensity,nsSpp3,mean)
EvenXsum<-aggregate(eveness~Livestockdensity,nsSpp3,mean)
BetaSXsum<-aggregate(beta.sor~Livestockdensity,nsSpp3beta,mean)

# Means + SE Shannon, evenness and species turnover 
# Means
ShanX<-aggregate(Shannon~Livestockdensity+Treatment+Season,nsSpp3ii,mean)
EvenX<-aggregate(eveness~Livestockdensity+Treatment+Season,nsSpp3ii,mean)
BetaSX<-aggregate(beta.sor~Livestockdensity+Treatment+Season,nsSpp3beta,mean)

# Errors bars
sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
ShanSEM<-aggregate(Shannon~Livestockdensity+Treatment+Season,nsSpp3ii,sem)
EvenSEM<-aggregate(eveness~Livestockdensity+Treatment+Season,nsSpp3ii,sem)
BetaSSEM<-aggregate(beta.sor~Livestockdensity+Treatment+Season,nsSpp3beta,sem)
#BetaSSEM<-aggregate(abs(Bcdist2)~Livestockdensity+Treatment+Season,reharv23,sem)

# Add SE terms to average summaries
ShanX$SEM<-ShanSEM$Shannon
EvenX$SEM<-EvenSEM$eveness
BetaSX$SEM<-BetaSSEM$beta.sor

# Create grouping variable - important for connecting points by a line
ShanX$Groups<-as.factor(with(ShanX, paste(Livestockdensity,Treatment,sep="_")))
EvenX$Groups<-as.factor(with(EvenX, paste(Livestockdensity,Treatment,sep="_")))
BetaSX$Groups<-as.factor(with(BetaSX, paste(Livestockdensity,Treatment,sep="_")))

# Reorder livestock treatment so low - medium - high (alphabetical high, low, medium)
ShanX$Livestockdensity<- factor(ShanX$Livestockdensity, levels = c("Low","Medium","High"))
EvenX$Livestockdensity<- factor(EvenX$Livestockdensity, levels = c("Low","Medium","High"))
BetaSX$Livestockdensity<- factor(BetaSX$Livestockdensity, levels = c("Low","Medium","High"))

# Rename Treatment names - Control = Open 
levels(ShanX$Treatment)<-c("Open","Exclosed")
levels(EvenX$Treatment)<-c("Open","Exclosed")
levels(BetaSX$Treatment)<-c("Open","Exclosed")

# Fill livestock and treatment control - this will mean we need to override the legend (see below)
ShanX$LivTrt<-as.factor(with(ShanX, paste(Livestockdensity, Treatment, sep="")))
EvenX$LivTrt<-as.factor(with(EvenX, paste(Livestockdensity, Treatment, sep="")))
BetaSX$LivTrt<-as.factor(with(BetaSX, paste(Livestockdensity, Treatment, sep="")))

# Plotting Shannon Diversity graph
pd <- position_dodge(0.5) # Dodge term - there are two options here either in ggplot arguement or seperate 
ShanP<-ggplot(ShanX, aes(x=Season, y=Shannon,colour=Livestockdensity,shape=Livestockdensity, fill=LivTrt, group=Groups, alpha=Treatment)) 
ShanP<-ShanP+geom_errorbar(aes(x = Season, ymin=Shannon-SEM,ymax=Shannon+SEM),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
ShanP<-ShanP+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
ShanP<-ShanP+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
ShanP<-ShanP+scale_shape_manual(values=c(21,24,22))
ShanP<-ShanP+scale_alpha_manual(values=c(1,1))
ShanP<-ShanP+scale_colour_manual(values=c("grey70","grey35","black"))
ShanP<-ShanP+scale_fill_manual(values=c("white","black","white","grey70","white","grey35"))
ShanP<-ShanP+ylab("Shannon Diversity")+xlab("")
ShanP<-ShanP+#theme_bw() +
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
        ,legend.background = element_rect(fill="transparent",colour=NA)
        ,legend.position = "right"
        ,legend.justification = "top"
        ,legend.direction="vertical"
        ,legend.key.width = unit(1.2,"cm"))
# Here we override the legend, it needs to be legend= T to work!
ShanP<- ShanP +guides(fill=F,shape=F, 
colour = guide_legend("Livestock density",override.aes = list(shape=c(21,24,22),size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1,linetype=NA)),
alpha = guide_legend(override.aes = list(shape=c(21,21),size=3.5,fill=c("white","grey50"),col=c("grey 50","grey 50"), stroke=1,linetype=NA)))
ShanP

# Plotting Eveness graph
pd <- position_dodge(0.5) # Dodge term - there are two options here either in ggplot arguement or seperate 
EvenP<-ggplot(EvenX, aes(x=Season, y=eveness,colour=Livestockdensity,shape=Livestockdensity, fill=LivTrt, group=Groups, alpha=Treatment)) 
EvenP<-EvenP+geom_errorbar(aes(x = Season, ymin=eveness-SEM,ymax=eveness+SEM),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
EvenP<-EvenP+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
EvenP<-EvenP+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
EvenP<-EvenP+scale_shape_manual(values=c(21,24,22))
EvenP<-EvenP+scale_alpha_manual(values=c(1,1))
EvenP<-EvenP+scale_colour_manual(values=c("grey70","grey35","black"))
EvenP<-EvenP+scale_fill_manual(values=c("white","black","white","grey70","white","grey35"))
EvenP<-EvenP+ylab("Evenness")
EvenP<-EvenP+#theme_bw() +
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
        ,legend.background = element_rect(fill="transparent",colour=NA)
        ,legend.position = "right"
        ,legend.justification = "top"
        ,legend.direction="vertical"
        ,legend.key.width = unit(1.2,"cm"))
# Here we override the legend, it needs to be legend= T to work!
EvenP<- EvenP +guides(fill=F,shape=F, 
                      colour = guide_legend("Livestock density",override.aes = list(shape=c(21,24,22),size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1,linetype=NA)),
                      alpha = guide_legend(override.aes = list(shape=c(21,21),size=3.5,fill=c("white","grey50"),col=c("grey 50","grey 50"), stroke=1,linetype=NA)))
EvenP

# Plotting Beta diversity graph
pd <- position_dodge(0.5) # Dodge term - there are two options here either in ggplot arguement or seperate 
BetaP<-ggplot(BetaSX, aes(x=Season, y=beta.sor,colour=Livestockdensity,shape=Livestockdensity, fill=LivTrt, group=Groups, alpha=Treatment)) 
BetaP<-BetaP+geom_errorbar(aes(x = Season, ymin=beta.sor-SEM,ymax=beta.sor+SEM),position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F)
BetaP<-BetaP+geom_line(position=pd,stat = "identity",size=.75,show.legend = T) 
BetaP<-BetaP+geom_point(position=pd,stat = "identity",size=3.5, stroke=1)
BetaP<-BetaP+scale_shape_manual(values=c(21,24,22))
BetaP<-BetaP+scale_alpha_manual(values=c(1,1))
BetaP<-BetaP+scale_colour_manual(values=c("grey70","grey35","black"))
BetaP<-BetaP+scale_fill_manual(values=c("white","black","white","grey70","white","grey35"))
BetaP<-BetaP+ylab((expression(italic(beta)~"- diversity")))+xlab("")
BetaP<-BetaP+#theme_bw() +
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
        ,legend.background = element_rect(fill="transparent",colour=NA)
        ,legend.position = "right"
        ,legend.justification = "top"
        ,legend.direction="vertical"
        ,legend.key.width = unit(1.2,"cm"))
# Here we override the legend, it needs to be legend= T to work!
BetaP
BetaP<- BetaP +guides(fill=F,shape=F, 
                      colour = guide_legend("Livestock density",override.aes = list(shape=c(21,24,22),size=3.5,fill=c("grey70","grey35","black"),col=c("grey70","grey35","black"), stroke=1,linetype=NA)),
                      alpha = guide_legend(override.aes = list(shape=c(21,21),size=3.5,fill=c("white","grey50"),col=c("grey 50","grey 50"), stroke=1,linetype=NA)))
BetaP

# Combine plots into panel
library(grid)
library(gridExtra)
library(egg)

# Extra legend from legend plot
library(ggpubr)
legend <- get_legend(ShanP)
arrangeGrob(legend)

p3 <- arrangeGrob(legend,p,dp,bp, ncol=3, nrow=2,widths=c(1,1,1), heights=c(1,1.4),layout_matrix = cbind(c(1,4), c(2,4),c(3,4))) #common.legend = T)
grid.arrange(p3)

egg::ggarrange(ShanP+ theme(legend.position="none"),EvenP+ theme(legend.position="none"),BetaP, ncol=3) #common.legend = T)

filename <- paste0("/Users/anotherswsmith/Documents/AfricanBioServices/Colloborators/Desalegn Wana /Pastoralism_in_NecSar/Pastoralism_in_NecSar/", "ShanEvenBeta", "_",Sys.Date(), ".jpeg" )
jpeg (filename, width=26.5, height=10, res=400, unit="cm")
p3
dev.off()

###################################################################################
##### Grasses #####
# Shannon index, species richness and Eveness
nsSpp3G$Shannon<-diversity(nsSpp3G[10:28], index = "shannon")
nsSpp3G$richness<-specnumber(nsSpp3G[10:28], MARGIN = 1)
nsSpp3G$eveness<-nsSpp3G$Shannon/log(nsSpp3G$richness)

# Grass  richness and biomass
plot(Shannon~TotalBiomass1,nsSpp3G)
abline(lm(Shannon~TotalBiomass1,nsSpp3G))
summary(lm(Shannon~TotalBiomass1,nsSpp3G)) # Highly significant
plot(richness~TotalBiomass1,nsSpp3G)
abline(lm(richness~TotalBiomass1,nsSpp3G))
summary(lm(richness~TotalBiomass1,nsSpp3G)) # NS
plot(eveness~TotalBiomass1,nsSpp3G)
abline(lm(eveness~TotalBiomass1,nsSpp3G))
summary(lm(eveness~TotalBiomass1,nsSpp3G)) # NS

# Drop Short I
nsSpp3Gi<-droplevels(nsSpp3G[nsSpp3G$Date!="21.10.2012",])

# Shannon diversity
ShanG<-lmer(Shannon~Livestockdensity+#Season+#Treatment+
               # Livestockdensity:Season+#Season:Treatment+
               # Livestockdensity:Treatment+  
                #Treatment:Livestockdensity:Season+
                (1 |Block), data=nsSpp3Gi)
summary(ShanG)
anova(ShanG) 
AIC(ShanG) #  177.3046
drop1(ShanG, test="Chi") # ALL NS
boxplot(Shannon~Livestockdensity,nsSpp3Gi) # Livestock lower diversity

# Update model
ShanGc1 <- update(ShanG, .~. -Livestockdensity)

anova(ShanG,ShanGc1)

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#ShanG    5 134.58 150.55 -62.291   124.58 7.4068      2    0.02464 ** #Livestock

# lsmeans
library(multcomp)
library(multcompView)
library(lsmeans)
library(lmerTest)
library(Hmisc)
library(pbkrtest)
boxplot(Shannon~Livestockdensity,nsSpp3Gi) # High has higher diversity
summary(glht(ShanG, mcp(Livestockdensity="Tukey")))  # High has higher diversity
lsShanG <-difflsmeans(ShanG, test.effs= "Livestockdensity:Season" )
lsShanG 

aggregate(Shannon~Livestockdensity+Season,nsSpp3Gi,mean)

with(nsSpp3Gi, {interaction.plot(Season,Livestockdensity,Shannon,
                                  xlab = "tree canopy",
                                  ylab = "Biomass Rainfall correlation",
                                  fun=mean)})
names(SerEbio5)

# Eveness of grasses
EveG<-lmer(eveness~Livestockdensity+Season+#rainmm+Treatment+
               #Livestockdensity:Season+Season:Treatment+
               #Livestockdensity:Treatment+   rainmm:Season+ 
               #rainmm:Livestockdensity+rainmm:Treatment+
               #Treatment:Livestockdensity:Season+
               #Treatment:Livestockdensity:rainmm+
               (1 |Block), data=nsSpp3G)
summary(EveG)
anova(EveG) 
AIC(EveG) # -30.63383
drop1(EveG, test="Chi") # ALL NS
bwplot(eveness~Livestockdensity|Season,nsSpp3G)

# Update model
EveGa <- update(EveG, .~. -Livestockdensity)
EveGb <- update(EveG, .~. -Season)

anova(EveG,EveGa)
anova(EveG,EveGb)

#      Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq) 
#EveG   7 -58.075 -33.286 36.037  -72.075 8.2332      2     0.0163 * #Livestockdensity
#EveG   7 -58.075 -33.286 36.037  -72.075 12.838      2   0.001631 ** #Season

# Model contrasts
boxplot(eveness~Livestockdensity,nsSpp3G)
boxplot(eveness~Season,nsSpp3G) 
summary(glht(EveG, mcp(Livestockdensity="Tukey")))  # Eveness higher in high livestock
summary(glht(EveG, mcp(Season="Tukey"))) # Increased from short I


# Analyse nested component (sne) - Grasses
nestGbt<-lmer(beta.sim~Livestockdensity+Treatment+Season+rainmm+
                #Livestockdensity:Season+Season:Treatment+
                Livestockdensity:Treatment+#rainmm:Season+
                #rainmm:Livestockdensity+rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3Gbeta)
summary(nestGbt)
anova(nestGbt) #
AIC(nestGbt) # -190.4819
drop1(nestGbt, test="Chi") # 
#Livestockdensity:Treatment 0.02193 *

# Update and remove factors # issues with interactions
nestGbta <- update(nestGbt, .~. -Livestockdensity:Treatment)
nestGbtb <- update(nestGbt, .~. -rainmm)
nestGbtc <- update(nestGbt, .~. -Season)
nestGbta2 <- update(nestGbta, .~. -Treatment)
nestGbta3 <- update(nestGbta, .~. -Livestockdensity)

anova(nestGbt,nestGbta)
anova(nestGbt,nestGbtb)
anova(nestGbt,nestGbtc)
anova(nestGbta,nestGbta2)
anova(nestGbta,nestGbta3)

#         Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
#nestGbt  10 -225.96 -194.03 122.98  -245.96 7.8661      2    0.01958 * #Livestockdensity:Treatment
#nestGbt  10 -225.96 -194.03 122.98  -245.96 7.3234      1   0.006806 ** # rain
#nestGbt  10 -225.96 -194.03 122.98  -245.96 7.3733      1    0.00662 ** # Season
#nestGbta   8 -222.10 -196.55 119.05  -238.10 1.261      1     0.2615 # Treatment
#nestGbta   8 -222.10 -196.55 119.05  -238.10 1.0472      2     0.5924 # Livestock

##### GRASSES #####
names(nsSpp3G)
#### Community diversity based on abundance ####
nsSpp3GLow<-droplevels(nsSpp3G[nsSpp3G$Livestockdensity=="Low",])
nsSpp3GMed<-droplevels(nsSpp3G[nsSpp3G$Livestockdensity=="Medium",])
nsSpp3GHigh<-droplevels(nsSpp3G[nsSpp3G$Livestockdensity=="High",])

OverallPresabsGL<-ifelse(nsSpp3GLow[,10:28]>0,1,0)
OverallPresabsGM<-ifelse(nsSpp3GMed[,10:28]>0,1,0)
OverallPresabsGH<-ifelse(nsSpp3GHigh[,10:28]>0,1,0)

beta.sample(OverallPresabsGL, index.family="sorensen", sites=nrow(nsSpp3GLow), samples = 1)
beta.sample(OverallPresabsGM, index.family="sorensen", sites=nrow(nsSpp3GMed), samples = 1)
beta.sample(OverallPresabsGH, index.family="sorensen", sites=nrow(nsSpp3GHigh), samples = 1)
# Low # 0.96444737
# Med # 0.93848338 # Lower diversity in medium - grasses
# High # 0.9558174

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

# Analyse turnover component (sor) - GRASSES
turnGSor<-lmer(beta.sor~Livestockdensity+Treatment+#Season+
                #Livestockdensity:Season+Season:Treatment+
                Livestockdensity:Treatment+
                #Treatment:Livestockdensity:Season+
                (1 |Block), data=nsSpp3Gbeta)
summary(turnGSor)
anova(turnGSor) 
AIC(turnGbt) # -55.17815
drop1(turnGSor, test="Chi") # ALL NS

# Update model
turnGSora <- update(turnGSor, .~. -Livestockdensity:Treatment)
turnGSora2 <- update(turnGSora, .~. -Livestockdensity)
turnGSora3 <- update(turnGSora, .~. -Treatment)

anova(turnGSor,turnGSora)
anova(turnGSora,turnGSora2)
anova(turnGSora,turnGSora3)

# Grass beta diversity
#          Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
#turnGSor   8 -76.635 -51.091 46.317  -92.635 6.9609      2    0.03079 * #Livestockdensity:Treatment
#turnGSora   6 -73.674 -54.516 42.837  -85.674 4.8194      2    0.08984 .#Livestockdensity
#turnGSora   6 -73.674 -54.516 42.837  -85.674 0.7815      1     0.3767 # Treatment

# Analyse turnover component (sim) - GRASSES
turnGbt<-lmer(beta.sim~Livestockdensity+Treatment+Season+rainmm+
                #Livestockdensity:Season+Season:Treatment+
               rainmm:Season+# Livestockdensity:Treatment+
                #rainmm:Livestockdensity+rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
               (1 |Block), data=nsSpp3Gbeta)
summary(turnGbt)
anova(turnGbt) 
AIC(turnGbt) # -55.17815
drop1(turnGbt, test="Chi") # ALL NS

# Analyse nested component (sne) - Herbs and climbers
nestGbt<-lmer(beta.sne~Livestockdensity+Treatment+Season+rainmm+
                #Livestockdensity:Season+Season:Treatment+
                 Livestockdensity:Treatment+#rainmm:Season+
                #rainmm:Livestockdensity+rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3Gbeta)
summary(nestGbt)
anova(nestGbt) #
AIC(nestGbt) # -175.3729
drop1(nestGbt, test="Chi") # 
#Livestockdensity:Treatment 0.02193 *

# Resid versus fit
E1 <- resid(nestGbt,type="pearson") #THIS IS FOR lme..NOT lme4
F1 <- fitted(nestGbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = F1, 
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)
# Update and remove factors # issues with interactions
nestGbta <- update(nestGbt, .~. -Livestockdensity:Treatment)
nestGbtb <- update(nestGbt, .~. -rainmm)
nestGbtc <- update(nestGbt, .~. -Season)
nestGbta2 <- update(nestGbta, .~. -Treatment)
nestGbta3 <- update(nestGbta, .~. -Livestockdensity)

anova(nestGbt,nestGbta)
anova(nestGbt,nestGbtb)
anova(nestGbt,nestGbtc)
anova(nestGbta,nestGbta2)
anova(nestGbta,nestGbta3)

#         Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
#nestGbt  10 -225.96 -194.03 122.98  -245.96 7.8661      2    0.01958 * #Livestockdensity:Treatment
#nestGbt  10 -225.96 -194.03 122.98  -245.96 7.3234      1   0.006806 ** # rain
#nestGbt  10 -225.96 -194.03 122.98  -245.96 7.3733      1    0.00662 ** # Season
#nestGbta   8 -222.10 -196.55 119.05  -238.10 1.261      1     0.2615 # Treatment
#nestGbta   8 -222.10 -196.55 119.05  -238.10 1.0472      2     0.5924 # Livestock

# Species graminoids number throughtime
nsSpp3Gbeta$Richness<-apply(nsSpp3Gbeta[,10:28]>0,1,sum)
GRich<-aggregate(Richness~Livestockdensity+Treatment+Season,nsSpp3Gbeta,mean)
GRich$Richness<-round(GRich$Richness, digits = 1)

# Exploring SNE nested component graphically
sem<-function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
nsSpp3GbetaMean<-aggregate(beta.sne~Livestockdensity+Treatment+Season,nsSpp3Gbeta,mean)
nsSpp3GbetaSD<-aggregate(beta.sne~Livestockdensity+Treatment+Season,nsSpp3Gbeta,sem)
nsSpp3GbetaMean<-cbind(nsSpp3GbetaMean,nsSpp3GbetaSD[4])
colnames(nsSpp3GbetaMean)[5]<-"sd"
pd <- position_dodge(0.5)
ggplot(nsSpp3GbetaMean,aes(x=Livestockdensity,y=sqrt(beta.sne),colour=Treatment, shape=Season))+
  geom_point(size=4.5, position=pd)+geom_errorbar(aes(x = Livestockdensity, ymin=sqrt(beta.sne)-sd,ymax=sqrt(beta.sne)+sd),
      position=pd,stat = "identity",linetype="solid",width=.2,show.legend=F) +
  geom_text(label = c(GRich$Richness), size=3.5,color="black",stat = "identity",position=pd,show.legend = F)+
  ylab(expression(paste((sqrt("beta[sne]")))))

# Plot Beta by livestock density and treatment
ggplot(nsSpp3Gbeta,aes(x=beta.sor,colour=Livestockdensity, linetype=Treatment))+geom_density()
ggplot(nsSpp3Gbeta,aes(x=beta.sim,colour=Livestockdensity, linetype=Treatment))+geom_density()
ggplot(nsSpp3Gbeta,aes(x=beta.sne,colour=Livestockdensity, linetype=Treatment))+geom_density()

# Comparison of the square root transformed bsim and bsne components of bsor b
with(nsSpp3Gbeta, plot(sqrt(beta.sim) ~ sqrt(beta.sne), 
                 type='n', ylab=expression(sqrt(beta[sim])), 
                 xlab=expression(sqrt(beta[sne]))))
with(nsSpp3Gbeta, points(y= sqrt(beta.sim), x=sqrt(beta.sne), pch=c(Treatment),col=c(Livestockdensity)))

# Nested - seems to be two plots that are outliers - are these low livestock?
plot(hclust(dist(nsSpp3Gbeta$beta.sim), method="single"),col=c(nsSpp3Gbeta$Livestockdensity),
     hang=-1, main='', sub='', xlab='')

names(nsSpp3Gbeta)
aggregate(beta.sne~Livestockdensity,nsSpp3Gbeta,median)
#  Livestockdensity   beta.sne
#1             High 0.10982804
#2              Low 0.10257937
#3           Medium 0.09832011


##### Woody plants #####

# Shannon index, species richness and Eveness
nsSpp3W$Shannon<-diversity(nsSpp3W[,10:19], index = "shannon")
nsSpp3W$richness<-specnumber(nsSpp3W[,10:19], MARGIN = 1)
nsSpp3W$eveness<-nsSpp3W$Shannon/log(nsSpp3W$richness)

aggregate(Shannon~Livestockdensity+Treatment,nsSpp3W,mean)

# Drop first date
nsSpp3Wi<-droplevels(nsSpp3W[nsSpp3W$Date!="21.10.2012",])

# Shannon diversity
ShanW<-lmer(Shannon~Livestockdensity+Season+rainmm+Treatment+
              #Livestockdensity:Season+Season:Treatment+
              Livestockdensity:Treatment+  # rainmm:Season+ 
              #rainmm:Livestockdensity+rainmm:Treatment+
              #Treatment:Livestockdensity:Season+
              #Treatment:Livestockdensity:rainmm+
              (1 |Block), data=nsSpp3Wi)
summary(ShanW)
anova(ShanW) 
AIC(ShanW) #  177.3046
drop1(ShanW, test="Chi") # ALL NS

# Update model
ShanWa <- update(ShanW, .~. -Livestockdensity:Treatment)
ShanWb <- update(ShanW, .~. -rainmm)
ShanWc <- update(ShanW, .~. -Season)
ShanWa1 <- update(ShanWa, .~. -Treatment)
ShanWa2 <- update(ShanWa, .~. -Livestockdensity)

anova(ShanW,ShanWa)
anova(ShanW,ShanWb)
anova(ShanW,ShanWc)
anova(ShanWa,ShanWa1)
anova(ShanWa,ShanWa2)

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq) 
#ShanW  10 100.71 132.64 -40.356   80.713 7.1883      2    0.02748 * # Livestockdensity:Treatment
#ShanW  10 100.71 132.64 -40.356   80.713 6.0508      1     0.0139 * #rainmm
#ShanW  10 100.71 132.64 -40.356   80.713 5.8807      1    0.01531 * #Season
#ShanWa   8 103.90 129.44 -43.951   87.901 1.7777      1     0.1824 # Treatment
#ShanWa   8 103.9 129.44 -43.951   87.901 0.7957      2     0.6718 #Livestockdensity

# Eveness of woody
EveW<-lmer(eveness~Livestockdensity+Treatment+#Season+rainmm+
             #Livestockdensity:Season+#Season:Treatment+
             Livestockdensity:Treatment+  # rainmm:Season+ 
             #rainmm:Livestockdensity+rainmm:Treatment+
             #Treatment:Livestockdensity:Season+
             #Treatment:Livestockdensity:rainmm+
             (1 |Block), data=nsSpp3W)
summary(EveW)
anova(EveW) 
AIC(EveW) # 140.0116
drop1(EveW, test="Chi") # ALL NS
bwplot(eveness~Livestockdensity|Treatment,nsSpp3W)

# Update model
EveWa <- update(EveW, .~. -Livestockdensity:Treatment)
EveWa1 <- update(EveWa, .~. -Livestockdensity)
EveWa2 <- update(EveWa, .~. -Treatment)

anova(EveW,EveWa)
anova(EveWa,EveWa1)
anova(EveWa,EveWa2)

#      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#EveW   8 123.25 147.34 -53.627   107.25 11.466      2   0.003237 **#Livestockdensity:Treatment
#EveWa   6 130.72 148.78 -59.360   118.72 0.6755      2     0.7134#Livestockdensity
#EveWa   6 130.72 148.78 -59.360   118.72 1.235      1     0.2664#Treatment

# Woody species beta part
# Remove rows with sum = zero
names(nsSpp3W[,10:19])

#### Community diversity based on abundance ####
nsSpp3WLow<-droplevels(nsSpp3W[nsSpp3W$Livestockdensity=="Low",])
nsSpp3WMed<-droplevels(nsSpp3W[nsSpp3W$Livestockdensity=="Medium",])
nsSpp3WHigh<-droplevels(nsSpp3W[nsSpp3W$Livestockdensity=="High",])

OverallPresabsWL<-ifelse(nsSpp3WLow[,10:19]>0,1,0)
OverallPresabsWM<-ifelse(nsSpp3WMed[,10:19]>0,1,0)
OverallPresabsWH<-ifelse(nsSpp3WHigh[,10:19]>0,1,0)

beta.sample(OverallPresabsWL, index.family="sorensen", sites=nrow(nsSpp3WLow), samples = 1)
beta.sample(OverallPresabsWM, index.family="sorensen", sites=nrow(nsSpp3WMed), samples = 1)
beta.sample(OverallPresabsWH, index.family="sorensen", sites=nrow(nsSpp3WHigh), samples = 1)
# Low # 0.97075948
# Med # 0.9563878  # Lower diversity in medium - grasses
# High #0.97089562

# Create a factor to define each plot shared across seasons
nsSpp3W$plot_codeII<-as.factor(with(nsSpp3W, paste(Transect,Block,Trtname,Replicate, sep="_")))

# Compare diversity through time
levels(nsSpp3W$Date)
nsSpp3WSeasonI<-droplevels(nsSpp3W[nsSpp3W$Date=="21.10.2012",])
nsSpp3WSeasonII<-droplevels(nsSpp3W[nsSpp3W$Date=="21.11.2013",])
nsSpp3WSeasonIII<-droplevels(nsSpp3W[nsSpp3W$Date=="21.6.2013",])

# Presence and absence of species
names(nsSpp3WSeasonI[,10:19])
names(nsSpp3WSeasonII[,10:19])
names(nsSpp3WSeasonIII[,10:19])
WpresabsI<-ifelse(nsSpp3WSeasonI[,10:19]>0,1,0)
WpresabsII<-ifelse(nsSpp3WSeasonII[,10:19]>0,1,0)
WpresabsIII<-ifelse(nsSpp3WSeasonIII[,10:19]>0,1,0)

# Convert to betapart objects
WpresabsI.core <- betapart.core(WpresabsI)
WpresabsII.core <- betapart.core(WpresabsII)
WpresabsIII.core <- betapart.core(WpresabsIII)

# Assign plot (without season) to each row
row.names(WpresabsI) <- paste(nsSpp3WSeasonI$plot_codeII, 1:nrow(WpresabsI), sep="")
row.names(WpresabsII) <- paste(nsSpp3WSeasonII$plot_codeII, 1:nrow(WpresabsII), sep="")
row.names(WpresabsIII) <- paste(nsSpp3WSeasonIII$plot_codeII, 1:nrow(WpresabsIII), sep="")

dimnames(WpresabsI) 
dimnames(WpresabsII)
dimnames(WpresabsIII)

# Beta.temp
WbtI<-beta.temp(WpresabsI,WpresabsII,index.family="sor") # this is soreson
WbtII<-beta.temp(WpresabsII,WpresabsIII,index.family="sor") 

nsSpp3WS2_3<-rbind(nsSpp3WSeasonII,nsSpp3WSeasonIII)
betaSeasonW<-rbind(WbtI,WbtII)
nsSpp3Wbeta<-cbind(nsSpp3WS2_3,betaSeasonW)

# Analyse turnover component (sor) - Woody species
turnWSor<-lmer(beta.sor~Livestockdensity+Treatment+Season+
                # Livestockdensity:Season+Season:Treatment+
                 Livestockdensity:Treatment+
                 #Treatment:Livestockdensity:Season+
                 (1 |Block), data=nsSpp3Wbeta)
summary(turnWSor)
anova(turnWSor) 
AIC(turnWbt) # 46.30523
drop1(turnWSor, test="Chi")

# Update model woody species
turnWSora <- update(turnWSor, .~. -Livestockdensity:Treatment)
turnWSora2 <- update(turnWSora, .~. -Livestockdensity)
turnWSora3 <- update(turnWSora, .~. -Treatment)

anova(turnWSor,turnWSora)
anova(turnWSora,turnWSora2)
anova(turnWSora,turnWSora3)

# Woody beta diversity
#          Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
#turnWSor   9 121.88 149.56 -51.939   103.88 15.99      2  0.0003371 *** #Livestockdensity:Treatment
#turnWSora   7 133.87 155.39 -59.934   119.87 2.8296      2      0.243#Livestockdensity
#turnWSora   7 133.87 155.39 -59.934   119.87 3.9743      1     0.0462 * # Treatment

# Analyse turnover component (sim) - Woody species
turnWbt<-lmer(beta.sim~Treatment+Season+rainmm+ Livestockdensity+
                Livestockdensity:Season+ #Season:Treatment+ 
                Livestockdensity:Treatment+ #rainmm:Season+
                rainmm:Livestockdensity+#rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3Wbeta)
summary(turnWbt)
anova(turnWbt) 
AIC(turnWbt) # 46.30523
drop1(turnWbt, test="Chi")

# Resid versus fit
Ew1 <- resid(turnWbt,type="pearson") 
Fw1 <- fitted(turnWbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = Fc1,  y = Ec1, xlab = "Fitted values",ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

# Update and remove factors # issues with interactions
turnWbt2<-lmer(beta.sim~Treatment+Season+rainmm+ Livestockdensity+(1 |Block), data=nsSpp3Wbeta)
turnWbta <- update(turnWbt, .~. -rainmm:Livestockdensity) 
turnWbtb <- update(turnWbt, .~. -Livestockdensity:Treatment) 
turnWbtc <- update(turnWbt, .~. -Livestockdensity:Season) 
turnWbt2a <- update(turnWbt2, .~. -rainmm)
turnWbt2b <- update(turnWbt2, .~. -Season)
turnWbt2c <- update(turnWbt2, .~. -Treatment)
turnWbt2d <- update(turnWbt2, .~. -Livestockdensity)

anova(turnWbt,turnWbta)
anova(turnWbt,turnWbtb)
anova(turnWbt,turnWbtc)
anova(turnWbt2,turnWbt2a)
anova(turnWbt2,turnWbt2b)
anova(turnWbt2,turnWbt2c)
anova(turnWbt2,turnWbt2d)

#         Df      AIC    BIC logLik deviance  Chisq Chi Df Pr(>Chisq) 
#turnWbt  14 -12.7871 27.572 20.393  -40.787 8.4276      2    0.01479 * # rainmm:Livestockdensity
#turnWbt  14 -12.7871 27.572 20.393  -40.787 7.2962      2    0.02604 * #Livestockdensity:Treatment
#turnWbt  14 -12.7871 27.572 20.393  -40.787 8.5324      2    0.01404 * # Livestockdensity:Season
#turnWbt2   8 -10.257 12.8052 13.129  -26.257 1.5478      1     0.2135 # rainmm
#turnWbt2   8 -10.257 12.8052 13.129  -26.257 1.5415      1     0.2144 #Season
#turnWbt2   8 -10.257 12.8052 13.129  -26.257  0.81      1     0.3681 #Treatment
#turnWbt2   8 -10.257 12.8052 13.129  -26.257 0.0349      2     0.9827 #Livestockdensity

# Analyse nested component (sne) - Herbs and climbers
nestWbt<-lmer(beta.sne~Livestockdensity+Season+rainmm+#Treatment+
                Livestockdensity:Season+
                rainmm:Livestockdensity+
                #+Season:Treatment+  
                #Livestockdensity:Treatment+#rainmm:Season+
                #rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3Wbeta)
summary(nestWbt)
anova(nestWbt) #
AIC(nestWbt) #-21.98277
drop1(nestWbt, test="Chi") # 

# Resid versus fit
EW1 <- resid(nestWbt,type="pearson") #THIS IS FOR lme..NOT lme4
FW1 <- fitted(nestWbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = FW1, 
     y = EW1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)


nestWbta<-lmer(beta.sne~Livestockdensity+Season+rainmm+#Treatment+
                 Livestockdensity:Season+
                 #rainmm:Livestockdensity+
                 #+Season:Treatment+  
                 #Livestockdensity:Treatment+#rainmm:Season+
                 #rainmm:Treatment+
                 #Treatment:Livestockdensity:Season+
                 #Treatment:Livestockdensity:rainmm+
                 (1 |Block), data=nsSpp3Wbeta)
  
# Update and remove factors # issues with interactions
nestWbtb <- update(nestWbt, .~. -Livestockdensity:Season)
nestWbt2 <- update(nestWbta, .~. -Livestockdensity:Season)
nestWbt2a <- update(nestWbt2, .~. -rainmm)
nestWbt2b <- update(nestWbt2, .~. -Season)
nestWbt2c <- update(nestWbt2, .~. -Livestockdensity)

anova(nestWbt,nestWbta)
anova(nestWbt,nestWbtb)
anova(nestWbt2,nestWbt2a)
anova(nestWbt2,nestWbt2b)
anova(nestWbt2,nestWbt2c)

#         Df     AIC     BIC logLik deviance Chisq Chi Df Pr(>Chisq)
#nestWbt  11 -72.821 -41.110 47.410  -94.821 7.3396      2    0.02548 *#rainmm:Livestockdensity 
#nestWbt  11 -72.821 -41.110 47.410  -94.821 7.6082      2    0.02228 * #Livestockdensity:Season
#nestWbt2   7 -69.225 -49.045 41.612  -83.225 0.1266      1      0.722 #rainmm
#nestWbt2   7 -69.225 -49.045 41.612  -83.225 0.1368      1     0.7115 #Season
#nestWbt2   7 -69.225 -49.045 41.612  -83.225 2.9739      2     0.2261 #Livestockdensity


##### Herbs and climbers  #####

# Shannon index, species richness and Eveness
nsSpp3C$Shannon<-diversity(nsSpp3C[,10:38], index = "shannon")
nsSpp3C$richness<-specnumber(nsSpp3C[,10:38], MARGIN = 1)
nsSpp3C$eveness<-nsSpp3C$Shannon/log(nsSpp3C$richness)

# Drop first date
nsSpp3Ci<-droplevels(nsSpp3C[nsSpp3C$Date!="21.10.2012",])

# Shannon diversity
ShanC<-lmer(Shannon~Livestockdensity+#Season+rainmm+Treatment+
              #Livestockdensity:Season+Season:Treatment+
              #Livestockdensity:Treatment+   #rainmm:Season+ 
              #rainmm:Livestockdensity+rainmm:Treatment+
              #Treatment:Livestockdensity:Season+
              #Treatment:Livestockdensity:rainmm+
              (1 |Block), data=nsSpp3Ci)
summary(ShanC)
anova(ShanC) 
AIC(ShanC) #  177.3046
drop1(ShanC, test="Chi") # ALL NS

# Update model
ShanCa <- update(ShanC, .~. -Livestockdensity)
anova(ShanC,ShanCa)

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#ShanC   5 228.13 244.09 -109.06   218.13 38.332      2  4.746e-09 ***

boxplot(Shannon~Livestockdensity,nsSpp3Ci) # Much lower in high>low>medium

# Evenness
EveC<-lmer(eveness~Livestockdensity+Season+rainmm+Treatment+
              Livestockdensity:Season+Season:Treatment+
              Livestockdensity:Treatment+   #rainmm:Season+ 
              rainmm:Livestockdensity+rainmm:Treatment+
              Treatment:Livestockdensity:Season+
              Treatment:Livestockdensity:rainmm+
              (1 |Block), data=nsSpp3Ci)
summary(EveC)
anova(EveC) 
AIC(EveC) #  177.3046
drop1(EveC, test="Chi") # ALL NS
bwplot(eveness~Livestockdensity|Treatment*Season,nsSpp3Ci)


# Update model
EveCa <- update(EveC, .~. -Livestockdensity:Season:Treatment)
EveCb <- update(EveC, .~. -Livestockdensity:rainmm:Treatment)
EveCc <- update(EveCa, .~. -Livestockdensity:rainmm:Treatment)
EveCc1 <- update(EveCc, .~. -rainmm:Treatment)
EveCc2 <- update(EveCc, .~. -rainmm:Livestockdensity)
EveCc3 <- update(EveCc, .~. -Livestockdensity:Treatment)
EveCc4 <- update(EveCc, .~. -Season:Treatment)
EveCc5 <- update(EveCc, .~. -Livestockdensity:Season)
EveCd<-lmer(eveness~Livestockdensity+Season+rainmm+Treatment+(1 |Block), data=nsSpp3Ci)
EveCd1 <- update(EveCd, .~. -Treatment)
EveCd2 <- update(EveCd, .~. -rainmm)
EveCd3 <- update(EveCd, .~. -Season)
EveCd4 <- update(EveCd, .~. -Livestockdensity)


anova(EveC,EveCa)

#       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
#EveC  20 69.154 127.84 -14.577   29.154 7.0043      2    0.03013 * #Livestockdensity:Season:Treatment

boxplot(Shannon~Livestockdensity,nsSpp3Ci) # Much lower in high>low>medium

# Remove rows with sum = zero
nsSpp3Cb<- nsSpp3C[ rowSums(nsSpp3C[,10:38])!=0, ] 
names(nsSpp3Cb[,10:38])
#### Community diversity based on abundance ####
nsSpp3CLow<-droplevels(nsSpp3C[nsSpp3C$Livestockdensity=="Low",])
nsSpp3CMed<-droplevels(nsSpp3C[nsSpp3C$Livestockdensity=="Medium",])
nsSpp3CHigh<-droplevels(nsSpp3C[nsSpp3C$Livestockdensity=="High",])

OverallPresabsCL<-ifelse(nsSpp3CLow[,10:38]>0,1,0)
OverallPresabsCM<-ifelse(nsSpp3CMed[,10:38]>0,1,0)
OverallPresabsCH<-ifelse(nsSpp3CHigh[,10:38]>0,1,0)

beta.sample(OverallPresabsCL, index.family="sorensen", sites=nrow(nsSpp3CLow), samples = 1)
beta.sample(OverallPresabsCM, index.family="sorensen", sites=nrow(nsSpp3CMed), samples = 1)
beta.sample(OverallPresabsCH, index.family="sorensen", sites=nrow(nsSpp3CHigh), samples = 1)
# Low # 0.97461346 
# Med # 0.97409506  # Lower diversity in medium - grasses
# High #0.97105263

# Create a factor to define each plot shared across seasons
nsSpp3C$plot_codeII<-as.factor(with(nsSpp3C, paste(Transect,Block,Trtname,Replicate, sep="_")))

# Compare diversity through time
levels(nsSpp3Cb$Date)
nsSpp3CSeasonI<-droplevels(nsSpp3C[nsSpp3C$Date=="21.10.2012",])
nsSpp3CSeasonII<-droplevels(nsSpp3C[nsSpp3C$Date=="21.11.2013",])
nsSpp3CSeasonIII<-droplevels(nsSpp3C[nsSpp3C$Date=="21.6.2013",])

# Presence and absence of species
names(nsSpp3CSeasonI[,10:38])
names(nsSpp3CSeasonII[,10:38])
names(nsSpp3CSeasonIII[,10:38])
CpresabsI<-ifelse(nsSpp3CSeasonI[,10:38]>0,1,0)
CpresabsII<-ifelse(nsSpp3CSeasonII[,10:38]>0,1,0)
CpresabsIII<-ifelse(nsSpp3CSeasonIII[,10:38]>0,1,0)

# Convert to betapart objects
CpresabsI.core <- betapart.core(CpresabsI)
CpresabsII.core <- betapart.core(CpresabsII)
CpresabsIII.core <- betapart.core(CpresabsIII)

# Assign plot (without season) to each row
row.names(CpresabsI) <- paste(nsSpp3CSeasonI$plot_codeII, 1:nrow(CpresabsI), sep="")
row.names(CpresabsII) <- paste(nsSpp3CSeasonII$plot_codeII, 1:nrow(CpresabsII), sep="")
row.names(CpresabsIII) <- paste(nsSpp3CSeasonIII$plot_codeII, 1:nrow(CpresabsIII), sep="")

dimnames(CpresabsI) 
dimnames(CpresabsII)
dimnames(CpresabsIII)

# Beta.temp
CbtI<-beta.temp(CpresabsI,CpresabsII,index.family="sor") # this is soreson
CbtII<-beta.temp(GpresabsII,GpresabsIII,index.family="sor") 

nsSpp3CS2_3<-rbind(nsSpp3CSeasonII,nsSpp3CSeasonIII)
betaSeasonC<-rbind(CbtI,CbtII)
nsSpp3Cbeta<-cbind(nsSpp3CS2_3,betaSeasonC)

# Analyse turnover component (sor) - Herbs and climbers
turnCSor<-lmer(beta.sor~Livestockdensity+Treatment+Season+
                  Livestockdensity:Season+Season:Treatment+
                 Livestockdensity:Treatment+
                 Treatment:Livestockdensity:Season+
                 (1 |Block), data=nsSpp3Cbeta)
summary(turnCSor)
anova(turnCSor) 
AIC(turnCSor) # 71.35126
drop1(turnCSor, test="Chi")

# Update model woody species
turnCSora <- update(turnCSor, .~. -Livestockdensity:Treatment:Season)
turnCSora2 <- update(turnCSora, .~. -Livestockdensity:Treatment)
turnCSora3 <- update(turnCSora, .~. -Season:Treatment)
turnCSora4 <- update(turnCSora, .~. - Livestockdensity:Season)

anova(turnCSor,turnCSora)
anova(turnWSora,turnWSora2)
anova(turnWSora,turnWSora3)

# Woody beta diversity
#          Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
#turnWSor   9 121.88 149.56 -51.939   103.88 15.99      2  0.0003371 *** #Livestockdensity:Treatment
#turnWSora   7 133.87 155.39 -59.934   119.87 2.8296      2      0.243#Livestockdensity
#turnWSora   7 133.87 155.39 -59.934   119.87 3.9743      1     0.0462 * # Treatment


# Analyse turnover component (sim) - Herbs and climbers
turnCbt<-lmer(beta.sim~Treatment+Season+#rainmm+ #Livestockdensity
                Season:Treatment+ #Livestockdensity:Season
                #rainmm:Season+ #Livestockdensity:Treatment+
                #rainmm:Livestockdensity+rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3Cbeta)
summary(turnCbt)
anova(turnCbt) 
AIC(turnCbt) # 83.73584
drop1(turnCbt, test="Chi") # ALL NS

# Resid versus fit
Ec1 <- resid(turnCbt,type="pearson") 
Fc1 <- fitted(turnCbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = Fc1,  y = Ec1, xlab = "Fitted values",ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

# Update and remove factors # issues with interactions
turnCbt2 <- update(turnCbt, .~. -Season:Treatment) # No three way interactions..
turnCbt2a <- update(turnCbt2, .~. -Treatment)
turnCbt2b <- update(turnCbt2, .~. -Season)

anova(turnCbt,turnCbt2)
anova(turnCbt2,turnCbt2a)
anova(turnCbt2,turnCbt2b)

#turnCbt: beta.sim ~ Treatment + Season + Season:Treatment + (1 | Block)
#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
#turnCbt   6 48.896 67.496 -18.448   36.896 7.0705      1   0.007836 ** # Season:Treatment
#turnCbt2   5 53.967 69.466 -21.983   43.967 7.9182      1   0.004894 ** #Treatment
#turnCbt2   5 53.967 69.466 -21.984   43.967 30.82      1  2.831e-08 *** # Season


# Analyse nested component (sne) - Herbs and climbers
nestCbt<-lmer(beta.sne~Livestockdensity+Treatment+Season+rainmm+
                +Season:Treatment+ # Livestockdensity:Season
                Livestockdensity:Treatment+#rainmm:Season+
                #rainmm:Livestockdensity+rainmm:Treatment+
                #Treatment:Livestockdensity:Season+
                #Treatment:Livestockdensity:rainmm+
                (1 |Block), data=nsSpp3Cbeta)
summary(nestCbt)
anova(nestCbt) #
AIC(nestCbt) #-110.4011
drop1(nestCbt, test="Chi") # 

# Resid versus fit
EC1 <- resid(nestCbt,type="pearson") #THIS IS FOR lme..NOT lme4
FC1 <- fitted(nestCbt)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = FC1, 
     y = EC1,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(v = 100, lwd = 2, col = 2) 
abline(h = 0, lty = 2, col = 1)

Livestockdensity+Treatment+Season+rainmm+
  +Season:Treatment+ # Livestockdensity:Season
  Livestockdensity:Treatment

# Update and remove factors # issues with interactions
nestCbta <- update(nestCbt, .~. -Livestockdensity:Treatment)
nestCbtb <- update(nestCbt, .~. -Season:Treatment)
nestCbt2 <- update(nestCbta, .~. -Season:Treatment)
nestCbt2a <- update(nestCbt2, .~. -rainmm)
nestCbt2b <- update(nestCbt2, .~. -Season)
nestCbt2c <- update(nestCbt2, .~. -Treatment)
nestCbt2d <- update(nestCbt2, .~. -Livestockdensity)

anova(nestCbt,nestCbta)
anova(nestCbt,nestCbtb)
anova(nestCbt2,nestCbt2a)
anova(nestCbt2,nestCbt2b)
anova(nestCbt2,nestCbt2c)
anova(nestCbt2,nestCbt2d)

#         Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq) 
#nestCbt  11 -162.98 -128.88 92.489  -184.98 18.835      2   8.13e-05 ***#Livestockdensity:Treatment
#nestCbt  11 -162.98 -128.88 92.489  -184.98 8.462      1   0.003626 ** #Season:Treatment
#nestCbt2   8 -143.72 -118.92 79.861  -159.72 7.804      1   0.005213 ** #rainmm
#nestCbt2   8 -143.72 -118.92 79.861  -159.72 7.2161      1   0.007225 ** #Season
#nestCbt2   8 -143.72 -118.92 79.861  -159.72 6.2618      1    0.01234 * #Treatment
#nestCbt2   8 -143.72 -118.92 79.861  -159.72 2.1265      2     0.3453 # Livestockdensity

##########################################################################################
#### END #### 
#########################################################################################