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

