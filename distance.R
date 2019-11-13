##_________##
#### TOP ####
##_________##
rm(list=ls())

#library(emmeans)
#library(lme4)
#library(ggplot2)
#library(ggpubr)
#install.packages("dplyr")
library(dplyr)
library(tidyr)

setwd("C:\\Users\\Kjirsten\\Documents\\R")
dist <- read.csv("distance.csv")
str(dist)

names(dist)

dist1 <- subset(dist, NEAR_FID == 0)
dist1$exclosure <- "1Ha"

dist2 <- subset(dist, NEAR_FID == 1)
dist2$exclosure <- "1Hb"

dist3 <- subset(dist, NEAR_FID == 2)
dist3$exclosure <- "1La"

dist4 <- subset(dist, NEAR_FID == 3)
dist4$exclosure <- "1Lb"

dist5 <- subset(dist, NEAR_FID == 4)
dist5$exclosure <- "1Ma"

dist6 <- subset(dist, NEAR_FID == 5)
dist6$exclosure <- "1Mb"

dist7 <- subset(dist, NEAR_FID == 6)
dist7$exclosure <- "2Ha"

dist8 <- subset(dist, NEAR_FID == 7)
dist8$exclosure <- "2Hb"

dist9 <- subset(dist, NEAR_FID == 8)
dist9$exclosure <- "2La"

dist9 <- subset(dist, NEAR_FID == 8)
dist9$exclosure <- "2La"

dist10 <- subset(dist, NEAR_FID == 9)
dist10$exclosure <- "2Lb"

dist11 <- subset(dist, NEAR_FID == 10)
dist11$exclosure <- "2Ma"

dist12 <- subset(dist, NEAR_FID == 11)
dist12$exclosure <- "2Mb"

dist13 <- subset(dist, NEAR_FID == 12)
dist13$exclosure <- "3Ha"

dist14 <- subset(dist, NEAR_FID == 13)
dist14$exclosure <- "3Hb"

dist15 <- subset(dist, NEAR_FID == 14)
dist15$exclosure <- "3La"

dist16 <- subset(dist, NEAR_FID == 15)
dist16$exclosure <- "3Lb"

dist17 <- subset(dist, NEAR_FID == 16)
dist17$exclosure <- "3Ma"

dist18 <- subset(dist, NEAR_FID == 17)
dist18$exclosure <- "3Mb"

distance <- rbind(dist1, dist2, dist3, dist4, dist5, dist6, dist7, dist8,
                  dist9, dist10, dist11, dist12, dist13, dist14, dist15,
                  dist16, dist17, dist18)
names(distance)

colnames(distance)[colnames(distance)=="IN_FID"] <- "herb.poly.id"
