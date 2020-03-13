rm(list=ls())
#quartz()
#Code for Alpine insect diversity, Glacier National Park
#Timothy Cline
#8/28/2019
library(dplyr)
library(tidyr)
library(lme4)
library(MuMIn)
library(merTools)
library(doParallel)
library(RColorBrewer)
library(vegan)

logit<-function(x){return(log(x/(1-x)))}
iLogit<-function(x){return(exp(x)/(1+exp(x)))}

setwd('~/Documents/AlpineBugs')
AlpineMaster<-read.csv('AlpineMasterData20191122_FewerGroups_GlacialSprings.csv',header=T,stringsAsFactors=FALSE)
IceAndSnow<-read.csv('LIA calcs_20191122_2.csv',header=T,stringsAsFactors=F)


SpeciesAbundance<-AlpineMaster[,which(colnames(AlpineMaster)=="Liodessus_affinis"):ncol(AlpineMaster)]
SpeciesAbundance<-SpeciesAbundance[,colSums(SpeciesAbundance)>0]

SpeciesPresence <- data.frame(matrix(as.numeric(SpeciesAbundance>0),nrow=nrow(SpeciesAbundance),ncol=ncol(SpeciesAbundance)))
colnames(SpeciesPresence) <- colnames(SpeciesAbundance)

#Calculate Richness and Diversity
AlpineMaster <- AlpineMaster %>% mutate(Richness = rowSums(SpeciesPresence), Diversity = apply(SpeciesAbundance,1,FUN=diversity))

AlpineMaster$DOY <- as.numeric(format(as.Date(AlpineMaster$Date),'%j'))

AlpineMaster$Pcnt_Just_Ice <- 100 * (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])
AlpineMaster$LIAprop <- IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]
AlpineMaster$PropLoss <- 1-((IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])/(IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]))#LIA$JG.prop.loss[LIAsiteInd]
AlpineMaster$PropLoss[is.na(AlpineMaster$PropLoss)]<-0 #zeros for sites that had no glaciers at LIA

AlpineMaster$PropDiff <- (IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]) - (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])

library(betapart)


betapair1<-beta.pair(SpeciesPresence,index.family = 'jaccard')
hist(betapair1$beta.jtu)
hist(betapair1$beta.jne)
hist(betapair1$beta.jac)

quartz(height=4,width=3)
par(mfrow=c(2,1),mar=c(3,3,1,1))
hist(betapair1$beta.jtu,xlim=c(0.1,1),main='',xlab='Species turnover among sites')
hist(betapair1$beta.jne/betapair1$beta.jac,xlim=c(0.1,1),main='',xlab='Species turnover * total diversity-1')

bmp<-beta.multi(x=SpeciesPresence,index.family = 'jaccard')
bmp$beta.JTU/bmp$beta.JAC
bma<-beta.multi.abund(x=SpeciesAbundance,index.family='bray')


