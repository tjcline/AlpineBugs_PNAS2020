rm(list=ls())
#quartz()
#Code for Alpine insect diversity, Glacier National Park
#Timothy Cline
#8/28/2019
library(dplyr)
library(lme4)
library(MuMIn)
library(merTools)
library(doParallel)
library(RColorBrewer)
library(vegan)

#Required functions for dealing with proportions
logit<-function(x){return(log(x/(1-x)))}
iLogit<-function(x){return(exp(x)/(1+exp(x)))}

setwd('~/Documents/AlpineBugs')

#load bug sampling and glacier/snow datasets
AlpineMaster<-read.csv('AlpineMasterData20191122_FewerGroups_GlacialSprings.csv',header=T,stringsAsFactors=FALSE)
IceAndSnow<-read.csv('LIA calcs_20191122_2.csv',header=T,stringsAsFactors=F)

#Create table of abundance by site and species
SpeciesAbundance<-AlpineMaster[,which(colnames(AlpineMaster)=="Liodessus_affinis"):ncol(AlpineMaster)]
SpeciesAbundance<-SpeciesAbundance[,colSums(SpeciesAbundance)>0]

#Create table of presence/absence by site and species
SpeciesPresence <- data.frame(matrix(as.numeric(SpeciesAbundance>0),nrow=nrow(SpeciesAbundance),ncol=ncol(SpeciesAbundance)))
colnames(SpeciesPresence) <- colnames(SpeciesAbundance)

#Calculate Richness and Diversity
AlpineMaster <- AlpineMaster %>% mutate(Richness = rowSums(SpeciesPresence), Diversity = apply(SpeciesAbundance,1,FUN=diversity))

#compute ice statistics
AlpineMaster$Pcnt_Just_Ice <- 100 * (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])
AlpineMaster$LIAprop <- IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]
AlpineMaster$PropLoss <- 1-((IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])/(IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]))#LIA$JG.prop.loss[LIAsiteInd]
AlpineMaster$PropLoss[is.na(AlpineMaster$PropLoss)]<-0 #zeros for sites that had no glaciers at LIA
AlpineMaster$PropDiff <- (IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]) - (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])


###Choose glacial index for analysis

#Quick reference for covariates
GCC <- AlpineMaster$Pcnt_Just_Ice
TEMP <- AlpineMaster$ts8_base
DIST <- AlpineMaster$Distance
ELEV <- AlpineMaster$Elevation

RICH <- AlpineMaster$Richness
DIV <- AlpineMaster$Diversity
ABUND <- rowSums(SpeciesAbundance)
LABUND <- log(ABUND)

#Covariate Transformations
LDIST<-log(DIST)
LGCC <- logit(GCC/100+0.01)
LRICH<-log(RICH)
LDIV<-log(DIV+0.1)

#Grouping variables
DRAIN<-AlpineMaster$DrainageCode
STREAM<-AlpineMaster$StreamNameCode

#### Calculating Beta Diversity

DRC<-matrix(unique(AlpineMaster$DrainageCode),ncol=1)
DrainageComposition_Abundance<-t(apply(DRC,1,FUN=function(x){colMeans(SpeciesAbundance[which(AlpineMaster$DrainageCode==x),])}))
DrainageComposition_Presence<-t(apply(DRC,1,FUN=function(x){colSums(SpeciesAbundance[which(AlpineMaster$DrainageCode==x),])>0}))
row.names(DrainageComposition_Abundance)<-DRC
row.names(DrainageComposition_Presence)<-DRC

RegionalDiversity_Abundance <- colMeans(DrainageComposition_Abundance)
RegionalDiversity_Presence <- (colSums(DrainageComposition_Presence)>0)

#calculate turnover (b-diversity) between sites and the regional species pool on abundance and presence/absence
for(i in 1:nrow(AlpineMaster)){
   BD_regional_abundance<-vegdist(rbind(SpeciesAbundance[i,],RegionalDiversity_Abundance),method='jaccard')
   AlpineMaster$BetaDiversityAbundance_Region[i]<-BD_regional_abundance
  
    BD_regional_presence<-vegdist(rbind(SpeciesPresence[i,],RegionalDiversity_Presence),method='jaccard')
    AlpineMaster$BetaDiversityPresence_Region[i]<-BD_regional_presence
}


#TestRandomEffect
lm.rich.full.RE1<-lmer(LRICH ~ LGCC*ELEV*TEMP*LDIST + (1|STREAM),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.rich.full.RE2<-lmer(LRICH ~ LGCC*ELEV*TEMP*LDIST + (1|DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.rich.full.fullRE<-lmer(LRICH ~ LGCC*ELEV*TEMP*LDIST + (1|STREAM/DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
AIC(lm.rich.full.RE1);AIC(lm.rich.full.RE2);AIC(lm.rich.full.fullRE)
#Drainage is the best grouping variable



### ALL REGRESSION ANALYSES FOR REGRESSIONS INCLUDED IN FIG. 2 (MUHLFELD ET AL. 2020, PNAS)
cl <- makeCluster(detectCores())
registerDoParallel(cl)

#Glacier prediction series for interval prediction
PredGlac<-seq(min(LGCC),max(LGCC),length=129)

RichRegBase<-lmer(LRICH ~ 1 + (1|DRAIN),REML=F)
RichReg<-lmer(LRICH~LGCC + (1|DRAIN),REML=F)
dAIC_RichReg<-round(AICc(RichReg)-AICc(RichRegBase),1)
RichRegPred<-predictInterval(RichReg,newdata=data.frame(LGCC=PredGlac),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

DivRegBase<-lmer(LDIV ~ 1 + (1|DRAIN),REML=F)
DivReg<-lmer(LDIV~LGCC + (1|DRAIN),REML=F)
DivRegPred<-predictInterval(DivReg,newdata=data.frame(LGCC=PredGlac),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)
dAIC_DivReg<-round(AICc(DivReg)-AICc(DivRegBase),1)

LBETA<-logit(AlpineMaster$BetaDiversityPresence_Region)
BetaRegBase<-lmer(LBETA ~ 1 + (1|DRAIN),REML=F)
BetaReg<-lmer(LBETA ~ LGCC + (1|DRAIN),REML=F)
BetaRegPred<-predictInterval(BetaReg,newdata=data.frame(LGCC=PredGlac),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)
dAIC_BetaReg<-round(AICc(BetaReg)-AICc(BetaRegBase),1)

#Glacier prediction series for interval prediction
PredTemp<-seq(min(TEMP),max(TEMP),length=129)
RichReg_Temp<-lmer(LRICH~TEMP + (1|DRAIN),REML=F)
RichRegPred_Temp<-predictInterval(RichReg_Temp,newdata=data.frame(TEMP=PredTemp),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

DivReg_Temp<-lmer(LDIV~TEMP + (1|DRAIN),REML=F)
DivRegPred_Temp<-predictInterval(DivReg_Temp,newdata=data.frame(TEMP=PredTemp),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

BetaReg_Temp<-lmer(LBETA~TEMP + (1|DRAIN),REML=F)
BetaRegPred_Temp<-predictInterval(BetaReg_Temp,newdata=data.frame(TEMP=PredTemp),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

stopCluster(cl)

save(list=ls(),file='AlpineBugsAnalysis_Short.Rdata')

#quick model selection tables
Richtab <- data.frame(logLik= c(round(logLik(RichReg),2),round(logLik(RichReg_Temp),2),round(logLik(RichRegBase),2)),
                       NumPars = c(3,3,2),
                       df = c(125,125,126),
                       AICc=c(round(AICc(RichReg),2),round(AICc(RichReg_Temp),2),round(AICc(RichRegBase),2)),
                       dAICc=c(round(AICc(RichReg)-AICc(RichRegBase),2),round(AICc(RichReg_Temp)-AICc(RichRegBase),2),'--'))
#quick model selection tables
Divtab <- data.frame(logLik= c(round(logLik(DivReg),2),round(logLik(DivReg_Temp),2),round(logLik(DivRegBase),2)),
                       NumPars = c(3,3,2),
                       df = c(125,125,126),
                       AICc=c(round(AICc(DivReg),2),round(AICc(DivReg_Temp),2),round(AICc(DivRegBase),2)),
                       dAICc=c(round(AICc(DivReg)-AICc(DivRegBase),2),round(AICc(DivReg_Temp)-AICc(DivRegBase),2),'--'))
#quick model selection tables
Betatab <- data.frame(logLik= c(round(logLik(BetaReg),2),round(logLik(BetaReg_Temp),2),round(logLik(BetaRegBase),2)),
                       NumPars = c(3,3,2),
                       df = c(125,125,126),
                       AICc=c(round(AICc(BetaReg),2),round(AICc(BetaReg_Temp),2),round(AICc(BetaRegBase),2)),
                       dAICc=c(round(AICc(BetaReg)-AICc(BetaRegBase),2),round(AICc(BetaReg_Temp)-AICc(BetaRegBase),2),'--'))
                       
write.table(Richtab, file='~/Desktop/RichTab.csv',sep=',')
write.table(Divtab, file='~/Desktop/DivTab.csv',sep=',')
write.table(Betatab, file='~/Desktop/BetaTab.csv',sep=',')


#### Analysis asking which taxa are found in only glacial streams
unidrain<-unique(AlpineMaster$Drainage)
SpeciesPresence_Drainage<-matrix(NA,nrow=length(unidrain),ncol=ncol(SpeciesPresence))
for(i in 1:length(unidrain)){
  SpeciesPresence_Drainage[i,]<-as.numeric(colSums(SpeciesPresence[which(AlpineMaster$Drainage==unidrain[i]),])>0)
}

#Which Taxa are only found in 1 drainage
colnames(SpeciesPresence)[which(colSums(SpeciesPresence_Drainage)==1)]
#Which Taxa are only found in 2 drainages
colnames(SpeciesPresence)[which(colSums(SpeciesPresence_Drainage)<=2)]

#Only taxa found in more than 2 sample sites
SpeciesPresence_Common<-SpeciesPresence[,which(colSums(SpeciesPresence)>2)]

#Species in Glacial but NOT non-glacial Streams
GNG_IND<-which(colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Just_Ice>0),])>0 & colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Just_Ice==0),])==0)
colnames(SpeciesPresence_Common)[GNG_IND]

