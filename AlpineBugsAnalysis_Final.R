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


###Choose glacial index for analysis
#GLACIER <- AlpineMaster$Glac_indx
#GLACIER <- AlpineMaster$Glac_indx_SnowIce

#GLACIER1 <- AlpineMaster$Glac_indx

#Quick reference for covariates
GCC <- AlpineMaster$Pcnt_Just_Ice#AlpineMaster$Pcnt_Total_Ice_Snow
TEMP <- AlpineMaster$ts8_base
DIST <- AlpineMaster$Distance
ELEV <- AlpineMaster$Elevation
ASPECT <- AlpineMaster$Aspect
DOY <- AlpineMaster$DOY/365
LAKE <- AlpineMaster$Pcnt_Lake





RICH<-AlpineMaster$Richness
DIV<-AlpineMaster$Diversity
ABUND<-rowSums(SpeciesAbundance)
LABUND<-log(ABUND)

#EVEN<-AlpineMaster$Evenness_e.H.S_TAXON
#DOMN<-AlpineMaster$Dominance_D_TAXON

#Grouping
DRAIN<-AlpineMaster$DrainageCode
STREAM<-AlpineMaster$StreamNameCode

#Covariate Transformations
LDIST<-log(DIST)
TASP<-sin(((360-ASPECT) +90)*pi/180)
LGCC <- logit(GCC/100+0.01)
LRICH<-log(RICH)
LDIV<-log(DIV+0.1)
LLAKE <- logit(LAKE/100+0.01)
LB<- as.integer(ifelse(LAKE>0,1,0))

#Looking at correlations between covariates
plot(ELEV,LDIST)
cor(ELEV,LDIST)
hist(LDIST)
hist(DIST)

plot(LDIST,TEMP)
cor(LDIST,TEMP)

plot(ELEV,TEMP)
cor(ELEV,TEMP)

plot(LGCC,TEMP)
cor(LGCC,TEMP)
summary(lm(LGCC~TEMP))

plot(TASP,TEMP)
cor(TASP,TEMP)

#QuickPCA to look at correlations
# PCvars<-prcomp(AlpineMaster %>% select(Elevation,Aspect,Distance,Watershed.area,Pcnt_Lake,Glac_indx,ts8_base),center=T,scale=T)
# PCvars
# summary(PCvars)
# autoplot(PCvars,loadings=T,loadings.label=T)



#TestRandomEffect
lm.rich.full.RE1<-lmer(LRICH ~ LGCC*ELEV*TEMP*LDIST + (1|STREAM),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.rich.full.RE2<-lmer(LRICH ~ LGCC*ELEV*TEMP*LDIST + (1|DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.rich.full.fullRE<-lmer(LRICH ~ LGCC*ELEV*TEMP*LDIST + (1|STREAM/DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))

AIC(lm.rich.full.RE1);AIC(lm.rich.full.RE2);AIC(lm.rich.full.fullRE)


lm.div.full.RE1<-lmer(LDIV ~ LGCC*ELEV*TEMP*LDIST + (1|STREAM),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.div.full.RE2<-lmer(LDIV ~ LGCC*ELEV*TEMP*LDIST + (1|DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.div.full.fullRE<-lmer(LDIV ~ LGCC*ELEV*TEMP*LDIST + (1|STREAM/DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
AIC(lm.div.full.RE1);AIC(lm.div.full.RE2);AIC(lm.div.full.fullRE)
#BestRandomEffectStructure is just drainage


#Drivers of Richness
lm.rich.full<-lmer(LRICH ~ LGCC*ELEV*LDIST*TEMP + (1|DRAIN),na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
rich.dredge<-dredge(lm.rich.full,extra=c('R2c'=function(x){r.squaredGLMM(x, null = nullmodel)[1,2]},'R2m'=function(x){r.squaredGLMM(x, null = nullmodel)[1,1]}))
rich.dredge

lm.div.full<-lmer(LDIV ~ LGCC*TEMP*LDIST*TEMP + (1|DRAIN),na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
div.dredge<-dredge(lm.div.full,extra=c('R2c'=function(x){r.squaredGLMM(x, null = nullmodel)[1,2]},'R2m'=function(x){r.squaredGLMM(x, null = nullmodel)[1,1]}))
div.dredge


#####Abundance
GRAND_TOTAL <- rowSums(SpeciesAbundance)
LGT<-log(GRAND_TOTAL)
lm.ta.full.RE1<-lmer(LGT ~ LGCC*TEMP*LDIST + (1|STREAM),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.ta.full.RE2<-lmer(LGT ~ LGCC*TEMP*LDIST + (1|DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
lm.ta.full.fullRE<-lmer(LGT ~ LGCC*TEMP*LDIST + (1|STREAM/DRAIN),na.action=na.fail,REML=T,control=lmerControl(optimizer="Nelder_Mead"))
AIC(lm.ta.full.RE1);AIC(lm.ta.full.RE2);AIC(lm.ta.full.fullRE)

# GlacBool<-as.numeric(AlpineAbund$Glac_indx>0)
# lm.total.abun.full<-lmer(log(Grand_Total) ~ GlacBool + GlacBool:Glac_indx + (1|DrainageCode),data=AlpineAbund,na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
# summary(lm.total.abun.full)


lm.total.abun.full<-lmer(LGT ~ LGCC*TEMP*LDIST*ELEV + (1|DRAIN),na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
ta.dredge<-dredge(lm.total.abun.full)
ta.dredge

# lm.EPT.abun.full<-lmer(log(EPT_total+1) ~ Glac_indx*ts8_base*logDist + (1|DrainageCode),data=AlpineAbund,na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
# EPTa.dredge<-dredge(lm.EPT.abun.full)
# EPTa.dredge

# lm.EPT.perc.full<-lmer(logit(EPT_pcnt+0.01) ~ Glac_indx*ts8_base*logDist + (1|DrainageCode),data=AlpineAbund,na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
# EPTp.dredge<-dredge(lm.EPT.perc.full)
# EPTp.dredge

#####AccumulatedSpeciesLoss

#sorted<-AlpineMaster %>% arrange(desc(Pcnt_Total_Ice_Snow),desc(Distance))#arrange(desc(Glac_indx),desc(Distance))
#sorted<-AlpineMaster %>% arrange(desc(ts8_base),desc(Distance))
#sorted<-SpeciesPresence[order(AlpineMaster$Elevation,decreasing=T),]
unidrain<-unique(AlpineMaster$Drainage)
SpeciesPresence_Drainage<-matrix(NA,nrow=length(unidrain),ncol=ncol(SpeciesPresence))
for(i in 1:length(unidrain)){
  SpeciesPresence_Drainage[i,]<-as.numeric(colSums(SpeciesPresence[which(AlpineMaster$Drainage==unidrain[i]),])>0)
}

#Which Taxa are only found in 1 drainage
colnames(SpeciesPresence)[which(colSums(SpeciesPresence_Drainage)==1)]
colnames(SpeciesPresence)[which(colSums(SpeciesPresence_Drainage)<=2)]

#Only taxa found in more than 2 drainages
#SpeciesPresence_Common<-SpeciesPresence[,which(colSums(SpeciesPresence_Drainage)>2)]

#Only taxa found in more than 2 sample sites
SpeciesPresence_Common<-SpeciesPresence[,which(colSums(SpeciesPresence)>2)]

#Species in Glacial but NOT non-glacial Streams
GNG_IND<-which(colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Just_Ice>0),])>0 & colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Just_Ice==0),])==0)
colnames(SpeciesPresence_Common)[GNG_IND]


#AlpineMaster[which(SpeciesPresence[,'Stygobromus_glacialis']>0),]

AlpineMaster$Pcnt_Bin<-NA
Breaks<-c(seq(90,10,by=-10),5,1,0) #Categories of GCC

#Assign sites to categories of Ice Cover
AlpineMaster$Pcnt_Bin[AlpineMaster$Pcnt_Just_Ice>Breaks[1]]<-1
for(i in 2:length(Breaks)){
  AlpineMaster$Pcnt_Bin[AlpineMaster$Pcnt_Just_Ice>=Breaks[i] & AlpineMaster$Pcnt_Just_Ice<Breaks[i-1]]<-i
}

#Tally up species occurences by categories of Ice Cover
MAT<-matrix(0,nrow=length(Breaks),ncol=ncol(SpeciesPresence_Common))
for(i in 1:length(Breaks)){
  if(length(which(AlpineMaster$Pcnt_Bin==i))>0){
    MAT[i,]<-as.numeric(colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Bin==i),])>0)
  }
}
rowSums(MAT) #Number of species occuring in each bin

colnames(SpeciesPresence_Common)[MAT[1,]==1] #Species only found at 90% or above Ice
colnames(SpeciesPresence_Common)[MAT[2,]==1] #Species only found at 80-90%
colnames(SpeciesPresence_Common)[MAT[3,]==1] #Species only found at 70-80%

Counts<-rep(NA,nrow(MAT))
Names<-list()
for(i in 1:(nrow(MAT)-1)){
  thisCount<-0
  thisSpecies<-c()
  Present<-MAT[i,]
  
  if(sum(Present)>0){
    
    if(i==1){Above<-rep(0,ncol(MAT))
    }else if(i==2){
      Above<-as.numeric(MAT[1,]>0)
    }else{Above<-as.numeric(colSums(MAT[1:(i-1),])>0)}
    
    Below<-as.numeric(colSums(rbind(rep(0,ncol(MAT)),MAT[(i+1):nrow(MAT),]))>0)
    
    for(j in 1:ncol(MAT)){
      if((Present[j]==1 | Above[j] == 1) & Below[j] == 0){
        thisCount<-thisCount+1
        thisSpecies<-c(thisSpecies,colnames(SpeciesPresence_Common)[j])
      }
    }
  }else{thisCount<-NA}
  Counts[i]<-thisCount
  Names[[i]]<-thisSpecies
}

OldNames<-''
for(i in 1:length(Names)){
  if(is.null(Names[[i]])){
    Names[[i]]<-''
  }else{
    Names[[i]] <- c(Names[[i]][!(Names[[i]] %in% OldNames)],'')
    OldNames<-c(OldNames,Names[[i]])
  }
}

Counts<-unlist(lapply(Names,FUN=length))-1

SpecLoss<-na.omit(cbind(Breaks[1:length(Counts)],Counts))

quartz()
plot(SpecLoss[,1],cumsum(SpecLoss[,2]),type='o',xlab='GCC (%)',ylab='SpeciesAbove_NotBelow',pch=16)
#points(lowess(SpecLoss[,1],SpecLoss[,2],f=2/3),type='l')
# text(61,1,Names[[5]])
# text(47,1.3,Names[[6]])
# for(i in 1:3){
#    text(28+i*2,3-(0.2*i),Names[[7]][i],adj=c(0,1))
# }
# for(i in 1:3){
#    text(18+i*2,3-(0.2*i),Names[[8]][i],adj=c(0,1))
# }
# for(i in 1:6){
#   text(2+2*i,6.3-(0.2*i),Names[[10]][i],adj=c(0,1))
# }

GstreamsSpecList<-colnames(SpeciesPresence_Common)[colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Just_Ice>0),])>0]
NGstreamsSpecList<-colnames(SpeciesPresence_Common)[colSums(SpeciesPresence_Common[which(AlpineMaster$Pcnt_Just_Ice==0),])>0]

GstreamsSpecList[which(!(GstreamsSpecList %in% NGstreamsSpecList))]

sum(is.na(match(GstreamsSpecList,NGstreamsSpecList)))
sum(!is.na(match(GstreamsSpecList,NGstreamsSpecList)))
GstreamsSpecList[is.na(match(GstreamsSpecList,NGstreamsSpecList))]



##### Calculating Beta Diversity

DRC<-matrix(unique(AlpineMaster$DrainageCode),ncol=1)
DrainageComposition_Abundance<-t(apply(DRC,1,FUN=function(x){colMeans(SpeciesAbundance[which(AlpineMaster$DrainageCode==x),])}))
DrainageComposition_Presence<-t(apply(DRC,1,FUN=function(x){colSums(SpeciesAbundance[which(AlpineMaster$DrainageCode==x),])>0}))
row.names(DrainageComposition_Abundance)<-DRC
row.names(DrainageComposition_Presence)<-DRC

RegionalDiversity_Abundance <- colMeans(DrainageComposition_Abundance)
RegionalDiversity_Presence <- (colSums(DrainageComposition_Presence)>0)



for(i in 1:nrow(AlpineMaster)){
  BD_basin_abundance<-vegdist(rbind(SpeciesAbundance[i,],DrainageComposition_Abundance[which(DRC==AlpineMaster$DrainageCode[i]),]),method='jaccard')
  BD_regional_abundance<-vegdist(rbind(SpeciesAbundance[i,],RegionalDiversity_Abundance),method='jaccard')
  AlpineMaster$BetaDiversityAbundance_Basin[i]<-BD_basin_abundance
  AlpineMaster$BetaDiversityAbundance_Region[i]<-BD_regional_abundance
  
  BD_basin_presence<-vegdist(rbind(SpeciesPresence[i,],DrainageComposition_Presence[which(DRC==AlpineMaster$DrainageCode[i]),]),method='jaccard')
  BD_regional_presence<-vegdist(rbind(SpeciesPresence[i,],RegionalDiversity_Presence),method='jaccard')
  AlpineMaster$BetaDiversityPresence_Basin[i]<-BD_basin_presence
  AlpineMaster$BetaDiversityPresence_Region[i]<-BD_regional_presence
}

#BETA_BASIN_ABUNDANCE
# plot(AlpineMaster$BetaDiversityAbundance_Basin~AlpineMaster$Pcnt_Just_Ice)
# points(lowess(AlpineMaster$BetaDiversityAbundance_Basin~AlpineMaster$Pcnt_Just_Ice),type='l')
# 
# plot(AlpineMaster$BetaDiversityAbundance_Basin~AlpineMaster$Distance)
# points(lowess(AlpineMaster$BetaDiversityAbundance_Basin~AlpineMaster$Distance),type='l')
# 
# plot(AlpineMaster$BetaDiversityAbundance_Basin~AlpineMaster$Elevation)
# points(lowess(AlpineMaster$BetaDiversityAbundance_Basin~AlpineMaster$Elevation),type='l')

#BETA_REGION_ABUNDANCE
plot(AlpineMaster$BetaDiversityAbundance_Region~AlpineMaster$Pcnt_Just_Ice)
points(lowess(AlpineMaster$BetaDiversityAbundance_Region~AlpineMaster$Pcnt_Just_Ice),type='l')

plot(AlpineMaster$BetaDiversityAbundance_Region~AlpineMaster$Distance)
points(lowess(AlpineMaster$BetaDiversityAbundance_Region~AlpineMaster$Distance),type='l')

plot(AlpineMaster$BetaDiversityAbundance_Region~AlpineMaster$Elevation)
points(lowess(AlpineMaster$BetaDiversityAbundance_Region~AlpineMaster$Elevation),type='l')


#BETA_BASIN_PRESENCE
plot(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Pcnt_Just_Ice,col=AlpineMaster$GradColor,pch=16,cex=2)
points(lowess(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Pcnt_Just_Ice),type='l')

plot(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Glac_indx)
points(lowess(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Glac_indx),type='l')

plot(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Distance)
points(lowess(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Distance),type='l')

plot(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Elevation)
points(lowess(AlpineMaster$BetaDiversityPresence_Basin~AlpineMaster$Elevation),type='l')

#BETA_REGION_PRESENCE
HistGrad<-hist(AlpineMaster$Distance,plot=F,breaks=c(0,50,100,200,400,800,1600,3200,6400))
BluePalGen<-colorRampPalette(brewer.pal(9,'Blues'))
BluePal<-BluePalGen(length(HistGrad$breaks))
for(i in 1:length(HistGrad$breaks)){
  AlpineMaster$GradColor[AlpineMaster$Distance>HistGrad$breaks[i]]<-rev(BluePal)[i]
}

# plot(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Pcnt_Just_Ice,col=AlpineMaster$GradColor,pch=16,cex=2)
# points(lowess(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Pcnt_Just_Ice),type='l')
# 
# plot(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Glac_indx)
# points(lowess(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Glac_indx),type='l')
# 
# plot(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Distance)
# points(lowess(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Distance),type='l')
# 
# plot(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Elevation)
# points(lowess(AlpineMaster$BetaDiversityPresence_Region~AlpineMaster$Elevation),type='l')
# summary(lm(logit(AlpineMaster$BetaDiversityPresence_Region)~AlpineMaster$Elevation))

#lm.full<-lmer(logit(AlpineMaster$BetaDiversityPresence_Region)~ELEV*TEMP*LGLACIER*DIST + (1|DRAIN),na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
#d1<-dredge(lm.full)




AlpineMaster$Glac_Y_N <- AlpineMaster$Pcnt_Just_Ice>0
GYN<-AlpineMaster$Glac_Y_N

cl <- makeCluster(detectCores())
registerDoParallel(cl)

#LPREDGLACIER<-logit((100-GLACIER)/100-0.01)
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


PredTemp<-seq(min(TEMP),max(TEMP),length=129)
RichReg_Temp<-lmer(LRICH~TEMP + (1|DRAIN),REML=F)
RichRegPred_Temp<-predictInterval(RichReg_Temp,newdata=data.frame(TEMP=PredTemp),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

DivReg_Temp<-lmer(LDIV~TEMP + (1|DRAIN),REML=F)
DivRegPred_Temp<-predictInterval(DivReg_Temp,newdata=data.frame(TEMP=PredTemp),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

BetaReg_Temp<-lmer(LBETA~TEMP + (1|DRAIN),REML=F)
BetaRegPred_Temp<-predictInterval(BetaReg_Temp,newdata=data.frame(TEMP=PredTemp),n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)

stopCluster(cl)

save(list=ls(),file='AlpineBugsAnalysis_Short.Rdata')

Richtab <- data.frame(logLik= c(round(logLik(RichReg),2),round(logLik(RichReg_Temp),2),round(logLik(RichRegBase),2)),
                       NumPars = c(3,3,2),
                       df = c(125,125,126),
                       AICc=c(round(AICc(RichReg),2),round(AICc(RichReg_Temp),2),round(AICc(RichRegBase),2)),
                       dAICc=c(round(AICc(RichReg)-AICc(RichRegBase),2),round(AICc(RichReg_Temp)-AICc(RichRegBase),2),'--'))

Divtab <- data.frame(logLik= c(round(logLik(DivReg),2),round(logLik(DivReg_Temp),2),round(logLik(DivRegBase),2)),
                       NumPars = c(3,3,2),
                       df = c(125,125,126),
                       AICc=c(round(AICc(DivReg),2),round(AICc(DivReg_Temp),2),round(AICc(DivRegBase),2)),
                       dAICc=c(round(AICc(DivReg)-AICc(DivRegBase),2),round(AICc(DivReg_Temp)-AICc(DivRegBase),2),'--'))

Betatab <- data.frame(logLik= c(round(logLik(BetaReg),2),round(logLik(BetaReg_Temp),2),round(logLik(BetaRegBase),2)),
                       NumPars = c(3,3,2),
                       df = c(125,125,126),
                       AICc=c(round(AICc(BetaReg),2),round(AICc(BetaReg_Temp),2),round(AICc(BetaRegBase),2)),
                       dAICc=c(round(AICc(BetaReg)-AICc(BetaRegBase),2),round(AICc(BetaReg_Temp)-AICc(BetaRegBase),2),'--'))
                       

write.table(Richtab, file='~/Desktop/RichTab.csv',sep=',')
write.table(Divtab, file='~/Desktop/DivTab.csv',sep=',')
write.table(Betatab, file='~/Desktop/BetaTab.csv',sep=',')

L20<-subset(AlpineMaster,AlpineMaster$LIAprop<0.1)
L60<-subset(AlpineMaster,AlpineMaster$LIAprop>=0.1 & AlpineMaster$LIAprop < 0.3)
L80<-subset(AlpineMaster,AlpineMaster$LIAprop>=0.3)

plot(L20$Pcnt_Just_Ice, L20$Richness,type='p',col='red',pch=16,xlim=c(0,80))

points(L60$Pcnt_Just_Ice, L60$Richness,type='p',col='green',pch=16)

points(L80$Pcnt_Just_Ice, L80$Richness,type='p',col='purple',pch=16)
LLIA<-logit((AlpineMaster$LIAprop+0.01))
LLIA[is.na(LLIA)]<-max(LLIA,na.rm=T)

summary(lm(log(AlpineMaster$Richness) ~ LGCC * LLIA))

plot(AlpineMaster$LIAprop, AlpineMaster$Pcnt_Just_Ice)

DiffIce<-AlpineMaster$LIAprop - (AlpineMaster$Pcnt_Just_Ice/100)
plot(AlpineMaster$Pcnt_Just_Ice,DiffIce)


summary(lm(Richne ~ LGCC*DiffIce))

lm1<-lm(LCC ~ LGCC)

plot(residuals(lm1) ~ AlpineMaster$LIAprop)

summary(lm(residuals(lm1) ~ AlpineMaster$LIAprop))
