rm(list=ls())
#Code for Alpine insect diversity, Glacier National Park
#Timothy Cline
#8/28/2019
library(dplyr)
library(topicmodels)
library(ldatuning)
library(slam)
library(tidytext)
library(ggplot2)
library(DirichletReg)
library(tidyr)
library(lme4)
library(MuMIn)
library(merTools)
library(doParallel)
library(parallel)
library(pscl)
library(RColorBrewer)
library(AICcmodavg)

#Required functions
logit<-function(x){return(log(x/(1-x)))}
iLogit<-function(x){return(exp(x)/(1+exp(x)))}
Zscore<-function(x){return((x-mean(x))/sd(x))}

setwd('~/Documents/AlpineBugs')
AlpineMaster<-read.csv('AlpineMasterData20191122_FewerGroups_GlacialSprings.csv',header=T,stringsAsFactors=FALSE)
IceAndSnow<-read.csv('LIA calcs_20191122_2.csv',header=T,stringsAsFactors=F)

SpeciesAbundance<-AlpineMaster[,which(colnames(AlpineMaster)=="Liodessus_affinis"):ncol(AlpineMaster)]
SpeciesAbundance<-SpeciesAbundance[,colSums(SpeciesAbundance)>0]

SpeciesPresence <- data.frame(matrix(as.numeric(SpeciesAbundance>0),nrow=nrow(SpeciesAbundance),ncol=ncol(SpeciesAbundance)))
colnames(SpeciesPresence) <- colnames(SpeciesAbundance)

AlpineMaster$Pcnt_Just_Ice <- 100 * (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])
AlpineMaster$LIAprop <- IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]
AlpineMaster$PropLoss <- 1-((IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])/(IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]))#LIA$JG.prop.loss[LIAsiteInd]
AlpineMaster$PropLoss[is.na(AlpineMaster$PropLoss)]<-0 #zeros for sites that had no glaciers at LIA

AlpineMaster$PropDiff <- (IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]) - (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])


#### Using Latent Dirichlet Allocation to separate communities.

#Function to create integer values
forceMatrixToInteger <- function(m){
    apply (m, c (1, 2), function (x) {
         (as.integer(x))
    })
}

#LDA takes counts. 
rndSpecAbund<-forceMatrixToInteger(round(SpeciesAbundance*1E3))
intSpeciesAbundance<-as.simple_triplet_matrix(rndSpecAbund)


#RUN LDA with 2 groups
SEED<-55
nComm<-2
VEM<-LDA(rndSpecAbund,k=nComm, control = list(seed = SEED,best=TRUE))

PerCommunity_PerSpecies_Prob<-tidy(VEM,matrix='beta')

topSpecies<-PerCommunity_PerSpecies_Prob %>%
  group_by(topic)%>%
  top_n(20,beta) %>%
  ungroup() %>%
  arrange(topic,-beta)

#plot topSpecies
quartz()
topSpecies %>%
  mutate(term = reorder(term, beta)) %>%
  ggplot(aes(term, beta, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free") +
  coord_flip()

z=posterior(VEM)
commun.plot=z$topics
commun.spp=z$terms


#creates a table of species that differ most between communities (topics)
beta_spread <- PerCommunity_PerSpecies_Prob %>%
  mutate(topic = paste0("topic", topic)) %>%
  spread(topic, beta) %>%
  filter(topic1 > .001 | topic2 > .001) %>%
  mutate(log_ratio = log2(topic2 / topic1))

#10 most different species from cold and warm
BetaRat.Cold<-beta_spread %>% arrange(log_ratio) %>% slice(1:10)
BetaRat.Warm<-beta_spread %>% arrange(desc(log_ratio)) %>% slice(1:10)

#Cleaning up names
Names.BetaRat.Cold<-c('Stygobromus glacialis','Gymnopais','Allomyia bifosa','Gonomyodes','Lednia tumana',
                 'Pseudokiefferiella','Prosimulium','Rhyacophila ebria','Diamesa','Thienemanniella')
Names.BetaRat.Warm<-c('Parapsyche elsis','Syndiamesa','Paraperla','Roederiodes','Diplocladius cultriger',
                 'Epeorus grandis','Pedicia','Drunella coloradensis','Epeorus deceptivus','Zapada glacier')

#ColdWater and WarmWater communities
CWC<-rev(sort(commun.spp[1,]))
WWC<-rev(sort(commun.spp[2,]))
Names.WWC<-c('Orthocladius','Tvetenia','Pagastia','Rheocricotopus','Sweltsa','Eukiefferiella','Diamesa','Megarcys',
             'Rhyacophila belona','Rhithrogena','Zapada_glacier','Zapada_columbiana','Diplocladius cultriger','Tokunagaia',
             'Cinygmula','Clinocera','Setvena bradleyi','Chaetocladius','Parorthocladius','Ameletus')
Names.CWC<-c('Diamesa','Lednia tumana','Prosimulium','Allomyia bifosa','Orthocladius','Tokunagaia','Rhyacophila ebria',
             'Polycelis','Thienemanniella','Gonomyodes','Chaetocladius','Gymnopais','Pseudodiamesa','Ameletus','Pseudokiefferiella',
             'Stilocladius','Tvetenia', 'Allomyia tripunctata','Corynoneura','Dicranota')

#Binary as to whether glacial stream or not
AlpineMaster$Glac_Y_N <- ifelse(AlpineMaster$Pcnt_Just_Ice>0,'Y','N')

#Add communities to data table
ColdComm<-commun.plot[,1]
AlpineMaster$ColdComm <- ColdComm
AlpineMaster$logitColdComm<-logit(ColdComm)
WarmComm<-commun.plot[,2]



LCC <-AlpineMaster$logitColdComm
DIST <- AlpineMaster$Distance
ELEV<-AlpineMaster$Elevation
LDIST <- log(DIST)
GCC <- AlpineMaster$Pcnt_Just_Ice
LGCC <- logit(GCC/100+0.01)
TEMP <- AlpineMaster$ts8_base
GYN <- as.numeric(AlpineMaster$Glac_Y_N=='Y')
DRAIN <- AlpineMaster$DrainageCode
STREAM <- AlpineMaster$Stream_Name
nL<-length(LCC)

#Regression
#1Distance
#2Elev
#3GCC
#4Temp

if(TRUE){
for(j in c(1,2,3,4)){
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  VAR <- switch(j,LDIST,ELEV,LGCC,TEMP)
  if(j!=3) lmIn<-lmer(LCC ~  VAR + GYN + VAR:GYN + (1|DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
  if(j==3) lmIn<-lmer(LCC ~  VAR + (1|DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
  VarDredge<-dredge(lmIn);VarDredge
  TopModels<-get.models(VarDredge,delta<2)
  TopModTab<-aictab(TopModels)
  
  pVAR <- seq(min(VAR),max(VAR),length=nL)
  pLGR<-seq(min(LGCC),max(LGCC),length=nL)
  pY <- rep(1,nL)
  pN <- rep(0,nL)
  dfY<-data.frame(VAR=pVAR,GYN=pY)
  dfN<-data.frame(VAR=pVAR,GYN=pN)
  
  #dfY<-data.frame(VAR=pVAR,LGCC=rep(mean(LGCC),nL))
  #dfN<-data.frame(VAR=pVAR,LGCC=rep(min(LGCC,nL)))
  
  thisPredG<-list()
  thisPredNG<-list()
  for(i in 1:length(TopModels)){
    thisPredG[[i]]<-predictInterval(TopModels[[i]],newdata=dfY,n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)
    thisPredNG[[i]]<-predictInterval(TopModels[[i]],newdata=dfN,n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE)
  }
  
  thisPredG_W<-matrix(NA,nrow=length(TopModels),ncol=length(thisPredG[[1]][,1]))
  thisPredNG_W<-matrix(NA,nrow=length(TopModels),ncol=length(thisPredG[[1]][,1]))
  
  thisUpperG_W<-matrix(NA,nrow=length(TopModels),ncol=length(thisPredG[[1]][,1]))
  thisUpperNG_W<-matrix(NA,nrow=length(TopModels),ncol=length(thisPredG[[1]][,1]))
  
  thisLowerG_W<-matrix(NA,nrow=length(TopModels),ncol=length(thisPredG[[1]][,1]))
  thisLowerNG_W<-matrix(NA,nrow=length(TopModels),ncol=length(thisPredG[[1]][,1]))
  
  for(i in 1:length(TopModels)){
    thisPredG_W[i,]<-thisPredG[[i]]$fit
    thisPredNG_W[i,]<-thisPredNG[[i]]$fit
    
    thisUpperG_W[i,]<-thisPredG[[i]]$upr
    thisUpperNG_W[i,]<-thisPredNG[[i]]$upr
    
    thisLowerG_W[i,]<-thisPredG[[i]]$lwr
    thisLowerNG_W[i,]<-thisPredNG[[i]]$lwr
  }
  
  FullPred_G<-iLogit(t(thisPredG_W) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1))
  FullPred_NG<-iLogit(t(thisPredNG_W) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1))
  
  FullLower_G<-iLogit(t(thisLowerG_W) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1))
  FullLower_NG<-iLogit(t(thisLowerNG_W) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1))
  
  FullUpper_G<-iLogit(t(thisUpperG_W) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1))
  FullUpper_NG<-iLogit(t(thisUpperNG_W) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1))
  
  assign(x=switch(j,'PredG_Dist','PredG_Elev','PredG_GCC','PredG_Temp'),value=cbind(FullLower_G,FullPred_G,FullUpper_G))
  assign(x=switch(j,'PredNG_Dist','PredNG_Elev','PredNG_GCC','PredNG_Temp'),value=cbind(FullLower_NG,FullPred_NG,FullUpper_NG))
  assign(x=switch(j,'TopModels_Dist','TopModels_Elev','TopModels_GCC','TopModels_Temp'),value=TopModels)
  assign(x=switch(j,'Lm_Dist','Lm_Elev','Lm_GCC','Lm_Temp'),value=lmIn)
  assign(x=switch(j,'Dredge_Dist','Dredge_Elev','Dredge_GCC','Dredge_Temp'),value=VarDredge)
  
  stopCluster(cl)
}
}


###### #Regression analysis for sites without glaciers

p0<-AlpineMaster %>% filter(LIAprop < 0.01) #Sites that had NO glacier at LIA
p100<-AlpineMaster %>% filter(LIAprop > 0.01 & Pcnt_Just_Ice==0) #Sites that LOST glacier since LIA
p00<-AlpineMaster %>% filter((LIAprop < 0.01 & Pcnt_Just_Ice==0) | (LIAprop > 0.01 & Pcnt_Just_Ice == 0)) #Sites that either lost or had no glacier (combination of previous 2)
pElse<-AlpineMaster %>% filter(Pcnt_Just_Ice>0) #Sites with Glaciers

cl <- makeCluster(detectCores())
registerDoParallel(cl)

G0set<-p00

P00LCC<-G0set$logitColdComm
P00ELEV<-G0set$Elevation/1000
P00DIST<-log(G0set$Distance)
P00SNOW<-logit((IceAndSnow$permanent_snow_area_m[match(G0set$Site_Name,IceAndSnow$Site_Name)]/IceAndSnow$watershed_area_m[match(G0set$Site_Name,IceAndSnow$Site_Name)])+0.01)
P00TEMP<-G0set$ts8_base
P00DRAIN<-G0set$DrainageCode
P00STREAM <- G0set$Stream_Name

#Test for random effects
p00full_DRAIN<-lmer(P00LCC ~ P00DIST+P00ELEV+P00TEMP+P00SNOW + (1|P00DRAIN),REML=T,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
p00full_STREAM<-lmer(P00LCC ~ P00DIST+P00ELEV+P00TEMP+P00SNOW + (1|P00STREAM),REML=T,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
AICc(p00full_DRAIN)
AICc(p00full_STREAM)
#DRAINAGE

lmBase<-lmer(P00LCC ~ 1 +(1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmBase)

lmElv<-lmer(P00LCC ~ P00ELEV + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmElv)
dfP00elev<-data.frame(P00ELEV=seq(min(P00ELEV),max(P00ELEV),length=nrow(G0set)))
P00PRED<-dfP00elev$P00ELEV
if(F){
  PP_ELV<-iLogit(predictInterval(lmElv,newdata=dfP00elev,n.sims=10000,level=0.95,which='fixed',.parallel=TRUE))
}
plot(G0set$Elevation,G0set$ColdComm,pch=16,cex=1.25)


lmDist<-lmer(P00LCC ~ P00DIST + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmDist)
dfP00dist<-data.frame(P00DIST=seq(min(P00DIST),max(P00DIST),length=nrow(G0set)))
P00PRED<-dfP00dist$P00DIST
if(FALSE){
  PP_DIST<-iLogit(predictInterval(lmDist,newdata=dfP00dist,n.sims=100000,level=0.50,which='fixed',.parallel=TRUE))
}
plot(G0set$Distance,G0set$ColdComm,pch=16,cex=1.25)

lmTemp<-lmer(P00LCC ~ P00TEMP + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmTemp)
dfP00temp<-data.frame(P00TEMP=seq(min(P00TEMP),max(P00TEMP),length=nrow(G0set)))
P00PRED<-dfP00temp$P00TEMP

PP_TEMP<-iLogit(predictInterval(lmTemp,newdata=dfP00temp,n.sims=10000,level=0.95,which='fixed',.parallel=TRUE))

stopCluster(cl)


save(list=ls(),file='AlpineBugs_LDA_Clean.Rdata')

