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


NewlyEmerged<-read.csv('RecentlyEmergedSites.csv',header=T,stringsAsFactors=F)
AlpineMaster$RecentSite <- NewlyEmerged$X[match(AlpineMaster$SiteCode,NewlyEmerged$SiteCode)]

#HistGrad<-hist(AlpineMaster$Distance,plot=F,breaks=c(0,50,100,200,400,800,1600,3200,6400))


#### Using Latent Dirichlet Allocation to separate communities.


forceMatrixToInteger <- function(m){
    apply (m, c (1, 2), function (x) {
         (as.integer(x))
    })
}


#sortedAbundance<-SpeciesAbundance[order(AlpineMaster$ts8_base),]
rndSpecAbund<-forceMatrixToInteger(round(SpeciesAbundance*1E3))
intSpeciesAbundance<-as.simple_triplet_matrix(rndSpecAbund)


# topicNumber <- FindTopicsNumber(rndSpecAbund,topics=seq(2,10),
#                                 metrics=c("Griffiths2004", "CaoJuan2009","Arun2010", "Deveaud2014"),
#                                 method='VEM',mc.cores=4L,verbose=T)
# FindTopicsNumber_plot(topicNumber)

SEED<-2005
nComm<-2
VEM<-LDA(rndSpecAbund,k=nComm, control = list(seed = SEED,best=TRUE))
#VEM3<-LDA(rndSpecAbund,k=3, control = list(seed = SEED,best=TRUE))


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



beta_spread <- PerCommunity_PerSpecies_Prob %>%
  mutate(topic = paste0("topic", topic)) %>%
  spread(topic, beta) %>%
  filter(topic1 > .001 | topic2 > .001) %>%
  mutate(log_ratio = log2(topic2 / topic1))

BetaRat.Cold<-beta_spread %>% arrange(log_ratio) %>% slice(1:10)
BetaRat.Warm<-beta_spread %>% arrange(desc(log_ratio)) %>% slice(1:10)

Names.BetaRat.Cold<-c('Stygobromus glacialis','Gymnopais','Allomyia bifosa','Gonomyodes','Lednia tumana',
                 'Pseudokiefferiella','Prosimulium','Rhyacophila ebria','Diamesa','Thienemanniella')

Names.BetaRat.Warm<-c('Parapsyche elsis','Syndiamesa','Paraperla','Roederiodes','Diplocladius cultriger',
                 'Epeorus grandis','Pedicia','Drunella coloradensis','Epeorus deceptivus','Zapada glacier')

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

#sum(AlpineMaster$Pcnt_Total_Ice_Snow<=3)

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
for(j in c(1,3)){
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


######

p0<-AlpineMaster %>% filter(LIAprop < 0.01) #HAD NO GLACIER AT LIA
p100<-AlpineMaster %>% filter(LIAprop > 0.01 & Pcnt_Just_Ice==0) #LOST GLACIER SINCE LIA
p00<-AlpineMaster %>% filter((LIAprop < 0.01 & Pcnt_Just_Ice==0) | (LIAprop > 0.01 & Pcnt_Just_Ice == 0))
pElse<-AlpineMaster %>% filter(Pcnt_Just_Ice>0)

nrow(p0)
nrow(p100)
nrow(p00)

cl <- makeCluster(detectCores())
registerDoParallel(cl)

#WhichDataSet
G0set<-p00

P00LCC<-G0set$logitColdComm
P00ELEV<-G0set$Elevation/1000
P00DIST<-log(G0set$Distance)
P00SNOW<-logit((IceAndSnow$permanent_snow_area_m[match(G0set$Site_Name,IceAndSnow$Site_Name)]/IceAndSnow$watershed_area_m[match(G0set$Site_Name,IceAndSnow$Site_Name)])+0.01)
P00TEMP<-G0set$ts8_base
P00DRAIN<-G0set$DrainageCode
P00STREAM <- G0set$Stream_Name

p00full_DRAIN<-lmer(P00LCC ~ P00DIST+P00ELEV+P00TEMP+P00SNOW + (1|P00DRAIN),REML=T,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
p00full_STREAM<-lmer(P00LCC ~ P00DIST+P00ELEV+P00TEMP+P00SNOW + (1|P00STREAM),REML=T,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
AICc(p00full_DRAIN)
AICc(p00full_STREAM)

nrow(p00)
#DRAINAGE

p00full<-lmer(P00LCC ~ P00DIST+P00ELEV+P00TEMP+P00SNOW + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
dredge(p00full)

lmBase<-lmer(P00LCC ~ 1 +(1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmBase)

lmElv<-lmer(P00LCC ~ P00ELEV + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmElv)
dfP00elev<-data.frame(P00ELEV=seq(min(P00ELEV),max(P00ELEV),length=nrow(G0set)))
P00PRED<-dfP00elev$P00ELEV
if(F){
  PP_ELV<-iLogit(predictInterval(lmElv,newdata=dfP00elev,n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE))
}
plot(G0set$Elevation,G0set$ColdComm,pch=16,cex=1.25)


lmDist<-lmer(P00LCC ~ P00DIST + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmDist)
dfP00dist<-data.frame(P00DIST=seq(min(P00DIST),max(P00DIST),length=nrow(G0set)))
P00PRED<-dfP00dist$P00DIST
if(FALSE){
  PP_DIST<-iLogit(predictInterval(lmDist,newdata=dfP00dist,n.sims=100000,level=0.50,which='fixed',.parallel=TRUE,include.resid.var = FALSE))
}
plot(G0set$Distance,G0set$ColdComm,pch=16,cex=1.25)

lmTemp<-lmer(P00LCC ~ P00TEMP + (1|P00DRAIN),REML=F,na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"))
summary(lmTemp)
dfP00temp<-data.frame(P00TEMP=seq(min(P00TEMP),max(P00TEMP),length=nrow(G0set)))
P00PRED<-dfP00temp$P00TEMP
#PP_TEMP <- predict(lmTemp,newdata=dfP00,se.fit=T)
PP_TEMP<-iLogit(predictInterval(lmTemp,newdata=dfP00temp,n.sims=10000,level=0.95,which='fixed',.parallel=TRUE,include.resid.var = FALSE))

# plot(p100$ts8_base,p100$ColdComm,type='p',pch=16,col=p100$GradColor)
# #points(p0$ts8_base,p0$ColdComm,type='p')
# polygon(c(P00PRED[1],P00PRED,rev(P00PRED)),c(PP_TEMP[1,3],PP_TEMP[,2],rev(PP_TEMP[,3])),col=paste0(BluePalGen(9)[5],44),border=paste0(BluePalGen(9)[5],44),xpd=T)#'#1e90ff33'
# points(P00PRED,PP_TEMP[,1],type='l',col=BluePalGen(9)[5],lwd=3,xpd=T)

stopCluster(cl)

#points(P00PRED,iLogit(PP_TEMP$fit+2*PP_TEMP$se.fit),type='l',lty=3,col=BluePal(9)[5],lwd=3,xpd=T)
#points(P00PRED,iLogit(PP_TEMP$fit-2*PP_TEMP$se.fit),type='l',lty=3,col=BluePal(9)[5],lwd=3,xpd=T)


#Connect points longitudinally with a single stream
# if(FALSE){
#   plot(p0$Distance,p0$ColdComm,pch=16)
#   colInd<-6
#   count<-1
#   Streams<-unique(p0$Stream_Name)
#   for(i in 1:length(Streams)){
#     s1<-p0 %>% filter(Stream_Name == Streams[i])
#     if(nrow(s1)>1){
#       cols<-c(brewer.pal(9,'Oranges')[colInd],brewer.pal(9,'Blues')[colInd],brewer.pal(9,'Greens')[colInd],brewer.pal(9,'Greys')[colInd])
#       points(s1$Distance,s1$ColdComm,type='l',col=cols[count],lwd=3)
#       count<-count+1
#     }
#   }
# }


save(list=ls(),file='AlpineBugs_LDA_Clean.Rdata')



###Test Dean's ideas
summary(lmer(LCC ~ LGCC + (1|DRAIN),REML=F))
summary(lmer(LCC ~ TEMP + (1|DRAIN),REML=F))
summary(lmer(LCC ~ LGCC*TEMP + (1|DRAIN),REML=F))


plot(AlpineMaster$LIAprop,(GCC/100),ylab='Modern GCC',xlab='LIA GCC')



plot(LCC ~ AlpineMaster$PropDiff)

summary(lmer(LCC ~ AlpineMaster$PropDiff + (1|DRAIN),REML=F))

plot(p100$logitColdComm~p100$Pcnt_Just_Ice)

# load('AlpineBugs_LDA_Clean.Rdata')
# 
# topSpeciesOUT<-PerCommunity_PerSpecies_Prob %>%
#   group_by(topic)%>%
#   ungroup() %>%
#   arrange(topic,-beta)
# 
# write.table(topSpeciesOUT,file='~/Desktop/TopSpeciesOut.csv',sep=',')
# Site_community_prob <- tidy(VEM,matrix='gamma')
# 
# pdf('Gamma_byGroup.pdf')
# for(i in 1:nComm){
#   Gamma<-Site_community_prob %>% filter(topic==i)
#   par(mfrow=c(4,1))
#   LG<-logit(Gamma$gamma)
#   EV<-AlpineMaster$Elevation
#   DT<-AlpineMaster$Distance
#   lm1<-lm(LG~EV)
#   lm2<-lm(LG~DT)
#   P1<-predict(lm1,newdata=list(EV=seq(min(AlpineMaster$Elevation),max(AlpineMaster$Elevation))))
#   P2<-predict(lm2,newdata=list(DT=seq(min(AlpineMaster$Distance),max(AlpineMaster$Distance))))
#   plot(Gamma$gamma ~ AlpineMaster$Elevation,pch=16,col=c('dodgerblue','darkorange','darkgreen')[AlpineMaster$NewSource])
#   points(iLogit(P1)~seq(min(AlpineMaster$Elevation),max(AlpineMaster$Elevation)),type='l',lwd=2)
#   plot(Gamma$gamma ~ (AlpineMaster$Pcnt_Total_Ice_Snow),pch=16,col=c('dodgerblue','darkorange','darkgreen')[AlpineMaster$NewSource])
#   plot(Gamma$gamma ~ (AlpineMaster$Distance),pch=16,col=c('dodgerblue','darkorange','darkgreen')[AlpineMaster$NewSource])
#   points(iLogit(P2)~seq(min(AlpineMaster$Distance),max(AlpineMaster$Distance)),type='l',lwd=2)
#   plot(Gamma$gamma ~ (AlpineMaster$ts8_base),pch=16,col=c('dodgerblue','darkorange','darkgreen')[AlpineMaster$NewSource])
# }
# dev.off()








# inputData_train$Y <- DR_data (inputData_train[,1:3])
# 
# 
# Gamma<-Site_community_prob %>% filter(topic==i)
# LG<-logit(Gamma$gamma)
# EV<-AlpineMaster$Elevation
# 
# 
# 
# library(lme4)
# summary(lmer(logit(Gamma$gamma)~AlpineMaster$DistanceAlpineMaster$Distance:factor(AlpineMaster$NewSource) + (1|AlpineMaster$DrainageCode)))
# 
# 
# 
# 
# z=posterior(VEM)
# commun.plot=z$topics
# commun.spp=z$terms
# 
# 
# #plot relative abundance of component communities for each sampling unit 1,...,1000 (i.e., along the gradient)
# par(mfrow=c(1,1))
# plot(NA,NA,xlim=c(0,nrow(rndSpecAbund)),ylim=c(0,1),xlab='Gradient/Sampling units',ylab='Relative abundance')
# for (i in 1:nGroups){
#   lines(1:nrow(rndSpecAbund),commun.plot[,i],col=i)
# }
# 
# 
# plot(commun.plot[,1]~AlpineMaster$Distance,type='p')
# 
# plot(commun.plot[,1]~AlpineMaster$ts8_base,type='p')
# 
# 
# plot(commun.plot[,1]~AlpineMaster$Pcnt_Total_Ice_Snow,type='p')
# 
# summary(lm(logit(commun.plot[,2])~logDist*AlpineMaster$Glac_indx_SnowIce*AlpineMaster$ts8_base))
# 
# points(commun.plot[,]~AlpineMaster$Glac_indx_SnowIce)
# 
# #plot relative abundance of species 1,...,2000 in component community
# par(mfrow=c(nGroups,1))
# for (i in 1:4){
#   plot(1:ncol(rndSpecAbund),commun.spp[i,],type='l',col=i,xlab='Species',ylab='Relative abundance')
# }
# 
# plot(1:ncol(rndSpecAbund),commun.spp[1,],type='p')
# 
# rev(sort(commun.spp[2,]))

hist(AlpineMaster$LIAprop)

L20<-subset(AlpineMaster,AlpineMaster$LIAprop<0.3)
L60<-subset(AlpineMaster,AlpineMaster$LIAprop>=0.3 & AlpineMaster$LIAprop < 0.8)
L80<-subset(AlpineMaster,AlpineMaster$LIAprop>=0.8)

plot(L20$Pcnt_Just_Ice, L20$ColdComm,type='p',col='red',pch=16,xlim=c(0,80))

plot(L60$Pcnt_Just_Ice, L60$ColdComm,type='p',col='blue',pch=16)

plot(L80$Pcnt_Just_Ice, L80$ColdComm,type='p',col='purple',pch=16)

plot(AlpineMaster$LIAprop, AlpineMaster$Pcnt_Just_Ice)

DiffIce<-AlpineMaster$LIAprop - (AlpineMaster$Pcnt_Just_Ice/100)
plot(AlpineMaster$Pcnt_Just_Ice,DiffIce)


summary(lm(LCC ~ LGCC*DiffIce))

lm1<-lm(LCC ~ LGCC)

plot(residuals(lm1) ~ AlpineMaster$LIAprop)

summary(lm(residuals(lm1) ~ AlpineMaster$LIAprop))


AlpineMaster$PcntLoss<-100*AlpineMaster$PropLoss
AlpineNo0<-filter(AlpineMaster,PcntLoss>0)
quartz(width=4,height=3.5)
par(mar=c(4,4,1,1))
plot(AlpineNo0$ColdComm ~ AlpineNo0$PcntLoss,xpd=T,xlab='Ice loss (%)',ylab='Relative abundance of \n cold-water community',pch=16,axes=F,ann=F,yaxs='i',xaxs='i')
axis(1)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0))
box(bty='l')
mtext('Relative abundance of \n cold-water community',2,line=2)
mtext('Ice Loss (%)',1,line=2)
#abline(lm(AlpineNo0$logitColdComm ~ AlpineNo0$PcntLoss))

plot(AlpineNo0$logitColdComm ~ AlpineNo0$PropLoss*100,xlab='Ice loss (%)')
abline(lm(AlpineNo0$logitColdComm ~ AlpineNo0$PropLoss))

summary(lmer(AlpineNo0$logitColdComm ~ AlpineNo0$PropLoss + (1|AlpineNo0$Drainage),REML=F))
