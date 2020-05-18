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
library(RColorBrewer)
library(AICcmodavg)
library(parallel)

#Required functions
logit<-function(x){return(log(x/(1-x)))}
iLogit<-function(x){return(exp(x)/(1+exp(x)))}
Zscore<-function(x){return((x-mean(x))/sd(x))}

setwd('~/Documents/AlpineBugs')

#load required LDA results
load('AlpineBugs_LDA_Clean.Rdata')

AllTemps<-read.csv('Alpine with all temp predictions 20190906.csv',header=T,stringsAsFactors=F)

#Future temperature predictions for all sites
AlpineMaster$ts8_45_75 <- AllTemps$ts8_45_75[match(AlpineMaster$SiteCode,AllTemps$SiteCode)]
AlpineMaster$ts8_85_75 <- AllTemps$ts8_85_75[match(AlpineMaster$SiteCode,AllTemps$SiteCode)]


CC<-AlpineMaster$ColdComm
LCC<- AlpineMaster$logitColdComm
DIST<- AlpineMaster$Distance
LDIST<-log(DIST)
TEMP <- AlpineMaster$ts8_base
LGCC <- logit((AlpineMaster$Pcnt_Just_Ice)/100 + 0.01)
ELEV<- AlpineMaster$Elevation/1000

DRAIN<-AlpineMaster$DrainageCode
STREAM<-AlpineMaster$Stream_Name

#Dredge all models
bigOL <- lmer(LCC ~ LDIST+TEMP+LGCC+ELEV + LDIST:LGCC + (1|STREAM),na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"),REML=F)
bigOLdredge <- dredge(bigOL)

#Only allow 2way interactions
TwoWayModels<-get.models(bigOLdredge,subset=T)
ModTabTwoWay<-aictab(TwoWayModels)

#Select top models with AICc < 4 different than the best
TopModelsIndex<-which(ModTabTwoWay$Delta_AICc<4)
TopModels<-TwoWayModels[TopModelsIndex]
TopModTab<-aictab(TopModels)

#write.table(TopModTab,file='~/Desktop/ProjectionsModelTabs.csv',row.names=FALSE,sep=',')

#Predict regressions with future scenarios of glacier and temperature
PRED_ELEV<-ELEV
#Scenarios
#t:current,rcp4.5_75,rcp8.5_75
#g:current, 0% cover

#predictions are made for all of the top models and weighted according to model weights

ScenarioPredictions <- list()
count<-1
for(g in 1:2){#1 is modern day GCC, 2 is 0 GCC
  PRED_DIST <- LDIST
  PRED_GCC <- switch(g,logit(AlpineMaster$Pcnt_Just_Ice/100 + 0.01),logit(rep(0.01,nrow(AlpineMaster))))
  for(t in 1:2){ #1 is modern temp, 2 is without glaciers
    if(g==1){
      PRED_TEMP <- AlpineMaster$ts8_base
    }else{
      PRED_TEMP <- switch(t,AlpineMaster$ts8_45_75,AlpineMaster$ts8_85_75)
    }
    
    #MakeWeightedPrediction
    thisPred<-matrix(NA,ncol=length(LCC),nrow=nrow(TopModTab))
    for(i in 1:length(TopModels)){
      thisPred[i,]<-(predict(TopModels[[i]],newdata=data.frame(LDIST=PRED_DIST,TEMP=PRED_TEMP,LGCC=PRED_GCC,ELEV=PRED_ELEV,STREAM=STREAM),re.form=~(1|STREAM)))
      #thisPred[i,]<-(predict(TopModels[[i]],newdata=data.frame(DIST=PRED_DIST,TEMP=PRED_TEMP,LGCC=PRED_GCC,ELEV=PRED_ELEV),re.form=~0))
    
      }
    WGHT_PRED<-iLogit(as.vector(t(thisPred) %*% matrix(TopModTab$AICcWt[1:length(TopModels)],ncol=1)))

    ScenarioPredictions[[count]] <- list(t=t,g=g,pred=WGHT_PRED)
    count<-count+1
  }
}

## Calculate cumulative habitat supporting differnt levels of the cold water community
xSeq<-seq(0,1,length=25)
Ymat<-matrix(NA,nrow=length(ScenarioPredictions),ncol=(length(xSeq)-1))
CCdist<-rep(NA,(length(xSeq)-1))
for(i in 1:length(ScenarioPredictions)){
  for(j in 1:(length(xSeq)-1)){
    CCdist[j]<-sum(CC>xSeq[j])/length(CC)
    P1<-ScenarioPredictions[[i]]$pred
    Ymat[i,j]<-sum(P1>xSeq[j])/length(P1)
  }
}