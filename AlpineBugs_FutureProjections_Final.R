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


logit<-function(x){return(log(x/(1-x)))}
iLogit<-function(x){return(exp(x)/(1+exp(x)))}
Zscore<-function(x){return((x-mean(x))/sd(x))}

setwd('~/Documents/AlpineBugs')

#source('AlpineBugs_LDA_Clean.R')
load('AlpineBugs_LDA_Clean.Rdata')

AllTemps<-read.csv('Alpine with all temp predictions 20190906.csv',header=T,stringsAsFactors=F)
head(AllTemps)

length(which(p0$ColdComm >=0.8)) / nrow(p0)
length(which(p100$ColdComm >=0.8)) / nrow(p100)


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


bigOL.stream <- lmer(LCC ~ LDIST+TEMP+LGCC+ELEV + LDIST:LGCC + (1|STREAM),na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"),REML=T)
bigOL.drain <- lmer(LCC ~ LDIST+TEMP+LGCC+ELEV + LDIST:LGCC + (1|DRAIN),na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"),REML=T)
bigOL.streamDrain <- lmer(LCC ~ LDIST+TEMP+LGCC+ELEV + LDIST:LGCC + (1|STREAM/DRAIN),na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"),REML=T)
AICc(bigOL.stream)
AICc(bigOL.drain)
AICc(bigOL.streamDrain)

bigOL <- lmer(LCC ~ LDIST+TEMP+LGCC+ELEV + LDIST:LGCC + (1|STREAM),na.action=na.fail,control=lmerControl(optimizer="Nelder_Mead"),REML=F)
bigOLdredge <- dredge(bigOL)

# TwoWayModels<-get.models(bigOLdredge,
#                       (is.na(bigOLdredge[,8])))

TwoWayModels<-get.models(bigOLdredge,subset=T)

ModTabTwoWay<-aictab(TwoWayModels)
TopModelsIndex<-which(ModTabTwoWay$Delta_AICc<4)

TopModels<-TwoWayModels[TopModelsIndex]
TopModTab<-aictab(TopModels)

#write.table(TopModTab,file='~/Desktop/ProjectionsModelTabs.csv',row.names=FALSE,sep=',')



PRED_ELEV<-ELEV
#Scenarios
#t:current,rcp4.5_75,rcp8.5_75
#g:0,50,100 Loss

ScenarioPredictions <- list()
count<-1
for(g in 1:2){#1 us modern day GCC, 2 is 0 GCC
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

##
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

quartz(width=4,height=4)
#par(mfrow=c(1,2))
##plot(AlpineMaster$ColdComm ~ ScenarioPredictions[[1]]$pred,type='p',pch=16,xlab='Predicted',ylab='Observed')
#plot(logit(AlpineMaster$ColdComm) ~ logit(ScenarioPredictions[[1]]$pred),type='p',xlim=c(-15,15),pch=16,ylim=c(-15,15),xlab='Predicted',ylab='Observed')
#abline(0,1)
#summary(lm(logit(AlpineMaster$ColdComm) ~ logit(ScenarioPredictions[[1]]$pred)))
par(mar=c(4,5,1,1))
plot(logit(AlpineMaster$ColdComm) ~ logit(ScenarioPredictions[[1]]$pred),type='p',xlim=c(-16,16),pch=16,ylim=c(-16,16),xlab='Predicted relative abundance \n of cold-water community',ylab='Observed relative abundance \n of cold-water community')

abline(0,1)
#abline(lm(logit(AlpineMaster$ColdComm) ~ logit(ScenarioPredictions[[1]]$pred)),type='l')

quartz(width=2.52,height=2.5)
par(mar=c(2.5,1.75,2,1.75))
par(lheight=0.75)
par(mgp=c(0,0.2,0))
plot(xSeq[2:length(xSeq)],Ymat[1,],type='l',col='black',lwd=3,xlim=c(0,1),ylim=c(0,1),ann=F,axes=F,xaxs='i',yaxs='i')
#points(xSeq[2:50],CCdist,type='l',col='darkorange',lwd=3)
for(i in c(3,4)){
  LTY<-1
  COL<-brewer.pal(3,'Blues')[2:3][ScenarioPredictions[[i]]$t]
  points(xSeq[2:length(xSeq)],Ymat[i,],type='l',col=COL,lwd=3,lty=LTY)
}

for(i in c(3,4)){
  LTY<-3
  COL<-brewer.pal(3,'Blues')[2:3][ScenarioPredictions[[i]]$t]
  points(xSeq[2:length(xSeq)],1-(Ymat[i,]/Ymat[1,]),type='l',col=COL,lwd=3,lty=LTY)
}

mtext('Proportion of sites',2,line=0.8,cex=8/12)
mtext('Relative abundance \n of cold-water community',1,line=1.25,cex=8/12)
mtext('% loss in future habitat',4,line=0.75,cex=8/12)
axis(1,cex.axis=7/12,tck=-0.05)
axis(2,cex.axis=7/12,tck=-0.05)
axis(4,at=seq(0,1,length=6),labels=seq(0,100,length=6),cex.axis=7/12,tck=-0.05)
box(bty='l')

par(xpd=T)
legend(0,1.3,cex=8/12,bty='n',legend=c('Current','RCP 4.5 2075 (0% glaciers)','RCP 8.5 2075 (0% glaciers)'),col=c('black',brewer.pal(3,'Blues')[2:3]),lty=c(1,1,1),lwd=c(3,3,3))



#legend(0.4,0.95,bty='n',legend=c('Current GCC','0% GCC'),col=c('black',brewer.pal(3,'Blues')[3]),lty=rep(1,1),lwd=c(2,2),cex=8/12)
#legend(0.4,0.75,bty='n',legend=c('RCP 4.5 2075','RCP 8.5 2075'),lty=c(2,3),lwd=c(2,2),cex=8/12)





# thisPred<-matrix(NA,ncol=length(LCC),nrow=length(TopModels))
# for(i in 1:length(TopModels)){
#   thisPred[i,]<-(predict(TopModels[[i]],newdata=data.frame(DIST=PRED_DIST,TEMP=PRED_TEMP,GCC=PRED_GCC,ELEV=PRED_ELEV),re.form=~(1|DRAIN)))
# }
# 
# WGHT_PRED<-iLogit(as.vector(t(thisPred) %*% matrix(bigOLdredge$weight[1:length(TopModels)],ncol=1)))


# plot(CC~DIST)
# points(WGHT_PRED~DIST)
# 
# mean(CC)
# mean(WGHT_PRED)
# 
# sum(CC>0.8)/length(CC)
# sum(WGHT_PRED>0.8)/length(CC)
# 
# sum(CC<0.2)/length(CC)
# sum(WGHT_PRED<0.2)/length(CC)
# 
# hist(CC,col='#ff8c0033',ylim=c(0,60))
# hist(WGHT_PRED,add=T,col='#1e90ff33')
# 
# plot(density(CC),col='darkorange',xlim=c(0,1),lwd=3,ylim=c(0,4))
# points(density(WGHT_PRED),col='dodgerblue',lwd=3,type='l')






