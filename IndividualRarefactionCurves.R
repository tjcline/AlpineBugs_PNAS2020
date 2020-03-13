logit<-function(x){return(log(x/(1-x)))}
RawCounts<-read.csv('~/Documents/AlpineBugs/RawCount_AggregateMidges.csv',header=T,stringsAsFactors = F)
RawCounts[is.na(RawCounts)]<-0
CountsOnly<-RawCounts[,-c(1,2)]

AlpineMaster<-read.csv('~/Documents/AlpineBugs/AlpineMasterData20191122_FewerGroups_GlacialSprings.csv',header=T,stringsAsFactors=FALSE)
IceAndSnow<-read.csv('~/Documents/AlpineBugs/LIA calcs_20191122_2.csv',header=T,stringsAsFactors=F)

AlpineMaster$Pcnt_Just_Ice <- 100 * (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])
AlpineMaster$LIAprop <- IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]
AlpineMaster$PropLoss <- 1-((IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])/(IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]))#LIA$JG.prop.loss[LIAsiteInd]
AlpineMaster$PropLoss[is.na(AlpineMaster$PropLoss)]<-0 #zeros for sites that had no glaciers at LIA

AlpineMaster$PropDiff <- (IceAndSnow$LIA_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]) - (IceAndSnow$just_ice_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)]/IceAndSnow$watershed_area_m[match(IceAndSnow$Site_Name,AlpineMaster$Site_Name)])



library(vegan)

xrare<-ceiling(CountsOnly)
raremax <- min(rowSums(xrare))
S<-specnumber(xrare)
Srare<-rarefy(x=xrare,sample=raremax)

sort(rowSums(xrare))
hist(rowSums(xrare))
mean(rowSums(xrare))

quantile(rowSums(xrare),c(0.1,0.2,0.3,0.4,0.5))

standardSample<-100
rarecurve(xrare,step=1,col='blue',cex=0.6)
abline(v=standardSample)
RarifiedRichness <- log(rarefy(xrare,sample=standardSample))

MatchingLGCC <- logit(AlpineMaster$Pcnt_Just_Ice[match(RawCounts[,1],AlpineMaster$Site_Name)]/100 + 0.01)
MatchingDRAIN<-AlpineMaster$Drainage[match(RawCounts[,1],AlpineMaster$Site_Name)]

library(lme4)

RarifiedRichnessModel<-lmer(RarifiedRichness ~ MatchingLGCC + (1|MatchingDRAIN),na.action=na.fail,REML=F,control=lmerControl(optimizer="Nelder_Mead"))
summary(RarifiedRichnessModel)


