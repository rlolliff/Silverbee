#Olliff Yang & Mesler 2018. The potential for phenological mismatch between a perennial herb and its ground nesting bee pollinator. AOB Plants 2018

#setup with overlapW/percentiles.r and MS- final regressions, cleaned script - july 2017.r
#####SET UP #####
#libraries:
library(dplyr)
library(ggplot2)
library(lme4)
library(colorRamps)  
library(grDevices)
#library(gridExtra) # for contour plot at the end (SI Fig 5)

#datasheets:
AggTemp<-read.csv("DuneTemperatures_Aggregations_tidyr.csv")
LathyTemp<-read.csv("Lathyrus_temp_tidyr.csv")
MSdata<-read.csv("All_FullPhenology.csv") #with new "UniquePlot" column, Unique plot <=17 are bee plots, >=18 are flower plots

#Functions: 
error.bar <- function(x, y, upper, lower=upper, length=0.05,...){  #for graphing error bars
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

se <- function(x) sqrt(var(x)/length(x)) #for clculating standard errors


overlapTOT<-function(fromb, fromf, tob, tof) { #overlap function- will calculate overlap - accounts for all scenarios! (ie. flowers completeley overlapped, vs bees)
  if(fromb <= fromf) { 
    if(tof > tob) overlap<-tob-fromf 
    else overlap<-tof-fromf 
  } 
  else { 
    if(tob > tof) overlap<-tof-fromb 
    else overlap<-tob - fromb 
  } 
  if(overlap < 0) overlap<-0 
  return(overlap) 
}   


#Data set up: 
df1<-MSdata %>%
  filter(NUMBER>=1) %>% #all dates with ACTIVE NESTS OR OPEN flowers recorded
  group_by(UniquePlot)  %>%  #don't need type or site actually if unique plot is used
  mutate(from = min(DOY), to = max(DOY))  %>% #from =start date, to= end date
  group_by(UniquePlot) %>% 
  filter(NUMBER==max(NUMBER))%>%
  mutate(peakDOY=DOY)  %>%  
  select(UniquePlot, from, to, peakDOY) 
df1.1<-merge(MSdata, df1, by="UniquePlot", all.y = TRUE, all.x = TRUE)#to add back all of the columns

df50<-MSdata %>%  
  filter(CUM_PROPORTION>=0.5)%>% #to calculate date of 50% active/flowering at each plot
  group_by(UniquePlot) %>% #SiteAsNumber not needed because UniquePlot is unique ID
  summarize(fifty= min(DOY)) 

df<-merge(df1.1, df50, by="UniquePlot", all.y = TRUE, all.x = TRUE)

AveTM<-AggTemp %>% #to Average points taken across a plot...
  group_by(SITE, PLOT, DOY)  %>% # select each date to average on each date
  summarise(ONEcm_T=mean(ONEcm_T, na.rm=TRUE),TENcm_T = mean(TENcm_T, na.rm=TRUE ), TWENTYFIVEcm_T= mean(TWENTYFIVEcm_T, na.rm=TRUE),  SOILMOIST= mean(SOILMOIST, na.rm=TRUE))


MinMaxTemp<-AggTemp %>% 
  group_by(SITE, PLOT)  %>% 
  summarise(minONEcm_T=min(ONEcm_T, na.rm=TRUE), minTENcm_T = min(TENcm_T, na.rm=TRUE ), 
            minTWENTYFIVEcm_T= min(TWENTYFIVEcm_T, na.rm=TRUE), minSOILMOIST= min(SOILMOIST, na.rm=TRUE), 
            maxONEcm_T=max(ONEcm_T, na.rm=TRUE),maxTENcm_T = max(TENcm_T, na.rm=TRUE ), 
            maxTWENTYFIVEcm_T= max(TWENTYFIVEcm_T, na.rm=TRUE), maxSOILMOIST= max(SOILMOIST, na.rm=TRUE))



TotAve<-AveTM%>%
  filter(DOY>=60 & DOY <=151)  %>% #Modify to specify month/dates - right now mar- May #January= 1-31, Feb= 32-59, March =60-90, Aril= 91-120, May= 121-151, Jun= 152-181 (in 2013, and all non leap years (in 2016, 2020, etc. add 1 to Feb end and all afterwards))
  group_by(SITE, PLOT) %>% 
  summarise(ONEcm_T=mean(ONEcm_T, na.rm=TRUE),TENcm_T = mean(TENcm_T, na.rm=TRUE ), 
            TWENTYFIVEcm_T= mean(TWENTYFIVEcm_T, na.rm=TRUE), SOILMOIST= mean(SOILMOIST, na.rm=TRUE)) %>%
  mutate(TYPE="Bee")

TotB<-merge(MinMaxTemp, TotAve, by=c("SITE", "PLOT"), all.y = TRUE, all.x = TRUE)


MinMaxTempL <-LathyTemp %>% #
  group_by(SITE, PLOT)  %>% 
  summarise(minONEcm_T=min(ONEcm_T, na.rm=TRUE), minTENcm_T = min(TENcm_T, na.rm=TRUE ),
            minTWENTYFIVEcm_T= min(TWENTYFIVEcm_T, na.rm=TRUE), minSOILMOIST= min(SOILMOIST, na.rm=TRUE),  
            maxONEcm_T=max(ONEcm_T, na.rm=TRUE), maxTENcm_T = max(TENcm_T, na.rm=TRUE ), 
            maxTWENTYFIVEcm_T= max(TWENTYFIVEcm_T, na.rm=TRUE), maxSOILMOIST= max(SOILMOIST, na.rm=TRUE))

TotAveL<-LathyTemp%>%
  filter(DOY>=60 & DOY <=151)  %>% #specify month/dates - right now Feb- May #January= 1-31, Feb= 32-59, March =60-90, Aril= 91-120, May= 121-151, Jun= 152-181 (in 2013, and all non leap years (in 2016, 2020, etc. add 1 to Feb end and all afterwards))
  group_by(SITE, PLOT) %>% 
  summarise(ONEcm_T=mean(ONEcm_T, na.rm=TRUE),TENcm_T = mean(TENcm_T, na.rm=TRUE ), TWENTYFIVEcm_T= mean(TWENTYFIVEcm_T, na.rm=TRUE), SOILMOIST= mean(SOILMOIST, na.rm=TRUE))%>%
  mutate(TYPE="Flower")

TotL<-merge(MinMaxTempL, TotAveL, by=c("SITE", "PLOT"), all.y = TRUE, all.x = TRUE)

FullTempsa<-merge(TotB, TotL, by=c("TYPE", "SITE", "PLOT", "minONEcm_T", "maxONEcm_T", "minTENcm_T", "maxTENcm_T", 
                                   "minTWENTYFIVEcm_T", "maxTWENTYFIVEcm_T", "minSOILMOIST", "maxSOILMOIST", 
                                   "ONEcm_T", "TENcm_T", "TWENTYFIVEcm_T",  "SOILMOIST"), all.y = TRUE, all.x = TRUE)

df_mean<-df %>%   #to calculate "mean date" at each plot
  group_by(SiteAsNumber, UniquePlot, TYPE) %>%
  summarise(from=mean(from), to=mean(to), fifty=mean(fifty))%>% #have to do this to make into single values
  group_by(SiteAsNumber, UniquePlot) %>%
  summarise(MeanFromTo= ((from+to)/2), fifty=fifty, TYPE=TYPE) %>%
  group_by(SiteAsNumber,UniquePlot)  %>% # have to repeat otherwise type gets left behind. Not sure why...
  summarise(MeanDate=((MeanFromTo + fifty)/2), TYPE)

df_A<-df_mean%>% 
  group_by(SiteAsNumber, TYPE) %>%
  summarize(SiteMean= mean(MeanDate))

head(df_A)
library("reshape2") #need reshape for the dcast function
test3<-df_A%>% 
  dcast(SiteAsNumber ~ TYPE, value.var = "SiteMean", fill = 0) #making a small data frame so that bee mean and flower mean can be in seperate columns, and filled in for entire site 
df_B<-merge(df_mean, test3,  by="SiteAsNumber", all.y = TRUE, all.x = TRUE) #merge should create the columns needed (bee mean and fl mean for each site in seperate coulmns)

#Now to subtract Mean date from opposite type Site Mean: 

df_ASYN<-df_B  %>% 
  mutate( diff= if_else(TYPE=="Bee", (MeanDate-Flower) ,if_else(TYPE=="Flower", (MeanDate-Bee), 999999))) %>%  # have to add in something for "false", make sure there are no 999999's in final version. If there are, something is wrong
  mutate(Asynchrony=abs(diff)) #%>%  # abs(diff) - taking absolute value of difference

dfAll<-full_join(df,df_ASYN, by=c("TYPE", "SiteAsNumber", "UniquePlot"),  all.y = TRUE, all.x = TRUE)
FullDF<-merge(dfAll, FullTempsa, by=c("TYPE", "SITE", "PLOT"),  all.y = TRUE, all.x = TRUE) 

#just bees: 
Bdf<-FullDF %>%
  filter(TYPE=="Bee")
#just flowers 
Fdf<-FullDF %>%
  filter(TYPE=="Flower")




#calculating the overlap####
BeeToFrom<-MSdata %>%
  filter(NUMBER>=1, TYPE=="Bee") %>% #all dates with ACTIVE NESTS OR OPEN flowers recorded
  group_by(SiteAsNumber, UniquePlot)  %>% #, PLOT
  summarise(fromb = min(DOY), tob = max(DOY))%>%  #1st and last dates
  mutate(totalBee=(tob-fromb))  #total # days active

SiteB<-BeeToFrom %>%
  group_by(SiteAsNumber)%>%
  mutate(startB=min(fromb), endB=max(tob)) %>%
  ungroup()%>%
  distinct(SiteAsNumber,startB, endB) #USE "DISTINCT" (VS. "SELECT")to remove repeated lines! numbers should not change

FlToFrom<-MSdata %>%
  filter(NUMBER>=1, TYPE=="Flower") %>% #all dates with ACTIVE NESTS OR OPEN flowers recorded
  group_by(SiteAsNumber, UniquePlot)  %>% #, 
  summarise(fromf = min(DOY), tof = max(DOY)) %>% #1st and last dates
  mutate(totalFl=(tof-fromf)) #total # days active - LENGTH! (CHANGE NAME?)

SiteF<-FlToFrom %>%
  group_by(SiteAsNumber)%>%
  mutate(startF=min(fromf), endF=max(tof))%>%
  ungroup()%>%
  distinct(SiteAsNumber,startF, endF) #USE "DISTINCT" (VS. "SELECT")to remove repeated lines! numbers should not change

testf<-merge(FlToFrom, SiteB, by=c("SiteAsNumber")) %>% #this is now set to calculate overlap of each plot with overall site timing of the other
  group_by(SiteAsNumber)

testb<-merge(BeeToFrom, SiteF, by=c("SiteAsNumber")) %>% 
  group_by(SiteAsNumber)

##for individual plot overlaps (based on surrounding site "resources")
bov<-testb%>%
  group_by(SiteAsNumber, UniquePlot)%>%
  summarise(doverlap=overlapTOT(fromb, startF, tob, endF), total=totalBee) %>% #overlapTOT function created at the top of script
  mutate(Overlap=(doverlap/total), percentOverlap=((doverlap/total*100)))  
head(bov)   #overlap of bees based on site

fov<-testf%>%
  group_by(SiteAsNumber, UniquePlot)%>%
  summarise(doverlap=overlapTOT(fromf, startB, tof, endB), total=totalFl) %>% 
  mutate(Overlap=(doverlap/total), percentOverlap=((doverlap/total*100)))  #overlap of flowers based on site 

df_OV<-full_join(bov,fov, by=c("SiteAsNumber", "UniquePlot", "doverlap", "total" , "Overlap", "percentOverlap" ))

df_OV<-df_OV %>%
  mutate(TYPE = ifelse(UniquePlot<=17, "Bee", ifelse(UniquePlot>17, "Flower","NA"))) #adding "type back into the dataframe

unique(df_OV$TYPE) #make sure there are no NA's - if do, ifelse above is incorrect


####Manuscript Figures #####
#Figure 3



#fig 3 panel 2
OVMean <-aggregate(df_OV$percentOverlap, by=list(df_OV$TYPE), 
                   FUN=mean, na.rm=TRUE)
OVSterr <-aggregate(df_OV$percentOverlap, by=list(df_OV$TYPE), 
                    FUN=se) # see top of code for se function

OVStdev <-aggregate(df_OV$percentOverlap, by=list(df_OV$TYPE), 
                    FUN=sd, na.rm=TRUE)

segX<-(OVMean$x - OVSterr$x * 2)
segY<-(OVMean$x + OVSterr$x * 2)

segments(segX, segY, x1=segX)

beex<-c(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7) #for adding points to bargraph - this is the x value modify for proper alignment with bargraph
flx<-c(1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9) #for adding points to bargraph -  this is the x value modify for proper alignment with bargraph


#fig 3 panel 3 - asynchrony
head(df_ASYN)
ASYNMean <-aggregate(df_ASYN$Asynchrony, by=list(df_ASYN$TYPE), 
                     FUN=mean, na.rm=TRUE)
ASYNSterr <-aggregate(df_ASYN$Asynchrony, by=list(df_ASYN$TYPE), 
                      FUN=se)# se function at top of script

ASYNStdev <-aggregate(df_ASYN$Asynchrony, by=list(df_ASYN$TYPE), 
                      FUN=sd, na.rm=TRUE)

segX<-(ASYNMean$x - ASYNSterr$x * 2)
segY<-(ASYNMean$x + ASYNSterr$x * 2)


#To make figure:
par(mfrow=c(1,3), tck=0.03, mar=c(5,5,2,1)+0.1)  #includes outward ticks, as requested by AOB plants
plot(df$PROPORTION~df$DOY, pch=c(20, 18)[as.numeric(df$TYPE)], col=c("dark grey", "dimgray")[as.numeric(df$TYPE)],xlim=c(85,185), xlab=c("Day of Year"),ylab=c("Proportion nesting/flowering")) #jitter(df$DOY,3) #jitter doesn't really help
legend(75,0.62, c("(a)"), box.lty = 0)

PLOT<-barplot(OVMean$x,, #score is # of buds and flowers and fruit (total flower production as of April 2015))
              xlab=" ", ylab="Percentage Overlap", col=c( "dark grey", "dimgray"),
              names.arg=c("Bees","Flowers"), ylim=c(0,115))
error.bar(PLOT, OVMean$x,OVSterr$x) #must have plot, Means & CI; Error bar function at top of code
points(beex,df_OV$percentOverlap[df_OV$TYPE=="Bee"], pch=20)
points(flx,jitter(df_OV$percentOverlap[df_OV$TYPE=="Flower"]), pch=20) # jittered to show most points occur at or near 100%
legend(0,113, c("(b)"), box.lty = 0)
abline(0,0, lwd=1)

PLOT<-barplot(ASYNMean$x,, #score is # of buds and flowers and fruit (total flower production as of April 2015))
              xlab=" ", ylab="Asynchrony", col=c("dark grey", "dimgray"), 
              names.arg=c("Bees","Flowers"), ylim=c(-1, 23))
error.bar(PLOT, ASYNMean$x,ASYNSterr$x) #must have plot, Means & CI, Error bar function at top of code
abline(0,0, lwd=1)
legend(0,21, c("(c)"), box.lty = 0)
points(beex,df_ASYN$Asynchrony[df_ASYN$TYPE=="Bee"], pch=20)
points(flx,df_ASYN$Asynchrony[df_ASYN$TYPE=="Flower"], pch=20)

#Figure 4 - 3x3 panel
par(mfrow = c(3,3), mar=c(4, 5, 1, 1), + 0.1)#it goes c(bottom, left, top, right) 

#25 max start
Bs25maxXSM<-(lmer(from~maxTWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Bdf))
summary(Bs25maxXSM)#extract info for tables
plot(Bdf$from~Bdf$maxTWENTYFIVEcm_T, pch=20,  xlab=c(""), ylab=c("Start Date"),cex.lab=1.5, col="darkgrey") 
# code from Meagan F. Oldfather - using predict to plot lines with marginalslope estimates
par.val<-seq(min(Bdf$maxTWENTYFIVEcm_T, na.rm=T), max(Bdf$maxTWENTYFIVEcm_T, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Bdf$SOILMOIST, na.rm=T)))
colnames(new.data)<-c("maxTWENTYFIVEcm_T", "SOILMOIST") # make to match names in models
new.data$predictions<-predict(object=Bs25maxXSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$maxTWENTYFIVEcm_T, new.data$predictions)

#25 ave start
Fs25SM<-(lmer(from~TWENTYFIVEcm_T + SOILMOIST+ (1|SITE), REML = F, data=Fdf)) 
summary(Fs25SM)
plot(Fdf$from~Fdf$TWENTYFIVEcm_T,1, pch=18, col="dim grey", xlab=c(""), ylab=c(""),cex.lab=1.5) 
par.val<-seq(min(Fdf$TWENTYFIVEcm_T, na.rm=T), max(Fdf$TWENTYFIVEcm_T, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Fdf$SOILMOIST, na.rm=T)))
colnames(new.data)<-c("TWENTYFIVEcm_T", "SOILMOIST") # make to match names in models
new.data$predictions<-predict(object=Fs25SM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$TWENTYFIVEcm_T, new.data$predictions)




#SM start


plot(Bdf$from~Bdf$SOILMOIST , pch=20, ylim=c(83, 147), xlim =c (0.054,0.1055), xlab=c(""), ylab=c(""), cex.lab=1.5, col="darkgrey")
points(Fdf$from~Fdf$SOILMOIST, pch=18, col="dim grey")



Fs25SM<-(lmer(from~TWENTYFIVEcm_T + SOILMOIST+ (1|SITE), REML = F, data=Fdf)) 
par.val<-seq(min(Fdf$SOILMOIST, na.rm=T), max(Fdf$SOILMOIST, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Fdf$TWENTYFIVEcm_T, na.rm=T)))
colnames(new.data)<-c("SOILMOIST", "TWENTYFIVEcm_T") # make to match names in models
summary(Fs25SM)
head(new.data)
new.data$predictions<-predict(object=Fs25SM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$SOILMOIST, new.data$predictions)

Bs25maxXSM<-(lmer(from~maxTWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Bdf)) 
par.val<-seq(min(Bdf$SOILMOIST, na.rm=T), max(Bdf$SOILMOIST, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Bdf$maxTWENTYFIVEcm_T, na.rm=T)))
colnames(new.data)<-c("SOILMOIST", "maxTWENTYFIVEcm_T") # make to match names in models
summary(Bs25maxXSM)
head(new.data)
new.data$predictions<-predict(object=Bs25maxXSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$SOILMOIST, new.data$predictions)








##peak dates: 
#25 max peak
plot(Bdf$peakDOY~Bdf$maxTWENTYFIVEcm_T, pch=20, xlab=c(""), cex.lab=1.5, col="darkgrey", ylab=c("Peak Date"), cex.main=1.5)

Bp25maxXSM<-(lmer(peakDOY~maxTWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Bdf))
summary(Bp25maxXSM)
par.val<-seq(min(Bdf$maxTWENTYFIVEcm_T, na.rm=T), max(Bdf$maxTWENTYFIVEcm_T, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Bdf$SOILMOIST, na.rm=T)))
colnames(new.data)<-c("maxTWENTYFIVEcm_T", "SOILMOIST") # make to match names in models
new.data$predictions<-predict(object=Bp25maxXSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$maxTWENTYFIVEcm_T, new.data$predictions)


#25 ave peak
plot(Fdf$peakDOY~Fdf$TWENTYFIVEcm_T,1, pch=18, col="dim grey", xlab=c(""),cex.lab=1.5, ylab=c("")) #"dark grey", "dimgray"
Fs25xSM<-(lmer(peakDOY~TWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Fdf)) 
par.val<-seq(min(Fdf$TWENTYFIVEcm_T, na.rm=T), max(Fdf$TWENTYFIVEcm_T, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Fdf$SOILMOIST, na.rm=T)))
colnames(new.data)<-c("TWENTYFIVEcm_T", "SOILMOIST") # make to match names in models
summary(Fs25xSM)
head(new.data)
new.data$predictions<-predict(object=Fs25xSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$TWENTYFIVEcm_T, new.data$predictions)

#SM peak
plot(Bdf$peakDOY~Bdf$SOILMOIST , pch=20, ylim=c(120, 161), xlim =c (0.054,0.1055), xlab=c(""), ylab=c(""),  col="darkgrey")
points(Fdf$peakDOY~Fdf$SOILMOIST, pch=18, col="dim grey")

Fp25xSM<-(lmer(peakDOY~TWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Fdf)) 
par.val<-seq(min(Fdf$SOILMOIST, na.rm=T), max(Fdf$SOILMOIST, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Fdf$TWENTYFIVEcm_T, na.rm=T)))
colnames(new.data)<-c("SOILMOIST", "TWENTYFIVEcm_T") # make to match names in models
summary(Fp25xSM)
head(new.data)
new.data$predictions<-predict(object=Fp25xSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$SOILMOIST, new.data$predictions)
head(Bdf)
Bp25maxXSM<-(lmer(peakDOY~maxTWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Bdf)) 
par.val<-seq(min(Bdf$SOILMOIST, na.rm=T), max(Bdf$SOILMOIST, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Bdf$maxTWENTYFIVEcm_T, na.rm=T)))
colnames(new.data)<-c("SOILMOIST", "maxTWENTYFIVEcm_T") # make to match names in models
summary(Bp25maxXSM)
head(new.data)
new.data$predictions<-predict(object=Bp25maxXSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$SOILMOIST, new.data$predictions)



#25 max end
plot(Bdf$to~Bdf$maxTWENTYFIVEcm_T, pch=20, xlab=c("25cm Max Temperature"),cex.lab=1.5, ylab=c("End Date"), col="darkgrey")

Be25maxXSM<-(lmer(to~maxTWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Bdf))
# code from Meagan- using predict!
par.val<-seq(min(Bdf$maxTWENTYFIVEcm_T, na.rm=T), max(Bdf$maxTWENTYFIVEcm_T, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Bdf$SOILMOIST, na.rm=T)))
colnames(new.data)<-c("maxTWENTYFIVEcm_T", "SOILMOIST") # make to match names in models
summary(Be25maxXSM)
head(new.data)
new.data$predictions<-predict(object=Be25maxXSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$maxTWENTYFIVEcm_T, new.data$predictions)

#25 ave end
plot(Fdf$to~Fdf$TWENTYFIVEcm_T,1, pch=18, col="dim grey", xlab=c("25cm Ave Temperature"),cex.lab=1.5, ylab=c("")) #"dark grey", "dimgray"

Fe25xSM<-(lmer(to~TWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Fdf)) 
par.val<-seq(min(Fdf$TWENTYFIVEcm_T, na.rm=T), max(Fdf$TWENTYFIVEcm_T, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Fdf$SOILMOIST, na.rm=T)))
colnames(new.data)<-c("TWENTYFIVEcm_T", "SOILMOIST") # make to match names in models
summary(Fe25xSM)
head(new.data)
new.data$predictions<-predict(object=Fe25xSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$TWENTYFIVEcm_T, new.data$predictions)



#SM end
plot(Bdf$to~Bdf$SOILMOIST , pch=20, ylim=c(150, 187), xlim =c (0.054,0.1055), xlab=c("Ave Soil Moist"), cex.lab=1.5, ylab=c(""), col="darkgrey")
points(Fdf$to~Fdf$SOILMOIST, pch=18, col="dim grey")
Fe25xSM<-(lmer(to~TWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Fdf)) 
par.val<-seq(min(Fdf$SOILMOIST, na.rm=T), max(Fdf$SOILMOIST, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Fdf$TWENTYFIVEcm_T, na.rm=T)))
colnames(new.data)<-c("SOILMOIST", "TWENTYFIVEcm_T") # make to match names in models
summary(Fe25xSM)
head(new.data)
new.data$predictions<-predict(object=Fe25xSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$SOILMOIST, new.data$predictions)
head(Bdf)
Be25maxXSM<-(lmer(to~maxTWENTYFIVEcm_T*SOILMOIST+ (1|SITE), REML = F, data=Bdf)) 
par.val<-seq(min(Bdf$SOILMOIST, na.rm=T), max(Bdf$SOILMOIST, na.rm=T), length.out=1000)
new.data<- data.frame(expand.grid(par.val,mean(Bdf$maxTWENTYFIVEcm_T, na.rm=T)))
colnames(new.data)<-c("SOILMOIST", "maxTWENTYFIVEcm_T") # make to match names in models
summary(Be25maxXSM)
head(new.data)
new.data$predictions<-predict(object=Be25maxXSM, newdata=new.data, re.form=NA, type = "response")
lines(new.data$SOILMOIST, new.data$predictions)


#Table 1 - Best abiotic variables - made in word
#Table 2 - Temp*Moist  #Linear regression coefficints / comparisons  - made in word


#SI2 - visitation 
vis<-read.csv("visitation_tot.csv")
EJapproved <- theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line=element_line(color="black"), axis.ticks=element_line(color="black"), legend.position="none", axis.title.y=element_text(vjust=0.3), axis.title.x=element_text(vjust=-0.5), axis.text.x = element_text(face='italic'))

# Make plot -
ggplot(vis, aes(x=Species, y=Number)) + geom_bar(stat="identity", fill="black") + ylab("Number of Visits") + EJapproved


#SI- Figure 5 - contour plot of interactions 
# call grid extra library : library(gridExtra)
# set up grid of (x,y) values
SM <- seq(0.085, 0.12, length.out=100) # from min & max values of the data
TEMP <- seq(20,30, length.out=100)
ggs <- data.frame(expand.grid(SOILMOIST=SM, maxTWENTYFIVEcm_T=TEMP))
ggp <- data.frame(expand.grid(SOILMOIST=SM, maxTWENTYFIVEcm_T=TEMP))
gge <- data.frame(expand.grid(SOILMOIST=SM, maxTWENTYFIVEcm_T=TEMP))

Bs25maxXSMf<-(lm(from~maxTWENTYFIVEcm_T*SOILMOIST, data=Bdf)) # models with random effects removed
Bp25maxXSMf<-(lm(peakDOY~maxTWENTYFIVEcm_T*SOILMOIST, data=Bdf))
Be25maxXSMf<-(lm(to~maxTWENTYFIVEcm_T*SOILMOIST, data=Bdf))

# prediction from the linear model
ggs$DOY <-predict(Bs25maxXSMf,newdata=ggs)
ggp$DOY <-predict(Bp25maxXSMf,newdata=ggp)
gge$DOY <-predict(Be25maxXSMf,newdata=gge)
head(ggs)
# contour plot - interactions 

jet.colors <- colorRampPalette(matlab.like(9))

Bs<-ggplot(ggs, aes(x=SOILMOIST, y=maxTWENTYFIVEcm_T, z=DOY))+
  stat_contour(aes(color=..level..),binwidth=10, size=2)+
  scale_color_gradientn(colours=jet.colors(8), limits=c(-50, 365)) +labs( x= " " ,  y= "25 cm Max Temperature", title = "  (a)                         START") + EJapproved + theme(legend.position="none")

Bp<-ggplot(ggp, aes(x=SOILMOIST, y=maxTWENTYFIVEcm_T, z=DOY))+
  stat_contour(aes(color=..level..),binwidth=10, size=2)+
  scale_color_gradientn(colours=jet.colors(8),limits=c(-50, 365)) + labs( x= " " ,  y= "", title = "  (b)                         PEAK") + EJapproved + theme(legend.position="none")

Be<-ggplot(gge, aes(x=SOILMOIST, y=maxTWENTYFIVEcm_T, z=DOY))+
  stat_contour(aes(color=..level..),binwidth=10, size=2)+
  scale_color_gradientn(colours=jet.colors(8),limits=c(-50, 365)) + labs( x= " " ,  y= "", colour = "DOY", title = "  (c)                         END") + EJapproved 

#grid.arrange(Bs, Bp, Be, nrow = 1)

#
SMf <- seq(0.05, .12, length.out=100)
TEMPf <- seq(15, 25, length.out=100)
ggfs <- data.frame(expand.grid(SOILMOIST=SMf, TWENTYFIVEcm_T=TEMPf))
ggfp <- data.frame(expand.grid(SOILMOIST=SMf, TWENTYFIVEcm_T=TEMPf))
ggfe <- data.frame(expand.grid(SOILMOIST=SMf, TWENTYFIVEcm_T=TEMPf))

Fs25XSMf<-(lm(from~TWENTYFIVEcm_T*SOILMOIST, data=Fdf))
Fp25XSMf<-(lm(peakDOY~TWENTYFIVEcm_T*SOILMOIST, data=Fdf))
Fe25XSMf<-(lm(to~TWENTYFIVEcm_T*SOILMOIST, data=Fdf))

# prediction from the linear model
ggfs$DOY <-predict(Fs25XSMf,newdata=ggfs)
ggfp$DOY <-predict(Fp25XSMf,newdata=ggfp)
ggfe$DOY <-predict(Fp25XSMf,newdata=ggfe)
# contour plot 

jet.colors <- colorRampPalette(matlab.like(9))
Fs<-ggplot(ggfs, aes(x=SOILMOIST, y=TWENTYFIVEcm_T, z=DOY))+
  stat_contour(aes(color=..level..),binwidth=10, size=2)+
  scale_color_gradientn(colours=jet.colors(8), limits=c(-50, 365)) +
  labs( x= "Ave Soil Moisture " ,  y= "25 cm Ave Temperature", title = "  (d)")+
  EJapproved + theme(legend.position="none") 


Fp<-ggplot(ggfp, aes(x=SOILMOIST, y=TWENTYFIVEcm_T, z=DOY))+
  stat_contour(aes(color=..level..),binwidth=10, size=2)+
  scale_color_gradientn(colours=jet.colors(8), limits=c(-50, 365)) + labs( x= "Ave Soil Moisture " ,  y= "", title = "  (e)") + EJapproved + theme(legend.position="none")

Fe<-ggplot(ggfe, aes(x=SOILMOIST, y=TWENTYFIVEcm_T, z=DOY))+
  stat_contour(aes(color=..level..),binwidth=10, size=2)+
  scale_color_gradientn(colours=jet.colors(8), limits=c(-50, 365)) + labs( x= "Ave Soil Moisture " ,  y= "", title = "  (f)", colour="DOY")  + EJapproved 


tiff("SI File 5.tiff", width = 17, height = 13, units = 'cm', res = 300, compression = 'none')

grid.arrange(Bs,  Bp, Be, Fs, Fp, Fe, nrow = 2) 

dev.off()








##citations:
R.Version()
citation()
citation("lme4")
citation("AICcmodavg")
citation("MuMIn")






