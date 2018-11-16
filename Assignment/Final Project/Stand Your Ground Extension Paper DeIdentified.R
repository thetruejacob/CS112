#######################################################
###            GOV 2001: Extension Paper            ###
###                   May 2, 2018                   ###
#######################################################

###################################
###  Table of Contents:         ###
###     I.   Replication        ###
###       a.    Table 1         ###
###       b.      Figure 1a     ###
###       c.      Figure 1b     ###
###       a.    Table 2         ###
###       b.      Figure 2a     ###
###       c.      Figure 2b     ###
###     II.  Extension          ###
###       a.    Table 1         ###
###       b.      Figure 1a     ###
###       c.      Figure 1b     ###
###       a.    Table 2         ###
###       b.      Figure 2a     ###
###       c.      Figure 2b     ###
###################################

##  SET-UP  ##

# set the working directory

library(foreign) ; library(tsModel) ; library(lmtest) ; library(Epi)
library(splines) ; library(grid); library(vcd)
library (dplyr) ; library (xtable)
library(sandwich) ; library(lmtest)

###########################
###                     ###
###     REPLICATION     ###
###                     ###
###########################

#######################
###                 ###
###     Table 1     ###
###                 ###
#######################

# Load full dataset

alldata <- read.csv("2001_Table1ALL.csv")
alldata$Rate <- alldata$Deaths/alldata$StdPop

####################################################################
##  MEAN MONTHLY COUNTS -- UNADJUSTED AND PER 100,000 POPULATION  ##
####################################################################

### HOMICIDE ###

# Subset Data
allhom <-subset(alldata, Cause=="AllHomicide")
FloridaHomOrig <- subset(allhom, Treatment ==1 & OrigData ==1)
ControlHomOrig <- subset(allhom, Treatment ==0 & OrigData ==1)

# Mean Monthly Counts - Original
mean(FloridaHomOrig$Deaths[FloridaHomOrig$Effective==0]) ##Florida before
mean(FloridaHomOrig$Deaths[FloridaHomOrig$Effective==1]) ##Florida after
mean(ControlHomOrig$Deaths[ControlHomOrig$Effective==0]) ##Control States before
mean(ControlHomOrig$Deaths[ControlHomOrig$Effective==1]) ##Control States after

# Mean Rate Counts 
mean(FloridaHomOrig$Rate[FloridaHomOrig$Effective==0]) #Florida Before
mean(FloridaHomOrig$Rate[FloridaHomOrig$Effective==1]) #Florida After
mean(ControlHomOrig$Rate[ControlHomOrig$Effective==0]) ##Control States before
mean(ControlHomOrig$Rate[ControlHomOrig$Effective==1]) ##Control States after


### FIREARM HOMICIDE ###

# Subset Data
firehom <- subset(alldata, Cause=="FirearmHomicide")
FFireHomOrig <- subset(firehom, Treatment ==1 & OrigData ==1)
CSFireHomOrig <- subset(firehom, Treatment ==0 & OrigData ==1)

# Mean Monthly Counts
mean(FFireHomOrig$Deaths[FFireHomOrig$Effective==0]) ##Florida before
mean(FFireHomOrig$Deaths[FFireHomOrig$Effective==1]) ##Florida after
mean(CSFireHomOrig$Deaths[CSFireHomOrig$Effective==0]) ##Control States before
mean(CSFireHomOrig$Deaths[CSFireHomOrig$Effective==1]) ##Control States after

# Mean Rate Counts
mean(FFireHomOrig$Rate[FFireHomOrig$Effective==0]) #Florida Before
mean(FFireHomOrig$Rate[FFireHomOrig$Effective==1]) # Florida After
mean(CSFireHomOrig$Rate[CSFireHomOrig$Effective==0]) ##Control States before
mean(CSFireHomOrig$Rate[CSFireHomOrig$Effective==1]) ##Control States after


### SUICIDES ###

# Subset Data
allsuicides <- subset(alldata, Cause=="Suicide")
Florida.SuiOrig <- subset(allsuicides, Treatment ==1 & OrigData ==1)
CS.SuiOrig <- subset(allsuicides, Treatment ==0 & OrigData ==1)

# Mean Monthly Counts - All Suicide
mean(Florida.SuiOrig$Deaths[Florida.SuiOrig$Effective==0]) ##Florida before
mean(Florida.SuiOrig$Deaths[Florida.SuiOrig$Effective==1]) ##Florida after
mean(CS.SuiOrig$Deaths[CS.SuiOrig$Effective==0]) ##Control States before
mean(CS.SuiOrig$Deaths[CS.SuiOrig$Effective==1]) ##Control States after

# Mean Rate Counts - All Suicide
mean(Florida.SuiOrig$Rate[Florida.SuiOrig$Effective==0]) #Florida Before
mean(Florida.SuiOrig$Rate[Florida.SuiOrig$Effective==1]) # Florida After
mean(CS.SuiOrig$Rate[CS.SuiOrig$Effective==0]) ##Control States before
mean(CS.SuiOrig$Rate[CS.SuiOrig$Effective==1]) ##Control States after


### FIREARM SUICIDES ###

# Subset Data
FASuic <- subset(alldata, Cause=="FirearmSuicide")
FFASuicOrig <- subset(FASuic, Treatment ==1 & OrigData ==1)
CSFASuicOrig <- subset(FASuic, Treatment ==0 & OrigData ==1)

# Mean Monthly Counts - Firearm Suicide
mean(FFASuicOrig$Deaths[FFASuicOrig$Effective==0]) ##Florida before
mean(FFASuicOrig$Deaths[FFASuicOrig$Effective==1]) ##Florida after
mean(CSFASuicOrig$Deaths[CSFASuicOrig$Effective==0]) ##Control States before
mean(CSFASuicOrig$Deaths[CSFASuicOrig$Effective==1]) ##Control States after

# Mean Rate Counts - Firearm Suicide
mean(FFASuicOrig$Rate[FFASuicOrig$Effective==0]) #Florida Before
mean(FFASuicOrig$Rate[FFASuicOrig$Effective==1]) # Florida After
mean(CSFASuicOrig$Rate[CSFASuicOrig$Effective==0]) ##Control States before
mean(CSFASuicOrig$Rate[CSFASuicOrig$Effective==1]) ##Control States after

########################################
##  TESTS FOR SERIAL AUTOCORRELATION  ##
########################################

# Homicide
BGtesthomF<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                 family=quasipoisson, FloridaHomOrig) 
bgtest(BGtesthomF)
bgtest(BGtesthomF, order=12)
##Result - Florida no
BGtesthomCS<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                  family=quasipoisson, ControlHomOrig) 
bgtest(BGtesthomCS)
bgtest(BGtesthomCS, order=12)
##Result - CS yes


# Homicide by Firearm
BGtestfirehomF<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, FFireHomOrig) 
bgtest(BGtestfirehomF)
BGtestfirehomCS<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                      family=quasipoisson, CSFireHomOrig) 
bgtest(BGtestfirehomCS)
##Result - CS and florida yes


# Suicide
BGtestsuicideF<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, Florida.SuiOrig) 
bgtest(BGtestsuicideF)
bgtest(BGtestsuicideF, order=12)
BGtestsuicideCS<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                      family=quasipoisson, CS.SuiOrig) 
bgtest(BGtestsuicideCS)
bgtest(BGtestsuicideCS, order=12)
##Result - CS and florida yes


## Suicide by Firearm
# Test for autocorrelation for Firearm Suicide
BGtestFsuicideF<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time, 
                      family=quasipoisson, FFASuicOrig) 
bgtest(BGtestFsuicideF, order=12)

BGtestFsuicideCS<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                       family=quasipoisson, CSFASuicOrig) 
bgtest(BGtestFsuicideCS, order=12)
##Result - CS yes but Florida yes


#######################################
##  REPLICATION FOR FLORIDA HOMICIDE ##
#######################################

Rep.Florida.Homicide <- rep(NA)
for (dataset in c("FloridaHomOrig")){
  if (dataset == "FloridaHomOrig"){
    temp.data <- FloridaHomOrig
  }
  Test<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
             family=quasipoisson, temp.data) 
  
  stepchangeestimate <- exp(ci.lin(Test, subset="Effective")[,"Estimate"])
  stepchangebottomCI <- exp(ci.lin(Test, subset="Effective")[,"2.5%"])
  stepchangetopCI <- exp(ci.lin(Test, subset="Effective")[,"97.5%"])
  stepchangep <- summary(Test)$coefficients["Effective",4]
  Step.Changeall <- rbind(stepchangebottomCI, stepchangeestimate, stepchangetopCI, stepchangep)
  trendestimate <- exp(ci.lin(Test, subset="Time")[,"Estimate"])
  trendbottomCI <- exp(ci.lin(Test, subset="Time")[,"2.5%"])
  trendtopCI <- exp(ci.lin(Test, subset="Time")[,"97.5%"])
  trendpval <- summary(Test)$coefficients["Time",4]
  trendall <- rbind(trendbottomCI, trendestimate, trendtopCI, trendpval)
  dta <- data.frame(Step.Changeall, trendall, type = dataset)
  Rep.Florida.Homicide <- cbind(Rep.Florida.Homicide, dta)
  rownames(Rep.Florida.Homicide) = c("Bottom CI","Estimate", "Top CI", "P Value")
}
Rep.Florida.Homicide

## FOR LOOP FOR EVERYTHING WITH ROBUST SE
Table1ModelsRSE <- NA
for (dataset in c("ControlHomOrig", "FFireHomOrig","CSFireHomOrig","Florida.SuiOrig","CS.SuiOrig","FFASuicOrig","CSFASuicOrig")){
  if (dataset == "ControlHomOrig"){
    temp.data <- ControlHomOrig
  }
  else if (dataset == "FFireHomOrig"){
    temp.data <- FFireHomOrig
  }
  else if (dataset == "CSFireHomOrig"){
    temp.data <- CSFireHomOrig
  }
  else if (dataset == "Florida.SuiOrig"){
    temp.data <- Florida.SuiOrig
  }
  else if (dataset == "CS.SuiOrig"){
    temp.data <- CS.SuiOrig
  }
  else if (dataset == "FFASuicOrig"){
    temp.data <- FFASuicOrig
  } 
  else if (dataset == "CSFASuicOrig"){
    temp.data <- CSFASuicOrig
  } 
  Harmonicmodel<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                      family=quasipoisson, temp.data) 
  coef <- coef(Harmonicmodel)["Effective"]
  se <- sqrt(vcovHAC(Harmonicmodel)["Effective","Effective"])
  Step.Change <-c(ll=exp(coef-qnorm(0.975)*se),RR=exp(coef), ul=exp(coef+qnorm(0.975)*se),
                  summary(Harmonicmodel)$coefficients["Effective",4] )
  trendestimate <- exp(ci.lin(Harmonicmodel, subset="Time")[,"Estimate"])
  trendpval <- summary(Harmonicmodel)$coefficients["Time",4]
  trendbottomCI <- exp(ci.lin(Harmonicmodel, subset="Time")[,"2.5%"])
  trendtopCI <- exp(ci.lin(Harmonicmodel, subset="Time")[,"97.5%"])
  trendall <- rbind(trendbottomCI, trendestimate, trendtopCI, trendpval)
  dta3 <- data.frame(Step.Change, trendall, type = dataset)
  Table1ModelsRSE <- cbind(Table1ModelsRSE, dta3)
  rownames(Table1ModelsRSE) = c("Bottom CI","Step Change", "Top CI", "P Value")
}
Table1ModelsRSE


#########################
##  INTERACTION MODELS ##
#########################

## Homicide
HomInteractOrig <- glm(Deaths ~ offset(log(StdPop)) + Effective*Treatment + Time*Treatment
                       + harmonic(MonthNo,2,12),family=quasipoisson, allhom, subset = OrigData ==1)
ci.lin(HomInteractOrig,subset="Effective:Treatment")

## All Suicides
SuicideInteract <- glm(Deaths ~ offset(log(StdPop)) + Effective*Treatment + Time*Treatment
                       + harmonic(MonthNo,2,12),family=quasipoisson, allsuicides, subset = OrigData ==1)
ci.lin(SuicideInteract,subset="Effective:Treatment")

## Firearm Homicide
FireHomInteract <- glm(Deaths ~ offset(log(StdPop)) + Effective*Treatment + Time*Treatment
                       + harmonic(MonthNo,2,12),family=quasipoisson, firehom, subset = OrigData ==1)
ci.lin(FireHomInteract,subset="Effective:Treatment")

## Firearm Suicides
FirearmSuicideInteract <- glm(Deaths ~ offset(log(StdPop)) + Effective*Treatment + Time*Treatment
                              + harmonic(MonthNo,2,12),family=quasipoisson, FASuic, subset = OrigData ==1)
ci.lin(FirearmSuicideInteract,subset="Effective:Treatment")

################
##  Figure 1A ##
################

repallhom <- subset(alldata, Cause=="AllHomicide" & Year < 2015)

##Sept 2001 was removed from the analysis but is added back in to graph
allhomtop<-repallhom[1:224,]
allhomtbottom <- repallhom[225:383,]

##Spet 2001 was given the mean monthly count and rate in the before period 
Sept2001 <- cbind("Sept.,2001", 2001,9,189.4, 0,0,46161275, 461.6128,30,0,1,"AllHomicide","Table1","All",.41)
colnames(Sept2001)<- colnames(repallhom)
allhomtop <- rbind(allhomtop, Sept2001)
repallhom <- rbind(allhomtop,allhomtbottom)
repallhom <- transform(repallhom, StdPop = as.numeric(StdPop), Deaths = as.numeric(Deaths), Year = as.numeric(Year),
                       Treatment = as.numeric(Treatment), Pop = as.numeric(Pop), Time = as.numeric(Time), 
                       Trend = as.numeric(Trend), Rate = as.numeric(Rate), MonthNo = as.numeric(MonthNo))

# Setup variables
repallhom$Level <- c(rep(0,81),rep(1,111),rep(0,81),rep(1,111))
repallhom$Time <- c(rep(c(1:192), 2))
repallhom$Trend <- rep(c(rep(0,81), 1:111),2)
repallhom$TXtime <- repallhom$Treatment * repallhom$Time
repallhom$TXlevel <- repallhom$Treatment * repallhom$Level
repallhom$TXtrend <- repallhom$Treatment * repallhom$Trend
repallhom_crappy <- glm(Rate ~ offset(log(StdPop)) + Time + Level + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data=repallhom)
repallhom.mod <- with(repallhom, Deaths/StdPop)
repallhom.datanew <- data.frame(StdPop=mean(repallhom$StdPop),
                                Level=rep(c(0,1), c(819,1101)),
                                Time = 1:1920/10,
                                MonthNo = rep(1:120/10, 16))

# Models
repallhom.pred1 <- predict(repallhom_crappy,type="response",repallhom.datanew)/mean(repallhom$StdPop)
repallhom.pred2 <- predict(repallhom_crappy,type="response",transform(repallhom,MonthNo=4.8))/mean(repallhom$StdPop)

repallhom_linear <- glm(Rate ~ offset(log(StdPop)) + Time + Treatment + Level +  
                          TXtime + TXlevel, 
                        family=quasipoisson, data=repallhom)
repallhom_harmonic <- glm(Rate ~ offset(log(StdPop)) + Time + Treatment + Level + 
                            TXtime + TXlevel + harmonic(MonthNo,2,12), 
                          family=quasipoisson, data=repallhom)

###Plot for replication all homicide 
plot(repfirehom$Time[1:192],repfirehom.mod[1:192],
     type="n",
     ylim=c(0,1),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:15*12+6,tick=T,labels=1999:2014) 
rect(81,0,192,1, col=grey(0.9),border=F)
points(repallhom$Time[1:192],repallhom.mod[1:192],
       col="darkorange2",
       pch=20)
points(repallhom$Time[193:384],repallhom.mod[193:384],
       col="dodgerblue4",
       pch=20)
lines(repallhom$Time[1:81], fitted(repallhom_linear)[1:81], col="darkorange",lty=5)
lines(repallhom$Time[1:81], fitted(repallhom_harmonic)[1:81], col="darkorange",lwd=2)
lines(repallhom$Time[82:192], fitted(repallhom_linear)[82:192], col="darkorange",lty=5)
lines(repallhom$Time[82:192], fitted(repallhom_harmonic)[82:192], col="darkorange",lwd=2)
lines(repallhom$Time[193:273], fitted(repallhom_linear)[193:273], col="dodgerblue4",lty=5)
lines(repallhom$Time[193:273], fitted(repallhom_harmonic)[193:273], col="dodgerblue4",lwd=2)
lines(repallhom$Time[274:384], fitted(repallhom_linear)[274:384], col="dodgerblue4",lty=5)
lines(repallhom$Time[274:384], fitted(repallhom_harmonic)[274:384], col="dodgerblue4",lwd=2)
legend(x=0, y=0.9, legend=c("Florida","Comparison States"), col=c("darkorange","dodgerblue4"),lwd = 2)
title(main="Replication: Homicide Rates \nFlorida and Comparison States 1999-2014")


################
##  Figure 1B ##
################

repfirehom <- subset(alldata, Cause=="FirearmHomicide" & Year < 2015)

# Setup variables
repfirehom$Level <- c(rep(0,81),rep(1,111),rep(0,81),rep(1,111))
repfirehom$TXtime <- repfirehom$Treatment * repfirehom$Time
repfirehom$TXlevel <- repfirehom$Treatment * repfirehom$Level
repfirehom$TXtrend <- repfirehom$Treatment * repfirehom$Trend

# Models
repfirehom_crappy <- glm(Rate ~ offset(log(StdPop)) + Time + Level + harmonic(MonthNo,2,12), 
                         family=quasipoisson, data=repfirehom)
repfirehom.mod <- with(repfirehom, Deaths/StdPop)
repfirehom.datanew <- data.frame(StdPop=mean(repfirehom$StdPop),
                                 Level=rep(c(0,1), c(819,1101)),
                                 Time = 1:1920/10,
                                 MonthNo = rep(1:120/10, 16))
repfirehom.pred1 <- predict(repfirehom_crappy,type="response",repfirehom.datanew)/mean(repfirehom$StdPop)
repfirehom.pred2 <- predict(repfirehom_crappy,type="response",transform(repfirehom.datanew,MonthNo=4.8))/mean(repfirehom$StdPop)

repfirehom_linear <- glm(Rate ~ offset(log(StdPop)) + Time + Treatment + Level +  
                           TXtime + TXlevel, 
                         family=quasipoisson, data=repfirehom)
repfirehom_harmonic <- glm(Rate ~ offset(log(StdPop)) + Time + Treatment + Level + 
                             TXtime + TXlevel + harmonic(MonthNo,2,12), 
                           family=quasipoisson, data=repfirehom)

# Replication Plot for replication firearm homicide 
plot(repfirehom$Time[1:192],repfirehom.mod[1:192],
     type="n",
     ylim=c(0,1),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:15*12+6,tick=T,labels=1999:2014) 
rect(81,0,192,1, col=grey(0.9),border=F)
points(repfirehom$Time[1:192],repfirehom.mod[1:192],
       col="darkorange2",
       pch=20)
points(repfirehom$Time[193:384],repfirehom.mod[193:384],
       col="dodgerblue4",
       pch=20)
lines(repfirehom$Time[1:81], fitted(repfirehom_linear)[1:81], col="darkorange",lty=5)
lines(repfirehom$Time[1:81], fitted(repfirehom_harmonic)[1:81], col="darkorange",lwd=2)
lines(repfirehom$Time[82:192], fitted(repfirehom_linear)[82:192], col="darkorange",lty=5)
lines(repfirehom$Time[82:192], fitted(repfirehom_harmonic)[82:192], col="darkorange",lwd=2)
lines(repfirehom$Time[193:273], fitted(repfirehom_linear)[193:273], col="dodgerblue4",lty=5)
lines(repfirehom$Time[193:273], fitted(repfirehom_harmonic)[193:273], col="dodgerblue4",lwd=2)
lines(repfirehom$Time[274:384], fitted(repfirehom_linear)[274:384], col="dodgerblue4",lty=5)
lines(repfirehom$Time[274:384], fitted(repfirehom_harmonic)[274:384], col="dodgerblue4",lwd=2)
legend(x=0, y=0.9, legend=c("Florida","Comparison States"), col=c("darkorange","dodgerblue4"),lwd = 2)
title(main="Replication: Firearm Homicide Rates \n Florida and Comparison States 1999-2014")

##########################################################################
###   EXTENSION FIGURE: All Firearm Related Deaths (Combo) 1999-2016   ###
##########################################################################

repcombo <- subset(alldata, Cause=="Combo")

# Setup variables
repcombo$Level <- c(rep(0,81),rep(1,135),rep(0,81),rep(1,135))
repcombo$TXtime <- repcombo$Treatment * repcombo$Time
repcombo$TXlevel <- repcombo$Treatment * repcombo$Level
repcombo$TXtrend <- repcombo$Treatment * repcombo$Trend

# Model
repcombo_linear <- glm(Rate ~ offset(log(StdPop)) + Time + Treatment + Level + Trend + 
                         TXtime + TXlevel + TXtrend, 
                       family=quasipoisson, data=repcombo)
repcombo_harmonic <- glm(Rate ~ offset(log(StdPop)) + Time + Treatment + Level + Trend + 
                           TXtime + TXlevel + TXtrend + harmonic(MonthNo,2,12), 
                         family=quasipoisson, data=repcombo) 

repcombo.mod <- with(repcombo, Deaths/StdPop)
repcombo.datanew <- data.frame(StdPop=repcombo$StdPop,Level=repcombo$Level, Times = repcombo$Time/10, MonthNo = repcombo$MonthNo/10)

repcombo.pred1 <- predict(repcombo_linear,type="response",repcombo)/repcombo$StdPop
repcombo.pred2 <- predict(repcombo_linear,type="response",transform(repcombo,MonthNo=4.8))/(repcombo$StdPop)

###Plot for combo extension 
plot(repcombo$Time[1:216],repcombo.mod[1:216],
     type="n",
     ylim=c(0,1),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:17*12+12,tick=T,labels=1999:2016) 
rect(81,0,216,1, col=grey(0.9),border=F)
points(repcombo$Time[1:216],repcombo.mod[1:216],
       col="darkorange2",
       pch=20)
points(repcombo$Time[217:432],repcombo.mod[217:432],
       col="dodgerblue4",
       pch=20)
lines(repcombo$Time[1:81], fitted(repcombo_linear)[1:81], col="darkorange",lty=5)
lines(repcombo$Time[1:81], fitted(repcombo_harmonic)[1:81], col="darkorange",lwd=2)
lines(repcombo$Time[82:216], fitted(repcombo_linear)[82:216], col="darkorange",lty=5)
lines(repcombo$Time[82:216], fitted(repcombo_harmonic)[82:216], col="darkorange",lwd=2)
lines(repcombo$Time[217:297], fitted(repcombo_linear)[217:297], col="dodgerblue4",lty=5)
lines(repcombo$Time[217:297], fitted(repcombo_harmonic)[217:297], col="dodgerblue4",lwd=2)
lines(repcombo$Time[298:432], fitted(repcombo_linear)[298:432], col="dodgerblue4",lty=5)
lines(repcombo$Time[298:432], fitted(repcombo_harmonic)[298:432], col="dodgerblue4",lwd=2)
legend(x=0, y=0.9, legend=c("Florida","Comparison States"), col=c("darkorange","dodgerblue4"),lwd = 2)
title(main="Extension: All Firearm Related Deaths, \n Florida and Comparison States 1999-2016")




#######################
###                 ###
###     Table 2     ###
###                 ###
#######################

# load data #
repdata <- read.csv("Table 2 Replication Data.csv")

# add columns
repdata$StdPop <- repdata$Population / 100000
repdata$Rate <- repdata$Deaths / repdata$StdPop
repdata$Treatment <- (rep(1, 8856))
repdata$Time <- c(rep(c(1:216), 41))
repdata$After <- rep(c(rep(0,81), rep(1,135)),41)
repdata$Trend <- rep(c(rep(0,81), 1:135),41)

##  SUBSET THE DATA  ##

# Homicide 
hom.white <- subset(repdata, Characteristic == "White" & Cause == "Homicide" & Year <2015)
hom.black <- subset(repdata, Characteristic == "Black or African American" & Cause == "Homicide" & Year <2015)
hom.male <- subset(repdata, Characteristic == "M" & Cause == "Homicide" & Year <2015)
hom.female <- subset(repdata, Characteristic == "F" & Cause == "Homicide" & Year <2015)
hom.20_34 <- subset(repdata, Characteristic == "20-34" & Cause == "Homicide" & Year <2015)
hom.35over <- subset(repdata, Characteristic == "35 and over" & Cause == "Homicide" & Year <2015)

# Suicide 
suic.white <- subset(repdata, Characteristic == "White" & Cause == "Suicide" & Year <2015)
suic.black <- subset(repdata, Characteristic == "Black or African American" & Cause == "Suicide" & Year <2015)
suic.male <- subset(repdata, Characteristic == "M" & Cause == "Suicide" & Year <2015)
suic.female <- subset(repdata, Characteristic == "F" & Cause == "Suicide" & Year <2015)
suic.20_34 <- subset(repdata, Characteristic == "20-34" & Cause == "Suicide" & Year <2015)
suic.35over <- subset(repdata, Characteristic == "35 and over" & Cause == "Suicide" & Year <2015)

# Firearm Homicide
FAhom.white <- subset(repdata, Characteristic == "White" & Cause == "Firearm Homicide" & Year <2015)
FAhom.black <- subset(repdata, Characteristic == "Black or African American" & Cause == "Firearm Homicide" & Year <2015)
FAhom.male <- subset(repdata, Characteristic == "M" & Cause == "Firearm Homicide" & Year <2015)
FAhom.20_34 <- subset(repdata, Characteristic == "20-34" & Cause == "Firearm Homicide" & Year <2015)
FAhom.35over <- subset(repdata, Characteristic == "35 and over" & Cause == "Firearm Homicide" & Year <2015)

# Firearm Suicide
FAsuic.white <- subset(repdata, Characteristic == "White" & Cause == "Firearm Suicide" & Year <2015)
FAsuic.black <- subset(repdata, Characteristic == "Black or African American" & Cause == "Firearm Suicide" & Year <2015)
FAsuic.male <- subset(repdata, Characteristic == "M" & Cause == "Firearm Suicide" & Year <2015)
FAsuic.female <- subset(repdata, Characteristic == "F" & Cause == "Firearm Suicide" & Year <2015)
FAsuic.20_34 <- subset(repdata, Characteristic == "20-34" & Cause == "Firearm Suicide" & Year <2015)
FAsuic.35over <- subset(repdata, Characteristic == "35 and over" & Cause == "Firearm Suicide" & Year <2015)

###########################
##  MEAN MONTHLY COUNTS  ##
###########################

# Homicide 

summary(hom.white$Deaths[repdata$After==0])
summary(hom.white$Deaths[repdata$After==1])

summary(hom.black$Deaths[repdata$After==0])
summary(hom.black$Deaths[repdata$After==1])

summary(hom.20_34$Deaths[repdata$After==0])
summary(hom.20_34$Deaths[repdata$After==1])

summary(hom.35over$Deaths[repdata$After==0])
summary(hom.35over$Deaths[repdata$After==1])

summary(hom.male$Deaths[repdata$After==0])
summary(hom.male$Deaths[repdata$After==1])

summary(hom.female$Deaths[repdata$After==0])
summary(hom.female$Deaths[repdata$After==1])

# Suicide

summary(suic.white$Deaths[repdata$After==0])
summary(suic.white$Deaths[repdata$After==1])

summary(suic.20_34$Deaths[repdata$After==0])
summary(suic.20_34$Deaths[repdata$After==1])

summary(suic.35over$Deaths[repdata$After==0])
summary(suic.35over$Deaths[repdata$After==1])

summary(suic.male$Deaths[repdata$After==0])
summary(suic.male$Deaths[repdata$After==1])

summary(suic.female$Deaths[repdata$After==0])
summary(suic.female$Deaths[repdata$After==1])

# Firearm Homicide

summary(FAhom.white$Deaths[repdata$After==0])
summary(FAhom.white$Deaths[repdata$After==1])

summary(FAhom.black$Deaths[repdata$After==0])
summary(FAhom.black$Deaths[repdata$After==1])

summary(FAhom.20_34$Deaths[repdata$After==0])
summary(FAhom.20_34$Deaths[repdata$After==1])

summary(FAhom.35over$Deaths[repdata$After==0])
summary(FAhom.35over$Deaths[repdata$After==1])

summary(FAhom.male$Deaths[repdata$After==0])
summary(FAhom.male$Deaths[repdata$After==1])

# Firearm Suicide

summary(FAsuic.white$Deaths[repdata$After==0])
summary(FAsuic.white$Deaths[repdata$After==1])

summary(FAsuic.20_34$Deaths[repdata$After==0])
summary(FAsuic.20_34$Deaths[repdata$After==1])

summary(FAsuic.35over$Deaths[repdata$After==0])
summary(FAsuic.35over$Deaths[repdata$After==1])

summary(FAsuic.male$Deaths[repdata$After==0])
summary(FAsuic.male$Deaths[repdata$After==1])

summary(FAsuic.female$Deaths[repdata$After==0])
summary(FAsuic.female$Deaths[repdata$After==1])


##################################################
##  MEAN MONTHLY COUNTS PER 100,000 POPULATION  ##
##################################################

# Homicide

summary(hom.white$Rate[repdata$After==0])
summary(hom.white$Rate[repdata$After==1])

summary(hom.black$Rate[repdata$After==0])
summary(hom.black$Rate[repdata$After==1])

summary(hom.20_34$Rate[repdata$After==0])
summary(hom.20_34$Rate[repdata$After==1])

summary(hom.35over$Rate[repdata$After==0])
summary(hom.35over$Rate[repdata$After==1])

summary(hom.male$Rate[repdata$After==0])
summary(hom.male$Rate[repdata$After==1])

summary(hom.female$Rate[repdata$After==0])
summary(hom.female$Rate[repdata$After==1])

# Suicide

summary(suic.white$Rate[repdata$After==0])
summary(suic.white$Rate[repdata$After==1])

summary(suic.20_34$Rate[repdata$After==0])
summary(suic.20_34$Rate[repdata$After==1])

summary(suic.35over$Rate[repdata$After==0])
summary(suic.35over$Rate[repdata$After==1])

summary(suic.male$Rate[repdata$After==0])
summary(suic.male$Rate[repdata$After==1])

summary(suic.female$Rate[repdata$After==0])
summary(suic.female$Rate[repdata$After==1])

# Firearm Homicide

summary(FAhom.white$Rate[repdata$After==0])
summary(FAhom.white$Rate[repdata$After==1])

summary(FAhom.black$Rate[repdata$After==0])
summary(FAhom.black$Rate[repdata$After==1])

summary(FAhom.20_34$Rate[repdata$After==0])
summary(FAhom.20_34$Rate[repdata$After==1])

summary(FAhom.35over$Rate[repdata$After==0])
summary(FAhom.35over$Rate[repdata$After==1])

summary(FAhom.male$Rate[repdata$After==0])
summary(FAhom.male$Rate[repdata$After==1])

# Firearm Suicide

summary(FAsuic.white$Rate[repdata$After==0])
summary(FAsuic.white$Rate[repdata$After==1])

summary(FAsuic.20_34$Rate[repdata$After==0])
summary(FAsuic.20_34$Rate[repdata$After==1])

summary(FAsuic.35over$Rate[repdata$After==0])
summary(FAsuic.35over$Rate[repdata$After==1])

summary(FAsuic.male$Rate[repdata$After==0])
summary(FAsuic.male$Rate[repdata$After==1])

summary(FAsuic.female$Rate[repdata$After==0])
summary(FAsuic.female$Rate[repdata$After==1])


###############################################
##  Step Change & Tests for Autocorrelation  ##
###############################################

###  HOMICIDE  ###

# WHITE #

## Non-Seasonal Model
m.hom.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                          family=quasipoisson, data = hom.white)
# Residual plot
hom.white.res1 <- residuals(m.hom.white,type="deviance")

plot(hom.white$Time,hom.white.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="White Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.hom.white)
bgtest(m.hom.white, order=12)

## Seasonal model
sm.hom.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                   family=quasipoisson, data = hom.white)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.hom.white)
bgtest(sm.hom.white, order=12)


# BLACK #

## Non-Seasonal Model
m.hom.black <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = hom.black)
# Residual plot
hom.black.res1 <- residuals(m.hom.black,type="deviance")

plot(hom.black$Time,hom.black.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="White Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.hom.black)
bgtest(m.hom.black, order=12)

## Seasonal model
sm.hom.black <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = hom.black)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.hom.black)
bgtest(sm.hom.black, order=12)


# YOUNGER #

## Non-Seasonal Model
m.hom.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = hom.20_34)
# Residual plot
hom.20_34.res1 <- residuals(m.hom.20_34,type="deviance")

plot(hom.20_34$Time,hom.20_34.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Younger Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.hom.20_34)            #p = 0.068
bgtest(m.hom.20_34, order=12)

## Seasonal model
sm.hom.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = hom.20_34)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.hom.20_34)
bgtest(sm.hom.20_34, order=12)


# OLDER #

## Non-Seasonal Model
m.hom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = hom.35over)
# Residual plot
hom.35over.res1 <- residuals(m.hom.35over,type="deviance")

plot(hom.35over$Time,hom.35over.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Older Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.hom.35over)  
bgtest(m.hom.35over, order=12)

## Seasonal model
sm.hom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = hom.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.hom.35over)
bgtest(sm.hom.35over, order=12)


# MALE #

## Non-Seasonal Model
m.hom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                    family=quasipoisson, data = hom.male)
# Residual plot
hom.male.res1 <- residuals(m.hom.male,type="deviance")

plot(hom.male$Time,hom.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.hom.male)  
bgtest(m.hom.male, order=12)

## Seasonal model
sm.hom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = hom.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.hom.male)
bgtest(sm.hom.male, order=12)
       

# FEMALE #

## Non-Seasonal Model
m.hom.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                  family=quasipoisson, data = hom.female)
# Residual plot
hom.female.res1 <- residuals(m.hom.female,type="deviance")

plot(hom.female$Time,hom.female.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Female Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.hom.female)  
bgtest(m.hom.female, order=12)

## Seasonal model
sm.hom.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                   family=quasipoisson, data = hom.female)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.hom.female)
bgtest(sm.hom.female, order=12)



###  SUICIDE  ###

# WHITE #

## Non-Seasonal Model
m.suic.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = suic.white)
# Residual plot
suic.white.res1 <- residuals(m.suic.white,type="deviance")

plot(suic.white$Time,suic.white.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="White Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.suic.white)            # p < 0.05
bgtest(m.suic.white, order=12)  # p < 0.01

## Seasonal model
sm.suic.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = suic.white)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.suic.white)
bgtest(sm.suic.white, order=12)  # p < 0.05


# YOUNGER #

## Non-Seasonal Model
m.suic.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = suic.20_34)
# Residual plot
suic.20_34.res1 <- residuals(m.suic.20_34,type="deviance")

plot(suic.20_34$Time,suic.20_34.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Younger Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.suic.20_34)  
bgtest(m.suic.20_34, order=12)

## Seasonal model
sm.suic.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = suic.20_34)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.suic.20_34)
bgtest(sm.suic.20_34, order=12)


# OLDER #

## Non-Seasonal Model
m.suic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                    family=quasipoisson, data = suic.35over)
# Residual plot
suic.35over.res1 <- residuals(m.suic.35over,type="deviance")

plot(suic.35over$Time,suic.35over.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Older Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.suic.35over)            # p < 0.05
bgtest(m.suic.35over, order=12)  # p < 0.01

## Seasonal model
sm.suic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = suic.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.suic.35over)
bgtest(sm.suic.35over, order=12)  # p = 0.06


# MALE #

## Non-Seasonal Model
m.suic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                  family=quasipoisson, data = suic.male)
# Residual plot
suic.male.res1 <- residuals(m.suic.male,type="deviance")

plot(suic.male$Time,suic.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.suic.male)            # p < 0.05  
bgtest(m.suic.male, order=12)  # p < 0.01

## Seasonal model
sm.suic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                   family=quasipoisson, data = suic.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.suic.male)
bgtest(sm.suic.male, order=12)  # p < 0.05


# FEMALE #

## Non-Seasonal Model
m.suic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                    family=quasipoisson, data = suic.female)
# Residual plot
suic.female.res1 <- residuals(m.suic.female,type="deviance")

plot(suic.female$Time,suic.female.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Female Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.suic.female)            # p < 0.01
bgtest(m.suic.female, order=12)  # p < 0.05

## Seasonal model
sm.suic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = suic.female)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.suic.female)            # p < 0.05
bgtest(sm.suic.female, order=12)



###  FIREARM HOMICIDE  ###

# WHITE #

## Non-Seasonal Model
m.FAhom.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = FAhom.white)
# Residual plot
FAhom.white.res1 <- residuals(m.FAhom.white,type="deviance")

plot(FAhom.white$Time,FAhom.white.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="White Firearm Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAhom.white)
bgtest(m.FAhom.white, order=12)

## Seasonal model
sm.FAhom.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = FAhom.white)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAhom.white)
bgtest(sm.FAhom.white, order=12)


# BLACK #

## Non-Seasonal Model
m.FAhom.black <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = FAhom.black)

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAhom.black)             # p < 0.05
bgtest(m.FAhom.black, order=12)

## Seasonal model
sm.FAhom.black <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = FAhom.black)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAhom.black)
bgtest(sm.FAhom.black, order=12)


# YOUNGER #

## Non-Seasonal Model
m.FAhom.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = FAhom.20_34)

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAhom.20_34)             # p < 0.05
bgtest(m.FAhom.20_34, order=12)

## Seasonal model
sm.FAhom.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = FAhom.20_34)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAhom.20_34)            # p < 0.05
bgtest(sm.FAhom.20_34, order=12)


# OLDER #

## Non-Seasonal Model
m.FAhom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                    family=quasipoisson, data = FAhom.35over)

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAhom.35over)            # p < 0.001
bgtest(m.FAhom.35over, order=12)

## Seasonal model
sm.FAhom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = FAhom.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAhom.35over)           # p < 0.001
bgtest(sm.FAhom.35over, order=12)


# MALE #

## Non-Seasonal Model
m.FAhom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                  family=quasipoisson, data = FAhom.male)
# Residual plot
FAhom.male.res1 <- residuals(m.FAhom.male,type="deviance")

plot(FAhom.male$Time,FAhom.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Firearm Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAhom.male)            # p < 0.01  
bgtest(m.FAhom.male, order=12)  # p ~ 0.05

## Seasonal model
sm.FAhom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                   family=quasipoisson, data = FAhom.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAhom.male)            # p < 0.01
bgtest(sm.FAhom.male, order=12)  # p < 0.05



###  FIREARM SUICIDE  ###

# WHITE #

## Non-Seasonal Model
m.FAsuic.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                    family=quasipoisson, data = FAsuic.white)
# Residual plot
FAsuic.white.res1 <- residuals(m.FAsuic.white,type="deviance")

plot(FAsuic.white$Time,FAsuic.white.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="White Firearm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAsuic.white)            
bgtest(m.FAsuic.white, order=12)  

## Seasonal model
sm.FAsuic.white <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = FAsuic.white)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAsuic.white)
bgtest(sm.FAsuic.white, order=12)


# YOUNGER #

## Non-Seasonal Model
m.FAsuic.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                    family=quasipoisson, data = FAsuic.20_34)

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAsuic.20_34)
bgtest(m.FAsuic.20_34, order=12)

## Seasonal model
sm.FAsuic.20_34 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = FAsuic.20_34)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAsuic.20_34)
bgtest(sm.FAsuic.20_34, order=12)


# OLDER #

## Non-Seasonal Model
m.FAsuic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                     family=quasipoisson, data = FAsuic.35over)
# Residual plot
FAsuic.35over.res1 <- residuals(m.FAsuic.35over,type="deviance")

plot(FAsuic.35over$Time,FAsuic.35over.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Older Firearm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAsuic.35over)           
bgtest(m.FAsuic.35over, order=12)  # p < 0.05

## Seasonal model
sm.FAsuic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                      family=quasipoisson, data = FAsuic.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAsuic.35over)
bgtest(sm.FAsuic.35over, order=12)  # p = 0.088


# MALE #

## Non-Seasonal Model
m.FAsuic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                   family=quasipoisson, data = FAsuic.male)
# Residual plot
FAsuic.male.res1 <- residuals(m.FAsuic.male,type="deviance")

plot(FAsuic.male$Time,FAsuic.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Firearm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAsuic.male)            
bgtest(m.FAsuic.male, order=12)  # p = 0.07

## Seasonal model
sm.FAsuic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = FAsuic.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAsuic.male)
bgtest(sm.FAsuic.male, order=12) 


# FEMALE #

## Non-Seasonal Model
m.FAsuic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time, 
                     family=quasipoisson, data = FAsuic.female)

# Breusch-Godfrey test for autocorrelation
bgtest(m.FAsuic.female)            # p < 0.005
bgtest(m.FAsuic.female, order=12)  # p < 0.05

## Seasonal model
sm.FAsuic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + harmonic(MonthNo,2,12), 
                      family=quasipoisson, data = FAsuic.female)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.FAsuic.female)            # p < 0.05
bgtest(sm.FAsuic.female, order=12)


############################
##  Relative Risk 95% CI  ##
############################

###  HOMICIDE  ###
# Breusch-Godrey tests did not indicate statistically significant autocorrelation
# so, we calculate and report normal standard errors for seasonal models of homicide

round(ci.lin(sm.hom.white,Exp=T)["After",5:7],2)
summary(sm.hom.white,Exp=T)$coefficients["After",4]

round(ci.lin(sm.hom.black,Exp=T)["After",5:7],2)
summary(sm.hom.black,Exp=T)$coefficients["After",4]

round(ci.lin(sm.hom.20_34,Exp=T)["After",5:7],2)
summary(sm.hom.20_34,Exp=T)$coefficients["After",4]

round(ci.lin(sm.hom.35over,Exp=T)["After",5:7],2)
summary(sm.hom.35over,Exp=T)$coefficients["After",4]

round(ci.lin(sm.hom.male,Exp=T)["After",5:7],2)
summary(sm.hom.male,Exp=T)$coefficients["After",4]

round(ci.lin(sm.hom.female,Exp=T)["After",5:7],2)
summary(sm.hom.female,Exp=T)$coefficients["After",4]


###  SUICIDE  ###
# Breusch-Godrey tests indicated statistically significant autocorrelation
# for white, older, male and female

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.suic.white )["After"]
se <- sqrt(vcovHAC(sm.suic.white )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.suic.white,Exp=T)$coefficients["After",4]

round(ci.lin(m.suic.20_34,Exp=T)["After",5:7],2)
summary(m.suic.20_34,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.suic.35over )["After"]
se <- sqrt(vcovHAC(sm.suic.35over )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.suic.35over,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.suic.male )["After"]
se <- sqrt(vcovHAC(sm.suic.male )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.suic.male,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.suic.female )["After"]
se <- sqrt(vcovHAC(sm.suic.female )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.suic.female,Exp=T)$coefficients["After",4]


###  FIREARM HOMICIDE  ###
# Breusch-Godrey tests indicated statistically significant autocorrelation
# for younger, older and male

round(ci.lin(sm.FAhom.white,Exp=T)["After",5:7],2)
summary(sm.FAhom.white,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.FAhom.black )["After"]
se <- sqrt(vcovHAC(sm.FAhom.black )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),3)
summary(sm.FAhom.black,Exp=T)$coefficients["After",4]

round(ci.lin(m.FAhom.20_34,Exp=T)["After",5:7],2)  # this gets closest to the numbers in the paper
summary(m.FAhom.20_34,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.FAhom.35over )["After"]
se <- sqrt(vcovHAC(sm.FAhom.35over )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.FAhom.35over,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.FAhom.male )["After"]
se <- sqrt(vcovHAC(sm.FAhom.male )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.FAhom.male,Exp=T)$coefficients["After",4]


###  FIREARM SUICIDE  ###
# Breusch-Godrey tests indicated statistically significant autocorrelation
# for older and female, however, the confidence intervals reported in the 
# paper appear to use normal standard errors

round(ci.lin(m.FAsuic.white,Exp=T)["After",5:7],2)
summary(m.FAsuic.white,Exp=T)$coefficients["After",4]

round(ci.lin(m.FAsuic.20_34,Exp=T)["After",5:7],2)
summary(m.FAsuic.20_34,Exp=T)$coefficients["After",4]

round(ci.lin(sm.FAsuic.35over,Exp=T)["After",5:7],2) #this produces the exact numbers from the paper
summary(sm.FAsuic.35over,Exp=T)$coefficients["After",4]

round(ci.lin(sm.FAsuic.male,Exp=T)["After",5:7],2)
summary(sm.FAsuic.male,Exp=T)$coefficients["After",4]

round(ci.lin(sm.FAsuic.female,Exp=T)["After",5:7],2) #this produces the exact numbers from the paper
summary(sm.FAsuic.female,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.FAsuic.35over )["After"]
se <- sqrt(vcovHAC(sm.FAsuic.35over )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.FAsuic.35over,Exp=T)$coefficients["After",4]

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.FAsuic.female )["After"]
se <- sqrt(vcovHAC(sm.FAsuic.female )["After","After"])
round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
summary(sm.FAsuic.female,Exp=T)$coefficients["After",4]


###############
##  FIGURES  ##
###############

###  Figure 2.A ###

# Data Set-Up
young <- subset(repdata, Characteristic == "20-34" & Cause == "Homicide" & Year < 2015| 
                  Characteristic == "20-34" & Cause == "Firearm Homicide" & Year < 2015)

# Setup variables
young$StdPop <- young$Population / 100000
young$Rate <- young$Deaths / young$StdPop
young$time <- rep(c(1:192), 2)
young$level <- c(rep(0,81),rep(1,111),rep(0,81),rep(1,111))
young$trend <- c(rep(0,81),1:111,rep(0,81),1:111)
young$firearms <- c(rep(1, 192), rep(0, 192))
young$FAtime <- young$firearms * young$time
young$FAlevel <- young$firearms * young$level
young$FAtrend <- young$firearms * young$trend

young_crappy <- glm(Rate ~ offset(log(StdPop)) + level + time + 
                      harmonic(MonthNo,2,12), family=quasipoisson, data=young)
young.mod <- with(young, Deaths/StdPop)
young.datanew <- data.frame(StdPop=mean(young$StdPop),level=rep(c(0,1),c(819,1101)),
                            time= 1:1920/10, MonthNo=rep(1:120/10,16))
young.crappy.pred1 <- predict(young_crappy,type="response",young.datanew)/mean(young$StdPop)
young.crappy.pred2 <- predict(young_crappy,type="response",transform(young.datanew,MonthNo=4.8))/mean(young$StdPop)

young_linear <- glm(Rate ~ offset(log(StdPop)) + time + firearms + FAtime + level + trend + FAlevel + FAtrend, 
                    family=quasipoisson, data=young)
young_harmonic <- glm(Rate ~ offset(log(StdPop)) + time + firearms + FAtime + level + trend + FAlevel + FAtrend + 
                        harmonic(MonthNo,2,12), family=quasipoisson, data=young)

pdf("RFig2A_HomvFAHomYoung.pdf")
plot(young$time[1:192],young.mod[1:192],
     type="n",
     ylim=c(0,2.5),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:15*12+6,tick=T,labels=1999:2014) 
rect(81,0,192,2.5, col=grey(0.9),border=F)
points(young$time[1:192],young.mod[1:192],
       col="dodgerblue4",
       pch=20)
points(young$time[193:384],young.mod[193:384],
       col="darkorange",
       pch=20)
lines(young$time[1:81], fitted(young_linear)[1:81], col="dodgerblue4",lty=5)
lines(young$time[1:81], fitted(young_harmonic)[1:81], col="dodgerblue4",lwd=2)
lines(young$time[82:192], fitted(young_linear)[82:192], col="dodgerblue4",lty=5)
lines(young$time[82:192], fitted(young_harmonic)[82:192], col="dodgerblue4",lwd=2)
lines(young$time[193:273], fitted(young_linear)[193:273], col="darkorange",lty=5)
lines(young$time[193:273], fitted(young_harmonic)[193:273], col="darkorange",lwd=2)
lines(young$time[274:384], fitted(young_linear)[274:384], col="darkorange",lty=5)
lines(young$time[274:384], fitted(young_harmonic)[274:384], col="darkorange",lwd=2)
legend(x=0, y=2.4, legend=c("Homicide","Homicide by firearm"), col=c("darkorange","dodgerblue4"),lwd = 2)
title(main="Change in homicide and all firearm deaths \n among people under age 35, 1999-2014")
graphics.off()

###  Figure 2.A ###

# Data Set-Up
white <- subset(repdata, Characteristic == "White" & Cause == "Homicide" & Year < 2015| 
                  Characteristic == "White" & Cause == "Firearm Homicide" & Year < 2015)

# Setup variables
white$StdPop <- white$Population / 100000
white$Rate <- white$Deaths / white$StdPop
white$time <- rep(c(1:192), 2)
white$level <- c(rep(0,81),rep(1,111),rep(0,81),rep(1,111))
white$trend <- c(rep(0,81),1:111,rep(0,81),1:111)
white$firearms <- c(rep(1, 192), rep(0, 192))
white$FAtime <- white$firearms * white$time
white$FAlevel <- white$firearms * white$level
white$FAtrend <- white$firearms * white$trend

white_crappy <- glm(Rate ~ offset(log(StdPop)) + level + time + 
                      harmonic(MonthNo,2,12), family=quasipoisson, data=white)
white.mod <- with(white, Deaths/StdPop)
white.datanew <- data.frame(StdPop=mean(white$StdPop),level=rep(c(0,1),c(819,1101)),
                            time= 1:1920/10, MonthNo=rep(1:120/10,16))
white.crappy.pred1 <- predict(young_crappy,type="response",white.datanew)/mean(white$StdPop)
white.crappy.pred2 <- predict(young_crappy,type="response",transform(white.datanew,MonthNo=4.8))/mean(white$StdPop)

white_linear <- glm(Rate ~ offset(log(StdPop)) + time + firearms + FAtime + level + trend + FAlevel + FAtrend, 
                    family=quasipoisson, data=white)
white_harmonic <- glm(Rate ~ offset(log(StdPop)) + time + firearms + FAtime + level + trend + FAlevel + FAtrend + 
                        harmonic(MonthNo,2,12), family=quasipoisson, data=white)

pdf("RFig2B_HomvFAHomWhite.pdf")
plot(white$time[1:192],white.mod[1:192],
     type="n",
     ylim=c(0,2.5),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:15*12+6,tick=T,labels=1999:2014) 
rect(81,0,192,2.5, col=grey(0.9),border=F)
points(white$time[1:192],white.mod[1:192],
       col="dodgerblue4",
       pch=20)
points(white$time[193:384],white.mod[193:384],
       col="darkorange",
       pch=20)
lines(white$time[1:81], fitted(white_linear)[1:81], col="dodgerblue4",lty=5)
lines(white$time[1:81], fitted(white_harmonic)[1:81], col="dodgerblue4",lwd=2)
lines(white$time[82:192], fitted(white_linear)[82:192], col="dodgerblue4",lty=5)
lines(white$time[82:192], fitted(white_harmonic)[82:192], col="dodgerblue4",lwd=2)
lines(white$time[193:273], fitted(white_linear)[193:273], col="darkorange",lty=5)
lines(white$time[193:273], fitted(white_harmonic)[193:273], col="darkorange",lwd=2)
lines(white$time[274:384], fitted(white_linear)[274:384], col="darkorange",lty=5)
lines(white$time[274:384], fitted(white_harmonic)[274:384], col="darkorange",lwd=2)
legend(x=0, y=2.4, legend=c("Homicide","Homicide by firearm"), col=c("darkorange","dodgerblue4"),lwd = 2)
title(main="Change in homicide and all firearm deaths among whites, 1999-2014")
graphics.off()


#########################
###                   ###
###     EXTENSION     ###
###                   ###
#########################

#######################
###                 ###
###     Table 1     ###
###                 ###
#######################

# Subsets
FallHom2016 <- subset(allhom, Treatment ==1)
CSallHom2016 <- subset(allhom, Treatment ==0) 
FFireHom2016 <- subset(firehom, Treatment ==1)
CSFireHom2016 <- subset(firehom, Treatment ==0) 
Florida.Sui2016 <- subset(allsuicides, Treatment ==1)
CS.Sui2016 <- subset(allsuicides, Treatment ==0)
FFASuic2016 <- subset(FASuic, Treatment ==1)
CSFASuic2016 <- subset(FASuic, Treatment ==0)

# NEW DATA - (Combo) ALL Firearm Related Deaths 
firecombo <- subset(alldata, Cause=="Combo")
FCombo2016 <- subset(firecombo, Treatment ==1)
CSCombo2016 <- subset(firecombo, Treatment ==0)

############################################
###     MEAN AND RATE MONTHLY COUNTS     ###
############################################

##2016 Homicide 
mean(FallHom2016$Deaths[FallHom2016$Effective==1]) ##Florida after
mean(CSallHom2016$Deaths[CSallHom2016$Effective==1]) ##Control States after
mean(FallHom2016$Rate[FallHom2016$Effective==1]) #Florida After
mean(CSallHom2016$Rate[CSallHom2016$Effective==1]) ##Control States after

##2016 - Firearm Homicide 
mean(FFireHom2016$Deaths[FFireHom2016$Effective==1]) ##Florida after
mean(CSFireHom2016$Deaths[CSFireHom2016$Effective==1]) ##Control States after
mean(FFireHom2016$Rate[FFireHom2016$Effective==1]) # Florida After
mean(CSFireHom2016$Rate[CSFireHom2016$Effective==1]) ##Control States after

##2016 - All Suicide
mean(Florida.Sui2016$Deaths[Florida.Sui2016$Effective==1]) ##Florida 2016
mean(CS.Sui2016$Deaths[CS.Sui2016$Effective==1]) ##Control States 2016
mean(Florida.Sui2016$Rate[Florida.Sui2016$Effective==1]) # Florida 2016
mean(CS.Sui2016$Rate[CS.Sui2016$Effective==1]) ##Control States 2016

##2016 - Firearm Suicide
mean(FFASuic2016$Deaths[FFASuic2016$Effective==1]) ##Florida after
mean(CSFASuic2016$Deaths[CSFASuic2016$Effective==1]) ##Control States after
mean(FFASuic2016$Rate[FFASuic2016$Effective==1]) # Florida After
mean(CSFASuic2016$Rate[CSFASuic2016$Effective==1]) ##Control States after

##ALL FIREARM RELATED DEATHS (Combo)
#Mean Monthly Counts - Combo
mean(FCombo2016$Deaths[FCombo2016$Effective==0]) ##Florida before
mean(FCombo2016$Deaths[FCombo2016$Effective==1]) ##Florida after
mean(CSCombo2016$Deaths[CSCombo2016$Effective==0]) ##Control States before
mean(CSCombo2016$Deaths[CSCombo2016$Effective==1]) ##Control States after

##Mean Rate Counts - Combo
mean(FCombo2016$Rate[FCombo2016$Effective==0]) #Florida Before
mean(FCombo2016$Rate[FCombo2016$Effective==1]) # Florida After
mean(CSCombo2016$Rate[CSCombo2016$Effective==0]) ##Control States before
mean(CSCombo2016$Rate[CSCombo2016$Effective==1]) ##Control States after


################################################
###     TESTS FOR SERIAL AUTOCORRELATION     ###
################################################

# Homicide
BGtesthomF2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                     family=quasipoisson, FallHom2016) 
bgtest(BGtesthomF2016, order=12)
BGtesthomCS2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                      family=quasipoisson, CSallHom2016) 
bgtest(BGtesthomCS2016, order=12)
##Result - yes both


# Homicide by Firearm
BGtestfirehomF2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                         family=quasipoisson, FFireHom2016) 
bgtest(BGtestfirehomF2016)
BGtestfirehomCS2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                          family=quasipoisson, CSFireHom2016) 
bgtest(BGtestfirehomCS2016)
##Result - yes both


# Suicide
BGtestsuicideF2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                         family=quasipoisson, Florida.Sui2016) 
bgtest(BGtestsuicideF2016, order=12)
BGtestsuicideCS2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                          family=quasipoisson, CS.Sui2016) 
bgtest(BGtestsuicideCS2016, order=12)
##Result - yes both


# Suicide by Firearm
BGtestFsuicideF20161<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time, 
                           family=quasipoisson, FFASuic2016) 
bgtest(BGtestFsuicideF20161)
bgtest(BGtestFsuicideF20161, order=12)
BGtestFsuicideF20162<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                           family=quasipoisson, FFASuic2016) 
bgtest(BGtestFsuicideF20162)
bgtest(BGtestFsuicideF20162, order=12)
BGtestFsuicideCS2016<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + harmonic(MonthNo,2,12), 
                           family=quasipoisson, CSFASuic2016) 
bgtest(BGtestFsuicideCS2016, order=12)
##Result - yes to Controls, yes to Florida 

# All firearm related deaths 
BGtestFcombo<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time, 
                   family=quasipoisson, FCombo2016) 
bgtest(BGtestFcombo)
BGtestCScombo<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time, 
                    family=quasipoisson, CSCombo2016) 
bgtest(BGtestCScombo)
#yes to serial autocorrelation 


###############################
##  EXTENSION WITH ROBUST SE ##
###############################

Table2ModelsRSE <- NA
for (dataset in c("FallHom2016", "CSallHom2016","FFireHom2016","CSFireHom2016","Florida.Sui2016","CS.Sui2016",
                  "FFASuic2016","CSFASuic2016","FCombo2016","CSCombo2016")){
  if (dataset == "FallHom2016"){
    temp.data <- FallHom2016
  }
  else if (dataset == "CSallHom2016"){
    temp.data <- CSallHom2016
  }
  else if (dataset == "FFireHom2016"){
    temp.data <- FFireHom2016
  }
  else if (dataset == "CSFireHom2016"){
    temp.data <- CSFireHom2016
  }
  else if (dataset == "Florida.Sui2016"){
    temp.data <- Florida.Sui2016
  }
  else if (dataset == "CS.Sui2016"){
    temp.data <- CS.Sui2016
  } 
  else if (dataset == "FFASuic2016"){
    temp.data <- FFASuic2016
  } 
  else if (dataset == "CSFASuic2016"){
    temp.data <- CSFASuic2016
  }
  else if (dataset == "FCombo2016"){
    temp.data <- FCombo2016
  }
  else if (dataset == "CSCombo2016"){
    temp.data <- CSCombo2016
  }
  Harmonicmodel<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Trend + harmonic(MonthNo,2,12), 
                      family=quasipoisson, temp.data) 
  coef <- coef(Harmonicmodel)["Effective"]
  se <- sqrt(vcovHAC(Harmonicmodel)["Effective","Effective"])
  Step.Change <-c(exp(coef-qnorm(0.975)*se), (exp(coef)-1)*100, exp(coef+qnorm(0.975)*se),
                  summary(Harmonicmodel)$coefficients["Effective",4] )
  linearmodel<- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Trend, 
                    family=quasipoisson, temp.data)
  timeestimate <- exp(ci.lin(linearmodel, subset="Time")[,"Estimate"])
  timepval <- summary(linearmodel)$coefficients["Time",4]
  timebottomCI <- exp(ci.lin(linearmodel, subset="Time")[,"2.5%"])
  timetopCI <- exp(ci.lin(linearmodel, subset="Time")[,"97.5%"])
  Before <- rbind(trendbottomCI, trendestimate, trendtopCI, trendpval)
  afterestimate <- exp(ci.lin(linearmodel, subset="Trend")[,"Estimate"])
  trendpval <- summary(linearmodel)$coefficients["Trend",4]
  afterbottomCI <- exp(ci.lin(linearmodel, subset="Trend")[,"2.5%"])
  aftertopCI <- exp(ci.lin(linearmodel, subset="Trend")[,"97.5%"])
  TrendChange <- rbind(afterbottomCI, afterestimate, aftertopCI, trendpval)
  dta3 <- data.frame(Step.Change, Before, TrendChange, type = dataset)
  Table2ModelsRSE <- cbind(Table2ModelsRSE, dta3)
  rownames(Table2ModelsRSE) = c("Bottom CI","Estimate", "Top CI", "P Value")
}
Table2ModelsRSE


#########################
##  INTERACTION MODELS ##
#########################

## ALL Homicide 2016
HomInteract2016 <- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Treatment + Trend + 
                         Effective*Treatment + Time*Treatment + Trend*Treatment + harmonic(MonthNo,2,12),
                       family=quasipoisson, allhom)
ci.lin(HomInteract2016, subset="Effective:Treatment")
ci.lin(HomInteract2016,subset="Time:Treatment")
ci.lin(HomInteract2016,subset="Treatment:Trend")

## Suicides 2016
SuicideInteract2016 <- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Treatment + Trend + 
                             Effective*Treatment + Time*Treatment + Trend*Treatment + harmonic(MonthNo,2,12),
                           family=quasipoisson, allsuicides)
ci.lin(SuicideInteract2016,subset="Effective:Treatment")
ci.lin(SuicideInteract2016,subset="Time:Treatment")
ci.lin(SuicideInteract2016,subset="Treatment:Trend")

## Firearm Homicide 2016
FireHomInteract2016 <- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Treatment + Trend + 
                             Effective*Treatment + Time*Treatment + Trend*Treatment + harmonic(MonthNo,2,12),
                           family=quasipoisson, firehom)
ci.lin(FireHomInteract2016,subset="Effective:Treatment")
ci.lin(FireHomInteract2016,subset="Time:Treatment")
ci.lin(FireHomInteract2016,subset="Treatment:Trend")

## Firearm Suicides 2016
FirearmSuicideInteract2016 <- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Treatment + Trend + 
                                    Effective*Treatment + Time*Treatment + Trend*Treatment + harmonic(MonthNo,2,12),
                                  family=quasipoisson, FASuic)
ci.lin(FirearmSuicideInteract2016,subset="Effective:Treatment")
ci.lin(FirearmSuicideInteract2016,subset="Time:Treatment")
ci.lin(FirearmSuicideInteract2016,subset="Treatment:Trend")

##ALL Firearm related deaths
InteractFireCombo <- glm(Deaths ~ offset(log(StdPop)) + Effective + Time + Treatment + Trend + 
                           Effective*Treatment + Time*Treatment + Trend*Treatment + harmonic(MonthNo,2,12),
                         family=quasipoisson, firecombo)

ci.lin(InteractFireCombo,subset="Effective:Treatment")
ci.lin(InteractFireCombo,subset="Time:Treatment")
ci.lin(InteractFireCombo,subset="Treatment:Trend")




#######################
###                 ###
###     Table 2     ###
###                 ###
#######################

# load data
extdemodata <- read.csv("Table 2 Extension Data.csv")

# add columns
extdemodata$StdPop <- extdemodata$Population / 100000
extdemodata$Rate <- extdemodata$Deaths / extdemodata$StdPop
extdemodata$Treatment <- (rep(1, 9072))
extdemodata$Time <- c(rep(c(1:216), 42))
extdemodata$After <- rep(c(rep(0,81), rep(1,135)),42)
extdemodata$Trend <- rep(c(rep(0,81), 1:135),42)

##  SUBSET THE DATA  ##

# Homicide
e.hom.NHwhite <- subset(extdemodata, Characteristic == "Non-Hispanic White" & Cause == "Homicide")
e.hom.NHblack <- subset(extdemodata, Characteristic == "Non-Hispanic Black" & Cause == "Homicide")
e.hom.hisp <- subset(extdemodata, Characteristic == "Hispanic Any Race" & Cause == "Homicide")
e.hom.under35 <- subset(extdemodata, Characteristic == "Under 35" & Cause == "Homicide")
e.hom.35over <- subset(extdemodata, Characteristic == "35 and over" & Cause == "Homicide")
e.hom.male <- subset(extdemodata, Characteristic == "M" & Cause == "Homicide")
e.hom.female <- subset(extdemodata, Characteristic == "F" & Cause == "Homicide")

# Suicide -- Non-Hispanic Black too small to report
e.suic.NHwhite <- subset(extdemodata, Characteristic == "Non-Hispanic White" & Cause == "Suicide")
e.suic.hisp <- subset(extdemodata, Characteristic == "Hispanic Any Race" & Cause == "Suicide")
e.suic.under35 <- subset(extdemodata, Characteristic == "Under 35" & Cause == "Suicide")
e.suic.35over <- subset(extdemodata, Characteristic == "35 and over" & Cause == "Suicide")
e.suic.male <- subset(extdemodata, Characteristic == "M" & Cause == "Suicide")
e.suic.female <- subset(extdemodata, Characteristic == "F" & Cause == "Suicide")

# Firearm Homicide -- Hispanic and Female too small to report
e.FAhom.NHwhite <- subset(extdemodata, Characteristic == "Non-Hispanic White" & Cause == "Firearm Homicide")
e.FAhom.NHblack <- subset(extdemodata, Characteristic == "Non-Hispanic Black" & Cause == "Firearm Homicide")
e.FAhom.male <- subset(extdemodata, Characteristic == "M" & Cause == "Firearm Homicide")
e.FAhom.under35 <- subset(extdemodata, Characteristic == "Under 35" & Cause == "Firearm Homicide")
e.FAhom.35over <- subset(extdemodata, Characteristic == "35 and over" & Cause == "Firearm Homicide")

# Firearm Suicide -- Non-Hispanic Black and Hispanic too small to report
e.FAsuic.NHwhite <- subset(extdemodata, Characteristic == "Non-Hispanic White" & Cause == "Firearm Suicide")
e.FAsuic.male <- subset(extdemodata, Characteristic == "M" & Cause == "Firearm Suicide")
e.FAsuic.female <- subset(extdemodata, Characteristic == "F" & Cause == "Firearm Suicide")
e.FAsuic.under35 <- subset(extdemodata, Characteristic == "Under 35" & Cause == "Firearm Suicide")
e.FAsuic.35over <- subset(extdemodata, Characteristic == "35 and over" & Cause == "Firearm Suicide")

# All Firearm Deaths
e.allFA.NHwhite <- subset(extdemodata, Characteristic == "Non-Hispanic White" & Cause == "All Firearm")
e.allFA.NHblack <- subset(extdemodata, Characteristic == "Non-Hispanic Black" & Cause == "All Firearm")
e.allFA.hisp <- subset(extdemodata, Characteristic == "Hispanic Any Race" & Cause == "All Firearm")
e.allFA.under35 <- subset(extdemodata, Characteristic == "Under 35" & Cause == "All Firearm")
e.allFA.35over <- subset(extdemodata, Characteristic == "35 and over" & Cause == "All Firearm")
e.allFA.male <- subset(extdemodata, Characteristic == "M" & Cause == "All Firearm")
e.allFA.female <- subset(extdemodata, Characteristic == "F" & Cause == "All Firearm")


#####################################################
##  MEAN MONTHLY COUNTS -- ALONE AND BY POPULATION ##
#####################################################

##  HOMICIDE ##

round(aggregate(e.hom.NHwhite[, "Deaths"], list(e.hom.NHwhite$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.NHwhite[, "Rate"], list(e.hom.NHwhite$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.hom.NHblack[, "Deaths"], list(e.hom.NHblack$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.NHblack[, "Rate"], list(e.hom.NHblack$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.hom.hisp[, "Deaths"], list(e.hom.hisp$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.hisp[, "Rate"], list(e.hom.hisp$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.hom.under35[, "Deaths"], list(e.hom.under35$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.under35[, "Rate"], list(e.hom.under35$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.hom.35over[, "Deaths"], list(e.hom.35over$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.35over[, "Rate"], list(e.hom.35over$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.hom.male[, "Deaths"], list(e.hom.male$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.male[, "Rate"], list(e.hom.male$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.hom.female[, "Deaths"], list(e.hom.female$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.hom.female[, "Rate"], list(e.hom.female$After), FUN = mean, na.rm = T), 2)

##  SUICIDE ##

round(aggregate(e.suic.NHwhite[, "Deaths"], list(e.suic.NHwhite$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.suic.NHwhite[, "Rate"], list(e.suic.NHwhite$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.suic.hisp[, "Deaths"], list(e.suic.hisp$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.suic.hisp[, "Rate"], list(e.suic.hisp$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.suic.under35[, "Deaths"], list(e.suic.under35$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.suic.under35[, "Rate"], list(e.suic.under35$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.suic.35over[, "Deaths"], list(e.suic.35over$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.suic.35over[, "Rate"], list(e.suic.35over$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.suic.male[, "Deaths"], list(e.suic.male$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.suic.male[, "Rate"], list(e.suic.male$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.suic.female[, "Deaths"], list(e.suic.female$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.suic.female[, "Rate"], list(e.suic.female$After), FUN = mean, na.rm = T), 2)

##  FIREARM HOMICIDE ##

round(aggregate(e.FAhom.NHwhite[, "Deaths"], list(e.FAhom.NHwhite$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAhom.NHwhite[, "Rate"], list(e.FAhom.NHwhite$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAhom.NHblack[, "Deaths"], list(e.FAhom.NHblack$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAhom.NHblack[, "Rate"], list(e.FAhom.NHblack$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAhom.hisp[, "Deaths"], list(e.FAhom.hisp$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAhom.hisp[, "Rate"], list(e.FAhom.hisp$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAhom.under35[, "Deaths"], list(e.FAhom.under35$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAhom.under35[, "Rate"], list(e.FAhom.under35$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAhom.35over[, "Deaths"], list(e.FAhom.35over$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAhom.35over[, "Rate"], list(e.FAhom.35over$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAhom.male[, "Deaths"], list(e.FAhom.male$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAhom.male[, "Rate"], list(e.FAhom.male$After), FUN = mean, na.rm = T), 2)

##  FIREARM SUICIDE ##

round(aggregate(e.FAsuic.NHwhite[, "Deaths"], list(e.FAsuic.NHwhite$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAsuic.NHwhite[, "Rate"], list(e.FAsuic.NHwhite$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAsuic.under35[, "Deaths"], list(e.FAsuic.under35$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAsuic.under35[, "Rate"], list(e.FAsuic.under35$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAsuic.35over[, "Deaths"], list(e.FAsuic.35over$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAsuic.35over[, "Rate"], list(e.FAsuic.35over$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAsuic.male[, "Deaths"], list(e.FAsuic.male$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAsuic.male[, "Rate"], list(e.FAsuic.male$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.FAsuic.female[, "Deaths"], list(e.FAsuic.female$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.FAsuic.female[, "Rate"], list(e.FAsuic.female$After), FUN = mean, na.rm = T), 2)


##  ALL FIREARM DEATHS ##

round(aggregate(e.allFA.NHwhite[, "Deaths"], list(e.allFA.NHwhite$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.allFA.NHwhite[, "Rate"], list(e.allFA.NHwhite$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.allFA.NHblack[, "Deaths"], list(e.allFA.NHblack$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.allFA.NHblack[, "Rate"], list(e.allFA.NHblack$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.allFA.under35[, "Deaths"], list(e.allFA.under35$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.allFA.under35[, "Rate"], list(e.allFA.under35$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.allFA.35over[, "Deaths"], list(e.allFA.35over$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.allFA.35over[, "Rate"], list(e.allFA.35over$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.allFA.male[, "Deaths"], list(e.allFA.male$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.allFA.male[, "Rate"], list(e.allFA.male$After), FUN = mean, na.rm = T), 2)

round(aggregate(e.allFA.female[, "Deaths"], list(e.allFA.female$After), FUN = mean, na.rm = T), 2)
round(aggregate(e.allFA.female[, "Rate"], list(e.allFA.female$After), FUN = mean, na.rm = T), 2)


##############################################################
##  Step Change, Trend Change  & Tests for Autocorrelation  ##
##############################################################

###  HOMICIDE  ###

  # WHITE #

## Non-Seasonal Model
m.e.hom.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                   family=quasipoisson, data = e.hom.NHwhite)
# Residual plot
e.hom.NHwhite.res1 <- residuals(m.e.hom.NHwhite,type="deviance")

plot(e.hom.NHwhite$Time,e.hom.NHwhite.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Non-Hispanic White Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.NHwhite)
bgtest(m.e.hom.NHwhite, order=12)

## Seasonal model
sm.e.hom.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                    family=quasipoisson, data = e.hom.NHwhite)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.NHwhite)
bgtest(sm.e.hom.NHwhite, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.NHwhite,Exp=T)["After",5:7],2)  
p <- round(summary(sm.e.hom.NHwhite)$coefficients["After",4],5)     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.NHwhite)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.NHwhite,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.NHwhite)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # BLACK #

## Non-Seasonal Model
m.e.hom.NHblack <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.hom.NHblack)
# Residual plot
e.hom.NHblack.res1 <- residuals(m.e.hom.NHblack,type="deviance")

plot(e.hom.NHblack$Time,e.hom.NHblack.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Non-Hispanic Black Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.NHblack)
bgtest(m.e.hom.NHblack, order=12)

## Seasonal model
sm.e.hom.NHblack <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.hom.NHblack)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.NHblack)
bgtest(sm.e.hom.NHblack, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.NHblack,Exp=T)["After",5:7],2)  
p <- round(summary(sm.e.hom.NHblack)$coefficients["After",4],5)     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.NHblack)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.NHblack,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.NHblack)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # HISPANIC #

## Non-Seasonal Model
m.e.hom.hisp <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.hom.hisp)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.hisp)
bgtest(m.e.hom.hisp, order=12)

## Seasonal model
sm.e.hom.hisp <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.hom.hisp)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.hisp)
bgtest(sm.e.hom.hisp, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.hisp,Exp=T)["After",5:7],2)  
p <- round(summary(sm.e.hom.hisp)$coefficients["After",4],5)     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.hisp)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.hisp,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.hisp)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # UNDER 35 #

## Non-Seasonal Model
m.e.hom.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                    family=quasipoisson, data = e.hom.under35)
# Residual plot
e.hom.under35.res1 <- residuals(m.e.hom.under35,type="deviance", na.rm = T)

plot(e.hom.under35$Time,e.hom.under35.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Under 35 Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.under35)
bgtest(m.e.hom.under35, order=12)

## Seasonal model
sm.e.hom.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = e.hom.under35)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.under35)
bgtest(sm.e.hom.under35, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.under35,Exp=T)["After",5:7],2)  
p <- summary(sm.e.hom.under35)$coefficients["After",4]     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.under35)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.under35,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.under35)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # 35 AND OVER #

## Non-Seasonal Model
m.e.hom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.hom.35over)
# Residual plot
e.hom.35over.res1 <- residuals(m.e.hom.35over,type="deviance", na.rm = T)

plot(e.hom.35over$Time,e.hom.35over.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="35 and Over Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.35over)
bgtest(m.e.hom.35over, order=12)

## Seasonal model
sm.e.hom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.hom.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.35over)
bgtest(sm.e.hom.35over, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.35over,Exp=T)["After",5:7],2)  
p <- summary(sm.e.hom.35over)$coefficients["After",4]     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.35over)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.35over,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.35over)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # MALE #

## Non-Seasonal Model
m.e.hom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                      family=quasipoisson, data = e.hom.male)
# Residual plot
e.hom.male.res1 <- residuals(m.e.hom.male,type="deviance", na.rm = T)

plot(e.hom.male$Time,e.hom.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.male)
bgtest(m.e.hom.male, order=12)

## Seasonal model
sm.e.hom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                       family=quasipoisson, data = e.hom.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.male)
bgtest(sm.e.hom.male, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.male,Exp=T)["After",5:7],2)  
p <- summary(sm.e.hom.male)$coefficients["After",4]     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.male)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.male,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.male)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # FEMALE #

## Non-Seasonal Model
m.e.hom.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                    family=quasipoisson, data = e.hom.female)
# Residual plot
e.hom.female.res1 <- residuals(m.e.hom.female,type="deviance", na.rm = T)

plot(e.hom.female$Time,e.hom.female.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Female Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.hom.female)
bgtest(m.e.hom.female, order=12)

## Seasonal model
sm.e.hom.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = e.hom.female)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.hom.female)
bgtest(sm.e.hom.female, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.hom.female,Exp=T)["After",5:7],2)  
p <- summary(sm.e.hom.female)$coefficients["After",4]    
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.hom.female)["Trend"]),4))
CI <- round(ci.lin(m.e.hom.female,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.hom.female)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))



###  SUICIDE  ###

  # WHITE #

## Non-Seasonal Model
m.e.suic.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.suic.NHwhite)
# Residual plot
e.suic.NHwhite.res1 <- residuals(m.e.suic.NHwhite,type="deviance")

plot(e.suic.NHwhite$Time,e.suic.NHwhite.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Non-Hispanic White Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.suic.NHwhite)
bgtest(m.e.suic.NHwhite, order=12)

## Seasonal model
sm.e.suic.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.suic.NHwhite)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.suic.NHwhite)
bgtest(sm.e.suic.NHwhite, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.suic.NHwhite )["After"]
se <- sqrt(vcovHAC(sm.e.suic.NHwhite )["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.suic.NHwhite)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.suic.NHwhite )["Trend"]
se <- sqrt(vcovHAC(m.e.suic.NHwhite )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.suic.NHwhite)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


  # HISPANIC #

## Non-Seasonal Model
m.e.suic.hisp <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                    family=quasipoisson, data = e.suic.hisp)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.suic.hisp)
bgtest(m.e.suic.hisp, order=12)

## Seasonal model
sm.e.suic.hisp <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = e.suic.hisp)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.suic.hisp)
bgtest(sm.e.suic.hisp, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.suic.hisp,Exp=T)["After",5:7],2)  
p <- round(summary(sm.e.suic.hisp)$coefficients["After",4],5)     
print(c(stepchange, p))

trendchange <- exp(coef(m.e.suic.hisp)["Trend"])
CI <- round(ci.lin(m.e.suic.hisp,Exp=T)["Trend",6:7],5)  
p <- summary(m.e.suic.hisp)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # UNDER 35 #

## Non-Seasonal Model
m.e.suic.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.suic.under35)
# Residual plot
e.suic.under35.res1 <- residuals(m.e.suic.under35,type="deviance", na.rm = T)

plot(e.suic.under35$Time,e.suic.under35.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Under 35 Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.suic.under35)
bgtest(m.e.suic.under35, order=12)

## Seasonal model
sm.e.suic.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.suic.under35)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.suic.under35)
bgtest(sm.e.suic.under35, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.suic.under35,Exp=T)["After",5:7],2)  
p <- summary(sm.e.suic.under35)$coefficients["After",4]     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.suic.under35)["Trend"]),4))
CI <- round(ci.lin(m.e.suic.under35,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.suic.under35)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # 35 AND OVER #

## Non-Seasonal Model
m.e.suic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                      family=quasipoisson, data = e.suic.35over)
# Residual plot
e.suic.35over.res1 <- residuals(m.e.suic.35over,type="deviance", na.rm = T)

plot(e.suic.35over$Time,e.suic.35over.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="35 and Over Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.suic.35over)
bgtest(m.e.suic.35over, order=12)

## Seasonal model
sm.e.suic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                       family=quasipoisson, data = e.suic.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.suic.35over)
bgtest(sm.e.suic.35over, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.suic.35over )["After"]
se <- sqrt(vcovHAC(sm.e.suic.35over)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.suic.35over)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.suic.35over )["Trend"]
se <- sqrt(vcovHAC(m.e.suic.35over )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.suic.35over)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


  # MALE #

## Non-Seasonal Model
m.e.suic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                    family=quasipoisson, data = e.suic.male)
# Residual plot
e.suic.male.res1 <- residuals(m.e.suic.male,type="deviance", na.rm = T)

plot(e.suic.male$Time,e.suic.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.suic.male)
bgtest(m.e.suic.male, order=12)

## Seasonal model
sm.e.suic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = e.suic.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.suic.male)
bgtest(sm.e.suic.male, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.suic.male )["After"]
se <- sqrt(vcovHAC(sm.e.suic.male)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.suic.male)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.suic.male )["Trend"]
se <- sqrt(vcovHAC(m.e.suic.male )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.suic.male)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


  # FEMALE #

## Non-Seasonal Model
m.e.suic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                      family=quasipoisson, data = e.suic.female)
# Residual plot
e.suic.female.res1 <- residuals(m.e.suic.female,type="deviance", na.rm = T)

plot(e.suic.female$Time,e.suic.female.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Female Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.suic.female)
bgtest(m.e.suic.female, order=12)

## Seasonal model
sm.e.suic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                       family=quasipoisson, data = e.suic.female)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.suic.female)
bgtest(sm.e.suic.female, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.suic.female )["After"]
se <- sqrt(vcovHAC(sm.e.suic.female)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.suic.female)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.suic.female )["Trend"]
se <- sqrt(vcovHAC(m.e.suic.female )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.suic.female)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))



###  FIREARM HOMICIDE  ###

  # WHITE #

## Non-Seasonal Model
m.e.FAhom.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.FAhom.NHwhite)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAhom.NHwhite)
bgtest(m.e.FAhom.NHwhite, order=12)

## Seasonal model
sm.e.FAHom.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.FAhom.NHwhite)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAHom.NHwhite)
bgtest(sm.e.FAHom.NHwhite, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.FAHom.NHwhite,Exp=T)["After",5:7],2)  
p <- round(summary(sm.e.FAHom.NHwhite)$coefficients["After",4],5)     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.FAhom.NHwhite)["Trend"]),4))
CI <- round(ci.lin(m.e.FAhom.NHwhite,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.FAhom.NHwhite)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # BLACK #

## Non-Seasonal Model
m.e.FAhom.NHblack <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.FAhom.NHblack)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAhom.NHblack)
bgtest(m.e.FAhom.NHblack, order=12)

## Seasonal model
sm.e.FAhom.NHblack <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.FAhom.NHblack)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAhom.NHblack)
bgtest(sm.e.FAhom.NHblack, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.FAhom.NHblack )["After"]
se <- sqrt(vcovHAC(sm.e.FAhom.NHblack)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.FAhom.NHblack)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.FAhom.NHblack )["Trend"]
se <- sqrt(vcovHAC(m.e.FAhom.NHblack )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),6)  
p <- summary(m.e.FAhom.NHblack)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


  # UNDER 35 #

## Non-Seasonal Model
m.e.FAhom.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.FAhom.under35)
# Residual plot
e.FAhom.under35.res1 <- residuals(m.e.FAhom.under35,type="deviance", na.rm = T)

plot(e.FAhom.under35$Time,e.FAhom.under35.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Under 35 Firearm Homicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAhom.under35)
bgtest(m.e.FAhom.under35, order=12)

## Seasonal model
sm.e.FAhom.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.FAhom.under35)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAhom.under35)
bgtest(sm.e.FAhom.under35, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.FAhom.under35 )["After"]
se <- sqrt(vcovHAC(sm.e.FAhom.under35)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
round(summary(sm.e.FAhom.under35)$coefficients["After", 4], 5)
p <- summary(sm.e.FAhom.under35)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.FAhom.under35 )["Trend"]
se <- sqrt(vcovHAC(m.e.FAhom.under35 )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),6)  
p <- summary(m.e.FAhom.under35)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


  # 35 AND OVER #

## Non-Seasonal Model
m.e.FAhom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                      family=quasipoisson, data = e.FAhom.35over)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAhom.35over)
bgtest(m.e.FAhom.35over, order=12)

## Seasonal model
sm.e.FAhom.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                       family=quasipoisson, data = e.FAhom.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAhom.35over)
bgtest(sm.e.FAhom.35over, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.FAhom.35over )["After"]
se <- sqrt(vcovHAC(sm.e.FAhom.35over)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.FAhom.35over)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.FAhom.35over )["Trend"]
se <- sqrt(vcovHAC(m.e.FAhom.35over )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),5)  
p <- summary(m.e.FAhom.35over)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


  # MALE #

## Non-Seasonal Model
m.e.FAhom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                    family=quasipoisson, data = e.FAhom.male)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAhom.male)
bgtest(m.e.FAhom.male, order=12)

## Seasonal model
sm.e.FAhom.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                     family=quasipoisson, data = e.FAhom.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAhom.male)
bgtest(sm.e.FAhom.male, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.FAhom.male )["After"]
se <- sqrt(vcovHAC(sm.e.FAhom.male)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.FAhom.male)$coefficients["After",4]
print(c(stepchange, p))

coef <- coef(m.e.FAhom.male )["Trend"]
se <- sqrt(vcovHAC(m.e.FAhom.male )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.FAhom.male)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))



###  FIREARM SUICIDE  ###

  # WHITE #

## Non-Seasonal Model
m.e.FAsuic.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                        family=quasipoisson, data = e.FAsuic.NHwhite)
# Residual plot
e.FAsuic.NHwhite.res1 <- residuals(m.e.FAsuic.NHwhite,type="deviance")

plot(e.FAsuic.NHwhite$Time,e.FAsuic.NHwhite.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Non-Hispanic White Firearm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAsuic.NHwhite)
bgtest(m.e.FAsuic.NHwhite, order=12)

## Seasonal model
sm.e.FAsuic.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                         family=quasipoisson, data = e.FAsuic.NHwhite)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAsuic.NHwhite)
bgtest(sm.e.FAsuic.NHwhite, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.FAsuic.NHwhite,Exp=T)["After",5:7],2)  
p <- round(summary(sm.e.FAsuic.NHwhite)$coefficients["After",4],5)     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.FAsuic.NHwhite)["Trend"]),4))
CI <- round(ci.lin(m.e.FAsuic.NHwhite,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.FAsuic.NHwhite)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # UNDER 35 #

## Non-Seasonal Model
m.e.FAsuic.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                        family=quasipoisson, data = e.FAsuic.under35)
# Residual plot
e.FAsuic.under35.res1 <- residuals(m.e.FAsuic.under35,type="deviance", na.rm = T)

plot(e.FAsuic.under35$Time,e.FAsuic.under35.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Under 35 Fierarm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAsuic.under35)
bgtest(m.e.FAsuic.under35, order=12)

## Seasonal model
sm.e.FAsuic.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                         family=quasipoisson, data = e.FAsuic.under35)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAsuic.under35)
bgtest(sm.e.FAsuic.under35, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.FAsuic.under35,Exp=T)["After",5:7],2)  
p <- summary(sm.e.FAsuic.under35)$coefficients["After",4]     
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.FAsuic.under35)["Trend"]),8))
CI <- round(ci.lin(m.e.FAsuic.under35,Exp=T)["Trend",6:7],8)  
p <- summary(m.e.FAsuic.under35)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # 35 AND OVER #

## Non-Seasonal Model
m.e.FAsuic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.FAsuic.35over)
# Residual plot
e.FAsuic.35over.res1 <- residuals(m.e.FAsuic.35over,type="deviance", na.rm = T)

plot(e.FAsuic.35over$Time,e.FAsuic.35over.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="35 and Over Firearm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAsuic.35over)
bgtest(m.e.FAsuic.35over, order=12)

## Seasonal model
sm.e.FAsuic.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.FAsuic.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAsuic.35over)
bgtest(sm.e.FAsuic.35over, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.FAsuic.35over,Exp=T)["After",5:7],2)  # step change w/ 95% CI
p <- summary(sm.e.FAsuic.35over)$coefficients["After",4]     # p-value on "After" = p-value on step change
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.FAsuic.35over)["Trend"]),4))
CI <- round(ci.lin(m.e.FAsuic.35over,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.FAsuic.35over)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # MALE #

## Non-Seasonal Model
m.e.FAsuic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                     family=quasipoisson, data = e.FAsuic.male)
# Residual plot
e.FAsuic.male.res1 <- residuals(m.e.FAsuic.male,type="deviance", na.rm = T)

plot(e.FAsuic.male$Time,e.FAsuic.male.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Male Firearm Suicide Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAsuic.male)
bgtest(m.e.FAsuic.male, order=12)

## Seasonal model
sm.e.FAsuic.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                      family=quasipoisson, data = e.FAsuic.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAsuic.male)
bgtest(sm.e.FAsuic.male, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.FAsuic.male,Exp=T)["After",5:7],2)  # step change w/ 95% CI
p <- summary(sm.e.FAsuic.male)$coefficients["After",4]     # p-value on "After" = p-value on step change
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.FAsuic.male)["Trend"]),4))
CI <- round(ci.lin(m.e.FAsuic.male,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.FAsuic.male)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))


  # FEMALE #

## Non-Seasonal Model
m.e.FAsuic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                       family=quasipoisson, data = e.FAsuic.female)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.FAsuic.female)
bgtest(m.e.FAsuic.female, order=12)

## Seasonal model
sm.e.FAsuic.female <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                        family=quasipoisson, data = e.FAsuic.female)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.FAsuic.female)
bgtest(sm.e.FAsuic.female, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model

# Robust standard errors calculated for the seasonal model
coef <- coef(sm.e.FAsuic.female )["After"]
se <- sqrt(vcovHAC(sm.e.FAsuic.female)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- summary(sm.e.FAsuic.female)$coefficients["After",4]
print(c(stepchange, p))

trendchange <- (round(exp(coef(m.e.FAsuic.female)["Trend"]),4))
CI <- round(ci.lin(m.e.FAsuic.female,Exp=T)["Trend",6:7],4)  
p <- summary(m.e.FAsuic.female)$coefficients["Trend",4]    
print(c((trendchange-1)*100, (CI-1)*100, p))



###  ALL FIREARM DEATHS  ###

# WHITE #

## Non-Seasonal Model
m.e.allFA.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                         family=quasipoisson, data = e.allFA.NHwhite)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.allFA.NHwhite)
bgtest(m.e.allFA.NHwhite, order=12)

## Seasonal model
sm.e.allFA.NHwhite <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                          family=quasipoisson, data = e.allFA.NHwhite)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.allFA.NHwhite)
bgtest(sm.e.allFA.NHwhite, order=12)

# No evidence of statistically significant serial autocorrelation

stepchange <- round(ci.lin(sm.e.allFA.NHwhite,Exp=T)["After",5:7],2) 
p <- round(summary(sm.e.allFA.NHwhite)$coefficients["After",4],7)     
print(c(stepchange, p))

trendchange <- round(ci.lin(m.e.allFA.NHwhite,Exp=T)["Trend",5:7],4)  
p <- round(summary(m.e.allFA.NHwhite)$coefficients["Trend",4],7)     
print(c((trendchange-1)*100, p))


# BLACK #

## Non-Seasonal Model
m.e.allFA.NHblack <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                         family=quasipoisson, data = e.allFA.NHblack)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.allFA.NHblack)
bgtest(m.e.allFA.NHblack, order=12)

## Seasonal model
sm.e.allFA.NHblack <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                          family=quasipoisson, data = e.allFA.NHblack)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.allFA.NHblack)
bgtest(sm.e.allFA.NHblack, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model
coef <- coef(sm.e.allFA.NHblack )["After"]
se <- sqrt(vcovHAC(sm.e.allFA.NHblack)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- round(summary(sm.e.allFA.NHblack)$coefficients["After",4],7)
print(c(stepchange, p))

coef <- coef(m.e.allFA.NHblack )["Trend"]
se <- sqrt(vcovHAC(m.e.allFA.NHblack )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.allFA.NHblack)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


# UNDER 35 #

## Non-Seasonal Model
m.e.allFA.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                         family=quasipoisson, data = e.allFA.under35)
# Residual plot
e.allFA.under35.res1 <- residuals(m.e.allFA.under35,type="deviance", na.rm = T)

plot(e.allFA.under35$Time,e.allFA.under35.res1,ylim=c(-10,10),pch=19,cex=0.7,col=grey(0.6),
     main="Under 35 All Firearm Deaths Residuals over Time",ylab="Deviance residuals",xlab="Date")
abline(h=0,lty=2,lwd=2, col = "red")

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.allFA.under35)
bgtest(m.e.allFA.under35, order=12)

## Seasonal model
sm.e.allFA.under35 <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                          family=quasipoisson, data = e.allFA.under35)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.allFA.under35)
bgtest(sm.e.allFA.under35, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal 

coef <- coef(sm.e.allFA.under35 )["After"]
se <- sqrt(vcovHAC(sm.e.allFA.under35)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- round(summary(sm.e.allFA.under35)$coefficients["After",4],7)
print(c(stepchange, p))

coef <- coef(m.e.allFA.under35 )["Trend"]
se <- sqrt(vcovHAC(m.e.allFA.under35 )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.allFA.under35)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


# 35 AND OVER #

## Non-Seasonal Model
m.e.allFA.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                        family=quasipoisson, data = e.allFA.35over)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.allFA.35over)
bgtest(m.e.allFA.35over, order=12)

## Seasonal model
sm.e.allFA.35over <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                         family=quasipoisson, data = e.allFA.35over)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.allFA.35over)
bgtest(sm.e.allFA.35over, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.allFA.35over )["After"]
se <- sqrt(vcovHAC(sm.e.allFA.35over)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- round(summary(sm.e.allFA.35over)$coefficients["After",4],7)
print(c(stepchange, p))

coef <- coef(m.e.allFA.35over )["Trend"]
se <- sqrt(vcovHAC(m.e.allFA.35over )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.allFA.35over)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))


# MALE #

## Non-Seasonal Model
m.e.allFA.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend, 
                      family=quasipoisson, data = e.allFA.male)

# Breusch-Godfrey test for autocorrelation
bgtest(m.e.allFA.male)
bgtest(m.e.allFA.male, order=12)

## Seasonal model
sm.e.allFA.male <- glm(Deaths ~ offset(log(StdPop)) + After + Time + Trend + harmonic(MonthNo,2,12), 
                       family=quasipoisson, data = e.allFA.male)

# Seasonal Breusch-Godfrey test for autocorrelation
bgtest(sm.e.allFA.male)
bgtest(sm.e.allFA.male, order=12)

# Statistically significant serial autocorrelation evident in seasonal and non-seasonal model
# Robust standard errors calculated for the seasonal model

coef <- coef(sm.e.allFA.male )["After"]
se <- sqrt(vcovHAC(sm.e.allFA.male)["After","After"])
stepchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),2)
p <- round(summary(sm.e.allFA.male)$coefficients["After",4],7)
print(c(stepchange, p))

coef <- coef(m.e.allFA.male )["Trend"]
se <- sqrt(vcovHAC(m.e.allFA.male )["Trend","Trend"])
trendchange <- round(c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se)),4)  
p <- summary(m.e.allFA.male)$coefficients["Trend",4]    
print(c((trendchange-1)*100, p))



#########################################################
###                     Figure 2.A                    ###
###   Homicide v. Suicide among Hispanic Floridians   ###
#########################################################

# Data Set-Up
e.men <- subset(extdemodata, Characteristic == "M" & Cause == "Suicide" | 
                      Characteristic == "M" & Cause == "Firearm Suicide" )

# Setup variables
e.men$Level <- c(rep(0,81),rep(1,135),rep(0,81),rep(1,135))
e.men$Firearm <- c(rep(1, 216), rep(0, 216))
e.men$FATime <- e.men$Firearm * e.men$Time
e.men$FALevel <- e.men$Firearm * e.men$Level
e.men$FATrend <- e.men$Firearm * e.men$Trend

# Models
e.men_linear <- glm(Rate ~ offset(log(StdPop)) + Time + Level + Trend + Firearm + FATime + FALevel + FATrend, 
                    family=quasipoisson, data=e.men)
e.men_harmonic <- glm(Rate ~ offset(log(StdPop)) + Time + Level + Trend + Firearm + FATime + FALevel + FATrend + 
                        harmonic(MonthNo,2,12), family=quasipoisson, data=e.men)
e.men.mod <- with(e.men, Deaths/StdPop)
e.men.datanew <- data.frame(StdPop=e.men$StdPop,Level=e.men$Level, Time = e.men$Time/10, MonthNo = e.men$MonthNo/10)

e.men.pred1 <- predict(e.men_linear,type="response",e.men)/e.men$StdPop
e.men.pred2 <- predict(e.men_linear,type="response",transform(e.men,MonthNo=4.8))/(e.men$StdPop)

# Plot
pdf("Fig2A_MaleSuicvFAsuic.pdf")
plot(e.men$Time[1:216],e.men.mod[1:216],
     type="n",
     ylim=c(0,3),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:17*12+12,tick=T,labels=1999:2016) 
rect(81,0,216,3, col=grey(0.9),border=F)
points(e.men$Time[1:216],e.men.mod[1:216],
       col="darkorange2",
       pch=20)
points(e.men$Time[217:432],e.men.mod[217:432],
       col="dodgerblue3",
       pch=20)
lines(e.men$Time[1:81], fitted(e.men_linear)[1:81], col="darkorange",lty=5)
lines(e.men$Time[1:81], fitted(e.men_harmonic)[1:81], col="darkorange",lwd=2)
lines(e.men$Time[82:216], fitted(e.men_linear)[82:216], col="darkorange",lty=5)
lines(e.men$Time[82:216], fitted(e.men_harmonic)[82:216], col="darkorange",lwd=2)
lines(e.men$Time[217:297], fitted(e.men_linear)[217:297], col="dodgerblue3",lty=5)
lines(e.men$Time[217:297], fitted(e.men_harmonic)[217:297], col="dodgerblue3",lwd=2)
lines(e.men$Time[298:432], fitted(e.men_linear)[298:432], col="dodgerblue3",lty=5)
lines(e.men$Time[298:432], fitted(e.men_harmonic)[298:432], col="dodgerblue3",lwd=2)
legend(x=0, y=2.9, legend=c("Suicide","Firearm Suicide"), col=c("dodgerblue3","darkorange"),lwd = 2)
title(main="Change in suicide and firearm suicide \n among male Floridians")
graphics.off()


####################################################################
###                         Figure 2.B                           ###
###   Suicide v. Firearm Suicide among Floridians Under Age 35   ###
####################################################################

# Data Set-Up
e.under35 <- subset(extdemodata, Characteristic == "Under 35" & Cause == "Suicide" | 
                   Characteristic == "Under 35" & Cause == "Firearm Suicide" )

# Setup variables
e.under35$Level <- c(rep(0,81),rep(1,135),rep(0,81),rep(1,135))
e.under35$Firearm <- c(rep(0, 216), rep(1, 216))
e.under35$FSTime <- e.under35$Firearm * e.under35$Time
e.under35$FSLevel <- e.under35$Firearm * e.under35$Level
e.under35$FSTrend <- e.under35$Firearm * e.under35$Trend

# Models
e.under35_linear <- glm(Rate ~ offset(log(StdPop)) + Time + Level + Trend + Firearm + FSTime + FSLevel + FSTrend, 
                     family=quasipoisson, data=e.under35)
e.under35_harmonic <- glm(Rate ~ offset(log(StdPop)) + Time + Level + Trend + Firearm + FSTime + FSLevel + FSTrend + 
                         harmonic(MonthNo,2,12), family=quasipoisson, data=e.under35)
e.under35.mod <- with(e.under35, Deaths/StdPop)
e.under35.datanew <- data.frame(StdPop=e.under35$StdPop,Level=e.under35$Level, Time = e.under35$Time/10, MonthNo = e.under35$MonthNo/10)

e.under35.pred1 <- predict(e.under35_linear,type="response",e.under35)/e.under35$StdPop
e.under35.pred2 <- predict(e.under35_linear,type="response",transform(e.under35,MonthNo=4.8))/(e.under35$StdPop)

# Plot
pdf("Under35HomvFAhom.pdf")
plot(e.under35$Time[1:216],e.under35.mod[1:216],
     type="n",
     ylim=c(0,1.5),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:17*12+12,tick=T,labels=1999:2016) 
rect(81,0,216,1.5, col=grey(0.9),border=F)
points(e.under35$Time[1:216],e.under35.mod[1:216],
       col="darkorange2",
       pch=20)
points(e.under35$Time[217:432],e.under35.mod[217:432],
       col="dodgerblue3",
       pch=20)
lines(e.under35$Time[1:81], fitted(e.under35_linear)[1:81], col="darkorange",lty=5)
lines(e.under35$Time[1:81], fitted(e.under35_harmonic)[1:81], col="darkorange",lwd=2)
lines(e.under35$Time[82:216], fitted(e.under35_linear)[82:216], col="darkorange",lty=5)
lines(e.under35$Time[82:216], fitted(e.under35_harmonic)[82:216], col="darkorange",lwd=2)
lines(e.under35$Time[217:297], fitted(e.under35_linear)[217:297], col="dodgerblue3",lty=5)
lines(e.under35$Time[217:297], fitted(e.under35_harmonic)[217:297], col="dodgerblue3",lwd=2)
lines(e.under35$Time[298:432], fitted(e.under35_linear)[298:432], col="dodgerblue3",lty=5)
lines(e.under35$Time[298:432], fitted(e.under35_harmonic)[298:432], col="dodgerblue3",lwd=2)
legend(x=0, y=1.4, legend=c("Suicide","Firearm Suicide"), col=c("darkorange","dodgerblue3"),lwd = 2)
title(main="Change in suicide and firearm suicide \n among Floridians under age 35")
graphics.off()


##########################################################
###                     Figure 2.C                     ###
###   Firearm Deaths among White v. Black Floridians   ###
##########################################################

# Data Set-Up
e.BWallFA <- subset(extdemodata, Characteristic == "Non-Hispanic White" & Cause == "All Firearm" | 
                      Characteristic == "Non-Hispanic Black" & Cause == "All Firearm" )

# Setup variables
e.BWallFA$Level <- c(rep(0,81),rep(1,135),rep(0,81),rep(1,135))
e.BWallFA$Black <- c(rep(1, 216), rep(0, 216))
e.BWallFA$BTime <- e.BWallFA$Black * e.BWallFA$Time
e.BWallFA$BLevel <- e.BWallFA$Black * e.BWallFA$Level
e.BWallFA$BTrend <- e.BWallFA$Black * e.BWallFA$Trend

# Models
e.BWallFA_linear <- glm(Rate ~ offset(log(StdPop)) + Time + Level + Trend + Black + BTime + BLevel + BTrend, 
                        family=quasipoisson, data=e.BWallFA)
e.BWallFA_harmonic <- glm(Rate ~ offset(log(StdPop)) + Time + Level + Trend + Black + BTime + BLevel + BTrend + 
                            harmonic(MonthNo,2,12), family=quasipoisson, data=e.BWallFA)
e.BWallFA.mod <- with(e.BWallFA, Deaths/StdPop)
e.BWallFA.datanew <- data.frame(StdPop=e.BWallFA$StdPop,Level=e.BWallFA$Level, Time = e.BWallFA$Time/10, MonthNo = e.BWallFA$MonthNo/10)

e.BWallFA.pred1 <- predict(e.BWallFA_linear,type="response",e.BWallFA)/e.BWallFA$StdPop
e.BWallFA.pred2 <- predict(e.BWallFA_linear,type="response",transform(e.BWallFA,MonthNo=4.8))/(e.BWallFA$StdPop)

# Plot
pdf("BlackWhiteallFA.pdf")
plot(e.BWallFA$Time[1:216],e.BWallFA.mod[1:216],
     type="n",
     ylim=c(0,2.5),
     frame.plot=F,
     ylab="Deaths per 100,000",
     xlab="Year",
     xaxt="n",
     las=2)
axis(1,at=0:17*12+12,tick=T,labels=1999:2016) 
rect(81,0,216,2.5, col=grey(0.9),border=F)
points(e.BWallFA$Time[1:216],e.BWallFA.mod[1:216],
       col="darkorange2",
       pch=20)
points(e.BWallFA$Time[217:432],e.BWallFA.mod[217:432],
       col="dodgerblue3",
       pch=20)
lines(e.BWallFA$Time[1:80], fitted(e.BWallFA_linear)[1:80], col="darkorange",lty=5)
lines(e.BWallFA$Time[1:80], fitted(e.BWallFA_harmonic)[1:80], col="darkorange",lwd=2)
lines(e.BWallFA$Time[81:215], fitted(e.BWallFA_linear)[81:215], col="darkorange",lty=5)
lines(e.BWallFA$Time[81:215], fitted(e.BWallFA_harmonic)[81:215], col="darkorange",lwd=2)
lines(e.BWallFA$Time[217:297], fitted(e.BWallFA_linear)[217:297], col="dodgerblue3",lty=5)
lines(e.BWallFA$Time[217:297], fitted(e.BWallFA_harmonic)[217:297], col="dodgerblue3",lwd=2)
lines(e.BWallFA$Time[298:432], fitted(e.BWallFA_linear)[298:432], col="dodgerblue3",lty=5)
lines(e.BWallFA$Time[298:432], fitted(e.BWallFA_harmonic)[298:432], col="dodgerblue3",lwd=2)
legend(x=0, y=2.4, legend=c("Black","White"), col=c("darkorange","dodgerblue3"),lwd = 2)
title(main="Change in firearm deaths \n among White v. Black Floridians")
graphics.off()
