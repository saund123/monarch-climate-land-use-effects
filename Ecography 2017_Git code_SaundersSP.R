#############################################################################
#############Data manipulation and model fitting code########################
##Local and cross-seasonal associations of climate and land use with 
##abundance of monarch butterflies Danaus plexippus
##Authors: Saunders, SP, Ries, L, Oberhauser, KS, Thogmartin, WE, Zipkin, EF
##Data sources are proprietary. Please contact authors for more information
##Note: only monarch observations were used from datasets of all butterfly
##species observations in Ohio and Illinois monitoring networks
############################################################################

#load libraries
library(jagsUI)
library(reshape)

#set working directory to where input files are located on user computer
setwd("C:/Users/")

#read in data 
#SurveyData file consists of file with observations by row and the following columns:
#Latitude, longitude, Site ID, Year, Week, Sum of Monarch Counts, Sum of Total Time
mon[1:10,]
names(mon)
str(mon)

#Look at data (checking all years/weeks/sites present)
uyear = sort(unique(mon$Year))
#94-13
uweek= sort(unique(mon$Week))
#5-35
usite = sort(unique(mon$SiteID))
#total sites: 262 (unpruned).

#Reshape the count data.
#NOTE: This assumes columns called Year, Week, SiteID and SumofMonarch count
#Does not need to be in any specific order
junk.melt=melt(mon,id.var=c("Year", "Week", "SiteID"), measure.var="SumOfMonarch.count")
monarchs=cast(junk.melt, SiteID ~ Week ~ Year, sum)
#This organizes into tables that are sums of counts according to Week by SiteID for each year

#need to put in NAs when site wasn't surveyed a given week of a given year
#If surveyed in a given week, insert count (including 0)
dates=read.csv("Weeks_IL94to14.csv", header=T,sep=',',na.strings=T) 
head(dates)
tail(dates)
names(dates)
dim(dates) #262 sites (matchs mon)

for (k in 1:length(uweek)){
  for (tt in 1:length(uyear)) {
    a=which(dates$Year==uyear[tt] & dates$Week==uweek[k]) 
    b= which(dates[a,]==0)
    bb=b-2
    monarchs[bb,k,tt]=NA
  }}
monarchs

#Include only spring/summer (e.g., weeks 5-29)
#Remove weeks 30-35 from the data
monarchs1=monarchs[,1:25,]
#keep weeks 5-29

########################
####DATA MANIP##########
########################
a=matrix(NA, ncol=126, nrow=20)

for (t in 1:length(uyear)) {
  b=monarchs[,14:25,t]
  cc= rep(0,262)
  
  for (j in 1:length(usite)) {
    cc[j]=length(which(b[j,]>(-1)))
  }
  
  dd=which(cc>0)
  ee= 126-length(dd)
  
  a[t,] = c(rownames(b[dd,]), rep(NA,ee))
}
#necessary for nested indexing to index over sites with fewer than 2 surveys/peak period

#need to pull from a[] how many sites surveyed per year and put in a vector with 20 spaces
for (t in 1:length(uyear)){
  ff[t]=sum(!is.na(a[t,])) 
}
usite.new<-ff #vector of how many sites surveyed at least once (during peak period) in given year

newmon = array(NA, c(126, length(uweek[14:25]),length(uyear)))

for (t in 1:length(uyear)){
  newmon[1:usite.new[t],,t] = monarchs[a[t,1:usite.new[t]],14:25,t]
}

#newmon is count data for re-indexed site only

mon3 = array(NA, c(126, length(uweek[14:25]),length(uyear)))
week3 = array(NA, c(126, length(uweek[14:25]),length(uyear)))
uweek.new = array(NA, c(126, length(uyear)))

for (t in 1:length(uyear)) {
  for (j in 1:usite.new[t]) {
    a = which(newmon[j,,t] > (-1))
    mon3[j,1:length(a),t] = newmon[j,a,t]
    week3[j,1:length(a),t] = a
    uweek.new[j,t] = length(a)
  }
}

#mon3 is left aligned. week3 is new week covariate (need to standardize below)
#uweek.new is new indexing of week to use in model loop

#####################################
###########Adding covariates#########
#####################################

#Site effects
coord=read.csv("SiteEffects_IL94to14.csv", header=T,sep=',',na.strings=T)
head(coord)
str(coord)
dim(coord) #262 rows (matches mon)
#this file is organized according to sites as rows and the following columns:
#longitude, latitude, site ID, percent open, avgerage GDD, Crop cover

#HABITAT COVARIATES
#Site effects--constant across years

#% open habitat (standardized)
amatch=pmatch(a[1,1:usite.new[1]],coord$SiteID)
openyr1=coord$X.Open[amatch]
#this reports a vector of the correct open values associated with surveyed sites for a given year
#repeat for 1 through 20 years

l<-list(openyr1,openyr2,openyr3,openyr4,openyr5,openyr6,openyr7,openyr8,openyr9,openyr10,openyr11,openyr12,openyr13,openyr14,openyr15,openyr16,openyr17,openyr18,openyr19,openyr20)

#creating matrix of open values for all years
n <- max(sapply(l, length)) #this returns length of longest vector
ll <- lapply(l, function(X) {
  c(as.numeric(X), rep(NA, times = n - length(X)))
}) #loop function fills in NA so that all vectors become equal length
newopen <- do.call(cbind, ll)

#then standardize
openmean=mean(newopen,na.rm=T)
opensd=sd(as.vector(newopen), na.rm=T)
newopen.st<-round((newopen-openmean)/opensd,3)
#replacing NA with 0 because covariate data can't have NAs
newopen.st[is.na(newopen.st)] <- 0

#final open covariate for model: newopen.st

#Average GDD
amatch=pmatch(a[1,1:usite.new[1]],coord$SiteID)
avggddyr1=coord$avgGDD10.28NEW[amatch]
#this reports a vector of the correct avggdd values associated with surveyed sites for a given year
#repeat for 1 through 20 years

l<-list(avggddyr1,avggddyr2,avggddyr3,avggddyr4,avggddyr5,avggddyr6,avggddyr7,avggddyr8,avggddyr9,avggddyr10,avggddyr11,avggddyr12,avggddyr13,avggddyr14,avggddyr15,avggddyr16,avggddyr17,avggddyr18,avggddyr19,avggddyr20)

#create matrix of avggdd values for all years
n <- max(sapply(l, length)) #this returns length of longest vector
ll <- lapply(l, function(X) {
  c(as.numeric(X), rep(NA, times = n - length(X)))
}) #loop function fills in NA so that all vectors become equal length
newavggdd <- do.call(cbind, ll) 

#then standardize
mavggdd=mean(newavggdd,na.rm=T)
savggdd=sd(as.vector(newavggdd), na.rm=T)
newavggdd.st<-round((newavggdd-mavggdd)/savggdd,3)
#replacing NA with 0
newavggdd.st[is.na(newavggdd.st)] <- 0

#final avggdd covariate to include in model: newavggdd.st
#also include quadratic term: newavggdd.st2

##Crop Cover
amatch=pmatch(a[1,1:usite.new[1]],coord$SiteID)
ccoveryr1=coord$CropCov[amatch]
#this reports a vector of the correct open values associated with surveyed sites for a given year
#repeat for 1 through 20 years

l<-list(ccoveryr1,ccoveryr2,ccoveryr3,ccoveryr4,ccoveryr5,ccoveryr6,ccoveryr7,ccoveryr8,ccoveryr9,ccoveryr10,ccoveryr11,ccoveryr12,ccoveryr13,ccoveryr14,ccoveryr15,ccoveryr16,ccoveryr17,ccoveryr18,ccoveryr19,ccoveryr20)

#creating matrix of crop cov values for all years
n <- max(sapply(l, length)) #this returns length of longest vector
ll <- lapply(l, function(X) {
  c(as.numeric(X), rep(NA, times = n - length(X)))
}) #loop function fills in NA so that all vectors become equal length
newcropcov <- do.call(cbind, ll) 

#then standardize
ccovmean=mean(newcropcov,na.rm=T)
ccovsd=sd(as.vector(newcropcov), na.rm=T)
newcropcov.st<-round((newcropcov-ccovmean)/ccovsd,3)
#replacing NA with 0 
newcropcov.st[is.na(newcropcov.st)] <- 0

#final crop cover covariate for model: newcropcov.st

#SITE BY YEAR effects
#NOTE: sites are in numerical order in this file
siteyear = read.csv("SiteYear_IL94to14_09max.csv", header=T,sep=',',na.strings=T)
str(siteyear)

#organize site by year
#report appropriate PDSI values according to siteID and year
#File is organized with each row as site-year combinations; thus, columns are:
#Site ID, Year, average PDSI, Glyphosate Application, Corn & Soy area, other crop area
sitepalm = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)) {
  for (j in 1:length(usite)) {
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    sitepalm[j,tt] = siteyear$AvgPDSI10.28[a] 
  }} 

newpalm=matrix(NA, nrow=126, ncol=20)

for (t in 1:length(uyear)) {
  junk=sitepalm[,t]
  names(junk) = usite
  match = pmatch(a[t,1:usite.new[t]],names(junk))
  newpalm[1:usite.new[t],t]= junk[match] 
}

#then standardize
newpalmmean=mean(newpalm,na.rm=T)
newpalmsd=sd(as.vector(newpalm), na.rm=T)
newpalm.st<-round((newpalm-newpalmmean)/newpalmsd,3)
#replacing NA with 0 
newpalm.st[is.na(newpalm.st)] <- 0

#final sitepalm covariate for model: newpalm.st

#report county-level app rates according to siteID and year
sitegly = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)){
  for (j in 1:length(usite)){
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    sitegly[j,tt] = siteyear$GlyApp[a]
  }
}

newgly=matrix(NA, nrow=126, ncol=20)

for (t in 1:length(uyear)) {
  junk=sitegly[,t]
  names(junk) = usite
  match = pmatch(a[t,1:usite.new[t]],names(junk))
  newgly[1:usite.new[t],t]= junk[match] 
}

#then standardize
newglymean=mean(newgly,na.rm=T)
newglysd=sd(as.vector(newgly), na.rm=T)
newgly.st<-round((newgly-newglymean)/newglysd,3)
#replacing NA with 0 
newgly.st[is.na(newgly.st)] <- 0

#final sitegly covariate for model: newgly.st

#report county-level app rates according to siteID and year. THIS IS % OF CROPS SPRAYED
sitegly3 = matrix(0,nrow=length(usite), ncol=length(uyear))
for (tt in 1:length(uyear)){
  for (j in 1:length(usite)){
    a=which(siteyear$SiteID == usite[j] & siteyear$Year == uyear[tt])
    sitegly3[j,tt] = siteyear$GlyApp3[a]
  }
}
#note:re run a from data manip above before running below code
newgly3=matrix(NA, nrow=126, ncol=20)

for (t in 1:length(uyear)) {
  junk=sitegly3[,t]
  names(junk) = usite
  match = pmatch(a[t,1:usite.new[t]],names(junk))
  newgly3[1:usite.new[t],t]= junk[match] 
}

#then standardize
newglymean3=mean(newgly3,na.rm=T)
newglysd3=sd(as.vector(newgly3), na.rm=T)
newgly.st3<-round((newgly3-newglymean3)/newglysd3,3)
#replacing NA with 0 
newgly.st3[is.na(newgly.st3)] <- 0
#note: correlation between crop cover (99-12) and newgly.st3(99-12) is 0.117. between crop and newgly.st3 (all years): 0.064

#YEAR effects (TX GDD and drought)
#NOTE: This file is first year of study (first row) to last year of study (bottom row) with columns:
#Year, spring GDD, spring precipitation, winter colony size
yeareffects=read.csv("YearEffects_IL94to14.csv", header=T,sep=',',na.strings=T)
str(yeareffects)
names(yeareffects)

spGDD=yeareffects$SprGDD4to9
spPrec=yeareffects$SprPrecipTX

#standardize spring GDD
mspGDD=mean(as.matrix(spGDD))
sdspGDD=sd(as.vector(spGDD))
spGDD1=as.vector((spGDD-mspGDD)/sdspGDD)

#standardize spring Precip
mspPrec=mean(as.matrix(spPrec))
sdspPrec=sd(as.vector(spPrec))
spPrec1=as.vector((spPrec-mspPrec)/sdspPrec)

#December colony size
earlywinter2=yeareffects$T2PrColSize

#standardize early winter colony size
mearlywinter2=mean(as.matrix(earlywinter2))
sdearlywinter2=sd(as.vector(earlywinter2))
earlywinter2.1=as.vector((earlywinter2-mearlywinter2)/sdearlywinter2)

#standardize week
week3mean=mean(week3,na.rm=T)
week3sd=sd(as.vector(week3), na.rm=T)
week3.st<-round((week3-week3mean)/week3sd,3)
#replacing NA with 0 
week3.st[is.na(week3.st)] <- 0

#use week3.st as final week covariate

###############################################
### Specify model in BUGS language#############
##############################################

sink("Ecogmodel.jags")
cat("
    model {
    a1  ~ dnorm(0,0.01)
    a2  ~ dnorm(0,0.01)
    a3  ~ dnorm(0,0.01)
    a4  ~ dnorm(0,0.01)
    a5  ~ dnorm(0,0.01)
    a6  ~ dnorm(0,0.01)
    a7  ~ dnorm(0,0.01)
    a8  ~ dnorm(0,0.01)
    a9  ~ dnorm(0,0.01)
    a10  ~ dnorm(0,0.01)
    a11  ~ dnorm(0,0.01)
    a12  ~ dnorm(0,0.01)
    a14  ~ dnorm(0,0.01)
    a15  ~ dnorm(0,0.01)
    
    rprec ~ dgamma(0.01,0.01)
    
    tau.a13 ~ dgamma(0.01, 0.01)
    
    for (m in 1:usite){
    a13[m] ~ dnorm(0, tau.a13)
    }
    
    
    for (t in 1:uyear){
    for (j in 1:usite.new[t]){
    for (k in 1:uweek.new[j,t]){
    
    y[j,k,t]  ~ dnegbin(p[j,k,t],rprec)
    
    p[j,k,t] <- rprec/(rprec + mu[j,k,t])
    
    log(mu[j,k,t]) <- a1 + a2*spPrec1[t] + a3*spPrec2[t] + 
    a4*spGDD1[t] + a5*spGDD2[t] + 
    a6*newavggdd.st[j,t] +
    a7*newpalm.st[j,t] +
    a8*newopen.st[j,t] + a9*earlywinter2.1[t] + 
    a10*newcropcov.st[j,t]+ a14*newgly.st3[j,t] + a15*newgly.st3[j,t]*newcropcov.st[j,t] +
    a11*week3.st[j,k,t] + a12*week3.st[j,k,t]*week3.st[j,k,t] + a13[siteRE2[j,t]]
    
    
    }
    }
    }
    
    }
    ", fill=TRUE)
sink()

#PREP JAGSUI DATA
bugsdata<-list(uyear=length(uyear), usite.new=usite.new, uweek.new=uweek.new, usite=length(usite),
               y=mon3, week3.st=week3.st, siteRE2=siteRE2,
               newcropcov.st=newcropcov.st,
               newpalm.st=newpalm.st, earlywinter2.1=earlywinter2.1,
               spGDD1=spGDD1, newopen.st=newopen.st, spPrec1=spPrec1, 
               spGDD2=spGDD1*spGDD1, newgly.st3=newgly.st3,
               spPrec2=spPrec1*spPrec1, newavggdd.st=newavggdd.st)

#inits for NEG BINOM model
inits<-function(){
  list(rprec=rgamma(1,1))
}

parameters<-c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 
              'a10', 'a11','a12', 'a14', 'a15', 'tau.a13')

############################
##RUN BUGS MODEL IN JAGSUI
############################
out<-jags(data = bugsdata,
          inits = inits,
          parameters.to.save = parameters,
          model.file = 'Ecogmodel.jags',
          n.chains = 3,
          n.adapt = 100,
          n.iter = 6000,
          n.burnin = 3000,
          n.thin = 3,
          parallel=TRUE
)
print(output,digits=4)
whiskerplot(output,parameters=c("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a14","a15"))

##################################################
##################################################
##Partial run: first 10 years: re-standardize data
##################################################
##################################################

#changing mon3 to be for first 10 years only (1994-2003)
dim(mon3)
mon3c<-mon3[,,1:10]

#changing uweek.new to be for first 10 years only
uweek.newc<-uweek.new[,1:10]

###COVARIATES###########
#Open covariate
#change newopen first
dim(newopen)
newopenc<-newopen[,1:10]
dim(newopenc)

#then standardize
openmeanc=mean(newopenc,na.rm=T)
opensdc=sd(as.vector(newopenc), na.rm=T)
newopen.stc<-round((newopenc-openmeanc)/opensdc,3)
#replacing NA with 0 because covariate data can't have NAs
newopen.stc[is.na(newopen.stc)] <- 0
#new covariate to use: newopen.stc

#Avg GDD
#change newavggdd first
dim(newavggdd)
newavggddc<-newavggdd[,1:10]
dim(newavggddc)

#then standardize
mavggddc=mean(newavggddc,na.rm=T)
savggddc=sd(as.vector(newavggddc), na.rm=T)
newavggdd.stc<-round((newavggddc-mavggddc)/savggddc,3)
#replacing NA with 0 (though shouldn't matter since won't iterate over those values)
newavggdd.stc[is.na(newavggdd.stc)] <- 0
#new covariate to use: newavggdd.stc

#Crop cover
dim(newcropcov)
newcropcovc<-newcropcov[,1:10]
dim(newcropcovc)

#then standardize
ccovmeanc=mean(newcropcovc,na.rm=T)
ccovsdc=sd(as.vector(newcropcovc), na.rm=T)
newcropcov.stc<-round((newcropcovc-ccovmeanc)/ccovsdc,3)
#replacing NA with 0 
newcropcov.stc[is.na(newcropcov.stc)] <- 0
#new covariate to use: newcropcov.stc

###SITE YEAR COVARIATES######
##PDI#########
dim(newpalm)
newpalmc<-newpalm[,1:10]
dim(newpalmc)

#then standardize
newpalmmeanc=mean(newpalmc,na.rm=T)
newpalmsdc=sd(as.vector(newpalmc), na.rm=T)
newpalm.stc<-round((newpalmc-newpalmmeanc)/newpalmsdc,3)
#replacing NA with 0 
newpalm.stc[is.na(newpalm.stc)] <- 0
#new covariate to use: newpalm.stc

#GLY 3 (% sprayed)##########
dim(newgly3)
newgly3c<-newgly3[,1:10]
dim(newgly3c)

#then standardize
newglymean3c=mean(newgly3c,na.rm=T)
newglysd3c=sd(as.vector(newgly3c), na.rm=T)
newgly.st3c<-round((newgly3c-newglymean3c)/newglysd3c,3)
#replacing NA with 0 
newgly.st3c[is.na(newgly.st3c)] <- 0
#new covariate to use: newgly.st3c

####YEAR EFFECTS############
#Spring precip and temps

#re-standardize spring GDD with only 10 years
mspGDDc=mean(as.matrix(spGDD[1:10]))
sdspGDDc=sd(as.vector(spGDD[1:10]))
spGDD1c=as.vector((spGDD[1:10]-mspGDDc)/sdspGDDc)
length(spGDD1c)
#new covariate to use: spGDD1c

#standardize spring Precip
mspPrecc=mean(as.matrix(spPrec[1:10]))
sdspPrecc=sd(as.vector(spPrec[1:10]))
spPrec1c=as.vector((spPrec[1:10]-mspPrecc)/sdspPrecc)
#new covariate to use: spPrec1c

##Early winter
#standardize early winter colony size
mearlywinter2c=mean(as.matrix(earlywinter2[1:10]))
sdearlywinter2c=sd(as.vector(earlywinter2[1:10]))
earlywinter2.1c=as.vector((earlywinter2[1:10]-mearlywinter2c)/sdearlywinter2c)
#new covariate to use: earlywinter2.1c

###WEEKS COVARIATE############
dim(week3)
week3c<-week3[,,1:10]
dim(week3c)

#then standardize
week3meanc=mean(week3c,na.rm=T)
week3sdc=sd(as.vector(week3c), na.rm=T)
week3.stc<-round((week3c-week3meanc)/week3sdc,3)
#replacing NA with 0 
week3.stc[is.na(week3.stc)] <- 0
dim(week3.stc)
#final covariate to use: week3.stc

##Adjusting random effect
dim(siteRE2)
siteRE2c<-siteRE2[,1:10]
dim(siteRE2c)

#########################
#PREP JAGSUI DATA
bugsdata<-list(uyear=length(uyear[1:10]), usite.new=usite.new[1:10], uweek.new=uweek.newc, usite=length(usite),
               y=mon3c, week3.st=week3.stc, siteRE2=siteRE2c,
               newcropcov.st=newcropcov.stc,  
               newpalm.st=newpalm.stc, earlywinter2.1=earlywinter2.1c,
               spGDD1=spGDD1c, newopen.st=newopen.stc, spPrec1=spPrec1c,
               spGDD2=spGDD1c*spGDD1c, newgly.st3=newgly.st3c, 
               spPrec2=spPrec1c*spPrec1c, newavggdd.st=newavggdd.stc)

#inits for NEG BINOM model
inits<-function(){
  list(rprec=rgamma(1,1))
}

parameters<-c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 
              'a10', 'a11','a12', 'a14', 'a15', 'tau.a13')

############################
##RUN BUGS MODEL IN JAGSUI
############################
output.first10<-jags(data = bugsdata,
                    inits = inits,
                    parameters.to.save = parameters,
                    model.file = 'Ecogmodel.jags',
                    n.chains = 3,
                    n.adapt = 100,
                    n.iter = 6000,
                    n.burnin = 3000,
                    n.thin = 3,
                    parallel=TRUE
)
print(output.first10,digits=4)

#################################################
#################################################
##Partial run: last 10 years: re-standardize data
#################################################
#################################################

#changing mon3 to be for last 10 years only (2004-2013)
dim(mon3)
mon3b<-mon3[,,11:20]

#changing uweek.new to be for last 10 years only
uweek.newb<-uweek.new[,11:20]

###COVARIATES#############
#Open covariate
#change newopen first
dim(newopen)
newopenb<-newopen[,11:20]
dim(newopenb)

#then standardize
openmeanb=mean(newopenb,na.rm=T)
opensdb=sd(as.vector(newopenb), na.rm=T)
newopen.stb<-round((newopenb-openmeanb)/opensdb,3)
#replacing NA with 0 because covariate data can't have NAs
newopen.stb[is.na(newopen.stb)] <- 0
#new covariate to use: newopen.stb

#Avg GDD
#change newavggdd first
dim(newavggdd)
newavggddb<-newavggdd[,11:20]
dim(newavggddb)

#then standardize
mavggddb=mean(newavggddb,na.rm=T)
savggddb=sd(as.vector(newavggddb), na.rm=T)
newavggdd.stb<-round((newavggddb-mavggddb)/savggddb,3)
#replacing NA with 0
newavggdd.stb[is.na(newavggdd.stb)] <- 0
#new covariate to use: newavggdd.stb

#Crop cover
dim(newcropcov)
newcropcovb<-newcropcov[,11:20]
dim(newcropcovb)

#then standardize
ccovmeanb=mean(newcropcovb,na.rm=T)
ccovsdb=sd(as.vector(newcropcovb), na.rm=T)
newcropcov.stb<-round((newcropcovb-ccovmeanb)/ccovsdb,3)
#replacing NA with 0 
newcropcov.stb[is.na(newcropcov.stb)] <- 0
#new covariate to use: newcropcov.stb

###SITE YEAR COVARIATES####################
##PDI#########
dim(newpalm)
newpalmb<-newpalm[,11:20]
dim(newpalmb)

#then standardize
newpalmmeanb=mean(newpalmb,na.rm=T)
newpalmsdb=sd(as.vector(newpalmb), na.rm=T)
newpalm.stb<-round((newpalmb-newpalmmeanb)/newpalmsdb,3)
#replacing NA with 0 
newpalm.stb[is.na(newpalm.stb)] <- 0
#new covariate to use: newpalm.stb

#GLY 3 (% sprayed)###################
dim(newgly3)
newgly3b<-newgly3[,11:20]
dim(newgly3b)

#then standardize
newglymean3b=mean(newgly3b,na.rm=T)
newglysd3b=sd(as.vector(newgly3b), na.rm=T)
newgly.st3b<-round((newgly3b-newglymean3b)/newglysd3b,3)
#replacing NA with 0 
newgly.st3b[is.na(newgly.st3b)] <- 0
#new covariate to use: newgly.st3b

####YEAR EFFECTS##########
#Spring precip and temps

#re-standardize spring GDD with only 10 years
mspGDDb=mean(as.matrix(spGDD[11:20]))
sdspGDDb=sd(as.vector(spGDD[11:20]))
spGDD1b=as.vector((spGDD[11:20]-mspGDDb)/sdspGDDb)
length(spGDD1b)
#new covariate to use: spGDD1b

#standardize spring Precip
mspPrecb=mean(as.matrix(spPrec[11:20]))
sdspPrecb=sd(as.vector(spPrec[11:20]))
spPrec1b=as.vector((spPrec[11:20]-mspPrecb)/sdspPrecb)
#new covariate to use: spPrec1b

##Early winter
#standardize early winter colony size
mearlywinter2b=mean(as.matrix(earlywinter2[11:20]))
sdearlywinter2b=sd(as.vector(earlywinter2[11:20]))
earlywinter2.1b=as.vector((earlywinter2[11:20]-mearlywinter2b)/sdearlywinter2b)
#new covariate to use: earlywinter2.1b

#Late winter (2005-2013)
latewinter<-c(11.12,1.62,6.28,3.51,0.89,3.4,0.77,2.95,1.19,0.89)
#standardize early winter colony size
mlatewinter=mean(as.matrix(latewinter))
sdlatewinter=sd(as.vector(latewinter))
latewinter1=as.vector((latewinter-mlatewinter)/sdlatewinter)
#new covariate to use: latewinter1

###WEEKS COVARIATE############
dim(week3)
week3b<-week3[,,11:20]
dim(week3b)

#then standardize
week3meanb=mean(week3b,na.rm=T)
week3sdb=sd(as.vector(week3b), na.rm=T)
week3.stb<-round((week3b-week3meanb)/week3sdb,3)
#replacing NA with 0 
week3.stb[is.na(week3.stb)] <- 0
dim(week3.stb)
#final covariate to use: week3.stb

##Adjusting random effect
dim(siteRE2)
siteRE2b<-siteRE2[,11:20]
dim(siteRE2b)

###########################
#PREP JAGSUI DATA
bugsdata<-list(uyear=length(uyear[11:20]), usite.new=usite.new[11:20], uweek.new=uweek.newb, usite=length(usite),
               y=mon3b, week3.st=week3.stb, siteRE2=siteRE2b,
               newcropcov.st=newcropcov.stb,
               newpalm.st=newpalm.stb, earlywinter2.1=earlywinter2.1b, #latewinter1=latewinter1, 
               spGDD1=spGDD1b, newopen.st=newopen.stb, spPrec1=spPrec1b,
               spGDD2=spGDD1b*spGDD1b, newgly.st3=newgly.st3b, 
               spPrec2=spPrec1b*spPrec1b, newavggdd.st=newavggdd.stb)

#inits for NEG BINOM model
inits<-function(){
  list(rprec=rgamma(1,1))
}

#NOTE: change # of parameters acc. to model structure!
parameters<-c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 
              'a10', 'a11','a12', 'a14', 'a15', 'tau.a13')


############################
##RUN BUGS MODEL IN JAGSUI
############################
output.last10<-jags(data = bugsdata,
                    inits = inits,
                    parameters.to.save = parameters,
                    model.file = 'Ecogmodel.jags',
                    n.chains = 3,
                    n.adapt = 100,
                    n.iter = 6000,
                    n.burnin = 3000,
                    n.thin = 3,
                    parallel=TRUE
)
print(output.last10,digits=4)

