#### PROJECT: Mimulus cardinalis northern transplant 2015-2016
#### PURPOSE: Perform model selection for each vital rate for subsequent use in IPMs
############# Vital rates include survival, growth, flowering, and fruit count
############# Fixed effects: size, site; Random effects: plot
#### AUTHOR: modified from Seema Sheth (Sheth & Angert 2018 PNAS)
#### DATE LAST MODIFIED: 20180224

# remove objects and clear workspace
rm(list = ls(all=TRUE))

# Install GLMMADMB package following instructions here: http://glmmadmb.r-forge.r-project.org/
# install.packages("R2admb")
# install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos",getOption("repos")),type="source")

#require packages
require(lme4)
require(glmmADMB)
require(MuMIn)
require(MASS)
require(pscl)
require(plyr)
require(dplyr)
require(tidyverse)


#*******************************************************************************
#### 1. Import and format data ###
#*******************************************************************************

data <- read_csv("Data/IPMData_transplant.csv")
head(data)

# log-transform size
data <- data %>% mutate(logSize = log(z), logSizeNext = log(z1), Fec0 = Repr, Fec1 = Fec)

# make sure plot and site are recognized as factors
data$PlotID = as.factor(data$PlotID)
data$SiteID = as.factor(data$SiteID)
data$Year = as.factor(data$Year)

# Variables are: 

# Surv: survival (1) or not (0) of individuals between time = t and time = t+1 
# logSize: total stem length of the individual at time t
# logSizeNext: total stem length of the individual at time t+1
# Repr (= Seema's Fec0): Flowering yes (1) or no (0)
# Fec (= Seema's Fec1): Total number of fruits per individual   
# SiteID: population
# PlotID: plot
# Region: within or beyond northern range edge
# ID: unique identifier for each individual
# NewPlot_13: don't know what this is
# NewPlot_14: ditto
# Year: annual transition 

# Not included in Matt's file but could bind in from other sources?

# Latitude: latitude of population
# Longitude: longitude of population
# Elevation: elevation of population
# SeedCt: mean seed count, rounded to the nearest integer, for each site
# ClassNext: stage class (juvenile, adult, dead, or NA) of plant at time = t+1 

# simple summaries for Supplementary Table 1
countNperplot <- data %>% 
  dplyr::group_by(PlotID) %>% 
  dplyr::summarize(count_per_plot = n_distinct(ID))
min(countNperplot$count_per_plot)
max(countNperplot$count_per_plot)
mean(countNperplot$count_per_plot)
sum(countNperplot$count_per_plot)

countNpersite <- data %>% 
  dplyr::group_by(SiteID) %>% 
  dplyr::summarize(count_per_site = n_distinct(ID))
min(countNpersite$count_per_site)
max(countNpersite$count_per_site)
mean(countNpersite$count_per_site)
sum(countNpersite$count_per_site)

#*******************************************************************************
#### 2. Survival ###
#*******************************************************************************

# fixed effects model w/ and w/out size
s1=glm(Surv~logSize,data=data,family=binomial)
s2=glm(Surv~1,data=data,family=binomial)
model.sel(s1,s2) # model w/ size is preferred

# site as fixed or random
sf=glm(Surv~SiteID,data=data,family=binomial)
sr=glmer(Surv~(1|SiteID),data=data,family=binomial)
model.sel(sf,sr) # fixed is preferred (is this a legit comparison??)

# A. interaction size x site; random intercepts & random slopes for Plot
s3=glmer(Surv~logSize*SiteID+(logSize|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# B. interaction size x site; random intercepts & constant slope for Plot
s4=glmer(Surv~logSize*SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# C. main effects size + site; random intercepts & random slopes for Plot
s5=glmer(Surv~logSize+SiteID+(logSize|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 
# warnings!

# D. main effects size + site; random intercepts & constant slope for Plot
s6=glmer(Surv~logSize+SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# E. interaction size x site; random intercepts & random slopes for Plot nested within Site
s7=glmer(Surv~logSize*SiteID+(z|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 
# warnings

# F. interaction size x site; random intercepts & constant slope for Plot nested within Site
s8=glmer(Surv~logSize*SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# G. main effects size + site; random intercepts & random slopes for Plot nested within Site
s9=glmer(Surv~logSize+SiteID+(logSize|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# H. main effects size + site; random intercepts & constant slope for Plot nested within Site
s10=glmer(Surv~logSize+SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# Compare models
anova(s3, s4, s5, s6, s7, s8, s9, s10)

AICc(s3, s4, s5, s6, s7, s8, s9, s10) 

model.sel(s3, s4, s5, s6, s7, s8, s9, s10) 

# PREFERRED MODEL IS s6

# Save top survival model to .rda file
save(s6, file='Robjects/surv.reg.rda')   

#*******************************************************************************
#### 3. Growth ###
#*******************************************************************************

# fixed effects model w/ and w/out size
g1=glm(logSizeNext~logSize,data=data[!is.na(data$logSize),])
g2=glm(logSizeNext~1,data=data[!is.na(data$logSize),])
model.sel(g1,g2) # model w/ size is preferred

# A. interaction Site x Size; random intercepts & random slopes for Plot
g3=lmer(logSizeNext~logSize*SiteID+(logSize|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# B. interaction Site x Size; random intercepts & random slopes for Plot
g4=lmer(logSizeNext~logSize*SiteID+(1|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# C. main effects Site + Size; random intercepts & random slopes for Plot
g5=lmer(logSizeNext~logSize+SiteID+(logSize|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# D. main effects Site + Size; random intercepts & constant slope for Plot
g6=lmer(logSizeNext~logSize+SiteID+(1|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# E. interaction Site x Size; random intercepts & random slopes for Plot nested within Site
g7=lmer(logSizeNext~logSize*SiteID+(logSize|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) # warnings

# F. interaction Site x Size; random intercepts & random slopes for Plot nested within Site
g8=lmer(logSizeNext~logSize*SiteID+(1|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# G. main effects Site + Size; random intercepts & random slopes for Plot nested within Site
g9=lmer(logSizeNext~logSize+SiteID+(logSize|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
# warnings

# H. main effects Site + Size; random intercepts & constant slope for Plot nested within Site
g10=lmer(logSizeNext~logSize+SiteID+(1|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 
# warnings

# Compare models
anova(g3, g4, g5, g6, g7, g8, g9, g10)

AICc(g3, g4, g5, g6, g7, g8, g9, g10) 

model.sel(g3, g4, g5, g6, g7, g8, g9, g10) 

# PREFERRED MODEL IS g6

# Save top growth model to .rda file
save(g6, file='Robjects/growth.reg.rda')   

#*******************************************************************************
#### 4. Flowering ###
#*******************************************************************************

# fixed effects model w/ and w/out size
fl1=glm(Fec0~logSize,data=data,family=binomial)
fl2=glm(Fec0~1,data=data,family=binomial)
model.sel(fl1,fl2) # model w/ size is preferred

# A. interaction Site x Size; random intercepts & random slopes for Plot
fl3=glmer(Fec0~logSize*SiteID+(logSize|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# B. interaction Site x Size; random intercepts & constant slope for Plot
fl4=glmer(Fec0~logSize*SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# C. main effects Site + Size; random intercepts & random slopes for Plot
fl5=glmer(Fec0~logSize+SiteID+(logSize|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# D. main effects Site + Size; random intercepts & constant slope for Plot
fl6=glmer(Fec0~logSize+SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# E. interaction Site x Size; random intercepts & random slopes for Plot nested within Site
fl7=glmer(Fec0~logSize*SiteID+(logSize|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# F. interaction Site x Size; random intercepts & constant slope for Plot nested within Site
fl8=glmer(Fec0~logSize*SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# G. main effects Site + Size; random intercepts & random slopes for Plot nested within Site
fl9=glmer(Fec0~logSize+SiteID+(logSize|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# H. main effects Site + Size; random intercepts & constant slope for Plot nested within Site
fl10=glmer(Fec0~logSize+SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 
# warning

# I. main effects Site + Size; random intercepts & constant slope for Plot; random intercepts and random slopes for Year
fl11=glmer(Fec0~logSize+SiteID+(1|PlotID)+(logSize|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# J. main effects Site + Size; random intercepts & constant slope for Plot; random intercepts and constant slope for Year
fl12=glmer(Fec0~logSize+SiteID+(1|PlotID)+(1|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# K. main effects Site + Size; random intercepts & random slopes for Plot; random intercepts and random slopes for Year
fl13=glmer(Fec0~logSize+SiteID+(logSize|PlotID)+(logSize|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# L. main effects Site + Size; random intercepts & random slopes for Plot; random intercepts and constant slope for Year
fl14=glmer(Fec0~logSize+SiteID+(logSize|PlotID)+(1|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# M. interaction Site x Size; random intercepts & constant slope for Plot; random intercepts and random slopes for Year
fl15=glmer(Fec0~logSize*SiteID+(1|PlotID)+(logSize|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# N. interaction Site x Size; random intercepts & constant slope for Plot; random intercepts and constant slope for Year
fl16=glmer(Fec0~logSize*SiteID+(1|PlotID)+(1|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# O. interaction Site x Size; random intercepts & random slopes for Plot; random intercepts and random slopes for Year
fl17=glmer(Fec0~logSize*SiteID+(logSize|PlotID)+(logSize|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# P. interaction Site x Size; random intercepts & random slopes for Plot; random intercepts and constant slope for Year
fl18=glmer(Fec0~logSize*SiteID+(logSize|PlotID)+(1|Year),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# Compare models
anova(fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10,fl11,fl12,fl13,fl14,fl15,fl16,fl17,fl18)

AICc(fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10,fl11,fl12,fl13,fl14,fl15,fl16,fl17,fl18) 

model.sel(fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10,fl11,fl12,fl13,fl14,fl15,fl16,fl17,fl18) 
model.sel(fl6, fl12) 

# PREFERRED MODEL IS fl6

# Save top flowering model to .rda file
save(fl6, file='Robjects/flowering.reg.rda')   

#*******************************************************************************
#### 5. Fruit number ###
#*******************************************************************************
  
	   #*******************************************************************************
	   ####5A. Fit fixed effects models (GLMs) only for initial exploratory model selection of variance structure (poisson vs. negative binomial)
	   #*******************************************************************************
	  
	   fr1=glm(Fec1~logSizeNext,data=data,na.action=na.omit,family=poisson) # poisson without 0-inflation 
	   fr2=glm.nb(Fec1~logSizeNext,data=data,na.action=na.omit) # negative binomial without 0-inflation
	   #fr3=zeroinfl(Fec~z,data=data,na.action=na.omit,dist="poisson") # poisson with 0-inflation
	   #fr4=zeroinfl(Fec~z,data=data,na.action=na.omit,dist="negbin") # negative binomial with 0-inflation
	   model.sel(fr1,fr2)#,fr3,fr4) 
	   # note: couldn't try 0-inflated models because min value of Fec = 1
	   
	   # PREFERRED MODEL IS fr2 (negative binomial w/out 0-inflation)
	   # fixed effects model w/ and w/out size
	  	fr5=glm.nb(Fec1~1,data=data,na.action=na.omit)
	  	model.sel(fr2,fr5) # model w/ size is preferred
	  	
		#*******************************************************************************
		####5B. Model selection of interaction and random effects structure
		#*******************************************************************************
		
		# A. interaction of Size x Site; random intercepts & random slopes for Plot
		fr6=glmmadmb(Fec1~logSize*SiteID+(logSize|PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
		
		# B. interaction of Size x Site; random intercepts & constant slope for Plot
		fr7=glmmadmb(Fec1~logSize*SiteID+(1|PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
		
		# C. main effects of Size + Site; random intercepts & random slopes for Plot
		fr8=glmmadmb(Fec1~logSize+SiteID+(logSize|PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
		
		# D. main effects of Size + Site; random intercepts & constant slope for Plot
		fr9=glmmadmb(Fec1~logSize+SiteID+(1|PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# E. interaction of Size x Site; random intercepts & random slopes for Plot nested within Site
		fr10=glmmadmb(Fec1~logSize*SiteID+(logSize|SiteID/PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
		
		# F. interaction of Size x Site; random intercepts & constant slope for Plot nested within Site
		fr11=glmmadmb(Fec1~logSize*SiteID+(1|SiteID/PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
		
		# G. main effects of Size + Site; random intercepts & random slopes for Plot nested within Site
		fr12=glmmadmb(Fec1~logSize+SiteID+(logSize|SiteID/PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log") 
		
		# H. main effects of Size + Site; random intercepts & constant slope for Plot nested within Site
		fr13=glmmadmb(Fec1~logSize+SiteID+(1|SiteID/PlotID),data=data[!is.na(data$Fec1),],family="nbinom",link="log")

		# I. main effects of Size + Site; random intercepts & constant slope for Plot; random intercepts and random slopes for Year
		fr14=glmmadmb(Fec1~logSize+SiteID+(1|PlotID)+(logSize|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# J. main effects of Size + Site; random intercepts & constant slope for Plot; random intercepts and constant slope for Year
		fr15=glmmadmb(Fec1~logSize+SiteID+(1|PlotID)+(1|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")

		# K. main effects of Size + Site; random intercepts & random slopes for Plot; random intercepts and random slopes for Year
		fr16=glmmadmb(Fec1~logSize+SiteID+(logSize|PlotID)+(logSize|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# L. main effects of Size + Site; random intercepts & random slopes for Plot; random intercepts and constant slope for Year
		fr17=glmmadmb(Fec1~logSize+SiteID+(logSize|PlotID)+(1|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")

		# M. interaction Size x Site; random intercepts & constant slope for Plot; random intercepts and random slopes for Year
		fr18=glmmadmb(Fec1~logSize*SiteID+(1|PlotID)+(logSize|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# N. interaction Size x Site; random intercepts & constant slope for Plot; random intercepts and constant slope for Year
		fr19=glmmadmb(Fec1~logSize*SiteID+(1|PlotID)+(1|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# O. interaction Size x Site; random intercepts & random slopes for Plot; random intercepts and random slopes for Year
		fr20=glmmadmb(Fec1~logSize*SiteID+(1|PlotID)+(logSize|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# P. intearction Size x Site; random intercepts & random slopes for Plot; random intercepts and constant slope for Year
		fr21=glmmadmb(Fec1~logSize*SiteID+(1|PlotID)+(1|Year),data=data[!is.na(data$Fec1),],family="nbinom",link="log")
		
		# Compare models	
		AICc(fr6,fr7,fr8,fr9,fr10,fr11,fr12,fr13,fr14,fr15,fr16,fr17,fr18,fr19,fr20,fr21) 
		
		model.sel(fr6,fr7,fr8,fr9,fr10,fr11,fr12,fr13,fr14,fr15,fr16,fr17,fr18,fr19,fr20,fr21)
		model.sel(fr9,fr15)
		
		# PREFERRED MODEL IS fr15
		
		# Save top fruit # model to .rda file because it takes a long time to run
		save(fr15, file='Robjects/fruit.reg.rda')   
		