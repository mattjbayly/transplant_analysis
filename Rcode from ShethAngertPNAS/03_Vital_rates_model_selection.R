#### PROJECT: Mimulus cardinalis northern transplant 2015-2016
#### PURPOSE: Perform model selection for each vital rate for subsequent use in IPMs
############# Vital rates include survival, growth, flowering, and fruit count
############# Fixed effects: size, site; Random effects: plot
#### AUTHOR: modified from Seema Sheth (Sheth & Angert 2018 PNAS)
#### DATE LAST MODIFIED: 20180220

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

# set working directory
# setwd("/Users/ssheth/Google Drive/demography_PNAS_November2017")

#*******************************************************************************
#### 1. Import and format data ###
#*******************************************************************************

data <- read_csv("Data/IPMData_transplant.csv")
head(data)

# log-transform size
data <- data %>% mutate(z = log(z), z1 = log(z1))

# make sure plot and site are recognized as factors
data$PlotID = as.factor(data$PlotID)
data$SiteID = as.factor(data$SiteID)

# Variables are: 

# Surv: survival (1) or not (0) of individuals between time = t and time = t+1 
# z (= Seema's logSize): total stem length of the individual
# z1 (= Seema's logSizeNext): same as "logSize" above, for t+1
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

#*******************************************************************************
#### 2. Survival ###
#*******************************************************************************

# fixed effects model w/ and w/out size
s1=glm(Surv~z,data=data,family=binomial)
s2=glm(Surv~1,data=data,family=binomial)
model.sel(s1,s2) # model w/ size is preferred

# site as fixed or random
sf=glm(Surv~SiteID,data=data,family=binomial)
sr=glmer(Surv~(1|SiteID),data=data,family=binomial)
model.sel(sf,sr) # fixed is preferred (is this a legit comparison??)

# A. interaction size x site; random intercepts & random slopes for Plot
s3=glmer(Surv~z*SiteID+(z|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# B. interaction size x site; random intercepts & constant slope for Plot
s4=glmer(Surv~z*SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# C. main effects size + site; random intercepts & random slopes for Plot
s5=glmer(Surv~z+SiteID+(z|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 
# warnings!

# D. main effects size + site; random intercepts & constant slope for Plot
s6=glmer(Surv~z+SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# E. interaction size x site; random intercepts & random slopes for Plot nested within Site
s7=glmer(Surv~z*SiteID+(z|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# F. interaction size x site; random intercepts & constant slope for Plot nested within Site
s8=glmer(Surv~z*SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# G. main effects size + site; random intercepts & random slopes for Plot nested within Site
s9=glmer(Surv~z+SiteID+(z|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# H. main effects size + site; random intercepts & constant slope for Plot nested within Site
s10=glmer(Surv~z+SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

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
g1=glm(z1~z,data=data[!is.na(data$z),])
g2=glm(z1~1,data=data[!is.na(data$z),])
model.sel(g1,g2) # model w/ size is preferred

# A. interaction Site x Size; random intercepts & random slopes for Plot
g3=lmer(z1~z*SiteID+(z|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# B. interaction Site x Size; random intercepts & random slopes for Plot
g4=lmer(z1~z*SiteID+(1|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# C. main effects Site + Size; random intercepts & random slopes for Plot
g5=lmer(z1~z+SiteID+(z|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# D. main effects Site + Size; random intercepts & constant slope for Plot
g6=lmer(z1~z+SiteID+(1|PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 

# E. interaction Site x Size; random intercepts & random slopes for Plot nested within Site
g7=lmer(z1~z*SiteID+(z|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) # warnings

# F. interaction Site x Size; random intercepts & random slopes for Plot nested within Site
g8=lmer(z1~z*SiteID+(1|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 
# G. main effects Site + Size; random intercepts & random slopes for Plot nested within Site
g9=lmer(z1~z+SiteID+(z|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
# warnings

# H. main effects Site + Size; random intercepts & constant slope for Plot nested within Site
g10=lmer(z1~z+SiteID+(1|SiteID/PlotID),data=data,control=lmerControl(optimizer = "bobyqa")) 
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
fl1=glm(Repr~z,data=data,family=binomial)
fl2=glm(Repr~1,data=data,family=binomial)
model.sel(fl1,fl2) # model w/ size is preferred

# A. interaction Site x Size; random intercepts & random slopes for Plot
fl3=glmer(Repr~z*SiteID+(z|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# B. interaction Site x Size; random intercepts & constant slope for Plot
fl4=glmer(Repr~z*SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# C. main effects Site + Size; random intercepts & random slopes for Plot
fl5=glmer(Repr~z+SiteID+(z|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# D. main effects Site + Size; random intercepts & constant slope for Plot
fl6=glmer(Repr~z+SiteID+(1|PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# E. interaction Site x Size; random intercepts & random slopes for Plot nested within Site
fl7=glmer(Repr~z*SiteID+(z|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# F. interaction Site x Size; random intercepts & constant slope for Plot nested within Site
fl8=glmer(Repr~z*SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 

# G. main effects Site + Size; random intercepts & random slopes for Plot nested within Site
fl9=glmer(Repr~z+SiteID+(z|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa")) 

# H. main effects Site + Size; random intercepts & constant slope for Plot nested within Site
fl10=glmer(Repr~z+SiteID+(1|SiteID/PlotID),data=data,family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000))) 
# warning

# Compare models
anova(fl3, fl4, fl5, fl6, fl7, fl8, fl9, fl10)

AICc(fl3, fl4, fl5, fl6, fl7, fl8, fl9, fl10) 

model.sel(fl3, fl4,fl5,fl6, fl7, fl8, fl9, fl10) 

# PREFERRED MODEL IS fl6

# Save top flowering model to .rda file
save(fl6, file='Robjects/flowering.reg.rda')   

#*******************************************************************************
#### 5. Fruit number ###
#*******************************************************************************
  
	   #*******************************************************************************
	   ####5A. Fit fixed effects models (GLMs) only for initial exploratory model selection of variance structure (poisson vs. negative binomial)
	   #*******************************************************************************
	  
	   fr1=glm(Fec~z,data=data,na.action=na.omit,family=poisson) # poisson without 0-inflation 
	   fr2=glm.nb(Fec~z,data=data,na.action=na.omit) # negative binomial without 0-inflation
	   #fr3=zeroinfl(Fec~z,data=data,na.action=na.omit,dist="poisson") # poisson with 0-inflation
	   #fr4=zeroinfl(Fec~z,data=data,na.action=na.omit,dist="negbin") # negative binomial with 0-inflation
	   model.sel(fr1,fr2)#,fr3,fr4) 
	   # note: couldn't try 0-inflated models because min value of Fec = 1
	   
	   # PREFERRED MODEL IS fr2 (negative binomial w/out 0-inflation)
	   # fixed effects model w/ and w/out size
	  	fr5=glm.nb(Fec~1,data=data,na.action=na.omit)
	  	model.sel(fr2,fr5) # model w/ size is preferred
	  	
		#*******************************************************************************
		####5B. Model selection of interaction and random effects structure
		#*******************************************************************************
		
		# A. interaction of Size x Site; random intercepts & random slopes for Plot
		fr6=glmmadmb(Fec~z*SiteID+(z|PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log") 
		
		# B. interaction of Size x Site; random intercepts & constant slope for Plot
		fr7=glmmadmb(Fec~z*SiteID+(1|PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log") 
		
		# C. main effects of Size + Site; random intercepts & random slopes for Plot
		fr8=glmmadmb(Fec~z+SiteID+(z|PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log") 
		
		# D. main effects of Size + Site; random intercepts & constant slope for Plot
		fr9=glmmadmb(Fec~z+SiteID+(1|PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log")
		
		# E. interaction of Size x Site; random intercepts & random slopes for Plot nested within Site
		fr10=glmmadmb(Fec~z*SiteID+(z|SiteID/PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log") 
		
		# F. interaction of Size x Site; random intercepts & constant slope for Plot nested within Site
		fr11=glmmadmb(Fec~z*SiteID+(1|SiteID/PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log") 
		
		# G. main effects of Size + Site; random intercepts & random slopes for Plot nested within Site
		fr12=glmmadmb(Fec~z+SiteID+(z|SiteID/PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log") 
		
		# H. main effects of Size + Site; random intercepts & constant slope for Plot nested within Site
		fr13=glmmadmb(Fec~z+SiteID+(1|SiteID/PlotID),data=data[!is.na(data$Fec),],family="nbinom",link="log")

				# Compare models	
		AICc(fr6,fr7,fr8,fr9,fr10,fr11,fr12,fr13) 
		
		model.sel(fr6,fr7,fr8,fr9,fr10,fr11,fr12,fr13)
		
		# PREFERRED MODEL IS fr9
		
		# Save top fruit # model to .rda file because it takes a long time to run
		save(fr9, file='Robjects/fruit.reg.rda')   
		