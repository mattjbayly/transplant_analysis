#### PROJECT: Mimulus cardinalis northern transplant
#### PURPOSE: Obtain bootstrapped vital rate coefficients
############# IPMs parameterized by bootstrapped coeffcients will then be run and lambdas will be obtained
############# Replicate bootstrap datasets will be used to obtain confidence intervals around lambda estimates for each site
#### AUTHOR: Seema Sheth
#### DATE LAST MODIFIED: 20180224

# remove objects and clear workspace
rm(list = ls(all=TRUE))

# require packages
require(plyr)
require(dplyr)
require(lme4)
require(glmmADMB)
require(doParallel)

# set working directory
# setwd("/Users/ssheth/Google Drive/demography_PNAS_November2017")

#*******************************************************************************
#### 1. Preliminaries ###
#*******************************************************************************

# Read in & examine bootstrapped data 
bootstrapped.data=readRDS("Robjects/Mcard_transplant_INDIV_BOOTSTRAP_data.rds")
str(bootstrapped.data)

# set # of cores to use for parallel processing
 registerDoParallel(cores=4)

# Set number of bootstrap replicate datasets
n.boot=2000

#*******************************************************************************
#### 2. Survival probability ###
#*******************************************************************************

# create function to run survival model
surv <- function(i) {
  glmer(Surv~logSize+SiteID+(1|PlotID),data=subset(bootstrapped.data,Replicate==i),family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))  
}

# run mixed effects model for survival on each bootstrapped dataset
surv.reg <- lapply(1:2000, surv)

# save model output
saveRDS(surv.reg,"Robjects/surv.reg_boot.rds") # 50+ convergence warnings

#*******************************************************************************
#### 3. Growth ###
#*******************************************************************************

# Create list to store regression objects
results=c()
growth.reg=list()

# parallel processing
results=foreach(i = 1:n.boot,.export=c('lmer'), .packages=c('lme4')) %dopar% {
  
  # Inspect model
  growth.reg[[i]]=lmer(logSizeNext~logSize+SiteID+(1|PlotID),data=subset(bootstrapped.data,Replicate==i),control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))   
} # End loop

# save model output
saveRDS(results,"Robjects/growth.reg_boot.rds")

#*******************************************************************************
#### 4. Flowering probability ###
#*******************************************************************************

# create function to run survival model
flowering <- function(i) {
  glmer(Fec0~logSize+SiteID+(1|PlotID),data=subset(bootstrapped.data,Replicate==i),family=binomial,control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))  
}

# run mixed effects model for survival on each bootstrapped dataset
flowering.reg <- lapply(1:2000, flowering)

# save model output
saveRDS(flowering.reg,"Robjects/flowering.reg_boot.rds") # 50+ convergence warnings

#*******************************************************************************
#### 5. Fruit # ###
#*******************************************************************************

# Create list to store regression objects
fruit.reg=list()
results=c()

# parallel processing
results=foreach(i = 1:n.boot,.export=c('glmmadmb'), .packages=c('glmmADMB')) %dopar% {
  
  # Inspect model
  fruit.reg[[i]]=glmmadmb(Fec1~logSize+SiteID+(1|PlotID),data=subset(bootstrapped.data,Replicate==i&!is.na(Fec1)),family="nbinom",link="log")   
} # End loop

# save model outputs
saveRDS(results,"Robjects/fruit.reg_boot.rds")




