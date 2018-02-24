#### PROJECT: Mimulus cardinalis northern transplant 2015-2016
#### PURPOSE: Create data frame of vital rate parameters and build integral projection models 
############# Obtain estimates of lambda for each transplant site
#### AUTHOR: modified from Seema Sheth (Sheth & Angert 2018 PNAS)
#### DATE LAST MODIFIED: 20171110

# remove objects and clear workspace
rm(list = ls(all=TRUE))

# require packages
require(lme4)
require(glmmADMB)
require(plyr)
require(dplyr)
require(tidyverse)
require(stringr)

#*******************************************************************************
#### 1. Import and format data ###
#*******************************************************************************
# these data are used to fill in the gap from fruits to recruits
# from Sheth and Angert 2018 from survey of natural populations

# individual-level data
# run data prep script to clean up and import data
source("Rcode/data_prep.R")
head(demo.data)
head(site_fruit_count_data)

# filter to northern sites only
demo.dat.north <- demo.data %>% 
  filter(Site=="Rock Creek"|Site=="Canton Creek"|Site=="Coast Fork of Williamette") %>% 
  droplevels()
head(demo.dat.north)

site_fruit_count_north <- site_fruit_count_data %>% 
  filter(Site=="Rock Creek"|Site=="Canton Creek"|Site=="Coast Fork of Williamette") %>% 
  droplevels()
head(site_fruit_count_north)

# seeds per fruit
seeds <- read_csv("Data/Amy_wild_demo_data/fec2.seed.per.fruit.2010.2011.2012.csv")
head(seeds)

# filter to northern sites only
seeds.north <- seeds %>% 
  filter(site=="Coast Fork of Williamette"|site=="Rock Creek"|site=="Canton Creek") %>% 
  droplevels()

# these are Matt's transplant data
# need them here to get site name vector
# also, they are used by source file 'integral_projection_model.R' below
matt_dat <- read_csv("Data/IPMData_transplant.csv")

# log-transform size
matt_dat <- matt_dat %>% mutate(logSize = log(z), logSizeNext = log(z1), Fec0 = Repr, Fec1 = Fec)

# make sure plot and site are recognized as factors
matt_dat$PlotID = as.factor(matt_dat$PlotID)
matt_dat$SiteID = as.factor(matt_dat$SiteID)

#*******************************************************************************
#### 2. Read in global survival, growth and fecundity models using data from all sites ###
#*******************************************************************************

# Create a vector of unique Site names for subsetting (n=8)
site=unique(matt_dat$SiteID)

# Set up data frame of model parameters
params=c()

#*******************************************************************************
  ### 3A. Survival ###
  #*******************************************************************************

  # Read in top survival model output (Formula: Surv ~ logSize + Site + (1 | Plot))
  surv.reg=load("Robjects/surv.reg.rda")

  # Get model coefficients
  fixef(s6)
  
  # Store model coefficients
  params$site=site
  params$surv.globint=fixef(s6)[1] 
  params$surv.siteint=c(0,fixef(s6)[3:9])
  params$surv.slope=fixef(s6)[2]
  params$surv.int = params$surv.globint + params$surv.siteint
    
  #*******************************************************************************
  ### 3B. Growth ###
  #*******************************************************************************
  
  # Read in top growth model output (Formula: logSizeNext ~ logSize + Site + (1 | Plot)
  growth.reg=load("Robjects/growth.reg.rda")
  
  # Get model coefficients
  fixef(g6)
  
  # Store model coefficients
  params$growth.globint=fixef(g6)[1] 
  params$growth.siteint=c(0,fixef(g6)[3:9])
  params$growth.slope=fixef(g6)[2] 
  params$growth.sd=rep(sigma(g6),times=length(site)) 
  params$growth.int = params$growth.globint + params$growth.siteint
  
  #*******************************************************************************
  ### 3C. Flowering ###
  #*******************************************************************************
  
  # Read in top flowering model output (Formula: Fec0 ~ logSize + Site + (1 | Plot))
  flowering.reg=load("Robjects/flowering.reg.rda")

  # Store model coefficients
  params$flowering.globint=fixef(fl6)[1] 
  params$flowering.siteint=c(0,fixef(fl6)[3:9])
  params$flowering.slope=fixef(fl6)[2] 
  params$flowering.int = params$flowering.globint + params$flowering.siteint
  
  #*******************************************************************************
  ### 3D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  
  # Read in top model output for fruit.reg (Formula: Fec0 ~ logSize + Site + (1 | Plot))   
  fruit.reg=load("Robjects/fruit.reg.rda")

  # Store model coefficients
  params$fruits.globint=fixef(fr9)[1] 
  params$fruits.siteint=c(0,fixef(fr9)[3:9])
  params$fruits.slope=fixef(fr9)[2] 
  params$fruits.int = params$fruits.globint + params$fruits.siteint
  
  #*******************************************************************************
  ### 3E. Create data frame of site-specific parameter estimates ###
  #*******************************************************************************
  
  params=data.frame(params)

  #*******************************************************************************
  ### 3F. Number of seeds per fruit (constant across sites) ###
  #*******************************************************************************
  
  seeds.north.mean <- seeds.north %>% 
    summarize(seeds.per.fruit = mean(newgrandmean)) 
  
  params$seeds.per.fruit = seeds.north.mean$seeds.per.fruit

  #*******************************************************************************
  ### 3G. Establishment probability (constant across sites) ###
  #*******************************************************************************
  
  # Obtain number of new recruits per 3 northern sites
  recruit.number.per.site=tapply(demo.dat.north$logSizeNext[is.na(demo.dat.north$logSize)],demo.dat.north$Site[is.na(demo.dat.north$logSize)],FUN="length") %>% data.frame()
  colnames(recruit.number.per.site)="recruit.number"
  
  # Obtain total fruit count per 3 northern sites 
  fruits.per.site=tapply(demo.dat.north$Fec1[!is.na(demo.dat.north$Fec1)],demo.dat.north$Site[!is.na(demo.dat.north$Fec1)],sum)

  # Obtain total seed count per site (= # fruits per site * # seeds per fruit per site)
  total.seeds.per.site=fruits.per.site*seeds.north$newgrandmean	
  
  # Estimate establishment probability as # of new recruits/# of seeds
  params$establishment.prob=mean(recruit.number.per.site$recruit.number/total.seeds.per.site)
  
  #*******************************************************************************
  ### 3H. Size distribution of recruits (constant across sites) ###
  #*******************************************************************************
  
  recruit.size.mean=mean(demo.dat.north$logSizeNext[is.na(demo.dat.north$logSize)])
  recruit.size.sd=sd(demo.dat.north$logSizeNext[is.na(demo.dat.north$logSize)])
  
  params$recruit.logSize.mean=recruit.size.mean  
  params$recruit.logSize.sd=recruit.size.sd  
  

#### Store parameters in .csv file for later use
write.csv(params,"Robjects/vital_rate_coefficients.csv",row.names=FALSE)
  
#*******************************************************************************
### 4. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
#*******************************************************************************

  #*******************************************************************************
  ### 4A. Subset data for site f
  #*******************************************************************************
  
  # create empty vectors for lambda and site to be filled
  lambda=c()
  Site=character()
  
  for (f in 1:length(site)) {
    data1=subset(matt_dat,SiteID==site[f])
    params1=subset(params,site==site[f])
    params1=subset(params1,select=-c(site,surv.globint, surv.siteint, growth.globint, growth.siteint, flowering.globint, flowering.siteint, fruits.globint, fruits.siteint))
    
    #*******************************************************************************
    ### 4B. Create survival, growth, and fecundity functions and build IPM by running integral_projection_model.R script
    #*******************************************************************************
    
    source("Rcode/integral_projection_model.R")
    
    #*******************************************************************************
    ### 4C. Obtain lambda estimate for site f
    #*******************************************************************************
    
    lambda[f] <- Re(eigen(K)$values[1])
    Site[f]=as.character(site[f])
    } # end loop to run IPMs and estimate lambdas for each site
    
    # make data frame of site and lambda
    site.lambda=data.frame(Site,lambda)
    
#*******************************************************************************
### 5. Merge site information with lambda estimates and save to .csv file
#*******************************************************************************

# Read in site info
site.info=read_csv("Data/raw_data/WunderGround/sites.csv")
    
# merge site info with lambda estimates
site.lambda=left_join(site.lambda, site.info, by=c("Site" = "ID2"))

plot(lambda ~ lat, data=site.lambda) # looks qualitatively the same as in Matt's thesis

# save to .csv file 
write.csv(site.lambda,"Robjects/site.lambda.csv",row.names=FALSE)

