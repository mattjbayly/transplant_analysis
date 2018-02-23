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

# (deleted 1. import and format data; we are using saved models)

#*******************************************************************************
#### 2. Create global survival, growth and fecundity models using data from all sites ###
#*******************************************************************************

# Create a vector of unique Site names for subsetting (n=8)
site=unique(data$SiteID)

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
  
  #*******************************************************************************
  ### 3C. Flowering ###
  #*******************************************************************************
  
  # Read in top flowering model output (Formula: Fec0 ~ logSize + Site + (1 | Plot))
  flowering.reg=load("Robjects/flowering.reg.rda")

  # Store model coefficients
  params$flowering.globint=fixef(fl6)[1] 
  params$flowering.siteint=c(0,fixef(fl6)[3:9])
  params$flowering.slope=fixef(fl6)[2] 
  
  #*******************************************************************************
  ### 3D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  
  # Read in top model output for fruit.reg (Formula: Fec0 ~ logSize + Site + (1 | Plot))   
  fruit.reg=load("Robjects/fruit.reg.rda")

  # Store model coefficients
  params$fruits.globint=fixef(fr9)[1] 
  params$fruits.siteint=c(0,fixef(fr9)[3:9])
  params$fruits.slope=fixef(fr9)[2] 
  
  #*******************************************************************************
  ### 3E. Create data frame of site-specific parameter estimates ###
  #*******************************************************************************
  
  params=data.frame(params)

  #*******************************************************************************
  ### 3F. Number of seeds per fruit (constant across sites) ###
  #*******************************************************************************
  
  # from demography data
  seeds <- read_csv("Data/Amy_wild_demo_data/fec2.seed.per.fruit.2010.2011.2012.csv")
  
  seeds.sources <- seeds %>% 
    filter(site=="Coast Fork of Williamette"|site=="Rock Creek"|site=="Canton Creek") %>% 
    summarize(seeds.per.fruit = mean(newgrandmean))  
  
  params$seeds.per.fruit = seeds.sources$seeds.per.fruit

  #*******************************************************************************
  ### 3G. Establishment probability (constant across sites) ###
  #*******************************************************************************
  
  # Obtain number of new recruits per site
  recruit.number=tapply(data$logSizeNext[is.na(data$logSize)],data$Site[is.na(data$logSize)],FUN="length") %>% data.frame()
  colnames(recruit.number)="recruit.number"
  
  # Obtain total fruit count per site 
  fruits.per.site=tapply(site_fruit_count_data$Fec1[!is.na(site_fruit_count_data$Fec1)],site_fruit_count_data$Site[!is.na(site_fruit_count_data$Fec1)],sum)
  
  # Obtain total seed count per site (= # fruits per site * # seeds per fruit per site)
  total.seeds.per.site=fruits.per.site*seeds.per.site$seed.ct	
  
  # Estimate establishment probability as # of new recruits/# of seeds
  params$establishment.prob=recruit.number$recruit.number/total.seeds.per.site
  
  # Set establishment probability as 0 for Hauser Creek (was calculated as NA because Hauser creek has 0 new recruits)
  params$establishment.prob[is.na(params$establishment.prob)]=0	

  #*******************************************************************************
  ### 3H. Size distribution of recruits (constant across sites) ###
  #*******************************************************************************
  
  # from the demography dataset 
  Rec_dist <- read_csv("Robjects/lnormFecKern.csv")
  
  N_Rec_dist <- Rec_dist %>% 
    filter(Reg=="N") %>% 
    summarize(recruit.size.mean=mean(meanlog),
              recruit.size.sd=mean(sdlog))
  
  params$recruit.logSize.mean=N_Rec_dist$recruit.size.mean
  params$recruit.logSize.sd=N_Rec_dist$recruit.size.sd


#### Store parameters in .csv file for later use
write.csv(params,"R_output/vital_rate_coefficients.csv",row.names=FALSE)
  
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
    data1=subset(data,Site==site[f])
    params1=subset(params,Site==site[f])
    params1=subset(params1,select=-Site)
    
    #*******************************************************************************
    ### 4B. Create survival, growth, and fecundity functions and build IPM by running integral_projection_model.R script
    #*******************************************************************************
    
    source("R_scripts/integral_projection_model.R")
    
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

# Create data frame of Site, Latitude, Longitude, Region, and Elevation for hypothesis testing
site.info=subset(data,select=c(Site,Latitude,Longitude,Elevation,Region,RegionRank)) %>% unique() %>% arrange(-Latitude)
    
# merge site info with lambda estimates
site.info=join(site.info,site.lambda)

# save to .csv file 
write.csv(site.info,"R_output/site.lambda.csv",row.names=FALSE)

