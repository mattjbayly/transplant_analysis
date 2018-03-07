#### PROJECT: Mimulus cardinalis northern transplant
#### PURPOSE: From each bootstrapped dataset, vital rate models are created and IPMs are run to obtain bootstrapped lambdas
#### AUTHOR: Seema Sheth
#### DATE LAST MODIFIED: 20180306

# remove objects and clear workspace
rm(list = ls(all=TRUE))

# require packages
require(plyr)
require(dplyr)
require(lme4)
require(glmmADMB)

# set working directory
# setwd("/Users/ssheth/Google Drive/demography_PNAS_November2017")

#*******************************************************************************
#### 1. Preliminaries ###
#*******************************************************************************

## Read in & examine bootstrapped data 
# Matt's individual data
bootstrapped.data=readRDS("Robjects/Mcard_transplant_INDIV_BOOTSTRAP_data.rds")
# str(bootstrapped.data)
bootstrapped.data$SiteID <- as.factor(bootstrapped.data$SiteID)
bootstrapped.data$PlotID <- as.factor(bootstrapped.data$PlotID)
bootstrapped.data$Year <- as.factor(bootstrapped.data$Year)

## Read in & examine bootstrapped coefficients from vital rate models
surv.reg_boot=readRDS("/Users/amyangert/Documents/GitClones/matt transplant bootstraps/surv.reg_boot.rds")
growth.reg_boot=readRDS("/Users/amyangert/Documents/GitClones/matt transplant bootstraps/growth.reg_boot.rds")
flowering.reg_boot=readRDS("/Users/amyangert/Documents/GitClones/matt transplant bootstraps/flowering.reg_boot.rds")
fruit.reg_boot=readRDS("/Users/amyangert/Documents/GitClones/matt transplant bootstraps/fruit.reg_boot.rds")

## Read in and examine other bootstrapped parameters (constants across sites)
# seeds per fruit
bootstrapped.seed.num=readRDS("Robjects/Mcard_transplant_SEEDS_BOOTSTRAP_data.rds")

# recruitment probabilty
bootstrapped.recruit.prob=as.data.frame(readRDS("Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_prob.rds"))
# add Replicate column
bootstrapped.recruit.prob$Replicate = seq(1:2000)
colnames(bootstrapped.recruit.prob) = c("recruit.prob", "Replicate")

# recruit size distribution
bootstrapped.recruit.dist=readRDS("Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_dist.rds")


## Create a vector of unique Site names for subsetting 
site=unique(bootstrapped.data$SiteID)

## Set number of bootstrap replicate datasets
n.boot=2000

#*******************************************************************************
#### 2. Obtain vital rate parameters across all sites for each replicate bootstrap dataset ###
#*******************************************************************************

# Create empty list to be filled in loop
params.boot=list()

# Begin loop to extract vital rate parameters for each bootstrapped dataset
for (k in 1:n.boot) {
  data.rep=subset(bootstrapped.data,Replicate==k) # select data from replicate k
  
  # Obtain total fruit count for each indivdiual at each site in each year, including monster plants
  #site_fruit_count_data=subset(data.rep,select=c(SiteID,Fec1)) 
  
  # Set up data frame of model parameters
  params=c()
  
  #*******************************************************************************
  ### 3A. Survival ###
  #*******************************************************************************
  
  # Read in model coefficients
  surv.reg=surv.reg_boot[[k]]
   
  # Store model coefficients
  params$surv.globint=fixef(surv.reg)[1]
  params$surv.siteint=c(0,fixef(surv.reg)[3:9]) 
  params$surv.slope=fixef(surv.reg)[2] 
  params$surv.int=params$surv.globint+params$surv.siteint
  params$Site=site
  
  # Make into data frame
  params=data.frame(params)
  
  #*******************************************************************************
  ### 3B. Growth ###
  #*******************************************************************************
  
  # Read in model coefficients
  growth.reg=growth.reg_boot[[k]]
  
  # Store model coefficients
  growth.params=c()
  growth.params$growth.globint=fixef(growth.reg)[1]  
  growth.params$growth.siteint=c(0,fixef(growth.reg)[3:9])
  growth.params$growth.slope=fixef(growth.reg)[2]
  growth.params$growth.int=growth.params$growth.globint+growth.params$growth.siteint
  growth.params$Site=site
  growth.params=data.frame(growth.params)
  params=join(params,growth.params,by="Site")
  params$growth.sd=rep(sigma(growth.reg),times=length(site)) 
  
  #*******************************************************************************
  ### 3C. Flowering ###
  #*******************************************************************************
  
  # Read in model coefficients
  flowering.reg=flowering.reg_boot[[k]]
  
  # Store model coefficients
  params$flowering.globint=fixef(flowering.reg)[1] 
  params$flowering.siteint=c(0,fixef(flowering.reg)[3:9]) 
  params$flowering.slope=fixef(flowering.reg)[2] 
  params$flowering.int=params$flowering.globint+params$flowering.siteint
  
  #*******************************************************************************
  ### 3D. Fruit number (untransformed) using negative binomial regression ###
  #*******************************************************************************
  
  # Read in model coefficients
  fruit.reg=fruit.reg_boot[[k]]
  
  # Store model coefficients
  params$fruits.globint=fixef(fruit.reg)[1] 
  params$fruits.siteint=c(0,fixef(fruit.reg)[3:9]) 
  params$fruits.int=params$fruits.globint+params$fruits.siteint
  params$fruits.slope=fixef(fruit.reg)[2] 
  
  #*******************************************************************************
  ### 3E. Size distribution of recruits ###
  #*******************************************************************************
 
  params$recruit.logSize.mean=unlist(rep(bootstrapped.recruit.dist[k,"recruit.size.mean"],times=length(site)))  
  params$recruit.logSize.sd=unlist(rep(bootstrapped.recruit.dist[k,"recruit.size.sd"],times=length(site)))  
  
  #*******************************************************************************
  ### 3F. Number of seeds per fruit ###
  #*******************************************************************************
  
  params$seeds.per.fruit=rep(bootstrapped.seed.num[k,"V1"],times=length(site))
  
  #*******************************************************************************
  ### 3G. Establishment probability ###
  #*******************************************************************************
  
  params$establishment.prob=rep(bootstrapped.recruit.prob[k,"bootstrapped.recruit.prob"],times=length(site))
  
  
  
  
  # Store data frame of parameter values for a given bootstrap replicate dataset into list
  params$Replicate=rep(k,each=nrow(params)) # create a column in data frame that corresponds to bootstrap replicate
  params.boot[[k]]=params
  print(k)
  } # end loop

# remove large objects
rm(surv.reg_boot)
rm(growth.reg_boot)
rm(flowering.reg_boot)
rm(fruit.reg_boot)
  
# Convert list of bootstrapped vital rate parameters to data frame
bootstrapped.params <- do.call(rbind, params.boot)
  
# Write bootstrapped parameter estimates to .csv file
write.csv(bootstrapped.params,"Robjects/Mcard_transplant_BOOTSTRAP_params.csv",row.names=FALSE)  
  
  #*******************************************************************************
  ### 4. Create site-specific IPMs parameterized by site-specific parameters derived from global vital rates models 
  #*******************************************************************************

# Create empty list to be filled in loop
lambda.boot=list()

# Begin loop to extract vital rate parameters for each bootstrapped dataset
for (k in 1:n.boot) {
  data.rep=subset(bootstrapped.data,Replicate==k) # select data from replicate k
  params=subset(bootstrapped.params,Replicate==k) # select parameters from replicate K
  
# Remove monster plants where individuals were not distinguished; do this again when estimating lambda to restrict size range
#### NOTE: these are plants that A. Angert noted as "not ok, definitely exclude from survival, growth, and fecundity but ok for seed input denominator for recruitment (history of lumping/splitting/relumping; redundant IDs)"
data.rep=subset(data.rep,NotAnIndividual!=1|is.na(NotAnIndividual))
  
  #*******************************************************************************
  ### 4A. Set up loop to create site-specific vital rate functions and IPMs
  #*******************************************************************************
  
  # create empty vector for lambda and siteID to be filled
  lambda=c()
  siteID=c()
  
  for (j in 1:length(site)) {
    data1=subset(data.rep,SiteID==site[j])
    params1=subset(params,Site==site[j])
    params1=subset(params1,selectID=-Site)
    
    # write if else statement so that if all individuals in bootstrap died, lambda is manually set equal to 0
    
    if(length(data1$logSizeNext[!is.na(data1$logSizeNext)])>0)
      
    {
      
      #*******************************************************************************
      ### 4B. Create survival, growth, and fecundity functions and build IPM by running integral_projection_model.R script
      #*******************************************************************************
      
      source("Rcode/integral_projection_model.R")
      
      #*******************************************************************************
      ### 4C. Obtain lambda
      #*******************************************************************************
      
      # obtain lambda estimate
      (lambda[j] <- Re(eigen(K)$values[1])) 
      
      # obtain site name corresponding to lambda estimate
      siteID[j]=as.character(unique(data1$Site))
      
    } # end if loop
    else{
      lambda[j]=0     
      siteID[j]=as.character(unique(data1$Site))} # end else loop
    
    # merge lambda estimate with site name
    lambda.site=data.frame(siteID,lambda)
    print(siteID[j]) # print site ID to keep track of what is currently running
  } # end among-site loop
  
  lambda.boot[[k]]=lambda.site
  print(k) # print replicate number to keep track of how many bootstraps are complete
} # end among-replicate loop

#*******************************************************************************
### 5. Write bootstrapped lambda and vital rate parameter estimates to .csv file
#*******************************************************************************

# Convert list of bootstrapped lambdas to data frame and arrange by Site
bootstrapped.lambda <- do.call(rbind, lambda.boot) %>% arrange(siteID)
bootstrapped.lambda$Replicate=rep(seq(1,n.boot,1),times=length(site))

# write bootstrapped lambda estimates to .csv
write.csv(bootstrapped.lambda,"R_output/Mcard_demog_INDIV_BOOTSTRAP_lambda_2010-2013.csv",row.names=FALSE) 