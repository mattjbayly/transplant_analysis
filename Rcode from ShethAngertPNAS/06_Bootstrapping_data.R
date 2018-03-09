#### PROJECT:  Mimulus cardinalis northern transplant 2015-2016
#### PURPOSE: Sample unique individuals from each site with replacement to create bootstrap datasets
############# From each bootstrapped dataset, vital rate models are created and IPMs are run to obtain bootstrapped lambdas
############# Replicate bootstrap datasets will be used to obtain confidence intervals around lambda estimates for each site
#### AUTHOR: Seema Sheth
#### DATE LAST MODIFIED: 20180305

# remove objects and clear workspace
rm(list = ls(all=TRUE))

# require packages
require(plyr)
require(dplyr)
require(tidyverse)

# set working directory
# setwd("/Users/ssheth/Google Drive/demography_PNAS_November2017")

#*******************************************************************************
#### 1. bring in Matt's cardinalis transplant data and demography data ###
#*******************************************************************************

matt_dat <- read_csv("Data/IPMData_transplant.csv")

# log-transform size and rename reproductive variables for consistency with Seema
matt_dat <- matt_dat %>% mutate(logSize = log(z), logSizeNext = log(z1), Fec0 = Repr, Fec1 = Fec)

# make sure plot and site are recognized as factors
matt_dat$PlotID = as.factor(matt_dat$PlotID)
matt_dat$SiteID = as.factor(matt_dat$SiteID)
matt_dat$Year = as.factor(matt_dat$Year)

# individual-level data
# run data prep script to clean up and import data
source("Rcode/data_prep.R")
# head(demo.data)
# head(site_fruit_count_data)

# filter to northern sites only
demo.dat.north <- demo.data %>% 
  filter(Site=="Rock Creek"|Site=="Canton Creek"|Site=="Coast Fork of Williamette") %>% 
  droplevels()
site_fruit_count_north <- site_fruit_count_data %>% 
  filter(Site=="Rock Creek"|Site=="Canton Creek"|Site=="Coast Fork of Williamette") %>% 
  droplevels()

# get new recruits from demo.dat.north
recruits <- demo.dat.north %>% filter(is.na(logSize))
 
# seeds per fruit
# don't have individual-level fruit data
# instead, draw from normal distribution with mean and sd of yearly averages
seeds <- read_csv("Data/Amy_wild_demo_data/fec2.seed.per.fruit.2010.2011.2012.csv")

# filter to northern sites only
# note: seed counts missing for Coast Fork (substituted Canton instead) so Coast Fork not included here
seeds.north <- seeds %>% 
  filter(site=="Rock Creek"|site=="Canton Creek") %>% 
  droplevels()

# note: 2011 and 2012 estimates missing for Canton (substituted 2010 instead) so set to NA here
seeds.north[seeds.north$site=="Canton Creek",]$mean2011 = NA
seeds.north[seeds.north$site=="Canton Creek",]$mean2012 = NA


#*******************************************************************************
#### 2. Create nested loop to obtain replicate bootstrap datasets for each of matt's sites, sampling with replacement ###
#*******************************************************************************

# Obtain a data frame of unique IDs from each Site
ID.by.Site=unique(matt_dat[,c("SiteID", "ID")])
# n=573 individuals

# Create a vector of unique Site names for subsetting; note this is sorted by decreasing latitude 
site=unique(matt_dat$SiteID)
# n=8 sites

# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set seed for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(site)) {
  data.site=subset(matt_dat,SiteID==site[i]) # select data from site i
  id.site=subset(ID.by.Site,SiteID==site[i]) # select list of unique individual IDs from site i
  id.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(id.site,size=nrow(id.site), replace = T)}) %>% ldply() # resample rows of site i's data with replacement and size=number of unique individuals in original dataset for each site and convert list to data frame
  
  id.boot$Replicate=rep(seq(1:n.boot),each=nrow(id.site)) # create a column in data frame that corresponds to bootstrap replicate
  data.boot=join(id.boot,data.site,type="left",match="all") # merge bootstrapped list of unique IDs to full dataset
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}

# Convert list to data frame
bootstrapped.data <- do.call(rbind, data.boot.rep) 

# Write bootstrapped datasets to .rds file
saveRDS(bootstrapped.data,"Robjects/Mcard_transplant_INDIV_BOOTSTRAP_data.rds")  


#*******************************************************************************
#### 3. Bootstrap the demography data used for fecundity estimates, sampling with replacement from three northern sites ###
#*******************************************************************************

### Required components:
# Number of seeds per fruit
  # draw from normal distribution
# Probability of recruitment 
  # for denominator, sample individuals with replacement from site_fruit_count_north
  # for numerator, sample individuals from demo.dat.north and tally recruits
# Size distribution of recruits
  # draw recruits with replacement and calculate size distribution

### Number of seeds per fruit
seeds.dist <- seeds.north %>% 
  summarize(seeds.mean = mean(c(mean2010, mean2011, mean2012), na.rm=T),
         seeds.sd = sd(c(mean2010, mean2011, mean2012), na.rm=T))
# but this mean is lower than that used for real estimates because of removing pseudoreplicated values for coast fork
# it could be causing bootstrapped lambdas to be consistently lower than real estimates
# use empirical sd, but same mean as used for real estimates

# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set seed for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
data.boot <- lapply(1:n.boot, function(j) { 
  set.seed(j+seed)
  rnorm(1, mean=1163.28, sd=seeds.dist$seeds.sd)}) %>% ldply() # 
  data.boot$Replicate=rep(seq(1:n.boot)) # create a column in data frame that corresponds to bootstrap replicate

# rename
bootstrapped.seeds <- data.boot
colnames(bootstrapped.seeds) <- c("seeds", "Replicate")
  
# Write bootstrapped datasets to .rds file
saveRDS(data.boot,"Robjects/Mcard_transplant_SEEDS_BOOTSTRAP_data.rds")  


### Probability of recruitment

## Total fruits produced at time t

# Create a vector of unique Site names for subsetting 
site=unique(site_fruit_count_north$Site)

# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set seed for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(site)) {
  data.site=subset(site_fruit_count_north,Site==site[i]) # select data from site i
  data.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(data.site,size=nrow(data.site), replace = T)}) %>% ldply() # resample rows of fruiting adults with replacement and size=original dataset for each site and convert list to data frame
  data.boot$Replicate=rep(seq(1:n.boot)) # create a column in data frame that corresponds to bootstrap replicate
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}

# Convert list to data frame
bootstrapped.fruits <- do.call(rbind, data.boot.rep) 

bootstrapped.fruit.sum <- bootstrapped.fruits %>% 
  group_by(Site, Replicate) %>% 
  summarize(total.fruits = sum(Fec1, na.rm=T)) %>% 
  ungroup

# Write bootstrapped datasets to .rds file
saveRDS(bootstrapped.fruits,"Robjects/Mcard_transplant_FRUITS_BOOTSTRAP_data.rds")  
saveRDS(bootstrapped.fruit.sum,"Robjects/Mcard_transplant_FRUITS_BOOTSTRAP_sum.rds") 

## Total recruits observed at time t+1

# Create a vector of unique Site names for subsetting 
# don't need to track IDs because each unique individual can only recruit once
site=unique(demo.dat.north$Site)
            
# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set seed for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(site)) {
  data.site=subset(demo.dat.north,Site==site[i]) # select data from site i
  data.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(data.site,size=nrow(data.site), replace = T)}) %>% ldply() # resample rows of recruits with replacement and size=original dataset for each site and convert list to data frame
  data.boot$Replicate=rep(seq(1:n.boot)) # create a column in data frame that corresponds to bootstrap replicate
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}

# Convert list to data frame
bootstrapped.indivs <- do.call(rbind, data.boot.rep) 

# Ok, so how many of these are recruits (as opposed to survivors)?
bootstrapped.recruit.num <- bootstrapped.indivs %>% 
  filter(is.na(logSize)) %>% 
  group_by(Site, Replicate) %>% 
  summarize(recruit.num = n()) %>% 
  ungroup

# Write bootstrapped datasets to .rds file
saveRDS(bootstrapped.indivs,"Robjects/Mcard_transplant_DEMOGINDIVS_BOOTSTRAP_data.rds")  
saveRDS(bootstrapped.recruit.num,"Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_num.rds")  

bootstrapped.indivs <- readRDS("Robjects/Mcard_transplant_DEMOGINDIVS_BOOTSTRAP_data.rds")
bootstrapped.recruit.num <- readRDS("Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_num.rds")
bootstrapped.seeds <- readRDS("Robjects/Mcard_transplant_SEEDS_BOOTSTRAP_data.rds")
colnames(bootstrapped.seeds) = c("seeds", "Replicate")

### Combine above into recruitment probability
df <- left_join(bootstrapped.recruit.num, bootstrapped.fruit.sum)
df2 <- left_join(df, bootstrapped.seeds)

bootstrapped.recruit.prob <- df2 %>% 
  mutate(recruit.prob = recruit.num/(total.fruits*seeds)) %>%
  group_by(Replicate) %>% 
  summarize(mean.recruit.prob = mean(recruit.prob))

saveRDS(bootstrapped.recruit.prob,"Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_prob.rds")  

### Size distribution of recruits

# Create a vector of unique Site names for subsetting 
# don't need to track IDs because each unique individual can only recruit once
site=unique(recruits$Site)

# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set seed for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(site)) {
  data.site=subset(recruits,Site==site[i]) # select data from site i
  id.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(data.site,size=nrow(data.site), replace = T)}) %>% ldply() # resample rows of recruits with replacement and size=original dataset for each site and convert list to data frame
  id.boot$Replicate=rep(seq(1:n.boot)) # create a column in data frame that corresponds to bootstrap replicate
  data.boot=join(id.boot,data.site,type="left",match="all") # merge bootstrapped list of unique IDs to full dataset
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}

# Convert list to data frame
bootstrapped.recruits <- do.call(rbind, data.boot.rep) 

bootstrapped.recruit.dist <- bootstrapped.recruits %>% 
  group_by(Replicate) %>% 
  summarize(recruit.size.mean = mean(logSizeNext),
            recruit.size.sd = sd(logSizeNext))

# Write bootstrapped datasets to .rds file
saveRDS(bootstrapped.recruits,"Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_data.rds")  
saveRDS(bootstrapped.recruit.dist,"Robjects/Mcard_transplant_RECRUITS_BOOTSTRAP_dist.rds")  
