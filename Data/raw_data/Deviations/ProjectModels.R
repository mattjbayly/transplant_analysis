###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   Read in saved model objects and thresholds 
##   Project models to proposed transplant sites
## last update:  2 June 2014
###################################################################################

# IMPORTANT: prior to running this script you must first check the file 'sites' 
# to make sure your selected sites have the right climate data & variables used in the model 
# In this files 'sites.csv' herbarium records have year specific climate avage data, but
# demography sites, transplant sites & Seema's sites all have climate data from 1980 - 2010. 
# 'sites_pred.csv' ~ final file with site predictions 

# load libraries now if desired; otherwise loaded below as needed
library(raster)
library(maptools)
library(gam)
library(randomForest)
library(dismo)
library(gbm)

#http://gis.stackexchange.com/questions/63577/joining-polygons
#http://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon

## set pathnames - Amy
path.root="/Users/amyangert/Google Drive/OccAmNat" 
path.dat = paste(path.root, "/data files", sep="")
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.dat, "/ecoregions.shp", sep="")
path.bio = paste(path.dat, "/wc0.5", sep="")
path.sta = paste(path.dat, "/gz_2010_us_040_00_500k", sep="")

## set pathnames - Matthew
path.root = "C:/Users/DW/Desktop/temp.sept.30" 
path.dat = paste(path.root, "/data files", sep="")
path.dat.fix = paste(path.root, "/data files", sep="") # older files relocated to another directory
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.obj, "/ecoregions.shp", sep="")
path.bio = paste(path.obj, "/wc0.5", sep="")
path.cod=paste(path.root, "/R code", sep="")
path.fig=paste(path.root, "/figures", sep="")
path.sta = paste(path.obj, "/gz_2010_us_040_00_500k", sep="")
path.gis=paste(path.root, "/path.gis", sep="")

################################################################################
######## START INITIALIZATION OF DATA STRUCTURES



## load SDM models 
setwd(path.obj)
for (i in 1:10) {
	mod.lr = get(load(paste("LR.mod2.",i,".pseudo11.Rda", sep="")))
	assign(paste("LR.mod2.",i, sep=""), mod.lr)
	mod.gam = get(load(paste("GAM.mod4.",i,".pseudo11.Rda", sep="")))
	assign(paste("GAM.mod4.",i, sep=""), mod.gam)
	mod.rf = get(load(paste("RF.mod1.",i,".pseudo11.Rda", sep="")))  
	assign(paste("RF.mod1.",i, sep=""), mod.rf)
	mod.brt = get(load(paste("BRT.mod4.",i,".pseudo11.Rda", sep="")))
	assign(paste("BRT.mod4.",i, sep=""), mod.brt)  
	mod.max = get(load(paste("MAX.mod1.",i,".pseudo11.Rda", sep="")))
	assign(paste("MAX.mod1.",i, sep=""), mod.max)  
	}

## read in replicate training datasets
library(sp)
library(rgdal)
setwd(path.dat) 
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
for (i in 1:10) {
	dat = read.csv(paste("dat",i,"c.csv", sep=""))
	dat = na.omit(dat)
	coordinates(dat) = ~Longitude + Latitude
	proj4string(dat) = CRS(prj.wgs)
	assign(paste("dat",i, sep=""), dat)
	dat.aea = spTransform(dat, CRS=CRS(prj.aea))
	dat.lcc = spTransform(dat, CRS=CRS(prj.lcc))
	assign(paste("dat",i,".aea", sep=""), dat.aea)
	assign(paste("dat",i,".lcc", sep=""), dat.lcc)
	}

## read in occupancy dataset
library(sp)
library(rgdal)
setwd(path.dat)
all = read.csv("all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
occ = all[all$DATASET=="occ",] #pull out occupancy dataset
dim(occ) #check size
occ$bio3 = log(occ$bio3+0.5) #make needed ln-transforms of predictors
occ$bio10 = log(occ$bio10+0.5)
occ$bio12 = log(occ$bio12+0.5)
occ$bio14 = log(occ$bio14+0.5)
coordinates(occ) = ~Longitude + Latitude
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
proj4string(occ) = CRS(prj.wgs)
occ.aea = spTransform(occ, CRS=CRS(prj.aea))
occ.lcc = spTransform(occ, CRS=CRS(prj.lcc))

## split by elevation 
L.thresh = 400
H.thresh = 1200
high <- occ[occ$Elevation > H.thresh,]
low <- occ[occ$Elevation < L.thresh,]
mid <- occ[occ$Elevation >= L.thresh & occ$Elevation <= H.thresh,]
high.aea = occ.aea[occ.aea$Elevation > H.thresh,]
low.aea <- occ.aea[occ.aea$Elevation < L.thresh,]
mid.aea <- occ.aea[occ.aea$Elevation >= L.thresh & occ.aea$Elevation <= H.thresh,]

high.lcc = occ.lcc[occ.lcc$Elevation > H.thresh,]
low.lcc <- occ.lcc[occ.lcc$Elevation < L.thresh,]
mid.lcc <- occ.lcc[occ.lcc$Elevation >= L.thresh & occ.lcc$Elevation <= H.thresh,]


## read in occupancy dataset split by ecoregions
setwd(path.dat)
north = read.csv("north.csv")
center = read.csv("center.csv")
south = read.csv("south.csv")
prj.wgs = "+proj=longlat +ellps=WGS84"

coordinates(north) = ~Longitude + Latitude
proj4string(north) = CRS(prj.wgs)
coordinates(center) = ~Longitude + Latitude
proj4string(center) = CRS(prj.wgs)
coordinates(south) = ~Longitude + Latitude
proj4string(south) = CRS(prj.wgs)

## reproject to albers equal area
library(rgdal)
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
occ.aea = spTransform(occ, CRS=CRS(prj.aea))
north.aea = spTransform(north, CRS=CRS(prj.aea))
center.aea = spTransform(center, CRS=CRS(prj.aea))
south.aea = spTransform(south, CRS=CRS(prj.aea))
north.lcc = spTransform(north, CRS=CRS(prj.lcc))
center.lcc = spTransform(center, CRS=CRS(prj.lcc))
south.lcc = spTransform(south, CRS=CRS(prj.lcc))

# get model cut points 
setwd(path.obj)
ext.accs.lr = get(load("LR.mod2.extaccs.pseudo11.Rda"))
ext.accs.gam = get(load("GAM.mod4.extaccs.pseudo11.Rda"))
ext.accs.rf = get(load("RF.mod1.extaccs.pseudo11.Rda"))
ext.accs.brt = get(load("BRT.mod4.extaccs.pseudo11.Rda"))
ext.accs.max = get(load("MAX.mod1.extaccs.pseudo11.Rda"))
lr.cuts = ext.accs.lr[ext.accs.lr$thresh=="SensSpec", "threshold"]
gam.cuts = ext.accs.gam[ext.accs.gam$thresh=="SensSpec", "threshold"]
rf.cuts = ext.accs.rf[ext.accs.rf$thresh=="SensSpec", "threshold"]
brt.cuts = ext.accs.brt[ext.accs.brt$thresh=="SensSpec", "threshold"]
max.cuts = ext.accs.max[ext.accs.max$thresh=="SensSpec", "threshold"]
cuts <- cbind(lr.cuts, gam.cuts, rf.cuts, brt.cuts, max.cuts)
cuts.avg <- colMeans(cuts)

#########################################################

### this seems to make predictions to sites based on short-term climate data
path.dev <- "C:/Users/DW/Desktop/transplant_analysis/Data/raw_data/Deviations"
# get sites climate data
setwd(path.dev)
sites <- read.csv(file="trans2015.csv")
	library(dismo)
	tmax <- sites[,c("Tmax01", "Tmax02", "Tmax03",    "Tmax04",   "Tmax05",    "Tmax06",    "Tmax07",    "Tmax08", "Tmax09" ,   "Tmax10"  ,  "Tmax11" ,   "Tmax12")]
	tmin <- sites[,c("Tmin01", "Tmin02", "Tmin03",    "Tmin04",   "Tmin05",    "Tmin06",    "Tmin07",    "Tmin08", "Tmin09" ,   "Tmin10"  ,  "Tmin11" ,   "Tmin12")]
	prec <- sites[,c("PPT01", "PPT02", "PPT03",    "PPT04",   "PPT05",    "PPT06",    "PPT07",    "PPT08", "PPT09" ,   "PPT10"  ,  "PPT11" ,   "PPT12")]
	tmax <- as.matrix(tmax)
	tmin <- as.matrix(tmin)
	prec <- as.matrix(prec)

	b <- biovars(prec, tmin, tmax)
	b <- data.frame(b)
	sites <- cbind(sites, b)

	# refine df to make model predictions 
	bio <- sites[,c("bio2", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15")]
	# in absence records sometimes climateWNA spits out negative precip for bio14 in super arid regions e.g. ~ -0.012 mm
	# need to fix here
	bio$bio14[bio$bio14 < 0] <- 0	

	# log transform predictors 
	bio$bio3 <- log(bio$bio3+0.5); bio$bio10 <- log(bio$bio10+0.5); bio$bio12 <- log(bio$bio12+0.5); bio$bio14 <- log(bio$bio14+0.5)
	pred.dom <- bio # will needed
	head(sites, 2)

	
### 
# let's use occ_site_preds_sept2014 as source of 1981-2010 climate variables for transplant sites. any differences between model predictions in this file and those we will output now should be due to changes in model objects during am nat revisions	
	
sites <- read.csv("/Users/amyangert/Documents/GitClones/transplant_analysis/Data/raw_data/occ_site_preds_sept2014.csv")
	
sites$bio3 <- log(sites$bio3+0.5)
sites$bio10 <- log(sites$bio10+0.5)
sites$bio12 <- log(sites$bio12+0.5)
sites$bio14 <- log(sites$bio14+0.5)

pred.dom = sites[,c("bio2", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15")]
#####################
# run model predictions 

## LR prediction to sites 
results = sites[,2:3]

for (i in 1:10) {
  setwd(path.obj); 
  mod = get(paste("LR.mod2.",i, sep="")); 
  modprob = predict(mod, pred.dom, type="response", fun=predict); 
  modprob <- data.frame(modprob); 
  results <- cbind(results, modprob)
  }
	# LR ave
	LRavg <- rowMeans(results[,3:13]); 
	final <- cbind(results, LRavg)
		
## GAM prediction to sites 
	library(gam); 
	library(dismo); 
	setwd(path.obj)
for (i in 1:10) {
  mod = get(paste("GAM.mod4.",i, sep="")); 
  gamprob = predict(mod, pred.dom, type="response", fun=predict); 
  gammy <- data.frame(gamprob); 
  results <- cbind(results, gammy)
  }  		
	GAMavg <- rowMeans(results[,14:23]); 
	final <- cbind(final, GAMavg)

## RF prediction and classification
	setwd(path.obj); 
	library(randomForest); 
	library(raster)
	for (i in 1:10) {
	  setwd(path.obj); 
	  mod = get(paste("RF.mod1.",i, sep="")); 
	  rfprob = predict(mod, pred.dom, type = "prob", fun=predict, index=2); 
	  rf_preddy <- data.frame(rfprob); 
	  results <- cbind(results, rf_preddy)
	  }
		RFavg <- rowMeans(results[,c(25,27,29,31,33,35,37,39,41,43)]); 
		final <- cbind(final, RFavg)
		
## BRT prediction and classification
	library(dismo); 
	library(gbm); 
	library(raster); 
	setwd(path.obj)
	for (i in 1:10) {
	  mod = get(paste("BRT.mod4.",i, sep="")); 
	  brprob = predict(mod, pred.dom, n.trees=mod$gbm.call$best.trees, type="response"); 
	  br_preddy <- data.frame(brprob); 
	  results <- cbind(results, br_preddy)
	  }
	BRTavg <- rowMeans(results[,44:53]); 
	final <- cbind(final, BRTavg)

## MAX prediction and classification
	library(dismo); 
	library(raster); 
	setwd(path.obj)
	for (i in 1:10) {
	  mod = get(paste("MAX.mod1.",i, sep="")); 
	  maxprob = predict(mod, pred.dom); 
	  max_preddy <- data.frame(maxprob); 
	  results <- cbind(results, max_preddy)
	  }
	  # MAX ave
		MAXavg <- rowMeans(results[,54:63]); 
		final <- cbind(final, MAXavg)

write.csv(final, file = "/Users/amyangert/Documents/GitClones/transplant_analysis/Data/site_preds_average.csv", row.names = FALSE)


#
##
###
####
#####
######
#######   -  examine & plot data
######
#####
####
###
##
#
