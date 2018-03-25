# GAM

library(PresenceAbsence)
library(DAAG)
library(gam)
library(car)
		# set directories to occupancy paper to extract data from (don't write here)
		path.root = "C:/Users/DW/Desktop/Main.Backup.15.03.25/PROJECT COMPONENTS/GIS PROJECTS/Hydrology" 
		path.ecocrop = paste(path.root, "/ecocrop/mergeOutput", sep="")
		path.KMLs <- paste0(path.root, "/Card Verified Recrods/separateKML")
		path.dat <- paste0(path.root, "/Model_Data")
		path.streamRast <- "E:/stream_gis_all"
		path.exDet <- paste0(path.root, "/ExDet")
		path.code <- paste0(path.root, "/R_code")
		path.fig <- paste0(path.root, "/Figures")
		path.obj <- paste0(path.root, "/objects")
		path.ecocrop = paste(path.root, "/ecocrop/mergeOutput", sep="")
		path.KMLs <- paste0(path.root, "/Card Verified Recrods/separateKML")
		path.dat <- paste0(path.root, "/Model_Data")
		path.objects <- paste0(path.root, "/objects")
		path.code <- paste0(path.root, "/R_code")

		# Projections: 
		prj.wgs = "+proj=longlat +ellps=WGS84"
		
		## load libraries 
		library(dismo)
		library(raster)
		library(rgdal)
		library(rgeos)

fileNames <- c("DataPseudo1.csv", "DataPseudo2.csv", "DataPseudo3.csv", "DataPseudo4.csv")
setNames <- c("DataPseudo1", "DataPseudo2", "DataPseudo3", "DataPseudo4")
		
		
######################################################
######################################################
######################################################
######################################################
######################################################
# Load in data
setwd(path.dat) 
	for (i in 1:4) {
		dat = read.csv(fileNames[i])
		dat <- dat[,c("PRESABS",
			"bio15", "SLOPE", "terrough20C", "bio12", "DrainAre")]
		dat <- dat[!is.na(dat$terrough20C),]
		dat$logbio12 <- log(dat$bio12 + 1)
		dat$logDrainAre <- log(dat$DrainAre + 1)
		drops <- c("bio12", "DrainAre")
		dat <- dat[,!(names(dat) %in% drops)]
		assign(setNames[i], dat)
		rm(dat)
		}
	head(DataPseudo1, 2)
	head(DataPseudo4, 2)
	dim(DataPseudo4)

	setwd(path.code)
	source("modforms.R")
	setwd(path.obj)
				
######################################################
######################################################
######################################################
######################################################
######################################################	
######## START COMPARATIVE GAM MODELLING
setwd(path.code)
source("modforms.R")
setwd(path.obj)
library(gam)

for (i in 1:4) {
	dat = get(setNames[i])
	mod1.GAM=gam(mod.form.4(dat,1,2), family=binomial, data=dat)
	assign(paste("GAM.mod1.",i, sep=""), mod1.GAM)
	mod2.GAM=gam(mod.form.3(dat,1,2), family=binomial, data=dat)
	assign(paste("GAM.mod2.",i, sep=""), mod2.GAM)
	mod3.GAM=gam(mod.form.2(dat,1,2), family=binomial, data=dat)
	assign(paste("GAM.mod3.",i, sep=""), mod3.GAM)
	mod4.GAM=step.gam(mod2.GAM,scope=list(
		"bio15"	=~1+ bio15 	+s(bio15,2) +s(bio15,3)	+s(bio15,4),
		"SLOPE"	=~1+ SLOPE 	+s(SLOPE,2) +s(SLOPE,3)	+s(SLOPE,4),
		"terrough20C"	=~1+ terrough20C 	+s(terrough20C,2) +s(terrough20C,3)	+s(terrough20C,4),
		"logbio12"	=~1+ logbio12 	+s(logbio12,2) +s(logbio12,3)	+s(logbio12,4),
		"logDrainAre"	=~1+ logDrainAre 	+s(logDrainAre,2) +s(logDrainAre,3)	+s(logDrainAre,4)),
		trace=F)
	assign(paste("GAM.mod4.",i, sep=""), mod4.GAM)
	setwd(path.obj)
	save(mod4.GAM, file=paste("GAM.mod4.",i,".pseudo11.Rda", sep=""))
	}
	
######## END MODEL COMPARATIVE GAM MODELLING
################################################################################

################################################################################
######## START RESUBSTITUTION ACCURACY COMPARISONS, MODEL = STEP GAM

library(PresenceAbsence)   
setwd(path.code)
source("accuracy.R")

accs = c()
for (i in 1:4) {
	dat = get(setNames[i])
	mod = get(paste("GAM.mod4.",i, sep=""))
	pred=predict(mod, type="response") # predict by model
	modl="mod4.GAM" # add var to keep track of model
	temp = accuracy(dat, pred, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "GAM.mod4"
	accs = rbind(accs, temp)
	}
setwd(path.obj)
save(accs, file="GAM.mod4.accs.pseudo11.Rda")

######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL= STEP GAM
################################################################################

################################################################################
######## START CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL=STEP GAM

library(PresenceAbsence)	
library(DAAG)          
setwd(path.code)
source("accuracy.R")

cv.accs = c()
for (i in 1:4) {
	## pull in replicate data and model
	mod = get(paste("GAM.mod4.",i, sep=""))
	dat = get(setNames[i])
	modl="mod4.GAM" # assign model to varname
	temp = cv.accuracy(mod, dat, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "GAM.mod4"
	cv.accs = rbind(cv.accs, temp)
	}
setwd(path.obj)
save(cv.accs, file="GAM.mod4.cvaccs.pseudo11.Rda")

######## END CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL= STEP GAM
################################################################################

################################################################################
######## START EXTERNAL VALIDATION CALCULATIONS, MODEL= STEP GAM
## predict to independent occupancy dataset
setwd(path.dat)
	ext = read.csv(file="OccSnapped.csv")
	ext <- ext[which(ext$lat > 35), ]
	ext$logbio12  <- log(ext$bio12 + 1)
	ext$logDrainAre  <- log(ext$DrainAre + 1)
	ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]

library(PresenceAbsence)
setwd(path.code)
source("accuracy.R")

ext.accs = c()
for (i in 1:4) {
	## pull in replicate model
	mod = get(paste("GAM.mod4.",i, sep=""))
	modl="mod4.GAM"  # assign model to varname
	temp = ext.accuracy(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "GAM.mod4"
	ext.accs = rbind(ext.accs, temp)
  }
setwd(path.obj)
save(ext.accs, file="GAM.mod4.extaccs.pseudo11.Rda")

##################################################
# set suitability to zero for all occupancy outside thermal envelope
		setwd(path.dat)
		ext = read.csv(file="OccSnapped.csv")
		ext$logbio12  <- log(ext$bio12 + 1)
		ext$logDrainAre  <- log(ext$DrainAre + 1)
		ext2 <- ext
		ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]
		setwd(path.ecocrop)
		library(raster)
		ecobin <- raster("EcoCropBin.tif")
		coordinates(ext2) = ~long + lat
		proj4string(ext2) = projection(ecobin)
		drop <- extract(ecobin, ext2)
		ext <- cbind(ext, drop)
		extPass <- ext[which(ext$drop=='1'),]
		extFail <- ext[which(ext$drop=='0'),]
		extFail$PRESABS <- 0
		ext <- rbind(extPass, extFail)
		ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]
		ext.accs2 = c()
		for (i in 1:4) {
			## pull in replicate model
			mod = get(paste("GAM.mod4.",i, sep=""))
			modl="mod4.GAM"  # assign model to varname
			temp = ext.accuracy(mod, ext, modl)
			temp$rep = i
			temp$thresh = c("SensSpec", "Kappa")
			temp$model = "GAM.mod4"
			ext.accs2 = rbind(ext.accs2, temp)
		  }
		setwd(path.obj)
		save(ext.accs2, file="GAM.mod4.extaccs.ecobin.Rda")


######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=STEP GAM
################################################################################


# Explore alternative thresholds 
library(SDMTools)

	
		for (i in 1:4) {
			## pull in replicate model
			mod = get(paste("GAM.mod4.",i, sep=""))
			preds <- predict(mod, newdata=ext, type="response") # model prediction
			threshsExt <- optim.thresh(ext$PRESABS, preds, threshold = 101)
			...
			....
			....
			......
			.......
			.........
			......
			
			
			
			modl="mod4.GAM"  # assign model to varname
			temp = ext.accuracy(mod, ext, modl)
			temp$rep = i
			temp$thresh = c("SensSpec", "Kappa")
			temp$model = "GAM.mod4"
			ext.accs2 = rbind(ext.accs2, temp)
		  }
