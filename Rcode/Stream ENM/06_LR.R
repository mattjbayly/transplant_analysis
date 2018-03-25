######################################################
# LOGISTIC REGRESSION 

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

for (i in 1:4) {
	## build initial model all variables: mod1.i.LR
	dat = get(setNames[i])
	mod1.LR=glm(mod.form.quad(dat,1,2),family=binomial,data=dat)
	assign(paste("LR.mod1.",i, sep=""), mod1.LR)
	#mod1.fit=100*(1-mod1.LR$deviance/mod1.LR$null.deviance)
	mod1.pred=predict(mod1.LR,type="response") # model prediction
	assign(paste("LR.mod1.",i,".pred", sep=""), mod1.pred)
	#summary(mod1.LR)        # full model summary stats
	save(mod1.LR, file=paste("LR.mod1.",setNames[i], ".rda", sep="")) # save model object

	## build parsimonious model w/ backwards variable reduction: mod2.i.LR
	mod2.LR=step(mod1.LR,trace=F) # backwards stepwise variable reduction
	assign(paste("LR.mod2.",i, sep=""), mod2.LR)
	#mod2.fit=100*(1-mod2.LR$deviance/mod2.LR$null.deviance) # model fit
	mod2.pred=predict(mod2.LR,type="response") # model prediction
	assign(paste("LR.mod2.",i,".pred", sep=""), mod2.pred)
	#summary(mod2.LR)  # reduced model summary	
	save(mod2.LR, file=paste("LR.mod2.",setNames[i], ".rda", sep="")) # save model object
	}

######## END BUILD FULL AND REDUCED MODELS
################################################################################




################################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
## requires pkg PresenceAbsence
## build testing dataframes using mod2.i predictions

	library(PresenceAbsence)   
	setwd(path.code)
	source("accuracy.R")

	accs = c()
	for (i in 1:4) {
		dat = get(setNames[i])
		pred = get(paste("LR.mod2.",i,".pred", sep=""))	
		modl="mod2.LR"
		temp = accuracy(dat, pred, modl)
		temp$rep = i
		temp$thresh = c("SensSpec", "Kappa")
		temp$model = "LR.mod2"
		accs = rbind(accs, temp)
		}
	setwd(path.obj)
	save(accs, file="LR.mod2.accs.pseudo11.Rda")

######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################


################################################################################
######## START CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
## perform  k-fold cross-validation; requires pkg DAAG

	library(PresenceAbsence)  
	library(DAAG)               
	setwd(path.code)
	source("accuracy.R")

	cv.accs = c()
	for (i in 1:4) {
		## pull in replicate data and model
		mod = get(paste("LR.mod2.",i, sep=""))
		dat = get(setNames[i])
		modl="mod2.LR" # assign model to varname
		temp = cv.accuracy(mod, dat, modl)
		temp$rep = i
		temp$thresh = c("SensSpec", "Kappa")
		temp$model = "LR.mod2"
		cv.accs = rbind(cv.accs, temp)
		}
	setwd(path.obj)
	save(cv.accs, file="LR.mod2.cvaccs.pseudo11.Rda")

######## END CROSS-VALIDATION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################


################################################################################
######## START EXTERNAL ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
## predict to independent occupancy dataset
setwd(path.dat)
	ext = read.csv(file="OccSnapped.csv")
	ext <- ext[which(ext$lat > 35), ]
	ext$logbio12  <- log(ext$bio12 + 1)
	ext$logDrainAre  <- log(ext$DrainAre + 1)
	ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]

	ext.accs = c()
	for (i in 1:4) {
		## pull in replicate model
		mod = get(paste("LR.mod2.",i, sep=""))
		modl="mod2.LR" # add var to keep track of model
		temp = ext.accuracy(mod, ext, modl)
		temp$rep = i
		temp$thresh = c("SensSpec", "Kappa")
		temp$model = "LR.mod2"
		ext.accs = rbind(ext.accs, temp)
		}
	setwd(path.obj)
	save(ext.accs, file="LR.mod2.extaccs.pseudo11.Rda")

######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################


	