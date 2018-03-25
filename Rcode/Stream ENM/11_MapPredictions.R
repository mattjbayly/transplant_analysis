##########################
# MAPS OF MODEL PREDICTIONS
##########################

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

setwd(path.dat)
	ext = read.csv(file="OccSnapped.csv")
	ext$logbio12  <- log(ext$bio12 + 1)
	ext$logDrainAre  <- log(ext$DrainAre + 1)
	ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]

# LOAD IN MODELS 
## load SDM models 
setwd(path.obj)
for (i in 1:4) {
	mod.lr = get(load(paste("LR.mod2.DataPseudo",i,".Rda", sep="")))
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

# LOAD IN RASTER LAYERS 
		path.streamRast <- "E:/stream_gis_all"
		setwd(path.streamRast)
		list.files(pattern = "\\.tif$")
		streamFiles <- c("ArbolateSu.tif", "bio12.tif", "bio13.tif", "bio14.tif", "bio15.tif", "bio16.tif", "bio17.tif",
			"DrainAre.tif", "SLOPE.tif", "stream.tif", "strem_el.tif", "terrough20C.tif", "TotDaSqKM.tif")  
		streamNames <- c("ArbolateSu", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
			"DrainAre", "SLOPE", "stream", "strem_el", "terrough20C", "TotDaSqKM")  
			
	# load in variables and make transformations
	setwd(path.streamRast)
	bio15 <- raster("bio15.tif")
	SLOPE <- raster("SLOPE.tif")
	terrough20C <- raster("terrough20C.tif")
	bio12 <- raster("bio12.tif")
	DrainAre <- raster("DrainAre.tif")
	logbio12 <- log(bio12+1)
	logDrainAre <- log(DrainAre+1)
	
	# create raster stack of bioclim predictors
	pred.dom=stack(bio15, SLOPE, terrough20C, logbio12, logDrainAre)            
	biovec = c("bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")
	names(pred.dom) = biovec

######################################################
######################################################
######################################################
######################################################
######################################################


## LR prediction and classification 
## LR prediction and classification 
setwd(path.obj)
for (i in 1:4) {
	mod = get(paste("LR.mod2.",i, sep=""))
	modprob = predict(pred.dom, mod, type="response", fun=predict)
	#assign(paste("LR.modprob.",i, sep=""), modprob)
	save(modprob, file=paste("LR.modprob.",i,".img", sep=""))
 	}		
	setwd(path.obj)
	writeRaster(LR.modprob.3, filename="LR.modprob.3.tif")
	writeRaster(LR.modprob.2, filename="LR.modprob.2.tif")
	writeRaster(LR.modprob.1, filename="LR.modprob.1.tif")
	writeRaster(LR.modprob.4, filename="LR.modprob.4.tif")
	
# FOR GAM 
	library(gam)
	library(dismo)
	setwd(path.obj)
	#for (i in 1:4) {
	i=3
	mod = get(paste("GAM.mod4.",i, sep=""))
	modprob = predict(pred.dom, mod, type="response", fun=predict)
	#assign(paste("GAM.modprob.",i, sep=""), modprob)
	writeRaster(modprob, filename=paste0("GAM.modprob.", i,".tif"))
	#} 	
	
# FOR RANDOM FOREST
	setwd(path.obj)
	library(randomForest) # load library
	library(raster)
for (i in 1:4) {
	mod = get(paste("RF.mod1.",i, sep=""))
	modprob = predict(pred.dom, mod, type = "prob", fun=predict, index=2)
	#assign(paste("RF.modprob.",i, sep=""), modprob)
	writeRaster(modprob, filename=paste0("RF.modprob.", i,".tif"))
	}
	

# BOOSTED REGRESSION TREES 
	library(dismo)
	library(gbm)
	library(raster)
	setwd(path.obj)
for (i in 1:4) {
	mod = get(paste("BRT.mod4.",i, sep=""))
	modprob = predict(pred.dom, mod, n.trees=mod$gbm.call$best.trees, type="response")
	#assign(paste("BRT.modprob.",i, sep=""), modprob)
	writeRaster(modprob, filename=paste0("BRT.modprob.", i,".tif"))
	}
	
## MAX prediction and classification
library(dismo) # load library
library(raster)
setwd(path.obj)
for (i in 1:10) {
	mod = get(paste("MAX.mod1.",i, sep=""))
	modprob = predict(mod, pred.dom)
	#assign(paste("MAX.modprob.",i, sep=""), modprob)
	writeRaster(modprob, filename=paste0("MAX.modprob.", i,".tif"))
	}	