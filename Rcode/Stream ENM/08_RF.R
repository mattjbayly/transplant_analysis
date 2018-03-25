##################
# RANDOM FOREST 
##################	

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
# RANDOM FOREST

library(randomForest) 
setwd(path.code)
source("modforms.R")


setwd(path.obj)
for (i in 1:4) {
	dat = get(setNames[i])
	mod1.RF = randomForest(mod.form.lin(dat,1,2), importance=T, keep.forest=T, data=dat)           
	assign(paste("RF.mod1.",i, sep=""), mod1.RF)
	save(mod1.RF, file=paste("RF.mod1.",i,".pseudo11.Rda", sep=""))	
	mod1.pred = predict(mod1.RF, type="prob")[,2] # predict from model
	assign(paste("RF.mod1.",i,".pred", sep=""), mod1.pred)
	}

	
################################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF
library(PresenceAbsence)   # PresenceAbsence for accuracy metrics
setwd(path.code)
source("accuracy.R")
accs = c()
for (i in 1:4) {
	dat = get(setNames[i])
	pred = get(paste("RF.mod1.",i,".pred", sep=""))
	modl = "mod1.RF"                     # add var to keep track of model
	temp = accuracy(dat, pred, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "RF.mod1"
	accs = rbind(accs, temp)
	}
setwd(path.obj)
save(accs, file="RF.mod1.accs.pseudo11.Rda")
######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF
################################################################################

############################################################################
######## START EXTERNAL VALIDATION, MODEL=RF
## predict to independent occupancy dataset

library(PresenceAbsence)  
setwd(path.code)
source("accuracy.R")

setwd(path.dat)
	ext = read.csv(file="OccSnapped.csv")
	ext <- ext[which(ext$lat > 35), ]
	ext$logbio12  <- log(ext$bio12 + 1)
	ext$logDrainAre  <- log(ext$DrainAre + 1)
	ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]

ext.accs = c()
for (i in 1:4) {
	mod = get(paste("RF.mod1.",i, sep=""))
	modl="mod1.RF"                     # add var to keep track of model
	temp = ext.accuracy.rf(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "RF.mod1"
	ext.accs = rbind(ext.accs, temp)
	}
setwd(path.obj)
save(ext.accs, file="RF.mod1.extaccs.pseudo11.Rda")
######## END EXTERNAL VALIDATION CALCULATIONS, MODEL=RF
################################################################################

