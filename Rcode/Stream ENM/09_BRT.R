##########################
# BOOSTED REGRESSION TREES
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
	ext <- ext[which(ext$lat > 35), ]
	ext$logbio12  <- log(ext$bio12 + 1)
	ext$logDrainAre  <- log(ext$DrainAre + 1)
	ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]

######################################################
######################################################
######################################################
######################################################
######################################################	
# BOOSTED REGRESSION TREES 

################################################################################
######## START BRT MODELS
library(gbm)    # load gbm package for BRT
library(dismo)  # for BRT calls per Elith et al (2008) JAnimalEcol 77:802-813

setwd(path.obj)
for (i in 1:4) {
	## call up replicate training data
	dat = get(setNames[i])
	resp=paste("as.factor(",colnames(dat[1]),")",sep="")  # assign response to column number
	n.col=ncol(dat)                                       # number of columns
	pred=2:n.col     				                       # assign predictors to column numbers
	## basic BRT model - LR adjusted down
	mod3.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	    tree.complexity=2, learning.rate=0.1, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=15)
 	assign(paste("BRT.mod3.",i, sep=""), mod3.BRT)
 	save(mod3.BRT, file=paste("BRT.mod3.",i,".pseudo11.Rda", sep=""))
	## basic BRT model - now bump TC up
	mod4.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	    tree.complexity=3, learning.rate=0.1, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
 	assign(paste("BRT.mod4.",i, sep=""), mod4.BRT)
 	save(mod4.BRT, file=paste("BRT.mod4.",i,".pseudo11.Rda", sep=""))
	## basic BRT model - keep TC up and bump LR back up
	mod5.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli",
	    tree.complexity=3, learning.rate=0.01, bag.fraction=0.75, n.folds=5,
	    plot.main=TRUE, keep.fold.fit=TRUE, step.size=5)
 	assign(paste("BRT.mod5.",i, sep=""), mod5.BRT)
 	save(mod5.BRT, file=paste("BRT.mod5.",i,".pseudo11.Rda", sep=""))
	}
	
# Any failed models? - YES 
		j=4
			for (i in 1:4) {
		mod = get(paste("BRT.mod",j,".",i, sep=""))
		print(paste("BRT.mod",j,".",i, sep=""))
		print(mod$cv.statistics[[1]])
		}
		
# re-run failed models (BRT.mod4 reps - )
		# remove step size argument for failed models
		failed <- c(1,2,3,4)
		for(i in 2:length(failed)){
		dat = get(paste(setNames[i], sep=""))
		resp=paste("as.factor(",colnames(dat[1]),")",sep="")
		n.col=ncol(dat)
		pred=2:n.col
		mod4.BRT=gbm.step(data=dat, gbm.x=pred, gbm.y=1, family="bernoulli", tree.complexity=3, learning.rate=0.01, bag.fraction=0.75, n.folds=5, plot.main=TRUE, keep.fold.fit=TRUE)
		assign(paste("BRT.mod4.",failed[i], sep=""), mod4.BRT)
		save(mod4.BRT, file=paste("BRT.mod4.",failed[i],".pseudo11.Rda", sep=""))
		}
		
#--------------------------------------------------------------------------------------#			
#--------------------------------------------------------------------------------------#
		## examine BRT output
		dev = matrix(ncol=6, nrow=11)
		ntrees = matrix(ncol=6, nrow=11)
		for (j in 3:5) {
			for (i in 1:4) {
				mod = get(paste("BRT.mod",j,".",i, sep=""))
						if((is.numeric(mod$cv.statistics[[1]]))==T){
						dev[i,j] = mod$cv.statistics[[1]]
						ntrees[i,j] = mod$n.trees
						} else {
						dev[i,j] = "NA" # failed model 
						}
		}}
		names(dev) = c("mod1", "mod2","mod3", "mod4", "mod5", "mod6")
		for (j in 3:5) {
			dev[11,j] = mean(as.numeric(as.character(dev[1:10,j])), na.rm=TRUE)}
		dev = as.data.frame(dev)
			dev
		ntrees = as.data.frame(ntrees)
		names(ntrees) = c("mod1", "mod2","mod3", "mod4", "mod5", "mod6")
			for (j in 3:5) {ntrees[11,j] = mean(as.numeric(as.character(ntrees[1:10,j])), na.rm=TRUE)}
			ntrees
		# mod4 wins in terms of lowest deviance and ntrees
		
		
################################################################################
######## START INTERNAL ACCURACY CALCULATIONS, MODEL=BRT
		library(PresenceAbsence)   
		setwd(path.code)
		source("accuracy.R")
		## be sure to adjust model #
		accs = c()
		for (i in 1:4) {
			dat = get(setNames[i])
			mod = get(paste("BRT.mod4.",i, sep=""))
			#pred = mod$fold.fit
				pred = mod$fitted
				modl="BRT.mod4"              
				temp = accuracy.brt(dat, pred, modl)
				temp$rep = i
				temp$thresh = c("SensSpec", "Kappa")
				temp$model = "BRT.mod4"
				accs = rbind(accs, temp)	
		}
		setwd(path.obj)
		save(accs, file="BRT.mod4.resubaccs.pseudo11.Rda")
######## END INTERNAL ACCURACY CALCULATIONS, MODEL=RF
################################################################################

################################################################################
######## START EXTERNAL ACCURACY CALCULATIONS, MODEL=BRT
	## predict to independent occupancy dataset
	library(PresenceAbsence)   
	setwd(path.code)
	source("accuracy.R")
	## be sure to adjust model #
	extaccs = c()
	for (i in 1:4) {
		mod = get(paste("BRT.mod4.",i, sep=""))
		modl = "BRT.mod4"
		temp = ext.accuracy.brt(mod, ext, modl)
		temp$rep = i
		temp$thresh = c("SensSpec", "Kappa")
		temp$model = "BRT.mod4"
		extaccs = rbind(extaccs, temp)	
		}
	setwd(path.obj)
	save(extaccs, file="BRT.mod4.extaccs.pseudo11.Rda")
######## END EXTERNAL ACCURACY CALCULATIONS, MODEL=BRT
################################################################################

		
		
		
		