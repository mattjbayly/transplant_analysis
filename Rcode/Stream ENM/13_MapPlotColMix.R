
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
	ext <- ext[which(ext$lat > 35), ]
	# drop out thermal envelop
	setwd(path.ecocrop)
	drop <- raster("EcoCropBin.tif")
	coordinates(ext)= ~long+lat
	projection(ext) = projection(drop)
	fins <- extract(drop, ext)
	
	ext <- data.frame(ext)
	ext$fins <- fins
	# only one case where a presence record is outside of the thermal envelope 
	ext <- ext[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]
	# drop out recrods under southern calibration limit 

	
	
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
##########################################################
#########################################################################
#########################################################################
#########################################################################
library(raster)
library(maptools)
library(gam)
library(randomForest)
library(dismo)
library(gbm)
library(biomod2) # general response plot methods for SDMs 

# Special fxn will use later. Obtained from http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
makeTransparent = function(..., alpha=0.5) {
					  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
					  alpha = floor(255*alpha)  
					  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
					  .makeTransparent = function(col, alpha) {
						rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
					  }
					  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
					  return(newColor)
					}

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

#############################################################
#########################################################################		
# GET MODEL PREDICTIONS FROM EXTERNAL TESTING DATASET & 
setwd(path.obj)
library(raster)
library(SDMTools)

				sets <- c("ext", "DataPseudo3")
	for(i in 1:length(sets)){
			# run model predictions 
				bio <- get(sets[i])
				pred.dom <- bio[,c("bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]
			## LR prediction to sites 
				setwd(path.obj)
				mod = LR.mod2.3
				modprob = predict(mod, pred.dom, type="response", fun=predict)
				modprob <- data.frame(modprob)
				bio <- cbind(bio, modprob)
			## GAM prediction to sites 
				library(gam); library(dismo); setwd(path.obj)
				mod = GAM.mod4.3
				gamprob = predict(mod, pred.dom, type="response", fun=predict)
				gammy <- data.frame(gamprob)
				bio <- cbind(bio, gammy)	
			## RF prediction and classification
				setwd(path.obj); library(randomForest); library(raster)
				mod = RF.mod1.3
				rfprob = predict(mod, pred.dom, type = "prob", fun=predict, index=2)
				rf_preddy <- data.frame(rfprob)
				bio <- cbind(bio, rf_preddy[1])		
			## BRT prediction and classification
				library(dismo); library(gbm); library(raster); setwd(path.obj)
				mod = BRT.mod4.3
				brprob = predict(mod, pred.dom, n.trees=mod$gbm.call$best.trees, type="response"); br_preddy <- data.frame(brprob)
				bio <- cbind(bio, br_preddy)
			## MAX prediction and classification
				library(dismo); library(raster); setwd(path.obj)
				mod = MAX.mod1.3
				maxprob = predict(mod, pred.dom)
				max_preddy <- data.frame(maxprob)
				bio <- cbind(bio, max_preddy)
				
			# NORMALIZE PREDICTED VALUES PRIOR TO ENSEMBLE
				colnames(bio) <- c("PRESABS","bio15","SLOPE","terrough20C","logbio12","logDrainAre","GLM","GAM","RF","BRT","MAX")				
				ENcolMeans <- colMeans(bio[,7:11])
				ENcolStDev <- apply(bio[,7:11], 2, sd)
				modNames <- c("GLM","GAM","RF","BRT","MAX")
				toEn <- bio[,7:11]
				for(j in 1:length(modNames)){
				toEn[,j] <- ((toEn[,j] - ENcolMeans[j])/ENcolStDev[j])
				}
				EN <- rowMeans(toEn)
				# set ensemble value
				bio$EN <- EN
				bio$EN <- rowMeans(bio[,7:11])

			# Thresholds	
				mm <- matrix(NA, nrow=3, ncol=6)
				threshy1 <- optim.thresh(bio$PRESABS, bio$GLM, threshold = 101)
				threshy2 <- optim.thresh(bio$PRESABS, bio$GAM, threshold = 101)
				threshy3 <- optim.thresh(bio$PRESABS, bio$RF, threshold = 101)
				threshy4 <- optim.thresh(bio$PRESABS, bio$BRT, threshold = 101)
				threshy5 <- optim.thresh(bio$PRESABS, bio$MAX, threshold = 101)
				threshy6 <- optim.thresh(bio$PRESABS, bio$EN, threshold = 101)
				mm[1,1] <- as.numeric(threshy1$maxKappa[1]); mm[2,1] <- as.numeric(max(unlist(threshy1[5])))
				mm[1,2] <- as.numeric(threshy2$maxKappa[1]); mm[2,2] <- as.numeric(max(unlist(threshy2[5])))
				mm[1,3] <- as.numeric(threshy3$maxKappa[1]); mm[2,3] <- as.numeric(max(unlist(threshy3[5])))
				mm[1,4] <- as.numeric(threshy4$maxKappa[1]); mm[2,4] <- as.numeric(max(unlist(threshy4[5])))
				mm[1,5] <- as.numeric(threshy5$maxKappa[1]); mm[2,5] <- as.numeric(max(unlist(threshy5[5])))
				mm[1,6] <- as.numeric(threshy6$maxKappa[1]); mm[2,6] <- as.numeric(max(unlist(threshy6[5])))
				mm[3,1] <- auc(bio$PRESABS, bio$GLM)
				mm[3,2] <- auc(bio$PRESABS, bio$GAM)
				mm[3,3] <- auc(bio$PRESABS, bio$RF)
				mm[3,4] <- auc(bio$PRESABS, bio$BRT)
				mm[3,5] <- auc(bio$PRESABS, bio$MAX)
				mm[3,6] <- auc(bio$PRESABS, bio$EN)
				mm <- data.frame(mm)
				mm <- cbind(data.frame(c("SensSpec", "MaxKappa", "rawAUC")),mm)
				colnames(mm) <- c("Thresh", "GLM","GAM","RF","BRT","MAX","EN")
				assign(paste0(sets[i], "_thresh"), mm)
				setwd(path.obj)
				write.csv(mm, file=paste0(sets[i], "_thresholds.csv"))
				assign(paste0(sets[i], "_accAuc"), bio)
				}
				
ext_thresh
DataPseudo3_thresh

accuracy(ext_accAuc$PRESABS, ext_accAuc$GLM, threshold = DataPseudo3_thresh[1,2])


# get AUC and accuracy values 
		setwd(path.obj)
		temp <- get(load("LR.mod2.accs.pseudo11.rda")); temp1 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp1
		temp <- get(load("LR.mod2.cvaccs.pseudo11.rda")); temp2 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp2
		temp <- get(load("LR.mod2.extaccs.pseudo11.rda")); temp3 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp3
		temp <- get(load("GAM.mod4.accs.pseudo11.rda")); temp4 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp4
		temp <- get(load("GAM.mod4.cvaccs.pseudo11.rda")); temp5 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp5
		temp <- get(load("GAM.mod4.extaccs.pseudo11.rda")); temp6 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp6
		temp <- get(load("RF.mod1.accs.pseudo11.rda")); temp7 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp7
		temp <- get(load("RF.mod1.extaccs.pseudo11.rda")); temp8 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp8
		temp <- get(load("BRT.mod4.resubaccs.pseudo11.rda")); temp9 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp9
		temp <- get(load("BRT.mod4.cvaccs.pseudo11.rda")); temp10 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp10
		temp <- get(load("BRT.mod4.extaccs.pseudo11.rda")); temp10 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp10
		temp <- get(load("MAX.mod1.accs.pseudo11.rda")); temp11 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp11
		temp <- get(load("MAX.mod1.cvaccs.pseudo11.rda")); temp12 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp12
		temp <- get(load("MAX.mod1.extaccs.pseudo11.rda")); temp13 <- temp[which(temp$rep==3 & temp$thresh=='SensSpec'),]; temp13


		setwd("C:/Users/DW/Desktop/temp.sept.30/R objects")
		temp <- get(load("LR.mod2.accs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("LR.mod2.cvaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("LR.mod2.extaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("GAM.mod4.accs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("GAM.mod4.cvaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("GAM.mod4.extaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("RF.mod1.accs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("RF.mod1.extaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("BRT.mod4.resubaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("BRT.mod4.extaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("MAX.mod1.accs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("MAX.mod1.cvaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		temp <- get(load("MAX.mod1.extaccs.pseudo11.rda")); temp <- temp[which(temp$thresh=='SensSpec'),]; summary(temp)
		setwd(path.obj)

		
#########################################################################
#########################################################################		
# Make Final plots
		
	setwd(path.obj)
	rglm <- raster("LR.modprob.3.tif")
		fun <- function(x) { x[x < 0.01] <- 0; x[x > 0.01] <- 1; return(x)} 
		rglm <- calc(rglm, fun)
		writeRaster(rglm, filename="binGLM.tif", overwrite=TRUE)
	rgam <- raster("GAM.modprob.3.tif")
		fun <- function(x) { x[x < 0.01] <- 0; x[x > 0.01] <- 1; return(x)} 
		rgam <- calc(rgam, fun)
		writeRaster(rgam, filename="binGAM.tif", overwrite=TRUE)
	rrf <- raster("RF.modprob.3.tif")
		fun <- function(x) { x[x < 0.01] <- 0; x[x > 0.01] <- 1; return(x)} 
		rrf <- calc(rrf, fun)
		writeRaster(rrf, filename="binRF.tif", overwrite=TRUE)
	rbrt <- raster("BRT.modprob.3.tif")
		fun <- function(x) { x[x < 0.12] <- 0; x[x > 0.12] <- 1; return(x)} 
		rbrt <- calc(rbrt, fun)
		writeRaster(rbrt, filename="binBRT.tif", overwrite=TRUE)		
	rmax <- raster("MAX.modprob.3.tif")
		fun <- function(x) { x[x < 0.35] <- 0; x[x > 0.35] <- 1; return(x)} 
		rmax <- calc(rmax, fun)
		writeRaster(rmax, filename="binMAX.tif", overwrite=TRUE)		
	# add together
	setwd(path.obj)
	rglm <- raster("binGLM.tif")
	rgam <- raster("binGAM.tif")
	rrf <- raster("binRF.tif")
	rbrt <- raster("binBRT.tif")
	rmax <- raster("binMAX.tif")
	alls <- stack(rglm, rgam, rrf, rbrt, rmax)
	comb <- calc(alls, fun=sum)
	
	comb10 <- aggregate(comb, fact=10, fun=max)
	writeRaster(comb10, filename="CombHyd10.tif", overwrite=TRUE)
	comb4 <- aggregate(comb, fact=4, fun=max)
	writeRaster(comb4, filename="CombHyd4.tif", overwrite=TRUE)
	# subtract out ecocrop envelope 
	setwd(path.ecocrop); dir()
	ecos <- raster("EcoCropBin.tif")
	ecosRE <- projectRaster(ecos, comb4, method="ngb")
	setwd(path.obj)
	writeRaster(ecosRE, filename="ecos4.tif", overwrite=TRUE)
	comb4 <- comb4*ecosRE
	writeRaster(comb4, filename="CombHyd4.tif", overwrite=TRUE)

	ecosRE <- projectRaster(ecos, comb10, method="ngb")
	setwd(path.obj)
	writeRaster(ecosRE, filename="ecos10.tif", overwrite=TRUE)
	comb10 <- comb10*ecosRE
	writeRaster(comb10, filename="CombHyd10.tif", overwrite=TRUE)

	
	#comb <- raster("binClimCOUNT10.tif")
	#comb10 <- aggregate(comb, fact=10, fun=max)
	#writeRaster(comb10, filename="CombClim10.tif", overwrite=TRUE)
	#comb <- raster("binClimCOUNT10.tif")
	#comb4 <- aggregate(comb, fact=4, fun=max)
	#writeRaster(comb4, filename="CombClim4.tif", overwrite=TRUE)

	# combination of predictions
	combclimb <- raster("binClimCOUNT10.tif")
	combclimb <- aggregate(combclimb, fact=10, fun=max)
	combhyd <- raster("CombHyd10.tif")
	combclimb <- crop(combclimb, combhyd)
	combclimb <- projectRaster(combclimb, combhyd, method="ngb")
	comb <- combclimb + combhyd
	writeRaster(comb, filename="CountCombined10.tif", overwrite=TRUE)
	comb <- raster("CountCombined10.tif")
	#vals <- unique(comb)
		
	
################################################################
#########################################################################
# MAKE COLOR SPACE
	setwd(path.code)
	source("val2col.R")
	blues <- c(1:6)
	reds <- c(seq(41,91, by=10))
	par(bg = "white")

#assign colors to grd levels
		#pal1 <- colorRampPalette(c("white", rgb(1,0,0)), space = "rgb")
		pal1 <- colorRampPalette(c("lightgrey", "red"), space = "rgb")
		col1 <- val2col(reds, col=pal1(6))
		pal2 <- colorRampPalette(c("lightgrey", "blue"), space = "rgb")
		col2 <- val2col(blues, col=pal2(6))
		par(mfrow=c(1,2)); plot(blues, col=col2, pch=19, cex=5); plot(reds, col=col1, pch=19, cex=5)
		col1 # red
		col2 # blue

	breakpoints <- c(0,1,2,3,4,5,10,11,12,13,14,15,
					20,21,22,23,24,25,30,31,32,33,34,35,
					40,41,42,43,44,45,50,51,52,53,54,55)
	breakpoints <- matrix(breakpoints, nrow=6, ncol=6)
		
		col3 <- matrix(NA*c(1:36), nrow=6, ncol=6)
		par(mfrow=c(1,1)); plot(1:7, col="white", xlab="i's", ylab="j's")
		for(j in 1:6){
			for(i in 1:6){
			coltmp <- (col2rgb(col1[j])/2) + (col2rgb(col2[i])/2)
			col3[i,j] <- rgb(coltmp[1], coltmp[2], coltmp[3], maxColorValue = 255)
			points(i, j, pch=15, cex=6, col=col3[i,j])
			text(i, j, label=breakpoints[i,j])
			text(i+0.2, j+0.2, label=col3[i,j], cex=0.6)
		}}
		
		col3B <- as.vector(col3)
		#plot(1:36, col=col3B, pch=19)
		col3B <- data.frame(col3B)

############################################################
############################################################
############################################################

		breakpoints <- c(0,1,2,3,4,5,10,11,12,13,14,15,
					20,21,22,23,24,25,30,31,32,33,34,35,
					40,41,42,43,44,45,50,51,52,53,54,55)

	col3B$val <- breakpoints
	colors1 <- as.character(col3B$col3B)	

# LEGEND PLOT 
		col3 <- matrix(NA*c(1:36), nrow=6, ncol=6)
		par(mfrow=c(1,1))
		plot(1:6, col="white", axes=FALSE, xlab="", ylab="", frame.plot=TRUE)
		for(j in 1:6){
			for(i in 1:6){
			coltmp <- (col2rgb(col1[j])/2) + (col2rgb(col2[i])/2)
			col3[i,j] <- rgb(coltmp[1], coltmp[2], coltmp[3], maxColorValue = 255)
			points(i, 7-j, pch=21, cex=3, bg=col3[i,j], col="lightgrey")
		}}	
		Axis(side=2, labels=c(0:5), at=6:1)
		Axis(side=3, labels=c(0:5), at=1:6)
		mtext("Climatic ENMs", side=2, line=2)
		mtext("Stream Habitat ENMs", side=3, line=1.5)
#.......................................	
# new color grid
	setwd(path.obj)
	cols <- read.csv(file="NewColgrid.csv", header=FALSE)
	cols <- as.matrix(cols)
	cols <- t(cols) # oops made the matrix backwards
	col3 <- cols
# LEGEND PLOT 
		par(mfrow=c(1,1))
		plot(1:6, col="white", axes=FALSE, xlab="", ylab="", frame.plot=TRUE)
		for(j in 1:6){
			for(i in 1:6){
			points(i, 7-j, pch=21, cex=3, bg=col3[i,j], col="lightgrey")
		}}	
		Axis(side=2, labels=c(0:5), at=6:1)
		Axis(side=3, labels=c(0:5), at=1:6)
		mtext("Climatic ENMs", side=2, line=2)
		mtext("Stream Habitat ENMs", side=3, line=1.5)
		
		colors1 <- as.vector(col3)
		col1 <- col3
		
############################################################
############################################################
############################################################	
# MAKE FINAL PLOTS
		setwd(path.code)
		source("plot_obj.R")
		
setwd(path.obj)	
# fix rasters only once
		e <- extent(-124.5642, -115, 34.99917, 47)
		e = as(e, "SpatialPolygons")
		proj4string(e) = CRS(prj.wgs)
		frame <- crop(countriesHigh, e)
		frame.aea <- spTransform(frame, CRS=CRS(prj.aea))
	# pres points
	pres = all[all$PRESABS=="1",] 	
	pres <- pres[which(pres$Latitude > 35),]
	coordinates(pres) = ~Longitude + Latitude
	proj4string(pres) = CRS(prj.wgs)
	pres <- spTransform(pres, CRS=CRS(prj.aea))

		
		comb <- raster("CountCombined10.tif")
		comb <- crop(comb, frame) # trim down 
		comb <- aggregate(comb, fact=2, fun=max, na.rm=TRUE)		
		comb.aea <- projectRaster(comb, crs=CRS(prj.aea), method = "ngb")
		comb.aea <- clip(comb.aea, frame.aea) # trim down  
		setwd(path.obj)
		writeRaster(comb.aea, filename="CountCombined10.aea.tif", overwrite=TRUE)

		comb <- raster("CombHyd10.tif")
		comb <- crop(comb, frame) # trim down  
		comb <- aggregate(comb, fact=2, fun=max, na.rm=TRUE)		
		comb.aea <- projectRaster(comb, crs=CRS(prj.aea), method = "ngb")
		comb.aea <- clip(comb.aea, frame.aea) # trim down  
		comb.aea[comb.aea<0] = 0
		setwd(path.obj)
		writeRaster(comb.aea, filename="CombHyd10.aea.tif", overwrite=TRUE)
	
		comb <- raster("CombClim10.tif")
		comb <- crop(comb, frame) # trim down  
		comb <- aggregate(comb, fact=2, fun=max, na.rm=TRUE)		
		comb.aea <- projectRaster(comb, crs=CRS(prj.aea), method = "ngb")
		comb.aea <- clip(comb.aea, frame.aea) # trim down  
		setwd(path.obj)
		writeRaster(comb.aea, filename="CombClim10.aea.tif", overwrite=TRUE)
			
		
setwd(path.fig)
#dir()

pdf(file = "13.1.combinePlotofPreds.pdf", height=(8.5 - (1.25 + 0.75)), width=(10-(0.75+0.75)), family="Times") #  - (0.75 + 0.75)
		par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
		par(mfrow=c(1,3), mar=c(0,0,0,0))

# STREAM HABITAT		
		plot(frame.aea, col="lightgrey")
		title(main="Stream habitat ENM", line=-6, cex=2)
		setwd(path.obj)
		comb <- raster("CombHyd10.aea.tif")	
		comb[comb<0] = 0
		plot(comb, col=col3[,1], add=TRUE, legend=FALSE)
		sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
		grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
		#points(pres, pch=1, cex=2, col="black")
		# add on legend
 		plot(comb, legend.only=TRUE, legend.width=2, legend.shrink=0.75, col=col3[,1], 
			axis.args=list(at=seq(0, 5, by=1), labels=seq(0, 5, by=1), cex.axis=0.8),
			smallplot=c(0.70,0.73, .1,.4), legend.args=list(text='Model Agreement', side=4, font=2, line=1.5, cex=0.8))

# CLIMATIC				
		plot(frame.aea, col="lightgrey")
		title(main="Climatic ENM", line=-6, cex=2)
		comb <- raster("CombClim10.aea.tif")
		comb[comb<0] = 0
		plot(comb, col=col3[1,], add=TRUE, legend=FALSE)
		sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
		grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
		#points(pres, pch=19, cex=0.3)
	 		plot(comb, legend.only=TRUE, legend.width=2, legend.shrink=0.75, col=col3[1,], 
			axis.args=list(at=seq(0, 50, by=10), labels=seq(0, 5, by=1), cex.axis=0.8),
			smallplot=c(0.70,0.73, .1,.4), legend.args=list(text='Model Agreement', side=4, font=2, line=1.5, cex=0.8))
	
# COMBINED PANNEL		
		plot(frame.aea, col="lightgrey")
		title(main="Combined", line=-6, cex=2)
		setwd(path.obj)
		comb <- raster("CountCombined10.aea.tif")
		comb[comb<0] = 0
		plot(comb, col=colors1, add=TRUE, legend=FALSE)
		sta.aea2 <- crop(sta.aea,frame.aea); plot(sta.aea2, add=T, border="black", lwd=0.2)
		grd.aea2 <- crop(grd.aea,frame.aea); plot(grd.aea2, add=T, lty="dashed", col="darkgrey", lwd=1)
		#points(pres, pch=19, cex=0.3)
		# reporject raster

		
		
#par(mai=c(0.2,0.1,5.3,0.5))
		
		
		dev.off()
		
		
		