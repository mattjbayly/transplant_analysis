##########################
# MAXENT
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
	


library(SDMTools)
accuracy.max = function(dat, mod, modl) {
	## get model predictions
	mod.pred = predict(mod,dat)
	
	## make dataframe with 
		## col1 = model label (required by optimal.thresholds)
		## col2 = presence/absence
		## col3 = model prediction	
	datx=cbind(modl,dat[1],mod.pred) # build dataframe w/model predictions

	## evaluate model
	mod.val=evaluate(dat[dat$PRESABS==1,2:ncol(dat)], dat[dat$PRESABS==0,2:ncol(dat)], mod.MAX) 
library(dismo)
	mod.val=evaluate(p=datx[datx$PRESABS==1,3],a=datx[datx$PRESABS==0,3]) 

	## determine thresholds
	mod.cut=threshold(mod.val) # view maxent thresholds
	mod.cut = as.data.frame(c(mod.cut$spec_sens, mod.cut$kappa))
	mod.cut$method = c("sensspec", "kappa")
	names(mod.cut)[1] = "pred"

	## generate confusion matrices
	mod.cfmat.specsens=table(datx[[2]], factor(as.numeric(datx$mod.pred>=mod.cut$pred[1])))
	mod.cfmat.kappa=table(datx[[2]], factor(as.numeric(datx$mod.pred>=mod.cut$pred[2])))

	## calculate model accuracies with standard deviation=F
	mod.acc=presence.absence.accuracy(datx, threshold=mod.cut$pred, st.dev=FALSE) 
	tss=mod.acc$sensitivity + mod.acc$specificity - 1 # code TSS metric
	mod.acc=cbind(mod.acc[1:7],tss) # bind all metrics	
}

cv.accuracy.max = function(dat, mod, modl, x.fold, n.col) {
	#assign to cross folds and predict
	dat.xf = sample(rep(c(1:5), length=nrow(dat))) # vector of random xfolds
	mod.predXF = rep(0, length=nrow(dat)) # empty vector of 0
	for (j in 1:x.fold) {
		tr = dat[dat.xf!=j, ] # training not eq. i
		te = dat[dat.xf==j, ] # test eq. i
		mx = maxent(tr[2:n.col], tr[1]) # maxent model on training
		mod.predXF[dat.xf==j] = predict(mx, te) # predict to test
		}
	dat <- cbind(dat, mod.predXF)

	## evaluate model, determine optimal threshold
	mod.val=evaluate(p=dat[dat$PRESABS==1,2:ncol(dat)], a=dat[dat$PRESABS==0,2:ncol(dat)], mod) 
	mod.cutXF = threshold(mod.val)
	mod.cutXF = as.data.frame(c(mod.cutXF$spec_sens, mod.cutXF$kappa))
	mod.cutXF$method = c("sensspec", "kappa")
	names(mod.cutXF)[1] = "cvpred"

	## build testing dataframe using model predictions
	dat2XF=cbind(modl, dat[1], mod.predXF) 

	## calculate model accuracies with standard deviation=F
	mod.accXF = presence.absence.accuracy(dat2XF, threshold=mod.cutXF$cvpred, st.dev=FALSE)
	tss = mod.accXF$sensitivity + mod.accXF$specificity - 1 # code TSS metric
	mod.accXF = cbind(mod.accXF[1:7], tss) # bind all metrics
}

ext.accuracy.max = function(mod, ext, modl) {
	## predict to new data
	mod.epred = predict(mod, ext) 

	## build testing dataframes using model predictions
	datx = cbind(modl, ext[1], mod.epred) # bind obs and predictions

	## evaluate model and find thresholds
	mod.eval = evaluate(p=ext[ext$PRESABS==1,c(2:ncol(ext))], a=ext[ext$PRESABS==0,c(2:ncol(ext))], mod) 
	mod.ecut = threshold(mod.eval)
	mod.ecut = as.data.frame(c(mod.ecut$spec_sens, mod.ecut$kappa))
	mod.ecut$method = c("sensspec", "kappa")
	names(mod.ecut)[1] = "epred"

	## generate confusion matrix
	mod.ecfmat.sensspec = table(datx[[2]], factor(as.numeric(datx$mod.epred >= mod.ecut$epred[1])))
	mod.ecfmat.kappa = table(datx[[2]], factor(as.numeric(datx$mod.epred >= mod.ecut$epred[2])))

	## calculate model accuracies with standard deviation=F
	mod.eacc = presence.absence.accuracy(datx, threshold=mod.ecut$epred, st.dev=FALSE) 
	tss = mod.eacc$sensitivity + mod.eacc$specificity - 1 # code TSS metric
	mod.eacc=cbind(mod.eacc[1:7], tss) # bind all metrics	
}	
		
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
# MAXENT
################################################################################
######## START MAXENT MODEL: TRUE PRESENCES & BUFFERED-POINTS PSEUDOABSENCES
library(dismo)                         	
setwd(path.obj)
for (i in 1:4) {
	dat = get(setNames[i])
	
	mod.MAX=maxent(dat[,c(2:6)],dat$PRESABS, args=c("betamultiplier=2.5", "threshold=FALSE")) 
	assign(paste("MAX.mod1.",i, sep=""), mod.MAX)
	save(mod.MAX, file=paste("MAX.mod1.",i,".pseudo11.Rda", sep=""))
	}
######## END MAXENT MODEL: TRUE PRESENCES & BUFFERED-POINTS PSEUDOABSENCES
################################################################################


################################################################################
######## START RESUBSITUTION ACCURACY, MODEL = MAXENT
		library(PresenceAbsence) 
		setwd(path.code)
		#source("accuracy.R")
		accs=c()
		for (i in 1:4) {
			dat = get(setNames[i])
			#dat <- dat[,c(2:ncol(dat))]
			mod = get(paste("MAX.mod1.",i, sep=""))
			modl = "mod.MAX" # add var to keep track of model
			temp = accuracy.max(dat, mod, modl)
			temp$rep = i
			temp$thresh = c("SensSpec", "Kappa")
			temp$model = "MAX.mod1"
			accs = rbind(accs, temp)	
		}	
		setwd(path.obj)
		save(accs, file="MAX.mod1.accs.pseudo11.Rda")	

		
		######## END INTERNAL ACCURACY & CLASSIFICATION 
################################################################################


################################################################################
######## START CROSS-VALIDATION METRICS FROM MAXENT
library(PresenceAbsence) 
setwd(path.code)
#source("accuracy.R")
cv.accs=c()
modl="mod.MAX" # var placeholder
x.fold=5
n.col=6
for (i in 1:4) {
	## call up replicate training data and predictions
	dat = get(setNames[i])
	mod = get(paste("MAX.mod1.",i, sep=""))
	temp = cv.accuracy.max(dat, mod, modl, x.fold, n.col)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "MAX.mod1"
	cv.accs=rbind(cv.accs, temp)
	}
setwd(path.obj)
save(cv.accs, file="MAX.mod1.cvaccs.pseudo11.Rda")
######## END CROSS-VALIDATION METRICS FROM MAXENT
################################################################################


################################################################################
######## START EXTERNAL ACCURACY METRICS FROM MAXENT
## predict to independent occupancy dataset
library(PresenceAbsence)   
setwd(path.code)
#source("accuracy.R")
ext.accs=c()
modl = "mod.MAX" # var placeholder
for (i in 1:4) {
	## pull in replicate model
	mod =  get(paste("MAX.mod1.",i, sep=""))
	temp = ext.accuracy.max(mod, ext, modl)
	temp$rep = i
	temp$thresh = c("SensSpec", "Kappa")
	temp$model = "MAX.mod1"
	ext.accs=rbind(ext.accs, temp)
	}
setwd(path.obj)
save(ext.accs, file="MAX.mod1.extaccs.pseudo11.Rda")
######## END EXTERNAL ACCURACY METRICS FROM MAXENT
################################################################################

