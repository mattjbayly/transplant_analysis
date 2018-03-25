
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



	# Make dataframes for each variable (sequence along, to make partial plots below)
			daty <- DataPseudo3
			myVars <- colnames(daty)[2:ncol(daty)]
			VarSummary <- summary(daty)

			for(i in 1:length(myVars)){
				this <- mean(daty[,c(myVars[i])])
				assign(paste("MeaN",myVars[i], sep=""), this)  
			}
			
			MeaNFrame <- data.frame(bio15=MeaNbio15,
				SLOPE=MeaNSLOPE, terrough20C=MeaNterrough20C,
				logbio12=MeaNlogbio12, logDrainAre=MeaNlogDrainAre)
			MeaNFrame <- MeaNFrame[rep(seq_len(nrow(MeaNFrame)), each=150),]
				
			for(i in 1:length(myVars)){
				VaryName <- myVars[i]
				VaryVal <- daty[,c(myVars[i])]	
				NewVals <- seq(min(VaryVal*0.8), max(VaryVal*1.2), length.out=150)		
				tempframe <- MeaNFrame
				tempframe[,c(VaryName)] <- NewVals
				assign(paste("Frame",myVars[i], sep=""), tempframe) 
				toSav <- get(paste0("Frame",myVars[i]))
				assign(paste0("Frame",myVars[i], "_", 3), toSav) 
			}
		myVars2 <- paste0("Frame", colnames(daty)[2:ncol(daty)])
		rm(daty)
		u=3
		
		
		VarImport <- matrix(NA, nrow=length(colnames(DataPseudo3)[2:6]), ncol=5)
		rownames(VarImport) <- colnames(DataPseudo3)[2:6]
		colnames(VarImport) <- c("GLM", "GAM", "RF", "BRT", "MAX")

#######################################
# Predict models to data sets & to partial plot matrices of vairables  			
	######################################
	# GLM
		dat = DataPseudo3
		mod.lr = get(paste0("LR.mod2.",3, sep=""))
		lr.pred=predict(mod.lr, type="response") # model prediction
				for(i in 1:length(colnames(dat2)[2:6])){
				dat2 <- dat
				dat2[,i+1] <- dat2[sample(nrow(dat2)),i+1]
				lr.pred2=predict(mod.lr, newdata=dat2, type="response") # model prediction
				thisCor <- cor(lr.pred2, lr.pred)
				thisCor <- 1-thisCor
				VarImport[i,1] <- thisCor}
		assign(paste("LR.",3, sep=""), lr.pred)
		dat <- cbind(dat, get(paste0("LR.", 3)))
		colnames(dat)[7] <- "LRpred"
		assign(paste("DataPseudo",3, sep=""), dat)
		
		for(i in 1:length(myVars2)){
			dat = get(paste0(myVars2[i], "_", 3))
			mod.lr = get(paste0("LR.mod2.",3, sep=""))
			lr.pred=predict(mod.lr, newdata=dat, type="response") # model prediction
			assign(paste("LR.",i, sep=""), lr.pred)
			dat <- cbind(dat, get(paste0("LR.", i)))
			colnames(dat)[6] <- "LRpred"
			assign(paste0(myVars2[i], "_",3), dat) 
		}	
		
	######################################
	# GAM 
		dat = DataPseudo3
		mod = get(paste0("GAM.mod4.",3, sep="")) 
		pred=predict(mod, dat, type="response") # model prediction
				for(i in 1:length(colnames(dat2)[2:6])){
				dat2 <- dat
				dat2[,i+1] <- dat2[sample(nrow(dat2)),i+1]
				pred2=predict(mod, newdata=dat2, type="response") # model prediction
				thisCor <- cor(pred2, pred)
				thisCor <- 1-thisCor
				VarImport[i,2] <- thisCor}

		assign(paste("GAM.",3, sep=""), pred)
		dat <- cbind(dat, get(paste0("GAM.", 3)))
		colnames(dat)[8] <- "GAMpred"
		assign(paste("DataPseudo",3, sep=""), dat)
		
		for(i in 1:length(myVars2)){
			dat = get(paste0(myVars2[i], "_", 3))
			mod = get(paste0("GAM.mod4.",3, sep=""))
			pred=predict(mod, newdata=dat, type="response") # model prediction
			dat <- cbind(dat, pred)
			names(dat)[names(dat)=="pred"] <- "GAMpred"
			assign(paste0(myVars2[i], "_",3), dat) 
	}
		
		
	######################################	
	# RANDOM FOREST
		dat =  DataPseudo3
		mod = get(paste0("RF.mod1.",3, sep="")) 
		pred=predict(mod, type="prob")[,2] # predict from model
				for(i in 1:length(colnames(dat2)[2:6])){
				dat2 <- dat
				dat2[,i+1] <- dat2[sample(nrow(dat2)),i+1]
				pred2=predict(mod, newdata=dat2, type="prob")[,2] # model prediction
				thisCor <- cor(pred2, pred)
				thisCor <- 1-thisCor
				VarImport[i,3] <- thisCor}		
	
		assign(paste("RF.",3, sep=""), pred)
		dat <- cbind(dat, get(paste0("RF.", 3)))
		colnames(dat)[9] <- "RFpred"
		assign(paste("DataPseudo",3, sep=""), dat)
	
	
		for(i in 1:length(myVars2)){
			dat = get(paste0(myVars2[i], "_", 3))	
			mod = get(paste0("RF.mod1.",3, sep=""))
			pred=predict(mod, newdata=dat, type="prob")[,2] # predict from model
			dat <- cbind(dat, pred)
			names(dat)[names(dat)=="pred"] <- "RFpred"
			assign(paste0(myVars2[i], "_",3), dat) 
	}
	
	######################################
	# BRT / GBM
		dat = DataPseudo3
		mod = get(paste0("BRT.mod4.",3, sep="")) 
		pred = predict(mod, newdata=dat, type="response", n.trees=50)	
				for(i in 1:length(colnames(dat2)[2:6])){
				dat2 <- dat
				dat2[,i+1] <- dat2[sample(nrow(dat2)),i+1]
				pred2=predict(mod, newdata=dat2, type="response", n.trees=50) # model prediction
				thisCor <- cor(pred2, pred)
				thisCor <- 1-thisCor
				VarImport[i,4] <- thisCor}		
		
		assign(paste("BRT.",3, sep=""), pred)
		dat <- cbind(dat, get(paste0("BRT.", 3)))
		colnames(dat)[10] <- "BRTpred"
		assign(paste("DataPseudo",3, sep=""), dat)
		
		
		for(i in 1:length(myVars2)){
			dat = get(paste0(myVars2[i], "_", 3))
			mod = get(paste0("BRT.mod4.",3, sep=""))
			pred = predict(mod, newdata=dat, type="response", n.trees=50)	
			dat <- cbind(dat, pred)
			names(dat)[names(dat)=="pred"] <- "BRTpred"
			assign(paste0(myVars2[i], "_",3), dat) 
	}

	######################################
	# MAXENT
		dat = DataPseudo3
		mod = get(paste0("MAX.mod1.",3, sep="")) 
		pred = predict(mod, dat)
				for(i in 1:length(colnames(dat2)[2:6])){
				dat2 <- dat
				dat2[,i+1] <- dat2[sample(nrow(dat2)),i+1]
				pred2=predict(mod, dat2) # model prediction
				thisCor <- cor(pred2, pred)
				thisCor <- 1-thisCor
				VarImport[i,5] <- thisCor}		

		assign(paste("MAX.",3, sep=""), pred)
		dat <- cbind(dat, get(paste0("MAX.", 3)))
		colnames(dat)[11] <- "MAXpred"
		assign(paste("DataPseudo",3, sep=""), dat)
		
		
		
		
		for(i in 1:length(myVars2)){
			dat = get(paste0(myVars2[i], "_", 3))
			mod = get(paste0("MAX.mod1.",3, sep=""))
			pred = predict(mod, dat)	
			dat <- cbind(dat, pred)
			names(dat)[names(dat)=="pred"] <- "MAXpred"
			assign(paste0(myVars2[i], "_",3), dat) 
	}	
	
	
#########################################################################	
#########################################################################
#########################################################################
# EXPLORE PARTIAL PLOTS AND PREDICTED RESPONSE IDEAS
	
	
# TAKE A LOOK AT CORRELATIONS BETWEEN MODEL TYPES 
	setwd(path.fig)	
	pdf(file = "Variable Response Plots2.pdf", width=(8.5 - (1.25 + 0.75)), height=(6.5), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(2,2))


	#			par(mfrow=c(2,2))
	#			plot(DataPseudo3$GAMpred, DataPseudo3$LRpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$GAMpred, DataPseudo3$RFpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$GAMpred, DataPseudo3$BRTpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$GAMpred, DataPseudo3$MAXpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			par(mfrow=c(2,2))
	#			plot(DataPseudo3$LRpred, DataPseudo3$GAMpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$LRpred, DataPseudo3$RFpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$LRpred, DataPseudo3$BRTpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$LRpred, DataPseudo3$MAXpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			par(mfrow=c(2,2))
	#			plot(DataPseudo3$RFpred, DataPseudo3$LRpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$RFpred, DataPseudo3$GAMpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$RFpred, DataPseudo3$BRTpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$RFpred, DataPseudo3$MAXpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			par(mfrow=c(2,2))
	#			plot(DataPseudo3$BRTpred, DataPseudo3$LRpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$BRTpred, DataPseudo3$GAMpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$BRTpred, DataPseudo3$RFpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$BRTpred, DataPseudo3$MAXpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			par(mfrow=c(2,2))
	#			plot(DataPseudo3$MAXpred, DataPseudo3$LRpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$MAXpred, DataPseudo3$GAMpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$MAXpred, DataPseudo3$BRTpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		
	#			plot(DataPseudo3$MAXpred, DataPseudo3$RFpred, pch=21, bg=as.factor(DataPseudo3$PRESABS))		

############################################################################################
############################################################################################
############################################################################################
############################################################################################

############################################################################################
############################################################################################
############################################################################################
############################################################################################

# ESTIMATES OF VARIABLE IMPORTANCE 
	# Idea taken from the Biomod2 package variables_importance  function {biomod2}
	# credit to Damien Georges Documentation reproduced from package biomod2, version 3.1-25. License: GPL-2
	#http://www.inside-r.org/packages/cran/biomod2/docs/variables_importance

	# It's more or less base on the same principle than randomForest variables
	# importance algorithm. The principle is to shuffle a single variable of the 
	# given data. Make model prediction with this 'shuffled' data.set. Then we 
	# compute a simple correlation (Pearson's by default) between references 
	# predictions and the 'shuffled' one. The return score is
	# 1-cor(pred_ref,pred_shuffled).
	# The highest the value, the more influence the variable has on the model.
	# A value of this 0 assumes no influence of that variable on the model.
	# Note that this technique does not account for interactions between the variables.
#####################################################################################
############################################################################################
############################################################################################
	VarImport
	
	
	
	
				
# TEMPLATE PLOT FOR A SINGLE VARIABLE 		
		plot(Framelogbio12_3$logbio12, Framelogbio12_3$MAXpred, col="black",
		type="n", main="Maxent Predictions",
		xlab=expression(paste("Bio12: " , italic("Total annual discharge"), " mm", sep="")),
		ylab=expression(paste("Predicted Suitability (", italic(" MAX"), ")", sep="")),
		ylim=c(0,1))
		palette(makeTransparent("firebrick3","black", alpha=0.5))
		points(DataPseudo3$logbio12, DataPseudo3$MAXpred, pch=c(19, 4)[as.numeric(as.factor(DataPseudo3$PRESABS))], col=as.factor(DataPseudo3$PRESABS), cex=0.8)		
		lines(Framelogbio12_3$logbio12, Framelogbio12_3$MAXpred, type="l", lwd=2.5, col="red")

# Loop through all variables and make a separate pdf for each model 		
	# Make a list of all the vairables to use for the x-axis labes 
		
	VarNames <- c(
		expression(paste("Bio15: " , italic("Discharge Seasonality"), " (coefficent of variation)", sep="")),
		expression(paste("Slope " ,  ""^o*" ", sep="")),
		expression(paste("TerRough: " , italic("Topographical Roughness"), sep="")),
		expression(paste("Bio12: log(" , italic("Total annual discharge"), ")", sep="")),
		expression(paste("DrainAre: log(" , italic("Upstream drainage area"), ")", sep=""))
	)
	
	# Start loop with model types (GLM,GAM, RF ect.)
	preddy <- c("LRpred",   "GAMpred",    "RFpred",   "BRTpred",   "MAXpred")
	modName <- c("GLM", "GAM", "RF", "BRT", "MAXENT")
	Varys <- c("bio15",  "SLOPE",  "terrough20C",  "logbio12", "logDrainAre")


	
	# MAKE RESPONSE PLOTS FOR ALL VARIABLES & MODELS 

	for(j in 1:length(preddy)){	
	    Thispreddy <- preddy[j]
			par(mfrow=c(3,2)) # 8 variables plus one variable importance plot 
			for(i in 1:length(Varys)){
				toPlot <- get(paste0(myVars2[i], "_", 3))	
				
				plot(toPlot[,c(Varys[i])], toPlot[,c(preddy[j])], col="black",
					type="n", main=paste0(modName[j], " Predictions"),
					xlab=VarNames[i],
					ylab="Predicted Suitability (0-1)",
					ylim=c(0,max(toPlot[,c(Thispreddy)])*1.05))
					
				palette(makeTransparent("black","firebrick3", alpha=0.4))
				
			#	points(DataPseudo3[,c(Varys[i])], DataPseudo3[,c(preddy[j])],
			#		pch=c(4, 19)[as.numeric(as.factor(DataPseudo3$PRESABS))], col=as.factor(DataPseudo3$PRESABS), cex=0.8)		
				toPlotB <- get(paste0(myVars2[i], "_", 3))	
				lines(toPlotB[,c(Varys[i])], toPlotB[,c(preddy[j])], type="l", lwd=2.5, col="red")
				
			}
			# need blank plot for now until we can figure out variable importance 
			 #plot.new() # remvoe this later
			
			barplot(VarImport[,j], horiz=TRUE,
				xlab="Variable Importance", col="darkred", border="black", lwd=1.5, xlim=c(0, 1),
				las=2)
		
			}
		dev.off()

		


	
	