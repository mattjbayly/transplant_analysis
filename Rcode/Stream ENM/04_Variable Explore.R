######################################################
# EXPLORE RELATIONSHIPS BETWEEN PRESENCE ABSENCE RECORDS & CORRELATION STRUCTURE OF VARIABLES,
# DO LOG TRANSFORMATION ECT. 
# 

## load libraries 
library(dismo)
library(raster)
library(rgdal)
library(rgeos)

## set pathnames - Matthew 
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


# Projections: 
prj.wgs = "+proj=longlat +ellps=WGS84"

###################
	# load files
	setwd(path.dat)
	dir()
	
	files <- c("DataPseudo1.csv", "DataPseudo2.csv", "DataPseudo3.csv", "DataPseudo4.csv")
	SetName <- c("DataPseudo1", "DataPseudo2", "DataPseudo3", "DataPseudo4")

	for(i in 1:length(files)){
		myf <- read.csv(file=files[i])
		#coordinates(myf) = ~long + lat
		#proj4string(myf) = CRS(projection(prj.wgs))
		myf$PRESABS <- as.factor(myf$PRESABS)
		myf$logArbolateSu <- log(myf$ArbolateSu+1)
		myf$logbio12 <- log(myf$bio12+1)
		myf$logbio13 <- log(myf$bio13+1)
		myf$logbio14 <- log(myf$bio14+1)
		#myf$logbio15 <- log(myf$bio15+1)
		myf$logbio16 <- log(myf$bio16+1)
		myf$logbio17 <- log(myf$bio17+1)
		myf$logDrainAre <- log(myf$DrainAre+1)
		#myf$logSLOPE <- log(myf$SLOPE+1)
		myf$logTotDaSqKM <- log(myf$TotDaSqKM+1)	
		assign(SetName[i], myf)	
	}
	head(DataPseudo2, 2)
 
	setwd(path.code)
	source('Xpairs.R')
	source('multiplot.R')

	streamNames <- c("ArbolateSu", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
		"DrainAre", "SLOPE", "stream", "strem_el", "terrough20C", "TotDaSqKM")  


########################################################################
# MAKE FIGURE SET FOR EACH OF THE PSEUDO DATASETS 


for(i in 1:length(SetName)){
		setwd(path.fig); pdf(file=paste0("04_Explor", SetName[i],".pdf"), width=11, height=8.5)
		CurrentSet <- get(SetName[i])
			###########################################################
			# MAKE HISTOGRAMS WHICH VARIABLES NEED LOG TRANSFORMATION
				#now make your lovely plot
				library(ggplot2)
				library(grid)

			g1 <- ggplot(CurrentSet, aes(log(ArbolateSu+1), fill = PRESABS)) + geom_density(alpha = 0.2) +
				ggtitle(SetName[i])
			g2 <- ggplot(CurrentSet, aes(log(bio12+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g3 <- ggplot(CurrentSet, aes(log(bio13+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g4 <- ggplot(CurrentSet, aes(log(bio14+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g5 <- ggplot(CurrentSet, aes(log(bio15+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g55 <- ggplot(CurrentSet, aes(bio15, fill = PRESABS)) + geom_density(alpha = 0.2)
			g6 <- ggplot(CurrentSet, aes(log(bio16+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g7 <- ggplot(CurrentSet, aes(log(bio17+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			#multiplot(g1, g2, g3, g4, g5, g55, g6, g7, cols=2)
			g1 <- ggplot(CurrentSet, aes(DrainAre, fill = PRESABS)) + geom_density(alpha = 0.2)
			g2 <- ggplot(CurrentSet, aes(SLOPE, fill = PRESABS)) + geom_density(alpha = 0.2)
			g3 <- ggplot(CurrentSet, aes(terrough20C, fill = PRESABS)) + geom_density(alpha = 0.2)
			g4 <- ggplot(CurrentSet, aes(TotDaSqKM, fill = PRESABS)) + geom_density(alpha = 0.2)
			g5 <- ggplot(CurrentSet, aes(log(DrainAre+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g6 <- ggplot(CurrentSet, aes(log(SLOPE+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g7 <- ggplot(CurrentSet, aes(log(terrough20C+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			g8 <- ggplot(CurrentSet, aes(log(TotDaSqKM+1), fill = PRESABS)) + geom_density(alpha = 0.2)
			#multiplot(g1, g2, g3, g4, g5, g6, g7, g8,cols=2)
			
		## Don't log transform these variables:
			# bio15
			# terrough20C
			# SLOPE	
			
		###########################################################
		# CORRELATION STRUCTURES WITH LATITUDE AND ELEVATION	
			setwd(path.code)
			source('Xpairs.R')
			CurrentSet <- CurrentSet[!is.na(CurrentSet$stream), ]; CurrentSet <- CurrentSet[!is.na(CurrentSet$bio12), ]; CurrentSet <- CurrentSet[!is.na(CurrentSet$strem_el), ]; CurrentSet <- CurrentSet[!is.na(CurrentSet$terrough20C), ]

			# Xpairs(CurrentSet[,c("logbio12", "logbio13", "logbio14", "bio15", "logbio16", "logbio17", "lat")])	
				moot <- as.matrix(CurrentSet[,c("logArbolateSu", "logDrainAre", "SLOPE", "terrough20C", "logbio12", "logbio13", "logbio14", "bio15", "logbio16", "logbio17")])
				M <- round(cor(moot, use = "pairwise.complete.obs", method = "spearman"), digits = 2)
				corrplot(M, method = "number")
				library(corrplot)
				dev.off()
				
				setwd(path.fig); pdf(file=paste0("04_ExplorB", SetName[i],".pdf"), width=11, height=8.5)
				set <- CurrentSet[,c("logArbolateSu", "logDrainAre", "SLOPE", "terrough20C", "logbio12", "logbio13", "logbio14", "bio15", "logbio16", "logbio17")]
				z <- princomp(na.omit(set), cor=TRUE)
				par(mfrow=c(1,1)); old.par <- par(mar = c(0, 0, 0, 0)); par(old.par)
				biplot(z)
				

			################################################################################
			######## START VARIABLE IMPORTANCE AMONG PREDICTORS
			## START function variable importance
			## Cheerfully borrowed and modified from:  Niklaus E. Zimmermann
			## This function evaluates the predictive power of individual
			##   predictor variables as available in sequence within a data frame

			#### START function variable importance
			varimp.glm=function(tr.spp,tr.var,pres,pf,pl) {
			  tmp.mat=matrix(ncol=2,nrow=(pl-pf+1))
			  for (i in pf:pl) {
				## option: linear+quadratic; linear only
				tmp=glm(tr.spp[,pres]~tr.var[,i]+I((tr.var[,i])^2),na.action=na.omit,family=binomial)
				#tmp=glm(tr.spp[,pres]~tr.var[,i],na.action=na.omit,family=binomial) # linear only glm
				tmp.mat[(i-pf+1),1]=tmp$aic
				tmp.mat[(i-pf+1),2]=(1-(tmp$deviance/tmp$null.deviance))
				}
				return(tmp.mat)
			  }
			#### END function variable importance

			## estimate VIP values => AIC & Adj deviance
			tr.vip=CurrentSet[,c("PRESABS", "bio15", "SLOPE", "terrough20C", "logbio12", "logbio13", "logbio14", "logbio16", "logbio17", "logDrainAre", "logTotDaSqKM")] # keep P/A & ALL predictors 
			pres=1                     # column for presence:absence
			v.start=2                  # column start predictor variables
			v.stop=ncol(tr.vip)        # last column predictor variables
			v.num=v.stop-1             # number predictor variables
			dev.fit=varimp.glm(tr.vip,tr.vip,pres,v.start,v.stop) # call VIP function
			x.labs=as.data.frame(names(tr.vip[2:v.stop]))
			dev.fit =cbind(as.data.frame(dev.fit), x.labs)                 # output matrix; col=1 AIC, col=2 Adj deviance
			names(dev.fit) = c("AIC", "AdjDev", "Var")
			dev.fit = dev.fit[order(dev.fit$AIC),]
			(dev.fit)
			setwd(path.obj)
			write.csv(dev.fit, file=paste0("Exp_Dev_", SetName[i], ".csv"))
			setwd(path.fig)
			
			
			## built basic barplot if desired
			d.max=ceiling(signif(max(dev.fit[,2]),2)*10)/10 # convoluted; max of y-axis
			ylim.r=range(0,d.max)                           # range y-axis
			x.labs=names(tr.vip[2:v.stop])                  # x-axis labels
			par(mfrow=c(1,1)); old.par <- par(mar = c(0, 0, 0, 0)); par(old.par)
			barplot(dev.fit[,2],col="darkgreen",ylim=ylim.r,
			  main=paste0("pred VIPs", SetName[i]),ylab="adj.D2",names=x.labs)
			abline(h=0); abline(mean(dev.fit[,2]),0,lt=3) # ref lines; dash=mean adj.dev

			######## END VARIABLE IMPORTANCE AMONG PREDICTORS	

		rm(dev.fit); rm(CurrentSet)
		dev.off()
}	# end loop 
		
	dev.off()
		
	
######################################################
######################################################
######################################################
######################################################
# GIVEN A CORRELATION THRESHOLD OF 0.80 (MAXIMUM CUTOFF) THE FOLLOWING 
# VARIABLE COMBINATIONS WORK 

		# pseduoset1 - bio15, SLOPE, terrough20C, logbio12, logDrainAre
		# pseudoset2 - bio15, SLOPE, terrough20C, logbio12, logDrainAre
		# pseudoset3 - bio15, SLOPE, terrough20C, logbio12, logDrainAre
		# pseudoset44 - bio15, SLOPE, terrough20C, logbio12, logDrainAre

		
# ~~~~~~~~ NUMBER OF PSEUDOS & DISTANCE DIDNT AFFECT VARIABLE SELECTION HERE 
		