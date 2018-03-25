######################################################
# RUN ECOCROP ON EXISTING PRES/ABSE RECORDS FROM OCCUPANCY 
# PAPER & DETERMINE IF A SUITABLE THRESHOLD EXISTS. 


## load libraryies 
library(dismo)
library(raster)

## set pathnames - Matthew 
# set directories to occupancy paper to extract data from (don't write here)
path.root = "C:/Users/DW/Desktop/temp.sept.30" 
path.dat = paste(path.root, "/data files", sep="")
path.dat.fix = paste(path.root, "/data files", sep="") # older files relocated to another directory
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.obj, "/ecoregions.shp", sep="")
path.bio = paste(path.obj, "/wc0.5", sep="")
path.cod=paste(path.root, "/R code", sep="")
path.fig=paste(path.root, "/figures", sep="")
path.sta = paste(path.obj, "/gz_2010_us_040_00_500k", sep="")

# load all records
setwd(path.dat)
all = read.csv("all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
pres =  all[all$PRESABS=="1",] 
abs  = all[all$PRESABS=="0",] 

# Will have to define a new species value in ecocrop for cardinalis 
head(ECOcrops, 3) # existing species in package

# Values from: 
		# Angert, Amy L., Seema N. Sheth, and John R. Paul. 
		# "Incorporating population-level variation in thermal performance into predictions 
		#of geographic range shifts." Integrative and comparative biology (2011): icr048.

		crop <- getCrop('Coffea excelsa L.')
		str(crop)
		crop@name <- "Scarlet Monkeyflower"
		crop@famname <- "Phrymaceae"
		crop@scientname <- "Mimulus cardinalis"
		crop@code <- "NA"
		crop@GMIN = 65
		crop@GMAX = 300 
		crop@KTMP = -4
		crop@TMIN = 2.58
		crop@TOPMN = 15.64
		crop@TOPMX = 35.05
		crop@TMAX = 48
		crop@RMIN = 0 
		crop@ROPMN = 0
		crop@ROPMX = 9900
		crop@RMAX = 9900
		LIG = NA; LIGR = NA; PP = NA; PPMIN = NA; PPMAX = NA; TEXT = NA; TEXTR = NA; DEP = NA; DEPR = NA; 
		DRA = NA; DRAR = NA; PHMIN = NA; PHOPMN = NA; PHOPMX = NA; PHMAX = NA; SAL = NA; SALR = NA; FER = NA; FERR = NA; LIMITS = NA
	
		# ECO CROP VALUES FOR PRESENCE AND ABSENCE RECORDS 
		pres$ecocrop <- NA
		for(i in 1:dim(pres)[1]){
			moo <- mean(ecocrop(crop, as.numeric(pres[i,c(21:32)]), as.numeric(pres[i,c(33:44)]), as.numeric(pres[i,c(45:56)]), rainfed=TRUE)@suitability)
			pres[i, c("ecocrop")] <- moo
		}
		
		abs$ecocrop <- NA
		for(i in 1:dim(abs)[1]){
			moo <- mean(ecocrop(crop, as.numeric(abs[i,c(21:32)]), as.numeric(abs[i,c(33:44)]), as.numeric(abs[i,c(45:56)]), rainfed=TRUE)@suitability)
			abs[i, c("ecocrop")] <- moo
			(print(i))
		}
		par(mfrow=c(2,1))
		hist(pres$ecocrop)
		hist(abs$ecocrop)
			
# Where should the cutoff be drawn?
		par(mfrow=c(1,1))
		hist(pres$ecocrop,  breaks = seq(0:0.15, by=0.01), xlim=c(0, 0.15))
		inspect <- head(pres[with(pres, order(ecocrop)), ], 20)
		#write.csv(inspect, file="checkThese.csv")
			
# Make latitude/elevation/ecopcrop heatmaps
	setwd("C:/Users/DW/Desktop/Main.Backup.15.03.25/PROJECT COMPONENTS/GIS PROJECTS/Hydrology/ecocrop/explore")
	pdf(file="presValues.pdf", width=11, height=8.5)
	library(akima) 
	library(fields) 		
	s <- interp(pres$Latitude,pres$Elevation,pres$ecocrop)
	image.plot(s)
	points(pres$Latitude, pres$Elevation, pch=21, cex=0.5, bg="white")
	text(pres$Latitude, pres$Elevation, labels=round(pres$ecocrop, 2), cex=0.5)
	und <- pres[which(pres$ecocrop < 0.05), ]
	points(und$Latitude, und$Elevation, pch=19, col="yellow")
	hist(pres$ecocrop)
	dev.off()

	#setwd(choose.dir())
	pdf(file="absValues.pdf", width=11, height=8.5)		
	s <- interp(abs$Latitude,abs$Elevation,abs$ecocrop)
	image.plot(s)
	und <- abs[which(abs$ecocrop < 0.05), ]
	points(abs$Latitude, abs$Elevation, pch=21, cex=0.5, bg="white")
	points(und$Latitude, und$Elevation, pch=19, col="yellow")
	#text(abs$Latitude, abs$Elevation, labels=round(abs$ecocrop, 2), cex=0.5)
	hist(abs$ecocrop)
	dev.off()

# Make a map 
	library(dismo); library(sp) 
	prj.wgs = "+proj=longlat +ellps=WGS84"
	coordinates(pres) = ~Longitude + Latitude
	proj4string(pres) = CRS(prj.wgs)
	rbPal <- colorRampPalette(c('lightgrey','yellow','red','purple', 'blue')) #Create a function to generate a continuous color palette
	pres$Col <- rbPal(100)[as.numeric(cut(pres$ecocrop, breaks = c(seq(0, 1, by=0.01))))]
	par(mfrow=c(1,2))
	pdf(file="pointvalues.pdf", width=11, height=8.5)		
	plot(pres, pch=21, cex=0.5, bg=pres$Col)
	coordinates(abs) = ~Longitude + Latitude
	proj4string(abs) = CRS(prj.wgs)
	rbPal <- colorRampPalette(c('lightgrey','yellow','red','purple', 'blue')) #Create a function to generate a continuous color palette
	abs$Col <- rbPal(100)[as.numeric(cut(abs$ecocrop, breaks = c(seq(0, 1, by=0.01))))]
	#plot(abs, pch=21, cex=0.5, bg=abs$Col)
		plot(pres, pch=21, bg=pres$Col, cex=0.5, main="threshold = 0.07", col="lightgrey")
		und <- abs[which(abs$ecocrop < 0.07), ]
		points(und, col="black", pch=19, cex=0.5)
		und <- pres[which(pres$ecocrop < 0.07), ]
		points(und, bg="red", pch=24, cex=1.2)
	plot(pres, pch=21, bg=pres$Col, cex=0.5, main="threshold = 0.06", col="lightgrey")
		und <- abs[which(abs$ecocrop < 0.06), ]
		points(und, col="black", pch=19, cex=0.5)
		und <- pres[which(pres$ecocrop < 0.06), ]
		points(und, bg="red", pch=24, cex=1.2)
	plot(pres, pch=21, bg=pres$Col, cex=0.5, main="threshold = 0.05", col="lightgrey")
		und <- abs[which(abs$ecocrop < 0.05), ]
		points(und, col="black", pch=19, cex=0.5)
		und <- pres[which(pres$ecocrop < 0.05), ]
		points(und, bg="red", pch=24, cex=1.2)
	plot(pres, pch=21, bg=pres$Col, cex=0.5, main="threshold = 0.04", col="lightgrey")
	und <- abs[which(abs$ecocrop < 0.04), ]
	points(und, col="black", pch=19, cex=0.5)
	und <- pres[which(pres$ecocrop < 0.04), ]
	points(und, bg="red", pch=24, cex=1.2)
		plot(pres, pch=21, bg=pres$Col, cex=0.5, main="threshold = 0.045", col="lightgrey")
		und <- abs[which(abs$ecocrop < 0.045), ]
		points(und, col="black", pch=19, cex=0.5)
		und <- pres[which(pres$ecocrop < 0.045), ]
		points(und, bg="red", pch=24, cex=1.2)
	
	dev.off()

	

##############################################################
# EVALUATE THRESHOLD BASED ON RASTERS FROM PRISM CLIMATE GROUP
##############################################################

	library(raster)
	path.main <- "C:/Users/DW/Desktop/Main.Backup.15.03.25/PROJECT COMPONENTS/GIS PROJECTS/Hydrology/"
	path.prism <- paste0(path.main, "PRISMrasters")
	setwd(path.prism); dir()[1:5]
	test <- raster('PRISM_ppt_30yr_normal_800mM2_01_asc.asc')

	e <- extent(-125.0208, -115, 24.0625, 49.9375) # extent object
	test <- crop(test, e)
	plot(test)
	
	# now do for all files
	files <- list.files(pattern = "\\.asc$") # $ sign means at the end of the file and \\. means nothing before (eg. .xlmasc)
	for(i in 1:length(files)){
		temp <- files[i]
		test <- raster(temp)
		test <- crop(test, e)
		newName <- strsplit(temp, "_asc.asc")[[1]]
		writeRaster(test, file=paste0(newName, ".tif"), bylayer = TRUE, 
			datatype = 'FLT4S', bandorder = 'tif', overwrite = TRUE)
	}
	
	
	# define monthly values 	
		files <- list.files(pattern = "\\.tif$") 
		tmins <- (files[grep(pattern = "tmin",files)])[1:12]
		taves <- (files[grep(pattern = "tmean",files)])[1:12]
		ppts <- (files[grep(pattern = "ppt",files)])[1:12]
	
		set <- c("tmins", "taves", "ppts")
		for(i in 1:length(set)){
			temp <- set[i]
			mySet <- get(temp)
				for(j in 1:length(mySet)){
					this <- mySet[j]
					moot <- raster(paste(this))
					assign(paste(set[i], j, sep="_"), moot)
				}
			}
		paste0("tmins_", 1:12)
		tmins <- stack(tmins_1, tmins_2, tmins_3, tmins_4, tmins_5, tmins_6, tmins_7, tmins_8, tmins_9, tmins_10, tmins_11, tmins_12)
		paste0("taves_", 1:12)
		taves <- stack(taves_1, taves_2, taves_3, taves_4, taves_5, taves_6, taves_7, taves_8, taves_9, taves_10, taves_11, taves_12)
		paste0("ppts_", 1:12)
		ppts <- stack(ppts_1, ppts_2, ppts_3, ppts_4, ppts_5, ppts_6, ppts_7, ppts_8, ppts_9, ppts_10, ppts_11, ppts_12)

	# run ecocrop (mean annual), crop is cardinalis 
		# will take a mintue or two 
		#ecoMap <- mean(ecocrop(crop, tmins, taves, ppts, rainfed=TRUE)@suitability)
		e2 <- extent(-121.30, -121.10, 39.60, 39.79) # extent object

	# For the function to run efficiently need to split main stack into tiles & process separately
	dir.lay <- "C:/Users/DW/Desktop/Main.Backup.15.03.25/PROJECT COMPONENTS/GIS PROJECTS/Hydrology/ecocrop/tempLayers"
	setwd(dir.lay)

	clipr <- seq(-125.020, -114.9958, by=((-114.9958 - -125.020)/800))
		
	for(j in 1:length(clipr)){	
		# make skinny strip
		e2 <- extent(clipr[j], clipr[j+1], 24.0625, 49.9375)
		tmins2 <- crop(tmins, e2)
		taves2 <- crop(taves, e2)
		ppts2 <- crop(ppts, e2)
				
		moo <- data.frame(Eco=NA, x=coordinates(tmins2)[,1], y=coordinates(tmins2)[,2])
		moo <- as.matrix(moo)
			for(i in 1:length(getValues(tmins2)[,1])){
				moo[i,1] <- mean(ecocrop(crop, getValues(tmins2)[i,], getValues(taves2)[i,], getValues(ppts2)[i,], rainfed=FALSE)@suitability)
				}
		update <- subset(ppts2, 2)
		update <- setValues(update,  moo[,c("Eco")])		
		writeRaster(update, file=paste0("Tile_", j,".tif"), bylayer = TRUE, 
			datatype = 'FLT4S', bandorder = 'tif', overwrite = TRUE)
		rm(moo); rm(e2); rm(update)
		}

#
##
###
####
##### # Need to paste together ecocrop tiles in QGIS - too slow in R. 
####
###
##
#
		

		
##############################################################
# Get ecocrop values from final raster
##############################################################

dir.ecoRast <- "C:/Users/DW/Desktop/Main.Backup.15.03.25/PROJECT COMPONENTS/GIS PROJECTS/Hydrology/ecocrop/mergeOutput"
setwd(dir.ecoRast) 
EcoCrop <- raster("EcoCropOut.tif")
# par(mfrow=c(1,1)); plot(EcoCrop); plot(pres, add=T)
EcoCrop
PresVals <- extract(EcoCrop, pres)
pres <- data.frame(pres)
pres <- cbind(pres, PresVals)
plot(pres$ecocrop, pres$PresVals, xlab="ClimateWNA", ylab="PRISM 800m")
names(pres)[names(pres) == 'ecocrop'] <- 'EcoClimWNA'
names(pres)[names(pres) == 'PresVals'] <- 'EcoPRISM'
z <- lm(pres$EcoPRISM ~ pres$EcoClimWNA)
pres$PRISM_Offset <- z$residuals 

	# Check residuals 
	rbPal <- colorRampPalette(c('yellow','red','blue')) #Create a function to generate a continuous color palette
	pres$Col <- rbPal(70)[as.numeric(cut(z$residuals, breaks = c(seq(-0.34350, 0.34970, by=0.01))))]
	plot(pres$EcoClimWNA, pres$EcoPRISM, xlab="ClimateWNA", ylab="PRISM 800m", col=pres$Col)

	# Check outliers in relation to elevation & latitude
	# 3D Scatterplot with Coloring and Vertical Drop Lines
		library(scatterplot3d) 
		attach(mtcars) 
		scatterplot3d(pres$Elevation,pres$Latitude,pres$PRISM_Offset, highlight.3d=TRUE, pch=16, main="3D Scatterplot")
		par(mfrow=c(1,2))
		plot(pres$Latitude, pres$PRISM_Offset, ylab="PRISM - ClimateWNA"); abline(h=0, col="red", lty = 3)
		plot(pres$Elevation, pres$PRISM_Offset, ylab="PRISM - ClimateWNA"); abline(h=0, col="red", lty = 3)

		
# Cut off for prism
	par(mfrow=c(1,1))		
	hist(pres$EcoPRISM, xlim=c(0, 0.2), breaks=50)

# Portion of point under set threshold 
	cutoffs <- seq(0.01, 0.17, by=0.005) 
	toFill <- matrix(NA, nrow=2, ncol=length(cutoffs))
	for(i in 1:length(cutoffs)){
		Val <- pres$EcoPRISM < cutoffs[i]
		toFill[1,i]<- length(Val[Val==TRUE])/length(pres$PRESABS)
		toFill[2,i] <- length(Val[Val==TRUE])
	}
	plot(cutoffs, toFill[1,], xlab="EcoCropPRISM cutoff", ylab="portion of excluded presences ") 
	text(cutoffs, (toFill[1,] + 0.01), labels=toFill[2,], cex=0.67)

#  Check values of demography and transplant sites 
	site <- read.csv("C:/Users/DW/Desktop/temp.sept.30/data files/sites.csv")
	coordinates(site) = ~Longitude + Latitude
	proj4string(site) = CRS(prj.wgs)
	SiteVals <- extract(EcoCrop, site)
	site <- data.frame(site)
	site <- cbind(site, SiteVals)
	SitePreds <- site[1:51,c("ID1", "ID2", "Latitude", "Longitude", "Elevation", "SiteVals")]	
	# 0.055 good cutoff
		
# Convert pres back to shapefile & save 
coordinates(pres) = ~Longitude + Latitude
proj4string(pres) = CRS(prj.wgs)	
library(rgdal)
writeOGR(pres, dir.ecoRast, "pres", driver="ESRI Shapefile")


################################################
# Final threshold selection: 

# From sites with 
	# Angert & Schemske 2005 THE EVOLUTION OF SPECIES' DISTRIBUTIONS: RECIPROCAL TRANSPLANTS ACROSS THE ELEVATION RANGES OF MIMULUS CARDINALIS AND M. LEWISII
		# Minimum = ecocrop 0.10 
	# Angert 2009: The niche, limits to species’ distributions, and spatiotemporal variation in demography across the elevation ranges of two monkeyflowers
		# Minimum = ecocrop 0.19 
	# Minimum from demography sites (mean=0.32)
		# Little Jameson Creek = 0.05040837
		# Mill Creek = 0.05059980
	# Matthew's transplant sites 
		# Range from 0.16 - 0.31 (warm) 
	
#
##
###
####
####
##### # FINAL ECOCROP THRESHOLD CUTOFF = 0.055 
####
###
##
#
	
	
	
