#########################################################
# EXDET TOOL R TUTORIAL 
# compute type 1 & type 2 novelty for sample points 
	
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

# http://stackoverflow.com/questions/11254524/omit-rows-containing-specific-column-of-na
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}		
		
# VARIABLE SET FOR HYDRO MODEL  (up to X5) 
		#  bio15, SLOPE, terrough20C, logbio12, logDrainAre
# VARIABLE SET FOR BIOCLIM (up to X8)
		# bio2, logbio3, bio4, logbio10, bio11, logbio12, logbio14, bio15
	

# IMPORT DATA AND FIX UP VARIABLES 
	setwd(path.exDet) 
	dir()
	bioclims <- read.csv(file="BG_bioclims.csv")
	hydros <- read.csv(file="BG_withStreams.csv")
	colnames(bioclims)
	colnames(hydros)
	bioclims <- bioclims[c("ID1", "ID2", "Latitude", "Longitude", "Elevation",
			"bio2", "bio3", "bio4", "bio10", "bio11", "bio12", "bio14", "bio15")]
	hydros <- hydros[c("ID1", "ID2", "lat", "long", "el",
			"bio15", "SLOPE", "terrough20C", "bio12", "DrainAre")]
# NEED TO MAKE NESSESARY TRANSFORMATIONS
	hydros$logbio12 <- log(hydros$bio12 + 1)
	hydros$logDrainAre <- log(hydros$DrainAre + 1)
	# drop untransformed columns
	drops <- c("bio12", "DrainAre")
	hydros <- hydros[,!(names(hydros) %in% drops)]
	
	bioclims$logbio3 <- log(bioclims$bio3 + 1)
	bioclims$logbio10 <- log(bioclims$bio10 + 1)
	bioclims$logbio12 <- log(bioclims$bio12 + 1)
	bioclims$logbio14 <- log(bioclims$bio14 + 1)
	drops <- c("bio3", "bio10", "bio12", "bio14")
	bioclims <- bioclims[,!(names(bioclims) %in% drops)]
	
	
########################################################
# IMPORT REFERENCE DATASET 
setwd(path.dat)
hydroref <- read.csv(file="DataPseudo4.csv")

	#missing <- hydroref[is.na(hydroref$terrough20C),]
	hydroref <- hydroref[!is.na(hydroref$terrough20C),]
	dim(missing); dim(hydroref)

	hydroref <- hydroref[c("name", "group", "lat", "long", "bio15", "SLOPE", "terrough20C", "bio12", "DrainAre")]
	hydroref$logbio12 <- log(hydroref$bio12 + 1)
	hydroref$logDrainAre <- log(hydroref$DrainAre + 1)
	# drop untransformed columns
	drops <- c("bio12", "DrainAre")
	hydroref <- hydroref[,!(names(hydroref) %in% drops)]

########################################################
# IMPORT REFERENCE DATASET 
		setwd(path.exDet)
		biosref <- read.csv(file="dat1c.csv")
		biosref <- biosref[,c("bio2", "bio4", "bio11", "bio15", "bio3", "bio10", "bio12", "bio14")] 
		names(biosref)[names(biosref)=="bio3"] <- "logbio3"
		names(biosref)[names(biosref)=="bio10"] <- "logbio10"
		names(biosref)[names(biosref)=="bio12"] <- "logbio12"
		names(biosref)[names(biosref)=="bio14"] <- "logbio14"
	
#############################
# old from tutorial 
#			library(raster) 
#			# direcotry where all asc files are stored 
#			setwd("C:/Users/DW/Desktop/Matts misc/ExDet/Reference data (southern Australia)")
#			# list of rasters 
#			ref <- c("AusBio13.asc", "AusBio14.asc", "AusBio5.asc", "AusBio6.asc")
#			pro <- c("SaBio13.asc", "SaBio14.asc", "SaBio5.asc", "SaBio6.asc")

#--------------------------------------------------------------------#		
# DATA IMPORT 
#--------------------------------------------------------------------# 
#############################
# old from tutorial 
#		    # Import the data in R using the read.asciigrid function 
#			# of the sp package. xy location is droped for now
#					(refdat <- prodat <- list()) # make empty list for refdat and prodat
#					for(i in 1:4){
#						refdat[[i]] <- read.asciigrid(fname=ref[i])@data   # make list of dataframes (with just data - no xy lat long)
#						prodat[[i]] <- read.asciigrid(fname=pro[i])@data   # make list of dataframes (with just data - no xy lat long)
#					}
#					refdat <- do.call(cbind, refdat)
#					prodat <- do.call(cbind, prodat)

# RUN FIRST FOR HYDRO DATA & SECOND FOR BIOCLIM DATA 
refdat <- hydroref[,c("bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]
prodat <- hydros[,c("bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre")]

#refdat <- biosref[,c("bio2", "bio4", "bio11", "bio15", "logbio3", "logbio10", "logbio12", "logbio14")]
#prodat <- bioclims[,c("bio2", "bio4", "bio11", "bio15", "logbio3", "logbio10", "logbio12", "logbio14")]

# bio2, bio4, bio11, bio15, logbio3, logbio10, logbio12, logbio14 
		
#--------------------------------------------------------------------#		
# COVARIANCE MATRIX 
#--------------------------------------------------------------------#
	# Calculate the average and covariance matrix of the variables 
	#  but just in the reference set. 
	ref.av  <- colMeans(refdat, na.rm=TRUE) # mean value of each environmental variable in reference matrix
	(ref.cov <- var(refdat, na.rm=TRUE)) # COVARIANCE MATRIX OF VARIABLES 

#--------------------------------------------------------------------#		
# NT2 - MAHALANOBIS DISTANCE 
#--------------------------------------------------------------------#	
	# Calculate the mahalanobis distance of each raster 
	# cell to the environmental center of the reference 
	# set for both the reference and the projection data 
	# set and calculate the ratio between the two.
	mah.ref   <- mahalanobis(x=refdat, center=ref.av, cov=ref.cov)
	mah.pro   <- mahalanobis(x=prodat, center=ref.av, cov=ref.cov)
	mah.max <- max(mah.ref[is.finite(mah.ref)])
	nt2 <- as.data.frame(mah.pro / mah.max)

#--------------------------------------------------------------------#		
# NT1 - UNIVARIATE EXTRAPOLATION
#--------------------------------------------------------------------#	
	# UD(ij) = min(	
	MaxMin <- matrix(NA, nrow=2, ncol=ncol(refdat))
	colnames(MaxMin) <- colnames(refdat); rownames(MaxMin) = c("max", "min")
	for(i in 1:ncol(refdat)){
		MaxMin[1,i]  <- max(refdat[,i], na.rm=TRUE) # max value of each variable in the reference matrix 
		MaxMin[2,i]  <- min(refdat[,i], na.rm=TRUE) # min value of each variable in the reference matrix 
	}
	UDs <- matrix(NA, nrow(prodat), ncol(prodat))
	for(j in 1:ncol(prodat)){
		for(i in 1:nrow(prodat)){
			UDs[i,j] <- (min(c(
				(prodat[i,j] - MaxMin[2,j]),(MaxMin[1,j] - prodat[i,j]), 0), na.rm=TRUE))/
				(MaxMin[1,j] - MaxMin[2,j])
			}
		}
	UDs <- data.frame(UDs)
	NT1 <- rowSums(UDs, na.rm=TRUE)
	
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	
# STORE VALUES TO PLOT 
# Create and plot the of layer NT1
		
		
		hydros$NT1 <- NT1
		library(sp)
		coordinates(hydros) = ~long + lat
		proj4string(hydros) = prj.wgs
			rbPal <- colorRampPalette(c('lightgrey','yellow','orange','red','purple')) #Create a function to generate a continuous color palette
			hydros$Col <- rbPal(100)[as.numeric(cut(hydros$NT1, breaks = c(seq(min(hydros$NT1), max(hydros$NT1), by=((max(hydros$NT1) - min(hydros$NT1))/100)))))]
		par(mfrow=c(2,2))
		plot(hydros, col=hydros$Col, pch=19, cex=0.2)	
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create and plot the of layer nt2
		hydros <- data.frame(hydros)
		names(nt2)[names(nt2)=="mah.pro/mah.max"] <- "nt2"
		nt2 <- nt2$nt2
		hydros$nt2 <- nt2
		#hist(hydros$nt2)
		hydros <- hydros[complete.cases(hydros),]

		coordinates(hydros) = ~long + lat
		proj4string(hydros) = prj.wgs
			
		rbPal <- colorRampPalette(c('lightgrey','yellow','orange','red','purple')) #Create a function to generate a continuous color palette
		hydros$Col2 <- rbPal(100)[as.numeric(cut(hydros$nt2, breaks = c(seq(min(hydros$nt2), max(hydros$nt2), by=((max(hydros$nt2) - min(hydros$nt2))/100)))))]	
		plot(hydros, col=hydros$Col2, pch=19, cex=0.2)	

		
#
##
###
####
#####
###### # THEN BIOCLIM DATA 
#####
####
###
##
#
		
refdat <- biosref[,c("bio2", "bio4", "bio11", "bio15", "logbio3", "logbio10", "logbio12", "logbio14")]
prodat <- bioclims[,c("bio2", "bio4", "bio11", "bio15", "logbio3", "logbio10", "logbio12", "logbio14")]		
	ref.av  <- colMeans(refdat, na.rm=TRUE) # mean value of each environmental variable in reference matrix
	(ref.cov <- var(refdat, na.rm=TRUE)) # COVARIANCE MATRIX OF VARIABLES 
	mah.ref   <- mahalanobis(x=refdat, center=ref.av, cov=ref.cov)
	mah.pro   <- mahalanobis(x=prodat, center=ref.av, cov=ref.cov)
	mah.max <- max(mah.ref[is.finite(mah.ref)])
	nt2 <- as.data.frame(mah.pro / mah.max)
	MaxMin <- matrix(NA, nrow=2, ncol=ncol(refdat))
	colnames(MaxMin) <- colnames(refdat); rownames(MaxMin) = c("max", "min")
	for(i in 1:ncol(refdat)){
		MaxMin[1,i]  <- max(refdat[,i], na.rm=TRUE) # max value of each variable in the reference matrix 
		MaxMin[2,i]  <- min(refdat[,i], na.rm=TRUE) # min value of each variable in the reference matrix 
	}
	UDs <- matrix(NA, nrow(prodat), ncol(prodat))
	for(j in 1:ncol(prodat)){
		for(i in 1:nrow(prodat)){
			UDs[i,j] <- (min(c(
				(prodat[i,j] - MaxMin[2,j]),(MaxMin[1,j] - prodat[i,j]), 0), na.rm=TRUE))/
				(MaxMin[1,j] - MaxMin[2,j])
			}
		}
	UDs <- data.frame(UDs)
	NT1 <- rowSums(UDs, na.rm=TRUE)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	
# Create and plot the of layer NT1
	bioclims <- data.frame(bioclims)
	bioclims$NT1 <- NT1
	
	#hist(bioclims$NT1)
	coordinates(bioclims) = ~Longitude + Latitude
	proj4string(bioclims) = prj.wgs
		
	rbPal <- colorRampPalette(c('lightgrey','yellow','orange','red','purple')) #Create a function to generate a continuous color palette
	bioclims$Col <- rbPal(100)[as.numeric(cut(bioclims$NT1, breaks = c(seq(min(bioclims$NT1), max(bioclims$NT1), by=((max(bioclims$NT1) - min(bioclims$NT1))/100)))))]
	
	plot(bioclims, col=bioclims$Col, pch=19, cex=0.2)	

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>	
# Create and plot the of layer nt2
		bioclims <- data.frame(bioclims)
		names(nt2)[names(nt2)=="mah.pro/mah.max"] <- "nt2"
		nt2 <- nt2$nt2
		bioclims$nt2 <- nt2

		#hist(bioclims$nt2)
		coordinates(bioclims) = ~Longitude + Latitude
		proj4string(bioclims) = prj.wgs
			
			rbPal <- colorRampPalette(c('lightgrey','yellow','orange','red','purple')) #Create a function to generate a continuous color palette
			bioclims$Col2 <- rbPal(100)[as.numeric(cut(bioclims$nt2, breaks = c(seq(min(bioclims$nt2), max(bioclims$nt2), by=((max(bioclims$nt2) - min(bioclims$nt2))/100)))))]
			
		plot(bioclims, col=bioclims$Col2, pch=19, cex=0.2)	

#
##
###
####
#####
###### # CONVERT TO RASTERS & PLOT VALUES  
#####
####
###
##
#			
	
par(mfrow=c(2,2))
		hist(bioclims$nt2, xlim=c(0,3))
		hist(hydros$nt2, xlim=c(0,3))
		hist(bioclims$NT1, xlim=c(-0.8, 0))
		hist(hydros$NT1, xlim=c(-0.8, 0))
			
# Import ecocrop prism raster for reference points to make spatial grid
	setwd(path.ecocrop)
	library(raster)
	Frame <- raster("EcoCropBin.tif")
	Frame
	projection(Frame) <- prj.wgs

	x.range <- as.numeric(c(xmin(Frame), xmax(Frame)))  # min/max longitude of the interpolation area
	y.range <- as.numeric(c(ymin(Frame), ymax(Frame)))  # min/max latitude of the interpolation area
	by.x <- (xmax(Frame) - xmin(Frame))/1203
	by.y <- (ymax(Frame) - ymin(Frame))/3105

	grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = by.x), y = seq(from = y.range[1], 
    to = y.range[2], by = by.y))  # expand points to grid
	coordinates(grd) <- ~x + y
	gridded(grd) <- TRUE
	projection(grd) <- prj.wgs

	hydros <- data.frame(hydros); head(hydros, 3)
	bioclims <- data.frame(bioclims); head(bioclims, 3)
		names(hydros)[names(hydros)=="nt2"] <- "nt2hyd"
		names(hydros)[names(hydros)=="NT1"] <- "nt1hyd"
		names(bioclims)[names(bioclims)=="nt2"] <- "nt2bio"
		names(bioclims)[names(bioclims)=="NT1"] <- "nt1bio"

	toBind <- hydros[,c("ID2", "bio15", "SLOPE", "terrough20C", "logbio12", "logDrainAre", "nt1hyd", "nt2hyd", "Col2")]

# merged dataframe
	prodat <- merge(bioclims, toBind, by.x="ID2", by.y="ID2", all.x=TRUE, all.y=TRUE)
# for spatial interpolation convert to spatial pixel dataframe	
	#g <- as(Frame, 'SpatialGridDataFrame')
# interpolate points to raster
	library(ggplot2)
	library(gstat)
	library(sp)
	library(maptools)
	# drop remaining NA's in prodat
	prodat <- prodat[!is.na(prodat$nt2hyd),]
	coordinates(prodat) = ~Longitude + Latitude
	proj4string(prodat) = prj.wgs
	setwd(path.exDet)
	library(rgdal)
	writeOGR(prodat, dsn=paste0(path.exDet, "/exDet.shp"), layer="exDet", driver="ESRI Shapefile") 
	
# IDW TO NEW RASTER	
#	idw_nt2bio <- idw(formula = nt2bio ~ 1, locations = prodat, 
#			newdata = grd)  # apply idw model for the data
			
library(rgdal)
ecoregs	= readOGR(dsn="C:/Users/DW/Desktop/temp.sept.30/R objects/ecoregions.shp/us_eco_l3_no_st.shp", layer="us_eco_l3_no_st") 
states = readOGR(dsn="C:/Users/DW/Desktop/temp.sept.30/R objects/gz_2010_us_040_00_500k/gz_2010_us_040_00_500k.shp", layer="gz_2010_us_040_00_500k")
AllPres = read.csv(file="C:/Users/DW/Desktop/temp.sept.30/data files/dat4c.csv")
		coordinates(AllPres) = ~Longitude + Latitude
		proj4string(AllPres ) = prj.wgs
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = projection(ecoregs)
prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
prodat.aea = spTransform(prodat, CRS=CRS(prj.aea))
states = spTransform(states, CRS=CRS(prj.aea))
AllPres = spTransform(AllPres, CRS=CRS(prj.aea))


# Make extent object to crop polygon 
	CP <- as(extent(prodat.aea), "SpatialPolygons")
	projection(CP) = prj.aea
	library(rgeos)
	#ecoregs2 <- gIntersection(ecoregs, CP, byid=TRUE)
	ecoregs2 <- crop(ecoregs, CP)
	states <- crop(states, CP)

#COMBINE 
	moot <- over(prodat.aea, ecoregs2)
		head(moot)
		moot <- data.frame(moot)
		prodat.aea <- data.frame(prodat.aea)
		binded <- cbind(prodat.aea, moot)

library(plyr)
cdata <- ddply(binded, "NA_L3CODE", summarise,
               nt1bio    = mean(nt1bio),
               nt2bio    = mean(nt2bio),
               nt1hyd    = mean(nt1hyd),
               nt2hyd    = mean(nt2hyd)
)

 
newobj <- merge(ecoregs2, cdata, by.x="NA_L3CODE", by.y="NA_L3CODE", all.x=TRUE, all.y=FALSE)
writeOGR(newobj, dsn=paste0(path.exDet, "/ecoreg.shp"), layer="ecoreg", driver="ESRI Shapefile") # save
writeOGR(states, dsn=paste0(path.exDet, "/states.shp"), layer="states", driver="ESRI Shapefile") # save

#########################################################################
#########################################################################
#########################################################################
#########################################################################
# MAKE EXDET PLOTS 
	setwd(path.code)
	source("legend.col.R")
	par(mfrow=c(2,2))
	setwd(path.fig)
	
pdf(file = "05_ExDet_NT2NT1_reDo.pdf", width=(8.5 - (1.25 + 0.75)), height=(5), family="Times") #  - (0.75 + 0.75)
	par(ps=10, mfrow=c(5,4), mar=c(3,3,2,1), mgp = c(2, 1, 0), oma=c(0,0,0,0)) # http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
	par(mfrow=c(2,3), mar=c(1,0,1,1))

	# layout of pointd
	plot(newobj,col = "lightgrey", sub="Univariate extrapolation: Bioclimatic model", border="darkgrey")
	plot(states, add=T, border="black", lwd=1)
	plot(AllPres, add=T, pch=21, col=as.factor(AllPres$PRESABS), bg="black", cex=0.45)
	mtext("Distribution of Presence",side=1,line=0)
	mtext("and absence records",side=1,line=1)
	mtext("from testing dataset",side=1,line=2)

	
# NT1		
	rbPal <- colorRampPalette(c('red', 'yellow', 'lightgrey')) #Create a function to generate a continuous color palette
	newobj$Colnt1bio <- rbPal(100)[as.numeric(cut(newobj$nt1bio, breaks = c(seq(0, -0.2, by=-0.002))))]
	plot(newobj,col = newobj$Colnt1bio, sub="Univariate extrapolation: Bioclimatic model", border="darkgrey")
	plot(states, add=T, border="black", lwd=1)
	palette(c("darkgray", "white"))
	#plot(AllPres, add=T, pch=21, col=as.factor(AllPres$PRESABS), bg="black", cex=0.45)
	legend.col(col = rev(rbPal(50)), lev = 50:1)
	mtext("Bioclimatic ENM",side=3,line=0)
	
	newobj$Colnt1hyd <- rbPal(100)[as.numeric(cut(newobj$nt1hyd, breaks = c(seq(0, -0.2, by=-0.002))))]
	plot(newobj,col = newobj$Colnt1hyd, sub="Univariate extrapolation: Hydrological model", border="darkgrey")
	plot(states, add=T, border="black", lwd=1)
	palette(c("darkgray", "white"))
	#plot(AllPres, add=T, pch=21, col=as.factor(AllPres$PRESABS), bg="black", cex=0.45)
	mtext("Stream Habitat ENM",side=3,line=0)
	mtext("Relative Univariate Novelty",side=2,line=1)

	plot.new()
# NT2		
	#par(mfrow=c(1,2))
	rbPal <- colorRampPalette(c('lightgrey','yellow','red')) #Create a function to generate a continuous color palette
	newobj$Colnt2bio <- rbPal(100)[as.numeric(cut(newobj$nt2bio, breaks = c(seq(0, 0.70, by=0.007))))]
	plot(newobj,col = newobj$Colnt2bio, sub="Multivariate extrapolation: Bioclimatic model", border="darkgrey")
	plot(states, add=T, border="black", lwd=1)
	palette(c("darkgray", "white"))
	#plot(AllPres, add=T, pch=21, col=as.factor(AllPres$PRESABS), bg="black", cex=0.45)
	legend.col(col = rbPal(50), lev = 50:1)
	mtext("Bioclimatic ENM",side=3,line=0)

	newobj$Colnt2hyd <- rbPal(100)[as.numeric(cut(newobj$nt2hyd, breaks = c(seq(0, 0.70, by=0.007))))]
	plot(newobj,col = newobj$Colnt2hyd, sub="Multivariate extrapolation: Hydrological model", border="darkgrey")
	plot(states, add=T, border="black", lwd=1)
	palette(c("darkgray", "white"))
	#plot(AllPres, add=T, pch=21, col=as.factor(AllPres$PRESABS), bg="black", cex=0.45)
	mtext("Stream Habitat ENM",side=3,line=0)
	mtext("Relative Multivariate Novelty",side=2,line=1)

	dev.off()
	
