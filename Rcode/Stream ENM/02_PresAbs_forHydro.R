######################################################
# GET PRESENCE AND GENERATE PSEUDO ABSENCE RECORDS TO USE IN HYDRO MODEL 
# - pseudo absence records are trimmed by ecocrop output. 

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
prj.wgs = "+proj=longlat +ellps=WGS84"
path.dat <- paste0(path.root, "/Model_Data")
path.ExDet <- paste0(path.root, "/ExDet")

#################################
#################################
#################################
# load ecocrop mask  & make binary 
	#ecocrop binary threshold 
	EcoThresh = 0.055
	setwd(path.ecocrop)
	EcoCrop <- raster("EcoCropOut.tif")
	par(mfrow=c(1,1)); plot(EcoCrop)
	EcoCrop
	
	m <- c(0, 0.055, 0,  0.055, 1, 1)# reclass frame
	rclmat <- matrix(m, ncol=3, byrow=TRUE)
	EcoCropBin <- reclassify(EcoCrop, rclmat) #
	# make ecocrop binary 
	writeRaster(EcoCropBin, filename="EcoCropBin.tif", datatype='INT4S', overwrite=TRUE)
	
	
#################################
#################################
# Load in adjusted cardinalis point records adjusted to stream 
	library(plotKML)
	library(rgdal)

	allherb <- readOGR(paste0(path.KMLs,'/14.04.20.all.herb.csv.kml'), "14.04.20.all.herb.csv")
	# need to get rid of all poorly geo-refer records and related errors 
			Frame <- matrix(NA, nrow=length(allherb$Name), ncol=1)
			length(allherb$Name)
			first <- strsplit(as.character(allherb$Description), "<br>PROBLEM:")	
			for(i in 1:length(allherb$Name)){
				second <- strsplit(first[[i]][2], "<br>LAT.streamstats")	
				third <- second[[1]][1]
				Frame[i,1] <- third	
			}
		unique(Frame) # view errors and issues 
	inat1 <- readOGR(paste0(path.KMLs,'/INaturalist.kml'), "INaturalist")
	inat2 <- readOGR(paste0(path.KMLs,'/INaturalistREDO.kml'), "INaturalistREDO")
	landsur <- readOGR(paste0(path.KMLs,'/LandSurveys.kml'), "LandSurveys")
	myplace <- readOGR(paste0(path.KMLs,'/My Places.kml'), "My Places")

	groups <- c("allherb", "inat1", "inat2", "landsur", "myplace")
		
	# convert to a formate where easy merging is possible 	
	for(i in 1:length(groups)){	
		this <- get(groups[i])
		myName <- this$Name 
		myGroup <- paste0(groups[i])
		myDesc <- this$Description
		mylong <- coordinates(this)[,1]	
		mylat <- coordinates(this)[,2]	
		moo <- data.frame(name = myName, long = mylong, lat = mylat, meta = myDesc)
		moo$group <- myGroup
		assign(paste0(groups[i]), moo)
		rm(moo)
		}
# need to drop bad records from allherb 
	Frame <- data.frame(Frame)
	unique(Frame)
	allherbB <- cbind(allherb, Frame)
	unique(allherbB$Frame) # all levels 
	# unacceptable levels 
	unaccepatable <- unique(allherbB$Frame)[c(2,4,5,6,8,9,10)]
	unaccepatable
	
	accepatable <- unique(allherbB$Frame)[c(1,3,7,11)]
	(accepatable  <- as.character(accepatable))
	
	allherbB <- allherbB[allherbB$Frame %in% accepatable,]
	allherbB$Frame <- as.factor(allherbB$Frame)
	unique(allherbB$Frame)
	allherb <- allherbB[,c("name", "long",  "lat",   "meta",  "group")] # overwrite existing file
	
# merge 
	allPRES <- rbind(allherb, inat1, inat2, landsur, myplace)	
	
# trim off records from SOUTH and TRANSVERSE REGION ~ will not be using them for our hydrology model 
	allPRES <- allPRES[which(allPRES$lat > 35), ]
	allPRES <- allPRES[which(allPRES$long < -10), ] # also records that didn't make it into the list or were in the middle of the atlantic (error coordinates) 

	coordinates(allPRES) = ~long + lat
		proj4string(allPRES) = CRS(projection(prj.wgs))
	plot(allPRES, pch=19, col=as.factor(allPRES$group))
	# ready to generate pseduo-absence records 
	
	
##################################################################
##################################################################	
##################################################################
##################################################################
# GENERATE PSEUDO-ABSENCE RECORD LOCATION

# Pseudo Set 1: 10\5:1 pseduoabses:presence ratio 
			  # (1.5km) minimum buffer
			  # (15km) maximum buffer
			  # exclusion area (unsuitable space from Ecocrop output).
			  n1=10
			  
# Pseudo Set 2: 100:1 pseduoabses:presence ratio 
			  # (1.5km) minimum buffer
			  # (15km) maximum buffer
			  # exclusion area (unsuitable space from Ecocrop output).
			  n2=100
			  dist1_2 <- 15000
			  
# Pseudo Set 3: 100:1 pseduoabses:presence ratio 
			  # (1.5km) minimum buffer
			  # (50km) maximum buffer
			  # exclusion area (unsuitable space from Ecocrop output).
			  n3=100
			  dist3 <- 50000
			  
# Pseudo Set 4: 15:1 pseduoabses:presence ratio 	  
			  # (1.5km) minimum buffer
			  # (50km) maximum buffer
			n4 = n1
			dist4 = dist3
		
	library(adehabitat)		  
	library(dismo)
	minCirc <- circles(allPRES, 1500, lonlat=TRUE); minCirc <- gUnaryUnion(minCirc@polygons); minCirc <- SpatialPolygonsDataFrame(minCirc,data=as.data.frame("buffer")); projection(minCirc) <- projection(prj.wgs)
	maxCirc15 <- circles(allPRES, dist1_2, lonlat=TRUE); maxCirc15 <- gUnaryUnion(maxCirc15@polygons); maxCirc15 <- SpatialPolygonsDataFrame(maxCirc15,data=as.data.frame("buffer")); projection(maxCirc15) <- projection(prj.wgs)
	maxCirc50 <- circles(allPRES, dist3, lonlat=TRUE); maxCirc50 <- gUnaryUnion(maxCirc50@polygons); maxCirc50 <- SpatialPolygonsDataFrame(maxCirc50,data=as.data.frame("buffer")); projection(maxCirc50) <- projection(prj.wgs)

	Pseu1_2 <-gDifference(maxCirc15, minCirc, byid=FALSE)
	Pseu3 <-gDifference(maxCirc50, minCirc, byid=FALSE)
	plot(Pseu1_2, col="green", bg="lightgrey")
	#writeOGR(Pseu1_2, dir.ecoRast, "Pseu1_2", driver="ESRI Shapefile")

#import stream raster set 
	stream <- raster("C:/Users/DW/Desktop/temp.sept.30/Occupancy_revisions/email_Feb23/stream_gis/stream.tif")
	
# Pseduo set 1 & 2
	test <- rasterize(Pseu1_2, EcoCropBin, field=1)
	test2 <- test + EcoCropBin
	test2[test2 < 2] = NA
	PseuoSpace1 <- test2
	PseuoSpace2 <- test2
# Pseduo set 3 & 4
	test <- rasterize(Pseu3, EcoCropBin, field=1)
	test2 <- test + EcoCropBin
	test2[test2 < 2] = NA
	PseuoSpace3 <- test2
	PseuoSpace4 <- test2

# set extent & sample points from limited area 
	myExtent <- extent(test2)
	# ~ 1hr to process (sorry)
	mainCut <- projectRaster(PseuoSpace3, stream) # takes ~ 1hr sorry 
	dir()
	writeRaster(mainCut, filename="mainCut.tif", datatype='INT4S', overwrite=TRUE)
	ToSamp <- mainCut*stream
	
	PseudoBin <- sampleRandom(ToSamp, 1200000, ext=myExtent, na.rm=TRUE, xy=TRUE, sp=FALSE)
	dim(PseudoBin) # removing NA values so in order to get a large set will need to keep repeating 

# take a giant background sample 
	setwd(path.ExDet)
	BigBin <- sampleRandom(stream, 5000000, ext=extent(stream), na.rm=TRUE, xy=TRUE, sp=TRUE)
	writeOGR(BigBin, path.ExDet, "StreamPoints", driver="ESRI Shapefile")

######################################
	PseudoBin <- data.frame(PseudoBin)
	coordinates(PseudoBin) = ~x + y
	proj4string(PseudoBin) = CRS(projection(prj.wgs))

# Pseuod1 	
	Vals1 <- extract(PseuoSpace1, PseudoBin)
	Pseudo1 <- PseudoBin
	Pseudo1$Vals <- Vals1
	Pseudo1 <- Pseudo1[which(Pseudo1$Vals == 2),]
	length(Pseudo1)
	Pseudo1 <- Pseudo1[sample(nrow(Pseudo1), length(allPRES)*n1), ]

# Pseuod2	
	Vals1 <- extract(PseuoSpace2, PseudoBin)
	Pseudo2 <- PseudoBin
	Pseudo2$Vals <- Vals1
	table(Pseudo2$Vals); length(allPRES)*n2 # will be cut short
	Pseudo2 <- Pseudo2[which(Pseudo2$Vals == 2),]
	length(Pseudo2); length(allPRES)*n2

# Pseuod3	
	Vals1 <- extract(PseuoSpace3, PseudoBin)
	Pseudo3 <- PseudoBin
	Pseudo3$Vals <- Vals1
	table(Pseudo3$Vals); length(allPRES)*n3 # will be cut short
	Pseudo3 <- Pseudo3[which(Pseudo3$Vals == 2),]
	Pseudo3 <- Pseudo3[sample(nrow(Pseudo3), length(allPRES)*n3), ]
	length(Pseudo3); length(allPRES)*n3

# Pseuod4	
	Vals1 <- extract(PseuoSpace4, PseudoBin)
	Pseudo4 <- PseudoBin
	Pseudo4$Vals <- Vals1
	table(Pseudo4$Vals); length(allPRES)*n4 # will be cut short
	Pseudo4 <- Pseudo4[which(Pseudo4$Vals == 2),]
	length(Pseudo4); length(allPRES)*n4
	Pseudo4 <- Pseudo4[sample(nrow(Pseudo4), length(allPRES)*n4), ]
	length(Pseudo4); length(allPRES)*n4
		
# View and reformat	
	plot(Pseudo3, pch=19, col="red", cex=0.15)
	plot(Pseudo2, pch=1, col="blue", cex=0.15, add=T)
	plot(Pseudo1, pch=3, col="yellow", cex=0.15, add=T)
	plot(Pseudo4, pch=3, col="green", cex=0.15, add=T)
		
length(Pseudo1); length(Pseudo2); length(Pseudo3); length(Pseudo4)
	
	Pseudo1 <- data.frame(Pseudo1); Pseudo1 <- Pseudo1[,c("x", "y")]; Pseudo1$name <- paste0("Pseudo", row(Pseudo1)[,1]); Pseudo1$meta <- "Set_Pseudo1" ; Pseudo1$group <- "Set_Pseudo1"; names(Pseudo1)[names(Pseudo1) == 'x'] <- 'long'; names(Pseudo1)[names(Pseudo1) == 'y'] <- 'lat'; Pseudo1$PRESABS <- 0; Pseudo1 <- Pseudo1[,c("name",  "long",  "lat",   "meta",  "group", "PRESABS")];
	Pseudo2 <- data.frame(Pseudo2); Pseudo2 <- Pseudo2[,c("x", "y")]; Pseudo2$name <- paste0("Pseudo", row(Pseudo2)[,1]); Pseudo2$meta <- "Set_Pseudo2" ; Pseudo2$group <- "Set_Pseudo2"; names(Pseudo2)[names(Pseudo2) == 'x'] <- 'long'; names(Pseudo2)[names(Pseudo2) == 'y'] <- 'lat'; Pseudo2$PRESABS <- 0; Pseudo2 <- Pseudo2[,c("name",  "long",  "lat",   "meta",  "group", "PRESABS")];
	Pseudo3 <- data.frame(Pseudo3); Pseudo3 <- Pseudo3[,c("x", "y")]; Pseudo3$name <- paste0("Pseudo", row(Pseudo3)[,1]); Pseudo3$meta <- "Set_Pseudo3" ; Pseudo3$group <- "Set_Pseudo3"; names(Pseudo3)[names(Pseudo3) == 'x'] <- 'long'; names(Pseudo3)[names(Pseudo3) == 'y'] <- 'lat'; Pseudo3$PRESABS <- 0; Pseudo3 <- Pseudo3[,c("name",  "long",  "lat",   "meta",  "group", "PRESABS")];
	Pseudo4 <- data.frame(Pseudo4); Pseudo4 <- Pseudo4[,c("x", "y")]; Pseudo4$name <- paste0("Pseudo", row(Pseudo4)[,1]); Pseudo4$meta <- "Set_Pseudo4" ; Pseudo4$group <- "Set_Pseudo4"; names(Pseudo4)[names(Pseudo4) == 'x'] <- 'long'; names(Pseudo4)[names(Pseudo4) == 'y'] <- 'lat'; Pseudo4$PRESABS <- 0; Pseudo4 <- Pseudo4[,c("name",  "long",  "lat",   "meta",  "group", "PRESABS")];

	head(Pseudo1, 2); head(Pseudo2, 2); head(Pseudo3, 2); head(Pseudo4, 2)
	
	
	moot <- extract(stream, allPRES)
	allPRES$moot <- moot
	setwd(path.ecocrop)
	writeOGR(allPRES, path.ecocrop, "AllpresSNAP", driver="ESRI Shapefile")
#
##
###
####
#####
###### # need to snap the few presence records manually only stream layer 
#####
####
###
#
#
	allPRESsnap <- readOGR(paste0(path.ecocrop, "/AllpresSNAP.shp"), layer="AllpresSNAP")
	moot <- extract(stream, allPRESsnap)
	allPRESsnap$moot <- moot
	allPRES <- allPRESsnap # overwrite 
		
	allPRES <- 	data.frame(allPRES)
	allPRES$PRESABS <- 1
	colnames(allPRES)

	allPRES <- allPRES[,c("name", "coords.x2", "coords.x1", "meta", "group", "PRESABS")]
	names(allPRES)[names(allPRES) == 'coords.x1'] <- "long"
	names(allPRES)[names(allPRES) == 'coords.x2'] <- "lat"
	head(allPRES)
	colnames(allPRES)
	
	DataPseudo1 <- rbind(allPRES, Pseudo1)	
	DataPseudo2 <- rbind(allPRES, Pseudo2)	
	DataPseudo3 <- rbind(allPRES, Pseudo3)	
	DataPseudo4 <- rbind(allPRES, Pseudo4)	

	# attach to allPRES and export dataframes
	path.dat <- paste0(path.root, "/Model_Data")
	setwd(path.dat)
	
	write.csv(DataPseudo1, file="DataPseudo1.csv", row.names=FALSE)
	write.csv(DataPseudo2, file="DataPseudo2.csv", row.names=FALSE)
	write.csv(DataPseudo3, file="DataPseudo3.csv", row.names=FALSE)
	write.csv(DataPseudo4, file="DataPseudo4.csv", row.names=FALSE)
	
	