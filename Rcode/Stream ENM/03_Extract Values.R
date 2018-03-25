######################################################
# Snap occurrence points and pseudo absence datasets to stream lines and extract values 
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
path.dat <- paste0(path.root, "/Model_Data")
path.streamRast <- "E:/stream_gis_all"
path.exDet <- paste0(path.root, "/ExDet")

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
		coordinates(myf) = ~long + lat
		proj4string(myf) = CRS(projection(prj.wgs))
		assign(SetName[i], myf)
	}
 head(DataPseudo2, 2)

 #########################
 # load stream variables individually and extract values, then remove 
 setwd(path.streamRast)
list.files(pattern = "\\.tif$")

streamFiles <- c("ArbolateSu.tif", "bio12.tif", "bio13.tif", "bio14.tif", "bio15.tif", "bio16.tif", "bio17.tif",
	"DrainAre.tif", "SLOPE.tif", "stream.tif", "strem_el.tif", "terrough20C.tif", "TotDaSqKM.tif")  
streamNames <- c("ArbolateSu", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
	"DrainAre", "SLOPE", "stream", "strem_el", "terrough20C", "TotDaSqKM")  

#GET ALL VALUES FOR STREAM DATA
	for(i in 1:length(SetName)){
		thisSet <- get(SetName[i])
		for(j in 1:length(streamNames)){
			setwd(path.streamRast)
			thisRast <- raster(paste0(streamFiles[j]))
			newVals <- extract(thisRast, thisSet)
			assign(streamNames[j], newVals)
			thisSet <- data.frame(thisSet)
			thisSet[,j+6] <- newVals
			coordinates(thisSet) = ~long + lat
			proj4string(thisSet) = CRS(projection(prj.wgs))
		}
		names(thisSet) <- c("name", "meta", "group", "PRESABS", streamNames)
		assign(SetName[i], thisSet)
		#rm(thiSet)
	}
	
	setwd(path.dat)
	write.csv(DataPseudo1, "DataPseudo1.csv")
	write.csv(DataPseudo2, "DataPseudo2.csv")
	write.csv(DataPseudo3, "DataPseudo3.csv")
	write.csv(DataPseudo4, "DataPseudo4.csv")

	
######################################################
# From exDet background layer extract values from stream/physical rasters
	path.exDet <- paste0(path.root, "/ExDet"); setwd(path.exDet); head(dir())
	exDetpts <- read.csv("BG_for_ClimWNA.csv")
	head(exDetpts, 2)
	#exDetpts
	coordinates(exDetpts) = ~long + lat
	proj4string(exDetpts) = CRS(projection(prj.wgs))
	saveName <- names(exDetpts)
	
	for(j in 1:length(streamNames)){
			setwd(path.streamRast)
			thisRast <- raster(paste0(streamFiles[j]))
			newVals <- extract(thisRast, exDetpts)
			assign(streamNames[j], newVals)
			exDetpts <- data.frame(exDetpts)
			exDetpts[,j+5] <- newVals
			coordinates(exDetpts) = ~long + lat
			proj4string(exDetpts) = CRS(projection(prj.wgs))
		}
	names(exDetpts) <- c(saveName, streamNames)
	exDetpts <- data.frame(exDetpts)
	head(exDetpts)
	summary(exDetpts)

	setwd(path.exDet)
	write.csv(exDetpts, "BG_withStreams.csv")

	# quick check 
	path.code <- paste0(path.root, "/R_code")
	setwd(path.code)
	source('Xpairs.R')
	Xpairs(exDetpts[,c("ArbolateSu", "bio12", "bio13", "bio14", "bio15",
		"bio16", "bio17", "DrainAre", "SLOPE", "strem_el","terrough20C", "TotDaSqKM")])

	Xpairs(exDetpts[,c("ArbolateSu", "bio12", "bio13", "bio14", "bio15",
		"bio16", "bio17", "DrainAre", "SLOPE", "strem_el","terrough20C", "TotDaSqKM")])

######################################################
	
	