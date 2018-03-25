###################################################################################
## Amy Angert
## Modified from Thomas C. Edwards, Jr., "Species Distribution Modelling Using R"
## PROGRAM FUNCTIONS: 
##   Create pretty map from model predictions
##   Project to region
## last update:  March 25, 2015 
###################################################################################

# load libraries now if desired; otherwise loaded below as needed
library(raster)
library(maptools)
library(gam)
library(randomForest)
library(dismo)
library(gbm)
library(rgdal)


#http://gis.stackexchange.com/questions/63577/joining-polygons
#http://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon

## set pathnames - Amy
path.root="/Users/amyangert/Desktop/Amy work August 2014" 
path.dat = paste(path.root, "/data files", sep="")
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.dat, "/ecoregions.shp", sep="")
path.bio = paste(path.dat, "/wc0.5", sep="")
path.sta = paste(path.dat, "/gz_2010_us_040_00_500k", sep="")
path.fig = paste(path.root, "/figures", sep="")

## set pathnames - Matthew
path.root = "C:/Users/DW/Desktop/temp.sept.30" 
path.dat = paste(path.root, "/data files", sep="")
path.dat.fix = paste(path.root, "/data files", sep="") # older files relocated to another directory
path.obj = paste(path.root, "/R objects", sep="")
path.eco = paste(path.obj, "/ecoregions.shp", sep="")
path.bio = paste(path.obj, "/wc0.5", sep="")
path.cod=paste(path.root, "/R code", sep="")
#path.fig=paste(path.root, "/figures", sep="")
path.sta = paste(path.obj, "/gz_2010_us_040_00_500k", sep="")


################################################################################
######## START INITIALIZATION OF DATA STRUCTURES

## read in replicate training datasets
library(sp)
library(rgdal)
setwd(path.dat) 
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

## read in occupancy dataset
library(sp)
library(rgdal)
setwd(path.dat)
all = read.csv("all.records.aug.31.csv") #includes occupancy dataset, cleaned herbarium records, and 20K pseudoabs drawn to match envir space of true absences
occ = all[all$DATASET=="occ",] #pull out occupancy dataset
dim(occ) #check size
occ$bio3 = log(occ$bio3+0.5) #make needed ln-transforms of predictors
occ$bio10 = log(occ$bio10+0.5)
occ$bio12 = log(occ$bio12+0.5)
occ$bio14 = log(occ$bio14+0.5)
coordinates(occ) = ~Longitude + Latitude
prj.wgs = "+proj=longlat +ellps=WGS84"
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
proj4string(occ) = CRS(prj.wgs)
occ.aea = spTransform(occ, CRS=CRS(prj.aea))
occ.lcc = spTransform(occ, CRS=CRS(prj.lcc))

## split by elevation 
L.thresh = 400
H.thresh = 1200
high <- occ[occ$Elevation > H.thresh,]
low <- occ[occ$Elevation < L.thresh,]
mid <- occ[occ$Elevation >= L.thresh & occ$Elevation <= H.thresh,]
high.aea = occ.aea[occ.aea$Elevation > H.thresh,]
low.aea <- occ.aea[occ.aea$Elevation < L.thresh,]
mid.aea <- occ.aea[occ.aea$Elevation >= L.thresh & occ.aea$Elevation <= H.thresh,]

high.lcc = occ.lcc[occ.lcc$Elevation > H.thresh,]
low.lcc <- occ.lcc[occ.lcc$Elevation < L.thresh,]
mid.lcc <- occ.lcc[occ.lcc$Elevation >= L.thresh & occ.lcc$Elevation <= H.thresh,]


## read in occupancy dataset split by ecoregions
setwd(path.dat)
north = read.csv("north.csv")
center = read.csv("center.csv")
south = read.csv("south.csv")
prj.wgs = "+proj=longlat +ellps=WGS84"

coordinates(north) = ~Longitude + Latitude
proj4string(north) = CRS(prj.wgs)
coordinates(center) = ~Longitude + Latitude
proj4string(center) = CRS(prj.wgs)
coordinates(south) = ~Longitude + Latitude
proj4string(south) = CRS(prj.wgs)

## reproject to albers equal area
library(rgdal)
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
occ.aea = spTransform(occ, CRS=CRS(prj.aea))
north.aea = spTransform(north, CRS=CRS(prj.aea))
center.aea = spTransform(center, CRS=CRS(prj.aea))
south.aea = spTransform(south, CRS=CRS(prj.aea))
north.lcc = spTransform(north, CRS=CRS(prj.lcc))
center.lcc = spTransform(center, CRS=CRS(prj.lcc))
south.lcc = spTransform(south, CRS=CRS(prj.lcc))

library(raster) 

###################
###################
###################
###################
###################
# 
# crop to buffered bounding box

occ2 = spTransform(occ, CRS=CRS(prj.lcc)) #ADDED
ymin = min(occ$Latitude) - 1.5
ymax = max(occ$Latitude) + 1.5
xmin = min(occ$Longitude) - 1.5
xmax = max(occ$Longitude) + 1.5
e = extent(xmin, xmax, ymin, ymax)

###################
###################
###################

setwd(path.eco)
# see 'help' for this function, you  need to set the directory to one folder level above where the destination file is saved 
ecoreg = readOGR(dsn="C:/Users/DW/Desktop/temp.sept.30/R objects/ecoregions.shp", layer="us_eco_l3_no_st") #ADDED

# project to match bioclim
library(rgdal)
prj.wgs = "+proj=longlat +ellps=WGS84"
ecoreg2 = spTransform(ecoreg, CRS(prj.lcc))
ecoreg.wgs = spTransform(ecoreg, CRS(prj.wgs))
ecoreg.aea = spTransform(ecoreg, CRS(prj.aea))


# crop polygon to buffered bounding box
# note: this loses the region names, so it doesn't work for subsetting data
library(rgeos)
bbox = as(e, "SpatialPolygons")
proj4string(bbox) = CRS(prj.wgs)
bbox.lcc = spTransform(bbox, CRS=CRS(prj.lcc))
bbox.aea = spTransform(bbox, CRS=CRS(prj.aea))
#ecoreg.crop.lcc = gIntersection(ecoreg2, bbox, byid=TRUE)
ecoreg.crop.lcc <- crop(ecoreg2, bbox.lcc) # ADDED

# pull out relevant ecoregions
eco.names = levels(ecoreg2$US_L3NAME)
poly.list = c(10,12,13,14,19,25,36,43,64,66,70,71,84)
for (i in 1:length(poly.list)) {
	poly = ecoreg2[ecoreg2$US_L3NAME==eco.names[poly.list[i]],]
	assign(paste("poly",poly.list[i], sep=""), poly)
	poly.aea = spTransform(poly, CRS=CRS(prj.aea))
	poly.lcc = spTransform(poly, CRS=CRS(prj.lcc))
	assign(paste("poly",poly.list[i],".lcc", sep=""), poly.lcc)
	assign(paste("poly",poly.list[i],".aea", sep=""), poly.aea)
	}

## world map polygons
library(rworldmap); data(countriesLow)
library(rworldxtra); data(countriesHigh)
library(rgdal)
prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
countriesLow.aea = spTransform(countriesLow, CRS=CRS(prj.aea))
countriesHigh.aea = spTransform(countriesHigh, CRS=CRS(prj.aea))
countriesLow.lcc = spTransform(countriesLow, CRS=CRS(prj.lcc))
countriesHigh.lcc = spTransform(countriesHigh, CRS=CRS(prj.lcc))

## state polygons
setwd(path.sta)
sta = readShapePoly("gz_2010_us_040_00_500k.shp")
projection(sta) = CRS(prj.wgs)
sta.aea = spTransform(sta, CRS=CRS(prj.aea))
sta.lcc = spTransform(sta, CRS=CRS(prj.lcc))

## gridlines
library(sp)
library(rgdal)#
#create unprojected gridlines
ecoreg.wgs = spTransform(ecoreg, CRS=CRS(prj.wgs))
grd.wgs = gridlines(ecoreg.wgs, ndiscr=200)
grd.wgs2 = gridlines(ecoreg.wgs, ndiscr=2000)


#project gridlines
grd.aea = spTransform(grd.wgs, CRS=CRS(prj.aea))
grd.lcc = spTransform(grd.wgs, CRS=CRS(prj.lcc))

#prepare labels for gridlines in unprojected space
gridat <- gridat(ecoreg.wgs, side="EN")
#project labels for gridlines
gridat.aea=spTransform(gridat, CRS=CRS(prj.aea), side="EN")
gridat.lcc=spTransform(gridat, CRS=CRS(prj.lcc), side="EN")

#slim down labels to fit on plots
lab.all = parse(text=as.character(gridat.aea$labels))
lab.cull = lab.all[c(2,4,6:9)]
coord.all = coordinates(gridat.aea)
coord.cull = coord.all[c(2,4,6:9),]

lab.all = parse(text=as.character(gridat.lcc$labels))
lab.cull = lab.all[c(2,4,6:9)]
coord.all = coordinates(gridat.lcc)
coord.cull = coord.all[c(2,4,6:9),]

######## END INITIALIZATION OF DATA STRUCTURES
################################################################################
  
#Clip function to trim probability rasters by bbox.lcc (extent obj)
clip<-function(raster,shape) {
          a1_crop<-crop(raster,shape)
          step1<-rasterize(shape,a1_crop)
          a1_crop*step1}

		  

#Make base plotting frame with gridlines
frame <- crop(countriesHigh.lcc, bbox.lcc); 
frame.aea <- crop(countriesHigh.aea, bbox.aea); 
frame.wgs = spTransform(frame, CRS=CRS(prj.wgs))
frame.grd = gridlines(frame.wgs, ndiscr=100)
frame.grd.lcc = spTransform(frame.grd, CRS=CRS(prj.lcc))
frame.grd.aea = spTransform(frame.grd, CRS=CRS(prj.aea))
plot(frame); plot(frame.grd.lcc, add=T) # check it out..OK?
plot(frame.aea); plot(frame.grd.aea, add=T) # check it out..OK?

## replicate herb + pseudo training sets
setwd(path.dat)
for (i in 1:10) {
	dat = read.csv(paste("dat",i,"c.csv", sep=""))
	dat = dat[,c(4:8,59:61,67:69,71:72)] #column indices changed 9/4/14
	assign(paste("dat",i, sep=""), dat)
	coordinates(dat) = ~Longitude + Latitude
	projection(dat) = prj.wgs
	dat = spTransform(dat, CRS=CRS(prj.aea))
	assign(paste("dat",i,".aea", sep=""), dat)
	}
rm(dat)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

		  
		  
		  
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
