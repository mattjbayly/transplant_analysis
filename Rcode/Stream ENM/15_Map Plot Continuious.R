
# Created: March 24, 2018
# Purpose: Produces plots on a continuious scale for stream habitat and climate ENMs of E. cardinalis.

# Overview 
# Load stream habitat ENM raster and take the mean of predicted values
# Load climate ENM raster and take the mean of predicted values
# Load additional spatial files - hillshade rasters; statelines; gridlines
# Create plots

# Options:
OPTION_redo_stream_habitat_raster_average = FALSE
OPTION_redo_aggregate_stream_habitat_raster_for_visability = FALSE
OPTION_redo_aggregate_climate_raster_for_visability = FALSE

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


## load libraries 
library(dismo)
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(sp)

# Use EPSG projection 3857 https://source.opennews.org/articles/choosing-right-map-projection/
prj_plotting <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"

#---------------------------------------------------
# Load additional plotting layers such as statelines
statelines <- readOGR(paste0(path.root, "/2018 Revisions/additional_plotting_layers/statelines.gpkg"), layer="statelines")

# Plotting extent -- 
  # Limit lower to 35 degrees and upper to 47.8
  plotting_extent <- extent(-126.3928, -115.08724, 35, 47.8)
  plotting_extent <- as(plotting_extent, "SpatialPolygons")
  projection(plotting_extent) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  plotting_extent <- spTransform(plotting_extent, prj_plotting)

# Hillshade
  hillshade <- raster(paste0(path.root,"/2018 Revisions/additional_plotting_layers/hillshade.tif"))

  
#------------------------------------------
# Load Transplant sites as points
  trans <- read.csv(paste0(path.root,"/2018 Revisions/Site Level ENM Preds.csv"))
  coordinates(trans) <- ~long+lat
  projection(trans) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  trans <- spTransform(trans, prj_plotting)
  
#-----------------------------------------
  # Load all occurence points 
  #Card Verified Recrods
  herb <- read.csv(paste0(path.root,"/Card Verified Recrods/14.04.20.all.herb.csv"))
  herb <- subset(herb, LAT.climate>35)
  coordinates(herb) <- ~LON.climate+LAT.climate
  projection(herb) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  herb <- spTransform(herb, prj_plotting)



#-----------------------------------------------------------------------
# Load stream habitat ENM raster and take the mean of predicted values
# Dont redo this if already done
if(OPTION_redo_stream_habitat_raster_average){
  mod_names <- c("LR", "GAM", "RF", "BRT", "MAX")
  rnames <- paste0(path.obj,"/", mod_names,".modprob.3.tif")
  stream_rasters <- lapply(rnames,raster)
  rstack <- stack(stream_rasters)
  # Very time consuming...
  overlay(rstack, fun=mean, filename=paste0(path.obj,"/mean_stream_pred.tif"))
}

# Pixel size too small to plot model predictions across pacific NW;
# instead aggrgate predictions to make visable at coarse scale
if(OPTION_redo_aggregate_stream_habitat_raster_for_visability){
  preds_stream <- raster(paste0(path.obj,"/mean_stream_pred.tif"))
  aggregate(preds_stream, fact=10, fun=max, expand=FALSE, na.rm=TRUE, format="GTiff",
            filename=paste0(path.obj,"/mean_stream_pred_agg.tif", overwrite=TRUE))
  
}
# Aggregate climate raster predictions to max for visibility
if(OPTION_redo_aggregate_climate_raster_for_visability){
  preds_clim <- raster(paste0(path.obj,"/mean_climate_pred.tif"))
  aggregate(preds_clim, fact=10, fun=max, expand=FALSE, na.rm=TRUE, format="GTiff",
            filename=paste0(path.obj,"/mean_climate_pred_agg.tif", overwrite=TRUE))
  
}


#--------------------------------------------------
# Load in stream and climate rasters for plotting
# Use EPSG projection 3857 https://source.opennews.org/articles/choosing-right-map-projection/
strmr <- raster(paste0(path.obj,"/mean_stream_pred_agg_3857.tif"))
climr <- raster(paste0(path.obj,"/mean_climate_pred_agg_3857.tif"))

# Crop down layers to plotting extent
strmrc <- crop(strmr, plotting_extent)
climrc <- crop(climr, plotting_extent)
# Fix zero bands from raster insets
climrc[climrc==0] <- NA
hillshade <- crop(hillshade, plotting_extent)
# Also clip statelines
statelinesc <- gIntersection(statelines, plotting_extent, byid=TRUE)

#
##
###
####
##### START PLOTS
####
###
##
#



#=======================================================
# START PLOT -- climate
#=======================================================

fname <- paste0(path.fig,"/spatial predictions/climate_enm_map.tif")
tiff(fname, width = 8, height = 12, units = "in", res=300, compression = c("none"),  bg = "transparent")


# Start base of plot with hillshade
par(mar=c(2.1,1.1,2.1,3.1), bg=NA)
plot(hillshade,
     axes=FALSE,
     box=FALSE,
     legend=FALSE,
     ext = extent(plotting_extent),
     col = "white"# create a color ramp of grey colors
)
#dev.off()


# Red color heat map for plotting -- from : http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=9
col_seq <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
colfunc <- colorRampPalette(col_seq)


#---------------------------------------
# Plot climate ENM as plot A to left
plot(climrc,
     breaks=c(seq(0, 1000, by=100)), 
     col=colfunc(10),
     add=TRUE,
     legend=FALSE)

#----------------------------------------
# Add on statelines
plot(statelinesc, col=NA, border="darkgrey", add=TRUE, lwd=1.8)

#------------------------------------------
# Add on native cardinalis records
plot(herb, pch=21, bg="white", col="black", add=TRUE, cex=0.8)

#----------------------------------------
# Add on Transplant site
plot(trans, add=TRUE, pch=1, col="black", cex=1.9, lwd=3.5)


dev.off()
# End of map
#-----------------------------------------------
# Start of legend

fname <- paste0(path.fig,"/spatial predictions/climate_enm_legend.png")
png(fname, width = 2.5, height = 7, units = "in", res=300, bg = "transparent")
par(bg=NA)


# Convert units -- stored as integer in raster
lbreaks=c(seq(0, 1000, by=200))/1000
# Legend colors
lcol=colfunc(10)

legend_image <- as.raster(matrix(rev(colfunc(10)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'ENM Prediction')
text(x=1.5, y = seq(0,1,l=length(lbreaks)), labels = lbreaks)
rasterImage(legend_image, 0, 0, 1,1)




dev.off()












#=======================================================
# START PLOT -- stream
#=======================================================

fname <- paste0(path.fig,"/spatial predictions/stream_enm_map.tif")
tiff(fname, width = 8, height = 12, units = "in", res=300, compression = c("none"), bg = "transparent")


# Start base of plot with hillshade
par(mar=c(2.1,1.1,2.1,3.1), bg=NA)
plot(hillshade,
     axes=FALSE,
     box=FALSE,
     legend=FALSE,
     ext = extent(plotting_extent),
     col = "white"# create a color ramp of grey colors
)
#dev.off()


# Red color heat map for plotting -- from : http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=9
col_seq <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
colfunc <- colorRampPalette(col_seq)


#---------------------------------------
# Plot climate ENM as plot A to left
plot(strmrc,
     breaks=c(seq(0, 0.2, by=0.02)), 
     col=colfunc(10),
     add=TRUE,
     legend=FALSE)

#----------------------------------------
# Add on statelines
plot(statelinesc, col=NA, border="darkgrey", add=TRUE, lwd=1.8)

#------------------------------------------
# Add on native cardinalis records
plot(herb, pch=21, bg="white", col="black", add=TRUE, cex=0.8)

#----------------------------------------
# Add on Transplant site
plot(trans, add=TRUE, pch=1, col="black", cex=1.9, lwd=3.5)


dev.off()
# End of map
#-----------------------------------------------
# Start of legend

fname <- paste0(path.fig,"/spatial predictions/stream_enm_legend.png")
png(fname, width = 2.5, height = 7, units = "in", res=300, bg = "transparent")

par(bg=NA)



# Convert units -- stored as integer in raster
lbreaks=c(seq(0, 0.2, by=0.04))
# Legend colors
lcol=colfunc(10)

legend_image <- as.raster(matrix(rev(colfunc(10)), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'ENM Prediction')
text(x=1.5, y = seq(0,1,l=length(lbreaks)), labels = lbreaks)
rasterImage(legend_image, 0, 0, 1,1)




dev.off()



#====================================================
# Extra legend for points
#------------------------------------------
fname <- paste0(path.fig,"/spatial predictions/point_legend.png")
png(fname, width = 5, height = 3, units = "in", res=300, bg = "transparent")
par(bg=NA)
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(1, 1, legend = c("Transplant Site", "Occurence Record"), 
       bg=c(NA, "white"),
       pch=c(1, 21),
       pt.lwd=c(3.5, 1),
       pt.cex=c(1.9, 0.8),
       col="black",
       cex=1.5,
       xjust=0.5, yjust=0.5)
dev.off()
