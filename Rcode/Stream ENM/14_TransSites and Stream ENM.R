# load in lambda estimates
lams <- read.csv("C:/Users/DW/Desktop/transplant_analysis/Robjects/siteLambdas.csv")

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
		library(raster)
		library(maptools)
		library(gam)
		library(randomForest)
		library(dismo)
		library(gbm)
		library(biomod2)		
		library(raster)
		library(SDMTools)
		library(sf)

		
	setwd(path.obj)
	trans <- readOGR(paste0(path.obj,'/transSites.shp'), "transSites")
	plot(trans)
	text(trans, labels=trans$site)
	
	# Load in transplant site zones encompassing all plots
	trans_zones <- st_read(paste0(path.root, "/2018 Revisions/transplant field sites/field_site_zones.shp"))
	trans_zones.sp <- as(trans_zones, "Spatial")


	# load models 
	m_glm <- raster("LR.modprob.3.tif")
	m_gam <- raster("GAM.modprob.3.tif")
	m_rf <- raster("RF.modprob.3.tif")
	m_brt <- raster("BRT.modprob.3.tif")
	m_max <- raster("MAX.modprob.3.tif")
	mods <- c(m_glm, m_gam, m_rf, m_brt, m_max)
	modsSimp <- c("glm", "gam", "rf", "brt", "max")

	m_ext <- function(rast) {
		vals <- extract(rast, trans)
		return(vals)
	}
	# attach to tran df
	transpreds <- lapply(mods, m_ext)
	trandf <- data.frame(trans)
	trandf$glm <- transpreds[[1]]
	trandf$gam <- transpreds[[2]]
	trandf$rf <- transpreds[[3]]
	trandf$brt <- transpreds[[4]]
	trandf$max <- transpreds[[5]]
	
	# Extact values across whole transplant site
	m_ext_zones <- function(rast) {
	  vals <- extract(rast, trans_zones.sp, fun=max, df=TRUE, na.rm=TRUE)
	  vals <- vals[,2]
	  return(vals)
	}
	transpreds_zones <- lapply(mods, m_ext_zones)
	modpreds <- matrix(unlist(transpreds_zones), nrow=8, ncol=5)
	modpreds <- as.data.frame(modpreds)
	colnames(modpreds) <- paste0(modsSimp, "_zone")
	modpreds$site <- trans_zones.sp$join_name
	# Merge to above df
	
	trandf <- merge(trandf, modpreds, by.x="site", by.y="site")
	

	dev.off()
	# Compare difference between point extraction and zone extraction
	# Go with zone -- more accurate and exactly reflects planting sites
	plot(trandf$glm, trandf$glm_zone)
	plot(trandf$gam, trandf$gam_zone)
	plot(trandf$rf, trandf$rf_zone)
	plot(trandf$brt, trandf$brt_zone)
	plot(trandf$max, trandf$max_zone)
	

	# join with site lambda value estimates
	lams <- lams[,c("V1", "V2")] # only usefull info
	trandf <- merge(trandf, lams, by.x="site", by.y="V1")
	names(trandf)[names(trandf)=="V2"] <- "lam"
	
	######################################################
	#plot ENM vs. Latitude -- for point
	dev.off()
	par(mfrow=c(2,3))
	for(i in 1:length(modsSimp)){
	temp <- modsSimp[i]
	plot(trandf$latitude, trandf[,7+i], pch=19, xlab="latitude", ylab=modsSimp[i])
	text((trandf$latitude + 0.1), trandf[,7+i], labels=trandf$site)
	}
	
	######################################################
	#plot ENM vs. Latitude -- for zone
	par(mfrow=c(2,3))
	for(i in 1:length(modsSimp)){
	  temp <- modsSimp[i]
	  plot(trandf$latitude, trandf[,paste0(temp, "_zone")], pch=19, xlab="latitude", ylab=modsSimp[i], main="Zone")
	  text((trandf$latitude + 0.1), trandf[,paste0(temp, "_zone")], labels=trandf$site)
	}
	
	
	#plot ENM vs. lambda
	par(mfrow=c(2,3))
	for(i in 1:length(modsSimp)){
	temp <- modsSimp[i]
	plot(trandf[,7+i], trandf$lam, pch=19, col="blue", ylab="lambda", xlab=modsSimp[i])
	text(trandf[,7+i], (trandf$lam + 0.1), labels=trandf$site)
	}
	
	dev.off()
	#plot ENM vs. lambda -- zones
	par(mfrow=c(2,3))
	for(i in 1:length(modsSimp)){
	  temp <- modsSimp[i]
	  plot(trandf[,paste0(temp, "_zone")], trandf$lam, pch=19, col="blue", ylab="lambda", xlab=modsSimp[i])
	  text(trandf[,paste0(temp, "_zone")], (trandf$lam + 0.1), labels=trandf$site)
	}
	

	# Mean predictions across all modesl
	mean_pred <- rowMeans(trandf[,c("glm_zone", "gam_zone", "rf_zone", "brt_zone", "max_zone")])
	trandf$mean_pred <- mean_pred
	
	# Check visually
	dev.off()
	plot(trandf$mean_pred, trandf$lam, pch=19, col="blue", ylab="lambda", xlab="mean pred")
	text(trandf$mean_pred, (trandf$lam + 0.1), labels=trandf$site)
	
	
	
	#==========================================
	# EXPORT PREDICTIONS
	
	head(trandf)
	exp <- trandf[,c("site", "glm_zone", "gam_zone", "rf_zone", "brt_zone", "max_zone", "mean_pred")]
	
	colnames(exp) <- c("site", "glm_stream", "gam_stream", "rf_stream", "brt_stream", "max_stream", "mean_pred_stream")
	
	write.csv(exp, paste0(path.root, "/2018 Revisions/transplant field sites/field_site_predictions.csv"), row.names=FALSE)
	