#required function 'val2col' from: http://www.menugget.blogspot.de/2011/09/converting-values-to-color-levels.html
val2col<-function(z, zlim, col = heat.colors(12), breaks){
											 if(!missing(breaks)){
											  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
											 }
											 if(missing(breaks) & !missing(zlim)){
											  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
											  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
											  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
											 }
											 if(missing(breaks) & missing(zlim)){
											  zlim <- range(z, na.rm=TRUE)
											  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
											  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
											  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
											 }
											 colorlevels <- col[((as.vector(z)-breaks[1])/(range(breaks)[2]-range(breaks)[1]))*(length(breaks)-1)+1] # assign colors to heights for each point
											 colorlevels
}