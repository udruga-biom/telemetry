library(raster)

gridmaker <- function(spatialobject, resolution = 500, extend = 0) {
  temp <- raster(ext = extent(spatialobject), resolution = resolution)
  return(extend(temp, extend))
}

mkspat <- function(tab, crs = NA) {
  coords <- cbind(tab$Longtitude, tab$Latitude)
  sptab <- SpatialPointsDataFrame(data = tab, coords = coords, proj4string = wgs84)
  if(!is.na(crs)){
    sptab <- spTransform(sptab, crs)
  }
  return(sptab)
}
