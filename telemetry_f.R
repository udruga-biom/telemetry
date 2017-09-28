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

# prostorni filter, za isključivanje gnjezdišta itd.
spatfilter <- function(input, excl_geom, radius = 0, filter = 'intersect') {
    if(radius > 0) {
      excl_geom <- gBuffer(excl_geom, width = radius)
    }
    izvan <- input[is.na(over(input, excl_geom)),]
    unutar <- input[!is.na(over(input, excl_geom)),]
    if(filter == 'intersect') {
      return(unutar)
    } else if (filter == 'difference') {
      return(izvan)
    }
}

# filter po vremenu unutar dana, timelow i timehigh su klase 'character' u UTC vremenskoj zoni
timefilter <- function(input, timelow, timehigh, invert = FALSE) {
  input <- input %>% mutate(GPSTimelow = parse_datetime(paste0(date(GPSTime), ' ', timelow)),
                            GPSTimehigh = parse_datetime(paste0(date(GPSTime), ' ', timehigh)))
  
  if(invert == FALSE) {
    output <- input %>% 
      filter(GPSTime >= GPSTimelow & GPSTime <= GPSTimehigh)
  } else {
    output <- input %>% 
      filter(GPSTime < GPSTimelow | GPSTime > GPSTimehigh)
  }
  output <- select(output, -GPSTimelow, -GPSTimehigh)
  return(output)
}

prepend <- function(string, path) {
  fConn <- file(path, 'r+') 
  Lines <- readLines(fConn) 
  writeLines(c(string, Lines), con = fConn) 
  close(fConn)
}
