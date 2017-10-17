library(raster)

get_ecotone_data <- function(year, month, project, username, password) {
  url = paste0('http://telemetry.ecotone.pl/', project, '/exports/positions/gps_pos_', year, sprintf("%02d", month), '.csv')
  auth = readLines("data/auth.txt")
  r <- GET(url, authenticate(username, password, type = "basic"))
  return(content(r))
}

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
spatfilter <- function(input, excl_geom, radius = 0) {
    if(radius > 0) {
      excl_geom <- gBuffer(excl_geom, width = radius)
    }
    return(!is.na(over(input, excl_geom)))
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
