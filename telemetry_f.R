library(raster)

# input podataka ----------------------------------------------------------

# dohvaća tablice direktno s ecotone servera
# zahtjeva paket 'httr'
get_ecotone_data <- function(year, month, project, username, password) {
  library(httr)
  url = paste0('http://telemetry.ecotone.pl/', project, '/exports/positions/gps_pos_', year, sprintf("%02d", month), '.csv')
  auth = readLines("data/auth.txt")
  r <- GET(url, authenticate(username, password, type = "basic"))
  return(content(r))
}

# obrada geospat tablica --------------------------------------------------

# radi raster s određenim parametrima
gridmaker <- function(spatialobject, resolution = 500, extend = 0) {
  temp <- raster(ext = extent(spatialobject), resolution = resolution)
  return(extend(temp, extend))
}

# običnu tablicu s 'Longitude' i 'Latitude' varijablama pretvara u prostornu
mkspat <- function(tab, crs = NA) {
  coords <- cbind(tab$Longitude, tab$Latitude)
  sptab <- SpatialPointsDataFrame(data = tab, coords = coords, proj4string = wgs84)
  if(!is.na(crs)){
    sptab <- spTransform(sptab, crs)
  }
  sptab$LatProj <- coordinates(sptab)[,2]
  sptab$LonProj <- coordinates(sptab)[,1]
  return(sptab)
}

# filtriranje -------------------------------------------------------------

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

# dodavanje prostornih varijabli ------------------------------------------

# funkcija za numeriranje odlazaka s kolonije
# ulaz: dataframe sa $AtColony stupcem (i.e. output spatfilter() funkcije)
foragingtrips <- function(sp_table, minpoints = 3){
  sp_table <- sp_table %>% 
    mutate(TripID = 1 + c(0, cumsum(abs(diff(sp_table$AtColony))))) %>% 
    filter(!AtColony) %>% 
    mutate(TripID = as.numeric(as.factor(TripID)))
  
  trippoints <- sp_table %>% 
    as.tibble() %>% 
    group_by(TripID) %>% 
    summarise(TripPoints = n())
  
  sp_table <- sp_table %>% 
    left_join(trippoints, by = "TripID") %>% 
    filter(TripPoints >= minpoints)
  # izlaz: dataframe sa TripID i TripPoints stupcima
  # TripID: redni broj za "foraging trip"
  # TripPoints: ukupni broj GPS lokacija za taj foraging trip (za potrebe filtriranja)
  return(sp_table)
}

# funkcija za izračun zračne udaljenosti između susjednih GPS točaka:
distance_p2p <- function(sptable) {
  output <- rep(0, nrow(sptable))
  for (i in 2:nrow(sptable)) {
    dx <- coordinates(sptable)[i,1] - coordinates(sptable)[i-1,1]
    dy <- coordinates(sptable)[i,2] - coordinates(sptable)[i-1,2]
    output[i] <- sqrt(dx^2 + dy^2)
  }
  return(output)
}

# funkcija za izračun vremenske razlike između susjednih GPS točaka:
time_p2p <- function(sptable, units = "secs") {
  td <- difftime(sptable$GPSTime[2:nrow(sptable)], sptable$GPSTime[1:nrow(sptable)-1], units = units)
  td <- c(0, td)
  return(td)
}

# funkcija koja to kombinira i dodaje tablici koju dobije kao input:
movestats <- function(sptable) {
  sptable$Distance <- distance_p2p(sptable)
  dx_nest <- coordinates(sptable)[,1] - coordinates(nest)[1]
  dy_nest <- coordinates(sptable)[,2] - coordinates(nest)[2]
  sptable$Dist_from_nest <- sqrt(dx_nest^2 + dy_nest^2)
  sptable$Timediff <- time_p2p(sptable, units = "secs")
  sptable$Speed <- sptable$Distance / sptable$Timediff
  for (i in 2:nrow(sptable)) {
    if (sptable$GpsID[i] != sptable$GpsID[i-1]) {
      sptable$Distance[i] <- 0
      sptable$Timediff[i] <- 0
      sptable$Speed[i] <- NA
    }
  }
  return(sptable)
}

# funkcija koja  kategorizira datetime varijablu u 4 kategorije
# cutpoints je vektor 4 cijela broja koji označavaju prijelomne točke u GMT vremenskoj zoni
timecat <- function(timevector, cutpoints, categories = c("prva", "druga", "treca", "cetvrta")) {
  require(lubridate)
  output <- rep(NA, length(timevector))
  output[hour(timevector) >= cutpoints[4] | hour(timevector) < cutpoints[1]] <- categories[1]
  output[hour(timevector) >= cutpoints[1] & hour(timevector) < cutpoints[2]] <- categories[2]
  output[hour(timevector) >= cutpoints[2] & hour(timevector) < cutpoints[3]] <- categories[3]
  output[hour(timevector) >= cutpoints[3] & hour(timevector) < cutpoints[4]] <- categories[4]
  return(output)
}


# pomoćne funkcije za plotanje --------------------------------------------

# dodatak za plot funkciju koji automatski zadaje područje istraživanja
plot.p <- function(...){
  plot(..., xlim = xlims, ylim = ylims)
}

# dodatak za plot funkciju koji automatski zadaje okolicu kolonije
plot.c <- function(..., radius = 10){
  plot(..., 
       xlim = c(xmin(nest) - radius * 1000, xmax(nest) + radius * 1000),
       ylim = c(ymin(nest) - radius * 1000, ymax(nest) + radius * 1000))
}