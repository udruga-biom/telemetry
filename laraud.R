#
# Analiza telemetrijskih podataka - Larus audouinii
# GPS logger: Ecotone SAKER L
# mzec 2017-10
#

verzija <- "0.7"

# setup -------------------------------------------------------------------

library(tidyverse)
library(raster)
library(sp)
library(spdplyr)
library(lubridate)
library(adehabitatHR)
library(adehabitatHS)
library(alphahull)
library(rgeos)
library(httr)
library(trip)
  
source(file = 'telemetry_f.R')
source(file = 'data/laraud_extent.R')

# postavke ----------------------------------------------------------------

# minimalni broj lokacija, foragingtrips() funkcija će filtrirati kraće tripove:
opt_minp <- 3
# radijus oko gnijezda za definiranje "kolonije" (u metrima):
opt_radius <- 500
# definiranje CRS-ova
wgs84 <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
htrs96 <- crs("+proj=tmerc +lat_0=0 +lon_0=16.5 +k=0.9999 \
              +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# definiranje mjesta gniježđenja (središte kolonije)
nest <- SpatialPoints(list(longitude = 549608, latitude = 4735945), proj4string = htrs96)

# pomoćne funkcije --------------------------------------------------------

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

# priprema podataka -------------------------------------------------------

# dohvaćanje podloga
# funkcija raster::getData() dohvaća s GADM online baze državne granice:
hr <- getData(name = "GADM", country = 'HRV', level = 0, path = 'data/')
it <- getData(name = "GADM", country = 'ITA', level = 0, path = 'data/')
al <- getData(name = "GADM", country = 'ALB', level = 0, path = 'data/')
mn <- getData(name = "GADM", country = 'BIH', level = 0, path = 'data/')
bh <- getData(name = "GADM", country = 'MNE', level = 0, path = 'data/')
sl <- getData(name = "GADM", country = 'SVN', level = 0, path = 'data/')


# spajanje pojedinačnih država u jedan sloj
drzave <- rbind(hr, it, al, mn, bh, sl) %>% spTransform(CRSobj = htrs96)

# username i šifra za pristup ecotone website-u
auth <- readLines("data/auth_laraud.txt")

# funkcija get_ecotone_data() dohvaća sa ecotone servera posljednju verziju podataka za određeni mjesec i godinu:
mj05 <- get_ecotone_data(2017, 5, "crogull", auth[1], auth[2])
mj06 <- get_ecotone_data(2017, 6, "crogull", auth[1], auth[2])
mj07 <- get_ecotone_data(2017, 7, "crogull", auth[1], auth[2])
mj08 <- get_ecotone_data(2017, 8, "crogull", auth[1], auth[2])

tab <- bind_rows(mj05, mj06, mj07, mj08)
tab$GpsID <- as.factor(tab$GpsDescription)

# definiranje spola
sextab <- data.frame(GpsID = c("CROG01", "CROG02", "CROG03", "CROG04", "CROG05"), 
                     Sex = as.factor(c("Female", "Female", "Male", "Male", "Female")))


# mkspat() fja od obične tablice radi prostornu tablicu, tj. spatialpoints dataframe
# movestats() fja dodaje point-to-point udaljenosti, vremensku razliku i brzinu
sptab <- mkspat(tab, crs = htrs96) %>% 
  movestats() %>% 
  left_join(sextab, by = "GpsID")

# dodavanje logičkog vektora za prisutnost na koloniji
sptab$AtColony <- spatfilter(input = sptab, excl_geom = nest, radius = opt_radius)
# pokazni plot
plot.c(drzave, radius = 1)
plot(sptab_g[sptab_g$AtColony,], add = TRUE)
plot(sptab_g[!sptab_g$AtColony,], add = TRUE, col = rgb(1, 0, 0))

movestats(sptab)

# definiranje perioda gniježđenja za svaku jedinku
period <- tibble(GpsID = as.factor(c("CROG01", "CROG02", "CROG03", "CROG04", "CROG05")),
                 # ovdje je definiran početak sezone, tj. za 2017. trenutak postavljanja trackera 
                 pocetak = as.POSIXct(c("2017-05-02 11:00:00",
                                        "2017-05-02 11:00:00",
                                        "2017-05-03 13:00:00",
                                        "2017-05-03 13:00:00",
                                        "2017-05-03 13:00:00"), 
                                      format = "%Y-%m-%d %H:%M:%S", 
                                      tz = "UTC"),
                 # ovdje je definiran kraj gniježđenja, tj. trenutak od kojeg se iz ponašanja ptice
                 # može zaključiti da više ne gnijezdi ("ručno" utvrđeno uvidom u podatke):
                 kraj = as.POSIXct(c("2017-06-04 03:00:00",
                                     "2017-07-09 12:00:00",
                                     "2017-06-30 10:00:00",
                                     "2017-06-08 07:30:00",
                                     "2017-06-25 07:00:00"), 
                                   format = "%Y-%m-%d %H:%M:%S", 
                                   tz = "UTC"))

# filtriranje prema periodu gniježđenja
sptab_g <- sptab %>%
  left_join(period, by = "GpsID") %>%
  filter(GPSTime > pocetak & GPSTime < kraj)

# kategoriziranje odlazaka s kolonije -------------------------------------

# numeriranje odlazaka s kolonije
g1 <- foragingtrips(filter(sptab_g, GpsID == "CROG01"), minpoints = opt_minp)
g2 <- foragingtrips(filter(sptab_g, GpsID == "CROG02"), minpoints = opt_minp)
g3 <- foragingtrips(filter(sptab_g, GpsID == "CROG03"), minpoints = opt_minp)
g4 <- foragingtrips(filter(sptab_g, GpsID == "CROG04"), minpoints = opt_minp)
g5 <- foragingtrips(filter(sptab_g, GpsID == "CROG05"), minpoints = opt_minp)

# generiranje trip objekata za sve track-ove
g1_trip <- trip(g1, TORnames = c("GPSTime", "TripID"))
g2_trip <- trip(g2, TORnames = c("GPSTime", "TripID"))
g3_trip <- trip(g3, TORnames = c("GPSTime", "TripID"))
g4_trip <- trip(g4, TORnames = c("GPSTime", "TripID"))
g5_trip <- trip(g5, TORnames = c("GPSTime", "TripID"))

# kreiranje "trajectory" objekta
g1_traj <- as.ltraj(xy = coordinates(g1), date = g1$GPSTime, id = g1$TripID)
g2_traj <- as.ltraj(xy = coordinates(g2), date = g2$GPSTime, id = g2$TripID)
g3_traj <- as.ltraj(xy = coordinates(g3), date = g3$GPSTime, id = g3$TripID)
g4_traj <- as.ltraj(xy = coordinates(g4), date = g4$GPSTime, id = g4$TripID)
g5_traj <- as.ltraj(xy = coordinates(g5), date = g5$GPSTime, id = g5$TripID)

# brzina, udaljenost ------------------------------------------------------

g <- rbind(g1, g2, g3, g4, g5)

gt <- g %>% 
  as.tibble() %>% 
  group_by(GpsID, TripID) %>% 
  summarise(pocetak = min(GPSTime), 
            kraj = max(GPSTime), 
            pocetak_mj = month(min(GPSTime), label = TRUE), 
            kraj_mj = month(max(GPSTime), label = TRUE),
            l_p2p_total = sum(Distance) / 1000, 
            l_nest_mean = mean(Dist_from_nest),
            l_nest_max = max(Dist_from_nest),
            nest_max_x = LonProj[[which.max(Dist_from_nest)]],
            nest_max_y = LatProj[[which.max(Dist_from_nest)]],
            # trajanje = sum(Timediff) / 3600, # stara metoda za izračun
            srednja_brzina = mean(Speed, na.rm = TRUE)) %>% 
  mutate(trajanje = difftime(kraj, pocetak, units = 'hours')) %>% 
  left_join(sextab, by = "GpsID")

write_csv(gt, path = 'output/laraud_trips.csv')


# jedinka, spol -----------------------------------------------------------

p <- ggplot(gt) + 
  labs(x = 'Jedinka')
  
p + geom_boxplot(aes(GpsID, l_p2p_total)) + labs(y = 'Ukupna duljina puta tijekom izlaska [m]')
p + geom_boxplot(aes(GpsID, srednja_brzina)) + labs(y = 'Srednja brzina tijekom izlaska [m/s]')
p + geom_boxplot(aes(GpsID, l_nest_mean)) + labs(y = 'Srednja zračna udaljenost od gn. tijekom izlaska [m]')
p + geom_boxplot(aes(GpsID, l_nest_max)) + labs(y = 'Maksimalna zračna udaljenost od gn. tijekom izlaska [m]')
p + geom_boxplot(aes(GpsID, as.numeric(trajanje))) + labs(y = 'Ukupno trajanje izlaska [h]')

p <- ggplot(gt) +
  labs(x = 'Spol')

p + geom_boxplot(aes(Sex, l_p2p_total)) + labs(y = 'Ukupna duljina puta tijekom izlaska [m]')
p + geom_boxplot(aes(Sex, srednja_brzina)) + labs(y = 'Srednja brzina tijekom izlaska [m/s]')
p + geom_boxplot(aes(Sex, l_nest_mean)) + labs(y = 'Srednja zračna udaljenost od gn. tijekom izlaska [m]')
p + geom_boxplot(aes(Sex, l_nest_max)) + labs(y = 'Maksimalna zračna udaljenost od gn. tijekom izlaska [m]')
p + geom_boxplot(aes(Sex, as.numeric(trajanje))) + labs(y = 'Ukupno trajanje izlaska [h]')

# kernel density estimates ------------------------------------------------

# g2points <- SpatialPoints(sptab_g2)
# g2grid <- gridmaker(sptab_g2, resolution = 500, extend = 50)
# g2gridsp <- as(g2grid, "SpatialPixels")
# 
# crs(g2grid) <- htrs96
# crs(g2points) <- htrs96
# crs(g2gridsp) <- htrs96
# 
# kda_g2 <- kernelUD(g2points, grid = g2gridsp)
# 
# kda_g2 <- kernelUD(spatfilter(g2points, excl_geom = hr_proj , radius = 0, filter = 'difference'), grid = g2gridsp)
# 
# plot(sptab_g2, pch = 16, col = rgb(0, 0, 0, 0.3))
# lines(hr_proj)
# lines(getverticeshr(kda_g2, percent = 95), col = rgb(1, 0.3, 0.3))
# lines(getverticeshr(kda_g2, percent = 90), col = rgb(0, 0.5, 0.5))
# lines(getverticeshr(kda_g2, percent = 80), col = rgb(0.2, 0.2, 1))
# lines(getverticeshr(kda_g2, percent = 75))
# 
# g3points <- SpatialPoints(sptab_g3)
# g3grid <- gridmaker(sptab_g3, resolution = 250, extend = 1000)
# g3gridsp <- as(g3grid, "SpatialPixels")
# 
# crs(g3grid) <- htrs96
# crs(g3points) <- htrs96
# crs(g3gridsp) <- htrs96
# 
# kda_g3 <- kernelUD(g3points, grid = g3gridsp)
# 
# ### potrebno izmijeniti, izmijenjen spatfilter()
# # kda_g3 <- kernelUD(spatfilter(g3points, excl_geom = hr_proj , radius = 0, filter = 'difference'), grid = g3gridsp)
# 
# plot(sptab_g3, pch = 16, col = rgb(0, 0, 0, 0.3))
# lines(hr_proj)
# lines(getverticeshr(kda_g3, percent = 95), col = rgb(1, 0.3, 0.3))
# lines(getverticeshr(kda_g3, percent = 90), col = rgb(0, 0.5, 0.5))
# lines(getverticeshr(kda_g3, percent = 80), col = rgb(0.2, 0.2, 1))
# lines(getverticeshr(kda_g3, percent = 30))
# 
# write.csv(sptab, file = "data/sptab.csv")
