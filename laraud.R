# install.packages('tidyverse')
# install.packages('ggplot2')
# install.packages('raster')
# install.packages('trip')
# install.packages('RNCEP')
# install.packages('spbabel')
# install.packages('spdplyr')

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

# priprema podloga
# funkcija raster::getData() dohvaća s GADM online baze državne granice:
hr <- getData(name = "GADM", country = 'HRV', level = 0, path = 'data/')
it <- getData(name = "GADM", country = 'ITA', level = 0, path = 'data/')
al <- getData(name = "GADM", country = 'ALB', level = 0, path = 'data/')
mn <- getData(name = "GADM", country = 'BIH', level = 0, path = 'data/')
bh <- getData(name = "GADM", country = 'MNE', level = 0, path = 'data/')
sl <- getData(name = "GADM", country = 'SVN', level = 0, path = 'data/')
drzave <- rbind(hr, it, al, mn, bh, sl) %>% spTransform(CRSobj = htrs96)

# definiranje CRS-ova
wgs84 <- crs(hr)
htrs96 <- crs("+proj=tmerc +lat_0=0 +lon_0=16.5 +k=0.9999 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# username i šifra za pristup ecotone website-u
auth <- readLines("data/auth_laraud.txt")

# funkcija get_ecotone_data() dohvaća sa ecotone servera posljednju verziju podataka za određeni mjesec i godinu:
mj05 <- get_ecotone_data(2017, 5, "crogull", auth[1], auth[2])
mj06 <- get_ecotone_data(2017, 6, "crogull", auth[1], auth[2])
mj07 <- get_ecotone_data(2017, 7, "crogull", auth[1], auth[2])
mj08 <- get_ecotone_data(2017, 8, "crogull", auth[1], auth[2])

tab <- bind_rows(mj05, mj06, mj07, mj08)
tab$GpsID <- as.factor(tab$GpsDescription)

# mkspat() funkcija od obične tablice radi prostornu, tj. spatialpoints dataframe:
sptab <- mkspat(tab, crs = htrs96)

# definiranje središnje točke kolonije i dodavanje logičkog vektora za prisutnost na koloniji
nest <- SpatialPoints(list(longitude = 549608, latitude = 4735945), proj4string = htrs96)
sptab$AtColony <- spatfilter(input = sptab, excl_geom = nest_g3, radius = 500)


# definiranje perioda gniježđenja za svaku jedinku
period <- tibble(GpsID = c("CROG01", "CROG02", "CROG03", "CROG04", "CROG05"),
                 # ovdje je definiran početak sezone, tj. za 2017. trenutak postavljanja trackera 
                 pocetak = as.POSIXct(c("2017-05-02 11:00:00",
                                        "2017-05-02 11:00:00",
                                        "2017-05-03 13:00:00",
                                        "2017-05-03 13:00:00",
                                        "2017-05-03 13:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 # ovdje je definiran kraj gniježđenja, tj. trenutak od kojeg se iz ponašanja ptice
                 # može zaključiti da više ne gnijezdi ("ručno" utvrđeno uvidom u podatke):
                 kraj = as.POSIXct(c("2017-06-04 03:00:00",
                                     "2017-07-09 12:00:00",
                                     "2017-06-30 10:00:00",
                                     "2017-06-08 07:30:00",
                                     "2017-06-25 07:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

sptab_g1 <- sptab %>%
  filter(GpsID == "CROG01") %>% 
  filter(GPSTime > period[period$GpsID == "CROG01",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG01",]$kraj)
sptab_g2 <- sptab %>%
  filter(GpsID == "CROG02") %>% 
  filter(GPSTime > period[period$GpsID == "CROG02",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG02",]$kraj)
sptab_g3 <- sptab %>%
  filter(GpsID == "CROG03") %>% 
  filter(GPSTime > period[period$GpsID == "CROG03",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG03",]$kraj)
sptab_g4 <- sptab %>%
  filter(GpsID == "CROG04") %>% 
  filter(GPSTime > period[period$GpsID == "CROG04",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG04",]$kraj)
sptab_g5 <- sptab %>%
  filter(GpsID == "CROG05") %>% 
  filter(GPSTime > period[period$GpsID == "CROG05",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG05",]$kraj)

# ponovno spajanje u jednu tablicu, filtriranu samo za gniježđenje
sptab_g <- rbind(sptab_g1, sptab_g2, sptab_g3, sptab_g4, sptab_g5)

# reprojiciranje državnih granica u htrs96, definiranje extenta za plot
xlims <- c(xmin(sptab_g), xmax(sptab_g))
ylims <- c(ymin(sptab_g), ymax(sptab_g))

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

# pokazni plot
plot.c(drzave, radius = 1)
plot(sptab_g[sptab_g$AtColony,], add = TRUE)
plot(sptab_g[!sptab_g$AtColony,], add = TRUE, col = rgb(1, 0, 0))

# funkcija za numeriranje odlazaka s kolonije
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
  
  return(sp_table)
}

# minimalni broj lokacija, foragingtrips() funkcija će filtrirati kraće tripove
minp <- 10

# numeriranje odlazaka s kolonije
g1 <- foragingtrips(sptab_g1, minpoints = minp)
g2 <- foragingtrips(sptab_g2, minpoints = minp)
g3 <- foragingtrips(sptab_g3, minpoints = minp)
g4 <- foragingtrips(sptab_g4, minpoints = minp)
g5 <- foragingtrips(sptab_g5, minpoints = minp)

# generiranje trip objekata za sve track-ove
g1_trip <- trip(g1, TORnames = c("GPSTime", "TripID"))
g2_trip <- trip(g2, TORnames = c("GPSTime", "TripID"))
g3_trip <- trip(g3, TORnames = c("GPSTime", "TripID"))
g4_trip <- trip(g4, TORnames = c("GPSTime", "TripID"))
g5_trip <- trip(g5, TORnames = c("GPSTime", "TripID"))

# dalje za igranje:
homedist()

# kreiranje "trajectory" objekta
g1_traj <- as.ltraj(xy = coordinates(g1), date = g1$GPSTime, id = g1$TripID)
g2_traj <- as.ltraj(xy = coordinates(g2), date = g2$GPSTime, id = g2$TripID)
g3_traj <- as.ltraj(xy = coordinates(g3), date = g3$GPSTime, id = g3$TripID)
g4_traj <- as.ltraj(xy = coordinates(g4), date = g4$GPSTime, id = g4$TripID)
g5_traj <- as.ltraj(xy = coordinates(g5), date = g5$GPSTime, id = g5$TripID)


#####
# kernel density estimates
#####

g2points <- SpatialPoints(sptab_g2)
g2grid <- gridmaker(sptab_g2, resolution = 500, extend = 50)
g2gridsp <- as(g2grid, "SpatialPixels")

crs(g2grid) <- htrs96
crs(g2points) <- htrs96
crs(g2gridsp) <- htrs96

kda_g2 <- kernelUD(g2points, grid = g2gridsp)

nest_g2 <- SpatialPoints(list(longitude = 549608, latitude = 4735945), proj4string = htrs96)

kda_g2 <- kernelUD(spatfilter(g2points, excl_geom = hr_proj , radius = 0, filter = 'difference'), grid = g2gridsp)

plot(sptab_g2, pch = 16, col = rgb(0, 0, 0, 0.3))
lines(hr_proj)
lines(getverticeshr(kda_g2, percent = 95), col = rgb(1, 0.3, 0.3))
lines(getverticeshr(kda_g2, percent = 90), col = rgb(0, 0.5, 0.5))
lines(getverticeshr(kda_g2, percent = 80), col = rgb(0.2, 0.2, 1))
lines(getverticeshr(kda_g2, percent = 75))

g3points <- SpatialPoints(sptab_g3)
g3grid <- gridmaker(sptab_g3, resolution = 250, extend = 1000)
g3gridsp <- as(g3grid, "SpatialPixels")

crs(g3grid) <- htrs96
crs(g3points) <- htrs96
crs(g3gridsp) <- htrs96

kda_g3 <- kernelUD(g3points, grid = g3gridsp)

### potrebno izmijeniti, izmijenjen spatfilter()
# kda_g3 <- kernelUD(spatfilter(g3points, excl_geom = hr_proj , radius = 0, filter = 'difference'), grid = g3gridsp)

plot(sptab_g3, pch = 16, col = rgb(0, 0, 0, 0.3))
lines(hr_proj)
lines(getverticeshr(kda_g3, percent = 95), col = rgb(1, 0.3, 0.3))
lines(getverticeshr(kda_g3, percent = 90), col = rgb(0, 0.5, 0.5))
lines(getverticeshr(kda_g3, percent = 80), col = rgb(0.2, 0.2, 1))
lines(getverticeshr(kda_g3, percent = 30))

write.csv(sptab, file = "data/sptab.csv")
