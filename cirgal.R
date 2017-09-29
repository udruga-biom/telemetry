#############################################################
# Obrada telemetrijskih podataka za objavu (Leaflet, Carto) #
# Udruga BIOM, 2017                                         #
#############################################################
ver = "0.6"

library(tidyverse)
library(raster)
library(sp)
library(spdplyr)
library(lubridate)
library(adehabitatHR)
library(adehabitatHS)
library(alphahull)
library(rgeos)
library(geojson)
library(geojsonio)
library(httr)

setwd('~/Documents/Code/telemetry/')
source('telemetry_f.R')

# priprema administrativnih podloga
hr <- getData(name = "GADM", country = 'HRV', level = 0, path = 'data/')
# drzave <- rbind(hr, it, al, mn, bh, sl)

# definiranje CRS-ova
wgs84 <- crs(hr)
htrs96 <- crs("+proj=tmerc +lat_0=0 +lon_0=16.5 +k=0.9999 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# reprojiciranje hr poligona u htrs96 koordinatni sustav
hr_proj <- spTransform(hr, htrs96)

# input podataka
# mj09 <- read_csv(file = "data/data-aquchr/gps_pos_201709.csv")
# mj10 <- read_csv(file = "data/data-aquchr/gps_pos_201710.csv")

# novi input podataka drito s ecotone web-a
auth = readLines("data/auth.txt")
mj09 <- get_ecotone_data(2017, 9, "croeagle", auth[1], auth[2])

tab <- bind_rows(mj09)

# dodavanje stringa za datum

tab$GPSTimeGMT <- tab$GPSTime
tab$GPSTime <- with_tz(tab$GPSTime, tzone = "Europe/Zagreb")
tab$Date <- paste0(sprintf("%02d", day(tab$GPSTime)), ".", sprintf("%02d", month(tab$GPSTime)), ".", year(tab$GPSTime), " ", 
                   sprintf("%02d", hour(tab$GPSTime)), ":", sprintf("%02d", minute(tab$GPSTime)))

# definiranje definiranje imena GPS-a
# tab$GpsID <- as.factor(tab$GpsNumber)
# levels(tab$GpsID)[levels(tab$GpsID) == "48505471292"] <- "CROE01"
# levels(tab$GpsID)[levels(tab$GpsID) == "48505476487"] <- "CROE02"

# definiranje perioda gniježđenja za svaku jedinku
period <- tibble(GpsDescription = c("CROE01"),
                 pocetak = as.POSIXct(c("2017-07-27 06:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Zagreb"),
                 kraj = as.POSIXct(c("2017-12-31 07:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "Europe/Zagreb"))

# tab %>% filter(GpsDescription == "CROG01") # za filtriranje po ptici

# odvojene tablice za jedinke
tab_e1 <- tab %>% filter(GpsDescription == "CROE01")

# generiranje prostornih dataframeova mkspat() funkcijom iz telemetry_f.R
sptab <- mkspat(tab, crs = wgs84)

sptab_e1 <- sptab %>% filter(GpsDescription == "CROE01")

# filtriranje po periodu
sptab_e1 <- sptab_e1 %>% 
  filter(GPSTime > period[period$GpsDescription == "CROE01",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsDescription == "CROE01",]$kraj)

# output filtriranog rastera na disk
write.csv(sptab_e1, file = "data/cirgal_sptab.csv")

output <- sptab_e1 %>% 
  transmute(gpstime = GPSTime,
            date = Date,
            latitude = Latitude,
            longitude = Longtitude,
            name = GpsDescription,
            temp = Temperature,
            smstime = SMSTime,
            battvol = BatteryVoltage,
            interval = GPSIntervals,
            gsmsignal = GSMSignal)

line <- Lines(list(Line(coordinates(output))), ID = "CROE01")
spline <- SpatialLines(list(line), proj4string = wgs84)

# geojson output
pointfile <- "/home/mzec/biom-cloud-it/telemetrija/points.json"
linefile <- "/home/mzec/biom-cloud-it/telemetrija/line.json"

output_gj <- as.geojson(output)
output_line_gj <- as.geojson(spline)
geo_write(output_gj, pointfile)
geo_write(output_line_gj, linefile)

