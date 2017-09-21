# install.packages('tidyverse')
# install.packages('ggplot2')
# install.packages('raster')
# install.packages('trip')
# install.packages('RNCEP')
# install.packages('spbabel')
# install.packages('spdplyr')
# install.packages('adehabitatHR')
# install.packages('adehabitatHS')
# install.packages('alphahull')

library(tidyverse)
library(raster)
library(sp)
library(spdplyr)
library(lubridate)
library(adehabitatHR)
library(adehabitatHS)
library(alphahull)
library(rgeos)

source('telemetry_f.R')

# priprema administrativnih podloga
hr <- getData(name = "GADM", country = 'HRV', level = 0, path = 'data/')
# it <- getData(name = "GADM", country = 'ITA', level = 0, path = 'data/')
# al <- getData(name = "GADM", country = 'ALB', level = 0, path = 'data/')
# mn <- getData(name = "GADM", country = 'BIH', level = 0, path = 'data/')
# bh <- getData(name = "GADM", country = 'MNE', level = 0, path = 'data/')
# sl <- getData(name = "GADM", country = 'SVN', level = 0, path = 'data/')
# drzave <- rbind(hr, it, al, mn, bh, sl)

# definiranje CRS-ova
wgs84 <- crs(hr)
htrs96 <- crs("+proj=tmerc +lat_0=0 +lon_0=16.5 +k=0.9999 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# reprojiciranje hr poligona u htrs96 koordinatni sustav
hr_proj <- spTransform(hr, htrs96)

# input podataka
mj4 <- read_csv(file = "data/data-aquchr/gps_pos_201704.csv")
mj5 <- read_csv(file = "data/data-aquchr/gps_pos_201705.csv")
mj6 <- read_csv(file = "data/data-aquchr/gps_pos_201706.csv")
mj7 <- read_csv(file = "data/data-aquchr/gps_pos_201707.csv")
mj8 <- read_csv(file = "data/data-aquchr/gps_pos_201708.csv")
mj9 <- read_csv(file = "data/data-aquchr/gps_pos_201709.csv")

tab <- bind_rows(mj4, mj5, mj6, mj7, mj8, mj9)
tab$GpsID <- as.factor(tab$GpsNumber)
levels(tab$GpsID)[levels(tab$GpsID) == "48505471292"] <- "CROE01"
levels(tab$GpsID)[levels(tab$GpsID) == "48505476487"] <- "CROE02"

# definiranje perioda gniježđenja za svaku jedinku
period <- tibble(GpsID = c("CROE01", "CROE02"),
                 pocetak = as.POSIXct(c("2017-04-29 06:00:00",
                                        "2017-06-29 06:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 kraj = as.POSIXct(c("2017-07-09 07:00:00",
                                     "2017-10-01 07:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

# tab %>% filter(GpsID == "CROG01") # za filtriranje po ptici

# odvojene tablice za jedinke
tab_e1 <- tab %>% filter(GpsID == "CROE01")
tab_e2 <- tab %>% filter(GpsID == "CROE02")

# generiranje prostornih dataframeova mkspat() funkcijom iz telemetry_f.R
sptab <- mkspat(tab, crs = htrs96)

sptab_e1 <- sptab %>% filter(GpsID == "CROE01")
sptab_e2 <- sptab %>% filter(GpsID == "CROE02")

# filtriranje po periodu
sptab_e1 <- sptab_e1 %>% 
  filter(GPSTime > period[period$GpsID == "CROE01",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROE01",]$kraj)
sptab_e2 <- sptab_e2 %>% 
  filter(GPSTime > period[period$GpsID == "CROE02",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROE02",]$kraj)

# output filtriranog rastera na disk
sptab_final <- rbind(sptab_e1, sptab_e2)
write.csv(sptab_final, file = "data/aquchr_sptab.csv")

# generiranje točaka bez atributa i praznih rastera za kde gridmaker() funkcijom iz telemetry_f.R
e1points <- SpatialPoints(sptab_e1)
e1grid <- gridmaker(sptab_e1, resolution = 500, extend = 50)
e1gridsp <- as(e1grid, "SpatialPixels")

e2points <- SpatialPoints(sptab_e2)
e2grid <- gridmaker(sptab_e2, resolution = 100, extend = 50)
e2gridsp <- as(e2grid, "SpatialPixels")

# unos lokacija gnjezdišta
nest_e2 <- SpatialPoints(list(longitude = 467457, latitude = 4886454), proj4string = htrs96)

##################################
#     minimum convex polygon     #
##################################

mcp_e1_95 <- mcp(e1points, percent = 95)
mcp_e1_90 <- mcp(e1points, percent = 90)
mcp_e1_80 <- mcp(e1points, percent = 80)
mcp_e1_50 <- mcp(e1points, percent = 50)

# ahull_e1_5 <- ahull(unique.array(coordinates(sptab_e1)), alpha = 5000)
# ahull_e1_25 <- ahull(unique.array(coordinates(sptab_e1)), alpha = 25000)
# ahull_e1_50 <- ahull(unique.array(coordinates(sptab_e1)), alpha = 50000)
# ahull_e1_75 <- ahull(unique.array(coordinates(sptab_e1)), alpha = 75000)

mcp_e2_95 <- mcp(SpatialPoints(sptab_e2), percent = 95)
mcp_e2_90 <- mcp(SpatialPoints(sptab_e2), percent = 90)
mcp_e2_80 <- mcp(SpatialPoints(sptab_e2), percent = 80)
mcp_e2_50 <- mcp(SpatialPoints(sptab_e2), percent = 50)

plot(sptab_e1)
lines(mcp_e1_95)
lines(mcp_e1_90)
lines(mcp_e1_80)
lines(mcp_e1_50)

plot(sptab_e2)
lines(mcp_e2_95)
lines(mcp_e2_90)
lines(mcp_e2_80)
lines(mcp_e2_50)

###################################
#     kernel density estimate     #
###################################

# izračun kernela
kda_e1 <- kernelUD(e1points, grid = e1gridsp)
kda_e2 <- kernelUD(e2points, grid = e2gridsp)

kda_e2ss <- kernelUD(spatfilter(e2points, excl_geom = nest_e2, radius = 400, filter = 'difference'), grid = e2gridsp)

# plotanje kernela
png(file = 'data/output-aquchr/croe01-kde.png', width = 800, height = 800, pointsize = 10)
plot(sptab_e1, pch = 16, col = rgb(0, 0, 0, 0.3))
lines(hr_proj)
lines(getverticeshr(kda_e1, percent = 95), col = rgb(1, 0.3, 0.3))
lines(getverticeshr(kda_e1, percent = 90), col = rgb(0, 0.5, 0.5))
lines(getverticeshr(kda_e1, percent = 80), col = rgb(0.2, 0.2, 1))
lines(getverticeshr(kda_e1, percent = 75))
dev.off()

png(file = 'data/output-aquchr/croe02-kde.png', width = 800, height = 800, pointsize = 10)
plot(sptab_e2, pch = 16, col = rgb(0, 0, 0, 0.3))
lines(hr_proj)
lines(getverticeshr(kda_e2, percent = 95), col = rgb(1, 0.3, 0.3))
lines(getverticeshr(kda_e2, percent = 90), col = rgb(0, 0.5, 0.5))
lines(getverticeshr(kda_e2, percent = 80), col = rgb(0.2, 0.2, 1))
lines(getverticeshr(kda_e2, percent = 75))
dev.off()

vert80 <- getverticeshr(kda_e1, percent = 80)



