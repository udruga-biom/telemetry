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

# priprema podloga
hr <- getData(name = "GADM", country = 'HRV', level = 0, path = 'data/')
it <- getData(name = "GADM", country = 'ITA', level = 0, path = 'data/')
al <- getData(name = "GADM", country = 'ALB', level = 0, path = 'data/')
mn <- getData(name = "GADM", country = 'BIH', level = 0, path = 'data/')
bh <- getData(name = "GADM", country = 'MNE', level = 0, path = 'data/')
sl <- getData(name = "GADM", country = 'SVN', level = 0, path = 'data/')
drzave <- rbind(hr, it, al, mn, bh, sl)

# definiranje CRS-ova
wgs84 <- crs(hr)
htrs96 <- crs("+proj=tmerc +lat_0=0 +lon_0=16.5 +k=0.9999 +x_0=500000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

mj5 <- read_csv(file = "data/data-ichaud/gps_pos_201705.csv")
mj6 <- read_csv(file = "data/data-ichaud/gps_pos_201706.csv")
mj7 <- read_csv(file = "data/data-ichaud/gps_pos_201707.csv")
mj8 <- read_csv(file = "data/data-ichaud/gps_pos_201708.csv")

tab <- bind_rows(mj5, mj6, mj7, mj8)
tab$GpsID <- as.factor(tab$GpsNumber)
levels(tab$GpsID)[levels(tab$GpsID) == "48505471249"] <- "CROG01"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471097"] <- "CROG02"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471049"] <- "CROG03"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471215"] <- "CROG04"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471291"] <- "CROG05"

sptab <- mkspat(tab, crs = htrs96)

# definiranje perioda gniježđenja za svaku jedinku
period <- tibble(GpsID = c("CROG01", "CROG02", "CROG03", "CROG04", "CROG05"),
                 pocetak = as.POSIXct(c("2017-05-02 11:00:00",
                                        "2017-05-02 11:00:00",
                                        "2017-05-03 13:00:00",
                                        "2017-05-03 13:00:00",
                                        "2017-05-03 13:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                 kraj = as.POSIXct(c("2017-06-04 03:00:00",
                                     "2017-07-09 12:00:00",
                                     "2017-06-30 10:00:00",
                                     "2017-06-08 07:30:00",
                                     "2017-06-25 07:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

sptab_g2 <- sptab %>% filter(GpsID == "CROG02")
sptab_g3 <- sptab %>% filter(GpsID == "CROG03")


#tab

# tab %>% filter(GpsID == "CROG01") # za filtriranje po ptici

# coords <- cbind(tab_g2$Latitude, tab_g2$Longtitude)
# sptab <- SpatialPointsDataFrame(data = tab_g2, coords = coords, proj4string = wgs84)

sptab_g2 <- sptab_g2 %>% 
  filter(GPSTime > period[period$GpsID == "CROG02",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG02",]$kraj)
sptab_g3 <- sptab_g3 %>% 
  filter(GPSTime > period[period$GpsID == "CROG03",]$pocetak) %>% 
  filter(GPSTime < period[period$GpsID == "CROG03",]$kraj)

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

nest_g3 <- SpatialPoints(list(longitude = 549608, latitude = 4735945), proj4string = htrs96)

kda_g3 <- kernelUD(spatfilter(g3points, excl_geom = hr_proj , radius = 0, filter = 'difference'), grid = g3gridsp)

plot(sptab_g3, pch = 16, col = rgb(0, 0, 0, 0.3))
lines(hr_proj)
lines(getverticeshr(kda_g3, percent = 95), col = rgb(1, 0.3, 0.3))
lines(getverticeshr(kda_g3, percent = 90), col = rgb(0, 0.5, 0.5))
lines(getverticeshr(kda_g3, percent = 80), col = rgb(0.2, 0.2, 1))
lines(getverticeshr(kda_g3, percent = 30))

write.csv(sptab, file = "data/sptab.csv")
