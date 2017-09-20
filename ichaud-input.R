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

mj5 <- read_csv(file = "data/gps_pos_201705.csv")
mj6 <- read_csv(file = "data/gps_pos_201706.csv")
mj7 <- read_csv(file = "data/gps_pos_201707.csv")
mj8 <- read_csv(file = "data/gps_pos_201708.csv")

tab <- bind_rows(mj5, mj6, mj7, mj8)
tab$GpsID <- as.factor(tab$GpsNumber)
levels(tab$GpsID)[levels(tab$GpsID) == "48505471249"] <- "CROG01"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471097"] <- "CROG02"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471049"] <- "CROG03"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471215"] <- "CROG04"
levels(tab$GpsID)[levels(tab$GpsID) == "48505471291"] <- "CROG05"

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

#tab

# tab %>% filter(GpsID == "CROG01") # za filtriranje po ptici

coords <- cbind(tab$Latitude, tab$Longtitude)
sptab <- SpatialPointsDataFrame(data = tab, coords = coords, proj4string = wgs84)

write.csv(sptab, file = "data/sptab.csv")
