# Script to calculate slopes from a 20 ft resolution DEM from:
# https://www.lib.ncsu.edu/gis/elevation
# using the whitebox tools package


library(raster)
library(sp)
library(rgdal)
library(ggplot2)
#install.packages("whitebox", repos="http://R-Forge.R-project.org")
# whitebox::install_whitebox()
library(whitebox)
library(tidyverse)

setwd('C:/Users/alice.carter/git/ghg_patterns_nhc/src/Amanda_code/')
nhc <- raster('orange/orange.asc')
sites <- read_csv('../data/site_data/NHCsite_metadata.csv')

plot(nhc)
nhc@crs

# breach depressions in the DEM:
whitebox::wbt_breach_depressions(dem = 'orange/orange.asc',
                                 output = 'orange/orange_breached.asc',
                                 fill_pits = TRUE)

nhc <- raster('orange/orange_breached.asc')
plot(nhc)


# calculate flowpath slope:
whitebox::wbt_average_flowpath_slope(dem = 'orange/orange.asc',
                                     output = 'orange/orange_slope.asc')

nhc <- raster('orange/orange_slope.asc')
plot(nhc)
points(sites$latitude, sites$longitude)
nhc_2<-aggregate(nhc, 6, FUN=conversion)
writeRaster(nhc_2, './nhc_agg.tif', format = 'GTiff', overwrite = TRUE)
#values don't seem to be changing


nhc_manip < -nhc
values(nhc_manip)<-values(nhc_manip/180)
plot(nhc_manip)
nhc_manip_agg<-aggregate(nhc_manip, 6, FUN=mean)


writeRaster(nhc_manip_agg, './nhc_ma.tif', format='GTiff', overwrite=TRUE)
nhc_ma<-raster('./nhc_ma.tif')
nhc_ma

#remove breach depressions first
wbt_breach_depressions('./nhc_ma.tif', './NHC_breachdepma.tif')
bd<-raster('./NHC_breachdepma.tif')
# wbt_breach_depressions('./nhc_agg.tif', './NHC_agg_breachdepma.tif')
# bd<-raster('./NHC_agg_breachdepma.tif')

#try flow accumulation where no branching is permitted
wbt_d8_flow_accumulation('./NHC_breachdepma.tif', './NHC_d8flowaccma.tif')
d8<-raster('./NHC_d8flowaccma.tif')
plot(d8)

#wbt_flow_accumulation_full_workflow('./nhc_2.tif', './nhc_dem2.tif', './nhc_pntr2.tif',
 #                                   './nhc_accum2.tif')


wbt_extract_streams('./NHC_d8flowaccma.tif', './streamsma.tif', threshold = 0.1, zero_background = TRUE)
streams<-raster('./streamsma.tif')
plot(streams)
#better with the dropped resolution!

#this prints in degrees - dont use
#slope('./nhc_ma.tif', './slope_ma.tif', units = 'degrees')
#slope<-raster('./slope_ma.tif')
#plot(slope)
#use this instead of the DEM for 2nd deriv?

#get the D8 flow pointer
wbt_d8_pointer('./NHC_breachdepma.tif', './NHC_d8flowpntrma.tif')
d8ptr<-raster('./NHC_d8flowpntrma.tif')
plot(d8ptr)


wbt_stream_slope_continuous('./NHC_d8flowpntrma.tif', './streamsma.tif',  './nhc_ma.tif',
                            './slope_contma.tif')
ssc<-raster('./slope_contma.tif')
plot(ssc)
ssc
writeRaster(ssc, './nhc_ssc.tif', format='GTiff', overwrite=TRUE)


#input slope to slope (2nd derivative)
wbt_stream_slope_continuous('./NHC_d8flowpntrma.tif', './streamsma.tif',  './slope_contma.tif',
                            './2deriv_contma.tif')
ssc2d<-raster('./2deriv_contma.tif')
plot(ssc2d)
writeRaster(ssc2d, './nhc_2d.tif', format='GTiff', overwrite=TRUE)

ssc2d
