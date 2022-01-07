#Script to calculate slopes - then import layer to Arc and pull the slopes right 
# at/upstream of each sampling point (so slope for 180m upstream)


library(tiff)
library(raster)
library(rgdal)
library(ggplot2)
#install.packages("whitebox", repos="http://R-Forge.R-project.org")
library(whitebox)

setwd('C://Users/Alice Carter/git/ghg_patterns_nhc/src/Amanda_code/')
#nhc<-raster('durham.asc')
#nhc_tif<-writeRaster(nhc, 'nhc_tif.tif', format='GTiff')
nhc<-raster('fivecounties.tif')

#they were 30m resolution grid cells
nhc
plot(nhc)
conversion<-function(x) {mean(x)/180}
nhc_2<-aggregate(nhc, 6, FUN=conversion)
writeRaster(nhc_2, './nhc_agg.tif', format = 'GTiff', overwrite = TRUE)
#values don't seem to be changing


nhc_manip<-nhc
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
