###########################################
# load in, summarize by sample/date and plot GHG data from NHC 2019-2020

# A carter
# 5/5/2020
setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/Projects/NHC Metabolism 2019-2020")

library(readr)
library(dplyr)

dat <- read_csv("../field data/2019/NHC gas data/NHC_2019-2020_processed_GHGdata.csv")

# set values below detection limits to 1/2 mdl
dat$CH4.ugL[dat$CH4.ugL<0]<- 0
dat$CO2.ugL[dat$CO2.ugL<0]<- 0
dat$N2O.ugL[dat$N2O.ugL<0]<- 0

dat.mean <- dat %>% group_by(Site, Date) %>%
  summarise(CH4.ugL.mean = mean(CH4.ugL, na.rm=T),
            CH4.ugL.sd = sd(CH4.ugL, na.rm=T),
            CO2.ugL.mean = mean(CO2.ugL, na.rm=T), 
            CO2.ugL.sd = sd(CO2.ugL, na.rm=T), 
            N2O.ugL.mean = mean(N2O.ugL, na.rm=T), 
            N2O.ugL.sd = sd(N2O.ugL, na.rm=T)
            )

 left_join(dat[,-c(7,8,9)], by=c("Site", "Date"))
