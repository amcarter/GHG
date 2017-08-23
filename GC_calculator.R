###########################################
## Calculating Dissolved gas Concentration from Field Headspace Equilibration Method
## Alice Carter
## Apr 2017
##
## Using gas headspace equilibration technique from ____?!
## The default for this method is 40 ml H2O and 20 ml headspace

library(readr)
library(tidyr) # spreading data
library(dplyr) # munging data
source("GC_functions.R")

# Read in data and standards
stds <- GCstandards()
dat <- read_csv("data/Raw/GCrun_19Jul17.csv")

# Compare standard curves:
stdsAll <- standard_curve_dates(dat)

# Calculate the ppm of the headspace using the standards from the respective date
  # calc_ppm returns a dataframe with samples only and removes the columns with areas

ddat<- calc_ppm(dat, stdsAll)

# Replace any values that are below minimum detection limits with 1/2 MDL
# MDLs: CO2 - 100 ppm     (these values are the concentration of the lowest standard)
#       CH4 - 0.3 ppm     Unless otherwise specified, replace_MDL uses these default values
#       N2O - 0.1 ppm   

ddat <- replace_MDL(ddat)

# Incorporate volume of N2 added and vol H2O sampled



# Convert headspace concentration to equilibrium with dissolved gas concentration

# Incorporate pH?



# Functions:

calc_water_conc<-function(dat, atm.ppres){
  CO2.ppres<- dat$CO2.ppm*dat$airpres.mmhg/760/10^6
  CH4.ppres<- dat$CH4.ppm*dat$airpres.mmhg/760/10^6
  N2O.ppres<- dat$N2O.ppm*dat$airpres.mmhg/760/10^6
  volH2O.ml <- dat$H2O.mass.g-dat$evac.mass.g
  
  # Calculate henry's solubility coefficients as a function of temp
  # Emperical equations from Hudson et al 2004 SOP, indirectly from a spreadsheet from AHelton
  # Constants are in units of atm/mol fraction
  CH4.henry<- 1/(exp(-365.183 + 18106.7/(dat$temp.water+273.15)+49.7554*log(dat$temp.water+273.15)-.00028503*(dat$temp.water+273.15))/1.98719)
  CO2.henry<- 1/(exp(-317.658 + 17371.2/(dat$temp.water+273.15)+43.0607*log(dat$temp.water+273.15)-.00021910*(dat$temp.water+273.15))/1.98719)
  N2O.henry<- 1/(exp(-180.95 + 13205.8/(dat$temp.water+273.15)+20.0399*log(dat$temp.water+273.15)+.0238544*(dat$temp.water+273.15))/1.98719)

  # Calculate the aqueous concentration in water after equilibration.  
  # CA = (n/V) * pg / H * MW *10^3 mg / g *10^3 ug / mg.  
  # Where n / V = 55 mol / L (molar concentration of water) and MW = molecular weight in g/mol of gas

  molmass<- c(44, 16, 44)
  atm.ug <- atm.ppres*molmass/22.4*vol.head*10^6

  CO2.aq.ugL <- 55.5*CO2.ppres*44*10^6
  CO2.air.ugL <- dat$volN2.ml/volH2O.ml*CO2.ppres*44/22.4*(273.15/(dat$temp.water+273.15)*10^6)
  dat$CO2.ugL <- (volH2O.ml*CO2.aq.ugL + CO2.air.ugL*dat$volN2.ml - atm.ug[1])/volH2O.ml
  CH4.aq.ugL <- 55.5*CH4.ppres*16*10^6
  CH4.air.ugL <- dat$volN2.ml/volH2O.ml*CH4.ppres*16/22.4*(273.15/(dat$temp.water+273.15)*10^6)
  dat$CH4.ugL <- (volH2O.ml*CH4.aq.ugL + CH4.air.ugL*dat$volN2.ml - atm.ug[2])/volH2O.ml
  N2O.aq.ugL <- 55.5*N2O.ppres*44*10^6
  N2O.air.ugL <- dat$volN2.ml/volH2O.ml*N2O.ppres*44/22.4*(273.15/(dat$temp.water+273.15)*10^6)
  dat$N2O.ugL <- (volH2O.ml*N2O.aq.ugL + N2O.air.ugL*dat$volN2.ml - atm.ug[3])/volH2O.ml

  dat <- select(dat, -volN2.ml, -evac.mass.g, -H2O.mass.g, -CO2.ppm, -CH4.ppm, -N2O.ppm)
  dat
  }

# Calculations:


dat$date.time <- as.POSIXct(dat$date.time)
#dat$solar.time <- convert_UTC_to_solartime(date.time=dat$date.time, longitude=-79, time.type="mean solar")
scurves<- standard_curve(dat, standards)

dat$CO2.ppm <- (dat$CO2area-scurves["CO2","intercept"])/scurves["CO2","slope"]
dat$CH4.ppm <- (dat$CH4area-scurves["CH4","intercept"])/scurves["CH4","slope"]
dat$N2O.ppm <- (dat$N2Oarea-scurves["N2O","intercept"])/scurves["N2O","slope"]
dat <-select(dat, -CO2area, -CH4area, -N2Oarea )



atm.ppres<-calc_atm_ppres(dat)
atm.ppres <- c(0,0,0)

dat<- calc_water_conc(dat, atm.ppres)

# dielsamples <- dat%>% filter(type=="water")%>%select(dat$date.time, dat$CO2.ugL, dat$CH4.ugL,dat$N2O.ugL )
# par(cex=1.1)
# plot(dielsamples$date.time, dielsamples$CO2.ugL/10^6, ylab = 'CO2 g/L')
# plot(dielsamples$date.time, dielsamples$CH4.ugL/10^3, ylab = "CH4 mg/L")
# plot(dielsamples$date.time, dielsamples$N2O.ugL/10^3, ylab = "N2O mg/L")
# 

sites16May <- matrix(NA,6,6)
rownames(sites16May)<- c("Eno","UEno","NHC","UNHC","Mud","Stony")
colnames(sites16May)<- c('CO2',"sd" ,"CH4","sd", "N2O","sd")
sites26May <- sites16May

tmp <- dat %>% filter( grepl("NHC", Sample, fixed = TRUE))
tmp<- tmp[1:3,]
sites26May["NHC",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("Eno", Sample, fixed = TRUE))
sites26May["Eno",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("UNHC", Sample, fixed = TRUE))
sites26May["UNHC",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("UENO", Sample, fixed = TRUE))
sites26May["UEno",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("Mud", Sample, fixed = TRUE))
sites26May["Mud",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("Stoney", Sample, fixed = TRUE))
sites26May["Stony",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )

dat<- dat16may


dat$date.time <- as.POSIXct(dat$date.time)
#dat$solar.time <- convert_UTC_to_solartime(date.time=dat$date.time, longitude=-79, time.type="mean solar")
scurves<- standard_curve(dat, standards)

dat$CO2.ppm <- (dat$CO2area-scurves["CO2","intercept"])/scurves["CO2","slope"]
dat$CH4.ppm <- (dat$CH4area-scurves["CH4","intercept"])/scurves["CH4","slope"]
dat$N2O.ppm <- (dat$N2Oarea-scurves["N2O","intercept"])/scurves["N2O","slope"]
dat <-select(dat, -CO2area, -CH4area, -N2Oarea )



atm.ppres<-calc_atm_ppres(dat)
atm.ppres <- c(0,0,0)

dat<- calc_water_conc(dat, atm.ppres)
 
tmp <- dat %>% filter( grepl("NHC", Sample, fixed = TRUE))
tmp<- tmp[1:3,]
sites16May["NHC",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("Eno", Sample, fixed = TRUE))
sites16May["Eno",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("UNHC", Sample, fixed = TRUE))
sites16May["UNHC",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("UENO", Sample, fixed = TRUE))
sites16May["UEno",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("Mud", Sample, fixed = TRUE))
sites16May["Mud",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )
tmp <- dat %>% filter( grepl("Stoney", Sample, fixed = TRUE))
sites16May["Stony",]<- c(mean(tmp$CO2.ugL),sd(tmp$CO2.ugL),mean(tmp$CH4.ugL),sd(tmp$CH4.ugL),mean(tmp$N2O.ugL),sd(tmp$N2O.ugL) )

       
par(mfrow = c(3,1))  
barplot(sites26May[,1]/10^6, ylab = "CO2 g/L")
barplot(sites26May[,3]/10^3, ylab = "CH4 mg/L")
barplot(sites26May[,5]/10^3, ylab = "N2O mg/L")

par(mfrow = c(3,1))  
barplot(sites16May[,1]/10^6, ylab = "CO2 g/L")
barplot(sites16May[,3]/10^3, ylab = "CH4 mg/L")
barplot(sites16May[,5]/10^3, ylab = "N2O mg/L")

par(mfrow = c(1,1))
barplot(height=c(sites[,"CO2"]/10^6, sites[,'CH4']/10^3, sites[,"N2O"]/10^3), 
        beside=TRUE, ylab = "conc")

barplot(bartable, beside = TRUE, legend = levels(unique(sites)))  ## plot 

#segments(sites[,1], (sites[,1]-sites[,2])/10^6,(sites[,1]+sites[,2])/10^6 )
#arrows(sites[,1], (sites[,1]-sites[,2])/10^6,sites[,1],(sites[,1]+sites[,2])/10^6 , angle=90)
                       
