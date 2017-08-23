##################################################
# GHG calculation functions and standards
# A Carter 24 July 2017
#
##################################################

# GHG standards from Bernhardt lab (C1-C5 standard mixes)
GCstandards <- function(){
  sample <- c("C1", "C2", "C3", "C4", "C5")
  CO2.ppm<- c(100,500,1000,2000,10000)
  CH4.ppm<- c(.3, 1, 5, 30, 200)
  N2O.ppm<- c(.1, .2, 1, 5, 25)
  stds <- data.frame(sample,CO2.ppm, CH4.ppm, N2O.ppm)
  stds
}

# Build a standard curve from a specified date:
standard_curve <- function(d, dte){
  tmp <- left_join(stds, d, by = "sample")
  w <- which(tmp$col.date==dte)
  tmp <- tmp[w,]
  scurves <- matrix(0,3,3)
  rownames(scurves)<- c('CO2', "CH4", "N2O")
  colnames(scurves)<- c( "intercept","slope", "r2")
  fitCO2 <- lm(tmp$CO2.area ~ tmp$CO2.ppm, na.rm=T)
  fitCH4 <- lm(tmp$CH4.area ~ tmp$CH4.ppm, na.rm=T)
  fitN2O <- lm(tmp$N2O.area ~ tmp$N2O.ppm, na.rm=T)
  scurves[1,1:2]<-fitCO2$coefficients
  scurves[2,1:2]<-fitCH4$coefficients
  scurves[3,1:2]<-fitN2O$coefficients
  scurves[1,3]<- summary(fitCO2)$r.squared
  scurves[2,3]<- summary(fitCH4)$r.squared
  scurves[3,3]<- summary(fitN2O)$r.squared
  
  
  par(mfrow = c(3,1), mar= c(2,5,1,1))
  plot(tmp$CO2.ppm, tmp$CO2.area, ylab = 'CO2 area')
  lines(tmp$CO2.ppm, scurves[1,2]*tmp$CO2.ppm + scurves[1,1])
  legend("topleft", paste0("r2 = ", round(scurves[1,3],3)), bty="n")
  plot(tmp$CH4.ppm, tmp$CH4.area, ylab = 'CH4 area')
  lines(tmp$CH4.ppm, scurves[2,2]*tmp$CH4.ppm + scurves[2,1])
  legend("topleft", paste0("r2 = ", round(scurves[2,3],3)), bty="n")
  plot(tmp$N2O.ppm, tmp$N2O.area, ylab = 'N2O area')
  lines(tmp$N2O.ppm, scurves[3,2]*tmp$N2O.ppm + scurves[3,1])
  legend("topleft", paste0("r2 = ", round(scurves[3,3],3)), bty="n")

  scurves
}

# Compare standard curves from different dates

standard_curve_dates <- function (d){
  ds <- d[which(d$type=="standard"),]
  dss <- left_join(stds, ds, by = "sample")
  dates <- unique(ds$col.date)
  stdsAll <- array(0, dim = c(3, 3, length(dates)), 
                   dimnames = list(c("CO2", "CH4", "N2O"), c("intercept","slope", "r2"), dates))
  for (i in dates){
    stdsAll[,,i] <- standard_curve(ds, i)
  }
  
  par(mfrow = c(length(dates),1))
  dcols <- rainbow(length(dates))
  g <- c("CO2", "CH4", "N2O")
  
  #CO2 plot:
    plot(dss$CO2.ppm, dss$CO2.area, main = "CO2", ylab = "area", xlab = "ppm")
    legend("topleft", dates, lty = rep(1, length(dates)), col = dcols)
    for(i in 1:3){
      abline(stdsAll["CO2", "intercept", dates[i]], stdsAll["CO2", "slope", dates[i]], col = dcols[i])
    }
  #CH4 plot:
    plot(dss$CH4.ppm, dss$CH4.area, main = "CH4", ylab = "area", xlab = "ppm")
    for(i in 1:3){
      abline(stdsAll["CH4", "intercept", dates[i]], stdsAll["CH4", "slope", dates[i]], col = dcols[i])
    }
  #N2O plot:
    plot(dss$N2O.ppm, dss$N2O.area, main = "N2O", ylab = "area", xlab = "ppm")
    for(i in 1:3){
      abline(stdsAll["N2O", "intercept", dates[i]], stdsAll["N2O", "slope", dates[i]], col = dcols[i])
    }
    stdsAll
}

# Calculate concentration of gases in ppm (normalize to standards only):
# this function uses the standards from the date matching the sample

calc_ppm <- function(d, stdsAll){
  dd <- d[which(d$type=="sample"),]
  dd$CO2.ppm <- 0
  dd$CH4.ppm <- 0
  dd$N2O.ppm <- 0
  dates <- unique(dd$col.date)
  for (i in dates){
    w<- which(dd$col.date==i)
    std <- stdsAll[,,i]
    dd$CO2.ppm[w] <- (dd$CO2.area[w]-std["CO2", "intercept"])/std["CO2","slope"]
    dd$CH4.ppm[w] <- (dd$CH4.area[w]-std["CH4", "intercept"])/std["CH4","slope"]
    dd$N2O.ppm[w] <- (dd$N2O.area[w]-std["N2O", "intercept"])/std["N2O","slope"]
  }
  dd %>% select(-CO2.area, -CH4.area, -N2O.area, -type)
}


# Replace values below minimum detection limits:
replace_MDL <- function(d, CO2 = 100, CH4=.3, N2O=.1){
  d$CO2.ppm[which(dat$CO2.ppm < CO2)]<- CO2/2
  d$CH4.ppm[which(dat$CH4.ppm < CH4)]<- CH4/2
  d$N2O.ppm[which(dat$N2O.ppm < N2O)]<- N2O/2
  d
}



calc_atm_ppres<- function(dat){
  tmp <- dat %>% filter(type=="atm") %>% select(-type)
  tmp$airpres.atm <- tmp$airpres.mmhg/760
  tmp$CO2.ppres <- tmp$CO2.ppm*tmp$airpres.atm/10^6
  tmp$CH4.ppres <- tmp$CH4.ppm*tmp$airpres.atm/10^6
  tmp$N2O.ppres <- tmp$N2O.ppm*tmp$airpres.atm/10^6
  
  atm <-c(mean(tmp$CO2.ppres), mean(tmp$CH4.ppres), mean(tmp$N2O.ppres))
  names(atm)<-c('CO2', "CH4", "N2O")
  atm
}

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


