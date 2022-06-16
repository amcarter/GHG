# gas tools

# convert from atm to concentration

convert_atm_to_mol <- function(atm, temp = 293) {
  R = 0.0821 # atm L / mol K
  conc = atm / (R * temp)
  return(conc)
}

HCO2<-0.00033 #Sander 2015, NEON guide
HN2O<-0.00024 #Sander 2015
HCH4<-0.00014 #Sander 2015


#calculate Henry's law constants and error (mol m-3 Pa-1)
#no longer have error in this because the temperatures dont have error
concs$HCO2_adj<-HCO2*exp(2400*((1/(concs$storageWaterTemp+273.15))-(1/298.15)))
mean(concs$HCO2_adj, na.rm=TRUE)

concs$HN2O_adj<-HN2O*exp(2700*((1/(concs$storageWaterTemp+273.15))-(1/298.15)))

concs$HCH4_adj<-HCH4*exp(1900*((1/(concs$storageWaterTemp+273.15))-(1/298.15)))

