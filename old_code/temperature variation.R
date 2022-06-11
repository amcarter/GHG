library(tidyverse)
library(lubridate)
setwd(metab_projdir)

sites <- c( "NHC","PM","CBP", "WB", "WBP", "UNHC" )
site <- "NHC"
dat <- readRDS(paste0("data/metabolism/modeled_prior_nreg_K/",site,".rds"))
  
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c("brown3","grey60"))

png(width=7, height=5, units='in', filename="figures/effective_dischargeNHC.png", type='cairo', res=300)
par(mfrow=c(3,2), mar=c(4,4,1,1))

for(i in sites){
  
  dat <- readRDS(paste0("data/metabolism/modeled_prior_nreg_K/",i,".rds"))
  dat <- dat$data %>%
    select(date, discharge, DO.obs, DO.sat, )%>%
    mutate(DO_pctsat = DO.obs/DO.sat) %>%
    left_join(dat$metab[,c(1,5)], by="date")
  dat$color <- rbPal(10)[as.numeric(cut(dat$ER,breaks = 10))]
  Qeff<- calc_QeffvsDO(dat)
  
  plot( dat$discharge,dat$DO_pctsat,pch=20, log="x", col=alpha(dat$color, .5),main=i, xlab="effective discharge",ylab="DO pct sat")
  lines(Qeff$Qeff, Qeff$ox, lwd=2)
  #plot( dat$metab$ER,dat$metab$K600,pch=19, log="x", main=i)
  
} 
dev.off()


  dat_daily <- dat$data %>% group_by(date)%>%
    summarize(temp.max=max(temp.water, na.rm=T),
              temp.min=min(temp.water, na.rm=T))
  dat_daily$temp.amp <- dat_daily$temp.max-dat_daily$temp.min
  
  dat_daily[which(is.infinite(dat_daily$temp.amp)),2:4]<- NA
  
  plot(dat$data$solar.time, dat$data$DO.obs, type="l", xlab="date", ylab="DO", col="grey60", lwd=.8, main=i)
  par(new=T)
  plot(dat$metab$date, dat$metab$discharge.m3s, type="l", axes=FALSE, col="blue", lwd=1.5, ylab="", xlab="", log="y")
  par(new=T)
  plot(dat_daily$date, dat_daily$temp.amp, pch=19, axes=F, xlab="", ylab="")
  axis(4)
  par(new=T)
}

avg <- dat$data %>% group_by(date)%>%
  summarise(max.temp = max(temp.water, na.rm=T),
            min.temp = min(temp.water, na.rm=T),
            mean.temp=mean(temp.water, na.rm=T))
avg$amp <- avg$max.temp-avg$min.temp
avg$month <- month(avg$date)
avg[is.infinite(avg$amp),2:5]<- NA
mons <- avg %>% group_by(month)%>%
  summarise(temp.amp=mean(amp, na.rm=T),
            temp.max = mean(max.temp, na.rm=T),
            temp.min=mean(min.temp, na.rm=T),
            temp.mean=mean(mean.temp, na.rm=T))

day <- "-17"
png(width=7, height=6, units='in', filename="night_DOsat_eachMonth.png", type='cairo', res=300)
par(mar=c(2,2,0,0), oma = c(2,2,3,1), mfrow = c(4,3))
for(i in 1:12){
  if(i %in% 1:2){
    date <- as.Date(paste0("2020-0",i,day))
  }else if(i %in% 3:9){
    date <- as.Date(paste0("2019-0",i,day))
  }else{
    date <- as.Date(paste0("2019-",i,day))
  }
  start <- ymd_hms(paste(date-1, "10:00:00"))
  end<- ymd_hms(paste(date,"15:00:00"))
  tmp <- dat$data[dat$data$solar.time >=start&dat$data$solar.time<=end,]
  ylims <- c(min(c(tmp$DO.obs, tmp$DO.sat), na.rm=T)-.5, max(c(tmp$DO.obs, tmp$DO.sat), na.rm=T)+.7)
  delT <- round(max(tmp$temp.water, na.rm=T)-min(tmp$temp.water, na.rm=T),1)
  plot(tmp$solar.time, tmp$DO.obs, type="l", ylim = ylims)
  lines(tmp$solar.time, tmp$DO.sat, lty=2)
  w<- which(tmp$light==0)
  polygon(tmp$solar.time[c(w[1], w[length(w)],w[length(w)],w[1])], c(0,0,15,15), col = alpha("grey70",.5), border=NA)
  mtext(as.character(month(i,label=TRUE, abbr=TRUE)), 3, line=-1.2, adj=.1, cex=.8)
  mtext(paste("delta T =",delT), 3, line=-1.2, adj=.9, cex=.8)
  #mtext(paste("average delta T =",round(mons$temp.amp[mons$month==i],1)), 1, line=-1.2, adj=.9, cex=.8)
}


par(new=T, mfrow=c(1,1), mar=c(2,4,0,1), oma=c(0,0,0,0))
plot(1,1,type="n", ylab="DO mg/L", xlab="", axes=F)
legend("top",c("DO observed", "DO saturation"), bty="n", lty=c(1,2), ncol=2, xpd=F )
dev.off()

#annual temp patterns

png(width=7, height=2, units="in", type="cairo", res=200, filename="tempswings_byMonth_NHC.png")
par(mar=c(4,4,1,4))
plot(mons$month, mons$temp.mean, type="n", axes=F, ylab="", xlab="", ylim = c(5,30))
polygon(c(mons$month, 12,1), c(mons$temp.mean,0,0), col=alpha("darkblue",.4),border=NA)
axis(4, col="darkblue", col.axis="darkblue")
mtext("mean temp C", side=4, line=2.2, col="darkblue")
par(new=T)
plot(mons$month, mons$temp.amp, type="l", lwd=2, axes=F, xlab="", ylab = "delta temp C", xpd=T, ylim = c(0,3))
axis(2)
axis(1,labels=F )
t <- seq(2,12, by=2)
text(x = t, y = -.15 , labels = as.character(month(t, label=T, abbr=T)), 
     pos = 1, xpd=NA,  cex=.9)
dev.off()



png(width=7, height=4,units="in", type="cairo", res=200, file="figures/tempswings_monthly_allsites.png")
par(mfrow = c(3,2),mar=c(4,4,1,4))
for(s in sites) {
  dat <- readRDS(paste0("data/metabolism/modeled_prior_nreg_K/",s,".rds"))
  
  avg <- dat$data %>% group_by(date)%>%
    summarise(max.temp = max(temp.water, na.rm=T),
              min.temp = min(temp.water, na.rm=T),
              mean.temp=mean(temp.water, na.rm=T))
  avg$amp <- avg$max.temp-avg$min.temp
  avg$month <- month(avg$date)
  avg[is.infinite(avg$amp),2:5]<- NA
  mons <- avg %>% group_by(month)%>%
    summarise(temp.amp=mean(amp, na.rm=T),
              temp.max = mean(max.temp, na.rm=T),
              temp.min=mean(min.temp, na.rm=T),
              temp.mean=mean(mean.temp, na.rm=T))
  
  plot(mons$month, mons$temp.mean, type="n", axes=F, ylab="", xlab="", ylim = c(5,30))
  polygon(c(mons$month, 12,1), c(mons$temp.mean,0,0), col=alpha("darkblue",.4),border=NA)
  axis(4, col="darkblue", col.axis="darkblue")
  mtext("mean temp C", side=4, line=2.2, col="darkblue", cex=.8)
  mtext(s,side=3, adj=.9, line=-1.5)
  par(new=T)
  plot(mons$month, mons$temp.amp, type="l", lwd=2, axes=F, xlab="", ylab = "delta temp C", xpd=T, ylim = c(1,4))
  axis(2)
  axis(1,labels=F )
  t <- seq(2,12, by=2)
  text(x = t, y = .8 , labels = as.character(month(t, label=T, abbr=T)), 
       pos = 1, xpd=NA,  cex=.9)
}

dev.off()

