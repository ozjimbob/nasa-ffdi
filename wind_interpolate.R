
library(RODBC)
library(mgcv)
library(fields)
source("R:\\SET\\PlantSci\\Staff\\Grant Williamson\\code\\space time met modelling\\elevation_functions.r")

biomass= odbcConnect("bio")
weather= odbcConnect("weather")

#pull in daily temp data for Jan 06 (only inside oz)

#Pull in list of dates.

dlist = seq(as.Date("1990-01-01"),as.Date("2010-12-31"),1)


stdata=sqlQuery(weather,"SELECT ele, stnum, name, lat, lon FROM weather_bom.combstats")
for(loc in 1:length(stdata$lon)){
if(is.na(stdata$ele[loc])) stdata$ele[loc] = ElePoint(stdata$lat[loc],stdata$lon[loc])
}
stdata$ele <- as.numeric(stdata$ele)
stdata$ele[is.na(stdata$ele)] <- 0
###

for(the_date in 1:length(dlist)){
rdate <- as.character(dlist[the_date])
year <- substr(rdate,1,4)
cat(rdate,"\n")

#Retrieve data for the day for which we do have temperature values, ie. model
thedata=sqlQuery(weather,paste("select * from weather_bom.combstats left join (select to_date(year || '-' || month || '-' || day,'YYYY-MM-DD') as date, mean_daily_wind_speed_in_km_h as temp, station_number as stnum from weather_bom.bom_daily_data_",year," where to_date(year || '-' || month || '-' || day,'YYYY-MM-DD') =  '",rdate,"') as goo on (weather_bom.combstats.stnum = goo.stnum);",sep=""))



cat("Got data\n")
#If met stations are missing elevation for the model set, fill in elevation from GEM

thedata$ele=0
for(loc in 1:length(thedata$stnum)){
thedata$ele[loc] = stdata$ele[stdata$stnum == thedata$stnum[loc]]
}

hasdata = subset(thedata,is.na(temp) == F)
nodata = subset(thedata,is.na(temp)== T)

hasdata$Modelled = "N"
nodata$Modelled = "Y"


cat("Got elevation\n")

gamMod=gam(temp~s(lon,lat,k=20)+s(ele,k=3),data=hasdata,family=gaussian(link="identity"))
pred=predict(gamMod,newdata=nodata,type='response')
nodata$mtemp <- pred
hasdata$mtemp <- hasdata$temp
cat("Run models\n")

mdata<-cbind(hasdata$stnum,rdate,hasdata$mtemp,hasdata$Modelled)
pdata<-cbind(nodata$stnum,rdate,nodata$mtemp,nodata$Modelled)
out <- data.frame(rbind(mdata,pdata))
names(out) = c("stnum","date","wind","Modelled")
write.table(out,"E:\\NASA_NERP\\HotspotFFDI\\datawind.csv",append=TRUE,sep=",",col.names = FALSE)
}




str(jan06)

# what family for GAM?
hist(jan06$temp)
#Gaussian eh?
# loop through each day and run a spatial gam, then predict the values on a grid and display
dates=names(table(jan06$dates))
for(i in 1:length(dates)){
test=jan06[jan06$dates==dates[i],]
attach(test)
gamMod=gam(temp~s(lon,lat),data=test)
res=100
xs=seq(min(test$lon),max(test$lon),len=res)
ys=seq(min(test$lat),max(test$lat),len=res)
gr=expand.grid(xs,ys)
pred=predict(gamMod,newdata=list(lon=gr[,1],lat=gr[,2]),se.fit=T,type='response')
image.plot(xs,ys,matrix(pred$fit,res,res),col=tim.colors(32),main=dates[i]) # or pred$se.fit
points(lon,lat)
detach(test)
}
