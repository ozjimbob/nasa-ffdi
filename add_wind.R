library(raster)
library(compiler)
library(sp)
library(zoo)
extract=cmpfun(extract)
raster=cmpfun(raster)

base::source("E:/NASA_NERP/HotspotFFDI/NASA_NERP_Hotspot FFDI/NASA_NERP_Hotspot FFDI/silo_functions.R")


wind_mean_root="G:\\wind\\daily_rasters\\"
wind_3pm_root="G:\\wind\\synoptic_rasters\\"
station_root="G:\\SILO_Data\\Patched_Point_Data\\"

stations=as.numeric(list.files(path = station_root))

mean.daily.list = list.files(path = wind_mean_root,pattern="daily")
syn.daily.list = list.files(path = wind_3pm_root,pattern="syn3")
mean.date.list=substr(mean.daily.list,7,16)
syn.date.list=substr(syn.daily.list,6,15)




station.list = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\metstations.csv")

for(station in 1:length(stations)){
	this_station = stations[station]
	print(paste("Station ",this_station," (",station,"/",length(stations),")",sep=""))
	this_lat = subset(station.list,stnum==this_station)$lat[1]
	this_lon = subset(station.list,stnum==this_station)$lon[1]
	cds <- SpatialPoints(cbind(this_lon,this_lat),proj4string=CRS("+proj=longlat +datum=WGS84"))
	station_data = SILO_Get_Station(this_station)
	station_data = subset(station_data,c_date >= as.Date("1900-01-02") & c_date <= as.Date("2011-04-30"))
	
	setwd(wind_mean_root)
	tot=c()
	for(year in 1900:2011){
	d= system(paste("C:\\OSGeo4W\\bin\\gdallocationinfo -geoloc stack",year,".tif ",this_lon," ",this_lat,sep=""),intern=T)
	d=d[grep("Value:",d,fixed=T)]
	d=unlist(strsplit(d,":"))
	nolist=grep("Value",d,fixed=T)
	nolist=nolist+1
	d=as.numeric(d[nolist])
	tot=c(tot,d)
	}
    to.merge=data.frame(c_date=as.Date(mean.date.list),wind_mean=tot)
	station_with_mean <- merge(station_data,to.merge,by="c_date",all.x=T)
	
	setwd(wind_3pm_root)
	tot=c()
	for(year in 1900:2011){
	d= system(paste("C:\\OSGeo4W\\bin\\gdallocationinfo -geoloc stack",year,".tif ",this_lon," ",this_lat,sep=""),intern=T)
	d=d[grep("Value:",d,fixed=T)]
	d=unlist(strsplit(d,":"))
	nolist=grep("Value",d,fixed=T)
	nolist=nolist+1
	d=as.numeric(d[nolist])
	tot=c(tot,d)
	}
    to.merge=data.frame(c_date=as.Date(syn.date.list),wind_3pm=tot)
	station.filled <- merge(station_with_mean,to.merge,by="c_date",all.x=T)
	
	
	write.csv(station.filled,paste("G:\\SILO_Data\\complete_stations\\",this_station,".csv",sep=""))
	
}