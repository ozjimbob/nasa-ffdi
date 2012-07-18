library(raster)
library(compiler)
library(sp)
extract=cmpfun(extract)
raster=cmpfun(raster)

test = stack("G:\\wind\\daily_rasters\\daily.vrt")

base::source("E:/NASA_NERP/HotspotFFDI/NASA_NERP_Hotspot FFDI/NASA_NERP_Hotspot FFDI/silo_functions.R")


wind_mean_root="G:\\wind\\daily_rasters\\"
station_root="G:\\SILO_Data\\Patched_Point_Data\\"

stations=as.numeric(list.files(path = station_root))
daily.list = list.files(path = wind_mean_root)
daily.fp.list = paste(wind_mean_root,daily.list,sep="")
sub.list = daily.fp.list[1:20]
st=stack(daily.fp.list)
writeRaster(st, 'G:\\wind\\daily_rasters\\daily.nc', overwrite=TRUE) 


rast=list()
for(day in daily.list){
	print(day)
	rast[day] = raster(paste(wind_mean_root,day,sep=""))
	}

station.list = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\metstations.csv")


for(station in stations){
	this_station = stations[station]
	this_lat = subset(station.list,stnum==this_station)$lat[1]
	this_lon = subset(station.list,stnum==this_station)$lon[1]
	cds <- SpatialPoints(cbind(this_lon,this_lat),proj4string=CRS("+proj=longlat +datum=WGS84"))
	station_data = SILO_Get_Station(this_station)
	station_data = subset(station_data,c_date >= as.Date("1900-01-02") & c_date <= as.Date("2011-04-30"))
	station_data$wind_daily=NA
	station_data$wind_3pm=NA
	for(line in 1:length(station_data$c_date)){
		this_date = as.character(station_data$c_date[line])
		print(this_date)
		daily.file = paste("daily_",this_date,".tif",sep="")
		if(daily.file %in% daily.list){
			daily = raster(paste(wind_mean_root,daily.file,sep=""))
			val = extract(daily,cds)
			station_data$wind_daily[line]=val
		}
		syn.file = paste("syn3_",this_date,".tif",sep="")
		if(syn.file %in% syn.list){
			syn = raster(paste(wind_3pm_root,syn.file,sep=""))
			val = extract(syn,cds)
			station_data$wind_3pm[line]=val
		}
	}
}