
base::source("E:/NASA_NERP/HotspotFFDI/NASA_NERP_Hotspot FFDI/NASA_NERP_Hotspot FFDI/fire_index_functions.R")
library(compiler)
index_KBDI_new = cmpfun(index_KBDI_new)
index_MA = cmpfun(index_MA)
index_CFWI = cmpfun(index_CFWI)
index_Fix = cmpfun(index_Fix)
index_FixD = cmpfun(index_FixD)
index_FOS = cmpfun(index_FOS)

station_root="G:\\SILO_Data\\complete_stations\\"
stations=list.files(path = station_root)
station.list = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\metstations.csv")

for(astation in 1:length(stations)){
	station = stations[astation]
	this_station = substr(station,start=1,stop=(nchar(station)-4))
	print(paste("Station ",this_station," (",astation,"/",length(stations),")",sep=""))
	this_lat = subset(station.list,stnum==this_station)$lat[1]
	station_data = read.csv(paste(station_root,station,sep=""))
	station_data$c_date= as.Date(as.character(station_data$c_date))
 	allrf = sum(station_data$rain)
	mean_rf = allrf / ((length(station_data$rain))/365.25)
	out_frame = DoFDIs(station_data$t_max,station_data$rain,station_data$rh_max,station_data$wind_mean,station_data$wind_3pm,station_data$c_date,this_lat,mean_rf)
	final_frame=cbind(station_data,out_frame)
	write.csv(final_frame,paste(station_root,station,sep=""))
}

#Demo functions
DoFDIs <- function(Temp,Rain,Dew,WS,WS3,DDate,latitude,meanrf){
	
	Dew = calc_dew(Temp,Dew)
	Days <- as.numeric(format(strptime(DDate,format="%Y-%m-%d"),"%j"))
	Months <- as.numeric(format(strptime(DDate,format="%Y-%m-%d"),"%m"))
	IDK <- index_KBDI_new(Temp,Rain,meanrf)
	
	IDMA <- index_MA(Temp,Rain,Dew,meanrf,WS, IDK)
	IDMA_3 <- index_MA(Temp,Rain,Dew,meanrf,WS3, IDK)
	
	ICFWI <- index_CFWI(Months,Days,Temp,Dew,WS,Rain,DayLengths(latitude))
	ICFWI_3 <- index_CFWI(Months,Days,Temp,Dew,WS3,Rain,DayLengths(latitude))
	
	ISharp = index_Fix(Temp,WS,Dew)
	ISharp_3 = index_Fix(Temp,WS3,Dew)
	
	ISharpD = index_FixD(Temp,WS,Dew,IDK)
	ISharpD_3 = index_FixD(Temp,WS3,Dew,IDK)
	
	IFos <- index_FOS(Temp,WS,Dew)
	IFos_3 <- index_FOS(Temp,WS3,Dew)
	
	data.frame(cbind(IDK,IDMA,IDMA_3,ICFWI,ICFWI_3,ISharp,ISharp_3,ISharpD,ISharpD_3,IFos,IFos_3 ))
}



