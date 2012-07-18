
library(RODBC)
library(mgcv)
library(fields)
source("R:\\SET\\PlantSci\\Staff\\Grant Williamson\\code\\space time met modelling\\elevation_functions.r")

biomass= odbcConnect("bio")
weather= odbcConnect("weather")


stdata=sqlQuery(weather,"SELECT ele, stnum, name, lat, lon FROM weather_bom.combstats")
for(loc in 1:length(stdata$lon)){
if(is.na(stdata$ele[loc])) stdata$ele[loc] = ElePoint(stdata$lat[loc],stdata$lon[loc])
}
stdata$ele <- as.numeric(stdata$ele)
stdata$ele[is.na(stdata$ele)] <- 0

stdata$FFDI=0
print(length(stdata.bak$lon))

adata = sqlQuery(biomass,"SELECT distinct(stnum) from meteorology.fire_weather_ibra")
stdata.bak = stdata
stdata = subset(stdata,stdata$stnum %in% adata$stnum)


for(station in 1:length(stdata$FFDI)){
	biomass= odbcConnect("bio")
	print(stdata$stnum[station])
	print(station / length(stdata$FFDI))
	alldays = sqlQuery(biomass,paste("SELECT \"IDMA\" from meteorology.fire_weather_ibra where stnum = ",stdata$stnum[station],";",sep=""))
	if(length(alldays$IDMA)==0){next}
	stdata$FFDI[station]=as.numeric(quantile(alldays$IDMA,.95))
	print("Done")
}

##### Try it backwards

alldays = sqlQuery(biomass,paste("SELECT stnum,\"IDMA\" from meteorology.fire_weather_ibra;",sep=""))

stationlist = unique(alldays$stnum)
bigframe = data.frame(stationlist,lat=0,lon=0,FFDI=0)

stdata=sqlQuery(weather,"SELECT ele, stnum, name, lat, lon FROM weather_bom.combstats")



for(station in 1:length(bigframe$stationlist)){
	print(station / length(bigframe$stationlist))
	lat = stdata$lat[stdata$stnum==bigframe$stationlist[station]]
	lon = stdata$lon[stdata$stnum==bigframe$stationlist[station]]
	ffdi = as.numeric(quantile(alldays$IDMA[alldays$stnum == bigframe$stationlist[station]],.95))
	bigframe$lat[station]=lat
	bigframe$lon[station]=lon
	bigframe$FFDI[station]=ffdi
	print("Done")
}


write.csv(bigframe,"E:\\tng\\envs\\ffdi95.csv")