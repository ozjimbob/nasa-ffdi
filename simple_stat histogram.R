library(RODBC)
biomass= odbcConnect("bio")

station=9574
data = sqlQuery(biomass,paste("select \"IDMA\" from meteorology.fire_weather_ibra where stnum=",station,";",sep=""))
hist(data$IDMA,xlim=c(0,50),main="Margaret River",xlab="FFDI")

