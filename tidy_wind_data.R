library("RODBC")
postgis = odbcConnect("postgis")

rootdir = "G:\\wind\\dwind\\"
statlist = read.csv("G:\\wind\\mws_stations.csv")
colnames(statlist)[2]="stnum"
colnames(statlist)[7]="lat"
colnames(statlist)[8]="lon"

files = list.files(path = rootdir,full.names = TRUE)
full_frame = data.frame(stnum=NA,full_date=NA,ws_mean=NA)

for(count_file in 1:length(files)){
this_file = files[count_file]
print(this_file)
file_data = read.csv(this_file)
file_data$full_date = as.Date(paste(file_data$Year,"-",file_data$Month,"-",file_data$Day,sep=""))
mini_frame=data.frame(stnum=as.numeric(file_data$Station.Number),full_date = as.Date(file_data$full_date),ws_mean=as.numeric(file_data$Mean.daily.wind.speed.in.km.h))
mini_frame = subset(mini_frame,!is.na(ws_mean))
full_frame=rbind(full_frame,mini_frame)
}

full_frame$lat=NA
full_frame$lon= NA

lat_vec = full_frame$lat
lon_vec = full_frame$lon
st_vec = full_frame$stnum

for(station in statlist$stnum){
	print(station)
	the_lat = subset(statlist,stnum==station)$lat[1]
	the_lon = subset(statlist,stnum==station)$lon[1]
	lat_vec[st_vec == station] = the_lat
	lon_vec[st_vec == station] = the_lon
}

full_frame$lat = lat_vec
full_frame$lon = lon_vec

full_frame$full_date = as.character(as.Date(full_frame$full_date,origin="1970-01-01"))
full_frame=subset(full_frame,!is.na(full_frame$full_date))

write.csv(full_frame,"G:\\wind\\daily.csv")

#####

rootdir = "G:\\wind\\synwind\\"
statlist = read.csv("G:\\wind\\ws_stations.csv")
colnames(statlist)[2]="stnum"
colnames(statlist)[7]="lat"
colnames(statlist)[8]="lon"

files = list.files(path = rootdir,full.names = TRUE)
full_frame = data.frame(stnum=NA,full_date=NA,ws_3pm=NA)


for(count_file in 1:length(files)){
this_file = files[count_file]
print(this_file)
file_data = read.csv(this_file)

file_sdata=subset(file_data,Hour>= 14 & Hour <=16)
if(length(file_sdata$Month)==0){
	file_sdata=subset(file_data,Hour>= 8 & Hour <=10)
}
if(length(file_sdata$Month)==0){
	file_sdata=subset(file_data,Hour>= 11 & Hour <=13)
}
if(length(file_sdata$Month)==0){
	file_sdata=subset(file_data,Hour==3)
}
file_data = file_sdata
file_data$full_date = as.Date(paste(file_data$Year,"-",file_data$Month,"-",file_data$Day,sep=""))
mini_frame=data.frame(stnum=as.numeric(file_data$Station.Number),full_date = as.Date(file_data$full_date),ws_3pm=as.numeric(file_data$Wind.speed.measured.in.km.h))
mini_frame = subset(mini_frame,!is.na(ws_3pm))
#full_frame=rbind(full_frame,mini_frame)

if(count_file==1){
	write.table(mini_frame,"G:\\wind\\synoptic_3pm.csv",col.names = T,row.names=F,sep=",")
}
else
{
	write.table(mini_frame,"G:\\wind\\synoptic_3pm.csv",append=T,col.names = F,row.names=F,sep=",")
}
}

full_frame = read.csv("G:\\wind\\synoptic_3pm.csv")
full_frame$full_date = as.Date(as.character(full_frame$full_date))



full_frame$lat=NA
full_frame$lon= NA

lat_vec = full_frame$lat
lon_vec = full_frame$lon
st_vec = full_frame$stnum

for(station in statlist$stnum){
	print(station)
	the_lat = subset(statlist,stnum==station)$lat[1]
	the_lon = subset(statlist,stnum==station)$lon[1]
	lat_vec[st_vec == station] = the_lat
	lon_vec[st_vec == station] = the_lon
}

full_frame$lat = lat_vec
full_frame$lon = lon_vec

full_frame$full_date = as.character(as.Date(full_frame$full_date,origin="1970-01-01"))
full_frame=subset(full_frame,!is.na(full_frame$full_date))

write.csv(full_frame,"G:\\wind\\synoptic_3pm.csv")
