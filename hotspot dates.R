data <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\hotspots\\hotspots.csv",stringsAsFactors=F)

datelist=c()
f=unlist(e)
v=rep(c(TRUE,FALSE),length(f)/2)
g=f[v==TRUE]
this=strsplit(g,"/")

td=function(v){
	paste(v[3],"-",v[2],"-",v[1],sep="")
}
this2=lapply(this,td)

this2=unlist(this2)
data$acq_date = this2

library(compiler)

datelist=data$acq_date
timelist=data$local_time

for(record in 1:length(datelist)){
	print(record)
	date=as.Date(datelist[record])
	time=timelist[record]
	if(time >= 2400){
		time=time-2400
		date=date+1
		datelist[record]=as.character(date)
		timelist[record]=time
	}
}

data$acq_date = datelist
data$local_time = timelist

data$pass=""
data$pass[data$pass == ""] = "D"



write.csv(data,"E:\\NASA_NERP\\HotspotFFDI\\data\\hotspots\\hotspots2.csv",row.names = F)

### More hotspot fixes///

data <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\hotspots\\hotspots.csv",stringsAsFactors=F)
