library(RODBC)
biomass= odbcConnect("bio")

met =sqlQuery(biomass,"select distinct(stnum,reg_code,lat,lon) from meteorology.fire_weather_ibra")

met$row = as.character(met$row)
for(a in 1:length(met$row)){
	tstr=met$row[a]
	h=substr(tstr,2,nchar(tstr)-1)
	stnum = strsplit(h,",")[[1]][1]
	ibra = strsplit(h,",")[[1]][2]
	lat = strsplit(h,",")[[1]][3]
	lon = strsplit(h,",")[[1]][4]
	met$stnum[a] = stnum
	met$ibra[a] = ibra
	met$lat[a] = lat
	met$lon[a] = lon	
}

write.csv(met,"E:\\NASA_NERP\\HotspotFFDI\\data\\metstations.csv")

met <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\metstations.csv")

brlist = unique(met$ibra)
bigstore=data.frame(IBRA=NA,FS=NA,DOY=NA,FFDI=NA)

for(thebr in 1:length(brlist)){
print(thebr)
totest = subset(met,ibra == brlist[thebr])
nostats = length(totest$ibra)
bframe = data.frame(date = as.Date(seq(as.Date("2002-01-01"),as.Date("2007-12-31"),1)), IDMA = 0, hs = 0)


for(stsel in 1:length(totest$ibra)){
station = totest$stnum[stsel]
print(station)
lat = totest$lat[stsel]
lon = totest$lon[stsel]

data = sqlQuery(biomass,paste("select \"Date\",\"IDK\",\"IDMA\",\"modis\" from meteorology.fire_weather_ibra where stnum=",station,";",sep=""))
hots = sqlQuery(biomass,paste("select count(the_geom),ACQ_DATE from spatial.modis_hs where ST_Distance(the_geom, ST_GeomFromText('POINT(",lon," ",lat,")',4326)) < .5 group by ACQ_DATE;",sep=""))

data$modis = 0
data$Date = as.Date(data$Date)
data = subset(data,Date >= "2002-01-01" & Date <= as.Date("2007-12-31"))

for(a in 1:length(hots$acq_date)){
	data$modis[data$Date == hots$acq_date[a]] = hots$count[a]
}
bframe$IDMA = bframe$IDMA + data$IDMA
bframe$hs = bframe$hs + as.numeric(data$modis)
}

bframe$IDMA = bframe$IDMA / nostats
bframe$hs = bframe$hs / nostats

#library(quantreg)
#mod = rq(hs ~ IDMA,tau=2:98/100,data=bframe)

bframe$jday = strptime(bframe$date, "%Y-%m-%d")$yday+1
k=tapply(bframe$hs,bframe$jday,sum)
s_year = data.frame(jday = as.numeric(names(k)),hs=as.vector(k))
hsplit = c(s_year$hs,s_year$hs,s_year$hs)
jsplit = c(s_year$jday-365,s_year$jday,s_year$jday+365)
s_year2 = data.frame(jday = jsplit, hs = hsplit)

library(MASS)
library(gam)

 mod = gam(hs ~ s(jday,df=20),data=s_year2)
ss = mod$smooth[366:(366+365)]
fs_start = which.min(ss)
hist(bframe$hs)

myvec = ss
myvec=c(myvec[fs_start:length(myvec)],myvec[1:fs_start-1])
cd_vec=c()
tot=0
for(a in 1:length(myvec)){
	tot=tot+myvec[a]
	cd_vec=c(cd_vec,tot)
	}


d=1:366
d_vec=diff(myvec)

peakday = which.max(d_vec)
doy = peakday + fs_start
doy = doy %% 366

l=tapply(bframe$IDMA,bframe$jday,mean)
peakffdi = l[doy]

plvec = c(as.character(brlist[thebr]),fs_start,doy,peakffdi)
bigstore=rbind(bigstore,plvec)
}


write.csv(bigstore,"E:\\NASA_NERP\\HotspotFFDI\\data\\fire_start_v2.csv")
