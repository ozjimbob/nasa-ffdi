library(raster)
library(sp)
library(fields) 
library(compiler)
raster=cmpfun(raster)
Tps=cmpfun(Tps)
interpolate=cmpfun(interpolate)
subset=cmpfun(subset)
Krig=cmpfun(Krig)

root="G:\\wind\\synoptic_rasters\\syn3_"
data = read.csv("G:\\wind\\synoptic_3pm.csv")

template = raster("G:\\wind\\template.txt")

alldates = as.Date(sort(as.numeric(as.Date(as.character(unique(data$full_date))))),origin="1970-01-01")

alldates = alldates[alldates>as.Date("1900-01-01")]

data$full_date = as.numeric(as.Date(as.character(data$full_date)))

for(date_number in 1:length(alldates)){
date_data = subset(data,full_date == as.numeric(alldates[date_number]))
this_date=alldates[date_number]

xy=data.frame(x=date_data$lon,y=date_data$lat)
tps <- Tps(xy, date_data$ws_3pm)
p <- interpolate(template, tps,file=paste(root,this_date,".tif",sep=""),format="GTiff",datatype="INT2S")
print(this_date)
#image(p,main=this_date,col=rainbow(20))
}


data <- read.csv("R:\\SET\\PlantSci\\Staff\\Grant Williamson\\Publications\\Tng Giant Eucs\\giants\\giants.csv")
eucs=subset(data,Fam=="Eucalyptus")
ang=subset(data,Fam=="Angiosperm")
coni = subset(data,Fam=="Conifer")
plot(eucs$pet_he_yr ~ eucs$prec_ann,pch=15,xlim=c(250,4000),ylim=c(500,2000),xlab="Annual Precipitation (mm)",ylab="Annual Potential Evapotranspiration(mm)")
points(ang$pet_he_yr ~ ang$prec_ann,pch="o")
points(coni$pet_he_yr ~ coni$prec_ann,pch=17)
text(I(data$pet_he_yr+50) ~ I(data$prec_ann+50),labels=data$MapNm)

