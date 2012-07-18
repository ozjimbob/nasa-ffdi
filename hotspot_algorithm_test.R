library(MASS)
library(gam)
library(compiler)
library(shape)
library(sp)
library(maptools)

wedge_geo=function(start,end,mid,s_sd,e_sd,m_sd,co){
	radi=1
	if(is.na(s_sd)){s_sd=1}
	if(is.na(e_sd)){e_sd=1}
	if(is.na(m_sd)){m_sd=1}
	start_h=(-(start-s_sd/2))+90
	start_l=(-(start+s_sd/2))+90
	start=(-(start))+90
	end_h=(-(end-e_sd/2))+90
	end_l=(-(end+e_sd/2))+90
	end=(-(end))+90
	mid1=(-(mid-m_sd/2))+90
	mid2=(-(mid+m_sd/2))+90
	mid=(-(mid))+90
	if(abs((start %% 360) - (end %% 360)) < 3){return()}
	#start=start %% 2*pi
	#end=end %% 2*pi
	#mid1 = mid1 %% 2*pi
	#mid2 = mid2 %% 2*pi
	#end_h = end_h %% 2*pi
	#end_l = end_l %% 2*pi
	#start_h = start_h %% 2*pi
	#end_h = end_h %% 2*pi
	filledellipse(rx1=0,rx2=radi,ry1=0,ry2=radi,mid=c(co[1],co[2]),dr=0.01,angle=0,0,2*pi,col="white",lcol="black",lwd=0.5)
	filledellipse(rx1=0,rx2=radi,ry1=0,ry2=radi,mid=c(co[1],co[2]),dr=0.01,angle=0,from=rad(start_l),to=rad(start_h),col="pink",lwd=0.5)
	filledellipse(rx1=0,rx2=radi,ry1=0,ry2=radi,mid=c(co[1],co[2]),dr=0.01,angle=0,from=rad(end_l),to=rad(end_h),col="pink",lwd=0.5)
	filledellipse(rx1=0,rx2=radi,ry1=0,ry2=radi,mid=c(co[1],co[2]),dr=0.01,angle=0,from=rad(end),to=rad(start),col=2,lwd=0.5)
	filledellipse(rx1=0,rx2=radi,ry1=0,ry2=radi,mid=c(co[1],co[2]),dr=0.01,angle=0,from=rad(mid2),to=rad(mid1),col=1,lwd=0.5)
	#k=getellipse(1.1,1.1,c(0,0),(2*pi)/12,-10,0,2*pi)
	#k=k[1:12,]
	#labs=c("Apr","Mar","Feb","Jan","Dec","Nov","Oct","Sep","Aug","Jul","Jun","May")
	#text(k,labs)
	#title(the_title)
	}
wedge=function(start,end,mid,s_sd,e_sd,m_sd,the_title){
	emptyplot(xlim = c(-1.2, 1.2), ylim = c(-1.2,1.2))
	start_h=(-(start-s_sd/2))+90
	start_l=(-(start+s_sd/2))+90
	start=(-(start))+90
	end_h=(-(end-e_sd/2))+90
	end_l=(-(end+e_sd/2))+90
	end=(-(end))+90
	mid1=(-(mid-m_sd/2))+90
	mid2=(-(mid+m_sd/2))+90
	mid=(-(mid))+90
	if(start+3 > end | start - 3 > end){return()}
	filledellipse(1,1,c(0,0),from=0,to=2*pi)
	filledellipse(rx1=0,rx2=1,ry1=0,ry2=1,mid=c(0,0),dr=0.01,angle=0,from=rad(start_l),to=rad(start_h),col="pink")
	filledellipse(rx1=0,rx2=1,ry1=0,ry2=1,mid=c(0,0),dr=0.01,angle=0,from=rad(end_l),to=rad(end_h),col="pink")
	filledellipse(rx1=0,rx2=1,ry1=0,ry2=1,mid=c(0,0),dr=0.01,angle=0,from=rad(end),to=rad(start),col=2)
	filledellipse(rx1=0,rx2=1,ry1=0,ry2=1,mid=c(0,0),dr=0.01,angle=0,from=rad(mid2),to=rad(mid1),col=1)
	k=getellipse(1.1,1.1,c(0,0),(2*pi)/12,-10,0,2*pi)
	k=k[1:12,]
	labs=c("Apr","Mar","Feb","Jan","Dec","Nov","Oct","Sep","Aug","Jul","Jun","May")
	text(k,labs)
	title(the_title)
}

search_year=cmpfun(function(hs_vec,start){
	fire_season_increase=0
	fire_season_increase_hs=0
	for(wind in start:(length(hs_vec)-3)){
		sub_vec=hs_vec[wind:(wind+2)]
		if(sub_vec[3] > sub_vec[2] & sub_vec[2] > sub_vec[1] & sub_vec[1] > 0){
			fire_season_increase=wind
			fire_season_increase_hs=sub_vec[3]
			break
		}
	}
	if(fire_season_increase == 0){
		return(c(0,0,0,0,0))
	}
	#pre_season_mean = mean(hs_vec[1:(fire_season_increase-1)])
	pre_season_mean = mean(c(hs_vec[1:14],hs_vec[350:365]))
	
	## Search for end of season
	fire_season_decrease=0
	if(fire_season_increase >= (365-18)){
			return(c(0,0,0,0,0))
		}
	for(wind in (fire_season_increase+3):(length(hs_vec)-18)){
			sub_vec=hs_vec[wind:(wind+14)]
			#test_sub = sum(as.numeric(sub_vec>=fire_season_increase_hs))
			test_sub = mean(sub_vec)
			if(test_sub<=pre_season_mean){
				fire_season_decrease=wind
				break
			}
		}
	if(fire_season_decrease == 0){
	for(wind in (length(hs_vec)-3):1){
		sub_vec=hs_vec[wind:(wind+2)]
		if(sub_vec[3] > 0 & sub_vec[3] < sub_vec[2] & sub_vec[2] < sub_vec[1]){
			fire_season_decrease=wind
			break
		}
	}
	}
	if(fire_season_decrease == 0){
		return(c(0,0,0,0,0))
		}
	## Search for peak day within range
	fire_season_peak = 0
	sub_vec=hs_vec[fire_season_increase:fire_season_decrease]
	peak_hs=max(sub_vec)
	peak_day=fire_season_increase+which.max(sub_vec)+1

	c(fire_season_increase,fire_season_increase_hs,fire_season_decrease,peak_day,peak_hs)
	
})

process_cluster=cmpfun(function(data,cluster_subset){
print(paste("Subsetting cluster:",cluster_subset,sep=""))
cluster_hotspots = subset(data,data$Cluster==cluster_subset)
hotspot_collapse=tapply(cluster_hotspots$Join_Count,cluster_hotspots$jd,sum)
hotspot_collapse = data.frame(jday = as.numeric(names(hotspot_collapse)),hs=as.vector(hotspot_collapse))
print("Filling blank days.")
filled=data.frame(jday = seq(1,365),hs=0)
for(rec in 1:length(filled$jday)){
	grab = subset(hotspot_collapse,hotspot_collapse$jday==rec)
	if(length(grab$jday)!=0){
		filled$hs[rec]=grab$hs[1]
	}
}
print("Creating merged cycle.")
hsplit = c(filled$hs,filled$hs,filled$hs)
jsplit = c(filled$jday-365,filled$jday,filled$jday+365)

print("Running spline.")
s_year2 = data.frame(jday = jsplit, hs = hsplit)
mod = gam(hs ~ s(jday,df=20),data=s_year2)
ss = mod$smooth
ss=ss[365:(365*2)]
fs_start = which.min(ss)

####
start_date = min(data$acq_date)
end_date = max(data$acq_date)

date_seq = seq(as.Date(start_date),as.Date(end_date),1)
full_hs = data.frame(date_seq,hs=0)

print("Filling individual days.")
for(the_date in 1:length(date_seq)){
	subf=subset(cluster_hotspots,cluster_hotspots$acq_date==date_seq[the_date])
	num = length(subf$acq_date)
	full_hs$hs[the_date] = num
}

full_hs$jday=as.numeric(format(full_hs$date_seq,format="%j"))
full_hs_t=tapply(full_hs$hs,full_hs$jday,sum)
if(fs_start==366){fs_start=365}

print("Identifying year passes.")
start_points=as.numeric(rownames(subset(full_hs,full_hs$jday==fs_start)))
end_points=start_points+365
count_points=length(end_points)
if(end_points[count_points] > length(full_hs$jday)){
	start_points=start_points[1:(count_points-1)]
	end_points=end_points[1:(count_points-1)]
}


cluster_summary=data.frame(clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")
for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	
	pre_season_mean = mean(c(hs_vec[1:14],hs_vec[350:365]))
	
	out_vec=search_year(hs_vec,1)
	newstart=out_vec[3]
	out_vec2=search_year(hs_vec,newstart)
	newstart=out_vec2[3]+1
	out_vec3=search_year(hs_vec,newstart)
	newstart=out_vec3[3]+1
	out_vec4=search_year(hs_vec,newstart)
	newstart=out_vec4[3]+1
	out_vec5=search_year(hs_vec,newstart)
	
	test_frame=data.frame(rbind(out_vec,out_vec2,out_vec3,out_vec4))
	test_frame$diff = test_frame$X3-test_frame$X1
	row=which.max(test_frame$diff)
	test_vec=test_frame[row,]
	
	cluster_summary$st_date[the_year] = as.character(start_date)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=pre_season_mean
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
print(plot(full_hs$hs ~ full_hs$jday))
print(abline(v=(cluster_summary$inc + fs_start) %% 365,col=3))
print(abline(v=(cluster_summary$dec + fs_start) %% 365,col=5))
print(abline(h=(cluster_summary$peak_hs),col=6))
print(abline(v=(cluster_summary$peak + fs_start) %% 365,col=6))
cluster_summary

})
rad=cmpfun(function(deg){deg*pi/180})
deg=cmpfun(function(rad){rad*180/pi})
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


data <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\hotspots.csv",stringsAsFactors=F)
data = subset(data,pass=="N",stringsAsFactors=F)

data$acq_date = as.Date(data$acq_date)
data$year=as.numeric(format(data$acq_date,format="%Y"))
data$jd=as.numeric(format(data$acq_date,format="%j"))

output_frame=data.frame(clust=NA,st_date=NA,fs_start=NA,inc=NA,inc_hs=NA,dec=NA,peak=NA,peak_hs=NA,len=NA,psm=NA,stringsAsFactors=F)

for(the_cluster in 1:61){
	frame=process_cluster(data,the_cluster)
	output_frame=rbind(output_frame,frame)
}

##################################################################


output_frame = subset(output_frame,inc > 0)
output_frame$inc = output_frame$inc + output_frame$fs_start
output_frame$dec = output_frame$dec + output_frame$fs_start
output_frame$peak = output_frame$peak + output_frame$fs_start
output_frame$inc = output_frame$inc %% 365
output_frame$dec = output_frame$dec %% 365
output_frame$peak = output_frame$peak %% 365

output_frame$inc_c = cos(rad(output_frame$inc*(360/366)))
output_frame$dec_c = cos(rad(output_frame$dec *(360/366)))
output_frame$peak_c = cos(rad(output_frame$peak *(360/366)))

output_frame$inc_s = sin(rad(output_frame$inc*(360/366)))
output_frame$dec_s = sin(rad(output_frame$dec *(360/366)))
output_frame$peak_s = sin(rad(output_frame$peak *(360/366)))

inc_c = as.numeric(tapply(output_frame$inc_c,output_frame$clust,mean))
dec_c = as.numeric(tapply(output_frame$dec_c,output_frame$clust,mean))
peak_c = as.numeric(tapply(output_frame$peak_c,output_frame$clust,mean))

inc_c_sd = as.numeric(tapply(output_frame$inc_c,output_frame$clust,sd))
dec_c_sd = as.numeric(tapply(output_frame$dec_c,output_frame$clust,sd))
peak_c_sd = as.numeric(tapply(output_frame$peak_c,output_frame$clust,sd))

inc_s = as.numeric(tapply(output_frame$inc_s,output_frame$clust,mean))
dec_s = as.numeric(tapply(output_frame$dec_s,output_frame$clust,mean))
peak_s = as.numeric(tapply(output_frame$peak_s,output_frame$clust,mean))

inc_s_sd = as.numeric(tapply(output_frame$inc_s,output_frame$clust,sd))
dec_s_sd = as.numeric(tapply(output_frame$dec_s,output_frame$clust,sd))
peak_s_sd = as.numeric(tapply(output_frame$peak_s,output_frame$clust,sd))

inc=(deg(atan2(inc_s,inc_c)) %% 360) /(360/365)
dec=(deg(atan2(dec_s,dec_c)) %% 360) /(360/365)
peak=(deg(atan2(peak_s,peak_c)) %% 360) /(360/365)

inc_sd=(deg(atan2(inc_s_sd,inc_c_sd)) %% 360) /(360/365)
dec_sd=(deg(atan2(dec_s_sd,dec_c_sd)) %% 360) /(360/365)
peak_sd=(deg(atan2(peak_s_sd,peak_c_sd)) %% 360) /(360/365)

rows = as.numeric(rownames(tapply(output_frame$fs_start,output_frame$clust,mean)))
fs_start = as.numeric(tapply(output_frame$fs_start,output_frame$clust,mean))
inc_hs = as.numeric(tapply(output_frame$inc_hs,output_frame$clust,mean))
peak_hs = as.numeric(tapply(output_frame$peak_hs,output_frame$clust,mean))
peak_hs_sd = as.numeric(tapply(output_frame$peak_hs,output_frame$clust,sd))
len = as.numeric(tapply(output_frame$len,output_frame$clust,mean))
len_sd = as.numeric(tapply(output_frame$len,output_frame$clust,sd))
psm = as.numeric(tapply(output_frame$psm,output_frame$clust,mean))


mean_frame=data.frame(clust=rows,fs_start,inc,inc_sd,inc_hs,dec,dec_sd,peak,peak_sd,peak_hs,peak_hs_sd,len,len_sd,psm)

map=readShapePoly("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.shp")
map.rec=data.frame(map)

write.csv(mean_frame,"E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial1\\mean_frame_N.csv",row.names=F)

m1 <- merge(map.rec, mean_frame, by.x = "GRIDCODE", by.y = "clust",all.x=T)
plot(map,col=m1$peak_hs)

X11(width = 11, height = 8.5)
plot(map)

for(record in 1:length(map.rec$GRIDCODE)){
	the_rec=map.rec$GRIDCODE[record]
	co=map@polygons[[record]]@labpt
	ss=subset(mean_frame,mean_frame$clust==the_rec)
	if(length(ss$clust)==0){next}
	wedge_geo(start=ss$inc,end=ss$dec,mid=ss$peak,s_sd=ss$inc_sd,e_sd=ss$dec_sd,m_sd=ss$peak_sd,co=co)

}