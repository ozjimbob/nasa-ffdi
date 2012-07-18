library(MASS)
library(gam)
library(compiler)
library(shape)
library(sp)
library(maptools)
library(signal)
library(fields)
library(KernSmooth)
library(glogis)
library(car)
###### http://cran.r-project.org/web/packages/biseVec/biseVec.pdf
library(biseVec)

search_year=cmpfun(function(hs_vec,start){
	hs.max = max(hs_vec)
	hs.min = min(hs_vec)
	pc10 = hs.min + (hs.max-hs.min)*0.1
	pc80 = hs.min + (hs.max-hs.min)*0.8
	peak_day=which.max(hs_vec)
	peak_hs=hs_vec[peak_day]
	
	## Find broad boundaries
	for(test in peak_day:1){
		if(hs_vec[test] <= pc10)
		{
			fire_season_increase = test
			fire_season_increase_hs = hs_vec[test]
			break
		}
	}
	for(test in peak_day:366){
		if(hs_vec[test] <= pc10)
		{
			fire_season_decrease = test
			break
		}
	}	
	for(test in peak_day:366){
		if(hs_vec[test] <= pc80)
		{
			po.peak = test
			break
		}
	}	
	for(test in peak_day:1){
		if(hs_vec[test] <= pc80)
		{
			pr.peak = test
			break
		}
	}
		
	#fire_season_increase=min(which(hs_vec>pc10))
	#fire_season_increase_hs=hs_vec[min(which(hs_vec>pc10))]
	pre_season_mean = mean(hs_vec[1:fire_season_increase-1])
	#fire_season_decrease=max(which(hs_vec>pc10))
	#pr.peak=min(which(hs_vec>pc80))
	#po.peak=max(which(hs_vec>pc80))

	
	c(fire_season_increase,fire_season_increase_hs,fire_season_decrease,peak_day,peak_hs,pre_season_mean,pr.peak,po.peak)
	
})

plot_year=cmpfun(function(data,cluster_subset,the_year,bandw=20){


data$acq_date = as.Date(data$acq_date)
data$year=as.numeric(format(data$acq_date,format="%Y"))
data$jd=as.numeric(format(data$acq_date,format="%j"))


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

print("Analyzing years...")

	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	print(start_date2)
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){return}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	out_vec=search_year(g)

	test_vec=out_vec
	
	A.start_date = as.character(start_date2)
	A.inc=test_vec[[1]]
	A.inc_hs=test_vec[[2]]
	A.dec=test_vec[[3]]
	A.peak=test_vec[[4]]
	A.peak_hs=test_vec[[5]]
	A.psm=test_vec[[6]]
	A.pr.peak=test_vec[[7]]
	A.po.peak=test_vec[[8]]
	A.points=hs_vec
	A.smooth=g


	data = subset(data,pass=="N",stringsAsFactors=F)
	
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

#print("Running spline.")
#s_year2 = data.frame(jday = jsplit, hs = hsplit)
#mod = gam(hs ~ s(jday,df=20),data=s_year2)
#ss = mod$smooth
#ss=ss[365:(365*2)]
#fs_start = which.min(ss)

####
#start_date = min(data$acq_date)
#end_date = max(data$acq_date)

#date_seq = seq(as.Date(start_date),as.Date(end_date),1)
#full_hs = data.frame(date_seq,hs=0)

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

print("Analyzing years...")

	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	print(start_date2)
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){return}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	out_vec=search_year(g)

	test_vec=out_vec
	
	N.start_date = as.character(start_date2)
	N.inc=test_vec[[1]]
	N.inc_hs=test_vec[[2]]
	N.dec=test_vec[[3]]
	N.peak=test_vec[[4]]
	N.peak_hs=test_vec[[5]]
	N.psm=test_vec[[6]]
	N.pr.peak=test_vec[[7]]
	N.po.peak=test_vec[[8]]
	N.points=hs_vec
	N.smooth=g
	
	### Make plot
clusters = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.csv")
cluster_name = as.character(subset(clusters,GRIDCODE==cluster_subset)$Name[1])
A.smooth=c(0,A.smooth,0)
N.smooth=c(0,N.smooth,0)
Xs=as.Date(as.Date(A.start_date)+seq_along(A.smooth))
print(A.start_date)
print(Xs)
Ys=format(Xs,format="%Y")
print(plot(rep(0,length(Xs))~Xs,xlim=c(0,366),ylim=c(0,max(A.points)),type="l",main=paste(cluster_name,".",Ys[1],sep=""),xaxt="n",xlab="Date",ylab="Hotspots"))
print(axis(side=1,at=seq_along(Xs),labels=months(Xs)))
print(polygon(A.smooth,col="grey80"))
print(polygon(N.smooth,col="grey20"))
print(points(A.points,cex=0.4,pch="x",col="red"))
print(points(N.points,cex=0.4,pch="x",col="blue"))

print(segments(A.inc+1,A.smooth[A.inc+1], A.dec+1,A.smooth[A.dec+1]))
print(segments(A.pr.peak+1,A.smooth[A.pr.peak+1], A.po.peak+1,A.smooth[A.po.peak+1]))
print(segments(A.peak+1,A.peak_hs,A.peak+1,0))

print(segments(N.inc+1,N.smooth[N.inc+1], N.dec+1,N.smooth[N.dec+1],col="white"))
print(segments(N.pr.peak+1,N.smooth[N.pr.peak+1], N.po.peak+1,N.smooth[N.po.peak+1],col="white"))
print(segments(N.peak+1,N.peak_hs,N.peak+1,0,col="white"))
})


data <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\hotspots.csv",stringsAsFactors=F)
plot_year(data,40,1,10)

plot_year(data,1,1,10)
