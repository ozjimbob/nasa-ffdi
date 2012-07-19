library(MASS)
library(gam)
library(compiler)
library(shape)
library(sp)
library(maptools)
library(lme4)
library(languageR)

#### For a given cluster, grab a 15-day window around the given date, and summarize into a series of data frames
get_met_data=cmpfun(function(cluster,the_date){
	start_interval = as.Date(as.Date(the_date)-7)
	end_interval = as.Date(as.Date(the_date)+7)
	label_seq = seq(-7,7,1)
	stations_to_get = subset(station_list,station_list$GRIDCODE==cluster)$no_
	if(length(stations_to_get)>20){
		stations_to_get = sample(stations_to_get,20)
	}
	end_frame = data.frame(stnum = NA, offset=NA, t_max=NA,t_min=NA,rain=NA,rh_max=NA,rh_min=NA,wind_mean=NA,wind_3pm=NA,IDK=NA,IDMA=NA,IDMA_3=NA,ICFWI=NA,ICFWI_3=NA,ISharp=NA,ISharp_3=NA,ISharpD=NA,ISharpD_3=NA,IFos=NA,IFos_3=NA)		
	for(idx in seq_along(stations_to_get)){
		station_data = read.csv(paste(station_root,stations_to_get[idx],".csv",sep=""))
		station_data$c_date = as.Date(station_data$c_date)
		station_data = subset(station_data,c_date >= start_interval & c_date <= end_interval)
		sub_frame = data.frame(stnum = station_data$stnum, offset=label_seq, t_max=station_data$t_max,t_min=station_data$t_min,rain=station_data$rain,rh_max=station_data$rh_max,rh_min=station_data$rh_min,wind_mean=station_data$wind_mean,wind_3pm=station_data$wind_3pm,IDK=station_data$IDK,IDMA=station_data$IDMA,IDMA_3=station_data$IDMA_3,ICFWI=station_data$ICFWI,ICFWI_3=station_data$ICFWI_3,ISharp=station_data$ISharp,ISharp_3=station_data$ISharp_3,ISharpD=station_data$ISharpD,ISharpD_3=station_data$ISharpD_3,IFos=station_data$IFos,IFos_3=station_data$IFos_3)
		end_frame = rbind(end_frame,sub_frame)
	}
	end_frame = subset(end_frame,!is.na(stnum))
	agg_frame=aggregate(end_frame,by=list(end_frame$offset),mean)
	pre_sub = subset(agg_frame,offset <= 0)
	post_sub = subset(agg_frame,offset >= 0)
	pre_mean = aggregate(pre_sub, by=list(pre_sub$stnum),mean)
	post_mean = aggregate(post_sub, by=list(post_sub$stnum),mean)
	day=subset(agg_frame,offset==0)
	output=list(full=end_frame,aggregate=agg_frame,pre_mean=pre_mean,post_mean=post_mean,day=day)
	output
})
subsamp <- cmpfun(function(de,ss){
de = de[rownames(de) %in% sample(rownames(de),ss),]
de
})
pv=cmpfun(pvals.fnc)
# For a windowed met-data set, run a t-test pre and post- the cutoff day, return a data frame with statistics for every variable
pre_post_t_test = cmpfun(function(full_met,mod.type="lme"){
full_frame = full_met
full_frame$slice = NA
full_frame$slice[full_frame$offset<0] = "Pre"
full_frame$slice[full_frame$offset>0] = "Post"
full_frame$slice = as.factor(full_frame$slice)
var.list=names(full_frame)[4:21]
out.frame = data.frame(varb=NA,est=NA,p=NA)
full_frame$c_date=as.character(full_frame$c_date)
for(idx in seq_along(var.list)){
	thevar=var.list[idx]
	if(mod.type=="glm"){
		formula = paste(thevar," ~ slice + c_date",sep="")
		model = glm(formula,data=full_frame)
		m.summary=summary(model)
		var.est = m.summary$coefficients[2,1]
		var.p = m.summary$coefficients[2,4]
	}else{
		formula = paste(thevar," ~ slice + (1|c_date)",sep="")
		formula.null = paste(thevar," ~ 1+ (1|c_date)",sep="")
		model <- lmer(formula, data=full_frame)
		model.null <- lmer(formula.null, data=full_frame)
		var.est = as.numeric(fixef(model)[2])
		aa=pv(model,nsim=1000,addPlot=F)
		var.est = aa$fixed[2,1]
		var.p = aa$fixed[2,6]
	}
	out_vec=c(thevar,var.est,var.p)
	out.frame=rbind(out.frame,out_vec)
}
out.frame=subset(out.frame,!is.na(varb))
out.frame
})

process_all_clusters=function(sset,t_var){
data = read.csv(paste(root_dir,input_file,sep=""))
data = subset(data,set==sset)
data$st_date = as.Date(data$st_date)
cluster.list = unique(data$clust)

for(idx.c in seq_along(cluster.list)){
	this.cluster = cluster.list[idx.c]
	this.data = subset(data,clust==this.cluster)
	test.first=1
	for(idx in seq_along(this.data$set)){

		row_data = this.data[idx,]
		if(is.na(row_data$st_date)){
			test.first = test.first + 1
			next
			}
    
		start_date = as.Date(row_data$st_date + as.numeric(row_data[t_var]))
		print(paste(idx/length(data$set)," - ",this.cluster," - ",start_date,sep=""))
		#end_date = as.Date(row_data$st_date + row_data$dec)
		#peak_date = as.Date(row_data$st_date + row_data$peak)
		#pre_date = as.Date(row_data$st_date + row_data$pr.peak)
		#post_date = as.Date(row_data$st_date + row_data$po.peak)
		cluster = row_data$clust
		
		print("Starting met")
		start_met = get_met_data(cluster,start_date)
		print("Got met")
		comb.met.d = start_met$aggregate
		mean.met.d = start_met$day
		mean.met.d$c_date = start_date
		comb.met.d$c_date = start_date
		if(idx==test.first){
			comb.met=comb.met.d
			}else{
			comb.met=rbind(comb.met,comb.met.d)
		}
		if(idx==test.first){
			mean.met=mean.met.d
			}else{
			mean.met=rbind(mean.met,mean.met.d)
		}
	}

	print("Running Statistics")
	var.summary = pre_post_t_test(comb.met)
	row_data[paste(var.summary$varb,".est",sep="")] = var.summary$est
	row_data[paste(var.summary$varb,".p",sep="")] = var.summary$p

	mean.met$Group.1 <- NULL
	mean.met$stnum <- NULL
	mean.met$offset <- NULL
	mean.met$c_date <- NULL

	this.data = subset(this.data,!is.na(st_date))
	new.data = cbind(this.data,mean.met)
	if(idx.c==1){
		new.yearly.data = new.data
		}else{
		new.yearly.data=rbind(new.yearly.data,new.data)
	}
	all.years.mean = mean(mean.met)

	row_data[paste(names(all.years.mean),".mean",sep="")] = as.numeric(all.years.mean)
	if(idx.c==1){
		new.cluster.data = row_data
		}else{
		new.cluster.data=rbind(new.cluster.data,row_data)
	}
  list(clust = new.cluster.data, year=new.yearly.data)
}}


root_dir="E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\"
out_root_dir="E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial3\\"
input_file="complete_frame.csv"
station_root = "G:\\SILO_Data\\complete_stations_min\\"
station_list = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\silo_stations.csv")

year.d=process_all_clusters("A","dec")

clust.d=new.cluster.data
year.n=new.yearly.data
clust.n=new.cluster.data
year.a = new.yearly.data
clust.a = new.cluster.data

cluster_year_frame = rbind(year.a,year.d,year.n)
cluster_frame = rbind(clust.a,clust.d,clust.n)

## Remove unwanted cluster_frame columns
cluster_frame["X"] = NULL
cluster_frame["st_date"] = NULL
cluster_frame["fs_start"] = NULL
cluster_frame["inc"] = NULL
cluster_frame["inc.hs"] = NULL
cluster_frame["dec"] = NULL
cluster_frame["peak"] = NULL
cluster_frame["peak_hs"] = NULL
cluster_frame["len"] = NULL
cluster_frame["psm"] = NULL
cluster_frame["pr.peak"] = NULL
cluster_frame["po.peak"] = NULL

cluster_frame=reshape(cluster_frame,varying=list(names(cluster_frame)[4:21],names(cluster_frame)[22:39],names(cluster_frame)[40:57]),v.names=c("est","p","mean"),direction="long",times=names(cluster_frame)[22:39])
cluster_frame["id"] = NULL
cluster_frame["inc_hs"] = NULL
names(cluster_frame)=c("set","clust","var","est","p","mean")

write.csv(cluster_frame,paste(out_root_dir,"cluster_start_stats.csv",sep=""))

### Remove unwanted cluster_year_frame columns
cluster_year_frame["X"] = NULL
ge=reshape(cluster_year_frame,varying=list(names(cluster_year_frame)[5:31]),v.names=c("mean"),direction="long",times=names(cluster_year_frame)[5:31])
write.csv(cluster_year_frame,paste(out_root_dir,"cluster_year_start_stats.csv",sep=""))




