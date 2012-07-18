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

cdf=function(g){
		h=g
	for(idx in 1:(length(g))){
		h[idx]=sum(g[1:idx])		
	}
	h
}
doublesigmoid <- function(res) {
        sigmoid <- function(res, maxvalue = 0, turnaround = FALSE) {
            days <- length(res)
            model <- vector(mode = "numeric", length = days)
            if (turnaround) {
                res <- res[days:1]
            }
            F <- mean(na.omit(res[1:(days/4)]))
            if (maxvalue > 0) {
                G <- maxvalue - F
            }
            else {
                res.order <- order(na.omit(res), decreasing = TRUE)
                meanmax <- mean(res[res.order[1:3]])
                G <- meanmax - F
            }
            model <- .C("sigmoid", rdays = as.integer(days), 
                ndvi = as.numeric(res), rF = as.numeric(F), rG = as.numeric(G), 
                model = as.numeric(model), PACKAGE = "biseVec")$model
            if (turnaround) {
                model <- model[length(model):1]
            }
            if (length(model) != days) {
                model <- rep(NA, days)
            }
            return(model)
        }
        days <- length(res)
        delay <- 10
        count <- 1
        while ((maxpos <- order(res, decreasing = TRUE)[count]) > 
            300) {
            count <- count + 1
            if ((count > 100) || (is.na(maxpos))) {
                return(rep(NA, days))
            }
        }
        modelfront <- sigmoid(res[1:(maxpos + delay)])
        if (length(na.omit(modelfront)) < 5) {
            return(rep(NA, days))
        }
        maxvalue <- max(na.omit(modelfront))
        modelback <- sigmoid(res[(maxpos - 1 + delay):days], 
            maxvalue = maxvalue, turnaround = TRUE)
        if ((length(modelback) < 5) || is.na(modelback[3])) {
            modelback <- rep(maxvalue, days - length(modelfront) + 
                2)
        }
        model <- vector(mode = "numeric", length = days)
        endfront <- length(modelfront)
        endback <- length(modelback)
        if (modelfront[endfront] != modelback[3]) {
            scal <- modelback[3]
            for (i in 1:endback) {
                modelback[i] <- (modelback[i]/scal) * modelfront[endfront]
            }
        }
        model <- c(modelfront, modelback[3:endback])
        return(model)
    }	
gf=function (ndvi, slidperiod = 40, globalthreshold = 0, asym = FALSE, resvec = c(0.55, 0.6, 0.65, 0.65, 0.7)) {
    days <- length(ndvi)
    model <- res <- correctedndvi <- vector(mode = "numeric", 
        length = days)
    #ndvi <- ndvi/10000
    #ndvi <- ifelse(ndvi <= 0 | ndvi > 1, NA, ndvi)
    #if (length(which(is.na(ndvi) == FALSE)) < 10) {
    #    return(rep(NA, 7))
    #}
    #if (ndvi[order(ndvi, decreasing = TRUE)[10]] < 0.1) {
    #    return(rep(NA, 7))
    #}
    ndvi <- ifelse(is.na(ndvi), 0, ndvi)
    res <- .C("bise", rdays = as.integer(days), ndvi = as.numeric(ndvi), 
        rslidperiod = as.integer(slidperiod), cndvi = as.numeric(correctedndvi), 
        PACKAGE = "biseVec")$cndvi
    res <- ifelse(res <= 0, NA, res)
    if (length(which(is.na(res) == FALSE)) < 5) {
        return(rep(NA, 7))
    }
    res <- rejectRapidIncrease(res)
    f <- approxfun(x = 1:(3 * days), y = c(res, res, res))
    temp <- f((days + 1):(2 * days))
    maxpos <- order(res, decreasing = TRUE)
    maxpos <- maxpos[which(!is.na(res[maxpos]))]
    count <- 1
    mu <- maxpos[count]
    while ((mu > 210) || (mu < 120)) {
        count <- count + 1
        if (count > length(maxpos)) {
            mu <- maxpos[1]
            break
        }
        mu <- maxpos[count]
    }
    sig <- sqrt(var(na.omit(temp)))
    base <- 0
    base <- mean(temp[1:(days/6)])
    if (asym) {
        res <- ifelse(is.na(res), -1, res)
        model <- .C("asymgauss", rdays = as.integer(days), ndvi = as.numeric(res), 
            mustart = as.integer(mu), sigstart = as.numeric(sig), 
            rbase = as.numeric(base), model = as.numeric(model), 
            PACKAGE = "biseVec")$model
    }
    else {
        ndvi <- res[which(is.na(res) == FALSE)]
        time <- which(is.na(res) == FALSE)
        G <- function(mu, sig, scal, base, time, days) {
            erg <- ((scal/(sig * sqrt(2 * pi))) * exp(-0.5 * 
                ((((time/days) - (mu/days))/sig)^2))) + base
            return(erg)
        }
        model.nls <- try(nls(ndvi ~ G(mu = mu, sig = sig, scal = scal, 
            base = base, time = time, days = days), start = list(sig = sig, 
            scal = 1, mu = mu), control = list(maxiter = 200)), 
            silent = TRUE)
        if (inherits(model.nls, "try-error") == FALSE) {
            model <- predict(model.nls, list(time = 1:days))
        }
        else {
            res <- ifelse(is.na(res), -1, res)
            model <- .C("gauss", rdays = as.integer(days), ndvi = as.numeric(res), 
                mustart = as.integer(mu), sigstart = as.numeric(sig), 
                rbase = as.numeric(base), model = as.numeric(model), 
                PACKAGE = "biseVec")$model
        }
    }
    if (length(which(is.na(model) == FALSE)) < 10) {
        return(rep(NA, 7))
    }
    sosvec <- vector(mode = "integer", length = 7)
    sosvec[1] <- localThresGreenup(model, resvec[1])
    sosvec[2] <- localThresGreenup(model, resvec[2])
    sosvec[3] <- localThresGreenup(model, resvec[3])
    sosvec[4] <- globalThresGreenup(model, resvec[4])
    sosvec[5] <- globalThresGreenup(model, resvec[5])
    sosvec[6] <- min(na.omit(model))
    sosvec[7] <- max(na.omit(model))
	return(model)
    #return(sosvec)
}
sf=function (ndvi, slidperiod = 40, mode = 0, resvec = c(0.55, 0.6,  0.65, 0.65, 0.7)) {
    days <- length(ndvi)
    res <- vector(mode = "numeric", length = days)
    correctedndvi <- vector(mode = "numeric", length = days)
    #ndvi <- ndvi/10000
    #ndvi <- ifelse(ndvi <= 0 | ndvi > 1, NA, ndvi)
    ndvi <- ifelse(is.na(ndvi), 0, ndvi)
    res <- .C("bise", rdays = as.integer(days), ndvi = as.numeric(ndvi), 
        rslidperiod = as.integer(slidperiod), cndvi = as.numeric(correctedndvi), 
        PACKAGE = "biseVec")$cndvi
    res <- ifelse(res <= 0, NA, res)
    if (length(which(is.na(res) == FALSE)) < 5) {
        return(rep(NA, 7))
    }
    res <- rejectRapidIncrease(res)
	print(plot(res))
    f <- splinefun(x = 1:(3 * days), y = c(res, res, res), method = "monoH.FC")
    res <- f((days + 1):(2 * days))
    res[which(res < 0)] <- 0
    model <- res
    sosvec <- vector(mode = "integer", length = 7)
    sosvec[1] <- localThresGreenup(model, resvec[1])
    sosvec[2] <- localThresGreenup(model, resvec[2])
    sosvec[3] <- localThresGreenup(model, resvec[3])
    sosvec[4] <- globalThresGreenup(model, resvec[4])
    sosvec[5] <- globalThresGreenup(model, resvec[5])
    sosvec[6] <- min(na.omit(model))
    sosvec[7] <- max(na.omit(model))
    return(model)
}
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
make_mean_frame = function(output_frame,the_set){
output_frame=subset(output_frame,set==the_set)
output_frame = subset(output_frame,inc > 0 & dec > 0 & pr.peak > 0 & po.peak > 0)
output_frame$inc = output_frame$inc + output_frame$fs_start
output_frame$dec = output_frame$dec + output_frame$fs_start
output_frame$peak = output_frame$peak + output_frame$fs_start
output_frame$pr.peak = output_frame$pr.peak + output_frame$fs_start
output_frame$po.peak = output_frame$po.peak + output_frame$fs_start

output_frame$inc = output_frame$inc %% 365
output_frame$dec = output_frame$dec %% 365
output_frame$peak = output_frame$peak %% 365
output_frame$pr.peak = output_frame$pr.peak %% 365
output_frame$po.peak = output_frame$po.peak %% 365


output_frame$inc_c = cos(rad(output_frame$inc*(360/366)))
output_frame$dec_c = cos(rad(output_frame$dec *(360/366)))
output_frame$peak_c = cos(rad(output_frame$peak *(360/366)))
output_frame$pr.peak_c = cos(rad(output_frame$pr.peak *(360/366)))
output_frame$po.peak_c = cos(rad(output_frame$po.peak *(360/366)))

output_frame$inc_s = sin(rad(output_frame$inc*(360/366)))
output_frame$dec_s = sin(rad(output_frame$dec *(360/366)))
output_frame$peak_s = sin(rad(output_frame$peak *(360/366)))
output_frame$pr.peak_s = sin(rad(output_frame$pr.peak *(360/366)))
output_frame$po.peak_s = sin(rad(output_frame$po.peak *(360/366)))

inc_c = as.numeric(tapply(output_frame$inc_c,output_frame$clust,mean))
dec_c = as.numeric(tapply(output_frame$dec_c,output_frame$clust,mean))
peak_c = as.numeric(tapply(output_frame$peak_c,output_frame$clust,mean))
pr.peak_c = as.numeric(tapply(output_frame$pr.peak_c,output_frame$clust,mean))
po.peak_c = as.numeric(tapply(output_frame$po.peak_c,output_frame$clust,mean))

inc_c_sd = as.numeric(tapply(output_frame$inc_c,output_frame$clust,sd))
dec_c_sd = as.numeric(tapply(output_frame$dec_c,output_frame$clust,sd))
peak_c_sd = as.numeric(tapply(output_frame$peak_c,output_frame$clust,sd))
pr.peak_c_sd = as.numeric(tapply(output_frame$pr.peak_c,output_frame$clust,sd))
po.peak_c_sd = as.numeric(tapply(output_frame$po.peak_c,output_frame$clust,sd))

inc_s = as.numeric(tapply(output_frame$inc_s,output_frame$clust,mean))
dec_s = as.numeric(tapply(output_frame$dec_s,output_frame$clust,mean))
peak_s = as.numeric(tapply(output_frame$peak_s,output_frame$clust,mean))
pr.peak_s = as.numeric(tapply(output_frame$pr.peak_s,output_frame$clust,mean))
po.peak_s = as.numeric(tapply(output_frame$po.peak_s,output_frame$clust,mean))

inc_s_sd = as.numeric(tapply(output_frame$inc_s,output_frame$clust,sd))
dec_s_sd = as.numeric(tapply(output_frame$dec_s,output_frame$clust,sd))
peak_s_sd = as.numeric(tapply(output_frame$peak_s,output_frame$clust,sd))
pr.peak_s_sd = as.numeric(tapply(output_frame$pr.peak_s,output_frame$clust,sd))
po.peak_s_sd = as.numeric(tapply(output_frame$po.peak_s,output_frame$clust,sd))

inc=(deg(atan2(inc_s,inc_c)) %% 360) /(360/365)
dec=(deg(atan2(dec_s,dec_c)) %% 360) /(360/365)
peak=(deg(atan2(peak_s,peak_c)) %% 360) /(360/365)
pr.peak=(deg(atan2(pr.peak_s,pr.peak_c)) %% 360) /(360/365)
po.peak=(deg(atan2(po.peak_s,po.peak_c)) %% 360) /(360/365)

inc_sd=(deg(atan2(inc_s_sd,inc_c_sd)) %% 360) /(360/365)
dec_sd=(deg(atan2(dec_s_sd,dec_c_sd)) %% 360) /(360/365)
peak_sd=(deg(atan2(peak_s_sd,peak_c_sd)) %% 360) /(360/365)
pr.peak_sd=(deg(atan2(pr.peak_s_sd,pr.peak_c_sd)) %% 360) /(360/365)
po.peak_sd=(deg(atan2(po.peak_s_sd,po.peak_c_sd)) %% 360) /(360/365)


rows = as.numeric(rownames(tapply(output_frame$fs_start,output_frame$clust,mean)))
fs_start = as.numeric(tapply(output_frame$fs_start,output_frame$clust,mean))
inc_hs = as.numeric(tapply(output_frame$inc_hs,output_frame$clust,mean))
peak_hs = as.numeric(tapply(output_frame$peak_hs,output_frame$clust,mean))
peak_hs_sd = as.numeric(tapply(output_frame$peak_hs,output_frame$clust,sd))
len = as.numeric(tapply(output_frame$len,output_frame$clust,mean))
len_sd = as.numeric(tapply(output_frame$len,output_frame$clust,sd))
psm = as.numeric(tapply(output_frame$psm,output_frame$clust,mean))


mean_frame=data.frame(set=the_set,clust=rows,fs_start,inc,inc_sd,inc_hs,dec,dec_sd,peak,peak_sd,pr.peak,pr.peak_sd,po.peak,po.peak_sd,peak_hs,peak_hs_sd,len,len_sd,psm)
mean_frame
}
search_year=cmpfun(function(hs_vec,start){
	fire_season_increase=0
	fire_season_increase_hs=0
	fire_season_decrease=0
	peak_day=0
	peak_hs=0
	pre_season_mean=0
	pr.peak=0
	po.peak=0
	
	hs.max = max(hs_vec)
	hs.min = min(hs_vec)
	pc10 = hs.min + (hs.max-hs.min)*0.1
	pc80 = hs.min + (hs.max-hs.min)*0.8
	peak_day=which.max(hs_vec)
	peak_hs=hs_vec[peak_day]
	
	## Find broad boundaries
	print("Increase")
	for(test in peak_day:1){
		if(hs_vec[test] <= pc10)
		{
			fire_season_increase = test
			fire_season_increase_hs = hs_vec[test]
			print("Found")
			break
		}
	}
	print("Decrease")
	for(test in peak_day:366){
		if(hs_vec[test] <= pc10)
		{
			fire_season_decrease = test
			print("Found")
			break
		}
	}	
	print("Post")
	for(test in peak_day:366){
		if(hs_vec[test] <= pc80)
		{
			po.peak = test
			print("Found")
			break
		}
	}	
	print("Pre")
	for(test in peak_day:1){
		if(hs_vec[test] <= pc80)
		{
			pr.peak = test
			print("Found")
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
process_cluster=cmpfun(function(data,cluster_subset,bandw=10){
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


cluster_summary=data.frame(set=rep("A",length(start_points)),clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),pr.peak=rep(0,length(start_points)),po.peak=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")

for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){next}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	#g1=locpoly(r_tm,r_vec,bandwidth=20,drv=2,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	#g1=g1$y[(length(hs_vec)):(length(hs_vec)*2)]
	print(plot(g))
#	h=cdf(hs_vec)
#	i=seq_along(h)
#	h1=h/max(h)
#	j1=lm(logit(h1) ~ i)

#	j=nls(h ~ theta1/(1 + exp(-(theta2 + theta3*i))),start=list(theta1=max(h),theta2=j1$coefficients[1],theta3=j1$coefficients[2]))

	#plot(g,type="l",lwd=2)
	#points(hs_vec)
	
	out_vec=search_year(g)

	test_vec=out_vec
	
	cluster_summary$st_date[the_year] = as.character(start_date2)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=test_vec[[6]]
	cluster_summary$pr.peak[the_year]=test_vec[[7]]
	cluster_summary$po.peak[the_year]=test_vec[[8]]
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
A.cs=cluster_summary


print(paste("Subsetting cluster:",cluster_subset,sep=""))
cluster_hotspots = subset(data,data$Cluster==cluster_subset & pass=="N")
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


cluster_summary=data.frame(set=rep("N",length(start_points)),clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),pr.peak=rep(0,length(start_points)),po.peak=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")

for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){next}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	#g1=locpoly(r_tm,r_vec,bandwidth=20,drv=2,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	#g1=g1$y[(length(hs_vec)):(length(hs_vec)*2)]
	
#	h=cdf(hs_vec)
#	i=seq_along(h)
#	h1=h/max(h)
#	j1=lm(logit(h1) ~ i)

#	j=nls(h ~ theta1/(1 + exp(-(theta2 + theta3*i))),start=list(theta1=max(h),theta2=j1$coefficients[1],theta3=j1$coefficients[2]))

	#plot(g,type="l",lwd=2)
	#points(hs_vec)
	
	out_vec=search_year(g)

	test_vec=out_vec
	
	cluster_summary$st_date[the_year] = as.character(start_date2)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=test_vec[[6]]
	cluster_summary$pr.peak[the_year]=test_vec[[7]]
	cluster_summary$po.peak[the_year]=test_vec[[8]]
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
N.cs=cluster_summary


print(paste("Subsetting cluster:",cluster_subset,sep=""))
cluster_hotspots = subset(data,data$Cluster==cluster_subset & pass!="N")
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


cluster_summary=data.frame(set=rep("D",length(start_points)),clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),pr.peak=rep(0,length(start_points)),po.peak=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")

for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){next}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	#g1=locpoly(r_tm,r_vec,bandwidth=20,drv=2,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	#g1=g1$y[(length(hs_vec)):(length(hs_vec)*2)]
	
#	h=cdf(hs_vec)
#	i=seq_along(h)
#	h1=h/max(h)
#	j1=lm(logit(h1) ~ i)

#	j=nls(h ~ theta1/(1 + exp(-(theta2 + theta3*i))),start=list(theta1=max(h),theta2=j1$coefficients[1],theta3=j1$coefficients[2]))

	#plot(g,type="l",lwd=2)
	#points(hs_vec)
	
	out_vec=search_year(g)

	test_vec=out_vec
	
	cluster_summary$st_date[the_year] = as.character(start_date2)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=test_vec[[6]]
	cluster_summary$pr.peak[the_year]=test_vec[[7]]
	cluster_summary$po.peak[the_year]=test_vec[[8]]
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
D.cs=cluster_summary

total_f=rbind(A.cs,N.cs,D.cs)
total_f

})
process_cluster_frp=cmpfun(function(data,cluster_subset,bandw=10){
print(paste("Subsetting cluster:",cluster_subset,sep=""))
cluster_hotspots = subset(data,data$Cluster==cluster_subset)
hotspot_collapse=tapply(cluster_hotspots$frp,cluster_hotspots$jd,sum)
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
	num = sum(subf$frp)
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


cluster_summary=data.frame(set=rep("A",length(start_points)),clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),pr.peak=rep(0,length(start_points)),po.peak=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")

for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){next}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	#g1=locpoly(r_tm,r_vec,bandwidth=20,drv=2,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	#g1=g1$y[(length(hs_vec)):(length(hs_vec)*2)]
	print(plot(g))
#	h=cdf(hs_vec)
#	i=seq_along(h)
#	h1=h/max(h)
#	j1=lm(logit(h1) ~ i)

#	j=nls(h ~ theta1/(1 + exp(-(theta2 + theta3*i))),start=list(theta1=max(h),theta2=j1$coefficients[1],theta3=j1$coefficients[2]))

	#plot(g,type="l",lwd=2)
	#points(hs_vec)
	
	out_vec=search_year(g)

	test_vec=out_vec
	
	cluster_summary$st_date[the_year] = as.character(start_date2)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=test_vec[[6]]
	cluster_summary$pr.peak[the_year]=test_vec[[7]]
	cluster_summary$po.peak[the_year]=test_vec[[8]]
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
A.cs=cluster_summary


print(paste("Subsetting cluster:",cluster_subset,sep=""))
cluster_hotspots = subset(data,data$Cluster==cluster_subset & pass=="N")
hotspot_collapse=tapply(cluster_hotspots$frp,cluster_hotspots$jd,sum)
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

print("Filling individual days.")
for(the_date in 1:length(date_seq)){
	subf=subset(cluster_hotspots,cluster_hotspots$acq_date==date_seq[the_date])
		num = sum(subf$frp)
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


cluster_summary=data.frame(set=rep("N",length(start_points)),clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),pr.peak=rep(0,length(start_points)),po.peak=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")

for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){next}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	#g1=locpoly(r_tm,r_vec,bandwidth=20,drv=2,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	#g1=g1$y[(length(hs_vec)):(length(hs_vec)*2)]
	
#	h=cdf(hs_vec)
#	i=seq_along(h)
#	h1=h/max(h)
#	j1=lm(logit(h1) ~ i)

#	j=nls(h ~ theta1/(1 + exp(-(theta2 + theta3*i))),start=list(theta1=max(h),theta2=j1$coefficients[1],theta3=j1$coefficients[2]))

	#plot(g,type="l",lwd=2)
	#points(hs_vec)
	
	out_vec=search_year(g)

	test_vec=out_vec
	
	cluster_summary$st_date[the_year] = as.character(start_date2)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=test_vec[[6]]
	cluster_summary$pr.peak[the_year]=test_vec[[7]]
	cluster_summary$po.peak[the_year]=test_vec[[8]]
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
N.cs=cluster_summary


print(paste("Subsetting cluster:",cluster_subset,sep=""))
cluster_hotspots = subset(data,data$Cluster==cluster_subset & pass!="N")
hotspot_collapse=tapply(cluster_hotspots$frp,cluster_hotspots$jd,sum)
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

print("Filling individual days.")
for(the_date in 1:length(date_seq)){
	subf=subset(cluster_hotspots,cluster_hotspots$acq_date==date_seq[the_date])
		num = sum(subf$frp)
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


cluster_summary=data.frame(set=rep("D",length(start_points)),clust=rep(cluster_subset,length(start_points)),st_date=rep("",length(start_points)),fs_start=rep(fs_start,length(start_points)),inc=rep(0,length(start_points)),inc_hs=rep(0,length(start_points)),dec=rep(0,length(start_points)),peak=rep(0,length(start_points)),peak_hs=rep(0,length(start_points)),psm=rep(0,length(start_points)),pr.peak=rep(0,length(start_points)),po.peak=rep(0,length(start_points)),stringsAsFactors=F)

print("Analyzing years...")

for(the_year in 1:length(start_points)){
	print(paste("Year: ",the_year,sep=""))
	this_data = subset(full_hs,as.numeric(rownames(full_hs))>= start_points[the_year] & as.numeric(rownames(full_hs))< end_points[the_year])
	start_date2=this_data$date_seq[1]
	### Search for 3-day increase
	hs_vec=this_data$hs
	if(mean(hs_vec)==0){next}
	
	#Filter
	hs_vec=sgolayfilt(hs_vec)

	r_vec = c(hs_vec,hs_vec,hs_vec)
	r_tm = seq_along(r_vec)
	
	g=locpoly(r_tm,r_vec,bandwidth=bandw,gridsize=length(r_vec))
	#g1=locpoly(r_tm,r_vec,bandwidth=20,drv=2,gridsize=length(r_vec))
	g=g$y[(length(hs_vec)):(length(hs_vec)*2)]
	#g1=g1$y[(length(hs_vec)):(length(hs_vec)*2)]
	
#	h=cdf(hs_vec)
#	i=seq_along(h)
#	h1=h/max(h)
#	j1=lm(logit(h1) ~ i)

#	j=nls(h ~ theta1/(1 + exp(-(theta2 + theta3*i))),start=list(theta1=max(h),theta2=j1$coefficients[1],theta3=j1$coefficients[2]))

	#plot(g,type="l",lwd=2)
	#points(hs_vec)
	
	out_vec=search_year(g)

	test_vec=out_vec
	
	cluster_summary$st_date[the_year] = as.character(start_date2)
	cluster_summary$inc[the_year]=test_vec[[1]]
	cluster_summary$inc_hs[the_year]=test_vec[[2]]
	cluster_summary$dec[the_year]=test_vec[[3]]
	cluster_summary$peak[the_year]=test_vec[[4]]
	cluster_summary$peak_hs[the_year]=test_vec[[5]]
	cluster_summary$psm[the_year]=test_vec[[6]]
	cluster_summary$pr.peak[the_year]=test_vec[[7]]
	cluster_summary$po.peak[the_year]=test_vec[[8]]
}

cluster_summary$len= cluster_summary$dec - cluster_summary$inc
D.cs=cluster_summary

total_f=rbind(A.cs,N.cs,D.cs)
total_f

})
rad=cmpfun(function(deg){deg*pi/180})
deg=cmpfun(function(rad){rad*180/pi})
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


data <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\hotspots.csv",stringsAsFactors=F)
#data = subset(data,pass!="N",stringsAsFactors=F)

data$acq_date = as.Date(data$acq_date)
data$year=as.numeric(format(data$acq_date,format="%Y"))
data$jd=as.numeric(format(data$acq_date,format="%j"))

output_frame=data.frame(set=NA,clust=NA,st_date=NA,fs_start=NA,inc=NA,inc_hs=NA,dec=NA,peak=NA,peak_hs=NA,len=NA,psm=NA,pr.peak=NA,po.peak=NA,stringsAsFactors=F)

for(the_cluster in 1:61){
	frame=process_cluster(data,the_cluster)
	output_frame=rbind(output_frame,frame)
}


output_frame = subset(output_frame,!is.na(st_date))
write.csv(output_frame,"E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\complete_frame.csv")

output_frame=read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\complete_frame_frp.csv")
##################################################################
big=output_frame

op=make_mean_frame(big,"D") 
write.csv(op,"E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\mean_frame_frp_D.csv",row.names=F)

################################################################
################################################################

mean_frame=read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\mean_frame_frp_N.csv")

map=readShapePoly("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.shp")
map.rec=data.frame(map)

#m1 <- merge(map.rec, mean_frame, by.x = "GRIDCODE", by.y = "clust",all.x=T)
#plot(map,col=m1$peak_hs)

X11(width = 11, height = 8.5)
plot(map,col="gray90")
title("Mean FRP / Night Hotspots / 10% Season")

#Wide map
for(record in 1:length(map.rec$GRIDCODE)){
	the_rec=map.rec$GRIDCODE[record]
	co=map@polygons[[record]]@labpt
	if(record==51){co[1] = co[1] - 2}
	if(record==16){co[1] = co[1] + 2}
	ss=subset(mean_frame,mean_frame$clust==the_rec)
	if(length(ss$clust)==0){next}
	wedge_geo(start=ss$inc,end=ss$dec,mid=ss$peak,s_sd=ss$inc_sd,e_sd=ss$dec_sd,m_sd=ss$peak_sd,co=co)

}

X11(width = 11, height = 8.5)
plot(map,col="gray90")
title("Mean FRP / Night Hotspots / 80% Season")
#Narrow map
for(record in 1:length(map.rec$GRIDCODE)){
	the_rec=map.rec$GRIDCODE[record]
	co=map@polygons[[record]]@labpt
	if(record==51){co[1] = co[1] - 2}
	if(record==16){co[1] = co[1] + 2}
	ss=subset(mean_frame,mean_frame$clust==the_rec)
	if(length(ss$clust)==0){next}
	wedge_geo(start=ss$pr.peak,end=ss$po.peak,mid=ss$peak,s_sd=ss$pr.peak_sd,e_sd=ss$po.peak_sd,m_sd=ss$peak_sd/2,co=co)

}

########  Sequence Plots - defined time

mean_frame=read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\mean_frame_A.csv")
map=readShapePoly("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.shp")
map.rec=data.frame(map)
clusters = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.csv")

do_rect=function(start,end,col,rnk){
	if(start < end){
		rect(start, rnk-0.4, end, rnk+0.4,col=col)
	}
	if(start > end){
		rect(start, rnk-0.4, 365, rnk+0.4,col=col)
		rect(0, rnk-0.4, end, rnk+0.4,col=col)
	}
}


mean_frame$lat=0
for(idx in seq_along(mean_frame$clust)){
	lat=map@polygons[[idx]]@labpt[2]
	mean_frame$lat[idx]=lat
}

month_vec=c("J","F","M","A","M","J","J","A","S","O","N","D")
loc_vec=c(0,28,59,89,120,150,181,212,242,273,303,334)
par(mai=c(0.8,1.9,0.1,0.1),las=1)
plot(0,0,type="n",xlim=c(-10,375),ylim=c(1,62),xaxt="n",yaxt="n",xlab="",ylab="")
abline(v=172,lwd=2)
axis(1,loc_vec,month_vec,cex=0.8)
for(sl_no in 1:61){
	slice=mean_frame[sl_no,]
	place = which(order(mean_frame$lat)==sl_no)
	cluster_name = as.character(subset(clusters,GRIDCODE==sl_no)$Name[1])
	axis(2,place,cluster_name)
	do_rect(slice$inc,slice$dec,"grey80",place)
	do_rect(slice$pr.peak,slice$po.peak,"grey50",place)
	do_rect(((slice$peak)-3)%%365,((slice$peak)+3)%%365,"grey20",place)
}

############# Day/Night Split ###########

mean_frame=read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\mean_frame_frp_A.csv")
mean_frame_n = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\results\\Trial2\\mean_frame_frp_N.csv")
map=readShapePoly("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.shp")
map.rec=data.frame(map)
clusters = read.csv("E:\\NASA_NERP\\HotspotFFDI\\data\\ready\\clusters.csv")

do_rect=function(start,end,col,rnk){
	if(start < end){
		rect(start, rnk-0.4, end, rnk,col=col)
	}
	if(start > end){
		rect(start, rnk-0.4, 365, rnk,col=col)
		rect(0, rnk-0.4, end, rnk,col=col)
	}
}


mean_frame$lat=0
for(idx in seq_along(mean_frame$clust)){
	lat=map@polygons[[idx]]@labpt[2]
	mean_frame$lat[idx]=lat
}

month_vec=c("J","F","M","A","M","J","J","A","S","O","N","D")
loc_vec=c(0,28,59,89,120,150,181,212,242,273,303,334)
par(mai=c(0.8,2.1,0.1,0.1),las=1,cex=0.8)
plot(0,0,type="n",xlim=c(-10,375),ylim=c(1,62),xaxt="n",yaxt="n",xlab="",ylab="")
abline(v=172,lwd=2)
abline(v=0,lwd=1)
abline(v=365,lwd=1)
axis(1,loc_vec,month_vec,cex=0.8)
for(sl_no in 1:61){
	slice=mean_frame[sl_no,]
	slice_n=mean_frame_n[sl_no,]
	place = which(order(mean_frame$lat)==sl_no)
	cluster_name = as.character(subset(clusters,GRIDCODE==sl_no)$Name[1])
	axis(2,place,cluster_name,cex=0.8)
	
	do_rect(slice$inc,slice$dec,"burlywood1",place)
	do_rect(slice$pr.peak,slice$po.peak,"burlywood3",place)
	do_rect(((slice$peak)-3)%%365,((slice$peak)+3)%%365,"chocolate4",place)
	
	do_rect(slice_n$inc,slice_n$dec,"aliceblue",place+0.4)
	do_rect(slice_n$pr.peak,slice_n$po.peak,"cadetblue2",place+0.4)
	do_rect(((slice_n$peak)-3)%%365,((slice_n$peak)+3)%%365,"blue",place+0.4)
}


