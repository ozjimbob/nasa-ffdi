require(raster)

root<<-"G:\\SILO_Data\\"

SILO_Get_Raster = function(variable,thedate,draw.plot=F){
	date_vec=strsplit(thedate, "-")[[1]]
	year = date_vec[1]
	the_file = paste(date_vec[1],date_vec[2],date_vec[3],"_",variable,".txt",sep="")
	fullpath = paste(root,variable,"\\",year,"\\",the_file,sep="")
	print(fullpath)
	inputraster = raster(fullpath)
	if(draw.plot){
		image(inputraster)
	}
	inputraster
}


SILO_Get_Value = function(variable,thedate,x,y){
	date_vec=strsplit(thedate, "-")[[1]]
	year = date_vec[1]
	the_file = paste(date_vec[1],date_vec[2],date_vec[3],"_",variable,".txt",sep="")
	fullpath = paste(root,variable,"\\",year,"\\",the_file,sep="")
	print(fullpath)
	inputraster = raster(fullpath)
	the_cell=cellFromXY(inputraster,c(x,y))
	the_value=extract(inputraster,the_cell)
	the_value
}


SILO_Get_Station=function(station){
	samplefile = paste("G:\\SILO_Data\\Patched_Point_Data\\",station,sep="")
	ppd_cols = c("stnum","year","month","day","t_max","t_max_s","t_min","t_min_s","rain","rain_s","pan","pan_s","rad","rad_s","vp","vp_s","rh_max","rh_min","fao65","m_lake","m_pot","m_act","m_wet","s_pan","pan_ws","pan_ws_s","pres","pres_s")
	col_sep=c(0,6,10,12,14,19,21,26,28,34,36,41,43,47,49,53,55,60,65,69,73,77,81,85,90,93,95,101,103)
	out_data=read.fwf(samplefile,widths=c(diff(col_sep)),header=F,skip=1,col.names=ppd_cols)
	out_data$c_date = as.Date(paste(out_data$year,out_data$month,out_data$day,sep="-"))
	out_data
}

