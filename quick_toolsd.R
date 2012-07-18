###  Create smaller met archive files, for quick loading

inf="G:\\SILO_Data\\complete_stations\\"
out="G:\\SILO_Data\\complete_stations_min\\"

files = list.files(inf)
itr=0
for(y in files){
		itr=itr+1
		print(y)
		print(itr)
		get.it = read.csv(paste(inf,y,sep=""))
		get.it=get.it[36361:40661,]
		write.csv(get.it,paste(out,y,sep=""))
	}