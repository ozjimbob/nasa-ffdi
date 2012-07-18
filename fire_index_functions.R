#######################
#######
#######   Fire Index Functions
#######
#######################



# Support Functions

#Load table needed for KBDI Lookup

#kbdilookup <- read.csv("E:\\NASA_NERP\\HotspotFFDI\\NASA_NERP_Hotspot FFDI\\kbdi.csv",header=T)

radians <- function(deg){
	deg*(pi/180)
}

degrees <- function(rad){
	rad*(180/pi)
}

#Calculate day length

sunaltitude <- function(latitude,dayofyear,time){
	decT = 2.0 * pi * ((dayofyear-1.0)/365.0)
	declination = 0.322003-22.971*cos(decT)-0.357898*cos(2*decT)-0.14398*cos(3*decT)+3.94638*sin(decT)+0.019334*sin(2*decT)+0.05928*sin(3*decT)
	hourangle = radians((15*(time-12.0)))
	altT = sin(radians(declination))*sin(radians(latitude))
	altU = cos(radians(declination))*cos(radians(latitude))*cos(hourangle)
	altcalc = asin(altT + altU)
	altitude = degrees(altcalc)
	altitude
}

dayintegral <- function(dayofyear,latitude){
	solarint = 0.0
	for(hour in 0:384){
		altitude = sunaltitude(latitude,dayofyear,hour/16)
		if(altitude > 0.0){
			solarint = solarint + 1/16
		}
	}
	solarint
}

hourwatt <- function(dayofyear,latitude,hour){
	altitude_deg=sunaltitude(latitude,dayofyear,hour)
	if(altitude_deg <= 0){
		radi = 0
	}else{
		GetAirMassRatio = (1/sin(radians(altitude_deg)))
		GetApparentExtraterrestrialFlux=1160 + (75 * sin((360/365) * (dayofyear - 275)))
		GetOpticalDepth = 0.174 + (0.035 * sin((360/365) * (dayofyear - 100)))
		radi = GetApparentExtraterrestrialFlux * exp(-1 * GetOpticalDepth * GetAirMassRatio)
	}
	radi
}


DayLengths <- function(latitude){
	daylist=c()
	for(a in 1:366){
		daylist = append(daylist,dayintegral(a,latitude))
	}
	daylist
}



calc_dew<-function(T,RH){
	RH[RH < 1] = 1
	dew_numer = 237.3*log(RH/100*6.112/6.1078*exp((17.67*T)/(T-0+243.5)))
	dew_denom = 17.27-log(RH/100*6.112/6.1078*exp((17.67*T)/(T-0+243.5)))
	H <- dew_numer/dew_denom
	H
}



calc_humid<-function(T,Dp){
	satvp <- 6.11*10^((7.5*T/(237.7+T)))
	actvp = 6.11*10^((7.5*Dp/(237.7+Dp)))
	H <- (actvp/satvp)*100
	H
}
      
#DF Lookup Function for KBDI

lookupDF<- function(DI,T,MAP,klook){
	#attach(kbdilookup)
	#Old Method
	if(T<1){T<-1}
	for(s in 1:length(klook$MinMAP)){
		if(T >= klook$MinTemp[s] & T < klook$MaxTemp[s]+1 & DI >= klook$MinDI[s] & DI < klook$MaxDI[s]+1 & MAP >= klook$MinMAP[s] & MAP < klook$MaxMAP[s]+1){
   			idx<-klook$DF[s]
			break
		}
	}
idx
}

lookupET<-function(pKBDI,T,MAP){
	ET = 0.001 * ((203.2-pKBDI)*(0.968*exp(0.0875*T+1.5552)-8.3) / (10.88*exp(-0.001736*MAP)+1))
	ET
}

#K Lookup Function

lookK <- function(Rain){
    if (Rain >= 0.0 & Rain < 0.1){
        K=1.0}
    else if (Rain >=0.1 & Rain <1.0){
        K=0.8}
    else if (Rain >=1.0 & Rain <3.0){
        K=0.6}
    else if (Rain >=3.0 & Rain <6.0){
        K=0.4}
    else if (Rain >=6.0 & Rain < 15.0){
        K=0.2}
    else if (Rain >=15.0 & Rain < 19.0){
        K=0.1}
    else{
        K = 0.0}
	K
}

# New calculation of DF, from rain and Soil Moisture Deficit


DF_new <- function(rain,SMD_vec){
	conseq_days_rain = rep(0,length(rain))
	conseq_sum_rain = rep(0,length(rain))
	size_of_largest=rep(0,length(rain))
	days_since_largest=rep(0,length(rain))
	rain[is.na(rain)]=0
	conseq = 0
	conseq_sum = 0

	rfvec <- rain

	for(d in 2:length(rfvec)){
		prev_rain = rfvec[d-1]
		cur_rain = rfvec[d]
		if(cur_rain > 2){
			conseq = conseq + 1
			conseq_sum = conseq_sum + cur_rain
		}
	
		if(cur_rain <= 2 & conseq != 0){
			conseq = 0
			conseq_sum = 0
		}

		biggest = max(rfvec[(d-conseq):d])
		which_day = which.max(rev(rfvec[(d-conseq):d]))-1
		conseq_days_rain[d] = conseq
		size_of_largest[d] = biggest

		if(conseq > 0 ){
			days_since_largest[d] = which_day
		}
		if(conseq == 0){
			days_since_largest[d] = days_since_largest[d-1]+1
		}

		if(conseq > 0 ){
			conseq_sum_rain[d] = conseq_sum
		}
		if(conseq == 0){
			conseq_sum_rain[d] = conseq_sum_rain[d-1]
		}	

	}

	look_frame = data.frame(rfvec,conseq_days_rain,conseq_sum_rain,size_of_largest,days_since_largest)
	look_frame$conseq_sum_rain[look_frame$days_since_larges > 20] = 0
	
	
	N_vec = look_frame$days_since_largest
	P_vec = look_frame$conseq_sum_rain

	x_vec = rep(0,length(P_vec))

	for(a in 1:length(N_vec)){
		N = N_vec[a]
		P = P_vec[a]

		if(N >= 1 & P >= 2){
			x = N^1.3/(N^1.3 + P - 2)
		}
		if(N == 0 & P >= 2){
			x = 0.7481988/(0.7481988 + P - 2)
		}
		if(P < 2){
			x = 1
		}
		x_vec[a] = x
	
	}
	
	xlim_vec = rep(0,length(SMD_vec))



	for(a in 1:length(SMD_vec)){
		SMD = SMD_vec[a]
		if(is.na(SMD)){
			next
		}
		if(SMD < 20){
			xlim = 1/(1+0.1135*SMD)
		}
		if(SMD >= 20){
			xlim = 75/(270.525-1.267*SMD)
		}
		if(x_vec[a] > xlim){
			x_vec[a] = xlim
		}
	}

	DF_vec = rep(0,length(SMD_vec))
	for(a in 1:length(DF_vec)){
		SMD = SMD_vec[a]
		x = x_vec[a]
		if(is.na(SMD)){
			next
		}
		DF = 10.5*(1-exp(-(SMD+30)/40))*(41*x^2+x)/(40*x^2+x+1)
		DF_vec[a] = DF
	}
	DF_vec
}


#McArthur Fire Index

index_MA <- function(Temperature, Rain, DewPoint,Wind, KBDI){
	McA<-0.0
	McA_list <- c()
	llength<-length(Temperature)
	DSinceRain<-0
	LastRain<-0.0
	for(day in 1:llength) {
		if(Rain[day] > 0.0){
			LastRain <- Rain[day]
			DSinceRain <- 0
		}else{
			DSinceRain <- DSinceRain + 1
		}


		satvp <- 6.11*10^((7.5*Temperature[day]/(237.7+Temperature[day])))
		actvp = 6.11*10^((7.5*DewPoint[day]/(237.7+DewPoint[day])))
		Humid <- (actvp/satvp)*100

		DF <- (0.191 * (KBDI[day] + 104.0)*(DSinceRain+1)**1.5) / ((3.52*(DSinceRain+1)**1.5)+LastRain-1.0)
		if(DF > 10.0){
			DF <- 10
		}

		#F <- 1.25 * DF * exp(((Temperature[day] - Humid)/30)+0.0234*Wind[day])
		F <- 2 * exp(-.45 + .987 * log(DF + .001) - .0345 * Humid + .0338 * Temperature[day] + .0234 * Wind[day]) 
		McA_list <- append(McA_list,F)
	}
	McA_list
}



#McArthur Fire Index From Drought Factor

index_MA_DF <- function(Temperature, DF, Humidity, Wind){
	McA<-0.0
	McA_list <- c()
	llength<-length(Temperature)
	for(day in 1:llength) {
		#F <- 1.25 * DF * exp(((Temperature[day] - Humidity[day])/30)+0.0234*Wind[day])
		F <- 2 * exp(-.45 + .987 * log(DF[day] + .001) - .0345 * Humidity[day] + .0338 * Temperature[day] + .0234 * Wind[day]) 
		McA_list <- append(McA_list,F)
	}
	McA_list
}


#McArthur Grassland Mk4

index_G4 <- function(Temperature, DewPoint,Wind){
	McA <- 0.0
	McA_list <- c()
	llength<-length(Temperature)
	for(day in 1:llength){
		satvp <- 6.11*10^((7.5*Temperature[day]/(237.7+Temperature[day])))
		actvp = 6.11*10^((7.5*DewPoint[day]/(237.7+DewPoint[day])))
		Humid <- (actvp/satvp)*100
		H <- Humid
		T <- Temperature[day]
		U <- Wind[day]
		C = 100
		Q = 4.5
		F = exp(-1.523 + 1.027 * log(Q)-0.009432*(100-C)^1.535+0.02764*T-0.2205*sqrt(H)+0.6422*sqrt(U))
		McA_list <- append(McA_list,F)
	}
	McA_list
}





# Keetch-Byran Drought Index

index_KBDI <- function(Temperature, Rain, MAP,klook){
	MAP <- MAP / 25.4
	Rain <- Rain / 25.4
	KBDI <- 400.0
	rainYD <- 0.0
	KBDI_list <- c()
	llength <- length(Temperature)
	for(day in 1:llength){
		pKBDI <- KBDI
		if(rainYD == 0.0){
			realrain <- Rain[day] - 0.2
		} else {
			realrain <- Rain[day]
		}
		rainYD <- Rain[day]
		if(realrain < 0.0){
			realrain <- 0.0
		}
		KBDI <- KBDI - (realrain * 100)
		if(KBDI < 0.0){
			KBDI <- 0.0
		}
		FTemp <- (Temperature[day] * (9.0/5.0)) + 32.0
		todayDF <- lookupDF(KBDI,FTemp,MAP,klook)
		todayET <- lookupET(pKBDI,Temperature[day],MAP*25.4)
		KBDI <- KBDI + todayDF
		#KBDI <- KBDI + todayET
		KBDI_list <- append(KBDI_list,KBDI)
	}
	KBDI_list
}


index_KBDI_new <- function(Temperature, Rain, MAP,klook){
		
	KBDI <- 200.0
	KBDI_list <- c()
	llength <- length(Temperature)
	for(day in 1:llength){
		pKBDI <- KBDI
		KBDI <- KBDI - Rain[day]
		if(KBDI < 0.0){
			KBDI <- 0.0
		}
		todayET <- lookupET(pKBDI,Temperature[day],MAP)
		KBDI <- KBDI + todayET
		KBDI_list <- append(KBDI_list,KBDI)
	}
	KBDI_list
}




# Nesterov Index

index_NI<-function(DewPoint,Temperature,Rain){
	Deficit <- Temperature - DewPoint
	G<-0.0
	G_list <- c()
	llength<-length(DewPoint)
	for(day in 1:llength){
		if (Rain[day] > 3.0){
			G <- 0.0
			G_list <- append(G_list,G)
			next
		}
		if(Temperature[day] > 0.0){
			G <- G + (Temperature[day] * Deficit[day])
		}
		G_list <- append(G_list,G)
	}
	G_list
}





# Zdenkho Index

index_Z <- function(DewPoint,Temperature,Rain){
	Deficit <- Temperature - DewPoint
	Z <- 0.0
	Z_list <- c()
	llength <- length(DewPoint)
	for(day in 1:llength){
		K <- lookK(Rain[day])
		Z <- (Z + Deficit[day]) * K
		Z_list <- append(Z_list,Z)
	}
	Z_list
}




#Modified Nesterov

index_MMI <- function(DewPoint,Temperature,Rain){
	Deficit <- Temperature - DewPoint
	M <- 0.0
	M_list <- c()
	llength <- length(DewPoint)
	for(day in 1:llength){
		if (Rain[day] > 3.0){
			M <- 0.0
			M_list <- append(M_list,M)
			next
		}
		if(Temperature[day] > 0.0){
			M <- M + (Temperature[day] * Deficit[day] * lookK(Rain[day]))
		}
		M_list <- append(M_list,M)
	}
	M_list
}




    
# Sharples Fire Weather Index - With and Without KBDI

index_Fix <- function(Temp,Wind,Dew){
	humid <- calc_humid(Temp,Dew)
	FFI_list <- c()
	llength <- length(Temp)
	for(day in 1:llength){
		FMI = 10 - 0.25*(Temp[day]-humid[day])
		f = max(Wind[day],1)/FMI
		FFI_list <- append(FFI_list,f)
	}
	FFI_list
}


index_FixD <- function(Temp,Wind,Dew,KBDI){
	humid <- calc_humid(Temp,Dew)
	FFI_list <- c()
	llength <- length(Temp)
	for(day in 1:llength){
		FMI = 10 - 0.25*(Temp[day]-humid[day])
		f = KBDI[day]*(max(Wind[day],1)/FMI)
		FFI_list <- append(FFI_list,f)
	}
	FFI_list
}




##########Fosberg

index_FOS <- function(Temp,Wind,Dew){
	humid <- calc_humid(Temp,Dew)
	Wind <- Wind / 1.6
	Temp = Temp * (9/5) + 32
	FFI_list <- c()
	llength <- length(Temp)
	for(day in 1:llength){
		if(is.na(humid[day])){
			FFI_list <- append(FFI_list,NA)
		} else {
			if(humid[day] < 10){
				m=0.03229 + 0.281073 * humid[day] - 0.000578 * humid[day] * Temp[day]
			}
			if(humid[day] >= 10 & humid[day] <= 50){
				m=2.22749 + 0.160107 * humid[day] - 0.01478 * Temp[day]
			}
			if(humid[day] > 50){
				m=21.0606 + 0.005565 * humid[day]^2 - 0.00035 * humid[day] * Temp[day] - 0.483199 * humid[day] 
			}
			n=1-2*(m/30)+1.5*(m/30)^2-0.5*(m/30)^3
			FFWI = n*((1+Wind[day]^2)^.5)/0.3002
			FFI_list <- append(FFI_list,FFWI)
		}
	}
	FFI_list
}





#Canadian Fire Weather Index

index_CFWI <- function(Month,Days,Temp,Dew,WS,Rain,daylist){
	Humid <- calc_humid(Temp,Dew)
      #print(Humid)
	#DLA <- c(6.4,5.0,2.4,0.4,-1.6,-1.6,-1.6,-1.6,-1.6,0.9,3.8,5.8)
	DLA <- c(0,0,0,0,0,0,0,0,0,0,0,0)
	CFWI <- 0.0
	F <- 85.0
	P <- 6.0
	D <- 15.0
	CFWI_list <- c()
	llength <- length(Temp)
	for(day in 1:llength){
		
		T = Temp[day]
		H = Humid[day]
		#print(day)
		if(H > 100){H=100}
		W = WS[day]
		if(W < 0){W=0}
		ro = Rain[day]
		if(ro < 0){ro = 0}
		Le = daylist[Days[day]]

		#Fine Fuel Moisture Code
		Fo = F
		#Calculate mo
		mo = 147.2 * (101 - Fo) / (59.5 + Fo)

		#Adjust mo for rain
		if(ro > 0.5){
			rf = ro - 0.5
			if(mo <= 150){
				mr = mo + 42.5 * rf * (exp(-100/(251-mo)))*(1-exp(-6.93/rf))
			} else {	
				mr = mo + 42.5 * rf * (exp(-100/(251-mo)))*(1-exp(-6.93/rf))+0.0015*((mo-150)^2)*(rf^0.5)
			}
			if (mr > 250){mr = 250}
			mo = mr
		}
		Ed = 0.942 * H^0.679 + 11 * exp((H-100)/10) + 0.18*(21.1-T)*(1-exp(-0.115 * H))
	
		if(mo > Ed){
			ko = 0.424 * (1-(H/100)^1.7) + 0.0694 * W^0.5 * (1-(H/100)^8)
			kd = ko * 0.581 * exp(0.0365 * T)
			m = Ed + (mo - Ed) * 10 ^ (-kd)
		}


		Ew = 0.618 * H^0.753 + 10 * exp((H-100)/10) + 0.18*(21.1-T)*(1-exp(-0.115 * H))

		if(mo < Ew & mo < Ed){
			kl = 0.424 * (1-((100-H)/100)^1.7) + 0.0694 * W ^ 0.5 * (1-((100/H)/100)^8)
			kw = kl * 0.581 * exp(0.0365 * T)
			m = Ew - (Ew - mo) * 10 ^ (-kw)
		}	


		if ( Ed >= mo & mo >= Ew){
			m = mo
		}

		F = 59.5 * (250 - m) / (147.2 + m)

		#Duff Moisture Code
		Po = P
		if (ro > 1.5){
			re = 0.92 * ro - 1.27
			Mo = 20 + exp(5.6348 - Po / 43.43)
			if(Po <= 33){
				b = 100 / (0.5 + 0.3 * Po)
			}
			if(Po > 33 & Po <= 65){
				b = 14 - 1.3 * log(Po)
			}
			if(Po > 65){
				b = 6.2 * log(Po) - 17.2
			}
			Mr = Mo + 1000 * re / (48.77 + b * re)
			Pr = 244.72 - 43.43 * log (Mr - 20)
			if (Pr < 0){Pr = 0}
			Po = Pr
		}
		
		if (T < -1.1){
			TT = -1.1
		} else {
			TT = T
		}

		K = 1.894 * (TT + 1.1) * (100 - H) * Le * 10E-6
		P = Po + 100 * K

		# Drought Code
		Do = D
		if (ro > 2.8){
			rd = 0.83 * ro - 1.27
			Qo = 800 * exp(-Do/400)
			Qr = Qo + 3.937 * rd
			Dr = 400 * log(800 / Qr)
			if (Dr < 0){Dr = 0}
			Do = Dr
		}
		
		if (T < 2.8){
			TT = 2.8
		} else {
			TT = T
		}
				
		V = 0.36 * (TT + 2.8) 
		D = Do + 0.5 * V

		# Initial Spread Index

		fW = exp(0.05039 * W)
		fF = 91.9 * exp(-0.1386 * m)*(1+m^5.31/(4.93 * 10E7))
		R = 0.208 * fW * fF


		# Buildup Index
		if ( P <= 0.4 * D){
			U = 0.8 * P * D / (P + 0.4 * D)
		}
		if ( P > 0.4 * D){
			U = P - (1-0.8*D/(P+0.4 * D))*(0.92 + (0.0114 * P)^1.7)
		}

		#Fire Weather Index
		if ( U <= 80){
			fD = 0.626 * U^0.809 + 2
		}
		if ( U > 80){
			fD = 1000/(25 + 108.64 * exp(-0.023 * U))
		}
		B = 0.1 * R * fD
		if(is.na(B)){
			S = 0 
		}else{
			if (B > 1){
				S = exp(2.72 * (0.434 * log(B))^0.647)
			}
			if (B <= 1){
				S = B
			}
		}
		
		
		CFWI_list = append(CFWI_list,S)
	}
	CFWI_list
}







# Calculation of Drought Factor from rain and soil moisture deficit
# DF_new <- function(Rain,SoilMoistureDeficit)

# McArthur Forest Fire Danger Index
# index_MA <- function(Temperature, Rain, DewPoint,Wind, KBDI)

# McArthur Forest Fire Danger Index with precalculated Drought Factor
# index_MA_DF <- function(Temperature, DF, Humidity, Wind)

# McArthur Grassland Fire Danger Index
# index_G4 <- function(Temperature, DewPoint,Wind)

# Keech-Byran Drought Index
# index_KBDI_new <- function(Temperature, Rain, MAP)

# Nesterov Index
# index_NI<-function(DewPoint,Temperature,Rain)

# Zdenkho Index
# index_Z <- function(DewPoint,Temperature,Rain)

# Modified Nesterov Index
# index_MMI <- function(DewPoint,Temperature,Rain)

# Sharples Index
# index_Fix <- function(Temp,Wind,Dew)

# Sharples Index with KBDI
# index_FixD <- function(Temp,Wind,Dew,KBDI)

# Fosberg Index
# index_FOS <- function(Temp,Wind,Dew)

# Canadian Fire Weather Index
# index_CFWI <- function(Month,Days,Temp,Dew,WS,Rain,daylist)




#Generate list of day lengths for this latitude





