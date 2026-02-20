source("FindPhases.R")

shortest = function(values,gamma){
# finds shortest segment of the vector values whose sum is at least the proportion gamma of total = sum(values)
# returns left and right extremities of such segment

	n = length(values)

	# cumul[i] = sum_{j < i} values[j] / total, for i=1,...,n+1
	cumul = c(0,cumsum(values)/sum(values))

	# maximum possible left extremity of segment
	maxleft = max(which(cumul[1:n] <= 1-gamma))

	# lengths of shortest segments containing at least proportion gamma of total with given left extremities
	len = sapply(1:maxleft,function(i) min(which(cumul[(i+1):(n+1)] >= cumul[i] + gamma)))

	left = which.min(len) # left extremity of the shortest possible segment

	return(list(left=left,right=left+len[left]-1))
}

shortest_date = function(dates,values,gamma=0.99){
# returns elements from vector dates corresponding to the extremities
# of the segment found by shortest(values,gamma)

	s = shortest(values,gamma)
	return(list(left=dates[s$left],right=dates[s$right]))
}

extract_beginning = function(values,gamma=0.1){
# extracts from the vector values a initial segment whose size is gamma times the size of the entire vector

	sizeseg = max(round(length(values)*gamma,0),1)
	return(values[1:sizeseg])
}

extract_peak = function(values,valuesref,gamma=0.1){
# extracts from the vector values a segment whose size is gamma times the size of the entire vector
# the segment starts at the peak of the vector valuesref

	sizeseg = max(round(length(values)*gamma,0),1)
	peak = which.max(valuesref)
	right = min(peak+sizeseg-1,length(values))
	return(values[peak:right])
}


# finds location of waves (shortest interval contained in the phase which contain at least 99% of deaths)

beginwave1 = as.Date(integer(0),origin = "1970-01-01")
endwave1 = as.Date(integer(0),origin = "1970-01-01")
beginwave2 = as.Date(integer(0),origin = "1970-01-01")
endwave2 = as.Date(integer(0),origin = "1970-01-01")
beginwave3 = as.Date(integer(0),origin = "1970-01-01")
endwave3 = as.Date(integer(0),origin = "1970-01-01")

for(country in countries)
{
	country_phase1 = dataset[dataset$location_name == country & dataset$date <= endphase1,]
	country_phase2 = dataset[dataset$location_name == country & dataset$date > endphase1 & dataset$date <= endphase2,]
	country_phase3 = dataset[dataset$location_name == country & dataset$date > endphase2,]

	wave1 = shortest_date(country_phase1$date,country_phase1$deathcovid)
	wave2 = shortest_date(country_phase2$date,country_phase2$deathcovid)
	wave3 = shortest_date(country_phase3$date,country_phase3$deathcovid)

	beginwave1 = c(beginwave1,wave1$left)
	endwave1 = c(endwave1,wave1$right)
	beginwave2 = c(beginwave2,wave2$left)
	endwave2 = c(endwave2,wave2$right)	
	beginwave3 = c(beginwave3,wave3$left)
	endwave3 = c(endwave3,wave3$right)	
}

# data frame containing one row per wave (3 waves x 24 countries = 72 rows)
dataperwave = data.frame(location_name=rep(countries,3),wave=rep(1:3,each=length(countries)),
			begin=c(beginwave1,beginwave2,beginwave3),end=c(endwave1,endwave2,endwave3))
dataperwave = dataperwave[order(dataperwave$location_name,dataperwave$wave),]

# adds columns to dataset with information about beginning and ending of each wave
# new column wave is the number of the wave corresponding to the given country and date (zero means interwave period)
datasetwaves = data.frame(location_name=countries,beginwave1,endwave1,beginwave2,endwave2,beginwave3,endwave3) |> left_join(dataset,y=_,by="location_name")
datasetwaves = datasetwaves |> mutate(wave = ifelse(date >= beginwave1 & date <= endwave1,1,0)+
		ifelse(date >= beginwave2 & date <= endwave2,2,0)+ifelse(date >= beginwave3 & date <= endwave3,3,0))

# computes mean Covid-19 mortality per million for each wave
dataperwave = datasetwaves |> group_by(location_name,wave) |> summarize(mean_death=mean(deathcovid),.groups="drop") |> left_join(dataperwave,y=_,by=c("location_name","wave"))


# finds beginning of waves (for index maskbeginwave), i.e., the decile of the wave starting at the beginning of the wave
# finds periods near peak of waves (for index maskpeakwave), i.e., the decile of the wave starting at the peak of mortality

end_beginning = as.Date(integer(0),origin = "1970-01-01")
begin_peak = as.Date(integer(0),origin = "1970-01-01")
end_peak = as.Date(integer(0),origin = "1970-01-01")
for(i in 1:nrow(dataperwave))
{
	dataset_countrywave = datasetwaves[datasetwaves$location_name == dataperwave$location_name[i] & datasetwaves$wave == dataperwave$wave[i],]
	end_beginning = c(end_beginning,max(extract_beginning(dataset_countrywave$date)))
	begin_peak = c(begin_peak,min(extract_peak(dataset_countrywave$date,dataset_countrywave$deathcovid)))
	end_peak = c(end_peak,max(extract_peak(dataset_countrywave$date,dataset_countrywave$deathcovid)))
}
dataperwave$end_beginning = end_beginning
dataperwave$begin_peak = begin_peak
dataperwave$end_peak = end_peak

if(TRUE) # you can replace "TRUE" with "FALSE" if you don't want to generate these graphs
{

# graphs (one per country) showing Covid-19 mortality per million, mask usage and waves annotated in aquamarine transparent rectangles

pdf("Waves.pdf",width=14.23,height=6.77)
for(i in 1:length(countries))
{
	print(ggplot(data=dataset[dataset$location_name == countries[i],])+
		geom_line(aes(x=date,y=deathcovid,col="deaths"))+
		geom_line(aes(x=date,y=mask_use_mean*10,col="masks"))+
		scale_x_continuous(breaks=xaxis_br,labels=xaxis_lab)+
		scale_color_manual(values=c("deaths"="black","masks"="purple"),labels=c("deaths"="Covid-19 deaths per million","masks"="mask usage times 10"),name=NULL)+
		labs(title=countries[i],y=NULL)+
		annotate("rect",xmin=beginwave1[i],xmax=endwave1[i],ymin=-Inf,ymax=Inf,fill="aquamarine",alpha=0.3)+
		annotate("rect",xmin=beginwave2[i],xmax=endwave2[i],ymin=-Inf,ymax=Inf,fill="aquamarine",alpha=0.3)+
		annotate("rect",xmin=beginwave3[i],xmax=endwave3[i],ymin=-Inf,ymax=Inf,fill="aquamarine",alpha=0.3))
}
dev.off()

# graphs like in Waves.pdf, but including yellow transparent rectangles indicating beginning of waves (for index maskbeginwave)
# and orange transparent rectangles indicating period near peak of waves (for index maskpeakwave)

pdf("WavesMarked.pdf",width=14.23,height=6.77)
for(i in 1:length(countries))
{
	dataperwave_country = dataperwave[dataperwave$location_name == countries[i],]
	endbeginning1 = dataperwave_country$end_beginning[dataperwave_country$wave == 1]
	endbeginning2 = dataperwave_country$end_beginning[dataperwave_country$wave == 2]
	endbeginning3 = dataperwave_country$end_beginning[dataperwave_country$wave == 3]
	beginpeak1 = dataperwave_country$begin_peak[dataperwave_country$wave == 1]
	beginpeak2 = dataperwave_country$begin_peak[dataperwave_country$wave == 2]
	beginpeak3 = dataperwave_country$begin_peak[dataperwave_country$wave == 3]
	endpeak1 = dataperwave_country$end_peak[dataperwave_country$wave == 1]
	endpeak2 = dataperwave_country$end_peak[dataperwave_country$wave == 2]
	endpeak3 = dataperwave_country$end_peak[dataperwave_country$wave == 3]

	print(ggplot(data=dataset[dataset$location_name == countries[i],])+
		geom_line(aes(x=date,y=deathcovid,col="deaths"))+
		geom_line(aes(x=date,y=mask_use_mean*10,col="masks"))+
		scale_x_continuous(breaks=xaxis_br,labels=xaxis_lab)+
		scale_color_manual(values=c("deaths"="black","masks"="purple"),labels=c("deaths"="Covid-19 deaths per million","masks"="mask usage times 10"),name=NULL)+
		labs(title=countries[i],y=NULL)+
		annotate("rect",xmin=beginwave1[i],xmax=endwave1[i],ymin=-Inf,ymax=Inf,fill="aquamarine",alpha=0.3)+
		annotate("rect",xmin=beginwave2[i],xmax=endwave2[i],ymin=-Inf,ymax=Inf,fill="aquamarine",alpha=0.3)+
		annotate("rect",xmin=beginwave3[i],xmax=endwave3[i],ymin=-Inf,ymax=Inf,fill="aquamarine",alpha=0.3)+
		annotate("rect",xmin=beginwave1[i],xmax=endbeginning1,ymin=-Inf,ymax=Inf,fill="yellow",alpha=0.2)+
		annotate("rect",xmin=beginwave2[i],xmax=endbeginning2,ymin=-Inf,ymax=Inf,fill="yellow",alpha=0.2)+
		annotate("rect",xmin=beginwave3[i],xmax=endbeginning3,ymin=-Inf,ymax=Inf,fill="yellow",alpha=0.2)+
		annotate("rect",xmin=beginpeak1,xmax=endpeak1,ymin=-Inf,ymax=Inf,fill="orange",alpha=0.2)+
		annotate("rect",xmin=beginpeak2,xmax=endpeak2,ymin=-Inf,ymax=Inf,fill="orange",alpha=0.2)+
		annotate("rect",xmin=beginpeak3,xmax=endpeak3,ymin=-Inf,ymax=Inf,fill="orange",alpha=0.2))
}
dev.off()
}
