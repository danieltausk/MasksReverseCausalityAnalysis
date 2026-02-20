source("ReadPrepareData.R")

library(ggplot2)

# normalizes covid deaths in each country by dividing by total covid deaths
dataset2 = dataset |> group_by(location_name) |>
		mutate(deathnorm = deathcovid/sum(deathcovid)) |> ungroup()

# for each date, sums normalized covid deaths over all countries
dataset2total = dataset2 |>
		group_by(date) |> summarize(totaldeathnorm = sum(deathnorm))

# arbitrary time intervals in which the global minimum of totaldeathnorm coincides with the relevant local minima
containmin1 = dataset2total[dataset2total$date > as.Date("2020-05-01","%Y-%m-%d") & dataset2total$date < as.Date("2020-12-01","%Y-%m-%d"),]
containmin2 = dataset2total[dataset2total$date > as.Date("2020-12-01","%Y-%m-%d"),]

# end dates of pandemic phases (phases are time periods that contain the waves)
endphase1 = containmin1$date[which.min(containmin1$totaldeathnorm)]
endphase2 = containmin2$date[which.min(containmin2$totaldeathnorm)]

# xaxis breaks and labels for graphs
xaxis_br = c(sapply(2:12,function(m) paste0("2020-",m,"-01")),sapply(1:12,function(m) paste0("2021-",m,"-01")))
xaxis_br = as.Date(xaxis_br,format="%Y-%m-%d")
xaxis_lab = c(sapply(2:12,function(m) paste0("20/",sprintf("%02d",m))),sapply(1:12,function(m) paste0("21/",sprintf("%02d",m))))

ggplot(data=dataset2,aes(x=date,y=deathnorm,col=location_name))+geom_line()+
	scale_x_continuous(breaks=xaxis_br,labels=xaxis_lab)+
	labs(y="normalized Covid-19 deaths",col="Country")+
	geom_vline(xintercept=c(endphase1,endphase2))
ggsave(filename="Phases.jpg",device="jpeg",width=14.23,height=6.77)

sink("Phases.txt")
cat("Phase 1: from ",format(min(dataset$date),"%Y/%m/%d")," to ",format(endphase1,"%Y/%m/%d"),"\n")
cat("Phase 2: from ",format(endphase1+1,"%Y/%m/%d")," to ",format(endphase2,"%Y/%m/%d"),"\n")
cat("Phase 3: from ",format(endphase2+1,"%Y/%m/%d")," to ",format(max(dataset$date),"%Y/%m/%d"),"\n")
sink()
