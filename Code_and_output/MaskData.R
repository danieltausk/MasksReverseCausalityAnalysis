source("FindWaves.R")

library(boot)

# computes indices maskbeginwave and maskpeakwave, one per wave (3 x 24 = 72 values)

mask_beginwave = numeric(0)
mask_peakwave = numeric(0)
for(i in 1:nrow(dataperwave))
{
	dataset_countrywave = datasetwaves[datasetwaves$location_name == dataperwave$location_name[i] & datasetwaves$wave == dataperwave$wave[i],]
	mask_beginwave = c(mask_beginwave,mean(extract_beginning(dataset_countrywave$mask_use_mean)))
	mask_peakwave = c(mask_peakwave,mean(extract_peak(dataset_countrywave$mask_use_mean,dataset_countrywave$deathcovid)))
}
dataperwave$mask_beginwave = mask_beginwave
dataperwave$mask_peakwave = mask_peakwave

write.csv(file="dataperwave.csv",x=dataperwave,row.names=FALSE)

# computes all five mask indices (one per country)

maskdata = dataset |> group_by(location_name) |> summarize(mask_all=mean(mask_use_mean)) |> as.data.frame()
maskdata = datasetwaves[datasetwaves$wave == 0,] |> group_by(location_name) |> summarize(mask_interwave=mean(mask_use_mean)) |> left_join(maskdata,y=_,by="location_name")
maskdata = datasetwaves[datasetwaves$wave != 0,] |> group_by(location_name) |> summarize(mask_inwave=mean(mask_use_mean)) |> left_join(maskdata,y=_,by="location_name")
mask_beginwave = numeric(0)
mask_peakwave = numeric(0)
for(country in maskdata$location_name)
{
	dataset_countrywave1 = datasetwaves[datasetwaves$location_name == country & datasetwaves$wave == 1,]
	dataset_countrywave2 = datasetwaves[datasetwaves$location_name == country & datasetwaves$wave == 2,]
	dataset_countrywave3 = datasetwaves[datasetwaves$location_name == country & datasetwaves$wave == 3,]
	mask_beginwave = c(mask_beginwave,mean(c(extract_beginning(dataset_countrywave1$mask_use_mean),
		extract_beginning(dataset_countrywave2$mask_use_mean),
		extract_beginning(dataset_countrywave3$mask_use_mean))))
	mask_peakwave = c(mask_peakwave,mean(c(extract_peak(dataset_countrywave1$mask_use_mean,dataset_countrywave1$deathcovid),
		extract_peak(dataset_countrywave2$mask_use_mean,dataset_countrywave2$deathcovid),
		extract_peak(dataset_countrywave3$mask_use_mean,dataset_countrywave3$deathcovid))))
}
maskdata$mask_beginwave = mask_beginwave
maskdata$mask_peakwave = mask_peakwave

write.csv(file="maskdata.csv",x=maskdata,row.names=FALSE)

# scatter plot of maskpeakwave against maskbeginwave, one point per wave
ggplot(data=dataperwave)+geom_point(aes(x=mask_beginwave,y=mask_peakwave,col=factor(wave)))+
	scale_x_continuous(lim=c(0,1),name="maskbeginwave")+scale_y_continuous(lim=c(0,1),name="maskpeakwave")+
	scale_color_manual(values=c("1"="black","2"="red","3"="blue"),labels=c("1"="first","2"="second","3"="third"),name="Wave",guide=guide_legend(order=1))+
	geom_smooth(aes(x=mask_beginwave,y=mask_peakwave,linetype="regression line"),method="lm",se=FALSE,formula=y~x,col="orange")+
	scale_linetype_manual(values="solid",name=NULL,guide=guide_legend(order=2))
ggsave(filename="mask_begin_peak_wave.jpg",device="jpeg",width=14.23,height=6.77)

dataperwave2 = mutate(dataperwave,quotient=mask_peakwave/mask_beginwave)

# scatter plot of quotient maskpeakwave/maskbeginwave against mean Covid-19 mortality per million during the wave
# one point per wave
ggplot(data=dataperwave2)+geom_point(aes(x=mean_death,y=quotient,col=factor(wave)))+
	labs(x="mean Covid-19 deaths per million during wave",y="maskpeakwave/maskbeginwave")+
	scale_y_continuous(lim=c(0,10),breaks=0:10)+
	scale_color_manual(values=c("1"="black","2"="red","3"="blue"),labels=c("1"="first","2"="second","3"="third"),name="Wave",guide=guide_legend(order=1))+
	geom_smooth(aes(x=mean_death,y=quotient,linetype="loess"),method="loess",se=FALSE,formula=y~x,col="orange")+
	scale_linetype_manual(values="solid",name=NULL,guide=guide_legend(order=2))
ggsave(filename="quotient_versus_mortality.jpg",device="jpeg",width=14.23,height=6.77)

maskdata2 = dataset |> group_by(location_name) |> summarize(total_death=sum(deathcovid)) |> left_join(maskdata,y=_,by="location_name")
maskdata2 = mutate(maskdata2,quotient=mask_inwave/mask_interwave)

# scatter plot of quotient maskinwave/maskinterwave against total Covid-19 mortality per million
# one point per country
ggplot(data=maskdata2)+geom_point(aes(x=total_death,y=quotient))+
	labs(x="total Covid-19 deaths per million",y="maskinwave/maskinterwave")+
	scale_y_continuous(breaks=0:7,lim=c(0,7))+
	geom_smooth(aes(x=total_death,y=quotient,linetype="loess"),method="loess",se=FALSE,formula=y~x,col="orange")+
	scale_linetype_manual(values="solid",name=NULL)
ggsave(filename="quotient_versus_mortality2.jpg",device="jpeg",width=14.23,height=6.77)

set.seed(123) # random seed for boostrap confidence intervals

cor_bootstrap_CI = function(x,y,method,conf=0.95,R=10^4){
# bootstrap confidence interval for correlation between x and y

	boots = boot(cbind(x,y),statistic=function(d,ind) cor(d[ind,1],d[ind,2],method=method),R=R)
	CI = boot.ci(boots,type="basic",conf=conf)$basic[4:5]
	return(list(est=cor(x,y,method=method),left=max(min(CI[1],1),-1),right=max(min(CI[2],1),-1)))
}

sink("VariousCorrelations.txt")

cat("maskbeginwave versus maskpeakwave (one datapoint per wave):\n")
CI = cor_bootstrap_CI(dataperwave$mask_beginwave,dataperwave$mask_peakwave,method="pearson")
cat("pearson correlation: ",round(CI$est,2)," 95% CI: [",round(CI$left,2),",",round(CI$right,2),"]\n")

cat("\n")

cat("maskpeakwave/maskbeginwave versus mean mortality (one datapoint per wave):\n")
CI = cor_bootstrap_CI(dataperwave2$quotient,dataperwave2$mean_death,method="spearman")
cat("spearman correlation: ",round(CI$est,2)," 95% CI: [",round(CI$left,2),",",round(CI$right,2),"]\n")

cat("\n")

cat("maskinwave/maskinterwave versus total Covid mortality (one datapoint per country):\n")
CI = cor_bootstrap_CI(maskdata2$quotient,maskdata2$total_death,method="spearman")
cat("spearman correlation: ",round(CI$est,2)," 95% CI: [",round(CI$left,2),",",round(CI$right,2),"]\n")

cat("\n")

cat("pearson correlations between mask indices (one datapoint per country):\n")

for(i in 2:(ncol(maskdata)-1))
{
	for(j in (i+1):ncol(maskdata))
	{
		CI = cor_bootstrap_CI(maskdata[,i],maskdata[,j],method="pearson")
		cat(names(maskdata)[i]," versus ",names(maskdata)[j],": ",round(CI$est,2)," 95% CI: [",round(CI$left,2),",",round(CI$right,2),"]\n")
	}
}

sink()

sink("MaskingComparisons.txt")

cat("Paired Wilcoxon test, H1: maskpeakwave > maskbeginwave (one datapoint per wave):\n")
print(wilcox.test(dataperwave$mask_peakwave,dataperwave$mask_beginwave,paired=TRUE,alternative="greater"))

cat("\n")

cat("maskpeakwave/maskbeginwave (one datapoint per wave):\n")
cat("median: ",round(median(dataperwave2$quotient),2),"\n")
cat("first quartile: ",round(quantile(dataperwave2$quotient,1/4),2),"\n")
cat("third quartile: ",round(quantile(dataperwave2$quotient,3/4),2),"\n")
cat("minimum: ",round(min(dataperwave2$quotient),2),"\n")
cat("maximum: ",round(max(dataperwave2$quotient),2),"\n")

cat("\n")

cat("Paired Wilcoxon test, H1: maskinwave > maskinterwave (one datapoint per country):\n")
print(wilcox.test(maskdata$mask_inwave,maskdata$mask_interwave,paired=TRUE,alternative="greater"))

cat("\n")

cat("maskinwave/maskinterwave (one datapoint per country):\n")
cat("median: ",round(median(maskdata2$quotient),2),"\n")
cat("first quartile: ",round(quantile(maskdata2$quotient,1/4),2),"\n")
cat("third quartile: ",round(quantile(maskdata2$quotient,3/4),2),"\n")
cat("minimum: ",round(min(maskdata2$quotient),2),"\n")
cat("maximum: ",round(max(maskdata2$quotient),2),"\n")

sink()
