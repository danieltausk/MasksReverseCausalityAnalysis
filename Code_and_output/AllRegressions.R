source("MaskData.R")

library(betaDelta)

output_reg = function(reg,level=0.95){
# generates output for linear regression stored in reg
# standardized coefficients, confidence intervals and p-values

	bd = BetaDelta(reg,type="mvn")
	ci = confint(bd,level=level)

	# uses p-values for non standardized coefficients --- null hypothesis is the same, but p-values are more accurate
	# because they are exact at least under normality assumptions
	output = data.frame(coefnames=names(bd$est),stdcoef=bd$est,
			CIleft=ci[,1],CIright=ci[,2],pval=summary(reg)$coefficients[-1,4])
	row.names(output) = NULL

	return(output)
}

# relevant dataset for regressions
data_reg = data_bsdt[,c("location","people_fully_vaccinated_per_hundred","human_development_index","cardlife","Per.Levitt.Age.Adjusted")]
data_reg = left_join(maskdata,data_reg,by=c("location_name"="location"))

reg_mask_all = lm(Per.Levitt.Age.Adjusted~mask_all+people_fully_vaccinated_per_hundred+human_development_index+cardlife,data=data_reg)
reg_mask_interwave = lm(Per.Levitt.Age.Adjusted~mask_interwave+people_fully_vaccinated_per_hundred+human_development_index+cardlife,data=data_reg)
reg_mask_inwave = lm(Per.Levitt.Age.Adjusted~mask_inwave+people_fully_vaccinated_per_hundred+human_development_index+cardlife,data=data_reg)
reg_mask_beginwave = lm(Per.Levitt.Age.Adjusted~mask_beginwave+people_fully_vaccinated_per_hundred+human_development_index+cardlife,data=data_reg)
reg_mask_peakwave = lm(Per.Levitt.Age.Adjusted~mask_peakwave+people_fully_vaccinated_per_hundred+human_development_index+cardlife,data=data_reg)

oldoptions = options()
sink("MainRegressions.txt")
options(width = 1000)

print(output_reg(reg_mask_all),row.names=FALSE)
cat("\n")
print(output_reg(reg_mask_interwave),row.names=FALSE)
cat("\n")
print(output_reg(reg_mask_inwave),row.names=FALSE)
cat("\n")
print(output_reg(reg_mask_beginwave),row.names=FALSE)
cat("\n")
print(output_reg(reg_mask_peakwave),row.names=FALSE)

options(oldoptions)
sink()

oldoptions = options()
sink("OtherRegressions.txt")
options(width = 1000)

cat("Regressing maskpeakwave on maskbeginwave (one datapoint per wave):\n\n")
reg = lm(mask_peakwave~mask_beginwave,data=dataperwave)
print(summary(reg))

cat("\n")

print(confint(reg))

cat("\n")

cat("Regressing maskinwave on maskinterwave (one datapoint per country):\n\n")
reg = lm(mask_inwave~mask_interwave,data=maskdata)
print(summary(reg))

cat("\n")

print(confint(reg))

cat("\n")

cat("Regressing maskpeakwave on maskbeginwave and mean mortality (one datapoint per wave):\n\n")
reg = lm(mask_peakwave~mask_beginwave+mean_death,data=dataperwave)
print(summary(reg))
cat("\n")
print(output_reg(reg),row.names=FALSE)

cat("\n")

cat("Regressing maskinwave on maskinterwave and total Covid deaths (one datapoint per country):\n\n")
reg = lm(mask_inwave~mask_interwave+total_death,data=maskdata2)
print(summary(reg))
cat("\n")
print(output_reg(reg),row.names=FALSE)

options(oldoptions)
sink()
