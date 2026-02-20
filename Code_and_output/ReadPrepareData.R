library(arrow)
library(dplyr)

# dataset files fulldataset.csv, masks20.parquet, masks21.parquet and owd.parquet should be stored in the parent directory of the code

# reads dataset from Tausk & Spira 2025
data_bsdt = read.csv(file="..\\fulldataset.csv",fileEncoding="UTF-8",stringsAsFactors=FALSE)
data_bsdt = data_bsdt[order(data_bsdt$location),]
countries = data_bsdt$location

# adds to the dataset the variable cardlife, which is the principal component of cardiovascular death rate and life expectancy
prcomp_cardlife = prcomp(data_bsdt[,c("cardiovasc_death_rate","life_expectancy")],scale=TRUE)
prcomp1_cardlife = prcomp_cardlife$x[,1] # extracts the principal component
data_bsdt$cardlife = prcomp1_cardlife

# reads datasets downloaded from https://github.com/csthiago/reanalysis_mask
masks20 = read_parquet("..\\masks20.parquet")
masks21 = read_parquet("..\\masks21.parquet")
owd = read_parquet("..\\owd.parquet")

masks20 = masks20[masks20$location_name %in% countries,c("location_name","date","mask_use_mean")]
masks21 = masks21[masks21$location_name %in% countries,c("location_name","date","mask_use_mean")]
owd = owd[owd$location %in% countries,c("location","date","new_deaths_smoothed_per_million")]

masks = rbind(masks20,masks21)
masks$date = as.Date(masks$date,format="%Y-%m-%d")
owd$date = as.Date(owd$date,format="%Y-%m-%d")

dataset = left_join(masks,owd,by=c("location_name"="location","date"="date"))
dataset = dataset[order(dataset$location_name,dataset$date),]
dataset = rename(dataset,deathcovid = new_deaths_smoothed_per_million)
