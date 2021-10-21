library(data.table)
library(lubridate)
library(ape)
library(FNN)
library(geoR)
library(EcoCountHelper)

####SppAccp Finalization####
SppAccp <- fread("./Data/Benchmark_Data/01_SppAccp1617.csv")
SpatialParams <- fread("./Data/Benchmark_Data/02_SiteSpatialParams.csv")

# Renaming non-intermediate columns and removing intermediate columns
NewNames <- c("Site", "LightCount", "BrightCount", "WaterDist", "Elev", "EastingUsAea", "NorthingUsAea", "PropCool", "ManualDevelPct", "ManualForestPct", "StrWeight")

setnames(x = SpatialParams, old = c("SiteCode", "LightCount", "BrightCount", "WaterDist", "Elev", "Lon.1", "Lat.1", "PropCool", "ManualDevelPct", "ManualForestPct", "StrWeight"),
         new = NewNames, skip_absent = T)

SpatialParams <- SpatialParams[,..NewNames]

# Filling NAs with zeros
for(i in 2:length(SpatialParams)){
  SpatialParams[[i]] <- nafill(SpatialParams[[i]], fill = 0)
}

# Merging bat detection data with site-level spatial parameters
SppAccpFinal <- merge(SppAccp, SpatialParams, by = "Site")

# Creating ordinal date vector
SppAccpFinal[,Yday := yday(SampleDate)]

# Creating year vector
SppAccpFinal[,Year := year(SampleDate)]

# Creating vectors with values for all numeric vectors standardized by two standard deviations
SppAccpFinal[,MoonScale := scale2(MoonPct)]
SppAccpFinal[,ForestScale := scale2(ManualForestPct)]
SppAccpFinal[,DevelScale := scale2(ManualDevelPct)]
SppAccpFinal[,WaterScale := scale2(WaterDist)]
SppAccpFinal[,ElevScale := scale2(Elev)]
SppAccpFinal[,YdayScale := scale2(Yday)]
SppAccpFinal[,BrightScale := scale2(BrightCount)]
SppAccpFinal[,CoolScale := scale2(PropCool)]
SppAccpFinal[,StrScale := scale2(StrWeight)]

# Writing prepped data to file
write.csv(SppAccpFinal, file = "./Data/Benchmark_Data/03_SppAccpFinalParams.csv", row.names = F)
