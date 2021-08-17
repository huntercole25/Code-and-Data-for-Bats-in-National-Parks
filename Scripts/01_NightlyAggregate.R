library(data.table)
library(lubridate)
library(tidyr)
library(lunar)

# Reading individual call sequence classifications table
SbFull1617 <- fread("./Data/RawData/RawCalls.csv", fill = T)

# Appending temporal metadata to each call sequence
SbFull1617[,DetectDt := ymd_hms(regmatches(Filename, regexpr("\\d{8}_\\d{6}", Filename)),
                                tz = "America/Denver")]

SbFull1617[,DetectDate := date(DetectDt)]
SbFull1617[hour(DetectDt) >= 0 & hour(DetectDt) <= 18, SampleDate := DetectDate - 1]
SbFull1617[hour(DetectDt) > 18, SampleDate := DetectDate]
SbFull1617[,SampleYear := year(SampleDate)]
SbFull1617[,Site := regmatches(Filename, regexpr("^[[:alnum:]]+", Filename))]

# Loading site coordinates
SiteLocs <- fread("./Data/RawData/SiteLocations.csv")

# # Consider removing - check script data dependencies
SiteList <- unique(SiteLocs$SiteCode)
# SbSiteList <- unique(SbFull1617$Site)
# SiteYears <- data.table(Site = SbSiteList)

#Aggregating positive species detections by site-night
SppAccp1617 <- aggregate(DetectDt ~ Site + SampleDate + SppAccp,
                         SbFull1617[!SppAccp == "",], FUN = length)
setnames(SppAccp1617, old = "DetectDt", new = "Count")

#Making aggregated data wide by species
SppAccp1617 <- as.data.table(spread(SppAccp1617, key = SppAccp, value = Count, fill = 0))

#Assigning spatial coordinates to each site-night
CoordAttr <- function(DataTable, SiteID){
  TmpData <- get("DataTable")
  SiteLocs <- get("SiteLocs")
  TmpData[Site == SiteID, Lat := SiteLocs$Lat[SiteLocs$SiteCode==SiteID]]
  TmpData[Site == SiteID, Lon := SiteLocs$Lon[SiteLocs$SiteCode==SiteID]]
  assign("DataTable", TmpData, pos = .GlobalEnv)
}

lapply(SiteList, CoordAttr, DataTable = SppAccp1617)

#Adding moon illumination data to each site-night
SppAccp1617[,MoonPct := lunar.illumination(SampleDate)]

#Write aggregated data to file
write.csv(x = SppAccp1617, file = "./Data/Benchmark_Data/01_SppAccp1617.csv", row.names = F)