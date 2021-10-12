library(data.table)
library(rgdal)
library(raster)
library(gdalUtils)
library(rgeos)
library(sf)
library(tidyr)
library(ape)

# Load detector locations
SiteLocations <- fread("./Data/RawData/SiteLocations.csv")

# Create spatial object with site names and coordinates
SitePoints <- SpatialPointsDataFrame(coords = SiteLocations[,c("Lon", "Lat")], data = SiteLocations,
                                     proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Write shapefile with native projection 
writeOGR(SitePoints, dsn = "./Data/GisData/SiteLocations", layer = "SitesWGS84", driver = "ESRI Shapefile", overwrite_layer = T)

# Project site points to USGS Albers Equal Area projection
SitePointsProj <- spTransform(SitePoints, CRSobj = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 
                                                       +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Write shapefile with USGS projection
writeOGR(SitePointsProj, dsn = "./Data/GisData/SiteLocations", layer = "SitesUsAea", driver = "ESRI Shapefile", overwrite_layer = T)

# Load GRTE Boundaries 
GrteBoundsRaw <- readOGR(dsn = "./Data/GisData/GrandTetonBoundaries", layer = "GrteBoundariesWGS84")

# Project GRTE boundaries to USGS Albers Equal Area
GrteBoundsUsAea <- spTransform(GrteBoundsRaw, CRSobj = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 
                                                       +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Project GRTE boundaries to NAD83 for later raster clipping
GrteBoundsNad83 <- spTransform(GrteBoundsRaw, CRSobj = crs("+proj=longlat +datum=NAD83 +no_defs"))

# Write shapefile with projected GRTE boundaries
writeOGR(GrteBoundsUsAea, dsn = "./Data/GisData/GrandTetonBoundaries", layer = "GrteBoundariesUsAea",
         driver = "ESRI Shapefile", overwrite_layer = T)

# Read in lakes and flowline shapefiles
RiversRaw <- readOGR(dsn = "./Data/GisData/Watershed/Shape", layer = "NHDFlowline")
LakesRaw <- readOGR(dsn = "./Data/GisData/Watershed/Shape", layer = "NHDWaterbody")

# Keep flowline Fcodes 46006 and 55800
RiversTrunc <- RiversRaw[RiversRaw$FCode == 46006 | RiversRaw$FCode == 55800,]

# Project shapefiles to USGS Albers Equal Area
RiversUsAea <- spTransform(RiversTrunc, CRSobj = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 
                                      +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
LakesUsAea <- spTransform(LakesRaw, CRSobj = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 
                                      +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Write river and lake shapefiles with projected geometries for future use
writeOGR(obj = RiversUsAea, dsn = "./Data/GisData/Watershed/Shape", layer = "RiversTruncUsAea", driver = "ESRI Shapefile", overwrite_layer = T)
writeOGR(obj = LakesUsAea, dsn = "./Data/GisData/Watershed/Shape", layer = "LakesUsAea", driver = "ESRI Shapefile", overwrite_layer = T)

# Use gDistance in rgeos to calculate distance from sites to body of water polygons and
# lines (NHD river and lake shapefiles)
SitePointsProj$RiverDist <- apply(gDistance(SitePointsProj, RiversUsAea, byid = T), 2, min)

#Calculate distance to water bodies
SitePointsProj$LakeDist <- apply(gDistance(SitePointsProj, LakesUsAea, byid = T), 2, min)
  
#Assign minimum distance vector to sites
SitePointsProj$WaterDist[SitePointsProj$LakeDist < SitePointsProj$RiverDist] <-
  SitePointsProj$LakeDist[SitePointsProj$LakeDist < SitePointsProj$RiverDist]

SitePointsProj$WaterDist[SitePointsProj$LakeDist > SitePointsProj$RiverDist] <-
  SitePointsProj$RiverDist[SitePointsProj$LakeDist > SitePointsProj$RiverDist]

# Write site shapefile with updated spatial metadata
writeOGR(SitePointsProj, dsn = "./Data/GisData/SiteLocations", layer = "SitesUsAea", driver = "ESRI Shapefile", overwrite_layer = T)

# Read in WGS84 light location points
LightSheet <- fread("./Data/RawData/RawLightLocs.csv")

# Make decimal degree longitude values negative to reflect apropriate hemisphere
LightSheet[,X := -X]

# Truncating vector names
setnames(LightSheet, old = c("Light ID", "Y", "X", "Position", "Brightness Index Score", "Color", "Height", "Directions Projected",
                             "Distance Cast", "Building #", "Switch Activated?", "Other Info"),
         new = c("LightID", "Lat", "Lon", "RelPosition", "BrightScore", "Color", "Height", "DirProject", "DistCast", "StructID",
                 "SwitchAct", "Notes"))

# Making light color category names consistent
LightSheet[Color == "Yellow white", Color := "Yellow-white"]
LightSheet[Color == "Blue white", Color := "Blue-white"]

# Create spatial object with light locations and light characteristic metadata
LightPointsWGS84 <- SpatialPointsDataFrame(coords = LightSheet[,c("Lon", "Lat")], data = LightSheet,
                       proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Project light locations spatial object to USGS Albers Equal Area projection
LightPointsUsAea <- spTransform(LightPointsWGS84, CRSobj = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96
                                      +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Write shapefile with projected light locations and light characteristic metadata
writeOGR(LightPointsUsAea, dsn = "./Data/GisData/Lights", layer = "LightsUsAea", driver = "ESRI Shapefile", overwrite_layer = T)

# Consider deleting this line - check data structures and vector names of Sp object vs shapefile
LightPointsUsAea <- readOGR(dsn = "./Data/GisData/Lights", layer = "LightsUsAea")

# Create 50m and 500m buffers around sites
Site50buff <- buffer(SitePointsProj, width = 50, dissolve = F)
Site500buff <- buffer(SitePointsProj, width = 500, dissolve = F)

# Write shapefiles with 50 meter and 500 meter buffers surrounding each site
writeOGR(Site50buff, dsn = "./Data/GisData/SiteBuffers", layer = "SiteBuf50UsAea", driver = "ESRI Shapefile", overwrite_layer = T)
writeOGR(Site500buff, dsn = "./Data/GisData/SiteBuffers", layer = "SiteBuf500UsAea", driver = "ESRI Shapefile", overwrite_layer = T)

writeOGR(SitePointsProj, dsn = "./Data/GisData/SiteLocations", layer = "SitesUsAea", driver = "ESRI Shapefile",
         overwrite_layer = T)

####Elevation####
# Compiles DEM tiles into single raster
TileList <- list.files("./Data/GisData/Tiles", full.names = T)

AllTiles <- lapply(list.files("./Data/GisData/Tiles", full.names = T), function(x)raster(x))
GrteRawDem <- do.call(merge, AllTiles)

#Writes raw full merged DEM
writeRaster(GrteRawDem, "./Data/GisData/Elevation/FullRawDem.tif", format = "GTiff", overwrite = T)

# Projects clipped GRTE raster to Alber's Equal Area Conic
GrteDemUsAea <- projectRaster(GrteDemRaw, crs = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96
                                      +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Writes clipped raster to file for later use
writeRaster(GrteDemUsAea, filename = "./Data/GisData/Elevation/GrteDemUsAea.tif", format = "GTiff", overwrite = T)

# GrteDemUsAea <- raster("./Data/GisData/Elevation/GrteDemUsAea.tif")

SitePointsProj$Elev <- raster::extract(GrteDemUsAea, SitePointsProj) # Extracts elevation for each site

####Bright Sum####
# Creates table of all lights with 50 meters of each site
LightSiteTable <- as.data.table(intersect(LightPointsUsAea, Site50buff))

# Calculates the sum of brightness index scores within 50 meters for each site
LightSumSite <- aggregate(LightSiteTable$BrghtSc ~ LightSiteTable$SiteCode, FUN = sum)

# Calculates the number of lights within 50 meters of each site
LightCountSite <- aggregate(LightSiteTable$BrghtSc ~ LightSiteTable$SiteCode, FUN = length)

# Renaming table columns
setnames(LightSumSite, c("SiteCode", "BrightCount"))
setnames(LightCountSite, c("SiteCode", "LightCount"))

# Assigning values calculated above to spatial object
SitePointsProj <- merge(LightSumSite, SitePointsProj, by = "SiteCode", all.y = T)
SitePointsProj <- merge(LightCountSite, SitePointsProj, by = "SiteCode", all.y = T)
SitePointsProj <- SpatialPointsDataFrame(SitePointsProj[,c("Lon.1", "Lat.1")], SitePointsProj)

# Writes updated site shapefile
writeOGR(SitePointsProj, dsn = "./Data/GisData/SiteLocations", layer = "SitesUsAea",
         driver = "ESRI Shapefile", overwrite_layer = T)

####Binary Color Temp####
# Assigning warm and cool color temperature classifications to qualitative color classifications
LightPointsUsAea$ColorTemp[LightPointsUsAea$Color == "Orange" | LightPointsUsAea$Color == "Yellow-white" | LightPointsUsAea$Color == "Yellow"] <- "Warm"
LightPointsUsAea$ColorTemp[LightPointsUsAea$Color == "Blue-white" | LightPointsUsAea$Color == "Green-white" | LightPointsUsAea$Color == "White"] <- "Cool"

# Making data table of light color classifications within 50 meters of sites
LightColTable <- as.data.table(intersect(LightPointsUsAea, Site50buff))
LightColTable[, SiteCode := as.character(SiteCode)]

# The function below tabulates the number and proportion of cool-colored lights within 50 meters of each site
CoolCalculator <- function(SiteLightIntersect, SiteCol, LightTempCol, SiteList, FinalTableName){
  require(data.table)
  DataList <- list()
  DataName <- deparse(substitute(SiteLightIntersect))
  
  CoolCalculatorSub <- function(z){
    FullData <- get(DataName, pos = .GlobalEnv)
    SiteData <- FullData[FullData[[SiteCol]] == z,]
    NTotal <- nrow(SiteData)
    NCool <- nrow(SiteData[SiteData[[LightTempCol]]=="Cool",])
    TmpData <- data.table(SiteCode = z, PropCool = (NCool/NTotal))
    TmpList <- get("DataList")
    DataList <<- c(list(TmpData), TmpList)
  }
  
  lapply(SiteList, CoolCalculatorSub)
  
  FinalData <- rbindlist(DataList)
  assign(deparse(substitute(FinalTableName)), FinalData, pos = .GlobalEnv)
}

SiteCodeList <- unique(LightColTable$SiteCode)
CoolCalculator(LightColTable, "SiteCode", "ColorTemp", SiteCodeList, LightColSite)

# Merging results of CoolCalculator call above with site shapefile
SitePointsProj <- sp::merge(SitePointsProj, LightColSite, by.x = "SiteCode", by.y = "SiteCode", all.x = T)

write.csv(as.data.table(SitePointsProj), file = "./Data/Benchmark_Data/02_SiteSpatialParams.csv", row.names = F)

####Manaual Landcover####
# Reading manually classified landcover vector data and truncating to 50 meters from each site
ManualLandcover <- readOGR(dsn = "./Data/GisData/ManualLandcover", layer = "ManualLandcoverClipped")
ManualLC <- intersect(ManualLandcover, Site50buff)

# Adding an area vector to the truncated manual land cover shapefile
ManualLC$Area <- area(ManualLC)

# Calculating the summed area of each land cover type present within 50 meters of each site
ManualLandAgg <- as.data.table(aggregate(Area ~ LCType + Site, data = ManualLC, FUN = sum))
ManualLandAgg[LCType == 1, LCType := "ManualDevelPct"]
ManualLandAgg[LCType == 2, LCType := "ManualForestPct"]

# Calculating the proportion of all areas within 50 meters of a site covered by each manually classified land cover type
ManualLandAgg$PctCoverage <- ManualLandAgg$Area/(50^2*pi)

# Making the data aggregated above wide-form
ManualLandWide <- as.data.table(pivot_wider(data = ManualLandAgg, id_cols = Site, names_from = LCType, values_from = PctCoverage, values_fill = list(Pctcoverage = 0)))

# Populating zeros for sites without any area covered by a given land cover type
ManualLandWide[is.na(ManualDevelPct), ManualDevelPct := 0]
ManualLandWide[is.na(ManualForestPct), ManualForestPct := 0]

# Merging proportional land cover data with the sites spatial object
SitePointsProj <- merge(as.data.table(SitePointsProj), ManualLandWide, by.x = "SiteCode", by.y = "Site", all.x = T)

write.csv(as.data.table(SitePointsProj), file = "./Data/Benchmark_Data/02_SiteSpatialParams.csv", row.names = F)


####Structures####
Structures <- fread("./Data/RawData/AggStructures.csv")

# Truncating to buildings with coordinates that were classified as suitable for bat inhabitation
Structures <- Structures[is.na(Lat)==F & is.na(Lon)==F & BatSuit == 1, c("Lat", "Lon", "BuildingNum", "BatSuit")]

# Creating a spatial object from the truncated table above
StrRaw <- SpatialPointsDataFrame(coords = Structures[,c("Lon", "Lat")], data = Structures,
                       proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Projecting the data above to albers Equal Area
StrProj <- spTransform(StrRaw, crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96
                                   +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

writeOGR(StrProj, dsn = "./Data/GisData/Structures", layer = "Structures_UsAea",
         driver = "ESRI Shapefile", overwrite_layer = T)

# The function below calculates the reciprocal of summed distances to structures within 1000 meters of a site
DistCounter <- function(x){
  y <- sort(x[x <= 1000])
  z <- 1/y
  sum(z)
}

SitePointsProj <- SpatialPointsDataFrame(coords = SitePointsProj[,c("Lon.1", "Lat.1")],
                                         data = SitePointsProj,
                                         proj4string = crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96
                                   +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

# Sums reciprocal distances to structures within 1000 meters of each site
SitePointsProj$StrWeight <- apply(gDistance(SitePointsProj, StrProj, byid = T), 2, DistCounter)
write.csv(as.data.table(SitePointsProj), file = "./Data/Benchmark_Data/02_SiteSpatialParams.csv", row.names = F)
