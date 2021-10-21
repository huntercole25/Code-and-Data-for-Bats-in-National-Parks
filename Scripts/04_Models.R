library(glmmTMB)
library(data.table)
library(lubridate)
library(plyr)
library(ggplot2)
library(performance)
library(DHARMa)
library(EcoCountHelper)

# Reading prepped data
SppAccpData <- fread("./Data/Benchmark_Data/03_SppAccpFinalParams.csv") 

# Converting SampleDate vector from character to Date class
SppAccpData[,SampleDate := ymd(SampleDate, tz = "America/Denver")]

# Making site and year factor class for categorical interpretation by models
SppAccpData[,Site := as.factor(Site)]
SppAccpData[,Year := as.factor(Year)]

# Creating list of species classified by SonoBat
UntruncSpecies <- grep("^[[:alpha:]]{4}$", names(SppAccpData), value = T)
PossibleSpecies <- UntruncSpecies[!UntruncSpecies == "Site" & !UntruncSpecies == "Elev" &
                                    !UntruncSpecies == "Yday" & !UntruncSpecies == "Year"]

# Truncating data to species that are present in 50 or more site-nights of data
ObsList <- list()

for(i in PossibleSpecies){
  TmpTable <- data.table(Species = i, PresentNights = nrow(SppAccpData[SppAccpData[[i]]>0,]))
  ObsList <- c(ObsList, list(TmpTable))
}

ObsTable <- rbindlist(ObsList)
TruncSpecies <- ObsTable[PresentNights >= 50, Species]

# Creating function to generate candidate models for each species
IdealModeller <- function(Data, Species){
  x <- get("x")
  
  BatData <- get(deparse(substitute(Data)))
  
  SpNb1 <- tryCatch({glmmTMB(BatData[[Species]] ~ StrScale + CoolScale + BrightScale + Year + YdayScale + ElevScale +
                     WaterScale + DevelScale + ForestScale + MoonScale + (1|Site),
                   data = BatData, family = nbinom1(link = "log"))},
                   error = function(cond){return(NA)},
                   warning = function(cond){return(NA)})
  x <- x+1
  setTxtProgressBar(pb, x)
  
  assign(paste0(Species, "Nb1"), SpNb1, pos = .GlobalEnv)
  
  SpNb2 <- tryCatch({glmmTMB(BatData[[Species]] ~ StrScale + CoolScale + BrightScale + Year + YdayScale + ElevScale +
                     WaterScale + DevelScale + ForestScale + MoonScale + (1|Site),
                   data = BatData, family = nbinom2(link = "log"))},
                   error = function(cond){return(NA)},
                   warning = function(cond){return(NA)})
  x <- x+1
  setTxtProgressBar(pb, x)
  
  assign(paste0(Species, "Nb2"), SpNb2, pos = .GlobalEnv)
  
  SpPoi <- tryCatch({glmmTMB(BatData[[Species]] ~ StrScale + CoolScale + BrightScale + Year + YdayScale + ElevScale +
                     WaterScale + DevelScale + ForestScale + MoonScale + (1|Site),
                   data = BatData, family = poisson(link = "log"))},
                   error = function(cond){return(NA)},
                   warning = function(cond){return(NA)})
  x <- x+1
  setTxtProgressBar(pb, x)
  
  assign(paste0(Species, "Poi"), SpPoi, pos = .GlobalEnv)
  
  SpZiNb1 <- tryCatch({glmmTMB(BatData[[Species]] ~ StrScale + CoolScale + BrightScale + Year + YdayScale + ElevScale +
                       WaterScale + DevelScale + ForestScale + MoonScale + (1|Site),
                   data = BatData, family = nbinom1(link = "log"), ziformula = ~YdayScale + Site)},
                   error = function(cond){return(NA)},
                   warning = function(cond){return(NA)})
  x <- x+1
  setTxtProgressBar(pb, x)
  
  assign(paste0(Species, "ZiNb1"), SpZiNb1, pos = .GlobalEnv)

  SpZiNb2 <- tryCatch({glmmTMB(BatData[[Species]] ~ StrScale + CoolScale + BrightScale + Year + YdayScale + ElevScale +
                       WaterScale + DevelScale + ForestScale + MoonScale + (1|Site),
                     data = BatData, family = nbinom2(link = "log"), ziformula = ~YdayScale + Site)},
                     error = function(cond){return(NA)},
                     warning = function(cond){return(NA)})
  x <- x+1
  setTxtProgressBar(pb, x)
  
  assign(paste0(Species, "ZiNb2"), SpZiNb2, pos = .GlobalEnv)
  
  SpZiPoi <- tryCatch({glmmTMB(BatData[[Species]] ~ StrScale + CoolScale + BrightScale + Year + YdayScale + ElevScale +
                       WaterScale + DevelScale + ForestScale + MoonScale + (1|Site),
                     data = BatData, family = poisson(link = "log"), ziformula = ~YdayScale + Site)},
                     error = function(cond){return(NA)},
                     warning = function(cond){return(NA)})
  x <- x+1
  setTxtProgressBar(pb, x)
  x <<- x
  
  assign(paste0(Species, "ZiPoi"), SpZiPoi, pos = .GlobalEnv)
}

# Starts progress bar
pb <- txtProgressBar(0, (6*length(TruncSpecies)), 0, style = 3)
x <- 0

lapply(TruncSpecies, IdealModeller, Data = SppAccpData)
close(pb)

# SAving glmmTMB objects to file
AllGlmmTMB <- Filter(function(x) inherits(get(x), "glmmTMB"), ls(pos = .GlobalEnv))

ModSaver <- function(ModList, SaveDir){
  ModSaverSub <- function(Model){
    saveRDS(get(Model, pos = .GlobalEnv), file = paste0(SaveDir, "/", Model, ".rds"))
  }
  lapply(ModList, FUN = ModSaverSub)
}

if(!"Models" %in% dir()){dir.create("./Models")}

ModSaver(AllGlmmTMB, "./Models")

# Creating AIC tables for each species
ModelCompare(TruncSpecies, TopModAll)

# Creating plots to visually assess mean-variance relationship for each species to inform
# error-distribution family selection
DistFitWide(c("Year", "Site"), SppAccpData, TruncSpecies)

# Checking VIFs for all top models
for(i in TopModAll$TopModel){
  TmpMod <- get(i)
  print(i)
  print(check_collinearity(TmpMod))
}

# Creating simulated residual diagnostic plots for checking goodness-of-fit
ResidPlotWide(SppAccpData, TopModAll$TopModel, "^[[:alpha:]]{4}", TestVals = F)

if(!"TopMods" %in% dir("./Models")){dir.create("./Models/TopMods")}

# Saving selected models for each species
for(i in TopModAll$TopModel){
  saveRDS(get(i), file = paste0("./Models/TopMods/", i, ".rds"))
}
