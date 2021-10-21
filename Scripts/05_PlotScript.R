library(glmmTMB)
library(ggplot2)
library(data.table)
library(EcoCountHelper)
library(lubridate)
library(cowplot)
library(kableExtra)
library(tidyr)

# Reading data used in models
SppAccpData <- fread("./Data/Benchmark_Data/03_SppAccpFinalParams.csv")

# Prepping data 
SppAccpData[,SampleDate := ymd(SampleDate, tz = "America/Denver")]
SppAccpData[,Site := as.factor(Site)]
SppAccpData[,Year := as.factor(Year)]

# Reading selected model for each species
Mods <- character()
for(i in list.files("./Models/TopMods", full.names = T)){
  TmpMod <- readRDS(i)
  ModName <- regmatches(basename(i), regexpr("^[[:alnum:]]+", basename(i)))
  assign(ModName, TmpMod)
  
  Mods <- c(Mods, ModName)
}

# Specifying labels in correct order for plots
Labs <- c("Intercept (Site)", "Structure Index", "Prop. Cool Lights", "Brightness Index", "Year (16-17)",
          "Ordinal Day", "Elevation", "Water Distance", "Prop. Developed", "Prop. Forest", "Moon Illum.")

# Generating effects plots for each species
EffectsPlotter(Mods, Labs,
               ConfInts = c(95, 85), ThemeBlack = F)

# Adding species title to each effects plot
EfPlots <- ls(pattern = "EffectsPlot")
for(i in EfPlots){
  TmpPlot <- get(i)
  TmpPlot <- TmpPlot + labs(title = substr(i, 1, 4))
  assign(i, TmpPlot)
}

# Creating multi-panel plot with all species effects plots
FullEffects <- DumbGrid(EpfuNb2EffectsPlot, LaciNb2EffectsPlot, LanoNb2EffectsPlot, MyevNb2EffectsPlot,
         MyluNb2EffectsPlot, MyvoNb2EffectsPlot, MyyuNb2EffectsPlot, Ncols = 3, FirstColWidth = 1.5)

# Saving multi-panel effects plot to file
if(!"Output" %in% dir()){dir.create("./Output")}

ggsave("./Output/FullEffects.png", FullEffects, width = 8, height = 8)


# Generating continuous data for each variable spanning the full extent of each
# continuous variable with all other variables held constant. These data will be
# used to visualize continuous model predictions for each numeric parameter.
Forest <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                      ForestScale = seq(min(SppAccpData$ForestScale), max(SppAccpData$ForestScale),
                                        length.out = nrow(SppAccpData)),
                      DevelScale = median(SppAccpData$DevelScale),
                      WaterScale = median(SppAccpData$WaterScale),
                      ElevScale = median(SppAccpData$ElevScale),
                      YdayScale = median(SppAccpData$YdayScale),
                      Year = "2016",
                      BrightScale = median(SppAccpData$BrightScale),
                      CoolScale = median(SppAccpData$CoolScale),
                      StrScale = median(SppAccpData$StrScale), Site = NA)

Devel <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                     ForestScale = median(SppAccpData$ForestScale),
                     DevelScale = seq(min(SppAccpData$DevelScale), max(SppAccpData$DevelScale),
                                      length.out = nrow(SppAccpData)),
                     WaterScale = median(SppAccpData$WaterScale),
                     ElevScale = median(SppAccpData$ElevScale),
                     YdayScale = median(SppAccpData$YdayScale),
                     Year = "2016",
                     BrightScale = median(SppAccpData$BrightScale),
                     CoolScale = median(SppAccpData$CoolScale),
                     StrScale = median(SppAccpData$StrScale), Site = NA)

Water <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                     ForestScale = median(SppAccpData$ForestScale),
                     DevelScale = median(SppAccpData$DevelScale),
                     WaterScale = seq(min(SppAccpData$WaterScale), max(SppAccpData$WaterScale),
                                      length.out = nrow(SppAccpData)),
                     ElevScale = median(SppAccpData$ElevScale),
                     YdayScale = median(SppAccpData$YdayScale),
                     Year = "2016",
                     BrightScale = median(SppAccpData$BrightScale),
                     CoolScale = median(SppAccpData$CoolScale),
                     StrScale = median(SppAccpData$StrScale), Site = NA)

Elev <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                    ForestScale = median(SppAccpData$ForestScale),
                    DevelScale = median(SppAccpData$DevelScale),
                    WaterScale = median(SppAccpData$WaterScale),
                    ElevScale = seq(min(SppAccpData$ElevScale), max(SppAccpData$ElevScale),
                                    length.out = nrow(SppAccpData)),
                    YdayScale = median(SppAccpData$YdayScale),
                    Year = "2016",
                    BrightScale = median(SppAccpData$BrightScale),
                    CoolScale = median(SppAccpData$CoolScale),
                    StrScale = median(SppAccpData$StrScale), Site = NA)

Bright <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                      ForestScale = median(SppAccpData$ForestScale),
                      DevelScale = median(SppAccpData$DevelScale),
                      WaterScale = median(SppAccpData$WaterScale),
                      ElevScale = median(SppAccpData$ElevScale),
                      YdayScale = median(SppAccpData$YdayScale),
                      Year = "2016",
                      BrightScale = seq(min(SppAccpData$BrightScale), max(SppAccpData$BrightScale),
                                        length.out = nrow(SppAccpData)),
                      CoolScale = median(SppAccpData$CoolScale),
                      StrScale = median(SppAccpData$StrScale), Site = NA)

Cool <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                    ForestScale = median(SppAccpData$ForestScale),
                    DevelScale = median(SppAccpData$DevelScale),
                    WaterScale = median(SppAccpData$WaterScale),
                    ElevScale = median(SppAccpData$ElevScale),
                    YdayScale = median(SppAccpData$YdayScale),
                    Year = "2016",
                    BrightScale = median(SppAccpData$BrightScale),
                    CoolScale = seq(min(SppAccpData$CoolScale), max(SppAccpData$CoolScale),
                                    length.out = nrow(SppAccpData)),
                    StrScale = median(SppAccpData$StrScale), Site = NA)

Structures <- expand.grid(MoonScale = median(SppAccpData$MoonScale),
                          ForestScale = median(SppAccpData$ForestScale),
                          DevelScale = median(SppAccpData$DevelScale),
                          WaterScale = median(SppAccpData$WaterScale),
                          ElevScale = median(SppAccpData$ElevScale),
                          YdayScale = median(SppAccpData$YdayScale),
                          Year = "2016",
                          BrightScale = median(SppAccpData$BrightScale),
                          CoolScale = median(SppAccpData$CoolScale),
                          StrScale = seq(min(SppAccpData$StrScale), max(SppAccpData$StrScale),
                                         length.out = nrow(SppAccpData)),
                          Site = NA)

# Generating table to iterate over with data table names, scaled term vector names,
# unscaled term vector names, and labels for each numeric vector.
PredMeta <- data.table(DfName = c("Forest", "Devel", "Water", "Elev", "Bright", "Cool", "Structures"),
                       ScaleTerm = c("ForestScale", "DevelScale", "WaterScale", "ElevScale", "BrightScale",
                                     "CoolScale", "StrScale"),
                       UsTerm = c("ForestUs", "DevelUs", "WaterUs", "ElevUs", "BrightUs", "CoolUs", "StrUs"),
                       UsCol = c("ManualForestPct", "ManualDevelPct", "WaterDist", "Elev", "BrightCount",
                                 "PropCool", "StrWeight"),
                       xlabel = c("Proportion Forested Area", "Proportion Developed Area",
                                  "Distance to Water (m)", "Elevation (m)", "Brightness Index Sum",
                                  "Proportion Cool Lights", "Structure Index"))

# Plotting each species' predicted response to changes in each continuous parameter
pb <- txtProgressBar(0, (nrow(PredMeta)*length(Mods)), style = 3)
x <- 0
ResList <- list()
TabList <- list()
for(h in 1:length(Mods)){
  i <- Mods[h]
  TmpMod <- get(i)
  Group <- substr(i, 1, 4)
  
  TmpDat <- as.data.table(summary(TmpMod)$coeff$cond, keep.rownames = T)
  TmpDat[, Mod := i]
  TabList <- c(TabList, list(TmpDat))
  
  GroupMean <- paste0(Group, "Mean")
  GroupSe <- paste0(Group, "Se")
  
  YmaxVect <- numeric()
  PlotVect <- character()
  PlotList <- list()
  
  for(j in 1:nrow(PredMeta)){
    TmpData <- get(PredMeta[j, DfName])
    
    TmpData[[GroupMean]] <-  predict(TmpMod, newdata = TmpData, se.fit = T, re.form = NA, type = "response")$fit
    TmpData[[GroupSe]] <- predict(TmpMod, newdata = TmpData, se.fit = T, re.form = NA, type = "response")$se.fit
    TmpData[[PredMeta[j, UsTerm]]] <- (TmpData[[PredMeta[j, ScaleTerm]]]*(2*sd(SppAccpData[[PredMeta[j, UsCol]]])))+
      mean(SppAccpData[[PredMeta[j, UsCol]]])

    TmpPlot <- ggplot(TmpData, aes(.data[[PredMeta[j, UsTerm]]], .data[[GroupMean]])) +
      geom_line() +
      geom_ribbon(aes(ymin=(.data[[GroupMean]]-1.96*.data[[GroupSe]]),
                  ymax=(.data[[GroupMean]]+1.96*.data[[GroupSe]])), alpha=0.3) +
      geom_ribbon(aes(ymin=(.data[[GroupMean]]-1.44*.data[[GroupSe]]),
                      ymax=(.data[[GroupMean]]+1.44*.data[[GroupSe]])), alpha=0.3) +
      labs(x = PredMeta[j, xlabel], y = "Count") +
      theme_light() +
      theme(plot.title = element_text(hjust = 0.5))

    ymax <- layer_scales(TmpPlot)$y$range$range[2]
    
    TmpPlot <- TmpPlot + coord_cartesian(ylim = c(0,ymax))
    
    YmaxVect <- c(YmaxVect, ymax)
    
    PlotName <- paste0(Group, PredMeta[j, DfName], "Plot")
    PlotVect <- c(PlotVect, PlotName)

    assign(PlotName, TmpPlot)
    x <- x+1
    setTxtProgressBar(pb, x)
  }
  Ylimit <- max(YmaxVect)
  
  for(k in 1:length(PlotVect)){
    TmpPlot <- get(PlotVect[k])
    assign(PlotVect[k], TmpPlot)
    
    PlotList[[k]] <- TmpPlot
  }
  TopPlot <- plot_grid(plotlist = PlotList[1:6], scale = 0.95, ncol = 3)
  BottomPlot <- plot_grid(NULL, PlotList[7][[1]], NULL, ncol = 3, scale = 0.95)
  GridPlot <- plot_grid(TopPlot, BottomPlot, ncol = 1, rel_heights = c(2,1))
  
  GridPlot <- plot_grid({ggdraw() + draw_label(substr(i, 1, 4),
                                               hjust = 0.5, fontface = "bold", size = 18)},
                        GridPlot,ncol = 1, rel_heights = c(0.05,1))
  
  GridPlot <- GridPlot + theme(plot.background = element_rect(color = "black"))
  
  GridName <- paste0(Group, "ParamGrid")
  assign(GridName, GridPlot)
}

# Creating multi-panel plot with predicted responses of each species to changes
# in continuous parameters
Top <- plot_grid(EpfuParamGrid, LaciParamGrid, LanoParamGrid, MyevParamGrid, MyluParamGrid, MyvoParamGrid,
                 ncol = 3)
Bottom <- plot_grid(NULL, MyyuParamGrid, NULL, ncol = 3)

Both <- plot_grid(Top, Bottom, ncol = 1, rel_heights = c(2,1))

ggsave("./Output/FullParamGrid.png", plot = Both, height = 18, width = 24, dpi = 200)

# Creating model summary table
ModSums <- rbindlist(TabList)
setnames(ModSums, c("Std. Error", "Pr(>|z|)"), c("SE", "p"))
ModSums[p < 0.001, pchar := "<0.001"]
ModSums[p >= 0.001, pchar := as.character(round(p, 3))]
ModSums[,p := NULL]
ModSums[,Estimate := as.character(round(Estimate, 3))]
ModSums[,SE := as.character(round(SE, 3))]
setnames(ModSums, "pchar", "p", skip_absent = T)
ModSums[,rn := rep(Labs, 7)]
ModSums[,PredNature := rep(c("AA", rep("Anthropogenic", 3), rep("Natural", 4),
                             "Anthropogenic", rep("Natural", 2)), 7)]
setorder(ModSums, PredNature, rn)
ModSums[,PredNature := NULL]
ModSums[,Species := substr(Mod, 1, 4)]
ModSums[,c("z value", "Mod") := NULL]
ModSumsLong <- pivot_longer(ModSums, c(Estimate, SE, p),
                            names_to = "Value", values_to = "Number")
ModSumsWide <- as.data.table(pivot_wider(ModSumsLong, names_from = rn, values_from = Number))

ModSumTab <- kable(ModSumsWide, digits = 3, align = "c") %>% 
  kable_styling(font_size = 14) %>% 
  add_header_above(c(" " = 3, "Anthropogenic" = 4, "Natural" = 6)) %>% 
  column_spec(3:length(ModSumsWide), width_min = "3cm", width_max = "3cm", border_left = T) %>% 
  collapse_rows(columns = 1) %>% 
  row_spec(seq(4, nrow(ModSumsWide), 6), background = "grey") %>% 
  row_spec(seq(5, nrow(ModSumsWide), 6), background = "grey") %>% 
  row_spec(seq(6, nrow(ModSumsWide), 6), background = "grey") %>% 
  row_spec(1:nrow(ModSumsWide), color = "black")
save_kable(ModSumTab, "./Output/ModSumTab.png")

if(!"Tables" %in% dir("./Output")){dir.create("./Output/Tables")}

fwrite(ModSumsWide, "./Output/Tables/ModelSummary.csv")

# Generating table with predicted response of each species to a specified change in
# each continuous parameter.
RealEf <- RealEffectTabWide(ls(pattern = "^[[:alpha:]]{4}Nb2$"), rownames(summary(EpfuNb2)$coeff$cond)[c(2:4, 7:10)],
                  c(0.05, 0.2, 5, 50, 250, 0.2, 0.2), ScaleSds = c(2,2,2,2,2,2,2),
                  PredVects = c("StrWeight", "PropCool", "BrightCount", "Elev", "WaterDist",
                                "ManualDevelPct", "ManualForestPct"), Data = SppAccpData)
RealEf[,Species := substr(Model, 1, 4)]
RealEf[DeltaPct >= 0, PctConf := paste0("+", round(DeltaPct, 2), "% (", round(LowerConf, 2), "% / ", round(UpperConf, 2), "%)")]
RealEf[DeltaPct < 0, PctConf := paste0(round(DeltaPct, 2), "% (", round(LowerConf, 2), "% / ", round(UpperConf, 2), "%)")]
RealEf[,UnitChange := as.character(UnitChange)]

RealEfWide <- as.data.table(pivot_wider(RealEf, id_cols = c(Predictor, UnitChange),
                                        names_from = Species, values_from = PctConf))
RealEfWide[,PredictorName := c("Structure Index", "Proportion Cool Lights", "Brightness Index",
                               "Elevation (m)", "Water Distance (m)", "Proportion Developed",
                               "Proportion Forested")]
RealEfWide[,PredNature := c(rep("Anthropogenic", 3), rep("Natural", 2), "Anthropogenic", "Natural")]
RealEfWide[,Predictor := PredictorName]
RealEfWide[,PredictorName := NULL]

setcolorder(RealEfWide, c("PredNature", "Predictor"))
setorder(RealEfWide, PredNature)
setnames(RealEfWide, old = c("PredNature", "UnitChange"),
         new = c("Predictor Classification", "Unit Increase"),
         skip_absent = T)

# Writing model predictions for a given change in continuous parameters to file
fwrite(RealEfWide, "./Output/Tables/RealEfWide.csv")

RealEfWide[,`Predictor Classification` := NULL]

# Formatting a kable table and writing it to file
RealEfTab <- kable(RealEfWide, align = "c") %>% 
  kable_styling(font_size = 24) %>% 
  row_spec(seq(2, nrow(RealEfWide), 2), background = "grey") %>% 
  row_spec(1:nrow(RealEfWide), color = "black") %>% 
  group_rows("Anthropogenic", 1, 4) %>%
  group_rows("Natural", 5, 7) %>%
  column_spec(3:length(RealEfWide), width_min = "4.2cm", width_max  = "4.2cm", border_left = T)

save_kable(RealEfTab, "./Output/RealEfTab.png")
