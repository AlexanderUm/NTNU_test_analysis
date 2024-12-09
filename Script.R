#-------------------------------------------------------------------------------
# A test statistical analysis ExMilk
#-------------------------------------------------------------------------------
set.seed(947354)

# Libraries 
Libs <- c("tidyverse", "lme4", "lmerTest", 
          "ComplexHeatmap", "broom.mixed")

for(i in Libs) {
  
  require(i, character.only = TRUE)
  
}

#-------------------------------------------------------------------------------
# Variables 

DataPath <- "data/Insulin data cross-over study.csv"

DirOut <- "out"

dir.create(DirOut, recursive = TRUE)


################################################################################
# Data overview
################################################################################
Data <- read.csv(DataPath, sep = ";")

# Fix values 
DataForm <- Data %>% 
              mutate(across(everything(), 
                            function(x){as.numeric(gsub(",", ".", x))})) %>% 
              mutate(Part_ID = as.factor(paste0("p", ID))) %>% 
              select(-ID)

#-------------------------------------------------------------------------------
# Spaghetti plot 
#-------------------------------------------------------------------------------
# Convert to long
DataLong <- DataForm %>% 
               pivot_longer(-Part_ID) %>% 
               mutate(Time = as.numeric(gsub(".*_", "", name)), 
                      Treatment = factor(gsub("_.*", "", name), 
                                    levels = c("REST", "MOD", "HIIT")))


DataLongSpagg <- bind_rows(mutate(DataLong, Trans = "None"),
                           mutate(DataLong, Trans = "Log", value = log(value))) %>% 
                    mutate(Trans = factor(Trans, levels = c("None", "Log")))

DataText <- DataLongSpagg %>% 
              filter(Time == max(Time))

OverviewPlot <- ggplot(DataLongSpagg, 
                       aes(y = value, 
                           x = Time, 
                           colour = Part_ID)) + 
                  geom_point(size = 1.5, alpha = 0.5) +
                  geom_line(aes(group = Part_ID), alpha = 0.5, linewidth = 0.75) + 
                  facet_grid(Trans ~ Treatment, scales = "free") + 
                  geom_text(data = DataText, 
                            aes(label = Part_ID, 
                                x = Time, y = value), hjust = -0.2) +
                  scale_x_continuous(breaks = sort(unique(DataLong$Time))) + 
                  theme_bw() + 
                  theme(legend.position = "none", 
                        panel.grid.minor = element_blank(), 
                        axis.title = element_blank(), 
                        strip.text = element_text(face = "bold")) + 
                  expand_limits(x = c(7, 16.5))

ggsave(filename = paste0(DirOut, "/overview.png"), 
       plot = OverviewPlot, width = 7, height = 6)  


#-------------------------------------------------------------------------------
# Heat map 
#-------------------------------------------------------------------------------
DataHeat <- DataForm %>% 
              column_to_rownames(var = "Part_ID")

ColSplitFac <- factor(gsub("_.*", "", colnames(DataHeat)), 
                      levels = c("REST", "MOD", "HIIT"))

ColLabels <- gsub(".*_", "", colnames(DataHeat))

HeatPlotRaw <- Heatmap(DataHeat, 
                        name = "None",
                        rect_gp = gpar(col = "gray40", 
                                       lwd = 0.5),
                        cluster_columns = FALSE, 
                        column_split = ColSplitFac, 
                        column_labels = ColLabels, 
                        row_title = "None", 
                        show_heatmap_legend = FALSE)

HeatPlotLog <-  Heatmap(log(DataHeat), 
                        name = "Log",
                        rect_gp = gpar(col = "gray40", 
                                       lwd = 0.5),
                        cluster_columns = FALSE, 
                        column_split = ColSplitFac, 
                        column_labels = ColLabels, 
                        row_title = "Log", 
                        show_heatmap_legend = FALSE)

HeatPlotComb <- HeatPlotRaw %v% HeatPlotLog

png(paste0(DirOut, "/heat_map.png"), 
    width = 6, height = 7, res = 300, units = "in")
draw(HeatPlotComb)
dev.off()

################################################################################
# Analysis type visualization 
################################################################################
#-------------------------------------------------------------------------------
# 1. Does overall dynamic of breastmilk insulin concentrations differ 
# between treatments?
#-------------------------------------------------------------------------------

Type1Plot <- DataLong %>% 
              mutate(Line_Group = paste0(Treatment, "_",Part_ID), 
                     Set = paste(sort(unique(Time)), collapse = " -> ")) %>% 
              ggplot(aes(y = log(value), 
                         x = Time, 
                         colour = Treatment)) + 
                    geom_point(size = 1.5, alpha = 0.5) +
                    geom_line(aes(group = Line_Group), alpha = 0.5, linewidth = 0.75) + 
                    scale_x_continuous(breaks = sort(unique(DataLong$Time))) + 
                    theme_bw() + 
                    theme(panel.grid.minor = element_blank(), 
                          strip.text = element_text(face = "bold"), 
                          axis.title = element_blank()) + 
                    scale_color_brewer(palette = "Dark2") + 
                    facet_wrap("Set")

ggsave(plot = Type1Plot, 
       filename = paste0(DirOut, "/model_type1.png"), 
       width = 5.5, height = 2)  


#-------------------------------------------------------------------------------
# 2. Is there a difference in rise and fall dynamic of breastmilk insulin 
# concentrations in response to meal consumption between treatments?
# 3. Are there differences in breastmilk insulin concentrations in 
# comparison with the baseline?
#-------------------------------------------------------------------------------
TimeSets <- list()

TimeSets[["Type2"]] <- list(c(7, 11), c(11, 12), c(12, 15))

TimeSets[["Type3"]] <- list(c(7, 11), c(7, 12), c(7, 15))

for(i in names(TimeSets)) {
  
  iDataSets <- NULL
  
  for(j in TimeSets[[i]]) {
    
    iDataSets <-  DataLong %>% 
                      filter(Time %in% j) %>% 
                      mutate(Set = paste(j, collapse = " -> ")) %>% 
                      droplevels() %>% 
                      bind_rows(iDataSets, .)
    
  }
  
  iTypePlot <- iDataSets %>% 
                        mutate(Line_Group = paste0(Treatment, "_",Part_ID)) %>%
                  ggplot(aes(y = log(value), 
                                       x = Time, 
                                       colour = Treatment)) + 
                    geom_point(size = 1.5, alpha = 0.5) +
                    geom_line(aes(group = Line_Group), alpha = 0.5, linewidth = 0.75) + 
                    facet_grid(.~Set, scales = "free_x") + 
                    scale_x_continuous(breaks = sort(unique(iDataSets$Time))) + 
                    theme_bw() + 
                    theme(panel.grid.minor = element_blank(), 
                          strip.text = element_text(face = "bold"), 
                          axis.title = element_blank()) + 
                    scale_color_brewer(palette = "Dark2")
  
  ggsave(plot = iTypePlot, 
         filename = paste0(DirOut, "/model_", i, ".png"), 
         width = 5.5, height = 2) 
  
}


################################################################################
# Statistical testing - LMM
################################################################################
# 1. Does overall dynamic of breast milk insulin concentrations differ 
# between treatments?
################################################################################
DataSetsLs <- list()

DataSetsLs[["All"]] <- DataLong

DataSetsLs[["Sub1"]] <- DataLong %>% 
                          filter(!Part_ID %in% c("p1", "p19")) %>% 
                          droplevels()

# Transformation functions
TransFun <- list()

TransFun[["None"]] <- function(x){x}

TransFun[["Log"]] <- function(x){log(x)}

# Parameters grid
LmmPrm <- expand.grid("Rand_effect" = c("(1|Part_ID)", 
                                        "(1+Time|Part_ID)", 
                                        "(1|Part_ID/Treatment)",
                                        "(1+Time|Part_ID/Treatment)"), 
                      "Samples_Set" = names(DataSetsLs), 
                      "Y_Tras" = names(TransFun), 
                      stringsAsFactors = FALSE)

FitRes <- NULL

ModelRes <- NULL

FullModelsLs <-list()

for(i in 1:nrow(LmmPrm)) {
  
  # Variables 
  iTransFun <- TransFun[[LmmPrm[i, "Y_Tras"]]]
  
  iData <- DataSetsLs[[LmmPrm[i, "Samples_Set"]]] %>% 
            mutate(across(value, iTransFun))
  
  iFormual <- paste0("value ~ Time*Treatment + ", LmmPrm[i, "Rand_effect"])
  
  # Model 
  Model <- lmerTest::lmer(as.formula(iFormual), data = iData)
  
  FullModelsLs[[paste(LmmPrm[i, ], collapse = "--")]] <- Model
  
  ModelRes <- Model %>% 
                broom.mixed::tidy() %>% 
                bind_cols(LmmPrm[i, ]) %>% 
                bind_rows(ModelRes, .)
  
  # Assumptions data 
  FitRes <- data.frame(Fitted = fitted(Model), 
                       Residuals = resid(Model)) %>% 
              bind_cols(LmmPrm[i, ]) %>% 
              bind_rows(FitRes, .)

}


#-------------------------------------------------------------------------------
# Visualize models assumptions
#-------------------------------------------------------------------------------
FitResPlot <- ggplot(FitRes, aes(y = Residuals, x = Fitted)) + 
                      geom_point(color = "steelblue", alpha = 0.5) + 
                      geom_abline(intercept = 0, slope = 0) + 
                      facet_wrap(c("Rand_effect", "Samples_Set", "Y_Tras"), 
                                 scales = "free", 
                                 ncol = 4) + 
                      theme_bw() 

ggsave(plot = FitResPlot, 
       filename = paste0(DirOut, "/Full_Fitted_vs_Residual.png"), 
       width = 12, 
       height = 9) 


QqResPlot <- ggplot(FitRes, aes(sample = Residuals)) + 
                    stat_qq(color = "steelblue", alpha = 0.5) + 
                    stat_qq_line() + 
                    facet_wrap(c("Rand_effect", "Samples_Set", "Y_Tras"), 
                               scales = "free", 
                               ncol = 4) + 
                    theme_bw()

ggsave(plot = QqResPlot, 
       filename = paste0(DirOut, "/Full_QQ_Residual.png"), 
       width = 12, 
       height = 9)

#-------------------------------------------------------------------------------
# Compare models 
#-------------------------------------------------------------------------------
# All samples log 
anova(FullModelsLs$`(1|Part_ID)--All--Log`, 
      FullModelsLs$`(1+Time|Part_ID)--All--Log`, 
      FullModelsLs$`(1|Part_ID/Treatment)--All--Log`, 
      FullModelsLs$`(1+Time|Part_ID/Treatment)--All--Log`) %>% 
  tidy() %>% 
  mutate(term = gsub(".*\\$|\\`", "", term), 
         across(where(is.numeric), function(x){round(x, 3)})) %>% 
  separate_wider_delim(term, 
                       delim =  "--", 
                       names = c("Random", "Samples", "Trans")) %>% 
  write.csv(paste0(DirOut, "/anova_all_samples.csv"), row.names = FALSE)

# Set1 samples log 
anova(FullModelsLs$`(1|Part_ID)--Sub1--Log`, 
      FullModelsLs$`(1+Time|Part_ID)--Sub1--Log`, 
      FullModelsLs$`(1|Part_ID/Treatment)--Sub1--Log`, 
      FullModelsLs$`(1+Time|Part_ID/Treatment)--Sub1--Log`) %>% 
  tidy() %>% 
  mutate(term = gsub(".*\\$|\\`", "", term), 
         across(where(is.numeric), function(x){round(x, 3)})) %>% 
  separate_wider_delim(term, 
                       delim =  "--", 
                       names = c("Random", "Samples", "Trans")) %>% 
  write.csv(paste0(DirOut, "/anova_sub_samples.csv"), row.names = FALSE)


# Select most promising models results 
ModelRes %>% 
   mutate(across(where(is.numeric), function(x) {round(x, 3)})) %>% 
   filter(Rand_effect %in% c("(1+Time|Part_ID/Treatment)",
                          "(1|Part_ID/Treatment)"), 
          Samples_Set == "All", 
          Y_Tras == "Log") %>% 
  write.csv(paste0(DirOut, "/full_best_mod_res.csv"), row.names = FALSE)
          

################################################################################
# Statistical testing - LMM
################################################################################
# 2. Is there a difference in rise and fall dynamic of breast milk insulin 
# concentrations in response to meal consumption between treatments?
# 3. Are there differences in breast milk insulin concentrations in 
# comparison with the baseline?
################################################################################

ModelSumLs <- list()

ModeLs <- list()

FitResLs <- list()


for (i in names(TimeSets)){
  
  iSubSet <- TimeSets[[i]] %>% 
              setNames(lapply(., function(x){paste(x, collapse = " -> ")}))
  
  iGrid <- expand.grid("Rand_effect" = c("(1|Part_ID)", 
                                         "(1+Time|Part_ID)", 
                                         "(1|Part_ID/Treatment)"), 
                       "Sub_Set" = names(iSubSet),
                       stringsAsFactors = FALSE)
  
  iModelSum <- NULL
  
  iModeLs <- list()
  
  iFitRes <- NULL
  
  for(j in 1:nrow(iGrid)) {
    
    jSubName <- iGrid[j, "Sub_Set"]
    
    jFormula <- paste0("value ~ Time*Treatment + ", 
                       iGrid[j, "Rand_effect"])
    
    jDataSets <-  DataLong %>% 
                      filter(Time %in% iSubSet[[jSubName]]) %>% 
                      mutate(value = log(value)) %>% 
                      droplevels() 
    
    jModel <- lmerTest::lmer(as.formula(jFormula), 
                            data = jDataSets)
    
    iModeLs[[jSubName]][[iGrid[j, "Rand_effect"]]] <- jModel
    
    iModelSum <- jModel %>% 
                   broom.mixed::tidy() %>% 
                   bind_cols(iGrid[j, ]) %>% 
                   mutate(Set_Type = i) %>% 
                   bind_rows(iModelSum, .)
    
    # Assumptions data 
    iFitRes <- data.frame(Fitted = fitted(jModel), 
                         Residuals = resid(jModel)) %>% 
                  bind_cols(iGrid[j, ]) %>% 
                  bind_rows(iFitRes, .)
    
  }
  
  ModelSumLs[[i]] <- iModelSum
  
  ModeLs[[i]] <- iModeLs
  
  #-------------------------------------------------------------------------------
  # Visualize models assumptions
  #-------------------------------------------------------------------------------
  iFitResPlot <- ggplot(iFitRes, aes(y = Residuals, x = Fitted)) + 
                    geom_point(color = "steelblue", alpha = 0.5) + 
                    geom_abline(intercept = 0, slope = 0) + 
                    facet_wrap(c("Rand_effect", "Sub_Set"), 
                               scales = "free", 
                               ncol = 3) + 
                    theme_bw() 
  
  ggsave(plot = iFitResPlot, 
         filename = paste0(DirOut, "/", i, "_Fitted_vs_Residual.png"), 
         width = 9, 
         height = 9) 
  
  
  iQqResPlot <- ggplot(iFitRes, aes(sample = Residuals)) + 
                      stat_qq(color = "steelblue", alpha = 0.5) + 
                      stat_qq_line() + 
                      facet_wrap(c("Rand_effect", "Sub_Set"), 
                                 scales = "free", 
                                 ncol = 3) + 
                      theme_bw() 
  
  ggsave(plot = iQqResPlot, 
         filename = paste0(DirOut, "/", i,"_QQ_Residual.png"), 
         width = 9, 
         height = 9)
  
}


#-------------------------------------------------------------------------------
# Compare models 
#-------------------------------------------------------------------------------
SetModComp <- NULL

for(i in names(ModeLs)) {
  
  for(j in names(ModeLs[[i]])) {
    
    SetModComp <- anova(ModeLs[[i]][[j]][[1]], 
                            ModeLs[[i]][[j]][[2]], 
                            ModeLs[[i]][[j]][[3]]) %>% 
                        tidy() %>% 
                        mutate(term = names(ModeLs[[i]][[j]]), 
                               Type = i, 
                               Set = j) %>% 
                        mutate(across(where(is.numeric), 
                                      function(x){round(x, 3)})) %>% 
                        add_row() %>% 
                        bind_rows(SetModComp, .)
  }
  
}

write.csv(SetModComp, paste0(DirOut, "/sub_mod_comp.csv"), 
          row.names = FALSE, na = "")


#-------------------------------------------------------------------------------
# Write out results 
#-------------------------------------------------------------------------------
for(i in names(ModelSumLs)) {
  
  ModelSumLs[[i]] %>% 
    filter(effect == "fixed") %>% 
    mutate(across(where(is.numeric), function(x){round(x, 3)})) %>% 
    select(-c(effect, group)) %>% 
    filter(grepl("Time:", term)) %>% 
  write.csv(paste0(DirOut, "/", i, "sub_mod_res.csv"), 
            row.names = FALSE, na = "")
  
}
