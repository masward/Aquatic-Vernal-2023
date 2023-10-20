#Mason Ward
##The Pennsylvania State University
###October 10, 2023
####Vernal Pond Environmental Variables

#Load in applicable libraries
library(dplyr)
library(ggplot2)
library(tibble)
library(vegan)
library(tidyverse)
library(reshape2)
library(viridis)
library(ggrepel)
library(rstatix)
library(RColorBrewer)
library(glmmTMB)
library(DataExplorer)

#bring in dataset from Github
##Username - masward
###Spongy-Moth repository
environmental_df <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/Vernal_Pond_Datasheet_Oct102023.csv")

#Change week, site, and wet/dry from characters to factors ----
environmental_df$WEEK_COL <- as.factor(environmental_df$WEEK_COL)
environmental_df$SITE_ID <- as.factor(environmental_df$SITE_ID)
environmental_df$WET_OR_DRY <- as.factor(environmental_df$WET_OR_DRY)



#Kruskal Wallis tests between SITE_ID and each variable on the raw data. This allows to determine ----
##if sites differ among variables statistically 
kruskal.test(PH ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 30.79, df = 7, p-value = 6.797e-05
kruskal.test(CONDUCTIVITY_uS ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 20.031, df = 7, p-value = 0.005502 *
kruskal.test(TEMPERATURE_C ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 11.162, df = 7, p-value = 0.1317
kruskal.test(DEPTH_CM ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 6.1182, df = 7, p-value = 0.526
kruskal.test(LENGTH_M ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 8.0268, df = 7, p-value = 0.3302
kruskal.test(WIDTH_M ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 27.645, df = 7, p-value = 0.0002549 *
kruskal.test(PCT_CANOPY_COVER ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 63.854, df = 7, p-value = 2.554e-11
kruskal.test(TOC_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 11.356, df = 6, p-value = 0.07798
kruskal.test(TN_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 11.137, df = 6, p-value = 0.08422
kruskal.test(DOC_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 11.28, df = 6, p-value = 0.08011
kruskal.test(TDN_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 7.0037, df = 6, p-value = 0.3205
kruskal.test(CHLORIDE_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 18.044, df = 6, p-value = 0.006124 *
kruskal.test(NITRATE_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 8.5485, df = 6, p-value = 0.2006
kruskal.test(SULFATE_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 10.889, df = 6, p-value = 0.09186
kruskal.test(PHOSPHATE_mg_L ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 7.8458, df = 5, p-value = 0.1649
kruskal.test(FRASS_g ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 10.032, df = 7, p-value = 0.1867
kruskal.test(LEAF_LITTER_g ~ SITE_ID, data = environmental_df) #Kruskal-Wallis chi-squared = 12.375, df = 7, p-value = 0.08889
##dunnTest test afterwards for significant results — on significant results


##Use "plot_histogram" to run histograms on all environmental variables to see if they are normally distributed ----
histogram_plot <- plot_histogram(environmental_df, binary_as_factor = TRUE, geom_histogram_args = list(bins = 30L),
  scale_x = "continuous", title = NULL, ggtheme = theme_gray(), theme_config = list(),
  nrow = 4L, ncol = 5L, parallel = FALSE) +
  fig_theme()


##Data is not normal — need to log transform data to normalize it (Except pH since it is already log)
env.log.data <- environmental_df %>% #Selecting which columns to  include
  group_by(WEEK_COL, SITE_ID, WET_OR_DRY, PH) %>% #grouping by variables
  transmute(log_conductivity_uS = log(CONDUCTIVITY_uS),
            log_temperature_C = log(TEMPERATURE_C), 
            log_depth_cm = log(DEPTH_CM),
            log_length_M = log(LENGTH_M),
            log_width_M = log(WIDTH_M),
            log_canopy_pct = log(PCT_CANOPY_COVER),
            log_TOC_mg_L = log(TOC_mg_L),
            log_TN_mg_L = log(TN_mg_L),
            log_DOC_mg_L = log(DOC_mg_L),
            log_TDN_mg_L = log(TDN_mg_L),
            log_chloride_mg_L = log(CHLORIDE_mg_L),
            log_nitrate_mg_L = log(NITRATE_mg_L),
            log_sulfate_mg_L = log(SULFATE_mg_L),
            log_phosphate_mg_L = log(PHOSPHATE_mg_L),
            log_frass_g = log(FRASS_g),
            log_leaf_g = log(LEAF_LITTER_g)) %>% 
  distinct()

log_hist <- plot_histogram(env.log.data, binary_as_factor = TRUE, geom_histogram_args = list(bins = 30L),
               scale_x = "continuous", title = NULL, ggtheme = theme_gray(), theme_config = list(),
               nrow = 4L, ncol = 5L, parallel = FALSE) +
  fig_theme()

##global model — based on output we can remove the variables, glmm — include time as a factor (blocking factor) ----










