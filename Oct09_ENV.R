#Mason S Ward
#2023-10-04
#Vernal pond dataset
#The Pennsylvania State University

#Load applicable libraries
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

#bring in dataset from Github - masward (username) under the "Spongy-Moth" repository
ENV.DF <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/Vernal_Pond_Datasheet_10042023.csv")
#Check dataframe — everything looks fine so far

#SEPARATE FRASS & LEAF-LITTER — these will be cumulative sums temporally, not averages ----

#AVERAGE & FIND STANDARD DEVS FOR ALL ENVIRONMENTAL VARIABLES BY SITE - EXCEPT FOR FRASS & LEAF LITTER ----
AVG_ENV <- ENV.DF %>%
  select(-c(WET_OR_DRY)) %>% #Selecting which column to not include
  group_by(SITE_ID) %>% #grouping by SITE_ID
  transmute(avg_PH = mean(PH, na.rm = TRUE), #Averaging values across different physical/chemical variables & not including NA values (pond dry)
            avg_CONDUCTIVITY_uS = mean(CONDUCTIVITY_uS, na.rm = TRUE),
            avg_TEMPERATURE_C = mean(TEMPERATURE_C, na.rm = TRUE),
            avg_DEPTH_CM = mean(DEPTH_CM, na.rm = TRUE), 
            avg_LENGTH_M = mean(LENGTH_M, na.rm = TRUE), 
            avg_WIDTH_M = mean(WIDTH_M, na.rm = TRUE),
            avg_CANOPY = mean(PCT_CANOPY_COVER, na.rm = TRUE),
            avg_TOC_mg_L = mean(TOC_mg_L, na.rm = TRUE),
            avg_TN_mg_L = mean(TN_mg_L, na.rm = TRUE),
            avg_DOC_mg_L = mean(DOC_mg_L, na.rm = TRUE), 
            avg_TDN_mg_L = mean(TDN_mg_L, na.rm = TRUE), 
            avg_CHLORIDE_mg_L = mean(CHLORIDE_mg_L, na.rm = TRUE), 
            avg_NITRATE_mg_L = mean(NITRATE_mg_L, na.rm = TRUE),
            avg_SULFATE_mg_L = mean(SULFATE_mg_L, na.rm = TRUE), 
            avg_PHOSPHATE_mg_L = mean(PHOSPHATE_mg_L, na.rm = TRUE),
            sd_PH = sd(PH, na.rm = TRUE),
            sd_CONDUCTIVITY_uS = sd(CONDUCTIVITY_uS, na.rm = TRUE),
            sd_TEMPERATURE_C = sd(TEMPERATURE_C, na.rm = TRUE),
            sd_DEPTH_CM = sd(DEPTH_CM, na.rm = TRUE), 
            sd_LENGTH_M = sd(LENGTH_M, na.rm = TRUE),
            sd_WIDTH_M = sd(WIDTH_M, na.rm = TRUE),
            sd_CANOPY = sd(PCT_CANOPY_COVER, na.rm = TRUE),
            sd_TOC_mg_L = sd(TOC_mg_L, na.rm = TRUE),
            sd_TN_mg_L = sd(TN_mg_L, na.rm = TRUE),
            sd_DOC_mg_L = sd(DOC_mg_L, na.rm = TRUE), 
            sd_TDN_mg_L = sd(TDN_mg_L, na.rm = TRUE), 
            sd_CHLORIDE_mg_L = sd(CHLORIDE_mg_L, na.rm = TRUE), 
            sd_NITRATE_mg_L = sd(NITRATE_mg_L, na.rm = TRUE),
            sd_SULFATE_mg_L = sd(SULFATE_mg_L, na.rm = TRUE), 
            sd_PHOSPHATE_mg_L = sd(PHOSPHATE_mg_L, na.rm = TRUE)) %>% 
  distinct()

ENV_long_data <- pivot_longer(
  data = AVG_ENV,
  cols = -SITE_ID,
  names_to = c(".value", "Variable"),
  names_pattern = "^(avg|sd)_(.*)")

#Facet wrap of environmental variables in ggplot ----
FACET <-ggplot(ENV_long_data, aes(x= SITE_ID, y= avg)) +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), 
                width=0.2, 
                position = position_dodge(0.05)) +
  geom_point() +
  facet_wrap(~ Variable, scales = "free_y") +
  fig_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#Figure theme from Christine - run before GGPLOT ----
fig_theme <- function(){
  theme_bw() +
    theme(text = element_text(family = "Times New Roman",
                              color = "black"),
          axis.text = element_text(color = "black",
                                   size = 12),
          axis.title = element_text(size = 16),
          panel.grid = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
}

#Test which variables are statistically different between sites
kruskal.test(avg_PH ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_CONDUCTIVITY_uS ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_TEMPERATURE_C ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_DEPTH_CM ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_LENGTH_M ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_WIDTH_M ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_CANOPY ~ SITE_ID, data = AVG_ENV) #P value is 0.4289
kruskal.test(avg_TOC_mg_L ~ SITE_ID, data = AVG_ENV) #P value is 0.4232
kruskal.test(avg_TN_mg_L ~ SITE_ID, data = AVG_ENV) #P value is 0.4232
kruskal.test(avg_DOC_mg_L ~ SITE_ID, data = AVG_ENV) #P value is 0.4232
kruskal.test(avg_TDN_mg_L ~ SITE_ID, data = AVG_ENV) #P value is 0.4232
kruskal.test(avg_CHLORIDE_mg_L ~ SITE_ID, data = AVG_ENV) #P value is 0.4232
kruskal.test(avg_NITRATE_mg_L ~ SITE_ID, data = AVG_ENV)#P value is 0.4232
kruskal.test(avg_SULFATE_mg_L ~ SITE_ID, data = AVG_ENV)#P value is 0.4232
kruskal.test(avg_SULFATE_mg_L ~ SITE_ID, data = AVG_ENV)#P value is 0.4232
kruskal.test(avg_PHOSPHATE_mg_L ~ SITE_ID, data = AVG_ENV) #P value is 0.4159



##PULL OUT CANOPY COVER & FRASS FOR LINEAR REGRESSION ----
LINEARREG_CC_FR <- ENV.DF %>%
  select(-c(WET_OR_DRY)) %>% #Selecting which column to not include
  group_by(SITE_ID) %>% #grouping by SITE_ID
  transmute(sum_FRASS_g = sum(FRASS_g, na.rm = TRUE),
            sum_LEAF_LITTER_g = sum(LEAF_LITTER_g, na.rm = TRUE), 
            avg_CANOPY = mean(PCT_CANOPY_COVER, na.rm = TRUE),
            sd_CANOPY = sd(PCT_CANOPY_COVER, na.rm = TRUE)) %>% 
  distinct()

linear_frass <- lm(sum_FRASS_g + sum_LEAF_LITER_g~ avg_CANOPY, data = LINEARREG_CC_FR)
summary(linear_frass) #P value 0.0622
linear_leaf <- lm(sum_LEAF_LITER_g ~ avg_CANOPY, data = LINEARREG_CC_FR)
summary(linear_leaf) #P value 0.029

LIN_FRA <- ggplot(LINEARREG_CC_FR, aes(x = avg_CANOPY, y = sum_FRASS_g, colour = SITE_ID)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE, color='black', aes(color=SITE_ID)) +
  theme_minimal() +
  labs(x='Average Percent Canopy Cover', y='Cumulative Frass (g)', title='Linear Regression Plot') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) 
LIN_LL <- ggplot(LINEARREG_CC_FR, aes(x = avg_CANOPY, y = sum_LEAF_LITER_g, colour = SITE_ID)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE, color='black', aes(color=SITE_ID)) +
  theme_minimal() +
  labs(x='Average Percent Canopy Cover', y='Cumulative Leaf Litter (g)', title='Linear Regression Plot') +
  theme(plot.title = element_text(hjust=0.5, size=20, face='bold')) 

linear_reg_plot <- ggplot(data=LINEARREG_CC_FR, aes(x=avg_CANOPY, y=sum_FRASS_g, )) + 
  geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(colour = "Frass")) +
  guides(color = guide_legend(title = "Terrestrial input")) +
  #Now add a second layer, with same x, but other y (and blue color for clarity)
  geom_point(aes(y = sum_LEAF_LITER_g)) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(y=sum_LEAF_LITER_g, color = "Leaf litter")) +
  scale_colour_manual(name="legend", values=c("#D95F02", "#56B4E9")) +
  xlab("Average canopy cover") +
  ylab("Cumlative terrestrial inputs (g)") +
  ggtitle("Relationship between terrestrial inputs and canopy cover") +
  labs(caption = "Figure 1: Relationship between terrestrial inputs (frass and leaf litter) and average canopy cover across vernal ponds (n = 8) (p-value = 0.0249).") +
  scale_linetype_discrete(labels=c("Leaf litter", "Frass")) +
  theme_bw() +
  theme(text = element_text(family = "Times New Roman", color = "black"),
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(size = 11),
        panel.grid = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        plot.caption.position = "plot", plot.caption = element_text(hjust =0, size = 11))
print(linear_reg_plot)

REG_LONG <- pivot_longer(
  data = LINEARREG_CC_FR,
  cols = -SITE_ID,
  names_to = c(".value", "Variable"),
  names_pattern = "^(avg|sum)_(.*)")



# IDENTIFY MULTICOLLINEARITY - identify the correlation between the outcome variables. ----
## A correlation above 0.9 is an indication of multicollinearity, which is problematic for MANOVA.
### 
library(corrplot)
colliniarity <- cor(AVG_ENV[,2:16], use="na.or.complete")
corrplot(colliniarity, method= "circle")
as.data.frame(colliniarity)
#Variables that are correlated 0.9 - Phosphate ~ Conductivity, Length ~ Width, Canopy ~ Sulfate, TOC ~ DOC, Chloride ~ Phosphate
##Variables that are on the verge of 0.9 (~0.8)- PH ~ Temp, PH ~ TDN, Conductivity ~ Chloride, Depth ~ DOC, Depth ~ TOC, Length ~ TOC, Length ~ DOC, TOC ~ Phosphate, TN ~ Phosphate, DOC ~ Phosphate, Nitrate ~ TDN,  