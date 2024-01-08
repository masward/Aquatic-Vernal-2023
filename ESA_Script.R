#Mason S Ward
#email: msw5688@psu.edu
#Co authors: Jon Sweetman (PhD), Sara Hermann (PhD), Alice Belskis (MS)
#Penn State University
#Department of Ecosystem Science & Management
#October 2023
#Research question: Does spongy moth influence macroinvertebrate diversity across vernal ponds in central Pennsylvania?

#Load all applicable libraries
library(ggplot2)
library(dplyr)
library(tibble)
library(vegan)
library(tidyverse)
library(reshape2)
library(viridis)
library(ggrepel)
library(rstatix)
library(adespatial)
library(betapart)
library(ecolTest)
library(flextable)
library(lme4)
library(RColorBrewer)
library(glmmTMB)
library(DataExplorer)
library(maps)
library(ggplot2) # used for plotting our maps
library(mapdata)# contains some higher-resolution outlines used in maps
library(tidyverse)
library(sf)
library(terra)
library(RColorBrewer)
library(LaplacesDemon)
library(AICcmodavg)
library(ggspatial)
library(pscl) #McFadden's R squared for glm

#bring in dataset from Github
##Username - masward -> Repository -> Spongy-Moth 
environmental_df <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/Vernal_Pond_Datasheet_Oct102023.csv")
macro.df <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/macroinvertebrate_datasheet_oct30.csv")

#Load in figure theme: Author (Christine Cornish) ----
fig.theme <- theme(text = element_text(family = "Times New Roman", color = "black"),
                   axis.text = element_text(color = "black", size = 10),
                   axis.title = element_text(size = 12),
                   panel.grid = element_blank(),
                   legend.text = element_text(size = 11),
                   legend.title = element_text(size = 11),
                   plot.caption.position = "plot", plot.caption = element_text(hjust =0, size = 11))


#Using shapefiles to create map of sites ----
##Color blind friendly palette 
color.blind.friendly <- c("#E5F5F9", "#1D91C0", "#67001F", "#F7FCFD", "#CB181D", "#78C679", "#F46D43", "#A6CEE3",
                          "#FD8D3C", "#A6D854", "#D4B9DA", "#6A51A3", "#7F0000", "#D9D9D9", "#FFF7BC", "#000000",
                          "#F0F0F0", "#C7EAE5", "#003C30", "#F16913", "#FFF7FB", "#8C6BB1", "#C7E9B4", "#762A83",
                          "#FC9272", "#AE017E", "#F7F7F7", "#DF65B0", "#EF3B2C", "#74C476")

##Create dateframe of vernal pond sites
vernal.map_df <- data.frame(state_forest  = c("Bald_Eagle", "Bald_Eagle","Bald_Eagle", "Moshannon", "Moshannon","Moshannon", "Sproul", "Sproul"),
                            SITE_ID  = c("BDEG_1", "BDEG_2", "BDEG_3", "MSHN_1", "MSHN_2", "MSHN_3", "SPRL_1", "SPRL_3"),
                            lat = c("40.799645", "40.80081", "40.80368", "40.84595", "40.93583", "40.94993", "41.10381", "41.15978"),
                            long = c("-77.512437", "-77.506195", "-77.499264", "-78.097878", "-78.032284", "-78.002096", "-77.988083", "-77.899456"))

##Use shapefile to create map!
pond_points <- st_as_sf(vernal.map_df, coords = c("long", "lat"), crs = 4269)  #4269 is the code for NAD83 obtained from here https://epsg.io/4269 


shp_forest <- st_read("/Users/masonward/Downloads/DCNR_BOF_StateForests202306/DCNR_BOF_StateForests202306.shp") #shape file of state forests
shp_forest #display information of the shapefile
st_crs(shp_forest)
st_bbox(shp_forest)
shp_boundaries <- st_read("/Users/masonward/Downloads/DCNR_BOF_Bndry_SFM201703/DCNR_BOF_Bndry_SFM201703.shp") #shape file of the boundaries of state forests
shp_boundaries #display information of the shapefile
st_crs(shp_boundaries)
st_bbox(shp_boundaries)

#Subset state forests by study area (both forests and boundaries)
shp_forest_subset <- subset(shp_forest, SF_Name %in% c('Moshannon', 'Bald Eagle','Sproul')) 
shp_boundaries_subset <- subset(shp_boundaries, DistrictNa %in% c('Moshannon', 'Bald Eagle', 'Sproul'))

#Plot and check how everything is looking
(ForestMap <- ggplot() +geom_sf(data = shp_forest, aes(fill = SF_Name)) +
    geom_sf(data = shp_boundaries) +
    ggtitle("Map of Pennsylvania State Forests") +
    ylab(expression("Latitude ("*degree*")" )) + 
    xlab(expression("Longitude ("*degree*")" )) +
    guides(fill=guide_legend(title="State Forests")) +
    fig.theme)

ggplot() +
  geom_sf(data = pond_points)

#Overlay study sites with the subset of state forests/boundaries
(VernalMap <- ggplot() + geom_sf(data = shp_forest_subset, aes(fill = SF_Name)) + 
    geom_sf(data = shp_boundaries_subset, aes(fill = DistrictNa)) + 
    geom_sf(data = pond_points) + 
    scale_fill_manual(values = c("#E5F5F9", "#D4B9DA", "#D9D9D9")) + 
    ggtitle("Map of vernal pond study sites") +
    theme_bw() + 
    ylab(expression("Latitude ("*degree*")" )) + 
    xlab(expression("Longitude ("*degree*")" )) +
    guides(fill=guide_legend(title="State Forests")) +
    coord_sf() +
    annotation_scale()) +
  fig.theme
#Edit map in Powerpoint to add on identification tags

#Start with macroinvertebrate data ----
##Alpha diversity 
###data is pulled from individual sample numbers, sum by week, site, and target taxon
## Sum by SITE_ID and like taxa
macro.abund <- macro.df %>%
  group_by(SITE_ID, TARGET_TAXON) %>% #Group by SITE_ID & TARGET_TAXON
  summarize(TAXA_ABUND = sum(TOTAL)) #Sum values together and create a TAXA_ABUND column

#Convert data to wide format for vegan analyses (vegan likes wide format for diversity indices)
macro.wide <- macro.abund %>%
  spread(key = "TARGET_TAXON",     # key = column to go on top
         value = "TAXA_ABUND") %>%     # value = column to spread under each new column
  mutate(across(everything(), ~replace_na(., 0))) # replaces all the NAs with 0s


#Barplot of Shannon Diversity for ESA
barplot_shannon <- data.frame(SITE_ID  = c("BDEG_2", "BDEG_3", "MSHN_1", "MSHN_2", "MSHN_3", "SPRL_1", "SPRL_3"),
                              Shannon_Diversity= c("2.04", "1.93", "2.12", "1.26", "1.38", "1.91", "1.64"))
barplot_shannon[2:3] <- lapply(barplot_shannon[2:3], as.numeric)

ggplot(data = barplot_shannon, aes(x = SITE_ID, y = Shannon_Diversity, fill = SITE_ID)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values = c("#E5F5F9", "#E5F5F9", "#D4B9DA","#D4B9DA", "#D4B9DA","#D9D9D9", "#D9D9D9")) +
  ylab("Shannon Diversity") +
  xlab("Site identification") +
  fig.theme

####Use Hutcheson's multiple t test to test alpha diversity site x site
alpha <- macro.wide #create separate dataframe for alpha diversity; KEEP macro.wide for beta diversity !!!!
 #remove first column since it is not important
alpha.pivot <- as.data.frame(t(alpha)) #Pivot dataframe for Hutcheson's t test
colnames(alpha.pivot) <- unlist(alpha.pivot[row.names(alpha.pivot)== 'SITE_ID',]) #Move Site_ID to column names
alpha.pivot <- alpha.pivot[-1,] #Remove the first row
alpha.pivot[] <- lapply(alpha.pivot, as.numeric) #make df numeric for H t test
results.ttest <- multiple_Hutcheson_t_test(alpha.pivot)
results.ttest <- as.data.frame(results.ttest) #There are significant differences between Sites for Shannon diversity values

##Shannon index is the same as what vegan calculated so it is correct. 
alpha.diversity <- data.frame(SITE_ID  = c("BDEG_2", "BDEG_3", "MSHN_1", "MSHN_2", "MSHN_3", "SPRL_1", "SPRL_3"),
                              Shannon_Diversity= c("2.04", "1.93", "2.12", "1.26", "1.38", "1.91", "1.64"))
alpha.diversity$Shannon_Diversity <- as.numeric(alpha.diversity$Shannon_Diversity)

###Can calulate Shannon diversity by week:
## Add back in WEEK_COL into a separate df; Sum by WEEK_COL, SITE_ID and like taxa
macro.abund_week <- macro.df %>%
  group_by(WEEK_COL, SITE_ID, TARGET_TAXON) %>% #Group by SITE_ID & TARGET_TAXON
  summarize(TAXA_ABUND = sum(TOTAL)) #Sum values together and create a TAXA_ABUND column

#Convert data to wide format for vegan analyses (vegan likes wide format for diversity indices)
macro.wide_week <- macro.abund_week %>%
  spread(key = "TARGET_TAXON",     # key = column to go on top
         value = "TAXA_ABUND") %>%     # value = column to spread under each new column
  mutate(across(everything(), ~replace_na(., 0)))     # replaces all the NAs with 0s

alpha.week <- macro.wide_week
#Convert data to wide format for vegan analyses (vegan likes wide format for diversity indices)
alpha.week <- macro.abund_week %>%
  spread(key = "TARGET_TAXON",     # key = column to go on top
         value = "TAXA_ABUND") %>%     # value = column to spread under each new column
  mutate(across(everything(), ~replace_na(., 0)))     # replaces all the NAs with 0s
###Use vegan for diversity index 
diversity.week <- diversity(alpha.week[3:37], index = "shannon")
print(diversity.week)

###Place values into new df
week.diversity <- alpha.week
week.diversity <- week.diversity[-c(3:37)]
week.diversity$shannon <- data.frame(shannon = c("1.2130076", "1.4750763", "2.1390063", "1.1262820", "1.6387253", "1.3321790", "1.6152233", "2.1547832",
                                                 "1.3293417", "0.7298429", "1.2698818", "1.0397208", "0.6931472", "1.7682782", "0.6551416",
                                                 "1.9112881", "1.6335430"))                                                
names(week.diversity)[3] <- "ShannonDiv"


###Beta diversity indices -> using a presence absenece matrix (calculate LCBD, replLCBD, richLCBD, SCBD) (code from Alice Belskis) ----
#First, create a presence absence matrix using macro.wide
macro.wide2 <- macro.wide

##Create presence/absence matrix to use for adespatial package with Jaccard's dissimilarity
macro.pre.ab <- macro.wide %>%
  mutate_if(is.numeric, ~1* (. !=0))
##make Site_ID the row names
as.data.frame(macro.pre.ab)
macro.pre.ab <- as.data.frame(macro.pre.ab)
rownames(macro.pre.ab) <- macro.pre.ab[,1] #Specify SITE_ID as rownames
macro.pre.ab <- macro.pre.ab[,-1] #Remove first column that has the old SITE_ID tags

##LCBD Analysis below 
#Run beta.div for LCBD values
full_LCBD <- beta.div(macro.pre.ab, sqrt.D = FALSE, nperm = 999)
LCBD <- full_LCBD$LCBD #pull out LCBD values into separate df
p.LCBD <- full_LCBD$p.LCBD #pull out p values for LCBD into separate df


lcbd <- cbind(LCBD, p.LCBD) #Make new dataframe
lcbd <- as.data.frame(lcbd) #Examine new df for any significant values
#No significant p values for LCBD (insert sad face)
##This means that sites are not ecologically unique among each other for their contributions to beta diversity
lcbd <- rownames_to_column(lcbd)
names(lcbd)[names(lcbd) == 'rowname'] <- 'SITE_ID'

##Beta replacement
###NO nestedness because my sites are not connected with each other
taxa_betas <- beta.div.comp(macro.pre.ab, coef = "J", quant = TRUE, save.abc = TRUE) #finding replacement and richness values
summary(taxa_betas)
taxa_repl <- taxa_betas$repl #pull out repl into different values
taxa_rich <- taxa_betas$rich #pull out rich into different values

taxa_rich_LCBD <- beta.div(taxa_rich, sqrt.D = FALSE, nperm = 999)
rich_LCBD <- taxa_rich_LCBD$LCBD
rich_p_LCBD <- taxa_rich_LCBD$p.LCBD

richness <- cbind(rich_LCBD, rich_p_LCBD) #Make new richness dataframe
richness <- as.data.frame(richness)  #Examine new df for any sig
#No significant site for richness 

taxa_repl_LCBD <- beta.div(taxa_repl, sqrt.D = FALSE, nperm = 999)
repl_LCBD <- taxa_repl_LCBD$LCBD
repl_p_LCBD <- taxa_repl_LCBD$p.LCBD

replacement <- cbind(repl_LCBD, repl_p_LCBD) #Make new replacement dataframe
replacement <- as.data.frame(replacement) #examine new df for sig
#replacement is not significant among sites (insert sad face)

#####Bar graph of replacement & richness
rich_row <- rownames_to_column(richness)
rep_row <- rownames_to_column(replacement)
repl.rich <- left_join(rich_row,rep_row, by= "rowname")
names(repl.rich)[names(repl.rich) == 'rowname'] <- 'SITE_ID'

##convert into three columns where you have a value and metric column -> did this in excel because it was easier. *See comment below regarding Github!*
#rr.plot <- read.csv("/Users/masonward/Desktop/Research/Spongy Moth/repl.rich. r stuff.csv") #check dataframe - delete three columns
#Added CSV onto my Github repository
rr.plot <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/repl_rich_oct30.csv")


#Calculate Jaccard dissimilarity — repl + rich - Legendre 2014 (supplemental material, pg 3)
repl.rich$jaccard <- repl.rich$rich_LCBD + repl.rich$repl_LCBD 

rr.plot %>% mutate(SITE_ID = fct_relevel(SITE_ID, "SPRL_3", "SPRL_1", "MSHN_3", "MSHN_2","MSHN_1", "BDEG_3", "BDEG_2")) %>%
  ggplot(aes(fill = metric, y = value, x = SITE_ID)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  coord_flip() +
  scale_fill_manual(labels = c("Replacement", "Richness difference"), values=c("#D95F02", "#56B4E9")) +
  ylab("Jaccard dissimilarity coefficient") +
  xlab("Site identification") +
  labs(fill = "") +
  fig.theme 

##SCBD Analysis 
SCBD <- as.data.frame(full_LCBD$SCBD)
SCBD$row_names <- row.names(SCBD)
SCBD <- SCBD[order(SCBD$`full_LCBD$SCBD`, decreasing = TRUE),]
names(SCBD)[names(SCBD) == 'full_LCBD$SCBD'] <- 'SCBD'
names(SCBD)[names(SCBD) == 'row_names'] <- 'Taxa'
rownames(SCBD) <- NULL
SCBD_list <- as.list(SCBD)
print(SCBD_list)
head(SCBD, 10)
##SCBD is useful to determine which species exhibit large variations across the study area.
###The sites where species with large SCBD values are abundant and dominate the community will normally also have large LCBD indices
#SCBD          Taxa
#1  0.04889426    Chauliodes
#2  0.04889426    Sminthurus
#3  0.04589514 Orthocladinae
#4  0.04563819   Limnephilus
#5  0.04475784     Liodessus
#6  0.04475784   Chaoboridae
#7  0.04194786       Lesteva
#8  0.04074522     Scirtidae
#9  0.04045957   Oligochaeta
#10 0.03814399        Gerris



SCBD %>% mutate(Taxa = fct_reorder(Taxa, SCBD)) %>% 
  ggplot(aes(SCBD, Taxa)) +
  geom_point() +
  ylab("Taxa") +
  xlab("SCBD values") +
  fig.theme
#fct_reorder makes the list go ascending to descending 
##It looks fine to me — maybe italize genera names

##Combine all diveristy indices into one dataframe
diversity_df <- left_join(repl.rich, lcbd, by= "SITE_ID")
diversity_df <- left_join(diversity_df, alpha.diversity, by = "SITE_ID")
#Remove p values because I do not need them in this dataframe
diversity_df <- diversity_df[-c(3,5,6,8)]




#Environmental data ----
##Use "plot_histogram" to run histograms on all environmental variables to see if they are normally distributed
histogram_plot <- plot_histogram(environmental_df, binary_as_factor = TRUE, geom_histogram_args = list(bins = 30L),
                                 scale_x = "continuous", title = NULL, ggtheme = theme_gray(), theme_config = list(),
                                 nrow = 4L, ncol = 5L, parallel = FALSE) 
#data is not normally distributed 
##Use methods from "Local contributions to beta diversity in urban pond networks: Implications for biodiversity conservation and management"
###Hill et al 2021 (Had one value per site so average values together; can look at individual weeks later)
avg_environmental <- environmental_df %>%
  select(-c(WEEK_COL, WET_OR_DRY)) %>% #Selecting which column to not include
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
            avg_FRASS_g = mean(FRASS_g, na.rm = TRUE),
            avg_LEAF_LITTER_g = mean(LEAF_LITTER_g, na.rm = TRUE))

#Remove duplicate rows
avg_environmental <- avg_environmental[!duplicated(avg_environmental),]

#Combine macroinvertebrate and environmental dataframes for GLM analysis ----
main_df <- left_join(diversity_df, avg_environmental, by = "SITE_ID") #woohoo, go me
#Transform diversity indices to have a normal distribution for GLM analysis
hist(main_df$Shannon_Diversity) #left-skewed data (sq rt)
hist(main_df$rich_LCBD) #I do not even know what skewed (log)
hist(main_df$repl_LCBD) #I do not know
hist(main_df$LCBD) #The hell is this?



#GLM analysis for shannon diversity metric x environmental variable ----
glm_sh_ph <- glm(Shannon_Diversity ~ avg_PH, data = main_df)
summary(glm_sh_ph) 
glm_sh_conductivity <- glm(Shannon_Diversity ~ avg_CONDUCTIVITY_uS, data = main_df)
summary(glm_sh_conductivity) 
glm_sh_temp <- glm(Shannon_Diversity ~ avg_TEMPERATURE_C, data = main_df)
summary(glm_sh_temp) 
glm_sh_depth <- glm(Shannon_Diversity ~ avg_DEPTH_CM, data = main_df)
summary(glm_sh_depth) 
glm_sh_length <- glm(Shannon_Diversity ~ avg_LENGTH_M, data = main_df)
summary(glm_sh_length) 
glm_sh_width <- glm(Shannon_Diversity ~ avg_WIDTH_M, data = main_df)
summary(glm_sh_width) 
glm_sh_canopy <- glm(Shannon_Diversity ~ avg_CANOPY, data = main_df)
summary(glm_sh_canopy) 
glm_sh_TOC <- glm(Shannon_Diversity ~ avg_TOC_mg_L, data = main_df)
summary(glm_sh_TOC) 
glm_sh_TN <- glm(Shannon_Diversity ~ avg_TN_mg_L, data = main_df)
summary(glm_sh_TN) 
glm_sh_DOC <- glm(Shannon_Diversity ~ avg_DOC_mg_L, data = main_df)
summary(glm_sh_DOC) 
glm_sh_TDN <- glm(Shannon_Diversity ~ avg_TDN_mg_L, data = main_df)
summary(glm_sh_TDN) 
glm_sh_chloride <- glm(Shannon_Diversity ~ avg_CHLORIDE_mg_L, data = main_df)
summary(glm_sh_chloride) 
glm_sh_nitrate <- glm(Shannon_Diversity ~ avg_NITRATE_mg_L, data = main_df)
summary(glm_sh_nitrate) 
glm_sh_sulfate <- glm(Shannon_Diversity ~ avg_SULFATE_mg_L, data = main_df)
summary(glm_sh_sulfate) 
glm_sh_phosphate <- glm(Shannon_Diversity ~ avg_PHOSPHATE_mg_L, data = main_df)
summary(glm_sh_phosphate) 
glm_sh_frass <- glm(Shannon_Diversity ~ avg_FRASS_g, data = main_df)
summary(glm_sh_frass) 
glm_sh_leaf <- glm(Shannon_Diversity ~ avg_LEAF_LITTER_g, data = main_df)
summary(glm_sh_leaf) 

Shannonmodels <- list(glm_sh_ph, glm_sh_conductivity, glm_sh_temp, glm_sh_depth, glm_sh_length, glm_sh_width, 
                      glm_sh_canopy, glm_sh_TOC, glm_sh_TN, glm_sh_DOC, glm_sh_TDN, glm_sh_chloride, glm_sh_nitrate, 
                      glm_sh_sulfate, glm_sh_phosphate, glm_sh_frass, glm_sh_leaf)
model.names <- c("PH", "Conductivity", "Temperature", "Depth", "Length", "Width", 
                      "Canopy", "TOC", "TN", "DOC", "TDN", "Chloride", "Nitrate", 
                      "Sulfate", "Phosphate", "Frass", "Leaf")
ShannonModelAICC <- as.data.frame(aictab(cand.set = Shannonmodels, modnames = model.names))
print(ShannonModelAICC)
#Modnames K     AICc Delta_AICc  ModelLik     AICcWt         LL    Cum.Wt
#14Sulfate 3 11.71022   0.000000 1.0000000 0.28283486  1.1448907 0.2828349
#Sulfate is the best model due to being 2 points below next variable in terms of AICc
pR2(glm_sh_sulfate)['McFadden'] #R2 value 0.8645898 

#GLM analysis for lcbd metric x environmental variable ----
glm_lcbd_ph <- glm(LCBD ~ avg_PH, data = main_df)
summary(glm_lcbd_ph) 
glm_lcbd_conductivity <- glm(LCBD ~ avg_CONDUCTIVITY_uS, data = main_df)
summary(glm_lcbd_conductivity) 
glm_lcbd_temp <- glm(LCBD ~ avg_TEMPERATURE_C, data = main_df)
summary(glm_lcbd_temp) 
glm_lcbd_depth <- glm(LCBD ~ avg_DEPTH_CM, data = main_df)
summary(glm_lcbd_depth) 
glm_lcbd_length <- glm(LCBD ~ avg_LENGTH_M, data = main_df)
summary(glm_lcbd_length) 
glm_lcbd_width <- glm(LCBD ~ avg_WIDTH_M, data = main_df)
summary(glm_lcbd_width) 
glm_lcbd_canopy <- glm(LCBD ~ avg_CANOPY, data = main_df)
summary(glm_lcbd_canopy) 
glm_lcbd_TOC <- glm(LCBD ~ avg_TOC_mg_L, data = main_df)
summary(glm_lcbd_TOC) 
glm_lcbd_TN <- glm(LCBD ~ avg_TN_mg_L, data = main_df)
summary(glm_lcbd_TN) 
glm_lcbd_DOC <- glm(LCBD ~ avg_DOC_mg_L, data = main_df)
summary(glm_lcbd_DOC) 
glm_lcbd_TDN <- glm(LCBD ~ avg_TDN_mg_L, data = main_df)
summary(glm_lcbd_TDN) 
glm_lcbd_chloride <- glm(LCBD ~ avg_CHLORIDE_mg_L, data = main_df)
summary(glm_lcbd_chloride) 
glm_lcbd_nitrate <- glm(LCBD ~ avg_NITRATE_mg_L, data = main_df)
summary(glm_lcbd_nitrate) 
glm_lcbd_sulfate <- glm(LCBD ~ avg_SULFATE_mg_L, data = main_df)
summary(glm_lcbd_sulfate) 
glm_lcbd_phosphate <- glm(LCBD ~ avg_PHOSPHATE_mg_L, data = main_df)
summary(glm_lcbd_phosphate) 
glm_lcbd_frass <- glm(LCBD ~ avg_FRASS_g, data = main_df)
summary(glm_lcbd_frass) 
glm_lcbd_leaf <- glm(LCBD ~ avg_LEAF_LITTER_g, data = main_df)
summary(glm_lcbd_leaf) 

LCBDmodels <- list(glm_lcbd_ph, glm_lcbd_conductivity, glm_lcbd_temp, glm_lcbd_depth, glm_lcbd_length, glm_lcbd_width, 
                      glm_lcbd_canopy, glm_lcbd_TOC, glm_lcbd_TN, glm_lcbd_DOC, glm_lcbd_TDN, glm_lcbd_chloride, glm_lcbd_nitrate, 
                      glm_lcbd_sulfate, glm_lcbd_phosphate, glm_lcbd_frass, glm_lcbd_leaf)
model.names <- c("PH", "Conductivity", "Temperature", "Depth", "Length", "Width", 
                 "Canopy", "TOC", "TN", "DOC", "TDN", "Chloride", "Nitrate", 
                 "Sulfate", "Phosphate", "Frass", "Leaf")
LCBDModelAICC <- as.data.frame(aictab(cand.set = LCBDmodels, modnames = model.names))
print(LCBDModelAICC)
#No significance 

#GLM analysis for repl metric x environmental variable ----
glm_repl_ph <- glm(repl_LCBD ~ avg_PH, data = main_df)
summary(glm_repl_ph) 
glm_repl_conductivity <- glm(repl_LCBD ~ avg_CONDUCTIVITY_uS, data = main_df)
summary(glm_repl_conductivity) 
glm_repl_temp <- glm(repl_LCBD ~ avg_TEMPERATURE_C, data = main_df)
summary(glm_repl_temp) 
glm_repl_depth <- glm(repl_LCBD ~ avg_DEPTH_CM, data = main_df)
summary(glm_repl_depth) 
glm_repl_length <- glm(repl_LCBD ~ avg_LENGTH_M, data = main_df)
summary(glm_repl_length) 
glm_repl_width <- glm(repl_LCBD ~ avg_WIDTH_M, data = main_df)
summary(glm_repl_width) 
glm_repl_canopy <- glm(repl_LCBD ~ avg_CANOPY, data = main_df)
summary(glm_repl_canopy) 
glm_repl_TOC <- glm(repl_LCBD ~ avg_TOC_mg_L, data = main_df)
summary(glm_repl_TOC) 
glm_repl_TN <- glm(repl_LCBD ~ avg_TN_mg_L, data = main_df)
summary(glm_repl_TN) 
glm_repl_DOC <- glm(repl_LCBD ~ avg_DOC_mg_L, data = main_df)
summary(glm_repl_DOC) 
glm_repl_TDN <- glm(repl_LCBD ~ avg_TDN_mg_L, data = main_df)
summary(glm_repl_TDN) 
glm_repl_chloride <- glm(repl_LCBD ~ avg_CHLORIDE_mg_L, data = main_df)
summary(glm_repl_chloride) 
glm_repl_nitrate <- glm(repl_LCBD ~ avg_NITRATE_mg_L, data = main_df)
summary(glm_repl_nitrate) 
glm_repl_sulfate <- glm(repl_LCBD ~ avg_SULFATE_mg_L, data = main_df)
summary(glm_repl_sulfate) 
glm_repl_phosphate <- glm(repl_LCBD ~ avg_PHOSPHATE_mg_L, data = main_df)
summary(glm_repl_phosphate) 
glm_repl_frass <- glm(repl_LCBD ~ avg_FRASS_g, data = main_df)
summary(glm_repl_frass) 
glm_repl_leaf <- glm(repl_LCBD ~ avg_LEAF_LITTER_g, data = main_df)
summary(glm_repl_leaf) 

REPLmodels <- list(glm_repl_ph, glm_repl_conductivity, glm_repl_temp, glm_repl_depth, glm_repl_length, glm_repl_width, 
                   glm_repl_canopy, glm_repl_TOC, glm_repl_TN, glm_repl_DOC, glm_repl_TDN, glm_repl_chloride, glm_repl_nitrate, 
                   glm_repl_sulfate, glm_repl_phosphate, glm_repl_frass, glm_repl_leaf)
model.names <- c("PH", "Conductivity", "Temperature", "Depth", "Length", "Width", 
                 "Canopy", "TOC", "TN", "DOC", "TDN", "Chloride", "Nitrate", 
                 "Sulfate", "Phosphate", "Frass", "Leaf")
ReplModelAICC <- as.data.frame(aictab(cand.set = REPLmodels, modnames = model.names))
print(ReplModelAICC)
#Modnames K       AICc Delta_AICc     ModelLik       AICcWt       LL    Cum.Wt
#1            PH 3 -26.02739   0.000000 1.0000000000 6.403231e-01 20.01369 0.6403231
#PH is the best model due to being 2 points below next variable in terms of AICc 
pR2(glm_repl_ph)['McFadden'] #R2 value is -0.3155121

#GLM analysis for rich metric x environmental variable ----
glm_rich_ph <- glm(rich_LCBD ~ avg_PH, data = main_df)
summary(glm_rich_ph) 
glm_rich_conductivity <- glm(rich_LCBD ~ avg_CONDUCTIVITY_uS, data = main_df)
summary(glm_rich_conductivity) 
glm_rich_temp <- glm(rich_LCBD ~ avg_TEMPERATURE_C, data = main_df)
summary(glm_rich_temp) 
glm_rich_depth <- glm(rich_LCBD ~ avg_DEPTH_CM, data = main_df)
summary(glm_rich_depth) 
glm_rich_length <- glm(rich_LCBD ~ avg_LENGTH_M, data = main_df)
summary(glm_rich_length) 
glm_rich_width <- glm(rich_LCBD ~ avg_WIDTH_M, data = main_df)
summary(glm_rich_width) 
glm_rich_canopy <- glm(rich_LCBD ~ avg_CANOPY, data = main_df)
summary(glm_rich_canopy) 
glm_rich_TOC <- glm(rich_LCBD ~ avg_TOC_mg_L, data = main_df)
summary(glm_rich_TOC) 
glm_rich_TN <- glm(rich_LCBD ~ avg_TN_mg_L, data = main_df)
summary(glm_rich_TN) 
glm_rich_DOC <- glm(rich_LCBD ~ avg_DOC_mg_L, data = main_df)
summary(glm_rich_DOC) 
glm_rich_TDN <- glm(rich_LCBD ~ avg_TDN_mg_L, data = main_df)
summary(glm_rich_TDN) 
glm_rich_chloride <- glm(rich_LCBD ~ avg_CHLORIDE_mg_L, data = main_df)
summary(glm_rich_chloride) 
glm_rich_nitrate <- glm(rich_LCBD ~ avg_NITRATE_mg_L, data = main_df)
summary(glm_rich_nitrate) 
glm_rich_sulfate <- glm(rich_LCBD ~ avg_SULFATE_mg_L, data = main_df)
summary(glm_rich_sulfate) 
glm_rich_phosphate <- glm(rich_LCBD ~ avg_PHOSPHATE_mg_L, data = main_df)
summary(glm_rich_phosphate) 
glm_rich_frass <- glm(rich_LCBD ~ avg_FRASS_g, data = main_df)
summary(glm_rich_frass) 
glm_rich_leaf <- glm(rich_LCBD ~ avg_LEAF_LITTER_g, data = main_df)
summary(glm_rich_leaf)

RICHmodels <- list(glm_rich_ph, glm_rich_conductivity, glm_rich_temp, glm_rich_depth, glm_rich_length, glm_rich_width, 
                   glm_rich_canopy, glm_rich_TOC, glm_rich_TN, glm_rich_DOC, glm_rich_TDN, glm_rich_chloride, glm_rich_nitrate, 
                   glm_rich_sulfate, glm_rich_phosphate, glm_rich_frass, glm_rich_leaf)
model.names <- c("PH", "Conductivity", "Temperature", "Depth", "Length", "Width", 
                 "Canopy", "TOC", "TN", "DOC", "TDN", "Chloride", "Nitrate", 
                 "Sulfate", "Phosphate", "Frass", "Leaf")
RichModelAICC <- as.data.frame(aictab(cand.set = RICHmodels, modnames = model.names))
print(RichModelAICC)
#Nothing is significant












#LCBD Maps of Pond sites ----
##Merge LCBD values with pond points
###Start here for shapefile November 1st
lcbd_map <- inner_join(lcbd, vernal.map_df, by = "SITE_ID")
lcbd_map <- st_as_sf(lcbd_map, coords = c("long", "lat"), crs = 4269)
repl_rich_map <- inner_join(repl.rich, vernal.map_df, by = "SITE_ID")
repl_rich_map <- st_as_sf(repl_rich_map, coords = c("long", "lat"), crs = 4269)


LCBDMap <- ggplot() + 
  geom_sf(data = shp_forest_subset, aes(fill = SF_Name)) + 
  geom_sf(data = shp_boundaries_subset, aes(fill = DistrictNa)) + 
  geom_sf(data = lcbd_map, aes(color = LCBD, size = 2)) + 
  scale_fill_manual(values = c("#E5F5F9", "#D4B9DA", "#D9D9D9")) + 
  ggtitle("Map of vernal pond study sites") +
  theme_bw() + 
  ylab(expression("Latitude ("*degree*")" )) + 
  xlab(expression("Longitude ("*degree*")" )) +
  guides(fill=guide_legend(title="State Forests")) +
  coord_sf() +
  annotation_scale() +
  fig.theme

REPLMap <- ggplot() + 
  geom_sf(data = shp_forest_subset, aes(fill = SF_Name)) + 
  geom_sf(data = shp_boundaries_subset, aes(fill = DistrictNa)) + 
  geom_sf(data = repl_rich_map, aes(color = repl_LCBD, size = 2)) + 
  scale_fill_manual(values = c("#E5F5F9", "#D4B9DA", "#D9D9D9")) + 
  ggtitle("Map of vernal pond study sites") +
  theme_bw() + 
  ylab(expression("Latitude ("*degree*")" )) + 
  xlab(expression("Longitude ("*degree*")" )) +
  guides(fill=guide_legend(title="State Forests")) +
  coord_sf() +
  annotation_scale() +
  fig.theme

RICHMap <- ggplot() + 
  geom_sf(data = shp_forest_subset, aes(fill = SF_Name)) + 
  geom_sf(data = shp_boundaries_subset, aes(fill = DistrictNa)) + 
  geom_sf(data = repl_rich_map, aes(color = rich_LCBD, size = 2)) + 
  scale_fill_manual(name="LCBD Richness",values = c("#E5F5F9", "#D4B9DA", "#D9D9D9")) + 
  ggtitle("Map of vernal pond study sites") +
  theme_bw() + 
  ylab(expression("Latitude ("*degree*")" )) + 
  xlab(expression("Longitude ("*degree*")" )) +
  guides(fill=guide_legend(title="State Forests")) +
  coord_sf() +
  annotation_scale() +
  fig.theme


