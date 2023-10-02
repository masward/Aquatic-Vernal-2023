#Mason S Ward
#2023-09-07
#Vernal pond dataset
#The Pennsylvania State University

#Load applicable libraries
library(dplyr)
library(tibble)
library(vegan)
library(tidyverse)
library(reshape2)
library(viridis)
library(ggrepel)

#bring in dataset from Github - masward (username)
vernal_df <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/Vernal_Pond_Datasheet_REVISED.csv")
##One row with no data - remove it
vernal_df <- vernal_df[-c(79),]

##Plotting data below - Canopy Cover is first
ggplot(vernal_df, aes(x = WEEK_COL, y = PCT_CANOPY_COVER, group = SITE_ID)) +
  geom_line(aes(colour = SITE_ID)) +
  geom_point() +
  scale_fill_viridis_d()

#Average all environmental variable values across each site
install.packages("tidyr")
library(tidyr)

colnames(vernal_df)

vernal_df$PH <- as.numeric(vernal_df$PH)

test <- vernal_df %>% 
  select(-WEEK_COL, -LATITUDE, -LONGITUDE, -WET_OR_DRY) %>% 
  group_by(SITE_ID) %>% 
  nest() 
# %>% 
#   summarise(avg_ph = mean(PH, na.rm = TRUE)) # this did not work 
as_tibble(test)

test <- vernal_df %>% 
  select(SITE_ID, PH) %>% 
  group_by(SITE_ID) %>% 
  mutate(test_mean = mean(PH)) %>% 
  print(n = Inf)

as_tibble(test)

test_2 <- vernal_df %>% 
  select(SITE_ID, PH)
as_tibble(test_2)

mean(test_2[test_2$SITE_ID == 'SPRL_1', 'PH'], na.rm = TRUE)
