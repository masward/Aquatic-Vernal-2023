#Mason S Ward
#2023-08-25
#Vernal pond macroinvertebrate dataset
#The Pennsylvania State University

#Load applicable libraries
library(dplyr)
library(tibble)
library(vegan)
library(tidyverse)
library(reshape2)
library(viridis)

#bring in dataset from Github repository - masward (username)
macro_df <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/Macroinvertebrate%20Datasheet%20-%20Revised.csv")

##Convert all of the dates to weeks using case_when function - allows for better accessibility and you can pull weeks faster
macro_week <- macro_df %>% 
  select(DATE_COL, SITE_ID, TARGET_TAXON, TOTAL, LIFE_STAGE) %>% # selecting the columns we want to keep in new dataframe
  mutate(
    DATE_COL = case_when( # case_when only works within mutate 
      DATE_COL == "2023-05-16" ~ "Week_02", # parameter 1 
      DATE_COL == "2023-05-17" ~ "Week_02",
      DATE_COL == "2023-05-22" ~ "Week_03", #parameter 2
      DATE_COL == "2023-05-23" ~ "Week_03",
      DATE_COL == "2023-05-24" ~ "Week_03", 
      DATE_COL == "2023-07-03" ~ "Week_09", #parameter 3
      DATE_COL == "2023-07-05" ~ "Week_09", 
      DATE_COL == "2023-07-06" ~ "Week_09",
      DATE_COL == "2023-07-10" ~ "Week_10", #Hester-Dendy samplers
      DATE_COL == "2023-07-11" ~ "Week_10",
      DATE_COL == "2023-07-12" ~ "Week_10",
      .default = as.character(DATE_COL) #Do not change anything if parameters are not met
    ) #%>% 
    #summarise()
  )
as_tibble(macro_week)
unique(macro_week$DATE_COL) #DATE_COL is in a row? delete it

macro_week2 <- macro_week[-c(389),] #deleting this row
unique(macro_week2$DATE_COL) #row is deleted

#rename column name from DATE_COL to WEEK_Col
colnames(macro_week2)[1] = "WEEK_COL"
print(macro_week2)
#change TOTAL from chr to num in order to summarize (must be numeric and not a character)
macro_week2$TOTAL <- as.numeric(macro_week2$TOTAL)

# sum abundance by site, taxa, and week
macro_abund <- macro_week2 %>%
  group_by(WEEK_COL, SITE_ID, TARGET_TAXON) %>% 
  summarize(TAXA_SUM = sum(TOTAL))

as_tibble(macro_abund)
#Remove the fist row
macro_abund <- macro_abund[-c(1),]

#Convert data to wide format for vegan analyses (vegan likes wide format for diversity indices)
taxa.wide <- macro_abund %>%
  spread(key = "TARGET_TAXON",     # key = column to go on top
         value = "TAXA_SUM") %>%     # value = column to spread under each new column
  mutate(across(everything(), ~replace_na(., 0)))     # replaces all the NAs with 0s

#two rows are taken from the same locality and same date; sum the two together (misspelled them)
## @ stat and coding group, kpg 8/26/2023 <- Katie Gundermann's code
taxa.wide[which(taxa.wide$WEEK_COL == "Week_02" & taxa.wide$SITE_ID == "MSHN_1"), c(3:ncol(taxa.wide))] <-as.list(colSums(
  taxa.wide[which(taxa.wide$WEEK_COL == "Week_02" & taxa.wide$SITE_ID %in% c("MSHN_1", "MHSN_1")),
                            c(3:ncol(taxa.wide))])) 

##delete the column with an error
taxa.wide <- taxa.wide[-which(taxa.wide$WEEK_COL == "Week_02" & taxa.wide$SITE_ID == "MHSN_1"),]
##Create main df for ease of access
MacroTaxa <- taxa.wide

##DIVERSITY INDICES BELOW USING VEGAN - specnumber (richness) & diversity (Shannon)

#Lets do some richness/diversity stuff
SpeciesRichness <- specnumber(MacroTaxa[3:40]) #calculated species richness across week and site; only select rows 3:40 because the first two columns are being calculated
MacroTaxa2 <- MacroTaxa # created new dataframe that includes richness values and not other data; do not want to confused R & better for plotting data
MacroTaxa2$SpeciesRichness <- SpeciesRichness #yay

ShannonDiversity <- diversity(MacroTaxa[3:40]) #calculated diversity across week and site; only select rows 3:40 because the first two columns are being calculated
MacroTaxa2$ShannonDiversity <- ShannonDiversity #add Shannon values to the new dataframe for ease of access

#Create new dataframe with only Week, Site, Richness, and Diversity and no other taxa. Taxa matrix will be used later on for ordinational plotting but not for boxplots
VernalPondMacro <- MacroTaxa2 %>% 
  select(WEEK_COL, SITE_ID, SpeciesRichness, ShannonDiversity)

#HESTER-DENDY should be in a new dataframe by itself
#remove data from WEEK_10 since it is HESTER-DENDY samples (inappropriate to compare among other weeks)
HESTER_DENDY <- VernalPondMacro[c(17:22),]
VernalPondMacro <- VernalPondMacro[-c(17:22),] #use code to remove Week_10 -> can add in later or create separate plot
unique(VernalPondMacro$WEEK_COL) #[1] "Week_02" "Week_03" "Week_09"



##PLOTTING DATA USING GENERALIZED GGPLOT CODE - YAY
##Plot data - suggested comment from lab meeting: Maybe do not use WEEK_COL as X-value

ggplot(data = VernalPondMacro, aes(x = WEEK_COL, y = ShannonDiversity, group = SITE_ID, fill = SITE_ID)) +
  geom_boxplot() +
  scale_fill_brewer() +
  theme_minimal() + 
  labs(x = "Week Collected", y = "Average Shannon Diversity (H)", title = "Average Shannon Diversity across state forests")
##Boxplot needs major help - maybe the x-value needs changed like above comment suggested?








#CODE BELOW MAY NOT BE USEFUL - ONLY USEFUL IF LOOKING ACROSS STATE FORESTS AND NOT ON A POND-BY-POND BASIS

#Think about combining all of the data - removing pond name and just look at richness/shannon across state forests
state_forest <- VernalPondMacro %>% #CHANGE ALL THE NAMES
  select(WEEK_COL, SITE_ID, SpeciesRichness, ShannonDiversity) %>% # selecting the columns we want 
  mutate(
    SITE_ID = case_when( # case_when only works within mutate 
      SITE_ID == "BDEG_1" ~ "BALD_EAGLE",
      SITE_ID == "BDEG_2" ~ "BALD_EAGLE",
      SITE_ID == "BDEG_3" ~ "BALD_EAGLE",
      SITE_ID == "MSHN_1" ~ "MOSHANNON", 
      SITE_ID == "MSHN_2" ~ "MOSHANNON",
      SITE_ID == "MSHN_3" ~ "MOSHANNON", 
      SITE_ID == "SPRL_1" ~ "SPROUL",
      SITE_ID == "SPRL_3" ~ "SPROUL",
      .default = as.character(SITE_ID)
      ))

##Do not use line of code below
##Do not use line of code below
#sproul does not have that much entries so let's remove her
#state_forest <- state_forest[-which(state_forest$SITE_ID == "SPROUL"),] ##ONLY USE IF REMOVING SPROUL - will keep in code just for future reference 09/08/2023

#Maybe do not average and leave values as is?
#Boxplot is not working so went with bar graph below since the average computes one value
#Code below is suggestive, maybe not use and remove averages among state forests

#use the group_by to average among sites (richness, diversity) - maybe do not average among the sites? talk to Jon about this 09/08/2023
AvgVernalPond <- state_forest %>%
  group_by(WEEK_COL, SITE_ID) %>% 
  summarize(Avg_SpeciesRichness = mean(SpeciesRichness),
            Avg_ShannonDiversity = mean(ShannonDiversity))

#It worked so let's plot some data
ggplot(data = AvgVernalPond, aes(x = WEEK_COL, y = Avg_ShannonDiversity, group = SITE_ID, fill = SITE_ID)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer() +
  theme_minimal() + 
  labs(x = "Week Collected", y = "Average Shannon Diversity (H)", title = "Average Shannon Diversity across state forests")
##Bar plot width on WEEK_2 not looking great - why is the width not changing to account for 0 value of Sproul?


