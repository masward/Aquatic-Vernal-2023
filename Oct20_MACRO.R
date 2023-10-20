#Mason S Ward
#2023-10-04
#Macroinvertebrate vernal pond dataset
#The Pennsylvania State University

#Load applicable libraries
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


#bring in dataset from Github - masward (username) under the "Spongy-Moth" repository
MACRO.DF <- read.csv("https://raw.githubusercontent.com/masward/Spongy-Moth/main/Macroinvertebrate_Datasheet_R%20copy.csv")
#check dataframe structure - looks good

# ALPHA DIVERSITY INDICES - start with main dataframe for alpha diveristy since it has abundance ----
## Sum by SITE_ID and like taxa
MACRO.ABUND <- MACRO.DF %>%
  group_by(SITE_ID, TARGET_TAXON) %>% #Group by SITE_ID & TARGET_TAXON
  summarize(TAXA_ABUND = sum(TOTAL)) #Sum values together and create a TAXA_ABUND column

#Convert data to wide format for vegan analyses (vegan likes wide format for diversity indices)
MACRO.WIDE <- MACRO.ABUND %>%
  spread(key = "TARGET_TAXON",     # key = column to go on top
         value = "TAXA_ABUND") %>%     # value = column to spread under each new column
  mutate(across(everything(), ~replace_na(., 0)))     # replaces all the NAs with 0s

#Alpha diversity below using the function "diversity" in the vegan package for Shannon Div. 
ALPHA.DIV <- MACRO.WIDE #create separate dataframe for alpha diverisity
ShannonDiversity <- diversity(MACRO.WIDE[2:39]) #calculated diversity across week and site; only select rows 3:40 because the first two columns are being calculated
ALPHA.DIV$ShannonDiversity <- ShannonDiversity #add Shannon values to the new dataframe for ease of access
##Remove all columns in ALPHA.DIV except for ShannonDiversity
ALPHA.DIV <- subset(ALPHA.DIV, select = -c(2:39))
###Plotting Shannon Diversity
ggplot(ALPHA.DIV, aes(x= SITE_ID, y = ShannonDiversity)) +
  geom_bar(stat= "identity",linewidth = 2) +
  ylab("Shannon's H'") + 
  xlab("Site ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
#Check for normality
hist(ALPHA.DIV$ShannonDiversity) #Not normal, run KW test
 
kruskal_result <- kruskal.test(ShannonDiversity ~ SITE_ID, ALPHA.DIV)
kruskal_result #Kruskal-Wallis chi-squared = 6, df = 6, p-value = 0.4232 (Not statistically different)






# BETA DIVERSITY INDICES - subset MACRO.ABUND dataframe to form presence/absence matrix ----
MACRO.WIDE2 <- MACRO.WIDE

##Create presence/absence matrix ----
MACRO.P.A <- MACRO.WIDE2 %>% 
  mutate_if(is.numeric, ~1* (. !=0))

##LCBD Analysis ----
#make Site_ID the row names
as.data.frame(MACRO.P.A)
MACRO.P.A <- as.data.frame(MACRO.P.A)
rownames(MACRO.P.A) <- MACRO.P.A[,1] #Specify SITE_ID as rownames
MACRO.P.A <- MACRO.P.A[,-1] #Remove first column 
#Run beta.div for LCBD values
full_LCBD <- beta.div(MACRO.P.A, sqrt.D = FALSE, nperm = 999)
LCBD <- full_LCBD$LCBD #pull out LCBD values into separate df
p.LCBD <- full_LCBD$p.LCBD #pull out p values for LCBD into separate df
qqplot(LCBD) #did not work ?

lcbd <- cbind(LCBD, p.LCBD) #Make new dataframe
lcbd <- as.data.frame(lcbd)
#No significant p values for LCBD (insert sad face)
##This means that sites are not ecologically unique among each other for their contributions to beta diversity

##Beta replacement ----
###NO nestedness because my sites are not connected with each other
taxa_betas <- beta.div.comp(MACRO.P.A, coef = "J", quant = TRUE, save.abc = TRUE) #finding replacement and richness values
summary(taxa_betas)
taxa_repl <- taxa_betas$repl #pull out repl into different values
taxa_rich <- taxa_betas$rich #pull out rich into different values

taxa_rich_LCBD <- beta.div(taxa_rich, sqrt.D = FALSE, nperm = 999)
rich_LCBD <- taxa_rich_LCBD$LCBD
rich_p_LCBD <- taxa_rich_LCBD$p.LCBD

richness <- cbind(rich_LCBD, rich_p_LCBD)
richness <- as.data.frame(richness) 
#One significant site for richness - MSHN 1 p value 0.037

taxa_repl_LCBD <- beta.div(taxa_repl, sqrt.D = FALSE, nperm = 999)
repl_LCBD <- taxa_repl_LCBD$LCBD
repl_p_LCBD <- taxa_repl_LCBD$p.LCBD

replacement <- cbind(repl_LCBD, repl_p_LCBD)
replacement <- as.data.frame(replacement)
#replacement is not significant (insert sad face)

##SCBD Analysis ----
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
#1  0.04591385     Chauliodes
#2  0.04591385     Sminthurus
#3  0.04171497 Limnephiloidea
#4  0.04150701      Liodessus
#5  0.04054906    Chaoboridae
#6  0.04014028        Lesteva
#7  0.03988088  Orthocladinae
#8  0.03944253      Scirtidae
#9  0.03756511    Oligochaeta
#10 0.03554032       Enochrus

##Create one dataframe for LCBD values
