#LIBRARIES

library(stringr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(wesanderson)
library(gridExtra)
library(ggpubr)

##EDITING VCF TABLE

#Set path
setwd("/home/mfcg/Descargas/covid/omicron/ver2gendiv")

#Read the identified variants mutations with frequencies and proportions and filtered metadata
fqnmeta <- read.csv(file = "newfqnmeta.csv",
                    header = TRUE)
#Metadata filtered alone
newfilmet <- read.delim(file = "newfilmet.tsv",
                        header = TRUE)
#Palettes
dj2pal <- wes_palette("Darjeeling2", 5, type = "discrete")
dj1pal <- wes_palette("Darjeeling1", 10, type = "continuous")


#PROVINCE
#Frequencies per province
table(newfilmet$Location)
#Frequencies per province and month
newfilmet %>% count(Location, Month)

#Adding variable that indicates province and its month
fqnmeta <- fqnmeta %>%
  mutate(promon = paste0(LOCAT, MONTH))

#Frequencies by province and month
prmnfq <- fqnmeta %>%
  group_by(POSIC, ALT1, LOCAT, MONTH, promon) %>% #Group by genomic position, alternative allele, province and month
  summarise(n=n()) %>%
  mutate(pmfq = case_when(
    promon == "Eastern CapeSeptember" ~ n/46,
    promon == "Eastern CapeOctober" ~ n/39,
    promon == "GautengSeptember" ~ n/47,
    promon == "GautengOctober" ~ n/25,
    promon == "GautengNovember" ~ n/70,
    promon == "KwaZulu-NatalSeptember" ~ n/12,
    promon == "KwaZulu-NatalOctober" ~ n/102,
    promon == "KwaZulu-NatalNovember" ~ n/102,
    promon == "MpumalangaSeptember" ~ n/5
  ))

ggplot(data = prmnfq, aes(x = pmfq, y = LOCAT, fill = MONTH)) + 
  geom_density_ridges(alpha = 0.6) +
  scale_fill_manual(values = dj2pal[c(2, 3, 4)]) +
  labs(x = "Polymorphic site frequency",
       y = "Locations") +
  theme_minimal()

ggplot(data = prmnfq, aes(x = pmfq, fill = MONTH)) + 
  geom_boxplot(varwidth = TRUE) +
  scale_fill_manual(values = dj1pal[1:9]) +
  facet_wrap(vars(LOCAT), nrow = 4) +
  theme_minimal()

ggplot(data = prmnfq, aes(x = pmfq, y = MONTH, fill = MONTH)) +
         geom_violin(alpha = 0.5, lwd = 0) +
         geom_boxplot(width = 0.05, lwd = 0.5, outlier.colour = NA) +
         scale_fill_manual(values = dj2pal[2:4]) +
         facet_wrap(vars(LOCAT), nrow = 4) +
  labs(x = "Genome-wide Alternative Allele Proportion",
       y = "") +
         theme_minimal() +
  theme(legend.position = "none")
