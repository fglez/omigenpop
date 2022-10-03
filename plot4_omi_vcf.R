#LIBRARIES

library(stringr)
library(tidyverse)
library(ggplot2)
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
dj1pal <- wes_palette("Darjeeling1", 11, type = "continuous")


#A function that shows the omicron updated mutations present in a dataframe
oldomi <- c("Q493R", "C25000T", "C25584T", "RG203KR", NA) #No longer part of B.1.1.529 lineage defining mutations
length(oldomi)
updomi <- function(fqndf){ #A dataframe with OMMU column that indicates omicron mutations
  uniom <- unique(fqndf$OMMU) #Unique omicron mutations added by hand
  uniom <- uniom[!uniom %in% oldomi] #Removing those that are no longer common to all B.1.1.529 lineage
  noomi <- length(uniom) #How many
  rsuom <- list(uniom, noomi)
  return(rsuom)
}


#MONTH
#Frequencies per month
table(newfilmet$Month) #Samples per month

#November
novfq <- fqnmeta %>%
  filter(MONTH == "November") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(timfq = n/172)
noval <- fqnmeta %>%
  filter(MONTH == "November") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD))#MEdian Locus Depth
novem <- inner_join(novfq, noval, by = c("POSIC", "ALT1"))
updomi(noval) #B.1.1.529 lineage defining mutations in November

#October
octfq <- fqnmeta %>%
  filter(MONTH == "October") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(timfq = n/166)
octal <- fqnmeta %>%
  filter(MONTH == "October") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD))#MEdian Locus Depth
octob <- inner_join(octfq, octal, by = c("POSIC", "ALT1"))
updomi(octal) #B.1.1.529 lineage defining mutations in October


#September
sepfq <- fqnmeta %>%
  filter(MONTH == "September") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(timfq = n/110)
sepal <- fqnmeta %>%
  filter(MONTH == "September") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD))#MEdian Locus Depth
septe <- inner_join(sepfq, sepal, by = c("POSIC", "ALT1"))
updomi(sepal) #B.1.1.529 lineage defining mutations in September


#All the data by months
meses <- rbind(septe, octob, novem)

#Violin plot with AAFs per month
mes_plotlin <- ggplot(data = meses,
                      aes(x = timfq, y = MONTH, fill = MONTH)) +
  geom_violin(alpha = 0.5, lwd = 0) +
  geom_boxplot(width = 0.05, lwd = 0.5, outlier.colour = NA)+
  scale_fill_manual(values = dj2pal[c(2, 3, 4)]) +
  labs(x = "Genome-wide Alternative Allele Proportion",
       y = "Month") +
  theme_minimal() +
  theme(text = element_text(size = 9))

#PROVINCE
#Frequencies per province
table(newfilmet$Location)

#Gauteng
gaufq <- fqnmeta %>%
  filter(LOCAT == "Gauteng") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(locfq = n/142)
gaual <- fqnmeta %>%
  filter(LOCAT ==  "Gauteng") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD))#MEdian Locus Depth
gaute <- inner_join(gaufq, gaual, by = c("POSIC", "ALT1"))
updomi(gaual) #Omicron mutations present in Gauteng

#KwaZulu-Natal
kwafq <- fqnmeta %>%
  filter(LOCAT == "KwaZulu-Natal") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(locfq = n/216)
kwaal <- fqnmeta %>%
  filter(LOCAT == "KwaZulu-Natal") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD))#MEdian Locus Depth
kwana <- inner_join(kwafq, kwaal, by = c("POSIC", "ALT1"))
updomi(kwaal) #Omicron mutations present in KwaZulu-Natal

#Eastern Cape
easfq <- fqnmeta %>%
  filter(LOCAT == "Eastern Cape") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(locfq = n/85)
easal <- fqnmeta %>%
  filter(LOCAT == "Eastern Cape") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD))#MEdian Locus Depth
easca <- inner_join(easfq, easal, by = c("POSIC", "ALT1"))
updomi(easal) #Omicron mutations present in Eastern Cape

#All provinces data together
provi <- rbind(easca, gaute, kwana)

#Violin plot with AAFs per month
loc_plotlin <- ggplot(data = provi,
                      aes(x = locfq, y = LOCAT, fill = LOCAT)) +
  geom_violin(alpha = 0.5, lwd = 0) +
  geom_boxplot(width = 0.05, lwd = 0.5, outlier.colour = NA)+
  scale_fill_manual(values = dj2pal[c(2, 3, 4)]) +
  labs(x = "Genome-wide Alternative Allele Proportion",
       y = "Location") +
  theme_minimal() +
  theme(text = element_text(size = 9))

#Full plot
ggarrange(mes_plotlin, loc_plotlin, nrow = 2, heights = c(1, 1), align = "v",
          labels = c("A", "B"), font.label = list(size = 9))


#Kolmogorov Smirnov test for months CDF
ks.test(x = novem$timfq, y = octob$timfq, B = 9999) #November vs October GW PSDs distribution
ks.test(x = novem$timfq, y = septe$timfq, B = 9999) #November vs September GW PSDs distribution
ks.test(x = septe$timfq, y = octob$timfq, B = 9999) #September vs October GW PSDs distribution

#Kolmogorov Smirnov test for locations CDF
ks.test(x = easca$locfq, y = gaute$locfq, B = 9999) #Eastern Cape vs Gauteng GW PSDSs distribution
ks.test(x = easca$locfq, y = kwana$locfq, B = 9999) #Eastern Cape vs KwaZulu-Natal GW PSDSs distribution
ks.test(x = kwana$locfq, y = gaute$locfq, B = 9999) #KwaZulu-Natal vs Gauteng GW PSDSs distribution
