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
fqnmeta <- read.csv(file = "bnewfqnmeta.csv",
                    header = TRUE)
#Metadata filtered alone
newfilmet <- read.delim(file = "newfilmet.tsv",
                        header = TRUE)
#Palettes

dj2pal <- wes_palette("Darjeeling2", 5, type = "discrete")
dj1pal <- wes_palette("Darjeeling1", 11, type = "continuous")

#Genetic regions SARS-CoV-2
sc2gg <- read.csv(file = "SC2-gene-posic.csv")
geneplo <- ggplot() +
  xlim(0, 30000) +
  geom_rect(data = sc2gg, #Displays in the background with colors the genes of SARS-CoV-2 by position
            aes(xmin = nuc_init, xmax = nuc_end, ymin= 0, ymax = 0.1, fill = Gen),
            alpha = 0.8,
            show.legend = FALSE) +
  scale_fill_manual(values = dj1pal[1:11]) +
  geom_text(data = sc2gg, 
            aes(x=c((nuc_init+nuc_end)/2), y=0.05, label=Gen),
            size = 2,
            angle = 90) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 3))



#A function that shows the omicron updated mutations present in a dataframe
oldomi <- c("Q493R", "D1146D", "T64T", "RG203KR") #No longer part of B.1.1.529 lineage defining mutations
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

#Months lollipop plot of genomic changes 

mupfift <- which(meses$timfq > 0.5)

mont_geno_lpplot <- ggplot() +
  ylim(0, 1) +
  xlim(0, 29903) +
  geom_point(data = meses,
             aes(x = POSIC, y = timfq, color = MEALP), #The alternative allele proportion is displayed by color code
             size = 4, 
             alpha = 1/14.9) +
  geom_segment(data = meses,
               aes(x = POSIC, xend = POSIC, y = 0, yend = timfq, color = MEALP)) +
  scale_colour_gradient(low = dj2pal[3], 
                        high = dj2pal[2],
                        limits = c(0.3, 1),
                        name = "Median alternative allele proportion")+
  geom_text(data = meses[mupfift, ],
            aes(x = POSIC, y = timfq, label = OMMU), #Display alternative alleles corresponding to omicron
            colour = "black",
            nudge_y = 0.02,
            size = 2,
            #angle = 30,
            check_overlap = TRUE) +
  geom_text(data = meses[mupfift, ],
            aes(x = POSIC, y = timfq, label = DEMMU), #Display alternative alleles corresponding to delta
            colour = "blue",
            nudge_y = 0.08,
            size = 2,
            #angle = 30,
            check_overlap = TRUE) +
  geom_text(data = meses[mupfift, ],
            aes(x = POSIC, y = timfq, label = ALMU), #Display alternative alleles corresponding to alpha
            colour = "green",
            nudge_y = 0.14,
            size = 2,
            #angle = 30,
            check_overlap = TRUE) +
  labs(x = "Genomic position", 
       y = "Polymorphic site frequency") +
  theme_minimal() +
  theme(text = element_text(size = 8))

mont_geno_lpplot <- mont_geno_lpplot +
  facet_wrap(vars(MONTH), nrow = 3) +
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.title = element_text(size = 6))

#Figure 2.  Mutations in the SARS-CoV-2 reads by month in four South African provinces.
ggarrange(mont_geno_lpplot, geneplo, ncol = 1, nrow = 2, heights = c(9, 1), align = "v")

#Which omicron mutations are present in >=50% of the samples in each month
updomi(septe[which(septe$timfq > 0.5), ]) #September
updomi(octob[which(octob$timfq > 0.5), ]) #October
updomi(novem[which(novem$timfq > 0.5), ]) #November

#Delta mutations present in each month
#September
unique(sepal$DEMMU)
length(unique(sepal$DEMMU))
unique(septe[which(septe$timfq > 0.5), "DEMMU"]) #>=50% of the samples
#October
unique(octal$DEMMU)
length(unique(octal$DEMMU))
unique(octob[which(octob$timfq > 0.5), "DEMMU"]) #>=50% of the samples
#November
unique(noval$DEMMU)
length(unique(noval$DEMMU))
unique(novem[which(novem$timfq > 0.5), "DEMMU"]) #>=50% of the samples

#ALpha mutations present in each month
#September
unique(sepal$ALMU)
length(unique(sepal$ALMU))
unique(septe[which(septe$timfq > 0.5), "ALMU"]) #>=50% of the samples
#October
unique(octal$ALMU)
length(unique(octal$ALMU))
unique(octob[which(octob$timfq > 0.5), "ALMU"]) #>=50% of the samples
#November
unique(noval$ALMU)
length(unique(noval$ALMU))
unique(novem[which(novem$timfq > 0.5), "ALMU"]) #>=50% of the samples

#Substitutions per site per year
(28-8)/29903/(3/12)

#No labels
mont_geno_lpplot <- ggplot() +
  ylim(0, 1) +
  xlim(0, 29903) +
  geom_point(data = meses,
             aes(x = POSIC, y = timfq, color = MEALP), #The alternative allele proportion is displayed by color code
             size = 4, 
             alpha = 1/14.9) +
  geom_segment(data = meses,
               aes(x = POSIC, xend = POSIC, y = 0, yend = timfq, color = MEALP)) +
  scale_colour_gradient(low = dj2pal[3], 
                        high = dj2pal[2],
                        limits = c(0.3, 1),
                        name = "Median alternative allele proportion")+
  labs(x = "Genomic position", 
       y = "Polymorphic site frequency") +
  theme_minimal() +
  theme(text = element_text(size = 8))

mont_geno_lpplot <- mont_geno_lpplot +
  facet_wrap(vars(MONTH), nrow = 3) +
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.title = element_text(size = 6))
