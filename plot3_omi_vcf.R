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
provi <- rbind(gaute, kwana)

#Location lollipop plot of genomic changes

lupfift <- which(provi$locfq > 0.5)

loc_geno_lpplot <- ggplot() +
  ylim(0, 1) +
  xlim(0, 30000) +
  geom_point(data = provi,
             aes(x = POSIC, y = locfq, color = MEALP), #The alternative allele proportion is displayed by color code
             size = 4, 
             alpha = 1/14.8) +
  geom_segment(data = provi,
               aes(x = POSIC, xend = POSIC, y = 0, yend = locfq, color = MEALP)) +
  scale_colour_gradient(low = dj2pal[3], 
                        high = dj2pal[2],
                        limits = c(0.3, 1),
                        name = "Median alternative allele proportion")+
  geom_text(data = provi[lupfift, ],
            aes(x = POSIC, y = locfq, label = OMMU), #Display alternative alleles corresponding to omicron
            colour = "black",
            nudge_y = 0.06,
            size = 3,
            check_overlap = TRUE) +
  geom_text(data = provi[lupfift, ],
            aes(x = POSIC, y = locfq, label = DEMMU), #Display alternative alleles corresponding to delta
            colour = "blue",
            nudge_y = 0.08,
            size = 3,
            check_overlap = TRUE) +
  geom_text(data = provi[lupfift, ],
            aes(x = POSIC, y = locfq, label = ALMU), #Display alternative alleles corresponding to alpha
            colour = "green",
            nudge_y = 0.1,
            size = 3,
            check_overlap = TRUE) +
  labs(x = "Genomic position", 
       y = "Polymorphic site frequency") +
  theme_minimal()

loc_geno_lpplot <- loc_geno_lpplot +
  facet_wrap(vars(LOCAT),
             nrow = 3) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        text = element_text(size = 9))

#Figure 3.  Mutations in the SARS-CoV-2 reads from South Africa during September-November 2022 by province. 
ggarrange(loc_geno_lpplot, geneplo, ncol = 1, nrow = 2, heights = c(9, 1), align = "v")


#Mpumalanga
#Which was not included in the plot due to the small sample size
mpufq <- fqnmeta %>%
  filter(LOCAT == "Mpumalanga") %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(locfq = n/5)
mpual <- fqnmeta %>%
  filter(LOCAT == "Mpumalanga") %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALternative allele Proportion
         MEAL = median(ALTPR), #MEdian ALternative Proportion
         MELD = median(LD)) #MEdian Locus Depth
mpuma <- inner_join(mpufq, mpual, by = c("POSIC", "ALT1"))
updomi(mpual)

#Omicron mutations in each province
updomi(easca)
updomi(gaute)
updomi(kwaal)
updomi(mpual)

#Which omicron mutations are present in >=50% of the samples in each province
updomi(easca[which(easca$locfq > 0.5), ]) #Eastern Cape
updomi(gaute[which(gaute$locfq > 0.5), ]) #Gauteng
updomi(kwana[which(kwana$locfq > 0.5), ]) #KwaZulu Natal
updomi(mpuma[which(mpuma$locfq > 0.5), ]) #Mpumalanga

#Delta mutations present in each province
#Eastern Cape
unique(easca$DEMMU)
length(unique(easca$DEMMU))
unique(easca[which(easca$locfq > 0.5), "DEMMU"]) #>=50% of the samples
#Gauteng
unique(gaute$DEMMU)
length(unique(gaute$DEMMU))
unique(gaute[which(gaute$locfq > 0.5), "DEMMU"]) #>=50% of the samples
#KwaZulu-Natal
unique(kwana$DEMMU)
length(unique(kwana$DEMMU))
unique(kwana[which(kwana$locfq > 0.5), "DEMMU"]) #>=50% of the samples
#Mpumalanga
unique(mpuma$DEMMU)
length(unique(mpuma$DEMMU))
unique(mpuma[which(mpuma$locfq > 0.5), "DEMMU"]) #>=50% of the samples

#Alpha mutations present in each province
#Eastern Cape
unique(easca$ALMU)
length(unique(easca$ALMU))
unique(easca[which(easca$locfq > 0.5), "ALMU"]) #>=50% of the samples
#Gauteng
unique(gaute$ALMU)
length(unique(gaute$ALMU))
unique(gaute[which(gaute$locfq > 0.5), "ALMU"]) #>=50% of the samples
#KwaZulu-Natal
unique(kwana$ALMU)
length(unique(kwana$ALMU))
unique(kwana[which(kwana$locfq > 0.5), "ALMU"]) #>=50% of the samples
#Mpumalanga
unique(mpuma$ALMU)
length(unique(mpuma$ALMU))
unique(mpuma[which(mpuma$locfq > 0.5), "ALMU"]) #>=50% of the samples
