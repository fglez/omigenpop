#LIBRARIES

library("stringr")
library("tidyverse")
library("ggplot2")
library("wesanderson")
library("gridExtra")
library("ggpubr")

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
            size = 3,
            angle = 90) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 2))

#FULL DATASET
#Polymorphic sites in all the samples that passed the filters
#Lollipop graph
upfift <- which(fqnmeta$freq > 0.5)

full_geno_lpplot <- ggplot() +
  xlim(0, 30000) +
  geom_point(data = fqnmeta,
             aes(x = POSIC, y = freq, color = MEALP), #The alternative allele proportion is displayed by color code
             size = 4, 
             alpha = 1/50.4) +
  geom_segment(data = fqnmeta,
               aes(x = POSIC, xend = POSIC, y = 0, yend = freq, color = MEALP)) +
  scale_colour_gradient(low = dj2pal[3], 
                        high = dj2pal[2],
                        limits = c(0.3, 1),
                        name = "Median alternative allele proportion")+
  geom_text(data = fqnmeta[upfift, ],
            aes(x = POSIC, y = freq, label = OMMU), #Display alternative alleles corresponding to omicron
            colour = "black", 
            nudge_y = 0.02,
            size = 3,
            check_overlap = TRUE) +
  geom_text(data = fqnmeta[upfift, ],
            aes(x = POSIC, y = freq, label = DEMMU), #Display alternative alleles corresponding to delta
            colour = "blue", 
            nudge_y = 0.08,
            size = 3,
            check_overlap = TRUE) +
  geom_text(data = fqnmeta[upfift, ],
            aes(x = POSIC, y = freq, label = ALMU), #Display alternative alleles corresponding to alpha
            colour = "green", 
            nudge_y = 0.06,
            size = 3,
            check_overlap = TRUE) +
  labs(x = "Genomic position", 
       y = "Alternative alleles frequency") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.box = "horizontal",
        text = element_text(size = 10))

##FIGURE 1. Genomic changes across the SARS-CoV-2 reads of South Africa during 2021. 
ggarrange(full_geno_lpplot, geneplo, ncol = 1, nrow = 2, heights = c(8, 2), align = "v", font.label = list(size = 12)) #Mutation labels


nolab_geno_lpplot <- ggplot() +
  xlim(0, 30000) +
  geom_point(data = fqnmeta,
             aes(x = POSIC, y = freq, color = MEALP), #The alternative allele proportion is displayed by color code
             size = 4, 
             alpha = 1/50.4) +
  geom_segment(data = fqnmeta,
               aes(x = POSIC, xend = POSIC, y = 0, yend = freq, color = MEALP)) +
  scale_colour_gradient(low = dj2pal[3], 
                        high = dj2pal[2],
                        #limits = c(0.5, 1),
                        name = "Median alternative allele proportion")+
  labs(x = "Genomic position", 
       y = "Alternative alleles frequency") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.box = "horizontal",
        text = element_text(size = 10))

ggarrange(nolab_geno_lpplot, geneplo, ncol = 1, nrow = 2, heights = c(8, 2), align = "v", font.label = list(size = 12)) #Without mutation labels

####Not in the article

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
updomi(fqnmeta[upfift, ])

#Barplots
#Subclonal or clonal mutations
#As defined in this article
#subclonal if it has a variant allele frequency (VAF) between 0.05 and 0.95 and is sup- ported by at least 5 variant reads
#clonal mutation, i.e. a mutation supported by at least 5 variants and a VAF of at least 0.95
#As all my variants are above 20 reads per variant 
#The only thing left is to classify them per alternative allele frequency/proportion "ALTPR"

fqnmeta <- fqnmeta %>%
  mutate(suboclo = ifelse(ALTPR < 0.95, "subclonal", "clonal"))

sbind <- fqnmeta %>%
  count(SAMPL, suboclo)
sbind <- spread(sbind, suboclo, n)
sbind[which(is.na(sbind$subclonal)), "subclonal"] <- 0
sbind <- sbind %>%
  mutate("Submut_pres" = ifelse(subclonal == 0, "Absent", "Present"))

table(sbind$Submut_pres)
104/448*100 #Absent percentage
344/448*100 #Present percentage
table(sbind$subclonal)
165/344*100 #Percentage of individuals with one subclonal mutation

sub_ind_plo <- ggplot(data = sbind) +
  geom_bar(aes(x = Submut_pres, y = stat(count)),
           fill = dj2pal[c(3,2)]) +
  labs(x = "Subclonal mutations",
       y = "Number of individuals") +
  theme_minimal() +
  theme(text = element_text(size = 9))

length(unique(sbind$subclonal))
sub_sam_plo <- ggplot(data = sbind) +
  geom_bar(aes(x = subclonal, y = stat(count)),
           fill = dj2pal[c(3, rep(2, 20))]) +
  labs(x = "Number of subclonal mutations",
       y = "Number of individuals") +
  theme_minimal() +
  theme(text = element_text(size = 9))

#Arranging the three plots in one for the final figure
ggarrange(ggarrange(full_geno_lpplot, geneplo, ncol = 1, nrow = 2, heights = c(9, 1), align = "v", labels = "A", font.label = list(size = 9)),
          ggarrange(sub_ind_plo, sub_sam_plo, ncol = 2, nrow = 1, widths = c(1, 3), labels = c("B", "C"), font.label = list(size = 9)),
          nrow = 2, ncol = 1, heights = c(0.75, 0.25))
  

