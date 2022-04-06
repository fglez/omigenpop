
#LIBRARIES

library(stringr)
library(tidyverse)
library(ggplot2)

##EDITING VCF TABLE

#Set path
setwd("/home/mfcg/Descargas/covid/omicron/vcfs")

rawvcftab <- read.delim(file = "sou_omi_SARS.tsv", #Read the tab separated vcf file
                        header = FALSE)

#Edit and separate columns 
edivcftab <- rawvcftab %>%
  separate(V1, into = c("FILE", "CHROM"), sep = "\\ ") %>% #Separate first column filename and reference delimited by a space
  mutate(SAMPL = str_remove_all(FILE, "trimmed_|.tsv")) %>% #Remove the prefix and sufix so the library or sample name can be retrieved
  mutate(POSIC = as.character(V2), #Position
         REFER = as.character(V4), #Reference allele
         ALTER = as.character(V5), #Alternative allele
         QUALI = as.numeric(V6)) %>% #Quality not in phred score
  separate(ALTER, into = c("ALT1", "ALT2", "ALT3"), sep = "\\,") %>% #Separate the alleles versions
  separate(V10, into = c("GT","DP", "AD", "RO", "QR", "AO", "QA", "GL"), sep = "\\:") %>% #Separate the info column
  separate(AD, into = c("refd", "alt1d", "alt2d", "alt3d"), sep = "\\,") %>% #Separate the allele depth values in reference's and alleles' depth
  mutate(RD = as.double(refd), #Reference Depth numeric
         AD1 = as.double(alt1d), #Allele Depth numeric
         AD2 = as.double(alt2d), #Allele 2 Depth numeric, when present
         AD3 = as.double(alt3d))

edivcftab$LD <- rowSums(edivcftab[,c("RD", "AD1", "AD2", "AD3")], na.rm = TRUE) #Locus depth as the sum of the reference's and alleles' depth

write_tsv(edivcftab, 
            file = "sou_SARS_edit.tsv")

#METADATA

#Set path
setwd("/home/mfcg/Descargas/covid/omicron/vcfs")

#Metadata from NCBI
metasa <- read.csv(file = "newSRATAB.txt",
                   header = TRUE)

gdmeta <- metasa %>%
  filter(Bases >= 40045840) %>% #From this number of bases or bigger there are >= 200k reads mapped to the reference genome as suggested by this paper https://www.sciencedirect.com/science/article/pii/S1198743X21001646
  mutate(Location = str_remove_all(geo_loc_name, "South Africa: "), #Location without country
         Month = strftime(as.Date(Collection_date), "%m")) %>% #Extract the month from the collection date
  select(Run, AvgSpotLen, Bases, Organism, geo_loc_name_country, BioProject,
         Collection_date, Month, Library.Name, Location, collected_by) #Choose only useful columns

##FILTERING DATA

#Set path
setwd("/home/mfcg/Descargas/covid/omicron/vcfs")

#Read the edited vcf table
edivcftab <- read.delim (file = "sou_SARS_edit.tsv",
                         header = TRUE)

#The data must meet the following criteria

filvcftab <- edivcftab %>%
  select(SAMPL, POSIC, REFER, ALT1, ALT2, ALT3, QUALI, LD, RD, AD1, AD2, AD3, GT) #Select only the useful columns

filvcftab <- filvcftab[edivcftab$SAMPL %in% gdmeta$Run, ] #Leave only in the vcf table those libraries with >= 200K reads mapped to the reference genome

filvcftab <- filvcftab %>%
  filter(QUALI >= 30, #Phred score of 30
         LD >= 20) %>% #Locus represented by 2000 reads
  mutate(BIALL = is.na(ALT2)) %>% #Is BIALLelic if there are only two alleles (REF and ALT1), more than one alternative alleles result in FALSE
  filter(BIALL == TRUE) %>% #Only biallelic
  select(SAMPL, POSIC, REFER, ALT1, QUALI, LD, RD, AD1, GT, BIALL)

write_tsv(filvcftab,
          file = "fil_sou_sars.tsv")

#Read filtered vcf table
filvcftab <- read.delim(file = "fil_sou_sars.tsv",
                        header = TRUE)
#Filtered metadata

filmeta <- gdmeta[gdmeta$Run %in% unique(filvcftab$SAMPL),]

write_tsv(filmeta,
          file = "fil_metadata.tsv")

#Filtered metadata
filmeta <- read.delim(file = "fil_metadata.tsv",
                        header = TRUE)

#Total number of libraries
libno <- length(unique(filvcftab$SAMPL))

###ANALYSES

##Frequency of position changes
fq_pos_var <- filvcftab %>%
  group_by(POSIC, ALT1) %>% #Group by position and by alternative allele
  summarise(n = n(),
            freq = n / libno)


#Omicron mutations
#Alternative alleles of B.1.1.529 common to the two lineages of omicron
#Nucleotide mutations
omimut <- read.csv(file = "B.1.1.529_nucmut.csv")

fq_pos_var <- fq_pos_var %>%
  rowwise() %>%
  mutate(OMMU = c(ifelse(POSIC %in% omimut$POSIC, "mutation", NA))) #Add mutation to all the mutations that are in the same position as those in the VOC omicron

#Save it in a tsv file
write_delim(x = fq_pos_var,
            file = "fqposvar.tsv")


#file edited by hand in order to add the specific name of the omicron mutation and to discard those that are in the same position but with a different allele
fqposvar <- read.table(file = "fqposvar.tsv",
                       header = TRUE)

#Delta mutations
#Adding the genomic-nucleotide changes
delmut <- read.csv(file = "B.1.617.2_nucmut.csv",
                   header = TRUE)

fqposvar <- fqposvar %>%
  rowwise() %>%
  mutate(DEMU = c(ifelse(POSIC %in% delmut$POSIC, "delta", NA))) #Add mutation to all the mutations that are in the same position of delta mutations

#Save the changes in a tsv file, addition of Delta mutations
write_delim(x = fqposvar, 
            file = "fqposvar.tsv")

#Updating the frequency file with the Delta mutations present
fqposvar <- read.table(file = "fqposvar.tsv",
                       header = TRUE)

#Alpha mutations
#Adding the genomic-nucleotide changes
alpmut <- read.csv(file = "B.1.1.7_nucmut.csv",
                   header = TRUE)

fqposvar <- fqposvar %>%
  rowwise() %>%
  mutate(ALMU = c(ifelse(POSIC %in% alpmut$POSIC, "alpha", NA))) #Add mutation to all the mutations that are in the same position of delta mutations

#Save the changes in a tsv file, addition of Delta mutations
write_delim(x = fqposvar, 
            file = "fqposvar.tsv")

#Updating the frequency file with the present Alpha mutations specified by hand
fqposvar <- read.table(file = "fqposvar.tsv",
                       header = TRUE)
  
#Adding frequencies
fqvcftab <- inner_join(fqposvar, filvcftab, by = c("POSIC", "ALT1")) %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD)) %>% #MEdian Locus Depth
  filter(freq >= 0.05) #Leave only the allele mutations that are present in 0.05 or more samples


#Set path
setwd("/home/mfcg/Descargas/covid/omicron/vcfs")

#Save changes in a tsv file
write_tsv(x = fqvcftab,
            file = "fqvcffil.tsv")
#Save the identified variants mutations with frequencies and proportions, filtered
fqvcftab <- read.delim(file = "fqvcffil.tsv",
           header = TRUE)

#Updating the metadata with only the samples that have mutations present in 0.05 samples
filmeta <- filmeta[filmeta$Run %in% unique(fqvcftab$SAMPL), ]
#Save changes in tsv file
write_tsv(filmeta,
          file = "fil_metadata.tsv")
#Filtered metadata
filmeta <- read.delim(file = "fil_metadata.tsv",
                      header = TRUE)

#Polymorphic sites in all the samples that passed the filters
#Lollipop graph
ggplot(fqvcftab, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  geom_point(size = 4, alpha = 1/15.6) +
  xlim(21563, 25384) + #Only spike region
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.02) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.1) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.18) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(title = "Genomic changes in SARS-CoV-2 in South Africa 2021",
       x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred", 
                        name = "Median alternative allele proportion") +
  theme_bw()

#Omicron mutations prior to split present in the dataset
length(unique(fqvcftab$OMMU))
#Mutations of the Delta constellation presen in the whole dataset
length(unique(fqvcftab$DEMU))
#Mutations of the Alpha constellation present in the whole dataset
unique(fqvcftab$ALMU)

#Adding some meadata
umeta <- filmeta %>%
  select(Run, Month, Location, BioProject)

colnames(umeta) <- c("SAMPL", "MONTH", "LOCAT", "PROJE")
umeta$MONTH <- as.character(umeta$MONTH)

fqnmeta <- inner_join(fqvcftab, umeta, by = "SAMPL")

#Only september data

sept <- fqnmeta %>%
  filter(MONTH == 9)

length(unique(sept$SAMPL))
length(unique(sept$OMMU))

ggplot(sept, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  #xlim(21563, 25384) +
  geom_point(size = 4, alpha = 1/4.5) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(title = "Genomic change in SARS-CoV-2 in September", 
       x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Only october data

octo <- fqnmeta %>%
  filter(MONTH == 10)

length(unique(octo$SAMPL))
length(unique(octo$OMMU))

ggplot(octo, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  #xlim(21563, 25384) +
  geom_point(size = 4, alpha = 1/5.5) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(title = "Genomic change in SARS-CoV-2 in October",
       x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Only November data

nove <- fqnmeta %>%
  filter(MONTH == 11)

length(unique(nove$SAMPL))
length(unique(nove$OMMU))

ggplot(nove, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  #xlim(21563, 25384) +
  geom_point(size = 4, alpha = 1/5.6) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.02) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(title = "Genomic change in SARS-CoV-2 in November",
       x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Only omicron data

onomi <- fqnmeta %>%
  filter(PROJE == "PRJNA784038")

length(unique(onomi$SAMPL))

onomi <- onomi %>%
  group_by(POSIC, ALT1) %>%
  summarise(n=n()) %>%
  mutate(omifq = n/59)

ggplot(onomi, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  geom_point(size = 4, alpha = 1/2.95) +
  #xlim(21563, 25384) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(title = "Genomic changes in omicron dataset in South Africa 2021",
       x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred") +
  theme_bw()

ggplot(data = onomi, aes(x = omifq, color = "PRJNA784038")) +
  geom_density() +
  labs(title ="Density of the polymorphic sites frequency of the omicron samples",
       x = "Polymorphic Site Frequency", y = "Density") +
  theme_bw()

#Genomic changes per province

#Gauteng

gaute <- fqnmeta %>%
  filter(LOCAT == "Gauteng")

length(unique(gaute$SAMPL))
length(unique(gaute$OMMU))

ggplot(gaute, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  geom_point(size = 4, alpha = 1/55) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Mpumalanga

mpuma <- fqnmeta %>%
  filter(LOCAT == "Mpumalanga")

length(unique(mpuma$SAMPL))
length(unique(mpuma$OMMU))

ggplot(mpuma, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  geom_point(size = 4, alpha = 1/8) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Eastern Cape

eacap <- fqnmeta %>%
  filter(LOCAT == "Eastern Cape")

length(unique(eacap$SAMPL))
length(unique(eacap$OMMU))

ggplot(eacap, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  geom_point(size = 4, alpha = 1/8) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Kwazulu-Natal

kwana <- fqnmeta %>%
  filter(LOCAT == "KwaZulu-Natal")

length(unique(kwana$SAMPL))
length(unique(kwana$OMMU))

ggplot(kwana, aes(x = POSIC, y = freq, color = MEALP)) + #The alternative allele proportion is displayed by color code
  geom_point(size = 4, alpha = 1/8) +
  geom_text(aes(label = OMMU), colour = "black", nudge_y = 0.025) +
  geom_text(aes(label = DEMU), colour = "deepskyblue4", nudge_y = 0.05) +
  geom_text(aes(label = ALMU), colour = "darkolivegreen", nudge_y = 0.075) +
  geom_segment(aes(x = POSIC, xend = POSIC, y = 0, yend = freq)) +
  labs(x = "Genomic position", y = "Polymorphic site frequency") +
  scale_colour_gradient(low = "darksalmon", high = "darkred",
                        name = "Median alternative allele proportion") +
  theme_bw()

#Frequencies taking into account the different sites

table(umeta$LOCAT) #number of samples for each location

altfq_popsite <- fqnmeta %>%
  group_by(POSIC, ALT1, LOCAT) %>%
  summarise(n=n()) %>%
  mutate(locfq = case_when(
    LOCAT == "Eastern Cape" ~ n/72,
    LOCAT == "Gauteng" ~ n/111,
    LOCAT == "KwaZulu-Natal" ~ n/124,
    LOCAT == "Mpumalanga" ~ n/5))

ggplot(data = altfq_popsite, aes(x = locfq, group = LOCAT, fill = LOCAT)) +
  geom_density(alpha = 0.25) +
  labs(title ="Density of the polymorphic sites frequency for each location",
       x = "Polymorphic Site Frequency", y = "Density") +
  theme_bw()

kruskal.test(locfq ~ LOCAT, data = altfq_popsite)

#Frequencies taking into account the different months

table(umeta$MONTH) #Samples per month

altfq_popmont <- fqnmeta %>%
  group_by(POSIC, ALT1, MONTH) %>%
  summarise(n=n()) %>%
  mutate(locfq = case_when(
    MONTH == 9 ~ n/90,
    MONTH == 10 ~ n/110,
    MONTH == 11 ~ n/112))

ggplot(data = altfq_popmont, aes(x = locfq, group = MONTH, fill = MONTH)) +
  geom_density(alpha = 0.3) +
  labs(title ="Density of the polymorphic sites frequency for each month",
       x = "Polymorphic Site Frequency", y = "Density") +
  theme_bw()

kruskal.test(locfq ~ MONTH, data = altfq_popmont)

#Frequencies change in time for the two provinces where omicron was first registered

table(umeta$LOCAT, umeta$MONTH)

#Gauteng
gautime <- fqnmeta %>%
  filter(LOCAT == "Gauteng") %>%
  group_by(POSIC, ALT1, MONTH) %>%
  summarise(n = n()) %>%
  mutate(timefq = case_when(
    MONTH == 9 ~ n/36,
    MONTH == 10 ~ n/15,
    MONTH == 11 ~ n/60)) %>%
  ggplot(aes(x = timefq, group = MONTH, fill = MONTH)) +
  geom_density(alpha = 0.33) +
  labs(title ="Gauteng density of the polymorphic sites frequency by month",
       x = "Polymorphic Site Frequency", y = "Density") +
  theme_bw()

#KwaZulu-Natal
kwatime <- fqnmeta %>%
  filter(LOCAT == "KwaZulu-Natal") %>%
  group_by(POSIC, ALT1, MONTH) %>%
  summarise(n = n()) %>%
  mutate(timefq = case_when(
    MONTH == 9 ~ n/8,
    MONTH == 10 ~ n/64,
    MONTH == 11 ~ n/52)) %>%
  ggplot(aes(x = timefq, group = MONTH, fill = MONTH)) +
  geom_density(alpha = 0.33) +
  labs(title ="KwaZulu-Natal density of the polymorphic sites frequency by month",
       x = "Polymorphic Site Frequency", y = "Density") +
  theme_bw()

#AMOVA entre poblaciones y dentro de poblaciones c/Provincia

ggplot(data = fqnmeta, aes(x = ALTPR, y = LOCAT, fill = LOCAT)) +
  xlim(0.995, 1) +
  geom_boxplot(width = 0.8, lwd = 0.2) +
  labs(x = "Genome-wide alternative allele proportion", y = "Location") +
  theme_bw()

kruskal.test(ALTPR~LOCAT, data = fqnmeta)

ggplot(data = fqnmeta, aes(x = ALTPR, y = MONTH, fill = MONTH)) +
  xlim(0.995, 1) +
  geom_boxplot(width = 0.8, lwd = 0.2) +
  labs(x = "Genome-wide alternative allele proportion", y = "Month") +
  theme_bw()

kruskal.test(ALTPR~MONTH, data = fqnmeta)

#Genetic regions SARS-CoV-2

sc2gg <- read.csv(file = "SC2-gene-posic.csv")
ggplot() +
  geom_rect(data = sc2gg, aes(xmin = nuc_init, xmax = nuc_end, ymin= 0, ymax = 1, fill = Gen), alpha = 0.25)

barplot(height = c(rep(1, 11)), 
        names = sc2gg$Gen, 
        width = sc2gg$nuc_end-sc2gg$nuc_init)
