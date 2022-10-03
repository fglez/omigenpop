##ACQUIRING INFORMATION FROM THE VARIANT CALLING FILES

#LIBRARIES

library(stringr)
library(tidyverse)

#EDITING VCF TABLE

#Set path
setwd("/home/mfcg/Descargas/covid/omicron/ver2gendiv")

rawvcftab <- read.delim(file = "sou_omi_SARS.tsv", #Read the tab separated vcf file
                        header = FALSE)

#Edit and separate columns 
edivcftab <- rawvcftab %>%
  separate(V1, into = c("FILE", "CHROM"), sep = "\\ ") %>% #Separate first column filename and reference delimited by a space
  mutate(SAMPL = str_remove_all(FILE, "trimmed_|.tsv")) %>% #Remove the prefix and sufix so the library or sample name can be retrieved
  mutate(POSIC = as.character(V2), #Position
         REFER = as.character(V4), #Reference allele
         ALTER = as.character(V5), #Alternative allele
         QUALI = as.numeric(V6)) %>% #Quality 
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

#Metadata from NCBI
metasa <- read.csv(file = "newSRATAB.txt",
                   header = TRUE)

gdmeta <- metasa %>%
  mutate(Lctn = str_remove_all(geo_loc_name, "South Africa: "), #Location without country
         Month = strftime(as.Date(Collection_date), "%m")) %>% #Extract the month from the collection date
  mutate(Location = str_replace_all(Lctn, "Kwazulu", "KwaZulu")) %>% #Homogenize names
  select(Run, AvgSpotLen, Bases, Organism, geo_loc_name_country, BioProject,
         Collection_date, Month, Library.Name, Location, collected_by) #Choose only useful columns

#EDITING MORE EXTRACTING INFO

edivcftab <- read.delim (file = "sou_SARS_edit.tsv",
                         header = TRUE)

vcfmrinf <- c("AB", "ABP", "AC", "AF", "AN", "AO", "CIGAR", "DP", "DPB", "DPRA", 
              "EPP", "EPPR", "GTI", "LEN", "MEANALT", "MQM", "MQMR", "NS",
              "NUMALT", "ODDS", "PAIRED", "PAIREDR", "PAO", "PQA", "PQR", "PRO", 
              "QA", "QR", "RO", "RPL", "RPP", "RPPR", "RPR", "RUN",
              "SAF", "SAP", "SAR", "SRF", "SRP", "SRR", "TYPE") #Info titles

newedivcf <- edivcftab %>%
  separate(V8, into = vcfmrinf, sep = ";")

newedivcf <- newedivcf %>%
  select(SAMPL, POSIC, REFER, ALT1, ALT2, ALT3, QUALI, LD, RD, AD1, AD2, AD3, GT, TYPE)

newedivcf <- newedivcf[newedivcf$SAMPL %in% gdmeta$Run, ]

#FILTERING DATA

newfilvcf <- newedivcf %>%
  filter(QUALI >= 30, #Phred score of 30
         LD >= 20) %>% #Locus represented by 2000 reads
  mutate(BIALL = is.na(ALT2)) %>% #Is BIALLelic if there are only two alleles (REF and ALT1), more than one alternative alleles result in FALSE
  filter(BIALL == TRUE) %>% #Only biallelic
  select(SAMPL, POSIC, REFER, ALT1, QUALI, LD, RD, AD1, GT, TYPE) %>%
  mutate(actal = paste0(REFER, POSIC, ALT1))

write_tsv(x = newfilvcf,
          file = "newfilvcf.tsv")

#FREQUENCIES

newfilvcf <- read.delim(file = "newfilvcf.tsv",
                        header = TRUE)

libno <- length(unique(newfilvcf$SAMPL))

newfqposv <- newfilvcf %>%
  group_by(POSIC, ALT1) %>% #Group by position and by alternative allele
  summarise(n = n(),
            freq = n / libno) %>%
  mutate(ALTAL = paste0(POSIC, ".", ALT1))
length(unique(newfqposv$ALTAL)) #Number of mutations before the polimorphism frequency

#Adding mutations
#Omicron mutations
#Alternative alleles of B.1.1.529 common to the two lineages of omicron
#Nucleotide mutations
omimut <- read.csv(file = "B.1.1.529_nucmut.csv") %>%
  mutate(ALTAL = paste0(POSIC, ".", ALTER))
newfqposv <- newfqposv %>%
  rowwise() %>%
  mutate(OMMU = c(ifelse(ALTAL %in% omimut$ALTAL, "mutation", NA))) #Add mutation to all the mutations that are in the same position as those in the VOC omicron
#Save it in a tsv file
write_delim(x = newfqposv,
            file = "newfqposv.tsv")
#file edited by hand in order to add the specific name of the omicron mutation and to discard those that are in the same position but with a different allele
newfqposv <- read.table(file = "newfqposv.tsv",
                       header = TRUE)
#Delta mutations
#Adding the genomic-nucleotide changes
delmut <- read.csv(file = "B.1.617.2_nucmut.csv",
                   header = TRUE) %>%
  mutate(ALTAL = paste0(POSIC, ".", ALTER))
newfqposv <- newfqposv %>%
  rowwise() %>%
  mutate(DEMMU = c(ifelse(ALTAL %in% delmut$ALTAL, "mutation", NA))) #Add mutation to all the mutations that are in the same position as those in the VOC DELTA
#Save it in a tsv file
write_delim(x = newfqposv,
            file = "newfqposv.tsv")
#file edited by hand in order to add the specific name of the omicron mutation and to discard those that are in the same position but with a different allele
newfqposv <- read.table(file = "newfqposv.tsv",
                        header = TRUE)
#Alpha mutations
#Adding the genomic-nucleotide changes
alpmut <- read.csv(file = "B.1.1.7_nucmut.csv",
                   header = TRUE) %>%
  mutate(ALTAL = paste0(POSIC, ".", ALTER))
newfqposv <- newfqposv %>%
  rowwise() %>%
  mutate(ALMU = ifelse(ALTAL %in% alpmut$ALTAL, "mutation", NA)) #Add mutation to all the mutations that are in the same position as those in the VOC ALPHA
#Save it in a tsv file
write_delim(x = newfqposv,
            file = "newfqposv.tsv")
#Updating the frequency file with the present Alpha mutations specified by hand
newfqposv <- read.table(file = "newfqposv.tsv",
                        header = TRUE)

newfqvcf <- inner_join(newfqposv, newfilvcf, by = c("POSIC", "ALT1")) %>%
  group_by(POSIC, ALT1) %>%
  mutate(ALTPR = AD1/LD, #ALTernative allele PRoportion
         MEALP = median(ALTPR), #MEdian ALTernative Proportion
         MELD = median(LD)) %>% #MEdian Locus Depth
  filter(freq >= 0.05) #Leave only the allele mutations that are present in 0.05 or more samples

#Save changes in a tsv file
write_tsv(x = newfqvcf, 
          file="bnewfqvcf.tsv")    

newfqvcf <- read.delim(file = "bnewfqvcf.tsv",
                       header = TRUE)
length(unique(newfqvcf$ALTAL)) #Number of mutations after the polimorphism frequency filtering

#FILTERING METADATA
newfilmet <- gdmeta[gdmeta$Run %in% unique(newfqvcf$SAMPL), ] #Only samples that passed the filters
#write_tsv(x = newfilmet,
#          file = "newfilmet.tsv")
newfilmet <- read.delim(file = "newfilmet.tsv",
                        header = TRUE)
length(unique(newfilmet$Run)) #Number of libraries after the filtering
newfilmet %>% 
  count(Location)

#JOINING METADATA A VARIANT DATA
umeta <- newfilmet %>%
  select(Run, Month, Location, BioProject)
colnames(umeta) <- c("SAMPL", "MONN", "LOCAT", "PROJE")
umeta <- umeta %>%
  mutate(MONTH = case_when( #Adding months name
    MONN == 9 ~ "September",
    MONN == 10 ~ "October",
    MONN == 11 ~ "November")
    )

#Put together metadata and vcf data
newfqnmt <- inner_join(newfqvcf, umeta, by = "SAMPL")

#Save it to a csv file
write_csv(x = newfqnmt,
          file = "bnewfqnmeta.csv")

#Read the new file
newfqnmt <- read.csv(file = "bnewfqnmeta.csv",
                     header = TRUE)

length(unique(newfqnmt$ALTAL)) #Different variants found
length(unique(newfqnmt$OMMU)) #Number of Omicron mutations in the dataset, including NA
uniomi <- unique(newfqnmt$OMMU)
oldomi <- c("Q493R", "D1146D", "T64T", "RG203KR") #No longer part of B.1.1.529 lineage defining mutations
uniomi <- uniomi[!uniomi %in% oldomi]
length(uniomi) #Updated number of mutations of B.1.1.529 lineage in dataset
#Delta mutations in dataset
length(unique(newfqnmt$DEMMU)) #Delta mutations in the dataset, including NA
unidel <- unique(newfqnmt$DEMMU)
#Mutation shared by Omicron prior to split and Delta
uniomi[which(uniomi %in% unidel)]
unidel[which(unidel %in% uniomi)]
#Alpha mutations in dataset
length(unique(newfqnmt$ALMU)) #Number of Alpha mutations in dataset, including NA
unial <- unique(newfqnmt$ALMU)
#Mutations shared by Omicron prior to split and Alpha
unial[which(unial %in% uniomi)]
uniomi[which(uniomi %in% unial)]

