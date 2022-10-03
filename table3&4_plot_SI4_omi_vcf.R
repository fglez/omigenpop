##GENETIC DIVERSITY MEASURES

#libraries
library("tidyverse")
library("poppr")
library("dplyr")
library("adegenet")
#library("pegas")
library("hierfstat")

#Path to files
setwd("/home/mfcg/Descargas/covid/omicron/ver2gendiv")

#Filtered data and metadata 
newfqnmt <- read.csv(file = "newfqnmeta.csv",
                     header = TRUE)
#Filtered metadata
newfilmet <- read.delim(file = "newfilmet.tsv",
                        header = TRUE)

#Provinces sample size
table(newfilmet$Location)
#Mpumalanga has a way smaller sample size of five

#Removing Mpumalanga from the variants information and metadata
thrfm <- newfqnmt %>%
  filter(LOCAT != "Mpumalanga") %>%
  select(SAMPL, POSIC, ALT1, ALTAL, REFER, GT, TYPE) %>% #Useful columns: sample name, position in genome, alternative allele, alternative allele name, reference allele, genotype
  mutate(REFAL = paste0(POSIC, ".", REFER)) #Column with reference allele name
#Removing Mpumalanga from the metadata
metthr <- newfilmet %>%
  filter(Location != "Mpumalanga")
#Samples of each province
kwasam <- newfilmet$Run[which(newfilmet$Location == "KwaZulu-Natal")] #KwaZulu-Natal
gausam <- newfilmet$Run[which(newfilmet$Location == "Gauteng")] #Gauteng
eacsam <- newfilmet$Run[which(newfilmet$Location == "Eastern Cape")] #Eastern Cape

#DOSAGE
#Table with genotypes of the individuals
#rows as samples and columns the alternative alleles at locus
# For alternatve alleles genotypes 1/1 and 0/1 will be 1 and NA will be 0
#Preparation objects
altlng <- thrfm %>%
  filter(TYPE == "TYPE=snp") %>%
  select(SAMPL, ALTAL, GT)

#Individuals, Alternative alleles and genotypes to sum up to ploidy, haploid = 1
altdf <- spread(altlng, ALTAL, GT)
alemat <- column_to_rownames(altdf, var = "SAMPL")
#Dosage all samples
alemat <- ifelse(is.na(alemat), 0, 1)

#Dosage per province
eacdos <- alemat[eacsam, ]
gaudos <- alemat[gausam, ]
kwados <- alemat[kwasam, ]

#MEASURES OF GENETIC DIVERSITY
#Obtained from dosage
#Pi per province
pi.dosage(dos = kwados, L=29903) #Kwazulu-Natal
pi.dosage(dos = gaudos, L=29903) #Gauteng
pi.dosage(dos = eacdos, L=29903) #Eastern Cape

#Tajima from dosage
TajimaD.dosage(dos = alemat) #From the full dataset
#Tajima's D per province
TajimaD.dosage(dos = kwados)
TajimaD.dosage(dos = gaudos)
TajimaD.dosage(dos = eacdos)

#Adding population
ppnd <- metthr %>%
  select(Run, Location) %>%
  mutate(LOCID = case_when(
    Location == "Eastern Cape" ~ 1,
    Location == "Gauteng" ~ 2,
    Location == "KwaZulu-Natal" ~ 3))
colnames(ppnd) <- c("SAMPL", "LOCATION", "LOCID")
sahf <- inner_join(altdf, ppnd, by = "SAMPL")
popvec <- sahf$LOCID

#FSTs
safst <- fs.dosage(dos = alemat,
          pop = popvec)
safst$Fs #Province specific Fsts and overall Fsts
safst$Fst2x2 #Pairwise Fsts between provinces

#For each site with a polymorphism
sites <- unique(alemat$POSIC)


#GENIND
#The main part of the genind requires a table with all the alleles of a loci
#And its count up to the ploidy (1 for haploid)
#There is already a table with the alternative alleles of type SNP
View(alemat)
#Making the one for reference alleles
#Genotype 1/1 means absence of reference alleles so is 0, 0/1 means presence thus 1
#Reference alleles for SNPs
reflng <- thrfm %>%
  filter(TYPE == "TYPE=snp") %>%
  select(SAMPL, REFAL, GT)
#Individuals, reference alleles and genotypes to sum up to ploidy, haploid = 1
refdf <- spread(reflng, REFAL, GT)
refmat <- column_to_rownames(refdf, var = "SAMPL")
refmat <- ifelse(refmat == "1/1", 0, 1)
#As NA means absence of the alternative allele thus it means the reference allele is present so is 1
refmat[is.na(refmat) == TRUE] <- 1

#Joining the reference mat and the alternative allele mat
allmat <- cbind(refmat, alemat)
allmat <- allmat[ , order(colnames(allmat))] #Putting together alleles of the same loci

#Loci reference to each allele
locre <- thrfm %>%
  filter(TYPE == "TYPE=snp") %>%
  select(POSIC, REFAL) %>%
  distinct() %>%
  mutate(LOCUS = as.character(POSIC))
local <- thrfm %>%
  filter(TYPE == "TYPE=snp") %>%
  select(POSIC, ALTAL) %>%
  distinct() %>%
  mutate(LOCUS = as.character(POSIC))
colnames(locre) <- c("POSIC", "ALNAM", "LOCUS")
colnames(local) <- c("POSIC", "ALNAM", "LOCUS")
localle <- rbind(locre, local)
localle <- localle[order(localle$ALNAM), -1]
#Factor that corresponds alleles with their loci
lcfc <- as.factor(localle$LOCUS)
#Alleles per locus
lcnll <- table(localle$LOCUS)
lcnm <- names(lcnll)
lcnll <- as.integer(lcnll)
names(lcnll) <- lcnm
#List with locus and their alleles
llnms <- split(localle$ALNAM, localle$LOCUS)

#Ploidy
plovec <- rep(1, 504)

#Population
popfc <- as.factor(ppnd$LOCATION)
popstr <- as.data.frame(popfc)
colnames(popstr) <- "Pop"

#Genind
sagenind <- genind(tab = allmat,
                   pop = popfc,
                   ploidy =  plovec, 
                   type = "codom",
                   strata = popstr)

sagenind@loc.fac <- lcfc
sagenind@loc.n.all <- lcnll
sagenind@all.names <- llnms

saamova <- poppr.amova(sagenind, ~Pop)

set.seed(100)
randtest(saamova)

