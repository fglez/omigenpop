
#libraries
#library(stringr)
library(tidyverse)
library("poppr")

#Path to files
setwd("/home/mfcg/Descargas/covid/omicron/vcfs")

#Open files with filtered identified variantes with frequencies and proportions
fqvcftab <- read.delim(file = "fqvcffil.tsv",
                       header = TRUE)
#Corresponding metadata
filmeta <- read.delim(file = "fil_metadata.tsv",
                      header = TRUE)

##Making of the genind file 

#Matrix of individuals, alleles and allele counts
#Retrieve in long format the individuals sample id, posicion of allele, allele depth
longiac <- fqvcftab[ , c("SAMPL", "POSIC", "AD1")]
#Transform long format to wide format, this is a dataframe
matriac <- spread(longiac, POSIC, AD1)
#Move the sample id comlumn to rownames
matriac <- column_to_rownames(matriac, var = "SAMPL")
#Change NAs to 0s so the counts make sense
matriac[apply(matriac, 2, is.na)] <- 0
#As matrix
matriac <- as.matrix(matriac)

#Table with the sample ids and population location
popind <- filmeta[ , c("Run", "Location")]
#Arrange in decreasing order by sample id
popind <- arrange(popind, Run)
#Retrieve only the location vector and change its class to factor
popind <- popind$Location
popind <- as.factor(popind)

#Strata data frame
stromi <- data.frame(Pop = popind) 

#Making an integer vector that indicates an haploid ploidy
vecplo <- as.integer(rep(1, 312))


#To geneind format
gentry <- genind(tab = matriac, 
                 pop = popind, 
                 ploidy = vecplo,
                 type = "codom", #This data are allele counts, so codominant is the type
                 strata = stromi)

#From the example by N Kamvar, SE Everhart and NJ Grünwald is followed
#it can be found here https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html
#To genclone
clogen <- as.genclone(gentry)
#Number of individuals per province population
table(strata(gentry, ~Pop))

#Amova analysis of the South Africa dataset from 2021 for the allele counts
saamova <- poppr.amova(clogen, ~Pop)

#Significance
set.seed(100)
sasigni <- randtest(saamova, nrepet = 1000)
