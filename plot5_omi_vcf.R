##GENETIC DIVERSITY MEASURES

#libraries
library("tidyverse")
library("poppr")
library("dplyr")
library("adegenet")
library("hierfstat")
library("wesanderson")
library("gridExtra")
library("ggpubr")


#Path to files
setwd("/home/mfcg/Descargas/covid/omicron/ver2gendiv")

#Filtered data and metadata 
newfqnmt <- read.csv(file = "newfqnmeta.csv",
                     header = TRUE)
#Filtered metadata
newfilmet <- read.delim(file = "newfilmet.tsv",
                        header = TRUE)

#Removing Mpumalanga
#Removing Mpumalanga from the variants information and metadata
thrfm <- newfqnmt %>%
  filter(LOCAT != "Mpumalanga") %>%
  select(SAMPL, POSIC, ALT1, ALTAL, REFER, GT, TYPE) %>% #Useful columns: sample name, position in genome, alternative allele, alternative allele name, reference allele, genotype
  mutate(REFAL = paste0(POSIC, ".", REFER)) %>% #Column with reference allele name
  filter(TYPE == "TYPE=snp")
#Removing Mpumalanga from the metadata
metthr <- newfilmet %>%
  filter(Location != "Mpumalanga")

#Genome positions with variants
poswvar <- unique(thrfm$POSIC)

#DOSAGE
#Table with genotypes of the individuals
#rows as samples and columns the alternative alleles at locus
# For alternatve alleles genotypes 1/1 and 0/1 will be 1 and NA will be 0
#Preparation objects
altlng <- thrfm %>%
  select(SAMPL, POSIC, GT)

#Individuals, Alternative alleles and genotypes to sum up to ploidy, haploid = 1
altdf <- spread(altlng, POSIC, GT)
alemat <- column_to_rownames(altdf, var = "SAMPL")
#Dosage all samples
alemat <- ifelse(is.na(alemat), 0, 1)


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

#List to save the fst valus for each snv
stfst <- vector(mode = "list", length = length(poswvar))

#Function to calculate the fst value for each snp or snv
for(i in 1:length(poswvar)){
  sitfst <- fs.dosage(dos = alemat[ , i], pop = popvec)
  solofs <- sitfst[["Fs"]][2, ]
  stfst[[i]] <- solofs
}

#Making a data frame to save the fst values for each snv
names(stfst) <- poswvar
stfstdf <- as.data.frame(t(as.data.frame(stfst)))

poswvar <- as.numeric(poswvar)
stfstdf <- cbind(stfstdf, poswvar)

#Plotting 
plot(x = stfstdf$poswvar, y = stfstdf$All, type = "b")

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

#Omicron mutation file
omimut <- read.csv(file = "B.1.1.529_nucmut.csv")
#Adding omicron mutations names to site fsts table
#sfsttab <- stfstdf %>%
#  rowwise() %>%
#  mutate(OMMU = c(ifelse(Position %in% omimut$POSIC, "mutation", NA)))
#write_delim(x = sfsttab,
#            file = "sfstab.tsv")
#Modified by hand to add the specific amino acid change of Omicron
sfsttab <- read.table(file = "sfstab.tsv",
                      header = TRUE)
omind <- which(sfsttab$OMMU != is.na(sfsttab$OMMU))

inomfst <- sfsttab %>%
  filter(OMMU != is.na(OMMU)) %>%
  filter(All > 0)

#Global FST per SNV
fstps_plo <- ggplot(data = sfsttab, aes(x = Position, y = All)) +
  ylim(-1, 1) +
  geom_point(color = dj2pal[4]) +
  geom_text(data = inomfst,
            aes(x = Position, y = All, label = OMMU),
            colour = "black",
            size = 3,
            angle = 45,
            check_overlap = TRUE) +
  labs(title = "Global", 
       x = "Genome position",
       y = "Fst value per SNV") +
  theme_minimal()

ggarrange(fstps_plo, geneplo, ncol = 1, nrow = 2, heights = c(2, 1), align = "v")

fstps_plo1 <- ggplot(data = sfsttab, aes(x = Position, y = EasternCape)) +
  ylim(-1.5, 1.5) +
  geom_point(color = dj2pal[1]) +
  geom_text(data = inomfst,
            aes(x = Position, y = EasternCape, label = OMMU),
            colour = "black",
            size = 3,
            angle = 45,
            check_overlap = TRUE) +
  labs(title = "Eastern Cape",
       x = "Genome position",
       y = "Fst value per SNV") +
  theme_minimal()

fstps_plo2 <- ggplot(data = sfsttab, aes(x = Position, y = Gauteng)) +
  ylim(-1.5, 1.5) +
  geom_point(color = dj2pal[2]) +
  geom_text(data = inomfst,
            aes(x = Position, y = Gauteng, label = OMMU),
            colour = "black",
            size = 3,
            angle = 45,
            check_overlap = TRUE) +
  labs(title = "Gauteng",
       x = "Genome position",
       y = "Fst value per SNV") +
  theme_minimal()

fstps_plo3 <- ggplot(data = sfsttab, aes(x = Position, y = KwaZuluNatal)) +
  ylim(-1.5, 1.5) +
  geom_point(color = dj2pal[3]) +
  geom_text(data = inomfst,
            aes(x = Position, y = KwaZuluNatal, label = OMMU),
            colour = "black",
            size = 3,
            angle = 45,
            check_overlap = TRUE) +
  labs(title =  "KwaZulu-Natal",
       x = "Genome position",
       y = "Fst value per SNV") +
  theme_minimal()

ggarrange(fstps_plo, geneplo, fstps_plo1, fstps_plo2, fstps_plo3, ncol = 1, nrow = 5, heights = c(3, 1, 2, 2, 2), align = "v",
          labels = c("A", "", "B", "C", "D"), font.label = list(size = 9))


