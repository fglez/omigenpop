#LIBRARIES

library(stringr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(wesanderson)
library(gridExtra)
library(ggpubr)

#Palettes

dj2pal <- wes_palette("Darjeeling2", 5, type = "discrete")
dj1pal <- wes_palette("Darjeeling1", 11, type = "continuous")

#Polymorphisms found vs locus depth

newfqvcf <- read.delim(file = "bnewfqvcf.tsv",
                       header = TRUE)

samdeps <- newfqvcf %>%
  count(LD)

plot(samdeps)

ggplot(data = samdeps, aes(x = LD, y = n)) +
  geom_point(color = dj2pal[2]) +
  labs(x = "Reads Depth",
       y = "Number of Nucleotide Variants Found") +
  theme_minimal()

cor.test(x = samdeps$LD, y = samdeps$n, method = "pearson")

#SNVs found vs locus depth

snvdeps <- newfqvcf %>%
  filter(TYPE == "TYPE=snp") %>%
  count(LD)

ggplot(data = snvdeps, aes(x = LD, y = n)) +
  geom_point(color = dj2pal[2]) +
  labs(x = "Reads Depth",
       y = "Number of SNVs Found") +
  theme_minimal()

plot(snvdeps)
cor.test(x = snvdeps$LD, y = snvdeps$n, method = "pearson")


