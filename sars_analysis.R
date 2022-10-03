library(tidyverse)

setwd("/Users/erik/SARS_refs/var_ana/ncovIllumina_sequenceAnalysis_callVariants/SNPS_pop/")

var_data_2 <- read_tsv("bajio_SARS_pop.tsv", col_names = F) %>% 
  rename(POS = X2, REF = X4, ALT = X5, qual = X6, format = X8, geno_info = X10) %>%
  separate(X1, into = c("ID_ind", "CHROM"), sep = "\\ ") %>%
  rowwise() %>%
  mutate(allele_depth = map_chr(
    geno_info, function(s) rev(strsplit(s, ":")[[1]])[6])) %>%
  separate(allele_depth, into = c("ref_d", "alt_d"), sep = "\\,") %>%
  rowwise() %>%
  mutate(POSI = as.character(POS),
         ref_dp = as.double(ref_d),
         alt_dp = as.double(alt_d),
         locus_dp = ref_dp + alt_dp,
         IND_IDs = str_sub(ID_ind, end = -5)) %>% 
  select(IND_IDs, POSI, REF, ALT, qual, locus_dp, ref_dp, alt_dp) %>% 
  filter(locus_dp >= 20, qual >= 30)
head(var_data_2)
dim(var_data_2)

additional_samples <- read_tsv("last_33_SARS_pop.tsv", col_names = F) %>% 
  rename(POS = X2, REF = X4, ALT = X5, qual = X6, format = X8, geno_info = X10) %>%
  separate(X1, into = c("ID_ind", "CHROM"), sep = "\\ ") %>%
  rowwise() %>%
  mutate(allele_depth = map_chr(
    geno_info, function(s) rev(strsplit(s, ":")[[1]])[6])) %>%
  separate(allele_depth, into = c("ref_d", "alt_d"), sep = "\\,") %>%
  rowwise() %>%
  mutate(POSI = as.character(POS),
         ref_dp = as.double(ref_d),
         alt_dp = as.double(alt_d),
         locus_dp = ref_dp + alt_dp,
         IND_IDs = str_sub(ID_ind, end = -5)) %>% 
  select(IND_IDs, POSI, REF, ALT, qual, locus_dp, ref_dp, alt_dp) %>% 
  filter(locus_dp >= 20, qual >= 30)
head(additional_samples)
dim(additional_samples)

whole_pop_sars <- bind_rows(var_data_2, additional_samples)
dim(whole_pop_sars)

uniquee_sites <- whole_pop_sars %>% 
  distinct(IND_IDs)
dim(uniquee_sites)

sample_IDs <- whole_pop_sars %>% 
  distinct(IND_IDs) %>%
  rowwise() %>% 
  mutate(IDS_sampls_sars = str_sub(
    IND_IDs, start = 9))
dim(sample_IDs)
head(sample_IDs)


pop_freq_var <- whole_pop_sars %>% 
  group_by(POSI) %>%
  summarise(n = n(),
            freq = n / 85)
head(pop_freq_var)
dim(pop_freq_var)
testito <- pop_freq_var %>% 
  distinct(POSI)
dim(testito)

#plots for distributions of genotype frequencies and median alt all freq


freq_whole <- inner_join(pop_freq_var, whole_pop_sars, by = "POSI") %>%
  rowwise() %>% 
  mutate(IDS_sampls_sars = str_sub(
    IND_IDs, start = 9),
    biallelic = if_else(str_detect(ALT, ","), "no", "yes")) %>% 
  group_by(POSI) %>% 
  mutate(alt_prop = alt_dp/locus_dp,
         median_alt_dp = median(alt_prop),
         median_locus_dp = median(locus_dp),
         phys_pos = as.double(POSI)) %>% 
  filter(freq >= 0.05, biallelic == "yes") %>% 
  arrange(phys_pos)
head(freq_whole)
dim(freq_whole)

uniquee_sites <- freq_whole %>% 
  distinct(POSI)
dim(uniquee_sites)
uniquee_samples <- freq_whole %>%
  ungroup() %>% 
  distinct(IDS_sampls_sars)
dim(uniquee_samples)


indels_sars <- freq_whole %>%
  ungroup() %>% 
  mutate(posit_len = str_length(ALT)) %>%
  select(posit_len, POSI, ALT) %>%
  distinct(POSI, posit_len, ALT)
head(indels_sars)
dim(indels_sars)

summary(indels_sars)


  distinct(POSI)
dim(testito_2)
head(testito_2)

testito_2 <- freq_whole %>%
  ungroup() %>% 
  distinct(POSI)
dim(testito_2)
head(testito_2)

metadata_sars <- read_tsv("sars_sample_ID_meta.tsv", col_names = T)
head(metadata_sars)
dim(metadata_sars)

final_whole_mat <- inner_join(metadata_sars, freq_whole, by = "IDS_sampls_sars") %>% 
  mutate(Symp = if_else(Symptoms == "SI", "Symptomatic", "Asymptomatic"))
head(final_whole_mat)
dim(final_whole_mat)
testito_3 <- final_whole_mat %>%
  ungroup() %>% 
  distinct(IDS_sampls_sars)
dim(testito_3)
head(testito_3)
write_tsv(testito_3, "final_lib_IDs_pop_gen.tsv")

location_samps <- final_whole_mat %>%
  distinct(IND_IDs, location) %>% 
  group_by(location) %>%
  summarise(n = n())
location_samps

freq_symp <- final_whole_mat %>%
  group_by(phys_pos, location) %>%
  summarise(n = n()) %>%
  mutate(freq = case_when(
    location == "GTO" ~ n/43,
    location == "SLP" ~ n/31,
    location == "JAL" ~ n/9)) %>%
  filter(freq >= 0.25) %>% 
  ggplot(aes(location, reorder(phys_pos, desc(phys_pos)), fill = freq)) +
  geom_tile() +
  ylab("Genomic \n position") +
  scale_fill_viridis_c(direction = -1, end = 1/2, alpha = 2/3) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 0),
        axis.text.x=element_text(size = 10, angle = 45, vjust = 1/2),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        axis.line = element_line())
freq_symp

head(freq_symp)

fixes_poly_sites <- final_whole_mat %>% 
  select(IDS_sampls_sars, POSI, freq) %>% 
  filter(freq >= 0.25) %>% 
  distinct(POSI)
head(fixes_poly_sites)
dim(fixes_poly_sites)

cov_mat_meta <- right_join(testito_3, final_whole_mat, by = "IDS_sampls_sars") %>% 
  distinct(IDS_sampls_sars, .keep_all = TRUE) %>% 
  select(IDS_sampls_sars, location, Symp)
head(cov_mat_meta)
dim(cov_mat_meta)


#test the discrete association of variants with the expression of symptoms
#
#
#
#
##waffle plot 
##
waff_data <- final_whole_mat %>% 
  waffle_iron(aes_d(group = phys_pos))

ggplot(waff_data, aes(x, y, fill = group)) + 
  geom_waffle()
  
##
##
##
##
disc_var_eff <- final_whole_mat %>%
  mutate(fisical_pos = as.factor(phys_pos)) %>% 
  filter(location == "GTO") %>% 
  group_by(fisical_pos, Symp) %>%
  summarise(n = n()) %>%
  filter(n() > 1) %>%
  mutate(testillo = rstatix::binom_test(n),
         dummy_var = if_else(fisical_pos > 1, "yes", "no")) %>% 
  filter(testillo$n >= 15) %>% 
  ggplot(aes(dummy_var, reorder(fisical_pos, desc(fisical_pos)),
             size = testillo$n, color = -log10(testillo$p))) +
  geom_point() +
  ylab("Genomic position") +
  scale_color_viridis_c(end = 1/2, direction = -1) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 0),
        axis.text.x=element_text(size = 0),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        axis.line = element_line())
disc_var_eff
    
head(disc_var_eff)


symp_gto <- final_whole_mat %>% 
  filter(location == "GTO",
         Symp == "Asymptomatic",
         median_alt_dp >= 0.8) %>% 
  distinct(POSI)
dim(symp_gto)

props_synt_gto <- final_whole_mat %>% 
  filter(location == "GTO", POSI == "28854") %>% 
  group_by(Symp) %>% 
  summarise(n_1 = n())
props_synt_gto

write_tsv(final_whole_mat, "SARS_central_mex_variants.tsv")


symp_test <- final_whole_mat %>%
  filter(n >= 25, location == "GTO") %>% 
  group_by(POSI) %>% 
  rstatix::anova_test(alt_prop ~ Symp)
symp_test

symp_wilc <- final_whole_mat %>%
  filter(location == "GTO") %>%
  ungroup() %>% 
  #group_by(POSI) %>% 
  rstatix::wilcox_test(alt_prop ~ Symp)
symp_wilc

eff_size_symp <- final_whole_mat %>%
  filter(location == "GTO") %>% 
  ungroup() %>% 
  rstatix::wilcox_effsize(alt_prop ~ Symp)
eff_size_symp


median_symp <- final_whole_mat %>%
  filter(location == "GTO") %>%
  group_by(Symp) %>% 
  summarise(mediana_s = median(alt_prop))
median_symp
  


symp_plot <- final_whole_mat %>%
  filter(location == "GTO") %>% 
  ggplot(aes(Symp, alt_prop, fill = Symp)) +
  geom_violin(color = NA) +
  geom_boxplot(width=0.075, lwd = 0.5, 
               fill = "#ffffff", outlier.colour = NA) +
  scale_fill_viridis_d(alpha = 2/3, end = 1/2, direction = -1) +
  coord_flip() +
  labs(y = "Genome-wide \n Alternative allele proportion") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.line = element_line(),
        legend.position = "none")
symp_plot

symp_locus <- final_whole_mat %>% 
  filter(location == "GTO", POSI == "28854") %>%
  ggplot(aes(Symp, alt_prop, fill = Symp)) +
  geom_violin(color = NA) +
  geom_boxplot(width=0.075, lwd = 0.5, 
               fill = "#ffffff", outlier.colour = NA) +
  scale_fill_viridis_d(alpha = 2/3, end = 1/2, direction = -1) +
  coord_flip() +
  labs(y = "Alternative allele proportion \n at locus S194L") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.line = element_line(),
        legend.position = "none")
symp_locus
  


order_plot <- c("JAL", "GTO", "SLP")
ridges_sars_freq <- final_whole_mat %>%
  group_by(phys_pos, location) %>%
  summarise(n = n()) %>%
  mutate(frequencies = case_when(
    location == "GTO" ~ n/43,
    location == "SLP" ~ n/31,
    location == "JAL" ~ n/9)) %>% 
  ggplot(aes(frequencies, factor(location, rev(order_plot)), fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis_c(name = "Frequency", begin = 1/6, end = 1/2, direction = -1) +
  labs(x = "Polymorphic site frequency", y = "Location") +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.line = element_line(),
        legend.position = "none")
ridges_sars_freq



ridges_sars <- final_whole_mat %>%
  ggplot(aes(alt_prop, factor(location, rev(order_plot)), fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3) +
  scale_fill_viridis_c(name = "Frequency", begin = 1/6, end = 1/2, direction = -1) +
  labs(x = "Alternative allele proportion", y = "Location") +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.line = element_line(),
        legend.position = "none")
ridges_sars



corr_net <- final_whole_mat %>%
  ungroup() %>% 
  select(freq, qual:alt_dp, alt_prop) %>% 
  corrr::correlate() %>% 
  corrr::network_plot(min_cor = 0.05) +
  scale_colour_viridis_c(end = 1/2)
corr_net

relevant_variants <- c("")



vars_plot <- final_whole_mat %>% 
  ggplot(aes(x = phys_pos, y = freq, color = median_alt_dp)) +
  geom_segment(aes(x = phys_pos, xend = phys_pos,
                   y = 0, yend = freq), size = 3/4) +
  geom_point(size = 6, alpha = 1/10) +
  scale_color_viridis_c(begin = 1/6, end = 1/2, direction = -1) +
  scale_x_continuous(breaks = seq(from = 0, to = 30000, by = 5000)) +
  labs(x = "Genomic position", y = "Polymorphic site \n frequency") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        axis.line = element_line())
vars_plot

fig_1_up <- cowplot::plot_grid(
  ridges_sars_freq, ridges_sars, symp_plot, 
  rel_widths = c(1, 1, 3/4),
  ncol = 3, labels = NA,
  label_size = 12)
fig_1_up

gene_models_sars <- read_tsv("visibleData.bed", col_names = T) %>% 
  select(X1:X4) %>%
  rowwise() %>% 
  mutate(y_start = if_else(X1 == "NC_045512.2", "-Inf", "No"),
         y_end = if_else(X1 == "NC_045512.2", "Inf", "No"))
head(gene_models_sars)

gene_models_draw <- ggplot() +
  geom_rect(data = gene_models_sars,
            mapping = aes(ymin = y_start, ymax = y_end,
                          xmin = X2, xmax = X3, fill = X4),
            alpha = 1/3, inherit.aes = FALSE, color = "#2f3542", size = 1/2) +
  scale_fill_viridis_d() +
  theme_void() +
  theme(legend.position="none")
gene_models_draw

fig_1_whole <- cowplot::plot_grid(
  vars_plot, gene_models_draw,
  rel_heights = c(1,1/6),
  nrow = 2, labels = NA, align = "H")
fig_1_whole

fig_1_up <- cowplot::plot_grid(
  ridges_sars_freq, ridges_sars, symp_plot, 
  rel_widths = c(1, 1, 3/4),
  ncol = 3, labels = NA,
  label_size = 12)
fig_1_up





summary(vars_plot)


#load admix data
#
#load bam list as data frame
bam_names <- read_tsv("bam.filelist", col_names = F) %>% 
  mutate(IND_IDs = str_sub(X1, end = -5),
         IDS_sampls_sars = str_sub(
           IND_IDs, start = 9)) %>% 
  select(IDS_sampls_sars)
head(bam_names)
dim(bam_names)
bams_to_remove <- full_join(testeable_bams, bam_names, by = "IDS_sampls_sars") %>%
  filter(is.na(IND_IDs))
head(bams_to_remove)
dim(bams_to_remove)

testeable_bams <- sample_IDs %>%
  mutate(IDS_sampls_sars = str_sub(
         IND_IDs, start = 9))
head(testeable_bams)
dim(testeable_bams)

whole_meta_bam <- inner_join(testeable_bams, bams_meta, by = "IDS_sampls_sars")
head(whole_meta_bam)
dim(whole_meta_bam)

head(metadata_sars)
dim(metadata_sars)

#load covariance matrix and set row and colnames
#
#lib_names <- bam_list$id_geno
head(bam_names)




lib_names <- bam_names$IDS_sampls_sars
head(lib_names)
dim(lib_names)
name <- "SARS_ANGSD.covMat"
m <- as.matrix(read.table(name))
rownames(m) <- lib_names
colnames(m) <- lib_names

aa <- as.data.frame(colnames(m) %in% bams_to_remove$IDS_sampls_sars)
aa
head(aa)
AA <- aa %>%
  mutate(sampID = rownames(aa)) %>% 
  filter(`colnames(m) %in% bams_to_remove$IDS_sampls_sars` == "TRUE")
AA
head(AA)
dim(AA)
test <- c(AA$sampID)
test

cov_mat <- m[-c(1,3,5,10,11,14,16,17,19,21,
               22,23,24,26,28,29,31,72,100),
             -c(1,3,5,10,11,14,16,17,19,21,
               22,23,24,26,28,29,31,72,100)]

par(mar=c(2,7,3,4)+.1)
heat_map_sars <- pheatmap(cov_mat,
                         border_color = F,
                         color = viridis(1000),
                         cellwidth = 6, cellheight = 6,
                         cluster_rows = T,
                         cluster_cols = T, scale = "row",
                         show_colnames = F, show_rownames = F)
heat_map_sars


PCA_SARS <- prcomp(cov_mat)
head(PCA_SARS)
summary(PCA_SARS)
PCadmix<-data.frame(PCA_SARS$x, Location = cov_mat_meta$location, symptoms = cov_mat_meta$Symp)
head(PCadmix)

anova_pca <- PCadmix %>% 
  select(PC1, Location, symptoms) %>% 
  rstatix::anova_test(PC1 ~ symptoms)
anova_pca

PCA_plot_cool <- PCadmix %>%
  ggplot(aes(PC1, PC2, color = Location)) +
  geom_point(size=6) +
  scale_color_viridis_d(1) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 8),
        axis.text.x=element_text(size = 0),
        axis.text.y=element_text(size = 8),
        axis.title.y=element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        axis.line = element_line())
PCA_plot_cool

longer_PCA <- PCadmix %>%
  select(Location, symptoms, PC1:PC3) %>%
  pivot_longer(cols = -c(Location, symptoms), names_to = "PC", values_to = "values")
head(longer_PCA)

violins_PC_loc <- longer_PCA %>% 
  ggplot(aes(factor(Location, rev(order_plot)), values, fill = Location)) +
  geom_violin() +
  geom_boxplot(width=0.1, lwd = 0.5, fill = "#ffffff", outlier.colour = NA) +
  coord_flip() +
  labs(y = "Score") +
  facet_wrap(~ PC, ncol = 1) +
  scale_fill_viridis_d(alpha = 2/3, end = 1/2) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.line = element_line(),
        legend.position = "none")
violins_PC_loc

violins_PC_sympt <- longer_PCA %>%
  filter(Location == "GTO") %>% 
  ggplot(aes(symptoms, values, fill = symptoms)) +
  geom_violin() +
  geom_boxplot(width=0.1, lwd = 0.5, fill = "#ffffff", outlier.colour = NA) +
  coord_flip() +
  labs(y = "Score") +
  facet_wrap(~ PC, ncol = 1) +
  scale_fill_viridis_d(alpha = 2/3, end = 1/2, direction = -1) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 10),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.y=element_text(size = 0),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        axis.line = element_line(),
        legend.position = "none")
violins_PC_sympt

head(violins_PC)
dim(violins_PC)
fig_1_whole

fig_1_1 <- cowplot::plot_grid(
  vars_plot, gene_models_draw,
  ncol = 1, rel_heights = c(5/6, 1/6), labels = "A")
fig_1_1

dt_plot_ar <- cowplot::plot_grid(
  disc_var_eff, NULL, ncol = 2,
  rel_widths = c(1/2, 1/2),
  labels = c("C", ""))
dt_plot_ar

fig_1_3_up_nor <- cowplot::plot_grid(
  ridges_sars_freq, dt_plot_ar,
  nrow = 2, labels = c("C", "D"),
  rel_heights = c(1/4, 3/4))
fig_1_3_up_nor

fig_1_3_up_nor_almost <- cowplot::plot_grid(
  symp_plot, symp_locus,
  nrow = 2, labels = c("E", "F"))
fig_1_3_up_nor_almost

fig_1_4 <- cowplot::plot_grid(
  freq_symp, NULL, fig_1_3_up_nor, NULL, fig_1_3_up_nor_almost,
  ncol = 5, labels = c("B", "", "", "", ""), rel_widths = c(2/10, 1/10, 3/10, 1/10, 3/10))
fig_1_4

fig_1_6 <- cowplot::plot_grid(
  fig_1_1, fig_1_4,
  ncol = 1, rel_heights = c(0.35, 0.65))
fig_1_6

fig_supp <- cowplot::plot_grid(
  violins_PC_loc, ridges_sars,
  ncol = 1, rel_heights = c(3/4, 1/4), labels = c("A", "B"))
fig_supp








#load admixture data for covid
#
admix_sars <- read_tsv("cluster_SARS.qopt", col_names = F) %>% 
  rename(pops = X1) %>% 
  separate(pops, into = c("pop_A", "pop_B", "pop_C"), sep = "\\ ")
head(admix_sars)
pop_admix_sars <- bind_cols(bam_names, admix_sars) %>% 
  filter(!IND_IDs %in% c("trimmed_LA_FB_05_0348_S60", "trimmed_UG_FR_P_0002_A7_S79"))
head(pop_admix_sars)
dim(pop_admix_sars)



SARS_admix <- pop_admix_sars %>% 
  rowwise() %>% 
  mutate(A_pop = as.double(pop_A),
         B_pop = as.double(pop_B),
         C_pop = as.double(pop_C)) %>%
  select(IND_IDs, A_pop:C_pop) %>%
  pivot_longer(A_pop:C_pop, names_to = "population", values_to = "proportion") %>%
  group_by(IND_IDs) %>% 
  mutate(likely_assignment = population[which.max(proportion)],
         assingment_prob = max(proportion)) %>%
  ungroup() %>% 
  arrange(likely_assignment, assingment_prob) %>%
  mutate(id_geno = forcats::fct_inorder(factor(IND_IDs)))
head(SARS_admix)
dim(SARS_admix)
distinct(SARS_admix$IND_IDs)

testillo <- inner_join(pop_admix_sars, SARS_admix, by = "IND_IDs") %>%
  select(IND_IDs, likely_assignment) %>% 
  distinct()
head(testillo)
dim(testillo)

A_pop_sars <- testillo %>% 
  filter(likely_assignment == "A_pop")
head(A_pop_sars)
dim(A_pop_sars)
write_tsv(A_pop_sars, "A_pop_sars.tsv", col_names = F)

B_pop_sars <- testillo %>% 
  filter(likely_assignment == "B_pop")
head(B_pop_sars)
dim(B_pop_sars)
write_tsv(B_pop_sars, "B_pop_sars.tsv", col_names = F)

C_pop_sars <- testillo %>% 
  filter(likely_assignment == "C_pop")
head(C_pop_sars)
dim(C_pop_sars)
write_tsv(C_pop_sars, "C_pop_sars.tsv", col_names = F)

plot_sars_admix <- SARS_admix %>%
  mutate(IDS_sam = str_sub(id_geno, start = 9)) %>% 
  ggplot(aes(IDS_sam, proportion, fill = population)) +
  geom_col(width = 1) +
  scale_fill_viridis_d(end = 2/3, alpha = 2/3) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 8),
        axis.text.x=element_text(size = 8, angle = 45),
        axis.text.y=element_text(size = 8),
        axis.title.y=element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        axis.line = element_line())
plot_sars_admix

summary(PCA_SARS)


PCadmix<-data.frame(PCA_SARS$x, Genotype = testillo$likely_assignment)
head(PCadmix)





PCA_plot_cool <- PCadmix %>%
  ggplot(aes(PC1, PC2, color = Genotype)) +
  geom_point(size=4) +
  scale_color_viridis_d(end = 2/3, alpha = 0.9) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x=element_text(size = 8),
        axis.text.x=element_text(size = 0),
        axis.text.y=element_text(size = 8),
        axis.title.y=element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        axis.line = element_line())
PCA_plot_cool

fig_2_final <- cowplot::plot_grid(PCA_plot_cool, NULL,
                                  ncol = 2, rel_widths = c(0.4, 0.6),
                                  labels = c("B", "C"))
fig_2_final

fig_2_final_2 <- cowplot::plot_grid(plot_sars_admix, fig_2_final,
                                    ncol = 1, rel_heights = c(2/3, 1),
                                    labels = c("A", ""))
fig_2_final_2






# comparison of distributions between asymptomatic and symptomatic


bam_cts <- read_tsv("bam_conts.tsv", col_names = F)
head(bam_cts)
dim(bam_cts)
meta_counts <- bind_cols(bam_names, bam_cts) %>% 
  rename(number_m_reads = X1)
head(meta_counts)
tail(meta_counts)

supp_3_meta <- inner_join(meta_counts, final_whole_mat, by = "IDS_sampls_sars")
dim(supp_3_meta)
head(supp_3_meta)

head(final_whole_mat)

write_tsv(supp_3_meta, "S3_pop_meta_cts.tsv")

summary(bam_cts)
%>%
  ggplot(aes(sqrt(X1))) +
  geom_density()
bam_cts
head(bam_cts)








