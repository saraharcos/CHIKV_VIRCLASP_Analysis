##Analyzing MaxQuant output
library(tidyverse)
library(here)
library(GGally)


mq_full <- read_tsv(here("MaxQuant", "proteinGroups.txt")) %>%
  filter(!grepl("REV_",`Protein IDs`) & 
        !grepl("CON_", `Protein IDs`)) %>%
  separate(`Majority protein IDs`, into = c("Trash", "Majority protein IDs", "ID"), sep = "\\|", extra = "drop") 

mq_coverage <- mq_full %>%
  select(ID, `Majority protein IDs`,contains("Coverage"), -contains("Razor"), -`Unique sequence coverage [%]`, -`Sequence coverage [%]`) %>%
  gather(Sample, Coverage, -`Majority protein IDs`) %>%
  mutate(Sample = gsub("no4SU", "Ctr", Sample),
         Sample = gsub("Sequence coverage ", "", Sample),
         Sample = gsub(" [//[]%", "", Sample),
         Sample = gsub("]", "", Sample))

mq_spectra <- mq_full %>%
  select(ID, `Majority protein IDs`, contains("MS/MS count ")) %>%
  gather(Sample, Spectra, -`Majority protein IDs`, -ID) %>%
  mutate(Sample = gsub("no4SU", "Ctr", Sample))

mq_intensities <- mq_full %>%
  select(ID, `Majority protein IDs`, contains("Intensity "), -contains("LFQ"))  %>%
  gather(Sample, Intensity, -`Majority protein IDs`, -ID) %>%
  mutate(Sample = gsub("no4SU", "Ctr", Sample))


#Spectral kernel density
ggplot(mq_spectra, aes(x = log2(Spectra), group = Sample, color = Sample)) +
  geom_density(alpha = 0.6) +
  scale_x_continuous(limits = c(-2, 10))

#Intensity kernel density
ggplot(mq_intensities, aes(x = log2(Intensity), group = Sample, color = Sample)) +
  geom_density(alpha = 0.6)

#Separate minus and plus, then add group based on control, not control
#Intensities
mq_intensities_minus <- mq_intensities %>%
  filter(grepl("minus", Sample)) %>%
  mutate(Sample = gsub("Intensity ", "", Sample),
         Group = substr(Sample, start = 1, stop = 3))

ggplot(mq_intensities_minus, aes(x = log2(Intensity), group = Sample, color = Group)) +
  geom_density(size = 1) +
  scale_x_continuous(limits = c(15, 38)) +
  theme_minimal() +
  labs(title = "Raw Intensities: Minus IFN") +
  theme(text = element_text(size = 20))

ggsave(filename = here("Figures", "intensities_minus_kde.pdf"), width = 8, height = 4)

ggplot(mq_intensities_minus, aes(y = log2(Intensity), x = Sample, color = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Raw Intensities: Minus IFN") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust= 0.5))

ggsave(filename = here("Figures", "intensities_minus_box.pdf"), width = 8, height = 4)

mq_intensities_plus <- mq_intensities %>%
  filter(grepl("plus", Sample)) %>%
  mutate(Sample = gsub("Intensity ", "", Sample),
         Group = substr(Sample, start = 1, stop = 3))

ggplot(mq_intensities_plus, aes(x = log2(Intensity), group = Sample, color = Group)) +
  geom_density(size = 1) +
  scale_x_continuous(limits = c(15, 38)) +
  theme_minimal() +
  labs(title = "Raw Intensities: Plus IFN") +
  theme(text = element_text(size = 20))

ggsave(filename = here("Figures", "intensities_plus_kde.pdf"), width = 8, height = 4)

ggplot(mq_intensities_plus, aes(y = log2(Intensity), x = Sample, color = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Raw Intensities: Plus IFN") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust= 0.5))

ggsave(filename = here("Figures", "intensities_plus_box.pdf"), width = 8, height = 4)

#Spectra scatterplots
mq_spectra_minus <- mq_spectra %>%
  filter(grepl("minus", Sample)) %>%
  mutate(Sample = gsub("MS/MS count ", "", Sample))%>%
  spread(Sample, Spectra, fill = 0) %>%
  filter(ID != "AHNK_HUMAN")

#library(GGally)
ggpairs(mq_spectra_minus[3:13], aes(alpha=0.4)) +
  theme_minimal()

#ggsave(filename = here("Figures", "spectra_minus.pdf"), width = 10, height = 10)

mq_spectra_plus <- mq_spectra %>%
  filter(grepl("plus", Sample)) %>%
  mutate(Sample = gsub("MS/MS count ", "", Sample))%>%
  spread(Sample, Spectra, fill = 0) %>%
  filter(ID != "AHNK_HUMAN")

ggpairs(mq_spectra_plus[3:13], aes(alpha=0.4)) +
  theme_minimal()

#ggsave(filename = here("Figures", "spectra_plus.pdf"), width = 10, height = 10)

#2-fold enrichment dataset
  
minus_enriched <- mq_intensities_minus %>%
  full_join(mq_coverage, by = c("Majority protein IDs", "Sample")) %>%
  mutate(Exists = 1) %>%
  filter(Intensity > 0) %>%
  #filter(Coverage >= 15 | Group == "Ctr") %>%
  select(-Coverage) %>%
  group_by(`Majority protein IDs`, Group) %>%
  mutate(Exists = sum(Exists)) %>%
  filter(Exists == 3 | Group == "Ctr") %>%
  select(-Exists) %>% 
  group_by(`Majority protein IDs`, Group) %>%
  mutate(ave_intensity = mean(Intensity)) %>%
  select(-Sample, -Intensity) %>%
  unique() %>%
  spread(Group, ave_intensity, fill = 0) %>%
  gather(Group, ave_intensity, -Ctr, -`Majority protein IDs`, -ID) %>%
  filter(ave_intensity > 0) %>%
  mutate(Enrichment = case_when(Ctr == 0 ~ ave_intensity,
                                TRUE ~ ave_intensity / Ctr)) %>%
  filter(Enrichment >= 5)
  
plus_enriched <- mq_intensities_plus %>%
  full_join(mq_coverage, by = c("Majority protein IDs", "Sample")) %>%
  mutate(Exists = 1) %>%
  filter(Intensity > 0) %>%
  #filter(Coverage >= 15 | Group == "Ctr") %>%
  select(-Coverage) %>%
  group_by(`Majority protein IDs`, Group) %>%
  mutate(Exists = sum(Exists)) %>%
  filter(Exists == 3 | Group == "Ctr") %>%
  select(-Exists) %>% 
  group_by(`Majority protein IDs`, Group) %>%
  mutate(ave_intensity = mean(Intensity)) %>%
  select(-Sample, -Intensity) %>%
  unique() %>%
  spread(Group, ave_intensity, fill = 0) %>%
  gather(Group, ave_intensity, -Ctr, -`Majority protein IDs`, -ID) %>%
  filter(ave_intensity > 0) %>%
  mutate(Enrichment = case_when(Ctr == 0 ~ ave_intensity,
                                TRUE ~ ave_intensity / Ctr)) %>%
  filter(Enrichment >= 5)
  
test <- tibble(name = unique(c(plus_enriched$ID, minus_enriched$ID))) %>%
  filter(!(name %in% crapome$UNIPROT_ID))


filter3_notest <- filter3 %>%
  filter(!(ID %in% test$name))

n_distinct(filter3_notest$Accession)

plus_ctr_0 <- plus_enriched %>% filter(Ctr == 0)
plus_ctr_val <- plus_enriched %>% filter(Ctr > 0)
minus_ctr_0 <- minus_enriched %>% filter(Ctr == 0)
minus_ctr_val <- minus_enriched %>% filter(Ctr > 0)


plus_ctr <- mq_intensities_plus %>%
  mutate(Exists = 1) %>%
  filter(Intensity > 0) %>%
  group_by(`Majority protein IDs`, Group) %>%
  mutate(Exists = sum(Exists)) %>%
  filter(Group == "Ctr" & Exists == 1)

plus_one_ctr <- plus_enriched %>%
  filter(ID %in% plus_ctr$ID)

plus_two_ctr <- plus_enriched %>%
  filter(!(ID %in% plus_ctr$ID))
    
minus_ctr <- mq_intensities_minus %>%
  mutate(Exists = 1) %>%
  filter(Intensity > 0) %>%
  group_by(`Majority protein IDs`, Group) %>%
  mutate(Exists = sum(Exists)) %>%
  filter(Group == "Ctr" & Exists == 1)

minus_one_ctr <- minus_enriched %>%
  filter(ID %in% plus_ctr$ID)

minus_two_ctr <- minus_enriched %>%
  filter(!(ID %in% plus_ctr$ID))



    
    
    



