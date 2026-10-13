setwd("/home/predatoma/Desktop/PBE/PF_collection/Classification/")
setwd("/media/paloma/TFM_Paloma/RyC_Collection/Classification/")

library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)
library(readr)
library(ggplot2)
library(tibble)
library(paletteer)
library(viridis)
library(cowplot)
library(GenomicRanges)



strain_order <- c("EC24", "EC06", "EC23", "EC25", "EC05", "EC14", "EC09", "EC07",
                  "EC02", "EC01", "EC19", "EC11", "EC10", "EC16", "EC21", "EC13",
                  "EC18", "EC22", "EC03", "EC08", "EC12", "EC15", "EC17", "EC04",
                  "EC20", "KPN04", "KPN03", "KPN05", "KPN01", "KPN07", "KPN10",
                  "KPN11", "KPN12", "KPN02", "KPN09", "KPN13", "KPN17", "KPN06",
                  "KPN15", "KPN18", "KPN08", "KPN16", "KPN19", "KPN14", "KPN20",
                  "KV01", "KQ03", "KQ01", "KQ02", "KQ04")


# FIG2
#Heatmap resistome per strain 
#RESISTANCE INFORMATION --> RESFINDER
resistance <- read.csv("Resistance_combined.csv", sep= "\t")
resistance$Strain <- sub("^PF_", "", resistance$Strain)
resistance <- resistance %>%
  separate_rows(RESISTANCE, sep = ";")

res_atb <- resistance %>%
  dplyr::slice(-1) %>%                 
  group_by(Strain, RESISTANCE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = rev(strain_order))) %>%
  ungroup()

family_resistance <- read.csv("ATB_family.csv", sep = ",") 
res_atb <- full_join(res_atb, family_resistance, by = "RESISTANCE")

res_atb <- res_atb%>%
  mutate(RESISTANCE = case_when(
    RESISTANCE == "" ~ "XXX",
    TRUE ~ RESISTANCE  
  ))
res_atb <- res_atb %>%
  arrange(Family, RESISTANCE) %>%
  mutate(RESISTANCE = factor(RESISTANCE, levels = unique(RESISTANCE)))

#Per species
ecoli <- res_atb[grepl("^EC", res_atb$Strain), ]

ecoli_resistance<- ggplot(ecoli, aes(x = RESISTANCE, y = Strain, fill = rep_count)) +
  geom_tile(color = "#5b176b", alpha = 0.6) +
  scale_fill_gradient(low = "#ffe5ec", high = "#9d3d7d", name = "Count") +
  theme_bw() +
  ylab("Predicted PTU") +
  xlab("Strain") +
  theme_bw() +
  theme (
    axis.text.x = element_text(angle = 90, hjust = 1, size = 17),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename="panel_coli_pang_resistome.svg", plot = ecoli_resistance, dpi =350)

klebsiella <- res_atb[grepl("^K", res_atb$Strain), ]

klebsiella_resistance <- ggplot(klebsiella, aes(x = RESISTANCE, y = Strain, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  ylab("Predicted PTU") +
  xlab("Strain") +
  theme_bw() +
  theme (
    axis.text.x = element_text(angle = 90, hjust = 1, size = 17),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename="panel_kleb_pang_resistome.svg", plot = klebsiella_resistance, dpi =350)


#we differenciate between AMR genes in plasmids vs in chromosomal genes --> Fig. Supp. 1
chromosome <- c("1", "contig_1_circular", "NODE_1_length_5105783_cov_0.998618",
                "contig_1", "NODE_1+_length_5132816_cov_1", "NODE_1+_length_5116136_cov_1",
                "NODE_1+_length_4664359_cov_1", "NODE_1+_length_5494081_cov_1",
                "NODE_1+_length_235159_cov_1", "NODE_edge_1+_length_5140923_cov_19",
                "NODE_1+_length_5246450_cov_1","NODE_1+_length_5059491_cov_1",
                "NODE_1+_length_5613265_cov_1"
)


resistome_mobilome <- resistance[!resistance$SEQUENCE %in% chromosome, ]
resistome_chromsome <- resistance[resistance$SEQUENCE %in% chromosome,]

resistome_mob1<- resistome_mobilome %>%
  dplyr::slice(-1) %>%                 
  group_by(Strain, RESISTANCE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = strain_order)) %>%
  ungroup()

resistome_mob1 <- full_join(resistome_mob1, family_resistance, by = "RESISTANCE")
resistome_mob1 <- resistome_mob1 %>%
  arrange(Family, RESISTANCE) %>% 
  mutate(RESISTANCE = factor(RESISTANCE, levels = rev(unique(RESISTANCE))))
resistome_mob1 <- resistome_mob1[!is.na(resistome_mob1$Strain), ]


resistance_plasm <- ggplot(resistome_mob1, aes(x = Strain, y = RESISTANCE, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

ggsave(filename="supp1_resistome_plasmidome.svg", plot = resistance_plasm, dpi =350)

resistome_chr1 <- resistome_chromsome %>%
  dplyr::slice(-1) %>%                 
  group_by(Strain, RESISTANCE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = strain_order)) %>%
  ungroup()
resistome_chr1 <- full_join(resistome_chr1, family_resistance, by = "RESISTANCE")
resistome_chr1 <- resistome_chr1 %>%
  arrange(Family, RESISTANCE) %>% 
  mutate(RESISTANCE = factor(RESISTANCE, levels = rev(unique(RESISTANCE))))
resistome_chr1 <- resistome_chr1[!is.na(resistome_chr1$Strain), ]

resistance_chr <- ggplot(resistome_chr1, aes(x = Strain, y = RESISTANCE, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

ggsave(filename="supp1_resistome_chromosome.svg", plot = resistance_chr, dpi =350)



#Supp. Fig. 2--> plotting the AMR gene
res_atb <- read.csv("Resistance_combined.csv", sep= "\t")
res_atb$Strain <- sub("^PF_", "", res_atb$Strain)

res_atb <- res_atb %>%  
  group_by(Strain, GENE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = rev(strain_order))) %>%
  ungroup()



plasma2 <- ggplot(res_atb, aes(x = GENE, y = Strain, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  scale_x_discrete(labels = function(x) {
    lapply(x, function(z) bquote(italic(.(z))))
  }) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename="genes_resistance.svg", plot =plasma2, dpi =350)



resistome_mobilome <- resistance[!resistance$SEQUENCE %in% chromosome, ]
res_mob <- resistome_mobilome %>%  
  group_by(Strain, GENE, RESISTANCE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = rev(strain_order))) %>%
  ungroup()

resistance_mob<- ggplot(res_mob, aes(x = Strain, y = GENE, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 



resistome_chromsome <- resistance[resistance$SEQUENCE %in% chromosome,]
res_chr <- resistome_chromsome %>%  
  group_by(Strain, GENE, RESISTANCE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = rev(strain_order))) %>%
  ungroup()

resistance_chr<- ggplot(res_chr, aes(x = Strain, y = GENE, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 


p2 <- resistance_chr

p1 <- resistance_mob + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank())


common_legend <- get_legend(resistance_chr + theme(legend.box.just = "left"))
resistance_location <- plot_grid(p1, p2, 
                        ncol = 1, align = "v", 
                        rel_heights = c( 1, 1)) 


#FIG3

###FIG3

# First plot --> number of plasmids and predicted mobilization 
#We delete the chromosome contig (first 1, count the )
#As frequency and as counts 
mobilization <- read.csv ("All_MOB_typer.csv", sep = "\t")
mobilization <- mobilization %>%
  group_by(Strain) %>%
  dplyr::slice(-1) %>%
  ungroup()

mob_colors <- c(
  conjugative       = "#a773b6ff",
  mobilizable       = "#e6d0a5ff",
  `non-mobilizable` = "#a4d0bdff")

# summing up number of mobilization items. 

mobilization_count <- mobilization%>%
  #dplyr::slice(-1) %>%                 
  group_by(Strain, predicted_mobility) %>%
  summarize(rep_count = n(), .groups = "drop") %>%  
  ungroup() %>%
  mutate(Strain = factor(Strain, levels = strain_order)) %>%
  arrange(Strain)

mobilization_count <- mobilization%>%
  #dplyr::slice(-1) %>%                 
  group_by(predicted_mobility) %>%
  summarize(rep_count = n(), .groups = "drop") %>%  
  ungroup() %>%
  mutate(Strain = factor(Strain, levels = strain_order)) %>%
  arrange(Strain)

counts <- ggplot(mobilization_count, aes(x = Strain, y = rep_count,
                                         fill = predicted_mobility)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = mob_colors) +
  scale_y_continuous(breaks = scales::pretty_breaks(),
                     labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  ylab("Counts") +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, size = 15),
    axis.text.x = element_text(angle = 90, size = 15, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave(filename="count_motility_MOB_typer.svg", plot = counts, dpi= 350)


# Plasmids per species filling with PCN, without ordering by family.
pcn <- read.csv("PCN.csv", sep = "\t")
pcn <- pcn %>%
  mutate(
    median_ratio = as.numeric(sub(",", ".", median_ratio)),
    mean_ratio = as.numeric(sub(",", ".", mean_ratio))
  )
PTU_prediction <- read.csv("combined_ptu_prediction.tsv", sep ="\t")


plasmids <- PTU_prediction %>%
  left_join(
    pcn %>% 
      select(Sample, p_contig, mean_ratio, median_ratio),
    by = c("Strain" = "Sample", "Contig" = "p_contig")
  )


plasmids <- plasmids %>%
  filter(!is.na(median_ratio))

plasmids$Strain <- factor(plasmids$Strain, levels = strain_order)
plasmids$X.Predicted_PTU <- factor(plasmids$X.Predicted_PTU, 
                                   levels = rev(sort(unique(plasmids$X.Predicted_PTU))))

plasmids_count1 <- plasmids %>%
  group_by(Strain) %>%
  summarize(rep_count = n(), .groups = "drop") %>%  
  ungroup() %>%
  arrange(Strain)

PTUS<- ggplot(plasmids, aes(x = Strain, y = X.Predicted_PTU, fill = median_ratio)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  ylab("Predicted PTU") +
  xlab("Strain") +
  theme_bw() +
  theme (
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename="PTU_prediction.svg", plot=PTUS,dpi = 350)

#Unify Fig3
p1 <- counts + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())

p2 <- PTUS 

common_legend <- get_legend(PTUS + theme(legend.box.just = "left"))
final_plot <- plot_grid(p1, p2, 
                        ncol = 1, align = "v", 
                        rel_heights = c(0.3, 1)) 


ggsave(filename="Fig3_v3_LALALAL.svg", plot = final_plot, dpi =350)


#Plasmids isdentified using plasmidfinder
plasfind <- read.csv("plasmids_abricate_plus_PCN.csv")
plasfind_counts <- plasfind %>%
  dplyr::slice(-1) %>%                 
  group_by(Sample, GENE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Sample = factor(Sample, levels = strain_order)) %>%
  ungroup()

plasfind_counts$GENE <- factor(plasfind_counts$GENE, 
                               levels = rev(sort(unique(plasfind_counts$GENE))))


plasfind_plot <- ggplot(plasfind_counts, aes(x = Sample, y = GENE, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Supp.Fig_3.svg", plot = plasfind_plot, dpi = 350)

## FIG 4+
### ISs from ISescan
ISs <- read.csv("ISescan_compilation.tsv", sep = "\t")
ISs$Strain <- sub("^PF_", "", ISs$Strain)

IS <- ISs %>%              
  group_by(Strain, family) %>%
  summarize(rep_counts = n(), .groups = "drop_last") %>%
  ungroup()


IS$Strain <- factor(IS$Strain, levels = strain_order)

IS$family <- factor(IS$family, 
                    levels = rev(sort(unique(IS$family))))

ISS_plot <- ggplot(IS, aes(x = Strain, y = family, fill = rep_counts)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  ylab("Strain") +
  xlab("IS family") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(filename = "ISs_clean.svg", plot = ISS_plot, dpi = 350)


prevalent_families <- IS %>%
  group_by(family) %>%
  summarize(total_rep = sum(rep_counts, na.rm = TRUE)) %>%
  slice_max(total_rep, n = 10) %>%
  pull(family)

IS <- IS %>%
  filter(family %in% prevalent_families)

ISS_prevalent <- ggplot(IS, aes(x = Strain, y = family, fill = rep_counts)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  ylab("IS family") +
  xlab("Strain") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(hjust = 1, size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(filename = "ISs_clean_filtered.svg", plot = ISS_prevalent, dpi = 350)


### Phages from PHASTEST
phages <- read_csv("Phages_Phastest.csv")
phages$Strain <- sub("^PF_", "", phages$Strain)
phages$Strain <- factor(phages$Strain, levels = strain_order)
all_strains <- unique(phages$Strain)

complete_phages <- phages[phages$Completeness=="intact", ]
complete_phages$summ_phages <- sapply(
  strsplit(complete_phages$`Most Common Phage`, "_"),
  function(x) paste(x[2], x[3], sep = "_")
)


phages_counts <- complete_phages %>%              
  group_by(Strain, summ_phages) %>%
  summarize(rep_counts = n(), .groups = "drop_last") %>%
  ungroup()


prevalent_phages <- phages_counts %>%
  group_by(summ_phages) %>%
  summarize(total_rep = sum(rep_counts, na.rm = TRUE)) %>%
  slice_max(total_rep, n = 10) %>%
  pull(summ_phages)

phages_counts <- phages_counts %>%
  filter(summ_phages %in% prevalent_phages)

template <- expand_grid(
  Strain = levels(all_strains),
  summ_phages  = prevalent_phages
)
phages_counts <- template %>%
  left_join(phages_counts, by = c("Strain", "summ_phages"))


phages_counts$summ_phages <- factor(phages_counts$summ_phages, 
                    levels = rev(sort(unique(phages_counts$summ_phages))))


phages_prevalent <- ggplot(phages_counts, aes(x = Strain, y = summ_phages, 
                                              fill = rep_counts)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count", 
                      na.value = "white") +
  theme_bw() +
  ylab("Phage Family") +
  xlab("Strain") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(hjust = 1, size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




### Integrons from Integron_finder
integrons <- read.csv("Integron_compilation.csv", sep = "\t")

integrons_format <- integrons %>%
  pivot_longer(cols = c(CALIN, complete, In0), 
               names_to = "Type", 
               values_to = "count") %>%
  group_by(STRAIN, Type) %>%
  summarise(rep_counts = sum(count, na.rm = TRUE), .groups = "drop") %>%
  rename(Strain = STRAIN) %>%
  mutate(rep_counts = na_if(rep_counts, 0)) 

integrons_format$Strain <- factor(integrons_format$Strain, levels = strain_order)
integrons_format$Type <- factor(integrons_format$Type, 
                                    levels = rev(sort(unique(integrons_format$Type))))

integron_plot <- ggplot(integrons_format, aes(x = Strain, y = Type, 
                                              fill = rep_counts)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count",
                      na.value = "white") +
  theme_bw() +
  ylab("Phage Family") +
  xlab("Strain") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 18),
    axis.text.y = element_text(hjust = 1, size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


p1 <-ISS_prevalent + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank())

p2 <- phages_prevalent + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank())

p3 <- integron_plot


common_legend <- get_legend(integron_plot + theme(legend.box.just = "left"))
final_plot4 <- plot_grid(p1, p2, p3,  
                         ncol = 1, align = "v", 
                         rel_heights = c(0.5, 0.5, 0.4)) 

ggsave(filename="combine_fig4.svg", plot = final_plot4, dpi =350)




###Fig 5
### Defense systems
defense <- read.csv("Defense_compilation.tsv", sep = "\t")

defense$General <- defense$type

ds_types <- c("DS-10", "DS-11", "DS-16", "DS-17", "DS-19", "DS-20", "DS-22", "DS-24", "DS-27", 
              "DS-28", "DS-3", "DS-34", "DS-36", "DS-37", "DS-39", "DS-4", "DS-40", "DS-6", 
              "DS-8", "DS-9")
defense <- defense %>%
  mutate(General = case_when(
    type %in% ds_types ~ "DS_type",
    TRUE ~ General
  ))
rm(ds_types)



abi_types <- c("AbiD", "AbiE", "AbiH", "AbiJ", "AbiP2", "AbiQ", "AbiU")
defense <- defense %>%
  mutate(General = case_when(
    type %in% abi_types ~ "Abi_type",
    TRUE ~ General
  ))
rm(abi_types)


pd_types <- c("PD-Lambda-1","PD-Lambda-3","PD-T4-1","PD-T4-3","PD-T4-5","PD-T4-6",
              "PD-T4-7","PD-T7-1","PD-T7-2","PD-T7-5")
defense <- defense %>%
  mutate(General = case_when(
    type %in% pd_types ~ "PD-T_type",
    TRUE ~ General
  ))
rm(pd_types)

gao_types <- c("Gao_Her", "Gao_Iet", "Gao_Ppl", "Gao_Qat", "Gao_Tmn")
defense <- defense %>%
  mutate(General = case_when(
    type %in% gao_types ~ "Gao_type",
    TRUE ~ General
  ))
rm(gao_types)

defense2  <- defense

#Differeciating between mechanism
defense2 <- defense %>%
  group_by(Strain, General, activity, genes_count) %>%
  summarize(genes_count = sum(genes_count), .groups = "drop_last")%>%
  ungroup()


defense2 <- defense2 %>% 
  mutate(General = factor(General, levels = unique(General[order(activity)])))

antidefense <- defense2[defense2$activity == "Antidefense", ]
antidefense$Strain <- sub("^PF_", "", antidefense$Strain)

antidefense$Strain <- factor(antidefense$Strain, levels = strain_order)
antidefense$General <- factor(antidefense$General, 
                    levels = sort(unique(antidefense$General)))



antidefense_plot <- ggplot(antidefense, aes(x = Strain, y = General, fill = genes_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  xlab("Strain") +
  ylab("Antidefense General Mechanism") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(hjust = 1, size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 


ggsave(filename = "Antidefense_plot.svg", plot = antidefense_plot, dpi = 350)


defense2 <- defense2[defense2$activity == "Defense", ]
defense2$Strain <- sub("^PF_", "", defense2$Strain)

prevalent_mechanisms <- defense2 %>%
  group_by(General) %>%
  summarize(total_rep = sum(genes_count, na.rm = TRUE)) %>%
  slice_max(total_rep, n = 10) %>%
  pull(General)


defense2 <- defense2 %>%
  filter(General %in% prevalent_mechanisms)
defense$Strain <- factor(defense$Strain, levels = strain_order)
defense2$General <- factor(defense2$General, 
                              levels = sort(unique(defense2$General)))


defense_plot <- ggplot(defense2, aes(x = Strain, y = General, fill = genes_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  xlab("Strain") +
  ylab("Defense General Mechanism") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(hjust = 1, size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

ggsave(filename = "Defense_plot_prevalent.svg", plot = defense_plot, dpi = 350)



###Unify figures IS, Defense Systems, Antidefense systems 
p1 <- defense_plot + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank())

p2 <- antidefense_plot 

common_legend <- get_legend(antidefense_plot + theme(legend.box.just = "left"))
final_plot4 <- plot_grid(p1, p2,  
                        ncol = 1, align = "v", 
                        rel_heights = c(0.5, 0.5)) 

ggsave(filename="combine_fig5.svg", plot = final_plot4, dpi =350)


##Supp. FIg3

AMR_Finder <- read.csv("Combined_resistance_AMRFinder.csv", sep= "\t")
AMR_Finder$Strain <- sub("^PF_", "", AMR_Finder$Strain)
res_atb <- AMR_Finder %>%  
  group_by(Strain, Element.symbol) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = rev(strain_order))) %>%
  ungroup()

res_atb$Element.symbol <- factor(res_atb$Element.symbol, 
                           levels = rev(sort(unique(res_atb$Element.symbol))))

AMRFinder_genes_plot <- ggplot(res_atb, aes(x = Element.symbol, y = Strain, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  theme_bw() +
  xlab("Strain") +
  ylab("Defense General Mechanism") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size =10),
    axis.text.y = element_text(hjust = 1, size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

 ggsave(filename="AMRFinder_genes_plot_1.svg", plot = AMRFinder_genes_plot, dpi =350)
 

res_atb <- res_atb %>%  
  group_by(Strain, GENE) %>%
  summarize(rep_count = n(), .groups = "drop_last") %>%
  mutate(Strain = factor(Strain, levels = rev(strain_order))) %>%
  ungroup()



AMRFinder_plot <- ggplot(res_atb, aes(x = GENE, y = Strain, fill = rep_count)) +
  geom_tile(color = "#03045e", alpha = 0.6) +
  scale_fill_gradient(low = "#a9d6e5", high = "#012a4a", name = "Count") +
  scale_x_discrete(labels = function(x) {
    lapply(x, function(z) bquote(italic(.(z))))
  }) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# --- Normalization function: remove trailing _number suffix ---
normalize_gene <- function(x) sub("_[0-9]+$", "", x)

# --- Per-strain comparison ---
strains <- union(unique(AMR_Finder$Strain), unique(resistance$Strain))

comparison_list <- lapply(strains, function(s) {
  amr_genes <- AMR_Finder %>%
    filter(Strain == s) %>%
    pull(Element.symbol) %>%
    normalize_gene() %>%
    unique()
  
  res_genes <- resistance %>%
    filter(Strain == s) %>%
    pull(GENE) %>%
    normalize_gene() %>%
    unique()
  
  shared      <- intersect(amr_genes, res_genes)
  only_amr    <- setdiff(amr_genes, res_genes)
  only_res    <- setdiff(res_genes, amr_genes)
  
  list(
    Strain             = s,
    shared             = shared,
    only_in_AMR_Finder = only_amr,
    only_in_resistance = only_res
  )
})

names(comparison_list) <- strains


# --- Summary table ---
summary_table <- bind_rows(lapply(comparison_list, function(x) {
  data.frame(
    Strain                 = x$Strain,
    N_shared               = length(x$shared),
    N_only_in_AMR_Finder   = length(x$only_in_AMR_Finder),
    N_only_in_resistance   = length(x$only_in_resistance),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    N_total        = N_shared + N_only_in_AMR_Finder + N_only_in_resistance,
    Pct_shared     = round(N_shared / N_total * 100, 1),
    Pct_only_AMR   = round(N_only_in_AMR_Finder / N_total * 100, 1),
    Pct_only_res   = round(N_only_in_resistance / N_total * 100, 1)
  )

print(summary_table)

totals_row <- summary_table %>%
  summarise(
    Strain               = "TOTAL",
    N_shared             = sum(N_shared),
    N_only_in_AMR_Finder = sum(N_only_in_AMR_Finder),
    N_only_in_resistance = sum(N_only_in_resistance),
    N_total              = sum(N_total)
  ) %>%
  mutate(
    Pct_shared   = round(N_shared / N_total * 100, 1),
    Pct_only_AMR = round(N_only_in_AMR_Finder / N_total * 100, 1),
    Pct_only_res = round(N_only_in_resistance / N_total * 100, 1)
  )



detail_table <- bind_rows(lapply(comparison_list, function(x) {
  bind_rows(
    if (length(x$shared) > 0)
      data.frame(Strain = x$Strain, GENE = x$shared,
                 Status = "Shared", stringsAsFactors = FALSE),
    if (length(x$only_in_AMR_Finder) > 0)
      data.frame(Strain = x$Strain, GENE = x$only_in_AMR_Finder,
                 Status = "Only in AMR_Finder", stringsAsFactors = FALSE),
    if (length(x$only_in_resistance) > 0)
      data.frame(Strain = x$Strain, GENE = x$only_in_resistance,
                 Status = "Only in resistance", stringsAsFactors = FALSE)
  )
}))

print(detail_table)
