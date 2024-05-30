#!/usr/bin/Rscript

require(data.table)
require(ggplot2)
require(tidyverse)

#====================GWAS enrichment in Epromoter====================

#Epromoter SNPs trait frequency
Epromoter_SNP_trait_num <- fread("Epromoter_LD_SNP_EFO-trait_parent-trait.trait-frequency.txt", header = F)
colnames(Epromoter_SNP_trait_num) <- c("Freq", "GWAS_trait")
all_SNP_trait_num <- fread("../../3.1GWAS_variants_list/all_ld_SNP_GWAS_trait.trait-frequency.unique.txt", header = F)
colnames(all_SNP_trait_num) <- c("Freq", "GWAS_trait")

#all GWAS LD SNPs trait frequency
all_SNP_trait_num_unique <- 
  all_SNP_trait_num %>%
  group_by(GWAS_trait) %>%
  summarise(Freq_unique = sum(Freq))

#merge Epromoter SNPs and all GWAS LD SNPs trait frequency
Epromoter_SNP_all_SNP_trait_num <- 
  merge(Epromoter_SNP_trait_num, all_SNP_trait_num_unique, 
        by.x = "GWAS_trait", by.y = "GWAS_trait", all.x = T, all.y = F)
colnames(Epromoter_SNP_all_SNP_trait_num) <- c("GWAS_trait", "Epromoter_SNP", "all_GWAS_LD_SNP")

#ratio: the percent of Epromoter SNPs in all GWAS LD SNPs
Epromoter_SNP_all_SNP_trait_num$enriched_ratio <- 
  Epromoter_SNP_all_SNP_trait_num$Epromoter_SNP/Epromoter_SNP_all_SNP_trait_num$all_GWAS_LD_SNP

#ranking by ratio
Epromoter_SNP_all_SNP_trait_ratio <- Epromoter_SNP_all_SNP_trait_num[order(Epromoter_SNP_all_SNP_trait_num$enriched_ratio, decreasing = T), ]




#====================GWAS enrichment in control====================

#Epromoter SNPs trait frequency
Epromoter_SNP_trait_num <- fread("control_promoter_LD_SNP_EFO-trait_parent-trait.trait-frequency.txt", header = F)
colnames(Epromoter_SNP_trait_num) <- c("Freq", "GWAS_trait")
all_SNP_trait_num <- fread("../../3.1GWAS_variants_list/all_ld_SNP_GWAS_trait.trait-frequency.unique.txt", header = F)
colnames(all_SNP_trait_num) <- c("Freq", "GWAS_trait")

#all GWAS LD SNPs trait frequency
all_SNP_trait_num_unique <- 
  all_SNP_trait_num %>%
  group_by(GWAS_trait) %>%
  summarise(Freq_unique = sum(Freq))

#merge Epromoter SNPs and all GWAS LD SNPs trait frequency
Epromoter_SNP_all_SNP_trait_num <- 
  merge(Epromoter_SNP_trait_num, all_SNP_trait_num_unique, 
        by.x = "GWAS_trait", by.y = "GWAS_trait", all.x = T, all.y = F)
colnames(Epromoter_SNP_all_SNP_trait_num) <- c("GWAS_trait", "Epromoter_SNP", "all_GWAS_LD_SNP")

#ratio: the percent of Epromoter SNPs in all GWAS LD SNPs
Epromoter_SNP_all_SNP_trait_num$enriched_ratio <- 
  Epromoter_SNP_all_SNP_trait_num$Epromoter_SNP/Epromoter_SNP_all_SNP_trait_num$all_GWAS_LD_SNP

#ranking by ratio
control_SNP_all_SNP_trait_ratio <- Epromoter_SNP_all_SNP_trait_num[order(Epromoter_SNP_all_SNP_trait_num$enriched_ratio, decreasing = T), ]





#=====================GWAS enrichment-difference-hypergeometric test========================================


GWAS_Epromoter_enrich <- Epromoter_SNP_all_SNP_trait_ratio
GWAS_control_enrich <- control_SNP_all_SNP_trait_ratio

GWAS_Epromoter_VS_control_enrich <- read.table("Epromoter_VS_control_GWAS_enrichment_bySNPs.txt",
                                               header = T, sep = "\t", quote = "")

#=====calculate hypergeometric_test_pvalue =======

#Epromoter
for (i in 1:nrow(GWAS_Epromoter_enrich)) {
  GWAS_Epromoter_enrich[i, "hypergeometric_test_pvalue"] <- 
    phyper(as.numeric(GWAS_Epromoter_enrich[i, "Epromoter_SNP"]), 
           as.numeric(GWAS_Epromoter_enrich[i, "all_GWAS_LD_SNP"]), 
           2445359-as.numeric(GWAS_Epromoter_enrich[i, "all_GWAS_LD_SNP"]),
           4330, lower.tail = F)
}

#control
for (i in 1:nrow(GWAS_control_enrich)) {
  GWAS_control_enrich[i, "hypergeometric_test_pvalue"] <- 
    phyper(as.numeric(GWAS_control_enrich[i, "Epromoter_SNP"]), 
           as.numeric(GWAS_control_enrich[i, "all_GWAS_LD_SNP"]), 
           2445359-as.numeric(GWAS_control_enrich[i, "all_GWAS_LD_SNP"]),
           4062, lower.tail = F)
}

colnames(GWAS_control_enrich)[2] <- c("control_SNP")
colnames(GWAS_Epromoter_VS_control_enrich)[1:2] <- c("GWAS_trait", "category")


GWAS_Epromoter_VS_control_enrich1 <- 
  merge(GWAS_Epromoter_VS_control_enrich, GWAS_Epromoter_enrich, 
        by = "GWAS_trait", all = T)

GWAS_Epromoter_VS_control_enrich2 <- 
  merge(GWAS_Epromoter_VS_control_enrich1, GWAS_control_enrich, 
        by = "GWAS_trait", all = T)

GWAS_Epromoter_VS_control_enrich3 <-
  GWAS_Epromoter_VS_control_enrich2[, c(1,2,9:16)]

#calculate ratio Epromoter_VS_control
GWAS_Epromoter_VS_control_enrich3$ratio.Ep_ctl <-
  GWAS_Epromoter_VS_control_enrich3$Epromoter_SNP/GWAS_Epromoter_VS_control_enrich3$control_SNP

#calculate pvalue.chisq.test Epromoter_VS_control
GWAS_Epromoter_VS_control_enrich3$pvalue.Ep_ctl <- NA

for (i in 1:nrow(GWAS_Epromoter_VS_control_enrich3)) {
  if (is.na(GWAS_Epromoter_VS_control_enrich3$Epromoter_SNP[i])==F & is.na(GWAS_Epromoter_VS_control_enrich3$control_SNP[i])==F) {
    pvalue.chisq.test<- 
      chisq.test(matrix(c(GWAS_Epromoter_VS_control_enrich3$Epromoter_SNP[i], 4330, GWAS_Epromoter_VS_control_enrich3$control_SNP[i], 4062), nrow = 2))
    GWAS_Epromoter_VS_control_enrich3$pvalue.Ep_ctl[i] <- pvalue.chisq.test$p.value
  }
}

write.table(GWAS_Epromoter_VS_control_enrich3,
            file = "GWAS_trait_enrichment_hypergeometric_in_Epromoter_VS_control.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")



#filter by difference 0.05 (GWAS_Epromoter_VS_control_enrich2$pvalue.x < 0.05)

GWAS_Epromoter_VS_control_enrich2_diff05 <-
  GWAS_Epromoter_VS_control_enrich2[GWAS_Epromoter_VS_control_enrich2$pvalue < 0.05, ]

GWAS_Epromoter_VS_control_enrich2_diff05 <-
  GWAS_Epromoter_VS_control_enrich2_diff05[!is.na(GWAS_Epromoter_VS_control_enrich2_diff05$GWAS_trait), ]

#filter by hypergeometric_test_pvalue < 0.001
GWAS_Epromoter_VS_control_enrich2_diff05_hyper001 <-
  GWAS_Epromoter_VS_control_enrich2_diff05[GWAS_Epromoter_VS_control_enrich2_diff05$hypergeometric_test_pvalue.x < 0.001 |
                                             GWAS_Epromoter_VS_control_enrich2_diff05$hypergeometric_test_pvalue.y < 0.001, ]



#filter by category
Ep_trait_plot <-
  GWAS_Epromoter_VS_control_enrich2_diff05_hyper001[!is.na(GWAS_Epromoter_VS_control_enrich2_diff05_hyper001$category) &
                                                      !(GWAS_Epromoter_VS_control_enrich2_diff05_hyper001$category %in% c("Body_measurement", "Other_measurement", "Other_trait", "Response_to_drug")), ]

Ep_trait_plot$GWAS_trait <- gsub("_", " ", Ep_trait_plot$GWAS_trait)
Ep_trait_plot$GWAS_trait <- factor(Ep_trait_plot$GWAS_trait, levels = Ep_trait_plot$GWAS_trait)
Ep_trait_plot$category <- factor(Ep_trait_plot$category)

Ep_trait_plot <- Ep_trait_plot[order(Ep_trait_plot$ratio, decreasing = T), ]

Ep_trait_plot$GWAS_trait <- as.character(Ep_trait_plot$GWAS_trait)
Ep_trait_plot$GWAS_trait <- factor(Ep_trait_plot$GWAS_trait, levels = Ep_trait_plot$GWAS_trait)


Ep_trait_plot_df <- Ep_trait_plot[, c("GWAS_trait", "category", "hypergeometric_test_pvalue.x", "hypergeometric_test_pvalue.y")]
colnames(Ep_trait_plot_df) <- c("GWAS_trait", "category", "Epromoter", "control")

Ep_trait_plot_df <- Ep_trait_plot_df %>%
  pivot_longer(
    cols = c(Epromoter, control),
    names_to = "Epromoter_control",
    values_to = "hypergeometric_test_pvalue"
  )

Ep_trait_plot_df$Epromoter_control <- factor(Ep_trait_plot_df$Epromoter_control, levels = c("Epromoter", "control"))

# Ep_trait_plot_df$pattern_fill <- ifelse(Ep_trait_plot_df$Epromoter_control == "Epromoter", "stripe", "none")


pdf("GWAS_trait_enrichment_hypergeometric_in_Epromoter-control_differential0.05.pdf", width = 10, height = 5)

ggplot(data = Ep_trait_plot_df, 
       mapping = aes(x = GWAS_trait, y = -log10(hypergeometric_test_pvalue), 
                     fill = Epromoter_control,
                     # fill = category, size = Epromoter_control, color = Epromoter_control, #pattern = pattern_fill, #
                     group = interaction(GWAS_trait, Epromoter_control)))+
  geom_bar(stat = "identity", width=0.8, position = "dodge")+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  scale_y_continuous(expand = c(0, 0))+
  labs(title = "GWAS trait enrichment in Epromoter/control (differential pvalue<0.05)", 
       # x = "", 
       y = "Enrichment   -log10(p-value)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(hjust = 1, angle = 45, size=12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        plot.margin = unit(c(1,1,0,4),"cm"),
        legend.position = c(0.8,0.6),
        legend.text = element_text(size=12),
        # legend.title = element_text(size=14),
        legend.title = element_blank())

dev.off()






