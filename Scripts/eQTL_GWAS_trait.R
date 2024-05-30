#!/usr/bin/Rscript

require(data.table)
require(tidyverse)
require(ggplot2)
require(ggsignif)


Epromoter_SNP_eQTL <- fread("Epromoter_GWAS_LD_SNP_eQTL_distance_merged.bed", header = F)
Epromoter_SNP_trait_num <- fread("Epromoter.trait-per-SNP.txt", header = F)

Epromoter_SNP_eQTL_trait_num <- merge(Epromoter_SNP_eQTL, Epromoter_SNP_trait_num, by.x = "V4", by.y = "V2", all.x = T)


#==========control promoter ==============

control_promoter_SNP_eQTL <- fread("control_promoter_GWAS_LD_SNP_eQTL_distance_merged.bed", header = F)
control_promoter_SNP_trait_num <- fread("control_promoter.trait-per-SNP.txt", header = F)

control_promoter_SNP_eQTL_trait_num <- merge(control_promoter_SNP_eQTL, control_promoter_SNP_trait_num, by.x = "V4", by.y = "V2", all.x = T)


#===merge Ep ctl======

Epromoter_eQTL_GWAS_trait <- 
  data.frame(GWAS_trait = Epromoter_SNP_eQTL_trait_num$V1.y, 
             eQTL_type = Epromoter_SNP_eQTL_trait_num$V41, 
             group = "Epromoter")

control_eQTL_GWAS_trait <- 
  data.frame(GWAS_trait = control_promoter_SNP_eQTL_trait_num$V1.y, 
             eQTL_type = control_promoter_SNP_eQTL_trait_num$V35, 
             group = "control")

Ep_ctl_eQTL_GWAS_trait <- rbind(Epromoter_eQTL_GWAS_trait, control_eQTL_GWAS_trait)

#remove NA, set 10 as maximal traits
Ep_ctl_eQTL_GWAS_trait <- 
  Ep_ctl_eQTL_GWAS_trait[!is.na(Ep_ctl_eQTL_GWAS_trait$GWAS_trait), ]

Ep_ctl_eQTL_GWAS_trait[Ep_ctl_eQTL_GWAS_trait$GWAS_trait > 10, "GWAS_trait"] <- 10

median_line <- Ep_ctl_eQTL_GWAS_trait[Ep_ctl_eQTL_GWAS_trait$eQTL_type %in% "proximal", ]
median_line <- median(median_line$GWAS_trait)

Ep_ctl_eQTL_GWAS_trait$eQTL_group <- paste(Ep_ctl_eQTL_GWAS_trait$eQTL_type, Ep_ctl_eQTL_GWAS_trait$group, sep = ".")
Ep_ctl_eQTL_GWAS_trait$eQTL_group <- factor(Ep_ctl_eQTL_GWAS_trait$eQTL_group, 
                                            levels = c("proximal.Epromoter", "proximal.control",
                                                       "distal.Epromoter", "distal.control", 
                                                       "distal;proximal.Epromoter", "distal;proximal.control"))

# Ep_ctl_eQTL_GWAS_trait$group <- factor(Ep_ctl_eQTL_GWAS_trait$group, levels = c("Epromoter", "control"))
# Ep_ctl_eQTL_GWAS_trait$eQTL_type <- factor(Ep_ctl_eQTL_GWAS_trait$eQTL_type, 
#                                            levels = c("proximal", "distal", "distal;proximal"))

pdf("Epromoter-control GWAS trait per eQTL by effect.pdf", width = 6, height = 4)

ggplot(Ep_ctl_eQTL_GWAS_trait, aes(x = eQTL_group, y = GWAS_trait, fill = eQTL_group))+
  # ggplot(Ep_ctl_eQTL_GWAS_trait, aes(x = eQTL_type, y = GWAS_trait, fill = group))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  # geom_violin(width = 1, position=position_dodge(0.9))+
  # geom_boxplot(width = 0.1, outlier.shape = NA, position=position_dodge(0.9))+
  scale_y_continuous(limits = c(1,10.5), breaks = c(1,2,3,5,10))+
  scale_x_discrete(labels = c("             proximal\n            Ep/ctl", "", "               distal\n                Ep/ctl", "", "               distal & proximal\n              Ep/ctl", ""))+
  # scale_fill_manual(values = c("#4DAF4A", "gray", "#377EB8", "gray", "#FFFF33", "gray"))+
  scale_fill_manual(values = c("#377EB8", "gray", "#377EB8", "gray", "#377EB8", "gray"))+
  geom_hline(yintercept = median_line, linetype = 2, color = "gray")+
  labs(title = "", y = "GWAS trait per eQTL")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("proximal.Epromoter", "proximal.control"),
                                 c("distal.Epromoter", "distal.control"), 
                                 c("distal;proximal.Epromoter", "distal;proximal.control")), 
              y_position = 9, textsize = 4, linetype = "blank")

dev.off()








