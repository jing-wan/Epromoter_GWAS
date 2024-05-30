

SNP_SELEX <- read.table("SNP-SELEX.txt", header = T)

Ep_GWAS <- read.table("Epromoters_ovTSS_ov_GWAS-LD-SNP_2023-02-15_gwas_info_merged.bed", header = F)
cp_GWAS <- read.table("control_promoters_5743merged_GWAS-LD-SNP_gwas_info_merged.bed", header = F)


Ep_GWAS_SELEX <- 
  merge(Ep_GWAS, SNP_SELEX, by.x = "V4", by.y = "rsID")

cp_GWAS_SELEX <- 
  merge(cp_GWAS, SNP_SELEX, by.x = "V4", by.y = "rsID")

Ep_GWAS_noSELEX <- Ep_GWAS[!(Ep_GWAS$V4 %in% Ep_GWAS_SELEX$V4), c("V4", "V17")]
cp_GWAS_noSELEX <- cp_GWAS[!(cp_GWAS$V4 %in% cp_GWAS_SELEX$V4), c("V4", "V14")]


#count GWAS trait per SNP
for (i in 1:nrow(Ep_GWAS_SELEX)) {
  Ep_GWAS_SELEX$GWAS_trait_num[i] <- length(unique(unlist(strsplit(Ep_GWAS_SELEX$V17[i], ";|,"))))
}

for (i in 1:nrow(cp_GWAS_SELEX)) {
  cp_GWAS_SELEX$GWAS_trait_num[i] <- length(unique(unlist(strsplit(cp_GWAS_SELEX$V14[i], ";|,"))))
}

for (i in 1:nrow(Ep_GWAS_noSELEX)) {
  Ep_GWAS_noSELEX$GWAS_trait_num[i] <- length(unique(unlist(strsplit(Ep_GWAS_noSELEX$V17[i], ";|,"))))
}

for (i in 1:nrow(cp_GWAS_noSELEX)) {
  cp_GWAS_noSELEX$GWAS_trait_num[i] <- length(unique(unlist(strsplit(cp_GWAS_noSELEX$V14[i], ";|,"))))
}



Ep_GWAS_SELEX$type <- "Ep_SNP_SELEX"
Ep_GWAS_noSELEX$type <- "Ep_SNP_noSELEX"
cp_GWAS_SELEX$type <- "ctl_SNP_SELEX"
cp_GWAS_noSELEX$type <- "ctl_SNP_noSELEX"


Ep_ctl_snp_GWAS_SELEX <- 
  rbind(Ep_GWAS_SELEX[, c("GWAS_trait_num", "type")], 
        Ep_GWAS_noSELEX[, c("GWAS_trait_num", "type")],
        cp_GWAS_SELEX[, c("GWAS_trait_num", "type")], 
        cp_GWAS_noSELEX[, c("GWAS_trait_num", "type")])

Ep_ctl_snp_GWAS_SELEX$type <- factor(Ep_ctl_snp_GWAS_SELEX$type, levels = c("Ep_SNP_SELEX", "Ep_SNP_noSELEX", "ctl_SNP_SELEX", "ctl_SNP_noSELEX"))

Ep_ctl_snp_GWAS_SELEX[Ep_ctl_snp_GWAS_SELEX$GWAS_trait_num > 10, "GWAS_trait_num"] <- 10


pdf("TF_binding_allelic_SNPs_GWAS_trait.v1.pdf", width = 5, height = 4)

ggplot(Ep_ctl_snp_GWAS_SELEX, aes(x = type, y = GWAS_trait_num, fill = type))+
  # geom_violin(width = 1, scale = "width")+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(limits = c(1,11), breaks = c(1,2,5,10))+
  scale_x_discrete(labels = c("significant\nSNP-SELEX", "other\nSNPs", "significant\nSNP-SELEX", "other\nSNPs"))+
  scale_fill_manual(values = c("#377EB8", "#377EB8", "gray", "gray"))+
  geom_hline(yintercept = 2, linetype = 2, color = "gray")+
  labs(title = "GWAS trait of TF binding allelic effect SNPs", 
       y = "GWAS trait per SNP",
       x = "          Epromoter                   control")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        # axis.title.x = element_text(size=16, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Ep_SNP_SELEX", "Ep_SNP_noSELEX"), c("ctl_SNP_SELEX", "ctl_SNP_noSELEX")), 
              y_position = 9, textsize = 4, linetype = "blank")+
  geom_signif(comparisons = list(c("Ep_SNP_SELEX", "ctl_SNP_SELEX"), c("Ep_SNP_noSELEX", "ctl_SNP_noSELEX")), 
              y_position = c(9.8, 10.5), textsize = 4)

dev.off()



