#!/usr/bin/Rscript


#========================merge all MPRA SNPs into one list============================

require(data.table)
require(dplyr)

snp_list_files <- list.files(path = "MPRA-SNP/SNP_list_files", pattern = ".txt")

#add file name as column resource
snp_list_all_df <- data.frame()

for (i in snp_list_files) {
  snp_list <- read.table(i, header = F)
  snp_list <- unique(snp_list)
  filename_base <- unlist(strsplit(i, ".txt"))
  # reference <- unlist(strsplit(filename_base, "_"))[1]
  # cell_line <- unlist(strsplit(filename_base, "_"))[2]
  snp_list_df <- data.frame(rsID = snp_list, resource = filename_base, stringsAsFactors = F)
  snp_list_all_df <- rbind(snp_list_all_df, snp_list_df)
}

colnames(snp_list_all_df) <- c("rsID", "resource")

#merge SNP from different resource
snp_list_all_df_merged <- 
  snp_list_all_df %>% 
  group_by(rsID) %>% 
  summarise(across(everything(), ~toString(.)))

write.table(snp_list_all_df_merged, file = "../published_MPRA_significant_allelic_SNP_list-2023-08.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")


#=================================================


mpra_snp <- fread("published_MPRA_significant_allelic_SNP_list-2023-08.txt", header = T)

Epromoter_snp <- fread("Epromoters_ovTSS_ov_GWAS-LD-SNP_2023-02-15_gwas_info_merged.bed", header = F)

#merge mpra snp list with Epromoter SNP

Epromoter_snp_mpra <- 
  merge(Epromoter_snp, mpra_snp, by.x = "V4", by.y = "rsID", all.x = T, all.y = F)

Epromoter_snp_mpra_rmNA <- Epromoter_snp_mpra[!is.na(Epromoter_snp_mpra$resource), ]

write.table(Epromoter_snp_mpra_rmNA, 
            file = "Epromoter_SNP_MPRA-validated_2023-08.txt", 
            col.names = F, row.names = F, quote = F, sep = "\t")

#========control========

control_snp <- fread("control_promoters_5743merged_GWAS-LD-SNP_gwas_info_merged.bed", header = F)

control_snp_mpra <- 
  merge(control_snp, mpra_snp, by.x = "V4", by.y = "rsID", all.x = T, all.y = F)


#===============GWAS trait per SNP: Epromoter/control MPRA/noMPRA=============

snp_GWAS_mpra <- Epromoter_snp_mpra[!is.na(Epromoter_snp_mpra$resource), c("V4", "V17")]
snp_GWAS_nompra <- Epromoter_snp_mpra[is.na(Epromoter_snp_mpra$resource), c("V4", "V17")]
control_snp_GWAS_mpra <- control_snp_mpra[!is.na(control_snp_mpra$resource), c("V4", "V14")]
control_snp_GWAS_nompra <- control_snp_mpra[is.na(control_snp_mpra$resource), c("V4", "V14")]

#count GWAS trait per SNP
for (i in 1:nrow(snp_GWAS_mpra)) {
  snp_GWAS_mpra$GWAS_trait_num[i] <- length(unique(unlist(strsplit(snp_GWAS_mpra$V17[i], ";|,"))))
}

for (i in 1:nrow(snp_GWAS_nompra)) {
  snp_GWAS_nompra$GWAS_trait_num[i] <- length(unique(unlist(strsplit(snp_GWAS_nompra$V17[i], ";|,"))))
}

for (i in 1:nrow(control_snp_GWAS_mpra)) {
  control_snp_GWAS_mpra$GWAS_trait_num[i] <- length(unique(unlist(strsplit(control_snp_GWAS_mpra$V14[i], ";|,"))))
}

for (i in 1:nrow(control_snp_GWAS_nompra)) {
  control_snp_GWAS_nompra$GWAS_trait_num[i] <- length(unique(unlist(strsplit(control_snp_GWAS_nompra$V14[i], ";|,"))))
}


snp_GWAS_mpra$type <- "Ep_SNP_MPRA"
snp_GWAS_nompra$type <- "Ep_SNP_noMPRA"
control_snp_GWAS_mpra$type <- "ctl_SNP_MPRA"
control_snp_GWAS_nompra$type <- "ctl_SNP_noMPRA"


Ep_ctl_snp_GWAS_mpra <- 
  rbind(snp_GWAS_mpra[, c("GWAS_trait_num", "type")], 
        snp_GWAS_nompra[, c("GWAS_trait_num", "type")],
        control_snp_GWAS_mpra[, c("GWAS_trait_num", "type")], 
        control_snp_GWAS_nompra[, c("GWAS_trait_num", "type")])

Ep_ctl_snp_GWAS_mpra$type <- factor(Ep_ctl_snp_GWAS_mpra$type, levels = c("Ep_SNP_MPRA", "Ep_SNP_noMPRA", "ctl_SNP_MPRA", "ctl_SNP_noMPRA"))

Ep_ctl_snp_GWAS_mpra[Ep_ctl_snp_GWAS_mpra$GWAS_trait_num > 10, "GWAS_trait_num"] <- 10


pdf("MPRA_allelic_SNPs_GWAS_trait.pdf", width = 5, height = 4)

ggplot(Ep_ctl_snp_GWAS_mpra, aes(x = type, y = GWAS_trait_num, fill = type))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(limits = c(1,11), breaks = c(1,2,5,10))+
  scale_x_discrete(labels = c("significant\nin MPRA", "other\nSNPs", "significant\nin MPRA", "other\nSNPs"))+
  scale_fill_manual(values = c("#377EB8", "#377EB8", "gray", "gray"))+
  geom_hline(yintercept = 2, linetype = 2, color = "gray")+
  labs(title = "GWAS trait of allelic effect SNPs", 
       y = "GWAS trait per SNP",
       x = "")+
  # x = "            Epromoter                 control")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        # axis.title.x = element_text(size=16, hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Ep_SNP_MPRA", "Ep_SNP_noMPRA"), c("ctl_SNP_MPRA", "ctl_SNP_noMPRA")), 
              y_position = 9, textsize = 4, linetype = "blank")+
  geom_signif(comparisons = list(c("Ep_SNP_MPRA", "ctl_SNP_MPRA"), c("Ep_SNP_noMPRA", "ctl_SNP_noMPRA")), 
              y_position = c(9.8, 10.5), textsize = 4)

dev.off()



