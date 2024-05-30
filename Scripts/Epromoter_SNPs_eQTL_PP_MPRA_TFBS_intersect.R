#!/usr/bin/Rscript

#====================Epromoter_distal_eQTL_PP_MPRA_TFBS_intersect===============


require(ggvenn)

Ep_GWAS <- fread("Epromoter_variants_annotation_SuppTable5.txt", header = T)

Ep_GWAS <- Ep_GWAS[, !c("TFs_affected_FABIAN")]

#count GWAS trait per SNP
for (i in 1:nrow(Ep_GWAS)) {
  Ep_GWAS$GWAS_trait_num[i] <- length(unique(unlist(strsplit(Ep_GWAS$GWAS_trait[i], ";|,"))))
}

eQTL_distal <- Ep_GWAS[Ep_GWAS$eQTL_effect %in% c("distal;proximal", "distal"), "rsID"]
PP_interaction <- Ep_GWAS[!is.na(Ep_GWAS$PP_target), "rsID"]
MPRA <- Ep_GWAS[!is.na(Ep_GWAS$MPRA_allelic_source), "rsID"]
TF_effect <- Ep_GWAS[!(is.na(Ep_GWAS$TFs_affected_ANANASTRA) & is.na(Ep_GWAS$TFs_affected_SNP2TFBS)), "rsID"]


#intersection
Epromoter_eQTL_PP_MPRA_TFBS <- list(eQTL_distal = unique(eQTL_distal$rsID), 
                                    PP_interaction = unique(PP_interaction$rsID),
                                    MPRA = unique(MPRA$rsID),
                                    TF_effect = unique(TF_effect$rsID))


pdf("Epromoter_SNPs_eQTL-PP-MPRA_TF_overlapping.pdf", width = 6, height = 5)

venn_ov_gg <- 
  ggvenn(Epromoter_eQTL_PP_MPRA_TFBS, 
         c("eQTL_distal", "PP_interaction", "MPRA", "TF_effect"), 
         show_percentage = F, fill_alpha = 0.8, stroke_size = 0.8, text_size = 6)+
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))

dev.off()

#get overlap list from venn diagram
intersect_SNPs <- Reduce(intersect, Epromoter_eQTL_PP_MPRA_TFBS)

Epromoter_eQTL_PP_MPRA_TFBS_intersect <- Ep_GWAS[Ep_GWAS$rsID %in% intersect_SNPs, ]

write.table(Epromoter_eQTL_PP_MPRA_TFBS_intersect, file = "Epromoter_distal_eQTL_PP_MPRA_TFBS_intersect_list.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")





