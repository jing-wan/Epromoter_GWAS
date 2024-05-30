#!/usr/bin/Rscript

require(data.table)
require(ggplot2)
require(ggsignif)

#========================all gene VS Epromoter VS control average gene expression================

gene_exp_order_Ep_cp_op <- fread("gene_exp_GTEx30_hcluster_order_Ep.txt", header = T)

gene_exp_order_Ep_cp_op$average_exp <- rowMeans(gene_exp_order_Ep_cp_op[, 5:34])

gene_exp_order_all <- gene_exp_order_Ep_cp_op
gene_exp_order_all$type <- "total_genes"

gene_exp_order_Ep_cp <- gene_exp_order_Ep_cp_op[gene_exp_order_Ep_cp_op$type %in% c("Epromoter", "control_promoter"), ]

gene_exp_order_Ep_cp[gene_exp_order_Ep_cp$type %in% "control_promoter", "type"] <- "control"

gene_exp_all_Ep_cp <- rbind(gene_exp_order_all, gene_exp_order_Ep_cp)

gene_exp_all_Ep_cp$type <- factor(gene_exp_all_Ep_cp$type, levels = c("total_genes", "Epromoter", "control"))

median_value_df <- gene_exp_all_Ep_cp[gene_exp_all_Ep_cp$type %in% "Epromoter", "average_exp"]


pdf("gene_expression_average_total_Epromoter-control_violinboxplot.pdf", width = 5, height = 4)

ggplot(gene_exp_all_Ep_cp, aes(x = type, y = average_exp, fill = type))+
  geom_violin(width = 0.8)+
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  # scale_fill_manual(values = c("#B2DF8A", "#1F78B4", "gray"))+#4DAF4A
  scale_fill_manual(values = c("#4DAF4A", "#1F78B4", "gray"))+
  geom_hline(yintercept = median(median_value_df$average_exp), linetype = 2, color = "gray0")+
  scale_x_discrete(labels = c("total", "Epromoter", "control"))+
  labs(title = "gene expression in average of 30 tissues", 
       x = "", y = "log2(fpkm+1)",
       fill = "gene type")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("total_genes", "Epromoter"), 
                                 c("Epromoter", "control")), 
              y_position = 9, textsize = 5, linetype = "blank")


dev.off()




#============Tissue specificity ======================

gene_exp_order_Ep_cp_op_RowMaxNorm <- gene_exp_all_Ep_cp
#get max expression between 30 tisssues for each gene
gene_exp_order_Ep_cp_op_RowMaxNorm$max_exp <- apply(gene_exp_order_Ep_cp_op_RowMaxNorm[, 5:34], 1, max)

#expression normalized between 30 tissues
gene_exp_order_Ep_cp_op_RowMaxNorm[, 5:34] <- gene_exp_order_Ep_cp_op_RowMaxNorm[, 5:34]/gene_exp_order_Ep_cp_op_RowMaxNorm$max_exp

#tissue_specificity_score
gene_exp_order_Ep_cp_op_RowMaxNorm$tissue_specificity_score <- rowSums(1-gene_exp_order_Ep_cp_op_RowMaxNorm[, 5:34])/(30-1)

#remove NA
gene_exp_order_Ep_cp_op_RowMaxNorm <- gene_exp_order_Ep_cp_op_RowMaxNorm[!is.na(gene_exp_order_Ep_cp_op_RowMaxNorm$tissue_specificity_score), ]

#median value
tissue_specificity_Ep_median_df <- gene_exp_order_Ep_cp_op_RowMaxNorm[gene_exp_order_Ep_cp_op_RowMaxNorm$type %in% "control", "tissue_specificity_score"]


pdf("tissue_specificity_score_total-Epromoter-control_violinboxplot.pdf", width = 5, height = 4)

ggplot(gene_exp_order_Ep_cp_op_RowMaxNorm, aes(x = type, y = tissue_specificity_score, fill = type))+
  geom_violin(width = 0.9)+
  geom_boxplot(width = 0.15)+
  scale_fill_manual(values = c("#4DAF4A", "#1F78B4", "gray"))+
  geom_hline(yintercept = median(tissue_specificity_Ep_median_df$tissue_specificity_score), linetype = 2, color = "gray0")+
  scale_y_continuous(limits = c(0,1.05), breaks = c(0.1, 0.25, 0.5, 0.75, 1))+
  scale_x_discrete(labels = c("total", "Epromoter", "control"))+
  labs(title = "tissue specificity in 30 tissues", 
       x = "", y = "tissue specificity score",
       fill = "gene type")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("total_genes", "Epromoter"), 
                                 c("Epromoter", "control")), 
              y_position = 0.95, textsize = 5, linetype = "blank")


dev.off()











