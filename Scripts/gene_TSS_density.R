#!/usr/bin/Rscript

hg38_promoters <- fread("../Ensembl_hg38_coding_transcript_TSS500_uniqTSS_rmTranscript.bed", header = F)

gene_exp_GTEx30_hcluster_order_Ep <- fread("../gene_exp_GTEx30_hcluster_order_Ep.txt", header = T)

all_Epromoters_geneID <- gene_exp_GTEx30_hcluster_order_Ep[gene_exp_GTEx30_hcluster_order_Ep$type %in% "Epromoter", "ensgid"]

control_promoter_geneID <- gene_exp_GTEx30_hcluster_order_Ep[gene_exp_GTEx30_hcluster_order_Ep$type %in% "control_promoter", "ensgid"]


Epromoter_gene_allTSS <- hg38_promoters[hg38_promoters$V5 %in% all_Epromoters_geneID$ensgid, ]

control_promoter_gene_allTSS <- hg38_promoters[hg38_promoters$V5 %in% control_promoter_geneID$ensgid, ]

Epromoter_gene_TSS_freq <- as.data.frame(table(Epromoter_gene_allTSS$V6))
Epromoter_gene_TSS_freq$category <- "Epromoter"
control_promoter_gene_TSS_freq <- as.data.frame(table(control_promoter_gene_allTSS$V6))
control_promoter_gene_TSS_freq$category <- "control"
hg38_promoter_gene_TSS_freq <- as.data.frame(table(hg38_promoters$V6))
hg38_promoter_gene_TSS_freq$category <- "total_genes"

gene_TSS_freq <- rbind(Epromoter_gene_TSS_freq[, 2:3], control_promoter_gene_TSS_freq[, 2:3], hg38_promoter_gene_TSS_freq[, 2:3])
colnames(gene_TSS_freq) <- c("TSS_number", "category")

gene_TSS_freq$category <- factor(gene_TSS_freq$category, levels = c("total_genes", "Epromoter", "control"))
gene_TSS_freq_Epromoter_median <- median(gene_TSS_freq[gene_TSS_freq$category %in% "Epromoter", ][, 1])
gene_TSS_freq_total_median <- median(gene_TSS_freq[gene_TSS_freq$category %in% "total_genes", ][, 1])


pdf("gene_TSS_number_total-Epromoter-control_violinboxplot.pdf", width = 5, height = 4)

ggplot(gene_TSS_freq[gene_TSS_freq$TSS_number<=20, ], aes(x = category, y = TSS_number, fill = category))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1)+
  scale_fill_manual(values = c("#B2DF8A", "#1F78B4", "gray"))+
  geom_hline(yintercept = gene_TSS_freq_total_median, linetype = 2, color = "gray0")+
  scale_y_continuous(limits = c(1,21), breaks = c(1,3,5,10,15,20))+
  # scale_x_discrete(labels = c("total genes\n19612", "Epromoter\n5333", "control\n5333"))+
  labs(title = "gene TSS number", 
       x = "", y = "TSS number per gene",
       fill = "gene type")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("total_genes", "Epromoter"),
                                 c("Epromoter", "control"), 
                                 c("control", "total_genes")), 
              y_position = c(19, 19, 20))


dev.off()


