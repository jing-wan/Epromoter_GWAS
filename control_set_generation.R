#!/usr/bin/Rscript

require(data.table)
require(stringr)
require(pheatmap)
require(reshape2)
require(ggpot2)

#============================control set generation===================================

#gene expression data download from Uhlen et al. 2016 (Table EV1)
gene_exp_62 <- fread("Table EV1 HPA and GTEx 62 tissues gene expression.txt", header = T)

gene_exp_HPA32 <- gene_exp_62[, c(1:4,7:38)]
gene_exp_GTEx30 <- gene_exp_62[, c(1,2,5,6,39:68)]

colnames(gene_exp_GTEx30) <- str_split(colnames(gene_exp_GTEx30), "\\.", 3, simplify = T)[, 3]

# remove those gene not in Ensembl_hg38_all_genes_TSS.bed (Epromoter annotation)
ensembl_hg38_gene <- fread("Ensembl_hg38_coding_transcript_TSS500_uniqTSS_rmTranscript.bed", header = F)
#some genes share same promoter also Epromoter 
ensembl_hg38_geneID <- unlist(strsplit(ensembl_hg38_gene$V5, ";"))

colnames(gene_exp_GTEx30)[1] <- "geneID"
gene_exp_GTEx30 <- gene_exp_GTEx30[gene_exp_GTEx30$geneID %in% ensembl_hg38_geneID, ]

#only keep fpkm value column
gene_exp_GTEx30_ma <- gene_exp_GTEx30[, c(5:34)]

gene_exp_GTEx30_ma <- log2(gene_exp_GTEx30_ma+1)


#=========== hierarchy cluster=================

gene_exp_GTEx30_ma_hcluster <- hclust(dist(gene_exp_GTEx30_ma, method = "euclidean"), method = "ward.D2")

#order by hcluster results
gene_exp_GTEx30_hcluster_order <- cbind(gene_exp_GTEx30[gene_exp_GTEx30_ma_hcluster$order, c(1:4)], 
                                        gene_exp_GTEx30_ma[gene_exp_GTEx30_ma_hcluster$order, ])

#calculate the nearest gene with Epromoter in cluster

all_Epromoters_list <- fread("Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15.bed", header = F)

all_Epromoters_gene <- all_Epromoters_list$V5

gene_alt_Ep_num <- as.data.frame(table(all_Epromoters_gene), stringsAsFactors = F)
colnames(gene_alt_Ep_num) <- c("geneID", "alt_Epromoter")


gene_exp_GTEx30_hcluster_order_Ep <- merge(gene_exp_GTEx30_hcluster_order, 
                                           gene_alt_Ep_num, 
                                           all.x = T, by.x = "geneID", by.y = "geneID", sort = F)


gene_exp_GTEx30_hcluster_order_Ep$type <- "NA"
gene_exp_GTEx30_hcluster_order_Ep[gene_exp_GTEx30_hcluster_order_Ep$alt_Epromoter >= 1, "type"] <- "Epromoter"
gene_exp_GTEx30_hcluster_order_Ep[gene_exp_GTEx30_hcluster_order_Ep$alt_Epromoter %in% NA, "type"] <- "other_promoter"


#--------find nearest expression genes as control -------------

Ep_index <- which(gene_exp_GTEx30_hcluster_order_Ep$type %in% "Epromoter")
Op_index <- which(gene_exp_GTEx30_hcluster_order_Ep$type %in% "other_promoter")

Op_nearest_index_all <- integer()
for (i in 1:length(Ep_index)) {
  #find nearest index by which.min(abs())
  Op_nearest_index <- Op_index[which.min(abs(Op_index-Ep_index[i]))]
  #store nearest index in vector
  Op_nearest_index_all <- c(Op_nearest_index_all, Op_nearest_index)
  #remove previous nearest index to avoid find duplicate value
  Op_index <- Op_index[-which.min(abs(Op_index-Ep_index[i]))]
}

gene_exp_GTEx30_hcluster_order_Ep[Op_nearest_index_all, "type"] <- "control_promoter"


write.table(gene_exp_GTEx30_hcluster_order_Ep, 
            file = "gene_exp_GTEx30_hcluster_order_Ep.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")


#===========================Epromoter VS nonEpromoter gene expression=========


gene_exp_GTEx30_hcluster_order_Ep <- fread("gene_exp_GTEx30_hcluster_order_Ep.txt", header = T)

gene_exp_Ep_ctl_op <- gene_exp_GTEx30_hcluster_order_Ep[, c(5:34, 36)]

gene_exp_Ep_ctl_op_long <- melt(gene_exp_Ep_ctl_op)

gene_exp_Ep_ctl_op_long[gene_exp_Ep_ctl_op_long$type == "control_promoter", "type"] <- "control"


#compare with total
gene_exp_total_long <- gene_exp_Ep_ctl_op_long
gene_exp_total_long$type <- "total_genes"

gene_exp_total_Ep_ctl_long <- rbind(gene_exp_Ep_ctl_op_long[gene_exp_Ep_ctl_op_long$type == "Epromoter", ],
                                    gene_exp_Ep_ctl_op_long[gene_exp_Ep_ctl_op_long$type == "control", ],
                                    gene_exp_total_long)


gene_exp_total_Ep_ctl_long$type <- factor(gene_exp_total_Ep_ctl_long$type, 
                                          levels = c("total_genes", "Epromoter", "control"))


pdf("gene_expression_30tissues_log2_total-Ep-cp_boxplot.pdf", width = 12, height = 6)


ggplot(gene_exp_total_Ep_ctl_long, aes(x = variable, y = value, fill = type))+
  geom_boxplot(lwd=0.1, outlier.shape = NA)+
  scale_fill_manual(values = c("#B2DF8A", "#1F78B4", "gray"))+
  ylim(0, 6)+
  labs(title = "gene expression in 30 tissues", 
       x = "", y = "log2(fpkm+1)",
       fill = "gene type")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=13, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = "top")+
  annotate("text", x = 1:30-0.1, y = 6, label = "***", size = 2)+
  annotate("text", x = 1:30+0.1, y = 5.5, label = "NS", size = 2)

dev.off()



#===============test significant by each group===================
require(dplyr)
options(scipen=2)

#keep two type, do group by tissue, do test by group, summarise pvalue
#all tissues p-values between Epromoters and control genes are >0.01, no significant difference
Ep_ctl_pvalue <-
  gene_exp_Ep_ctl_op_long[!(gene_exp_Ep_ctl_op_long$type %in% "other_promoter"), ] %>%
  group_by(variable) %>%
  do(w = wilcox.test(value~type, data=., paired=FALSE)) %>%
  summarise(variable, p_value = w$p.value)

#all tissues p-values between Epromoters and other genes are < 2.2e-16, significant difference
Ep_op_pvalue <-
  gene_exp_Ep_ctl_op_long[!(gene_exp_Ep_ctl_op_long$type %in% "control_promoter"), ] %>%
  group_by(variable) %>%
  do(w = wilcox.test(value~type, data=., paired=FALSE)) %>%
  summarise(variable, p_value = w$p.value)

Ep_op_pvalue$p_value <- 2.2e-16

Ep_ctl_op_pvalue <- data.frame(tissue = Ep_ctl_pvalue$variable,
                               Epromoter_control_pvalue = Ep_ctl_pvalue$p_value,
                               Epromoter_other_pvalue = Ep_op_pvalue$p_value,
                               stringsAsFactors = F)

write.table(Ep_ctl_op_pvalue, file = "Epromoter_control_other_gene_expression_compare_pvalue.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")



#=================get control gene pool TSS500===================

hg38_promoters <- fread("Ensembl_hg38_coding_transcript_TSS500_uniqTSS_rmTranscript.bed", header = F)
all_Epromoters_list <- fread("Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15.bed", header = F)

all_Epromoters_geneID <- unique(all_Epromoters_list$V5)

control_promoter_geneID <- gene_exp_GTEx30_hcluster_order_Ep[gene_exp_GTEx30_hcluster_order_Ep$type %in% "control_promoter", "ensgid"]

hg38_other_geneID <-
  hg38_promoters[!(hg38_promoters$V5 %in% all_Epromoters_geneID | 
                     hg38_promoters$V5 %in% control_promoter_geneID$ensgid), "V5"]

set.seed(2023)
control_geneID2 <- sample(unique(hg38_other_geneID$V5), 213)

control_promoter_gene_pool <-
  hg38_promoters[hg38_promoters$V5 %in% control_promoter_geneID$ensgid | 
                   hg38_promoters$V5 %in% control_geneID2, ]

write.table(control_promoter_gene_pool,
            file = "control_gene_pool_all_promoters_TSS500.bed",
            col.names = F, row.names = F, quote = F, sep = "\t")


#================take 5743 control promoters================

control_promoter_gene_pool_merged <- fread("control_gene_pool_all_promoters_TSS500_merged.bed", header = F)

#at least keep 1 TSS for each control gene
control_promoter_merged_set1 <- control_promoter_gene_pool_merged[!duplicated(control_promoter_gene_pool_merged$V5), ]

control_promoter_gene_pool_merged_rest <- control_promoter_gene_pool_merged[!(control_promoter_gene_pool_merged$V2 %in% control_promoter_merged_set1$V2), ]

control_promoter_merged_set2 <- control_promoter_gene_pool_merged_rest[sample(1:nrow(control_promoter_gene_pool_merged_rest), 191), ]

#total merged control promoters=5743
control_promoter_merged_total <- rbind(control_promoter_merged_set1, control_promoter_merged_set2)

control_promoter_merged_total <- control_promoter_merged_total[order(control_promoter_merged_total$V1, control_promoter_merged_total$V2), ]

write.table(control_promoter_merged_total, 
            file = "control_promoters.hg38_coding_gene.5743merged.bed", 
            col.names = F, row.names = F, quote = F, sep = "\t")



#===================examples to show the control of Epromoters=================

require(data.table)
require(tidyverse)

gene_exp_hcluster_order_Ep <- fread("gene_exp_GTEx30_hcluster_order_Ep.txt", header = T)

Ep_example1 <- gene_exp_hcluster_order_Ep[gene_exp_hcluster_order_Ep$gene_name %in% "RPS3", ]
cp_example1 <- gene_exp_hcluster_order_Ep[gene_exp_hcluster_order_Ep$gene_name %in% "RPS8", ]

Ep_example2 <- gene_exp_hcluster_order_Ep[gene_exp_hcluster_order_Ep$gene_name %in% "HBG2", ]
cp_example2 <- gene_exp_hcluster_order_Ep[gene_exp_hcluster_order_Ep$gene_name %in% "HBG1", ]


example <- rbind(Ep_example1, cp_example1, Ep_example2, cp_example2)

example[example$type %in% "control_promoter", "type"] <- "control"

example_long <-
  example[, c(2,5:34,36)] %>%
  pivot_longer(cols = !c(gene_name,type),
               names_to = "tissue", values_to = "expression")


#order rows by factor levels
example_long$gene_name <- factor(example_long$gene_name, levels = c("HBG1", "HBG2", "RPS8", "RPS3"))
example_long$type <- factor(example_long$type, levels = c("Epromoter", "control"))


pdf("Epromoter-control gene 2 examples expression heatmap.pdf", width = 9, height = 3)

ggplot(example_long)+
  geom_tile(aes(x = tissue, y = gene_name, alpha = expression, fill = type), color = "gray")+
  scale_fill_manual(values = c("#377EB8", "gray50"))+
  scale_alpha(range = c(0.1, 1))+
  labs(title = "2 examples: Epromoter-control genes")+
  theme_classic()+ 
  theme(#plot.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(angle = 45, hjust=1, vjust = 1, size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line= element_blank(),
        axis.ticks = element_blank(),
        # legend.position = "top",
        # legend.title = element_blank(),
        legend.text = element_text(size = 16))

dev.off()



