#!/usr/bin/Rscript

require(data.table)
require(tidyverse)


Ep_remap_TFs <- fread("eprom_nr_flag_bilan.tsv", header = F)
cp_remap_TFs <- fread("control_nr_flag_bilan.tsv", header = F)


Ep_TFs <- data.frame(promoterID = paste(Ep_remap_TFs$V1, Ep_remap_TFs$V2, Ep_remap_TFs$V3, sep = "."), 
                     gene = Ep_remap_TFs$V6, 
                     TFs = Ep_remap_TFs$V11,
                     type = "Epromoter")

cp_TFs <- data.frame(promoterID = paste(cp_remap_TFs$V1, cp_remap_TFs$V2, cp_remap_TFs$V3, sep = "."), 
                     gene = cp_remap_TFs$V6, 
                     TFs = cp_remap_TFs$V8,
                     type = "control")

#merge together
Ep_cp_TFs <- rbind(Ep_TFs, cp_TFs)
Ep_cp_TFs$promoterID_gene <- paste(Ep_cp_TFs$promoterID, Ep_cp_TFs$gene, sep = "_")

#split merged rows
Ep_cp_TFs_split <- 
  Ep_cp_TFs %>%
  separate_rows(TFs) %>%
  unique()

#TF exitst: 1
Ep_cp_TFs_split$True <- 1
Ep_cp_TFs_split <- Ep_cp_TFs_split[, c("promoterID_gene", "type", "TFs", "True")]

#convert into wide matrix
Ep_cp_TFs_ma <- 
  pivot_wider(Ep_cp_TFs_split, 
              names_from = TFs, values_from = True,
              values_fill = 0)

#=====add promoters without TFs binding======
Ep_list <- fread("../raw/Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15_merged.bed", header = F)
cp_list <- fread("../raw/control_promoters.hg38_coding_gene.5743merged.bed", header = F)
Ep_list$promoterID_gene <- paste(paste(Ep_list$V1, Ep_list$V2, Ep_list$V3, sep = "."), Ep_list$V6, sep = "_")
cp_list$promoterID_gene <- paste(paste(cp_list$V1, cp_list$V2, cp_list$V3, sep = "."), cp_list$V6, sep = "_")

Ep_TF0 <- data.frame(promoterID_gene = Ep_list$promoterID_gene[!(Ep_list$promoterID_gene %in% Ep_cp_TFs_ma$promoterID_gene)])
cp_TF0 <- data.frame(promoterID_gene = cp_list$promoterID_gene[!(cp_list$promoterID_gene %in% Ep_cp_TFs_ma$promoterID_gene)])

Ep_cp_TFs_ma1 <- merge(Ep_cp_TFs_ma, Ep_TF0, by = "promoterID_gene", all = T)
Ep_cp_TFs_ma1 <- merge(Ep_cp_TFs_ma1, cp_TF0, by = "promoterID_gene", all = T)

Ep_cp_TFs_ma1[Ep_cp_TFs_ma1$promoterID_gene %in% Ep_TF0$promoterID_gene, 2] <- "Epromoter"
Ep_cp_TFs_ma1[Ep_cp_TFs_ma1$promoterID_gene %in% cp_TF0$promoterID_gene, 2] <- "control"

Ep_cp_TFs_ma1[Ep_cp_TFs_ma1$promoterID_gene %in% Ep_TF0$promoterID_gene, 3:1198] <- 0
Ep_cp_TFs_ma1[Ep_cp_TFs_ma1$promoterID_gene %in% cp_TF0$promoterID_gene, 3:1198] <- 0


#========UMAP cluster===========

# Ep_TFs_umap_res <- umap(Ep_TFs_ma[, 2:1191])
# Ep_TFs_umap_res <- umap(Ep_TFs_ma[, 2:1194])

Ep_cp_TFs_ma1_res <- umap(Ep_cp_TFs_ma1[, 3:1198])

Ep_cp_TFs_umap_res_layout <- as.data.frame(Ep_cp_TFs_ma1_res$layout)

Ep_cp_TFs_umap_res_layout$promoterID_gene <- Ep_cp_TFs_ma1$promoterID_gene
Ep_cp_TFs_umap_res_layout$type <- Ep_cp_TFs_ma1$type

write.table(Ep_cp_TFs_umap_res_layout, file = "Ep_cp_TFs_umap_res.txt", 
            col.names = T, row.names = F, sep = "\t", quote = F)



#========UMAP scatterplot================

setwd("C:/work/sync/工作/PhD project/Epromoter_variants.v5.1/12.Genomic_features/Epromoter_ReMap/UMAP_cluster_byTFs")

Ep_cp_TFs_umap_res_layout <- read.table("Ep_cp_TFs_umap_res.txt", header = T)

Ep_cp_TFs_umap_res_layout$type <- factor(Ep_cp_TFs_umap_res_layout$type, levels = c("Epromoter", "control"))


# pdf("UMAP cluster promoters by TFs binding.pdf", width = 6, height = 6)
pdf("UMAP cluster promoters by TFs binding 3groups.pdf", width = 6, height = 6)

ggplot(Ep_cp_TFs_umap_res_layout, 
       aes(x = Ep_cp_TFs_umap_res_layout$V1, y = Ep_cp_TFs_umap_res_layout$V2, color = type))+
  geom_point(size = 0.3)+
  scale_x_continuous(breaks = c(-10, -7, 0, 10))+
  scale_color_manual(values = c("#377EB8", "gray"))+
  geom_vline(xintercept = c(-7, 0), linetype = 2, color = "gray")+
  labs(title = "UMAP cluster promoters by all TFs(ReMap) binding",
       x = "UMAP1", y = "UMAP2")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        # axis.ticks = element_blank(),
        axis.title = element_text(size=16),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  guides(color = guide_legend(override.aes = list(size=3)))

dev.off()



#==========add ReMap TFs binding count for each Epromoter/control================


Ep_remap_TFs <- fread("eprom_nr_flag_bilan.tsv", header = F)
cp_remap_TFs <- fread("control_nr_flag_bilan.tsv", header = F)

Ep_TFs_num <- data.frame(promoterID = paste(Ep_remap_TFs$V1, Ep_remap_TFs$V2, Ep_remap_TFs$V3, sep = "."), 
                         gene = Ep_remap_TFs$V6, 
                         TFs_num = Ep_remap_TFs$V10,
                         type = "Epromoter")

cp_TFs_num <- data.frame(promoterID = paste(cp_remap_TFs$V1, cp_remap_TFs$V2, cp_remap_TFs$V3, sep = "."), 
                         gene = cp_remap_TFs$V6, 
                         TFs_num = cp_remap_TFs$V7,
                         type = "control")

#merge together
Ep_cp_TFs_num <- rbind(Ep_TFs_num, cp_TFs_num)
Ep_cp_TFs_num$promoterID_gene <- paste(Ep_cp_TFs_num$promoterID, Ep_cp_TFs_num$gene, sep = "_")

#add ReMap TFs num
Ep_cp_TFs_umap_res_TFnum <- merge(Ep_cp_TFs_umap_res_layout, Ep_cp_TFs_num[, c("promoterID_gene", "TFs_num")], 
                                  by = "promoterID_gene", all.x = T)

Ep_cp_TFs_umap_res_TFnum[is.na(Ep_cp_TFs_umap_res_TFnum$TFs_num), "TFs_num"] <- 0

Ep_cp_TFs_umap_res_TFnum$type <- factor(Ep_cp_TFs_umap_res_TFnum$type, levels = c("Epromoter", "control"))

Ep_cp_TFs_umap_res_TFnum1 <-
  Ep_cp_TFs_umap_res_TFnum %>% filter(type == "control")
Ep_cp_TFs_umap_res_TFnum2 <-
  Ep_cp_TFs_umap_res_TFnum %>% filter(type == "Epromoter")


pdf("UMAP cluster promoters by TFs binding.pdf", width = 6, height = 6)

ggplot()+
  geom_point(Ep_cp_TFs_umap_res_TFnum1, aes(x = V1, y = V2, color = TFs_num), size = 0.8)+
  scale_color_gradient(low = "white", high = "gray")+
  # new_scale_color()+
  geom_point(Ep_cp_TFs_umap_res_TFnum2, aes(x = V1, y = V2, color = TFs_num), size = 0.8)+
  scale_color_gradient(values = c("white", "#377EB8"))+
  labs(title = "UMAP cluster promoters by TFs binding",
       x = "UMAP1", y = "UMAP2")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "top")

dev.off()



# pdf("UMAP cluster Epromoters by TFs binding.pdf", width = 6, height = 6)
pdf("UMAP cluster Epromoters by TFs binding density.pdf", width = 6, height = 6)

ggplot(Ep_cp_TFs_umap_res_TFnum, 
       aes(x = Ep_cp_TFs_umap_res_TFnum$V1, y = Ep_cp_TFs_umap_res_TFnum$V2, color = TFs_num))+
  geom_point(size = 0.3)+
  scale_color_gradient(low = "blue", high = "red")+
  labs(title = "UMAP cluster promoters by all TFs(ReMap) binding",
       x = "UMAP1", y = "UMAP2", color = "TFs binding\n density")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        # axis.ticks = element_blank(),
        axis.title = element_text(size=16),
        legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size=10))
# guides(color = guide_legend(override.aes = list(size=3)))

dev.off()





