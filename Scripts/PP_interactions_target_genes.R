#!/usr/bin/Rscript

#==================count target genes of PP interactions================


all_target <- read.table("total_coding_gene_promoters_PP_interactions_target-number.txt", header = F)
Ep_target <- read.table("3datasets_PP_interactions_Ep_target-number.txt", header = F)
cp_target <- read.table("3datasets_PP_interactions_cp_target-number.txt", header = F)

#include 0 target
all_target[all_target$V10 == ".", "V10"] <- 0

#count target genes number after 3datasets merged, split by ,|;
for (i in 1:nrow(all_target)) {
  if (all_target$V10[i] == "0") {
    all_target$target_genes_num[i] <- 0
  } 
  else {
    target_genes <- unlist(strsplit(all_target$V10[i], ",|;"))
    all_target$target_genes_num[i] <- length(unique(target_genes))
    all_target$target_genes[i] <- paste(unique(target_genes), collapse = ";")
  }
}

for (i in 1:nrow(Ep_target)) {
  target_genes <- unlist(strsplit(Ep_target$V4[i], ",|;"))
  Ep_target$target_genes_num[i] <- length(unique(target_genes))
  Ep_target$target_genes[i] <- paste(unique(target_genes), collapse = ";")
}

for (i in 1:nrow(cp_target)) {
  target_genes <- unlist(strsplit(cp_target$V4[i], ",|;"))
  cp_target$target_genes_num[i] <- length(unique(target_genes))
  cp_target$target_genes[i] <- paste(unique(target_genes), collapse = ";")
}


Ep_target$type <- "Epromoter"
cp_target$type <- "control"
all_target$type <- "total"

Ep_target0 <- rbind(Ep_target[, c("target_genes_num", "type")], 
                    data.frame(target_genes_num = rep(0, 5743-nrow(Ep_target)),
                               type = rep("Epromoter", 5743-nrow(Ep_target))))
cp_target0 <- rbind(cp_target[, c("target_genes_num", "type")],
                    data.frame(target_genes_num = rep(0, 5743-nrow(cp_target)),
                               type = rep("control", 5743-nrow(cp_target))))
all_target0 <- all_target[, c("target_genes_num", "type")]

Ep_cp_all_target <- rbind(all_target0, Ep_target0, cp_target0)


Ep_cp_all_target[Ep_cp_all_target$target_genes_num > 30, "target_genes_num"] <- 30
Ep_cp_all_target <- Ep_cp_all_target[!is.na(Ep_cp_all_target$target_genes_num), ]
Ep_cp_all_target$type <- factor(Ep_cp_all_target$type, levels = c("total", "Epromoter", "control"))


pdf("3datastes_PP_interactions_Ep-cp-total_target_genes_0.pdf", width = 5, height = 4)

ggplot(Ep_cp_all_target, aes(x = type, y = target_genes_num, fill = type))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(limits = c(0, 30), breaks = c(0,5,8,10,15,20,25,30))+
  scale_fill_manual(values = c("#4DAF4A", "#377EB8", "gray"))+
  # scale_fill_manual(values = c("#377EB8", "gray", "#4DAF4A"))+
  geom_hline(yintercept = 8, linetype = 2, color = "gray")+
  labs(title = "promoter-promoter interactions", y = "PP target genes")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "total"),
                                 c("Epromoter", "control")),
              y_position = c(26, 26),
              textsize = 4, linetype = "blank")

dev.off()






