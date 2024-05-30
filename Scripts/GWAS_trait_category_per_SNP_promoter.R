#!/usr/bin/Rscript

require(ggplot2)
require(ggsignif)
require(RColorBrewer)

#========GWAS trait per SNP==============

Ep_trait_per_SNP <- read.table("Epromoter.trait-per-SNP.txt", header = F)

op_trait_per_SNP <- read.table("control_promoter.trait-per-SNP.txt", header = F)

trait_per_SNP1 <- data.frame(trait_per_SNP = Ep_trait_per_SNP$V1, label = "Epromoter")
trait_per_SNP2 <- data.frame(trait_per_SNP = op_trait_per_SNP$V1, label = "control_promoter")

trait_per_SNP <- rbind(trait_per_SNP1, trait_per_SNP2)

trait_per_SNP[trait_per_SNP$trait_per_SNP >= 10, 1] <- 10

mean_value<- mean(trait_per_SNP[trait_per_SNP$label == "Epromoter", 1])

Ep_mean <- mean(trait_per_SNP[trait_per_SNP$label == "Epromoter", 1])
op_mean <- mean(trait_per_SNP[trait_per_SNP$label == "control_promoter", 1])

trait_per_SNP$label <- factor(trait_per_SNP$label, levels = c("Epromoter", "control_promoter"))

pdf("GWAS trait per SNP.pdf", width = 4.5, height = 4)

ggplot(trait_per_SNP, aes(x = label, y = trait_per_SNP, fill = label))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.shape = NA)+
  scale_y_continuous(breaks = c(1,2,5,10))+
  # scale_y_continuous(breaks = c(1,2,5,10,15))+
  # scale_x_discrete(labels = c("control_promoter\n3.24", "Epromoter\n3.42"))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = 2, linetype = 2, color = "gray")+
  labs(y = "GWAS trait per SNP")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control_promoter")), y_position = 9,
              textsize = 4, linetype = "blank")

dev.off()




#==============GWAS trait per promoter=================


Ep_trait_per_promoter <- read.table("GWAS-trait-per-Epromoter.bed", header = F)
op_trait_per_promoter <- read.table("GWAS-trait-per-control-promoter.bed", header = F)

trait_per_promoter1 <- data.frame(trait_per_promoter = Ep_trait_per_promoter$V6, 
                                  label = "Epromoter")
trait_per_promoter2 <- data.frame(trait_per_promoter = op_trait_per_promoter$V6, 
                                  label = "control_promoter")

trait_per_promoter <- rbind(trait_per_promoter1, trait_per_promoter2)

trait_per_promoter[trait_per_promoter$trait_per_promoter > 10, 1] <- 10

trait_per_promoter$label <- factor(trait_per_promoter$label, levels = c("Epromoter", "control_promoter"))

pdf("GWAS trait per promoter.pdf", width = 4.5, height = 4)

ggplot(trait_per_promoter, aes(x = label, y = trait_per_promoter, fill = label))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.shape = NA)+
  scale_y_continuous(breaks = c(0,2,3,5,10))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = 2, linetype = 2, color = "gray")+
  labs(y = "GWAS trait per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control_promoter")), y_position = 9, 
              textsize = 4, linetype = "blank")

dev.off()



#==============GWAS category per SNP=================

Ep_trait_category <- fread("Epromoter_LD_SNP_EFO-trait_parent-trait.txt", header = F)
cp_trait_category <- fread("control_promoter_LD_SNP_EFO-trait_parent-trait.txt", header = F)


colnames(Ep_trait_category) <- c("rsID", "trait", "traitID", "category", "categoryID")
colnames(cp_trait_category) <- c("rsID", "trait", "traitID", "category", "categoryID")

Ep_trait_category[is.na(Ep_trait_category$category), "trait"] <- "other"
cp_trait_category[is.na(cp_trait_category$category), "trait"] <- "other"
Ep_trait_category[is.na(Ep_trait_category$category), "category"] <- "other"
cp_trait_category[is.na(cp_trait_category$category), "category"] <- "other"


# Count unique traits and categories for each SNP
Ep_trait_category_count <- 
  aggregate(cbind(Traits = trait, Categories = category) ~ rsID, 
            data = Ep_trait_category, 
            FUN = function(x) length(unique(x)))

cp_trait_category_count <- 
  aggregate(cbind(Traits = trait, Categories = category) ~ rsID, 
            data = cp_trait_category, 
            FUN = function(x) length(unique(x)))

Ep_trait_category_count <- Ep_trait_category_count[Ep_trait_category_count$Traits >= 2, ]
cp_trait_category_count <- cp_trait_category_count[cp_trait_category_count$Traits >= 2, ]


trait_per_promoter1 <- data.frame(trait_per_promoter = Ep_trait_category_count$Categories, 
                                  label = "Epromoter")
trait_per_promoter2 <- data.frame(trait_per_promoter = cp_trait_category_count$Categories, 
                                  label = "control_promoter")

trait_per_promoter <- rbind(trait_per_promoter1, trait_per_promoter2)

trait_per_promoter$label <- factor(trait_per_promoter$label, levels = c("Epromoter", "control_promoter"))

pdf("GWAS Categories per SNP.pdf", width = 4.5, height = 4)

ggplot(trait_per_promoter, aes(x = label, y = trait_per_promoter, fill = label))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.shape = NA)+
  scale_y_continuous(breaks = c(1,2,3,5,10,15,18))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = 2, linetype = 2, color = "gray")+
  labs(y = "GWAS Categories per SNP")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control_promoter")), y_position = 12, 
              textsize = 4, linetype = "blank")

dev.off()




#=====================GWAS category per promoter==========================

require(data.table)

Ep_snp_gwas_info <- fread("Epromoter_GWAS_EFO_unique.bed", header = F)
efo_trait <- fread("trait_mappings", header = T)

Ep_snp_gwas_info <- unique(Ep_snp_gwas_info)
Ep_snp_gwas_info <- Ep_snp_gwas_info[order(Ep_snp_gwas_info$V5), ]
efo_trait <- unique(efo_trait[, 2:5])
efo_trait <- efo_trait[order(efo_trait$'EFO URI'), ]

#replace space in Parent term with "_"
efo_trait$'EFO term' <- gsub(" ", "_", efo_trait$'EFO term')
efo_trait$'Parent term' <- gsub(" ", "_", efo_trait$'Parent term')

#merge with parent trait
Ep_snp_gwas_info_efo_parent <- merge(Ep_snp_gwas_info, efo_trait, by.x = "V5", by.y = "EFO URI", all.x = T)
#order
Ep_snp_gwas_info_efo_parent_order <- Ep_snp_gwas_info_efo_parent[order(V1,V2), c(2,3,4,5,1,6,7,8)]

write.table(Ep_snp_gwas_info_efo_parent_order, file = "Epromoter_EFO-trait_parent-trait_byEp.bed", 
            sep = "\t", col.names = F, row.names = F, quote = F)


#control

Ep_snp_gwas_info <- fread("control_GWAS_EFO_unique.bed", header = F)
efo_trait <- fread("trait_mappings", header = T)

Ep_snp_gwas_info <- unique(Ep_snp_gwas_info)
Ep_snp_gwas_info <- Ep_snp_gwas_info[order(Ep_snp_gwas_info$V5), ]
efo_trait <- unique(efo_trait[, 2:5])
efo_trait <- efo_trait[order(efo_trait$'EFO URI'), ]

#replace space in Parent term with "_"
efo_trait$'EFO term' <- gsub(" ", "_", efo_trait$'EFO term')
efo_trait$'Parent term' <- gsub(" ", "_", efo_trait$'Parent term')

#merge with parent trait
Ep_snp_gwas_info_efo_parent <- merge(Ep_snp_gwas_info, efo_trait, by.x = "V5", by.y = "EFO URI", all.x = T)
#order
Ep_snp_gwas_info_efo_parent_order <- Ep_snp_gwas_info_efo_parent[order(V1,V2), c(2,3,4,5,1,6,7,8)]

write.table(Ep_snp_gwas_info_efo_parent_order, file = "control_EFO-trait_parent-trait_bycp.bed", 
            sep = "\t", col.names = F, row.names = F, quote = F)


#=====================GWAS category per promoter============

byEp_trait_category <- read.table("Epromoter_EFO-trait_parent-trait_byEp.bed", header = F)
bycp_trait_category <- read.table("control_EFO-trait_parent-trait_bycp.bed", header = F)

byEp_trait_category$promoterID <-
  paste(byEp_trait_category$V1, byEp_trait_category$V2, byEp_trait_category$V3, sep = "_")

bycp_trait_category$promoterID <-
  paste(bycp_trait_category$V1, bycp_trait_category$V2, bycp_trait_category$V3, sep = "_")


colnames(byEp_trait_category)[5:9] <- c("trait", "traitID", "category", "categoryID", "promoterID")
colnames(bycp_trait_category)[5:9] <- c("trait", "traitID", "category", "categoryID", "promoterID")

byEp_trait_category[is.na(byEp_trait_category$category), "trait"] <- "other"
bycp_trait_category[is.na(bycp_trait_category$category), "trait"] <- "other"
byEp_trait_category[is.na(byEp_trait_category$category), "category"] <- "other"
bycp_trait_category[is.na(bycp_trait_category$category), "category"] <- "other"


# Count unique traits and categories for each SNP
byEp_trait_category_count <- 
  aggregate(cbind(Traits = trait, Categories = category) ~ promoterID, 
            data = byEp_trait_category, 
            FUN = function(x) length(unique(x)))

bycp_trait_category_count <- 
  aggregate(cbind(Traits = trait, Categories = category) ~ promoterID, 
            data = bycp_trait_category, 
            FUN = function(x) length(unique(x)))


byEp_trait_category_count <- byEp_trait_category_count[byEp_trait_category_count$Traits >= 2, ]
bycp_trait_category_count <- bycp_trait_category_count[bycp_trait_category_count$Traits >= 2, ]


trait_per_promoter1 <- data.frame(trait_per_promoter = byEp_trait_category_count$Categories, 
                                  label = "Epromoter")
trait_per_promoter2 <- data.frame(trait_per_promoter = bycp_trait_category_count$Categories, 
                                  label = "control_promoter")

trait_per_promoter <- rbind(trait_per_promoter1, trait_per_promoter2)


trait_per_promoter$label <- factor(trait_per_promoter$label, levels = c("Epromoter", "control_promoter"))

pdf("GWAS Categories per promoter.pdf", width = 4.5, height = 4)

ggplot(trait_per_promoter, aes(x = label, y = trait_per_promoter, fill = label))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.05, outlier.shape = NA)+
  scale_y_continuous(breaks = c(1,2,3,5,10,15,18))+
  scale_x_discrete(labels = c("Epromoter", "control"))+
  scale_fill_manual(values = c("#377EB8", "gray"))+
  geom_hline(yintercept = 3, linetype = 2, color = "gray")+
  labs(y = "GWAS Categories per promoter")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Epromoter", "control_promoter")), y_position = 12, 
              textsize = 4, linetype = "blank")

dev.off()



