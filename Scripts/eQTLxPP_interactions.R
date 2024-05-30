

#======================Epromoter VS control eQTL x PP=============================

#==========Epromoter==============

test1 <- fread("Epromoter-variants_GWAS-eQTL-PP-MPRA-TF-CRISPR_list-2023-08_common-target.txt", header = T)

eQTL_PP_common_target_df <- test1[test1$eQTL_effect != "proximal" & test1$Epromoter_gene != test1$eQTL_PP_common_target, ]
eQTL_PP_common_target_df1 <- test1[test1$eQTL_effect == "distal" & test1$Epromoter_gene != test1$eQTL_PP_common_target, ]
eQTL_PP_common_target_df2 <- test1[test1$eQTL_effect == "distal;proximal" & test1$Epromoter_gene != test1$eQTL_PP_common_target, ]

eQTL_PP_common_target_df$type <- "Ep.common_target"
eQTL_PP_common_target_df1$type <- "Ep.common_target1"
eQTL_PP_common_target_df2$type <- "Ep.common_target2"

test1_eQTL <- test1[!is.na(test1$eQTL_effect),]
test1_PP <- test1[!is.na(test1$PP_target),]

eQTL_target_only_GWAS <- test1_eQTL[!(test1_eQTL$rsID %in% eQTL_PP_common_target_df$rsID), "GWAS_trait"]
eQTL_target_only_GWAS$type <- "Ep.eQTL_target_only"

PP_target_only_GWAS <- test1_PP[!(test1_PP$rsID %in% eQTL_PP_common_target_df$rsID), "GWAS_trait"]
PP_target_only_GWAS$type <- "Ep.PP_target_only"

#also some intersection between eQTL_target_only_GWAS and PP_target_only_GWAS
#means have both eQTL and PP, but target different genes

#choose one
#common_target: eQTL distal,distal&proximal x PP
# eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df[, c("GWAS_trait", "type")])
#common_target: eQTL distal x PP
eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df1[, c("GWAS_trait", "type")])
#common_target: eQTL distal&proximal x PP
eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df2[, c("GWAS_trait", "type")])
#common_target: eQTL distal,distal&proximal x PP
eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df1[, c("GWAS_trait", "type")], eQTL_PP_common_target_df2[, c("GWAS_trait", "type")])

#count GWAS trait number
for (i in 1:nrow(eQTL_PP_intersect_GWAS)) {
  eQTL_PP_intersect_GWAS$GWAS_trait_num[i] <- length(unlist(strsplit(eQTL_PP_intersect_GWAS$GWAS_trait[i], ",|;")))
}

eQTL_PP_intersect_GWAS[eQTL_PP_intersect_GWAS$GWAS_trait_num > 10, "GWAS_trait_num"] <- 10

Ep.eQTL_PP_intersect_GWAS <- eQTL_PP_intersect_GWAS


# eQTL_PP_intersect_GWAS$type <- factor(eQTL_PP_intersect_GWAS$type, levels = c("eQTL_target_only", "common_target", "PP_target_only"))
# eQTL_PP_intersect_GWAS$type <- factor(eQTL_PP_intersect_GWAS$type, levels = c("eQTL_target_only", "PP_target_only", "common_target"))


#==========control==============

eQTL_PP_common_target_df <- cp_GWAS_eQTL_PP[cp_GWAS_eQTL_PP$eQTL_effect != "proximal" & cp_GWAS_eQTL_PP$V12 != cp_GWAS_eQTL_PP$eQTL_PP_common_target, ]
eQTL_PP_common_target_df1 <- cp_GWAS_eQTL_PP[cp_GWAS_eQTL_PP$eQTL_effect == "distal" & cp_GWAS_eQTL_PP$V12 != cp_GWAS_eQTL_PP$eQTL_PP_common_target, ]
eQTL_PP_common_target_df2 <- cp_GWAS_eQTL_PP[cp_GWAS_eQTL_PP$eQTL_effect == "distal;proximal" & cp_GWAS_eQTL_PP$V12 != cp_GWAS_eQTL_PP$eQTL_PP_common_target, ]

eQTL_PP_common_target_df$type <- "cp.common_target"
eQTL_PP_common_target_df1$type <- "cp.common_target1"
eQTL_PP_common_target_df2$type <- "cp.common_target2"


test1_eQTL <- cp_GWAS_eQTL_PP[!is.na(cp_GWAS_eQTL_PP$eQTL_effect),]
test1_PP <- cp_GWAS_eQTL_PP[!is.na(cp_GWAS_eQTL_PP$PP_target),]

eQTL_target_only_GWAS <- test1_eQTL[!(test1_eQTL$V4 %in% eQTL_PP_common_target_df$V4), "V13"]
eQTL_target_only_GWAS$type <- "cp.eQTL_target_only"

PP_target_only_GWAS <- test1_PP[!(test1_PP$V4 %in% eQTL_PP_common_target_df$V4), "V13"]
PP_target_only_GWAS$type <- "cp.PP_target_only"

#also some intersection between eQTL_target_only_GWAS and PP_target_only_GWAS
#means have both eQTL and PP, but target different genes

#choose one
#common_target: eQTL distal,distal&proximal x PP
# eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df[, c("V13", "type")])
#common_target: eQTL distal x PP
eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df1[, c("V13", "type")])
#common_target: eQTL distal&proximal x PP
eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df2[, c("V13", "type")])
#common_target: eQTL distal,distal&proximal x PP
eQTL_PP_intersect_GWAS <- rbind(eQTL_target_only_GWAS, PP_target_only_GWAS, eQTL_PP_common_target_df1[, c("V13", "type")], eQTL_PP_common_target_df2[, c("V13", "type")])

colnames(eQTL_PP_intersect_GWAS) <- c("GWAS_trait", "type")

#count GWAS trait number
for (i in 1:nrow(eQTL_PP_intersect_GWAS)) {
  eQTL_PP_intersect_GWAS$GWAS_trait_num[i] <- length(unlist(strsplit(eQTL_PP_intersect_GWAS$GWAS_trait[i], ",|;")))
}

eQTL_PP_intersect_GWAS[eQTL_PP_intersect_GWAS$GWAS_trait_num > 10, "GWAS_trait_num"] <- 10

cp.eQTL_PP_intersect_GWAS <- eQTL_PP_intersect_GWAS


#====merge Ep cp=======
Ep.cp.eQTL_PP_intersect_GWAS <- rbind(Ep.eQTL_PP_intersect_GWAS, cp.eQTL_PP_intersect_GWAS)

Ep.cp.eQTL_PP_intersect_GWAS$type <- factor(Ep.cp.eQTL_PP_intersect_GWAS$type, 
                                            levels = c("Ep.eQTL_target_only", "cp.eQTL_target_only",
                                                       "Ep.PP_target_only", "cp.PP_target_only",
                                                       "Ep.common_target1", "cp.common_target1",
                                                       "Ep.common_target2", "cp.common_target2"))



pdf("Epromoter-control plieotropy eQTL X PP interactions.pdf", width = 8, height = 4)

ggplot(Ep.cp.eQTL_PP_intersect_GWAS, aes(x = type, y = GWAS_trait_num, fill = type))+
  # ggplot(Ep_ctl_eQTL_GWAS_trait, aes(x = eQTL_type, y = GWAS_trait, fill = group))+
  geom_violin(width = 0.9)+
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  scale_y_continuous(limits = c(1,10.5), breaks = c(1,2,3,5,10))+
  scale_x_discrete(labels = c("             eQTL only\n            Ep/ctl", "", 
                              "               PP only\n                Ep/ctl", "", 
                              "               eQTLdistal & PP\n              Ep/ctl", "",
                              "               eQTLdistal+proximal & PP\n              Ep/ctl", ""))+
  scale_fill_manual(values = c("#377EB8", "gray", "#377EB8", "gray", "#377EB8", "gray", "#377EB8", "gray"))+
  geom_hline(yintercept = median_line, linetype = 2, color = "gray")+
  labs(title = "eQTL X PP interactions", y = "GWAS trait")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size=12, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16))+
  guides(fill = FALSE)+
  geom_signif(comparisons = list(c("Ep.eQTL_target_only", "cp.eQTL_target_only"),
                                 c("Ep.PP_target_only", "cp.PP_target_only"), 
                                 c("Ep.common_target1", "cp.common_target1"),
                                 c("Ep.common_target2", "cp.common_target2")), 
              y_position = 9, textsize = 4, linetype = "blank")

dev.off()



