require(ggplot2)
require(reshape2)
require(tidyverse)


LDSC_res_files <- list.files(pattern = ".results")

Ep_cp_heritability_all <- data.frame()

for (i in LDSC_res_files) {
  # i <- "PASS_ADHD_Demontis2018_Epromoter_control.results"
  # i <- "UKB_460K.biochemistry_Albumin_Epromoter_control.results"
  LDSC_res <- read.table(i, header = T)
  GWAS_study <- strsplit(i, "PASS_|_Epromoter_control.results|UKB_460K.")[[1]][2]
  
  Ep_cp_heritability <- data.frame(GWAS_sumstats_study = GWAS_study,
                                   Epromoter_LDSC_Prop.SNPs = LDSC_res[LDSC_res$Category %in% "L2_0", "Prop._SNPs"],
                                   Epromoter_LDSC_Prop.h2 = LDSC_res[LDSC_res$Category %in% "L2_0", "Prop._h2"],
                                   Epromoter_LDSC_enrichment = LDSC_res[LDSC_res$Category %in% "L2_0", "Enrichment"],
                                   control_promoter_LDSC_Prop.SNPs = LDSC_res[LDSC_res$Category %in% "L2_1", "Prop._SNPs"],
                                   control_promoter_LDSC_Prop.h2 = LDSC_res[LDSC_res$Category %in% "L2_1", "Prop._h2"],
                                   control_promoter_LDSC_enrichment = LDSC_res[LDSC_res$Category %in% "L2_1", "Enrichment"],
                                   Enhancer_LDSC_Prop.SNPs = LDSC_res[LDSC_res$Category %in% "Enhancer_AnderssonL2_2", "Prop._SNPs"],
                                   Enhancer_LDSC_Prop.h2 = LDSC_res[LDSC_res$Category %in% "Enhancer_AnderssonL2_2", "Prop._h2"],
                                   Enhancer_LDSC_enrichment = LDSC_res[LDSC_res$Category %in% "Enhancer_AnderssonL2_2", "Enrichment"],
                                   Promoter_UCSC_LDSC_Prop.SNPs = LDSC_res[LDSC_res$Category %in% "Promoter_UCSCL2_2", "Prop._SNPs"],
                                   Promoter_UCSC_LDSC_Prop.h2 = LDSC_res[LDSC_res$Category %in% "Promoter_UCSCL2_2", "Prop._h2"],
                                   Promoter_UCSC_LDSC_enrichment = LDSC_res[LDSC_res$Category %in% "Promoter_UCSCL2_2", "Enrichment"],
                                   Coding_LDSC_Prop.SNPs = LDSC_res[LDSC_res$Category %in% "Coding_UCSCL2_2", "Prop._SNPs"],
                                   Coding_LDSC_Prop.h2 = LDSC_res[LDSC_res$Category %in% "Coding_UCSCL2_2", "Prop._h2"],
                                   Coding_LDSC_enrichment = LDSC_res[LDSC_res$Category %in% "Coding_UCSCL2_2", "Enrichment"])
  
  Ep_cp_heritability_all <- rbind(Ep_cp_heritability_all, Ep_cp_heritability)
  
}

Ep_cp_heritability_all_enrich_ranking <-
  Ep_cp_heritability_all[order(Ep_cp_heritability_all$Epromoter_LDSC_enrichment, decreasing = T), ]

write.table(Ep_cp_heritability_all_enrich_ranking, file = "../GWAS_sumstats_heritability_of_Epromoter_control_Enhancer_Promoter_Coding.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")



#======================GWAS heritability heatmap for each top10=================

require(pheatmap)
require(grid)

Ep_cp_heritability_all_c6 <- Ep_cp_heritability_all[, c("GWAS_sumstats_study", "Epromoter_LDSC_enrichment", "control_promoter_LDSC_enrichment", "Enhancer_LDSC_enrichment", "Promoter_UCSC_LDSC_enrichment", "Coding_LDSC_enrichment")]
colnames(Ep_cp_heritability_all_c6) <- c("GWAS_sumstats_study", "Epromoter", "control_promoter", "Enhancer_Andersson", "Promoter_UCSC", "Coding_UCSC")

heritability_top10_c7 <-
  rbind(Ep_cp_heritability_all_c6 %>% arrange(desc(Epromoter)) %>% slice(1:10) %>% mutate(category = "Epromoter"),
        Ep_cp_heritability_all_c6 %>% arrange(desc(control_promoter)) %>% slice(1:10) %>% mutate(category = "control_promoter"),
        Ep_cp_heritability_all_c6 %>% arrange(desc(Enhancer_Andersson)) %>% slice(1:10) %>% mutate(category = "Enhancer_Andersson"))

#merge common GWAS category
heritability_top10_merged <- 
  heritability_top10_c7 %>%
  group_by(GWAS_sumstats_study) %>%
  summarise(across(c("Epromoter", "control_promoter", "Enhancer_Andersson", "Promoter_UCSC", "Coding_UCSC"), ~unique(.)), 
            category = paste(category, collapse = ",")) %>%
  arrange(category, case_when(category == "Epromoter" ~ desc(Epromoter),
                              category == "control_promoter" ~ desc(control_promoter),
                              category == "Enhancer_Andersson" ~ desc(Enhancer_Andersson))) %>%
  arrange(match(category, c("Epromoter", "control_promoter", "Enhancer_Andersson", "Epromoter,Enhancer_Andersson", "control_promoter,Enhancer_Andersson")))


# heritability_top10_uniq <- heritability_top10_c7[!duplicated(heritability_top10_c7[, 1:6]), ]

heritability_top10_uniq <- as.data.frame(heritability_top10_merged)
rownames(heritability_top10_uniq) <- 1:nrow(heritability_top10_uniq)

#scale
heritability_top10_uniq[, 2:6][heritability_top10_uniq[, 2:6] > 100] <- 100
heritability_top10_uniq[, 2:6][heritability_top10_uniq[, 2:6] < -20] <- -20

heritability_top10_uniq[, 2:6] <- round(heritability_top10_uniq[, 2:6], digits = 1)

heritability_top10_uniq$category <- factor(heritability_top10_uniq$category, levels = unique(heritability_top10_uniq$category))

heritability_top10_uniq_pht <-
  pheatmap(heritability_top10_uniq[, 2:6], 
           cluster_rows = F, cluster_cols = F, show_rownames = T,
           breaks = c(seq(-20, 30, length.out=50), seq(31, 100, length.out=50)),
           legend_breaks = c(0, 50, 100),
           annotation_row = data.frame(category = as.factor(heritability_top10_uniq$category)),
           annotation_names_row = F, annotation_legend = T,
           gaps_row = c(6, 11, 15, 17, 20, 21),
           labels_row = heritability_top10_uniq$GWAS_sumstats_study,
           angle_col = 45, cellwidth = 40, cellheight = 20,
           main = "Top10 GWAS enriched in Epromoter, Control, Enhancer")


save_pheatmap_pdf <- function(x, filename, width=12, height=9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  grid.text("heritability=Prop.h2/Prop.SNPs", x=0.6, y=0.95, gp=gpar(fontsize=10))
  dev.off()
}

save_pheatmap_pdf(heritability_top10_uniq_pht, "../Top10_GWAS_heritability_in_Epromoter-Control-Enhancer_heatmap1.pdf")




