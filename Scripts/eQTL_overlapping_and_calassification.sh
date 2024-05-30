#!/usr/bin/bash

#overlap eQTLs and Epromoters
bedtools intersect \
-a Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15_merged.bed \
-b ../../4.1eQTL_list/eQTL_catalog_credible_sets_all_tissue_merged_sorted.bed \
-wa -wb >Epromoter_ov_eQTL_credible_sets.bed &

awk 'OFS="\t"{print $10,$11,$12,$13,$14,$15,$16,$17,$1,$2,$3,$4,$5,$6,$7,$8,$9}' \
Epromoter_ov_eQTL_credible_sets.bed \
| sort -k1,1V -k2,2n \
>Epromoter_ov_eQTL_credible_sets_byeQTL.bed


#merge eQTL and target gene TSS
#!/usr/bin/Rscript
# require(data.table)
# ensembl_tss <- fread("../Ensembl_hg38_all_genes_TSS.bed", header = F)
# eQTL_Epromoter <- fread("Epromoter_ov_eQTL_credible_sets_byeQTL.bed", header = F)
# #all.x = T, no.dups=F, allow.cartesian=T (keep duplicate record)
# eQTL_Epromoter_target_TSS <- merge(eQTL_Epromoter, ensembl_tss, by.x = "V4", by.y = "V4", all.x = T, no.dups=F, allow.cartesian=T)
# write.table(eQTL_Epromoter_target_TSS, file = "Epromoter_ov_eQTL_credible_sets_byeQTL_all_target_genes.bed", col.names=F, row.names=F, quote=F, sep = "\t")

#calculate distance between eQTL and target gene TSS 
#take <2kb as proximal, >=2kb as distal
awk 'OFS="\t"{print $2,$3,$4,$5,$18,$19,$20,$21,$22,$1,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' \
Epromoter_ov_eQTL_credible_sets_byeQTL_all_target_genes.bed \
| awk 'OFS="\t"{print $0,$2-$6}' \
| awk 'OFS="\t"{if (sqrt($23*$23)<2000) {print $0,$24="proximal"} else if (sqrt($23*$23)>=2000) {print $0,$24="distal"}}' \
| sort -k1,1V -k2,2n \
| awk '$5!~"NA"' \
>Epromoter_ov_eQTL_credible_sets_all-target-TSS-distance.bed


#eQTL-target_gene merge
bedtools merge -i Epromoter_ov_eQTL_credible_sets_all-target-TSS-distance.bed \
-c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 -o distinct -delim ";" -d -1 \
>Epromoter_ov_eQTL_credible_sets_all-target-TSS-distance_merged.bed









