#!/usr/bin/bash

#===========Partitioned heritability==========

#Epromoter
for i in {1..22}
do
##Step 1: Creating an annot file
python ~/tools/ldsc/make_annot.py \
--bed-file Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15_merged.bed \
--bimfile ../../GRCh38/plink_files/1000G.EUR.hg38.${i}.bim \
--annot-file Epromoter/Epromoter_chr${i}.annot.gz

##Step 2: Computing LD scores with an annot file
python ~/tools/ldsc/ldsc.py \
--l2 \
--bfile ../../GRCh38/plink_files/1000G.EUR.hg38.${i} \
--ld-wind-cm 1 \
--annot Epromoter/Epromoter_chr${i}.annot.gz \
--thin-annot \
--out Epromoter/Epromoter_chr${i} \
--print-snps ../../GRCh38/hapmap3_snps/hapmap3.${i}.snp
done


#control promoter
for i in {1..22}
do

##Step 1: Creating an annot file
python ~/tools/ldsc/make_annot.py \
--bed-file control_promoters.hg38_coding_gene.5743merged.bed \
--bimfile ../../GRCh38/plink_files/1000G.EUR.hg38.${i}.bim \
--annot-file control_promoter/control_promoter_chr${i}.annot.gz

##Step 2: Computing LD scores with an annot file
python ~/tools/ldsc/ldsc.py \
--l2 \
--bfile ../../GRCh38/plink_files/1000G.EUR.hg38.${i} \
--ld-wind-cm 1 \
--annot control_promoter/control_promoter_chr${i}.annot.gz \
--thin-annot \
--out control_promoter/control_promoter_chr${i} \
--print-snps ../../GRCh38/hapmap3_snps/hapmap3.${i}.snp
done


#LDSC partitioned heritability 
for file in ../../GWAS_sumstats/all_sumstats/*.sumstats
do
sample=$(basename $file .sumstats)
python  ~/tools/ldsc/ldsc.py \
--h2 $file \
--ref-ld-chr Epromoter/Epromoter_chr@,control_promoter/control_promoter_chr@,../../GRCh38/baseline_v1.2/baseline.@ \
--w-ld-chr ../../GRCh38/weights/weights.hm3_noMHC. \
--overlap-annot \
--frqfile-chr ../../GRCh38/plink_files/1000G.EUR.hg38. \
--out LDSC_res_Ep_cp_baseline/$sample"_Epromoter_control"
done


