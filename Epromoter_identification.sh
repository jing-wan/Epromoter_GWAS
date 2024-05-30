#!/bin/bash

#The shell command for Epromoter identification
#input files: enhancers collected from different datasets

#promoters (500bp) overlap with enhancers either 50% each other
bedtools intersect \
-a Ensembl_hg38_coding_transcript_TSS500_uniqTSS_rmTranscript.bed \
-b 28_STARR-seq_datasets_all_enhancers_list_merged.bed \
-wa -wb -f 0.5 -F 0.5 -e \
| cut -f 1-6,10-12 \
>Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15.bed

bedtools merge -i Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15.bed \
-c 4,5,6,7,8,9 -o distinct -delim "," \
| sort -k1,1V -k2,2n \
>Epromoters_ovTSS_Enh-TSS500-50pct_2023-02-15_merged.bed

