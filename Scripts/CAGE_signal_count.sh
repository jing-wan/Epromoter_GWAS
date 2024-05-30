#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --nodelist=bigmemorix
#SBATCH --mem=16000
#SBATCH -o myjob.o
#SBATCH -e myjob.e

cd $PWD

awk '{print $1,$2,$3}' Epromoters_5743merged.rm_bidirection.TSS_up500down500.bed \
| while read -r chr start end; do \
cmd="bigWigToWig ../../ctssTotalCounts.fwd.bw -chrom=$chr -start=$start -end=$end stdout >>ctssTotalCounts.fwd.Epromoters.wig" 
eval "$cmd"
done

awk '{print $1,$2,$3}' Epromoters_5743merged.rm_bidirection.TSS_up500down500.bed \
| while read -r chr start end; do \
cmd="bigWigToWig ../../ctssTotalCounts.rev.bw -chrom=$chr -start=$start -end=$end stdout >>ctssTotalCounts.rev.Epromoters.wig" 
eval "$cmd"
done

awk '{print $1,$2,$3}' control_promoters.5743merged.rm_bidirection.TSS_up500down500.bed \
| while read -r chr start end; do \
cmd="bigWigToWig ../../ctssTotalCounts.fwd.bw -chrom=$chr -start=$start -end=$end stdout >>ctssTotalCounts.fwd.control_promoters.wig"  
eval "$cmd"
done

awk '{print $1,$2,$3}' control_promoters.5743merged.rm_bidirection.TSS_up500down500.bed \
| while read -r chr start end; do \
cmd="bigWigToWig ../../ctssTotalCounts.rev.bw -chrom=$chr -start=$start -end=$end stdout >>ctssTotalCounts.rev.control_promoters.wig"  
eval "$cmd"
done

wig2bed < ctssTotalCounts.fwd.Epromoters.wig > ctssTotalCounts.fwd.Epromoters.bed
wig2bed < ctssTotalCounts.rev.Epromoters.wig > ctssTotalCounts.rev.Epromoters.bed
wig2bed < ctssTotalCounts.fwd.control_promoters.wig > ctssTotalCounts.fwd.control_promoters.bed
wig2bed < ctssTotalCounts.rev.control_promoters.wig > ctssTotalCounts.rev.control_promoters.bed


#calculate CAGE fwd/rev count in each promoter
bedtools intersect -a Epromoters_5743merged.rm_bidirection.TSS_up500down500.bed \
-b ctssTotalCounts.fwd.Epromoters.bed -wa -wb \
| bedtools merge -i - -d -1000 -c 4,5,6,7,8,9,14 \
-o distinct,distinct,distinct,distinct,distinct,distinct,sum \
>Epromoters_CAGE_fwd_total_count.bed

bedtools intersect -a Epromoters_5743merged.rm_bidirection.TSS_up500down500.bed \
-b ctssTotalCounts.rev.Epromoters.bed -wa -wb \
| bedtools merge -i - -d -1000 -c 4,5,6,7,8,9,14 \
-o distinct,distinct,distinct,distinct,distinct,distinct,sum \
>Epromoters_CAGE_rev_total_count.bed

bedtools intersect -a control_promoters.5743merged.rm_bidirection.TSS_up500down500.bed \
-b ctssTotalCounts.fwd.control_promoters.bed -wa -wb \
| bedtools merge -i - -d -1000 -c 4,5,6,11 \
-o distinct,distinct,distinct,sum \
>control_promoters_CAGE_fwd_total_count.bed

bedtools intersect -a control_promoters.5743merged.rm_bidirection.TSS_up500down500.bed \
-b ctssTotalCounts.rev.control_promoters.bed -wa -wb \
| bedtools merge -i - -d -1000 -c 4,5,6,11 \
-o distinct,distinct,distinct,sum \
>control_promoters_CAGE_rev_total_count.bed












