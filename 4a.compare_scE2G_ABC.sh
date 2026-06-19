
######################### OVERLAPPING REGULATORY REGIONS BETWEEN scE2G and ABC

module load bedtools/2.31.0

pairToPair -a data/a.compare_scE2G_ABC/abc_0hr_hg38.tsv \
-b data/a.compare_scE2G_ABC/scE2G_0hr_hg38.tsv -type both | bedtools overlap -i stdin -cols 2,3,12,13 > data/a.compare_scE2G_ABC/gene_link_tnf_0hr_overlap.txt

pairToPair -a data/a.compare_scE2G_ABC/abc_4hr_hg38.tsv \
-b data/a.compare_scE2G_ABC/scE2G_4hr_hg38.tsv -type both | bedtools overlap -i stdin -cols 2,3,12,13 > data/a.compare_scE2G_ABC/gene_link_tnf_4hr_overlap.txt

pairToPair -a data/a.compare_scE2G_ABC/abc_24hr_hg38.tsv \
-b data/a.compare_scE2G_ABC/scE2G_24hr_hg38.tsv -type both | bedtools overlap -i stdin -cols 2,3,12,13 > data/a.compare_scE2G_ABC/gene_link_tnf_24hr_overlap.txt


######################### MERGE TSV FILES
#I merge overlapping regulatory regions across conditions (TNF24hr,TNF4hr,NT) and report the linked genes by those regions
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/scE2G_only.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o collapse > data/a.compare_scE2G_ABC/merged_bed/scE2G_only_merged_w_gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} $12 != "promoter" {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/scE2G_only.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o collapse > data/a.compare_scE2G_ABC/merged_bed/scE2G_only_no_promoter_merged_w_gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/ABC_only.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o collapse  > data/a.compare_scE2G_ABC/merged_bed/ABC_only_merged_w_gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/scE2G_in_ABC.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o collapse  > data/a.compare_scE2G_ABC/merged_bed/scE2G_in_ABC_merged_w_gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/ABC_in_scE2G.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o collapse  > data/a.compare_scE2G_ABC/merged_bed/ABC_in_scE2G_merged_w_gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/scE2G.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o collapse  > data/a.compare_scE2G_ABC/merged_bed/scE2G_merged_w_gene.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");print $1,$2,$3,a[1]}' data/a.compare_scE2G_ABC/ABC.tsv |sort -k1,1 -k2,2n | bedtools merge -c 4 -o collapse  > data/a.compare_scE2G_ABC/merged_bed/ABC_merged_w_gene.bed
cat data/a.compare_scE2G_ABC/merged_bed/scE2G_in_ABC_merged_w_gene.bed data/a.compare_scE2G_ABC/merged_bed/ABC_in_scE2G_merged_w_gene.bed  |sort -k1,1 -k2,2n | bedtools merge  -c 4 -o collapse > data/a.compare_scE2G_ABC/merged_bed/ABC_scE2G_merged_w_gene.bed

#I merge overlapping regulatory regions across conditions (TNF24hr,TNF4hr,NT) and report the category of the region (promoter,intergenic,genic)
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");split(a[2],b,"|");print $1,$2,$3,b[1]}' data/a.compare_scE2G_ABC/scE2G_only.tsv | sort -k1,1 -k2,2n| bedtools merge -c 4 -o distinct > data/a.compare_scE2G_ABC/merged_bed/scE2G_only_merged_w_category.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");split(a[2],b,"|");print $1,$2,$3,b[1]}' data/a.compare_scE2G_ABC/ABC_only.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o distinct  > data/a.compare_scE2G_ABC/merged_bed/ABC_only_merged_w_category.bed
awk -F"\t" 'BEGIN{OFS="\t"} {split($7, a, "_");split(a[2],b,"|");print $1,$2,$3,b[1]}' data/a.compare_scE2G_ABC/scE2G_in_ABC.tsv |sort -k1,1 -k2,2n| bedtools merge -c 4 -o distinct  > data/a.compare_scE2G_ABC/merged_bed/scE2G_in_ABC_merged_w_category.bed #for ABE-scE2G

######################### MERGE BED FILES
#I merge overlapping regulatory regions across conditions (TNF24hr,TNF4hr,NT)
sort -k1,1 -k2,2n data/a.compare_scE2G_ABC/scE2G_only.bed | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/scE2G_only_merged.bed
sort -k1,1 -k2,2n data/a.compare_scE2G_ABC/ABC_only.bed | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/ABC_only_merged.bed
sort -k1,1 -k2,2n data/a.compare_scE2G_ABC/scE2G_in_ABC.bed | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/scE2G_in_ABC_merged.bed
sort -k1,1 -k2,2n data/a.compare_scE2G_ABC/ABC_in_scE2G.bed | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/ABC_in_scE2G_merged.bed
sort -k1,1 -k2,2n data/a.compare_scE2G_ABC/scE2G.bed | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/scE2G_merged.bed
sort -k1,1 -k2,2n data/a.compare_scE2G_ABC/ABC.bed | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/ABC_merged.bed
cat data/a.compare_scE2G_ABC/merged_bed/scE2G_in_ABC_merged.bed data/a.compare_scE2G_ABC/merged_bed/ABC_in_scE2G_merged.bed  |sort -k1,1 -k2,2n | bedtools merge > data/a.compare_scE2G_ABC/merged_bed/ABC_scE2G_merged.bed 
awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,$3}' data/a.compare_scE2G_ABC/merged_bed/scE2G_only_no_promoter_merged_w_gene.bed |sort -k1,1 -k2,2n| bedtools merge  > data/a.compare_scE2G_ABC/merged_bed/scE2G_only_no_promoter_merged.bed

##################### OVERLAP REGULATORY REGIONS WITH CAD AND DBP FINEMAPPED VARIANTS
#variants files format: chr    start   stop    id>PIP
module load bedtools/2.31.0
reg_region=()
while IFS= read -r line; do
    reg_region+=("$line")
done < data/a.compare_scE2G_ABC/regulatory_regions_list.txt

for r in "${reg_region[@]}"
do
bedtools intersect -a data/a.compare_scE2G_ABC/merged_bed/$r  \
-b data/a.compare_scE2G_ABC/finemapped_CAD_snps_PIP_hg38.bed -wa -wb -loj > data/a.compare_scE2G_ABC/overlap_variants/CAD_SNPS_${r}
done

for r in "${reg_region[@]}"
do
bedtools intersect -a data/a.compare_scE2G_ABC/merged_bed/$r  \
-b data/a.compare_scE2G_ABC/finemapped_DBP_snps_PIP_hg38.bed -wa -wb -loj > data/a.compare_scE2G_ABC/overlap_variants/DBP_SNPS_${r}
done

############################# REGULATORY REGIONS OVERLAPPING CAD FINEMAPPED VARIANTS + LINKED GENES
module load bedtools/2.31.0
bedtools intersect -a data/a.compare_scE2G_ABC/merged_bed/scE2G_only_merged_w_gene.bed -b data/a.compare_scE2G_ABC/finemapped_CAD_snps_PIP_hg38.bed -wa -wb > data/a.compare_scE2G_ABC/overlap_variants/scE2G_only_w_gene_CAD_SNPS.tsv
bedtools intersect -a data/a.compare_scE2G_ABC/merged_bed/scE2G_only_no_promoter_merged_w_gene.bed -b data/a.compare_scE2G_ABC/finemapped_CAD_snps_PIP_hg38.bed -wa -wb > data/a.compare_scE2G_ABC/overlap_variants/scE2G_only_no_promoter_w_gene_CAD_SNPS.tsv
bedtools intersect -a data/a.compare_scE2G_ABC/merged_bed/ABC_only_merged_w_gene.bed -b data/a.compare_scE2G_ABC/finemapped_CAD_snps_PIP_hg38.bed -wa -wb > data/a.compare_scE2G_ABC/overlap_variants/ABC_only_w_gene_CAD_SNPS.tsv 
bedtools intersect -a data/a.compare_scE2G_ABC/merged_bed/ABC_scE2G_merged_w_gene.bed -b data/a.compare_scE2G_ABC/finemapped_CAD_snps_PIP_hg38.bed -wa -wb > data/a.compare_scE2G_ABC/overlap_variants/ABC_scE2G_w_gene_CAD_SNPS.tsv 
