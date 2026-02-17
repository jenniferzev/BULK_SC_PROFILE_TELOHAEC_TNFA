library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(VennDiagram)
library(rsample)
library(ggridges)

############################### COMPARE OVERLAPPING AND NON OVERLAPPING scE2G AND ABC REGIONS

#see section OVERLAPPING REGULATORY REGIONS BETWEEN scE2G and ABC in 4a.compare_scE2G_ABC.sh

abc_0hr <- read.table("data/compare_scE2G_ABC/abc_0hr_hg38.tsv")
abc_0hr$V7 <- paste0(abc_0hr$V7,"_0hr")
abc_4hr <- read.table("data/compare_scE2G_ABC/abc_4hr_hg38.tsv")
abc_4hr$V7 <- paste0(abc_4hr$V7,"_4hr")
abc_24hr <- read.table("data/compare_scE2G_ABC/abc_24hr_hg38.tsv")
abc_24hr$V7 <- paste0(abc_24hr$V7,"_24hr")
abc <- rbind(abc_0hr,abc_4hr,abc_24hr)
colnames(abc) <- c("region_chr","region_start","region_end","gene_chr","gene_start","gene_end","gene_info","score","region_strand","gene_strand")
abc$condition <- ifelse(grepl("_0hr",abc$gene_info),"teloHAEC_TNFa_NT", ifelse(grepl("_4hr",abc$gene_info),"teloHAEC_TNFa_4hr","teloHAEC_TNFa_24hr"))

scE2G_0hr <- read.table("data/compare_scE2G_ABC/scE2G_0hr_hg38.tsv")
scE2G_0hr$V7 <- paste0(scE2G_0hr$V7,"_0hr")
scE2G_4hr <- read.table("data/compare_scE2G_ABC/scE2G_4hr_hg38.tsv")
scE2G_4hr$V7 <- paste0(scE2G_4hr$V7,"_4hr")
scE2G_24hr <- read.table("data/compare_scE2G_ABC/scE2G_24hr_hg38.tsv")
scE2G_24hr$V7 <- paste0(scE2G_24hr$V7,"_24hr")
scE2G <- rbind(scE2G_0hr,scE2G_4hr,scE2G_24hr)
colnames(scE2G) <- c("region_chr","region_start","region_end","gene_chr","gene_start","gene_end","gene_info","score","region_strand","gene_strand")
scE2G$condition <- ifelse(grepl("_0hr",scE2G$gene_info),"teloHAEC_TNFa_NT", ifelse(grepl("_4hr",scE2G$gene_info),"teloHAEC_TNFa_4hr","teloHAEC_TNFa_24hr"))

link_overlap_0hr <- read.table("data/compare_scE2G_ABC/gene_link_tnf_0hr_overlap.txt")
link_overlap_0hr$V7 <- paste0(link_overlap_0hr$V7,"_0hr")
link_overlap_0hr$V17 <- paste0(link_overlap_0hr$V17,"_0hr")
link_overlap_4hr <- read.table("data/compare_scE2G_ABC/gene_link_tnf_4hr_overlap.txt")
link_overlap_4hr$V7 <- paste0(link_overlap_4hr$V7,"_4hr")
link_overlap_4hr$V17 <- paste0(link_overlap_4hr$V17,"_4hr")
link_overlap_24hr <- read.table("data/compare_scE2G_ABC/gene_link_tnf_24hr_overlap.txt")
link_overlap_24hr$V7 <- paste0(link_overlap_24hr$V7,"_24hr")
link_overlap_24hr$V17 <- paste0(link_overlap_24hr$V17,"_24hr")
link_overlap <- rbind(link_overlap_0hr,link_overlap_4hr,link_overlap_24hr)
colnames(link_overlap) <- c(paste("abc",colnames(abc)[1:length(colnames(abc))-1], sep="_"),paste("scE2G",colnames(scE2G)[1:length(colnames(scE2G))-1], sep="_"),"overlap")
link_overlap$condition <- ifelse(grepl("_0hr",link_overlap$scE2G_gene_info),"teloHAEC_TNFa_NT", ifelse(grepl("_4hr",link_overlap$scE2G_gene_info),"teloHAEC_TNFa_4hr","teloHAEC_TNFa_24hr"))


scE2G$category <- as.factor(sub(".*_(.*?)\\|.*", "\\1",scE2G$gene_info))
abc$category <- as.factor(sub(".*_(.*?)\\|.*", "\\1",abc$gene_info))

#Filter abc-scE2G reg regions pairs
link_overlap <- link_overlap[link_overlap$overlap >= 250,]

abc$in_overlap <- ifelse(abc$gene_info %in% link_overlap$abc_gene_info, 1, 0)
scE2G$in_overlap <- ifelse(scE2G$gene_info %in% link_overlap$scE2G_gene_info, 1, 0)

#Count proportion of links in overlap between in scE2G and abc for each treatment time
abc_overlap_count <- table(abc$in_overlap,abc$condition)
scE2G_overlap_count <- table(scE2G$in_overlap,scE2G$condition)
overlap_count <-  data.frame(
  dataset = c(rep("abc",6),rep("scE2G",6)),
  treatment=c(rep("24h",2),rep("4h",2),rep("NT",2),rep("24h",2),rep("4h",2),rep("NT",2)),
  overlap = rep(c("0", "1"),6),
  count = c(as.vector(abc_overlap_count),as.vector(scE2G_overlap_count))
)
overlap_count$overlap <- factor(overlap_count$overlap, levels = rev(c("0", "1")))

svg("figures/compare_scE2G_ABC/scE2G_ABC_in_overlap_count.svg")
ggplot(overlap_count , aes(x = dataset, y = count, fill = overlap)) +
  facet_wrap(~ treatment) +
  geom_bar(stat = "identity") +
  labs(
    title = "Overlap counts per dataset",
    x = "Dataset",
    y = "Count",
    fill = "Overlap"
  ) +
  theme_minimal()
dev.off()


#Overlapping and specific regulatory regions for each treatment condition
scE2G_0hr$in_overlap <- ifelse(scE2G_0hr$V7 %in% link_overlap$scE2G_gene_info, 1, 0)
scE2G_4hr$in_overlap <- ifelse(scE2G_4hr$V7 %in% link_overlap$scE2G_gene_info, 1, 0)
scE2G_24hr$in_overlap <- ifelse(scE2G_24hr$V7 %in% link_overlap$scE2G_gene_info, 1, 0)

abc_0hr$in_overlap <- ifelse(abc_0hr$V7 %in% link_overlap$abc_gene_info, 1, 0)
abc_4hr$in_overlap <- ifelse(abc_4hr$V7 %in% link_overlap$abc_gene_info, 1, 0)
abc_24hr$in_overlap <- ifelse(abc_24hr$V7 %in% link_overlap$abc_gene_info, 1, 0)


#Write overlapping and specific regulatory regions
scE2G_only <- scE2G[!scE2G$gene_info %in% link_overlap$scE2G_gene_info, ]
scE2G_in_abc <- scE2G[scE2G$gene_info %in% link_overlap$scE2G_gene_info, ]
abc_only <- abc[!abc$gene_info %in% link_overlap$abc_gene_info, ]
abc_in_scE2G <- abc[abc$gene_info %in% link_overlap$abc_gene_info, ]

write.table(x=scE2G_only,file="data/compare_scE2G_ABC/scE2G_only.tsv",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=scE2G_in_abc,file="data/compare_scE2G_ABC/scE2G_in_ABC.tsv",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=abc_only,file="data/compare_scE2G_ABC/ABC_only.tsv",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=abc_in_scE2G,file="data/compare_scE2G_ABC/ABC_in_scE2G.tsv",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=scE2G,file="data/compare_scE2G_ABC/scE2G.tsv",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=abc,file="data/compare_scE2G_ABC/ABC.tsv",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)

write.table(x=scE2G_only[,1:3],file="data/compare_scE2G_ABC/scE2G_only.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=scE2G_in_abc[,1:3],file="data/compare_scE2G_ABC/scE2G_in_ABC.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=abc_only[,1:3],file="data/compare_scE2G_ABC/ABC_only.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=abc_in_scE2G[,1:3],file="data/compare_scE2G_ABC/ABC_in_scE2G.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=scE2G[,1:3],file="data/compare_scE2G_ABC/scE2G.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(x=abc[,1:3],file="data/compare_scE2G_ABC/ABC.bed",sep="\t",row.names=FALSE,quote=FALSE,col.names=FALSE)

############################### COMPARE scE2G SCORE DISTRIBUTION PER CATEGORY

svg("figures/compare_scE2G_ABC/scE2G_score_distrib_per_category.svg")
ggplot(scE2G, aes(x = score, y = category, fill = category)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(
    x="scE2G Score") +
  theme_minimal()
dev.off()

scE2G$distance_to_tss <- pmin(abs(scE2G$region_start - scE2G$gene_start), abs(scE2G$region_end - scE2G$gene_start))
svg("figures/compare_scE2G_ABC/scE2G_dist_to_tss_per_category-log10.svg")
ggplot(scE2G, aes(x = log10(distance_to_tss), y = category, fill = category)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(
  x="Distance to TSS (in bp) (log10)") +
  theme_minimal()
dev.off()

############################### PLOT THE NUMBER OF REGULATOREY REGION AFTER MERGING REDUNDANT REGULATORY REGIONS FOUND IN DIFFERENT TREATMENT TIME AND BETWEEN SCE2G AND ABC

#see section MERGE TSV FILES in 4a.compare_scE2G_ABC.sh

scE2G_only_merged <- read.table("data/compare_scE2G_ABC/merged_bed/scE2G_only_merged_w_gene.bed")
ABC_only_merged <- read.table("data/compare_scE2G_ABC/merged_bed/ABC_only_merged_w_gene.bed")
ABC_scE2G_merged <- read.table("data/compare_scE2G_ABC/merged_bed/ABC_scE2G_merged_w_gene.bed")

region_count <- data.frame(label=c("scE2G_only_merged","ABC_only_merged","ABC_scE2G_merged"),n_regions=c(nrow(scE2G_only_merged),nrow(ABC_only_merged),nrow(ABC_scE2G_merged)))
svg("figures/compare_scE2G_ABC/scE2G_ABC_in_overlap_count_merged.svg")
ggplot(region_count, aes(x = label, y = n_regions)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of regions",
    x = "Dataset",
    y = "Count"
  ) +
  theme_minimal()
dev.off()

############################## PLOT THE NUMBER OF EACH REGULATORY REGIONS CATEGORY (PROMOTER, GENIC, INTERGENIC), FOR scE2G-only ONLY

#see section MERGE TSV FILES in 4a.compare_scE2G_ABC.sh

scE2G_only_merged <- read.table("data/a.compare_scE2G_ABC/merged_bed/scE2G_only_merged_w_category.bed")
scE2G_only_category_count <- as.data.frame(table(scE2G_only_merged$V4))
scE2G_only_category_count <- scE2G_only_category_count[scE2G_only_category_count$Var1 %in% c("promoter","genic","intergenic"), ]
scE2G_only_category_count$scE2G_only <- "scE2G-only"
svg("figures/compare_scE2G_ABC/scE2G_only_regions_category_distrib.svg")
ggplot(scE2G_only_category_count, aes(x = scE2G_only, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of regions",
    x = "Category",
    y = "Count"
  ) +
  theme_minimal()
dev.off()


########################################## CAUSAL SIGNAL DENSITY, CAD

#see section OVERLAP REGULATORY REGIONS WITH CAD AND DBP FINEMAPPED VARIANTS in 4a.compare_scE2G_ABC.sh

r_regions <- c('ABC_merged.bed','ABC_only_merged.bed','ABC_scE2G_merged.bed','scE2G_merged.bed','scE2G_only_merged.bed')

ci_list <- list()
estimate_list <- list()
for(r in r_regions){
  cad_snps_reg_regions <- read.table(paste0("data/compare_scE2G_ABC/overlap_variants/CAD_SNPS_",r))
  cad_snps_reg_regions <- cad_snps_reg_regions[,c(1:3,7)]
  colnames(cad_snps_reg_regions) <- c("chr","start","end","snp_pip")
  cad_snps_reg_regions$snp_pip <- as.numeric(gsub(".*>","",cad_snps_reg_regions$snp_pip))

  cad_snps_reg_regions <- cad_snps_reg_regions %>%
    group_by(chr, start, end) %>%
    summarise(sum_pip = sum(snp_pip, na.rm = TRUE), .groups = "drop")

  cad_snps_reg_regions$coverage <- cad_snps_reg_regions$end - cad_snps_reg_regions$start

  genomic_coverage <- sum(cad_snps_reg_regions$coverage) / 3299210039 * 100

  data_sample <- bootstraps(cad_snps_reg_regions,times=1000)
  data_sample <- as.data.frame(data_sample$splits)
  data_sample <- data_sample[,grepl("sum_pip",colnames(as.data.frame(data_sample)))]
  sum_pip_estimates <- colSums(data_sample)
  csd_estimates <- sum_pip_estimates/genomic_coverage
  ci <- quantile(csd_estimates, c(0.025, 0.975))
  print(r)
  print(ci)
  ci_list[[r]] <- ci
  estimate_list[[r]] <- sum(cad_snps_reg_regions$sum_pip)/genomic_coverage
}

df <- data.frame(
  regions = r_regions,
  csd     = unlist(estimate_list[r_regions]),
  lower   = sapply(ci_list[r_regions], `[`, 1),
  upper   = sapply(ci_list[r_regions], `[`, 2)
)


svg("figures/compare_scE2G_ABC/cad_snps_causal_signal_density_CI95.svg")
ggplot(df, aes(x = reorder(regions,-csd), y = csd)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2)+
  labs(x="Regulatory regions", y = "sum(PIP) / Genome Coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(limits = c(NA, 42))
dev.off()

########################################## CAUSAL SIGNAL DENSITY, DBP

#see section OVERLAP REGULATORY REGIONS WITH CAD AND DBP FINEMAPPED VARIANTS in 4a.compare_scE2G_ABC.sh

r_regions <- c('ABC_merged.bed','ABC_only_merged.bed','ABC_scE2G_merged.bed','scE2G_merged.bed','scE2G_only_merged.bed')

ci_list <- list()
estimate_list <- list()
for(r in r_regions){
  cad_snps_reg_regions <- read.table(paste0("data/compare_scE2G_ABC/overlap_variants/DBP_SNPS_",r))
  cad_snps_reg_regions <- cad_snps_reg_regions[,c(1:3,7)]
  colnames(cad_snps_reg_regions) <- c("chr","start","end","snp_pip")
  cad_snps_reg_regions$snp_pip <- as.numeric(gsub(".*>","",cad_snps_reg_regions$snp_pip))

  cad_snps_reg_regions <- cad_snps_reg_regions %>%
    group_by(chr, start, end) %>%
    summarise(sum_pip = sum(snp_pip, na.rm = TRUE), .groups = "drop")

  cad_snps_reg_regions$coverage <- cad_snps_reg_regions$end - cad_snps_reg_regions$start

  genomic_coverage <- sum(cad_snps_reg_regions$coverage) / 3299210039 * 100

  data_sample <- bootstraps(cad_snps_reg_regions,times=1000)
  data_sample <- as.data.frame(data_sample$splits)
  data_sample <- data_sample[,grepl("sum_pip",colnames(as.data.frame(data_sample)))]
  sum_pip_estimates <- colSums(data_sample)
  csd_estimates <- sum_pip_estimates/genomic_coverage
  ci <- quantile(csd_estimates, c(0.025, 0.975))
  print(r)
  print(ci)
  ci_list[[r]] <- ci
  estimate_list[[r]] <- sum(cad_snps_reg_regions$sum_pip)/genomic_coverage
}

df <- data.frame(
  regions = r_regions,
  csd     = unlist(estimate_list[r_regions]),
  lower   = sapply(ci_list[r_regions], `[`, 1),
  upper   = sapply(ci_list[r_regions], `[`, 2)
)
mean(df$csd)#9.008526

svg("figures/compare_scE2G_ABC/dbp_snps_causal_signal_density_CI95.svg")
ggplot(df, aes(x = reorder(regions,-csd), y = csd)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2)+
  labs(x="Regulatory regions", y = "sum(PIP) / Genome Coverage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(limits = c(NA, 42))
dev.off()


########################################## EXPLORE OPEN TARGET GENES LINKED BY REGULATORY REGIONS OVERLAPPING CAD SNPS
open_target <- read.table("data/compare_scE2G_ABC/OT-EFO_0001645-associated-targets-06_06_2025-v25_03.tsv",header=TRUE,sep="\t")
open_target <- unique(open_target)  #remove duplicated genes
open_target <- open_target[!duplicated(open_target$symbol),] 

scE2G_only_region_genes <- read.table("data/compare_scE2G_ABC/overlap_variants/scE2G_only_w_gene_CAD_SNPS.tsv",header=FALSE)
scE2G_only_no_promoter_region_genes <- read.table("data/compare_scE2G_ABC/overlap_variants/scE2G_only_no_promoter_w_gene_CAD_SNPS.tsv",header=FALSE)
ABC_only_region_genes <- read.table("data/compare_scE2G_ABC/overlap_variants/ABC_only_w_gene_CAD_SNPS.tsv",header=FALSE)
ABC_scE2G_region_genes  <- read.table("data/compare_scE2G_ABC/overlap_variants/ABC_scE2G_w_gene_CAD_SNPS.tsv",header=FALSE)

colnames(scE2G_only_region_genes) <- c("region_chr","region_start","region_end","gene","snp_chr","snp_start","snp_stop","snp_id")
colnames(scE2G_only_no_promoter_region_genes) <- c("region_chr","region_start","region_end","gene","snp_chr","snp_start","snp_stop","snp_id")
colnames(ABC_only_region_genes) <- c("region_chr","region_start","region_end","gene","snp_chr","snp_start","snp_stop","snp_id")
colnames(ABC_scE2G_region_genes) <- c("region_chr","region_start","region_end","gene","snp_chr","snp_start","snp_stop","snp_id")

scE2G_only_region_genes$region_id <- paste(scE2G_only_region_genes$region_chr,scE2G_only_region_genes$region_start,scE2G_only_region_genes$region_end,sep="_")
scE2G_only_no_promoter_region_genes$region_id <- paste(scE2G_only_no_promoter_region_genes$region_chr,scE2G_only_no_promoter_region_genes$region_start,scE2G_only_no_promoter_region_genes$region_end,sep="_")
ABC_only_region_genes$region_id <- paste(ABC_only_region_genes$region_chr,ABC_only_region_genes$region_start,ABC_only_region_genes$region_end,sep="_")
ABC_scE2G_region_genes$region_id <- paste(ABC_scE2G_region_genes$region_chr,ABC_scE2G_region_genes$region_start,ABC_scE2G_region_genes$region_end,sep="_")

scE2G_only_region_genes<- scE2G_only_region_genes %>%
  separate_rows(gene, sep = ",")
scE2G_only_region_genes <- unique(scE2G_only_region_genes)
scE2G_only_region_genes$OT <- ifelse(scE2G_only_region_genes$gene %in% open_target$symbol,1,0)

scE2G_only_no_promoter_region_genes <- scE2G_only_no_promoter_region_genes %>%
  separate_rows(gene, sep = ",")
scE2G_only_no_promoter_region_genes <- unique(scE2G_only_no_promoter_region_genes)
scE2G_only_no_promoter_region_genes$OT <- ifelse(scE2G_only_no_promoter_region_genes$gene %in% open_target$symbol,1,0)

ABC_only_region_genes <- ABC_only_region_genes %>%
  separate_rows(gene, sep = ",")
ABC_only_region_genes<- unique(ABC_only_region_genes)
ABC_only_region_genes$OT <- ifelse(ABC_only_region_genes$gene %in% open_target$symbol,1,0)

ABC_scE2G_region_genes <- ABC_scE2G_region_genes  %>%
  separate_rows(gene, sep = ",")
ABC_scE2G_region_genes  <- unique(ABC_scE2G_region_genes )
ABC_scE2G_region_genes$OT <- ifelse(ABC_scE2G_region_genes$gene %in% open_target$symbol,1,0)


#ALL OT and non OT genes
all_genes <- data.frame(gene = unique(c(as.character(scE2G_only_region_genes$gene),as.character(scE2G_only_no_promoter_region_genes$gene), as.character(ABC_only_region_genes$gene), as.character(ABC_scE2G_region_genes$gene))))
all_genes$scE2G_ABC_region <- ifelse(all_genes$gene %in% ABC_scE2G_region_genes$gene,1,0)
all_genes$scE2G_only_no_promoter_region <- ifelse(all_genes$gene %in% scE2G_only_no_promoter_region_genes$gene,1,0)
all_genes$scE2G_only_region <- ifelse(all_genes$gene %in% scE2G_only_region_genes$gene,1,0)
all_genes$ABC_only_region <- ifelse(all_genes$gene %in% ABC_only_region_genes$gene,1,0)
all_genes <- merge(all_genes,open_target[,c("symbol","globalScore")],by.x="gene",by.y="symbol",all.x=TRUE) 

#number of regions overlapping a CAD snps and linking an OT gene or non OT gene
region_overlap_cad_snp_count_OT_gene <- data.frame(n_region = c(length(unique(scE2G_only_region_genes$region_id[scE2G_only_region_genes$OT == 1])),length(unique(scE2G_only_no_promoter_region_genes$region_id[scE2G_only_no_promoter_region_genes$OT == 1])),length(unique(ABC_only_region_genes$region_id[ABC_only_region_genes$OT == 1])),length(unique((ABC_scE2G_region_genes$region_id[ABC_scE2G_region_genes$OT == 1])))),
region = c("scE2G-only","scE2G-only-no-promoter","ABC-only","scE2G-ABC"))

#   n_region                 region
# 1      137             scE2G-only
# 2       91 scE2G-only-no-promoter
# 3       42               ABC-only
# 4       53              scE2G-ABC

region_overlap_cad_snp_count_non_OT_gene <- data.frame(n_region = c(length(unique(scE2G_only_region_genes$region_id[scE2G_only_region_genes$OT == 0])),length(unique(scE2G_only_no_promoter_region_genes$region_id[scE2G_only_no_promoter_region_genes$OT == 0])),length(unique(ABC_only_region_genes$region_id[ABC_only_region_genes$OT == 0])),length(unique((ABC_scE2G_region_genes$region_id[ABC_scE2G_region_genes$OT == 0])))),
region = c("scE2G-only","scE2G-only-no-promoter","ABC-only","scE2G-ABC"))

#   n_region                 region
# 1       30             scE2G-only
# 2       27 scE2G-only-no-promoter
# 3       38               ABC-only
# 4        6              scE2G-ABC

#number of linked genes (OT or non OT)
linked_genes_count <- data.frame(n_genes = c(sum(all_genes$scE2G_ABC_region[is.na(all_genes$globalScore)]),sum(all_genes$scE2G_only_region[is.na(all_genes$globalScore)]),sum(all_genes$scE2G_only_no_promoter_region[is.na(all_genes$globalScore)]),sum(all_genes$ABC_only_region[is.na(all_genes$globalScore)]),sum(all_genes$scE2G_ABC_region[!is.na(all_genes$globalScore)]),sum(all_genes$scE2G_only_region[!is.na(all_genes$globalScore)]),sum(all_genes$scE2G_only_no_promoter_region[!is.na(all_genes$globalScore)]),sum(all_genes$ABC_only_region[!is.na(all_genes$globalScore)])),
category=c(rep("non_OT",4),rep("OT",4)), region = rep(c("scE2G_ABC_region","scE2G_only_region","scE2G_only_no_promoter_region","ABC_only_region"),2))

linked_genes_count$category <- factor(linked_genes_count$category, levels = rev(c("non_OT", "OT")))

svg("figures/compare_scE2G_ABC/OT_genes_count.svg")
ggplot(linked_genes_count  , aes(x = region, y = n_genes, fill = category)) +
  geom_bar(stat = "identity") +
#  geom_hline(yintercept = nrow(open_target), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Number of genes linked per dataset",
    x = "Dataset",
    y = "Count",
    fill = "OpenTarget genes"
  ) +
  theme_minimal()
dev.off()


#ALL OT genes only
all_genes <- data.frame(gene = open_target$symbol)
all_genes$scE2G_ABC_region <- ifelse(all_genes$gene %in% ABC_scE2G_region_genes$gene,1,0)
all_genes$scE2G_only_no_promoter_region <- ifelse(all_genes$gene %in% scE2G_only_no_promoter_region_genes$gene,1,0)
all_genes$scE2G_only_region <- ifelse(all_genes$gene %in% scE2G_only_region_genes$gene,1,0)
all_genes$ABC_only_region <- ifelse(all_genes$gene %in% ABC_only_region_genes$gene,1,0)
all_genes <- merge(all_genes,open_target[,c("symbol","globalScore")],by.x="gene",by.y="symbol",all.x=TRUE) 
all_genes$all <-1

svg("figures/compare_scE2G_ABC/venn_diagram_opentarget_cad_genes.svg", width = 6, height = 6)
a <- venn.diagram(
  x = list( all_genes$gene[all_genes$scE2G_ABC_region == 1 & !is.na(all_genes$globalScore)],all_genes$gene[all_genes$scE2G_only_region == 1  & !is.na(all_genes$globalScore)],all_genes$gene[all_genes$ABC_only_region == 1  & !is.na(all_genes$globalScore)] ), # The two lists
  category.names = c("scE2G_ABC CAD genes", "scE2G_only CAD genes","ABC_only CAD genes"), # Labels for the lists
  filename = NULL, imagetype = "svg", output = FALSE,
  fill = c("lightgreen", "skyblue","orange"), # Colors for the circles
  alpha = 0.5, # Transparency of the circles
  cex = 1.5, # Size of the text
  fontface = "bold", # Text style
  cat.cex = 1.5, # Size of the category text
  cat.fontface = "bold", # Text style for category labels
  cat.pos = c(0, 120, 240) # Position of the category labels
)
 grid.draw(a)
 dev.off()



all_genes_long <- all_genes %>%
  pivot_longer(
    cols = c(scE2G_ABC_region, scE2G_only_region,scE2G_only_no_promoter_region,ABC_only_region,all),
    names_to = "region_type",
    values_to = "value"
  )
all_genes_long <- all_genes_long[!all_genes_long$value == 0,]

kruskal.test(globalScore ~ region_type, data = all_genes_long[!is.na(all_genes_long$globalScore),])  # non parametric anova, p-value = 0.005
wilcox_res <- as.data.frame(pairwise.wilcox.test(all_genes_long$globalScore,
                     all_genes_long$region_type,
                     p.adjust.method = "BH")$p.value)
wilcox_res[] <- lapply(wilcox_res, function(x) {
      if (is.numeric(x)) format(x, scientific = TRUE, digits = 3) else x
        })

wilcox_res

svg("figures/compare_scE2G_ABC/OT_genes_score_with_fullOT_category.svg")
ggplot(all_genes_long[!is.na(all_genes_long$globalScore),], aes(x=region_type, y = globalScore)) +
  geom_boxplot(fill = "steelblue") +
  theme_minimal() +
  labs(title = "OT gene score distrib", x = "Category", y = "Frequency")
dev.off()