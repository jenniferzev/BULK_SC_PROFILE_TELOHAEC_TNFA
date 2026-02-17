# COMPARE BULK AND SINGLE CELL DATA FOR TELOHAEC +/- TNFA

## PACKAGES
library(dplyr)
library(Seurat) #for RNAseq analysis
library(patchwork)
library(ggplot2)
library(presto) ####
library(Signac) #fot ATACseq analysis
library(clustree)
library(RColorBrewer)
library(DESeq2)
library(scDblFinder)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationHub)
library(biovizBase)
library(glmGamPoi)
library("openxlsx")
library(ggrepel)
library(VennDiagram)
library(tidyr)
library(dplyr)
library(UpSetR)
library(ggrastr)
library(ComplexUpset)
library(ggridges)


############################ PART 1
###### COMPARE ATAC PEAK LENGTHS BETWEEN SC AND BULK

# SC ATAC PEAKS
load(file="data/compare_sc_bulk_atac/sc_telohaec_tnfa_dop_results_26nov25.Rdata")
atac_dop_results$peak_id <- atac_dop_results$peak
atac_dop_results <- separate(atac_dop_results, col = peak, into = c("chr", "start","stop"), sep = "-")
atac_dop_results$stop <- as.numeric(atac_dop_results$stop)
atac_dop_results$start <- as.numeric(atac_dop_results$start)


atac_dop_results$peak_length <- atac_dop_results$stop - atac_dop_results$start
svg("figures/compare_sc_bulk_atac/sc_telohaec_tnfa_atac_peaks_length.svg")
hist(atac_dop_results$peak_length,breaks=50, main=NULL, xlab = "sc Peak length (bp)", ylab = "Number of peaks")
dev.off()

# BULK ATAC PEAKS
lalonde_atac_deseq <- read.table(file="data/compare_sc_bulk_atac/atacpeaks_DE_LFC.txt",header=TRUE)
lalonde_atac_deseq$peak_id <- paste(lalonde_atac_deseq$chr,lalonde_atac_deseq$start,lalonde_atac_deseq$stop, sep="_")
lalonde_atac_peak_hg38 <- read.table(file="data/compare_sc_bulk_atac/atacpeaks_DE_hg38.bed")  #liftOver results
colnames(lalonde_atac_peak_hg38) <- c("chr_hg38","start_hg38","stop_hg38","peak_id_hg19")
lalonde_atac_deseq <- merge(lalonde_atac_deseq,lalonde_atac_peak_hg38,by.x="peak_id",by.y="peak_id_hg19",all.x=TRUE) #merge the coordinates


lalonde_atac_deseq$peak_length <- lalonde_atac_deseq$stop_hg38 - lalonde_atac_deseq$start_hg38
svg("figures/compare_sc_bulk_atac/lalonde_atac_peaks_length.svg")
hist(lalonde_atac_deseq$peak_length, main=NULL, xlab = "bulk Peak length (bp)", ylab = "Number of peaks")
dev.off()


###### EXPLORE OVERLAP BETWEEN BULK AND SC ATAC PEAKS

# BED FILE WITH SC ATAC PEAK COORDINATES(hg18)
write.table(x=atac_dop_results[,c("chr","start","stop","peak_id")], file="data/compare_sc_bulk_atac/sc_telohaec_tnf_atac_peaks.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# SEE SECTION "OVERLAP BETWEEN SC AND BULK ATAC PEAK" in 3.compare_sc_bulk_atac.sh

atac_peak_overlap_amount <- read.table(file="data/compare_sc_bulk_atac/atac_peak_overlap_amount.bed")
colnames(atac_peak_overlap_amount) <- c("bulk_chr_hg38","bulk_start_hg38","bulk_stop_hg38","bulk_id_hg19","sc_chr_hg38","sc_start_hg38","sc_stop_hg38","sc_id_hg38","overlap")

svg("figures/compare_sc_bulk_atac/atac_peak_overlap_distrib.svg")
hist(atac_peak_overlap_amount$overlap,main = NULL, xlab = "Overlap length (bp)", ylab = "Number of overlapping peaks pairs")
dev.off()

###### DEFINE ATAC PEAK PAIRS

# Filter overlap by overlap size
atac_peak_pairs <- atac_peak_overlap_amount %>% 
  filter(overlap >= 250)
write.table(x=atac_peak_pairs, file="data/compare_sc_bulk_atac/atac_peak_pair.bed", quote=FALSE, row.names=FALSE)


############################ PART 2
###### COMPARE PEAK LFC BETWEEN SC AND BULK CHROMATIN ACCESSIBILITY

atac_peak_pairs <- read.table(file="data/compare_sc_bulk_atac/atac_peak_pair.bed", header=TRUE)

load(file="data/compare_sc_bulk_atac/sc_telohaec_tnfa_dop_results_26nov25.Rdata")

lalonde_atac_deseq <- read.table(file="data/compare_sc_bulk_atac/atacpeaks_DE_LFC.txt",header=TRUE)
lalonde_atac_deseq$peak_id <- paste(lalonde_atac_deseq$chr,lalonde_atac_deseq$start,lalonde_atac_deseq$stop, sep="_")

colnames(atac_dop_results) <- paste("sc_",colnames(atac_dop_results),sep="")
atac_peak_pairs <- merge(atac_peak_pairs,atac_dop_results[,c("sc_log10FC_4hr_0hr","sc_log10FC_24hr_0hr","sc_log10FC_24hr_4hr","sc_peak")], by.x="sc_id_hg38", by.y="sc_peak", all.x=TRUE)
atac_peak_pairs <- merge(atac_peak_pairs,lalonde_atac_deseq[,c("log10FC_wt_4h","log10FC_wt_24h","log10FC_4h_24h","peak_id")], by.x="bulk_id_hg19", by.y="peak_id", all.x=TRUE)

dop_compare_plot <- function(time){
  if(time == "4hr_0hr"){
    select_col <- c("sc_log10FC_4hr_0hr","log10FC_wt_4h")
  }else if(time == "24hr_0hr"){
    select_col <- c("sc_log10FC_24hr_0hr","log10FC_wt_24h")
  }else if(time == "24hr_4hr" ){
    select_col <- c("sc_log10FC_24hr_4hr","log10FC_4h_24h")
  }

  dop_plot_data <-atac_peak_pairs[,select_col]
  colnames(dop_plot_data) <- c("LFC_scATAC_seq","LFC_bulk_ATACseq")

  model <- lm(LFC_scATAC_seq ~ LFC_bulk_ATACseq, data = dop_plot_data)
  slope <- coef(model)[2]

  ggplot(dop_plot_data, aes(x = LFC_scATAC_seq, y = LFC_bulk_ATACseq)) +
  geom_point_rast(size = 2, alpha = 0.7, color="blue") + 
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  labs(
       x = "scATACseq LFC",
       y = "bulk ATACseq LFC",
       ) +

  theme_classic()
  ggsave(paste("figures/compare_sc_bulk_atac/lfc_compare_plot_",time,"_raster.svg",sep=""),device = "svg", width = 6, height = 4, dpi = 300)

  dop_plot_data <- dop_plot_data[complete.cases(dop_plot_data[, c("LFC_scATAC_seq", "LFC_bulk_ATACseq")]), ]
  spearman_cor <- cor.test(dop_plot_data[["LFC_scATAC_seq"]], dop_plot_data[["LFC_bulk_ATACseq"]], method = "spearman",use = "pairwise.complete.obs")
  pearson_cor<- cor.test(dop_plot_data[["LFC_scATAC_seq"]], dop_plot_data[["LFC_bulk_ATACseq"]], method = "pearson",use = "pairwise.complete.obs")
  correlation_df <- data.frame(
  spearman_cor = as.numeric(spearman_cor$estimate),
  pearson_cor = as.numeric(pearson_cor$estimate),
  spearman_pval = format.pval(spearman_cor$p.value, digits = 3, eps = .Machine$double.eps),
  pearson_pval = format.pval(pearson_cor$p.value, digits = 3, eps = .Machine$double.eps)
  )
   colnames(correlation_df) <- paste(
    c("spearman_cor", "pearson_cor", "spearman_pval", "pearson_pval"),
    time, sep = "_"
    )
  return(correlation_df)
}

cor_4hr_0hr <- dop_compare_plot("4hr_0hr") 
cor_24hr_0hr <- dop_compare_plot("24hr_0hr") 
cor_24hr_4hr <- dop_compare_plot("24hr_4hr")


############################ PART 3
################################## COMPARE PEAK DOP STATUS (FOR EACH COMPARISON)

# time="dop_4hr_0hr"
# sc_lfc="log10FC_4hr_0hr"
# lfc="LFC_NT_4h"

#OR

# time="dop_24hr_0hr"
# sc_lfc="log10FC_24hr_0hr"
# lfc="LFC_NT_24h"

#OR
# time="dop_24hr_4hr"
# sc_lfc="log10FC_24hr_4hr"
# lfc="LFC_4h_24h"


#filter single cell atac dop
load(file="data/compare_sc_bulk_atac/sc_telohaec_tnfa_dop_results_26nov25.Rdata")
filtered_dop <- atac_dop_results[atac_dop_results[[time]] == 1, ]

#dop from the article
original_lalonde_dop <- read.xlsx("data/compare_sc_bulk_atac/dop_lalonde_2019.xlsx",sheet= 1,startRow = 3)
original_lalonde_dop$peak_id <- paste(original_lalonde_dop$CHR,original_lalonde_dop$START,original_lalonde_dop$STOP, sep="_")
lalonde_dop <- original_lalonde_dop[original_lalonde_dop[[lfc]] > 0.3, ]

#overlapping pairs info
atac_peak_overlap_amount <- read.table(file="data/compare_sc_bulk_atac/atac_peak_overlap_amount.bed")
colnames(atac_peak_overlap_amount) <- c("bulk_chr_hg38","bulk_start_hg38","bulk_stop_hg38","bulk_id_hg19","sc_chr_hg38","sc_start_hg38","sc_stop_hg38","sc_id_hg38","overlap")
atac_peak_overlap_amount <- atac_peak_overlap_amount[atac_peak_overlap_amount$overlap >= 250, ]
atac_peak_pairs <- atac_peak_overlap_amount 

#define dop in the overlapping pairs
overlapping_dop <- atac_peak_pairs
overlapping_dop$sc_dop <- ifelse(overlapping_dop$sc_id_hg38 %in% filtered_dop$peak, 1, 0)
overlapping_dop <- merge(overlapping_dop,atac_dop_results[,c("peak",sc_lfc)],by.x="sc_id_hg38",by.y="peak", all.x=TRUE)
colnames(overlapping_dop)[which(colnames(overlapping_dop) == sc_lfc)] <- "sc_lfc"
overlapping_dop$bulk_dop <- ifelse(overlapping_dop$bulk_id_hg19 %in% lalonde_dop$peak_id, 1, 0)
overlapping_dop <- merge(overlapping_dop,original_lalonde_dop[,c("peak_id",lfc)],by.x="bulk_id_hg19",by.y="peak_id", all.x=TRUE)
colnames(overlapping_dop)[which(colnames(overlapping_dop) == lfc)] <- "bulk_lfc"

pair_dop_status <- function(sc_dop,bulk_dop){
  if(sc_dop==1 & bulk_dop==1){
    return('sc-bulk DOPs')
  } else if(sc_dop==1 & bulk_dop==0){
    return('sc-only DOPs')
  } else if(sc_dop==0 & bulk_dop==1){
    return('bulk-only DOPs')
  } else {
    return('not DOPs')
  }
}


overlapping_dop$pair_dop_status <- apply(overlapping_dop, 1, function(x) pair_dop_status(x['sc_dop'], x['bulk_dop']))
overlapping_dop$pair_dop_status <- factor(overlapping_dop$pair_dop_status, levels=c('bulk-only DOPs','not DOPs','sc-bulk DOPs','sc-only DOPs'))

svg(paste("figures/compare_sc_bulk_atac/dop_sc_lfc_ridgeplot_",time,".svg",sep=""))
ggplot(overlapping_dop[!overlapping_dop$pair_dop_status %in% c("not DOPs"),], aes(x = abs(sc_lfc), y = pair_dop_status, fill = pair_dop_status)) +
  geom_density_ridges() +
  theme_ridges() + 
  labs(x="log10FC ; scATACseq ",y="DOP annotation") +
  theme_minimal()
dev.off()


upset_plot_input <- data.frame(
  `sc DOP` = c(1,0,1,0),
  `bulk DOP` = c(1,0,0,1),
  `not sc DOP` = c(0,1,0,1),
  `not bulk DOP` = c(0,1,1,0),
  count = c(summary(as.factor(overlapping_dop$pair_dop_status))[3],
            summary(as.factor(overlapping_dop$pair_dop_status))[2],
            summary(as.factor(overlapping_dop$pair_dop_status))[4],
            summary(as.factor(overlapping_dop$pair_dop_status))[1])
)
upset_plot_input_expanded <- upset_plot_input[rep(1:nrow(upset_plot_input), upset_plot_input$count), -5]

svg(paste("figures/compare_sc_bulk_atac/compare_dop_upset",time,".svg"), width = 14, height = 10)
ComplexUpset::upset(
  upset_plot_input_expanded,
  intersect = c("sc.DOP", "bulk.DOP", "not.sc.DOP", "not.bulk.DOP"),
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
)+
  labs(x = "Number of peak pairs")+
  theme(
    strip.text.y = element_text(size = 20, face = "bold"),  # Increases set label text
    axis.text.x = element_blank(),                  # Increases intersection labels
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),              # Increases axis titles
    plot.title = element_text(size = 20, face = "bold")     # If you use ggtitle()
  )
dev.off()