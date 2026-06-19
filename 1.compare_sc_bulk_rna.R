# COMPARE BULK AND SINGLE CELL DATA FOR TELOHAEC +/- TNFA

## PACKAGES
library(dplyr)
library(Seurat) #for RNAseq analysis
library(patchwork)
library(ggplot2)
library(glmGamPoi)
library("openxlsx")
library(ggrepel)
library(VennDiagram)
library(tidyr)
library(svglite)
library(ggridges)
library(ggrastr)


################################## COMPARE DEG 

#SC DIFFERENTIAL EXPRESSION ANALYSIS
load(file="data/compare_sc_bulk_rna/sc_telohaec_tnfa_deg_results_gene_symbols.Rdata") #rna_deg_results
filtered_rna_deg_results <- rna_deg_results[rna_deg_results$is_deg == 1,]


#BULK DIFFERENTIAL EXPRESSION ANALYSIS
lalonde_rna_deseq <- read.table("data/compare_sc_bulk_rna/bulk_deg_deseq2_results.tsv",header=TRUE,sep="\t")
lalonde_rna_deseq$FDR <- as.numeric(lalonde_rna_deseq$FDR) 
lalonde_rna_deseq$deg_NT_4h <- ifelse(abs(lalonde_rna_deseq$LFC_NT_4h) > 0.3, 1, 0)
lalonde_rna_deseq$deg_NT_24h <- ifelse(abs(lalonde_rna_deseq$LFC_NT_24h) > 0.3, 1, 0)
lalonde_rna_deseq$deg_4h_24h <- ifelse(abs(lalonde_rna_deseq$LFC_4h_24h) > 0.3, 1, 0)
filtered_lalonde_rna_deseq <- lalonde_rna_deseq[(abs(lalonde_rna_deseq$LFC_NT_4h) > 0.3 | abs(lalonde_rna_deseq$LFC_NT_24h) > 0.3 | abs(lalonde_rna_deseq$LFC_4h_24h) > 0.3) & lalonde_rna_deseq$FDR < 0.001,]

#GENES THAT WERE PRESENTS IN BOTH DATASETS
intersect_genes <- intersect(rna_deg_results$gene, lalonde_rna_deseq$gene_name)
length(intersect_genes) # 11221 genes are in both analysis


## DEG VENN DIAGRAM
svg("figures/compare_sc_bulk_rna/compare_deg.svg", width = 6, height = 6)
a <- venn.diagram(
  x = list( filtered_rna_deg_results$gene[filtered_rna_deg_results$gene %in% intersect_genes ],filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$gene_name %in% intersect_genes ]), # The two lists
  category.names = c("DEGs from sc analysis", "DEGs from bulk analysis"), # Labels for the lists
  filename = NULL, imagetype = "svg", output = FALSE,
  fill = c("lightgreen", "skyblue"), # Colors for the circles
  alpha = 0.5, # Transparency of the circles
  cex = 1.5, # Size of the text
  fontface = "bold", # Text style
  cat.cex = 1.5, # Size of the category text
  cat.fontface = "bold", # Text style for category labels
  cat.pos = c(0, 170) # Position of the category labels
)
 grid.draw(a)
 dev.off()


## DIFFERENTIAL EXPRESSION ANALYSIS LOG10FC COMPARISON
deg_compare_plot <- function(time){
  if(time == "4hr_0hr"){
    rna_deg_results_col <- c("log10FC_4hr_0hr","deg_4hr_0hr")
    lalonde_rna_deseq_col <- c("LFC_NT_4h","deg_NT_4h")
  }else if(time == "24hr_0hr"){
    rna_deg_results_col <- c("log10FC_24hr_0hr","deg_24hr_0hr")
    lalonde_rna_deseq_col <- c("LFC_NT_24h","deg_NT_24h")
  }else if(time == "24hr_4hr" ){
    rna_deg_results_col <- c("log10FC_24hr_4hr","deg_24hr_4hr")
    lalonde_rna_deseq_col <- c("LFC_4h_24h","deg_4h_24h")
  }

  deg_plot_data <- merge(rna_deg_results[,c("gene",rna_deg_results_col)], lalonde_rna_deseq[,c("gene_name",lalonde_rna_deseq_col)], by.x="gene", by.y="gene_name") #only compare genes that are in both datasets
  colnames(deg_plot_data) <- c("gene","LFC_scRNA_seq","sc_deg","LFC_bulk_RNAseq","bulk_deg")
  deg_plot_data$deg_category <- ifelse(deg_plot_data$sc_deg == 1 & deg_plot_data$bulk_deg == 1,"both",
                                    ifelse(deg_plot_data$sc_deg == 1, "sc",
                                    ifelse(deg_plot_data$bulk_deg == 1, "bulk","none")))
  deg_plot_data$deg_category <- factor( deg_plot_data$deg_category, levels = c("none", "bulk", "sc","both"))
                                  

  model <- lm(LFC_scRNA_seq ~ LFC_bulk_RNAseq, data = deg_plot_data)
  slope <- coef(model)[2]

  ggplot(deg_plot_data[order(deg_plot_data$deg_category),], aes(x = LFC_scRNA_seq, y = LFC_bulk_RNAseq,color=deg_category)) +
  geom_point_rast(
    data = subset(deg_plot_data, deg_category == "none"),
    aes(LFC_scRNA_seq, LFC_bulk_RNAseq),
    color = "gray",
    alpha = 1,
    size = 0.5
  ) +
  geom_point_rast(
    data = subset(deg_plot_data, deg_category != "none"),
    aes(LFC_scRNA_seq, LFC_bulk_RNAseq,
        color = deg_category),
    alpha = 0.3,
    size = 0.8
  )+
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color="red")
  labs(
       x = "scRNAseq LFC",
       y = "bulk RNAseq LFC",
       ) +
  theme_classic()
  ggsave(paste("figures/compare_sc_bulk_rna/lfc_compare_plot_",time,".svg",sep=""),device = "svg", width = 6, height = 4, dpi = 300)
  deg_plot_data <- deg_plot_data[complete.cases(deg_plot_data[, c("LFC_scRNA_seq", "LFC_bulk_RNAseq")]), ]
  spearman_cor <- cor.test(deg_plot_data[["LFC_scRNA_seq"]], deg_plot_data[["LFC_bulk_RNAseq"]], method = "spearman",use = "pairwise.complete.obs")
  pearson_cor<- cor.test(deg_plot_data[["LFC_scRNA_seq"]], deg_plot_data[["LFC_bulk_RNAseq"]], method = "pearson",use = "pairwise.complete.obs")
 correlation_df <- data.frame(
  spearman_cor = as.numeric(spearman_cor$estimate),
  pearson_cor = as.numeric(pearson_cor$estimate),
  spearman_pval = format.pval(spearman_cor$p.value, digits = 3, eps = .Machine$double.eps),
  pearson_pval = format.pval(pearson_cor$p.value, digits = 3, eps = .Machine$double.eps)
)
   colnames(correlation_df) <- paste(
    c("spearman_cor", " pearson_cor", "spearman_pval", "pearson_pval"),
    time, sep = "_"
  )
  return(correlation_df)
}

cor_4hr_0hr <- deg_compare_plot("4hr_0hr") 
cor_24hr_0hr <- deg_compare_plot("24hr_0hr") 
cor_24hr_4hr <- deg_compare_plot("24hr_4hr")



## COMPARE SC-ONLY,BULK-ONLY,SC-BULK DEG -- TNF 4hr vs TNF 0hr

sc_bulk_deg_NT_4h <- intersect(filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_4hr_0hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes],
                              filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_NT_4h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes])
sc_only_deg_NT_4h <- setdiff(filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_4hr_0hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes],
                           filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_NT_4h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes])
bulk_only_deg_NT_4h <- setdiff(filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_NT_4h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes],
                             filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_4hr_0hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes])

deg_compare_NT_4h <- merge(rna_deg_results[,c("gene","log10FC_4hr_0hr","p_val_adj_4hr_0hr","deg_4hr_0hr")],
                          lalonde_rna_deseq[,c("gene_name","LFC_NT_4h","FDR","deg_NT_4h")],
                          by.x="gene",by.y="gene_name") 

deg_compare_NT_4h$deg <- ifelse(deg_compare_NT_4h$gene %in% sc_bulk_deg_NT_4h , "sc_bulk_NT_4h",
                            ifelse(deg_compare_NT_4h$gene %in% sc_only_deg_NT_4h, "sc_only_NT_4h",
                                  ifelse(deg_compare_NT_4h$gene %in% bulk_only_deg_NT_4h, "bulk_only_NT_4h","none")))
                                  
ggplot(deg_compare_NT_4h[!(deg_compare_NT_4h$deg %in% c("none")),], aes(x = abs(log10FC_4hr_0hr), y =deg,color=deg)) +
geom_density_ridges() +
labs(
      x = "abs(scRNAseq LFC TNFa 4h vs 0h)"
      ) +
theme_classic()
ggsave(paste("figures/compare_sc_bulk_rna/lfc_NT_4hr_DEG_sc_ridgeplot.svg",sep=""),device = "svg", width = 10, height = 4, dpi = 300)


## COMPARE SC-ONLY,BULK-ONLY,SC-BULK DEG -- TNF 24hr vs TNF 0hr


sc_bulk_deg_NT_24h <- intersect(filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_24hr_0hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes],
                              filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_NT_24h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes])
sc_only_deg_NT_24h <- setdiff(filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_24hr_0hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes],
                           filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_NT_24h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes])
bulk_only_deg_NT_24h <- setdiff(filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_NT_24h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes],
                             filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_24hr_0hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes])

deg_compare_NT_24h <- merge(rna_deg_results[,c("gene","log10FC_24hr_0hr","p_val_adj_24hr_0hr","deg_24hr_0hr")],
                          lalonde_rna_deseq[,c("gene_name","LFC_NT_24h","FDR","deg_NT_24h")],
                          by.x="gene",by.y="gene_name") 

deg_compare_NT_24h$deg <- ifelse(deg_compare_NT_24h$gene %in% sc_bulk_deg_NT_24h , "sc_bulk_NT_24h",
                             ifelse(deg_compare_NT_24h$gene %in% sc_only_deg_NT_24h, "sc_only_NT_24h",
                                    ifelse(deg_compare_NT_24h$gene %in% bulk_only_deg_NT_24h, "bulk_only_NT_24h","none")))

ggplot(deg_compare_NT_24h[!(deg_compare_NT_24h$deg %in% c("none")),], aes(x = abs(log10FC_24hr_0hr), y =deg,color=deg)) +
geom_density_ridges() +
labs(
      x = "abs(scRNAseq LFC TNFa 24h vs 0h)"
      ) +
theme_classic()
ggsave(paste("figures/compare_sc_bulk_rna/lfc_NT_24hr_DEG_sc_ridgeplot.svg",sep=""),device = "svg", width = 10, height = 4, dpi = 300)


## COMPARE SC-ONLY,BULK-ONLY,SC-BULK DEG -- TNF 24hr vs TNF 4hr

sc_bulk_deg_4h_24h <- intersect(filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_24hr_4hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes],
                            filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_4h_24h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes])
sc_only_deg_4h_24h <- setdiff(filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_24hr_4hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes],
                          filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_4h_24h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes])
bulk_only_deg_4h_24h <- setdiff(filtered_lalonde_rna_deseq$gene_name[filtered_lalonde_rna_deseq$deg_4h_24h == 1 & filtered_lalonde_rna_deseq$gene_name %in% intersect_genes],
                            filtered_rna_deg_results$gene[filtered_rna_deg_results$deg_24hr_4hr == 1 & filtered_rna_deg_results$gene %in% intersect_genes])

deg_compare_4h_24h <- merge(rna_deg_results[,c("gene","log10FC_24hr_4hr","p_val_adj_24hr_4hr","deg_24hr_4hr")],
                        lalonde_rna_deseq[,c("gene_name","LFC_4h_24h","FDR","deg_4h_24h")],
                        by.x="gene",by.y="gene_name") 

deg_compare_4h_24h$deg <- ifelse(deg_compare_4h_24h$gene %in% sc_bulk_deg_4h_24h , "sc_bulk_4h_24h",
                            ifelse(deg_compare_4h_24h$gene %in% sc_only_deg_4h_24h, "sc_only_4h_24h",
                                  ifelse(deg_compare_4h_24h$gene %in% bulk_only_deg_4h_24h, "bulk_only_4h_24h","none")))

                                      ggplot(deg_compare_4h_24h[!(deg_compare_4h_24h$deg %in% c("none")),], aes(x = abs(log10FC_24hr_4hr), y =deg,color=deg)) +
geom_density_ridges() +
labs(
      x = "abs(scRNAseq LFC TNFa 24h vs 4h)"
      ) +
theme_classic()
ggsave(paste("figures/compare_sc_bulk_rna/lfc_4h_24hr_DEG_sc_ridgeplot.svg",sep=""),device = "svg", width = 10, height = 4, dpi = 300)
