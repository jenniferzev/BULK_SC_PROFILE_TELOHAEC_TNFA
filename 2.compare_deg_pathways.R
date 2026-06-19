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


pathways_res <- read.table("data/compare_deg_pathways/sc_telohaec_pathway.txt",header=TRUE,sep="\t")
filtered_pathway_res <- pathways_res[pathways_res$PValue < 0.05, c("Term","Count","PValue","Genes")]
filtered_pathway_res$ID <- gsub(":.*","",filtered_pathway_res$Term)

lalonde_pathways <- read.table("data/compare_deg_pathways/new_lalonde_pathways.txt",header=TRUE,sep="\t")
lalonde_pathways<- lalonde_pathways[lalonde_pathways$PValue < 0.05, c("Term","Count","PValue","Genes")]
lalonde_pathways$ID <- gsub(":.*","",lalonde_pathways$Term)

svg("figures/compare_deg_pathway/venn_diagram_pathway.svg", width = 6, height = 6)
b <- venn.diagram(
  x = list(filtered_pathway_res$Term,lalonde_pathways$Term), # The two lists
  category.names = c("sc DEGs enriched KEGG Pathways", "bulk DEGs enriched KEGG Pathways"), # Labels for the lists
  filename = NULL, imagetype = "svg", # You can specify a file to save the image (e.g., "venn.png")
  fill = c("lightgreen", "skyblue"), # Colors for the circles
  alpha = 0.5, # Transparency of the circles
  cex = 1.5, # Size of the text
  fontface = "bold", # Text style
  cat.cex = 1.5, # Size of the category text
  cat.fontface = "bold", # Text style for category labels
  cat.pos = c(0, 180) # Position of the category labels
)
grid.draw(b)
dev.off()

#I keep the old bulk pathway analysis for now
pathway_bulk_sc <- intersect(filtered_pathway_res$Term,lalonde_pathways$Term)
pathway_sc_only <- setdiff(filtered_pathway_res$Term,lalonde_pathways$Term)
pathway_bulk_only <- setdiff(lalonde_pathways$Term,filtered_pathway_res$Term)


#compare pvalue
filtered_pathway_res$experience <- "single cell"
lalonde_pathways$experience <- "bulk"

pathway_data <- rbind(filtered_pathway_res,lalonde_pathways)
pathway_data$log_pval <- -(log10(pathway_data$PValue))
pathway_data <- pathway_data[order(pathway_data$log_pval,decreasing = TRUE),]


ggplot(pathway_data[pathway_data$Term %in% pathway_bulk_sc, ], aes(x = experience, y = Term, size = log_pval, color=log_pval)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Pathway", size = "-Log10 P", color = "-Log10 P") +
    theme_bw()
ggsave("figures/compare_deg_pathway/pathways_sc_bulk.svg")

ggplot(pathway_data[pathway_data$Term %in% pathway_sc_only , ], aes(x = experience, y = Term, size = log_pval, color=log_pval)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Pathway", size = "-Log10 P", color = "-Log10 P") +
    theme_bw()
ggsave("figures/compare_deg_pathway/pathways_sc_only.svg")

ggplot(pathway_data[pathway_data$Term %in% pathway_bulk_only, ], aes(x = experience, y = Term, size = log_pval, color=log_pval)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Pathway", size = "-Log10 P", color = "-Log10 P") +
    theme_bw()
ggsave("figures/compare_deg_pathway/pathways_bulk_only.svg")


sc_only_deg_pathways_res <- read.table("/home/jzev/projects/def-glettre/jzev/BULK_SC_PROFILE_TELOHAEC_TNFA/data/compare_deg_pathways/sc_only_deg_telohaec_pathway.txt",header=TRUE,sep="\t")
sc_only_deg_pathways_res <- sc_only_deg_pathways_res[sc_only_deg_pathways_res$PValue < 0.05, c("Term","Count","PValue","Genes","Category")]
sc_only_deg_pathways_res$ID <- gsub(":.*","",sc_only_deg_pathways_res$Term)
sc_only_deg_pathways_res$log_pval <- -(log10(sc_only_deg_pathways_res$PValue))

ggplot(sc_only_deg_pathways_res, aes(x = Category, y = Term, size = log_pval, color=log_pval)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Pathway", size = "-Log10 P", color = "-Log10 P") +
    theme_bw()
ggsave("figures/compare_deg_pathway/pathways_sc_deg_only.svg")

bulk_only_deg_pathways_res <- read.table("/home/jzev/projects/def-glettre/jzev/BULK_SC_PROFILE_TELOHAEC_TNFA/data/compare_deg_pathways/bulk_only_deg_telohaec_pathway.txt",header=TRUE,sep="\t")
bulk_only_deg_pathways_res <- bulk_only_deg_pathways_res[bulk_only_deg_pathways_res$PValue < 0.05, c("Term","Count","PValue","Genes","Category")]
bulk_only_deg_pathways_res$ID <- gsub(":.*","",bulk_only_deg_pathways_res$Term)
bulk_only_deg_pathways_res$log_pval <- -(log10(bulk_only_deg_pathways_res$PValue))
bulk_only_deg_pathways_res <- bulk_only_deg_pathways_res[order(bulk_only_deg_pathways_res$log_pval,decreasing = TRUE),]

ggplot(bulk_only_deg_pathways_res, aes(x = Category, y = Term, size = log_pval, color=log_pval)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", y = "Pathway", size = "-Log10 P", color = "-Log10 P") +
    theme_bw()
ggsave("figures/compare_deg_pathway/pathways_bulk_deg_only.svg")

