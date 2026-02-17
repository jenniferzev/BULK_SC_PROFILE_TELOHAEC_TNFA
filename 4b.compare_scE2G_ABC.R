library(ggplot2)


reformat_result <- function(h2_results){
    h2_results$Enrichment_fdr <- p.adjust(h2_results$Enrichment_p, method = "fdr") 
    h2_results$Coefficient_p <-  2 * (1 - pnorm(abs(h2_results$Coefficient_z.score)))
    h2_results$Coefficient_fdr <- p.adjust(h2_results$Coefficient_p, method = "fdr")
    if (any(h2_results$Enrichment_fdr <= 0.05, na.rm = TRUE)) {
        max_enrich_p <- max(h2_results$Enrichment_p[h2_results$Enrichment_fdr <= 0.05], na.rm = TRUE)
        } else {
        max_enrich_p <- NA
        }
    if (any(h2_results$Coefficient_fdr <= 0.05, na.rm = TRUE)) {
        max_coef_p <- max(h2_results$Coefficient_p[h2_results$Coefficient_fdr <= 0.05], na.rm = TRUE)
        } else {
        max_coef_p <- NA
        }
    return(list(df=h2_results,enrichment_fdr_thr=max_enrich_p,coefficient_fdr_thr=max_coef_p))
}

h2_enrichement_figure <- function(h2_results,enrichment_threshold,label){

    p <- ggplot(h2_results , aes(x = Enrichment, y = -log10(Enrichment_p))) +
    geom_point(color='black') +                             # Scatter plot
    theme_minimal() +
    labs(title = label, x = "Enrichment", y = "-log10(Enrichment_p)")+
    geom_text(aes(label = Phenotype), vjust = -1, size = 4)+
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "red")+
    geom_hline(yintercept = -log10(enrichment_threshold), color = "blue", linetype = "dotted") +
    theme(plot.background = element_rect(fill = "white", color = NA))+
    coord_cartesian(
    xlim = c(5, 12),
    ylim = c(0, 4)
    )
    #coord_cartesian(clip = "off") 
    ggsave(paste("figures/b.compare_scE2G_ABC/",label,"_enrich_plot.svg",sep=""), p, width = 16, height = 8)

}

dbp_h2 <- read.table("data/b.compare_scE2G_ABC/DBP_h2.txt",header=TRUE)
dbp_h2 <- reformat_result(dbp_h2)
h2_enrichement_figure(dbp_h2$df,dbp_h2$enrichment_fdr_thr,"dbp_h2")

cad_h2 <- read.table("data/b.compare_scE2G_ABC/CAD_h2.txt",header=TRUE)
cad_h2 <- reformat_result(cad_h2)
h2_enrichement_figure(cad_h2$df,cad_h2$enrichment_fdr_thr,"cad_h2")