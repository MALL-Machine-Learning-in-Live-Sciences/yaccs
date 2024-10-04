# Plot heatmap of expression ordered by plage score and divided by cell cycle
source("requirements.R")

# Input paths
tcga_inputpath <- "data/pp/counts_norm_coad_patients.rds"
tcga_annot <- "data/pathway-enrichment/tcga_yaccs_annot.rds"
cc_signatures_inputpath <- "extdata/signatures/cell_cycle.rds"
coad_signatures_inputpath <- "extdata/signatures/signatures_coad.rds"

# Load data
tcga <- readRDS(tcga_inputpath)
tcga_annot <- readRDS(tcga_annot)
cc_signatures <- readRDS(cc_signatures_inputpath)

# Plot heatmap
meta <-
    tcga_annot %>%
    arrange(plageScore)

# Annotations
ha_col <- HeatmapAnnotation(
    df = meta
)

ha <- list()
for (i in seq_along(cc_signatures)) {

    common_genes <- intersect(cc_signatures[[i]], colnames(tcga))

    ha_row <- rowAnnotation(
        Signature = rep(names(cc_signatures)[i], length(common_genes)))

    # Plot
    mat <- scale(tcga)
    toplot <- t(mat)
    ha[[i]] <- Heatmap(
        toplot[common_genes, rownames(meta)],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        row_title = "Signatures",
        column_title = "Samples",
        top_annotation = ha_col,
        right_annotation = ha_row
        )


    print(names(cc_signatures)[i])

}