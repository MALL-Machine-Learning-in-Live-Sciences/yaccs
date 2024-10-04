# Network plot from HTRAnalyzeR2 results
source("requirements.R")

# Inputs
inputpath <- "data/pathway-enrichment/tcga-plage-reactome_high_plage.rds"

# Load data
data <- readRDS(inputpath)

data@para$pValueCutoff <- 0.05
data@result$GSEA.results$Reactome$Gene.Set.Term <- tolower(gsub('REACTOME_', '', data@result$GSEA.results$Reactome$Gene.Set.Term))
data@result$GSEA.results$Reactome$Gene.Set.Term = gsub('_', ' ', data@result$GSEA.results$Reactome$Gene.Set.Term)

# Plot network
# ===
viewEnrichMap(
    data,
    resultName = "GSEA.results",
    gscs = "Reactome",
    allSig = TRUE,
    gsNameType = "term")

# Plot GSEA
# ====
topGS <- getTopGeneSets(
    data,
    resultName = "GSEA.results",
    gscs = "Reactome",
    allSig = TRUE)

viewGSEA(data, gscName = 'Reactome', gsName = topGS[["Reactome"]][grep("ESTABLISHMENT_OF_SISTER_CHROMATID_COHESION", topGS[["Reactome"]])])
viewGSEA(data, gscName = 'Reactome', gsName = topGS[["Reactome"]][grep("MITOTIC_PROMETAPHASE", topGS[["Reactome"]])])
viewGSEA(data, gscName = 'Reactome', gsName = topGS[["Reactome"]][grep("M_PHASE", topGS[["Reactome"]])])
viewGSEA(data, gscName = 'Reactome', gsName = topGS[["Reactome"]][grep("CELL_CYCLE_MITOTIC", topGS[["Reactome"]])])

