# Plot correlation plot of coad signatures and cell cycle signatures
source("requirements.R")

# Input paths
cc_signatures_inputpath <- "extdata/signatures/cell_cycle.rds"
tcga_inputpath <- "data/pp/counts_norm_coad_patients.rds"
coad_signatures_inputpath <- "extdata/signatures/signatures_coad.rds"

# Load data
train <- readRDS(tcga_inputpath)
cc_signatures <- readRDS(cc_signatures_inputpath)
coad_signatures <- readRDS(coad_signatures_inputpath)

# Define signatures
signatures <- c(cc_signatures, coad_signatures)

# Calculate signature´s activity
gsva <- gsva(
    t(train),
    signatures,
    method = "plage",
    kcdf = "Gaussian")
gsva <- as.data.frame(t(gsva))

signatures_ordered <- c(
    "G0.EarlyG1",
    "G1Phase",
    "G1S.transition",
    "S.phase",
    "G2.phase",
    "G2M.transition",
    "M.prophase",
    "M.prometaphase",
    "M.metaphase.anaphase",
    "M.telophase.cytokinesis",
    "ours",
    "mda114",
    "zhang",
    "coloGuideEx",
    "coloPrint",
    "coloGuidePro"
    )

gsva <- gsva[, match(signatures_ordered, names(gsva))]

# Plotting
corrplot::corrplot(
    cor(gsva),
    method = "number",
    type = "lower",
    col = viridis(16),
    tl.col = "black")
