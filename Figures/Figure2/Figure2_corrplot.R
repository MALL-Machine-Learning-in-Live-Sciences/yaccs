# Corrplot literature signatures
# =====
source("requirements.R")

# Inputpaths
signatures_inputpath <- "extdata/signatures/signatures_coad.rds"
counts_inputpath <- "data/pp/counts_norm_coad_patients.rds"

# Load data
coad_signatures <- readRDS(signatures_inputpath)
counts <- readRDS(counts_inputpath) %>% t()

# Arguments
method <- "plage"

# Convert to
signatures_coad_symbol <- list()
for (i in seq_along(coad_signatures)) {
  signatures_coad_symbol[[i]] <- AnnotationDbi::mapIds(
                                                    org.Hs.eg.db,
                                                    keys = coad_signatures[[i]],
                                                    column = "SYMBOL",
                                                    keytype = "ENSEMBL",
                                                    multiVals = "first"
                                                    )
}
names(signatures_coad_symbol) <- names(coad_signatures)

for (s in seq_along(signatures_coad_symbol)) {
    signatures_coad_symbol[[s]] <- names(signatures_coad_symbol[[s]])
}

# Run GSVA
gsva <- gsva(
    counts,
    signatures_coad_symbol,
    method = method,
    kcdf = "Gaussian"
    )
gsva <- as.data.frame(t(gsva))
names(gsva)[grep("ours", names(gsva))] <- "YACCS"

# order by SigCheck significance
s <- readRDS("data/sigcheck/scKnownLiterature_ours.rds")
sort <- names(sort(c(YACCS = s$survivalPval, s$survivalPvalsKnown)))

gsva <- gsva[, match(sort, colnames(gsva))]

# Plotting
corrplot(
    cor(gsva),
    method = "circle",
    type = "lower",
    col = viridis(10),
    tl.col = "black"
    )
