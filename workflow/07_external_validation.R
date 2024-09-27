# External validation
# ----------------
source("requirements.R")

# Inputpaths
coad_signatures_inputpath <- "extdata/signatures/signatures_coad.rds"

# Load data
coad_signatures <- readRDS(coad_signatures_inputpath)

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

# 
names(signatures_coad_symbol)
