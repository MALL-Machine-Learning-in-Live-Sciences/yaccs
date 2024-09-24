# Feature selection with FCBF
source("R/utils.R")
source("requirements.R")

# Inputpaths
metadata_inputpath <- "data/data_partitions/metadata_train.rds"
counts_inputpath <- "data/data_partitions/counts_train.rds"

# Outputpaths
outputpath <- "data/yaccs/yaccs_annot.rds"

# Arguments
min_su <- 0.0025

# Load data
counts_train <- readRDS(counts_inputpath)
metadata_train <- readRDS(metadata_inputpath)

# Run FCBF
data_fcbf <- fast.cor.FS(
    x = counts_train,
    y = metadata_train$vital_status,
    min_su = min_su
)

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
yaccs <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = names(data_fcbf),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL",
  multiVals = "first"
)

yaccs <- yaccs[-which(yaccs$SYMBOL == "LGALS7B"), ]
rownames(yaccs) <- 1:nrow(yaccs)

# Store signature
saveRDS(yaccs, file = file.path(outputpath))
