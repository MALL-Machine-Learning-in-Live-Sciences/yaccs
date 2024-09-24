source("R/utils.R")
source("requirements.R")

# Inputpaths
metadata_inputpath <- "data/pp/metadata_coad_patients.rds"
counts_inputpath <- "data/pp/counts_norm_coad_patients.rds"

# Outputpaths
outputdir <- "data/data_partitions/"

# Load data
metadata <- readRDS(metadata_inputpath)
counts <- readRDS(counts_inputpath)

# Define partitions
smpSize <- floor(0.9 * nrow(metadata))
set.seed(55)
trainIdx <- sample(seq_len(nrow(metadata)), size = smpSize)

# Split into train and test
metadata_train <- metadata[trainIdx, ]
metadata_test <- metadata[-trainIdx, ]

counts_train <- counts[trainIdx, ]
counts_test <- counts[-trainIdx, ]

# Save partitions
saveRDS(metadata_train, file.path(outputdir, "metadata_train.rds"))
saveRDS(metadata_test, file.path(outputdir, "metadata_test.rds"))
saveRDS(counts_train, file.path(outputdir, "counts_train.rds"))
saveRDS(counts_test, file.path(outputdir, "counts_test.rds"))
