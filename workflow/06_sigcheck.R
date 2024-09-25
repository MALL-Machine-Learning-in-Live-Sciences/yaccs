# Runing SigCheck experiment
source("requirements.R")

# Inputpaths
counts_train_inputpath <- "data/data_partitions/counts_train.rds"
counts_test_inputpath <- "data/data_partitions/counts_test.rds"
metadata_train_inputpath <- "data/data_partitions/metadata_train.rds"
metadata_test_inputpath <- "data/data_partitions/metadata_test.rds"

signatures_inputpath <- "extdata/signatures/signatures_coad.rds"
inputdir <- "extdata/signatures/sigcheck/"

# Outpaths
outputdir <- "data/sigcheck/"

# Arguments
it <- 10000
# CoresToUse <- 7
# mcp <- MulticoreParam(workers = CoresToUse)
# register(mcp, default = TRUE)

# Load data
## TCGA data
counts_train <- readRDS(counts_train_inputpath)
counts_test <- readRDS(counts_test_inputpath)
counts <- rbind.data.frame(counts_train, counts_test)
counts <- t(counts)
metadata_train <- readRDS(metadata_train_inputpath)
metadata_test <- readRDS(metadata_test_inputpath)
metadata <- rbind.data.frame(metadata_train, metadata_test)
rownames(metadata) <- metadata$barcode

## Coad signatures
coad_signatures <- readRDS(signatures_inputpath)

## Signatures to compare (MSigdb + literature)
c1 <- read.gmt(file.path(inputdir, "c1.all.v7.2.symbols.gmt"))
c2 <- read.gmt(file.path(inputdir, "c2.all.v7.2.symbols.gmt"))
c3 <- read.gmt(file.path(inputdir, "c3.all.v7.2.symbols.gmt"))
c4 <- read.gmt(file.path(inputdir, "c4.all.v7.2.symbols.gmt"))
c5 <- read.gmt(file.path(inputdir, "c5.all.v7.2.symbols.gmt"))
c6 <- read.gmt(file.path(inputdir, "c6.all.v7.2.symbols.gmt"))
c7 <- read.gmt(file.path(inputdir, "c7.all.v7.2.symbols.gmt"))
c8 <- read.gmt(file.path(inputdir, "c8.all.v7.2.symbols.gmt"))
literature <- readRDS(file.path(inputdir, "literatureSignatures.rds"))

# Convert ENSEMBL to SYMBOL
genes <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

counts <-
  counts %>%
  as_tibble() %>%
  mutate(
    symbol = genes
  ) %>%
  filter(
    !is.na(symbol)
  ) %>%
  group_by(symbol) %>%
  slice(1) %>%
  tibble::column_to_rownames("symbol") %>%
  as.matrix()

# Create eset object
pdata <- AnnotatedDataFrame(metadata)
eset <- ExpressionSet(
    assayData = counts,
    phenoData = pdata
)

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

# Create sigcheck objects
sigChecks <- list()
for (i in seq_along(signatures_coad_symbol)) {
    sigChecks[[i]] <- sigCheck(
        eset,
        classes = "vital_status",
        survival = "SurvTime",
        signature = signatures_coad_symbol[[i]],
        validationSamples = seq(293, 325, 1),
        scoreMethod = "classifier"
    )
}

names(sigChecks) <- names(signatures_coad_symbol)


# Run SigChecks!
for (i in seq_along(sigChecks)) {

  # retrieve name
  name <- names(sigChecks)[i]

  random <- sigCheckRandom(sigChecks[[i]], iterations = it)
  saveRDS(random, file = paste0(outputdir, "sigCheckRandom_", name, ".rds"))

  perFeat <- sigCheckPermuted(sigChecks[[i]], toPermute = "features", iterations = it)
  saveRDS(perFeat, file = paste0(outputdir, "sigCheckPermutedFeatures_", name, ".rds"))

  perSurv <- sigCheckPermuted(sigChecks[[i]], toPermute = "survival", iterations = it)
  saveRDS(perSurv, file = paste0(outputdir, "sigCheckPermutedSurvival_", name, ".rds"))

  scKnown.lit <- sigCheckKnown(sigChecks[[i]], literature)
  saveRDS(scKnown.lit, file = paste0(outputdir, "scKnownLiterature_", name, ".rds"))

  msigDB <- c(c1, c2, c3, c4, c5, c6, c7, c8)
  scKnown.msigDB <- sigCheckKnown(sigCheck, msigDB)
  saveRDS(scKnown.msigDB, file = paste0(outputdir, "scKnownMSigDB_", name, ".rds"))

}