# Preprocess counts data
# =====
source("requirements.R")
source("R/utils.R")

# Directories
inputpath <- "extdata/TCGA-COAD-RNASeq_SumExp.rds"
extclinpath <- "00_pp/data/clinical.rds"
outputdir <- "00_pp/data"

# Load objects
clinical <- readRDS(extclinpath)
omic <- readRDS(inputpath)

# Get counts data
omic <- assay(omic)
omic <- as.data.frame(t(omic))

# Retrieve only omic data of patient with clinical data available
patIDs <- clinical$barcode
omic <- omic[match(patIDs, rownames(omic)), ]

names(omic) <- sapply(strsplit(names(omic), ".", fixed = TRUE), "[", 1)

# Filter omic data by protein coding only
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
convert <- getBM(filters = c('ensembl_gene_id'), 
                 attributes = c('ensembl_gene_id', 'gene_biotype'),
                 values = names(omic),
                 mart = mart)
protCod <- convert[which(convert$gene_biotype == 'protein_coding'), ]
protCodIDs <- protCod$ensembl_gene_id
omic <- omic[, match(protCodIDs, names(omic))]

# Remove NZV
percZeros <- colSums(omic == 0)/nrow(omic)*100
keepGenes <- names(which(percZeros < 80))
omic <- omic[, match(keepGenes, names(omic))]

# Filter those genes not included Affymetrix GPL570 for validating
platID570 <- "GPL570"

gpl570 <- getGEO(platID570)
gpl570 <- gpl570@dataTable@table
annot570 <- gpl570[, c('ID', 'Gene Symbol', 'Gene Title', 'ENTREZ_GENE_ID')]

dict570 <- getBM(filters = c('affy_hg_u133_plus_2'),
                attributes = c('ensembl_gene_id', 'affy_hg_u133_plus_2'),
                values = annot570$ID,
                mart = mart)

genes570 <- unique(dict570$ensembl_gene_id)
omic570 <- omic[, intersect(names(omic), genes570)]

# Normalize data
counts <- apply(t(omic570), 2, function(x) as.vector(as.integer(x)))
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(condition = clinical$vital_status),
    design = ~ condition
    )
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::estimateDispersions(dds, fitType = "parametric")
vst_data <- DESeq2::vst(dds, blind = TRUE)
vst_matrix <- assay(vst_data)
rownames(vst_matrix) <- colnames(omic570)

# Create list for counts and clinical
gpl570 <- list(omic = as.data.frame(t(vst_matrix)), clinical = clinical)

# Split train and test
smpSize <- floor(0.9 * nrow(clinical))
set.seed(55)
trainIdx <- sample(seq_len(nrow(clinical)), size = smpSize)

# Add to the list
gpl570$omicTrain <- gpl570$omic[trainIdx, ]
gpl570$omicTest <- gpl570$omic[-trainIdx, ]

saveRDS(gpl570, file = "00_pp/data/TCGA-gpl570.rds")
