# Pathway enrichment
# ----------

# Inputpath
tcga_inputpath <- "data/pp/counts_norm_coad_patients.rds"
signatures_inputpath <- "extdata/signatures/signatures_coad.rds"
tcga_counts_inputpath <- "extdata/tcga-coad/TCGA-COAD-RNASeq_counts.rds"

# Outputpath
outputpath <- "data/pathway-enrichment/yaccs_pathway_enrichment.rds"

# Arguments

# Load data
signatures <- readRDS(signatures_inputpath)
train <- readRDS(tcga_inputpath)
counts <- readRDS(tcga_counts_inputpath)
counts <- assays(counts)[[1]]  # unstranded counts to perform DESeq analysis
counts <- t(counts)
colnames(counts) <- sapply(strsplit(colnames(counts), ".", fixed = TRUE), "[", 1)
common_genes <- intersect(colnames(train), colnames(counts))

# Calculate pathway activity
act <- gsva(
    t(train),
    list(signatures[["ours"]]),
    method = "plage",
    kcdf = "Gaussian"
) %>%
t() %>% 
as.data.frame() %>%
rename("yaccs" = "V1")

# Cluster patients based on yaccs score
mclustBIC(act$yaccs) # identify optimal number of clusters

clust <- Mclust(
    act$yaccs,
    modelNames = "E",
    na.rm = TRUE,
    G = 2
)

act$classif <- clust$classification
act$classif <- ifelse(act$classif == 2, "high", "low")

# Save
saveRDS(act, file = outputpath)

# Differential expression analysis
counts[1:5, 1:5]
head(act)
ddse <- DESeqDataSetFromMatrix(
    countData = t(counts),
    colData = act,
    design = ~ classif
)
dds <- DESeq(ddse)
res <- results(dds)

# Pathway enrichment
score <- -log10(res$padj) * sign(res$log2FoldChange)
names(score) <- rownames(res)
phenotype <- score

phenotype_id <- mapIds(
    org.Hs.eg.db,
    keys = names(phenotype),
    column = "ENTREZID",
    keytype = "ENSEMBL"
    )

names(phenotype) <- phenotype_id

# Define Gene Set Collections
reactome <- MSigDBGeneSets(
    species = "Hs",
    collection = "C2",
    subcategory = "REACTOME"
    )
ListGSC <- list(Reactome = reactome)

# Initiate GSCA object
gsca <- GSCA(
    listOfGeneSetCollections = ListGSC,
    geneList = phenotype
    )

# preprocess
gsca1 <- preprocess(
    gsca,
    species = "Hs",
    initialIDs = "ENTREZID",
    keepMultipleMappings = TRUE,
    duplicateRemoverMethod = "max",
    orderAbsValue = FALSE
    )

# analysis
if (requireNamespace("doParallel", quietly = TRUE)) {
    doParallel::registerDoParallel(cores = 6)
}

gsca2 <- analyze(
    gsca1,
    para = list(
        pValueCutoff = 0.05,
        pAdjustMethod = "BH",
        nPermutations = 10000,
        minGeneSetSize = 10,
        exponent = 1),
    doGSOA = FALSE,
    doGSEA = TRUE
    )

# append gene sets terms
gsca3 <- appendGSTerms(gsca2, msigdbGSCs = c("Reactome"))

# Save results
saveRDS(
    list(
        Cluster = act,
        DiffExpression = res,
        PathwayEnrichment = gsca3
        ),
    file = outputpath)
