# 
source("requirements.R")

# Inputpaths
tcga_train_counts_inputpath <- "data/data_partitions/counts_train.rds"
tcga_train_meta_inputpath <- "data/data_partitions/metadata_train.rds"
signatures_coad_inputpath <- "extdata/signatures/signatures_coad.rds"
ccle_gex_inputpath <- "extdata/repurposing/CCLE_expression.csv"
ccle_info_inputpath <- "extdata/repurposing/sample_info.csv"

# Outputpath
outputpath <- "data/repurposing/signatures_scores.rds"

# Load data
tcga_train <- readRDS(tcga_train_counts_inputpath)
tcga_train_meta <- readRDS(tcga_train_meta_inputpath)

signatures_coad <- readRDS(signatures_coad_inputpath)

ccle_info <- read.delim2(
    ccle_info_inputpath,
    header = TRUE,
    row.names = 1,
    sep = ","
    )

ccle_gex <- data.table::fread(
    ccle_gex_inputpath,
    header = TRUE,
    sep = ",") %>%
    tibble::column_to_rownames("V1")

# select coad cell lines
ccle_info <-
    ccle_info %>%
    filter(
        primary_disease == "Colon/Colorectal Cancer"
    )
ccle_gex <- ccle_gex[match(rownames(ccle_info), rownames(ccle_gex)), ]
ccle_gex <- ccle_gex[complete.cases(ccle_gex), ]

colnames(ccle_gex) <- sapply(strsplit(colnames(ccle_gex), " (", fixed = TRUE), "[", 1)


# Loop
i <- 4
res <- list()
for (i in seq_along(signatures_coad)) {

    signature_name <- names(signatures_coad)[i]
    signature_symbol <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = signatures_coad[[i]],
        column = c("SYMBOL"),
        keytype = "ENSEMBL",
        multiVals = "first"
    )
    signature_symbol <- make.names(signature_symbol)
    names(signature_symbol) <- signatures_coad[[i]]

    diff_genes <- setdiff(signature_symbol, colnames(ccle_gex))
    if (length(diff_genes) > 0) {

        to_remove <- which(signature_symbol %in% diff_genes)
        signature_symbol <- signature_symbol[-to_remove]
        signatures_coad[[i]] <- signatures_coad[[i]][-to_remove]

    }

    ccle_gex_f <- ccle_gex[, signature_symbol]
    ccle_gex_f$cohort <- "ccle"

    tcga_train_f <- tcga_train[, signatures_coad[[i]]]
    colnames(tcga_train_f) <- signature_symbol
    tcga_train_f <- scale(tcga_train_f, scale = FALSE) %>% as.data.frame()
    tcga_train_f$cohort <- "tcga"

    # Correct by combat
    df <- rbind.data.frame(tcga_train_f, ccle_gex_f)
    mat <- df %>% dplyr::select(-c(cohort))
    meta <- df %>% dplyr::select(c(cohort))

    modcombat <- model.matrix(~ 1, data = meta)
    combat <- ComBat(
        dat = t(mat),
        batch = meta$cohort,
        mod = modcombat,
        ref.batch = "tcga"
        )

    combat <- cbind.data.frame(t(combat), cohort = meta$cohort)

    train <- combat[which(combat$cohort == "tcga"), ]
    train <- subset(train, select = -c(cohort))
    test <- combat[which(combat$cohort == "ccle"), ]
    test <- subset(test, select = -c(cohort))

    # Train cox model
    train <- cbind.data.frame(
        tcga_train_meta %>% 
            dplyr::select(c(vital_status, SurvTime)) %>%
            mutate(
                vital_status = ifelse(vital_status == "Dead", TRUE, FALSE)
            ),
        train
            )

    cox <- coxph(Surv(SurvTime, vital_status) ~ ., data = train, x = TRUE)
    pred <- predict(cox, test)

    res[[i]] <- data.frame(
        cell_lines = names(pred),
        score = pred,
        signature = signature_name
    )

    print(paste0("Predicted using ", signature_name))

}

# Bind results
res_df <- data.table::rbindlist(res)

# Save
saveRDS(res_df, file = outputpath)